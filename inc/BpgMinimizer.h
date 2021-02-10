#ifndef _BPG_MINIMIZER_GUARD
#define _BPG_MINIMIZER_GUARD

#include <stdlib.h>
#include <stdexcept>
#include <cmath>

#include "fftw3.h"
#include "FieldProvider.h"
#include "LbFunctionalCalculator.h"

/* class uses the BPG minimizer algorithm described by Jiang et. al. (2020), arXiv:2002.09898, to minimize the free-energy
 * of a given intial condition with respect to a given free energy functional.
 *
 * TFunctionalCalculator represents the free-energy functional we are minimizing over
 */
template<typename TFunctionalCalculator>
class BpgMinimizer {
private:
    // toggle algorithm components
    bool   m_useNesterov;         // use Nesterov acceleration
    bool   m_adaptiveTimeSteps;   // use adaptive time-steps
    double m_errorTolerance;      // used for convergence test
    int    m_maxIterations;       // used to abort loops
    int    m_iterator;            // used to track total iterations
    int    m_fieldIterator;       // no. of iterations used by field minimizer
    int    m_periodIterator;      // no. of iterations used by period optimizer

public:
    // constructor
    BpgMinimizer(bool useNesterov, bool useAdaptiveTimeSteps, double errorTolerance, int maxIterations)
    : m_useNesterov{useNesterov},
      m_adaptiveTimeSteps{useAdaptiveTimeSteps},
      m_errorTolerance{errorTolerance},
      m_maxIterations{maxIterations}
      {
        // initialize iterator at zero !
        m_iterator = 0;
        m_fieldIterator = 0;
        m_periodIterator = 0;
      }

      // getter / setter for variables
      void setNesterov(bool useNesterov) { m_useNesterov = useNesterov; }
      bool getNesterov() { return m_useNesterov; }

      void setATS(bool useATS) { m_adaptiveTimeSteps = useATS; }
      bool getATS() { return m_adaptiveTimeSteps; }

      void setErrorTol(double errorTol) { m_errorTolerance = errorTol; }
      double getErrorTol() { return m_errorTolerance; }

      void setMaxIterations(int maxIterations) { m_maxIterations = maxIterations; }
      int getMaxIterations() { return m_maxIterations; }

      void setIterator(int iterator) { m_iterator = iterator; }
      int getIterator() { return m_iterator; }


    // minimizer function - uses SIS minimizer with (optional) Nesterov acceleration technique and
    // (optional) adaptive time-stepping (NOTE: adaptive time-stepping has not been implemented yet)
    // to find field configuration that minimizes free-energy represented by functional calculator
    void minimizeField(
      FieldProvider &field,
      TFunctionalCalculator &calculator)
    {
      m_fieldIterator = 0;

      // get information about field
      fftw_complex* cplxFieldData = field.getCplxDataPointer();
      fftw_complex* realFieldData = field.getRealDataPointer();
      const int dimension         = field.getDimension();
      int* gridSizes              = field.getGridSizes();
      const int numFieldElements  = field.getNumFieldElements();
      double* dx 		  = field.getDx();
      double* dq                  = field.getDq();
      const int phaseDimension    = field.getPhaseDimension();

      
      // initialize nesterov field - copy field data
      int* nestGridSizes = (int*) malloc(sizeof(int) * dimension);
      memcpy(nestGridSizes, gridSizes, sizeof(int) * dimension);

      double* nestDq = (double*) malloc(sizeof(double) * dimension);
      memcpy(nestDq, dq, sizeof(double) * dimension); 

      FieldProvider nestField(
        cplxFieldData, //initNestData,
        dimension,
        nestGridSizes,
        nestDq,
        false,
        phaseDimension);
      fftw_complex* realNestData = nestField.getRealDataPointer();
      fftw_complex* cplxNestData = nestField.getCplxDataPointer();

      // NL derivative field - compute Nl derivative and create field provider
      fftw_complex* initNlDerivData = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numFieldElements);
      double deltaNLVal = 0.0;
      calculator.nlDeriv(realNestData, initNlDerivData, numFieldElements, deltaNLVal);

      /*
      if (dimension == 1) {
        // initialize fftw stuff
        fftw_complex *in, *out;
        in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numFieldElements);
        initNlDerivData = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numFieldElements);

        for (int index = 0; index < numFieldElements; index++) {
          in[index][0] = 0;
          in[index][1] = 0;
          initNlDerivData[index][0] = 0;
          initNlDerivData[index][1] = 0;
        }

        fftw_plan p = fftw_plan_dft_1d(numFieldElements, in, initNlDerivData, FFTW_FORWARD,FFTW_ESTIMATE);

        // load data into in
        double temp;
        calculator.nlDeriv(realNestData, in, numFieldElements, temp);

        // execute
        fftw_execute(p);

        // delete everything
        fftw_destroy_plan(p);
        fftw_free(in);
      }
      */

      int* nlDerivGridSizes = (int*) malloc(sizeof(int) * dimension);
      memcpy(nlDerivGridSizes, gridSizes, sizeof(int) * dimension);

      double* nlDerivDx = (double*) malloc(sizeof(double) * dimension);
      memcpy(nlDerivDx, dx, sizeof(double) * dimension);

      FieldProvider nlDerivField(
        initNlDerivData,
        dimension,
        nlDerivGridSizes,
        nlDerivDx,    // note for real-field initialization we need real grid spacing
        false,
        phaseDimension);
      fftw_free(initNlDerivData);
      fftw_complex* realNlDerivData = nlDerivField.getRealDataPointer();
      fftw_complex* cplxNlDerivData = nlDerivField.getCplxDataPointer();

      /*
      // timestep field
      fftw_complex* initTimeStepData = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numFieldElements);
      memcpy(initTimeStepData, cplxFieldData, sizeof(fftw_complex) * numFieldElements);
      FieldProvider timeStepField(
        initTimeStepData,
        dimension,
        gridSizes,
        dq,
        false,
        phaseDimension);
      fftw_complex* timeStepData = timeStepField.getCplxDataPointer();
      */

      // allocate memory for laplacian and initialize
      double* laplacian = (double*) malloc(sizeof(double) * numFieldElements);
      field.laplacian(laplacian);

      // allocate memory for quadratic coefficient array and initialize
      double* D = (double*) malloc(sizeof(double) * numFieldElements);
      for (int index = 0; index < numFieldElements; index++)
        D[index] = calculator.quadraticCoeff(-laplacian[index]);

      // allocate memory for storing old nesterov field data
      fftw_complex* oldNestField = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numFieldElements);

      // compute free-energy
      double f0 = calculator.f(cplxFieldData, realFieldData, laplacian, numFieldElements);
      double fOld = f0;

      // initialize iteration variables

      // timestep and adaptive time-stepping variables
      const double tMin = 0.1;      // adaptive - never allow timestep to get smaller than this
      const double tMax = 2.0;      // adaptive - never allow timestep to get larger than this
      const double t0 = 1.0;        // adaptive ? initial timestep : timestep
      double tStep = t0;
      //const double rho    = 0.5 * (sqrt(5) - 1.0); // update coefficient (rho < 1)
      //const double delta  = 0.01;                              // used for timestep update inequality test

      // Nesterov step variables
      double       theta  = 1.0;
      const double thetaQ = 0.0;
      double       beta   = 0.0;

      // error
      double currentError = 1.0;

      // stop criterion: stop loop when error is below tolerance OR max iterations reached
      bool stopCriterion = false;

      // while loop 

      while(!stopCriterion)
      {
        // increment iterator
        m_fieldIterator++;

	// update field :
        // double deltaNest = 0.0; // needed for time-step estimation
        for (int index = 0; index < numFieldElements; index++) {

          // store current field values
          double reCplxFieldOldVal = cplxFieldData[index][0];
          double imCplxFieldOldVal = cplxFieldData[index][1];

          // extract Nl deriv values:
          double reCplxNlDerivVal  = cplxNlDerivData[index][0] / numFieldElements;
          double imCplxNlDerivVal  = cplxNlDerivData[index][1] / numFieldElements;

          // matrix element for update step
          double A = 1.0 - tStep * laplacian[index] * D[index];

          // update cplx field
          cplxFieldData[index][0] = (cplxNestData[index][0] + tStep * laplacian[index] * reCplxNlDerivVal) / A;
          cplxFieldData[index][1] = (cplxNestData[index][1] + tStep * laplacian[index] * imCplxNlDerivVal) / A;

          // store current Nesterov field values
          double reNestFieldOldVal = cplxNestData[index][0];
          double imNestFieldOldVal = cplxNestData[index][1];
          oldNestField[index][0] = realNestData[index][0];

          // update cplx Nesterov field
          double reNestFieldVal, imNestFieldVal;
          if (m_useNesterov) {
            reNestFieldVal = (1.0 + beta) * cplxFieldData[index][0] - beta * reCplxFieldOldVal;
            imNestFieldVal = (1.0 + beta) * cplxFieldData[index][1] - beta * imCplxFieldOldVal;
          } else {
            reNestFieldVal = cplxFieldData[index][0];
            imNestFieldVal = cplxFieldData[index][1];
          }
          cplxNestData[index][0] = reNestFieldVal;
          cplxNestData[index][1] = imNestFieldVal;

          // calculate change in nest field (used for time-step estimation)
          //deltaNest += (reNestFieldVal - reNestFieldOldVal) * (reNestFieldVal - reNestFieldOldVal) +
          //             (imNestFieldVal - imNestFieldOldVal) * (imNestFieldVal - imNestFieldOldVal);
        }

        // transform field and Nesterov field back into real space
        field.transformC2R();
        nestField.transformC2R();

        // update NL derivative field in real space and transform
	calculator.nlDeriv(realNestData, realNlDerivData, numFieldElements, deltaNLVal);
	nlDerivField.transformR2C();

        // compute new free-energy
        double f     = calculator.f(cplxFieldData, realFieldData, laplacian, numFieldElements);
        //double fNest = calculator.f(cplxNestData, realNestData, laplacian, numFieldElements);

        // compute error
        currentError = fOld - f;
        fOld = f;

        // update Nesterov coefficient beta
	// reset condition - free-energy increases:
        if (currentError < 0 && m_useNesterov) {
          theta = 1.0;
        }
        double theta2   = theta * theta;
        double newTheta = -0.5 * (theta2 - thetaQ) + sqrt(0.25 * (theta2 - thetaQ) * (theta2 - thetaQ) + theta2);
        beta            = theta * (1 - theta) / (theta2 + newTheta);
        theta           = newTheta;

        // update stopCriterion
        if (m_fieldIterator >= m_maxIterations || std::abs(currentError) < m_errorTolerance)
	  stopCriterion = true;

      } // end while loop

      // check no. of iterations
      if (m_fieldIterator == m_maxIterations)
        throw std::runtime_error("Maximum iterations reached in field optimization");

      free(laplacian);
      free(D);
      fftw_free(oldNestField);
    } // end minimizeField method

    /*
     * Collection of methods for period optimization
     */

    // method computes optimal lattice vector size(s) using conjugate gradient method
    // and uses them to update the real and cplx lattice spacing of field
    void optimizePeriods(
      FieldProvider &field,
      TFunctionalCalculator &calculator)
    {
      m_periodIterator = 0;

      // unpack variables in field provider a little bit
      fftw_complex* cplxFieldData = field.getCplxDataPointer(); // cplx coefficients phi_n

      const int d = field.getDimension();

      const int numFieldElements = field.getNumFieldElements();
      const int* N = field.getGridSizes();

      // initial reciprocal lattice spacing
      double *dQ = field.getDq();

      // initialize arrays to hold b_k and b_{k-1} :
      double* bOld = (double*) calloc(3, sizeof(double));
      double* bNew = (double*) calloc(3, sizeof(double));
      // reorder elements of lattice spacing into 3-component vector b_k:
      bNew[0]   = dQ[d - 1];
      if (d > 1) {
        bNew[1]   = dQ[d - 2];
        if (d > 2)
          bNew[2] = dQ[d - 3];
      }

      // initialize arrays to hold residuals and step direction
      double* rOld = (double*) calloc(3, sizeof(double));
      double* rNew = (double*) calloc(3, sizeof(double));
      double* s    = (double*) calloc(3, sizeof(double));

      // compute residual and magnitude:
      double gradMagnitude = computeResidual(rNew, bNew, numFieldElements, N, d, calculator, cplxFieldData);

      // initialize conjugate coefficient
      double beta = 0.0;

      // iteration variables
      double f0 = calculator.fQuad(field), f; // used to verify free-energy decrease every loop
      bool posFlag = false;
      bool stopCriterion = gradMagnitude < m_errorTolerance;  // we may already be at optimal unit cell sizing!

      while(!stopCriterion) {
        m_periodIterator++;

        // compute line search direction s:
        s[0] = rNew[0] + beta * rOld[0];
        s[1] = rNew[1] + beta * rOld[1];
        s[2] = rNew[2] + beta * rOld[2];

        // compute optimal step-size, alpha
        double alpha = 1.0;
        findAlpha(alpha, calculator, field, bOld, s, -gradMagnitude, posFlag);

        // store step
        bOld[0] = bNew[0];
        bOld[1] = bNew[1];
        bOld[2] = bNew[2];

        // update b
        bNew[0] = bOld[0] + alpha * s[0];
        bNew[1] = bOld[1] + alpha * s[1];
        bNew[2] = bOld[2] + alpha * s[2];

        // store residuals:
        rOld[0] = rNew[0];
        rOld[1] = rNew[1];
        rOld[2] = rNew[2];

        // update residuals:
        double gradMagnitudeOld = gradMagnitude;
        gradMagnitude = computeResidual(rNew, bNew, numFieldElements, N, d, calculator, cplxFieldData);

        // compute conjugate coefficient for next time (Fletcher-Reeves) - use restart to prevent oscillation:
        if (gradMagnitude > gradMagnitudeOld) beta = 0.0;
        else beta = gradMagnitude / gradMagnitudeOld;

        // compute change in free-energy
        double f = calculator.fQuadB(field, bNew);
        double err = f - f0;
        if (err > 0)
          posFlag = true;
        f0 = f;

        // update stop criterion
        if (abs(err) < m_errorTolerance || m_periodIterator == m_maxIterations)
          stopCriterion = true;
      }

      // verify max iterations have not been reached
      if (m_periodIterator == m_maxIterations)
        throw std::runtime_error("Maximum iterations reached in period optimization");

      // repackage into format that field-provider likes
      double* newDq = (double*) malloc(sizeof(double) * d);
      newDq[d - 1] = bNew[0];
      if (d > 1) {
        newDq[d - 2] = bNew[1];
        if (d > 2)
          newDq[d - 3] = bNew[2];
      }

      // update field provider lattice spacings:
      field.setDq(newDq);

      free(bOld);
      free(bNew);

      free(rOld);
      free(rNew);
      free(s);
    } // end period optimization method

    double computeResidual(
      double* r,
      double* b,
      const int numFieldElements,
      const int* N,
      const int d,
      TFunctionalCalculator &calculator,
      fftw_complex* cplxFieldData)
    {
      r[0] = 0.0;
      r[1] = 0.0;
      r[2] = 0.0;

      // loop through all reciprocal lattice points
      for (int index = 0; index < numFieldElements; index++) {

        // determine x, y, z coordinates of lattice point:
        int i{0}, j{0}, k{0};
        int iMod{0}, jMod{0}, kMod{0};

        int Nx = N[d-1];
        i = index % Nx;
        iMod = i < Nx / 2 ? i : Nx - i;

        if (d > 1) {
          int Ny = N[d-2];
          j = (index - i) / Nx % Ny;
          jMod = j < Ny / 2 ? j : Ny - j;

          if (d > 2) {
            int Nz = N[d-3];
            k = ((index - i) / Nx - j) / Ny;
            kMod = k < Nz / 2 ? k : Nz - k;
          }
        }

        // calculate reciprocal lattice vector at this point
        double q2 =   (iMod * iMod) * (b[0] * b[0])
                    + (jMod * jMod) * (b[1] * b[1])
                    + (kMod * kMod) * (b[2] * b[2]);

        // derivative of Gamma(q2)
        double gammaPrime = calculator.quadraticCoeffDeriv(q2);

        // amplitude of fourier peak:
        double phi2 = cplxFieldData[index][0] * cplxFieldData[index][0] + cplxFieldData[index][1] * cplxFieldData[index][1];

        r[0] += -gammaPrime * phi2 * iMod * iMod * b[0];
        r[1] += -gammaPrime * phi2 * jMod * jMod * b[1];
        r[2] += -gammaPrime * phi2 * kMod * kMod * b[2];
      } // end loop over reciprocal lattice points

      return (r[0] * r[0]) + (r[1] * r[1]) + (r[2] * r[2]);
    } // end computeResiduals method

    // returns step size found via backtracking line search method:
    void findAlpha(
      double &alpha,
      TFunctionalCalculator &calculator,
      FieldProvider &field,
      double* bk,
      const double* sk,
      const double m,
      bool posFlag)
    {
      // set alpha to start
      alpha = 1.0;  // max value of alpha

      // initial free-energy:
      double fB0 = calculator.fQuadB(field,bk);

      // algorithm variables:
      const double c = 0.2;
      const double rho = 0.9;

      double alphaMin = 1e-1;
      if (m_periodIterator > 50 || posFlag) alphaMin = 1e-2;

      double* b = (double*) malloc(sizeof(double) * 3);
      
      bool stopCriterion = false;
      
      while (!stopCriterion) {
        // update b vectors:
        b[0] = bk[0] + alpha * sk[0];
        b[1] = bk[1] + alpha * sk[1];
        b[2] = bk[2] + alpha * sk[2];

        double lhs = fB0 - calculator.fQuadB(field,bk);
        double rhs = - alpha * c * m;

        if (lhs >= rhs || alpha < alphaMin) {
          stopCriterion = true;
          alpha = alphaMin;
        }
        else alpha *= rho;
      }
      free(b);
    } // end findAlpha method

    // optimize field config and period
    void findMinimum(
      FieldProvider &field,
      TFunctionalCalculator &calculator)
    {
      // set iterator to zero
      m_iterator = 0;

      // loop until error is smaller than tolerance or max no. of iterations are reached
      bool stopCriterion = false;

      // error - computed using change in free-energy:
      double f0 = calculator.f(field); // intial free-energy
      double error = 1.0;

      while (!stopCriterion)
      {
        m_iterator++;

        // minimize field
        minimizeField(field, calculator);

        // minimize periods
        optimizePeriods(field, calculator);

        // recompute free-energy
        double f = calculator.f(field);

        // update error
        error = f - f0;

        // store free-energy of current config
        f0 = f;

        // update stop criterion
        if (abs(error) < m_errorTolerance || m_iterator == m_maxIterations)
          stopCriterion = true;
      }

      // check if max iterations have been reached
      if (m_iterator == m_maxIterations)
        throw std::runtime_error("Maximum iterations reached");
    } // end find minimum method
};
#endif
