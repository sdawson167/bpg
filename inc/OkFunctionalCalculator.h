#ifndef _OK_CALCULATOR_GUARD
#define _OK_CALCULATOR_GUARD

#include "FieldProvider.h"
#include "fftw3.h"

/*
 * Class contains information about the landau-brazovskii free-energy functional
 */

class OkFunctionalCalculator {

private:
    double  m_tau;
    double  m_gamma;

    const int m_paramNum = 2;

public:
    // constructor
    OkFunctionalCalculator(double tau, double gamma)
    : m_tau{tau},
      m_gamma{gamma}
    { }

    // getter and setter functions
    double getTau() { return m_tau; }
    void setTau(double tau) { m_tau = tau; }

    double getGamma() { return m_gamma; }
    void setGamma(double gamma) { m_gamma = gamma; }

    int getParamNum() { return m_paramNum; }

    // find quadratic potential coefficients, D(q)
    double quadraticCoeff(double q2) 
    { 
      if (q2 == 0) return 0.0;
      else return q2 + 1.0/q2 - 2; 
    }

    // compute quadratic free-energy
    double fQuad(FieldProvider &field);
    double fQuad(fftw_complex* cplxFieldData, double* laplacian, int numFieldElements);

    // compute non-linear free-energy
    double fNL(FieldProvider &field);
    double fNL(fftw_complex* realFieldData, int numFieldElements);

    // compute full free-energy
    double f(FieldProvider &field) {return fNL(field) + fQuad(field);}
    double f(fftw_complex* cplxFieldData, fftw_complex* realFieldData, double* laplacian, int numFieldElements)
    {
      return fQuad(cplxFieldData, laplacian, numFieldElements) + fNL(realFieldData, numFieldElements);
    }

    // compute derivative of NL free-energy and store in field
    void nlDeriv(fftw_complex* realFieldData, fftw_complex* realNLFieldData, int numFieldElements, double &deltaNLDeriv);

    // period optimization method
    double quadraticCoeffDeriv(double q2) { return 1 - 1.0/(q2 * q2); }

    // compute quadratic free-energy for given lattice spacings
    double fQuadB(FieldProvider &field, double* b);

};
#endif
