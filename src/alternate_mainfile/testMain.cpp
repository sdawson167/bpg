#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <utility>
#include <cstring>
#include <fstream>
#include <cmath>
#include <chrono>

#include "utils.h"
#include "fftw3.h"

#include "BpgMinimizer.h"
#include "FieldProvider.h"
#include "GenericProvider.h"

#include "PwFunctionalCalculator.h"

int main(int argc, char** argv)
{
  /*
   * ===============================================================================
   *                     Parse input args - get tau val
   * ===============================================================================
   */
   
   double tau = 0.0;
   if (argc > 1) { 
     try {
       tau = std::atof(argc[1]);
     } catch (...) {
       std::cout << "invalid tau value" << std::endl;
       return 1;
     }
   }

   int phaseID = 1;
   if (argc > 2) {
     try { 
       phaseID = std::atoi(argc[2]); 
     } catch (...) {
       std::cout << "invalid phaseID" << std::endl;
       return 1;
     }
   }
  
  

  //bool pointInStableRegion(double tau, double gamma, int phaseID);

  /*
   * ==============================================================================
   *               Start optimization for all points listed in input
   * ==============================================================================
   */

  // minimization parameters (max iterations and error tolerance)
  const double errorTol = 1e-9;
  const int    maxIter  = 15;
  
  // ensure we don't print more digits of precion than we know (based on errorTol):
  std::cout << std::fixed << std::setprecision(8);

  // create minimizer (object that handles implementation of minimization algorithm) 
  BpgMinimizer<PwFunctionalCalculator> minimizer(errorTol, maxIter);

  // create field-provider object (represents the order parameter in real and fourier space) 
  // initialize field-provider in chosen phase
  GenericPhaseProvider provider(phaseID); 
  FieldProvider *field = provider.generateInitialCondition();

  // get 
  double gamma
    // loop through phase points (tau, gamma, ... vals) - compute free-energy at each point
    for (paramList::const_iterator it = phasePoints.begin(); it != phasePoints.end(); it++)
    {
      // get phasePoint
      std::vector<double> otherParams = *it;

      // get tau/gamma and lam1/lam2 values
      double tau   = otherParams[0];
      double gamma = otherParams[1];
      double lam1  = otherParams[2];
      double lam2  = otherParams[3];

      // create 'functional calculator' object - handles free-energy calculations
      PwFunctionalCalculator calculator(tau, gamma, lam1, lam2);

      // minimize field, print result (if success) or '100' (if failure)
      try {
        minimizer.minimize(*field, calculator);
        //std::cout << calculator.f(*field) << std::endl;
	fP.push_back(calculator.f(*field));
      } catch (...) {
        //std::cout << 100 << std::endl;
	fP.push_back(100);
      }

      // reset initial condition
      if (resetFlag) 
        provider.resetCondition(*field);

    } // end loop over phase points

    // add list of free-energies to vector
    freeEnergies.push_back(fP);
    
    // clear FieldProvider object for this phase
    delete(field);

  } // end loop over phases
  
  // print free-energies
  for (size_t i = 0; i < phasePoints.size(); i++) {
    std::vector<double> x = phasePoints[i];
    std::cout << x[0] << ", " << x[1];

    for (size_t j = 0; j < phaseIDList.size(); j++) { 
      doubleVec fP = freeEnergies[j];
      std::cout << ", " << fP[i];
    }
    std::cout << std::endl;
  }

  return 0;
}
