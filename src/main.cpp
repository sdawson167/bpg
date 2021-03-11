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
   *             Parse input args - determine which calculation to do
   * ===============================================================================
   */
  
  int inMode;
  paramList phasePoints;
  std::vector<int> phaseIDList;

  // attempt to extract list of calculations to perform 
  try {
    parseInputArgs(argc, argv, inMode, phasePoints, phaseIDList);

  } catch (std::string errorMessage){
    std::cout << errorMessage << std::endl;
    return 1;
  }

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

  // loop through phases
  for (std::vector<int>::const_iterator idIter = phaseIDList.begin(); idIter != phaseIDList.end(); idIter++) { 
   
    // get phaseID - integer that indicates which phase we are optimizing
    int phaseID = *idIter;

    // create field-provider object (represents the order parameter in real and fourier space) 
    // initialize field-provider in chosen phase
    GenericPhaseProvider provider(phaseID); 
    FieldProvider *field = provider.generateInitialCondition();

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
        std::cout << calculator.f(*field) << std::endl;
      } catch (...) {
        std::cout << 100 << std::endl;
      }

      // reset initial condition
      if (inMode == 1 || inMode == 2) {
        std::cout << "resetting provider" << std::endl;
        provider.resetCondition(*field);
      }

    } // end loop over phase points
    
    // clear FieldProvider object for this phase
    delete(field);

  } // end loop over phases
  
  return 0;
}
