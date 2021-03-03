#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <tuple>
#include <utility>
#include <cstring>
#include <fstream>
#include <cmath>
#include <chrono>

#include "fftw3.h"

#include "BpgMinimizer.h"
#include "FieldProvider.h"
#include "GenericProvider.h"

#include "LbFunctionalCalculator.h"
#include "OkFunctionalCalculator.h"
#include "PwFunctionalCalculator.h"

std::vector<std::string> divideWords(std::string str);

int main(int argc, char** argv)
{

  /*
   * =============================================================
   *     Parse input args - determine which calculation to do
   * =============================================================
   */

  // choose model - LB, OK, PW
  // if model is not chosen as input arg, give user option to enter it:
  std::string model;
  if (argc > 1) model = argv[1];
  else {
    std::cout << "choose model: LB, OK, or PW" << std::endl;
    std::getline(std::cin, model);
    std::cin.clear();
  }

  // verify model choice is valid
  if ( model != "LB" && model != "OK" && model != "PW") {
    std::cout << "first input arg. is invalid - user must choose LB, OK, or PW" << std::endl;
    return 1;
  }

  // depending on chosen model paramNumber checks that we are given enough parametsr
  int paramNumber = 0;
  if (model == "LB") paramNumber = 2;      // LB model: tau, gamma
  else if (model == "OK") paramNumber = 2; // OK model: tau, gamma
  else if (model == "PW") paramNumber = 4; // PW model: tau, gamma, lam1, lam2

  // determine form of input arguments - input file (1) or command line (2):
  int inMode;
  if (argc > 2) inMode = std::atoi(argv[2]);
  else
  {
    std::string input;
    std::cout << "choose input mode: 1 (input file) or 2 (command line input)" << std::endl;
    std::getline(std::cin, input);
    
    std::stringstream stream(input);
    if (!(stream >> inMode)) {
      std::cout << "second input arg is invalid - choose 1 or 2" << std::endl;
      return 1;
    }
    std::cin.clear();
  }
  // verify choice of input mode is valid
  if (inMode != 1 && inMode != 2) {
    std::cout << "second input arg is invalid - choose 1 or 2" << std::endl;
    return 1;
  }

  /* initialize params from remaining input arguments
   * how we do this depends on choice of input mode:
   * 	1 (input file) - need name of file from cmd line
   * 	2 (command line) - need remaining params from cmd line
   */

  // define variable type to hold phase point data
  typedef std::vector<double> param;  
  typedef std::vector<param> paramList;

  // initialize all points on which calculation will be done
  paramList phasePoints;
  std::vector<int> phaseIdList;

  // if we chose to initialize from file:
  if (inMode == 1)
  {
    // get name of file
    std::string inString, inFileName;
    if (argc < 4) {
      std::cout << "enter name of input file" << std::endl;
      std::getline(std::cin, inString);
      std::cin.clear();
    } else 
      inString = argv[3];
      
    inFileName = "/home/sdawson/projects/def-shi/sdawson/BPG/" + inString;

    // open file
    std::ifstream inFile;
    inFile.open(inFileName); 
    // verify open:
    if (!inFile) {
      std::cout << "could not read input file " << inFileName << std::endl;
      return 1;
    }
    
    // read input params from file and store in paramList object
    std::string line;
    while (std::getline(inFile,line)) 
    {
      std::vector<std::string> brokenLine = divideWords(line);
      if ((int) brokenLine.size() != paramNumber) {
        std::cout << "invalid parameter file for " << model << " model, file must provide " << paramNumber << " input params per line" << std::endl;
        return 1;
      }
      param otherParams;
      for (size_t index = 0; index < brokenLine.size(); index++) 
        otherParams.push_back(stof(brokenLine[index]));
      
      phasePoints.push_back(otherParams);
    }

    // now get phaseIDs - list of phases to optimize:
    if (argc < 5 ) {
      std::cout << "choose phases (1 - 7) to optimize, enter any non-numerical character to finish:" << std::endl;

      int phaseID;
      std::string input;
      bool readInFlag = true;
      while (std::getline(std::cin, input)) {
    
        std::vector<std::string> dividedInput = divideWords(input);
    
        for (size_t index = 0; index < dividedInput.size(); index++) {
          std::stringstream stream(dividedInput[index]);
                                                                                                               
          if (stream >> phaseID) {
            if (phaseID > 0 && phaseID < 8)
              phaseIdList.push_back(phaseID);
          } else 
            readInFlag = false;
        }
        if (!readInFlag)
          break;
      }
  
      // make sure we have at least one phase:  
      if (phaseIdList.size() == 0) {
        std::cout << "user must choose 1 or more phases to optimize" << std::endl;
        return 1;
      }
      std::cin.clear();
                                                                                                               
    } else {
      for (int p = 4; p < argc; p++) 
        phaseIdList.push_back(std::atoi(argv[p]));
    }

  // end inMode = 1 case 

  // inMode = 2: read input params (tau, gamma, ...) and list of phaseIDs from command line:
  } else if (inMode == 2) {
    
    // first, get input params (tau, gamma, ...):
    param otherParams;
    if (argc < (3 + paramNumber) ) {
      std::cout << "enter " << paramNumber << " input params required by " << model << " model:" << std::endl;
      
      double x;
      std::string input;
      int n = 0;
      while (std::getline(std::cin, input)) {
        
	std::vector<std::string> dividedInput = divideWords(input);
	
	for (size_t index = 0; index < dividedInput.size(); index++) {
          std::stringstream stream(dividedInput[index]);
	  
	  if (stream >> x) {
	    if (n < paramNumber) {
	      otherParams.push_back(x);
	      n++;
	    } 
	  
	  } 
	}

	if (n == paramNumber)
          break;
      }
      
      std::cin.clear();
    
    } else {
      for (int index = 3; index < (3 + paramNumber); index++) 
      otherParams.push_back(std::atof(argv[index])); 
    }
    phasePoints.push_back(otherParams);

    // now get phaseIDs - list of phases to optimize:
    if (argc < (3 + paramNumber + 1) ) {
      std::cout << "choose phases (1 - 7) to optimize, enter any non-numerical character to finish:" << std::endl;
    
      int phaseID;
      std::string input;
      bool readInFlag = true;
      while (std::getline(std::cin, input)) {
        
	std::vector<std::string> dividedInput = divideWords(input);
	
	for (size_t index = 0; index < dividedInput.size(); index++) {
	  std::stringstream stream(dividedInput[index]);

	  if (stream >> phaseID) {
	    if (phaseID > 0 && phaseID < 8)
              phaseIdList.push_back(phaseID);
	  } else 
	    readInFlag = false;
	}
	if (!readInFlag)
	  break;
      }
      
      // make sure we have at least one phase:  
      if (phaseIdList.size() == 0) {
        std::cout << "user must choose 1 or more phases to optimize" << std::endl;
        return 1;
      }
      std::cin.clear();

    } else {
      for (int p = (3 + paramNumber); p < argc; p++)
        phaseIdList.push_back(std::atoi(argv[p]));
    }
  } // end inMode = 2 case

  // loop through phases
  for (std::vector<int>::const_iterator idIter = phaseIdList.begin(); idIter != phaseIdList.end(); idIter++) { 
    // get phaseID
    int phaseID = *idIter;

    // create field-provider object and initialize phase
    GenericPhaseProvider provider(phaseID);
    FieldProvider *field = provider.generateInitialCondition();

    // loop through phase points - compute free-energy at each point
    for (paramList::const_iterator it = phasePoints.begin(); it != phasePoints.end(); it++)
    {
      // get phasePoint
      std::vector<double> otherParams = *it;

      // get tau/gamma values - used for all models!
      double tau   = otherParams[0];
      double gamma = otherParams[1];

      // minimization parameters - used for all models
      double errorTol = 1e-8;
      int    maxIter  = 15;
      // ensure we don't print more digits of precision than we know:
      std::cout << std::fixed;
      std::cout << std::setprecision(8);

      // code timing stuff
      // auto start = std::chrono::steady_clock::now();

      // this part depends on the choice of model
      
      if (model == "LB") {
        // create calculator object
	LbFunctionalCalculator calculator(tau, gamma);

	// create minimizer object
	BpgMinimizer<LbFunctionalCalculator> minimizer(errorTol, maxIter);

	// minimize field, print result (if success) or '100' (if failure)
	try {
	  minimizer.minimize(*field, calculator);
	  std::cout << calculator.f(*field) << std::endl;
	} catch (...) {
	  std::cout << 100 << std::endl;
	}
      
      } else if (model == "OK") {
        // create calculator object
	OkFunctionalCalculator calculator(tau, gamma);

	// create minimizer object
	BpgMinimizer<OkFunctionalCalculator> minimizer(errorTol, maxIter);

	// minimize field, print result (if success) or '100' (if failure)
	try {
	  minimizer.minimize(*field, calculator);
	  std::cout << calculator.f(*field) << std::endl;
	} catch (...) {
	  std::cout << 100 << std::endl;
	}

      } else if (model == "PW") {
        // read remaining input params
	double lam1 = otherParams[2];
	double lam2 = otherParams[3];

	// create calculator object
	PwFunctionalCalculator calculator(tau, gamma, lam1, lam2);

	/*
	// create list of q-values
	std::vector<double> q;
	double dq = 1e-2;
	int N = 401;
	for (int index = 0; index < N; index++) {
	  q.push_back(dq * index);
	  std::cout << index << ": " << dq * index << std::endl;
	}
	*/

	// create minimizer object
	BpgMinimizer<PwFunctionalCalculator> minimizer(errorTol, maxIter);

	// minimize field, print result (if success) or '100' (if failure)
	try {
	  minimizer.minimize(*field, calculator);
	  std::cout << calculator.f(*field) << std::endl;
	} catch (...) {
	  std::cout << 100 << std::endl;
	}

      } 

      // timing stuff:
      /*
      auto end = std::chrono::steady_clock::now();
      std::cout << "elapsed time: " 
	        << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
		<< std::endl;                                                                 
		*/
      
      // reset initial condition
      provider.resetCondition(*field);

    } // end loop over phase points
    
    // clear FieldProvider object for this phase
    delete(field);

  } // end loop over phases
  
  return 0;
}

/*
 * ======================================================
 *            method to parse input files
 * ======================================================
 *
 *   method takes string, divides it into words 
 *   (separated by spaces) and returns a vector of 
 *   those words 
 *
 */

std::vector<std::string> divideWords(std::string str)
{
  using namespace std;
  string w = "";  // current word
  vector<string> strVec;

  // parse each char in the string:
  for (auto c : str)
  {
    // if we've found a space, divide the sentence here
    if (c == ' ') {
      strVec.push_back(w);
      w = "";
    // otherwise add character to current word
    } else {
      w = w + c;
    }
  }
  strVec.push_back(w);
  return strVec;
}

