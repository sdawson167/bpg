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
#include "OkFunctionalCalculator.h"

std::vector<std::string> divideWords(std::string str)
{
  using namespace std;
  string w = "";
  vector<string> strVec;
  for (auto c : str)
  {
    if (c == ' ') {
      strVec.push_back(w);
      w = "";
    } else {
      w = w + c;
    }
  }
  strVec.push_back(w);
  return strVec;
}

int main(int argc, char** argv)
{
  // determine form of input arguments:
  int inMode;
  if (argc > 1) inMode = std::atoi(argv[1]);
  else
  {
    std::cout << "please choose input mode: 1 (input file) or 2 (command line input)" << std::endl;
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
  if (inMode == 1)
  {
    // get name of file
    std::string inString, inFileName;
    if (argc < 3) {
      std::cout << "please enter name of input file" << std::endl;
      return 1;
    } else {
      inString = argv[2];
      inFileName = "/home/sdawson/projects/def-shi/sdawson/BPG/" + inString;
    }

    // open file
    std::ifstream inFile;
    inFile.open(inFileName); 
    if (!inFile) {
      std::cout << "could not read input file " << inFileName << std::endl;
      return 1;
    }

    // get phaseID from remaining input args
    if (argc < 4) {
      std::cout << "please choose one or more phases (1-7) to compute" << std::endl;
      return 1;
    }

    for (int p = 3; p < argc; p++) 
      phaseIdList.push_back(std::atoi(argv[p]));
    
    // read input params from file and store in paramList object
    std::string line;
    while (std::getline(inFile,line)) 
    {
      std::vector<std::string> brokenLine = divideWords(line);
      param otherParams;
      for (size_t index = 0; index < brokenLine.size(); index++) 
	otherParams.push_back(stof(brokenLine[index]));
      
      phasePoints.push_back(otherParams);
    }

  } else if (inMode == 2) {

    if (argc < 3) {
      std::cout << "please enter input params" << std::endl;
      return 1;
    }

    if (argc < 5) {
      std::cout << "incorrect no. of input params for Ok model, need: tau, gamma, phaseID(s)" << std::endl;
      return 1;
    }

    param otherParams;
    for (int index = 2; index < 4; index++) 
      otherParams.push_back(std::atof(argv[index]));
    
    phasePoints.push_back(otherParams);

    for (int index = 4; index < argc; index++)
      phaseIdList.push_back(std::atoi(argv[index]));

  } else {
    std::cout << "input mode must be 1 (read from file) or 2 (read from command line)" << std::endl;
    return 1;
  }
  
  // initialize minimizer object
  BpgMinimizer<OkFunctionalCalculator> minimizer(1e-8, 15);

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

      // initialize calculator
      double tau   = otherParams[0];
      double gamma = otherParams[1];
      OkFunctionalCalculator calculator(tau, gamma);

      std::cout << std::setprecision(10);
      //auto start = std::chrono::steady_clock::now();
      
      try {
        minimizer.minimize(*field, calculator);
        std::cout << calculator.f(*field) << std::endl;
   
      } catch (...) {
        std::cout << 100 << std::endl;
      } 

      /*
      auto end = std::chrono::steady_clock::now();
      std::cout << "elapsed time: " 
	        << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
		<< std::endl;                                                                 
		*/
      // reset initial condition
      provider.resetCondition(*field);

    } // end loop over phase points
    
    delete(field);
  } // end loop over phases

  return 0;
}
