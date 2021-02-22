#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <tuple>
#include <utility>
#include <cstring>
#include <fstream>
#include <sys/time.h>
#include <cmath>

#include "fftw3.h"

#include "BpgMinimizer.h"
#include "FieldProvider.h"
#include "GenericProvider.h"
#include "LamellarPhaseProvider.h"
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

    // read input params from file and store in paramList object
    std::string line;
    while (std::getline(inFile,line)) 
    {
      std::vector<std::string> brokenLine = divideWords(line);
      phaseIdList.push_back(std::stoi(brokenLine[0]));
      param otherParams;
      for (size_t index = 1; index < brokenLine.size(); index++) 
	otherParams.push_back(stof(brokenLine[index]));
      
      phasePoints.push_back(otherParams);
    }

  } else if (inMode == 2) {

    if (argc < 3) {
      std::cout << "please enter input params" << std::endl;
      return 1;
    }

    if (argc != 5) {
      std::cout << "incorrect no. of input params for Ok model" << std::endl;
      return 1;
    }

    phaseIdList.push_back(std::atoi(argv[2]));
    param otherParams;
    for (int index = 3; index < 5; index++) 
      otherParams.push_back(std::atof(argv[index]));
    
    phasePoints.push_back(otherParams);

  } else {
    std::cout << "input mode must be 1 (read from file) or 2 (read from command line)" << std::endl;
    return 1;
  }
  
  // initialize minimizer object
  BpgMinimizer<OkFunctionalCalculator> minimizer(1e-8, 1000);

  // loop through phases
  for (std::vector<int>::const_iterator idIter = phaseIdList.begin(); idIter != phaseIdList.end(); idIter++) { 
    // get phaseID
    int phaseID = *idIter;

    // initialize phase
    GenericPhaseProvider provider(phaseID);
    FieldProvider *field = provider.generateInitialCondition();

    // need to reset initialization if convergence fails
    bool convergenceFlag = true;

    // loop through phase points - compute free-energy at each point
    for (paramList::const_iterator it = phasePoints.begin(); it != phasePoints.end(); it++)
    {
      // check if convergence was achieved at last point - if not, need to reset initial condition
      /*
      if (!convergenceFlag){
	delete(field);
        *field = provider.generateInitialCondition();
      }
      */

      // get phasePoint
      std::vector<double> otherParams = *it;

      // initialize calculator
      double tau   = otherParams[0];
      double gamma = otherParams[1];
      OkFunctionalCalculator calculator(tau, gamma);

      std::cout << std::setprecision(10);
      
      try {
        minimizer.minimize(*field, calculator);
        std::cout << calculator.f(*field) << std::endl;
	convergenceFlag = true;
      } catch (...) {
        std::cout << 100 << std::endl;
	convergenceFlag = false;
      } 
    } // end loop over phase points
    delete(field);
  } // end loop over phases

  return 0;
}
