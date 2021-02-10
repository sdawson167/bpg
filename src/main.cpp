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
  typedef std::pair<int, std::vector<double>> param;  
  typedef std::vector<param> paramList;

  // initialize all points on which calculation will be done
  paramList phasePoints;
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
      int phaseID = std::stoi(brokenLine[0]);
      std::vector<double> otherParams;
      for (size_t index = 1; index < brokenLine.size(); index++) 
	otherParams.push_back(stof(brokenLine[index]));
      
    param phasePoint(phaseID, otherParams);
    phasePoints.push_back(phasePoint);
    }

  } else if (inMode == 2) {

    if (argc < 3) {
      std::cout << "please enter input params" << std::endl;
      return 1;
    }

    if (argc != 5) {
      std::cout << "incorrect no. of input params for Lb model" << std::endl;
      return 1;
    }

    int phaseID = std::atoi(argv[2]);
    std::vector<double> otherParams;
    for (int index = 3; index < 5; index++) 
      otherParams.push_back(std::atof(argv[index]));
    
    param phasePoint(phaseID, otherParams);
    phasePoints.push_back(phasePoint);

  } else {
    std::cout << "input mode must be 1 (read from file) or 2 (read from command line)" << std::endl;
    return 1;
  }
  
  // initialize minimizer object
  BpgMinimizer<OkFunctionalCalculator> minimizer(true, false, 1e-8, 1000);

  // loop through phase points - compute free-energy at each point
  for (paramList::const_iterator it = phasePoints.begin(); it != phasePoints.end(); it++)
  {
    // get phasePoint
    int phaseID = (*it).first;
    std::vector<double> otherParams = (*it).second;

    // initialize field from phaseID
    GenericPhaseProvider provider(phaseID);
    FieldProvider *field = provider.generateInitialCondition();

    // initialize calculator
    double tau   = otherParams[0];
    double gamma = otherParams[1];
    OkFunctionalCalculator calculator(tau, gamma);

    // minimize field
    std::cout << std::setprecision(10);
    try {
      minimizer.findMinimum(*field, calculator);
      std::cout << calculator.f(*field) << std::endl;
    } catch (...) {
      //std::cout << 100 << std::endl;
    }
    
    delete(field);
  } // end loop over phase points

  return 0;
}
