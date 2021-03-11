#ifndef _UTILS_GUARD
#define _UTILS_GUARD

#include <vector>
#include <string>

// header file contains a couple useful methods used in the main code

// method to break string up into list of words (separated by spaces) and return
// vector of words
std::vector<std::string> divideWords(std::string str);

// method to read arguments from command line and extract list of phase points
// (tau, gamma, ... values) and phaseIDs (integer indicators of phases to
// minimize)
typedef std::vector<double> param;  
typedef std::vector<param> paramList;
void parseInputArgs(int argc, char** argv, int &inMode, paramList &phasePoints, std::vector<int> &phaseIDList);

#endif
