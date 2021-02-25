#include "BccPhaseProvider.h"

#define _USE_MATH_DEFINES
#include <math.h>

FieldProvider BccPhaseProvider::generateInitialCondition(int gridSize) {

  // initialize field provider variables:

  const int N = gridSize;

  int* gridSizes = (int*) malloc(3 * sizeof(int));
  gridSizes[0] = N; gridSizes[1] = N; gridSizes[2] = N;

  const int numFieldElements = N * N * N;

  // grid spacing based on period (box size)
  const double dq = 2 * M_PI / m_period;
  double* dqVec = (double*) malloc(3 * sizeof(double));
  dqVec[0] = dq; dqVec[1] = dq; dqVec[2] = dq;

  // initialize field in complex space:
  const bool real = false;

  fftw_complex* data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numFieldElements);
  populateDataArray(data, numFieldElements, gridSizes);
  
  // create field provider object
  FieldProvider initialCondition{
    data,
    m_dimension,
    gridSizes,
    dqVec,
    real,
    m_phaseID};

  fftw_free(data);

  return initialCondition;
}

void BccPhaseProvider::resetCondition(FieldProvider &field)
{
  // verify phaseID
  const int phaseID = field.getPhaseID();
  if (phaseID != m_phaseID)
    throw std::runtime_error("Could not reset BCC phase - incorrect phase ID");

  // unpack field provider
  int* gridSizes         = field.getGridSizes();
  int  numFieldElements  = field.getNumFieldElements();
  fftw_complex* cplxData = field.getCplxDataPointer();
  
  // set values of cplxData
  populateDataArray(cplxData, numFieldElements, gridSizes);

  // update real data
  field.transformC2R();
}

void BccPhaseProvider::populateDataArray(fftw_complex* cplxData, int numFieldElements, int* gridSizes)
{
  // initialize values to zero
  for (int index = 0; index < numFieldElements; index++) {
    cplxData[index][0] = 0.0;
    cplxData[index][1] = 0.0;
  }

  // set of fourier peaks: 
  const double amp = m_amplitude;
  std::vector<point> initVals;
  initVals.push_back( makePoint( {  0,  0,  0}, m_avDensity));
  initVals.push_back( makePoint( {  1,  1,  0}, amp/12 ));
  initVals.push_back( makePoint( { -1,  1,  0}, amp/12 ));
  initVals.push_back( makePoint( {  1, -1,  0}, amp/12 ));
  initVals.push_back( makePoint( { -1, -1,  0}, amp/12 ));
  initVals.push_back( makePoint( {  1,  0,  1}, amp/12 ));
  initVals.push_back( makePoint( { -1,  0,  1}, amp/12 ));
  initVals.push_back( makePoint( {  1,  0, -1}, amp/12 ));
  initVals.push_back( makePoint( { -1,  0, -1}, amp/12 ));
  initVals.push_back( makePoint( {  0,  1,  1}, amp/12 ));
  initVals.push_back( makePoint( {  0, -1,  1}, amp/12 ));
  initVals.push_back( makePoint( {  0,  1, -1}, amp/12 ));
  initVals.push_back( makePoint( {  0, -1, -1}, amp/12 ));

  // populate appropriate array vals
  for (std::vector<point>::iterator it = initVals.begin(); it != initVals.end(); it++) {
    std::vector<int> k = std::get<0>(*it);
    double	     u = std::get<1>(*it);
                                                                                         
    int kz = k[2] < 0 ? k[2] + gridSizes[0] : k[2];
    int ky = k[1] < 0 ? k[1] + gridSizes[1] : k[1];
    int kx = k[0] < 0 ? k[0] + gridSizes[2] : k[0];
                                                                                         
    int index = kx + (gridSizes[2] * ky) + (gridSizes[2] * gridSizes[1] * kz);
                                                                                         
    cplxData[index][0] = u;
  }
} // end populateDataArray method
