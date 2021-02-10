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

  // initialize data -
  // set of fourier peaks from Kai's code:
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

  fftw_complex* data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numFieldElements);
  // initialize values to zero
  for (int index = 0; index < numFieldElements; index++) {
    data[index][0] = 0.0;
    data[index][1] = 0.0;
  }

  for (std::vector<point>::iterator it = initVals.begin(); it != initVals.end(); it++) {
    std::vector<int> k = std::get<0>(*it);
    double	     u = std::get<1>(*it);

    int kz = k[2] < 0 ? k[2] + N : k[2];
    int ky = k[1] < 0 ? k[1] + N : k[1];
    int kx = k[0] < 0 ? k[0] + N : k[0];

    int index = kx + (N * ky) + (N * N * kz);

    data[index][0] = u;
  }

  // create field provider object
  FieldProvider initialCondition{
    data,
    m_dimension,
    gridSizes,
    dqVec,
    real,
    m_phaseDimension};

  fftw_free(data);

  return initialCondition;
}
