#include "CylindricalHexagonalPhaseProvider.h"

#define _USE_MATH_DEFINES
#include <math.h>

FieldProvider CylindricalHexagonalPhaseProvider::generateInitialCondition(int gridSize) {

  // initialize field provider variables:

  // Ny = 2 * Nx for hexagonal phase
  const int N = gridSize;

  int* gridSizes = (int*) malloc(2 * sizeof(int));
  gridSizes[0] = 2 * N; gridSizes[1] = N;

  const int numFieldElements = 2 * N * N;

  // grid spacing based on period (box size)
  const double dqx = 2 * M_PI / m_period;
  const double dqy = 2 * M_PI / (2 * sqrt(3) * m_period);
  double* dq = (double*) malloc(2 * sizeof(double));
  dq[0] = dqy;  dq[1] = dqx;

  // initialize field in complex space:
  const bool real = false;

  // initialize data -
  // set of fourier peaks from Kai's code:
  std::vector<point> initVals;
  initVals.push_back( makePoint({ 1,  2}, 1.0/6 ));
  initVals.push_back( makePoint({-1,  2}, 1.0/6 ));
  initVals.push_back( makePoint({ 1, -2}, 1.0/6 ));
  initVals.push_back( makePoint({-1, -2}, 1.0/6 ));
  initVals.push_back( makePoint({ 0,  4}, 1.0/6 ));
  initVals.push_back( makePoint({ 0, -4}, 1.0/6 ));

  fftw_complex* data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numFieldElements);
  // initialize values to zero
  for (int index = 0; index < numFieldElements; index++) {
    data[index][0] = 0.0;
    data[index][1] = 0.0;
  }
  for (std::vector<point>::iterator it = initVals.begin(); it != initVals.end(); it++) {
    std::vector<int> k = std::get<0>(*it);	  
    double	     u = std::get<1>(*it);

    int ky = k[1] < 0 ? k[1] + 2 * N : k[1];
    int kx = k[0] < 0 ? k[0] + N : k[0];

    int index = kx + (N * ky);

    data[index][0] = u;
  }

  // create field provider object
  FieldProvider initialCondition{
    data,
    m_dimension,
    gridSizes,
    dq,
    real,
    m_phaseID};

  fftw_free(data);

  return initialCondition;
}

void CylindricalHexagonalPhaseProvider::resetCondition(FieldProvider field)
{
}
