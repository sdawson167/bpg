#include "LamellarPhaseProvider.h"

#define _USE_MATH_DEFINES
#include <math.h>

FieldProvider LamellarPhaseProvider::generateInitialCondition(int gridSize) {

  // initialize field provider variables:

  int* gridSizes = (int*) malloc(sizeof(int));
  gridSizes[0] = gridSize;

  // grid spacing based on period (box size)
  const double gridSpacing = 2 * M_PI / m_period;
  double* dq = (double*) malloc(sizeof(double));
  dq[0] = gridSpacing;

  // initialize field in complex space:
  const bool real = false;

  // initialize data - 1D cosine function:
  fftw_complex* data;
  data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * gridSize);
  for (int index = 0; index < gridSize; index++) {
    data[index][0] = 0.0;
    data[index][1] = 0.0;
  }
  data[0][0] = m_avDensity;
  data[1][0] = m_amplitude;

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

void LamellarPhaseProvider::resetCondition(FieldProvider field)
{
}
