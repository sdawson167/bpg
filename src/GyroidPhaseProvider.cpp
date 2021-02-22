#include "GyroidPhaseProvider.h"

#define _USE_MATH_DEFINES
#include <math.h>

typedef std::vector<int>    intPoint;
typedef std::tuple<intPoint, double> point;
point makePoint(intPoint coords, double amp) { return point(coords, amp); }

FieldProvider GyroidPhaseProvider::generateInitialCondition(int gridSize) {

  // initialize field provider variables:

  // Ny = 2 * Nx for hexagonal phase
  int* gridSizes = (int*) malloc(3 * sizeof(int));
  gridSizes[0] = gridSize;  gridSizes[1] = gridSize;  gridSizes[2] = gridSize;

  const int numFieldElements = gridSize * gridSize * gridSize;

  // grid spacing based on period (box size)
  const double dx = m_period / gridSize;
  double* dxVec = (double*) malloc(3 * sizeof(double));
  dxVec[0] = dx;  dxVec[1] = dx;  dxVec[2] = dx;

  // initialize field in real space:
  const bool real = true;

  // initialize data -
  fftw_complex* data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numFieldElements);
  // initialize values to zero
  for (int index = 0; index < numFieldElements; index++) {
    data[index][1] = 0.0;
  }
  double w = 2 * M_PI / m_period;
  double sum = 0.0;
  for (int k = 0, index = 0; k < gridSize; k++) {
    double z = k * dx;

    for (int j = 0; j < gridSize; j++) {
      double y = j * dx;

      for (int i = 0; i < gridSize; i++, index++) {
        double x = i * dx;

        double val = sin(w * x) * cos(w * y) +
                     sin(w * y) * cos(w * z) +
                     sin(w * z) * cos(w * x);

        if (val > 1 || val < -1)
          data[index][0] = m_amplitude;
        else
          data[index][0] = 0.0;

        sum += data[index][0];
      }
    }
  }

  double avDensity = sum / numFieldElements;
  for(int index = 0; index < numFieldElements; index++)
    data[index][0] += (m_avDensity - avDensity);

  // create field provider object
  FieldProvider initialCondition{
    data,
    m_dimension,
    gridSizes,
    dxVec,
    real,
    m_phaseID};

  fftw_free(data);

  return initialCondition;
}

void GyroidPhaseProvider::resetCondition(FieldProvider field)
{
}
