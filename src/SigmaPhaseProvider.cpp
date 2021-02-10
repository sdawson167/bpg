#include "SigmaPhaseProvider.h"

#define _USE_MATH_DEFINES
#include <math.h>

FieldProvider SigmaPhaseProvider::generateInitialCondition(int gridSize) {

  // initialize field provider variables:

  const int N  = 2 * gridSize;
  const int Nz = gridSize;
  int* gridSizes = (int*) malloc(3 * sizeof(int));
  gridSizes[0] = Nz;  gridSizes[1] = N;  gridSizes[2] = N;
  const int numFieldElements = Nz * N * N;

  // grid spacing based on period (box size)
  const double dqx = 2 * M_PI/m_periodX;
  const double dqz = 2 * M_PI/m_periodZ;
  double* dqVec = (double*) malloc(3 * sizeof(double));
  dqVec[0] = dqz; dqVec[1] = dqx; dqVec[2] = dqx;

  // initialize field in complex space:
  const bool real = false;

  // initialize data -
  // set of fourier peaks from Kai's code:
  typedef std::vector<int>    intPoint;
  typedef std::tuple<intPoint, double> point;
  std::vector<point> initVals;
  initVals.push_back(makePoint( {   0,   0,  0 },  m_avDensity));
  initVals.push_back(makePoint( {   4,   0,  0 }, -0.0167031));
  initVals.push_back(makePoint( { 124,   0,  0 }, -0.0167031));
  initVals.push_back(makePoint( {   3,   1,  0 }, -0.0152148));
  initVals.push_back(makePoint( {   4,   1,  0 },  0.1859610));
  initVals.push_back(makePoint( {   5,   1,  0 },  0.0155903));
  initVals.push_back(makePoint( { 123,   1,  0 },  0.0155903));
  initVals.push_back(makePoint( { 124,   1,  0 }, -0.1859610));
  initVals.push_back(makePoint( { 125,   1,  0 }, -0.0152148));
  initVals.push_back(makePoint( {   3,   2,  0 },  0.0133784));
  initVals.push_back(makePoint( { 125,   2,  0 }, -0.0133784));
  initVals.push_back(makePoint( {   1,   3,  0 }, -0.0152148));
  initVals.push_back(makePoint( {   2,   3,  0 },  0.0133784));
  initVals.push_back(makePoint( {   3,   3,  0 },  0.1906830));
  initVals.push_back(makePoint( {   5,   3,  0 },  0.0105081));
  initVals.push_back(makePoint( { 123,   3,  0 },  0.0105081));
  initVals.push_back(makePoint( { 125,   3,  0 },  0.1906830));
  initVals.push_back(makePoint( { 126,   3,  0 }, -0.0133784));
  initVals.push_back(makePoint( { 127,   3,  0 }, -0.0152148));
  initVals.push_back(makePoint( {   0,   4,  0 }, -0.0167031));
  initVals.push_back(makePoint( {   1,   4,  0 },  0.1859610));
  initVals.push_back(makePoint( { 127,   4,  0 }, -0.1859610));
  initVals.push_back(makePoint( {   1,   5,  0 },  0.0155903));
  initVals.push_back(makePoint( {   3,   5,  0 },  0.0105081));
  initVals.push_back(makePoint( {   5,   5,  0 },  0.0107489));
  initVals.push_back(makePoint( { 123,   5,  0 },  0.0107489));
  initVals.push_back(makePoint( { 125,   5,  0 },  0.0105081));
  initVals.push_back(makePoint( { 127,   5,  0 },  0.0155903));
  initVals.push_back(makePoint( {   1, 123,  0 },  0.0155903));
  initVals.push_back(makePoint( {   3, 123,  0 },  0.0105081));
  initVals.push_back(makePoint( {   5, 123,  0 },  0.0107489));
  initVals.push_back(makePoint( { 123, 123,  0 },  0.0107489));
  initVals.push_back(makePoint( { 125, 123,  0 },  0.0105081));
  initVals.push_back(makePoint( { 127, 123,  0 },  0.0155903));
  initVals.push_back(makePoint( {   0, 124,  0 }, -0.0167031));
  initVals.push_back(makePoint( {   1, 124,  0 }, -0.1859610));
  initVals.push_back(makePoint( { 127, 124,  0 },  0.1859610));
  initVals.push_back(makePoint( {   1, 125,  0 }, -0.0152148));
  initVals.push_back(makePoint( {   2, 125,  0 }, -0.0133784));
  initVals.push_back(makePoint( {   3, 125,  0 },  0.1906830));
  initVals.push_back(makePoint( {   5, 125,  0 },  0.0105081));
  initVals.push_back(makePoint( { 123, 125,  0 },  0.0105081));
  initVals.push_back(makePoint( { 125, 125,  0 },  0.1906830));
  initVals.push_back(makePoint( { 126, 125,  0 },  0.0133784));
  initVals.push_back(makePoint( { 127, 125,  0 }, -0.0152148));
  initVals.push_back(makePoint( {   3, 126,  0 }, -0.0133784));
  initVals.push_back(makePoint( { 125, 126,  0 },  0.0133784));
  initVals.push_back(makePoint( {   3, 127,  0 }, -0.0152148));
  initVals.push_back(makePoint( {   4, 127,  0 }, -0.1859610));
  initVals.push_back(makePoint( {   5, 127,  0 },  0.0155903));
  initVals.push_back(makePoint( { 123, 127,  0 },  0.0155903));
  initVals.push_back(makePoint( { 124, 127,  0 },  0.1859610));
  initVals.push_back(makePoint( { 125, 127,  0 }, -0.0152148));
  initVals.push_back(makePoint( {   3,   0,  1 }, -0.0127882));
  initVals.push_back(makePoint( { 125,   0,  1 }, -0.0127882));
  initVals.push_back(makePoint( {   3,   1,  1 }, -0.0340823));
  initVals.push_back(makePoint( {   4,   1,  1 },  0.1662350));
  initVals.push_back(makePoint( {   5,   1,  1 },  0.0224927));
  initVals.push_back(makePoint( { 123,   1,  1 }, -0.0224927));
  initVals.push_back(makePoint( { 124,   1,  1 },  0.1662350));
  initVals.push_back(makePoint( { 125,   1,  1 },  0.0340823));
  initVals.push_back(makePoint( {   2,   2,  1 }, -0.0142726));
  initVals.push_back(makePoint( {   3,   2,  1 }, -0.0145262));
  initVals.push_back(makePoint( {   4,   2,  1 },  0.0103089));
  initVals.push_back(makePoint( {   5,   2,  1 },  0.0102342));
  initVals.push_back(makePoint( { 123,   2,  1 },  0.0102342));
  initVals.push_back(makePoint( { 124,   2,  1 }, -0.0103089));
  initVals.push_back(makePoint( { 125,   2,  1 }, -0.0145262));
  initVals.push_back(makePoint( { 126,   2,  1 },  0.0142726));
  initVals.push_back(makePoint( {   0,   3,  1 }, -0.0127882));
  initVals.push_back(makePoint( {   1,   3,  1 }, -0.0340823));
  initVals.push_back(makePoint( {   2,   3,  1 }, -0.0145262));
  initVals.push_back(makePoint( {   3,   3,  1 }, -0.1606770));
  initVals.push_back(makePoint( {   4,   3,  1 },  0.0191507));
  initVals.push_back(makePoint( { 124,   3,  1 },  0.0191507));
  initVals.push_back(makePoint( { 125,   3,  1 },  0.1606770));
  initVals.push_back(makePoint( { 126,   3,  1 }, -0.0145262));
  initVals.push_back(makePoint( { 127,   3,  1 },  0.0340823));
  initVals.push_back(makePoint( {   1,   4,  1 },  0.1662350));
  initVals.push_back(makePoint( {   2,   4,  1 },  0.0103089));
  initVals.push_back(makePoint( {   3,   4,  1 },  0.0191507));
  initVals.push_back(makePoint( { 125,   4,  1 },  0.0191507));
  initVals.push_back(makePoint( { 126,   4,  1 }, -0.0103089));
  initVals.push_back(makePoint( { 127,   4,  1 },  0.1662350));
  initVals.push_back(makePoint( {   1,   5,  1 },  0.0224927));
  initVals.push_back(makePoint( {   2,   5,  1 },  0.0102342));
  initVals.push_back(makePoint( { 126,   5,  1 },  0.0102342));
  initVals.push_back(makePoint( { 127,   5,  1 }, -0.0224927));
  initVals.push_back(makePoint( {   1, 123,  1 }, -0.0224927));
  initVals.push_back(makePoint( {   2, 123,  1 },  0.0102342));
  initVals.push_back(makePoint( { 126, 123,  1 },  0.0102342));
  initVals.push_back(makePoint( { 127, 123,  1 },  0.0224927));
  initVals.push_back(makePoint( {   1, 124,  1 },  0.1662350));
  initVals.push_back(makePoint( {   2, 124,  1 }, -0.0103089));
  initVals.push_back(makePoint( {   3, 124,  1 },  0.0191507));
  initVals.push_back(makePoint( { 125, 124,  1 },  0.0191507));
  initVals.push_back(makePoint( { 126, 124,  1 },  0.0103089));
  initVals.push_back(makePoint( { 127, 124,  1 },  0.1662350));
  initVals.push_back(makePoint( {   0, 125,  1 }, -0.0127882));
  initVals.push_back(makePoint( {   1, 125,  1 },  0.0340823));
  initVals.push_back(makePoint( {   2, 125,  1 }, -0.0145262));
  initVals.push_back(makePoint( {   3, 125,  1 },  0.1606770));
  initVals.push_back(makePoint( {   4, 125,  1 },  0.0191507));
  initVals.push_back(makePoint( { 124, 125,  1 },  0.0191507));
  initVals.push_back(makePoint( { 125, 125,  1 }, -0.1606770));
  initVals.push_back(makePoint( { 126, 125,  1 }, -0.0145262));
  initVals.push_back(makePoint( { 127, 125,  1 }, -0.0340823));
  initVals.push_back(makePoint( {   2, 126,  1 },  0.0142726));
  initVals.push_back(makePoint( {   3, 126,  1 }, -0.0145262));
  initVals.push_back(makePoint( {   4, 126,  1 }, -0.0103089));
  initVals.push_back(makePoint( {   5, 126,  1 },  0.0102342));
  initVals.push_back(makePoint( { 123, 126,  1 },  0.0102342));
  initVals.push_back(makePoint( { 124, 126,  1 },  0.0103089));
  initVals.push_back(makePoint( { 125, 126,  1 }, -0.0145262));
  initVals.push_back(makePoint( { 126, 126,  1 }, -0.0142726));
  initVals.push_back(makePoint( {   3, 127,  1 },  0.0340823));
  initVals.push_back(makePoint( {   4, 127,  1 },  0.1662350));
  initVals.push_back(makePoint( {   5, 127,  1 }, -0.0224927));
  initVals.push_back(makePoint( { 123, 127,  1 },  0.0224927));
  initVals.push_back(makePoint( { 124, 127,  1 },  0.1662350));
  initVals.push_back(makePoint( { 125, 127,  1 }, -0.0340823));
  initVals.push_back(makePoint( {   0,   0,  2 },  0.1558300));
  initVals.push_back(makePoint( {   2,   0,  2 },  0.1327290));
  initVals.push_back(makePoint( { 126,   0,  2 },  0.1327290));
  initVals.push_back(makePoint( {   1,   1,  2 }, -0.0358995));
  initVals.push_back(makePoint( {   2,   1,  2 },  0.1359910));
  initVals.push_back(makePoint( {   3,   1,  2 },  0.0529617));
  initVals.push_back(makePoint( { 125,   1,  2 },  0.0529617));
  initVals.push_back(makePoint( { 126,   1,  2 }, -0.1359910));
  initVals.push_back(makePoint( { 127,   1,  2 }, -0.0358995));
  initVals.push_back(makePoint( {   0,   2,  2 },  0.1327290));
  initVals.push_back(makePoint( {   1,   2,  2 },  0.1359910));
  initVals.push_back(makePoint( {   2,   2,  2 }, -0.0743355));
  initVals.push_back(makePoint( {   3,   2,  2 }, -0.0204280));
  initVals.push_back(makePoint( { 125,   2,  2 },  0.0204280));
  initVals.push_back(makePoint( { 126,   2,  2 }, -0.0743355));
  initVals.push_back(makePoint( { 127,   2,  2 }, -0.1359910));
  initVals.push_back(makePoint( {   1,   3,  2 },  0.0529617));
  initVals.push_back(makePoint( {   2,   3,  2 }, -0.0204280));
  initVals.push_back(makePoint( {   5,   3,  2 },  0.0120951));
  initVals.push_back(makePoint( { 123,   3,  2 },  0.0120951));
  initVals.push_back(makePoint( { 126,   3,  2 },  0.0204280));
  initVals.push_back(makePoint( { 127,   3,  2 },  0.0529617));
  initVals.push_back(makePoint( {   3,   5,  2 },  0.0120951));
  initVals.push_back(makePoint( { 125,   5,  2 },  0.0120951));
  initVals.push_back(makePoint( {   3, 123,  2 },  0.0120951));
  initVals.push_back(makePoint( { 125, 123,  2 },  0.0120951));
  initVals.push_back(makePoint( {   1, 125,  2 },  0.0529617));
  initVals.push_back(makePoint( {   2, 125,  2 },  0.0204280));
  initVals.push_back(makePoint( {   5, 125,  2 },  0.0120951));
  initVals.push_back(makePoint( { 123, 125,  2 },  0.0120951));
  initVals.push_back(makePoint( { 126, 125,  2 }, -0.0204280));
  initVals.push_back(makePoint( { 127, 125,  2 },  0.0529617));
  initVals.push_back(makePoint( {   0, 126,  2 },  0.1327290));
  initVals.push_back(makePoint( {   1, 126,  2 }, -0.1359910));
  initVals.push_back(makePoint( {   2, 126,  2 }, -0.0743355));
  initVals.push_back(makePoint( {   3, 126,  2 },  0.0204280));
  initVals.push_back(makePoint( { 125, 126,  2 }, -0.0204280));
  initVals.push_back(makePoint( { 126, 126,  2 }, -0.0743355));
  initVals.push_back(makePoint( { 127, 126,  2 },  0.1359910));
  initVals.push_back(makePoint( {   1, 127,  2 }, -0.0358995));
  initVals.push_back(makePoint( {   2, 127,  2 }, -0.1359910));
  initVals.push_back(makePoint( {   3, 127,  2 },  0.0529617));
  initVals.push_back(makePoint( { 125, 127,  2 },  0.0529617));
  initVals.push_back(makePoint( { 126, 127,  2 },  0.1359910));
  initVals.push_back(makePoint( { 127, 127,  2 }, -0.0358995));
  initVals.push_back(makePoint( {   0,   0, 62 },  0.1558300));
  initVals.push_back(makePoint( {   2,   0, 62 },  0.1327290));
  initVals.push_back(makePoint( { 126,   0, 62 },  0.1327290));
  initVals.push_back(makePoint( {   1,   1, 62 }, -0.0358995));
  initVals.push_back(makePoint( {   2,   1, 62 },  0.1359910));
  initVals.push_back(makePoint( {   3,   1, 62 },  0.0529617));
  initVals.push_back(makePoint( { 125,   1, 62 },  0.0529617));
  initVals.push_back(makePoint( { 126,   1, 62 }, -0.1359910));
  initVals.push_back(makePoint( { 127,   1, 62 }, -0.0358995));
  initVals.push_back(makePoint( {   0,   2, 62 },  0.1327290));
  initVals.push_back(makePoint( {   1,   2, 62 },  0.1359910));
  initVals.push_back(makePoint( {   2,   2, 62 }, -0.0743355));
  initVals.push_back(makePoint( {   3,   2, 62 }, -0.0204280));
  initVals.push_back(makePoint( { 125,   2, 62 },  0.0204280));
  initVals.push_back(makePoint( { 126,   2, 62 }, -0.0743355));
  initVals.push_back(makePoint( { 127,   2, 62 }, -0.1359910));
  initVals.push_back(makePoint( {   1,   3, 62 },  0.0529617));
  initVals.push_back(makePoint( {   2,   3, 62 }, -0.0204280));
  initVals.push_back(makePoint( {   5,   3, 62 },  0.0120951));
  initVals.push_back(makePoint( { 123,   3, 62 },  0.0120951));
  initVals.push_back(makePoint( { 126,   3, 62 },  0.0204280));
  initVals.push_back(makePoint( { 127,   3, 62 },  0.0529617));
  initVals.push_back(makePoint( {   3,   5, 62 },  0.0120951));
  initVals.push_back(makePoint( { 125,   5, 62 },  0.0120951));
  initVals.push_back(makePoint( {   3, 123, 62 },  0.0120951));
  initVals.push_back(makePoint( { 125, 123, 62 },  0.0120951));
  initVals.push_back(makePoint( {   1, 125, 62 },  0.0529617));
  initVals.push_back(makePoint( {   2, 125, 62 },  0.0204280));
  initVals.push_back(makePoint( {   5, 125, 62 },  0.0120951));
  initVals.push_back(makePoint( { 123, 125, 62 },  0.0120951));
  initVals.push_back(makePoint( { 126, 125, 62 }, -0.0204280));
  initVals.push_back(makePoint( { 127, 125, 62 },  0.0529617));
  initVals.push_back(makePoint( {   0, 126, 62 },  0.1327290));
  initVals.push_back(makePoint( {   1, 126, 62 }, -0.1359910));
  initVals.push_back(makePoint( {   2, 126, 62 }, -0.0743355));
  initVals.push_back(makePoint( {   3, 126, 62 },  0.0204280));
  initVals.push_back(makePoint( { 125, 126, 62 }, -0.0204280));
  initVals.push_back(makePoint( { 126, 126, 62 }, -0.0743355));
  initVals.push_back(makePoint( { 127, 126, 62 },  0.1359910));
  initVals.push_back(makePoint( {   1, 127, 62 }, -0.0358995));
  initVals.push_back(makePoint( {   2, 127, 62 }, -0.1359910));
  initVals.push_back(makePoint( {   3, 127, 62 },  0.0529617));
  initVals.push_back(makePoint( { 125, 127, 62 },  0.0529617));
  initVals.push_back(makePoint( { 126, 127, 62 },  0.1359910));
  initVals.push_back(makePoint( { 127, 127, 62 }, -0.0358995));
  initVals.push_back(makePoint( {   3,   0, 63 }, -0.0127882));
  initVals.push_back(makePoint( { 125,   0, 63 }, -0.0127882));
  initVals.push_back(makePoint( {   3,   1, 63 }, -0.0340823));
  initVals.push_back(makePoint( {   4,   1, 63 },  0.1662350));
  initVals.push_back(makePoint( {   5,   1, 63 },  0.0224927));
  initVals.push_back(makePoint( { 123,   1, 63 }, -0.0224927));
  initVals.push_back(makePoint( { 124,   1, 63 },  0.1662350));
  initVals.push_back(makePoint( { 125,   1, 63 },  0.0340823));
  initVals.push_back(makePoint( {   2,   2, 63 }, -0.0142726));
  initVals.push_back(makePoint( {   3,   2, 63 }, -0.0145262));
  initVals.push_back(makePoint( {   4,   2, 63 },  0.0103089));
  initVals.push_back(makePoint( {   5,   2, 63 },  0.0102342));
  initVals.push_back(makePoint( { 123,   2, 63 },  0.0102342));
  initVals.push_back(makePoint( { 124,   2, 63 }, -0.0103089));
  initVals.push_back(makePoint( { 125,   2, 63 }, -0.0145262));
  initVals.push_back(makePoint( { 126,   2, 63 },  0.0142726));
  initVals.push_back(makePoint( {   0,   3, 63 }, -0.0127882));
  initVals.push_back(makePoint( {   1,   3, 63 }, -0.0340823));
  initVals.push_back(makePoint( {   2,   3, 63 }, -0.0145262));
  initVals.push_back(makePoint( {   3,   3, 63 }, -0.1606770));
  initVals.push_back(makePoint( {   4,   3, 63 },  0.0191507));
  initVals.push_back(makePoint( { 124,   3, 63 },  0.0191507));
  initVals.push_back(makePoint( { 125,   3, 63 },  0.1606770));
  initVals.push_back(makePoint( { 126,   3, 63 }, -0.0145262));
  initVals.push_back(makePoint( { 127,   3, 63 },  0.0340823));
  initVals.push_back(makePoint( {   1,   4, 63 },  0.1662350));
  initVals.push_back(makePoint( {   2,   4, 63 },  0.0103089));
  initVals.push_back(makePoint( {   3,   4, 63 },  0.0191507));
  initVals.push_back(makePoint( { 125,   4, 63 },  0.0191507));
  initVals.push_back(makePoint( { 126,   4, 63 }, -0.0103089));
  initVals.push_back(makePoint( { 127,   4, 63 },  0.1662350));
  initVals.push_back(makePoint( {   1,   5, 63 },  0.0224927));
  initVals.push_back(makePoint( {   2,   5, 63 },  0.0102342));
  initVals.push_back(makePoint( { 126,   5, 63 },  0.0102342));
  initVals.push_back(makePoint( { 127,   5, 63 }, -0.0224927));
  initVals.push_back(makePoint( {   1, 123, 63 }, -0.0224927));
  initVals.push_back(makePoint( {   2, 123, 63 },  0.0102342));
  initVals.push_back(makePoint( { 126, 123, 63 },  0.0102342));
  initVals.push_back(makePoint( { 127, 123, 63 },  0.0224927));
  initVals.push_back(makePoint( {   1, 124, 63 },  0.1662350));
  initVals.push_back(makePoint( {   2, 124, 63 }, -0.0103089));
  initVals.push_back(makePoint( {   3, 124, 63 },  0.0191507));
  initVals.push_back(makePoint( { 125, 124, 63 },  0.0191507));
  initVals.push_back(makePoint( { 126, 124, 63 },  0.0103089));
  initVals.push_back(makePoint( { 127, 124, 63 },  0.1662350));
  initVals.push_back(makePoint( {   0, 125, 63 }, -0.0127882));
  initVals.push_back(makePoint( {   1, 125, 63 },  0.0340823));
  initVals.push_back(makePoint( {   2, 125, 63 }, -0.0145262));
  initVals.push_back(makePoint( {   3, 125, 63 },  0.1606770));
  initVals.push_back(makePoint( {   4, 125, 63 },  0.0191507));
  initVals.push_back(makePoint( { 124, 125, 63 },  0.0191507));
  initVals.push_back(makePoint( { 125, 125, 63 }, -0.1606770));
  initVals.push_back(makePoint( { 126, 125, 63 }, -0.0145262));
  initVals.push_back(makePoint( { 127, 125, 63 }, -0.0340823));
  initVals.push_back(makePoint( {   2, 126, 63 },  0.0142726));
  initVals.push_back(makePoint( {   3, 126, 63 }, -0.0145262));
  initVals.push_back(makePoint( {   4, 126, 63 }, -0.0103089));
  initVals.push_back(makePoint( {   5, 126, 63 },  0.0102342));
  initVals.push_back(makePoint( { 123, 126, 63 },  0.0102342));
  initVals.push_back(makePoint( { 124, 126, 63 },  0.0103089));
  initVals.push_back(makePoint( { 125, 126, 63 }, -0.0145262));
  initVals.push_back(makePoint( { 126, 126, 63 }, -0.0142726));
  initVals.push_back(makePoint( {   3, 127, 63 },  0.0340823));
  initVals.push_back(makePoint( {   4, 127, 63 },  0.1662350));
  initVals.push_back(makePoint( {   5, 127, 63 }, -0.0224927));
  initVals.push_back(makePoint( { 123, 127, 63 },  0.0224927));
  initVals.push_back(makePoint( { 124, 127, 63 },  0.1662350));
  initVals.push_back(makePoint( { 125, 127, 63 }, -0.0340823));

  fftw_complex* data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numFieldElements);
  // initialize values to zero
  for (int index = 0; index < numFieldElements; index++) {
    data[index][0] = 0.0;
    data[index][1] = 0.0;
  }
  for (std::vector<point>::iterator it = initVals.begin(); it != initVals.end(); it++) {
    std::vector<int> k = std::get<0>(*it);	  
    double	     u = std::get<1>(*it);

    int kz = k[2];
    int ky = k[1];
    int kx = k[0];

    int index = kx + (N * ky) + (N * N * kz);

    data[index][0] = u * m_amplitude;
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