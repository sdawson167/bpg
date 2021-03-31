#include "A15PhaseProvider.h"

#define _USE_MATH_DEFINES
#include <math.h>

FieldProvider A15PhaseProvider::generateInitialCondition(int gridSize) {

  // initialize field provider variables:

  const int N = gridSize;

  int* gridSizes = (int*) malloc(3 * sizeof(int));
  gridSizes[0] = N; gridSizes[1] = N; gridSizes[2] = N;

  const int numFieldElements = N * N * N;

  // grid spacing based on period (box size)
  double* dqVec = (double*) malloc(3 * sizeof(double));

  // initialize field in complex space:
  const bool real = false;

  // initialize cplx data array:
  fftw_complex* data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numFieldElements);
  populateDataArray(dqVec, data, numFieldElements, gridSizes);

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

void A15PhaseProvider::resetCondition(FieldProvider &field) 
{
  // verify phaseID
  const int phaseID = field.getPhaseID();
  if (phaseID != m_phaseID)
	  throw std::runtime_error("Could not reset A15 phase - incorrect phase ID");

  // unpack field provider
  int* gridSizes         = field.getGridSizes();
  int  numFieldElements  = field.getNumFieldElements();
  double* dqVec 	 = field.getDq(); 
  fftw_complex* cplxData = field.getCplxDataPointer();
  
  // set values of cplxData
  populateDataArray(dqVec, cplxData, numFieldElements, gridSizes);
  field.updateDx();

  // update real data
  field.transformC2R();
}

void A15PhaseProvider::populateDataArray(double* dqVec, fftw_complex* cplxData, int numFieldElements, int* gridSizes)
{
  // reset grid spacing:
  double dq = 2 * M_PI / m_period;
  dqVec[0] = dq;
  dqVec[1] = dq;
  dqVec[2] = dq;

  // set all values in cplxData to zero
  for (int index = 0; index < numFieldElements; index++) {
    cplxData[index][0] = 0.0;
    cplxData[index][1] = 0.0;
  }

  // Fourier peaks describing A15 phase
  std::vector<point> initVals;
  initVals.push_back( makePoint( {  0,  0,  0 },    m_avDensity));          

  initVals.push_back( makePoint( { -2, -1,  0 },    0.0338277052131978 / 2)); 
  initVals.push_back( makePoint( {  2,  1,  0 },    0.0338277052131978 / 2)); 

  initVals.push_back( makePoint( {  2, -1,  0 },    0.0338277052131977 / 2)); 
  initVals.push_back( makePoint( { -2,  1,  0 },    0.0338277052131977 / 2)); 

  initVals.push_back( makePoint( { -1,  2,  0 },   -0.0338277052131977 / 2)); 
  initVals.push_back( makePoint( {  1, -2,  0 },   -0.0338277052131977 / 2)); 
  
  initVals.push_back( makePoint( {  1,  2,  0 },   -0.0338277052131977 / 2)); 
  initVals.push_back( makePoint( { -1, -2,  0 },   -0.0338277052131977 / 2)); 
  
  initVals.push_back( makePoint( {  0, -2,  1 },    0.0338228603092057 / 2)); 
  initVals.push_back( makePoint( {  0,  2, -1 },    0.0338228603092057 / 2)); // new
  
  initVals.push_back( makePoint( { -2,  0,  1 },   -0.0338228603092056 / 2)); 
  initVals.push_back( makePoint( {  2,  0, -1 },   -0.0338228603092056 / 2)); // new  
  
  initVals.push_back( makePoint( {  0,  2,  1 },    0.0338228603092055 / 2)); 
  initVals.push_back( makePoint( {  0, -2, -1 },    0.0338228603092055 / 2)); // new 
 
  initVals.push_back( makePoint( {  2,  0,  1 },   -0.0338228603092055 / 2)); 
  initVals.push_back( makePoint( { -2,  0, -1 },   -0.0338228603092055 / 2)); // new 

  initVals.push_back( makePoint( { -1,  0,  2 },    0.0338139694707322 / 2)); 
  initVals.push_back( makePoint( {  1,  0, -2 },    0.0338139694707322 / 2)); // new 
 
  initVals.push_back( makePoint( {  0, -1,  2 },   -0.0338139694707320 / 2)); 
  initVals.push_back( makePoint( {  0,  1, -2 },   -0.0338139694707320 / 2)); // new 

  initVals.push_back( makePoint( {  0,  1,  2 },   -0.0338139694707319 / 2)); 
  initVals.push_back( makePoint( {  0, -1, -2 },   -0.0338139694707319 / 2)); // new 
 
  initVals.push_back( makePoint( {  1,  0,  2 },    0.0338139694707319 / 2)); 
  initVals.push_back( makePoint( { -1,  0, -2 },    0.0338139694707319 / 2)); // new 
  
  initVals.push_back( makePoint( {  0,  0,  2 },    0.0331853928814813 / 2)); 
  initVals.push_back( makePoint( {  0,  0, -2 },    0.0331853928814813 / 2)); // new 
  
  initVals.push_back( makePoint( {  0,  2,  0 },    0.0331847214985317 / 2)); 
  initVals.push_back( makePoint( {  0, -2,  0 },    0.0331847214985317 / 2)); 
 
  initVals.push_back( makePoint( {  2,  0,  0 },    0.0331847214985317 / 2)); 
  initVals.push_back( makePoint( { -2,  0,  0 },    0.0331847214985317 / 2)); 
  
  initVals.push_back( makePoint( {  1, -1,  2 },    0.0297069055330425 / 2)); 
  initVals.push_back( makePoint( { -1,  1, -2 },    0.0297069055330425 / 2)); // new 
  
  initVals.push_back( makePoint( { -1, -1,  2 },    0.0297069055330425 / 2)); 
  initVals.push_back( makePoint( {  1,  1, -2 },    0.0297069055330425 / 2)); // new 
  
  initVals.push_back( makePoint( {  1,  1,  2 },    0.0297069055330424 / 2)); 
  initVals.push_back( makePoint( { -1, -1, -2 },    0.0297069055330424 / 2)); // new 
  
  initVals.push_back( makePoint( { -1,  1,  2 },    0.0297069055330424 / 2)); 
  initVals.push_back( makePoint( {  1, -1, -2 },    0.0297069055330424 / 2)); // new 
  
  initVals.push_back( makePoint( { -1,  2,  1 },    0.0297065456817734 / 2)); 
  initVals.push_back( makePoint( {  1, -2, -1 },    0.0297065456817734 / 2)); // new 
  
  initVals.push_back( makePoint( { -2, -1,  1 },    0.0297065456817733 / 2)); 
  initVals.push_back( makePoint( {  2,  1, -1 },    0.0297065456817733 / 2)); // new 
  
  initVals.push_back( makePoint( { -1, -2,  1 },    0.0297065456817733 / 2)); 
  initVals.push_back( makePoint( {  1,  2, -1 },    0.0297065456817733 / 2)); // new 
  
  initVals.push_back( makePoint( {  2, -1,  1 },    0.0297065456817733 / 2)); 
  initVals.push_back( makePoint( { -2,  1, -1 },    0.0297065456817733 / 2)); // new 
  
  initVals.push_back( makePoint( { -2,  1,  1 },    0.0297065456817732 / 2)); 
  initVals.push_back( makePoint( {  2, -1, -1 },    0.0297065456817732 / 2)); // new 
  
  initVals.push_back( makePoint( {  2,  1,  1 },    0.0297065456817732 / 2)); 
  initVals.push_back( makePoint( { -2, -1, -1 },    0.0297065456817732 / 2)); // new 
  
  initVals.push_back( makePoint( {  1,  2,  1 },    0.0297065456817732 / 2)); 
  initVals.push_back( makePoint( { -1, -2, -1 },    0.0297065456817732 / 2)); // new 
  
  initVals.push_back( makePoint( {  1, -2,  1 },    0.0297065456817731 / 2)); 
  initVals.push_back( makePoint( { -1,  2, -1 },    0.0297065456817731 / 2)); // new 

  // load non-zero peaks into cplx data array
  for (std::vector<point>::iterator it = initVals.begin(); it != initVals.end(); it++) {
    std::vector<int> k = std::get<0>(*it);                                              
    double	     u = std::get<1>(*it);
                                                                                         
    int kz = k[2] < 0 ? k[2] + gridSizes[0] : k[2];
    int ky = k[1] < 0 ? k[1] + gridSizes[1] : k[1];
    int kx = k[0] < 0 ? k[0] + gridSizes[2] : k[0];
                                                                                         
    int index = kx + (gridSizes[2] * ky) + (gridSizes[2] * gridSizes[1] * kz);
                                                                                         
    cplxData[index][0] = u * m_amplitude;
  }
}
