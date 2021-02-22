#include "GenericProvider.h"

#include "LamellarPhaseProvider.h"
#include "GyroidPhaseProvider.h"
#include "CylindricalHexagonalPhaseProvider.h"
#include "BccPhaseProvider.h"
#include "FccPhaseProvider.h"
#include "A15PhaseProvider.h"
#include "SigmaPhaseProvider.h"


FieldProvider* GenericPhaseProvider::generateInitialCondition()
{
  switch(m_phaseID) {
    case 1: {   // lam
      LamellarPhaseProvider lamProvider{ m_lamPeriod, m_avgDensity, m_amplitude };
      FieldProvider* initialFieldProvider = new FieldProvider(lamProvider.generateInitialCondition(64));
      return initialFieldProvider;
    }
    case 2: {   // gyr
      GyroidPhaseProvider gyrProvider{ m_gyrPeriod, m_avgDensity, m_amplitude };
      FieldProvider* initialFieldProvider = new FieldProvider(gyrProvider.generateInitialCondition(64));
      return initialFieldProvider;
    }
    case 3: {   // hex
      CylindricalHexagonalPhaseProvider hexProvider{ m_hexPeriod, m_avgDensity, m_amplitude };
      FieldProvider* initialFieldProvider = new FieldProvider(hexProvider.generateInitialCondition(64));
      return initialFieldProvider;
    }
    case 4: {   // bcc
      BccPhaseProvider bccProvider{ m_bccPeriod, m_avgDensity, m_amplitude};
      FieldProvider* initialFieldProvider = new FieldProvider(bccProvider.generateInitialCondition(64));
      return initialFieldProvider;
    }
    case 5: {   // fcc
      FccPhaseProvider fccProvider{ m_fccPeriod, m_avgDensity, m_amplitude};
      FieldProvider* initialFieldProvider = new FieldProvider(fccProvider.generateInitialCondition(64));
      return initialFieldProvider;
    }
    case 6: {   // a15
      A15PhaseProvider a15Provider{ m_a15Period, m_avgDensity, m_amplitude};
      FieldProvider* initialFieldProvider = new FieldProvider(a15Provider.generateInitialCondition(64));
      return initialFieldProvider;
    }
    case 7: {
      SigmaPhaseProvider sigProvider{ m_sigPeriodX, m_sigPeriodZ, m_avgDensity, m_amplitude};
      FieldProvider* initialFieldProvider = new FieldProvider(sigProvider.generateInitialCondition(64));
      return initialFieldProvider;
    }
    default: {  // dis
      int* gridSize = (int*) malloc(sizeof(int));
      gridSize[0] = 1;      

      double* dr = (double*) malloc(sizeof(double));
      dr[0] = 0.0;

      fftw_complex* data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex));
      FieldProvider* initialFieldProvider = 
      new FieldProvider{
        data, 
        1,
        gridSize,
        dr,
        false,
        0};
      fftw_free(data);
      
      return initialFieldProvider;
    }
  } // end switch
}
