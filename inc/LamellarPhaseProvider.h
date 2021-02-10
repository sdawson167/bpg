#ifndef _LAMELLAR_GUARD
#define _LAMELLAR_GUARD

#include "Phase.h"
#include "FieldProvider.h"

/*
 * class object faciliates creation of field provider in lamellar phase
 */

class LamellarPhaseProvider {
private:
    double m_period;
    double m_avDensity;
    double m_amplitude;

    const int m_dimension = 1;
    const int m_phaseDimension = 1;

public:
    LamellarPhaseProvider(double period, double avDensity, double amplitude)
    : m_period{period},
      m_avDensity{avDensity},
      m_amplitude{amplitude}
    {};

    /*
     * ======================================
     *          getters and setters
     * ======================================
     */
    double getPeriod(){ return m_period; }

    double getAverageDensity() { return m_avDensity; }

    double getAmplitude() { return m_amplitude; }

    Phase getPhase()
    {
      return Phase::lam;
    }

    /*
     * ======================================
     *      initialize field provider
     * ======================================
     */
    FieldProvider generateInitialCondition(int gridSize);

};
#endif
