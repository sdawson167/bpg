#ifndef _FCC_GUARD
#define _FCC_GUARD

#include "Phase.h"
#include "FieldProvider.h"

/*
 * class object facilitates creation of field provider in face-centred cubic phase
 */

class FccPhaseProvider {
private:
    double m_period;
    double m_avDensity;
    double m_amplitude;
    const int m_dimension = 3;
    const int m_phaseDimension = 1;

    typedef std::vector<int>    intPoint;
    typedef std::tuple<intPoint, double> point;
    point makePoint(intPoint coords, double amp) { return point(coords, amp); }

public:
    FccPhaseProvider(double period, double avDensity, double amplitude)
      : m_period{ period },
        m_avDensity{ avDensity },
        m_amplitude{ amplitude }
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
      return Phase::fcc;
    }

    /*
     * ======================================
     *      initialize field provider
     * ======================================
     */
    FieldProvider generateInitialCondition(int gridSize);

};
#endif
