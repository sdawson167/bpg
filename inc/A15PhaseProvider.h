#ifndef _A15_GUARD
#define _A15_GUARD

#include "Phase.h"
#include "FieldProvider.h"

/*
 * class object facilitates creation of field provider in Frank-Kasper A15 phase
 */

class A15PhaseProvider {
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
    A15PhaseProvider(double period, double avDensity, double amplitude)
      : m_period{ period },
        m_avDensity{ avDensity },
        m_amplitude{ amplitude }
    { };

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
      return Phase::a15;
    }

    /*
     * ======================================
     *      initialize field provider
     * ======================================
     */
    FieldProvider generateInitialCondition(int gridSize);

};
#endif
