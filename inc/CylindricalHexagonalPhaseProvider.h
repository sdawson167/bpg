#ifndef _HEX_GUARD
#define _HEX_GUARD

#include "Phase.h"
#include "FieldProvider.h"

/*
 * class object facilitates creation of field provider in cylindrical hexagonal phase
 */

class CylindricalHexagonalPhaseProvider {
private:
    double m_period;
    double m_avDensity;
    double m_amplitude;
    const int m_dimension = 2;
    const int m_phaseDimension = 2;

    typedef std::vector<int>    intPoint;                                      	
    typedef std::tuple<intPoint, double> point;
    point makePoint(intPoint coords, double amp) { return point(coords, amp); }

public:
    CylindricalHexagonalPhaseProvider(double period, double avDensity, double amplitude)
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
      return Phase::hex;
    }

    /*
     * ======================================
     *      initialize field provider
     * ======================================
     */
    FieldProvider generateInitialCondition(int gridSize);
};
#endif