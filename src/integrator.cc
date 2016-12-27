//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/integrator.h"
#include "apfel/tools.h"
#include <cmath>

namespace apfel
{
  //_________________________________________________________________________
  Integrator::Integrator():
    _cst(1E-25)
  {
    // prepare the dgauss contants.
    _w = {
               0.1012285362903762591525313543e0,
               0.2223810344533744705443559944e0,
               0.3137066458778872873379622020e0,
               0.3626837833783619829651504493e0,
               0.0271524594117540948517805725e0,
               0.0622535239386478928628438370e0,
               0.0951585116824927848099251076e0,
               0.1246289712555338720524762822e0,
               0.1495959888165767320815017305e0,
               0.1691565193950025381893120790e0,
               0.1826034150449235888667636680e0,
               0.1894506104550684962853967232e0
    };
    _x = {
               0.9602898564975362316835608686e0,
               0.7966664774136267395915539365e0,
               0.5255324099163289858177390492e0,
               0.1834346424956498049394761424e0,
               0.9894009349916499325961541735e0,
               0.9445750230732325760779884155e0,
               0.8656312023878317438804678977e0,
               0.7554044083550030338951011948e0,
               0.6178762444026437484466717640e0,
               0.4580167776572273863424194430e0,
               0.2816035507792589132304605015e0,
               0.0950125098376374401853193354e0
    };
  }

  //_________________________________________________________________________
  double Integrator::integrate(const double &xmin, const double &xmax, const double &eps) const
  {    
    return dgauss(xmin, xmax, eps);
  }

  //_________________________________________________________________________
  double Integrator::dgauss(const double &a, const double &b, const double &eps) const
  {
    int i;
    double aa, bb, y, c1, c2, s8, s16, u;

    const double delta = _cst*fabs(a-b);
    double dgauss = 0;
    aa = a;

    goto5:
    y = b -aa;
    if (fabs(y) <= delta) return dgauss;

    goto2:
    bb = aa+y;
    c1 = 0.5*(aa+bb);
    c2 = c1-aa;
    s8 =0;
    s16 = 0;

    for (i = 1; i <= 4; i++)
      {
        u = _x[i-1]*c2;
        s8 += _w[i-1]*(integrand(c1+u)+integrand(c1-u));
      }

    for (i = 5; i <= 12 ; i++)
      {
        u = _x[i-1]*c2;
        s16 += _w[i-1]*(integrand(c1+u)+integrand(c1-u));
      }

    s8 *= c2;
    s16 *= c2;
    if (fabs(s16-s8) > eps*(1.0+fabs(s16))) goto goto4;
    dgauss += s16;
    aa = bb;
    goto goto5;

    goto4:
    y = 0.5*y;
    if (fabs(y) > delta) goto goto2;

    throw runtime_exception("Integrator::dgauss", "too high accuracy required");

    return dgauss;
  }
}
