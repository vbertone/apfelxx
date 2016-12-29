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
  Integrator::Integrator()
  {
  }

  //_________________________________________________________________________
  double Integrator::integrate(double const& a, double const& b, double const& eps) const
  {
    const double delta = cst * fabs(a-b);
    double dgauss = 0;
    double aa = a;

    goto5:
    double y = b - aa;
    if (fabs(y) <= delta) return dgauss;

    goto2:
    double bb  = aa + y;
    double c1  = 0.5 * ( aa + bb );
    double c2  = c1 - aa;

    double s8  = 0;
    for (auto i = 0; i < 4; i++)
      {
        double u = gq_x[1][i] * c2;
        s8 += gq_w[1][i] * ( integrand(c1+u) + integrand(c1-u) );
      }

    double s16 = 0;
    for (auto i = 0; i < 8; i++)
      {
        double u = gq_x[2][i] * c2;
        s16 += gq_w[2][i] * ( integrand(c1+u) + integrand(c1-u) );
      }

    s8  *= c2;
    s16 *= c2;
    if (fabs(s16-s8) > eps*(1+fabs(s16))) goto goto4;

    dgauss += s16;
    aa = bb;
    goto goto5;

    goto4:
    y = 0.5 * y;
    if (fabs(y) > delta) goto goto2;

    throw runtime_exception("Integrator::dgauss", "too high accuracy required");

    return dgauss;
  }

  //_________________________________________________________________________
  double Integrator::integrate(double const& a, double const& b, int const& m) const
  {
    double I  = 0;
    double k1 = 0.5 * ( b - a );
    double k2 = k1 + a;
    for(auto i = 0; i < (int)gq_x[m].size(); i++) 
      {
	double delta = gq_x[m][i] * k1;
	I += gq_w[m][i] * ( integrand( k2 + delta ) + integrand( k2 - delta ) );
      }
    I *= k1;
    return I;
  }

}
