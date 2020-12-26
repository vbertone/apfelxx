//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/integrator.h"
#include "apfel/constants.h"
#include "apfel/messages.h"

#include <array>
#include <algorithm>

namespace apfel
{
  /**
    * Parameters of the Gauss-Legendre integration with 8 and 16-point integration.
    */
  const std::array<std::vector<double>, 2> gl_x =
  {
    {
      std::vector<double>{0.1834346424956498, 0.5255324099163289, 0.7966664774136267, 0.9602898564975362},
      std::vector<double>{
        0.0950125098376374, 0.2816035507792589, 0.4580167776572273, 0.6178762444026437,
        0.7554044083550030, 0.8656312023878317, 0.9445750230732325, 0.9894009349916499
      }
    }
  };
  const std::array<std::vector<double>, 2> gl_w =
  {
    {
      std::vector<double>{0.3626837833783619, 0.3137066458778872, 0.2223810344533744, 0.1012285362903762},
      std::vector<double>{
        0.1894506104550684, 0.1826034150449235, 0.1691565193950025, 0.1495959888165767,
        0.1246289712555338, 0.0951585116824927, 0.0622535239386478, 0.0271524594117540
      }
    }
  };

  /**
    * Parameters of the Gauss-Kronrod integration with 7 and 15-point integration.
    */
  const std::array<std::vector<double>, 2> gk_x =
  {
    {
      std::vector<double>{0.0000000000000000e+00, 4.0584515137739717e-01, 7.4153118559939444e-01, 9.4910791234275852e-01},
      std::vector<double>{
        0.0000000000000000e+00, 2.0778495500789847e-01, 4.0584515137739717e-01, 5.8608723546769113e-01,
        7.4153118559939444e-01, 8.6486442335976907e-01, 9.4910791234275852e-01, 9.9145537112081264e-01
      }
    }
  };
  const std::array<std::vector<double>, 2> gk_w =
  {
    {
      std::vector<double>{4.1795918367346939e-01, 3.8183005050511894e-01, 2.7970539148927667e-01, 1.2948496616886969e-01},
      std::vector<double>{
        2.0948214108472783e-01, 2.0443294007529889e-01, 1.9035057806478541e-01, 1.6900472663926790e-01,
        1.4065325971552592e-01, 1.0479001032225018e-01, 6.3092092629978553e-02, 2.2935322010529225e-02
      }
    }
  };

  //_________________________________________________________________________
  Integrator::Integrator(std::function<double(double const&)> const& func, IntegrationMethod const& method):
    _func(func),
    _method(method)
  {
  }

  //_________________________________________________________________________
  double Integrator::integrate(double const& a, double const& b, double const& eps) const
  {
    switch (_method)
      {
      case GAUSS_LEGENDRE:
        return integrateGL(a, b, eps);
      case GAUSS_KRONROD:
        return integrateGK(a, b, eps);
      default:
        throw std::runtime_error(error("Integrator::integrate", "Unknown integration method"));
      }
  }

  //_________________________________________________________________________
  double Integrator::integrateGL(double const& a, double const& b, double const& eps) const
  {
    const double delta = eps25 * std::abs(a - b);
    double dgauss = 0;
    double aa = a;

goto5:
    double y = b - aa;
    if (std::abs(y) <= delta)
      return dgauss;

goto2:
    const double bb = aa + y;
    const double c1 = 0.5 * ( aa + bb );
    const double c2 = c1 - aa;

    double s8 = 0;
    for (int i = 0; i < 4; i++)
      {
        const double u = gl_x[0][i] * c2;
        s8 += gl_w[0][i] * ( integrand(c1+u) + integrand(c1-u) );
      }

    double s16 = 0;
    for (int i = 0; i < 8; i++)
      {
        const double u = gl_x[1][i] * c2;
        s16 += gl_w[1][i] * ( integrand(c1+u) + integrand(c1-u) );
      }

    s8  *= c2;
    s16 *= c2;
    if (std::abs(s16 - s8) > eps * ( 1 + std::abs(s16) ))
      goto goto4;

    dgauss += s16;
    aa = bb;
    goto goto5;

goto4:
    y = 0.5 * y;
    if (std::abs(y) > delta)
      goto goto2;

    throw std::runtime_error(error("Integrator::integrateGL", "Too high accuracy required."));
  }

  //_________________________________________________________________________
  double Integrator::integrateGK(double const& a, double const& b, double const& eps) const
  {
    const double delta = eps25 * std::abs(a - b);
    double dgauss = 0;
    double aa = a;

goto5:
    double y = b - aa;
    if (std::abs(y) <= delta)
      return dgauss;

goto2:
    double bb = aa + y;
    double c1 = 0.5 * ( aa + bb );
    double c2 = c1 - aa;

    const double f = integrand(c1);
    double s7  = gk_w[0][0] * f;
    double s15 = gk_w[1][0] * f;
    for (int i = 1; i < 4; i++)
      {
        const double u  = gk_x[0][i] * c2;
        const double fu = integrand(c1 + u);
        const double fd = integrand(c1 - u);
        s7 += gk_w[0][i] * ( fu + fd );

        const double u1 = gk_x[1][2 * i - 1] * c2;
        s15 += gk_w[1][2 * i - 1] * ( integrand(c1 + u1) + integrand(c1 - u1) ) + gk_w[1][2 * i] * ( fu + fd );
      }
    double w = gk_x[1][7] * c2;
    s15 += gk_w[1][7] * ( integrand(c1 + w) + integrand(c1 - w) );

    s7  *= c2;
    s15 *= c2;
    if (std::abs(s15 - s7) > eps * ( 1 + std::abs(s15) ))
      goto goto4;

    dgauss += s15;
    aa = bb;
    goto goto5;

goto4:
    y = 0.5 * y;
    if (std::abs(y) > delta)
      goto goto2;

    throw std::runtime_error(error("Integrator::integrateGK", "Too high accuracy required."));
  }

  //_________________________________________________________________________
  double Integrator::integrate(double const& a, double const& b, std::vector<double> const& FixPts, double const& eps) const
  {
    // Create vector of fixed points including the integration bounds
    // in ordered sequence and removing the points outside the
    // integration range.
    const double ap = (a < b ? a : b);
    const double bp = (a < b ? b : a);
    std::vector<double> bounds{ap, bp};
    for (auto const& ib : FixPts)
      if(ib > ap && ib < bp)
        bounds.push_back(ib);

    // Sort vector according on how the integration bounds are
    // ordered.
    std::sort(bounds.begin(), bounds.end(), [a, b] (double const& x1, double const& x2) -> bool { return (b > a ? x1 < x2 : x1 > x2); });

    // Now compute sub-integrals and sum them up
    double dgauss = 0;
    for (int i = 1; i < (int) bounds.size(); i++)
      dgauss += integrate(bounds[i-1], bounds[i], eps);

    return dgauss;
  }
}
