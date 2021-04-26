//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Fulvio Piacenza: fulvio.piacenza01@universitadipavia.it
//

#include "apfel/twobodyphasespace.h"

#include <math.h>
#include <algorithm>
#include <apfel/integrator.h>

namespace apfel
{
  //_________________________________________________________________________
  TwoBodyPhaseSpace::TwoBodyPhaseSpace(double const& pTmin1, double const& pTmin2, double const& etamin, double const& etamax, double const& eps):
    _pTmin1(pTmin1),
    _pTmin2(pTmin2),
    _etamin(etamin),
    _etamax(etamax),
    _eps(eps)
  {
  }

  //_________________________________________________________________________
  TwoBodyPhaseSpace::TwoBodyPhaseSpace(double const& pTmin, double const& etamin, double const& etamax, double const& eps):
    TwoBodyPhaseSpace{pTmin, pTmin, etamin, etamax, eps}
  {
  }

  //_________________________________________________________________________
  double TwoBodyPhaseSpace::PhaseSpaceReduction(double const& Q, double const& y, double const& qT)
  {
    // Return automatically zero if "y" is larger that "_etamax"
    if (y <= _etamin || y >= _etamax)
      return 0;

    // Useful definitions
    const double Q2      = Q * Q;
    const double qT2     = qT * qT;
    const double qT4     = qT2 * qT2;
    const double pTmin22 = _pTmin2 * _pTmin2;
    const double M2      = Q2 + qT2;
    const double M       = sqrt(M2);
    const double ctghmax = 1 / tanh(y - _etamax);
    const double ctghmin = 1 / tanh(y - _etamin);

    // Integrand function
    const auto IntegrandP = [&] (double const& eta) -> double
    {
      // More useful definitions
      const double ch    = cosh(eta - y);
      const double sh2   = ch * ch - 1;
      const double Eq    = M * ch;
      const double Eq2   = Eq * Eq;
      const double Eq4   = Eq2 * Eq2;
      const double EmqT2 = Eq2 - qT2;
      const double EmqT4 = EmqT2 * EmqT2;
      const double EmqT  = sqrt(EmqT2);

      // Auxiliary functions
      const double f2    = ( 2 * _pTmin1 * Eq - Q2 ) / 2 / _pTmin1 / qT;
      const double f3max = Eq / qT - Q2 * ( sinh(eta - y) * ctghmax + ch ) / 2 / qT / M;
      const double f3min = Eq / qT - Q2 * ( sinh(eta - y) * ctghmin + ch ) / 2 / qT / M;
      const double f4    = ( Eq * ( Q2 - 2 * pTmin22 + 2 * qT2 ) - Q2 * sqrt( Eq2 - M2 + pTmin22 ) ) / 2 / qT / ( M2 - pTmin22 );

      // Integration limits in cos(phi)
      const double x1 = std::max(f2, -1.);
      const double x2 = std::min(std::min(f4, std::min(f3max, f3min)), 1.);

      // Primitive function of the intregration in cos(phi) (up to a
      // factor that can be compute externally) when contractin the
      // leptonic tensor with g^{\mu\nu}_\perp.
      const auto Fbar = [&] (double const& x) -> double
      {
        const double xs = x * x;
        const double xp = sqrt( 1 - xs );
        const double atanfact = Eq * ( atan( ( qT - x * Eq ) / EmqT / xp ) - atan( ( qT + x * Eq ) / EmqT / xp ) ) / EmqT;
        const double Fi = qT2 * x * xp / ( x2 * qT2 - Eq2 ) - atanfact;
        const double Gi = Q2 * sh2
        * ( qT2 * xp * ( ( ( 11 * Eq2 * qT2 + 4 * qT4 ) * x2 + 3 * Eq * qT * ( 9 * Eq2 + qT2 ) * x + 18 * Eq4 - 5 * Eq2 * qT2 + 2 * qT4 ) / pow(x * qT + Eq, 3) +
                         ( ( 11 * Eq2 * qT2 + 4 * qT4 ) * x2 - 3 * Eq * qT * ( 9 * Eq2 + qT2 ) * x + 18 * Eq4 - 5 * Eq2 * qT2 + 2 * qT4 ) / pow(x * qT - Eq, 3) )
            - 6 * ( 2 * Eq2 + 3 * qT2 ) * atanfact ) / EmqT4 / 4;

        return 3 * Fi + Gi;
      };

      // Construct integrand
      if (x2 > x1)
        return ( Fbar(x2) - Fbar(x1) ) / EmqT2;
      else
	return 0;
    };

    // Return integral (symmetrise with respect to the center of the
    // rapidity range to make sure that the result is symmetric). This
    // "misbehaviour" is due to the fact that the integrand function
    // is pieceswise and thus the numerical integrator struggles.
    const apfel::Integrator Ieta{IntegrandP};
    return Q2 * ( y > ( _etamin + _etamax ) / 2 ? Ieta.integrate(_etamin, _etamax, _eps) : - Ieta.integrate(_etamax, _etamin, _eps) ) / 16 / M_PI;
  }

  //_________________________________________________________________________
  double TwoBodyPhaseSpace::DerivePhaseSpaceReduction(double const& Q, double const& y, double const& qT)
  {
    const double eps = 1e-3;
    return ( PhaseSpaceReduction(Q, y, qT * ( 1 + eps )) - PhaseSpaceReduction(Q, y, qT * ( 1 - eps )) ) / 2 / eps / qT;
  }

  //_________________________________________________________________________
  double TwoBodyPhaseSpace::ParityViolatingPhaseSpaceReduction(double const& Q, double const& y, double const& qT)
  {
    // Return automatically zero if "y" is larger that "_etamax"
    if (y <= _etamin || y >= _etamax)
      return 0;

    // Useful definitions
    const double Q2      = Q * Q;
    const double qT2     = qT * qT;
    const double pTmin22 = _pTmin2 * _pTmin2;
    const double M2      = Q2 + qT2;
    const double M       = sqrt(M2);
    const double ctghmax = 1 / tanh(y - _etamax);
    const double ctghmin = 1 / tanh(y - _etamin);

    // Integrand function
    const auto IntegrandP = [&] (double const& eta) -> double
    {
      // More useful definitions
      const double ch    = cosh(eta - y);
      const double sh    = sinh(y - eta);
      const double Eq    = M * ch;
      const double Eq2   = Eq * Eq;
      const double EmqT2 = Eq2 - qT2;
      const double EmqT  = sqrt(EmqT2);

      // Auxiliary functions
      const double f2    = ( 2 * _pTmin1 * Eq - Q2 ) / 2 / _pTmin1 / qT;
      const double f3max = Eq / qT - Q2 * ( sinh(eta - y) * ctghmax + ch ) / 2 / qT / M;
      const double f3min = Eq / qT - Q2 * ( sinh(eta - y) * ctghmin + ch ) / 2 / qT / M;
      const double f4    = ( Eq * ( Q2 - 2 * pTmin22 + 2 * qT2 ) - Q2 * sqrt( Eq2 - M2 + pTmin22 ) ) / 2 / qT / ( M2 - pTmin22 );

      // Integration limits in cos(phi)
      const double x1 = std::max(f2, -1.);
      const double x2 = std::min(std::min(f4, std::min(f3max, f3min)), 1.);

      // Primitive function of the intregration in cos(phi) (up to a
      // factor that can be compute externally) when contractin the
      // leptonic tensor with g^{\mu\nu}_\perp.
      const auto Hbar = [&] (double const& x) -> double
      {
        const double xs = x * x;
        const double xp = sqrt( 1 - xs );
        const double Hi = sh
        * ( qT2 * xp * ( ( 3 * Eq * qT * x - ( 4 * Eq2 - qT2 ) ) / pow(x * qT - Eq, 2) +
                         ( 3 * Eq * qT * x + ( 4 * Eq2 - qT2 ) ) / pow(x * qT + Eq, 2) )
            - ( 2 * Eq2 + qT2 )
            * ( atan( ( qT - x * Eq ) / EmqT / xp ) - atan( ( qT + x * Eq ) / EmqT / xp ) ) / EmqT ) / EmqT2 / EmqT2;
        return Hi;
      };

      // Construct integrand
      if (x2 > x1)
        return ( Hbar(x2) - Hbar(x1) ) / EmqT2;
      else
	return 0;
    };

    // Return integral (symmetrise with respect to the center of the
    // rapidity range to make sure that the result is symmetric). This
    // "misbehaviour" is due to the fact that the integrand function
    // is pieceswise and thus the numerical integrator struggles.
    const apfel::Integrator Ieta{IntegrandP};
    return 3 * Q2 * Q * ( y > ( _etamin + _etamax ) / 2 ? Ieta.integrate(_etamin, _etamax, _eps) : - Ieta.integrate(_etamax, _etamin, _eps) ) / 128 / M_PI;
  }
}
