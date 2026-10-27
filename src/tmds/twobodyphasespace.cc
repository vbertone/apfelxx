//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Fulvio Piacenza: fulvio.piacenza01@universitadipavia.it
//

#include "apfel/twobodyphasespace.h"

#include <math.h>
#include <algorithm>
#include <iostream>
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
  double TwoBodyPhaseSpace::ZeroqTPhaseSpaceReduction(double Q, double y) const
  {
    // Note: no support at the moment for asymemtric cuts in _pTmin1, _pTmin2
    if (Q <= 2. * _pTmin1) return 0.;

    const double delta = std::min(y - _etamin, _etamax - y);
    if (delta <= 0.) return 0.;

    const double rho = 2. * _pTmin1 / Q;
    const double t   = std::min(std::tanh(delta), std::sqrt(1. - rho*rho));
    return t * (3. + t * t) / 4.;
  }

  //_________________________________________________________________________
  double TwoBodyPhaseSpace::_integrand(double eta, double qT, double y, double MT) const
  {
    const double c  = std::cosh(y - eta);
    auto [x_l, x_u] = _get_x_bounds(qT, MT, c, y, eta);
    if (x_l >= x_u) return 0.;
    return _primitive(x_u, qT, MT, c) - _primitive(x_l, qT, MT, c);
  }

  //_________________________________________________________________________
  std::pair<double, double> TwoBodyPhaseSpace::_get_x_bounds(double qT, double MT, double c, double y,
                                                             double eta) const
  {
    const double f1 = (MT * c - (MT * MT - qT * qT) / (2. * _pTmin1)) / qT;
    const double xbar = [&]() -> double
    {
      const double den = 2. * qT * (MT * MT - _pTmin1 * _pTmin1);
      const double num = c * MT * (MT * MT + qT * qT - 2 * _pTmin1 * _pTmin1) +
      (qT * qT - MT * MT) * std::sqrt(_pTmin1 * _pTmin1 + MT * MT * (c * c - 1));
      return num / den;
    }();

    auto Rinv = [&](double eta_bar)
    {
      const double den = 2. * qT * MT * std::sinh(y - eta_bar);
      const double num =
        MT * MT * std::sinh(2 * y - eta - eta_bar) + qT * qT * std::sinh(eta - eta_bar);
      return num / den;
    };

    return {std::max(-1., f1),
            std::min({1., xbar, Rinv(_etamin), Rinv(_etamax)})};
  }

  //_________________________________________________________________________
  double TwoBodyPhaseSpace::_primitive(double x, double qT, double MT, double c) const
  {
    // Without 3. * sq(Q) / (8. * M_PI) prefactor
    const double qT2 = qT*qT;
    const double u = MT * c;
    const double u2 = u*u;
    const double d = (1 - c * c) * (qT2 - MT * MT);
    const double z = qT * x - u;
    const double z2 = z*z;
    const double k = u * u - qT2;
    const double k2 = k*k;

    const double sqrt_1_m_x2 = std::sqrt((1. - x) * (1. + x));
    const double sqrt_pref   = (d / z2 - 5. * u * d / (2. * k * z) + 2. * qT2 * d / k2 +
                                11 * d * u2 / (2. * k2) + 3) /
                               (3. * k * z);

    const double arctan_pref =
      u * (2. * k2 + d * (3. * qT2 + 2. * u2)) / std::pow(u2 - qT2, 3.5);
    const double arctan = -std::atan2(  //
                            std::sqrt((u + qT) * (1. - x)), //
                            std::sqrt((u - qT) * (1. + x))  //
                          );
    return qT * sqrt_1_m_x2 * sqrt_pref + arctan * arctan_pref;
  }

  //_________________________________________________________________________
  std::vector<double> TwoBodyPhaseSpace::eta_breakpoints(double Q, double y, double qT) const
  {
    const double Q2 = Q*Q;
    const double MT = std::hypot(Q, qT);

    std::vector<double> points{_etamin, _etamax};

    // No open region if MT <= 2h
    if (MT <= 2. * _pTmin1) return points;

    auto add = [&](double eta)
    {
      if (std::isfinite(eta) && eta > _etamin && eta < _etamax) points.push_back(eta);
    };

    auto add_partner = [&](double eta_known, double pT_known)
    {
      const double delta = eta_known - y;
      const double plus  = MT - pT_known * std::exp(delta);
      const double minus = MT - pT_known * std::exp(-delta);

      if (plus > 0. && minus > 0.) add(y + 0.5 * (std::log(plus / minus)));
    };

    add(y);

    // lT1 = -qTCut at x = +/-1.
    // Exchanging leptons also supplies lT2 = -qTCut at x = +/-1.
    for (double sign :
         {
           -1., 1.
         })
      {
        const double c = (Q2 / (2. * _pTmin1) + sign * qT) / MT;
        if (c < 1.) continue;

        const double delta = std::acosh(c);

        for (double side :
             {
               -1., 1.
             })
          {
            const double eta = y + side * delta;
            add(eta);
            add_partner(eta, _pTmin1);
          }
      }

    // Both transverse-momentum cuts saturated:
    // eta1 + eta2 = 2y and MT = 2h cosh(eta1 - z).
    const double delta_both = std::acosh(MT / (2. * _pTmin1));
    add(y - delta_both);
    add(y + delta_both);

    for (double eta_bar :
         {
           _etamin, _etamax
         })
      {
        // pT1 = h and eta2 = eta_bar.
        add(eta_bar + std::asinh((MT / _pTmin1) * std::sinh(y - eta_bar)));

        // pT2 = h and eta2 = eta_bar.
        add_partner(eta_bar, _pTmin1);

        // eta2 = eta_bar at a collinear transverse configuration.
        // Parametrize the second lepton with cos(phi2) = +/-1.
        const double c = std::cosh(eta_bar - y);

        for (double sign :
             {
               -1., 1.
             })
          {
            const double pT2 = Q2 / (2. * (MT * c - sign * qT));
            add_partner(eta_bar, pT2);
          }
      }

    std::sort(points.begin(), points.end());

    // Different transition formulas can produce numerically close breakpoints
    constexpr double eta_merge_tol = 1.0e-12;

    std::vector<double> out_points;
    out_points.reserve(points.size());
    out_points.push_back(_etamin);

    for (const double eta : points)
      {
        if (eta <= _etamin || eta >= _etamax) continue;

        if (eta - out_points.back() <= eta_merge_tol || _etamax - eta <= eta_merge_tol)
          {
            continue;
          }

        out_points.push_back(eta);
      }

    out_points.push_back(_etamax);

    return out_points;
  }

  //_________________________________________________________________________
  double TwoBodyPhaseSpace::PhaseSpaceReduction(double const& Q, double const& y, double const& qTin)
  {
    // Return automatically zero if "y" is larger that "_etamax"
    if (y <= _etamin || y >= _etamax) return 0.;

    constexpr double qT_rel_epsilon = 1.e-8;
    if (qTin <= qT_rel_epsilon * Q) return ZeroqTPhaseSpaceReduction(Q, y);

    const double pref = 3. * Q*Q / (8. * M_PI);
    const double MT   = std::hypot(Q, qTin);


    auto fnc = [&](double const& eta) -> double { return _integrand(eta, qTin, y, MT); };
    const auto points = eta_breakpoints(Q, y, qTin);

    const apfel::Integrator Ieta{fnc};

    double result = 0.;

    for (std::size_t i = 1; i < points.size(); i++)
      {
        try
          {
            result += Ieta.integrate(points[i - 1], points[i], _eps);
          }
        catch (...)
          {
            std::cout << "TwoBodyPhaseSpace integral between " << points[i-1]
               << " and " << points[i] << " failed to reach accuracy.\n";
          }
      }

    return pref * result;
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
    const std::function<double(double const&)> IntegrandP = [&] (double const& eta) -> double
    {
      // More definitions
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
      const std::function<double(double const&)> Hbar = [&] (double const& x) -> double
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
