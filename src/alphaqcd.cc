//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/alphaqcd.h"
#include "apfel/tools.h"

using namespace std;

namespace apfel {

  //_________________________________________________________________________________
  AlphaQCD::AlphaQCD(double const& AlphaRef, double const& MuRef, vector<double> const& Masses, int const& pt):
    Coupling(AlphaRef, MuRef, Masses),
    _pt(pt)
  {}

  //_________________________________________________________________________________
  double AlphaQCD::Coup(int const& nf, double const& as0, double const& mu02, double const& mu2) const
  {
    // Return immediately "as0" if "mu02" and "mu2" are equal
    if (mu02 == mu2) return as0;

    // Numerical solution of the evolution equation with fourth-order Runge-Kutta.
    // Use "nstep" steps for the evolution.
    const auto lrrat = log(mu2/mu02);
    const auto nstep = 10;
    const auto dlr   = lrrat / nstep;
    auto as          = as0 / FourPi;
    for (auto k = 0; k < nstep; k++)
      {
	auto xk0 = dlr * fbeta(as            ,nf);
	auto xk1 = dlr * fbeta(as + 0.5 * xk0,nf);
	auto xk2 = dlr * fbeta(as + 0.5 * xk1,nf);
	auto xk3 = dlr * fbeta(as +       xk2,nf);
	as      += ( xk0 + 2 * xk1 + 2 * xk2 + xk3 ) / 6.;
      }
    return FourPi * as;
  }

  //_________________________________________________________________________________
  double AlphaQCD::MatchCoupling(bool const& Up, double const& Coup, double const& LogKth) const
  {
    const auto sgn = ( Up ? 1 : -1);    
    const auto ep = Coup / FourPi;
    const double c[] = { 1, sgn * 2. / 3. * LogKth, 4. / 9. * pow(LogKth,2) + sgn *  38. / 3. * LogKth + sgn * 14. / 3. };
    double match = 0;
    for (auto i = 0; i <= _pt; i++) match += c[i] * pow(ep,i);
    return Coup * match;
  }

  //_________________________________________________________________________________
  double AlphaQCD::betaQCD(int const& pt, int const& nf) const
  {
    if      ( pt == 0 ) return ( 33. - 2. * nf ) / 3.;
    else if ( pt == 1 ) return 102. - 38. / 3. * nf;
    else if ( pt == 2 ) return 2857. / 2. - 5033. / 18. * nf + 325. / 54. * pow(nf,2);
    else                return 0;
  }

  //_________________________________________________________________________________
  double AlphaQCD::fbeta(double const& as, int const& nf) const
  {
    double bt = 0;
    for (auto i = 0; i <= _pt; i++) bt -= pow(as,i+2) * betaQCD(i,nf);
    return bt;
  }

}
