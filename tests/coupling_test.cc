//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <apfel/coupling.h>
#include <apfel/tools.h>
#include <cmath>
#include <vector>
#include <iostream>

using namespace apfel;
using namespace std;

class AlphaQCD: public Coupling
{
public:
  AlphaQCD(double const& AlphaRef, double const& MuRef, vector<double> const& Masses, int const& pt = 0):
    Coupling(AlphaRef, MuRef, Masses),
    _pt(pt)
  {}

  double Coup(int const& nf, double const& as0, double const& mu02, double const& mu2) const
  {
    auto lrrat = log(mu2/mu02);
    // Analytical solution at leading order
    if ( _pt == 0 )
      return as0 / ( 1 + beta0(nf) * as0 * lrrat );
    else
      {
	const int nstep = 10;
	auto as  = as0;
	auto dlr = lrrat / nstep;
	// Numerical solution of the evolution equation
	// with fourth-order Runge-Kutta beyond leading order
	for (auto k1 = 0; k1 < nstep; k1++)
	  {
            auto xk0 = dlr * fbeta(as            ,nf);
            auto xk1 = dlr * fbeta(as + 0.5 * xk0,nf);
            auto xk2 = dlr * fbeta(as + 0.5 * xk1,nf);
            auto xk3 = dlr * fbeta(as +       xk2,nf);
            as      += ( xk0 + 2 * xk1 + 2 * xk2 + xk3 ) / 6.;
	  }
	return as;
      }
  }

  double MatchCoupling(bool const& Up, double const& Coup, double const& LogKth) const
  {
    if (_pt == 0) return Coup;
    else if (_pt == 1)
      {
	double c1;
	if (Up)
	  c1 = 2. / 3. * LogKth;
	else
	  c1 = - 2. / 3. * LogKth;
	return Coup * ( 1 + c1 * Coup / FourPi );
      }
    else if (_pt == 2)
      {
	double c1, c2;
	if (Up)
	  {
	    c1 = 2. / 3. * LogKth;
	    c2 = 4. / 9. * pow(LogKth,2) + 38. / 3. * LogKth + 14. / 3.;
	  }
	else
	  {
	    c1 = - 2. / 3. * LogKth;
	    c2 = 4. / 9. * pow(LogKth,2) - 38. / 3. * LogKth - 14. / 3.;
	  }
	return Coup * ( 1 + c1 * Coup / FourPi + c2 * pow(Coup / FourPi,2) );
      }
    else return 0;
  }

  double beta0(int const& nf) const { return ( 33. - 2. * nf ) / 3. / FourPi; }
  double beta1(int const& nf) const { return ( 102. - 38. / 3. * nf ) / pow(FourPi,2); }
  double beta2(int const& nf) const { return ( 2857. / 2. - 5033. / 18. * nf + 325. / 54. * pow(nf,2) ) / pow(FourPi,3); }
  double fbeta(double const& as, int const& nf) const
  {
    auto bt = - pow(as,2) * beta0(nf);
    if (_pt >= 1) bt -= pow(as,3) * beta1(nf);
    if (_pt >= 2) bt -= pow(as,4) * beta2(nf);
    return bt;
  }

private:
  int const _pt;
};

int main()
{
  const double         AlphaRef = 0.35;
  const double         MuRef    = sqrt(2);
  const vector<double> Masses   = {0, 0, 0, sqrt(2), 4.5, 175};
  const AlphaQCD asLO(AlphaRef, MuRef, Masses);
  const AlphaQCD asNLO(AlphaRef, MuRef, Masses, 1);
  const AlphaQCD asNNLO(AlphaRef, MuRef, Masses, 2);

  const double Mu = 100;

  cout << "\nLO:    alpha_s(Mu = " << Mu << " GeV) = " << asLO.GetCoupling(Mu) << endl;
  cout << "NLO:   alpha_s(Mu = " << Mu << " GeV) = " << asNLO.GetCoupling(Mu)  << endl;
  cout << "NNLO:  alpha_s(Mu = " << Mu << " GeV) = " << asNNLO.GetCoupling(Mu) << "\n" << endl;

  return 0;
}
