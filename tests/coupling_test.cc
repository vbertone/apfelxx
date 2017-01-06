//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <apfel/coupling.h>
#include <cmath>
#include <vector>
#include <iostream>

using namespace apfel;
using namespace std;

class AlphaQCD: public Coupling
{
public:
  AlphaQCD(double const& AlphaRef, double const& MuRef, vector<double> const& Masses): Coupling(AlphaRef, MuRef, Masses) {}
  double Coup(int const& nf, double const& as0, double const& mu02, double const& mu2) const { return as0 / ( 1 + beta0(nf) * as0 * log(mu2/mu02) ); }
  double MatchCoupling(bool const& Up, double const& Coup, double const& LogKth)       const { if(Up) {}; return Coup * ( 1 + Coup * LogKth ); }
  double beta0(int const& nf) const { return ( 33. - 2. * nf ) / 3. / ( 4. * M_PI ); }
};

int main()
{
  const double         AlphaRef = 0.35;
  const double         MuRef    = sqrt(2);
  const vector<double> Masses   = {0, 0, 0, sqrt(2), 4.5, 175};
  const AlphaQCD as(AlphaRef, MuRef, Masses);

  const double Mu = 100;
  cout << "\nalpha_s(Mu = " << Mu << " GeV) = " << as.GetCoupling(Mu) << "\n" << endl;

  return 0;
}
