//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <iostream>
#include <cmath>
#include <iomanip>

#include <apfel/qgrid.h>
#include <apfel/alphaqcd.h>
#include <apfel/timer.h>

using namespace apfel;
using namespace std;

/**
 * @brief Alpha_s on a QGrid
 */
class GridAlphaQCD: public QGrid
{
public:
  GridAlphaQCD(AlphaQCD const& as,
	       int const& nQ,
	       double const& QMin,
	       double const& QMax,
	       int const& InterDegree):
    QGrid(nQ, QMin, QMax, InterDegree, as.GetThresholds())
  {
    for (auto const& iQ : _Qg)
      _GridValues.push_back(as.GetCoupling(iQ));
  }
};

int main()
{
  // Constructor of QGrid
  const QGrid qg{50, 1, 1000, 3, {0, 0, 0, sqrt(2), 4.5, 175.}};

  cout << qg << endl;

  // Interpolate AlphaQCD on a QGrid
  const AlphaQCD as{0.35, sqrt(2), {0, 0, 0, sqrt(2), 4.5, 175}, 2};
  const GridAlphaQCD gas{as, 50, 1, 1000, 3};

  // Precision test
  auto nQ   = 20;
  auto Qmin = 1.;
  auto Qmax = 1000.;
  auto Step = exp( log( Qmax / Qmin ) / ( nQ - 1 ) );
  auto Q = Qmin;
  cout << setprecision(8);
  cout << "Q       \t\tDirect  \t\tInterpolated\t\tRatio" << endl;
  for (auto iQ = 0; iQ < nQ; iQ++)
    {
      cout << Q << "\t\t" << as.GetCoupling(Q) << "\t\t" << gas.Evaluate(Q) << "\t\t" << as.GetCoupling(Q) / gas.Evaluate(Q) << endl;
      Q *=Step;
    }

  // Speed test
  Timer t;
  nQ   = 1000000;
  Qmin = 1.1;
  Qmax = 999.;
  Step = exp( log( Qmax / Qmin ) / ( nQ - 1 ) );

  cout << "\nDirect calculation of " << nQ << " points ..." << endl;
  Q = Qmin;
  t.start();
  for (auto iQ = 0; iQ < nQ; iQ++)
    {
      as.GetCoupling(Q);
      Q *=Step;
    }
  t.printTime(t.stop());

  cout << "\nInterpolated calculation of " << nQ << " points ..." << endl;
  Q = Qmin;
  t.start();
  for (auto iQ = 0; iQ < nQ; iQ++)
    {
      gas.Evaluate(Q);
      Q *=Step;
    }
  t.printTime(t.stop());

  return 0;
}
