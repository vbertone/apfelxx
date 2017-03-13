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
#include <apfel/gridalphaqcd.h>
#include <apfel/tabulateobject.h>
#include <apfel/timer.h>

using namespace apfel;
using namespace std;

int main()
{
  // Constructor of QGrid
  const QGrid<double> qg{50, 1, 1000, 3, {0, 0, 0, sqrt(2), 4.5, 175.}};
  cout << qg << endl;

  // Direct AlphaQCD
  const AlphaQCD as{0.35, sqrt(2), {0, 0, 0, sqrt(2), 4.5, 175}, 2};
  MatchedEvolution<double> *tas = new AlphaQCD{0.35, sqrt(2), {0, 0, 0, sqrt(2), 4.5, 175}, 2};

  // Tabulate AlphaQCD on a QGrid
  //const GridAlphaQCD gas{0.35, sqrt(2), {0, 0, 0, sqrt(2), 4.5, 175}, 2, 50, 1, 1000, 3};
  const TabulateObject<double> gas{tas, 50, 1, 1000, 3};

  cout << "Precision test ..." << endl;
  auto nQ   = 20;
  auto Qmin = 1.1;
  auto Qmax = 999.;
  auto Step = exp( log( Qmax / Qmin ) / ( nQ - 1 ) );
  auto Q = Qmin;
  cout << setprecision(8) << scientific;
  cout << "Q       \t\tDirect  \t\tInterpolated\t\tRatio" << endl;
  for (auto iQ = 0; iQ < nQ; iQ++)
    {
      cout << Q << "\t\t" << as.Evaluate(Q) << "\t\t" << gas.Evaluate(Q) << "\t\t" << as.Evaluate(Q) / gas.Evaluate(Q) << endl;
      //cout << Q << "\t\t" << as.Evaluate(Q) << "\t\t" << tas.Evaluate(Q) << "\t\t" << as.Evaluate(Q) / gas.Evaluate(Q) << endl;
      Q *=Step;
    }

  cout << "\nSpeed test ..." << endl;
  Timer t;
  nQ   = 1000000;
  Step = exp( log( Qmax / Qmin ) / ( nQ - 1 ) );

  t.start();
  cout << "Direct calculation of " << nQ << " points ..." << endl;
  Q = Qmin;
  for (auto iQ = 0; iQ < nQ; iQ++)
    {
      as.Evaluate(Q);
      Q *=Step;
    }
  t.printTime(t.stop());

  t.start();
  cout << "Interpolated calculation of " << nQ << " points ..." << endl;
  Q = Qmin;
  for (auto iQ = 0; iQ < nQ; iQ++)
    {
      gas.Evaluate(Q);
      Q *=Step;
    }
  t.printTime(t.stop());

  return 0;
}
