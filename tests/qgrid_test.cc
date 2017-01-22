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
 * @brief Alpha_s on a QGrid: (version 1)
 * This class fills out the _GridValues vector by a direct calculation
 * of the coupling starting alwyes from the reference values given in
 * the constructore of AlphaQCD. This is the simplest way best not the
 * most performing because is is assumed that alpha_s is computed with
 * a good accuracy for any distance between initial and final scales.
 * The consequence is that one has to do many steps (~10) in the Runge-
 * Kutta algorithm.
 */
class GridAlphaQCD1: public QGrid<double>
{
public:
  GridAlphaQCD1(AlphaQCD const& as,
		int      const& nQ,
		double   const& QMin,
		double   const& QMax,
		int      const& InterDegree):
    QGrid<double>(nQ, QMin, QMax, InterDegree, as.GetThresholds())
  {
    for (auto const& iQ : _Qg)
      _GridValues.push_back(as.GetCoupling(iQ));
  }
};

/**
 * @brief Alpha_s on a QGrid: (version 2)
 * This class fills out the _GridValues vector step by step as the
 * calculation of the next step is always computed by starting from
 * the previous step. This method is more performing because the steps
 * are typically small and thus it is sufficient to make one single
 * step in the Runge-Kutta algorithm.
 */
class GridAlphaQCD2: public QGrid<double>
{
public:
  GridAlphaQCD2(double         const& AlphaRef,
		double         const& MuRef,
		vector<double> const& Masses,
		int            const& pt,
		int            const& nQ,
		double         const& QMin,
		double         const& QMax,
		int            const& InterDegree):
    QGrid<double>(nQ, QMin, QMax, InterDegree, Masses)
  {
    // Initialze alphas
    AlphaQCD as{AlphaRef, MuRef, Masses, pt, 1};

    // Find the point on the QGrid right below MuRef
    const auto tQ = lower_bound(_Qg.begin(), _Qg.end(), MuRef) - _Qg.begin() - 1;

    // Resize container
    _GridValues.resize(_Qg.size());

    // Loop on "_Qg" below "MuRef"
    for (auto iQ = tQ; iQ >= 0; iQ--)
      {
	auto a = as.GetCoupling(_Qg[iQ]);
	_GridValues[iQ] = a;
	as.SetAlphaRef(a);
	as.SetMuRef(_Qg[iQ]);
      }

    // Loop on "_Qg" above "MuRef"
    as.SetAlphaRef(AlphaRef);
    as.SetMuRef(MuRef);
    for (auto iQ = tQ + 1; iQ < (int) _Qg.size(); iQ++)
      {
	auto a = as.GetCoupling(_Qg[iQ]);
	_GridValues[iQ] = a;
	as.SetAlphaRef(a);
	as.SetMuRef(_Qg[iQ]);
      }
  }
};

int main()
{
  // Constructor of QGrid
  const QGrid<double> qg{50, 1, 1000, 3, {0, 0, 0, sqrt(2), 4.5, 175.}};
  cout << qg << endl;

  // Tabulate AlphaQCD on a QGrid (version 1)
  const AlphaQCD as{0.35, sqrt(2), {0, 0, 0, sqrt(2), 4.5, 175}, 2};
  const GridAlphaQCD1 gas1{as, 50, 1, 1000, 3};

  // Tabulate AlphaQCD on a QGrid (version 2)
  const GridAlphaQCD2 gas2{0.35, sqrt(2), {0, 0, 0, sqrt(2), 4.5, 175}, 2, 50, 1, 1000, 3};

  cout << "Precision test ..." << endl;
  auto nQ   = 20;
  auto Qmin = 1.1;
  auto Qmax = 999.;
  auto Step = exp( log( Qmax / Qmin ) / ( nQ - 1 ) );
  auto Q = Qmin;
  cout << setprecision(8) << scientific;
  cout << "Q       \t\tDirect  \t\tInterpolated(1)\t\tInterpolated(2)\t\tRatio(1)   \t\tRatio(2)" << endl;
  for (auto iQ = 0; iQ < nQ; iQ++)
    {
      cout << Q << "\t\t" << as.GetCoupling(Q) << "\t\t" << gas1.Evaluate(Q) << "\t\t" << gas2.Evaluate(Q) << "\t\t" << as.GetCoupling(Q) / gas1.Evaluate(Q) << "\t\t" << as.GetCoupling(Q) / gas2.Evaluate(Q) << endl;
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
      as.GetCoupling(Q);
      Q *=Step;
    }
  t.printTime(t.stop());

  t.start();
  cout << "Interpolated(1) calculation of " << nQ << " points ..." << endl;
  Q = Qmin;
  for (auto iQ = 0; iQ < nQ; iQ++)
    {
      gas1.Evaluate(Q);
      Q *=Step;
    }
  t.printTime(t.stop());

  t.start();

  cout << "Interpolated(2) calculation of " << nQ << " points ..." << endl;
  Q = Qmin;
  for (auto iQ = 0; iQ < nQ; iQ++)
    {
      gas2.Evaluate(Q);
      Q *=Step;
    }
  t.printTime(t.stop());
  return 0;
}
