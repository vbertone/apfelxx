//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <cmath>
#include <map>
#include <iomanip>

#include <apfel/dglapbuilder.h>
#include <apfel/structurefunctionbuilder.h>
#include <apfel/grid.h>
#include <apfel/timer.h>
#include <apfel/tools.h>
#include <apfel/alphaqcd.h>
#include <apfel/tabulateobject.h>

using namespace apfel;
using namespace std;

// LH Toy PDFs
double xupv(double const& x)  { return 5.107200 * pow(x,0.8) * pow((1-x),3); }
double xdnv(double const& x)  { return 3.064320 * pow(x,0.8) * pow((1-x),4); }
double xglu(double const& x)  { return 1.7 * pow(x,-0.1) * pow((1-x),5); }
double xdbar(double const& x) { return 0.1939875 * pow(x,-0.1) * pow((1-x),6); }
double xubar(double const& x) { return xdbar(x) * (1-x); }
double xsbar(double const& x) { return 0.2 * ( xdbar(x) + xubar(x) ); }
unordered_map<int,double> LHToyPDFs(double const& x, double const&)
{
  // Call all functions once.
  const double upv  = xupv (x);
  const double dnv  = xdnv (x);
  const double glu  = xglu (x);
  const double dbar = xdbar(x);
  const double ubar = xubar(x);
  const double sbar = xsbar(x);

  // Construct QCD evolution basis conbinations.
  double const Gluon   = glu;
  double const Singlet = dnv + 2 * dbar + upv + 2 * ubar + 2 * sbar;
  double const T3      = upv + 2 * ubar - dnv - 2 * dbar;
  double const T8      = upv + 2 * ubar + dnv + 2 * dbar - 4 * sbar;
  double const Valence = upv + dnv;
  double const V3      = upv - dnv;

  // Fill in map in the QCD evolution basis.
  unordered_map<int,double> QCDEvMap;
  QCDEvMap.insert({0 , Gluon});
  QCDEvMap.insert({1 , Singlet});
  QCDEvMap.insert({2 , Valence});
  QCDEvMap.insert({3 , T3});
  QCDEvMap.insert({4 , V3});
  QCDEvMap.insert({5 , T8});
  QCDEvMap.insert({6 , Valence});
  QCDEvMap.insert({7 , Singlet});
  QCDEvMap.insert({8 , Valence});
  QCDEvMap.insert({9 , Singlet});
  QCDEvMap.insert({10, Valence});
  QCDEvMap.insert({11, Singlet});
  QCDEvMap.insert({12, Valence});

  return QCDEvMap;
}

int main()
{
  // x-space grid
  const Grid g{{SubGrid{100,1e-5,3}, SubGrid{60,1e-1,3}, SubGrid{50,6e-1,3}, SubGrid{50,8e-1,3}}};

  // Initial scale
  const double mu0 = sqrt(2);

  // Vectors of masses and thresholds
  const vector<double> Masses = {0, 0, 0, sqrt(2), 4.5, 175}; // Check in the level above that they are ordered
  const vector<double> Thresholds = Masses;

  // Perturbative order
  const int PerturbativeOrder = 2;

  // Running coupling
  const double AlphaQCDRef = 0.35;
  const double MuAlphaQCDRef = sqrt(2);
  AlphaQCD a{AlphaQCDRef, MuAlphaQCDRef, Masses, PerturbativeOrder};
  const TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // Charges
  function<vector<double>(double const&)> fBq = [Thresholds] (double const& Q) -> vector<double>
    {
      vector<double> Bq;
      for (auto i = 0; i < (int) Thresholds.size(); i++)
	Bq.push_back((Q > Thresholds[i] ? QCh2[i] : 0));
      return Bq;
    };
  function<vector<double>(double const&)> fDq = [Thresholds] (double const&) -> vector<double>{ return {0, 0, 0, 0, 0, 0}; };

  // Initialize DGLAP evolution
  auto EvolvedPDFs = DglapBuildQCD(g, LHToyPDFs, mu0, Masses, Thresholds, PerturbativeOrder, as);

  // Tabulate PDFs
  const TabulateObject<Set<Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3};

  // Evolved PDFs
  const auto PDFs = [&] (double const& x, double const& Q) -> unordered_map<int,double>{ return TabulatedPDFs.EvaluateMapxQ(x,Q); };

  // Initialize structure functions
  const auto F2 = F2NCBuildZM(g, PDFs, Thresholds, PerturbativeOrder, as, fBq);
  const auto FL = FLNCBuildZM(g, PDFs, Thresholds, PerturbativeOrder, as, fBq);
  const auto F3 = F3NCBuildZM(g, PDFs, Thresholds, PerturbativeOrder, as, fDq);

  const TabulateObject<Distribution> F2total {[&] (double const& Q) -> Distribution{ return F2.at(0).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const TabulateObject<Distribution> F2light {[&] (double const& Q) -> Distribution{ return F2.at(1).Evaluate(Q) + F2.at(2).Evaluate(Q) + F2.at(3).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const TabulateObject<Distribution> F2charm {[&] (double const& Q) -> Distribution{ return F2.at(4).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const TabulateObject<Distribution> F2bottom{[&] (double const& Q) -> Distribution{ return F2.at(5).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};

  Timer t;
  t.start();

  // Final scale
  const auto Q = 100;

  cout << scientific << endl;
  cout << "Alphas(Q) = " << as(Q) << endl;
  cout << endl;

  const vector<double> xlha = { 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1 };

  cout << "    x   "
       << "  F2light   "
       << "  F2charm   "
       << "  F2bottom  "
       << "  F2total   "
       << endl;
  for (auto i = 2; i < (int) xlha.size(); i++)
    cout << setprecision(1) << xlha[i] << "  " << setprecision(4)
	 << (F2.at(1).Evaluate(Q)+F2.at(2).Evaluate(Q)+F2.at(3).Evaluate(Q)).Evaluate(xlha[i]) << "  "
	 << F2.at(4).Evaluate(Q).Evaluate(xlha[i]) << "  "
	 << F2.at(5).Evaluate(Q).Evaluate(xlha[i]) << "  "
	 << F2.at(0).Evaluate(Q).Evaluate(xlha[i]) << "  "
	 << endl;
  cout << endl;

  for (auto i = 2; i < (int) xlha.size(); i++)
    cout << setprecision(1) << xlha[i] << "  " << setprecision(4)
	 << F2light.EvaluatexQ(xlha[i],Q) << "  "
	 << F2charm.EvaluatexQ(xlha[i],Q) << "  "
	 << F2bottom.EvaluatexQ(xlha[i],Q) << "  "
	 << F2total.EvaluatexQ(xlha[i],Q) << "  "
	 << endl;
  cout << endl;

  cout << "    x   "
       << "  FLlight   "
       << "  FLcharm   "
       << "  FLbottom  "
       << "  FLtotal   "
       << endl;
  for (auto i = 2; i < (int) xlha.size(); i++)
    cout << setprecision(1) << xlha[i] << "  " << setprecision(4)
	 << (FL.at(1).Evaluate(Q)+FL.at(2).Evaluate(Q)+FL.at(3).Evaluate(Q)).Evaluate(xlha[i]) << "  "
	 << FL.at(4).Evaluate(Q).Evaluate(xlha[i]) << "  "
	 << FL.at(5).Evaluate(Q).Evaluate(xlha[i]) << "  "
	 << FL.at(0).Evaluate(Q).Evaluate(xlha[i]) << "  "
	 << endl;
  cout << endl;

  cout << "    x   "
       << "  F3light   "
       << "  F3charm   "
       << "  F3bottom  "
       << "  F3total   "
       << endl;
  for (auto i = 2; i < (int) xlha.size(); i++)
    cout << setprecision(1) << xlha[i] << "  " << setprecision(4)
	 << (F3.at(1).Evaluate(Q)+F3.at(2).Evaluate(Q)+F3.at(3).Evaluate(Q)).Evaluate(xlha[i]) << "  "
	 << F3.at(4).Evaluate(Q).Evaluate(xlha[i]) << "  "
	 << F3.at(5).Evaluate(Q).Evaluate(xlha[i]) << "  "
	 << F3.at(0).Evaluate(Q).Evaluate(xlha[i]) << "  "
	 << endl;
  cout << endl;

  t.stop();

  const int k = 1000000;
  cout << "Interpolating " << k << " times F2 on the grid... ";
  t.start();
  for (auto i = 0; i < k; i++)
    F2total.EvaluatexQ(0.05,Q);
  t.stop();

  return 0;
}
