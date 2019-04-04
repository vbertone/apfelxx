//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/dglapbuilder.h>
#include <apfel/structurefunctionbuilder.h>
#include <apfel/grid.h>
#include <apfel/timer.h>
#include <apfel/constants.h>
#include <apfel/alphaqcd.h>
#include <apfel/tabulateobject.h>
#include <apfel/lhtoypdfs.h>

#include <cmath>
#include <map>
#include <iomanip>

int main()
{
  //SetVerbosityLevel(0);
  apfel::Banner();
  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{100,1e-5,3}, apfel::SubGrid{60,1e-1,3}, apfel::SubGrid{50,6e-1,3}, apfel::SubGrid{50,8e-1,3}}};

  // Initial scale
  const double mu0 = sqrt(2);

  // Vectors of masses and thresholds
  const std::vector<double> Masses = {0, 0, 0, sqrt(2), 4.5, 175}; // Check in the level above that they are ordered
  const std::vector<double> Thresholds = Masses;

  // Perturbative order
  const int PerturbativeOrder = 2;

  // Running coupling
  const double AlphaQCDRef = 0.35;
  const double MuAlphaQCDRef = sqrt(2);
  apfel::AlphaQCD a{AlphaQCDRef, MuAlphaQCDRef, Masses, PerturbativeOrder};
  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // Effective charges.
  std::function<std::vector<double>(double const&)> fBq = [=] (double const&) -> std::vector<double>{ return apfel::QCh2; };
  std::function<std::vector<double>(double const&)> fDq = [=] (double const&) -> std::vector<double>{ return {0, 0, 0, 0, 0, 0}; };

  // Initialize QCD evolution objects
  const auto DglapObj = InitializeDglapObjectsQCD(g, Masses, Thresholds);

  // Construct the DGLAP object
  auto EvolvedPDFs = BuildDglap(DglapObj, apfel::LHToyPDFs, mu0, PerturbativeOrder, as);

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3};

  // Evolved PDFs
  const auto PDFs = [&] (double const& x, double const& Q) -> std::map<int,double>{ return TabulatedPDFs.EvaluateMapxQ(x,Q); };

  // Initialize coefficient functions
  const auto F2Obj = InitializeF2NCObjectsZM(g, Thresholds);
  const auto FLObj = InitializeFLNCObjectsZM(g, Thresholds);
  const auto F3Obj = InitializeF3NCObjectsZM(g, Thresholds);

  // Initialize structure functions
  const auto F2 = BuildStructureFunctions(F2Obj, PDFs, PerturbativeOrder, as, fBq);
  const auto FL = BuildStructureFunctions(FLObj, PDFs, PerturbativeOrder, as, fBq);
  const auto F3 = BuildStructureFunctions(F3Obj, PDFs, PerturbativeOrder, as, fDq);

  const apfel::TabulateObject<apfel::Distribution> F2total {[&] (double const& Q) -> apfel::Distribution{ return F2.at(0).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F2light {[&] (double const& Q) -> apfel::Distribution{ return F2.at(1).Evaluate(Q) + F2.at(2).Evaluate(Q) + F2.at(3).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F2charm {[&] (double const& Q) -> apfel::Distribution{ return F2.at(4).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F2bottom{[&] (double const& Q) -> apfel::Distribution{ return F2.at(5).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};

  const apfel::TabulateObject<apfel::Distribution> FLtotal {[&] (double const& Q) -> apfel::Distribution{ return FL.at(0).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> FLlight {[&] (double const& Q) -> apfel::Distribution{ return FL.at(1).Evaluate(Q) + FL.at(2).Evaluate(Q) + FL.at(3).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> FLcharm {[&] (double const& Q) -> apfel::Distribution{ return FL.at(4).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> FLbottom{[&] (double const& Q) -> apfel::Distribution{ return FL.at(5).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};

  const apfel::TabulateObject<apfel::Distribution> F3total {[&] (double const& Q) -> apfel::Distribution{ return F3.at(0).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F3light {[&] (double const& Q) -> apfel::Distribution{ return F3.at(1).Evaluate(Q) + F3.at(2).Evaluate(Q) + F3.at(3).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F3charm {[&] (double const& Q) -> apfel::Distribution{ return F3.at(4).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F3bottom{[&] (double const& Q) -> apfel::Distribution{ return F3.at(5).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};

  apfel::Timer t;

  // Final scale
  const double Q = 100;

  std::cout << std::scientific << std::endl;
  std::cout << "Alphas(Q) = " << as(Q) << std::endl;
  std::cout << std::endl;

  const std::vector<double> xlha = { 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1 };

  std::cout << "    x   "
       << "  F2light   "
       << "  F2charm   "
       << "  F2bottom  "
       << "  F2total   "
       << std::endl;
  for (int i = 2; i < (int) xlha.size(); i++)
    std::cout << std::setprecision(1) << xlha[i] << "  " << std::setprecision(4)
	 << F2light.EvaluatexQ(xlha[i],Q)  << "  "
	 << F2charm.EvaluatexQ(xlha[i],Q)  << "  "
	 << F2bottom.EvaluatexQ(xlha[i],Q) << "  "
	 << F2total.EvaluatexQ(xlha[i],Q)  << "  "
	 << std::endl;
  std::cout << std::endl;

  std::cout << "    x   "
       << "  FLlight   "
       << "  FLcharm   "
       << "  FLbottom  "
       << "  FLtotal   "
       << std::endl;
  for (int i = 2; i < (int) xlha.size(); i++)
    std::cout << std::setprecision(1) << xlha[i] << "  " << std::setprecision(4)
	 << FLlight.EvaluatexQ(xlha[i],Q)  << "  "
	 << FLcharm.EvaluatexQ(xlha[i],Q)  << "  "
	 << FLbottom.EvaluatexQ(xlha[i],Q) << "  "
	 << FLtotal.EvaluatexQ(xlha[i],Q)  << "  "
	 << std::endl;
  std::cout << std::endl;

  std::cout << "    x   "
       << "  F3light   "
       << "  F3charm   "
       << "  F3bottom  "
       << "  F3total   "
       << std::endl;
  for (int i = 2; i < (int) xlha.size(); i++)
    std::cout << std::setprecision(1) << xlha[i] << "  " << std::setprecision(4)
	 << F3light.EvaluatexQ(xlha[i],Q)  << "  "
	 << F3charm.EvaluatexQ(xlha[i],Q)  << "  "
	 << F3bottom.EvaluatexQ(xlha[i],Q) << "  "
	 << F3total.EvaluatexQ(xlha[i],Q)  << "  "
	 << std::endl;
  std::cout << std::endl;

  t.stop();

  const int k = 1000000;
  std::cout << "Interpolating " << k << " times F2 on the grid... ";
  t.start();
  for (int i = 0; i < k; i++)
    F2total.EvaluatexQ(0.05,Q);
  t.stop();

  return 0;
}
