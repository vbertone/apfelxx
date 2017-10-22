//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <apfel/grid.h>
#include <apfel/timer.h>
#include <apfel/alphaqcd.h>
#include <apfel/tabulateobject.h>
#include <apfel/evolutionbasisqcd.h>
#include <apfel/matchingbasisqcd.h>
#include <apfel/dglap.h>
#include <apfel/dglapbuilder.h>
#include <apfel/lhtoypdfs.h>

#include <functional>

using namespace apfel;
using namespace std;

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

  // Initialize QCD evolution objects
  const auto DglapObj = InitializeDglapObjectsQCD(g, Masses, Thresholds);

  // Construct the DGLAP object
  auto EvolvedPDFs = BuildDglap(DglapObj, LHToyPDFs, mu0, PerturbativeOrder, as);

  // Tabulate PDFs
  const TabulateObject<Set<Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3};

  // Final scale
  double mu = 100;

  // Print results
  cout << scientific;

  // Evolve PDFs to the final Scale
  cout << "Direct evolution (4th order Runge-Kutta) from Q0 = " << mu0 << " GeV to Q = " << mu << " GeV... ";
  Timer t;
  auto pdfs = EvolvedPDFs->Evaluate(mu);
  t.stop();

  cout << "Interpolation of the tabulated PDFs... ";
  t.start();
  auto tpdfs = TabulatedPDFs.Evaluate(mu);
  t.stop();

  double xlha[] = {1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};

  cout << "\nAlphaQCD(Q) = " << Alphas.Evaluate(mu) << endl;
  cout << "\n   x    "
       << "   u-ubar   "
       << "   d-dbar   "
       << " 2(ubr+dbr) "
       << "   c+cbar   "
       << "    gluon   "
       << endl;

  cout << "Direct Evolution:" << endl;
  for (auto i = 2; i < 11; i++)
    {
      cout.precision(1);
      cout << xlha[i];
      cout.precision(4);
      cout << "  " <<
	pdfs.at(2).Evaluate(xlha[i])  / 6  +
	pdfs.at(4).Evaluate(xlha[i])  / 2  +
	pdfs.at(6).Evaluate(xlha[i])  / 6  +
	pdfs.at(8).Evaluate(xlha[i])  / 12 +
	pdfs.at(10).Evaluate(xlha[i]) / 20 +
	pdfs.at(12).Evaluate(xlha[i]) / 30
	   << "  " <<
	pdfs.at(2).Evaluate(xlha[i])  / 6  -
	pdfs.at(4).Evaluate(xlha[i])  / 2  +
	pdfs.at(6).Evaluate(xlha[i])  / 6  +
	pdfs.at(8).Evaluate(xlha[i])  / 12 +
	pdfs.at(10).Evaluate(xlha[i]) / 20 +
	pdfs.at(12).Evaluate(xlha[i]) / 30
	   << "  " <<
	( pdfs.at(1).Evaluate(xlha[i])  - pdfs.at(2).Evaluate(xlha[i])  ) / 3  +
	( pdfs.at(5).Evaluate(xlha[i])  - pdfs.at(6).Evaluate(xlha[i])  ) / 3  +
	( pdfs.at(7).Evaluate(xlha[i])  - pdfs.at(8).Evaluate(xlha[i])  ) / 6  +
	( pdfs.at(9).Evaluate(xlha[i])  - pdfs.at(10).Evaluate(xlha[i]) ) / 10 +
	( pdfs.at(11).Evaluate(xlha[i]) - pdfs.at(12).Evaluate(xlha[i]) ) / 15
	   << "  " <<
	pdfs.at(1).Evaluate(xlha[i])  / 6  -
	pdfs.at(7).Evaluate(xlha[i])  / 4  +
	pdfs.at(9).Evaluate(xlha[i])  / 20 +
	pdfs.at(11).Evaluate(xlha[i]) / 30
	   << "  " <<
	pdfs.at(0).Evaluate(xlha[i]) << "  "
	   << endl;
    }
  cout << "      " << endl;

  cout << "Interpolation on the PDF table (all x for each Q):" << endl;
  for (auto i = 2; i < 11; i++)
    {
      cout.precision(1);
      cout << xlha[i];
      cout.precision(4);
      cout << "  " <<
	tpdfs.at(2).Evaluate(xlha[i])  / 6  +
	tpdfs.at(4).Evaluate(xlha[i])  / 2  +
	tpdfs.at(6).Evaluate(xlha[i])  / 6  +
	tpdfs.at(8).Evaluate(xlha[i])  / 12 +
	tpdfs.at(10).Evaluate(xlha[i]) / 20 +
	tpdfs.at(12).Evaluate(xlha[i]) / 30
	   << "  " <<
	tpdfs.at(2).Evaluate(xlha[i])  / 6  -
	tpdfs.at(4).Evaluate(xlha[i])  / 2  +
	tpdfs.at(6).Evaluate(xlha[i])  / 6  +
	tpdfs.at(8).Evaluate(xlha[i])  / 12 +
	tpdfs.at(10).Evaluate(xlha[i]) / 20 +
	tpdfs.at(12).Evaluate(xlha[i]) / 30
	   << "  " <<
	( tpdfs.at(1).Evaluate(xlha[i])  - tpdfs.at(2).Evaluate(xlha[i])  ) / 3  +
	( tpdfs.at(5).Evaluate(xlha[i])  - tpdfs.at(6).Evaluate(xlha[i])  ) / 3  +
	( tpdfs.at(7).Evaluate(xlha[i])  - tpdfs.at(8).Evaluate(xlha[i])  ) / 6  +
	( tpdfs.at(9).Evaluate(xlha[i])  - tpdfs.at(10).Evaluate(xlha[i]) ) / 10 +
	( tpdfs.at(11).Evaluate(xlha[i]) - tpdfs.at(12).Evaluate(xlha[i]) ) / 15
	   << "  " <<
	tpdfs.at(1).Evaluate(xlha[i])  / 6  -
	tpdfs.at(7).Evaluate(xlha[i])  / 4  +
	tpdfs.at(9).Evaluate(xlha[i])  / 20 +
	tpdfs.at(11).Evaluate(xlha[i]) / 30
	   << "  " <<
	tpdfs.at(0).Evaluate(xlha[i]) << "  "
	   << endl;
    }
  cout << "      " << endl;

  cout << "Interpolation on the PDF table (x and Q independently):" << endl;
  for (auto i = 2; i < 11; i++)
    {
      cout.precision(1);
      cout << xlha[i];
      cout.precision(4);
      cout << "  " <<
	TabulatedPDFs.EvaluatexQ(2,xlha[i],mu)  / 6  +
	TabulatedPDFs.EvaluatexQ(4,xlha[i],mu)  / 2  +
	TabulatedPDFs.EvaluatexQ(6,xlha[i],mu)  / 6  +
	TabulatedPDFs.EvaluatexQ(8,xlha[i],mu)  / 12 +
	TabulatedPDFs.EvaluatexQ(10,xlha[i],mu) / 20 +
	TabulatedPDFs.EvaluatexQ(12,xlha[i],mu) / 30
	   << "  " <<
	TabulatedPDFs.EvaluatexQ(2,xlha[i],mu)  / 6  -
	TabulatedPDFs.EvaluatexQ(4,xlha[i],mu)  / 2  +
	TabulatedPDFs.EvaluatexQ(6,xlha[i],mu)  / 6  +
	TabulatedPDFs.EvaluatexQ(8,xlha[i],mu)  / 12 +
	TabulatedPDFs.EvaluatexQ(10,xlha[i],mu) / 20 +
	TabulatedPDFs.EvaluatexQ(12,xlha[i],mu) / 30
	   << "  " <<
	( TabulatedPDFs.EvaluatexQ(1,xlha[i],mu)  - TabulatedPDFs.EvaluatexQ(2,xlha[i],mu)  ) / 3  +
	( TabulatedPDFs.EvaluatexQ(5,xlha[i],mu)  - TabulatedPDFs.EvaluatexQ(6,xlha[i],mu)  ) / 3  +
	( TabulatedPDFs.EvaluatexQ(7,xlha[i],mu)  - TabulatedPDFs.EvaluatexQ(8,xlha[i],mu)  ) / 6  +
	( TabulatedPDFs.EvaluatexQ(9,xlha[i],mu)  - TabulatedPDFs.EvaluatexQ(10,xlha[i],mu) ) / 10 +
	( TabulatedPDFs.EvaluatexQ(11,xlha[i],mu) - TabulatedPDFs.EvaluatexQ(12,xlha[i],mu) ) / 15
	   << "  " <<
	TabulatedPDFs.EvaluatexQ(1,xlha[i],mu)  / 6  -
	TabulatedPDFs.EvaluatexQ(7,xlha[i],mu)  / 4  +
	TabulatedPDFs.EvaluatexQ(9,xlha[i],mu)  / 20 +
	TabulatedPDFs.EvaluatexQ(11,xlha[i],mu) / 30
	   << "  " <<
	TabulatedPDFs.EvaluatexQ(0,xlha[i],mu) << "  "
	   << endl;
    }
  cout << "      " << endl;

  cout << "Interpolation on the PDF table as a map (x and Q independently):" << endl;
  for (auto i = 2; i < 11; i++)
    {
      const auto DistMap = TabulatedPDFs.EvaluateMapxQ(xlha[i],mu);
      cout.precision(1);
      cout << xlha[i];
      cout.precision(4);
      cout << "  " <<
	DistMap.at(2)  / 6  +
	DistMap.at(4)  / 2  +
	DistMap.at(6)  / 6  +
	DistMap.at(8)  / 12 +
	DistMap.at(10) / 20 +
	DistMap.at(12) / 30
	   << "  " <<
	DistMap.at(2)  / 6  -
	DistMap.at(4)  / 2  +
	DistMap.at(6)  / 6  +
	DistMap.at(8)  / 12 +
	DistMap.at(10) / 20 +
	DistMap.at(12) / 30
	   << "  " <<
	( DistMap.at(1)  - DistMap.at(2)  ) / 3  +
	( DistMap.at(5)  - DistMap.at(6)  ) / 3  +
	( DistMap.at(7)  - DistMap.at(8)  ) / 6  +
	( DistMap.at(9)  - DistMap.at(10) ) / 10 +
	( DistMap.at(11) - DistMap.at(12) ) / 15
	   << "  " <<
	DistMap.at(1)  / 6  -
	DistMap.at(7)  / 4  +
	DistMap.at(9)  / 20 +
	DistMap.at(11) / 30
	   << "  " <<
	DistMap.at(0) << "  "
	   << endl;
    }
  cout << "      " << endl;

  int k = 1000000;
  cout << "Interpolating " << k << " times a single PDF on the (x,Q) grid... ";
  t.start();
  for (auto i = 0; i < k; i++)
    TabulatedPDFs.EvaluatexQ(0,0.05,mu);
  t.stop();

  k = 100000;
  cout << "Interpolating " << k << " times an map of PDFs on the (x,Q) grid... ";
  t.start();
  for (auto i = 0; i < k; i++)
    TabulatedPDFs.EvaluateMapxQ(0.05,mu);
  t.stop();

  return 0;
}
