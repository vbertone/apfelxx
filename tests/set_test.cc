//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <apfel/distribution.h>
#include <apfel/grid.h>
#include <apfel/subgrid.h>
#include <apfel/operator.h>
#include <apfel/expression.h>
#include <apfel/timer.h>
#include <apfel/tools.h>
#include <apfel/alphaqcd.h>
#include <apfel/set.h>
#include <apfel/evolutionbasisqcd.h>
#include <apfel/splittingfunctions.h>

#include <cmath>
#include <map>
#include <functional>

using namespace apfel;
using namespace std;

// LH Toy PDFs
double xupv(double const& x)  { return 5.107200 * pow(x,0.8) * pow((1-x),3); }
double xdnv(double const& x)  { return 3.064320 * pow(x,0.8) * pow((1-x),4); }
double xglu(double const& x)  { return 1.7 * pow(x,-0.1) * pow((1-x),5); }
double xdbar(double const& x) { return 0.1939875 * pow(x,-0.1) * pow((1-x),6); }
double xubar(double const& x) { return xdbar(x) * (1-x); }
double xsbar(double const& x) { return 0.2 * ( xdbar(x) + xubar(x) ); }
double LHToyPDFs(int const& i, double const& x)
{
  // Gluon
  if      (i == EvolutionBasisQCD::GLUON    ) return xglu(x);
  // Singlet, T15, T24, T35
  else if (i == EvolutionBasisQCD::SIGMA   ||
	   i == EvolutionBasisQCD::T15     ||
	   i == EvolutionBasisQCD::T24     ||
	   i == EvolutionBasisQCD::T35      ) return xdnv(x) + 2 * xdbar(x) + xupv(x) + 2 * xubar(x) + 2 * xsbar(x);
  // T3
  else if (i == EvolutionBasisQCD::T3       ) return xupv(x) + 2 * xubar(x) - xdnv(x) - 2 * xdbar(x);
  // T8
  else if (i == EvolutionBasisQCD::T8       ) return xupv(x) + 2 * xubar(x) + xdnv(x) + 2 * xdbar(x) - 4 * xsbar(x);
  // Valence, V8, V15, V24, V35
  else if (i == EvolutionBasisQCD::VALENCE ||
	   i == EvolutionBasisQCD::V8      ||
	   i == EvolutionBasisQCD::V15     ||
	   i == EvolutionBasisQCD::V24     ||
	   i == EvolutionBasisQCD::V35      ) return xupv(x) + xdnv(x);
  // V3
  else if (i == EvolutionBasisQCD::V3       )  return xupv(x) - xdnv(x);
  else              return 0;
}

int main()
{
  // Time counter
  Timer t,ttot;
  ttot.start();

  // Grid
  const Grid g{{SubGrid{80,1e-5,3}, SubGrid{50,1e-1,5}, SubGrid{40,8e-1,5}}};

  // ===============================================================
  // Allocate LO splitting functions operators
  cout << "Initializing operators ..." << endl;
  t.start();
  map<int,map<int,Operator>> OpMap;
  const Operator O0ns{g, P0ns{}};
  const Operator O0qg{g, P0qg{}};
  const Operator O0gq{g, P0gq{}};
  for (int nf = 3; nf <= 6; nf++)
    {
      const Operator O0gg{g, P0gg{nf}};
      const Operator O0qgnf = nf * O0qg;
      map<int,Operator> OM;
      OM.insert({EvolutionBasisQCD::PNSP,O0ns});
      OM.insert({EvolutionBasisQCD::PNSM,O0ns});
      OM.insert({EvolutionBasisQCD::PNSV,O0ns});
      OM.insert({EvolutionBasisQCD::PQQ, O0ns});
      OM.insert({EvolutionBasisQCD::PQG, O0qgnf});
      OM.insert({EvolutionBasisQCD::PGQ, O0gq});
      OM.insert({EvolutionBasisQCD::PGG, O0gg});
      OpMap.insert({nf,OM});
    }
  t.stop();

  // Allocate distributions
  cout << "Initializing distributions ..." << endl;
  t.start();
  map<int,Distribution> DistMap;
  for (int i = EvolutionBasisQCD::GLUON; i <= EvolutionBasisQCD::V35; i++)
    DistMap.insert({i,Distribution{g, LHToyPDFs, i}});
  t.stop();

  cout << "Initializing set of operators and distributions ..." << endl;
  t.start();
  // Allocate maps
  map<int,EvolutionBasisQCD> basis;
  for (int nf = 3; nf <= 6; nf++)
    basis.insert({nf,EvolutionBasisQCD{nf}});

  // Allocate set of operators
  map<int,Set<Operator>> Splittings;
  for (int nf = 3; nf <= 6; nf++)
    Splittings.insert({nf,Set<Operator>{basis.at(nf), OpMap.at(nf)}});

  // Allocate set of initial distributions
  Set<Distribution> PDFs{basis.at(5), DistMap};

  t.stop();

  // ===============================================================
  // Test products
  cout << "\nTesting products ..." << endl;
  t.start();
  auto Product = Splittings.at(5) * PDFs;
  cout << "(Splitting * PDFs)[GLUON](x=0.1) = "
       << Product.at(EvolutionBasisQCD::GLUON).Evaluate(0.1) << endl;

  auto Product2 = 2 * Product;
  cout << "(2 * Splitting * PDFs)[GLUON](x=0.1) = "
       << Product2.at(EvolutionBasisQCD::GLUON).Evaluate(0.1) << endl;

  auto Sum = Product.at(EvolutionBasisQCD::GLUON) + Product.at(EvolutionBasisQCD::GLUON);
  cout << "[(Splitting * PDFs)[GLUON] + (Splitting * PDFs)[GLUON]](x=0.1) = "
       << Sum.Evaluate(0.1) << endl;
  t.stop();

  // ===============================================================
  cout << "\nFull computation time:" << endl;
  ttot.stop();

  return 0;
}
