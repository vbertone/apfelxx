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
#include <apfel/lhtoypdfs.h>

#include <cmath>
#include <map>
#include <functional>

using namespace apfel;
using namespace std;

int main()
{
  // Time counter
  Timer t,ttot;

  // Grid
  const Grid g{{SubGrid{80,1e-5,3}, SubGrid{50,1e-1,5}, SubGrid{40,8e-1,5}}};

  // ===============================================================
  // Allocate LO splitting functions operators
  cout << "Initializing operators... ";
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
  cout << "Initializing distributions... ";
  t.start();
  map<int,Distribution> DistMap = DistributionMap(g, LHToyPDFs, 0);
  t.stop();

  cout << "Initializing set of operators and distributions... ";
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

  // Test products
  cout << "\nTesting products..." << endl;
  t.start();
  auto Product = Splittings.at(5) * PDFs;
  cout << "(Splitting * PDFs)[GLUON](x=0.1) = "
       << Product.at(EvolutionBasisQCD::GLUON).Evaluate(0.1) << endl;

  auto Product2 = 3 * Splittings.at(5) * PDFs;
  cout << "(3 * Splitting * PDFs)[GLUON](x=0.1) = "
       << Product2.at(EvolutionBasisQCD::GLUON).Evaluate(0.1) << endl;

  auto Sum = ( Splittings.at(5) + 2 * Splittings.at(5) ) * PDFs;
  cout << "[(Splitting * PDFs)[GLUON] + 2 * (Splitting * PDFs)[GLUON]](x=0.1) = "
       << Sum.at(EvolutionBasisQCD::GLUON).Evaluate(0.1) << endl;
  t.stop();

  cout << "\nTotal ";
  ttot.stop();

  return 0;
}
