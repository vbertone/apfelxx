//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
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

int main()
{
  // Time counter
  apfel::Timer t,ttot;

  // Grid
  const apfel::Grid g{{apfel::SubGrid{80,1e-5,3}, apfel::SubGrid{50,1e-1,5}, apfel::SubGrid{40,8e-1,5}}};

  // ===============================================================
  // Allocate LO splitting functions operators
  std::cout << "Initializing operators... ";
  t.start();
  std::map<int, std::map<int, apfel::Operator>> OpMap;
  const apfel::Operator O0ns{g, apfel::P0ns{}};
  const apfel::Operator O0qg{g, apfel::P0qg{}};
  const apfel::Operator O0gq{g, apfel::P0gq{}};
  for (int nf = 3; nf <= 6; nf++)
    {
      const apfel::Operator O0gg{g, apfel::P0gg{nf}};
      const apfel::Operator O0qgnf = nf * O0qg;
      std::map<int, apfel::Operator> OM;
      OM.insert({apfel::EvolutionBasisQCD::PNSP,O0ns});
      OM.insert({apfel::EvolutionBasisQCD::PNSM,O0ns});
      OM.insert({apfel::EvolutionBasisQCD::PNSV,O0ns});
      OM.insert({apfel::EvolutionBasisQCD::PQQ, O0ns});
      OM.insert({apfel::EvolutionBasisQCD::PQG, O0qgnf});
      OM.insert({apfel::EvolutionBasisQCD::PGQ, O0gq});
      OM.insert({apfel::EvolutionBasisQCD::PGG, O0gg});
      OpMap.insert({nf,OM});
    }
  t.stop();

  // Allocate distributions
  std::cout << "Initializing distributions... ";
  t.start();
  std::map<int, apfel::Distribution> DistMap = DistributionMap(g, apfel::LHToyPDFs, 0);
  t.stop();

  std::cout << "Initializing set of operators and distributions... ";
  t.start();
  // Allocate maps
  std::map<int, apfel::EvolutionBasisQCD> basis;
  for (int nf = 3; nf <= 6; nf++)
    basis.insert({nf, apfel::EvolutionBasisQCD{nf}});

  // Allocate set of operators
  std::map<int, apfel::Set<apfel::Operator>> Splittings;
  for (int nf = 3; nf <= 6; nf++)
    Splittings.insert({nf, apfel::Set<apfel::Operator>{basis.at(nf), OpMap.at(nf)}});

  // Allocate set of initial distributions
  apfel::Set<apfel::Distribution> PDFs{basis.at(5), DistMap};
  t.stop();

  // Test products
  std::cout << "\nTesting products..." << std::endl;
  t.start();
  const apfel::Set<apfel::Distribution> Product = Splittings.at(5) * PDFs;
  std::cout << "(Splitting * PDFs)[GLUON](x=0.1) = "
       << Product.at(apfel::EvolutionBasisQCD::GLUON).Evaluate(0.1) << std::endl;

  const apfel::Set<apfel::Distribution> Product2 = 3 * Splittings.at(5) * PDFs;
  std::cout << "(3 * Splitting * PDFs)[GLUON](x=0.1) = "
       << Product2.at(apfel::EvolutionBasisQCD::GLUON).Evaluate(0.1) << std::endl;

  const apfel::Set<apfel::Distribution> Sum = ( Splittings.at(5) + 2 * Splittings.at(5) ) * PDFs;
  std::cout << "[(Splitting * PDFs)[GLUON] + 2 * (Splitting * PDFs)[GLUON]](x=0.1) = "
       << Sum.at(apfel::EvolutionBasisQCD::GLUON).Evaluate(0.1) << std::endl;
  t.stop();

  std::cout << "\nTotal ";
  ttot.stop();

  return 0;
}
