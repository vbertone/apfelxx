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
#include <apfel/timer.h>
#include <apfel/tools.h>
#include <apfel/alphaqcd.h>
#include <apfel/set.h>
#include <apfel/dglap.h>
#include <apfel/evolutionbasis.h>
#include <apfel/splittingfunctions.h>
#include <apfel/matchingconditions.h>

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
  if      (i == EvolutionBasis::GLUON    ) return xglu(x);
  // Singlet, T15, T24, T35
  else if (i == EvolutionBasis::SIGMA   ||
	   i == EvolutionBasis::T15     ||
	   i == EvolutionBasis::T24     ||
	   i == EvolutionBasis::T35      ) return xdnv(x) + 2 * xdbar(x) + xupv(x) + 2 * xubar(x) + 2 * xsbar(x);
  // T3
  else if (i == EvolutionBasis::T3       ) return xupv(x) + 2 * xubar(x) - xdnv(x) - 2 * xdbar(x);
  // T8
  else if (i == EvolutionBasis::T8       ) return xupv(x) + 2 * xubar(x) + xdnv(x) + 2 * xdbar(x) - 4 * xsbar(x);
  // Valence, V8, V15, V24, V35
  else if (i == EvolutionBasis::VALENCE ||
	   i == EvolutionBasis::V8      ||
	   i == EvolutionBasis::V15     ||
	   i == EvolutionBasis::V24     ||
	   i == EvolutionBasis::V35      ) return xupv(x) + xdnv(x);
  // V3
  else if (i == EvolutionBasis::V3       )  return xupv(x) - xdnv(x);
  else              return 0;
}

// The PDF class
class PDF: public Distribution
{
public:
  // Standard constructor
  PDF(Grid const& gr, function<double(int,double)> const& inPDFs, int const& ipdf): Distribution(gr)
  {
    for (auto const& ix: _grid.GetJointGrid().GetGrid())
      if (ix < 1) _distributionJointGrid.push_back(inPDFs(ipdf,ix));
      else        _distributionJointGrid.push_back(0);

    for (auto ig=0; ig<_grid.nGrids(); ig++)
      {
        vector<double> sg;
        for (auto const& ix: _grid.GetSubGrid(ig).GetGrid())
          if (ix < 1) sg.push_back(inPDFs(ipdf,ix));
          else        sg.push_back(0);
        _distributionSubGrid.push_back(sg);
      }
  }
};

int main()
{
  Timer t;

  cout << "Initialization ..." << endl;
  t.start();
  // Allocate grid
  const Grid g{{SubGrid{80,1e-5,3}, SubGrid{50,1e-1,3}, SubGrid{40,8e-1,3}}};

  // LO Matching conditions
  const Operator Id{g, Identity{}};
  const Operator Zero{g, Null{}};
  unordered_map<int,Operator> Match;
  Match.insert({EvolutionBasis::PNSP,Id});
  Match.insert({EvolutionBasis::PNSM,Id});
  Match.insert({EvolutionBasis::PNSV,Id});
  Match.insert({EvolutionBasis::PQQ, Id});
  Match.insert({EvolutionBasis::PQG, Zero});
  Match.insert({EvolutionBasis::PGQ, Zero});
  Match.insert({EvolutionBasis::PGG, Id});

  // ===============================================================
  // Allocate LO splitting functions operators
  unordered_map<int,unordered_map<int,Operator>> OpMapLO;
  const Operator O0ns{g, P0ns{}};
  const Operator O0qg{g, P0qg{}};
  const Operator O0gq{g, P0gq{}};
  for (int nf = 3; nf <= 6; nf++)
    {
      const Operator O0gg{g, P0gg{nf}};
      const Operator O0qgnf = nf * O0qg;
      unordered_map<int,Operator> OM;
      OM.insert({EvolutionBasis::PNSP,O0ns});
      OM.insert({EvolutionBasis::PNSM,O0ns});
      OM.insert({EvolutionBasis::PNSV,O0ns});
      OM.insert({EvolutionBasis::PQQ, O0ns});
      OM.insert({EvolutionBasis::PQG, O0qgnf});
      OM.insert({EvolutionBasis::PGQ, O0gq});
      OM.insert({EvolutionBasis::PGG, O0gg});
      OpMapLO.insert({nf,OM});
    }

  // ===============================================================
  // Allocate NLO splitting functions operators
  unordered_map<int,unordered_map<int,Operator>> OpMapNLO;
  for (int nf = 3; nf <= 6; nf++)
    {
      const Operator O1nsp{g, P1nsp{nf}};
      const Operator O1nsm{g, P1nsm{nf}};
      const Operator O1ps{g, P1ps{nf}};
      const Operator O1qq = O1nsp + O1ps;
      const Operator O1qg{g, P1qg{nf}};
      const Operator O1gq{g, P1gq{nf}};
      const Operator O1gg{g, P1gg{nf}};
      unordered_map<int,Operator> OM;
      OM.insert({EvolutionBasis::PNSP,O1nsp});
      OM.insert({EvolutionBasis::PNSM,O1nsm});
      OM.insert({EvolutionBasis::PNSV,O1nsm});
      OM.insert({EvolutionBasis::PQQ, O1qq});
      OM.insert({EvolutionBasis::PQG, O1qg});
      OM.insert({EvolutionBasis::PGQ, O1gq});
      OM.insert({EvolutionBasis::PGG, O1gg});
      OpMapNLO.insert({nf,OM});
    }

  // ===============================================================
  // Allocate NNLO splitting functions operators
  unordered_map<int,unordered_map<int,Operator>> OpMapNNLO;
  for (int nf = 3; nf <= 6; nf++)
    {
      const Operator O2nsp{g, P2nsp{nf}};
      const Operator O2nsm{g, P2nsm{nf}};
      const Operator O2nss{g, P2nss{nf}};
      const Operator O2nsv = O2nsm + O2nss;
      const Operator O2ps{g, P2ps{nf}};
      const Operator O2qq = O2nsp + O2ps;
      const Operator O2qg{g, P2qg{nf}};
      const Operator O2gq{g, P2gq{nf}};
      const Operator O2gg{g, P2gg{nf}};
      unordered_map<int,Operator> OM;
      OM.insert({EvolutionBasis::PNSP,O2nsp});
      OM.insert({EvolutionBasis::PNSM,O2nsm});
      OM.insert({EvolutionBasis::PNSV,O2nsv});
      OM.insert({EvolutionBasis::PQQ, O2qq});
      OM.insert({EvolutionBasis::PQG, O2qg});
      OM.insert({EvolutionBasis::PGQ, O2gq});
      OM.insert({EvolutionBasis::PGG, O2gg});
      OpMapNNLO.insert({nf,OM});
    }

  // Allocate distributions
  unordered_map<int,Distribution> DistMap;
  for (int i = EvolutionBasis::GLUON; i <= EvolutionBasis::V35; i++)
    DistMap.insert({i,PDF{g, LHToyPDFs, i}});

  // Allocate maps
  unordered_map<int,EvolutionBasis> basis;
  for (int nf = 3; nf <= 6; nf++)
    basis.insert({nf,EvolutionBasis{nf}});

  // Allocate set of operators
  unordered_map<int,Set<Operator>> SplittingsLO;
  unordered_map<int,Set<Operator>> SplittingsNLO;
  unordered_map<int,Set<Operator>> SplittingsNNLO;
  unordered_map<int,Set<Operator>> Matching;
  for (int nf = 3; nf <= 6; nf++)
    {
      SplittingsLO.insert({nf,Set<Operator>{basis.at(nf), OpMapLO.at(nf)}});
      SplittingsNLO.insert({nf,Set<Operator>{basis.at(nf), OpMapNLO.at(nf)}});
      SplittingsNNLO.insert({nf,Set<Operator>{basis.at(nf), OpMapNNLO.at(nf)}});
      Matching.insert({nf,Set<Operator>{basis.at(nf), Match}});
    }

  // Allocate set of initial distributions
  Set<Distribution> PDFs{basis.at(3), DistMap};

  // Coupling
  AlphaQCD as{0.35, sqrt(2), {0, 0, 0, sqrt(2), 4.5, 175}, 1};

  int nsteps = 10;
  DGLAP evolution{
    [&] (int const& nf, double const& mu) -> Set<Operator> {
      const double cp = as.GetObject(mu) / FourPi;
      auto LO   = cp * SplittingsLO.at(nf);
      auto NLO  = cp * cp * SplittingsNLO.at(nf);
      auto NNLO = cp * cp * cp * SplittingsNNLO.at(nf);
      return LO + NLO + NNLO;
    },
      [&] (bool const&, int const& nf, double const&) -> Set<Operator>{ return Matching.at(nf); },
	PDFs, sqrt(2), {0, 0, 0, sqrt(2), 4.5, 175}, nsteps};
  t.printTime(t.stop());

  double Q = 100;
  cout << scientific;
  cout << "\nEvolution (4th order Runge-Kutta with " << nsteps << " steps) to Q = " << Q << " GeV ..." << endl;
  t.start();
  auto pdfs = evolution.GetObject(Q);
  t.printTime(t.stop());

  double xlha[] = {1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2,
		   1e-1, 3e-1, 5e-1, 7e-1, 9e-1};

  cout << "\nalpha_QCD(Q) = " << as.GetObject(Q) << endl;
  cout << "Standard evolution:" << endl;
  cout << "   x    "
       << "   u-ubar   "
       << "   d-dbar   "
       << " 2(ubr+dbr) "
       << "   c+cbar   "
       << "    gluon   "
       << endl;

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

  return 0;
}
