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

  // ===============================================================
  // Allocate LO Matching conditions
  unordered_map<int,Operator> MatchLO;
  const Operator Id{g, Identity{}};
  const Operator Zero{g, Null{}};
  MatchLO.insert({EvolutionBasis::PNSP,Id});
  MatchLO.insert({EvolutionBasis::PNSM,Id});
  MatchLO.insert({EvolutionBasis::PNSV,Id});
  MatchLO.insert({EvolutionBasis::PQQ, Id});
  MatchLO.insert({EvolutionBasis::PQG, Zero});
  MatchLO.insert({EvolutionBasis::PGQ, Zero});
  MatchLO.insert({EvolutionBasis::PGG, Id});
  MatchLO.insert({EvolutionBasis::PT3Q, Id});
  MatchLO.insert({EvolutionBasis::PT3G, Zero});
  MatchLO.insert({EvolutionBasis::PT8Q, Id});
  MatchLO.insert({EvolutionBasis::PT8G, Zero});
  MatchLO.insert({EvolutionBasis::PT15Q,Id});
  MatchLO.insert({EvolutionBasis::PT15G,Zero});
  MatchLO.insert({EvolutionBasis::PT24Q,Id});
  MatchLO.insert({EvolutionBasis::PT24G,Zero});
  MatchLO.insert({EvolutionBasis::PT35Q,Id});
  MatchLO.insert({EvolutionBasis::PT35G,Zero});

  // ===============================================================
  // Allocate NNLO Matching conditions
  unordered_map<int,unordered_map<int,Operator>> MatchNNLO;
  const Operator APS2Hq{g, APS2Hq_0{}};
  const Operator ANS2qqH{g, ANS2qqH_0{}};
  const Operator AS2qqH = ANS2qqH + APS2Hq;
  const Operator AS2Hg{g, AS2Hg_0{}};
  const Operator AS2gqH{g, AS2gqH_0{}};
  const Operator AS2ggH{g, AS2ggH_0{}};
  for (int nf = 3; nf <= 6; nf++)
    {
      unordered_map<int,Operator> OM;
      OM.insert({EvolutionBasis::PNSP, ANS2qqH});
      OM.insert({EvolutionBasis::PNSM, ANS2qqH});
      OM.insert({EvolutionBasis::PNSV, ANS2qqH});
      OM.insert({EvolutionBasis::PQQ,  AS2qqH});
      OM.insert({EvolutionBasis::PQG,  AS2Hg});
      OM.insert({EvolutionBasis::PGQ,  AS2gqH});
      OM.insert({EvolutionBasis::PGG,  AS2ggH});
      OM.insert({EvolutionBasis::PT3Q, ANS2qqH});
      OM.insert({EvolutionBasis::PT3G, Zero});
      OM.insert({EvolutionBasis::PT8Q, ANS2qqH});
      OM.insert({EvolutionBasis::PT8G, Zero});
      if ( nf == 3 )
	{
	  OM.insert({EvolutionBasis::PT15Q,ANS2qqH-3*APS2Hq});
	  OM.insert({EvolutionBasis::PT15G,-3*AS2Hg});
	  OM.insert({EvolutionBasis::PT24Q,AS2qqH});
	  OM.insert({EvolutionBasis::PT24G,AS2Hg});
	  OM.insert({EvolutionBasis::PT35Q,AS2qqH});
	  OM.insert({EvolutionBasis::PT35G,AS2Hg});
	}
      else if ( nf == 4 )
	{
	  OM.insert({EvolutionBasis::PT15Q,ANS2qqH});
	  OM.insert({EvolutionBasis::PT15G,Zero});
	  OM.insert({EvolutionBasis::PT24Q,ANS2qqH-4*APS2Hq});
	  OM.insert({EvolutionBasis::PT24G,-4*AS2Hg});
	  OM.insert({EvolutionBasis::PT35Q,AS2qqH});
	  OM.insert({EvolutionBasis::PT35G,AS2Hg});
	}
      else if ( nf == 5 )
	{
	  OM.insert({EvolutionBasis::PT15Q,ANS2qqH});
	  OM.insert({EvolutionBasis::PT15G,Zero});
	  OM.insert({EvolutionBasis::PT24Q,ANS2qqH});
	  OM.insert({EvolutionBasis::PT24G,Zero});
	  OM.insert({EvolutionBasis::PT35Q,ANS2qqH-5*APS2Hq});
	  OM.insert({EvolutionBasis::PT35G,-5*AS2Hg});
	}
      else
	{
	  OM.insert({EvolutionBasis::PT15Q,ANS2qqH});
	  OM.insert({EvolutionBasis::PT15G,Zero});
	  OM.insert({EvolutionBasis::PT24Q,ANS2qqH});
	  OM.insert({EvolutionBasis::PT24G,Zero});
	  OM.insert({EvolutionBasis::PT35Q,ANS2qqH});
	  OM.insert({EvolutionBasis::PT35G,Zero});
	}
      MatchNNLO.insert({nf,OM});
    }

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
      OM.insert({EvolutionBasis::PNSP, O0ns});
      OM.insert({EvolutionBasis::PNSM, O0ns});
      OM.insert({EvolutionBasis::PNSV, O0ns});
      OM.insert({EvolutionBasis::PQQ,  O0ns});
      OM.insert({EvolutionBasis::PQG,  O0qgnf});
      OM.insert({EvolutionBasis::PGQ,  O0gq});
      OM.insert({EvolutionBasis::PGG,  O0gg});
      OM.insert({EvolutionBasis::PT3Q, O0ns});
      OM.insert({EvolutionBasis::PT3G, O0qgnf});
      OM.insert({EvolutionBasis::PT8Q, O0ns});
      OM.insert({EvolutionBasis::PT8G, O0qgnf});
      OM.insert({EvolutionBasis::PT15Q,O0ns});
      OM.insert({EvolutionBasis::PT15G,O0qgnf});
      OM.insert({EvolutionBasis::PT24Q,O0ns});
      OM.insert({EvolutionBasis::PT24G,O0qgnf});
      OM.insert({EvolutionBasis::PT35Q,O0ns});
      OM.insert({EvolutionBasis::PT35G,O0qgnf});
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
      OM.insert({EvolutionBasis::PNSP, O1nsp});
      OM.insert({EvolutionBasis::PNSM, O1nsm});
      OM.insert({EvolutionBasis::PNSV, O1nsm});
      OM.insert({EvolutionBasis::PQQ,  O1qq});
      OM.insert({EvolutionBasis::PQG,  O1qg});
      OM.insert({EvolutionBasis::PGQ,  O1gq});
      OM.insert({EvolutionBasis::PGG,  O1gg});
      OM.insert({EvolutionBasis::PT3Q, O1qq});
      OM.insert({EvolutionBasis::PT3G, O1qg});
      OM.insert({EvolutionBasis::PT8Q, O1qq});
      OM.insert({EvolutionBasis::PT8G, O1qg});
      OM.insert({EvolutionBasis::PT15Q,O1qq});
      OM.insert({EvolutionBasis::PT15G,O1qg});
      OM.insert({EvolutionBasis::PT24Q,O1qq});
      OM.insert({EvolutionBasis::PT24G,O1qg});
      OM.insert({EvolutionBasis::PT35Q,O1qq});
      OM.insert({EvolutionBasis::PT35G,O1qg});
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
      OM.insert({EvolutionBasis::PNSP, O2nsp});
      OM.insert({EvolutionBasis::PNSM, O2nsm});
      OM.insert({EvolutionBasis::PNSV, O2nsv});
      OM.insert({EvolutionBasis::PQQ,  O2qq});
      OM.insert({EvolutionBasis::PQG,  O2qg});
      OM.insert({EvolutionBasis::PGQ,  O2gq});
      OM.insert({EvolutionBasis::PGG,  O2gg});
      OM.insert({EvolutionBasis::PT3Q, O2qq});
      OM.insert({EvolutionBasis::PT3G, O2qg});
      OM.insert({EvolutionBasis::PT8Q, O2qq});
      OM.insert({EvolutionBasis::PT8G, O2qg});
      OM.insert({EvolutionBasis::PT15Q,O2qq});
      OM.insert({EvolutionBasis::PT15G,O2qg});
      OM.insert({EvolutionBasis::PT24Q,O2qq});
      OM.insert({EvolutionBasis::PT24G,O2qg});
      OM.insert({EvolutionBasis::PT35Q,O2qq});
      OM.insert({EvolutionBasis::PT35G,O2qg});
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
  unordered_map<int,Set<Operator>> MatchingLO;
  unordered_map<int,Set<Operator>> MatchingNNLO;
  for (int nf = 3; nf <= 6; nf++)
    {
      SplittingsLO.insert({nf,Set<Operator>{basis.at(nf), OpMapLO.at(nf)}});
      SplittingsNLO.insert({nf,Set<Operator>{basis.at(nf), OpMapNLO.at(nf)}});
      SplittingsNNLO.insert({nf,Set<Operator>{basis.at(nf), OpMapNNLO.at(nf)}});
      MatchingLO.insert({nf,Set<Operator>{basis.at(nf), MatchLO}});
      MatchingNNLO.insert({nf,Set<Operator>{basis.at(nf), MatchNNLO.at(nf)}});
    }

  // Allocate set of initial distributions
  Set<Distribution> PDFs{basis.at(3), DistMap};

  // Coupling
  AlphaQCD as{0.35, sqrt(2), {0, 0, 0, sqrt(2), 4.5, 175}, 2};

  // Compute alpha_s right above the thresholds
  const auto eps = 1e-8;
  const double asth[] = {as.GetObject(sqrt(2)+eps),as.GetObject(4.5+eps),as.GetObject(175+eps)};

  int nsteps = 10;
  DGLAP evolution{
    [&] (int const& nf, double const& mu) -> Set<Operator>{
      const auto cp = as.GetObject(mu) / FourPi;
      const auto LO   = cp * SplittingsLO.at(nf);
      const auto NLO  = cp * cp * SplittingsNLO.at(nf);
      const auto NNLO = cp * cp * cp * SplittingsNNLO.at(nf);
      return LO + NLO + NNLO;
    },
      [&] (bool const&, int const& nf, double const&) -> Set<Operator>{
	const auto cp = asth[nf-3] / FourPi;
	const auto LO = MatchingLO.at(nf);
	const auto NNLO = cp * cp * MatchingNNLO.at(nf);
	return LO + NNLO;
      },
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
