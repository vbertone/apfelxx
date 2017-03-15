//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/grid.h"
#include "apfel/dglapbuilder.h"
#include "apfel/dglap.h"
#include "apfel/operator.h"
#include "apfel/set.h"
#include "apfel/timer.h"
#include "apfel/tools.h"
#include "apfel/evolutionbasisqcd.h"
#include "apfel/matchingbasisqcd.h"
#include "apfel/splittingfunctions.h"
#include "apfel/matchingconditions.h"

#include <map>

using namespace std;

namespace apfel {

  /**
   * @brief The PDF class
   *
   * Helper class for the construction of PDFs from function
   */
  class PDF: public Distribution
      {
        public:
        PDF(Grid                                        const& g,
            function<double(int const&, double const&)> const& InPDFsFunc,
            int                                         const& ipdf):
          Distribution(g)
        {
          for (auto const& ix: _grid.GetJointGrid().GetGrid())
            if (ix < 1)
              _distributionJointGrid.push_back(InPDFsFunc(ipdf,ix));
            else
              _distributionJointGrid.push_back(0);

          for (auto ig=0; ig<_grid.nGrids(); ig++)
            {
              vector<double> sg;
              for (auto const& ix: _grid.GetSubGrid(ig).GetGrid())
                if (ix < 1)
                  sg.push_back(InPDFsFunc(ipdf,ix));
                else
                  sg.push_back(0);
              _distributionSubGrid.push_back(sg);
            }
        }
  };

  //_____________________________________________________________________________
  Dglap DglapBuildQCD(Grid                                        const& g,
                      function<double(int const&, double const&)> const& InPDFsFunc,
                      double                                      const& MuRef,
                      vector<double>                              const& Masses,
                      vector<double>                              const& Thresholds,
                      int                                         const& PerturbativeOrder,
                      function<double(double const&)>             const& Alphas,
                      double                                      const& IntEps,
                      int                                         const& nsteps)
  {
    cout << "Initialization... ";
    Timer t;
    t.start();

    // Compute initial and final number of active flavours 
    // according to the vector of thresholds (it assumes that
    // the thresholds vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if ( v <= 0 )
	nfi++;

    // Compute AlphaQCD above and below the thresholds
    unordered_map<int,double> asThUp;
    unordered_map<int,double> asThDown;
    for (auto nf = nfi + 1; nf <= nff; nf++)
      {
	asThDown.insert({nf,Alphas(Thresholds[nf-1])/FourPi});
	asThUp.insert({nf,Alphas(Thresholds[nf-1]+eps8)/FourPi});
      }

    // Allocate convolution maps for the evolution and matching
    unordered_map<int,EvolutionBasisQCD> evbasis;
    unordered_map<int,MatchingBasisQCD>  matchbasis;
    for (int nf = nfi; nf <= nff; nf++)
      {
	evbasis.insert({nf,EvolutionBasisQCD{nf}});
	matchbasis.insert({nf,MatchingBasisQCD{nf}});
      }

    // Allocate initial scale distributions
    unordered_map<int,Distribution> DistMap;
    for (int i = EvolutionBasisQCD::GLUON; i <= EvolutionBasisQCD::V35; i++)
      DistMap.insert({i,PDF{g, InPDFsFunc, i}});

    // Compute numeber of active flavouts the the PDF initial scale and 
    // allocate set of initial distributions.
    int nf0 = 0;
    for (auto const& v : Thresholds)
      if (MuRef > v)
	nf0++;
      else
	break;

    // Create set of initial distributions
    // (assumed to be in the QCD evolution basis).
    Set<Distribution> InPDFs{evbasis.at(nf0), DistMap};

    // Allocate needed operators (matching conditions and splitting functions).
    // By now the code is fast enough to precompute everything at all available
    // perturbative orders and the current perturbative order is accounted for
    // only when the actual splitting functions and matching conditions (lambda)
    // functions are defined.
    // ===============================================================
    // LO Matching conditions
    unordered_map<int,Operator> MatchLO;
    const Operator Id{g, Identity{}, IntEps};
    const Operator Zero{g, Null{}, IntEps};
    MatchLO.insert({MatchingBasisQCD::PNSP, Id});
    MatchLO.insert({MatchingBasisQCD::PNSM, Id});
    MatchLO.insert({MatchingBasisQCD::PNSV, Id});
    MatchLO.insert({MatchingBasisQCD::PQQ,  Id});
    MatchLO.insert({MatchingBasisQCD::PQG,  Zero});
    MatchLO.insert({MatchingBasisQCD::PGQ,  Zero});
    MatchLO.insert({MatchingBasisQCD::PGG,  Id});
    for (int i = MatchingBasisQCD::PT3Q; i <= MatchingBasisQCD::PT35Q; i++)
      MatchLO.insert({i, Id});
    for (int i = MatchingBasisQCD::PT3G; i <= MatchingBasisQCD::PT35G; i++)
      MatchLO.insert({i, Zero});

    // ===============================================================
    // LO splitting functions operators
    unordered_map<int,unordered_map<int,Operator>> OpMapLO;
    const Operator O0ns{g, P0ns{}, IntEps};
    const Operator O0qg{g, P0qg{}, IntEps};
    const Operator O0gq{g, P0gq{}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
	const Operator O0gg{g, P0gg{nf}, IntEps};
	const Operator O0qgnf = nf * O0qg;
	unordered_map<int,Operator> OM;
	OM.insert({EvolutionBasisQCD::PNSP, O0ns});
	OM.insert({EvolutionBasisQCD::PNSM, O0ns});
	OM.insert({EvolutionBasisQCD::PNSV, O0ns});
	OM.insert({EvolutionBasisQCD::PQQ,  O0ns});
	OM.insert({EvolutionBasisQCD::PQG,  O0qgnf});
	OM.insert({EvolutionBasisQCD::PGQ,  O0gq});
	OM.insert({EvolutionBasisQCD::PGG,  O0gg});
	OpMapLO.insert({nf,OM});
      }

    // ===============================================================
    // NLO splitting functions operators
    unordered_map<int,unordered_map<int,Operator>> OpMapNLO;
    for (int nf = nfi; nf <= nff; nf++)
      {
	const Operator O1nsp{g, P1nsp{nf}, IntEps};
	const Operator O1nsm{g, P1nsm{nf}, IntEps};
	const Operator O1ps{g, P1ps{nf}, IntEps};
	const Operator O1qq = O1nsp + O1ps;
	const Operator O1qg{g, P1qg{nf}, IntEps};
	const Operator O1gq{g, P1gq{nf}, IntEps};
	const Operator O1gg{g, P1gg{nf}, IntEps};
	unordered_map<int,Operator> OM;
	OM.insert({EvolutionBasisQCD::PNSP, O1nsp});
	OM.insert({EvolutionBasisQCD::PNSM, O1nsm});
	OM.insert({EvolutionBasisQCD::PNSV, O1nsm});
	OM.insert({EvolutionBasisQCD::PQQ,  O1qq});
	OM.insert({EvolutionBasisQCD::PQG,  O1qg});
	OM.insert({EvolutionBasisQCD::PGQ,  O1gq});
	OM.insert({EvolutionBasisQCD::PGG,  O1gg});
	OpMapNLO.insert({nf,OM});
      }

    // ===============================================================
    // Allocate NNLO Matching conditions
    unordered_map<int,unordered_map<int,Operator>> MatchNNLO;  
    const Operator APS2Hq{g, APS2Hq_0{}, IntEps};
    const Operator ANS2qqH{g, ANS2qqH_0{}, IntEps};
    const Operator AS2qqH = ANS2qqH + APS2Hq;
    const Operator AS2Hg{g, AS2Hg_0{}, IntEps};
    const Operator AS2gqH{g, AS2gqH_0{}, IntEps};
    const Operator AS2ggH{g, AS2ggH_0{}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
	unordered_map<int,Operator> OM;
	OM.insert({MatchingBasisQCD::PNSP, ANS2qqH});
	OM.insert({MatchingBasisQCD::PNSM, ANS2qqH});
	OM.insert({MatchingBasisQCD::PNSV, ANS2qqH});
	OM.insert({MatchingBasisQCD::PQQ,  AS2qqH});
	OM.insert({MatchingBasisQCD::PQG,  AS2Hg});
	OM.insert({MatchingBasisQCD::PGQ,  AS2gqH});
	OM.insert({MatchingBasisQCD::PGG,  AS2ggH});
	const Operator AS2TqH = ANS2qqH - nf * APS2Hq;
	const Operator AS2Tg  = - nf * AS2Hg;
	for (int i = MatchingBasisQCD::PT3Q; i <= MatchingBasisQCD::PT35Q; i++)
	  if (i > MatchingBasisQCD::PT3Q + nf - 1)
	    OM.insert({i, AS2qqH});
	  else if (i == MatchingBasisQCD::PT3Q + nf - 1)
	    OM.insert({i, AS2TqH});
	  else
	    OM.insert({i, Zero});
	for (int i = MatchingBasisQCD::PT3G; i <= MatchingBasisQCD::PT35G; i++)
	  if (i > MatchingBasisQCD::PT3G + nf - 1)
	    OM.insert({i, AS2Hg});
	  else if (i == MatchingBasisQCD::PT3G + nf - 1)
	    OM.insert({i, AS2Tg});
	  else
	    OM.insert({i, Zero});
	MatchNNLO.insert({nf,OM});
      }

    // ===============================================================
    // Allocate NNLO splitting functions operators
    unordered_map<int,unordered_map<int,Operator>> OpMapNNLO;
    for (int nf = nfi; nf <= nff; nf++)
      {
	const Operator O2nsp{g, P2nsp{nf}, IntEps};
	const Operator O2nsm{g, P2nsm{nf}, IntEps};
	const Operator O2nss{g, P2nss{nf}, IntEps};
	const Operator O2nsv = O2nsm + O2nss;
	const Operator O2ps{g, P2ps{nf}, IntEps};
	const Operator O2qq = O2nsp + O2ps;
	const Operator O2qg{g, P2qg{nf}, IntEps};
	const Operator O2gq{g, P2gq{nf}, IntEps};
	const Operator O2gg{g, P2gg{nf}, IntEps};
	unordered_map<int,Operator> OM;
	OM.insert({EvolutionBasisQCD::PNSP, O2nsp});
	OM.insert({EvolutionBasisQCD::PNSM, O2nsm});
	OM.insert({EvolutionBasisQCD::PNSV, O2nsv});
	OM.insert({EvolutionBasisQCD::PQQ,  O2qq});
	OM.insert({EvolutionBasisQCD::PQG,  O2qg});
	OM.insert({EvolutionBasisQCD::PGQ,  O2gq});
	OM.insert({EvolutionBasisQCD::PGG,  O2gg});
	OpMapNNLO.insert({nf,OM});
      }

    // Allocate set of operators
    unordered_map<int,Set<Operator>> P0;
    unordered_map<int,Set<Operator>> P1;
    unordered_map<int,Set<Operator>> P2;
    unordered_map<int,Set<Operator>> M0;
    unordered_map<int,Set<Operator>> M2;
    for (int nf = nfi; nf <= nff; nf++)
      {
	P0.insert({nf,Set<Operator>{evbasis.at(nf), OpMapLO.at(nf)}});
	P1.insert({nf,Set<Operator>{evbasis.at(nf), OpMapNLO.at(nf)}});
	P2.insert({nf,Set<Operator>{evbasis.at(nf), OpMapNNLO.at(nf)}});
	M0.insert({nf,Set<Operator>{matchbasis.at(nf), MatchLO}});
	M2.insert({nf,Set<Operator>{matchbasis.at(nf), MatchNNLO.at(nf)}});
      }

    // Create splitting functions and matching conditions lambda functions
    // according to the requested perturbative order.
    function<Set<Operator>(int const&, double const&)>              SplittingFunctions;
    function<Set<Operator>(bool const&, int const&, double const&)> MatchingConditions;
    if (PerturbativeOrder == 0)
      {
        SplittingFunctions = [=] (int const& nf, double const& mu) -> Set<Operator>
	  { const auto cp = Alphas(mu)/FourPi; return cp * P0.at(nf); };
	MatchingConditions = [=] (bool const&, int const& nf, double const&) -> Set<Operator>
	  { return M0.at(nf); };
      }
    else if (PerturbativeOrder == 1)
      {
        SplittingFunctions = [=] (int const& nf, double const& mu) -> Set<Operator>
	  { const auto cp = Alphas(mu)/FourPi; return cp * ( P0.at(nf) + cp * P1.at(nf) ); };
	MatchingConditions = [=] (bool const&, int const& nf, double const&) -> Set<Operator>
	  { return M0.at(nf); };
      }
    else if (PerturbativeOrder == 2)
      {
        SplittingFunctions = [=] (int const& nf, double const& mu) -> Set<Operator>
	  { const auto cp = Alphas(mu)/FourPi; return cp * ( P0.at(nf) + cp * ( P1.at(nf) + cp * P2.at(nf) ) ); };
	MatchingConditions = [=] (bool const& Up, int const& nf, double const&) -> Set<Operator>
	  { const auto cp = asThUp.at(nf+1); return M0.at(nf) + ( Up ? 1 : -1) * cp * cp * M2.at(nf); };
      }

    t.stop();

    // Initialize DGLAP evolution
    return Dglap{SplittingFunctions, MatchingConditions, InPDFs, MuRef, Masses, Thresholds, nsteps};
  }
}
