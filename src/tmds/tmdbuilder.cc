//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/grid.h"
#include "apfel/tmdbuilder.h"
#include "apfel/timer.h"
#include "apfel/matchingfunctions.h"
#include "apfel/evolutionbasisqcd.h"
#include "apfel/tools.h"

using namespace std;

namespace apfel {
  /*
   * For all expressions see arXiv:1604.07869.
   */
  //_____________________________________________________________________________
  map<int,TmdObjects> InitializeTmdObjects(Grid           const& g,
					   vector<double> const& Thresholds,
					   double         const& IntEps)
  {
    cout << "Initializing TMD objects for matching and evolution... ";
    Timer t;
    t.start();

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the thresholds
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
	nfi++;

    // ===============================================================
    // LO matching functions operators.
    map<int,Operator> MatchLO;
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};
    MatchLO.insert({EvolutionBasisQCD::PNSP, Id});
    MatchLO.insert({EvolutionBasisQCD::PNSM, Id});
    MatchLO.insert({EvolutionBasisQCD::PNSV, Id});
    MatchLO.insert({EvolutionBasisQCD::PQQ,  Id});
    MatchLO.insert({EvolutionBasisQCD::PQG,  Zero});
    MatchLO.insert({EvolutionBasisQCD::PGQ,  Zero});
    MatchLO.insert({EvolutionBasisQCD::PGG,  Id});

    // ===============================================================
    // NLO matching functions operators.
    map<int,Operator> MatchNLO;
    const Operator O1ns{g, C1ns{}, IntEps};
    const Operator O1qg{g, C1qg{}, IntEps};
    const Operator O1gq{g, C1gq{}, IntEps};
    const Operator O1gg{g, C1gg{}, IntEps};
    MatchNLO.insert({EvolutionBasisQCD::PNSP, O1ns});
    MatchNLO.insert({EvolutionBasisQCD::PNSM, O1ns});
    MatchNLO.insert({EvolutionBasisQCD::PNSV, O1ns});
    MatchNLO.insert({EvolutionBasisQCD::PQQ,  O1ns});
    MatchNLO.insert({EvolutionBasisQCD::PQG,  O1qg});
    MatchNLO.insert({EvolutionBasisQCD::PGQ,  O1gq});
    MatchNLO.insert({EvolutionBasisQCD::PGG,  O1gg});

    // Define object of the structure containing the TmdObjects
    map<int,TmdObjects> TmdObj;

    // Construcuct sets of operators for each perturbative order for
    // the matching functions. Initialize also coefficients of: beta
    // function, GammaCusp, gammaV, and Collins-Soper anomalous
    // dimensions.
    map<int,map<int,Set<Operator>>> MatchingFunctions;
    for (auto nf = nfi; nf <= nff; nf++)
      {
	TmdObjects obj;

	// Beta function
	obj.Beta.insert({0, beta0(nf)});
	obj.Beta.insert({1, beta1(nf)});
	obj.Beta.insert({2, beta2(nf)});

	// GammaCusp
	obj.GammaCuspq.insert({0, CF * GammaCusp0()});
	obj.GammaCuspq.insert({1, CF * GammaCusp1(nf)});
	obj.GammaCuspq.insert({2, CF * GammaCusp2(nf)});
	obj.GammaCuspg.insert({0, CA * GammaCusp0()});
	obj.GammaCuspg.insert({1, CA * GammaCusp1(nf)});
	obj.GammaCuspg.insert({2, CA * GammaCusp2(nf)});

	// GammaV
	obj.GammaVq.insert({0, gammaVq0()});
	obj.GammaVq.insert({1, gammaVq1(nf)});
	obj.GammaVq.insert({2, gammaVq2(nf)});
	obj.GammaVg.insert({0, gammaVg0(nf)});
	obj.GammaVg.insert({1, gammaVg1(nf)});
	obj.GammaVg.insert({2, gammaVg2(nf)});

	// CSd
	valarray<double> d1{CSd10(), CSd11()};
	valarray<double> d2{CSd20(nf), CSd21(nf), CSd22(nf)};
	valarray<double> d3{CSd30(nf), CSd31(nf), CSd32(nf), CSd33(nf)};
	obj.CSdq.insert({0, CF * d1});
	obj.CSdq.insert({1, CF * d2});
	obj.CSdq.insert({2, CF * d3});
	obj.CSdg.insert({0, CA * d1});
	obj.CSdg.insert({1, CA * d2});
	obj.CSdg.insert({2, CA * d3});

	// Matching functions.
	const EvolutionBasisQCD evb{nf};
	obj.MatchingFunctions.insert({0, Set<Operator>{evb, MatchLO}});
	obj.MatchingFunctions.insert({1, Set<Operator>{evb, MatchNLO}});
	//obj.MatchingFunctions.insert({2, Set<Operator>{evb, MatchNNLO}});

	TmdObj.insert({nf, obj});
      }
    t.stop();

    return TmdObj;
  }
}
