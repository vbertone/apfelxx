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
#include "apfel/betaqcd.h"
#include "apfel/gammacusp.h"
#include "apfel/gammav.h"
#include "apfel/gammacs.h"
#include "apfel/logzeta.h"
#include "apfel/tools.h"
#include "apfel/constants.h"

using namespace std;

namespace apfel {
  //_____________________________________________________________________________
  map<int,TmdObjects> InitializeTmdObjects(Grid           const& g,
					   vector<double> const& Thresholds,
					   double         const& IntEps)
  {
    report("Initializing TMD objects for matching and evolution... ");
    Timer t;

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
    map<int,map<int,Operator>> MatchPDFsNLO;
    const Operator O1ns{g, C1ns{}, IntEps};
    const Operator O1qg{g, C1qg{}, IntEps};
    const Operator O1gq{g, C1gq{}, IntEps};
    const Operator O1gg{g, C1gg{}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
	const Operator O1qgnf = nf * O1qg;
	map<int,Operator> OM;
	OM.insert({EvolutionBasisQCD::PNSP, O1ns});
	OM.insert({EvolutionBasisQCD::PNSM, O1ns});
	OM.insert({EvolutionBasisQCD::PNSV, O1ns});
	OM.insert({EvolutionBasisQCD::PQQ,  O1ns});
	OM.insert({EvolutionBasisQCD::PQG,  O1qgnf});
	OM.insert({EvolutionBasisQCD::PGQ,  O1gq});
	OM.insert({EvolutionBasisQCD::PGG,  O1gg});
	MatchPDFsNLO.insert({nf,OM});
      }

    // ===============================================================
    // NNLO matching functions operators.
    map<int,map<int,Operator>> MatchPDFsNNLO;
    const Operator O2Vqqb{g, C2Vqqb{}, IntEps};
    const Operator O2ps{g, C2ps{}, IntEps};
    const Operator O2qg{g, C2qg{}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
	const Operator O2Vqq{g, C2Vqq{nf}, IntEps};
	const Operator O2qgnf = nf * O2qg;
	const Operator O2gq{g, C2gq{nf}, IntEps};
	const Operator O2gg{g, C2gg{nf}, IntEps};
	const Operator O2nsp = O2Vqq + O2Vqqb;
	const Operator O2nsm = O2Vqq - O2Vqqb;
	const Operator O2qq  = O2nsp + nf * O2ps;
	map<int,Operator> OM;
	OM.insert({EvolutionBasisQCD::PNSP, O2nsp});
	OM.insert({EvolutionBasisQCD::PNSM, O2nsm});
	OM.insert({EvolutionBasisQCD::PNSV, O2nsm});
	OM.insert({EvolutionBasisQCD::PQQ,  O2qq});
	OM.insert({EvolutionBasisQCD::PQG,  O2qgnf});
	OM.insert({EvolutionBasisQCD::PGQ,  O2gq});
	OM.insert({EvolutionBasisQCD::PGG,  O2gg});
	MatchPDFsNNLO.insert({nf,OM});
      }

    // Define map containing the TmdObjects for each nf.
    map<int,TmdObjects> TmdObj;

    // Construct sets of operators for each perturbative order for the
    // matching functions. Initialize also coefficients of: beta
    // function, GammaCusp, gammaV, and Collins-Soper anomalous
    // dimensions.
    map<int,map<int,Set<Operator>>> MatchingFunctionsPDFs;
    for (auto nf = nfi; nf <= nff; nf++)
      {
	TmdObjects obj;

	// Threshold
	obj.Threshold = Thresholds[nf-1];

	// Beta function
	obj.Beta.insert({0, beta0(nf)});
	obj.Beta.insert({1, beta1(nf)});
	obj.Beta.insert({2, beta2(nf)});

	// GammaCusp
	obj.GammaCuspq.insert({0, 4 * CF * GammaCusp0()});
	obj.GammaCuspq.insert({1, 4 * CF * GammaCusp1(nf)});
	obj.GammaCuspq.insert({2, 4 * CF * GammaCusp2(nf)});
	obj.GammaCuspg.insert({0, 4 * CA * GammaCusp0()});
	obj.GammaCuspg.insert({1, 4 * CA * GammaCusp1(nf)});
	obj.GammaCuspg.insert({2, 4 * CA * GammaCusp2(nf)});

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

	// Log(zeta) coefficients (for implementing the zeta-
	// prescription)
	obj.Lzetaq.insert({0, {Lzetaq10(), Lzetaq11()}});
	obj.Lzetaq.insert({1, {Lzetaq20(nf), 0, Lzetaq22(nf)}});
	obj.Lzetag.insert({0, {Lzetag10(nf), Lzetag11()}});
	obj.Lzetag.insert({1, {Lzetag20(nf), 0, Lzetag22(nf)}});

	// Matching functions.
	const EvolutionBasisQCD evb{nf};
	obj.MatchingFunctionsPDFs.insert({0, Set<Operator>{evb, MatchLO}});
	obj.MatchingFunctionsPDFs.insert({1, Set<Operator>{evb, MatchPDFsNLO.at(nf)}});
	obj.MatchingFunctionsPDFs.insert({2, Set<Operator>{evb, MatchPDFsNNLO.at(nf)}});

	//obj.MatchingFunctionsFFs.insert({0, Set<Operator>{evb, MatchLO}});
	//obj.MatchingFunctionsFFs.insert({1, Set<Operator>{evb, MatchFFsNLO}});
	//obj.MatchingFunctionsFFs.insert({2, Set<Operator>{evb, MatchFFsNNLO}});

	TmdObj.insert({nf, obj});
      }
    t.stop();

    return TmdObj;
  }

  //_____________________________________________________________________________
  function<Set<Distribution>(double const&, double const&, double const&)> BuildTmdPDFs(map<int,TmdObjects>                            const& TmdObj,
											map<int,DglapObjects>                          const& DglapObj,
											function<Set<Distribution>(double const&)>     const& CollPDFs,
											function<double(double const&, double const&)> const& fNP,
											function<double(double const&)>                const& Mu0b,
											function<double(double const&)>                const& Mub,
											int                                            const& PerturbativeOrder,
											function<double(double const&)>                const& Alphas,
											double                                         const& IntEps)
  {
    // Computed TMDs at the initial scale by convoluting PDFs,
    // matching functions, and non-perturbative function.
    const auto MatchedTmdPDFs = MatchTmdPDFs(TmdObj, DglapObj, CollPDFs, fNP, Mub, PerturbativeOrder, Alphas);

    // Compute TMD evolution factors.
    const auto EvolFactors = EvolutionFactors(TmdObj, Mu0b, Mub, PerturbativeOrder, Alphas, IntEps);

    // Construct the function that returns the product.
    const auto EvolvedTMDs = [=] (double const& b, double const& muf, double const& zetaf) -> Set<Distribution>
      {
	return EvolFactors(b, muf, zetaf) * MatchedTmdPDFs(b);
      };

    return EvolvedTMDs;
  }

  //_____________________________________________________________________________
  function<Set<Distribution>(double const&)> MatchTmdPDFs(map<int,TmdObjects>                            const& TmdObj,
							  map<int,DglapObjects>                          const& DglapObj,
							  function<Set<Distribution>(double const&)>     const& CollPDFs,
							  function<double(double const&, double const&)> const& fNP,
							  function<double(double const&)>                const& Mub,
							  int                                            const& PerturbativeOrder,
							  function<double(double const&)>                const& Alphas)
  {
    // Retrieve thresholds from "TmdObj".
    vector<double> thrs;
    for(auto const& obj : TmdObj)
      {
	const int    nf  = obj.first;
	const double thr = obj.second.Threshold;
	if ((int) thrs.size() < nf)
	  thrs.resize(nf);
	thrs[nf-1] = thr;
      }

    // Define the LX fuction as in eq. (4.9) of arXiv:1604.07869.
    const double C0 = 2 * exp(- emc);
    const auto LX = [C0] (double const& mu, double const& b) -> double{ return 2 * log( b * mu / C0 ); };

    // Matching functions as functions of the absolute value of the
    // impact parameter b.
    function<Set<Operator>(double const&)> MatchFunc;
    if (PerturbativeOrder == 0)
      MatchFunc = [=] (double const& b) -> Set<Operator>
	{
	  return TmdObj.at(NF(Mub(b),thrs)).MatchingFunctionsPDFs.at(0);
	};
    else if (PerturbativeOrder == 1)
      MatchFunc = [=] (double const& b) -> Set<Operator>
	{
	  const double mu   = Mub(b);
	  const int nf      = NF(mu,thrs);
	  const auto mf     = TmdObj.at(nf).MatchingFunctionsPDFs;
	  const auto sf     = DglapObj.at(nf).SplittingFunctions;
	  const double Lmu  = LX(mu,b);
	  const double coup = Alphas(mu) / FourPi;
	  return mf.at(0) + coup * ( - Lmu * sf.at(0) + mf.at(1) );
	};
    else if (PerturbativeOrder == 2)
      {
	// Precompute set of operators on the O(as^2) bit that are
	// proportional to the different powers of Lmu (see eq. (2.37)
	// of https://arxiv.org/pdf/1706.01473.pdf).
	map<int,Set<Operator>> SetLcoef;
	map<int,Set<Operator>> SetL2coef;
	for (auto const& to : TmdObj)
	  {
	    const int nf    = to.first;
	    const double b0 = to.second.Beta.at(0);
	    const auto mf   = to.second.MatchingFunctionsPDFs;
	    const auto sf   = DglapObj.at(nf).SplittingFunctions;

	    // Construct matricial product of P0 * P0 and C1 * P0 (see
	    // eq. (B.15) https://arxiv.org/pdf/1706.01473.pdf).
	    const auto P0 = sf.at(0);
	    const auto P1 = sf.at(1);
	    const auto C1 = mf.at(1);
	    const auto ConvMap = P0.GetMap();

	    map<int,Operator> MapL2coef;
	    MapL2coef.insert({EvolutionBasisQCD::PNSP, ( P0.at(0) * P0.at(0)                       - b0 * P0.at(0) ) / 2});
	    MapL2coef.insert({EvolutionBasisQCD::PNSM, ( P0.at(1) * P0.at(1)                       - b0 * P0.at(1) ) / 2});
	    MapL2coef.insert({EvolutionBasisQCD::PNSV, ( P0.at(2) * P0.at(2)                       - b0 * P0.at(2) ) / 2});
	    MapL2coef.insert({EvolutionBasisQCD::PQQ,  ( P0.at(3) * P0.at(3) + P0.at(4) * P0.at(5) - b0 * P0.at(3) ) / 2});
	    MapL2coef.insert({EvolutionBasisQCD::PQG,  ( P0.at(3) * P0.at(4) + P0.at(4) * P0.at(6) - b0 * P0.at(4) ) / 2});
	    MapL2coef.insert({EvolutionBasisQCD::PGQ,  ( P0.at(5) * P0.at(3) + P0.at(6) * P0.at(5) - b0 * P0.at(5) ) / 2});
	    MapL2coef.insert({EvolutionBasisQCD::PGG,  ( P0.at(5) * P0.at(4) + P0.at(6) * P0.at(6) - b0 * P0.at(6) ) / 2});

	    map<int,Operator> MapLcoef;
	    MapLcoef.insert({EvolutionBasisQCD::PNSP, P1.at(0) + C1.at(0) * P0.at(0)                       - b0 * C1.at(0)});
	    MapLcoef.insert({EvolutionBasisQCD::PNSM, P1.at(1) + C1.at(1) * P0.at(1)                       - b0 * C1.at(1)});
	    MapLcoef.insert({EvolutionBasisQCD::PNSV, P1.at(2) + C1.at(2) * P0.at(2)                       - b0 * C1.at(2)});
	    MapLcoef.insert({EvolutionBasisQCD::PQQ,  P1.at(3) + C1.at(3) * P0.at(3) + C1.at(4) * P0.at(5) - b0 * C1.at(3)});
	    MapLcoef.insert({EvolutionBasisQCD::PQG,  P1.at(4) + C1.at(3) * P0.at(4) + C1.at(4) * P0.at(6) - b0 * C1.at(4)});
	    MapLcoef.insert({EvolutionBasisQCD::PGQ,  P1.at(5) + C1.at(5) * P0.at(3) + C1.at(6) * P0.at(5) - b0 * C1.at(5)});
	    MapLcoef.insert({EvolutionBasisQCD::PGG,  P1.at(6) + C1.at(5) * P0.at(4) + C1.at(6) * P0.at(6) - b0 * C1.at(6)});

	    SetL2coef.insert({nf,Set<Operator>{ConvMap, MapL2coef}});
	    SetLcoef.insert({nf,Set<Operator>{ConvMap, MapLcoef}});
	  }

	// Now contruct the actual matching-function function.
	MatchFunc = [=] (double const& b) -> Set<Operator>
	  {
	    const double mu   = Mub(b);
	    const int nf      = NF(mu,thrs);
	    const auto mf     = TmdObj.at(nf).MatchingFunctionsPDFs;
	    const auto sf     = DglapObj.at(nf).SplittingFunctions;
	    const double Lmu  = LX(mu,b);
	    const double coup = Alphas(mu) / FourPi;
	    const auto nlo    = mf.at(1) - Lmu * sf.at(0);
	    const auto nnlo   = mf.at(2) - Lmu * ( SetLcoef.at(nf) - Lmu * SetL2coef.at(nf) );
	    return mf.at(0) + coup * ( nlo + coup * nnlo );
	  };
      }

    // Construct function that returns the product of: matching
    // functions, PDFs, NP function.
    const auto MatchedTMDs = [=] (double const& b) -> Set<Distribution>
      {
	return MatchFunc(b) * ( [&] (double const& x) -> double{ return fNP(x, b); } * CollPDFs(Mub(b)) );
      };

    return MatchedTMDs;
  }

  //_____________________________________________________________________________
  function<vector<double>(double const&, double const&, double const&)> EvolutionFactors(map<int,TmdObjects>             const& TmdObj,
											 function<double(double const&)> const& Mu0b,
											 function<double(double const&)> const& Mub,
											 int                             const& PerturbativeOrder,
											 function<double(double const&)> const& Alphas,
											 double                          const& IntEps)
  {
    // Retrieve thresholds from "TmdObj".
    vector<double> thrs;
    for(auto const& obj : TmdObj)
      {
	const int    nf  = obj.first;
	const double thr = obj.second.Threshold;
	if ((int) thrs.size() < nf)
	  thrs.resize(nf);
	thrs[nf-1] = thr;
      }

    // Define the LX fuction as in eq. (4.9) of arXiv:1604.07869.
    const double C0 = 2 * exp(- emc);
    const auto LX = [C0] (double const& mu, double const& b) -> double{ return 2 * log( b * mu / C0 ); };

    // Create functions needed for the TMD evolution.
    function<double(double const&)> gammaVq;
    function<double(double const&)> gammaVg;
    function<double(double const&)> GammaCuspq;
    function<double(double const&)> GammaCuspg;
    function<double(double const&, double const&)> DCSq;
    function<double(double const&, double const&)> DCSg;
    function<double(double const&, double const&)> zetaq;
    function<double(double const&, double const&)> zetag;
    // LL
    if (PerturbativeOrder == 0)
      {
	gammaVq = [=] (double const&) -> double{ return 0; };
	gammaVg = [=] (double const&) -> double{ return 0; };
	GammaCuspq = [=] (double const& mu) -> double
	  {
	    const double coup = Alphas(mu) / FourPi;
	    return coup * TmdObj.at(NF(mu,thrs)).GammaCuspq.at(0);
	  };
	GammaCuspg = [=] (double const& mu) -> double
	  {
	    const double coup = Alphas(mu) / FourPi;
	    return coup * TmdObj.at(NF(mu,thrs)).GammaCuspg.at(0);
	  };
	//DCSq = [=] (double const&, double const&) -> double{ return 0; };
	//DCSg = [=] (double const&, double const&) -> double{ return 0; };
	DCSq = [=] (double const& mu, double const& b) -> double
	  {
	    const auto d      = TmdObj.at(NF(mu,thrs)).CSdq;
	    const double coup = Alphas(mu) / FourPi;
	    const double Lmu  = LX(mu,b);
	    const double lo   = d.at(0)[0] + Lmu * d.at(0)[1];
	    return coup * lo;
	  };
	DCSg = [=] (double const& mu, double const& b) -> double
	  {
	    const auto d      = TmdObj.at(NF(mu,thrs)).CSdg;
	    const double coup = Alphas(mu) / FourPi;
	    const double Lmu  = LX(mu,b);
	    const double lo   = d.at(0)[0] + Lmu * d.at(0)[1];
	    return coup * lo;
	  };
	zetaq = [=] (double const& mu, double const&) -> double{ return mu * mu; };
	zetag = [=] (double const& mu, double const&) -> double{ return mu * mu; };
      }
    // NLL
    else if (PerturbativeOrder == 1)
      {
	gammaVq = [=] (double const& mu) -> double
	  {
	    const double coup = Alphas(mu) / FourPi;
	    return coup * TmdObj.at(NF(mu,thrs)).GammaVq.at(0);
	  };
	gammaVg = [=] (double const& mu) -> double
	  {
	    const double coup = Alphas(mu) / FourPi;
	    return coup * TmdObj.at(NF(mu,thrs)).GammaVg.at(0);
	  };
	GammaCuspq = [=] (double const& mu) -> double
	  {
	    const auto gc     = TmdObj.at(NF(mu,thrs)).GammaCuspq;
	    const double coup = Alphas(mu) / FourPi;
	    return coup * ( gc.at(0) + coup * gc.at(1) );
	  };
	GammaCuspg = [=] (double const& mu) -> double
	  {
	    const auto gc     = TmdObj.at(NF(mu,thrs)).GammaCuspg;
	    const double coup = Alphas(mu) / FourPi;
	    return coup * ( gc.at(0) + coup * gc.at(1) );
	  };
	//DCSq = [=] (double const& mu, double const& b) -> double
	//  {
	//    const auto d      = TmdObj.at(NF(mu,thrs)).CSdq;
	//    const double coup = Alphas(mu) / FourPi;
	//    const double Lmu  = LX(mu,b);
	//   const double lo   = d.at(0)[0] + Lmu * d.at(0)[1];
	//    return coup * lo;
	//  };
	//DCSg = [=] (double const& mu, double const& b) -> double
	//  {
	//    const auto d      = TmdObj.at(NF(mu,thrs)).CSdg;
	//    const double coup = Alphas(mu) / FourPi;
	//    const double Lmu  = LX(mu,b);
	//    const double lo   = d.at(0)[0] + Lmu * d.at(0)[1];
	//    return coup * lo;
	//  };
	DCSq = [=] (double const& mu, double const& b) -> double
	  {
	    const auto d      = TmdObj.at(NF(mu,thrs)).CSdq;
	    const double coup = Alphas(mu) / FourPi;
	    const double Lmu  = LX(mu,b);
	    const double lo   = d.at(0)[0] + Lmu * d.at(0)[1];
	    const double nlo  = d.at(1)[0] + Lmu * ( d.at(1)[1] + Lmu * d.at(1)[2] );
	    return coup * ( lo + coup * nlo );
	  };
	DCSg = [=] (double const& mu, double const& b) -> double
	  {
	    const auto d      = TmdObj.at(NF(mu,thrs)).CSdg;
	    const double coup = Alphas(mu) / FourPi;
	    const double Lmu  = LX(mu,b);
	    const double lo   = d.at(0)[0] + Lmu * d.at(0)[1];
	    const double nlo  = d.at(1)[0] + Lmu * ( d.at(1)[1] + Lmu * d.at(1)[2] );
	    return coup * ( lo + coup * nlo );
	  };
	zetaq = [=] (double const& mu, double const& b) -> double
	  {
	    const auto lz      = TmdObj.at(NF(mu,thrs)).Lzetaq;
	    const double Lmu   = LX(mu,b);
	    const double lo    = lz.at(0)[0] + Lmu * lz.at(0)[1];
	    const double lzeta = lo;
	    return 4 * exp( - lzeta + Lmu - 2 * emc ) / b / b;
	  };
	zetag = [=] (double const& mu, double const& b) -> double
	  {
	    const auto lz      = TmdObj.at(NF(mu,thrs)).Lzetag;
	    const double Lmu   = LX(mu,b);
	    const double lo    = lz.at(0)[0] + Lmu * lz.at(0)[1];
	    const double lzeta = lo;
	    return 4 * exp( - lzeta + Lmu - 2 * emc ) / b / b;
	  };
      }
    // NNLL
    else if (PerturbativeOrder == 2)
      {
	gammaVq = [=] (double const& mu) -> double
	  {
	    const auto gv     = TmdObj.at(NF(mu,thrs)).GammaVq;
	    const double coup = Alphas(mu) / FourPi;
	    return coup * ( gv.at(0) + coup * gv.at(1) );
	  };
	gammaVg = [=] (double const& mu) -> double
	  {
	    const auto gv     = TmdObj.at(NF(mu,thrs)).GammaVg;
	    const double coup = Alphas(mu) / FourPi;
	    return coup * ( gv.at(0) + coup * gv.at(1) );
	  };
	GammaCuspq = [=] (double const& mu) -> double
	  {
	    const auto gc     = TmdObj.at(NF(mu,thrs)).GammaCuspq;
	    const double coup = Alphas(mu) / FourPi;
	    return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * gc.at(2) ) );
	  };
	GammaCuspg = [=] (double const& mu) -> double
	  {
	    const auto gc     = TmdObj.at(NF(mu,thrs)).GammaCuspg;
	    const double coup = Alphas(mu) / FourPi;
	    return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * gc.at(2) ) );
	  };
	//DCSq = [=] (double const& mu, double const& b) -> double
	//  {
	//    const auto d      = TmdObj.at(NF(mu,thrs)).CSdq;
	//    const double coup = Alphas(mu) / FourPi;
	//    const double Lmu  = LX(mu,b);
	//    const double lo   = d.at(0)[0] + Lmu * d.at(0)[1];
	//    const double nlo  = d.at(1)[0] + Lmu * ( d.at(1)[1] + Lmu * d.at(1)[2] );
	//    return coup * ( lo + coup * nlo );
	//  };
	//DCSg = [=] (double const& mu, double const& b) -> double
	//  {
	//    const auto d      = TmdObj.at(NF(mu,thrs)).CSdg;
	//    const double coup = Alphas(mu) / FourPi;
	//    const double Lmu  = LX(mu,b);
	//    const double lo   = d.at(0)[0] + Lmu * d.at(0)[1];
	//    const double nlo  = d.at(1)[0] + Lmu * ( d.at(1)[1] + Lmu * d.at(1)[2] );
	//    return coup * ( lo + coup * nlo );
	//  };
	DCSq = [=] (double const& mu, double const& b) -> double
	  {
	    const auto d      = TmdObj.at(NF(mu,thrs)).CSdq;
	    const double coup = Alphas(mu) / FourPi;
	    const double Lmu  = LX(mu,b);
	    const double lo   = d.at(0)[0] + Lmu * d.at(0)[1];
	    const double nlo  = d.at(1)[0] + Lmu * ( d.at(1)[1] + Lmu * d.at(1)[2] );
	    const double nnlo = d.at(2)[0] + Lmu * ( d.at(2)[1] + Lmu * ( d.at(2)[2] + Lmu * d.at(2)[3] ) );
	    return coup * ( lo + coup * ( nlo + coup * nnlo ) );
	  };
	DCSg = [=] (double const& mu, double const& b) -> double
	  {
	    const auto d      = TmdObj.at(NF(mu,thrs)).CSdg;
	    const double coup = Alphas(mu) / FourPi;
	    const double Lmu  = LX(mu,b);
	    const double lo   = d.at(0)[0] + Lmu * d.at(0)[1];
	    const double nlo  = d.at(1)[0] + Lmu * ( d.at(1)[1] + Lmu * d.at(1)[2] );
	    const double nnlo = d.at(2)[0] + Lmu * ( d.at(2)[1] + Lmu * ( d.at(2)[2] + Lmu * d.at(2)[3] ) );
	    return coup * ( lo + coup * ( nlo + coup * nnlo ) );
	  };
	zetaq = [=] (double const& mu, double const& b) -> double
	  {
	    const auto lz      = TmdObj.at(NF(mu,thrs)).Lzetaq;
	    const double coup  = Alphas(mu) / FourPi;
	    const double Lmu   = LX(mu,b);
	    const double lo    = lz.at(0)[0] + Lmu * lz.at(0)[1];
	    const double nlo   = lz.at(1)[0] + Lmu * Lmu * lz.at(1)[2];
	    const double lzeta = lo + coup * nlo;
	    return 4 * exp( - lzeta + Lmu - 2 * emc ) / b / b;
	  };
	zetag = [=] (double const& mu, double const& b) -> double
	  {
	    const auto lz      = TmdObj.at(NF(mu,thrs)).Lzetag;
	    const double coup  = Alphas(mu) / FourPi;
	    const double Lmu   = LX(mu,b);
	    const double lo    = lz.at(0)[0] + Lmu * lz.at(0)[1];
	    const double nlo   = lz.at(1)[0] + Lmu * Lmu * lz.at(1)[2];
	    const double lzeta = lo + coup * nlo;
	    return 4 * exp( - lzeta + Lmu - 2 * emc ) / b / b;
	  };
      }

    // Define the integrands.
    Integrator I1q{[=] (double const& mu) -> double{ return ( - gammaVq(mu) + 2 * log(mu) * GammaCuspq(mu) ) / mu; }};
    Integrator I1g{[=] (double const& mu) -> double{ return ( - gammaVg(mu) + 2 * log(mu) * GammaCuspg(mu) ) / mu; }};
    Integrator I2q{[=] (double const& mu) -> double{ return GammaCuspq(mu) / mu; }};
    Integrator I2g{[=] (double const& mu) -> double{ return GammaCuspg(mu) / mu; }};

    // Construct function that returns the product of: matching
    // functions, PDFs, NP function, and evolution factors, i.e. the
    // set of evolved TMDs (times x).
    const auto EvolFactors = [=] (double const& b, double const& muf, double const& zetaf) -> vector<double>
      {
	// Define relevant scales
	const double mu0    = Mu0b(b);
	const double mui    = Mub(b);
	const double zetaiq = zetaq(mui,b);
	const double zetaig = zetag(mui,b);

	// Compute argument of the exponent of the evolution factors.
	const double LRq = I1q.integrate(mui, muf, thrs, IntEps)
	- I2q.integrate(mu0, muf, thrs, IntEps) * log(zetaf)
	+ I2q.integrate(mu0, mui, thrs, IntEps) * log(zetaiq);
	const double LRg = I1g.integrate(mui, muf, thrs, IntEps)
	- I2g.integrate(mu0, muf, thrs, IntEps) * log(zetaf)
	+ I2g.integrate(mu0, mui, thrs, IntEps) * log(zetaig);

	// Compute the evolution factors.
	const double Rq = exp( LRq - DCSq(mu0,b) * log( zetaf / zetaiq ) );
	const double Rg = exp( LRg - DCSg(mu0,b) * log( zetaf / zetaig ) );

	// Return vector of evolution factors.
	return vector<double>{Rg, Rq, Rq, Rq, Rq, Rq, Rq, Rq, Rq, Rq, Rq, Rq, Rq};
      };

    return EvolFactors;
  }
}
