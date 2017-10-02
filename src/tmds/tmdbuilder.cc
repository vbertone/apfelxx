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
    map<int,Operator> MatchPDFsNLO;
    const Operator O1ns{g, C1ns{}, IntEps};
    const Operator O1qg{g, C1qg{}, IntEps};
    const Operator O1gq{g, C1gq{}, IntEps};
    const Operator O1gg{g, C1gg{}, IntEps};
    MatchPDFsNLO.insert({EvolutionBasisQCD::PNSP, O1ns});
    MatchPDFsNLO.insert({EvolutionBasisQCD::PNSM, O1ns});
    MatchPDFsNLO.insert({EvolutionBasisQCD::PNSV, O1ns});
    MatchPDFsNLO.insert({EvolutionBasisQCD::PQQ,  O1ns});
    MatchPDFsNLO.insert({EvolutionBasisQCD::PQG,  O1qg});
    MatchPDFsNLO.insert({EvolutionBasisQCD::PGQ,  O1gq});
    MatchPDFsNLO.insert({EvolutionBasisQCD::PGG,  O1gg});

    // Define object of the structure containing the TmdObjects
    map<int,TmdObjects> TmdObj;

    // Construcuct sets of operators for each perturbative order for
    // the matching functions. Initialize also coefficients of: beta
    // function, GammaCusp, gammaV, and Collins-Soper anomalous
    // dimensions.
    map<int,map<int,Set<Operator>>> MatchingFunctionsPDFs;
    for (auto nf = nfi; nf <= nff; nf++)
      {
	TmdObjects obj;

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
	obj.MatchingFunctionsPDFs.insert({1, Set<Operator>{evb, MatchPDFsNLO}});
	//obj.MatchingFunctionsPDFs.insert({2, Set<Operator>{evb, MatchPDFsNNLO}});

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
											TabulateObject<Set<Distribution>>              const& CollPDFs,
											function<double(double const&, double const&)> const& fNP,
											function<double(double const&)>                const& Mu0b,
											function<double(double const&)>                const& Mub,
											int                                            const& PerturbativeOrder,
											function<double(double const&)>                const& Alphas,
											double                                         const& IntEps)
  {
    // Get thresholds from the tabulated PDFs.
    const vector<double> thrs = CollPDFs.GetThresholds();

    // Define the LX fuction as in eq. (4.9) of arXiv:1604.07869.
    const double C0 = 2 * exp(- emc);
    const auto LX = [C0] (double const& mu, double const& b) -> double{ return 2 * log( b * mu / C0 ); };

    // =================================================================
    // TMD matching on collinear distributions. What follows assume the
    // use of the zeta-prescription. As consequence the matching
    // functions have a simplified form and zeta depends on mu in a
    // little convoluted way that depends on the perturbative order and
    // on whether it enters the evolution of quarks or gluon.
    // =================================================================
    // Matching functions as functions of the absolute value of the
    // impact parameter b.
    function<Set<Operator>(double const&, double const&)> MatchFunc;
    if (PerturbativeOrder == 0)
      MatchFunc = [=] (double const& mu, double const&) -> Set<Operator>
	{
	  return TmdObj.at(NF(mu,thrs)).MatchingFunctionsPDFs.at(0);
	};
    else if (PerturbativeOrder == 1)
      MatchFunc = [=] (double const& mu, double const& b) -> Set<Operator>
	{
	  const double Lmu  = LX(mu,b);
	  const double coup = Alphas(mu) / FourPi;
	  const double nf   = NF(mu,thrs);
	  return TmdObj.at(nf).MatchingFunctionsPDFs.at(0) + coup * ( - Lmu * DglapObj.at(nf).SplittingFunctions.at(0) + TmdObj.at(nf).MatchingFunctionsPDFs.at(1) );
	};
    else if (PerturbativeOrder == 2)
      MatchFunc = [=] (double const& mu, double const& b) -> Set<Operator>
	{
	  const double Lmu  = LX(mu,b);
	  const double coup = Alphas(mu) / FourPi;
	  const double nf   = NF(mu,thrs);
	  return TmdObj.at(nf).MatchingFunctionsPDFs.at(0) + coup * ( - Lmu * DglapObj.at(nf).SplittingFunctions.at(0) + TmdObj.at(nf).MatchingFunctionsPDFs.at(1) );
	};

    // =================================================================
    // TMD evolution
    // =================================================================
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
	gammaVq = [=] (double const&)->double{ return 0; };
	gammaVg = [=] (double const&)->double{ return 0; };
	GammaCuspq = [=] (double const& mu)->double
	  {
	    const double coup = Alphas(mu) / FourPi;
	    return coup * TmdObj.at(NF(mu,thrs)).GammaCuspq.at(0);
	  };
	GammaCuspg = [=] (double const& mu)->double
	  {
	    const double coup = Alphas(mu) / FourPi;
	    return coup * TmdObj.at(NF(mu,thrs)).GammaCuspg.at(0);
	  };
	DCSq = [=] (double const&, double const&)->double{ return 0; };
	DCSg = [=] (double const&, double const&)->double{ return 0; };
	zetaq = [=] (double const& mu, double const&)->double{ return mu * mu; };
	zetag = [=] (double const& mu, double const&)->double{ return mu * mu; };
      }
    // NLL
    else if (PerturbativeOrder == 1)
      {
	gammaVq = [=] (double const& mu)->double
	  {
	    const double coup = Alphas(mu) / FourPi;
	    return coup * TmdObj.at(NF(mu,thrs)).GammaVq.at(0);
	  };
	gammaVg = [=] (double const& mu)->double
	  {
	    const double coup = Alphas(mu) / FourPi;
	    return coup * TmdObj.at(NF(mu,thrs)).GammaVg.at(0);
	  };
	GammaCuspq = [=] (double const& mu)->double
	  {
	    const auto gc     = TmdObj.at(NF(mu,thrs)).GammaCuspq;
	    const double coup = Alphas(mu) / FourPi;
	    return coup * ( gc.at(0) + coup * gc.at(1) );
	  };
	GammaCuspg = [=] (double const& mu)->double
	  {
	    const auto gc     = TmdObj.at(NF(mu,thrs)).GammaCuspg;
	    const double coup = Alphas(mu) / FourPi;
	    return coup * ( gc.at(0) + coup * gc.at(1) );
	  };
	DCSq = [=] (double const& mu, double const& b)->double
	  {
	    const auto d      = TmdObj.at(NF(mu,thrs)).CSdq;
	    const double coup = Alphas(mu) / FourPi;
	    const double Lmu  = LX(mu,b);
	    const double lo   = d.at(0)[0] + Lmu * d.at(0)[1];
	    return coup * lo;
	  };
	DCSg = [=] (double const& mu, double const& b)->double
	  {
	    const auto d      = TmdObj.at(NF(mu,thrs)).CSdg;
	    const double coup = Alphas(mu) / FourPi;
	    const double Lmu  = LX(mu,b);
	    const double lo   = d.at(0)[0] + Lmu * d.at(0)[1];
	    return coup * lo;
	  };
	zetaq = [=] (double const& mu, double const& b)->double
	  {
	    const auto lz      = TmdObj.at(NF(mu,thrs)).Lzetaq;
	    const double Lmu   = LX(mu,b);
	    const double lo    = lz.at(0)[0] + Lmu * lz.at(0)[1];
	    const double lzeta = lo;
	    return 4 * exp( - lzeta + Lmu - 2 * emc ) / b / b;
	  };
	zetag = [=] (double const& mu, double const& b)->double
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
	gammaVq = [=] (double const& mu)->double
	  {
	    const auto gv     = TmdObj.at(NF(mu,thrs)).GammaVq;
	    const double coup = Alphas(mu) / FourPi;
	    return coup * ( gv.at(0) + coup * gv.at(1) );
	  };
	gammaVg = [=] (double const& mu)->double
	  {
	    const auto gv     = TmdObj.at(NF(mu,thrs)).GammaVg;
	    const double coup = Alphas(mu) / FourPi;
	    return coup * ( gv.at(0) + coup * gv.at(1) );
	  };
	GammaCuspq = [=] (double const& mu)->double
	  {
	    const auto gc     = TmdObj.at(NF(mu,thrs)).GammaCuspq;
	    const double coup = Alphas(mu) / FourPi;
	    return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * gc.at(2) ) );
	  };
	GammaCuspg = [=] (double const& mu)->double
	  {
	    const auto gc     = TmdObj.at(NF(mu,thrs)).GammaCuspg;
	    const double coup = Alphas(mu) / FourPi;
	    return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * gc.at(2) ) );
	  };
	DCSq = [=] (double const& mu, double const& b)->double
	  {
	    const auto d      = TmdObj.at(NF(mu,thrs)).CSdq;
	    const double coup = Alphas(mu) / FourPi;
	    const double Lmu  = LX(mu,b);
	    const double lo   = d.at(0)[0] + Lmu * d.at(0)[1];
	    const double nlo  = d.at(1)[0] + Lmu * ( d.at(1)[1] + Lmu * d.at(1)[2] );
	    return coup * ( lo + coup * nlo );
	  };
	DCSg = [=] (double const& mu, double const& b)->double
	  {
	    const auto d      = TmdObj.at(NF(mu,thrs)).CSdg;
	    const double coup = Alphas(mu) / FourPi;
	    const double Lmu  = LX(mu,b);
	    const double lo   = d.at(0)[0] + Lmu * d.at(0)[1];
	    const double nlo  = d.at(1)[0] + Lmu * ( d.at(1)[1] + Lmu * d.at(1)[2] );
	    return coup * ( lo + coup * nlo );
	  };
	zetaq = [=] (double const& mu, double const& b)->double
	  {
	    const auto lz      = TmdObj.at(NF(mu,thrs)).Lzetaq;
	    const double coup  = Alphas(mu) / FourPi;
	    const double Lmu   = LX(mu,b);
	    const double lo    = lz.at(0)[0] + Lmu * lz.at(0)[1];
	    const double nlo   = lz.at(1)[0] + Lmu * Lmu * lz.at(1)[2];
	    const double lzeta = lo + coup * nlo;
	    return 4 * exp( - lzeta + Lmu - 2 * emc ) / b / b;
	  };
	zetag = [=] (double const& mu, double const& b)->double
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
    const auto EvolvedTMDs = [=] (double const& b, double const& muf, double const& zetaf)->Set<Distribution>
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

	// Compute the actual evolution factors
	const double Rq = exp( LRq - DCSq(mu0,b) * log( zetaf / zetaiq ) );
	const double Rg = exp( LRg - DCSg(mu0,b) * log( zetaf / zetaig ) );

	// Multiply PDFs by the non-perturbative function and the
	// evolution factors.
	const auto CollPDFsAtMu = CollPDFs.Evaluate(mui);
	const map<int,Distribution> PertPDFs = CollPDFsAtMu.GetObjects();
	map<int,Distribution> DistMap;
	for (auto const& id : PertPDFs)
	  DistMap.insert({id.first, (id.first == 0 ? Rg : Rq) * id.second * [=] (double const& x)->double{ return fNP(x, b); }});

	Set<Distribution> EvNonPertPDFs{CollPDFsAtMu.GetMap(), DistMap};

	return MatchFunc(mui,b) * EvNonPertPDFs;
      };

    return EvolvedTMDs;
  }
}
