//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/gtmdbuilder.h"
#include "apfel/tmdbuilder.h"
#include "apfel/gpdbuilder.h"
#include "apfel/timer.h"
#include "apfel/gmatchingfunctions.h"
#include "apfel/evolutionbasisqcd.h"
#include "apfel/betaqcd.h"
#include "apfel/gammak.h"
#include "apfel/gammaf.h"
#include "apfel/kcs.h"
#include "apfel/tools.h"
#include "apfel/integrator.h"

namespace apfel
{
  //_____________________________________________________________________________
  std::map<int, GtmdObjects> InitializeGtmdObjectsEven(Grid                const& g,
                                                       std::vector<double> const& Thresholds,
                                                       double              const& xi,
                                                       double              const& IntEps,
                                                       Polarization        const& pol)
  {
    // Initialise GPD splitting functions on the grid required to
    // compute the log terms of the matching functions.
    const std::map<int, DglapObjects> DglapObj = pol == Polarization::U ? InitializeGpdObjects(g, Thresholds, xi, false, IntEps) :
                                                 pol == Polarization::L ? InitializeGpdObjectsPol(g, Thresholds, xi, false, IntEps) :
                                                                          InitializeGpdObjectsTrans(g, Thresholds, xi, false, IntEps);

    report("Initializing parity-even unpolarised GTMD objects for matching and evolution... ");
    Timer t;

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the threshold
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
        nfi++;

    // ===============================================================
    // LO matching functions operators
    std::map<int, std::map<int, Operator>> C00;
    const Operator Id  {g, Identity{}, IntEps, true};
    const Operator Zero{g, Null{},     IntEps, true};
    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, Id});
        OM.insert({EvolutionBasisQCD::PNSM, Id});
        OM.insert({EvolutionBasisQCD::PNSV, Id});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * Id});
        OM.insert({EvolutionBasisQCD::PQG,  Zero});
        OM.insert({EvolutionBasisQCD::PGQ,  Zero});
        OM.insert({EvolutionBasisQCD::PGG,  Id});
        C00.insert({nf, OM});
      }

    // ===============================================================
    // NLO matching functions operators
    std::map<int, std::map<int, Operator>> C10;
    const Expression tmpns = pol == Polarization::U ? static_cast<Expression>(Cgtmd1nseUU{xi}): 
                             pol == Polarization::L ? static_cast<Expression>(Cgtmd1nseLL{xi}): 
                                                      static_cast<Expression>(Null{});
    const Operator O1ns{g, tmpns, IntEps, true};

    const Expression tmpqq = pol == Polarization::U ? static_cast<Expression>(Cgtmd1qqeUU{xi}) : 
                             pol == Polarization::L ? static_cast<Expression>(Cgtmd1qqeLL{xi}) : 
                                                      static_cast<Expression>(Null{});
    const Operator O1qq{g, tmpqq, IntEps, true};
    
    const Expression tmpqg = pol == Polarization::U ? static_cast<Expression>(Cgtmd1qgeUU{xi}) : 
                             pol == Polarization::L ? static_cast<Expression>(Cgtmd1qgeLL{xi}) : 
                                                      static_cast<Expression>(Null{});
    const Operator O1qg{g, tmpqg, IntEps, true};

    const Expression tmpgq = pol == Polarization::U ? static_cast<Expression>(Cgtmd1gqeUU{xi}) : 
                             pol == Polarization::L ? static_cast<Expression>(Cgtmd1gqeLL{xi}) : 
                                                      static_cast<Expression>(Null{});
    const Operator O1gq{g, tmpgq, IntEps, true};

    const Expression tmpgg = pol == Polarization::U ? static_cast<Expression>(Cgtmd1ggeUU{xi}) : 
                             pol == Polarization::L ? static_cast<Expression>(Cgtmd1ggeLL{xi}) : 
                                                      static_cast<Expression>(Cgtmd1ggeTT{xi});
    const Operator O1gg{g, tmpgg, IntEps, true};

    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O1ns});
        OM.insert({EvolutionBasisQCD::PNSM, O1ns});
        OM.insert({EvolutionBasisQCD::PNSV, O1ns});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O1qq});
        OM.insert({EvolutionBasisQCD::PQG,           nf * O1qg});
        OM.insert({EvolutionBasisQCD::PGQ,  ( nf / 6. ) * O1gq});
        OM.insert({EvolutionBasisQCD::PGG,                O1gg});
        C10.insert({nf, OM});
      }

    // Terms proportion to one power of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C11;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O11gmVq = gammaFq0() * Id;
        const Operator O11gmVg = gammaFg0(nf) * Id;
        const auto P0 = DglapObj.at(nf).SplittingFunctions.at(0);
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O11gmVq - 2 * P0.at(0)});
        OM.insert({EvolutionBasisQCD::PNSM, O11gmVq - 2 * P0.at(1)});
        OM.insert({EvolutionBasisQCD::PNSV, O11gmVq - 2 * P0.at(2)});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * ( O11gmVq - 2 * P0.at(3) )});
        OM.insert({EvolutionBasisQCD::PQG,                          - 2 * P0.at(4)});
        OM.insert({EvolutionBasisQCD::PGQ,                - ( nf / 3. ) * P0.at(5)});
        OM.insert({EvolutionBasisQCD::PGG,                  O11gmVg - 2 * P0.at(6)});
        C11.insert({nf, OM});
      }

    // Terms proportion to two powers of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C12;
    const Operator O12gmKq = - CF * gammaK0() / 2 * Id;
    const Operator O12gmKg = - CA * gammaK0() / 2 * Id;
    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP,               O12gmKq});
        OM.insert({EvolutionBasisQCD::PNSM,               O12gmKq});
        OM.insert({EvolutionBasisQCD::PNSV,               O12gmKq});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O12gmKq});
        OM.insert({EvolutionBasisQCD::PQG,                   Zero});
        OM.insert({EvolutionBasisQCD::PGQ,                   Zero});
        OM.insert({EvolutionBasisQCD::PGG,                O12gmKg});
        C12.insert({nf, OM});
      }

    // Define map containing the GtmdObjects for each nf
    std::map<int, GtmdObjects> GtmdObj;

    // Construct a set of operators for each perturbative order for
    // the matching functions. Initialize also coefficients of: beta
    // function, gammaK, gammaF, and Collins-Soper anomalous
    // dimensions.
    for (int nf = nfi; nf <= nff; nf++)
      {
        GtmdObjects obj;

        // Threshold
        obj.Threshold = Thresholds[nf-1];

        // Skewness
        obj.xi = xi;

        // Beta function
        obj.Beta.insert({0, beta0qcd(nf)});
        obj.Beta.insert({1, beta1qcd(nf)});

        // GammaF quark
        obj.GammaFq.insert({0, gammaFq0()});
        obj.GammaFq.insert({1, gammaFq1(nf)});

        // GammaF gluon
        obj.GammaFg.insert({0, gammaFg0(nf)});
        obj.GammaFg.insert({1, gammaFg1(nf)});

        // gammaK (multiply by CF for quarks and by CA for gluons)
        obj.GammaK.insert({0, gammaK0()});
        obj.GammaK.insert({1, gammaK1(nf)});
        obj.GammaK.insert({2, gammaK2(nf)});

        // Collins-Soper anomalous dimensions (multiply by CF for
        // quarks and by CA for gluons).
        obj.KCS.insert({0, {KCS00(),   KCS01()}});
        obj.KCS.insert({1, {KCS10(nf), KCS11(nf), KCS12(nf)}});

        // Matching functions
        const EvolutionBasisQCD evb{nf};
        obj.MatchingFunctions.insert({0, {{evb, C00.at(nf)}}});
        obj.MatchingFunctions.insert({1, {{evb, C10.at(nf)}, {evb, C11.at(nf)}, {evb, C12.at(nf)}}});

        // Insert full object
        GtmdObj.insert({nf, obj});
      }
    t.stop();

    return GtmdObj;
  }


  //_____________________________________________________________________________
  std::map<int, GtmdObjects> InitializeGtmdObjectsOdd(Grid                const& g,
                                                      std::vector<double> const& Thresholds,
                                                      double              const& xi,
                                                      double              const& IntEps,
                                                      Polarization        const& pol)
  {
    report("Initializing parity-odd unpolarised GTMD objects for matching and evolution... ");
    Timer t;

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the threshold
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
        nfi++;

    // Construct map of null operators to be used below
    const Operator Zero{g, Null{}, IntEps, true};
    std::map<int, Operator> ZeroSet;
    ZeroSet.insert({EvolutionBasisQCD::PNSP, Zero});
    ZeroSet.insert({EvolutionBasisQCD::PNSM, Zero});
    ZeroSet.insert({EvolutionBasisQCD::PNSV, Zero});
    ZeroSet.insert({EvolutionBasisQCD::PQQ,  Zero});
    ZeroSet.insert({EvolutionBasisQCD::PQG,  Zero});
    ZeroSet.insert({EvolutionBasisQCD::PGQ,  Zero});
    ZeroSet.insert({EvolutionBasisQCD::PGG,  Zero});

    // ===============================================================
    // NLO matching functions operators
    std::map<int, std::map<int, Operator>> C10;
    
    const Expression tmpqg = pol == Polarization::U ? static_cast<Expression>(Cgtmd1qgoUU{xi}) : 
                             pol == Polarization::L ? static_cast<Expression>(Cgtmd1qgoLL{xi}) : 
                                                      static_cast<Expression>(Null{});
    const Operator O1qg{g, tmpqg, IntEps, true};

    const Expression tmpgg = pol == Polarization::U ? static_cast<Expression>(Cgtmd1ggoUU{xi}) : 
                             pol == Polarization::L ? static_cast<Expression>(Cgtmd1ggoLL{xi}) : 
                                                      static_cast<Expression>(Cgtmd1ggoTT{xi});
    const Operator O1gg{g, tmpgg, IntEps, true};

    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP,      Zero});
        OM.insert({EvolutionBasisQCD::PNSM,      Zero});
        OM.insert({EvolutionBasisQCD::PNSV,      Zero});
        OM.insert({EvolutionBasisQCD::PQQ,       Zero});
        OM.insert({EvolutionBasisQCD::PQG,  nf * O1qg});
        OM.insert({EvolutionBasisQCD::PGQ,       Zero});
        OM.insert({EvolutionBasisQCD::PGG,       O1gg});
        C10.insert({nf, OM});
      }

    // Define map containing the GtmdObjects for each nf
    std::map<int, GtmdObjects> GtmdObj;

    // Construct a set of operators for each perturbative order for
    // the matching functions. Initialize also coefficients of: beta
    // function, gammaK, gammaF, and Collins-Soper anomalous
    // dimensions.
    for (int nf = nfi; nf <= nff; nf++)
      {
        GtmdObjects obj;

        // Threshold
        obj.Threshold = Thresholds[nf-1];

        // Skewness
        obj.xi = xi;

        // Beta function
        obj.Beta.insert({0, beta0qcd(nf)});
        obj.Beta.insert({1, beta1qcd(nf)});

        // GammaF quark
        obj.GammaFq.insert({0, gammaFq0()});
        obj.GammaFq.insert({1, gammaFq1(nf)});

        // GammaF gluon
        obj.GammaFg.insert({0, gammaFg0(nf)});
        obj.GammaFg.insert({1, gammaFg1(nf)});

        // gammaK (multiply by CF for quarks and by CA for gluons)
        obj.GammaK.insert({0, gammaK0()});
        obj.GammaK.insert({1, gammaK1(nf)});
        obj.GammaK.insert({2, gammaK2(nf)});

        // Collins-Soper anomalous dimensions (multiply by CF for
        // quarks and by CA for gluons).
        obj.KCS.insert({0, {KCS00(),   KCS01()}});
        obj.KCS.insert({1, {KCS10(nf), KCS11(nf), KCS12(nf)}});

        // Matching functions
        const EvolutionBasisQCD evb{nf};
        obj.MatchingFunctions.insert({0, {{evb, ZeroSet}}});
        obj.MatchingFunctions.insert({1, {{evb, C10.at(nf)}, {evb, ZeroSet}, {evb, ZeroSet}}});

        // Insert full object
        GtmdObj.insert({nf, obj});
      }
    t.stop();

    return GtmdObj;
  }

  //_____________________________________________________________________________
  std::map<int, GtmdObjects> InitializeGtmdObjectsYT(Grid                const& g,
                                                     std::vector<double> const& Thresholds,
                                                     double              const& xi,
                                                     double              const& IntEps,
                                                     Polarization        const& pol,
                                                     bool                const& even)
  {
    report("Initializing parity-odd unpolarised GTMD objects for matching and evolution... ");
    if(pol == Polarization::T) throw std::runtime_error(error("InitializeGtmdObjectsEvenYT", "Cannot be invoked for diagonal TT channel."));
    Timer t;

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the threshold
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
        nfi++;

    // Construct map of null operators to be used below
    const Operator Zero{g, Null{}, IntEps, true};
    std::map<int, Operator> ZeroSet;
    ZeroSet.insert({EvolutionBasisQCD::PNSP, Zero});
    ZeroSet.insert({EvolutionBasisQCD::PNSM, Zero});
    ZeroSet.insert({EvolutionBasisQCD::PNSV, Zero});
    ZeroSet.insert({EvolutionBasisQCD::PQQ,  Zero});
    ZeroSet.insert({EvolutionBasisQCD::PQG,  Zero});
    ZeroSet.insert({EvolutionBasisQCD::PGQ,  Zero});
    ZeroSet.insert({EvolutionBasisQCD::PGG,  Zero});

    // ===============================================================
    // NLO matching functions operators
    // Only q->g and g->g are non-vanishing
    std::map<int, std::map<int, Operator>> C10;
    
    const Expression tmpqg = pol == Polarization::U ? (even ? static_cast<Expression>(Cgtmd1qgeUT{xi}) : static_cast<Expression>(Cgtmd1qgoUT{xi})) : 
                                                      (even ? static_cast<Expression>(Cgtmd1qgeLT{xi}) : static_cast<Expression>(Cgtmd1qgoLT{xi}));
    const Operator O1qg{g, tmpqg, IntEps, true};

    const Expression tmpgg = pol == Polarization::U ? (even ? static_cast<Expression>(Cgtmd1ggeUT{xi}) : static_cast<Expression>(Cgtmd1ggoUT{xi})) : 
                                                      (even ? static_cast<Expression>(Cgtmd1ggeLT{xi}) : static_cast<Expression>(Cgtmd1ggoLT{xi}));
    const Operator O1gg{g, tmpgg, IntEps, true};

    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP,      Zero});
        OM.insert({EvolutionBasisQCD::PNSM,      Zero});
        OM.insert({EvolutionBasisQCD::PNSV,      Zero});
        OM.insert({EvolutionBasisQCD::PQQ,       Zero});
        OM.insert({EvolutionBasisQCD::PQG,  nf * O1qg});
        OM.insert({EvolutionBasisQCD::PGQ,       Zero});
        OM.insert({EvolutionBasisQCD::PGG,       O1gg});
        C10.insert({nf, OM});
      }

    // Define map containing the GtmdObjects for each nf
    std::map<int, GtmdObjects> GtmdObj;

    // Construct a set of operators for each perturbative order for
    // the matching functions. Initialize also coefficients of: beta
    // function, gammaK, gammaF, and Collins-Soper anomalous
    // dimensions.
    for (int nf = nfi; nf <= nff; nf++)
      {
        GtmdObjects obj;

        // Threshold
        obj.Threshold = Thresholds[nf-1];

        // Skewness
        obj.xi = xi;

        // Beta function
        obj.Beta.insert({0, beta0qcd(nf)});
        obj.Beta.insert({1, beta1qcd(nf)});

        // GammaF quark
        obj.GammaFq.insert({0, gammaFq0()});
        obj.GammaFq.insert({1, gammaFq1(nf)});

        // GammaF gluon
        obj.GammaFg.insert({0, gammaFg0(nf)});
        obj.GammaFg.insert({1, gammaFg1(nf)});

        // gammaK (multiply by CF for quarks and by CA for gluons)
        obj.GammaK.insert({0, gammaK0()});
        obj.GammaK.insert({1, gammaK1(nf)});
        obj.GammaK.insert({2, gammaK2(nf)});

        // Collins-Soper anomalous dimensions (multiply by CF for
        // quarks and by CA for gluons).
        obj.KCS.insert({0, {KCS00(),   KCS01()}});
        obj.KCS.insert({1, {KCS10(nf), KCS11(nf), KCS12(nf)}});

        // Matching functions
        const EvolutionBasisQCD evb{nf};
        obj.MatchingFunctions.insert({0, {{evb, ZeroSet}}});
        obj.MatchingFunctions.insert({1, {{evb, C10.at(nf)}, {evb, ZeroSet}, {evb, ZeroSet}}});

        // Insert full object
        GtmdObj.insert({nf, obj});
      }
    t.stop();

    return GtmdObj;
  }

  //_____________________________________________________________________________
  std::map<int, GtmdObjects> InitializeGtmdObjectsTY(Grid                const& g,
                                                     std::vector<double> const& Thresholds,
                                                     double              const& xi,
                                                     double              const& IntEps,
                                                     Polarization        const& pol,
                                                     bool                const& even)
  {
    report("Initializing parity-odd unpolarised GTMD objects for matching and evolution... ");
    if(pol == Polarization::T) throw std::runtime_error(error("InitializeGtmdObjectsEvenYT", "Cannot be invoked for diagonal TT channel."));
    Timer t;

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the threshold
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
        nfi++;

    // Construct map of null operators to be used below
    const Operator Zero{g, Null{}, IntEps, true};
    std::map<int, Operator> ZeroSet;
    ZeroSet.insert({EvolutionBasisQCD::PNSP, Zero});
    ZeroSet.insert({EvolutionBasisQCD::PNSM, Zero});
    ZeroSet.insert({EvolutionBasisQCD::PNSV, Zero});
    ZeroSet.insert({EvolutionBasisQCD::PQQ,  Zero});
    ZeroSet.insert({EvolutionBasisQCD::PQG,  Zero});
    ZeroSet.insert({EvolutionBasisQCD::PGQ,  Zero});
    ZeroSet.insert({EvolutionBasisQCD::PGG,  Zero});

    // ===============================================================
    // NLO matching functions operators
    // Only g->q and g->g are non-vanishing
    std::map<int, std::map<int, Operator>> C10;
    
    const Expression tmpgq = pol == Polarization::U ? (even ? static_cast<Expression>(Cgtmd1gqeTU{xi}) : static_cast<Expression>(Null{})) : 
                                                      (even ? static_cast<Expression>(Cgtmd1gqeTL{xi}) : static_cast<Expression>(Null{}));
    const Operator O1gq{g, tmpgq, IntEps, true};

    const Expression tmpgg = pol == Polarization::U ? (even ? static_cast<Expression>(Cgtmd1ggeTU{xi}) : static_cast<Expression>(Cgtmd1ggoTU{xi})) : 
                                                      (even ? static_cast<Expression>(Cgtmd1ggeTL{xi}) : static_cast<Expression>(Cgtmd1ggoTL{xi}));
    const Operator O1gg{g, tmpgg, IntEps, true};

    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP,               Zero});
        OM.insert({EvolutionBasisQCD::PNSM,               Zero});
        OM.insert({EvolutionBasisQCD::PNSV,               Zero});
        OM.insert({EvolutionBasisQCD::PQQ,                Zero});
        OM.insert({EvolutionBasisQCD::PQG,                Zero});
        OM.insert({EvolutionBasisQCD::PGQ,  ( nf / 6. ) * O1gq});
        OM.insert({EvolutionBasisQCD::PGG,                O1gg});
        C10.insert({nf, OM});
      }

    // Define map containing the GtmdObjects for each nf
    std::map<int, GtmdObjects> GtmdObj;

    // Construct a set of operators for each perturbative order for
    // the matching functions. Initialize also coefficients of: beta
    // function, gammaK, gammaF, and Collins-Soper anomalous
    // dimensions.
    for (int nf = nfi; nf <= nff; nf++)
      {
        GtmdObjects obj;

        // Threshold
        obj.Threshold = Thresholds[nf-1];

        // Skewness
        obj.xi = xi;

        // Beta function
        obj.Beta.insert({0, beta0qcd(nf)});
        obj.Beta.insert({1, beta1qcd(nf)});

        // GammaF quark
        obj.GammaFq.insert({0, gammaFq0()});
        obj.GammaFq.insert({1, gammaFq1(nf)});

        // GammaF gluon
        obj.GammaFg.insert({0, gammaFg0(nf)});
        obj.GammaFg.insert({1, gammaFg1(nf)});

        // gammaK (multiply by CF for quarks and by CA for gluons)
        obj.GammaK.insert({0, gammaK0()});
        obj.GammaK.insert({1, gammaK1(nf)});
        obj.GammaK.insert({2, gammaK2(nf)});

        // Collins-Soper anomalous dimensions (multiply by CF for
        // quarks and by CA for gluons).
        obj.KCS.insert({0, {KCS00(),   KCS01()}});
        obj.KCS.insert({1, {KCS10(nf), KCS11(nf), KCS12(nf)}});

        // Matching functions
        const EvolutionBasisQCD evb{nf};
        obj.MatchingFunctions.insert({0, {{evb, ZeroSet}}});
        obj.MatchingFunctions.insert({1, {{evb, C10.at(nf)}, {evb, ZeroSet}, {evb, ZeroSet}}});

        // Insert full object
        GtmdObj.insert({nf, obj});
      }
    t.stop();

    return GtmdObj;
  }

  //_____________________________________________________________________________
  std::function<Set<Distribution>(double const&, double const&, double const&)> BuildGtmds(std::map<int, GtmdObjects>                      const& GtmdObj,
                                                                                           std::function<Set<Distribution>(double const&)> const& CollGPDs,
                                                                                           std::function<double(double const&)>            const& Alphas,
                                                                                           int                                             const& PerturbativeOrder,
                                                                                           double                                          const& Ci,
                                                                                           double                                          const& IntEps)
  {
    // Match GTMDs onto collinear GPDs
    const std::function<Set<Distribution>(double const&)> MatchedGtmds = MatchGtmds(GtmdObj, CollGPDs, Alphas, PerturbativeOrder, Ci);

    // Compute GTMD evolution factors
    const std::function<std::vector<double>(double const&, double const&, double const&, double const&)> EvolFactors = EvolutionFactors(GtmdObj, Alphas, PerturbativeOrder, Ci, IntEps);

    // Compute GTMDs at the final scale by multiplying the initial
    // scale GTMDs by the evolution factor.
    const auto EvolvedGTMDs = [=] (double const& b, double const& muf, double const& zetaf) -> Set<Distribution>
    {
      return [=] (double const& x) -> std::vector<double> { return EvolFactors(x, b, muf, zetaf); } * MatchedGtmds(b);
    };

    return EvolvedGTMDs;
  }

  //_____________________________________________________________________________
  std::function<Set<Distribution>(double const&)> MatchGtmds(std::map<int, GtmdObjects>                      const& GtmdObj,
                                                             std::function<Set<Distribution>(double const&)> const& CollGPDs,
                                                             std::function<double(double const&)>            const& Alphas,
                                                             int                                             const& PerturbativeOrder,
                                                             double                                          const& Ci)
  {
    // Get matching functions
    const std::function<Set<Operator>(double const&)> MatchFunc = MatchingFunctions(GtmdObj, Alphas, PerturbativeOrder, Ci);

    // Construct function that returns the product of matching
    // functions and collinear GPDs.
    const auto MatchedGTMDs = [=] (double const& b) -> Set<Distribution>
    {
      // Define lower scales
      const double mu0 = Ci * 2 * exp(- emc) / b;

      // Convolute matching functions with the collinear GPDs and
      // return.
      return MatchFunc(mu0) * CollGPDs(mu0);
    };

    return MatchedGTMDs;
  }

  //_____________________________________________________________________________
  std::function<Set<Operator>(double const&)> MatchingFunctions(std::map<int, GtmdObjects>           const& GtmdObj,
                                                                std::function<double(double const&)> const& Alphas,
                                                                int                                  const& PerturbativeOrder,
                                                                double                               const& Ci)
  {
    // Retrieve thresholds from "GtmdObj"
    std::vector<double> thrs;
    for(auto const& obj : GtmdObj)
      {
        const int    nf  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) thrs.size() < nf)
          thrs.resize(nf);
        thrs[nf-1] = thr;
      }

    // Define the log(Ci) to assess scale variations
    const double Lmu = log(Ci);

    // Matching functions as functions of the absolute value of the
    // impact parameter b.
    std::function<Set<Operator>(double const&)> MatchFunc;
    if (PerturbativeOrder == LL || PerturbativeOrder == NLL)
      MatchFunc = [=] (double const& mu) -> Set<Operator>
      {
        return GtmdObj.at(NF(mu, thrs)).MatchingFunctions.at(0)[0];
      };
    else if (PerturbativeOrder == NNLL || PerturbativeOrder == NLLp)
      MatchFunc = [=] (double const& mu) -> Set<Operator>
      {
        const double coup = Alphas(mu) / FourPi;
        const auto& mf = GtmdObj.at(NF(mu, thrs)).MatchingFunctions;
        const auto c0  = mf.at(0);
        const auto c1  = mf.at(1);
        const auto lo  = c0[0];
        const auto nlo = c1[0] + Lmu * ( c1[1] + Lmu * c1[2] );
        return lo + coup * nlo;
      };

    // Return matching functions
    return MatchFunc;
  }

  //_____________________________________________________________________________
  std::function<std::vector<double>(double const&, double const&, double const&, double const&)> EvolutionFactors(std::map<int, GtmdObjects>           const& GtmdObj,
                                                                                                                  std::function<double(double const&)> const& Alphas,
                                                                                                                  int                                  const& PerturbativeOrder,
                                                                                                                  double                               const& Ci,
                                                                                                                  double                               const& IntEps)
  {
    // Retrieve thresholds from "GtmdObj"
    std::vector<double> thrs;
    for(auto const& obj : GtmdObj)
      {
        const int    nf  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) thrs.size() < nf)
          thrs.resize(nf);
        thrs[nf-1] = thr;
      }

    // Get skewness from the first element of GtmdObj
    const double xi = GtmdObj.begin()->second.xi;

    // Define the log(Ci) to assess scale variations
    const double Lmu = log(Ci);

    // Create functions needed for the GTMD evolution
    std::function<double(double const&)> gammaFq;
    std::function<double(double const&)> gammaFg;
    std::function<double(double const&)> gammaK;
    std::function<double(double const&)> K;
    // LL
    if (PerturbativeOrder == LL)
      {
        gammaFq = [=] (double const&) -> double{ return 0; };
        gammaFg = [=] (double const&) -> double{ return 0; };
        gammaK  = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * GtmdObj.at(NF(mu, thrs)).GammaK.at(0);
        };
        K = [=] (double const&) -> double{ return 0; };
      }
    // NLL
    else if (PerturbativeOrder == NLL || PerturbativeOrder == NLLp)
      {
        gammaFq = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * GtmdObj.at(NF(mu, thrs)).GammaFq.at(0);
        };
        gammaFg = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * GtmdObj.at(NF(mu, thrs)).GammaFg.at(0);
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = GtmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * gc.at(1) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = GtmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          return coup * lo;
        };
      }
    // NNLL
    else if (PerturbativeOrder == NNLL || PerturbativeOrder == NNLLp)
      {
        gammaFq = [=] (double const& mu) -> double
        {
          const auto& gv    = GtmdObj.at(NF(mu, thrs)).GammaFq;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * gv.at(1) );
        };
        gammaFg = [=] (double const& mu) -> double
        {
          const auto& gv    = GtmdObj.at(NF(mu, thrs)).GammaFg;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * gv.at(1) );
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = GtmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * gc.at(2) ) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = GtmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const std::vector<double> d1 = d.at(1);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          const double nlo  = d1[0] + Lmu * ( d1[1] + Lmu * d1[2] );
          return coup * ( lo + coup * nlo );
        };
      }

    // Define the integrands
    const Integrator I1q{[=] (double const& mu) -> double{ return gammaFq(mu) / mu; }};
    const Integrator I1g{[=] (double const& mu) -> double{ return gammaFg(mu) / mu; }};
    const Integrator I2 {[=] (double const& mu) -> double{ return gammaK(mu) / mu; }};
    const Integrator I3 {[=] (double const& mu) -> double{ return gammaK(mu) * log(mu) / mu; }};

    // Construct function that returns the perturbative evolution
    // kernel.
    const auto EvolFactors = [=] (double const& x, double const& b, double const& muf, double const& zetaf) -> std::vector<double>
    {
      // Define lower scales
      const double mu0   = Ci * 2 * exp(- emc) / b;
      const double zeta0 = mu0 * mu0;
      const double omk2 = std::abs(1 - pow(xi / x, 2));

      // Compute argument of the exponent of the evolution factors
      const double IntI1q = I1q.integrate(mu0, muf, thrs, IntEps);
      const double IntI1g = I1g.integrate(mu0, muf, thrs, IntEps);
      const double IntI2  = I2.integrate(mu0, muf, thrs, IntEps) * log(omk2 * zetaf);
      const double IntI3  = I3.integrate(mu0, muf, thrs, IntEps);

      // Compute the evolution factors
      const double Klz = ( K(mu0) * log( omk2 * zetaf / zeta0 ) - IntI2 ) / 2 + IntI3;
      const double Rq  = exp( CF * Klz + IntI1q );
      const double Rg  = exp( CA * Klz + IntI1g );

      // Return vector of evolution factors
      return std::vector<double>{Rg, Rq, Rq, Rq, Rq, Rq, Rq, Rq, Rq, Rq, Rq, Rq, Rq};
    };

    return EvolFactors;
  }

  //_____________________________________________________________________________
  std::function<double(double const&, double const&, double const&, double const&)> QuarkEvolutionFactor(std::map<int, GtmdObjects>           const& GtmdObj,
                                                                                                         std::function<double(double const&)> const& Alphas,
                                                                                                         int                                  const& PerturbativeOrder,
                                                                                                         double                               const& Ci,
                                                                                                         double                               const& IntEps)
  {
    // Retrieve thresholds from "GtmdObj"
    std::vector<double> thrs;
    for(auto const& obj : GtmdObj)
      {
        const int    nf  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) thrs.size() < nf)
          thrs.resize(nf);
        thrs[nf-1] = thr;
      }

    // Get skewness from the first element of GtmdObj
    const double xi = GtmdObj.begin()->second.xi;

    // Define the log(Ci) to assess scale variations
    const double Lmu = log(Ci);

    // Create functions needed for the GTMD evolution
    std::function<double(double const&)> gammaFq;
    std::function<double(double const&)> gammaK;
    std::function<double(double const&)> K;
    // LL
    if (PerturbativeOrder == LL)
      {
        gammaFq = [=] (double const&) -> double{ return 0; };
        gammaK  = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * GtmdObj.at(NF(mu, thrs)).GammaK.at(0);
        };
        K = [=] (double const&) -> double{ return 0; };
      }
    // NLL
    else if (PerturbativeOrder == NLL || PerturbativeOrder == NLLp)
      {
        gammaFq = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * GtmdObj.at(NF(mu, thrs)).GammaFq.at(0);
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = GtmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * gc.at(1) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = GtmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          return coup * lo;
        };
      }
    // NNLL
    else if (PerturbativeOrder == NNLL || PerturbativeOrder == NNLLp)
      {
        gammaFq = [=] (double const& mu) -> double
        {
          const auto& gv    = GtmdObj.at(NF(mu, thrs)).GammaFq;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * gv.at(1) );
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = GtmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * gc.at(2) ) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = GtmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const std::vector<double> d1 = d.at(1);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          const double nlo  = d1[0] + Lmu * ( d1[1] + Lmu * d1[2] );
          return coup * ( lo + coup * nlo );
        };
      }

    // Define the integrands
    const Integrator I1{[=] (double const& mu) -> double{ return gammaFq(mu) / mu; }};
    const Integrator I2{[=] (double const& mu) -> double{ return gammaK(mu) / mu; }};
    const Integrator I3{[=] (double const& mu) -> double{ return gammaK(mu) * log(mu) / mu; }};

    // Construct function that returns the perturbative evolution
    // kernel.
    const auto EvolFactor = [=] (double const& x, double const& b, double const& muf, double const& zetaf) -> double
    {
      // Define lower scales
      const double mu0   = Ci * 2 * exp(- emc) / b;
      const double zeta0 = mu0 * mu0;
      const double omk2 = std::abs(1 - pow(xi / x, 2));

      // Compute argument of the exponent of the evolution factors
      const double IntI1 = I1.integrate(mu0, muf, thrs, IntEps);
      const double IntI2 = I2.integrate(mu0, muf, thrs, IntEps) * log(omk2 * zetaf);
      const double IntI3 = I3.integrate(mu0, muf, thrs, IntEps);

      // Compute the evolution factors
      const double Klz = ( K(mu0) * log( omk2 * zetaf / zeta0 ) - IntI2 ) / 2 + IntI3;
      const double Rq  = exp( CF * Klz + IntI1 );

      // Return the evolution factor
      return Rq;
    };

    return EvolFactor;
  }

  //_____________________________________________________________________________
  std::function<double(double const&, double const&, double const&, double const&)> GluonEvolutionFactor(std::map<int, GtmdObjects>           const& GtmdObj,
                                                                                                         std::function<double(double const&)> const& Alphas,
                                                                                                         int                                  const& PerturbativeOrder,
                                                                                                         double                               const& Ci,
                                                                                                         double                               const& IntEps)
  {
    // Retrieve thresholds from "GtmdObj"
    std::vector<double> thrs;
    for(auto const& obj : GtmdObj)
      {
        const int    nf  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) thrs.size() < nf)
          thrs.resize(nf);
        thrs[nf-1] = thr;
      }

    // Get skewness from the first element of GtmdObj
    const double xi = GtmdObj.begin()->second.xi;

    // Define the log(Ci) to assess scale variations
    const double Lmu = log(Ci);

    // Create functions needed for the GTMD evolution
    std::function<double(double const&)> gammaFg;
    std::function<double(double const&)> gammaK;
    std::function<double(double const&)> K;
    // LL
    if (PerturbativeOrder == LL)
      {
        gammaFg = [=] (double const&) -> double{ return 0; };
        gammaK  = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * GtmdObj.at(NF(mu, thrs)).GammaK.at(0);
        };
        K = [=] (double const&) -> double{ return 0; };
      }
    // NLL
    else if (PerturbativeOrder == NLL || PerturbativeOrder == NLLp)
      {
        gammaFg = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * GtmdObj.at(NF(mu, thrs)).GammaFg.at(0);
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = GtmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * gc.at(1) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = GtmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          return coup * lo;
        };
      }
    // NNLL
    else if (PerturbativeOrder == NNLL || PerturbativeOrder == NNLLp)
      {
        gammaFg = [=] (double const& mu) -> double
        {
          const auto& gv    = GtmdObj.at(NF(mu, thrs)).GammaFg;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * gv.at(1) );
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = GtmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * gc.at(2) ) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = GtmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const std::vector<double> d1 = d.at(1);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          const double nlo  = d1[0] + Lmu * ( d1[1] + Lmu * d1[2] );
          return coup * ( lo + coup * nlo );
        };
      }

    // Define the integrands
    const Integrator I1{[=] (double const& mu) -> double{ return gammaFg(mu) / mu; }};
    const Integrator I2{[=] (double const& mu) -> double{ return gammaK(mu) / mu; }};
    const Integrator I3{[=] (double const& mu) -> double{ return gammaK(mu) * log(mu) / mu; }};

    // Construct function that returns the perturbative evolution
    // kernel.
    const auto EvolFactor = [=] (double const& x, double const& b, double const& muf, double const& zetaf) -> double
    {
      // Define lower scales
      const double mu0   = Ci * 2 * exp(- emc) / b;
      const double zeta0 = mu0 * mu0;
      const double omk2 = std::abs(1 - pow(xi / x, 2));

      // Compute argument of the exponent of the evolution factors
      const double IntI1 = I1.integrate(mu0, muf, thrs, IntEps);
      const double IntI2 = I2.integrate(mu0, muf, thrs, IntEps) * log(omk2 * zetaf);
      const double IntI3 = I3.integrate(mu0, muf, thrs, IntEps);

      // Compute the evolution factors
      const double Klz = ( K(mu0) * log( omk2 * zetaf / zeta0 ) - IntI2 ) / 2 + IntI3;
      const double Rg  = exp( CA * Klz + IntI1 );

      // Return the factor
      return Rg;
    };

    return EvolFactor;
  }

  //_____________________________________________________________________________
  std::function<double(double const&, double const&)> CollinsSoperKernel(std::map<int, GtmdObjects>           const& GtmdObj,
                                                                         std::function<double(double const&)> const& Alphas,
                                                                         int                                  const& PerturbativeOrder,
                                                                         double                               const& Ci,
                                                                         double                               const& IntEps)
  {
    // Retrieve thresholds from "GtmdObj"
    std::vector<double> thrs;
    for(auto const& obj : GtmdObj)
      {
        const int    nf  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) thrs.size() < nf)
          thrs.resize(nf);
        thrs[nf-1] = thr;
      }

    // Define the log(Ci) to assess scale variations
    const double Lmu = log(Ci);

    // Create functions needed for the TMD evolution
    std::function<double(double const&)> gammaK;
    std::function<double(double const&)> K;
    // LL
    if (PerturbativeOrder == LL)
      {
        gammaK  = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * GtmdObj.at(NF(mu, thrs)).GammaK.at(0);
        };
        K = [=] (double const&) -> double{ return 0; };
      }
    // NLL
    else if (PerturbativeOrder == NLL || PerturbativeOrder == NLLp)
      {
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = GtmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * gc.at(1) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = GtmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          return coup * lo;
        };
      }
    // NNLL
    else if (PerturbativeOrder == NNLL || PerturbativeOrder == NNLLp)
      {
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = GtmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * gc.at(2) ) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = GtmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const std::vector<double> d1 = d.at(1);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          const double nlo  = d1[0] + Lmu * ( d1[1] + Lmu * d1[2] );
          return coup * ( lo + coup * nlo );
        };
      }
    else
      throw std::runtime_error(error("CollinsSoperKernel", "Perturbative order not available."));

    // Define the integrands
    const Integrator I2{[=] (double const& mu) -> double{ return gammaK(mu) / mu; }};

    // Construct function that returns the perturbative evolution
    // kernel.
    const auto CSKernel = [=] (double const& b, double const& muf) -> double
    {
      // Define lower scales
      const double mu0 = Ci * 2 * exp(- emc) / b;

      // Compute argument of the exponent of the evolution factors
      const double IntI2 = I2.integrate(mu0, muf, thrs, IntEps);

      // Compute the evolution factors
      const double Klz = CF * ( K(mu0) - IntI2 );

      // Return the evolution factor
      return Klz;
    };

    return CSKernel;
  }
}
