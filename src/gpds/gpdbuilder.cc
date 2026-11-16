//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/gpdbuilder.h"
#include "apfel/dglap.h"
#include "apfel/set.h"
#include "apfel/timer.h"
#include "apfel/constants.h"
#include "apfel/tools.h"
#include "apfel/messages.h"
#include "apfel/gpdsplittingfunctionsunp_sl.h"
#include "apfel/gpdsplittingfunctionspol_sl.h"
#include "apfel/gpdsplittingfunctionstrans_sl.h"
#include "apfel/gpdmatchingconditionsunp_sl.h"
#include "apfel/matchingconditions_sl.h"
#include "apfel/evolutionbasisqcd.h"
#include "apfel/matchingbasisqcd.h"

namespace apfel
{
  //_____________________________________________________________________________
  std::map<int, DglapObjects> InitializeGpdObjects(Grid                const& g,
                                                   std::vector<double> const& Masses,
                                                   std::vector<double> const& Thresholds,
                                                   double              const& xi,
                                                   bool                const& OpEvol,
                                                   double              const& IntEps)
  {
    report("Initializing DglapObjects for GPD unpolarised evolution at xi = " + std::to_string(xi) + "... ");
    Timer t;

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the threshold
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
        nfi++;

    // Compute logs of muth2 / m2 needed for the matching
    // conditions. Push a zero at last to extend the vector but that
    // entry will never be effectively used in the evolution.
    std::vector<double> LogKth;
    for (int im = 0; im < (int) Thresholds.size(); im++)
      if (Thresholds[im] < eps12 || Masses[im] < eps12)
        LogKth.push_back(0);
      else
        LogKth.push_back(2 * log( Thresholds[im] / Masses[im] ));
    LogKth.push_back(0);

    // If xi happens to coincide with one of the nodes of the x-space
    // grid, numerical instabilities may prevent the integrals from
    // converging. In this eventuality, displace xi by a small amount.
    double xid = xi;
    for (double x : g.GetJointGrid().GetGrid())
      if (std::abs(x - xi) < eps12)
        {
          xid += eps10;
          break;
        }

    // Allocate needed operators (matching conditions and splitting
    // functions). By now the code is fast enough to precompute
    // everything at all available perturbative orders and the current
    // perturbative order is accounted for only when the actual
    // splitting functions and matching conditions (lambda) functions
    // are defined.
    // ===============================================================
    // LO Matching conditions
    std::map<int, Operator> MatchLO;
    const Operator Id  {g, Identity{}, IntEps, true};
    const Operator Zero{g, Null{},     IntEps, true};
    MatchLO.insert({MatchingBasisQCD::M0,  Id});
    MatchLO.insert({MatchingBasisQCD::M1,  Zero});
    MatchLO.insert({MatchingBasisQCD::M2,  Zero});
    MatchLO.insert({MatchingBasisQCD::M3,  Zero});
    MatchLO.insert({MatchingBasisQCD::M4,  Zero});
    MatchLO.insert({MatchingBasisQCD::M5,  Zero});
    MatchLO.insert({MatchingBasisQCD::M6,  Zero});
    MatchLO.insert({MatchingBasisQCD::M7,  Zero});
    MatchLO.insert({MatchingBasisQCD::M8,  Zero});
    MatchLO.insert({MatchingBasisQCD::M9,  Zero});
    MatchLO.insert({MatchingBasisQCD::M10, Zero});
    MatchLO.insert({MatchingBasisQCD::M11, Zero});

    // ===============================================================
    // LO splitting function operators
    std::map<int, std::map<int, Operator>> OpMapLO;
    const Operator O0ns{g, Pgpd0ns{xid}, IntEps, true};
    const Operator O0qq{g, Pgpd0qq{xid}, IntEps, true};
    const Operator O0gq{g, Pgpd0gq{xid}, IntEps, true};
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O0qg{g, Pgpd0qg{nf, xid}, IntEps, true};
        const Operator O0gg{g, Pgpd0gg{nf, xid}, IntEps, true};
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O0qq});
        OM.insert({EvolutionBasisQCD::PNSM, O0ns});
        OM.insert({EvolutionBasisQCD::PNSV, O0ns});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O0qq});
        OM.insert({EvolutionBasisQCD::PQG,                O0qg});
        OM.insert({EvolutionBasisQCD::PGQ,  ( nf / 6. ) * O0gq});
        OM.insert({EvolutionBasisQCD::PGG,                O0gg});
        OpMapLO.insert({nf, OM});
      }

    // ===============================================================
    // NLO matching conditions. The terms not multiplying
    // ln(muth2 / m2) of M_gh and M_hh are only known for parton-level
    // skewness smaller than one and are set to zero above.
    std::map<int, std::map<int, Operator>> MatchNLO;
    const Operator AS1HgL {g, Pgpd0qg{1, xid}, IntEps, true};
    const Operator AS1ggHL{g, AS1ggH_L{},      IntEps, true};
    const Operator AS1gH0 {g, AS1gpdgH_0{xid}, IntEps, true};
    const Operator AS1HH0 {g, AS1gpdHH_0{xid}, IntEps, true};
    const Operator AS1HHL {g, AS1gpdHH_L{xid}, IntEps, true};
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator AS1Hg  =          LogKth[nf] * AS1HgL;
        const Operator AS1ggH =          LogKth[nf] * AS1ggHL;
        const Operator AS1gH  = AS1gH0 + LogKth[nf] * O0gq;
        const Operator AS1HH  = AS1HH0 + LogKth[nf] * AS1HHL;
        std::map<int, Operator> OM;
        OM.insert({MatchingBasisQCD::M0,  Zero});
        OM.insert({MatchingBasisQCD::M1,  AS1ggH});
        OM.insert({MatchingBasisQCD::M2,  AS1gH});
        OM.insert({MatchingBasisQCD::M3,  AS1gH});
        OM.insert({MatchingBasisQCD::M4,  AS1Hg});
        OM.insert({MatchingBasisQCD::M5,  AS1HH});
        OM.insert({MatchingBasisQCD::M6,  AS1HH});
        OM.insert({MatchingBasisQCD::M7,  Zero});
        OM.insert({MatchingBasisQCD::M8,  AS1Hg});
        OM.insert({MatchingBasisQCD::M9,  AS1HH});
        OM.insert({MatchingBasisQCD::M10, Zero});
        OM.insert({MatchingBasisQCD::M11, Zero});
        MatchNLO.insert({nf, OM});
      }

    // ===============================================================
    // NLO splitting function operators. Only the non-singlet kernels
    // are available. The NLO singlet entries (pure-singlet qq, qg,
    // gq, gg) are set to zero, so that singlet and gluon evolve at LO
    // only.
    std::map<int, std::map<int, Operator>> OpMapNLO;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O1nsp{g, Pgpd1nsp{nf, xid}, IntEps, true};
        const Operator O1nsm{g, Pgpd1nsm{nf, xid}, IntEps, true};
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O1nsp});
        OM.insert({EvolutionBasisQCD::PNSM, O1nsm});
        OM.insert({EvolutionBasisQCD::PNSV, O1nsm});
        OM.insert({EvolutionBasisQCD::PQQ,  Zero});
        OM.insert({EvolutionBasisQCD::PQG,  Zero});
        OM.insert({EvolutionBasisQCD::PGQ,  Zero});
        OM.insert({EvolutionBasisQCD::PGG,  Zero});
        OpMapNLO.insert({nf, OM});
      }

    // Define object of the structure containing the DglapObjects
    std::map<int, DglapObjects> DglapObj;

    // Allocate convolution maps for evolution and matching, and set
    // of operators.
    for (int nf = nfi; nf <= nff; nf++)
      {
        DglapObjects obj;
        obj.Threshold = (nf > 0 ? Thresholds[nf - 1] : 0);
        if (OpEvol)
          {
            std::map<int, Operator> MapUnity;
            for (auto const& coord : Gkj)
              MapUnity.insert({coord.second, (coord.first.first == coord.first.second ? Id : Zero)});
            obj.UnitySet = Set<Operator> {EvolutionOperatorBasisQCD{nf}, MapUnity};
            obj.SplittingFunctions.insert({0, Set<Operator>{EvolutionOperatorBasisQCD{nf}, OpMapLO.at(nf)}});
            obj.SplittingFunctions.insert({1, Set<Operator>{EvolutionOperatorBasisQCD{nf}, OpMapNLO.at(nf)}});
            obj.MatchingConditions.insert({0, Set<Operator>{MatchingOperatorBasisQCD{nf},  MatchLO}});
            obj.MatchingConditions.insert({1, Set<Operator>{MatchingOperatorBasisQCD{nf},  MatchNLO.at(nf)}});
          }
        else
          {
            obj.SplittingFunctions.insert({0, Set<Operator>{EvolutionBasisQCD{nf}, OpMapLO.at(nf)}});
            obj.SplittingFunctions.insert({1, Set<Operator>{EvolutionBasisQCD{nf}, OpMapNLO.at(nf)}});
            obj.MatchingConditions.insert({0, Set<Operator>{MatchingBasisQCD{nf},  MatchLO}});
            obj.MatchingConditions.insert({1, Set<Operator>{MatchingBasisQCD{nf},  MatchNLO.at(nf)}});
          }
        DglapObj.insert({nf,obj});
      }
    t.stop();

    return DglapObj;
  }

  //_____________________________________________________________________________
  std::map<int, DglapObjects> InitializeGpdObjects(Grid                const& g,
                                                   std::vector<double> const& Thresholds,
                                                   double              const& xi,
                                                   bool                const& OpEvol,
                                                   double              const& IntEps)
  {
    return InitializeGpdObjects(g, Thresholds, Thresholds, xi, OpEvol, IntEps);
  }

  //_____________________________________________________________________________
  std::map<int, DglapObjects> InitializeGpdObjectsPol(Grid                const& g,
                                                      std::vector<double> const&, // masses not needed at LO
                                                      std::vector<double> const& Thresholds,
                                                      double              const& xi,
                                                      bool                const& OpEvol,
                                                      double              const& IntEps)
  {
    report("Initializing DglapObjects for GPD longitudinally polarised evolution at xi = " + std::to_string(xi) + "... ");
    Timer t;

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the threshold
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
        nfi++;

    // If xi happens to coincide with one of the nodes of the x-space
    // grid, numerical instabilities may prevent the integrals from
    // converging. In this eventuality, displace xi by a small amount.
    double xid = xi;
    for (double x : g.GetJointGrid().GetGrid())
      if (std::abs(x - xi) < eps12)
        {
          xid += eps10;
          break;
        }

    // Allocate needed operators (matching conditions and splitting
    // functions). By now the code is fast enough to precompute
    // everything at all available perturbative orders and the current
    // perturbative order is accounted for only when the actual
    // splitting functions and matching conditions (lambda) functions
    // are defined.
    // ===============================================================
    // LO Matching conditions
    std::map<int, Operator> MatchLO;
    const Operator Id  {g, Identity{}, IntEps, true};
    const Operator Zero{g, Null{},     IntEps, true};
    MatchLO.insert({MatchingBasisQCD::M0,  Id});
    MatchLO.insert({MatchingBasisQCD::M1,  Zero});
    MatchLO.insert({MatchingBasisQCD::M2,  Zero});
    MatchLO.insert({MatchingBasisQCD::M3,  Zero});
    MatchLO.insert({MatchingBasisQCD::M4,  Zero});
    MatchLO.insert({MatchingBasisQCD::M5,  Zero});
    MatchLO.insert({MatchingBasisQCD::M6,  Zero});
    MatchLO.insert({MatchingBasisQCD::M7,  Zero});
    MatchLO.insert({MatchingBasisQCD::M8,  Zero});
    MatchLO.insert({MatchingBasisQCD::M9,  Zero});
    MatchLO.insert({MatchingBasisQCD::M10, Zero});
    MatchLO.insert({MatchingBasisQCD::M11, Zero});

    // ===============================================================
    // LO splitting function operators
    std::map<int, std::map<int, Operator>> OpMapLO;
    const Operator O0ns{g, Pgpd0polns{xid}, IntEps, true};
    const Operator O0qq{g, Pgpd0polqq{xid}, IntEps, true};
    const Operator O0gq{g, Pgpd0polgq{xid}, IntEps, true};
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O0qg{g, Pgpd0polqg{nf, xid}, IntEps, true};
        const Operator O0gg{g, Pgpd0polgg{nf, xid}, IntEps, true};
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O0qq});
        OM.insert({EvolutionBasisQCD::PNSM, O0ns});
        OM.insert({EvolutionBasisQCD::PNSV, O0ns});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O0qq});
        OM.insert({EvolutionBasisQCD::PQG,                O0qg});
        OM.insert({EvolutionBasisQCD::PGQ,  ( nf / 6. ) * O0gq});
        OM.insert({EvolutionBasisQCD::PGG,                O0gg});
        OpMapLO.insert({nf, OM});
      }

    // Define object of the structure containing the DglapObjects
    std::map<int, DglapObjects> DglapObj;

    // Allocate convolution maps for evolution and matching, and set
    // of operators.
    for (int nf = nfi; nf <= nff; nf++)
      {
        DglapObjects obj;
        obj.Threshold = (nf > 0 ? Thresholds[nf - 1] : 0);
        if (OpEvol)
          {
            std::map<int, Operator> MapUnity;
            for (auto const& coord : Gkj)
              MapUnity.insert({coord.second, (coord.first.first == coord.first.second ? Id : Zero)});
            obj.UnitySet = Set<Operator> {EvolutionOperatorBasisQCD{nf}, MapUnity};
            obj.SplittingFunctions.insert({0, Set<Operator>{EvolutionOperatorBasisQCD{nf}, OpMapLO.at(nf)}});
            obj.MatchingConditions.insert({0, Set<Operator>{MatchingOperatorBasisQCD{nf},  MatchLO}});
          }
        else
          {
            obj.SplittingFunctions.insert({0, Set<Operator>{EvolutionBasisQCD{nf}, OpMapLO.at(nf)}});
            obj.MatchingConditions.insert({0, Set<Operator>{MatchingBasisQCD{nf},  MatchLO}});
          }
        DglapObj.insert({nf,obj});
      }
    t.stop();

    return DglapObj;
  }

  //_____________________________________________________________________________
  std::map<int, DglapObjects> InitializeGpdObjectsPol(Grid                const& g,
                                                      std::vector<double> const& Thresholds,
                                                      double              const& xi,
                                                      bool                const& OpEvol,
                                                      double              const& IntEps)
  {
    return InitializeGpdObjectsPol(g, Thresholds, Thresholds, xi, OpEvol, IntEps);
  }

  //_____________________________________________________________________________
  std::map<int, DglapObjects> InitializeGpdObjectsTrans(Grid                const& g,
                                                        std::vector<double> const&, // masses not needed at LO
                                                        std::vector<double> const& Thresholds,
                                                        double              const& xi,
                                                        bool                const& OpEvol,
                                                        double              const& IntEps)
  {
    report("Initializing DglapObjects for GPD transversely polarised evolution at xi = " + std::to_string(xi) + "... ");
    Timer t;

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the threshold
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
        nfi++;

    // If xi happens to coincide with one of the nodes of the x-space
    // grid, numerical instabilities may prevent the integrals from
    // converging. In this eventuality, displace xi by a small amount.
    double xid = xi;
    for (double x : g.GetJointGrid().GetGrid())
      if (std::abs(x - xi) < eps12)
        {
          xid += eps10;
          break;
        }

    // Allocate needed operators (matching conditions and splitting
    // functions). By now the code is fast enough to precompute
    // everything at all available perturbative orders and the current
    // perturbative order is accounted for only when the actual
    // splitting functions and matching conditions (lambda) functions
    // are defined.
    // ===============================================================
    // LO Matching conditions
    std::map<int, Operator> MatchLO;
    const Operator Id  {g, Identity{}, IntEps, true};
    const Operator Zero{g, Null{},     IntEps, true};
    MatchLO.insert({MatchingBasisQCD::M0,  Id});
    MatchLO.insert({MatchingBasisQCD::M1,  Zero});
    MatchLO.insert({MatchingBasisQCD::M2,  Zero});
    MatchLO.insert({MatchingBasisQCD::M3,  Zero});
    MatchLO.insert({MatchingBasisQCD::M4,  Zero});
    MatchLO.insert({MatchingBasisQCD::M5,  Zero});
    MatchLO.insert({MatchingBasisQCD::M6,  Zero});
    MatchLO.insert({MatchingBasisQCD::M7,  Zero});
    MatchLO.insert({MatchingBasisQCD::M8,  Zero});
    MatchLO.insert({MatchingBasisQCD::M9,  Zero});
    MatchLO.insert({MatchingBasisQCD::M10, Zero});
    MatchLO.insert({MatchingBasisQCD::M11, Zero});

    // ===============================================================
    // LO splitting function operators
    std::map<int, std::map<int, Operator>> OpMapLO;
    const Operator O0ns{g, Pgpd0transns{xid}, IntEps, true};
    const Operator O0qq{g, Pgpd0transqq{xid}, IntEps, true};
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O0gg{g, Pgpd0transgg{nf, xid}, IntEps, true};
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O0qq});
        OM.insert({EvolutionBasisQCD::PNSM, O0ns});
        OM.insert({EvolutionBasisQCD::PNSV, O0ns});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O0qq});
        OM.insert({EvolutionBasisQCD::PQG,                Zero});
        OM.insert({EvolutionBasisQCD::PGQ,                Zero});
        OM.insert({EvolutionBasisQCD::PGG,                O0gg});
        OpMapLO.insert({nf, OM});
      }

    // Define object of the structure containing the DglapObjects
    std::map<int, DglapObjects> DglapObj;

    // Allocate convolution maps for evolution and matching, and set
    // of operators.
    for (int nf = nfi; nf <= nff; nf++)
      {
        DglapObjects obj;
        obj.Threshold = (nf > 0 ? Thresholds[nf - 1] : 0);
        if (OpEvol)
          {
            std::map<int, Operator> MapUnity;
            for (auto const& coord : Gkj)
              MapUnity.insert({coord.second, (coord.first.first == coord.first.second ? Id : Zero)});
            obj.UnitySet = Set<Operator> {EvolutionOperatorBasisQCD{nf}, MapUnity};
            obj.SplittingFunctions.insert({0, Set<Operator>{EvolutionOperatorBasisQCD{nf}, OpMapLO.at(nf)}});
            obj.MatchingConditions.insert({0, Set<Operator>{MatchingOperatorBasisQCD{nf},  MatchLO}});
          }
        else
          {
            obj.SplittingFunctions.insert({0, Set<Operator>{EvolutionBasisQCD{nf}, OpMapLO.at(nf)}});
            obj.MatchingConditions.insert({0, Set<Operator>{MatchingBasisQCD{nf},  MatchLO}});
          }
        DglapObj.insert({nf,obj});
      }
    t.stop();

    return DglapObj;
  }

  //_____________________________________________________________________________
  std::map<int, DglapObjects> InitializeGpdObjectsTrans(Grid                const& g,
                                                        std::vector<double> const& Thresholds,
                                                        double              const& xi,
                                                        bool                const& OpEvol,
                                                        double              const& IntEps)
  {
    return InitializeGpdObjectsTrans(g, Thresholds, Thresholds, xi, OpEvol, IntEps);
  }
}
