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
#include "apfel/evolutionbasisqcd.h"
#include "apfel/matchingbasisqcd.h"

namespace apfel
{
  //_____________________________________________________________________________
  std::map<int, DglapObjects> InitializeGpdObjects(Grid                const& g,
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

    // If xi happens to coincide with one of the nodes of the x-space
    // grid, numerical instabilities may prevent the integrals from
    // converging. In this eventuality, displace xi by a small amount.
    double xid = xi;
    for (double x : g.GetJointGrid().GetGrid())
      if (std::abs(x - xi) > eps12)
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

    // Define object of the structure containing the DglapObjects
    std::map<int, DglapObjects> DglapObj;

    // Allocate convolution maps for evolution and matching, and set
    // of operators.
    for (int nf = nfi; nf <= nff; nf++)
      {
        DglapObjects obj;
        obj.Threshold = Thresholds[nf-1];
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
      if (std::abs(x - xi) > eps12)
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
        obj.Threshold = Thresholds[nf-1];
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
      if (std::abs(x - xi) > eps12)
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
        obj.Threshold = Thresholds[nf-1];
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
}
