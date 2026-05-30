//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/factorisationschemes.h"
#include "apfel/factorisationschemekernels.h"
#include "apfel/evolutionbasisqcd.h"
#include "apfel/matchingbasisqcd.h"
#include "apfel/betaqcd.h"
#include "apfel/messages.h"
#include "apfel/timer.h"

namespace apfel
{
  //_____________________________________________________________________________
  std::map<int, std::map<int, Set<Operator>>> InitializeSchemeChangeKernelsKrk(Grid                const& g,
                                                                               std::vector<double> const& Thresholds,
                                                                               double              const& IntEps)
  {
    report("Initializing scheme-change kernels from MSbar to Krk for space-like QCD unpolarised evolution... ");
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
    // NLO scheme-change-kernel operators
    std::map<int, std::map<int, Operator>> KNLO;
    const Operator K1ns{g, K1qqKrk{}, IntEps};
    const Operator K1qg{g, K1qgKrk{}, IntEps};
    const Operator K1gq{g, K1gqKrk{}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator K1gg{g, K1ggKrk{nf}, IntEps};
        std::map<int, Operator> OK;
        OK.insert({EvolutionBasisQCD::PNSP, K1ns});
        OK.insert({EvolutionBasisQCD::PNSM, K1ns});
        OK.insert({EvolutionBasisQCD::PNSV, K1ns});
        OK.insert({EvolutionBasisQCD::PQQ,  K1ns});
        OK.insert({EvolutionBasisQCD::PQG,  nf * K1qg});
        OK.insert({EvolutionBasisQCD::PGQ,  K1gq});
        OK.insert({EvolutionBasisQCD::PGG,  K1gg});
        KNLO.insert({nf, OK});
      }

    // Define object of the structure containing the kernels
    std::map<int, std::map<int, Set<Operator>>> KObj;

    // Run over active flavours
    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Set<Operator>> Knf;
        Knf.insert({1, Set<Operator>{EvolutionBasisQCD{nf}, KNLO.at(nf)}});
        KObj.insert({nf, Knf});
      }
    t.stop();

    return KObj;
  }

  //_____________________________________________________________________________
  std::map<int, DglapObjects> ChangeFactorisationSchemeMSbarToK(std::map<int, DglapObjects>                 const& DOMSbar,
                                                                std::map<int, std::map<int, Set<Operator>>> const& K)
  {
    // Set output equal to input
    std::map<int, DglapObjects> DOK = DOMSbar;

    // Adjust NLO splitting functions and matching conditions
    for (auto const& Knf : K)
      {
        const int nf = Knf.first;
        const std::map<int, Operator> KnfNLO = Knf.second.at(1).GetObjects();
        const std::map<int, Operator> PMSLO  = DOMSbar.at(nf).SplittingFunctions.at(0).GetObjects();
        const std::map<int, Operator> PMSNLO = DOMSbar.at(nf).SplittingFunctions.at(1).GetObjects();
        const std::map<int, Operator> MMSNLO = DOMSbar.at(nf).MatchingConditions.at(1).GetObjects();

        // NLO splitting functions tranformed
        std::map<int, Operator> PKNLO;
        PKNLO.insert({EvolutionBasisQCD::PNSP, PMSNLO.at(EvolutionBasisQCD::PNSP) - beta0qcd(nf) * KnfNLO.at(EvolutionBasisQCD::PNSP)});
        PKNLO.insert({EvolutionBasisQCD::PNSM, PMSNLO.at(EvolutionBasisQCD::PNSM) - beta0qcd(nf) * KnfNLO.at(EvolutionBasisQCD::PNSM)});
        PKNLO.insert({EvolutionBasisQCD::PNSV, PMSNLO.at(EvolutionBasisQCD::PNSV) - beta0qcd(nf) * KnfNLO.at(EvolutionBasisQCD::PNSV)});
        PKNLO.insert({EvolutionBasisQCD::PQQ, PMSNLO.at(EvolutionBasisQCD::PQQ) - beta0qcd(nf) * KnfNLO.at(EvolutionBasisQCD::PQQ)
                      + KnfNLO.at(EvolutionBasisQCD::PQQ) * PMSLO.at(EvolutionBasisQCD::PQQ) + KnfNLO.at(EvolutionBasisQCD::PQG) * PMSLO.at(EvolutionBasisQCD::PGQ)
                      - PMSLO.at(EvolutionBasisQCD::PQQ) * KnfNLO.at(EvolutionBasisQCD::PQQ) + PMSLO.at(EvolutionBasisQCD::PQG) * KnfNLO.at(EvolutionBasisQCD::PGQ)});
        PKNLO.insert({EvolutionBasisQCD::PQG, PMSNLO.at(EvolutionBasisQCD::PQG) - beta0qcd(nf) * KnfNLO.at(EvolutionBasisQCD::PQG)
                      + KnfNLO.at(EvolutionBasisQCD::PQQ) * PMSLO.at(EvolutionBasisQCD::PQG) + KnfNLO.at(EvolutionBasisQCD::PQG) * PMSLO.at(EvolutionBasisQCD::PGG)
                      - PMSLO.at(EvolutionBasisQCD::PQQ) * KnfNLO.at(EvolutionBasisQCD::PQG) + PMSLO.at(EvolutionBasisQCD::PQG) * KnfNLO.at(EvolutionBasisQCD::PGG)});
        PKNLO.insert({EvolutionBasisQCD::PGQ, PMSNLO.at(EvolutionBasisQCD::PGQ) - beta0qcd(nf) * KnfNLO.at(EvolutionBasisQCD::PGQ)
                      + KnfNLO.at(EvolutionBasisQCD::PGQ) * PMSLO.at(EvolutionBasisQCD::PQQ) + KnfNLO.at(EvolutionBasisQCD::PGG) * PMSLO.at(EvolutionBasisQCD::PGQ)
                      - PMSLO.at(EvolutionBasisQCD::PGQ) * KnfNLO.at(EvolutionBasisQCD::PQQ) + PMSLO.at(EvolutionBasisQCD::PGG) * KnfNLO.at(EvolutionBasisQCD::PGQ)});
        PKNLO.insert({EvolutionBasisQCD::PGG, PMSNLO.at(EvolutionBasisQCD::PGG) - beta0qcd(nf) * KnfNLO.at(EvolutionBasisQCD::PGG)
                      + KnfNLO.at(EvolutionBasisQCD::PGQ) * PMSLO.at(EvolutionBasisQCD::PQG) + KnfNLO.at(EvolutionBasisQCD::PGG) * PMSLO.at(EvolutionBasisQCD::PGG)
                      - PMSLO.at(EvolutionBasisQCD::PGQ) * KnfNLO.at(EvolutionBasisQCD::PQG) + PMSLO.at(EvolutionBasisQCD::PGG) * KnfNLO.at(EvolutionBasisQCD::PGG)});

        // Replace NLO splitting functions
        DOK.at(nf).SplittingFunctions.at(1).SetObjects(PKNLO);

        // NLO matching conditions tranformed
        std::map<int, Operator> MKNLO;
        if (K.contains(nf+1))
          {
            MKNLO.insert({MatchingBasisQCD::M0,  MMSNLO.at(MatchingBasisQCD::M0)});
            MKNLO.insert({MatchingBasisQCD::M1,  MMSNLO.at(MatchingBasisQCD::M1) + K.at(nf+1).at(1).at(EvolutionBasisQCD::PGG) - KnfNLO.at(EvolutionBasisQCD::PGG)});
            MKNLO.insert({MatchingBasisQCD::M2,  MMSNLO.at(MatchingBasisQCD::M2)});
            MKNLO.insert({MatchingBasisQCD::M3,  MMSNLO.at(MatchingBasisQCD::M3)});
            MKNLO.insert({MatchingBasisQCD::M4,  MMSNLO.at(MatchingBasisQCD::M4)});
            MKNLO.insert({MatchingBasisQCD::M5,  MMSNLO.at(MatchingBasisQCD::M5)});
            MKNLO.insert({MatchingBasisQCD::M6,  MMSNLO.at(MatchingBasisQCD::M6)});
            MKNLO.insert({MatchingBasisQCD::M7,  MMSNLO.at(MatchingBasisQCD::M7)});
            MKNLO.insert({MatchingBasisQCD::M8,  MMSNLO.at(MatchingBasisQCD::M8)});
            MKNLO.insert({MatchingBasisQCD::M9,  MMSNLO.at(MatchingBasisQCD::M9)});
            MKNLO.insert({MatchingBasisQCD::M10, MMSNLO.at(MatchingBasisQCD::M10)});
            MKNLO.insert({MatchingBasisQCD::M11, MMSNLO.at(MatchingBasisQCD::M11)});
          }

        // Replace NLO matching conditions
        DOK.at(nf).MatchingConditions.at(1).SetObjects(MKNLO);
      }

    return DOK;
  }
}
