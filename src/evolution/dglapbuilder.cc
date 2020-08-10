//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/dglapbuilder.h"
#include "apfel/timer.h"
#include "apfel/constants.h"
#include "apfel/tools.h"
#include "apfel/messages.h"
#include "apfel/splittingfunctions.h"
#include "apfel/splittingfunctions_tl.h"
#include "apfel/matchingconditions.h"
#include "apfel/matchingconditions_tl.h"
#include "apfel/evolutionbasisqcd.h"
#include "apfel/matchingbasisqcd.h"

namespace apfel
{
  //_____________________________________________________________________________
  std::map<int, DglapObjects> InitializeDglapObjectsQCD(Grid                const& g,
                                                        std::vector<double> const& Masses,
                                                        std::vector<double> const& Thresholds,
                                                        bool                const& OpEvol,
                                                        double              const& IntEps)
  {
    report("Initializing DglapObjects for space-like QCD unpolarised evolution... ");
    Timer t;

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the thresholds
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
        nfi++;

    // Compute logs of muth2 / m2 needed for the the matching
    // conditions.
    std::vector<double> LogKth;
    for (int im = 0; im < (int) Thresholds.size(); im++)
      if (Thresholds[im] < eps12 || Masses[im] < eps12)
        LogKth.push_back(0);
      else
        LogKth.push_back(2 * log( Thresholds[im] / Masses[im] ));

    // Allocate needed operators (matching conditions and splitting
    // functions). By now the code is fast enough to precompute
    // everything at all available perturbative orders and the current
    // perturbative order is accounted for only when the actual
    // splitting functions and matching conditions (lambda) functions
    // are defined.
    // ===============================================================
    // LO Matching conditions.
    std::map<int, Operator> MatchLO;
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};
    MatchLO.insert({MatchingBasisQCD::M0, Id});
    MatchLO.insert({MatchingBasisQCD::M1, Zero});
    MatchLO.insert({MatchingBasisQCD::M2, Zero});
    MatchLO.insert({MatchingBasisQCD::M3, Zero});
    MatchLO.insert({MatchingBasisQCD::M4, Zero});
    MatchLO.insert({MatchingBasisQCD::M5, Zero});
    MatchLO.insert({MatchingBasisQCD::M6, Zero});
    MatchLO.insert({MatchingBasisQCD::M7, Zero});

    // ===============================================================
    // LO splitting function operators.
    std::map<int, std::map<int, Operator>> OpMapLO;
    const Operator O0ns{g, P0ns{}, IntEps};
    const Operator O0qg{g, P0qg{}, IntEps};
    const Operator O0gq{g, P0gq{}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O0gg{g, P0gg{nf}, IntEps};
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O0ns});
        OM.insert({EvolutionBasisQCD::PNSM, O0ns});
        OM.insert({EvolutionBasisQCD::PNSV, O0ns});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O0ns});
        OM.insert({EvolutionBasisQCD::PQG,           nf * O0qg});
        OM.insert({EvolutionBasisQCD::PGQ,  ( nf / 6. ) * O0gq});
        OM.insert({EvolutionBasisQCD::PGG,                O0gg});
        OpMapLO.insert({nf, OM});
      }

    // ===============================================================
    // NLO Matching conditions.
    std::map<int, std::map<int, Operator>> MatchNLO;
    const Operator AS1HgL {g, AS1Hg_L{},  IntEps};
    const Operator AS1ggHL{g, AS1ggH_L{}, IntEps};
    const Operator AS1gH0 {g, AS1gH_0{},  IntEps};
    const Operator AS1gHL {g, AS1gH_L{},  IntEps};
    const Operator AS1HH0 {g, AS1HH_0{},  IntEps};
    const Operator AS1HHL {g, AS1HH_L{},  IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator AS1Hg  =          LogKth[nf-1] * AS1HgL;
        const Operator AS1ggH =          LogKth[nf-1] * AS1ggHL;
        const Operator AS1gH  = AS1gH0 + LogKth[nf-1] * AS1gHL;
        const Operator AS1HH  = AS1HH0 + LogKth[nf-1] * AS1HHL;
        std::map<int, Operator> OM;
        OM.insert({MatchingBasisQCD::M0, Zero});
        OM.insert({MatchingBasisQCD::M1, AS1ggH});
        OM.insert({MatchingBasisQCD::M2, AS1gH});
        OM.insert({MatchingBasisQCD::M3, AS1gH});
        OM.insert({MatchingBasisQCD::M4, AS1Hg});
        OM.insert({MatchingBasisQCD::M5, AS1HH});
        OM.insert({MatchingBasisQCD::M6, AS1HH});
        OM.insert({MatchingBasisQCD::M7, Zero});
        MatchNLO.insert({nf, OM});
      }

    // ===============================================================
    // NLO splitting function operators.
    std::map<int, std::map<int, Operator>> OpMapNLO;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O1nsp{g, P1nsp{nf}, IntEps};
        const Operator O1nsm{g, P1nsm{nf}, IntEps};
        const Operator O1ps {g, P1ps{nf},  IntEps};
        const Operator O1qg {g, P1qg{nf},  IntEps};
        const Operator O1gq {g, P1gq{nf},  IntEps};
        const Operator O1gg {g, P1gg{nf},  IntEps};
        const Operator O1qq = O1nsp + O1ps;
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O1nsp});
        OM.insert({EvolutionBasisQCD::PNSM, O1nsm});
        OM.insert({EvolutionBasisQCD::PNSV, O1nsm});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O1qq});
        OM.insert({EvolutionBasisQCD::PQG,                O1qg});
        OM.insert({EvolutionBasisQCD::PGQ,  ( nf / 6. ) * O1gq});
        OM.insert({EvolutionBasisQCD::PGG,                O1gg});
        OpMapNLO.insert({nf, OM});
      }

    // ===============================================================
    // NNLO Matching conditions.
    std::map<int, std::map<int, Operator>> MatchNNLO;
    const Operator APS2Hq0  {g, APS2Hq_0{},   IntEps};
    const Operator APS2HqL  {g, APS2Hq_L{},   IntEps};
    const Operator APS2HqL2 {g, APS2Hq_L2{},  IntEps};
    const Operator ANS2qqH0 {g, ANS2qqH_0{},  IntEps};
    const Operator ANS2qqHL {g, ANS2qqH_L{},  IntEps};
    const Operator ANS2qqHL2{g, ANS2qqH_L2{}, IntEps};
    const Operator AS2Hg0   {g, AS2Hg_0{},    IntEps};
    const Operator AS2HgL   {g, AS2Hg_L{},    IntEps};
    const Operator AS2HgL2  {g, AS2Hg_L2{},   IntEps};
    const Operator AS2gqH0  {g, AS2gqH_0{},   IntEps};
    const Operator AS2gqHL  {g, AS2gqH_L{},   IntEps};
    const Operator AS2gqHL2 {g, AS2gqH_L2{},  IntEps};
    const Operator AS2ggH0  {g, AS2ggH_0{},   IntEps};
    const Operator AS2ggHL  {g, AS2ggH_L{},   IntEps};
    const Operator AS2ggHL2 {g, AS2ggH_L2{},  IntEps};
    const Operator AS2qqH0  = ANS2qqH0  + APS2Hq0;
    const Operator AS2qqHL  = ANS2qqHL  + APS2HqL;
    const Operator AS2qqHL2 = ANS2qqHL2 + APS2HqL2;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const double lnk  = LogKth[nf-1];
        const double lnk2 = lnk * lnk;
        const Operator APS2Hq  = APS2Hq0  + lnk * APS2HqL  + lnk2 * APS2HqL2;
        const Operator ANS2qqH = ANS2qqH0 + lnk * ANS2qqHL + lnk2 * ANS2qqHL2;
        const Operator AS2Hg   = AS2Hg0   + lnk * AS2HgL   + lnk2 * AS2HgL2;
        const Operator AS2gqH  = AS2gqH0  + lnk * AS2gqHL  + lnk2 * AS2gqHL2;
        const Operator AS2ggH  = AS2ggH0  + lnk * AS2ggHL  + lnk2 * AS2ggHL2;
        const Operator AS2qqH  = AS2qqH0  + lnk * AS2qqHL  + lnk2 * AS2qqHL2;
        std::map<int, Operator> OM;
        OM.insert({MatchingBasisQCD::M0, Zero});
        OM.insert({MatchingBasisQCD::M1, AS2ggH});
        OM.insert({MatchingBasisQCD::M2, nf * AS2gqH});
        OM.insert({MatchingBasisQCD::M3, (-1) * AS2gqH});
        OM.insert({MatchingBasisQCD::M4, AS2Hg});
        OM.insert({MatchingBasisQCD::M5, nf * AS2qqH});
        OM.insert({MatchingBasisQCD::M6, (-1) * AS2qqH});
        OM.insert({MatchingBasisQCD::M7, ANS2qqH});
        MatchNNLO.insert({nf, OM});
      }

    // Auxiliary NNLO contributions to be used for backward
    // matching. The are essentially the square of the NLO matching
    // matrix. They are labelled with perturbative order -2.
    std::map<int, std::map<int, Operator>> MatchNNLOb;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator AS1Hg  =          LogKth[nf-1] * AS1HgL;
        const Operator AS1ggH =          LogKth[nf-1] * AS1ggHL;
        const Operator AS1gH  = AS1gH0 + LogKth[nf-1] * AS1gHL;
        const Operator AS1HH  = AS1HH0 + LogKth[nf-1] * AS1HHL;
        const Operator AS1Hg2  = AS1Hg  * AS1ggH + AS1gH * AS1HH;;
        const Operator AS1ggH2 = AS1ggH * AS1ggH + AS1gH * AS1Hg;
        const Operator AS1gH2  = AS1ggH * AS1gH  + AS1gH * AS1HH;
        const Operator AS1HH2  = AS1Hg  * AS1gH  + AS1HH * AS1HH;
        std::map<int, Operator> OM;
        OM.insert({MatchingBasisQCD::M0, Zero});
        OM.insert({MatchingBasisQCD::M1, AS1ggH2});
        OM.insert({MatchingBasisQCD::M2, AS1gH2});
        OM.insert({MatchingBasisQCD::M3, AS1gH2});
        OM.insert({MatchingBasisQCD::M4, AS1Hg2});
        OM.insert({MatchingBasisQCD::M5, AS1HH2});
        OM.insert({MatchingBasisQCD::M6, AS1HH2});
        OM.insert({MatchingBasisQCD::M7, Zero});
        MatchNNLOb.insert({nf, OM});
      }

    // ===============================================================
    // NNLO splitting function operators.
    std::map<int, std::map<int, Operator>> OpMapNNLO;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O2nsp{g, P2nsp{nf}, IntEps};
        const Operator O2nsm{g, P2nsm{nf}, IntEps};
        const Operator O2nss{g, P2nss{nf}, IntEps};
        const Operator O2ps {g, P2ps{nf},  IntEps};
        const Operator O2qg {g, P2qg{nf},  IntEps};
        const Operator O2gq {g, P2gq{nf},  IntEps};
        const Operator O2gg {g, P2gg{nf},  IntEps};
        const Operator O2qq  = O2nsp + O2ps;
        const Operator O2nsv = O2nsm + O2nss;
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O2nsp});
        OM.insert({EvolutionBasisQCD::PNSM, O2nsm});
        OM.insert({EvolutionBasisQCD::PNSV, O2nsv});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O2qq});
        OM.insert({EvolutionBasisQCD::PQG,                O2qg});
        OM.insert({EvolutionBasisQCD::PGQ,  ( nf / 6. ) * O2gq});
        OM.insert({EvolutionBasisQCD::PGG,                O2gg});
        OpMapNNLO.insert({nf, OM});
      }

    // ===============================================================
    // NNNLO splitting function operators. For now only the
    // non-singlet splitting functions have been computed to leading
    // colour even though the subleading colour part is estimated
    // through an approximate parameterisation.
    std::map<int, std::map<int, Operator>> OpMapNNNLO;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O3nsp{g, P3nsp{nf}, IntEps};
        const Operator O3nsm{g, P3nsm{nf}, IntEps};
        const Operator O3nss{g, P3nss{nf}, IntEps};
        const Operator O3nsv = O3nsm + O3nss;
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O3nsp});
        OM.insert({EvolutionBasisQCD::PNSM, O3nsm});
        OM.insert({EvolutionBasisQCD::PNSV, O3nsv});
        OM.insert({EvolutionBasisQCD::PQQ,  Zero});
        OM.insert({EvolutionBasisQCD::PQG,  Zero});
        OM.insert({EvolutionBasisQCD::PGQ,  Zero});
        OM.insert({EvolutionBasisQCD::PGG,  Zero});
        OpMapNNNLO.insert({nf, OM});
      }

    // Define object of the structure containing the DglapObjects.
    std::map<int, DglapObjects> DglapObj;

    // Allocate convolution maps for evolution and matching, and set
    // of operators.
    for (int nf = nfi; nf <= nff; nf++)
      {
        DglapObjects obj;
        obj.Threshold = Thresholds[nf-1];
        if (OpEvol)
          {
            obj.SplittingFunctions.insert({ 0, Set<Operator>{EvolutionOperatorBasisQCD{nf}, OpMapLO.at(nf)}});
            obj.SplittingFunctions.insert({ 1, Set<Operator>{EvolutionOperatorBasisQCD{nf}, OpMapNLO.at(nf)}});
            obj.SplittingFunctions.insert({ 2, Set<Operator>{EvolutionOperatorBasisQCD{nf}, OpMapNNLO.at(nf)}});
            obj.SplittingFunctions.insert({ 3, Set<Operator>{EvolutionOperatorBasisQCD{nf}, OpMapNNNLO.at(nf)}});
            obj.MatchingConditions.insert({ 0, Set<Operator>{MatchingOperatorBasisQCD{nf},  MatchLO}});
            obj.MatchingConditions.insert({ 1, Set<Operator>{MatchingOperatorBasisQCD{nf},  MatchNLO.at(nf)}});
            obj.MatchingConditions.insert({ 2, Set<Operator>{MatchingOperatorBasisQCD{nf},  MatchNNLO.at(nf)}});
            obj.MatchingConditions.insert({-2, Set<Operator>{MatchingOperatorBasisQCD{nf},  MatchNNLOb.at(nf)}});
          }
        else
          {
            obj.SplittingFunctions.insert({ 0, Set<Operator>{EvolutionBasisQCD{nf}, OpMapLO.at(nf)}});
            obj.SplittingFunctions.insert({ 1, Set<Operator>{EvolutionBasisQCD{nf}, OpMapNLO.at(nf)}});
            obj.SplittingFunctions.insert({ 2, Set<Operator>{EvolutionBasisQCD{nf}, OpMapNNLO.at(nf)}});
            obj.SplittingFunctions.insert({ 3, Set<Operator>{EvolutionBasisQCD{nf}, OpMapNNNLO.at(nf)}});
            obj.MatchingConditions.insert({ 0, Set<Operator>{MatchingBasisQCD{nf},  MatchLO}});
            obj.MatchingConditions.insert({ 1, Set<Operator>{MatchingBasisQCD{nf},  MatchNLO.at(nf)}});
            obj.MatchingConditions.insert({ 2, Set<Operator>{MatchingBasisQCD{nf},  MatchNNLO.at(nf)}});
            obj.MatchingConditions.insert({-2, Set<Operator>{MatchingBasisQCD{nf},  MatchNNLOb.at(nf)}});
          }
        DglapObj.insert({nf, obj});
      }
    t.stop();

    return DglapObj;
  }

  //_____________________________________________________________________________
  std::map<int, DglapObjects> InitializeDglapObjectsQCD(Grid                const& g,
                                                        std::vector<double> const& Thresholds,
                                                        bool                const& OpEvol,
                                                        double              const& IntEps)
  {
    return InitializeDglapObjectsQCD(g, Thresholds, Thresholds, OpEvol, IntEps);
  }

  //_____________________________________________________________________________
  std::map<int, DglapObjects> InitializeDglapObjectsQCDT(Grid                const& g,
                                                         std::vector<double> const& Masses,
                                                         std::vector<double> const& Thresholds,
                                                         bool                const& OpEvol,
                                                         double              const& IntEps)
  {
    report("Initializing DglapObjects for time-like QCD unpolarised evolution... ");
    Timer t;

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the thresholds
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
        nfi++;

    // Compute logs of muth2 / m2 needed for the the matching
    // conditions.
    std::vector<double> LogKth;
    for (int im = 0; im < (int) Thresholds.size(); im++)
      if (Thresholds[im] < eps12 || Masses[im] < eps12)
        LogKth.push_back(0);
      else
        LogKth.push_back(2 * log( Thresholds[im] / Masses[im] ));

    // Allocate needed operators (matching conditions and splitting
    // functions). By now the code is fast enough to precompute
    // everything at all available perturbative orders and the current
    // perturbative order is accounted for only when the actual
    // splitting functions and matching conditions (lambda) functions
    // are defined.
    // ===============================================================
    // LO Matching conditions.
    std::map<int, Operator> MatchLO;
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};
    MatchLO.insert({MatchingBasisQCD::M0, Id});
    MatchLO.insert({MatchingBasisQCD::M1, Zero});
    MatchLO.insert({MatchingBasisQCD::M2, Zero});
    MatchLO.insert({MatchingBasisQCD::M3, Zero});
    MatchLO.insert({MatchingBasisQCD::M4, Zero});
    MatchLO.insert({MatchingBasisQCD::M5, Zero});
    MatchLO.insert({MatchingBasisQCD::M6, Zero});
    MatchLO.insert({MatchingBasisQCD::M7, Zero});

    // ===============================================================
    // LO splitting function operators.
    std::map<int, std::map<int, Operator>> OpMapLO;
    const Operator O0ns{g, P0Tns{}, IntEps};
    const Operator O0qg{g, P0Tqg{}, IntEps};
    const Operator O0gq{g, P0Tgq{}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O0gg{g, P0Tgg{nf}, IntEps};
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O0ns});
        OM.insert({EvolutionBasisQCD::PNSM, O0ns});
        OM.insert({EvolutionBasisQCD::PNSV, O0ns});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O0ns});
        OM.insert({EvolutionBasisQCD::PQG,           nf * O0qg});
        OM.insert({EvolutionBasisQCD::PGQ,  ( nf / 6. ) * O0gq});
        OM.insert({EvolutionBasisQCD::PGG,                O0gg});
        OpMapLO.insert({nf, OM});
      }

    // ===============================================================
    // NLO Matching conditions.
    std::map<int, std::map<int, Operator>> MatchNLO;
    const Operator AS1Hg0 {g, ATS1Hg_0{},  IntEps};
    const Operator AS1HgL {g, ATS1Hg_L{},  IntEps};
    const Operator AS1ggHL{g, ATS1ggH_L{}, IntEps};
    const Operator AS1gHL {g, ATS1gH_L{},  IntEps};
    const Operator AS1HH0 {g, ATS1HH_0{},  IntEps};
    const Operator AS1HHL {g, ATS1HH_L{},  IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator AS1Hg  = AS1Hg0 + LogKth[nf-1] * AS1HgL;
        const Operator AS1ggH =          LogKth[nf-1] * AS1ggHL;
        const Operator AS1gH  =          LogKth[nf-1] * AS1gHL;
        const Operator AS1HH  = AS1HH0 + LogKth[nf-1] * AS1HHL;
        std::map<int, Operator> OM;
        OM.insert({MatchingBasisQCD::M0, Zero});
        OM.insert({MatchingBasisQCD::M1, AS1ggH});
        OM.insert({MatchingBasisQCD::M2, AS1gH});
        OM.insert({MatchingBasisQCD::M3, AS1gH});
        OM.insert({MatchingBasisQCD::M4, AS1Hg});
        OM.insert({MatchingBasisQCD::M5, AS1HH});
        OM.insert({MatchingBasisQCD::M6, AS1HH});
        OM.insert({MatchingBasisQCD::M7, Zero});
        MatchNLO.insert({nf, OM});
      }

    // ===============================================================
    // NLO splitting function operators.
    std::map<int, std::map<int, Operator>> OpMapNLO;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O1nsp{g, P1Tnsp{nf}, IntEps};
        const Operator O1nsm{g, P1Tnsm{nf}, IntEps};
        const Operator O1ps {g, P1Tps{nf},  IntEps};
        const Operator O1qg {g, P1Tqg{nf},  IntEps};
        const Operator O1gq {g, P1Tgq{nf},  IntEps};
        const Operator O1gg {g, P1Tgg{nf},  IntEps};
        const Operator O1qq = O1nsp + O1ps;
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O1nsp});
        OM.insert({EvolutionBasisQCD::PNSM, O1nsm});
        OM.insert({EvolutionBasisQCD::PNSV, O1nsm});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O1qq});
        OM.insert({EvolutionBasisQCD::PQG,                O1qg});
        OM.insert({EvolutionBasisQCD::PGQ,  ( nf / 6. ) * O1gq});
        OM.insert({EvolutionBasisQCD::PGG,                O1gg});
        OpMapNLO.insert({nf, OM});
      }

    // ===============================================================
    // NNLO Matching conditions. Set to zero for now because they are
    // not know yet.
    std::map<int, Operator> MatchNNLO;
    for (int i = MatchingBasisQCD::M0; i <= MatchingBasisQCD::M7; i++)
      MatchNNLO.insert({i, Zero});

    // ===============================================================
    // NNLO splitting function operators.
    std::map<int, std::map<int, Operator>> OpMapNNLO;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O2nsp{g, P2Tnsp{nf}, IntEps};
        const Operator O2nsm{g, P2Tnsm{nf}, IntEps};
        const Operator O2nss{g, P2Tnss{nf}, IntEps};
        const Operator O2ps {g, P2Tps{nf},  IntEps};
        const Operator O2qg {g, P2Tqg{nf},  IntEps};
        const Operator O2gq {g, P2Tgq{nf},  IntEps};
        const Operator O2gg {g, P2Tgg{nf},  IntEps};
        const Operator O2qq  = O2nsp + O2ps;
        const Operator O2nsv = O2nsm + O2nss;
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O2nsp});
        OM.insert({EvolutionBasisQCD::PNSM, O2nsm});
        OM.insert({EvolutionBasisQCD::PNSV, O2nsv});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O2qq});
        OM.insert({EvolutionBasisQCD::PQG,                O2qg});
        OM.insert({EvolutionBasisQCD::PGQ,  ( nf / 6. ) * O2gq});
        OM.insert({EvolutionBasisQCD::PGG,                O2gg});
        OpMapNNLO.insert({nf, OM});
      }

    // Define object of the structure containing the DglapObjects.
    std::map<int, DglapObjects> DglapObj;

    // Allocate convolution maps for evolution and matching, and set
    // of operators.
    for (int nf = nfi; nf <= nff; nf++)
      {
        DglapObjects obj;
        obj.Threshold = Thresholds[nf-1];
        if (OpEvol)
          {
            obj.SplittingFunctions.insert({ 0, Set<Operator>{EvolutionOperatorBasisQCD{nf}, OpMapLO.at(nf)}});
            obj.SplittingFunctions.insert({ 1, Set<Operator>{EvolutionOperatorBasisQCD{nf}, OpMapNLO.at(nf)}});
            obj.SplittingFunctions.insert({ 2, Set<Operator>{EvolutionOperatorBasisQCD{nf}, OpMapNNLO.at(nf)}});
            obj.MatchingConditions.insert({ 0, Set<Operator>{MatchingOperatorBasisQCD{nf},  MatchLO}});
            obj.MatchingConditions.insert({ 1, Set<Operator>{MatchingOperatorBasisQCD{nf},  MatchNLO.at(nf)}});
            obj.MatchingConditions.insert({ 2, Set<Operator>{MatchingOperatorBasisQCD{nf},  MatchNNLO}});
            obj.MatchingConditions.insert({-2, Set<Operator>{MatchingOperatorBasisQCD{nf},  MatchNNLO}});
          }
        else
          {
            obj.SplittingFunctions.insert({ 0, Set<Operator>{EvolutionBasisQCD{nf}, OpMapLO.at(nf)}});
            obj.SplittingFunctions.insert({ 1, Set<Operator>{EvolutionBasisQCD{nf}, OpMapNLO.at(nf)}});
            obj.SplittingFunctions.insert({ 2, Set<Operator>{EvolutionBasisQCD{nf}, OpMapNNLO.at(nf)}});
            obj.MatchingConditions.insert({ 0, Set<Operator>{MatchingBasisQCD{nf},  MatchLO}});
            obj.MatchingConditions.insert({ 1, Set<Operator>{MatchingBasisQCD{nf},  MatchNLO.at(nf)}});
            obj.MatchingConditions.insert({ 2, Set<Operator>{MatchingBasisQCD{nf},  MatchNNLO}});
            obj.MatchingConditions.insert({-2, Set<Operator>{MatchingBasisQCD{nf},  MatchNNLO}});
          }
        DglapObj.insert({nf,obj});
      }
    t.stop();

    return DglapObj;
  }

  //_____________________________________________________________________________
  std::map<int, DglapObjects> InitializeDglapObjectsQCDT(Grid                const& g,
                                                         std::vector<double> const& Thresholds,
                                                         bool                const& OpEvol,
                                                         double              const& IntEps)
  {
    return InitializeDglapObjectsQCDT(g, Thresholds, Thresholds, OpEvol, IntEps);
  }

  //_____________________________________________________________________________
  std::map<int, DglapObjects> InitializeDglapObjectsQCDtrans(Grid                const& g,
                                                             std::vector<double> const& Masses,
                                                             std::vector<double> const& Thresholds,
                                                             bool                const& OpEvol,
                                                             double              const& IntEps)
  {
    report("Initializing DglapObjects for space-like QCD transversely polarised evolution... ");
    Timer t;

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the thresholds
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
        nfi++;

    // Compute logs of muth2 / m2 needed for the the matching
    // conditions (not possible for transversities yet).
    std::vector<double> LogKth;
    for (int im = 0; im < (int) Thresholds.size(); im++)
      if (Thresholds[im] < eps12 || Masses[im] < eps12)
        LogKth.push_back(0);
      else
        LogKth.push_back(2 * log( Thresholds[im] / Masses[im] ));

    // Allocate needed operators (matching conditions and splitting
    // functions). By now the code is fast enough to precompute
    // everything at all available perturbative orders and the current
    // perturbative order is accounted for only when the actual
    // splitting functions and matching conditions (lambda) functions
    // are defined.
    // ===============================================================
    // LO Matching conditions.
    std::map<int, Operator> MatchLO;
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};
    MatchLO.insert({MatchingBasisQCD::M0, Id});
    MatchLO.insert({MatchingBasisQCD::M1, Zero});
    MatchLO.insert({MatchingBasisQCD::M2, Zero});
    MatchLO.insert({MatchingBasisQCD::M3, Zero});
    MatchLO.insert({MatchingBasisQCD::M4, Zero});
    MatchLO.insert({MatchingBasisQCD::M5, Zero});
    MatchLO.insert({MatchingBasisQCD::M6, Zero});
    MatchLO.insert({MatchingBasisQCD::M7, Zero});

    // ===============================================================
    // LO splitting function operators.
    std::map<int, std::map<int, Operator>> OpMapLO;
    const Operator O0ns{g, P0transns{}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O0ns});
        OM.insert({EvolutionBasisQCD::PNSM, O0ns});
        OM.insert({EvolutionBasisQCD::PNSV, O0ns});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O0ns});
        OM.insert({EvolutionBasisQCD::PQG,                Zero});
        OM.insert({EvolutionBasisQCD::PGQ,                Zero});
        OM.insert({EvolutionBasisQCD::PGG,                Zero});
        OpMapLO.insert({nf, OM});
      }

    // ===============================================================
    // NLO Matching conditions (Null).
    std::map<int, Operator> MatchNLO;
    for (int i = MatchingBasisQCD::M0; i <= MatchingBasisQCD::M7; i++)
      MatchNLO.insert({i, Zero});

    // ===============================================================
    // NLO splitting function operators.
    std::map<int, std::map<int, Operator>> OpMapNLO;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O1nsp{g, P1transnsp{nf}, IntEps};
        const Operator O1nsm{g, P1transnsm{nf}, IntEps};
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O1nsp});
        OM.insert({EvolutionBasisQCD::PNSM, O1nsm});
        OM.insert({EvolutionBasisQCD::PNSV, O1nsm});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O1nsp});
        OM.insert({EvolutionBasisQCD::PQG,                Zero});
        OM.insert({EvolutionBasisQCD::PGQ,                Zero});
        OM.insert({EvolutionBasisQCD::PGG,                Zero});
        OpMapNLO.insert({nf, OM});
      }

    // Define object of the structure containing the DglapObjects.
    std::map<int, DglapObjects> DglapObj;

    // Allocate convolution maps for evolution and matching, and set
    // of operators.
    for (int nf = nfi; nf <= nff; nf++)
      {
        DglapObjects obj;
        obj.Threshold = Thresholds[nf-1];
        if (OpEvol)
          {
            obj.SplittingFunctions.insert({0, Set<Operator>{EvolutionOperatorBasisQCD{nf}, OpMapLO.at(nf)}});
            obj.SplittingFunctions.insert({1, Set<Operator>{EvolutionOperatorBasisQCD{nf}, OpMapNLO.at(nf)}});
            obj.MatchingConditions.insert({0, Set<Operator>{MatchingOperatorBasisQCD{nf},  MatchLO}});
            obj.MatchingConditions.insert({1, Set<Operator>{MatchingOperatorBasisQCD{nf},  MatchNLO}});
          }
        else
          {
            obj.SplittingFunctions.insert({0, Set<Operator>{EvolutionBasisQCD{nf}, OpMapLO.at(nf)}});
            obj.SplittingFunctions.insert({1, Set<Operator>{EvolutionBasisQCD{nf}, OpMapNLO.at(nf)}});
            obj.MatchingConditions.insert({0, Set<Operator>{MatchingBasisQCD{nf},  MatchLO}});
            obj.MatchingConditions.insert({1, Set<Operator>{MatchingBasisQCD{nf},  MatchNLO}});
          }
        DglapObj.insert({nf,obj});
      }
    t.stop();

    return DglapObj;
  }

  //_____________________________________________________________________________
  std::map<int, DglapObjects> InitializeDglapObjectsQCDtrans(Grid                const& g,
                                                             std::vector<double> const& Thresholds,
                                                             bool                const& OpEvol,
                                                             double              const& IntEps)
  {
    return InitializeDglapObjectsQCDtrans(g, Thresholds, Thresholds, OpEvol, IntEps);
  }

  //_____________________________________________________________________________
  std::map<int, DglapObjects> InitializeDglapObjectsQCDTtrans(Grid                const& g,
                                                              std::vector<double> const& Masses,
                                                              std::vector<double> const& Thresholds,
                                                              bool                const& OpEvol,
                                                              double              const& IntEps)
  {
    report("Initializing DglapObjects for time-like QCD transversely polarised evolution... ");
    Timer t;

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the thresholds
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
        nfi++;

    // Compute logs of muth2 / m2 needed for the the matching
    // conditions (not possible for transversities yet).
    std::vector<double> LogKth;
    for (int im = 0; im < (int) Thresholds.size(); im++)
      if (Thresholds[im] < eps12 || Masses[im] < eps12)
        LogKth.push_back(0);
      else
        LogKth.push_back(2 * log( Thresholds[im] / Masses[im] ));

    // Allocate needed operators (matching conditions and splitting
    // functions). By now the code is fast enough to precompute
    // everything at all available perturbative orders and the current
    // perturbative order is accounted for only when the actual
    // splitting functions and matching conditions (lambda) functions
    // are defined.
    // ===============================================================
    // LO Matching conditions.
    std::map<int, Operator> MatchLO;
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};
    MatchLO.insert({MatchingBasisQCD::M0, Id});
    MatchLO.insert({MatchingBasisQCD::M1, Zero});
    MatchLO.insert({MatchingBasisQCD::M2, Zero});
    MatchLO.insert({MatchingBasisQCD::M3, Zero});
    MatchLO.insert({MatchingBasisQCD::M4, Zero});
    MatchLO.insert({MatchingBasisQCD::M5, Zero});
    MatchLO.insert({MatchingBasisQCD::M6, Zero});
    MatchLO.insert({MatchingBasisQCD::M7, Zero});

    // ===============================================================
    // LO splitting function operators.
    std::map<int, std::map<int, Operator>> OpMapLO;
    const Operator O0ns{g, P0Ttransns{}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O0ns});
        OM.insert({EvolutionBasisQCD::PNSM, O0ns});
        OM.insert({EvolutionBasisQCD::PNSV, O0ns});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O0ns});
        OM.insert({EvolutionBasisQCD::PQG,                Zero});
        OM.insert({EvolutionBasisQCD::PGQ,                Zero});
        OM.insert({EvolutionBasisQCD::PGG,                Zero});
        OpMapLO.insert({nf, OM});
      }

    // ===============================================================
    // NLO Matching conditions (Null).
    std::map<int, Operator> MatchNLO;
    for (int i = MatchingBasisQCD::M0; i <= MatchingBasisQCD::M7; i++)
      MatchNLO.insert({i, Zero});

    // ===============================================================
    // NLO splitting function operators.
    std::map<int, std::map<int, Operator>> OpMapNLO;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O1nsp{g, P1Ttransnsp{nf}, IntEps};
        const Operator O1nsm{g, P1Ttransnsm{nf}, IntEps};
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O1nsp});
        OM.insert({EvolutionBasisQCD::PNSM, O1nsm});
        OM.insert({EvolutionBasisQCD::PNSV, O1nsm});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O1nsp});
        OM.insert({EvolutionBasisQCD::PQG,                Zero});
        OM.insert({EvolutionBasisQCD::PGQ,                Zero});
        OM.insert({EvolutionBasisQCD::PGG,                Zero});
        OpMapNLO.insert({nf, OM});
      }

    // Define object of the structure containing the DglapObjects.
    std::map<int, DglapObjects> DglapObj;

    // Allocate convolution maps for evolution and matching, and set
    // of operators.
    for (int nf = nfi; nf <= nff; nf++)
      {
        DglapObjects obj;
        obj.Threshold = Thresholds[nf-1];
        if (OpEvol)
          {
            obj.SplittingFunctions.insert({0, Set<Operator>{EvolutionOperatorBasisQCD{nf}, OpMapLO.at(nf)}});
            obj.SplittingFunctions.insert({1, Set<Operator>{EvolutionOperatorBasisQCD{nf}, OpMapNLO.at(nf)}});
            obj.MatchingConditions.insert({0, Set<Operator>{MatchingOperatorBasisQCD{nf},  MatchLO}});
            obj.MatchingConditions.insert({1, Set<Operator>{MatchingOperatorBasisQCD{nf},  MatchNLO}});
          }
        else
          {
            obj.SplittingFunctions.insert({0, Set<Operator>{EvolutionBasisQCD{nf}, OpMapLO.at(nf)}});
            obj.SplittingFunctions.insert({1, Set<Operator>{EvolutionBasisQCD{nf}, OpMapNLO.at(nf)}});
            obj.MatchingConditions.insert({0, Set<Operator>{MatchingBasisQCD{nf},  MatchLO}});
            obj.MatchingConditions.insert({1, Set<Operator>{MatchingBasisQCD{nf},  MatchNLO}});
          }
        DglapObj.insert({nf,obj});
      }
    t.stop();

    return DglapObj;
  }

  //_____________________________________________________________________________
  std::map<int, DglapObjects> InitializeDglapObjectsQCDTtrans(Grid                const& g,
                                                              std::vector<double> const& Thresholds,
                                                              bool                const& OpEvol,
                                                              double              const& IntEps)
  {
    return InitializeDglapObjectsQCDTtrans(g, Thresholds, Thresholds, OpEvol, IntEps);
  }

  //_____________________________________________________________________________
  std::function<Set<Operator>(int const&, double const&)> SplittingFunctions(std::map<int, DglapObjects>          const& DglapObj,
                                                                             int                                  const& PerturbativeOrder,
                                                                             std::function<double(double const&)> const& Alphas)
  {
    if (PerturbativeOrder == 0)
      return [=] (int const& nf, double const& mu) -> Set<Operator>
      {
        const double cp = Alphas(mu) / FourPi;
        return cp * DglapObj.at(nf).SplittingFunctions.at(0);
      };
    else if (PerturbativeOrder == 1)
      return [=] (int const& nf, double const& mu) -> Set<Operator>
      {
        const double cp = Alphas(mu) / FourPi;
        const auto sf = DglapObj.at(nf).SplittingFunctions;
        return cp * ( sf.at(0) + cp * sf.at(1) );
      };
    else if (PerturbativeOrder == 2)
      return [=] (int const& nf, double const& mu) -> Set<Operator>
      {
        const double cp = Alphas(mu) / FourPi;
        const auto sf = DglapObj.at(nf).SplittingFunctions;
        return cp * ( sf.at(0) + cp * ( sf.at(1) + cp * sf.at(2) ) );
      };
    else if (PerturbativeOrder == 3)
      return [=] (int const& nf, double const& mu) -> Set<Operator>
      {
        const double cp = Alphas(mu) / FourPi;
        const auto sf = DglapObj.at(nf).SplittingFunctions;
        return cp * ( sf.at(0) + cp * ( sf.at(1) + cp * ( sf.at(2) + cp * sf.at(3) ) ) );
      };
    else
      throw std::runtime_error(error("SplittingFunctions","Perturbative order not allowed."));
  }

  //_____________________________________________________________________________
  std::function<Set<Operator>(bool const&, int const&)> MatchingConditions(std::map<int, DglapObjects>          const& DglapObj,
                                                                           int                                  const& PerturbativeOrder,
                                                                           std::function<double(double const&)> const& Alphas)
  {
    // Compute coupling above and below the thresholds.
    std::map<int, double> asThUp;
    std::map<int, double> asThDown;
    for (auto obj = std::next(DglapObj.begin()); obj != DglapObj.end(); ++obj)
      {
        const int    nf  = obj->first;
        const double thr = obj->second.Threshold;
        asThDown.insert({nf, Alphas(thr * ( 1 - eps8 ) ) / FourPi});
        asThUp.insert({nf, Alphas(thr * ( 1 + eps8 ) ) / FourPi});
      }

    if (PerturbativeOrder == 0)
      return [=] (bool const&, int const& nf) -> Set<Operator>
      {
        return DglapObj.at(nf).MatchingConditions.at(0);
      };
    else if (PerturbativeOrder == 1)
      return [=] (bool const& Up, int const& nf) -> Set<Operator>
      {
        const double cp = (Up ? asThUp.at(nf+1) : asThDown.at(nf));
        const auto mc = DglapObj.at(nf).MatchingConditions;
        return mc.at(0) + (Up ? 1 : -1) * cp * mc.at(1);
      };
    else if (PerturbativeOrder == 2)
      return [=] (bool const& Up, int const& nf) -> Set<Operator>
      {
        const double cp = (Up ? asThUp.at(nf+1) : asThDown.at(nf));
        const auto mc = DglapObj.at(nf).MatchingConditions;
        return mc.at(0) + (Up ? 1 : -1) * cp * ( mc.at(1) + cp * ( mc.at(2) - (Up ? 0 : 1) * mc.at(-2) ) );
      };
    else if (PerturbativeOrder == 3)
      return [=] (bool const& Up, int const& nf) -> Set<Operator>
      {
        const double cp = (Up ? asThUp.at(nf+1) : asThDown.at(nf));
        const auto mc = DglapObj.at(nf).MatchingConditions;
        return mc.at(0) + (Up ? 1 : -1) * cp * ( mc.at(1) + cp * ( mc.at(2) - (Up ? 0 : 1) * mc.at(-2) ) );
      };
    else
      throw std::runtime_error(error("MatchingConditions","Perturbative order not allowed."));
  }

  //_____________________________________________________________________________
  std::unique_ptr<Dglap<Distribution>> BuildDglap(std::map<int, DglapObjects>                                        const& DglapObj,
                                                  std::function<std::map<int, double>(double const&, double const&)> const& InDistFunc,
                                                  double                                                             const& MuRef,
                                                  int                                                                const& PerturbativeOrder,
                                                  std::function<double(double const&)>                               const& Alphas,
                                                  int                                                                const& nsteps)
  {
    // Collect thresholds.
    std::vector<double> Thresholds;
    for (auto const& obj : DglapObj)
      {
        const int    nf  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) Thresholds.size() < nf)
          Thresholds.resize(nf);
        Thresholds[nf-1] = thr;
      }

    // Create set of initial distributions.
    const Set<Distribution> InPDFs{DglapObj.at(NF(MuRef, Thresholds)).SplittingFunctions.at(0).GetMap(),
                                   DistributionMap(DglapObj.begin()->second.SplittingFunctions.at(0).at(0).GetGrid(), InDistFunc, MuRef)};

    // Initialize DGLAP evolution.
    return std::unique_ptr<Dglap<Distribution>>(new Dglap<Distribution> {SplittingFunctions(DglapObj, PerturbativeOrder, Alphas),
                                                                         MatchingConditions(DglapObj, PerturbativeOrder, Alphas), InPDFs, MuRef, Thresholds, nsteps
                                                                        });
  }

  //_____________________________________________________________________________
  std::unique_ptr<Dglap<Operator>> BuildDglap(std::map<int, DglapObjects>          const& DglapObj,
                                              double                               const& MuRef,
                                              int                                  const& PerturbativeOrder,
                                              std::function<double(double const&)> const& Alphas,
                                              int                                  const& nsteps)
  {
    // Collect thresholds.
    std::vector<double> Thresholds;
    for (auto const& obj : DglapObj)
      {
        const int    nf  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) Thresholds.size() < nf)
          Thresholds.resize(nf);
        Thresholds[nf-1] = thr;
      }

    // Allocate Identity and Zero operators.
    const Operator One{DglapObj.begin()->second.SplittingFunctions.at(0).at(0).GetGrid(), Identity{}};
    const Operator Zero{DglapObj.begin()->second.SplittingFunctions.at(0).at(0).GetGrid(), Null{}};

    // Create set of initial operators that represent the unity set of
    // operators.
    std::map<int, Operator> MapUnity;
    for (int i = 0; i < 55; i++)
      if (i == 0  || i == 8  || i == 14 || i == 17 || i == 22 || i == 26 || i == 30 ||
          i == 35 || i == 38 || i == 44 || i == 46 || i == 53 || i == 54)
        MapUnity.insert({i, One});
      else
        MapUnity.insert({i, Zero});
    Set<Operator> Unity{EvolutionOperatorBasisQCD{NF(MuRef, Thresholds)}, MapUnity};

    // Initialize DGLAP evolution.
    return std::unique_ptr<Dglap<Operator>>(new Dglap<Operator> {SplittingFunctions(DglapObj, PerturbativeOrder, Alphas),
                                                                 MatchingConditions(DglapObj, PerturbativeOrder, Alphas), Unity, MuRef, Thresholds, nsteps
                                                                });
  }
}
