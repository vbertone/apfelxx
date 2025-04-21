//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/dglapbuilderqcdqed.h"
#include "apfel/timer.h"
#include "apfel/messages.h"
#include "apfel/tools.h"
#include "apfel/matchingbasisqcdqed.h"
#include "apfel/evolutionbasisqcdqed.h"
#include "apfel/splittingfunctionsunp_sl.h"
#include "apfel/splittingfunctionsunp_sl_qed.h"
#include "apfel/matchingconditions_sl.h"

namespace apfel
{
  //_____________________________________________________________________________
  std::map<int, DglapObjectsQCDQED> InitializeDglapObjectsQCDQED(Grid                const& g,
                                                                 std::vector<double> const& QuarkThresholds,
                                                                 std::vector<double> const& LeptonThresholds,
                                                                 bool                const& OpEvol,
                                                                 double              const& IntEps,
                                                                 bool                const& n3lo,
                                                                 std::vector<int>    const& IMod)
  {
    if (!n3lo)
      warning("InitializeDglapObjectsQCDQED", "N3LO corrections will not be initialised.");
    report("Initializing DglapObjectsQCDQED for space-like QCDxQED unpolarised evolution... ");
    Timer t;

    // Assembly total set of threshold by concatenating and sorting
    // quark and lepton thresholds.
    const std::vector<double> Thresholds = ConcatenateAndSortVectors(QuarkThresholds, LeptonThresholds);

    // Number of down and up quarks depending on nf
    const std::vector<int> NDW = {0, 0, 1, 2, 2, 3, 3};
    const std::vector<int> NUP = {0, 1, 1, 1, 2, 2, 3};

    // Determine nd, nd, and nl given the total number of active
    // partons.
    std::vector<std::vector<int>> NDUL;
    for (int nt = 0; nt <= (int) Thresholds.size(); nt++)
      {
        const int nl = NF((nt > 0 ? Thresholds[nt - 1] : 0) + eps8, LeptonThresholds);
        const int nf = NF((nt > 0 ? Thresholds[nt - 1] : 0) + eps8, QuarkThresholds);
        const int nd = NDW[nf];
        const int nu = NUP[nf];
        NDUL.push_back({nd, nu, nl});
      }

    // Determine initial and final number of active flavours according
    // to the vector of thresholds.
    int nti = 0;
    int ntf = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
        nti++;

    // Determine parton species give the number of active partons
    std::map<int, PartonSpecies> Species;
    for (int nt = nti; nt <= ntf; nt++)
      if (NDUL[nt][0] - NDUL[std::min(nt + 1, ntf)][0] != 0)
        Species.insert({nt, PartonSpecies::DOWNTYPEQUARK});
      else if (NDUL[nt][1] - NDUL[std::min(nt + 1, ntf)][1] != 0)
        Species.insert({nt, PartonSpecies::UPTYPEQUARK});
      else if (NDUL[nt][2] - NDUL[std::min(nt + 1, ntf)][2] != 0)
        Species.insert({nt, PartonSpecies::CHARGEDLEPTON});
      else
        Species.insert({nt, PartonSpecies::UNKNOWN});

    // Compute logs of muth2 / m2 needed for the matching
    // conditions. Push a zero at last to extend the vector but that
    // entry will never be effectively used in the evolution. Set to
    // zero for now.
    std::vector<double> LogKth;
    for (int im = 0; im < (int) Thresholds.size(); im++)
      LogKth.push_back(0);
    LogKth.push_back(0);

    // Allocate needed operators (matching conditions and splitting
    // functions). By now the code is fast enough to precompute
    // everything at all available perturbative orders and the current
    // perturbative order is accounted for only when the actual
    // splitting functions and matching conditions (lambda) functions
    // are defined.
    // ===============================================================
    // LO matching conditions
    std::map<int, Operator> Match00;
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};
    Match00.insert({MatchingBasisQCDQED::ONE,    Id});
    Match00.insert({MatchingBasisQCDQED::KQg,    Zero});
    Match00.insert({MatchingBasisQCDQED::KXgm,   Zero});
    Match00.insert({MatchingBasisQCDQED::KQqp,   Zero});
    Match00.insert({MatchingBasisQCDQED::KXX,    Zero});
    Match00.insert({MatchingBasisQCDQED::Kqg,    Zero});
    Match00.insert({MatchingBasisQCDQED::KNSq,   Zero});
    Match00.insert({MatchingBasisQCDQED::Kqqp,   Zero});
    Match00.insert({MatchingBasisQCDQED::Kgg,    Zero});
    Match00.insert({MatchingBasisQCDQED::Kgq,    Zero});
    Match00.insert({MatchingBasisQCDQED::KgQ,    Zero});
    Match00.insert({MatchingBasisQCDQED::KXgmgm, Zero});
    Match00.insert({MatchingBasisQCDQED::KgmX,   Zero});

    // ===============================================================
    // O(as) splitting function operators
    std::map<int, std::map<int, Operator>> OpMap10;
    const Operator O10ns{g, P0ns{},  IntEps};
    const Operator O10gq{g, P0gq{},  IntEps};
    const Operator O10qg{g, P0qg{1}, IntEps};
    for (int nt = nti; nt <= ntf; nt++)
      {
        // Determine number of active quarks
        const int nf = NDUL[nt][0] + NDUL[nt][1];
        const Operator O10gg{g, P0gg{nf}, IntEps};
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCDQED::PPDD, O10ns});
        OM.insert({EvolutionBasisQCDQED::PPUU, O10ns});
        OM.insert({EvolutionBasisQCDQED::PMDD, O10ns});
        OM.insert({EvolutionBasisQCDQED::PMUU, O10ns});
        OM.insert({EvolutionBasisQCDQED::PDg,  O10qg});
        OM.insert({EvolutionBasisQCDQED::PUg,  O10qg});
        OM.insert({EvolutionBasisQCDQED::PgD,  O10gq});
        OM.insert({EvolutionBasisQCDQED::PgU,  O10gq});
        OM.insert({EvolutionBasisQCDQED::Pgg,  O10gg});
        // Insert Zero in the remaining slots
        for (int i = EvolutionBasisQCDQED::PPDD; i <= EvolutionBasisQCDQED::Pgmgm; i++)
          OM.insert({i, Zero});
        OpMap10.insert({nt, OM});
      }

    // ===============================================================
    // O(a) splitting function operators
    std::map<int, std::map<int, Operator>> OpMap01;
    const Operator O01ns  {g, P01qedns{},   IntEps};
    const Operator O01gmq {g, P01qedgmq{},  IntEps};
    const Operator O01qgm {g, P01qedqgm{},  IntEps};
    const Operator O01gmgm{g, P01qedgmgm{}, IntEps};
    for (int nt = nti; nt <= ntf; nt++)
      {
        // Determine number of active quarks and leptons
        const int nf = NDUL[nt][0] + NDUL[nt][1];
        const int nl = NDUL[nt][2];
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCDQED::PPDD, ed2 * O01ns});
        OM.insert({EvolutionBasisQCDQED::PPUU, eu2 * O01ns});
        OM.insert({EvolutionBasisQCDQED::PPLL,       O01ns});
        OM.insert({EvolutionBasisQCDQED::PMDD, ed2 * O01ns});
        OM.insert({EvolutionBasisQCDQED::PMUU, eu2 * O01ns});
        OM.insert({EvolutionBasisQCDQED::PMLL,       O01ns});
        OM.insert({EvolutionBasisQCDQED::PDgm, ( NC * ed2 ) * O01qgm});
        OM.insert({EvolutionBasisQCDQED::PUgm, ( NC * eu2 ) * O01qgm});
        OM.insert({EvolutionBasisQCDQED::PLgm,                O01qgm});
        OM.insert({EvolutionBasisQCDQED::PgmD, ed2 * O01gmq});
        OM.insert({EvolutionBasisQCDQED::PgmU, eu2 * O01gmq});
        OM.insert({EvolutionBasisQCDQED::PgmL,       O01gmq});
        OM.insert({EvolutionBasisQCDQED::Pgmgm, ( NC * SumCh2[nf] + nl ) * O01gmgm});
        // Insert Zero in the remaining slots
        for (int i = EvolutionBasisQCDQED::PPDD; i <= EvolutionBasisQCDQED::Pgmgm; i++)
          OM.insert({i, Zero});
        OpMap01.insert({nt, OM});
      }

    // ===============================================================
    // O(as) matching conditions
    std::map<int, std::map<int, Operator>> Match10;
    const Operator AS1HgL {g, AS1Hg_L{},  IntEps};
    const Operator AS1ggHL{g, AS1ggH_L{}, IntEps};
    const Operator AS1gH0 {g, AS1gH_0{},  IntEps};
    const Operator AS1gHL {g, AS1gH_L{},  IntEps};
    const Operator AS1HH0 {g, AS1HH_0{},  IntEps};
    const Operator AS1HHL {g, AS1HH_L{},  IntEps};
    for (int nt = nti; nt <= ntf; nt++)
      {
        const Operator AS1Hg  =          LogKth[nt] * AS1HgL;
        const Operator AS1ggH =          LogKth[nt] * AS1ggHL;
        const Operator AS1gH  = AS1gH0 + LogKth[nt] * AS1gHL;
        const Operator AS1HH  = AS1HH0 + LogKth[nt] * AS1HHL;
        std::map<int, Operator> OM;
        OM.insert({MatchingBasisQCDQED::ONE,    Zero});
        OM.insert({MatchingBasisQCDQED::KXgm,   Zero});
        OM.insert({MatchingBasisQCDQED::KQqp,   Zero});
        OM.insert({MatchingBasisQCDQED::Kqg,    Zero});
        OM.insert({MatchingBasisQCDQED::KNSq,   Zero});
        OM.insert({MatchingBasisQCDQED::Kqqp,   Zero});
        OM.insert({MatchingBasisQCDQED::Kgq,    Zero});
        OM.insert({MatchingBasisQCDQED::KXgmgm, Zero});
        OM.insert({MatchingBasisQCDQED::KgmX,   Zero});
        if (Species.at(nt) == PartonSpecies::DOWNTYPEQUARK || Species.at(nt) == PartonSpecies::UPTYPEQUARK)
          {
            OM.insert({MatchingBasisQCDQED::KQg, AS1Hg});
            OM.insert({MatchingBasisQCDQED::KXX, AS1HH});
            OM.insert({MatchingBasisQCDQED::Kgg, AS1ggH});
            OM.insert({MatchingBasisQCDQED::KgQ, AS1gH});
          }
        else
          {
            OM.insert({MatchingBasisQCDQED::KQg, Zero});
            OM.insert({MatchingBasisQCDQED::KXX, Zero});
            OM.insert({MatchingBasisQCDQED::Kgg, Zero});
            OM.insert({MatchingBasisQCDQED::KgQ, Zero});
          }
        Match10.insert({nt, OM});
      }

    // ===============================================================
    // O(a) matching conditions
    std::map<int, std::map<int, Operator>> Match01;
    const Operator AS1qedHgmL  {g, AS1qedHgm_L{},   IntEps};
    const Operator AS1qedgmgmHL{g, AS1qedgmgmH_L{}, IntEps};
    const Operator AS1qedgmH0  {g, AS1qedgmH_0{},   IntEps};
    const Operator AS1qedgmHL  {g, AS1qedgmH_L{},   IntEps};
    const Operator AS1qedHH0   {g, AS1qedHH_0{},    IntEps};
    const Operator AS1qedHHL   {g, AS1qedHH_L{},    IntEps};
    for (int nt = nti; nt <= ntf; nt++)
      {
        const Operator AS1qedHgm   =              LogKth[nt] * AS1qedHgmL;
        const Operator AS1qedgmgmH =              LogKth[nt] * AS1qedgmgmHL;
        const Operator AS1qedgmH   = AS1qedgmH0 + LogKth[nt] * AS1qedgmHL;
        const Operator AS1qedHH    = AS1qedHH0  + LogKth[nt] * AS1qedHHL;
        std::map<int, Operator> OM;
        OM.insert({MatchingBasisQCDQED::ONE,  Zero});
        OM.insert({MatchingBasisQCDQED::KQg,  Zero});
        OM.insert({MatchingBasisQCDQED::KQqp, Zero});
        OM.insert({MatchingBasisQCDQED::Kqg,  Zero});
        OM.insert({MatchingBasisQCDQED::KNSq, Zero});
        OM.insert({MatchingBasisQCDQED::Kqqp, Zero});
        OM.insert({MatchingBasisQCDQED::Kgq,  Zero});
        OM.insert({MatchingBasisQCDQED::Kgg,  Zero});
        OM.insert({MatchingBasisQCDQED::KgQ,  Zero});
        if (Species.at(nt) == PartonSpecies::CHARGEDLEPTON)
          {
            OM.insert({MatchingBasisQCDQED::KXgm,   AS1qedHgm});
            OM.insert({MatchingBasisQCDQED::KXX,    AS1qedHH});
            OM.insert({MatchingBasisQCDQED::KXgmgm, AS1qedgmgmH});
            OM.insert({MatchingBasisQCDQED::KgmX,   AS1qedgmH});
          }
        else
          {
            OM.insert({MatchingBasisQCDQED::KXgm,   Zero});
            OM.insert({MatchingBasisQCDQED::KXX,    Zero});
            OM.insert({MatchingBasisQCDQED::KXgmgm, Zero});
            OM.insert({MatchingBasisQCDQED::KgmX,   Zero});
          }
        Match01.insert({nt, OM});
      }

    // ===============================================================
    // O(as^2) splitting function operators
    std::map<int, std::map<int, Operator>> OpMap20;
    const Operator O20qg{g, P1qg{1}, IntEps};
    const Operator O20ps{g, P1ps{1}, IntEps};
    for (int nt = nti; nt <= ntf; nt++)
      {
        // Determine number of active quarks
        const int nf = NDUL[nt][0] + NDUL[nt][1];
        const Operator O20nsp{g, P1nsp{nf}, IntEps};
        const Operator O20nsm{g, P1nsm{nf}, IntEps};
        const Operator O20gq {g, P1gq{nf},  IntEps};
        const Operator O20gg {g, P1gg{nf},  IntEps};
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCDQED::PPDD,  O20nsp});
        OM.insert({EvolutionBasisQCDQED::PPUU,  O20nsp});
        OM.insert({EvolutionBasisQCDQED::PMDD,  O20nsm});
        OM.insert({EvolutionBasisQCDQED::PMUU,  O20nsm});
        OM.insert({EvolutionBasisQCDQED::PPSDD, O20ps});
        OM.insert({EvolutionBasisQCDQED::PPSDU, O20ps});
        OM.insert({EvolutionBasisQCDQED::PPSUD, O20ps});
        OM.insert({EvolutionBasisQCDQED::PPSUU, O20ps});
        OM.insert({EvolutionBasisQCDQED::PDg,   O20qg});
        OM.insert({EvolutionBasisQCDQED::PUg,   O20qg});
        OM.insert({EvolutionBasisQCDQED::PgD,   O20gq});
        OM.insert({EvolutionBasisQCDQED::PgU,   O20gq});
        OM.insert({EvolutionBasisQCDQED::Pgg,   O20gg});
        // Insert Zero in the remaining slots
        for (int i = EvolutionBasisQCDQED::PPDD; i <= EvolutionBasisQCDQED::Pgmgm; i++)
          OM.insert({i, Zero});
        OpMap20.insert({nt, OM});
      }

    // ===============================================================
    // O(as*a) splitting function operators
    std::map<int, std::map<int, Operator>> OpMap11;

    const Operator O11nsp {g, P11qednsp{},  IntEps};
    const Operator O11nsm {g, P11qednsm{},  IntEps};
    const Operator O11qg  {g, P11qedqg{},   IntEps};
    const Operator O11qgm {g, P11qedqgm{},  IntEps};
    const Operator O11gq  {g, P11qedgq{},   IntEps};
    const Operator O11gmq {g, P11qedgmq{},  IntEps};
    const Operator O11gg  {g, P11qedgg{},   IntEps};
    const Operator O11ggm {g, P11qedggm{},  IntEps};
    const Operator O11gmg {g, P11qedgmg{},  IntEps};
    const Operator O11gmgm{g, P11qedgmgm{}, IntEps};
    for (int nt = nti; nt <= ntf; nt++)
      {
        // Determine number of active quarks
        const int nd = NDUL[nt][0];
        const int nu = NDUL[nt][1];
        const double bbf = nd * ed2 + nu * eu2;
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCDQED::PPDD,  ed2 * O11nsp});
        OM.insert({EvolutionBasisQCDQED::PPUU,  eu2 * O11nsp});
        OM.insert({EvolutionBasisQCDQED::PMDD,  ed2 * O11nsm});
        OM.insert({EvolutionBasisQCDQED::PMUU,  eu2 * O11nsm});
        OM.insert({EvolutionBasisQCDQED::PDg,   ed2 * O11qg});
        OM.insert({EvolutionBasisQCDQED::PUg,   eu2 * O11qg});
        OM.insert({EvolutionBasisQCDQED::PDgm,  ( NC * ed2 ) * O11qgm});
        OM.insert({EvolutionBasisQCDQED::PUgm,  ( NC * eu2 ) * O11qgm});
        OM.insert({EvolutionBasisQCDQED::PgD,   ed2 * O11gq});
        OM.insert({EvolutionBasisQCDQED::PgU,   eu2 * O11gq});
        OM.insert({EvolutionBasisQCDQED::PgmD,  ed2 * O11gmq});
        OM.insert({EvolutionBasisQCDQED::PgmU,  eu2 * O11gmq});
        OM.insert({EvolutionBasisQCDQED::Pgg,          bbf   * O11gg});
        OM.insert({EvolutionBasisQCDQED::Pggm,  ( NC * bbf ) * O11ggm});
        OM.insert({EvolutionBasisQCDQED::Pgmg,         bbf   * O11gmg});
        OM.insert({EvolutionBasisQCDQED::Pgmgm, ( NC * bbf ) * O11gmgm});
        // Insert Zero in the remaining slots
        for (int i = EvolutionBasisQCDQED::PPDD; i <= EvolutionBasisQCDQED::Pgmgm; i++)
          OM.insert({i, Zero});
        OpMap11.insert({nt, OM});
      }

    // ===============================================================
    // O(a^2) splitting function operators
    std::map<int, std::map<int, Operator>> OpMap02;
    const Operator O02qgm {g, P02qedqgm{},  IntEps};
    const Operator O02ps  {g, P02qedps{},   IntEps};
    const Operator O02gmgm{g, P02qedgmgm{}, IntEps};
    for (int nt = nti; nt <= ntf; nt++)
      {
        // Determine number of active quarks and leptons
        const int nf = NDUL[nt][0] + NDUL[nt][1];
        const int nl = NDUL[nt][2];
        const double eSigma2 = ( NC * SumCh2[nf] + nl );
        const double ed4     = ed2 * ed2;
        const double eu4     = eu2 * eu2;
        const Operator O02nspD{g, P02qednsm{eSigma2 / ed2}, IntEps};
        const Operator O02nsmD{g, P02qednsm{eSigma2 / ed2}, IntEps};
        const Operator O02gmqD{g, P02qedgmq{eSigma2 / ed2}, IntEps};
        const Operator O02nspU{g, P02qednsm{eSigma2 / eu2}, IntEps};
        const Operator O02nsmU{g, P02qednsm{eSigma2 / eu2}, IntEps};
        const Operator O02gmqU{g, P02qedgmq{eSigma2 / eu2}, IntEps};
        const Operator O02nspL{g, P02qednsm{eSigma2},       IntEps};
        const Operator O02nsmL{g, P02qednsm{eSigma2},       IntEps};
        const Operator O02gmqL{g, P02qedgmq{eSigma2},       IntEps};
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCDQED::PPDD, ed4 * O02nspD});
        OM.insert({EvolutionBasisQCDQED::PPUU, eu4 * O02nspU});
        OM.insert({EvolutionBasisQCDQED::PPLL,       O02nspL});
        OM.insert({EvolutionBasisQCDQED::PMDD, ed4 * O02nsmD});
        OM.insert({EvolutionBasisQCDQED::PMUU, eu4 * O02nsmU});
        OM.insert({EvolutionBasisQCDQED::PMLL,       O02nsmL});
        OM.insert({EvolutionBasisQCDQED::PPSDD, ( NC * ed2 * ed2 ) * O02ps});
        OM.insert({EvolutionBasisQCDQED::PPSDU, ( NC * ed2 * eu2 ) * O02ps});
        OM.insert({EvolutionBasisQCDQED::PPSDL,   NC * ed2         * O02ps});
        OM.insert({EvolutionBasisQCDQED::PPSUD, ( NC * eu2 * ed2 ) * O02ps});
        OM.insert({EvolutionBasisQCDQED::PPSUU, ( NC * eu2 * eu2 ) * O02ps});
        OM.insert({EvolutionBasisQCDQED::PPSUL,   NC * eu2         * O02ps});
        OM.insert({EvolutionBasisQCDQED::PPSLD,              ed2   * O02ps});
        OM.insert({EvolutionBasisQCDQED::PPSLU,              eu2   * O02ps});
        OM.insert({EvolutionBasisQCDQED::PPSLL,                      O02ps});
        OM.insert({EvolutionBasisQCDQED::PDgm, ( NC * ed4 ) * O02qgm});
        OM.insert({EvolutionBasisQCDQED::PUgm, ( NC * eu4 ) * O02qgm});
        OM.insert({EvolutionBasisQCDQED::PLgm,                O02qgm});
        OM.insert({EvolutionBasisQCDQED::PgmD, ed4 * O02gmqU});
        OM.insert({EvolutionBasisQCDQED::PgmU, eu4 * O02gmqD});
        OM.insert({EvolutionBasisQCDQED::PgmL,       O02gmqL});
        OM.insert({EvolutionBasisQCDQED::Pgmgm, ( NC * SumCh4[nf] + nl ) * O02gmgm});
        // Insert Zero in the remaining slots
        for (int i = EvolutionBasisQCDQED::PPDD; i <= EvolutionBasisQCDQED::Pgmgm; i++)
          OM.insert({i, Zero});
        OpMap02.insert({nt, OM});
      }

    // Define object of the structure containing the DglapObjects
    std::map<int, DglapObjectsQCDQED> DglapObj;

    // Allocate convolution maps for evolution and matching, and set
    // of operators.
    for (int nt = nti; nt <= ntf; nt++)
      {
        DglapObjectsQCDQED obj;
        const int nd = NDUL[nt][0];
        const int nu = NDUL[nt][1];
        const int nl = NDUL[nt][2];
        obj.Threshold = (nt > 0 ? Thresholds[nt - 1] : 0);
        obj.Species = Species.at(nt);
        if (OpEvol)
          {
            std::map<int, Operator> MapUnity;
            for (auto const& coord : GkjQCDQED)
              MapUnity.insert({coord.second, (coord.first.first == coord.first.second ? Id : Zero)});
            obj.UnitySet = Set<Operator> {EvolutionOperatorBasisQCDQED{nd, nu, nl}, MapUnity};
            obj.SplittingFunctions.insert({{ 1, 0}, Set<Operator>{EvolutionOperatorBasisQCDQED{nd, nu, nl}, OpMap10.at(nt)}});
            obj.SplittingFunctions.insert({{ 0, 1}, Set<Operator>{EvolutionOperatorBasisQCDQED{nd, nu, nl}, OpMap01.at(nt)}});
            obj.SplittingFunctions.insert({{ 2, 0}, Set<Operator>{EvolutionOperatorBasisQCDQED{nd, nu, nl}, OpMap20.at(nt)}});
            obj.SplittingFunctions.insert({{ 1, 1}, Set<Operator>{EvolutionOperatorBasisQCDQED{nd, nu, nl}, OpMap11.at(nt)}});
            obj.SplittingFunctions.insert({{ 0, 2}, Set<Operator>{EvolutionOperatorBasisQCDQED{nd, nu, nl}, OpMap02.at(nt)}});
            obj.MatchingConditions.insert({{ 0, 0}, Set<Operator>{MatchingOperatorBasisQCDQED{nd, nu, nl, obj.Species}, Match00}});
            obj.MatchingConditions.insert({{ 1, 0}, Set<Operator>{MatchingOperatorBasisQCDQED{nd, nu, nl, obj.Species}, Match10.at(nt)}});
            obj.MatchingConditions.insert({{ 0, 1}, Set<Operator>{MatchingOperatorBasisQCDQED{nd, nu, nl, obj.Species}, Match01.at(nt)}});
          }
        else
          {
            obj.SplittingFunctions.insert({{ 1, 0}, Set<Operator>{EvolutionBasisQCDQED{nd, nu, nl}, OpMap10.at(nt)}});
            obj.SplittingFunctions.insert({{ 0, 1}, Set<Operator>{EvolutionBasisQCDQED{nd, nu, nl}, OpMap01.at(nt)}});
            obj.SplittingFunctions.insert({{ 2, 0}, Set<Operator>{EvolutionBasisQCDQED{nd, nu, nl}, OpMap20.at(nt)}});
            obj.SplittingFunctions.insert({{ 1, 1}, Set<Operator>{EvolutionBasisQCDQED{nd, nu, nl}, OpMap11.at(nt)}});
            obj.SplittingFunctions.insert({{ 0, 2}, Set<Operator>{EvolutionBasisQCDQED{nd, nu, nl}, OpMap02.at(nt)}});
            obj.MatchingConditions.insert({{ 0, 0}, Set<Operator>{MatchingBasisQCDQED{nd, nu, nl, obj.Species}, Match00}});
            obj.MatchingConditions.insert({{ 1, 0}, Set<Operator>{MatchingBasisQCDQED{nd, nu, nl, obj.Species}, Match10.at(nt)}});
            obj.MatchingConditions.insert({{ 0, 1}, Set<Operator>{MatchingBasisQCDQED{nd, nu, nl, obj.Species}, Match01.at(nt)}});
          }
        DglapObj.insert({nt, obj});
      }
    t.stop();

    return DglapObj;
  }

  //_____________________________________________________________________________
  std::function<Set<Operator>(int const&, double const&)> SplittingFunctionsQCDQED(std::map<int, DglapObjectsQCDQED>    const& DglapObj,
                                                                                   int                                  const& PerturbativeOrder,
                                                                                   std::function<double(double const&)> const& Alphas,
                                                                                   std::function<double(double const&)> const& Alphaem)
  {
    if (PerturbativeOrder == 0)
      return [=] (int const& nt, double const& t) -> Set<Operator>
      {
        return ( Alphas(exp(t / 2)) / FourPi )  * DglapObj.at(nt).SplittingFunctions.at({1, 0})
        + ( Alphaem(exp(t / 2)) / FourPi ) * DglapObj.at(nt).SplittingFunctions.at({0, 1});
      };
    else if (PerturbativeOrder == 1)
      return [=] (int const& nt, double const& t) -> Set<Operator>
      {
        const double cp10 = Alphas(exp(t / 2)) / FourPi;
        const double cp01 = Alphaem(exp(t / 2)) / FourPi;
        const double cp20 = cp10 * cp10;
        const double cp11 = cp10 * cp01;
        const double cp02 = cp01 * cp01;
        const auto sf = DglapObj.at(nt).SplittingFunctions;
        return cp10 * sf.at({1, 0}) + cp01 * sf.at({0, 1})
        + cp20 * sf.at({2, 0}) + cp11 * sf.at({1, 1}) + cp02 * sf.at({0, 2});
      };
    else if (PerturbativeOrder == 2)
      return [=] (int const& nt, double const& t) -> Set<Operator>
      {
        const double cp10 = Alphas(exp(t / 2)) / FourPi;
        const double cp01 = Alphaem(exp(t / 2)) / FourPi;
        const double cp20 = cp10 * cp10;
        const double cp11 = cp10 * cp01;
        const double cp02 = cp01 * cp01;
        const double cp30 = cp10 * cp20;
        const auto sf = DglapObj.at(nt).SplittingFunctions;
        return cp10 * sf.at({1, 0}) + cp01 * sf.at({0, 1})
        + cp20 * sf.at({2, 0}) + cp11 * sf.at({1, 1}) + cp02 * sf.at({0, 2})
        + cp30 * sf.at({3, 0});
      };
    else if (PerturbativeOrder == 3)
      return [=] (int const& nt, double const& t) -> Set<Operator>
      {
        const double cp10 = Alphas(exp(t / 2)) / FourPi;
        const double cp01 = Alphaem(exp(t / 2)) / FourPi;
        const double cp20 = cp10 * cp10;
        const double cp11 = cp10 * cp01;
        const double cp02 = cp01 * cp01;
        const double cp30 = cp10 * cp20;
        const double cp40 = cp10 * cp30;
        const auto sf = DglapObj.at(nt).SplittingFunctions;
        return cp10 * sf.at({1, 0}) + cp01 * sf.at({0, 1})
        + cp20 * sf.at({2, 0}) + cp11 * sf.at({1, 1}) + cp02 * sf.at({0, 2})
        + cp30 * sf.at({3, 0})
        + cp40 * sf.at({4, 0});
      };
    else
      throw std::runtime_error(error("SplittingFunctionsQCDQED","Perturbative order not allowed."));
  }

  //_____________________________________________________________________________
  std::function<Set<Operator>(bool const&, int const&)> MatchingConditionsQCDQED(std::map<int, DglapObjectsQCDQED>        const& DglapObj,
                                                                                 int                                      const& PerturbativeOrder,
                                                                                 std::map<int, std::pair<double, double>> const& AlphasTh,
                                                                                 std::map<int, std::pair<double, double>> const& AlphaemTh)
  {
    if (PerturbativeOrder == 0)
      return [=] (bool const&, int const& nt) -> Set<Operator>
      {
        return DglapObj.at(nt).MatchingConditions.at({0, 0});
      };
    else if (PerturbativeOrder == 1)
      return [=] (bool const& Up, int const& nt) -> Set<Operator>
      {
        const double cp10 = (Up ? AlphasTh.at(nt+1).second  : AlphasTh.at(nt+1).first)  / FourPi;
        const double cp01 = (Up ? AlphaemTh.at(nt+1).second : AlphaemTh.at(nt+1).first) / FourPi;
        const auto mc = DglapObj.at(nt).MatchingConditions;
        return mc.at({0, 0})
        + (Up ? 1 : -1) * cp10 * mc.at({1, 0}) + (Up ? 1 : -1) * cp01 * mc.at({0, 1});
      };
    else if (PerturbativeOrder == 2)
      return [=] (bool const& Up, int const& nt) -> Set<Operator>
      {
        const double cp10 = (Up ? AlphasTh.at(nt+1).second  : AlphasTh.at(nt+1).first)  / FourPi;
        const double cp01 = (Up ? AlphaemTh.at(nt+1).second : AlphaemTh.at(nt+1).first) / FourPi;
        const double cp20 = cp10 * cp10;
        const auto mc = DglapObj.at(nt).MatchingConditions;
        return mc.at({0, 0})
        + (Up ? 1 : -1) * cp10 * mc.at({1, 0}) + (Up ? 1 : -1) * cp01 * mc.at({0, 1});
        + (Up ? 1 : -1) * cp20 * ( mc.at({2, 0}) - (Up ? 0 : 1) * mc.at({-2, 0}) );
      };
    else if (PerturbativeOrder == 3)
      return [=] (bool const& Up, int const& nt) -> Set<Operator>
      {
        const double cp10 = (Up ? AlphasTh.at(nt+1).second  : AlphasTh.at(nt+1).first)  / FourPi;
        const double cp01 = (Up ? AlphaemTh.at(nt+1).second : AlphaemTh.at(nt+1).first) / FourPi;
        const double cp20 = cp10 * cp10;
        const double cp30 = cp10 * cp20;
        const auto mc = DglapObj.at(nt).MatchingConditions;
        return mc.at({0, 0})
        + (Up ? 1 : -1) * cp10 * mc.at({1, 0}) + (Up ? 1 : -1) * cp01 * mc.at({0, 1});
        + (Up ? 1 : -1) * cp20 * ( mc.at({2, 0}) - (Up ? 0 : 1) * mc.at({-2, 0}) );
        + (Up ? 1 : -1) * cp30 * mc.at({3, 0});
      };
    else
      throw std::runtime_error(error("MatchingConditionsQCDQED","Perturbative order not allowed."));
  }


  //_____________________________________________________________________________
  std::unique_ptr<Dglap<Distribution>> BuildDglap(std::map<int, DglapObjectsQCDQED>                                  const& DglapObj,
                                                  std::function<std::map<int, double>(double const&, double const&)> const& InDistFunc,
                                                  double                                                             const& MuRef,
                                                  int                                                                const& PerturbativeOrder,
                                                  std::function<double(double const&)>                               const& Alphas,
                                                  std::function<double(double const&)>                               const& Alphaem,
                                                  int                                                                const& nsteps)
  {
    // Collect thresholds and coupling above and below them
    std::vector<double> Thresholds;
    std::map<int, std::pair<double, double>> AlphasTh;
    std::map<int, std::pair<double, double>> AlphaemTh;
    for (auto const& obj : DglapObj)
      {
        const int    nt  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) Thresholds.size() < nt)
          Thresholds.resize(nt);
        if (nt > 0)
          Thresholds[nt-1] = thr;
        AlphasTh.insert({nt, std::make_pair(Alphas(thr * ( 1 - eps8 )), Alphas(thr * ( 1 + eps8 )))});
        AlphaemTh.insert({nt, std::make_pair(Alphaem(thr * ( 1 - eps8 )), Alphaem(thr * ( 1 + eps8 )))});
      }

    // Create set of initial distributions
    const Set<Distribution> InPDFs{DglapObj.at(NF(MuRef, Thresholds)).SplittingFunctions.begin()->second.GetMap(),
                                   DistributionMap(DglapObj.begin()->second.SplittingFunctions.begin()->second.at(0).GetGrid(), InDistFunc, MuRef)};

    // Initialize DGLAP evolution
    return std::unique_ptr<Dglap<Distribution>>(new Dglap<Distribution> {SplittingFunctionsQCDQED(DglapObj, PerturbativeOrder, Alphas, Alphaem),
                                                                         MatchingConditionsQCDQED(DglapObj, PerturbativeOrder, AlphasTh, AlphaemTh),
                                                                         nullptr, InPDFs, MuRef, Thresholds, nsteps
                                                                        });
  }

  //_____________________________________________________________________________
  std::unique_ptr<Dglap<Operator>> BuildDglap(std::map<int, DglapObjectsQCDQED>    const& DglapObj,
                                              double                               const& MuRef,
                                              int                                  const& PerturbativeOrder,
                                              std::function<double(double const&)> const& Alphas,
                                              std::function<double(double const&)> const& Alphaem,
                                              int                                  const& nsteps)
  {
    // Collect thresholds and coupling above and below them
    std::vector<double> Thresholds;
    std::map<int, std::pair<double, double>> AlphasTh;
    std::map<int, std::pair<double, double>> AlphaemTh;
    for (auto const& obj : DglapObj)
      {
        const int    nt  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) Thresholds.size() < nt)
          Thresholds.resize(nt);
        if (nt > 0)
          Thresholds[nt-1] = thr;
        AlphasTh.insert({nt, std::make_pair(Alphas(thr * ( 1 - eps8 )), Alphas(thr * ( 1 + eps8 )))});
        AlphaemTh.insert({nt, std::make_pair(Alphaem(thr * ( 1 - eps8 )), Alphaem(thr * ( 1 + eps8 )))});
      }

    // Initialize DGLAP evolution. When computing evolution operators,
    // no inhomogeneous terms are allowed because their presence would
    // prevent wrinting the DGLAP evolution equations in terms of the
    // evolution operators. In other words, evolution operators can be
    // computed in the homogeneous case only. Set InhomogeneousTerms
    // to nullptr.
    // Initialize DGLAP evolution
    return std::unique_ptr<Dglap<Operator>>(new Dglap<Operator> {SplittingFunctionsQCDQED(DglapObj, PerturbativeOrder, Alphas, Alphaem),
                                                                 MatchingConditionsQCDQED(DglapObj, PerturbativeOrder, AlphasTh, AlphaemTh),
                                                                 nullptr, DglapObj.at(NF(MuRef, Thresholds)).UnitySet, MuRef, Thresholds, nsteps
                                                                });
  }
}
