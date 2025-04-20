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
        // Determine numer of active quarks
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
    const Operator O01ns{g, P0qedns{},   IntEps};
    const Operator O01gq{g, P0qedgmq{},  IntEps};
    const Operator O01qg{g, P0qedqgm{},  IntEps};
    const Operator O01gg{g, P0qedgmgm{}, IntEps};
    for (int nt = nti; nt <= ntf; nt++)
      {
        // Determine numer of active quarks and leptons
        const int nf = NDUL[nt][0] + NDUL[nt][1];
        const int nl = NDUL[nt][2];
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCDQED::PPDD, ed2 * O01ns});
        OM.insert({EvolutionBasisQCDQED::PPUU, eu2 * O01ns});
        OM.insert({EvolutionBasisQCDQED::PPLL,       O01ns});
        OM.insert({EvolutionBasisQCDQED::PMDD, ed2 * O01ns});
        OM.insert({EvolutionBasisQCDQED::PMUU, eu2 * O01ns});
        OM.insert({EvolutionBasisQCDQED::PMLL,       O01ns});
        OM.insert({EvolutionBasisQCDQED::PDgm, ( NC * ed2 ) * O01qg});
        OM.insert({EvolutionBasisQCDQED::PUgm, ( NC * eu2 ) * O01qg});
        OM.insert({EvolutionBasisQCDQED::PLgm,                O01qg});
        OM.insert({EvolutionBasisQCDQED::PgmD, ed2 * O01gq});
        OM.insert({EvolutionBasisQCDQED::PgmU, eu2 * O01gq});
        OM.insert({EvolutionBasisQCDQED::PgmL,       O01gq});
        OM.insert({EvolutionBasisQCDQED::Pgmgm, ( NC * SumCh2[nf] + nl ) * O01gg});
        // Insert Zero in the remaining slots
        for (int i = EvolutionBasisQCDQED::PPDD; i <= EvolutionBasisQCDQED::Pgmgm; i++)
          OM.insert({i, Zero});
        OpMap01.insert({nt, OM});
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
        if (nd - NDUL[std::min(nt + 1, ntf)][0] != 0)
          obj.Species = PartonSpecies::DOWNTYPEQUARK;
        else if (nu - NDUL[std::min(nt + 1, ntf)][1] != 0)
          obj.Species = PartonSpecies::UPTYPEQUARK;
        else if (nl - NDUL[std::min(nt + 1, ntf)][2] != 0)
          obj.Species = PartonSpecies::CHARGEDLEPTON;
        if (OpEvol)
          {
            std::map<int, Operator> MapUnity;
            for (auto const& coord : GkjQCDQED)
              MapUnity.insert({coord.second, (coord.first.first == coord.first.second ? Id : Zero)});
            obj.UnitySet = Set<Operator> {EvolutionOperatorBasisQCDQED{nd, nu, nl}, MapUnity};
            obj.SplittingFunctions.insert({{ 1, 0}, Set<Operator>{EvolutionOperatorBasisQCDQED{nd, nu, nl}, OpMap10.at(nt)}});
            obj.SplittingFunctions.insert({{ 0, 1}, Set<Operator>{EvolutionOperatorBasisQCDQED{nd, nu, nl}, OpMap01.at(nt)}});

            obj.MatchingConditions.insert({{ 0, 0}, Set<Operator>{MatchingOperatorBasisQCDQED{nd, nu, nl, obj.Species}, Match00}});
          }
        else
          {
            obj.SplittingFunctions.insert({{ 1, 0}, Set<Operator>{EvolutionBasisQCDQED{nd, nu, nl}, OpMap10.at(nt)}});
            obj.SplittingFunctions.insert({{ 0, 1}, Set<Operator>{EvolutionBasisQCDQED{nd, nu, nl}, OpMap01.at(nt)}});

            obj.MatchingConditions.insert({{ 0, 0}, Set<Operator>{MatchingBasisQCDQED{nd, nu, nl, obj.Species}, Match00}});
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
