// APFEL++ 2017

#pragma once

#include "apfel/doubleoperator.h"
#include "apfel/doubledistribution.h"
#include "apfel/distribution.h"
#include "apfel/grid.h"

#include <functional>
#include <map>
#include <string>
#include <vector>

namespace apfel
{
  /**
   * @brief Structure containing pre-computed DoubleOperator objects for
   * all SIDIS coefficient functions up to exact NNLO.
   *
   * Operators are organised in three maps by structure function type:
   * - @c FT: unpolarised transverse (= F<SUB>1</SUB>) operators
   * - @c FL: unpolarised longitudinal operators
   * - @c G1: longitudinally polarised G<SUB>1</SUB> operators
   *
   * Within each map, entries are keyed by the name of the corresponding
   * DoubleExpression (i.e. matching DoubleExpression::GetName()), e.g.:
   * - @c "DoubleIdentity"  → LO (FT and G1 only)
   * - @c "C1TQ2Q"          → NLO FT qq channel
   * - @c "C2TQ2QNS_nf3"    → NNLO FT non-singlet at nf = 3
   * - @c "DC2Q2G_nf5"      → NNLO G1 gluon-in-quark at nf = 5
   *
   * @note The DoubleOperator objects hold references to the Grid objects
   * passed to InitializeSidisObjects(). Those grids must remain alive
   * for the lifetime of this structure.
   */
  struct SidisNNLOObjects
  {
    std::map<std::string, DoubleOperator> FT;
    std::map<std::string, DoubleOperator> FL;
    std::map<std::string, DoubleOperator> G1;
  };

  /**
   * @brief Initialises all SIDIS DoubleOperator objects (LO + NLO + exact NNLO)
   * for the unpolarised (FT, FL) and polarised (G1) structure functions.
   *
   * This is the exact NNLO implementation. It uses the DoubleExpression
   * subclasses from sidiscoefficientfunctionsunp.h and sidiscoefficientfunctionspol.h
   * directly.
   *
   * Coefficient function references:
   * - NLO:  hep-ph/9711387 Appendix C
   * - NNLO unpolarised:  arXiv:2401.16281
   * - NNLO polarised:    arXiv:2404.08597
   *
   * @param gx: x-space grid (must outlive the returned object)
   * @param gz: z-space grid (must outlive the returned object)
   * @param Thresholds: heavy-quark thresholds determining the active-nf range
   * @param IntEps: integration accuracy for DoubleOperator construction (default: 1e-3)
   * @return SidisNNLOObjects with all operators indexed by channel name
   */
  SidisNNLOObjects InitializeSidisObjects(Grid                const& gx,
                                          Grid                const& gz,
                                          std::vector<double> const& Thresholds,
                                          double              const& IntEps = 1e-3);

  /**
   * @brief Builds the unpolarised transverse (F<SUB>T</SUB> = F<SUB>1</SUB>) SIDIS
   * cross-section function from pre-computed operators.
   *
   * The returned function maps a scale Q to a DoubleDistribution f(x, z) such
   * that the physical cross section is obtained via
   * @code
   *   auto FT = BuildSidisUnpFT(obj, InPDFs, InFFs, Alphas, Thresholds);
   *   double xsec = FT(Q).Evaluate1(x).Integrate(zmin, zmax);
   * @endcode
   * The 1/(xz) kinematic prefactor is included.
   *
   * @param obj: pre-computed operators from InitializeSidisObjects (must outlive the returned function)
   * @param InPDFs: unpolarised PDFs as xf(x, Q), mapping flavour id to value
   * @param InFFs: fragmentation functions as zD(z, Q), mapping flavour id to value
   * @param Alphas: strong coupling α<SUB>s</SUB>(Q)
   * @param Thresholds: heavy-quark thresholds (used to determine nf at each Q)
   * @param PerturbativeOrder: 0 = LO, 1 = NLO, 2 = NNLO (default: 2); only used when channel = -1
   * @param channel: partonic channel selector (default: -1 = full sum controlled by PerturbativeOrder).
   *   Individual channels follow the benchmark convention:
   *   0 = LO ns; 1 = LO+NLO ns; 2 = NLO qg; 3 = NLO gq;
   *   4 = ns+ps through NNLO; 5 = qg through NNLO; 6 = gq through NNLO;
   *   7 = NNLO gg; 8 = NNLO qbq; 9 = NNLO qpq.
   * @return function Q → F<SUB>T</SUB> DoubleDistribution
   */
  std::function<DoubleDistribution(double const&)> BuildSidisUnpFT(
    SidisNNLOObjects                                                          const& obj,
    std::function<std::map<int, double>(double const&, double const&)> const& InPDFs,
    std::function<std::map<int, double>(double const&, double const&)> const& InFFs,
    std::function<double(double const&)>                               const& Alphas,
    std::vector<double>                                                const& Thresholds,
    int                                                                const& PerturbativeOrder = 2,
    int                                                                const& channel = -1);

  /**
   * @brief Builds the unpolarised longitudinal (F<SUB>L</SUB>) SIDIS
   * cross-section function from pre-computed operators.
   *
   * Identical API to BuildSidisUnpFT but for F<SUB>L</SUB>.
   * F<SUB>L</SUB> has no LO contribution; PerturbativeOrder = 0 returns zero.
   *
   * @param obj: pre-computed operators from InitializeSidisObjects (must outlive the returned function)
   * @param InPDFs: unpolarised PDFs as xf(x, Q)
   * @param InFFs: fragmentation functions as zD(z, Q)
   * @param Alphas: strong coupling α<SUB>s</SUB>(Q)
   * @param Thresholds: heavy-quark thresholds
   * @param PerturbativeOrder: 0 = 0 (no LO), 1 = NLO, 2 = NNLO (default: 2)
   * @return function Q → F<SUB>L</SUB> DoubleDistribution
   */
  std::function<DoubleDistribution(double const&)> BuildSidisUnpFL(
    SidisNNLOObjects                                                          const& obj,
    std::function<std::map<int, double>(double const&, double const&)> const& InPDFs,
    std::function<std::map<int, double>(double const&, double const&)> const& InFFs,
    std::function<double(double const&)>                               const& Alphas,
    std::vector<double>                                                const& Thresholds,
    int                                                                const& PerturbativeOrder = 2);

  /**
   * @brief Builds the longitudinally polarised G<SUB>1</SUB> SIDIS
   * cross-section function from pre-computed operators.
   *
   * Identical API to BuildSidisUnpFT but for G<SUB>1</SUB>.
   * @c InPDFs should provide the polarised (helicity) PDFs Δf(x, Q).
   *
   * @param obj: pre-computed operators from InitializeSidisObjects (must outlive the returned function)
   * @param InPDFs: polarised PDFs as xΔf(x, Q)
   * @param InFFs: fragmentation functions as zD(z, Q)
   * @param Alphas: strong coupling α<SUB>s</SUB>(Q)
   * @param Thresholds: heavy-quark thresholds
   * @param PerturbativeOrder: 0 = LO, 1 = NLO, 2 = NNLO (default: 2); only used when channel = -1
   * @param channel: partonic channel selector (default: -1 = full sum); see BuildSidisUnpFT for convention
   * @return function Q → G<SUB>1</SUB> DoubleDistribution
   */
  std::function<DoubleDistribution(double const&)> BuildSidisG1(
    SidisNNLOObjects                                                          const& obj,
    std::function<std::map<int, double>(double const&, double const&)> const& InPDFs,
    std::function<std::map<int, double>(double const&, double const&)> const& InFFs,
    std::function<double(double const&)>                               const& Alphas,
    std::vector<double>                                                const& Thresholds,
    int                                                                const& PerturbativeOrder = 2,
    int                                                                const& channel = -1);
}
