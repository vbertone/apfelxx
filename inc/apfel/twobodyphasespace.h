//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Fulvio Piacenza: fulvio.piacenza01@universitadipavia.it
//

#pragma once

namespace apfel
{
  /**
   * @brief Class for the calculation of the phase-space reduction
   * factor due to cuts on the single outgoing lepton in Drell-Yan
   * production. The relevant process is:<br>
   * 		&gamma;(q) &rarr; l<SUP>+</SUP>(k<SUB>1</SUB>) + l<SUP>-</SUP>(k<SUB>2</SUB>) <br>
   *  with:<br>
   * 		k<SUB>T,1(2)</SUB> > p<SUB>T,min</SUB><br>
   * 		|&eta;<SUB>1(2)</SUB>| < &eta;<SUB>max</SUB><br>
   */
  class TwoBodyPhaseSpace
  {
  public:
    /**
     * @brief The "TwoBodyPhaseSpace" constructor for asymmetric cuts.
     * @param pTmin1: the minimum cut in the p<SUB>T</SUB> of the first single lepton
     * @param pTmin2: the minimum cut in the p<SUB>T</SUB> of the second single lepton
     * @param etamin: the minimum cut in the &eta; of the single lepton
     * @param etamax: the maximum cut in the &eta; of the single lepton
     * @param eps: the integration accuracy
     */
    TwoBodyPhaseSpace(double const& pTmin1, double const& pTmin2, double const& etamin, double const& etamax, double const& eps);

    /**
     * @brief The "TwoBodyPhaseSpace" constructor for symmetric cuts.
     * @param pTmin: the minimum cut in the p<SUB>T</SUB> of the single lepton
     * @param etamin: the minimum cut in the &eta; of the single lepton
     * @param etamax: the maximum cut in the &eta; of the single lepton
     * @param eps: the integration accuracy (default = 10<SUP>-9</SUP>)
     */
    TwoBodyPhaseSpace(double const& pTmin, double const& etamin, double const& etamax, double const& eps = 1e-9);

    /**
     * @brief Function that returns the phase-space reduction factor.
     * @param Q: invariant mass of the leptonic pair
     * @param y: rapidity of the leptonic pair
     * @param qT: transverse momentum of the leptonic pair
     * @return the phase-space reduction factor as a function of the
     * invariant mass, transverse momentum and rapidity of the lepton
     * pair.
     */
    double PhaseSpaceReduction(double const& Q, double const& y, double const& qT);

    /**
     * @brief Function that returns the derivative w.r.t. qT of the
     * phase-space reduction factor.
     * @param Q: invariant mass of the leptonic pair
     * @param y: rapidity of the leptonic pair
     * @param qT: transverse momentum of the leptonic pair
     * @return the derivative of the phase-space reduction factor as a
     * function of the invariant mass, transverse momentum and rapidity
     * of the lepton pair.
     */
    double DerivePhaseSpaceReduction(double const& Q, double const& y, double const& qT);

    /**
     * @brief Function that returns the phase-space reduction factor
     * associated to the parity violating contribution.
     * @param Q: invariant mass of the leptonic pair
     * @param y: rapidity of the leptonic pair
     * @param qT: transverse momentum of the leptonic pair
     * @return the phase-space reduction factor as a function of the
     * invariant mass, transverse momentum and rapidity of the lepton
     * pair.
     */
    double ParityViolatingPhaseSpaceReduction(double const& Q, double const& y, double const& qT);

  private:
    /**
     * @name Cut variables
     * Cut Variables that define the cuts on the single leptons
     */
    ///@{
    double const _pTmin1;
    double const _pTmin2;
    double const _etamin;
    double const _etamax;
    double const _eps;
    ///@}
  };
}
