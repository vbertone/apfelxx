//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/expression.h"

namespace apfel
{
  /**
   * @defgroup PhysSF NLO unpolarised splitting functions in the PHYS scheme
   * Two-loop (O(alpha_s^2)) space-like splitting functions in the PHYS
   * factorisation scheme, following arXiv:2506.XXXXX (Delorme et al.).
   *
   * Extracted from PrecomputedSplittingFunctions.m via Pqiqi[PHYS][2],
   * Pqiqk[PHYS][2], Pqiqbi[PHYS][2], Pqig[PHYS][2], Pgq[PHYS][2],
   * Pgg[PHYS][2].
   *
   * Note: Pqiqbk[PHYS][2] == Pqiqk[PHYS][2]; use P1qiqkPHYS for both.
   *
   * Convention: apfelxx uses alpha_s/(4*pi); the source file uses
   * alpha_s/(2*pi). All coefficients here are multiplied by 4.
   */
  ///@{

  /**
   * @brief P^{PHYS}_{q_i q_i} at NLO (diagonal non-singlet plus).
   * Has plus-distribution and delta(1-x) contributions.
   */
  class P1qiqiPHYS: public Expression
  {
  public:
    P1qiqiPHYS(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int    const _nf;
    double       _a0;    ///< coefficient of [1/(1-x)]_+
    double       _a1;    ///< coefficient of [log(1-x)/(1-x)]_+
    double       _delta; ///< coefficient of delta(1-x)
  };

  /**
   * @brief P^{PHYS}_{q_i q_k} at NLO (k != i; also equals P^{PHYS}_{q_i qbar_k}).
   * Purely regular (no distributions, no nf).
   */
  class P1qiqkPHYS: public Expression
  {
  public:
    P1qiqkPHYS();
    double Regular(double const& x) const;
  };

  /**
   * @brief P^{PHYS}_{q_i qbar_i} at NLO (diagonal non-singlet minus contribution).
   * Purely regular (no distributions, no nf).
   */
  class P1qiqbiPHYS: public Expression
  {
  public:
    P1qiqbiPHYS();
    double Regular(double const& x) const;
  };

  /**
   * @brief P^{PHYS}_{q_i g} at NLO (quark-gluon singlet channel, per flavour).
   * Purely regular. Uses Li_2(x) and Li_2(x^2).
   */
  class P1qigPHYS: public Expression
  {
  public:
    P1qigPHYS(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief P^{PHYS}_{gq} at NLO (gluon-quark channel).
   * Purely regular. Uses Li_2(x) and Li_2(x^2).
   */
  class P1gqPHYS: public Expression
  {
  public:
    P1gqPHYS(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief P^{PHYS}_{gg} at NLO (gluon-gluon channel).
   * Has plus-distribution and delta(1-x) contributions.
   */
  class P1ggPHYS: public Expression
  {
  public:
    P1ggPHYS(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int    const _nf;
    double       _a0;    ///< coefficient of [1/(1-x)]_+
    double       _a1;    ///< coefficient of [log(1-x)/(1-x)]_+
    double       _delta; ///< coefficient of delta(1-x)
  };
  ///@}
}
