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
   * @defgroup KrkSF NLO unpolarised splitting functions in the Krk scheme
   * Two-loop (O(alpha_s^2)) space-like splitting functions in the Krk
   * factorisation scheme, following arXiv:2506.XXXXX (Delorme et al.).
   *
   * Extracted from PrecomputedSplittingFunctions.m via Pqiqi[Krk][2],
   * Pqiqk[Krk][2], Pqiqbi[Krk][2], Pqig[Krk][2], Pgq[Krk][2],
   * Pgg[Krk][2].
   *
   * Note: Pqiqbk[Krk][2] == Pqiqk[Krk][2]; use P1qiqkKrk for both.
   *
   * Convention: apfelxx uses alpha_s/(4*pi); the source file uses
   * alpha_s/(2*pi). All coefficients here are multiplied by 4.
   */
  ///@{

  /**
   * @brief P^{Krk}_{q_i q_i} at NLO (diagonal non-singlet plus).
   * Has plus-distribution and delta(1-x) contributions.
   */
  class P1qiqiKrk: public Expression
  {
  public:
    P1qiqiKrk(int const& nf);
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
   * @brief P^{Krk}_{q_i q_k} at NLO (k != i; also equals P^{Krk}_{q_i qbar_k}).
   * Purely regular (no distributions, no nf).
   */
  class P1qiqkKrk: public Expression
  {
  public:
    P1qiqkKrk();
    double Regular(double const& x) const;
  };

  /**
   * @brief P^{Krk}_{q_i qbar_i} at NLO (diagonal non-singlet minus contribution).
   * Purely regular (no distributions, no nf).
   */
  class P1qiqbiKrk: public Expression
  {
  public:
    P1qiqbiKrk();
    double Regular(double const& x) const;
  };

  /**
   * @brief P^{Krk}_{q_i g} at NLO (quark-gluon singlet channel, per flavour).
   * Purely regular.
   */
  class P1qigKrk: public Expression
  {
  public:
    P1qigKrk(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief P^{Krk}_{gq} at NLO (gluon-quark channel).
   * Purely regular.
   */
  class P1gqKrk: public Expression
  {
  public:
    P1gqKrk(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief P^{Krk}_{gg} at NLO (gluon-gluon channel).
   * Has plus-distribution and delta(1-x) contributions.
   */
  class P1ggKrk: public Expression
  {
  public:
    P1ggKrk(int const& nf);
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
