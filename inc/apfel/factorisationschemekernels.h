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
   * @defgroup FSKernels Factorisation-scheme transition kernels
   * NLO (O(alpha_s)) kernels K_{ij}(x) for transitions between
   * factorisation schemes, following arXiv:2501.18289, eqs. (69)-(72).
   *
   * Each kernel is decomposed as
   *   K(x) = [a0/(1-x)]_+ + [a1*ln(1-x)/(1-x)]_+ + b*ln(1-x)
   *           + c*ln(x) + P + Delta*delta(1-x)
   *
   * mapped to apfel::Expression methods as:
   *   Regular(x)  = b(x)*ln(1-x) + c(x)*ln(x) + P(x)
   *   Singular(x) = (a0 + a1*ln(1-x)) / (1-x)
   *   Local(x)    = Delta + a0*ln(1-x) + (a1/2)*ln^2(1-x)
   *
   * Channels are indexed as: 1=qq, 2=qg, 3=gq, 4=gg.
   * Classes that vanish identically use apfel::Null.
   */
  ///@{
  // ===========================================================================
  // POS scheme
  // ===========================================================================
  ///@{
  /// @brief K^POS_{qg} at NLO (qq and gg vanish: use Null)
  class K1qgPOS: public Expression
  {
  public:
    K1qgPOS();
    double Regular(double const& x) const;
  };

  /// @brief K^POS_{gq} at NLO
  class K1gqPOS: public Expression
  {
  public:
    K1gqPOS();
    double Regular(double const& x) const;
  };
  ///@}

  // ===========================================================================
  // MPOS scheme
  // ===========================================================================
  ///@{
  /// @brief K^MPOS_{qq} at NLO
  class K1qqMPOS: public Expression
  {
  public:
    K1qqMPOS();
    double Regular(double const& x) const;
  };

  /// @brief K^MPOS_{qg} at NLO (same b,c,P as POS)
  class K1qgMPOS: public Expression
  {
  public:
    K1qgMPOS();
    double Regular(double const& x) const;
  };

  /// @brief K^MPOS_{gq} at NLO (same b,c,P as POS)
  class K1gqMPOS: public Expression
  {
  public:
    K1gqMPOS();
    double Regular(double const& x) const;
  };

  /// @brief K^MPOS_{gg} at NLO
  class K1ggMPOS: public Expression
  {
  public:
    K1ggMPOS(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };
  ///@}

  // ===========================================================================
  // MPOSDelta scheme
  // ===========================================================================
  ///@{
  /// @brief K^{MPOSDelta}_{qq} at NLO (pure delta contribution)
  class K1qqMPOSDelta: public Expression
  {
  public:
    K1qqMPOSDelta();
    double Local(double const&) const;
  };

  /// @brief K^{MPOSDelta}_{qg} at NLO (same b,c,P as POS)
  class K1qgMPOSDelta: public Expression
  {
  public:
    K1qgMPOSDelta();
    double Regular(double const& x) const;
  };

  /// @brief K^{MPOSDelta}_{gq} at NLO (same b,c,P as POS)
  class K1gqMPOSDelta: public Expression
  {
  public:
    K1gqMPOSDelta();
    double Regular(double const& x) const;
  };

  /// @brief K^{MPOSDelta}_{gg} at NLO (pure delta contribution)
  class K1ggMPOSDelta: public Expression
  {
  public:
    K1ggMPOSDelta(int const& nf);
    double Local(double const&) const;
  private:
    int const _nf;
  };
  ///@}

  // ===========================================================================
  // Krk scheme
  // ===========================================================================
  ///@{
  /// @brief K^{Krk}_{qq} at NLO
  class K1qqKrk: public Expression
  {
  public:
    K1qqKrk();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /// @brief K^{Krk}_{qg} at NLO
  class K1qgKrk: public Expression
  {
  public:
    K1qgKrk();
    double Regular(double const& x) const;
  };

  /// @brief K^{Krk}_{gq} at NLO
  class K1gqKrk: public Expression
  {
  public:
    K1gqKrk();
    double Regular(double const& x) const;
  };

  /// @brief K^{Krk}_{gg} at NLO
  class K1ggKrk: public Expression
  {
  public:
    K1ggKrk(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };
  ///@}

  // ===========================================================================
  // Aversa scheme
  // ===========================================================================
  ///@{
  /// @brief K^{Aversa}_{qq} at NLO
  class K1qqAversa: public Expression
  {
  public:
    K1qqAversa();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /// @brief K^{Aversa}_{qg} at NLO
  class K1qgAversa: public Expression
  {
  public:
    K1qgAversa();
    double Regular(double const& x) const;
  };

  /// @brief K^{Aversa}_{gq} at NLO
  class K1gqAversa: public Expression
  {
  public:
    K1gqAversa();
    double Regular(double const& x) const;
  };

  /// @brief K^{Aversa}_{gg} at NLO
  class K1ggAversa: public Expression
  {
  public:
    K1ggAversa(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };
  ///@}

  // ===========================================================================
  // DIS scheme
  // ===========================================================================
  ///@{
  /// @brief K^{DIS}_{qq} at NLO
  class K1qqDIS: public Expression
  {
  public:
    K1qqDIS();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /// @brief K^{DIS}_{qg} at NLO
  class K1qgDIS: public Expression
  {
  public:
    K1qgDIS();
    double Regular(double const& x) const;
  };

  /// @brief K^{DIS}_{gq} at NLO
  class K1gqDIS: public Expression
  {
  public:
    K1gqDIS();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /// @brief K^{DIS}_{gg} at NLO (gg vanishes for DIS: a0=a1=Delta=0; only b,c,P nonzero)
  class K1ggDIS: public Expression
  {
  public:
    K1ggDIS(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };
  ///@}

  // ===========================================================================
  // PHYS scheme
  // ===========================================================================
  ///@{
  /// @brief K^{PHYS}_{qq} at NLO
  class K1qqPHYS: public Expression
  {
  public:
    K1qqPHYS();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /// @brief K^{PHYS}_{qg} at NLO
  class K1qgPHYS: public Expression
  {
  public:
    K1qgPHYS();
    double Regular(double const& x) const;
  };

  /// @brief K^{PHYS}_{gq} at NLO
  class K1gqPHYS: public Expression
  {
  public:
    K1gqPHYS();
    double Regular(double const& x) const;
  };

  /// @brief K^{PHYS}_{gg} at NLO
  class K1ggPHYS: public Expression
  {
  public:
    K1ggPHYS(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };
  ///@}
  ///@}
}
