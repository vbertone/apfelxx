//
// APFEL++ 2017
//
// Author:Valerio Bertone:valerio.bertone@cern.ch
//

#pragma once

#include "apfel/doubleexpression.h"

namespace apfel
{
  /**
   * @defgroup NCMasslessSIDISpol Zero-mass coefficient functions for longitudinally polarised SIDIS
   * Collection of the zero-mass coefficient functions for the
   * structure function G<SUB>1</SUB> in longitudinally polarised
   * SIDIS up to O(&alpha;<SUB>s</SUB><SUP>2</SUP>).
   * @note Expressions are extracted from the following references:
   * https://arxiv.org/pdf/hep-ph/9711387 and
   * https://arxiv.org/pdf/2404.08597.
   */
  ///@{
  /**
   * @defgroup NLOzmSIDISpol NLO zero-mass coefficient functions for longitudinally polarised SIDIS
   * @ingroup NCMasslessSIDISpol
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet quark-in-quark
   * coefficient function for G1
   */
  class DC1Q2Q: public DoubleExpression
  {
  public:
    DC1Q2Q();
    std::string GetName() const override { return "DC1Q2Q"; }
    double LocalLocal(double const& x, double const& z) const override;
    double LocalSingular(double const& x, double const& z) const override;
    double LocalRegular(double const& x, double const& z) const override;
    double SingularLocal(double const& x, double const& z) const override;
    double SingularSingular(double const& x, double const& z) const override;
    double SingularRegular(double const& x, double const& z) const override;
    double RegularLocal(double const& x, double const& z) const override;
    double RegularSingular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon-in-quark coefficient function
   * for G1
   */
  class DC1Q2G: public DoubleExpression
  {
  public:
    DC1Q2G();
    std::string GetName() const override { return "DC1Q2G"; }
    double LocalRegular(double const& x, double const& z) const override;
    double SingularRegular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) quark-in-gluon coefficient function
   * for G1
   */
  class DC1G2Q: public DoubleExpression
  {
  public:
    DC1G2Q();
    std::string GetName() const override { return "DC1G2Q"; }
    double RegularLocal(double const& x, double const& z) const override;
    double RegularSingular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };
  ///@}

  /**
   * @defgroup NNLOzmSIDISpol NNLO zero-mass coefficient functions for longitudinally polarised SIDIS
   * @ingroup NCMasslessSIDISpol
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon-in-quark
   * coefficient function for G1.
   */
  class DC2Q2G: public DoubleExpression
  {
  public:
    DC2Q2G(int const& nf);
    std::string GetName() const override { return "DC2Q2G"; }
    double LocalRegular(double const& x, double const& z) const override;
    double SingularRegular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  private:
    int const _nf;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon-in-gluon
   * coefficient function for G1.
   */
  class DC2G2G: public DoubleExpression
  {
  public:
    DC2G2G();
    std::string GetName() const override { return "DC2G2G"; }
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) quark-in-gluon
   * coefficient function for G1.
   */
  class DC2G2Q: public DoubleExpression
  {
  public:
    DC2G2Q(int const& nf);
    std::string GetName() const override { return "DC2G2Q"; }
    double RegularLocal(double const& x, double const& z) const override;
    double RegularSingular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  private:
    int const _nf;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet
   * quark-in-quark coefficient function for G1.
   */
  class DC2Q2QNS: public DoubleExpression
  {
  public:
    DC2Q2QNS(int const& nf);
    std::string GetName() const override { return "DC2Q2QNS"; }
    double LocalLocal(double const& x, double const& z) const override;
    double LocalSingular(double const& x, double const& z) const override;
    double LocalRegular(double const& x, double const& z) const override;
    double SingularLocal(double const& x, double const& z) const override;
    double SingularSingular(double const& x, double const& z) const override;
    double SingularRegular(double const& x, double const& z) const override;
    double RegularLocal(double const& x, double const& z) const override;
    double RegularSingular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  private:
    int const _nf;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * quark-in-quark coefficient function for G1.
   */
  class DC2Q2QPS: public DoubleExpression
  {
  public:
    DC2Q2QPS();
    std::string GetName() const override { return "DC2Q2QPS"; }
    double RegularLocal(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * antiquark-in-quark coefficient function for G1.
   */
  class DC2Q2QB: public DoubleExpression
  {
  public:
    DC2Q2QB();
    std::string GetName() const override { return "DC2Q2QB"; }
    double LocalRegular(double const& x, double const& z) const override;
    double SingularRegular(double const& x, double const& z) const override;
    double RegularLocal(double const& x, double const& z) const override;
    double RegularSingular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * quark(prime)-in-quark (1) coefficient function for G1.
   */
  class DC2Q2QP1: public DoubleExpression
  {
  public:
    DC2Q2QP1();
    std::string GetName() const override { return "DC2Q2QP1"; }
    double LocalRegular(double const& x, double const& z) const override;
    double SingularRegular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * quark(prime)-in-quark (2) coefficient function for G1.
   */
  class DC2Q2QP2: public DoubleExpression
  {
  public:
    DC2Q2QP2();
    std::string GetName() const override { return "DC2Q2QP2"; }
    double RegularLocal(double const& x, double const& z) const override;
    double RegularSingular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * quark(prime)-in-quark (1) coefficient function for G1.
   */
  class DC2Q2QP3: public DoubleExpression
  {
  public:
    DC2Q2QP3();
    std::string GetName() const override { return "DC2Q2QP3"; }
    double RegularRegular(double const& x, double const& z) const override;
  };
  ///@}
  ///@}
}
