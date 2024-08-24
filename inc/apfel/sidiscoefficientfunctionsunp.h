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
   * @defgroup NCMasslessSIDIS Zero-mass coefficient functions for unpolarised SIDIS
   * Collection of the zero-mass coefficient functions for the
   * structure functions F<SUB>L</SUB> and F<SUB>T</SUB> in
   * unpolarised SIDIS up to O(&alpha;<SUB>s</SUB><SUP>2</SUP>).
   * @note Expressions are extracted from the following references:
   * https://arxiv.org/pdf/hep-ph/9711387 and
   * https://arxiv.org/pdf/2401.16281.
   */
  ///@{
  /**
   * @defgroup NLOzmSIDIS NLO zero-mass coefficient functions for unpolarised SIDIS
   * @ingroup NCMasslessSIDIS
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet quark-in-quak
   * coefficient function for FL
   */
  class C1LQ2Q: public DoubleExpression
  {
  public:
    C1LQ2Q();
    std::string GetName() const override { return "C1LQ2Q"; }
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon-in-quak coefficient function
   * for FL
   */
  class C1LQ2G: public DoubleExpression
  {
  public:
    C1LQ2G();
    std::string GetName() const override { return "C1LQ2G"; }
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) quark-in-gluon coefficient function
   * for FL
   */
  class C1LG2Q: public DoubleExpression
  {
  public:
    C1LG2Q();
    std::string GetName() const override { return "C1LG2Q"; }
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet quark-in-quark
   * coefficient function for FT
   */
  class C1TQ2Q: public DoubleExpression
  {
  public:
    C1TQ2Q();
    std::string GetName() const override { return "C1TQ2Q"; }
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
   * for FT
   */
  class C1TQ2G: public DoubleExpression
  {
  public:
    C1TQ2G();
    std::string GetName() const override { return "C1TQ2G"; }
    double LocalRegular(double const& x, double const& z) const override;
    double SingularRegular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) quark-in-gluon coefficient function
   * for FT
   */
  class C1TG2Q: public DoubleExpression
  {
  public:
    C1TG2Q();
    std::string GetName() const override { return "C1TG2Q"; }
    double RegularLocal(double const& x, double const& z) const override;
    double RegularSingular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };
  ///@}

  /**
   * @defgroup NNLOzmSIDIS NNLO zero-mass coefficient functions for unpolarised SIDIS
   * @ingroup NCMasslessSIDIS
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon-in-quark
   * coefficient function for FL.
   */
  class C2LQ2G: public DoubleExpression
  {
  public:
    C2LQ2G();
    std::string GetName() const override { return "C2LQ2G"; }
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon-in-gluon
   * coefficient function for FL.
   */
  class C2LG2G: public DoubleExpression
  {
  public:
    C2LG2G();
    std::string GetName() const override { return "C2LG2G"; }
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) quark-in-gluon
   * coefficient function for FL.
   */
  class C2LG2Q: public DoubleExpression
  {
  public:
    C2LG2Q();
    std::string GetName() const override { return "C2LG2Q"; }
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet
   * quark-in-quark coefficient function for FL.
   */
  class C2LQ2QNS: public DoubleExpression
  {
  public:
    C2LQ2QNS(int const& nf);
    std::string GetName() const override { return "C2LQ2QNS"; }
    double RegularRegular(double const& x, double const& z) const override;
  private:
    int const _nf;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * quark-in-quark coefficient function for FL.
   */
  class C2LQ2QPS: public DoubleExpression
  {
  public:
    C2LQ2QPS();
    std::string GetName() const override { return "C2LQ2QPS"; }
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * antiquark-in-quark coefficient function for FL.
   */
  class C2LQ2QB: public DoubleExpression
  {
  public:
    C2LQ2QB();
    std::string GetName() const override { return "C2LQ2QB"; }
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * quark(prime)-in-quark (1) coefficient function for FL.
   */
  class C2LQ2QP1: public DoubleExpression
  {
  public:
    C2LQ2QP1();
    std::string GetName() const override { return "C2LQ2QP1"; }
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * quark(prime)-in-quark (2) coefficient function for FL.
   */
  class C2LQ2QP2: public DoubleExpression
  {
  public:
    C2LQ2QP2();
    std::string GetName() const override { return "C2LQ2QP2"; }
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * quark(prime)-in-quark (1) coefficient function for FL.
   */
  class C2LQ2QP3: public DoubleExpression
  {
  public:
    C2LQ2QP3();
    std::string GetName() const override { return "C2LQ2QP3"; }
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon-in-quark
   * coefficient function for FT.
   */
  class C2TQ2G: public DoubleExpression
  {
  public:
    C2TQ2G();
    std::string GetName() const override { return "C2TQ2G"; }
    double LocalRegular(double const& x, double const& z) const override;
    double SingularRegular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon-in-gluon
   * coefficient function for FT.
   */
  class C2TG2G: public DoubleExpression
  {
  public:
    C2TG2G();
    std::string GetName() const override { return "C2TG2G"; }
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) quark-in-gluon
   * coefficient function for FT.
   */
  class C2TG2Q: public DoubleExpression
  {
  public:
    C2TG2Q();
    std::string GetName() const override { return "C2TG2Q"; }
    double RegularLocal(double const& x, double const& z) const override;
    double RegularSingular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet
   * quark-in-quark coefficient function for FT.
   */
  class C2TQ2QNS: public DoubleExpression
  {
  public:
    C2TQ2QNS(int const& nf);
    std::string GetName() const override { return "C2TQ2QNS"; }
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
   * quark-in-quark coefficient function for FT.
   */
  class C2TQ2QPS: public DoubleExpression
  {
  public:
    C2TQ2QPS();
    std::string GetName() const override { return "C2TQ2QPS"; }
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * antiquark-in-quark coefficient function for FT.
   */
  class C2TQ2QB: public DoubleExpression
  {
  public:
    C2TQ2QB();
    std::string GetName() const override { return "C2TQ2QB"; }
    double LocalRegular(double const& x, double const& z) const override;
    double SingularRegular(double const& x, double const& z) const override;
    double RegularLocal(double const& x, double const& z) const override;
    double RegularSingular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * quark(prime)-in-quark (1) coefficient function for FT.
   */
  class C2TQ2QP1: public DoubleExpression
  {
  public:
    C2TQ2QP1();
    std::string GetName() const override { return "C2TQ2QP1"; }
    double LocalRegular(double const& x, double const& z) const override;
    double SingularRegular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * quark(prime)-in-quark (2) coefficient function for FT.
   */
  class C2TQ2QP2: public DoubleExpression
  {
  public:
    C2TQ2QP2();
    std::string GetName() const override { return "C2TQ2QP2"; }
    double RegularLocal(double const& x, double const& z) const override;
    double RegularSingular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * quark(prime)-in-quark (1) coefficient function for FT.
   */
  class C2TQ2QP3: public DoubleExpression
  {
  public:
    C2TQ2QP3();
    std::string GetName() const override { return "C2TQ2QP3"; }
    double RegularRegular(double const& x, double const& z) const override;
  };
  ///@}
  ///@}
}
