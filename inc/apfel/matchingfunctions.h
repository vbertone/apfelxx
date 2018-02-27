//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/expression.h"

namespace apfel
{
  /**
   * @defgroup TMDMatchingFunctions TMD matching functions
   * @note The perturbative matching functions for PDFs and FFs
   * (references: https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf). Notice that the
   * expressions for FFs are implicitly multiplied by a factor
   * x<SUP>2</SUP>. This is because the convolution with FFs has the
   * form:
   *
   * \f$D(x) = C(x) \otimes d(x) / x^2 = (1/x^2) \int_x^1 (dy/y) [y^2 C(y)] d(x/y)\f$
   *
   * For the implementation of the O(&alpha;<SUB>s</SUB><SUP>2</SUP>)
   * expressions I have used the parameterization taken from the
   * arTeMiDe code but only for the regular part. The exact
   * coefficicients of the rest are implemented explicitly as reported
   * in eq. (B.3) of https://arxiv.org/pdf/1706.01473.pdf.
   */
  ///@{
  ///@}
  /**
   * @defgroup SLMatchFunc Space-like matching functions
   * @ingroup TMDMatchingFunctions
   */
  ///@{
  ///@}
  /**
   * @defgroup TLMatchFunc Time-like matching functions
   * @ingroup TMDMatchingFunctions
   */
  ///@{
  ///@}
  /**
   * @defgroup NLO NLO matching functions for PDFs
   * @ingroup SLMatchFunc
   */
  ///@{
  /**
   * @brief The O(&alpha;<SUB>s</SUB>) non-singlet matching function
   * for PDFs (references: https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C1ns: public Expression
  {
  public:
    C1ns();
    double Regular(double const& x) const;
    double Local(double const&)     const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) quark-gluon matching function
   * for PDFs (references: https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C1qg: public Expression
  {
  public:
    C1qg();
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) gluon-quark matching function
   * for PDFs (references: https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C1gq: public Expression
  {
  public:
    C1gq();
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) gluon-gluon matching function
   * for PDFs (references: https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C1gg: public Expression
  {
  public:
    C1gg();
    double Local(double const&) const;
  };
  ///@}

  /**
   * @defgroup NNLO NNLO matching functions for PDFs
   * @ingroup SLMatchFunc
   */
  ///@{
  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>2</SUP>) quark-quark
   * matching function for PDFs (references:
   * https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C2Vqq: public Expression
  {
  public:
    C2Vqq(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  protected:
    int const _nf;
    double    _A2;
    double    _A3;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>2</SUP>) quark-antiquark
   * matching function for PDFs (references:
   * https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C2Vqqb: public Expression
  {
  public:
    C2Vqqb();
    double Regular(double const& x)  const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * matching function for PDFs (references:
   * https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C2ps: public Expression
  {
  public:
    C2ps();
    double Regular(double const& x)  const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>2</SUP>) quark-gluon
   * matching function for PDFs (references:
   * https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C2qg: public Expression
  {
  public:
    C2qg();
    double Regular(double const& x)  const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon-quark
   * matching function for PDFs (references:
   * https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C2gq: public Expression
  {
  public:
    C2gq(int const& nf);
    double Regular(double const& x)  const;
  private:
    int const _nf;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon-gluon
   * matching function for PDFs (references:
   * https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C2gg: public Expression
  {
  public:
    C2gg(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
    double    _A2;
    double    _A3;
  };
  ///@}

  /**
   * @defgroup NLOff NLO matching functions for FFs
   * @ingroup TLMatchFunc
   */
  ///@{
  /**
   * @brief The O(&alpha;<SUB>s</SUB>) non-singlet matching function
   * for FFs (references: https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C1nsff: public Expression
  {
  public:
    C1nsff();
    double Regular(double const& x) const;
    double Local(double const&)     const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) quark-gluon matching function
   * for FFs (references: https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C1qgff: public Expression
  {
  public:
    C1qgff();
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) gluon-quark matching function
   * for FFs (references: https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C1gqff: public Expression
  {
  public:
    C1gqff();
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) gluon-gluon matching function
   * for FFs (references: https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C1ggff: public Expression
  {
  public:
    C1ggff();
    double Regular(double const& x) const;
    double Local(double const&)     const;
  };
  ///@}

  /**
   * @defgroup NNLOff NNLO matching functions for FFs
   * @ingroup TLMatchFunc
   */
  ///@{
  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>2</SUP>) quark-quark
   * matching function for FFs (references:
   * https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C2Vqqff: public Expression
  {
  public:
    C2Vqqff(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  protected:
    int const _nf;
    double    _A2;
    double    _A3;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>2</SUP>) quark-antiquark
   * matching function for FFs (references:
   * https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C2Vqqbff: public Expression
  {
  public:
    C2Vqqbff();
    double Regular(double const& x)  const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * matching function for FFs (references:
   * https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C2psff: public Expression
  {
  public:
    C2psff();
    double Regular(double const& x)  const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>2</SUP>) quark-gluon
   * matching function for FFs (references:
   * https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C2qgff: public Expression
  {
  public:
    C2qgff();
    double Regular(double const& x)  const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon-quark
   * matching function for FFs (references:
   * https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C2gqff: public Expression
  {
  public:
    C2gqff(int const& nf);
    double Regular(double const& x)  const;
  private:
    int const _nf;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon-gluon
   * matching function for FFs (references:
   * https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C2ggff: public Expression
  {
  public:
    C2ggff(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
    double    _A2;
    double    _A3;
  };
  ///@}
}
