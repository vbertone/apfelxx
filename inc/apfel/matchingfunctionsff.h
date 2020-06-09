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
   * @defgroup TLMatchFunc Time-like matching functions
   * TMD mathing functions for FFs (time-like hard scales)
   * @ingroup TMDMatchingFunctions
   */
  ///@{
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
    double Regular(double const& x) const;
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
    double Regular(double const& x) const;
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
    double Regular(double const& x) const;
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
    double Regular(double const& x) const;
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
