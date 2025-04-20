//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/expression.h"

#include <vector>

namespace apfel
{
  /**
   * @defgroup UnpSFQED Unpolarised splitting functions for QCDxQED evolution
   * @ingroup SLSplittings
   */
  ///@{
  ///@}
  /**
   * @defgroup LOunpsfQED LO splitting functions for QCDxQED
   * @ingroup UnpSFQED
   */
  ///@{
  /**
   * @brief Space-like O(&alpha;) non-singlet unpolarised splitting
   * function.
   */
  class P01qedns: public Expression
  {
  public:
    P01qedns();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief Space-like O(&alpha;) quark-gamma unpolarised splitting
   * function.
   */
  class P01qedqgm: public Expression
  {
  public:
    P01qedqgm();
    double Regular(double const& x) const;
  };

  /**
   * @brief Space-like O(&alpha;) gamma-quark unpolarised splitting
   * function.
   */
  class P01qedgmq: public Expression
  {
  public:
    P01qedgmq();
    double Regular(double const& x) const;
  };

  /**
   * @brief Space-like O(&alpha;) gamma-gamma unpolarised splitting
   * function.
   */
  class P01qedgmgm: public Expression
  {
  public:
    P01qedgmgm();
    double Local(double const& x) const;
  };
  ///@}

  /**
   * @defgroup NLOunpsfQED NLO splitting functions for QCDxQED
   * @ingroup UnpSFQED
   */
  ///@{
  /**
   * @brief Space-like O(&alpha;<SUP>2</SUP>) non-singlet-plus
   * unpolarised splitting function.
   */
  class P02qednsp: public Expression
  {
  public:
    P02qednsp(double const& crat2);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  protected:
    double const _crat2;
    double       _a2;
  };

  /**
   * @brief Space-like O(&alpha;<SUP>2</SUP>) non-singlet-minus
   * unpolarised splitting function.
   */
  class P02qednsm: public Expression
  {
  public:
    P02qednsm(double const& crat2);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  protected:
    double const _crat2;
    double       _a2;
  };

  /**
   * @brief Space-like O(&alpha;<SUP>2</SUP>) pure-singlet
   * unpolarised splitting function.
   */
  class P02qedps: public Expression
  {
  public:
    P02qedps();
    double Regular(double const& x) const;
  };

  /**
   * @brief Space-like O(&alpha;<SUP>2</SUP>) quark-photon unpolarised
   * splitting function.
   */
  class P02qedqgm: public Expression
  {
  public:
    P02qedqgm();
    double Regular(double const& x) const;
  };

  /**
   * @brief Space-like O(&alpha;<SUP>2</SUP>) photon-quark unpolarised
   * splitting function.
   */
  class P02qedgmq: public Expression
  {
  public:
    P02qedgmq(double const& crat2);
    double Regular(double const& x) const;
  private:
    int const _crat2;
  };

  /**
   * @brief Space-like O(&alpha;<SUP>2</SUP>) photon-photon unpolarised
   * splitting function.
   */
  class P02qedgmgm: public Expression
  {
  public:
    P02qedgmgm();
    double Regular(double const& x) const;
    double Local(double const& x)   const;
  };
  ///@}
}
