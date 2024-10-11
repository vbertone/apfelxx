//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/expression.h"
#include "apfel/splittingfunctionsunp_tl.h"

namespace apfel
{
  /**
   * @defgroup PolSFtl Longitudinally polarised splitting functions
   * @ingroup TLSplittings
   */
  ///@{
  ///@}
  /**
   * @defgroup LOpolsf LO splitting functions
   * @ingroup PolSF
   */
  ///@{
  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB>) non-singlet
   * longitudinally polarised splitting function. This is equal to the
   * non-singlet unpolarised splitting function.
   */
  class P0Tpolns: public P0Tns
  {
  public:
    P0Tpolns();
  };

  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB>) quark-gluon
   * longitudinally polarised splitting function.
   */
  class P0Tpolqg: public Expression
  {
  public:
    P0Tpolqg(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB>) gluon-quark
   * longitudinally polarised splitting function.
   */
  class P0Tpolgq: public Expression
  {
  public:
    P0Tpolgq();
    double Regular(double const& x) const;
  };

  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB>) gluon-gluon
   * longitudinally polarised splitting function.
   */
  class P0Tpolgg: public Expression
  {
  public:
    P0Tpolgg(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };
  ///@}
}
