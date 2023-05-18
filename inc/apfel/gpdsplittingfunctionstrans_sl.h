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
   * @defgroup GPDTransSF Transversely polarised evolution kernels
   * @ingroup GPDEvKernels
   */
  ///@{
  ///@}
  /**
   * @defgroup LOtransevk LO evolution kernels
   * @ingroup GPDTransSF
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet polarised evolution
   * kernel.
   */
  class Pgpd0transns: public Expression
  {
  public:
    Pgpd0transns(double const& xi);
    double Regular(double const& y)  const;
    double Singular(double const& y) const;
    double Local(double const& y)    const;
    double LocalPP(double const& y)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) quark-quark transversely polarised
   * splitting function.
   */
  class Pgpd0transqq: public Expression
  {
  public:
    Pgpd0transqq(double const& xi);
    double Regular(double const& y)  const;
    double Singular(double const& y) const;
    double Local(double const& y)    const;
    double LocalPP(double const& y)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon-gluon linearly polarised
   * splitting function.
   */
  class Pgpd0transgg: public Expression
  {
  public:
    Pgpd0transgg(int const& nf, double const& xi);
    double Regular(double const& y)    const;
    double Singular(double const& y)   const;
    double Local(double const& y)      const;
    double LocalPP(double const& y)    const;
  private:
    int const _nf;
  };
  ///@}
}
