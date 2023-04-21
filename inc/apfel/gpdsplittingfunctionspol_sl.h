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
   * @defgroup GPDPolSF Polarised evolution kernels
   * @ingroup GPDEvKernels
   */
  ///@{
  ///@}
  /**
   * @defgroup LOpolevk LO evolution kernels
   * @ingroup GPDPolSF
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet polarised evolution
   * kernel.
   */
  class Pgpd0polns: public Expression
  {
  public:
    Pgpd0polns(double const& xi);
    double Regular(double const& y)  const;
    double Singular(double const& y) const;
    double Local(double const& y)    const;
    double LocalPP(double const& y)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) quark-quark polarised splitting
   * function.
   */
  class Pgpd0polqq: public Expression
  {
  public:
    Pgpd0polqq(double const& xi);
    double Regular(double const& y)  const;
    double Singular(double const& y) const;
    double Local(double const& y)    const;
    double LocalPP(double const& y)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) quark-gluon polarised splitting
   * function.
   */
  class Pgpd0polqg: public Expression
  {
  public:
    Pgpd0polqg(int const& nf, double const& xi);
    double Regular(double const& y) const;
  private:
    int const _nf;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon-quark polarised splitting
   * function.
   */
  class Pgpd0polgq: public Expression
  {
  public:
    Pgpd0polgq(double const& xi);
    double Regular(double const& y) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon-gluon polarised splitting
   * function.
   */
  class Pgpd0polgg: public Expression
  {
  public:
    Pgpd0polgg(int const& nf, double const& xi);
    double Regular(double const& y)    const;
    double Singular(double const& y)   const;
    double Local(double const& y)      const;
    double LocalPP(double const& y)    const;
    double SingularPV(double const& y) const;
    double LocalPV(double const& y)    const;
  private:
    int const _nf;
  };
  ///@}
}
