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
   * @defgroup GPDEvKernels GPD evolution kernels
   * Collection of the MSbar evolution kernels for the evolution of
   * GPDs. They are split into two categories that apply to the DGLAP
   * and ERBL regions respectively.
   */
  ///@{
  ///@}
  /**
   * @defgroup GPDUnpSF Unpolarised evolution kernels
   * @ingroup GPDEvKernels
   */
  ///@{
  ///@}
  /**
   * @defgroup LOunpevk LO evolution kernels
   * @ingroup GPDUnpSF
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet unpolarised evolution
   * kernel.
   */
  class Pgpd0ns: public Expression
  {
  public:
    Pgpd0ns(double const& xi);
    double Regular(double const& y)  const;
    double Singular(double const& y) const;
    double Local(double const& y)    const;
    double LocalPV(double const& y)  const;
  private:
    double const _xi;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) quark-quark unpolarised splitting
   * function.
   */
  class Pgpd0qq: public Expression
  {
  public:
    Pgpd0qq(double const& xi);
    double Regular(double const& y)  const;
    double Singular(double const& y) const;
    double Local(double const& y)    const;
    double LocalPV(double const& y)  const;
  private:
    double const _xi;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) quark-gluon unpolarised splitting
   * function.
   */
  class Pgpd0qg: public Expression
  {
  public:
    Pgpd0qg(int const& nf, double const& xi);
    double Regular(double const& y) const;
  private:
    int    const _nf;
    double const _xi;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon-quark unpolarised splitting
   * function.
   */
  class Pgpd0gq: public Expression
  {
  public:
    Pgpd0gq(double const& xi);
    double Regular(double const& y) const;
  private:
    double const _xi;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon-gluon unpolarised splitting
   * function.
   */
  class Pgpd0gg: public Expression
  {
  public:
    Pgpd0gg(int const& nf, double const& xi);
    double Regular(double const& y)  const;
    double Singular(double const& y) const;
    double Local(double const& y)    const;
    double LocalPV(double const& y)  const;
  private:
    int    const _nf;
    double const _xi;
  };
  ///@}
}
