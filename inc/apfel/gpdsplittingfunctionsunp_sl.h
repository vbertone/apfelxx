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
    double LocalPP(double const& c, double const& d) const;
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
    double LocalPP(double const& c, double const& d) const;
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
    int const _nf;
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
    double LocalPP(double const& c, double const& d) const;
  private:
    int const _nf;
  };
  ///@}

  /**
   * @defgroup NLOunpevk NLO evolution kernels
   * @ingroup GPDUnpSF
   * Two-loop non-singlet kernels in the HDR (hybrid) scheme. With
   * P(y, kappa) the quark-in-quark kernel P_qq^[1](x = y, xi = kappa y),
   * the "plus" combination is P(y, kappa) - P(-y, kappa) = P_qq^V + P_qqbar^V
   * and the "minus" (valence) combination is P(y, kappa) + P(-y, kappa) =
   * P_qq^V - P_qqbar^V. Pure-singlet contributions are not included.
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet-plus
   * unpolarised evolution kernel.
   */
  class Pgpd1nsp: public Expression
  {
  public:
    Pgpd1nsp(int const& nf, double const& xi);
    double Regular(double const& y)  const;
    double Singular(double const& y) const;
    double Local(double const& y)    const;
    double LocalPP(double const& c, double const& d) const;
  private:
    int const _nf;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet-minus
   * (valence) unpolarised evolution kernel.
   */
  class Pgpd1nsm: public Expression
  {
  public:
    Pgpd1nsm(int const& nf, double const& xi);
    double Regular(double const& y)  const;
    double Singular(double const& y) const;
    double Local(double const& y)    const;
    double LocalPP(double const& c, double const& d) const;
  private:
    int const _nf;
  };
  ///@}
}
