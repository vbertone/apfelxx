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
   * kernel for the DGLAP region.
   */
  class Pgpd0nsDGLAP: public Expression
  {
  public:
    Pgpd0nsDGLAP(double const& xi);
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    double const _xi;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet unpolarised evolution
   * kernel for the ERBL region.
   */
  class Pgpd0nsERBL: public Expression
  {
  public:
    Pgpd0nsERBL(double const& xi);
    double Singular(double const& x) const;
  private:
    double const _xi;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) quark-gluon unpolarised splitting
   * function for the DGLAP region.
   */
  class Pgpd0qgDGLAP: public Expression
  {
  public:
    Pgpd0qgDGLAP(double const& xi);
    double Regular(double const& x) const;
  private:
    double const _xi;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) quark-gluon unpolarised splitting
   * function for the ERBL region.
   */
  class Pgpd0qgERBL: public Expression
  {
  public:
    Pgpd0qgERBL(double const& xi);
    double Regular(double const& x) const;
  private:
    double const _xi;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon-quark unpolarised splitting
   * function for the DGLAP region.
   */
  class Pgpd0gqDGLAP: public Expression
  {
  public:
    Pgpd0gqDGLAP(double const& xi);
    double Regular(double const& x) const;
  private:
    double const _xi;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon-quark unpolarised splitting
   * function for the ERBL region.
   */
  class Pgpd0gqERBL: public Expression
  {
  public:
    Pgpd0gqERBL(double const& xi);
    double Regular(double const& x) const;
  private:
    double const _xi;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon-gluon unpolarised splitting
   * function for the DGLAP region.
   */
  class Pgpd0ggDGLAP: public Expression
  {
  public:
    Pgpd0ggDGLAP(int const& nf, double const& xi);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int    const _nf;
    double const _xi;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon-gluon unpolarised splitting
   * function for the ERBL region.
   */
  class Pgpd0ggERBL: public Expression
  {
  public:
    Pgpd0ggERBL(int const& nf, double const& xi);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int    const _nf;
    double const _xi;
  };
  ///@}
}
