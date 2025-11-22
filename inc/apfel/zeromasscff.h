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
   * @defgroup CFFMassless Zero-mass coefficient functions for the Compton form factors
   * Collection of the MSbar coefficient functions for the compton
   * form factors.
   */
  ///@{
  ///@}

  /**
   * @defgroup CFFLOzm LO zero-mass coefficient functions
   * @ingroup CFFMassless
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>0</SUP>) non-singlet real part of the
   * coefficient function for F1.
   */
  class ReCFF10ns: public Expression
  {
  public:
    ReCFF10ns();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double LocalPP(double const& x)  const;
  };
  ///@}

  /**
   * @defgroup CFFNLOzm NLO zero-mass coefficient functions
   * @ingroup CFFMassless
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet imaginary part of the
   * coefficient function for F1.
   */
  class ImCFF11ns: public Expression
  {
  public:
    ImCFF11ns();
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet real part of the
   * coefficient function for F1.
   */
  class ReCFF11ns: public Expression
  {
  public:
    ReCFF11ns();
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon imaginary part of the
   * coefficient function for F1.
   */
  class ImCFF11g: public Expression
  {
  public:
    ImCFF11g();
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon real part of the
   * coefficient function for F1.
   */
  class ReCFF11g: public Expression
  {
  public:
    ReCFF11g();
  };
  ///@}
}
