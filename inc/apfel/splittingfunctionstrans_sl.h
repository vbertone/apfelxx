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
   * @defgroup TransSF Transversely polarised splitting functions
   * @ingroup SLSplittings
   * @note Reference https://arxiv.org/pdf/hep-ph/9706511v2.pdf.
   */
  ///@{
  ///@}
  /**
   * @defgroup LOtranssf LO splitting functions
   * @ingroup TransSF
   */
  ///@{
  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB>) non-singlet transversely polarised
   * splitting function.
   */
  class P0transns: public Expression
  {
  public:
    P0transns();
    double Regular(double const&)    const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };
  ///@}

  /**
   * @defgroup NLOtranssf NLO splitting functions
   * @ingroup TransSF
   */
  ///@{
  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet-plus
   * transversely polarised splitting function.
   */
  class P1transnsp: public Expression
  {
  public:
    P1transnsp(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  protected:
    int const _nf;
    double    _a2;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet-minus
   * transversely polarised splitting function.
   */
  class P1transnsm: public P1transnsp
  {
  public:
    P1transnsm(int const& nf);
    double Regular(double const& x) const;
  };
  ///@}
}
