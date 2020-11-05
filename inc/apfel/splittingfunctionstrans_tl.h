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
   * @defgroup TransSFtl Transversely polarised splitting functions
   * @ingroup TLSplittings
   * @note Reference https://arxiv.org/pdf/hep-ph/0108241v1.pdf.
   */
  ///@{
  ///@}
  /**
   * @defgroup LOtranssf LO splitting functions
   * @ingroup TransSFtl
   */
  ///@{
  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB>) non-singlet transversely
   * polarised splitting function.
   */
  class P0Ttransns: public Expression
  {
  public:
    P0Ttransns();
    double Regular(double const&)    const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };
  ///@}

  /**
   * @defgroup NLOtranssf NLO splitting functions
   * @ingroup TransSFtl
   */
  ///@{
  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>)
   * non-singlet-plus transversely polarised splitting function.
   */
  class P1Ttransnsp: public Expression
  {
  public:
    P1Ttransnsp(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  protected:
    int const _nf;
    double    _a2;
  };

  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>)
   * non-singlet-minus transversely polarised splitting function.
   */
  class P1Ttransnsm: public P1Ttransnsp
  {
  public:
    P1Ttransnsm(int const& nf);
    double Regular(double const& x) const;
  };
  ///@}
}
