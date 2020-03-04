//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/expression.h"

namespace apfel
{
  /**
   * @defgroup MatchCondTL Time-like matching conditions
   * @notes The expressions are taken from:
   * https://arxiv.org/pdf/hep-ph/0504192.pdf.
   */
  ///@{
  ///@}
  /**
   * @defgroup NLOMC NLO matching conditions
   * @ingroup MatchCondTL
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB>) constant term of eq. (15) of
   * https://arxiv.org/pdf/hep-ph/0504192.pdf.
   */
  class ATS1Hg_0: public Expression
  {
  public:
    ATS1Hg_0();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>) of eq. (15) of
   * https://arxiv.org/pdf/hep-ph/0504192.pdf.
   */
  class ATS1Hg_L: public Expression
  {
  public:
    ATS1Hg_L();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>) of eq. (22) of
   * https://arxiv.org/pdf/hep-ph/0504192.pdf.
   */
  class ATS1ggH_L: public Expression
  {
  public:
    ATS1ggH_L();
    double Local(double const& x) const;
  };
  ///@}
}
