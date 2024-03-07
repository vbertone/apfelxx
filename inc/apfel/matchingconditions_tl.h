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
   * @defgroup MatchCondTL Time-like matching conditions
   * @note The expressions are taken from:
   * https://arxiv.org/pdf/hep-ph/0504192.pdf.
   */
  ///@{
  ///@}
  /**
   * @defgroup NLOMC NLO unpolarised matching conditions
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

  /**
   * @brief O(&alpha;<SUB>s</SUB>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>) for the HH matching. This is the QCD adaptation of
   * Eq. (4.121) of https://arxiv.org/pdf/1909.03886.pdf.
   */
  class ATS1HH_L: public Expression
  {
  public:
    ATS1HH_L();
    double Singular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) constant term for the HH
   * matching. This is the QCD adaptation of Eq. (4.121) of
   * https://arxiv.org/pdf/1909.03886.pdf.
   */
  class ATS1HH_0: public Expression
  {
  public:
    ATS1HH_0();
    double Singular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>) of eq. (B.2) of
   * https://arxiv.org/pdf/hep-ph/9612398.pdf.
   */
  class ATS1gH_L: public Expression
  {
  public:
    ATS1gH_L();
    double Regular(double const& x) const;
  };
  ///@}
}
