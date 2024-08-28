//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/expression.h"

#include <vector>

namespace apfel
{
  /**
   * @defgroup MatchCond Space-like matching conditions
   * @note The expressions are taken from:
   * https://arxiv.org/pdf/hep-ph/9612398.pdf. Note that in these
   * expressions ln(m<SUP>2</SUP>/&mu;<SUP>2</SUP>) appears while we
   * need ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>), so we need to include a
   * minus sign in front of every term linear in this log.
   */
  ///@{
  /**
   * @defgroup NLOMC NLO unpolarised matching conditions
   * @ingroup MatchCond
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>) of eq. (B.2) of
   * https://arxiv.org/pdf/hep-ph/9612398.pdf.
   */
  class AS1Hg_L: public Expression
  {
  public:
    AS1Hg_L();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>) of eq (B.6) of
   * https://arxiv.org/pdf/hep-ph/9612398.pdf.
   */
  class AS1ggH_L: public Expression
  {
  public:
    AS1ggH_L();
    double Local(double const&) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>) for the HH matching. This is the QCD adaptation of
   * Eq. (4.121) of https://arxiv.org/pdf/1909.03886.pdf.
   */
  class AS1HH_L: public Expression
  {
  public:
    AS1HH_L();
    double Singular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) constant term for the HH
   * matching. This is the QCD adaptation of Eq. (4.121) of
   * https://arxiv.org/pdf/1909.03886.pdf.
   */
  class AS1HH_0: public Expression
  {
  public:
    AS1HH_0();
    double Singular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>) for the gH matching. This is
   * the QCD adaptation of Eq. (4.189) of
   * https://arxiv.org/pdf/1909.03886.pdf.
   */
  class AS1gH_L: public Expression
  {
  public:
    AS1gH_L();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) constant term for the gH
   * matching. This is the QCD adaptation of Eq. (4.189) of
   * https://arxiv.org/pdf/1909.03886.pdf.
   */
  class AS1gH_0: public Expression
  {
  public:
    AS1gH_0();
    double Regular(double const& x) const;
  };
  ///@}

  /**
   * @defgroup NNLOMC NNLO unpolarised matching conditions
   * @ingroup MatchCond
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) constant term of eq
   * (B.1) of https://arxiv.org/pdf/hep-ph/9612398.pdf.
   */
  class APS2Hq_0: public Expression
  {
  public:
    APS2Hq_0();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>) of eq (B.1) of
   * https://arxiv.org/pdf/hep-ph/9612398.pdf.
   */
  class APS2Hq_L: public Expression
  {
  public:
    APS2Hq_L();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln<SUP>2</SUP>(&mu;<SUP>2</SUP>/m<SUP>2</SUP>) of eq (B.1) of
   * https://arxiv.org/pdf/hep-ph/9612398.pdf.
   */
  class APS2Hq_L2: public Expression
  {
  public:
    APS2Hq_L2();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) constant term of eq
   * (B.3) of https://arxiv.org/pdf/hep-ph/9612398.pdf.
   */
  class AS2Hg_0: public Expression
  {
  public:
    AS2Hg_0();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>) of eq (B.3) of
   * https://arxiv.org/pdf/hep-ph/9612398.pdf.
   */
  class AS2Hg_L: public Expression
  {
  public:
    AS2Hg_L();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln<SUP>2</SUP>(&mu;<SUP>2</SUP>/m<SUP>2</SUP>) of eq (B.3) of
   * https://arxiv.org/pdf/hep-ph/9612398.pdf.
   */
  class AS2Hg_L2: public Expression
  {
  public:
    AS2Hg_L2();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) constant term of eq
   * (B.4) of https://arxiv.org/pdf/hep-ph/9612398.pdf.
   */
  class ANS2qqH_0: public Expression
  {
  public:
    ANS2qqH_0();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>) of eq (B.4) of
   * https://arxiv.org/pdf/hep-ph/9612398.pdf.
   */
  class ANS2qqH_L: public Expression
  {
  public:
    ANS2qqH_L();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln<SUP>2</SUP>(&mu;<SUP>2</SUP>/m<SUP>2</SUP>) of eq (B.4) of
   * https://arxiv.org/pdf/hep-ph/9612398.pdf.
   */
  class ANS2qqH_L2: public Expression
  {
  public:
    ANS2qqH_L2();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) constant term of eq
   * (B.5) of https://arxiv.org/pdf/hep-ph/9612398.pdf.
   */
  class AS2gqH_0: public Expression
  {
  public:
    AS2gqH_0();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>) of eq (B.5) of
   * https://arxiv.org/pdf/hep-ph/9612398.pdf.
   */
  class AS2gqH_L: public Expression
  {
  public:
    AS2gqH_L();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln<SUP>2</SUP>(&mu;<SUP>2</SUP>/m<SUP>2</SUP>) of eq (B.5) of
   * https://arxiv.org/pdf/hep-ph/9612398.pdf.
   */
  class AS2gqH_L2: public Expression
  {
  public:
    AS2gqH_L2();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) constant term of eq
   * (B.7) of https://arxiv.org/pdf/hep-ph/9612398.pdf.
   */
  class AS2ggH_0: public Expression
  {
  public:
    AS2ggH_0();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>) of eq (B.7) of
   * https://arxiv.org/pdf/hep-ph/9612398.pdf.
   */
  class AS2ggH_L: public Expression
  {
  public:
    AS2ggH_L();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln<SUP>2</SUP>(&mu;<SUP>2</SUP>/m<SUP>2</SUP>) of eq (B.7) of
   * https://arxiv.org/pdf/hep-ph/9612398.pdf.
   */
  class AS2ggH_L2: public Expression
  {
  public:
    AS2ggH_L2();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };
  ///@}

  /**
   * @defgroup NNNLOMC NNNLO unpolarised matching conditions
   * @note Approximated expressions from
   * https://github.com/MSHTPDF/N3LO_additions. Details to be found in
   * https://arxiv.org/pdf/2207.04739.pdf. Logarithmic terms currently
   * unavailable.
   * @ingroup MatchCond
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) constant term.
   */
  class APS3Hq_0: public Expression
  {
  public:
    APS3Hq_0();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) constant term.
   */
  class AS3Hg_0: public Expression
  {
  public:
    AS3Hg_0(double const& rho = 12214.000);
    double Regular(double const& x) const;
  private:
    double        const _rho;
    std::vector<double> _C;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) constant term.
   */
  class ANS3qqH_0: public Expression
  {
  public:
    ANS3qqH_0(double const& rho = -64.411);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    double        const _rho;
    std::vector<double> _C;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) constant term.
   */
  class AS3gqH_0: public Expression
  {
  public:
    AS3gqH_0();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) constant term.
   */
  class AS3ggH_0: public Expression
  {
  public:
    AS3ggH_0(double const& rho = -1951.600);
    double Regular(double const& x)  const;
  private:
    double const        _rho;
    std::vector<double> _C;
  };
  ///@}

  /**
   * @defgroup NLOMCpol NLO longitudinally polarised matching conditions
   * @ingroup MatchCond
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS1polHg_L: public Expression
  {
  public:
    AS1polHg_L();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS1polggH_L: public Expression
  {
  public:
    AS1polggH_L();
    double Local(double const&) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>) for the HH matching.
   */
  class AS1polHH_L: public Expression
  {
  public:
    AS1polHH_L();
    double Singular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) constant term for the HH
   * matching.
   */
  class AS1polHH_0: public Expression
  {
  public:
    AS1polHH_0();
    double Singular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>) for the gH matching.
   */
  class AS1polgH_L: public Expression
  {
  public:
    AS1polgH_L();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) constant term for the gH
   * matching.
   */
  class AS1polgH_0: public Expression
  {
  public:
    AS1polgH_0();
    double Regular(double const& x) const;
  };
  ///@}

  /**
   * @defgroup NNLOMCpol NNLO longitudinally polarised matching conditions
   * @ingroup MatchCond
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) constant term.
   */
  class APS2polHq_0: public Expression
  {
  public:
    APS2polHq_0();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class APS2polHq_L: public Expression
  {
  public:
    APS2polHq_L();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln<SUP>2</SUP>(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class APS2polHq_L2: public Expression
  {
  public:
    APS2polHq_L2();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) constant term.
   */
  class AS2polHg_0: public Expression
  {
  public:
    AS2polHg_0();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS2polHg_L: public Expression
  {
  public:
    AS2polHg_L();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln<SUP>2</SUP>(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS2polHg_L2: public Expression
  {
  public:
    AS2polHg_L2();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) constant term.
   */
  class ANS2polqqH_0: public Expression
  {
  public:
    ANS2polqqH_0();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class ANS2polqqH_L: public Expression
  {
  public:
    ANS2polqqH_L();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln<SUP>2</SUP>(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class ANS2polqqH_L2: public Expression
  {
  public:
    ANS2polqqH_L2();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) constant term.
   */
  class AS2polgqH_0: public Expression
  {
  public:
    AS2polgqH_0();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS2polgqH_L: public Expression
  {
  public:
    AS2polgqH_L();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln<SUP>2</SUP>(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS2polgqH_L2: public Expression
  {
  public:
    AS2polgqH_L2();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) constant term.
   */
  class AS2polggH_0: public Expression
  {
  public:
    AS2polggH_0();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS2polggH_L: public Expression
  {
  public:
    AS2polggH_L();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln<SUP>2</SUP>(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS2polggH_L2: public Expression
  {
  public:
    AS2polggH_L2();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };
  ///@}
  ///@}
}
