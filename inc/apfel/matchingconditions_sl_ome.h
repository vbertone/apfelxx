//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/expression.h"
#include "apfel/ome.h"

namespace apfel
{

  /**
   * @defgroup MatchCond Space-like matching conditions from libome
   * @note This is a set of classes whic wrap the functions of
   * libome. The original code can be found at:
   * https://gitlab.com/libome/libome and the reference publication
   * is:
   * J. Ablinger, A. Behring, J. Blümlein, A. De Freitas, A. von
   * Manteuffel, C. Schneider, and K. Schönwald, "The Single-Mass
   * Variable Flavor Number Scheme at Three-Loop Order",
   * https://arxiv.org/pdf/2510.02175 (DESY 24-037).
   */
  ///@{
  /**
   * @defgroup NLOMCOME NLO unpolarised matching conditions from libome
   * @ingroup MatchCond
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS1Hg_L_ome: public Expression
  {
  public:
    AS1Hg_L_ome();
    double Regular(double const& x) const;
  private:
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS1ggH_L_ome: public Expression
  {
  public:
    AS1ggH_L_ome();
    double Local(double const&) const;
  private:
    ome::ome_nf_const<double> const& _ome_l;
  };
  ///@}

  /**
   * @defgroup NNLOMCOME NNLO unpolarised matching conditions from libome
   * @ingroup MatchCond
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) constant term.
   */
  class APS2Hq_0_ome: public Expression
  {
  public:
    APS2Hq_0_ome();
    double Regular(double const& x) const;
  private:
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class APS2Hq_L_ome: public Expression
  {
  public:
    APS2Hq_L_ome();
    double Regular(double const& x) const;
  private:
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln<SUP>2</SUP>(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class APS2Hq_L2_ome: public Expression
  {
  public:
    APS2Hq_L2_ome();
    double Regular(double const& x) const;
  private:
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) constant term.
   */
  class AS2Hg_0_ome: public Expression
  {
  public:
    AS2Hg_0_ome();
    double Regular(double const& x) const;
  private:
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS2Hg_L_ome: public Expression
  {
  public:
    AS2Hg_L_ome();
    double Regular(double const& x) const;
  private:
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln<SUP>2</SUP>(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS2Hg_L2_ome: public Expression
  {
  public:
    AS2Hg_L2_ome();
    double Regular(double const& x) const;
  private:
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) constant term.
   */
  class ANS2qqH_0_ome: public Expression
  {
  public:
    ANS2qqH_0_ome();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    ome::ome_nf<double>       const& _ome_r;
    ome::ome_nf_plus<double>  const& _ome_s;
    ome::ome_nf_const<double> const& _ome_l;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class ANS2qqH_L_ome: public Expression
  {
  public:
    ANS2qqH_L_ome();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    ome::ome_nf<double>       const& _ome_r;
    ome::ome_nf_plus<double>  const& _ome_s;
    ome::ome_nf_const<double> const& _ome_l;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln<SUP>2</SUP>(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class ANS2qqH_L2_ome: public Expression
  {
  public:
    ANS2qqH_L2_ome();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    ome::ome_nf<double>       const& _ome_r;
    ome::ome_nf_plus<double>  const& _ome_s;
    ome::ome_nf_const<double> const& _ome_l;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) constant term.
   */
  class AS2gqH_0_ome: public Expression
  {
  public:
    AS2gqH_0_ome();
    double Regular(double const& x) const;
  private:
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS2gqH_L_ome: public Expression
  {
  public:
    AS2gqH_L_ome();
    double Regular(double const& x) const;
  private:
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln<SUP>2</SUP>(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS2gqH_L2_ome: public Expression
  {
  public:
    AS2gqH_L2_ome();
    double Regular(double const& x) const;
  private:
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) constant term.
   */
  class AS2ggH_0_ome: public Expression
  {
  public:
    AS2ggH_0_ome();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    ome::ome_nf<double>       const& _ome_r;
    ome::ome_nf_plus<double>  const& _ome_s;
    ome::ome_nf_const<double> const& _ome_l;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS2ggH_L_ome: public Expression
  {
  public:
    AS2ggH_L_ome();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    ome::ome_nf<double>       const& _ome_r;
    ome::ome_nf_plus<double>  const& _ome_s;
    ome::ome_nf_const<double> const& _ome_l;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln<SUP>2</SUP>(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS2ggH_L2_ome: public Expression
  {
  public:
    AS2ggH_L2_ome();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    ome::ome_nf<double>       const& _ome_r;
    ome::ome_nf_plus<double>  const& _ome_s;
    ome::ome_nf_const<double> const& _ome_l;
  };
  ///@}

  /**
   * @defgroup NNNLOMCOME NNNLO unpolarised matching conditions from libome
   * @ingroup MatchCond
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) constant term.
   */
  class APS3Hq_0_ome: public Expression
  {
  public:
    APS3Hq_0_ome(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) constant term.
   */
  class AS3Hg_0_ome: public Expression
  {
  public:
    AS3Hg_0_ome(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) constant term.
   */
  class ANS3qqH_0_ome: public Expression
  {
  public:
    ANS3qqH_0_ome(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
    ome::ome_nf<double>       const& _ome_r;
    ome::ome_nf_plus<double>  const& _ome_s;
    ome::ome_nf_const<double> const& _ome_l;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) constant term.
   */
  class ANS3qqHm_0_ome: public Expression
  {
  public:
    ANS3qqHm_0_ome(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
    ome::ome_nf<double>       const& _ome_r;
    ome::ome_nf_plus<double>  const& _ome_s;
    ome::ome_nf_const<double> const& _ome_l;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) constant term.
   */
  class AS3gqH_0_ome: public Expression
  {
  public:
    AS3gqH_0_ome(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) constant term.
   */
  class AS3ggH_0_ome: public Expression
  {
  public:
    AS3ggH_0_ome(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
    ome::ome_nf<double>       const& _ome_r;
    ome::ome_nf_plus<double>  const& _ome_s;
    ome::ome_nf_const<double> const& _ome_l;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) constant term.
   */
  class AS3qgQ_0_ome: public Expression
  {
  public:
    AS3qgQ_0_ome(int const& nf);
    double Regular(double const& x)  const;
  private:
    int const _nf;
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) constant term.
   */
  class APS3qqQ_0_ome: public Expression
  {
  public:
    APS3qqQ_0_ome(int const& nf);
    double Regular(double const& x)  const;
  private:
    int const _nf;
    ome::ome_nf<double> const& _ome_r;
  };
  ///@}

  /**
   * @defgroup NLOMCpolOME NLO longitudinally polarised matching conditions from libome
   * @ingroup MatchCond
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS1polHg_L_ome: public Expression
  {
  public:
    AS1polHg_L_ome();
    double Regular(double const& x) const;
  private:
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS1polggH_L_ome: public Expression
  {
  public:
    AS1polggH_L_ome();
    double Local(double const&) const;
  private:
    ome::ome_nf_const<double> const& _ome_l;
  };
  ///@}

  /**
   * @defgroup NNLOMCpolOME NNLO longitudinally polarised matching conditions from libome
   * @ingroup MatchCond
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) constant term.
   */
  class APS2polHq_0_ome: public Expression
  {
  public:
    APS2polHq_0_ome();
    double Regular(double const& x) const;
  private:
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class APS2polHq_L_ome: public Expression
  {
  public:
    APS2polHq_L_ome();
    double Regular(double const& x) const;
  private:
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln<SUP>2</SUP>(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class APS2polHq_L2_ome: public Expression
  {
  public:
    APS2polHq_L2_ome();
    double Regular(double const& x) const;
  private:
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) constant term.
   */
  class AS2polHg_0_ome: public Expression
  {
  public:
    AS2polHg_0_ome();
    double Regular(double const& x) const;
  private:
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS2polHg_L_ome: public Expression
  {
  public:
    AS2polHg_L_ome();
    double Regular(double const& x) const;
  private:
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln<SUP>2</SUP>(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS2polHg_L2_ome: public Expression
  {
  public:
    AS2polHg_L2_ome();
    double Regular(double const& x) const;
  private:
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) constant term.
   */
  class ANS2polqqH_0_ome: public Expression
  {
  public:
    ANS2polqqH_0_ome();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    ome::ome_nf<double>       const& _ome_r;
    ome::ome_nf_plus<double>  const& _ome_s;
    ome::ome_nf_const<double> const& _ome_l;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class ANS2polqqH_L_ome: public Expression
  {
  public:
    ANS2polqqH_L_ome();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    ome::ome_nf<double>       const& _ome_r;
    ome::ome_nf_plus<double>  const& _ome_s;
    ome::ome_nf_const<double> const& _ome_l;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln<SUP>2</SUP>(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class ANS2polqqH_L2_ome: public Expression
  {
  public:
    ANS2polqqH_L2_ome();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    ome::ome_nf<double>       const& _ome_r;
    ome::ome_nf_plus<double>  const& _ome_s;
    ome::ome_nf_const<double> const& _ome_l;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) constant term.
   */
  class AS2polgqH_0_ome: public Expression
  {
  public:
    AS2polgqH_0_ome();
    double Regular(double const& x) const;
  private:
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS2polgqH_L_ome: public Expression
  {
  public:
    AS2polgqH_L_ome();
    double Regular(double const& x) const;
  private:
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln<SUP>2</SUP>(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS2polgqH_L2_ome: public Expression
  {
  public:
    AS2polgqH_L2_ome();
    double Regular(double const& x) const;
  private:
    ome::ome_nf<double> const& _ome_r;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) constant term.
   */
  class AS2polggH_0_ome: public Expression
  {
  public:
    AS2polggH_0_ome();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    ome::ome_nf<double>       const& _ome_r;
    ome::ome_nf_plus<double>  const& _ome_s;
    ome::ome_nf_const<double> const& _ome_l;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS2polggH_L_ome: public Expression
  {
  public:
    AS2polggH_L_ome();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    ome::ome_nf<double>       const& _ome_r;
    ome::ome_nf_plus<double>  const& _ome_s;
    ome::ome_nf_const<double> const& _ome_l;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) term propotional to
   * ln<SUP>2</SUP>(&mu;<SUP>2</SUP>/m<SUP>2</SUP>).
   */
  class AS2polggH_L2_ome: public Expression
  {
  public:
    AS2polggH_L2_ome();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    ome::ome_nf<double>       const& _ome_r;
    ome::ome_nf_plus<double>  const& _ome_s;
    ome::ome_nf_const<double> const& _ome_l;
  };
  ///@}
  ///@}
}
