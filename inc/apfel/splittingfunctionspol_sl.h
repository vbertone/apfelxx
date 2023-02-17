//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/expression.h"
#include "apfel/splittingfunctionsunp_sl.h"

namespace apfel
{
  /**
   * @defgroup PolSF Longitudinally polarised splitting functions
   * @ingroup SLSplittings
   * @note Reference: https://arxiv.org/pdf/hep-ph/9603366.pdf
   */
  ///@{
  ///@}
  /**
   * @defgroup LOpolsf LO splitting functions
   * @ingroup PolSF
   */
  ///@{
  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB>) non-singlet
   * longitudinally polarised splitting function. This is equal to the
   * non-singlet unpolarised splitting function.
   */
  class P0polns: public P0ns
  {
  public:
    P0polns();
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB>) quark-gluon
   * longitudinally polarised splitting function.
   */
  class P0polqg: public Expression
  {
  public:
    P0polqg(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB>) gluon-quark
   * longitudinally polarised splitting function.
   */
  class P0polgq: public Expression
  {
  public:
    P0polgq();
    double Regular(double const& x) const;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB>) gluon-gluon
   * longitudinally polarised splitting function.
   */
  class P0polgg: public Expression
  {
  public:
    P0polgg(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };
  ///@}

  /**
   * @defgroup NLOpolsf NLO splitting functions
   * @ingroup PolSF
   */
  ///@{
  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>)
   * non-singlet-plus longitudinally polarised splitting
   * function. This is equal to the non-singlet-minus unpolarised
   * splitting function.
   */
  class P1polnsp: public P1nsm
  {
  public:
    P1polnsp(int const& nf);
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>)
   * non-singlet-minus longitudinally polarised splitting
   * function. This is equal to the non-singlet-plus unpolarised
   * splitting function.
   */
  class P1polnsm: public P1nsp
  {
  public:
    P1polnsm(int const& nf);
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * longitudinally polarised splitting function.
   */
  class P1polps: public Expression
  {
  public:
    P1polps(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>) quark-gluon
   * longitudinally polarised splitting function.
   */
  class P1polqg: public Expression
  {
  public:
    P1polqg(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon-quark
   * longitudinally polarised splitting function.
   */
  class P1polgq: public Expression
  {
  public:
    P1polgq(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon-gluon
   * longitudinally polarised splitting function.
   */
  class P1polgg: public Expression
  {
  public:
    P1polgg(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
    double    _a2g;
  };
  ///@}

  /**
   * @defgroup NNLOpolsf NNLO splitting functions
   * @ingroup PolSF
   */
  ///@{
  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>3</SUP>)
   * non-singlet-plus longitudinally polarised splitting
   * function. This is equal to the non-singlet-minus unpolarised
   * splitting function.
   */
  class P2polnsp: public P2nsm
  {
  public:
    P2polnsp(int const& nf);
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>3</SUP>)
   * non-singlet-minus longitudinally polarised splitting
   * function. This is equal to the non-singlet-plus unpolarised
   * splitting function.
   */
  class P2polnsm: public P2nsp
  {
  public:
    P2polnsm(int const& nf);
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>3</SUP>)
   * non-singlet-valence longitudinally polarised splitting function
   * minus non-singlet-minus longitudinally polarised splitting
   * function.
   * @note This has been computed in https://arxiv.org/pdf/1506.04517.pdf.
   */
  class P2polnss: public Expression
  {
  public:
    P2polnss(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>3</SUP>) pure-singlet
   * longitudinally polarised splitting function.
   */
  class P2polps: public Expression
  {
  public:
    P2polps(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>3</SUP>) quark-gluon
   * longitudinally polarised splitting function.
   */
  class P2polqg: public Expression
  {
  public:
    P2polqg(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>3</SUP>) gluon-quark
   * longitudinally polarised splitting function.
   */
  class P2polgq: public Expression
  {
  public:
    P2polgq(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>3</SUP>) gluon-gluon
   * longitudinally polarised splitting function.
   */
  class P2polgg: public Expression
  {
  public:
    P2polgg(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };
  ///@}
}
