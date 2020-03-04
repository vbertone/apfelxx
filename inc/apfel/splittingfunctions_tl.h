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
   * @defgroup TLSplittings Time-like splitting function
   * Collection of the MSbar time-like splitting functions up to the
   * highest order currently known for unpolarised, polarised (not
   * yet!), and transversity evolution.
   * @note While for the O(&alpha;<SUB>s</SUB>) and
   * O(&alpha;<SUB>s</SUB><SUP>2</SUP>) splitting functions exact
   * expressions are used, a fast parameterisation for the
   * O(&alpha;<SUB>s</SUB><SUP>3</SUP>) ones is used. See
   * https://www.liverpool.ac.uk/~avogt/split.html for more details.
   */
  ///@{
  ///@}
  /**
   * @defgroup UnpSFtl Unpolarised splitting functions
   * @ingroup TLSplittings
   */
  ///@{
  ///@}
  /**
   * @defgroup LOunpsftl LO splitting functions
   * @ingroup UnpSFtl
   */
  ///@{
  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB>) non-singlet unpolarised
   * splitting function.
   */
  class P0Tns: public Expression
  {
  public:
    P0Tns();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB>) quark-gluon unpolarised
   * splitting function.
   */
  class P0Tqg: public Expression
  {
  public:
    P0Tqg();
    double Regular(double const& x) const;
  };

  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB>) gluon-quark unpolarised
   * splitting function.
   */
  class P0Tgq: public Expression
  {
  public:
    P0Tgq();
    double Regular(double const& x) const;
  };

  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB>) gluon-gluon unpolarised
   * splitting function.
   */
  class P0Tgg: public Expression
  {
  public:
    P0Tgg(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };
  ///@}

  /**
   * @defgroup NLOunpsftl NLO splitting functions
   * @ingroup UnpSFtl
   */
  ///@{
  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>)
   * non-singlet-plus unpolarised splitting function.
   */
  class P1Tnsp: public Expression
  {
  public:
    P1Tnsp(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  protected:
    int const _nf;
    double    _a2;
  };

  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>)
   * non-singlet-minus unpolarised splitting function.
   */
  class P1Tnsm: public P1Tnsp
  {
  public:
    P1Tnsm(int const& nf);
    double Regular(double const& x) const;
  };

  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * unpolarised splitting function.
   */
  class P1Tps: public Expression
  {
  public:
    P1Tps(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>) quark-gluon
   * unpolarised splitting function.
   */
  class P1Tqg: public Expression
  {
  public:
    P1Tqg(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon-quark
   * unpolarised splitting function.
   */
  class P1Tgq: public Expression
  {
  public:
    P1Tgq(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon-gluon
   * unpolarised splitting function.
   */
  class P1Tgg: public Expression
  {
  public:
    P1Tgg(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
    double    _a2g;
  };
  ///@}

  /**
   * @defgroup NNLOunpsftl NNLO splitting functions
   * @ingroup UnpSFtl
   */
  ///@{
  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB><SUP>3</SUP>)
   * non-singlet-plus unpolarised splitting function.
   */
  class P2Tnsp: public Expression
  {
  public:
    P2Tnsp(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };

  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB><SUP>3</SUP>)
   * non-singlet-minus unpolarised splitting function.
   */
  class P2Tnsm: public Expression
  {
  public:
    P2Tnsm(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };

  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB><SUP>3</SUP>)
   * non-singlet-valence unpolarised splitting function minus
   * non-singlet-minus unpolarised splitting function.
   */
  class P2Tnss: public Expression
  {
  public:
    P2Tnss(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB><SUP>3</SUP>) pure-singlet
   * unpolarised splitting function.
   */
  class P2Tps: public Expression
  {
  public:
    P2Tps(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB><SUP>3</SUP>) quark-gluon
   * unpolarised splitting function.
   */
  class P2Tqg: public Expression
  {
  public:
    P2Tqg(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB><SUP>3</SUP>) gluon-quark
   * unpolarised splitting function.
   */
  class P2Tgq: public Expression
  {
  public:
    P2Tgq(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Time-like O(&alpha;<SUB>s</SUB><SUP>3</SUP>) gluon-gluon
   * unpolarised splitting function.
   */
  class P2Tgg: public Expression
  {
  public:
    P2Tgg(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };
  ///@}

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
