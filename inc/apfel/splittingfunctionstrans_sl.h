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
   * @note References: https://arxiv.org/pdf/hep-ph/9706511v2.pdf,
   * https://lib-extopc.kek.jp/preprints/PDF/2000/0032/0032201.pdf,
   * https://arxiv.org/pdf/hep-ph/9805295v1.pdf
   *
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

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB>) gluon-gluon linearly
   * polarised splitting function.
   * @note
   */
  class P0transgg: public Expression
  {
  public:
    P0transgg(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
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

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon-gluon linearly
   * polarised splitting function.
   * @note
   */
  class P1transgg: public Expression
  {
  public:
    P1transgg(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
    double    _a2;
  };
  ///@}
}
