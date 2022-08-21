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
   * @defgroup NCMasslesspol Zero-mass coefficient functions for longitudinally polarised DIS
   * Collection of the Zero-mass coefficient functions up to
   * O(&alpha;<SUB>s</SUB>) for g<SUB>4</SUB> and g<SUB>L</SUB>, and
   * up to O(&alpha;<SUB>s</SUB><SUP>2</SUP>) for g<SUB>1</SUB>.
   * @note The expression for g1 are taken from:
   * https://inspirehep.net/literature/353973.
   */
  ///@{
  /**
   * @defgroup NLOzmpol NLO zero-mass coefficient functions
   * @ingroup NCMasslesspol
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet coefficient function
   * for g<SUB>4</SUB>.
   */
  class G41ns: public Expression
  {
  public:
    G41ns();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet coefficient function
   * for g<SUB>L</SUB>.
   */
  class GL1ns: public Expression
  {
  public:
    GL1ns();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet coefficient function
   * for g<SUB>1</SUB>.
   */
  class G11ns: public Expression
  {
  public:
    G11ns();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon coefficient function for for
   * g<SUB>1</SUB>.
   */
  class G11g: public Expression
  {
  public:
    G11g();
    double Regular(double const& x) const;
  };
  ///@}

  /**
   * @defgroup NNLOzmpol NNLO zero-mass coefficient functions
   * @ingroup NCMasslesspol
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet-plus
   * coefficient function for F2.
   */
  class G12nsp: public Expression
  {
  public:
    G12nsp(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * coefficient function for F2.
   */
  class G12ps: public Expression
  {
  public:
    G12ps();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon coefficient
   * function for F2.
   */
  class G12g: public Expression
  {
  public:
    G12g();
    double Regular(double const& x) const;
  };
  ///@}
  ///@}
}
