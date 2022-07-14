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
   * @defgroup GTMDMatchingFunctionsUnp Unpolarised GTMD matching functions
   * The perturbative matching functions for unpolarised GTMDs.
   */
  ///@{
  ///@}
  /**
   * @defgroup NLOmatchUnp NLO matching functions for GTMDs
   * @ingroup GTMDMatchingFunctionsUnp
   */
  ///@{
  /**
   * @brief The O(&alpha;<SUB>s</SUB>) non-singlet unpolarised
   * matching function for GTMDs.
   */
  class Cgtmd1ns: public Expression
  {
  public:
    Cgtmd1ns(double const& xi);
    double Regular(double const& x) const;
    double Local(double const&)     const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) quark-quark unpolarised
   * matching function for GTMDs.
   */
  class Cgtmd1qq: public Expression
  {
  public:
    Cgtmd1qq(double const& xi);
    double Regular(double const& x) const;
    double Local(double const&)     const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) quark-gluon unpolarised
   * matching function for GTMDs.
   */
  class Cgtmd1qg: public Expression
  {
  public:
    Cgtmd1qg(double const& xi);
    double Regular(double const& x)    const;
    double SingularPV(double const& x) const;
    double LocalPV(double const& x)    const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) gluon-quark unpolarised
   * matching function for GTMDs.
   */
  class Cgtmd1gq: public Expression
  {
  public:
    Cgtmd1gq(double const& xi);
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) gluon-gluon unpolarised
   * matching function for GTMDs.
   */
  class Cgtmd1gg: public Expression
  {
  public:
    Cgtmd1gg(double const& xi);
    double Regular(double const& x)    const;
    double Local(double const&)        const;
    double SingularPV(double const& x) const;
    double LocalPV(double const& x)    const;
  };
  ///@}
}
