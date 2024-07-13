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
   * @defgroup GTMDMatchingFunctions GTMD matching functions
   * The perturbative matching functions for GTMDs.
   */
  ///@{
  ///@}
  /**
   * @defgroup NLOmatchUnp NLO unpolarised matching functions for GTMDs
   * @ingroup GTMDMatchingFunctions
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

  /**
   * @defgroup NLOmatchPol NLO longitudinally polarised matching functions for GTMDs
   * @ingroup GTMDMatchingFunctions
   */
  ///@{
  /**
   * @brief The O(&alpha;<SUB>s</SUB>) non-singlet longitudinally
   * polarised matching function for GTMDs.
   */
  class Cgtmd1polns: public Expression
  {
  public:
    Cgtmd1polns(double const& xi);
    double Regular(double const& x) const;
    double Local(double const&)     const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) quark-quark longitudinally
   * polarised matching function for GTMDs.
   */
  class Cgtmd1polqq: public Expression
  {
  public:
    Cgtmd1polqq(double const& xi);
    double Regular(double const& x) const;
    double Local(double const&)     const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) quark-gluon longitudinally
   * polarised matching function for GTMDs.
   */
  class Cgtmd1polqg: public Expression
  {
  public:
    Cgtmd1polqg(double const& xi);
    double Regular(double const& x)    const;
    double SingularPV(double const& x) const;
    double LocalPV(double const& x)    const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) gluon-quark longitudinally
   * polarised matching function for GTMDs.
   */
  class Cgtmd1polgq: public Expression
  {
  public:
    Cgtmd1polgq(double const& xi);
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) gluon-gluon longitudinally
   * polarised matching function for GTMDs.
   */
  class Cgtmd1polgg: public Expression
  {
  public:
    Cgtmd1polgg(double const& xi);
    double Regular(double const& x)    const;
    double Local(double const&)        const;
    double SingularPV(double const& x) const;
    double LocalPV(double const& x)    const;
  };
  ///@}

  /**
   * @defgroup NLOmatchTrans NLO transversely/linearly polarised
   * matching functions for GTMDs
   * @ingroup GTMDMatchingFunctions
   */
  ///@{
  /**
   * @brief The O(&alpha;<SUB>s</SUB>) gluon-gluon linearly polarised
   * matching function for GTMDs.
   */
  class Cgtmd1lingg: public Expression
  {
  public:
    Cgtmd1lingg(double const& xi);
    double Regular(double const& x)    const;
    double Local(double const&)        const;
    double SingularPV(double const& x) const;
    double LocalPV(double const& x)    const;
  };
  ///@}

  /**
   * @defgroup NLOmatchLinUnp NLO matching functions for GTMDs for linealy-polarised gluons onto unpolarised partons
   * @ingroup GTMDMatchingFunctions
   */
  ///@{
  /**
   * @brief The O(&alpha;<SUB>s</SUB>) gluon(lin)-quark(unpol)
   * matching function for GTMDs.
   */
  class Cgtmd1linunpgq: public Expression
  {
  public:
    Cgtmd1linunpgq(double const& xi);
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) gluon(lin)-gluon(unpol)
   * matching function for GTMDs.
   */
  class Cgtmd1linunpgg: public Expression
  {
  public:
    Cgtmd1linunpgg(double const& xi);
    double Regular(double const& x)    const;
    double Local(double const&)        const;
    double SingularPV(double const& x) const;
    double LocalPV(double const& x)    const;
  };
  ///@}

  /**
   * @defgroup NLOmatchLinPol NLO matching functions for GTMDs for linealy-polarised gluons onto longitudinally polarised partons
   * @ingroup GTMDMatchingFunctions
   */
  ///@{
  /**
   * @brief The O(&alpha;<SUB>s</SUB>) gluon(lin)-quark(long pol)
   * matching function for GTMDs.
   */
  class Cgtmd1linpolgq: public Expression
  {
  public:
    Cgtmd1linpolgq(double const& xi);
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) gluon(lin)-gluon(long pol)
   * matching function for GTMDs.
   */
  class Cgtmd1linpolgg: public Expression
  {
  public:
    Cgtmd1linpolgg(double const& xi);
    double Regular(double const& x)    const;
    double Local(double const&)        const;
    double SingularPV(double const& x) const;
    double LocalPV(double const& x)    const;
  };
  ///@}
}
