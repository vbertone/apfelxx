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
   * @brief The O(&alpha;<SUB>s</SUB>) parity even non-singlet
   * unpolarised matching function for GTMDs.
   */
  class Cgtmd1nse: public Expression
  {
  public:
    Cgtmd1nse(double const& xi);
    double Regular(double const& x) const;
    double Local(double const&)     const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even quark-quark
   * unpolarised matching function for GTMDs.
   */
  class Cgtmd1qqe: public Expression
  {
  public:
    Cgtmd1qqe(double const& xi);
    double Regular(double const& x) const;
    double Local(double const&)     const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even quark-gluon
   * unpolarised matching function for GTMDs.
   */
  class Cgtmd1qge: public Expression
  {
  public:
    Cgtmd1qge(double const& xi);
    double Regular(double const& x)    const;
    double SingularPV(double const& x) const;
    double LocalLogPV(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even gluon-quark
   * unpolarised matching function for GTMDs.
   */
  class Cgtmd1gqe: public Expression
  {
  public:
    Cgtmd1gqe(double const& xi);
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even gluon-gluon
   * unpolarised matching function for GTMDs.
   */
  class Cgtmd1gge: public Expression
  {
  public:
    Cgtmd1gge(double const& xi);
    double Regular(double const& x)    const;
    double Local(double const&)        const;
    double SingularPV(double const& x) const;
    double LocalLogPV(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity odd quark-gluon
   * unpolarised matching function for GTMDs.
   */
  class Cgtmd1qgo: public Expression
  {
  public:
    Cgtmd1qgo(double const& xi);
    double LocalPV() const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity odd gluon-gluon
   * unpolarised matching function for GTMDs.
   */
  class Cgtmd1ggo: public Expression
  {
  public:
    Cgtmd1ggo(double const& xi);
    double LocalPV() const;
  };
  ///@}

  /**
   * @defgroup NLOmatchPol NLO longitudinally polarised matching functions for GTMDs
   * @ingroup GTMDMatchingFunctions
   */
  ///@{
  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even non-singlet
   * longitudinally polarised matching function for GTMDs.
   */
  class Cgtmd1polnse: public Expression
  {
  public:
    Cgtmd1polnse(double const& xi);
    double Regular(double const& x) const;
    double Local(double const&)     const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even quark-quark
   * longitudinally polarised matching function for GTMDs.
   */
  class Cgtmd1polqqe: public Expression
  {
  public:
    Cgtmd1polqqe(double const& xi);
    double Regular(double const& x) const;
    double Local(double const&)     const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even quark-gluon
   * longitudinally polarised matching function for GTMDs.
   */
  class Cgtmd1polqge: public Expression
  {
  public:
    Cgtmd1polqge(double const& xi);
    double Regular(double const& x)    const;
    double SingularPV(double const& x) const;
    double LocalLogPV(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even gluon-quark
   * longitudinally polarised matching function for GTMDs.
   */
  class Cgtmd1polgqe: public Expression
  {
  public:
    Cgtmd1polgqe(double const& xi);
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even gluon-gluon
   * longitudinally polarised matching function for GTMDs.
   */
  class Cgtmd1polgge: public Expression
  {
  public:
    Cgtmd1polgge(double const& xi);
    double Regular(double const& x)    const;
    double Local(double const&)        const;
    double SingularPV(double const& x) const;
    double LocalLogPV(double const& x) const;
  };
  ///@}

  /**
   * @defgroup NLOmatchTrans NLO transversely/linearly polarised
   * matching functions for GTMDs
   * @ingroup GTMDMatchingFunctions
   */
  ///@{
  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even gluon-gluon
   * linearly polarised matching function for GTMDs.
   */
  class Cgtmd1lingge: public Expression
  {
  public:
    Cgtmd1lingge(double const& xi);
    double Regular(double const& x)    const;
    double Local(double const&)        const;
    double SingularPV(double const& x) const;
    double LocalLogPV(double const& x) const;
  };
  ///@}

  /**
   * @defgroup NLOmatchLinUnp NLO matching functions for GTMDs for linealy-polarised gluons onto unpolarised partons
   * @ingroup GTMDMatchingFunctions
   */
  ///@{
  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even
   * gluon(lin)-quark(unpol) matching function for GTMDs.
   */
  class Cgtmd1linunpgqe: public Expression
  {
  public:
    Cgtmd1linunpgqe(double const& xi);
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even
   * gluon(lin)-gluon(unpol) matching function for GTMDs.
   */
  class Cgtmd1linunpgge: public Expression
  {
  public:
    Cgtmd1linunpgge(double const& xi);
    double Regular(double const& x)    const;
    double Local(double const&)        const;
    double SingularPV(double const& x) const;
    double LocalLogPV(double const& x) const;
  };
  ///@}

  /**
   * @defgroup NLOmatchLinPol NLO matching functions for GTMDs for linealy-polarised gluons onto longitudinally polarised partons
   * @ingroup GTMDMatchingFunctions
   */
  ///@{
  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even
   * gluon(lin)-quark(long pol) matching function for GTMDs.
   */
  class Cgtmd1linpolgqe: public Expression
  {
  public:
    Cgtmd1linpolgqe(double const& xi);
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even
   * gluon(lin)-gluon(long pol) matching function for GTMDs.
   */
  class Cgtmd1linpolgge: public Expression
  {
  public:
    Cgtmd1linpolgge(double const& xi);
    double Regular(double const& x)    const;
    double Local(double const&)        const;
    double SingularPV(double const& x) const;
    double LocalLogPV(double const& x) const;
  };
  ///@}
}
