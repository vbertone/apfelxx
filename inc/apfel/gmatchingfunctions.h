//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//         Simone Rodini: rodini.simone.luigi@gmail.com
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
  class Cgtmd1nseUU: public Expression
  {
  public:
    Cgtmd1nseUU(double const& xi);
    double Regular(double const& x) const;
    double Local(double const&)     const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even quark-quark
   * unpolarised matching function for GTMDs.
   */
  class Cgtmd1qqeUU: public Expression
  {
  public:
    Cgtmd1qqeUU(double const& xi);
    double Regular(double const& x) const;
    double Local(double const&)     const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even quark-gluon
   * unpolarised matching function for GTMDs.
   */
  class Cgtmd1qgeUU: public Expression
  {
  public:
    Cgtmd1qgeUU(double const& xi);
    double Regular(double const& x)    const;
    double SingularPV(double const& x) const;
    double LocalLogPV(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even gluon-quark
   * unpolarised matching function for GTMDs.
   */
  class Cgtmd1gqeUU: public Expression
  {
  public:
    Cgtmd1gqeUU(double const& xi);
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even gluon-gluon
   * unpolarised matching function for GTMDs.
   */
  class Cgtmd1ggeUU: public Expression
  {
  public:
    Cgtmd1ggeUU(double const& xi);
    double Regular(double const& x)    const;
    double Local(double const&)        const;
    double SingularPV(double const& x) const;
    double LocalLogPV(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity odd quark-gluon
   * unpolarised matching function for GTMDs.
   */
  class Cgtmd1qgoUU: public Expression
  {
  public:
    Cgtmd1qgoUU(double const& xi);
    double LocalPV() const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity odd gluon-gluon
   * unpolarised matching function for GTMDs.
   */
  class Cgtmd1ggoUU: public Expression
  {
  public:
    Cgtmd1ggoUU(double const& xi);
    double LocalPV() const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even quark-gluon
   * unpolarised to transversely polarized matching function for
   * GTMDs.
   */
  class Cgtmd1qgeUT : public Expression
  {
   public:
     Cgtmd1qgeUT(double const &xi);
     double Regular(double const &x) const;
     double SingularPV(double const &x) const;
     double LocalLogPV(double const &x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity odd quark-gluon
   * unpolarised to transversely polarized matching function for
   * GTMDs.
   */
  class Cgtmd1qgoUT : public Expression
  {
   public:
     Cgtmd1qgoUT(double const &xi);
     double LocalPV() const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even gluon-gluon
   * unpolarised to transversely polarized matching function for
   * GTMDs.
   */
  class Cgtmd1ggeUT : public Expression
  {
   public:
     Cgtmd1ggeUT(double const &xi);
     double Regular(double const &x) const;
     double SingularPV(double const &x) const;
     double LocalLogPV(double const &x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity odd gluon-gluon
   * unpolarised to transversely polarized matching function for
   * GTMDs.
   */
  class Cgtmd1ggoUT : public Expression
  {
   public:
     Cgtmd1ggoUT(double const &xi);
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
  class Cgtmd1nseLL: public Expression
  {
  public:
    Cgtmd1nseLL(double const& xi);
    double Regular(double const& x) const;
    double Local(double const&)     const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even quark-quark
   * longitudinally polarised matching function for GTMDs.
   */
  class Cgtmd1qqeLL: public Expression
  {
  public:
    Cgtmd1qqeLL(double const& xi);
    double Regular(double const& x) const;
    double Local(double const&)     const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even quark-gluon
   * longitudinally polarised matching function for GTMDs.
   */
  class Cgtmd1qgeLL: public Expression
  {
  public:
    Cgtmd1qgeLL(double const& xi);
    double Regular(double const& x)    const;
    double SingularPV(double const& x) const;
    double LocalLogPV(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity odd quark-gluon
   * longitudinally polarised matching function for GTMDs.
   */
  class Cgtmd1qgoLL : public Expression
  {
   public:
     Cgtmd1qgoLL(double const &xi);
     double LocalPV() const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even gluon-quark
   * longitudinally polarised matching function for GTMDs.
   */
  class Cgtmd1gqeLL: public Expression
  {
  public:
    Cgtmd1gqeLL(double const& xi);
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even gluon-gluon
   * longitudinally polarised matching function for GTMDs.
   */
  class Cgtmd1ggeLL: public Expression
  {
  public:
    Cgtmd1ggeLL(double const& xi);
    double Regular(double const& x)    const;
    double Local(double const&)        const;
    double SingularPV(double const& x) const;
    double LocalLogPV(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity odd gluon-gluon
   * longitudinally polarised matching function for GTMDs.
   */
  class Cgtmd1ggoLL : public Expression
  {
   public:
     Cgtmd1ggoLL(double const &xi);
     double LocalPV() const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even quark-gluon
   * longitudinally polarised to transversely polarized matching
   * function for GTMDs.
   */
  class Cgtmd1qgeLT : public Expression
  {
   public:
     Cgtmd1qgeLT(double const &xi);
     double Regular(double const &x) const;
     double SingularPV(double const &x) const;
     double LocalLogPV(double const &x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity odd quark-gluon
   * longitudinally polarised to transversely polarized matching
   * function for GTMDs.
   */
  class Cgtmd1qgoLT : public Expression
  {
   public:
     Cgtmd1qgoLT(double const &xi);
     double LocalPV() const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even gluon-gluon
   * longitudinally polarised to transversely polarized matching
   * function for GTMDs.
   */
  class Cgtmd1ggeLT : public Expression
  {
   public:
     Cgtmd1ggeLT(double const &xi);
     double Regular(double const &x) const;
     double SingularPV(double const &x) const;
     double LocalLogPV(double const &x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity odd gluon-gluon
   * longitudinally polarised to transversely polarized matching
   * function for GTMDs.
   */
  class Cgtmd1ggoLT : public Expression
  {
   public:
     Cgtmd1ggoLT(double const &xi);
     double LocalPV() const;
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
  class Cgtmd1ggeTT: public Expression
  {
  public:
    Cgtmd1ggeTT(double const& xi);
    double Regular(double const& x)    const;
    double Local(double const&)        const;
    double SingularPV(double const& x) const;
    double LocalLogPV(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity odd gluon-gluon
   * linearly polarised matching function for GTMDs.
   */
  class Cgtmd1ggoTT: public Expression
  {
  public:
    Cgtmd1ggoTT(double const& xi);
    double LocalPV() const;
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
  class Cgtmd1gqeTU: public Expression
  {
  public:
    Cgtmd1gqeTU(double const& xi);
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even
   * gluon(lin)-gluon(unpol) matching function for GTMDs.
   */
  class Cgtmd1ggeTU: public Expression
  {
  public:
    Cgtmd1ggeTU(double const& xi);
    double Regular(double const& x)    const;
    double SingularPV(double const& x) const;
    double LocalLogPV(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity odd
   * gluon(lin)-gluon(unpol) matching function for GTMDs.
   */
  class Cgtmd1ggoTU: public Expression
  {
  public:
    Cgtmd1ggoTU(double const& xi);
    double LocalPV() const;;
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
  class Cgtmd1gqeTL: public Expression
  {
  public:
    Cgtmd1gqeTL(double const& xi);
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity even
   * gluon(lin)-gluon(long pol) matching function for GTMDs.
   */
  class Cgtmd1ggeTL: public Expression
  {
  public:
    Cgtmd1ggeTL(double const& xi);
    double Regular(double const& x)    const;
    double SingularPV(double const& x) const;
    double LocalLogPV(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) parity odd
   * gluon(lin)-gluon(long pol) matching function for GTMDs.
   */
  class Cgtmd1ggoTL: public Expression
  {
  public:
    Cgtmd1ggoTL(double const& xi);
    double LocalPV() const;;
  };
  ///@}
}
