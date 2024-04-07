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
   * @defgroup SLSplittings Space-like splitting functions
   * Collection of the MSbar space-like splitting functions up to the
   * highest order currently known for unpolarised, polarised, and
   * transversity evolution.

   * @note While for the O(&alpha;<SUB>s</SUB>) and
   * O(&alpha;<SUB>s</SUB><SUP>2</SUP>) splitting functions exact
   * expressions are used, a fast parameterisation for the
   * O(&alpha;<SUB>s</SUB><SUP>3</SUP>) (and
   * O(&alpha;<SUB>s</SUB><SUP>4</SUP>) when available) ones is
   * used. See https://www.liverpool.ac.uk/~avogt/split.html for more
   * details. Approximated O(&alpha;<SUB>s</SUB><SUP>4</SUP>)
   * splitting functions are taken from:
   * https://github.com/MSHTPDF/N3LO_additions and best fit parameters
   * taken from Table 8 of https://arxiv.org/pdf/2207.04739.pdf.
   */
  ///@{
  ///@}
  /**
   * @defgroup UnpSF Unpolarised splitting functions
   * @ingroup SLSplittings
   */
  ///@{
  ///@}
  /**
   * @defgroup LOunpsf LO splitting functions
   * @ingroup UnpSF
   */
  ///@{
  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB>) non-singlet unpolarised splitting
   * function.
   */
  class P0ns: public Expression
  {
  public:
    P0ns();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB>) quark-gluon unpolarised splitting
   * function.
   */
  class P0qg: public Expression
  {
  public:
    P0qg(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB>) gluon-quark unpolarised splitting
   * function.
   */
  class P0gq: public Expression
  {
  public:
    P0gq();
    double Regular(double const& x) const;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB>) gluon-gluon unpolarised splitting
   * function.
   */
  class P0gg: public Expression
  {
  public:
    P0gg(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };
  ///@}

  /**
   * @defgroup NLOunpsf NLO splitting functions
   * @ingroup UnpSF
   */
  ///@{
  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet-plus
   * unpolarised splitting function.
   */
  class P1nsp: public Expression
  {
  public:
    P1nsp(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  protected:
    int const _nf;
    double    _a2;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet-minus
   * unpolarised splitting function.
   */
  class P1nsm: public P1nsp
  {
  public:
    P1nsm(int const& nf);
    double Regular(double const& x) const;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * unpolarised splitting function.
   */
  class P1ps: public Expression
  {
  public:
    P1ps(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>) quark-gluon unpolarised
   * splitting function.
   */
  class P1qg: public Expression
  {
  public:
    P1qg(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon-quark unpolarised
   * splitting function.
   */
  class P1gq: public Expression
  {
  public:
    P1gq(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon-gluon unpolarised
   * splitting function.
   */
  class P1gg: public Expression
  {
  public:
    P1gg(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
    double    _a2g;
  };
  ///@}

  /**
   * @defgroup NNLOunpsf NNLO splitting functions
   * @ingroup UnpSF
   */
  ///@{
  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>3</SUP>) non-singlet-plus
   * unpolarised splitting function.
   */
  class P2nsp: public Expression
  {
  public:
    P2nsp(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>3</SUP>) non-singlet-minus
   * unpolarised splitting function.
   */
  class P2nsm: public Expression
  {
  public:
    P2nsm(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>3</SUP>) non-singlet-valence
   * unpolarised splitting function minus non-singlet-minus
   * unpolarised splitting function.
   */
  class P2nss: public Expression
  {
  public:
    P2nss(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>3</SUP>) pure-singlet
   * unpolarised splitting function.
   */
  class P2ps: public Expression
  {
  public:
    P2ps(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>3</SUP>) quark-gluon unpolarised
   * splitting function.
   */
  class P2qg: public Expression
  {
  public:
    P2qg(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>3</SUP>) gluon-quark unpolarised
   * splitting function.
   */
  class P2gq: public Expression
  {
  public:
    P2gq(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>3</SUP>) gluon-gluon unpolarised
   * splitting function.
   */
  class P2gg: public Expression
  {
  public:
    P2gg(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };
  ///@}

  /**
   * @defgroup NNNLOunpsf NNNLO splitting functions
   * @ingroup UnpSF
   * @note For now only leading-color plus, minus, and valence
   * contributions have been computed and parameterised. The singlet
   * ones are also parameterised and using the first Mellin moments
   * and the small- and large-x asymptotic behaviours.
   */
  ///@{
  /* /\** */
  /*  * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>4</SUP>) non-singlet-plus */
  /*  * unpolarised splitting function. */
  /*  *\/ */
  /* class P3nsp: public Expression */
  /* { */
  /* public: */
  /*   P3nsp(int const& nf, int const& imod = 0, double const& rho = 0.007); */
  /*   double Regular(double const& x)  const; */
  /*   double Singular(double const& x) const; */
  /*   double Local(double const& x)    const; */
  /* private: */
  /*   int           const _nf; */
  /*   int           const _imod; */
  /*   double        const _rho; */
  /*   std::vector<double> _C; */
  /* }; */

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>4</SUP>)
   * non-singlet-plus unpolarised splitting function. Parameterisation
   * determined in https://arxiv.org/pdf/1707.08315.pdf
   */
  class P3nsp: public Expression
  {
  public:
    P3nsp(int const& nf, int const& imod = 0);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
    int const _imod;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>4</SUP>)
   * non-singlet-minus unpolarised splitting
   * function. Parameterisation determined in
   * https://arxiv.org/pdf/1707.08315.pdf
   */
  class P3nsm: public Expression
  {
  public:
    P3nsm(int const& nf, int const& imod = 0);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
    int const _imod;
  };

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>4</SUP>)
   * non-singlet-valence unpolarised splitting
   * function. Parameterisation determined in
   * https://arxiv.org/pdf/1707.08315.pdf
   */
  class P3nss: public Expression
  {
  public:
    P3nss(int const& nf, int const& imod = 0);
    double Regular(double const& x) const;
  private:
    int const _nf;
    int const _imod;
  };

  /* /\** */
  /*  * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>4</SUP>) pure-singlet */
  /*  * unpolarised splitting function. */
  /*  *\/ */
  /* class P3ps: public Expression */
  /* { */
  /* public: */
  /*   P3ps(int const& nf, int const& imod = 0); */
  /*   double Regular(double const& x) const; */
  /* private: */
  /*   int const _nf; */
  /*   int const _imod; */
  /* }; */

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>4</SUP>) pure-singlet
   * unpolarised splitting function. Parameterisation determined in
   * https://arxiv.org/pdf/2302.07593.pdf
   */
  class P3ps: public Expression
  {
  public:
    P3ps(int const& nf, int const& imod = 0);
    double Regular(double const& x) const;
  private:
    int const _nf;
    int const _imod;
  };

  /* /\** */
  /*  * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>4</SUP>) quark-gluon unpolarised */
  /*  * splitting function. */
  /*  *\/ */
  /* class P3qg: public Expression */
  /* { */
  /* public: */
  /*   P3qg(int const& nf, double const& rho = -1.754); */
  /*   double Regular(double const& x) const; */
  /* private: */
  /*   int           const _nf; */
  /*   double        const _rho; */
  /*   std::vector<double> _C; */
  /* }; */

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>4</SUP>) quark-gluon
   * unpolarised splitting function. Parameterisation determined in
   * https://arxiv.org/pdf/2307.04158.pdf
   */
  class P3qg: public Expression
  {
  public:
    P3qg(int const& nf, int const& imod = 0);
    double Regular(double const& x) const;
  private:
    int const _nf;
    int const _imod;
  };

  /* /\** */
  /*  * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>4</SUP>) gluon-quark unpolarised */
  /*  * splitting function. */
  /*  *\/ */
  /* class P3gq: public Expression */
  /* { */
  /* public: */
  /*   P3gq(double const& rho = -1.784); */
  /*   double Regular(double const& x) const; */
  /* private: */
  /*   double        const _rho; */
  /*   std::vector<double> _C; */
  /* }; */

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>4</SUP>) gluon-quark
   * unpolarised splitting function. Parameterisation determined in
   * https://arxiv.org/pdf/2310.05744.pdf
   */
  class P3gq: public Expression
  {
  public:
    P3gq(int const& nf, int const& imod = 0);
    double Regular(double const& x) const;
  private:
    int const _nf;
    int const _imod;
  };

  /* /\** */
  /*  * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>4</SUP>) gluon-gluon unpolarised */
  /*  * splitting function. */
  /*  *\/ */
  /* class P3gg: public Expression */
  /* { */
  /* public: */
  /*   P3gg(double const& rho = 19.245); */
  /*   double Regular(double const& x)  const; */
  /* private: */
  /*   double        const _rho; */
  /*   std::vector<double> _C; */
  /* }; */

  /**
   * @brief Space-like O(&alpha;<SUB>s</SUB><SUP>4</SUP>) gluon-gluon
   * unpolarised splitting function. Parameterisation determined in
   * https://arxiv.org/pdf/2310.05744.pdf
   */
  class P3gg: public Expression
  {
  public:
    P3gg(int const& nf, int const& imod = 0);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
    int const _imod;
    double _A4gluon;
  };
  ///@}
}
