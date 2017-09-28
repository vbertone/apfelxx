//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <apfel/expression.h>

namespace apfel
{
  /**
   * @brief The NLO matching function classes
   * @brief Reference: arXiv:1604.07869
   */
  //_________________________________________________________________________________
  class C1ns: public Expression
  {
  public:
    C1ns();
    double Regular(double const& x) const;
    double Local(double const&)     const;
  };

  //_________________________________________________________________________________
  class C1qg: public Expression
  {
  public:
    C1qg();
    double Regular(double const& x) const;
  };

  //_________________________________________________________________________________
  class C1gq: public Expression
  {
  public:
    C1gq();
    double Regular(double const& x) const;
  };

  //_________________________________________________________________________________
  class C1gg: public Expression
  {
  public:
    C1gg();
    double Local(double const&) const;
  };

  /* /\** */
  /*  * @brief The NLO splitting function classes */
  /*  *\/ */
  /* //_________________________________________________________________________________ */
  /* class P1nsp: public Expression */
  /* { */
  /* public: */
  /*   P1nsp(int const& nf); */
  /*   double Regular(double const& x)  const; */
  /*   double Singular(double const& x) const; */
  /*   double Local(double const& x)    const; */
  /* protected: */
  /*   int const _nf; */
  /*   double    _a2; */
  /* }; */

  /* //_________________________________________________________________________________ */
  /* class P1nsm: public P1nsp */
  /* { */
  /* public: */
  /*   P1nsm(int const& nf); */
  /*   double Regular(double const& x)  const; */
  /* }; */

  /* //_________________________________________________________________________________ */
  /* class P1ps: public Expression */
  /* { */
  /* public: */
  /*   P1ps(int const& nf); */
  /*   double Regular(double const& x)  const; */
  /* private: */
  /*   int const _nf; */
  /* }; */

  /* //_________________________________________________________________________________ */
  /* class P1qg: public Expression */
  /* { */
  /* public: */
  /*   P1qg(int const& nf); */
  /*   double Regular(double const& x)  const; */
  /* private: */
  /*   int const _nf; */
  /* }; */

  /* //_________________________________________________________________________________ */
  /* class P1gq: public Expression */
  /* { */
  /* public: */
  /*   P1gq(int const& nf); */
  /*   double Regular(double const& x)  const; */
  /* private: */
  /*   int const _nf; */
  /* }; */

  /* //_________________________________________________________________________________ */
  /* class P1gg: public Expression */
  /* { */
  /* public: */
  /*   P1gg(int const& nf); */
  /*   double Regular(double const& x)  const; */
  /*   double Singular(double const& x) const; */
  /*   double Local(double const& x)    const; */
  /* private: */
  /*   int const _nf; */
  /*   double    _a2g; */
  /* }; */

  /* /\** */
  /*  * @brief The NNLO splitting function classes (parametrized) */
  /*  *\/ */
  /* //_________________________________________________________________________________ */
  /* class P2nsp: public Expression */
  /* { */
  /* public: */
  /*   P2nsp(int const& nf); */
  /*   double Regular(double const& x)  const; */
  /*   double Singular(double const& x) const; */
  /*   double Local(double const& x)    const; */
  /* private: */
  /*   int const _nf; */
  /* }; */

  /* //_________________________________________________________________________________ */
  /* class P2nsm: public Expression */
  /* { */
  /* public: */
  /*   P2nsm(int const& nf); */
  /*   double Regular(double const& x)  const; */
  /*   double Singular(double const& x) const; */
  /*   double Local(double const& x)    const; */
  /* private: */
  /*   int const _nf; */
  /* }; */

  /* //_________________________________________________________________________________ */
  /* class P2nss: public Expression */
  /* { */
  /* public: */
  /*   P2nss(int const& nf); */
  /*   double Regular(double const& x)  const; */
  /* private: */
  /*   int const _nf; */
  /* }; */

  /* //_________________________________________________________________________________ */
  /* class P2ps: public Expression */
  /* { */
  /* public: */
  /*   P2ps(int const& nf); */
  /*   double Regular(double const& x)  const; */
  /* private: */
  /*   int const _nf; */
  /* }; */

  /* //_________________________________________________________________________________ */
  /* class P2qg: public Expression */
  /* { */
  /* public: */
  /*   P2qg(int const& nf); */
  /*   double Regular(double const& x)  const; */
  /* private: */
  /*   int const _nf; */
  /* }; */

  /* //_________________________________________________________________________________ */
  /* class P2gq: public Expression */
  /* { */
  /* public: */
  /*   P2gq(int const& nf); */
  /*   double Regular(double const& x)  const; */
  /* private: */
  /*   int const _nf; */
  /* }; */

  /* //_________________________________________________________________________________ */
  /* class P2gg: public Expression */
  /* { */
  /* public: */
  /*   P2gg(int const& nf); */
  /*   double Regular(double const& x)  const; */
  /*   double Singular(double const& x) const; */
  /*   double Local(double const& x)    const; */
  /* private: */
  /*   int const _nf; */
  /* }; */
}
