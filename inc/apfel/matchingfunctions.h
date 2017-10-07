//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/expression.h"

namespace apfel
{
  /**
   * @brief The NLO matching function classes
   * @brief References: arXiv:1604.07869 and arXiv:1706.01473.
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

  /**
   * @brief The NNLO matching function classes. In the implementation
   * of these expressions I used the parameterization only for the
   * regular part. The exact coefficicients of the rest are
   * implemented explicitly as reported in appendix B.3 of
   * https://arxiv.org/pdf/1706.01473.pdf.
   */
  //_________________________________________________________________________________
  class C2Vqq: public Expression
  {
  public:
    C2Vqq(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  protected:
    int const _nf;
    double    _A2;
    double    _A3;
  };

  //_________________________________________________________________________________
  class C2Vqqb: public Expression
  {
  public:
    C2Vqqb();
    double Regular(double const& x)  const;
  };

  //_________________________________________________________________________________
  class C2ps: public Expression
  {
  public:
    C2ps();
    double Regular(double const& x)  const;
  };

  //_________________________________________________________________________________
  class C2qg: public Expression
  {
  public:
    C2qg();
    double Regular(double const& x)  const;
  };

  //_________________________________________________________________________________
  class C2gq: public Expression
  {
  public:
    C2gq(int const& nf);
    double Regular(double const& x)  const;
  private:
    int const _nf;
  };

  //_________________________________________________________________________________
  class C2gg: public Expression
  {
  public:
    C2gg(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
    double    _A2;
    double    _A3;
  };
}
