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
   * @brief The LO time-like splitting function classes
   */
  //_________________________________________________________________________________
  class P0Tns: public Expression
  {
  public:
    P0Tns();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  //_________________________________________________________________________________
  class P0Tqg: public Expression
  {
  public:
    P0Tqg();
    double Regular(double const& x)  const;
  };

  //_________________________________________________________________________________
  class P0Tgq: public Expression
  {
  public:
    P0Tgq();
    double Regular(double const& x)  const;
  };

  //_________________________________________________________________________________
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

  /**
   * @brief The NLO time-like splitting function classes
   */
  //_________________________________________________________________________________
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

  //_________________________________________________________________________________
  class P1Tnsm: public P1Tnsp
  {
  public:
    P1Tnsm(int const& nf);
    double Regular(double const& x)  const;
  };

  //_________________________________________________________________________________
  class P1Tps: public Expression
  {
  public:
    P1Tps(int const& nf);
    double Regular(double const& x)  const;
  private:
    int const _nf;
  };

  //_________________________________________________________________________________
  class P1Tqg: public Expression
  {
  public:
    P1Tqg(int const& nf);
    double Regular(double const& x)  const;
  private:
    int const _nf;
  };

  //_________________________________________________________________________________
  class P1Tgq: public Expression
  {
  public:
    P1Tgq(int const& nf);
    double Regular(double const& x)  const;
  private:
    int const _nf;
  };

  //_________________________________________________________________________________
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

  /**
   * @brief The NNLO time-like splitting function classes (parametrized)
   */
  //_________________________________________________________________________________
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

  //_________________________________________________________________________________
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

  //_________________________________________________________________________________
  class P2Tnss: public Expression
  {
  public:
    P2Tnss(int const& nf);
    double Regular(double const& x)  const;
  private:
    int const _nf;
  };

  //_________________________________________________________________________________
  class P2Tps: public Expression
  {
  public:
    P2Tps(int const& nf);
    double Regular(double const& x)  const;
  private:
    int const _nf;
  };

  //_________________________________________________________________________________
  class P2Tqg: public Expression
  {
  public:
    P2Tqg(int const& nf);
    double Regular(double const& x)  const;
  private:
    int const _nf;
  };

  //_________________________________________________________________________________
  class P2Tgq: public Expression
  {
  public:
    P2Tgq(int const& nf);
    double Regular(double const& x)  const;
  private:
    int const _nf;
  };

  //_________________________________________________________________________________
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
}
