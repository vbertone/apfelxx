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
   * @brief O(as) coefficient function classes for F2
   */
  //_________________________________________________________________________________
  class C21ns: public Expression
  {
  public:
    C21ns();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  //_________________________________________________________________________________
  class C21g: public Expression
  {
  public:
    C21g();
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(as^2) coefficient function classes for F2 (parametrized)
   */
  //_________________________________________________________________________________
  class C22nsp: public Expression
  {
  public:
    C22nsp(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };

  //_________________________________________________________________________________
  class C22nsm: public Expression
  {
  public:
    C22nsm(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };

  //_________________________________________________________________________________
  class C22ps: public Expression
  {
  public:
    C22ps();
    double Regular(double const& x)  const;
  };

  //_________________________________________________________________________________
  class C22g: public Expression
  {
  public:
    C22g();
    double Regular(double const& x)  const;
    double Local(double const& x)    const;
  };

  /**
   * @brief O(as) coefficient function classes for FL
   */
  //_________________________________________________________________________________
  class CL1ns: public Expression
  {
  public:
    CL1ns();
    double Regular(double const& x)  const;
  };

  //_________________________________________________________________________________
  class CL1g: public Expression
  {
  public:
    CL1g();
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(as^2) coefficient function classes for FL (parametrized)
   */
  //_________________________________________________________________________________
  class CL2nsp: public Expression
  {
  public:
    CL2nsp(int const& nf);
    double Regular(double const& x)  const;
    double Local(double const&)      const;
  private:
    int const _nf;
  };

  //_________________________________________________________________________________
  class CL2nsm: public Expression
  {
  public:
    CL2nsm(int const& nf);
    double Regular(double const& x)  const;
    double Local(double const&)      const;
  private:
    int const _nf;
  };

  //_________________________________________________________________________________
  class CL2ps: public Expression
  {
  public:
    CL2ps();
    double Regular(double const& x)  const;
  };

  //_________________________________________________________________________________
  class CL2g: public Expression
  {
  public:
    CL2g();
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(as) coefficient function classes for F3
   */
  //_________________________________________________________________________________
  class C31ns: public Expression
  {
  public:
    C31ns();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief O(as^2) coefficient function classes for F3 (parametrized)
   */
  //_________________________________________________________________________________
  class C32nsp: public Expression
  {
  public:
    C32nsp(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };

  //_________________________________________________________________________________
  class C32nsm: public Expression
  {
  public:
    C32nsm(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };
}
