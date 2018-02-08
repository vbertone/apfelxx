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
   * @brief The LO space-like splitting function classes
   */
  //_________________________________________________________________________________
  class P0ns: public Expression
  {
  public:
    P0ns();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  //_________________________________________________________________________________
  class P0qg: public Expression
  {
  public:
    P0qg();
    double Regular(double const& x)  const;
  };

  //_________________________________________________________________________________
  class P0gq: public Expression
  {
  public:
    P0gq();
    double Regular(double const& x)  const;
  };

  //_________________________________________________________________________________
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

  /**
   * @brief The NLO space-like splitting function classes
   */
  //_________________________________________________________________________________
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

  //_________________________________________________________________________________
  class P1nsm: public P1nsp
  {
  public:
    P1nsm(int const& nf);
    double Regular(double const& x)  const;
  };

  //_________________________________________________________________________________
  class P1ps: public Expression
  {
  public:
    P1ps(int const& nf);
    double Regular(double const& x)  const;
  private:
    int const _nf;
  };

  //_________________________________________________________________________________
  class P1qg: public Expression
  {
  public:
    P1qg(int const& nf);
    double Regular(double const& x)  const;
  private:
    int const _nf;
  };

  //_________________________________________________________________________________
  class P1gq: public Expression
  {
  public:
    P1gq(int const& nf);
    double Regular(double const& x)  const;
  private:
    int const _nf;
  };

  //_________________________________________________________________________________
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

  /**
   * @brief The NNLO space-like splitting function classes (parametrized)
   */
  //_________________________________________________________________________________
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

  //_________________________________________________________________________________
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

  //_________________________________________________________________________________
  class P2nss: public Expression
  {
  public:
    P2nss(int const& nf);
    double Regular(double const& x)  const;
  private:
    int const _nf;
  };

  //_________________________________________________________________________________
  class P2ps: public Expression
  {
  public:
    P2ps(int const& nf);
    double Regular(double const& x)  const;
  private:
    int const _nf;
  };

  //_________________________________________________________________________________
  class P2qg: public Expression
  {
  public:
    P2qg(int const& nf);
    double Regular(double const& x)  const;
  private:
    int const _nf;
  };

  //_________________________________________________________________________________
  class P2gq: public Expression
  {
  public:
    P2gq(int const& nf);
    double Regular(double const& x)  const;
  private:
    int const _nf;
  };

  //_________________________________________________________________________________
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

  /**
   * @brief The NNNLO space-like splitting function classes
   * (parametrized and leading color). Only the +, -, and valence
   * contributions have been computed so far.
   *
   */
  //_________________________________________________________________________________
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

  //_________________________________________________________________________________
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

  //_________________________________________________________________________________
  class P3nss: public Expression
  {
  public:
    P3nss(int const& nf, int const& imod = 0);
    double Regular(double const& x)  const;
  private:
    int const _nf;
    int const _imod;
  };

  /**
   * @brief The LO space-like splitting function for tranversely
   * polarised PDFs. Reference
   * https://arxiv.org/pdf/hep-ph/9706511v2.pdf.
   */
  //_________________________________________________________________________________
  class P0transns: public Expression
  {
  public:
    P0transns();
    double Regular(double const&)    const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief The NLO space-like splitting function for tranversely
   * polarised PDFs. Reference
   * https://arxiv.org/pdf/hep-ph/9706511v2.pdf.
   */
  //_________________________________________________________________________________
  class P1transnsp: public Expression
  {
  public:
    P1transnsp(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  protected:
    int const _nf;
    double    _a2;
  };

  //_________________________________________________________________________________
  class P1transnsm: public P1transnsp
  {
  public:
    P1transnsm(int const& nf);
    double Regular(double const& x)  const;
  };
}
