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
   * Fortran functions for the O(as^2) coefficient functions from
   * 'src/dis/hqcoef.f'. Use the analytic expressions whenever
   * available.
   */
  extern"C" {
    double c2log_(double *wr,double *xi);
    double cllog_(double *wr,double *xi);
    double d2nloq_(double *wr,double *xi);
    double dlnloq_(double *wr,double *xi);
    double c2nlog_(double *wr,double *xi);
    double clnlog_(double *wr,double *xi);
    double c2nloq_(double *wr,double *xi);
    double clnloq_(double *wr,double *xi);
    double c2nlobarg_(double *wr,double *xi);
    double clnlobarg_(double *wr,double *xi);
    double c2nlobarq_(double *wr,double *xi);
    double clnlobarq_(double *wr,double *xi);
  }

  /**
   * Massive Netral current coefficient functions. In the following
   * 'xi' indicates the ratio Q^2 / mh^2.
   */

  /**
   * @brief O(as) gluon coefficient function for F2. See eq. (53) of
   * https://arxiv.org/pdf/1001.2312.pdf. Or it uses the fortran
   * routines in 'src/dis/hqcoef.f'
   */
  //_________________________________________________________________________________
  class Cm21gNC: public Expression
  {
  public:
    Cm21gNC(double const& eta);
    double Regular(double const& x) const;
  };

  /**
   * @brief O(as) gluon coefficient function for FL. Analytic
   * expression. Or it uses the fortran routines in 'src/dis/hqcoef.f'
   */
  //_________________________________________________________________________________
  class CmL1gNC: public Expression
  {
  public:
    CmL1gNC(double const& eta);
    double Regular(double const& x) const;
  };

  /**
   * @brief O(as^2) non-singlet coefficient function for F2. See
   * Appendix A of https://arxiv.org/pdf/hep-ph/9601302.pdf. Or it
   * uses the fortran routines in 'src/dis/hqcoef.f'
   */
  //_________________________________________________________________________________
  class Cm22nsNC: public Expression
  {
  public:
    Cm22nsNC(double const& eta);
    double Regular(double const& x) const;
    double Local(double const&)     const;
  private:
    double _adler;
  };

  /**
   * @brief O(as^2) non-singlet coefficient function for FL. See
   * Appendix A of https://arxiv.org/pdf/hep-ph/9601302.pdf. Or it
   * uses the fortran routines in 'src/dis/hqcoef.f'
   */
  //_________________________________________________________________________________
  class CmL2nsNC: public Expression
  {
  public:
    CmL2nsNC(double const& eta);
    double Regular(double const& x) const;
  };

  /**
   * @brief O(as^2) gluon coefficient function for F2. Uses the
   * fortran routines in 'src/dis/hqcoef.f'
   */
  //_________________________________________________________________________________
  class Cm22gNC: public Expression
  {
  public:
    Cm22gNC(double const& eta);
    double Regular(double const& x) const;
  };

  /**
   * @brief O(as^2) gluon coefficient function for FL. Uses the
   * fortran routines in 'src/dis/hqcoef.f'
   */
  //_________________________________________________________________________________
  class CmL2gNC: public Expression
  {
  public:
    CmL2gNC(double const& eta);
    double Regular(double const& x) const;
  };

  /**
   * @brief O(as^2) pure-singlet coefficient function for F2. Uses the
   * fortran routines in 'src/dis/hqcoef.f'
   */
  //_________________________________________________________________________________
  class Cm22psNC: public Expression
  {
  public:
    Cm22psNC(double const& eta);
    double Regular(double const& x) const;
  };

  /**
   * @brief O(as^2) pure-singlet coefficient function for FL. Uses the
   * fortran routines in 'src/dis/hqcoef.f'
   */
  //_________________________________________________________________________________
  class CmL2psNC: public Expression
  {
  public:
    CmL2psNC(double const& eta);
    double Regular(double const& x) const;
  };

  /**
   * @brief O(as^2) gluon coefficient function proportional to
   * ln(q2/m2) for F2. Uses the fortran routines in 'src/dis/hqcoef.f'
   */
  //_________________________________________________________________________________
  class Cm22bargNC: public Expression
  {
  public:
    Cm22bargNC(double const& eta);
    double Regular(double const& x) const;
  };

  /**
   * @brief O(as^2) gluon coefficient function proportional to
   * ln(q2/m2) for FL. Uses the fortran routines in 'src/dis/hqcoef.f'
   */
  //_________________________________________________________________________________
  class CmL2bargNC: public Expression
  {
  public:
    CmL2bargNC(double const& eta);
    double Regular(double const& x) const;
  };

  /**
   * @brief O(as^2) pure-singlet coefficient function proportional to
   * ln(q2/m2) for F2. Uses the fortran routines in 'src/dis/hqcoef.f'
   */
  //_________________________________________________________________________________
  class Cm22barpsNC: public Expression
  {
  public:
    Cm22barpsNC(double const& eta);
    double Regular(double const& x) const;
  };

  /**
   * @brief O(as^2) pure-singlet coefficient function proportional to
   * ln(q2/m2) for FL. Uses the fortran routines in 'src/dis/hqcoef.f'
   */
  //_________________________________________________________________________________
  class CmL2barpsNC: public Expression
  {
  public:
    CmL2barpsNC(double const& eta);
    double Regular(double const& x) const;
  };
}
