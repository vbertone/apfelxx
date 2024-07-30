//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/expression.h"

namespace apfel
{
  /// @cond UNNECESSARY
  /**
   * @name Fortran massive coefficient functions
   * Fortran functions for the O(&alpha;<SUB>s</SUB><SUP>2</SUP>)
   * coefficient functions from 'src/dis/hqcoef.f'. Use the analytic
   * expressions whenever available.
   */
  ///@{
  extern"C"
  {
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
  ///@}
  /// \endcond

  /**
   * @defgroup NCMassive Massive neutral-current coefficient functions
   * Collection of the neutral-current massive coefficient functions
   * for F<SUB>2</SUB> and F<SUB>L</SUB> up to O(&alpha;<SUB>s</SUB>).
   * @note In the following 'xi' indicates the ratio Q<SUP>2</SUP> /
   * M<SUP>2</SUP>.
   */
  ///@{
  /**
   * @defgroup LO LO massive coefficient functions
   * @ingroup NCMassive
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon coefficient function for
   * F2. See eq. (53) of https://arxiv.org/pdf/1001.2312.pdf. Or it
   * uses the fortran routines in 'src/dis/hqcoef.f'
   */
  class Cm21gNC: public Expression
  {
  public:
    Cm21gNC(double const& eta);
    double Regular(double const& x) const;
  private:
    double const _eta;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon coefficient function for FL.
   */
  class CmL1gNC: public Expression
  {
  public:
    CmL1gNC(double const& eta);
    double Regular(double const& x) const;
  private:
    double const _eta;
  };
  ///@}

  /**
   * @defgroup NLO NLO massive coefficient functions
   * @ingroup NCMassive
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet coefficient
   * function for F2. See Appendix A of
   * https://arxiv.org/pdf/hep-ph/9601302.pdf. Or it uses the fortran
   * routines in 'src/dis/hqcoef.f'
   */
  class Cm22nsNC: public Expression
  {
  public:
    Cm22nsNC(double const& eta);
    double Regular(double const& x) const;
    double Local(double const&)     const;
  private:
    double const _eta;
    double       _adler;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet coefficient
   * function for FL. See Appendix A of
   * https://arxiv.org/pdf/hep-ph/9601302.pdf. Or it uses the fortran
   * routines in 'src/dis/hqcoef.f'
   */
  class CmL2nsNC: public Expression
  {
  public:
    CmL2nsNC(double const& eta);
    double Regular(double const& x) const;
  private:
    double const _eta;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon coefficient
   * function for F2. Uses the fortran routines in 'src/dis/hqcoef.f'
   */
  class Cm22gNC: public Expression
  {
  public:
    Cm22gNC(double const& eta);
    double Regular(double const& x) const;
  private:
    double const _eta;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon coefficient
   * function for FL. Uses the fortran routines in 'src/dis/hqcoef.f'
   */
  class CmL2gNC: public Expression
  {
  public:
    CmL2gNC(double const& eta);
    double Regular(double const& x) const;
  private:
    double const _eta;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * coefficient function for F2. Uses the fortran routines in
   * 'src/dis/hqcoef.f'
   */
  class Cm22psNC: public Expression
  {
  public:
    Cm22psNC(double const& eta);
    double Regular(double const& x) const;
  private:
    double const _eta;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * coefficient function for FL. Uses the fortran routines in
   * 'src/dis/hqcoef.f'
   */
  class CmL2psNC: public Expression
  {
  public:
    CmL2psNC(double const& eta);
    double Regular(double const& x) const;
  private:
    double const _eta;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon coefficient
   * function proportional to ln(Q<SUP>2</SUP>/M<SUP>2</SUP>) for
   * F2. Uses the fortran routines in 'src/dis/hqcoef.f'
   */
  class Cm22bargNC: public Expression
  {
  public:
    Cm22bargNC(double const& eta);
    double Regular(double const& x) const;
  private:
    double const _eta;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon coefficient
   * function proportional to ln(Q<SUP>2</SUP>/M<SUP>2</SUP>) for
   * FL. Uses the fortran routines in 'src/dis/hqcoef.f'
   */
  class CmL2bargNC: public Expression
  {
  public:
    CmL2bargNC(double const& eta);
    double Regular(double const& x) const;
  private:
    double const _eta;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * coefficient function proportional to
   * ln(Q<SUP>2</SUP>/M<SUP>2</SUP>) for F2. Uses the fortran routines
   * in 'src/dis/hqcoef.f'
   */
  class Cm22barpsNC: public Expression
  {
  public:
    Cm22barpsNC(double const& eta);
    double Regular(double const& x) const;
  private:
    double const _eta;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * coefficient function proportional to
   * ln(Q<SUP>2</SUP>/M<SUP>2</SUP>) for FL. Uses the fortran routines
   * in 'src/dis/hqcoef.f'
   */
  class CmL2barpsNC: public Expression
  {
  public:
    CmL2barpsNC(double const& eta);
    double Regular(double const& x) const;
  private:
    double const _eta;
  };
  ///@}

  /**
   * @defgroup NNLOthr Approximated NNLO massive coefficient functions near threshold
   * @ingroup NCMassive
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) gluon coefficient
   * function for F2 near
   * threshold. Reference:https://arxiv.org/pdf/1205.5727.
   */
  class CmTh23gNC: public Expression
  {
  public:
    CmTh23gNC(int const& nf, double const& eta);
    double Regular(double const& x) const;
  private:
    int     const _nf;
    double  const _eta;
    Cm21gNC const _c21g;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) gluon coefficient
   * function for FL near
   * threshold. Reference:https://arxiv.org/pdf/1205.5727.
   */
  class CmThL3gNC: public Expression
  {
  public:
    CmThL3gNC(int const& nf, double const& eta);
    double Regular(double const& x) const;
  private:
    int     const _nf;
    double  const _eta;
    CmL1gNC const _cL1g;
  };
  ///@}

  /**
   * @defgroup NNLOsx Approximated NNLO massive coefficient functions at small x
   * @ingroup NCMassive
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) gluon coefficient
   * function for F2 at small
   * x. Reference:https://arxiv.org/pdf/1205.5727.
   */
  class Cmsx23gNC: public Expression
  {
  public:
    Cmsx23gNC(int const& nf, double const& eta);
    double Regular(double const& x) const;
  private:
    int    const _nf;
    double const _eta;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) pure-singlet
   * coefficient function for F2 at small
   * x. Reference:https://arxiv.org/pdf/1205.5727.
   */
  class Cmsx23psNC: public Expression
  {
  public:
    Cmsx23psNC(int const& nf, double const& eta);
    double Regular(double const& x) const;
  private:
    Cmsx23gNC const _c23g;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) gluon coefficient
   * function for FL at small
   * x. Reference:https://arxiv.org/pdf/1205.5727.
   */
  class CmsxL3gNC: public Expression
  {
  public:
    CmsxL3gNC(int const& nf, double const& eta);
    double Regular(double const& x) const;
  private:
    int    const _nf;
    double const _eta;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) pure-singlet
   * coefficient function for FL at small
   * x. Reference:https://arxiv.org/pdf/1205.5727.
   */
  class CmsxL3psNC: public Expression
  {
  public:
    CmsxL3psNC(int const& nf, double const& eta);
    double Regular(double const& x) const;
  private:
    CmsxL3gNC const _cL3g;
  };
  ///@}

  /**
   * @defgroup NNLOm0sx Approximated NNLO massive coefficient functions at small x and Q >> m
   * @ingroup NCMassive
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) gluon coefficient
   * function for F2 at Q >> m and small
   * x. Reference:https://arxiv.org/pdf/1205.5727.
   */
  class Cm0sx23gNC: public Expression
  {
  public:
    Cm0sx23gNC(int const& nf, double const& eta);
    double Regular(double const& x) const;
  private:
    int    const _nf;
    double const _eta;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) pure-singlet
   * coefficient function for F2 at Q >> m and small
   * x. Reference:https://arxiv.org/pdf/1205.5727.
   */
  class Cm0sx23psNC: public Expression
  {
  public:
    Cm0sx23psNC(int const& nf, double const& eta);
    double Regular(double const& x) const;
  private:
    Cm0sx23gNC const _c23g;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) gluon coefficient
   * function for FL at Q >> m and small
   * x. Reference:https://arxiv.org/pdf/1205.5727.
   */
  class Cm0sxL3gNC: public Expression
  {
  public:
    Cm0sxL3gNC(int const& nf, double const& eta);
    double Regular(double const& x) const;
  private:
    int    const _nf;
    double const _eta;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) pure-singlet
   * coefficient function for FL at Q >> m and small
   * x. Reference:https://arxiv.org/pdf/1205.5727.
   */
  class Cm0sxL3psNC: public Expression
  {
  public:
    Cm0sxL3psNC(int const& nf, double const& eta);
    double Regular(double const& x) const;
  private:
    Cm0sxL3gNC const _cL3g;
  };
  ///@}

  /**
   * @defgroup NLOhq NLO massive heavy-quark-initiated coefficient functions
   * Collection of the massive coefficient functions for processes
   * with a massive heavy quark in the intial state. These coefficient
   * functions can be used both in the neutral-current and in the
   * charged-current cases.
   * @note All the expressions of the one-loop coefficient functions
   * are extracted from https://arxiv.org/pdf/hep-ph/9805233.pdf. The
   * coefficient function C is witten in terms of the private _R
   * function such that:
   *
   * C(x) = ( _R(x) - _R(1) ) / ( 1 - x ) + _R(x) / ( 1 - x )_+ + ( L + _R(1) ) delta( 1- x )
   *
   * @ingroup NCMassive
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet coefficient function
   * for 2xF1. The relevant function is in Eq. (C4) of
   * https://arxiv.org/pdf/hep-ph/9805233.pdf.
   */
  class Cm11ns: public Expression
  {
  public:
    Cm11ns(double const& m1, double const& m2, double const& Q, double const& Splus, double const& Sminus);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    double _R(double const& x) const;
    double const _m1;
    double const _m2;
    double const _Splus;
    double const _Sminus;
    double const _m12;
    double const _m22;
    double const _Q2;
    double const _Del;
    double const _Del2;
    double const _Spp;
    double const _Spm;
    double const _Smp;
    double const _fact1;
    double _S1;
    double _V1;
    double _R1;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet coefficient function
   * for F2. The relevant function is in Eq. (C4) of
   * https://arxiv.org/pdf/hep-ph/9805233.pdf.
   */
  class Cm21ns: public Expression
  {
  public:
    Cm21ns(double const& m1, double const& m2, double const& Q, double const& Splus, double const& Sminus);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    double _R(double const& x) const;
    double const _m1;
    double const _m2;
    double const _Splus;
    double const _Sminus;
    double const _m12;
    double const _m22;
    double const _Q2;
    double const _Del;
    double const _Del2;
    double const _Spp;
    double const _Spm;
    double const _Smp;
    double const _fact2;
    double _S2;
    double _V2;
    double _R1;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet coefficient function
   * for F3. The relevant function is in Eq. (C4) of
   * https://arxiv.org/pdf/hep-ph/9805233.pdf.
   */
  class Cm31ns: public Expression
  {
  public:
    Cm31ns(double const& m1, double const& m2, double const& Q, double const& Rplus, double const& Rminus);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    double _R(double const& x) const;
    double const _m1;
    double const _m2;
    double const _Rplus;
    double const _Rminus;
    double const _m12;
    double const _m22;
    double const _Q2;
    double const _Del;
    double const _Del2;
    double const _Spp;
    double const _Spm;
    double const _Smp;
    double const _fact3;
    double _S3;
    double _V3;
    double _R1;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet coefficient function
   * for FL = F2 - 2xF1. The relevant function is in Eq. (C4) of
   * https://arxiv.org/pdf/hep-ph/9805233.pdf.
   */
  class CmL1ns: public Expression
  {
  public:
    CmL1ns(double const& m1, double const& m2, double const& Q, double const& Splus, double const& Sminus);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    Cm11ns const _C1;
    Cm21ns const _C2;
    double const _factL;
  };
  ///@}
  ///@}

  /**
   * @defgroup CCMassive Massive charged-current coefficient functions
   * Collection of the charged-current massive coefficient functions
   * for F<SUB>2</SUB>, F<SUB>2</SUB>, and xF<SUB>3</SUB> to
   * O(&alpha;<SUB>s</SUB>). Expressions taken from
   * https://arxiv.org/pdf/hep-ph/9603304.pdf.
   * @note In the following 'xi' indicates the ratio Q<SUP>2</SUP> /
   * M<SUP>2</SUP>.
   */
  ///@{
  /**
   * @defgroup NLO NLO massive coefficient functions
   * @ingroup CCMassive
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB>) quark coefficient function for F2.
   */
  class Cm21qCC: public Expression
  {
  public:
    Cm21qCC(double const& lambda);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    double const _lambda;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon coefficient function for F2.
   */
  class Cm21gCC: public Expression
  {
  public:
    Cm21gCC(double const& lambda);
    double Regular(double const& x) const;
  private:
    double const _lambda;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) quark coefficient function for FL.
   */
  class CmL1qCC: public Expression
  {
  public:
    CmL1qCC(double const& lambda);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    double const _lambda;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon coefficient function for FL.
   */
  class CmL1gCC: public Expression
  {
  public:
    CmL1gCC(double const& lambda);
    double Regular(double const& x) const;
  private:
    double const _lambda;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) quark coefficient function for xF3.
   */
  class Cm31qCC: public Expression
  {
  public:
    Cm31qCC(double const& lambda);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    double const _lambda;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon coefficient function for xF3.
   */
  class Cm31gCC: public Expression
  {
  public:
    Cm31gCC(double const& lambda);
    double Regular(double const& x) const;
  private:
    double const _lambda;
  };
  ///@}
}
