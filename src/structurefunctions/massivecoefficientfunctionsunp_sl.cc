//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/massivecoefficientfunctionsunp_sl.h"
#include "apfel/constants.h"
#include "apfel/specialfunctions.h"
#include "apfel/integrator.h"
#include "apfel/tools.h"

namespace apfel
{
  //_________________________________________________________________________________
  Cm21gNC::Cm21gNC(double const& eta):
    Expression(),
    _eta(eta)
  {
  }
  double Cm21gNC::Regular(double const& x) const
  {
    if (x >= 1)
      return 0;
    /*
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi  = 4 * eta / ( 1 - eta );
    double wr  = xi * ( 1 / z - 1 ) / 4 - 1;
    const double cm21g = xi * c2log_(&wr, &xi) / z / M_PI;
    return _eta * cm21g;
    */
    const double eta   = this->_eta;
    const double z     = eta * x;
    const double z2    = z * z;
    const double epsi  = ( 1 - eta ) / eta / 4;
    const double epsi2 = epsi * epsi;
    const double v     = sqrt( 1 - 4 * epsi * z / ( 1 - z ));
    return _eta * 4 * TR * ( ( z2 + ( 1 - z ) * ( 1 - z ) + 4 * epsi * z * ( 1 - 3 * z )
                               - 8 * epsi2 * z2 ) * log( ( 1 + v ) / ( 1 - v ) )
                             + ( 8 * z * ( 1 - z ) - 1 - 4 * epsi * z * ( 1 - z ) ) * v );
  }

  //_________________________________________________________________________________
  CmL1gNC::CmL1gNC(double const& eta):
    Expression(),
    _eta(eta)
  {
  }
  double CmL1gNC::Regular(double const& x) const
  {
    if (x >= 1)
      return 0;
    /*
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi  = 4 * eta / ( 1 - eta );
    double wr  = xi * ( 1 / z - 1 ) / 4 - 1;
    const double cml1g = xi * cllog_(&wr, &xi) / z / M_PI;
    return _eta * cml1g;
    */
    const double eta  = this->_eta;
    const double z    = eta * x;
    const double z2   = z * z;
    const double epsi = ( 1 - eta ) / eta / 4;
    const double v    = sqrt( 1 - 4 * z / ( 1 - z ) * epsi );
    return _eta * 4 * TR * ( - 8 * epsi * z2 * log( ( 1 + v ) / ( 1 - v ) )
                             + 4 * v * z * ( 1 - z ) );
  }

  //_________________________________________________________________________________
  Cm22nsNC::Cm22nsNC(double const& eta):
    Expression(),
    _eta(eta)
  {
    // Compute integral needed to enforce the Adler sum rule. See
    // eqs. (71) and (97) of https://arxiv.org/pdf/1001.2312.pdf.
    const Integrator Integrand{[&] (double const& y) -> double{ return Regular(y); }};
    _adler = - Integrand.integrate(0, 1, eps5);
  }
  double Cm22nsNC::Regular(double const& x) const
  {
    if (x >= 1)
      return 0;
    /*
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi  = 4 * eta / ( 1 - eta );
    double wr  = xi * ( 1 / z - 1 ) / 4 - 1;
    const double cm22ns = 16 * M_PI * xi * d2nloq_(&wr, &xi) / z;
    return _eta * cm22ns;
    */
    const double eta = this->_eta;
    const double xi  = 4 * eta / ( 1 - eta );
    const double z   = eta * x;
    const double z2  = z * z;
    const double z3  = z * z2;
    const double omz = 1 - z;
    const double zr2  = z2 / xi / xi;
    const double qr   = z / omz / xi;
    const double qr2  = qr * qr;
    const double sq1  = sqrt( 1 - 4 * z / xi / omz );
    const double sq2  = sqrt( 1 - 4 * z / xi );
    if (sq1 == sq2)
      return 0;
    const double l1   = log( ( 1 + sq1 ) / ( 1 - sq1 ) );
    const double l2   = log( ( 1 + sq2 ) / ( 1 - sq2 ) );
    const double l3   = log( ( sq2 + sq1 ) / ( sq2 - sq1 ) );
    const double dil1 = dilog( omz * ( 1 + sq1 ) / ( 1 + sq2 ) );
    const double dil2 = dilog( ( 1 - sq2 ) / ( 1 + sq1 ) );
    const double dil3 = dilog( ( 1 - sq1 ) / ( 1 + sq2 ) );
    const double dil4 = dilog( ( 1 + sq1 ) / ( 1 + sq2 ) );
    return _eta * CF * TR * ( ( 4. / 3 * ( 1 + z2 ) / omz - 16 / omz * zr2 * ( 1 - 9 * z + 9 * z2 ) )
                              * ( log( omz / z2 ) * l1 + l1 * l2 + 2 * ( - dil1 + dil2 + dil3 - dil4 ) )
                              + ( - 8. / 3 + 4 / omz + qr2 * ( 128. - 432 * z + 288 * z2 - 8 / omz ) ) * l1
                              + ( 88. / 9 + 136. / 9 * z - 152. / 9 / omz
                                  + qr * ( 464. / 9 - 512. / 3 * z + 2048. / 9 * z2 )
                                  + qr2 * ( - 832. / 9 + 6208. / 9 * z - 11392. / 9 * z2 + 6016. / 9 * z3 ) ) * l3 / sq2
                              + ( - 272. / 27 - 1244. / 27 * z + 718. / 27 / omz
                                  + qr * ( - 3424. / 27 + 15608. / 27 * z - 4304. / 9 * z2 + 20. / 27 / omz ) ) * sq1 );
  }
  double Cm22nsNC::Local(double const&) const
  {
    return _eta * _adler;
  }

  //_________________________________________________________________________________
  CmL2nsNC::CmL2nsNC(double const& eta):
    Expression(),
    _eta(eta)
  {
  }
  double CmL2nsNC::Regular(double const& x) const
  {
    if (x >= 1)
      return 0;
    /*
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi  = 4 * eta / ( 1 - eta );
    double wr  = xi * ( 1 / z - 1 ) / 4 - 1;
    const double cml2ns = 16 * M_PI * xi * dlnloq_(&wr, &xi) / z;
    return _eta * cml2ns;
    */
    const double eta = this->_eta;
    const double xi  = 4 * eta / ( 1 - eta );
    const double z   = eta * x;
    const double z2  = z * z;
    const double omz = 1 - z;
    const double zr2 = z2 / xi / xi;
    const double qr  = z / omz / xi;
    const double qr2 = qr * qr;
    const double sq1 = sqrt( 1 - 4 * z / xi / omz );
    const double sq2 = sqrt( 1 - 4 * z / xi );
    if (sq1 == sq2)
      return 0;
    const double l1   = log( ( 1 + sq1 ) / ( 1 - sq1 ) );
    const double l2   = log( ( 1 + sq2 ) / ( 1 - sq2 ) );
    const double l3   = log( ( sq2 + sq1 ) / ( sq2 - sq1 ) );
    const double dil1 = dilog( omz * ( 1 + sq1 ) / ( 1 + sq2 ) );
    const double dil2 = dilog( ( 1 - sq2 ) / ( 1 + sq1 ) );
    const double dil3 = dilog( ( 1 - sq1 ) / ( 1 + sq2 ) );
    const double dil4 = dilog( ( 1 + sq1 ) / ( 1 + sq2 ) );
    return _eta * CF * TR * ( 96 * z * zr2 * ( log( omz / z2 ) * l1 + l1 * l2 + 2 * ( - dil1 + dil2 + dil3 - dil4 ) )
                              + qr2 * ( 64 - 288 * z + 192 * z2 ) * l1
                              + z * ( 16. / 3 - 416 * z / 3 / xi + 1408 * z2 / 3 / xi / xi ) * l3 / sq2
                              + ( 16. / 3 - 400 * z / 18
                                  + z * ( - 160. / 3 + 3824 * z / 9 - 992 * z2 / 3 ) / omz / xi ) * sq1);
  }

  //_________________________________________________________________________________
  Cm22gNC::Cm22gNC(double const& eta):
    Expression(),
    _eta(eta)
  {
  }
  double Cm22gNC::Regular(double const& x) const
  {
    if (x >= 1)
      return 0;
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi = 4 * eta / ( 1 - eta );
    double wr = xi * ( 1 / z - 1 ) / 4 - 1;
    return _eta * 16 * M_PI * xi * c2nlog_(&wr, &xi) / z;
  }

  //_________________________________________________________________________________
  CmL2gNC::CmL2gNC(double const& eta):
    Expression(),
    _eta(eta)
  {
  }
  double CmL2gNC::Regular(double const& x) const
  {
    if (x >= 1)
      return 0;
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi = 4 * eta / ( 1 - eta );
    double wr = xi * ( 1 / z - 1 ) / 4 - 1;
    return _eta * 16 * M_PI * xi * clnlog_(&wr, &xi) / z;
  }

  //_________________________________________________________________________________
  Cm22psNC::Cm22psNC(double const& eta):
    Expression(),
    _eta(eta)
  {
  }
  double Cm22psNC::Regular(double const& x) const
  {
    if (x >= 1)
      return 0;
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi = 4 * eta / ( 1 - eta );
    double wr = xi * ( 1 / z - 1 ) / 4 - 1;
    return _eta * 16 * M_PI * xi * c2nloq_(&wr, &xi) / z;
  }

  //_________________________________________________________________________________
  CmL2psNC::CmL2psNC(double const& eta):
    Expression(),
    _eta(eta)
  {
  }
  double CmL2psNC::Regular(double const& x) const
  {
    if (x >= 1)
      return 0;
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi = 4 * eta / ( 1 - eta );
    double wr = xi * ( 1 / z - 1 ) / 4 - 1;
    return _eta * 16 * M_PI * xi * clnloq_(&wr, &xi) / z;
  }

  //_________________________________________________________________________________
  Cm22bargNC::Cm22bargNC(double const& eta):
    Expression(),
    _eta(eta)
  {
  }
  double Cm22bargNC::Regular(double const& x) const
  {
    if (x >= 1)
      return 0;
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi = 4 * eta / ( 1 - eta );
    double wr = xi * ( 1 / z - 1 ) / 4 - 1;
    return _eta * 16 * M_PI * xi * c2nlobarg_(&wr, &xi) / z;
  }

  //_________________________________________________________________________________
  CmL2bargNC::CmL2bargNC(double const& eta):
    Expression(),
    _eta(eta)
  {
  }
  double CmL2bargNC::Regular(double const& x) const
  {
    if (x >= 1)
      return 0;
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi = 4 * eta / ( 1 - eta );
    double wr = xi * ( 1 / z - 1 ) / 4 - 1;
    return _eta * 16 * M_PI * xi * clnlobarg_(&wr, &xi) / z;
  }

  //_________________________________________________________________________________
  Cm22barpsNC::Cm22barpsNC(double const& eta):
    Expression(),
    _eta(eta)
  {
  }
  double Cm22barpsNC::Regular(double const& x) const
  {
    if (x >= 1)
      return 0;
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi = 4 * eta / ( 1 - eta );
    double wr = xi * ( 1 / z - 1 ) / 4 - 1;
    return _eta * 16 * M_PI * xi * c2nlobarq_(&wr, &xi) / z;
  }

  //_________________________________________________________________________________
  CmL2barpsNC::CmL2barpsNC(double const& eta):
    Expression(),
    _eta(eta)
  {
  }
  double CmL2barpsNC::Regular(double const& x) const
  {
    if (x >= 1)
      return 0;
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi = 4 * eta / ( 1 - eta );
    double wr = xi * ( 1 / z - 1 ) / 4 - 1;
    return _eta * 16 * M_PI * xi * clnlobarq_(&wr, & xi) / z;
  }

  //_________________________________________________________________________________
  Cm11ns::Cm11ns(double const& m1, double const& m2, double const& Q, double const& Splus, double const& Sminus):
    Expression(2 * pow(Q, 2) / ( pow(Q, 2) + pow(m2, 2) - pow(m1, 2) + DeltaFun(pow(m1, 2), pow(m2, 2), - pow(Q, 2)))),
    _m1(m1),
    _m2(m2),
    _Splus(Splus),
    _Sminus(Sminus),
    _m12(_m1 * _m1),
    _m22(_m2 * _m2),
    _Q2(Q * Q),
    _Del(DeltaFun(_m12, _m22, - _Q2)),
    _Del2(_Del * _Del),
    _Spp(_Q2 + _m22 + _m12),
    _Spm(_Q2 + _m22 - _m12),
    _Smp(_Q2 - _m22 + _m12),
    _fact1(2 * CF * ( _Spp - 2 * _m1 * _m2 * _Sminus / _Splus ) / _Del)
  {
    const double I1    = log( ( _Spp + _Del ) / ( _Spp - _Del ) ) / _Del;
    const double Cplus = 2 * m1 * m2 * I1;
    const double CRm   = ( _Del2 / 2 / _Q2 + _Spp * ( 1 + log( _Q2 / _Del ) ) ) * I1
                         + ( _m22 - _m12 ) / 2 / _Q2 * log( _m12 / _m22 )
                         - log( _Q2 / _m12 ) - log( _Q2 / _m22 ) - 4
                         + _Spp / _Del * ( + pow(log( dabs( ( _Del - _Spm ) / 2 / _Q2 ) ), 2) / 2
                                           + pow(log( dabs( ( _Del - _Smp ) / 2 / _Q2 ) ), 2) / 2
                                           - pow(log( dabs( ( _Del + _Spm ) / 2 / _Q2 ) ), 2) / 2
                                           - pow(log( dabs( ( _Del + _Smp ) / 2 / _Q2 ) ), 2) / 2
                                           - dilog( ( _Del - _Spm ) / 2 / _Del ) - dilog( ( _Del - _Smp ) / 2 / _Del )
                                           + dilog( ( _Del + _Spm ) / 2 / _Del ) + dilog( ( _Del + _Smp ) / 2 / _Del ) );

    // Soft and virtual contributions
    _S1 = 2 + _Spp / _Del * ( _Del * I1 + dilog( 2 * _Del / ( _Del - _Spp ) ) - dilog( 2 * _Del / ( _Del + _Spp ) ) )
          + log( _Del2 / _m22 / _Q2 ) * ( - 2 + _Spp * I1 );
    _V1 = CRm + ( _Sminus * _Spp - 2 * _Splus * _m1 * _m2 ) / ( _Splus * _Spp - 2 * _Sminus * _m1 * _m2 ) * Cplus;

    // Do not compute _R exactly at one to avoid numerical problems
    _R1 = _R(1 - eps15);
  }

  //_________________________________________________________________________________
  double Cm11ns::_R(double const& x) const
  {
    const double s1h   = ( 1 - x ) * ( ( _Del - _Spm ) * x + _Del + _Spm ) / 2 / x;
    const double s1h2  = s1h * s1h;
    const double Delp  = DeltaFun(_m12, s1h + _m22, - _Q2);
    const double Delp2 = Delp * Delp;
    const double Lxi   = log( ( _Spp + s1h - Delp ) / ( _Spp + s1h + Delp ) );
    const double Ixi   = ( ( s1h + 2 * _m22 ) / s1h2 + ( s1h + _m22 ) / Delp / s1h2 * _Spp * Lxi );
    const double N1    = ( _Splus * _Spp - 2 * _m1 * _m2 * _Sminus ) / 2 / _Del;
    const double f1hat = 8 / Delp2 * ( - _Del2 * ( _Splus * _Spp - 2 * _m1 * _m2 * _Sminus ) * Ixi + 2 * _m1 * _m2 * _Sminus
                                       * ( 1 / s1h * ( Delp2 + 4 * _m22 * _Spm ) + 2 * _Spm - _Smp + (_Spp + s1h) / 2 + ( s1h + _m22 ) / Delp / s1h
                                           * ( Delp2 + 2 * _Spm * _Spp + ( _m22 + _Q2 ) * s1h ) * Lxi )
                                       + _Splus * ( ( - _m22 * _Spp ) / ( ( s1h + _m22 ) * s1h ) * ( _Del2 + 4 * _m22 * _Spm) - 1. / 4 / ( s1h + _m22 )
                                                    * ( 3 * pow(_Spp, 2) * _Smp + 4 * _m22 * (10 * _Spp * _Spm - _Spm * _Smp - _m12 * _Spp )
                                                        + s1h * ( - 7 * _Spp * _Smp + 18 * _Del2 - 4 * _m12 * ( 7 * _Q2 - 4 * _m22 + 7 * _m12 ) )
                                                        + 3 * s1h2 * ( _Spm - 2 * _m12 ) - pow(s1h, 3) ) + ( s1h + _m22 ) / 2 / Delp
                                                    * ( - 2 / s1h * _Spp * ( _Del2 + 2 * _Spm * _Spp )
                                                        + ( 4 * _m12 * _m22 - 7 * _Spm * _Spp ) - 4 * _Spm * s1h - pow(s1h, 2) ) * Lxi ) );

    return ( 1 - x ) * s1h * f1hat / N1 / 8 / ( s1h + _m22 );
  }

  //_________________________________________________________________________________
  double Cm11ns::Regular(double const& x) const
  {
    return _fact1 * ( _R(x) - _R1 ) / ( 1 - x );
  }

  //_________________________________________________________________________________
  double Cm11ns::Singular(double const& x) const
  {
    return _fact1 *_R1 / ( 1 - x );
  }

  //_________________________________________________________________________________
  double Cm11ns::Local(double const& x) const
  {
    return _fact1 * ( _S1 + _V1 + _R1 * log(1 - x) );
  }

  //_________________________________________________________________________________
  Cm21ns::Cm21ns(double const& m1, double const& m2, double const& Q, double const& Splus, double const& Sminus):
    Expression(2 * pow(Q, 2) / ( pow(Q, 2) + pow(m2, 2) - pow(m1, 2) + DeltaFun(pow(m1, 2), pow(m2, 2), - pow(Q, 2)))),
    _m1(m1),
    _m2(m2),
    _Splus(Splus),
    _Sminus(Sminus),
    _m12(_m1 * _m1),
    _m22(_m2 * _m2),
    _Q2(Q * Q),
    _Del(DeltaFun(_m12, _m22, - _Q2)),
    _Del2(_Del * _Del),
    _Spp(_Q2 + _m22 + _m12),
    _Spm(_Q2 + _m22 - _m12),
    _Smp(_Q2 - _m22 + _m12),
    _fact2(2 * CF * _Del / _Q2)
  {
    const double I1    = log( ( _Spp + _Del ) / ( _Spp - _Del ) ) / _Del;
    const double Cplus = 2 * m1 * m2 * I1;
    const double C1m   = - ( _Spm * I1 + log( _m12 / _m22 ) ) / _Q2;
    const double C1p   = - ( _Smp * I1 - log( _m12 / _m22 ) ) / _Q2;
    const double CRm   = ( _Del2 / 2 / _Q2 + _Spp * ( 1 + log( _Q2 / _Del ) ) ) * I1
                         + ( _m22 - _m12 ) / 2 / _Q2 * log( _m12 / _m22 )
                         - log( _Q2 / _m12 ) - log( _Q2 / _m22 ) - 4
                         + _Spp / _Del * ( + pow(log( dabs( ( _Del - _Spm ) / 2 / _Q2 ) ), 2) / 2
                                           + pow(log( dabs( ( _Del - _Smp ) / 2 / _Q2 ) ), 2) / 2
                                           - pow(log( dabs( ( _Del + _Spm ) / 2 / _Q2 ) ), 2) / 2
                                           - pow(log( dabs( ( _Del + _Smp ) / 2 / _Q2 ) ), 2) / 2
                                           - dilog( ( _Del - _Spm ) / 2 / _Del ) - dilog( ( _Del - _Smp ) / 2 / _Del )
                                           + dilog( ( _Del + _Spm ) / 2 / _Del ) + dilog( ( _Del + _Smp ) / 2 / _Del ) );

    // Soft and virtual contributions
    _S2 = 2 + _Spp / _Del * ( _Del * I1 + dilog( 2 * _Del / ( _Del - _Spp ) ) - dilog( 2 * _Del / ( _Del + _Spp ) ) )
          + log( _Del2 / _m22 / _Q2 ) * ( - 2 + _Spp * I1 );
    _V2 = CRm + ( _m12 * C1p + _m22 * C1m ) / 2 + _Sminus / _Splus * ( Cplus + m1 * m2 / 2 * ( C1p + C1m ) );

    // Do not compute _R exactly at one to avoid numerical problems
    _R1 = _R(1 - eps15);
  }

  //_________________________________________________________________________________
  double Cm21ns::_R(double const& x) const
  {
    const double s1h   = ( 1 - x ) * ( ( _Del - _Spm ) * x + _Del + _Spm ) / 2 / x;
    const double s1h2  = s1h * s1h;
    const double Delp  = DeltaFun(_m12, s1h + _m22, - _Q2);
    const double Delp2 = Delp * Delp;
    const double Lxi   = log( ( _Spp + s1h - Delp ) / ( _Spp + s1h + Delp ) );
    const double Ixi   = ( ( s1h + 2 * _m22 ) / s1h2 + ( s1h + _m22 ) / Delp / s1h2 * _Spp * Lxi );
    const double N2    = 2 * _Splus * _Del / Delp2;
    const double f2hat = 16 / pow(Delp, 4) * ( - 2 * pow(_Del, 4) * _Splus * Ixi + 2 * _m1 * _m2 * _Sminus *
                                               ( ( ( s1h + _m22 ) / Delp ) * ( Delp2 - 6 * _m12 * _Q2 ) * Lxi
                                                 - Delp2 * ( s1h + _Spp ) / 2 / ( s1h + _m22 ) + ( 2 * Delp2 - 3 * _Q2 * ( s1h + _Spp) ) )
                                               + _Splus * ( - 2 * ( _Del2 - 6 * _m12 * _Q2 ) * ( s1h + _m22 ) - 2 * ( _m12 + _m22 ) * s1h2
                                                            - 9 * _m22 * pow(_Spm, 2) + _Del2 * ( 2 * _Spp - _m22 ) + 2 * s1h * ( 2 * _Del2 + ( _m12 - 5 * _m22 ) * _Spm )
                                                            + ( Delp2 - 6 * _Q2 * ( _m22 + s1h ) ) * _Spp * ( s1h + _Spp ) / 2 / ( s1h + _m22 )
                                                            - 2 * _Del2 / s1h * ( _Del2 + 2 * ( 2 * _m22 + s1h ) * _Spm ) + ( s1h + _m22 ) / Delp
                                                            * ( - 2 / s1h * _Del2 * ( _Del2 + 2 * _Spm * _Spp ) - 2 * s1h * ( _Del2 - 6 * _m12 * _Q2 )
                                                                - ( Delp2 - 18 * _m12 * _Q2 ) * _Spp - 2 * _Del2 * ( _Spp + 2 * _Spm) ) * Lxi ) );

    return ( 1 - x ) * s1h * f2hat / N2 / 8 / ( s1h + _m22 );
  }

  //_________________________________________________________________________________
  double Cm21ns::Regular(double const& x) const
  {
    return _fact2 * ( _R(x) - _R1 ) / ( 1 - x );
  }

  //_________________________________________________________________________________
  double Cm21ns::Singular(double const& x) const
  {
    return _fact2 * _R1 / ( 1 - x );
  }

  //_________________________________________________________________________________
  double Cm21ns::Local(double const& x) const
  {
    return _fact2 * ( _S2 + _V2 + _R1 * log(1 - x) );
  }

  //_________________________________________________________________________________
  Cm31ns::Cm31ns(double const& m1, double const& m2, double const& Q, double const& Rplus, double const& Rminus):
    Expression(2 * pow(Q, 2) / ( pow(Q, 2) + pow(m2, 2) - pow(m1, 2) + DeltaFun(pow(m1, 2), pow(m2, 2), - pow(Q, 2)))),
    _m1(m1),
    _m2(m2),
    _Rplus(Rplus),
    _Rminus(Rminus),
    _m12(_m1 * _m1),
    _m22(_m2 * _m2),
    _Q2(Q * Q),
    _Del(DeltaFun(_m12, _m22, - _Q2)),
    _Del2(_Del * _Del),
    _Spp(_Q2 + _m22 + _m12),
    _Spm(_Q2 + _m22 - _m12),
    _Smp(_Q2 - _m22 + _m12),
    _fact3(2 * CF)
  {
    const double I1    = log( ( _Spp + _Del ) / ( _Spp - _Del ) ) / _Del;
    const double Cplus = 2 * m1 * m2 * I1;
    const double CRm   = ( _Del2 / 2 / _Q2 + _Spp * ( 1 + log( _Q2 / _Del ) ) ) * I1
                         + ( _m22 - _m12 ) / 2 / _Q2 * log( _m12 / _m22 )
                         - log( _Q2 / _m12 ) - log( _Q2 / _m22 ) - 4
                         + _Spp / _Del * ( + pow(log( dabs( ( _Del - _Spm ) / 2 / _Q2 ) ), 2) / 2
                                           + pow(log( dabs( ( _Del - _Smp ) / 2 / _Q2 ) ), 2) / 2
                                           - pow(log( dabs( ( _Del + _Spm ) / 2 / _Q2 ) ), 2) / 2
                                           - pow(log( dabs( ( _Del + _Smp ) / 2 / _Q2 ) ), 2) / 2
                                           - dilog( ( _Del - _Spm ) / 2 / _Del ) - dilog( ( _Del - _Smp ) / 2 / _Del )
                                           + dilog( ( _Del + _Spm ) / 2 / _Del ) + dilog( ( _Del + _Smp ) / 2 / _Del ) );

    // Soft and virtual contributions
    _S3 = 2 + _Spp / _Del * ( _Del * I1 + dilog( 2 * _Del / ( _Del - _Spp ) ) - dilog( 2 * _Del / ( _Del + _Spp ) ) )
          + log( _Del2 / _m22 / _Q2 ) * ( - 2 + _Spp * I1 );
    _V3 = CRm + _Rminus / _Rplus * Cplus;

    // Do not compute _R exactly at one to avoid numerical problems
    _R1 = _R(1 - eps15);
  }

  //_________________________________________________________________________________
  double Cm31ns::_R(double const& x) const
  {
    const double s1h   = ( 1 - x ) * ( ( _Del - _Spm ) * x + _Del + _Spm ) / 2 / x;
    const double s1h2  = s1h * s1h;
    const double Delp  = DeltaFun(_m12, s1h + _m22, - _Q2);
    const double Delp2 = Delp * Delp;
    const double Lxi   = log( ( _Spp + s1h - Delp ) / ( _Spp + s1h + Delp ) );
    const double Ixi   = ( ( s1h + 2 * _m22 ) / s1h2 + ( s1h + _m22 ) / Delp / s1h2 * _Spp * Lxi );
    const double N3    = 2 * _Rplus / Delp;
    const double f3hat = 16 / Delp2 * ( - 2 * _Del2 * _Rplus * Ixi + 2 * _m1 * _m2 * _Rminus
                                        * ( 1 - _Smp / s1h + ( s1h + _m22 ) * ( s1h + _Spm ) / Delp / s1h * Lxi )
                                        + _Rplus * ( _Smp - 3 * _Spm - 2 / s1h * ( _Del2 + 2 * _m22 * _Spm )
                                                     - ( s1h - _Smp ) * ( s1h + _Spp ) / 2 / ( s1h + _m22 )
                                                     + ( s1h + _m22 ) / Delp / s1h
                                                     * ( - s1h2 + 4 * ( _m12 * _Smp - _Del2 ) - 3 * s1h * _Spm ) * Lxi ) );

    return ( 1 - x ) * s1h * f3hat / N3 / 8 / ( s1h + _m22 );
  }

  //_________________________________________________________________________________
  double Cm31ns::Regular(double const& x) const
  {
    return _fact3 * ( _R(x) - _R1 ) / ( 1 - x );
  }

  //_________________________________________________________________________________
  double Cm31ns::Singular(double const& x) const
  {
    return _fact3 * _R1 / ( 1 - x );
  }

  //_________________________________________________________________________________
  double Cm31ns::Local(double const& x) const
  {
    return _fact3 * ( _S3 + _V3 + _R1 * log(1 - x) );
  }

  //_________________________________________________________________________________
  CmL1ns::CmL1ns(double const& m1, double const& m2, double const& Q, double const& Splus, double const& Sminus):
    Expression(2 * pow(Q, 2) / ( pow(Q, 2) + pow(m2, 2) - pow(m1, 2) + DeltaFun(pow(m1, 2), pow(m2, 2), - pow(Q, 2)))),
    _C1(Cm11ns{m1, m2, Q, Splus, Sminus}),
    _C2(Cm21ns{m1, m2, Q, Splus, Sminus}),
    _factL(2 * Splus / ( Splus + Sminus ))
  {
  }

  //_________________________________________________________________________________
  double CmL1ns::Regular(double const& x) const
  {
    return _factL * ( _C2.Regular(x) - _C1.Regular(x) );
  }

  //_________________________________________________________________________________
  double CmL1ns::Singular(double const& x) const
  {
    return _factL * ( _C2.Singular(x) - _C1.Singular(x) );
  }

  //_________________________________________________________________________________
  double CmL1ns::Local(double const& x) const
  {
    return _factL * ( _C2.Local(x) - _C1.Local(x) );
  }

  //_________________________________________________________________________________
  Cm21qCC::Cm21qCC(double const& lambda):
    Expression(),
    _lambda(lambda)
  {
  }
  double Cm21qCC::Regular(double const& z) const
  {
    return 2 * CF * ( - ( 1 + pow(z, 2) ) * log(z) / ( 1 - z )
                      + ( 2 + log(_lambda) - 2 * log(1 - z) + log(1 - _lambda * z) ) * ( 1 + z )
                      + 1 / _lambda );
  }
  double Cm21qCC::Singular(double const& z) const
  {
    return 2 * CF * ( 2 * ( 2 * log(1 - z) - log(1 - _lambda * z) ) / ( 1 - z )
                      + 2 * ( - 1 - log(_lambda) ) / ( 1 - z )
                      + ( 2 * pow(_lambda, 2) - _lambda - 1 ) / _lambda / ( 1 - _lambda * z )
                      + ( 1 - z ) / pow(1 - _lambda * z, 2) / 2 );
  }
  double Cm21qCC::Local(double const& z) const
  {
    const double KA = ( 1 - _lambda ) * log(1 - _lambda) / _lambda;
    const double Rx = - dilog(1 / ( 1 - _lambda )) + dilog(( 1 - z * _lambda ) / ( 1 - _lambda ))
                      + log(_lambda * ( 1 - z ) / ( 1 - _lambda )) * log(1 - _lambda * z);
    return 2 * CF * ( - 4 - 1. / 2. / _lambda - 2 * zeta2 - ( 1 + 3 * _lambda ) * KA / 2 / _lambda - 3 * log(_lambda) / 2
                      + 2 * pow(log(1 - z), 2) - 2 * Rx + 2 * ( - 1 - log(_lambda) ) * log(1 - z)
                      + ( 2 * pow(_lambda, 2) - _lambda - 1 ) * log(1 - _lambda * z) / pow(_lambda, 2)
                      + log(1 - _lambda * z) / 2 / pow(_lambda, 2)
                      + ( 1 - _lambda ) * z / 2 / _lambda / ( 1 - _lambda * z ) );
  }

  //_________________________________________________________________________________
  Cm21gCC::Cm21gCC(double const& lambda):
    Expression(),
    _lambda(lambda)
  {
  }
  double Cm21gCC::Regular(double const& z) const
  {
    return 4 * TR * ( ( pow(z, 2) + pow(1 - z, 2) ) * ( log(( 1 - z ) / z) - log(1 - _lambda) / 2 - log(_lambda) / 2 )
                      + 8 * z * ( 1 - z ) - 1
                      + ( 1 - _lambda ) * ( - 6 * ( 1 + 2 * _lambda ) * z * ( 1 - z ) + 1 / ( 1 - _lambda * z )
                                            + 6 * _lambda * z * ( 1 - 2 * _lambda * z ) * log(( 1 - _lambda * z ) / ( 1 - _lambda ) / z) ) );
  }

  //_________________________________________________________________________________
  CmL1qCC::CmL1qCC(double const& lambda):
    Expression(),
    _lambda(lambda)
  {
  }
  double CmL1qCC::Regular(double const& z) const
  {
    return 2 * CF * ( 1 - _lambda ) * ( - ( 1 + pow(z, 2) ) * log(z) / ( 1 - z )
                                        + ( log(_lambda) - 2 * log(1 - z) + log(1 - _lambda * z) ) * ( 1 + z )
                                        + 3 )
           + 2 * CF * ( 1 + _lambda ) * z;
  }
  double CmL1qCC::Singular(double const& z) const
  {
    return 2 * CF * ( 1 - _lambda ) * ( 2 * ( 2 * log(1 - z) - log(1 - _lambda * z) ) / ( 1 - z )
                                        + 2 * ( - 1 - log(_lambda) ) / ( 1 - z )
                                        - 2 / ( 1 - _lambda * z )
                                        + ( 1 - z ) / pow(1 - _lambda * z, 2) / 2 );
  }
  double CmL1qCC::Local(double const& z) const
  {
    const double KA = ( 1 - _lambda ) * log(1 - _lambda) / _lambda;
    const double Rx = - dilog(1 / ( 1 - _lambda )) + dilog(( 1 - z * _lambda ) / ( 1 - _lambda ))
                      + log(_lambda * ( 1 - z ) / ( 1 - _lambda )) * log(1 - _lambda * z);
    return 2 * CF * _lambda * KA
           + 2 * CF * ( 1 - _lambda ) * ( 2 * pow(log(1 - z), 2) - 2 * Rx + 2 * ( - 1 - log(_lambda) ) * log(1 - z)
                                          - 2 * log(1 - _lambda * z) / _lambda
                                          + log(1 - _lambda * z) / 2 / pow(_lambda, 2)
                                          + ( 1 - _lambda ) * z / 2 / _lambda / ( 1 - _lambda * z ) );
  }

  //_________________________________________________________________________________
  CmL1gCC::CmL1gCC(double const& lambda):
    Expression(),
    _lambda(lambda)
  {
  }
  double CmL1gCC::Regular(double const& z) const
  {
    return 4 * TR * ( ( 1 - _lambda ) * ( pow(z, 2) + pow(1 - z, 2) ) * ( log(( 1 - z ) / z) - log(1 - _lambda) / 2 - log(_lambda) / 2 )
                      + 4 * ( 2 - _lambda ) * z * ( 1 - z )
                      + ( 1 - _lambda ) * ( - 2 * ( 3 + 4 * _lambda ) * z * ( 1 - z )
                                            + 4 * _lambda * z * ( 1 - 2 * _lambda * z ) * log(( 1 - _lambda * z ) / ( 1 - _lambda ) / z) ) );
  }

  //_________________________________________________________________________________
  Cm31qCC::Cm31qCC(double const& lambda):
    Expression(),
    _lambda(lambda)
  {
  }
  double Cm31qCC::Regular(double const& z) const
  {
    return 2 * CF * ( - ( 1 + pow(z, 2) ) * log(z) / ( 1 - z )
                      + ( 1 + log(_lambda) - 2 * log(1 - z) + log(1 - _lambda * z) ) * ( 1 + z )
                      + 1 / _lambda );
  }
  double Cm31qCC::Singular(double const& z) const
  {
    return 2 * CF * ( 2 * ( 2 * log(1 - z) - log(1 - _lambda * z) ) / ( 1 - z )
                      + 2 * ( - 1 - log(_lambda) ) / ( 1 - z )
                      + ( _lambda - 1 ) / _lambda / ( 1 - _lambda * z )
                      + ( 1 - z ) / pow(1 - _lambda * z, 2) / 2 );
  }
  double Cm31qCC::Local(double const& z) const
  {
    const double KA = ( 1 - _lambda ) * log(1 - _lambda) / _lambda;
    const double Rx = - dilog(1 / ( 1 - _lambda )) + dilog(( 1 - z * _lambda ) / ( 1 - _lambda ))
                      + log(_lambda * ( 1 - z ) / ( 1 - _lambda )) * log(1 - _lambda * z);
    return 2 * CF * ( - 4 - 1. / 2. / _lambda - 2 * zeta2 - ( 1 + 3 * _lambda ) * KA / 2 / _lambda - 3 * log(_lambda) / 2
                      + 2 * pow(log(1 - z), 2) - 2 * Rx + 2 * ( - 1 - log(_lambda) ) * log(1 - z)
                      + ( _lambda - 1 ) * log(1 - _lambda * z) / pow(_lambda, 2)
                      + log(1 - _lambda * z) / 2 / pow(_lambda, 2)
                      + ( 1 - _lambda ) * z / 2 / _lambda / ( 1 - _lambda * z ) );
  }

  //_________________________________________________________________________________
  Cm31gCC::Cm31gCC(double const& lambda):
    Expression(),
    _lambda(lambda)
  {
  }
  double Cm31gCC::Regular(double const& z) const
  {
    return 4 * TR * ( ( pow(z, 2) + pow(1 - z, 2) ) * ( log(( 1 - z ) / ( 1 - _lambda * z )) + log(1 - _lambda) / 2 - log(_lambda) / 2 )
                      + ( 1 - _lambda ) * ( 2 * z * ( 1 - z ) - 2 * z * ( 1 - ( 1 + _lambda ) * z ) * log(( 1 - _lambda * z ) / ( 1 - _lambda ) / z) ) );
  }
}
