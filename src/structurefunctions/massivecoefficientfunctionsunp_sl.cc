//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/massivecoefficientfunctionsunp_sl.h"
#include "apfel/betaqcd.h"
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
    return _eta * 16 * M_PI * xi * clnlobarq_(&wr, &xi) / z;
  }

  //_________________________________________________________________________________
  // Helper functions for the approximated O(as^3) coefficient
  // functions. The correspond to Eqs. (3.10) and (3.11) of
  // https://arxiv.org/pdf/1205.5727.
  //_________________________________________________________________________________
  double c0Th(double const& xi)
  {
    const double y  = sqrt(1 + 4 / xi);
    const double L1 = log(1 + xi / 2);
    const double L2 = log(2 + xi / 2);
    const double L3 = log(sqrt(xi) * ( y - 1 ) / 2);
    const double xip2 = 2 + xi;
    const double xip4 = 4 + xi;
    const double Li2 = dilog(-2 / xip2);
    return CA * ( 50 - Pi2 + 12 * L3 / y + 4 * L3 * L3 + L1 * L1 + 6 * L2
                  - 4 * L2 * L2 + 2 * Li2 + 48 / xip2 - 4 * L2 / xip2
                  + 64 * L2 / xip2 / xip2 - 128 * L2 / ( xip2 * xip2 * xip4 )
                  - 160 / xip2 / xip4 - 64 * L2 / xip2 / xip4
                  + 128 / ( xip2 * xip4 * xip4 ) - 12 * ( 4 + zeta2 ) / xip4
                  - 8 * L3 * L3 / xip4 + 64 / xip4 / xip4 )
           + CF * ( -18. - 2. / 3. * Pi2 - 24. * L3 / y - 8. * L3 * L3
                    + 2. * L1 * L1 - 6. * L2 + 4. * Li2 - 48. / xip2
                    + 8. * L2 / xip2 + 360. / xip2 / xip4
                    + 128. * L2 / xip2 / xip4 - 544. / (xip2 * xip4 * xip4)
                    + 48. * L3 * L3 / xip4 - 8. * L1 * L1 / xip4
                    + (44. + 40. * zeta2) / xip4 - 120. * L2 / xip2 / xip2
                    + 256. * L2 / (xip2 * xip2 * xip4) - 16. * Li2 / xip4
                    - 272. / xip4 / xip4 );
  }
  double c0Thbar(double const& xi)
  {
    return 4 * CA * ( 2 + log(1 + xi / 4) ) - 4 * TR / 3;
  }
  double ThExpansiong3(int const& nf, double const& xi, double const& z)
  {
    const double beta = sqrt(1 - 4 * z / xi / ( 1 - z ));
    const double l  = log(beta);
    const double l2 = l * l;
    const double l3 = l * l2;
    const double l4 = l * l3;
    const double Lm  = - log(xi);
    const double Lm2 = Lm * Lm;
    const double ln2  = log(2);
    const double ln22 = ln2 * ln2;
    const double ln23 = ln2 * ln22;
    const double c_log4 = 128 * CA * CA;
    const double c_log3 = ( 768 * ln2 - 6464 / 9. ) * CA * CA + 128. / 9. * CA * nf + 128 * CA * CA * Lm;
    const double c_log2 = ( 1728 * ln22 - 3232 * ln2 - 208 / 3. * Pi2 + 15520 / 9. ) * CA * CA
                          + ( 64. * ln2 - 640 / 9. ) * CA * nf + 16 * CA * c0Th(xi) + 32 * CA * ( CF - CA / 2 ) * Pi2 / beta
                          - ( ( - 512 * ln2 + 1136 / 3. ) * CA * CA - 32 / 3. * CA * nf + 16 * CA * c0Thbar(xi) ) * Lm
                          + 32 * CA * CA * Lm2;
    const double c_log_const = ( 1728 * ln23 - 4848 * ln22 + 15520 / 3. * ln2 - 208 * Pi2 * ln2
                                 + 936 * zeta3 + 608. / 3. * Pi2 - 88856. / 27.) * CA * CA
                               + ( 96 * ln22 - 640 / 3. * ln2 - 16 / 3. * Pi2 + 4592 / 27. ) * CA * nf
                               - 32 * CF * ( CF - CA / 2 ) * Pi2 + ( 48 * ln2 - 40 ) * CA * c0Th(xi);
    const double c_log_fracbeta = ( ( - 92 / 3. + 32 * ln2 ) * CA + 8 / 3. * nf) * ( CF - CA / 2 ) * Pi2;
    const double c_log_Lm = - ( ( - 672 * ln22 + 976 * ln2 + 104 / 3. * Pi2 - 4160 / 9. ) * CA * CA
                                + ( - 32 * ln2 + 320 / 9. ) * CA * nf + ( 48  * ln2 - 40 ) * CA * c0Thbar(xi)
                                - 8 * CA * c0Th(xi) - 16 * CA * ( CF - CA / 2 ) * Pi2 / beta );
    const double c_log_Lm2 = ( 64 * ln2 - 44 / 3. ) * CA * CA + 8 / 3. * CA * nf - 8 * CA * c0Thbar(xi);
    const double c_log = c_log_const + c_log_fracbeta / beta + c_log_Lm * Lm + c_log_Lm2 * Lm2;
    const double c_fracbeta = ( ( 8 * ln22 - 68 / 3. * ln2 + 8 / 3. * Pi2 - 658 / 9. ) * CA
                                + ( 8 / 3. * ln2 - 20 / 9.) * nf + 2 * c0Th(xi)
                                + ( 26 / 3. * CA + 4 / 3. * nf - 2 * c0Thbar(xi) ) * Lm ) * ( CF - CA / 2 ) * Pi2;
    const double c_fracbeta2 = 4 / 3. * ( CF - CA / 2 ) * ( CF - CA / 2 ) * Pi2 * Pi2;
    return c_log4 * l4 + c_log3 * l3 + c_log2 * l2 + c_log * l + c_fracbeta / beta + c_fracbeta2 / beta / beta;
  }
  double ThExpansiong3const(double const& xi)
  {
    return pow(c0Th(xi) + 36 * CA * pow(log(2), 2) - 60 * CA * log(2) - log(xi) * ( 8 * CA * log(2) - c0Thbar(xi) ), 2);
  }

  //_________________________________________________________________________________
  CmTh23gNC::CmTh23gNC(int const& nf, double const& eta):
    Expression(),
    _nf(nf),
    _eta(eta),
    _c21g(Cm21gNC{eta})
  {
  }
  double CmTh23gNC::Regular(double const& x) const
  {
    if (x >= 1)
      return 0;
    const double eta = this->_eta;
    const double z   = eta * x;
    const double xi  = 4 * eta / ( 1 - eta );
    return _c21g.Regular(x) * ( ThExpansiong3(_nf, xi, z) + ThExpansiong3const(xi) );
  }

  //_________________________________________________________________________________
  CmThL3gNC::CmThL3gNC(int const& nf, double const& eta):
    Expression(),
    _nf(nf),
    _eta(eta),
    _cL1g(CmL1gNC{eta})
  {
  }
  double CmThL3gNC::Regular(double const& x) const
  {
    if (x >= 1)
      return 0;
    const double eta = this->_eta;
    const double z   = eta * x;
    const double xi  = 4 * eta / ( 1 - eta );
    return _cL1g.Regular(x) * ( ThExpansiong3(_nf, xi, z) + ThExpansiong3const(xi) );
  }

  //_________________________________________________________________________________
  Cmsx23gNC::Cmsx23gNC(int const& nf, double const& eta):
    Expression(),
    _nf(nf),
    _eta(eta)
  {
  }
  double Cmsx23gNC::Regular(double const& x) const
  {
    if (x >= 1)
      return 0;
    const double eta  = this->_eta;
    const double xi   = 4 * eta / ( 1 - eta );
    const double Lmu  = - log(xi);
    const double Lmu2 = Lmu * Lmu;
    const double a11  = CA;
    const double a21  = _nf * ( 26 * CF - 23 * CA ) / 36;
    const double a10  = - ( 11 * CA + 2 * _nf * ( 1 - 2 * CF / CA ) ) / 12;
    const double bt0  = beta0qcd(_nf);
    const double z    = sqrt(1 / ( 1 + 4 / xi ));
    const double L    = log(( 1 + z ) / ( 1 - z ));

    // Call polylogs
    double wx = z;
    int nw = 3;
    int n1 = -1;
    int n2 = 1;
    int sz = n2 - n1 + 1;
    double *Hr1 = new double[sz];
    double *Hr2 = new double[sz * sz];
    double *Hr3 = new double[sz * sz * sz];
    double *Hr4 = new double[sz * sz * sz * sz];
    double *Hr5 = new double[sz * sz * sz * sz * sz];
    apf_hplog_(&wx, &nw, Hr1, Hr2, Hr3, Hr4, Hr5, &n1, &n2);

    // weight 2
    const double Hm1m1 = Hr2[0];
    const double H1m1  = Hr2[2];
    const double Hm11  = Hr2[6];
    const double H11   = Hr2[8];

    // weight 3
    const double Hm1m1m1 = Hr3[0];
    const double H1m1m1  = Hr3[2];
    const double Hm11m1  = Hr3[6];
    const double H11m1   = Hr3[8];
    const double Hm1m11  = Hr3[18];
    const double H1m11   = Hr3[20];
    const double Hm111   = Hr3[24];
    const double H111    = Hr3[26];

    // Deallocate pointers
    delete[] Hr1;
    delete[] Hr2;
    delete[] Hr3;
    delete[] Hr4;
    delete[] Hr5;

    const double Hmp  = H11 + H1m1 - Hm11 - Hm1m1;
    const double Hmpm = H111 - H11m1 + H1m11 - H1m1m1 - Hm111 + Hm11m1 - Hm1m11 + Hm1m1m1;
    const double I    = 4 * z * Hmp;
    const double J    = 4 * z * L;
    const double K    = 4 * z * Hmpm;
    const double Lnxi = log(1 + xi / 4);
    //const double ln2  = log(2);
    const double ll   = CA * CA * ( - 1472 / 27. - 8 / 3. * K * ( - 1 + 1 / xi ) + 8 / 27. * J * ( - 71 + 92 / xi )
                                    + I * ( 8 / 3. * Lnxi * ( - 1 + 1 / xi ) + 8 / 9. * ( - 13 + 10 / xi ) )
                                    + ( - 160 / 9. + 16 / 3. * I * ( - 1 + 1 / xi ) + 8 / 9. * J * ( - 13 + 10 / xi ) ) * Lmu
                                    + ( - 16 / 3. + 8 / 3. * J * ( - 1 + 1 / xi ) ) * Lmu2 ) * log(x) / x;
    const double nllc =  ( a21 * ( 160 / 9. - 16 / 3. * I * ( - 1 + 1 / xi ) - 8 / 9. * J * ( - 13 + 10 * 1 / xi ) )
                           + a10 * a11 * ( 2944 / 27. + 16 / 3. * K * ( - 1 + 1 / xi ) - 16 / 27. * J * ( - 71 + 92 * 1 / xi )
                                           + I * ( - 16 / 3. * Lnxi * ( - 1 + 1 / xi ) - 16 / 9. * ( - 13 + 10 * 1 / xi ) ) )
                           + a11 * bt0 * ( - 1472 / 27. - 8 / 3. * K * ( - 1 + 1 / xi ) + 8 / 27. * J * ( - 71 + 92 * 1 / xi )
                                           + I * ( 8 / 3. * Lnxi * ( - 1 + 1 / xi ) + 8 / 9. * ( - 13 + 10 * 1 / xi ) ) )
                           + ( a21 * ( 32 / 3. - 16 / 3. * J * ( - 1 + 1 / xi ) )
                               + a10 * a11 * ( 320 / 9. - 32 / 3. * I * ( - 1 + 1 / xi ) - 16 / 9. * J * ( - 13 + 10 * 1 / xi ) )
                               + a11 * bt0 * ( - 160 / 9. + 16 / 3. * I * ( - 1 + 1 / xi ) + 8 / 9. * J * ( - 13 + 10 * 1 / xi ) ) ) * Lmu
                           + ( a10 * a11 * ( 32 / 3. - 16 / 3. * J * ( - 1 + 1 / xi ) )
                               + a11 * bt0 * ( - 16 / 3. + 8 / 3. * J * ( - 1 + 1 / xi ) ) ) * Lmu2 ) / x;
    /*
    // Variation of nllc
    const double nllv = ( a10 * a11 * ( 2944 / 27. + 16 / 3. * K * ( - 1 + 1 / xi ) - 16 / 27. * J * ( - 71 + 92 * 1 / xi )
                                        + I * ( - 16 / 3. * Lnxi * ( - 1 + 1 / xi ) - 16 / 9. * ( - 13 + 10 * 1 / xi ) ) )
                          + ( a10 * a11 * ( 32 / 3. - 16 / 3. * J * ( - 1 + 1 / xi ) )
                              + a11 * bt0 * ( - 16 / 3. + 8 / 3. * J * ( - 1 + 1 / xi ) ) ) * Lmu2
                          + a11 * bt0 * ( - 1472 / 27. - 8 / 3. * K * ( - 1 + 1 / xi ) - 640 * ln2 / 9. + 140 * zeta3 / 3
                                          + I * ( 8 / 3. * Lnxi * ( - 1 + 1 / xi ) + 8 / 9. * ( - 13 + 10 * 1 / xi )
                                                  + 64 / 3. * ( - 1 + 1 / xi ) * ln2 - 14 * ( - 1 + 1 / xi ) * zeta3 )
                                          + J * ( 8 / 27. * ( - 71 + 92 * 1 / xi ) + 32 / 9. * ( - 13 + 10 * 1 / xi ) * ln2
                                                  - 7 / 3. * ( - 13 + 10 * 1 / xi ) * zeta3) )
                          + Lmu * ( a10 * a11 * ( 320 / 9. - 32 / 3. * I * ( - 1 + 1 / xi ) - 16 / 9. * J * ( - 13 + 10 * 1 / xi ) )
                                    + a11 * bt0 * ( - 160 / 9. + 16 / 3. * I * ( - 1 + 1 / xi ) - 128 * ln2 / 3 + 28 * zeta3
                                                    + J * ( 8 / 9. * ( - 13 + 10 * 1 / xi ) + 64 / 3. * ( - 1 + 1 / xi ) * ln2
                                                            - 14 * ( - 1 + 1 / xi ) * zeta3) ) ) ) / x;
    */
    return ll + nllc;
  }

  //_________________________________________________________________________________
  Cmsx23psNC::Cmsx23psNC(int const& nf, double const& eta):
    Expression(),
    _c23g(Cmsx23gNC{nf, eta})
  {
  }
  double Cmsx23psNC::Regular(double const& x) const
  {
    return CF / CA * _c23g.Regular(x);
  }

  //_________________________________________________________________________________
  CmsxL3gNC::CmsxL3gNC(int const& nf, double const& eta):
    Expression(),
    _nf(nf),
    _eta(eta)
  {
  }
  double CmsxL3gNC::Regular(double const& x) const
  {
    if (x >= 1)
      return 0;
    const double eta  = this->_eta;
    const double xi   = 4 * eta / ( 1 - eta );
    const double Lmu  = - log(xi);
    const double Lmu2 = Lmu * Lmu;
    const double a11  = CA;
    const double a21  = _nf * ( 26 * CF - 23 * CA ) / 36;
    const double a10  = - ( 11 * CA + 2 * _nf * ( 1 - 2 * CF / CA ) ) / 12;
    const double bt0  = beta0qcd(_nf);
    const double z    = sqrt(1 / ( 1 + 4 / xi ));
    const double L    = log(( 1 + z ) / ( 1 - z ));

    // Call polylogs
    double wx = z;
    int nw = 3;
    int n1 = -1;
    int n2 = 1;
    int sz = n2 - n1 + 1;
    double *Hr1 = new double[sz];
    double *Hr2 = new double[sz * sz];
    double *Hr3 = new double[sz * sz * sz];
    double *Hr4 = new double[sz * sz * sz * sz];
    double *Hr5 = new double[sz * sz * sz * sz * sz];
    apf_hplog_(&wx, &nw, Hr1, Hr2, Hr3, Hr4, Hr5, &n1, &n2);

    // weight 2
    const double Hm1m1 = Hr2[0];
    const double H1m1  = Hr2[2];
    const double Hm11  = Hr2[6];
    const double H11   = Hr2[8];

    // weight 3
    const double Hm1m1m1 = Hr3[0];
    const double H1m1m1  = Hr3[2];
    const double Hm11m1  = Hr3[6];
    const double H11m1   = Hr3[8];
    const double Hm1m11  = Hr3[18];
    const double H1m11   = Hr3[20];
    const double Hm111   = Hr3[24];
    const double H111    = Hr3[26];

    // Deallocate pointers
    delete[] Hr1;
    delete[] Hr2;
    delete[] Hr3;
    delete[] Hr4;
    delete[] Hr5;

    const double Hmp  = H11 + H1m1 - Hm11 - Hm1m1;
    const double Hmpm = H111 - H11m1 + H1m11 - H1m1m1 - Hm111 + Hm11m1 - Hm1m11 + Hm1m1m1;
    const double I    = 4 * z * Hmp;
    const double J    = 4 * z * L;
    const double K    = 4 * z * Hmpm;
    const double Lnxi = log(1 + xi / 4);
    //const double ln2  = log(2);
    const double ll   = a11 * a11 * ( - 32 / 3. * K / xi * ( 1 + 3 / xi ) - 128 / 27. * ( 17 + 120 / xi )
                                      + 16 / 27. * J * (3 + 136 / xi + 480 / xi / xi )
                                      + I * ( 32 / 3. * Lnxi / xi * ( 1 + 3 / xi ) + 16 / 9. * ( - 3 - 4 / xi + 24 / xi / xi ) )
                                      + ( 64 / 3. * I / xi * ( 1 + 3 / xi ) - 64 / 9. * ( - 1 + 12 / xi )
                                          + 16 / 9. * J * ( - 3 - 4 / xi + 24 / xi / xi ) ) * Lmu
                                      + ( 32 / 3. * J / xi * ( 1 + 3 / xi ) - 32 / 3. * ( 1 + 6 / xi ) ) * Lmu2) * log(x) / x / ( 1 + 4 / xi );
    const double nllc = ( a21 * ( - 64 / 3. * I / xi * ( 1 + 3 / xi ) + 64 / 9. * ( - 1 + 12 / xi )
                                  - 16 / 9. * J * ( - 3 - 4 / xi + 24 / xi / xi ) )
                          + a10 * a11 * ( 64 / 3. * K / xi * ( 1 + 3 / xi ) + 256 / 27. * ( 17 + 120 / xi )
                                          - 32 / 27. * J * (3 + 136 / xi + 480 / xi / xi )
                                          + I * ( - 64 / 3. * Lnxi / xi * ( 1 + 3 / xi ) - 32 / 9. * ( - 3 - 4 / xi + 24 / xi / xi ) ) )
                          + a11 * bt0 * ( - 32 / 3. * K / xi * ( 1 + 3 / xi ) - 128 / 27. * ( 17 + 120 / xi )
                                          + 16 / 27. * J * (3 + 136 / xi + 480 / xi / xi )
                                          + I * ( 32 / 3. * Lnxi / xi * ( 1 + 3 / xi ) + 16 / 9. * ( - 3 - 4 / xi + 24 / xi / xi ) ) )
                          + ( a21 * ( - 64 / 3. * J / xi * ( 1 + 3 / xi ) + 64 / 3. * ( 1 + 6 / xi ) )
                              + a10 * a11 * ( - 128 / 3. * I / xi * ( 1 + 3 / xi ) + 128 / 9. * ( - 1 + 12 / xi )
                                              - 32 / 9. * J * ( - 3 - 4 / xi + 24 / xi / xi ) )
                              + a11 * bt0 * ( 64 / 3. * I / xi * ( 1 + 3 / xi ) - 64 / 9. * ( - 1 + 12 / xi )
                                              + 16 / 9. * J * ( - 3 - 4 / xi + 24 / xi / xi ) ) ) * Lmu
                          + ( a11 * bt0 * ( 32 / 3. * J / xi * ( 1 + 3 / xi ) - 32 / 3. * ( 1 + 6 / xi ) )
                              + a10 * a11 * ( - 64 / 3. * J / xi * ( 1 + 3 / xi ) + 64 / 3. * ( 1 + 6 / xi ) ) ) * Lmu2 ) / x / ( 1 + 4 / xi );
    /*
    // Variation of nllc
    const double nllv = ( a10 * a11 * ( 64 / 3. * K / xi * ( 1 + 3 / xi ) + 256 / 27. * ( 17 + 120 / xi )
                                        - 32 / 27. * J * (3 + 136 / xi + 480 / xi / xi )
                                        + I * ( - 64 / 3. * Lnxi / xi * ( 1 + 3 / xi ) - 32 / 9. * ( - 3 - 4 / xi + 24 / xi / xi ) ) )
                          + ( a11 * bt0 * ( 32 / 3. * J / xi * ( 1 + 3 / xi ) - 32 / 3. * ( 1 + 6 / xi ) )
                              + a10 * a11 * ( - 64 / 3. * J / xi * ( 1 + 3 / xi ) + 64 / 3. * ( 1 + 6 / xi ) ) ) * Lmu2
                          + a11 * bt0 * ( - 32 / 3. * K / xi * ( 1 + 3 / xi ) - 128 / 27. * ( 17 + 120 / xi )
                                          - 256 / 9. * ( - 1 + 12 / xi ) * ln2 + 56 / 3. * ( - 1 + 12 / xi ) * zeta3
                                          + I * ( 32 / 3. * Lnxi / xi * ( 1 + 3 / xi ) + 16 / 9. * ( - 3 - 4 / xi + 24 / xi / xi )
                                                  + 256 / 3. / xi * ( 1 + 3 / xi ) * ln2 - 56 / xi * ( 1 + 3 / xi ) * zeta3 )
                                          + J * ( 16 / 27. * (3 + 136 / xi + 480 / xi / xi ) + 64 / 9. * ( - 3 - 4 / xi + 24 / xi / xi ) * ln2
                                                  - 14 / 3. * ( - 3 - 4 / xi + 24 / xi / xi ) * zeta3) )
                          + Lmu * ( a10 * a11 * ( - 128 / 3. * I / xi * ( 1 + 3 / xi ) + 128 / 9. * ( - 1 + 12 / xi )
                                                  - 32 / 9. * J * ( - 3 - 4 / xi + 24 / xi / xi ) )
                                    + a11 * bt0 * ( 64 / 3. * I / xi * ( 1 + 3 / xi ) - 64 / 9. * ( - 1 + 12 / xi )
                                                    - 256 / 3. * ( 1 + 6 / xi ) * ln2 + 56 * ( 1 + 6 / xi ) * zeta3
                                                    + J * ( 16 / 9. * ( - 3 - 4 / xi + 24 / xi / xi ) + 256 / 3. / xi * ( 1 + 3 / xi ) * ln2
                                                            - 56 / xi * ( 1 + 3 / xi ) * zeta3) ) ) ) / x / ( 1 + 4 / xi );
    */
    return ll + nllc;
  }

  //_________________________________________________________________________________
  CmsxL3psNC::CmsxL3psNC(int const& nf, double const& eta):
    Expression(),
    _cL3g(CmsxL3gNC{nf, eta})
  {
  }
  double CmsxL3psNC::Regular(double const& x) const
  {
    return CF / CA * _cL3g.Regular(x);
  }

  //_________________________________________________________________________________
  Cm0sx23gNC::Cm0sx23gNC(int const& nf, double const& eta):
    Expression(),
    _nf(nf),
    _eta(eta)
  {
  }
  double Cm0sx23gNC::Regular(double const& x) const
  {
    if (x >= 1)
      return 0;
    const double eta  = this->_eta;
    const double xi   = 4 * eta / ( 1 - eta );
    const double Lmu  = - log(xi);
    const double Lmu2 = Lmu * Lmu;
    const double LQ   = Lmu;
    const double LQ2  = LQ * LQ;
    const double LQ3  = LQ * LQ2;
    const double a11  = CA;
    const double a21  = _nf * ( 26 * CF - 23 * CA ) / 36;
    const double a10  = - ( 11 * CA + 2 * _nf * ( 1 - 2 * CF / CA ) ) / 12;
    const double bt0  = beta0qcd(_nf);
    const double ll   = CA * CA * log(x) * ( - 32 / 27. * ( - 71 + 18 * zeta2 ) * LQ - 208 / 9. * LQ2 + 32 / 9. * LQ3
                                             + Lmu2 * ( - 16 / 3. + 32 / 3. * LQ )
                                             + Lmu * ( 32 / 9. * ( - 5 + 6 * zeta2 ) + 416 / 9. * LQ - 32 / 3. * LQ2 )
                                             + 16 / 27. * ( - 92 + 78 * zeta2 - 72 * zeta3 ) ) / x;
    const double nllc = ( - 32 / 9. * a21 * ( - 5 + Pi2 )
                          + ( - 416 * a21 / 9 + 64 / 27. * a10 * a11 * ( - 71 + 3 * Pi2 ) - 32 / 27. * a11 * bt0 * ( - 71 + 3 * Pi2 ) ) * LQ
                          + ( 416 * a10 * a11 / 9 + 32 * a21 / 3 - 208 * a11 * bt0 / 9 ) * LQ2
                          + ( - 64 * a10 * a11 / 9 + 32 * a11 * bt0 / 9 ) * LQ3
                          + Lmu2 * ( 32 * a10 * a11 / 3 - 16 * a11 * bt0 / 3 + ( - 64 * a10 * a11 / 3 + 32 * a11 * bt0 / 3 ) * LQ )
                          + Lmu * ( 32 * a21 / 3 - 64 / 9. * a10 * a11 * ( - 5 + Pi2 ) + 32 / 9. * a11 * bt0 * ( - 5 + Pi2 )
                                    + ( - 832 * a10 * a11 / 9 - 64 * a21 / 3 + 416 * a11 * bt0 / 9 ) * LQ
                                    + ( 64 * a10 * a11 / 3 - 32 * a11 * bt0 / 3 ) * LQ2 )
                          - 32 / 27. * a10 * a11 * ( - 92 + 13 * Pi2 - 72 * zeta3 )
                          + 16 / 27. * a11 * bt0 * ( - 92 + 13 * Pi2 - 72 * zeta3 ) ) / x;
    /*
    // Variation of nllc
    const double nllv = ( ( - 64 * a10 * a11 / 9 + 32 * a11 * bt0 / 9 ) * LQ3
                          + Lmu2 * ( 32 * a10 * a11 / 3 - 16 * a11 * bt0 / 3 + ( - 64 * a10 * a11 / 3 + 32 * a11 * bt0 / 3 ) * LQ )
                          + LQ2 * ( 416 * a10 * a11 / 9 - 4 / 9. * a11 * bt0 * ( 52 + 96 * log(2) - 63 * zeta3 ) )
                          - 32 / 27. * a10 * a11 * ( - 92 + 13 * Pi2 - 72 * zeta3 )
                          + 4 / 27. * a11 * bt0 * ( - 368 + 52 * Pi2 - 480 * log(2) + 96 * Pi2 * log(2) + 27 * zeta3 - 63 * Pi2 * zeta3 )
                          + Lmu * ( - 64 / 9. * a10 * a11 * ( - 5 + Pi2 ) + ( 64 * a10 * a11 / 3 - 32 * a11 * bt0 / 3 ) * LQ2
                                    + LQ * ( - 832 * a10 * a11 / 9 + 8 / 9. * a11 * bt0 * ( 52 + 96 * log(2) - 63 * zeta3 ) )
                                    + 4 / 9. * a11 * bt0 * ( - 40 + 8 * Pi2 - 96 * log(2) + 63 * zeta3 ) )
                          + LQ * ( 64 / 27. * a10 * a11 * ( - 71 + 3 * Pi2 )
                                   - 4. / 27 * a11 * bt0 * ( - 568 + 24 * Pi2 - 1248 * log(2) + 819 * zeta3 ) ) ) / x;
    */
    return ll + nllc;
  }

  //_________________________________________________________________________________
  Cm0sx23psNC::Cm0sx23psNC(int const& nf, double const& eta):
    Expression(),
    _c23g(Cm0sx23gNC{nf, eta})
  {
  }
  double Cm0sx23psNC::Regular(double const& x) const
  {
    return CF / CA * _c23g.Regular(x);
  }

  //_________________________________________________________________________________
  Cm0sxL3gNC::Cm0sxL3gNC(int const& nf, double const& eta):
    Expression(),
    _nf(nf),
    _eta(eta)
  {
  }
  double Cm0sxL3gNC::Regular(double const& x) const
  {
    if (x >= 1)
      return 0;
    const double eta  = this->_eta;
    const double xi   = 4 * eta / ( 1 - eta );
    const double Lmu  = - log(xi);
    const double Lmu2 = Lmu * Lmu;
    const double LQ   = Lmu;
    const double LQ2  = LQ * LQ;
    const double a11  = CA;
    const double a21  = _nf * ( 26 * CF - 23 * CA ) / 36;
    const double a10  = - ( 11 * CA + 2 * _nf * ( 1 - 2 * CF / CA ) ) / 12;
    const double bt0  = beta0qcd(_nf);
    const double ll   = CA * CA * ( 32 / 27. * ( - 68 + 18 * zeta2 ) - 32 / 3. * Lmu2 - 64 / 9. * LQ
                                    - 32 / 3. * LQ2 + Lmu * ( 64 / 9. + 64 / 3. * LQ ) ) * log(x) / x;
    const double nllc = ( - 64 * a21 / 9 - 64 / 27. * a10 * a11 * ( - 68 + 3 * Pi2 )
                          + 32 / 27. * a11 * bt0 * ( - 68 + 3 * Pi2 )
                          + ( 64 * a10 * a11 / 3 - 32 * a11 * bt0 / 3 ) * Lmu2
                          + ( 128 * a10 * a11 / 9 - 64 * a21 / 3 - 64 * a11 * bt0 / 9 ) * LQ
                          + ( 64 * a10 * a11 / 3 - 32 * a11 * bt0 / 3 ) * LQ2
                          + ( - 128 * a10 * a11 / 9 + 64 * a21 / 3 + 64 * a11 * bt0 / 9
                              + ( - 128 * a10 * a11 / 3 + 64 * a11 * bt0 / 3 ) * LQ ) * Lmu ) / x;
    /*
    // Variation of nllc
    const double nllv = ( - 64 / 27 * a10 * a11 * ( - 68 + 3 * Pi2 )
                          + ( 64 * a10 * a11 / 3 - 32 * a11 * bt0 / 3 ) * Lmu2
                          + ( 64 * a10 * a11 / 3 - 32 * a11 * bt0 / 3 ) * LQ2
                          + ( - 128 * a10 * a11 / 9 + ( - 128 * a10 * a11 / 3 + 64 * a11 * bt0 / 3 ) * LQ
                              - 8 / 9 * a11 * bt0 * ( - 8 + 96 * log(2) - 63 * zeta3 ) ) * Lmu
                          + ( 128 * a10 * a11 / 9 + 8. / 9 * a11 * bt0 * ( - 8 + 96 * log(2) - 63 * zeta3 ) ) * LQ
                          + 8 / 27. * a11 * bt0 * ( - 272 + 12 * Pi2 + 96 * log(2) - 63 * zeta3 ) ) / x;
    */
    return ll + nllc;
  }

  //_________________________________________________________________________________
  Cm0sxL3psNC::Cm0sxL3psNC(int const& nf, double const& eta):
    Expression(),
    _cL3g(Cm0sxL3gNC{nf, eta})
  {
  }
  double Cm0sxL3psNC::Regular(double const& x) const
  {
    return CF / CA * _cL3g.Regular(x);
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
