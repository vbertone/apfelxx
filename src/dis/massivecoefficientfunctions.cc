//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/massivecoefficientfunctions.h"
#include "apfel/constants.h"
#include "apfel/specialfunctions.h"
#include "apfel/integrator.h"

namespace apfel
{
  //_________________________________________________________________________________
  Cm21gNC::Cm21gNC(double const& eta):
    Expression(eta)
  {
  }
  double Cm21gNC::Regular(double const& x) const
  {
    /*
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi  = 4 * eta / ( 1 - eta );
    double wr  = xi * ( 1 / z - 1 ) / 4 - 1;
    const double cm21g = xi * c2log_(&wr,&xi) / z / M_PI;
    return cm21g;
    */
    const double eta   = this->_eta;
    const double z     = eta * x;
    const double z2    = z * z;
    const double epsi  = ( 1 - eta ) / eta / 4;
    const double epsi2 = epsi * epsi;
    const double v     = sqrt( 1 - 4 * epsi * z / ( 1 - z ));
    const double cm21g = 4 * TR * ( ( z2 + ( 1 - z ) * ( 1 - z ) + 4 * epsi * z * ( 1 - 3 * z )
                                      - 8 * epsi2 * z2 ) * log( ( 1 + v ) / ( 1 - v ) )
                                    + ( 8 * z * ( 1 - z ) - 1 - 4 * epsi * z * ( 1 - z ) ) * v );
    return cm21g;
  }

  //_________________________________________________________________________________
  CmL1gNC::CmL1gNC(double const& eta):
    Expression(eta)
  {
  }
  double CmL1gNC::Regular(double const& x) const
  {
    /*
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi  = 4 * eta / ( 1 - eta );
    double wr  = xi * ( 1 / z - 1 ) / 4 - 1;
    const double cml1g = xi * cllog_(&wr,&xi) / z / M_PI;
    return cml1g;
    */
    const double eta   = this->_eta;
    const double z     = eta * x;
    const double z2    = z * z;
    const double epsi  = ( 1 - eta ) / eta / 4;
    const double v     = sqrt( 1 - 4 * z / ( 1 - z ) * epsi );
    const double cml1g = 4 * TR * ( - 8 * epsi * z2 * log( ( 1 + v ) / ( 1 - v ) )
                                    + 4 * v * z * ( 1 - z ) );
    return cml1g;
  }

  //_________________________________________________________________________________
  Cm22nsNC::Cm22nsNC(double const& eta):
    Expression(eta)
  {
    // Compute integral needed to enforce the Adler sum rule. See
    // eqs. (71) and (97) of https://arxiv.org/pdf/1001.2312.pdf.
    const Integrator Integrand{[&] (double const& y)->double{ return Regular(y); }};
    _adler = - Integrand.integrate(0,1,eps5);
  }
  double Cm22nsNC::Regular(double const& x) const
  {
    /*
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi  = 4 * eta / ( 1 - eta );
    double wr  = xi * ( 1 / z - 1 ) / 4 - 1;
    const double cm22ns = 16 * M_PI * xi * d2nloq_(&wr,&xi) / z;
    return cm22ns;
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
    const double cm22ns = CF * TR * ( ( 4. / 3 * ( 1 + z2 ) / omz - 16 / omz * zr2 * ( 1 - 9 * z + 9 * z2 ) )
                                      * ( log( omz / z2 ) * l1 + l1 * l2 + 2 * ( - dil1 + dil2 + dil3 - dil4 ) )
                                      + ( - 8. / 3 + 4 / omz + qr2 * ( 128. - 432 * z + 288 * z2 - 8 / omz ) ) * l1
                                      + ( 88. / 9 + 136. / 9 * z - 152. / 9 / omz
                                          + qr * ( 464. / 9 - 512. / 3 * z + 2048. / 9 * z2 )
                                          + qr2 * ( - 832. / 9 + 6208. / 9 * z - 11392. / 9 * z2 + 6016. / 9 * z3 ) ) * l3 / sq2
                                      + ( - 272. / 27 - 1244. / 27 * z + 718. / 27 / omz
                                          + qr * ( - 3424. / 27 + 15608. / 27 * z - 4304. / 9 * z2 + 20. / 27 / omz ) ) * sq1 );
    return cm22ns;
  }
  double Cm22nsNC::Local(double const&) const
  {
    return _adler;
  }

  //_________________________________________________________________________________
  CmL2nsNC::CmL2nsNC(double const& eta):
    Expression(eta)
  {
  }
  double CmL2nsNC::Regular(double const& x) const
  {
    /*
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi  = 4 * eta / ( 1 - eta );
    double wr  = xi * ( 1 / z - 1 ) / 4 - 1;
    const double cml2ns = 16 * M_PI * xi * dlnloq_(&wr,&xi) / z;
    return cml2ns;
    */
    const double eta = this->_eta;
    const double xi  = 4 * eta / ( 1 - eta );
    const double z   = eta * x;
    const double z2  = z * z;
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
    const double cml2ns = CF * TR * ( 96 * z * zr2 * ( log( omz / z2 ) * l1 + l1 * l2 + 2 * ( - dil1 + dil2 + dil3 - dil4 ) )
                                      + qr2 * ( 64 - 288 * z + 192 * z2 ) * l1
                                      + z * ( 16. / 3 - 416 * z / 3 / xi + 1408 * z2 / 3 / xi / xi ) * l3 / sq2
                                      + ( 16. / 3 - 400 * z / 18
                                          + z * ( - 160. / 3 + 3824 * z / 9 - 992 * z2 / 3 ) / omz / xi ) * sq1);
    return cml2ns;
  }

  //_________________________________________________________________________________
  Cm22gNC::Cm22gNC(double const& eta):
    Expression(eta)
  {
  }
  double Cm22gNC::Regular(double const& x) const
  {
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi  = 4 * eta / ( 1 - eta );
    double wr  = xi * ( 1 / z - 1 ) / 4 - 1;
    const double cm22g = 16 * M_PI * xi * c2nlog_(&wr,&xi) / z;
    return cm22g;
  }

  //_________________________________________________________________________________
  CmL2gNC::CmL2gNC(double const& eta):
    Expression(eta)
  {
  }
  double CmL2gNC::Regular(double const& x) const
  {
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi  = 4 * eta / ( 1 - eta );
    double wr  = xi * ( 1 / z - 1 ) / 4 - 1;
    const double cml2g = 16 * M_PI * xi * clnlog_(&wr,&xi) / z;
    return cml2g;
  }

  //_________________________________________________________________________________
  Cm22psNC::Cm22psNC(double const& eta):
    Expression(eta)
  {
  }
  double Cm22psNC::Regular(double const& x) const
  {
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi  = 4 * eta / ( 1 - eta );
    double wr  = xi * ( 1 / z - 1 ) / 4 - 1;
    const double cm22ps = 16 * M_PI * xi * c2nloq_(&wr,&xi) / z;
    return cm22ps;
  }

  //_________________________________________________________________________________
  CmL2psNC::CmL2psNC(double const& eta):
    Expression(eta)
  {
  }
  double CmL2psNC::Regular(double const& x) const
  {
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi  = 4 * eta / ( 1 - eta );
    double wr  = xi * ( 1 / z - 1 ) / 4 - 1;
    const double cml2ps = 16 * M_PI * xi * clnloq_(&wr,&xi) / z;
    return cml2ps;
  }

  //_________________________________________________________________________________
  Cm22bargNC::Cm22bargNC(double const& eta):
    Expression(eta)
  {
  }
  double Cm22bargNC::Regular(double const& x) const
  {
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi  = 4 * eta / ( 1 - eta );
    double wr  = xi * ( 1 / z - 1 ) / 4 - 1;
    const double cm22barg = 16 * M_PI * xi * c2nlobarg_(&wr,&xi) / z;
    return cm22barg;
  }

  //_________________________________________________________________________________
  CmL2bargNC::CmL2bargNC(double const& eta):
    Expression(eta)
  {
  }
  double CmL2bargNC::Regular(double const& x) const
  {
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi  = 4 * eta / ( 1 - eta );
    double wr  = xi * ( 1 / z - 1 ) / 4 - 1;
    const double cml2barg = 16 * M_PI * xi * clnlobarg_(&wr,&xi) / z;
    return cml2barg;
  }

  //_________________________________________________________________________________
  Cm22barpsNC::Cm22barpsNC(double const& eta):
    Expression(eta)
  {
  }
  double Cm22barpsNC::Regular(double const& x) const
  {
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi  = 4 * eta / ( 1 - eta );
    double wr  = xi * ( 1 / z - 1 ) / 4 - 1;
    const double cm22barps = 16 * M_PI * xi * c2nlobarq_(&wr,&xi) / z;
    return cm22barps;
  }

  //_________________________________________________________________________________
  CmL2barpsNC::CmL2barpsNC(double const& eta):
    Expression(eta)
  {
  }
  double CmL2barpsNC::Regular(double const& x) const
  {
    const double eta = this->_eta;
    const double z   = eta * x;
    double xi  = 4 * eta / ( 1 - eta );
    double wr  = xi * ( 1 / z - 1 ) / 4 - 1;
    const double cml2barps = 16 * M_PI * xi * clnlobarq_(&wr,&xi) / z;
    return cml2barps;
  }
}
