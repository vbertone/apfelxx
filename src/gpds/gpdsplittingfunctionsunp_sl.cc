//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/gpdsplittingfunctionsunp_sl.h"
#include "apfel/constants.h"
#include "apfel/specialfunctions.h"

#include <cmath>

namespace apfel
{
  //_________________________________________________________________________________
  Pgpd0ns::Pgpd0ns(double const& xi):
    Expression(1 / xi)
  {
  }
  double Pgpd0ns::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? - 2 * CF * ( 1 + y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * ( 1 + ( 1 + kappa ) * y + ( 1 + kappa - pow(kappa, 2) ) * pow(y, 2) ) / ( 1 + y ) / ( 1 - ky2 ) : 0);
  }
  double Pgpd0ns::Singular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( (y <= 1 ? 2 : 0) - (kappa > 1 ? 1 : 0) ) / ( 1 - y );
  }
  double Pgpd0ns::Local(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( 2 * log(1 - y) + 3. / 2. - log(std::abs(1 - pow(kappa, 2))) );
  }
  double Pgpd0ns::LocalPP(double const& c, double const& d) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( - (kappa > 1 ? log(1 - c) - (d > 1 ? log(d - 1) : 0) : 0) );
  }

  //_________________________________________________________________________________
  Pgpd0qq::Pgpd0qq(double const& xi):
    Expression(1 / xi)
  {
  }
  double Pgpd0qq::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? - 2 * CF * ( 1 + y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * ( 1 + y + kappa * y + pow(kappa, 3) * pow(y, 2) ) / kappa / ( 1 + y ) / ( 1 - ky2 ) : 0);
  }
  double Pgpd0qq::Singular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( (y <= 1 ? 2 : 0) - (kappa > 1 ? 1 : 0) ) / ( 1 - y );
  }
  double Pgpd0qq::Local(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( 2 * log(1 - y) + 3. / 2. - log(std::abs(1 - pow(kappa, 2))) );
  }
  double Pgpd0qq::LocalPP(double const& c, double const& d) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( - (kappa > 1 ? log(1 - c) - (d > 1 ? log(d - 1) : 0) : 0) );
  }

  //_________________________________________________________________________________
  Pgpd0qg::Pgpd0qg(int const& nf, double const& xi):
    Expression(1 / xi),
    _nf(nf)
  {
  }
  double Pgpd0qg::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 4 * _nf * TR * ( pow(y, 2) + pow(1 - y, 2) - ky2 ) / pow(1 - ky2, 2) : 0)
           + (kappa > 1 ? 4 * _nf * TR * ( 1 - kappa ) * ( 1 - kappa * ( kappa + 2 ) * pow(y, 2) ) / kappa / pow(1 - ky2, 2) : 0);
  }

  //_________________________________________________________________________________
  Pgpd0gq::Pgpd0gq(double const& xi):
    Expression(1 / xi)
  {
  }
  double Pgpd0gq::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 2 * CF * ( 1 + pow(1 - y, 2) - ky2 ) / y / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? - 2 * CF * pow(1 - kappa, 2) / kappa / ( 1 - ky2 ) : 0);
  }

  //_________________________________________________________________________________
  Pgpd0gg::Pgpd0gg(int const& nf, double const& xi):
    Expression(1 / xi),
    _nf(nf)
  {
  }
  double Pgpd0gg::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    const double k2    = kappa * kappa;
    return 4 * CA * (y <= 1 ? - ( 1 + k2 * y ) / ( 1 - ky2 ) + ( ( 1 - y ) / y + y * ( 1 - y ) ) / pow(1 - ky2, 2) : 0)
           + 2 * CA * (kappa > 1 ? 2 * ( 1 - kappa ) * ( 1 + pow(y, 2) ) / pow(1 - ky2, 2) + k2 * ( 1 + y ) / ( 1 - ky2 )
                       + ( 1 - k2 ) * ( 2 - 1 / kappa - 1 / ( 1 + y ) ) / ( 1 - ky2 ) : 0);
  }
  double Pgpd0gg::Singular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CA * ( (y <= 1 ? 2 : 0) - (kappa > 1 ? 1 : 0) ) / ( 1 - y );
  }
  double Pgpd0gg::Local(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 4 * CA * log(1 - y) + ( 11 * CA - 4 * _nf * TR ) / 3 - 2 * CA * log(std::abs(1 - pow(kappa, 2)));
  }
  double Pgpd0gg::LocalPP(double const& c, double const& d) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CA * ( - (kappa > 1 ? log(1 - c) - (d > 1 ? log(d - 1) : 0) : 0) );
  }

  //_________________________________________________________________________________
  // NLO non-singlet kernels. Agpd1ns: coefficient A of the 1/(1-y)
  // singularity (cusp, kappa-independent). Cgpd1ns: constant c0 of
  // the delta(1-y) term, c(kappa) = c0 - (A/2) ln|1 - kappa^2|, fixed
  // by the vector sum rule.
  namespace
  {
    double Agpd1ns(int const& nf)
    {
      return CF * CA * 4 * ( 67 - 3 * Pi2 ) / 9 - CF * TR * nf * 80. / 9;
    }
    double Cgpd1ns(int const& nf)
    {
      return CF * CF * ( 3. / 2 - 2 * Pi2 + 24 * zeta3 )
             + CF * CA * ( 17. / 6 + 22 * Pi2 / 9 - 12 * zeta3 )
             - CF * TR * nf * ( 2. / 3 + 8 * Pi2 / 9 );
    }
  }

  namespace
  {
    // Regular part of the two-loop (O(as^2), as = alpha_s/(4 pi))
    // quark-in-quark GPD kernel P_qq^[1](x, xi), in parton-in-parton
    // variables (parent momentum fraction 1, skewness xi > 0, x !=
    // 1), summed over colour factors.  Regions: xi < 1: DGLAP xi<x<1,
    // ERBL -xi<x<xi, q->qbar -1<x<-xi; xi > 1: -1<x<1, 1<x<xi,
    // -xi<x<-1.  The singular part A/(1-x) (xi < 1) or (A/2)/|1-x|
    // (xi > 1) is included here and subtracted in Regular().
    double Pgpd1qq(double const& x, double const& xi, int const& nf)
    {
      const double x2  = x * x;
      const double x3  = x2 * x;
      const double x4  = x2 * x2;
      const double xi2 = xi * xi;
      const double xi3 = xi2 * xi;
      const double ln2 = log(2.);
      if (xi < 1)
        {
          if (x >= 1 || x <= -1)
            return 0;
          // small xi (and xi << |x|): series in xi, avoids the 1/xi cancellations of the closed forms
          const bool ser = ( xi < 1e-4 && xi < 1e-3 * std::abs(x) );
          if (x > xi && ser)
            {
              const double s1 = log(1 - x);
              const double s2 = log(x);
              return CF * CF * ( (20*((-1 + x) * (-1 + x)) - 2*(s2 * s2)*(-1 + x2) + 4*s2*(3 + 2*x - 2*x2 + 2*s1*(1 + x2)))/(-1 + x) + (2*(-10 + 64*x - 117*x2 - 9*(s2 * s2)*(-1 + x2)*x2 + 18*s1*(-1 + x2)*(1 + x2 + 2*s2*x2) - 44*x3 + s2*(-3 + 12*x + 72*x2 - 45*x4) + 107*x4)*xi2)/(9.*(-1 + x)*x2) )
                     + CF * CA * ( -0.1111111111111111*(374 - 480*x + 374*x2 - 6*Pi2*(1 + x2) + 18*(s2 * s2)*(1 + x2) + 6*s2*(17 + 5*x2))/(-1 + x) + ((-15 + 105*x + (-65 + 6*Pi2)*x2 - 18*(s2 * s2)*(1 + x)*x2 + (-425 + 6*Pi2)*x3 - 6*s2*(1 - 3*x + 2*x2 + 6*x3))*xi2)/(9.*x2) )
                     + CF * TR * nf * ( (88 - 96*x + 88*x2 + 24*s2*(1 + x2))/(9.*(-1 + x)) + (4*(1 - 3*x + 7*x2 + 6*s2*(1 + x)*x2 + 27*x3)*xi2)/(9.*x2) );
            }
          else if (x > xi)
            {
              const double l1 = log(1 - x);
              const double l2 = log(1 - xi);
              const double l3 = log(x - xi);
              const double l4 = log(1 + xi);
              const double l5 = log(x + xi);
              const double d1 = dilog((2*xi)/(1 + xi));
              const double d2 = dilog((2*xi)/(x + xi));
              return CF * CF * ( (4*(-1 + x))/((-1 + xi)*(1 + xi)) - (2*l4*l5*(x - xi))/(xi*(1 + xi)) - (2*l2*l3*(x + xi))/((-1 + xi)*xi) - (4*l2*l4*(-1 + x - 2*xi)*(x + xi))/((-1 + x)*xi*(1 + xi)) + (4*l2*l5*(-1 + x - 2*xi)*(x + xi))/((-1 + x)*xi*(1 + xi)) + (4*l3*l4*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) - (4*l3*l5*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) - ((l2 * l2)*(-x + x2 - 3*xi - x*xi - 4*xi2))/((-1 + x)*xi*(1 + xi)) + ((l5 * l5)*(-3*x + 3*x2 + xi + 3*x*xi - 4*xi2))/((-1 + x)*(-1 + xi)*xi) + (2*l4*(x + xi)*(2 - 2*x + 2*xi + x*xi - 3*xi2))/((-1 + x)*(-1 + xi)*xi*(1 + xi)) - (2*l5*(x + xi)*(2 - 2*x + 2*xi + x*xi - 3*xi2))/((-1 + x)*(-1 + xi)*xi*(1 + xi)) + (4*l1*l2*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (4*l1*l3*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) + (4*l1*l4*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (4*l1*l5*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (8*d1*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (8*d2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (2*l2*(x - xi)*(-2 + 2*x + 2*xi + x*xi + 3*xi2))/((-1 + x)*(-1 + xi)*xi*(1 + xi)) - (2*l3*(x - xi)*(-2 + 2*x + 2*xi + x*xi + 3*xi2))/((-1 + x)*(-1 + xi)*xi*(1 + xi)) + ((l3 * l3)*(x - xi)*(-1 + x + xi + 3*x*xi + 4*xi2))/((-1 + x)*(-1 + xi)*xi*(1 + xi)) - ((l4 * l4)*(-5*x + 5*x2 + 3*xi + x2*xi - 5*xi2 + 5*x*xi2 - 4*xi3))/((-1 + x)*(-1 + xi)*xi*(1 + xi)) )
                     + CF * CA * ( (-11*l4*(-1 + x - 2*xi)*(x + xi))/(3.*(-1 + x)*xi*(1 + xi)) + (2*l2*l4*(-1 + x - 2*xi)*(x + xi))/((-1 + x)*xi*(1 + xi)) + (11*l5*(-1 + x - 2*xi)*(x + xi))/(3.*(-1 + x)*xi*(1 + xi)) - (2*l2*l5*(-1 + x - 2*xi)*(x + xi))/((-1 + x)*xi*(1 + xi)) - (11*l2*(x - xi)*(-1 + x + 2*xi))/(3.*(-1 + x)*(-1 + xi)*xi) + (11*l3*(x - xi)*(-1 + x + 2*xi))/(3.*(-1 + x)*(-1 + xi)*xi) - (2*l3*l4*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) + (2*l3*l5*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) + (4*(59 - 51*x + 59*x2 - 67*xi2))/(9.*(-1 + x)*(-1 + xi)*(1 + xi)) - (2*Pi2*(1 + x2 - 2*xi2))/(3.*(-1 + x)*(-1 + xi)*(1 + xi)) + (4*d1*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (4*d2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (2*(l4 * l4)*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (2*(l5 * l5)*(x + xi2))/((-1 + xi)*xi*(1 + xi)) )
                     + CF * TR * nf * ( (4*l4*(-1 + x - 2*xi)*(x + xi))/(3.*(-1 + x)*xi*(1 + xi)) - (4*l5*(-1 + x - 2*xi)*(x + xi))/(3.*(-1 + x)*xi*(1 + xi)) + (4*l2*(x - xi)*(-1 + x + 2*xi))/(3.*(-1 + x)*(-1 + xi)*xi) - (4*l3*(x - xi)*(-1 + x + 2*xi))/(3.*(-1 + x)*(-1 + xi)*xi) - (16*(4 - 3*x + 4*x2 - 5*xi2))/(9.*(-1 + x)*(-1 + xi)*(1 + xi)) );
            }
          else if (x > -xi)
            {
              const double l1 = log(1 - x);
              const double l2 = log(xi);
              const double l3 = log(1 + xi);
              const double l4 = log(-x + xi);
              const double l5 = log(x + xi);
              const double d1 = dilog((x + xi)/(2.*xi));
              const double d2 = dilog((x + xi)/(1 + xi));
              return CF * CF * ( (l4*(x - xi))/(xi*(1 + xi)) + ((l4 * l4)*(x - xi))/(xi*(1 + xi)) - (2*l2*l4*(x + xi))/((-1 + xi)*xi) - (2*l4*ln2*(x + xi))/((-1 + xi)*xi) + (2*(x + xi))/(xi*(1 + xi)) - (3*l5*(-1 + x - 2*xi)*(x + xi))/((-1 + x)*xi*(1 + xi)) + (2*(l5 * l5)*(-1 + x - 2*xi)*(x + xi))/((-1 + x)*xi*(1 + xi)) - (4*Pi2*(-1 + x - 2*xi)*(x + xi))/(3.*(-1 + x)*xi*(1 + xi)) + (4*l3*l4*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) + ((l3 * l3)*(x + xi)*(1 - x + xi + 3*x*xi - 4*xi2))/((-1 + x)*(-1 + xi)*xi*(1 + xi)) + (2*l3*(x + xi)*(2 - 2*x + 2*xi + x*xi - 3*xi2))/((-1 + x)*(-1 + xi)*xi*(1 + xi)) + (8*d2*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (4*l1*l3*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) + (4*l1*l5*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (4*l4*l5*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (8*d1*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (2*l2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (2*(l2 * l2)*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (2*ln2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (4*l2*ln2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (2*(ln2 * ln2)*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (2*l2*l5*(3*x + xi + x*xi + 3*xi2))/((-1 + xi)*xi*(1 + xi)) + (2*l5*ln2*(3*x + xi + x*xi + 3*xi2))/((-1 + xi)*xi*(1 + xi)) - (2*l3*l5*(x - xi)*(-1 + x + xi + 3*x*xi + 4*xi2))/((-1 + x)*(-1 + xi)*xi*(1 + xi)) )
                     + CF * CA * ( (2*(-59 + 59*x - 67*xi)*(x + xi))/(9.*(-1 + x)*xi*(1 + xi)) - (11*l3*(-1 + x - 2*xi)*(x + xi))/(3.*(-1 + x)*xi*(1 + xi)) + (11*l5*(-1 + x - 2*xi)*(x + xi))/(3.*(-1 + x)*xi*(1 + xi)) + (Pi2*(-1 + x - 2*xi)*(x + xi))/(3.*(-1 + x)*xi*(1 + xi)) - (2*l3*l4*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) + (2*l3*l5*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) + (2*l4*l5*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) - (4*d2*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) + (4*l1*l3*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (2*(l3 * l3)*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (4*l1*l5*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) + (4*d1*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (2*(l2 * l2)*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (4*l2*l5*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (4*l2*ln2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (4*l5*ln2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (2*(ln2 * ln2)*(x + xi2))/((-1 + xi)*xi*(1 + xi)) )
                     + CF * TR * nf * ( (-8*(-4 + 4*x - 5*xi)*(x + xi))/(9.*(-1 + x)*xi*(1 + xi)) + (4*l3*(-1 + x - 2*xi)*(x + xi))/(3.*(-1 + x)*xi*(1 + xi)) - (4*l5*(-1 + x - 2*xi)*(x + xi))/(3.*(-1 + x)*xi*(1 + xi)) );
            }
          else if (ser)
            {
              const double s1 = log(1 - x);
              const double s2 = log(-x);
              const double s3 = log(x/(-1 + x));
              const double s4 = dilog(1/(1 - x));
              const double s5 = dilog(x/(-1 + x));
              return CF * CF * ( (4*(4 + 2*s2 - 2*s1*s2 + (s2 * s2) - 4*s2*x - 4*x2 + 2*s2*x2 - 2*s1*s2*x2 + (s2 * s2)*x2 - 2*s4*(1 + x2) + 2*s5*(1 + x2)))/(-1 + x) - (4*(2 + 3*s2 - 38*x - 12*s2*x - 36*x2 + 9*s1*x2 + 27*s2*x2 - 18*s1*s2*x2 + 9*(s2 * s2)*x2 + 9*s3*x2 + 18*s4*(-1 + x2)*x2 - 18*s5*(-1 + x2)*x2 + 38*x3 + 24*s2*x3 + 34*x4 + 9*s1*x4 - 24*s2*x4 + 18*s1*s2*x4 - 9*(s2 * s2)*x4 + 9*s3*x4)*xi2)/(9.*(-1 + x)*x2) )
                     + CF * CA * ( (-2*(4 + 2*s2 - 2*s1*s2 + (s2 * s2) - 4*s2*x - 4*x2 + 2*s2*x2 - 2*s1*s2*x2 + (s2 * s2)*x2 - 2*s4*(1 + x2) + 2*s5*(1 + x2)))/(-1 + x) + (2*(2 + 3*s2 - 38*x - 12*s2*x - 36*x2 + 9*s1*x2 + 27*s2*x2 - 18*s1*s2*x2 + 9*(s2 * s2)*x2 + 9*s3*x2 + 18*s4*(-1 + x2)*x2 - 18*s5*(-1 + x2)*x2 + 38*x3 + 24*s2*x3 + 34*x4 + 9*s1*x4 - 24*s2*x4 + 18*s1*s2*x4 - 9*(s2 * s2)*x4 + 9*s3*x4)*xi2)/(9.*(-1 + x)*x2) )
                     + CF * TR * nf * ( 0 );
            }
          else
            {
              const double l1 = log(1 - x);
              const double l2 = log(1 - xi);
              const double l3 = log(-x - xi);
              const double l4 = log(1 + xi);
              const double l5 = log(-x + xi);
              const double d1 = dilog((-1 + xi)/(-1 + x));
              const double d2 = dilog((2*xi)/(1 + xi));
              const double d3 = dilog((2*xi)/(-x + xi));
              const double d4 = dilog((x + xi)/(-1 + x));
              return CF * CF * ( (4*l2*l4*(-1 + x - 2*xi)*(x + xi))/((-1 + x)*xi*(1 + xi)) - (4*l3*l5*(-1 + x - 2*xi)*(x + xi))/((-1 + x)*xi*(1 + xi)) - (4*l3*l4*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) + (4*l2*l5*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) + (8*d1*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (8*d4*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (8*l1*l2*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) + (8*l1*l3*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) + (8*d2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (8*d3*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (4*(l4 * l4)*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (4*(l5 * l5)*(x + xi2))/((-1 + xi)*xi*(1 + xi)) )
                     + CF * CA * ( (-2*l2*l4*(-1 + x - 2*xi)*(x + xi))/((-1 + x)*xi*(1 + xi)) + (2*l3*l5*(-1 + x - 2*xi)*(x + xi))/((-1 + x)*xi*(1 + xi)) + (2*l3*l4*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) - (2*l2*l5*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) - (4*d1*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) + (4*d4*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) + (4*l1*l2*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (4*l1*l3*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (4*d2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (4*d3*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (2*(l4 * l4)*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (2*(l5 * l5)*(x + xi2))/((-1 + xi)*xi*(1 + xi)) )
                     + CF * TR * nf * ( 0 );
            }
        }
      else
        {
          if (x >= xi || x <= -xi)
            return 0;
          if (x > 1)
            {
              const double l1 = log(-1 + x);
              const double l2 = log(-1 + xi);
              const double l3 = log(xi);
              const double l4 = log(1 + xi);
              const double l5 = log(-x + xi);
              const double l6 = log(x + xi);
              const double d1 = dilog((-1 + xi)/(2.*xi));
              const double d2 = dilog((1 + xi)/(x + xi));
              return CF * CF * ( (-2*(x - xi))/((-1 + xi)*xi) + (2*l3*l6*(x - xi))/(xi*(1 + xi)) + (2*l6*ln2*(x - xi))/(xi*(1 + xi)) + (2*l2*l5*(x + xi))/((-1 + xi)*xi) - (2*l3*l5*(x + xi))/((-1 + xi)*xi) - (l6*(x + xi))/((-1 + xi)*xi) - (2*l5*ln2*(x + xi))/((-1 + xi)*xi) - (4*l2*l6*(-1 + x - 2*xi)*(x + xi))/((-1 + x)*xi*(1 + xi)) + (4*l4*l6*(-1 + x - 2*xi)*(x + xi))/((-1 + x)*xi*(1 + xi)) + (4*l2*l4*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) + (3*l5*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) - (2*(l5 * l5)*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) + ((l2 * l2)*(-x + x2 - 3*xi - x*xi - 4*xi2))/((-1 + x)*xi*(1 + xi)) - (8*d2*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (4*l1*l2*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (8*l1*l4*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) + (4*l1*l5*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) + (8*l1*l6*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) + (4*Pi2*(1 + x2 - 2*xi2))/(3.*(-1 + x)*(-1 + xi)*(1 + xi)) + (8*d1*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (2*l3*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (8*l2*l3*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (6*(l3 * l3)*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (4*l5*l6*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (2*ln2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (8*l2*ln2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (12*l3*ln2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (6*(ln2 * ln2)*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (2*l2*(x - xi)*(-2 + 2*x + 2*xi + x*xi + 3*xi2))/((-1 + x)*(-1 + xi)*xi*(1 + xi)) - ((l6 * l6)*(-x + x2 + 3*xi + 5*x2*xi - xi2 + x*xi2 - 8*xi3))/((-1 + x)*(-1 + xi)*xi*(1 + xi)) )
                     + CF * CA * ( (2*l2*l6*(-1 + x - 2*xi)*(x + xi))/((-1 + x)*xi*(1 + xi)) - (2*l4*l6*(-1 + x - 2*xi)*(x + xi))/((-1 + x)*xi*(1 + xi)) - (Pi2*(-1 + x - 2*xi)*(x + xi))/(3.*(-1 + x)*xi*(1 + xi)) + (11*l2*(x - xi)*(-1 + x + 2*xi))/(3.*(-1 + x)*(-1 + xi)*xi) - (2*l2*l4*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) - (11*l5*(x - xi)*(-1 + x + 2*xi))/(3.*(-1 + x)*(-1 + xi)*xi) - (2*(x - xi)*(-59 + 59*x + 67*xi))/(9.*(-1 + x)*(-1 + xi)*xi) + (4*d2*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) + (4*l1*l4*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (4*l1*l6*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) + (2*(l6 * l6)*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (4*d1*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (4*l2*l3*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (2*(l3 * l3)*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (4*l2*ln2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (4*l3*ln2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (2*(ln2 * ln2)*(x + xi2))/((-1 + xi)*xi*(1 + xi)) )
                     + CF * TR * nf * ( (-4*l2*(x - xi)*(-1 + x + 2*xi))/(3.*(-1 + x)*(-1 + xi)*xi) + (4*l5*(x - xi)*(-1 + x + 2*xi))/(3.*(-1 + x)*(-1 + xi)*xi) + (8*(x - xi)*(-4 + 4*x + 5*xi))/(9.*(-1 + x)*(-1 + xi)*xi) );
            }
          else if (x > -1)
            {
              const double l1 = log(1 - x);
              const double l2 = log(xi);
              const double l3 = log(1 + xi);
              const double l4 = log(-x + xi);
              const double l5 = log(x + xi);
              const double d1 = dilog((-1 + x)/(x + xi));
              const double d2 = dilog((x + xi)/(2.*xi));
              return CF * CF * ( (l4*(x - xi))/(xi*(1 + xi)) + ((l4 * l4)*(x - xi))/(xi*(1 + xi)) - (2*l2*l4*(x + xi))/((-1 + xi)*xi) - (2*l4*ln2*(x + xi))/((-1 + xi)*xi) + (2*(x + xi))/(xi*(1 + xi)) - (3*l5*(-1 + x - 2*xi)*(x + xi))/((-1 + x)*xi*(1 + xi)) + (4*l3*l4*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) - ((l3 * l3)*(-x + x2 + 3*xi + x*xi - 4*xi2))/((-1 + x)*(-1 + xi)*xi) + (2*l3*(x + xi)*(2 - 2*x + 2*xi + x*xi - 3*xi2))/((-1 + x)*(-1 + xi)*xi*(1 + xi)) + (8*d1*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) + (4*l1*l3*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (4*l1*l5*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (4*l4*l5*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (8*d2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (2*l2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (2*(l2 * l2)*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (2*ln2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (4*l2*ln2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (2*(ln2 * ln2)*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (4*Pi2*(x + xi2))/(3.*(-1 + xi)*xi*(1 + xi)) + (2*l2*l5*(3*x + xi + x*xi + 3*xi2))/((-1 + xi)*xi*(1 + xi)) + (2*l5*ln2*(3*x + xi + x*xi + 3*xi2))/((-1 + xi)*xi*(1 + xi)) - (2*l3*l5*(x - xi)*(-1 + x + xi + 3*x*xi + 4*xi2))/((-1 + x)*(-1 + xi)*xi*(1 + xi)) + (2*(l5 * l5)*(x - x2 + 3*xi + 3*x2*xi + xi2 - x*xi2 - 6*xi3))/((-1 + x)*(-1 + xi)*xi*(1 + xi)) )
                     + CF * CA * ( (2*(-59 + 59*x - 67*xi)*(x + xi))/(9.*(-1 + x)*xi*(1 + xi)) - (11*l3*(-1 + x - 2*xi)*(x + xi))/(3.*(-1 + x)*xi*(1 + xi)) + (11*l5*(-1 + x - 2*xi)*(x + xi))/(3.*(-1 + x)*xi*(1 + xi)) - (2*l3*l4*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) + (2*l3*l5*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) + (2*l4*l5*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) - (Pi2*(x - xi)*(-1 + x + 2*xi))/(3.*(-1 + x)*(-1 + xi)*xi) - (4*d1*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (2*(l5 * l5)*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) + (4*d2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (2*(l2 * l2)*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (4*l2*l5*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (4*l2*ln2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (4*l5*ln2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (2*(ln2 * ln2)*(x + xi2))/((-1 + xi)*xi*(1 + xi)) )
                     + CF * TR * nf * ( (-8*(-4 + 4*x - 5*xi)*(x + xi))/(9.*(-1 + x)*xi*(1 + xi)) + (4*l3*(-1 + x - 2*xi)*(x + xi))/(3.*(-1 + x)*xi*(1 + xi)) - (4*l5*(-1 + x - 2*xi)*(x + xi))/(3.*(-1 + x)*xi*(1 + xi)) );
            }
          else
            {
              const double l1 = log(1 - x);
              const double l2 = log(-1 + xi);
              const double l3 = log(xi);
              const double l4 = log(1 + xi);
              const double l5 = log(-x + xi);
              const double l6 = log(x + xi);
              const double d1 = dilog((-1 + xi)/(2.*xi));
              const double d2 = dilog((-1 + xi)/(-x + xi));
              return CF * CF * ( (l5*(x - xi))/(xi*(1 + xi)) + (2*l3*l6*(x - xi))/(xi*(1 + xi)) - (2*l4*l6*(x - xi))/(xi*(1 + xi)) + (2*l6*ln2*(x - xi))/(xi*(1 + xi)) - (2*l3*l5*(x + xi))/((-1 + xi)*xi) - (2*l5*ln2*(x + xi))/((-1 + xi)*xi) + (2*(x + xi))/(xi*(1 + xi)) - (3*l6*(-1 + x - 2*xi)*(x + xi))/((-1 + x)*xi*(1 + xi)) + (2*(l6 * l6)*(-1 + x - 2*xi)*(x + xi))/((-1 + x)*xi*(1 + xi)) - (4*Pi2*(-1 + x - 2*xi)*(x + xi))/(3.*(-1 + x)*xi*(1 + xi)) - (4*l2*l4*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) - (4*l2*l5*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) + (4*l4*l5*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) - ((l4 * l4)*(-x + x2 + 3*xi + x*xi - 4*xi2))/((-1 + x)*(-1 + xi)*xi) + (2*l4*(x + xi)*(2 - 2*x + 2*xi + x*xi - 3*xi2))/((-1 + x)*(-1 + xi)*xi*(1 + xi)) + (8*d2*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) + (8*l1*l2*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) + (4*l1*l4*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (8*l1*l5*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (4*l1*l6*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (8*d1*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (2*l3*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (8*l2*l3*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (2*(l3 * l3)*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (4*l5*l6*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (2*ln2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (8*l2*ln2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (4*l3*ln2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (2*(ln2 * ln2)*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + ((l5 * l5)*(x - x2 + 3*xi + 5*x2*xi + xi2 - x*xi2 - 8*xi3))/((-1 + x)*(-1 + xi)*xi*(1 + xi)) )
                     + CF * CA * ( (2*(-59 + 59*x - 67*xi)*(x + xi))/(9.*(-1 + x)*xi*(1 + xi)) - (11*l4*(-1 + x - 2*xi)*(x + xi))/(3.*(-1 + x)*xi*(1 + xi)) + (11*l6*(-1 + x - 2*xi)*(x + xi))/(3.*(-1 + x)*xi*(1 + xi)) + (Pi2*(-1 + x - 2*xi)*(x + xi))/(3.*(-1 + x)*xi*(1 + xi)) + (2*l2*l4*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) + (2*l2*l5*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) - (2*l4*l5*(x - xi)*(-1 + x + 2*xi))/((-1 + x)*(-1 + xi)*xi) - (4*d2*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (4*l1*l2*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) + (4*l1*l5*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) - (2*(l5 * l5)*(1 + x2 - 2*xi2))/((-1 + x)*(-1 + xi)*(1 + xi)) + (4*d1*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (4*l2*l3*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (2*(l3 * l3)*(x + xi2))/((-1 + xi)*xi*(1 + xi)) - (4*l2*ln2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (4*l3*ln2*(x + xi2))/((-1 + xi)*xi*(1 + xi)) + (2*(ln2 * ln2)*(x + xi2))/((-1 + xi)*xi*(1 + xi)) )
                     + CF * TR * nf * ( (-8*(-4 + 4*x - 5*xi)*(x + xi))/(9.*(-1 + x)*xi*(1 + xi)) + (4*l4*(-1 + x - 2*xi)*(x + xi))/(3.*(-1 + x)*xi*(1 + xi)) - (4*l6*(-1 + x - 2*xi)*(x + xi))/(3.*(-1 + x)*xi*(1 + xi)) );
            }
        }
    }
  }

  //_________________________________________________________________________________
  Pgpd1nsp::Pgpd1nsp(int const& nf, double const& xi):
    Expression(1 / xi),
    _nf(nf)
  {
  }
  double Pgpd1nsp::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return Pgpd1qq(y, kappa * y, _nf) - Pgpd1qq(- y, kappa * y, _nf) - Singular(y);
  }
  double Pgpd1nsp::Singular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return Agpd1ns(_nf) * ( (y <= 1 ? 1 : 0) - (kappa > 1 ? 0.5 : 0) ) / ( 1 - y );
  }
  double Pgpd1nsp::Local(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double A     = Agpd1ns(_nf);
    return A * log(1 - y) + Cgpd1ns(_nf) - A / 2 * log(std::abs(1 - kappa * kappa));
  }
  double Pgpd1nsp::LocalPP(double const& c, double const& d) const
  {
    const double kappa = 1 / _eta / _extvar;
    return ( kappa > 1 ? - Agpd1ns(_nf) / 2 * ( log(1 - c) - (d > 1 ? log(d - 1) : 0) ) : 0 );
  }

  //_________________________________________________________________________________
  Pgpd1nsm::Pgpd1nsm(int const& nf, double const& xi):
    Expression(1 / xi),
    _nf(nf)
  {
  }
  double Pgpd1nsm::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return Pgpd1qq(y, kappa * y, _nf) + Pgpd1qq(- y, kappa * y, _nf) - Singular(y);
  }
  double Pgpd1nsm::Singular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return Agpd1ns(_nf) * ( (y <= 1 ? 1 : 0) - (kappa > 1 ? 0.5 : 0) ) / ( 1 - y );
  }
  double Pgpd1nsm::Local(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double A     = Agpd1ns(_nf);
    return A * log(1 - y) + Cgpd1ns(_nf) - A / 2 * log(std::abs(1 - kappa * kappa));
  }
  double Pgpd1nsm::LocalPP(double const& c, double const& d) const
  {
    const double kappa = 1 / _eta / _extvar;
    return ( kappa > 1 ? - Agpd1ns(_nf) / 2 * ( log(1 - c) - (d > 1 ? log(d - 1) : 0) ) : 0 );
  }
}
