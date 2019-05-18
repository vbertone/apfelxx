//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/zeromasscoefficientfunctions.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________________
  C21ns::C21ns():
    Expression()
  {
  }
  double C21ns::Regular(double const& x) const
  {
    return 2 * CF * ( - ( 1 + x ) * log( 1 - x ) - ( 1 + pow(x,2) ) * log(x) / ( 1 - x ) + 3 + 2 * x );
  }
  double C21ns::Singular(double const& x) const
  {
    return 2 * CF * ( 2 * log( 1 - x ) - 3 / 2. ) / ( 1 - x );
  }
  double C21ns::Local(double const& x) const
  {
    return 2 * CF * ( pow(log(1-x),2) - 3 * log( 1 - x ) / 2 - ( 2 * zeta2 + 9 / 2. ) );
  }

  //_________________________________________________________________________________
  C21g::C21g():
    Expression()
  {
  }
  double C21g::Regular(double const& x) const
  {
    return 4 * TR * ( ( pow((1-x),2) + pow(x,2) ) * log( ( 1 - x ) / x ) - 8 * x * ( x - 1 ) - 1 );
  }

  //_________________________________________________________________________________
  C22nsp::C22nsp(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C22nsp::Regular(double const& x) const
  {
    double const dl      = log(x);
    double const dl_2    = dl * dl;
    double const dl_3    = dl_2 * dl;
    double const dl1     = log(1-x);
    double const dl1_2   = dl1 * dl1;
    double const dl1_3   = dl1_2 * dl1;
    double const c2nsp2a =
      - 69.59 - 1008 * x - 2.835 * dl_3 - 17.08 * dl_2 + 5.986 * dl - 17.19 * dl1_3 + 71.08 * dl1_2 - 660.7 * dl1 - 174.8 * dl * dl1_2 + 95.09 * dl_2 * dl1
      + _nf * ( - 5.691 - 37.91 * x + 2.244 * dl_2 + 5.770 * dl - 1.707 * dl1_2  + 22.95 * dl1 + 3.036 * dl_2 * dl1 + 17.97 * dl * dl1 );
    return c2nsp2a;
  }
  double C22nsp::Singular(double const& x) const
  {
    double const dl1    = log(1-x);
    double const dl1_2  = dl1 * dl1;
    double const dl1_3  = dl1_2 * dl1;
    double const c2ns2b =
      + 14.2222 * dl1_3 - 61.3333 * dl1_2 - 31.105 * dl1 + 188.64
      + _nf * ( 1.77778 * dl1_2 - 8.5926 * dl1 + 6.3489 );
    return c2ns2b / ( 1 - x );
  }
  double C22nsp::Local(double const& x) const
  {
    double const dl1     = log(1-x);
    double const dl1_2   = dl1 * dl1;
    double const dl1_3   = dl1_2 * dl1;
    double const dl1_4   = dl1_3 * dl1;
    double const c2nsp2c =
      + 3.55555 * dl1_4 - 20.4444 * dl1_3 - 15.5525 * dl1_2 + 188.64 * dl1 - 338.531 + 0.485
      + _nf * ( 0.592593 * dl1_3 - 4.2963 * dl1_2 + 6.3489 * dl1 + 46.844 - 0.0035 );
    return c2nsp2c;
  }

  //_________________________________________________________________________________
  C22nsm::C22nsm(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C22nsm::Regular(double const& x) const
  {
    double const dl      = log(x);
    double const dl_2    = dl * dl;
    double const dl_3    = dl_2 * dl;
    double const dl1     = log(1-x);
    double const dl1_2   = dl1 * dl1;
    double const dl1_3   = dl1_2 * dl1;
    double const c2nsm2a =
      - 84.18 - 1010 * x - 3.748 * dl_3 - 19.56 * dl_2 - 1.235 * dl - 17.19 * dl1_3 + 71.08 * dl1_2 - 663.0 * dl1 - 192.4 * dl * dl1_2 + 80.41 * dl_2 * dl1
      + _nf * ( - 5.691 - 37.91 * x + 2.244 * dl_2 + 5.770 * dl - 1.707 * dl1_2  + 22.95 * dl1 + 3.036 * dl_2 * dl1 + 17.97 * dl * dl1 );
    return c2nsm2a;
  }
  double C22nsm::Singular(double const& x) const
  {
    double const dl1    = log(1-x);
    double const dl1_2  = dl1 * dl1;
    double const dl1_3  = dl1_2 * dl1;
    double const c2ns2b =
      + 14.2222 * dl1_3 - 61.3333 * dl1_2 - 31.105 * dl1 + 188.64
      + _nf * ( 1.77778 * dl1_2 - 8.5926 * dl1 + 6.3489 );
    return c2ns2b / ( 1 - x );
  }
  double C22nsm::Local(double const& x) const
  {
    double const dl1     = log(1-x);
    double const dl1_2   = dl1 * dl1;
    double const dl1_3   = dl1_2 * dl1;
    double const dl1_4   = dl1_3 * dl1;
    double const c2nsm2c =
      + 3.55555 * dl1_4 - 20.4444 * dl1_3 - 15.5525 * dl1_2 + 188.64 * dl1 - 338.531 + 0.537
      + _nf * ( 0.592593 * dl1_3 - 4.2963 * dl1_2 + 6.3489 * dl1 + 46.844 - 0.0035 );
    return c2nsm2c;
  }

  //_________________________________________________________________________________
  C22ps::C22ps():
    Expression()
  {
  }
  double C22ps::Regular(double const& x) const
  {
    double const dl     = log(x);
    double const dl_2   = dl * dl;
    double const dl_3   = dl_2 * dl;
    double const dl1    = log(1-x);
    double const dl1_2  = dl1 * dl1;
    double const dl1_3  = dl1_2 * dl1;
    double const c2ps2a =
      5.290 * ( 1 / x - 1 ) + 4.310 * dl_3 - 2.086 * dl_2 + 39.78 * dl - 0.101 * ( 1 - x ) * dl1_3 - ( 24.75 - 13.80 * x ) * dl_2 * dl1 + 30.23 * dl * dl1;
    return c2ps2a;
  }

  //_________________________________________________________________________________
  C22g::C22g():
    Expression()
  {
  }
  double C22g::Regular(double const& x) const
  {
    double const dl    = log(x);
    double const dl_2  = dl * dl;
    double const dl_3  = dl_2 * dl;
    double const dl1   = log(1-x);
    double const dl1_2 = dl1 * dl1;
    double const dl1_3 = dl1_2 * dl1;
    double const c2g2a =
      1 / x * ( 11.90 + 1494.* dl1 ) + 5.319 * dl_3 - 59.48 * dl_2 - 284.8 * dl + 392.4 - 1483 * dl1
      + ( 6.445 + 209.4 * ( 1 - x ) ) * dl1_3 - 24.00 * dl1_2 - 724.1 * dl_2 * dl1 - 871.8 * dl * dl1_2;
    return c2g2a;
  }
  double C22g::Local(double const&) const
  {
    return - 0.28;
  }

  //_________________________________________________________________________________
  CL1ns::CL1ns():
    Expression()
  {
  }
  double CL1ns::Regular(double const& x) const
  {
    return 4 * CF * x;
  }

  //_________________________________________________________________________________
  CL1g::CL1g():
    Expression()
  {
  }
  double CL1g::Regular(double const& x) const
  {
    return 16 * TR * x * ( 1 - x );
  }

  //_________________________________________________________________________________
  CL2nsp::CL2nsp(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double CL2nsp::Regular(double const& x) const
  {
    double const dl      = log(x);
    double const dl_2    = dl * dl;
    double const dl1     = log(1-x);
    double const dl1_2   = dl1 * dl1;
    double const clnsp2a =
      - 40.41 + 97.48 * x + ( 26.56 * x - 0.031 ) * dl_2 - 14.85 * dl + 13.62 * dl1_2 - 55.79 * dl1 - 150.5 * dl * dl1
      + _nf * 16 / 27. * ( 6 * x * dl1 - 12 * x * dl - 25 * x + 6 );
    return clnsp2a;
  }
  double CL2nsp::Local(double const&) const
  {
    return - 0.164;
  }

  //_________________________________________________________________________________
  CL2nsm::CL2nsm(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double CL2nsm::Regular(double const& x) const
  {
    double const dl      = log(x);
    double const dl_2    = dl * dl;
    double const dl1     = log(1-x);
    double const dl1_2   = dl1 * dl1;
    double const clnsm2a =
      - 52.27 + 100.8 * x + ( 23.29 * x - 0.043 ) * dl_2 - 22.21 * dl + 13.30 * dl1_2 - 59.12 * dl1 - 141.7 * dl * dl1
      + _nf * 16 / 27. * ( 6 * x * dl1 - 12 * x * dl - 25 * x + 6 );
    return clnsm2a;
  }
  double CL2nsm::Local(double const&) const
  {
    return - 0.150;
  }

  //_________________________________________________________________________________
  CL2ps::CL2ps():
    Expression()
  {
  }
  double CL2ps::Regular(double const& x) const
  {
    double const dl     = log(x);
    double const dl_2   = dl * dl;
    double const dl1    = log(1-x);
    double const omx    = 1 - x;
    double const omx2   = omx * omx;
    double const omx3   = omx2 * omx;
    double const clps2a =
      ( 15.94 - 5.212 * x ) * omx2 * dl1 + ( 0.421 + 1.520 * x ) * dl_2 + 28.09 * omx * dl - ( 2.370 / x - 19.27 ) * omx3;
    return clps2a;
  }

  //_________________________________________________________________________________
  CL2g::CL2g():
    Expression()
  {
  }
  double CL2g::Regular(double const& x) const
  {
    double const dl    = log(x);
    double const dl_2  = dl * dl;
    double const dl1   = log(1-x);
    double const dl1_2 = dl1 * dl1;
    double const omx   = 1 - x;
    double const clg2a =
      ( 94.74 - 49.20 * x ) * omx * dl1_2 + 864.8 * omx * dl1 + 1161 * x * dl * dl1 + 60.06 * x * dl_2 + 39.66 * omx * dl - 5.333 * ( 1 / x - 1 );
    return clg2a;
  }

  //_________________________________________________________________________________
  C31ns::C31ns():
    Expression()
  {
  }
  double C31ns::Regular(double const& x) const
  {
    return 2 * CF * ( - ( 1 + x ) * log( 1 - x ) - ( 1 + pow(x,2) ) * log(x) / ( 1 - x ) + 2 + x );
  }
  double C31ns::Singular(double const& x) const
  {
    return 2 * CF * ( 2 * log( 1 - x ) - 3 / 2. ) / ( 1 - x );
  }
  double C31ns::Local(double const& x) const
  {
    return 2 * CF * ( pow(log(1-x),2) - 3 * log( 1 - x ) / 2 - ( 2 * zeta2 + 9 / 2. ) );
  }

  //_________________________________________________________________________________
  C32nsp::C32nsp(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C32nsp::Regular(double const& x) const
  {
    double const dl      = log(x);
    double const dl_2    = dl * dl;
    double const dl_3    = dl_2 * dl;
    double const dl1     = log(1-x);
    double const dl1_2   = dl1 * dl1;
    double const dl1_3   = dl1_2 * dl1;
    double const c3nsp2a =
      - 242.9 - 467.2 * x - 3.049 * dl_3 - 30.14 * dl_2 - 79.14 * dl - 15.20 * dl1_3 + 94.61 * dl1_2 - 396.1 * dl1 - 92.43 * dl * dl1_2
      + _nf * ( - 6.337 - 14.97 * x + 2.207 * dl_2 + 8.683 * dl + 0.042 * dl1_3 - 0.808 * dl1_2  + 25.00 * dl1 + 9.684 * dl * dl1 );
    return c3nsp2a;
  }
  double C32nsp::Singular(double const& x) const
  {
    double const dl1    = log(1-x);
    double const dl1_2  = dl1 * dl1;
    double const dl1_3  = dl1_2 * dl1;
    double const c3ns2b =
      + 14.2222 * dl1_3 - 61.3333 * dl1_2 - 31.105 * dl1 + 188.64
      + _nf * ( 1.77778 * dl1_2 - 8.5926 * dl1 + 6.3489 );
    return c3ns2b / ( 1 - x );
  }
  double C32nsp::Local(double const& x) const
  {
    double const dl1     = log(1-x);
    double const dl1_2   = dl1 * dl1;
    double const dl1_3   = dl1_2 * dl1;
    double const dl1_4   = dl1_3 * dl1;
    double const c3nsp2c =
      + 3.55555 * dl1_4 - 20.4444 * dl1_3 - 15.5525 * dl1_2 + 188.64 * dl1 - 338.531 - 0.152
      + _nf * ( 0.592593 * dl1_3 - 4.2963 * dl1_2 + 6.3489 * dl1 + 46.844 + 0.013 );
    return c3nsp2c;
  }

  //_________________________________________________________________________________
  C32nsm::C32nsm(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C32nsm::Regular(double const& x) const
  {
    double const dl      = log(x);
    double const dl_2    = dl * dl;
    double const dl_3    = dl_2 * dl;
    double const dl1     = log(1-x);
    double const dl1_2   = dl1 * dl1;
    double const dl1_3   = dl1_2 * dl1;
    double const c3nsm2a =
      - 206.1 - 576.8 * x - 3.922 * dl_3 - 33.31 * dl_2 - 67.60 * dl - 15.20 * dl1_3 + 94.61 * dl1_2 - 409.6 * dl1 - 147.9 * dl * dl1_2
      + _nf * ( - 6.337 - 14.97 * x + 2.207 * dl_2 + 8.683 * dl + 0.042 * dl1_3 - 0.808 * dl1_2 + 25.00 * dl1 + 9.684 * dl * dl1 );
    return c3nsm2a;
  }
  double C32nsm::Singular(double const& x) const
  {
    double const dl1    = log(1-x);
    double const dl1_2  = dl1 * dl1;
    double const dl1_3  = dl1_2 * dl1;
    double const c3ns2b =
      + 14.2222 * dl1_3 - 61.3333 * dl1_2 - 31.105 * dl1 + 188.64
      + _nf * ( 1.77778 * dl1_2 - 8.5926 * dl1 + 6.3489 );
    return c3ns2b / ( 1 - x );
  }
  double C32nsm::Local(double const& x) const
  {
    double const dl1     = log(1-x);
    double const dl1_2   = dl1 * dl1;
    double const dl1_3   = dl1_2 * dl1;
    double const dl1_4   = dl1_3 * dl1;
    double const c3nsm2c =
      + 3.55555 * dl1_4 - 20.4444 * dl1_3 - 15.5525 * dl1_2 + 188.64 * dl1 - 338.531 - 0.104
      + _nf * ( 0.592593 * dl1_3 - 4.2963 * dl1_2 + 6.3489 * dl1 + 46.844 + 0.013 );
    return c3nsm2c;
  }
}
