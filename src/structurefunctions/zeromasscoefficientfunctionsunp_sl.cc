//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/zeromasscoefficientfunctionsunp_sl.h"
#include "apfel/constants.h"
#include "apfel/specialfunctions.h"

namespace apfel
{
  //_________________________________________________________________________________
  C21ns::C21ns():
    Expression()
  {
  }
  double C21ns::Regular(double const& x) const
  {
    return 2 * CF * ( - ( 1 + x ) * log(1 - x) - ( 1 + pow(x, 2) ) * log(x) / ( 1 - x ) + 3 + 2 * x );
  }
  double C21ns::Singular(double const& x) const
  {
    return 2 * CF * ( 2 * log(1 - x) - 3 / 2. ) / ( 1 - x );
  }
  double C21ns::Local(double const& x) const
  {
    return 2 * CF * ( pow(log(1 - x), 2) - 3 * log(1 - x) / 2 - ( 2 * zeta2 + 9 / 2. ) );
  }

  //_________________________________________________________________________________
  C21g::C21g():
    Expression()
  {
  }
  double C21g::Regular(double const& x) const
  {
    return 4 * TR * ( ( pow(1 - x, 2) + pow(x, 2) ) * log( ( 1 - x ) / x ) - 8 * x * ( x - 1 ) - 1 );
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
  C31ns::C31ns():
    Expression()
  {
  }
  double C31ns::Regular(double const& x) const
  {
    return 2 * CF * ( - ( 1 + x ) * log(1 - x) - ( 1 + pow(x, 2) ) * log(x) / ( 1 - x ) + 2 + x );
  }
  double C31ns::Singular(double const& x) const
  {
    return 2 * CF * ( 2 * log(1 - x) - 3 / 2. ) / ( 1 - x );
  }
  double C31ns::Local(double const& x) const
  {
    return 2 * CF * ( pow(log(1 - x), 2) - 3 * log(1 - x) / 2 - ( 2 * zeta2 + 9 / 2. ) );
  }

  //_________________________________________________________________________________
  C22nsp::C22nsp(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C22nsp::Regular(double const& x) const
  {
    double const dl    = log(x);
    double const dl_2  = dl * dl;
    double const dl_3  = dl_2 * dl;
    double const dl1   = log(1-x);
    double const dl1_2 = dl1 * dl1;
    double const dl1_3 = dl1_2 * dl1;
    return
      - 69.59 - 1008 * x - 2.835 * dl_3 - 17.08 * dl_2 + 5.986 * dl - 17.19 * dl1_3 + 71.08 * dl1_2 - 660.7 * dl1 - 174.8 * dl * dl1_2 + 95.09 * dl_2 * dl1
      + _nf * ( - 5.691 - 37.91 * x + 2.244 * dl_2 + 5.770 * dl - 1.707 * dl1_2  + 22.95 * dl1 + 3.036 * dl_2 * dl1 + 17.97 * dl * dl1 );
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
    double const dl1   = log(1-x);
    double const dl1_2 = dl1 * dl1;
    double const dl1_3 = dl1_2 * dl1;
    double const dl1_4 = dl1_3 * dl1;
    return
      + 3.55555 * dl1_4 - 20.4444 * dl1_3 - 15.5525 * dl1_2 + 188.64 * dl1 - 338.531 + 0.485
      + _nf * ( 0.592593 * dl1_3 - 4.2963 * dl1_2 + 6.3489 * dl1 + 46.844 - 0.0035 );
  }

  //_________________________________________________________________________________
  C22nsm::C22nsm(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C22nsm::Regular(double const& x) const
  {
    double const dl    = log(x);
    double const dl_2  = dl * dl;
    double const dl_3  = dl_2 * dl;
    double const dl1   = log(1-x);
    double const dl1_2 = dl1 * dl1;
    double const dl1_3 = dl1_2 * dl1;
    return
      - 84.18 - 1010 * x - 3.748 * dl_3 - 19.56 * dl_2 - 1.235 * dl - 17.19 * dl1_3 + 71.08 * dl1_2 - 663.0 * dl1 - 192.4 * dl * dl1_2 + 80.41 * dl_2 * dl1
      + _nf * ( - 5.691 - 37.91 * x + 2.244 * dl_2 + 5.770 * dl - 1.707 * dl1_2  + 22.95 * dl1 + 3.036 * dl_2 * dl1 + 17.97 * dl * dl1 );
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
    double const dl1   = log(1-x);
    double const dl1_2 = dl1 * dl1;
    double const dl1_3 = dl1_2 * dl1;
    double const dl1_4 = dl1_3 * dl1;
    return
      + 3.55555 * dl1_4 - 20.4444 * dl1_3 - 15.5525 * dl1_2 + 188.64 * dl1 - 338.531 + 0.537
      + _nf * ( 0.592593 * dl1_3 - 4.2963 * dl1_2 + 6.3489 * dl1 + 46.844 - 0.0035 );
  }

  //_________________________________________________________________________________
  C22ps::C22ps():
    Expression()
  {
  }
  double C22ps::Regular(double const& x) const
  {
    double const dl    = log(x);
    double const dl_2  = dl * dl;
    double const dl_3  = dl_2 * dl;
    double const dl1   = log(1-x);
    double const dl1_2 = dl1 * dl1;
    double const dl1_3 = dl1_2 * dl1;
    return
      5.290 * ( 1 / x - 1 ) + 4.310 * dl_3 - 2.086 * dl_2 + 39.78 * dl - 0.101 * ( 1 - x ) * dl1_3 - ( 24.75 - 13.80 * x ) * dl_2 * dl1 + 30.23 * dl * dl1;
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
    return
      1 / x * ( 11.90 + 1494.* dl1 ) + 5.319 * dl_3 - 59.48 * dl_2 - 284.8 * dl + 392.4 - 1483 * dl1
      + ( 6.445 + 209.4 * ( 1 - x ) ) * dl1_3 - 24.00 * dl1_2 - 724.1 * dl_2 * dl1 - 871.8 * dl * dl1_2;
  }
  double C22g::Local(double const&) const
  {
    return - 0.28;
  }

  //_________________________________________________________________________________
  CL2nsp::CL2nsp(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double CL2nsp::Regular(double const& x) const
  {
    double const dl    = log(x);
    double const dl_2  = dl * dl;
    double const dl1   = log(1-x);
    double const dl1_2 = dl1 * dl1;
    return
      - 40.41 + 97.48 * x + ( 26.56 * x - 0.031 ) * dl_2 - 14.85 * dl + 13.62 * dl1_2 - 55.79 * dl1 - 150.5 * dl * dl1
      + _nf * 16 / 27. * ( 6 * x * dl1 - 12 * x * dl - 25 * x + 6 );
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
    double const dl    = log(x);
    double const dl_2  = dl * dl;
    double const dl1   = log(1-x);
    double const dl1_2 = dl1 * dl1;
    return
      - 52.27 + 100.8 * x + ( 23.29 * x - 0.043 ) * dl_2 - 22.21 * dl + 13.30 * dl1_2 - 59.12 * dl1 - 141.7 * dl * dl1
      + _nf * 16 / 27. * ( 6 * x * dl1 - 12 * x * dl - 25 * x + 6 );
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
    double const dl   = log(x);
    double const dl_2 = dl * dl;
    double const dl1  = log(1-x);
    double const omx  = 1 - x;
    double const omx2 = omx * omx;
    double const omx3 = omx2 * omx;
    return
      ( 15.94 - 5.212 * x ) * omx2 * dl1 + ( 0.421 + 1.520 * x ) * dl_2 + 28.09 * omx * dl - ( 2.370 / x - 19.27 ) * omx3;
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
    return
      ( 94.74 - 49.20 * x ) * omx * dl1_2 + 864.8 * omx * dl1 + 1161 * x * dl * dl1 + 60.06 * x * dl_2 + 39.66 * omx * dl - 5.333 * ( 1 / x - 1 );
  }

  //_________________________________________________________________________________
  C32nsp::C32nsp(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C32nsp::Regular(double const& x) const
  {
    double const dl    = log(x);
    double const dl_2  = dl * dl;
    double const dl_3  = dl_2 * dl;
    double const dl1   = log(1-x);
    double const dl1_2 = dl1 * dl1;
    double const dl1_3 = dl1_2 * dl1;
    return
      - 242.9 - 467.2 * x - 3.049 * dl_3 - 30.14 * dl_2 - 79.14 * dl - 15.20 * dl1_3 + 94.61 * dl1_2 - 396.1 * dl1 - 92.43 * dl * dl1_2
      + _nf * ( - 6.337 - 14.97 * x + 2.207 * dl_2 + 8.683 * dl + 0.042 * dl1_3 - 0.808 * dl1_2  + 25.00 * dl1 + 9.684 * dl * dl1 );
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
    double const dl1   = log(1-x);
    double const dl1_2 = dl1 * dl1;
    double const dl1_3 = dl1_2 * dl1;
    double const dl1_4 = dl1_3 * dl1;
    return
      + 3.55555 * dl1_4 - 20.4444 * dl1_3 - 15.5525 * dl1_2 + 188.64 * dl1 - 338.531 - 0.152
      + _nf * ( 0.592593 * dl1_3 - 4.2963 * dl1_2 + 6.3489 * dl1 + 46.844 + 0.013 );
  }

  //_________________________________________________________________________________
  C32nsm::C32nsm(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C32nsm::Regular(double const& x) const
  {
    double const dl    = log(x);
    double const dl_2  = dl * dl;
    double const dl_3  = dl_2 * dl;
    double const dl1   = log(1-x);
    double const dl1_2 = dl1 * dl1;
    double const dl1_3 = dl1_2 * dl1;
    return
      - 206.1 - 576.8 * x - 3.922 * dl_3 - 33.31 * dl_2 - 67.60 * dl - 15.20 * dl1_3 + 94.61 * dl1_2 - 409.6 * dl1 - 147.9 * dl * dl1_2
      + _nf * ( - 6.337 - 14.97 * x + 2.207 * dl_2 + 8.683 * dl + 0.042 * dl1_3 - 0.808 * dl1_2 + 25.00 * dl1 + 9.684 * dl * dl1 );
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
    double const dl1   = log(1-x);
    double const dl1_2 = dl1 * dl1;
    double const dl1_3 = dl1_2 * dl1;
    double const dl1_4 = dl1_3 * dl1;
    return
      + 3.55555 * dl1_4 - 20.4444 * dl1_3 - 15.5525 * dl1_2 + 188.64 * dl1 - 338.531 - 0.104
      + _nf * ( 0.592593 * dl1_3 - 4.2963 * dl1_2 + 6.3489 * dl1 + 46.844 + 0.013 );
  }

  //_________________________________________________________________________________
  C23nsp::C23nsp(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C23nsp::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double x3   = x * x2;
    const double dl   = log(x);
    const double dl2  = dl * dl;
    const double dl3  = dl * dl2;
    const double dl4  = dl * dl3;
    const double dl5  = dl * dl4;
    const double x1   = 1 - x;
    const double dl1  = log(x1);
    const double dl12 = dl1 * dl1;
    const double dl13 = dl1 * dl12;
    const double dl14 = dl1 * dl13;
    const double dl15 = dl1 * dl14;
    const double d27  = 1. / 27;
    const double d243 = 1. / 243;
    return
      - 4926. + 7725. * x + 57256. * x2 + 12898. * x3 - 32. * d27 * dl5 - 8796. * d243 * dl4 - 309.1 * dl3
      - 899.6 * dl2 - 775.8 * dl + 4.719 * x * dl5 - 512. * d27 * dl15 + 6336. * d27 * dl14
      - 3368. * dl13 - 2978. * dl12 + 18832. * dl1 - 56000. * x1 * dl12 - dl*dl1 * ( 6158. + 1836. * dl )
      + _nf * ( 831.6 - 6752. * x - 2778. * x2 + 728. * d243 * dl4 + 12224. * d243 * dl3 + 187.3 * dl2 + 275.6 * dl
                + 4.102 * x * dl4 - 1920. * d243 * dl14 + 153.5 * dl13 - 828.7 * dl12 - 501.1 * dl1 + 171.0 * x1 * dl14
                + dl * dl1 * ( 4365. + 716.2 * dl - 5983. * dl1 ) )
      + _nf * _nf * ( 129.2 * x + 102.5 * x2 - 368. * d243 * dl3 - 1984. * d243 * dl2 - 8.042 * dl
                      - 192. * d243 * dl13 + 18.21 * dl12 - 19.09 * dl1 + dl * dl1 * ( - 96.07 - 12.46 * dl + 85.88 * dl1 ) )
      + fl11ns[_nf-1] * _nf * ( ( 126.42 - 50.29 * x - 50.15 * x2 ) * x1 - 26.717 - 960. * d243 * dl2 * ( dl + 5 ) + 59.59 * dl
                                - x * dl2 * ( 101.8 + 34.79 * dl + 3.070 * dl2 ) - 9.075 * x * x1 * dl1 ) * x;
  }
  double C23nsp::Singular(double const& x) const
  {
    const double dl1  = log(1-x);
    const double dl12 = dl1 * dl1;
    const double dl13 = dl1 * dl12;
    const double dl14 = dl1 * dl13;
    const double dl15 = dl1 * dl14;
    const double d81  = 1. / 81;
    const double c2ns3b =
      + 1536. * d81 * dl15 - 16320. * d81 * dl14 + 5.01099e+2 * dl13 + 1.17154e+3 * dl12 - 7.32845e+3 * dl1 + 4.44276e+3
      + _nf * ( 640. * d81 * dl14 - 6592. * d81 * dl13 + 220.573 * dl12 + 294.906 * dl1 - 729.359 )
      + _nf * _nf * ( 64. * d81 * dl13 - 464. * d81 * dl12 + 7.67505 * dl1 + 1.00830 );
    return c2ns3b / ( 1 - x );
  }
  double C23nsp::Local(double const& x) const
  {
    const double dl1  = log(1-x);
    const double dl12 = dl1 * dl1;
    const double dl13 = dl1 * dl12;
    const double dl14 = dl1 * dl13;
    const double dl15 = dl1 * dl14;
    const double dl16 = dl1 * dl15;
    const double d3   = 1. / 3;
    const double d81  = 1. / 81;
    return
      + 256. * d81 * dl16 - 3264. * d81 * dl15 + 1.252745e+2 * dl14 + 3.905133e+2 * dl13 - 3.664225e+3 * dl12 + 4.44276e+3  * dl1 - 9195.48 + 25.10
      + _nf * ( 128. * d81 * dl15 - 1648. * d81 * dl14 + 220.573 * d3 * dl13 + 147.453 * dl12 - 729.359 * dl1 + 2575.074 - 0.387 )
      + _nf * _nf * ( 16. * d81 * dl14 - 464. * d81 * d3 * dl13 + 7.67505 * 0.5 * dl12 + 1.0083 * dl1 - 103.2521 + 0.0155 )
      - fl11ns[_nf-1] * _nf * 11.8880;
  }

  //_________________________________________________________________________________
  C23nsm::C23nsm(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C23nsm::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double x3   = x * x2;
    const double dl   = log(x);
    const double dl2  = dl * dl;
    const double dl3  = dl * dl2;
    const double dl4  = dl * dl3;
    const double dl5  = dl * dl4;
    const double x1   = 1 - x;
    const double dl1  = log(x1);
    const double dl12 = dl1 * dl1;
    const double dl13 = dl1 * dl12;
    const double dl14 = dl1 * dl13;
    const double dl15 = dl1 * dl14;
    const double d27  = 1. / 27;
    const double d243 = 1. / 243;

    // Plus piece
    const double C23nsp =
      - 4926. + 7725. * x + 57256. * x2 + 12898. * x3 - 32. * d27 * dl5 - 8796. * d243 * dl4 - 309.1 * dl3
      - 899.6 * dl2 - 775.8 * dl + 4.719 * x * dl5 - 512. * d27 * dl15 + 6336. * d27 * dl14
      - 3368. * dl13 - 2978. * dl12 + 18832. * dl1 - 56000. * x1 * dl12 - dl*dl1 * ( 6158. + 1836. * dl )
      + _nf * ( 831.6 - 6752. * x - 2778. * x2 + 728. * d243 * dl4 + 12224. * d243 * dl3 + 187.3 * dl2 + 275.6 * dl
                + 4.102 * x * dl4 - 1920. * d243 * dl14 + 153.5 * dl13 - 828.7 * dl12 - 501.1 * dl1 + 171.0 * x1 * dl14
                + dl * dl1 * ( 4365. + 716.2 * dl - 5983. * dl1 ) )
      + _nf * _nf * ( 129.2 * x + 102.5 * x2 - 368. * d243 * dl3 - 1984. * d243 * dl2 - 8.042 * dl
                      - 192. * d243 * dl13 + 18.21 * dl12 - 19.09 * dl1 + dl * dl1 * ( - 96.07 - 12.46 * dl + 85.88 * dl1 ) )
      + fl11ns[_nf-1] * _nf * ( ( 126.42 - 50.29 * x - 50.15 * x2 ) * x1 - 26.717 - 960. * d243 * dl2 * ( dl + 5 ) + 59.59 * dl
                                - x * dl2 * ( 101.8 + 34.79 * dl + 3.070 * dl2 ) - 9.075 * x * x1 * dl1 ) * x;

    // Compute and include difference to get the minus piece
    const double c2q30a = ( 54.478 * dl12 + 304.60 * dl1 + 691.68 * x ) * x1 + 179.14 * dl - 0.1826 * dl3;
    const double c2q30b = - ( 13.378 * dl12 + 97.60 * dl1 + 118.12 * x ) * x1 - 91.196 * dl2 - 0.4644 * dl5;
    const double c2q31a = ( 20.822 * x2 - 282.10 * ( 1. + x / 2. ) ) * x1 - 285.58 * x * dl - 112.30 * dl + 3.587 * dl3;
    const double c2q31b = ( 4.522 * dl1 + 447.88 * ( 1. + x / 2. ) ) * x1 + 514.02 * x * dl + 147.05 * dl + 7.386 * dl2;
    return C23nsp - 0.5 * ( c2q30a + c2q30b + _nf * ( c2q31a + c2q31b ) );
  }
  double C23nsm::Singular(double const& x) const
  {
    const double dl1  = log(1-x);
    const double dl12 = dl1 * dl1;
    const double dl13 = dl1 * dl12;
    const double dl14 = dl1 * dl13;
    const double dl15 = dl1 * dl14;
    const double d81  = 1. / 81;
    const double c2ns3b =
      + 1536. * d81 * dl15 - 16320. * d81 * dl14 + 5.01099e+2 * dl13 + 1.17154e+3 * dl12 - 7.32845e+3 * dl1 + 4.44276e+3
      + _nf * ( 640. * d81 * dl14 - 6592. * d81 * dl13 + 220.573 * dl12 + 294.906 * dl1 - 729.359 )
      + _nf * _nf * ( 64. * d81 * dl13 - 464. * d81 * dl12 + 7.67505 * dl1 + 1.00830 );
    return c2ns3b / ( 1 - x );
  }
  double C23nsm::Local(double const& x) const
  {
    const double dl1  = log(1-x);
    const double dl12 = dl1 * dl1;
    const double dl13 = dl1 * dl12;
    const double dl14 = dl1 * dl13;
    const double dl15 = dl1 * dl14;
    const double dl16 = dl1 * dl15;
    const double d3   = 1. / 3;
    const double d81  = 1. / 81;
    return
      + 256. * d81 * dl16 - 3264. * d81 * dl15 + 1.252745e+2 * dl14 + 3.905133e+2 * dl13 - 3.664225e+3 * dl12 + 4.44276e+3  * dl1 - 9195.48 + 25.10
      + _nf * ( 128. * d81 * dl15 - 1648. * d81 * dl14 + 220.573 * d3 * dl13 + 147.453 * dl12 - 729.359 * dl1 + 2575.074 - 0.387 )
      + _nf * _nf * ( 16. * d81 * dl14 - 464. * d81 * d3 * dl13 + 7.67505 * 0.5 * dl12 + 1.0083 * dl1 - 103.2521 + 0.0155 )
      - fl11ns[_nf-1] * _nf * 11.8880;
  }

  //_________________________________________________________________________________
  C23ps::C23ps(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C23ps::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double dl    = log(x);
    const double dl2   = dl * dl;
    const double dl3   = dl * dl2;
    const double dl4   = dl * dl3;
    const double dl5   = dl * dl4;
    const double x1    = 1 - x;
    const double dl1   = log(x1);
    const double dl12  = dl1 * dl1;
    const double dl13  = dl1 * dl12;
    const double dl14  = dl1 * dl13;
    const double d9    = 1. / 9;
    const double d81   = 1. / 81;
    const double c2s31 = ( 856. * d81 * dl14 - 6032. * d81 * dl13 + 130.57 * dl12 - 542.0 * dl1 + 8501. - 4714. * x + 61.50 * x2 ) * x1
                         + dl * dl1 * ( 8831. * dl + 4162. * x1 ) - 15.44 * x * dl5 + 3333. * x * dl2 + 1615. * dl + 1208. * dl2  - 333.73 * dl3
                         + 4244. * d81 * dl4 - 40. * d9 * dl5 - 2731.82 * x1 / x - 414.262 * dl / x;
    const double c2s32 = ( - 64. * d81 * dl13 + 208. * d81 * dl12 + 23.09 * dl1 - 220.27 + 59.80 * x - 177.6 * x2) * x1
                         - dl * dl1 * ( 160.3 * dl + 135.4 * x1 ) - 24.14 * x * dl3 - 215.4 * x * dl2 - 209.8 * dl - 90.38 * dl2
                         - 3568. / 243. * dl3 - 184. * d81 * dl4 + 40.2426 * x1 / x;
    const double c2s3f = ( ( 126.42 - 50.29 * x - 50.15 * x2) * x1 - 26.717 - 320. * d81 * dl2 * ( dl + 5. ) + 59.59 * dl
                           - x * dl2 * ( 101.8 + 34.79 * dl + 3.070 * dl2 ) - 9.075 * x * x1 * dl1 ) * x;
    return c2s31 + _nf * c2s32 + ( fl11sg[_nf-1] - fl11ns[_nf-1] ) * c2s3f;
  }
  double C23ps::Local(double const&) const
  {
    return - ( fl11sg[_nf-1] - fl11ns[_nf-1] ) * 11.8880;
  }

  //_________________________________________________________________________________
  C23g::C23g(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C23g::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double x3    = x * x2;
    const double dl    = log(x);
    const double dl2   = dl * dl;
    const double dl3   = dl * dl2;
    const double dl4   = dl * dl3;
    const double dl5   = dl * dl4;
    const double x1    = 1 - x;
    const double dl1   = log(x1);
    const double dl12  = dl1 * dl1;
    const double dl13  = dl1 * dl12;
    const double dl14  = dl1 * dl13;
    const double dl15  = dl1 * dl14;
    const double d9    = 1. / 9;
    const double d81   = 1. / 81;
    const double c2g31 =
      966. * d81 * dl15 - 935.5 * d9 * dl14 + 89.31 * dl13 + 979.2 * dl12 - 2405. * dl1 + 1372. * x1 * dl14 - 15729. - 310510. * x
      + 331570. * x2 - 244150. * x * dl2 - 253.3 * x * dl5 + dl * dl1 * ( 138230. - 237010. * dl ) - 11860. * dl - 700.8 * dl2
      - 1440. * dl3 + 2480.5 * d81 * dl4 - 134. * d9 * dl5 - 6362.54 / x - 932.089 * dl / x;
    const double c2g32 =
      131. * d81 * dl14 - 14.72 * dl13 + 3.607 * dl12 - 226.1 * dl1 + 4.762 - 190.0 * x - 818.4 * x2 - 4019. * x * dl2
      - dl * dl1 * ( 791.5 + 4646 * dl ) + 739.0 * dl + 418.0 * dl2 + 104.3 * dl3 + 809. * d81 * dl4 + 12. * d9 * dl5 + 84.423 / x;
    const double c2g3f =
      3.211 * dl12 + 19.04 * x * dl1 + 0.623 * x1 * dl13 - 64.47 * x + 121.6 * x2 - 45.82 * x3 - x * dl * dl1 * ( 31.68 + 37.24 * dl )
      - x * dl * ( 82.40 + 16.08 * dl ) + x * dl3 * ( 520. * d81 + 11.27 * x ) + 60. * d81 * x * dl4;
    return c2g31 + _nf * ( c2g32 + fl11sg[_nf-1] * c2g3f );
  }
  double C23g::Local(double const&) const
  {
    return 0.625;
  }

  //_________________________________________________________________________________
  CL3nsp::CL3nsp(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double CL3nsp::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double dl    = log(x);
    const double dl2   = dl * dl;
    const double dl3   = dl * dl2;
    const double dl1   = log(1-x);
    const double dl12  = dl1 * dl1;
    const double dl13  = dl1 * dl12;
    const double dl14  = dl1 * dl13;
    const double d81   = 1. / 81;
    return
      - 2220.5 - 7884. * x + 4168. * x2 - 1280. * d81 * dl3 - 7456. / 27. * dl2 - 1355.7 * dl + 512. / 27. * dl14 - 177.40 * dl13 + 650.6 * dl12
      - 2729. * dl1 + 208.3 * x * dl3 - dl13 * ( 1. - x ) * ( 125.3 - 195.6 * dl1) - dl* dl1 * ( 844.7 * dl + 517.3 * dl1 )
      + _nf * ( 408.4 - 9.345 * x - 919.3 * x2 + 1728. * d81 * dl2 + 200.73 * dl - 1792. * d81 * x * dl3 + 1024. * d81 * dl13 - 112.35 * dl12 + 344.1 * dl1
                + ( 1. - x ) * dl12 * ( 239.7 + 20.63 * dl1 ) + dl* dl1 * ( 887.3 + 294.5 * dl - 59.14 * dl1 ) )
      + _nf * _nf * ( - 19. + ( 317. / 6. - 12. * zeta2 ) * x + 9. * x * dl2 + dl * ( - 6. + 50. * x ) + 3. * x * dl12 + dl1 * ( 6. - 25. * x )
                      - 6. * x * dl* dl1 + 6. * x * dilog(x) ) * 64. * d81
      + fl11ns[_nf-1] * _nf * ( ( 107.0 + 321.05 * x - 54.62 * x2 ) * ( 1. - x ) - 26.717 - 320. * d81 * dl3 - 640. * d81 * dl2
                                + 9.773 * dl + x * dl * ( 363.8 + 68.32 * dl ) ) * x;
  }
  double CL3nsp::Local(double const&) const
  {
    return 0.113 + _nf * 0.006;
  }

  //_________________________________________________________________________________
  CL3nsm::CL3nsm(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double CL3nsm::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double x3    = x * x2;
    const double dl    = log(x);
    const double dl2   = dl * dl;
    const double dl3   = dl * dl2;
    const double x1    = 1 - x;
    const double x12   = x1 * x1;
    const double dl1   = log(x1);
    const double dl12  = dl1 * dl1;
    const double dl13  = dl1 * dl12;
    const double dl14  = dl1 * dl13;
    const double d81   = 1. / 81;

    // Plus piece
    const double CL3nsp =
      - 2220.5 - 7884. * x + 4168. * x2 - 1280. * d81 * dl3 - 7456. / 27. * dl2 - 1355.7 * dl + 512. / 27. * dl14 - 177.40 * dl13 + 650.6 * dl12
      - 2729. * dl1 + 208.3 * x * dl3 - dl13 * ( 1. - x ) * ( 125.3 - 195.6 * dl1) - dl* dl1 * ( 844.7 * dl + 517.3 * dl1 )
      + _nf * ( 408.4 - 9.345 * x - 919.3 * x2 + 1728. * d81 * dl2 + 200.73 * dl - 1792. * d81 * x * dl3 + 1024. * d81 * dl13 - 112.35 * dl12 + 344.1 * dl1
                + ( 1. - x ) * dl12 * ( 239.7 + 20.63 * dl1 ) + dl* dl1 * ( 887.3 + 294.5 * dl - 59.14 * dl1 ) )
      + _nf * _nf * ( - 19. + ( 317. / 6. - 12. * zeta2 ) * x + 9. * x * dl2 + dl * ( - 6. + 50. * x ) + 3. * x * dl12 + dl1 * ( 6. - 25. * x )
                      - 6. * x * dl* dl1 + 6. * x * dilog(x) ) * 64. * d81
      + fl11ns[_nf-1] * _nf * ( ( 107.0 + 321.05 * x - 54.62 * x2 ) * ( 1. - x ) - 26.717 - 320. * d81 * dl3 - 640. * d81 * dl2
                                + 9.773 * dl + x * dl * ( 363.8 + 68.32 * dl ) ) * x;

    // Compute and include difference to get the minus piece
    const double clq30a = - ( 495.49 * x2 + 906.86 ) * x12 - 983.23 * x * x1 * dl + 53.706 * dl2 + 5.3059 * dl3;
    const double clq30b = ( 78.306 * dl1 + 6.3838 * x ) * x12 + 20.809 * x * x1 * dl - 114.47 * dl2 - 22.222 * dl3;
    const double clq31a = ( 29.95 * x3 - 59.087 * x2 + 379.91 ) * x12 - 273.042 * x * dl2 + 71.482 * x1 * dl;
    const double clq31b = ( 12.532 * dl1 + 141.99 * x2 - 250.62 * x ) * x12 - ( 153.586 * x - 0.6569 ) * x1 * dl;
    return CL3nsp - 0.5 * ( clq30a + clq30b + _nf * ( clq31a + clq31b ) );
  }
  double CL3nsm::Local(double const&) const
  {
    return 0.113 + _nf * 0.006;
  }

  //_________________________________________________________________________________
  CL3ps::CL3ps(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double CL3ps::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double dl    = log(x);
    const double dl2   = dl * dl;
    const double dl3   = dl * dl2;
    const double x1    = 1 - x;
    const double x12   = x1 * x1;
    const double dl1   = log(x1);
    const double dl12  = dl1 * dl1;
    const double dl13  = dl1 * dl12;
    const double d27   = 1. / 27;
    const double d81   = 1. / 81;
    const double cls31 = ( 1568. * d27 * dl13 - 11904. * d27 * dl12 + 5124. * dl1 ) * x12 + dl * dl1 * ( 2184. * dl + 6059. * x1 )
                         - ( 795.6 + 1036. * x ) * x12 - 143.6 * dl * x1 + 8544. * d27 * dl2 - 1600. * d27 * dl3 - 885.53 / x * x12 - 182. * dl / x * x1;
    const double cls32 = ( - 96. * d27 * dl12 + 29.52 * dl1) * x12 + dl * dl1 * ( 35.18 * dl + 73.06 * x1 ) - ( 14.16 - 69.84 * x ) * x12
                         - 35.24 * x * dl2 - 69.41 * dl * x1 - 384. * d27 * dl2 + 40.239 / x * x12;
    const double cls3f = ( ( 107.0 + 321.05 * x - 54.62 * x2 ) * ( 1 - x ) - 26.717 - 320. * d81 * dl3 - 640. * d81 * dl2
                           + 9.773 * dl + x * dl * ( 363.8 + 68.32 * dl ) ) * x;
    return cls31 + _nf * cls32 + ( fl11sg[_nf-1] - fl11ns[_nf-1] ) * cls3f;
  }

  //_________________________________________________________________________________
  CL3g::CL3g(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double CL3g::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double dl    = log(x);
    const double dl2   = dl * dl;
    const double dl3   = dl * dl2;
    const double dl4   = dl * dl3;
    const double x1    = 1 - x;
    const double dl1   = log(x1);
    const double dl12  = dl1 * dl1;
    const double dl13  = dl1 * dl12;
    const double dl14  = dl1 * dl13;
    const double d27   = 1. / 27;
    const double clg31 = ( 144. * dl14 - 47024. * d27 * dl13 + 6319. * dl12 + 53160. * dl1 ) * x1 + dl * dl1 * ( 72549. + 88238. * dl )
                         + ( 3709. - 33514. * x - 9533. * x2 ) * x1 + 66773. * x * dl2 - 1117. * dl + 45.37 * dl2 - 5360. * d27 * dl3 - 2044.70 / x * x1 - 409.506 * dl / x;
    const double clg32 = ( 288. * d27 * dl13 - 3648. * d27 * dl12 - 592.3 * dl1 + 1511. * x * dl1 ) * x1 + dl * dl1 * ( 311.3 + 14.24 * dl )
                         + ( 577.3 - 729.0 * x ) * x1 + 30.78 * x * dl3 + 366.0 * dl + 3000. * d27 * dl2 + 480. * d27 * dl3 + 88.5037 / x * x1;
    const double clg3f = ( - 0.0105 * dl13 + 1.550 * dl12 + 19.72 * x * dl1 - 66.745 * x + 0.615 * x2 ) * x1 + 20. * d27 * x * dl4
                         + ( 280. / 81. + 2.260 * x) * x * dl3 - ( 15.40 - 2.201 * x ) * x * dl2 - ( 71.66 - 0.121 * x ) * x * dl;
    return clg31 + _nf * ( clg32 + fl11sg[_nf-1] * clg3f );
  }

  //_________________________________________________________________________________
  C33nsp::C33nsp(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C33nsp::Regular(double const& x) const
  {
    const double fl02 = 1;
    const double dl   = log(x);
    const double dl2  = dl * dl;
    const double dl3  = dl * dl2;
    const double dl4  = dl * dl3;
    const double dl5  = dl * dl4;
    const double dl6  = dl * dl5;
    const double x1   = 1 - x;
    const double dl1  = log(x1);
    const double dl12 = dl1 * dl1;
    const double dl13 = dl1 * dl12;
    const double dl14 = dl1 * dl13;
    const double dl15 = dl1 * dl14;
    const double d27  = 1. / 27;
    const double d243 = 1. / 243;

    // Minus piece
    const double CL3nsm =
      - 1853. - 5709. * x + x * x1 * ( 5600. - 1432. * x ) - 536. / 405. * dl5 - 4036. / 81. * dl4 - 496.95 * dl3 - 1488. * dl2 - 293.3 * dl
      - 512. * d27 * dl15 + 8896. * d27 * dl14 - 1396. * dl13 + 3990. * dl12 + 14363. * dl1 - 0.463 * x * dl6 - dl* dl1 * ( 4007. + 1312. * dl )
      + _nf * ( 516.1 - 465.2 * x + x * x1 * ( 635.3 + 310.4 * x ) + 304. / 81. * dl4 + 48512. / 729. * dl3 + 305.32 * dl2 + 366.9 * dl - 1.200 * x * dl4
                - 640. / 81. * dl14 + 32576. / 243. * dl13  - 660.7 * dl12 + 959.1 * dl1 + 31.95 * x1 * dl14 + dl * dl1 * ( 1496. + 270.1 * dl - 1191. * dl1 ) )
      + _nf * _nf * ( 11.32 + 51.94 * x - x * x1 * ( 44.52 + 11.05 * x ) - 368. * d243* dl3 - 2848. / 243. * dl2 - 16.00 * dl
                      - 64. / 81. * dl13 + 992. / 81. * dl12 - 49.65 * dl1 - dl* dl1 * ( 39.99 + 5.103 * dl - 16.30 * dl1 ) + 0.0647 * x * dl4 )
      + fl02 * _nf * ( 48.79 - ( 242.4 - 150.7 * x ) * x1 - 16. / 27. * dl5 + 17.26* dl3 - 113.4 * dl2 - 477.0 * dl + 2.147 * dl12 - 24.57 * dl1
                       + x * dl * ( 218.1 + 82.27 * dl2 ) - dl * dl1 * ( 81.70 + 9.412 * dl1 ) ) * x1;

    // Compute and include difference to get the minus piece
    const double c3q30a = - ( 46.72 * dl12 + 267.26 * dl1 + 719.49 * x ) * x1 - 171.98 * dl + 9.470 * dl3;
    const double c3q30b = ( 3.216 * dl12 + 44.50 * dl1 - 34.588 ) * x1 + 98.719 * dl2 + 2.6208 * dl5;
    const double c3q31a = ( 0.8489 * dl1 + 67.928 * ( 1. + 0.5 * x ) ) * x1 + 97.922 * x * dl - 17.070 * dl2 - 3.132 * dl3;
    const double c3q31b = - ( 0.186 * dl1 + 61.102 * ( 1. + x ) ) * x1 - 122.51 * x * dl + 10.914 * dl2 + 2.748 * dl3;
    return CL3nsm - 0.5 * ( c3q30a + c3q30b + _nf * ( c3q31a + c3q31b ) );
  }
  double C33nsp::Singular(double const& x) const
  {
    const double dl1  = log(1 - x);
    const double dl12 = dl1 * dl1;
    const double dl13 = dl1 * dl12;
    const double dl14 = dl1 * dl13;
    const double dl15 = dl1 * dl14;
    const double d81  = 1. / 81;
    const double c3ns3b =
      + 1536. * d81 * dl15 - 16320. * d81 * dl14 + 5.01099e+2 * dl13 + 1.17154e+3 * dl12 - 7.32845e+3 * dl1 + 4.44276e+3
      + _nf * ( 640. * d81 * dl14 - 6592. * d81 * dl13  + 220.573 * dl12 + 294.906 * dl1 - 729.359 )
      + _nf * _nf * ( 64. * d81 * dl13 - 464. * d81 * dl12 + 7.67505 * dl1 + 1.00830 );
    return c3ns3b / ( 1 - x );
  }
  double C33nsp::Local(double const& x) const
  {
    const double dl1  = log(1 - x);
    const double dl12 = dl1 * dl1;
    const double dl13 = dl1 * dl12;
    const double dl14 = dl1 * dl13;
    const double dl15 = dl1 * dl14;
    const double dl16 = dl1 * dl15;
    const double d3   = 1. / 3;
    const double d81  = 1. / 81;
    return
      + 256. * d81 * dl16 - 3264. * d81 * dl15 + 1.252745e+2 * dl14 + 3.905133e+2 * dl13 - 3.664225e+3 * dl12 + 4.44276e+3  * dl1 - 9195.48 + 22.80
      + _nf * ( 128. * d81 * dl15 - 1648. * d81 * dl14 + 220.573 * d3 * dl13 + 147.453 * dl12 - 729.359 * dl1 + 2575.074 + 0.386 )
      + _nf * _nf * ( 16. * d81 * dl14 - 464. * d81* d3 * dl13 + 7.67505 * 0.5 * dl12 + 1.0083 * dl1 - 103.2521 - 0.0081 );
  }

  //_________________________________________________________________________________
  C33nsm::C33nsm(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C33nsm::Regular(double const& x) const
  {
    const double fl02 = 1;
    const double dl   = log(x);
    const double dl2  = dl * dl;
    const double dl3  = dl * dl2;
    const double dl4  = dl * dl3;
    const double dl5  = dl * dl4;
    const double dl6  = dl * dl5;
    const double x1   = 1 - x;
    const double dl1  = log(x1);
    const double dl12 = dl1 * dl1;
    const double dl13 = dl1 * dl12;
    const double dl14 = dl1 * dl13;
    const double dl15 = dl1 * dl14;
    const double d27  = 1. / 27;
    const double d243 = 1. / 243;
    return
      - 1853. - 5709. * x + x * x1 * ( 5600. - 1432. * x ) - 536. / 405. * dl5 - 4036. / 81. * dl4 - 496.95 * dl3 - 1488. * dl2 - 293.3 * dl
      - 512. * d27 * dl15 + 8896. * d27 * dl14 - 1396. * dl13 + 3990. * dl12 + 14363. * dl1 - 0.463 * x * dl6 - dl* dl1 * ( 4007. + 1312. * dl )
      + _nf * ( 516.1 - 465.2 * x + x * x1 * ( 635.3 + 310.4 * x ) + 304. / 81. * dl4 + 48512. / 729. * dl3 + 305.32 * dl2 + 366.9 * dl - 1.200 * x * dl4
                - 640. / 81. * dl14 + 32576. / 243. * dl13  - 660.7 * dl12 + 959.1 * dl1 + 31.95 * x1 * dl14 + dl * dl1 * ( 1496. + 270.1 * dl - 1191. * dl1 ) )
      + _nf * _nf * ( 11.32 + 51.94 * x - x * x1 * ( 44.52 + 11.05 * x ) - 368. * d243* dl3 - 2848. / 243. * dl2 - 16.00 * dl
                      - 64. / 81. * dl13 + 992. / 81. * dl12 - 49.65 * dl1 - dl* dl1 * ( 39.99 + 5.103 * dl - 16.30 * dl1 ) + 0.0647 * x * dl4 )
      + fl02 * _nf * ( 48.79 - ( 242.4 - 150.7 * x ) * x1 - 16. / 27. * dl5 + 17.26* dl3 - 113.4 * dl2 - 477.0 * dl + 2.147 * dl12 - 24.57 * dl1
                       + x * dl * ( 218.1 + 82.27 * dl2 ) - dl * dl1 * ( 81.70 + 9.412 * dl1 ) ) * x1;
  }
  double C33nsm::Singular(double const& x) const
  {
    const double dl1  = log(1 - x);
    const double dl12 = dl1 * dl1;
    const double dl13 = dl1 * dl12;
    const double dl14 = dl1 * dl13;
    const double dl15 = dl1 * dl14;
    const double d81  = 1. / 81;
    const double c3ns3b =
      + 1536. * d81 * dl15 - 16320. * d81 * dl14 + 5.01099e+2 * dl13 + 1.17154e+3 * dl12 - 7.32845e+3 * dl1 + 4.44276e+3
      + _nf * ( 640. * d81 * dl14 - 6592. * d81 * dl13  + 220.573 * dl12 + 294.906 * dl1 - 729.359 )
      + _nf * _nf * ( 64. * d81 * dl13 - 464. * d81 * dl12 + 7.67505 * dl1 + 1.00830 );
    return c3ns3b / ( 1 - x );
  }
  double C33nsm::Local(double const& x) const
  {
    const double dl1  = log(1 - x);
    const double dl12 = dl1 * dl1;
    const double dl13 = dl1 * dl12;
    const double dl14 = dl1 * dl13;
    const double dl15 = dl1 * dl14;
    const double dl16 = dl1 * dl15;
    const double d3   = 1. / 3;
    const double d81  = 1. / 81;
    return
      + 256. * d81 * dl16 - 3264. * d81 * dl15 + 1.252745e+2 * dl14 + 3.905133e+2 * dl13 - 3.664225e+3 * dl12 + 4.44276e+3  * dl1 - 9195.48 + 22.80
      + _nf * ( 128. * d81 * dl15 - 1648. * d81 * dl14 + 220.573 * d3 * dl13 + 147.453 * dl12 - 729.359 * dl1 + 2575.074 + 0.386 )
      + _nf * _nf * ( 16. * d81 * dl14 - 464. * d81* d3 * dl13 + 7.67505 * 0.5 * dl12 + 1.0083 * dl1 - 103.2521 - 0.0081 );
  }
}
