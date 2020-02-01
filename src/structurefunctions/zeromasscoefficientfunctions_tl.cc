//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/zeromasscoefficientfunctions_tl.h"
#include "apfel/constants.h"
#include "apfel/specialfunctions.h"

namespace apfel
{
  //_________________________________________________________________________________
  C21Tns::C21Tns():
    Expression()
  {
  }
  double C21Tns::Regular(double const& x) const
  {
    return 2 * CF * ( - ( 1 + x ) * log( 1 - x ) + 2 * ( 1 + pow(x, 2) ) * log(x) / ( 1 - x ) + 5. / 2. - 3 * x / 2 );
  }
  double C21Tns::Singular(double const& x) const
  {
    return 2 * CF * ( 2 * log( 1 - x ) - 3 / 2. ) / ( 1 - x );
  }
  double C21Tns::Local(double const& x) const
  {
    return 2 * CF * ( pow(log( 1 - x ), 2) - 3 * log( 1 - x ) / 2 + ( 4 * zeta2 - 9 / 2. ) );
  }

  //_________________________________________________________________________________
  C21Tg::C21Tg():
    Expression()
  {
  }
  double C21Tg::Regular(double const& x) const
  {
    return 4 * CF * ( ( 1 + pow(( 1 - x ), 2) ) * log( pow(x, 2) * ( 1 - x ) ) / x );
  }

  //_________________________________________________________________________________
  C22Tnsp::C22Tnsp(int const& nf):
    Expression(),
    _nf(nf)
  {
    const double CF2 =  CF * CF;
    _a3 = 8.*CF2;
    _a2 = - 22./3.*CA*CF
          - 18.*CF2
          + 4./3.*CF*_nf;
    _a1 = - 8.*zeta2*CA*CF
          + 16.*zeta2*CF2
          + 367./9.*CA*CF
          - 27.*CF2
          - 58./9.*CF*_nf;
    _a0 = + 44./3.*zeta2*CA*CF
          + 40.*zeta3*CA*CF
          - 8.*zeta3*CF2
          - 3155./54.*CA*CF
          + 51./2.*CF2
          + 247./27.*CF*_nf
          - 8./3.*zeta2*CF*_nf;
  }
  double C22Tnsp::Regular(double const& x) const
  {
    // Useful definitions
    const double CF2 =  CF * CF;
    const double x2  = x * x;
    const double x3  = x * x2;
    const double dx  = 1 / x;
    const double dx2 = dx * dx;;
    const double dm  = 1 / ( 1 - x );
    const double dp  = 1 / ( 1 + x );
    const double dl1 = log( 1 - x );

    // Allocate pointers for the harmonic polylogs
    double wx = x;
    int    nw = 3;
    int    n1 = -1;
    int    n2 = 1;
    int    sz = n2 - n1 + 1;
    double *Hr1 = new double[sz];
    double *Hr2 = new double[sz*sz];
    double *Hr3 = new double[sz*sz*sz];
    double *Hr4 = new double[sz*sz*sz*sz];

    // Call polylogs
    hplog_(&wx, &nw, Hr1, Hr2, Hr3, Hr4, &n1, &n2);

    const double Hr1m1 = Hr1[0];
    const double Hr10  = Hr1[1];
    const double Hr11  = Hr1[2];

    //const double Hr2m1m1 = Hr2[0];
    //const double Hr20m1  = Hr2[1];
    //const double Hr21m1  = Hr2[2];
    const double Hr2m10  = Hr2[3];
    const double Hr200   = Hr2[4];
    const double Hr210   = Hr2[5];
    //const double Hr2m11  = Hr2[6];
    const double Hr201   = Hr2[7];
    const double Hr211   = Hr2[8];

    //const double Hr3m1m1m1 = Hr3[0];
    //const double Hr30m1m1  = Hr3[1];
    //const double Hr31m1m1  = Hr3[2];
    //const double Hr3m10m1  = Hr3[3];
    //const double Hr300m1   = Hr3[4];
    //const double Hr310m1   = Hr3[5];
    //const double Hr3m11m1  = Hr3[6];
    //const double Hr301m1   = Hr3[7];
    //const double Hr311m1   = Hr3[8];
    const double Hr3m1m10  = Hr3[9];
    const double Hr30m10   = Hr3[10];
    //const double Hr31m10   = Hr3[11];
    const double Hr3m100   = Hr3[12];
    const double Hr3000    = Hr3[13];
    const double Hr3100    = Hr3[14];
    //const double Hr3m110   = Hr3[15];
    const double Hr3010    = Hr3[16];
    const double Hr3110    = Hr3[17];
    //const double Hr3m1m11  = Hr3[18];
    //const double Hr30m11   = Hr3[19];
    //const double Hr31m11   = Hr3[20];
    const double Hr3m101   = Hr3[21];
    const double Hr3001    = Hr3[22];
    const double Hr3101    = Hr3[23];
    //const double Hr3m111   = Hr3[24];
    const double Hr3011    = Hr3[25];
    const double Hr3111    = Hr3[26];

    const double C2eq2 =
      // Transverse contribution
      + CF*CA * ( 1271./270. + 4829./270.*x - 24./5.*x2
                  + 24./5.*dx - 3155./54.*dm - 36.*zeta3*x - 28.*zeta3*
                  dp + 28.*zeta3*dm - 4.*zeta2*x + 24./5.*zeta2*x3 - 12.*
                  Hr1m1*zeta2 + 20.*Hr1m1*zeta2*x + 32.*Hr1m1*zeta2*dp - 359.
                  /15.*Hr10 - 143./5.*Hr10*x - 24./5.*Hr10*x2
                  - 24./5.*Hr10*dx + 8.*Hr10*dp + 206./3.*Hr10
                  *dm + 12.*Hr10*zeta2 + 20.*Hr10*zeta2*x + 8.*Hr10*zeta2*dp
                  - 32.*Hr10*zeta2*dm - 25./9.*Hr11 + 311./9.*Hr11
                  *x - 367./9.*Hr11*dm + 8.*Hr11*zeta2 - 8.*Hr11*zeta2
                  *dm - 4.*Hr2m10 - 4.*Hr2m10*x + 24./5.*Hr2m10
                  *x3 + 24./5.*Hr2m10*dx2 + 11./3.*Hr200 +
                  23./3.*Hr200*x - 24./5.*Hr200*x3 - 22./3.
                  *Hr200*dm - 22./3.*Hr201 - 22./3.*Hr201*x +
                  44./3.*Hr201*dm + 22./3.*Hr211 + 22./3.*
                  Hr211*x - 44./3.*Hr211*dm - 8.*Hr3m1m10 + 24.
                  *Hr3m1m10*x )
      + CF*CA * ( 32.*Hr3m1m10*dp - 16.*Hr3m100
                  + 8.*Hr3m100*x + 24.*Hr3m100*dp + 8.*Hr3m101
                  - 8.*Hr3m101*x - 16.*Hr3m101*dp + 16.*Hr30m10
                  *x + 8.*Hr30m10*dp - 24.*Hr30m10*dm - 36.
                  *Hr3000*x - 36.*Hr3000*dp + 36.*Hr3000*dm - 4.
                  *Hr3001 + 4.*Hr3001*x + 8.*Hr3001*dp + 4.
                  *Hr3010 + 4.*Hr3010*x - 8.*Hr3010*dm + 4.*
                  Hr3100 + 12.*Hr3100*x - 16.*Hr3100*dm - 4.*
                  Hr3101 - 4.*Hr3101*x + 8.*Hr3101*dm + 4.*
                  Hr3110 + 4.*Hr3110*x - 8.*Hr3110*dm )
      + CF2 * ( 279./10. - 279./10.*x + 48./5.
                *x2 - 48./5.*dx + 51./2.*dm + 56.*zeta3 + 128.*
                zeta3*x + 56.*zeta3*dp - 152.*zeta3*dm - 24.*zeta2 - 48./5.*zeta2*
                x3 + 12.*zeta2*dm + 24.*Hr1m1*zeta2 - 40.*Hr1m1*zeta2*x -
                64.*Hr1m1*zeta2*dp + 376./5.*Hr10 + 166./5.*Hr10
                *x + 48./5.*Hr10*x2 + 48./5.*Hr10*dx - 16.*
                Hr10*dp - 106.*Hr10*dm - 4.*Hr10*zeta2 - 20.*Hr10*
                zeta2*x - 16.*Hr10*zeta2*dp + 40.*Hr10*zeta2*dm + 13.*Hr11
                - 51.*Hr11*x + 27.*Hr11*dm - 8.*Hr11*zeta2 + 8.*
                Hr11*zeta2*x + 8.*Hr2m10 + 8.*Hr2m10*x - 48./5.*
                Hr2m10*x3 - 48./5.*Hr2m10*dx2 - 66.*Hr200
                - 30.*Hr200*x + 48./5.*Hr200*x3 + 66.*Hr200
                *dm + 12.*Hr201 - 4.*Hr201*x + 12.*Hr201*dm
                - 28.*Hr210 + 28.*Hr210*x + 24.*Hr210*dm + 16.
                *Hr211 + 8.*Hr211*x - 36.*Hr211*dm + 16.*
                Hr3m1m10 )
      + CF2 * (  - 48.*Hr3m1m10*x - 64.*Hr3m1m10
                 *dp + 32.*Hr3m100 - 16.*Hr3m100*x - 48.*
                 Hr3m100*dp - 16.*Hr3m101 + 16.*Hr3m101*x + 32.
                 *Hr3m101*dp - 32.*Hr30m10*x - 16.*Hr30m10*
                 dp + 48.*Hr30m10*dm + 66.*Hr3000 + 138.*Hr3000
                 *x + 72.*Hr3000*dp - 160.*Hr3000*dm - 8.*
                 Hr3001 - 24.*Hr3001*x - 16.*Hr3001*dp + 8.*
                 Hr3001*dm + 36.*Hr3010 + 36.*Hr3010*x - 72.
                 *Hr3010*dm - 16.*Hr3011 - 16.*Hr3011*x + 40.
                 *Hr3011*dm - 12.*Hr3100 - 28.*Hr3100*x + 40.
                 *Hr3100*dm - 16.*Hr3101 - 16.*Hr3101*x + 32.
                 *Hr3101*dm - 24.*Hr3110 - 24.*Hr3110*x + 48.
                 *Hr3110*dm + 24.*Hr3111 + 24.*Hr3111*x - 48.
                 *Hr3111*dm )
      + _nf*CF * ( - 59./27. - 17./27.*x + 247./
                   27.*dm + 10./3.*Hr10 + 6.*Hr10*x - 32./3.*
                   Hr10*dm - 14./9.*Hr11 - 26./9.*Hr11*x + 58./9.
                   *Hr11*dm - 2./3.*Hr200 - 2./3.*Hr200*x + 4.
                   /3.*Hr200*dm + 4./3.*Hr201 + 4./3.*Hr201*
                   x - 8./3.*Hr201*dm - 4./3.*Hr211 - 4./3.*
                   Hr211*x + 8./3.*Hr211*dm )
      // Longitudinal contribution
      + CF * CA * ( 1729./45. - 98./15.*x - 16./5.*x2 -
                    24./5.*dx + 16.*zeta2*x + 16./5.*zeta2*x3 - 8.*Hr1m1*zeta2
                    - 146./15.*Hr10 + 8./5.*Hr10*x - 16./5.*
                    Hr10*x2 + 24./5.*Hr10*dx + 46./3.*Hr11 - 8.
                    *Hr11*zeta2 + 8.*Hr2m10 + 16.*Hr2m10*x + 16./5.*
                    Hr2m10*x3 - 24./5.*Hr2m10*dx2 - 16.*Hr200*
                    x - 16./5.*Hr200*x3 - 16.*Hr3m1m10 + 8.*Hr3m100
                    + 16.*Hr30m10 + 8.*Hr3100 )
      + CF * CF * ( - 147./5. - 18./5.*x + 32./5.
                    *x2 + 48./5.*dx + 4.*zeta2 - 32.*zeta2*x - 32./5.*
                    zeta2*x3 + 16.*Hr1m1*zeta2 + 34./5.*Hr10 + 24./5.*
                    Hr10*x + 32./5.*Hr10*x2 - 48./5.*Hr10*dx - 14.
                    *Hr11 - 4.*Hr11*x + 16.*Hr11*zeta2 - 16.*Hr2m10
                    - 32.*Hr2m10*x - 32./5.*Hr2m10*x3 + 48./5.
                    *Hr2m10*dx2 - 12.*Hr200 + 32.*Hr200*x + 32./
                    5.*Hr200*x3 - 4.*Hr201 - 16.*Hr210 + 8.*
                    Hr211 + 32.*Hr3m1m10 - 16.*Hr3m100 - 32.*
                    Hr30m10 - 16.*Hr3100 )
      + _nf * CF * ( - 50./9. + 4./3.*x + 4./3.
                     *Hr10 - 4./3.*Hr11 );

    // Deallocate pointers for the harmonic polylogs
    delete[] Hr4;
    delete[] Hr3;
    delete[] Hr2;
    delete[] Hr1;

    // Singular term to be subtracted
    const double C2eq2L = dm * ( pow(dl1, 3) * _a3 + pow(dl1, 2) * _a2 + dl1 * _a1 + _a0 );

    // Return regular
    return C2eq2 - C2eq2L;
  }
  double C22Tnsp::Singular(double const& x) const
  {
    const double dl1 = log( 1 - x );
    return ( pow(dl1, 3) * _a3 + pow(dl1, 2) * _a2 + dl1 * _a1 + _a0 ) / ( 1 - x );;
  }
  double C22Tnsp::Local(double const& x) const
  {
    const double c2delt =
      + CA*CF * ( - 5465./72. + 140./3.*zeta3 + 215./3.*zeta2
                  - 49./5.*pow(zeta2, 2) )
      + CF*CF * ( 331./8. - 78.*zeta3 - 39.*zeta2 + 30.*pow(zeta2, 2) )
      + CF*_nf * ( 457./36. + 4./3.*zeta3 - 38./3.*zeta2 );
    const double dl1 = log( 1 - x );
    return pow(dl1, 4) * _a3 / 4 + pow(dl1, 3) * _a2 / 3 + pow(dl1, 2) * _a1 / 2 + dl1 * _a0 + c2delt;
  }

  //_________________________________________________________________________________
  C22Tps::C22Tps():
    Expression()
  {
  }
  double C22Tps::Regular(double const& x) const
  {
    // Useful definitions
    const double x2 = x * x;

    // Allocate pointers for the harmonic polylogs
    double wx = x;
    int    nw = 3;
    int    n1 = -1;
    int    n2 = 1;
    int    sz = n2 - n1 + 1;
    double *Hr1 = new double[sz];
    double *Hr2 = new double[sz*sz];
    double *Hr3 = new double[sz*sz*sz];
    double *Hr4 = new double[sz*sz*sz*sz];

    // Call polylogs
    hplog_(&wx, &nw, Hr1, Hr2, Hr3, Hr4, &n1, &n2);

    const double Hr10 = Hr1[1];
    const double Hr11 = Hr1[2];

    const double Hr2m10 = Hr2[3];
    const double Hr200  = Hr2[4];
    const double Hr201  = Hr2[7];
    const double Hr211  = Hr2[8];

    const double Hr3000 = Hr3[13];
    const double Hr3001 = Hr3[22];
    const double Hr3011 = Hr3[25];

    const double cTeqps2 =
      // Transverse contribution
      + CF * ( - 118./3. + 70./3.*x + 512./27.*x2
               - 80./27./x + 16*zeta3 + 16*zeta3*x - 8*zeta2 - 32*zeta2*x - 200./
               3.*Hr10 - 104./3.*Hr10*x - 128./9.*Hr10*x2
               - 16./3.*Hr10/x + 16*Hr10*zeta2 + 16*Hr10*zeta2*x + 92.
               /3.*Hr11 - 68./3.*Hr11*x - 32./3.*Hr11*x2
               + 8./3.*Hr11/x - 16*Hr2m10 - 16*Hr2m10*x - 16.
               /3.*Hr2m10*x2 - 16./3.*Hr2m10/x - 14*Hr200
               - 14*Hr200*x + 16./3.*Hr200*x2 + 64./3.*Hr200/x
               + 4*Hr201 + 20*Hr201*x + 16./3.*Hr201*
               x2 - 32./3.*Hr201/x + 4*Hr211 - 4*Hr211*x -
               16./3.*Hr211*x2 + 16./3.*Hr211/x + 44*Hr3000
               + 44*Hr3000*x - 24*Hr3001 - 24*Hr3001*x + 8*
               Hr3011 + 8*Hr3011*x )
      // Longitudinal contribution
      + CF * ( - 56./3. + 104./3.*x - 8*x2 - 8/x + 8*
               zeta2 - 16*Hr10 - 16*Hr10*x + 8./3.*Hr10*x2 + 32./
               3.*Hr10/x + 8*Hr11*x - 8./3.*Hr11*x2 - 16./3.
               *Hr11/x + 24*Hr200 - 8*Hr201 );

    // Deallocate pointers for the harmonic polylogs
    delete[] Hr4;
    delete[] Hr3;
    delete[] Hr2;
    delete[] Hr1;
    return cTeqps2;
  }

  //_________________________________________________________________________________
  C22Tg::C22Tg():
    Expression()
  {
  }
  double C22Tg::Regular(double const& x) const
  {
    // Useful definitions
    const double x2  = x * x;
    const double x3  = x * x2;
    const double dx  = 1 / x;
    const double dx2 = dx * dx;

    // Allocate pointers for the harmonic polylogs
    double wx = x;
    int    nw = 3;
    int    n1 = -1;
    int    n2 = 1;
    int    sz = n2 - n1 + 1;
    double *Hr1 = new double[sz];
    double *Hr2 = new double[sz*sz];
    double *Hr3 = new double[sz*sz*sz];
    double *Hr4 = new double[sz*sz*sz*sz];

    // Call polylogs
    hplog_(&wx, &nw, Hr1, Hr2, Hr3, Hr4, &n1, &n2);

    const double Hr1m1 = Hr1[0];
    const double Hr10  = Hr1[1];
    const double Hr11  = Hr1[2];

    const double Hr2m10 = Hr2[3];
    const double Hr200  = Hr2[4];
    const double Hr210  = Hr2[5];
    const double Hr201  = Hr2[7];
    const double Hr211  = Hr2[8];

    const double Hr3m1m10 = Hr3[9];
    const double Hr30m10  = Hr3[10];
    const double Hr3m100  = Hr3[12];
    const double Hr3000   = Hr3[13];
    const double Hr3100   = Hr3[14];
    const double Hr3010   = Hr3[16];
    const double Hr3110   = Hr3[17];
    const double Hr3m101  = Hr3[21];
    const double Hr3001   = Hr3[22];
    const double Hr3101   = Hr3[23];
    const double Hr3011   = Hr3[25];
    const double Hr3111   = Hr3[26];

    const double cTeg2 =
      // Transverse contribution
      + CF * CA * ( - 36 - 106*x - 928./27.*x2 + 4438./27.*
                    dx + 64*zeta3 - 136*zeta3*x - 240*zeta3*dx + 32*zeta2 + 64*zeta2*x - 56*zeta2*
                    dx + 16*Hr1m1*zeta2 + 8*Hr1m1*zeta2*x + 32*Hr1m1*zeta2*dx + 772.
                    /3.*Hr10 + 172./3.*Hr10*x + 256./9.*Hr10*
                    x2 + 496./3.*Hr10*dx - 128*Hr10*zeta2 - 16*Hr10*zeta2*x
                    + 32*Hr10*zeta2*dx + 236./3.*Hr11 + 4./3.*Hr11*x
                    + 32./3.*Hr11*x2 - 356./3.*Hr11*dx - 48*Hr11
                    *zeta2 + 24*Hr11*zeta2*x + 32*Hr11*zeta2*dx + 80*Hr2m10 + 56*
                    Hr2m10*x + 32./3.*Hr2m10*x2 + 80./3.*Hr2m10*dx
                    - 32*Hr200 + 4*Hr200*x - 32./3.*Hr200*x2
                    - 464./3.*Hr200*dx - 96*Hr201 - 16*Hr201*x - 32.
                    /3.*Hr201*x2 + 496./3.*Hr201*dx - 64*Hr210
                    + 8*Hr210*x + 64*Hr210*dx + 96*Hr211 - 16*Hr211*
                    x + 32./3.*Hr211*x2 - 344./3.*Hr211*dx - 32*
                    Hr3m1m10 - 16*Hr3m1m10*x + 96*Hr3m100 + 48*Hr3m100*x
                    + 80*Hr3m100*dx - 32*Hr3m101 - 16*Hr3m101
                    *x - 32*Hr3m101*dx + 64*Hr30m10 + 32*Hr30m10*x +
                    64*Hr30m10*dx - 176*Hr3000 - 248*Hr3000*x - 320*
                    Hr3000*dx + 128*Hr3001 + 96*Hr3001*x + 96*Hr3001*dx
                    + 96*Hr3010 - 48*Hr3010*x - 96*Hr3010*dx
                    - 48*Hr3011 - 24*Hr3011*x - 16*Hr3011*dx + 64*
                    Hr3100 - 32*Hr3100*x - 48*Hr3100*dx - 32*Hr3101
                    + 16*Hr3101*x + 32*Hr3101*dx - 64*Hr3110 + 32*
                    Hr3110*x + 64*Hr3110*dx + 16*Hr3111 - 8*Hr3111*x
                    - 16*Hr3111*dx )
      + CF * CF * ( - 604./5. + 154./5.*x - 16./
                    5.*x2 + 316./5.*dx + 32*zeta3 - 80*zeta3*x - 64*zeta3*dx - 32*
                    zeta2 - 72*zeta2*x + 16./5.*zeta2*x3 + 64*Hr1m1*zeta2 + 32*Hr1m1*zeta2*x
                    + 32*Hr1m1*zeta2*dx + 418./5.*Hr10 - 262./5.*
                    Hr10*x - 16./5.*Hr10*x2 + 144./5.*Hr10*dx +
                    32*Hr10*zeta2 - 16*Hr10*zeta2*x + 24*Hr11*x - 8*Hr11*dx +
                    80*Hr11*zeta2 - 40*Hr11*zeta2*x - 48*Hr11*zeta2*dx - 64*Hr2m10
                    - 96*Hr2m10*x + 16./5.*Hr2m10*x3 - 64./5.*
                    Hr2m10*dx2 - 64*Hr200 + 166*Hr200*x - 16./5.*
                    Hr200*x3 - 80*Hr201 + 12*Hr201*x + 96*Hr201*dx
                    - 16*Hr210*x + 112*Hr211 - 28*Hr211*x - 96*Hr211
                    *dx + 128*Hr3m1m10 + 64*Hr3m1m10*x + 64*Hr3m1m10*
                    dx - 64*Hr3m100 - 32*Hr3m100*x - 32*Hr3m100*dx -
                    128*Hr30m10 + 88*Hr3000 - 44*Hr3000*x + 16*Hr3001
                    - 8*Hr3001*x - 64*Hr3001*dx + 64*Hr3010 - 32
                    *Hr3010*x - 64*Hr3010*dx - 80*Hr3011 + 40*Hr3011*x
                    + 96*Hr3011*dx - 128*Hr3100 + 64*Hr3100*x
                    + 96*Hr3100*dx - 48*Hr3101 + 24*Hr3101*x + 48*
                    Hr3101*dx - 16*Hr3110 + 8*Hr3110*x + 16*Hr3110*dx
                    + 80*Hr3111 - 40*Hr3111*x - 80*Hr3111*dx )
      // Longitudinal contribution
      + CF * CA * ( - 320./3. - 160./3.*x + 32./3.*x2 +
                    448./3.*dx - 64*zeta2 + 32*zeta2*dx + 112*Hr10 + 32*Hr10*x
                    - 16./3.*Hr10*x2 - 352./3.*Hr10*dx - 144*Hr11
                    - 16*Hr11*x + 16./3.*Hr11*x2 + 464./3.*Hr11
                    *dx + 32*Hr2m10 + 32*Hr2m10*dx - 96*Hr200 - 128*Hr200*dx
                    + 64*Hr201 + 64*Hr210 - 64*Hr210*dx - 32*
                    Hr211 + 32*Hr211*dx )
      + CF * CF * ( 24./5. + 248./15.*x - 32./15.
                    *x2 - 96./5.*dx + 16*zeta2 + 32./15.*zeta2*x3 - 8./5.
                    *Hr10 - 224./15.*Hr10*x - 32./15.*Hr10*x2
                    + 96./5.*Hr10*dx + 24*Hr11 + 8*Hr11*x - 32*Hr11*
                    dx - 32./3.*Hr2m10 + 32./15.*Hr2m10*x3 + 64.
                    /5.*Hr2m10*dx2 + 48*Hr200 - 32./15.*Hr200*
                    x3 - 16*Hr201 );

    // Deallocate pointers for the harmonic polylogs
    delete[] Hr4;
    delete[] Hr3;
    delete[] Hr2;
    delete[] Hr1;
    return cTeg2;
  }

  //_________________________________________________________________________________
  CL1Tns::CL1Tns():
    Expression()
  {
  }
  double CL1Tns::Regular(double const&) const
  {
    return 2 * CF;
  }

  //_________________________________________________________________________________
  CL1Tg::CL1Tg():
    Expression()
  {
  }
  double CL1Tg::Regular(double const& x) const
  {
    return 8 * CF * ( 1 - x ) / x;
  }

  //_________________________________________________________________________________
  CL2Tnsp::CL2Tnsp(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double CL2Tnsp::Regular(double const& x) const
  {
    // Useful definitions
    const double x2  = x * x;
    const double x3  = x * x2;
    const double dx  = 1 / x;
    const double dx2 = dx * dx;

    // Allocate pointers for the harmonic polylogs
    double wx = x;
    int    nw = 3;
    int    n1 = -1;
    int    n2 = 1;
    int    sz = n2 - n1 + 1;
    double *Hr1 = new double[sz];
    double *Hr2 = new double[sz*sz];
    double *Hr3 = new double[sz*sz*sz];
    double *Hr4 = new double[sz*sz*sz*sz];

    // Call polylogs
    hplog_(&wx, &nw, Hr1, Hr2, Hr3, Hr4, &n1, &n2);

    const double Hr1m1 = Hr1[0];
    const double Hr10  = Hr1[1];
    const double Hr11  = Hr1[2];

    const double Hr2m10 = Hr2[1];
    const double Hr200  = Hr2[4];
    const double Hr201  = Hr2[5];
    const double Hr210  = Hr2[7];
    const double Hr211  = Hr2[8];

    const double Hr3m1m10 = Hr3[9];
    const double Hr30m10  = Hr3[10];
    const double Hr3m100  = Hr3[12];
    const double Hr3100   = Hr3[14];

    const double CLeq2 =
      + CF * CA * ( 1729./45. - 98./15.*x - 16./5.*x2 -
                    24./5.*dx + 16.*zeta2*x + 16./5.*zeta2*x3 - 8.*Hr1m1*zeta2
                    - 146./15.*Hr10 + 8./5.*Hr10*x - 16./5.*
                    Hr10*x2 + 24./5.*Hr10*dx + 46./3.*Hr11 - 8.
                    *Hr11*zeta2 + 8.*Hr2m10 + 16.*Hr2m10*x + 16./5.*
                    Hr2m10*x3 - 24./5.*Hr2m10*dx2 - 16.*Hr200*
                    x - 16./5.*Hr200*x3 - 16.*Hr3m1m10 + 8.*Hr3m100
                    + 16.*Hr30m10 + 8.*Hr3100 )
      + CF * CF * ( - 147./5. - 18./5.*x + 32./5.
                    *x2 + 48./5.*dx + 4.*zeta2 - 32.*zeta2*x - 32./5.*
                    zeta2*x3 + 16.*Hr1m1*zeta2 + 34./5.*Hr10 + 24./5.*
                    Hr10*x + 32./5.*Hr10*x2 - 48./5.*Hr10*dx - 14.
                    *Hr11 - 4.*Hr11*x + 16.*Hr11*zeta2 - 16.*Hr2m10
                    - 32.*Hr2m10*x - 32./5.*Hr2m10*x3 + 48./5.
                    *Hr2m10*dx2 - 12.*Hr200 + 32.*Hr200*x + 32./
                    5.*Hr200*x3 - 4.*Hr201 - 16.*Hr210 + 8.*
                    Hr211 + 32.*Hr3m1m10 - 16.*Hr3m100 - 32.*
                    Hr30m10 - 16.*Hr3100 )
      + _nf * CF * ( - 50./9. + 4./3.*x + 4./3.
                     *Hr10 - 4./3.*Hr11 );

    // Deallocate pointers for the harmonic polylogs
    delete[] Hr4;
    delete[] Hr3;
    delete[] Hr2;
    delete[] Hr1;
    return CLeq2;
  }

  //_________________________________________________________________________________
  CL2Tps::CL2Tps():
    Expression()
  {
  }
  double CL2Tps::Regular(double const& x) const
  {
    // Useful definitions
    const double x2 = x * x;

    // Allocate pointers for the harmonic polylogs
    double wx = x;
    int    nw = 3;
    int    n1 = -1;
    int    n2 = 1;
    int    sz = n2 - n1 + 1;
    double *Hr1 = new double[sz];
    double *Hr2 = new double[sz*sz];
    double *Hr3 = new double[sz*sz*sz];
    double *Hr4 = new double[sz*sz*sz*sz];

    // Call polylogs
    hplog_(&wx, &nw, Hr1, Hr2, Hr3, Hr4, &n1, &n2);

    const double Hr10  = Hr1[1];
    const double Hr11  = Hr1[2];

    const double Hr200  = Hr2[4];
    const double Hr201  = Hr2[7];

    const double cLeqps2 =
      + CF * ( - 56./3. + 104./3.*x - 8*x2 - 8/x + 8*
               zeta2 - 16*Hr10 - 16*Hr10*x + 8./3.*Hr10*x2 + 32./
               3.*Hr10/x + 8*Hr11*x - 8./3.*Hr11*x2 - 16./3.
               *Hr11/x + 24*Hr200 - 8*Hr201 );

    // Deallocate pointers for the harmonic polylogs
    delete[] Hr4;
    delete[] Hr3;
    delete[] Hr2;
    delete[] Hr1;
    return cLeqps2;
  }

  //_________________________________________________________________________________
  CL2Tg::CL2Tg():
    Expression()
  {
  }
  double CL2Tg::Regular(double const& x) const
  {
    // Useful definitions
    const double x2  = x * x;
    const double x3  = x * x2;
    const double dx  = 1 / x;
    const double dx2 = dx * dx;

    // Allocate pointers for the harmonic polylogs
    double wx = x;
    int    nw = 3;
    int    n1 = -1;
    int    n2 = 1;
    int    sz = n2 - n1 + 1;
    double *Hr1 = new double[sz];
    double *Hr2 = new double[sz*sz];
    double *Hr3 = new double[sz*sz*sz];
    double *Hr4 = new double[sz*sz*sz*sz];

    // Call polylogs
    hplog_(&wx, &nw, Hr1, Hr2, Hr3, Hr4, &n1, &n2);

    const double Hr10  = Hr1[1];
    const double Hr11  = Hr1[2];

    const double Hr2m10 = Hr2[3];
    const double Hr200  = Hr2[4];
    const double Hr210  = Hr2[5];
    const double Hr201  = Hr2[7];
    const double Hr211  = Hr2[8];

    const double cLeg2 =
      + CF * CA * ( - 320./3. - 160./3.*x + 32./3.*x2 +
                    448./3.*dx - 64*zeta2 + 32*zeta2*dx + 112*Hr10 + 32*Hr10*x
                    - 16./3.*Hr10*x2 - 352./3.*Hr10*dx - 144*Hr11
                    - 16*Hr11*x + 16./3.*Hr11*x2 + 464./3.*Hr11
                    *dx + 32*Hr2m10 + 32*Hr2m10*dx - 96*Hr200 - 128*Hr200*dx
                    + 64*Hr201 + 64*Hr210 - 64*Hr210*dx - 32*
                    Hr211 + 32*Hr211*dx )
      + CF * CF * ( 24./5. + 248./15.*x - 32./15.
                    *x2 - 96./5.*dx + 16*zeta2 + 32./15.*zeta2*x3 - 8./5.
                    *Hr10 - 224./15.*Hr10*x - 32./15.*Hr10*x2
                    + 96./5.*Hr10*dx + 24*Hr11 + 8*Hr11*x - 32*Hr11*
                    dx - 32./3.*Hr2m10 + 32./15.*Hr2m10*x3 + 64.
                    /5.*Hr2m10*dx2 + 48*Hr200 - 32./15.*Hr200*
                    x3 - 16*Hr201 );

    // Deallocate pointers for the harmonic polylogs
    delete[] Hr4;
    delete[] Hr3;
    delete[] Hr2;
    delete[] Hr1;
    return cLeg2;
  }

  //_________________________________________________________________________________
  C31Tns::C31Tns():
    Expression()
  {
  }
  double C31Tns::Regular(double const& x) const
  {
    return 2 * CF * ( - ( 1 + x ) * log( 1 - x ) - 2 * ( 1 + pow(x, 2) ) * log(x) / ( 1 - x ) + 1. / 2. - x / 2 );
  }
  double C31Tns::Singular(double const& x) const
  {
    return 2 * CF * ( 2 * log( 1 - x ) - 3 / 2. ) / ( 1 - x );
  }
  double C31Tns::Local(double const& x) const
  {
    return 2 * CF * ( pow(log( 1 - x ), 2) - 3 * log( 1 - x ) / 2 + ( 4 * zeta2 - 9 / 2. ) );
  }

  //_________________________________________________________________________________
  C32Tnsp::C32Tnsp(int const& nf):
    Expression(),
    _nf(nf)
  {
    const double CF2 = CF * CF;
    _a3 = + 8.*CF2;
    _a2 = - 22./3.*CA*CF
          - 18.*CF2
          + 4./3.*CF*_nf;
    _a1 = - 8.*zeta2*CA*CF
          + 16.*zeta2*CF2
          + 367./9.*CA*CF
          - 27.*CF2
          - 58./9.*CF*_nf;
    _a0 = + 44./3.*zeta2*CA*CF
          + 40.*zeta3*CA*CF
          - 8.*zeta3*CF2
          - 3155./54.*CA*CF
          + 51./2.*CF2
          + 247./27.*CF*_nf
          - 8./3.*zeta2*CF*_nf;
  }
  double C32Tnsp::Regular(double const& x) const
  {
    // Useful definitions
    const double CF2 =  CF * CF;
    const double x2  = x * x;
    const double dx  = 1 / x;
    const double dm  = 1 / ( 1 - x );
    const double dp  = 1 / ( 1 + x );
    const double dl1 = log( 1 - x );

    // Allocate pointers for the harmonic polylogs
    double wx = x;
    int    nw = 3;
    int    n1 = -1;
    int    n2 = 1;
    int    sz = n2 - n1 + 1;
    double *Hr1 = new double[sz];
    double *Hr2 = new double[sz*sz];
    double *Hr3 = new double[sz*sz*sz];
    double *Hr4 = new double[sz*sz*sz*sz];

    // Call polylogs
    hplog_(&wx, &nw, Hr1, Hr2, Hr3, Hr4, &n1, &n2);

    const double Hr1m1 = Hr1[0];
    const double Hr10  = Hr1[1];
    const double Hr11  = Hr1[2];

    const double Hr2m10 = Hr2[3];
    const double Hr200  = Hr2[4];
    const double Hr210  = Hr2[5];
    const double Hr201  = Hr2[7];
    const double Hr211  = Hr2[8];

    const double Hr3m1m10 = Hr3[9];
    const double Hr30m10  = Hr3[10];
    const double Hr3m100  = Hr3[12];
    const double Hr3000   = Hr3[13];
    const double Hr3100   = Hr3[14];
    const double Hr3010   = Hr3[16];
    const double Hr3110   = Hr3[17];
    const double Hr3m101  = Hr3[21];
    const double Hr3001   = Hr3[22];
    const double Hr3101   = Hr3[23];
    const double Hr3011   = Hr3[25];
    const double Hr3111   = Hr3[26];

    const double C3eq2 =
      + CF * CA * ( 325./54. + 895./54.*x - 3155./54.*dm -
                    36.*zeta3 + 28.*zeta3*dp + 28.*zeta3*dm + 12.*zeta2 + 8.*zeta2*x
                    + 8.*zeta2*x2 + 20.*Hr1m1*zeta2 - 12.*Hr1m1*zeta2*x - 32.
                    *Hr1m1*zeta2*dp + 27.*Hr10 - 193./3.*Hr10*x - 8.*
                    Hr10*dp + 206./3.*Hr10*dm + 20.*Hr10*zeta2 + 12.*
                    Hr10*zeta2*x - 8.*Hr10*zeta2*dp - 32.*Hr10*zeta2*dm - 19./
                    9.*Hr11 + 305./9.*Hr11*x - 367./9.*Hr11*dm +
                    8.*Hr11*zeta2*x - 8.*Hr11*zeta2*dm + 4.*Hr2m10 + 4.*
                    Hr2m10*x + 8.*Hr2m10*x2 + 8.*Hr2m10*dx + 59./
                    3.*Hr200 + 71./3.*Hr200*x - 8.*Hr200*x2 -
                    22./3.*Hr200*dm - 46./3.*Hr201 - 46./3.*
                    Hr201*x + 44./3.*Hr201*dm + 22./3.*Hr211 +
                    22./3.*Hr211*x - 44./3.*Hr211*dm + 24.*Hr3m1m10
                    - 8.*Hr3m1m10*x - 32.*Hr3m1m10*dp + 8.*
                    Hr3m100 - 16.*Hr3m100*x - 24.*Hr3m100*dp - 8.
                    *Hr3m101 )
      + CF * CA * ( 8.*Hr3m101*x + 16.*Hr3m101*
                    dp + 16.*Hr30m10 - 8.*Hr30m10*dp - 24.*Hr30m10*dm
                    - 36.*Hr3000 + 36.*Hr3000*dp + 36.*Hr3000*dm
                    + 4.*Hr3001 - 4.*Hr3001*x - 8.*Hr3001*dp
                    + 4.*Hr3010 + 4.*Hr3010*x - 8.*Hr3010*dm
                    + 12.*Hr3100 + 4.*Hr3100*x - 16.*Hr3100*dm
                    - 4.*Hr3101 - 4.*Hr3101*x + 8.*Hr3101*
                    dm + 4.*Hr3110 + 4.*Hr3110*x - 8.*Hr3110*dm )
      + CF2 * (  - 19./2. + 19./2.*x + 51./2.
                 *dm + 128.*zeta3 + 56.*zeta3*x - 56.*zeta3*dp - 152.*zeta3*dm -
                 52.*zeta2 - 20.*zeta2*x - 16.*zeta2*x2 + 12.*zeta2*dm - 40.*
                 Hr1m1*zeta2 + 24.*Hr1m1*zeta2*x + 64.*Hr1m1*zeta2*dp - 2.*
                 Hr10 + 92.*Hr10*x + 16.*Hr10*dp - 106.*Hr10*dm
                 - 20.*Hr10*zeta2 - 4.*Hr10*zeta2*x + 16.*Hr10*zeta2*dp +
                 40.*Hr10*zeta2*dm - 5.*Hr11 - 33.*Hr11*x + 27.*
                 Hr11*dm + 8.*Hr11*zeta2 - 8.*Hr11*zeta2*x - 8.*Hr2m10
                 - 8.*Hr2m10*x - 16.*Hr2m10*x2 - 16.*Hr2m10
                 *dx - 86.*Hr200 - 74.*Hr200*x + 16.*Hr200*x2
                 + 66.*Hr200*dm + 32.*Hr201 + 8.*Hr201*x + 12.
                 *Hr201*dm - 12.*Hr210 + 12.*Hr210*x + 24.*
                 Hr210*dm + 8.*Hr211 + 16.*Hr211*x - 36.*Hr211*dm
                 - 48.*Hr3m1m10 + 16.*Hr3m1m10*x + 64.*
                 Hr3m1m10*dp - 16.*Hr3m100 + 32.*Hr3m100*x +
                 48.*Hr3m100*dp )
      + CF2 * ( 16.*Hr3m101 - 16.*Hr3m101*x
                - 32.*Hr3m101*dp - 32.*Hr30m10 + 16.*Hr30m10*dp
                + 48.*Hr30m10*dm + 138.*Hr3000 + 66.*Hr3000*x
                - 72.*Hr3000*dp - 160.*Hr3000*dm - 24.
                *Hr3001 - 8.*Hr3001*x + 16.*Hr3001*dp + 8.*
                Hr3001*dm + 36.*Hr3010 + 36.*Hr3010*x - 72.
                *Hr3010*dm - 16.*Hr3011 - 16.*Hr3011*x + 40.
                *Hr3011*dm - 28.*Hr3100 - 12.*Hr3100*x + 40.
                *Hr3100*dm - 16.*Hr3101 - 16.*Hr3101*x + 32.
                *Hr3101*dm - 24.*Hr3110 - 24.*Hr3110*x + 48.
                *Hr3110*dm + 24.*Hr3111 + 24.*Hr3111*x - 48.
                *Hr3111*dm )
      + _nf * CF * ( 55./27. - 131./27.*x + 247./
                     27.*dm + 2.*Hr10 + 22./3.*Hr10*x - 32./3.*
                     Hr10*dm - 2./9.*Hr11 - 38./9.*Hr11*x + 58./9.
                     *Hr11*dm - 2./3.*Hr200 - 2./3.*Hr200*x + 4.
                     /3.*Hr200*dm + 4./3.*Hr201 + 4./3.*Hr201*x
                     - 8./3.*Hr201*dm - 4./3.*Hr211 - 4./3.*
                     Hr211*x + 8./3.*Hr211*dm );

    // Deallocate pointers for the harmonic polylogs
    delete[] Hr4;
    delete[] Hr3;
    delete[] Hr2;
    delete[] Hr1;

    // Singular term term to be subtracted
    const double C3eq2L = dm * ( pow(dl1, 3) * _a3 + pow(dl1, 2) * _a2 + dl1 * _a1 + _a0 );

    // Return regular
    return C3eq2 - C3eq2L;
  }
  double C32Tnsp::Singular(double const& x) const
  {
    const double dl1 = log( 1 - x );
    return ( pow(dl1, 3) * _a3 + pow(dl1, 2) * _a2 + dl1 * _a1 + _a0 ) / ( 1 - x );
  }
  double C32Tnsp::Local(double const& x) const
  {
    const double c3delt =
      + CA * CF * ( - 5465./72. + 140./3.*zeta3 + 215./3.*zeta2
                    - 49./5.*pow(zeta2, 2) )
      + CF * CF * ( 331./8. - 78.*zeta3 - 39.*zeta2 + 30.*pow(zeta2, 2) )
      + CF * _nf * ( 457./36. + 4./3.*zeta3 - 38./3.*zeta2 );

    const double dl1 = log( 1 - x );
    return pow(dl1, 4) * _a3/4. + pow(dl1, 3) * _a2/3. + pow(dl1, 2) * _a1/2. + dl1 * _a0 + c3delt;
  }
}
