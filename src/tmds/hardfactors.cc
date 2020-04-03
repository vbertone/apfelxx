//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/hardfactors.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________
  double H1DY()
  {
    return 2 * CF * ( - 8 + 7 * Pi2 / 6 );
  }

  //_________________________________________________________________________
  double H2DY(int const& nf)
  {
    return 2 * CF * ( CF * ( 511. / 8. - 83 * Pi2 / 6 - 30 * zeta3 + 67 * Pi2 * Pi2 / 60 ) +
                      CA * ( - 51157. / 648. + 1061 * Pi2 / 108 + 313 * zeta3 / 9 - 4 * Pi2 * Pi2 / 45 ) +
                      TR * nf * ( 4085. / 162. - 91 * Pi2 / 27 + 4 * zeta3 / 9 ) );
  }

  //_________________________________________________________________________
  double H3DY(int const& nf)
  {
    const int nf2       = nf * nf;
    const double CF2    = CF * CF;
    const double CF3    = CF * CF2;
    const double CA2    = CA * CA;
    const double Pi4    = Pi2 * Pi2;
    const double Pi6    = Pi2 * Pi4;
    const double zeta32 = zeta3 * zeta3;
    return CF*nf2*(-29.100899253162627 + (1612*Pi2)/243. + (86*Pi4)/1215. - (832*zeta3)/243.)
           + CA*CF2*(2544.0771604938273 - (1478*Pi6)/1701. - (406507*Pi2)/972. + (92237*Pi4)/2430.
                     - (52564*zeta3)/27. + (1690*Pi2*zeta3)/9. + (592*zeta32)/3. - (5512*zeta5)/9.)
           + CA2*CF*(-1946.4519509221154 + (4784*Pi6)/25515. + (596513*Pi2)/2187. - (4303*Pi4)/4860.
                     + (505087*zeta3)/243. - (1168*Pi2*zeta3)/9. - (2272*zeta32)/9. - (868*zeta5)/9.)
           + CF2*nf*(-117.20781893004116 + (13705*Pi2)/243. - (1463*Pi4)/243. + (26080*zeta3)/81.
                     - (148*Pi2*zeta3)/9. - (832*zeta5)/9.)
           + CA*CF*nf*(518.2658131382411 - (201749*Pi2)/2187. - (35*Pi4)/243. - (8576*zeta3)/27.
                       + (148*Pi2*zeta3)/9. - (8*zeta5)/3.)
           + CF3*(-933.1666666666666 + (27403*Pi6)/17010. + (4339*Pi2)/36. - (346*Pi4)/15.
                  - 460*zeta3 - (140*Pi2*zeta3)/3. + 32*zeta32 + 1328*zeta5);
  }

  //_________________________________________________________________________
  double H1SIDIS()
  {
    return 2 * CF * ( - 8 + zeta2 );
  }

  //_________________________________________________________________________
  double H2SIDIS(int const& nf)
  {
    return 2 * CF * ( CF * ( 511. / 8. + 13 * zeta2 - 30 * zeta3 + 39 * zeta4 / 2 ) +
                      CA * ( - 51157. / 648. - 337 * zeta2 / 18 + 313 * zeta3 / 9 + 22 * zeta4 ) +
                      TR * nf * ( 4085. / 162. + 46 * zeta2 / 9 + 4 * zeta3 / 9 ) );
  }

  //_________________________________________________________________________
  double H3SIDIS(int const& nf)
  {
    const int nf2       = nf * nf;
    const double CF2    = CF * CF;
    const double CF3    = CF * CF2;
    const double CA2    = CA * CA;
    const double Pi4    = Pi2 * Pi2;
    const double Pi6    = Pi2 * Pi4;
    const double zeta32 = zeta3 * zeta3;
    return CF*nf2*(-29.100899253162627 - (824*Pi2)/243. - (94*Pi4)/1215. - (832*zeta3)/243.)
           + CA*CF2*(2544.0771604938273 + (292367*Pi2)/972. - (14503*Pi4)/2430. - (2476*Pi6)/8505.
                     - (52564*zeta3)/27. - (382*Pi2*zeta3)/3. + (592*zeta32)/3. - (5512*zeta5)/9.)
           + CA2*CF*(-1946.4519509221154 - (412315*Pi2)/2187. + (22157*Pi4)/4860. - (1538*Pi6)/5103.
                     + (505087*zeta3)/243. + (416*Pi2*zeta3)/9. - (2272*zeta32)/9. - (868*zeta5)/9.)
           + CF2*nf*(-117.20781893004116 - (8567*Pi2)/243. - (131*Pi4)/243. + (26080*zeta3)/81.
                     - (4*Pi2*zeta3)/3. - (832*zeta5)/9.)
           + CA*CF*nf*(518.2658131382411 + (115555*Pi2)/2187. + Pi4/243. - (8576*zeta3)/27.
                       + (4*Pi2*zeta3)/9. - (8*zeta5)/3.)
           + CF3*(-933.1666666666666 - (4859*Pi2)/36. + (4*Pi4)/15. + (1625*Pi6)/3402.
                  - 460*zeta3 + (220*Pi2*zeta3)/3. + 32*zeta32 + 1328*zeta5);
  }

  //_________________________________________________________________________
  double H3Ch()
  {
    const int nc     = 3;
    const double Pi4 = Pi2 * Pi2;
    return CF*(-32/nc + 8*nc - (40*Pi2)/(3.*nc) + (10*nc*Pi2)/3. + (4*Pi4)/(45.*nc)
               - (nc*Pi4)/45. - (112*zeta3)/(3.*nc) + (28*nc*zeta3)/3. + (640*zeta5)/(3.*nc)
               - (160*nc*zeta5)/3.);
  }

  //_________________________________________________________________________
  double H1ggH()
  {
    return 2 * ( CA * ( 5 + 7 * Pi2 / 6 ) - 3 * CF );
  }

  //_________________________________________________________________________
  double H2ggH(int const& nf)
  {
    const double mH = 125.1;
    const double mt = 172.9;
    const double Lt = 2 * log( mH / mt );
    return ( - 135 * CA + 23827 * CA * CA - 216 * CF - 15660 * CA * CF + 5832 * CF * CF
             + 2268 * CA * CA * Lt - 3564 * CA * CF * Lt - 4510 * CA * nf - 4428 * CF * nf
             + 1296 * CF * nf * Lt + 6795 * CA * CA * Pi2 - 2268 * CA * CF * Pi2
             - 450 * Pi2 * CA * nf + 333 * Pi2 * Pi2 * CA * CA - 5148 * CA * CA * zeta3
             - 1656 * zeta3 * CA * nf + 2592 * zeta3 * CF * nf ) / 162;
  }
}
