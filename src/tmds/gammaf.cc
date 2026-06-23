//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/gammaf.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________
  double gammaFq0()
  {
    return 6 * CF;
  }

  //_________________________________________________________________________
  double gammaFq1(int const& nf)
  {
    return - ( CF * CF * ( - 3 + 4 * Pi2 - 48 * zeta3 )
               + CF * CA * ( - 961. / 27. - 11 * Pi2 / 3 + 52 * zeta3 )
               + CF * TR * nf * (  260. / 27. + 4 * Pi2 / 3 ) );
  }

  //_________________________________________________________________________
  double gammaFq2(int const& nf)
  {
    return - ( pow(CF, 3) * ( - 29 - 6 * Pi2 -  16 * Pi2 * Pi2 / 5 - 136 * zeta3 +  32 * Pi2 * zeta3 / 3  + 480 * zeta5 )
               + CF * CF * CA * ( - 151. / 2. +  410 * Pi2 / 9 + 494 * Pi2 * Pi2 / 135 - 1688 * zeta3 / 3 - 16 * Pi2 * zeta3 / 3 - 240 * zeta5 )
               + CF * CA * CA * ( - 139345. / 1458. - 7163 * Pi2 / 243 - 83 * Pi2 * Pi2 / 45 + 7052 * zeta3 / 9 - 88 * Pi2 * zeta3 / 9 - 272 * zeta5 )
               + CF * CF * TR * nf * ( 5906. / 27. - 52 * Pi2 / 9 - 56 * Pi2 * Pi2 / 27 + 1024 * zeta3 / 9 )
               + CF * CA * TR * nf * ( - 34636. / 729. + 5188 * Pi2 / 243 + 44 * Pi2 * Pi2 / 45 - 3856 * zeta3 / 27 )
               + CF * TR * TR * nf * nf * ( 19336. / 729. - 80 * Pi2 / 27 - 64 * zeta3 / 27 ) );
  }

  //_________________________________________________________________________
  double gammaFq3(int const& nf)
  {
    const double H      = - 0.7015802723647;
    const double d4FFNF = 5. / 36.;
    const double d4AFNF = 5. / 2.;
    const double nf2    = nf * nf;
    const double nf3    = nf * nf2;
    const double zeta22 = zeta2 * zeta2;
    const double zeta23 = zeta2 * zeta22;
    const double zeta32 = zeta3 * zeta3;
    const double CF2    = CF * CF;
    const double CF3    = CF * CF2;
    const double CF4    = CF * CF3;
    const double CA2    = CA * CA;
    const double CA3    = CA * CA2;
    return CF*nf3*(-37382/6561 + (16*zeta2)/27 +
                   (128*zeta22)/135 + (1424*zeta3)/243) +
           (d4FFNF*nf)*(-384 + (4544*zeta2)/3 - (320*zeta22)/3 +
                        (9472*zeta23)/315 - (5312*zeta3)/9 + 128*zeta2*zeta3 +
                        (1216*zeta32)/3 - (21760*zeta5)/9) + CA*CF2*nf*
           (-1092511/972 + (673*zeta2)/27 + (105488*zeta22)/135 +
            (5744*zeta23)/35 - (23518*zeta3)/81 + (3904*zeta2*zeta3)/9 -
            (3400*zeta32)/3 - (4472*zeta5)/3) + CA*CF*nf2*
           (-97189/17496 + (41579*zeta2)/729 + (152*zeta22)/15 +
            (14872*zeta3)/243 + (256*zeta2*zeta3)/9 - (1184*zeta5)/9) +
           CF2*nf2*(9965/486 + (1972*zeta2)/27 - (8032*zeta22)/135 -
                    (4232*zeta3)/81 - (224*zeta2*zeta3)/9 + (1040*zeta5)/9) +
           CA2*CF*nf*(326863/1944 - (445117*zeta2)/729 -
                      (17164*zeta22)/45 + (24184*zeta23)/315 + (140632*zeta3)/243 -
                      (3584*zeta2*zeta3)/9 + (6916*zeta32)/9 + (6088*zeta5)/27) +
           CF3*nf*(27949/108 + 322*zeta2 - (668*zeta22)/5 -
                   (117344*zeta23)/315 - (1120*zeta3)/9 - (512*zeta2*zeta3)/3 +
                   368*zeta32 + (3872*zeta5)/3) + CA*CF3*
           (-2085/2 + 320*H + 2334*zeta2 + (8188*zeta22)/5 +
            (527336*zeta23)/315 - 9400*zeta3 + (1784*zeta2*zeta3)/3 +
            (4512*zeta22*zeta3)/5 + 3240*zeta32 + 6048*zeta5 +
            3328*zeta2*zeta5 - 25060*zeta7) + CA3*CF*
           (7179083/26244 + (560*H)/3 + (1062149*zeta2)/729 +
            (179182*zeta22)/135 - (139592*zeta23)/315 -
            (2159464*zeta3)/243 + (25480*zeta2*zeta3)/9 +
            (956*zeta22*zeta3)/5 - (11674*zeta32)/9 +
            (301166*zeta5)/27 + (248*zeta2*zeta5)/3 - (18927*zeta7)/2) +
           d4AFNF*(192 - 640*H - (2176*zeta2)/3 + (3104*zeta22)/15 +
                   (241888*zeta23)/315 + (44032*zeta3)/9 - 5632*zeta2*zeta3 -
                   (8736*zeta22*zeta3)/5 + (15856*zeta32)/3 -
                   (145840*zeta5)/9 + 2624*zeta2*zeta5 + 9924*zeta7) +
           CF4*(4873/12 - 900*zeta2 - (1368*zeta22)/5 -
                (33776*zeta23)/35 + 4008*zeta3 - 240*zeta2*zeta3 +
                (256*zeta22*zeta3)/5 - 2304*zeta32 - 5040*zeta5 -
                768*zeta2*zeta5 + 11760*zeta7) + CA2*CF2*
           (29639/18 - 480*H - (93542*zeta2)/27 - (44792*zeta22)/27 -
            (26136*zeta23)/35 + (375964*zeta3)/27 - (21728*zeta2*zeta3)/
            9 - (6128*zeta22*zeta3)/5 + (196*zeta32)/3 -
            (97292*zeta5)/9 - 3008*zeta2*zeta5 + 22050*zeta7);
  }

  //_________________________________________________________________________
  double gammaFg0(int const& nf)
  {
    return - ( - 22 * CA / 3 + 8 * TR * nf / 3 );
  }

  //_________________________________________________________________________
  double gammaFg1(int const& nf)
  {
    return - ( CA * CA * ( - 1384. / 27. + 11 * Pi2 / 9 + 4 * zeta3 )
               + CA * TR * nf * ( 512. / 27. - 4 * Pi2 / 9 ) + 8 * CF * TR * nf );
  }

  //_________________________________________________________________________
  double gammaFg2(int const& nf)
  {
    return - ( 2 * pow(CA, 3) * ( - 97186. / 729. + 6109 * Pi2 / 486 - 319 * Pi2 * Pi2 / 270 + 122 * zeta3 / 3 - 20 * Pi2 * zeta3 / 9 - 16 * zeta5 )
               + 2 * CA * CA * TR * nf * ( 30715. / 729. - 1198 * Pi2 / 243 + 82 * Pi2 * Pi2 / 135 + 712 * zeta3 / 27 )
               + 2 * CA * CF * TR * nf * ( 2434. / 27. - 2 * Pi2 / 3 - 8 * Pi2 * Pi2 / 45  - 304 * zeta3 / 9 )
               - 4 * CF * CF * TR * nf
               + 2 * CA * TR * TR * nf * nf * ( - 538. / 729. + 40 * Pi2 / 81 - 224 * zeta3 / 27 )
               - 88 * CF * TR * TR * nf * nf / 9 );
  }

  //_________________________________________________________________________
  double gammaFg3(int const& nf)
  {
    const double H      = - 0.7015802723647;
    const double d4FFNA = 5. / 96.;
    const double d4AFNA = 15. / 16.;
    const double d4AANA = 135. / 8.;
    const double nf2    = nf * nf;
    const double nf3    = nf * nf2;
    const double zeta22 = zeta2 * zeta2;
    const double zeta23 = zeta2 * zeta22;
    const double zeta32 = zeta3 * zeta3;
    const double CF2    = CF * CF;
    const double CF3    = CF * CF2;
    const double CA2    = CA * CA;
    const double CA3    = CA * CA2;
    const double CA4    = CA * CA3;
    return 46*CF3*nf + (308*CF*nf3)/243 +
           CF2*nf2*(676/27 - (352*zeta3)/9) +
           CA*nf3*(-15890/6561 - (16*zeta2)/81 + (256*zeta22)/135 -
                   (400*zeta3)/243) + (d4FFNA*nf2)*
           (-1408/9 + (1024*zeta3)/3) + (d4AFNA*nf)*
           (448/9 - 64*zeta2 + (2464*zeta22)/15 - (14464*zeta23)/315 +
            (2560*zeta3)/9 + 1216*zeta2*zeta3 + (1216*zeta32)/3 -
            (30880*zeta5)/9) + CA*CF2*nf*(685/12 - 2*zeta2 +
                                          (148*zeta22)/5 - (320*zeta23)/7 + (1592*zeta3)/3 -
                                          80*zeta32 - (1600*zeta5)/3) + CA2*nf2*
           (611939/17496 - (13483*zeta2)/729 + (3128*zeta22)/135 +
            (37354*zeta3)/243 - 32*zeta2*zeta3 - (1024*zeta5)/9) +
           CA2*CF*nf*(-903983/972 + (3023*zeta2)/9 -
                      (1196*zeta22)/45 + (5632*zeta23)/315 + (29606*zeta3)/81 -
                      176*zeta2*zeta3 + 152*zeta32 + (8*zeta5)/9) +
           CA*CF*nf2*(1199/18 - (172*zeta2)/9 + (128*zeta22)/45 -
                      (1688*zeta3)/81 + (32*zeta2*zeta3)/3 + (304*zeta5)/9) +
           CA3*nf*(-421325/1944 + (155273*zeta2)/729 -
                   (69502*zeta22)/135 + (148976*zeta23)/945 -
                   (260822*zeta3)/243 + 148*zeta2*zeta3 - (596*zeta32)/9 +
                   (16066*zeta5)/27) + CA4*(10672040/6561 + (80*H)/3 -
                                            (1051411*zeta2)/729 + (248368*zeta22)/135 -
                                            (100208*zeta23)/135 - (21940*zeta3)/243 + (2068*zeta2*zeta3)/
                                            3 - (404*zeta22*zeta3)/5 - (2686*zeta32)/9 +
                                            (37232*zeta5)/27 - (1096*zeta2*zeta5)/3 - (1427*zeta7)/2) +
           d4AANA*(128/9 - 640*H + 64*zeta2 + (1072*zeta22)/15 +
                   (253856*zeta23)/315 + (39328*zeta3)/9 - 6176*zeta2*zeta3 -
                   (8736*zeta22*zeta3)/5 + (15856*zeta32)/3 -
                   (141280*zeta5)/9 + 2624*zeta2*zeta5 + 9924*zeta7);
  }
}
