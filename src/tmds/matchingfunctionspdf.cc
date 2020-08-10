//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/matchingfunctionspdf.h"
#include "apfel/constants.h"
#include "apfel/specialfunctions.h"

#include <numeric>

namespace apfel
{
  //_________________________________________________________________________________
  C1nspdf::C1nspdf():
    Expression()
  {
  }
  double C1nspdf::Regular(double const& x) const
  {
    return 2 * CF * ( 1 - x );
  }
  double C1nspdf::Local(double const&) const
  {
    return - CF * zeta2;
  }

  //_________________________________________________________________________________
  C1qgpdf::C1qgpdf():
    Expression()
  {
  }
  double C1qgpdf::Regular(double const& x) const
  {
    return 8 * TR * x * ( 1 - x );
  }

  //_________________________________________________________________________________
  C1gqpdf::C1gqpdf():
    Expression()
  {
  }
  double C1gqpdf::Regular(double const& x) const
  {
    return 2 * CF * x;
  }

  //_________________________________________________________________________________
  C1ggpdf::C1ggpdf():
    Expression()
  {
  }
  double C1ggpdf::Local(double const&) const
  {
    return - CA * zeta2;
  }

  //_________________________________________________________________________________
  C2Vqqpdf::C2Vqqpdf(int const& nf):
    Expression(),
    _nf(nf)
  {
    /*
        _A2 = 224 * CF * _nf * TR / 27.
          + CA * CF * ( - 29.925925925925927 + 28 * zeta3 );
    */
    _A2 = - 3232. / 27. + 112 * zeta3 + 448. * _nf / 81.;
  }
  double C2Vqqpdf::Regular(double const& x) const
  {
    /*
        // Polylogs
        double *xx  = new double{x};
        int    *nw  = new int{3};
        int    *wn1 = new int{-1};
        int    *wn2 = new int{1};
        double *H   = new double[363];
        apf_hplog_(xx, nw, &H[0], &H[3], &H[12], &H[39], &H[120], wn1, wn2);

        // Defintions
        const double x2 = x * x;
        const double H0  = H[HPLogMap({0})];
        const double H1  = H[HPLogMap({1})];

        const double H00 = H[HPLogMap({0,0})];
        const double H2  = H[HPLogMap({2})];
        const double H10 = H[HPLogMap({1,0})];

        const double H000 = H[HPLogMap({0,0,0})];
        const double H12  = H[HPLogMap({1,2})];
        const double H20  = H[HPLogMap({2,0})];
        const double H21  = H[HPLogMap({2,1})];
        const double H100 = H[HPLogMap({1,0,0})];
        const double H110 = H[HPLogMap({1,1,0})];

        // Delete pointers
        delete xx;
        delete nw;
        delete wn1;
        delete wn2;
        delete[] H;

        return CF*_nf*TR*((-4*(37 + 19*x))/27. + (20*(1 + x2)*H0)/(9.*(1 - x)) + (4*(1 + x2)*H00)/(3.*(1 - x)))
          + pow(CF,2)*(-22*(1 - x) + (2*(5 - 13*x + 16*x2)*H0)/(1 - x) + 2*x*H1 - (2*(-3 - 2*x + 2*x2)*H00)/(1 - x)
    		   + 2*(1 + x)*H000 + 2*(1 - x)*(2*H2 + 6*H10 + 3*zeta2)
    		   + (4*(1 + x2)*(2*H12 + 2*H20 + H21 - H100 + 2*H110 + 6*zeta3))/(1 - x))
          + CA*CF*((8*(100 + x))/27. - (2*(29 - 36*x + 83*x2)*H0)/(9.*(1 - x)) - 2*x*H1
    	       + ((-11 - 12*x + x2)*H00)/(3.*(1 - x)) - 2*(1 - x)*(2*H10 + 3*zeta2)
    	       + (-2*(1 + x2)*(2*H12 + 2*H20 + H000 + 2*H110) + 2*(-13 + x2)*zeta3)/(1 - x));
    */
    const double x2   = x * x;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx2 * lx;
    const double l1x  = log(1-x);
    const double l1x2 = l1x * l1x;
    const double l1x3 = l1x2 * l1x;
    const std::vector<double> fReg
    {
      l1x,
      l1x2,
      l1x3,
      1 / x,
      lx / x,
      lx,
      lx2,
      lx3,
      1,
      x,
      x2,
      x * lx / ( 1 - x ),
      x * lx,
      x2 *lx,
      x * lx2 / ( 1 - x ),
      x * lx2,
      ( lx / ( 1 - x ) + 1 ) * l1x,
      lx * l1x,
      x * lx * l1x,
      ( 1 - x ) * l1x / x,
      ( 1 - x ) * l1x,
      ( 1 - x ) * ( 1 - x ) * l1x,
      ( 1 - x ) * l1x2};

    const std::vector<double> CoeffReg
    {
      200. / 9.,
      - 64. / 9.,
      0.,
      0.,
      0.,
      - 8. + 40. * _nf/27.,
      - 2. + 4. * _nf/9.,
      - 20. / 27.,
      - 296. * _nf / 81. + 1076.6744297226016,
      - 152. * _nf / 81. + 7792.719665777814,
      111.49810429898287,
      80. * _nf / 27. + 8980.334190376141,
      - 40. * _nf / 27. - 3795.008745809993,
      82.30795871692112,
      8. * _nf / 9. - 201.0129463471822,
      - 4. * _nf / 9. + 206.75145891009598,
      5603.371344939401,
      - 526.1352578350587,
      1382.8610999663256,
      1092.9256332669593,
      2547.784733022028,
      - 147.17479558391307,
      3.564983084988843};

    return std::inner_product(fReg.begin(), fReg.end(), CoeffReg.begin(), 0.);
  }
  double C2Vqqpdf::Singular(double const& x) const
  {
    return _A2 / ( 1 - x );
  }
  double C2Vqqpdf::Local(double const& x) const
  {
    const double ln1mx = log( 1 - x );
    const double A1    = - 2416. / 81. - 134 * zeta2 / 3 + 448 * zeta3 / 9 + 200 * zeta4 / 9
                         + _nf * ( 352. / 243. + 20 * zeta2 / 9 + 56 * zeta3 / 27 );

    return A1 + _A2 * ln1mx;
  }

  //_________________________________________________________________________________
  C2Vqqbpdf::C2Vqqbpdf():
    Expression()
  {
  }
  double C2Vqqbpdf::Regular(double const& x) const
  {
    /*
        // Polylogs
        double *xx  = new double{x};
        int    *nw  = new int{3};
        int    *wn1 = new int{-1};
        int    *wn2 = new int{1};
        double *H   = new double[363];
        apf_hplog_(xx, nw, &H[0], &H[3], &H[12], &H[39], &H[120], wn1, wn2);

        // Defintions
        const double Hm1 = H[HPLogMap({-1})];
        const double H0  = H[HPLogMap({0})];

        const double Hm10 = H[HPLogMap({-1,0})];
        const double H10  = H[HPLogMap({1,0})];

        const double Hm20   = H[HPLogMap({-2,0})];
        const double H000   = H[HPLogMap({0,0,0})];
        const double H20    = H[HPLogMap({2,0})];
        const double Hm1m10 = H[HPLogMap({-1,-1,0})];
        const double Hm100  = H[HPLogMap({-1,0,0})];

        // Delete pointers
        delete xx;
        delete nw;
        delete wn1;
        delete wn2;
        delete[] H;

        return (CA - 2*CF)*CF*(-15*(1 - x) + (-3 - 11*x)*H0 + 4*(1 + x)*Hm10
    			   + 4*(1 - x)*H10 - 2*(-3 + x)*zeta2
    			   - (2*(1 + pow(x,2))*(4*Hm20 - 2*H20 - 4*Hm1m10 + 2*Hm100
    						- H000 - 2*Hm1*zeta2 + zeta3))/(1 + x));
    */
    const double x2   = x * x;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx2 * lx;
    const double l1x  = log(1-x);
    const double l1x2 = l1x * l1x;
    const double l1x3 = l1x2 * l1x;
    const std::vector<double> fReg
    {
      l1x,
      l1x2,
      l1x3,
      1 / x,
      lx / x,
      lx,
      lx2,
      lx3,
      1,
      x,
      x2,
      x * lx / ( 1 - x ),
      x * lx,
      x2 *lx,
      x * lx2 / ( 1 - x ),
      x * lx2,
      ( lx / ( 1 - x ) + 1 ) * l1x,
      lx * l1x,
      x * lx * l1x,
      ( 1 - x ) * l1x / x,
      ( 1 - x ) * l1x,
      ( 1 - x ) * ( 1 - x ) * l1x,
      ( 1 - x ) * l1x2};

    const std::vector<double> CoeffReg
    {
      0.,
      0.,
      0.,
      0.,
      0.,
      - 4. / 3.,
      0.,
      4. / 27.,
      540.6779037661292,
      4561.881499265996,
      330.3186826846845,
      5432.878085716809,
      - 2563.28806902233,
      - 17.78991469587654,
      - 76.36755190995123,
      78.87763768206408,
      3443.1424448947796,
      - 599.7983485750569,
      839.4963238597323,
      544.0265746128981,
      1417.3624064262308,
      - 113.40180701234924,
      0.009338040302161874};

    return std::inner_product(fReg.begin(), fReg.end(), CoeffReg.begin(), 0.);
  }

  //_________________________________________________________________________________
  C2pspdf::C2pspdf():
    Expression()
  {
  }
  double C2pspdf::Regular(double const& x) const
  {
    /*
        // Polylogs
        double *xx  = new double{x};
        int    *nw  = new int{3};
        int    *wn1 = new int{-1};
        int    *wn2 = new int{1};
        double *H   = new double[363];
        apf_hplog_(xx, nw, &H[0], &H[3], &H[12], &H[39], &H[120], wn1, wn2);

        // Defintions
        const double x2 = x * x;
        const double H0   = H[HPLogMap({0})];
        const double H00  = H[HPLogMap({0,0})];
        const double H10  = H[HPLogMap({1,0})];
        const double H000 = H[HPLogMap({0,0,0})];

        // Delete pointers
        delete xx;
        delete nw;
        delete wn1;
        delete wn2;
        delete[] H;

        return 2*CF*TR*((2*(1 - x)*(172 - 143*x + 136*x2))/(27.*x)
    		    + (4*(21 - 30*x + 32*x2)*H0)/9. - (2*(3 + 3*x + 8*x2)*H00)/3.
    		    + 4*(1 + x)*H000 - (8*(1 - x)*(2 - x + 2*x2)*(H10 + zeta2))/(3.*x));
    */
    /*
      return 2*CF*TR*((2*(1 - x)*(172 - 143*x + 136*x2))/(27.*x)
      + (4*(21 - 30*x + 32*x2)*log(x))/9. - ((3 + 3*x + 8*x2)*pow(log(x),2))/3.
      + (2*(1 + x)*pow(log(x),3))/3. - (8*(1 - x)*(2 - x + 2*x2)*(-(log(1 - x)*log(x))
      - dilog(x) + zeta2))/(3.*x));
    */
    const double x2   = x * x;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx2 * lx;
    const double l1x  = log(1-x);
    const double l1x2 = l1x * l1x;
    const double l1x3 = l1x2 * l1x;
    const std::vector<double> fReg
    {
      l1x,
      l1x2,
      l1x3,
      1 / x,
      lx / x,
      lx,
      lx2,
      lx3,
      1,
      x,
      x2,
      x * lx / ( 1 - x ),
      x * lx,
      x2 *lx,
      x * lx2 / ( 1 - x ),
      x * lx2,
      ( lx / ( 1 - x ) + 1 ) * l1x,
      lx * l1x,
      x * lx * l1x,
      ( 1 - x ) * l1x / x,
      ( 1 - x ) * l1x,
      ( 1 - x ) * ( 1 - x ) * l1x,
      ( 1 - x ) * l1x2};

    const std::vector<double> CoeffReg
    {
      0.,
      0.,
      0.,
      688. / 81. - 32. * zeta2 / 9.,
      0.,
      8. / 3.,
      - 2. / 3.,
      4. / 9.,
      - 603.924227035499,
      - 4636.485211211568,
      - 49.76555465398209,
      - 5287.52982020046,
      2269.612280502602,
      - 58.06494427244658,
      119.75596348356056,
      - 129.79997791500733,
      - 3369.724995744234,
      427.8946185110229,
      - 812.4665998422748,
      - 600.6972087253562,
      - 1469.0061919285804,
      92.73615445019001,
      - 0.021528986881156446};

    return 2 * std::inner_product(fReg.begin(), fReg.end(), CoeffReg.begin(), 0.);
  }

  //_________________________________________________________________________________
  C2qgpdf::C2qgpdf():
    Expression()
  {
  }
  double C2qgpdf::Regular(double const& x) const
  {
    /*
        // Polylogs
        double *xx  = new double{x};
        int    *nw  = new int{3};
        int    *wn1 = new int{-1};
        int    *wn2 = new int{1};
        double *H   = new double[363];
        apf_hplog_(xx, nw, &H[0], &H[3], &H[12], &H[39], &H[120], wn1, wn2);

        // Defintions
        const double x2 = x * x;
        const double x3 = x * x2;
        const double Hm1 = H[HPLogMap({-1})];
        const double H0  = H[HPLogMap({0})];
        const double H1  = H[HPLogMap({1})];

        const double Hm10 = H[HPLogMap({-1,0})];
        const double H00  = H[HPLogMap({0,0})];
        const double H2   = H[HPLogMap({2})];
        const double H10  = H[HPLogMap({1,0})];
        const double H11  = H[HPLogMap({1,1})];

        const double Hm20   = H[HPLogMap({-2,0})];
        const double H000   = H[HPLogMap({0,0,0})];
        const double H12    = H[HPLogMap({1,2})];
        const double H20    = H[HPLogMap({2,0})];
        const double H21    = H[HPLogMap({2,1})];
        const double H100   = H[HPLogMap({1,0,0})];
        const double H110   = H[HPLogMap({1,1,0})];
        const double Hm1m10 = H[HPLogMap({-1,-1,0})];
        const double Hm100  = H[HPLogMap({-1,0,0})];
        const double H111   = H[HPLogMap({1,1,1})];

        // Delete pointers
        delete xx;
        delete nw;
        delete wn1;
        delete wn2;
        delete[] H;

        return 2*(CA*TR*((-2*(-172 + 315*x - 387*x2 + 298*x3))/(27.*x) + (4*(21 - 30*x + 68*x2)*H0)/9.
    		     + 2*x*(-3 + 4*x)*H1 + 8*x*(1 + x)*Hm10 - (2*(3 - 12*x + 44*x2)*H00)/3.
    		     - (8*(1 - x)*(2 - x + 11*x2)*H10)/(3.*x) - 8*(1 - x)*x*H11 + 4*(1 + 2*x)*H000
    		     + 4*(1 - 2*x + 2*x2)*(H12 + H110 - H111) + (8*(-2 + 3*x - 9*x2 + 11*x3)*zeta2)/(3.*x)
    		     + 4*(1 + 2*x + 2*x2)*(2*Hm20 - 2*Hm1m10 + Hm100 - Hm1*zeta2) - 8*x*(2*H20 - zeta3))
    	      + CF*TR*(-13 + 75*x - 72*x2 + (8 + 15*x - 8*x2)*H0 - 2*x*(-3 + 4*x)*H1 + (1 + 12*x - 8*x2)*H00
    		       - 2*(1 - 2*x + 4*x2)*H000 + 4*(1 - x)*x*(2*H2 + 2*H10 + 2*H11 - 3*zeta2)
    		       + 4*(1 - 2*x + 2*x2)*(H21 - H100 + H111 + 7*zeta3)));
    */
    const double x2   = x * x;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx2 * lx;
    const double l1x  = log(1-x);
    const double l1x2 = l1x * l1x;
    const double l1x3 = l1x2 * l1x;
    const std::vector<double> fReg
    {
      l1x,
      l1x2,
      l1x3,
      1 / x,
      lx / x,
      lx,
      lx2,
      lx3,
      1,
      x,
      x2,
      x * lx / ( 1 - x ),
      x * lx,
      x2 *lx,
      x * lx2 / ( 1 - x ),
      x * lx2,
      ( lx / ( 1 - x ) + 1 ) * l1x,
      lx * l1x,
      x * lx * l1x,
      ( 1 - x ) * l1x / x,
      ( 1 - x ) * l1x,
      ( 1 - x ) * ( 1 - x ) * l1x,
      ( 1 - x ) * l1x2};

    const std::vector<double> CoeffReg
    {
      - 5. / 3.,
        0.,
        5. / 9.,
        172. / 9. - 8. * zeta2,
        0.,
        34. / 3.,
        - 7. / 6.,
        7. / 9.,
        11804.917158263002,
        43420.91122559322,
        - 3972.3493167408856,
        51274.90217179903,
        - 11322.606629412014,
        1490.8267572457225,
        365.38010578193774,
        - 402.6162297779727,
        21027.101170409886,
        19269.1080641264,
        11357.814388977074,
        11798.406220601846,
        29476.38189858057,
        - 227.4273932698467,
        15.408922558612357};

    return 2 * std::inner_product(fReg.begin(), fReg.end(), CoeffReg.begin(), 0.);
  }

  //_________________________________________________________________________________
  C2gqpdf::C2gqpdf(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C2gqpdf::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx2 * lx;
    const double l1x  = log(1-x);
    const double l1x2 = l1x * l1x;
    const double l1x3 = l1x2 * l1x;
    const std::vector<double> fReg
    {
      l1x,
      l1x2,
      l1x3,
      1 / x,
      lx / x,
      lx,
      lx2,
      lx3,
      1,
      x,
      x2,
      x * lx / ( 1 - x ),
      x * lx,
      x2 *lx,
      x * lx2 / ( 1 - x ),
      x * lx2,
      ( lx / ( 1 - x ) + 1 ) * l1x,
      lx * l1x,
      x * lx * l1x,
      ( 1 - x ) * l1x / x,
      ( 1 - x ) * l1x,
      ( 1 - x ) * ( 1 - x ) * l1x,
      ( 1 - x ) * l1x2};

    const std::vector<double> CoeffReg
    {
      - 184. / 9. + 32. * _nf / 27.,
        - 44. / 9. + 8. * _nf / 9.,
        - 40. / 27.,
        - 12640. / 27. + 896. * _nf / 81. + 352. * zeta2 / 3. + 192. * zeta3,
        0.,
        - 200. / 3.,
        112. / 9.,
        - 112. / 27.,
        25387.766684267783 - 155.5078210719438 * _nf,
        128641.80191589122 - 529.7740230354465 * _nf,
        - 4691.90206636766 + 33.2493657381115 * _nf,
        149304.47081486747 - 643.5386512087849 * _nf,
        - 49319.02516115282 + 160.4570922525814 * _nf,
        3077.798040737399 - 15.139140368040819 * _nf,
        - 1297.381684493681 + 2.266359003377411 * _nf,
        1379.187692679483 - 2.2865882004537834 * _nf,
        79157.51026314059 - 313.45913506070985 * _nf,
        20294.80059647656 - 169.33020796238313 * _nf,
        27727.02916158908 - 115.18065110064754 * _nf,
        25431.533693585556 - 138.52016675095615 * _nf,
        62304.61943144707 - 297.70931938249043 * _nf,
        - 1802.796656932475 + 4.76862827955143 * _nf,
        5.352770139788647 + 0.9080131334622416 * _nf};

    return std::inner_product(fReg.begin(), fReg.end(), CoeffReg.begin(), 0.);
  }

  //_________________________________________________________________________________
  C2ggpdf::C2ggpdf(int const& nf):
    Expression(),
    _nf(nf)
  {
    _A2 = - 808. / 3. + 252. * zeta3 + 112. * _nf / 9.;
  }
  double C2ggpdf::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx2 * lx;
    const double l1x  = log(1-x);
    const double l1x2 = l1x * l1x;
    const double l1x3 = l1x2 * l1x;
    const std::vector<double> fReg
    {
      l1x,
      l1x2,
      l1x3,
      1 / x,
      lx / x,
      lx,
      lx2,
      lx3,
      1,
      x,
      x2,
      x * lx / ( 1 - x ),
      x * lx,
      x2 *lx,
      x * lx2 / ( 1 - x ),
      x * lx2,
      ( lx / ( 1 - x ) + 1 ) * l1x,
      lx * l1x,
      x * lx * l1x,
      ( 1 - x ) * l1x / x,
      ( 1 - x ) * l1x,
      ( 1 - x ) * ( 1 - x ) * l1x,
      ( 1 - x ) * l1x2};

    const std::vector<double> CoeffReg
    {
      6. - 2. * _nf,
      - 36.,
      0.,
      - 3160. / 3. + 226. * _nf / 9. + 264. * zeta2 + 432. * zeta3,
      0.,
      - 293. + 74. * _nf / 3.,
      3. + 6. * _nf,
      - 12. + 8. * _nf / 9.,
      28350.905925013587 - 1165.597999487027 * _nf,
      204106.38400427837 - 9112.918558788755 * _nf,
      - 3654.296204682914 - 202.82515955533333 * _nf,
      228606.3613205384 - 10439.786162275559 * _nf,
      - 92515.7827039022 + 4572.107887712472 * _nf,
      4598.901824247715 - 82.57482307730146 * _nf,
      -5961.66732755631 + 242.87579888233213 * _nf,
      6430.443013955799 - 258.29777371318767 * _nf,
      144229.2381370356 - 6707.548773597168 * _nf,
      - 11622.452033932963 + 958.7180437906424 * _nf,
      35791.91490687339 - 1542.1541360841081 * _nf,
      28788.30003557479 - 1171.3757772648048 * _nf,
      67492.96763028382 - 2764.3635575309518 * _nf,
      - 3881.2814494298837 + 176.34274530829245 * _nf,
      18.352295609756112 - 0.04234942977167983 * _nf};

    return std::inner_product(fReg.begin(), fReg.end(), CoeffReg.begin(), 0.);
  }
  double C2ggpdf::Singular(double const& x) const
  {
    return _A2 / ( 1 - x );
  }
  double C2ggpdf::Local(double const& x) const
  {
    const double ln1mx = log(1-x);
    const double A1    = - 112. - 201. * zeta2 / 2. + 154. * zeta3 + 225. * zeta4 / 4.
                         + _nf * ( 548. / 27. + 5. * zeta2 - 28. * zeta3 / 3. )
                         - 56. * _nf * _nf  / 81.;

    return A1 + _A2 * ln1mx;
  }

  //_________________________________________________________________________________
  C3Vqqpdf::C3Vqqpdf(int const& nf):
    Expression(),
    _nf(nf)
  {
    const double CF2 = CF * CF;
    const double CA2 = CA * CA;
    const double TR2 = TR * TR;
    const double nf2 = _nf * _nf;
    _A2 = CF*nf2*TR2*(-10.183813443072703 - (128*zeta3)/9.)
          + CA*CF2*((808*zeta2)/27. - 28*zeta2*zeta3)
          + CF2*_nf*TR*(126.74074074074075 - (224*zeta2)/27. - (608*zeta3)/9. - 32*zeta4)
          + CA*CF*_nf*TR*(171.81344307270234 - (1648*zeta2)/81. - (1808*zeta3)/27. + (40*zeta4)/3.)
          + CA2*CF*(-407.4471879286694 + (6392*zeta2)/81. + (12328*zeta3)/27.
                    - (176*zeta2*zeta3)/3. + (154*zeta4)/3. - 192*zeta5);
  }
  double C3Vqqpdf::Regular(double const& x) const
  {
    // Polylogs
    double *xx  = new double{x};
    int    *nw  = new int{5};
    int    *wn1 = new int{-1};
    int    *wn2 = new int{1};
    double *H   = new double[363];
    apf_hplog_(xx, nw, &H[0], &H[3], &H[12], &H[39], &H[120], wn1, wn2);

    // Defintions
    const double x2  = x * x;
    const double x3  = x * x2;
    const double CF2 = CF * CF;
    const double CF3 = CF * CF2;
    const double CA2 = CA * CA;
    const double TR2 = TR * TR;
    const double nf2 = _nf * _nf;

    const double Hm1 = H[HPLogMap({-1})];
    const double H0  = H[HPLogMap({0})];
    const double H1  = H[HPLogMap({1})];

    const double Hm10  = H[HPLogMap({-1,0})];
    const double H00   = H[HPLogMap({0,0})];
    const double H2    = H[HPLogMap({2})];
    const double H10   = H[HPLogMap({1,0})];
    const double H11   = H[HPLogMap({1,1})];
    const double Hm2   = H[HPLogMap({-2})];
    const double Hm1m1 = H[HPLogMap({-1,-1})];

    const double Hm20   = H[HPLogMap({-2,0})];
    const double H000   = H[HPLogMap({0,0,0})];
    const double H12    = H[HPLogMap({1,2})];
    const double H20    = H[HPLogMap({2,0})];
    const double H21    = H[HPLogMap({2,1})];
    const double H100   = H[HPLogMap({1,0,0})];
    const double H110   = H[HPLogMap({1,1,0})];
    const double Hm1m10 = H[HPLogMap({-1,-1,0})];
    const double Hm100  = H[HPLogMap({-1,0,0})];
    const double H111   = H[HPLogMap({1,1,1})];
    const double H3     = H[HPLogMap({3})];
    const double Hm12   = H[HPLogMap({-1,2})];
    const double Hm3    = H[HPLogMap({-3})];
    const double Hm2m1  = H[HPLogMap({-2,-1})];
    const double H1m2   = H[HPLogMap({1,-2})];

    const double H4       = H[HPLogMap({4})];
    const double H22      = H[HPLogMap({2,2})];
    const double H30      = H[HPLogMap({3,0})];
    const double H210     = H[HPLogMap({2,1,0})];
    const double H13      = H[HPLogMap({1,3})];
    const double H31      = H[HPLogMap({3,1})];
    const double H112     = H[HPLogMap({1,1,2})];
    const double H120     = H[HPLogMap({1,2,0})];
    const double H121     = H[HPLogMap({1,2,1})];
    const double H200     = H[HPLogMap({2,0,0})];
    const double H0000    = H[HPLogMap({0,0,0,0})];
    const double H1000    = H[HPLogMap({1,0,0,0})];
    const double H1100    = H[HPLogMap({1,1,0,0})];
    const double H1110    = H[HPLogMap({1,1,1,0})];
    const double H211     = H[HPLogMap({2,1,1})];
    const double Hm22     = H[HPLogMap({-2,2})];
    const double Hm2m10   = H[HPLogMap({-2,-1,0})];
    const double Hm200    = H[HPLogMap({-2,0,0})];
    const double Hm13     = H[HPLogMap({-1,3})];
    const double Hm1m20   = H[HPLogMap({-1,-2,0})];
    const double Hm1m12   = H[HPLogMap({-1,-1,2})];
    const double Hm120    = H[HPLogMap({-1,2,0})];
    const double Hm121    = H[HPLogMap({-1,2,1})];
    const double Hm1m1m10 = H[HPLogMap({-1,-1,-1,0})];
    const double Hm1m100  = H[HPLogMap({-1,-1,0,0})];
    const double Hm1000   = H[HPLogMap({-1,0,0,0})];
    const double Hm30     = H[HPLogMap({-3,0})];
    const double H1m20    = H[HPLogMap({1,-2,0})];
    const double Hm130    = H[HPLogMap({-1,3,0})];
    const double Hm1200   = H[HPLogMap({-1,2,0,0})];
    const double Hm1210   = H[HPLogMap({-1,2,1,0})];

    const double H5       = H[HPLogMap({5})];
    const double Hm40     = H[HPLogMap({-4,0})];
    const double H23      = H[HPLogMap({2,3})];
    const double H41      = H[HPLogMap({4,1})];
    const double Hm300    = H[HPLogMap({-3,0,0})];
    const double H2m20    = H[HPLogMap({2,-2,0})];
    const double H220     = H[HPLogMap({2,2,0})];
    const double H300     = H[HPLogMap({3,0,0})];
    const double H2000    = H[HPLogMap({2,0,0,0})];
    const double H2100    = H[HPLogMap({2,1,0,0})];
    const double H2110    = H[HPLogMap({2,1,1,0})];
    const double H00000   = H[HPLogMap({0,0,0,0,0})];
    const double H212     = H[HPLogMap({2,1,2})];
    const double Hm32     = H[HPLogMap({-3,2})];
    const double Hm23     = H[HPLogMap({-2,3})];
    const double H14      = H[HPLogMap({1,4})];
    const double H32      = H[HPLogMap({3,2})];
    const double H40      = H[HPLogMap({4,0})];
    const double Hm3m10   = H[HPLogMap({-3,-1,0})];
    const double Hm2m20   = H[HPLogMap({-2,-2,0})];
    const double Hm2m12   = H[HPLogMap({-2,-1,2})];
    const double Hm220    = H[HPLogMap({-2,2,0})];
    const double Hm221    = H[HPLogMap({-2,2,1})];
    const double H1m30    = H[HPLogMap({1,-3,0})];
    const double H1m22    = H[HPLogMap({1,-2,2})];
    const double H113     = H[HPLogMap({1,1,3})];
    const double H122     = H[HPLogMap({1,2,2})];
    const double H130     = H[HPLogMap({1,3,0})];
    const double H131     = H[HPLogMap({1,3,1})];
    const double H221     = H[HPLogMap({2,2,1})];
    const double H310     = H[HPLogMap({3,1,0})];
    const double Hm2m1m10 = H[HPLogMap({-2,-1,-1,0})];
    const double Hm2m100  = H[HPLogMap({-2,-1,0,0})];
    const double Hm2000   = H[HPLogMap({-2,0,0,0})];
    const double H1m2m10  = H[HPLogMap({1,-2,-1,0})];
    const double H1m200   = H[HPLogMap({1,-2,0,0})];
    const double H11m20   = H[HPLogMap({1,1,-2,0})];
    const double H1112    = H[HPLogMap({1,1,1,2})];
    const double H1120    = H[HPLogMap({1,1,2,0})];
    const double H1121    = H[HPLogMap({1,1,2,1})];
    const double H1200    = H[HPLogMap({1,2,0,0})];
    const double H1210    = H[HPLogMap({1,2,1,0})];
    const double H10000   = H[HPLogMap({1,0,0,0,0})];
    const double H11000   = H[HPLogMap({1,1,0,0,0})];
    const double H11100   = H[HPLogMap({1,1,1,0,0})];
    const double H311     = H[HPLogMap({3,1,1})];

    // Delete pointers
    delete xx;
    delete nw;
    delete wn1;
    delete wn2;
    delete[] H;

    return CF*nf2*TR2*((64*(-41 + 157*x))/729. - (32*(43 - 34*x + 43*x2)*H0)/(81.*(1 - x))
                       - (32*(37 - 24*x + 37*x2)*H00)/(81.*(1 - x)) - (160*(1 + x2)*H000)/(27.*(1 - x))
                       + (64*(1 + x)*zeta3)/9.)
           + CA*CF*_nf*TR*((-346*(-89 + 451*x))/729. + (80*(168 - 121*x + 157*x2)*H0)/(81.*(1 - x))
                           + (64*(1 + 5*x)*H1)/27. - (8*H11)/3. - (16*(1 + x)*H21)/9. + (8*(-25 + 18*x)*H1*zeta2)/9.
                           + (-24*(-25 + 27*x + x2)*H2 + 8*(1247 - 882*x + 1322*x2)*H00 + 144*(4 - 6*x + x2)*H10
                              + 8*(175 - 231*x + 47*x2)*zeta2)/(81.*(1 - x))
                           + (8*(28 - 33*x + 16*x2)*H3 + 24*(37 - 27*x + 30*x2)*H12 + 24*(26 - 17*x + 34*x2)*H20
                              + 16*(98 + 21*x + 41*x2)*H000 + 96*(4 - 13*x + 4*x2)*H100 + 24*(9 + 15*x + 16*x2)*H110
                              - 8*(20 + 3*x + 8*x2)*H0*zeta2 + 8*(391 - 51*x + 87*x2)*zeta3)/(27.*(1 - x))
                           + (16*(4 + 3*x + x2)*(H4 + H22) + 48*(3 - x + 4*x2)*H30 + 16*(2 - 3*x + 5*x2)*H210
                              + 16*(-2 - 3*x + x2)*H2*zeta2 - 16*(5 + 3*x + 2*x2)*H00*zeta2 + 16*(13 + 9*x + 4*x2)*H0*zeta3
                              + 16*(1 + x2)*(9*H13 - H31 - 4*H112 + 3*H120 - 3*H121 + H200 + 8*H0000 - 8*H1000 + 5*H1100 - 8*H1110
                                             - 7*H10*zeta2 - 4*H11*zeta2 - 5*H1*zeta3) + 12*(-3 - 14*x + 21*x2)*zeta4)/(9.*(1 - x)))
           + CF2*_nf*TR*((-2389 - 1033*x)/27. - (4*(1470 - 1628*x + 1127*x2)*H0)/(81.*(1 - x))
                         - (104*(2 + x)*H1)/27. + (32*x*H11)/9. + (224*(1 - x)*H1*zeta2)/9.
                         + (8*(191 + 174*x + 146*x2)*H2 - 16*(503 - 474*x + 349*x2)*H00 - 16*(53 - 312*x + 35*x2)*H10
                            - 24*(82 - 156*x + 95*x2)*zeta2)/(81.*(1 - x))
                         + (16*(29 + 12*x + 8*x2)*H3 - 96*(13 - 6*x + 13*x2)*H12 - 16*(73 - 72*x + 88*x2)*H20
                            - 192*(3 - x + 3*x2)*H21 + 8*(-169 - 126*x + 191*x2)*H000 + 16*(1 + 138*x + x2)*H100
                            - 192*(3 + 4*x + 3*x2)*H110 - 4*(163 - 96*x + 79*x2)*H0*zeta2
                            - 8*(421 - 234*x + 667*x2)*zeta3)/(27.*(1 - x))
                         + (16*(-5 + 17*x2)*H0000 + 4*(1 + x2)*(24*H4 - 40*H13 - 32*H22 - 56*H30 - 16*H31 + 48*H112
                                                                - 40*H120 + 32*H121 + 8*H200 - 24*H210 + 16*H211 + 104*H1000
                                                                - 32*H1100 + 64*H1110 + 16*H2*zeta2 - 19*H00*zeta2 + 24*H10*zeta2
                                                                + 16*H11*zeta2 - 176*H0*zeta3 - 8*H1*zeta3)
                            - 8*(23 + 59*x2)*zeta4)/(9.*(1 - x)))
           + CA2*CF*((-542581 + 1136639*x)/1458. + (2*(-15080 - 12422*x + 12863*x2)*H0)/(81.*(1 - x))
                     - (2*(-4526 + 5369*x)*H1)/27. - (32*(4 + 7*x)*Hm22)/3. + (4*(1 + x)*(3 + 64*x + 6*x2)*Hm10)/(3.*x)
                     + (52*x*H11)/3. + (4*(14 + 11*x)*H21)/9. - (32*(-2 + 9*x)*Hm2m10)/3. - (160*(1 - x)*Hm200)/3.
                     + (80*(2 + x)*Hm2*zeta2)/3. - (2*(-431 + 579*x)*H1*zeta2)/9.
                     + (-6*(-1492 - 909*x + 2512*x2)*H2 + 6*(-2974 + 5026*x - 6333*x2 + 108*x3)*H00
                        + 18*(86 - 339*x + 185*x2)*H10 - 2*(3605 + 4821*x - 8471*x2 + 324*x3)*zeta2)/(81.*(1 - x))
                     - (8*(1 + x)*(14*Hm12 - 44*Hm1m10 + 45*Hm100 - 36*Hm1*zeta2))/3.
                     + (-2*(-100 + 492*x + 281*x2)*H3 + 36*(27 - 136*x + 127*x2)*Hm20 - 6*(527 - 642*x + 651*x2)*H12
                        - 6*(232 - 406*x + 623*x2)*H20 - 30*(103 - 64*x + 149*x2)*H000 - 12*(349 - 847*x + 340*x2)*H100
                        - 6*(309 - 18*x + 245*x2)*H110 + 2*(448 - 264*x + 1225*x2)*H0*zeta2
                        + 2*(-2183 - 4506*x + 2598*x2)*zeta3)/(27.*(1 - x))
                     - (32*(1 + x)*(3*Hm13 + 2*Hm1m20 - 7*Hm1m12 - 2*Hm120 - Hm121 + Hm1m1m10 - 2*Hm1m100 + 4*Hm1000
                                    + (15*Hm1m1*zeta2)/2. - 3*Hm10*zeta2 - 6*Hm1*zeta3))/3.
                     + (-4*(68 - 39*x + 59*x2)*H4 - 24*(-8 - 25*x + 18*x2)*Hm30 + 12*(40 - 113*x + 31*x2)*H13
                        - 8*(49 - 45*x + 49*x2)*H22 + 4*(23 - 39*x + 56*x2)*H31 + 48*(-1 + 2*x)*(-14 + 11*x)*H1m20
                        - 4*(-14 - 75*x + x2)*H112 - 24*(2 + 3*x + 9*x2)*H120 + 12*(9 - x + 14*x2)*H121
                        - 4*(64 + 72*x + 67*x2)*H0000 + 4*(37 + 84*x + 19*x2)*H1000 - 4*(121 - 120*x + 109*x2)*H1100
                        + 352*(1 + x2)*H1110 - 4*(43 - 66*x + 13*x2)*H10*zeta2 + 8*(31 - 27*x + 40*x2)*H11*zeta2
                        + 4*(133 - 330*x + 145*x2)*H1*zeta3)/(9.*(1 - x))
                     + (-12*(15 + 65*x + 54*x2 + 10*x3)*H30 + 4*(-83 - 110*x + 124*x2 + 115*x3)*H200
                        - 4*(85 + 76*x + 31*x2 + 4*x3)*H210 + 8*(50 + 23*x - 28*x2 + 17*x3)*H2*zeta2
                        + 4*(106 + 37*x - 23*x2 + 28*x3)*H00*zeta2 - 4*(-127 - 202*x + 59*x2 + 98*x3)*H0*zeta3
                        - 3*(-240 + 313*x + 680*x2 + 115*x3)*zeta4)/(9.*(1 - x)*(1 + x))
                     + 2*(1 + x)*(2*Hm130 + 4*Hm1200 - 4*Hm1210 - 4*Hm12*zeta2 + 2*Hm100*zeta2 - 4*Hm10*zeta3 - Hm1*zeta4)
                     + (-16*(5 + x2)*H5 - 32*(2 + x2)*Hm40 - 4*(13 + 31*x2)*H23 + 4*(7 + x2)*H41 - 16*(7 + 5*x2)*Hm300
                        - 16*(9 + 7*x2)*H2m20 + 16*(1 + 2*x2)*H220 + 4*(13 + 15*x2)*H300 + 4*(35 + 33*x2)*H2000
                        + 4*(31 + 49*x2)*H2100 + 4*(-17 + 13*x2)*H2110 + 12*(3 + 5*x2)*H00000 + 8*(21 + 17*x2)*H3*zeta2
                        + 4*(9 + 35*x2)*H20*zeta2 + 64*x2*H21*zeta2 + 4*(-17 + x2)*(H212 - H000*zeta2) - 20*(13 + 11*x2)*H2*zeta3
                        - 4*(11 + 37*x2)*H00*zeta3 + 4*(95 + 49*x2)*zeta2*zeta3 - 4*(42 + 31*x2)*H0*zeta4
                        - 4*(1 + x2)*(64*Hm32 + 24*Hm23 + 6*H14 + 12*H32 - 12*H40 + 48*Hm3m10 + 16*Hm2m20 - 56*Hm2m12
                                      - 16*Hm220 - 8*Hm221 + 72*H1m30 + 8*H1m22 - 32*H113 + 16*H122 - 44*H130 - 16*H131
                                      + 2*H221 + 6*H310 + 8*Hm2m1m10 - 16*Hm2m100 + 32*Hm2000 - 24*H1m2m10 + 64*H1m200
                                      - 40*H11m20 + 16*H1112 - 12*H1120 + 8*H1121 - 30*H1200 + 2*H1210 - 32*H10000
                                      - 28*H11000 - 16*H11100 - 40*Hm3*zeta2 + 60*Hm2m1*zeta2 - 24*Hm20*zeta2
                                      - 20*H1m2*zeta2 - 40*H12*zeta2 + 2*H100*zeta2 - 12*H110*zeta2 - 16*H111*zeta2
                                      - 48*Hm2*zeta3 + 38*H10*zeta3 - 12*H11*zeta3 + 21*H1*zeta4) - 4*(-1 + 75*x2)*zeta5)/(3.*(1 - x)))
           + CA*CF2*((944*(1 - x))/3. + ((30003 + 68912*x - 76790*x2)*H0)/(81.*(1 - x)) + (4*(-6709 + 7333*x)*H1)/27.
                     + (32*(11 + 14*x)*Hm22)/3. - (4*(1 + x)*(8 + x)*(1 + 16*x)*Hm10)/(3.*x) - (4*(-15 + 91*x)*H11)/9.
                     + (32*(-5 + 28*x)*Hm2m10)/3. - 144*Hm2*zeta2
                     + (2*(-15685 - 13731*x + 23474*x2)*H2 - 2*(-12035 + 29454*x - 29056*x2 + 864*x3)*H00
                        - 2*(-704 + 4155*x + 1069*x2)*H10 + 3*(7297 + 7044*x - 13969*x2 + 576*x3)*zeta2)/(81.*(1 - x))
                     + (8*(1 + x)*(44*Hm12 - 134*Hm1m10 + 143*Hm100 - 111*Hm1*zeta2))/3.
                     + (-2*(1382 - 3354*x + 1217*x2)*H3 - 36*(51 - 370*x + 345*x2)*Hm20 + 6*(1079 - 1686*x + 1631*x2)*H12
                        + 2*(1348 - 4554*x + 4387*x2)*H20 + 12*(168 - 137*x + 225*x2)*H21 + 2*(2081 - 3024*x + 3305*x2)*H000
                        + 4*(1837 - 6279*x + 1810*x2)*H100 + 6*(1209 - 1010*x + 873*x2)*H110 - 2*(-373 + 798*x + 2042*x2)*H0*zeta2
                        - 6*(671 - 2023*x + 1304*x2)*H1*zeta2 + 2*(-3424 + 6840*x + 4145*x2)*zeta3)/(27.*(1 - x))
                     + (16*(1 + x)*(14*Hm13 + 8*Hm1m20 - 32*Hm1m12 - 8*Hm120 - 4*Hm121 + 8*Hm1m1m10 - 8*Hm1m100 + 22*Hm1000
                                    + 36*Hm1m1*zeta2 - 13*Hm10*zeta2 - 29*Hm1*zeta3))/3.
                     + (12*(16 - 66*x + 45*x2)*H4 + 48*(-8 - 31*x + 24*x2)*Hm30 - 4*(517 - 1116*x + 469*x2)*H13
                        + 8*(125 - 120*x + 101*x2)*H22 - 16*(-2 - 18*x + 25*x2)*H31 + 48*(31 - 66*x + 41*x2)*Hm200
                        - 144*(15 - 38*x + 21*x2)*H1m20 + 12*(-11 - 54*x + x2)*H112 + 4*(71 + 96*x + 161*x2)*H120
                        - 4*(61 - 12*x + 91*x2)*H121 - 176*(1 + x2)*H211 - 4*(-52 - 162*x + 157*x2)*H0000
                        - 8*(35 + 288*x + 8*x2)*H1000 + 8*(62 - 3*x + 47*x2)*H1100 - 8*(70 + 33*x + 73*x2)*H1110
                        + 12*(59 - 44*x + 37*x2)*H10*zeta2 - 16*(20 - 9*x + 29*x2)*H11*zeta2
                        - 4*(239 - 1008*x + 239*x2)*H1*zeta3)/(9.*(1 - x))
                     + (-4*(-82 - 628*x - 307*x2 + 167*x3)*H30 - 4*(-83 - 389*x + 406*x2 + 568*x3)*H200
                        - 12*(-99 - 139*x - 32*x2 + 56*x3)*H210 - 4*(269 + 35*x - 235*x2 + 143*x3)*H2*zeta2
                        + (-607 + 629*x + 593*x2 - 355*x3)*H00*zeta2 + 4*(43 - 479*x + 355*x2 + 733*x3)*H0*zeta3
                        + (281 + 3257*x + 4106*x2 + 986*x3)*zeta4)/(9.*(1 - x)*(1 + x))
                     + (8*(1 + x)*(16*H5 + 4*H32 - 6*Hm130 - 12*Hm1200 + 12*Hm1210 + 12*Hm12*zeta2
                                   - 6*Hm100*zeta2 + 12*Hm10*zeta3 + 3*Hm1*zeta4))/3.
                     + (-64*(-2 + x2)*Hm40 - 8*(23 + 3*x2)*H23 - 16*(11 + 17*x2)*H40 + 32*(-1 + 2*x2)*H41
                        + 16*(13 + 5*x2)*Hm300 + 32*(13 + 11*x2)*H2m20 - 48*(-1 + 5*x2)*H212 - 8*(29 + 45*x2)*H220
                        - 8*(35 + 43*x2)*H300 - 8*(47 + 27*x2)*H2000 - 352*(1 + 2*x2)*H2100 - 8*(-13 + 35*x2)*H2110
                        - 8*(3 + 7*x2)*(H310 + 6*H00000) - 48*(7 + 5*x2)*H3*zeta2 - 4*(-37 + 43*x2)*H20*zeta2
                        - 16*(-1 + 7*x2)*H21*zeta2 + 18*(-5 + 11*x2)*H000*zeta2 + 8*(107 + 79*x2)*H2*zeta3
                        + 16*(15 + 29*x2)*H00*zeta3 - 6*(159 + 181*x2)*zeta2*zeta3 + 4*(181 + 153*x2)*H0*zeta4
                        + 4*(1 + x2)*(156*Hm32 + 56*Hm23 - 22*H14 + 180*Hm3m10 + 32*Hm2m20 - 128*Hm2m12 - 32*Hm220
                                      - 16*Hm221 + 152*H1m30 + 40*H1m22 - 188*H113 - 8*H122 - 176*H130 - 72*H131
                                      - 14*H221 + 32*Hm2m1m10 - 32*Hm2m100 + 88*Hm2000 - 40*H1m2m10 + 184*H1m200
                                      - 120*H11m20 + 16*H1112 - 72*H1120 + 8*H1121 - 134*H1200 - 42*H1210
                                      - 100*H10000 - 88*H11100 - 66*Hm3*zeta2 + 144*Hm2m1*zeta2 - 52*Hm20*zeta2
                                      - 60*H1m2*zeta2 - 97*H12*zeta2 + 38*H100*zeta2 + 23*H110*zeta2 - 16*H111*zeta2
                                      - 116*Hm2*zeta3 + 78*H10*zeta3 - 24*H11*zeta3 + 111*H1*zeta4)
                        + 8*(325 + 289*x2)*zeta5)/(3.*(1 - x)))
           + CF3*(964*(1 - x) + ((471 - 648*x + 250*x2)*H0)/(3.*(1 - x)) - (2*(-1204 + 1245*x)*H1)/3. + 64*x*Hm30
                  - 64*Hm22 + (8*(1 + x)*(2 + x + 4*x2)*Hm10)/(3.*x) - (8*(-13 + 7*x)*H11)/3. - (64*(-1 + 10*x)*Hm2m10)/3.
                  - (32*(-17 + 19*x)*H1m20)/3. + (4*(-173 + 133*x)*H110)/3. + (8*(-13 + 15*x)*H112)/3. + (32*(-2 + x)*H121)/3.
                  - (8*(-29 + 35*x)*H1100)/3. + (8*(-9 + 11*x)*H1110)/3. - (32*(-7 + 10*x)*Hm2*zeta2)/3.
                  + (-2*(-199 - 624*x + 840*x2)*H2 + 8*(-27 + 93*x - 95*x2 + 4*x3)*H00 - 16*(39 - 67*x + 34*x2)*H10
                     - 2*(520 + 30*x - 535*x2 + 16*x3)*zeta2)/(3.*(1 - x))
                  - (16*(1 + x)*(16*Hm12 - 46*Hm1m10 + 53*Hm100 - 39*Hm1*zeta2))/3.
                  + (8*(29 - 130*x + 89*x2)*H3 + 8*(-3 - 98*x + 91*x2)*Hm20 - 4*(95 - 264*x + 161*x2)*H12
                     - 24*(1 - 24*x + 7*x2)*H20 - 4*(5 - 60*x + 51*x2)*H21 - 4*(36 - 126*x + 119*x2)*H000
                     - 4*(91 - 314*x + 91*x2)*H100 - 2*(131 - 303*x + 44*x2)*H0*zeta2 + 2*(46 - 259*x + 197*x2)*H1*zeta2
                     - 4*(-83 - 176*x + 193*x2)*zeta3)/(3.*(1 - x))
                  - (32*(1 + x)*(2*Hm13 - 4*Hm1m12 + 4*Hm1m1m10 + 6*Hm1000 + 6*Hm1m1*zeta2 - Hm10*zeta2 - 5*Hm1*zeta3))/3.
                  + (-8*(-3 - 26*x + 17*x2)*H4 + 8*(95 - 164*x + 93*x2)*H13 + 24*(-6 - 2*x + 7*x2)*H22 + 8*(3 - 22*x + 25*x2)*H31
                     - 32*(11 - 26*x + 21*x2)*Hm200 + 16*(3 - 26*x + 2*x2)*H120 - 12*(-1 - 10*x + 5*x2)*H0000
                     - 16*(10 - 50*x + 13*x2)*H1000 - 12*(23 - 14*x + 19*x2)*H10*zeta2 - 16*(11 + 28*x + 15*x2)*H1*zeta3)/(3.*(1 - x))
                  + (16*(9 - 73*x - 12*x2 + 58*x3)*H30 + 48*(7 - 23*x + 13*x2 + 35*x3)*H200 + 16*(-25 - 77*x - 23*x2 + 53*x3)*H210
                     + 8*(43 + 17*x - 73*x2 + x3)*H2*zeta2 + 4*(27 - 139*x - 24*x2 + 94*x3)*H00*zeta2
                     - 16*(61 + x + 2*x2 + 38*x3)*H0*zeta3 + (-1573 - 1987*x - 143*x2 + 367*x3)*zeta4)/(6.*(1 - x)*(1 + x))
                  + 8*(1 + x)*(2*Hm130 + 4*Hm1200 - 4*Hm1210 - 4*Hm12*zeta2 + 2*Hm100*zeta2 - 4*Hm10*zeta3 - Hm1*zeta4)
                  + (-16*(5 + 7*x2)*H5 + 256*x2*Hm40 + 16*(23 + 17*x2)*H23 + 32*(5 + 9*x2)*H32 + 128*(1 + 2*x2)*H40
                     + 32*(1 + 5*x2)*Hm300 + 32*(7 + 13*x2)*H212 + 32*(15 + 23*x2)*H220 + 32*(9 + 10*x2)*H300 + 96*(2 + 3*x2)*H310
                     - 96*(-1 + 2*x2)*H2000 + 32*(9 + 19*x2)*H2100 + 8*(-29 + 7*x2)*H20*zeta2 + 4*(-23 + x2)*H21*zeta2
                     + 2*(29 + 35*x2)*H000*zeta2 + 16*(1 + 3*x2)*(H41 + 9*H2110 - 2*H3*zeta2) + 16*(-5 + 19*x2)*H2*zeta3
                     + 48*(-1 + 3*x)*(1 + 3*x)*H00*zeta3 + 8*(73 + 85*x2)*zeta2*zeta3 - 4*(139 + 17*x2)*H0*zeta4
                     - 4*(1 + x2)*(56*Hm32 + 16*Hm23 + 168*Hm3m10 - 32*Hm2m12 + 16*H1m30 + 48*H1m22 - 192*H113 - 84*H122
                                   - 152*H130 - 104*H131 + 64*H2m20 - 64*H221 - 48*H311 + 32*Hm2m1m10 + 48*Hm2000 + 16*H1m2m10
                                   + 112*H1m200 - 80*H11m20 - 120*H1120 - 136*H1200 - 104*H1210 - 24*H00000 + 112*H11000
                                   - 96*H11100 + 28*Hm3*zeta2 + 48*Hm2m1*zeta2 - 8*Hm20*zeta2 - 40*H1m2*zeta2 - 30*H12*zeta2
                                   - 3*H100*zeta2 + 38*H110*zeta2 - 40*Hm2*zeta3 - 136*H10*zeta3 - 72*H11*zeta3 + 18*H1*zeta4)
                     - 8*(381 + 355*x2)*zeta5)/(3.*(1 - x)));
  }
  double C3Vqqpdf::Singular(double const& x) const
  {
    return _A2 / ( 1 - x );
  }
  double C3Vqqpdf::Local(double const& x) const
  {
    const double CF2 = CF * CF;
    const double CF3 = CF * CF2;
    const double CA2 = CA * CA;
    const double TR2 = TR * TR;
    const double nf2 = _nf * _nf;
    const double Pi6 = Pi2 * Pi2 * Pi2;
    const double ln1mx = log( 1 - x );
    const double A1    = -(CF3*Pi6)/1296.
                         + CF*nf2*TR2*(-0.07803688462124676 - (272*zeta2)/27. - (1120*zeta3)/243. - (88*zeta4)/27.)
                         + CA*CF2*(-Pi6/108. - (1214*zeta2)/81. + (77*zeta2*zeta3)/9. + (335*zeta4)/12.)
                         + CA*CF*_nf*TR*(-62.91190367322054 + (74530*zeta2)/729. + (8152*zeta3)/81.
                                         + (40*zeta2*zeta3)/9. - (416*zeta4)/27. - (184*zeta5)/3.)
                         + CF2*_nf*TR*(-87.9156378600823 + (2803*zeta2)/81. + (3488*zeta3)/81.
                                       - (268*zeta2*zeta3)/9. + (77*zeta4)/9. + (224*zeta5)/9.)
                         + CA2*CF*(198.59583142813597 - (1543*Pi6)/25515. - (297481*zeta2)/1458.
                                   - (75566*zeta3)/243. + (550*zeta2*zeta3)/9. + (464*pow(zeta3,2))/9.
                                   + (3649*zeta4)/54. + (902*zeta5)/9.);// + (3./4.) * _A2;

    return A1 + _A2 * ln1mx;
  }

  //_________________________________________________________________________________
  C3Vqqbpdf::C3Vqqbpdf(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C3Vqqbpdf::Regular(double const& x) const
  {
    // Polylogs
    double *xx  = new double{x};
    int    *nw  = new int{5};
    int    *wn1 = new int{-1};
    int    *wn2 = new int{1};
    double *H   = new double[363];
    apf_hplog_(xx, nw, &H[0], &H[3], &H[12], &H[39], &H[120], wn1, wn2);

    // Defintions
    const double x2  = x * x;
    const double x3  = x * x2;
    const double x4  = x * x3;

    const double Hm1 = H[HPLogMap({-1})];
    const double H0  = H[HPLogMap({0})];
    const double H1  = H[HPLogMap({1})];

    const double Hm10  = H[HPLogMap({-1,0})];
    const double H00   = H[HPLogMap({0,0})];
    const double H2    = H[HPLogMap({2})];
    const double H10   = H[HPLogMap({1,0})];
    const double H11   = H[HPLogMap({1,1})];
    const double Hm2   = H[HPLogMap({-2})];
    const double Hm1m1 = H[HPLogMap({-1,-1})];

    const double Hm20    = H[HPLogMap({-2,0})];
    const double H000    = H[HPLogMap({0,0,0})];
    const double H12     = H[HPLogMap({1,2})];
    const double H20     = H[HPLogMap({2,0})];
    const double H21     = H[HPLogMap({2,1})];
    const double H100    = H[HPLogMap({1,0,0})];
    const double H110    = H[HPLogMap({1,1,0})];
    const double Hm1m10  = H[HPLogMap({-1,-1,0})];
    const double Hm100   = H[HPLogMap({-1,0,0})];
    const double H3      = H[HPLogMap({3})];
    const double Hm12    = H[HPLogMap({-1,2})];
    const double Hm3     = H[HPLogMap({-3})];
    const double Hm2m1   = H[HPLogMap({-2,-1})];
    const double Hm1m2   = H[HPLogMap({-1,-2})];
    const double Hm1m1m1 = H[HPLogMap({-1,-1,-1})];

    const double H4       = H[HPLogMap({4})];
    const double H22      = H[HPLogMap({2,2})];
    const double H30      = H[HPLogMap({3,0})];
    const double H210     = H[HPLogMap({2,1,0})];
    const double H13      = H[HPLogMap({1,3})];
    const double H31      = H[HPLogMap({3,1})];
    const double H112     = H[HPLogMap({1,1,2})];
    const double H120     = H[HPLogMap({1,2,0})];
    const double H200     = H[HPLogMap({2,0,0})];
    const double H0000    = H[HPLogMap({0,0,0,0})];
    const double H1000    = H[HPLogMap({1,0,0,0})];
    const double H1100    = H[HPLogMap({1,1,0,0})];
    const double H1110    = H[HPLogMap({1,1,1,0})];
    const double Hm22     = H[HPLogMap({-2,2})];
    const double Hm2m10   = H[HPLogMap({-2,-1,0})];
    const double Hm200    = H[HPLogMap({-2,0,0})];
    const double Hm13     = H[HPLogMap({-1,3})];
    const double Hm1m20   = H[HPLogMap({-1,-2,0})];
    const double Hm1m12   = H[HPLogMap({-1,-1,2})];
    const double Hm120    = H[HPLogMap({-1,2,0})];
    const double Hm121    = H[HPLogMap({-1,2,1})];
    const double Hm1m1m10 = H[HPLogMap({-1,-1,-1,0})];
    const double Hm1m100  = H[HPLogMap({-1,-1,0,0})];
    const double Hm1000   = H[HPLogMap({-1,0,0,0})];
    const double Hm30     = H[HPLogMap({-3,0})];
    const double H1m20    = H[HPLogMap({1,-2,0})];
    const double Hm130    = H[HPLogMap({-1,3,0})];
    const double Hm1200   = H[HPLogMap({-1,2,0,0})];
    const double Hm1210   = H[HPLogMap({-1,2,1,0})];

    const double H5        = H[HPLogMap({5})];
    const double Hm40      = H[HPLogMap({-4,0})];
    const double H23       = H[HPLogMap({2,3})];
    const double H41       = H[HPLogMap({4,1})];
    const double Hm300     = H[HPLogMap({-3,0,0})];
    const double H2m20     = H[HPLogMap({2,-2,0})];
    const double H220      = H[HPLogMap({2,2,0})];
    const double H300      = H[HPLogMap({3,0,0})];
    const double H2000     = H[HPLogMap({2,0,0,0})];
    const double H2100     = H[HPLogMap({2,1,0,0})];
    const double H2110     = H[HPLogMap({2,1,1,0})];
    const double H00000    = H[HPLogMap({0,0,0,0,0})];
    const double H212      = H[HPLogMap({2,1,2})];
    const double Hm32      = H[HPLogMap({-3,2})];
    const double Hm23      = H[HPLogMap({-2,3})];
    const double H32       = H[HPLogMap({3,2})];
    const double H40       = H[HPLogMap({4,0})];
    const double Hm3m10    = H[HPLogMap({-3,-1,0})];
    const double Hm2m20    = H[HPLogMap({-2,-2,0})];
    const double Hm2m12    = H[HPLogMap({-2,-1,2})];
    const double Hm220     = H[HPLogMap({-2,2,0})];
    const double Hm221     = H[HPLogMap({-2,2,1})];
    const double H310      = H[HPLogMap({3,1,0})];
    const double Hm2m1m10  = H[HPLogMap({-2,-1,-1,0})];
    const double Hm2m100   = H[HPLogMap({-2,-1,0,0})];
    const double Hm2000    = H[HPLogMap({-2,0,0,0})];
    const double Hm14      = H[HPLogMap({-1,4})];
    const double Hm1m30    = H[HPLogMap({-1,-3,0})];
    const double Hm1m22    = H[HPLogMap({-1,-2,2})];
    const double Hm1m13    = H[HPLogMap({-1,-1,3})];
    const double Hm122     = H[HPLogMap({-1,2,2})];
    const double Hm131     = H[HPLogMap({-1,3,1})];
    const double Hm1m2m10  = H[HPLogMap({-1,-2,-1,0})];
    const double Hm1m200   = H[HPLogMap({-1,-2,0,0})];
    const double Hm1m1m20  = H[HPLogMap({-1,-1,-2,0})];
    const double Hm1m1m12  = H[HPLogMap({-1,-1,-1,2})];
    const double Hm1m120   = H[HPLogMap({-1,-1,2,0})];
    const double Hm1m121   = H[HPLogMap({-1,-1,2,1})];
    const double Hm1m1m100 = H[HPLogMap({-1,-1,-1,0,0})];
    const double Hm1m1000  = H[HPLogMap({-1,-1,0,0,0})];
    const double Hm10000   = H[HPLogMap({-1,0,0,0,0})];

    // Delete pointers
    delete xx;
    delete nw;
    delete wn1;
    delete wn2;
    delete[] H;

    return
      (CA - 2*CF)*CF*(_nf*TR*((392*(1 - x))/81. + (32*(11 + 47*x)*H0)/81. - (592*(1 - x)*H1)/27.
                              - (16*(-7 + 29*x)*H2)/27. - (64*(1 - x)*(H10 + H11))/9.
                              + (448*(1 - 6*x + x2)*Hm10 + 8*(-193 + 54*x + 23*x2)*H00
                                 + 16*(-43 - 18*x + 137*x2)*zeta2)/(81.*(1 + x))
                              + (16*(19 + 18*x + 19*x2)*H3 + 32*(35 + 18*x + 53*x2)*Hm20
                                 - 32*(7 - 6*x + 7*x2)*Hm12 - 48*(13 + 6*x + 13*x2)*H20
                                 - 192*(7 + 4*x + 7*x2)*Hm1m10 + 16*(49 + 18*x + 49*x2)*Hm100
                                 - 16*(43 + 36*x + 43*x2)*H000 - 64*(7 + 9*x + 7*x2)*Hm1*zeta2
                                 - 16*(23 + 18*x + 5*x2)*H0*zeta2)/(27.*(1 + x))
                              + (16*(1 - x)*(2*H12 - 3*H100 + 4*H110 + 2*H1*zeta2))/3.
                              - (32*(1 + x)*(H21 - zeta3))/9.
                              + (16*(1 + x2)*(5*H4 + 20*Hm30 - 14*Hm22 - 8*Hm13 + 6*H22 - 9*H30
                                              - 2*H31 - 20*Hm2m10 + 17*Hm200 - 20*Hm1m20 + 12*Hm1m12
                                              + 6*Hm120 + 4*Hm121 - 9*H200 + 12*H210 + 16*Hm1m1m10
                                              - 14*Hm1m100 + 11*Hm1000 - 8*H0000 + 4*Hm2*zeta2
                                              + 6*H2*zeta2 - 4*Hm1m1*zeta2 + 4*Hm10*zeta2
                                              - 4*H00*zeta2 + 2*Hm1*zeta3 - H0*zeta3 + 6*zeta4))/(9.*(1 + x)))
                      + CA*((22292*(1 - x))/81. - (4*(-859 + 2287*x + 3065*x2)*H0)/(81.*(1 + x))
                            + (10832*(1 - x)*H1)/27. + (4*(439 + 1744*x)*H2)/27.
                            - (4*(81 + 1457*x - 723*x2 + 1619*x3 + 162*x4)*Hm10)/(81.*x*(1 + x))
                            - (8*(1 - x)*(17*H10 - 70*H11))/9. + (280*(1 + x)*H21)/9. + (16*(1 - x)*Hm2m1m10)/3.
                            + (2*(2183 + 4590*x + 6287*x2 + 324*x3)*H00
                               - 2*(6701 + 12294*x + 9473*x2 + 324*x3)*zeta2)/(81.*(1 + x))
                            - (4*(1 - x)*(86*H12 + 41*H100 + 67*H110 - 47*H1*zeta2))/3.
                            + (4*(-296 + 81*x + 73*x2)*H3 - 4*(983 + 486*x + 1415*x2)*Hm20
                               - 8*(259 + 822*x + 259*x2)*Hm12 + 12*(116 + 33*x + 185*x2)*H20
                               + 24*(229 + 166*x + 205*x2)*Hm1m10 - 8*(709 + 864*x + 709*x2)*Hm100
                               + 4*(470 + 972*x + 1307*x2)*H000 + 4*(1205 + 2142*x + 1133*x2)*Hm1*zeta2
                               - 4*(-22 + 738*x + 689*x2)*H0*zeta2 - 6*(515 + 592*x + 41*x2)*zeta3)/(27.*(1 + x))
                            - (16*(1 - x)*(11*H13 + 8*H1m20 - 13*H112 - 13*H120 + 11*H1000 - H1100 - 20*H1110
                                           - 20*H10*zeta2 - 9*H11*zeta2 - 22*H1*zeta3))/3.
                            + (-4*(103 + 186*x + 193*x2)*H4 - 8*(59 - 39*x + 95*x2)*Hm30
                               + 8*(101 + 261*x + 314*x2)*Hm22 + 64*(10 + 9*x + 10*x2)*Hm13
                               - 24*(27 + 58*x + 53*x2)*H22 - 12*(-15 + 62*x + 11*x2)*H30
                               + 8*(23 + 48*x + 47*x2)*H31 + 8*(161 + 9*x + 95*x2)*Hm2m10
                               - 4*(184 - 168*x + 31*x2)*Hm200 + 8*(98 + 39*x + 134*x2)*Hm1m20
                               - 48*(23 + 24*x + 23*x2)*Hm1m12 - 24*(6 + 8*x + 15*x2)*Hm120
                               - 16*(35 + 48*x + 35*x2)*Hm121 - 12*(11 + 68*x + 15*x2)*H200
                               - 12*(51 + 112*x + 131*x2)*H210 - 8*(76 - 21*x + 79*x2)*Hm1m1m10
                               + 44*(14 + 3*x + 17*x2)*Hm1m100 - 8*(50 - 30*x + 59*x2)*Hm1000
                               + 4*(64 - 72*x + 67*x2)*H0000 - 4*(41 + 513*x + 533*x2)*Hm2*zeta2
                               - 12*(-3 - 28*x + x2)*H2*zeta2 + 4*(200 + 309*x + 197*x2)*Hm1m1*zeta2
                               - 4*(62 + 15*x + 32*x2)*Hm10*zeta2 + 4*(65 + 129*x + 134*x2)*H00*zeta2
                               - 4*(67 + 141*x + 91*x2)*Hm1*zeta3 + 4*(-19 + 21*x + 35*x2)*H0*zeta3
                               + 3*(279 + 249*x + 145*x2)*zeta4)/(9.*(1 + x))
                            + (-4*(25 + 23*x2)*H5 + 8*(7 + 13*x2)*Hm40 + 8*(73 + 71*x2)*Hm32
                               + 4*(145 + 147*x2)*Hm23 - 8*(17 + 15*x2)*H32 + 16*(15 + 11*x2)*Hm3m10
                               + 4*(41 + 59*x2)*Hm300 + 8*(11 + 9*x2)*Hm2m20 - 8*(151 + 149*x2)*Hm2m12
                               - 8*(7 + 9*x2)*Hm220 - 4*(35 + 29*x2)*H300 - 4*(111 + 133*x2)*Hm2m100
                               + 8*(25 + 33*x2)*Hm2000 - 12*(3 + 5*x2)*H00000 - 16*(29 + 30*x2)*Hm3*zeta2
                               + 32*(38 + 37*x2)*Hm2m1*zeta2 - 4*(11 + 21*x2)*(H310 + H3*zeta2)
                               - 8*(107 + 105*x2)*Hm2*zeta3 + 4*(11 + 37*x2)*H00*zeta3
                               + 16*(5 + 8*x2)*zeta2*zeta3 + (231 + 197*x2)*H0*zeta4
                               + 8*(1 + x2)*(36*Hm14 - 22*H23 - 6*H40 + 4*H41 - 16*Hm221 - 20*Hm1m30
                                             - 150*Hm1m22 - 84*Hm1m13 + 58*Hm122 + 31*Hm130 - 16*Hm131
                                             - 16*H2m20 + 26*H212 + 26*H220 - 24*Hm1m2m10 - 57*Hm1m200
                                             - 18*Hm1m1m20 + 168*Hm1m1m12 + 14*Hm1m120 + 32*Hm1m121
                                             + 34*Hm1200 + 56*Hm1210 - 22*H2000 + 2*H2100 + 40*H2110
                                             + 60*Hm1m1m100 - 38*Hm1m1000 + 16*Hm10000 - 94*Hm20*zeta2
                                             + 138*Hm1m2*zeta2 + 2*Hm12*zeta2 + 40*H20*zeta2 + 18*H21*zeta2
                                             - 168*Hm1m1m1*zeta2 + 95*Hm1m10*zeta2 - 34*Hm100*zeta2
                                             + 14*H000*zeta2 + 44*H2*zeta3 + 113*Hm1m1*zeta3 - 22*Hm10*zeta3
                                             - 34*Hm1*zeta4) + 8*(30 + 11*x2)*zeta5)/(3.*(1 + x)))
                      + CF*((1294*(1 - x))/3. + (4*(62 + 703*x + 638*x2)*H0)/(3.*(1 + x)) - (1532*(1 - x)*H1)/3.
                            - (4*(103 + 205*x)*H2)/3. + (8*(1 - 29*x - 65*x2 - 27*x3 + 2*x4)*Hm10)/(3.*x*(1 + x))
                            - 80*(1 - x)*(H10 + H11) - (4*(17 + 56*x)*H20)/3. - (8*(19 + 11*x)*H21)/3.
                            - (16*(1 - x)*Hm220)/3. - (8*(51 + 43*x)*Hm1m10)/3. + (176*(1 + x)*Hm120)/3.
                            + (16*(25 + 13*x)*H200)/3. + (8*(11 + 9*x)*Hm1m1m10)/3.
                            + (-2*(19 + 670*x + 671*x2 + 8*x3)*H00 + (151 + 968*x + 857*x2 + 16*x3)*zeta2)/(3.*(1 + x))
                            + (4*(1 - x)*(101*H12 + 174*H100 + 15*H110 - 129*H1*zeta2))/3.
                            + (-4*(-30 + 87*x + 85*x2)*H3 + 4*(27 + 162*x + 167*x2)*Hm20 + 8*(49 + 130*x + 49*x2)*Hm12
                               + 16*(45 + 82*x + 45*x2)*Hm100 - 4*(33 + 177*x + 187*x2)*H000
                               - 4*(149 + 354*x + 141*x2)*Hm1*zeta2 + (-111 + 646*x + 673*x2)*H0*zeta2
                               + 2*(-123 + 368*x + 395*x2)*zeta3)/(3.*(1 + x))
                            + (4*(1 - x)*(96*H13 + 8*H1m20 - 76*H112 - 52*H120 + 84*H1000 + 24*H1100 - 128*H1110
                                          - 155*H10*zeta2 - 64*H11*zeta2 - 188*H1*zeta3))/3.
                            + (8*(28 + 50*x + 43*x2)*H4 - 8*(-18 + 22*x + 13*x2)*Hm30 - 8*(28 + 118*x + 135*x2)*Hm22
                               - 16*(11 + 10*x + 11*x2)*Hm13 + 4*(73 + 156*x + 125*x2)*H22 + 4*(-12 + 44*x + 29*x2)*H30
                               - 8*(13 + 32*x + 25*x2)*H31 + 8*(-47 - 16*x + x2)*Hm2m10 - 4*(-61 + 98*x + 108*x2)*Hm200
                               - 8*(37 + 56*x + 37*x2)*Hm1m20 + 144*(3 + 4*x + 3*x2)*Hm1m12 + 48*(5 + 8*x + 5*x2)*Hm121
                               + 4*(27 + 128*x + 131*x2)*H210 - 4*(35 + 36*x + 37*x2)*Hm1m100 + 8*(10 - 4*x + 13*x2)*Hm1000
                               + 12*(-5 + 26*x + 16*x2)*H0000 + 4*(9 + 220*x + 271*x2)*Hm2*zeta2
                               - 4*(47 + 68*x + 27*x2)*H2*zeta2 - 4*(97 + 124*x + 99*x2)*Hm1m1*zeta2
                               + 8*(15 + 13*x + 13*x2)*Hm10*zeta2 - 4*(44 + 104*x + 93*x2)*H00*zeta2
                               + 8*(26 + 32*x + 27*x2)*Hm1*zeta3 - 4*(49 + 140*x + 109*x2)*H0*zeta3
                               + (-1307 - 1050*x - 97*x2)*zeta4)/(3.*(1 + x))
                            + (8*(31 + 29*x2)*H5 - 16*(7 + 29*x2)*Hm40 - 32*(32 + 31*x2)*Hm32 - 16*(68 + 69*x2)*Hm23
                               + 16*(17 + 15*x2)*H32 + 16*(5 + 11*x2)*H40 - 16*(38 + 59*x2)*Hm300 + 96*(1 + 2*x2)*Hm2m20
                               + 32*(59 + 58*x2)*Hm2m12 + 56*(7 + 9*x2)*H300 + 8*(13 + 11*x2)*H310 + 272*(4 + 5*x2)*Hm2m100
                               - 24*(25 + 33*x2)*Hm2000 + 8*(1 + 3*x2)*(8*Hm3m10 + 9*H00000) + 32*(33 + 34*x2)*Hm3*zeta2
                               + 16*H3*zeta2 - 16*(121 + 119*x2)*Hm2m1*zeta2 + 8*(171 + 179*x2)*Hm20*zeta2
                               - 8*(31 + 57*x2)*H00*zeta3 + 2*(119 + 127*x2)*zeta2*zeta3 - 2*(267 + 353*x2)*H0*zeta4
                               - 2*(1 + x2)*(264*Hm14 - 192*H23 + 48*H41 - 112*Hm221 - 224*Hm1m30 - 936*Hm1m22
                                             - 592*Hm1m13 + 312*Hm122 + 148*Hm130 - 128*Hm131 - 16*H2m20 + 152*H212
                                             + 104*H220 + 48*Hm2m1m10 - 620*Hm1m200 + 24*Hm1m1m20 + 1056*Hm1m1m12
                                             + 40*Hm1m120 + 192*Hm1m121 + 304*Hm1200 + 256*Hm1210 - 168*H2000
                                             - 48*H2100 + 256*H2110 + 624*Hm1m1m100 - 416*Hm1m1000 + 216*Hm10000
                                             + 936*Hm1m2*zeta2 - 40*Hm12*zeta2 + 310*H20*zeta2 + 128*H21*zeta2
                                             - 1056*Hm1m1m1*zeta2 + 744*Hm1m10*zeta2 - 294*Hm100*zeta2 + 131*H000*zeta2
                                             - 676*Hm2*zeta3 + 376*H2*zeta3 + 708*Hm1m1*zeta3 - 360*Hm10*zeta3 - 233*Hm1*zeta4)
                               - 4*(231 + 313*x2)*zeta5)/(3.*(1 + x))));
  }

  //_________________________________________________________________________________
  C3pspdf::C3pspdf(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C3pspdf::Regular(double const& x) const
  {
    // Polylogs
    double *xx  = new double{x};
    int    *nw  = new int{5};
    int    *wn1 = new int{-1};
    int    *wn2 = new int{1};
    double *H   = new double[363];
    apf_hplog_(xx, nw, &H[0], &H[3], &H[12], &H[39], &H[120], wn1, wn2);

    // Defintions
    const double x2  = x * x;
    const double x3  = x * x2;
    const double x4  = x * x3;
    const double CF2 = CF * CF;
    const double TR2 = TR * TR;

    const double Hm1 = H[HPLogMap({-1})];
    const double H0  = H[HPLogMap({0})];
    const double H1  = H[HPLogMap({1})];

    const double Hm10  = H[HPLogMap({-1,0})];
    const double H00   = H[HPLogMap({0,0})];
    const double H2    = H[HPLogMap({2})];
    const double H10   = H[HPLogMap({1,0})];
    const double H11   = H[HPLogMap({1,1})];
    const double Hm2   = H[HPLogMap({-2})];
    const double Hm1m1 = H[HPLogMap({-1,-1})];

    const double Hm20    = H[HPLogMap({-2,0})];
    const double H000    = H[HPLogMap({0,0,0})];
    const double H12     = H[HPLogMap({1,2})];
    const double H20     = H[HPLogMap({2,0})];
    const double H21     = H[HPLogMap({2,1})];
    const double H100    = H[HPLogMap({1,0,0})];
    const double H110    = H[HPLogMap({1,1,0})];
    const double Hm1m10  = H[HPLogMap({-1,-1,0})];
    const double Hm100   = H[HPLogMap({-1,0,0})];
    const double H3      = H[HPLogMap({3})];
    const double Hm12    = H[HPLogMap({-1,2})];
    const double Hm3     = H[HPLogMap({-3})];
    const double H111    = H[HPLogMap({1,1,1})];

    const double H4       = H[HPLogMap({4})];
    const double H22      = H[HPLogMap({2,2})];
    const double H30      = H[HPLogMap({3,0})];
    const double H210     = H[HPLogMap({2,1,0})];
    const double H13      = H[HPLogMap({1,3})];
    const double H31      = H[HPLogMap({3,1})];
    const double H112     = H[HPLogMap({1,1,2})];
    const double H120     = H[HPLogMap({1,2,0})];
    const double H200     = H[HPLogMap({2,0,0})];
    const double H0000    = H[HPLogMap({0,0,0,0})];
    const double H1000    = H[HPLogMap({1,0,0,0})];
    const double H1100    = H[HPLogMap({1,1,0,0})];
    const double H1110    = H[HPLogMap({1,1,1,0})];
    const double Hm22     = H[HPLogMap({-2,2})];
    const double Hm2m10   = H[HPLogMap({-2,-1,0})];
    const double Hm200    = H[HPLogMap({-2,0,0})];
    const double Hm13     = H[HPLogMap({-1,3})];
    const double Hm1m20   = H[HPLogMap({-1,-2,0})];
    const double Hm1m12   = H[HPLogMap({-1,-1,2})];
    const double Hm120    = H[HPLogMap({-1,2,0})];
    const double Hm121    = H[HPLogMap({-1,2,1})];
    const double Hm1m1m10 = H[HPLogMap({-1,-1,-1,0})];
    const double Hm1m100  = H[HPLogMap({-1,-1,0,0})];
    const double Hm1000   = H[HPLogMap({-1,0,0,0})];
    const double Hm30     = H[HPLogMap({-3,0})];
    const double H1m20    = H[HPLogMap({1,-2,0})];
    const double Hm130    = H[HPLogMap({-1,3,0})];
    const double Hm1200   = H[HPLogMap({-1,2,0,0})];
    const double Hm1210   = H[HPLogMap({-1,2,1,0})];
    const double H211     = H[HPLogMap({2,1,1})];
    const double H1111    = H[HPLogMap({1,1,1,1})];

    const double H5        = H[HPLogMap({5})];
    const double Hm40      = H[HPLogMap({-4,0})];
    const double H23       = H[HPLogMap({2,3})];
    const double H41       = H[HPLogMap({4,1})];
    const double Hm300     = H[HPLogMap({-3,0,0})];
    const double H2m20     = H[HPLogMap({2,-2,0})];
    const double H220      = H[HPLogMap({2,2,0})];
    const double H300      = H[HPLogMap({3,0,0})];
    const double H2000     = H[HPLogMap({2,0,0,0})];
    const double H2100     = H[HPLogMap({2,1,0,0})];
    const double H2110     = H[HPLogMap({2,1,1,0})];
    const double H00000    = H[HPLogMap({0,0,0,0,0})];
    const double H212      = H[HPLogMap({2,1,2})];
    const double Hm32      = H[HPLogMap({-3,2})];
    const double Hm23      = H[HPLogMap({-2,3})];
    const double H32       = H[HPLogMap({3,2})];
    const double H40       = H[HPLogMap({4,0})];
    const double Hm3m10    = H[HPLogMap({-3,-1,0})];
    const double Hm2m20    = H[HPLogMap({-2,-2,0})];
    const double Hm2m12    = H[HPLogMap({-2,-1,2})];
    const double Hm220     = H[HPLogMap({-2,2,0})];
    const double H310      = H[HPLogMap({3,1,0})];
    const double Hm2m1m10  = H[HPLogMap({-2,-1,-1,0})];
    const double Hm2m100   = H[HPLogMap({-2,-1,0,0})];
    const double Hm2000    = H[HPLogMap({-2,0,0,0})];
    const double H311      = H[HPLogMap({3,1,1})];
    const double H2111     = H[HPLogMap({2,1,1,1})];

    // Delete pointers
    delete xx;
    delete nw;
    delete wn1;
    delete wn2;
    delete[] H;

    return 2*CF*_nf*TR2*((32*(1 - x)*(1058 - 2689*x + 1571*x2))/(729.*x) + (16*(-661 - 337*x + 534*x2)*H0)/243.
                         + (32*(1 - x)*(75 + 19*x + 75*x2)*H1)/(81.*x) + (32*(100 - 26*x + 57*x2)*H2)/81.
                         - (32*(137 + 38*x + 18*x2)*H00)/81.
                         + ((1 - x)*(-576*(1 - x + x2)*H10 + 16*(38 - 79*x + 38*x2)*H11))/(81.*x)
                         + (16*(-49 - 49*x + 24*x2)*H000)/27. - (32*(18 + 64*x + 10*x2 + 39*x3)*zeta2)/(81.*x)
                         - (32*(-5 - 5*x + 6*x2)*(2*H3 - H21 - 2*H0*zeta2))/27.
                         + ((1 - x)*(16*(4 + 7*x + 4*x2)*H111 - 192*(2 - x + 2*x2)*(H12 - H100 + 2*H110 + H1*zeta2)))/(27.*x)
                         + (32*(-12 + 13*x - 23*x2 + 18*x3)*zeta3)/(27.*x)
                         + (16*(1 + x)*(8*H4 - 4*H31 + 2*H211 - 10*H0000 - 8*H00*zeta2 - 4*H0*zeta3 - 9*zeta4))/9.)
           + 2*CF2*TR*((-2*(1 - x)*(2421 + 56*x + 7512*x2))/(81.*x) + (2*(2707 - 8308*x + 5720*x2)*H0)/81.
                       - (2*(1 - x)*(4343 - 5805*x + 7088*x2)*H1)/(81.*x) - (4*(-1524 + 75*x + 2882*x2)*H2)/81.
                       + (8*(229 - 55*x + 104*x2)*H3)/9. - (16*(-24 - 27*x + 40*x2)*H4)/9. - (2*(2733 + 9273*x + 5864*x2)*H00)/81.
                       + ((1 - x)*(8*(101 + 161*x + 362*x2)*H10 - 4*(361 + 2557*x + 406*x2)*H11))/(81.*x)
                       + (16*(-15 + 123*x + 68*x2)*H20)/27. - (4*(612 + 423*x + 184*x2)*H21)/27. + (16*(11 + 10*x)*H22)/3.
                       - (16*(12 + 45*x + 4*x2)*H30)/9. + (16*(-6 + 9*x + 16*x2)*H31)/9. + (4*(257 - 167*x + 376*x2)*H000)/9.
                       + (16*(-2 + x)*(-5 + 4*x)*H200)/3. - (16*(7 + 2*x + 4*x2)*H210)/3. + (16*(9 + 21*x + 4*x2)*H211)/9.
                       - (8*(66 - 15*x + 140*x2)*H0000)/9. + (2*(-1812 - 1815*x + 537*x2 + 5956*x3)*zeta2)/(81.*x)
                       - (4*(617 - 134*x + 272*x2)*H0*zeta2)/9. - (16*(14 + 7*x + 4*x2)*H2*zeta2)/3.
                       + (2*(-111 - 279*x + 440*x2)*H00*zeta2)/9.
                       + ((1 - x)*(16*(20 + 203*x + 20*x2)*H12 - 16*(-1 - 160*x + 26*x2)*H100 - 8*(8 + 209*x + 8*x2)*H110
                                   + 12*(4 + 47*x + 40*x2)*H111 + 8*(8 - 493*x + 8*x2)*H1*zeta2))/(27.*x)
                       - (4*(-768 - 720*x - 999*x2 + 1528*x3)*zeta3)/(27.*x) + (32*(42 + 81*x + 38*x2)*H0*zeta3)/9.
                       + ((1 - x)*(-16*(20 - 19*x + 20*x2)*H13 + 64*(5 + 2*x + 5*x2)*H112 - 8*(52 - 17*x + 52*x2)*H1000
                                   + 32*(2 + 17*x + 2*x2)*H1100 + 16*(32 - 25*x + 32*x2)*H1110 - 16*(4 + 7*x + 4*x2)*H1111
                                   + 8*(70 - 53*x + 70*x2)*H10*zeta2 + 16*(8 - 13*x + 8*x2)*(H120 + 2*H11*zeta2)
                                   + 96*(4 - 11*x + 4*x2)*H1*zeta3))/(9.*x) - (4*(-78 + 531*x + 66*x2 + 82*x3)*zeta4)/(9.*x)
                       + (4*(1 + x)*(56*H5 + 8*H23 + 16*H40 - 24*H41 + 16*H212 - 8*H220 + 16*H300 - 8*H310 - 4*H311
                                     - 4*H2000 + 32*H2100 - 8*H2110 - 8*H2111 + 48*H00000 - 8*H20*zeta2 - 16*H21*zeta2
                                     - 59*H000*zeta2 - 48*H2*zeta3 - 72*H00*zeta3 + 28*zeta2*zeta3 - 56*H0*zeta4 - 144*zeta5))/3.)
           + 2*CA*CF*TR*((-4*(1 - x)*(228133 - 36824*x + 356788*x2))/(729.*x)
                         - (2*(17152 + 143678*x + 34691*x2 + 214384*x3)*H0)/(243.*x)
                         + (2*(1 - x)*(14944 - 3245*x + 15520*x2)*H1)/(243.*x) + (4*(-724 - 451*x + 224*x2)*H3)/27.
                         + (4*(-109 - 91*x + 8*x2)*H4)/9. + (16*(-25 + 13*x)*Hm40)/3. + (8*(87 - 99*x + 112*x2)*Hm30)/9.
                         - (8*(576 + 93*x + 544*x2)*Hm20)/27. + (8*(27 + x + 16*x2)*Hm22)/3.
                         + (4*(1 + x)*(906 - 77*x + 828*x2)*Hm10)/(27.*x) + (4*(14795 + 3176*x + 16868*x2)*H00)/81.
                         + ((1 - x)*(4*(4481 + 272*x + 6272*x2)*H10 + 4*(463 + 2098*x + 184*x2)*H11))/(81.*x)
                         + (4*(59 + 35*x + 8*x2)*H31)/9. - (16*(5 + 3*x)*H32)/3. - 64*x*H40 + (8*(-37 + 29*x)*Hm300)/3.
                         - (64*(-2 + x)*x*Hm2m10)/3. + (4*(159 - 315*x + 184*x2)*Hm200)/9. + 8*(1 + x)*Hm120
                         - (4*(1261 - 587*x + 1758*x2)*H000)/27.
                         + (48*(-4 - 6*x - 7*x2 + 4*x3)*H22 - 8*(-4 + 11*x + 35*x2 + 4*x3)*H211)/(9.*x)
                         - (16*(-6 + 11*x)*H300)/3. + (8*(-11 + 3*x)*H310)/3. + (8*(121 + 25*x)*H0000)/9.
                         + 16*(-4 + 3*x)*H00000 - (8*(27 - 7*x + 20*x2)*Hm2*zeta2)/3. - 20*(1 - x)*H10*zeta2
                         + (4*(1256 - 2068*x + 2969*x2 + 1000*x3)*H2 - 8*(-3002 + 1099*x - 2936*x2 + 3932*x3)*zeta2)/(81.*x)
                         + ((1 + x)*(-8*(116 - 101*x + 116*x2)*Hm12 + 16*(20 + 31*x + 74*x2)*Hm1m10
                                     - 4*(484 - 25*x + 592*x2)*Hm100 + 16*(68 - 35*x + 95*x2)*Hm1*zeta2))/(27.*x)
                         + ((1 - x)*(8*(-5 + 2*x)*(-16 + 67*x)*H12 - 24*(249 - 2*x + 249*x2)*H100
                                     + 16*(230 + 11*x + 257*x2)*H110 - 4*(-32 + 103*x + 76*x2)*H111
                                     + 4*(748 + 475*x + 748*x2)*H1*zeta2))/(27.*x)
                         + (12*(-24 + 26*x + 47*x2)*H30 + 4*(159 + 294*x + 265*x2 + 112*x3)*H210
                            - 4*(-37 - 56*x - 2*x2 + 8*x3)*H00*zeta2)/(9.*(1 + x)) - (208*(2 + x)*H00*zeta3)/3.
                         - (8*(11 + 27*x)*zeta2*zeta3)/3.
                         + (24*(52 + 144*x + 123*x2 + 124*x3)*H20 - 4*(-116 - 371*x - 488*x2 + 40*x3)*H21
                            + 8*(156 + 776*x + 731*x2 + 254*x3)*H0*zeta2
                            - 4*(-1620 + 2203*x - 4100*x2 + 700*x3)*zeta3)/(27.*x)
                         - (8*(1 - x)*(H5 - 2*Hm32 + Hm23 - 8*Hm3m10 - 10*Hm2m20 + 2*Hm2m12 - 2*Hm220 + 4*Hm2m1m10
                                       - 23*Hm2m100 + 16*Hm2000 - 2*Hm3*zeta2 + 7*H3*zeta2 - 8*Hm20*zeta2 - 2*Hm2*zeta3))/3.
                         + ((1 + x)*(12*(8 + x + 8*x2)*Hm13 - 8*(16 - 37*x + 16*x2)*Hm1m20 - 8*(16 + 17*x + 16*x2)*Hm1m12
                                     + 16*(4 - 7*x + 4*x2)*Hm1m1m10 - 36*(8 - 19*x + 8*x2)*Hm1m100
                                     + 96*(-2 + x)*(-1 + 2*x)*Hm1000 - 16*(14 - 11*x + 14*x2)*Hm10*zeta2
                                     - 16*(2 + x + 2*x2)*(2*Hm121 - 5*Hm1m1*zeta2) - 24*(4 - x + 4*x2)*Hm1*zeta3))/(9.*x)
                         + ((1 - x)*(-4*(64 - 5*x + 64*x2)*H13 + 16*(8 + 5*x + 8*x2)*H1m20 - 4*(32 + 29*x + 32*x2)*H112
                                     - 24*(16 + 13*x + 16*x2)*H120 + 4*(104 + 119*x + 104*x2)*H1000
                                     - 28*(16 + 19*x + 16*x2)*H1100 + 16*(4 + 7*x + 4*x2)*H1111
                                     - 32*(1 - 5*x + x2)*H11*zeta2 - 4*(32 - 25*x + 32*x2)*(H1110 - H1*zeta3)))/(9.*x)
                         + (2*(-201 + 5*x)*H0*zeta4)/3. + (-8*(72 + 141*x + 141*x2 + 109*x3 + 28*x4)*H200
                                                           + 12*(16 + 79*x + 120*x2 + 87*x3 + 24*x4)*H2*zeta2
                                                           + 8*(24 - 14*x - 136*x2 + 47*x3 + 136*x4)*H0*zeta3
                                                           + (1840 + 2551*x + 4956*x2 + 4935*x3 + 672*x4)*zeta4)/(9.*x*(1 + x))
                         - (8*(1 + x)*(3*H23 + H41 - 3*Hm130 - 4*H2m20 + 5*H212 + 14*H220 - 2*H311 - 6*Hm1200 + 6*Hm1210
                                       - 19*H2000 + 21*H2100 - H2110 - 4*H2111 + 6*Hm12*zeta2 + 5*H20*zeta2 - 4*H21*zeta2
                                       - 3*Hm100*zeta2 - H000*zeta2 + H2*zeta3 + 6*Hm10*zeta3 + (3*Hm1*zeta4)/2.))/3.
                         + 32*(-1 + 9*x)*zeta5);
  }

  //_________________________________________________________________________________
  C3pvpdf::C3pvpdf():
    Expression()
  {
  }
  double C3pvpdf::Regular(double const& x) const
  {
    // Polylogs
    double *xx  = new double{x};
    int    *nw  = new int{5};
    int    *wn1 = new int{-1};
    int    *wn2 = new int{1};
    double *H   = new double[363];
    apf_hplog_(xx, nw, &H[0], &H[3], &H[12], &H[39], &H[120], wn1, wn2);

    // Defintions
    const double x2    = x * x;
    const double x3    = x * x2;
    const double nc    = 3;
    const double nc2   = nc * nc;
    const double dabc2 = ( nc2 - 1 ) * ( nc2 - 4 ) / nc;

    const double Hm1 = H[HPLogMap({-1})];
    const double H0  = H[HPLogMap({0})];
    const double H1  = H[HPLogMap({1})];

    const double Hm10  = H[HPLogMap({-1,0})];
    const double H00   = H[HPLogMap({0,0})];
    const double H2    = H[HPLogMap({2})];
    const double H10   = H[HPLogMap({1,0})];
    const double H11   = H[HPLogMap({1,1})];
    const double Hm2   = H[HPLogMap({-2})];
    const double Hm1m1 = H[HPLogMap({-1,-1})];

    const double H000    = H[HPLogMap({0,0,0})];
    const double H12     = H[HPLogMap({1,2})];
    const double H20     = H[HPLogMap({2,0})];
    const double H21     = H[HPLogMap({2,1})];
    const double H100    = H[HPLogMap({1,0,0})];
    const double H110    = H[HPLogMap({1,1,0})];
    const double Hm1m10  = H[HPLogMap({-1,-1,0})];
    const double Hm100   = H[HPLogMap({-1,0,0})];
    const double H3      = H[HPLogMap({3})];
    const double Hm12    = H[HPLogMap({-1,2})];
    const double Hm3     = H[HPLogMap({-3})];
    const double Hm2m1   = H[HPLogMap({-2,-1})];
    const double Hm20    = H[HPLogMap({-2,0})];

    const double H4       = H[HPLogMap({4})];
    const double H22      = H[HPLogMap({2,2})];
    const double H30      = H[HPLogMap({3,0})];
    const double H210     = H[HPLogMap({2,1,0})];
    const double H13      = H[HPLogMap({1,3})];
    const double H31      = H[HPLogMap({3,1})];
    const double H112     = H[HPLogMap({1,1,2})];
    const double H120     = H[HPLogMap({1,2,0})];
    const double H200     = H[HPLogMap({2,0,0})];
    const double H0000    = H[HPLogMap({0,0,0,0})];
    const double H1000    = H[HPLogMap({1,0,0,0})];
    const double H1100    = H[HPLogMap({1,1,0,0})];
    const double H1110    = H[HPLogMap({1,1,1,0})];
    const double Hm22     = H[HPLogMap({-2,2})];
    const double Hm2m10   = H[HPLogMap({-2,-1,0})];
    const double Hm200    = H[HPLogMap({-2,0,0})];
    const double Hm13     = H[HPLogMap({-1,3})];
    const double Hm1m20   = H[HPLogMap({-1,-2,0})];
    const double Hm1m12   = H[HPLogMap({-1,-1,2})];
    const double Hm120    = H[HPLogMap({-1,2,0})];
    const double Hm121    = H[HPLogMap({-1,2,1})];
    const double Hm1m1m10 = H[HPLogMap({-1,-1,-1,0})];
    const double Hm1m100  = H[HPLogMap({-1,-1,0,0})];
    const double Hm1000   = H[HPLogMap({-1,0,0,0})];
    const double Hm30     = H[HPLogMap({-3,0})];
    const double H1m20    = H[HPLogMap({1,-2,0})];
    const double Hm130    = H[HPLogMap({-1,3,0})];
    const double Hm1200   = H[HPLogMap({-1,2,0,0})];
    const double Hm1210   = H[HPLogMap({-1,2,1,0})];

    const double H5        = H[HPLogMap({5})];
    const double Hm40      = H[HPLogMap({-4,0})];
    const double H23       = H[HPLogMap({2,3})];
    const double H41       = H[HPLogMap({4,1})];
    const double Hm300     = H[HPLogMap({-3,0,0})];
    const double H2m20     = H[HPLogMap({2,-2,0})];
    const double H220      = H[HPLogMap({2,2,0})];
    const double H300      = H[HPLogMap({3,0,0})];
    const double H2000     = H[HPLogMap({2,0,0,0})];
    const double H2100     = H[HPLogMap({2,1,0,0})];
    const double H2110     = H[HPLogMap({2,1,1,0})];
    const double H00000    = H[HPLogMap({0,0,0,0,0})];
    const double H212      = H[HPLogMap({2,1,2})];
    const double Hm32      = H[HPLogMap({-3,2})];
    const double Hm23      = H[HPLogMap({-2,3})];
    const double H32       = H[HPLogMap({3,2})];
    const double Hm3m10    = H[HPLogMap({-3,-1,0})];
    const double Hm2m20    = H[HPLogMap({-2,-2,0})];
    const double Hm2m12    = H[HPLogMap({-2,-1,2})];
    const double Hm220     = H[HPLogMap({-2,2,0})];
    const double H310      = H[HPLogMap({3,1,0})];
    const double Hm2m1m10  = H[HPLogMap({-2,-1,-1,0})];
    const double Hm2m100   = H[HPLogMap({-2,-1,0,0})];
    const double Hm2000    = H[HPLogMap({-2,0,0,0})];

    // Delete pointers
    delete xx;
    delete nw;
    delete wn1;
    delete wn2;
    delete[] H;

    return
      (dabc2*TR*((171520*(1 - x))/27. + (2816*(-2 + 9*x)*H0)/27. + (256*(1 - x)*(-36 + 629*x)*H1)/(27.*x)
                 + (64*(355 + 1669*x)*H2)/27. - (64*(435 + 627*x + 512*x2)*H3)/27. + (256*(15 + 15*x + 28*x2)*H4)/9.
                 - (256*(7 + 9*x)*H5)/3. - (512*(-1 + 5*x)*Hm40)/3. - (512*(21 - 9*x + 16*x2)*Hm30)/9.
                 + (128*(207 + 369*x + 256*x2)*Hm20)/27. - (256*(27 + x + 16*x2)*Hm22)/3. - (34048*(1 + x)*Hm10)/27.
                 + (64*(-794 + 1171*x)*H00)/27. - (128*(1 - x)*(153*H10 - 91*H11))/9. + (512*(-5 + 4*x)*H20)/3.
                 + (5248*(1 + x)*H21)/9. + (256*(6 + 3*x + 8*x2)*H22)/3. - (128*(-3 + 27*x + 16*x2)*H31)/9.
                 - (256*(-5 + 13*x)*Hm300)/3. + (512*(1 - x + 4*x2)*Hm2m10)/3. - (128*(93 - 105*x + 88*x2)*Hm200)/9.
                 - 256*(1 + x)*(Hm1m20 + Hm120) - (128*(345 + 261*x + 256*x2)*H000)/27. - 512*(1 - x)*H1m20
                 - (64*(1 - x)*(108*H12 + 193*H100 + 288*H110))/9. - (512*(-1 + 2*x)*H300)/3. + (1024*(1 + 2*x2)*H0000)/3.
                 - (64*(1715 + 841*x)*zeta2)/27. + (1024*(7 + 5*x2)*Hm2*zeta2)/3. + (64*(489 + 825*x + 1024*x2)*H0*zeta2)/27.
                 - (1024*(1 - x)*(16 + 13*x + 16*x2)*H1*zeta2)/(27.*x) + (256*(-1 + 9*x)*H3*zeta2)/3.
                 - (1024*pow(1 + x,3)*Hm1m1*zeta2)/(3.*x)
                 + ((1 + x)*(256*(128 - 95*x + 128*x2)*Hm12 - 256*(128 + 67*x + 128*x2)*Hm1m10
                             + 256*(128 + 25*x + 128*x2)*Hm100 - 384*(128 - 41*x + 128*x2)*Hm1*zeta2))/(27.*x)
                 + (256*(540 - 141*x + 320*x2)*zeta3)/27. - (256*(-7 + 5*x)*zeta2*zeta3)/3.
                 - (256*(1 - x)*(2*Hm32 - Hm23 - 2*H32 + 8*Hm3m10 + 2*Hm2m20 - 2*Hm2m12 + 2*Hm220 + 5*H310
                                 + 4*Hm2m1m10 + 11*Hm2m100 - 8*Hm2000 + 2*Hm3*zeta2 + 4*Hm2m1*zeta2 - 2*Hm2*zeta3))/3.
                 + ((1 + x)*(-384*(8 + x + 8*x2)*Hm13 + 256*(16 + 17*x + 16*x2)*Hm1m12 + 512*(4 - 7*x + 4*x2)*Hm1m1m10
                             + 384*(8 - 29*x + 8*x2)*Hm1m100 - 2048*(1 - 4*x + x2)*Hm1000
                             + 512*(2 + x + 2*x2)*(2*Hm121 + 3*Hm10*zeta2) + 256*(4 + 11*x + 4*x2)*Hm1*zeta3))/(9.*x)
                 + ((1 - x)*(128*(32 + 65*x + 32*x2)*H13 - 384*(16 + 19*x + 16*x2)*H112 - 768*(8 - x + 8*x2)*H120
                             + 128*(32 - 7*x + 32*x2)*H1000 - 128*(16 + 73*x + 16*x2)*H1100
                             - 384*(32 + 29*x + 32*x2)*H1110 - 128*(80 + 77*x + 80*x2)*H10*zeta2
                             - 512*(14 + 11*x + 14*x2)*H11*zeta2 - 640*(16 + x + 16*x2)*H1*zeta3))/(9.*x)
                 - (512*(3*H00000 - 13*H00*zeta3))/3. - (64*(39 + 5*x)*H0*zeta4)/3.
                 + (384*(12 + 8*x + 9*x2 + 16*x3)*H30 + 128*(-15 - 6*x + 7*x2 + 16*x3)*H200
                    + 384*(15 + 12*x + 35*x2 + 32*x3)*H210 + 384*(7 + 2*x + 9*x2 + 8*x3)*H2*zeta2
                    - 384*(26 + 42*x + 43*x2 + 24*x3)*H00*zeta2 - 256*(-9 - 18*x + 26*x2 + 44*x3)*H0*zeta3
                    + 64*(15 - 243*x - 173*x2 + 76*x3)*zeta4)/(9.*(1 + x))
                 + (256*(1 + x)*(9*H23 + 3*H41 + 3*Hm130 - 4*H2m20 - 9*H212 - 2*H220 + 6*Hm1200 - 6*Hm1210
                                 + H2000 - 9*H2100 - 15*H2110 - 6*Hm12*zeta2 - 13*H20*zeta2 - 8*H21*zeta2
                                 + 3*Hm100*zeta2 + 9*H000*zeta2 - 5*H2*zeta3 - 6*Hm10*zeta3 - (3*Hm1*zeta4)/2.))/3.
                 - (256*(53 + 15*x)*zeta5)/3.))/(64.*nc);
  }

  //_________________________________________________________________________________
  C3qgpdf::C3qgpdf(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C3qgpdf::Regular(double const& x) const
  {
    // Polylogs
    double *xx  = new double{x};
    int    *nw  = new int{5};
    int    *wn1 = new int{-1};
    int    *wn2 = new int{1};
    double *H   = new double[363];
    apf_hplog_(xx, nw, &H[0], &H[3], &H[12], &H[39], &H[120], wn1, wn2);

    // Defintions
    const double x2   = x * x;
    const double x3   = x * x2;
    const double x4   = x * x3;
    const double x5   = x * x4;
    const double opx2 = ( 1 + x ) * ( 1 + x );
    const double CF2  = CF * CF;
    const double CA2  = CA * CA;
    const double TR2  = TR * TR;

    const double Hm1 = H[HPLogMap({-1})];
    const double H0  = H[HPLogMap({0})];
    const double H1  = H[HPLogMap({1})];

    const double Hm10  = H[HPLogMap({-1,0})];
    const double H00   = H[HPLogMap({0,0})];
    const double H2    = H[HPLogMap({2})];
    const double H10   = H[HPLogMap({1,0})];
    const double H11   = H[HPLogMap({1,1})];
    const double Hm2   = H[HPLogMap({-2})];
    const double Hm1m1 = H[HPLogMap({-1,-1})];

    const double Hm20    = H[HPLogMap({-2,0})];
    const double H000    = H[HPLogMap({0,0,0})];
    const double H12     = H[HPLogMap({1,2})];
    const double H20     = H[HPLogMap({2,0})];
    const double H21     = H[HPLogMap({2,1})];
    const double H100    = H[HPLogMap({1,0,0})];
    const double H110    = H[HPLogMap({1,1,0})];
    const double Hm1m10  = H[HPLogMap({-1,-1,0})];
    const double Hm100   = H[HPLogMap({-1,0,0})];
    const double H111    = H[HPLogMap({1,1,1})];
    const double H3      = H[HPLogMap({3})];
    const double Hm12    = H[HPLogMap({-1,2})];
    const double Hm3     = H[HPLogMap({-3})];
    const double Hm2m1   = H[HPLogMap({-2,-1})];
    const double H1m2    = H[HPLogMap({1,-2})];
    const double Hm1m2	 = H[HPLogMap({-1,-2})];
    const double Hm1m1m1 = H[HPLogMap({-1,-1,-1})];

    const double H4       = H[HPLogMap({4})];
    const double H22      = H[HPLogMap({2,2})];
    const double H30      = H[HPLogMap({3,0})];
    const double H210     = H[HPLogMap({2,1,0})];
    const double H13      = H[HPLogMap({1,3})];
    const double H31      = H[HPLogMap({3,1})];
    const double H112     = H[HPLogMap({1,1,2})];
    const double H120     = H[HPLogMap({1,2,0})];
    const double H121     = H[HPLogMap({1,2,1})];
    const double H200     = H[HPLogMap({2,0,0})];
    const double H0000    = H[HPLogMap({0,0,0,0})];
    const double H1000    = H[HPLogMap({1,0,0,0})];
    const double H1100    = H[HPLogMap({1,1,0,0})];
    const double H1110    = H[HPLogMap({1,1,1,0})];
    const double H211     = H[HPLogMap({2,1,1})];
    const double Hm22     = H[HPLogMap({-2,2})];
    const double Hm2m10   = H[HPLogMap({-2,-1,0})];
    const double Hm200    = H[HPLogMap({-2,0,0})];
    const double Hm13     = H[HPLogMap({-1,3})];
    const double Hm1m20   = H[HPLogMap({-1,-2,0})];
    const double Hm1m12   = H[HPLogMap({-1,-1,2})];
    const double Hm120    = H[HPLogMap({-1,2,0})];
    const double Hm121    = H[HPLogMap({-1,2,1})];
    const double Hm1m1m10 = H[HPLogMap({-1,-1,-1,0})];
    const double Hm1m100  = H[HPLogMap({-1,-1,0,0})];
    const double Hm1000   = H[HPLogMap({-1,0,0,0})];
    const double Hm30     = H[HPLogMap({-3,0})];
    const double H1m20    = H[HPLogMap({1,-2,0})];
    const double Hm130    = H[HPLogMap({-1,3,0})];
    const double Hm1200   = H[HPLogMap({-1,2,0,0})];
    const double Hm1210   = H[HPLogMap({-1,2,1,0})];
    const double H1111 	  = H[HPLogMap({1,1,1,1})];

    const double H5         = H[HPLogMap({5})];
    const double Hm40       = H[HPLogMap({-4,0})];
    const double H23        = H[HPLogMap({2,3})];
    const double H41        = H[HPLogMap({4,1})];
    const double Hm300      = H[HPLogMap({-3,0,0})];
    const double H2m20      = H[HPLogMap({2,-2,0})];
    const double H220       = H[HPLogMap({2,2,0})];
    const double H300       = H[HPLogMap({3,0,0})];
    const double H2000      = H[HPLogMap({2,0,0,0})];
    const double H2100      = H[HPLogMap({2,1,0,0})];
    const double H2110      = H[HPLogMap({2,1,1,0})];
    const double H00000     = H[HPLogMap({0,0,0,0,0})];
    const double H212       = H[HPLogMap({2,1,2})];
    const double Hm32       = H[HPLogMap({-3,2})];
    const double Hm23       = H[HPLogMap({-2,3})];
    const double H14        = H[HPLogMap({1,4})];
    const double H32        = H[HPLogMap({3,2})];
    const double H40        = H[HPLogMap({4,0})];
    const double Hm3m10     = H[HPLogMap({-3,-1,0})];
    const double Hm2m20     = H[HPLogMap({-2,-2,0})];
    const double Hm2m12     = H[HPLogMap({-2,-1,2})];
    const double Hm220      = H[HPLogMap({-2,2,0})];
    const double Hm221      = H[HPLogMap({-2,2,1})];
    const double H1m30      = H[HPLogMap({1,-3,0})];
    const double H1m22      = H[HPLogMap({1,-2,2})];
    const double H113       = H[HPLogMap({1,1,3})];
    const double H122       = H[HPLogMap({1,2,2})];
    const double H130       = H[HPLogMap({1,3,0})];
    const double H131       = H[HPLogMap({1,3,1})];
    const double H221       = H[HPLogMap({2,2,1})];
    const double H310       = H[HPLogMap({3,1,0})];
    const double Hm2m1m10   = H[HPLogMap({-2,-1,-1,0})];
    const double Hm2m100    = H[HPLogMap({-2,-1,0,0})];
    const double Hm2000     = H[HPLogMap({-2,0,0,0})];
    const double H1m2m10    = H[HPLogMap({1,-2,-1,0})];
    const double H1m200     = H[HPLogMap({1,-2,0,0})];
    const double H11m20     = H[HPLogMap({1,1,-2,0})];
    const double H1112      = H[HPLogMap({1,1,1,2})];
    const double H1120      = H[HPLogMap({1,1,2,0})];
    const double H1121      = H[HPLogMap({1,1,2,1})];
    const double H1200      = H[HPLogMap({1,2,0,0})];
    const double H1210      = H[HPLogMap({1,2,1,0})];
    const double H10000     = H[HPLogMap({1,0,0,0,0})];
    const double H11000     = H[HPLogMap({1,1,0,0,0})];
    const double H11100     = H[HPLogMap({1,1,1,0,0})];
    const double H311       = H[HPLogMap({3,1,1})];
    const double H2111	    = H[HPLogMap({2,1,1,1})];
    const double H1211	    = H[HPLogMap({1,2,1,1})];
    const double H11110	    = H[HPLogMap({1,1,1,1,0})];
    const double H11111	    = H[HPLogMap({1,1,1,1,1})];
    const double Hm1m121    = H[HPLogMap({-1,-1,2,1})];
    const double Hm10000    = H[HPLogMap({-1,0,0,0,0})];
    const double Hm1m1m20   = H[HPLogMap({-1,-1,-2,0})];
    const double Hm14	    = H[HPLogMap({-1,4})];
    const double Hm1m30	    = H[HPLogMap({-1,-3,0})];
    const double Hm1m22	    = H[HPLogMap({-1,-2,2})];
    const double Hm1m13	    = H[HPLogMap({-1,-1,3})];
    const double Hm122	    = H[HPLogMap({-1,2,2})];
    const double Hm131	    = H[HPLogMap({-1,3,1})];
    const double Hm1m2m10   = H[HPLogMap({-1,-2,-1,0})];
    const double Hm1m200    = H[HPLogMap({-1,-2,0,0})];
    const double Hm1m1m12   = H[HPLogMap({-1,-1,-1,2})];
    const double Hm1m120    = H[HPLogMap({-1,-1,2,0})];
    const double Hm1211	    = H[HPLogMap({-1,2,1,1})];
    const double Hm1m1m1m10 = H[HPLogMap({-1,-1,-1,-1,0})];
    const double Hm1m1m100  = H[HPLogMap({-1,-1,-1,0,0})];
    const double Hm1m1000   = H[HPLogMap({-1,-1,0,0,0})];

    // Delete pointers
    delete xx;
    delete nw;
    delete wn1;
    delete wn2;
    delete[] H;

    return
      2*(CA*_nf*TR2*((4*(-1582 - 8697*x - 24198*x2 + 43324*x3))/(729.*x) + (8*(-523 - 8119*x + 2280*x2)*H0)/243.
                     + (8*(-36 + 103*x - 1064*x2 + 668*x3)*H1)/(243.*x) + (16*(-3 - 23*x + 3*x2)*H2)/27.
                     - (16*(10 + 21*x + 20*x2)*Hm20)/9. - (8*(85 + 290*x + 299*x2)*Hm10)/81.
                     - (8*(32 + 1407*x + 391*x2)*H00)/81. - (8*(85 - 326*x + 290*x2)*H11)/81.
                     + (8*(12 + 65*x + 8*x2)*H20)/9. - (8*x*(5*H3 - 2*H21))/9. - (16*(25 + 62*x + 62*x2)*Hm100)/27.
                     + (32*(-4 - 74*x + 23*x2)*H000)/27. - (32*(10 - 23*x + 23*x2)*(H12 - H111))/27.
                     + (16*(1 + 6*x + 2*x2)*H200)/3. - (32*(-1 + 12*x)*H0000)/9. + (8*(12 + 27*x + 8*x2)*H0*zeta2)/9.
                     + (-24*(-40 + 3*x - 123*x2 + 154*x3)*H10 - 8*(-120 - 18*x - 172*x2 + 435*x3)*zeta2)/(81.*x)
                     + (80*(1 + 2*x + 2*x2)*(2*Hm1m10 + Hm1*zeta2))/9.
                     + (-16*(-24 + 31*x - 128*x2 + 116*x3)*H100 + 32*(-6 + 4*x - 29*x2 + 26*x3)*H110
                        + 16*(-12 + 23*x - 88*x2 + 82*x3)*H1*zeta2 + 48*(-4 + 18*x + 7*x2 + 30*x3)*zeta3)/(27.*x)
                     - (32*(1 + 2*x + 2*x2)*(3*Hm30 - 3*Hm2m10 + (9*Hm200)/2. - 3*Hm1m20 + 3*Hm1m1m10
                                             - (9*Hm1m100)/2. + 4*Hm1000 - (3*Hm2*zeta2)/2.
                                             + (3*Hm1m1*zeta2)/2. - (3*Hm10*zeta2)/2. - (3*Hm1*zeta3)/2.))/9.
                     - (16*(1 - 2*x + 2*x2)*(2*H13 - 8*H112 - 2*H120 - 5*H121 - 3*H1000 + 4*H1100 - 4*H1110
                                             + 5*H1111 + H10*zeta2 + 2*H11*zeta2 - 10*H1*zeta3))/9.
                     + (16*x*(12*H30 - 12*H210 - 12*H2*zeta2 - 19*zeta4))/9.)
         + CF*_nf*TR2*(-(-160736 - 1294539*x + 1012650*x2 + 417704*x3)/(1458.*x)
                       + (8*(20752 + 7693*x + 17172*x2)*H0)/243. + (8*(-355 + 155*x + 268*x2)*H1)/243.
                       + (224*(4 - 11*x + 11*x2)*H2)/81. - (4*(-6023 + 6559*x + 7424*x2)*H00)/81.
                       + (8*(112 - 353*x + 308*x2)*H11)/81. + (4*(1805 - 1282*x + 1344*x2)*H000)/27.
                       - (32*(10 - 23*x + 23*x2)*(H21 + H111))/27. - (16*(-47 + 94*x + 44*x2)*H0000)/9.
                       + (80*(1 - 2*x + 2*x2)*(H211 + 2*H1000 + H1111))/9. - 64*(-1 + 2*x)*H00000
                       + (32*(-96 + 131*x - 403*x2 + 382*x3)*H10 + 8*(-384 + 356*x - 961*x2 + 877*x3)*zeta2)/(81.*x)
                       - (32*x*(-9 + 10*x)*(H20 + H0*zeta2))/9. - (32*(1 - x)*(4 + 13*x + 22*x2)*(H110 + H1*zeta2))/(9.*x)
                       + (-16*(-48 - 133*x - 28*x2 + 184*x3)*H100 - 16*(24 + 154*x - 287*x2 + 221*x3)*zeta3)/(27.*x)
                       + (32*(3 + 6*x + 4*x2)*(H30 + 2*H200 - H210 - H2*zeta2 + H00*zeta2 + H0*zeta3))/3.
                       - (8*(73 + 106*x + 104*x2)*zeta4)/9.)
         + CF2*TR*((2127 + 2406*x - 5020*x2)/24. + ((235 - 734*x - 708*x2)*H0)/3. - (2*(-505 + 191*x + 354*x2)*H1)/3.
                   - (4*(-70 - 343*x + 151*x2)*H2)/3. + (4*(49 - 11*x + 24*x2)*H3)/3. - (4*(5 - 136*x + 50*x2)*H4)/3.
                   + (64*(7 + 16*x2)*Hm40)/3. + (16*(47 + 46*x + 53*x2)*Hm30)/3. + (8*(135 + 98*x + 63*x2)*Hm20)/3.
                   + (4*(395 + 634*x + 251*x2)*Hm10)/3. + ((15 - 857*x - 1012*x2)*H00)/3. - (4*(77 - 26*x + 2*x2)*H10)/3.
                   - (2*(110 - 353*x + 302*x2)*H11)/3. + (16*(2 - 9*x + 6*x2)*H12)/3. - (8*(-7 - 12*x + 25*x2)*H13)/3.
                   - (4*(-7 + 108*x + 4*x2)*H20)/3. + (4*(35 - 47*x + 27*x2)*H21)/3. + (8*(2 + 13*x2)*H22)/3.
                   - (32*(3 - 8*x + 3*x2)*H23)/3. + (8*(7 - 98*x + 33*x2)*H30)/3. - (8*(-3 - 8*x + 6*x2)*H31)/3.
                   + (64*(1 - 2*x + 3*x2)*H32)/3. + (32*(11 + 18*x + 26*x2)*Hm300)/3. - (64*(3 + 4*x + 8*x2)*Hm2m20)/3.
                   + (8*(55 + 138*x + 112*x2)*Hm200)/3. - (112*(1 + x)*(3 + 7*x)*Hm1m20)/3.
                   + (4*(183 + 384*x + 209*x2)*Hm100)/3. + 64*x*(1 + x)*Hm120 + ((65 - 694*x - 792*x2)*H000)/3.
                   + (32*(17 - 22*x + 2*x2)*H1m20)/3. + (4*(1 + x)*(-21 + 11*x)*H100)/3. - (8*(-6 - 5*x + 2*x2)*H110)/3.
                   + (4*(15 - 38*x + 27*x2)*H111)/3. + (8*(8 - 12*x + 13*x2)*H112)/3. + (8*(23 - 68*x + 57*x2)*H120)/3.
                   - (8*(-5 - 4*x + 6*x2)*H121)/3. + (64*(3 + 2*x2)*H2m20)/3. + (16*(2 - 35*x + 4*x2)*H200)/3.
                   - (8*(-13 + 19*x2)*H210)/3. - (4*(-1 + 2*x)*(3 + 8*x)*H211)/3. + (16*(5 - 6*x + 2*x2)*H212)/3.
                   + (32*(5 - 12*x + 17*x2)*H220)/3. + (16*(7 - 14*x + 16*x2)*H221)/3. - (32*(1 + 14*x + 16*x2)*H300)/3.
                   + (16*(13 - 2*x + 48*x2)*H310)/3. + (16*(13 - 26*x + 28*x2)*H311)/3. - (32*(9 + 10*x + 24*x2)*Hm2m100)/3.
                   + (64*(1 + x + 3*x2)*Hm2000)/3. - (64*(1 + x)*(9 + 14*x)*Hm1m100)/3. + (8*(27 + 82*x + 61*x2)*Hm1000)/3.
                   - (8*(-4 + 21*x + 94*x2)*H0000)/3. - (8*(17 - 56*x + 33*x2)*H1000)/3. + (8*(45 - 50*x + 8*x2)*H1100)/3.
                   - (8*(26 - 42*x + 19*x2)*H1110)/3. - (4*(1 - 20*x + 16*x2)*H1111)/3. - (16*(11 - 26*x + 40*x2)*H2000)/3.
                   + (32*(1 + 2*x + 4*x2)*H2100)/3. - (16*(-5 + 14*x + 4*x2)*H2110)/3. + (8*(31 - 62*x + 60*x2)*H2111)/3.
                   - (16*(3 + 4*x + 24*x2)*H00000)/3. + ((-241 - 81*x + 1768*x2)*zeta2)/3.
                   + ((-220 - 177*x + 252*x2)*H0*zeta2)/3. + (2*(40 + 27*x)*H1*zeta2)/3. - (8*(-4 + 7*x + 7*x2)*H2*zeta2)/3.
                   + (64*(1 + 2*x + 7*x2)*H3*zeta2)/3. - (64*x*(2 + 3*x)*Hm20*zeta2)/3. - (8*(-13 + 2*x + 21*x2)*Hm10*zeta2)/3.
                   + ((17 - 708*x + 656*x2)*H00*zeta2)/3. + (8*(30 - 73*x + 61*x2)*H10*zeta2)/3.
                   - (8*(17 - 15*x + 7*x2)*H11*zeta2)/3. + (32*(6 - 10*x + 9*x2)*H20*zeta2)/3.
                   - (4*(-1 + 18*x + 30*x2)*H21*zeta2)/3. - (32*(5 + 6*x + 10*x2)*(2*Hm3m10 + Hm3*zeta2))/3.
                   - (8*(21 + 38*x + 36*x2)*(2*Hm2m10 + Hm2*zeta2))/3. - (4*(41 + 96*x + 63*x2)*(2*Hm1m10 + Hm1*zeta2))/3.
                   + 32*(1 + x)*(1 + 3*x)*(2*Hm1m1m10 + Hm1m1*zeta2) + (4*(317 - 758*x + 757*x2)*zeta3)/3.
                   - 32*(1 + x)*(1 + 2*x)*Hm1*zeta3 + (4*(57 - 344*x + 398*x2)*H0*zeta3)/3.
                   - (8*(1 - x)*(-73 + 175*x)*H1*zeta3)/3. + (32*(20 - 34*x + 39*x2)*H2*zeta3)/3.
                   - (4*(197 - 122*x + 410*x2)*zeta2*zeta3)/3. + (32*(3 + 6*x + 8*x2)*(2*Hm2m1m10 + Hm2m1*zeta2 - Hm2*zeta3))/3.
                   - (2*(1 - 2*x + 4*x2)*(52*H5 + 16*H40 - 16*H41 - 55*H000*zeta2 - 172*H00*zeta3))/3.
                   + ((688 + 907*x + 379*x2)*zeta4)/3. + (4*(109 - 118*x + 348*x2)*H0*zeta4)/3.
                   - (16*(1 + 2*x + 2*x2)*(22*Hm1m30 - 6*Hm130 - 20*Hm1m2m10 + 30*Hm1m200 - 16*Hm1m1m20 - 22*Hm1200 + 12*Hm1210
                                           + 16*Hm1m1m1m10 - 24*Hm1m1m100 + 6*Hm1m1000 - 5*Hm10000 - 10*Hm1m2*zeta2 + 12*Hm12*zeta2
                                           + 8*Hm1m1m1*zeta2 - 6*Hm1m10*zeta2 - 8*Hm1m1*zeta3 + 4*Hm10*zeta3 + 22*Hm1*zeta4))/3.
                   - 4*(1 - 2*x + 2*x2)*(16*H14 - (40*H1m30)/3. + (8*H113)/3. - 8*H122 - (8*H130)/3. - 4*H131 + (16*H1m200)/3.
                                         - (16*H11m20)/3. - (4*H1112)/3. - (68*H1120)/3. - (28*H1121)/3. - (20*H1200)/3.
                                         - (52*H1210)/3. - (56*H1211)/3. + (28*H10000)/3. + (80*H11000)/3. - 4*H11100
                                         + (4*H11110)/3. - 20*H11111 - (8*H12*zeta2)/3. - 17*H100*zeta2 - (32*H110*zeta2)/3.
                                         + 5*H111*zeta2 - 52*H10*zeta3 - 52*H11*zeta3 - (100*H1*zeta4)/3.)
                   - (8*(63 - 454*x + 84*x2)*zeta5)/3.)
         + CA2*TR*((2*(-470494 + 450714*x - 2175249*x2 + 2153647*x3))/(729.*x)
                   + (-2*(17152 + 167620*x + 62809*x2 + 747175*x3)*H0 - 2*(-6528 - 50692*x + 24401*x2 + 30151*x3)*H1)/(243.*x)
                   + (4*(-321 - 483*x + 1460*x2)*H3)/27. - (8*(51 - 48*x + 137*x2)*H4)/9. - (4*(9 - 58*x + 72*x2)*H5)/3.
                   + (8*(-41 + 74*x + 8*x2)*Hm40)/3. + (8*(216 + 180*x + 817*x2)*Hm30)/9. + (8*(29 + 46*x + 56*x2)*Hm32)/3.
                   - (4*(-213 - 30*x + 2426*x2)*Hm20)/27. + (8*(31 + 37*x + 96*x2)*Hm22)/3. + (4*(49 + 110*x + 100*x2)*Hm23)/3.
                   + (2*(20278 + 29235*x + 126201*x2)*H00)/81. - (4*(29 + 66*x + 28*x2)*H23)/3. + (4*(51 - 135*x + 317*x2)*H31)/9.
                   - (8*(19 + 50*x + 20*x2)*H32)/3. - 192*x*H40 + 4*(-3 - 14*x + 8*x2)*H41 - (16*(3 + 22*x + 10*x2)*Hm3m10)/3.
                   + (4*(-35 + 214*x + 52*x2)*Hm300)/3. + (8*(13 - 26*x + 8*x2)*Hm2m20)/3. - (8*(14 + 37*x + 78*x2)*Hm2m10)/3.
                   - (8*(39 + 66*x + 76*x2)*Hm2m12)/3. + (4*(297 + 246*x + 1522*x2)*Hm200)/9. - (8*(3 + 18*x + 8*x2)*Hm220)/3.
                   + (56*(1 + x)*Hm120)/3.
                   + ((1 + x)*(12*(8 + 9*x + 60*x2)*Hm13 - 8*(16 + 53*x + 112*x2)*Hm1m12 - 32*(2 + 4*x + 17*x2)*Hm121))/(9.*x)
                   + 4*(3 + 6*x + 4*x2)*Hm130 - (4*(1679 - 1712*x + 7385*x2)*H000)/27. + (16*(7 + 10*x + 8*x2)*H2m20)/3.
                   - (4*(9 - 62*x + 16*x2)*H212)/3. + (8*(-1 - 62*x + 24*x2)*H220)/3. - (8*(1 - 26*x + 2*x2)*H221)/3.
                   + (8*(17 - 100*x + 6*x2)*H300)/3. + (4*(-25 - 86*x + 4*x2)*H310)/3. + (16*(3 + 8*x)*H311)/3.
                   - 32*(1 + x2)*Hm2m1m10 - (4*(5 + 254*x + 80*x2)*Hm2m100)/3. + (8*(-7 + 78*x + 10*x2)*Hm2000)/3.
                   + (8*(-1 - 2*x + 4*x2)*Hm1210)/3. - (8*(-55 - 345*x + 108*x2)*H0000)/9. - (4*(-11 - 118*x + 52*x2)*H2000)/3.
                   - (4*(43 + 198*x + 20*x2)*H2100)/3. + (4*(27 + 94*x + 20*x2)*H2110)/3. + (16*(1 + 4*x)*H2111)/3.
                   + 16*(-4 + 15*x)*H00000 - (16*(16 + 34*x + 33*x2)*Hm3*zeta2)/3. - (4*(76 + 111*x + 270*x2)*Hm2*zeta2)/3.
                   + (4*(-13 + 18*x + 20*x2)*H3*zeta2)/3. + (8*(33 + 66*x + 70*x2)*Hm2m1*zeta2)/3.
                   - (8*(17 + 90*x + 50*x2)*Hm20*zeta2)/3. + (4*(47 - 2*x + 88*x2)*H20*zeta2)/3.
                   + (32*(6 - x + 8*x2)*H21*zeta2)/3. - 4*(7 + 14*x + 16*x2)*Hm100*zeta2 + 12*(-1 + 2*x)*(-1 + 4*x)*H000*zeta2
                   + (-4*(-1256 - 702*x - 12231*x2 + 98*x3)*H2 + 2*(5436 + 17579*x + 33484*x2 + 23245*x3)*Hm10
                      - 2*(-9490 + 9435*x - 43401*x2 + 43978*x3)*H10 + 2*(-918 + 7796*x - 12541*x2 + 6748*x3)*H11
                      - 2*(-12536 + 11376*x - 58657*x2 + 49619*x3)*zeta2)/(81.*x)
                   + (-8*(116 + 96*x + 393*x2 + 422*x3)*Hm12 - 4*(-208 + 1463*x - 1678*x2 + 92*x3)*H12
                      - 4*(-116 - 111*x - 861*x2 + 445*x3)*H21 + 8*(40 + 54*x + 276*x2 + 7*x3)*Hm1m10
                      - 4*(484 + 124*x + 431*x2 + 393*x3)*Hm100 - 4*(40 + 175*x - 194*x2 + 337*x3)*H111
                      + 4*(272 + 246*x + 1062*x2 + 851*x3)*Hm1*zeta2)/(27.*x) + 8*(5 + 10*x + 8*x2)*(Hm1200 - Hm12*zeta2)
                   + (36*(-8 - 89*x - 138*x2 - 37*x3 + 19*x4)*H30 + 4*(84 + 393*x + 1024*x2 + 1133*x3 + 436*x4)*H210
                      + 4*(30 - 321*x - 323*x2 + 473*x3 + 436*x4)*H00*zeta2)/(9.*opx2) - (16*(8 + 22*x + 19*x2)*Hm2*zeta3)/3.
                   + 4*(27 + 30*x + 52*x2)*H2*zeta3 - (16*(14 + 28*x + 25*x2)*Hm10*zeta3)/3. - (16*(27 + 41*x)*H00*zeta3)/3.
                   - (4*(93 + 242*x + 136*x2)*zeta2*zeta3)/3.
                   + (6*(208 + 556*x + 1118*x2 + 2811*x3 + 2059*x4)*H20 + 4*(-1494 - 235*x - 6050*x2 + 123*x3 + 7486*x4)*H100
                      - 2*(-1984 + 805*x - 6983*x2 - 1361*x3 + 8519*x4)*H110
                      + 2*(624 + 2238*x + 5946*x2 + 7511*x3 + 3233*x4)*H0*zeta2
                      - 2*(-1544 - 1467*x - 6989*x2 + 1403*x3 + 8577*x4)*H1*zeta2
                      - 2*(-3288 + 2178*x - 8142*x2 - 13207*x3 + 509*x4)*zeta3)/(27.*x*(1 + x))
                   + (4*(-112 + 217*x - 761*x2 + 624*x3)*H13 + 12*(-16 - 59*x - 172*x2 + 97*x3)*H22
                      - 8*(16 + 39*x + 195*x2 + 205*x3)*Hm1m20 - 16*(-8 - 60*x + 69*x2 + 8*x3)*H1m20
                      - 4*(-64 + 220*x - 539*x2 + 471*x3)*H112 + 4*(-48 - 25*x - 352*x2 + 457*x3)*H120
                      - 4*(-48 + 124*x - 404*x2 + 383*x3)*H121 - 16*(-2 + 6*x2 + 53*x3)*H211
                      + 8*(8 + 51*x + 84*x2 + 74*x3)*Hm1m1m10 - 12*(24 + 16*x + 245*x2 + 286*x3)*Hm1m100
                      + 4*(48 + x + 494*x2 + 647*x3)*Hm1000 - 4*(-8 - 147*x - 63*x2 + 233*x3)*H1000
                      + 4*(-112 - 4*x - 685*x2 + 845*x3)*H1100 - 4*(-112 + 233*x - 811*x2 + 734*x3)*H1110
                      - 4*(-8 - 61*x + 14*x2)*H1111 + 4*(40 + 189*x + 414*x2 + 298*x3)*Hm1m1*zeta2
                      - 4*(56 - 30*x + 264*x2 + 401*x3)*Hm10*zeta2 - 4*(-96 + 97*x - 389*x2 + 287*x3)*H10*zeta2
                      - 8*(-20 - 5*x - 98*x2 + 112*x3)*H11*zeta2 - 12*(8 + 33*x + 72*x2 + 58*x3)*Hm1*zeta3
                      - 4*(-176 + 38*x - 307*x2 + 555*x3)*H1*zeta3)/(9.*x)
                   + (8*(1 + 2*x + 2*x2)*(10*Hm14 - 16*Hm221 - 24*Hm1m30 - 24*Hm1m22 - 11*Hm1m13 + 10*Hm122 - 10*Hm131
                                          + 8*Hm1m2m10 - 16*Hm1m200 + 2*Hm1m1m20 + 18*Hm1m1m12 + 16*Hm1m121 + 4*Hm1211
                                          + 7*Hm1m1m100 - 3*Hm1m1000 + 13*Hm10000 + 28*Hm1m2*zeta2 - 18*Hm1m1m1*zeta2
                                          + 18*Hm1m10*zeta2 + 4*Hm1m1*zeta3))/3. - (4*(83 + 166*x + 163*x2)*Hm1*zeta4)/3.
                   + ((-425 - 38*x + 136*x2)*H0*zeta4)/3.
                   + (-4*(144 + 426*x + 1830*x2 + 3167*x3 + 1756*x4 + 155*x5)*H200
                      + 12*(16 + 103*x + 398*x2 + 671*x3 + 456*x4 + 102*x5)*H2*zeta2
                      + 4*(48 - 42*x - 1470*x2 - 811*x3 + 2308*x4 + 1757*x5)*H0*zeta3
                      + (1840 + 5012*x + 17187*x2 + 31855*x3 + 22925*x4 + 5103*x5)*zeta4)/(9.*x*opx2)
                   - (2*(1 - 2*x + 2*x2)*(4*H14 - 56*H1m30 + 16*H1m22 - 64*H113 - 28*H122 - 76*H130 - 20*H131 - 16*H1m2m10
                                          + 140*H1112 - 8*H1120 + 80*H1121 - 116*H1200 - 4*H1210 + 40*H1211 + 44*H10000 + 76*H11000
                                          - 52*H11100 + 172*H11110 - 120*H11111 - 24*H1m2*zeta2 - 4*H12*zeta2 - 44*H100*zeta2
                                          + 32*H110*zeta2 + 16*H111*zeta2 - 120*H10*zeta3 + 172*H11*zeta3 + 57*H1*zeta4))/3.
                   - (4*(77 - 1206*x + 136*x2)*zeta5)/3.)
         + CA*CF*TR*((-116208 - 834293*x + 59486*x2 + 803580*x3)/(1944.*x) + ((-13165 - 32212*x - 26970*x2)*H0)/243.
                     + (2*(-2285 - 60229*x + 15479*x2 + 48075*x3)*H1)/(243.*x) - (8*(2576 + 3701*x + 2672*x2)*H2)/81.
                     + (4*(134 - 278*x + 179*x2)*H3)/9. - (4*(-6 + 78*x + 631*x2)*H4)/9. + (16*(11 + 16*x + 8*x2)*H5)/3.
                     - (16*(25 - 14*x + 28*x2)*Hm40)/3. - (8*(139 + 84*x + 63*x2)*Hm30)/3. + (16*(23 + 10*x + 44*x2)*Hm32)/3.
                     + (4*(-339 + 194*x + 70*x2)*Hm20)/3. + 8*(15 + 28*x + 16*x2)*Hm22 + (16*(23 + 30*x + 47*x2)*Hm23)/3.
                     + (2*(-929 - 824*x + 125*x2)*Hm10)/3. + (8*(55 + 82*x + 28*x2)*Hm12)/3. + 8*(1 + x)*(10 + 13*x)*Hm13
                     + ((1745 + 26249*x - 47946*x2)*H00)/81. + (4*(-2 - 881*x + 694*x2)*H21)/27. - (8*(-18 - 7*x + 14*x2)*H22)/3.
                     + (16*(9 - 40*x + 3*x2)*H23)/3. + (4*(3 + 186*x + 172*x2)*H31)/9. + (16*(1 - 14*x + 4*x2)*H32)/3.
                     + (64*(1 - x)*H40)/3. + (8*(-3 - 10*x + 4*x2)*H41)/3. - (16*(1 - 10*x + 10*x2)*Hm3m10)/3.
                     + (16*(4 + 8*x + 25*x2)*Hm300)/3. + (16*(3 - 22*x + 8*x2)*Hm2m20)/3. + 8*(6 + 20*x + x2)*Hm2m10
                     - (32*(17 + 14*x + 33*x2)*Hm2m12)/3. - (4*(79 + 45*x2)*Hm200)/3. + (16*(3 + 22*x + 4*x2)*Hm220)/3.
                     + (128*x*Hm221)/3. + (8*(53 + 92*x + 27*x2)*Hm1m20)/3. - (8*(36 + 64*x + 35*x2)*Hm1m10)/3.
                     - 16*(1 + x)*(11 + 8*x)*Hm1m12 + (4*(-133 - 98*x + 27*x2)*Hm100)/3. + (8*(-14 - 18*x + 5*x2)*Hm120)/3.
                     + (16*(-7 + x)*(1 + x)*Hm121)/3. - 16*(1 + 2*x + 3*x2)*Hm130 + ((299 + 11222*x + 1836*x2)*H000)/27.
                     - (32*(22 - 34*x + 9*x2)*H1m20)/3.
                     + (4*(-16 + 1056*x - 1290*x2 + 349*x3)*H12 + 8*(46 + 20*x + 74*x2 + 21*x3)*H111)/(27.*x)
                     - (32*(5 - 4*x + 2*x2)*H2m20)/3. + (4*(-37 + 128*x + 176*x2)*H211)/9. + (8*(35 - 34*x + 76*x2)*H212)/3.
                     - (32*(1 + 7*x + 2*x2)*H220)/3. + (8*(11 - 46*x + 18*x2)*H221)/3. - (8*(17 - 22*x + 10*x2)*H300)/3.
                     - (16*(9 + 4*x + 31*x2)*H310)/3. - (32*(2 + 3*x + 2*x2)*H311)/3. - (16*(5 + 2*x + 18*x2)*Hm2m1m10)/3.
                     - (8*(27 + 70*x + 64*x2)*Hm2m100)/3. + (8*(37 + 70*x + 90*x2)*Hm2000)/3. - 8*(10 + 16*x + x2)*Hm1m1m10
                     + (4*(8 + 5*x)*(10 + 9*x)*Hm1m100)/3. - (4*(-9 - 26*x + 13*x2)*Hm1000)/3. - (32*(2 + 4*x + 7*x2)*Hm1200)/3.
                     + 32*(2 + 4*x + 5*x2)*Hm1210 - (4*(121 - 230*x + 948*x2)*H0000)/9. + (16*(-2 - 24*x + 7*x2)*H2000)/3.
                     - (16*(-3 - 14*x + 7*x2)*H2100)/3. + (16*(10 - 14*x + 47*x2)*H2110)/3. - 8*(11 - 18*x + 20*x2)*H2111
                     + (32*(3 + 8*x)*H00000)/3. - (8*(47 + 10*x + 98*x2)*Hm3*zeta2)/3. - 4*(24 + 36*x + 31*x2)*Hm2*zeta2
                     - (4*(146 + 228*x + 91*x2)*Hm1*zeta2)/3. - 16*(5 - 10*x + 17*x2)*H3*zeta2 + 8*(21 + 18*x + 38*x2)*Hm2m1*zeta2
                     - (16*(31 + 30*x + 54*x2)*Hm20*zeta2)/3. + 4*(34 + 60*x + 31*x2)*Hm1m1*zeta2
                     - (4*(143 + 260*x + 99*x2)*Hm10*zeta2)/3. + 16*(3 + 6*x + 8*x2)*Hm12*zeta2
                     - (8*(25 - 92*x + 8*x2)*H20*zeta2)/3. - (32*(6 - 11*x + x2)*H21*zeta2)/3.
                     - (4*(49 + 98*x + 110*x2)*Hm100*zeta2)/3. - (4*(47 + 62*x + 32*x2)*H000*zeta2)/3.
                     + (-2*(-404 - 1007*x - 5051*x2 + 6190*x3)*H10 + 4*(661 - 1894*x + 905*x2 + 564*x3)*H11
                        + 2*(-1812 + 13193*x - 10027*x2 + 13814*x3)*zeta2)/(81.*x)
                     + (2*(708 + 180*x + 1765*x2 + 2239*x3)*H20 - 6*(310 - 1272*x - 1097*x2 + 503*x3)*H0*zeta2)/(27.*(1 + x))
                     + ((1 - x)*(-96*(2 + 17*x2)*H121 + 4*(16 + 49*x + 25*x2)*H11*zeta2))/(9.*x)
                     - (16*(23 + 18*x + 42*x2)*Hm2*zeta3)/3. - (8*(32 + 66*x + 37*x2)*Hm1*zeta3)/3.
                     + (8*(-31 - 38*x + 20*x2)*H2*zeta3)/3. - (8*(41 + 82*x + 70*x2)*Hm10*zeta3)/3.
                     + (8*(-33 - 66*x + 4*x2)*H00*zeta3)/3. + (8*(77 - 9*x + 128*x2)*zeta2*zeta3)/3.
                     + (-4*(92 + 736*x - 575*x2 - 563*x3 + 710*x4)*H100 + 2*(-224 + 655*x - 1863*x2 - 763*x3 + 2087*x4)*H110
                        - 2*(-32 + 1999*x + 432*x2 - 1324*x3 + 167*x4)*H1*zeta2
                        - 2*(-1536 - 5458*x - 15548*x2 - 4303*x3 + 7215*x4)*zeta3)/(27.*x*(1 + x))
                     + (-4*(84 - 546*x - 758*x2 + 494*x3 + 613*x4)*H30 - 12*(72 + 40*x + 45*x2 + 282*x3 + 199*x4)*H200
                        + 12*(-11 - 48*x - 24*x2 + 76*x3 + 57*x4)*H210 + 12*(-39 - 104*x - 72*x2 + 36*x3 + 37*x4)*H2*zeta2
                        + 2*(-3 + 438*x + 1847*x2 + 2296*x3 + 908*x4)*H00*zeta2
                        + 4*(255 + 2046*x + 4405*x2 + 3764*x3 + 1132*x4)*H0*zeta3)/(9.*opx2)
                     + (4*(-32 - 81*x - 72*x2 + 329*x3)*H13 + 8*(-8 + 69*x - 87*x2 + 17*x3)*H112
                        - 4*(16 + 114*x - 204*x2 + 173*x3)*H120 - 4*(8 + 239*x - 340*x2 + 203*x3)*H1000
                        - 4*(-16 + 174*x - 522*x2 + 373*x3)*H1100 + 4*(-16 + 273*x - 552*x2 + 277*x3)*H1110
                        + 8*(-4 - 29*x - 23*x2 + 24*x3)*H1111 - 4*(-44 + 177*x - 564*x2 + 683*x3)*H10*zeta2
                        - 24*(8 + 97*x - 187*x2 + 88*x3)*H1*zeta3)/(9.*x)
                     + (8*(1 + 2*x + 2*x2)*(20*Hm14 + 32*Hm1m30 - 60*Hm1m22 - 53*Hm1m13 + 10*Hm122 - 6*Hm131 - 24*Hm1m2m10
                                            - 22*Hm1m200 - 10*Hm1m1m20 + 78*Hm1m1m12 - 4*Hm1m120 - 4*Hm1211 + 32*Hm1m1m1m10
                                            + 41*Hm1m1m100 - 51*Hm1m1000 + 15*Hm10000 + 48*Hm1m2*zeta2 - 62*Hm1m1m1*zeta2
                                            + 61*Hm1m10*zeta2 + 44*Hm1m1*zeta3))/3.
                     + ((312 - 2291*x - 10332*x2 - 15038*x3 - 9732*x4 - 2441*x5)*zeta4)/(9.*x*opx2)
                     + (8*(59 + 118*x + 121*x2)*Hm1*zeta4)/3. - (2*(167 - 130*x + 148*x2)*H0*zeta4)/3.
                     + (2*(1 - 2*x + 2*x2)*(44*H14 - 8*H1m30 - 16*H1m22 + 24*H113 + 44*H122 - 52*H130 + 52*H131
                                            - 80*H1m2m10 + 48*H1m200 - 32*H11m20 + 132*H1112 - 24*H1120 + 24*H1121
                                            - 92*H1200 + 12*H1210 - 72*H1211 - 36*H10000 + 68*H11000 - 28*H11100
                                            + 180*H11110 - 240*H11111 - 24*H1m2*zeta2 - 74*H12*zeta2 - 68*H100*zeta2
                                            - 38*H110*zeta2 + 46*H111*zeta2 - 80*H10*zeta3 - 20*H11*zeta3 + H1*zeta4))/3.
                     + (4*(-49 - 974*x + 304*x2)*zeta5)/3.));
  }

  //_________________________________________________________________________________
  C1gqpdfBM::C1gqpdfBM():
    Expression()
  {
  }
  double C1gqpdfBM::Regular(double const& x) const
  {
    return - 4 * CF * ( 1 - x ) / x;
  }

  //_________________________________________________________________________________
  C1ggpdfBM::C1ggpdfBM():
    Expression()
  {
  }
  double C1ggpdfBM::Regular(double const& x) const
  {
    return - 4 * CA * ( 1 - x ) / x;
  }
}
