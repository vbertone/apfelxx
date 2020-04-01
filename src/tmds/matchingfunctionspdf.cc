//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/matchingfunctionspdf.h"
#include "apfel/constants.h"

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
    _A2 = - 3232. / 27. + 112 * zeta3 + 448. * _nf / 81.;
  }
  double C2Vqqpdf::Regular(double const& x) const
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
    const double ln1mx = log(1-x);
    const double A1    = - 2416. / 81. - 134 * zeta2 / 3 + 448 * zeta3 / 9 + 200 * zeta4 / 9
                         + _nf * ( 352. / 243. + 20 * zeta2 / 9 + 56 * zeta3 / 27 );

    return A1 + _A2 * ln1mx;
  }

  //_________________________________________________________________________________
  C2Vqqbpdf::C2Vqqbpdf()
  {
  }
  double C2Vqqbpdf::Regular(double const& x) const
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
