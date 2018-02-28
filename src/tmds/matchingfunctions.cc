//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/matchingfunctions.h"
#include "apfel/constants.h"

#include <numeric>

namespace apfel
{
  //_________________________________________________________________________________
  C1ns::C1ns():
    Expression()
  {
  }
  double C1ns::Regular(double const& x) const
  {
    return 2 * CF * ( 1 - x );
  }
  double C1ns::Local(double const&) const
  {
    return - CF * zeta2;
  }
  
  //_________________________________________________________________________________
  C1qg::C1qg():
    Expression()
  {
  }
  double C1qg::Regular(double const& x) const
  {
    return 8 * TR * x * ( 1 - x );
  }

  //_________________________________________________________________________________
  C1gq::C1gq():
    Expression()
  {
  }
  double C1gq::Regular(double const& x) const
  {
    return 2 * CF * x;
  }

  //_________________________________________________________________________________
  C1gg::C1gg():
    Expression()
  {
  }
  double C1gg::Local(double const&) const
  {
    return - CA * zeta2;
  }

  //_________________________________________________________________________________
  C2Vqq::C2Vqq(int const& nf):
    Expression(),
    _nf(nf)
  {
    _A2 = - 3232. / 27. + 112 * zeta3 + 448. * _nf / 81.;
    _A3 = 0.;
  }
  double C2Vqq::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx2 * lx;
    const double l1x  = log(1-x);
    const double l1x2 = l1x * l1x;
    const double l1x3 = l1x2 * l1x;
    const std::vector<double> fReg{
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

    const std::vector<double> CoeffReg{
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
  double C2Vqq::Singular(double const& x) const
  {
    const double fA2 = 1 / ( 1 - x );
    const double fA3 = fA2 * log(1-x);

    return _A2 * fA2 + _A3 * fA3;
  }
  double C2Vqq::Local(double const& x) const
  {
    const double ln1mx  = log(1-x);
    const double ln1mx2 = ln1mx * ln1mx;
    const double A1     = - 2416. / 81. - 134 * zeta2 / 3 + 448 * zeta3 / 9 + 200 * zeta4 / 9
      + _nf * ( 352. / 243. + 20 * zeta2 / 9 + 56 * zeta3 / 27 );

    return A1 + _A2 * ln1mx + _A3 * ln1mx2 / 2;
  }

  //_________________________________________________________________________________
  C2Vqqb::C2Vqqb()
  {
  }
  double C2Vqqb::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx2 * lx;
    const double l1x  = log(1-x);
    const double l1x2 = l1x * l1x;
    const double l1x3 = l1x2 * l1x;
    const std::vector<double> fReg{
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

    const std::vector<double> CoeffReg{
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
  C2ps::C2ps():
    Expression()
  {
  }
  double C2ps::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx2 * lx;
    const double l1x  = log(1-x);
    const double l1x2 = l1x * l1x;
    const double l1x3 = l1x2 * l1x;
    const std::vector<double> fReg{
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

    const std::vector<double> CoeffReg{
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
  C2qg::C2qg():
    Expression()
  {
  }
  double C2qg::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx2 * lx;
    const double l1x  = log(1-x);
    const double l1x2 = l1x * l1x;
    const double l1x3 = l1x2 * l1x;
    const std::vector<double> fReg{
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

    const std::vector<double> CoeffReg{
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
  C2gq::C2gq(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C2gq::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx2 * lx;
    const double l1x  = log(1-x);
    const double l1x2 = l1x * l1x;
    const double l1x3 = l1x2 * l1x;
    const std::vector<double> fReg{
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

    const std::vector<double> CoeffReg{
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
  C2gg::C2gg(int const& nf):
    Expression(),
    _nf(nf)
  {
    _A2 = - 808. / 3. + 252. * zeta3 + 112. * _nf / 9.;
    _A3 = 0.;
  }
  double C2gg::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx2 * lx;
    const double l1x  = log(1-x);
    const double l1x2 = l1x * l1x;
    const double l1x3 = l1x2 * l1x;
    const std::vector<double> fReg{
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

    const std::vector<double> CoeffReg{
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
  double C2gg::Singular(double const& x) const
  {
    const double fA2 = 1 / ( 1 - x );
    const double fA3 = fA2 * log(1-x);

    return _A2 * fA2 + _A3 * fA3;
  }
  double C2gg::Local(double const& x) const
  {
    const double ln1mx  = log(1-x);
    const double ln1mx2 = ln1mx * ln1mx;
    const double A1     = - 112. - 201. * zeta2 / 2. + 154. * zeta3 + 225. * zeta4 / 4.
      + _nf * ( 548. / 27. + 5. * zeta2 - 28. * zeta3 / 3. )
      - 56. * _nf * _nf  / 81.;

    return A1 + _A2 * ln1mx + _A3 * ln1mx2 / 2;
  }

  //_________________________________________________________________________________
  C1nsff::C1nsff():
    Expression()
  {
  }
  double C1nsff::Regular(double const& x) const
  {
    return CF * ( 2 * ( 1 - x ) + 4 * ( 1 + x * x ) * log(x) / ( 1 - x ) );
  }
  double C1nsff::Local(double const&) const
  {
    return - CF * zeta2;
  }

  //_________________________________________________________________________________
  C1qgff::C1qgff():
    Expression()
  {
  }
  double C1qgff::Regular(double const& x) const
  {
    return 2 * CF * ( 2 * x + 4 * ( 1 + pow(1 - x, 2) ) * log(x) / x );
  }

  //_________________________________________________________________________________
  C1gqff::C1gqff():
    Expression()
  {
  }
  double C1gqff::Regular(double const& x) const
  {
    return TR * ( 4 * x * ( 1 - x ) + 4 * ( 1 - 2 * x * ( 1 - x ) ) * log(x) );
  }

  //_________________________________________________________________________________
  C1ggff::C1ggff():
    Expression()
  {
  }
  double C1ggff::Regular(double const& x) const
  {
    return 8 * CA * pow(1 - x * ( 1 - x ), 2) * log(x) / x / ( 1 - x );
  }
  double C1ggff::Local(double const&) const
  {
    return - CA * zeta2;
  }

  //_________________________________________________________________________________
  C2Vqqff::C2Vqqff(int const& nf):
    Expression(),
    _nf(nf)
  {
    _A2 = - 3232. / 27. + 112 * zeta3 + 448. * _nf / 81.;
    _A3 = 0.;
  }
  double C2Vqqff::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx2 * lx;
    const double l1x  = log(1-x);
    const double l1x2 = l1x * l1x;
    const double l1x3 = l1x2 * l1x;
    const std::vector<double> fReg{
      l1x,
	l1x2,
	l1x3,
	1 / x,
	lx / x,
	lx2 / x,
	lx3 / x,
	lx,
	lx2,
	lx3,
	1,
	x,
	x2,
	x * lx / ( 1 - x ),
	x * lx,
	x2 * lx,
	x * lx2 / ( 1 - x ),
	x * lx2,
	x * lx3,
	( lx / ( 1 - x ) + 1 ) * l1x,
	lx * l1x,
	x * lx * l1x,
	( 1 - x ) * l1x / x,
	( 1 - x ) * l1x,
	( 1 - x ) * l1x2};

    const std::vector<double> CoeffReg{
      - 200. / 9.,
	64. / 9.,
	0.,
	0.,
	0.,
	0.,
	0.,
	1496. / 9. - 8. * _nf - 560. * zeta2 / 9.,
	- 130. / 9. + 4. * _nf / 9.,
	- 140. / 27.,
	- 296. * _nf / 81. - 301.03247439776976,
	- 152. * _nf / 81. - 989.2167272286393,
	82.59056818996338,
	- 80. * _nf / 9. - 1063.98482846164,
	8. * _nf + 206.28577290245227,
	-18.651271690975136,
	8. * _nf / 9. + 83.00296625888389,
	- 4. * _nf / 9. - 70.52319745715631,
	- 4.911975080877064,
	- 1105.7693500845382,
	327.376932797376,
	- 109.45363459015105,
	- 174.6471655693391,
	- 112.83673919345797,
	- 3.5294575557396084};

    return std::inner_product(fReg.begin(), fReg.end(), CoeffReg.begin(), 0.);
  }
  double C2Vqqff::Singular(double const& x) const
  {
    const double fA2 = 1 / ( 1 - x );
    const double fA3 = fA2 * log(1-x);

    return _A2 * fA2 + _A3 * fA3;
  }
  double C2Vqqff::Local(double const& x) const
  {
    const double ln1mx  = log(1-x);
    const double ln1mx2 = ln1mx * ln1mx;
    const double A1     = - 2416. / 81. - 134 * zeta2 / 3 + 448 * zeta3 / 9 + 2360 * zeta4 / 9
      + _nf * ( 352. / 243. + 20. * zeta2 / 9 + 56 * zeta3 / 27 );

    return A1 + _A2 * ln1mx + _A3 * ln1mx2 / 2;
  }

  //_________________________________________________________________________________
  C2Vqqbff::C2Vqqbff()
  {
  }
  double C2Vqqbff::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx2 * lx;
    const double l1x  = log(1-x);
    const double l1x2 = l1x * l1x;
    const double l1x3 = l1x2 * l1x;
    const std::vector<double> fReg{
      l1x,
	l1x2,
	l1x3,
	1 / x,
	lx / x,
	lx2 / x,
	lx3 / x,
	lx,
	lx2,
	lx3,
	1,
	x,
	x2,
	x * lx / ( 1 - x ),
	x * lx,
	x2 * lx,
	x * lx2 / ( 1 - x ),
	x * lx2,
	x * lx3,
	( lx / ( 1 - x ) + 1 ) * l1x,
	lx * l1x,
	x * lx * l1x,
	( 1 - x ) * l1x / x,
	( 1 - x ) * l1x,
	( 1 - x ) * l1x2};

    const std::vector<double> CoeffReg{
      0.,
	0.,
	0.,
	0.,
	0.,
	0.,
	0.,
	- 76. / 9.,
	- 32. / 9.,
	- 4. / 3.,
	- 448.1996646927744,
	- 5647.490349524448,
	- 257.33492226817515,
	- 6353.024936485396,
	3394.3630072013384,
	84.62666948172131,
	483.57763340201967,
	- 518.3011843768235,
	- 0.3446236997186985,
	- 6484.394475423013,
	3796.7810643424796,
	339.26036997443003,
	- 436.86534933603923,
	1330.4397468097895,
	- 0.02218996661410778};

    return std::inner_product(fReg.begin(), fReg.end(), CoeffReg.begin(), 0.);
  }

  //_________________________________________________________________________________
  C2psff::C2psff():
    Expression()
  {
  }
  double C2psff::Regular(double const& x) const
  {    const double x2   = x * x;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx2 * lx;
    const double l1x  = log(1-x);
    const double l1x2 = l1x * l1x;
    const double l1x3 = l1x2 * l1x;
    const std::vector<double> fReg{
      l1x,
	l1x2,
	l1x3,
	1 / x,
	lx / x,
	lx2 / x,
	lx3 / x,
	lx,
	lx2,
	lx3,
	1,
	x,
	x2,
	x * lx / ( 1 - x ),
	x * lx,
	x2 * lx,
	x * lx2 / ( 1 - x ),
	x * lx2,
	x * lx3,
	( lx / ( 1 - x ) + 1 ) * l1x,
	lx * l1x,
	x * lx * l1x,
	( 1 - x ) * l1x / x,
	( 1 - x ) * l1x,
	( 1 - x ) * l1x2};

    const std::vector<double> CoeffReg{
      0.,
	0.,
	0.,
	- 592. / 81. + 32. * zeta2 / 9.,
	32. / 9.,
	64. / 9.,
	0.,
	- 48.,
	- 22. / 3.,
	44. / 9.,
	- 28.813571629909163,
	206.17553030550255,
	76.0334481394452,
	251.93541929963473,
	- 169.05906218222225,
	- 32.29013619101719,
	- 10.685799944808078,
	3.7898590887852626,
	4.909581801691148,
	188.58291934528984,
	- 90.34188300607897,
	- 1.4634823045099683,
	18.62607672661471,
	- 16.127886782439663,
	0.0009797967055855182};

    return 2 * std::inner_product(fReg.begin(), fReg.end(), CoeffReg.begin(), 0.);
  }

  //_________________________________________________________________________________
  C2qgff::C2qgff():
    Expression()
  {
  }
  double C2qgff::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx2 * lx;
    const double l1x  = log(1-x);
    const double l1x2 = l1x * l1x;
    const double l1x3 = l1x2 * l1x;
    const std::vector<double> fReg{
      l1x,
	l1x2,
	l1x3,
	1 / x,
	lx / x,
	lx2 / x,
	lx3 / x,
	lx,
	lx2,
	lx3,
	1,
	x,
	x2,
	x * lx / ( 1 - x ),
	x * lx,
	x2 * lx,
	x * lx2 / ( 1 - x ),
	x * lx2,
	x * lx3,
	( lx / ( 1 - x ) + 1 ) * l1x,
	lx * l1x,
	x * lx * l1x,
	( 1 - x ) * l1x / x,
	( 1 - x ) * l1x,
	( 1 - x ) * l1x2};

    const std::vector<double> CoeffReg{
      - 40./9. - 80. * zeta2 / 3.,
	- 40. /9.,
	40. / 27.,
	12512. / 27. - 352. * zeta2 / 3. - 352. * zeta3,
	400. / 3.,
	- 848. / 3.,
	- 320. / 3.,
	1144. / 3.,
	64. / 3.,
	- 1232. / 27.,
	- 11519.346897372414,
	- 34241.466106099186,
	5326.770104932414,
	- 40601.87518176106,
	4178.463030904903,
	- 1705.6350033291087,
	- 966.5754106411847,
	1144.6267544753136,
	- 84.66732541780037,
	- 11393.035115581788,
	- 25857.49295712562,
	- 10601.55795204891,
	- 12368.214954397781,
	- 29991.60756399795,
	- 5.282535747460972};

    return 2 * std::inner_product(fReg.begin(), fReg.end(), CoeffReg.begin(), 0.);
  }

  //_________________________________________________________________________________
  C2gqff::C2gqff(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C2gqff::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx2 * lx;
    const double l1x  = log(1-x);
    const double l1x2 = l1x * l1x;
    const double l1x3 = l1x2 * l1x;
    const std::vector<double> fReg{
      l1x,
	l1x2,
	l1x3,
	1 / x,
	lx / x,
	lx2 / x,
	lx3 / x,
	lx,
	lx2,
	lx3,
	1,
	x,
	x2,
	x * lx / ( 1 - x ),
	x * lx,
	x2 * lx,
	x * lx2 / ( 1 - x ),
	x * lx2,
	x * lx3,
	( lx / ( 1 - x ) + 1 ) * l1x,
	lx * l1x,
	x * lx * l1x,
	( 1 - x ) * l1x / x,
	( 1 - x ) * l1x,
	( 1 - x ) * l1x2};

    const std::vector<double> CoeffReg{
      - 44. / 3. + 10. * _nf / 9. + 10. * zeta2,
	- 7. / 2. + _nf / 3.,
	- 5. / 9.,
	- 148. / 9. + 8. * zeta2,
	8.,
	16.,
	0.,
	- 157. / 3. - 10. * _nf / 3. - 100. * zeta2 / 3.,
	- 73. / 3. + _nf / 3.,
	77. / 9.,
	- 5042.153939120427 - 66.62886623235278 * _nf,
	303.86339171965324 - 317.3462434009907 * _nf,
	3703.409480057778 + 21.346846003738534 * _nf,
	- 1059.0639607339292 - 361.41246956998253 * _nf,
	- 8193.44131612091 + 116.14537890481989 * _nf,
	- 1046.4745476064809 - 18.841069748305667 * _nf,
	- 1458.9834262434276 + 5.214301688311585 * _nf,
	1590.9730212100676 - 6.201059346425027 * _nf,
	73.07631646309086 - 0.014948523947447324 * _nf,
	10420.770382128236 - 200.7638366136303 * _nf,
	- 20512.88777168424 - 49.38978687382705 * _nf,
	- 4683.056024325436 - 58.30447555159878 * _nf,
	- 4929.425576232434 - 65.41307217273041 * _nf,
	- 15198.304463044708 - 144.13788267277235 * _nf,
	- 11.384648952243019 - 0.6611076391799198 * _nf};

    return std::inner_product(fReg.begin(), fReg.end(), CoeffReg.begin(), 0.);
  }

  //_________________________________________________________________________________
  C2ggff::C2ggff(int const& nf):
    Expression(),
    _nf(nf)
  {
    _A2 = - 808. / 3. + 252 * zeta3 + 112. * _nf / 9.;
    _A3 = 0.;
  }
  double C2ggff::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx2 * lx;
    const double l1x  = log(1-x);
    const double l1x2 = l1x * l1x;
    const double l1x3 = l1x2 * l1x;
    const std::vector<double> fReg{
      l1x,
	l1x2,
	l1x3,
	1 / x,
	lx / x,
	lx2 / x,
	lx3 / x,
	lx,
	lx2,
	lx3,
	1,
	x,
	x2,
	x * lx / ( 1 - x ),
	x * lx,
	x2 * lx,
	x * lx2 / ( 1 - x ),
	x * lx2,
	x * lx3,
	( lx / ( 1 - x ) + 1 ) * l1x,
	lx * l1x,
	x * lx * l1x,
	( 1 - x ) * l1x / x,
	( 1 - x ) * l1x,
	( 1 - x ) * l1x2};

    const std::vector<double> CoeffReg{
      - 6. + 2. * _nf,
	36.,
	0.,
	3268. / 3. + 278. * _nf / 81. - 264. * zeta2 - 792. * zeta3,
	536. - 244. * _nf / 9. - 360. * zeta2,
	- 660. - 8. * _nf / 9.,
	- 240.,
	819. - 50. * _nf + 360. * zeta2,
	33. + 22. * _nf / 3.,
	- 132. + 88. * _nf / 9.,
	- 4536.42238065816 - 163.95224892156722 * _nf,
	18080.936974539778 - 625.9193713983447 * _nf,
	14139.810571485807 - 184.19497231902477 * _nf,
	27937.822349639453 - 975.5233827623935 * _nf,
	- 35293.57801079519 + 498.8898228158924 * _nf,
	- 5939.483006875214 + 81.15085632871367 * _nf,
	- 6775.737946774227 + 30.343871838863492 * _nf,
	7209.990192354363 - 51.50813149633906 * _nf,
	- 720.1068548927861 + 9.691339026103627 * _nf,
	62244.278348304506 - 725.7040599397907 * _nf,
	- 73593.14246325135 + 340.92281843959506 * _nf,
	- 19274.85766344364 + 2.9241023778057107 * _nf,
	- 6776.731748013898 - 73.730026699345 * _nf,
	- 54745.95834168829 + 54.675824792462755 * _nf,
	- 17.510572839038723 - 0.003977517619150036 * _nf};

    return std::inner_product(fReg.begin(), fReg.end(), CoeffReg.begin(), 0.);
  }
  double C2ggff::Singular(double const& x) const
  {
    const double fA2 = 1 / ( 1 - x );
    const double fA3 = fA2 * log(1-x);

    return _A2 * fA2 + _A3 * fA3;
  }
  double C2ggff::Local(double const& x) const
  {
    const double ln1mx  = log(1-x);
    const double ln1mx2 = ln1mx * ln1mx;
    const double A1     = - 112. - 201 * zeta2 / 2 + 154 * zeta3 + 2385 * zeta4 / 4
      + _nf * ( 548. / 27. + 5 * zeta2 - 28 * zeta3 / 3 )
      - 56.  *_nf * _nf / 81.;

    return A1 + _A2 * ln1mx + _A3 * ln1mx2 / 2;
  }
}
