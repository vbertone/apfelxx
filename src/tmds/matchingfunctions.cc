//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/matchingfunctions.h"
#include "apfel/tools.h"

#include <cmath>
#include <numeric>

using std::inner_product;

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
    const double fA4 = log(1-x);
    const double fA5 = fA4 * fA4;
    const double fA6 = fA5 * fA4;

    const double fB1 = log(x);
    const double fB2 = fB1 * fB1;
    const double fB3 = fB2 * fB1;
    const double fB4 = 1 / x;
    const double fB5 = fB1 * fB4;

    const double fc1 = 1;
    const double fc2 = x;
    const double fc3 = fc2 * fc2;
    const double fc4 = fc3 * fc2;
    const double fc5 = fA4 * fB1;
    const double fc6 = fc5 * fB1;

    const vector<double> fReg{fc1, fc2, fc3, fc4, fA4, fA5, fA6, fB1, fB2, fB3, fc5, fc6, fB4, fB5};
    const vector<double> CoeffReg{- 15.591812729537994 - 3.6403764478352265 * _nf,
	- 26.19077732354684 - 4.035868117241033 * _nf,
	61.362896852454675 - 0.6181383623814728 * _nf,
	- 18.603456697689623 - 0.19612093929051017 * _nf,
	200. / 9.,
	- 64. / 9.,
	0.,
	- 8. + 40. * _nf / 27.,
	- 2. + 4. * _nf / 9.,
	- 20. / 27.,
	32.70163574653397 - 0.069575874398736 * _nf,
	12.62916203247169 - 0.17334865938082786 * _nf,
	0.,
	0.};

    return inner_product(fReg.begin(), fReg.end(), CoeffReg.begin(), 0.);
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
    const double A1     = - 2416. / 81. - 67 * Pi2 / 9 + 448 * zeta3 / 9 + 20 * Pi2 * Pi2 / 81
      + _nf * ( 352. / 243. + 10 * Pi2 / 27 + 56 * zeta3 / 27 );

    return A1 + _A2 * ln1mx + _A3 * ln1mx2 / 2;
  }

  //_________________________________________________________________________________
  C2Vqqb::C2Vqqb()
  {
  }
  double C2Vqqb::Regular(double const& x) const
  {
    const double fA4 = log(1-x);
    const double fA5 = fA4 * fA4;
    const double fA6 = fA5 * fA4;

    const double fB1 = log(x);
    const double fB2 = fB1 * fB1;
    const double fB3 = fB2 * fB1;

    const double fc1 = 1;
    const double fc2 = x;
    const double fc3 = fc2 * fc2;
    const double fc4 = fc3 * fc2;
    const double fc5 = fA4 * fB1;
    const double fc6 = fc5 * fB1;

    const vector<double> fReg{fc1, fc2, fc3, fc4, fA4, fA5, fA6, fB1, fB2, fB3, fc5, fc6};

    const vector<double> CoeffReg{- 3.329144380448291,
	5.852825050121677,
	- 4.796244777535505,
	2.2557153775946435,
	0.,
	0.,
	0.,
	- 4. / 3.,
	0.,
	4. / 27.,
	0.505971782717726,
	0.08759292525297521};

    return inner_product(fReg.begin(), fReg.end(), CoeffReg.begin(), 0.);
  }

  //_________________________________________________________________________________
  C2ps::C2ps():
    Expression()
  {
  }
  double C2ps::Regular(double const& x) const
  {
    const double fA4 = log(1-x);
    const double fA5 = fA4 * fA4;
    const double fA6 = fA5 * fA4;

    const double fB1 = log(x);
    const double fB2 = fB1 * fB1;
    const double fB3 = fB2 * fB1;
    const double fB4 = 1 / x;
    const double fB5 = fB1 * fB4;

    const double fc1 = 1;
    const double fc2 = x;
    const double fc3 = fc2 * fc2;
    const double fc4 = fc3 * fc2;
    const double fc5 = fA4 * fB1;
    const double fc6 = fc5 * fB1;

    const vector<double> fReg{fc1, fc2, fc3, fc4, fA4, fA5, fA6, fB1, fB2, fB3, fc5, fc6, fB4, fB5};
    const vector<double> CoeffReg{- 3.54792360390173,
	4.949994905325898,
	- 11.357230229521692,
	7.2245197835090105,
	0.,
	0.,
	0.,
	8. / 3.,
	- 2. / 3.,
	4. / 9.,
	2.255041649827519,
	0.7663091186642608,
	688. / 81. - 32 * zeta2 / 9,
	0.};

    return inner_product(fReg.begin(), fReg.end(), CoeffReg.begin(), 0.);
  }

  //_________________________________________________________________________________
  C2qg::C2qg():
    Expression()
  {
  }
  double C2qg::Regular(double const& x) const
  {
    const double fA4 = log(1-x);
    const double fA5 = fA4 * fA4;
    const double fA6 = fA5 * fA4;

    const double fB1 = log(x);
    const double fB2 = fB1 * fB1;
    const double fB3 = fB2 * fB1;
    const double fB4 = 1 / x;
    const double fB5 = fB1 * fB4;

    const double fc1 = 1;
    const double fc2 = x;
    const double fc3 = fc2 * fc2;
    const double fc4 = fc3 * fc2;
    const double fc5 = fA4 * fB1;
    const double fc6 = fc5 * fB1;

    const vector<double> fReg{fc1, fc2, fc3, fc4, fA4, fA5, fA6, fB1, fB2, fB3, fc5, fc6, fB4, fB5};
    const vector<double> CoeffReg{5.32456247529479,
	- 37.8039790325823,
	55.939541231907775,
	- 44.592546819513025,
	- 5. / 3.,
	0.,
	5. / 9.,
	34. / 3.,
	- 7. / 6.,
	7. / 9.,
	- 10.961736574029366,
	- 6.522774793663304,
	172. / 9. - 4 * Pi2 / 3,
	0.};

    return 2 * inner_product(fReg.begin(), fReg.end(), CoeffReg.begin(), 0.);
  }

  //_________________________________________________________________________________
  C2gq::C2gq(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C2gq::Regular(double const& x) const
  {
    const double fA4 = log(1-x);
    const double fA5 = fA4 * fA4;
    const double fA6 = fA5 * fA4;

    const double fB1 = log(x);
    const double fB2 = fB1 * fB1;
    const double fB3 = fB2 * fB1;
    const double fB4 = 1 / x;
    const double fB5 = fB1 * fB4;

    const double fc1 = 1;
    const double fc2 = x;
    const double fc3 = fc2 * fc2;
    const double fc4 = fc3 * fc2;
    const double fc5 = fA4 * fB1;
    const double fc6 = fc5 * fB1;

    const vector<double> fReg{fc1, fc2, fc3, fc4, fA4, fA5, fA6, fB1, fB2, fB3, fc5, fc6, fB4, fB5};
    const vector<double> CoeffReg{- 43.40171996161091 - 16.984343881161237 * _nf,
	173.18374973344203 + 9.595722778694926 * _nf,
	- 88.4323955556243 - 2.1972189030825384 * _nf,
	- 8.718852731902624 + 1.1549244843223576 * _nf,
	- 184. / 9. + 32. * _nf / 27.,
	- 44. / 9. + 8. * _nf / 9.,
	- 40. / 27.,
	- 200. / 3., // !!!!!! - 200. / 9. in the paper  !!!!!!!
	112. / 9.,
	- 112. / 27.,
	- 40.28067793284879 - 0.552346106244215 * _nf,
	- 35.36530647197912 - 0.0627184829983808 * _nf,
	- 12640. / 27. + 352 * zeta2 / 3 + 192 * zeta3 + 896. * _nf / 81.,
	  0.};

    return inner_product(fReg.begin(), fReg.end(), CoeffReg.begin(), 0.);
  }

  //_________________________________________________________________________________
  C2gg::C2gg(int const& nf):
    Expression(),
    _nf(nf)
  {
    _A2 = - 808. / 3. + 252 * zeta3 + 112. * _nf / 9.;
    _A3 = 0.;
  }
  double C2gg::Regular(double const& x) const
  {
    const double fA4 = log(1-x);
    const double fA5 = fA4 * fA4;
    const double fA6 = fA5 * fA4;

    const double fB1 = log(x);
    const double fB2 = fB1 * fB1;
    const double fB3 = fB2 * fB1;
    const double fB4 = 1 / x;
    const double fB5 = fB1 * fB4;

    const double fc1 = 1;
    const double fc2 = x;
    const double fc3 = fc2 * fc2;
    const double fc4 = fc3 * fc2;
    const double fc5 = fA4 * fB1;
    const double fc6 = fc5 * fB1;

    const vector<double> fReg{fc1, fc2, fc3, fc4, fA4, fA5, fA6, fB1, fB2, fB3, fc5, fc6, fB4, fB5};
    const vector<double> CoeffReg{160.3174084388 + 1.8798838185664 * _nf,
	235.25958439812186 - 13.729196847420008 * _nf,
	- 77.67188779845978 - 0.6784851744787738 * _nf,
	51.61266721115288 - 1.3169890096131303 * _nf,
	2. - 2. * _nf / 3.,
	- 12.,
	0.,
	- 293. / 3. + 74. * _nf / 9.,
	1. + 2. * _nf,
	- 4. + 8. * _nf / 27.,
	24.137288673188998 - 1.0158393615503878 * _nf,
	- 22.406633639406444 + 1.090527996849969 * _nf,
	- 3160. / 9. + 264 * zeta2 / 3 + 432 * zeta3 / 3 + 226. * _nf / 27.,
	0.};

    return 3 * inner_product(fReg.begin(), fReg.end(), CoeffReg.begin(), 0.);
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
    const double A1     = - 112. - 201 * zeta2 / 2 + 154 * zeta3 + 225 * zeta4 / 4
      + _nf * ( 548. / 27. + 5 * zeta2 - 28 * zeta3 / 3 )
      - 56.  *_nf * _nf / 81.;

    return A1 + _A2 * ln1mx + _A3 * ln1mx2 / 2;
  }
}
