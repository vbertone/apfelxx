//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/zeromasscoefficientfunctions_tl.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________________
  C21Tns::C21Tns():
    Expression()
  {
  }
  double C21Tns::Regular(double const& x) const
  {
    return 2 * CF * ( - ( 1 + x ) * log( 1 - x ) + 2 * ( 1 + pow(x,2) ) * log(x) / ( 1 - x ) + 5. / 2. - 3 * x / 2 );
  }
  double C21Tns::Singular(double const& x) const
  {
    return 2 * CF * ( 2 * log( 1 - x ) - 3 / 2. ) / ( 1 - x );
  }
  double C21Tns::Local(double const& x) const
  {
    return 2 * CF * ( pow(log(1-x),2) - 3 * log( 1 - x ) / 2 + ( 4 * zeta2 - 9 / 2. ) );
  }

  //_________________________________________________________________________________
  C21Tg::C21Tg():
    Expression()
  {
  }
  double C21Tg::Regular(double const& x) const
  {
    return 4 * CF * ( ( 1 + pow((1-x),2) ) * log( pow(x, 2) * ( 1 - x ) ) / x );
  }

  //_________________________________________________________________________________
  C22Tnsp::C22Tnsp(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C22Tnsp::Regular(double const&) const
  {
    return 0 * _nf;
  }
  double C22Tnsp::Singular(double const&) const
  {
    return 0;
  }
  double C22Tnsp::Local(double const&) const
  {
    return 0;
  }

  //_________________________________________________________________________________
  C22Tnsm::C22Tnsm(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C22Tnsm::Regular(double const&) const
  {
    return 0 * _nf;
  }
  double C22Tnsm::Singular(double const&) const
  {
    return 0;
  }
  double C22Tnsm::Local(double const&) const
  {
    return 0;
  }

  //_________________________________________________________________________________
  C22Tps::C22Tps():
    Expression()
  {
  }
  double C22Tps::Regular(double const&) const
  {
    return 0;
  }

  //_________________________________________________________________________________
  C22Tg::C22Tg():
    Expression()
  {
  }
  double C22Tg::Regular(double const&) const
  {
    return 0;
  }
  double C22Tg::Local(double const&) const
  {
    return 0;
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
  double CL2Tnsp::Regular(double const&) const
  {
    return 0 * _nf;
  }
  double CL2Tnsp::Local(double const&) const
  {
    return 0;
  }

  //_________________________________________________________________________________
  CL2Tnsm::CL2Tnsm(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double CL2Tnsm::Regular(double const&) const
  {
    return 0 * _nf;
  }
  double CL2Tnsm::Local(double const&) const
  {
    return 0;
  }

  //_________________________________________________________________________________
  CL2Tps::CL2Tps():
    Expression()
  {
  }
  double CL2Tps::Regular(double const&) const
  {
    return 0;
  }

  //_________________________________________________________________________________
  CL2Tg::CL2Tg():
    Expression()
  {
  }
  double CL2Tg::Regular(double const&) const
  {
    return 0;
  }

  //_________________________________________________________________________________
  C31Tns::C31Tns():
    Expression()
  {
  }
  double C31Tns::Regular(double const& x) const
  {
    return 2 * CF * ( - ( 1 + x ) * log( 1 - x ) - 2 * ( 1 + pow(x,2) ) * log(x) / ( 1 - x ) + 1. / 2. - x / 2 );
  }
  double C31Tns::Singular(double const& x) const
  {
    return 2 * CF * ( 2 * log( 1 - x ) - 3 / 2. ) / ( 1 - x );
  }
  double C31Tns::Local(double const& x) const
  {
    return 2 * CF * ( pow(log(1-x),2) - 3 * log( 1 - x ) / 2 + ( 4 * zeta2 - 9 / 2. ) );
  }

  //_________________________________________________________________________________
  C32Tnsp::C32Tnsp(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C32Tnsp::Regular(double const&) const
  {
    return 0 * _nf;
  }
  double C32Tnsp::Singular(double const&) const
  {
    return 0;
  }
  double C32Tnsp::Local(double const&) const
  {
    return 0;
  }

  //_________________________________________________________________________________
  C32Tnsm::C32Tnsm(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C32Tnsm::Regular(double const&) const
  {
    return 0 * _nf;
  }
  double C32Tnsm::Singular(double const&) const
  {
    return 0;
  }
  double C32Tnsm::Local(double const&) const
  {
    return 0;
  }
}
