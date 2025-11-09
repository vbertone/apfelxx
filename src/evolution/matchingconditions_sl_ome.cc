//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/matchingconditions_sl_ome.h"
#include "apfel/integrator.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________________
  AS1Hg_L_ome::AS1Hg_L_ome():
    Expression(),
    _ome_r(ome::AQg_reg[1][1])
  {
  }
  double AS1Hg_L_ome::Regular(double const& x) const
  {
    return - _ome_r(0, x);
  }

  //_________________________________________________________________________________
  AS1ggH_L_ome::AS1ggH_L_ome():
    Expression(),
    _ome_l(ome::AggQ_delta[1][1])
  {
  }
  double AS1ggH_L_ome::Local(double const&) const
  {
    return - _ome_l(0);
  }

  //_________________________________________________________________________________
  APS2Hq_0_ome::APS2Hq_0_ome():
    Expression(),
    _ome_r(ome::AQqPS_reg[2][0])
  {
  }
  double APS2Hq_0_ome::Regular(double const& x) const
  {
    return _ome_r(0, x);
  }

  //_________________________________________________________________________________
  APS2Hq_L_ome::APS2Hq_L_ome():
    Expression(),
    _ome_r(ome::AQqPS_reg[2][1])
  {
  }
  double APS2Hq_L_ome::Regular(double const& x) const
  {
    return - _ome_r(0, x);
  }

  //_________________________________________________________________________________
  APS2Hq_L2_ome::APS2Hq_L2_ome():
    Expression(),
    _ome_r(ome::AQqPS_reg[2][2])
  {
  }
  double APS2Hq_L2_ome::Regular(double const& x) const
  {
    return _ome_r(0, x);
  }

  //_________________________________________________________________________________
  AS2Hg_0_ome::AS2Hg_0_ome():
    Expression(),
    _ome_r(ome::AQg_reg[2][0])
  {
  }
  double AS2Hg_0_ome::Regular(double const& x) const
  {
    return _ome_r(0, x);
  }

  //_________________________________________________________________________________
  AS2Hg_L_ome::AS2Hg_L_ome():
    Expression(),
    _ome_r(ome::AQg_reg[2][1])
  {
  }
  double AS2Hg_L_ome::Regular(double const& x) const
  {
    return - _ome_r(0, x);
  }

  //_________________________________________________________________________________
  AS2Hg_L2_ome::AS2Hg_L2_ome():
    Expression(),
    _ome_r(ome::AQg_reg[2][2])
  {
  }
  double AS2Hg_L2_ome::Regular(double const& x) const
  {
    return _ome_r(0, x);
  }

  //_________________________________________________________________________________
  ANS2qqH_0_ome::ANS2qqH_0_ome():
    Expression(),
    _ome_r(ome::AqqQNSEven_reg[2][0]),
    _ome_s(ome::AqqQNSEven_plus[2][0]),
    _ome_l(ome::AqqQNSEven_delta[2][0])
  {
  }
  double ANS2qqH_0_ome::Regular(double const& x) const
  {
    return _ome_r(0, x);
  }
  double ANS2qqH_0_ome::Singular(double const& x) const
  {
    return _ome_s(0, x);
  }
  double ANS2qqH_0_ome::Local(double const& x) const
  {
    return - Integrator{[=] (double const& y) -> double{ return _ome_s(0, y); }}.integrate(0, x, eps7) + _ome_l(0);
  }

  //_________________________________________________________________________________
  ANS2qqH_L_ome::ANS2qqH_L_ome():
    Expression(),
    _ome_r(ome::AqqQNSEven_reg[2][1]),
    _ome_s(ome::AqqQNSEven_plus[2][1]),
    _ome_l(ome::AqqQNSEven_delta[2][1])
  {
  }
  double ANS2qqH_L_ome::Regular(double const& x) const
  {
    return - _ome_r(0, x);
  }
  double ANS2qqH_L_ome::Singular(double const& x) const
  {
    return - _ome_s(0, x);
  }
  double ANS2qqH_L_ome::Local(double const& x) const
  {
    return Integrator{[=] (double const& y) -> double{ return _ome_s(0, y); }}.integrate(0, x, eps7) - _ome_l(0);
  }

  //_________________________________________________________________________________
  ANS2qqH_L2_ome::ANS2qqH_L2_ome():
    Expression(),
    _ome_r(ome::AqqQNSEven_reg[2][2]),
    _ome_s(ome::AqqQNSEven_plus[2][2]),
    _ome_l(ome::AqqQNSEven_delta[2][2])
  {
  }
  double ANS2qqH_L2_ome::Regular(double const& x) const
  {
    return _ome_r(0, x);
  }
  double ANS2qqH_L2_ome::Singular(double const& x) const
  {
    return _ome_s(0, x);
  }
  double ANS2qqH_L2_ome::Local(double const& x) const
  {
    return - Integrator{[=] (double const& y) -> double{ return _ome_s(0, y); }}.integrate(0, x, eps7) + _ome_l(0);
  }

  //_________________________________________________________________________________
  AS2gqH_0_ome::AS2gqH_0_ome():
    Expression(),
    _ome_r(ome::AgqQ_reg[2][0])
  {
  }
  double AS2gqH_0_ome::Regular(double const& x) const
  {
    return _ome_r(0, x);
  }

  //_________________________________________________________________________________
  AS2gqH_L_ome::AS2gqH_L_ome():
    Expression(),
    _ome_r(ome::AgqQ_reg[2][1])
  {
  }
  double AS2gqH_L_ome::Regular(double const& x) const
  {
    return - _ome_r(0, x);
  }

  //_________________________________________________________________________________
  AS2gqH_L2_ome::AS2gqH_L2_ome():
    Expression(),
    _ome_r(ome::AgqQ_reg[2][2])
  {
  }
  double AS2gqH_L2_ome::Regular(double const& x) const
  {
    return _ome_r(0, x);
  }

  //_________________________________________________________________________________
  AS2ggH_0_ome::AS2ggH_0_ome():
    Expression(),
    _ome_r(ome::AggQ_reg[2][0]),
    _ome_s(ome::AggQ_plus[2][0]),
    _ome_l(ome::AggQ_delta[2][0])
  {
  }
  double AS2ggH_0_ome::Regular(double const& x) const
  {
    return _ome_r(0, x);
  }
  double AS2ggH_0_ome::Singular(double const& x) const
  {
    return _ome_s(0, x);
  }
  double AS2ggH_0_ome::Local(double const& x) const
  {
    return - Integrator{[=] (double const& y) -> double{ return _ome_s(0, y); }}.integrate(0, x, eps7) + _ome_l(0);
  }

  //_________________________________________________________________________________
  AS2ggH_L_ome::AS2ggH_L_ome():
    Expression(),
    _ome_r(ome::AggQ_reg[2][1]),
    _ome_s(ome::AggQ_plus[2][1]),
    _ome_l(ome::AggQ_delta[2][1])
  {
  }
  double AS2ggH_L_ome::Regular(double const& x) const
  {
    return - _ome_r(0, x);
  }
  double AS2ggH_L_ome::Singular(double const& x) const
  {
    return - _ome_s(0, x);
  }
  double AS2ggH_L_ome::Local(double const& x) const
  {
    return Integrator{[=] (double const& y) -> double{ return _ome_s(0, y); }}.integrate(0, x, eps7) - _ome_l(0);
  }

  //_________________________________________________________________________________
  AS2ggH_L2_ome::AS2ggH_L2_ome():
    Expression(),
    _ome_r(ome::AggQ_reg[2][2]),
    _ome_s(ome::AggQ_plus[2][2]),
    _ome_l(ome::AggQ_delta[2][2])
  {
  }
  double AS2ggH_L2_ome::Regular(double const& x) const
  {
    return _ome_r(0, x);
  }
  double AS2ggH_L2_ome::Singular(double const& x) const
  {
    return _ome_s(0, x);
  }
  double AS2ggH_L2_ome::Local(double const& x) const
  {
    return - Integrator{[=] (double const& y) -> double{ return _ome_s(0, y); }}.integrate(0, x, eps7) + _ome_l(0);
  }

  //_________________________________________________________________________________
  APS3Hq_0_ome::APS3Hq_0_ome(int const& nf):
    Expression(),
    _nf(nf),
    _ome_r(ome::AQqPS_reg[3][0])
  {
  }
  double APS3Hq_0_ome::Regular(double const& x) const
  {
    return _ome_r(_nf, x);
  }

  //_________________________________________________________________________________
  AS3Hg_0_ome::AS3Hg_0_ome(int const& nf):
    Expression(),
    _nf(nf),
    _ome_r(ome::AQg_reg[3][0])
  {
  }
  double AS3Hg_0_ome::Regular(double const& x) const
  {
    return _ome_r(_nf, x);
  }

  //_________________________________________________________________________________
  ANS3qqH_0_ome::ANS3qqH_0_ome(int const& nf):
    Expression(),
    _nf(nf),
    _ome_r(ome::AqqQNSEven_reg[3][0]),
    _ome_s(ome::AqqQNSEven_plus[3][0]),
    _ome_l(ome::AqqQNSEven_delta[3][0])
  {
  }
  double ANS3qqH_0_ome::Regular(double const& x) const
  {
    return _ome_r(_nf, x);
  }
  double ANS3qqH_0_ome::Singular(double const& x) const
  {
    return _ome_s(_nf, x);
  }
  double ANS3qqH_0_ome::Local(double const& x) const
  {
    return - Integrator{[=] (double const& y) -> double{ return _ome_s(_nf, y); }}.integrate(0, x, eps7) + _ome_l(_nf);
  }

  //_________________________________________________________________________________
  ANS3qqHm_0_ome::ANS3qqHm_0_ome(int const& nf):
    Expression(),
    _nf(nf),
    _ome_r(ome::AqqQNSOdd_reg[3][0]),
    _ome_s(ome::AqqQNSOdd_plus[3][0]),
    _ome_l(ome::AqqQNSOdd_delta[3][0])
  {
  }
  double ANS3qqHm_0_ome::Regular(double const& x) const
  {
    return _ome_r(_nf, x);
  }
  double ANS3qqHm_0_ome::Singular(double const& x) const
  {
    return _ome_s(_nf, x);
  }
  double ANS3qqHm_0_ome::Local(double const& x) const
  {
    return - Integrator{[=] (double const& y) -> double{ return _ome_s(_nf, y); }}.integrate(0, x, eps7) + _ome_l(_nf);
  }

  //_________________________________________________________________________________
  AS3gqH_0_ome::AS3gqH_0_ome(int const& nf):
    Expression(),
    _nf(nf),
    _ome_r(ome::AgqQ_reg[3][0])
  {
  }
  double AS3gqH_0_ome::Regular(double const& x) const
  {
    return _ome_r(_nf, x);
  }

  //_________________________________________________________________________________
  AS3ggH_0_ome::AS3ggH_0_ome(int const& nf):
    Expression(),
    _nf(nf),
    _ome_r(ome::AggQ_reg[3][0]),
    _ome_s(ome::AggQ_plus[3][0]),
    _ome_l(ome::AggQ_delta[3][0])
  {
  }
  double AS3ggH_0_ome::Regular(double const& x) const
  {
    return _ome_r(_nf, x);
  }
  double AS3ggH_0_ome::Singular(double const& x) const
  {
    return _ome_s(_nf, x);
  }
  double AS3ggH_0_ome::Local(double const& x) const
  {
    return - Integrator{[=] (double const& y) -> double{ return _ome_s(_nf, y); }}.integrate(0, x, eps7) + _ome_l(_nf);
  }

  //_________________________________________________________________________________
  AS3qgQ_0_ome::AS3qgQ_0_ome(int const& nf):
    Expression(),
    _nf(nf),
    _ome_r(ome::AqgQ_reg[3][0])
  {
  }
  double AS3qgQ_0_ome::Regular(double const& x) const
  {
    return _ome_r(_nf, x);
  }

  //_________________________________________________________________________________
  APS3qqQ_0_ome::APS3qqQ_0_ome(int const& nf):
    Expression(),
    _nf(nf),
    _ome_r(ome::AqqQPS_reg[3][0])
  {
  }
  double APS3qqQ_0_ome::Regular(double const& x) const
  {
    return _ome_r(_nf, x);
  }

  //_________________________________________________________________________________
  AS1polHg_L_ome::AS1polHg_L_ome():
    Expression(),
    _ome_r(ome::polAQg_reg[1][1])
  {
  }
  double AS1polHg_L_ome::Regular(double const& x) const
  {
    return - _ome_r(0, x);
  }

  //_________________________________________________________________________________
  AS1polggH_L_ome::AS1polggH_L_ome():
    Expression(),
    _ome_l(ome::polAggQ_delta[1][1])
  {
  }
  double AS1polggH_L_ome::Local(double const&) const
  {
    return - _ome_l(0);
  }

  //_________________________________________________________________________________
  APS2polHq_0_ome::APS2polHq_0_ome():
    Expression(),
    _ome_r(ome::polAQqPS_reg[2][0])
  {
  }
  double APS2polHq_0_ome::Regular(double const& x) const
  {
    return _ome_r(0, x);
  }

  //_________________________________________________________________________________
  APS2polHq_L_ome::APS2polHq_L_ome():
    Expression(),
    _ome_r(ome::polAQqPS_reg[2][1])
  {
  }
  double APS2polHq_L_ome::Regular(double const& x) const
  {
    return - _ome_r(0, x);
  }

  //_________________________________________________________________________________
  APS2polHq_L2_ome::APS2polHq_L2_ome():
    Expression(),
    _ome_r(ome::polAQqPS_reg[2][2])
  {
  }
  double APS2polHq_L2_ome::Regular(double const& x) const
  {
    return _ome_r(0, x);
  }

  //_________________________________________________________________________________
  AS2polHg_0_ome::AS2polHg_0_ome():
    Expression(),
    _ome_r(ome::polAQg_reg[2][0])
  {
  }
  double AS2polHg_0_ome::Regular(double const& x) const
  {
    return _ome_r(0, x);
  }

  //_________________________________________________________________________________
  AS2polHg_L_ome::AS2polHg_L_ome():
    Expression(),
    _ome_r(ome::polAQg_reg[2][1])
  {
  }
  double AS2polHg_L_ome::Regular(double const& x) const
  {
    return - _ome_r(0, x);
  }

  //_________________________________________________________________________________
  AS2polHg_L2_ome::AS2polHg_L2_ome():
    Expression(),
    _ome_r(ome::polAQg_reg[2][2])
  {
  }
  double AS2polHg_L2_ome::Regular(double const& x) const
  {
    return _ome_r(0, x);
  }

  //_________________________________________________________________________________
  ANS2polqqH_0_ome::ANS2polqqH_0_ome():
    Expression(),
    _ome_r(ome::polAqqQNSEven_reg[2][0]),
    _ome_s(ome::polAqqQNSEven_plus[2][0]),
    _ome_l(ome::polAqqQNSEven_delta[2][0])
  {
  }
  double ANS2polqqH_0_ome::Regular(double const& x) const
  {
    return _ome_r(0, x);
  }
  double ANS2polqqH_0_ome::Singular(double const& x) const
  {
    return _ome_s(0, x);
  }
  double ANS2polqqH_0_ome::Local(double const& x) const
  {
    return - Integrator{[=] (double const& y) -> double{ return _ome_s(0, y); }}.integrate(0, x, eps7) + _ome_l(0);
  }

  //_________________________________________________________________________________
  ANS2polqqH_L_ome::ANS2polqqH_L_ome():
    Expression(),
    _ome_r(ome::polAqqQNSEven_reg[2][1]),
    _ome_s(ome::polAqqQNSEven_plus[2][1]),
    _ome_l(ome::polAqqQNSEven_delta[2][1])
  {
  }
  double ANS2polqqH_L_ome::Regular(double const& x) const
  {
    return - _ome_r(0, x);
  }
  double ANS2polqqH_L_ome::Singular(double const& x) const
  {
    return - _ome_s(0, x);
  }
  double ANS2polqqH_L_ome::Local(double const& x) const
  {
    return Integrator{[=] (double const& y) -> double{ return _ome_s(0, y); }}.integrate(0, x, eps7) - _ome_l(0);
  }

  //_________________________________________________________________________________
  ANS2polqqH_L2_ome::ANS2polqqH_L2_ome():
    Expression(),
    _ome_r(ome::polAqqQNSEven_reg[2][2]),
    _ome_s(ome::polAqqQNSEven_plus[2][2]),
    _ome_l(ome::polAqqQNSEven_delta[2][2])
  {
  }
  double ANS2polqqH_L2_ome::Regular(double const& x) const
  {
    return _ome_r(0, x);
  }
  double ANS2polqqH_L2_ome::Singular(double const& x) const
  {
    return _ome_s(0, x);
  }
  double ANS2polqqH_L2_ome::Local(double const& x) const
  {
    return - Integrator{[=] (double const& y) -> double{ return _ome_s(0, y); }}.integrate(0, x, eps7) + _ome_l(0);
  }

  //_________________________________________________________________________________
  AS2polgqH_0_ome::AS2polgqH_0_ome():
    Expression(),
    _ome_r(ome::polAgqQ_reg[2][0])
  {
  }
  double AS2polgqH_0_ome::Regular(double const& x) const
  {
    return _ome_r(0, x);
  }

  //_________________________________________________________________________________
  AS2polgqH_L_ome::AS2polgqH_L_ome():
    Expression(),
    _ome_r(ome::polAgqQ_reg[2][1])
  {
  }
  double AS2polgqH_L_ome::Regular(double const& x) const
  {
    return - _ome_r(0, x);
  }

  //_________________________________________________________________________________
  AS2polgqH_L2_ome::AS2polgqH_L2_ome():
    Expression(),
    _ome_r(ome::polAgqQ_reg[2][2])
  {
  }
  double AS2polgqH_L2_ome::Regular(double const& x) const
  {
    return _ome_r(0, x);
  }

  //_________________________________________________________________________________
  AS2polggH_0_ome::AS2polggH_0_ome():
    Expression(),
    _ome_r(ome::polAggQ_reg[2][0]),
    _ome_s(ome::polAggQ_plus[2][0]),
    _ome_l(ome::polAggQ_delta[2][0])
  {
  }
  double AS2polggH_0_ome::Regular(double const& x) const
  {
    return _ome_r(0, x);
  }
  double AS2polggH_0_ome::Singular(double const& x) const
  {
    return _ome_s(0, x);;
  }
  double AS2polggH_0_ome::Local(double const& x) const
  {
    return - Integrator{[=] (double const& y) -> double{ return _ome_s(0, y); }}.integrate(0, x, eps7) + _ome_l(0);
  }

  //_________________________________________________________________________________
  AS2polggH_L_ome::AS2polggH_L_ome():
    Expression(),
    _ome_r(ome::polAggQ_reg[2][1]),
    _ome_s(ome::polAggQ_plus[2][1]),
    _ome_l(ome::polAggQ_delta[2][1])
  {
  }
  double AS2polggH_L_ome::Regular(double const& x) const
  {
    return - _ome_r(0, x);
  }
  double AS2polggH_L_ome::Singular(double const& x) const
  {
    return - _ome_s(0, x);
  }
  double AS2polggH_L_ome::Local(double const& x) const
  {
    return Integrator{[=] (double const& y) -> double{ return _ome_s(0, y); }}.integrate(0, x, eps7) - _ome_l(0);
  }

  //_________________________________________________________________________________
  AS2polggH_L2_ome::AS2polggH_L2_ome():
    Expression(),
    _ome_r(ome::polAggQ_reg[2][2]),
    _ome_s(ome::polAggQ_plus[2][2]),
    _ome_l(ome::polAggQ_delta[2][2])
  {
  }
  double AS2polggH_L2_ome::Regular(double const& x) const
  {
    return _ome_r(0, x);
  }
  double AS2polggH_L2_ome::Singular(double const& x) const
  {
    return _ome_s(0, x);
  }
  double AS2polggH_L2_ome::Local(double const& x) const
  {
    return - Integrator{[=] (double const& y) -> double{ return _ome_s(0, y); }}.integrate(0, x, eps7) + _ome_l(0);
  }
}
