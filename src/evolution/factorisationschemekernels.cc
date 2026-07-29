//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/factorisationschemekernels.h"
#include "apfel/constants.h"

#include <cmath>

namespace apfel
{
  //_________________________________________________________________________________
  K1qgPOS::K1qgPOS():
    Expression()
  {
  }
  double K1qgPOS::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double ln1mx = log(1 - x);
    const double lnx   = log(x);
    return 2 * ( 2 * TR * ( 2 * x2 - 2 * x + 1 ) * ln1mx
                 + TR * ( -2 * x2 + 2 * x - 1 ) * lnx
                 + TR * ( -2 * x2 + 2 * x - 1 ) );
  }

  //_________________________________________________________________________________
  K1gqPOS::K1gqPOS():
    Expression()
  {
  }
  double K1gqPOS::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double ln1mx = log(1 - x);
    const double lnx   = log(x);
    return 2 * ( 2 * CF * ( x2 - 2 * x + 2 ) / x * ln1mx
                 + CF * ( - x2 + 2 * x - 2 ) / x * lnx
                 + CF * ( - x2 + 2 * x - 2 ) / x );
  }

  //_________________________________________________________________________________
  K1qqMPOS::K1qqMPOS():
    Expression()
  {
  }
  double K1qqMPOS::Regular(double const& x) const
  {
    const double xm1 = x - 1;
    return 2 * ( ( 350. / 3. ) * CF * x * x * xm1 * xm1 );
  }

  //_________________________________________________________________________________
  K1qgMPOS::K1qgMPOS():
    Expression()
  {
  }
  double K1qgMPOS::Regular(double const& x) const
  {
    return K1qgPOS{}.Regular(x);
  }

  //_________________________________________________________________________________
  K1gqMPOS::K1gqMPOS():
    Expression()
  {
  }
  double K1gqMPOS::Regular(double const& x) const
  {
    return K1gqPOS{}.Regular(x);
  }

  //_________________________________________________________________________________
  K1ggMPOS::K1ggMPOS(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double K1ggMPOS::Regular(double const& x) const
  {
    const double xm1 = x - 1;
    return 2 * ( ( 475. / 3. ) * TR * _nf * x * x * xm1 * xm1 );
  }

  //_________________________________________________________________________________
  K1qqMPOSDelta::K1qqMPOSDelta():
    Expression()
  {
  }
  double K1qqMPOSDelta::Local(double const&) const
  {
    return 2 * ( ( 35. / 18. ) * CF );
  }

  //_________________________________________________________________________________
  K1qgMPOSDelta::K1qgMPOSDelta():
    Expression()
  {
  }
  double K1qgMPOSDelta::Regular(double const& x) const
  {
    return K1qgPOS{}.Regular(x);
  }

  //_________________________________________________________________________________
  K1gqMPOSDelta::K1gqMPOSDelta():
    Expression()
  {
  }
  double K1gqMPOSDelta::Regular(double const& x) const
  {
    return K1gqPOS{}.Regular(x);
  }

  //_________________________________________________________________________________
  K1ggMPOSDelta::K1ggMPOSDelta(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double K1ggMPOSDelta::Local(double const&) const
  {
    return 2 * ( ( 95. / 36. ) * TR * _nf );
  }

  //_________________________________________________________________________________
  K1qqKrk::K1qqKrk():
    Expression()
  {
  }
  double K1qqKrk::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double ln1mx = log(1 - x);
    const double lnx   = log(x);
    return 2 * ( 2 * CF * ( - x - 1 ) * ln1mx
                 + CF * ( x2 + 1 ) / ( x - 1 ) * lnx
                 + CF * ( 1 - x ) );
  }
  double K1qqKrk::Singular(double const& x) const
  {
    return 2 * 4 * CF * log(1 - x) / ( 1 - x );
  }
  double K1qqKrk::Local(double const& x) const
  {
    const double ln1mx = log(1 - x);
    return 2 * ( - CF * ( Pi2 / 3 + 17. / 4. ) + 2 * CF * ln1mx * ln1mx );
  }

  //_________________________________________________________________________________
  K1qgKrk::K1qgKrk():
    Expression()
  {
  }
  double K1qgKrk::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double ln1mx = log(1 - x);
    const double lnx   = log(x);
    return 2 * ( 2 * TR * ( 2 * x2 - 2 * x + 1 ) * ln1mx
                 + TR * ( -2 * x2 + 2 * x - 1 ) * lnx
                 + 2 * TR * x * ( 1 - x ) );
  }

  //_________________________________________________________________________________
  K1gqKrk::K1gqKrk():
    Expression()
  {
  }
  double K1gqKrk::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double ln1mx = log(1 - x);
    const double lnx   = log(x);
    return 2 * ( 2 * CF * ( x2 - 2 * x + 2 ) / x * ln1mx
                 + CF * ( - x2 + 2 * x - 2 ) / x * lnx
                 + CF * x );
  }

  //_________________________________________________________________________________
  K1ggKrk::K1ggKrk(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double K1ggKrk::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double x3    = x * x2;
    const double xm1   = x - 1;
    const double ln1mx = log(1 - x);
    const double lnx   = log(x);
    const double b = 4 * CA * ( - x3 + x2 - 2 * x + 1 ) / x;
    const double c = 2 * CA * ( x * xm1 * ( x2 - x + 2 ) + 1 ) / ( x * xm1 );
    return 2 * ( b * ln1mx + c * lnx );
  }
  double K1ggKrk::Singular(double const& x) const
  {
    return 2 * 4 * CA * log(1 - x) / ( 1 - x );
  }
  double K1ggKrk::Local(double const& x) const
  {
    const double ln1mx = log(1 - x);
    return 2 * ( - CA * ( Pi2 / 3 + 341. / 72. ) + ( 59. / 36. ) * TR * _nf
                 + 2 * CA * ln1mx * ln1mx );
  }

  //_________________________________________________________________________________
  K1qqAversa::K1qqAversa():
    Expression()
  {
  }
  double K1qqAversa::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double ln1mx = log(1 - x);
    const double lnx   = log(x);
    return 2 * ( CF * ( - x - 1 ) * ln1mx
                 + CF * ( x2 + 1 ) / ( x - 1 ) * lnx
                 + CF * ( 2 * x + 3 ) );
  }
  double K1qqAversa::Singular(double const& x) const
  {
    return 2 * ( ( -3. / 2. * CF + 2 * CF * log(1 - x) ) / ( 1 - x ) );
  }
  double K1qqAversa::Local(double const& x) const
  {
    const double ln1mx = log(1 - x);
    return 2 * ( - CF * ( Pi2 / 3 + 9. / 2. ) - 3. / 2. * CF * ln1mx + CF * ln1mx * ln1mx );
  }

  //_________________________________________________________________________________
  K1qgAversa::K1qgAversa():
    Expression()
  {
  }
  double K1qgAversa::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double ln1mx = log(1 - x);
    const double lnx   = log(x);
    return 2 * ( TR * ( 2 * x2 - 2 * x + 1 ) * ln1mx
                 + TR * ( -2 * x2 + 2 * x - 1 ) * lnx );
  }

  //_________________________________________________________________________________
  K1gqAversa::K1gqAversa():
    Expression()
  {
  }
  double K1gqAversa::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double ln1mx = log(1 - x);
    const double lnx   = log(x);
    return 2 * ( CF * ( x2 - 2 * x + 2 ) / x * ln1mx
                 + CF * ( - x2 + 2 * x - 2 ) / x * lnx
                 - 4. / 3. * CF );
  }

  //_________________________________________________________________________________
  K1ggAversa::K1ggAversa(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double K1ggAversa::Regular(double const& x) const
  {
    const double ln1mx = log(1 - x);
    const double lnx   = log(x);
    return 2 * ( - 2 * CA * ln1mx + 2 * CA * x / ( x - 1 ) * lnx );
  }
  double K1ggAversa::Singular(double const& x) const
  {
    return 2 * ( 2 * CA * log(1 - x) / ( 1 - x ) );
  }
  double K1ggAversa::Local(double const& x) const
  {
    const double ln1mx = log(1 - x);
    return 2 * ( - CA * ( 1 + Pi2 / 3 ) + 5. / 6. * TR * _nf + CA * ln1mx * ln1mx );
  }

  //_________________________________________________________________________________
  K1qqDIS::K1qqDIS():
    Expression()
  {
  }
  double K1qqDIS::Regular(double const& x) const
  {
    return K1qqAversa{}.Regular(x);
  }
  double K1qqDIS::Singular(double const& x) const
  {
    return K1qqAversa{}.Singular(x);
  }
  double K1qqDIS::Local(double const& x) const
  {
    return K1qqAversa{}.Local(x);
  }

  //_________________________________________________________________________________
  K1qgDIS::K1qgDIS():
    Expression()
  {
  }
  double K1qgDIS::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double ln1mx = log(1 - x);
    const double lnx   = log(x);
    return 2 * ( TR * ( 2 * x2 - 2 * x + 1 ) * ln1mx
                 + TR * ( -2 * x2 + 2 * x - 1 ) * lnx
                 + TR * ( -8 * x2 + 8 * x - 1 ) );
  }

  //_________________________________________________________________________________
  K1gqDIS::K1gqDIS():
    Expression()
  {
  }
  double K1gqDIS::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double ln1mx = log(1 - x);
    const double lnx   = log(x);
    return 2 * ( CF * ( x + 1 ) * ln1mx
                 + CF * ( - x2 - 1 ) / ( x - 1 ) * lnx
                 + CF * ( -2 * x - 3 ) );
  }
  double K1gqDIS::Singular(double const& x) const
  {
    return 2 * ( ( 3. / 2. * CF - 2 * CF * log(1 - x) ) / ( 1 - x ) );
  }
  double K1gqDIS::Local(double const& x) const
  {
    const double ln1mx = log(1 - x);
    return 2 * ( CF * ( 9. / 2. + Pi2 / 3 ) + 3. / 2. * CF * ln1mx - CF * ln1mx * ln1mx );
  }

  //_________________________________________________________________________________
  K1ggDIS::K1ggDIS(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double K1ggDIS::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double ln1mx = log(1 - x);
    const double lnx   = log(x);
    return 2 * ( 2 * TR * _nf * ( -2 * x2 + 2 * x - 1 ) * ln1mx
                 + 2 * TR * _nf * (  2 * x2 - 2 * x + 1 ) * lnx
                 + 2 * TR * _nf * (  8 * x2 - 8 * x + 1 ) );
  }

  //_________________________________________________________________________________
  K1qqPHYS::K1qqPHYS():
    Expression()
  {
  }
  double K1qqPHYS::Regular(double const& x) const
  {
    const double ln1mx = log(1 - x);
    return 2 * ( CF * ( - x - 1 ) * ln1mx + CF * ( 1 - x ) );
  }
  double K1qqPHYS::Singular(double const& x) const
  {
    return 2 * ( 2 * CF * log(1 - x) / ( 1 - x ) );
  }
  double K1qqPHYS::Local(double const& x) const
  {
    const double ln1mx = log(1 - x);
    return 2 * ( - 11. / 4. * CF + CF * ln1mx * ln1mx );
  }

  //_________________________________________________________________________________
  K1qgPHYS::K1qgPHYS():
    Expression()
  {
  }
  double K1qgPHYS::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double ln1mx = log(1 - x);
    return 2 * ( TR * ( 2 * x2 - 2 * x + 1 ) * ln1mx + 2 * TR * x * ( 1 - x ) );
  }

  //_________________________________________________________________________________
  K1gqPHYS::K1gqPHYS():
    Expression()
  {
  }
  double K1gqPHYS::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double ln1mx = log(1 - x);
    return 2 * ( CF * ( x2 - 2 * x + 2 ) / x * ln1mx + CF * x );
  }

  //_________________________________________________________________________________
  K1ggPHYS::K1ggPHYS(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double K1ggPHYS::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double x3    = x * x2;
    const double ln1mx = log(1 - x);
    return 2 * ( 2 * CA * ( - x3 + x2 - 2 * x + 1 ) / x * ln1mx );
  }
  double K1ggPHYS::Singular(double const& x) const
  {
    return 2 * ( 2 * CA * log(1 - x) / ( 1 - x ) );
  }
  double K1ggPHYS::Local(double const& x) const
  {
    const double ln1mx = log(1 - x);
    return 2 * ( - CA * ( 203. / 72. ) + ( 29. / 36. ) * TR * _nf + CA * ln1mx * ln1mx );
  }
}
