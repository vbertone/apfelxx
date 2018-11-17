//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/lhtoypdfs.h"

#include <cmath>

namespace apfel
{
  //_________________________________________________________________________________
  double xupv(double const& x)
  {
    return 5.107200 * pow(x, 0.8) * pow(1 - x, 3);
  }

  //_________________________________________________________________________________
  double xdnv(double const& x)
  {
    return 3.064320 * pow(x, 0.8) * pow(1 - x, 4);
  }

  //_________________________________________________________________________________
  double xglu(double const& x)
  {
    return 1.7 * pow(x, - 0.1) * pow(1 - x, 5);
  }

  //_________________________________________________________________________________
  double xdbar(double const& x)
  {
    return 0.1939875 * pow(x, - 0.1) * pow(1 - x, 6);
  }

  //_________________________________________________________________________________
  double xubar(double const& x)
  {
    return xdbar(x) * ( 1 - x );
  }

  //_________________________________________________________________________________
  double xsbar(double const& x)
  {
    return 0.2 * ( xdbar(x) + xubar(x) );
  }

  //_________________________________________________________________________________
  std::map<int,double> LHToyPDFs(double const& x, double const&)
  {
    // Call all functions only once.
    const double upv  = xupv(x);
    const double dnv  = xdnv(x);
    const double glu  = xglu(x);
    const double dbar = xdbar(x);
    const double ubar = xubar(x);
    const double sbar = xsbar(x);

    // Construct QCD evolution basis conbinations.
    double const Gluon   = glu;
    double const Singlet = dnv + 2 * dbar + upv + 2 * ubar + 2 * sbar;
    double const T3      = upv + 2 * ubar - dnv - 2 * dbar;
    double const T8      = upv + 2 * ubar + dnv + 2 * dbar - 4 * sbar;
    double const Valence = upv + dnv;
    double const V3      = upv - dnv;

    // Fill in map in the QCD evolution basis.
    std::map<int,double> QCDEvMap;
    QCDEvMap[0]  = Gluon;
    QCDEvMap[1]  = Singlet;
    QCDEvMap[2]  = Valence;
    QCDEvMap[3]  = T3;
    QCDEvMap[4]  = V3;
    QCDEvMap[5]  = T8;
    QCDEvMap[6]  = Valence;
    QCDEvMap[7]  = Singlet;
    QCDEvMap[8]  = Valence;
    QCDEvMap[9]  = Singlet;
    QCDEvMap[10] = Valence;
    QCDEvMap[11] = Singlet;
    QCDEvMap[12] = Valence;

    return QCDEvMap;
  }
}
