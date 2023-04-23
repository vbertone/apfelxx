//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/lhtoypdfs.h"
#include "apfel/rotations.h"
#include "apfel/messages.h"

#include <cmath>

namespace apfel
{
  //_________________________________________________________________________________
  std::map<int, double> LHToyPDFs(double const& x, double const&)
  {
    return Provider(x, 0);
  }

  //_________________________________________________________________________________
  std::map<int, double> LHToyPDFsPhys(double const& x, double const&)
  {
    return QCDEvToPhys(Provider(x, 0));
  }

  //_________________________________________________________________________________
  std::map<int, double> LHToyPDFsPol(double const& x, double const&)
  {
    return Provider(x, 1);
  }

  //_________________________________________________________________________________
  std::map<int, double> LHToyFFs(double const& x, double const&)
  {
    return Provider(x, 2);
  }

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
  double xupvpol(double const& x)
  {
    return 1.3 * pow(x, 0.7) * pow(1 - x, 3) * ( 1 + 3 * x );
  }

  //_________________________________________________________________________________
  double xdnvpol(double const& x)
  {
    return - 0.5 * pow(x, 0.7) * pow(1 - x, 4) * ( 1 + 4 * x );
  }

  //_________________________________________________________________________________
  double xglupol(double const& x)
  {
    return 1.5 * pow(x, 0.5) * pow(1 - x, 5);
  }

  //_________________________________________________________________________________
  double xdbarpol(double const& x)
  {
    return - 0.05 * pow(x, 0.3) * pow(1 - x, 7);
  }

  //_________________________________________________________________________________
  double xubarpol(double const& x)
  {
    return xdbarpol(x);
  }

  //_________________________________________________________________________________
  double xsbarpol(double const& x)
  {
    return 0.5 * xdbarpol(x);
  }

  //_________________________________________________________________________________
  double xsff(double const& x)
  {
    const double as = 0.718;
    const double bs = 6.266;
    const double Ns = 0.094 / ( std::tgamma(as + 2) * std::tgamma(bs + 1) / std::tgamma(as + bs + 3) );
    return Ns * pow(x, as + 1) * pow(1 - x, bs);
  }

  //_________________________________________________________________________________
  double xlff(double const& x)
  {
    const double al = -0.963;
    const double bl = 1.37;
    const double Nl = 0.401 / ( std::tgamma(al + 2) * std::tgamma(bl + 1) / std::tgamma(al + bl + 3) );
    return Nl * pow(x, al + 1) * pow(1 - x, bl);
  }

  //_________________________________________________________________________________
  double xgff(double const& x)
  {
    const double ag = 1.943;
    const double bg = 8.000;
    const double Ng = 0.238 / ( std::tgamma(ag + 2) * std::tgamma(bg + 1) / std::tgamma(ag + bg + 3) );
    return Ng * pow(x, ag + 1) * pow(1 - x, bg);
  }

  //_________________________________________________________________________________
  std::map<int, double> Provider(double const& x, int const& isel)
  {
    // Initialise distributions
    double upv  = 0;
    double dnv  = 0;
    double glu  = 0;
    double dbar = 0;
    double ubar = 0;
    double sbar = 0;
    // Unpolarised PDFs
    if (isel == 0)
      {
        upv  = xupv(x);
        dnv  = xdnv(x);
        glu  = xglu(x);
        dbar = xdbar(x);
        ubar = xubar(x);
        sbar = xsbar(x);
      }
    // Polarised PDFs
    else if (isel == 1)
      {
        upv  = xupvpol(x);
        dnv  = xdnvpol(x);
        glu  = xglupol(x);
        dbar = xdbarpol(x);
        ubar = xubarpol(x);
        sbar = xsbarpol(x);
      }
    // Unpolarised FFs
    else if (isel == 2)
      {
        glu  = xgff(x);
        dbar = xlff(x);
        sbar = xsff(x);
        upv  = dbar - sbar;
        dnv  = sbar - dbar;
        ubar = sbar;
      }
    else
      throw std::runtime_error(error("Provider", "Unknown selector."));

    // Construct QCD evolution basis conbinations
    double const Gluon   = glu;
    double const Singlet = dnv + 2 * dbar + upv + 2 * ubar + 2 * sbar;
    double const T3      = upv + 2 * ubar - dnv - 2 * dbar;
    double const T8      = upv + 2 * ubar + dnv + 2 * dbar - 4 * sbar;
    double const Valence = upv + dnv;
    double const V3      = upv - dnv;

    // Fill in map in the QCD evolution basis
    std::map<int, double> QCDEvMap;
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
