//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/gmatchingfunctions.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________________
  Cgtmd1ns::Cgtmd1ns(double const& xi):
    Expression(),
    _xi(xi)
  {
  }
  double Cgtmd1ns::Regular(double const& y) const
  {
    const double kappa = _xi / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 2 * CF * ( 1 - y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * ( 1 - kappa ) * y / ( 1 - ky2 ) : 0);
  }
  double Cgtmd1ns::Local(double const&) const
  {
    return - CF * zeta2;
  }

  //_________________________________________________________________________________
  Cgtmd1qq::Cgtmd1qq(double const& xi):
    Expression(),
    _xi(xi)
  {
  }
  double Cgtmd1qq::Regular(double const& y) const
  {
    const double kappa = _xi / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 2 * CF * ( 1 - y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * ( 1 - kappa ) / kappa / ( 1 - ky2 ) : 0);
  }
  double Cgtmd1qq::Local(double const&) const
  {
    return - CF * zeta2;
  }

  //_________________________________________________________________________________
  Cgtmd1qg::Cgtmd1qg(double const& xi):
    Expression(),
    _xi(xi)
  {
  }
  double Cgtmd1qg::Regular(double const& y) const
  {
    const double kappa = _xi / _extvar;
    const double ky2   = pow(kappa * y, 2);
    //return (y <= 1 ? 8 * TR * y * ( 1 - y ) / pow(1 - ky2, 2) : 0)
    //       + (kappa > 1 ? 8 * TR * ( 1 - kappa ) * pow(y, 2) / pow(1 - ky2, 2) : 0);

    // I multiplied this matching function by ( 1 - kappa * y )
    // because otherwise the divergence at y = 1/kappa would not
    // cancel and the integrals would diverge. This is probably a
    // signal of the fact that the specific definition of the gluon
    // GPD that is being used here vanishes at y = xi. This is visible
    // from Eq. (116) of https://arxiv.org/pdf/2206.01412.pdf.
    // Thefore, the additional factor ( 1 - kappa * y ) makes this
    // matching function useable with any gluon GPD that does not
    // vanish at y = xi and effectively amounts to replacing in the
    // convolution integral:
    //
    //  G(y, xi) -> y / ( y - xi ) G(y, xi)
    //
    // that is assumed to be finite at y = xi. The same applies to Cgg
    // below.
    return ( 1 - kappa * y ) * ( (y <= 1 ? 8 * TR * y * ( 1 - y ) / pow(1 - ky2, 2) : 0)
                                 + (kappa > 1 ? 8 * TR * ( 1 - kappa ) * pow(y, 2) / pow(1 - ky2, 2) : 0) );
  }

  //_________________________________________________________________________________
  Cgtmd1gq::Cgtmd1gq(double const& xi):
    Expression(),
    _xi(xi)
  {
  }
  double Cgtmd1gq::Regular(double const& y) const
  {
    const double kappa = _xi / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 2 * CF * ( 1 - pow(kappa, 2) ) * y / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? - 2 * CF * ( 1 - pow(kappa, 2) ) / kappa / ( 1 - ky2 ) : 0);
  }

  //_________________________________________________________________________________
  Cgtmd1gg::Cgtmd1gg(double const& xi):
    Expression(),
    _xi(xi)
  {
  }
  double Cgtmd1gg::Regular(double const& y) const
  {
    const double kappa = _xi / _extvar;
    const double ky2   = pow(kappa * y, 2);
    //return (y <= 1 ? 8 * CA * pow(kappa, 2) * y * ( 1 - y ) / pow(1 - ky2, 2) : 0)
    //       + (kappa > 1 ? CA * ( 1 - kappa ) * ( 1 + kappa - ( 1 - 7 * kappa ) * ky2 ) / kappa / pow(1 - ky2, 2) : 0);

    // See Cqg above
    return ( 1 - kappa * y ) * ( (y <= 1 ? 8 * CA * pow(kappa, 2) * y * ( 1 - y ) / pow(1 - ky2, 2) : 0)
                                 + (kappa > 1 ? CA * ( 1 - kappa ) * ( 1 + kappa - ( 1 - 7 * kappa ) * ky2 ) / kappa / pow(1 - ky2, 2) : 0) );
  }
  double Cgtmd1gg::Local(double const&) const
  {
    return - CA * zeta2;
  }
}
