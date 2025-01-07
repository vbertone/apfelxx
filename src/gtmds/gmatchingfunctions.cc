//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//         Simone Rodini: rodini.simone.luigi@gmail.com
//

#include "apfel/gmatchingfunctions.h"
#include "apfel/constants.h"

/*
* Structure of the matching functions. Columns correspond to GPDs and
* rows to GTMDs:
*
* |     | q,U | g,U | q,L | g,L | q,T | g,T |
* | --- | --- | --- | --- | --- | --- | --- |
* | q,U |  v  |  v  |     |     |     |  v  |
* | g,U |  v  |  v  |     |     |     |  v  |
* | q,L |     |     |  v  |  v  |     |  v  |
* | g,L |     |     |  v  |  v  |     |  v  |
* | q,T |     |     |     |     |     |     |
* | g,T |  v  |  v  |  v  |  v  |     |  v  |
*/

namespace apfel
{
  //_________________________________________________________________________________
  Cgtmd1nseUU::Cgtmd1nseUU(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1nseUU::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 2 * CF * ( 1 - y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * ( 1 - kappa ) * y / ( 1 - ky2 ) : 0);
  }
  double Cgtmd1nseUU::Local(double const&) const
  {
    return - CF * zeta2;
  }

  //_________________________________________________________________________________
  Cgtmd1qqeUU::Cgtmd1qqeUU(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1qqeUU::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 2 * CF * ( 1 - y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * ( 1 - kappa ) / kappa / ( 1 - ky2 ) : 0);
  }
  double Cgtmd1qqeUU::Local(double const&) const
  {
    return - CF * zeta2;
  }

  //_________________________________________________________________________________
  Cgtmd1qgeUU::Cgtmd1qgeUU(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1qgeUU::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? 8 * TR * y * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? - 2 * TR * ( 1 - ky ) / kappa / pow(1 + ky, 2) : 8 * TR * ( 1 - kappa ) * pow(y, 2) / pow(1 - ky2, 2));
  }
  double Cgtmd1qgeUU::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return 2 * TR / kappa / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1qgeUU::LocalLogPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return 2 * TR / pow(kappa, 2) * log(1 - kappa * y);
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1gqeUU::Cgtmd1gqeUU(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1gqeUU::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 2 * CF * ( 1 - pow(kappa, 2) ) * y / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? - 2 * CF * ( 1 - pow(kappa, 2) ) / kappa / ( 1 - ky2 ) : 0);
  }

  //_________________________________________________________________________________
  Cgtmd1ggeUU::Cgtmd1ggeUU(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1ggeUU::Regular(double const& y) const // R2 and Rhat in different but equivalent form compared to the notes
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? - 8 * CA * pow(kappa, 2) * y * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? - CA * ( 1 - 3 * pow(kappa, 2) + ( 1 + pow(kappa, 2) ) * ky ) / kappa / pow(1 + ky, 2) : - 2 * CA * ( 1 - kappa ) * ( 1 + kappa - ( 1 - 3 * kappa ) * ky2 ) / kappa / pow(1 - ky2, 2));
  }
  double Cgtmd1ggeUU::Local(double const&) const
  {
    return - CA * zeta2;
  }
  double Cgtmd1ggeUU::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return - CA * ( 1 + pow(kappa, 2) ) / kappa / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1ggeUU::LocalLogPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return - CA * ( 1 + pow(kappa, 2) ) / pow(kappa, 2) * log(1 - kappa * y);
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1qgoUU::Cgtmd1qgoUU(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1qgoUU::LocalPV() const // From (3.11) master formula this includes: -\pi/\kappa. 'is' is kept externally
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1)
      return -2 * TR / pow(kappa, 2) * M_PI;
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1ggoUU::Cgtmd1ggoUU(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1ggoUU::LocalPV() const // From (3.11) master formula this includes: -\pi/\kappa. 'is' is kept externally
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1)
      return CA * ( 1 + pow(kappa, 2) ) / pow(kappa, 2) * M_PI;
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1qgeUT::Cgtmd1qgeUT(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1qgeUT::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? 8 * TR * y * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? - 2 * TR * ( 1 - ky ) / kappa / pow(1 + ky, 2) : 8 * TR * ( 1 - kappa ) * pow(y, 2) / pow(1 - ky2, 2));
  }
  double Cgtmd1qgeUT::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return 2 * TR / kappa / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1qgeUT::LocalLogPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return 2 * TR / pow(kappa, 2) * log(1 - kappa * y);
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1qgoUT::Cgtmd1qgoUT(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1qgoUT::LocalPV() const // From (3.11) master formula this includes: -\pi/\kappa. 'is' is kept externally
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1)
      return -2 * TR / pow(kappa, 2) * M_PI;
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1ggeUT::Cgtmd1ggeUT(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1ggeUT::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? - 4 * CA * (1 + pow(kappa, 2)) * y * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? CA * (1 + pow(kappa, 2)) * (1 - ky) / kappa / pow(1 + ky, 2) : - 4 * CA * pow(y, 2) * (1 - kappa) * (1 + pow(kappa, 2)) / pow(1 - ky2, 2));
  }
  double Cgtmd1ggeUT::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return - CA * ( 1 + pow(kappa, 2) ) / kappa / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1ggeUT::LocalLogPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return - CA * ( 1 + pow(kappa, 2) ) / pow(kappa, 2) * log(1 - kappa * y);
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1ggoUT::Cgtmd1ggoUT(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1ggoUT::LocalPV() const // From (3.11) master formula this includes: -\pi/\kappa. 'is' is kept externally
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1)
      return CA * ( 1 + pow(kappa, 2) ) / pow(kappa, 2) * M_PI;
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1nseLL::Cgtmd1nseLL(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1nseLL::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 2 * CF * ( 1 - y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * ( 1 - kappa ) / kappa / ( 1 - ky2 ) : 0);
  }
  double Cgtmd1nseLL::Local(double const&) const
  {
    return - CF * zeta2;
  }

  //_________________________________________________________________________________
  Cgtmd1qqeLL::Cgtmd1qqeLL(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1qqeLL::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 2 * CF * ( 1 - y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * ( 1 - kappa ) * y / ( 1 - ky2 ) : 0);
  }
  double Cgtmd1qqeLL::Local(double const&) const
  {
    return - CF * zeta2;
  }

  //_________________________________________________________________________________
  Cgtmd1qgeLL::Cgtmd1qgeLL(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1qgeLL::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? 8 * TR * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? 2 * TR * ( 3 + ky ) / pow(1 + ky, 2) : 8 * TR * ( 1 - kappa ) * y / pow(1 - ky2, 2));
  }
  double Cgtmd1qgeLL::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return 2 * TR / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1qgeLL::LocalLogPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return 2 * TR * log(1 - kappa * y) / kappa;
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1gqeLL::Cgtmd1gqeLL(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1gqeLL::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? - 4 * CF * ( 1 - y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? - 4 * CF * ( 1 - kappa ) * y / ( 1 - ky2 ) : 0);
  }

  //_________________________________________________________________________________
  Cgtmd1ggeLL::Cgtmd1ggeLL(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1ggeLL::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? - 8 * CA * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? - 2 * CA * ( 3 + ky ) / pow(1 + ky, 2) : - 8 * CA * ( 1 - kappa ) * y / pow(1 - ky2, 2));
  }
  double Cgtmd1ggeLL::Local(double const&) const
  {
    return - CA * zeta2;
  }
  double Cgtmd1ggeLL::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return - 2 * CA / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1ggeLL::LocalLogPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return - 2 * CA * log(1 - kappa * y) / kappa;
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1qgoLL::Cgtmd1qgoLL(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1qgoLL::LocalPV() const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1)
      return - 2 * TR / kappa * M_PI;
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1ggoLL::Cgtmd1ggoLL(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1ggoLL::LocalPV() const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1)
      return 2 * CA / kappa * M_PI;
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1qgeLT::Cgtmd1qgeLT(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1qgeLT::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? 8 * TR * ky * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? - 2 * TR * ( 1 - ky ) / kappa / pow(1 + ky, 2) : 8 * TR * ( 1 - kappa ) * y / pow(1 - ky2, 2));
  }
  double Cgtmd1qgeLT::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return 2 * TR / kappa / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1qgeLT::LocalLogPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return 2 * TR * log(1 - kappa * y) / pow(kappa, 2);
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1ggeLT::Cgtmd1ggeLT(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1ggeLT::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? - 8 * CA * ky * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? 2 * CA * ( 1 - ky ) / kappa / pow(1 + ky, 2) : - 8 * CA * ( 1 - kappa ) * y / pow(1 - ky2, 2));
  }
  double Cgtmd1ggeLT::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return - 2 * CA / kappa / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1ggeLT::LocalLogPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return - 2 * CA * log(1 - kappa * y) / pow(kappa, 2);
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1qgoLT::Cgtmd1qgoLT(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1qgoLT::LocalPV() const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1)
      return - 2 * TR / pow(kappa, 2) * M_PI;
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1ggoLT::Cgtmd1ggoLT(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1ggoLT::LocalPV() const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1)
      return 2 * CA / pow(kappa, 2) * M_PI;
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1ggeTT::Cgtmd1ggeTT(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1ggeTT::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? - 8 * CA * pow(kappa, 2) * y * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? 2 * CA * kappa * ( 1 - ky ) / pow(1 + ky, 2) : - 8 * CA * ky2 * ( 1 - kappa ) / pow(1 - ky2, 2));
  }
  double Cgtmd1ggeTT::Local(double const&) const
  {
    return - CA * zeta2;
  }
  double Cgtmd1ggeTT::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return - 2 * CA * kappa / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1ggeTT::LocalLogPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return - 2 * CA * log(1 - kappa * y);
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1ggoTT::Cgtmd1ggoTT(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1ggoTT::LocalPV() const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1)
      return 2 * CA * M_PI;
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1gqeTU::Cgtmd1gqeTU(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1gqeTU::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? - 4 * CF * ( 1 - y ) / y / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? - 4 * CF * ( 1 - kappa ) / ( 1 - ky2 ) : 0);
  }

  //_________________________________________________________________________________
  Cgtmd1ggeTU::Cgtmd1ggeTU(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1ggeTU::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? - 4 * CA * ( 1 - y ) * ( 1 + ky2 ) / y / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? CA * ( 2 * kappa / (1 + ky) + 4 * kappa / pow(1 + ky, 2) - 4 / y ) : - 4 * CA * ( 1 - kappa ) * ( 1 + ky2 ) / pow(1 - ky2, 2));
  }
  double Cgtmd1ggeTU::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return - 2 * CA * kappa / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1ggeTU::LocalLogPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return - 2 * CA * log(1 - kappa * y);
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1ggoTU::Cgtmd1ggoTU(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1ggoTU::LocalPV() const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1)
      return 2 * CA * M_PI;
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1gqeTL::Cgtmd1gqeTL(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1gqeTL::Regular(double const& y) const // We changed the definitions from very first iteration, see tables in draft.
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 4 * CF * ( 1 - y ) * kappa / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 4 * CF * ( 1 - kappa ) / ( 1 - ky2 ) : 0);
  }

  //_________________________________________________________________________________
  Cgtmd1ggeTL::Cgtmd1ggeTL(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1ggeTL::Regular(double const& y) const // same as for gqeTL
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? 8 * CA * kappa * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? 2 * CA * ( 3 + ky ) / pow(1 + ky, 2) : 8 * CA * ( 1 - kappa ) / pow(1 - ky2, 2));
  }
  double Cgtmd1ggeTL::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return 2 * CA / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1ggeTL::LocalLogPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return 2 * CA * log(1 - kappa * y) / kappa;
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1ggoTL::Cgtmd1ggoTL(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1ggoTL::LocalPV() const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1)
      return  - 2 * CA / kappa * M_PI;
    else
      return 0;
  }
}
