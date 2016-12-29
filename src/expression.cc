//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/expression.h"

namespace apfel
{
  //_________________________________________________________________________
  Expression::Expression():
    _MassIndex(1),
    _FlavNumb(0)
  {
  }

  //_________________________________________________________________________
  Expression::Expression(double const& MassIndex):
    _MassIndex(MassIndex),
    _FlavNumb(0)
  {
  }

  //_________________________________________________________________________
  Expression::Expression(int const& FlavNumb):
    _MassIndex(1),
    _FlavNumb(FlavNumb)
  {
  }

}
