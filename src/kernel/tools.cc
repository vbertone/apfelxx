//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/tools.h"

#include <math.h>

namespace apfel
{
  //_________________________________________________________________________
  int NF(double const& Q, std::vector<double> const& Thresholds)
  {
    // Compute number of active flavours the the PDF initial scale
    int nf = 0;
    for (auto const& v : Thresholds)
      if (Q > v)
	nf++;
      else
	break;
    return nf;
  }

  //_________________________________________________________________________
  double DeltaFun(double const& a, double const& b, double const& c)
  {
    return sqrt( a * a + b * b + c * c - 2 * ( a * b + b * c + c * a ) );
  }
}
