//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/tools.h"
#include "apfel/constants.h"

#include <math.h>
#include <algorithm>

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

  //_________________________________________________________________________
  std::vector<double> ElectroWeakCharges(double const& Q, bool const& virt, int const& sel)
  {
    // Relevant constants
    const double Q2    = Q * Q;
    const double MZ2   = apfel::ZMass * apfel::ZMass;
    const double GmZ2  = apfel::GammaZ * apfel::GammaZ;
    const double S2ThW = apfel::Sin2ThetaW;
    const double VD    = - 0.5 + 2 * S2ThW / 3;
    const double VU    = + 0.5 - 4 * S2ThW / 3;
    const double AD    = - 0.5;
    const double AU    = + 0.5;
    const double Ve    = - 0.5 + 2 * S2ThW;
    const double Ae    = - 0.5;
    const std::vector<double> Vq = {VD, VU, VD, VU, VD, VU};
    const std::vector<double> Aq = {AD, AU, AD, AU, AD, AU};

    // Propagator and its square
    double PZ;
    double PZ2;
    if (virt)
      {
	PZ  = Q2 * ( Q2 -  MZ2 ) / ( pow(Q2 - MZ2,2) + MZ2 * GmZ2 ) / ( 4 * S2ThW * ( 1 - S2ThW ) );
	PZ2 = pow(Q2,2) / ( pow(Q2 - MZ2,2) + MZ2 * GmZ2 ) / pow(4 * S2ThW * ( 1 - S2ThW ),2);
      }
    else
      {
	PZ  = Q2 / ( Q2 + MZ2 ) / ( 4 * S2ThW * ( 1 - S2ThW ) );
	PZ2 = PZ * PZ;
      }

    // Build electroweak charges
    std::vector<double> Charges(6, 0.);
    if (sel < 0 || sel > 5)
      for (auto i = 0; i < 6; i++)
	Charges[i] = apfel::QCh2[i]
	  - 2 * apfel::QCh[i] * Vq[i] * Ve * PZ
	  + ( Ve * Ve + Ae * Ae ) * ( Vq[i] * Vq[i] + Aq[i] * Aq[i] ) * PZ2;
    else
      Charges[sel] = apfel::QCh2[sel]
	- 2 * apfel::QCh[sel] * Vq[sel] * Ve * PZ
	+ ( Ve * Ve + Ae * Ae ) * ( Vq[sel] * Vq[sel] + Aq[sel] * Aq[sel] ) * PZ2;

    return Charges;
  }

  //_________________________________________________________________________________
  std::vector<double> ConcatenateAndSortVectors(std::vector<double> const& v1, std::vector<double> const& v2)
  {
    std::vector<double> v12 = v1;
    v12.insert(v12.end(), v2.begin(), v2.end());
    std::sort(v12.begin(), v12.end());
    return v12;
  }
}
