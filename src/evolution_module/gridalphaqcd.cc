//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/gridalphaqcd.h"
#include "apfel/alphaqcd.h"

using namespace std;

namespace apfel {

  //_________________________________________________________________________________
  GridAlphaQCD::GridAlphaQCD(double         const& AlphaRef,
			     double         const& MuRef,
			     vector<double> const& Masses,
			     vector<double> const& Thresholds,
			     int            const& pt,
			     int            const& nQ,
			     double         const& QMin,
			     double         const& QMax,
			     int            const& InterDegree):
    QGrid<double>(nQ, QMin, QMax, InterDegree, Thresholds)
  {
    // Initialze alphas
    AlphaQCD as{AlphaRef, MuRef, Masses, Thresholds, pt};

    // Set number of steps of the RK algorith to 1
    as.SetNumberOfSteps(1);

    // Find the point on the QGrid right below MuRef
    const auto tQ = lower_bound(_Qg.begin(), _Qg.end(), as.GetMuRef()) - _Qg.begin() - 1;

    // Resize container
    _GridValues.resize(_Qg.size());

    // Loop on "_Qg" below "MuRef"
    for (auto iQ = tQ; iQ >= 0; iQ--)
      {
	auto a = as.Evaluate(_Qg[iQ]);
	_GridValues[iQ] = a;
	as.SetObjectRef(a);
	as.SetMuRef(_Qg[iQ]);
      }

    // Loop on "_Qg" above "MuRef"
    as.SetObjectRef(AlphaRef);
    as.SetMuRef(MuRef);
    for (auto iQ = tQ + 1; iQ < (int) _Qg.size(); iQ++)
      {
	auto a = as.Evaluate(_Qg[iQ]);
	_GridValues[iQ] = a;
	as.SetObjectRef(a);
	as.SetMuRef(_Qg[iQ]);
      }
  }

  //_________________________________________________________________________________
  GridAlphaQCD::GridAlphaQCD(double         const& AlphaRef,
			     double         const& MuRef,
			     vector<double> const& Masses,
			     int            const& pt,
			     int            const& nQ,
			     double         const& QMin,
			     double         const& QMax,
			     int            const& InterDegree):
  GridAlphaQCD::GridAlphaQCD(AlphaRef, MuRef, Masses, Masses, pt, nQ, QMin, QMax, InterDegree)
  {
  }

}
