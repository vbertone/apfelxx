//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/tabulateobject.h"
#include "apfel/distribution.h"
#include "apfel/set.h"

#include <algorithm>

using namespace std;

namespace apfel {

  //_________________________________________________________________________________
  template<class T>
  TabulateObject<T>::TabulateObject(MatchedEvolution<T> &Object,
				    int                 const& nQ,
				    double              const& QMin,
				    double              const& QMax,
				    int                 const& InterDegree):
    QGrid<T>(nQ, QMin, QMax, InterDegree, Object.GetThresholds())
  {
    // Save initial conditions
    const auto nsteps = Object.GetNumberOfSteps();
    const auto ObjRef = Object.GetObjectRef();
    const auto MuRef  = Object.GetMuRef();

    // Set number of steps of the RK algorith to 1
    Object.SetNumberOfSteps(1);

    // Find the point on the QGrid right below MuRef
    const auto tQ = lower_bound(this->_Qg.begin(), this->_Qg.end(), MuRef) - this->_Qg.begin() - 1;

    // Loop on "_Qg" below "MuRef"
    for (auto iQ = tQ; iQ >= 0; iQ--)
      {
	auto o = Object.Evaluate(this->_Qg[iQ]);
	this->_GridValues.push_back(o);
	Object.SetObjectRef(o);
	Object.SetMuRef(this->_Qg[iQ]);
      }

    // Reverse order of the elements
    reverse(this->_GridValues.begin(),this->_GridValues.end());

    // Loop on "_Qg" above "MuRef"
    Object.SetObjectRef(ObjRef);
    Object.SetMuRef(MuRef);
    for (auto iQ = tQ + 1; iQ < (int) this->_Qg.size(); iQ++)
      {
	auto o = Object.Evaluate(this->_Qg[iQ]);
	this->_GridValues.push_back(o);
	Object.SetObjectRef(o);
	Object.SetMuRef(this->_Qg[iQ]);
      }

    // Reset initial conditions
    Object.SetNumberOfSteps(nsteps);
    Object.SetObjectRef(ObjRef);
    Object.SetMuRef(MuRef);
  }

  template class TabulateObject<double>;
  template class TabulateObject<Distribution>;
  template class TabulateObject<Set<Distribution>>;

}
