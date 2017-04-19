//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/tabulateobject.h"
#include "apfel/distribution.h"
#include "apfel/set.h"
#include "apfel/tools.h"

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

  //_________________________________________________________________________________
  template<class T>
  TabulateObject<T>::TabulateObject(function<T(double)> const& Object,
				    int                 const& nQ,
				    double              const& QMin,
				    double              const& QMax,
				    int                 const& InterDegree,
				    vector<double>      const& Thresholds):
    QGrid<T>(nQ, QMin, QMax, InterDegree, Thresholds)
  {
    for (auto const& iQ : this->_Qg)
      this->_GridValues.push_back(Object(iQ));
  }

  // Specializations
  //_________________________________________________________________________________
  template class TabulateObject<double>;
  template class TabulateObject<Distribution>;
  template class TabulateObject<Set<Distribution>>;

  //_________________________________________________________________________________
  template<>
  double TabulateObject<double>::EvaluatexQ(double const&, double const&) const
  {
    throw runtime_exception("TabulateObject::EvaluatexQ(x,Q)","This function can't be used for the specialization 'double' of the TabulateObject class.");
  }

  template<>
  double TabulateObject<Distribution>::EvaluatexQ(double const& x, double const& Q) const
  {
    cout << ", x = " << x << ", Q = " << Q << endl;
    return 0;
  };

  template<>
  double TabulateObject<Set<Distribution>>::EvaluatexQ(double const&, double const&) const
  {
    throw runtime_exception("TabulateObject::EvaluatexQ(x,Q)","This function can't be used for the specialization 'Set<Distribution>' of the TabulateObject class.");
  }

  //_________________________________________________________________________________
  template<>
  double TabulateObject<double>::EvaluatexQ(int const&, double const&, double const&) const
  {
    throw runtime_exception("TabulateObject::EvaluatexQ(i,x,Q)","This function can't be used for the specialization 'double' of the TabulateObject class.");
  }

  template<>
  double TabulateObject<Distribution>::EvaluatexQ(int const&, double const&, double const&) const
  {
    throw runtime_exception("TabulateObject::EvaluatexQ(i,x,Q)","This function can't be used for the specialization 'Distribution' of the TabulateObject class.");
  }

  template<>
  double TabulateObject<Set<Distribution>>::EvaluatexQ(int const& i, double const& x, double const& Q) const
  {
    const auto ll2ql  = log( 2 * log( Q / this->_Lambda ) );
    const auto bounds = this->SumBounds(Q);

    // Loop over the nodes
    double result = 0;
    for (auto tau = get<1>(bounds); tau < get<2>(bounds); tau++)
      result += Interpolant(get<0>(bounds), tau, ll2ql) * this->_GridValues[tau].at(i).Evaluate(x);

    return result;
  }

}
