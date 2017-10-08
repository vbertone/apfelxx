//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/tabulateobject.h"
#include "apfel/distribution.h"
#include "apfel/operator.h"
#include "apfel/doubleobject.h"
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
				    int                 const& InterDegree,
				    double              const& Lambda):
    QGrid<T>(nQ, QMin, QMax, InterDegree, Object.GetThresholds(), Lambda)
  {
    // Save initial conditions.
    const auto nsteps = Object.GetNumberOfSteps();
    const auto ObjRef = Object.GetObjectRef();
    const auto MuRef  = Object.GetMuRef();

    // Set number of steps of the RK algorith to 1.
    Object.SetNumberOfSteps(1);

    // Find the point on the QGrid right below MuRef.
    const auto tQ = lower_bound(this->_Qg.begin(), this->_Qg.end(), MuRef) - this->_Qg.begin() - 1;

    // Loop on "_Qg" below "MuRef".
    for (auto iQ = tQ; iQ >= 0; iQ--)
      {
	auto o = Object.Evaluate(this->_Qg[iQ]);
	this->_GridValues.push_back(o);
	Object.SetObjectRef(o);
	Object.SetMuRef(this->_Qg[iQ]);
      }

    // Reverse order of the elements.
    reverse(this->_GridValues.begin(),this->_GridValues.end());

    // Loop on "_Qg" above "MuRef".
    Object.SetObjectRef(ObjRef);
    Object.SetMuRef(MuRef);
    for (auto iQ = tQ + 1; iQ < (int) this->_Qg.size(); iQ++)
      {
	auto o = Object.Evaluate(this->_Qg[iQ]);
	this->_GridValues.push_back(o);
	Object.SetObjectRef(o);
	Object.SetMuRef(this->_Qg[iQ]);
      }

    // Reset initial conditions.
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
				    vector<double>      const& Thresholds,
				    double              const& Lambda):
    QGrid<T>(nQ, QMin, QMax, InterDegree, Thresholds, Lambda)
  {
    // Fill in Qgrid with the object.
    for (auto const& iQ : this->_Qg)
      this->_GridValues.push_back(Object(iQ));
  }

  //_________________________________________________________________________________
  template<class T>
  TabulateObject<T>::TabulateObject(function<T(double)>             const& Object,
				    int                             const& nQ,
				    double                          const& QMin,
				    double                          const& QMax,
				    int                             const& InterDegree,
				    vector<double>                  const& Thresholds,
				    function<double(double const&)> const& TabFunc,
				    function<double(double const&)> const& InvTabFunc):
    QGrid<T>(nQ, QMin, QMax, InterDegree, Thresholds, TabFunc, InvTabFunc)
  {
    // Fill in Qgrid with the object.
    for (auto const& iQ : this->_Qg)
      this->_GridValues.push_back(Object(iQ));
  }

  // Specializations
  //_________________________________________________________________________________
  template class TabulateObject<double>;
  template class TabulateObject<Distribution>;
  template class TabulateObject<Set<Distribution>>;
  template class TabulateObject<DoubleObject<Distribution>>;
  template class TabulateObject<Operator>;
  template class TabulateObject<Set<Operator>>;

  //_________________________________________________________________________________
  template<>
  double TabulateObject<Distribution>::EvaluatexQ(double const& x, double const& Q) const
  {
    const auto fq     = this->_TabFunc(Q);
    const auto bounds = this->SumBounds(Q);

    // Loop over the nodes.
    double result = 0;
    for (auto tau = get<1>(bounds); tau < get<2>(bounds); tau++)
      result += Interpolant(get<0>(bounds), tau, fq) * this->_GridValues[tau].Evaluate(x);

    return result;
  };

  //_________________________________________________________________________________
  template<>
  double TabulateObject<Set<Distribution>>::EvaluatexQ(int const& i, double const& x, double const& Q) const
  {
    const auto fq     = this->_TabFunc(Q);
    const auto bounds = this->SumBounds(Q);

    // Loop over the nodes.
    double result = 0;
    for (auto tau = get<1>(bounds); tau < get<2>(bounds); tau++)
      result += Interpolant(get<0>(bounds), tau, fq) * this->_GridValues[tau].at(i).Evaluate(x);

    return result;
  }

  //_________________________________________________________________________________
  template<>
  double TabulateObject<DoubleObject<Distribution>>::EvaluatexzQ(double const& x, double const& z, double const& Q) const
  {
    const auto fq     = this->_TabFunc(Q);
    const auto bounds = this->SumBounds(Q);

    // Loop over the nodes.
    double result = 0;
    for (auto tau = get<1>(bounds); tau < get<2>(bounds); tau++)
      result += Interpolant(get<0>(bounds), tau, fq) * this->_GridValues[tau].Evaluate(x,z);

    return result;
  }

  //_________________________________________________________________________________
  template<>
  map<int,double> TabulateObject<Set<Distribution>>::EvaluateMapxQ(double const& x, double const& Q) const
  {
    const auto fq     = this->_TabFunc(Q);
    const auto bounds = this->SumBounds(Q);
    const int cp      = get<0>(bounds);
    const int lower   = get<1>(bounds);
    const int upper   = get<2>(bounds);

    // Fill in map.
    map<int,double> result;
    for (auto tau = lower; tau < upper; tau++)
      {
	const auto& obj = this->_GridValues[tau].GetObjects();
	const double w = Interpolant(cp, tau, fq);
	for (auto it = obj.begin(); it != obj.end(); ++it)
	  result[it->first] += w * it->second.Evaluate(x);
      }
    return result;
  }
}
