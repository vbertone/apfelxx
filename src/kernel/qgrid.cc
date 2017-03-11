//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/qgrid.h"
#include "apfel/tools.h"
#include "apfel/distribution.h"
#include "apfel/set.h"

#include <algorithm>

namespace apfel
{
  //_________________________________________________________________________________
  template<class T>
  QGrid<T>::QGrid(int const& nQ, double const& QMin, double const& QMax, int const& InterDegree, vector<double> const& Thresholds, double const& Lambda):
    _nQ(nQ),
    _InterDegree(InterDegree),
    _QMin(QMin),
    _QMax(QMax),
    _Lambda(Lambda),
    _Thresholds(Thresholds)
  {
    // Check that QMin is actually smaller than QMax
    if (QMax <= QMin)
      throw runtime_exception("QGrid::QGrid","QMax must be larger than QMin");

    // Find initial and final number of flavours
    const auto nfin = ( _Thresholds.empty() ? 0 : lower_bound(_Thresholds.begin()+1, _Thresholds.end(), _QMin) - _Thresholds.begin() );
    const auto nffi = ( _Thresholds.empty() ? 0 : lower_bound(_Thresholds.begin()+1, _Thresholds.end(), _QMax) - _Thresholds.begin() );

    // Compute a temporary grid constant in ln(ln(Q^2/Lambda^2)
    // without taking into account the threholds
    vector<double> llql = { log( 2 * log( _QMin / _Lambda ) ) };
    auto Step = log( log( _QMax / _Lambda ) / log( _QMin / _Lambda ) ) / _nQ;
    for (auto iq = 1; iq <= _nQ; iq++) llql.push_back(llql.back()+Step);

    // Identify the indices of the grid points that fall right below the thresholds.
    // This will define the number of points that fall in each interval defined by
    // the grid bounds and the thresholds.
    _nQg.push_back(0);
    vector<double> llThL = {log( 2 * log( _QMin / _Lambda ) )};
    for (auto isg = nfin+1; isg <= nffi; isg++)
      {
	llThL.push_back(log( 2 * log( _Thresholds[isg-1] / _Lambda ) ));
	_nQg.push_back(lower_bound(llql.begin()+1, llql.end(), llThL.back()) - llql.begin());
      }
    _nQg.push_back(_nQ);
    llThL.push_back(log( 2 * log( _QMax / _Lambda ) ));

    // Check that all intervals have at least two points.
    // If not, adjust it.
    // Check also whether there is any subgrid whose number of points
    // is smaller than the interpolation degree plus one.
    // If so, adjust the interpolation degree.
    for (auto isg = 0; isg < (int) _nQg.size() - 1; isg++)
      {
	if (_nQg[isg+1] - _nQg[isg] < 2)                _nQg[isg+1]  = _nQg[isg] + 2;
	if (_nQg[isg+1] - _nQg[isg] < _InterDegree + 2) _InterDegree = _nQg[isg+1] - _nQg[isg] - 1;
      }

    // Adjust _nQ if needed
    if (_nQ != _nQg.back()) _nQ = _nQg.back();

    // Now construct the actual grid in such a way that the threshold concides
    // with nodes of the grid.
    _llQ2g.push_back(log( 2 * log( _QMin / _Lambda ) ));
    for (auto isg = 0; isg < (int) _nQg.size() - 1; isg++)
      {
	Step = ( llThL[isg+1] - llThL[isg] ) / ( _nQg[isg+1] - _nQg[isg] - 1 );
	for (auto iq = _nQg[isg] + 1; iq < _nQg[isg+1]; iq++) _llQ2g.push_back(_llQ2g.back()+Step);
	_llQ2g.push_back(_llQ2g.back());
      }

    // Displace slightly the values below and above the thresholds
    for (auto isg = 1; isg < (int) _nQg.size() - 1; isg++)
      {
	_llQ2g[_nQg[isg]-1] *= 1 - eps15;
	_llQ2g[_nQg[isg]]   *= 1 + eps15;
      }

    // Now compute grid in Q.
    for (auto const& lq : _llQ2g) _Qg.push_back(_Lambda * exp( exp(lq) / 2 ));
  }

  //_________________________________________________________________________________
  template<class T>
  double QGrid<T>::Interpolant(int const& tQ, int const& tau, double const& lnln2ql) const
  {
    // Return immediately 1 if "Q" coincides with "_Qg[tau]" unless tau coincides
    // with a threshold index. In that case return zero.
    if (fabs(lnln2ql / _llQ2g[tau] - 1) < eps11)
      return 1;

    // Define the lower bound of the interpolation range.
    int bound = tau + tQ - _InterDegree;
    if (_InterDegree > tau + tQ) bound = 0;

    // Initialize interpolant
    double w_int = 1;

    // Find the the neighbors of "Q" on the grid
    int j;
    for (j = 0; j <= tau+tQ-bound; j++)
      if (lnln2ql >= _llQ2g[tau+tQ-j] && lnln2ql < _llQ2g[tau+tQ-j+1])
	break;

    // Compute the interpolant
    for (auto delta = 0; delta <= _InterDegree; delta++)
      if(delta != j)
	w_int *= ( lnln2ql - _llQ2g[tau-j+delta] ) / ( _llQ2g[tau] - _llQ2g[tau-j+delta] );

    return w_int;
  }

  //_________________________________________________________________________________
  template<class T>
  tuple<int,int,int> QGrid<T>::SumBounds(double const& Q) const
  {
    tuple<int,int,int> bounds (0,0,0);

    // Return if "Q" is outside the grid range (no sum will be performed)
    if (Q < _Qg.front() - eps12 || Q > _Qg.back() + eps12)
      return bounds;

    // If Q falls in the tiny gap at the thresholds assume the node below
    for (auto iQ = 1; iQ < (int) _nQg.size() - 1; iQ++)
      if (Q > _Qg[_nQg[iQ]-1] && Q <= _Qg[_nQg[iQ]])
	{
	  get<1>(bounds) = _nQg[iQ] - 1;
	  get<2>(bounds) = _nQg[iQ];
	  return bounds;
	}
      
    // Identify the subgrid in which Q falls
    int iQ;
    for (iQ = 0; iQ < (int) _nQg.size() - 1; iQ++)
      if (Q > _Qg[_nQg[iQ]] && Q <= _Qg[_nQg[iQ+1]])
	break;

    // Determine the control parameter and put it in the first entry of the tuple
    if (Q > _Qg[_nQg[iQ+1]-_InterDegree])
      {
	int id;
	for (id = 2; id <= _InterDegree; id++)
	  if (Q > _Qg[_nQg[iQ+1]-id])
	    break;
	get<0>(bounds) = _InterDegree - id + 1;
      }

    // Determine the actual bounds
    const auto low = lower_bound(_Qg.begin()+1, _Qg.end(), Q) - _Qg.begin();
    get<1>(bounds) = low;
    get<2>(bounds) = low;

    if (fabs(Q / _Qg[low] - 1) <= eps12)
          get<2>(bounds) += 1;
    else
      {
        get<1>(bounds) += - get<0>(bounds) - 1;
        get<2>(bounds) += - get<0>(bounds) + _InterDegree;
      }

    return bounds;
  }

  //_________________________________________________________________________________
  template<class T>
  T QGrid<T>::Evaluate(double const& Q) const
  {
    auto const bounds = SumBounds(Q);
    auto const ll2ql  = log( 2 * log( Q / _Lambda ) );

    // first create a copy of template object with the first component
    auto tau = get<1>(bounds);
    T result = Interpolant(get<0>(bounds), tau, ll2ql) * _GridValues[tau];

    // then loop and add the extra terms
    for (tau = tau+1; tau < get<2>(bounds); tau++)
      result += Interpolant(get<0>(bounds), tau, ll2ql) * _GridValues[tau];

    return result;
  }

  //_________________________________________________________________________________
  template<class T>
  bool QGrid<T>::operator == (QGrid const& Qg) const
  {
    if(_nQ != Qg._nQ)                   return false;
    if(_QMin != Qg._QMin)               return false;
    if(_QMax != Qg._QMax)               return false;
    if(_InterDegree != Qg._InterDegree) return false;
    if(_Thresholds != Qg._Thresholds)   return false;
    if(_Lambda != Qg._Lambda)           return false;
    return true;
  }

  //_________________________________________________________________________________
  template<class T>
  bool QGrid<T>::operator != (QGrid const& Qg) const
  {
    if (*this == Qg) return false;
    else return true;
  }

  // template fixed types
  template class QGrid<double>;
  template class QGrid<Distribution>;
  template class QGrid<Set<Distribution>>;
}
