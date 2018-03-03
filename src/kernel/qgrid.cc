//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/qgrid.h"
#include "apfel/constants.h"
#include "apfel/messages.h"
#include "apfel/operator.h"
#include "apfel/set.h"
#include "apfel/doubleobject.h"

#include <stdexcept>

namespace apfel
{
  //_________________________________________________________________________________
  template<class T>
  QGrid<T>::QGrid(int                                  const& nQ,
		  double                               const& QMin,
		  double                               const& QMax,
		  int                                  const& InterDegree,
		  std::vector<double>                  const& Thresholds,
		  std::function<double(double const&)> const& TabFunc,
		  std::function<double(double const&)> const& InvTabFunc):
    _nQ(nQ),
    _QMin(QMin),
    _QMax(QMax),
    _InterDegree(InterDegree),
    _Thresholds(Thresholds),
    _TabFunc(TabFunc)
  {
    // Check that QMin is actually smaller than QMax.
    if (QMax <= QMin)
      throw std::runtime_error(error("QGrid::QGrid","QMax must be larger than QMin"));

    // Check that "TabFunc" and "InvTabFunc" are actually the inverse
    // function of each other. The check is done at grid bounds and in
    // the middle.

    const std::vector<double> TestPoints{_QMin, ( _QMax +  _QMin ) / 2, _QMax};
    for (auto const& p : TestPoints)
      {
	const double reldiff = abs(InvTabFunc( TabFunc(p) ) / p - 1);
	if (reldiff > eps8)
	  throw std::runtime_error(error("QGrid::QGrid","TabFunc and InvTabFunc are not the inverse of each other."));
      }

    // Find initial and final number of flavours.
    const int nfin = NF(_QMin, _Thresholds);
    const int nffi = NF(_QMax, _Thresholds);

    // Compute a temporary grid constant in ln(ln(Q^2/Lambda^2)
    // without taking into account the threholds.
    std::vector<double> fq = {_TabFunc(_QMin)};
    double Step = ( _TabFunc(_QMax) - _TabFunc(_QMin) ) / _nQ;
    for (int iq = 1; iq <= _nQ; iq++)
      fq.push_back(fq.back()+Step);

    // Identify the indices of the grid points that fall right below
    // the thresholds.  This will define the number of points that
    // fall in each interval defined by the grid bounds and the
    // thresholds.
    _nQg.push_back(0);
    std::vector<double> fqTh = {_TabFunc(_QMin)};
    for (int isg = nfin+1; isg <= nffi; isg++)
      {
	fqTh.push_back(_TabFunc(_Thresholds[isg-1]));
	_nQg.push_back(lower_bound(fq.begin()+1, fq.end(), fqTh.back()) - fq.begin());
      }
    _nQg.push_back(_nQ);
    fqTh.push_back(_TabFunc(_QMax));

    // Check that all intervals have at least two points. If not,
    // adjust it. Check also whether there is any subgrid whose number
    // of points is smaller than the interpolation degree plus one.
    // If so, adjust the interpolation degree.
    for (int isg = 0; isg < (int) _nQg.size() - 1; isg++)
      {
	if (_nQg[isg+1] - _nQg[isg] < 2)
	  _nQg[isg+1]  = _nQg[isg] + 2;
	if (_nQg[isg+1] - _nQg[isg] < _InterDegree + 2)
	  _InterDegree = _nQg[isg+1] - _nQg[isg] - 1;
      }

    // Adjust _nQ if needed.
    if (_nQ != _nQg.back()) _nQ = _nQg.back();

    // Now construct the actual grid in such a way that the threshold
    // concides with nodes of the grid.
    _fQg.push_back(_TabFunc(_QMin));
    for (int isg = 0; isg < (int) _nQg.size() - 1; isg++)
      {
	Step = ( fqTh[isg+1] - fqTh[isg] ) / ( _nQg[isg+1] - _nQg[isg] - 1 );
	for (int iq = _nQg[isg] + 1; iq < _nQg[isg+1]; iq++) _fQg.push_back(_fQg.back()+Step);
	_fQg.push_back(_fQg.back());
      }

    // Now compute grid in Q.
    for (auto const& lq : _fQg)
      _Qg.push_back(InvTabFunc(lq));

    // Displace slightly the values below and above the thresholds.
    for (int isg = 1; isg < (int) _nQg.size() - 1; isg++)
      {
	_Qg[_nQg[isg]-1] *= 1 - eps12;
	_Qg[_nQg[isg]]   *= 1 + eps12;
	_fQg[_nQg[isg]-1] = TabFunc(_Qg[_nQg[isg]-1]);
	_fQg[_nQg[isg]]   = TabFunc(_Qg[_nQg[isg]]);
      }
  }

  //_________________________________________________________________________________
  template<class T>
  QGrid<T>::QGrid(int                 const& nQ,
		  double              const& QMin,
		  double              const& QMax,
		  int                 const& InterDegree,
		  std::vector<double> const& Thresholds,
		  double              const& Lambda):
    QGrid<T>(nQ, QMin, QMax, InterDegree, Thresholds,
	     [Lambda] (double const& Q)->double{ return log( 2 * log( Q / Lambda ) ); },
	     [Lambda] (double const& fQ)->double{ return Lambda * exp( exp(fQ) / 2 ); })
  {
  }

  //_________________________________________________________________________________
  template<class T>
  double QGrid<T>::Interpolant(int const& tQ, int const& tau, double const& fq) const
  {
    // Return immediately 1 if "Q" coincides with "_Qg[tau]" unless
    // tau coincides with a threshold index. In that case return zero.
    if (abs(fq / _fQg[tau] - 1) < eps11)
      return 1;

    // Define the lower bound of the interpolation range.
    int bound = tau + tQ - _InterDegree;
    if (_InterDegree > tau + tQ)
      bound = 0;
    //if (fq < _fQg[bound] || fq >= _fQg[tau+tQ+1])
    //  return 0;

    // Initialize interpolant
    double w_int = 1;

    // Find the the neighbors of "Q" on the grid
    int j;
    for (j = tau+tQ-bound; j >=0; j--)
      if (fq < _fQg[tau+tQ-j+1])
	break;

    // Compute the interpolant
    for (int delta = tau-j; delta <= tau-j+_InterDegree; delta++)
      if (delta != tau)
	w_int *= ( fq - _fQg[delta] ) / ( _fQg[tau] - _fQg[delta] );

    return w_int;
  }

  //_________________________________________________________________________________
  template<class T>
  std::tuple<int,int,int> QGrid<T>::SumBounds(double const& Q) const
  {
    std::tuple<int,int,int> bounds{0, 0, 0};

    // Return if "Q" is outside the grid range (no sum will be
    // performed).
    if (Q < _Qg.front() - eps12 || Q > _Qg.back() + eps12)
      return bounds;

    // If Q falls in the tiny gap at the thresholds assume the node
    // below.
    for (int iQ = 1; iQ < (int) _nQg.size() - 1; iQ++)
      if (Q > _Qg[_nQg[iQ]-1] && Q <= _Qg[_nQg[iQ]])
	{
	  std::get<1>(bounds) = _nQg[iQ] - 1;
	  std::get<2>(bounds) = _nQg[iQ];
	  return bounds;
	}
      
    // Identify the subgrid in which Q falls.
    int iQ;
    for (iQ = 0; iQ < (int) _nQg.size() - 1; iQ++)
      if (Q > _Qg[_nQg[iQ]] && Q <= _Qg[_nQg[iQ+1]])
	break;

    if (iQ == (int) _nQg.size()-1)
      iQ--;

    // Determine the control parameter and put it in the first entry
    // of the tuple.
    if (Q > _Qg[_nQg[iQ+1] - _InterDegree])
      {
	int id;
	for (id = 2; id <= _InterDegree; id++)
	  if (Q > _Qg[_nQg[iQ+1] - id])
	    break;
	std::get<0>(bounds) = _InterDegree - id + 1;
      }

    // Determine the actual bounds.
    const int low       = lower_bound(_Qg.begin()+1, _Qg.end(), Q) - _Qg.begin();
    std::get<1>(bounds) = low;
    std::get<2>(bounds) = low;

    if (abs(Q / _Qg[low] - 1) <= eps12)
      std::get<2>(bounds) += 1;
    else
      {
        std::get<1>(bounds) += - std::get<0>(bounds) - 1;
        std::get<2>(bounds) += - std::get<0>(bounds) + _InterDegree;
      }

    return bounds;
  }

  //_________________________________________________________________________________
  template<class T>
  T QGrid<T>::Evaluate(double const& Q) const
  {
    const std::tuple<int,int,int> bounds = SumBounds(Q);
    const double                  fq     = _TabFunc(Q);

    // First create a copy of template object with the first
    // component...
    int tau  = std::get<1>(bounds);
    T result = Interpolant(std::get<0>(bounds), tau, fq) * _GridValues[tau];

    // ...then loop and add the extra terms.
    for (tau = tau+1; tau < std::get<2>(bounds); tau++)
      result += Interpolant(std::get<0>(bounds), tau, fq) * _GridValues[tau];

    return result;
  }

  //_________________________________________________________________________________
  template<class T>
  bool QGrid<T>::operator == (QGrid const& Qg) const
  {
    if (_nQ != Qg._nQ)
      return false;
    if (_QMin != Qg._QMin)
      return false;
    if (_QMax != Qg._QMax)
      return false;
    if (_InterDegree != Qg._InterDegree)
      return false;
    if (_Thresholds != Qg._Thresholds)
      return false;

    return true;
  }

  //_________________________________________________________________________________
  template<class T>
  bool QGrid<T>::operator != (QGrid const& Qg) const
  {
    if (*this == Qg)
      return false;
    else
      return true;
  }

  // Fixed template types.
  template class QGrid<double>;
  template class QGrid<Distribution>;
  template class QGrid<Set<Distribution>>;
  template class QGrid<DoubleObject<Distribution>>;
  template class QGrid<Operator>;
  template class QGrid<Set<Operator>>;
}
