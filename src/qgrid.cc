//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/qgrid.h"
#include "apfel/tools.h"

namespace apfel
{

  //_________________________________________________________________________________
  QGrid::QGrid(int const& nQ, double const& QMin, double const& QMax, int const& InterDegree, vector<double> const& Thresholds, double const& Lambda):
    _nQ(nQ),
    _InterDegree(InterDegree),
    _QMin(QMin),
    _QMax(QMax),
    _Thresholds(Thresholds),
    _Lambda(Lambda)
  {
    // Check that QMin is actually smaller than QMax
    if (QMax <= QMin)
      throw runtime_exception("QGrid::QGrid","QMax must be larger than QMin");

    // Find initial and final number of flavours
    const auto nfin = lower_bound(_Thresholds.begin()+1, _Thresholds.end(), _QMin) - _Thresholds.begin();
    const auto nffi = lower_bound(_Thresholds.begin()+1, _Thresholds.end(), _QMax) - _Thresholds.begin();

    // Check that the number of points of the grid is at least equal
    // to the number of thresholds falling in the grid interval plus
    // two given by the grid intervals.
    if (_nQ < nffi - nfin + 2)
      throw runtime_exception("QGrid::QGrid","Number of point of the Q-grid too small.");

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
    for (auto isg = 0; isg < (int) _nQg.size() - 1; isg++)
	if (_nQg[isg] == _nQg[isg+1]) _nQg[isg]--;

    // Now construct the actual grid in such a qay that the threshold concides
    // with nodes of the grid.
    _llQ2g.push_back(log( 2 * log( _QMin / _Lambda ) ));
    for (auto isg = 0; isg < (int) _nQg.size() - 1; isg++)
      {
	Step = ( llThL[isg+1] - llThL[isg] ) / ( _nQg[isg+1] - _nQg[isg] );
	for (auto iq = _nQg[isg]; iq < _nQg[isg+1]; iq++) _llQ2g.push_back(_llQ2g.back()+Step);
      }

    // Now compute grid in Q.
    for (auto const& lq : _llQ2g) _Qg.push_back(_Lambda * exp( exp(lq) / 2 ));
  }

  //_________________________________________________________________________________
  bool QGrid::operator == (QGrid const& Qg) const
  {
    if(_nQ != Qg._nQ)                   return false;
    if(_QMin != Qg._QMin)               return false;
    if(_QMax != Qg._QMax)               return false;
    if(_InterDegree != Qg._InterDegree) return false;
    if(_Thresholds != Qg._Thresholds)   return false;
    if(_Lambda != Qg._Lambda)   return false;
    return true;
  }

  //_________________________________________________________________________________
  bool QGrid::operator != (QGrid const& Qg) const
  {
    if (*this == Qg) return false;
    else return true;
  }

  //_________________________________________________________________________________
  std::ostream& operator<<(std::ostream& os, QGrid const& Qg)
  {
    os << "QGrid: " << &Qg << "\n";
    os << "nQ          = " << Qg._nQ << "\n";
    os << "QMin        = " << Qg._QMin << "\n";
    os << "QMax        = " << Qg._QMax << "\n";
    os << "InterDegree = " << Qg._InterDegree << "\n";
    os << "Lambda      = " << Qg._Lambda << "\n";
    os << "Thresholds  = ";
    for (const auto &v: Qg._Thresholds) os << v << " ";
    os << "\n";
    os << "Qg          = ";
    for (const auto &v: Qg._Qg) os << v << " ";
    os << "\n\n";
    return os;
  }

}
