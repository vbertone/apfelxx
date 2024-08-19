//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/qgrid.h"
#include "apfel/constants.h"
#include "apfel/messages.h"
#include "apfel/operator.h"
#include "apfel/set.h"
#include "apfel/doubleobject.h"
#include "apfel/doubledistribution.h"
#include "apfel/doubleoperator.h"

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
    // Check that QMin is actually smaller than QMax
    if (QMax <= QMin)
      throw std::runtime_error(error("QGrid::QGrid", "QMax must be larger than QMin"));

    // Check that "TabFunc" and "InvTabFunc" are actually the inverse
    // function of each other. The check is done at grid bounds and in
    // the middle.
    const std::vector<double> TestPoints{_QMin, ( _QMax +  _QMin ) / 2, _QMax};
    for (auto const& p : TestPoints)
      {
        const double reldiff = std::abs(InvTabFunc( TabFunc(p) ) / p - 1);
        if (reldiff > eps8)
          throw std::runtime_error(error("QGrid::QGrid", "TabFunc and InvTabFunc are not the inverse of each other."));
      }

    // Find initial and final number of flavours
    const int nfin = NF(_QMin, _Thresholds);
    const int nffi = NF(_QMax, _Thresholds);

    // Compute a temporary grid constant in '_TabFunc(Q)' without
    // taking into account the threholds.
    std::vector<double> fq = {_TabFunc(_QMin)};
    double Step = ( _TabFunc(_QMax) - _TabFunc(_QMin) ) / _nQ;
    for (int iq = 1; iq <= _nQ; iq++)
      fq.push_back(fq.back() + Step);

    // Identify the indices of the grid points that fall right below
    // the thresholds. This will define the number of points that fall
    // in each interval defined by the grid bounds and the thresholds.
    _nQg.push_back(0);
    std::vector<double> fqTh = {_TabFunc(_QMin)};
    for (int isg = nfin + 1; isg <= nffi; isg++)
      {
        fqTh.push_back(_TabFunc(_Thresholds[isg-1]));
        _nQg.push_back(std::lower_bound(fq.begin()+1, fq.end(), fqTh.back()) - fq.begin());
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
          _nQg[isg+1] = _nQg[isg] + 2;
        if (_nQg[isg+1] - _nQg[isg] < _InterDegree + 2)
          _InterDegree = _nQg[isg+1] - _nQg[isg] - 1;
      }

    // Adjust _nQ if needed
    if (_nQ != _nQg.back())
      _nQ = _nQg.back();

    // Now construct the actual grid in such a way that the threshold
    // concides with nodes of the grid.
    _fQg.push_back(_TabFunc(_QMin));
    for (int isg = 0; isg < (int) _nQg.size() - 1; isg++)
      {
        Step = ( fqTh[isg+1] - fqTh[isg] ) / ( _nQg[isg+1] - _nQg[isg] - 1 );
        for (int iq = _nQg[isg] + 1; iq < _nQg[isg+1]; iq++)
          _fQg.push_back(_fQg.back()+Step);
        _fQg.push_back(_fQg.back());
      }

    // Now compute grid in Q
    for (auto const& lq : _fQg)
      _Qg.push_back(InvTabFunc(lq));

    // Displace slightly the values below and above the thresholds
    for (int isg = 1; isg < (int) _nQg.size() - 1; isg++)
      {
        _Qg[_nQg[isg] - 1]  *= 1 - eps12;
        _Qg[_nQg[isg]]      *= 1 + eps12;
        _fQg[_nQg[isg] - 1]  = TabFunc(_Qg[_nQg[isg] - 1]);
        _fQg[_nQg[isg]]      = TabFunc(_Qg[_nQg[isg]]);
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
             [Lambda] (double const& Q) -> double{ return log( 2 * log( Q / Lambda ) ); },
             [Lambda] (double const& fQ) -> double{ return Lambda * exp( exp(fQ) / 2 ); })
  {
  }

  //_________________________________________________________________________________
  template<class T>
  QGrid<T>::QGrid(std::vector<double> const& Qg,
                  int                 const& InterDegree):
    _nQ(Qg.size() - 1),
    _QMin(Qg.front()),
    _QMax(Qg.back()),
    _InterDegree(InterDegree),
    _Thresholds{},
    _TabFunc([] (double const& Q) -> double{ return Q; }),
    _Qg(Qg),
    _fQg(Qg),
    _nQg{0, _nQ}
  {
  }

  //_________________________________________________________________________________
  template<class T>
  double QGrid<T>::Interpolant(int const& tQ, int const& tau, double const& fq) const
  {
    // Return immediately 1 if "Q" coincides with "_Qg[tau]"
    if (std::abs(fq / _fQg[tau] - 1) < eps11)
      return 1;

    // Define the lower bound of the interpolation range
    const int bound = std::max(tau + tQ - _InterDegree, 0);

    // Return zero if fq is outside the allowed range
    if (fq < _fQg[bound] || fq >= _fQg[std::min(tau + tQ + 1, _nQ)])
      return 0;

    // Find the the neighbours of "Q" on the grid
    int j;
    for (j = tau + tQ - bound; j >= 0; j--)
      if (fq < _fQg[tau+tQ-j+1])
        break;

    // Compute the interpolant
    double w_int = 1;
    for (int delta = tau - j; delta <= tau - j + _InterDegree; delta++)
      if (delta != tau)
        w_int *= ( fq - _fQg[delta] ) / ( _fQg[tau] - _fQg[delta] );

    return w_int;
  }

  //_________________________________________________________________________________
  template<class T>
  double QGrid<T>::DerInterpolant(int const& tQ, int const& tau, double const& Q) const
  {
    // Define the lower bound of the interpolation range
    const int bound = std::max(tau + tQ - _InterDegree, 0);

    // Return zero if Q is outside the allowed range
    if (Q < _Qg[bound] || Q >= _Qg[std::min(tau + tQ + 1, _nQ)])
      return 0;

    // Find the the neighbours of "Q" on the grid
    int j;
    for (j = tau + tQ - bound; j >= 0; j--)
      if (Q < _Qg[tau+tQ-j+1])
        break;

    // Compute the interpolant
    double dw_int = 0;
    for (int gamma = tau - j; gamma <= tau - j + _InterDegree; gamma++)
      {
        double w = 1;
        for (int delta = tau - j; delta <= tau - j + _InterDegree; delta++)
          if (delta != tau && delta != gamma)
            w *= ( Q - _Qg[delta] ) / ( _Qg[tau] - _Qg[delta] );
        if (gamma != tau)
          {
            w /= _Qg[tau] - _Qg[gamma];
            dw_int += w;
          }
      }
    return dw_int;
  }

  //_________________________________________________________________________________
  template<class T>
  double QGrid<T>::IntInterpolant(int const& tQ, int const& tau, double const& Qa, double const& Qb) const
  {
    // Return 0 if "Qa" and "Qb" are outside the range in which the
    // interpolant is different from zero.
    if (Qa > _Qg[tau + tQ + 1] || Qb < _Qg[std::max(tau + tQ - _InterDegree, 0)])
      return 0;

    // Construct interpolant
    double iw_int = 0;
    for (int i = 0; i <= std::min(_InterDegree, tau); i++)
      {
        if (_Qg[tau + tQ - i] > Qb || _Qg[tau + tQ - i + 1] < Qa)
          continue;

        // Product of denominators
        double dp = 1;
        std::vector<double> r(_InterDegree);
        int j = 0;
        for (int m = 0; m <= _InterDegree; m++)
          if(m != i)
            {
              dp /= _Qg[tau] - _Qg[tau-i+m];
              r[j++] = _Qg[tau-i+m];
            }

        // Expansion coefficients
        const std::vector<double> p = ProductExpansion(r);

        // Integration bounds
        const double Qab = std::max(Qa, _Qg[tau + tQ - i]);
        const double Qbb = std::min(Qb, _Qg[tau + tQ - i + 1]);

        // Sum of the integrals
        double sum = 0;
        for (int n = 0; n <= _InterDegree; n++)
          sum += pow(-1, n) * p[n] * ( pow(Qbb, _InterDegree - n + 1) - pow(Qab, _InterDegree - n + 1) ) / ( _InterDegree - n + 1 );

        iw_int += dp * sum;
      }
    return iw_int;
  }

  //_________________________________________________________________________________
  template<class T>
  std::tuple<int, int, int> QGrid<T>::SumBounds(double const& Q) const
  {
    // Initialise output tuple
    std::tuple<int, int, int> bounds{0, 0, 0};

    // Return if "Q" is outside the grid range (no sum will be
    // performed).
    if (Q < _Qg.front() - eps12 || Q > _Qg.back() + eps12)
      return bounds;

    // If Q falls in the tiny gap at the thresholds assume the node
    // below.
    for (int iQ = 1; iQ < (int) _nQg.size() - 1; iQ++)
      if (Q > _Qg[_nQg[iQ] - 1] && Q <= _Qg[_nQg[iQ]])
        {
          std::get<1>(bounds) = _nQg[iQ] - 1;
          std::get<2>(bounds) = _nQg[iQ];
          return bounds;
        }

    // Identify the subgrid in which Q falls
    int iQ;
    for (iQ = 0; iQ < (int) _nQg.size() - 1; iQ++)
      if (Q > _Qg[_nQg[iQ]] && Q <= _Qg[_nQg[iQ + 1]])
        break;

    if (iQ == (int) _nQg.size() - 1)
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

    // Determine the actual bounds
    const int low       = std::lower_bound(_Qg.begin() + 1, _Qg.end(), Q) - _Qg.begin();
    std::get<1>(bounds) = low;
    std::get<2>(bounds) = low;

    if (std::abs(Q / _Qg[low] - 1) <= eps12)
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
    // Get summation bounds
    const std::tuple<int, int, int> bounds = SumBounds(Q);
    const double                    fq     = _TabFunc(Q);

    // First create a copy of template object with the first
    // component...
    int tau  = std::get<1>(bounds);
    T result = Interpolant(std::get<0>(bounds), tau, fq) * _GridValues[tau];

    // ...then loop and add the extra terms
    for (tau = tau + 1; tau < std::get<2>(bounds); tau++)
      result += Interpolant(std::get<0>(bounds), tau, fq) * _GridValues[tau];

    return result;
  }

  //_________________________________________________________________________________
  template<class T>
  T QGrid<T>::Derive(double const& Q) const
  {
    // Get summation bounds
    const std::tuple<int, int, int> bounds = SumBounds(Q);

    // First create a copy of template object with the first
    // component...
    int tau  = std::get<1>(bounds);
    T result = DerInterpolant(std::get<0>(bounds), tau, Q) * _GridValues[tau];

    // ...then loop and add the extra terms
    for (tau = tau + 1; tau < std::get<2>(bounds); tau++)
      result += DerInterpolant(std::get<0>(bounds), tau, Q) * _GridValues[tau];

    return result;
  }

  //_________________________________________________________________________________
  template<class T>
  T QGrid<T>::Integrate(double const& Qa, double const& Qb) const
  {
    // Order integration bounds and adjust sign if necessary
    double Qao = std::min(Qa, Qb);
    double Qbo = std::max(Qa, Qb);
    int    sgn = (Qb > Qa ? 1 : -1);

    // Sum bounds and control parameters at the integral boundaries
    const std::tuple<int, int, int> boundsa = SumBounds(Qao);
    const std::tuple<int, int, int> boundsb = SumBounds(Qbo);

    // Initialise result to zero
    T result = T{0*_GridValues[0]};

    // First term
    for (int tau = std::get<1>(boundsa); tau < std::get<2>(boundsa); tau++)
      result += IntInterpolant(std::get<0>(boundsa), tau, Qao, _Qg[std::get<1>(boundsa) + std::get<0>(boundsa) + 1]) * _GridValues[tau];

    // Second term
    for (int gamma = std::get<1>(boundsa) + std::get<0>(boundsa) + 1; gamma <= std::get<1>(boundsb) + std::get<0>(boundsb); gamma++)
      {
        // Skip the tiny interval at the thresholds
        if (std::abs(_Qg[gamma+1] - _Qg[gamma]) < eps8)
          continue;
        const std::tuple<int, int, int> bounds = SumBounds(_Qg[gamma] * ( 1 + eps8 ));
        for (int tau = std::get<1>(bounds); tau < std::get<2>(bounds); tau++)
          result += IntInterpolant(std::get<0>(bounds), tau, _Qg[gamma], _Qg[gamma+1]) * _GridValues[tau];
      }

    // Third term
    for (int tau = std::get<1>(boundsb); tau < std::get<2>(boundsb); tau++)
      result -= IntInterpolant(std::get<0>(boundsb), tau, Qbo, _Qg[std::get<1>(boundsb) + std::get<0>(boundsb) + 1]) * _GridValues[tau];

    return sgn * result;
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
  template class QGrid<matrix<double>>;
  template class QGrid<Distribution>;
  template class QGrid<Set<Distribution>>;
  template class QGrid<DoubleObject<Distribution>>;
  template class QGrid<Operator>;
  template class QGrid<Set<Operator>>;
  template class QGrid<DoubleObject<Operator>>;
  template class QGrid<DoubleObject<Distribution, Operator>>;
  template class QGrid<DoubleObject<Operator, Distribution>>;
  template class QGrid<Set<DoubleObject<Distribution, Operator>> >;
  template class QGrid<Set<DoubleObject<Operator, Distribution>> >;
  template class QGrid<DoubleDistribution>;
  template class QGrid<DoubleOperator>;
}
