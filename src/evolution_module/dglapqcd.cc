//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/dglapqcd.h"
#include "apfel/ode.h"
#include "apfel/subgrid.h"
#include "apfel/grid.h"

using namespace std;

namespace apfel {

  //_________________________________________________________________________________
  DglapQCD::DglapQCD(Grid                         const& g,
		     function<double(int,double)> const& InPDFs,
		     double                       const& MuDistRef,
		     double                       const& AlphaRef,
		     double                       const& MuAlphaRef,
		     vector<double>               const& Masses,
		     vector<double>               const& Thresholds,
		     int                          const& PertOrder,
		     int                          const& nstep,
		     double                       const& xi):
    MatchedEvolution(Set<Distribution>{EvolutionBasis{}, {}}, MuDistRef, Masses, Thresholds),
    _InPDFs(InPDFs),
    _PertOrder(PertOrder),
    _nstep(nstep),
    _as(AlphaRef, MuAlphaRef, Masses, Thresholds, PertOrder, nstep)
  {
    // Allocate splitting functions operators
    unordered_map<int,unordered_map<int,Operator>> OpMap;
    const Operator O0ns{g, P0ns{}};
    const Operator O0qg{g, P0qg{}};
    const Operator O0gq{g, P0gq{}};
    for (int nf = 3; nf <= 6; nf++)
      {
	const Operator O0gg{g, P0gg{nf}};
	const Operator O0qgnf = nf * O0qg;
	unordered_map<int,Operator> OM;
	OM.insert({EvolutionBasis::PNSP,O0ns});
	OM.insert({EvolutionBasis::PNSM,O0ns});
	OM.insert({EvolutionBasis::PNSV,O0ns});
	OM.insert({EvolutionBasis::PQQ, O0ns});
	OM.insert({EvolutionBasis::PQG, O0qgnf});
	OM.insert({EvolutionBasis::PGQ, O0gq});
	OM.insert({EvolutionBasis::PGG, O0gg});
	OpMap.insert({nf,OM});
      }

    // Allocate initial scale distributions
    unordered_map<int,Distribution> DistMap;
    for (int i = EvolutionBasis::GLUON; i <= EvolutionBasis::V35; i++)
      DistMap.insert({i,PDF{g, InPDFs, i}});

    // Allocate maps
    unordered_map<int,EvolutionBasis> Bases;
    for (int nf = 3; nf <= 6; nf++)
      Bases.insert({nf,EvolutionBasis{nf}});

    // Allocate set of operators
    for (int nf = 3; nf <= 6; nf++)
      _SplittingFunctions.insert({nf,Set<Operator>{Bases.at(nf), OpMap.at(nf)}});

    // Find number of active flavour at the initial scale
    const auto nfi = lower_bound(Thresholds.begin()+1, Thresholds.end(), MuDistRef) - Thresholds.begin();

    // Allocate set of initial distributions.
    //SetObjectRef(Set<Distribution>{Bases.at(nfi), DistMap});
    _ObjRef = Set<Distribution>{Bases.at(nfi), DistMap};

    ///////////////////////////
    _SplittingFunctions.at(nfi) * _ObjRef;
    ///////////////////////////
  }

  //_________________________________________________________________________________
  Set<Distribution> DglapQCD::EvolveObject(int const& nf, Set<Distribution> const& sd0, double const& mu02, double const& mu2) const
  {
    // Return immediately "sd0" if "mu02" and "mu2" are equal
    if (mu02 == mu2)
       return sd0;

    // Numerical solution of the evolution equation with fourth-order Runge-Kutta.
    // Use "_nstep" steps for the evolution.
    auto       lrrat = log(mu2/mu02);
    const auto dlr   = lrrat / _nstep;
    auto       f     = sd0;

    const auto df = rk4setd([&](double const& mu, Set<Distribution> const& f)->Set<Distribution>{ return Derivative(nf, _as.GetObject(mu), f); });

    for (auto k = 0; k < _nstep; k++)
      {
	f += df(lrrat, f, dlr);
	lrrat += dlr;
      }

    return f;
  }

  //_________________________________________________________________________________
  Set<Distribution> DglapQCD::MatchObject(bool const& Up, Set<Distribution> const& sd, double const& LogKth) const
  {
    /*
    if(_PertOrder <= 0)
      return sd;

    const auto sgn = ( Up ? 1 : -1);
    //const auto ep = Coup / FourPi;
    //const double c[] = { 1, sgn * 2. / 3. * LogKth, 4. / 9. * pow(LogKth,2) + sgn *  38. / 3. * LogKth + sgn * 14. / 3. };
    //double match = 0;
    //for (auto i = 0; i <= _PertOrder; i++) match += c[i] * pow(ep,i);
    */
    return sd;
  }

  //_________________________________________________________________________________
  Set<Distribution> DglapQCD::Derivative(int const& nf, double const& as, Set<Distribution> const& f) const
  {
    // LO
    auto df = as * _SplittingFunctions.at(nf) * f;

    /*
    // NLO
    if (_PertOrder > 0)
      df += as * as * f;

    // NNLO
    if (_PertOrder > 0)
      df += as * as * as * f;
    */
    return df;
  }

  //_________________________________________________________________________________
  DglapQCD::P0ns::P0ns():
    Expression()
  {
  }
  double DglapQCD::P0ns::Regular(double const& x) const
  {
    return - 2 * CF * ( 1 + x );
  }
  double DglapQCD::P0ns::Singular(double const& x) const
  {
    return 4 * CF / ( 1 - x );
  }
  double DglapQCD::P0ns::Local(double const& x) const
  {
    return 4 * CF * log( 1 - x ) + 3 * CF;
  }

  //_________________________________________________________________________________
  DglapQCD::P0qg::P0qg():
    Expression()
  {
  }
  double DglapQCD::P0qg::Regular(double const& x) const
  {
    return 2 * ( 1 - 2 * x + 2 * x * x );
  }

  //_________________________________________________________________________________
  DglapQCD::P0gq::P0gq():
    Expression()
  {
  }
  double DglapQCD::P0gq::Regular(double const& x) const
  {
    return 4 * CF * ( - 1 + 0.5 * x + 1 / x );
  }

  //_________________________________________________________________________________
  DglapQCD::P0gg::P0gg(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double DglapQCD::P0gg::Regular(double const& x) const
  {
    return 4 * CA * ( - 2 + x - x * x + 1 / x );
  }
  double DglapQCD::P0gg::Singular(double const& x) const
  {
    return 4 * CA / ( 1 - x );
  }
  double DglapQCD::P0gg::Local(double const& x) const
  {
    return 4 * CA * log( 1 - x ) - 2 / 3. * _nf + 11 / 3. * CA;
  }

  //_________________________________________________________________________________
  DglapQCD::EvolutionBasis::EvolutionBasis():
    ConvolutionMap("EvolutionBasis")
  {
  }

  //_________________________________________________________________________________
  DglapQCD::EvolutionBasis::EvolutionBasis(int const& nf):
    ConvolutionMap("EvolutionBasis")
  {
    // dg = Pgg * g + Pgq * Sigma
    _rules[GLUON] = { {PGG, GLUON, +1}, {PGQ, SIGMA, +1} };
    // dSigma = Pqg * g + ( Pnsp + Pps ) * Sigma
    _rules[SIGMA] = { {PQG, GLUON, +1}, {PQQ, SIGMA, +1} };
    // dV = Pnsv * V
    _rules[VALENCE] = { {PNSV, VALENCE, +1} };
    // d{T,V}3 = Pnsp * {T,V}3
    if (nf > 1)
      {
	_rules[T3] = { {PNSP, T3, +1} };
	_rules[V3] = { {PNSM, V3, +1} };
      }
    else
      {
	_rules[T3] = _rules[SIGMA];
	_rules[V3] = _rules[VALENCE];
      }
    // d{T,V}8 = Pnsp * {T,V}8
    if (nf > 2)
      {
	_rules[T8] = { {PNSP, T8, +1} };
	_rules[V8] = { {PNSM, V8, +1} };
      }
    else
      {
	_rules[T8] = _rules[SIGMA];
	_rules[V8] = _rules[VALENCE];
      }
    // d{T,V}15 = Pnsp * {T,V}15
    if (nf > 3)
      {
	_rules[T15] = { {PNSP, T15, +1} };
	_rules[V15] = { {PNSM, V15, +1} };
      }
    else
      {
	_rules[T15] = _rules[SIGMA];
	_rules[V15] = _rules[VALENCE];
      }
    // d{T,V}24 = Pnsp * {T,V}24
    if (nf > 4)
      {
	_rules[T24] = { {PNSP, T24, +1} };
	_rules[V24] = { {PNSM, V24, +1} };
      }
    else
      {
	_rules[T24] = _rules[SIGMA];
	_rules[V24] = _rules[VALENCE];
      }
    // d{T,V}35 = Pnsp * {T,V}35
    if (nf > 5)
      {
	_rules[T35] = { {PNSP, T35, +1} };
	_rules[V35] = { {PNSM, V35, +1} };
      }
    else
      {
	_rules[T35] = _rules[SIGMA];
	_rules[V35] = _rules[VALENCE];
      }
  };

  //_________________________________________________________________________________
  DglapQCD::PDF::PDF(Grid const& gr, function<double(int,double)> const& InPDFs, int const& ipdf):
    Distribution(gr)
  {
    for (auto const& ix: _grid.GetJointGrid().GetGrid())
      if (ix < 1) _distributionJointGrid.push_back(InPDFs(ipdf,ix));
      else        _distributionJointGrid.push_back(0);

    for (auto ig=0; ig<_grid.nGrids(); ig++)
      {
        vector<double> sg;
        for (auto const& ix: _grid.GetSubGrid(ig).GetGrid())
          if (ix < 1) sg.push_back(InPDFs(ipdf,ix));
          else        sg.push_back(0);
        _distributionSubGrid.push_back(sg);
      }
  };

}

