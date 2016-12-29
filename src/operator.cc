//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <iostream>
#include <cmath>

#include "apfel/operator.h"
#include "apfel/grid.h"
#include "apfel/subgrid.h"
#include "apfel/expression.h"
#include "apfel/interpolator.h"

namespace apfel
{
  //_________________________________________________________________________
  Operator::Operator(Grid const& gr, Expression const& expr, Interpolator const& inter, double const& eps):
    _grid(gr),
    _expr(expr),
    _inter(inter),
    _eps(eps)
  {
    // Number of grids
    auto ng = _grid.nGrids();

    // Loop over the subgrids
    _Operator = new double**[ng];
    for (_ig = 0; _ig < ng; _ig++)
      {
	// Define the global subgrid
	auto sg = _grid.GetSubGrid(_ig);

	// Get vector with the grid nodes
	auto xg = sg.GetGrid();

	// Number of grid points
	auto nx = sg.nx();

	// Interpolation degree
	auto id = sg.InterDegree();

	// Limit the loop over "_beta" according to whether "sg" is external
	int gbound = ( sg.IsExternal() ? nx : 1 );

	_Operator[_ig] = new double*[gbound];
	for (_beta = 0; _beta < gbound; _beta++)
	  {
	    auto xbeta  = xg[_beta];
	    auto lxbeta = log(xbeta);
	    _Operator[_ig][_beta] = new double[nx];
	    for (_alpha = _beta; _alpha < nx; _alpha++)
	      {
		// Weight of the subtraction term (independent of x)
		_ws = _inter.Interpolant(_alpha, lxbeta, sg);

		/*
		double c = fmax(xbeta, xbeta / xg[_alpha+1]);
		double d = fmin(1, xbeta / xg[fmax(_alpha-id,_beta)]);

		// Compute the integral
		double I = this->integrate(c, d, _eps);
		*/

		// Given that the interpolation functions have discontinuos derivative
		// on the nodes and are wiggly, it turns out that it is conveniente to
		// split the integrals into (id+1) intervals on each of which the integrand
		// is smooth. This way, even though more integrals have to be computed, the
		// integrator converges faster.

		// Number of grid intervals we need to integrate over.
		auto nint = fmin(id,_alpha)+1;

		// Integral
		double I = 0;
		for(int jint=0; jint<nint; jint++)
		  {
		    // Define integration bounds of the first iteration
		    double c = xbeta / xg[_alpha-jint+1];
		    double d = xbeta / xg[_alpha-jint];

		    // Compute the integral
		    I += this->integrate(c, d, _eps);
		  }
		_Operator[_ig][_beta][_alpha] = I;
	      }
	    // Add local part
	    _Operator[_ig][_beta][_beta] += _expr.Local(xbeta/xg[_beta+1]);
	  }
      }
  }

  //_________________________________________________________________________
  double Operator::integrand(double const& x) const
  {
    double wr = _inter.Interpolant(_alpha, log(_grid.GetSubGrid(_ig).GetGrid()[_beta] / x), _grid.GetSubGrid(_ig));
    return _expr.Regular(x) * wr + _expr.Singular(x) * ( wr - _ws );
  }

}