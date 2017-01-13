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
#include "apfel/distribution.h"
#include "apfel/tools.h"

namespace apfel
{
  //_________________________________________________________________________
  Operator::Operator(Grid const& gr, Expression const& expr, double const& eps):
    Integrator{},
    LagrangeInterpolator{gr},
    _grid(gr),
    _expr(&expr),
    _eps(eps)
  {
    // Number of grids
    const int ng = _grid.nGrids();

    // Loop over the subgrids    
    _Operator.resize(ng);
    for (_ig = 0; _ig < ng; _ig++)
      {
	// Define the global subgrid
	const auto& sg = _grid.GetSubGrid(_ig);

	// Get vector with the grid nodes
	const auto& xg = sg.GetGrid();

	// Number of grid points
	auto nx = sg.nx();

	// Interpolation degree
	const int id = sg.InterDegree();

	// Limit the loop over "_beta" according to whether "sg" is external
	const int gbound = ( sg.IsExternal() ? nx : 0 );

	_Operator[_ig].resize(gbound+1, nx+1, 0);
	for (_beta = 0; _beta <= gbound; _beta++)
	  {
	    const auto xbeta = xg[_beta];
	    for (_alpha = _beta; _alpha <= nx; _alpha++)
	      {
		// Weight of the subtraction term (independent of x)
		_ws = ( _alpha == _beta ? 1 : 0 );

		// Given that the interpolation functions have discontinuos derivative
		// on the nodes and are wiggly, it turns out that it is conveniente to
		// split the integrals into (id+1) intervals on each of which the integrand
		// is smooth. This way, even though more integrals have to be computed, the
		// integrator converges faster.

		// Number of grid intervals we need to integrate over.
		int nmax = fmin(id,_alpha-_beta) + 1;
		int nmin = fmax(0,_alpha+1-nx);

		// Integral
		double I = 0;
		for(auto jint = nmin; jint < nmax; jint++)
		  {
		    // Define integration bounds of the first iteration
		    const double c = xbeta / xg[_alpha-jint+1];
		    const double d = xbeta / xg[_alpha-jint];

		    // Compute the integral
		    I += integrate(c, d, _eps);
		  }
		_Operator[_ig](_beta,_alpha) = I;
	      }
	    // Add local part
	    _Operator[_ig](_beta,_beta) += _expr->Local(xbeta/xg[_beta+1]);
	  }
      }
  }

  //_________________________________________________________________________
  Operator::Operator(Operator const& obj):
    Integrator{},
    LagrangeInterpolator{obj._grid},
    _grid(obj._grid),
    _expr(nullptr),
    _eps(obj._eps),
    _Operator(obj._Operator)
  {
  }

  //_________________________________________________________________________
  double Operator::integrand(double const& x) const
  {
    double res = 0;
    try
      {
        const double wr = Interpolant(_alpha, log(_grid.GetSubGrid(_ig).GetGrid()[_beta] / x), _grid.GetSubGrid(_ig));
        res = _expr->Regular(x) * wr + _expr->Singular(x) * ( wr - _ws );
      }
    catch (std::exception &e)
      {
        throw runtime_exception("Operator::integrand",  "Operator Expression not defined.");
      }

    return res;
  }

  //_________________________________________________________________________
  Distribution Operator::operator*=(Distribution const& d) const
  {
    // Fast method to check that we are using the same Grid
    if (&this->_grid != &d.GetGrid())
      throw runtime_exception("Operator::operator*=", "Operator and Distribution grids do not match");

    // Compute the the distribution on the subgrids
    const auto& sg = d.GetDistributionSubGrid();
    vector<vector<double>> s(sg);

    int const ng = _grid.nGrids(); //sg.size();
    for (auto ig = 0; ig < ng; ig++)
      {
        int const nx = this->_grid.GetSubGrid(ig).nx();

	// If the grid is external the product between the operator and the distribution
	// has to be done in a standard way.
	if (this->_grid.GetSubGrid(ig).IsExternal())
	  {
	    for (auto alpha = 0; alpha <= nx; alpha++)
	      {
		s[ig][alpha] = 0;
		for (auto beta = alpha; beta <= nx; beta++)
		  s[ig][alpha] += _Operator[ig](alpha,beta) * sg[ig][beta];
	      }
	  }
	// If the grid is internal the product between the operator and the distribution
	// has to be done exploiting the symmetry of the operator.
	else
	  {
	    for (auto alpha = 0; alpha <= nx; alpha++)
	      {
		s[ig][alpha] = 0;
		for (auto beta = alpha; beta <= nx; beta++)
		  s[ig][alpha] += _Operator[ig](0,beta-alpha) * sg[ig][beta];
	      }
	  }

	// Set to zero the values above one
	for (auto alpha = nx + 1; alpha < this->_grid.GetSubGrid(ig).InterDegree() + nx + 1; alpha++)
	  s[ig][alpha] = 0;
      }

    // Compute the the distribution on the joint grid
    vector<double> j;
    for(auto ig = 0; ig < ng; ig++)
      {
        int const nx = this->_grid.GetSubGrid(ig).nx();

        double xtrans;
        if(ig < ng-1) xtrans = this->_grid.GetSubGrid(ig+1).xMin();
        else          xtrans = 1 + 2 * eps12;

        for(auto alpha = 0; alpha <= nx; alpha++)
          {
            double const x = this->_grid.GetSubGrid(ig).GetGrid()[alpha];
            if(xtrans - x < eps12) break;
            j.push_back(s[ig][alpha]);
          }
      }

    // Set to zero the values above one
    for (auto alpha = 0; alpha < this->_grid.GetJointGrid().InterDegree(); alpha++)
      j.push_back(0);

    return Distribution{d, s, j};
  }

  //_________________________________________________________________________
  Operator& Operator::operator*=(Operator const& o)
  {
    // fast method to check that we are using the same Grid
    if (&this->_grid != &o.GetGrid())
      throw runtime_exception("Operator::operator*=", "Operators grid does not match");

    auto v = _Operator;

    int const ng = _grid.nGrids(); //sg.size();
    for (auto ig = 0; ig < ng; ig++)
      {
        int const nx = this->_grid.GetSubGrid(ig).nx();

	// If the grid is external the product between the operators
	// has to be done in a standard way.
	if (this->_grid.GetSubGrid(ig).IsExternal())
	  {
	    for (auto alpha = 0; alpha <= nx; alpha++)
	      for (auto beta = alpha; beta <= nx; beta++)
		{
		  _Operator[ig](alpha,beta) = 0;
		  for (auto gamma = alpha; gamma <= beta; gamma++)
		    _Operator[ig](alpha,beta) += v[ig](alpha,gamma) * o._Operator[ig](gamma,beta);
		}
	  }
	// If the grid is internal the product between the operators
	// has to be done exploiting the symmetry of the operators.
	else
	  {
	    _Operator[ig].resize(nx+1, nx+1, 0);
	    for (auto alpha = 0; alpha <= nx; alpha++)
	      for (auto beta = alpha; beta <= nx; beta++)
		{
		  _Operator[ig](alpha,beta) = 0;
		  for (auto gamma = alpha; gamma <= beta; gamma++)
		    _Operator[ig](alpha,beta) += v[ig](0,gamma-alpha) * o._Operator[ig](0,beta-gamma);
		}
	  }
      }
    return *this;
  }

  //_________________________________________________________________________
  Operator& Operator::operator*=(double const& s)
  {
    for (size_t ig = 0; ig < _Operator.size(); ig++)
      for (size_t alpha = 0; alpha < _Operator[ig].size().first; alpha++)
        for (size_t beta = alpha; beta < _Operator[ig].size().second; beta++)
          _Operator[ig](alpha,beta) *= s;

    return *this;
  }

  //_________________________________________________________________________
  Operator& Operator::operator+=(Operator const& o)
  {
    // fast method to check that we are using the same Grid
    if (&this->_grid != &o.GetGrid())
      throw runtime_exception("Operator::operator+=", "Operators grid does not match");

    for (size_t ig = 0; ig < _Operator.size(); ig++)
      for (size_t alpha = 0; alpha < _Operator[ig].size().first; alpha++)
        for (size_t beta = alpha; beta < _Operator[ig].size().second; beta++)
          _Operator[ig](alpha,beta) += o._Operator[ig](alpha,beta);

    return *this;
  }

  //_________________________________________________________________________
  Operator& Operator::operator-=(Operator const& o)
  {
    // fast method to check that we are using the same Grid
    if (&this->_grid != &o.GetGrid())
      throw runtime_exception("Operator::operator+=", "Operators grid does not match");

    for (size_t ig = 0; ig < _Operator.size(); ig++)
      for (size_t alpha = 0; alpha < _Operator[ig].size().first; alpha++)
        for (size_t beta = alpha; beta < _Operator[ig].size().second; beta++)
          _Operator[ig](alpha,beta) -= o._Operator[ig](alpha,beta);

    return *this;
  }

  //_________________________________________________________________________
  Distribution operator*(Operator lhs, Distribution const& rhs)
  {
    return lhs *= rhs;
  }

  //_________________________________________________________________________
  Operator operator*(Operator lhs, Operator const& rhs)
  {
    return lhs *= rhs;
  }

  //_________________________________________________________________________
  Operator operator*(double const& s, Operator rhs)
  {
    return rhs *= s;
  }

  //_________________________________________________________________________
  Operator operator*(Operator lhs, double const& s)
  {
    return lhs *= s;
  }

  //_________________________________________________________________________
  Operator operator+(Operator lhs, Operator const& rhs)
  {
    return lhs += rhs;
  }

  //_________________________________________________________________________
  Operator operator-(Operator lhs, Operator const& rhs)
  {
    return lhs -= rhs;
  }
}
