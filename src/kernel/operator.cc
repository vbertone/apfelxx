//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/operator.h"
#include "apfel/tools.h"

#include <iostream>
#include <cmath>

namespace apfel
{
  //_________________________________________________________________________
  Operator::Operator(Grid const& gr, Expression const& expr, accuracy const& eps):
    Integrator{},
    LagrangeInterpolator{gr},
    _grid(gr),
    _expr(&expr)
  {
    // Scaling factor
    _eta = _expr->eta();

    // Number of grids.
    const int ng = _grid.nGrids();

    // Loop over the subgrids.
    _Operator.resize(ng);
    for (_ig = 0; _ig < ng; _ig++)
      {
	// Define the global subgrid.
	const SubGrid& sg = _grid.GetSubGrid(_ig);

	// Get vector with the grid nodes.
	const vector<double>& xg = sg.GetGrid();

	// Number of grid points.
	const int nx = sg.nx();

	// Interpolation degree.
	const int id = sg.InterDegree();

	// Limit the loop over "beta" according to whether "sg" is
	// external.
	const int gbound = ( sg.IsExternal() ? nx : 0 );

	_Operator[_ig].resize(gbound+1, nx+1, 0);
	for (auto beta = 0; beta <= gbound; beta++)
	  {
	    _xbeta = xg[beta];
	    // Local function. Can be computed outside the 'alpha'
	    // loop.
	    const double L = _eta * _expr->Local(_xbeta / xg[beta+1] / _eta);
	    for (_alpha = beta; _alpha <= nx; _alpha++)
	      {
		// Weight of the subtraction term (it is a delta if
		// _eta = 1).
		_ws = Interpolant(_alpha, log(_xbeta / _eta), sg);

		// Given that the interpolation functions have
		// discontinuos derivative on the nodes and are
		// wiggly, it turns out that it is convenient to split
		// the integrals into (id+1) intervals on each of
		// which the integrand is smooth. This way, even
		// though more integrals have to be computed, the
		// integration converges faster.

		// Number of grid intervals we need to integrate over.
		const int nmax = fmin(id,_alpha-beta) + 1;
		const int nmin = fmax(0,_alpha+1-nx);

		// Integral.
		double I = 0;
		for (auto jint = nmin; jint < nmax; jint++)
		  {
		    // Define integration bounds of the first
		    // iteration.
		    const double c = _xbeta / xg[_alpha-jint+1];
		    const double d = _xbeta / xg[_alpha-jint];

		    // Compute the integral.
		    I += integrate(c, d, eps);
		  }
		// Add the local part.
		_Operator[_ig](beta,_alpha) = I + L * _ws;
	      }
	  }
      }
  }

  //_________________________________________________________________________
  double Operator::integrand(double const& x) const
  {
    const double z = x / _eta;
    if (z >= 1)
      return 0;
    const double wr = Interpolant(_alpha, log(_xbeta / x), _grid.GetSubGrid(_ig));
    return _expr->Regular(z) * wr + _expr->Singular(z) * ( wr - _ws );
  }

  //_________________________________________________________________________
  Operator& Operator::operator = (Operator const& o)
  {
    if (this != &o)
      *this = o;
    return *this;
  }

  //_________________________________________________________________________
  Distribution Operator::operator *= (Distribution const& d) const
  {
    // Fast method to check that we are using the same Grid.
    if (&this->_grid != &d.GetGrid())
      throw runtime_exception("Operator::operator*=", "Operator and Distribution grids do not match");

    const vector<vector<double>>& sg = d.GetDistributionSubGrid();
    vector<vector<double>> s(sg);
    vector<double> j;

    // Compute the the distribution on the subgrids.
    int const ng = _grid.nGrids(); //sg.size();
    for (auto ig = 0; ig < ng; ig++)
      {
        int const nx = _grid.GetSubGrid(ig).nx();

	// If the grid is external the product between the operator
	// and the distribution has to be done in a standard way.
	if (this->_grid.GetSubGrid(ig).IsExternal())
	  {
	    for (auto alpha = 0; alpha <= nx; alpha++)
	      {
		s[ig][alpha] = 0;
		for (auto beta = alpha; beta <= nx; beta++)
		  s[ig][alpha] += _Operator[ig](alpha,beta) * sg[ig][beta];
	      }
	  }
	// If the grid is internal the product between the operator
	// and the distribution has to be done exploiting the symmetry
	// of the operator.
	else
	  {
	    for (auto alpha = 0; alpha <= nx; alpha++)
	      {
		s[ig][alpha] = 0;
		for (auto beta = alpha; beta <= nx; beta++)
		  s[ig][alpha] += _Operator[ig](0,beta-alpha) * sg[ig][beta];
	      }
	  }

	// Set to zero the values above one.
	for (auto alpha = nx + 1; alpha < this->_grid.GetSubGrid(ig).InterDegree() + nx + 1; alpha++)
	  s[ig][alpha] = 0;

	// Compute the the distribution on the joint grid.
        double xtrans;
        if (ig < ng-1)
	  xtrans = this->_grid.GetSubGrid(ig+1).xMin();
        else
          xtrans = 1 + 2 * eps12;

        for (auto alpha = 0; alpha <= nx; alpha++)
          {
            const double x = this->_grid.GetSubGrid(ig).GetGrid()[alpha];
            if (xtrans - x < eps12)
	      break;
            j.push_back(s[ig][alpha]);
          }
      }

    // Set to zero the values above one.
    for (auto alpha = 0; alpha < this->_grid.GetJointGrid().InterDegree(); alpha++)
      j.push_back(0);

    return Distribution{d, s, j};
  }

  //_________________________________________________________________________
  Operator& Operator::operator *= (Operator const& o)
  {
    // fast method to check that we are using the same Grid.
    if (&this->_grid != &o.GetGrid())
      throw runtime_exception("Operator::operator*=", "Operators grid does not match");

    const auto v = _Operator;

    const int ng = _grid.nGrids(); //sg.size();
    for (auto ig = 0; ig < ng; ig++)
      {
        const int nx = this->_grid.GetSubGrid(ig).nx();

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
  Operator& Operator::operator *= (double const& s)
  {
    for (size_t ig = 0; ig < _Operator.size(); ig++)
      for (size_t alpha = 0; alpha < _Operator[ig].size(0); alpha++)
        for (size_t beta = alpha; beta < _Operator[ig].size(1); beta++)
          _Operator[ig](alpha,beta) *= s;

    return *this;
  }


  //_________________________________________________________________________
  Operator& Operator::operator *= (function<double(double const&)> f)
  {
    if (!_grid.ExtGrids())
      throw runtime_exception("Operator::operator*=", "Multiplication by a function not allowed on internal grids");

    for (size_t ig = 0; ig < _Operator.size(); ig++)
      {
	// Get ig-th subgrid
	const auto& sg = _grid.GetSubGrid(ig).GetGrid();
	for (size_t alpha = 0; alpha < _Operator[ig].size(0); alpha++)
	  for (size_t beta = alpha; beta < _Operator[ig].size(1); beta++)
	    _Operator[ig](alpha,beta) *= f(sg[alpha]);
      }
    return *this;
  }

  //_________________________________________________________________________
  Operator& Operator::operator /= (double const& s)
  {
    const double r = 1 / s;
    for (size_t ig = 0; ig < _Operator.size(); ig++)
      for (size_t alpha = 0; alpha < _Operator[ig].size(0); alpha++)
        for (size_t beta = alpha; beta < _Operator[ig].size(1); beta++)
          _Operator[ig](alpha,beta) *= r;

    return *this;
  }

  //_________________________________________________________________________
  Operator& Operator::operator += (Operator const& o)
  {
    // fast method to check that we are using the same Grid.
    if (&this->_grid != &o.GetGrid())
      throw runtime_exception("Operator::operator+=", "Operators grid does not match");

    for (size_t ig = 0; ig < _Operator.size(); ig++)
      for (size_t alpha = 0; alpha < _Operator[ig].size(0); alpha++)
        for (size_t beta = alpha; beta < _Operator[ig].size(1); beta++)
          _Operator[ig](alpha,beta) += o._Operator[ig](alpha,beta);

    return *this;
  }

  //_________________________________________________________________________
  Operator& Operator::operator -= (Operator const& o)
  {
    // fast method to check that we are using the same Grid.
    if (&this->_grid != &o.GetGrid())
      throw runtime_exception("Operator::operator+=", "Operators grid does not match");

    for (size_t ig = 0; ig < _Operator.size(); ig++)
      for (size_t alpha = 0; alpha < _Operator[ig].size(0); alpha++)
        for (size_t beta = alpha; beta < _Operator[ig].size(1); beta++)
          _Operator[ig](alpha,beta) -= o._Operator[ig](alpha,beta);

    return *this;
  }

  //_________________________________________________________________________
  Distribution operator * (Operator lhs, Distribution const& rhs)
  {
    return lhs *= rhs;
  }

  //_________________________________________________________________________
  Operator operator * (Operator lhs, Operator const& rhs)
  {
    return lhs *= rhs;
  }

  //_________________________________________________________________________
  Operator operator * (double const& s, Operator rhs)
  {
    return rhs *= s;
  }

  //_________________________________________________________________________
  Operator operator * (Operator lhs, double const& s)
  {
    return lhs *= s;
  }

  //_________________________________________________________________________
  Operator operator * (function<double(double const&)> f, Operator rhs)
  {
    return rhs *= f;
  }

  //_________________________________________________________________________
  Operator operator * (Operator lhs, function<double(double const&)> f)
  {
    return lhs *= f;
  }

  //_________________________________________________________________________
  Operator operator / (Operator lhs, double const& s)
  {
    return lhs /= s;
  }

  //_________________________________________________________________________
  Operator operator + (Operator lhs, Operator const& rhs)
  {
    return lhs += rhs;
  }

  //_________________________________________________________________________
  Operator operator - (Operator lhs, Operator const& rhs)
  {
    return lhs -= rhs;
  }
}
