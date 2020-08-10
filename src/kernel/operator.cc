//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/operator.h"
#include "apfel/lagrangeinterpolator.h"
#include "apfel/integrator.h"
#include "apfel/constants.h"
#include "apfel/messages.h"

namespace apfel
{
  //_________________________________________________________________________
  Operator::Operator(Grid const& gr):
    _grid(gr)
  {
  }

  //_________________________________________________________________________
  Operator::Operator(Grid const& gr, Expression const& expr, double const& eps):
    _grid(gr)
  {
    // Interpolator object for the interpolating functions
    const LagrangeInterpolator li{gr};

    // Scaling factor
    const double eta = expr.eta();

    // Number of grids.
    const int ng = _grid.nGrids();

    // Loop over the subgrids
    _Operator.resize(ng);
    for (int ig = 0; ig < ng; ig++)
      {
        // Define the global subgrid
        const SubGrid& sg = _grid.GetSubGrid(ig);

        // Get vector with the grid nodes
        const std::vector<double>& xg = sg.GetGrid();

        // Number of grid points
        const int nx = sg.nx();

        // Interpolation degree.
        const int id = sg.InterDegree();

        // Limit the loop over "beta" according to whether "sg" is
        // external.
        const int gbound = ( sg.IsExternal() ? nx : 0 );

        _Operator[ig].resize(gbound + 1, nx + 1, 0);
        for (int beta = 0; beta <= gbound; beta++)
          {
            const double xbeta = xg[beta];

            // Set xbeta as external variable in case needed for the
            // computation of the operator (this is in general not
            // needed but in the case of GPD evolution it is).
            expr.SetExternalVariable(xbeta);

            // Local function. Can be computed outside the 'alpha'
            // loop.
            const double L   = eta * expr.Local(xbeta / xg[beta+1] / eta);
            const double lxe = log(xbeta / eta);
            for (int alpha = beta; alpha <= nx; alpha++)
              {
                // Weight of the subtraction term (it is a delta if
                // eta = 1).
                const double ws = li.InterpolantLog(alpha, lxe, sg);

                // Given that the interpolation functions have
                // discontinuos derivative on the nodes and are
                // wiggly, it turns out that it is convenient to split
                // the integrals into (id+1) intervals on each of
                // which the integrand is smooth. This way, even
                // though more integrals have to be computed, the
                // integration converges faster and is more accurate.

                // Number of grid intervals we need to integrate over
                const int nmin = fmax(0, alpha + 1 - nx);
                const int nmax = fmin(id, alpha - beta) + 1;

                // Integral
                double I = 0;
                for (int jint = nmin; jint < nmax; jint++)
                  {
                    // Define integration bounds of the first
                    // iteration.
                    const double c = xbeta / xg[alpha-jint+1];
                    const double d = xbeta / xg[alpha-jint];

                    // Define integrand and the corresponding
                    // "Integrator" object.
                    const std::function<double(double const&)> integrand = [&] (double const& x) -> double
                    {
                      const double z = x / eta;
                      if (z >= 1)
                        return 0;
                      const double wr = li.InterpolantLog(alpha, log(xbeta / x), _grid.GetSubGrid(ig));
                      return expr.Regular(z) * wr + expr.Singular(z) * ( wr - ws );
                    };
                    const Integrator Io{integrand};

                    // Compute the integral
                    I += Io.integrate(c, d, eps);
                  }
                // Add the local part
                _Operator[ig](beta, alpha) = I + L * ws;
              }
          }
      }
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
    // Fast method to check that we are using the same Grid
    if (&this->_grid != &d.GetGrid())
      throw std::runtime_error(error("Operator::operator*=", "Operator and Distribution grids do not match"));

    const std::vector<std::vector<double>>& sg = d.GetDistributionSubGrid();
    std::vector<std::vector<double>> s(sg);
    std::vector<double> j;

    // Compute the the distribution on the subgrids
    int const ng = _grid.nGrids();
    for (int ig = 0; ig < ng; ig++)
      {
        int const nx = _grid.GetSubGrid(ig).nx();

        // If the grid is external the product between the operator
        // and the distribution has to be done in a standard way.
        if (this->_grid.GetSubGrid(ig).IsExternal())
          for (int alpha = 0; alpha <= nx; alpha++)
            {
              s[ig][alpha] = 0;
              for (int beta = alpha; beta <= nx; beta++)
                s[ig][alpha] += _Operator[ig](alpha, beta) * sg[ig][beta];
            }
        // If the grid is internal the product between the operator
        // and the distribution has to be done exploiting the symmetry
        // of the operator.
        else
          for (int alpha = 0; alpha <= nx; alpha++)
            {
              s[ig][alpha] = 0;
              for (int beta = alpha; beta <= nx; beta++)
                s[ig][alpha] += _Operator[ig](0,beta-alpha) * sg[ig][beta];
            }

        // Set to zero the values above one
        for (int alpha = nx + 1; alpha < this->_grid.GetSubGrid(ig).InterDegree() + nx + 1; alpha++)
          s[ig][alpha] = 0;

        // Compute the distribution on the joint grid
        double xtrans;
        if (ig < ng-1)
          xtrans = this->_grid.GetSubGrid(ig+1).xMin();
        else
          xtrans = 1 + 2 * eps12;

        for (int alpha = 0; alpha <= nx; alpha++)
          {
            const double x = this->_grid.GetSubGrid(ig).GetGrid()[alpha];
            if (xtrans - x < eps12)
              break;
            j.push_back(s[ig][alpha]);
          }
      }

    // Set to zero the values above one
    for (int alpha = 0; alpha < this->_grid.GetJointGrid().InterDegree(); alpha++)
      j.push_back(0);

    return Distribution{d, s, j};
  }

  //_________________________________________________________________________
  Operator& Operator::operator *= (Operator const& o)
  {
    // Fast method to check that we are using the same Grid
    if (&this->_grid != &o.GetGrid())
      throw std::runtime_error(error("Operator::operator*=", "Operators grid does not match"));

    const std::vector<matrix<double>> v = _Operator;

    const int ng = _grid.nGrids();
    for (int ig = 0; ig < ng; ig++)
      {
        const int nx = this->_grid.GetSubGrid(ig).nx();

        // If the grid is external the product between the operators
        // has to be done in a standard way.
        if (this->_grid.GetSubGrid(ig).IsExternal())
          {
            for (int alpha = 0; alpha <= nx; alpha++)
              for (int beta = alpha; beta <= nx; beta++)
                {
                  _Operator[ig](alpha, beta) = 0;
                  for (int gamma = alpha; gamma <= beta; gamma++)
                    _Operator[ig](alpha, beta) += v[ig](alpha, gamma) * o._Operator[ig](gamma, beta);
                }
          }
        // If the grid is internal the product between the operators
        // can be done exploiting the symmetry of the operators.
        else
          {
            for (int beta = 0; beta <= nx; beta++)
              {
                _Operator[ig](0, beta) = 0;
                for (int gamma = 0; gamma <= beta; gamma++)
                  _Operator[ig](0, beta) += v[ig](0, gamma) * o._Operator[ig](0, beta - gamma);
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
          _Operator[ig](alpha, beta) *= s;

    return *this;
  }

  //_________________________________________________________________________
  Operator& Operator::operator *= (std::function<double(double const&)> f)
  {
    if (!_grid.ExtGrids())
      throw std::runtime_error(error("Operator::operator*=", "Multiplication by a function not allowed on internal grids"));

    for (size_t ig = 0; ig < _Operator.size(); ig++)
      {
        // Get ig-th subgrid
        const std::vector<double>& sg = _grid.GetSubGrid(ig).GetGrid();
        for (size_t alpha = 0; alpha < _Operator[ig].size(0); alpha++)
          for (size_t beta = alpha; beta < _Operator[ig].size(1); beta++)
            _Operator[ig](alpha, beta) *= f(sg[alpha]);
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
          _Operator[ig](alpha, beta) *= r;

    return *this;
  }

  //_________________________________________________________________________
  Operator& Operator::operator += (Operator const& o)
  {
    // Fast method to check that we are using the same Grid
    if (&this->_grid != &o.GetGrid())
      throw std::runtime_error(error("Operator::operator+=", "Operators grid does not match"));

    for (size_t ig = 0; ig < _Operator.size(); ig++)
      for (size_t alpha = 0; alpha < _Operator[ig].size(0); alpha++)
        for (size_t beta = alpha; beta < _Operator[ig].size(1); beta++)
          _Operator[ig](alpha, beta) += o._Operator[ig](alpha, beta);

    return *this;
  }

  //_________________________________________________________________________
  Operator& Operator::operator -= (Operator const& o)
  {
    // Fast method to check that we are using the same Grid
    if (&this->_grid != &o.GetGrid())
      throw std::runtime_error(error("Operator::operator+=", "Operators grid does not match"));

    for (size_t ig = 0; ig < _Operator.size(); ig++)
      for (size_t alpha = 0; alpha < _Operator[ig].size(0); alpha++)
        for (size_t beta = alpha; beta < _Operator[ig].size(1); beta++)
          _Operator[ig](alpha, beta) -= o._Operator[ig](alpha, beta);

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
  Operator operator * (std::function<double(double const&)> f, Operator rhs)
  {
    return rhs *= f;
  }

  //_________________________________________________________________________
  Operator operator * (Operator lhs, std::function<double(double const&)> f)
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

  //_________________________________________________________________________________
  std::ostream& operator << (std::ostream& os, Operator const& op)
  {
    const std::vector<matrix<double>> om = op.GetOperator();
    os << "Operator: " << &op << "\n";
    os << "Operator on the first SubGrid:" << "\n";
    const std::ostringstream default_format;
    os << std::scientific;
    os.precision(1);
    for (int i = 0; i < (int) om[0].size(0); i++)
      {
        os << "{i: " << i << ", [";
        for (int j = 0; j < (int) om[0].size(1); j++)
          os << om[0](i, j) << " ";
        os << "\b]";
        if (i != (int) om[0].size(0) - 1)
          os << "}\n";
      }
    os.copyfmt(default_format);
    return os;
  }
}
