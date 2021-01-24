//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/operator.h"
#include "apfel/integrator.h"
#include "apfel/messages.h"

#include <cmath>

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
    const LagrangeInterpolator li{_grid};

    // Scaling factor
    const double eta = expr.eta();

    // Number of grids.
    const int ng = _grid.nGrids();

    // Loop over the subgrids
    _Operator.resize(ng);
    for (int ig = 0; ig < ng; ig++)
      {
        // Get current subgrid
        const SubGrid& sg = _grid.GetSubGrid(ig);

        // Get vector with the grid nodes
        const std::vector<double>& xg = sg.GetGrid();

        // Get number of grid intervals in the definition range
        const int nx = sg.nx();

        // Get interpolation degree
        const int id = sg.InterDegree();

        // Get logarithmic step
        const double s = sg.Step();

        // Local function
        const double L = eta * expr.Local(1 / exp(s) / eta);

        // Initialise operator
        _Operator[ig].resize(1, nx);

        // Loop over the index beta. In fact beta = 0 because the size
        // of the first dimension of "_Operator" is one.
        for (int beta = 0; beta < (int) _Operator[ig].size(0); beta++)
          {
            // Log of the lower bound
            const double lxe = ( beta - nx ) * s - log(eta);

            // Loop over the index alpha
            for (int alpha = beta; alpha < (int) _Operator[ig].size(1); alpha++)
              {
                // Weight of the subtraction term (it is a delta if
                // eta = 1).
                const double ws = li.InterpolantLog(alpha, lxe, sg);

                // Given that the interpolation functions have
                // discontinuos derivative at the nodes, it turns out
                // that it is convenient to split the integrals into
                // id + 1 intervals on each of which the integrand is
                // smooth. Despite more integrals have to be computed,
                // the integration converges faster and is more
                // accurate.

                // Initialise integral
                double I = 0;

                // Run over the grid intervals over which the
                // interpolating function is different from zero.
                for (int j = 0; j < std::min(id, alpha - beta) + 1; j++)
                  {
                    // Define "Integrator" object.
                    const Integrator Ij{[&] (double const& y) -> double
                      {
                        const double z  = y / eta;
                        const double wr = li.InterpolantLog(alpha, log(xg[beta] / y), sg);
                        return expr.Regular(z) * wr + expr.Singular(z) * ( wr - ws );
                      }};
                    // Compute the integral
                    I += Ij.integrate(exp((beta - alpha + j - 1) * s), exp((beta - alpha + j) * s), eps);
                  }
                // Add the local part
                _Operator[ig](beta, alpha) = I + L * ws;
              }
          }
      }
  }

  //_________________________________________________________________________
  Distribution Operator::operator *= (Distribution const& d) const
  {
    // Fast method to check that we are using the same Grid
    if (&_grid != &d.GetGrid())
      throw std::runtime_error(error("Operator::operator *=", "Operator and Distribution grids do not match"));

    // Get number of subgrids
    const int ng = _grid.nGrids();

    // Get map of indices map from joint to subgrids
    const std::vector<std::vector<int>>& jsmap = _grid.JointToSubMap();

    // Get map of indices map from sub to joint to grids
    const std::vector<std::pair<int, int>>& sjmap = _grid.SubToJointMap();

    // Get joint distribution
    const std::vector<double>& dj = d.GetDistributionJointGrid();

    // Initialise output vectors
    std::vector<double> j(d.GetDistributionJointGrid().size(), 0);
    std::vector<std::vector<double>> s(ng);

    // Construct joint distribution first. The product between the
    // operator and the distribution is done exploiting the symmetry
    // of the operator.
    for (int alpha = 0; alpha < _grid.GetJointGrid().nx(); alpha++)
      {
        const std::pair<int, int> m = sjmap[alpha];
        for (int beta = m.second; beta < _grid.GetSubGrid(m.first).nx(); beta++)
          j[alpha] += _Operator[m.first](0, beta - m.second) * dj[jsmap[m.first][beta]];
      }

    // Compute the the distribution on the subgrids
    for (int ig = 0; ig < ng; ig++)
      {
        // Resize output on the "ig"-th subgrid
        s[ig].resize(d.GetDistributionSubGrid()[ig].size(), 0);
        for (int alpha = 0; alpha < _grid.GetSubGrid(ig).nx(); alpha++)
          s[ig][alpha] += j[jsmap[ig][alpha]];
      }

    // Return distribution object
    return Distribution{_grid, s, j};
  }

  //_________________________________________________________________________
  Operator& Operator::operator *= (Operator const& o)
  {
    // Fast method to check that we are using the same Grid
    if (&_grid != &o.GetGrid())
      throw std::runtime_error(error("Operator::operator *=", "Operators grid does not match"));

    const std::vector<matrix<double>> v = _Operator;
    for (int ig = 0; ig < (int) v.size(); ig++)
      {
        // Set operator entries to zero
        _Operator[ig].set(0);

        // The product between the operators is done exploiting the
        // symmetry of the operators.
        for (int alpha = 0; alpha < (int) _Operator[ig].size(0); alpha++)
          for (int beta = 0; beta < (int) _Operator[ig].size(1); beta++)
            for (int gamma = 0; gamma <= beta; gamma++)
              _Operator[ig](alpha, beta) += v[ig](alpha, gamma) * o._Operator[ig](0, beta - gamma);
      }
    return *this;
  }

  //_________________________________________________________________________
  Operator& Operator::operator = (Operator const& o)
  {
    if (this != &o)
      *this = o;
    return *this;
  }

  //_________________________________________________________________________
  Operator& Operator::operator *= (double const& s)
  {
    for (int ig = 0; ig < (int) _Operator.size(); ig++)
      for (int alpha = 0; alpha < (int) _Operator[ig].size(0); alpha++)
        for (int beta = 0; beta < (int) _Operator[ig].size(1); beta++)
          _Operator[ig](alpha, beta) *= s;

    return *this;
  }

  //_________________________________________________________________________
  Operator& Operator::operator *= (std::function<double(double const&)> f)
  {
    for (int ig = 0; ig < (int) _Operator.size(); ig++)
      {
        // Get ig-th subgrid
        const std::vector<double>& sg = _grid.GetSubGrid(ig).GetGrid();
        for (int alpha = 0; alpha < (int) _Operator[ig].size(0); alpha++)
          for (int beta = 0; beta < (int) _Operator[ig].size(1); beta++)
            _Operator[ig](alpha, beta) *= f(sg[beta]);
      }
    return *this;
  }

  //_________________________________________________________________________
  Operator& Operator::operator /= (double const& s)
  {
    for (int ig = 0; ig < (int) _Operator.size(); ig++)
      for (int alpha = 0; alpha < (int) _Operator[ig].size(0); alpha++)
        for (int beta = 0; beta < (int) _Operator[ig].size(1); beta++)
          _Operator[ig](alpha, beta) /= s;

    return *this;
  }

  //_________________________________________________________________________
  Operator& Operator::operator += (Operator const& o)
  {
    // Fast method to check that we are using the same Grid
    if (&_grid != &o.GetGrid())
      throw std::runtime_error(error("Operator::operator +=", "Operators grid does not match"));

    for (int ig = 0; ig < (int) _Operator.size(); ig++)
      for (int alpha = 0; alpha < (int) _Operator[ig].size(0); alpha++)
        for (int beta = 0; beta < (int) _Operator[ig].size(1); beta++)
          _Operator[ig](alpha, beta) += o._Operator[ig](alpha, beta);

    return *this;
  }

  //_________________________________________________________________________
  Operator& Operator::operator -= (Operator const& o)
  {
    // Fast method to check that we are using the same Grid
    if (&_grid != &o.GetGrid())
      throw std::runtime_error(error("Operator::operator +=", "Operators grid does not match"));

    for (int ig = 0; ig < (int) _Operator.size(); ig++)
      for (int alpha = 0; alpha < (int) _Operator[ig].size(0); alpha++)
        for (int beta = 0; beta < (int) _Operator[ig].size(1); beta++)
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
    os << "Operator on the SubGrids:" << "\n";
    const std::ostringstream default_format;
    os << std::scientific;
    os.precision(2);
    for (int i = 0; i < (int) om.size(); i++)
      {
        os << "O[" << i << "]: [";
        for (int alpha = 0; alpha < (int) om[i].size(0); alpha++)
          for (int beta = 0; beta < (int) om[i].size(1); beta++)
            os << "{(" << alpha << ", " << beta << ") : " << om[i](alpha, beta) << "} ";
        os << "]\n";
      }
    os.copyfmt(default_format);
    return os;
  }
}
