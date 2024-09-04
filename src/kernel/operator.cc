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
  Operator::Operator(Grid const& gr, Expression const& expr, double const& eps, bool const& gpd):
    _grid(gr),
    _expr(expr),
    _eps(eps),
    _gpd(gpd)
  {
    if (_gpd)
      BuildOperatorGPD();
    else
      BuildOperatorDGLAP();
  }

  //_________________________________________________________________________
  void Operator::BuildOperatorDGLAP()
  {
    // Interpolator object for the interpolating functions
    const LagrangeInterpolator li{_grid};

    // Scaling factor
    const double eta = _expr.eta();

    // Number of grids
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
        const double L = eta * _expr.Local(1 / exp(s) / eta);

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
                // id + 1 intervals over each of which the integrand
                // is smooth. Despite more integrals have to be
                // computed, the integration converges faster and is
                // more accurate.

                // Initialise integral
                double I = 0;

                // Run over the grid intervals over which the
                // interpolating function is different from zero.
                for (int j = 0; j < std::min(id, alpha - beta) + 1; j++)
                  {
                    // Define "Integrator" object
                    const Integrator Ij{[&] (double const& y) -> double
                      {
                        const double z  = y / eta;
                        const double wr = li.InterpolantLog(alpha, log(xg[beta] / y), sg);
                        return _expr.Regular(z) * wr + _expr.Singular(z) * ( wr - ws );
                      }};
                    // Compute the integral
                    I += Ij.integrate(exp((beta - alpha + j - 1) * s), exp((beta - alpha + j) * s), _eps);
                  }
                // Add the local part
                _Operator[ig](beta, alpha) = I + L * ws;
              }
          }
      }
  }

  //_________________________________________________________________________
  void Operator::BuildOperatorGPD()
  {
    // Interpolator object for the interpolating functions
    const LagrangeInterpolator li{_grid};

    // Get skewness from the expression
    const double xi = 1 / _expr.eta();

    // Loop over the subgrids
    _Operator.resize(1);

    // Get joint grid
    const SubGrid& jg = _grid.GetJointGrid();

    // Get vector with the grid nodes
    const std::vector<double>& xg = jg.GetGrid();

    // Get number of grid intervals in the definition range
    const int nx = jg.nx();

    // Get interpolation degree
    const int k = jg.InterDegree();

    // Initialise operator
    _Operator[0].resize(nx, nx, 0);

    // Loop over the index beta. In fact beta = 0 because the size
    // of the first dimension of "_Operator" is one.
    for (int beta = 0; beta < (int) _Operator[0].size(0); beta++)
      {
        // Set xg[beta] as external variable for the computation
        // of the operator (this is in general not needed but in
        // the case of GPD evolution).
        _expr.SetExternalVariable(xg[beta]);

        // Define "kappa" variable
        const double kappa = xi / xg[beta];

        // Loop over the index alpha
        for (int alpha = 0; alpha < (int) _Operator[0].size(1); alpha++)
          {
            // Weight of the subtraction term: this is
            // \delta_{\alpha\beta}
            const double ws = (beta == alpha ? 1 : 0);

            // Weight of the PV subtraction terms at y = 1 / kappa
            const double wspv = li.Interpolant(alpha, xi, jg);

            // Given that the interpolation functions have
            // discontinuos derivative at the nodes, it turns out that
            // it is convenient to split the integrals into k + 1
            // intervals over each of which the integrand is
            // smooth. Despite more integrals have to be computed, the
            // integration converges faster and is more accurate.

            // Run over the grid intervals over which the
            // interpolating function is different from zero.
            for (int j = 0; j <= std::min(alpha, k); j++)
              {
                // Define "Integrator" object. IMPORTANT: the
                // particular form of the subtraction term only
                // applies to singular terms that behave as
                // 1/(1-y). If other singular functions are used, this
                // term has to be adjusted. The integral also contains
                // a possible principal-valued contribution with
                // singularity at y = 1 / kappa integrated between x
                // and 1. Also in this case the procedure is specific
                // for 1 / ( 1 - kappa y).
                const Integrator Ij{[&] (double const& y) -> double
                  {
                    const double wr = li.Interpolant(alpha, xg[beta] / y, jg);
                    const double ky = kappa * y;
                    return _expr.Regular(y) * wr
                                + _expr.Singular(y) * ( wr - ws * ( 1 + (y > 1 ? ( 1 - y ) / y : 0) ) )
                                + _expr.SingularPV(y) * ( wr - wspv * ( 1 + (ky > 1 ? ( 1 - ky ) / ky : 0) ) );
                  }};
                // Compute the integral
                _Operator[0](beta, alpha) += Ij.integrate(xg[beta] / xg[alpha - j + 1], xg[beta] / xg[alpha - j], _eps);
              }
            // Add PV local part from the y = 1 / kappa singularity
            _Operator[0](beta, alpha) += wspv * ( _expr.LocalPV(xg[beta] / xg[alpha + 1]) - _expr.LocalPV(xg[std::max(0, alpha - k)] / xg[beta] / pow(kappa, 2)) );
          }
        // Add the local parts: that from standard +-prescripted terms
        // ("Local") and that deriving from principal-valued integrals
        // at x = 1, i.e. from the ++-prescription ("LocalPP").
        _Operator[0](beta, beta) += _expr.Local(xg[beta] / xg[beta + 1])
                                    + _expr.LocalPP(xg[beta] / xg[beta + 1])
                                    - (beta == 0 ? 0 : _expr.LocalPP(xg[std::max(0, beta - k)] / xg[beta]));
      }
  }

  //_________________________________________________________________________
  Distribution Operator::Evaluate(double const& x) const
  {
    // Initialise distribution
    apfel::Distribution d{_grid};

    // Number of sub-grids
    const int ng = _grid.nGrids();

    // Initialise vector of vectors containing the distribution on the
    // sub-grids.
    std::vector<std::vector<double>> dsg(ng);

    // Fill in sub-grids
    for (int ig = 0; ig < ng; ig++)
      {
        // Number of nodes
        const int nx = _grid.GetSubGrid(ig).nx();

        // Resize distribution subgrid
        dsg[ig].resize(_grid.GetSubGrid(ig).GetGrid().size(), 0.);

        // Get summation bounds for the interpolation
        const std::array<int, 2> bounds = d.SumBounds(x, _grid.GetSubGrid(ig));

        // Run over the second index of the operator
        for (int alpha = 0; alpha < nx; alpha++)
          {
            // If the operator has one single line, this means that
            // one can (must) use the symmetry _Operator[ig](beta,
            // alpha) = _Operator[ig](0, alpha - beta). Otherwise, the
            // operator matrix is assumed to be a full nx^2 matrix.
            if (_Operator[ig].size(0) == 1)
              for (int beta = bounds[0]; beta < std::min(bounds[1], alpha + 1); beta++)
                dsg[ig][alpha] += d.Interpolant(beta, x, _grid.GetSubGrid(ig)) * _Operator[ig](0, alpha - beta);
            else
              for (int beta = bounds[0]; beta < std::min(bounds[1], nx + 1); beta++)
                dsg[ig][alpha] += d.Interpolant(beta, x, _grid.GetSubGrid(ig)) * _Operator[ig](beta, alpha);
          }
      }

    // Initialise vector containing the distribution on the joint
    // grid.
    std::vector<double> djg(_grid.GetJointGrid().GetGrid().size(), 0.);

    // Determine sub-grid to use to fill in the joint distribution
    int jg;
    for (jg = ng - 1; jg >= 0; jg--)
      if (x >= _grid.GetSubGrid(jg).xMin())
        break;

    // Determine the weight of the previous grid that enforces a
    // smoother transition from the (jg-1)-th and the jg-th subgrid.
    // Do it only from the second grid on.
    double w2 = 0;
    if (jg > 0)
      {
        // Get interpolation degree on the jg-th subgrid
        const int id = _grid.GetSubGrid(jg).InterDegree();

        // Determine the closest bottom node on the jg-th subgrid
        const int ndb = d.SumBounds(x, _grid.GetSubGrid(jg))[0];

        // If "ndb" is larger than the interpolation degree of the jg-th
        // subgrid "id", use only this subgrid to compute the distribution
        // on the joint grid. Otherwise compute the joint grid as a
        // weighted average of (jg-1)-th and jg-th subgrids with weight
        // depending on the distance from the bottom node of the gj-th
        // subgrid. This should ensure a smoother transition from the
        // (jg-1)-th and the jg-th subgrid.
        w2 = std::max(id - ndb, 0) * pow(id + 1, -1);

        // Get map of indices from joint to the sub-grid preceding the
        // selected one.
        const std::vector<int>& jsmapp = _grid.JointToSubMap()[jg-1];

        // Fill in joint grid
        for (int alpha = 0; alpha < (int) jsmapp.size(); alpha++)
          djg[jsmapp[alpha]] += w2 * dsg[jg-1][alpha];
      }

    // Determine weight on the current grid according to the weight on
    // the previous grid.
    const double w1 = 1 - w2;

    // Get map of indices map from joint to the selected sub-grid
    const std::vector<int>& jsmap = _grid.JointToSubMap()[jg];

    // Fill in joint grid
    for (int alpha = 0; alpha < (int) jsmap.size(); alpha++)
      djg[jsmap[alpha]] += w1 * dsg[jg][alpha];

    // Set sub-grids and joint grid
    d.SetSubGrids(dsg);
    d.SetJointGrid(djg);

    // Return distribution
    return d;
  }

  //_________________________________________________________________________
  Distribution Operator::operator *= (Distribution const& d) const
  {
    // Fast method to check that we are using the same Grid
    if (&_grid != &d.GetGrid())
      throw std::runtime_error(error("Operator::operator *=", "Operator and Distribution grids do not match."));

    // Get number of subgrids
    const int ng = _grid.nGrids();

    // Get map of indices map from joint to subgrids
    const std::vector<std::vector<int>>& jsmap = _grid.JointToSubMap();

    // Get joint distribution
    const std::vector<double>& dj = d.GetDistributionJointGrid();

    // Initialise output vectors
    std::vector<double> j(dj.size(), 0);
    std::vector<std::vector<double>> s(ng);

    // Get number of grid intervals in the definition range of the
    // joint grid.
    const int nx = _grid.GetJointGrid().nx();

    // Get map of indices map from sub to joint to grids
    const std::vector<std::pair<int, int>>& sjmap = _grid.SubToJointMap();

    // Construct joint distribution first. The product between the
    // operator and the distribution is done exploiting the symmetry
    // of the operator if the first operator has one line
    // only. Otherwise the product is done in a standard way. This
    // should be enough to distinguish between DGLAP- and GPD-like
    // operators.
    for (int beta = 0; beta < nx; beta++)
      if (_Operator[0].size(0) == 1)
        {
          const std::pair<int, int> m = sjmap[beta];
          for (int alpha = m.second; alpha < _grid.GetSubGrid(m.first).nx(); alpha++)
            j[beta] += _Operator[m.first](0, alpha - m.second) * dj[jsmap[m.first][alpha]];
        }
      else
        for (int alpha = 0; alpha < nx; alpha++)
          j[beta] += _Operator[0](beta, alpha) * dj[alpha];

    // Compute the the distribution on the subgrids
    for (int ig = 0; ig < ng; ig++)
      {
        // Resize output on the "ig"-th subgrid
        s[ig].resize(d.GetDistributionSubGrid()[ig].size(), 0.);
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
      throw std::runtime_error(error("Operator::operator *=", "Grids do not match."));

    const std::vector<matrix<double>> v = _Operator;
    for (int ig = 0; ig < (int) v.size(); ig++)
      {
        // Set operator entries to zero
        _Operator[ig].set(0);

        // The product between the operators is done exploiting the
        // symmetry of the operators if the operator matrix only
        // contains one line. Otherwise the product is done in a
        // standard way. This should be enough to distinguish between
        // DGLAP- and GPD-like operators.
        if (_Operator[ig].size(0) == 1)
          for (int beta = 0; beta < (int) _Operator[ig].size(1); beta++)
            for (int gamma = 0; gamma <= beta; gamma++)
              _Operator[ig](0, beta) += v[ig](0, gamma) * o._Operator[ig](0, beta - gamma);
        else
          for (int alpha = 0; alpha < (int) _Operator[ig].size(0); alpha++)
            for (int beta = 0; beta < (int) _Operator[ig].size(1); beta++)
              for (int gamma = 0; gamma < (int) _Operator[ig].size(1); gamma++)
                _Operator[ig](alpha, beta) += v[ig](alpha, gamma) * o._Operator[ig](gamma, beta);
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
      throw std::runtime_error(error("Operator::operator +=", "Grids do not match."));

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
      throw std::runtime_error(error("Operator::operator +=", "Grids do not match."));

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
