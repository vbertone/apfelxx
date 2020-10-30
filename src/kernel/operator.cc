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
  Operator::Operator(Grid const& gr, Expression const& expr, double const& eps, bool const& erbl):
    _grid(gr),
    _erbl(erbl)
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
        _Operator[ig].resize((_erbl ? 2 : 1 ) * nx + 1, 0, (_erbl ? - nx : 0 ));

        // The index beta is set to zero corresponding to the first
        // line of the operator matrix. We keep it explicit to make
        // contact with the notes.
        const int beta = 0;

        // Set xg[beta] as external variable for the computation of
        // the operator (this is in general not needed but in the case
        // of GPD evolution).
        expr.SetExternalVariable(xg[beta]);

        // Log of the lower bound
        const double lxe = ( beta - nx ) * s - log(eta);
        for (int alpha = beta; alpha <= nx; alpha++)
          {
            // Weight of the subtraction term (it is a delta if
            // eta = 1).
            const double ws = li.InterpolantLog(alpha, lxe, sg);

            // Given that the interpolation functions have
            // discontinuos derivative at the nodes, it turns out that
            // it is convenient to split the integrals into id + 1
            // intervals on each of which the integrand is smooth.
            // Even though more integrals have to be computed, the
            // integration converges faster and is more accurate.

            // Initialise integral
            double I = 0;
            for (int j = std::max(0, alpha + 1 - nx); j < (_erbl ? id : std::min(id, alpha - beta) ) + 1; j++)
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
            _Operator[ig][alpha] = I + L * ws;

            // If the operator is ERB-like compute the negative
            // components of the operator vector. This is done in the
            // exact same way as above but running over the first
            // column of the operator matrix rather than over the
            // first line. In practice, I've taken the code above and
            // exchanged alpha and beta. Therfore, setting beta to
            // zero and running over alpha is equivalent to computing
            // the first column.
            if (_erbl)
              // Start from one because the alpha = 0 has already been
              // computed above.
              for (int alpha = 1; alpha <= nx; alpha++)
                {
                  expr.SetExternalVariable(xg[alpha]);
                  const double ws = li.InterpolantLog(beta, ( alpha - nx ) * s - log(eta), sg);
                  double I = 0;
                  for (int j = std::max(0, beta + 1 - nx); j < id + 1; j++)
                    {
                      const Integrator Ij{[&] (double const& y) -> double
                        {
                          const double z  = y / eta;
                          const double wr = li.InterpolantLog(beta, log(xg[alpha] / y), sg);
                          return expr.Regular(z) * wr + expr.Singular(z) * ( wr - ws );
                        }};
                      I += Ij.integrate(exp((alpha - beta + j - 1) * s), exp((alpha - beta + j) * s), eps);
                    }
                  _Operator[ig][-alpha] = I + L * ws;
                }
          }
      }
  }

  //_________________________________________________________________________
  Operator::Operator(Grid const& gr, Expression const& expr, bool const& erbl, double const& eps):
    Operator(gr, expr, eps, erbl)
  {
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
    if (&_grid != &d.GetGrid())
      throw std::runtime_error(error("Operator::operator*=", "Operator and Distribution grids do not match"));

    // Get vector of distributions on the subgrids
    const std::vector<std::vector<double>>& dsg = d.GetDistributionSubGrid();
    const int ng = _grid.nGrids();

    // Initialise output vectors
    std::vector<std::vector<double>> s(ng);
    std::vector<double> j(d.GetDistributionJointGrid().size(), 0);

    // Compute the the distribution on the subgrids
    int i = 0;
    for (int ig = 0; ig < ng; ig++)
      {
        // Resize output on the "ig"-th subgrid
        s[ig].resize(dsg[ig].size(), 0);

        // Get current subgrid
        const SubGrid& sg = _grid.GetSubGrid(ig);

        // Number of interval in the definition region
        int const nx = sg.nx();

        // The product between the operator and the distribution is
        // done exploiting the symmetry of the operator.
        for (int alpha = 0; alpha <= nx; alpha++)
          for (int beta = (_erbl ? 0 : alpha); beta <= nx; beta++)
            s[ig][alpha] += _Operator[ig][beta - alpha] * dsg[ig][beta];

        // Compute the distribution on the joint grid
        const double xtrans = ( ig < ng - 1 ? _grid.GetSubGrid(ig + 1).xMin() : 1 + 2 * eps12 );
        for (int alpha = 0; alpha <= nx && xtrans - sg.GetGrid()[alpha] > eps12; alpha++)
          j[i++] = s[ig][alpha];
      }
    return Distribution{d, s, j};
  }

  //_________________________________________________________________________
  Operator& Operator::operator *= (Operator const& o)
  {
    // Fast method to check that we are using the same Grid
    if (&_grid != &o.GetGrid())
      throw std::runtime_error(error("Operator::operator*=", "Operators grid does not match"));

    const std::vector<ExtendedVector<double>> v = _Operator;
    for (int ig = 0; ig < (int) v.size(); ig++)
      {
        // Set operator entries to zero
        std::fill(_Operator[ig].begin(), _Operator[ig].end(), 0);

        // The product between the operators is done exploiting the
        // symmetry of the operators.
        int const nx = _grid.GetSubGrid(ig).nx();
        for (int beta = _Operator[ig].min(); beta < _Operator[ig].max(); beta++)
          for (int gamma = std::min(0, beta) ; gamma <= (_erbl ? nx : beta) + std::min(0, beta); gamma++)
            _Operator[ig][beta] += v[ig][gamma] * o._Operator[ig][beta - gamma];
      }
    return *this;
  }

  //_________________________________________________________________________
  Operator& Operator::operator *= (double const& s)
  {
    for (int ig = 0; ig < (int) _Operator.size(); ig++)
      for (int alpha = _Operator[ig].min(); alpha < _Operator[ig].max(); alpha++)
        _Operator[ig][alpha] *= s;

    return *this;
  }

  //_________________________________________________________________________
  Operator& Operator::operator *= (std::function<double(double const&)> f)
  {
    for (int ig = 0; ig < (int) _Operator.size(); ig++)
      {
        // Get ig-th subgrid
        const std::vector<double>& sg = _grid.GetSubGrid(ig).GetGrid();
        for (int alpha = _Operator[ig].min(); alpha < _Operator[ig].max(); alpha++)
          _Operator[ig][alpha] *= f(sg[alpha]);
      }
    return *this;
  }

  //_________________________________________________________________________
  Operator& Operator::operator /= (double const& s)
  {
    for (int ig = 0; ig < (int) _Operator.size(); ig++)
      for (int alpha = _Operator[ig].min(); alpha < _Operator[ig].max(); alpha++)
        _Operator[ig][alpha] /= s;

    return *this;
  }

  //_________________________________________________________________________
  Operator& Operator::operator += (Operator const& o)
  {
    // Fast method to check that we are using the same Grid
    if (&_grid != &o.GetGrid())
      throw std::runtime_error(error("Operator::operator+=", "Operators grid does not match"));

    for (int ig = 0; ig < (int) _Operator.size(); ig++)
      for (int alpha = _Operator[ig].min(); alpha < _Operator[ig].max(); alpha++)
        _Operator[ig][alpha] += o._Operator[ig][alpha];

    return *this;
  }

  //_________________________________________________________________________
  Operator& Operator::operator -= (Operator const& o)
  {
    // Fast method to check that we are using the same Grid
    if (&_grid != &o.GetGrid())
      throw std::runtime_error(error("Operator::operator+=", "Operators grid does not match"));

    for (int ig = 0; ig < (int) _Operator.size(); ig++)
      for (int alpha = _Operator[ig].min(); alpha < _Operator[ig].max(); alpha++)
        _Operator[ig][alpha] -= o._Operator[ig][alpha];

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
    const std::vector<ExtendedVector<double>> om = op.GetOperator();
    os << "Operator: " << &op << "\n";
    os << "Operator on the SubGrids:" << "\n";
    const std::ostringstream default_format;
    os << std::scientific;
    os.precision(2);
    for (int i = 0; i < (int) om.size(); i++)
      {
        os << om[i].min() << "\n";
        os << "O[" << i << "]: [";
        for (int alpha = om[i].min(); alpha < (int) om[i].max(); alpha++)
          os << "{" << alpha << " : " << om[i][alpha] << "} ";
        os << "]\n";
      }
    os.copyfmt(default_format);
    return os;
  }
}
