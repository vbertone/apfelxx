//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/doubleoperator.h"
#include "apfel/operator.h"
#include "apfel/integrator.h"
#include "apfel/messages.h"
#include "apfel/constants.h"

#include <cmath>

namespace apfel
{
  //_________________________________________________________________________
  DoubleOperator::DoubleOperator(Grid const& gr1, Grid const& gr2, DoubleExpression const& dexpr, double const& eps):
    _grid1(gr1),
    _grid2(gr2),
    _eps(eps),
    _dexprName(dexpr.GetName())
  {
    // Interpolator objects for the interpolating functions
    const LagrangeInterpolator li1{_grid1};
    const LagrangeInterpolator li2{_grid2};

    // Number of sub grids
    const int ng1 = _grid1.nGrids();
    const int ng2 = _grid2.nGrids();

    // Resize operator w.r.t. the first grid
    _dOperator.resize(ng1);

    // Loop over the subgrids of the first grid
    for (int ig1 = 0; ig1 < ng1; ig1++)
      {
        // Get current subgrid parameters on the first grid
        const SubGrid& sg1 = _grid1.GetSubGrid(ig1);
        const std::vector<double>& xg1 = sg1.GetGrid();
        const int nx1 = sg1.nx();
        const int id1 = sg1.InterDegree();
        const double dx1 = sg1.Step();
        const double edx1 = exp(-dx1);

        // Resize operator w.r.t. the second grid
        _dOperator[ig1].resize(ng2);

        // Loop over the subgrids of the second grid
        for (int ig2 = 0; ig2 < ng2; ig2++)
          {
            // Get current subgrid parameters on the second grid
            const SubGrid& sg2 = _grid2.GetSubGrid(ig2);
            const std::vector<double>& xg2 = sg2.GetGrid();
            const int nx2 = sg2.nx();
            const int id2 = sg2.InterDegree();
            const double dx2 = sg2.Step();
            const double edx2 = exp(-dx2);

            // Evaluate Local-Local contribution
            const double cLL = dexpr.LocalLocal(edx1, edx2);

            // Evaluate Local-Singular and Local-Regular making the
            // delta-function of the first variable act and
            // integrating over the second variable.
            matrix<double> cLS_LR{1, (size_t) nx2};
            for (int delta = 0; delta < (int) cLS_LR.size(0); delta++)
              for (int gamma = delta; gamma < (int) cLS_LR.size(1); gamma++)
                {
                  const double ws2 = (delta == gamma ? 1 : 0);
                  for (int k = 0; k < std::min(id2, gamma - delta) + 1; k++)
                    {
                      const double b2 = exp((delta - gamma + k) * dx2);
                      const double a2 = b2 * edx2;
                      const Integrator Ik{[&] (double const& y2) -> double
                        {
                          const double wr2 = li2.InterpolantLog(gamma, log(xg2[delta] / y2), sg2);
                          return dexpr.LocalSingular(edx1, y2) * ( wr2 - ws2 ) + dexpr.LocalRegular(edx1, y2) * wr2;
                        }};
                      cLS_LR(delta, gamma) += Ik.integrate(a2, b2, _eps);
                    }
                }

            // Evaluate Singular-Local and Regular-Local making the
            // delta-function of the second variable act and
            // integrating over the first variable.
            matrix<double> cSL_RL{1, (size_t) nx1};
            for (int beta = 0; beta < (int) cSL_RL.size(0); beta++)
              for (int alpha = beta; alpha < (int) cSL_RL.size(1); alpha++)
                {
                  const double ws1 = (beta == alpha ? 1 : 0);
                  for (int j = 0; j < std::min(id1, alpha - beta) + 1; j++)
                    {
                      const double b1 = exp((beta - alpha + j) * dx1);
                      const double a1 = b1 * edx1;
                      const Integrator Ij{[&] (double const& y1) -> double
                        {
                          const double wr1 = li1.InterpolantLog(alpha, log(xg1[beta] / y1), sg1);
                          return ( wr1 - ws1 ) * dexpr.SingularLocal(y1, edx2) + wr1 * dexpr.RegularLocal(y1, edx2);
                        }};
                      cSL_RL(beta, alpha) += Ij.integrate(a1, b1, _eps);
                    }
                }

            // Finally move to compute Singular-Singular,
            // Singular-Regular, Regular-Singular, and Regular-Regular
            // that require integrating over both y1 and y2. Include
            // Local-Local, Local-Singular, Local-Regular,
            // Singular-Local, and Regular-Local along the way.

            // Resize operator container w.r.t. the first variable
            _dOperator[ig1][ig2].resize(1, nx1);

            // Loop over the index beta. In fact, beta = 0 because the
            // size of the first dimension of "_dOperator" is one.
            for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(0); beta++)
              {
                // Useful precomputation
                const double lxg1beta = log(xg1[beta]);

                // Loop over the index alpha
                for (int alpha = beta; alpha < (int) _dOperator[ig1][ig2].size(1); alpha++)
                  {
                    // Resize operator w.r.t. the second variable
                    _dOperator[ig1][ig2](beta, alpha).resize(1, nx2);

                    // Weight of the subtraction term w.r.t. the first
                    // grid
                    const double ws1 = (beta == alpha ? 1 : 0);

                    // Starting value for the integration in y1
                    const double b10 = exp(( beta - alpha ) * dx1);

                    // Loop over the index delta. In fact, delta = 0
                    // because the size of the first dimension of
                    // "_dOperator" is one.
                    for (int delta = 0; delta < (int) _dOperator[ig1][ig2](beta, alpha).size(0); delta++)
                      {
                        // Useful precomputation
                        const double lxg2delta = log(xg2[delta]);

                        // Loop over the index gamma
                        for (int gamma = delta; gamma < (int) _dOperator[ig1][ig2](beta, alpha).size(1); gamma++)
                          {
                            // Weight of the subtraction term
                            // w.r.t. the second grid.
                            const double ws2 = (delta == gamma ? 1 : 0);

                            // Starting value for the integration in y2
                            const double b20 = exp(( delta - gamma ) * dx2);

                            // Initialise integral
                            double cSS_SR_RS_RR = 0;

                            // Run over the grid intervals that
                            // concern the integral.
                            double b1 = b10;
                            for (int j = 0; j < std::min(id1, alpha - beta) + 1; j++)
                              {
                                double b2 = b20;
                                for (int k = 0; k < std::min(id2, gamma - delta) + 1; k++)
                                  {
                                    const Integrator Iy1{[&] (double const& y1) -> double
                                      {
                                        if (y1 > 1 - eps8)
                                          return 0;
                                        const double wr1 = li1.InterpolantLog(alpha, lxg1beta - log(y1), sg1);
                                        const Integrator Iy1y2{
                                          [&] (double const& y2) -> double
                                          {
                                            if (y2 > 1 - eps8)
                                              return 0;
                                            const double wr2 = li2.InterpolantLog(gamma, lxg2delta - log(y2), sg2);
                                            return ( wr1 - ws1 ) * dexpr.SingularSingular(y1, y2) * ( wr2 - ws2 )
                                                                        + ( wr1 - ws1 ) * dexpr.SingularRegular(y1, y2) * wr2
                                                                        + wr1 * dexpr.RegularSingular(y1, y2) * ( wr2 - ws2 )
                                                                        + wr1 * dexpr.RegularRegular(y1, y2) * wr2;
                                          }};
                                        return Iy1y2.integrate(b2 * edx2, b2, _eps);
                                      }};
                                    cSS_SR_RS_RR += Iy1.integrate(b1 * edx1, b1, _eps);
                                    b2 /= edx2;
                                  }
                                b1 /= edx1;
                              }
                            _dOperator[ig1][ig2](beta, alpha)(delta, gamma) = ws1 * cLL * ws2 + ws1 * cLS_LR(delta, gamma) + cSL_RL(beta, alpha) * ws2 + cSS_SR_RS_RR;
                          }
                      }
                  }
              }
          }
      }
  }

  //_________________________________________________________________________
  DoubleOperator::DoubleOperator(Operator const& O1, Operator const& O2):
    _grid1(O1.GetGrid()),
    _grid2(O2.GetGrid()),
    _eps(std::max(O1.GetIntegrationAccuracy(), O2.GetIntegrationAccuracy())),
    _dexprName("None")
  {
    // Stop the computation if any of the two operators is of GPD
    // type.
    if (O1.IsGPD() || O2.IsGPD())
      throw std::runtime_error(error("DoubleOperator::DoubleOperator", "GPD double operators cannot be handled yet."));

    // Get single operators
    const std::vector<matrix<double>> o1 = O1.GetOperator();
    const std::vector<matrix<double>> o2 = O2.GetOperator();

    // Get number of grids
    const int ng1 = _grid1.nGrids();
    const int ng2 = _grid2.nGrids();

    // Fill in double operator
    _dOperator.resize(ng1);
    for (int ig1 = 0; ig1 < ng1; ig1++)
      {
        const int nx1 = _grid1.GetSubGrid(ig1).nx();
        _dOperator[ig1].resize(ng2);
        for (int ig2 = 0; ig2 < ng2; ig2++)
          {
            const int nx2 = _grid2.GetSubGrid(ig2).nx();
            _dOperator[ig1][ig2].resize(1, nx1);
            for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(0); beta++)
              for (int alpha = beta; alpha < (int) _dOperator[ig1][ig2].size(1); alpha++)
                {
                  _dOperator[ig1][ig2](beta, alpha).resize(1, nx2);
                  for (int delta = 0; delta < (int) _dOperator[ig1][ig2](alpha, beta).size(0); delta++)
                    for (int gamma = delta; gamma < (int) _dOperator[ig1][ig2](alpha, beta).size(1); gamma++)
                      _dOperator[ig1][ig2](beta, alpha)(delta, gamma) = o1[ig1](beta, alpha) * o2[ig2](delta, gamma);
                }
          }
      }
  }

  //_________________________________________________________________________
  DoubleOperator::DoubleOperator(DoubleObject<Operator> const& DObj):
    _grid1(DObj.GetTerms()[0].object1.GetGrid()),
    _grid2(DObj.GetTerms()[0].object2.GetGrid()),
    _eps(std::max(DObj.GetTerms()[0].object1.GetIntegrationAccuracy(), DObj.GetTerms()[0].object2.GetIntegrationAccuracy())),
    _dexprName("None")
  {
    // Stop the computation if any of the two operators in the first
    // term of the double objects is of GPD type.
    if (DObj.GetTerms()[0].object1.IsGPD() || DObj.GetTerms()[0].object2.IsGPD())
      throw std::runtime_error(error("DoubleOperator::DoubleOperator", "GPD double operators cannot be handled yet."));

    // Number of subgrids
    const int ng1 = _grid1.nGrids();
    const int ng2 = _grid2.nGrids();

    // Resize matrix and set it to zero
    _dOperator.resize(ng1, std::vector<matrix<matrix<double>>>(ng2));
    for (int ig1 = 0; ig1 < ng1; ig1++)
      for (int ig2 = 0; ig2 < ng2; ig2++)
        {
          const int nx1 = _grid1.GetSubGrid(ig1).nx();
          const int nx2 = _grid2.GetSubGrid(ig2).nx();
          _dOperator[ig1][ig2].resize(1, nx1);
          for (int alpha = 0; alpha < nx1; alpha++)
            _dOperator[ig1][ig2](0, alpha).resize(1, nx2);
        }

    // Fill in container
    for (auto const& term : DObj.GetTerms())
      {
        // Get coefficient and objects of the current term
        const double coef = term.coefficient;
        const std::vector<matrix<double>> O1 = term.object1.GetOperator();
        const std::vector<matrix<double>> O2 = term.object2.GetOperator();

        // Fill in matrix
        for (int ig1 = 0; ig1 < ng1; ig1++)
          for (int ig2 = 0; ig2 < ng2; ig2++)
            for(int beta = 0; beta < (int) _dOperator[ig1][ig2].size(0); beta++)
              for(int alpha = 0; alpha < (int) _dOperator[ig1][ig2].size(1); alpha++)
                for(int delta = 0; delta < (int) _dOperator[ig1][ig2](beta, alpha).size(0); delta++)
                  for(int gamma = 0; gamma < (int) _dOperator[ig1][ig2](beta, alpha).size(1); gamma++)
                    _dOperator[ig1][ig2](beta, alpha)(delta, gamma) += coef * O1[ig1](beta, alpha) * O2[ig2](delta, gamma);
      }
  }

#if WITH_YAML_CPP == 1
  //_________________________________________________________________________
  DoubleOperator::DoubleOperator(YAML::Node const& Node, Grid const& gr1, Grid const& gr2, DoubleExpression const& dexpr):
    _grid1(gr1),
    _grid2(gr2),
    _eps(Node["Integration accuracy"].as<double>()),
    _dexprName(dexpr.GetName())
  {
    // Check that the DoubleExpression stored in the node matches with
    // the input one by checking the names.
    if (_dexprName != Node["DoubleExpression"].as<std::string>())
      throw std::runtime_error(error("DoubleOperator::DoubleOperator", "Input DoubleExpression and DoubleExpression name stored in YAML::Node do not match."));

    // First construct interpolation grids using the parameters give
    // in the node ...
    std::vector<apfel::SubGrid> vsg1;
    for (auto const &sg1 : Node["FirstGrid"])
      vsg1.push_back({sg1[0].as<int>(), sg1[1].as<double>(), sg1[2].as<int>()});
    std::vector<apfel::SubGrid> vsg2;
    for (auto const &sg2 : Node["SecondGrid"])
      vsg2.push_back({sg2[0].as<int>(), sg2[1].as<double>(), sg2[2].as<int>()});
    const apfel::Grid g1{vsg1};
    const apfel::Grid g2{vsg2};

    // ... then check that they match with those given as input.
    if (g1 != _grid1 || g2 != _grid2)
      throw std::runtime_error(error("DoubleOperator::DoubleOperator", "Input grids and grids stored in YAML::Node do not match."));

    // Get double operator values from the node
    const std::vector<std::vector<std::vector<std::vector<double>>>> NodeOperator = Node["DoubleOperator"].as<std::vector<std::vector<std::vector<std::vector<double>>>>>();

    // Get number of subgrids
    const int ng1 = _grid1.nGrids();
    const int ng2 = _grid2.nGrids();

    // Resize matrix and set it to zero
    _dOperator.resize(ng1, std::vector<matrix<matrix<double>>>(ng2));
    for (int ig1 = 0; ig1 < ng1; ig1++)
      for (int ig2 = 0; ig2 < ng2; ig2++)
        {
          const int nx1 = _grid1.GetSubGrid(ig1).nx();
          const int nx2 = _grid2.GetSubGrid(ig2).nx();
          _dOperator[ig1][ig2].resize(1, nx1);
          for (int alpha = 0; alpha < nx1; alpha++)
            _dOperator[ig1][ig2](0, alpha) = matrix<double> {1, (size_t) nx2, NodeOperator[ig1][ig2][alpha]};
        }
  }

  //_________________________________________________________________________
  std::string DoubleOperator::EmitDoubleOperator() const
  {
    YAML::Emitter DOTab;

    // Dump DoubleOperator object on a YAML::Emitter object
    DOTab.SetDoublePrecision(8);
    DOTab << YAML::BeginMap;
    DOTab << YAML::Key << "DoubleExpression" << YAML::Value << _dexprName;
    DOTab << YAML::Key << "Integration accuracy" << YAML::Value << _eps;
    DOTab << YAML::Key << "FirstGrid";
    DOTab << YAML::Value << YAML::BeginSeq;
    for (auto const& sg1 : _grid1.GetSubGrids())
      DOTab << YAML::Flow << YAML::BeginSeq << sg1.nx() << sg1.xMin() << sg1.InterDegree() << YAML::EndSeq;
    DOTab << YAML::EndSeq;
    DOTab << YAML::Key << "SecondGrid";
    DOTab << YAML::Value << YAML::BeginSeq;
    for (auto const& sg2 : _grid2.GetSubGrids())
      DOTab << YAML::Flow << YAML::BeginSeq << sg2.nx() << sg2.xMin() << sg2.InterDegree() << YAML::EndSeq;
    DOTab << YAML::EndSeq;
    DOTab << YAML::Key << "DoubleOperator";
    DOTab << YAML::Value << YAML::Flow << YAML::BeginSeq;
    for (auto const& og1 : _dOperator)
      {
        DOTab << YAML::Flow << YAML::BeginSeq;
        for (auto const& og2 : og1)
          {
            DOTab << YAML::Flow << YAML::BeginSeq;
            for (auto const& m1 : og2.data())
              {
                DOTab << YAML::Flow << YAML::BeginSeq;
                for (auto const& m2 : m1.data())
                  DOTab << YAML::Flow << m2;
                DOTab << YAML::EndSeq;
              }
            DOTab << YAML::EndSeq;
          }
        DOTab << YAML::EndSeq;
      }
    DOTab << YAML::EndSeq;
    DOTab << YAML::EndMap;

    // Return DoubleOperator as a string
    return DOTab.c_str();
  }
#endif

  //_________________________________________________________________________
  DoubleDistribution DoubleOperator::operator *= (DoubleDistribution const& d) const
  {
    // Fast method to check that we are using the same grids
    if (&_grid1 != &d.GetFirstGrid() || &_grid2 != &d.GetSecondGrid())
      throw std::runtime_error(error("DoubleOperator::operator *=", "DoubleOperator and DoubleDistribution grids do not match."));

    // Get numbers of subgrids
    const int ng1 = _grid1.nGrids();
    const int ng2 = _grid2.nGrids();

    // Get number of grid intervals in the definition range of the
    // joint grids.
    const int nx1 = _grid1.GetJointGrid().nx();
    const int nx2 = _grid2.GetJointGrid().nx();

    // Get maps of indices map from joint to subgrids
    const std::vector<std::vector<int>>& jsmap1 = _grid1.JointToSubMap();
    const std::vector<std::vector<int>>& jsmap2 = _grid2.JointToSubMap();

    // Get map of indices map from sub to joint grids
    const std::vector<std::pair<int, int>>& sjmap1 = _grid1.SubToJointMap();
    const std::vector<std::pair<int, int>>& sjmap2 = _grid2.SubToJointMap();

    // Get joint distribution of the input distribution
    const matrix<double>& dj = d.GetDistributionJointGrid();

    // Initialise output vectors
    matrix<double> j(dj.size(0), dj.size(1));
    std::vector<std::vector<matrix<double>>> s(ng1, std::vector<matrix<double>>(ng2));

    // Construct joint distributions first. The product between the
    // operator and the distribution is done exploiting the symmetry
    // of the operator, which implies that it has a single line for
    // each pair of indices.
    for (int beta = 0; beta < nx1; beta++)
      {
        const std::pair<int, int> m1 = sjmap1[beta];
        for (int delta = 0; delta < nx2; delta++)
          {
            const std::pair<int, int> m2 = sjmap2[delta];
            for (int alpha = m1.second; alpha < _grid1.GetSubGrid(m1.first).nx(); alpha++)
              for (int gamma = m2.second; gamma < _grid2.GetSubGrid(m2.first).nx(); gamma++)
                j(beta, delta) += _dOperator[m1.first][m2.first](0, alpha - m1.second)(0, gamma - m2.second) * dj(jsmap1[m1.first][alpha], jsmap2[m2.first][gamma]);
          }
      }

    // Compute the the distribution on the subgrids using the joint
    // grids.
    for (int ig1 = 0; ig1 < ng1; ig1++)
      for (int ig2 = 0; ig2 < ng2; ig2++)
        {
          // Resize output on the "ig1-ig2"-th subgrid
          s[ig1][ig2].resize(d.GetDistributionSubGrid()[ig1][ig2].size(0), d.GetDistributionSubGrid()[ig1][ig2].size(1));
          for (int alpha = 0; alpha < _grid1.GetSubGrid(ig1).nx(); alpha++)
            for (int gamma = 0; gamma < _grid2.GetSubGrid(ig2).nx(); gamma++)
              s[ig1][ig2](alpha, gamma) += j(jsmap1[ig1][alpha], jsmap2[ig2][gamma]);
        }

    // Return double distribution object
    return DoubleDistribution{_grid1, _grid2, s, j};
  }

  //_________________________________________________________________________
  DoubleOperator& DoubleOperator::operator *= (double const& s)
  {
    for (int ig1 = 0; ig1 < (int) _dOperator.size(); ig1++)
      for (int ig2 = 0; ig2 < (int) _dOperator[ig1].size(); ig2++)
        for (int alpha = 0; alpha < (int) _dOperator[ig1][ig2].size(0); alpha++)
          for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(1); beta++)
            for (int gamma = 0; gamma < (int) _dOperator[ig1][ig2](alpha, beta).size(0); gamma++)
              for (int delta = 0; delta < (int) _dOperator[ig1][ig2](alpha, beta).size(1); delta++)
                _dOperator[ig1][ig2](alpha, beta)(gamma, delta) *= s;

    return *this;
  }

  //_________________________________________________________________________
  DoubleOperator& DoubleOperator::operator /= (double const& s)
  {
    for (int ig1 = 0; ig1 < (int) _dOperator.size(); ig1++)
      for (int ig2 = 0; ig2 < (int) _dOperator[ig1].size(); ig2++)
        for (int alpha = 0; alpha < (int) _dOperator[ig1][ig2].size(0); alpha++)
          for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(1); beta++)
            for (int gamma = 0; gamma < (int) _dOperator[ig1][ig2](alpha, beta).size(0); gamma++)
              for (int delta = 0; delta < (int) _dOperator[ig1][ig2](alpha, beta).size(1); delta++)
                _dOperator[ig1][ig2](alpha, beta)(gamma, delta) /= s;

    return *this;
  }

  //_________________________________________________________________________
  DoubleOperator& DoubleOperator::operator += (DoubleOperator const& o)
  {
    // Fast method to check that we are using the same grids
    if (&_grid1 != &o.GetFirstGrid() || &_grid2 != &o.GetSecondGrid())
      throw std::runtime_error(error("DoubleOperator::operator +=", "Grids do not match."));

    for (int ig1 = 0; ig1 < (int) _dOperator.size(); ig1++)
      for (int ig2 = 0; ig2 < (int) _dOperator[ig1].size(); ig2++)
        for (int alpha = 0; alpha < (int) _dOperator[ig1][ig2].size(0); alpha++)
          for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(1); beta++)
            for (int gamma = 0; gamma < (int) _dOperator[ig1][ig2](alpha, beta).size(0); gamma++)
              for (int delta = 0; delta < (int) _dOperator[ig1][ig2](alpha, beta).size(1); delta++)
                _dOperator[ig1][ig2](alpha, beta)(gamma, delta) += o._dOperator[ig1][ig2](alpha, beta)(gamma, delta);

    return *this;
  }

  //_________________________________________________________________________
  DoubleOperator& DoubleOperator::operator -= (DoubleOperator const& o)
  {
    // Fast method to check that we are using the same grids
    if (&_grid1 != &o.GetFirstGrid() || &_grid2 != &o.GetSecondGrid())
      throw std::runtime_error(error("DoubleOperator::operator +=", "Grids do not match."));

    for (int ig1 = 0; ig1 < (int) _dOperator.size(); ig1++)
      for (int ig2 = 0; ig2 < (int) _dOperator[ig1].size(); ig2++)
        for (int alpha = 0; alpha < (int) _dOperator[ig1][ig2].size(0); alpha++)
          for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(1); beta++)
            for (int gamma = 0; gamma < (int) _dOperator[ig1][ig2](alpha, beta).size(0); gamma++)
              for (int delta = 0; delta < (int) _dOperator[ig1][ig2](alpha, beta).size(1); delta++)
                _dOperator[ig1][ig2](alpha, beta)(gamma, delta) -= o._dOperator[ig1][ig2](alpha, beta)(gamma, delta);

    return *this;
  }

  //_________________________________________________________________________
  DoubleOperator& DoubleOperator::operator *= (DoubleOperator const& o)
  {
    // Fast method to check that we are using the same grids
    if (&_grid1 != &o.GetFirstGrid() || &_grid2 != &o.GetSecondGrid())
      throw std::runtime_error(error("DoubleOperator::operator *=", "Grids do not match."));

    const std::vector<std::vector<matrix<matrix<double>>>> v = _dOperator;
    for (int ig1 = 0; ig1 < (int) v.size(); ig1++)
      for (int ig2 = 0; ig2 < (int) v[ig1].size(); ig2++)
        {
          // Set operator entries to zero
          _dOperator[ig1][ig2].set(0);

          // The product between the operators is done ussuming that
          // the first dimension of both pairs of indices is one,
          // which allows us to exploit the symmetry of the operators
          // deriving from the logarithmic distribution of the nodes.
          for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(1); beta++)
            for (int delta = 0; delta < (int) _dOperator[ig1][ig2](0, beta).size(1); delta++)
              for (int alpha = 0; alpha <= beta; alpha++)
                for (int gamma = 0; gamma <= delta; gamma++)
                  _dOperator[ig1][ig2](0, beta)(0, delta) += v[ig1][ig2](0, alpha)(0, gamma) * o._dOperator[ig1][ig2](0, beta - alpha)(0, delta - gamma);
        }
    return *this;
  }

  //_________________________________________________________________________
  DoubleOperator& DoubleOperator::operator = (DoubleOperator const& o)
  {
    if (this != &o)
      *this = o;

    return *this;
  }

  //_________________________________________________________________________
  DoubleOperator& DoubleOperator::operator *= (std::function<double(double const&, double const&)> f)
  {
    for (int ig1 = 0; ig1 < (int) _dOperator.size(); ig1++)
      {
        const std::vector<double>& sg1 = _grid1.GetSubGrid(ig1).GetGrid();
        for (int ig2 = 0; ig2 < (int) _dOperator[ig1].size(); ig2++)
          {
            const std::vector<double>& sg2 = _grid2.GetSubGrid(ig2).GetGrid();
            for (int alpha = 0; alpha < (int) _dOperator[ig1][ig2].size(0); alpha++)
              for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(1); beta++)
                for (int gamma = 0; gamma < (int) _dOperator[ig1][ig2](alpha, beta).size(0); gamma++)
                  for (int delta = 0; delta < (int) _dOperator[ig1][ig2](alpha, beta).size(1); delta++)
                    _dOperator[ig1][ig2](alpha, beta)(gamma, delta) *= f(sg1[beta], sg2[delta]);
          }
      }
    return *this;
  }

  //_________________________________________________________________________
  DoubleOperator& DoubleOperator::operator *= (std::function<double(double const&)> f)
  {
    for (int ig1 = 0; ig1 < (int) _dOperator.size(); ig1++)
      {
        const std::vector<double>& sg1 = _grid1.GetSubGrid(ig1).GetGrid();
        for (int ig2 = 0; ig2 < (int) _dOperator[ig1].size(); ig2++)
          {
            const std::vector<double>& sg2 = _grid2.GetSubGrid(ig2).GetGrid();
            for (int alpha = 0; alpha < (int) _dOperator[ig1][ig2].size(0); alpha++)
              for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(1); beta++)
                for (int gamma = 0; gamma < (int) _dOperator[ig1][ig2](alpha, beta).size(0); gamma++)
                  for (int delta = 0; delta < (int) _dOperator[ig1][ig2](alpha, beta).size(1); delta++)
                    _dOperator[ig1][ig2](alpha, beta)(gamma, delta) *= f(sg1[beta]) * f(sg2[delta]);
          }
      }
    return *this;
  }

  //_________________________________________________________________________
  DistributionOperator DoubleOperator::MultiplyFirstBy(Distribution const& d) const
  {
    // Fast method to check that we are using the same grid for the
    // first variable
    if (&_grid1 != &d.GetGrid())
      throw std::runtime_error(error("DoubleOperator::MultiplyFirstBy", "first grid of the DoubleOperator and and grid of the Distribution do not match."));

    // Get numbers of subgrids
    const int ng1 = _grid1.nGrids();
    const int ng2 = _grid2.nGrids();

    // Get distributions on the subgrids
    const std::vector<std::vector<double>>& sg = d.GetDistributionSubGrid();

    // Initialise output vector
    std::vector<std::vector<std::vector<matrix<double>>>> s(ng1, std::vector<std::vector<matrix<double>>>(ng2));

    // Compute the the distribution on the subgrids using the joint
    // grids.
    for (int ig1 = 0; ig1 < ng1; ig1++)
      {
        // Number of nodes of the ig1-th subgrid
        const int nx1 = _grid1.GetSubGrid(ig1).nx();
        for (int ig2 = 0; ig2 < ng2; ig2++)
          {
            // Number of nodes of the ig1-th subgrid
            const int nx2 = _grid2.GetSubGrid(ig2).nx();

            // Resize container w.r.t the first grid
            s[ig1][ig2].resize(nx1);

            for (int alpha = 0; alpha < nx1; alpha++)
              {
                // Resize container w.r.t the second grid
                s[ig1][ig2][alpha].resize(0, nx2);
                for (int delta = 0; delta < nx2; delta++)
                  // Compute product between first dimension of the
                  // double operator and distribution.
                  for (int beta = 0; beta < nx1 - alpha; beta++)
                    s[ig1][ig2][alpha](0, delta) += _dOperator[ig1][ig2](0, beta)(0, delta) * sg[ig1][beta + alpha];
              }
          }
      }

    // Return DistributionOperator object
    return DistributionOperator{_grid1, _grid2, s};
  }

  //_________________________________________________________________________
  OperatorDistribution DoubleOperator::MultiplySecondBy(Distribution const& d) const
  {
    // Fast method to check that we are using the same grid for the
    // first variable
    if (&_grid2 != &d.GetGrid())
      throw std::runtime_error(error("DoubleOperator::MultiplySecondBy", "second grid of the DoubleOperator and and grid of the Distribution do not match."));

    // Get numbers of subgrids
    const int ng1 = _grid1.nGrids();
    const int ng2 = _grid2.nGrids();

    // Get distributions on the subgrids
    const std::vector<std::vector<double>>& sg = d.GetDistributionSubGrid();

    // Initialise output vector
    std::vector<std::vector<matrix<std::vector<double>>>> s(ng1, std::vector<matrix<std::vector<double>>>(ng2));

    // Compute the the distribution on the subgrids using the joint
    // grids.
    for (int ig1 = 0; ig1 < ng1; ig1++)
      {
        // Number of nodes of the ig1-th subgrid
        const int nx1 = _grid1.GetSubGrid(ig1).nx();
        for (int ig2 = 0; ig2 < ng2; ig2++)
          {
            // Number of nodes of the ig1-th subgrid
            const int nx2 = _grid2.GetSubGrid(ig2).nx();

            // Resize container w.r.t the firt and second grids
            s[ig1][ig2].resize(0, nx1, std::vector<double>(nx2));

            for (int beta = 0; beta < nx1; beta++)
              for (int gamma = 0; gamma < nx2; gamma++)
                // Compute product between second dimension of the
                // double operator and distribution
                for (int delta = 0; delta < nx2 - gamma; delta++)
                  s[ig1][ig2](0, beta)[gamma] += _dOperator[ig1][ig2](0, beta)(0, delta) * sg[ig2][delta + gamma];
          }
      }

    // Return OperatorDistribution object
    return OperatorDistribution{_grid1, _grid2, s};
  }

  //_________________________________________________________________________
  DoubleDistribution operator * (DoubleOperator const& lhs, DoubleDistribution const& rhs)
  {
    return lhs *= rhs;
  }

  //_________________________________________________________________________
  DistributionOperator operator * (Distribution const& lhs, DoubleOperator const& rhs)
  {
    return rhs.MultiplyFirstBy(lhs);
  }
  //_________________________________________________________________________
  OperatorDistribution operator * (DoubleOperator const& lhs, Distribution const& rhs)
  {
    return lhs.MultiplySecondBy(rhs);
  }

  //_________________________________________________________________________
  DoubleOperator operator * (DoubleOperator lhs, DoubleOperator const& rhs)
  {
    return lhs *= rhs;
  }

  //_________________________________________________________________________
  DoubleOperator operator * (double const& s, DoubleOperator rhs)
  {
    return rhs *= s;
  }

  //_________________________________________________________________________
  DoubleOperator operator * (DoubleOperator lhs, double const& s)
  {
    return lhs *= s;
  }

  //_________________________________________________________________________
  DoubleOperator operator * (std::function<double(double const&, double const&)> f, DoubleOperator rhs)
  {
    return rhs *= f;
  }

  //_________________________________________________________________________
  DoubleOperator operator * (DoubleOperator lhs, std::function<double(double const&, double const&)> f)
  {
    return lhs *= f;
  }

  //_________________________________________________________________________
  DoubleOperator operator * (std::function<double(double const&)> f, DoubleOperator rhs)
  {
    return rhs *= f;
  }

  //_________________________________________________________________________
  DoubleOperator operator * (DoubleOperator lhs, std::function<double(double const&)> f)
  {
    return lhs *= f;
  }

  //_________________________________________________________________________
  DoubleOperator operator / (DoubleOperator lhs, double const& s)
  {
    return lhs /= s;
  }

  //_________________________________________________________________________
  DoubleOperator operator + (DoubleOperator lhs, DoubleOperator const& rhs)
  {
    return lhs += rhs;
  }

  //_________________________________________________________________________
  DoubleOperator operator - (DoubleOperator lhs, DoubleOperator const& rhs)
  {
    return lhs -= rhs;
  }

  //_________________________________________________________________________________
  std::ostream& operator << (std::ostream& os, DoubleOperator const& dop)
  {
    const std::vector<std::vector<matrix<matrix<double>>>> om = dop.GetDoubleOperator();
    os << "Operator: " << &dop << "\n";
    os << "Operator on the SubGrids:" << "\n";
    const std::ostringstream default_format;
    os << std::scientific;
    os.precision(2);
    for (int i = 0; i < (int) om.size(); i++)
      for (int j = 0; j < (int) om[i].size(); j++)
        {
          os << "O[" << i << "][" << j << "]: [";
          for (int alpha = 0; alpha < (int) om[i][j].size(0); alpha++)
            for (int beta = 0; beta < (int) om[i][j].size(1); beta++)
              for (int gamma = 0; gamma < (int) om[i][j](alpha, beta).size(0); gamma++)
                for (int delta = 0; delta < (int) om[i][j](alpha, beta).size(1); delta++)
                  os << "{(" << alpha << ", " << beta << ")(" << gamma << ", " << delta << ") : " << om[i][j](alpha, beta)(gamma, delta) << "} ";
          os << "]\n";
        }
    os.copyfmt(default_format);
    return os;
  }
}
