/*
 * APFEL++ 2017
 *
 * Authors: Valerio Bertone: valerio.bertone@cern.ch
 *          Stefano Carrazza: stefano.carrazza@cern.ch
 */

#include <iostream>
#include <iomanip>
#include <apfel/grid.h>
#include <apfel/subgrid.h>
#include <apfel/interpolator.h>
using namespace apfel;
using namespace std;

/**
 * Prototype function for testing purproses.
 */
double xg(double const& x)
{
  return x*(1-x);
}

/**
 * @brief The Parton class
 */
class Parton: public LagrangeInterpolator
{
public:
  /**
   * Allocate the langrage interpolation and fill the inherited
   * \c _distribution object with the jointed grid.
   */
  Parton(Grid const& gr): LagrangeInterpolator(gr)
  {
    for (auto const& ix: _grid.GetJointGrid().GetGrid())
      _distribution.push_back(xg(ix));
  }
};

int main()
{
  cout << setprecision(15) << scientific;

  const Grid g{
    {SubGrid{10,1e-5,3}, SubGrid{20,1e-1,3}}
  };

  const Parton xgluon(g);

  vector<double> x = {1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7};
  cout << "x, original function, interpolated function, ratio" << endl;
  for (auto const& ix: x)
    {
      const auto original = xg(ix);
      const auto interpol = xgluon.Evaluate(ix);
      cout << ix << " "
           << original << " "
           << interpol << " "
           << original/interpol<< endl;
    }

  return 0;
}
