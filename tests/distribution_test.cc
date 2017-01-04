//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <apfel/distribution.h>
#include <apfel/grid.h>
#include <apfel/subgrid.h>
#include <apfel/operator.h>
#include <apfel/expression.h>
#include <apfel/timer.h>
#include <apfel/tools.h>
#include <cmath>

using namespace apfel;
using namespace std;

/**
 * Prototype function for testing purproses.
 */
double xg(double const& x)
{
  return 1 - x;
}

/**
 * @brief The Parton class
 */
class myPDF: public Distribution
{
public:
  /**
   * Allocate the langrage interpolation and fill the inherited
   * \c _distribution object with the jointed grid.
   */
  myPDF(Grid const& gr): Distribution(gr)
  {
    for (auto const& ix: _grid.GetJointGrid().GetGrid())
      if (ix < 1) _distributionJointGrid.push_back(xg(ix));
      else        _distributionJointGrid.push_back(0);

    for (auto ig=0; ig<_grid.nGrids(); ig++)
      {
        vector<double> sg;
        for (auto const& ix: _grid.GetSubGrid(ig).GetGrid())
          if (ix < 1) sg.push_back(xg(ix));
	  else        sg.push_back(0);
        _distributionSubGrid.push_back(sg);
      }
  }
};

class p0qq: public Expression
{
public:
  p0qq(): Expression() {}
  double Regular(double const& x)  const { return - 2 * CF * ( 1 + x ); }
  double Singular(double const& x) const { return 4 * CF / ( 1 - x ); }
  double Local(double const& x)    const { return 4 * CF * log( 1 - x ) + 3 * CF; }
};


int main()
{
  // Grid
  const Grid g{{SubGrid{80,1e-5,3}, SubGrid{50,1e-1,3}, SubGrid{40,8e-1,5}}};

  // Distribution
  const myPDF d{g};

  // Expression
  const p0qq p;

  // Operator
  const Operator op{g, p};

  // Multiply operator by the distribution to create a new distribution
  auto new_d = op * d;

  // Check the numerical accuracy of "new_d" by comparing with the analytical result
  cout << scientific;
  for (auto ix = 0; ix <= g.GetJointGrid().nx(); ix++)
    {
      double x = g.GetJointGrid().GetGrid()[ix];
      // Analytic result for x \int_x^1 dy Pqq(y) ( 1 - x / y )
      double Ix = CF * ( - 2 * ( 3. / 2. - x - pow(x,2) / 2. ) + 4 * ( 1 - x ) * log( 1 -  x ) + 3 * ( 1 - x ) + 2 * x * ( log(x) + 1 - x ) );
      cout << x << "\t\t" << new_d.Evaluate(x) << "\t\t" << Ix << "\t\t" << new_d.Evaluate(x) / Ix << endl;
    }

  // Multiply operator by itself to create a new operator
  auto new_o = op * op;

  return 0;
}
