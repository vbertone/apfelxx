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
  if(x < 1) return x * ( 1 - x );
  else      return 0;
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
      _distributionJointGrid.push_back(xg(ix));

    for (auto ig=0; ig<_grid.nGrids(); ig++)
      {
        vector<double> sg;
        for (auto const& ix: _grid.GetSubGrid(ig).GetGrid())
          sg.push_back(xg(ix));
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
  Timer t;
  t.start();

  // Grid
  const Grid g{{SubGrid{80,1e-5,3}, SubGrid{50,1e-1,5}, SubGrid{40,8e-1,5}}, false};

  // Expression
  const p0qq p;

  const myPDF d{g};
  const Operator op(g, p);

  auto new_d = op*d;
  auto new_o = op*op;

  t.printTime(t.stop());

  return 0;
}
