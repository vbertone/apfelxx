//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <apfel/constants.h>
#include <apfel/expression.h>
#include <apfel/lagrangeinterpolator.h>
#include <apfel/grid.h>
#include <apfel/subgrid.h>
#include <apfel/operator.h>
#include <apfel/timer.h>

#include <cmath>

using namespace apfel;
using namespace std;

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
  const Grid g{{SubGrid{80,1e-5,3}, SubGrid{50,1e-1,5}, SubGrid{40,8e-1,5}}, false};

  // Expression
  const p0qq p;

  Timer t;
  // Construct the operator
  const Operator op{g, p};
  t.stop();

  return 0;
}
