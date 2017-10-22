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
#include <apfel/constants.h>

#include <cmath>

using namespace apfel;
using namespace std;

// Class to define the analytical expression of LO splitting function P0qq
class p0qq: public Expression
{
public:
  p0qq(): Expression() {}
  double Regular(double const& x)  const { return - 2 * CF * ( 1 + x ); }
  double Singular(double const& x) const { return 4 * CF / ( 1 - x ); }
  double Local(double const& x)    const { return 4 * CF * log( 1 - x ) + 3 * CF; }
};

// Class to define the analytical expression of the squared LO splitting function P0qq
class p0qq2: public Expression
{
public:
  p0qq2(): Expression() {}
  double Regular(double const& x)  const { return 4 * CF * CF * ( - 4 * log(x) / ( 1 - x ) - 4 * ( 1 + x ) * log( 1 - x ) + 3 * ( 1 + x ) * log(x) - ( x + 5 ) ); }
  double Singular(double const& x) const { return 4 * CF * CF * ( 8 * log( 1 - x ) + 6 ) / ( 1 - x ); }
  double Local(double const& x)    const { return 4 * CF * CF * ( 4 * pow(log( 1 - x ),2) + 6 * log( 1 - x ) + 9. / 4. - 4 * pow(M_PI,2) / 6. ) ;}
};

int main()
{
  Timer t;

  // Grid
  const Grid g{{SubGrid{80,1e-5,3}, SubGrid{50,1e-1,3}, SubGrid{40,8e-1,3}}, false};

  // Distribution
  const Distribution d{g, [&] (double const& x)->double{ return 1 - x; }};

  // Expression
  const p0qq p;

  // Operator
  cout << "\nInitialization ..." << endl;
  t.start();
  const Operator O{g, p};
  t.stop();

  // Multiply operator by the distribution to create a new distribution
  cout << "\nConvolution between operator and distribution (O * d) ..." << endl;
  t.start();
  auto Od = O * d;
  t.stop();

  // Multiply operator by itself to create a new operator
  cout << "\nConvolution between two operators (O * O) ..." << endl;
  t.start();
  auto OO = O * O;
  t.stop();

  // Check the numerical accuracy of "Od" by comparing with the analytical result
  cout << "\nChecking the numerical accuracy of O * d ... " << endl;
  cout << scientific;
  for (auto ix = 0; ix <= g.GetJointGrid().nx(); ix++)
    {
      double x = g.GetJointGrid().GetGrid()[ix];
      // Analytic result for x \int_x^1 dy Pqq(y) ( 1 - x / y )
      double Ix = CF * ( - 2 * ( 3. / 2. - x - pow(x,2) / 2. ) + 4 * ( 1 - x ) * log( 1 -  x ) + 3 * ( 1 - x ) + 2 * x * ( log(x) + 1 - x ) );
      cout << x << "\t\t" << Od.Evaluate(x) << "\t\t" << Ix << "\t\t" << Od.Evaluate(x) / Ix << endl;
    }

  // Check the numerical accuracy of "Od" by comparing with the analytical result
  // Analytical expression of P0qq \otimes P0qq
  const p0qq2 p2;
  const Operator O2{g, p2};

  // Now convolute both "OO" and "O2" with the test distribution "d" and check the result
  auto OOd = OO * d;
  auto O2d = O2 * d;

  cout << "\nChecking the numerical accuracy of O * O ... " << endl;
  cout << scientific;
  for (auto ix = 0; ix <= g.GetJointGrid().nx(); ix++)
    {
      double x = g.GetJointGrid().GetGrid()[ix];
      cout << x << "\t\t" << OOd.Evaluate(x) << "\t\t" << O2d.Evaluate(x) << "\t\t" << OOd.Evaluate(x) / O2d.Evaluate(x) << endl;
    }

  return 0;
}
