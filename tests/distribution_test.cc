//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

#include <cmath>

// Class to define the analytical expression of LO splitting function P0qq
class p0qq: public apfel::Expression
{
public:
  p0qq(): Expression() {}
  double Regular(double const& x)  const
  {
    return - 2 * apfel::CF * ( 1 + x );
  }
  double Singular(double const& x) const
  {
    return 4 * apfel::CF / ( 1 - x );
  }
  double Local(double const& x)    const
  {
    return 4 * apfel::CF * log( 1 - x ) + 3 * apfel::CF;
  }
};

// Class to define the analytical expression of the squared LO splitting function P0qq
class p0qq2: public apfel::Expression
{
public:
  p0qq2(): Expression() {}
  double Regular(double const& x)  const
  {
    return 4 * apfel::CF * apfel::CF * ( - 4 * log(x) / ( 1 - x ) - 4 * ( 1 + x ) * log( 1 - x ) + 3 * ( 1 + x ) * log(x) - ( x + 5 ) );
  }
  double Singular(double const& x) const
  {
    return 4 * apfel::CF * apfel::CF * ( 8 * log( 1 - x ) + 6 ) / ( 1 - x );
  }
  double Local(double const& x)    const
  {
    return 4 * apfel::CF * apfel::CF * ( 4 * pow(log( 1 - x ),2) + 6 * log( 1 - x ) + 9. / 4. - 4 * pow(M_PI,2) / 6. ) ;
  }
};

int main()
{
  apfel::Timer t;

  // Grid
  const apfel::Grid g{{apfel::SubGrid{80,1e-5,3}, apfel::SubGrid{50,1e-1,3}, apfel::SubGrid{40,8e-1,3}}};

  // Distribution
  const apfel::Distribution d{g, [&] (double const& x) -> double{ return 1 - x; }};

  // Print distribution
  std::cout << d;

  // Expression
  const p0qq p;

  // Operator
  std::cout << "\nInitialization ..." << std::endl;
  t.start();
  const apfel::Operator O{g, p};
  t.stop();

  // Multiply operator by the distribution to create a new distribution
  std::cout << "\nConvolution between operator and distribution (O * d) ..." << std::endl;
  t.start();
  const apfel::Distribution Od = O * d;
  t.stop();

  // Multiply operator by itself to create a new operator
  std::cout << "\nConvolution between two operators (O * O) ..." << std::endl;
  t.start();
  const apfel::Operator OO = O * O;
  t.stop();

  // Check the numerical accuracy of "Od" by comparing with the analytical result
  std::cout << "\nChecking the numerical accuracy of O * d ... " << std::endl;
  std::cout << std::scientific;
  for (int ix = 0; ix <= g.GetJointGrid().nx(); ix++)
    {
      double x = g.GetJointGrid().GetGrid()[ix];
      // Analytic result for x \int_x^1 dy Pqq(y) ( 1 - x / y )
      double Ix = apfel::CF * ( - 2 * ( 3. / 2. - x - pow(x,2) / 2. ) + 4 * ( 1 - x ) * log( 1 -  x ) + 3 * ( 1 - x ) + 2 * x * ( log(x) + 1 - x ) );
      std::cout << x << "\t\t" << Od.Evaluate(x) << "\t\t" << Ix << "\t\t" << Od.Evaluate(x) / Ix << std::endl;
    }

  // Check the numerical accuracy of "Od" by comparing with the analytical result
  // Analytical expression of P0qq \otimes P0qq
  const p0qq2 p2;
  const apfel::Operator O2{g, p2};

  // Now convolute both "OO" and "O2" with the test distribution "d" and check the result
  const apfel::Distribution OOd = OO * d;
  const apfel::Distribution O2d = O2 * d;

  std::cout << "\nChecking the numerical accuracy of O * O ... " << std::endl;
  for (int ix = 0; ix <= g.GetJointGrid().nx(); ix++)
    {
      double x = g.GetJointGrid().GetGrid()[ix];
      std::cout << x << "\t\t" << OOd.Evaluate(x) << "\t\t" << O2d.Evaluate(x) << "\t\t" << OOd.Evaluate(x) / O2d.Evaluate(x) << std::endl;
    }

  // Now define a double object with O and d and print it.
  const apfel::DoubleObject<apfel::Operator, apfel::Distribution> DObj{{{1, O, d}}};
  std::cout << DObj << std::endl;

  return 0;
}
