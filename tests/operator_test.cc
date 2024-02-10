//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

// Test expression (the LO Pqq splitting function)
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

int main()
{
  // Grid
  const apfel::Grid g{{apfel::SubGrid{80, 1e-5, 3}, apfel::SubGrid{50, 1e-1, 5}, apfel::SubGrid{40, 8e-1, 5}}};

  // Expression
  const p0qq p;

  apfel::Timer t;
  // Construct the operator
  const apfel::Operator op{g, p};
  t.stop();

  // Print operator
  std::cout << op << std::endl;;

  return 0;
}
