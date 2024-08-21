//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

#include <iomanip>

int main()
{
  std::cout << std::setprecision(15) << std::scientific;

  // Allocate vector of subgrids
  std::vector<apfel::SubGrid> sgs = {apfel::SubGrid{100, 1e-5, 3}, apfel::SubGrid{100, 1e-3, 2}, apfel::SubGrid{100, 1e-1, 2}};

  // Construct grid
  apfel::Grid g{sgs};

  // Print grid
  std::cout << g << std::endl;

  return 0;
}
