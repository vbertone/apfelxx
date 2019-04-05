//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

#include <iostream>
#include <iomanip>

int main()
{
  std::cout << std::setprecision(15) << std::scientific;

  // Allocate vector of subgrids
  std::vector<apfel::SubGrid> sgs = {apfel::SubGrid{10,1e-5,3}, apfel::SubGrid{20,1e-3,2}, apfel::SubGrid{30,1e-1,2}};

  // Construct grid with and without locking the subgrids
  apfel::Grid gl{sgs, true};
  apfel::Grid g{sgs, false};

  // Print grids
  std::cout << g << std::endl;
  std::cout << gl << std::endl;

  if (g == gl)
    return 0;

  return 0;
}
