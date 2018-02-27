//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <apfel/grid.h>
#include <apfel/subgrid.h>

#include <iostream>
#include <iomanip>

using namespace apfel;
using namespace std;

int main()
{
  cout << setprecision(15) << scientific;

  // Allocate vector of subgrids
  vector<SubGrid> sgs = {
    SubGrid{10,1e-5,3},
    SubGrid{20,1e-3,2},
    SubGrid{30,1e-1,2}
  };

  // Construct grid with and without locking the subgrids
  Grid gl{sgs, true};
  Grid g{sgs, false};

  // Print grids
  cout << g << endl;
  cout << gl << endl;

  if (g == gl)
    return 0;

  return 0;
}
