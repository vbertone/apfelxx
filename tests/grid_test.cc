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
using namespace apfel;
using namespace std;

int main()
{
  cout << setprecision(15) << scientific;

  vector<SubGrid> sgs = {
    SubGrid{10,1e-5,3},
    SubGrid{20,1e-3,2},
    SubGrid{30,1e-1,2}
  };

  Grid g{sgs, false};
  Grid gl{sgs, true};

  cout << g << endl;
  cout << gl << endl;

  if (g == gl) return false;

  return 0;
}
