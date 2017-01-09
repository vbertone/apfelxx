//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <iostream>
#include <cmath>
#include <iomanip>
#include <apfel/qgrid.h>

using namespace apfel;
using namespace std;

int main()
{
  // Constructor tests
  QGrid qg(50, 0.5, 1000, 3, {0, 0, 0, sqrt(2), 4.5, 175.});

//  cout << qg << endl;

  return 0;
}
