/*
 * APFEL++ 2017
 *
 * Authors: Valerio Bertone: valerio.bertone@vu.nl
 *          Stefano Carrazza: stefano.carrazza@cern.ch
 */

#include <iostream>
#include <memory>
#include <iomanip>
#include <apfel/subgrid.h>
using namespace apfel;
using namespace std;

int main()
{
  cout << setprecision(15) << scientific;

  // Constructor tests
  SubGrid a{50, 1e-5, 3};
  cout << a << endl;

  // Copy constructor
  SubGrid b = a;
  cout << b << endl;

  // Dynamic allocation through copy constructor
  unique_ptr<SubGrid> c{new SubGrid{a}};
  cout << *c << endl;

  // Copy assigment
  SubGrid d{100, 1e-3, 4};
  cout << d << endl;
  d = a;
  cout << d << endl;

  // External constructor
  vector<double> v = {1e-5, 1e-4, 1e-3, 1e-2, 1};
  SubGrid e{v, 3};
  cout << e << endl;

  // Check equality operator
  if ( (a == b) == false  ||
       (a == *c) == false ||
       (a == d) == false  ||
       (a == e) == true
     )
    return -1;

  return 0;
}
