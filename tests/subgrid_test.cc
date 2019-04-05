//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

#include <iostream>
#include <memory>
#include <iomanip>

int main()
{
  std::cout << std::setprecision(12) << std::scientific;

  // Construct test subgrid
  apfel::SubGrid a{50, 1e-5, 3};
  std::cout << a << std::endl;

  // Copy constructor
  apfel::SubGrid b = a;
  std::cout << b << std::endl;

  // Dynamic allocation through copy constructor
  std::unique_ptr<apfel::SubGrid> c{new apfel::SubGrid{a}};
  std::cout << *c << std::endl;

  // Copy assigment
  apfel::SubGrid d{100, 1e-3, 4};
  std::cout << d << std::endl;
  d = a;
  std::cout << d << std::endl;

  // External constructor
  std::vector<double> v = {1e-5, 1e-4, 1e-3, 1e-2, 1};
  apfel::SubGrid e{v, 3};
  std::cout << e << std::endl;

  // Check equality operator
  if ( (a == b)  == false ||
       (a == *c) == false ||
       (a == d)  == false ||
       (a == e)  == true )
    return -1;

  return 0;
}
