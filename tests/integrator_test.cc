//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <apfel/integrator.h>

#include <iostream>
#include <cmath>

using namespace apfel;
using namespace std;

int main()
{ 
  const Integrator f{[&] (double const& x)->double{ return log(x); }};
  // Integrate using dgauss with a given accuracy
  const double res1 = f.integrate(0, 2, 1e-5);
  // integrate using dgauss with a given number of nodes
  const double res2 = f.integrate(0, 2, 4);

  cout << res1 << "  " << res2 << "  " << res1 / res2 << endl;

  return 0;
}
