//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <iostream>
#include <cmath>
#include <apfel/integrator.h>
using namespace apfel;
using namespace std;

class MyFunction: public Integrator
{
public:
  MyFunction(): Integrator() {}
  double integrand(double const& x) const
  {
    return log(x);
  }

};

int main()
{ 
  const MyFunction f;
  // Integrate using dgauss with a given accuracy
  const double res1 = f.integrate(0, 2, 1e-5);
  // integrate using dgauss with a given number of nodes
  const double res2 = f.integrate(0, 2, 4);

  cout << res1 << "  " << res2 << endl;

  return 0;
}
