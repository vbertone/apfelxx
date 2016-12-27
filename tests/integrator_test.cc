//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <iostream>
#include <apfel/integrator.h>
using namespace apfel;
using namespace std;

class MyFunction: public Integrator
{
public:
  MyFunction(): Integrator() {}
  double integrand(double const& x) const
  {
    return x*x*x;
  }

};

int main()
{ 
  const MyFunction f;
  const double res = f.integrate(0, 1);

  if (res != 1/3.)
    return -1;

  return 0;
}
