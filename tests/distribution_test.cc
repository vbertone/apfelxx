//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <apfel/distribution.h>
#include <apfel/timer.h>

using namespace apfel;
using namespace std;

class myPDF: public Distribution
{
public:
  myPDF(): Distribution() {}
  double GetPDF(int const& id, double const& x, double const& q, int const& n = 0) const
  {
    if(x < 1) return x * ( 1 - x );
    else      return 0;
  }
};

int main()
{
  Timer t;
  t.start();

  const myPDF p;

  t.printTime(t.stop());

  return 0;
}
