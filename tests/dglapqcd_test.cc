//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <apfel/dglapqcd.h>
#include <apfel/grid.h>
#include <apfel/subgrid.h>
#include <apfel/timer.h>
#include <cmath>
#include <map>
#include <functional>

using namespace apfel;
using namespace std;

// LH Toy PDFs
double xupv(double const& x)  { return 5.107200 * pow(x,0.8) * pow((1-x),3); }
double xdnv(double const& x)  { return 3.064320 * pow(x,0.8) * pow((1-x),4); }
double xglu(double const& x)  { return 1.7 * pow(x,-0.1) * pow((1-x),5); }
double xdbar(double const& x) { return 0.1939875 * pow(x,-0.1) * pow((1-x),6); }
double xubar(double const& x) { return xdbar(x) * (1-x); }
double xsbar(double const& x) { return 0.2 * ( xdbar(x) + xubar(x) ); }
double LHToyPDFs(int const& i, double const& x)
{
  // Gluon
  if      (i == 0  ) return xglu(x);
  // Singlet, T15, T24, T35
  else if (i == 1  ||
	   i == 7  ||
	   i == 9  ||
	   i == 11 ) return xdnv(x) + 2 * xdbar(x) + xupv(x) + 2 * xubar(x) + 2 * xsbar(x);
  // T3
  else if (i == 3  ) return xupv(x) + 2 * xubar(x) - xdnv(x) - 2 * xdbar(x);
  // T8
  else if (i == 5  ) return xupv(x) + 2 * xubar(x) + xdnv(x) + 2 * xdbar(x) - 4 * xsbar(x);
  // Valence, V8, V15, V24, V35
  else if (i == 2  ||
	   i == 6  ||
	   i == 8  ||
	   i == 10 ||
	   i == 12 ) return xupv(x) + xdnv(x);
  // V3
  else if (i == 4  ) return xupv(x) - xdnv(x);
  else               return 0;
}

int main()
{
  // Time counter
  Timer t;
  t.start();

  // Grid
  const Grid g{{SubGrid{80,1e-5,3}, SubGrid{50,1e-1,5}, SubGrid{40,8e-1,5}}};

  // Initialize QCD DGLAP evolution
  const DglapQCD evolution{g, LHToyPDFs, sqrt(2), 0.35, sqrt(2), {0,0,0,sqrt(2),4.5,175}, {0,0,0,sqrt(2),4.5,175}, 0};

  t.printTime(t.stop());

  auto PDFs = evolution.GetObject(100);

  //cout << PDFs.GetObjects().size() << endl;

  //for ( auto it = PDFs.GetObjects().begin(); it != PDFs.GetObjects().end(); ++it )
  //  std::cout << " " << it->first << ":" << endl;


  //cout << PDFs.at(0).Evaluate(0.0001) << endl;

  return 0;
}
