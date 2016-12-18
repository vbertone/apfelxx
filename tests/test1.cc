#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <iomanip>

#include <apfel/subgrid.h>
#include <apfel/grid.h>
#include <apfel/evolsetup.h>
#include <apfel/evolinit.h>
#include <apfel/distset.h>
#include <apfel/opset.h>
#include <apfel/splittings.h>
#include <apfel/utils.h>

using namespace std;

using namespace apfel;

// Test functions
double func0(double const& x) { return 1;}
double func1(double const& x) { return x;}
double func2(double const& x) { return x * x;}
double func3(double const& x) { return exp(x);}

void LHToyPDFsPhys(double const& x, double const& Q, double* physxpdf);
void LHToyPDFsEvol(double const& x, double const& Q, double* evolxpdf)
{
  double physxpdf[14];
  LHToyPDFsPhys(x, Q, physxpdf);
  Phys2Evol(physxpdf, evolxpdf);
}

int main() {

  // Initialize 3 different setups: setup, setup1, setup2
  evolsetup setup;
  setup.SetQLimits(2,100);
  setup.SetGaussPoints(0);
  evolsetup setup1;
  evolsetup setup2 = setup;

  // Print the number of Gauss points for the three setups
  cout << setup.GaussPoints() << endl;
  cout << setup1.GaussPoints() << endl;
  cout << setup2.GaussPoints() << endl;

  // Initialize the evolutions
  evolinit initdefault;         // Default settings
  evolinit initdefault1;        // Default settings
  evolinit initcustom2(setup2); // Settings in setup2
  evolinit initcustom1(setup1); // Settings in setup1
  evolinit initcustom(setup);   // Settings in setup

  // Retrieve global grid from "initcustom"
  grid ggrid = initcustom.GlobalGrid();

  // Define a set of "distribution" objects of the "ggrid" grid
  distset *testdist = new distset(ggrid, 1.4142, 14, LHToyPDFsEvol);
  cout << testdist->Scale()[0] << endl;

  // Define a set of "operators" objects on the "ggrid" grid
  opset *testop = new QCD_space_unpol(3, ggrid);
  for(int imass=3; imass<=6; imass++) testop->CreateOperators(imass);
  cout << testop->Map(2).size() << " " << testop->Map(1).size() << " " << testop->Map(0).size() << " " << testop->Map(3).size() << endl;
  cout << testop->Map(2)[1] << " " << testop->Map(1)[1] << " " << testop->Map(3)[0] << endl;
  cout << testop->Operators().size() << "  " << testop->NumberOfMembers() << endl;
  cout << testop->Operators().size() / testop->NumberOfMembers() << endl;

}

void LHToyPDFsPhys(double const& x, double const& Q, double* xfx)
{
  // Parameters
  double N_uv = 5.107200;
  double auv  = 0.8;
  double buv  = 3;
  double N_dv = 3.064320;
  double adv  = 0.8;
  double bdv  = 4;
  double N_g  = 1.7;
  double ag   = -0.1;
  double bg   = 5;
  double N_db = 0.1939875;
  double adb  = -0.1;
  double bdb  = 6;
  double fs   = 0.2;

  // Construct PDFs
  double xuv   = N_uv * pow(x,auv) * pow(1-x,buv);
  double xdv   = N_dv * pow(x,adv) * pow(1-x,bdv);
  double xg    = N_g  * pow(x,ag)  * pow(1-x,bg);
  double xdbar = N_db * pow(x,adb) * pow(1-x,bdb);
  double xubar = xdbar * ( 1 - x );
  double xs    = fs * ( xdbar + xubar );
  double xsbar = xs;

  xfx[0]  = 0;
  xfx[1]  = 0;
  xfx[2]  = 0;
  xfx[3]  = xs;
  xfx[4]  = xuv + xubar;
  xfx[5]  = xdv + xdbar;
  xfx[6]  = xg;
  xfx[7]  = xdbar;
  xfx[8]  = xubar;
  xfx[9]  = xsbar;
  xfx[10] = 0;
  xfx[11] = 0;
  xfx[12] = 0;
  xfx[13] = 0;
}
