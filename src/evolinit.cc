/*
  evolinit.cc:

  Author: Valerio Bertone
*/

#include <iostream>
#include <stdexcept>
#include <ctime>

#include "apfel/evolinit.h"
#include "apfel/splittings.h"

using namespace std;

namespace apfel {

  // ================================================================================
  // Initializer function
  void evolinit::Initializer(evolsetup Setup_)
  {
    // Timer
    std::clock_t start = clock();

    // Create subgrids
    for(int ig=0; ig<Setup_.nGrids(); ig++) {
      subgrid sg(Setup_.nx(ig), Setup_.xMin(ig), Setup_.InterDegree(ig));
      _GlobalGrid.AddGrid(sg);
    }

    // Lock subgrids if required
    if(Setup_.Locked()) _GlobalGrid.LockGrids();

    // Create global x-space grid
    _GlobalGrid.CreateGrid();

    // Allocate QCD splitting functions
    _sfQCD = new QCD_space_unpol(Setup_.GaussPoints(), _GlobalGrid);

    // Create QCD operators
    for(int imass=3; imass<=6; imass++) _sfQCD->CreateOperators(imass);

    // Print time elapsed
    double dt = ( clock() - start ) / (double) CLOCKS_PER_SEC;
    cout << "Time elapsed: " << dt << " s" << endl;
  }

}
