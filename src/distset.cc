/*
  distset.cc:

  Author: Valerio Bertone
*/

#include <cmath>
#include <stdexcept>

#include "apfel/distset.h"

using namespace std;

namespace apfel {

  // ================================================================================
  // Constructor
  distset::distset(grid   const& GlobalGrid_,
		   double const& Scale_,
		   int    const& NDistributions_,
		   void   (*DistributionFunction_)(const double& x, const double& Scale, double* xfx)):
    _GlobalGrid(GlobalGrid_),
    _DistributionFunction(DistributionFunction_)
  {
    // Push back scale
    _Scale.push_back(Scale_);

    // Push back number of distributions
    _NDistributions.push_back(NDistributions_);

    // Number of SubGrids
    int ng = _GlobalGrid.nGrids();

    // Check that the grid has been initialized
    if(ng == -1) throw runtime_error("The global grid has not been initialized.");

    // Allocate distributions
    _Distributions.push_back(new double**[NDistributions_]);
    for(int ifl=0; ifl<NDistributions_; ifl++) {
      _Distributions[0][ifl] = new double*[ng];
      for(int ig=0; ig<ng; ig++) {
        int nx = _GlobalGrid.GetSubGrid(ig).nx();
	_Distributions[0][ifl][ig] = new double[nx];
	for(int alpha=0; alpha<nx; alpha++) _Distributions[0][ifl][ig][alpha] = 0;
      }
    }

    // Loop over the SubGrids
    for(int ig=0; ig<ng; ig++) {
      // Number of grid points
      int nx = _GlobalGrid.GetSubGrid(ig).nx();
      // Loop over the grid points
      for(int alpha=0; alpha<nx; alpha++) {
        double x = _GlobalGrid.GetSubGrid(ig).xg(alpha);
	double xfx[NDistributions_];
	DistributionFunction_(x, Scale_, xfx);
	for(int ifl=0; ifl<NDistributions_; ifl++) _Distributions[0][ifl][ig][alpha] = xfx[ifl];
      }
    }
  }

  // ================================================================================
  // Destructor
  distset::~distset()
  {
    int ng = _GlobalGrid.nGrids();
    for(int id=0; id<_Distributions.size(); id++) {
      for(int ifl=0; ifl<_NDistributions[id]; ifl++) {
	for(int ig=0; ig<ng; ig++) delete _Distributions[0][ifl][ig];
	delete _Distributions[0][ifl];
      }
      delete _Distributions[0];
    }
    _Scale.clear();
    _NDistributions.clear();
    _Distributions.clear();
  }

  // ================================================================================
  // Include a new distribution
  distset distset::operator+=(distset const& dist)
  {
    // Check that **dist** is defined on the same global grid.
    if( !(dist.GlobalGrid() == _GlobalGrid) )
      throw runtime_error("To add a distset object it must be defined of the same global grid.");

    for(int id=0; id<dist.Distributions().size(); id++) {
      // Push back scale
      _Scale.push_back(dist.Scale()[id]);

      // Push back number of distributions
      _NDistributions.push_back(dist.NDistributions()[id]);

      // Push back distributions
      _Distributions.push_back(dist.Distributions()[id]);
    }
  }

  // ================================================================================
  // Multiply distribution for and operator
  distset::distset(opset const* op, distset const* dist)
  {
    // Check that **op** and **dist** are defined on the same global grid.
    if( !(op->GlobalGrid() == dist->GlobalGrid()) )
      throw runtime_error("To multiply a opset object with a distset object they must be defined of the same global grid.");
    _GlobalGrid = dist->GlobalGrid();

    // Number of members of the operator
    int NMembers = op->NumberOfMembers();

    // Number of mass scales of the operator. This is given by the total number
    // of operators divided by the number of members.
    int NScale = op->Operators().size() / NMembers;

    // Number of SubGrids
    int ng = _GlobalGrid.nGrids();
    /*
    // Loop over the distribution scales
    for(int iscd=0; iscd<dist->Scale().size(); iscd++) {

      // Number of distributions at isc
      int NDist = dist->NDistributions()[iscd];

      // Loop over the operator mass scales
      for(int isco=0; isco<NScales; isco++) {

	// Push back scale
	_Scale.push_back(Scale()[iscd]);

	// Push back number of distributions
	_NDistributions.push_back(NDistributions()[iscd]);

	// Number of final distribution scales
	int iscf = _Scale.size();

	// Push back distribution
	_Distributions.push_back(new double**[NDist]);

	// Loop over the distributions
	for(int idist=0; idist<NDist; idist++) {

	  _Distributions[iscf][idist] = new double*[ng];

	  // Loop over the SubGrids
	  for(int ig=0; ig<ng; ig++) {

	    int nx = _GlobalGrid.SubGrid(ig).nx();
	    _Distributions[iscf][idist][ig] = new double[nx];

	    // Loop over the nodes of the SubGrid
	    for(int alpha=0; alpha<nx; alpha++) {

	  // Loop over the map to construct the derived distributions
	  for(int imap=0; imap<op->Map(idist).size(); imap++) {

	    int jm = NMembers * isc + imap;



	      _Distributions[iscf][idist][ig][alpha] += op->Operators()[jm][ig][alpha][gamma] dist->Distributions()[iscd][idist][ig][gamma];
		}



	  }
	}
      }
    }

    */
  }

}
