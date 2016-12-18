/*
  opset.cc:

  Author: Valerio Bertone
*/

#include <iostream>
#include <stdexcept>

#include <apfel/opset.h>
#include <apfel/dgauss.h>

using namespace std;

namespace apfel {

  // ================================================================================
  // Create the the convolution operator by integrating the operator
  // with interpolation functions.
  void opset::CreateOperators(int const& imass)
  {
    // Get number of subgrids
    int ng = _GlobalGrid.nGrids();

    // Check that the grid has been initialized
    if(ng == -1) throw runtime_error("The global grid has not been initialized.");

    // Loop over the member functions
    for(int im=0; im<NumberOfMembers(); im++) {

      // Loop over the subgrids
      _Operators.push_back(new double**[ng]);
      for(int ig=0; ig<ng; ig++) {
	subgrid sg = _GlobalGrid.SubGrid(ig);

	// Interpolation degree
	int id = sg.InterDegree();

	// Number of grid points
	int nx = sg.nx();

	int gbound = 1;
	if(sg.IsExt()) gbound = nx;

	// Loop over the grid points
	_Operators[im][ig] = new double*[gbound];
	for(int beta=0; beta<gbound; beta++) {

	  double xbeta = sg.xg(beta);
	  _Operators[im][ig][beta] = new double[nx];
	  for(int alpha=beta; alpha<nx; alpha++) {

	    // Delta function
	    double wL = 0;
	    if(alpha == beta) wL = 1;

	    // Integral
	    double I = 0;

	    // Given that the interpolation functions are discontinuos and
	    // wiggly, split the integrals into (id+1) intervals on each of
	    // which the integrand is smooth.

	    // How many grid intervals do we need to integrate over?
	    int nint = min(id+1,alpha);

	    for(int jint=0; jint<nint; jint++) {

	      // Define integration bounds of the first iteration
	      double c = xbeta / sg.xg(alpha-jint+1);
	      double d = xbeta / sg.xg(alpha-jint);

	      // Gauss integration
	      double k1 = ( d - c ) / 2;
	      double k2 = ( d + c ) / 2;
	      int m = pow(2,_GaussPoints+1);
	      for(int i=0; i<m; i++) {

		double xp = k2 + gq_x[_GaussPoints][i] * k1;
		double xm = k2 - gq_x[_GaussPoints][i] * k1;

		double wp = sg.Interpolant(alpha,xbeta/xp);
		double wm = sg.Interpolant(alpha,xbeta/xm);

		// Regular part
		I += OperatorFunction(im,REGULAR,imass,xp) * wp + OperatorFunction(im,REGULAR,imass,xm) * wm;

		// Singular part
		I += OperatorFunction(im,SINGULAR,imass,xp) * ( wp - wL ) + OperatorFunction(im,SINGULAR,imass,xm) * ( wm - wL );

		// Multiply by the weight
		I *= gq_w[_GaussPoints][i];
	      }
	      // Overall normalization
	      I *= k1;
	    }
	    _Operators[im][ig][beta][alpha] = I;
	  }
	  // Include local part
	  _Operators[im][ig][beta][beta] += OperatorFunction(im,LOCAL,imass,xbeta/sg.xg(beta+1));
	}
      }
    }
  }


  // ================================================================================
  // Multiplication constructor
  opset::opset(int const& im1, opset const* op1, int const& im2, opset const* op2)
  {
    // Check that **op1** and **op2** are defined on the same global grid.
    if( !(op1->GlobalGrid() == op2->GlobalGrid()) )
      throw runtime_error("opset objects cannot be added because the respective global grids are different.");
    else
      _GlobalGrid = op1->GlobalGrid();

    // Check that **op1** and **op2** have been computed using the same number of points for the Gauss integration.
    // (This is not strictly necessary but it's safer to add this requirement).
    if( op1->GaussPoints() != op2->GaussPoints() )
      throw runtime_error("opset objects cannot be added because the respective numbers of Gauss points are different.");
    else
      _GaussPoints = op1->GaussPoints();

    // Start the multiplication.
    // Get number of subgrids
    int ng = _GlobalGrid.nGrids();

    // Allocate operators
    _Operators.push_back(new double**[ng]);

    // Loop over the subgrids
    for(int ig=0; ig<ng; ig++) {
      subgrid sg = _GlobalGrid.SubGrid(ig);

      // Number of grid points
      int nx = sg.nx();

      int gbound = 1;
      if(sg.IsExt()) gbound = nx;

      // Loop over the grid points
      _Operators[0][ig] = new double*[gbound];
      for(int beta=0; beta<gbound; beta++) {

	_Operators[0][ig][beta] = new double[nx];
	for(int alpha=beta; alpha<nx; alpha++) {

	  _Operators[0][ig][beta][alpha] = 0;
	  for(int gamma=beta; gamma<alpha; gamma++) {

	    _Operators[0][ig][beta][alpha] +=
	      op1->Operators()[im1][ig][beta][gamma] *
	      op2->Operators()[im2][ig][gamma][alpha];
	  }
	}
      }
    }

  }

  // ================================================================================
  // Destructor
  opset::~opset()
  {
    int ng = _GlobalGrid.nGrids();
    for(int im=0; im<_Operators.size(); im++) {
      for(int ig=0; ig<ng; ig++) {
	subgrid sg = _GlobalGrid.SubGrid(ig);
	int nx = sg.nx();
	int gbound = 1;
	if(sg.IsExt()) gbound = nx;
	for(int beta=0; beta<gbound; beta++) delete _Operators[im][ig][beta];
	delete _Operators[im][ig];
      }
      delete _Operators[im];
    }
    _Operators.clear();
  }


  // ================================================================================
  // Constructor to join two operators
  opset::opset(opset const* op1, opset const* op2)
  {
    // Check that **op1** and **op2** are defined on the same global grid.
    if( !(op1->GlobalGrid() == op2->GlobalGrid()) )
      throw runtime_error("To join two opset objects they must be defined on the same global grid.");

    // Check that **op1** and **op2** have been computed using the same number of points for the Gauss integration.
    // (This is not strictly necessary but it's safer to add this requirement).
    if( op1->GaussPoints() != op2->GaussPoints() )
      throw runtime_error("To join two opset objects they must be computed using the same number of Gauss points.");

    // Push back operators from **op1**
    for(int im=0; im<op1->Operators().size(); im++) _Operators.push_back(op1->Operators()[im]);

    // Push back operators from **op2**
    for(int im=0; im<op2->Operators().size(); im++) _Operators.push_back(op2->Operators()[im]);
  }

}
