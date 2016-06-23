/*
  distset.hh:

  Author: Valerio Bertone
*/

#pragma once

#include <string>

#include "grid.hh"
#include "opset.hh"

using namespace std;

namespace apfel {

  /**
   * class distset:
   * Class for the management of the distributions (PDFs, FFs, atc.).
   **/
  class distset {
  private:
    // Attributes
    grid              _GlobalGrid;
    vector<double>    _Scale;
    vector<int>       _NDistributions;
    void              (*_DistributionFunction)(double const&, double const&, double*);
    vector<double***> _Distributions;

  public:
    // Constructors
    distset(grid   const& GlobalGrid_,
	    double const& Scale_,
	    int    const& NDistributions_,
	    void   (*DistributionFunction_)(const double& x, const double& Scale, double* xfx));

    distset(opset const* op, distset const* dist);

    // Destructor
    ~distset();

    // Getters
    grid GlobalGrid()                 const { return _GlobalGrid; }
    vector<double> Scale()            const { return _Scale; }
    vector<int> NDistributions()      const { return _NDistributions; }
    vector<double***> Distributions() const { return _Distributions; }

    // Include a new distribution
    distset operator+=(distset const& dist);
  };

}
