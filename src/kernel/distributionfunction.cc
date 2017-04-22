//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/distributionfunction.h"

#include <functional>

using std::function;

namespace apfel {

  //_________________________________________________________________________
  DistributionFunction::DistributionFunction(Grid                                        const& g,
					     function<double(int const&, double const&)> const& InPDFsFunc,
					     int                                         const& ipdf):
    Distribution(g)
  {
    for (auto const& ix: _grid.GetJointGrid().GetGrid())
      if (ix < 1)
	_distributionJointGrid.push_back(InPDFsFunc(ipdf,ix));
      else
	_distributionJointGrid.push_back(0);

    for (auto ig=0; ig<_grid.nGrids(); ig++)
      {
	vector<double> sg;
	for (auto const& ix: _grid.GetSubGrid(ig).GetGrid())
	  if (ix < 1)
	    sg.push_back(InPDFsFunc(ipdf,ix));
	  else
	    sg.push_back(0);
	_distributionSubGrid.push_back(sg);
      }
  }

  //_________________________________________________________________________
  DistributionFunction::DistributionFunction(Grid                                                       const& g,
					     function<double(int const&, double const&, double const&)> const& InPDFsFunc,
					     int                                                        const& ipdf,
					     double                                                     const& Q):
    Distribution(g)
  {
    for (auto const& ix: _grid.GetJointGrid().GetGrid())
      if (ix < 1)
	_distributionJointGrid.push_back(InPDFsFunc(ipdf,ix,Q));
      else
	_distributionJointGrid.push_back(0);

    for (auto ig=0; ig<_grid.nGrids(); ig++)
      {
	vector<double> sg;
	for (auto const& ix: _grid.GetSubGrid(ig).GetGrid())
	  if (ix < 1)
	    sg.push_back(InPDFsFunc(ipdf,ix,Q));
	  else
	    sg.push_back(0);
	_distributionSubGrid.push_back(sg);
      }
  }

}
