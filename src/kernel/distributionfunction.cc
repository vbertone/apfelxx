//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/distributionfunction.h"

namespace apfel {

  //_________________________________________________________________________
  DistributionFunction::DistributionFunction(Grid                                        const& g,
					     function<double(int const&, double const&)> const& InDistFunc,
					     int                                         const& ipdf):
    Distribution(g)
  {
    for (auto const& ix: _grid.GetJointGrid().GetGrid())
      _distributionJointGrid.push_back(InDistFunc(ipdf,ix < 1 ? ix : 1));

    for (auto ig = 0; ig < _grid.nGrids(); ig++)
      {
	vector<double> sg;
	for (auto const& ix: _grid.GetSubGrid(ig).GetGrid())
	  sg.push_back(InDistFunc(ipdf,ix < 1 ? ix : 1));

	_distributionSubGrid.push_back(sg);
      }
  }

  //_________________________________________________________________________
  DistributionFunction::DistributionFunction(Grid                                                       const& g,
					     function<double(int const&, double const&, double const&)> const& InDistFunc,
					     int                                                        const& ipdf,
					     double                                                     const& Q):
    Distribution(g)
  {
    for (auto const& ix: _grid.GetJointGrid().GetGrid())
      _distributionJointGrid.push_back(InDistFunc(ipdf,ix < 1 ? ix : 1,Q));

    for (auto ig = 0; ig < _grid.nGrids(); ig++)
      {
	vector<double> sg;
	for (auto const& ix: _grid.GetSubGrid(ig).GetGrid())
	  sg.push_back(InDistFunc(ipdf,ix < 1 ? ix : 1,Q));

	_distributionSubGrid.push_back(sg);
      }
  }

  //_________________________________________________________________________
  DistributionFunction::DistributionFunction(Grid                   const& g,
					     vector<double>         const& DistJointGrid,
					     vector<vector<double>> const& DistSubGrid):
    Distribution(g)
  {
    _distributionJointGrid = DistJointGrid;
    _distributionSubGrid   = DistSubGrid;
  }

}
