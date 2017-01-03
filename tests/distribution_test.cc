//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <apfel/distribution.h>
#include <apfel/grid.h>
#include <apfel/subgrid.h>
#include <apfel/timer.h>

using namespace apfel;
using namespace std;

/**
 * Prototype function for testing purproses.
 */
double xg(double const& x)
{
  if(x < 1) return x * ( 1 - x );
  else      return 0;
}

/**
 * @brief The Parton class
 */
class myPDF: public Distribution
{
public:
  /**
   * Allocate the langrage interpolation and fill the inherited
   * \c _distribution object with the jointed grid.
   */
  myPDF(Grid const& gr): Distribution(gr)
  {
    for (auto const& ix: _grid.GetJointGrid().GetGrid())
      _distributionJointGrid.push_back(xg(ix));

    for (auto ig=0; ig<_grid.nGrids(); ig++)
      {
        vector<double> sg;
        for (auto const& ix: _grid.GetSubGrid(ig).GetGrid())
          sg.push_back(xg(ix));
        _distributionSubGrid.push_back(sg);
      }
  }
};

int main()
{
  Timer t;
  t.start();

  const Grid g{
    {SubGrid{80,1e-5,3}, SubGrid{50,1e-1,5}, SubGrid{40,8e-1,5}}, false
  };
  const myPDF p{g};

  t.printTime(t.stop());

  return 0;
}
