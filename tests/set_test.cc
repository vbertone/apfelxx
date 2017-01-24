//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <apfel/distribution.h>
#include <apfel/grid.h>
#include <apfel/subgrid.h>
#include <apfel/operator.h>
#include <apfel/expression.h>
#include <apfel/timer.h>
#include <apfel/tools.h>
#include <apfel/alphaqcd.h>
#include <apfel/set.h>
#include <cmath>
#include <map>

using namespace apfel;
using namespace std;

// The PDF class
class PDF: public Distribution
{
public:
  // Standard constructor
  PDF(int const& ipdf, Grid const& gr): Distribution(gr)
  {
    for (auto const& ix: _grid.GetJointGrid().GetGrid())
      if (ix < 1) _distributionJointGrid.push_back(LHToyPDFs(ipdf,ix));
      else        _distributionJointGrid.push_back(0);

    for (auto ig=0; ig<_grid.nGrids(); ig++)
      {
        vector<double> sg;
        for (auto const& ix: _grid.GetSubGrid(ig).GetGrid())
          if (ix < 1) sg.push_back(LHToyPDFs(ipdf,ix));
          else        sg.push_back(0);
        _distributionSubGrid.push_back(sg);
      }
  }  

  // LH Toy PDFs
  double xupv(double const& x)  { return 5.107200 * pow(x,0.8) * pow((1-x),3); }
  double xdnv(double const& x)  { return 3.064320 * pow(x,0.8) * pow((1-x),4); }
  double xglu(double const& x)  { return 1.7 * pow(x,-0.1) * pow((1-x),5); }
  double xdbar(double const& x) { return 0.1939875 * pow(x,-0.1) * pow((1-x),6); }
  double xubar(double const& x) { return xdbar(x) * (1-x); }
  double xsbar(double const& x) { return 0.2 * ( xdbar(x) + xubar(x) ); }
  double LHToyPDFs(int const& i, double const& x)
  {
    if      (i == 3)  return xsbar(x);
    else if (i == 2)  return xupv(x) + xubar(x);
    else if (i == 1)  return xdnv(x) + xdbar(x);
    else if (i == 0)  return xglu(x);
    else if (i == -1) return xdbar(x);
    else if (i == -2) return xubar(x);
    else if (i == -3) return xsbar(x);
    else              return 0;
  }

};

/**
 * @brief The p0qq class
 */
class p0qq: public Expression
{
public:
  p0qq(): Expression() {}
  double Regular(double const& x)  const { return - 2 * CF * ( 1 + x ); }
  double Singular(double const& x) const { return 4 * CF / ( 1 - x ); }
  double Local(double const& x)    const { return 4 * CF * log( 1 - x ) + 3 * CF; }
};

/**
 * @brief The p0gg class
 */
class p0gg: public Expression
{
public:
  p0gg(): Expression() {}
  double Regular(double const& x)  const { return 4 * CA * (-2 + x - x*x + 1./x); }
  double Singular(double const& x) const { return 4 * CA / ( 1 - x ); }
  double Local(double const& x)    const { return 4 * CA * log(1-x) - 2/3.*5 + 11/3.*CA; }
};

/**
 * @brief The p0qg class
 */
class p0qg: public Expression
{
public:
  p0qg(): Expression() {}
  double Regular(double const& x)  const { return 2 * 5 * (1-2*x + 2*x*x); }
  double Singular(double const& x) const { return 0*x; }
  double Local(double const& x) const { return 0*x; }
};

/**
 * @brief The p0gq class
 */
class p0gq: public Expression
{
public:
  p0gq(): Expression() {}
  double Regular(double const& x)  const { return 4 * CF * (-1 + 0.5*x + 1./x); }
  double Singular(double const& x) const { return 0*x; }
  double Local(double const& x) const { return 0*x; }
};

int main()
{
  // Time counter
  Timer t;
  t.start();

  // Grid
  const Grid grid{{SubGrid{80,1e-5,3}, SubGrid{50,1e-1,5}, SubGrid{40,8e-1,5}}, false};

  // load PDFs in flavor basis
  const PDF sbar{-3, grid};
  const PDF ubar{-2, grid};
  const PDF dbar{-1, grid};
  const PDF g{0, grid};
  const PDF d{1, grid};
  const PDF u{2, grid};
  const PDF s{3, grid};

  const Operator pqq{grid, p0qq{}};
  const Operator pgg{grid, p0gg{}};
  const Operator pgq{grid, p0gq{}};
  const Operator pqg{grid, p0qg{}};

  unordered_map<int,Operator> omap = {
                                      {FlvrBasis::PGG, pgg},
                                      {FlvrBasis::PQQ, pqq},
                                      {FlvrBasis::PQG, pqg},
                                      {FlvrBasis::PGQ, pgq}
                                     };

  unordered_map<int,Distribution> dmap = {
                                          {FlvrBasis::SBAR, sbar},
                                          {FlvrBasis::UBAR, ubar},
                                          {FlvrBasis::DBAR, dbar},
                                          {FlvrBasis::GLU,  g},
                                          {FlvrBasis::D,    d},
                                          {FlvrBasis::U,    u},
                                          {FlvrBasis::S,    s}
                                        };

  FlvrBasis basis;

  // allocating operators following the flavor basis
  Set<Operator> splittings(basis, omap);

  // allocating PDFs following the flavor basis
  Set<Distribution> pdfs(basis, dmap);

  // testing product
  auto product = splittings * pdfs;

  // getting new distribution
  cout << "(Splitting * PDFs)[GLUON](x=0.1) = "
       << product.at(FlvrBasis::GLU).Evaluate(0.1) << endl;

  t.printTime(t.stop());

  return 0;
}
