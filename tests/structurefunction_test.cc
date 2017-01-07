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
#include <cmath>

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
  if      (i == 3)  return xsbar(x);
  else if (i == 2)  return xupv(x) + xubar(x);
  else if (i == 1)  return xdnv(x) + xdbar(x);
  else if (i == 0)  return xglu(x);
  else if (i == -1) return xdbar(x);
  else if (i == -2) return xubar(x);
  else if (i == -3) return xsbar(x);
  else              return 0;
}

// Zero mass NLO coefficient functions for F2
double C2qR(double const& x) { return 2. * CF * ( - ( 1. + x ) * log( 1. - x ) - ( 1. + pow(x,2) ) * log(x) / ( 1. - x ) + 3. + 2. * x ); }
double C2qS(double const& x) { return 2. * CF * ( 2. * log( 1 - x ) - 3. / 2. ) / ( 1. - x ); }
double C2qL(double const& x) { return 2. * CF * ( pow(log( 1. - x ),2) - 3. * log( 1. - x ) / 2. - ( 2. * zeta2 + 9. / 2. ) ); }
double C2gR(double const& x) { return 4. * TR * ( ( pow(( 1. - x ),2) + pow(x,2) ) * log( ( 1. - x ) / x ) - 8. * pow(x,2) + 8. * x - 1. ); }

// The PDF class
class PDF: public Distribution
{
public:
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
};

// Enumerator to identify the cf. component
enum comp { Q = 0, G = 1 };

// The coefficient function class
class f2cf: public Expression
{
public:
  f2cf(comp const& cmp, int const& pt):
    Expression(),
    _cmp(cmp),
    _pt(pt)
  {}
  double Regular(double const& x) const
  {
    if (_cmp == Q)
	if(_pt == 0) return 0;
	else         return C2qR(x);
    else if (_cmp == G)
	if(_pt == 0) return 0;
	else         return C2gR(x);
    else
      return 0;
  }
  double Singular(double const& x) const
  {
    if (_cmp == Q)
	if(_pt == 0) return 0;
	else         return C2qS(x);
    else
      return 0;
  }
  double Local(double const& x)    const
  {
    if (_cmp == Q)
	if(_pt == 0) return 1;
	else         return C2qL(x);
    else
      return 0;
  }
private:
  comp const _cmp;
  int  const _pt;
};


int main()
{
  Timer t;
  t.start();

  // Define x and scale
  const auto x  = 0.1;
  const auto Mu = sqrt(2);

  // Grid
  const Grid g{{SubGrid{80,1e-5,3}, SubGrid{50,1e-1,3}, SubGrid{40,8e-1,3}}, false};

  // PDFs
  vector<PDF> pdfs;
  for (auto i = -6; i < 7; i++) 
    {
      PDF d{i,g};
      pdfs.push_back(d);
    }

  // Alphas
  const AlphaQCD Coup(0.35, sqrt(2), {0, 0, 0, sqrt(2), 4.5, 175}, 1);
  const auto as = Coup.GetCoupling(Mu) / FourPi;

  // Coefficient functions
  const f2cf C2qLO(Q, 0);
  const f2cf C2qNLO(Q, 1);
  const f2cf C2gNLO(G, 1);

  // Construct operators
  const Operator OqLO{g, C2qLO};
  const Operator OqNLO{g, C2qNLO};
  const Operator OgNLO{g, C2gNLO};

  // Combine operators with alphas
  const auto Oq = OqLO + as * OqNLO;
  const auto Og =        as * OgNLO;

  // Construct the structure function as a combination of distributions and operators
  const auto eu2 = 4. / 9.;
  const auto ed2 = 1. / 9.;
  const auto F2 = Oq * ( ed2 * ( pdfs[3] + pdfs[9] + pdfs[5] + pdfs[7] ) + eu2 * ( pdfs[4] + pdfs[8] ) ) + 2 * ( ed2 + eu2 ) * Og * pdfs[6];

  cout << "\nF2(x = " << x << ", Q = " << Mu << ") = " << F2.Evaluate(x) << "\n" << endl;

  t.printTime(t.stop());

  return 0;
}
