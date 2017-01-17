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
#include <map>

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
  // Void constructor
  PDF(Grid const& gr): Distribution(gr)
  {
    _distributionJointGrid.resize(_grid.GetJointGrid().GetGrid().size());
    for (auto ig=0; ig<_grid.nGrids(); ig++)
      {
        vector<double> sg;
	sg.resize(_grid.GetSubGrid(ig).GetGrid().size());
        _distributionSubGrid.push_back(sg);
      }
  }
};

// The structure function class
class StructureFunction: public Distribution
{
public:
  StructureFunction(Grid const& gr): Distribution(gr)
  {
    _distributionJointGrid.resize(_grid.GetJointGrid().GetGrid().size());

    for (auto ig=0; ig<_grid.nGrids(); ig++)
      {
        vector<double> sg(_grid.GetSubGrid(ig).GetGrid().size(),0);
        _distributionSubGrid.push_back(sg);
      }
  }
};

// Enumerator to identify the cf. component
enum comp { G = 0, PS = 1, NSP = 2, NSM = 3, V = 4 };

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
    if (_cmp == NSP)
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
    if (_cmp == NSP)
      if(_pt == 0) return 0;
      else         return C2qS(x);
    else
      return 0;
  }
  double Local(double const& x)    const
  {
    if (_cmp == NSP)
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
  // Time counter
  Timer t;

  // Define the scale
  const auto Q = sqrt(2);

  // ========== Initialization phase ==========
  t.start();

  // Grid
  const Grid g{{SubGrid{80,1e-5,3}, SubGrid{50,1e-1,5}, SubGrid{40,8e-1,5}}, false};

  // PDFs
  vector<PDF> pdfs;
  for (auto i = -6; i < 7; i++) 
    {
      PDF d{i,g};
      pdfs.push_back(d);
    }

  // Alphas
  const AlphaQCD Coup{0.35, sqrt(2), {0, 0, 0, sqrt(2), 4.5, 175}, 1};
  const auto as = Coup.GetCoupling(Q) / FourPi;

  // Coefficient functions
  const f2cf C2qLO{NSP, 0};
  const f2cf C2qNLO{NSP, 1};
  const f2cf C2gNLO{G, 1};

  // Construct operators
  const Operator OqLO{g, C2qLO};
  const Operator OqNLO{g, C2qNLO};
  const Operator OgNLO{g, C2gNLO};

  // Combine operators with alphas and put it in a map
  const map<comp, Operator> Om = { { G, as * OgNLO }, { NSP, OqLO + as * OqNLO } };

  // Combine operators with alphas and put it in a vector
  const vector<Operator> Ov = { OqLO + as * OqNLO, as * OgNLO };

  cout << "Initialization ..." << endl;
  t.printTime(t.stop());

  // Charges
  const auto eu2 = 4. / 9.;
  const auto ed2 = 1. / 9.;
  const auto eg2 = 2 * ed2 + 2 * eu2;

  // Construct the structure function as a combination of distributions and operators
  // ========== Computation phase ==========
  t.start();

  // Define map
  map <comp, map<int,double>> Map;
  Map[G][6]   = eg2;
  Map[NSP][3] = Map[NSP][9] = Map[NSP][5] = Map[NSP][7] = ed2;
  Map[NSP][4] = Map[NSP][8] = eu2;

  // Combine distributions and operators according to the map a produce a StructureFunction object
  StructureFunction F2map{g};
  for (auto const& iO : Map)
    for (auto const& id : iO.second)
      F2map += id.second * Om.at(iO.first) * pdfs[id.first];

  // Print the results
  cout << scientific;
  const vector<double> xlha = { 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1 };
  for (auto i = 2; i < (int) xlha.size(); i++)
    cout << "F2(x = " << xlha[i] << ", Q = " << Q << " GeV) = " << F2map.Evaluate(xlha[i]) << endl;

  cout << "Computation ..." << endl;
  t.printTime(t.stop());

  // Alterantive (and more optimal) way to combine Operators and distributions
  // ========== Computation phase ==========
  t.start();

  // Define multimap bewteen operators and distributions
  multimap<int,int> MMap;
  MMap.insert(pair<int,int>(0,3)); // sbar
  MMap.insert(pair<int,int>(0,4)); // ubar
  MMap.insert(pair<int,int>(0,5)); // dbar
  MMap.insert(pair<int,int>(1,6)); // gluon
  MMap.insert(pair<int,int>(0,7)); // d
  MMap.insert(pair<int,int>(0,8)); // u
  MMap.insert(pair<int,int>(0,9)); // s

  // Charges
  vector<double> ch2 = { eu2, ed2, eu2, ed2, eu2, ed2, eg2, ed2, eu2, ed2, eu2, ed2, eu2 };

  // Construct sturcture functions
  StructureFunction F2vec{g};
  for (auto i = 0; i< (int) Ov.size(); i++)
    {
      // Combine distributions according to charges
      PDF pdfcomb{g};
      pair <multimap<int,int>::iterator, multimap<int,int>::iterator> ret;
      ret = MMap.equal_range(i);
      for (multimap<int,int>::iterator it = ret.first; it != ret.second; ++it)
	pdfcomb += ch2[it->second] * pdfs[it->second];

      // Multiply combined distribution with the appropriate operator and sum it to the structure function
      F2vec += Ov[i] * pdfcomb;
    }

  for (auto i = 2; i < (int) xlha.size(); i++)
    cout << "F2(x = " << xlha[i] << ", Q = " << Q << " GeV) = " << F2vec.Evaluate(xlha[i]) << endl;

  cout << "Computation ..." << endl;
  t.printTime(t.stop());

  return 0;
}
