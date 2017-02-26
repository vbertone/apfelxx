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
#include <functional>

using namespace apfel;
using namespace std;

/**
 * @brief A very simple example of ConvolutionMap derivation.
 *
 * This class, following the derivation procedure from ConvolutionMap
 * implements the Basis enumerator with custom tags for the objects.
 */
class EvolutionBasis: public ConvolutionMap
{
public:
  /**
   * @brief The map enums
   */
  enum Operand: int {PNSP, PNSM, PNSV, PQQ, PQG, PGQ, PGG};
  enum Object:  int {GLUON, SIGMA, VALENCE, T3, V3, T8, V8, T15, V15, T24, V24, T35, V35};

  /**
   * @brief The class constructor
   */
  EvolutionBasis(int const& nf):
    ConvolutionMap{"EvolutionBasis"}
  {
    // dg = Pgg * g + Pgq * Sigma
    _rules[GLUON] = { {PGG, GLUON, +1}, {PGQ, SIGMA, +1} };

    // dSigma = Pqg * g + ( Pnsp + Pps ) * Sigma
    _rules[SIGMA] = { {PQG, GLUON, +1}, {PQQ, SIGMA, +1} };

    // dV = Pnsv * V
    _rules[VALENCE] = { {PNSV, VALENCE, +1} };

    // d{T,V}3 = Pnsp * {T,V}3
    if (nf > 1)
      {
	_rules[T3] = { {PNSP, T3, +1} };
	_rules[V3] = { {PNSM, V3, +1} };
      }
    else
      {
	_rules[T3] = _rules[SIGMA];
	_rules[V3] = _rules[VALENCE];
      }

    // d{T,V}8 = Pnsp * {T,V}8
    if (nf > 2)
      {
	_rules[T8] = { {PNSP, T8, +1} };
	_rules[V8] = { {PNSM, V8, +1} };
      }
    else
      {
	_rules[T8] = _rules[SIGMA];
	_rules[V8] = _rules[VALENCE];
      }

    // d{T,V}15 = Pnsp * {T,V}15
    if (nf > 3)
      {
	_rules[T15] = { {PNSP, T15, +1} };
	_rules[V15] = { {PNSM, V15, +1} };
      }
    else
      {
	_rules[T15] = _rules[SIGMA];
	_rules[V15] = _rules[VALENCE];
      }

    // d{T,V}24 = Pnsp * {T,V}24
    if (nf > 4)
      {
	_rules[T24] = { {PNSP, T24, +1} };
	_rules[V24] = { {PNSM, V24, +1} };
      }
    else
      {
	_rules[T24] = _rules[SIGMA];
	_rules[V24] = _rules[VALENCE];
      }

    // d{T,V}35 = Pnsp * {T,V}35
    if (nf > 5)
      {
	_rules[T35] = { {PNSP, T35, +1} };
	_rules[V35] = { {PNSM, V35, +1} };
      }
    else
      {
	_rules[T35] = _rules[SIGMA];
	_rules[V35] = _rules[VALENCE];
      }
  };
};

// LH Toy PDFs
double xupv(double const& x)  { return 5.107200 * pow(x,0.8) * pow((1-x),3); }
double xdnv(double const& x)  { return 3.064320 * pow(x,0.8) * pow((1-x),4); }
double xglu(double const& x)  { return 1.7 * pow(x,-0.1) * pow((1-x),5); }
double xdbar(double const& x) { return 0.1939875 * pow(x,-0.1) * pow((1-x),6); }
double xubar(double const& x) { return xdbar(x) * (1-x); }
double xsbar(double const& x) { return 0.2 * ( xdbar(x) + xubar(x) ); }
double LHToyPDFs(int const& i, double const& x)
{
  // Gluon
  if      (i == EvolutionBasis::GLUON    ) return xglu(x);
  // Singlet, T15, T24, T35
  else if (i == EvolutionBasis::SIGMA   ||
	   i == EvolutionBasis::T15     ||
	   i == EvolutionBasis::T24     ||
	   i == EvolutionBasis::T35      ) return xdnv(x) + 2 * xdbar(x) + xupv(x) + 2 * xubar(x) + 2 * xsbar(x);
  // T3
  else if (i == EvolutionBasis::T3       ) return xupv(x) + 2 * xubar(x) - xdnv(x) - 2 * xdbar(x);
  // T8
  else if (i == EvolutionBasis::T8       ) return xupv(x) + 2 * xubar(x) + xdnv(x) + 2 * xdbar(x) - 4 * xsbar(x);
  // Valence, V8, V15, V24, V35
  else if (i == EvolutionBasis::VALENCE ||
	   i == EvolutionBasis::V8      ||
	   i == EvolutionBasis::V15     ||
	   i == EvolutionBasis::V24     ||
	   i == EvolutionBasis::V35      ) return xupv(x) + xdnv(x);
  // V3
  else if (i == EvolutionBasis::V3       )  return xupv(x) - xdnv(x);
  else              return 0;
}

// The PDF class
class PDF: public Distribution
{
public:
  // Standard constructor
  PDF(Grid const& gr, function<double(int,double)> const& inPDFs, int const& ipdf): Distribution(gr)
  {
    for (auto const& ix: _grid.GetJointGrid().GetGrid())
      if (ix < 1) _distributionJointGrid.push_back(inPDFs(ipdf,ix));
      else        _distributionJointGrid.push_back(0);

    for (auto ig=0; ig<_grid.nGrids(); ig++)
      {
        vector<double> sg;
        for (auto const& ix: _grid.GetSubGrid(ig).GetGrid())
          if (ix < 1) sg.push_back(inPDFs(ipdf,ix));
          else        sg.push_back(0);
        _distributionSubGrid.push_back(sg);
      }
  }
};

/**
 * @brief The LO splitting function classes
 */
class P0ns: public Expression
{
public:
  P0ns(): Expression() { }
  double Regular(double const& x)  const { return - 2 * CF * ( 1 + x ); }
  double Singular(double const& x) const { return 4 * CF / ( 1 - x ); }
  double Local(double const& x)    const { return 4 * CF * log( 1 - x ) + 3 * CF; }
};

class P0qg: public Expression
{
public:
  P0qg(): Expression() { }
  double Regular(double const& x)  const { return 2 * ( 1 - 2 * x + 2 * x * x ); }
};

class P0gq: public Expression
{
public:
  P0gq(): Expression() { }
  double Regular(double const& x)  const { return 4 * CF * ( - 1 + 0.5 * x + 1 / x ); }
};

class P0gg: public Expression
{
public:
  P0gg(int const& nf): Expression(), _nf(nf) { }
  double Regular(double const& x)  const { return 4 * CA * ( - 2 + x - x * x + 1 / x ); }
  double Singular(double const& x) const { return 4 * CA / ( 1 - x ); }
  double Local(double const& x)    const { return 4 * CA * log( 1 - x ) - 2 / 3. * _nf + 11 / 3. * CA; }
private:
  int const _nf;
};

int main()
{
  // Time counter
  Timer t,ttot;
  ttot.start();

  // Grid
  const Grid g{{SubGrid{80,1e-5,3}, SubGrid{50,1e-1,5}, SubGrid{40,8e-1,5}}};

  // ===============================================================
  // Allocate LO splitting functions operators
  cout << "Initializing operators ..." << endl;
  t.start();
  unordered_map<int,unordered_map<int,Operator>> OpMap;
  const Operator O0ns{g, P0ns{}};
  const Operator O0qg{g, P0qg{}};
  const Operator O0gq{g, P0gq{}};
  for (int nf = 3; nf <= 6; nf++)
    {
      const Operator O0gg{g, P0gg{nf}};
      const Operator O0qgnf = nf * O0qg;
      unordered_map<int,Operator> OM;
      OM.insert({EvolutionBasis::PNSP,O0ns});
      OM.insert({EvolutionBasis::PNSM,O0ns});
      OM.insert({EvolutionBasis::PNSV,O0ns});
      OM.insert({EvolutionBasis::PQQ, O0ns});
      OM.insert({EvolutionBasis::PQG, O0qgnf});
      OM.insert({EvolutionBasis::PGQ, O0gq});
      OM.insert({EvolutionBasis::PGG, O0gg});
      OpMap.insert({nf,OM});
    }
  t.printTime(t.stop());

  // Allocate distributions
  cout << "Initializing distributions ..." << endl;
  t.start();
  unordered_map<int,Distribution> DistMap;
  for (int i = EvolutionBasis::GLUON; i <= EvolutionBasis::V35; i++)
    DistMap.insert({i,PDF{g, LHToyPDFs, i}});
  t.printTime(t.stop());

  cout << "Initializing set of operators and distributions ..." << endl;
  t.start();
  // Allocate maps
  unordered_map<int,EvolutionBasis> basis;
  for (int nf = 3; nf <= 6; nf++)
    basis.insert({nf,EvolutionBasis{nf}});

  // Allocate set of operators
  unordered_map<int,Set<Operator>> Splittings;
  for (int nf = 3; nf <= 6; nf++)
    Splittings.insert({nf,Set<Operator>{basis.at(nf), OpMap.at(nf)}});

  // Allocate set of initial distributions
  Set<Distribution> PDFs{basis.at(5), DistMap};

  t.printTime(t.stop());

  // ===============================================================
  // Test products
  cout << "\nTesting products ..." << endl;
  t.start();
  auto Product = Splittings.at(5) * PDFs;
  cout << "(Splitting * PDFs)[GLUON](x=0.1) = "
       << Product.at(EvolutionBasis::GLUON).Evaluate(0.1) << endl;

  auto Product2 = 2 * Product;
  cout << "(2 * Splitting * PDFs)[GLUON](x=0.1) = "
       << Product2.at(EvolutionBasis::GLUON).Evaluate(0.1) << endl;

  auto Sum = Product.at(EvolutionBasis::GLUON) + Product.at(EvolutionBasis::GLUON);
  cout << "[(Splitting * PDFs)[GLUON] + (Splitting * PDFs)[GLUON]](x=0.1) = "
       << Sum.Evaluate(0.1) << endl;
  t.printTime(t.stop());

  // ===============================================================
  cout << "\nFull computation time:" << endl;
  ttot.printTime(ttot.stop());

  return 0;
}
