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

/**
 * @brief A very simple example of BasisMap derivation.
 *
 * This class, following the derivation procedure from BasisMap
 * implements the Basis enumerator with custom tags for the objects.
 */
class EvolutionMap: public BasisMap
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
  EvolutionMap(int const& nf):
    BasisMap{"EvolutionMap"}
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
  if      (i == EvolutionMap::GLUON    ) return xglu(x);
  // Singlet, T15, T24, T35
  else if (i == EvolutionMap::SIGMA   ||
	   i == EvolutionMap::T15     ||
	   i == EvolutionMap::T24     ||
	   i == EvolutionMap::T35      ) return xdnv(x) + 2 * xdbar(x) + xupv(x) + 2 * xubar(x) + 2 * xsbar(x);
  // T3
  else if (i == EvolutionMap::T3       ) return xupv(x) + 2 * xubar(x) - xdnv(x) - 2 * xdbar(x);
  // T8
  else if (i == EvolutionMap::T8       ) return xupv(x) + 2 * xubar(x) + xdnv(x) + 2 * xdbar(x) - 4 * xsbar(x);
  // Valence, V8, V15, V24, V35
  else if (i == EvolutionMap::VALENCE ||
	   i == EvolutionMap::V8      ||
	   i == EvolutionMap::V15     ||
	   i == EvolutionMap::V24     ||
	   i == EvolutionMap::V35      ) return xupv(x) + xdnv(x);
  // V3
  else if (i == EvolutionMap::V3       )  return xupv(x) - xdnv(x);
  else              return 0;
}

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
};

/**
 * @brief The LO splitting function class
 */
class P0: public Expression
{
public:
  P0(int const& iop, int const& nf):
    Expression(),
    _iop(iop),
    _nf(nf)
  {}
  double Regular(double const& x) const
  {
    switch(_iop)
      {
	case 0:
	  return - 2 * CF * ( 1 + x );
	case 1:
	  return - 2 * CF * ( 1 + x );
	case 2:
	  return - 2 * CF * ( 1 + x );
	case 3:
	  return - 2 * CF * ( 1 + x );
	case 4:
	  return 2 * _nf * ( 1 - 2 * x + 2 * x * x );
	case 5:
	  return 4 * CF * ( - 1 + 0.5 * x + 1 / x );
	case 6:
	  return 4 * CA * ( - 2 + x - x * x + 1 / x );
	default:
	  return 0.;
      }
  }
  double Singular(double const& x) const
  {
    switch(_iop)
      {
	case 0:
	  return 4 * CF / ( 1 - x );
	case 1:
	  return 4 * CF / ( 1 - x );
	case 2:
	  return 4 * CF / ( 1 - x );
	case 3:
	  return 4 * CF / ( 1 - x );
	case 6:
	  return 4 * CA / ( 1 - x );
	default:
	  return 0.;
      }
  }
  double Local(double const& x) const
  {
    switch(_iop)
      {
	case 0:
	  return 4 * CF * log( 1 - x ) + 3 * CF;
	case 1:
	  return 4 * CF * log( 1 - x ) + 3 * CF;
	case 2:
	  return 4 * CF * log( 1 - x ) + 3 * CF;
	case 3:
	  return 4 * CF * log( 1 - x ) + 3 * CF;
	case 6:
	  return 4 * CA * log( 1 - x ) - 2 / 3. * _nf + 11 / 3. * CA;
	default:
	  return 0.;
      }
  }
private:
  int const _iop;
  int const _nf;
};

int main()
{
  // Time counter
  Timer t,ttot;
  ttot.start();

  // Grid
  const Grid g{{SubGrid{80,1e-5,3}, SubGrid{50,1e-1,5}, SubGrid{40,8e-1,5}}, false};

  // Allocate map
  EvolutionMap basis(6);

  // Allocate operators
  // Brute force: compute all operators for all nf's.
  // At LO, this is not optimal because most of the operators do not
  // depend on nf. In addition most of the operators are equal.
  // This can be optimized.
  cout << "Initializing operators ..." << endl;
  t.start();
  vector<unordered_map<int,Operator>> OpMap;
  for (auto nf = 3; nf <= 6; nf++)
    {
      unordered_map<int,Operator> OM;
      for (int i = EvolutionMap::PNSP; i <= EvolutionMap::PGG; i++)
	{
	  const Operator O{g, P0{i, nf}};
	  OM.insert({i,O});
	}
      OpMap.push_back(OM);
    }
  t.printTime(t.stop());

  // Allocate distributions
  cout << "Initializing distributions ..." << endl;
  t.start();
  unordered_map<int,Distribution> DistMap;
  for (int i = EvolutionMap::GLUON; i <= EvolutionMap::V35; i++)
    {
      const PDF f{i, g};
      DistMap.insert({i,f});
    }
  t.printTime(t.stop());

  cout << "Initializing set of operators and distributions ..." << endl;
  t.start();
  // Allocate set of operators
  Set<Operator> Splittings{basis, OpMap[2]};

  // Allocate set of operators
  Set<Distribution> PDFs{basis, DistMap};
  t.printTime(t.stop());

  // Test products
  cout << "\nTesting products ..." << endl;
  t.start();
  auto Product = Splittings * PDFs;
  cout << "(Splitting * PDFs)[GLUON](x=0.1) = "
       << Product.at(EvolutionMap::GLUON).Evaluate(0.1) << endl;

  auto Product2 = 2 * Product;
  cout << "(2 * Splitting * PDFs)[GLUON](x=0.1) = "
       << Product2.at(EvolutionMap::GLUON).Evaluate(0.1) << endl;

  auto Sum = Product.at(EvolutionMap::GLUON) + Product.at(EvolutionMap::GLUON);
  cout << "[(Splitting * PDFs)[GLUON] + (Splitting * PDFs)[GLUON]](x=0.1) = "
       << Sum.Evaluate(0.1) << endl;
  t.printTime(t.stop());

  cout << "\nFull computation time:" << endl;
  ttot.printTime(ttot.stop());

  return 0;
}
