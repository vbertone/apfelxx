//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <cmath>
#include <map>
#include <iomanip>

#include <apfel/zeromasscoefficientfunctions.h>
#include <apfel/disbasis.h>
#include <apfel/expression.h>
#include "apfel/operator.h"
#include "apfel/set.h"
#include <apfel/distribution.h>
#include <apfel/grid.h>
#include <apfel/subgrid.h>
#include <apfel/operator.h>
#include <apfel/expression.h>
#include <apfel/timer.h>
#include <apfel/tools.h>
#include <apfel/alphaqcd.h>

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
  // Gluon
  if      (i == 0)
    return xglu(x);
  // Singlet, T15, T24, T35
  else if (i == 1 || i == 7 || i == 9 || i == 11 )
    return xdnv(x) + 2 * xdbar(x) + xupv(x) + 2 * xubar(x) + 2 * xsbar(x);
  // T3
  else if (i == 3)
    return xupv(x) + 2 * xubar(x) - xdnv(x) - 2 * xdbar(x);
  // T8
  else if (i == 5)
    return xupv(x) + 2 * xubar(x) + xdnv(x) + 2 * xdbar(x) - 4 * xsbar(x);
  // Valence, V8, V15, V24, V35
  else if (i == 2 || i == 6 || i == 8 || i == 10 || i == 12)
    return xupv(x) + xdnv(x);
  // V3
  else if (i == 4)
    return xupv(x) - xdnv(x);
  else
    return 0;
}

/**
 * @brief The PDF class
 */
class PDF: public Distribution
{
public:
  PDF(Grid                                        const& g,
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
};

int main()
{
  // Time counter
  Timer t;

  // Define the scale and number of active flavours
  const auto Q = sqrt(2);
  const auto nf = 3;

  // Grid
  const Grid g{{SubGrid{80,1e-5,3}, SubGrid{50,1e-1,5}, SubGrid{40,8e-1,5}}, false};

  // Alphas
  const AlphaQCD Coup{0.35, sqrt(2), {0, 0, 0, sqrt(2), 4.5, 175}, 2};
  const auto as  = Coup.Evaluate(Q) / FourPi;
  const auto as2 = as * as;

  // Charges
  const auto eu2 = 4. / 9.;
  const auto ed2 = 1. / 9.;
  vector<double> Bq = { eu2, ed2, ed2, eu2, ed2, eu2 };
  vector<double> Dq = { 0, 0, 0, 0, 0, 0 };

  // Integration accuracy
  const auto IntEps = 1e-5;

  t.start();
  cout << "Initialization ..." << endl;

  // ===============================================================
  // Allocate convolution maps for all structure functions
  unordered_map<int,DISNCBasis> cmap;
  for (int k = 1; k <= 6; k++)
    cmap.insert({k,DISNCBasis{k,nf}});

  // Allocate distributions for F2
  unordered_map<int,Distribution> F2Map;
  F2Map.insert({0,PDF{g, LHToyPDFs, 0 }});
  F2Map.insert({1,PDF{g, LHToyPDFs, 1 }});
  F2Map.insert({2,PDF{g, LHToyPDFs, 3 }});
  F2Map.insert({3,PDF{g, LHToyPDFs, 5 }});
  F2Map.insert({4,PDF{g, LHToyPDFs, 7 }});
  F2Map.insert({5,PDF{g, LHToyPDFs, 9 }});
  F2Map.insert({6,PDF{g, LHToyPDFs, 11}});

  // Allocate distributions for FL
  unordered_map<int,Distribution> FLMap = F2Map;

  // Allocate distributions for F3
  unordered_map<int,Distribution> F3Map;
  F3Map.insert({0,PDF{g, LHToyPDFs, 0 }});
  F3Map.insert({1,PDF{g, LHToyPDFs, 2 }});
  F3Map.insert({2,PDF{g, LHToyPDFs, 4 }});
  F3Map.insert({3,PDF{g, LHToyPDFs, 6 }});
  F3Map.insert({4,PDF{g, LHToyPDFs, 8 }});
  F3Map.insert({5,PDF{g, LHToyPDFs, 10}});
  F3Map.insert({6,PDF{g, LHToyPDFs, 12}});

  // Allocate set of distributions
  unordered_map<int,Set<Distribution>> sPDFsF2;
  unordered_map<int,Set<Distribution>> sPDFsFL;
  unordered_map<int,Set<Distribution>> sPDFsF3;
  for (int k = 1; k <= 6; k++)
    {
      sPDFsF2.insert({k,Set<Distribution>{cmap.at(k), F2Map}});
      sPDFsFL.insert({k,Set<Distribution>{cmap.at(k), FLMap}});
      sPDFsF3.insert({k,Set<Distribution>{cmap.at(k), F3Map}});
    }

  // ===============================================================
  const Operator Id{g, Identity{}, IntEps};
  const Operator Zero{g, Null{}, IntEps};

  // Coefficient functions for F2
  unordered_map<int,Operator> C2;
  const Operator O21ns{g, C21ns{}, IntEps};
  const Operator O21g{g, C21g{}, IntEps};
  const Operator O22nsp{g, C22nsp{nf}, IntEps};
  const Operator O22ps{g, C22ps{}, IntEps};
  const Operator O22g{g, C22g{}, IntEps};
  const Operator C2ns = Id   + as * O21ns + as2 * O22nsp;
  const Operator C2s  = Id   + as * O21ns + as2 * ( O22nsp + nf * O22ps );
  const Operator C2g  =        as * O21g  + as2 * O22g;
  C2.insert({DISNCBasis::CNS, C2ns});
  C2.insert({DISNCBasis::CT,  C2s});
  C2.insert({DISNCBasis::CG,  C2g});

  // Coefficient functions for FL
  unordered_map<int,Operator> CL;
  const Operator OL1ns{g, CL1ns{}, IntEps};
  const Operator OL1g{g, CL1g{}, IntEps};
  const Operator OL2nsp{g, CL2nsp{nf}, IntEps};
  const Operator OL2ps{g, CL2ps{}, IntEps};
  const Operator OL2g{g, CL2g{}, IntEps};
  const Operator CLns = as * OL1ns + as2 * OL2nsp;
  const Operator CLs  = as * OL1ns + as2 * ( OL2nsp + nf * OL2ps );
  const Operator CLg  = as * OL1g  + as2 * OL2g;
  CL.insert({DISNCBasis::CNS, CLns});
  CL.insert({DISNCBasis::CT,  CLs});
  CL.insert({DISNCBasis::CG,  CLg});

  // Coefficient functions for F3
  unordered_map<int,Operator> C3;
  const Operator O31ns{g, C31ns{}, IntEps};
  const Operator O32nsm{g, C32nsm{nf}, IntEps};
  const Operator C3ns = Id + as * O31ns + as2 * O32nsm;
  C3.insert({DISNCBasis::CNS, C3ns});
  C3.insert({DISNCBasis::CT,  C3ns});
  C3.insert({DISNCBasis::CG,  Zero});

  // Allocate set of operators
  unordered_map<int,Set<Operator>> sC2;
  unordered_map<int,Set<Operator>> sCL;
  unordered_map<int,Set<Operator>> sC3;
  for (int k = 1; k <= 6; k++)
    {
      sC2.insert({k,Set<Operator>{cmap.at(k), C2}});
      sCL.insert({k,Set<Operator>{cmap.at(k), CL}});
      sC3.insert({k,Set<Operator>{cmap.at(k), C3}});
    }

  // Compute structure functions
  unordered_map<int,Distribution> F2;
  unordered_map<int,Distribution> FL;
  unordered_map<int,Distribution> F3;
  for (int k = 1; k <= 6; k++)
    {
      Set<Distribution> sF2 = sC2.at(k) * sPDFsF2.at(k);
      Set<Distribution> sFL = sCL.at(k) * sPDFsFL.at(k);
      Set<Distribution> sF3 = sC3.at(k) * sPDFsF3.at(k);
      Distribution F2k = Bq[k-1] * sF2.Combine();
      Distribution FLk = Bq[k-1] * sFL.Combine();
      Distribution F3k = Dq[k-1] * sF3.Combine();
      F2.insert({k, F2k});
      FL.insert({k, FLk});
      F3.insert({k, F3k});
    }

  // Compute total structure functions
  Distribution TotF2 = F2.at(1);
  Distribution TotFL = FL.at(1);
  Distribution TotF3 = F3.at(1);
  for (int k = 2; k <= nf; k++)
    {
      TotF2 += F2.at(k);
      TotFL += FL.at(k);
      TotF3 += F3.at(k);
    }
  F2.insert({0, TotF2});
  FL.insert({0, TotFL});
  F3.insert({0, TotF3});
  t.stop();

  // Print results
  t.start();
  cout << scientific;

  cout << "Alphas(Q) = " << as * FourPi << endl;
  cout << endl;

  const vector<double> xlha = { 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1 };

  cout << "    x   "
       << "  F2light   "
       << "  F2charm   "
       << "  F2bottom  "
       << "  F2total   "
       << endl;
  for (auto i = 2; i < (int) xlha.size(); i++)
    cout << setprecision(1) << xlha[i] << "  " << setprecision(4)
	 << (F2.at(1)+F2.at(2)+F2.at(3)).Evaluate(xlha[i]) << "  "
	 << F2.at(4).Evaluate(xlha[i]) << "  "
	 << F2.at(5).Evaluate(xlha[i]) << "  "
	 << F2.at(0).Evaluate(xlha[i]) << "  "
	 << endl;
  cout << endl;

  cout << "    x   "
       << "  FLlight   "
       << "  FLcharm   "
       << "  FLbottom  "
       << "  FLtotal   "
       << endl;
  for (auto i = 2; i < (int) xlha.size(); i++)
    cout << setprecision(1) << xlha[i] << "  " << setprecision(4)
	 << (FL.at(1)+FL.at(2)+FL.at(3)).Evaluate(xlha[i]) << "  "
	 << FL.at(4).Evaluate(xlha[i]) << "  "
	 << FL.at(5).Evaluate(xlha[i]) << "  "
	 << FL.at(0).Evaluate(xlha[i]) << "  "
	 << endl;
  cout << endl;

  cout << "    x   "
       << "  F3light   "
       << "  F3charm   "
       << "  F3bottom  "
       << "  F3total   "
       << endl;
  for (auto i = 2; i < (int) xlha.size(); i++)
    cout << setprecision(1) << xlha[i] << "  " << setprecision(4)
	 << (F3.at(1)+F3.at(2)+F3.at(3)).Evaluate(xlha[i]) << "  "
	 << F3.at(4).Evaluate(xlha[i]) << "  "
	 << F3.at(5).Evaluate(xlha[i]) << "  "
	 << F3.at(0).Evaluate(xlha[i]) << "  "
	 << endl;
  cout << endl;

  return 0;
}
