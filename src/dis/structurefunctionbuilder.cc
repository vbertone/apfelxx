//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/structurefunctionbuilder.h"
#include "apfel/grid.h"
#include "apfel/operator.h"
#include "apfel/set.h"
#include "apfel/timer.h"
#include "apfel/zeromasscoefficientfunctions.h"

#include <map>

using namespace std;

namespace apfel {

  //_____________________________________________________________________________
  StructureFunctionObjects InitializeF2ObjectsZM(Grid const& g, double const& IntEps)
  {
    cout << "Initializing StructureFunctionObjects for F2 Zero Mass... ";
    Timer t;
    t.start();

    // Define object of the structure containing the DglapObjects
    StructureFunctionObjects F2Obj;

    // Set initial and final number of active flavours. Precompute
    // objects for all numbers of flavours because it cheap enough.
    int nfi = 1;
    int nff = 6;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // Coefficient functions for F2.
    // LO
    unordered_map<int,Operator> C2LO;
    C2LO.insert({DISNCBasis::CNS, Id});
    C2LO.insert({DISNCBasis::CS,  Id});
    C2LO.insert({DISNCBasis::CG,  Zero});

    // NLO
    unordered_map<int,Operator> C2NLO;
    const Operator O21ns{g, C21ns{}, IntEps};
    const Operator O21g {g, C21g{},  IntEps};
    C2NLO.insert({DISNCBasis::CNS, O21ns});
    C2NLO.insert({DISNCBasis::CS,  O21ns});
    C2NLO.insert({DISNCBasis::CG,  O21g});

    // NNLO
    unordered_map<int,unordered_map<int,Operator>> C2NNLO;
    const Operator O22ps{g, C22ps{}, IntEps};
    const Operator O22g {g, C22g{},  IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
	const Operator O22nsp{g, C22nsp{nf}, IntEps};
	const Operator O22t = O22nsp + 6 * O22ps;
	unordered_map<int,Operator> C2NNLOnf;
	C2NNLOnf.insert({DISNCBasis::CNS, O22nsp});
	C2NNLOnf.insert({DISNCBasis::CS,  O22t});
	C2NNLOnf.insert({DISNCBasis::CG,  O22g});
	C2NNLO.insert({nf,C2NNLOnf});
      }

    // Allocate set of operators.
    F2Obj.skip = {2,4,6,8,10,12};
    for (int k = 1; k <= 6; k++)
      {
	F2Obj.ConvBasis.insert({k,DISNCBasis{k}});
	F2Obj.C0.insert({k,Set<Operator>{F2Obj.ConvBasis.at(k), C2LO}});
	F2Obj.C1.insert({k,Set<Operator>{F2Obj.ConvBasis.at(k), C2NLO}});
	unordered_map<int,Set<Operator>> NNLO;
	for (int nf = nfi; nf <= nff; nf++)
	  NNLO.insert({nf,Set<Operator>{F2Obj.ConvBasis.at(k), C2NNLO.at(nf)}});
	F2Obj.C2.insert({k,NNLO});
      }
    t.stop();

    return F2Obj;
  }

  //_____________________________________________________________________________
  StructureFunctionObjects InitializeFLObjectsZM(Grid const& g, double const& IntEps)
  {
    cout << "Initializing StructureFunctionObjects for FL Zero Mass... ";
    Timer t;
    t.start();

    // Define object of the structure containing the DglapObjects
    StructureFunctionObjects FLObj;

    // Set initial and final number of active flavours. Precompute
    // objects for all numbers of flavours because it cheap enough.
    int nfi = 1;
    int nff = 6;

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps};

    // Coefficient functions for FL.
    // LO
    unordered_map<int,Operator> CLLO;
    CLLO.insert({DISNCBasis::CNS, Zero});
    CLLO.insert({DISNCBasis::CS,  Zero});
    CLLO.insert({DISNCBasis::CG,  Zero});

    // NLO
    unordered_map<int,Operator> CLNLO;
    const Operator OL1ns{g, CL1ns{}, IntEps};
    const Operator OL1g {g, CL1g{},  IntEps};
    CLNLO.insert({DISNCBasis::CNS, OL1ns});
    CLNLO.insert({DISNCBasis::CS,  OL1ns});
    CLNLO.insert({DISNCBasis::CG,  OL1g});

    // NNLO
    unordered_map<int,unordered_map<int,Operator>> CLNNLO;
    const Operator OL2ps{g, CL2ps{}, IntEps};
    const Operator OL2g {g, CL2g{},  IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
	const Operator OL2nsp{g, CL2nsp{nf}, IntEps};
	const Operator OL2t = OL2nsp + 6 * OL2ps;
	unordered_map<int,Operator> CLNNLOnf;
	CLNNLOnf.insert({DISNCBasis::CNS, OL2nsp});
	CLNNLOnf.insert({DISNCBasis::CS,  OL2t});
	CLNNLOnf.insert({DISNCBasis::CG,  OL2g});
	CLNNLO.insert({nf,CLNNLOnf});
      }

    // Allocate set of operators.
    FLObj.skip = {2,4,6,8,10,12};
    for (int k = 1; k <= 6; k++)
      {
	FLObj.ConvBasis.insert({k,DISNCBasis{k}});
	FLObj.C0.insert({k,Set<Operator>{FLObj.ConvBasis.at(k), CLLO}});
	FLObj.C1.insert({k,Set<Operator>{FLObj.ConvBasis.at(k), CLNLO}});
	unordered_map<int,Set<Operator>> NNLO;
	for (int nf = nfi; nf <= nff; nf++)
	  NNLO.insert({nf,Set<Operator>{FLObj.ConvBasis.at(k), CLNNLO.at(nf)}});
	FLObj.C2.insert({k,NNLO});
      }
    t.stop();

    return FLObj;
  }

  //_____________________________________________________________________________
  StructureFunctionObjects InitializeF3ObjectsZM(Grid const& g, double const& IntEps)
  {
    cout << "Initializing StructureFunctionObjects for F3 Zero Mass... ";
    Timer t;
    t.start();

    // Define object of the structure containing the DglapObjects
    StructureFunctionObjects F3Obj;

    // Set initial and final number of active flavours. Precompute
    // objects for all numbers of flavours because it cheap enough.
    int nfi = 1;
    int nff = 6;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // Coefficient functions for F3.
    // LO
    unordered_map<int,Operator> C3LO;
    C3LO.insert({DISNCBasis::CNS, Id});
    C3LO.insert({DISNCBasis::CS,  Id});
    C3LO.insert({DISNCBasis::CG,  Zero});

    // NLO
    unordered_map<int,Operator> C3NLO;
    const Operator O31ns{g, C31ns{}, IntEps};
    C3NLO.insert({DISNCBasis::CNS, O31ns});
    C3NLO.insert({DISNCBasis::CS,  O31ns});
    C3NLO.insert({DISNCBasis::CG,  Zero});

    // NNLO
    unordered_map<int,unordered_map<int,Operator>> C3NNLO;
    for (int nf = nfi; nf <= nff; nf++)
      {
	const Operator O32nsm{g, C32nsm{nf}, IntEps};
	const Operator O32t = O32nsm;
	unordered_map<int,Operator> C3NNLOnf;
	C3NNLOnf.insert({DISNCBasis::CNS, O32nsm});
	C3NNLOnf.insert({DISNCBasis::CS,  O32t});
	C3NNLOnf.insert({DISNCBasis::CG,  Zero});
	C3NNLO.insert({nf,C3NNLOnf});
      }

    // Allocate set of operators.
    F3Obj.skip = {1,3,5,7,9,11};
    for (int k = 1; k <= 6; k++)
      {
	F3Obj.ConvBasis.insert({k,DISNCBasis{k}});
	F3Obj.C0.insert({k,Set<Operator>{F3Obj.ConvBasis.at(k), C3LO}});
	F3Obj.C1.insert({k,Set<Operator>{F3Obj.ConvBasis.at(k), C3NLO}});
	unordered_map<int,Set<Operator>> NNLO;
	for (int nf = nfi; nf <= nff; nf++)
	  NNLO.insert({nf,Set<Operator>{F3Obj.ConvBasis.at(k), C3NNLO.at(nf)}});
	F3Obj.C2.insert({k,NNLO});
      }
    t.stop();

    return F3Obj;
  }

  //_____________________________________________________________________________
  unordered_map<int,Observable> StructureFunctionBuildNC(StructureFunctionObjects                                          const& FObj,
							 function<unordered_map<int,double>(double const&, double const&)> const& InDistFunc,
							 vector<double>                                                    const& Thresholds,
							 int                                                               const& PerturbativeOrder,
							 function<double(double const&)>                                   const& Alphas,
							 function<vector<double>(double const&)>                           const& Charges)
  {
    // Get grid.
    Grid const& g = FObj.C0.at(1).at(0).GetGrid();

    // Compute observables for structure function.
    unordered_map<int,Observable> F;
    for (int k = 1; k <= 6; k++)
      {
	// Define coefficient function functions.
	function<Set<Operator>(double const&)> Cf;
	if (PerturbativeOrder == 0)
	  {
	    Cf = [=] (double const& Q) -> Set<Operator>
	      {
		return Charges(Q)[k-1] * FObj.C0.at(k);
	      };
	  }
	else if (PerturbativeOrder == 1)
	  {
	    Cf = [=] (double const& Q) -> Set<Operator>
	      {
		const auto cp = Alphas(Q) / FourPi;
		return Charges(Q)[k-1] * ( FObj.C0.at(k) + cp * FObj.C1.at(k) );
	      };
	  }
	else if (PerturbativeOrder == 2)
	  {
	    Cf = [=] (double const& Q) -> Set<Operator>
	      {
		const auto cp = Alphas(Q) / FourPi;
		const auto nf = NF(Q, Thresholds);
		return Charges(Q)[k-1] * ( FObj.C0.at(k) + cp * ( FObj.C1.at(k) +  + cp * FObj.C2.at(k).at(nf) ) );
	      };
	  }

	// Define distribution function functions.
	const auto DistF = [=,&g] (double const& Q) -> Set<Distribution>
	  {
	    return Set<Distribution>{FObj.ConvBasis.at(k), DistributionMap(g, InDistFunc, Q, FObj.skip)};
	  };
	// Initialize "Observable".
	F.insert({k,Observable{Cf, DistF}});
      }

    // Total structure function.
    function<Set<Operator>(double const&)> Cf;
    if (PerturbativeOrder == 0)
      {
	Cf = [=] (double const& Q) -> Set<Operator>
	  {
	    const DISNCBasis basis{Charges(Q)};
	    Set<Operator> LO{basis, FObj.C0.at(1).GetObjects()};
	    return LO;
	  };
      }
    else if (PerturbativeOrder == 1)
      {
	Cf = [=] (double const& Q) -> Set<Operator>
	  {
	    const auto cp = Alphas(Q) / FourPi;
	    const DISNCBasis basis{Charges(Q)};
	    Set<Operator> LO {basis, FObj.C0.at(1).GetObjects()};
	    Set<Operator> NLO{basis, FObj.C1.at(1).GetObjects()};
	    return LO + cp * NLO;
	  };
      }
    else if (PerturbativeOrder == 2)
      {
	Cf = [=] (double const& Q) -> Set<Operator>
	  {
	    const auto cp = Alphas(Q) / FourPi;
	    const auto nf = NF(Q, Thresholds);
	    const DISNCBasis basis{Charges(Q)};
	    Set<Operator> LO  {basis, FObj.C0.at(1).GetObjects()};
	    Set<Operator> NLO {basis, FObj.C1.at(1).GetObjects()};
	    Set<Operator> NNLO{basis, FObj.C2.at(1).at(nf).GetObjects()};
	    return LO + cp * ( NLO +  + cp * NNLO );
	  };
      }

    // Define distribution function functions.
    const auto DistF = [=,&g] (double const& Q) -> Set<Distribution>
      {
	return Set<Distribution>{DISNCBasis{Charges(Q)}, DistributionMap(g, InDistFunc, Q, FObj.skip)};
      };
    F.insert({0,Observable{Cf, DistF}});

    return F;
  }

  //_____________________________________________________________________________
  unordered_map<int,Observable> StructureFunctionBuildNC(StructureFunctionObjects                                   const& FObj,
							 function<double(int const&, double const&, double const&)> const& InDistFunc,
							 vector<double>                                             const& Thresholds,
							 int                                                        const& PerturbativeOrder,
							 function<double(double const&)>                            const& Alphas,
							 function<vector<double>(double const&)>                    const& Charges)
  {
    const auto InDistFuncMap = [=] (double const& x, double const& Q) -> unordered_map<int,double>
      {
	unordered_map<int,double> DistMap;
	for (int i = 0; i <= 12; i++)
	  DistMap.insert({i,InDistFunc(i, x, Q)});
	return DistMap;
      };
    return StructureFunctionBuildNC(FObj, InDistFuncMap, Thresholds, PerturbativeOrder, Alphas, Charges);
  }


}
