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
#include "apfel/tools.h"
#include "apfel/constants.h"
#include "apfel/zeromasscoefficientfunctions.h"
#include "apfel/massivecoefficientfunctions.h"
#include "apfel/massivezerocoefficientfunctions.h"
#include "apfel/tabulateobject.h"

using namespace std;

namespace apfel {
  //_____________________________________________________________________________
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeF2NCObjectsZM(Grid           const& g,
												   vector<double> const& Thresholds,
												   double         const& IntEps)
  {
    cout << "Initializing StructureFunctionObjects for F2 NC Zero Mass... ";
    Timer t;
    t.start();

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // LO
    map<int,Operator> C2LO;
    C2LO.insert({CNS, Id});
    C2LO.insert({CS,  Id});
    C2LO.insert({CG,  Zero});

    // NLO
    map<int,Operator> C2NLO;
    const Operator O21ns{g, C21ns{}, IntEps};
    const Operator O21g {g, C21g{},  IntEps};
    C2NLO.insert({CNS, O21ns});
    C2NLO.insert({CS,  O21ns});
    C2NLO.insert({CG,  O21g});

    // NNLO
    map<int,map<int,Operator>> C2NNLO;
    const Operator O22ps{g, C22ps{}, IntEps};
    const Operator O22g {g, C22g{},  IntEps};
    for (int nf = 1; nf <= 6; nf++)
      {
	const Operator O22nsp{g, C22nsp{nf}, IntEps};
	const Operator O22t = O22nsp + 6 * O22ps;
	map<int,Operator> C2NNLOnf;
	C2NNLOnf.insert({CNS, O22nsp});
	C2NNLOnf.insert({CS,  O22t});
	C2NNLOnf.insert({CG,  O22g});
	C2NNLO.insert({nf,C2NNLOnf});
      }

    // Vector of distributions to skip
    const vector<int> skip = {2,4,6,8,10,12};

    // Define object of the structure containing the DglapObjects
    const auto F2Obj = [=] (double const& Q, vector<double> const& Ch) -> StructureFunctionObjects
      {
	// Determine number of active flavours.
	const int nf = NF(Q, Thresholds);

	// Effective charges. The charges of the components with mh >
	// Q are set to zero.
	vector<double> EffCh;
	for (int k = 1; k <= 6; k++)
	  EffCh.push_back(( k > nf ? 0 : Ch[k-1]));

	// Fill in structure function object
	StructureFunctionObjects FObj;
	FObj.skip = skip;
	// Single structure function components.
	for (int k = 0; k <= 6; k++)
	  {
	    FObj.ConvBasis.insert({k,(k == 0 ? DISNCBasis{EffCh} : DISNCBasis{k, EffCh[k-1]})});
	    FObj.C0.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C2LO}});
	    FObj.C1.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C2NLO}});
	    FObj.C2.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C2NNLO.at(nf)}});
	  }
	return FObj;
      };
    t.stop();

    return F2Obj;
  }

  //_____________________________________________________________________________
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeFLNCObjectsZM(Grid           const& g,
												   vector<double> const& Thresholds,
												   double         const& IntEps)
  {
    cout << "Initializing StructureFunctionObjects for FL NC Zero Mass... ";
    Timer t;
    t.start();

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps};

    // LO
    map<int,Operator> CLLO;
    CLLO.insert({CNS, Zero});
    CLLO.insert({CS,  Zero});
    CLLO.insert({CG,  Zero});

    // NLO
    map<int,Operator> CLNLO;
    const Operator OL1ns{g, CL1ns{}, IntEps};
    const Operator OL1g {g, CL1g{},  IntEps};
    CLNLO.insert({CNS, OL1ns});
    CLNLO.insert({CS,  OL1ns});
    CLNLO.insert({CG,  OL1g});

    // NNLO (for nf from 1 to 6)
    map<int,map<int,Operator>> CLNNLO;
    const Operator OL2ps{g, CL2ps{}, IntEps};
    const Operator OL2g {g, CL2g{},  IntEps};
    for (int nf = 1; nf <= 6; nf++)
      {
	const Operator OL2nsp{g, CL2nsp{nf}, IntEps};
	const Operator OL2t = OL2nsp + 6 * OL2ps;
	map<int,Operator> CLNNLOnf;
	CLNNLOnf.insert({CNS, OL2nsp});
	CLNNLOnf.insert({CS,  OL2t});
	CLNNLOnf.insert({CG,  OL2g});
	CLNNLO.insert({nf,CLNNLOnf});
      }

    // Vector of distributions to skip
    const vector<int> skip = {2,4,6,8,10,12};

    // Define object of the structure containing the DglapObjects
    const auto FLObj = [=] (double const& Q, vector<double> const& Ch) -> StructureFunctionObjects
      {
	// Determine number of active flavours.
	const int nf = NF(Q, Thresholds);

	// Effective charges. The charges of the components with mh >
	// Q are set to zero.
	vector<double> EffCh;
	for (int k = 1; k <= 6; k++)
	  EffCh.push_back(( k > nf ? 0 : Ch[k-1]));

	// Fill in structure function object
	StructureFunctionObjects FObj;
	FObj.skip = skip;
	// Single structure function components.
	for (int k = 0; k <= 6; k++)
	  {
	    FObj.ConvBasis.insert({k,(k == 0 ? DISNCBasis{EffCh} : DISNCBasis{k, EffCh[k-1]})});
	    FObj.C0.insert({k,Set<Operator>{FObj.ConvBasis.at(k), CLLO}});
	    FObj.C1.insert({k,Set<Operator>{FObj.ConvBasis.at(k), CLNLO}});
	    FObj.C2.insert({k,Set<Operator>{FObj.ConvBasis.at(k), CLNNLO.at(nf)}});
	  }
	return FObj;
      };
    t.stop();

    return FLObj;
  }

  //_____________________________________________________________________________
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeF3NCObjectsZM(Grid           const& g,
												   vector<double> const& Thresholds,
												   double         const& IntEps)
  {
    cout << "Initializing StructureFunctionObjects for F3 NC Zero Mass... ";
    Timer t;
    t.start();

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // LO
    map<int,Operator> C3LO;
    C3LO.insert({CNS, Id});
    C3LO.insert({CS,  Id});
    C3LO.insert({CG,  Zero});

    // NLO
    map<int,Operator> C3NLO;
    const Operator O31ns{g, C31ns{}, IntEps};
    C3NLO.insert({CNS, O31ns});
    C3NLO.insert({CS,  O31ns});
    C3NLO.insert({CG,  Zero});

    // NNLO
    map<int,map<int,Operator>> C3NNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
	const Operator O32nsm{g, C32nsm{nf}, IntEps};
	const Operator O32t = O32nsm;
	map<int,Operator> C3NNLOnf;
	C3NNLOnf.insert({CNS, O32nsm});
	C3NNLOnf.insert({CS,  O32t});
	C3NNLOnf.insert({CG,  Zero});
	C3NNLO.insert({nf,C3NNLOnf});
      }

    // Vector of distributions to skip
    const vector<int> skip = {1,3,5,7,9,11};

    // Define object of the structure containing the DglapObjects
    const auto F3Obj = [=] (double const& Q, vector<double> const& Ch) -> StructureFunctionObjects
      {
	// Determine number of active flavours.
	const int nf = NF(Q, Thresholds);

	// Effective charges. The charges of the components with mh >
	// Q are set to zero.
	vector<double> EffCh;
	for (int k = 1; k <= 6; k++)
	  EffCh.push_back(( k > nf ? 0 : Ch[k-1]));

	// Fill in structure function object
	StructureFunctionObjects FObj;
	FObj.skip = skip;
	// Single structure function components.
	for (int k = 0; k <= 6; k++)
	  {
	    FObj.ConvBasis.insert({k,(k == 0 ? DISNCBasis{EffCh} : DISNCBasis{k, EffCh[k-1]})});
	    FObj.C0.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C3LO}});
	    FObj.C1.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C3NLO}});
	    FObj.C2.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C3NNLO.at(nf)}});
	  }
	return FObj;
      };
    t.stop();

    return F3Obj;
  }

  //_____________________________________________________________________________
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeF2CCPlusObjectsZM(Grid           const& g,
												       vector<double> const& Thresholds,
												       double         const& IntEps)
  {
    cout << "Initializing StructureFunctionObjects for ( F2(nu) + F2(nubar) ) / 2 Zero Mass... ";
    Timer t;
    t.start();

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // LO
    map<int,Operator> C2LO;
    C2LO.insert({CNS, Id});
    C2LO.insert({CS,  Id});
    C2LO.insert({CG,  Zero});

    // NLO
    map<int,Operator> C2NLO;
    const Operator O21ns{g, C21ns{}, IntEps};
    const Operator O21g {g, C21g{},  IntEps};
    C2NLO.insert({CNS, O21ns});
    C2NLO.insert({CS,  O21ns});
    C2NLO.insert({CG,  O21g});

    // NNLO
    map<int,map<int,Operator>> C2NNLO;
    const Operator O22ps{g, C22ps{}, IntEps};
    const Operator O22g {g, C22g{},  IntEps};
    for (int nf = 1; nf <= 6; nf++)
      {
	const Operator O22nsp{g, C22nsp{nf}, IntEps};
	const Operator O22t = O22nsp + 6 * O22ps;
	map<int,Operator> C2NNLOnf;
	C2NNLOnf.insert({CNS, O22nsp});
	C2NNLOnf.insert({CS,  O22t});
	C2NNLOnf.insert({CG,  O22g});
	C2NNLO.insert({nf,C2NNLOnf});
      }

    // Vector of distributions to skip
    const vector<int> skip = {2,4,6,8,10,12};

    // Define object of the structure containing the DglapObjects
    const auto F2Obj = [=] (double const& Q, vector<double> const& Ch) -> StructureFunctionObjects
      {
	// Determine number of active flavours.
	const int nf = NF(Q, Thresholds);

	// Effective charges.
	vector<double> EffCh;
	if (nf <= 3)
	  EffCh = {Ch[0], Ch[1], 0, 0, 0, 0, 0, 0, 0};
	else if (nf == 4)
	  EffCh = {Ch[0], Ch[1], 0, Ch[3], Ch[4], 0, 0, 0, 0};
	else if (nf == 5)
	  EffCh = {Ch[0], Ch[1], Ch[2], Ch[3], Ch[4], Ch[5], 0, 0, 0};
	else
	  EffCh = Ch;

	// Fill in structure function object
	StructureFunctionObjects FObj;
	FObj.skip = skip;
	// Single structure function components.
	for (int k = 0; k <= 6; k++)
	  {
	    FObj.ConvBasis.insert({k,(k == 0 ? DISCCBasis{EffCh, false} : DISCCBasis{k, false, EffCh[k-1]})});
	    FObj.C0.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C2LO}});
	    FObj.C1.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C2NLO}});
	    FObj.C2.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C2NNLO.at(nf)}});
	  }
	return FObj;
      };
    t.stop();

    return F2Obj;
  }

  //_____________________________________________________________________________
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeF2CCMinusObjectsZM(Grid           const& g,
													vector<double> const& Thresholds,
													double         const& IntEps)
  {
    cout << "Initializing StructureFunctionObjects for ( F2(nu) - F2(nubar) ) / 2 Zero Mass... ";
    Timer t;
    t.start();

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // LO
    map<int,Operator> C2LO;
    C2LO.insert({CNS, Id});
    C2LO.insert({CS,  Zero});
    C2LO.insert({CG,  Zero});

    // NLO
    map<int,Operator> C2NLO;
    const Operator O21ns{g, C21ns{}, IntEps};
    C2NLO.insert({CNS, O21ns});
    C2NLO.insert({CS,  Zero});
    C2NLO.insert({CG,  Zero});

    // NNLO
    map<int,map<int,Operator>> C2NNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
	const Operator O22nsm{g, C22nsm{nf}, IntEps};
	map<int,Operator> C2NNLOnf;
	C2NNLOnf.insert({CNS, O22nsm});
	C2NNLOnf.insert({CS,  Zero});
	C2NNLOnf.insert({CG,  Zero});
	C2NNLO.insert({nf,C2NNLOnf});
      }

    // Vector of distributions to skip
    const vector<int> skip = {0,1,2,3,5,7,9,11};

    // Define object of the structure containing the DglapObjects
    const auto F2Obj = [=] (double const& Q, vector<double> const& Ch) -> StructureFunctionObjects
      {
	// Determine number of active flavours.
	const int nf = NF(Q, Thresholds);

	// Effective charges.
	vector<double> EffCh;
	if (nf <= 3)
	  EffCh = {Ch[0], Ch[1], 0, 0, 0, 0, 0, 0, 0};
	else if (nf == 4)
	  EffCh = {Ch[0], Ch[1], 0, Ch[3], Ch[4], 0, 0, 0, 0};
	else if (nf == 5)
	  EffCh = {Ch[0], Ch[1], Ch[2], Ch[3], Ch[4], Ch[5], 0, 0, 0};
	else
	  EffCh = Ch;

	// Fill in structure function object
	StructureFunctionObjects FObj;
	FObj.skip = skip;
	// Single structure function components.
	for (int k = 0; k <= 6; k++)
	  {
	    FObj.ConvBasis.insert({k,(k == 0 ? DISCCBasis{EffCh, false} : DISCCBasis{k, false, EffCh[k-1]})});
	    FObj.C0.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C2LO}});
	    FObj.C1.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C2NLO}});
	    FObj.C2.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C2NNLO.at(nf)}});
	  }
	return FObj;
      };
    t.stop();

    return F2Obj;
  }

  //_____________________________________________________________________________
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeFLCCPlusObjectsZM(Grid           const& g,
												       vector<double> const& Thresholds,
												       double         const& IntEps)
  {
    cout << "Initializing StructureFunctionObjects for ( FL(nu) + FL(nubar) ) / 2 Zero Mass... ";
    Timer t;
    t.start();

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps};

    // LO
    map<int,Operator> CLLO;
    CLLO.insert({CNS, Zero});
    CLLO.insert({CS,  Zero});
    CLLO.insert({CG,  Zero});

    // NLO
    map<int,Operator> CLNLO;
    const Operator OL1ns{g, CL1ns{}, IntEps};
    const Operator OL1g {g, CL1g{},  IntEps};
    CLNLO.insert({CNS, OL1ns});
    CLNLO.insert({CS,  OL1ns});
    CLNLO.insert({CG,  OL1g});

    // NNLO
    map<int,map<int,Operator>> CLNNLO;
    const Operator OL2ps{g, CL2ps{}, IntEps};
    const Operator OL2g {g, CL2g{},  IntEps};
    for (int nf = 1; nf <= 6; nf++)
      {
	const Operator OL2nsp{g, CL2nsp{nf}, IntEps};
	const Operator OL2t = OL2nsp + 6 * OL2ps;
	map<int,Operator> CLNNLOnf;
	CLNNLOnf.insert({CNS, OL2nsp});
	CLNNLOnf.insert({CS,  OL2t});
	CLNNLOnf.insert({CG,  OL2g});
	CLNNLO.insert({nf,CLNNLOnf});
      }

    // Vector of distributions to skip
    const vector<int> skip = {2,4,6,8,10,12};

    // Define object of the structure containing the DglapObjects
    const auto FLObj = [=] (double const& Q, vector<double> const& Ch) -> StructureFunctionObjects
      {
	// Determine number of active flavours.
	const int nf = NF(Q, Thresholds);

	// Effective charges.
	vector<double> EffCh;
	if (nf <= 3)
	  EffCh = {Ch[0], Ch[1], 0, 0, 0, 0, 0, 0, 0};
	else if (nf == 4)
	  EffCh = {Ch[0], Ch[1], 0, Ch[3], Ch[4], 0, 0, 0, 0};
	else if (nf == 5)
	  EffCh = {Ch[0], Ch[1], Ch[2], Ch[3], Ch[4], Ch[5], 0, 0, 0};
	else
	  EffCh = Ch;

	// Fill in structure function object
	StructureFunctionObjects FObj;
	FObj.skip = skip;
	// Single structure function components.
	for (int k = 0; k <= 6; k++)
	  {
	    FObj.ConvBasis.insert({k,(k == 0 ? DISCCBasis{EffCh, false} : DISCCBasis{k, false, EffCh[k-1]})});
	    FObj.C0.insert({k,Set<Operator>{FObj.ConvBasis.at(k), CLLO}});
	    FObj.C1.insert({k,Set<Operator>{FObj.ConvBasis.at(k), CLNLO}});
	    FObj.C2.insert({k,Set<Operator>{FObj.ConvBasis.at(k), CLNNLO.at(nf)}});
	  }
	return FObj;
      };
    t.stop();

    return FLObj;
  }

  //_____________________________________________________________________________
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeFLCCMinusObjectsZM(Grid           const& g,
												       vector<double> const& Thresholds,
												       double         const& IntEps)
  {
    cout << "Initializing StructureFunctionObjects for ( FL(nu) - FL(nubar) ) / 2 Zero Mass... ";
    Timer t;
    t.start();

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps};

    // LO
    map<int,Operator> CLLO;
    CLLO.insert({CNS, Zero});
    CLLO.insert({CS,  Zero});
    CLLO.insert({CG,  Zero});

    // NLO
    map<int,Operator> CLNLO;
    const Operator OL1ns{g, CL1ns{}, IntEps};
    CLNLO.insert({CNS, OL1ns});
    CLNLO.insert({CS,  Zero});
    CLNLO.insert({CG,  Zero});

    // NNLO
    map<int,map<int,Operator>> CLNNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
	const Operator OL2nsm{g, CL2nsm{nf}, IntEps};
	map<int,Operator> CLNNLOnf;
	CLNNLOnf.insert({CNS, OL2nsm});
	CLNNLOnf.insert({CS,  Zero});
	CLNNLOnf.insert({CG,  Zero});
	CLNNLO.insert({nf,CLNNLOnf});
      }

    // Vector of distributions to skip
    const vector<int> skip = {0,1,2,3,5,7,9,11};

    // Define object of the structure containing the DglapObjects
    const auto FLObj = [=] (double const& Q, vector<double> const& Ch) -> StructureFunctionObjects
      {
	// Determine number of active flavours.
	const int nf = NF(Q, Thresholds);

	// Effective charges.
	vector<double> EffCh;
	if (nf <= 3)
	  EffCh = {Ch[0], Ch[1], 0, 0, 0, 0, 0, 0, 0};
	else if (nf == 4)
	  EffCh = {Ch[0], Ch[1], 0, Ch[3], Ch[4], 0, 0, 0, 0};
	else if (nf == 5)
	  EffCh = {Ch[0], Ch[1], Ch[2], Ch[3], Ch[4], Ch[5], 0, 0, 0};
	else
	  EffCh = Ch;

	// Fill in structure function object
	StructureFunctionObjects FObj;
	FObj.skip = skip;
	// Single structure function components.
	for (int k = 0; k <= 6; k++)
	  {
	    FObj.ConvBasis.insert({k,(k == 0 ? DISCCBasis{EffCh, false} : DISCCBasis{k, false, EffCh[k-1]})});
	    FObj.C0.insert({k,Set<Operator>{FObj.ConvBasis.at(k), CLLO}});
	    FObj.C1.insert({k,Set<Operator>{FObj.ConvBasis.at(k), CLNLO}});
	    FObj.C2.insert({k,Set<Operator>{FObj.ConvBasis.at(k), CLNNLO.at(nf)}});
	  }
	return FObj;
      };
    t.stop();

    return FLObj;
  }

  //_____________________________________________________________________________
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeF3CCPlusObjectsZM(Grid           const& g,
												       vector<double> const& Thresholds,
												       double         const& IntEps)
  {
    cout << "Initializing StructureFunctionObjects for ( F3(nu) + F3(nubar) ) / 2 Zero Mass... ";
    Timer t;
    t.start();

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // LO
    map<int,Operator> C3LO;
    C3LO.insert({CNS, Id});
    C3LO.insert({CS,  Id});
    C3LO.insert({CG,  Zero});

    // NLO
    map<int,Operator> C3NLO;
    const Operator O31ns{g, C31ns{}, IntEps};
    C3NLO.insert({CNS, O31ns});
    C3NLO.insert({CS,  O31ns});
    C3NLO.insert({CG,  Zero});

    // NNLO
    map<int,map<int,Operator>> C3NNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
	const Operator O32nsp{g, C32nsp{nf}, IntEps};
	const Operator O32t = O32nsp;
	map<int,Operator> C3NNLOnf;
	C3NNLOnf.insert({CNS, O32nsp});
	C3NNLOnf.insert({CS,  O32t});
	C3NNLOnf.insert({CG,  Zero});
	C3NNLO.insert({nf,C3NNLOnf});
      }

    // Vector of distributions to skip
    const vector<int> skip = {0,1,2,4,6,8,10,12};

    // Define object of the structure containing the DglapObjects
    const auto F3Obj = [=] (double const& Q, vector<double> const& Ch) -> StructureFunctionObjects
      {
	// Determine number of active flavours.
	const int nf = NF(Q, Thresholds);

	// Effective charges.
	vector<double> EffCh;
	if (nf <= 3)
	  EffCh = {Ch[0], Ch[1], 0, 0, 0, 0, 0, 0, 0};
	else if (nf == 4)
	  EffCh = {Ch[0], Ch[1], 0, Ch[3], Ch[4], 0, 0, 0, 0};
	else if (nf == 5)
	  EffCh = {Ch[0], Ch[1], Ch[2], Ch[3], Ch[4], Ch[5], 0, 0, 0};
	else
	  EffCh = Ch;

	// Fill in structure function object
	StructureFunctionObjects FObj;
	FObj.skip = skip;
	// Single structure function components.
	for (int k = 0; k <= 6; k++)
	  {
	    FObj.ConvBasis.insert({k,(k == 0 ? DISCCBasis{EffCh, true} : DISCCBasis{k, true, EffCh[k-1]})});
	    FObj.C0.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C3LO}});
	    FObj.C1.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C3NLO}});
	    FObj.C2.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C3NNLO.at(nf)}});
	  }
	return FObj;
      };
    t.stop();

    return F3Obj;
  }

  //_____________________________________________________________________________
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeF3CCMinusObjectsZM(Grid           const& g,
												       vector<double> const& Thresholds,
												       double         const& IntEps)
  {
    cout << "Initializing StructureFunctionObjects for ( F3(nu) - F3(nubar) ) / 2 Zero Mass... ";
    Timer t;
    t.start();

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // LO
    map<int,Operator> C3LO;
    C3LO.insert({CNS, Id});
    C3LO.insert({CS,  Id});
    C3LO.insert({CG,  Zero});

    // NLO
    map<int,Operator> C3NLO;
    const Operator O31ns{g, C31ns{}, IntEps};
    C3NLO.insert({CNS, O31ns});
    C3NLO.insert({CS,  O31ns});
    C3NLO.insert({CG,  Zero});

    // NNLO
    map<int,map<int,Operator>> C3NNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
	const Operator O32nsm{g, C32nsm{nf}, IntEps};
	map<int,Operator> C3NNLOnf;
	C3NNLOnf.insert({CNS, O32nsm});
	C3NNLOnf.insert({CS,  O32nsm});
	C3NNLOnf.insert({CG,  Zero});
	C3NNLO.insert({nf,C3NNLOnf});
      }

    // Vector of distributions to skip
    const vector<int> skip = {0,1,3,5,7,9,11};

    // Define object of the structure containing the DglapObjects
    const auto F3Obj = [=] (double const& Q, vector<double> const& Ch) -> StructureFunctionObjects
      {
	// Determine number of active flavours.
	const int nf = NF(Q, Thresholds);

	// Effective charges.
	vector<double> EffCh;
	if (nf <= 3)
	  EffCh = {Ch[0], Ch[1], 0, 0, 0, 0, 0, 0, 0};
	else if (nf == 4)
	  EffCh = {Ch[0], Ch[1], 0, Ch[3], Ch[4], 0, 0, 0, 0};
	else if (nf == 5)
	  EffCh = {Ch[0], Ch[1], Ch[2], Ch[3], Ch[4], Ch[5], 0, 0, 0};
	else
	  EffCh = Ch;

	// Fill in structure function object
	StructureFunctionObjects FObj;
	FObj.skip = skip;
	// Single structure function components.
	for (int k = 0; k <= 6; k++)
	  {
	    FObj.ConvBasis.insert({k,(k == 0 ? DISCCBasis{EffCh, true} : DISCCBasis{k, true, EffCh[k-1]})});
	    FObj.C0.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C3LO}});
	    FObj.C1.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C3NLO}});
	    FObj.C2.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C3NNLO.at(nf)}});
	  }
	return FObj;
      };
    t.stop();

    return F3Obj;
  }

  //_____________________________________________________________________________
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeF2NCObjectsMassive(Grid           const& g,
													vector<double> const& Masses,
													double         const& IntEps,
													int            const& nxi,
													double         const& ximin,
													double         const& ximax,
													int            const& intdeg,
													double         const& lambda)
  {
    Timer t;
    t.start();

    // Determine number of active flavours
    int actnf = 0;
    for (auto const& m : Masses)
      if (m < eps8)
	actnf++;

    cout << "Initializing StructureFunctionObjects for F2 NC Massive with " << actnf << " active flavours... \n";

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // Zero Mass coefficient functions
    // LO
    map<int,Operator> C2LO;
    C2LO.insert({CNS, Id});
    C2LO.insert({CS,  Id});
    C2LO.insert({CG,  Zero});

    // NLO
    map<int,Operator> C2NLO;
    const Operator O21ns{g, C21ns{}, IntEps};
    const Operator O21g {g, C21g{},  IntEps};
    C2NLO.insert({CNS, O21ns});
    C2NLO.insert({CS,  O21ns});
    C2NLO.insert({CG,  O21g});

    // NNLO
    const Operator O22ps {g, C22ps{}, IntEps};
    const Operator O22g  {g, C22g{},  IntEps};
    const Operator O22nsp{g, C22nsp{actnf}, IntEps};
    const Operator O22t = O22nsp + 6 * O22ps;

    // Massive coefficient functions
    // Null set of operator needed for the LO coefficient functions.
    map<int,Operator> Not;
    Not.insert({CNS, Zero});
    Not.insert({CS,  Zero});
    Not.insert({CG,  Zero});

    // NLO
    const TabulateObject<Operator> TabO21g{[=,&g] (double const& xi) -> Operator
	{
	  const double teta = 1 / ( 1 + 4 / xi );
	  return Operator{g, Cm21gNC{teta}, IntEps};
	}, nxi, ximin, ximax, intdeg, {}, lambda};

    // NNLO
    const TabulateObject<Operator> TabO22ns{[=,&g] (double const& xi) -> Operator
	{
	  const double teta = 1 / ( 1 + 4 / xi );
	  return Operator{g, Cm22nsNC{teta}, IntEps};
	}, nxi, ximin, ximax, intdeg, {}, lambda};
    const auto fO22s = [=,&g] (double const& xi) -> Operator
      {
	const double teta = 1 / ( 1 + 4 / xi );
	const Operator O22psc{g, Cm22psNC{teta}, IntEps};
	const Operator O22psl{g, Cm22barpsNC{teta}, IntEps};
	return 6 * ( O22psc + log(xi) * O22psl );
      };
    const TabulateObject<Operator> TabO22s{fO22s, nxi, ximin, ximax, intdeg, {}, lambda};
    const auto fO22g = [=,&g] (double const& xi) -> Operator
      {
	const double teta = 1 / ( 1 + 4 / xi );
	const Operator O22gc{g, Cm22gNC{teta}, IntEps};
	const Operator O22gl{g, Cm22bargNC{teta}, IntEps};
	return O22gc + log(xi) * O22gl;
      };
    const TabulateObject<Operator> TabO22g{fO22g, nxi, ximin, ximax, intdeg, {}, lambda};

    // Vector of distributions to skip
    const vector<int> skip = {2,4,6,8,10,12};

    // Define object of the structure containing the DglapObjects
    const auto F2Obj = [=] (double const& Q, vector<double> const& Ch) -> StructureFunctionObjects
      {
	// Fill in structure function object
	StructureFunctionObjects FObj;
	FObj.skip = skip;

	// Start with the heavy quark structure function components.
	const double Q2  = Q * Q;
	Operator lNNLOns = O22nsp;
	for (int k = actnf + 1; k <= 6; k++)
	  {
	    // Compute value of xi.
	    const double M2 = Masses[k-1] * Masses[k-1];
	    const double xi = Q2 / M2;

	    // Convolution Basis
	    FObj.ConvBasis.insert({k,DISNCBasis{k, Ch[k-1]}});

	    // Set LO non-singlet and singlet coefficient
	    // functions to zero (gluon is already zero).
	    FObj.C0.insert({k,Set<Operator>{FObj.ConvBasis.at(k),Not}});

	    // Now insert NLO
	    map<int,Operator> NLO;
	    NLO.insert({CNS, Zero});
	    NLO.insert({CS,  Zero});
	    NLO.insert({CG,  TabO21g.Evaluate(xi)});
	    FObj.C1.insert({k,Set<Operator>{FObj.ConvBasis.at(k),NLO}});

	    // Now insert NNLO
	    map<int,Operator> NNLO;
	    NNLO.insert({CNS, Zero});
	    NNLO.insert({CS,  TabO22s.Evaluate(xi)});
	    NNLO.insert({CG,  TabO22g.Evaluate(xi)});
	    FObj.C2.insert({k,Set<Operator>{FObj.ConvBasis.at(k),NNLO}});

	    // Include the massive bit to the non-singlet coefficient
	    // function needed for the light structure functions.
	    for (int i = 0; i <= actnf; i++)
	      lNNLOns += TabO22ns.Evaluate(xi);
	  }

	// Now fill in the light components. To be updated in the
	// loop over the heavy component.
	map<int,Operator> C2NNLO;
	C2NNLO.insert({CNS, lNNLOns});
	C2NNLO.insert({CS,  O22t});
	C2NNLO.insert({CG,  O22g});
	for (int k = 1; k <= actnf; k++)
	  {
	    // Convolution Basis
	    FObj.ConvBasis.insert({k,DISNCBasis{k, Ch[k-1]}});
	    FObj.C0.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C2LO}});
	    FObj.C1.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C2NLO}});
	    FObj.C2.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C2NNLO}});
	  }

	// Total structure function set to zero for now. Need to
	// construct it a posteriori as a sum of light and heavy
	// components.
	FObj.ConvBasis.insert({0, DISNCBasis{Ch}});
	FObj.C0.insert({0,Set<Operator>{FObj.ConvBasis.at(0), Not}});
	FObj.C1.insert({0,Set<Operator>{FObj.ConvBasis.at(0), Not}});
	FObj.C2.insert({0,Set<Operator>{FObj.ConvBasis.at(0), Not}});
	return FObj;
      };
    t.stop();

    return F2Obj;
  }

  //_____________________________________________________________________________
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeFLNCObjectsMassive(Grid           const& g,
													vector<double> const& Masses,
													double         const& IntEps,
													int            const& nxi,
													double         const& ximin,
													double         const& ximax,
													int            const& intdeg,
													double         const& lambda)
  {
    Timer t;
    t.start();

    // Determine number of active flavours
    int actnf = 0;
    for (auto const& m : Masses)
      if (m < eps8)
	actnf++;

    cout << "Initializing StructureFunctionObjects for FL NC Massive with " << actnf << " active flavours... \n";

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps};

    // Zero Mass coefficient functions
    // LO
    map<int,Operator> CLLO;
    CLLO.insert({CNS, Zero});
    CLLO.insert({CS,  Zero});
    CLLO.insert({CG,  Zero});

    // NLO
    map<int,Operator> CLNLO;
    const Operator OL1ns{g, CL1ns{}, IntEps};
    const Operator OL1g {g, CL1g{},  IntEps};
    CLNLO.insert({CNS, OL1ns});
    CLNLO.insert({CS,  OL1ns});
    CLNLO.insert({CG,  OL1g});

    // NNLO
    const Operator OL2ps {g, CL2ps{}, IntEps};
    const Operator OL2g  {g, CL2g{},  IntEps};
    const Operator OL2nsp{g, CL2nsp{actnf}, IntEps};
    const Operator OL2t = OL2nsp + 6 * OL2ps;

    // Massive coefficient functions
    // Null set of operator needed for the LO coefficient functions.
    map<int,Operator> Not;
    Not.insert({CNS, Zero});
    Not.insert({CS,  Zero});
    Not.insert({CG,  Zero});

    // NLO
    const TabulateObject<Operator> TabOL1g{[=,&g] (double const& xi) -> Operator
	{
	  const double teta = 1 / ( 1 + 4 / xi );
	  return Operator{g, CmL1gNC{teta}, IntEps};
	}, nxi, ximin, ximax, intdeg, {}, lambda};

    // NNLO
    const TabulateObject<Operator> TabOL2ns{[=,&g] (double const& xi) -> Operator
	{
	  const double teta = 1 / ( 1 + 4 / xi );
	  return Operator{g, CmL2nsNC{teta}, IntEps};
	}, nxi, ximin, ximax, intdeg, {}, lambda};
    const auto fOL2s = [=,&g] (double const& xi) -> Operator
      {
	const double teta = 1 / ( 1 + 4 / xi );
	const Operator OL2psc{g, CmL2psNC{teta}, IntEps};
	const Operator OL2psl{g, CmL2barpsNC{teta}, IntEps};
	return 6 * ( OL2psc + log(xi) * OL2psl );
      };
    const TabulateObject<Operator> TabOL2s{fOL2s, nxi, ximin, ximax, intdeg, {}, lambda};
    const auto fOL2g = [=,&g] (double const& xi) -> Operator
      {
	const double teta = 1 / ( 1 + 4 / xi );
	const Operator OL2gc{g, CmL2gNC{teta}, IntEps};
	const Operator OL2gl{g, CmL2bargNC{teta}, IntEps};
	return OL2gc + log(xi) * OL2gl;
      };
    const TabulateObject<Operator> TabOL2g{fOL2g, nxi, ximin, ximax, intdeg, {}, lambda};

    // Vector of distributions to skip
    const vector<int> skip = {2,4,6,8,10,12};

    // Define object of the structure containing the DglapObjects
    const auto FLObj = [=] (double const& Q, vector<double> const& Ch) -> StructureFunctionObjects
      {
	// Fill in structure function object
	StructureFunctionObjects FObj;
	FObj.skip = skip;

	// Start with the heavy quark structure function components.
	const double Q2  = Q * Q;
	Operator lNNLOns = OL2nsp;
	for (int k = actnf + 1; k <= 6; k++)
	  {
	    // Compute value of xi.
	    const double M2 = Masses[k-1] * Masses[k-1];
	    const double xi = Q2 / M2;

	    // Convolution Basis
	    FObj.ConvBasis.insert({k,DISNCBasis{k, Ch[k-1]}});

	    // Set LO non-singlet and singlet coefficient
	    // functions to zero (gluon is already zero).
	    FObj.C0.insert({k,Set<Operator>{FObj.ConvBasis.at(k),Not}});

	    // Now insert NLO
	    map<int,Operator> NLO;
	    NLO.insert({CNS, Zero});
	    NLO.insert({CS,  Zero});
	    NLO.insert({CG,  TabOL1g.Evaluate(xi)});
	    FObj.C1.insert({k,Set<Operator>{FObj.ConvBasis.at(k),NLO}});

	    // Now insert NNLO
	    map<int,Operator> NNLO;
	    NNLO.insert({CNS, Zero});
	    NNLO.insert({CS,  TabOL2s.Evaluate(xi)});
	    NNLO.insert({CG,  TabOL2g.Evaluate(xi)});
	    FObj.C2.insert({k,Set<Operator>{FObj.ConvBasis.at(k),NNLO}});

	    // Include the massive bit to the non-singlet coefficient
	    // function needed for the light structure functions.
	    for (int i = 0; i <= actnf; i++)
	      lNNLOns += TabOL2ns.Evaluate(xi);
	  }

	// Now fill in the light components. To be updated in the
	// loop over the heavy component.
	map<int,Operator> CLNNLO;
	CLNNLO.insert({CNS, lNNLOns});
	CLNNLO.insert({CS,  OL2t});
	CLNNLO.insert({CG,  OL2g});
	for (int k = 1; k <= actnf; k++)
	  {
	    // Convolution Basis
	    FObj.ConvBasis.insert({k, DISNCBasis{k, Ch[k-1]}});
	    FObj.C0.insert({k,Set<Operator>{FObj.ConvBasis.at(k), CLLO}});
	    FObj.C1.insert({k,Set<Operator>{FObj.ConvBasis.at(k), CLNLO}});
	    FObj.C2.insert({k,Set<Operator>{FObj.ConvBasis.at(k), CLNNLO}});
	  }

	// Total structure function set to zero for now. Need to
	// construct it a posteriori as a sum of light and heavy
	// components.
	FObj.ConvBasis.insert({0, DISNCBasis{Ch}});
	FObj.C0.insert({0,Set<Operator>{FObj.ConvBasis.at(0), Not}});
	FObj.C1.insert({0,Set<Operator>{FObj.ConvBasis.at(0), Not}});
	FObj.C2.insert({0,Set<Operator>{FObj.ConvBasis.at(0), Not}});
	return FObj;
      };
    t.stop();

    return FLObj;
  }

  //_____________________________________________________________________________
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeF2NCObjectsMassiveZero(Grid           const& g,
													    vector<double> const& Masses,
													    double         const& IntEps,
													    int            const& nxi,
													    double         const& ximin,
													    double         const& ximax,
													    int            const& intdeg,
													    double         const& lambda)
  {
    Timer t;
    t.start();

    // Determine number of active flavours
    int actnf = 0;
    for (auto const& m : Masses)
      if (m < eps8)
	actnf++;

    cout << "Initializing StructureFunctionObjects for F2 NC Massive Zero with " << actnf << " active flavours... \n";

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // Zero Mass coefficient functions
    // LO
    map<int,Operator> C2LO;
    C2LO.insert({CNS, Id});
    C2LO.insert({CS,  Id});
    C2LO.insert({CG,  Zero});

    // NLO
    map<int,Operator> C2NLO;
    const Operator O21ns{g, C21ns{}, IntEps};
    const Operator O21g {g, C21g{},  IntEps};
    C2NLO.insert({CNS, O21ns});
    C2NLO.insert({CS,  O21ns});
    C2NLO.insert({CG,  O21g});

    // NNLO
    const Operator O22ps {g, C22ps{}, IntEps};
    const Operator O22g  {g, C22g{},  IntEps};
    const Operator O22nsp{g, C22nsp{actnf}, IntEps};
    const Operator O22t = O22nsp + 6 * O22ps;

    // Massive zero coefficient functions
    // Null set of operator needed for the LO coefficient functions.
    map<int,Operator> Not;
    Not.insert({CNS, Zero});
    Not.insert({CS,  Zero});
    Not.insert({CG,  Zero});

    // Initialize massive zero coefficient functions
    const Operator Om021gc  {g, Cm021gNC_c{},   IntEps};
    const Operator Om021gl  {g, Cm021gNC_l{},   IntEps};
    const Operator Om022nsc {g, Cm022nsNC_c{},  IntEps};
    const Operator Om022nsl {g, Cm022nsNC_l{},  IntEps};
    const Operator Om022nsl2{g, Cm022nsNC_l2{}, IntEps};
    const Operator Om022psc {g, Cm022psNC_c{},  IntEps};
    const Operator Om022psl {g, Cm022psNC_l{},  IntEps};
    const Operator Om022psl2{g, Cm022psNC_l2{}, IntEps};
    const Operator Om022psf {g, Cm022psNC_f{},  IntEps};
    const Operator Om022pslf{g, Cm022psNC_lf{}, IntEps};
    const Operator Om022gc  {g, Cm022gNC_c{},   IntEps};
    const Operator Om022gl  {g, Cm022gNC_l{},   IntEps};
    const Operator Om022gl2 {g, Cm022gNC_l2{},  IntEps};
    const Operator Om022gf  {g, Cm022gNC_f{},   IntEps};
    const Operator Om022glf {g, Cm022gNC_lf{},  IntEps};

    // NLO
    const TabulateObject<Operator> TabO21g{[=,&g] (double const& xi) -> Operator
	{
	  const double lxi = log(xi);
	  return Om021gc + lxi * Om021gl;
	}, nxi, ximin, ximax, intdeg, {}, lambda};

    // NNLO
    const TabulateObject<Operator> TabO22ns{[=,&g] (double const& xi) -> Operator
	{
	  const double lxi  = log(xi);
	  const double lxi2 = lxi * lxi;
	  return Om022nsc + lxi * Om022nsl + lxi2 * Om022nsl2;
	}, nxi, ximin, ximax, intdeg, {}, lambda};
    const auto fO22s = [=,&g] (double const& xi) -> Operator
      {
	const double lxi  = log(xi);
	const double lxi2 = lxi * lxi;
	const double lxiF = - lxi;
	return 6 * ( Om022psc + lxi * Om022psl + lxi2 * Om022psl2 + lxiF * Om022psf + lxi * lxiF * Om022pslf );
      };
    const TabulateObject<Operator> TabO22s{fO22s, nxi, ximin, ximax, intdeg, {}, lambda};
    const auto fO22g = [=,&g] (double const& xi) -> Operator
      {
	const double lxi  = log(xi);
	const double lxi2 = lxi * lxi;
	const double lxiF = - lxi;
	return Om022gc  + lxi * Om022gl  + lxi2 * Om022gl2 + lxiF * Om022gf + lxi * lxiF * Om022glf;
      };
    const TabulateObject<Operator> TabO22g{fO22g, nxi, ximin, ximax, intdeg, {}, lambda};

    // Vector of distributions to skip
    const vector<int> skip = {2,4,6,8,10,12};

    // Define object of the structure containing the DglapObjects
    const auto F2Obj = [=] (double const& Q, vector<double> const& Ch) -> StructureFunctionObjects
      {
	// Fill in structure function object
	StructureFunctionObjects FObj;
	FObj.skip = skip;

	// Start with the heavy quark structure function components.
	const double Q2  = Q * Q;
	Operator lNNLOns = O22nsp;
	for (int k = actnf + 1; k <= 6; k++)
	  {
	    // Compute value of xi.
	    const double M2 = Masses[k-1] * Masses[k-1];
	    const double xi = Q2 / M2;

	    // Convolution Basis
	    FObj.ConvBasis.insert({k,DISNCBasis{k, Ch[k-1]}});

	    // Set LO non-singlet and singlet coefficient
	    // functions to zero (gluon is already zero).
	    FObj.C0.insert({k,Set<Operator>{FObj.ConvBasis.at(k),Not}});

	    // Now insert NLO
	    map<int,Operator> NLO;
	    NLO.insert({CNS, Zero});
	    NLO.insert({CS,  Zero});
	    NLO.insert({CG,  TabO21g.Evaluate(xi)});
	    FObj.C1.insert({k,Set<Operator>{FObj.ConvBasis.at(k),NLO}});

	    // Now insert NNLO
	    map<int,Operator> NNLO;
	    NNLO.insert({CNS, Zero});
	    NNLO.insert({CS,  TabO22s.Evaluate(xi)});
	    NNLO.insert({CG,  TabO22g.Evaluate(xi)});
	    FObj.C2.insert({k,Set<Operator>{FObj.ConvBasis.at(k),NNLO}});

	    // Include the massive bit to the non-singlet coefficient
	    // function needed for the light structure functions.
	    for (int i = 0; i <= actnf; i++)
	      lNNLOns += TabO22ns.Evaluate(xi);
	  }

	// Now fill in the light components. To be updated in the
	// loop over the heavy component.
	map<int,Operator> C2NNLO;
	C2NNLO.insert({CNS, lNNLOns});
	C2NNLO.insert({CS,  O22t});
	C2NNLO.insert({CG,  O22g});
	for (int k = 1; k <= actnf; k++)
	  {
	    // Convolution Basis
	    FObj.ConvBasis.insert({k, DISNCBasis{k, Ch[k-1]}});
	    FObj.C0.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C2LO}});
	    FObj.C1.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C2NLO}});
	    FObj.C2.insert({k,Set<Operator>{FObj.ConvBasis.at(k), C2NNLO}});
	  }

	// Total structure function set to zero for now. Need to
	// construct it a posteriori as a sum of light and heavy
	// components.
	FObj.ConvBasis.insert({0, DISNCBasis{Ch}});
	FObj.C0.insert({0,Set<Operator>{FObj.ConvBasis.at(0), Not}});
	FObj.C1.insert({0,Set<Operator>{FObj.ConvBasis.at(0), Not}});
	FObj.C2.insert({0,Set<Operator>{FObj.ConvBasis.at(0), Not}});
	return FObj;
      };
    t.stop();

    return F2Obj;
  }

  //_____________________________________________________________________________
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeFLNCObjectsMassiveZero(Grid           const& g,
													    vector<double> const& Masses,
													    double         const& IntEps,
													    int            const& nxi,
													    double         const& ximin,
													    double         const& ximax,
													    int            const& intdeg,
													    double         const& lambda)
  {
    Timer t;
    t.start();

    // Determine number of active flavours
    int actnf = 0;
    for (auto const& m : Masses)
      if (m < eps8)
	actnf++;

    cout << "Initializing StructureFunctionObjects for FL NC Massive Zero with " << actnf << " active flavours... \n";

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps};

    // Zero Mass coefficient functions
    // LO
    map<int,Operator> CLLO;
    CLLO.insert({CNS, Zero});
    CLLO.insert({CS,  Zero});
    CLLO.insert({CG,  Zero});

    // NLO
    map<int,Operator> CLNLO;
    const Operator OL1ns{g, CL1ns{}, IntEps};
    const Operator OL1g {g, CL1g{},  IntEps};
    CLNLO.insert({CNS, OL1ns});
    CLNLO.insert({CS,  OL1ns});
    CLNLO.insert({CG,  OL1g});

    // NNLO
    const Operator OL2ps {g, CL2ps{}, IntEps};
    const Operator OL2g  {g, CL2g{},  IntEps};
    const Operator OL2nsp{g, CL2nsp{actnf}, IntEps};
    const Operator OL2t = OL2nsp + 6 * OL2ps;

    // Massive zero coefficient functions
    // Null set of operator needed for the LO coefficient functions.
    map<int,Operator> Not;
    Not.insert({CNS, Zero});
    Not.insert({CS,  Zero});
    Not.insert({CG,  Zero});

    // Initialize massive zero coefficient functions
    const Operator Om0L1g  {g, Cm0L1gNC_c{},  IntEps};
    const Operator Om0L2nsc{g, Cm0L2nsNC_c{}, IntEps};
    const Operator Om0L2nsl{g, Cm0L2nsNC_l{}, IntEps};
    const Operator Om0L2psc{g, Cm0L2psNC_c{}, IntEps};
    const Operator Om0L2psl{g, Cm0L2psNC_l{}, IntEps};
    const Operator Om0L2psf{g, Cm0L2psNC_f{}, IntEps};
    const Operator Om0L2gc {g, Cm0L2gNC_c{},  IntEps};
    const Operator Om0L2gl {g, Cm0L2gNC_l{},  IntEps};
    const Operator Om0L2gf {g, Cm0L2gNC_f{},  IntEps};

    // NNLO
    const TabulateObject<Operator> TabOL2ns{[=,&g] (double const& xi) -> Operator
	{
	  const double lxi = log(xi);
	  return Om0L2nsc + lxi * Om0L2nsl;
	}, nxi, ximin, ximax, intdeg, {}, lambda};
    const auto fOL2s = [=,&g] (double const& xi) -> Operator
      {
	const double lxi  = log(xi);
	const double lxiF = - lxi;
	return 6 * ( Om0L2psc + lxi * Om0L2psl + lxiF * Om0L2psf );
      };
    const TabulateObject<Operator> TabOL2s{fOL2s, nxi, ximin, ximax, intdeg, {}, lambda};
    const auto fOL2g = [=,&g] (double const& xi) -> Operator
      {
	const double lxi  = log(xi);
	const double lxiF = - lxi;
	return Om0L2gc  + lxi * Om0L2gl + lxiF * Om0L2gf;
      };
    const TabulateObject<Operator> TabOL2g{fOL2g, nxi, ximin, ximax, intdeg, {}, lambda};

    // Vector of distributions to skip
    const vector<int> skip = {2,4,6,8,10,12};

    // Define object of the structure containing the DglapObjects
    const auto FLObj = [=] (double const& Q, vector<double> const& Ch) -> StructureFunctionObjects
      {
	// Fill in structure function object
	StructureFunctionObjects FObj;
	FObj.skip = skip;

	// Start with the heavy quark structure function components.
	const double Q2  = Q * Q;
	Operator lNNLOns = OL2nsp;
	for (int k = actnf + 1; k <= 6; k++)
	  {
	    // Compute value of xi.
	    const double M2 = Masses[k-1] * Masses[k-1];
	    const double xi = Q2 / M2;

	    // Convolution Basis
	    FObj.ConvBasis.insert({k,DISNCBasis{k, Ch[k-1]}});

	    // Set LO non-singlet and singlet coefficient
	    // functions to zero (gluon is already zero).
	    FObj.C0.insert({k,Set<Operator>{FObj.ConvBasis.at(k),Not}});

	    // Now insert NLO
	    map<int,Operator> NLO;
	    NLO.insert({CNS, Zero});
	    NLO.insert({CS,  Zero});
	    NLO.insert({CG,  Om0L1g});
	    FObj.C1.insert({k,Set<Operator>{FObj.ConvBasis.at(k),NLO}});

	    // Now insert NNLO
	    map<int,Operator> NNLO;
	    NNLO.insert({CNS, Zero});
	    NNLO.insert({CS,  TabOL2s.Evaluate(xi)});
	    NNLO.insert({CG,  TabOL2g.Evaluate(xi)});
	    FObj.C2.insert({k,Set<Operator>{FObj.ConvBasis.at(k),NNLO}});

	    // Include the massive bit to the non-singlet coefficient
	    // function needed for the light structure functions.
	    for (int i = 0; i <= actnf; i++)
	      lNNLOns += TabOL2ns.Evaluate(xi);
	  }

	// Now fill in the light components. To be updated in the
	// loop over the heavy component.
	map<int,Operator> CLNNLO;
	CLNNLO.insert({CNS, lNNLOns});
	CLNNLO.insert({CS,  OL2t});
	CLNNLO.insert({CG,  OL2g});
	for (int k = 1; k <= actnf; k++)
	  {
	    // Convolution Basis
	    FObj.ConvBasis.insert({k, DISNCBasis{k, Ch[k-1]}});
	    FObj.C0.insert({k,Set<Operator>{FObj.ConvBasis.at(k), CLLO}});
	    FObj.C1.insert({k,Set<Operator>{FObj.ConvBasis.at(k), CLNLO}});
	    FObj.C2.insert({k,Set<Operator>{FObj.ConvBasis.at(k), CLNNLO}});
	  }

	// Total structure function set to zero for now. Need to
	// construct it a posteriori as a sum of light and heavy
	// components.
	FObj.ConvBasis.insert({0, DISNCBasis{Ch}});
	FObj.C0.insert({0,Set<Operator>{FObj.ConvBasis.at(0), Not}});
	FObj.C1.insert({0,Set<Operator>{FObj.ConvBasis.at(0), Not}});
	FObj.C2.insert({0,Set<Operator>{FObj.ConvBasis.at(0), Not}});
	return FObj;
      };
    t.stop();

    return FLObj;
  }

  //_____________________________________________________________________________
  map<int,Observable> BuildStructureFunctions(function<StructureFunctionObjects(double const&, vector<double> const&)> const& FObj,
					      function<map<int,double>(double const&, double const&)>                  const& InDistFunc,
					      int                                                                      const& PerturbativeOrder,
					      function<double(double const&)>                                          const& Alphas,
					      function<vector<double>(double const&)>                                  const& Couplings)
  {
    // Call FObj at energy 1 to use it for those quantities that do
    // not depend on Q.
    const StructureFunctionObjects FObj1 = FObj(1, Couplings(1));

    // Get grid.
    Grid const& g = FObj1.C0.at(1).at(0).GetGrid();

    // Get skip vector.
    const vector<int> skip = FObj1.skip;

    // Cycle over the key of the convolution basis map.
    map<int,Observable> F;
    for (auto it = FObj1.ConvBasis.begin(); it != FObj1.ConvBasis.end(); ++it)
      {
	// Structure function index.
	const int k = it->first;
	// Define coefficient function functions.
	const auto Cf = [=] (double const& Q) -> Set<Operator>
	  {
	    const double cp = Alphas(Q) / FourPi;
	    const StructureFunctionObjects FObjQ = FObj(Q,Couplings(Q));
	    Set<Operator> CoefFuncs = FObjQ.C0.at(k);
	    if (PerturbativeOrder > 0)
	      CoefFuncs += cp * FObjQ.C1.at(k);
	    if (PerturbativeOrder > 1)
	      CoefFuncs += cp * cp * FObjQ.C2.at(k);
	    return CoefFuncs;
	  };

	// Define distribution function functions.
	const auto DistF = [=,&g] (double const& Q) -> Set<Distribution>{ return Set<Distribution>{FObj(Q,Couplings(Q)).ConvBasis.at(k), DistributionMap(g, InDistFunc, Q, skip)}; };

	// Initialize "Observable".
	F.insert({k,Observable{Cf, DistF}});
      }
    return F;
  }

  //_____________________________________________________________________________
  map<int,Observable> BuildStructureFunctions(function<StructureFunctionObjects(double const&, vector<double> const&)> const& FObj,
					      function<double(int const&, double const&, double const&)>               const& InDistFunc,
					      int                                                                      const& PerturbativeOrder,
					      function<double(double const&)>                                          const& Alphas,
					      function<vector<double>(double const&)>                                  const& Couplings)
  {
    const auto InDistFuncMap = [=] (double const& x, double const& Q) -> map<int,double>
      {
	map<int,double> DistMap;
	for (int i = 0; i <= 12; i++)
	  DistMap.insert({i,InDistFunc(i, x, Q)});
	return DistMap;
      };
    return BuildStructureFunctions(FObj, InDistFuncMap, PerturbativeOrder, Alphas, Couplings);
  }
}
