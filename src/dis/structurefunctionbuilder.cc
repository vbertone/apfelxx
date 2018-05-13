//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/structurefunctionbuilder.h"
#include "apfel/timer.h"
#include "apfel/tools.h"
#include "apfel/constants.h"
#include "apfel/zeromasscoefficientfunctions.h"
#include "apfel/massivecoefficientfunctions.h"
#include "apfel/massivezerocoefficientfunctions.h"
#include "apfel/tabulateobject.h"

namespace apfel {
  //_____________________________________________________________________________
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF2NCObjectsZM(Grid                const& g,
													     std::vector<double> const& Thresholds,
													     double              const& IntEps)
  {
    report("Initializing StructureFunctionObjects for F2 NC Zero Mass... ");
    Timer t;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // LO
    std::map<int,Operator> C2LO;
    C2LO.insert({DISNCBasis::CNS, Id});
    C2LO.insert({DISNCBasis::CS,  Id});
    C2LO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int,Operator> C2NLO;
    const Operator O21ns{g, C21ns{}, IntEps};
    const Operator O21g {g, C21g{},  IntEps};
    C2NLO.insert({DISNCBasis::CNS, O21ns});
    C2NLO.insert({DISNCBasis::CS,  O21ns});
    C2NLO.insert({DISNCBasis::CG,  O21g});

    // NNLO
    std::map<int,std::map<int,Operator>> C2NNLO;
    const Operator O22ps{g, C22ps{}, IntEps};
    const Operator O22g {g, C22g{},  IntEps};
    for (int nf = 1; nf <= 6; nf++)
      {
	const Operator O22nsp{g, C22nsp{nf}, IntEps};
	const Operator O22t = O22nsp + 6 * O22ps;
	std::map<int,Operator> C2NNLOnf;
	C2NNLOnf.insert({DISNCBasis::CNS, O22nsp});
	C2NNLOnf.insert({DISNCBasis::CS,  O22t});
	C2NNLOnf.insert({DISNCBasis::CG,  O22g});
	C2NNLO.insert({nf,C2NNLOnf});
      }

    // Vector of distributions to skip
    const std::vector<int> skip = {2,4,6,8,10,12};

    // Define object of the structure containing the DglapObjects
    const auto F2Obj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
      {
	// Determine number of active flavours.
	const int nf = NF(Q, Thresholds);

	// Effective charges. The charges of the components with mh >
	// Q are set to zero.
	std::vector<double> EffCh;
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
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeFLNCObjectsZM(Grid                const& g,
													     std::vector<double> const& Thresholds,
													     double              const& IntEps)
  {
    report("Initializing StructureFunctionObjects for FL NC Zero Mass... ");
    Timer t;

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps};

    // LO
    std::map<int,Operator> CLLO;
    CLLO.insert({DISNCBasis::CNS, Zero});
    CLLO.insert({DISNCBasis::CS,  Zero});
    CLLO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int,Operator> CLNLO;
    const Operator OL1ns{g, CL1ns{}, IntEps};
    const Operator OL1g {g, CL1g{},  IntEps};
    CLNLO.insert({DISNCBasis::CNS, OL1ns});
    CLNLO.insert({DISNCBasis::CS,  OL1ns});
    CLNLO.insert({DISNCBasis::CG,  OL1g});

    // NNLO (for nf from 1 to 6)
    std::map<int,std::map<int,Operator>> CLNNLO;
    const Operator OL2ps{g, CL2ps{}, IntEps};
    const Operator OL2g {g, CL2g{},  IntEps};
    for (int nf = 1; nf <= 6; nf++)
      {
	const Operator OL2nsp{g, CL2nsp{nf}, IntEps};
	const Operator OL2t = OL2nsp + 6 * OL2ps;
	std::map<int,Operator> CLNNLOnf;
	CLNNLOnf.insert({DISNCBasis::CNS, OL2nsp});
	CLNNLOnf.insert({DISNCBasis::CS,  OL2t});
	CLNNLOnf.insert({DISNCBasis::CG,  OL2g});
	CLNNLO.insert({nf,CLNNLOnf});
      }

    // Vector of distributions to skip
    const std::vector<int> skip = {2,4,6,8,10,12};

    // Define object of the structure containing the DglapObjects
    const auto FLObj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
      {
	// Determine number of active flavours.
	const int nf = NF(Q, Thresholds);

	// Effective charges. The charges of the components with mh >
	// Q are set to zero.
	std::vector<double> EffCh;
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
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF3NCObjectsZM(Grid                const& g,
													     std::vector<double> const& Thresholds,
													     double              const& IntEps)
  {
    report("Initializing StructureFunctionObjects for F3 NC Zero Mass... ");
    Timer t;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // LO
    std::map<int,Operator> C3LO;
    C3LO.insert({DISNCBasis::CNS, Id});
    C3LO.insert({DISNCBasis::CS,  Id});
    C3LO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int,Operator> C3NLO;
    const Operator O31ns{g, C31ns{}, IntEps};
    C3NLO.insert({DISNCBasis::CNS, O31ns});
    C3NLO.insert({DISNCBasis::CS,  O31ns});
    C3NLO.insert({DISNCBasis::CG,  Zero});

    // NNLO
    std::map<int,std::map<int,Operator>> C3NNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
	const Operator O32nsm{g, C32nsm{nf}, IntEps};
	const Operator O32t = O32nsm;
	std::map<int,Operator> C3NNLOnf;
	C3NNLOnf.insert({DISNCBasis::CNS, O32nsm});
	C3NNLOnf.insert({DISNCBasis::CS,  O32t});
	C3NNLOnf.insert({DISNCBasis::CG,  Zero});
	C3NNLO.insert({nf,C3NNLOnf});
      }

    // Vector of distributions to skip
    const std::vector<int> skip = {1,3,5,7,9,11};

    // Define object of the structure containing the DglapObjects
    const auto F3Obj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
      {
	// Determine number of active flavours.
	const int nf = NF(Q, Thresholds);

	// Effective charges. The charges of the components with mh >
	// Q are set to zero.
	std::vector<double> EffCh;
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
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF2CCPlusObjectsZM(Grid                const& g,
														 std::vector<double> const& Thresholds,
														 double              const& IntEps)
  {
    report("Initializing StructureFunctionObjects for ( F2(nu) + F2(nubar) ) / 2 Zero Mass... ");
    Timer t;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // LO
    std::map<int,Operator> C2LO;
    C2LO.insert({DISCCBasis::CNS, Id});
    C2LO.insert({DISCCBasis::CS,  Id});
    C2LO.insert({DISCCBasis::CG,  Zero});

    // NLO
    std::map<int,Operator> C2NLO;
    const Operator O21ns{g, C21ns{}, IntEps};
    const Operator O21g {g, C21g{},  IntEps};
    C2NLO.insert({DISCCBasis::CNS, O21ns});
    C2NLO.insert({DISCCBasis::CS,  O21ns});
    C2NLO.insert({DISCCBasis::CG,  O21g});

    // NNLO
    std::map<int,std::map<int,Operator>> C2NNLO;
    const Operator O22ps{g, C22ps{}, IntEps};
    const Operator O22g {g, C22g{},  IntEps};
    for (int nf = 1; nf <= 6; nf++)
      {
	const Operator O22nsp{g, C22nsp{nf}, IntEps};
	const Operator O22t = O22nsp + 6 * O22ps;
	std::map<int,Operator> C2NNLOnf;
	C2NNLOnf.insert({DISCCBasis::CNS, O22nsp});
	C2NNLOnf.insert({DISCCBasis::CS,  O22t});
	C2NNLOnf.insert({DISCCBasis::CG,  O22g});
	C2NNLO.insert({nf,C2NNLOnf});
      }

    // Vector of distributions to skip
    const std::vector<int> skip = {2,4,6,8,10,12};

    // Define object of the structure containing the DglapObjects
    const auto F2Obj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
      {
	// Determine number of active flavours.
	const int nf = NF(Q, Thresholds);

	// Effective charges.
	std::vector<double> EffCh;
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
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF2CCMinusObjectsZM(Grid                const& g,
														  std::vector<double> const& Thresholds,
														  double              const& IntEps)
  {
    report("Initializing StructureFunctionObjects for ( F2(nu) - F2(nubar) ) / 2 Zero Mass... ");
    Timer t;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // LO
    std::map<int,Operator> C2LO;
    C2LO.insert({DISCCBasis::CNS, Id});
    C2LO.insert({DISCCBasis::CS,  Zero});
    C2LO.insert({DISCCBasis::CG,  Zero});

    // NLO
    std::map<int,Operator> C2NLO;
    const Operator O21ns{g, C21ns{}, IntEps};
    C2NLO.insert({DISCCBasis::CNS, O21ns});
    C2NLO.insert({DISCCBasis::CS,  Zero});
    C2NLO.insert({DISCCBasis::CG,  Zero});

    // NNLO
    std::map<int,std::map<int,Operator>> C2NNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
	const Operator O22nsm{g, C22nsm{nf}, IntEps};
	std::map<int,Operator> C2NNLOnf;
	C2NNLOnf.insert({DISCCBasis::CNS, O22nsm});
	C2NNLOnf.insert({DISCCBasis::CS,  Zero});
	C2NNLOnf.insert({DISCCBasis::CG,  Zero});
	C2NNLO.insert({nf,C2NNLOnf});
      }

    // Vector of distributions to skip
    const std::vector<int> skip = {0,1,2,3,5,7,9,11};

    // Define object of the structure containing the DglapObjects
    const auto F2Obj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
      {
	// Determine number of active flavours.
	const int nf = NF(Q, Thresholds);

	// Effective charges.
	std::vector<double> EffCh;
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
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeFLCCPlusObjectsZM(Grid                const& g,
														 std::vector<double> const& Thresholds,
														 double              const& IntEps)
  {
    report("Initializing StructureFunctionObjects for ( FL(nu) + FL(nubar) ) / 2 Zero Mass... ");
    Timer t;

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps};

    // LO
    std::map<int,Operator> CLLO;
    CLLO.insert({DISCCBasis::CNS, Zero});
    CLLO.insert({DISCCBasis::CS,  Zero});
    CLLO.insert({DISCCBasis::CG,  Zero});

    // NLO
    std::map<int,Operator> CLNLO;
    const Operator OL1ns{g, CL1ns{}, IntEps};
    const Operator OL1g {g, CL1g{},  IntEps};
    CLNLO.insert({DISCCBasis::CNS, OL1ns});
    CLNLO.insert({DISCCBasis::CS,  OL1ns});
    CLNLO.insert({DISCCBasis::CG,  OL1g});

    // NNLO
    std::map<int,std::map<int,Operator>> CLNNLO;
    const Operator OL2ps{g, CL2ps{}, IntEps};
    const Operator OL2g {g, CL2g{},  IntEps};
    for (int nf = 1; nf <= 6; nf++)
      {
	const Operator OL2nsp{g, CL2nsp{nf}, IntEps};
	const Operator OL2t = OL2nsp + 6 * OL2ps;
	std::map<int,Operator> CLNNLOnf;
	CLNNLOnf.insert({DISCCBasis::CNS, OL2nsp});
	CLNNLOnf.insert({DISCCBasis::CS,  OL2t});
	CLNNLOnf.insert({DISCCBasis::CG,  OL2g});
	CLNNLO.insert({nf,CLNNLOnf});
      }

    // Vector of distributions to skip
    const std::vector<int> skip = {2,4,6,8,10,12};

    // Define object of the structure containing the DglapObjects
    const auto FLObj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
      {
	// Determine number of active flavours.
	const int nf = NF(Q, Thresholds);

	// Effective charges.
	std::vector<double> EffCh;
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
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeFLCCMinusObjectsZM(Grid                const& g,
														  std::vector<double> const& Thresholds,
														  double              const& IntEps)
  {
    report("Initializing StructureFunctionObjects for ( FL(nu) - FL(nubar) ) / 2 Zero Mass... ");
    Timer t;

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps};

    // LO
    std::map<int,Operator> CLLO;
    CLLO.insert({DISCCBasis::CNS, Zero});
    CLLO.insert({DISCCBasis::CS,  Zero});
    CLLO.insert({DISCCBasis::CG,  Zero});

    // NLO
    std::map<int,Operator> CLNLO;
    const Operator OL1ns{g, CL1ns{}, IntEps};
    CLNLO.insert({DISCCBasis::CNS, OL1ns});
    CLNLO.insert({DISCCBasis::CS,  Zero});
    CLNLO.insert({DISCCBasis::CG,  Zero});

    // NNLO
    std::map<int,std::map<int,Operator>> CLNNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
	const Operator OL2nsm{g, CL2nsm{nf}, IntEps};
	std::map<int,Operator> CLNNLOnf;
	CLNNLOnf.insert({DISCCBasis::CNS, OL2nsm});
	CLNNLOnf.insert({DISCCBasis::CS,  Zero});
	CLNNLOnf.insert({DISCCBasis::CG,  Zero});
	CLNNLO.insert({nf,CLNNLOnf});
      }

    // Vector of distributions to skip
    const std::vector<int> skip = {0,1,2,3,5,7,9,11};

    // Define object of the structure containing the DglapObjects
    const auto FLObj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
      {
	// Determine number of active flavours.
	const int nf = NF(Q, Thresholds);

	// Effective charges.
	std::vector<double> EffCh;
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
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF3CCPlusObjectsZM(Grid                const& g,
														 std::vector<double> const& Thresholds,
														 double              const& IntEps)
  {
    report("Initializing StructureFunctionObjects for ( F3(nu) + F3(nubar) ) / 2 Zero Mass... ");
    Timer t;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // LO
    std::map<int,Operator> C3LO;
    C3LO.insert({DISCCBasis::CNS, Id});
    C3LO.insert({DISCCBasis::CS,  Id});
    C3LO.insert({DISCCBasis::CG,  Zero});

    // NLO
    std::map<int,Operator> C3NLO;
    const Operator O31ns{g, C31ns{}, IntEps};
    C3NLO.insert({DISCCBasis::CNS, O31ns});
    C3NLO.insert({DISCCBasis::CS,  O31ns});
    C3NLO.insert({DISCCBasis::CG,  Zero});

    // NNLO
    std::map<int,std::map<int,Operator>> C3NNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
	const Operator O32nsp{g, C32nsp{nf}, IntEps};
	const Operator O32t = O32nsp;
	std::map<int,Operator> C3NNLOnf;
	C3NNLOnf.insert({DISCCBasis::CNS, O32nsp});
	C3NNLOnf.insert({DISCCBasis::CS,  O32t});
	C3NNLOnf.insert({DISCCBasis::CG,  Zero});
	C3NNLO.insert({nf,C3NNLOnf});
      }

    // Vector of distributions to skip
    const std::vector<int> skip = {0,1,2,4,6,8,10,12};

    // Define object of the structure containing the DglapObjects
    const auto F3Obj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
      {
	// Determine number of active flavours.
	const int nf = NF(Q, Thresholds);

	// Effective charges.
	std::vector<double> EffCh;
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
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF3CCMinusObjectsZM(Grid                const& g,
														  std::vector<double> const& Thresholds,
														  double              const& IntEps)
  {
    report("Initializing StructureFunctionObjects for ( F3(nu) - F3(nubar) ) / 2 Zero Mass... ");
    Timer t;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // LO
    std::map<int,Operator> C3LO;
    C3LO.insert({DISCCBasis::CNS, Id});
    C3LO.insert({DISCCBasis::CS,  Id});
    C3LO.insert({DISCCBasis::CG,  Zero});

    // NLO
    std::map<int,Operator> C3NLO;
    const Operator O31ns{g, C31ns{}, IntEps};
    C3NLO.insert({DISCCBasis::CNS, O31ns});
    C3NLO.insert({DISCCBasis::CS,  O31ns});
    C3NLO.insert({DISCCBasis::CG,  Zero});

    // NNLO
    std::map<int,std::map<int,Operator>> C3NNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
	const Operator O32nsm{g, C32nsm{nf}, IntEps};
	std::map<int,Operator> C3NNLOnf;
	C3NNLOnf.insert({DISCCBasis::CNS, O32nsm});
	C3NNLOnf.insert({DISCCBasis::CS,  O32nsm});
	C3NNLOnf.insert({DISCCBasis::CG,  Zero});
	C3NNLO.insert({nf,C3NNLOnf});
      }

    // Vector of distributions to skip
    const std::vector<int> skip = {0,1,3,5,7,9,11};

    // Define object of the structure containing the DglapObjects
    const auto F3Obj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
      {
	// Determine number of active flavours.
	const int nf = NF(Q, Thresholds);

	// Effective charges.
	std::vector<double> EffCh;
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
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF2NCObjectsMassive(Grid                const& g,
														  std::vector<double> const& Masses,
														  double              const& IntEps,
														  int                 const& nxi,
														  double              const& ximin,
														  double              const& ximax,
														  int                 const& intdeg,
														  double              const& lambda)
  {
    Timer t;

    // Determine number of active flavours
    int actnf = 0;
    for (auto const& m : Masses)
      if (m < eps8)
	actnf++;

    report("Initializing StructureFunctionObjects for F2 NC Massive with " + std::to_string(actnf) + " active flavours... \n");

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // Zero Mass coefficient functions
    // LO
    std::map<int,Operator> C2LO;
    C2LO.insert({DISNCBasis::CNS, Id});
    C2LO.insert({DISNCBasis::CS,  Id});
    C2LO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int,Operator> C2NLO;
    const Operator O21ns{g, C21ns{}, IntEps};
    const Operator O21g {g, C21g{},  IntEps};
    C2NLO.insert({DISNCBasis::CNS, O21ns});
    C2NLO.insert({DISNCBasis::CS,  O21ns});
    C2NLO.insert({DISNCBasis::CG,  O21g});

    // NNLO
    const Operator O22ps {g, C22ps{}, IntEps};
    const Operator O22g  {g, C22g{},  IntEps};
    const Operator O22nsp{g, C22nsp{actnf}, IntEps};
    const Operator O22t = O22nsp + 6 * O22ps;

    // Massive coefficient functions
    // Null set of operator needed for the LO coefficient functions.
    std::map<int,Operator> Not;
    Not.insert({DISNCBasis::CNS, Zero});
    Not.insert({DISNCBasis::CS,  Zero});
    Not.insert({DISNCBasis::CG,  Zero});

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
    const std::vector<int> skip = {2,4,6,8,10,12};

    // Define object of the structure containing the DglapObjects
    const auto F2Obj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
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
	    std::map<int,Operator> NLO;
	    NLO.insert({DISNCBasis::CNS, Zero});
	    NLO.insert({DISNCBasis::CS,  Zero});
	    NLO.insert({DISNCBasis::CG,  TabO21g.Evaluate(xi)});
	    FObj.C1.insert({k,Set<Operator>{FObj.ConvBasis.at(k),NLO}});

	    // Now insert NNLO
	    std::map<int,Operator> NNLO;
	    NNLO.insert({DISNCBasis::CNS, Zero});
	    NNLO.insert({DISNCBasis::CS,  TabO22s.Evaluate(xi)});
	    NNLO.insert({DISNCBasis::CG,  TabO22g.Evaluate(xi)});
	    FObj.C2.insert({k,Set<Operator>{FObj.ConvBasis.at(k),NNLO}});

	    // Include the massive bit to the non-singlet coefficient
	    // function needed for the light structure functions.
	    for (int i = 0; i <= actnf; i++)
	      lNNLOns += TabO22ns.Evaluate(xi);
	  }

	// Now fill in the light components. To be updated in the
	// loop over the heavy component.
	std::map<int,Operator> C2NNLO;
	C2NNLO.insert({DISNCBasis::CNS, lNNLOns});
	C2NNLO.insert({DISNCBasis::CS,  O22t});
	C2NNLO.insert({DISNCBasis::CG,  O22g});
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
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeFLNCObjectsMassive(Grid                const& g,
														  std::vector<double> const& Masses,
														  double              const& IntEps,
														  int                 const& nxi,
														  double              const& ximin,
														  double              const& ximax,
														  int                 const& intdeg,
														  double              const& lambda)
  {
    Timer t;

    // Determine number of active flavours
    int actnf = 0;
    for (auto const& m : Masses)
      if (m < eps8)
	actnf++;

    report("Initializing StructureFunctionObjects for FL NC Massive with " + std::to_string(actnf) + " active flavours... \n");

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps};

    // Zero Mass coefficient functions
    // LO
    std::map<int,Operator> CLLO;
    CLLO.insert({DISNCBasis::CNS, Zero});
    CLLO.insert({DISNCBasis::CS,  Zero});
    CLLO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int,Operator> CLNLO;
    const Operator OL1ns{g, CL1ns{}, IntEps};
    const Operator OL1g {g, CL1g{},  IntEps};
    CLNLO.insert({DISNCBasis::CNS, OL1ns});
    CLNLO.insert({DISNCBasis::CS,  OL1ns});
    CLNLO.insert({DISNCBasis::CG,  OL1g});

    // NNLO
    const Operator OL2ps {g, CL2ps{}, IntEps};
    const Operator OL2g  {g, CL2g{},  IntEps};
    const Operator OL2nsp{g, CL2nsp{actnf}, IntEps};
    const Operator OL2t = OL2nsp + 6 * OL2ps;

    // Massive coefficient functions
    // Null set of operator needed for the LO coefficient functions.
    std::map<int,Operator> Not;
    Not.insert({DISNCBasis::CNS, Zero});
    Not.insert({DISNCBasis::CS,  Zero});
    Not.insert({DISNCBasis::CG,  Zero});

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
    const std::vector<int> skip = {2,4,6,8,10,12};

    // Define object of the structure containing the DglapObjects
    const auto FLObj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
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
	    std::map<int,Operator> NLO;
	    NLO.insert({DISNCBasis::CNS, Zero});
	    NLO.insert({DISNCBasis::CS,  Zero});
	    NLO.insert({DISNCBasis::CG,  TabOL1g.Evaluate(xi)});
	    FObj.C1.insert({k,Set<Operator>{FObj.ConvBasis.at(k),NLO}});

	    // Now insert NNLO
	    std::map<int,Operator> NNLO;
	    NNLO.insert({DISNCBasis::CNS, Zero});
	    NNLO.insert({DISNCBasis::CS,  TabOL2s.Evaluate(xi)});
	    NNLO.insert({DISNCBasis::CG,  TabOL2g.Evaluate(xi)});
	    FObj.C2.insert({k,Set<Operator>{FObj.ConvBasis.at(k),NNLO}});

	    // Include the massive bit to the non-singlet coefficient
	    // function needed for the light structure functions.
	    for (int i = 0; i <= actnf; i++)
	      lNNLOns += TabOL2ns.Evaluate(xi);
	  }

	// Now fill in the light components. To be updated in the
	// loop over the heavy component.
	std::map<int,Operator> CLNNLO;
	CLNNLO.insert({DISNCBasis::CNS, lNNLOns});
	CLNNLO.insert({DISNCBasis::CS,  OL2t});
	CLNNLO.insert({DISNCBasis::CG,  OL2g});
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
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF2NCObjectsMassiveZero(Grid                const& g,
														      std::vector<double> const& Masses,
														      double              const& IntEps,
														      int                 const& nxi,
														      double              const& ximin,
														      double              const& ximax,
														      int                 const& intdeg,
														      double              const& lambda)
  {
    Timer t;

    // Determine number of active flavours
    int actnf = 0;
    for (auto const& m : Masses)
      if (m < eps8)
	actnf++;

    report("Initializing StructureFunctionObjects for F2 NC Massive Zero with " + std::to_string(actnf) + " active flavours... \n");

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // Zero Mass coefficient functions
    // LO
    std::map<int,Operator> C2LO;
    C2LO.insert({DISNCBasis::CNS, Id});
    C2LO.insert({DISNCBasis::CS,  Id});
    C2LO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int,Operator> C2NLO;
    const Operator O21ns{g, C21ns{}, IntEps};
    const Operator O21g {g, C21g{},  IntEps};
    C2NLO.insert({DISNCBasis::CNS, O21ns});
    C2NLO.insert({DISNCBasis::CS,  O21ns});
    C2NLO.insert({DISNCBasis::CG,  O21g});

    // NNLO
    const Operator O22ps {g, C22ps{}, IntEps};
    const Operator O22g  {g, C22g{},  IntEps};
    const Operator O22nsp{g, C22nsp{actnf}, IntEps};
    const Operator O22t = O22nsp + 6 * O22ps;

    // Massive zero coefficient functions
    // Null set of operator needed for the LO coefficient functions.
    std::map<int,Operator> Not;
    Not.insert({DISNCBasis::CNS, Zero});
    Not.insert({DISNCBasis::CS,  Zero});
    Not.insert({DISNCBasis::CG,  Zero});

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
    const TabulateObject<Operator> TabO21g{[=] (double const& xi) -> Operator
	{
	  const double lxi = log(xi);
	  return Om021gc + lxi * Om021gl;
	}, nxi, ximin, ximax, intdeg, {}, lambda};

    // NNLO
    const TabulateObject<Operator> TabO22ns{[=] (double const& xi) -> Operator
	{
	  const double lxi  = log(xi);
	  const double lxi2 = lxi * lxi;
	  return Om022nsc + lxi * Om022nsl + lxi2 * Om022nsl2;
	}, nxi, ximin, ximax, intdeg, {}, lambda};
    const auto fO22s = [=] (double const& xi) -> Operator
      {
	const double lxi  = log(xi);
	const double lxi2 = lxi * lxi;
	const double lxiF = - lxi;
	return 6 * ( Om022psc + lxi * Om022psl + lxi2 * Om022psl2 + lxiF * Om022psf + lxi * lxiF * Om022pslf );
      };
    const TabulateObject<Operator> TabO22s{fO22s, nxi, ximin, ximax, intdeg, {}, lambda};
    const auto fO22g = [=] (double const& xi) -> Operator
      {
	const double lxi  = log(xi);
	const double lxi2 = lxi * lxi;
	const double lxiF = - lxi;
	return Om022gc  + lxi * Om022gl  + lxi2 * Om022gl2 + lxiF * Om022gf + lxi * lxiF * Om022glf;
      };
    const TabulateObject<Operator> TabO22g{fO22g, nxi, ximin, ximax, intdeg, {}, lambda};

    // Vector of distributions to skip
    const std::vector<int> skip = {2,4,6,8,10,12};

    // Define object of the structure containing the DglapObjects
    const auto F2Obj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
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
	    std::map<int,Operator> NLO;
	    NLO.insert({DISNCBasis::CNS, Zero});
	    NLO.insert({DISNCBasis::CS,  Zero});
	    NLO.insert({DISNCBasis::CG,  TabO21g.Evaluate(xi)});
	    FObj.C1.insert({k,Set<Operator>{FObj.ConvBasis.at(k),NLO}});

	    // Now insert NNLO
	    std::map<int,Operator> NNLO;
	    NNLO.insert({DISNCBasis::CNS, Zero});
	    NNLO.insert({DISNCBasis::CS,  TabO22s.Evaluate(xi)});
	    NNLO.insert({DISNCBasis::CG,  TabO22g.Evaluate(xi)});
	    FObj.C2.insert({k,Set<Operator>{FObj.ConvBasis.at(k),NNLO}});

	    // Include the massive bit to the non-singlet coefficient
	    // function needed for the light structure functions.
	    for (int i = 0; i <= actnf; i++)
	      lNNLOns += TabO22ns.Evaluate(xi);
	  }

	// Now fill in the light components. To be updated in the
	// loop over the heavy component.
	std::map<int,Operator> C2NNLO;
	C2NNLO.insert({DISNCBasis::CNS, lNNLOns});
	C2NNLO.insert({DISNCBasis::CS,  O22t});
	C2NNLO.insert({DISNCBasis::CG,  O22g});
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
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeFLNCObjectsMassiveZero(Grid                const& g,
														      std::vector<double> const& Masses,
														      double              const& IntEps,
														      int                 const& nxi,
														      double              const& ximin,
														      double              const& ximax,
														      int                 const& intdeg,
														      double              const& lambda)
  {
    Timer t;

    // Determine number of active flavours
    int actnf = 0;
    for (auto const& m : Masses)
      if (m < eps8)
	actnf++;

    report("Initializing StructureFunctionObjects for FL NC Massive Zero with " + std::to_string(actnf) + " active flavours... \n");

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps};

    // Zero Mass coefficient functions
    // LO
    std::map<int,Operator> CLLO;
    CLLO.insert({DISNCBasis::CNS, Zero});
    CLLO.insert({DISNCBasis::CS,  Zero});
    CLLO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int,Operator> CLNLO;
    const Operator OL1ns{g, CL1ns{}, IntEps};
    const Operator OL1g {g, CL1g{},  IntEps};
    CLNLO.insert({DISNCBasis::CNS, OL1ns});
    CLNLO.insert({DISNCBasis::CS,  OL1ns});
    CLNLO.insert({DISNCBasis::CG,  OL1g});

    // NNLO
    const Operator OL2ps {g, CL2ps{}, IntEps};
    const Operator OL2g  {g, CL2g{},  IntEps};
    const Operator OL2nsp{g, CL2nsp{actnf}, IntEps};
    const Operator OL2t = OL2nsp + 6 * OL2ps;

    // Massive zero coefficient functions
    // Null set of operator needed for the LO coefficient functions.
    std::map<int,Operator> Not;
    Not.insert({DISNCBasis::CNS, Zero});
    Not.insert({DISNCBasis::CS,  Zero});
    Not.insert({DISNCBasis::CG,  Zero});

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
    const TabulateObject<Operator> TabOL2ns{[=] (double const& xi) -> Operator
	{
	  const double lxi = log(xi);
	  return Om0L2nsc + lxi * Om0L2nsl;
	}, nxi, ximin, ximax, intdeg, {}, lambda};
    const auto fOL2s = [=] (double const& xi) -> Operator
      {
	const double lxi  = log(xi);
	const double lxiF = - lxi;
	return 6 * ( Om0L2psc + lxi * Om0L2psl + lxiF * Om0L2psf );
      };
    const TabulateObject<Operator> TabOL2s{fOL2s, nxi, ximin, ximax, intdeg, {}, lambda};
    const auto fOL2g = [=] (double const& xi) -> Operator
      {
	const double lxi  = log(xi);
	const double lxiF = - lxi;
	return Om0L2gc  + lxi * Om0L2gl + lxiF * Om0L2gf;
      };
    const TabulateObject<Operator> TabOL2g{fOL2g, nxi, ximin, ximax, intdeg, {}, lambda};

    // Vector of distributions to skip
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

    // Define object of the structure containing the DglapObjects
    const auto FLObj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
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
	    std::map<int,Operator> NLO;
	    NLO.insert({DISNCBasis::CNS, Zero});
	    NLO.insert({DISNCBasis::CS,  Zero});
	    NLO.insert({DISNCBasis::CG,  Om0L1g});
	    FObj.C1.insert({k,Set<Operator>{FObj.ConvBasis.at(k),NLO}});

	    // Now insert NNLO
	    std::map<int,Operator> NNLO;
	    NNLO.insert({DISNCBasis::CNS, Zero});
	    NNLO.insert({DISNCBasis::CS,  TabOL2s.Evaluate(xi)});
	    NNLO.insert({DISNCBasis::CG,  TabOL2g.Evaluate(xi)});
	    FObj.C2.insert({k,Set<Operator>{FObj.ConvBasis.at(k),NNLO}});

	    // Include the massive bit to the non-singlet coefficient
	    // function needed for the light structure functions.
	    for (int i = 0; i <= actnf; i++)
	      lNNLOns += TabOL2ns.Evaluate(xi);
	  }

	// Now fill in the light components. To be updated in the
	// loop over the heavy component.
	std::map<int,Operator> CLNNLO;
	CLNNLO.insert({DISNCBasis::CNS, lNNLOns});
	CLNNLO.insert({DISNCBasis::CS,  OL2t});
	CLNNLO.insert({DISNCBasis::CG,  OL2g});
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
  std::map<int,Observable> BuildStructureFunctions(std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> const& FObj,
						   std::function<std::map<int,double>(double const&, double const&)>                  const& InDistFunc,
						   int                                                                                const& PerturbativeOrder,
						   std::function<double(double const&)>                                               const& Alphas,
						   std::function<std::vector<double>(double const&)>                                  const& Couplings)
  {
    // Call FObj at energy 1 to use it for those quantities that do
    // not depend on Q.
    const StructureFunctionObjects FObj1 = FObj(1, Couplings(1));

    // Get grid.
    Grid const& g = FObj1.C0.at(1).at(0).GetGrid();

    // Get skip vector.
    const std::vector<int> skip = FObj1.skip;

    // Cycle over the key of the convolution basis map.
    std::map<int,Observable> F;
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
	const auto DistF = [=,&g] (double const& Q) -> Set<Distribution>
	  {
	    return Set<Distribution>{FObj(Q,Couplings(Q)).ConvBasis.at(k), DistributionMap(g, InDistFunc, Q, skip)};
	  };

	// Initialize "Observable".
	F.insert({k,Observable{Cf, DistF}});
      }
    return F;
  }

  //_____________________________________________________________________________
  std::map<int,Observable> BuildStructureFunctions(std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> const& FObj,
						   std::function<double(int const&, double const&, double const&)>                    const& InDistFunc,
						   int                                                                                const& PerturbativeOrder,
						   std::function<double(double const&)>                                               const& Alphas,
						   std::function<std::vector<double>(double const&)>                                  const& Couplings)
  {
    const auto InDistFuncMap = [=] (double const& x, double const& Q) -> std::map<int,double>
      {
	std::map<int,double> DistMap;
	for (int i = 0; i <= 12; i++)
	  DistMap.insert({i,InDistFunc(i, x, Q)});
	return DistMap;
      };
    return BuildStructureFunctions(FObj, InDistFuncMap, PerturbativeOrder, Alphas, Couplings);
  }
}
