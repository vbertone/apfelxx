//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/structurefunctionbuilder.h"
#include "apfel/grid.h"
#include "apfel/distributionfunction.h"
#include "apfel/operator.h"
#include "apfel/set.h"
#include "apfel/timer.h"
#include "apfel/tools.h"
#include "apfel/disbasis.h"
#include "apfel/zeromasscoefficientfunctions.h"
#include "apfel/rotations.h"

#include <map>

using namespace std;

namespace apfel {

  //_____________________________________________________________________________
  unordered_map<int,Observable> F2NCBuildZM(Grid                                                       const& g,
					    function<double(int const&, double const&, double const&)> const& InDistFunc,
					    vector<double>                                             const& Thresholds,
					    int                                                        const& PerturbativeOrder,
					    function<double(double const&)>                            const& Alphas,
					    function<vector<double>(double const&)>                    const& Charges,
					    bool                                                       const& RotateInput,
					    double                                                     const& IntEps)
  {
    cout << "Initializing F2NCBuildZM... ";
    Timer t;
    t.start();

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the thresholds
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
	nfi++;

    // Rotate input distributions into the QCD evolution basis if
    // required.
    function<double(int const&, double const&, double const&)> QCDEvPDFsFunc;
    if (RotateInput)
      QCDEvPDFsFunc = [=] (int const& i, double const& x, double const& Q) -> double
	{ return QCDEvToPhys(i, x, Q, InDistFunc); };
    else
      QCDEvPDFsFunc = InDistFunc;

    // Allocate distributions.
    function<unordered_map<int,Distribution>(double const&)> fF2Map = [&g,QCDEvPDFsFunc] (double const& Q) -> unordered_map<int,Distribution>
      {
	unordered_map<int,Distribution> F2Map;
	F2Map.insert({0, DistributionFunction{g, QCDEvPDFsFunc, 0, Q}});
	for (int k = 1; k <= 6; k++)
	  F2Map.insert({k, DistributionFunction{g, QCDEvPDFsFunc, 2 * k - 1, Q}});

	// Change sign to T3 to exchange "up" with "down".
	F2Map.at(2) *= -1;
	return F2Map;
      };

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // Coefficient functions for F2.
    // LO
    unordered_map<int,Operator> C2LO;
    C2LO.insert({DISNCBasis::CNS, Id});
    C2LO.insert({DISNCBasis::CT,  Id});
    C2LO.insert({DISNCBasis::CG,  Zero});

    // NLO
    unordered_map<int,Operator> C2NLO;
    const Operator O21ns{g, C21ns{}, IntEps};
    const Operator O21g {g, C21g{},  IntEps};
    C2NLO.insert({DISNCBasis::CNS, O21ns});
    C2NLO.insert({DISNCBasis::CT,  O21ns});
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
	C2NNLOnf.insert({DISNCBasis::CT,  O22t});
	C2NNLOnf.insert({DISNCBasis::CG,  O22g});
	C2NNLO.insert({nf,C2NNLOnf});
      }

    // Compute observables for F2.
    unordered_map<int,Observable> F2;
    for (int k = 1; k <= 6; k++)
      {
	// Convolution basis.
	const DISNCBasis basis{k};

	// Define sets of operators.
	Set<Operator> LO {basis, C2LO};
	Set<Operator> NLO{basis, C2NLO};
	unordered_map<int,Set<Operator>> NNLO;
	for (int nf = nfi; nf <= nff; nf++)
	  NNLO.insert({nf,Set<Operator>{basis, C2NNLO.at(nf)}});

	// Define coefficient function functions.
	function<Set<Operator>(double const&)> C2f;
	if (PerturbativeOrder == 0)
	  {
	    C2f = [=] (double const& Q) -> Set<Operator>
	      {
		return Charges(Q)[k-1] * LO;
	      };
	  }
	else if (PerturbativeOrder == 1)
	  {
	    C2f = [=] (double const& Q) -> Set<Operator>
	      {
		const auto cp = Alphas(Q) / FourPi;
		return Charges(Q)[k-1] * ( LO + cp * NLO );
	      };
	  }
	else if (PerturbativeOrder == 2)
	  {
	    C2f = [=] (double const& Q) -> Set<Operator>
	      {
		const auto cp = Alphas(Q) / FourPi;
		const auto nf = NF(Q, Thresholds);
		return Charges(Q)[k-1] * ( LO + cp * ( NLO +  + cp * NNLO.at(nf) ) );
	      };
	  }
	// Define distribution function functions.
	function<Set<Distribution>(double const&)> DistF2 = [=] (double const& Q) -> Set<Distribution>
	  {
	    return Set<Distribution>{basis, fF2Map(Q)};
	  };
	// Initialize "Observable".
	F2.insert({k,Observable{C2f, DistF2}});
      }

    // Total structure function.
    function<Set<Operator>(double const&)> C2f;
    if (PerturbativeOrder == 0)
      {
	C2f = [=] (double const& Q) -> Set<Operator>
	  {
	    const DISNCBasis basis{Charges(Q)};
	    Set<Operator> LO{basis, C2LO};
	    return LO;
	  };
      }
    else if (PerturbativeOrder == 1)
      {
	C2f = [=] (double const& Q) -> Set<Operator>
	  {
	    const auto cp = Alphas(Q) / FourPi;
	    const DISNCBasis basis{Charges(Q)};
	    Set<Operator> LO {basis, C2LO};
	    Set<Operator> NLO{basis, C2NLO};
	    return LO + cp * NLO;
	  };
      }
    else if (PerturbativeOrder == 2)
      {
	C2f = [=] (double const& Q) -> Set<Operator>
	  {
	    const auto cp = Alphas(Q) / FourPi;
	    const auto nf = NF(Q, Thresholds);
	    const DISNCBasis basis{Charges(Q)};
	    Set<Operator> LO  {basis, C2LO};
	    Set<Operator> NLO {basis, C2NLO};
	    Set<Operator> NNLO{basis, C2NNLO.at(nf)};
	    return LO + cp * ( NLO +  + cp * NNLO );
	  };
      }
    // Define distribution function functions.
    function<Set<Distribution>(double const&)> DistF2 = [=] (double const& Q) -> Set<Distribution>
      {
	return Set<Distribution>{DISNCBasis{Charges(Q)}, fF2Map(Q)};
      };
    F2.insert({0,Observable{C2f, DistF2}});

    t.stop();

    return F2;
  }

  //_____________________________________________________________________________
  unordered_map<int,Observable> FLNCBuildZM(Grid                                                       const& g,
					    function<double(int const&, double const&, double const&)> const& InDistFunc,
					    vector<double>                                             const& Thresholds,
					    int                                                        const& PerturbativeOrder,
					    function<double(double const&)>                            const& Alphas,
					    function<vector<double>(double const&)>                    const& Charges,
					    bool                                                       const& RotateInput,
					    double                                                     const& IntEps)
  {
    cout << "Initializing FLNCBuildZM... ";
    Timer t;
    t.start();

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the thresholds
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
	nfi++;

    // Rotate input distributions into the QCD evolution basis if
    // required.
    function<double(int const&, double const&, double const&)> QCDEvPDFsFunc;
    if (RotateInput)
      QCDEvPDFsFunc = [=] (int const& i, double const& x, double const& Q) -> double
	{ return QCDEvToPhys(i, x, Q, InDistFunc); };
    else
      QCDEvPDFsFunc = InDistFunc;

    // Allocate distributions.
    function<unordered_map<int,Distribution>(double const&)> fFLMap = [&g,QCDEvPDFsFunc] (double const& Q) -> unordered_map<int,Distribution>
      {
	unordered_map<int,Distribution> FLMap;
	FLMap.insert({0, DistributionFunction{g, QCDEvPDFsFunc, 0, Q}});
	for (int k = 1; k <= 6; k++)
	  FLMap.insert({k, DistributionFunction{g, QCDEvPDFsFunc, 2 * k - 1, Q}});

	// Change sign to T3 to exchange "up" with "down".
	FLMap.at(2) *= -1;
	return FLMap;
      };

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps};

    // Coefficient functions for FL.
    // LO
    unordered_map<int,Operator> CLLO;
    CLLO.insert({DISNCBasis::CNS, Zero});
    CLLO.insert({DISNCBasis::CT,  Zero});
    CLLO.insert({DISNCBasis::CG,  Zero});

    // NLO
    unordered_map<int,Operator> CLNLO;
    const Operator OL1ns{g, CL1ns{}, IntEps};
    const Operator OL1g {g, CL1g{},  IntEps};
    CLNLO.insert({DISNCBasis::CNS, OL1ns});
    CLNLO.insert({DISNCBasis::CT,  OL1ns});
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
	CLNNLOnf.insert({DISNCBasis::CT,  OL2t});
	CLNNLOnf.insert({DISNCBasis::CG,  OL2g});
	CLNNLO.insert({nf,CLNNLOnf});
      }

    // Compute observables for FL.
    unordered_map<int,Observable> FL;
    for (int k = 1; k <= 6; k++)
      {
	// Convolution basis.
	const DISNCBasis basis{k};

	// Define sets of operators.
	Set<Operator> LO {basis, CLLO};
	Set<Operator> NLO{basis, CLNLO};
	unordered_map<int,Set<Operator>> NNLO;
	for (int nf = nfi; nf <= nff; nf++)
	  NNLO.insert({nf,Set<Operator>{basis, CLNNLO.at(nf)}});

	// Define coefficient function functions.
	function<Set<Operator>(double const&)> CLf;
	if (PerturbativeOrder == 0)
	  {
	    CLf = [=] (double const& Q) -> Set<Operator>
	      {
		return Charges(Q)[k-1] * LO;
	      };
	  }
	else if (PerturbativeOrder == 1)
	  {
	    CLf = [=] (double const& Q) -> Set<Operator>
	      {
		const auto cp = Alphas(Q) / FourPi;
		return Charges(Q)[k-1] * ( LO + cp * NLO );
	      };
	  }
	else if (PerturbativeOrder == 2)
	  {
	    CLf = [=] (double const& Q) -> Set<Operator>
	      {
		const auto cp = Alphas(Q) / FourPi;
		const auto nf = NF(Q, Thresholds);
		return Charges(Q)[k-1] * ( LO + cp * ( NLO +  + cp * NNLO.at(nf) ) );
	      };
	  }
	// Define distribution function functions.
	function<Set<Distribution>(double const&)> DistFL = [=] (double const& Q) -> Set<Distribution>
	  {
	    return Set<Distribution>{basis, fFLMap(Q)};
	  };
	// Initialize "Observable".
	FL.insert({k,Observable{CLf, DistFL}});
      }

    // Total structure function.
    function<Set<Operator>(double const&)> CLf;
    if (PerturbativeOrder == 0)
      {
	CLf = [=] (double const& Q) -> Set<Operator>
	  {
	    const DISNCBasis basis{Charges(Q)};
	    Set<Operator> LO{basis, CLLO};
	    return LO;
	  };
      }
    else if (PerturbativeOrder == 1)
      {
	CLf = [=] (double const& Q) -> Set<Operator>
	  {
	    const auto cp = Alphas(Q) / FourPi;
	    const DISNCBasis basis{Charges(Q)};
	    Set<Operator> LO {basis, CLLO};
	    Set<Operator> NLO{basis, CLNLO};
	    return LO + cp * NLO;
	  };
      }
    else if (PerturbativeOrder == 2)
      {
	CLf = [=] (double const& Q) -> Set<Operator>
	  {
	    const auto cp = Alphas(Q) / FourPi;
	    const auto nf = NF(Q, Thresholds);
	    const DISNCBasis basis{Charges(Q)};
	    Set<Operator> LO  {basis, CLLO};
	    Set<Operator> NLO {basis, CLNLO};
	    Set<Operator> NNLO{basis, CLNNLO.at(nf)};
	    return LO + cp * ( NLO +  + cp * NNLO );
	  };
      }
    // Define distribution function functions.
    function<Set<Distribution>(double const&)> DistFL = [=] (double const& Q) -> Set<Distribution>
      {
	return Set<Distribution>{DISNCBasis{Charges(Q)}, fFLMap(Q)};
      };
    FL.insert({0,Observable{CLf, DistFL}});

    t.stop();

    return FL;
  }

  //_____________________________________________________________________________
  unordered_map<int,Observable> F3NCBuildZM(Grid                                                       const& g,
					    function<double(int const&, double const&, double const&)> const& InDistFunc,
					    vector<double>                                             const& Thresholds,
					    int                                                        const& PerturbativeOrder,
					    function<double(double const&)>                            const& Alphas,
					    function<vector<double>(double const&)>                    const& Charges,
					    bool                                                       const& RotateInput,
					    double                                                     const& IntEps)
  {
    cout << "Initializing F3NCBuildZM... ";
    Timer t;
    t.start();

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the thresholds
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
	nfi++;

    // Rotate input distributions into the QCD evolution basis if
    // required.
    function<double(int const&, double const&, double const&)> QCDEvPDFsFunc;
    if (RotateInput)
      QCDEvPDFsFunc = [=] (int const& i, double const& x, double const& Q) -> double
	{ return QCDEvToPhys(i, x, Q, InDistFunc); };
    else
      QCDEvPDFsFunc = InDistFunc;

    // Allocate distributions.
    function<unordered_map<int,Distribution>(double const&)> fF3Map = [&g,QCDEvPDFsFunc] (double const& Q) -> unordered_map<int,Distribution>
      {
	unordered_map<int,Distribution> F3Map;
	F3Map.insert({0, DistributionFunction{g, QCDEvPDFsFunc, 0, Q}});
	for (int k = 1; k <= 6; k++)
	  F3Map.insert({k, DistributionFunction{g, QCDEvPDFsFunc, 2 * k, Q}});

	// Change sign to T3 to exchange "up" with "down".
	F3Map.at(2) *= -1;
	return F3Map;
      };

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // Coefficient functions for F3.
    // LO
    unordered_map<int,Operator> C3LO;
    C3LO.insert({DISNCBasis::CNS, Id});
    C3LO.insert({DISNCBasis::CT,  Id});
    C3LO.insert({DISNCBasis::CG,  Zero});

    // NLO
    unordered_map<int,Operator> C3NLO;
    const Operator O31ns{g, C31ns{}, IntEps};
    C3NLO.insert({DISNCBasis::CNS, O31ns});
    C3NLO.insert({DISNCBasis::CT,  O31ns});
    C3NLO.insert({DISNCBasis::CG,  Zero});

    // NNLO
    unordered_map<int,unordered_map<int,Operator>> C3NNLO;
    for (int nf = nfi; nf <= nff; nf++)
      {
	const Operator O32nsm{g, C32nsm{nf}, IntEps};
	const Operator O32t = O32nsm;
	unordered_map<int,Operator> C3NNLOnf;
	C3NNLOnf.insert({DISNCBasis::CNS, O32nsm});
	C3NNLOnf.insert({DISNCBasis::CT,  O32t});
	C3NNLOnf.insert({DISNCBasis::CG,  Zero});
	C3NNLO.insert({nf,C3NNLOnf});
      }

    // Compute observables for F3.
    unordered_map<int,Observable> F3;
    for (int k = 1; k <= 6; k++)
      {
	// Convolution basis.
	const DISNCBasis basis{k};

	// Define sets of operators.
	Set<Operator> LO {basis, C3LO};
	Set<Operator> NLO{basis, C3NLO};
	unordered_map<int,Set<Operator>> NNLO;
	for (int nf = nfi; nf <= nff; nf++)
	  NNLO.insert({nf,Set<Operator>{basis, C3NNLO.at(nf)}});

	// Define coefficient function functions.
	function<Set<Operator>(double const&)> C3f;
	if (PerturbativeOrder == 0)
	  {
	    C3f = [=] (double const& Q) -> Set<Operator>
	      {
		return Charges(Q)[k-1] * LO;
	      };
	  }
	else if (PerturbativeOrder == 1)
	  {
	    C3f = [=] (double const& Q) -> Set<Operator>
	      {
		const auto cp = Alphas(Q) / FourPi;
		return Charges(Q)[k-1] * ( LO + cp * NLO );
	      };
	  }
	else if (PerturbativeOrder == 2)
	  {
	    C3f = [=] (double const& Q) -> Set<Operator>
	      {
		const auto cp = Alphas(Q) / FourPi;
		const auto nf = NF(Q, Thresholds);
		return Charges(Q)[k-1] * ( LO + cp * ( NLO +  + cp * NNLO.at(nf) ) );
	      };
	  }
	// Define distribution function functions.
	function<Set<Distribution>(double const&)> DistF3 = [=] (double const& Q) -> Set<Distribution>
	  {
	    return Set<Distribution>{basis, fF3Map(Q)};
	  };
	// Initialize "Observable".
	F3.insert({k,Observable{C3f, DistF3}});
      }

    // Total structure function.
    function<Set<Operator>(double const&)> C3f;
    if (PerturbativeOrder == 0)
      {
	C3f = [=] (double const& Q) -> Set<Operator>
	  {
	    const DISNCBasis basis{Charges(Q)};
	    Set<Operator> LO{basis, C3LO};
	    return LO;
	  };
      }
    else if (PerturbativeOrder == 1)
      {
	C3f = [=] (double const& Q) -> Set<Operator>
	  {
	    const auto cp = Alphas(Q) / FourPi;
	    const DISNCBasis basis{Charges(Q)};
	    Set<Operator> LO {basis, C3LO};
	    Set<Operator> NLO{basis, C3NLO};
	    return LO + cp * NLO;
	  };
      }
    else if (PerturbativeOrder == 2)
      {
	C3f = [=] (double const& Q) -> Set<Operator>
	  {
	    const auto cp = Alphas(Q) / FourPi;
	    const auto nf = NF(Q, Thresholds);
	    const DISNCBasis basis{Charges(Q)};
	    Set<Operator> LO  {basis, C3LO};
	    Set<Operator> NLO {basis, C3NLO};
	    Set<Operator> NNLO{basis, C3NNLO.at(nf)};
	    return LO + cp * ( NLO +  + cp * NNLO );
	  };
      }
    // Define distribution function functions.
    function<Set<Distribution>(double const&)> DistF3 = [=] (double const& Q) -> Set<Distribution>
      {
	return Set<Distribution>{DISNCBasis{Charges(Q)}, fF3Map(Q)};
      };
    F3.insert({0,Observable{C3f, DistF3}});

    t.stop();

    return F3;
  }

}
