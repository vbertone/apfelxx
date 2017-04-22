//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/structurefunctionbuilder.h"
#include <apfel/grid.h>
#include <apfel/distributionfunction.h>
#include "apfel/operator.h"
#include "apfel/set.h"
#include "apfel/timer.h"
#include "apfel/tools.h"
#include "apfel/disbasis.h"
#include <apfel/zeromasscoefficientfunctions.h>

#include <map>

using namespace std;

namespace apfel {

  //_____________________________________________________________________________
  unordered_map<int,Observable> F2BuildZM(Grid                                        const& g,
					  function<double(int const&, double const&)> const& InDistFunc,
					  vector<double>                              const& Thresholds,
					  int                                         const& PerturbativeOrder,
					  function<double(double const&)>             const& Alphas,
					  function<vector<double>(double const&)>     const& Charges,
					  double                                      const& IntEps)
  {
    cout << "Initializing F2BuildZM... ";
    Timer t;
    t.start();

    // Compute initial and final number of active flavours 
    // according to the vector of thresholds (it assumes that
    // the thresholds vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
	nfi++;

    // Allocate distributions
    unordered_map<int,Distribution> F2Map;
    F2Map.insert({0, DistributionFunction{g, InDistFunc, 0 }});
    F2Map.insert({1, DistributionFunction{g, InDistFunc, 1 }});
    F2Map.insert({2, DistributionFunction{g, InDistFunc, 3 }});
    F2Map.insert({3, DistributionFunction{g, InDistFunc, 5 }});
    F2Map.insert({4, DistributionFunction{g, InDistFunc, 7 }});
    F2Map.insert({5, DistributionFunction{g, InDistFunc, 9 }});
    F2Map.insert({6, DistributionFunction{g, InDistFunc, 11}});

    // Change sign to T3 to exchange "up" with "down"
    F2Map.at(2) *= -1;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // Coefficient functions for F2
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
	const Operator O22t = O22nsp + nf * O22ps;
	unordered_map<int,Operator> C2NNLOnf;
	C2NNLOnf.insert({DISNCBasis::CNS, O22nsp});
	C2NNLOnf.insert({DISNCBasis::CT,  O22t});
	C2NNLOnf.insert({DISNCBasis::CG,  O22g});
	C2NNLO.insert({nf,C2NNLOnf});
      }

    // Compute observables for F2
    unordered_map<int,Observable> F2;
    for (int k = 1; k <= 6; k++)
      {
	// Define coefficient function functions
	function<Set<Operator>(double const&)> C2f;
	if (PerturbativeOrder == 0)
	  {
	    C2f = [=] (double const& Q) -> Set<Operator>
	      {
		const auto nf = NF(Q, Thresholds);
		Set<Operator> LO{DISNCBasis{k,nf}, C2LO};
		return Charges(Q)[k-1] * LO;
	      };
	  }
	else if (PerturbativeOrder == 1)
	  {
	    C2f = [=] (double const& Q) -> Set<Operator>
	      {
		const auto cp = Alphas(Q) / FourPi;
		const auto nf = NF(Q, Thresholds);
		Set<Operator> LO {DISNCBasis{k,nf}, C2LO};
		Set<Operator> NLO{DISNCBasis{k,nf}, C2NLO};
		return Charges(Q)[k-1] * ( LO + cp * NLO );
	      };
	  }
	else if (PerturbativeOrder == 2)
	  {
	    C2f = [=] (double const& Q) -> Set<Operator>
	      {
		const auto cp = Alphas(Q) / FourPi;
		const auto nf = NF(Q, Thresholds);
		Set<Operator> LO  {DISNCBasis{k,nf}, C2LO};
		Set<Operator> NLO {DISNCBasis{k,nf}, C2NLO};
		Set<Operator> NNLO{DISNCBasis{k,nf}, C2NNLO.at(nf)};
		return Charges(Q)[k-1] * ( LO + cp * ( NLO +  + cp * NNLO ) );
	      };
	  }
	// Define distribution function functions
	function<Set<Distribution>(double const&)> DistF2 = [=] (double const& Q) -> Set<Distribution>
	  {
	    const auto nf = NF(Q, Thresholds);
	    return Set<Distribution>{DISNCBasis{k,nf}, F2Map};
	  };
	// Initialize "Observable"
	F2.insert({k,Observable{C2f, DistF2}});
      }

    // Total structure function
    function<Set<Operator>(double const&)> C2f;
    if (PerturbativeOrder == 0)
      {
	C2f = [=] (double const& Q) -> Set<Operator>
	  {
	    const auto nf = NF(Q, Thresholds);
	    Set<Operator> LO{DISNCBasis{Charges(Q),nf}, C2LO};
	    return LO;
	  };
      }
    else if (PerturbativeOrder == 1)
      {
	C2f = [=] (double const& Q) -> Set<Operator>
	  {
	    const auto cp = Alphas(Q) / FourPi;
	    const auto nf = NF(Q, Thresholds);
	    Set<Operator> LO {DISNCBasis{Charges(Q),nf}, C2LO};
	    Set<Operator> NLO{DISNCBasis{Charges(Q),nf}, C2NLO};
	    return LO + cp * NLO;
	  };
      }
    else if (PerturbativeOrder == 2)
      {
	C2f = [=] (double const& Q) -> Set<Operator>
	  {
	    const auto cp = Alphas(Q) / FourPi;
	    const auto nf = NF(Q, Thresholds);
	    Set<Operator> LO  {DISNCBasis{Charges(Q),nf}, C2LO};
	    Set<Operator> NLO {DISNCBasis{Charges(Q),nf}, C2NLO};
	    Set<Operator> NNLO{DISNCBasis{Charges(Q),nf}, C2NNLO.at(nf)};
	    return LO + cp * ( NLO +  + cp * NNLO );
	  };
      }
    // Define distribution function functions
    function<Set<Distribution>(double const&)> DistF2 = [=] (double const& Q) -> Set<Distribution>
      {
	const auto nf = NF(Q, Thresholds);
	return Set<Distribution>{DISNCBasis{Charges(Q), nf}, F2Map};
      };
    F2.insert({0,Observable{C2f, DistF2}});

    t.stop();

    return F2;
  }

  //_____________________________________________________________________________
  unordered_map<int,Observable> FLBuildZM(Grid                                        const& g,
					  function<double(int const&, double const&)> const& InDistFunc,
					  vector<double>                              const& Thresholds,
					  int                                         const& PerturbativeOrder,
					  function<double(double const&)>             const& Alphas,
					  function<vector<double>(double const&)>     const& Charges,
					  double                                      const& IntEps)
  {
    cout << "Initializing FLBuildZM... ";
    Timer t;
    t.start();

    // Compute initial and final number of active flavours 
    // according to the vector of thresholds (it assumes that
    // the thresholds vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
	nfi++;

    // Allocate distributions
    unordered_map<int,Distribution> FLMap;
    FLMap.insert({0, DistributionFunction{g, InDistFunc, 0 }});
    FLMap.insert({1, DistributionFunction{g, InDistFunc, 1 }});
    FLMap.insert({2, DistributionFunction{g, InDistFunc, 3 }});
    FLMap.insert({3, DistributionFunction{g, InDistFunc, 5 }});
    FLMap.insert({4, DistributionFunction{g, InDistFunc, 7 }});
    FLMap.insert({5, DistributionFunction{g, InDistFunc, 9 }});
    FLMap.insert({6, DistributionFunction{g, InDistFunc, 11}});

    // Change sign to T3 to exchange "up" with "down"
    FLMap.at(2) *= -1;

    // ===============================================================
    const Operator Zero{g, Null{},     IntEps};

    // Coefficient functions for FL
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
	const Operator OL2t = OL2nsp + nf * OL2ps;
	unordered_map<int,Operator> CLNNLOnf;
	CLNNLOnf.insert({DISNCBasis::CNS, OL2nsp});
	CLNNLOnf.insert({DISNCBasis::CT,  OL2t});
	CLNNLOnf.insert({DISNCBasis::CG,  OL2g});
	CLNNLO.insert({nf,CLNNLOnf});
      }

    // Compute observables for FL
    unordered_map<int,Observable> FL;
    for (int k = 1; k <= 6; k++)
      {
	// Define coefficient function functions
	function<Set<Operator>(double const&)> CLf;
	if (PerturbativeOrder == 0)
	  {
	    CLf = [=] (double const& Q) -> Set<Operator>
	      {
		const auto nf = NF(Q, Thresholds);
		Set<Operator> LO{DISNCBasis{k,nf}, CLLO};
		return Charges(Q)[k-1] * LO;
	      };
	  }
	else if (PerturbativeOrder == 1)
	  {
	    CLf = [=] (double const& Q) -> Set<Operator>
	      {
		const auto cp = Alphas(Q) / FourPi;
		const auto nf = NF(Q, Thresholds);
		Set<Operator> NLO{DISNCBasis{k,nf}, CLNLO};
		return Charges(Q)[k-1] * cp * NLO;
	      };
	  }
	else if (PerturbativeOrder == 2)
	  {
	    CLf = [=] (double const& Q) -> Set<Operator>
	      {
		const auto cp = Alphas(Q) / FourPi;
		const auto nf = NF(Q, Thresholds);
		Set<Operator> NLO {DISNCBasis{k,nf}, CLNLO};
		Set<Operator> NNLO{DISNCBasis{k,nf}, CLNNLO.at(nf)};
		return Charges(Q)[k-1] * cp * ( NLO +  + cp * NNLO );
	      };
	  }
	// Define distribution function functions
	function<Set<Distribution>(double const&)> DistFL = [=] (double const& Q) -> Set<Distribution>
	  {
	    const auto nf = NF(Q, Thresholds);
	    return Set<Distribution>{DISNCBasis{k,nf}, FLMap};
	  };
	// Initialize "Observable"
	FL.insert({k,Observable{CLf, DistFL}});
      }

    // Total structure function
    function<Set<Operator>(double const&)> CLf;
    if (PerturbativeOrder == 0)
      {
	CLf = [=] (double const& Q) -> Set<Operator>
	  {
	    const auto nf = NF(Q, Thresholds);
	    Set<Operator> LO{DISNCBasis{Charges(Q),nf}, CLLO};
	    return LO;
	  };
      }
    else if (PerturbativeOrder == 1)
      {
	CLf = [=] (double const& Q) -> Set<Operator>
	  {
	    const auto cp = Alphas(Q) / FourPi;
	    const auto nf = NF(Q, Thresholds);
	    Set<Operator> NLO{DISNCBasis{Charges(Q),nf}, CLNLO};
	    return cp * NLO;
	  };
      }
    else if (PerturbativeOrder == 2)
      {
	CLf = [=] (double const& Q) -> Set<Operator>
	  {
	    const auto cp = Alphas(Q) / FourPi;
	    const auto nf = NF(Q, Thresholds);
	    Set<Operator> NLO {DISNCBasis{Charges(Q),nf}, CLNLO};
	    Set<Operator> NNLO{DISNCBasis{Charges(Q),nf}, CLNNLO.at(nf)};
	    return cp * ( NLO +  + cp * NNLO );
	  };
      }
    // Define distribution function functions
    function<Set<Distribution>(double const&)> DistFL = [=] (double const& Q) -> Set<Distribution>
      {
	const auto nf = NF(Q, Thresholds);
	return Set<Distribution>{DISNCBasis{Charges(Q), nf}, FLMap};
      };
    FL.insert({0,Observable{CLf, DistFL}});

    t.stop();

    return FL;
  }

  //_____________________________________________________________________________
  unordered_map<int,Observable> F3BuildZM(Grid                                        const& g,
					  function<double(int const&, double const&)> const& InDistFunc,
					  vector<double>                              const& Thresholds,
					  int                                         const& PerturbativeOrder,
					  function<double(double const&)>             const& Alphas,
					  function<vector<double>(double const&)>     const& Charges,
					  double                                      const& IntEps)
  {
    cout << "Initializing F3BuildZM... ";
    Timer t;
    t.start();

    // Compute initial and final number of active flavours 
    // according to the vector of thresholds (it assumes that
    // the thresholds vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
	nfi++;

    // Allocate distributions
    unordered_map<int,Distribution> F3Map;
    F3Map.insert({0, DistributionFunction{g, InDistFunc, 0 }});
    F3Map.insert({1, DistributionFunction{g, InDistFunc, 2 }});
    F3Map.insert({2, DistributionFunction{g, InDistFunc, 4 }});
    F3Map.insert({3, DistributionFunction{g, InDistFunc, 6 }});
    F3Map.insert({4, DistributionFunction{g, InDistFunc, 8 }});
    F3Map.insert({5, DistributionFunction{g, InDistFunc, 10}});
    F3Map.insert({6, DistributionFunction{g, InDistFunc, 12}});

    // Change sign to V3 to exchange "up" with "down"
    F3Map.at(2) *= -1;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // Coefficient functions for F3
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

    // Compute observables for F3
    unordered_map<int,Observable> F3;
    for (int k = 1; k <= 6; k++)
      {
	// Define coefficient function functions
	function<Set<Operator>(double const&)> C3f;
	if (PerturbativeOrder == 0)
	  {
	    C3f = [=] (double const& Q) -> Set<Operator>
	      {
		const auto nf = NF(Q, Thresholds);
		Set<Operator> LO{DISNCBasis{k,nf}, C3LO};
		return Charges(Q)[k-1] * LO;
	      };
	  }
	else if (PerturbativeOrder == 1)
	  {
	    C3f = [=] (double const& Q) -> Set<Operator>
	      {
		const auto cp = Alphas(Q) / FourPi;
		const auto nf = NF(Q, Thresholds);
		Set<Operator> NLO{DISNCBasis{k,nf}, C3NLO};
		return Charges(Q)[k-1] * cp * NLO;
	      };
	  }
	else if (PerturbativeOrder == 2)
	  {
	    C3f = [=] (double const& Q) -> Set<Operator>
	      {
		const auto cp = Alphas(Q) / FourPi;
		const auto nf = NF(Q, Thresholds);
		Set<Operator> NLO {DISNCBasis{k,nf}, C3NLO};
		Set<Operator> NNLO{DISNCBasis{k,nf}, C3NNLO.at(nf)};
		return Charges(Q)[k-1] * cp * ( NLO +  + cp * NNLO );
	      };
	  }
	// Define distribution function functions
	function<Set<Distribution>(double const&)> DistF3 = [=] (double const& Q) -> Set<Distribution>
	  {
	    const auto nf = NF(Q, Thresholds);
	    return Set<Distribution>{DISNCBasis{k,nf}, F3Map};
	  };
	// Initialize "Observable"
	F3.insert({k,Observable{C3f, DistF3}});
      }

    // Total structure function
    function<Set<Operator>(double const&)> C3f;
    if (PerturbativeOrder == 0)
      {
	C3f = [=] (double const& Q) -> Set<Operator>
	  {
	    const auto nf = NF(Q, Thresholds);
	    Set<Operator> LO{DISNCBasis{Charges(Q),nf}, C3LO};
	    return LO;
	  };
      }
    else if (PerturbativeOrder == 1)
      {
	C3f = [=] (double const& Q) -> Set<Operator>
	  {
	    const auto cp = Alphas(Q) / FourPi;
	    const auto nf = NF(Q, Thresholds);
	    Set<Operator> NLO{DISNCBasis{Charges(Q),nf}, C3NLO};
	    return cp * NLO;
	  };
      }
    else if (PerturbativeOrder == 2)
      {
	C3f = [=] (double const& Q) -> Set<Operator>
	  {
	    const auto cp = Alphas(Q) / FourPi;
	    const auto nf = NF(Q, Thresholds);
	    Set<Operator> NLO {DISNCBasis{Charges(Q),nf}, C3NLO};
	    Set<Operator> NNLO{DISNCBasis{Charges(Q),nf}, C3NNLO.at(nf)};
	    return cp * ( NLO +  + cp * NNLO );
	  };
      }
    // Define distribution function functions
    function<Set<Distribution>(double const&)> DistF3 = [=] (double const& Q) -> Set<Distribution>
      {
	const auto nf = NF(Q, Thresholds);
	return Set<Distribution>{DISNCBasis{Charges(Q), nf}, F3Map};
      };
    F3.insert({0,Observable{C3f, DistF3}});

    t.stop();

    return F3;
  }

}
