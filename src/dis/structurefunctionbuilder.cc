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
#include "apfel/massivecoefficientfunctions.h"
#include "apfel/tabulateobject.h"

#include <map>

using namespace std;

namespace apfel {

  //_____________________________________________________________________________
  function<StructureFunctionObjects(double const&)> InitializeF2NCObjectsMassive(Grid           const& g,
										 vector<double> const& Masses,
										 double         const& IntEps,
										 int            const& nxi,
										 double         const& ximin,
										 double         const& ximax,
										 int            const& intdeg,
										 double         const& lambda)
  {
    cout << "Initializing StructureFunctionObjects for F2 NC Massive... \n";
    Timer t;
    t.start();

    // ===============================================================
    // Massive coefficient functions
    // NLO
    const TabulateObject<Operator> TabO21g{[=,&g] (double const& xi) -> Operator
	{
	  const double teta = 1 / ( 1 + 4 / xi );
	  return Operator{g, Cm21gNC{teta}};
	}, nxi, ximin, ximax, intdeg, {}, lambda};

    // NNLO
    const TabulateObject<Operator> TabO22ns{[=,&g] (double const& xi) -> Operator
	{
	  const double teta = 1 / ( 1 + 4 / xi );
	  return Operator{g, Cm22nsNC{teta}};
	}, nxi, ximin, ximax, intdeg, {}, lambda};
    const auto fO22s = [=,&g] (double const& xi) -> Operator
      {
	const double teta = 1 / ( 1 + 4 / xi );
	const Operator O22psc{g, Cm22psNC{teta}};
	const Operator O22psl{g, Cm22barpsNC{teta}};
	return 6 * ( O22psc + log(xi) * O22psl );
      };
    const TabulateObject<Operator> TabO22s{fO22s, nxi, ximin, ximax, intdeg, {}, lambda};
    const auto fO22g = [=,&g] (double const& xi) -> Operator
      {
	const double teta = 1 / ( 1 + 4 / xi );
	const Operator O22gc{g, Cm22gNC{teta}};
	const Operator O22gl{g, Cm22bargNC{teta}};
	return O22gc + log(xi) * O22gl;
      };
    const TabulateObject<Operator> TabO22g{fO22g, nxi, ximin, ximax, intdeg, {}, lambda};

    // Zero-mass coefficient functions have to be initialized anyway.
    StructureFunctionObjects FObjZM = InitializeF2NCObjectsZM(g,IntEps);

    // Null set of operators for each component as an utility.
    const Operator Zero{g, Null{}, IntEps};
    map<int,Operator> None;
    None.insert({CNS, Zero});
    None.insert({CS,  Zero});
    None.insert({CG,  Zero});
    map<int,Set<Operator>> Not;
    for (auto it = FObjZM.ConvBasis.begin(); it != FObjZM.ConvBasis.end(); ++it)
      Not.insert({it->first,Set<Operator>{it->second, None}});

    // Size of the mass vector
    const int nm = Masses.size();

    // Define object of the structure containing the DglapObjects
    auto F2Obj = [=] (double const& Q) -> StructureFunctionObjects
      {
	const double Q2  = Q * Q;
	StructureFunctionObjects FObj;
	FObj.skip         = FObjZM.skip;
	FObj.ConvBasis    = FObjZM.ConvBasis;
	FObj.ConvBasisTot = FObjZM.ConvBasisTot;
	for (auto it = FObjZM.ConvBasis.begin(); it != FObjZM.ConvBasis.end(); ++it)
	  {
	    const int k = it->first;
	    // Now include the massive bits whenever needed, that is
	    // when there is a corresponding entry in the mass vector
	    // and if the respective value is biggere than zero.
	    if (k <= nm && Masses[k-1] > 0)
	      {
		// Compute value of xi.
		const double M2 = Masses[k-1] * Masses[k-1];
		const double xi = Q2 / M2;

		// Set LO non-singlet and singlet coefficient
		// functions to zero (gluon is already zero).
		FObj.C0.insert({k,Not.at(k)});

		// Now insert NLO
		map<int,Operator> NLO;
		NLO.insert({CNS, Zero});
		NLO.insert({CS,  Zero});
		NLO.insert({CG,  TabO21g.Evaluate(xi)});
		FObj.C1.insert({k,Set<Operator>{it->second,NLO}});

		// Now insert NNLO
		map<int,Set<Operator>> NNLO;
		for (auto inf = FObjZM.ConvBasis.begin(); inf != FObjZM.ConvBasis.end(); ++inf)
		  {
		    map<int,Operator> mNNLO;
		    mNNLO.insert({CNS, Zero});
		    mNNLO.insert({CS,  TabO22s.Evaluate(xi)});
		    mNNLO.insert({CG,  TabO22g.Evaluate(xi)});
		    NNLO.insert({inf->first,Set<Operator>{it->second, mNNLO}});
		  }
		FObj.C2.insert({k,NNLO});
	      }
	    else
	      {
		FObj.C0.insert({k,FObjZM.C0.at(k)});
		FObj.C1.insert({k,FObjZM.C1.at(k)});
		FObj.C2.insert({k,FObjZM.C2.at(k)});
	      }
	  }
	return FObj;
      };
    t.stop();

    return F2Obj;
  }

  //_____________________________________________________________________________
  StructureFunctionObjects InitializeF2NCObjectsZM(Grid const& g, double const& IntEps)
  {
    cout << "Initializing StructureFunctionObjects for F2 NC Zero Mass... ";
    Timer t;
    t.start();

    // Set initial and final number of active flavours. Precompute
    // objects for all numbers of flavours because it cheap enough.
    int nfi = 1;
    int nff = 6;

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
    for (int nf = nfi; nf <= nff; nf++)
      {
	const Operator O22nsp{g, C22nsp{nf}, IntEps};
	const Operator O22t = O22nsp + 6 * O22ps;
	map<int,Operator> C2NNLOnf;
	C2NNLOnf.insert({CNS, O22nsp});
	C2NNLOnf.insert({CS,  O22t});
	C2NNLOnf.insert({CG,  O22g});
	C2NNLO.insert({nf,C2NNLOnf});
      }

    // Define object of the structure containing the DglapObjects
    StructureFunctionObjects F2Obj;
    F2Obj.skip = {2,4,6,8,10,12};
    for (int k = nfi; k <= nff; k++)
      {
	F2Obj.ConvBasis.insert({k,DISNCBasis{k}});
	F2Obj.C0.insert({k,Set<Operator>{F2Obj.ConvBasis.at(k), C2LO}});
	F2Obj.C1.insert({k,Set<Operator>{F2Obj.ConvBasis.at(k), C2NLO}});
	map<int,Set<Operator>> NNLO;
	for (int nf = nfi; nf <= nff; nf++)
	  NNLO.insert({nf,Set<Operator>{F2Obj.ConvBasis.at(k), C2NNLO.at(nf)}});
	F2Obj.C2.insert({k,NNLO});
      }
    F2Obj.ConvBasisTot = [=] (vector<double> const& Ch) -> ConvolutionMap{ return DISNCBasis{Ch}; };
    t.stop();

    return F2Obj;
  }

  //_____________________________________________________________________________
  StructureFunctionObjects InitializeFLNCObjectsZM(Grid const& g, double const& IntEps)
  {
    cout << "Initializing StructureFunctionObjects for FL NC Zero Mass... ";
    Timer t;
    t.start();

    // Set initial and final number of active flavours. Precompute
    // objects for all numbers of flavours because it cheap enough.
    int nfi = 1;
    int nff = 6;

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
    for (int nf = nfi; nf <= nff; nf++)
      {
	const Operator OL2nsp{g, CL2nsp{nf}, IntEps};
	const Operator OL2t = OL2nsp + 6 * OL2ps;
	map<int,Operator> CLNNLOnf;
	CLNNLOnf.insert({CNS, OL2nsp});
	CLNNLOnf.insert({CS,  OL2t});
	CLNNLOnf.insert({CG,  OL2g});
	CLNNLO.insert({nf,CLNNLOnf});
      }

    // Define object of the structure containing the DglapObjects
    StructureFunctionObjects FLObj;
    FLObj.skip = {2,4,6,8,10,12};
    for (int k = nfi; k <= nff; k++)
      {
	FLObj.ConvBasis.insert({k,DISNCBasis{k}});
	FLObj.C0.insert({k,Set<Operator>{FLObj.ConvBasis.at(k), CLLO}});
	FLObj.C1.insert({k,Set<Operator>{FLObj.ConvBasis.at(k), CLNLO}});
	map<int,Set<Operator>> NNLO;
	for (int nf = nfi; nf <= nff; nf++)
	  NNLO.insert({nf,Set<Operator>{FLObj.ConvBasis.at(k), CLNNLO.at(nf)}});
	FLObj.C2.insert({k,NNLO});
      }
    FLObj.ConvBasisTot = [=] (vector<double> const& Ch) -> ConvolutionMap{ return DISNCBasis{Ch}; };
    t.stop();

    return FLObj;
  }

  //_____________________________________________________________________________
  StructureFunctionObjects InitializeF3NCObjectsZM(Grid const& g, double const& IntEps)
  {
    cout << "Initializing StructureFunctionObjects for F3 NC Zero Mass... ";
    Timer t;
    t.start();

    // Set initial and final number of active flavours. Precompute
    // objects for all numbers of flavours because it cheap enough.
    int nfi = 1;
    int nff = 6;

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
    for (int nf = nfi; nf <= nff; nf++)
      {
	const Operator O32nsm{g, C32nsm{nf}, IntEps};
	const Operator O32t = O32nsm;
	map<int,Operator> C3NNLOnf;
	C3NNLOnf.insert({CNS, O32nsm});
	C3NNLOnf.insert({CS,  O32t});
	C3NNLOnf.insert({CG,  Zero});
	C3NNLO.insert({nf,C3NNLOnf});
      }

    // Define object of the structure containing the DglapObjects
    StructureFunctionObjects F3Obj;
    F3Obj.skip = {1,3,5,7,9,11};
    for (int k = nfi; k <= nff; k++)
      {
	F3Obj.ConvBasis.insert({k,DISNCBasis{k}});
	F3Obj.C0.insert({k,Set<Operator>{F3Obj.ConvBasis.at(k), C3LO}});
	F3Obj.C1.insert({k,Set<Operator>{F3Obj.ConvBasis.at(k), C3NLO}});
	map<int,Set<Operator>> NNLO;
	for (int nf = nfi; nf <= nff; nf++)
	  NNLO.insert({nf,Set<Operator>{F3Obj.ConvBasis.at(k), C3NNLO.at(nf)}});
	F3Obj.C2.insert({k,NNLO});
      }
    F3Obj.ConvBasisTot = [=] (vector<double> const& Ch) -> ConvolutionMap{ return DISNCBasis{Ch}; };
    t.stop();

    return F3Obj;
  }

  //_____________________________________________________________________________
  StructureFunctionObjects InitializeF2CCPlusObjectsZM(Grid const& g, double const& IntEps)
  {
    cout << "Initializing StructureFunctionObjects for ( F2(nu) + F2(nubar) ) / 2 Zero Mass... ";
    Timer t;
    t.start();

    // Set initial and final number of active channels.
    int chi = 1;
    int chf = 9;

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
    for (int ch = chi; ch <= chf; ch++)
      {
	const Operator O22nsp{g, C22nsp{ch}, IntEps};
	const Operator O22t = O22nsp + 6 * O22ps;
	map<int,Operator> C2NNLOch;
	C2NNLOch.insert({CNS, O22nsp});
	C2NNLOch.insert({CS,  O22t});
	C2NNLOch.insert({CG,  O22g});
	C2NNLO.insert({ch,C2NNLOch});
      }

    // Define object of the structure containing the DglapObjects
    StructureFunctionObjects F2Obj;
    F2Obj.skip = {2,4,6,8,10,12};
    for (int k = 1; k <= 6; k++)
      {
	F2Obj.ConvBasis.insert({k,DISCCBasis{k,false}});
	F2Obj.C0.insert({k,Set<Operator>{F2Obj.ConvBasis.at(k), C2LO}});
	F2Obj.C1.insert({k,Set<Operator>{F2Obj.ConvBasis.at(k), C2NLO}});
	map<int,Set<Operator>> NNLO;
	for (int ch = chi; ch <= chf; ch++)
	  NNLO.insert({ch,Set<Operator>{F2Obj.ConvBasis.at(k), C2NNLO.at(ch)}});
	F2Obj.C2.insert({k,NNLO});
      }
    F2Obj.ConvBasisTot = [=] (vector<double> const& CKM) -> ConvolutionMap{ return DISCCBasis{CKM,false}; };
    t.stop();

    return F2Obj;
  }

  //_____________________________________________________________________________
  StructureFunctionObjects InitializeF2CCMinusObjectsZM(Grid const& g, double const& IntEps)
  {
    cout << "Initializing StructureFunctionObjects for ( F2(nu) - F2(nubar) ) / 2 Zero Mass... ";
    Timer t;
    t.start();

    // Set initial and final number of active channels.
    int chi = 1;
    int chf = 9;

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
    for (int ch = chi; ch <= chf; ch++)
      {
	const Operator O22nsm{g, C22nsm{ch}, IntEps};
	map<int,Operator> C2NNLOch;
	C2NNLOch.insert({CNS, O22nsm});
	C2NNLOch.insert({CS,  Zero});
	C2NNLOch.insert({CG,  Zero});
	C2NNLO.insert({ch,C2NNLOch});
      }

    // Define object of the structure containing the DglapObjects
    StructureFunctionObjects F2Obj;
    F2Obj.skip = {0,1,2,3,5,7,9,11};
    for (int k = 1; k <= 6; k++)
      {
	F2Obj.ConvBasis.insert({k,DISCCBasis{k,false}});
	F2Obj.C0.insert({k,Set<Operator>{F2Obj.ConvBasis.at(k), C2LO}});
	F2Obj.C1.insert({k,Set<Operator>{F2Obj.ConvBasis.at(k), C2NLO}});
	map<int,Set<Operator>> NNLO;
	for (int ch = chi; ch <= chf; ch++)
	  NNLO.insert({ch,Set<Operator>{F2Obj.ConvBasis.at(k), C2NNLO.at(ch)}});
	F2Obj.C2.insert({k,NNLO});
      }
    F2Obj.ConvBasisTot = [=] (vector<double> const& CKM) -> ConvolutionMap{ return DISCCBasis{CKM,false}; };
    t.stop();

    return F2Obj;
  }

  //_____________________________________________________________________________
  StructureFunctionObjects InitializeFLCCPlusObjectsZM(Grid const& g, double const& IntEps)
  {
    cout << "Initializing StructureFunctionObjects for ( FL(nu) + FL(nubar) ) / 2 Zero Mass... ";
    Timer t;
    t.start();

    // Set initial and final number of active channels.
    int chi = 1;
    int chf = 9;

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
    for (int ch = chi; ch <= chf; ch++)
      {
	const Operator OL2nsp{g, CL2nsp{ch}, IntEps};
	const Operator OL2t = OL2nsp + 6 * OL2ps;
	map<int,Operator> CLNNLOch;
	CLNNLOch.insert({CNS, OL2nsp});
	CLNNLOch.insert({CS,  OL2t});
	CLNNLOch.insert({CG,  OL2g});
	CLNNLO.insert({ch,CLNNLOch});
      }

    // Define object of the structure containing the DglapObjects
    StructureFunctionObjects FLObj;
    FLObj.skip = {2,4,6,8,10,12};
    for (int k = 1; k <= 6; k++)
      {
	FLObj.ConvBasis.insert({k,DISCCBasis{k,false}});
	FLObj.C0.insert({k,Set<Operator>{FLObj.ConvBasis.at(k), CLLO}});
	FLObj.C1.insert({k,Set<Operator>{FLObj.ConvBasis.at(k), CLNLO}});
	map<int,Set<Operator>> NNLO;
	for (int ch = chi; ch <= chf; ch++)
	  NNLO.insert({ch,Set<Operator>{FLObj.ConvBasis.at(k), CLNNLO.at(ch)}});
	FLObj.C2.insert({k,NNLO});
      }
    FLObj.ConvBasisTot = [=] (vector<double> const& CKM) -> ConvolutionMap{ return DISCCBasis{CKM,false}; };
    t.stop();

    return FLObj;
  }

  //_____________________________________________________________________________
  StructureFunctionObjects InitializeFLCCMinusObjectsZM(Grid const& g, double const& IntEps)
  {
    cout << "Initializing StructureFunctionObjects for ( FL(nu) - FL(nubar) ) / 2 Zero Mass... ";
    Timer t;
    t.start();

    // Set initial and final number of active channels.
    int chi = 1;
    int chf = 9;

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
    for (int ch = chi; ch <= chf; ch++)
      {
	const Operator OL2nsm{g, CL2nsm{ch}, IntEps};
	map<int,Operator> CLNNLOch;
	CLNNLOch.insert({CNS, OL2nsm});
	CLNNLOch.insert({CS,  Zero});
	CLNNLOch.insert({CG,  Zero});
	CLNNLO.insert({ch,CLNNLOch});
      }

    // Define object of the structure containing the DglapObjects
    StructureFunctionObjects FLObj;
    FLObj.skip = {0,1,2,3,5,7,9,11};
    for (int k = 1; k <= 6; k++)
      {
	FLObj.ConvBasis.insert({k,DISCCBasis{k,false}});
	FLObj.C0.insert({k,Set<Operator>{FLObj.ConvBasis.at(k), CLLO}});
	FLObj.C1.insert({k,Set<Operator>{FLObj.ConvBasis.at(k), CLNLO}});
	map<int,Set<Operator>> NNLO;
	for (int ch = chi; ch <= chf; ch++)
	  NNLO.insert({ch,Set<Operator>{FLObj.ConvBasis.at(k), CLNNLO.at(ch)}});
	FLObj.C2.insert({k,NNLO});
      }
    FLObj.ConvBasisTot = [=] (vector<double> const& CKM) -> ConvolutionMap{ return DISCCBasis{CKM,false}; };
    t.stop();

    return FLObj;
  }

  //_____________________________________________________________________________
  StructureFunctionObjects InitializeF3CCPlusObjectsZM(Grid const& g, double const& IntEps)
  {
    cout << "Initializing StructureFunctionObjects for ( F3(nu) + F3(nubar) ) / 2 Zero Mass... ";
    Timer t;
    t.start();

    // Set initial and final number of active channels.
    int chi = 1;
    int chf = 9;

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
    for (int ch = chi; ch <= chf; ch++)
      {
	const Operator O32nsp{g, C32nsp{ch}, IntEps};
	const Operator O32t = O32nsp;
	map<int,Operator> C3NNLOch;
	C3NNLOch.insert({CNS, O32nsp});
	C3NNLOch.insert({CS,  O32t});
	C3NNLOch.insert({CG,  Zero});
	C3NNLO.insert({ch,C3NNLOch});
      }

    // Define object of the structure containing the DglapObjects
    StructureFunctionObjects F3Obj;
    F3Obj.skip = {0,1,2,4,6,8,10,12};
    for (int k = 1; k <= 6; k++)
      {
	F3Obj.ConvBasis.insert({k,DISCCBasis{k,true}});
	F3Obj.C0.insert({k,Set<Operator>{F3Obj.ConvBasis.at(k), C3LO}});
	F3Obj.C1.insert({k,Set<Operator>{F3Obj.ConvBasis.at(k), C3NLO}});
	map<int,Set<Operator>> NNLO;
	for (int ch = chi; ch <= chf; ch++)
	  NNLO.insert({ch,Set<Operator>{F3Obj.ConvBasis.at(k), C3NNLO.at(ch)}});
	F3Obj.C2.insert({k,NNLO});
      }
    F3Obj.ConvBasisTot = [=] (vector<double> const& CKM) -> ConvolutionMap{ return DISCCBasis{CKM,true}; };
    t.stop();

    return F3Obj;
  }

  //_____________________________________________________________________________
  StructureFunctionObjects InitializeF3CCMinusObjectsZM(Grid const& g, double const& IntEps)
  {
    cout << "Initializing StructureFunctionObjects for ( F3(nu) - F3(nubar) ) / 2 Zero Mass... ";
    Timer t;
    t.start();

    // Set initial and final number of active channels.
    int chi = 1;
    int chf = 9;

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
    for (int ch = chi; ch <= chf; ch++)
      {
	const Operator O32nsm{g, C32nsm{ch}, IntEps};
	map<int,Operator> C3NNLOch;
	C3NNLOch.insert({CNS, O32nsm});
	C3NNLOch.insert({CS,  O32nsm});
	C3NNLOch.insert({CG,  Zero});
	C3NNLO.insert({ch,C3NNLOch});
      }

    // Define object of the structure containing the DglapObjects
    StructureFunctionObjects F3Obj;
    F3Obj.skip = {0,1,3,5,7,9,11};
    for (int k = 1; k <= 6; k++)
      {
	F3Obj.ConvBasis.insert({k,DISCCBasis{k,true}});
	F3Obj.C0.insert({k,Set<Operator>{F3Obj.ConvBasis.at(k), C3LO}});
	F3Obj.C1.insert({k,Set<Operator>{F3Obj.ConvBasis.at(k), C3NLO}});
	map<int,Set<Operator>> NNLO;
	for (int ch = chi; ch <= chf; ch++)
	  NNLO.insert({ch,Set<Operator>{F3Obj.ConvBasis.at(k), C3NNLO.at(ch)}});
	F3Obj.C2.insert({k,NNLO});
      }
    F3Obj.ConvBasisTot = [=] (vector<double> const& CKM) -> ConvolutionMap{ return DISCCBasis{CKM,true}; };
    t.stop();

    return F3Obj;
  }

  //_____________________________________________________________________________
  map<int,Observable> BuildStructureFunctions(StructureFunctionObjects                                const& FObj,
					      function<map<int,double>(double const&, double const&)> const& InDistFunc,
					      vector<double>                                          const& Thresholds,
					      int                                                     const& PerturbativeOrder,
					      function<double(double const&)>                         const& Alphas,
					      function<vector<double>(double const&)>                 const& Couplings)
  {
    // Get grid.
    Grid const& g = FObj.C0.at(1).at(0).GetGrid();

    // Cycle over the key of the convolution basis map.
    map<int,Observable> F;
    for (auto it = FObj.ConvBasis.begin(); it != FObj.ConvBasis.end(); ++it)
      {
	// Structure function index.
	const int k = it->first;

	// Define coefficient function functions.
	function<Set<Operator>(double const&)> Cf;
	if (PerturbativeOrder == 0)
	  {
	    Cf = [=] (double const& Q) -> Set<Operator>
	      {
		return Couplings(Q)[k-1] * FObj.C0.at(k);
	      };
	  }
	else if (PerturbativeOrder == 1)
	  {
	    Cf = [=] (double const& Q) -> Set<Operator>
	      {
		const auto cp = Alphas(Q) / FourPi;
		return Couplings(Q)[k-1] * ( FObj.C0.at(k) + cp * FObj.C1.at(k) );
	      };
	  }
	else if (PerturbativeOrder == 2)
	  {
	    Cf = [=] (double const& Q) -> Set<Operator>
	      {
		const auto cp = Alphas(Q) / FourPi;
		const auto nf = NF(Q, Thresholds);
		return Couplings(Q)[k-1] * ( FObj.C0.at(k) + cp * ( FObj.C1.at(k) +  + cp * FObj.C2.at(k).at(nf) ) );
	      };
	  }

	// Define distribution function functions.
	const auto DistF = [=,&g] (double const& Q) -> Set<Distribution>
	  {
	    return Set<Distribution>{it->second, DistributionMap(g, InDistFunc, Q, FObj.skip)};
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
	    const auto basis = FObj.ConvBasisTot(Couplings(Q));
	    Set<Operator> LO{basis, FObj.C0.at(1).GetObjects()};
	    return LO;
	  };
      }
    else if (PerturbativeOrder == 1)
      {
	Cf = [=] (double const& Q) -> Set<Operator>
	  {
	    const auto cp = Alphas(Q) / FourPi;
	    const auto basis = FObj.ConvBasisTot(Couplings(Q));
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
	    const auto basis = FObj.ConvBasisTot(Couplings(Q));
	    Set<Operator> LO  {basis, FObj.C0.at(1).GetObjects()};
	    Set<Operator> NLO {basis, FObj.C1.at(1).GetObjects()};
	    Set<Operator> NNLO{basis, FObj.C2.at(1).at(nf).GetObjects()};
	    return LO + cp * ( NLO +  + cp * NNLO );
	  };
      }

    // Define distribution function functions.
    const auto DistF = [=,&g] (double const& Q) -> Set<Distribution>
      {
	return Set<Distribution>{FObj.ConvBasisTot(Couplings(Q)), DistributionMap(g, InDistFunc, Q, FObj.skip)};
      };
    F.insert({0,Observable{Cf, DistF}});

    return F;
  }

  //_____________________________________________________________________________
  map<int,Observable> BuildStructureFunctions(StructureFunctionObjects                                   const& FObj,
					      function<double(int const&, double const&, double const&)> const& InDistFunc,
					      vector<double>                                             const& Thresholds,
					      int                                                        const& PerturbativeOrder,
					      function<double(double const&)>                            const& Alphas,
					      function<vector<double>(double const&)>                    const& Couplings)
  {
    const auto InDistFuncMap = [=] (double const& x, double const& Q) -> map<int,double>
      {
	map<int,double> DistMap;
	for (int i = 0; i <= 12; i++)
	  DistMap.insert({i,InDistFunc(i, x, Q)});
	return DistMap;
      };
    return BuildStructureFunctions(FObj, InDistFuncMap, Thresholds, PerturbativeOrder, Alphas, Couplings);
  }

  //_____________________________________________________________________________
  map<int,Observable> BuildStructureFunctions(function<StructureFunctionObjects(double const&)>       const& FObj,
					      function<map<int,double>(double const&, double const&)> const& InDistFunc,
					      vector<double>                                          const& Thresholds,
					      int                                                     const& PerturbativeOrder,
					      function<double(double const&)>                         const& Alphas,
					      function<vector<double>(double const&)>                 const& Couplings)
  {
    // Call FObj at zero energy to use it for those quantities that do
    // not depend on Q.
    const StructureFunctionObjects FObj0 = FObj(0);

    // Get grid.
    Grid const& g = FObj0.C0.at(1).at(0).GetGrid();

    // Cycle over the key of the convolution basis map.
    map<int,Observable> F;
    for (auto it = FObj0.ConvBasis.begin(); it != FObj0.ConvBasis.end(); ++it)
      {
	// Structure function index.
	const int k = it->first;

	// Define coefficient function functions.
	function<Set<Operator>(double const&)> Cf;
	if (PerturbativeOrder == 0)
	  {
	    Cf = [=] (double const& Q) -> Set<Operator>
	      {
		return Couplings(Q)[k-1] * FObj(Q).C0.at(k);
	      };
	  }
	else if (PerturbativeOrder == 1)
	  {
	    Cf = [=] (double const& Q) -> Set<Operator>
	      {
		const auto cp = Alphas(Q) / FourPi;
		const StructureFunctionObjects FObjQ = FObj(Q);
		return Couplings(Q)[k-1] * ( FObjQ.C0.at(k) + cp * FObjQ.C1.at(k) );
	      };
	  }
	else if (PerturbativeOrder == 2)
	  {
	    Cf = [=] (double const& Q) -> Set<Operator>
	      {
		const auto cp = Alphas(Q) / FourPi;
		const auto nf = NF(Q, Thresholds);
		const StructureFunctionObjects FObjQ = FObj(Q);
		return Couplings(Q)[k-1] * ( FObjQ.C0.at(k) + cp * ( FObjQ.C1.at(k) +  + cp * FObjQ.C2.at(k).at(nf) ) );
	      };
	  }

	// Define distribution function functions.
	const auto DistF = [=,&g] (double const& Q) -> Set<Distribution>
	  {
	    const StructureFunctionObjects FObjQ = FObj(Q);
	    return Set<Distribution>{FObjQ.ConvBasis.at(k), DistributionMap(g, InDistFunc, Q, FObjQ.skip)};
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
	    const StructureFunctionObjects FObjQ = FObj(Q);
	    const auto basis = FObjQ.ConvBasisTot(Couplings(Q));
	    Set<Operator> LO{basis, FObjQ.C0.at(1).GetObjects()};
	    return LO;
	  };
      }
    else if (PerturbativeOrder == 1)
      {
	Cf = [=] (double const& Q) -> Set<Operator>
	  {
	    const auto cp = Alphas(Q) / FourPi;
	    const StructureFunctionObjects FObjQ = FObj(Q);
	    const auto basis = FObjQ.ConvBasisTot(Couplings(Q));
	    Set<Operator> LO {basis, FObjQ.C0.at(1).GetObjects()};
	    Set<Operator> NLO{basis, FObjQ.C1.at(1).GetObjects()};
	    return LO + cp * NLO;
	  };
      }
    else if (PerturbativeOrder == 2)
      {
	Cf = [=] (double const& Q) -> Set<Operator>
	  {
	    const auto cp = Alphas(Q) / FourPi;
	    const auto nf = NF(Q, Thresholds);
	    const StructureFunctionObjects FObjQ = FObj(Q);
	    const auto basis = FObjQ.ConvBasisTot(Couplings(Q));
	    Set<Operator> LO  {basis, FObjQ.C0.at(1).GetObjects()};
	    Set<Operator> NLO {basis, FObjQ.C1.at(1).GetObjects()};
	    Set<Operator> NNLO{basis, FObjQ.C2.at(1).at(nf).GetObjects()};
	    return LO + cp * ( NLO +  + cp * NNLO );
	  };
      }

    // Define distribution function functions.
    const auto DistF = [=,&g] (double const& Q) -> Set<Distribution>
      {
	const StructureFunctionObjects FObjQ = FObj(Q);
	return Set<Distribution>{FObjQ.ConvBasisTot(Couplings(Q)), DistributionMap(g, InDistFunc, Q, FObjQ.skip)};
      };
    F.insert({0,Observable{Cf, DistF}});

    return F;
  }

  //_____________________________________________________________________________
  map<int,Observable> BuildStructureFunctions(function<StructureFunctionObjects(double const&)>          const& FObj,
					      function<double(int const&, double const&, double const&)> const& InDistFunc,
					      vector<double>                                             const& Thresholds,
					      int                                                        const& PerturbativeOrder,
					      function<double(double const&)>                            const& Alphas,
					      function<vector<double>(double const&)>                    const& Couplings)
  {
    const auto InDistFuncMap = [=] (double const& x, double const& Q) -> map<int,double>
      {
	map<int,double> DistMap;
	for (int i = 0; i <= 12; i++)
	  DistMap.insert({i,InDistFunc(i, x, Q)});
	return DistMap;
      };
    return BuildStructureFunctions(FObj, InDistFuncMap, Thresholds, PerturbativeOrder, Alphas, Couplings);
  }
}
