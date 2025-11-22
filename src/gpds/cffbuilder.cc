//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/cffbuilder.h"
#include "apfel/gpdbuilder.h"
#include "apfel/timer.h"
#include "apfel/tools.h"
#include "apfel/zeromasscff.h"

namespace apfel
{
  //_____________________________________________________________________________
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeImCFF1NCObjectsZM(Grid                const& g,
                                                                                                                 std::vector<double> const& Thresholds,
                                                                                                                 double              const& xi,
                                                                                                                 double              const& IntEps)
  {
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = InitializeGpdObjects(g, Thresholds, xi);

    report("Initializing StructureFunctionObjects for the imaginary part of Compton form factor F1 Zero Mass at xi = " + std::to_string(xi) + "... ");
    Timer t;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps, true};
    const Operator Zero{g, Null{},     IntEps, true};

    // LO
    std::map<int, Operator> C1LO;
    C1LO.insert({DISNCBasis::CNS, M_PI * Id});
    C1LO.insert({DISNCBasis::CS,  M_PI * Id});
    C1LO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> C1NLO;
    const Operator O11ns{g, ImCFF11ns{}, IntEps, true};
    const Operator O11g {g, ImCFF11g{},  IntEps, true};
    C1NLO.insert({DISNCBasis::CNS, O11ns});
    C1NLO.insert({DISNCBasis::CS,  O11ns});
    C1NLO.insert({DISNCBasis::CG,  O11g});

    // Vector of distributions to skip
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

    // Define object of the structure containing the DglapObjects
    const auto F1Obj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
    {
      // Determine number of active flavours.
      const int nf = NF(Q, Thresholds);

      // Effective charges. The charges of the components with mh >
      // Q are set to zero.
      std::vector<double> EffCh;
      for (int k = 1; k <= 6; k++)
        EffCh.push_back((k > nf ? 0 : Ch[k-1]));

      // Fill in structure function object
      StructureFunctionObjects FObj;
      FObj.nf   = nf;
      FObj.P    = PDFObj.at(nf);
      FObj.skip = skip;
      // Single structure function components.
      for (int k = 0; k <= 6; k++)
        {
// *INDENT-OFF*
          FObj.ConvBasis.insert({k, (k == 0 ? DISNCBasis{EffCh} : DISNCBasis{k, EffCh[k-1]})});
// *INDENT-ON*
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C1LO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C1NLO}});
        }
      return FObj;
    };
    t.stop();

    return F1Obj;
  }

  //_____________________________________________________________________________
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeReCFF1NCObjectsZM(Grid                const& g,
                                                                                                                 std::vector<double> const& Thresholds,
                                                                                                                 double              const& xi,
                                                                                                                 double              const& IntEps)
  {
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = InitializeGpdObjects(g, Thresholds, xi);

    report("Initializing StructureFunctionObjects for the real part of Compton form factor F1 Zero Mass at xi = " + std::to_string(xi) + "... ");
    Timer t;

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps, true};

    // LO
    std::map<int, Operator> C1LO;
    const Operator O10ns{g, ReCFF10ns{}, IntEps, true};
    C1LO.insert({DISNCBasis::CNS, O10ns});
    C1LO.insert({DISNCBasis::CS,  O10ns});
    C1LO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> C1NLO;
    const Operator O11ns{g, ReCFF11ns{}, IntEps, true};
    const Operator O11g {g, ReCFF11g{},  IntEps, true};
    C1NLO.insert({DISNCBasis::CNS, O11ns});
    C1NLO.insert({DISNCBasis::CS,  O11ns});
    C1NLO.insert({DISNCBasis::CG,  O11g});

    // Vector of distributions to skip
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

    // Define object of the structure containing the DglapObjects
    const auto F1Obj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
    {
      // Determine number of active flavours.
      const int nf = NF(Q, Thresholds);

      // Effective charges. The charges of the components with mh >
      // Q are set to zero.
      std::vector<double> EffCh;
      for (int k = 1; k <= 6; k++)
        EffCh.push_back((k > nf ? 0 : Ch[k-1]));

      // Fill in structure function object
      StructureFunctionObjects FObj;
      FObj.nf   = nf;
      FObj.P    = PDFObj.at(nf);
      FObj.skip = skip;
      // Single structure function components.
      for (int k = 0; k <= 6; k++)
        {
// *INDENT-OFF*
          FObj.ConvBasis.insert({k, (k == 0 ? DISNCBasis{EffCh} : DISNCBasis{k, EffCh[k-1]})});
// *INDENT-ON*
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C1LO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C1NLO}});
        }
      return FObj;
    };
    t.stop();

    return F1Obj;
  }
}
