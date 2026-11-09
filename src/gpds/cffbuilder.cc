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
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeImCFF2NCObjectsZMSpaceLike(Grid                const& g,
                                                                                                                          std::vector<double> const& Thresholds,
                                                                                                                          double              const& xi,
                                                                                                                          double              const& IntEps)
  {
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = InitializeGpdObjects(g, Thresholds, xi);

    report("Initializing StructureFunctionObjects for the imaginary part of Compton form factor F2 Zero Mass at xi = " + std::to_string(xi) + " (Q2 > 0)... ");
    Timer t;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps, true};
    const Operator Zero{g, Null{},     IntEps, true};

    // LO
    std::map<int, Operator> C2LO;
    C2LO.insert({DISNCBasis::CNS, 2 * M_PI * Id});
    C2LO.insert({DISNCBasis::CS,  2 * M_PI * Id});
    C2LO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> C2NLO;
    C2NLO.insert({DISNCBasis::CNS, Zero});
    C2NLO.insert({DISNCBasis::CS,  Zero});
    C2NLO.insert({DISNCBasis::CG,  Zero});

    // Vector of distributions to skip
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

    // Define object of the structure containing the DglapObjects
    const auto F2Obj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
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
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2LO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2NLO}});
        }
      return FObj;
    };
    t.stop();

    return F2Obj;
  }

  //_____________________________________________________________________________
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeImCFF2NCObjectsZMTimeLike(Grid                const& g,
                                                                                                                         std::vector<double> const& Thresholds,
                                                                                                                         double              const& xi,
                                                                                                                         double              const& IntEps)
  {
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = InitializeGpdObjects(g, Thresholds, xi);

    report("Initializing StructureFunctionObjects for the imaginary part of Compton form factor F2 Zero Mass at xi = " + std::to_string(xi) + " (Q2 < 0)... ");
    Timer t;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps, true};
    const Operator Zero{g, Null{},     IntEps, true};

    // LO
    std::map<int, Operator> C2LO;
    C2LO.insert({DISNCBasis::CNS, - 2 * M_PI * Id});
    C2LO.insert({DISNCBasis::CS,  - 2 * M_PI * Id});
    C2LO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> C2NLO;
    C2NLO.insert({DISNCBasis::CNS, Zero});
    C2NLO.insert({DISNCBasis::CS,  Zero});
    C2NLO.insert({DISNCBasis::CG,  Zero});

    // Vector of distributions to skip
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

    // Define object of the structure containing the DglapObjects
    const auto F2Obj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
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
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2LO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2NLO}});
        }
      return FObj;
    };
    t.stop();

    return F2Obj;
  }

  //_____________________________________________________________________________
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeReCFF2NCObjectsZMSpaceLike(Grid                const& g,
                                                                                                                          std::vector<double> const& Thresholds,
                                                                                                                          double              const& xi,
                                                                                                                          double              const& IntEps)
  {
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = InitializeGpdObjects(g, Thresholds, xi);

    report("Initializing StructureFunctionObjects for the real part of Compton form factor F2 Zero Mass at xi = " + std::to_string(xi) + " (Q2 > 0)... ");
    Timer t;

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps, true};
    const Operator O20ns{g, ReCFF20ns{}, IntEps, true};

    // LO
    std::map<int, Operator> C2LO;
    C2LO.insert({DISNCBasis::CNS, O20ns});
    C2LO.insert({DISNCBasis::CS,  O20ns});
    C2LO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> C2NLO;
    C2NLO.insert({DISNCBasis::CNS, Zero});
    C2NLO.insert({DISNCBasis::CS,  Zero});
    C2NLO.insert({DISNCBasis::CG,  Zero});

    // Vector of distributions to skip
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

    // Define object of the structure containing the DglapObjects
    const auto F2Obj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
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
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2LO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2NLO}});
        }
      return FObj;
    };
    t.stop();

    return F2Obj;
  }

  //_____________________________________________________________________________
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeReCFF2NCObjectsZMTimeLike(Grid                const& g,
                                                                                                                         std::vector<double> const& Thresholds,
                                                                                                                         double              const& xi,
                                                                                                                         double              const& IntEps)
  {
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = InitializeGpdObjects(g, Thresholds, xi);

    report("Initializing StructureFunctionObjects for the real part of Compton form factor F2 Zero Mass at xi = " + std::to_string(xi) + " (Q2 < 0)... ");
    Timer t;

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps, true};
    const Operator O20ns{g, ReCFF20ns{}, IntEps, true};

    // LO
    std::map<int, Operator> C2LO;
    C2LO.insert({DISNCBasis::CNS, O20ns});
    C2LO.insert({DISNCBasis::CS,  O20ns});
    C2LO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> C2NLO;
    C2NLO.insert({DISNCBasis::CNS, Zero});
    C2NLO.insert({DISNCBasis::CS,  Zero});
    C2NLO.insert({DISNCBasis::CG,  Zero});

    // Vector of distributions to skip
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

    // Define object of the structure containing the DglapObjects
    const auto F2Obj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
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
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2LO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2NLO}});
        }
      return FObj;
    };
    t.stop();

    return F2Obj;
  }

  //_____________________________________________________________________________
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeImCFFLNCObjectsZMSpaceLike(Grid                const& g,
                                                                                                                          std::vector<double> const& Thresholds,
                                                                                                                          double              const& xi,
                                                                                                                          double              const& IntEps)
  {
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = InitializeGpdObjects(g, Thresholds, xi);

    report("Initializing StructureFunctionObjects for the imaginary part of Compton form factor FL Zero Mass at xi = " + std::to_string(xi) + " (Q2 > 0)... ");
    Timer t;

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps, true};

    // LO
    std::map<int, Operator> CLLO;
    CLLO.insert({DISNCBasis::CNS, Zero});
    CLLO.insert({DISNCBasis::CS,  Zero});
    CLLO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> CLNLO;
    CLNLO.insert({DISNCBasis::CNS, Zero});
    CLNLO.insert({DISNCBasis::CS,  Zero});
    CLNLO.insert({DISNCBasis::CG,  Zero});

    // Vector of distributions to skip
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

    // Define object of the structure containing the DglapObjects
    const auto FLObj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
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
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLLO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLNLO}});
        }
      return FObj;
    };
    t.stop();

    return FLObj;
  }

  //_____________________________________________________________________________
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeImCFFLNCObjectsZMTimeLike(Grid                const& g,
                                                                                                                         std::vector<double> const& Thresholds,
                                                                                                                         double              const& xi,
                                                                                                                         double              const& IntEps)
  {
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = InitializeGpdObjects(g, Thresholds, xi);

    report("Initializing StructureFunctionObjects for the imaginary part of Compton form factor FL Zero Mass at xi = " + std::to_string(xi) + " (Q2 < 0)... ");
    Timer t;

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps, true};

    // LO
    std::map<int, Operator> CLLO;
    CLLO.insert({DISNCBasis::CNS, Zero});
    CLLO.insert({DISNCBasis::CS,  Zero});
    CLLO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> CLNLO;
    CLNLO.insert({DISNCBasis::CNS, Zero});
    CLNLO.insert({DISNCBasis::CS,  Zero});
    CLNLO.insert({DISNCBasis::CG,  Zero});

    // Vector of distributions to skip
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

    // Define object of the structure containing the DglapObjects
    const auto FLObj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
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
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLLO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLNLO}});
        }
      return FObj;
    };
    t.stop();

    return FLObj;
  }

  //_____________________________________________________________________________
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeReCFFLNCObjectsZMSpaceLike(Grid                const& g,
                                                                                                                          std::vector<double> const& Thresholds,
                                                                                                                          double              const& xi,
                                                                                                                          double              const& IntEps)
  {
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = InitializeGpdObjects(g, Thresholds, xi);

    report("Initializing StructureFunctionObjects for the real part of Compton form factor FL Zero Mass at xi = " + std::to_string(xi) + " (Q2 > 0)... ");
    Timer t;

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps, true};

    // LO
    std::map<int, Operator> CLLO;
    CLLO.insert({DISNCBasis::CNS, Zero});
    CLLO.insert({DISNCBasis::CS,  Zero});
    CLLO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> CLNLO;
    CLNLO.insert({DISNCBasis::CNS, Zero});
    CLNLO.insert({DISNCBasis::CS,  Zero});
    CLNLO.insert({DISNCBasis::CG,  Zero});

    // Vector of distributions to skip
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

    // Define object of the structure containing the DglapObjects
    const auto FLObj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
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
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLLO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLNLO}});
        }
      return FObj;
    };
    t.stop();

    return FLObj;
  }

  //_____________________________________________________________________________
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeReCFFLNCObjectsZMTimeLike(Grid                const& g,
                                                                                                                         std::vector<double> const& Thresholds,
                                                                                                                         double              const& xi,
                                                                                                                         double              const& IntEps)
  {
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = InitializeGpdObjects(g, Thresholds, xi);

    report("Initializing StructureFunctionObjects for the real part of Compton form factor FL Zero Mass at xi = " + std::to_string(xi) + " (Q2 < 0)... ");
    Timer t;

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps, true};

    // LO
    std::map<int, Operator> CLLO;
    CLLO.insert({DISNCBasis::CNS, Zero});
    CLLO.insert({DISNCBasis::CS,  Zero});
    CLLO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> CLNLO;
    CLNLO.insert({DISNCBasis::CNS, Zero});
    CLNLO.insert({DISNCBasis::CS,  Zero});
    CLNLO.insert({DISNCBasis::CG,  Zero});

    // Vector of distributions to skip
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

    // Define object of the structure containing the DglapObjects
    const auto FLObj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
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
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLLO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLNLO}});
        }
      return FObj;
    };
    t.stop();

    return FLObj;
  }
}
