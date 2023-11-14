//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/structurefunctionbuilder.h"
#include "apfel/timer.h"
#include "apfel/tools.h"
#include "apfel/constants.h"
#include "apfel/zeromasscoefficientfunctionsunp_sl.h"
#include "apfel/zeromasscoefficientfunctionspol_sl.h"
#include "apfel/massivecoefficientfunctionsunp_sl.h"
#include "apfel/massivezerocoefficientfunctionsunp_sl.h"
#include "apfel/zeromasscoefficientfunctionsunp_tl.h"
#include "apfel/tabulateobject.h"
#include "apfel/betaqcd.h"

#include <numeric>

namespace apfel
{
  //_____________________________________________________________________________
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF2NCObjectsZM(Grid                const& g,
                                                                                                             std::vector<double> const& Thresholds,
                                                                                                             double              const& IntEps)
  {
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = apfel::InitializeDglapObjectsQCD(g, Thresholds);

    report("Initializing StructureFunctionObjects for F2 NC Zero Mass... ");
    Timer t;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // LO
    std::map<int, Operator> C2LO;
    C2LO.insert({DISNCBasis::CNS, Id});
    C2LO.insert({DISNCBasis::CS,  Id});
    C2LO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> C2NLO;
    const Operator O21ns{g, C21ns{}, IntEps};
    const Operator O21g {g, C21g{},  IntEps};
    C2NLO.insert({DISNCBasis::CNS, O21ns});
    C2NLO.insert({DISNCBasis::CS,  O21ns});
    C2NLO.insert({DISNCBasis::CG,  O21g});

    // NNLO
    std::map<int, std::map<int, Operator>> C2NNLO;
    const Operator O22ps{g, C22ps{}, IntEps};
    const Operator O22g {g, C22g{},  IntEps};
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator O22nsp{g, C22nsp{nf}, IntEps};
        const Operator O22t = O22nsp + 6 * O22ps;
        std::map<int, Operator> C2NNLOnf;
        C2NNLOnf.insert({DISNCBasis::CNS, O22nsp});
        C2NNLOnf.insert({DISNCBasis::CS,  O22t});
        C2NNLOnf.insert({DISNCBasis::CG,  O22g});
        C2NNLO.insert({nf, C2NNLOnf});
      }

    // NNNLO
    std::map<int, std::map<int, Operator>> C2NNNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator O23ps{g, C23ps{nf}, IntEps};
        const Operator O23g {g, C23g{nf},  IntEps};
        const Operator O23nsp{g, C23nsp{nf}, IntEps};
        const Operator O23t = O23nsp + 6 * O23ps;
        std::map<int, Operator> C2NNNLOnf;
        C2NNNLOnf.insert({DISNCBasis::CNS, O23nsp});
        C2NNNLOnf.insert({DISNCBasis::CS,  O23t});
        C2NNNLOnf.insert({DISNCBasis::CG,  O23g});
        C2NNNLO.insert({nf, C2NNNLOnf});
      }

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
          FObj.ConvBasis.insert({k, (k == 0 ? DISNCBasis{EffCh} : DISNCBasis{k, EffCh[k-1]})});
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2LO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2NLO}});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2NNLO.at(nf)}});
          FObj.C3.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2NNNLO.at(nf)}});
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
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = apfel::InitializeDglapObjectsQCD(g, Thresholds);

    report("Initializing StructureFunctionObjects for FL NC Zero Mass... ");
    Timer t;

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps};

    // LO
    std::map<int, Operator> CLLO;
    CLLO.insert({DISNCBasis::CNS, Zero});
    CLLO.insert({DISNCBasis::CS,  Zero});
    CLLO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> CLNLO;
    const Operator OL1ns{g, CL1ns{}, IntEps};
    const Operator OL1g {g, CL1g{},  IntEps};
    CLNLO.insert({DISNCBasis::CNS, OL1ns});
    CLNLO.insert({DISNCBasis::CS,  OL1ns});
    CLNLO.insert({DISNCBasis::CG,  OL1g});

    // NNLO
    std::map<int, std::map<int, Operator>> CLNNLO;
    const Operator OL2ps{g, CL2ps{}, IntEps};
    const Operator OL2g {g, CL2g{},  IntEps};
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator OL2nsp{g, CL2nsp{nf}, IntEps};
        const Operator OL2t = OL2nsp + 6 * OL2ps;
        std::map<int, Operator> CLNNLOnf;
        CLNNLOnf.insert({DISNCBasis::CNS, OL2nsp});
        CLNNLOnf.insert({DISNCBasis::CS,  OL2t});
        CLNNLOnf.insert({DISNCBasis::CG,  OL2g});
        CLNNLO.insert({nf, CLNNLOnf});
      }

    // NNNLO
    std::map<int, std::map<int, Operator>> CLNNNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator OL3ps{g, CL3ps{nf}, IntEps};
        const Operator OL3g {g, CL3g{nf},  IntEps};
        const Operator OL3nsp{g, CL3nsp{nf}, IntEps};
        const Operator OL3t = OL3nsp + 6 * OL3ps;
        std::map<int, Operator> CLNNNLOnf;
        CLNNNLOnf.insert({DISNCBasis::CNS, OL3nsp});
        CLNNNLOnf.insert({DISNCBasis::CS,  OL3t});
        CLNNNLOnf.insert({DISNCBasis::CG,  OL3g});
        CLNNNLO.insert({nf, CLNNNLOnf});
      }

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
      FObj.skip = skip;
      FObj.nf   = nf;
      FObj.P    = PDFObj.at(nf);
      // Single structure function components.
      for (int k = 0; k <= 6; k++)
        {
          FObj.ConvBasis.insert({k, (k == 0 ? DISNCBasis{EffCh} : DISNCBasis{k, EffCh[k-1]})});
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLLO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLNLO}});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLNNLO.at(nf)}});
          FObj.C3.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLNNNLO.at(nf)}});
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
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = apfel::InitializeDglapObjectsQCD(g, Thresholds);

    report("Initializing StructureFunctionObjects for F3 NC Zero Mass... ");
    Timer t;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // LO
    std::map<int, Operator> C3LO;
    C3LO.insert({DISNCBasis::CNS, Id});
    C3LO.insert({DISNCBasis::CS,  Id});
    C3LO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> C3NLO;
    const Operator O31ns{g, C31ns{}, IntEps};
    C3NLO.insert({DISNCBasis::CNS, O31ns});
    C3NLO.insert({DISNCBasis::CS,  O31ns});
    C3NLO.insert({DISNCBasis::CG,  Zero});

    // NNLO
    std::map<int, std::map<int, Operator>> C3NNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator O32nsm{g, C32nsm{nf}, IntEps};
        std::map<int, Operator> C3NNLOnf;
        C3NNLOnf.insert({DISNCBasis::CNS, O32nsm});
        C3NNLOnf.insert({DISNCBasis::CS,  O32nsm});
        C3NNLOnf.insert({DISNCBasis::CG,  Zero});
        C3NNLO.insert({nf, C3NNLOnf});
      }

    // NNNLO
    std::map<int, std::map<int, Operator>> C3NNNLO;
    const Operator O33nsv{g, C33nsv{}, IntEps};
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator O33nsm{g, C33nsm{nf}, IntEps};
        const Operator O33t = O33nsm + 6 * O33nsv;
        std::map<int, Operator> C3NNNLOnf;
        C3NNNLOnf.insert({DISNCBasis::CNS, O33nsm});
        C3NNNLOnf.insert({DISNCBasis::CS,  O33t});
        C3NNNLOnf.insert({DISNCBasis::CG,  Zero});
        C3NNNLO.insert({nf, C3NNNLOnf});
      }

    // Vector of distributions to skip
    const std::vector<int> skip = {1, 3, 5, 7, 9, 11};

    // Define object of the structure containing the DglapObjects
    const auto F3Obj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
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
          FObj.ConvBasis.insert({k, (k == 0 ? DISNCBasis{EffCh} : DISNCBasis{k, EffCh[k-1]})});
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C3LO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C3NLO}});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C3NNLO.at(nf)}});
          FObj.C3.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C3NNNLO.at(nf)}});
        }
      return FObj;
    };
    t.stop();

    return F3Obj;
  }

  //_____________________________________________________________________________
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> Initializeg4NCObjectsZM(Grid                const& g,
                                                                                                             std::vector<double> const& Thresholds,
                                                                                                             double              const& IntEps)
  {
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = apfel::InitializeDglapObjectsQCDpol(g, Thresholds);

    report("Initializing StructureFunctionObjects for g4 NC Zero Mass... ");
    Timer t;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // LO
    std::map<int, Operator> G4LO;
    G4LO.insert({DISNCBasis::CNS, Id});
    G4LO.insert({DISNCBasis::CS,  Id});
    G4LO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> G4NLO;
    const Operator O41ns{g, G41ns{}, IntEps};
    G4NLO.insert({DISNCBasis::CNS, O41ns});
    G4NLO.insert({DISNCBasis::CS,  O41ns});
    G4NLO.insert({DISNCBasis::CG,  Zero});

    // NNLO
    // According to Eq. (19) of https://arxiv.org/pdf/2210.12014.pdf
    // the NNLO corrections to g4 coincide with the NNLO non-singlet
    // corrections to F2.
    std::map<int, std::map<int, Operator>> G4NNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator O42nsp{g, C22nsp{nf}, IntEps};
        std::map<int, Operator> G4NNLOnf;
        G4NNLOnf.insert({DISNCBasis::CNS, O42nsp});
        G4NNLOnf.insert({DISNCBasis::CS,  O42nsp});
        G4NNLOnf.insert({DISNCBasis::CG,  Zero});
        G4NNLO.insert({nf, G4NNLOnf});
      }

    // Vector of distributions to skip
    const std::vector<int> skip = {1, 3, 5, 7, 9, 11};

    // Define object of the structure containing the DglapObjects
    const auto g4Obj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
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
          FObj.ConvBasis.insert({k, (k == 0 ? DISNCBasis{EffCh} : DISNCBasis{k, EffCh[k-1]})});
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), G4LO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), G4NLO}});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), G4NNLO.at(nf)}});
        }
      return FObj;
    };
    t.stop();

    return g4Obj;
  }

  //_____________________________________________________________________________
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializegLNCObjectsZM(Grid                const& g,
                                                                                                             std::vector<double> const& Thresholds,
                                                                                                             double              const& IntEps)
  {
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = apfel::InitializeDglapObjectsQCDpol(g, Thresholds);

    report("Initializing StructureFunctionObjects for gL NC Zero Mass... ");
    Timer t;

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps};

    // LO
    std::map<int, Operator> GLLO;
    GLLO.insert({DISNCBasis::CNS, Zero});
    GLLO.insert({DISNCBasis::CS,  Zero});
    GLLO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> GLNLO;
    const Operator OL1ns{g, GL1ns{}, IntEps};
    GLNLO.insert({DISNCBasis::CNS, OL1ns});
    GLNLO.insert({DISNCBasis::CS,  OL1ns});
    GLNLO.insert({DISNCBasis::CG,  Zero});

    // NNLO
    // According to Eq. (19) of https://arxiv.org/pdf/2210.12014.pdf
    // the NNLO corrections to gL coincide with the NNLO non-singlet
    // corrections to FL.
    std::map<int, std::map<int, Operator>> GLNNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator OL2nsp{g, CL2nsp{nf}, IntEps};
        std::map<int, Operator> GLNNLOnf;
        GLNNLOnf.insert({DISNCBasis::CNS, OL2nsp});
        GLNNLOnf.insert({DISNCBasis::CS,  OL2nsp});
        GLNNLOnf.insert({DISNCBasis::CG,  Zero});
        GLNNLO.insert({nf, GLNNLOnf});
      }

    // Vector of distributions to skip
    const std::vector<int> skip = {1, 3, 5, 7, 9, 11};

    // Define object of the structure containing the DglapObjects
    const auto gLObj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
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
          FObj.ConvBasis.insert({k, (k == 0 ? DISNCBasis{EffCh} : DISNCBasis{k, EffCh[k-1]})});
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), GLLO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), GLNLO}});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), GLNNLO.at(nf)}});
        }
      return FObj;
    };
    t.stop();

    return gLObj;
  }

  //_____________________________________________________________________________
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> Initializeg1NCObjectsZM(Grid                const& g,
                                                                                                             std::vector<double> const& Thresholds,
                                                                                                             double              const& IntEps)
  {
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = apfel::InitializeDglapObjectsQCDpol(g, Thresholds);

    report("Initializing StructureFunctionObjects for g1 NC Zero Mass... ");
    Timer t;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // LO
    std::map<int, Operator> G1LO;
    G1LO.insert({DISNCBasis::CNS, Id});
    G1LO.insert({DISNCBasis::CS,  Id});
    G1LO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> G1NLO;
    const Operator O11ns{g, G11ns{}, IntEps};
    const Operator O11g {g, G11g{},  IntEps};
    G1NLO.insert({DISNCBasis::CNS, O11ns});
    G1NLO.insert({DISNCBasis::CS,  O11ns});
    G1NLO.insert({DISNCBasis::CG,  O11g});

    // NNLO
    std::map<int, std::map<int, Operator>> G1NNLO;
    const Operator O12ps{g, G12ps{}, IntEps};
    const Operator O12g {g, G12g{},  IntEps};
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator O12nsp{g, G12nsp{nf}, IntEps};
        const Operator O12t = O12nsp + 6 * O12ps;
        std::map<int, Operator> G1NNLOnf;
        G1NNLOnf.insert({DISNCBasis::CNS, O12nsp});
        G1NNLOnf.insert({DISNCBasis::CS,  O12t});
        G1NNLOnf.insert({DISNCBasis::CG,  O12g});
        G1NNLO.insert({nf, G1NNLOnf});
      }

    // Vector of distributions to skip
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

    // Define object of the structure containing the DglapObjects
    const auto g1Obj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
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
          FObj.ConvBasis.insert({k, (k == 0 ? DISNCBasis{EffCh} : DISNCBasis{k, EffCh[k-1]})});
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), G1LO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), G1NLO}});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), G1NNLO.at(nf)}});
        }
      return FObj;
    };
    t.stop();

    return g1Obj;
  }

  //_____________________________________________________________________________
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF2CCPlusObjectsZM(Grid                const& g,
                                                                                                                 std::vector<double> const& Thresholds,
                                                                                                                 double              const& IntEps)
  {
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = apfel::InitializeDglapObjectsQCD(g, Thresholds);

    report("Initializing StructureFunctionObjects for ( F2(nu) + F2(nubar) ) / 2 Zero Mass... ");
    Timer t;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // LO
    std::map<int, Operator> C2LO;
    C2LO.insert({DISCCBasis::CNS, Id});
    C2LO.insert({DISCCBasis::CS,  Id});
    C2LO.insert({DISCCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> C2NLO;
    const Operator O21ns{g, C21ns{}, IntEps};
    const Operator O21g {g, C21g{},  IntEps};
    C2NLO.insert({DISCCBasis::CNS, O21ns});
    C2NLO.insert({DISCCBasis::CS,  O21ns});
    C2NLO.insert({DISCCBasis::CG,  O21g});

    // NNLO
    std::map<int, std::map<int, Operator>> C2NNLO;
    const Operator O22ps{g, C22ps{}, IntEps};
    const Operator O22g {g, C22g{},  IntEps};
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator O22nsp{g, C22nsp{nf}, IntEps};
        const Operator O22t = O22nsp + 6 * O22ps;
        std::map<int, Operator> C2NNLOnf;
        C2NNLOnf.insert({DISCCBasis::CNS, O22nsp});
        C2NNLOnf.insert({DISCCBasis::CS,  O22t});
        C2NNLOnf.insert({DISCCBasis::CG,  O22g});
        C2NNLO.insert({nf, C2NNLOnf});
      }

    // NNNLO
    std::map<int, std::map<int, Operator>> C2NNNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator O23ps{g, C23ps{nf, false}, IntEps};
        const Operator O23g {g, C23g{nf, false},  IntEps};
        const Operator O23nsp{g, C23nsp{nf, false}, IntEps};
        const Operator O23t = O23nsp + 6 * O23ps;
        std::map<int, Operator> C2NNNLOnf;
        C2NNNLOnf.insert({DISCCBasis::CNS, O23nsp});
        C2NNNLOnf.insert({DISCCBasis::CS,  O23t});
        C2NNNLOnf.insert({DISCCBasis::CG,  O23g});
        C2NNNLO.insert({nf, C2NNNLOnf});
      }

    // Vector of distributions to skip
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

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
      FObj.nf   = nf;
      FObj.P    = PDFObj.at(nf);
      FObj.skip = skip;
      // Single structure function components.
      for (int k = 0; k <= 6; k++)
        {
          FObj.ConvBasis.insert({k, (k == 0 ? DISCCBasis{EffCh, false} : DISCCBasis{k, false, EffCh[k-1]})});
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2LO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2NLO}});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2NNLO.at(nf)}});
          FObj.C3.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2NNNLO.at(nf)}});
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
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = apfel::InitializeDglapObjectsQCD(g, Thresholds);

    report("Initializing StructureFunctionObjects for ( F2(nu) - F2(nubar) ) / 2 Zero Mass... ");
    Timer t;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // LO
    std::map<int, Operator> C2LO;
    C2LO.insert({DISCCBasis::CNS, Id});
    C2LO.insert({DISCCBasis::CS,  Zero});
    C2LO.insert({DISCCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> C2NLO;
    const Operator O21ns{g, C21ns{}, IntEps};
    C2NLO.insert({DISCCBasis::CNS, O21ns});
    C2NLO.insert({DISCCBasis::CS,  Zero});
    C2NLO.insert({DISCCBasis::CG,  Zero});

    // NNLO
    std::map<int, std::map<int, Operator>> C2NNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator O22nsm{g, C22nsm{nf}, IntEps};
        std::map<int, Operator> C2NNLOnf;
        C2NNLOnf.insert({DISCCBasis::CNS, O22nsm});
        C2NNLOnf.insert({DISCCBasis::CS,  Zero});
        C2NNLOnf.insert({DISCCBasis::CG,  Zero});
        C2NNLO.insert({nf, C2NNLOnf});
      }

    // NNNLO
    std::map<int, std::map<int, Operator>> C2NNNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator O23nsm{g, C23nsm{nf, false}, IntEps};
        std::map<int, Operator> C2NNNLOnf;
        C2NNNLOnf.insert({DISCCBasis::CNS, O23nsm});
        C2NNNLOnf.insert({DISCCBasis::CS,  Zero});
        C2NNNLOnf.insert({DISCCBasis::CG,  Zero});
        C2NNNLO.insert({nf, C2NNNLOnf});
      }

    // Vector of distributions to skip
    const std::vector<int> skip = {0, 1, 2, 3, 5, 7, 9, 11};

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
      FObj.nf   = nf;
      FObj.P    = PDFObj.at(nf);
      FObj.skip = skip;
      // Single structure function components.
      for (int k = 0; k <= 6; k++)
        {
          FObj.ConvBasis.insert({k, (k == 0 ? DISCCBasis{EffCh, false} : DISCCBasis{k, false, EffCh[k-1]})});
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2LO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2NLO}});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2NNLO.at(nf)}});
          FObj.C3.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2NNNLO.at(nf)}});
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
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = apfel::InitializeDglapObjectsQCD(g, Thresholds);

    report("Initializing StructureFunctionObjects for ( FL(nu) + FL(nubar) ) / 2 Zero Mass... ");
    Timer t;

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps};

    // LO
    std::map<int, Operator> CLLO;
    CLLO.insert({DISCCBasis::CNS, Zero});
    CLLO.insert({DISCCBasis::CS,  Zero});
    CLLO.insert({DISCCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> CLNLO;
    const Operator OL1ns{g, CL1ns{}, IntEps};
    const Operator OL1g {g, CL1g{},  IntEps};
    CLNLO.insert({DISCCBasis::CNS, OL1ns});
    CLNLO.insert({DISCCBasis::CS,  OL1ns});
    CLNLO.insert({DISCCBasis::CG,  OL1g});

    // NNLO
    std::map<int, std::map<int, Operator>> CLNNLO;
    const Operator OL2ps{g, CL2ps{}, IntEps};
    const Operator OL2g {g, CL2g{},  IntEps};
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator OL2nsp{g, CL2nsp{nf}, IntEps};
        const Operator OL2t = OL2nsp + 6 * OL2ps;
        std::map<int, Operator> CLNNLOnf;
        CLNNLOnf.insert({DISCCBasis::CNS, OL2nsp});
        CLNNLOnf.insert({DISCCBasis::CS,  OL2t});
        CLNNLOnf.insert({DISCCBasis::CG,  OL2g});
        CLNNLO.insert({nf, CLNNLOnf});
      }

    // NNNLO
    std::map<int, std::map<int, Operator>> CLNNNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator OL3ps{g, CL3ps{nf, false}, IntEps};
        const Operator OL3g {g, CL3g{nf, false},  IntEps};
        const Operator OL3nsp{g, CL3nsp{nf, false}, IntEps};
        const Operator OL3t = OL3nsp + 6 * OL3ps;
        std::map<int, Operator> CLNNNLOnf;
        CLNNNLOnf.insert({DISCCBasis::CNS, OL3nsp});
        CLNNNLOnf.insert({DISCCBasis::CS,  OL3t});
        CLNNNLOnf.insert({DISCCBasis::CG,  OL3g});
        CLNNNLO.insert({nf, CLNNNLOnf});
      }

    // Vector of distributions to skip
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

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
      FObj.nf   = nf;
      FObj.P    = PDFObj.at(nf);
      FObj.skip = skip;
      // Single structure function components.
      for (int k = 0; k <= 6; k++)
        {
          FObj.ConvBasis.insert({k, (k == 0 ? DISCCBasis{EffCh, false} : DISCCBasis{k, false, EffCh[k-1]})});
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLLO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLNLO}});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLNNLO.at(nf)}});
          FObj.C3.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLNNNLO.at(nf)}});
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
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = apfel::InitializeDglapObjectsQCD(g, Thresholds);

    report("Initializing StructureFunctionObjects for ( FL(nu) - FL(nubar) ) / 2 Zero Mass... ");
    Timer t;

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps};

    // LO
    std::map<int, Operator> CLLO;
    CLLO.insert({DISCCBasis::CNS, Zero});
    CLLO.insert({DISCCBasis::CS,  Zero});
    CLLO.insert({DISCCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> CLNLO;
    const Operator OL1ns{g, CL1ns{}, IntEps};
    CLNLO.insert({DISCCBasis::CNS, OL1ns});
    CLNLO.insert({DISCCBasis::CS,  Zero});
    CLNLO.insert({DISCCBasis::CG,  Zero});

    // NNLO
    std::map<int, std::map<int, Operator>> CLNNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator OL2nsm{g, CL2nsm{nf}, IntEps};
        std::map<int, Operator> CLNNLOnf;
        CLNNLOnf.insert({DISCCBasis::CNS, OL2nsm});
        CLNNLOnf.insert({DISCCBasis::CS,  Zero});
        CLNNLOnf.insert({DISCCBasis::CG,  Zero});
        CLNNLO.insert({nf, CLNNLOnf});
      }

    // NNNLO
    std::map<int, std::map<int, Operator>> CLNNNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator OL3nsm{g, CL3nsm{nf, false}, IntEps};
        std::map<int, Operator> CLNNNLOnf;
        CLNNNLOnf.insert({DISCCBasis::CNS, OL3nsm});
        CLNNNLOnf.insert({DISCCBasis::CS,  Zero});
        CLNNNLOnf.insert({DISCCBasis::CG,  Zero});
        CLNNNLO.insert({nf, CLNNNLOnf});
      }

    // Vector of distributions to skip
    const std::vector<int> skip = {0, 1, 2, 3, 5, 7, 9, 11};

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
      FObj.nf   = nf;
      FObj.P    = PDFObj.at(nf);
      FObj.skip = skip;
      // Single structure function components.
      for (int k = 0; k <= 6; k++)
        {
          FObj.ConvBasis.insert({k, (k == 0 ? DISCCBasis{EffCh, false} : DISCCBasis{k, false, EffCh[k-1]})});
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLLO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLNLO}});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLNNLO.at(nf)}});
          FObj.C3.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLNNNLO.at(nf)}});
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
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = apfel::InitializeDglapObjectsQCD(g, Thresholds);

    report("Initializing StructureFunctionObjects for ( F3(nu) + F3(nubar) ) / 2 Zero Mass... ");
    Timer t;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // LO
    std::map<int, Operator> C3LO;
    C3LO.insert({DISCCBasis::CNS, Id});
    C3LO.insert({DISCCBasis::CS,  Id});
    C3LO.insert({DISCCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> C3NLO;
    const Operator O31ns{g, C31ns{}, IntEps};
    C3NLO.insert({DISCCBasis::CNS, O31ns});
    C3NLO.insert({DISCCBasis::CS,  O31ns});
    C3NLO.insert({DISCCBasis::CG,  Zero});

    // NNLO
    std::map<int, std::map<int, Operator>> C3NNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator O32nsp{g, C32nsp{nf}, IntEps};
        std::map<int, Operator> C3NNLOnf;
        C3NNLOnf.insert({DISCCBasis::CNS, O32nsp});
        C3NNLOnf.insert({DISCCBasis::CS,  O32nsp});
        C3NNLOnf.insert({DISCCBasis::CG,  Zero});
        C3NNLO.insert({nf, C3NNLOnf});
      }

    // NNNLO
    std::map<int, std::map<int, Operator>> C3NNNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator O33nsp{g, C33nsp{nf}, IntEps};
        std::map<int, Operator> C3NNNLOnf;
        C3NNNLOnf.insert({DISCCBasis::CNS, O33nsp});
        C3NNNLOnf.insert({DISCCBasis::CS,  O33nsp});
        C3NNNLOnf.insert({DISCCBasis::CG,  Zero});
        C3NNNLO.insert({nf, C3NNNLOnf});
      }

    // Vector of distributions to skip
    const std::vector<int> skip = {0, 1, 2, 4, 6, 8, 10, 12};

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
      FObj.nf   = nf;
      FObj.P    = PDFObj.at(nf);
      FObj.skip = skip;
      // Single structure function components.
      for (int k = 0; k <= 6; k++)
        {
          FObj.ConvBasis.insert({k, (k == 0 ? DISCCBasis{EffCh, true} : DISCCBasis{k, true, EffCh[k-1]})});
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C3LO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C3NLO}});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C3NNLO.at(nf)}});
          FObj.C3.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C3NNNLO.at(nf)}});
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
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = apfel::InitializeDglapObjectsQCD(g, Thresholds);

    report("Initializing StructureFunctionObjects for ( F3(nu) - F3(nubar) ) / 2 Zero Mass... ");
    Timer t;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // LO
    std::map<int, Operator> C3LO;
    C3LO.insert({DISCCBasis::CNS, Id});
    C3LO.insert({DISCCBasis::CS,  Id});
    C3LO.insert({DISCCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> C3NLO;
    const Operator O31ns{g, C31ns{}, IntEps};
    C3NLO.insert({DISCCBasis::CNS, O31ns});
    C3NLO.insert({DISCCBasis::CS,  O31ns});
    C3NLO.insert({DISCCBasis::CG,  Zero});

    // NNLO
    std::map<int, std::map<int, Operator>> C3NNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator O32nsm{g, C32nsm{nf}, IntEps};
        std::map<int, Operator> C3NNLOnf;
        C3NNLOnf.insert({DISCCBasis::CNS, O32nsm});
        C3NNLOnf.insert({DISCCBasis::CS,  O32nsm});
        C3NNLOnf.insert({DISCCBasis::CG,  Zero});
        C3NNLO.insert({nf, C3NNLOnf});
      }

    // NNNLO
    std::map<int, std::map<int, Operator>> C3NNNLO;
    const Operator O33nsv{g, C33nsv{}, IntEps};
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator O33nsm{g, C33nsm{nf}, IntEps};
        const Operator O33t = O33nsm + 6 * O33nsv;
        std::map<int, Operator> C3NNNLOnf;
        C3NNNLOnf.insert({DISCCBasis::CNS, O33nsm});
        C3NNNLOnf.insert({DISCCBasis::CS,  O33t});
        C3NNNLOnf.insert({DISCCBasis::CG,  Zero});
        C3NNNLO.insert({nf, C3NNNLOnf});
      }

    // Vector of distributions to skip
    const std::vector<int> skip = {0, 1, 3, 5, 7, 9, 11};

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
      FObj.nf   = nf;
      FObj.P    = PDFObj.at(nf);
      FObj.skip = skip;
      // Single structure function components.
      for (int k = 0; k <= 6; k++)
        {
          FObj.ConvBasis.insert({k, (k == 0 ? DISCCBasis{EffCh, true} : DISCCBasis{k, true, EffCh[k-1]})});
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C3LO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C3NLO}});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C3NNLO.at(nf)}});
          FObj.C3.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C3NNNLO.at(nf)}});
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
    // Make sure that the vector of masses contains all the 6 masses.
    if (Masses.size() != 6)
      throw std::runtime_error(error("InitializeF2NCObjectsMassive", "The vector of masses has to contain exactly 6 ordered masses."));

    // Determine number of active flavours
    int actnf = 0;
    for (auto const& m : Masses)
      if (m < eps8)
        actnf++;

    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = apfel::InitializeDglapObjectsQCD(g, std::vector<double>(actnf, 0.));

    report("Initializing StructureFunctionObjects for F2 NC Massive with " + std::to_string(actnf) + " active flavours... ");
    Timer t;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // Zero Mass coefficient functions
    // LO
    std::map<int, Operator> C2LO;
    C2LO.insert({DISNCBasis::CNS, Id});
    C2LO.insert({DISNCBasis::CS,  Id});
    C2LO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> C2NLO;
    const Operator O21ns{g, C21ns{}, IntEps};
    const Operator O21g {g, C21g{},  IntEps};
    C2NLO.insert({DISNCBasis::CNS, O21ns});
    C2NLO.insert({DISNCBasis::CS,  O21ns});
    C2NLO.insert({DISNCBasis::CG,  O21g});

    // NNLO
    const Operator O22ps {g, C22ps{},       IntEps};
    const Operator O22g  {g, C22g{},        IntEps};
    const Operator O22nsp{g, C22nsp{actnf}, IntEps};
    const Operator O22t = O22nsp + 6 * O22ps;

    // Massive coefficient functions
    // Null set of operator needed for the LO coefficient functions.
    std::map<int, Operator> Not;
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
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

    // Define object of the structure containing the DglapObjects
    const auto F2Obj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
    {
      // Fill in structure function object
      StructureFunctionObjects FObj;
      FObj.nf   = actnf;
      FObj.P    = PDFObj.at(actnf);
      FObj.skip = skip;

      // Start with the heavy quark structure function components.
      const double Q2  = Q * Q;
      Operator lNNLOns = O22nsp;
      for (int k = actnf + 1; k <= 6; k++)
        {
          // Compute value of xi.
          const double M2 = Masses[k-1] * Masses[k-1];
          const double xi = Q2 / M2;

          // Boolean that defines whether xi is withing the allowed
          // interpolation range
          const bool xiin = xi >= ximin && xi <= ximax;

          // Convolution Basis
          FObj.ConvBasis.insert({k, DISNCBasis{k, Ch[k-1]}});

          // Set LO non-singlet and singlet coefficient
          // functions to zero (gluon is already zero).
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), Not}});

          // Now insert NLO
          std::map<int, Operator> NLO;
          NLO.insert({DISNCBasis::CNS, Zero});
          NLO.insert({DISNCBasis::CS,  Zero});
          NLO.insert({DISNCBasis::CG,  (xiin ? TabO21g.Evaluate(xi) : Zero)});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), NLO}});

          // Now insert NNLO
          std::map<int, Operator> NNLO;
          NNLO.insert({DISNCBasis::CNS, Zero});
          NNLO.insert({DISNCBasis::CS,  (xiin ? TabO22s.Evaluate(xi) : Zero)});
          NNLO.insert({DISNCBasis::CG,  (xiin ? TabO22g.Evaluate(xi) : Zero)});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), NNLO}});

          // Include the massive bit to the non-singlet coefficient
          // function needed for the light structure functions.
          if (xiin)
            lNNLOns += TabO22ns.Evaluate(xi);
        }

      // Now fill in the light components. To be updated in the
      // loop over the heavy component.
      std::map<int, Operator> C2NNLO;
      C2NNLO.insert({DISNCBasis::CNS, lNNLOns});
      C2NNLO.insert({DISNCBasis::CS,  O22t});
      C2NNLO.insert({DISNCBasis::CG,  O22g});
      for (int k = 1; k <= actnf; k++)
        {
          // Convolution Basis
          FObj.ConvBasis.insert({k, DISNCBasis{k, Ch[k-1]}});
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2LO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2NLO}});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2NNLO}});
        }

      // Total structure function. This simple construction is based
      // on the fact that up to O(as^2) there are no non-singlet
      // contributions to the massive components.
      std::vector<double> ChL(6, 0.);
      std::map<int, Operator> C2NLOtot  = C2NLO;
      std::map<int, Operator> C2NNLOtot = C2NNLO;
      for (int k = 0; k < actnf; k++)
        ChL[k] += Ch[k];
      const double SumChL = accumulate(ChL.begin(), ChL.end(), 0.);
      for (int k = actnf + 1; k <= 6; k++)
        {
          const double wch = Ch[k-1] / SumChL;
          C2NLOtot.at(DISNCBasis::CG)  += wch * FObj.C1.at(k).GetObjects().at(DISNCBasis::CG);
          C2NNLOtot.at(DISNCBasis::CS) += wch * FObj.C2.at(k).GetObjects().at(DISNCBasis::CS);
          C2NNLOtot.at(DISNCBasis::CG) += wch * FObj.C2.at(k).GetObjects().at(DISNCBasis::CG);
        }

      FObj.ConvBasis.insert({0, DISNCBasis{ChL}});
      FObj.C0.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C2LO}});
      FObj.C1.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C2NLOtot}});
      FObj.C2.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C2NNLOtot}});
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
    // Make sure that the vector of masses contains all the 6 masses.
    if (Masses.size() != 6)
      throw std::runtime_error(error("InitializeFLNCObjectsMassive", "The vector of masses has to contain exactly 6 ordered masses."));

    // Determine number of active flavours
    int actnf = 0;
    for (auto const& m : Masses)
      if (m < eps8)
        actnf++;

    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = apfel::InitializeDglapObjectsQCD(g, std::vector<double>(actnf, 0.));

    report("Initializing StructureFunctionObjects for FL NC Massive with " + std::to_string(actnf) + " active flavours... ");
    Timer t;

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps};

    // Zero Mass coefficient functions
    // LO
    std::map<int, Operator> CLLO;
    CLLO.insert({DISNCBasis::CNS, Zero});
    CLLO.insert({DISNCBasis::CS,  Zero});
    CLLO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> CLNLO;
    const Operator OL1ns{g, CL1ns{}, IntEps};
    const Operator OL1g {g, CL1g{},  IntEps};
    CLNLO.insert({DISNCBasis::CNS, OL1ns});
    CLNLO.insert({DISNCBasis::CS,  OL1ns});
    CLNLO.insert({DISNCBasis::CG,  OL1g});

    // NNLO
    const Operator OL2ps {g, CL2ps{},       IntEps};
    const Operator OL2g  {g, CL2g{},        IntEps};
    const Operator OL2nsp{g, CL2nsp{actnf}, IntEps};
    const Operator OL2t = OL2nsp + 6 * OL2ps;

    // Massive coefficient functions
    // Null set of operator needed for the LO coefficient functions.
    std::map<int, Operator> Not;
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
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

    // Define object of the structure containing the DglapObjects
    const auto FLObj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
    {
      // Fill in structure function object
      StructureFunctionObjects FObj;
      FObj.nf   = actnf;
      FObj.P    = PDFObj.at(actnf);
      FObj.skip = skip;

      // Start with the heavy quark structure function components.
      const double Q2  = Q * Q;
      Operator lNNLOns = OL2nsp;
      for (int k = actnf + 1; k <= 6; k++)
        {
          // Compute value of xi.
          const double M2 = Masses[k-1] * Masses[k-1];
          const double xi = Q2 / M2;

          // Boolean that defines whether xi is withing the allowed
          // interpolation range
          const bool xiin = xi >= ximin && xi <= ximax;

          // Convolution Basis
          FObj.ConvBasis.insert({k, DISNCBasis{k, Ch[k-1]}});

          // Set LO non-singlet and singlet coefficient
          // functions to zero (gluon is already zero).
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), Not}});

          // Now insert NLO
          std::map<int, Operator> NLO;
          NLO.insert({DISNCBasis::CNS, Zero});
          NLO.insert({DISNCBasis::CS,  Zero});
          NLO.insert({DISNCBasis::CG,  (xiin ? TabOL1g.Evaluate(xi) : Zero)});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), NLO}});

          // Now insert NNLO
          std::map<int, Operator> NNLO;
          NNLO.insert({DISNCBasis::CNS, Zero});
          NNLO.insert({DISNCBasis::CS,  (xiin ? TabOL2s.Evaluate(xi) : Zero)});
          NNLO.insert({DISNCBasis::CG,  (xiin ? TabOL2g.Evaluate(xi) : Zero)});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), NNLO}});

          // Include the massive bit to the non-singlet coefficient
          // function needed for the light structure functions.
          if (xiin)
            lNNLOns += TabOL2ns.Evaluate(xi);
        }

      // Now fill in the light components. To be updated in the
      // loop over the heavy component.
      std::map<int, Operator> CLNNLO;
      CLNNLO.insert({DISNCBasis::CNS, lNNLOns});
      CLNNLO.insert({DISNCBasis::CS,  OL2t});
      CLNNLO.insert({DISNCBasis::CG,  OL2g});
      for (int k = 1; k <= actnf; k++)
        {
          // Convolution Basis
          FObj.ConvBasis.insert({k, DISNCBasis{k, Ch[k-1]}});
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLLO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLNLO}});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLNNLO}});
        }

      // Total structure function. This simple construction is based
      // on the fact that up to O(as^2) there are no non-singlet
      // contributions to the massive components.
      std::vector<double> ChL(6, 0.);
      std::map<int, Operator> CLNLOtot  = CLNLO;
      std::map<int, Operator> CLNNLOtot = CLNNLO;
      for (int k = 0; k < actnf; k++)
        ChL[k] += Ch[k];
      const double SumChL = accumulate(ChL.begin(), ChL.end(), 0.);
      for (int k = actnf + 1; k <= 6; k++)
        {
          const double wch = Ch[k-1] / SumChL;
          CLNLOtot.at(DISNCBasis::CG)  += wch * FObj.C1.at(k).GetObjects().at(DISNCBasis::CG);
          CLNNLOtot.at(DISNCBasis::CS) += wch * FObj.C2.at(k).GetObjects().at(DISNCBasis::CS);
          CLNNLOtot.at(DISNCBasis::CG) += wch * FObj.C2.at(k).GetObjects().at(DISNCBasis::CG);
        }

      FObj.ConvBasis.insert({0, DISNCBasis{ChL}});
      FObj.C0.insert({0, Set<Operator>{FObj.ConvBasis.at(0), CLLO}});
      FObj.C1.insert({0, Set<Operator>{FObj.ConvBasis.at(0), CLNLOtot}});
      FObj.C2.insert({0, Set<Operator>{FObj.ConvBasis.at(0), CLNNLOtot}});
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
    // Make sure that the vector of masses contains all the 6 masses.
    if (Masses.size() != 6)
      throw std::runtime_error(error("InitializeF2NCObjectsMassiveZero", "The vector of masses has to contain exactly 6 ordered masses."));

    // Determine number of active flavours
    int actnf = 0;
    for (auto const& m : Masses)
      if (m < eps8)
        actnf++;

    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = apfel::InitializeDglapObjectsQCD(g, std::vector<double>(actnf, 0.));

    report("Initializing StructureFunctionObjects for F2 NC Massive Zero with " + std::to_string(actnf) + " active flavours... ");
    Timer t;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // Zero Mass coefficient functions
    // LO
    std::map<int, Operator> C2LO;
    C2LO.insert({DISNCBasis::CNS, Id});
    C2LO.insert({DISNCBasis::CS,  Id});
    C2LO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> C2NLO;
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
    std::map<int, Operator> Not;
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
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

    // Define object of the structure containing the DglapObjects
    const auto F2Obj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
    {
      // Fill in structure function object
      StructureFunctionObjects FObj;
      FObj.nf   = actnf;
      FObj.P    = PDFObj.at(actnf);
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
          FObj.ConvBasis.insert({k, DISNCBasis{k, Ch[k-1]}});

          // Set LO non-singlet and singlet coefficient
          // functions to zero (gluon is already zero).
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), Not}});

          // Now insert NLO
          std::map<int, Operator> NLO;
          NLO.insert({DISNCBasis::CNS, Zero});
          NLO.insert({DISNCBasis::CS,  Zero});
          NLO.insert({DISNCBasis::CG,  TabO21g.Evaluate(xi)});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), NLO}});

          // Now insert NNLO
          std::map<int, Operator> NNLO;
          NNLO.insert({DISNCBasis::CNS, Zero});
          NNLO.insert({DISNCBasis::CS,  TabO22s.Evaluate(xi)});
          NNLO.insert({DISNCBasis::CG,  TabO22g.Evaluate(xi)});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), NNLO}});

          // Include the massive bit to the non-singlet coefficient
          // function needed for the light structure functions.
          for (int i = 0; i <= actnf; i++)
            lNNLOns += TabO22ns.Evaluate(xi);
        }

      // Now fill in the light components. To be updated in the
      // loop over the heavy component.
      std::map<int, Operator> C2NNLO;
      C2NNLO.insert({DISNCBasis::CNS, lNNLOns});
      C2NNLO.insert({DISNCBasis::CS,  O22t});
      C2NNLO.insert({DISNCBasis::CG,  O22g});
      for (int k = 1; k <= actnf; k++)
        {
          // Convolution Basis
          FObj.ConvBasis.insert({k, DISNCBasis{k, Ch[k-1]}});
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2LO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2NLO}});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2NNLO}});
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
    // Make sure that the vector of masses contains all the 6 masses.
    if (Masses.size() != 6)
      throw std::runtime_error(error("InitializeFLNCObjectsMassiveZero", "The vector of masses has to contain exactly 6 ordered masses."));

    // Determine number of active flavours
    int actnf = 0;
    for (auto const& m : Masses)
      if (m < eps8)
        actnf++;

    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = apfel::InitializeDglapObjectsQCD(g, std::vector<double>(actnf, 0.));

    report("Initializing StructureFunctionObjects for FL NC Massive Zero with " + std::to_string(actnf) + " active flavours... ");
    Timer t;

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps};

    // Zero Mass coefficient functions
    // LO
    std::map<int, Operator> CLLO;
    CLLO.insert({DISNCBasis::CNS, Zero});
    CLLO.insert({DISNCBasis::CS,  Zero});
    CLLO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> CLNLO;
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
    std::map<int, Operator> Not;
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
      FObj.nf   = actnf;
      FObj.P    = PDFObj.at(actnf);
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
          FObj.ConvBasis.insert({k, DISNCBasis{k, Ch[k-1]}});

          // Set LO non-singlet and singlet coefficient
          // functions to zero (gluon is already zero).
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), Not}});

          // Now insert NLO
          std::map<int, Operator> NLO;
          NLO.insert({DISNCBasis::CNS, Zero});
          NLO.insert({DISNCBasis::CS,  Zero});
          NLO.insert({DISNCBasis::CG,  Om0L1g});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), NLO}});

          // Now insert NNLO
          std::map<int, Operator> NNLO;
          NNLO.insert({DISNCBasis::CNS, Zero});
          NNLO.insert({DISNCBasis::CS,  TabOL2s.Evaluate(xi)});
          NNLO.insert({DISNCBasis::CG,  TabOL2g.Evaluate(xi)});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), NNLO}});

          // Include the massive bit to the non-singlet coefficient
          // function needed for the light structure functions.
          for (int i = 0; i <= actnf; i++)
            lNNLOns += TabOL2ns.Evaluate(xi);
        }

      // Now fill in the light components. To be updated in the
      // loop over the heavy component.
      std::map<int, Operator> CLNNLO;
      CLNNLO.insert({DISNCBasis::CNS, lNNLOns});
      CLNNLO.insert({DISNCBasis::CS,  OL2t});
      CLNNLO.insert({DISNCBasis::CG,  OL2g});
      for (int k = 1; k <= actnf; k++)
        {
          // Convolution Basis
          FObj.ConvBasis.insert({k, DISNCBasis{k, Ch[k-1]}});
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLLO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLNLO}});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLNNLO}});
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
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF2NCObjectsZMT(Grid                const& g,
                                                                                                              std::vector<double> const& Thresholds,
                                                                                                              double              const& IntEps)
  {
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = apfel::InitializeDglapObjectsQCDT(g, Thresholds);

    report("Initializing StructureFunctionObjects for F2 Zero Mass for SIA... ");
    Timer t;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // LO
    std::map<int, Operator> C2LO;
    C2LO.insert({DISNCBasis::CNS, Id});
    C2LO.insert({DISNCBasis::CS,  Id});
    C2LO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> C2NLO;
    const Operator O21ns{g, C21Tns{}, IntEps};
    const Operator O21g {g, C21Tg{},  IntEps};
    C2NLO.insert({DISNCBasis::CNS, O21ns});
    C2NLO.insert({DISNCBasis::CS,  O21ns});
    C2NLO.insert({DISNCBasis::CG,  O21g});

    // NNLO
    std::map<int, std::map<int, Operator>> C2NNLO;
    const Operator O22ps{g, C22Tps{}, IntEps};
    const Operator O22g {g, C22Tg{},  IntEps};
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator O22nsp{g, C22Tnsp{nf}, IntEps};
        const Operator O22t = O22nsp + 6 * O22ps;
        std::map<int, Operator> C2NNLOnf;
        C2NNLOnf.insert({DISNCBasis::CNS, O22nsp});
        C2NNLOnf.insert({DISNCBasis::CS,  O22t});
        C2NNLOnf.insert({DISNCBasis::CG,  O22g});
        C2NNLO.insert({nf, C2NNLOnf});
      }

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
          FObj.ConvBasis.insert({k, (k == 0 ? DISNCBasis{EffCh} : DISNCBasis{k, EffCh[k-1]})});
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2LO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2NLO}});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C2NNLO.at(nf)}});
        }
      return FObj;
    };
    t.stop();

    return F2Obj;
  }

  //_____________________________________________________________________________
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeFLNCObjectsZMT(Grid                const& g,
                                                                                                              std::vector<double> const& Thresholds,
                                                                                                              double              const& IntEps)
  {
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = apfel::InitializeDglapObjectsQCDT(g, Thresholds);

    report("Initializing StructureFunctionObjects for FL Zero Mass for SIA... ");
    Timer t;

    // ===============================================================
    const Operator Zero{g, Null{}, IntEps};

    // LO
    std::map<int, Operator> CLLO;
    CLLO.insert({DISNCBasis::CNS, Zero});
    CLLO.insert({DISNCBasis::CS,  Zero});
    CLLO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> CLNLO;
    const Operator OL1ns{g, CL1Tns{}, IntEps};
    const Operator OL1g {g, CL1Tg{},  IntEps};
    CLNLO.insert({DISNCBasis::CNS, OL1ns});
    CLNLO.insert({DISNCBasis::CS,  OL1ns});
    CLNLO.insert({DISNCBasis::CG,  OL1g});

    // NNLO (for nf from 1 to 6)
    std::map<int, std::map<int, Operator>> CLNNLO;
    const Operator OL2ps{g, CL2Tps{}, IntEps};
    const Operator OL2g {g, CL2Tg{},  IntEps};
    for (int nf = 1; nf <= 6; nf++)
      {
        const Operator OL2nsp{g, CL2Tnsp{nf}, IntEps};
        const Operator OL2t = OL2nsp + 6 * OL2ps;
        std::map<int, Operator> CLNNLOnf;
        CLNNLOnf.insert({DISNCBasis::CNS, OL2nsp});
        CLNNLOnf.insert({DISNCBasis::CS,  OL2t});
        CLNNLOnf.insert({DISNCBasis::CG,  OL2g});
        CLNNLO.insert({nf, CLNNLOnf});
      }

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
          FObj.ConvBasis.insert({k, (k == 0 ? DISNCBasis{EffCh} : DISNCBasis{k, EffCh[k-1]})});
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLLO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLNLO}});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), CLNNLO.at(nf)}});
        }
      return FObj;
    };
    t.stop();

    return FLObj;
  }

  //_____________________________________________________________________________
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF3NCObjectsZMT(Grid                const& g,
                                                                                                              std::vector<double> const& Thresholds,
                                                                                                              double              const& IntEps)
  {
    // Initalise DGLAP objects need for scale variations
    const auto PDFObj = apfel::InitializeDglapObjectsQCDT(g, Thresholds);

    report("Initializing StructureFunctionObjects for F3 Zero Mass for SIA...\n");
    warning("InitializeF3NCObjectsZMT", "NNLO corrections currently unavailable");
    Timer t;

    // ===============================================================
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // LO
    std::map<int, Operator> C3LO;
    C3LO.insert({DISNCBasis::CNS, Id});
    C3LO.insert({DISNCBasis::CS,  Id});
    C3LO.insert({DISNCBasis::CG,  Zero});

    // NLO
    std::map<int, Operator> C3NLO;
    const Operator O31ns{g, C31Tns{}, IntEps};
    C3NLO.insert({DISNCBasis::CNS, O31ns});
    C3NLO.insert({DISNCBasis::CS,  O31ns});
    C3NLO.insert({DISNCBasis::CG,  Zero});

    // NNLO
    std::map<int, std::map<int, Operator>> C3NNLO;
    for (int nf = 1; nf <= 6; nf++)
      {
        std::map<int, Operator> C3NNLOnf;
        C3NNLOnf.insert({DISNCBasis::CNS, Zero});
        C3NNLOnf.insert({DISNCBasis::CS,  Zero});
        C3NNLOnf.insert({DISNCBasis::CG,  Zero});
        C3NNLO.insert({nf, C3NNLOnf});
      }
    /*
        // NNLO
        std::map<int, std::map<int, Operator>> C3NNLO;
        for (int nf = 1; nf <= 6; nf++)
          {
            const Operator O32nsm{g, C32Tnsm{nf}, IntEps};
            const Operator O32t = O32nsm;
            std::map<int, Operator> C3NNLOnf;
            C3NNLOnf.insert({DISNCBasis::CNS, O32nsm});
            C3NNLOnf.insert({DISNCBasis::CS,  O32t});
            C3NNLOnf.insert({DISNCBasis::CG,  Zero});
            C3NNLO.insert({nf, C3NNLOnf});
          }
    */
    // Vector of distributions to skip
    const std::vector<int> skip = {1, 3, 5, 7, 9, 11};

    // Define object of the structure containing the DglapObjects
    const auto F3Obj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects
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
          FObj.ConvBasis.insert({k, (k == 0 ? DISNCBasis{EffCh} : DISNCBasis{k, EffCh[k-1]})});
          FObj.C0.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C3LO}});
          FObj.C1.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C3NLO}});
          FObj.C2.insert({k, Set<Operator>{FObj.ConvBasis.at(k), C3NNLO.at(nf)}});
        }
      return FObj;
    };
    t.stop();

    return F3Obj;
  }

  //_____________________________________________________________________________
  std::map<int, Observable<>> BuildStructureFunctions(std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> const& FObj,
                                                      std::function<std::map<int, double>(double const&, double const&)>                 const& InDistFunc,
                                                      int                                                                                const& PerturbativeOrder,
                                                      std::function<double(double const&)>                                               const& Alphas,
                                                      std::function<std::vector<double>(double const&)>                                  const& Couplings,
                                                      double                                                                             const& xiR,
                                                      double                                                                             const& xiF)
  {
    // Scale variation factors
    const double tR = 2 * log(xiR);
    const double tF = 2 * log(xiF);

    // Call FObj at energy 1 to use it for those quantities that do
    // not depend on Q.
    const StructureFunctionObjects FObj1 = FObj(1, Couplings(1));

    // Get grid.
    Grid const& g = FObj1.C0.at(1).at(0).GetGrid();

    // Get skip vector.
    const std::vector<int> skip = FObj1.skip;

    // Cycle over the key of the convolution basis map.
    std::map<int, Observable<>> F;
    for (auto it = FObj1.ConvBasis.begin(); it != FObj1.ConvBasis.end(); ++it)
      {
        // Structure function index.
        const int k = it->first;

        // Create Observable
        Observable<> Obs{};

        // Split the computation into different cases to avoid
        // computing more convolutions than necessary.
        if (xiF == 1)
          {
            // Declare and then allocate coefficient function
            // according to whether renormalisation-scale-variation
            // terms have to be included or not.
            std::function<Set<Operator>(double const&)> Cf;

            if (xiR == 1)
              Cf = [=] (double const& Q) -> Set<Operator>
              {
                const StructureFunctionObjects FObjQ = FObj(Q, Couplings(Q));
                const double cp  = Alphas(Q) / FourPi;
                const double cp2 = cp * cp;
                const double cp3 = cp * cp2;
                Set<Operator> CoefFuncs = FObjQ.C0.at(k);
                if (PerturbativeOrder > 0)
                  CoefFuncs += cp * FObjQ.C1.at(k);
                if (PerturbativeOrder > 1)
                  CoefFuncs += cp2 * FObjQ.C2.at(k);
                if (PerturbativeOrder > 2)
                  CoefFuncs += cp3 * FObjQ.C3.at(k);
                return CoefFuncs;
              };
            else
              Cf = [=] (double const& Q) -> Set<Operator>
              {
                const StructureFunctionObjects FObjQ = FObj(Q, Couplings(Q));
                const double cp  = Alphas(xiR * Q) / FourPi;
                const double cp2 = cp * cp;
                const double cp3 = cp * cp2;
                const double b0  = beta0qcd(FObjQ.nf);
                const double b1  = beta1qcd(FObjQ.nf);
                Set<Operator> CoefFuncs = FObjQ.C0.at(k);
                if (PerturbativeOrder > 0)
                  CoefFuncs += cp * FObjQ.C1.at(k);
                if (PerturbativeOrder > 1)
                  CoefFuncs += cp2 * ( FObjQ.C2.at(k) + ( tR * b0 ) * FObjQ.C1.at(k) );
                if (PerturbativeOrder > 2)
                  CoefFuncs += cp3 * ( FObjQ.C3.at(k) + ( 2 * tR * b0 ) * FObjQ.C2.at(k) + ( tR * ( b1 + pow(b0, 2) * tR ) ) * FObjQ.C1.at(k) );
                return CoefFuncs;
              };

            // Define distribution-function functions
            const auto DistF = [=, &g] (double const& Q) -> Set<Distribution>
            {
              return Set<Distribution>{FObj(Q, Couplings(Q)).ConvBasis.at(k), DistributionMap(g, InDistFunc, xiF * Q, skip)};
            };

            // Add pair to the observable
            Obs.AddConvolutionPair(Cf, DistF);
          }
        else
          {
            // Declare vectors of coefficient functions and PDFs, each
            // one having as many entried a perturbative orders
            // needed.
            std::vector<std::function<Set<Operator>(double const&)>> Cfv(PerturbativeOrder + 1);
            std::vector<std::function<Set<Distribution>(double const&)>> DistFv(PerturbativeOrder + 1);

            // LO
            Cfv[0] = [=] (double const& Q) -> Set<Operator> { return FObj(Q, Couplings(Q)).C0.at(k); };
            DistFv[0] = [=, &g] (double const& Q) -> Set<Distribution>
            {
              return Set<Distribution>{FObj(Q, Couplings(Q)).ConvBasis.at(k), DistributionMap(g, InDistFunc, xiF * Q, skip)};
            };

            // NLO
            if (PerturbativeOrder > 0)
              {
                Cfv[1] = [=] (double const& Q) -> Set<Operator> { return Alphas(xiR * Q) / FourPi * FObj(Q, Couplings(Q)).C1.at(k); };
                DistFv[1] = [=, &g] (double const& Q) -> Set<Distribution>
                {
                  const StructureFunctionObjects FObjQ = FObj(Q, Couplings(Q));
                  const double cp = Alphas(xiR * Q) / FourPi;
                  const ConvolutionMap CMap = FObjQ.P.SplittingFunctions.at(0).GetMap();
                  const Set<Distribution> F{CMap, DistributionMap(g, InDistFunc, xiF * Q)};
                  const Set<Distribution> P0F = FObjQ.P.SplittingFunctions.at(0) * F;
                  const Set<Distribution> D1  = ( - cp * tF ) * P0F;
                  return Set<Distribution>{FObjQ.ConvBasis.at(k), D1.GetObjects()};
                };
              }

            // NNLO
            if (PerturbativeOrder > 1)
              {
                Cfv[2] = [=] (double const& Q) -> Set<Operator>
                {
                  const StructureFunctionObjects FObjQ = FObj(Q, Couplings(Q));
                  return pow(Alphas(xiR * Q) / FourPi, 2) * ( FObjQ.C2.at(k) + ( tR * beta0qcd(FObjQ.nf) ) * FObjQ.C1.at(k) );
                };
                DistFv[2] = [=, &g] (double const& Q) -> Set<Distribution>
                {
                  const StructureFunctionObjects FObjQ = FObj(Q, Couplings(Q));
                  const double cp2 = pow(Alphas(xiR * Q) / FourPi, 2);
                  const double b0  = beta0qcd(FObjQ.nf);
                  const ConvolutionMap CMap = FObjQ.P.SplittingFunctions.at(0).GetMap();
                  const Set<Distribution> F{CMap, DistributionMap(g, InDistFunc, xiF * Q)};
                  const Set<Distribution> P0F   = FObjQ.P.SplittingFunctions.at(0) * F;
                  const Set<Distribution> P1F   = FObjQ.P.SplittingFunctions.at(1) * F;
                  const Set<Distribution> P0P0F = FObjQ.P.SplittingFunctions.at(0) * P0F;
                  const Set<Distribution> D2    = cp2 * ( ( - tF ) * P1F + ( pow(tF, 2) / 2 ) * P0P0F + ( - tF * tR * b0 + pow(tF, 2) * b0 / 2 ) * P0F );
                  return Set<Distribution>{FObjQ.ConvBasis.at(k), D2.GetObjects()};
                };
              }

            // NNNLO
            if (PerturbativeOrder > 2)
              {
                Cfv[3] = [=] (double const& Q) -> Set<Operator>
                {
                  const StructureFunctionObjects FObjQ = FObj(Q, Couplings(Q));
                  return pow(Alphas(xiR * Q) / FourPi, 3) * ( FObjQ.C3.at(k) + ( 2 * tR * beta0qcd(FObjQ.nf) ) * FObjQ.C2.at(k)
                                                              + ( tR * ( beta1qcd(FObjQ.nf) + pow(beta0qcd(FObjQ.nf), 2) * tR ) ) * FObjQ.C1.at(k) );
                };
                DistFv[3] = [=, &g] (double const& Q) -> Set<Distribution>
                {
                  const StructureFunctionObjects FObjQ = FObj(Q, Couplings(Q));
                  const double cp3 = pow(Alphas(xiR * Q) / FourPi, 3);
                  const double b0  = beta0qcd(FObjQ.nf);
                  const double b1  = beta1qcd(FObjQ.nf);
                  const ConvolutionMap CMap = FObjQ.P.SplittingFunctions.at(0).GetMap();
                  const Set<Distribution> F{CMap, DistributionMap(g, InDistFunc, xiF * Q)};
                  const Set<Distribution> P0F     = FObjQ.P.SplittingFunctions.at(0) * F;
                  const Set<Distribution> P1F     = FObjQ.P.SplittingFunctions.at(1) * F;
                  const Set<Distribution> P2F     = FObjQ.P.SplittingFunctions.at(2) * F;
                  const Set<Distribution> P0P0F   = FObjQ.P.SplittingFunctions.at(0) * P0F;
                  const Set<Distribution> P0P0P0F = FObjQ.P.SplittingFunctions.at(0) * P0P0F;
                  const Set<Distribution> P0P1F   = 0.5 * ( FObjQ.P.SplittingFunctions.at(0) * P1F + FObjQ.P.SplittingFunctions.at(1) * P0F );
                  const Set<Distribution> D3      = cp3 * ( ( pow(tF, 2) * b1 / 2 - tF * tR * b1 - pow(tF, 3) * pow(b0, 2) / 3 + pow(tF, 2) * tR * pow(b0, 2) - tF * pow(tR, 2) * pow(b0, 2) ) * P0F
                                                            + ( pow(tF, 2) * b0 - 2 * tF * tR * b0 ) * P1F
                                                            + ( - tF ) * P2F
                                                            + ( - pow(tF, 3) * b0 / 2 + pow(tF, 2) * tR * b0 ) * P0P0F
                                                            + ( - pow(tF, 3) / 6 ) * P0P0P0F
                                                            + pow(tF, 2) * P0P1F );
                  return Set<Distribution>{FObjQ.ConvBasis.at(k), D3.GetObjects()};
                };
              }

            // Combine coefficient functions and PDFs
            for (int i = 0; i <= PerturbativeOrder; i++)
              for (int j = 0; j <= PerturbativeOrder - i; j++)
                Obs.AddConvolutionPair(Cfv[i], DistFv[j]);
          }

        // Finally insert Observable
        F.insert({k, Obs});
      }
    return F;
  }

  //_____________________________________________________________________________
  std::map<int, Observable<>> BuildStructureFunctions(std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> const& FObj,
                                                      std::function<double(int const&, double const&, double const&)>                    const& InDistFunc,
                                                      int                                                                                const& PerturbativeOrder,
                                                      std::function<double(double const&)>                                               const& Alphas,
                                                      std::function<std::vector<double>(double const&)>                                  const& Couplings,
                                                      double                                                                             const& xiR,
                                                      double                                                                             const& xiF)
  {
    const auto InDistFuncMap = [=] (double const& x, double const& Q) -> std::map<int, double>
    {
      std::map<int, double> DistMap;
      for (int i = 0; i <= 12; i++)
        DistMap.insert({i,InDistFunc(i, x, Q)});
      return DistMap;
    };
    return BuildStructureFunctions(FObj, InDistFuncMap, PerturbativeOrder, Alphas, Couplings, xiR, xiF);
  }

  //_____________________________________________________________________________
  Distribution BuildStructureFunctions(StructureFunctionObjects    const& FObjQ,
                                       std::map<int, Distribution> const& InDistFuncQ,
                                       int                         const& PerturbativeOrder,
                                       double                      const& AlphasQ,
                                       int                         const& k,
                                       double                      const& xiR,
                                       double                      const& xiF)
  {
    // Scale variation factors
    const double tR = 2 * log(xiR);
    const double tF = 2 * log(xiF);

    // Gather coupling constant and its powers, beta function
    // coefficients, and splitting functions convoluted with the input
    // PDFs.
    const double cp  = AlphasQ / FourPi;
    const double cp2 = cp * cp;
    const double cp3 = cp * cp2;
    const double b0  = beta0qcd(FObjQ.nf);
    const double b1  = beta1qcd(FObjQ.nf);

    const ConvolutionMap CMap = FObjQ.P.SplittingFunctions.at(0).GetMap();
    const Set<Distribution> F{CMap, InDistFuncQ};
    const Set<Distribution> P0F     = FObjQ.P.SplittingFunctions.at(0) * F;
    const Set<Distribution> P1F     = FObjQ.P.SplittingFunctions.at(1) * F;
    const Set<Distribution> P2F     = FObjQ.P.SplittingFunctions.at(2) * F;
    const Set<Distribution> P0P0F   = FObjQ.P.SplittingFunctions.at(0) * P0F;
    const Set<Distribution> P0P1F   = 0.5 * ( FObjQ.P.SplittingFunctions.at(0) * P1F + FObjQ.P.SplittingFunctions.at(1) * P0F );
    const Set<Distribution> P0P0P0F = FObjQ.P.SplittingFunctions.at(0) * P0P0F;

    // Allocate vectors of coeffient functions and PDFs
    std::vector<Set<Operator>> Cfv(PerturbativeOrder + 1);
    std::vector<Set<Distribution>> DistFv(PerturbativeOrder + 1);

    // LO
    Cfv[0] = FObjQ.C0.at(k);
    DistFv[0] = Set<Distribution> {FObjQ.ConvBasis.at(k), InDistFuncQ};

    // NLO
    if (PerturbativeOrder > 0)
      {
        Cfv[1] = cp * FObjQ.C1.at(k);
        DistFv[1] = Set<Distribution> {FObjQ.ConvBasis.at(k), (( - cp * tF ) * P0F).GetObjects()};
      }

    // NNLO
    if (PerturbativeOrder > 1)
      {
        Cfv[2] = cp2 * ( FObjQ.C2.at(k) + ( tR * b0 ) * FObjQ.C1.at(k) );
        DistFv[2] = Set<Distribution> {FObjQ.ConvBasis.at(k), (cp2 * ( ( - tF ) * P1F + ( pow(tF, 2) / 2 ) * P0P0F + ( - tF * tR * b0 + pow(tF, 2) * b0 / 2 ) * P0F )).GetObjects()};
      }

    // NNNLO
    if (PerturbativeOrder > 2)
      {
        Cfv[3] = cp3 * ( FObjQ.C3.at(k) + ( 2 * tR * b0 ) * FObjQ.C2.at(k) + ( tR * ( b1 + pow(b0, 2) * tR ) ) * FObjQ.C1.at(k) );
        DistFv[3] = Set<Distribution> {FObjQ.ConvBasis.at(k), (cp3 * ( ( pow(tF, 2) * b1 / 2 - tF * tR * b1 - pow(tF, 3) * pow(b0, 2) / 3 + pow(tF, 2) * tR * pow(b0, 2) - tF * pow(tR, 2) * pow(b0, 2) ) * P0F
                                                                       + ( pow(tF, 2) * b0 - 2 * tF * tR * b0 ) * P1F
                                                                       + ( - tF ) * P2F
                                                                       + ( - pow(tF, 3) * b0 / 2 + pow(tF, 2) * tR * b0 ) * P0P0F
                                                                       + ( - pow(tF, 3) / 6 ) * P0P0P0F
                                                                       + pow(tF, 2) * P0P1F )).GetObjects()
                                      };
      }

    // Convolute coefficient function with set of distributions
    Set<Distribution> SF = Cfv[0] * DistFv[0];
    for (int i = 0; i <= PerturbativeOrder; i++)
      for (int j = 0; j <= PerturbativeOrder - i; j++)
        if (i != 0 && j != 0)
          SF += Cfv[i] * DistFv[j];

    // Combine set of distributions and return
    return SF.Combine();
  }

  //_____________________________________________________________________________
  std::map<int, Distribution> BuildStructureFunctions(StructureFunctionObjects    const& FObjQ,
                                                      std::map<int, Distribution> const& InDistFuncQ,
                                                      int                         const& PerturbativeOrder,
                                                      double                      const& AlphasQ,
                                                      double                      const& xiR,
                                                      double                      const& xiF)
  {
    // Cycle over the key of the convolution basis map.
    std::map<int, Distribution> F;
    for (auto it = FObjQ.ConvBasis.begin(); it != FObjQ.ConvBasis.end(); ++it)
      // Push "Distribution".
      F.insert({it->first, BuildStructureFunctions(FObjQ, InDistFuncQ, PerturbativeOrder, AlphasQ, it->first, xiR, xiF)});

    return F;
  }
}
