// The **evolsetup** class contains all the parameters needed to perform the DGLAP evolution.
// The constructor initializes a set of default parameters which correspond to the LH
// benchmark parameters at NNLO.

// Include only once
#pragma once

// Standard libraries
#include <string>
#include <cmath>
#include <vector>

using namespace std;

// All classes are namespaced in **apfel::**
namespace apfel {

  // Enumerators to facilitate the interpretation

  // Flavour scheme
  enum flavourscheme {
    VFNS,
    FFNS
  };

  // Evolution theory
  enum theory {
    QCD,
    QED,
    QUniD
  };

  // Logarithmic accuracy od the small-x resummed evolution (inactive is the small-x evolution is not enabled)
  enum logaccuracy {
    LL,
    NLL
  };

  // Solution of the coupling RG equations
  enum alphaevolution {
    exact,
    expanded,
    lambda
  };

  // Solution of the DGLAP equation
  enum pdfevolution {
    exactmu,
    exactalpha,
    exapandalpha,
    truncated
  };

  // Renormalization scheme used for the heavy quark masses
  enum massrenscheme {
    POLE,
    MSBAR
  };

  // Structure for the **SubGrid** parameters
  struct gridparams {
    int    nx;      // Number of nodes
    double xmin;    // Minimum value of x
    int    id;      // Interpolation degree
    double *xgext;  // External grid (if any)
  };

  // Begin of the class **evolsetup**
  class evolsetup {
  private:
    // Attributes of the class
    bool               _WelcomeMessage;      // Welcome message (ON, OFF)
    double             _Qmin;                // Lower bound of the evolution range
    double             _Qmax;                // Upper bound of the evolution range
    int                _PerturbativeOrder;   // Perturbative order of the evolution (LO, NLO, NNLO)
    flavourscheme      _FlavourScheme;       // Flavour scheme (FFNS, VFNS)
    int                _Nf_FF;               // Numeber of active flavours in the FFNS
    theory             _Theory;              // Evolution theory (QCD, QED, QUniD)
    bool               _FastEvolution;       // Fast evolution (ON, OFF)
    bool               _TimeLikeEvolution;   // Time-like evolution (ON, OFF)
    bool               _PolarizedEvolution;  // Polarized evolution (ON, OFF)
    bool               _SmallxResummation;   // Small-x resummation (ON, OFF)
    logaccuracy        _LogAccuracy;         // Logarithmic accuracy of the small-x resummation (LL, NLL)
    double             _AlphaQCDRef;         // Reference value of the QCD coupling
    double             _QQCDRef;             // reference scale of the QCD coupling
    double             _AlphaQEDRef;         // Reference value of the QED coupling
    double             _QQEDRef;             // Reference scale of the QCD coupling
    double             _LambdaQCDRef;        // Value of LambdaQCD
    int                _nLambdaQCDRef;       // Number of flavours of LambdaQCD
    double             _EpsilonTruncation;   // Value of the (small) parameter to compute the numerical derivatives
    alphaevolution     _AlphaEvolution;      // Solution of the RGE of the couplings
    pdfevolution       _PDFEvolution;        // Solution of the DGLAP equations
    double             _RenFacRatio;         // ratio muF / muR
    massrenscheme      _MassRenScheme;       // Renormalization scheme for the heavy quark masses
    double             _MCharm;              // Reference value of the charm mass
    double             _MBottom;             // Reference value of the bottom mass
    double             _MTop;                // Reference value of the top mass
    double             _kThCharm;            // Displacement of the charm threshold
    double             _kThBottom;           // Displacement of the bottom threshold
    double             _kThTop;              // Displacement of the top threshold
    double             _QRefMCharm;          // Reference scale of the charm mass
    double             _QRefMBottom;         // Reference scale of the bottom mass
    double             _QRefMTop;            // Reference scale of the top mass
    double             _TauMass;             // Mass of the tau lepton
    bool               _MassRunning;         // Running of the heavy quark masses (ON, OFF)
    int                _MaxFlavourPDFs;      // Maximum number of flavours allowed in the PDF evolution
    int                _MaxFlavourAlpha;     // Maximum number of flavours allowed in the coupling evolution
    string             _PDFSet;              // Name of the input PDF set
    int                _Replica;             // PDF replica/member
    bool               _EvolutionOperator;   // Compute evolution operator (ON, OFF)
    bool               _LeptonEvolution;     // Lepton evolution (ON, OFF)
    double             _Q2minLHA;            // Parametes of the LHAPDF (x,Q^2) grid
    double             _Q2maxLHA;
    double             _xminLHA;
    double             _xmLHA;
    double             _xmaxLHA;
    int                _nxLHA;
    int                _nxmLHA;
    int                _nQ2LHA;
    double             _Lambda2LHA;
    int                _nQ2g;                // Number of Q^2 intervals for the PDF chaching
    int                _InterDegreeQ;        // Interpolation degree on the Q^2 grid
    bool               _Locked;              // Subgrid locked (ON, OFF)
    vector<gridparams> _GridParams;          // Vector of the parameters of the SubGrids
    int                _GaussPoints;         // Number of points used for the Gauss integrations

  public:
    // Constructor
    evolsetup():
      _WelcomeMessage(true),
      _Qmin(1), _Qmax(10000),
      _PerturbativeOrder(0),
      _FlavourScheme(VFNS),
      _Nf_FF(3),
      _Theory(QCD),
      _FastEvolution(true),
      _TimeLikeEvolution(false),
      _PolarizedEvolution(false),
      _SmallxResummation(false), _LogAccuracy(NLL),
      _AlphaQCDRef(0.35), _QQCDRef(sqrt(2)),
      _AlphaQEDRef(7.496252e-3), _QQEDRef(1.777),
      _LambdaQCDRef(0.220), _nLambdaQCDRef(5),
      _EpsilonTruncation(1e-2),
      _AlphaEvolution(exact),
      _PDFEvolution(exactalpha),
      _RenFacRatio(1),
      _MassRenScheme(POLE),
      _MCharm(sqrt(2)), _MBottom(4.5), _MTop(175),
      _kThCharm(1), _kThBottom(1), _kThTop(1),
      _QRefMCharm(_MCharm), _QRefMBottom(_MTop), _QRefMTop(_MTop),
      _TauMass(1.777),
      _MassRunning(true),
      _MaxFlavourPDFs(6),
      _MaxFlavourAlpha(6),
      _PDFSet("ToyLH"),
      _Replica(0),
      _EvolutionOperator(false),
      _LeptonEvolution(false),
      _Q2minLHA(1), _Q2maxLHA(1e10),
      _xminLHA(1e-9), _xmLHA(1e-1), _xmaxLHA(1),
      _nxLHA(100), _nxmLHA(50), _nQ2LHA(50),
      _Lambda2LHA(0.0625),
      _nQ2g(70), _InterDegreeQ(3),
      _Locked(true),
      _GaussPoints(3)
    {
      // Assign parameters to the SubGrids
      gridparams gp1 = {80, 1e-5, 3, NULL};
      gridparams gp2 = {50, 1e-1, 5, NULL};
      gridparams gp3 = {40, 8e-1, 5, NULL};
      _GridParams.push_back(gp1);
      _GridParams.push_back(gp2);
      _GridParams.push_back(gp3);
    }

    // Setter functions
    void EnableWelcomeMessage(bool const& WelcomeMessage_);
    void SetQLimits(double const& Qmin_, double const& Qmax_);
    void SetPerturbativeOrder(int const& PerturbativeOrder_);
    void SetVFNS();
    void SetFFNS(int const& Nf_FF_);
    void SetTheory(theory const& Theory_);
    void SetFastEvolution(bool const& FastEvolution_);
    void SetTimeLikeEvolution(bool const& TimeLikeEvolution_);
    void SetPolarizedEvolution(bool const& PolarizedEvolution_);
    void SetSmallxResummation(bool const& SmallxResummation_, logaccuracy const& LogAccuracy_);
    void SetAlphaQCDRef(double const& AlphaQCDRef_, double const& QQCDRef_);
    void SetAlphaQEDRef(double const& AlphaQEDRef_, double const& QQEDRef_);
    void SetLambdaQCDRef(double const& LambdaQCDRef_, int const& nLambdaQCDRef_);
    void SetEpsilonTruncation(double const& EpsilonTruncation_);
    void SetAlphaEvolution(alphaevolution const& AlphaEvolution_);
    void SetPDFEvolution(pdfevolution const& PDFEvolution_);
    void SetRenFacRatio(double const& RenFacRatio_);
    void SetPoleMasses(double const& MCharm_, double const& MBottom_, double const& MTop_);
    void SetMSBarMasses(double const& MCharm_, double const& MBottom_, double const& MTop_);
    void SetMassMatchingScales(double const& kThCharm_, double const& kThBottom_, double const& kThTop_);
    void SetMassScaleReference(double const& MCharm_, double const& MBottom_, double const& MTop_);
    void SetTauMass(double const& TauMass_);
    void EnableMassRunning(bool const& MassRunning_);
    void SetMaxFlavourPDFs(int const& MaxFlavourPDFs_);
    void SetMaxFlavourAlpha(int const& MaxFlavourAlpha_);
    void SetPDFSet(string const& PDFSet_);
    void SetReplica(int const& Replica_);
    void EnableEvolutionOperator(bool const& EvolutionOperator_);
    void EnableLeptonEvolution(bool const& LeptonEvolution_);
    void SetLHgridParameters(int const& nxLHA_, int const& nxmLHA_, double const& xminLHA_, double const& xmLHA_, double const& xmaxLHA_,
			     int const& nQ2LHA_, double const& Q2minLHA_, double const& Q2maxLHA_);
    void SetQGridParameters(int const& nQ2g_, int const& InterDegreeQ_);
    void LockGrids(bool const& Locked_);
    void SetGridParameters(int const& nx_, int const& id_, double const& xmin_);
    void SetGridParameters(int const& nx_, int const& id_, double *xgext_);
    void SetGaussPoints(int const& GaussPoints_);

    // Getter functions
    bool           WelcomeMessage()            const { return _WelcomeMessage; }
    double         Qmin()                      const { return _Qmin; }
    double         Qmax()                      const { return _Qmax; }
    int            PerturbativeOrder()         const { return _PerturbativeOrder; }
    flavourscheme  FlavourScheme()             const { return _FlavourScheme; }
    int            Nf_FF()                     const { return _Nf_FF; }
    theory         Theory()                    const { return _Theory; }
    bool           FastEvolution()             const { return _FastEvolution; }
    bool           TimeLikeEvolution()         const { return _TimeLikeEvolution; }
    bool           PolarizedEvolution()        const { return _PolarizedEvolution; }
    bool           SmallxResummation()         const { return _SmallxResummation; }
    logaccuracy    LogAccuracy()               const { return _LogAccuracy; }
    double         AlphaQCDRef()               const { return _AlphaQCDRef; }
    double         QQCDRef()                   const { return _QQCDRef; }
    double         AlphaQEDRef()               const { return _AlphaQEDRef; }
    double         QQEDRef()                   const { return _QQEDRef; }
    double         LambdaQCDRef()              const { return _LambdaQCDRef; }
    int            nLambdaQCDRef()             const { return _nLambdaQCDRef; }
    double         EpsilonTruncation()         const { return _EpsilonTruncation; }
    alphaevolution AlphaEvolution()            const { return _AlphaEvolution; }
    pdfevolution   PDFEvolution()              const { return _PDFEvolution; }
    double         RenFacRatio()               const { return _RenFacRatio; }
    massrenscheme  MassRenScheme()             const { return _MassRenScheme; }
    double         MCharm()                    const { return _MCharm; }
    double         MBottom()                   const { return _MBottom; }
    double         MTop()                      const { return _MTop; }
    double         MuCharm()                   const { return _kThCharm * _MCharm; }
    double         MuBottom()                  const { return _kThBottom * _MBottom; }
    double         MuTop()                     const { return _kThTop * _MTop; }
    double         QRefMCharm()                const { return _QRefMCharm; }
    double         QRefMBottom()               const { return _QRefMBottom; }
    double         QRefMTop()                  const { return _QRefMTop; }
    double         TauMass()                   const { return _TauMass; }
    bool           MassRunning()               const { return _MassRunning; }
    int            MaxFlavourPDFs()            const { return _MaxFlavourPDFs; }
    int            MaxFlavourAlpha()           const { return _MaxFlavourAlpha; }
    string         PDFSet()                    const { return _PDFSet; }
    int            Replica()                   const { return _Replica; }
    bool           EvolutionOperator()         const { return _EvolutionOperator; }
    bool           LeptonEvolution()           const { return _LeptonEvolution; }
    int            nxLHA()                     const { return _nxLHA; }
    int            nxmLHA()                    const { return _nxmLHA; }
    double         xminLHA()                   const { return _xminLHA; }
    double         xmLHA()                     const { return _xmLHA; }
    double         xmaxLHA()                   const { return _xmaxLHA; }
    int            nQ2LHA()                    const { return _nQ2LHA; }
    double         Q2minLHA()                  const { return _Q2minLHA; }
    double         Q2maxLHA()                  const { return _Q2maxLHA; }
    int            nQ2g()                      const { return _nQ2g; }
    int            InterDegreeQ()              const { return _InterDegreeQ; }
    bool           Locked()                    const { return _Locked; }
    int            nGrids()                    const { return _GridParams.size(); }
    int            GaussPoints()               const { return _GaussPoints; }
    int            nx(int const& ig_)          const;
    int            InterDegree(int const& ig_) const;
    double         xMin(int const& ig_)        const;
    int            NfFinPDF()                  const;
    int            NfIniPDF()                  const;
    int            NfFinAlpha()                const;
    int            NfIniAlpha()                const;

  };

}
