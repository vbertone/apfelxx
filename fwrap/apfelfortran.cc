//
// apfelxx Fortran wrapper
//
// C-linkage shim called by the apfel_fortran Fortran module
// (fwrap/apfel_fortran.f90) via bind(C) interfaces. Mirrors the
// procedural, stateful interface of the old Fortran-native APFEL
// (APFELfwevol.h) on top of apfelxx's C++ evolution machinery.
//

#include "apfel/alphaqcd.h"
#include "apfel/config.h"
#include "apfel/dglapbuilder.h"
#include "apfel/grid.h"
#include "apfel/lhtoypdfs.h"
#include "apfel/rotations.h"
#include "apfel/set.h"
#include "apfel/subgrid.h"
#include "apfel/tabulateobject.h"

#ifdef WITH_LHAPDF
#include <LHAPDF/LHAPDF.h>
#endif

#include <cstdlib>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace
{
  // Parameters of one x-space subgrid, staged by SetGridParameters and
  // only turned into an apfel::SubGrid (which has no default
  // constructor) once InitializeAPFEL is called.
  struct GridParams
  {
    int nx = 0;
    double xmin = 0;
    int id = 0;
  };

  // Mirrors the old Fortran commons blocks: configuration staged by the
  // SetXxx calls, consumed by InitializeAPFEL.
  struct FortranState
  {
    // Defaults mirror apfel::EvolutionSetup's own defaults
    // (inc/apfel/evolutionsetup.h), so InitializeAPFEL/EvolveAPFEL work
    // out of the box without any SetXxx calls.
    int PerturbativeOrder = 2;
    std::vector<double> QuarkMasses = {0, 0, 0, std::sqrt(2.), 4.5, 175};
    std::vector<double> QuarkThresholds = {0, 0, 0, std::sqrt(2.), 4.5, 175};
    double AlphaQCDRef = 0.35;
    double QRef = std::sqrt(2.);
    double Qmin = 1;
    double Qmax = 1000;
    int nQg = 50;
    int InterDegreeQ = 3;
    std::vector<GridParams> subgridParams =
    {
      {100, 1e-5, 3}, {100, 1e-1, 3}, {100, 6e-1, 3}, {80, 8.5e-1, 5}
    };

    // Renormalisation/factorisation scale-variation ratio xi = muR/muF,
    // passed straight through to BuildDglap (which defaults to 1 itself).
    double xi = 1;

    // LHAPDF member index, consulted by SetPDFSet when it next builds an
    // LHAPDF-backed InSet. Has no effect on the built-in "toyLH" set.
    int replica = 0;

    // Built by InitializeAPFEL (grid + evolution-kernel integrals only,
    // no PDFs involved yet, exactly like the old InitializeAPFEL.f)
    std::unique_ptr<const apfel::Grid> g;
    std::function<double(double const&)> as;
    std::map<int, apfel::DglapObjects> DglapObj;
    bool initialised = false;

    // Set by SetPDFSet: the initial-scale condition, in evolution basis
    // (BuildDglap's InDistFunc contract). Defaults to the built-in toy
    // set so EvolveAPFEL/CachePDFsAPFEL work without calling SetPDFSet.
    std::function<std::map<int, double>(double const&, double const&)> InSet = apfel::LHToyPDFs;

    // Built lazily by EvolveAPFEL/CachePDFsAPFEL from InSet; rebuilt
    // whenever SetPDFSet is called again or Q0 changes.
    std::unique_ptr<apfel::Dglap<apfel::Distribution>> dglap;
    double dglapQ0 = -1;

    // Result of the last EvolveAPFEL (always a direct, un-cached
    // dglap->Evaluate), rotated to physical basis, ready for xPDF.
    std::map<int, apfel::Distribution> cachedPhys;
    bool evolved = false;

    // Built by CachePDFsAPFEL: wraps `dglap` in a TabulateObject so that
    // xPDFxQ can look up any (x, Q) by interpolation instead of a fresh
    // ODE solve. Independent of EvolveAPFEL/xPDF.
    std::unique_ptr<const apfel::TabulateObject<apfel::Set<apfel::Distribution>>> tab;
  };

  FortranState& State()
  {
    static FortranState state;
    return state;
  }

  // Mirrors the old API's "print message and call exit(-10)" convention
  // for fatal configuration errors.
  [[noreturn]] void FatalError(std::string const& where, std::string const& message)
  {
    std::cerr << "In " << where << ":\n" << message << std::endl;
    std::exit(-10);
  }

  // Builds (or reuses) the Dglap evolution object for the given Q0.
  // Shared by EvolveAPFEL and CachePDFsAPFEL, which otherwise both need
  // the same "solve the DGLAP equations from this Q0" starting point.
  void EnsureDglap(FortranState& state, std::string const& where, double q0)
  {
    if (!state.initialised)
      FatalError(where, "APFEL has not been initialised. Call InitializeAPFEL first.");
    if (!state.InSet)
      FatalError(where, "No PDF set has been defined. Call SetPDFSet first.");

    if (!state.dglap || state.dglapQ0 != q0)
      {
        state.dglap = apfel::BuildDglap(state.DglapObj, state.InSet, q0, state.PerturbativeOrder, state.as, state.xi);
        state.dglapQ0 = q0;
      }
  }
}

extern "C"
{
  bool apfelxxf_checkapfel()
  {
    return true;
  }

  void apfelxxf_setperturbativeorder(int pto)
  {
    State().PerturbativeOrder = pto;
  }

  void apfelxxf_setpolemasses(double mc, double mb, double mt)
  {
    State().QuarkMasses = {0, 0, 0, mc, mb, mt};
    State().QuarkThresholds = {0, 0, 0, mc, mb, mt};
  }

  void apfelxxf_setalphaqcdref(double alpharef, double qref)
  {
    State().AlphaQCDRef = alpharef;
    State().QRef = qref;
  }

  void apfelxxf_setqlimits(double qmin, double qmax)
  {
    State().Qmin = qmin;
    State().Qmax = qmax;
  }

  void apfelxxf_setnumberofgrids(int n)
  {
    if (n <= 0)
      FatalError("SetNumberOfGrids", "The number of grids must be positive.");
    State().subgridParams.assign(n, GridParams{});
  }

  void apfelxxf_setgridparameters(int i, int np, int deg, double x)
  {
    auto& sgp = State().subgridParams;
    if (i < 1 || i > (int) sgp.size())
      FatalError("SetGridParameters", "Grid index out of range. Call SetNumberOfGrids first.");
    sgp[i - 1] = GridParams{np, x, deg};
  }

  void apfelxxf_setqgridparameters(int npq, int degq)
  {
    State().nQg = npq;
    State().InterDegreeQ = degq;
  }

  void apfelxxf_initializeapfel()
  {
    FortranState& state = State();

    if (state.subgridParams.empty())
      FatalError("InitializeAPFEL", "No x-space grid defined. Call SetNumberOfGrids and SetGridParameters first.");

    // Build the x-space grid
    std::vector<apfel::SubGrid> sg;
    for (auto const& gp : state.subgridParams)
      sg.emplace_back(gp.nx, gp.xmin, gp.id);
    state.g = std::unique_ptr<const apfel::Grid>(new apfel::Grid{sg});

    // Build and tabulate the strong coupling
    apfel::AlphaQCD a{state.AlphaQCDRef, state.QRef, state.QuarkMasses, state.QuarkThresholds, state.PerturbativeOrder};
    const auto Alphas = std::make_shared<const apfel::TabulateObject<double>>(a, 2 * state.nQg, state.Qmin - 0.1, state.Qmax + 1, state.InterDegreeQ);
    state.as = [Alphas] (double const& mu) -> double{ return Alphas->Evaluate(mu); };

    // Precompute the perturbative coefficients of the splitting
    // functions and matching conditions (no PDFs involved yet)
    state.DglapObj = apfel::InitializeDglapObjectsQCDOme(*state.g, state.QuarkMasses, state.QuarkThresholds);

    state.initialised = true;
  }

  double apfelxxf_alphaqcd(double q)
  {
    if (!State().initialised)
      FatalError("alphaQCD", "APFEL has not been initialised. Call InitializeAPFEL first.");
    return State().as(q);
  }

  int apfelxxf_getperturbativeorder()
  {
    return State().PerturbativeOrder;
  }

  void apfelxxf_setrenfacratio(double xi)
  {
    FortranState& state = State();
    state.xi = xi;
    // The dglap object was built with the previous xi: force a rebuild.
    state.dglap.reset();
    state.tab.reset();
    state.evolved = false;
  }

  void apfelxxf_setreplica(int irep)
  {
    State().replica = irep;
  }

  void apfelxxf_setpdfset(char const* name, int length)
  {
    FortranState& state = State();
    std::string const pdfset(name, length);

    if (pdfset == "toyLH")
      state.InSet = apfel::LHToyPDFs;
    else
      {
#ifdef WITH_LHAPDF
        auto pdf = std::shared_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(pdfset, state.replica));
        state.InSet = [pdf] (double const& x, double const& Q0) -> std::map<int, double>
        {
          return apfel::PhysToQCDEv(pdf->xfxQ(x, Q0));
        };
#else
        FatalError("SetPDFSet", "Unknown PDF set \"" + pdfset + "\": LHAPDF support was not built into apfelxx, only \"toyLH\" is available.");
#endif
      }

    // Invalidate anything that depended on the previous PDF set
    state.dglap.reset();
    state.tab.reset();
    state.evolved = false;
  }

  void apfelxxf_evolveapfel(double q0, double q)
  {
    FortranState& state = State();
    EnsureDglap(state, "EvolveAPFEL", q0);

    state.cachedPhys = apfel::QCDEvToPhys(state.dglap->Evaluate(q).GetObjects());
    state.evolved = true;
  }

  void apfelxxf_cachepdfsapfel(double q0)
  {
    FortranState& state = State();
    EnsureDglap(state, "CachePDFsAPFEL", q0);

    state.tab = std::unique_ptr<const apfel::TabulateObject<apfel::Set<apfel::Distribution>>>(new apfel::TabulateObject<apfel::Set<apfel::Distribution>> {*state.dglap, state.nQg, state.Qmin, state.Qmax, state.InterDegreeQ});
  }

  double apfelxxf_xpdf(int i, double x)
  {
    FortranState& state = State();
    if (i < -6 || i > 6)
      FatalError("xPDF", "Invalid PDF index, i = " + std::to_string(i));
    if (!state.evolved)
      FatalError("xPDF", "PDFs have not been evolved. Call EvolveAPFEL first.");
    return state.cachedPhys.at(i).Evaluate(x);
  }

  double apfelxxf_xpdfxq(int i, double x, double q)
  {
    FortranState& state = State();
    if (i < -6 || i > 6)
      FatalError("xPDFxQ", "Invalid PDF index, i = " + std::to_string(i));
    if (!state.tab)
      FatalError("xPDFxQ", "PDFs have not been cached. Call CachePDFsAPFEL first.");
    return apfel::QCDEvToPhys(state.tab->EvaluateMapxQ(x, q)).at(i);
  }
}
