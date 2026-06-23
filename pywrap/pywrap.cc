#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>

#include <apfel/apfelxx.h>
#include <apfel/betaqcd.h>
#include <apfel/gammak.h>
#include <apfel/gammaf.h>
#include <apfel/kcs.h>
#include <apfel/massivecoefficientfunctionsunp_sl.h>

namespace py = pybind11;
using namespace pybind11::literals;

// Helper that registers the common interface of the QGrid<T> template
// for a given specialisation. Specialisation-specific methods (if any)
// can be chained onto the returned object at the call site.
template <typename T>
py::class_<apfel::QGrid<T>> bind_qgrid(py::module_& m, char const* name)
{
  return py::class_<apfel::QGrid<T>>(m, name, "Tabulation grid in the scale Q with interpolation (QGrid<T> specialisation).")
         .def(py::init<int const&, double const&, double const&, int const&, std::vector<double> const&, std::function<double(double const&)> const&, std::function<double(double const&)>>(), "Construct a Q-grid with custom tabulation function and its inverse.", "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "TabFunc"_a, "InvTabFunc"_a)
         .def(py::init<int const&, double const&, double const&, int const&, std::vector<double> const&, double const&>(), "Construct a Q-grid using the default log-log tabulation function with scale Lambda.", "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "Lambda"_a = 0.25)
         .def(py::init<std::vector<double> const&, int const&>(), "Construct a Q-grid from an explicit list of Q nodes.", "Qg"_a, "InterDegree"_a)
         .def("Evaluate", &apfel::QGrid<T>::Evaluate, "Evaluate the tabulated object at Q.", "Q"_a)
         .def("Derive", &apfel::QGrid<T>::Derive, "Evaluate the derivative of the tabulated object at Q.", "Q"_a)
         .def("Integrate", &apfel::QGrid<T>::Integrate, "Integrate the tabulated object between Qa and Qb.", "Qa"_a, "Qb"_a)
         .def("nQ", &apfel::QGrid<T>::nQ, "Return the number of Q nodes.")
         .def("InterDegree", &apfel::QGrid<T>::InterDegree, "Return the interpolation degree.")
         .def("QMin", &apfel::QGrid<T>::QMin, "Return the lower bound in Q.")
         .def("QMax", &apfel::QGrid<T>::QMax, "Return the upper bound in Q.")
         .def("TabFunc", &apfel::QGrid<T>::TabFunc, "Return the tabulation function.")
         .def("GetThresholds", &apfel::QGrid<T>::GetThresholds, "Return the heavy-quark thresholds.")
         .def("GetQGrid", &apfel::QGrid<T>::GetQGrid, "Return the Q-grid nodes.")
         .def("GetFQGrid", &apfel::QGrid<T>::GetFQGrid, "Return the tabulated values on the Q-grid.")
         .def("GetThesholdIndices", &apfel::QGrid<T>::GetThesholdIndices, "Return the indices of the thresholds on the grid.")
         .def("GetQGridValues", &apfel::QGrid<T>::GetQGridValues, "Return the values of the tabulation function on the Q-grid.")
         .def("Interpolant", &apfel::QGrid<T>::Interpolant, "Interpolating function (interpolation weights) in Q.", "tQ"_a, "tau"_a, "fq"_a)
         .def("DerInterpolant", &apfel::QGrid<T>::DerInterpolant, "Derivative of the interpolating function in Q.", "tQ"_a, "tau"_a, "Q"_a)
         .def("IntInterpolant", &apfel::QGrid<T>::IntInterpolant, "Integral of the interpolating function in Q.", "tQ"_a, "tau"_a, "Qa"_a, "Qb"_a)
         .def("Print", &apfel::QGrid<T>::Print, "Print the QGrid object.")
         .def(py::self == py::self)
         .def(py::self != py::self);
}

// Helper that registers the common interface of the TabulateObject<T>
// template (its four constructors) for a given specialisation. The
// interpolation/evaluation methods are inherited from the QGrid<T> base,
// while specialisation-specific evaluators (EvaluatexQ, ...) can be
// chained onto the returned object at the call site.
template <typename T>
py::class_<apfel::TabulateObject<T>, apfel::QGrid<T>> bind_tabulateobject(py::module_& m, char const* name, char const* doc)
{
  return py::class_<apfel::TabulateObject<T>, apfel::QGrid<T>>(m, name, doc)
         .def(py::init<apfel::MatchedEvolution<T>&, int const&, double const&, double const&, int const&, double const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Lambda"_a = 0.25)
         .def(py::init<std::function<T(double const&)> const&, int const&, double const&, double const&, int const&, std::vector<double> const&, double const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "Lambda"_a = 0.25)
         .def(py::init<std::function<T(double const&)> const&, int const&, double const&, double const&, int const&, std::vector<double> const&, std::function<double(double const&)> const&, std::function<double(double const&)> const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "TabFunc"_a, "InvTabFunc"_a)
         .def(py::init<std::function<T(double const&)> const&, std::vector<double> const&, int const&>(), "Object"_a, "Qg"_a, "InterDegree"_a);
}

// Trampoline class template for the virtual MatchedEvolution<T> base,
// enabling Python subclasses to override the evolution/matching methods.
template <typename T>
class PyMatchedEvolution: public apfel::MatchedEvolution<T>
{
public:
  using apfel::MatchedEvolution<T>::MatchedEvolution;
  T EvolveObject(int const& nf, double const& mu02, double const& mu2, T const& Obj0) const override
  {
    PYBIND11_OVERRIDE(T, apfel::MatchedEvolution<T>, EvolveObject, nf, mu02, mu2, Obj0);
  };
  T MatchObject(bool const& Up, int const& nf, T const& Obj) const override
  {
    PYBIND11_OVERRIDE_PURE(T, apfel::MatchedEvolution<T>, MatchObject, Up, nf, Obj);
  };
  T Derivative(int const& nf, double const& Mu, T const& Obj) const override
  {
    PYBIND11_OVERRIDE_PURE(T, apfel::MatchedEvolution<T>, Derivative, nf, Mu, Obj);
  };
};

// Helper that registers the common interface of the MatchedEvolution<T>
// template for a given specialisation, using the trampoline above.
template <typename T>
py::class_<apfel::MatchedEvolution<T>, PyMatchedEvolution<T>> bind_matchedevolution(py::module_& m, char const* name, char const* doc)
{
  return py::class_<apfel::MatchedEvolution<T>, PyMatchedEvolution<T>>(m, name, doc)
         .def(py::init<T const&, double const&, std::vector<double> const&, int const&>(), "ObjRef"_a, "MuRef"_a, "Thresholds"_a, "nsteps"_a = 10)
         .def("EvolveObject", &apfel::MatchedEvolution<T>::EvolveObject, "Evolve the object with nf flavours from scale mu02 to mu2.", "nf"_a, "mu02"_a, "mu2"_a, "Obj0"_a)
         .def("MatchObject", &apfel::MatchedEvolution<T>::MatchObject, "Match the object across a threshold (Up = increasing the number of flavours).", "Up"_a, "nf"_a, "Obj"_a)
         .def("Derivative", &apfel::MatchedEvolution<T>::Derivative, "Return the right-hand side of the evolution equation at scale Mu.", "nf"_a, "Mu"_a, "Obj"_a)
         .def("Evaluate", &apfel::MatchedEvolution<T>::Evaluate, "Return the evolved object at scale mu.", "mu"_a)
         .def("GetObjectRef", &apfel::MatchedEvolution<T>::GetObjectRef, "Return the reference object.")
         .def("GetMuRef", &apfel::MatchedEvolution<T>::GetMuRef, "Return the reference scale.")
         .def("GetThresholds", &apfel::MatchedEvolution<T>::GetThresholds, "Return the matching thresholds.")
         .def("GetNumberOfSteps", &apfel::MatchedEvolution<T>::GetNumberOfSteps, "Return the number of Runge-Kutta steps.")
         .def("SetObjectRef", &apfel::MatchedEvolution<T>::SetObjectRef, "Set the reference object.", "ObjRef"_a)
         .def("SetMuRef", &apfel::MatchedEvolution<T>::SetMuRef, "Set the reference scale.", "MuRef"_a)
         .def("SetNumberOfSteps", &apfel::MatchedEvolution<T>::SetNumberOfSteps, "Set the number of Runge-Kutta steps.", "nsteps"_a);
}

PYBIND11_MODULE(apfelpy, m)
{
  // Documentation
  m.doc() = "Python wrapper of APFEL++";

  // Constants
  py::module_ _constants = m.def_submodule("constants", "Numerical constants");

  // Utility functions
  py::module_ _utilities = m.def_submodule("utilities", "Utility functions");

  // Initializers
  py::module_ _initializers = m.def_submodule("initializers", "Initialisers");

  // Builders
  py::module_ _builders = m.def_submodule("builders", "Builders");

  // betaQCD
  py::module_ _betaQCD = m.def_submodule("betaQCD", "Coefficients of the QCD beta function");

  // Gamma cusp
  py::module_ _gammaK = m.def_submodule("gammaK", "Coefficients of the QCD gamma cusp");

  // Gamma non-cusp
  py::module_ _gammaF = m.def_submodule("gammaF", "Coefficients of the QCD gamma non-cusp");

  // Collins-Soper kernel
  py::module_ _KCS = m.def_submodule("KCS", "Coefficients of the QCD Collins-Soper kernel");

  // Wrappers of "messages.h"
  m.def("SetVerbosityLevel", &apfel::SetVerbosityLevel, "Set the global verbosity level.", "vl"_a);
  m.def("GetVerbosityLevel", &apfel::GetVerbosityLevel, "Return the global verbosity level.");
  m.def("Banner",            &apfel::Banner, "Print the APFEL++ banner.");

  // Wrappers of "constants.h"
  _constants.attr("Pi2")        = apfel::Pi2;
  _constants.attr("FourPi")     = apfel::FourPi;
  _constants.attr("emc")        = apfel::emc;
  _constants.attr("zeta2")      = apfel::zeta2;
  _constants.attr("zeta3")      = apfel::zeta3;
  _constants.attr("zeta4")      = apfel::zeta4;
  _constants.attr("zeta5")      = apfel::zeta5;
  _constants.attr("zeta6")      = apfel::zeta6;
  _constants.attr("TR")         = apfel::TR;
  _constants.attr("CF")         = apfel::CF;
  _constants.attr("CA")         = apfel::CA;
  _constants.attr("NC")         = apfel::NC;
  _constants.attr("ed")         = apfel::ed;
  _constants.attr("eu")         = apfel::eu;
  _constants.attr("ed2")        = apfel::ed2;
  _constants.attr("eu2")        = apfel::eu2;
  _constants.attr("QCh")        = apfel::QCh;
  _constants.attr("QCh2")       = apfel::QCh2;
  _constants.attr("fl11ns")     = apfel::fl11ns;
  _constants.attr("fl11sg")     = apfel::fl11sg;
  _constants.attr("ConvFact")   = apfel::ConvFact;
  _constants.attr("ZMass")      = apfel::ZMass;
  _constants.attr("GammaZ")     = apfel::GammaZ;
  _constants.attr("WMass")      = apfel::WMass;
  _constants.attr("GammaW")     = apfel::GammaW;
  _constants.attr("ProtonMass") = apfel::ProtonMass;
  _constants.attr("Sin2ThetaW") = apfel::Sin2ThetaW;
  _constants.attr("GFermi")     = apfel::GFermi;
  _constants.attr("alphaem")    = apfel::alphaem;
  _constants.attr("Vud")        = apfel::Vud;
  _constants.attr("Vus")        = apfel::Vus;
  _constants.attr("Vub")        = apfel::Vub;
  _constants.attr("Vcd")        = apfel::Vcd;
  _constants.attr("Vcs")        = apfel::Vcs;
  _constants.attr("Vcb")        = apfel::Vcb;
  _constants.attr("Vtd")        = apfel::Vtd;
  _constants.attr("Vts")        = apfel::Vts;
  _constants.attr("Vtb")        = apfel::Vtb;
  _constants.attr("Vud2")       = apfel::Vud2;
  _constants.attr("Vus2")       = apfel::Vus2;
  _constants.attr("Vub2")       = apfel::Vub2;
  _constants.attr("Vcd2")       = apfel::Vcd2;
  _constants.attr("Vcs2")       = apfel::Vcs2;
  _constants.attr("Vcb2")       = apfel::Vcb2;
  _constants.attr("Vtd2")       = apfel::Vtd2;
  _constants.attr("Vts2")       = apfel::Vts2;
  _constants.attr("Vtb2")       = apfel::Vtb2;
  _constants.attr("CKM")        = apfel::CKM;
  _constants.attr("CKM2")       = apfel::CKM2;

  // Wrappers of "betaqcd.h"
  _betaQCD.def("beta0qcd", &apfel::beta0qcd, "nf"_a);
  _betaQCD.def("beta1qcd", &apfel::beta1qcd, "nf"_a);
  _betaQCD.def("beta2qcd", &apfel::beta2qcd, "nf"_a);
  _betaQCD.def("beta3qcd", &apfel::beta3qcd, "nf"_a);

  // Wrappers of "gammak.h"
  _gammaK.def("gammaK0",    &apfel::gammaK0);
  _gammaK.def("gammaK1",    &apfel::gammaK1,    "nf"_a);
  _gammaK.def("gammaK2",    &apfel::gammaK2,    "nf"_a);
  _gammaK.def("gammaK3",    &apfel::gammaK3,    "nf"_a);
  _gammaK.def("gammaK3gmq", &apfel::gammaK3gmq, "nf"_a);

  // Wrappers of "gammaf.h"
  _gammaF.def("gammaFq0", &apfel::gammaFq0);
  _gammaF.def("gammaFq1", &apfel::gammaFq1, "nf"_a);
  _gammaF.def("gammaFq2", &apfel::gammaFq2, "nf"_a);
  _gammaF.def("gammaFg0", &apfel::gammaFg0, "nf"_a);
  _gammaF.def("gammaFg1", &apfel::gammaFg1, "nf"_a);
  _gammaF.def("gammaFg2", &apfel::gammaFg2, "nf"_a);

  // Wrappers of "kcs.h"
  _KCS.def("KCS00", &apfel::KCS00);
  _KCS.def("KCS01", &apfel::KCS01);
  _KCS.def("KCS10", &apfel::KCS10, "nf"_a);
  _KCS.def("KCS11", &apfel::KCS11, "nf"_a);
  _KCS.def("KCS12", &apfel::KCS12, "nf"_a);
  _KCS.def("KCS20", &apfel::KCS20, "nf"_a);
  _KCS.def("KCS21", &apfel::KCS21, "nf"_a);
  _KCS.def("KCS22", &apfel::KCS22, "nf"_a);
  _KCS.def("KCS23", &apfel::KCS23, "nf"_a);

  // Wrappers of "lhtoypdfs.h"
  _utilities.def("LHToyPDFs", &apfel::LHToyPDFs, "x"_a, "Q"_a);
  _utilities.def("LHToyPDFsPhys", &apfel::LHToyPDFsPhys, "x"_a, "Q"_a);
  _utilities.def("LHToyPDFsPol", &apfel::LHToyPDFsPol, "x"_a, "Q"_a);
  _utilities.def("LHToyFFs", &apfel::LHToyFFs, "x"_a, "Q"_a);

  // Wrappers of "tools.h"
  py::enum_<apfel::QuarkFlavour>(_utilities, "QuarkFlavour")
  .value("TOTAL",   apfel::QuarkFlavour::TOTAL)
  .value("DOWN",    apfel::QuarkFlavour::DOWN)
  .value("UP",      apfel::QuarkFlavour::UP)
  .value("STRANGE", apfel::QuarkFlavour::STRANGE)
  .value("CHARM",   apfel::QuarkFlavour::CHARM)
  .value("BOTTOM",  apfel::QuarkFlavour::BOTTOM)
  .value("TOP",     apfel::QuarkFlavour::TOP);
  _utilities.def("NF", &apfel::NF, "Q"_a, "Thresholds"_a);
  _utilities.def("ElectroWeakCharges", &apfel::ElectroWeakCharges, "Q"_a, "virt"_a, "Comp"_a = apfel::QuarkFlavour::TOTAL);
  _utilities.def("ParityViolatingElectroWeakCharges", &apfel::ParityViolatingElectroWeakCharges, "Q"_a, "virt"_a, "Comp"_a = apfel::QuarkFlavour::TOTAL);
  _utilities.def("ElectroWeakChargesNWA", &apfel::ElectroWeakChargesNWA);
  _utilities.def("ProductExpansion", &apfel::ProductExpansion, "r"_a);
  _utilities.def("GetSIATotalCrossSection", &apfel::GetSIATotalCrossSection,
                 "PerturbativeOrder"_a,
                 "Q"_a,
                 "AlphaQCD"_a,
                 "AlphaQED"_a,
                 "Thresholds"_a,
                 "Comp"_a = apfel::QuarkFlavour::TOTAL,
                 "NoCharges"_a = false);

  // Wrappers of "rotations.h"
  _utilities.def("PhysToQCDEv", py::overload_cast<std::map<int, double> const&>(&apfel::PhysToQCDEv),              "InPhysMap"_a);
  //_utilities.def("PhysToQCDEv", py::overload_cast<apfel::Set<apfel::Distribution> const&>(&apfel::PhysToQCDEv),    "InPhysMap"_a, "nf"_a);
  //_utilities.def("PhysToQCDEv", py::overload_cast<apfel::Set<apfel::Operator> const&>(&apfel::PhysToQCDEv),        "InPhysMap"_a, "nf"_a);
  _utilities.def("QCDEvToPhys", py::overload_cast<std::map<int, double> const&>(&apfel::QCDEvToPhys),              "QCDEvMap"_a);
  _utilities.def("QCDEvToPhys", py::overload_cast<std::map<int, apfel::Distribution> const&>(&apfel::QCDEvToPhys), "QCDEvMap"_a);
  _utilities.def("QCDEvToPhys", py::overload_cast<std::map<int, apfel::Operator> const&>(&apfel::QCDEvToPhys),     "QCDEvMap"_a);

  // Wrappers of "timer.h"
  py::class_<apfel::Timer>(m, "Timer", "Computes the time elapsed between start and stop.")
  .def(py::init<>())
  .def("start", &apfel::Timer::start, "Start the timer.")
  .def("stop",  &apfel::Timer::stop, "Stop the timer and report the elapsed time.", "ForceDisplay"_a = false);

  // Wrappers of "subgrid.h"
  py::class_<apfel::SubGrid>(m, "SubGrid", "Class for the x-space interpolation subgrids.")
  .def(py::init<int const&, double const&, int const&>(), "Construct a logarithmically-spaced subgrid from the number of nodes, the lower bound and the interpolation degree.", "nx"_a, "xMin"_a, "InterDegree"_a)
  .def(py::init<std::vector<double> const&, int const&>(), "Construct a subgrid from an explicit list of nodes and the interpolation degree.", "xsg"_a, "InterDegree"_a)
  .def("nx", &apfel::SubGrid::nx, "Return the number of x points.")
  .def("InterDegree", &apfel::SubGrid::InterDegree, "Return the interpolation degree.")
  .def("xMin", &apfel::SubGrid::xMin, "Return the minimum node value.")
  .def("xMax", &apfel::SubGrid::xMax, "Return the maximum node value.")
  .def("Step", &apfel::SubGrid::Step, "Return the step size of the log grid.")
  .def("GetGrid", &apfel::SubGrid::GetGrid, "Return the grid.")
  .def("GetLogGrid", &apfel::SubGrid::GetLogGrid, "Return the log-grid.")
  .def("Print", &apfel::SubGrid::Print, "Print the SubGrid object.")
  .def(py::self == py::self)
  .def(py::self != py::self);

  // Wrappers of "grid.h"
  py::class_<apfel::Grid>(m, "Grid", "Collection of SubGrids defining the full x-space interpolation grid.")
  .def(py::init<std::vector<apfel::SubGrid> const&>(), "Construct a Grid from a vector of subgrids.", "grs"_a)
  .def("nGrids", &apfel::Grid::nGrids, "Return the number of subgrids.")
  .def("SubToJointMap", &apfel::Grid::SubToJointMap, "Return the map of indices from the joint grid to the subgrids.")
  .def("JointToSubMap", &apfel::Grid::JointToSubMap, "Return the map of indices from the subgrids to the joint grid.")
  .def("TransitionPoints", &apfel::Grid::TransitionPoints, "Return the vector of transition indices on the joint grid.")
  .def("GetSubGrids", &apfel::Grid::GetSubGrids, "Return the vector of subgrids.")
  .def("GetSubGrid", &apfel::Grid::GetSubGrid, "Return the ig-th subgrid.", "ig"_a)
  .def("GetJointGrid", &apfel::Grid::GetJointGrid, "Return the joint subgrid.")
  .def("Print", &apfel::Grid::Print, "Print the Grid object.")
  .def(py::self == py::self)
  .def(py::self != py::self);

  // Wrappers of "interpolator.h"
  // Trampoline class for virtual class
  class PyInterpolator: public apfel::Interpolator
  {
  public:
    using Interpolator::Interpolator;
    double InterpolantLog(int const& beta, double const& lnx, apfel::SubGrid const& sg) const override
    {
      PYBIND11_OVERRIDE_PURE(double, Interpolator, InterpolantLog, beta, lnx, sg);
    };
    double Interpolant(int const& beta, double const& x, apfel::SubGrid const& sg) const override
    {
      PYBIND11_OVERRIDE_PURE(double, Interpolator, Interpolant, beta, x, sg);
    };
    std::array<int, 2> SumBounds(double const& x, apfel::SubGrid const& sg) const override
    {
      PYBIND11_OVERRIDE_PURE(PYBIND11_TYPE(std::array<int, 2>), Interpolator, SumBounds, x, sg);
    };
  };
  py::class_<apfel::Interpolator, PyInterpolator>(m, "Interpolator", "Mother class for the x-space interpolation.")
  .def(py::init<apfel::Grid const&>(), "Construct an interpolator on the given grid.", "gr"_a)
  .def(py::init<apfel::Grid const&, std::vector<std::vector<double>> const&, std::vector<double> const&>(), "Construct an interpolator on the given grid from explicit subgrid and joint-grid values.", "gr"_a, "distsubgrid"_a, "distjointgrid"_a)
  .def("Evaluate", py::overload_cast<double const&>(&apfel::Interpolator::Evaluate, py::const_), "Evaluate the interpolated function on the joint grid.", "x"_a)
  .def("Evaluate", py::overload_cast<double const&, int const&>(&apfel::Interpolator::Evaluate, py::const_), "Evaluate the interpolated function on a given subgrid.", "x"_a, "ig"_a)
  .def("Derive", &apfel::Interpolator::Derive, "Evaluate the derivative of the interpolated function.", "x"_a)
  .def("Integrate", &apfel::Interpolator::Integrate, "Evaluate the integral of the interpolated function between a and b.", "a"_a, "b"_a)
  .def("Squash", &apfel::Interpolator::Squash, "Sum all entries of the distribution and return the squashed distribution.")
  .def("InterpolantLog", &apfel::Interpolator::InterpolantLog, "Interpolating function in log space (interpolation weights).", "beta"_a, "lnx"_a, "_asg"_a)
  .def("Interpolant", &apfel::Interpolator::Interpolant, "Interpolating function (interpolation weights).", "beta"_a, "x"_a, "sg"_a)
  .def("DerInterpolant", &apfel::Interpolator::DerInterpolant, "Derivative of the interpolating function.")
  .def("IntInterpolant", &apfel::Interpolator::IntInterpolant, "Integral of the interpolating function.")
  .def("SumBounds", &apfel::Interpolator::SumBounds, "Lower and upper bounds of the grid index entering the interpolation sum.")
  .def("GetGrid", &apfel::Interpolator::GetGrid, "Return the Grid object.")
  .def("GetDistributionSubGrid", &apfel::Interpolator::GetDistributionSubGrid, "Return the distribution values on the subgrids.")
  .def("GetDistributionJointGrid", &apfel::Interpolator::GetDistributionJointGrid, "Return the distribution values on the joint grid.")
  .def("Print", &apfel::Interpolator::Print, "Print the Interpolator object.");

  // Wrappers of "lagrangeinterpolator.h"
  py::class_<apfel::LagrangeInterpolator, apfel::Interpolator>(m, "LagrangeInterpolator", "Specialisation of the Interpolator class using Lagrange interpolation.")
  .def(py::init<apfel::Grid const&>(), "Construct a Lagrange interpolator on the given grid.", "gr"_a)
  .def(py::init<apfel::Grid const&, std::vector<std::vector<double>> const&, std::vector<double> const&>(), "Construct a Lagrange interpolator on the given grid from explicit subgrid and joint-grid values.", "gr"_a, "distsubgrid"_a, "distjointgrid"_a)
  .def("InterpolantLog", &apfel::LagrangeInterpolator::InterpolantLog, "Interpolating function in log space (interpolation weights).", "beta"_a, "lnx"_a, "_asg"_a)
  .def("Interpolant", &apfel::LagrangeInterpolator::Interpolant, "Interpolating function (interpolation weights).", "beta"_a, "x"_a, "sg"_a)
  .def("DerInterpolant", &apfel::LagrangeInterpolator::DerInterpolant, "Derivative of the interpolating function.", "beta"_a, "x"_a, "sg"_a)
  .def("IntInterpolant", &apfel::LagrangeInterpolator::IntInterpolant, "Integral of the interpolating function.", "beta"_a, "a"_a, "b"_a, "sg"_a)
  .def("SumBounds", &apfel::LagrangeInterpolator::SumBounds, "Lower and upper bounds of the grid index entering the interpolation sum.", "x"_a, "sg"_a);

  // Wrappers of "distribution.h"
  py::class_<apfel::Distribution, apfel::LagrangeInterpolator>(m, "Distribution", "One of the basic objects of the library: a function discretised on an x-space grid.")
  .def(py::init<apfel::Grid const&>(), "g"_a)
  .def(py::init<apfel::Distribution const&>(), "obj"_a)
  .def(py::init<apfel::Distribution const&, std::vector<std::vector<double>> const&, std::vector<double> const&>(), "obj"_a, "distsubgrid"_a, "distjointgrid"_a)
  .def(py::init<apfel::Grid const&, std::vector<std::vector<double>> const&, std::vector<double> const&>(), "g"_a, "distsubgrid"_a, "distjointgrid"_a)
  .def(py::init<apfel::Grid const&, std::function<double(double const&)>>(), "g"_a, "InDistFunc"_a)
  .def(py::init<apfel::Grid const&, std::function<double(double const&, double const&)>, double const&>(), "g"_a, "InDistFunc"_a, "Q"_a)
  .def(py::init<apfel::Grid const&, std::function<double(int const&, double const&)>, int const&>(), "g"_a, "InDistFunc"_a, "ipdf"_a)
  .def(py::init<apfel::Grid const&, std::function<double(int const&, double const&, double const&)>, int const&, double const&>(), "g"_a, "InDistFunc"_a, "ipdf"_a, "Q"_a)
  .def("SetJointGrid", static_cast<void (apfel::Distribution::*)(int const&, double const&)>(&apfel::Distribution::SetJointGrid), "Set the value of the joint grid at node ix.", "ix"_a, "x"_a)
  .def("SetSubGrid", static_cast<void (apfel::Distribution::*)(int const&, int const&, double const&)>(&apfel::Distribution::SetSubGrid), "Set the value of subgrid ig at node ix.", "ig"_a, "ix"_a, "x"_a)
  .def("SetJointGrid", static_cast<void (apfel::Distribution::*)(std::vector<double> const&)>(&apfel::Distribution::SetJointGrid), "Set all the joint-grid values at once.", "jg"_a)
  .def("SetSubGrid", static_cast<void (apfel::Distribution::*)(int const&, std::vector<double> const&)>(&apfel::Distribution::SetSubGrid), "Set all the values of subgrid ig at once.", "ig"_a, "sg"_a)
  .def("Derivative", &apfel::Distribution::Derivative, "Return the derivative of the Distribution as a new Distribution.")
  //.def(py::self = py::self)  // DOES NOT WORK!
  .def(py::self *= double())
  .def(py::self *= std::function<double(double const&)>())
  .def(py::self /= double())
  .def(py::self *= py::self)
  .def(py::self += py::self)
  //.def(py::self -= py::self)
  .def(py::self * double())
  .def(double() * py::self)
  .def(py::self * std::function<double(double const&)>())
  .def(std::function<double(double const&)>() * py::self)
  .def(py::self / double())
  .def(py::self + py::self)
  .def(py::self - py::self)
  .def(py::self * py::self);

  // Wrappers of "expression.h"
  // Trampoline class for virtual class
  class PyExpression: public apfel::Expression
  {
  public:
    using Expression::Expression;
    double Regular(double const& x) const override
    {
      PYBIND11_OVERRIDE(double, Expression, Regular, x);
    };
    double Singular(double const& x) const override
    {
      PYBIND11_OVERRIDE(double, Expression, Singular, x);
    };
    double Local(double const& x) const override
    {
      PYBIND11_OVERRIDE(double, Expression, Local, x);
    };
    double LocalPP(double const& x) const override
    {
      PYBIND11_OVERRIDE(double, Expression, LocalPP, x);
    };
    double SingularPV(double const& x) const override
    {
      PYBIND11_OVERRIDE(double, Expression, SingularPV, x);
    };
    double LocalPV() const override
    {
      PYBIND11_OVERRIDE(double, Expression, LocalPV);
    };
    double LocalLogPV(double const& x) const override
    {
      PYBIND11_OVERRIDE(double, Expression, LocalLogPV, x);
    };
  };
  py::class_<apfel::Expression, PyExpression>(m, "Expression", "Encapsulates a convolution kernel in terms of its regular, singular and local parts.")
  .def(py::init<double const&>(), "eta"_a = 1)
  .def("Regular", &apfel::Expression::Regular, "Regular term at x.", "x"_a)
  .def("Singular", &apfel::Expression::Singular, "Singular term at x.", "x"_a)
  .def("Local", &apfel::Expression::Local, "Local term at x.", "x"_a)
  .def("LocalPP", &apfel::Expression::LocalPP, "Local term for ++-prescribed distributions at x.", "x"_a)
  .def("SingularPV", &apfel::Expression::SingularPV, "Singular term for principal-valued distributions at x.", "x"_a)
  .def("LocalPV", &apfel::Expression::LocalPV, "Local term for principal-valued distributions.")
  .def("LocalLogPV", &apfel::Expression::LocalLogPV, "Log-dependent local term for principal-valued distributions at x.", "x"_a)
  .def("SetExternalVariable", &apfel::Expression::SetExternalVariable, "Set the value of a possible external variable.", "extvar"_a)
  .def("eta", &apfel::Expression::eta, "Return the value of the scaling parameter eta.");

  py::class_<apfel::Identity, apfel::Expression>(m, "Identity", "Expression implementing the identity operator.")
  .def(py::init<>())
  .def("Local", &apfel::Expression::Local, "Local term at x.", "x"_a);

  py::class_<apfel::Null, apfel::Expression>(m, "Null", "Expression implementing the null operator.")
  .def(py::init<>());

  // Wrappers of "matrix.h" (this is a template class and needs a
  // wrapper for any specialisation)
  py::class_<apfel::matrix<double>>(m, "matrixDouble", "Row-major matrix container with double entries.")
  .def(py::init<size_t const&, size_t const&>(), "row"_a = 0, "col"_a = 0)
  .def("resize", &apfel::matrix<double>::resize, "row"_a, "col"_a, "v"_a = 0)
  .def("set", &apfel::matrix<double>::set, "v"_a)
  .def("size", py::overload_cast<>(&apfel::matrix<double>::size, py::const_))
  .def("size", py::overload_cast<size_t const&>(&apfel::matrix<double>::size, py::const_), "dim"_a)
  .def("__call__", [] (apfel::matrix<double>& c, size_t const& i, size_t const& j)
  {
    return c(i, j);
  })
  .def("__call__", [] (apfel::matrix<double> const& c, size_t const& i, size_t const& j)
  {
    return c(i, j);
  });

  py::class_<apfel::matrix<std::vector<int>>>(m, "matrixVectorInt", "Row-major matrix container whose entries are vectors of int.")
  .def(py::init<size_t const&, size_t const&>(), "row"_a = 0, "col"_a = 0)
  .def("resize", &apfel::matrix<std::vector<int>>::resize, "row"_a, "col"_a, "v"_a = 0)
  .def("set", &apfel::matrix<std::vector<int>>::set, "v"_a)
  .def("size", py::overload_cast<>(&apfel::matrix<std::vector<int>>::size, py::const_))
  .def("size", py::overload_cast<size_t const&>(&apfel::matrix<std::vector<int>>::size, py::const_), "dim"_a)
  .def("__call__", [] (apfel::matrix<std::vector<int>>& c, size_t const& i, size_t const& j)
  {
    return c(i, j);
  })
  .def("__call__", [] (apfel::matrix<std::vector<int>> const& c, size_t const& i, size_t const& j)
  {
    return c(i, j);
  });

  py::class_<apfel::matrix<std::vector<double>>>(m, "matrixVectorDouble", "Row-major matrix container whose entries are vectors of double.")
  .def(py::init<size_t const&, size_t const&>(), "row"_a = 0, "col"_a = 0)
  .def("resize", &apfel::matrix<std::vector<double>>::resize, "row"_a, "col"_a, "v"_a = 0)
  .def("set", &apfel::matrix<std::vector<double>>::set, "v"_a)
  .def("size", py::overload_cast<>(&apfel::matrix<std::vector<double>>::size, py::const_))
  .def("size", py::overload_cast<size_t const&>(&apfel::matrix<std::vector<double>>::size, py::const_), "dim"_a)
  .def("__call__", [] (apfel::matrix<std::vector<double>>& c, size_t const& i, size_t const& j)
  {
    return c(i, j);
  })
  .def("__call__", [] (apfel::matrix<std::vector<double>> const& c, size_t const& i, size_t const& j)
  {
    return c(i, j);
  });

  // Wrappers of "operator.h"
  py::class_<apfel::Operator>(m, "Operator", "Discretised convolution kernel (operator) acting on distributions over the grid.")
  .def(py::init<apfel::Operator const&>(), "g"_a)
  //.def(py::init<apfel::Grid const&>(), "g"_a)
  .def(py::init<apfel::Grid const&, apfel::Expression const&, double const&>(), "Construct an operator on the grid from an expression, with integration accuracy eps.", "g"_a, "expr"_a, "eps"_a = 1e-5)
  .def("Evaluate", &apfel::Operator::Evaluate, "Interpolate the operator over its first variable at x.", "x"_a)
  .def("GetGrid", &apfel::Operator::GetGrid, "Return the Grid object associated with the operator.")
  .def("GetOperator", &apfel::Operator::GetOperator, "Return the operator container (the raw operator values).")
  .def("Print", &apfel::Operator::Print, "Print the Operator object.")
  .def(py::self *= py::self)
  //.def(py::self = py::self)  // DOES NOT WORK!
  .def(py::self *= double())
  .def(py::self *= std::function<double(double const&)>())
  .def(py::self /= double())
  .def(py::self += py::self)
  //.def(py::self -= py::self)
  .def(py::self * apfel::Distribution(apfel::Grid{{apfel::SubGrid{10, 1e-5, 3}}}))
  .def(py::self * py::self)
  .def(py::self * double())
  .def(double() * py::self)
  .def(py::self * std::function<double(double const&)>())
  .def(std::function<double(double const&)>() * py::self)
  .def(py::self / double())
  .def(py::self + py::self)
  .def(py::self - py::self);

  // Wrappers of "integrator.h"
  py::enum_<apfel::Integrator::IntegrationMethod>(m, "IntegrationMethod", "Numerical integration method.")
  .value("GAUSS_LEGENDRE", apfel::Integrator::IntegrationMethod::GAUSS_LEGENDRE)
  .value("GAUSS_KRONROD",  apfel::Integrator::IntegrationMethod::GAUSS_KRONROD);
  py::class_<apfel::Integrator>(m, "Integrator", "Performs one-dimensional numerical integration of a function.")
  .def(py::init<std::function<double(double const&)> const&, apfel::Integrator::IntegrationMethod const&>(), "Construct an integrator from the integrand and the integration method.", "func"_a, "method"_a = apfel::Integrator::IntegrationMethod::GAUSS_KRONROD)
  .def("integrate", py::overload_cast<double const&, double const&, double const&>(&apfel::Integrator::integrate, py::const_), "Adaptive integration between xmin and xmax to relative accuracy eps.", "xmin"_a, "xmax"_a, "eps"_a)
  .def("integrate", py::overload_cast<double const&, double const&, std::vector<double> const&, double const&>(&apfel::Integrator::integrate, py::const_), "Adaptive integration between xmin and xmax with fixed points, to relative accuracy eps.", "xmin"_a, "xmax"_a, "FixPts"_a, "eps"_a)
  .def("integrate", py::overload_cast<std::vector<double> const&, double const&>(&apfel::Integrator::integrate, py::const_), "Adaptive integration over the interval spanned by the fixed points, to relative accuracy eps.", "FixPts"_a, "eps"_a)
  .def("integrate", py::overload_cast<double const&, double const&, int const&>(&apfel::Integrator::integrate, py::const_), "Fixed-order integration between xmin and xmax using n subintervals.", "xmin"_a, "xmax"_a, "n"_a = 1)
  .def("integrate", py::overload_cast<double const&, double const&, std::vector<double> const&, int const&>(&apfel::Integrator::integrate, py::const_), "Fixed-order integration between xmin and xmax with fixed points, using n subintervals.", "xmin"_a, "xmax"_a, "FixPts"_a, "n"_a = 1)
  .def("integrate", py::overload_cast<std::vector<double> const&, int const&>(&apfel::Integrator::integrate, py::const_), "Fixed-order integration over the interval spanned by the fixed points, using n subintervals.", "FixPts"_a, "n"_a = 1)
  .def("integrand", &apfel::Integrator::integrand, "Return the integrand evaluated at x.", "x"_a)
  .def("Method", &apfel::Integrator::Method, "Return the integration method.");

  // Wrappers of "integrator2d.h"
  py::class_<apfel::Integrator2D>(m, "Integrator2D", "Performs two-dimensional numerical integration of a function.")
  .def(py::init<std::function<double(double const&, double const&)> const&, apfel::Integrator::IntegrationMethod const&>(), "Construct a 2D integrator from the integrand and the integration method.", "func"_a, "method"_a = apfel::Integrator::IntegrationMethod::GAUSS_KRONROD)
  .def("integrate", &apfel::Integrator2D::integrate, "Integrate over the rectangle [xmin, xmax] x [ymin, ymax] to relative accuracy eps.", "xmin"_a, "xmax"_a, "ymin"_a, "ymax"_a, "eps"_a)
  .def("integrand", &apfel::Integrator2D::integrand, "Return the integrand evaluated at (x, y).", "x"_a, "y"_a);

  // Wrappers of "doubleobject.h"
  py::class_<apfel::term<apfel::Distribution>>(m, "termD", "A single term of a DoubleObject (a coefficient and a pair of single objects).");
  py::class_<apfel::term<apfel::Operator>>(m, "termO", "A single term of a DoubleObject (a coefficient and a pair of single objects).");
  py::class_<apfel::term<apfel::Distribution, apfel::Operator>>(m, "termDO", "A single term of a DoubleObject (a coefficient and a pair of single objects).");
  py::class_<apfel::term<apfel::Operator, apfel::Distribution>>(m, "termOD", "A single term of a DoubleObject (a coefficient and a pair of single objects).");

  py::class_<apfel::DoubleObject<apfel::Distribution>>(m, "DoubleObjectD", "Collection of pairs of single objects representing a double (two-variable) object.")
  .def(py::init<>())
  .def(py::init<std::vector<apfel::term<apfel::Distribution>>>(), "terms"_a)
  .def("AddTerm", &apfel::DoubleObject<apfel::Distribution>::AddTerm, "Add a term to the double object.", "newterm"_a)
  .def("GetTerms", &apfel::DoubleObject<apfel::Distribution>::GetTerms, "Return the list of terms.")
  .def("Evaluate", py::overload_cast<double const&, double const&>(&apfel::DoubleObject<apfel::Distribution>::Evaluate, py::const_), "Evaluate the double object at (x, z).", "x"_a, "z"_a)
  .def("Evaluate", py::overload_cast<int const&, double const&, double const&>(&apfel::DoubleObject<apfel::Distribution>::Evaluate, py::const_), "Evaluate the i-th term at (x, z).", "i"_a, "x"_a, "z"_a)
  .def("Evaluate1", &apfel::DoubleObject<apfel::Distribution>::Evaluate1, "Evaluate the first single object at x.", "x"_a)
  .def("Evaluate2", &apfel::DoubleObject<apfel::Distribution>::Evaluate2, "Evaluate the second single object at z.", "z"_a)
  .def("Derive", &apfel::DoubleObject<apfel::Distribution>::Derive, "Evaluate the derivative at (x, z).", "x"_a, "z"_a)
  .def("Derive1", &apfel::DoubleObject<apfel::Distribution>::Derive1, "Evaluate the derivative of the first single object at x.", "x"_a)
  .def("Derive2", &apfel::DoubleObject<apfel::Distribution>::Derive2, "Evaluate the derivative of the second single object at z.", "z"_a)
  .def("Integrate", py::overload_cast<double const&, double const&, double const&, double const&>(&apfel::DoubleObject<apfel::Distribution>::Integrate, py::const_), "Integrate over the rectangle [xl, xu] x [zl, zu].", "xl"_a, "xu"_a, "zl"_a, "zu"_a)
  .def("Integrate1", &apfel::DoubleObject<apfel::Distribution>::Integrate1, "Integrate the first single object between xl and xu.", "xl"_a, "xu"_a)
  .def("Integrate2", &apfel::DoubleObject<apfel::Distribution>::Integrate2, "Integrate the second single object between zl and zu.", "zl"_a, "zu"_a)
  .def("Integrate", py::overload_cast<double const&, double const&, std::function<double(double const&)>, std::function<double(double const&)>>(&apfel::DoubleObject<apfel::Distribution>::Integrate, py::const_), "Integrate with z-bounds given as functions of x.", "xl"_a, "xu"_a, "zlx"_a, "zux"_a)
  .def("Integrate", py::overload_cast<std::function<double(double const&)>, std::function<double(double const&)>, double const&, double const&>(&apfel::DoubleObject<apfel::Distribution>::Integrate, py::const_), "Integrate with x-bounds given as functions of z.", "xlz"_a, "xuz"_a, "zl"_a, "zu"_a)
  .def("MultiplyBy", &apfel::DoubleObject<apfel::Distribution>::MultiplyBy, "Multiply the single objects by functions of their respective variables.", "fx"_a, "fz"_a)
  .def("Print", &apfel::DoubleObject<apfel::Distribution>::Print, "Print the DoubleObject.")
  .def(py::self *= double())
  //.def(py::self *= py::self)
  .def(py::self *= std::function<double(double const&)>())
  .def(py::self /= double())
  .def(py::self += py::self)
  //.def(py::self -= py::self)
  .def(py::self * double())
  .def(double() * py::self)
  .def(py::self / double())
  //.def(py::self * py::self)
  .def(py::self + py::self)
  .def(py::self - py::self);

  py::class_<apfel::DoubleObject<apfel::Operator>>(m, "DoubleObjectO", "Collection of pairs of single objects representing a double (two-variable) object.")
  .def(py::init<>())
  .def(py::init<std::vector<apfel::term<apfel::Operator>>>(), "terms"_a)
  .def("AddTerm", &apfel::DoubleObject<apfel::Operator>::AddTerm, "newterm"_a)
  .def("GetTerms", &apfel::DoubleObject<apfel::Operator>::GetTerms)
  .def("Print", &apfel::DoubleObject<apfel::Operator>::Print)
  .def(py::self *= double())
  //.def(py::self *= py::self)
  .def(py::self *= std::function<double(double const&)>())
  .def(py::self /= double())
  .def(py::self += py::self)
  //.def(py::self -= py::self)
  .def(py::self * double())
  .def(double() * py::self)
  .def(py::self / double())
  //.def(py::self * py::self)
  .def(py::self + py::self)
  .def(py::self - py::self);

  py::class_<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>(m, "DoubleObjectDO", "Collection of pairs of single objects representing a double (two-variable) object.")
  .def(py::init<>())
  .def(py::init<std::vector<apfel::term<apfel::Distribution, apfel::Operator>>>(), "terms"_a)
  .def("AddTerm", &apfel::DoubleObject<apfel::Distribution, apfel::Operator>::AddTerm, "newterm"_a)
  .def("GetTerms", &apfel::DoubleObject<apfel::Distribution, apfel::Operator>::GetTerms)
  .def("Print", &apfel::DoubleObject<apfel::Distribution, apfel::Operator>::Print)
  .def(py::self *= double())
  //.def(py::self *= py::self)
  .def(py::self *= std::function<double(double const&)>())
  .def(py::self /= double())
  .def(py::self += py::self)
  //.def(py::self -= py::self)
  .def(py::self * double())
  .def(double() * py::self)
  .def(py::self / double())
  //.def(py::self * py::self)
  .def(py::self + py::self)
  .def(py::self - py::self);

  py::class_<apfel::DoubleObject<apfel::Operator, apfel::Distribution>>(m, "DoubleObjectOD", "Collection of pairs of single objects representing a double (two-variable) object.")
  .def(py::init<>())
  .def(py::init<std::vector<apfel::term<apfel::Operator, apfel::Distribution>>>(), "terms"_a)
  .def("AddTerm", &apfel::DoubleObject<apfel::Operator, apfel::Distribution>::AddTerm, "newterm"_a)
  .def("GetTerms", &apfel::DoubleObject<apfel::Operator, apfel::Distribution>::GetTerms)
  .def("Print", &apfel::DoubleObject<apfel::Operator, apfel::Distribution>::Print)
  .def(py::self *= double())
  //.def(py::self *= py::self)
  .def(py::self *= std::function<double(double const&)>())
  .def(py::self /= double())
  .def(py::self += py::self)
  //.def(py::self -= py::self)
  .def(py::self * double())
  .def(double() * py::self)
  .def(py::self / double())
  .def(py::self + py::self)
  .def(py::self - py::self);

  // Wrappers of "convolutionmap.h"
  py::class_<apfel::ConvolutionMap::rule>(m, "rule", "A single rule of a convolution map (operand, object and coefficient).")
  .def_readwrite("operand", &apfel::ConvolutionMap::rule::operand)
  .def_readwrite("object", &apfel::ConvolutionMap::rule::object)
  .def_readwrite("coefficient", &apfel::ConvolutionMap::rule::coefficient)
  .def(py::self == py::self);

  py::class_<apfel::ConvolutionMap>(m, "ConvolutionMap", "Encapsulates the rules describing how operators act on objects in a convolution.")
  .def(py::init<std::string const&>(), "Construct a (named) convolution map.", "name"_a)
  .def("SetRules", &apfel::ConvolutionMap::SetRules, "Set the rules of the convolution map.", "rules"_a)
  .def("GetName", &apfel::ConvolutionMap::GetName, "Return the name of the map.")
  .def("GetRules", &apfel::ConvolutionMap::GetRules, "Return the full set of rules.")
  .def("GetRuleMatrix", &apfel::ConvolutionMap::GetRuleMatrix, "Return the rule matrix of coefficients.")
  .def("GetRuleIndices", &apfel::ConvolutionMap::GetRuleIndices, "Return the operand indices of the rules.")
  .def("Print", &apfel::ConvolutionMap::Print, "Print the ConvolutionMap object.");

  py::class_<apfel::DiagonalBasis, apfel::ConvolutionMap>(m, "DiagonalBasis", "Diagonal convolution map (each object convolved with a single operator).")
  .def(py::init<int const&, int const&>(), "nf"_a, "offset"_a = 0);

  // Wrappers of "disbasis.h"
  py::class_<apfel::DISNCBasis, apfel::ConvolutionMap>(m, "DISNCBasis", "Convolution map for neutral-current DIS structure functions.")
  .def(py::init<int const&, double const&>(), "k"_a, "fact"_a = 1)
  .def(py::init<std::vector<double> const&>(), "Ch"_a);

  py::class_<apfel::DISCCBasis, apfel::ConvolutionMap>(m, "DISCCBasis", "Convolution map for charged-current DIS structure functions.")
  .def(py::init<int const&, bool const&, double const&>(), "l"_a, "Is3"_a, "fact"_a = 1)
  .def(py::init<std::vector<double> const&, bool const&>(), "CKM"_a, "Is3"_a);

  // Wrappers of "evolutionbasis.h"
  py::class_<apfel::EvolutionBasisQCD, apfel::ConvolutionMap>(m, "EvolutionBasisQCD", "QCD evolution basis: convolution map for the DGLAP evolution of distributions.")
  .def(py::init<int const&>(), "nf"_a);

  py::class_<apfel::EvolutionOperatorBasisQCD, apfel::ConvolutionMap>(m, "EvolutionOperatorBasisQCD", "QCD evolution basis: convolution map for the DGLAP evolution of operators.")
  .def(py::init<int const&>(), "nf"_a);

  py::class_<apfel::EvolveDistributionsBasisQCD, apfel::ConvolutionMap>(m, "EvolveDistributionsBasisQCD", "Convolution map to evolve a full set of distributions in the QCD basis.")
  .def(py::init<>());

  // Wrappers of "matchingbasisqcd.h"
  py::class_<apfel::MatchingBasisQCD, apfel::ConvolutionMap>(m, "MatchingBasisQCD", "Convolution map for the QCD matching of distributions at heavy-quark thresholds.")
  .def(py::init<int const&>(), "nf"_a);

  py::class_<apfel::MatchingOperatorBasisQCD, apfel::ConvolutionMap>(m, "MatchingOperatorBasisQCD", "Convolution map for the QCD matching of operators at heavy-quark thresholds.")
  .def(py::init<int const&>(), "nf"_a);

  // Wrappers of "set.h"
  py::class_<apfel::Set<apfel::Distribution>>(m, "SetD", "Collection of objects (distributions) together with a convolution map.")
  .def(py::init<apfel::ConvolutionMap const&, std::map<int, apfel::Distribution> const&>(), "Map"_a = apfel::ConvolutionMap{"UNDEFINED"}, "in"_a = std::map<int, apfel::Distribution> {})
  .def(py::init<std::map<int, apfel::Distribution> const&>(), "in"_a)
  .def("at", &apfel::Set<apfel::Distribution>::at, "Return the object with the given ID.", "id"_a)
  .def("GetMap", &apfel::Set<apfel::Distribution>::GetMap, "Return the convolution map.")
  .def("GetObjects", &apfel::Set<apfel::Distribution>::GetObjects, "Return the full map of objects.")
  .def("SetMap", &apfel::Set<apfel::Distribution>::SetMap, "(Re)set the convolution map.", "map"_a)
  .def("SetObjects", &apfel::Set<apfel::Distribution>::SetObjects, "(Re)set the map of objects.", "objects"_a)
  .def("Combine", py::overload_cast<>(&apfel::Set<apfel::Distribution>::Combine, py::const_), "Sum all objects of the set into a single one.")
  .def("Combine", py::overload_cast<std::vector<double> const&>(&apfel::Set<apfel::Distribution>::Combine, py::const_), "Sum all objects of the set weighted by the given coefficients.", "weights"_a)
  .def("Print", &apfel::Set<apfel::Distribution>::Print, "Print the Set object.")
  .def(py::self *= double())
  //.def(py::self *= py::self)
  .def(py::self *= std::function<double(double const&)>())
  .def(py::self *= std::vector<double>())
  .def(py::self *= std::map<int, double>())
  .def(py::self /= double())
  .def(py::self += py::self)
  //.def(py::self -= py::self)
  .def(py::self * double())
  .def(double() * py::self)
  .def(py::self * std::function<double(double const&)>())
  .def(std::function<double(double const&)>() * py::self)
  .def(py::self * std::vector<double>())
  .def(std::vector<double>() * py::self)
  .def(py::self * std::map<int, double>())
  .def(std::map<int, double>() * py::self)
  .def(py::self / double())
  .def(py::self + py::self)
  .def(py::self - py::self);

  py::class_<apfel::Set<apfel::Operator>>(m, "SetO", "Collection of objects (operators) together with a convolution map.")
  .def(py::init<apfel::ConvolutionMap const&, std::map<int, apfel::Operator> const&>(), "Map"_a = apfel::ConvolutionMap{"UNDEFINED"}, "in"_a = std::map<int, apfel::Operator> {})
  .def(py::init<std::map<int, apfel::Operator> const&>(), "in"_a)
  .def("at", &apfel::Set<apfel::Operator>::at, "id"_a)
  .def("GetMap", &apfel::Set<apfel::Operator>::GetMap)
  .def("GetObjects", &apfel::Set<apfel::Operator>::GetObjects)
  .def("SetMap", &apfel::Set<apfel::Operator>::SetMap, "map"_a)
  .def("SetObjects", &apfel::Set<apfel::Operator>::SetObjects, "objects"_a)
  .def("Combine", py::overload_cast<>(&apfel::Set<apfel::Operator>::Combine, py::const_))
  .def("Combine", py::overload_cast<std::vector<double> const&>(&apfel::Set<apfel::Operator>::Combine, py::const_), "weights"_a)
  .def("Print", &apfel::Set<apfel::Operator>::Print)
  .def(py::self *= double())
  //.def(py::self *= py::self)
  .def(py::self *= std::function<double(double const&)>())
  .def(py::self *= std::vector<double>())
  .def(py::self *= std::map<int, double>())
  .def(py::self /= double())
  .def(py::self += py::self)
  //.def(py::self -= py::self)
  .def(py::self * double())
  .def(double() * py::self)
  .def(py::self * std::function<double(double const&)>())
  .def(std::function<double(double const&)>() * py::self)
  .def(py::self * std::vector<double>())
  .def(std::vector<double>() * py::self)
  .def(py::self * std::map<int, double>())
  .def(std::map<int, double>() * py::self)
  .def(py::self / double())
  .def(py::self + py::self)
  .def(py::self - py::self);

  py::class_<apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>>(m, "SetDO", "Collection of objects (distribution-operator double objects) together with a convolution map.")
  .def(py::init<apfel::ConvolutionMap const&, std::map<int, apfel::DoubleObject<apfel::Distribution, apfel::Operator>> const&>(), "Map"_a = apfel::ConvolutionMap{"UNDEFINED"}, "in"_a = std::map<int, apfel::DoubleObject<apfel::Distribution, apfel::Operator>> {})
  .def(py::init<std::map<int, apfel::DoubleObject<apfel::Distribution, apfel::Operator>> const&>(), "in"_a)
  .def("at", &apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::at, "id"_a)
  .def("GetMap", &apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::GetMap)
  .def("GetObjects", &apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::GetObjects)
  .def("SetMap", &apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::SetMap, "map"_a)
  .def("SetObjects", &apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::SetObjects, "objects"_a)
  .def("Combine", py::overload_cast<>(&apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::Combine, py::const_))
  .def("Combine", py::overload_cast<std::vector<double> const&>(&apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::Combine, py::const_), "weights"_a)
  .def("Print", &apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::Print)
  .def(py::self *= double())
  //.def(py::self *= py::self)
  .def(py::self *= std::function<double(double const&)>())
  .def(py::self *= std::vector<double>())
  .def(py::self *= std::map<int, double>())
  .def(py::self /= double())
  .def(py::self += py::self)
  //.def(py::self -= py::self)
  .def(py::self * double())
  .def(double() * py::self)
  .def(py::self * std::function<double(double const&)>())
  .def(std::function<double(double const&)>() * py::self)
  .def(py::self * std::vector<double>())
  .def(std::vector<double>() * py::self)
  .def(py::self * std::map<int, double>())
  .def(std::map<int, double>() * py::self)
  .def(py::self / double())
  .def(py::self + py::self)
  .def(py::self - py::self);

  // Wrappers of "observable.h"
  py::class_<apfel::Observable<apfel::Distribution>::ConvolutionPair>(m, "ConvolutionPairD", "A pair of sets: coefficient functions and the objects they are convolved with.")
  .def(py::init<std::function<apfel::Set<apfel::Operator>(double const&)> const&, std::function<apfel::Set<apfel::Distribution>(double const&)> const&>(), "C"_a, "O"_a)
  .def_readwrite("CoefficientFunctions", &apfel::Observable<apfel::Distribution>::ConvolutionPair::CoefficientFunctions)
  .def_readwrite("Objects", &apfel::Observable<apfel::Distribution>::ConvolutionPair::Objects);

  py::class_<apfel::Observable<apfel::Operator>::ConvolutionPair>(m, "ConvolutionPairO", "A pair of sets: coefficient functions and the objects they are convolved with.")
  .def(py::init<std::function<apfel::Set<apfel::Operator>(double const&)> const&, std::function<apfel::Set<apfel::Operator>(double const&)> const&>(), "C"_a, "O"_a)
  .def_readwrite("CoefficientFunctions", &apfel::Observable<apfel::Operator>::ConvolutionPair::CoefficientFunctions)
  .def_readwrite("Objects", &apfel::Observable<apfel::Operator>::ConvolutionPair::Objects);

  py::class_<apfel::Observable<apfel::Distribution>>(m, "ObservableD", "Encapsulates sets of coefficient functions and objects to build a physical observable.")
  .def(py::init<std::vector<apfel::Observable<apfel::Distribution>::ConvolutionPair>>(), "ConvPair"_a)
  .def(py::init<std::function<apfel::Set<apfel::Operator>(double const&)> const&, std::function<apfel::Set<apfel::Distribution>(double const&)>>(), "CoefficientFunctions"_a, "Objects"_a)
  .def("AddConvolutionPair", &apfel::Observable<apfel::Distribution>::AddConvolutionPair, "Add a convolution pair (coefficient functions and objects).", "CoefficientFunctions"_a, "Objects"_a)
  .def("Evaluate", py::overload_cast<double const&>(&apfel::Observable<apfel::Distribution>::Evaluate, py::const_), "Return the observable as a distribution at scale Q.", "Q"_a)
  .def("Evaluate", py::overload_cast<double const&, double const&>(&apfel::Observable<apfel::Distribution>::Evaluate, py::const_), "Return the observable interpolated in x at scale Q.", "x"_a, "Q"_a)
  .def("SetObjects", &apfel::Observable<apfel::Distribution>::SetObjects, "Set the objects, keeping the same coefficient functions.", "Objects"_a, "ip"_a = 0)
  .def("GetCoefficientFunctions", &apfel::Observable<apfel::Distribution>::GetCoefficientFunctions, "Return the coefficient functions.", "ip"_a = 0);

  py::class_<apfel::Observable<apfel::Operator>>(m, "ObservableO", "Encapsulates sets of coefficient functions and objects to build a physical observable.")
  .def(py::init<std::vector<apfel::Observable<apfel::Operator>::ConvolutionPair>>(), "ConvPair"_a)
  .def(py::init<std::function<apfel::Set<apfel::Operator>(double const&)> const&, std::function<apfel::Set<apfel::Operator>(double const&)>>(), "CoefficientFunctions"_a, "Objects"_a)
  .def("AddConvolutionPair", &apfel::Observable<apfel::Operator>::AddConvolutionPair, "Add a convolution pair (coefficient functions and objects).", "CoefficientFunctions"_a, "Objects"_a)
  .def("SetObjects", &apfel::Observable<apfel::Operator>::SetObjects, "Set the objects, keeping the same coefficient functions.", "Objects"_a, "ip"_a = 0)
  .def("GetCoefficientFunctions", &apfel::Observable<apfel::Operator>::GetCoefficientFunctions, "Return the coefficient functions.", "ip"_a = 0);

  // Wrapers of "qgrid.h"
  bind_qgrid<double>(m, "QGrid");

  bind_qgrid<apfel::matrix<double>>(m, "QGridMatrix");

  bind_qgrid<apfel::Distribution>(m, "QGridD");

  bind_qgrid<apfel::Operator>(m, "QGridO");

  bind_qgrid<apfel::Set<apfel::Distribution>>(m, "QGridSetD");

  bind_qgrid<apfel::Set<apfel::Operator>>(m, "QGridSetO");

  bind_qgrid<apfel::DoubleObject<apfel::Distribution>>(m, "QGridDD");

  bind_qgrid<apfel::DoubleObject<apfel::Operator>>(m, "QGridOO");

  bind_qgrid<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>(m, "QGridDO");

  bind_qgrid<apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>>(m, "QGridSetDO");

  // Wrappers of "matchedevolution.h"
  constexpr char const* medoc = "Object evolved across heavy-quark thresholds via evolution and matching rules.";
  bind_matchedevolution<double>(m, "MatchedEvolution", "Mother class for objects evolved across heavy-quark thresholds via evolution and matching rules.");
  bind_matchedevolution<apfel::matrix<double>>(m, "MatchedEvolutionMatrix", medoc);
  bind_matchedevolution<apfel::Distribution>(m, "MatchedEvolutionD", medoc);
  bind_matchedevolution<apfel::Set<apfel::Distribution>>(m, "MatchedEvolutionSetD", medoc);
  bind_matchedevolution<apfel::DoubleObject<apfel::Distribution>>(m, "MatchedEvolutionDD", medoc);
  bind_matchedevolution<apfel::Operator>(m, "MatchedEvolutionO", medoc);
  bind_matchedevolution<apfel::Set<apfel::Operator>>(m, "MatchedEvolutionSetO", medoc);
  bind_matchedevolution<apfel::DoubleObject<apfel::Operator>>(m, "MatchedEvolutionOO", medoc);

  // Wrappers of "dglap.h"
  py::class_<apfel::Dglap<apfel::Distribution>, apfel::MatchedEvolution<apfel::Set<apfel::Distribution>>>(m, "DglapD", "DGLAP evolution of a set of distributions across heavy-quark thresholds.")
  .def(py::init<std::function<apfel::Set<apfel::Operator>(int const&, double const&)> const&, std::function<apfel::Set<apfel::Operator>(bool const&, int const&)> const&, std::function<apfel::Set<apfel::Distribution>(int const&, double const&)> const&, apfel::Set<apfel::Distribution>, double const&, std::vector<double>const&, int const&>(), "SplittingFunctions"_a, "MatchingConditions"_a, "InhomogeneousTerms"_a, "ObjRef"_a, "MuRef"_a, "Thresholds"_a, "nsteps"_a = 10)
  .def("MatchObject", &apfel::Dglap<apfel::Distribution>::MatchObject, "Match the set of distributions across a threshold.", "Up"_a, "nf"_a, "sd"_a)
  .def("Derivative", &apfel::Dglap<apfel::Distribution>::Derivative, "Return the right-hand side of the DGLAP equation.", "nf"_a, "mu"_a, "f"_a)
  .def("SetInitialDistributions", py::overload_cast<std::function<double(int const&, double const&)> const&>(&apfel::Dglap<apfel::Distribution>::SetInitialDistributions), "Set the initial-scale distributions from a function of (flavour, x).", "InDistFunc"_a)
  .def("SetInitialDistributions", py::overload_cast<std::function<std::map<int, double>(double const&)> const&>(&apfel::Dglap<apfel::Distribution>::SetInitialDistributions), "Set the initial-scale distributions from a function of x returning a flavour map.", "InDistFunc"_a)
  .def("SetInitialDistributions", py::overload_cast<std::function<std::map<int, double>(double const&, double const&)> const&, double const&>(&apfel::Dglap<apfel::Distribution>::SetInitialDistributions), "Set the initial-scale distributions from a function of (x, mu) at scale mu.", "InDistFunc"_a, "mu"_a);

  py::class_<apfel::Dglap<apfel::Operator>, apfel::MatchedEvolution<apfel::Set<apfel::Operator>>>(m, "DglapO", "DGLAP evolution of a set of operators across heavy-quark thresholds.")
  .def(py::init<std::function<apfel::Set<apfel::Operator>(int const&, double const&)> const&, std::function<apfel::Set<apfel::Operator>(bool const&, int const&)> const&, std::function<apfel::Set<apfel::Operator>(int const&, double const&)> const&, apfel::Set<apfel::Operator>, double const&, std::vector<double>const&, int const&>(), "SplittingFunctions"_a, "MatchingConditions"_a, "InhomogeneousTerms"_a, "ObjRef"_a, "MuRef"_a, "Thresholds"_a, "nsteps"_a = 10)
  .def("MatchObject", &apfel::Dglap<apfel::Operator>::MatchObject, "Match the set of operators across a threshold.", "Up"_a, "nf"_a, "sd"_a)
  .def("Derivative", &apfel::Dglap<apfel::Operator>::Derivative, "Return the right-hand side of the DGLAP equation.", "nf"_a, "mu"_a, "f"_a);

  // Wrappers of "alphaqcd.h"
  py::class_<apfel::AlphaQCD, apfel::MatchedEvolution<double>>(m, "AlphaQCD", "Running of the QCD coupling alpha_s across heavy-quark thresholds.")
  .def(py::init<double const&, double const&, std::vector<double> const&, std::vector<double> const&, int const&, int const&>(), "AlphaRef"_a, "MuRef"_a, "Masses"_a, "Thresholds"_a, "pt"_a, "nsteps"_a = 10)
  .def(py::init<double const&, double const&, std::vector<double> const&, int const&, int const&>(), "AlphaRef"_a, "MuRef"_a, "Masses"_a, "pt"_a, "nsteps"_a = 10)
  .def("MatchObject", &apfel::AlphaQCD::MatchObject, "Match the coupling across a threshold.", "Up"_a, "nf"_a, "Coup"_a)
  .def("Derivative", &apfel::AlphaQCD::Derivative, "Return the beta function (right-hand side of the running equation).", "nf"_a, "void"_a, "as"_a);

  // Wrappers of "alphaqcdg.h"
  py::class_<apfel::AlphaQCDg, apfel::MatchedEvolution<double>>(m, "AlphaQCDg", "Running of the QCD coupling alpha_s using the analytic 'g-function' solution.")
  .def(py::init<double const&, double const&, std::vector<double> const&, std::vector<double> const&, int const&, double const&>(), "AlphaRef"_a, "MuRef"_a, "Masses"_a, "Thresholds"_a, "pt"_a, "kappa"_a = 1)
  .def(py::init<double const&, double const&, std::vector<double> const&, int const&, double const&>(), "AlphaRef"_a, "MuRef"_a, "Masses"_a, "pt"_a, "kappa"_a = 1)
  .def("MatchObject", &apfel::AlphaQCDg::MatchObject, "Match the coupling across a threshold.", "Up"_a, "nf"_a, "Coup"_a)
  .def("EvolveObject", &apfel::AlphaQCDg::EvolveObject, "Evolve the coupling with nf flavours between the log-scales lnmu02 and lnmu2.", "nf"_a, "lnmu02"_a, "lnmu2"_a, "as0"_a);

  // Wrappers of "alphaqcdxi.h"
  py::class_<apfel::AlphaQCDxi, apfel::MatchedEvolution<double>>(m, "AlphaQCDxi", "Running of the QCD coupling alpha_s with renormalisation-scale variation parameter xi.")
  .def(py::init<double const&, double const&, std::vector<double> const&, std::vector<double> const&, int const&, double const&, int const&>(), "AlphaRef"_a, "MuRef"_a, "Masses"_a, "Thresholds"_a, "pt"_a, "xi"_a = 1, "nsteps"_a = 10)
  .def(py::init<double const&, double const&, std::vector<double> const&, int const&, double const&, int const&>(), "AlphaRef"_a, "MuRef"_a, "Masses"_a, "pt"_a, "xi"_a = 1, "nsteps"_a = 10)
  .def("MatchObject", &apfel::AlphaQCDxi::MatchObject, "Match the coupling across a threshold.", "Up"_a, "nf"_a, "Coup"_a)
  .def("Derivative", &apfel::AlphaQCDxi::Derivative, "Return the beta function (right-hand side of the running equation).", "nf"_a, "void"_a, "as"_a);

  // Wrappers of "alphaqed.h"
  py::class_<apfel::AlphaQED, apfel::MatchedEvolution<double>>(m, "AlphaQED", "Running of the QED coupling alpha across lepton and quark thresholds.")
  .def(py::init<double const&, double const&, std::vector<double> const&, std::vector<double> const&, int const&, int const&>(), "AlphaRef"_a, "MuRef"_a, "LeptThresholds"_a, "QuarkThresholds"_a, "pt"_a, "nsteps"_a = 10)
  .def("MatchObject", &apfel::AlphaQED::MatchObject, "Match the coupling across a threshold.", "Up"_a, "nf"_a, "Coup"_a)
  .def("Derivative", &apfel::AlphaQED::Derivative, "Return the beta function (right-hand side of the running equation).", "nfl"_a, "void"_a, "a"_a);

  // Wrappers of "alphaqcdqed.h"
  py::class_<apfel::AlphaQCDQED, apfel::MatchedEvolution<apfel::matrix<double>>>(m, "AlphaQCDQED", "Combined running of the QCD and QED couplings (2x2 coupling matrix).")
  .def(py::init<double const&, double const&, double const&, std::vector<double> const&, std::vector<double> const&, int const&, int const&>(), "AlphaQCDRef"_a, "AlphaQEDRef"_a, "MuRef"_a, "LeptThresholds"_a, "QuarkThresholds"_a, "pt"_a, "nsteps"_a = 10)
  .def("MatchObject", &apfel::AlphaQCDQED::MatchObject, "Match the couplings across a threshold.", "Up"_a, "nf"_a, "Coup"_a)
  .def("Derivative", &apfel::AlphaQCDQED::Derivative, "Return the coupled beta functions (right-hand side of the running equations).", "nfl"_a, "void"_a, "a"_a);

  // Wrappers of "dglapbuilder.h"
  py::class_<apfel::DglapObjects>(m, "DglapObjects", "Container of the operators (splitting functions and matching conditions) needed to build the DGLAP evolution.")
  .def_readwrite("Threshold", &apfel::DglapObjects::Threshold)
  .def_readwrite("SplittingFunctions", &apfel::DglapObjects::SplittingFunctions)
  .def_readwrite("MatchingConditions", &apfel::DglapObjects::MatchingConditions);

  _builders.def("BuildDglap", py::overload_cast<std::map<int, apfel::DglapObjects> const&, std::function<std::map<int, double>(double const&, double const&)> const&, double const&, int const&, std::function<double(double const&)> const&, double const&, int const&>(&apfel::BuildDglap), "Build the DGLAP evolution object from the DGLAP objects and the initial-scale distributions.", "DglapObj"_a, "InDistFunc"_a, "MuRef"_a, "PerturbativeOrder"_a, "Alphas"_a, "xi"_a = 1, "nsteps"_a = 10);
  _builders.def("BuildDglap", py::overload_cast<std::map<int, apfel::DglapObjects> const&, double const&, int const&, std::function<double(double const&)> const&, double const&, int const&>(&apfel::BuildDglap), "Build the DGLAP evolution object from the DGLAP objects (no initial distributions set).", "DglapObj"_a, "MuRef"_a, "PerturbativeOrder"_a, "Alphas"_a, "xi"_a = 1, "nsteps"_a = 10);
  _builders.def("BuildDglap", py::overload_cast<std::function<apfel::DglapObjects(double const&)> const&, std::vector<double> const&, std::function<std::map<int, double>(double const&, double const&)> const&, double const&, int const&, std::function<double(double const&)> const&, int const&>(&apfel::BuildDglap), "Build the DGLAP evolution object from a DGLAP-objects function tabulated on the thresholds.", "DglapObj"_a, "Thresholds"_a, "InDistFunc"_a, "MuRef"_a, "PerturbativeOrder"_a, "Alphas"_a, "nsteps"_a = 10);

  _initializers.def("InitializeDglapObjectsQCD", py::overload_cast<apfel::Grid const&, std::vector<double> const&, std::vector<double> const&, bool const&, double const&, bool const&, std::vector<int> const&>(&apfel::InitializeDglapObjectsQCD), "g"_a, "Masses"_a, "Thresholds"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5, "n3lo"_a = false, "IMod"_a = std::vector<int> {0, 0, 0, 0, 0, 0, 0});
  _initializers.def("InitializeDglapObjectsQCD", py::overload_cast<apfel::Grid const&, std::vector<double> const&, bool const&, double const&, bool const&, std::vector<int> const&>(&apfel::InitializeDglapObjectsQCD), "g"_a, "Thresholds"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5, "n3lo"_a = false, "IMod"_a = std::vector<int> {0, 0, 0, 0, 0, 0, 0});
  _initializers.def("InitializeDglapObjectsQCDpol", py::overload_cast<apfel::Grid const&, std::vector<double> const&, std::vector<double> const&, bool const&, double const&>(&apfel::InitializeDglapObjectsQCDpol), "g"_a, "Masses"_a, "Thresholds"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);
  _initializers.def("InitializeDglapObjectsQCDpol", py::overload_cast<apfel::Grid const&, std::vector<double> const&, bool const&, double const&>(&apfel::InitializeDglapObjectsQCDpol), "g"_a, "Thresholds"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);
  _initializers.def("InitializeDglapObjectsQCDT", py::overload_cast<apfel::Grid const&, std::vector<double> const&, std::vector<double> const&, bool const&, double const&>(&apfel::InitializeDglapObjectsQCDT), "g"_a, "Masses"_a, "Thresholds"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);
  _initializers.def("InitializeDglapObjectsQCDT", py::overload_cast<apfel::Grid const&, std::vector<double> const&, bool const&, double const&>(&apfel::InitializeDglapObjectsQCDT), "g"_a, "Thresholds"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);
  _initializers.def("InitializeDglapObjectsQCDtrans", py::overload_cast<apfel::Grid const&, std::vector<double> const&, std::vector<double> const&, bool const&, double const&>(&apfel::InitializeDglapObjectsQCDtrans), "g"_a, "Masses"_a, "Thresholds"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);
  _initializers.def("InitializeDglapObjectsQCDtrans", py::overload_cast<apfel::Grid const&, std::vector<double> const&, bool const&, double const&>(&apfel::InitializeDglapObjectsQCDtrans), "g"_a, "Thresholds"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);
  _initializers.def("InitializeDglapObjectsQCDTtrans", py::overload_cast<apfel::Grid const&, std::vector<double> const&, std::vector<double> const&, bool const&, double const&>(&apfel::InitializeDglapObjectsQCDTtrans), "g"_a, "Masses"_a, "Thresholds"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);
  _initializers.def("InitializeDglapObjectsQCDTtrans", py::overload_cast<apfel::Grid const&, std::vector<double> const&, bool const&, double const&>(&apfel::InitializeDglapObjectsQCDTtrans), "g"_a, "Thresholds"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);

  // Wrappers of "tabulateobject.h"
  constexpr char const* tabdoc = "Tabulates an object on a Q-grid for fast interpolation (specialisation of QGrid<T>).";
  bind_tabulateobject<double>(m, "TabulateObject", tabdoc);
  bind_tabulateobject<apfel::matrix<double>>(m, "TabulateObjectMatrix", tabdoc);
  bind_tabulateobject<apfel::Distribution>(m, "TabulateObjectD", tabdoc)
  .def("EvaluatexQ", py::overload_cast<double const&, double const&>(&apfel::TabulateObject<apfel::Distribution>::EvaluatexQ, py::const_), "Evaluate the tabulated distribution at (x, Q).", "x"_a, "Q"_a);
  bind_tabulateobject<apfel::Set<apfel::Distribution>>(m, "TabulateObjectSetD", tabdoc)
  .def("EvaluatexQ", py::overload_cast<int const&, double const&, double const&>(&apfel::TabulateObject<apfel::Set<apfel::Distribution>>::EvaluatexQ, py::const_), "Evaluate the i-th tabulated distribution at (x, Q).", "i"_a, "x"_a, "Q"_a)
  .def("EvaluateMapxQ", py::overload_cast<double const&, double const&>(&apfel::TabulateObject<apfel::Set<apfel::Distribution>>::EvaluateMapxQ, py::const_), "Evaluate the full flavour map at (x, Q).", "x"_a, "Q"_a);
  bind_tabulateobject<apfel::DoubleObject<apfel::Distribution>>(m, "TabulateObjectDD", tabdoc)
  .def("EvaluatexzQ", py::overload_cast<double const&, double const&, double const&>(&apfel::TabulateObject<apfel::DoubleObject<apfel::Distribution>>::EvaluatexzQ, py::const_), "Evaluate the tabulated double distribution at (x, z, Q).", "x"_a, "z"_a, "Q"_a);
  bind_tabulateobject<apfel::Operator>(m, "TabulateObjectO", tabdoc);
  bind_tabulateobject<apfel::Set<apfel::Operator>>(m, "TabulateObjectSetO", tabdoc);
  bind_tabulateobject<apfel::DoubleObject<apfel::Operator>>(m, "TabulateObjectOO", tabdoc);

  // Wrappers of "structurefunctionbuilder.h"
  py::class_<apfel::StructureFunctionObjects>(m, "StructureFunctionObjects", "Container of the coefficient functions and convolution basis needed to build a structure function.")
  .def_readwrite("nf", &apfel::StructureFunctionObjects::nf)
  .def_readwrite("P", &apfel::StructureFunctionObjects::P)
  .def_readwrite("skip", &apfel::StructureFunctionObjects::skip)
  .def_readwrite("ConvBasis", &apfel::StructureFunctionObjects::ConvBasis)
  .def_readwrite("C0", &apfel::StructureFunctionObjects::C0)
  .def_readwrite("C1", &apfel::StructureFunctionObjects::C1)
  .def_readwrite("C2", &apfel::StructureFunctionObjects::C2)
  .def_readwrite("C3", &apfel::StructureFunctionObjects::C3);

  _builders.def("BuildStructureFunctions", py::overload_cast<std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> const&, std::function<std::map<int, double>(double const&, double const&)> const&, int const&, std::function<double(double const&)> const&, std::function<std::vector<double>(double const&)> const&, double const&, double const&>(&apfel::BuildStructureFunctions), "Build the structure functions as a map of observables from the structure-function objects.", "FObj"_a, "InDistFunc"_a, "PerturbativeOrder"_a, "Alphas"_a, "Couplings"_a, "xiR"_a = 1, "xiF"_a = 1);
  _builders.def("BuildStructureFunctions", py::overload_cast<std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> const&, std::function<double(int const&, double const&, double const&)> const&, int const&, std::function<double(double const&)> const&, std::function<std::vector<double>(double const&)> const&, double const&, double const&>(&apfel::BuildStructureFunctions), "Build the structure functions as a map of observables (distributions given as a (flavour, x, Q) function).", "FObj"_a, "InDistFunc"_a, "PerturbativeOrder"_a, "Alphas"_a, "Couplings"_a, "xiR"_a = 1, "xiF"_a = 1);
  _builders.def("BuildStructureFunctions", py::overload_cast<apfel::StructureFunctionObjects const&, std::map<int, apfel::Distribution> const&, int const&, double const&, int const&, double const&, double const&>(&apfel::BuildStructureFunctions), "Build a single structure function (index k) at fixed scale from precomputed distributions.", "FObjQ"_a, "InDistFuncQ"_a, "PerturbativeOrder"_a, "AlphasQ"_a, "k"_a, "xiR"_a = 1, "xiF"_a = 1);
  _builders.def("BuildStructureFunctions", py::overload_cast<apfel::StructureFunctionObjects const&, std::map<int, apfel::Distribution> const&, int const&, double const&, double const&, double const&>(&apfel::BuildStructureFunctions), "Build the structure functions at fixed scale from precomputed distributions.", "FObjQ"_a, "InDistFuncQ"_a, "PerturbativeOrder"_a, "AlphasQ"_a, "xiR"_a = 1, "xiF"_a = 1);

  _initializers.def("InitializeF2NCObjectsZM",      py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&>(&apfel::InitializeF2NCObjectsZM),      "g"_a, "Thresholds"_a, "IntEps"_a = 1e-5);
  _initializers.def("InitializeFLNCObjectsZM",      py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&>(&apfel::InitializeFLNCObjectsZM),      "g"_a, "Thresholds"_a, "IntEps"_a = 1e-5);
  _initializers.def("InitializeF3NCObjectsZM",      py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&>(&apfel::InitializeF3NCObjectsZM),      "g"_a, "Thresholds"_a, "IntEps"_a = 1e-5);
  _initializers.def("Initializeg4NCObjectsZM",      py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&>(&apfel::Initializeg4NCObjectsZM),      "g"_a, "Thresholds"_a, "IntEps"_a = 1e-5);
  _initializers.def("InitializegLNCObjectsZM",      py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&>(&apfel::InitializegLNCObjectsZM),      "g"_a, "Thresholds"_a, "IntEps"_a = 1e-5);
  _initializers.def("Initializeg1NCObjectsZM",      py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&>(&apfel::Initializeg1NCObjectsZM),      "g"_a, "Thresholds"_a, "IntEps"_a = 1e-5);
  _initializers.def("InitializeF2CCPlusObjectsZM",  py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&>(&apfel::InitializeF2CCPlusObjectsZM),  "g"_a, "Thresholds"_a, "IntEps"_a = 1e-5);
  _initializers.def("InitializeF2CCMinusObjectsZM", py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&>(&apfel::InitializeF2CCMinusObjectsZM), "g"_a, "Thresholds"_a, "IntEps"_a = 1e-5);
  _initializers.def("InitializeFLCCPlusObjectsZM",  py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&>(&apfel::InitializeFLCCPlusObjectsZM),  "g"_a, "Thresholds"_a, "IntEps"_a = 1e-5);
  _initializers.def("InitializeFLCCMinusObjectsZM", py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&>(&apfel::InitializeFLCCMinusObjectsZM), "g"_a, "Thresholds"_a, "IntEps"_a = 1e-5);
  _initializers.def("InitializeF3CCPlusObjectsZM",  py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&>(&apfel::InitializeF3CCPlusObjectsZM),  "g"_a, "Thresholds"_a, "IntEps"_a = 1e-5);
  _initializers.def("InitializeF3CCMinusObjectsZM", py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&>(&apfel::InitializeF3CCMinusObjectsZM), "g"_a, "Thresholds"_a, "IntEps"_a = 1e-5);
  _initializers.def("InitializeF2NCObjectsZMT",     py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&>(&apfel::InitializeF2NCObjectsZMT),     "g"_a, "Thresholds"_a, "IntEps"_a = 1e-5);
  _initializers.def("InitializeFLNCObjectsZMT",     py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&>(&apfel::InitializeFLNCObjectsZMT),     "g"_a, "Thresholds"_a, "IntEps"_a = 1e-5);
  _initializers.def("InitializeF3NCObjectsZMT",     py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&>(&apfel::InitializeF3NCObjectsZMT),     "g"_a, "Thresholds"_a, "IntEps"_a = 1e-5);
  _initializers.def("InitializeF2NCObjectsMassive", py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&, int const&, double const&, double const&, int const&, double const&, double const&, double const&, double const&>(&apfel::InitializeF2NCObjectsMassive), "g"_a, "Masses"_a, "IntEps"_a = 1e-5, "nxi"_a = 150, "ximin"_a = 0.01, "ximax"_a = 10000, "intdeg"_a = 3, "lambda"_a = 0.0005, "IModThr"_a = 0, "IModAsy"_a = 0, "DampPow"_a = 1);
  _initializers.def("InitializeFLNCObjectsMassive", py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&, int const&, double const&, double const&, int const&, double const&, double const&, double const&, double const&>(&apfel::InitializeFLNCObjectsMassive), "g"_a, "Masses"_a, "IntEps"_a = 1e-5, "nxi"_a = 150, "ximin"_a = 0.01, "ximax"_a = 10000, "intdeg"_a = 3, "lambda"_a = 0.0005, "IModThr"_a = 0, "IModAsy"_a = 0, "DampPow"_a = 1);
  _initializers.def("InitializeF2NCObjectsMassiveZero", py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&, int const&, double const&, double const&, int const&, double const&>(&apfel::InitializeF2NCObjectsMassiveZero), "g"_a, "Masses"_a, "IntEps"_a = 1e-5, "nxi"_a = 150, "ximin"_a = 0.01, "ximax"_a = 10000, "intdeg"_a = 3, "lambda"_a = 0.0005);
  _initializers.def("InitializeFLNCObjectsMassiveZero", py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&, int const&, double const&, double const&, int const&, double const&>(&apfel::InitializeFLNCObjectsMassiveZero), "g"_a, "Masses"_a, "IntEps"_a = 1e-5, "nxi"_a = 150, "ximin"_a = 0.01, "ximax"_a = 10000, "intdeg"_a = 3, "lambda"_a = 0.0005);

  // Wrappers of "cffbuilder.h"
  _initializers.def("InitializeImCFF1NCObjectsZM", py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&, double const&>(&apfel::InitializeImCFF1NCObjectsZM),"g"_a, "Thresholds"_a, "xi"_a, "IntEps"_a = 1e-5);
  _initializers.def("InitializeReCFF1NCObjectsZM", py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&, double const&>(&apfel::InitializeReCFF1NCObjectsZM),"g"_a, "Thresholds"_a, "xi"_a, "IntEps"_a = 1e-5);

  // Wrappers of "dglapbuilder.h"
  _initializers.def("InitializeGpdObjects", py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&, bool const&, double const&>(&apfel::InitializeGpdObjects), "g"_a, "Thresholds"_a, "xi"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);
  _initializers.def("InitializeGpdObjectsPol", py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&, bool const&, double const&>(&apfel::InitializeGpdObjectsPol), "g"_a, "Thresholds"_a, "xi"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);
  _initializers.def("InitializeGpdObjectsTrans", py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&, bool const&, double const&>(&apfel::InitializeGpdObjectsTrans), "g"_a, "Thresholds"_a, "xi"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);

  // Wrappers of "massivecoefficientfunctionsunp_sl.h"
  py::class_<apfel::Cm21gNC, apfel::Expression>(m, "Cm21gNC", "O(as) gluon massive neutral-current coefficient function for F2.")
  .def(py::init<double const&>(), "eta"_a)
  .def("Regular", &apfel::Cm21gNC::Regular, "Regular term at x.", "x"_a);
  py::class_<apfel::CmL1gNC, apfel::Expression>(m, "CmL1gNC", "O(as) gluon massive neutral-current coefficient function for FL.")
  .def(py::init<double const&>(), "eta"_a)
  .def("Regular", &apfel::CmL1gNC::Regular, "Regular term at x.", "x"_a);
  py::class_<apfel::Cm22gNC, apfel::Expression>(m, "Cm22gNC", "O(as^2) gluon massive neutral-current coefficient function for F2.")
  .def(py::init<double const&>(), "eta"_a)
  .def("Regular", &apfel::Cm22gNC::Regular, "Regular term at x.", "x"_a);
  py::class_<apfel::CmL2gNC, apfel::Expression>(m, "CmL2gNC", "O(as^2) gluon massive neutral-current coefficient function for FL.")
  .def(py::init<double const&>(), "eta"_a)
  .def("Regular", &apfel::CmL2gNC::Regular, "Regular term at x.", "x"_a);
  py::class_<apfel::Cm22psNC, apfel::Expression>(m, "Cm22psNC", "O(as^2) pure-singlet massive neutral-current coefficient function for F2.")
  .def(py::init<double const&>(), "eta"_a)
  .def("Regular", &apfel::Cm22psNC::Regular, "Regular term at x.", "x"_a);
  py::class_<apfel::CmL2psNC, apfel::Expression>(m, "CmL2psNC", "O(as^2) pure-singlet massive neutral-current coefficient function for FL.")
  .def(py::init<double const&>(), "eta"_a)
  .def("Regular", &apfel::CmL2psNC::Regular, "Regular term at x.", "x"_a);
  py::class_<apfel::Cm22nsNC, apfel::Expression>(m, "Cm22nsNC", "O(as^2) non-singlet massive neutral-current coefficient function for F2.")
  .def(py::init<double const&>(), "eta"_a)
  .def("Regular", &apfel::Cm22nsNC::Regular, "Regular term at x.", "x"_a);
  py::class_<apfel::CmL2nsNC, apfel::Expression>(m, "CmL2nsNC", "O(as^2) non-singlet massive neutral-current coefficient function for FL.")
  .def(py::init<double const&>(), "eta"_a)
  .def("Regular", &apfel::CmL2nsNC::Regular, "Regular term at x.", "x"_a);
  py::class_<apfel::Cm22bargNC, apfel::Expression>(m, "Cm22bargNC", "O(as^2) gluon massive neutral-current coefficient function for F2, term proportional to ln(Q^2/M^2).")
  .def(py::init<double const&>(), "eta"_a)
  .def("Regular", &apfel::Cm22bargNC::Regular, "Regular term at x.", "x"_a);
  py::class_<apfel::CmL2bargNC, apfel::Expression>(m, "CmL2bargNC", "O(as^2) gluon massive neutral-current coefficient function for FL, term proportional to ln(Q^2/M^2).")
  .def(py::init<double const&>(), "eta"_a)
  .def("Regular", &apfel::CmL2bargNC::Regular, "Regular term at x.", "x"_a);
  py::class_<apfel::Cm22barpsNC, apfel::Expression>(m, "Cm22barpsNC", "O(as^2) pure-singlet massive neutral-current coefficient function for F2, term proportional to ln(Q^2/M^2).")
  .def(py::init<double const&>(), "eta"_a)
  .def("Regular", &apfel::Cm22barpsNC::Regular, "Regular term at x.", "x"_a);
  py::class_<apfel::CmL2barpsNC, apfel::Expression>(m, "CmL2barpsNC", "O(as^2) pure-singlet massive neutral-current coefficient function for FL, term proportional to ln(Q^2/M^2).")
  .def(py::init<double const&>(), "eta"_a)
  .def("Regular", &apfel::CmL2barpsNC::Regular, "Regular term at x.", "x"_a);
  py::class_<apfel::Cm2a3gNC, apfel::Expression>(m, "Cm2a3gNC", "Approximated O(as^3) gluon massive neutral-current coefficient function for F2.")
  .def(py::init<int const&, double const&, int const&, bool const&>(), "nf"_a, "eta"_a, "imod"_a = 0, "muterms"_a = true)
  .def("Regular", &apfel::Cm2a3gNC::Regular, "Regular term at x.", "x"_a);
  py::class_<apfel::Cm2a3psNC, apfel::Expression>(m, "Cm2a3psNC", "Approximated O(as^3) pure-singlet massive neutral-current coefficient function for F2.")
  .def(py::init<int const&, double const&, int const&, bool const&>(), "nf"_a, "eta"_a, "imod"_a = 0, "muterms"_a = true)
  .def("Regular", &apfel::Cm2a3psNC::Regular, "Regular term at x.", "x"_a);
  py::class_<apfel::CmLa3gNC, apfel::Expression>(m, "CmLa3gNC", "Approximated O(as^3) gluon massive neutral-current coefficient function for FL.")
  .def(py::init<int const&, double const&, int const&, bool const&>(), "nf"_a, "eta"_a, "imod"_a = 0, "muterms"_a = true)
  .def("Regular", &apfel::CmLa3gNC::Regular, "Regular term at x.", "x"_a);
  py::class_<apfel::CmLa3psNC, apfel::Expression>(m, "CmLa3psNC", "Approximated O(as^3) pure-singlet massive neutral-current coefficient function for FL.")
  .def(py::init<int const&, double const&, int const&, bool const&>(), "nf"_a, "eta"_a, "imod"_a = 0, "muterms"_a = true)
  .def("Regular", &apfel::CmLa3psNC::Regular, "Regular term at x.", "x"_a);
}
