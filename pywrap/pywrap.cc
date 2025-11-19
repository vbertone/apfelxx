#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>
#include <apfel/ogataquadrature.h>
#include <apfel/hardfactors.h>
#include <apfel/tmdbuilder.h>
#include <apfel/apfelxx.h>
#include <apfel/betaqcd.h>
#include <apfel/gammak.h>
#include <apfel/gammaf.h>
#include <apfel/kcs.h>

namespace py = pybind11;
using namespace pybind11::literals;

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

  // Hard factors
  py::module_ _hardFactors = m.def_submodule("hardFactors", "Hard factors for TMD calculations");

  // TMD evolution and matching functions
  py::module_ _tmd = m.def_submodule("tmd", "TMD evolution and matching functions");  
  
  // Ogara quadrature
  py::module_ _ogata = m.def_submodule("ogata", "Ogata quadrature for Hankel transforms");

  // Wrappers of "messages.h"
  m.def("SetVerbosityLevel", &apfel::SetVerbosityLevel, "vl"_a);
  m.def("GetVerbosityLevel", &apfel::GetVerbosityLevel);
  m.def("Banner",            &apfel::Banner);

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
  _utilities.def("LHToyFFs", &apfel::LHToyPDFsPol, "x"_a, "Q"_a);

  // Wrapper of "hardfactors.h"
  _hardFactors.def("H1DY", &apfel::H1DY, "One-loop hard factor for Drell-Yan");
  _hardFactors.def("H2DY", &apfel::H2DY, "nf"_a, "Two-loop hard factor for Drell-Yan");
  _hardFactors.def("H3DY", &apfel::H3DY, "nf"_a, "Three-loop hard factor for Drell-Yan");

  _hardFactors.def("H1SIDIS", &apfel::H1SIDIS, "One-loop hard factor for SIDIS");
  _hardFactors.def("H2SIDIS", &apfel::H2SIDIS, "nf"_a, "Two-loop hard factor for SIDIS");
  _hardFactors.def("H3SIDIS", &apfel::H3SIDIS, "nf"_a, "Three-loop hard factor for SIDIS");

  _hardFactors.def("H3Ch", &apfel::H3Ch, "Three-loop hard factor for Higgs production");
  _hardFactors.def("H1ggH", &apfel::H1ggH, "One-loop hard factor for gluon-gluon fusion Higgs production");
  _hardFactors.def("H2ggH", &apfel::H2ggH, "nf"_a, "Two-loop hard factor for gluon-gluon fusion Higgs production");

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
  py::class_<apfel::Timer>(m, "Timer")
  .def(py::init<>())
  .def("start", &apfel::Timer::start)
  .def("stop",  &apfel::Timer::stop, "ForceDisplay"_a = false);

  // Wrappers of "subgrid.h"
  py::class_<apfel::SubGrid>(m, "SubGrid")
  .def(py::init<int const&, double const&, int const&>(), "nx"_a, "xMin"_a, "InterDegree"_a)
  .def(py::init<std::vector<double> const&, int const&>(), "xsg"_a, "InterDegree"_a)
  .def("nx", &apfel::SubGrid::nx)
  .def("InterDegree", &apfel::SubGrid::InterDegree)
  .def("xMin", &apfel::SubGrid::xMin)
  .def("xMax", &apfel::SubGrid::xMax)
  .def("Step", &apfel::SubGrid::Step)
  .def("GetGrid", &apfel::SubGrid::GetGrid)
  .def("GetLogGrid", &apfel::SubGrid::GetLogGrid)
  .def("Print", &apfel::SubGrid::Print)
  .def(py::self == py::self)
  .def(py::self != py::self);

  // Wrappers of "grid.h"
  py::class_<apfel::Grid>(m, "Grid")
  .def(py::init<std::vector<apfel::SubGrid> const&>(), "grs"_a)
  .def("nGrids", &apfel::Grid::nGrids)
  .def("SubToJointMap", &apfel::Grid::SubToJointMap)
  .def("JointToSubMap", &apfel::Grid::JointToSubMap)
  .def("TransitionPoints", &apfel::Grid::TransitionPoints)
  .def("GetSubGrids", &apfel::Grid::GetSubGrids)
  .def("GetSubGrid", &apfel::Grid::GetSubGrid, "ig"_a)
  .def("GetJointGrid", &apfel::Grid::GetJointGrid)
  .def("Print", &apfel::Grid::Print)
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
  py::class_<apfel::Interpolator, PyInterpolator>(m, "Interpolator")
  .def(py::init<apfel::Grid const&>(), "gr"_a)
  .def(py::init<apfel::Grid const&, std::vector<std::vector<double>> const&, std::vector<double> const&>(), "gr"_a, "distsubgrid"_a, "distjointgrid"_a)
  .def("Evaluate", py::overload_cast<double const&>(&apfel::Interpolator::Evaluate, py::const_), "x"_a)
  .def("Evaluate", py::overload_cast<double const&, int const&>(&apfel::Interpolator::Evaluate, py::const_), "x"_a, "ig"_a)
  .def("Derive", &apfel::Interpolator::Derive, "x"_a)
  .def("Integrate", &apfel::Interpolator::Integrate, "a"_a, "b"_a)
  .def("InterpolantLog", &apfel::Interpolator::InterpolantLog, "beta"_a, "lnx"_a, "_asg"_a)
  .def("Interpolant", &apfel::Interpolator::Interpolant, "beta"_a, "x"_a, "sg"_a)
  .def("DerInterpolant", &apfel::Interpolator::DerInterpolant)
  .def("IntInterpolant", &apfel::Interpolator::IntInterpolant)
  .def("SumBounds", &apfel::Interpolator::SumBounds)
  .def("GetGrid", &apfel::Interpolator::GetGrid)
  .def("GetDistributionSubGrid", &apfel::Interpolator::GetDistributionSubGrid)
  .def("GetDistributionJointGrid", &apfel::Interpolator::GetDistributionJointGrid)
  .def("Print", &apfel::Interpolator::Print);

  // Wrappers of "lagrangeinterpolator.h"
  py::class_<apfel::LagrangeInterpolator, apfel::Interpolator>(m, "LagrangeInterpolator")
  .def(py::init<apfel::Grid const&>(), "gr"_a)
  .def(py::init<apfel::Grid const&, std::vector<std::vector<double>> const&, std::vector<double> const&>(), "gr"_a, "distsubgrid"_a, "distjointgrid"_a)
  .def("InterpolantLog", &apfel::LagrangeInterpolator::InterpolantLog, "beta"_a, "lnx"_a, "_asg"_a)
  .def("Interpolant", &apfel::LagrangeInterpolator::Interpolant, "beta"_a, "x"_a, "sg"_a)
  .def("DerInterpolant", &apfel::LagrangeInterpolator::DerInterpolant, "beta"_a, "x"_a, "sg"_a)
  .def("IntInterpolant", &apfel::LagrangeInterpolator::IntInterpolant, "beta"_a, "a"_a, "b"_a, "sg"_a)
  .def("SumBounds", &apfel::LagrangeInterpolator::SumBounds, "x"_a, "sg"_a);

  // Wrappers of "distribution.h"
  py::class_<apfel::Distribution, apfel::LagrangeInterpolator>(m, "Distribution")
  .def(py::init<apfel::Grid const&>(), "g"_a)
  .def(py::init<apfel::Distribution const&>(), "obj"_a)
  .def(py::init<apfel::Distribution const&, std::vector<std::vector<double>> const&, std::vector<double> const&>(), "obj"_a, "distsubgrid"_a, "distjointgrid"_a)
  .def(py::init<apfel::Grid const&, std::vector<std::vector<double>> const&, std::vector<double> const&>(), "g"_a, "distsubgrid"_a, "distjointgrid"_a)
  .def(py::init<apfel::Grid const&, std::function<double(double const&)>>(), "g"_a, "InDistFunc"_a)
  .def(py::init<apfel::Grid const&, std::function<double(double const&, double const&)>, double const&>(), "g"_a, "InDistFunc"_a, "Q"_a)
  .def(py::init<apfel::Grid const&, std::function<double(int const&, double const&)>, int const&>(), "g"_a, "InDistFunc"_a, "ipdf"_a)
  .def(py::init<apfel::Grid const&, std::function<double(int const&, double const&, double const&)>, int const&, double const&>(), "g"_a, "InDistFunc"_a, "ipdf"_a, "Q"_a)
  .def("SetJointGrid", static_cast<void (apfel::Distribution::*)(int const&, double const&)>(&apfel::Distribution::SetJointGrid), "ix"_a, "x"_a)
  .def("SetSubGrid", static_cast<void (apfel::Distribution::*)(int const&, int const&, double const&)>(&apfel::Distribution::SetSubGrid), "ig"_a, "ix"_a, "x"_a)
  .def("SetJointGrid", static_cast<void (apfel::Distribution::*)(std::vector<double> const&)>(&apfel::Distribution::SetJointGrid), "jg"_a)
  .def("SetSubGrid", static_cast<void (apfel::Distribution::*)(int const&, std::vector<double> const&)>(&apfel::Distribution::SetSubGrid), "ig"_a, "sg"_a)
  .def("Derivative", &apfel::Distribution::Derivative)
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
    double LocalPV(double const& x) const override
    {
      PYBIND11_OVERRIDE(double, Expression, LocalPV, x);
    };
  };
  py::class_<apfel::Expression, PyExpression>(m, "Expression")
  .def(py::init<double const&>(), "eta"_a = 1)
  .def("Regular", &apfel::Expression::Regular)
  .def("Singular", &apfel::Expression::Singular)
  .def("Local", &apfel::Expression::Local)
  .def("LocalPV", &apfel::Expression::LocalPV)
  .def("SetExternalVariable", &apfel::Expression::SetExternalVariable, "extvar"_a)
  .def("eta", &apfel::Expression::eta);

  py::class_<apfel::Identity, apfel::Expression>(m, "Identity")
  .def(py::init<>())
  .def("Local", &apfel::Expression::Local);

  py::class_<apfel::Null, apfel::Expression>(m, "Null")
  .def(py::init<>());

  // Wrappers of "matrix.h" (this is a template class and needs a
  // wrapper for any specialisation)
  py::class_<apfel::matrix<double>>(m, "matrixDouble")
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

  py::class_<apfel::matrix<std::vector<int>>>(m, "matrixVectorInt")
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

  py::class_<apfel::matrix<std::vector<double>>>(m, "matrixVectorDouble")
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
  py::class_<apfel::Operator>(m, "Operator")
  .def(py::init<apfel::Operator const&>(), "g"_a)
  //.def(py::init<apfel::Grid const&>(), "g"_a)
  .def(py::init<apfel::Grid const&, apfel::Expression const&, double const&>(), "g"_a, "expr"_a, "eps"_a = 1e-5)
  .def("GetGrid", &apfel::Operator::GetGrid)
  .def("GetOperator", &apfel::Operator::GetOperator)
  .def("Print", &apfel::Operator::Print)
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
  py::enum_<apfel::Integrator::IntegrationMethod>(m, "IntegrationMethod")
  .value("GAUSS_LEGENDRE", apfel::Integrator::IntegrationMethod::GAUSS_LEGENDRE)
  .value("GAUSS_KRONROD",  apfel::Integrator::IntegrationMethod::GAUSS_KRONROD);
  py::class_<apfel::Integrator>(m, "Integrator")
  .def(py::init<std::function<double(double const&)> const&, apfel::Integrator::IntegrationMethod const&>(), "func"_a, "method"_a = apfel::Integrator::IntegrationMethod::GAUSS_KRONROD)
  .def("integrate", py::overload_cast<double const&, double const&, double const&>(&apfel::Integrator::integrate, py::const_), "xmin"_a, "xmax"_a, "eps"_a)
  .def("integrate", py::overload_cast<double const&, double const&, std::vector<double> const&, double const&>(&apfel::Integrator::integrate, py::const_), "xmin"_a, "xmax"_a, "FixPts"_a, "eps"_a)
  .def("integrate", py::overload_cast<std::vector<double> const&, double const&>(&apfel::Integrator::integrate, py::const_), "FixPts"_a, "eps"_a)
  .def("integrate", py::overload_cast<double const&, double const&, int const&>(&apfel::Integrator::integrate, py::const_), "xmin"_a, "xmax"_a, "n"_a = 1)
  .def("integrate", py::overload_cast<double const&, double const&, std::vector<double> const&, int const&>(&apfel::Integrator::integrate, py::const_), "xmin"_a, "xmax"_a, "FixPts"_a, "n"_a = 1)
  .def("integrate", py::overload_cast<std::vector<double> const&, int const&>(&apfel::Integrator::integrate, py::const_), "FixPts"_a, "n"_a = 1)
  .def("integrand", &apfel::Integrator::integrand, "x"_a)
  .def("Method", &apfel::Integrator::Method);

  // Wrappers of "integrator2d.h"
  py::class_<apfel::Integrator2D>(m, "Integrator2D")
  .def(py::init<std::function<double(double const&, double const&)> const&, apfel::Integrator::IntegrationMethod const&>(), "func"_a, "method"_a = apfel::Integrator::IntegrationMethod::GAUSS_KRONROD)
  .def("integrate", &apfel::Integrator2D::integrate, "xmin"_a, "xmax"_a, "ymin"_a, "ymax"_a, "eps"_a)
  .def("integrand", &apfel::Integrator2D::integrand, "x"_a, "y"_a);

  // Wrappers of "doubleobject.h"
  py::class_<apfel::term<apfel::Distribution>>(m, "termD");
  py::class_<apfel::term<apfel::Operator>>(m, "termO");
  py::class_<apfel::term<apfel::Distribution, apfel::Operator>>(m, "termDO");
  py::class_<apfel::term<apfel::Operator, apfel::Distribution>>(m, "termOD");

  py::class_<apfel::DoubleObject<apfel::Distribution>>(m, "DoubleObjectD")
  .def(py::init<>())
  .def(py::init<std::vector<apfel::term<apfel::Distribution>>>(), "terms"_a)
  .def("AddTerm", &apfel::DoubleObject<apfel::Distribution>::AddTerm, "newterm"_a)
  .def("GetTerms", &apfel::DoubleObject<apfel::Distribution>::GetTerms)
  .def("Evaluate", &apfel::DoubleObject<apfel::Distribution>::Evaluate, "x"_a, "z"_a)
  .def("Evaluate1", &apfel::DoubleObject<apfel::Distribution>::Evaluate1, "x"_a)
  .def("Evaluate2", &apfel::DoubleObject<apfel::Distribution>::Evaluate2, "z"_a)
  .def("Derive", &apfel::DoubleObject<apfel::Distribution>::Derive, "x"_a, "z"_a)
  .def("Derive1", &apfel::DoubleObject<apfel::Distribution>::Derive1, "x"_a)
  .def("Derive2", &apfel::DoubleObject<apfel::Distribution>::Derive2, "z"_a)
  .def("Integrate", py::overload_cast<double const&, double const&, double const&, double const&>(&apfel::DoubleObject<apfel::Distribution>::Integrate, py::const_), "xl"_a, "xu"_a, "zl"_a, "zu"_a)
  .def("Integrate1", &apfel::DoubleObject<apfel::Distribution>::Integrate1, "xl"_a, "xu"_a)
  .def("Integrate2", &apfel::DoubleObject<apfel::Distribution>::Integrate2, "zl"_a, "zu"_a)
  .def("Integrate", py::overload_cast<double const&, double const&, std::function<double(double const&)>, std::function<double(double const&)>>(&apfel::DoubleObject<apfel::Distribution>::Integrate, py::const_), "xl"_a, "xu"_a, "zlx"_a, "zux"_a)
  .def("Integrate", py::overload_cast<std::function<double(double const&)>, std::function<double(double const&)>, double const&, double const&>(&apfel::DoubleObject<apfel::Distribution>::Integrate, py::const_), "xlz"_a, "xuz"_a, "zl"_a, "zu"_a)
  .def("MultiplyBy", &apfel::DoubleObject<apfel::Distribution>::MultiplyBy, "fx"_a, "fz"_a)
  .def("Print", &apfel::DoubleObject<apfel::Distribution>::Print)
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

  py::class_<apfel::DoubleObject<apfel::Operator>>(m, "DoubleObjectO")
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

  py::class_<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>(m, "DoubleObjectDO")
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

  py::class_<apfel::DoubleObject<apfel::Operator, apfel::Distribution>>(m, "DoubleObjectOD")
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
  py::class_<apfel::ConvolutionMap::rule>(m, "rule")
  .def_readwrite("operand", &apfel::ConvolutionMap::rule::operand)
  .def_readwrite("object", &apfel::ConvolutionMap::rule::object)
  .def_readwrite("coefficient", &apfel::ConvolutionMap::rule::coefficient)
  .def(py::self == py::self);

  py::class_<apfel::ConvolutionMap>(m, "ConvolutionMap")
  .def(py::init<std::string const&>(), "name"_a)
  .def("SetRules", &apfel::ConvolutionMap::SetRules, "rules"_a)
  .def("GetName", &apfel::ConvolutionMap::GetName)
  .def("GetRules", &apfel::ConvolutionMap::GetRules)
  .def("GetRuleMatrix", &apfel::ConvolutionMap::GetRuleMatrix)
  .def("GetRuleIndices", &apfel::ConvolutionMap::GetRuleIndices)
  .def("Print", &apfel::ConvolutionMap::Print);

  py::class_<apfel::DiagonalBasis, apfel::ConvolutionMap>(m, "DiagonalBasis")
  .def(py::init<int const&, int const&>(), "nf"_a, "offset"_a = 0);

  // Wrappers of "disbasis.h"
  py::class_<apfel::DISNCBasis, apfel::ConvolutionMap>(m, "DISNCBasis")
  .def(py::init<int const&, double const&>(), "k"_a, "fact"_a = 1)
  .def(py::init<std::vector<double> const&>(), "Ch"_a);

  py::class_<apfel::DISCCBasis, apfel::ConvolutionMap>(m, "DISCCBasis")
  .def(py::init<int const&, bool const&, double const&>(), "l"_a, "Is3"_a, "fact"_a = 1)
  .def(py::init<std::vector<double> const&, bool const&>(), "CKM"_a, "Is3"_a);

  // Wrappers of "evolutionbasis.h"
  py::class_<apfel::EvolutionBasisQCD, apfel::ConvolutionMap>(m, "EvolutionBasisQCD")
  .def(py::init<int const&>(), "nf"_a);

  py::class_<apfel::EvolutionOperatorBasisQCD, apfel::ConvolutionMap>(m, "EvolutionOperatorBasisQCD")
  .def(py::init<int const&>(), "nf"_a);

  py::class_<apfel::EvolveDistributionsBasisQCD, apfel::ConvolutionMap>(m, "EvolveDistributionsBasisQCD")
  .def(py::init<>());

  // Wrappers of "matchingbasisqcd.h"
  py::class_<apfel::MatchingBasisQCD, apfel::ConvolutionMap>(m, "MatchingBasisQCD")
  .def(py::init<int const&>(), "nf"_a);

  py::class_<apfel::MatchingOperatorBasisQCD, apfel::ConvolutionMap>(m, "MatchingOperatorBasisQCD")
  .def(py::init<int const&>(), "nf"_a);

  // Wrappers of "set.h"
  py::class_<apfel::Set<apfel::Distribution>>(m, "SetD")
  .def(py::init<apfel::ConvolutionMap const&, std::map<int, apfel::Distribution> const&>(), "Map"_a = apfel::ConvolutionMap{"UNDEFINED"}, "in"_a = std::map<int, apfel::Distribution> {})
  .def(py::init<std::map<int, apfel::Distribution> const&>(), "in"_a)
  .def("at", &apfel::Set<apfel::Distribution>::at, "id"_a)
  .def("GetMap", &apfel::Set<apfel::Distribution>::GetMap)
  .def("GetObjects", &apfel::Set<apfel::Distribution>::GetObjects)
  .def("SetMap", &apfel::Set<apfel::Distribution>::SetMap, "map"_a)
  .def("SetObjects", &apfel::Set<apfel::Distribution>::SetObjects, "objects"_a)
  .def("Combine", py::overload_cast<>(&apfel::Set<apfel::Distribution>::Combine, py::const_))
  .def("Combine", py::overload_cast<std::vector<double> const&>(&apfel::Set<apfel::Distribution>::Combine, py::const_), "weights"_a)
  .def("Print", &apfel::Set<apfel::Distribution>::Print)
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
  .def(py::self - py::self)
  // For convenience, add __len__:
  .def("__len__", [](const apfel::Set<apfel::Distribution>& self) {
    return self.GetObjects().size(); // or use self.size() if available
  });

  py::class_<apfel::Set<apfel::Operator>>(m, "SetO")
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

  py::class_<apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>>(m, "SetDO")
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
  py::class_<apfel::Observable<apfel::Distribution>::ConvolutionPair>(m, "ConvolutionPairD")
  .def(py::init<std::function<apfel::Set<apfel::Operator>(double const&)> const&, std::function<apfel::Set<apfel::Distribution>(double const&)> const&>(), "C"_a, "O"_a)
  .def_readwrite("CoefficientFunctions", &apfel::Observable<apfel::Distribution>::ConvolutionPair::CoefficientFunctions)
  .def_readwrite("Objects", &apfel::Observable<apfel::Distribution>::ConvolutionPair::Objects);

  py::class_<apfel::Observable<apfel::Operator>::ConvolutionPair>(m, "ConvolutionPairO")
  .def(py::init<std::function<apfel::Set<apfel::Operator>(double const&)> const&, std::function<apfel::Set<apfel::Operator>(double const&)> const&>(), "C"_a, "O"_a)
  .def_readwrite("CoefficientFunctions", &apfel::Observable<apfel::Operator>::ConvolutionPair::CoefficientFunctions)
  .def_readwrite("Objects", &apfel::Observable<apfel::Operator>::ConvolutionPair::Objects);

  py::class_<apfel::Observable<apfel::Distribution>>(m, "ObservableD")
  .def(py::init<std::vector<apfel::Observable<apfel::Distribution>::ConvolutionPair>>(), "ConvPair"_a)
  .def(py::init<std::function<apfel::Set<apfel::Operator>(double const&)> const&, std::function<apfel::Set<apfel::Distribution>(double const&)>>(), "CoefficientFunctions"_a, "Objects"_a)
  .def("AddConvolutionPair", &apfel::Observable<apfel::Distribution>::AddConvolutionPair, "CoefficientFunctions"_a, "Objects"_a)
  .def("Evaluate", py::overload_cast<double const&>(&apfel::Observable<apfel::Distribution>::Evaluate, py::const_), "Q"_a)
  .def("Evaluate", py::overload_cast<double const&, double const&>(&apfel::Observable<apfel::Distribution>::Evaluate, py::const_), "x"_a, "Q"_a)
  .def("SetObjects", &apfel::Observable<apfel::Distribution>::SetObjects, "Objects"_a, "ip"_a = 0)
  .def("GetCoefficientFunctions", &apfel::Observable<apfel::Distribution>::GetCoefficientFunctions, "ip"_a = 0);

  py::class_<apfel::Observable<apfel::Operator>>(m, "ObservableO")
  .def(py::init<std::vector<apfel::Observable<apfel::Operator>::ConvolutionPair>>(), "ConvPair"_a)
  .def(py::init<std::function<apfel::Set<apfel::Operator>(double const&)> const&, std::function<apfel::Set<apfel::Operator>(double const&)>>(), "CoefficientFunctions"_a, "Objects"_a)
  .def("AddConvolutionPair", &apfel::Observable<apfel::Operator>::AddConvolutionPair, "CoefficientFunctions"_a, "Objects"_a)
  .def("SetObjects", &apfel::Observable<apfel::Operator>::SetObjects, "Objects"_a, "ip"_a = 0)
  .def("GetCoefficientFunctions", &apfel::Observable<apfel::Operator>::GetCoefficientFunctions, "ip"_a = 0);

  // Wrapers of "qgrid.h"
  py::class_<apfel::QGrid<double>>(m, "QGrid")
  .def(py::init<int const&, double const&, double const&, int const&, std::vector<double> const&, std::function<double(double const&)> const&, std::function<double(double const&)>>(), "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "TabFunc"_a, "InvTabFunc"_a)
  .def(py::init<int const&, double const&, double const&, int const&, std::vector<double> const&, double const&>(), "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "Lambda"_a = 0.25)
  .def(py::init<std::vector<double> const&, int const&>(), "Qg"_a, "InterDegree"_a)
  .def("Evaluate", &apfel::QGrid<double>::Evaluate, "Q"_a)
  .def("Derive", &apfel::QGrid<double>::Derive, "Q"_a)
  .def("Integrate", &apfel::QGrid<double>::Integrate, "Qa"_a, "Qb"_a)
  .def("nQ", &apfel::QGrid<double>::nQ)
  .def("InterDegree", &apfel::QGrid<double>::InterDegree)
  .def("QMin", &apfel::QGrid<double>::QMin)
  .def("QMax", &apfel::QGrid<double>::QMax)
  .def("TabFunc", &apfel::QGrid<double>::TabFunc)
  .def("GetThresholds", &apfel::QGrid<double>::GetThresholds)
  .def("GetQGrid", &apfel::QGrid<double>::GetQGrid)
  .def("GetFQGrid", &apfel::QGrid<double>::GetFQGrid)
  .def("GetThesholdIndices", &apfel::QGrid<double>::GetThesholdIndices)
  .def("GetQGridValues", &apfel::QGrid<double>::GetQGridValues)
  .def("Interpolant", &apfel::QGrid<double>::Interpolant, "tQ"_a, "tau"_a, "fq"_a)
  .def("DerInterpolant", &apfel::QGrid<double>::DerInterpolant, "tQ"_a, "tau"_a, "Q"_a)
  .def("IntInterpolant", &apfel::QGrid<double>::IntInterpolant, "tQ"_a, "tau"_a, "Qa"_a, "Qb"_a)
  .def("Print", &apfel::QGrid<double>::Print)
  .def(py::self == py::self)
  .def(py::self != py::self);

  py::class_<apfel::QGrid<apfel::matrix<double>>>(m, "QGridMatrix")
  .def(py::init<int const&, double const&, double const&, int const&, std::vector<double> const&, std::function<double(double const&)> const&, std::function<double(double const&)>>(), "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "TabFunc"_a, "InvTabFunc"_a)
  .def(py::init<int const&, double const&, double const&, int const&, std::vector<double> const&, double const&>(), "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "Lambda"_a = 0.25)
  .def(py::init<std::vector<double> const&, int const&>(), "Qg"_a, "InterDegree"_a)
  .def("Evaluate", &apfel::QGrid<apfel::matrix<double>>::Evaluate, "Q"_a)
  .def("Derive", &apfel::QGrid<apfel::matrix<double>>::Derive, "Q"_a)
  .def("Integrate", &apfel::QGrid<apfel::matrix<double>>::Integrate, "Qa"_a, "Qb"_a)
  .def("nQ", &apfel::QGrid<apfel::matrix<double>>::nQ)
  .def("InterDegree", &apfel::QGrid<apfel::matrix<double>>::InterDegree)
  .def("QMin", &apfel::QGrid<apfel::matrix<double>>::QMin)
  .def("QMax", &apfel::QGrid<apfel::matrix<double>>::QMax)
  .def("TabFunc", &apfel::QGrid<apfel::matrix<double>>::TabFunc)
  .def("GetThresholds", &apfel::QGrid<apfel::matrix<double>>::GetThresholds)
  .def("GetQGrid", &apfel::QGrid<apfel::matrix<double>>::GetQGrid)
  .def("GetFQGrid", &apfel::QGrid<apfel::matrix<double>>::GetFQGrid)
  .def("GetThesholdIndices", &apfel::QGrid<apfel::matrix<double>>::GetThesholdIndices)
  .def("GetQGridValues", &apfel::QGrid<apfel::matrix<double>>::GetQGridValues)
  .def("Interpolant", &apfel::QGrid<apfel::matrix<double>>::Interpolant, "tQ"_a, "tau"_a, "fq"_a)
  .def("DerInterpolant", &apfel::QGrid<apfel::matrix<double>>::DerInterpolant, "tQ"_a, "tau"_a, "Q"_a)
  .def("IntInterpolant", &apfel::QGrid<apfel::matrix<double>>::IntInterpolant, "tQ"_a, "tau"_a, "Qa"_a, "Qb"_a)
  .def("Print", &apfel::QGrid<apfel::matrix<double>>::Print)
  .def(py::self == py::self)
  .def(py::self != py::self);

  py::class_<apfel::QGrid<apfel::Distribution>>(m, "QGridD")
  .def(py::init<int const&, double const&, double const&, int const&, std::vector<double> const&, std::function<double(double const&)> const&, std::function<double(double const&)>>(), "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "TabFunc"_a, "InvTabFunc"_a)
  .def(py::init<int const&, double const&, double const&, int const&, std::vector<double> const&, double const&>(), "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "Lambda"_a = 0.25)
  .def(py::init<std::vector<double> const&, int const&>(), "Qg"_a, "InterDegree"_a)
  .def("Evaluate", &apfel::QGrid<apfel::Distribution>::Evaluate, "Q"_a)
  .def("Derive", &apfel::QGrid<apfel::Distribution>::Derive, "Q"_a)
  .def("Integrate", &apfel::QGrid<apfel::Distribution>::Integrate, "Qa"_a, "Qb"_a)
  .def("nQ", &apfel::QGrid<apfel::Distribution>::nQ)
  .def("InterDegree", &apfel::QGrid<apfel::Distribution>::InterDegree)
  .def("QMin", &apfel::QGrid<apfel::Distribution>::QMin)
  .def("QMax", &apfel::QGrid<apfel::Distribution>::QMax)
  .def("TabFunc", &apfel::QGrid<apfel::Distribution>::TabFunc)
  .def("GetThresholds", &apfel::QGrid<apfel::Distribution>::GetThresholds)
  .def("GetQGrid", &apfel::QGrid<apfel::Distribution>::GetQGrid)
  .def("GetFQGrid", &apfel::QGrid<apfel::Distribution>::GetFQGrid)
  .def("GetThesholdIndices", &apfel::QGrid<apfel::Distribution>::GetThesholdIndices)
  .def("GetQGridValues", &apfel::QGrid<apfel::Distribution>::GetQGridValues)
  .def("Interpolant", &apfel::QGrid<apfel::Distribution>::Interpolant, "tQ"_a, "tau"_a, "fq"_a)
  .def("DerInterpolant", &apfel::QGrid<apfel::Distribution>::DerInterpolant, "tQ"_a, "tau"_a, "Q"_a)
  .def("IntInterpolant", &apfel::QGrid<apfel::Distribution>::IntInterpolant, "tQ"_a, "tau"_a, "Qa"_a, "Qb"_a)
  .def("Print", &apfel::QGrid<apfel::Distribution>::Print)
  .def(py::self == py::self)
  .def(py::self != py::self);

  py::class_<apfel::QGrid<apfel::Operator>>(m, "QGridO")
  .def(py::init<int const&, double const&, double const&, int const&, std::vector<double> const&, std::function<double(double const&)> const&, std::function<double(double const&)>>(), "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "TabFunc"_a, "InvTabFunc"_a)
  .def(py::init<int const&, double const&, double const&, int const&, std::vector<double> const&, double const&>(), "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "Lambda"_a = 0.25)
  .def(py::init<std::vector<double> const&, int const&>(), "Qg"_a, "InterDegree"_a)
  .def("Evaluate", &apfel::QGrid<apfel::Operator>::Evaluate, "Q"_a)
  .def("Derive", &apfel::QGrid<apfel::Operator>::Derive, "Q"_a)
  .def("Integrate", &apfel::QGrid<apfel::Operator>::Integrate, "Qa"_a, "Qb"_a)
  .def("nQ", &apfel::QGrid<apfel::Operator>::nQ)
  .def("InterDegree", &apfel::QGrid<apfel::Operator>::InterDegree)
  .def("QMin", &apfel::QGrid<apfel::Operator>::QMin)
  .def("QMax", &apfel::QGrid<apfel::Operator>::QMax)
  .def("TabFunc", &apfel::QGrid<apfel::Operator>::TabFunc)
  .def("GetThresholds", &apfel::QGrid<apfel::Operator>::GetThresholds)
  .def("GetQGrid", &apfel::QGrid<apfel::Operator>::GetQGrid)
  .def("GetFQGrid", &apfel::QGrid<apfel::Operator>::GetFQGrid)
  .def("GetThesholdIndices", &apfel::QGrid<apfel::Operator>::GetThesholdIndices)
  .def("GetQGridValues", &apfel::QGrid<apfel::Operator>::GetQGridValues)
  .def("Interpolant", &apfel::QGrid<apfel::Operator>::Interpolant, "tQ"_a, "tau"_a, "fq"_a)
  .def("DerInterpolant", &apfel::QGrid<apfel::Operator>::DerInterpolant, "tQ"_a, "tau"_a, "Q"_a)
  .def("IntInterpolant", &apfel::QGrid<apfel::Operator>::IntInterpolant, "tQ"_a, "tau"_a, "Qa"_a, "Qb"_a)
  .def("Print", &apfel::QGrid<apfel::Operator>::Print)
  .def(py::self == py::self)
  .def(py::self != py::self);

  py::class_<apfel::QGrid<apfel::Set<apfel::Distribution>>>(m, "QGridSetD")
  .def(py::init<int const&, double const&, double const&, int const&, std::vector<double> const&, std::function<double(double const&)> const&, std::function<double(double const&)>>(), "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "TabFunc"_a, "InvTabFunc"_a)
  .def(py::init<int const&, double const&, double const&, int const&, std::vector<double> const&, double const&>(), "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "Lambda"_a = 0.25)
  .def(py::init<std::vector<double> const&, int const&>(), "Qg"_a, "InterDegree"_a)
  .def("Evaluate", &apfel::QGrid<apfel::Set<apfel::Distribution>>::Evaluate, "Q"_a)
  .def("Derive", &apfel::QGrid<apfel::Set<apfel::Distribution>>::Derive, "Q"_a)
  .def("Integrate", &apfel::QGrid<apfel::Set<apfel::Distribution>>::Integrate, "Qa"_a, "Qb"_a)
  .def("nQ", &apfel::QGrid<apfel::Set<apfel::Distribution>>::nQ)
  .def("InterDegree", &apfel::QGrid<apfel::Set<apfel::Distribution>>::InterDegree)
  .def("QMin", &apfel::QGrid<apfel::Set<apfel::Distribution>>::QMin)
  .def("QMax", &apfel::QGrid<apfel::Set<apfel::Distribution>>::QMax)
  .def("TabFunc", &apfel::QGrid<apfel::Set<apfel::Distribution>>::TabFunc)
  .def("GetThresholds", &apfel::QGrid<apfel::Set<apfel::Distribution>>::GetThresholds)
  .def("GetQGrid", &apfel::QGrid<apfel::Set<apfel::Distribution>>::GetQGrid)
  .def("GetFQGrid", &apfel::QGrid<apfel::Set<apfel::Distribution>>::GetFQGrid)
  .def("GetThesholdIndices", &apfel::QGrid<apfel::Set<apfel::Distribution>>::GetThesholdIndices)
  .def("GetQGridValues", &apfel::QGrid<apfel::Set<apfel::Distribution>>::GetQGridValues)
  .def("Interpolant", &apfel::QGrid<apfel::Set<apfel::Distribution>>::Interpolant, "tQ"_a, "tau"_a, "fq"_a)
  .def("DerInterpolant", &apfel::QGrid<apfel::Set<apfel::Distribution>>::DerInterpolant, "tQ"_a, "tau"_a, "Q"_a)
  .def("IntInterpolant", &apfel::QGrid<apfel::Set<apfel::Distribution>>::IntInterpolant, "tQ"_a, "tau"_a, "Qa"_a, "Qb"_a)
  .def("Print", &apfel::QGrid<apfel::Set<apfel::Distribution>>::Print)
  .def(py::self == py::self)
  .def(py::self != py::self);

  py::class_<apfel::QGrid<apfel::Set<apfel::Operator>>>(m, "QGridSetO")
  .def(py::init<int const&, double const&, double const&, int const&, std::vector<double> const&, std::function<double(double const&)> const&, std::function<double(double const&)>>(), "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "TabFunc"_a, "InvTabFunc"_a)
  .def(py::init<int const&, double const&, double const&, int const&, std::vector<double> const&, double const&>(), "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "Lambda"_a = 0.25)
  .def(py::init<std::vector<double> const&, int const&>(), "Qg"_a, "InterDegree"_a)
  .def("Evaluate", &apfel::QGrid<apfel::Set<apfel::Operator>>::Evaluate, "Q"_a)
  .def("Derive", &apfel::QGrid<apfel::Set<apfel::Operator>>::Derive, "Q"_a)
  .def("Integrate", &apfel::QGrid<apfel::Set<apfel::Operator>>::Integrate, "Qa"_a, "Qb"_a)
  .def("nQ", &apfel::QGrid<apfel::Set<apfel::Operator>>::nQ)
  .def("InterDegree", &apfel::QGrid<apfel::Set<apfel::Operator>>::InterDegree)
  .def("QMin", &apfel::QGrid<apfel::Set<apfel::Operator>>::QMin)
  .def("QMax", &apfel::QGrid<apfel::Set<apfel::Operator>>::QMax)
  .def("TabFunc", &apfel::QGrid<apfel::Set<apfel::Operator>>::TabFunc)
  .def("GetThresholds", &apfel::QGrid<apfel::Set<apfel::Operator>>::GetThresholds)
  .def("GetQGrid", &apfel::QGrid<apfel::Set<apfel::Operator>>::GetQGrid)
  .def("GetFQGrid", &apfel::QGrid<apfel::Set<apfel::Operator>>::GetFQGrid)
  .def("GetThesholdIndices", &apfel::QGrid<apfel::Set<apfel::Operator>>::GetThesholdIndices)
  .def("GetQGridValues", &apfel::QGrid<apfel::Set<apfel::Operator>>::GetQGridValues)
  .def("Interpolant", &apfel::QGrid<apfel::Set<apfel::Operator>>::Interpolant, "tQ"_a, "tau"_a, "fq"_a)
  .def("DerInterpolant", &apfel::QGrid<apfel::Set<apfel::Operator>>::DerInterpolant, "tQ"_a, "tau"_a, "Q"_a)
  .def("IntInterpolant", &apfel::QGrid<apfel::Set<apfel::Operator>>::IntInterpolant, "tQ"_a, "tau"_a, "Qa"_a, "Qb"_a)
  .def("Print", &apfel::QGrid<apfel::Set<apfel::Operator>>::Print)
  .def(py::self == py::self)
  .def(py::self != py::self);

  py::class_<apfel::QGrid<apfel::DoubleObject<apfel::Distribution>>>(m, "QGridDD")
  .def(py::init<int const&, double const&, double const&, int const&, std::vector<double> const&, std::function<double(double const&)> const&, std::function<double(double const&)>>(), "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "TabFunc"_a, "InvTabFunc"_a)
  .def(py::init<int const&, double const&, double const&, int const&, std::vector<double> const&, double const&>(), "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "Lambda"_a = 0.25)
  .def(py::init<std::vector<double> const&, int const&>(), "Qg"_a, "InterDegree"_a)
  .def("Evaluate", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution>>::Evaluate, "Q"_a)
  .def("Derive", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution>>::Derive, "Q"_a)
  .def("Integrate", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution>>::Integrate, "Qa"_a, "Qb"_a)
  .def("nQ", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution>>::nQ)
  .def("InterDegree", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution>>::InterDegree)
  .def("QMin", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution>>::QMin)
  .def("QMax", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution>>::QMax)
  .def("TabFunc", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution>>::TabFunc)
  .def("GetThresholds", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution>>::GetThresholds)
  .def("GetQGrid", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution>>::GetQGrid)
  .def("GetFQGrid", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution>>::GetFQGrid)
  .def("GetThesholdIndices", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution>>::GetThesholdIndices)
  .def("GetQGridValues", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution>>::GetQGridValues)
  .def("Interpolant", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution>>::Interpolant, "tQ"_a, "tau"_a, "fq"_a)
  .def("DerInterpolant", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution>>::DerInterpolant, "tQ"_a, "tau"_a, "Q"_a)
  .def("IntInterpolant", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution>>::IntInterpolant, "tQ"_a, "tau"_a, "Qa"_a, "Qb"_a)
  .def("Print", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution>>::Print)
  .def(py::self == py::self)
  .def(py::self != py::self);

  py::class_<apfel::QGrid<apfel::DoubleObject<apfel::Operator>>>(m, "QGridOO")
  .def(py::init<int const&, double const&, double const&, int const&, std::vector<double> const&, std::function<double(double const&)> const&, std::function<double(double const&)>>(), "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "TabFunc"_a, "InvTabFunc"_a)
  .def(py::init<int const&, double const&, double const&, int const&, std::vector<double> const&, double const&>(), "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "Lambda"_a = 0.25)
  .def(py::init<std::vector<double> const&, int const&>(), "Qg"_a, "InterDegree"_a)
  .def("Evaluate", &apfel::QGrid<apfel::DoubleObject<apfel::Operator>>::Evaluate, "Q"_a)
  .def("Derive", &apfel::QGrid<apfel::DoubleObject<apfel::Operator>>::Derive, "Q"_a)
  .def("Integrate", &apfel::QGrid<apfel::DoubleObject<apfel::Operator>>::Integrate, "Qa"_a, "Qb"_a)
  .def("nQ", &apfel::QGrid<apfel::DoubleObject<apfel::Operator>>::nQ)
  .def("InterDegree", &apfel::QGrid<apfel::DoubleObject<apfel::Operator>>::InterDegree)
  .def("QMin", &apfel::QGrid<apfel::DoubleObject<apfel::Operator>>::QMin)
  .def("QMax", &apfel::QGrid<apfel::DoubleObject<apfel::Operator>>::QMax)
  .def("TabFunc", &apfel::QGrid<apfel::DoubleObject<apfel::Operator>>::TabFunc)
  .def("GetThresholds", &apfel::QGrid<apfel::DoubleObject<apfel::Operator>>::GetThresholds)
  .def("GetQGrid", &apfel::QGrid<apfel::DoubleObject<apfel::Operator>>::GetQGrid)
  .def("GetFQGrid", &apfel::QGrid<apfel::DoubleObject<apfel::Operator>>::GetFQGrid)
  .def("GetThesholdIndices", &apfel::QGrid<apfel::DoubleObject<apfel::Operator>>::GetThesholdIndices)
  .def("GetQGridValues", &apfel::QGrid<apfel::DoubleObject<apfel::Operator>>::GetQGridValues)
  .def("Interpolant", &apfel::QGrid<apfel::DoubleObject<apfel::Operator>>::Interpolant, "tQ"_a, "tau"_a, "fq"_a)
  .def("DerInterpolant", &apfel::QGrid<apfel::DoubleObject<apfel::Operator>>::DerInterpolant, "tQ"_a, "tau"_a, "Q"_a)
  .def("IntInterpolant", &apfel::QGrid<apfel::DoubleObject<apfel::Operator>>::IntInterpolant, "tQ"_a, "tau"_a, "Qa"_a, "Qb"_a)
  .def("Print", &apfel::QGrid<apfel::DoubleObject<apfel::Operator>>::Print)
  .def(py::self == py::self)
  .def(py::self != py::self);

  py::class_<apfel::QGrid<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>>(m, "QGridDO")
  .def(py::init<int const&, double const&, double const&, int const&, std::vector<double> const&, std::function<double(double const&)> const&, std::function<double(double const&)>>(), "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "TabFunc"_a, "InvTabFunc"_a)
  .def(py::init<int const&, double const&, double const&, int const&, std::vector<double> const&, double const&>(), "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "Lambda"_a = 0.25)
  .def(py::init<std::vector<double> const&, int const&>(), "Qg"_a, "InterDegree"_a)
  .def("Evaluate", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::Evaluate, "Q"_a)
  .def("Derive", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::Derive, "Q"_a)
  .def("Integrate", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::Integrate, "Qa"_a, "Qb"_a)
  .def("nQ", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::nQ)
  .def("InterDegree", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::InterDegree)
  .def("QMin", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::QMin)
  .def("QMax", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::QMax)
  .def("TabFunc", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::TabFunc)
  .def("GetThresholds", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::GetThresholds)
  .def("GetQGrid", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::GetQGrid)
  .def("GetFQGrid", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::GetFQGrid)
  .def("GetThesholdIndices", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::GetThesholdIndices)
  .def("GetQGridValues", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::GetQGridValues)
  .def("Interpolant", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::Interpolant, "tQ"_a, "tau"_a, "fq"_a)
  .def("DerInterpolant", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::DerInterpolant, "tQ"_a, "tau"_a, "Q"_a)
  .def("IntInterpolant", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::IntInterpolant, "tQ"_a, "tau"_a, "Qa"_a, "Qb"_a)
  .def("Print", &apfel::QGrid<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>::Print)
  .def(py::self == py::self)
  .def(py::self != py::self);

  py::class_<apfel::QGrid<apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>>>(m, "QGridSetDO")
  .def(py::init<int const&, double const&, double const&, int const&, std::vector<double> const&, std::function<double(double const&)> const&, std::function<double(double const&)>>(), "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "TabFunc"_a, "InvTabFunc"_a)
  .def(py::init<int const&, double const&, double const&, int const&, std::vector<double> const&, double const&>(), "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "Lambda"_a = 0.25)
  .def(py::init<std::vector<double> const&, int const&>(), "Qg"_a, "InterDegree"_a)
  .def("Evaluate", &apfel::QGrid<apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>>::Evaluate, "Q"_a)
  .def("Derive", &apfel::QGrid<apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>>::Derive, "Q"_a)
  .def("Integrate", &apfel::QGrid<apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>>::Integrate, "Qa"_a, "Qb"_a)
  .def("nQ", &apfel::QGrid<apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>>::nQ)
  .def("InterDegree", &apfel::QGrid<apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>>::InterDegree)
  .def("QMin", &apfel::QGrid<apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>>::QMin)
  .def("QMax", &apfel::QGrid<apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>>::QMax)
  .def("TabFunc", &apfel::QGrid<apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>>::TabFunc)
  .def("GetThresholds", &apfel::QGrid<apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>>::GetThresholds)
  .def("GetQGrid", &apfel::QGrid<apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>>::GetQGrid)
  .def("GetFQGrid", &apfel::QGrid<apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>>::GetFQGrid)
  .def("GetThesholdIndices", &apfel::QGrid<apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>>::GetThesholdIndices)
  .def("GetQGridValues", &apfel::QGrid<apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>>::GetQGridValues)
  .def("Interpolant", &apfel::QGrid<apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>>::Interpolant, "tQ"_a, "tau"_a, "fq"_a)
  .def("DerInterpolant", &apfel::QGrid<apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>>::DerInterpolant, "tQ"_a, "tau"_a, "Q"_a)
  .def("IntInterpolant", &apfel::QGrid<apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>>::IntInterpolant, "tQ"_a, "tau"_a, "Qa"_a, "Qb"_a)
  .def("Print", &apfel::QGrid<apfel::Set<apfel::DoubleObject<apfel::Distribution, apfel::Operator>>>::Print)
  .def(py::self == py::self)
  .def(py::self != py::self);

  // Wrappers of "matchedevolution.h"
  // Trampoline class for virtual class
  class PyMatchedEvolution: public apfel::MatchedEvolution<double>
  {
  public:
    using MatchedEvolution::MatchedEvolution;
    double EvolveObject(int const& nf, double const& mu02, double const& mu2, double const& Obj0) const override
    {
      PYBIND11_OVERRIDE(double, MatchedEvolution<double>, EvolveObject, nf, mu02, mu2, Obj0);
    };
    double MatchObject(bool const& Up, int const& nf, double const& Obj) const override
    {
      PYBIND11_OVERRIDE_PURE(double, MatchedEvolution<double>, MatchObject, Up, nf, Obj);
    };
    double Derivative(int const& nf, double const& Mu, double const& Obj) const override
    {
      PYBIND11_OVERRIDE_PURE(double, MatchedEvolution<double>, Derivative, nf, Mu, Obj);
    };
  };
  py::class_<apfel::MatchedEvolution<double>, PyMatchedEvolution>(m, "MatchedEvolution")
  .def(py::init<double const&, double const&, std::vector<double> const&, int const&>(), "ObjRef"_a, "MuRef"_a, "Thresholds"_a, "nsteps"_a = 10)
  .def("EvolveObject", &apfel::MatchedEvolution<double>::EvolveObject, "nf"_a, "mu02"_a, "mu2"_a, "Obj0"_a)
  .def("MatchObject", &apfel::MatchedEvolution<double>::MatchObject, "Up"_a, "nf"_a, "Obj"_a)
  .def("Derivative", &apfel::MatchedEvolution<double>::Derivative, "nf"_a, "Mu"_a, "Obj"_a)
  .def("Evaluate", &apfel::MatchedEvolution<double>::Evaluate, "mu"_a)
  .def("GetObjectRef", &apfel::MatchedEvolution<double>::GetObjectRef)
  .def("GetMuRef", &apfel::MatchedEvolution<double>::GetMuRef)
  .def("GetThresholds", &apfel::MatchedEvolution<double>::GetThresholds)
  .def("GetNumberOfSteps", &apfel::MatchedEvolution<double>::GetNumberOfSteps)
  .def("SetObjectRef", &apfel::MatchedEvolution<double>::SetObjectRef, "ObjRef"_a)
  .def("SetMuRef", &apfel::MatchedEvolution<double>::SetMuRef, "MuRef"_a)
  .def("SetNumberOfSteps", &apfel::MatchedEvolution<double>::SetNumberOfSteps, "nsteps"_a);

  // Trampoline class for virtual class
  class PyMatchedEvolutionMatrix: public apfel::MatchedEvolution<apfel::matrix<double>>
  {
  public:
    using MatchedEvolution::MatchedEvolution;
    apfel::matrix<double> EvolveObject(int const& nf, double const& mu02, double const& mu2, apfel::matrix<double> const& Obj0) const override
    {
      PYBIND11_OVERRIDE(apfel::matrix<double>, MatchedEvolution<apfel::matrix<double>>, EvolveObject, nf, mu02, mu2, Obj0);
    };
    apfel::matrix<double> MatchObject(bool const& Up, int const& nf, apfel::matrix<double> const& Obj) const override
    {
      PYBIND11_OVERRIDE_PURE(apfel::matrix<double>, MatchedEvolution<apfel::matrix<double>>, MatchObject, Up, nf, Obj);
    };
    apfel::matrix<double> Derivative(int const& nf, double const& Mu, apfel::matrix<double> const& Obj) const override
    {
      PYBIND11_OVERRIDE_PURE(apfel::matrix<double>, MatchedEvolution<apfel::matrix<double>>, Derivative, nf, Mu, Obj);
    };
  };
  py::class_<apfel::MatchedEvolution<apfel::matrix<double>>, PyMatchedEvolutionMatrix>(m, "MatchedEvolutionMatrix")
  .def(py::init<apfel::matrix<double> const&, double const&, std::vector<double> const&, int const&>(), "ObjRef"_a, "MuRef"_a, "Thresholds"_a, "nsteps"_a = 10)
  .def("EvolveObject", &apfel::MatchedEvolution<apfel::matrix<double>>::EvolveObject, "nf"_a, "mu02"_a, "mu2"_a, "Obj0"_a)
  .def("MatchObject", &apfel::MatchedEvolution<apfel::matrix<double>>::MatchObject, "Up"_a, "nf"_a, "Obj"_a)
  .def("Derivative", &apfel::MatchedEvolution<apfel::matrix<double>>::Derivative, "nf"_a, "Mu"_a, "Obj"_a)
  .def("Evaluate", &apfel::MatchedEvolution<apfel::matrix<double>>::Evaluate, "mu"_a)
  .def("GetObjectRef", &apfel::MatchedEvolution<apfel::matrix<double>>::GetObjectRef)
  .def("GetMuRef", &apfel::MatchedEvolution<apfel::matrix<double>>::GetMuRef)
  .def("GetThresholds", &apfel::MatchedEvolution<apfel::matrix<double>>::GetThresholds)
  .def("GetNumberOfSteps", &apfel::MatchedEvolution<apfel::matrix<double>>::GetNumberOfSteps)
  .def("SetObjectRef", &apfel::MatchedEvolution<apfel::matrix<double>>::SetObjectRef, "ObjRef"_a)
  .def("SetMuRef", &apfel::MatchedEvolution<apfel::matrix<double>>::SetMuRef, "MuRef"_a)
  .def("SetNumberOfSteps", &apfel::MatchedEvolution<apfel::matrix<double>>::SetNumberOfSteps, "nsteps"_a);

  // Trampoline class for virtual class
  class PyMatchedEvolutionD: public apfel::MatchedEvolution<apfel::Distribution>
  {
  public:
    using MatchedEvolution::MatchedEvolution;
    apfel::Distribution EvolveObject(int const& nf, double const& mu02, double const& mu2, apfel::Distribution const& Obj0) const override
    {
      PYBIND11_OVERRIDE(apfel::Distribution, MatchedEvolution<apfel::Distribution>, EvolveObject, nf, mu02, mu2, Obj0);
    };
    apfel::Distribution MatchObject(bool const& Up, int const& nf, apfel::Distribution const& Obj) const override
    {
      PYBIND11_OVERRIDE_PURE(apfel::Distribution, MatchedEvolution<apfel::Distribution>, MatchObject, Up, nf, Obj);
    };
    apfel::Distribution Derivative(int const& nf, double const& Mu, apfel::Distribution const& Obj) const override
    {
      PYBIND11_OVERRIDE_PURE(apfel::Distribution, MatchedEvolution<apfel::Distribution>, Derivative, nf, Mu, Obj);
    };
  };
  py::class_<apfel::MatchedEvolution<apfel::Distribution>, PyMatchedEvolutionD>(m, "MatchedEvolutionD")
  .def(py::init<apfel::Distribution const&, double const&, std::vector<double> const&, int const&>(), "ObjRef"_a, "MuRef"_a, "Thresholds"_a, "nsteps"_a = 10)
  .def("EvolveObject", &apfel::MatchedEvolution<apfel::Distribution>::EvolveObject, "nf"_a, "mu02"_a, "mu2"_a, "Obj0"_a)
  .def("MatchObject", &apfel::MatchedEvolution<apfel::Distribution>::MatchObject, "Up"_a, "nf"_a, "Obj"_a)
  .def("Derivative", &apfel::MatchedEvolution<apfel::Distribution>::Derivative, "nf"_a, "Mu"_a, "Obj"_a)
  .def("Evaluate", &apfel::MatchedEvolution<apfel::Distribution>::Evaluate, "mu"_a)
  .def("GetObjectRef", &apfel::MatchedEvolution<apfel::Distribution>::GetObjectRef)
  .def("GetMuRef", &apfel::MatchedEvolution<apfel::Distribution>::GetMuRef)
  .def("GetThresholds", &apfel::MatchedEvolution<apfel::Distribution>::GetThresholds)
  .def("GetNumberOfSteps", &apfel::MatchedEvolution<apfel::Distribution>::GetNumberOfSteps)
  .def("SetObjectRef", &apfel::MatchedEvolution<apfel::Distribution>::SetObjectRef, "ObjRef"_a)
  .def("SetMuRef", &apfel::MatchedEvolution<apfel::Distribution>::SetMuRef, "MuRef"_a)
  .def("SetNumberOfSteps", &apfel::MatchedEvolution<apfel::Distribution>::SetNumberOfSteps, "nsteps"_a);

  // Trampoline class for virtual class
  class PyMatchedEvolutionSetD: public apfel::MatchedEvolution<apfel::Set<apfel::Distribution>>
  {
  public:
    using MatchedEvolution::MatchedEvolution;
    apfel::Set<apfel::Distribution> EvolveObject(int const& nf, double const& mu02, double const& mu2, apfel::Set<apfel::Distribution> const& Obj0) const override
    {
      PYBIND11_OVERRIDE(apfel::Set<apfel::Distribution>, MatchedEvolution<apfel::Set<apfel::Distribution>>, EvolveObject, nf, mu02, mu2, Obj0);
    };
    apfel::Set<apfel::Distribution> MatchObject(bool const& Up, int const& nf, apfel::Set<apfel::Distribution> const& Obj) const override
    {
      PYBIND11_OVERRIDE_PURE(apfel::Set<apfel::Distribution>, MatchedEvolution<apfel::Set<apfel::Distribution>>, MatchObject, Up, nf, Obj);
    };
    apfel::Set<apfel::Distribution> Derivative(int const& nf, double const& Mu, apfel::Set<apfel::Distribution> const& Obj) const override
    {
      PYBIND11_OVERRIDE_PURE(apfel::Set<apfel::Distribution>, MatchedEvolution<apfel::Set<apfel::Distribution>>, Derivative, nf, Mu, Obj);
    };
  };
  py::class_<apfel::MatchedEvolution<apfel::Set<apfel::Distribution>>, PyMatchedEvolutionSetD>(m, "MatchedEvolutionSetD")
  .def(py::init<apfel::Set<apfel::Distribution> const&, double const&, std::vector<double> const&, int const&>(), "ObjRef"_a, "MuRef"_a, "Thresholds"_a, "nsteps"_a = 10)
  .def("EvolveObject", &apfel::MatchedEvolution<apfel::Set<apfel::Distribution>>::EvolveObject, "nf"_a, "mu02"_a, "mu2"_a, "Obj0"_a)
  .def("MatchObject", &apfel::MatchedEvolution<apfel::Set<apfel::Distribution>>::MatchObject, "Up"_a, "nf"_a, "Obj"_a)
  .def("Derivative", &apfel::MatchedEvolution<apfel::Set<apfel::Distribution>>::Derivative, "nf"_a, "Mu"_a, "Obj"_a)
  .def("Evaluate", &apfel::MatchedEvolution<apfel::Set<apfel::Distribution>>::Evaluate, "mu"_a)
  .def("GetObjectRef", &apfel::MatchedEvolution<apfel::Set<apfel::Distribution>>::GetObjectRef)
  .def("GetMuRef", &apfel::MatchedEvolution<apfel::Set<apfel::Distribution>>::GetMuRef)
  .def("GetThresholds", &apfel::MatchedEvolution<apfel::Set<apfel::Distribution>>::GetThresholds)
  .def("GetNumberOfSteps", &apfel::MatchedEvolution<apfel::Set<apfel::Distribution>>::GetNumberOfSteps)
  .def("SetObjectRef", &apfel::MatchedEvolution<apfel::Set<apfel::Distribution>>::SetObjectRef, "ObjRef"_a)
  .def("SetMuRef", &apfel::MatchedEvolution<apfel::Set<apfel::Distribution>>::SetMuRef, "MuRef"_a)
  .def("SetNumberOfSteps", &apfel::MatchedEvolution<apfel::Set<apfel::Distribution>>::SetNumberOfSteps, "nsteps"_a);

  // Trampoline class for virtual class
  class PyMatchedEvolutionDD: public apfel::MatchedEvolution<apfel::DoubleObject<apfel::Distribution>>
  {
  public:
    using MatchedEvolution::MatchedEvolution;
    apfel::DoubleObject<apfel::Distribution> EvolveObject(int const& nf, double const& mu02, double const& mu2, apfel::DoubleObject<apfel::Distribution> const& Obj0) const override
    {
      PYBIND11_OVERRIDE(apfel::DoubleObject<apfel::Distribution>, MatchedEvolution<apfel::DoubleObject<apfel::Distribution>>, EvolveObject, nf, mu02, mu2, Obj0);
    };
    apfel::DoubleObject<apfel::Distribution> MatchObject(bool const& Up, int const& nf, apfel::DoubleObject<apfel::Distribution> const& Obj) const override
    {
      PYBIND11_OVERRIDE_PURE(apfel::DoubleObject<apfel::Distribution>, MatchedEvolution<apfel::DoubleObject<apfel::Distribution>>, MatchObject, Up, nf, Obj);
    };
    apfel::DoubleObject<apfel::Distribution> Derivative(int const& nf, double const& Mu, apfel::DoubleObject<apfel::Distribution> const& Obj) const override
    {
      PYBIND11_OVERRIDE_PURE(apfel::DoubleObject<apfel::Distribution>, MatchedEvolution<apfel::DoubleObject<apfel::Distribution>>, Derivative, nf, Mu, Obj);
    };
  };
  py::class_<apfel::MatchedEvolution<apfel::DoubleObject<apfel::Distribution>>, PyMatchedEvolutionDD>(m, "MatchedEvolutionDD")
  .def(py::init<apfel::DoubleObject<apfel::Distribution> const&, double const&, std::vector<double> const&, int const&>(), "ObjRef"_a, "MuRef"_a, "Thresholds"_a, "nsteps"_a = 10)
  .def("EvolveObject", &apfel::MatchedEvolution<apfel::DoubleObject<apfel::Distribution>>::EvolveObject, "nf"_a, "mu02"_a, "mu2"_a, "Obj0"_a)
  .def("MatchObject", &apfel::MatchedEvolution<apfel::DoubleObject<apfel::Distribution>>::MatchObject, "Up"_a, "nf"_a, "Obj"_a)
  .def("Derivative", &apfel::MatchedEvolution<apfel::DoubleObject<apfel::Distribution>>::Derivative, "nf"_a, "Mu"_a, "Obj"_a)
  .def("Evaluate", &apfel::MatchedEvolution<apfel::DoubleObject<apfel::Distribution>>::Evaluate, "mu"_a)
  .def("GetObjectRef", &apfel::MatchedEvolution<apfel::DoubleObject<apfel::Distribution>>::GetObjectRef)
  .def("GetMuRef", &apfel::MatchedEvolution<apfel::DoubleObject<apfel::Distribution>>::GetMuRef)
  .def("GetThresholds", &apfel::MatchedEvolution<apfel::DoubleObject<apfel::Distribution>>::GetThresholds)
  .def("GetNumberOfSteps", &apfel::MatchedEvolution<apfel::DoubleObject<apfel::Distribution>>::GetNumberOfSteps)
  .def("SetObjectRef", &apfel::MatchedEvolution<apfel::DoubleObject<apfel::Distribution>>::SetObjectRef, "ObjRef"_a)
  .def("SetMuRef", &apfel::MatchedEvolution<apfel::DoubleObject<apfel::Distribution>>::SetMuRef, "MuRef"_a)
  .def("SetNumberOfSteps", &apfel::MatchedEvolution<apfel::DoubleObject<apfel::Distribution>>::SetNumberOfSteps, "nsteps"_a);

  // Trampoline class for virtual class
  class PyMatchedEvolutionO: public apfel::MatchedEvolution<apfel::Operator>
  {
  public:
    using MatchedEvolution::MatchedEvolution;
    apfel::Operator EvolveObject(int const& nf, double const& mu02, double const& mu2, apfel::Operator const& Obj0) const override
    {
      PYBIND11_OVERRIDE(apfel::Operator, MatchedEvolution<apfel::Operator>, EvolveObject, nf, mu02, mu2, Obj0);
    };
    apfel::Operator MatchObject(bool const& Up, int const& nf, apfel::Operator const& Obj) const override
    {
      PYBIND11_OVERRIDE_PURE(apfel::Operator, MatchedEvolution<apfel::Operator>, MatchObject, Up, nf, Obj);
    };
    apfel::Operator Derivative(int const& nf, double const& Mu, apfel::Operator const& Obj) const override
    {
      PYBIND11_OVERRIDE_PURE(apfel::Operator, MatchedEvolution<apfel::Operator>, Derivative, nf, Mu, Obj);
    };
  };
  py::class_<apfel::MatchedEvolution<apfel::Operator>, PyMatchedEvolutionO>(m, "MatchedEvolutionO")
  .def(py::init<apfel::Operator const&, double const&, std::vector<double> const&, int const&>(), "ObjRef"_a, "MuRef"_a, "Thresholds"_a, "nsteps"_a = 10)
  .def("EvolveObject", &apfel::MatchedEvolution<apfel::Operator>::EvolveObject, "nf"_a, "mu02"_a, "mu2"_a, "Obj0"_a)
  .def("MatchObject", &apfel::MatchedEvolution<apfel::Operator>::MatchObject, "Up"_a, "nf"_a, "Obj"_a)
  .def("Derivative", &apfel::MatchedEvolution<apfel::Operator>::Derivative, "nf"_a, "Mu"_a, "Obj"_a)
  .def("Evaluate", &apfel::MatchedEvolution<apfel::Operator>::Evaluate, "mu"_a)
  .def("GetObjectRef", &apfel::MatchedEvolution<apfel::Operator>::GetObjectRef)
  .def("GetMuRef", &apfel::MatchedEvolution<apfel::Operator>::GetMuRef)
  .def("GetThresholds", &apfel::MatchedEvolution<apfel::Operator>::GetThresholds)
  .def("GetNumberOfSteps", &apfel::MatchedEvolution<apfel::Operator>::GetNumberOfSteps)
  .def("SetObjectRef", &apfel::MatchedEvolution<apfel::Operator>::SetObjectRef, "ObjRef"_a)
  .def("SetMuRef", &apfel::MatchedEvolution<apfel::Operator>::SetMuRef, "MuRef"_a)
  .def("SetNumberOfSteps", &apfel::MatchedEvolution<apfel::Operator>::SetNumberOfSteps, "nsteps"_a);

  // Trampoline class for virtual class
  class PyMatchedEvolutionSetO: public apfel::MatchedEvolution<apfel::Set<apfel::Operator>>
  {
  public:
    using MatchedEvolution::MatchedEvolution;
    apfel::Set<apfel::Operator> EvolveObject(int const& nf, double const& mu02, double const& mu2, apfel::Set<apfel::Operator> const& Obj0) const override
    {
      PYBIND11_OVERRIDE(apfel::Set<apfel::Operator>, MatchedEvolution<apfel::Set<apfel::Operator>>, EvolveObject, nf, mu02, mu2, Obj0);
    };
    apfel::Set<apfel::Operator> MatchObject(bool const& Up, int const& nf, apfel::Set<apfel::Operator> const& Obj) const override
    {
      PYBIND11_OVERRIDE_PURE(apfel::Set<apfel::Operator>, MatchedEvolution<apfel::Set<apfel::Operator>>, MatchObject, Up, nf, Obj);
    };
    apfel::Set<apfel::Operator> Derivative(int const& nf, double const& Mu, apfel::Set<apfel::Operator> const& Obj) const override
    {
      PYBIND11_OVERRIDE_PURE(apfel::Set<apfel::Operator>, MatchedEvolution<apfel::Set<apfel::Operator>>, Derivative, nf, Mu, Obj);
    };
  };
  py::class_<apfel::MatchedEvolution<apfel::Set<apfel::Operator>>, PyMatchedEvolutionSetO>(m, "MatchedEvolutionSetO")
  .def(py::init<apfel::Set<apfel::Operator> const&, double const&, std::vector<double> const&, int const&>(), "ObjRef"_a, "MuRef"_a, "Thresholds"_a, "nsteps"_a = 10)
  .def("EvolveObject", &apfel::MatchedEvolution<apfel::Set<apfel::Operator>>::EvolveObject, "nf"_a, "mu02"_a, "mu2"_a, "Obj0"_a)
  .def("MatchObject", &apfel::MatchedEvolution<apfel::Set<apfel::Operator>>::MatchObject, "Up"_a, "nf"_a, "Obj"_a)
  .def("Derivative", &apfel::MatchedEvolution<apfel::Set<apfel::Operator>>::Derivative, "nf"_a, "Mu"_a, "Obj"_a)
  .def("Evaluate", &apfel::MatchedEvolution<apfel::Set<apfel::Operator>>::Evaluate, "mu"_a)
  .def("GetObjectRef", &apfel::MatchedEvolution<apfel::Set<apfel::Operator>>::GetObjectRef)
  .def("GetMuRef", &apfel::MatchedEvolution<apfel::Set<apfel::Operator>>::GetMuRef)
  .def("GetThresholds", &apfel::MatchedEvolution<apfel::Set<apfel::Operator>>::GetThresholds)
  .def("GetNumberOfSteps", &apfel::MatchedEvolution<apfel::Set<apfel::Operator>>::GetNumberOfSteps)
  .def("SetObjectRef", &apfel::MatchedEvolution<apfel::Set<apfel::Operator>>::SetObjectRef, "ObjRef"_a)
  .def("SetMuRef", &apfel::MatchedEvolution<apfel::Set<apfel::Operator>>::SetMuRef, "MuRef"_a)
  .def("SetNumberOfSteps", &apfel::MatchedEvolution<apfel::Set<apfel::Operator>>::SetNumberOfSteps, "nsteps"_a);

  // Trampoline class for virtual class
  class PyMatchedEvolutionOO: public apfel::MatchedEvolution<apfel::DoubleObject<apfel::Operator>>
  {
  public:
    using MatchedEvolution::MatchedEvolution;
    apfel::DoubleObject<apfel::Operator> EvolveObject(int const& nf, double const& mu02, double const& mu2, apfel::DoubleObject<apfel::Operator> const& Obj0) const override
    {
      PYBIND11_OVERRIDE(apfel::DoubleObject<apfel::Operator>, MatchedEvolution<apfel::DoubleObject<apfel::Operator>>, EvolveObject, nf, mu02, mu2, Obj0);
    };
    apfel::DoubleObject<apfel::Operator> MatchObject(bool const& Up, int const& nf, apfel::DoubleObject<apfel::Operator> const& Obj) const override
    {
      PYBIND11_OVERRIDE_PURE(apfel::DoubleObject<apfel::Operator>, MatchedEvolution<apfel::DoubleObject<apfel::Operator>>, MatchObject, Up, nf, Obj);
    };
    apfel::DoubleObject<apfel::Operator> Derivative(int const& nf, double const& Mu, apfel::DoubleObject<apfel::Operator> const& Obj) const override
    {
      PYBIND11_OVERRIDE_PURE(apfel::DoubleObject<apfel::Operator>, MatchedEvolution<apfel::DoubleObject<apfel::Operator>>, Derivative, nf, Mu, Obj);
    };
  };
  py::class_<apfel::MatchedEvolution<apfel::DoubleObject<apfel::Operator>>, PyMatchedEvolutionOO>(m, "MatchedEvolutionOO")
  .def(py::init<apfel::DoubleObject<apfel::Operator> const&, double const&, std::vector<double> const&, int const&>(), "ObjRef"_a, "MuRef"_a, "Thresholds"_a, "nsteps"_a = 10)
  .def("EvolveObject", &apfel::MatchedEvolution<apfel::DoubleObject<apfel::Operator>>::EvolveObject, "nf"_a, "mu02"_a, "mu2"_a, "Obj0"_a)
  .def("MatchObject", &apfel::MatchedEvolution<apfel::DoubleObject<apfel::Operator>>::MatchObject, "Up"_a, "nf"_a, "Obj"_a)
  .def("Derivative", &apfel::MatchedEvolution<apfel::DoubleObject<apfel::Operator>>::Derivative, "nf"_a, "Mu"_a, "Obj"_a)
  .def("Evaluate", &apfel::MatchedEvolution<apfel::DoubleObject<apfel::Operator>>::Evaluate, "mu"_a)
  .def("GetObjectRef", &apfel::MatchedEvolution<apfel::DoubleObject<apfel::Operator>>::GetObjectRef)
  .def("GetMuRef", &apfel::MatchedEvolution<apfel::DoubleObject<apfel::Operator>>::GetMuRef)
  .def("GetThresholds", &apfel::MatchedEvolution<apfel::DoubleObject<apfel::Operator>>::GetThresholds)
  .def("GetNumberOfSteps", &apfel::MatchedEvolution<apfel::DoubleObject<apfel::Operator>>::GetNumberOfSteps)
  .def("SetObjectRef", &apfel::MatchedEvolution<apfel::DoubleObject<apfel::Operator>>::SetObjectRef, "ObjRef"_a)
  .def("SetMuRef", &apfel::MatchedEvolution<apfel::DoubleObject<apfel::Operator>>::SetMuRef, "MuRef"_a)
  .def("SetNumberOfSteps", &apfel::MatchedEvolution<apfel::DoubleObject<apfel::Operator>>::SetNumberOfSteps, "nsteps"_a);

  // Wrappers of "dglap.h"
  py::class_<apfel::Dglap<apfel::Distribution>, apfel::MatchedEvolution<apfel::Set<apfel::Distribution>>>(m, "DglapD")
  .def(py::init<std::function<apfel::Set<apfel::Operator>(int const&, double const&)> const&, std::function<apfel::Set<apfel::Operator>(bool const&, int const&)> const&, apfel::Set<apfel::Distribution>, double const&, std::vector<double>const&, int const&>(), "SplittingFunctions"_a, "MatchingConditions"_a, "ObjRef"_a, "MuRef"_a, "Thresholds"_a, "nsteps"_a = 10)
  .def("MatchObject", &apfel::Dglap<apfel::Distribution>::MatchObject, "Up"_a, "nf"_a, "sd"_a)
  .def("Derivative", &apfel::Dglap<apfel::Distribution>::Derivative, "nf"_a, "mu"_a, "f"_a)
  .def("SetInitialDistributions", py::overload_cast<std::function<double(int const&, double const&)> const&>(&apfel::Dglap<apfel::Distribution>::SetInitialDistributions), "InDistFunc"_a)
  .def("SetInitialDistributions", py::overload_cast<std::function<std::map<int, double>(double const&)> const&>(&apfel::Dglap<apfel::Distribution>::SetInitialDistributions), "InDistFunc"_a)
  .def("SetInitialDistributions", py::overload_cast<std::function<std::map<int, double>(double const&, double const&)> const&, double const&>(&apfel::Dglap<apfel::Distribution>::SetInitialDistributions), "InDistFunc"_a, "mu"_a);

  py::class_<apfel::Dglap<apfel::Operator>, apfel::MatchedEvolution<apfel::Set<apfel::Operator>>>(m, "DglapO")
  .def(py::init<std::function<apfel::Set<apfel::Operator>(int const&, double const&)> const&, std::function<apfel::Set<apfel::Operator>(bool const&, int const&)> const&, apfel::Set<apfel::Operator>, double const&, std::vector<double>const&, int const&>(), "SplittingFunctions"_a, "MatchingConditions"_a, "ObjRef"_a, "MuRef"_a, "Thresholds"_a, "nsteps"_a = 10)
  .def("MatchObject", &apfel::Dglap<apfel::Operator>::MatchObject, "Up"_a, "nf"_a, "sd"_a)
  .def("Derivative", &apfel::Dglap<apfel::Operator>::Derivative, "nf"_a, "mu"_a, "f"_a);

  // Wrappers of "alphaqcd.h"
  py::class_<apfel::AlphaQCD, apfel::MatchedEvolution<double>>(m, "AlphaQCD")
  .def(py::init<double const&, double const&, std::vector<double> const&, std::vector<double> const&, int const&, int const&>(), "AlphaRef"_a, "MuRef"_a, "Masses"_a, "Thresholds"_a, "pt"_a, "nsteps"_a = 10)
  .def(py::init<double const&, double const&, std::vector<double> const&, int const&, int const&>(), "AlphaRef"_a, "MuRef"_a, "Masses"_a, "pt"_a, "nsteps"_a = 10)
  .def("MatchObject", &apfel::AlphaQCD::MatchObject, "Up"_a, "nf"_a, "Coup"_a)
  .def("Derivative", &apfel::AlphaQCD::Derivative, "nf"_a, "void"_a, "as"_a);

  // Wrappers of "alphaqcdg.h"
  py::class_<apfel::AlphaQCDg, apfel::MatchedEvolution<double>>(m, "AlphaQCDg")
  .def(py::init<double const&, double const&, std::vector<double> const&, std::vector<double> const&, int const&, double const&>(), "AlphaRef"_a, "MuRef"_a, "Masses"_a, "Thresholds"_a, "pt"_a, "kappa"_a = 1)
  .def(py::init<double const&, double const&, std::vector<double> const&, int const&, double const&>(), "AlphaRef"_a, "MuRef"_a, "Masses"_a, "pt"_a, "kappa"_a = 1)
  .def("MatchObject", &apfel::AlphaQCDg::MatchObject, "Up"_a, "nf"_a, "Coup"_a)
  .def("EvolveObject", &apfel::AlphaQCDg::EvolveObject, "nf"_a, "lnmu02"_a, "lnmu2"_a, "as0"_a);

  // Wrappers of "alphaqcdxi.h"
  py::class_<apfel::AlphaQCDxi, apfel::MatchedEvolution<double>>(m, "AlphaQCDxi")
  .def(py::init<double const&, double const&, std::vector<double> const&, std::vector<double> const&, int const&, double const&, int const&>(), "AlphaRef"_a, "MuRef"_a, "Masses"_a, "Thresholds"_a, "pt"_a, "xi"_a = 1, "nsteps"_a = 10)
  .def(py::init<double const&, double const&, std::vector<double> const&, int const&, double const&, int const&>(), "AlphaRef"_a, "MuRef"_a, "Masses"_a, "pt"_a, "xi"_a = 1, "nsteps"_a = 10)
  .def("MatchObject", &apfel::AlphaQCDxi::MatchObject, "Up"_a, "nf"_a, "Coup"_a)
  .def("Derivative", &apfel::AlphaQCDxi::Derivative, "nf"_a, "void"_a, "as"_a);

  // Wrappers of "alphaqed.h"
  py::class_<apfel::AlphaQED, apfel::MatchedEvolution<double>>(m, "AlphaQED")
  .def(py::init<double const&, double const&, std::vector<double> const&, std::vector<double> const&, int const&, int const&>(), "AlphaRef"_a, "MuRef"_a, "LeptThresholds"_a, "QuarkThresholds"_a, "pt"_a, "nsteps"_a = 10)
  .def("MatchObject", &apfel::AlphaQED::MatchObject, "Up"_a, "nf"_a, "Coup"_a)
  .def("Derivative", &apfel::AlphaQED::Derivative, "nfl"_a, "void"_a, "a"_a);

  // Wrappers of "alphaqcdqed.h"
  py::class_<apfel::AlphaQCDQED, apfel::MatchedEvolution<apfel::matrix<double>>>(m, "AlphaQCDQED")
  .def(py::init<double const&, double const&, double const&, std::vector<double> const&, std::vector<double> const&, int const&, int const&>(), "AlphaQCDRef"_a, "AlphaQEDRef"_a, "MuRef"_a, "LeptThresholds"_a, "QuarkThresholds"_a, "pt"_a, "nsteps"_a = 10)
  .def("MatchObject", &apfel::AlphaQCDQED::MatchObject, "Up"_a, "nf"_a, "Coup"_a)
  .def("Derivative", &apfel::AlphaQCDQED::Derivative, "nfl"_a, "void"_a, "a"_a);

  // Wrappers of "dglapbuilder.h"
  py::class_<apfel::DglapObjects>(m, "DglapObjects")
  .def_readwrite("Threshold", &apfel::DglapObjects::Threshold)
  .def_readwrite("SplittingFunctions", &apfel::DglapObjects::SplittingFunctions)
  .def_readwrite("MatchingConditions", &apfel::DglapObjects::MatchingConditions);

  _builders.def("BuildDglap", py::overload_cast<std::map<int, apfel::DglapObjects> const&, std::function<std::map<int, double>(double const&, double const&)> const&, double const&, int const&, std::function<double(double const&)> const&, double const&, int const&>(&apfel::BuildDglap), "DglapObj"_a, "InDistFunc"_a, "MuRef"_a, "PerturbativeOrder"_a, "Alphas"_a, "xi"_a = 1, "nsteps"_a = 10);
  _builders.def("BuildDglap", py::overload_cast<std::map<int, apfel::DglapObjects> const&, double const&, int const&, std::function<double(double const&)> const&, double const&, int const&>(&apfel::BuildDglap), "DglapObj"_a, "MuRef"_a, "PerturbativeOrder"_a, "Alphas"_a, "xi"_a = 1, "nsteps"_a = 10);
  _builders.def("BuildDglap", py::overload_cast<std::function<apfel::DglapObjects(double const&)> const&, std::vector<double> const&, std::function<std::map<int, double>(double const&, double const&)> const&, double const&, int const&, std::function<double(double const&)> const&, int const&>(&apfel::BuildDglap), "DglapObj"_a, "Thresholds"_a, "InDistFunc"_a, "MuRef"_a, "PerturbativeOrder"_a, "Alphas"_a, "nsteps"_a = 10);

  _initializers.def("InitializeDglapObjectsQCD", py::overload_cast<apfel::Grid const&, std::vector<double> const&, std::vector<double> const&, bool const&, double const&>(&apfel::InitializeDglapObjectsQCD), "g"_a, "Masses"_a, "Thresholds"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);
  _initializers.def("InitializeDglapObjectsQCD", py::overload_cast<apfel::Grid const&, std::vector<double> const&, bool const&, double const&>(&apfel::InitializeDglapObjectsQCD), "g"_a, "Thresholds"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);
  _initializers.def("InitializeDglapObjectsQCDpol", py::overload_cast<apfel::Grid const&, std::vector<double> const&, std::vector<double> const&, bool const&, double const&>(&apfel::InitializeDglapObjectsQCDpol), "g"_a, "Masses"_a, "Thresholds"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);
  _initializers.def("InitializeDglapObjectsQCDpol", py::overload_cast<apfel::Grid const&, std::vector<double> const&, bool const&, double const&>(&apfel::InitializeDglapObjectsQCDpol), "g"_a, "Thresholds"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);
  _initializers.def("InitializeDglapObjectsQCDT", py::overload_cast<apfel::Grid const&, std::vector<double> const&, std::vector<double> const&, bool const&, double const&>(&apfel::InitializeDglapObjectsQCDT), "g"_a, "Masses"_a, "Thresholds"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);
  _initializers.def("InitializeDglapObjectsQCDT", py::overload_cast<apfel::Grid const&, std::vector<double> const&, bool const&, double const&>(&apfel::InitializeDglapObjectsQCDT), "g"_a, "Thresholds"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);
  _initializers.def("InitializeDglapObjectsQCDtrans", py::overload_cast<apfel::Grid const&, std::vector<double> const&, std::vector<double> const&, bool const&, double const&>(&apfel::InitializeDglapObjectsQCDtrans), "g"_a, "Masses"_a, "Thresholds"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);
  _initializers.def("InitializeDglapObjectsQCDtrans", py::overload_cast<apfel::Grid const&, std::vector<double> const&, bool const&, double const&>(&apfel::InitializeDglapObjectsQCDtrans), "g"_a, "Thresholds"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);
  _initializers.def("InitializeDglapObjectsQCDTtrans", py::overload_cast<apfel::Grid const&, std::vector<double> const&, std::vector<double> const&, bool const&, double const&>(&apfel::InitializeDglapObjectsQCDTtrans), "g"_a, "Masses"_a, "Thresholds"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);
  _initializers.def("InitializeDglapObjectsQCDTtrans", py::overload_cast<apfel::Grid const&, std::vector<double> const&, bool const&, double const&>(&apfel::InitializeDglapObjectsQCDTtrans), "g"_a, "Thresholds"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);

  // Wrappers of "tabulateobject.h"
  py::class_<apfel::TabulateObject<double>, apfel::QGrid<double>>(m, "TabulateObject")
  .def(py::init<apfel::MatchedEvolution<double>&, int const&, double const&, double const&, int const&, double const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Lambda"_a = 0.25)
  .def(py::init<std::function<double(double const&)> const&, int const&, double const&, double const&, int const&, std::vector<double> const&, double const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "Lambda"_a = 0.25)
  .def(py::init<std::function<double(double const&)> const&, int const&, double const&, double const&, int const&, std::vector<double> const&, std::function<double(double const&)> const&, std::function<double(double const&)> const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "TabFunc"_a, "InvTabFunc"_a)
  .def(py::init<std::function<double(double const&)> const&, std::vector<double> const&, int const&>(), "Object"_a, "Qg"_a, "InterDegree"_a);

  py::class_<apfel::TabulateObject<apfel::matrix<double>>, apfel::QGrid<apfel::matrix<double>>>(m, "TabulateObjectMatrix")
  .def(py::init<apfel::MatchedEvolution<apfel::matrix<double>>&, int const&, double const&, double const&, int const&, double const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Lambda"_a = 0.25)
  .def(py::init<std::function<apfel::matrix<double>(double const&)> const&, int const&, double const&, double const&, int const&, std::vector<double> const&, double const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "Lambda"_a = 0.25)
  .def(py::init<std::function<apfel::matrix<double>(double const&)> const&, int const&, double const&, double const&, int const&, std::vector<double> const&, std::function<double(double const&)> const&, std::function<double(double const&)> const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "TabFunc"_a, "InvTabFunc"_a)
  .def(py::init<std::function<apfel::matrix<double>(double const&)> const&, std::vector<double> const&, int const&>(), "Object"_a, "Qg"_a, "InterDegree"_a);

  py::class_<apfel::TabulateObject<apfel::Distribution>, apfel::QGrid<apfel::Distribution>>(m, "TabulateObjectD")
  .def(py::init<apfel::MatchedEvolution<apfel::Distribution>&, int const&, double const&, double const&, int const&, double const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Lambda"_a = 0.25)
  .def(py::init<std::function<apfel::Distribution(double const&)> const&, int const&, double const&, double const&, int const&, std::vector<double> const&, double const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "Lambda"_a = 0.25)
  .def(py::init<std::function<apfel::Distribution(double const&)> const&, int const&, double const&, double const&, int const&, std::vector<double> const&, std::function<double(double const&)> const&, std::function<double(double const&)> const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "TabFunc"_a, "InvTabFunc"_a)
  .def(py::init<std::function<apfel::Distribution(double const&)> const&, std::vector<double> const&, int const&>(), "Object"_a, "Qg"_a, "InterDegree"_a)
  .def("EvaluatexQ", py::overload_cast<double const&, double const&>(&apfel::TabulateObject<apfel::Distribution>::EvaluatexQ, py::const_),"x"_a, "Q"_a);

  py::class_<apfel::TabulateObject<apfel::Set<apfel::Distribution>>, apfel::QGrid<apfel::Set<apfel::Distribution>>>(m, "TabulateObjectSetD")
  .def(py::init<apfel::MatchedEvolution<apfel::Set<apfel::Distribution>>&, int const&, double const&, double const&, int const&, double const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Lambda"_a = 0.25)
  .def(py::init<std::function<apfel::Set<apfel::Distribution>(double const&)> const&, int const&, double const&, double const&, int const&, std::vector<double> const&, double const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "Lambda"_a = 0.25)
  .def(py::init<std::function<apfel::Set<apfel::Distribution>(double const&)> const&, int const&, double const&, double const&, int const&, std::vector<double> const&, std::function<double(double const&)> const&, std::function<double(double const&)> const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "TabFunc"_a, "InvTabFunc"_a)
  .def(py::init<std::function<apfel::Set<apfel::Distribution>(double const&)> const&, std::vector<double> const&, int const&>(), "Object"_a, "Qg"_a, "InterDegree"_a)
  .def("EvaluatexQ", py::overload_cast<int const&, double const&, double const&>(&apfel::TabulateObject<apfel::Set<apfel::Distribution>>::EvaluatexQ, py::const_), "i"_a, "x"_a, "Q"_a)
  .def("EvaluateMapxQ", py::overload_cast<double const&, double const&>(&apfel::TabulateObject<apfel::Set<apfel::Distribution>>::EvaluateMapxQ, py::const_), "x"_a, "Q"_a);

  py::class_<apfel::TabulateObject<apfel::DoubleObject<apfel::Distribution>>, apfel::QGrid<apfel::DoubleObject<apfel::Distribution>>>(m, "TabulateObjectDD")
  .def(py::init<apfel::MatchedEvolution<apfel::DoubleObject<apfel::Distribution>>&, int const&, double const&, double const&, int const&, double const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Lambda"_a = 0.25)
  .def(py::init<std::function<apfel::DoubleObject<apfel::Distribution>(double const&)> const&, int const&, double const&, double const&, int const&, std::vector<double> const&, double const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "Lambda"_a = 0.25)
  .def(py::init<std::function<apfel::DoubleObject<apfel::Distribution>(double const&)> const&, int const&, double const&, double const&, int const&, std::vector<double> const&, std::function<double(double const&)> const&, std::function<double(double const&)> const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "TabFunc"_a, "InvTabFunc"_a)
  .def(py::init<std::function<apfel::DoubleObject<apfel::Distribution>(double const&)> const&, std::vector<double> const&, int const&>(), "Object"_a, "Qg"_a, "InterDegree"_a)
  .def("EvaluatexzQ", py::overload_cast<double const&, double const&, double const&>(&apfel::TabulateObject<apfel::DoubleObject<apfel::Distribution>>::EvaluatexzQ, py::const_), "x"_a, "z"_a, "Q"_a);

  py::class_<apfel::TabulateObject<apfel::Operator>, apfel::QGrid<apfel::Operator>>(m, "TabulateObjectO")
  .def(py::init<apfel::MatchedEvolution<apfel::Operator>&, int const&, double const&, double const&, int const&, double const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Lambda"_a = 0.25)
  .def(py::init<std::function<apfel::Operator(double const&)> const&, int const&, double const&, double const&, int const&, std::vector<double> const&, double const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "Lambda"_a = 0.25)
  .def(py::init<std::function<apfel::Operator(double const&)> const&, int const&, double const&, double const&, int const&, std::vector<double> const&, std::function<double(double const&)> const&, std::function<double(double const&)> const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "TabFunc"_a, "InvTabFunc"_a)
  .def(py::init<std::function<apfel::Operator(double const&)> const&, std::vector<double> const&, int const&>(), "Object"_a, "Qg"_a, "InterDegree"_a);

  py::class_<apfel::TabulateObject<apfel::Set<apfel::Operator>>, apfel::QGrid<apfel::Set<apfel::Operator>>>(m, "TabulateObjectSetO")
  .def(py::init<apfel::MatchedEvolution<apfel::Set<apfel::Operator>>&, int const&, double const&, double const&, int const&, double const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Lambda"_a = 0.25)
  .def(py::init<std::function<apfel::Set<apfel::Operator>(double const&)> const&, int const&, double const&, double const&, int const&, std::vector<double> const&, double const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "Lambda"_a = 0.25)
  .def(py::init<std::function<apfel::Set<apfel::Operator>(double const&)> const&, int const&, double const&, double const&, int const&, std::vector<double> const&, std::function<double(double const&)> const&, std::function<double(double const&)> const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "TabFunc"_a, "InvTabFunc"_a)
  .def(py::init<std::function<apfel::Set<apfel::Operator>(double const&)> const&, std::vector<double> const&, int const&>(), "Object"_a, "Qg"_a, "InterDegree"_a);

  py::class_<apfel::TabulateObject<apfel::DoubleObject<apfel::Operator>>, apfel::QGrid<apfel::DoubleObject<apfel::Operator>>>(m, "TabulateObjectOO")
  .def(py::init<apfel::MatchedEvolution<apfel::DoubleObject<apfel::Operator>>&, int const&, double const&, double const&, int const&, double const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Lambda"_a = 0.25)
  .def(py::init<std::function<apfel::DoubleObject<apfel::Operator>(double const&)> const&, int const&, double const&, double const&, int const&, std::vector<double> const&, double const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "Lambda"_a = 0.25)
  .def(py::init<std::function<apfel::DoubleObject<apfel::Operator>(double const&)> const&, int const&, double const&, double const&, int const&, std::vector<double> const&, std::function<double(double const&)> const&, std::function<double(double const&)> const&>(), "Object"_a, "nQ"_a, "QMin"_a, "QMax"_a, "InterDegree"_a, "Thresholds"_a, "TabFunc"_a, "InvTabFunc"_a)
  .def(py::init<std::function<apfel::DoubleObject<apfel::Operator>(double const&)> const&, std::vector<double> const&, int const&>(), "Object"_a, "Qg"_a, "InterDegree"_a);

  // Wrappers of "structurefunctionbuilder.h"
  py::class_<apfel::StructureFunctionObjects>(m, "StructureFunctionObjects")
  .def_readwrite("nf", &apfel::StructureFunctionObjects::nf)
  .def_readwrite("P", &apfel::StructureFunctionObjects::P)
  .def_readwrite("skip", &apfel::StructureFunctionObjects::skip)
  .def_readwrite("ConvBasis", &apfel::StructureFunctionObjects::ConvBasis)
  .def_readwrite("C0", &apfel::StructureFunctionObjects::C0)
  .def_readwrite("C1", &apfel::StructureFunctionObjects::C1)
  .def_readwrite("C2", &apfel::StructureFunctionObjects::C2)
  .def_readwrite("C3", &apfel::StructureFunctionObjects::C3);

  _builders.def("BuildStructureFunctions", py::overload_cast<std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> const&, std::function<std::map<int, double>(double const&, double const&)> const&, int const&, std::function<double(double const&)> const&, std::function<std::vector<double>(double const&)> const&, double const&, double const&>(&apfel::BuildStructureFunctions), "FObj"_a, "InDistFunc"_a, "PerturbativeOrder"_a, "Alphas"_a, "Couplings"_a, "xiR"_a = 1, "xiF"_a = 1);
  _builders.def("BuildStructureFunctions", py::overload_cast<std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> const&, std::function<double(int const&, double const&, double const&)> const&, int const&, std::function<double(double const&)> const&, std::function<std::vector<double>(double const&)> const&, double const&, double const&>(&apfel::BuildStructureFunctions), "FObj"_a, "InDistFunc"_a, "PerturbativeOrder"_a, "Alphas"_a, "Couplings"_a, "xiR"_a = 1, "xiF"_a = 1);
  _builders.def("BuildStructureFunctions", py::overload_cast<apfel::StructureFunctionObjects const&, std::map<int, apfel::Distribution> const&, int const&, double const&, int const&, double const&, double const&>(&apfel::BuildStructureFunctions), "FObjQ"_a, "InDistFuncQ"_a, "PerturbativeOrder"_a, "AlphasQ"_a, "k"_a, "xiR"_a = 1, "xiF"_a = 1);
  _builders.def("BuildStructureFunctions", py::overload_cast<apfel::StructureFunctionObjects const&, std::map<int, apfel::Distribution> const&, int const&, double const&, double const&, double const&>(&apfel::BuildStructureFunctions), "FObjQ"_a, "InDistFuncQ"_a, "PerturbativeOrder"_a, "AlphasQ"_a, "xiR"_a = 1, "xiF"_a = 1);

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
  _initializers.def("InitializeF2NCObjectsMassive", py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&, int const&, double const&, double const&, int const&, double const&>(&apfel::InitializeF2NCObjectsMassive), "g"_a, "Masses"_a, "IntEps"_a = 1e-5, "nxi"_a = 150, "ximin"_a = 0.01, "ximax"_a = 10000, "intdeg"_a = 3, "lambda"_a = 0.0005);
  _initializers.def("InitializeFLNCObjectsMassive", py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&, int const&, double const&, double const&, int const&, double const&>(&apfel::InitializeFLNCObjectsMassive), "g"_a, "Masses"_a, "IntEps"_a = 1e-5, "nxi"_a = 150, "ximin"_a = 0.01, "ximax"_a = 10000, "intdeg"_a = 3, "lambda"_a = 0.0005);
  _initializers.def("InitializeF2NCObjectsMassiveZero", py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&, int const&, double const&, double const&, int const&, double const&>(&apfel::InitializeF2NCObjectsMassiveZero), "g"_a, "Masses"_a, "IntEps"_a = 1e-5, "nxi"_a = 150, "ximin"_a = 0.01, "ximax"_a = 10000, "intdeg"_a = 3, "lambda"_a = 0.0005);
  _initializers.def("InitializeFLNCObjectsMassiveZero", py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&, int const&, double const&, double const&, int const&, double const&>(&apfel::InitializeFLNCObjectsMassiveZero), "g"_a, "Masses"_a, "IntEps"_a = 1e-5, "nxi"_a = 150, "ximin"_a = 0.01, "ximax"_a = 10000, "intdeg"_a = 3, "lambda"_a = 0.0005);

  // Wrappers of "dglapbuilder.h"
  _initializers.def("InitializeGpdObjects", py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&, bool const&, double const&>(&apfel::InitializeGpdObjects), "g"_a, "Thresholds"_a, "xi"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);
  _initializers.def("InitializeGpdObjectsPol", py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&, bool const&, double const&>(&apfel::InitializeGpdObjectsPol), "g"_a, "Thresholds"_a, "xi"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);
  _initializers.def("InitializeGpdObjectsTrans", py::overload_cast<apfel::Grid const&, std::vector<double> const&, double const&, bool const&, double const&>(&apfel::InitializeGpdObjectsTrans), "g"_a, "Thresholds"_a, "xi"_a, "OpEvol"_a = false, "IntEps"_a = 1e-5);

  // Wrappers of "tmdbuilder.h"
  // Define the TmdObjects structure
  py::class_<apfel::TmdObjects>(_tmd, "TmdObjects")
  .def_readwrite("Threshold", &apfel::TmdObjects::Threshold)
  .def_readwrite("Beta",      &apfel::TmdObjects::Beta)
  .def_readwrite("GammaFq",   &apfel::TmdObjects::GammaFq)
  .def_readwrite("GammaFg",   &apfel::TmdObjects::GammaFg)
  .def_readwrite("GammaK",    &apfel::TmdObjects::GammaK)
  .def_readwrite("KCS",       &apfel::TmdObjects::KCS)
  .def_readwrite("MatchingFunctionsPDFs", &apfel::TmdObjects::MatchingFunctionsPDFs)
  .def_readwrite("MatchingFunctionsFFs",  &apfel::TmdObjects::MatchingFunctionsFFs)
  .def_readwrite("HardFactors", &apfel::TmdObjects::HardFactors);

  // Add TMD object initializers
  _tmd.def("InitializeTmdObjects",            &apfel::InitializeTmdObjects,            "g"_a, "Thresholds"_a, "IntEps"_a = 1e-5, "Initialize TMD objects for evolution and matching of TMD PDFs and FFs");
  _tmd.def("InitializeTmdObjectsDYResScheme", &apfel::InitializeTmdObjectsDYResScheme, "g"_a, "Thresholds"_a, "IntEps"_a = 1e-5, "Initialize TMD objects with resummation-scheme transformation");
  _tmd.def("InitializeTmdObjectsBM",          &apfel::InitializeTmdObjectsBM,          "g"_a, "Thresholds"_a, "IntEps"_a = 1e-5, "Initialize TMD objects for Boer-Mulders TMD PDF");
  _tmd.def("InitializeTmdObjectsSivers",      &apfel::InitializeTmdObjectsSivers,      "g"_a, "Thresholds"_a, "IntEps"_a = 1e-5, "Initialize TMD objects for Sivers TMD PDF");
  _tmd.def("InitializeTmdObjectsg1",          &apfel::InitializeTmdObjectsg1,          "g"_a, "Thresholds"_a, "IntEps"_a = 1e-5, "Initialize TMD objects for g1 TMD PDFs");

  // Add TMD builders
  _tmd.def("BuildTmdPDFs", &apfel::BuildTmdPDFs, "TmdObj"_a, "CollPDFs"_a,          "Alphas"_a, "PerturbativeOrder"_a,             "Ci"_a = 1, "IntEps"_a = 1e-7, "Build matched and evolved TMD PDFs");
  _tmd.def("BuildTmdFFs",  &apfel::BuildTmdFFs,  "TmdObj"_a, "CollFFs"_a,           "Alphas"_a, "PerturbativeOrder"_a,             "Ci"_a = 1, "IntEps"_a = 1e-7, "Build matched and evolved TMD FFs");
  _tmd.def("BuildTmdJet",  &apfel::BuildTmdJet,  "TmdObj"_a, "JetAlgo"_a, "JetR"_a, "Alphas"_a, "PerturbativeOrder"_a, "CJ"_a = 1, "Ci"_a = 1, "IntEps"_a = 1e-7, "Build TMD of a jet in b-space");
  _tmd.def("MatchTmdPDFs", &apfel::MatchTmdPDFs, "TmdObj"_a, "CollPDFs"_a,          "Alphas"_a, "PerturbativeOrder"_a, "Ci"_a = 1,                                "Match TMD PDFs in b-space");
  _tmd.def("MatchTmdFFs",  &apfel::MatchTmdFFs,  "TmdObj"_a, "CollFFs"_a,           "Alphas"_a, "PerturbativeOrder"_a, "Ci"_a = 1,                                "Match TMD FFs in b-space");
  _tmd.def("MatchTmdJet",  &apfel::MatchTmdJet,  "TmdObj"_a, "JetAlgo"_a,  "tR"_a,  "Alphas"_a, "PerturbativeOrder"_a, "CJ"_a = 1, "Ci"_a = 1, "IntEps"_a = 1e-7, "Match jet TMD in b-space at initial scale");
  _tmd.def("MatchingFunctionsPDFs", &apfel::MatchingFunctionsPDFs, "TmdObj"_a, "Alphas"_a, "PerturbativeOrder"_a, "Ci"_a = 1, "Get matching functions for TMD PDFs");
  _tmd.def("MatchingFunctionsFFs",  &apfel::MatchingFunctionsFFs,  "TmdObj"_a, "Alphas"_a, "PerturbativeOrder"_a, "Ci"_a = 1, "Get matching functions for TMD FFs");

  // Add evolution factor functions
  _tmd.def("EvolutionFactors",       static_cast<std::function<std::vector<double>(double const&, double const&, double const&)>(*)(std::map<int, apfel::TmdObjects> const&, std::function<double(double const&)> const&, int const&, double const&, double const&)>(&apfel::EvolutionFactors), "TmdObj"_a, "Alphas"_a, "PerturbativeOrder"_a, "Ci"_a = 1.0, "IntEps"_a = 1e-7, "Return evolution factors for gluon and quarks as a function of (bT, mu, zeta).\n" "The returned callable accepts three doubles and returns a vector of doubles.\n" "Index 0: gluon factor, indices 1-12: quark factors."); 
  _tmd.def("EvolutionFactorsK",      &apfel::EvolutionFactorsK, "TmdObj"_a, "Alphas"_a, "PerturbativeOrder"_a, "Ci"_a = 1, "IntEps"_a = 1e-7, "Get evolution factors for gluon and quarks with double logs isolated");
  _tmd.def("QuarkEvolutionFactor",   static_cast<std::function<double(double const&, double const&, double const&)> (*)(std::map<int, apfel::TmdObjects> const&, std::function<double(double const&)> const&, int const&, double const&, double const&)>(&apfel::QuarkEvolutionFactor), "TmdObj"_a, "Alphas"_a, "PerturbativeOrder"_a, "Ci"_a = 1, "IntEps"_a = 1e-7, "Get evolution factor for quarks");
  _tmd.def("QuarkEvolutionFactorxi", &apfel::QuarkEvolutionFactorxi, "TmdObj"_a, "Alphas"_a, "PerturbativeOrder"_a, "xi"_a = 1, "Ci"_a = 1, "IntEps"_a = 1e-7, "Get evolution factor for quarks with resummation-scale parameter");
  _tmd.def("GluonEvolutionFactor",   static_cast<std::function<double(double const&, double const&, double const&)>(*)(std::map<int, apfel::TmdObjects> const&, std::function<double(double const&)> const&, int const&, double const&, double const&)>(&apfel::GluonEvolutionFactor), "TmdObj"_a, "Alphas"_a, "PerturbativeOrder"_a, "Ci"_a = 1, "IntEps"_a = 1e-7, "Get evolution factor for gluon");  

  // Add analytic evolution and kernel functions
  _tmd.def("QuarkAnalyticEvolutionFactor", &apfel::QuarkAnalyticEvolutionFactor, "TmdObj"_a, "mu0"_a, "Alphas0"_a, "kappa"_a, "kappa0"_a, "PerturbativeOrder"_a, "Get analytic evolution factor for quark TMD");
  _tmd.def("GluonAnalyticEvolutionFactor", &apfel::GluonAnalyticEvolutionFactor, "TmdObj"_a, "mu0"_a, "Alphas0"_a, "kappa"_a, "kappa0"_a, "PerturbativeOrder"_a, "Get analytic evolution factor for gluon TMD");
  _tmd.def("CollinsSoperKernel",           &apfel::CollinsSoperKernel,           "TmdObj"_a, "Alphas"_a, "PerturbativeOrder"_a, "Ci"_a = 1, "IntEps"_a = 1e-7,   "Get perturbative part of Collins-Soper kernel for quarks");
  _tmd.def("HardFactor",                   &apfel::HardFactor,      "Process"_a, "TmdObj"_a, "Alphas"_a, "PerturbativeOrder"_a, "Cf"_a = 1,                      "Get hard factor for specified process");

  // Wrappers of "ogataquadradure.h"
  py::class_<apfel::OgataQuadrature>(_ogata, "OgataQuadrature")
  .def(py::init<int const&, double const&, double const&, int const&>(), "nu"_a = 0, "CutOff"_a = 1e-5, "h"_a = 0.001, "nZeroMax"_a = 1000, "The Integrator constructor.\n Parameters:\n nu: the order of the Bessel function (default: 0)\n CutOff: the accuracy computed as a cutoff on the size of the last computed term relative to the total (default: 10^-5)\n h: internal variable of the algorithm (default: 0.001)\n nZeroMax: maximum number of terms in the Ogata quadrature (default: 1000)")
  
  .def("InitialiseWeights", &apfel::OgataQuadrature::InitialiseWeights, "nZeroMax"_a, "Function that initialises the coordinates and the weights for the Ogata quadrature.")
  // transform is a templated function in cpp. Here the tamplate instantiation is explicitly chosen to be double.
  .def("transform", [](const apfel::OgataQuadrature &self, std::function<double(double const&)> const& func, double const& qT, bool const& Dynh, int const& nmax, int const& period) {return self.transform<double>(func, qT, Dynh, nmax, period);}, "func"_a, "qT"_a, "Dynh"_a = true, "nmax"_a = 1000, "period"_a = 10, "Function that transform the input function.\n" "Parameters:\n" "  func: function to be transformed\n" "  qT: value of qT in which to compute the transform\n" "  Dynh: switch to compute the step parameter _h dynamically (default: true)\n" "  nmax: maximum number of terms in the Ogata quadrature (default: 1000)\n" "  period: interval across which the integral is checked to determine where the sum is truncated (default: 10)")
  .def("GetBesselOrder",    &apfel::OgataQuadrature::GetBesselOrder,    "Function that returns the Bessel order.")
  .def("GetCutOff",         &apfel::OgataQuadrature::GetCutOff,         "Function that returns the Ogata cut-off parameter.")
  .def("GetStepParameter",  &apfel::OgataQuadrature::GetStepParameter,  "Function that returns the Ogata step parameter.")
  .def("GetCoordinates",    &apfel::OgataQuadrature::GetCoordinates,    "Function that returns the unscaled coordinates used in the Ogata quadrature.")
  .def("GetWeights",        &apfel::OgataQuadrature::GetWeights,        "Function that returns the weights used in the Ogata quadrature.")
  .def("SetStepParameter",  &apfel::OgataQuadrature::SetStepParameter,  "h"_a,  "Function that sets the Ogata step parameter.")
  .def("JnuZerosGenerator", &apfel::OgataQuadrature::JnuZerosGenerator, "nu"_a, "Function that writes on screen the first 1000 zeros of the Bessel function J0.");

}