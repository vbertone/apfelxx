#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>
#include <apfel/apfelxx.h>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(apfelpy, m) {
  // Documentation
  m.doc() = "Python wrapper of APFEL++";

  // Wrappers of "messages.h"
  m.def("SetVerbosityLevel", &apfel::SetVerbosityLevel, "vl"_a);
  m.def("GetVerbosityLevel", &apfel::GetVerbosityLevel);
  m.def("Banner",            &apfel::Banner);

  // Wrappers of "constants.h"
  py::module_ _constants = m.def_submodule("constants", "Numerical constants");
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
  _constants.attr("ConvFact")   = apfel::ConvFact;
  _constants.attr("ZMass")      = apfel::ZMass;
  _constants.attr("GammaZ")     = apfel::GammaZ;
  _constants.attr("WMass")      = apfel::WMass;
  _constants.attr("GammaW")     = apfel::GammaW;
  _constants.attr("ProtonMass") = apfel::ProtonMass;
  _constants.attr("Sin2ThetaW") = apfel::Sin2ThetaW;
  _constants.attr("GFermi")     = apfel::GFermi;
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

  // Utility functions
  py::module_ _utilities = m.def_submodule("utilities", "Utility functions");

  // Wrappers of "lhtoypdfs.h"
  _utilities.def("xupv",      &apfel::xupv,      "x"_a);
  _utilities.def("xdnv",      &apfel::xdnv,      "x"_a);
  _utilities.def("xglu",      &apfel::xglu,      "x"_a);
  _utilities.def("xdbar",     &apfel::xdbar,     "x"_a);
  _utilities.def("xubar",     &apfel::xubar,     "x"_a);
  _utilities.def("xsbar",     &apfel::xsbar,     "x"_a);
  _utilities.def("LHToyPDFs", &apfel::LHToyPDFs, "x"_a, "Q"_a);

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
  _utilities.def("PhysToQCDEv", &apfel::PhysToQCDEv,              "InPhysMap"_a);
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

  // Wrappers form "interpolator.h"
  // Trampoline class for virtual class
  class PyInterpolator: public apfel::Interpolator
  {
  public:
    using Interpolator::Interpolator;
    double InterpolantLog(int const& beta, double const& lnx, apfel::SubGrid const& sg) const override { PYBIND11_OVERRIDE_PURE(double, Interpolator, InterpolantLog, beta, lnx, sg); };
    double Interpolant(int const& beta, double const& x, apfel::SubGrid const& sg) const override { PYBIND11_OVERRIDE_PURE(double, Interpolator, Interpolant, beta, x, sg); };
    std::array<int, 2> SumBounds(double const& x, apfel::SubGrid const& sg) const override { PYBIND11_OVERRIDE_PURE(PYBIND11_TYPE(std::array<int, 2>), Interpolator, SumBounds, x, sg); };
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

  // Wrappers form "lagrangeinterpolator.h"
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
    .def("SetJointGrid", &apfel::Distribution::SetJointGrid, "ix"_a, "x"_a)
    .def("SetSubGrid", &apfel::Distribution::SetSubGrid, "ig"_a, "ix"_a, "x"_a)
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
    double Regular(double const& x) const override { PYBIND11_OVERRIDE(double, Expression, Regular, x); };
    double Singular(double const& x) const override { PYBIND11_OVERRIDE(double, Expression, Singular, x); };
    double Local(double const& x) const override { PYBIND11_OVERRIDE(double, Expression, Local, x); };
    double LocalPV(double const& x) const override { PYBIND11_OVERRIDE(double, Expression, LocalPV, x); };
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
  py::class_<apfel::matrix<size_t>>(m, "matrixSize_t")
    .def(py::init<size_t const&, size_t const&>(), "row"_a = 0, "col"_a = 0)
    .def("resize", &apfel::matrix<size_t>::resize, "row"_a, "col"_a, "v"_a = 0)
    .def("set", &apfel::matrix<size_t>::set, "v"_a)
    .def("size", &apfel::matrix<size_t>::size, "dim"_a)
    .def("__call__", [] (apfel::matrix<size_t>& c, size_t const& i, size_t const& j){ return c(i, j); })
    .def("__call__", [] (apfel::matrix<size_t> const& c, size_t const& i, size_t const& j){ return c(i, j); });

  py::class_<apfel::matrix<int>>(m, "matrixInt")
    .def(py::init<size_t const&, size_t const&>(), "row"_a = 0, "col"_a = 0)
    .def("resize", &apfel::matrix<int>::resize, "row"_a, "col"_a, "v"_a = 0)
    .def("set", &apfel::matrix<int>::set, "v"_a)
    .def("size", &apfel::matrix<int>::size, "dim"_a)
    .def("__call__", [] (apfel::matrix<int>& c, size_t const& i, size_t const& j){ return c(i, j); })
    .def("__call__", [] (apfel::matrix<int> const& c, size_t const& i, size_t const& j){ return c(i, j); });

  py::class_<apfel::matrix<float>>(m, "matrixFloat")
    .def(py::init<size_t const&, size_t const&>(), "row"_a = 0, "col"_a = 0)
    .def("resize", &apfel::matrix<float>::resize, "row"_a, "col"_a, "v"_a = 0)
    .def("set", &apfel::matrix<float>::set, "v"_a)
    .def("size", &apfel::matrix<float>::size, "dim"_a)
    .def("__call__", [] (apfel::matrix<float>& c, size_t const& i, size_t const& j){ return c(i, j); })
    .def("__call__", [] (apfel::matrix<float> const& c, size_t const& i, size_t const& j){ return c(i, j); });

  py::class_<apfel::matrix<double>>(m, "matrixDouble")
    .def(py::init<size_t const&, size_t const&>(), "row"_a = 0, "col"_a = 0)
    .def("resize", &apfel::matrix<double>::resize, "row"_a, "col"_a, "v"_a = 0)
    .def("set", &apfel::matrix<double>::set, "v"_a)
    .def("size", &apfel::matrix<double>::size, "dim"_a)
    .def("__call__", [] (apfel::matrix<double>& c, size_t const& i, size_t const& j){ return c(i, j); })
    .def("__call__", [] (apfel::matrix<double> const& c, size_t const& i, size_t const& j){ return c(i, j); });

  py::class_<apfel::matrix<std::vector<int>>>(m, "matrixVectorInt")
    .def(py::init<size_t const&, size_t const&>(), "row"_a = 0, "col"_a = 0)
    .def("resize", &apfel::matrix<std::vector<int>>::resize, "row"_a, "col"_a, "v"_a = 0)
    .def("set", &apfel::matrix<std::vector<int>>::set, "v"_a)
    .def("size", &apfel::matrix<std::vector<int>>::size, "dim"_a)
    .def("__call__", [] (apfel::matrix<std::vector<int>>& c, size_t const& i, size_t const& j){ return c(i, j); })
    .def("__call__", [] (apfel::matrix<std::vector<int>> const& c, size_t const& i, size_t const& j){ return c(i, j); });

  py::class_<apfel::matrix<std::vector<double>>>(m, "matrixVectorDouble")
    .def(py::init<size_t const&, size_t const&>(), "row"_a = 0, "col"_a = 0)
    .def("resize", &apfel::matrix<std::vector<double>>::resize, "row"_a, "col"_a, "v"_a = 0)
    .def("set", &apfel::matrix<std::vector<double>>::set, "v"_a)
    .def("size", &apfel::matrix<std::vector<double>>::size, "dim"_a)
    .def("__call__", [] (apfel::matrix<std::vector<double>>& c, size_t const& i, size_t const& j){ return c(i, j); })
    .def("__call__", [] (apfel::matrix<std::vector<double>> const& c, size_t const& i, size_t const& j){ return c(i, j); });

  // Wrappers of "operator.h"
  py::class_<apfel::Operator>(m, "Operator")
    .def(py::init<apfel::Operator const&>(), "g"_a)
    .def(py::init<apfel::Grid const&>(), "g"_a)
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
    .def(py::self *= py::self)
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
    .def("Evaluate", &apfel::DoubleObject<apfel::Operator>::Evaluate, "x"_a, "z"_a)
    .def("Evaluate1", &apfel::DoubleObject<apfel::Operator>::Evaluate1, "x"_a)
    .def("Evaluate2", &apfel::DoubleObject<apfel::Operator>::Evaluate2, "z"_a)
    .def("Derive", &apfel::DoubleObject<apfel::Operator>::Derive, "x"_a, "z"_a)
    .def("Derive1", &apfel::DoubleObject<apfel::Operator>::Derive1, "x"_a)
    .def("Derive2", &apfel::DoubleObject<apfel::Operator>::Derive2, "z"_a)
    .def("Integrate", py::overload_cast<double const&, double const&, double const&, double const&>(&apfel::DoubleObject<apfel::Operator>::Integrate, py::const_), "xl"_a, "xu"_a, "zl"_a, "zu"_a)
    .def("Integrate1", &apfel::DoubleObject<apfel::Operator>::Integrate1, "xl"_a, "xu"_a)
    .def("Integrate2", &apfel::DoubleObject<apfel::Operator>::Integrate2, "zl"_a, "zu"_a)
    .def("Integrate", py::overload_cast<double const&, double const&, std::function<double(double const&)>, std::function<double(double const&)>>(&apfel::DoubleObject<apfel::Operator>::Integrate, py::const_), "xl"_a, "xu"_a, "zlx"_a, "zux"_a)
    .def("Integrate", py::overload_cast<std::function<double(double const&)>, std::function<double(double const&)>, double const&, double const&>(&apfel::DoubleObject<apfel::Operator>::Integrate, py::const_), "xlz"_a, "xuz"_a, "zl"_a, "zu"_a)
    .def("MultiplyBy", &apfel::DoubleObject<apfel::Operator>::MultiplyBy, "fx"_a, "fz"_a)
    .def("Print", &apfel::DoubleObject<apfel::Operator>::Print)
    .def(py::self *= double())
    .def(py::self *= py::self)
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
    .def("Evaluate", &apfel::DoubleObject<apfel::Distribution, apfel::Operator>::Evaluate, "x"_a, "z"_a)
    .def("Evaluate1", &apfel::DoubleObject<apfel::Distribution, apfel::Operator>::Evaluate1, "x"_a)
    .def("Evaluate2", &apfel::DoubleObject<apfel::Distribution, apfel::Operator>::Evaluate2, "z"_a)
    .def("Derive", &apfel::DoubleObject<apfel::Distribution, apfel::Operator>::Derive, "x"_a, "z"_a)
    .def("Derive1", &apfel::DoubleObject<apfel::Distribution, apfel::Operator>::Derive1, "x"_a)
    .def("Derive2", &apfel::DoubleObject<apfel::Distribution, apfel::Operator>::Derive2, "z"_a)
    .def("Integrate", py::overload_cast<double const&, double const&, double const&, double const&>(&apfel::DoubleObject<apfel::Distribution, apfel::Operator>::Integrate, py::const_), "xl"_a, "xu"_a, "zl"_a, "zu"_a)
    .def("Integrate1", &apfel::DoubleObject<apfel::Distribution, apfel::Operator>::Integrate1, "xl"_a, "xu"_a)
    .def("Integrate2", &apfel::DoubleObject<apfel::Distribution, apfel::Operator>::Integrate2, "zl"_a, "zu"_a)
    .def("Integrate", py::overload_cast<double const&, double const&, std::function<double(double const&)>, std::function<double(double const&)>>(&apfel::DoubleObject<apfel::Distribution, apfel::Operator>::Integrate, py::const_), "xl"_a, "xu"_a, "zlx"_a, "zux"_a)
    .def("Integrate", py::overload_cast<std::function<double(double const&)>, std::function<double(double const&)>, double const&, double const&>(&apfel::DoubleObject<apfel::Distribution, apfel::Operator>::Integrate, py::const_), "xlz"_a, "xuz"_a, "zl"_a, "zu"_a)
    .def("MultiplyBy", &apfel::DoubleObject<apfel::Distribution, apfel::Operator>::MultiplyBy, "fx"_a, "fz"_a)
    .def("Print", &apfel::DoubleObject<apfel::Distribution, apfel::Operator>::Print)
    .def(py::self *= double())
    .def(py::self *= py::self)
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
    .def("Evaluate", &apfel::DoubleObject<apfel::Operator, apfel::Distribution>::Evaluate, "x"_a, "z"_a)
    .def("Evaluate1", &apfel::DoubleObject<apfel::Operator, apfel::Distribution>::Evaluate1, "x"_a)
    .def("Evaluate2", &apfel::DoubleObject<apfel::Operator, apfel::Distribution>::Evaluate2, "z"_a)
    .def("Derive", &apfel::DoubleObject<apfel::Operator, apfel::Distribution>::Derive, "x"_a, "z"_a)
    .def("Derive1", &apfel::DoubleObject<apfel::Operator, apfel::Distribution>::Derive1, "x"_a)
    .def("Derive2", &apfel::DoubleObject<apfel::Operator, apfel::Distribution>::Derive2, "z"_a)
    .def("Integrate", py::overload_cast<double const&, double const&, double const&, double const&>(&apfel::DoubleObject<apfel::Operator, apfel::Distribution>::Integrate, py::const_), "xl"_a, "xu"_a, "zl"_a, "zu"_a)
    .def("Integrate1", &apfel::DoubleObject<apfel::Operator, apfel::Distribution>::Integrate1, "xl"_a, "xu"_a)
    .def("Integrate2", &apfel::DoubleObject<apfel::Operator, apfel::Distribution>::Integrate2, "zl"_a, "zu"_a)
    .def("Integrate", py::overload_cast<double const&, double const&, std::function<double(double const&)>, std::function<double(double const&)>>(&apfel::DoubleObject<apfel::Operator, apfel::Distribution>::Integrate, py::const_), "xl"_a, "xu"_a, "zlx"_a, "zux"_a)
    .def("Integrate", py::overload_cast<std::function<double(double const&)>, std::function<double(double const&)>, double const&, double const&>(&apfel::DoubleObject<apfel::Operator, apfel::Distribution>::Integrate, py::const_), "xlz"_a, "xuz"_a, "zl"_a, "zu"_a)
    .def("MultiplyBy", &apfel::DoubleObject<apfel::Operator, apfel::Distribution>::MultiplyBy, "fx"_a, "fz"_a)
    .def("Print", &apfel::DoubleObject<apfel::Operator, apfel::Distribution>::Print)
    .def(py::self *= double())
    .def(py::self *= py::self)
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

  // Wrappers of "convolutionmap.h"

  // Wrappers of "set.h"
  
}
