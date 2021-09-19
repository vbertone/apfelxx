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
  _utilities.def("PhysToQCDEv", py::overload_cast<std::map<int, double> const&>(&apfel::PhysToQCDEv),              "InPhysMap"_a);
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
  py::class_<apfel::Expression>(m, "Expression")
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

  // Wrappers of "matrix.h"
  py::class_<apfel::matrix<T>>(m, "matrix");
  //.def(py::init<apfel::Operator const&>(), "g"_a)

/*  
    matrix(size_t const& row = 0, size_t const& col = 0);
    void resize(size_t const& row, size_t const& col, T const& v = 0);
    void set(T const& v);
    size_t const& size(size_t const& dim) const { return _size[dim]; }
    T&       operator()(size_t const& i, size_t const& j)       { return _data[i*_size[1]+j]; }
    T const& operator()(size_t const& i, size_t const& j) const { return _data[i*_size[1]+j]; }
*/

  // Wrappers of "operator.h"
  py::class_<apfel::Operator>(m, "Operator")
    .def(py::init<apfel::Operator const&>(), "g"_a)
    .def(py::init<apfel::Grid const&>(), "g"_a)
    .def(py::init<apfel::Grid const&, apfel::Expression const&, double const&>(), "g"_a, "expr"_a, "eps"_a)
;
/*
    Operator(Operator const&) = default;
    Operator(Grid const& gr);
    Operator(Grid const& gr, Expression const& expr, double const& eps = 1e-5);
    Distribution operator *= (Distribution const& d) const;         //!< this *= Distribution
    Operator& operator *= (Operator const& o);                      //!< this *= Operator
    Operator& operator  = (Operator const& o);                      //!< this  = Operator
    Operator& operator *= (double const& s);                        //!< this *= Scalar
    Operator& operator *= (std::function<double(double const&)> f); //!< This *= Function
    Operator& operator /= (double const& s);                        //!< this /= Scalar
    Operator& operator += (Operator const& o);                      //!< this += Operator
    Operator& operator -= (Operator const& o);                      //!< this -= Operator
    Grid const& GetGrid() const { return _grid; }
    std::vector<matrix<double>> GetOperator() const { return _Operator; }
  Distribution operator * (Operator lhs, Distribution const& rhs);                //!< Operator*Distribution
  Operator     operator * (Operator lhs, Operator const& rhs);                    //!< Operator*Operator
  Operator     operator * (double const& s, Operator rhs);                        //!< Scalar*Operator
  Operator     operator * (Operator lhs, double const& s);                        //!< Operator*Scalar
  Operator     operator * (std::function<double(double const&)> f, Operator rhs); //!< function*Operator
  Operator     operator * (Operator lhs, std::function<double(double const&)> f); //!< Operator*function
  Operator     operator / (Operator lhs, double const& s);                        //!< Operator/Scalar
  Operator     operator + (Operator lhs, Operator const& rhs);                    //!< Operator+Operator
  Operator     operator - (Operator lhs, Operator const& rhs);                    //!< Operator-Operator
*/
}
