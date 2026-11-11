//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>
#include <apfel/zeromasscff.h>

#include <iomanip>

int main()
{
  apfel::Banner();

  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{30, 1e-5, 3}, apfel::SubGrid{20, 1e-1, 3}, apfel::SubGrid{10, 9e-1, 3}}};

  // Initial scale
  const double mu0 = sqrt(2);

  // Vector of thresholds
  const std::vector<double> Thresholds = {0, 0, 0, sqrt(2), 4.5, 175};

  // Perturbative order of the CFFs (GPD evolution is only available
  // at LO)
  const int PerturbativeOrder = 1;

  // Skewness
  const double xi = 0.27;

  // Running coupling
  const double AlphaQCDRef = 0.35;
  const double MuAlphaQCDRef = sqrt(2);
  apfel::AlphaQCD a{AlphaQCDRef, MuAlphaQCDRef, Thresholds, PerturbativeOrder};
  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // Effective charges
  std::function<std::vector<double>(double const&)> fBq = [=] (double const&) -> std::vector<double> { return apfel::QCh2; };

  // Initialize GPD evolution objects
  const auto GpdObj = InitializeGpdObjects(g, Thresholds, xi);

  // Define set of toy GPDs
  const auto ToyGPDs = [=] (double const& x, double const&) -> std::map<int, double>
  {
    const double gl  = x * ( 1 - x * x ) / x;
    const double up  = x * x * ( 1 - x * x );
    const double upb = - x * x * x * ( 1 - x * x );
    const double dn  = x * ( 1 - x * x ) * std::exp(- pow(x - xi, 2));
    const double dnb = - x * ( 1 - x * x ) * std::exp(- pow(- x - xi, 2));

    // Construct QCD evolution basis conbinations
    double const Gluon   = gl;
    double const Singlet = dn + dnb + up + upb;
    double const T3      = up + upb - dn - dnb;
    double const Valence = dn - dnb + up - upb;
    double const V3      = up - upb - dn + dnb;

    // Fill in map in the QCD evolution basis
    std::map<int, double> QCDEvMap;
    QCDEvMap[0]  = Gluon;
    QCDEvMap[1]  = Singlet;
    QCDEvMap[2]  = Valence;
    QCDEvMap[3]  = T3;
    QCDEvMap[4]  = V3;
    QCDEvMap[5]  = Singlet;
    QCDEvMap[6]  = Valence;
    QCDEvMap[7]  = Singlet;
    QCDEvMap[8]  = Valence;
    QCDEvMap[9]  = Singlet;
    QCDEvMap[10] = Valence;
    QCDEvMap[11] = Singlet;
    QCDEvMap[12] = Valence;

    return QCDEvMap;
  };

  // Construct the DGLAP objects
  const auto EvolvedGPDs = BuildDglap(GpdObj, ToyGPDs, mu0, 0, as);

  // Tabulate GPDs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedGPDs{*EvolvedGPDs, 30, 1, 100, 3};

  // Evolved GPDs
  const auto GPDs = [&] (double const& x, double const& Q) -> std::map<int, double> { return TabulatedGPDs.EvaluateMapxQ(x, Q); };

  // Initialize coefficient functions for CFFs
  const auto ImF1SLObj = InitializeImCFF1NCObjectsZMSpaceLike(g, Thresholds, xi);
  const auto ReF1SLObj = InitializeReCFF1NCObjectsZMSpaceLike(g, Thresholds, xi);
  const auto ImFLSLObj = InitializeImCFFLNCObjectsZMSpaceLike(g, Thresholds, xi);
  const auto ReFLSLObj = InitializeReCFFLNCObjectsZMSpaceLike(g, Thresholds, xi);
  const auto ImF1TLObj = InitializeImCFF1NCObjectsZMTimeLike(g, Thresholds, xi);
  const auto ReF1TLObj = InitializeReCFF1NCObjectsZMTimeLike(g, Thresholds, xi);
  const auto ImFLTLObj = InitializeImCFFLNCObjectsZMTimeLike(g, Thresholds, xi);
  const auto ReFLTLObj = InitializeReCFFLNCObjectsZMTimeLike(g, Thresholds, xi);

  // Initialize form factors
  const auto ImF1SL = BuildStructureFunctions(ImF1SLObj, GPDs, PerturbativeOrder, as, fBq);
  const auto ReF1SL = BuildStructureFunctions(ReF1SLObj, GPDs, PerturbativeOrder, as, fBq);
  const auto ImFLSL = BuildStructureFunctions(ImFLSLObj, GPDs, PerturbativeOrder, as, fBq);
  const auto ReFLSL = BuildStructureFunctions(ReFLSLObj, GPDs, PerturbativeOrder, as, fBq);
  const auto ImF1TL = BuildStructureFunctions(ImF1TLObj, GPDs, PerturbativeOrder, as, fBq);
  const auto ReF1TL = BuildStructureFunctions(ReF1TLObj, GPDs, PerturbativeOrder, as, fBq);
  const auto ImFLTL = BuildStructureFunctions(ImFLTLObj, GPDs, PerturbativeOrder, as, fBq);
  const auto ReFLTL = BuildStructureFunctions(ReFLTLObj, GPDs, PerturbativeOrder, as, fBq);

  // Tabulate Structure functions
  const apfel::TabulateObject<apfel::Distribution> ImF1SLtotal {[&] (double const& Q) -> apfel::Distribution{ return ImF1SL.at(0).Evaluate(Q); }, 50, 1, 99, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> ReF1SLtotal {[&] (double const& Q) -> apfel::Distribution{ return ReF1SL.at(0).Evaluate(Q); }, 50, 1, 99, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> ImFLSLtotal {[&] (double const& Q) -> apfel::Distribution{ return ImFLSL.at(0).Evaluate(Q); }, 50, 1, 99, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> ReFLSLtotal {[&] (double const& Q) -> apfel::Distribution{ return ReFLSL.at(0).Evaluate(Q); }, 50, 1, 99, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> ImF1TLtotal {[&] (double const& Q) -> apfel::Distribution{ return ImF1TL.at(0).Evaluate(Q); }, 50, 1, 99, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> ReF1TLtotal {[&] (double const& Q) -> apfel::Distribution{ return ReF1TL.at(0).Evaluate(Q); }, 50, 1, 99, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> ImFLTLtotal {[&] (double const& Q) -> apfel::Distribution{ return ImFLTL.at(0).Evaluate(Q); }, 50, 1, 99, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> ReFLTLtotal {[&] (double const& Q) -> apfel::Distribution{ return ReFLTL.at(0).Evaluate(Q); }, 50, 1, 99, 3, Thresholds};

  apfel::Timer t;

  // Final scale
  const double Q = 10;

  std::cout << std::scientific << std::endl;
  std::cout << "xi = " << xi << std::endl;
  std::cout << "Alphas(Q) = " << as(Q) << std::endl;
  std::cout << std::endl;

  // Print results
  const std::vector<double> xlha = {1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};

  std::cout << "Q2 > 0:" << std::endl;
  std::cout << "     x    "
            << "     ReF1tot    "
            << "     ImF1tot    "
            << "     ReFLtot    "
            << "     ImFLtot    "
            << std::endl;
  for (double const& x : xlha)
    std::cout << std::setprecision(4) << x << "  " << std::setprecision(8)
              << ReF1SLtotal.EvaluatexQ(x, Q) / x << "  "
              << ImF1SLtotal.EvaluatexQ(x, Q) / x << "  "
              << ReFLSLtotal.EvaluatexQ(x, Q) / x << "  "
              << ImFLSLtotal.EvaluatexQ(x, Q) / x << "  "
              << std::endl;
  std::cout << std::endl;

  std::cout << "Q2 < 0:" << std::endl;
  std::cout << "     x    "
            << "     ReF1tot    "
            << "     ImF1tot    "
            << "     ReFLtot    "
            << "     ImFLtot    "
            << std::endl;
  for (double const& x : xlha)
    std::cout << std::setprecision(4) << x << "  " << std::setprecision(8)
              << ReF1TLtotal.EvaluatexQ(x, Q) / x << "  "
              << ImF1TLtotal.EvaluatexQ(x, Q) / x << "  "
              << ReFLTLtotal.EvaluatexQ(x, Q) / x << "  "
              << ImFLTLtotal.EvaluatexQ(x, Q) / x << "  "
              << std::endl;
  std::cout << std::endl;

  // Check against DIS at xi = 0, where Im F1 = pi (F2 - FL) and Im
  // FL = 2 pi FL, using the same distributions
  const auto ImF1SL0Obj = InitializeImCFF1NCObjectsZMSpaceLike(g, Thresholds, 0);
  const auto ImFLSL0Obj = InitializeImCFFLNCObjectsZMSpaceLike(g, Thresholds, 0);
  const auto F2Obj      = InitializeF2NCObjectsZM(g, Thresholds);
  const auto FLObj      = InitializeFLNCObjectsZM(g, Thresholds);
  const apfel::Distribution ImF1SL0 = BuildStructureFunctions(ImF1SL0Obj, ToyGPDs, PerturbativeOrder, as, fBq).at(0).Evaluate(Q);
  const apfel::Distribution ImFLSL0 = BuildStructureFunctions(ImFLSL0Obj, ToyGPDs, PerturbativeOrder, as, fBq).at(0).Evaluate(Q);
  const apfel::Distribution F2      = BuildStructureFunctions(F2Obj,      ToyGPDs, PerturbativeOrder, as, fBq).at(0).Evaluate(Q);
  const apfel::Distribution FL      = BuildStructureFunctions(FLObj,      ToyGPDs, PerturbativeOrder, as, fBq).at(0).Evaluate(Q);

  std::cout << "Q2 > 0, xi = 0 vs DIS:" << std::endl;
  std::cout << "     x    "
            << " ImF1/pi(F2-FL) "
            << "   ImFL/2piFL   "
            << std::endl;
  for (double const& x : xlha)
    std::cout << std::setprecision(4) << x << "  " << std::setprecision(8)
              << ImF1SL0.Evaluate(x) / ( M_PI * ( F2.Evaluate(x) - FL.Evaluate(x) ) ) << "  "
              << ImFLSL0.Evaluate(x) / ( 2 * M_PI * FL.Evaluate(x) ) << "  "
              << std::endl;
  std::cout << std::endl;

  // DVCS at x = xi, with hard scale Q2 = - q1^2
  const apfel::Distribution ImF1DVCS = BuildStructureFunctions(InitializeImCFF1NCObjectsZMDVCS(g, Thresholds, xi), GPDs, PerturbativeOrder, as, fBq).at(0).Evaluate(Q);
  const apfel::Distribution ReF1DVCS = BuildStructureFunctions(InitializeReCFF1NCObjectsZMDVCS(g, Thresholds, xi), GPDs, PerturbativeOrder, as, fBq).at(0).Evaluate(Q);

  std::cout << "DVCS (x = xi):" << std::endl;
  std::cout << "     ReF1tot    "
            << "     ImF1tot    "
            << std::endl;
  std::cout << std::setprecision(8)
            << ReF1DVCS.Evaluate(xi) / xi << "  "
            << ImF1DVCS.Evaluate(xi) / xi << "  "
            << std::endl;
  std::cout << std::endl;

  // Compare the NLO DVCS kernels with the kappa -> 1 limit of the
  // DDVCS. The DDVCS kernels assume mu2 = Q2 / 2, hence the shift -
  // ln(2) C0 x P0.
  std::vector<double> nodes(g.GetJointGrid().GetGrid().begin(), g.GetJointGrid().GetGrid().begin() + g.GetJointGrid().nx() + 1);
  nodes.push_back(xi);
  std::sort(nodes.begin(), nodes.end());
  const apfel::Grid gx{{apfel::SubGrid{nodes, 3}}};
  const double xid = xi * ( 1 + apfel::eps6 );
  const apfel::Distribution qd{gx, [&] (double const& x) -> double { return ToyGPDs(x, 0).at(1); }};
  const apfel::Distribution gd{gx, [&] (double const& x) -> double { return ToyGPDs(x, 0).at(0); }};
  const auto P0 = InitializeGpdObjects(gx, Thresholds, xi).at(3).SplittingFunctions.at(0).GetObjects();
  const apfel::Distribution P0q = P0.at(apfel::EvolutionBasisQCD::PNSP) * qd;
  const apfel::Distribution P0g = ( 1. / 3 ) * ( P0.at(apfel::EvolutionBasisQCD::PQG) * gd );
  const apfel::Operator ReC0{gx, apfel::ReCFF10ns{}, apfel::eps5, true};
  const std::vector<double> ReqDVCS  = {( apfel::Operator{gx, apfel::ReCFF11nsDVCS{}, apfel::eps5, true} * qd ).Evaluate(xi) / xi,
                                        ( apfel::Operator{gx, apfel::ReCFF11nsSpaceLike{xid}, apfel::eps5, true} * qd ).Evaluate(xi) / xi - log(2) * ( ReC0 * P0q ).Evaluate(xi) / xi
                                       };
  const std::vector<double> ImqDVCS  = {( apfel::Operator{gx, apfel::ImCFF11nsDVCS{}, apfel::eps5, true} * qd ).Evaluate(xi) / xi,
                                        ( apfel::Operator{gx, apfel::ImCFF11nsSpaceLike{xid}, apfel::eps5, true} * qd ).Evaluate(xi) / xi - log(2) * M_PI * P0q.Evaluate(xi) / xi
                                       };
  const std::vector<double> RegDVCS  = {( apfel::Operator{gx, apfel::ReCFF11gDVCS{}, apfel::eps5, true} * gd ).Evaluate(xi) / xi,
                                        ( apfel::Operator{gx, apfel::ReCFF11gSpaceLike{xid}, apfel::eps5, true} * gd ).Evaluate(xi) / xi - log(2) * ( ReC0 * P0g ).Evaluate(xi) / xi
                                       };
  const std::vector<double> ImgDVCS  = {( apfel::Operator{gx, apfel::ImCFF11gDVCS{}, apfel::eps5, true} * gd ).Evaluate(xi) / xi,
                                        ( apfel::Operator{gx, apfel::ImCFF11gSpaceLike{xid}, apfel::eps5, true} * gd ).Evaluate(xi) / xi - log(2) * M_PI * P0g.Evaluate(xi) / xi
                                       };

  std::cout << "NLO DVCS vs kappa -> 1 limit of DDVCS (x = xi):" << std::endl;
  std::cout << "               DVCS          DDVCS(k->1)       ratio" << std::endl;
  std::cout << std::setprecision(8)
            << "Re quark: " << ReqDVCS[0] << "  " << ReqDVCS[1] << "  " << ReqDVCS[0] / ReqDVCS[1] << std::endl
            << "Im quark: " << ImqDVCS[0] << "  " << ImqDVCS[1] << "  " << ImqDVCS[0] / ImqDVCS[1] << std::endl
            << "Re gluon: " << RegDVCS[0] << "  " << RegDVCS[1] << "  " << RegDVCS[0] / RegDVCS[1] << std::endl
            << "Im gluon: " << ImgDVCS[0] << "  " << ImgDVCS[1] << "  " << ImgDVCS[0] / ImgDVCS[1] << std::endl;
  std::cout << std::endl;

  // Compare F1 at x = xi from the DVCS and DDVCS initialisers. Since
  // the DDVCS hard scale is - q1^2 / 2, DDVCS is computed at Q /
  // sqrt(2) with renormalisation and factorisation scales sqrt(2)
  // times larger, i.e. Q, so that the two coincide at fixed order.
  const apfel::Distribution ReF1DVCSx  = BuildStructureFunctions(InitializeReCFF1NCObjectsZMDVCS(gx, Thresholds, xi), GPDs, PerturbativeOrder, as, fBq).at(0).Evaluate(Q);
  const apfel::Distribution ImF1DVCSx  = BuildStructureFunctions(InitializeImCFF1NCObjectsZMDVCS(gx, Thresholds, xi), GPDs, PerturbativeOrder, as, fBq).at(0).Evaluate(Q);
  const apfel::Distribution ReF1DDVCSx = BuildStructureFunctions(InitializeReCFF1NCObjectsZMSpaceLike(gx, Thresholds, xi), GPDs, PerturbativeOrder, as, fBq, sqrt(2), sqrt(2)).at(0).Evaluate(Q / sqrt(2));
  const apfel::Distribution ImF1DDVCSx = BuildStructureFunctions(InitializeImCFF1NCObjectsZMSpaceLike(gx, Thresholds, xi), GPDs, PerturbativeOrder, as, fBq, sqrt(2), sqrt(2)).at(0).Evaluate(Q / sqrt(2));

  std::cout << "F1 from DVCS vs DDVCS (x = xi):" << std::endl;
  std::cout << "               DVCS          DDVCS(k->1)       ratio" << std::endl;
  std::cout << std::setprecision(8)
            << "Re F1:    " << ReF1DVCSx.Evaluate(xi) / xi << "  " << ReF1DDVCSx.Evaluate(xi) / xi << "  " << ReF1DVCSx.Evaluate(xi) / ReF1DDVCSx.Evaluate(xi) << std::endl
            << "Im F1:    " << ImF1DVCSx.Evaluate(xi) / xi << "  " << ImF1DDVCSx.Evaluate(xi) / xi << "  " << ImF1DVCSx.Evaluate(xi) / ImF1DDVCSx.Evaluate(xi) << std::endl;
  std::cout << std::endl;

  t.stop();

  const int k = 1000000;
  std::cout << "Interpolating " << k << " times F1 on the grid... ";
  t.start();
  for (int i = 0; i < k; i++)
    ReF1SLtotal.EvaluatexQ(0.05, Q);
  t.stop();

  return 0;
}
