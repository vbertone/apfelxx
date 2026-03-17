// Test for `BuildSidisUnpFT` and `BuildSidisG1`.
//
// Uses the new `InitializeSidisObjects` + Build* API from `sidisbuilder.h`.
// This test uses toy PDFs/FFs and ensures internal consistency.

#include <apfel/apfelxx.h>

#include <iomanip>
#include <iostream>
#include <cmath>
#include <map>
#include <string>

int main()
{
  apfel::SetVerbosityLevel(1);

  // Use a very coarse grid for a fast test.
  const apfel::Grid gx{{apfel::SubGrid{5, 1e-3, 3}}};
  const apfel::Grid gz{{apfel::SubGrid{5, 1e-2, 3}}};

  // Heavy-quark thresholds
  const std::vector<double> Thresholds{0, 0, 0, 100, 200};

  // Running coupling (NLO)
  const double as_ref = 0.35;
  const double mu_ref = sqrt(2);
  apfel::AlphaQCD as_obj{as_ref, mu_ref, Thresholds, 1};
  const auto as = [&] (double const& mu)
  {
    return as_obj.Evaluate(mu);
  };

  // Toy PDF / FF sets (at mu_ref)
  const auto InPDFs  = [&] (double const& x, double const& Q)
  {
    return apfel::QCDEvToPhys(apfel::LHToyPDFs(x, Q));
  };
  const auto InpPDFs = [&] (double const& x, double const& Q)
  {
    return apfel::QCDEvToPhys(apfel::LHToyPDFsPol(x, Q));
  };
  const auto InFFs   = [&] (double const& x, double const& Q)
  {
    return apfel::QCDEvToPhys(apfel::LHToyFFs(x, Q));
  };

  // Pre-compute all NNLO SIDIS operators (done once for all channels)
  const apfel::SidisNNLOObjects obj = apfel::InitializeSidisObjects(gx, gz, Thresholds, 1e-1);

  // Kinematics
  const double Q = 10;
  const double x = 0.1;
  const double zmin = 0.2, zmax = 0.8;

  std::cout << std::scientific << std::setprecision(8);
  std::cout << "Testing BuildSidis* at Q = " << Q << " GeV, x = " << x << "\n\n";

  // 1. FT unpolarised (NLO sum)
  const auto FT_full = apfel::BuildSidisUnpFT(obj, InPDFs, InFFs, as, Thresholds, 1, -1);
  const double ft_total = FT_full(Q).Evaluate1(x).Integrate(zmin, zmax);

  // Individual channels (LO ns=0, NLO ns=1, NLO qg=2, NLO gq=3)
  const auto FT0 = apfel::BuildSidisUnpFT(obj, InPDFs, InFFs, as, Thresholds, 1, 0);
  const auto FT2 = apfel::BuildSidisUnpFT(obj, InPDFs, InFFs, as, Thresholds, 1, 2);
  const auto FT3 = apfel::BuildSidisUnpFT(obj, InPDFs, InFFs, as, Thresholds, 1, 3);
  const auto FT1 = apfel::BuildSidisUnpFT(obj, InPDFs, InFFs, as, Thresholds, 1, 1);

  const double ft0 = FT0(Q).Evaluate1(x).Integrate(zmin, zmax);
  const double ft1 = FT1(Q).Evaluate1(x).Integrate(zmin, zmax);
  const double ft2 = FT2(Q).Evaluate1(x).Integrate(zmin, zmax);
  const double ft3 = FT3(Q).Evaluate1(x).Integrate(zmin, zmax);
  const double ft_sum = ft1 + ft2 + ft3;

  std::cout << "FT (NLO sum): " << ft_total << "\n";
  std::cout << "FT (channel sum): " << ft_sum << " (LO_ns: " << ft0 << ", NLO_ns: " << ft1-ft0 << ", NLO_qg: " << ft2 << ", NLO_gq: " << ft3 << ")\n";
  std::cout << "FT Ratio: " << ft_total / ft_sum << "\n\n";

  // 2. G1 polarised (LO sum)
  const auto G1_full = apfel::BuildSidisG1(obj, InpPDFs, InFFs, as, Thresholds, 0, -1);
  const double g1_total = G1_full(Q).Evaluate1(x).Integrate(zmin, zmax);

  const auto G1_ch0 = apfel::BuildSidisG1(obj, InpPDFs, InFFs, as, Thresholds, 0, 0);
  const double g1_ch0 = G1_ch0(Q).Evaluate1(x).Integrate(zmin, zmax);

  std::cout << "G1 (LO sum): " << g1_total << "\n";
  std::cout << "G1 (LO ns): " << g1_ch0 << "\n";
  std::cout << "G1 Ratio: " << g1_total / g1_ch0 << "\n\n";

  const double tol = 1e-8;
  bool ok = (std::abs(ft_total / ft_sum - 1.0) < tol) && (std::abs(g1_total / g1_ch0 - 1.0) < tol);

  if (ok)
    {
      std::cout << "ALL PASSED\n";
      return 0;
    }
  else
    {
      std::cout << "SOME FAILED\n";
      return 1;
    }
}
