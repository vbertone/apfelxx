// Validate `BuildSidisUnpFT` and `BuildSidisG1` against `sidisbuilder_test.dat`
// (channel -1) from the repository https://github.com/vbertone/SIDIS_NNLO_tests.
//
// Uses the new `InitializeSidisObjects` + Build* API from `sidisbuilder.h`.
// Results are printed to stdout and dumped to `sidisbuilder_test_output.dat`.

#include <LHAPDF/LHAPDF.h>
#include <apfel/apfelxx.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <string>

int main()
{
  apfel::SetVerbosityLevel(1);

  // Use a coarse grid for a quick benchmark.
  const apfel::Grid gx{{apfel::SubGrid{30, 1e-3, 3}, apfel::SubGrid{20, 1e-1, 3}, apfel::SubGrid{15, 8e-1, 3}}};
  const apfel::Grid gz{{apfel::SubGrid{30, 1e-2, 3}, apfel::SubGrid{20, 2e-1, 3}, apfel::SubGrid{15, 8e-1, 3}}};

  // PDF / FF sets
  LHAPDF::PDF* PDFs  = LHAPDF::mkPDF("MSHT20nnlo_as118");
  LHAPDF::PDF* pPDFs = LHAPDF::mkPDF("BDSSV24pol_NNLO");
  LHAPDF::PDF* FFs   = LHAPDF::mkPDF("BDSSV_NNLO_PionPlus");

  const std::vector<double> Thresholds{0, 0, 0, PDFs->quarkThreshold(4), PDFs->quarkThreshold(5)};

  const auto as      = [&] (double const& mu) { return PDFs->alphasQ(mu); };
  auto       InPDFs  = [&] (double const& x, double const& Q) { return PDFs->xfxQ(x, Q); };
  auto       InpPDFs = [&] (double const& x, double const& Q) { return pPDFs->xfxQ(x, Q); };
  auto       InFFs   = [&] (double const& x, double const& Q) { return FFs->xfxQ(x, Q); };

  // Pre-compute all NNLO SIDIS operators
  const apfel::SidisNNLOObjects obj = apfel::InitializeSidisObjects(gx, gz, Thresholds, 1e-2);

  // Build full NNLO cross-section functions
  const auto FT = apfel::BuildSidisUnpFT(obj, InPDFs,  InFFs, as, Thresholds, 2);
  const auto G1 = apfel::BuildSidisG1(obj, InpPDFs, InFFs, as, Thresholds, 2);

  // Kinematics matching benchmark.dat
  const std::vector<double> xs {0.0052, 0.0079, 0.0142, 0.0245, 0.0346, 0.0487,
                                 0.0765, 0.121,  0.172,  0.24,   0.341,  0.48};
  const std::vector<double> Q2s{1.25001, 1.46,  2.12,  3.22,  4.36,  5.97,
                                 8.96,   13.8,  19.6,  27.6,  40.1,  55.6};
  const double zmin = 0.2, zmax = 0.85;

  // Read APFEL++ reference values from benchmark.dat (channel -1 = full NNLO).
  // Columns: i Q2 x zmin zmax G1_APFEL G1_Zurich F1_APFEL F1_Zurich
  std::vector<double> g1ref, F1ref;
  {
    std::ifstream f("sidisbuilder_test.dat");
    if (!f) { std::cerr << "Cannot open sidisbuilder_test.dat\n"; return 1; }
    bool in_ch = false;
    std::string line;
    while (std::getline(f, line))
      {
        if (line.find("# Computing channel -1") != std::string::npos) { in_ch = true; continue; }
        if (in_ch && line.find("# Computing channel") != std::string::npos) break;
        if (!in_ch || line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        int idx; double Q2, x, z1, z2, g1a, g1z, F1a, F1z;
        ss >> idx >> Q2 >> x >> z1 >> z2 >> g1a >> g1z >> F1a >> F1z;
        g1ref.push_back(g1a);
        F1ref.push_back(F1a);
      }
  }
  if ((int) g1ref.size() != (int) xs.size())
    { std::cerr << "Reference data size mismatch\n"; return 1; }

  std::ofstream fout("sidisbuilder_test_output.dat");
  fout << std::scientific;
  std::cout << std::scientific;

  const std::string header =
    "#i\tQ2\t\t\tx\t\t\tzmin\t\tzmax\t\tG1(builder)\tG1(benchmark)\tFT(builder)\tFT(benchmark)";
  std::cout << header << "\n";
  fout  << header << "\n";

  bool all_ok = true;
  const double tol = 0.10; // 10% — coarse grid

  for (int i = 0; i < (int) xs.size(); i++)
    {
      const double Q   = std::sqrt(Q2s[i]);
      const double x   = xs[i];
      const double g1v = G1(Q).Evaluate1(x).Integrate(zmin, zmax);
      const double F1v = FT(Q).Evaluate1(x).Integrate(zmin, zmax);

      const double rg1 = (g1ref[i] != 0.) ? g1v / g1ref[i] - 1. : g1v;
      const double rF1 = (F1ref[i] != 0.) ? F1v / F1ref[i] - 1. : F1v;
      if (std::abs(rg1) >= tol || std::abs(rF1) >= tol) all_ok = false;

      const std::string row = " " + std::to_string(i)
        + "\t" + std::to_string(Q2s[i])
        + "\t" + std::to_string(x)
        + "\t" + std::to_string(zmin)
        + "\t" + std::to_string(zmax)
        + "\t" + std::to_string(g1v)
        + "\t" + std::to_string(g1ref[i])
        + "\t" + std::to_string(F1v)
        + "\t" + std::to_string(F1ref[i]);

      std::cout << row << "\n";
      fout  << row << "\n";
    }

  const std::string verdict = all_ok ? "# ALL PASSED" : "# SOME FAILED";
  std::cout << "\n" << verdict << "\n";
  fout  << "\n" << verdict << "\n";
  fout.close();

  return all_ok ? 0 : 1;
}
