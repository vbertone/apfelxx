//
// APFEL++ 2017
//
// Author: Emanuele R. Nocera: emanuele.roberto.nocera@gmail.com
//

#include <apfel/apfelxx.h>

#include <cmath>
#include <map>
#include <iomanip>

int main()
{
  // Test scales
  std::vector<double> Mu = {sqrt(2), 4, 50, 200};

  // Test perturbative order
  int pto = 2;

  //Reference values of the strong coupling and HQ thresholds
  const double AlphaQCDRef = 0.35;
  const double MuQCDRef    = sqrt(2);
  const std::vector<double> QuarkThresholds = {0, 0, 0, sqrt(2), 4.5, 175};
  const apfel::AlphaQCD asQCD{AlphaQCDRef, MuQCDRef, QuarkThresholds, pto};

  // Reference value of the QED coupling and HQ thresholds
  const double AlphaQEDRef = 1. / 128.;
  const double MuQEDRef    = 91.2;
  const std::vector<double> LeptThresholds = {0, 0, 1.777};
  const apfel::AlphaQED asQED{AlphaQEDRef, MuQEDRef, QuarkThresholds, LeptThresholds, 0};

  // Compute SIA total cross section
  std::cout << std::setw(10) << "Q2"
            << std::setw(18) << "SIA total Xsec"
            << std::endl;
  for (double mu : Mu)
    std::cout << std::setprecision(4)
              << std::setw(10) << mu
              << std::setprecision(4)
              << std::setw(10) << apfel::GetSIATotalCrossSection(pto, mu, asQCD.Evaluate(mu), asQED.Evaluate(mu), QuarkThresholds)
              << " [nbarn]" << std::endl;

  return 0;
}
