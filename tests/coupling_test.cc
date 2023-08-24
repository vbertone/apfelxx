//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

#include <cmath>
#include <iostream>

int main()
{
  // Test scale
  double Mu = 1;

  // Reference value of the strong coupling and heavy-quark
  // thresholds.
  const double AlphaQCDRef = 0.118;
  const double MuQCDRef    = 91.1876;
  const std::vector<double> QuarkThresholds = {0, 0, 0, sqrt(2), 4.5, 175};

  // Iniatialize the running of the coupling at all available
  // perturbative orders.
  apfel::AlphaQCD asLO{AlphaQCDRef, MuQCDRef, QuarkThresholds, 0};
  apfel::AlphaQCD asNLO{AlphaQCDRef, MuQCDRef, QuarkThresholds, 1};
  apfel::AlphaQCD asNNLO{AlphaQCDRef, MuQCDRef, QuarkThresholds, 2};
  apfel::AlphaQCD asNNNLO{AlphaQCDRef, MuQCDRef, QuarkThresholds, 3};
  apfel::AlphaQCD asNNNNLO{AlphaQCDRef, MuQCDRef, QuarkThresholds, 4};

  // Compute and print values at Mu.
  std::cout << "\nNumerical evolution of the strong coupling:" << std::endl;
  std::cout << "LO:     alpha_s(Mu = " << Mu << " GeV) = " << asLO.Evaluate(Mu) << std::endl;
  std::cout << "NLO:    alpha_s(Mu = " << Mu << " GeV) = " << asNLO.Evaluate(Mu) << " (NLO/LO       = " << 100 * asNLO.Evaluate(Mu) / asLO.Evaluate(Mu)<< "%)" << std::endl;
  std::cout << "NNLO:   alpha_s(Mu = " << Mu << " GeV) = " << asNNLO.Evaluate(Mu) << " (NNLO/NLO     = " << 100 * asNNLO.Evaluate(Mu) / asNLO.Evaluate(Mu)<< "%)" << std::endl;
  std::cout << "NNNLO:  alpha_s(Mu = " << Mu << " GeV) = " << asNNNLO.Evaluate(Mu) << " (NNNLO/NNLO   = " << 100 * asNNNLO.Evaluate(Mu) / asNNLO.Evaluate(Mu)<< "%)" << std::endl;
  std::cout << "NNNNLO: alpha_s(Mu = " << Mu << " GeV) = " << asNNNNLO.Evaluate(Mu) << " (NNNNLO/NNNLO = " << 100 * asNNNNLO.Evaluate(Mu) / asNNNLO.Evaluate(Mu)<< "%)" << std::endl;

  // Iniatialize the running of the coupling at all available
  // perturbative orders.
  apfel::AlphaQCDg asLOg{AlphaQCDRef, MuQCDRef, QuarkThresholds, 0};
  apfel::AlphaQCDg asNLOg{AlphaQCDRef, MuQCDRef, QuarkThresholds, 1};
  apfel::AlphaQCDg asNNLOg{AlphaQCDRef, MuQCDRef, QuarkThresholds, 2};
  apfel::AlphaQCDg asNNNLOg{AlphaQCDRef, MuQCDRef, QuarkThresholds, 3};

  // Compute and print values at Mu.
  std::cout << "\nAnalytic evolution of the strong coupling:" << std::endl;
  std::cout << "LO:    alpha_s(Mu = " << Mu << " GeV) = " << asLOg.Evaluate(Mu) << std::endl;
  std::cout << "NLO:   alpha_s(Mu = " << Mu << " GeV) = " << asNLOg.Evaluate(Mu) << " (NLO/LO     = " << 100 * asNLOg.Evaluate(Mu) / asLOg.Evaluate(Mu)<< "%)" << std::endl;
  std::cout << "NNLO:  alpha_s(Mu = " << Mu << " GeV) = " << asNNLOg.Evaluate(Mu) << " (NNLO/NLO   = " << 100 * asNNLOg.Evaluate(Mu) / asNLOg.Evaluate(Mu)<< "%)" << std::endl;
  std::cout << "NNNLO: alpha_s(Mu = " << Mu << " GeV) = " << asNNNLOg.Evaluate(Mu) << " (NNNLO/NNLO = " << 100 * asNNNLOg.Evaluate(Mu) / asNNLOg.Evaluate(Mu)<< "%)" << std::endl;

  // Iniatialize the running of the coupling at all available
  // perturbative orders.
  const double xi = 1;
  apfel::AlphaQCDxi asLOxi{AlphaQCDRef, MuQCDRef, QuarkThresholds, 0, xi};
  apfel::AlphaQCDxi asNLOxi{AlphaQCDRef, MuQCDRef, QuarkThresholds, 1, xi};
  apfel::AlphaQCDxi asNNLOxi{AlphaQCDRef, MuQCDRef, QuarkThresholds, 2, xi};
  apfel::AlphaQCDxi asNNNLOxi{AlphaQCDRef, MuQCDRef, QuarkThresholds, 3, xi};

  // Compute and print values at Mu.
  std::cout << "\nNumerical evolution of the strong coupling with resummation scale (xi = " << xi << "):" << std::endl;
  std::cout << "LO:    alpha_s(Mu = " << Mu << " GeV) = " << asLOxi.Evaluate(Mu) << std::endl;
  std::cout << "NLO:   alpha_s(Mu = " << Mu << " GeV) = " << asNLOxi.Evaluate(Mu) << " (NLO/LO     = " << 100 * asNLOxi.Evaluate(Mu) / asLOxi.Evaluate(Mu)<< "%)" << std::endl;
  std::cout << "NNLO:  alpha_s(Mu = " << Mu << " GeV) = " << asNNLOxi.Evaluate(Mu) << " (NNLO/NLO   = " << 100 * asNNLOxi.Evaluate(Mu) / asNLOxi.Evaluate(Mu)<< "%)" << std::endl;
  std::cout << "NNNLO: alpha_s(Mu = " << Mu << " GeV) = " << asNNNLOxi.Evaluate(Mu) << " (NNNLO/NNLO = " << 100 * asNNNLOxi.Evaluate(Mu) / asNNLOxi.Evaluate(Mu)<< "%)" << std::endl;

  // Ratio between numerical and analytic solutions order by order
  std::cout << "\nRatio between numerical and analytic solutions:" << std::endl;
  std::cout << "LO:    " << 100 * asLO.Evaluate(Mu) / asLOg.Evaluate(Mu) << "%"<< std::endl;
  std::cout << "NLO:   " << 100 * asNLO.Evaluate(Mu) / asNLOg.Evaluate(Mu)<< "%" << std::endl;
  std::cout << "NNLO:  " << 100 * asNNLO.Evaluate(Mu) / asNNLOg.Evaluate(Mu)<< "%" << std::endl;
  std::cout << "NNNLO: " << 100 * asNNNLO.Evaluate(Mu) / asNNNLOg.Evaluate(Mu)<< "%" << std::endl;

  // Reference value of the QED coupling and heavy-quark
  // thresholds.
  const double AlphaQEDRef = 1. / 128.;
  const double MuQEDRef    = MuQCDRef;
  const std::vector<double> LeptThresholds = {0, 0, 1.777};

  // Iniatialize the running of the QED coupling at all available
  // perturbative orders.
  apfel::AlphaQED aLO{AlphaQEDRef, MuQEDRef, QuarkThresholds, LeptThresholds, 0};
  apfel::AlphaQED aNLO{AlphaQEDRef, MuQEDRef, QuarkThresholds, LeptThresholds, 1};

  // Compute and print values at Mu.
  Mu = 1e10;
  std::cout << "\nNumeric evolution of the electromagnetic coupling:" << std::endl;
  std::cout << "LO:    alpha_em(Mu = " << Mu << " GeV) = " << aLO.Evaluate(Mu) << std::endl;
  std::cout << "NLO:   alpha_em(Mu = " << Mu << " GeV) = " << aNLO.Evaluate(Mu)  << " (NLO/LO = " << 100 * aNLO.Evaluate(Mu) / aLO.Evaluate(Mu)<< "%)" << std::endl;

  // Compute mixed evolution
  apfel::AlphaQCDQED aLOmix{AlphaQCDRef, AlphaQEDRef, MuQEDRef, QuarkThresholds, LeptThresholds, 0};
  apfel::AlphaQCDQED aNLOmix{AlphaQCDRef, AlphaQEDRef, MuQEDRef, QuarkThresholds, LeptThresholds, 1};

  // Compute and print values at Mu.
  Mu = 1;
  std::cout << "\nCoupled numerical evolution of strong and electromagnetic couplings:" << std::endl;
  std::cout << "LO:    alpha_s(Mu = " << Mu << " GeV)[coup]  = " << aLOmix.Evaluate(Mu)(0, 0)
            << ", alpha_s(Mu = " << Mu << " GeV)[dec]  = " << asLO.Evaluate(Mu)
            << ", (ratio = " << aLOmix.Evaluate(Mu)(0, 0) / asLO.Evaluate(Mu) << ")" <<std::endl;
  std::cout << "LO:    alpha_em(Mu = " << Mu << " GeV)[coup] = " << aLOmix.Evaluate(Mu)(1, 0)
            << ", alpha_em(Mu = " << Mu << " GeV)[dec] = " << aLO.Evaluate(Mu)
            << ", (ratio = " << aLOmix.Evaluate(Mu)(1, 0) / aLO.Evaluate(Mu) << ")" <<std::endl;
  std::cout << "NLO:   alpha_s(Mu = " << Mu << " GeV)[coup]  = " << aNLOmix.Evaluate(Mu)(0, 0)
            << ", alpha_s(Mu = " << Mu << " GeV)[dec]  = " << asNLO.Evaluate(Mu)
            << ", (ratio = " << aNLOmix.Evaluate(Mu)(0, 0) / asNLO.Evaluate(Mu) << ")" <<std::endl;
  std::cout << "NLO:   alpha_em(Mu = " << Mu << " GeV)[coup] = " << aNLOmix.Evaluate(Mu)(1, 0)
            << ", alpha_em(Mu = " << Mu << " GeV)[dec] = " << aNLO.Evaluate(Mu)
            << ", (ratio = " << aNLOmix.Evaluate(Mu)(1, 0) / aNLO.Evaluate(Mu) << ")\n" <<std::endl;

  return 0;
}
