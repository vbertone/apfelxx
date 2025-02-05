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
  // Reference value of the strong coupling and heavy-quark
  // thresholds.
  const double AlphaRef = 0.35;
  const double MuRef    = sqrt(2);
  const std::vector<double> Masses = {0, 0, 0, 1.25, 4.2, 160};

  // Iniatialize the running of the coupling at all available
  // perturbative orders.
  apfel::AlphaQCDMSbarMass AlphasLO{AlphaRef, MuRef, Masses, 0};
  apfel::AlphaQCDMSbarMass AlphasNLO{AlphaRef, MuRef, Masses, 1};
  apfel::AlphaQCDMSbarMass AlphasNNLO{AlphaRef, MuRef, Masses, 2};
  const auto asLO = [&] (double const& mu) -> double{ return AlphasLO.Evaluate(mu); };
  const auto asNLO = [&] (double const& mu) -> double{ return AlphasNLO.Evaluate(mu); };
  const auto asNNLO = [&] (double const& mu) -> double{ return AlphasNNLO.Evaluate(mu); };

  // Initialise running mass
  apfel::MSbarMass mcLO{Masses[3], Masses[3], Masses, 0, asLO};
  apfel::MSbarMass mcNLO{Masses[3], Masses[3], Masses, 1, asNLO};
  apfel::MSbarMass mcNNLO{Masses[3], Masses[3], Masses, 2, asNNLO};

  // Test running of the mass
  const int nmu = 20;
  const double mumin = 1;
  const double mumax = 100;
  const double mustp = exp(log(mumax / mumin) / ( nmu - 1 ));

  // Compute and print values of the coupling at mumax.
  std::cout << "\nLO:   alpha_s(Mu = " << mumax << " GeV) = " << asLO(mumax) << std::endl;
  std::cout << "NLO:  alpha_s(Mu = " << mumax << " GeV) = " << asNLO(mumax) << " (NLO/LO   = " << 100 * asNLO(mumax) / asLO(mumax)<< "%)" << std::endl;
  std::cout << "NNLO: alpha_s(Mu = " << mumax << " GeV) = " << asNNLO(mumax) << " (NNLO/NLO = " << 100 * asNNLO(mumax) / asNLO(mumax)<< "%)" << std::endl;

  // Now tabulate masses
  std::cout << std::scientific << "\nmu [GeV]\tmc(mu)@LO\tmc(mu)@NLO\tmc(mu)@NNLO\n";
  for (double mu = mumin; mu <= 1.00001 * mumax; mu *= mustp)
    std::cout << mu << "\t" << mcLO.Evaluate(mu) << "\t" << mcNLO.Evaluate(mu) << "\t" << mcNNLO.Evaluate(mu) << std::endl;
  std::cout << "\n";

  return 0;
}
