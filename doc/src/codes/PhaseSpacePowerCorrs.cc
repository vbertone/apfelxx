/*
 * Author: Valerio Bertone
 */

#include <iostream>

// NangaParbat libs
#include <apfel/apfelxx.h>

int main() {
  // Define lepton cuts
  const double pTmin = 20;
  const double etamax = 2.4;
  const double etamin = - etamax;

  // Phase-space reduction factor
  apfel::TwoBodyPhaseSpace ps{pTmin, etamin, etamax, 1e-20};

  // Vector boson kinematics
  const double Q = 91;
  const double y = 0;

  // Analytic solution at qT = y = 0
  double psapprox = 0;
  if (Q >= 2 * pTmin && Q < 2 * pTmin * cosh(etamax))
    psapprox = ( 1 - pow(pTmin / Q, 2) ) * sqrt( 1 - pow(2 * pTmin / Q, 2) );
  else if (Q >= 2 * pTmin * cosh(etamax))
    psapprox = tanh(etamax) * ( 1 - 1 / pow(2 * cosh(etamax), 2 ) );

  const std::vector<double> qTv{0, 0.01, 0.1, 1, 5, 10, 15, 20, 30, 40, 50};

  for (auto const& qT : qTv)
    std::cout << std::scientific << qT << "\t" << ps.PhaseSpaceReduction(Q, y, qT) / psapprox << std::endl;
  std::cout << "\n";

  return 0;
}
