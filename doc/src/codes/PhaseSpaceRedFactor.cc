/*
 * Author: Valerio Bertone
 */

#include <iostream>

// NangaParbat libs
#include <NangaParbat/twobodyphasespace.h>

int main() {
  // Define lepton cuts
  const double pTmin = 20;
  const double etamin = -2.4;
  const double etamax = 2.4;

  // Phase-space reduction factor
  NangaParbat::TwoBodyPhaseSpace ps{pTmin, etamin, etamax};

  // Vector boson kinematics
  const double Q  = 91;
  double qT;

  const int ny = 10000;
  const double ymin = etamin + 1e-2;
  const double ymax = etamax - 1e-2;
  const double step = ( ymax - ymin ) / ( ny - 1 );

  qT = 1;
  for (double y = ymin; y < ymax; y += step)
    std::cout << std::scientific << qT << "\t" << y << "\t" <<ps.ParityViolatingPhaseSpaceReduction(Q, y, qT) / ps.PhaseSpaceReduction(Q, y, qT) << std::endl;
  std::cout << "\n\n";

  qT = 3;
  for (double y = ymin; y < ymax; y += step)
    std::cout << std::scientific << qT << "\t" << y << "\t" <<ps.ParityViolatingPhaseSpaceReduction(Q, y, qT) / ps.PhaseSpaceReduction(Q, y, qT) << std::endl;
  std::cout << "\n\n";

  qT = 10;
  for (double y = ymin; y < ymax; y += step)
    std::cout << std::scientific << qT << "\t" << y << "\t" <<ps.ParityViolatingPhaseSpaceReduction(Q, y, qT) / ps.PhaseSpaceReduction(Q, y, qT) << std::endl;
  std::cout << "\n\n";
  return 0;
}
