//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <apfel/alphaqcd.h>

#include <cmath>
#include <iostream>

using namespace apfel;
using namespace std;

int main()
{
  // Reference value of the strong coupling and heavy-quark
  // thresholds.
  const double         AlphaRef = 0.35;
  const double         MuRef    = sqrt(2);
  const vector<double> Masses   = {0, 0, 0, sqrt(2), 4.5, 175};

  // Iniatialize the running of the coupling at all available
  // perturbative orders.
  AlphaQCD asLO{AlphaRef, MuRef, Masses, 0};
  AlphaQCD asNLO{AlphaRef, MuRef, Masses, 1};
  AlphaQCD asNNLO{AlphaRef, MuRef, Masses, 2};
  AlphaQCD asNNNLO{AlphaRef, MuRef, Masses, 3};

  // Compute and print values at Mu.
  const auto Mu = 100.;
  cout << "\nLO:    alpha_s(Mu = " << Mu << " GeV) = " << asLO.Evaluate(Mu) << endl;
  cout << "NLO:   alpha_s(Mu = " << Mu << " GeV) = " << asNLO.Evaluate(Mu)  << " (NLO/LO     = " << 100 * asNLO.Evaluate(Mu) / asLO.Evaluate(Mu)<< "%)" << endl;
  cout << "NNLO:  alpha_s(Mu = " << Mu << " GeV) = " << asNNLO.Evaluate(Mu) << " (NNLO/NLO   = " << 100 * asNNLO.Evaluate(Mu) / asNLO.Evaluate(Mu)<< "%)" << endl;
  cout << "NNNLO: alpha_s(Mu = " << Mu << " GeV) = " << asNNNLO.Evaluate(Mu) << " (NNNLO/NNLO = " << 100 * asNNNLO.Evaluate(Mu) / asNNLO.Evaluate(Mu)<< "%)\n" << endl;

  return 0;
}
