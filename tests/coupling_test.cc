//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/alphaqcd.h>
#include <apfel/alphaqed.h>

#include <cmath>
#include <iostream>

using namespace apfel;
using namespace std;

int main()
{
  // Test scale
  double Mu = 100.;

  // Reference value of the strong coupling and heavy-quark
  // thresholds.
  const double         AlphaQCDRef     = 0.35;
  const double         MuQCDRef        = sqrt(2);
  const vector<double> QuarkThresholds = {0, 0, 0, sqrt(2), 4.5, 175};

  // Iniatialize the running of the coupling at all available
  // perturbative orders.
  AlphaQCD asLO{AlphaQCDRef, MuQCDRef, QuarkThresholds, 0};
  AlphaQCD asNLO{AlphaQCDRef, MuQCDRef, QuarkThresholds, 1};
  AlphaQCD asNNLO{AlphaQCDRef, MuQCDRef, QuarkThresholds, 2};
  AlphaQCD asNNNLO{AlphaQCDRef, MuQCDRef, QuarkThresholds, 3};

  // Compute and print values at Mu.
  cout << "\nLO:    alpha_s(Mu = " << Mu << " GeV) = " << asLO.Evaluate(Mu) << endl;
  cout << "NLO:   alpha_s(Mu = " << Mu << " GeV) = " << asNLO.Evaluate(Mu)  << " (NLO/LO     = " << 100 * asNLO.Evaluate(Mu) / asLO.Evaluate(Mu)<< "%)" << endl;
  cout << "NNLO:  alpha_s(Mu = " << Mu << " GeV) = " << asNNLO.Evaluate(Mu) << " (NNLO/NLO   = " << 100 * asNNLO.Evaluate(Mu) / asNLO.Evaluate(Mu)<< "%)" << endl;
  cout << "NNNLO: alpha_s(Mu = " << Mu << " GeV) = " << asNNNLO.Evaluate(Mu) << " (NNNLO/NNLO = " << 100 * asNNNLO.Evaluate(Mu) / asNNLO.Evaluate(Mu)<< "%)" << endl;

  // Reference value of the QED coupling and heavy-quark
  // thresholds.
  const double         AlphaQEDRef    = 1. / 128.;
  const double         MuQEDRef       = 91.2;
  const vector<double> LeptThresholds = {0, 0, 1.777};

  // Iniatialize the running of the QED coupling at all available
  // perturbative orders.
  AlphaQED aLO{AlphaQEDRef, MuQEDRef, QuarkThresholds, LeptThresholds, 0};
  AlphaQED aNLO{AlphaQEDRef, MuQEDRef, QuarkThresholds, LeptThresholds, 1};

  // Compute and print values at Mu.
  Mu = 1e10;
  cout << "\nLO:    alpha_em(Mu = " << Mu << " GeV) = " << aLO.Evaluate(Mu) << endl;
  cout << "NLO:   alpha_em(Mu = " << Mu << " GeV) = " << aNLO.Evaluate(Mu)  << " (NLO/LO = " << 100 * aNLO.Evaluate(Mu) / aLO.Evaluate(Mu)<< "%)\n" << endl;

  return 0;
}
