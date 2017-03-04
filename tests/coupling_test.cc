//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <cmath>
#include <iostream>

#include <apfel/alphaqcd.h>

using namespace apfel;
using namespace std;

int main()
{
  const auto           AlphaRef = 0.35;
  const auto           MuRef    = sqrt(2);
  const vector<double> Masses   = {0, 0, 0, sqrt(2), 4.5, 175};

  AlphaQCD asLO{AlphaRef, MuRef, Masses, 0};
  AlphaQCD asNLO{AlphaRef, MuRef, Masses, 1};
  AlphaQCD asNNLO{AlphaRef, MuRef, Masses, 2};

  const auto Mu = 100.;
  cout << "\nLO:    alpha_s(Mu = " << Mu << " GeV) = " << asLO.GetObject(Mu) << endl;
  cout << "NLO:   alpha_s(Mu = " << Mu << " GeV) = " << asNLO.GetObject(Mu)  << endl;
  cout << "NNLO:  alpha_s(Mu = " << Mu << " GeV) = " << asNNLO.GetObject(Mu) << "\n" << endl;

  return 0;
}
