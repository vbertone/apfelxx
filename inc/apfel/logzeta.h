//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

namespace apfel
{
  // Coeffiecient of log(zeta) in the zeta-prescription scheme.
  double Lzetaq10();
  double Lzetaq11();
  double Lzetaq20(int const& nf);
  double Lzetaq22(int const& nf);
  double Lzetag10(int const& nf);
  double Lzetag11();
  double Lzetag20(int const& nf);
  double Lzetag22(int const& nf);
}
