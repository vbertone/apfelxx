//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include <map>

namespace apfel
{
  /// @cond UNNECESSARY
  /**
   * @name Les Houches PDFs
   * Les Houches toy parameterisation at Q = &radic;2 GeV. Used for
   * test purposes. This also includes a set of longitudinally
   * polarised PDFs and unpolarised FFs.
   */
  ///@{
  std::map<int, double> Provider(double const& x, int const& isel);
  std::map<int, double> LHToyPDFs(double const& x, double const&);
  std::map<int, double> LHToyPDFsPhys(double const& x, double const&);
  std::map<int, double> LHToyPDFsPol(double const& x, double const&);
  std::map<int, double> LHToyFFs(double const& x, double const&);
  ///@}
  /// @endcond
}
