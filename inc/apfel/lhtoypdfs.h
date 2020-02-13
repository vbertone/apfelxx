//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include <map>

namespace apfel
{
  /// @cond UNNECESSARY
  /**
   * @name Les Houches PDFs
   * Les Houches toy parameterisation at Q = &radic;2 GeV. Used for
   * test purposes.
   */
  ///@{
  double xupv(double const& x);
  double xdnv(double const& x);
  double xglu(double const& x);
  double xdbar(double const& x);
  double xubar(double const& x);
  double xsbar(double const& x);
  std::map<int,double> LHToyPDFs(double const& x, double const&);
  ///@}
  /// @endcond
}
