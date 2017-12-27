//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/expression.h"

namespace apfel
{
  /**
   * @notes Expressions taken from hep-ph/0504192.
   */

  /**
   * @brief O(alpha_s) matching conditions.
   */
  class ATS1Hg_0: public Expression
  {
  public:
    ATS1Hg_0();
    double Regular(double const& x)  const;
  };

  class ATS1Hg_L: public Expression
  {
  public:
    ATS1Hg_L();
    double Regular(double const& x)  const;
  };

  //_________________________________________________________________________________
  class ATS1ggH_L: public Expression
  {
  public:
    ATS1ggH_L();
    double Local(double const& x)  const;
  };
}
