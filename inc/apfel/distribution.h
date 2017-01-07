//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <apfel/lagrangeinterpolator.h>

namespace apfel
{
  /**
   * @brief The Distribution class for PDFs.
   *
   * This class provides methods to inherit a custom PDF distribution
   * in a generic basis.
   */
  class Distribution: public LagrangeInterpolator
  {
  public:
    /**
     * @brief Distribution constructors.
     * @param gr the Grid object
     */
    Distribution(Grid const& gr);

    /**
     * @brief Distribution constructors.
     * @param gr the Grid object
     * @param distsubgrid a 2d vector with the distribution values for each subgrid.
     * @param distjointgrid a vector with the distribution values on the joint grid.
     */
    Distribution(Distribution const& obj, vector<vector<double>> const& distsubgrid, vector<double> const& distjointgrid);

    // operators
    Distribution& operator*=(double const& s);       //!< this *= Scalar
    Distribution& operator+=(Distribution const& d); //!< this += Distribution
    Distribution& operator-=(Distribution const& d); //!< this -= Distribution
  };

  // Extra operation definitions where Distribution is at the left hand side (lhs).
  Distribution operator*(double const& s, Distribution rhs);         //!< Scalar*Distribution
  Distribution operator*(Distribution lhs, double const& s);         //!< Distribution*Scalar
  Distribution operator+(Distribution lhs, Distribution const& rhs); //!< Distribution+Distribution
  Distribution operator-(Distribution lhs, Distribution const& rhs); //!< Distribution-Distribution
}
