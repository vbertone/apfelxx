//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/lagrangeinterpolator.h"

#include <functional>
#include <map>

using std::function;
using std::map;

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
     *
     * @param gr the Grid object
     */
    Distribution(Grid const& gr);

    /**
     * @brief Distribution constructors.
     *
     * @param obj Distribution object
     * @param distsubgrid a 2d vector with the distribution values for each subgrid.
     * @param distjointgrid a vector with the distribution values on the joint grid.
     */
    Distribution(Distribution           const& obj,
		 vector<vector<double>> const& distsubgrid,
		 vector<double>         const& distjointgrid);

    /**
     * @brief Distribution constructors.
     *
     * @param gr the Grid object
     * @param distsubgrid a 2d vector with the distribution values for each subgrid.
     * @param distjointgrid a vector with the distribution values on the joint grid.
     */
    Distribution(Grid                   const& g,
		 vector<vector<double>> const& distsubgrid,
		 vector<double>         const& distjointgrid);

    /**
     * @brief Distribution constructors.
     *
     * @param gr the Grid object
     * @param InDistFunc function of ipdf and x to be tabulated.
     * @param ipdf int to be fed to InDistFunc.
     */
    Distribution(Grid                            const& g,
		 function<double(double const&)> const& InDistFunc);

    /**
     * @brief Distribution constructors.
     *
     * @param gr the Grid object
     * @param InDistFunc function of ipdf, x, and Q to be tabulated.
     * @param Q double to be fed to InDistFunc.
     */
    Distribution(Grid                                           const& g,
		 function<double(double const&, double const&)> const& InDistFunc,
		 double                                         const& Q);

    /**
     * @brief Distribution constructors.
     *
     * @param gr the Grid object
     * @param InDistFunc function of ipdf and x to be tabulated.
     * @param ipdf int to be fed to InDistFunc.
     */
    Distribution(Grid                                        const& g,
		 function<double(int const&, double const&)> const& InDistFunc,
		 int                                         const& ipdf);

    /**
     * @brief Distribution constructors.
     *
     * @param gr the Grid object
     * @param InDistFunc function of ipdf, x, and Q to be tabulated.
     * @param ipdf int to be fed to InDistFunc.
     * @param Q double to be fed to InDistFunc.
     */
    Distribution(Grid                                                       const& g,
		 function<double(int const&, double const&, double const&)> const& InDistFunc,
		 int                                                        const& ipdf,
		 double                                                     const& Q);

    /**
     * @brief Function to push back the values of the joint grid.
     */
    void PushJointGrid(double const& xi);

    /**
     * @brief Function to push back the values of the subgrid.
     */
    void PushSubGrid(double const& xi, bool const& next);

    // Operators
    Distribution& operator  = (Distribution const& d);                    //!< this  = Distribution
    Distribution& operator *= (double const& s);                          //!< this *= Scalar
    Distribution& operator *= (function<double(double const&)> const& f); //!< this *= Function of the integration variable
    Distribution& operator /= (double const& s);                          //!< this /= Scalar
    Distribution& operator *= (Distribution const& d);                    //!< this *= Distribution
    Distribution& operator += (Distribution const& d);                    //!< this += Distribution
    Distribution& operator -= (Distribution const& d);                    //!< this -= Distribution
  };

  // Extra operation definitions where Distribution is at the left hand side (lhs).
  Distribution operator * (double const& s, Distribution rhs);                          //!< Scalar*Distribution
  Distribution operator * (Distribution lhs, double const& s);                          //!< Distribution*Scalar
  Distribution operator * (function<double(double const&)> const& f, Distribution rhs); //!< Function*Distribution
  Distribution operator * (Distribution lhs, function<double(double const&)> const& f); //!< Distribution*Function
  Distribution operator / (Distribution lhs, double const& s);                          //!< Distribution/Scalar
  Distribution operator + (Distribution lhs, Distribution const& rhs);                  //!< Distribution+Distribution
  Distribution operator - (Distribution lhs, Distribution const& rhs);                  //!< Distribution-Distribution
  Distribution operator * (Distribution lhs, Distribution const& rhs);                  //!< Distribution*Distribution

  // Fill in an undordered_map of distributions from a map of distributions.
  map<int,Distribution> DistributionMap(Grid                                                    const& g,
					function<map<int,double>(double const&, double const&)> const& InDistFunc,
					double                                                  const& Q,
					vector<int>                                             const& skip = {});

  // Fill in an undordered_map of distributions from a map of distributions.
  map<int,Distribution> DistributionMap(Grid                                     const& g,
					function<map<int,double>(double const&)> const& InDistFunc,
					vector<int>                              const& skip = {});
}
