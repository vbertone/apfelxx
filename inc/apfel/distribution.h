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
   * @brief The Distribution class defines one of the basic objects of
   * APFEL++. This is essentially the discretisation of a function
   * that can be conveniently used for convolutions.
   */
  class Distribution: public LagrangeInterpolator
  {
  public:
    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    /**
     * @brief The Distribution constructors.
     * @param gr: the Grid object that defines the x-space interpolation grid
     */
    Distribution(Grid const& gr);

    /**
     * @brief The Distribution constructors.
     * @param obj: a reference distribution from wich the grid is extracted
     * @param distsubgrid: the vector of the distribution on the subgrids
     * @param distjointgrid: the vector of the distribution on the joint grid
     */
    Distribution(Distribution           const& obj,
		 vector<vector<double>> const& distsubgrid,
		 vector<double>         const& distjointgrid);

    /**
     * @brief The Distribution constructors.
     * @param gr: the Grid object that defines the x-space interpolation grid
     * @param distsubgrid: the vector of the distribution on the subgrids
     * @param distjointgrid: the vector of the distribution on the joint grid
     */
    Distribution(Grid                   const& g,
		 vector<vector<double>> const& distsubgrid,
		 vector<double>         const& distjointgrid);

    /**
     * @brief The Distribution constructors.
     * @param gr: the Grid object that defines the x-space interpolation grid
     * @param InDistFunc: a function of x to be tabulated on the grid in x
     */
    Distribution(Grid                            const& g,
		 function<double(double const&)> const& InDistFunc);

    /**
     * @brief The Distribution constructors.
     * @param gr: the Grid object that defines the x-space interpolation grid
     * @param InDistFunc: a function of x and Q to be tabulated on the grid in x
     * @param Q: the value of Q in which InDistFunc has to be tabulated
     */
    Distribution(Grid                                           const& g,
		 function<double(double const&, double const&)> const& InDistFunc,
		 double                                         const& Q);

    /**
     * @brief The Distribution constructors.
     * @param gr: the Grid object that defines the x-space interpolation grid
     * @param InDistFunc: a function of ipdf and x to be tabulated on the grid in x
     * @param ipdf: the value of ipdf in which InDistFunc has to be tabulated
     */
    Distribution(Grid                                        const& g,
		 function<double(int const&, double const&)> const& InDistFunc,
		 int                                         const& ipdf);

    /**
     * @brief The Distribution constructors.
     * @param gr: the Grid object that defines the x-space interpolation grid
     * @param InDistFunc: a function of ipdf, x, and Q to be tabulated on the grid in x
     * @param ipdf: the value of ipdf in which InDistFunc has to be tabulated
     * @param Q: the value of Q in which InDistFunc has to be tabulated
     */
    Distribution(Grid                                                       const& g,
		 function<double(int const&, double const&, double const&)> const& InDistFunc,
		 int                                                        const& ipdf,
		 double                                                     const& Q);
    ///@}

    /**
     * @brief Function to push back the values of the joint grid.
     * @param xi: value of of the distribution to be appended to the distribution vector on the joint grid
     */
    void PushJointGrid(double const& xi);

    /**
     * @brief Function to push back the values of the subgrid.
     * @param xi: value of of the distribution to be appended to the distribution vector on one of the subgrids
     * @param next: switch to tell the function whether it should start filling in the next subgrid
     */
    void PushSubGrid(double const& xi, bool const& next);

    /**
     * @name Binary operators
     */
    ///@{
    Distribution& operator  = (Distribution const& d);                    //!< this  = Distribution
    Distribution& operator *= (double const& s);                          //!< this *= Scalar
    Distribution& operator *= (function<double(double const&)> const& f); //!< this *= Function of the integration variable
    Distribution& operator /= (double const& s);                          //!< this /= Scalar
    Distribution& operator *= (Distribution const& d);                    //!< this *= Distribution
    Distribution& operator += (Distribution const& d);                    //!< this += Distribution
    Distribution& operator -= (Distribution const& d);                    //!< this -= Distribution
    ///@}
  };

  /**
   * @name Ternary operators
   */
  ///@{
  Distribution operator * (double const& s, Distribution rhs);                          //!< Scalar*Distribution
  Distribution operator * (Distribution lhs, double const& s);                          //!< Distribution*Scalar
  Distribution operator * (function<double(double const&)> const& f, Distribution rhs); //!< Function*Distribution
  Distribution operator * (Distribution lhs, function<double(double const&)> const& f); //!< Distribution*Function
  Distribution operator / (Distribution lhs, double const& s);                          //!< Distribution/Scalar
  Distribution operator + (Distribution lhs, Distribution const& rhs);                  //!< Distribution+Distribution
  Distribution operator - (Distribution lhs, Distribution const& rhs);                  //!< Distribution-Distribution
  Distribution operator * (Distribution lhs, Distribution const& rhs);                  //!< Distribution*Distribution
  ///@}

  /**
   * @name Map of Distribution functions
   * Function that return maps pf distributions.
   */
  ///@{
  /**
   * @brief Function that fills in a map of distributions from a
   * map-valued function.
   * @param g: Grid object
   * @param InDistFunc: map-valued function dependent on x and a scale Q.
   * @param Q: the value of Q in which InDistFunc has to be tabulated
   * @param skip: vector of map indices to be skipped in the tabulation.
   */
  map<int,Distribution> DistributionMap(Grid                                                    const& g,
					function<map<int,double>(double const&, double const&)> const& InDistFunc,
					double                                                  const& Q,
					vector<int>                                             const& skip = {});

  /**
   * @brief Function that fills in a map of distributions from a
   * map-valued function.
   * @param g: Grid object
   * @param InDistFunc: map-valued function dependent on x
   * @param skip: vector of map indices to be skipped in the tabulation.
   */
  map<int,Distribution> DistributionMap(Grid                                     const& g,
					function<map<int,double>(double const&)> const& InDistFunc,
					vector<int>                              const& skip = {});
  ///@}
}
