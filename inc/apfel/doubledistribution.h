//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/distribution.h"
#include "apfel/matrix.h"
#include "apfel/lagrangeinterpolator.h"
#include "apfel/doubleobject.h"

namespace apfel
{
  /**
   * @brief The DoubleDistribution class defines one of the basic
   * objects that provides a discretisation of a two-dimensional
   * function and that can be used for convolutions.
   */
  class DoubleDistribution
  {
  public:
    DoubleDistribution() = delete;

    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    /**
     * @brief The DoubleDistribution constructor
     * @param g1: the Grid object that defines the first interpolation grid
     * @param g2: the Grid object that defines the second interpolation grid
     * @param InDistFunc: the function of "x1" and "x2" to be tabulated
     */
    DoubleDistribution(Grid                                                const& g1,
                       Grid                                                const& g2,
                       std::function<double(double const&, double const&)> const& InDistFunc = [] (double const&, double const&) -> double{ return 0; });

    /**
     * @brief The DoubleDistribution copy constructor
     * @param obj: the object to be copied
     */
    DoubleDistribution(DoubleDistribution const& obj);

    /**
     * @brief The DoubleDistribution constructor
     * @param d1: the first single distribution
     * @param d2: the second single distribution
     */
    DoubleDistribution(Distribution const& d1, Distribution const& d2);

    /**
     * @brief The DoubleDistribution constructor
     * @param g1: the Grid object that defines the first interpolation grid
     * @param g2: the Grid object that defines the second interpolation grid
     * @param distsubgrid: the distribution on the subgrids
     * @param distjointgrid: the distribution on the joint grid
     */
    DoubleDistribution(Grid                                     const& g1,
                       Grid                                     const& g2,
                       std::vector<std::vector<matrix<double>>> const& distsubgrid,
                       matrix<double>                           const& distjointgrid);

    /**
     * @brief The DoubleDistribution constructor
     * @param DObj: the DoubleObject of distributions
     */
    DoubleDistribution(DoubleObject<Distribution> const& DObj);
    ///@}

    /**
     * @name Evaluate functions
     * List of functions that perform the interpolation on the x-space
     * grids.
     */
    ///@{
    /**
     * @brief Function that evaluates the interpolated function on the joint grids.
     * @param x1: the value in x1 to be interpolated
     * @param x2: the value in x2 to be interpolated
     * @return the interpolated result
     */
    double Evaluate(double const& x1, double const& x2) const;

    /**
     * @brief Function that evaluates the derivative of the
     * interpolated function on the joint grids.
     * @param x1: the value in x1 where the derivative has to be computed
     * @param x2: the value in x2 where the derivative has to be computed
     * @return the derivative of the interpolated function
     */
    double Derive(double const& x1, double const& x2) const;

    /**
     * @brief Function that evaluates the integral of the interpolated
     * function in the interval x1 in [a1,b1] and x2 in [a2,b2] on the
     * joint grids.
     * @param a1: the lower integration bound for the first variable
     * @param b1: the upper integration bound for the first variable
     * @param a2: the lower integration bound for the second variable
     * @param b2: the upper integration bound for the second variable
     * @return the integral of the interpolated function
     */
    double Integrate(double const& a1, double const& b1, double const& a2, double const& b2) const;

    /**
     * @brief Function that evaluates the DoubleDistribution in the
     * first variable leaving the second undetermined. This produces a
     * single distribution.
     * @param x1: the value of the first variable
     */
    Distribution Evaluate1(double const& x1) const;

    /**
     * @brief Function that evaluates the DoubleDistribution in the
     * second variable leaving the first undetermined. This produces a
     * single distribution.
     * @param x2: the value of the second variable
     */
    Distribution Evaluate2(double const& x2) const;

    /**
     * @brief Function that evaluates the derivative of a
     * DoubleDistribution in the first variable leaving the second
     * undetermined. This produces a single distribution.
     * @param x1: the value of the first variable
     */
    Distribution Derive1(double const& x1) const;

    /**
     * @brief Function that evaluates the derivative of a
     * DoubleDistribution in the second variable leaving the first
     * undetermined. This produces a single distribution.
     * @param x2: the value of the second variable
     */
    Distribution Derive2(double const& x2) const;

    /**
     * @brief Function that evaluates the integral of the interpolated
     * function in the interval [a1,b1] for the first variable leaving
     * the second undetermined.
     * @param a1: the lower integration bound for the first variable
     * @param b1: the upper integration bound for the first variable
     * @return the integral of the interpolated function
     */
    Distribution Integrate1(double const& a1, double const& b1) const;

    /**
     * @brief Function that evaluates the integral of the interpolated
     * function in the interval [a2,b2] for the second variable leaving
     * the first undetermined.
     * @param a2: the lower integration bound for the second variable
     * @param b2: the upper integration bound for the second variable
     * @return the integral of the interpolated function
     */
    Distribution Integrate2(double const& a2, double const& b2) const;
    ///@}

    /**
     * @brief Function that returns the derivative of the
     * DoubleDistribution in the form of a DoubleDistribution object.
     */
    DoubleDistribution Derivative() const;

    /**
     * @name Getters
     */
    ///@{
    Grid                                     const& GetFirstGrid()             const { return _g1; }
    Grid                                     const& GetSecondGrid()            const { return _g2; }
    LagrangeInterpolator                     const  GetFirstInterpolator()     const { return _li1; }
    LagrangeInterpolator                     const  GetSecondInterpolator()    const { return _li2; }
    std::vector<std::vector<matrix<double>>> const  GetDistributionSubGrid()   const { return _dDSubGrid; }
    matrix<double>                           const  GetDistributionJointGrid() const { return _dDJointGrid; }
    ///@}

    /**
     * @name Binary operators
     */
    ///@{
    DoubleDistribution& operator  = (DoubleDistribution const& d);                                  //!< this  = Distribution
    DoubleDistribution& operator *= (double const& s);                                              //!< this *= Scalar
    DoubleDistribution& operator *= (std::function<double(double const&, double const&)> const& f); //!< this *= Function of both the integration variables
    DoubleDistribution& operator *= (std::function<double(double const&)> const& f);                //!< this *= Function of one variable used for both the integration variables
    DoubleDistribution& operator /= (double const& s);                                              //!< this /= Scalar
    DoubleDistribution& operator *= (DoubleDistribution const& d);                                  //!< this *= Distribution
    DoubleDistribution& operator += (DoubleDistribution const& d);                                  //!< this += Distribution
    DoubleDistribution& operator -= (DoubleDistribution const& d);                                  //!< this -= Distribution
    ///@}

    /**
     * @brief Print the DoubleDistribution object
     */
    void Print() const { std::cout << *this << std::endl; }

  private:
    Grid                              const& _g1;           //!< The first interpolation grid
    Grid                              const& _g2;           //!< The second interpolation grid
    LagrangeInterpolator              const  _li1;          //!< The first Lagrange interpolator
    LagrangeInterpolator              const  _li2;          //!< The second Lagrange interpolator
    std::vector<std::vector<matrix<double>>> _dDSubGrid;    //!< The array with the double distribution values on the subgrids
    matrix<double>                           _dDJointGrid;  //!< The array with the double distribution values on the joint grids

    friend std::ostream& operator << (std::ostream& os, DoubleDistribution const& sg);
  };

  /**
   * @name Ternary operators
   */
  ///@{
  DoubleDistribution operator * (double const& s, DoubleDistribution rhs);                                              //!< Scalar*Distribution
  DoubleDistribution operator * (DoubleDistribution lhs, double const& s);                                              //!< Distribution*Scalar
  DoubleDistribution operator * (std::function<double(double const&, double const&)> const& f, DoubleDistribution rhs); //!< Function*Distribution
  DoubleDistribution operator * (DoubleDistribution lhs, std::function<double(double const&, double const&)> const& f); //!< Distribution*Function
  DoubleDistribution operator * (std::function<double(double const&)> const& f, DoubleDistribution rhs);                //!< Function*Distribution
  DoubleDistribution operator * (DoubleDistribution lhs, std::function<double(double const&)> const& f);                //!< Distribution*Function
  DoubleDistribution operator / (DoubleDistribution lhs, double const& s);                                              //!< Distribution/Scalar
  DoubleDistribution operator + (DoubleDistribution lhs, DoubleDistribution const& rhs);                                //!< Distribution+Distribution
  DoubleDistribution operator - (DoubleDistribution lhs, DoubleDistribution const& rhs);                                //!< Distribution-Distribution
  DoubleDistribution operator * (DoubleDistribution lhs, DoubleDistribution const& rhs);                                //!< Distribution*Distribution
  ///@}

  /**
   * @brief Method which prints DoubleDistribution with cout <<.
   */
  std::ostream& operator << (std::ostream& os, DoubleDistribution const& in);
}
