//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/qgrid.h"
#include "apfel/matchedevolution.h"

#include <map>

namespace apfel
{
  /**
   * @brief The template TabulateObject class is a derived of the
   * QGrid class that tabulates on object of type T (it can be a
   * double, a Distribution, an Operator, Set<Distribution>, a
   * Set<Operator>) over a grid in Q, taking into account the possible
   * presence of thresholds, and provides the method to evaluate the
   * tabulated object at any generic value of Q.
   */
  template<class T>
  class TabulateObject: public QGrid<T>
  {
  public:
    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    /**
     * @brief The TabulateObject default constructor for an "evolving"
     * object (MatchedEvolution).
     * @param Object: the MatchedEvolution type object to be tabulated in Q
     * @param nQ: the number of on nodes of the grid in Q
     * @param QMin: the lower bound of the grid in Q
     * @param QMax: the upper bound of the grid in Q
     * @param InterDegree: the interpolation degree on the grid in Q
     * @param Lambda: the value of the parameter in the function ln(ln(Q<SUP>2</SUP>/&Lambda;<SUP>2</SUP>)) used for the tabulation (default: 0.25)
     */
    TabulateObject(MatchedEvolution<T> &Object,
		   int                 const& nQ,
		   double              const& QMin,
		   double              const& QMax,
		   int                 const& InterDegree,
		   double              const& Lambda = 0.25);

    /**
     * @brief The TabulateObject default constructor for a Q dependent 
     * object.
     * @param Object: the T-valued function to be tabulated in Q
     * @param nQ: the number of on nodes of the grid in Q
     * @param QMin: the lower bound of the grid in Q
     * @param QMax: the upper bound of the grid in Q
     * @param InterDegree: the interpolation degree on the grid in Q
     * @param Lambda: the value of the parameter in the function ln(ln(Q<SUP>2</SUP>/&Lambda;<SUP>2</SUP>)) used for the tabulation (default: 0.25)
     */
    TabulateObject(std::function<T(double)> const& Object,
		   int                      const& nQ,
		   double                   const& QMin,
		   double                   const& QMax,
		   int                      const& InterDegree,
		   std::vector<double>      const& Thresholds,
		   double                   const& Lambda = 0.25);

    /**
     * @brief The TabulateObject default constructor for a Q dependent 
     * object tabulated on a used given distribution in Q.
     * @param Object: the T-valued function to be tabulated in Q
     * @param nQ: the number of on nodes of the grid in Q
     * @param QMin: the lower bound of the grid in Q
     * @param QMax: the upper bound of the grid in Q
     * @param InterDegree: the interpolation degree on the grid in Q
     * @param TabFunc: the function to be used for the tabulation in Q
     * @param InvTabFunc: the inverse function of TabFunc (it has to be provided analytically)
     */
    TabulateObject(std::function<T(double)>             const& Object,
		   int                                  const& nQ,
		   double                               const& QMin,
		   double                               const& QMax,
		   int                                  const& InterDegree,
		   std::vector<double>                  const& Thresholds,
		   std::function<double(double const&)> const& TabFunc,
		   std::function<double(double const&)> const& InvTabFunc);
    ///@}

    /**
     * @name Evaluate functions
     * List of functions to interpolate the tabulated object.
     */
    ///@{
    /**
     * @brief This function interpolates in x and Q. It is used for T
     * = Distribution.
     */
    double EvaluatexQ(double const& x, double const& Q) const;
    /**
     * @brief This function interpolates in x and Q for a given
     * element i. It is used for T = Set<Distribution>.
     */
    double EvaluatexQ(int const& i, double const& x, double const& Q) const;
    /**
     * @brief This function interpolates in x, z, and Q. It is used
     * for T = DoubleOject<Distribution>.
     */
    double EvaluatexzQ(double const& x, double const& z, double const& Q) const;
    /**
     * @brief This function interpolates in x and Q and returns a
     * map. It is used for T = Set<Distribution>.
     */
    std::map<int,double> EvaluateMapxQ(double const& x, double const& Q) const;
    ///@}
  };
}
