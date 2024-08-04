//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/matchedevolution.h"
#include "apfel/set.h"
#include "apfel/operator.h"

namespace apfel
{
  /**
   * @brief The Dglap class is specialization class of the
   * MatchedEvolution class for the computation of the DGLAP
   * evolution.
   */
  template<class T>
  class Dglap: public MatchedEvolution<Set<T>>
  {
  public:
    Dglap() =  delete;

    /**
     * @brief Dglap constructor.
     * @param SplittingFunctions: set of splitting functions
     * @param MatchingConditions: set of matching conditions
     * @param InhomogeneousTerms: inhomogeneous term to be added to the derivative
     * @param ObjRef: reference object to be evolved
     * @param MuRef: reference scale from which the evolution starts
     * @param Thresholds: vector of the heavy-quark thresholds
     * @param nsteps: number of steps of the ODE solver (default: 10)
     */
    Dglap(std::function<Set<Operator>(int const&, double const&)> const& SplittingFunctions,
          std::function<Set<Operator>(bool const&, int const&)>   const& MatchingConditions,
          std::function<Set<T>(int const&, double const&)>        const& InhomogeneousTerms,
          Set<T>                                                  const& ObjRef,
          double                                                  const& MuRef,
          std::vector<double>                                     const& Thresholds,
          int                                                     const& nsteps = 10);

    /**
     * @brief Function that matches the evolved object at the thresholds.
     * @param Up: whether the matching has to be done upward or backward
     * @param nf: number of active flavours on this side of the threshold
     * @param sd: object on this side of the threshold
     * @return The matched object on the other side of the threshold.
     */
    Set<T> MatchObject(bool const& Up, int const& nf, Set<T> const& sd) const;

    /**
     * @brief This function returns the r.h.s. of the DGLAP equation,
     * i.e. the convolution between splitting functions and
     * distributions.
     * @param nf: number of active flavours
     * @param t: value of the log of the factorisation scale squared
     * @param f: set of distributions at the scale mu
     * @return The r.h.s. of the DGLAP equation.
     */
    Set<T> Derivative(int const& nf, double const& t, Set<T> const& f) const;

    /**
     * @name Setters
     */
    ///@{
    /**
     * @brief Function that sets the reference object at the reference
     * scale using a function of the index and x.
     * @param InDistFunc: the function that returns the distributions.
     */
    void SetInitialDistributions(std::function<double(int const&, double const&)> const& InDistFunc);

    /**
     * @brief Function that sets the reference distribution at the
     * reference scale using a map of the distribution as function of
     * x.
     * @param InDistFunc: the function that returns the distributions.
     */
    void SetInitialDistributions(std::function<std::map<int, double>(double const&)> const& InDistFunc);

    /**
     * @brief Function that sets the reference distribution at the
     * reference scale using a map of the distribution as function of
     * x and Q.
     * @param InDistFunc: the function that returns the distributions.
     * @param mu: the scale at which distributions have to be computed
     */
    void SetInitialDistributions(std::function<std::map<int, double>(double const&, double const&)> const& InDistFunc, double const& mu);
    ///@}
  private:
    std::function<Set<Operator>(int const&, double const&)> const _SplittingFunctions;
    std::function<Set<Operator>(bool const&, int const&)>   const _MatchingConditions;
    std::function<Set<T>(int const&, double const&)>        const _InhomogeneousTerms;
  };
}
