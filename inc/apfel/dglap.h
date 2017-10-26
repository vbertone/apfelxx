//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/matchedevolution.h"
#include "apfel/set.h"
#include "apfel/distribution.h"
#include "apfel/operator.h"

#include <functional>
#include <map>

using std::function;
using std::map;

namespace apfel
{
  /**
   * @brief The Dglap class.
   *
   * A specialization class of the MatchedEvolution class for the
   * computation of the DGLAP evolution.
   */
  template<class T>
  class Dglap: public MatchedEvolution<Set<T>>
  {
  public:

    Dglap() =  delete;

    /**
     * @brief Dglap default constructor.
     */
    Dglap(function<Set<Operator>(int const&,double const&)> const& SplittingFunctions,
	  function<Set<Operator>(bool const&,int const&)>   const& MatchingConditions,
	  Set<T>                                            const& ObjRef,
	  double                                            const& MuRef,
	  vector<double>                                    const& Thresholds,
	  int                                               const& nsteps = 10);

    /**
     * @brief
     * @param
     * @return
     */
    Set<T> MatchObject(bool const& Up, int const& nf, Set<T> const& sd) const;

    /**
     * @brief
     * @param
     * @return
     */
    Set<T> Derivative(int const& nf, double const& mu, Set<T> const& f) const;

    /**
     * @brief Function that sets the reference object at the reference
     * scale using a function of the index and x.
     */
    void SetInitialDistributions(function<double(int const&, double const&)> const& InDistFunc);

    /**
     * @brief Function that sets the reference distribution at the
     * reference scale using a map of the distribution as function of
     * x.
     */
    void SetInitialDistributions(function<map<int,double>(double const&)> const& InDistFunc);

  private:
    function<Set<Operator>(int const&,double const&)> _SplittingFunctions;
    function<Set<Operator>(bool const&,int const&)>   _MatchingConditions;
  };
}
