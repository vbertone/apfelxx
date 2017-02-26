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
#include "apfel/expression.h"
#include "apfel/alphaqcd.h"
#include "apfel/tools.h"

#include <array>
#include <functional>

using std::array;
using std::function;

namespace apfel
{
  /**
   * @brief The DglapQCD class.
   *
   * A specialization class of the MatchedEvolution class
   * for the computation of the DglapQCD evolution.
   */
  class DglapQCD: public MatchedEvolution<Set<Distribution>>
  {
  public:

    DglapQCD() =  delete;

    /**
     * @brief DglapQCD default constructor.
     */
    DglapQCD(Grid                         const& g,
	     function<double(int,double)> const& InPDFs,
	     double                       const& MuDistRef,
	     double                       const& AlphaRef,
	     double                       const& MuAlphaRef,
	     vector<double>               const& Masses,
	     vector<double>               const& Thresholds,
	     int                          const& PertOrder,
	     int                          const& nstep = 10,
	     double                       const& xi = 1);

    /**
     * @brief Function for the computation of the coupling given nf. This function can be overriden.
     * @param nf number of active flavours.
     * @param sd0 starting value of the coupling.
     * @param mu02 initial squared scale.
     * @param mu2 final squared scale.
     * @return value of the coupling at mu2.
     */
    Set<Distribution> EvolveObject(int const& nf, Set<Distribution> const& sd0, double const& mu02, double const& mu2) const;

    /**
     * @brief Function for the computation of the matching. This function can be overriden.
     * @param Up tells whether the matching is upward or not (downward).
     * @param Coup value of the coupling to be matched.
     * @param LogKth value of ln( muth2 / m2 ).
     * @return the matched value of the coupling.
     */
    Set<Distribution> MatchObject(bool const& Up, Set<Distribution> const& sd, double const& LogKth) const;

    /**
     * @brief Function for the computation of the full QCD beta function.
     * @param as value of the coupling.
     * @param nf number of active flavours.
     * @return the value of the beta function.
     */
    Set<Distribution> Derivative(int const& nf, double const& as, Set<Distribution> const& f) const;

    /**
     * @brief Function that returns the perturbative order.
     */
    int const& GetPerturbativeOrder() const { return _PertOrder; }

    /**
     * @brief Function that returns the number of steps.
     */
    int const& GetNumberOfSteps()     const { return _nstep; }

    /**
     * @brief Function that returns AlphaQCD.
     */
    AlphaQCD const& GetAlphaQCD()     const { return _as; }

    //_________________________________________________________________________________
    /**
     * @brief The QCD splitting function classes
     */
    class P0ns: public Expression
    {
    public:
      P0ns();
      double Regular(double const& x)  const;
      double Singular(double const& x) const;
      double Local(double const& x)    const;
    };

    class P0qg: public Expression
    {
    public:
      P0qg();
      double Regular(double const& x)  const;
    };

    class P0gq: public Expression
    {
    public:
      P0gq();
      double Regular(double const& x)  const;
    };

    class P0gg: public Expression
    {
    public:
      P0gg(int const& nf);
      double Regular(double const& x)  const;
      double Singular(double const& x) const;
      double Local(double const& x)    const;
    private:
      int const _nf;
    };

    /**
     * @brief A very simple example of ConvolutionMap derivation.
     *
     * This class, following the derivation procedure from ConvolutionMap
     * implements the Basis enumerator with custom tags for the objects.
     */
    class EvolutionBasis: public ConvolutionMap
    {
    public:
      /**
       * @brief The map enums
       */
      enum Operand: int {PNSP, PNSM, PNSV, PQQ, PQG, PGQ, PGG};
      enum Object:  int {GLUON, SIGMA, VALENCE, T3, V3, T8, V8, T15, V15, T24, V24, T35, V35};

      /**
       * @brief The class constructors
       */
    EvolutionBasis();
    EvolutionBasis(int const& nf);
    };

    /**
     * @brief Class for the distributions
     */
    class PDF: public Distribution
    {
    public:
      PDF(Grid const& gr, function<double(int,double)> const& InPDFs, int const& ipdf);
    };

  private:
    function<double(int,double)>     const _InPDFs;
    int                              const _PertOrder;
    int                              const _nstep;
    AlphaQCD                         const _as;
    unordered_map<int,Set<Operator>>       _SplittingFunctions;
  };

}
