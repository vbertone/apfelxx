//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/evolutionsetup.h"
#include "apfel/tabulateobject.h"
#include "apfel/dglapbuilder.h"

namespace apfel
{
  /**
   * @brief The LHKnotArray structure emulates the KnotArray1F class
   * of LHAPDF and contains the grids in x, Q2 (only a given subgrid),
   * and one single distribution tabulated on the (x,Q) bidimensional
   * grid.
   */
  struct LHKnotArray
  {
    std::vector<double> xs;     //!< List of x knots
    std::vector<double> q2s;    //!< List of Q2 knots
    std::vector<double> xfs;    //!< List of xf values across the 2D knot array, stored as a strided [ix][iQ2] 1D array
  };

  /**
   * @brief The InitialiseEvolution performs all the operations to
   * initialise a DGLAP evolution using an EvolutionSetup object to
   * retrieve the relevant information. This class also provides the
   * necessary functions to access the evolved distributions,
   * coupling(s), and masses.
   */
  class InitialiseEvolution
  {
  public:
    /**
     * @brief The InitialiseEvolution constructor.
     * @param setup: the EvolutionSetup data structure encapsulate the evolution parameters
     */
    InitialiseEvolution(EvolutionSetup const& setup);

    /**
     * @brief The CheckSetup function checks that the input setup is
     * meaningful and compatible with the current capabilities of the
     * code.
     * @return true or false according to whether the check is succesful or not
     * @note Of course, not all possible checks are implemented. The
     * user has to be careful when modifying the evolution setup and
     * make sure that the settings are reasonable.
     */
    bool CheckSetup() const;

    /**
     * @brief The ReportSetup function reports the parameters in a
     * human readable format.
     */
    void ReportSetup() const;

    /**
     * @brief The InitialiseCouplings function intialises and
     * tabulates the running coupling(s).
     */
    void InitialiseCouplings();

    /**
     * @brief The InitialiseDglapObject intialises the relevant
     * objects for the DGLAP evolution and constructs a Dglap object
     * to be used for the evolution.
     */
    void InitialiseDglapObject();

    /**
     * @brief The TabulateEvolution computes the DGLAP evolution and
     * tabulates the distributions over and (x,Q2) grid. The tabulated
     * distributions can be accessed via the KnotArray() array
     * function.
     */
    void TabulateEvolution();

    /**
     * @brief The ResetInitialDistributions resets the initial set of
     * distributions without the need of generating a new
     * EvolutionSetup object.
     */
    void ResetInitialDistributions(std::function<std::map<int,double>(double const&, double const&)> InSet) { _setup.InSet = InSet; }

    /**
     * @brief Function that returns the evolved strong coupling.
     * @param mu: the value of the renormalisation scale &mu;<SUB>R</SUB>
     * @return the evolved strong coupling
     */
    double Alphas(double const& mu) const { return _as(mu); }

    /**
     * @brief Function that returns the full set of distributions tabulated of the (x,Q2) grid.
     * @return the _KnotArray attribute
     */
    std::map<double, std::map<int, LHKnotArray>> KnotArray() const { return _KnotArray; }

  private:
    EvolutionSetup                               _setup;        //!< Evolution setup object
    std::unique_ptr<const Grid>                  _g;            //!< x-space grid
    std::vector<double>                          _ThAlpha;      //!< Vector of quark thresholds for the couplings
    std::vector<double>                          _MAlpha;       //!< Vector of quark masses for the couplings
    std::vector<double>                          _ThDist;       //!< Vector of quark thresholds for the distributions
    std::vector<double>                          _MDist;        //!< Vector of quark masses for the distributions
    std::function<double(double const&)>         _as;           //!< Strong coupling function
    std::map<int, DglapObjects>                  _DglapObj;     //!< Dglap evolution objects
    std::map<double, std::map<int, LHKnotArray>> _KnotArray;    //!< Object that emulates the KnotArray of LHAPDF to be fet to LHAPDF itself
  };
}
