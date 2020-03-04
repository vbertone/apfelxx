//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
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
    InitialiseEvolution(EvolutionSetup const& setup, bool const& WriteGrid = false);

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
     * @brief The InitialiseDglapObject function intialises the
     * relevant objects for the DGLAP evolution and constructs a Dglap
     * object to be used for the evolution.
     */
    void InitialiseDglapObject();

    /**
     * @brief The TabulateEvolution function computes the DGLAP
     * evolution and tabulates the distributions over and (x,Q2)
     * grid. The tabulated distributions can be accessed via the
     * KnotArray() array function.
     * @param InSet: the input set of distributions
     */
    void TabulateEvolution(std::function<std::map<int, double>(double const&, double const&)> const& InSet);

    /**
     * @brief The WriteGridInfo function creates the folder and writes
     * the info file of the LHAPDF grid.
     */
    void WriteGridInfo();

    /**
     * @brief The WriteGrid function dumps to file in the LHAPDF
     * format the actual PDF grid.
     * @note If the writing of the grid is enabled, this function is
     * called every time that the TabulateEvolution function is
     * called. Therefore this cab be used to compute more members of a
     * given set without reinitialising the evolution.
     */
    void WriteGrid();

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
    bool                                         _WriteGrid;    //!< Switch to write LHAPDF grids
    std::unique_ptr<const Grid>                  _g;            //!< x-space grid
    std::function<double(double const&)>         _as;           //!< Strong coupling function
    std::map<int, DglapObjects>                  _DglapObj;     //!< Dglap evolution objects
    std::map<double, std::map<int, LHKnotArray>> _KnotArray;    //!< Object that emulates the KnotArray of LHAPDF to be fed to LHAPDF itself
  };
}
