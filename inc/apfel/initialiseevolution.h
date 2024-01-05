//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/evolutionsetup.h"
#include "apfel/tabulateobject.h"
#include "apfel/dglapbuilder.h"
#include "apfel/config.h"

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
     * @param WriteGrid: switch to enable the writing of grids in the LHAPDF format (default: false)
     * @param GridHeader: part of the LHAPDF grid header that can be set externally (default: empty = use default)
     */
    InitialiseEvolution(EvolutionSetup const& setup, bool const& WriteGrid = false, std::string const& GridHeader = "");

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
     * @brief Function that returns the full set of distributions
     * tabulated on the (x,Q2) grid.
     * @return the _KnotArray attribute
     */
    std::map<double, std::map<int, LHKnotArray>> KnotArray() const { return _KnotArray; }

    /**
     * @brief Function that returns the full set of distributions as a
     * TabulateObject object.
     * @return the _KnotArray attribute
     */
    TabulateObject<Set<Distribution>> TabulatedDistributions() const { return *_TabulatedDists; }

    /**
     * @brief Function that returns the EvolutionSetup object.
     * @return the the _setup object
     */
    EvolutionSetup GetEvolutionSetup() const { return _setup; }

  private:
    EvolutionSetup                                           _setup;          //!< Evolution setup object
    bool                                                     _WriteGrid;      //!< Switch to write LHAPDF grids
    std::string                                              _GridHeader;     //!< Part of the LHAPDF grid header that can be set externally (the format is resposibility of the user)
    std::unique_ptr<const Grid>                              _g;              //!< x-space grid
    std::function<double(double const&)>                     _as;             //!< Strong coupling function
    std::map<int, DglapObjects>                              _DglapObj;       //!< Dglap evolution objects
    std::map<double, std::map<int, LHKnotArray>>             _KnotArray;      //!< Object that emulates the KnotArray of LHAPDF to be fed to LHAPDF itself
    std::unique_ptr<const TabulateObject<Set<Distribution>>> _TabulatedDists; //!< Tabulated distributions
  };
}

#if WITH_LHAPDF == 1

#include <LHAPDF/LHAPDF.h>
#include <LHAPDF/GridPDF.h>

/**
 * @brief Function that constructs a pointer to an LHAPDF::PDF object
 * using an apfel::InitialiseEvolution object as an input in the place
 * of a standard LHAPDF grid.
 * @param ev: APFEL++ InitialiseEvolution object
 * @return a pointer to an LHAPDF::PDF object
 */
LHAPDF::PDF* mkPDF(apfel::InitialiseEvolution const& ev)
{
  // Get knot array from APFEL
  const std::map<double, std::map<int, apfel::LHKnotArray>> ka = ev.KnotArray();

  // Knot array to be fed to LHAPDF
  LHAPDF::KnotArray data;

  // Loop over Q subgrids to accumulate x-grid, q2-grid, and flavour
  // IDs.
  for (auto const& sg : ka)
    {
      // Get x-grid from first flavour of the first subgrid
      if (data.xs().empty())
        data.setxknots() = sg.second.begin()->second.xs;

      // Accumulate q2 grid using the first flavour
      data.setq2knots().insert(data.setq2knots().end(), sg.second.begin()->second.q2s.begin(), sg.second.begin()->second.q2s.end());

      // Fill in the flavour ID vector (use 21 for the gluon)
      if (data.setPids().empty())
        for (auto const& id : sg.second)
          data.setPids().push_back((id.first == 0 ? 21 : id.first));
    }

  // Set up the knots of the Knotarray
  data.setShape() = std::vector<size_t> {data.xs().size(), data.q2s().size(), data.setPids().size()};
  data.fillLogKnots();
  data.initPidLookup();

  // Set size of data vector
  data.setGrid().resize(data.shape(0) * data.shape(1) * data.shape(2), 0.);

  // Loop over Q subgrids to accumulate data
  int nqprev = 0;
  for (auto const& sg : ka)
    {
      // Initialise flavour index
      int ifl = 0;

      // Get q2-subgrid size
      const int q2size = (int) sg.second.begin()->second.xfs.size() / data.shape(0);

      // Loop over flavours
      for (auto const& id : sg.second)
        {
          // Get vector of distributions
          const std::vector<double> f = id.second.xfs;

          // Loop over x-grid
          for (int ix = 0; ix < (int) data.shape(0); ix++)
            {
              int iq = nqprev;
              // Loop over q2-sugrid
              for (int iqs = 0; iqs < q2size; iqs++)
                data.setGrid()[ix * data.shape(2) * data.shape(1) + iq++ * data.shape(2) + ifl] = f[ix * q2size + iqs];
            }
          ifl++;
        }
      nqprev += q2size;
    }

  // Fill in alpha_s vector
  std::vector<double> als;
  for (auto const& q2 : data.q2s())
    als.push_back(ev.Alphas(sqrt(q2)));

  // Construct alpha_s object
  LHAPDF::AlphaS_Ipol* as = new LHAPDF::AlphaS_Ipol{};
  as->setQ2Values(data.q2s());
  as->setAlphaSValues(als);

  // Define an LHAPDF GridPDF object
  LHAPDF::GridPDF* dist = new LHAPDF::GridPDF{};

  // Pass knot arrays to LHAPDF
  dist->Data() = data;

  // Set interpolator
  dist->setInterpolator((std::string) "logcubic");

  // Set extrapolator
  dist->setExtrapolator((std::string) "continuation");

  // Set flavours
  dist->setFlavors(data.setPids());

  // Set alpha_s
  dist->setAlphaS(as);

  // Set scale bounds
  dist->info().set_entry("QMin", sqrt(data.q2s().front()));
  dist->info().set_entry("QMax", sqrt(data.q2s().back()));

  // Set x bounds
  dist->info().set_entry("XMin", data.xs().front());
  dist->info().set_entry("XMax", data.xs().back());

  // Set quark masses and thresholds
  const std::vector<std::string> Qnames{"Down", "Up", "Strange", "Charm", "Bottom", "Top"};
  const std::vector<double> trhs = ev.GetEvolutionSetup().Thresholds;
  for (int iq = 0; iq < (int) Qnames.size(); iq++)
    {
      dist->info().set_entry("M" + Qnames[iq], (iq < (int) trhs.size() ? trhs[iq] : 1e8 + iq));
      dist->info().set_entry("Threshold" + Qnames[iq], (iq < (int) trhs.size() ? trhs[iq] : 1e8 + iq));
    }

  // Return object
  return dist;
}

#endif
