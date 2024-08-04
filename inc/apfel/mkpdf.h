//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/config.h"

#if WITH_LHAPDF == 1

#include "apfel/initialiseevolution.h"

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
