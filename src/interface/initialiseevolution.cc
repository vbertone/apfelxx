//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/initialiseevolution.h"
#include "apfel/grid.h"
#include "apfel/messages.h"
#include "apfel/alphaqcd.h"
#include "apfel/tabulateobject.h"
#include "apfel/constants.h"
#include "apfel/rotations.h"

#include <iostream>
#include <algorithm>

namespace apfel {
  //_________________________________________________________________________________
  InitialiseEvolution::InitialiseEvolution(EvolutionSetup const& setup):
    _setup(setup)
  {
    // First check the sanity of the setup object
    if (!CheckSetup())
      throw std::runtime_error(error("InitialiseEvolution::InitialiseEvolution","Terminating program due to inconsistent setup."));

    // Now report the evolution set up (if allowed by the vebosity
    // level).
    ReportSetup();

    // Construct vectors of thresholds, masses, and mass reference
    // scales assuming that the light-quark masses are zero.
    const std::vector<double> Thresholds{0, 0, 0, _setup.ThCharm, _setup.ThBottom, _setup.ThTop};
    const std::vector<double> Masses{0, 0, 0, _setup.MCharm, _setup.MBottom, _setup.MTop};

    // Resize mass vectors according to the maximun number of active
    // flavours for the couplings.
    _ThAlpha = Thresholds;
    _MAlpha  = Masses;
    _ThAlpha.resize(_setup.MaxFlavourAlpha);
    _MAlpha.resize(_setup.MaxFlavourAlpha);

    // Resize mass vectors according to the maximun number of active
    // flavours for the couplings.
    _ThDist = Thresholds;
    _MDist  = Masses;
    _ThDist.resize(_setup.MaxFlavourPDFs);
    _MDist.resize(_setup.MaxFlavourPDFs);

    // Initialise the couplings.
    InitialiseCouplings();

    // Initialize Dglap objects.
    InitialiseDglapObject();

    // Make a tabulation just in case the knot array is called straight away.
    TabulateEvolution();    
  }

  //_________________________________________________________________________________
  void InitialiseEvolution::InitialiseCouplings()
  {
    if (_setup.Theory == EvolutionSetup::QCD)
      if (_setup.MassRenScheme == EvolutionSetup::POLE)
	{
	  AlphaQCD a{_setup.AlphaQCDRef, _setup.QQCDRef, _MAlpha, _ThAlpha, _setup.PerturbativeOrder};
	  const TabulateObject<double> Alphas{a, 2 * _setup.nQg, _setup.Qmin - 0.1, _setup.Qmax + 1, _setup.InterDegreeQ};
	  _as = [=] (double const& mu) -> double{ return Alphas.Evaluate(mu); };
	}
  }

  //_________________________________________________________________________________
  void InitialiseEvolution::InitialiseDglapObject()
  {
    // Construct vector of subgrids.
    std::vector<SubGrid> sg;
    for (auto const& gp : _setup.GridParameters)
      if(gp.xgext.empty())
	sg.push_back(SubGrid{gp.nx, gp.xmin, gp.id});
      else
	sg.push_back(SubGrid{gp.xgext, gp.id});

    // Intialise the x-space grid.
    _g = std::unique_ptr<const Grid>(new Grid{sg, _setup.Locked});

    // Integration accuracy.
    const double IntEps = _setup.GaussAccuracy;

    if (_setup.Virtuality == EvolutionSetup::SPACE)
      {
	if (_setup.EvolPolarisation == EvolutionSetup::UNP)
	  _DglapObj = InitializeDglapObjectsQCD(*_g, _MDist, _ThDist, false, IntEps);
	else if (_setup.EvolPolarisation == EvolutionSetup::TRANS)
	  _DglapObj = InitializeDglapObjectsQCDtrans(*_g, _MDist, _ThDist, false, IntEps);
      }
    else if (_setup.Virtuality == EvolutionSetup::TIME)
      {
	if (_setup.EvolPolarisation == EvolutionSetup::UNP)
	  _DglapObj = InitializeDglapObjectsQCDT(*_g, _MDist, _ThDist, false, IntEps);
	else if (_setup.EvolPolarisation == EvolutionSetup::TRANS)
	  _DglapObj = InitializeDglapObjectsQCDTtrans(*_g, _MDist, _ThDist, false, IntEps);
      }
  }

  //_________________________________________________________________________________
  void InitialiseEvolution::TabulateEvolution()
  {
    // Construct the Dglap object.
    std::unique_ptr<Dglap<Distribution>> EvolvedDists = BuildDglap(_DglapObj, _setup.InSet, _setup.Q0, _setup.PerturbativeOrder, _as);

    // Tabulate distributions
    const TabulateObject<Set<Distribution>> TabulatedDists{*EvolvedDists, _setup.nQg, _setup.Qmin, _setup.Qmax, _setup.InterDegreeQ};

    // Get Q-grid from the tabulated object.
    const std::vector<double> qg = TabulatedDists.GetQGrid();

    // Get threshold indices.
    const std::vector<int> tind = TabulatedDists.GetThesholdIndices();

    // Get Set of distributions on the Q-grid
    const std::vector<Set<Distribution>> xfg = TabulatedDists.GetQGridValues();

    // Run over the threshold indices.
    for (int i = 0; i < (int) tind.size(); i++)
      {
	// Threshold index
	const int ti = tind[i];

	// Retrieve distributions at the threshold.
	const std::map<int, Distribution> tdist = xfg[ti].GetObjects();

	// Run over distributions.
	std::map<int, LHKnotArray> LHKnotArrayNF;
	for (auto const& d : tdist)
	  {
	    LHKnotArray ka;

	    // x-space grid
	    ka.xs = d.second.GetGrid().GetJointGrid().GetGrid();

	    // Remove nodes above one.
	    ka.xs.resize(ka.xs.size()-d.second.GetGrid().GetJointGrid().InterDegree());

	    // Q2-space (sub)grid
	    std::vector<double> q2;
	    for (int iq = ti; iq < tind[i] + 1;iq++)
	      q2.push_back(qg[iq] * qg[iq]);
	    ka.q2s = q2;

	    // Get sizes of both x and q2 grids.
	    const int xsize = ka.xs.size();
	    const int q2size = ka.q2s.size();

	    // Now fill in vector of distributions.
	    std::vector<double> xf(xsize * q2size);
	    for (int iq = ti; iq < tind[i] + 1;iq++)
	      {
		const std::vector<double> xfx = xfg[iq].GetObjects().at(d.first).GetDistributionJointGrid();
		for (int ix = 0; ix < (int) ka.xs.size(); ix++)
		  xf[ix * q2size + iq - ti] = xfx[ix];
	      }
	    ka.xfs = xf;

	    LHKnotArrayNF.insert({d.first, ka});
	  }
	_KnotArray.insert({qg[ti] * qg[ti], LHKnotArrayNF});
      }
  }

  //_________________________________________________________________________________
  bool InitialiseEvolution::CheckSetup() const
  {
    // Initialise switch to true.
    bool passed = true;

    // Check initial scale.
    if (_setup.Q0 <= 0.5)
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "The initial evolution scale must be larger than 0.5 GeV.") << std::endl;
	passed = false;
      }

    // Check the x-space grid.
    for (auto const& gp : _setup.GridParameters)
      {
	// Interpolation degree.
	if (gp.id <= 0)
	  {
	    std::cout << error("InitialiseEvolution::CheckSetup", "The interpolation degree of each subgrid must be positive.") << std::endl;
	    passed = false;
	  }

	// Internal grid.
	if(gp.xgext.empty())
	  {
	    if (gp.nx < gp.id + 1)
	      {
		std::cout << error("InitialiseEvolution::CheckSetup", "The number of nodes of each subgrid must be bigger than the interpolation degree plus one.") << std::endl;
		passed = false;
	      }

	    if (gp.xmin <= 0 || gp.xmin >= 1)
	      {
		std::cout << error("InitialiseEvolution::CheckSetup", "The lower bound xmin of each subgrid must be such that 0 < xmin < 1.") << std::endl;
		passed = false;	
	      }
	  }
	else
	  {
	    // Esternal grid.
	    const std::vector<double> v = gp.xgext;

	    if (!std::is_sorted(v.begin(),v.end()))
	      {
		std::cout << error("InitialiseEvolution::CheckSetup", "The grid vector of a subgrid is not sorted.") << std::endl;
		passed = false;
	      }

	    if (v.front() <= 0 || v.front() >= 1)
	      {
		std::cout << error("InitialiseEvolution::CheckSetup", "The first element of each subgrid vector must be between 0 and 1.") << std::endl;
		passed = false;
	      }

	    if (std::abs(v.back() - 1) > eps10)
	      {
		std::cout << error("InitialiseEvolution::CheckSetup", "The last element of each subgrid vector must be equal to 1.") << std::endl;
		passed = false;
	      }
	    
	    if ((int) v.size() < gp.id+1)
	      {
		std::cout << error("InitialiseEvolution::CheckSetup", "The number of nodes of each subgrid must be bigger than the interpolation degree plus one.") << std::endl;
		passed = false;
	      }
	  }
      }

    // Check the Q-space grid.
    if (_setup.InterDegreeQ <= 0)
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "The interpolation degree of the grid in Q must be positive.") << std::endl;
	passed = false;
      }

    if (_setup.nQg < _setup.InterDegreeQ + 1)
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "The number of nodes of the grid in Q must be bigger than the interpolation degree plus one.") << std::endl;
	passed = false;
      }

    if (_setup.Lambda <= 0)
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "The Lambda parameter of the grid in Q must be positive.") << std::endl;
	passed = false;
      }

    if (_setup.Qmin <= _setup.Lambda)
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "The lower bound of the grid in Q must be larger than Lambda.") << std::endl;
	passed = false;
      }

    if (_setup.Qmax <= _setup.Qmin)
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "The upper bound of the grid in Q must larger than the lower bound.") << std::endl;
	passed = false;
      }

    // Check flavour scheme
    if (_setup.FlavourScheme != EvolutionSetup::VFNS &&
	_setup.FlavourScheme != EvolutionSetup::FFNS)
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "Unknown flavour scheme.") << std::endl;
	passed = false;
      }

    // Check number of flavours in the FFNS (if relevant)
    if (_setup.FlavourScheme == EvolutionSetup::FFNS)
      if (_setup.Nf_FF < 3 || _setup.Nf_FF > 6)
	{
	  std::cout << error("InitialiseEvolution::CheckSetup", "Number of active flavours in the FFNS out of range.") << std::endl;
	  passed = false;
	}

    // Check theory.
    if (_setup.Theory != EvolutionSetup::QCD)
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "Unknown theory.") << std::endl;
	passed = false;
      }

    // Check virtuality
    if (_setup.Virtuality != EvolutionSetup::SPACE &&
	_setup.Virtuality != EvolutionSetup::TIME)
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "Unknown virtuality.") << std::endl;
	passed = false;
      }

    // Check evolution polarisation.
    if (_setup.EvolPolarisation != EvolutionSetup::UNP &&
	_setup.EvolPolarisation != EvolutionSetup::POL &&
	_setup.EvolPolarisation != EvolutionSetup::TRANS)
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "Unknown evolution polarisation.") << std::endl;
	passed = false;
      }

    // Check perturbative order.
    if (_setup.PerturbativeOrder < 0 || _setup.PerturbativeOrder > 3)
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "Perturbative order out of range.") << std::endl;
	passed = false;
      }

    // Make sure that for each vituality and polarisation, the maximum
    // perturbative order in not exceeded.
    if (_setup.Virtuality == EvolutionSetup::SPACE)
      {
	if (_setup.EvolPolarisation == EvolutionSetup::POL)
	  {
	    std::cout << error("InitialiseEvolution::CheckSetup", "Space-like polarised evolution not implemented yet.") << std::endl;
	    passed = false;
	    //if (_setup.PerturbativeOrder > 2)
	    //  {
	    //     std::cout << error("InitialiseEvolution::CheckSetup", "Space-like polarised evolution implemented up to NNLO.") << std::endl;
	    //     passed = false;
	    //  }
	  }
	else if (_setup.EvolPolarisation == EvolutionSetup::TRANS)
	  {
	    if (_setup.PerturbativeOrder > 1)
	      {
		std::cout << error("InitialiseEvolution::CheckSetup", "Space-like transverse evolution implemented up to NLO.") << std::endl;
		passed = false;
	      }
	  }
      }
    else if (_setup.Virtuality == EvolutionSetup::TIME)
      {
	if (_setup.EvolPolarisation == EvolutionSetup::UNP)
	  {
	    if (_setup.PerturbativeOrder > 2)
	      {
		std::cout << error("InitialiseEvolution::CheckSetup", "Time-like polarised evolution implemented up to NNLO.") << std::endl;
		passed = false;
	      }
	  }
	if (_setup.EvolPolarisation == EvolutionSetup::POL)
	  {
	    std::cout << error("InitialiseEvolution::CheckSetup", "Time-like polarised evolution not implemented yet.") << std::endl;
	  }
	else if (_setup.EvolPolarisation == EvolutionSetup::TRANS)
	  {
	    if (_setup.PerturbativeOrder > 1)
	      {
		std::cout << error("InitialiseEvolution::CheckSetup", "Time-like transverse evolution implemented up to NLO.") << std::endl;
		passed = false;
	      }
	  }
      }

    // Check coupling evolution.
    if (_setup.CouplingEvolution != EvolutionSetup::exact &&
	_setup.CouplingEvolution != EvolutionSetup::expanded)
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "Unknown coupling evolution.") << std::endl;
	passed = false;
      }

    // Check coupling evolution.
    if (_setup.PDFEvolution != EvolutionSetup::exactmu &&
	_setup.PDFEvolution != EvolutionSetup::exactalpha &&
	_setup.PDFEvolution != EvolutionSetup::expandalpha &&
	_setup.PDFEvolution != EvolutionSetup::truncated)
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "Unknown PDF evolution.") << std::endl;
	passed = false;
      }

    // Check maximum number of active flavours in the couplings.
    if (_setup.MaxFlavourAlpha < 3 || _setup.MaxFlavourAlpha > 6)
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "Maximum number of active flavours allowed in the coupling evolution out of range.") << std::endl;
	passed = false;
      }
    // Check maximum number of active flavours in PDFs.
    if (_setup.MaxFlavourPDFs < 3 || _setup.MaxFlavourPDFs > 6)
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "Maximum number of active flavours allowed in the coupling evolution out of range.") << std::endl;
	passed = false;
      }

    // Check mass renormalisation scheme.
    if (_setup.MassRenScheme != EvolutionSetup::POLE &&
	_setup.MassRenScheme != EvolutionSetup::MSBAR)
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "Unknown mass renormalisation scheme.") << std::endl;
	passed = false;
      }

    // Check Gauss integration accuracy.
    if (_setup.GaussAccuracy <= 0 || _setup.GaussAccuracy > 1)
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "The dGauss integration accuracy must me in the (0:1] range.") << std::endl;
	passed = false;
      }

    // Check that heavy-quark thresholds, masses, and reference scales
    // are sorted.
    const std::vector<double> Thresholds{_setup.ThCharm, _setup.ThBottom, _setup.ThTop};
    if (!std::is_sorted(Thresholds.begin(),Thresholds.end()))
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "The heavy-quark thresholds are not sorted.") << std::endl;
	passed = false;
      }

    const std::vector<double> Masses{_setup.MCharm, _setup.MBottom, _setup.MTop};
    if (!std::is_sorted(Masses.begin(),Masses.end()))
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "The heavy-quark masses are not sorted.") << std::endl;
	passed = false;
      }

    // Check renormalisarion / factorisation scale ratio.
    if (_setup.RenFacRatio <= 0)
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "The ratio between renormalisation and factorisation scales cannot be negative.") << std::endl;
	passed = false;
      }

    return passed;
  }

  //_________________________________________________________________________________
  void InitialiseEvolution::ReportSetup() const
  {
    // Apfel banner.
    Banner();

    // String with the report of parameters
    std::string report = "\n\nCurrent evolution setup:\n";

    // Starting scale.
    report += "- Evolution starting scale: " + std::to_string(_setup.Q0) + " GeV\n";

    // Perturbative order.
    report += "- Perturbative order of the evolution: N" + std::to_string(_setup.PerturbativeOrder) + "LO\n";

    // Flavour scheme.
    report += "- Flavous scheme: ";
    if (_setup.FlavourScheme == EvolutionSetup::VFNS)
      report += "VFNS\n";
    else if (_setup.FlavourScheme == EvolutionSetup::FFNS)
      report += "FFNS with " + std::to_string(_setup.Nf_FF) + " active flavours\n";

    // Evolution theory.
    report += "- Evolution theory: ";
    if (_setup.Theory == EvolutionSetup::QCD)
      report += "QCD\n";
    else if (_setup.Theory == EvolutionSetup::QCD_QED)
      report += "QCD + QED\n";

    // Virtuality.
    report += "- Evolution virtuality: ";
    if (_setup.Virtuality == EvolutionSetup::SPACE)
      report += "Space-like (PDFs)\n";
    else if (_setup.Virtuality == EvolutionSetup::TIME)
      report += "Time-like (FFs)\n";

    // Polarisation.
    report += "- Evolution polarisation: ";
    if (_setup.EvolPolarisation == EvolutionSetup::UNP)
      report += "Unpolarised\n";
    else if (_setup.EvolPolarisation == EvolutionSetup::POL)
      report += "longitudinally polarised\n";
    else if (_setup.EvolPolarisation == EvolutionSetup::TRANS)
      report += "transversely polarised\n";

    // AlphaQCD reference value.
    report += "- QCD coupling reference value: AlphaQCD(" + std::to_string(_setup.QQCDRef) + " GeV) = " + std::to_string(_setup.AlphaQCDRef) + "\n";

    // AlphaQED reference value.
    if (_setup.Theory == EvolutionSetup::QCD_QED)
      report += "- QED coupling reference value: AlphaQED(" + std::to_string(_setup.QQEDRef) + " GeV) = " + std::to_string(_setup.AlphaQEDRef) + "\n";

    // Coupling evolution.
    report += "- Coupling evolution: ";
    if (_setup.CouplingEvolution == EvolutionSetup::exact)
      report += "Exact";
    else if (_setup.CouplingEvolution == EvolutionSetup::expanded)
      report += "Expanded";
    report += " with maximum " + std::to_string(_setup.MaxFlavourAlpha) + " active flavours\n";

    // PDF evolution.
    report += "- PDF evolution: ";
    if (_setup.PDFEvolution == EvolutionSetup::exactmu)
      report += "Exact in mu";
    else if (_setup.PDFEvolution == EvolutionSetup::exactalpha)
      report += "Exact in alpha";
    else if (_setup.PDFEvolution == EvolutionSetup::expandalpha)
      report += "Expanded in alpha";
    else if (_setup.PDFEvolution == EvolutionSetup::truncated)
      report += "Truncated";
    report += " with maximum " + std::to_string(_setup.MaxFlavourPDFs) + " active flavours\n";

    // Ren. / fact. scales ratio.
    report += "- muR / mF: " + std::to_string(_setup.RenFacRatio) + "\n";

    // Heavy-quark masses and thresholds.
    if (_setup.MassRenScheme == EvolutionSetup::POLE)
      {
	report += "- Pole heavy-quark masses:\n";
	report += "  + Mc = " + std::to_string(_setup.MCharm)  + " GeV\n";
	report += "  + Mb = " + std::to_string(_setup.MBottom) + " GeV\n";
	report += "  + Mt = " + std::to_string(_setup.MTop)    + " GeV\n";
	report += "- Pole heavy-quark threholds:\n";
	report += "  + muc = " + std::to_string(_setup.ThCharm)  + " GeV\n";
	report += "  + mub = " + std::to_string(_setup.ThBottom) + " GeV\n";
	report += "  + mut = " + std::to_string(_setup.ThTop)    + " GeV\n";
      }
    if (_setup.MassRenScheme == EvolutionSetup::MSBAR)
      {
	report += "- MSbar heavy-quark masses:\n";
	report += "  + mc(mc) = " + std::to_string(_setup.MCharm)  + " GeV\n";
	report += "  + mb(mb) = " + std::to_string(_setup.MBottom) + " GeV\n";
	report += "  + mt(mt) = " + std::to_string(_setup.MTop)    + " GeV\n";
	report += "- MSbar heavy-quark threholds:\n";
	report += "  + muc = " + std::to_string(_setup.ThCharm)  + " GeV\n";
	report += "  + mub = " + std::to_string(_setup.ThBottom) + " GeV\n";
	report += "  + mut = " + std::to_string(_setup.ThTop)    + " GeV\n";
      }

    // Tau mass.
    report += "- Tau lepton mass: " + std::to_string(_setup.TauMass) + " GeV\n";

    // Integration accuracy.
    report += "- Relative Gauss integration accuracy: " + std::to_string(_setup.GaussAccuracy) + "\n\n";

    // Report x-grid parameters.
    report += "- Grid in x: " + std::to_string(_setup.GridParameters.size()) + " subgrids in x found with the following parameters:\n";
    for (auto const& gp : _setup.GridParameters)
      if (gp.xgext.empty())
	report += "   + internal grid with " + std::to_string(gp.nx)
	  + " nodes in the range [" + std::to_string(gp.xmin)
	  + ":1] with interpolation degree " + std::to_string(gp.id) + "\n";
      else
	report += "   + ixternal grid with " + std::to_string(gp.xgext.size())
	  + " nodes in the range [" + std::to_string(gp.xgext.front())
	  + ":" + std::to_string(gp.xgext.back())
	  + "] with interpolation degree " + std::to_string(gp.id) + "\n";

    if (_setup.Locked)
      report += "   The sugrids are locked\n\n";

    // Report Q-grid parameters.
    report += "- Grid in Q: " + std::to_string(_setup.nQg)
      + " nodes in the range [" + std::to_string(_setup.Qmin)
      + ":" + std::to_string(_setup.Qmax) + "] GeV"
      + " with interpolation degree " + std::to_string(_setup.InterDegreeQ) + "\n";

    // Print report.
    info("InitialiseEvolution::ReportSetup", report);
  }
}
