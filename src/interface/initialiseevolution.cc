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
#include <sys/stat.h>
#include <fstream>

namespace apfel {
  //_________________________________________________________________________________
  InitialiseEvolution::InitialiseEvolution(EvolutionSetup const& setup, bool const& WriteGrid):
    _setup(setup),
    _WriteGrid(WriteGrid)
  {
    // Check the sanity of the setup object.
    if (!CheckSetup())
      throw std::runtime_error(error("InitialiseEvolution::InitialiseEvolution","Terminating program due to inconsistent setup."));

    // Report the evolution set up (if allowed by the vebosity level).
    ReportSetup();

    // Initialise the couplings.
    InitialiseCouplings();

    // Initialize Dglap objects.
    InitialiseDglapObject();

    // If the production of grids has been enabled, create a folded
    // with the name of the setup and write the info file.
    if (WriteGrid)
      WriteGridInfo();

    // Do the tabulation.
    TabulateEvolution(setup.InSet);
  }

  //_________________________________________________________________________________
  void InitialiseEvolution::InitialiseCouplings()
  {
    if (_setup.Theory == EvolutionSetup::QCD)
      if (_setup.MassRenScheme == EvolutionSetup::POLE)
	{
	  AlphaQCD a{_setup.AlphaQCDRef, _setup.QQCDRef, _setup.Masses, _setup.Thresholds, _setup.PerturbativeOrder};
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
	  _DglapObj = InitializeDglapObjectsQCD(*_g, _setup.Masses, _setup.Thresholds, false, IntEps);
	else if (_setup.EvolPolarisation == EvolutionSetup::TRANS)
	  _DglapObj = InitializeDglapObjectsQCDtrans(*_g, _setup.Masses, _setup.Thresholds, false, IntEps);
      }
    else if (_setup.Virtuality == EvolutionSetup::TIME)
      {
	if (_setup.EvolPolarisation == EvolutionSetup::UNP)
	  _DglapObj = InitializeDglapObjectsQCDT(*_g, _setup.Masses, _setup.Thresholds, false, IntEps);
	else if (_setup.EvolPolarisation == EvolutionSetup::TRANS)
	  _DglapObj = InitializeDglapObjectsQCDTtrans(*_g, _setup.Masses, _setup.Thresholds, false, IntEps);
      }
  }

  //_________________________________________________________________________________
  void InitialiseEvolution::TabulateEvolution(std::function<std::map<int,double>(double const&, double const&)> const& InSet)
  {
    // Clear the knot array map before filling it (in case this
    // function was already called before).
    _KnotArray.clear();

    // Construct the Dglap object.
    std::unique_ptr<Dglap<Distribution>> EvolvedDists = BuildDglap(_DglapObj, InSet, _setup.Q0, _setup.PerturbativeOrder, _as);

    // Tabulate distributions.
    const TabulateObject<Set<Distribution>> TabulatedDists{*EvolvedDists, _setup.nQg, _setup.Qmin, _setup.Qmax, _setup.InterDegreeQ};

    // Get Q-grid from the tabulated object.
    const std::vector<double> qg = TabulatedDists.GetQGrid();

    // Get threshold indices.
    const std::vector<int> tind = TabulatedDists.GetThesholdIndices();

    // Get Set of distributions on the Q-grid.
    const std::vector<Set<Distribution>> xfg = TabulatedDists.GetQGridValues();

    // Run over the threshold indices. Skipe the last because it is
    // the last point of the grid in Q.
    for (int i = 0; i < (int) tind.size() - 1; i++)
      {
	// Threshold index.
	const int ti = tind[i];

	// Retrieve distributions at the threshold and rotate
	// distribution into the physical basis.
	const std::map<int, Distribution> tdist = QCDEvToPhys(xfg[ti].GetObjects());

	// Run over distributions.
	std::map<int, LHKnotArray> LHKnotArrayNF;
	for (auto const& d : tdist)
	  {
	    LHKnotArray ka;

	    // x-space grid.
	    ka.xs = d.second.GetGrid().GetJointGrid().GetGrid();

	    // Remove nodes above one.
	    ka.xs.resize(ka.xs.size() - d.second.GetGrid().GetJointGrid().InterDegree());

	    // Q2-space (sub)grid.
	    std::vector<double> q2;
	    for (int iq = ti; iq < tind[i+1]; iq++)
	      q2.push_back(qg[iq] * qg[iq]);
	    ka.q2s = q2;

	    // Get sizes of both x and q2 grids.
	    const int xsize  = ka.xs.size();
	    const int q2size = ka.q2s.size();

	    // Now fill in vector of distributions.
	    std::vector<double> xf(xsize * q2size);
	    for (int iq = ti; iq < tind[i+1]; iq++)
	      {
		const std::vector<double> xfx = QCDEvToPhys(xfg[iq].GetObjects()).at(d.first).GetDistributionJointGrid();
		for (int ix = 0; ix < (int) ka.xs.size(); ix++)
		  xf[ix * q2size + iq - ti] = xfx[ix];
	      }
	    ka.xfs = xf;

	    LHKnotArrayNF.insert({d.first, ka});
	  }
	_KnotArray.insert({qg[ti] * qg[ti], LHKnotArrayNF});
      }

    // Write grid if required.
    if (_WriteGrid)
      WriteGrid();
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

    // Check flavour scheme.
    if (_setup.FlavourScheme != EvolutionSetup::VFNS &&
	_setup.FlavourScheme != EvolutionSetup::FFNS)
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "Unknown flavour scheme.") << std::endl;
	passed = false;
      }

    // Check number of flavours in the FFNS (if relevant).
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

    // Check virtuality.
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

    // Check distribution evolution.
    if (_setup.PDFEvolution != EvolutionSetup::exactmu &&
	_setup.PDFEvolution != EvolutionSetup::exactalpha &&
	_setup.PDFEvolution != EvolutionSetup::expandalpha &&
	_setup.PDFEvolution != EvolutionSetup::truncated)
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "Unknown PDF evolution.") << std::endl;
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

    // Check that heavy-quark thresholds and masses are sorted.
    if (!std::is_sorted(_setup.Thresholds.begin(), _setup.Thresholds.end()))
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "The heavy-quark thresholds are not sorted.") << std::endl;
	passed = false;
      }

    if (!std::is_sorted(_setup.Masses.begin(), _setup.Masses.end()))
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "The heavy-quark masses are not sorted.") << std::endl;
	passed = false;
      }

    // Check maximum number of active flavours.
    if (_setup.Thresholds.size() < 3 || _setup.Thresholds.size() > 6)
      {
	std::cout << error("InitialiseEvolution::CheckSetup", "Maximum number of active flavours out of range.") << std::endl;
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

    // String with the report of parameters.
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
    report += " with maximum " + std::to_string(_setup.Thresholds.size()) + " active flavours\n";

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
    report += " with maximum " + std::to_string(_setup.Thresholds.size()) + " active flavours\n";

    // Ren. / fact. scales ratio.
    report += "- muR / mF: " + std::to_string(_setup.RenFacRatio) + "\n";

    // Heavy-quark masses and thresholds.
    std::vector<std::string> Thq{"mud", "muu", "mus", "muc", "mub", "mut"};
    if (_setup.MassRenScheme == EvolutionSetup::POLE)
      {
	std::vector<std::string> Mq{"Md", "Mu", "Ms", "Mc", "Mb", "Mt"};
	report += "- Pole heavy-quark masses:\n";
	for (int i = 0; i < (int) _setup.Masses.size(); i++)
	  if (_setup.Masses[i] > 0)
	    report += "  + " + Mq[i] + " = " + std::to_string(_setup.Masses[i]) + " GeV\n";
	report += "- Pole heavy-quark threholds:\n";
	for (int i = 0; i < (int) _setup.Thresholds.size(); i++)
	  if (_setup.Thresholds[i] > 0)
	    report += "  + " + Thq[i] + " = " + std::to_string(_setup.Thresholds[i]) + " GeV\n";
      }
    if (_setup.MassRenScheme == EvolutionSetup::MSBAR)
      {
	std::vector<std::string> Mq{"md(md)", "mu(mu)", "ms(ms)", "mc(mc)", "mb(mb)", "mt(mt)"};
	report += "- MSbar heavy-quark masses:\n";
	for (int i = 0; i < (int) _setup.Masses.size(); i++)
	  if (_setup.Masses[i] > 0)
	    report += "  + " + Mq[i] + " = " + std::to_string(_setup.Masses[i]) + " GeV\n";
	report += "- Pole heavy-quark threholds:\n";
	for (int i = 0; i < (int) _setup.Thresholds.size(); i++)
	  if (_setup.Thresholds[i] > 0)
	    report += "  + " + Thq[i] + " = " + std::to_string(_setup.Thresholds[i]) + " GeV\n";
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

  //_________________________________________________________________________________
  void InitialiseEvolution::WriteGridInfo()
  {
    // Create folder for the info file and the grid(s).
    if (mkdir(_setup.name.c_str(), 0777) != 0)
      throw std::runtime_error(error("InitialiseEvolution::WriteGridInfo", "Cannot create folder for the LHAPDF set."));

    // Information string.
    std::string info;

    // Write information.
    info += "SetDesc: 'Set generated with APFEL++. Name: " + _setup.name + "'\n";
    info += "Authors: V. Bertone, S. Carrazza\n";
    info += "Reference: arXiv:1708.00911\n";
    info += "Format: lhagrid1\n";
    info += "DataVersion: 1\n";
    info += "NumMembers: 0\n";
    info += "Particle: 2212\n";

    info += "Flavors: [";
    for (int i = - (int) _setup.Thresholds.size(); i <= (int) _setup.Thresholds.size(); i++)
      info += std::to_string((i == 0 ? 21 : i)) + ", ";
    info = info.substr(0, info.size() - 2);
    info += "]\n";

    info += "OrderQCD: " + std::to_string(_setup.PerturbativeOrder) + "\n";
    info += "FlavorScheme: variable\n";
    info += "NumFlavors: " + std::to_string(_setup.Thresholds.size()) + "\n";

    double xmin = 1;
    for (auto const& gp : _setup.GridParameters)
      xmin = std::min(xmin, gp.xmin);
    info += "XMin: " + std::to_string(xmin) + "\n";
    info += "XMax: 1\n";

    info += "QMin: " + std::to_string(_setup.Qmin) + "\n";
    info += "QMax: " + std::to_string(_setup.Qmax) + "\n";

    info += "MZ: " + std::to_string(_setup.QQCDRef) + "\n";

    const std::vector<std::string> ms{"MUp: ", "MDown: ", "MStrange: ", "MCharm: ", "MBottom: ", "MTop: "};
    for (int i = 0; i < (int) _setup.Masses.size(); i++)
      info += ms[i] + std::to_string(_setup.Masses[i]) + "\n";

    const std::vector<std::string> th{"ThresholdUp: ", "ThresholdDown: ", "ThresholdStrange: ", "ThresholdCharm: ", "ThresholdBottom: ", "ThresholdTop: "};
    for (int i = 0; i < (int) _setup.Thresholds.size(); i++)
      info += th[i] + std::to_string(_setup.Thresholds[i]) + "\n";

    info += "AlphaS_MZ: " + std::to_string(_setup.AlphaQCDRef) + "\n";
    info += "AlphaS_OrderQCD: " + std::to_string(_setup.PerturbativeOrder) + "\n";
    info += "AlphaS_Type: ipol\n";

    const QGrid<double> qg{_setup.nQg, _setup.Qmin, _setup.Qmax, _setup.InterDegreeQ, _setup.Thresholds};
    info += "AlphaS_Qs: [";
    for (auto const& q : qg.GetQGrid())
      info += std::to_string(q) + ", ";
    info = info.substr(0, info.size() - 2);
    info += "]\n";

    info += "AlphaS_Vals: [";
    for (auto const& q : qg.GetQGrid())
      info += std::to_string(_as(q)) + ", ";
    info = info.substr(0, info.size() - 2);
    info += "]\n";

    // Open info file and print the information in it.
    std::ofstream out(_setup.name + "/" + _setup.name + ".info");
    out << info;
    out.close();
  }

  //_________________________________________________________________________________
  void InitialiseEvolution::WriteGrid()
  {
    // Check if the destination folder exists.
    struct stat info;
    if (stat(_setup.name.c_str(), &info) != 0)
      throw std::runtime_error(error("InitialiseEvolution::WriteGrid", "Folder for the LHAPDF set does not exist."));

    // Determine the first available member index and constuct the
    // open file name.
    std::ofstream out;
    std::string filename = _setup.name + "/" + _setup.name + "_0000";
    const int maxmem = 10000;
    for (int i = 1; i < maxmem; i++)
      {
	if (!(bool) std::ifstream((filename + ".dat").c_str()))
	  {
	    out.open(filename + ".dat");
	    break;
	  }
	filename = filename.substr(0, filename.size() - (int) log10(i) - 1) + std::to_string(i);
      }

    out << std::scientific;
    out << "PdfType: replica\n";
    out << "Format: lhagrid1\n";

    for(auto const& ka : _KnotArray)
      {
	out << "---\n";

	// Write x grid.
	for (auto const& x : ka.second.begin()->second.xs)
	  out << x << " ";
	out << "\n";

	// Write Q grid.
	for (auto const& q2 : ka.second.begin()->second.q2s)
	  out << sqrt(q2) << " ";
	out << "\n";

	// Write flavour indices.
	for (int i = - (int) _setup.Thresholds.size(); i <= (int) _setup.Thresholds.size(); i++)
	  out << (i == 0 ? 21 : i) << " ";
	out << "\n";

	// Array of distributions.
	const int nd = 2 * _setup.Thresholds.size() + 1;
	const int np = ka.second.begin()->second.xs.size() * ka.second.begin()->second.q2s.size();
	double dist[np][nd];

	// Gather tabulated distributions.
	int id = 0;
	for (auto const& d : ka.second)
	  if (std::abs(d.first) <= (int) _setup.Thresholds.size())
	    {
	      const std::vector<double> xf = d.second.xfs;
	      for (int jp = 0; jp < (int) xf.size(); jp++)
		dist[jp][id] = (std::abs(xf[jp]) > eps15 ? xf[jp] : 0);
	      id++;
	    }

	// Print distributions.
	for (int jp = 0; jp < np; jp++)
	  {
	    for (int jd = 0; jd < nd; jd++)
	      out << dist[jp][jd] << "  ";
	    out << "\n";
	  }
      }

    out << "---\n";
    out.close();
  }
}
