//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/tmdcrosssections.h"
#include "apfel/tabulateobject.h"
#include "apfel/rotations.h"
#include "apfel/ogataquadrature.h"
#include "apfel/integrator.h"

using namespace std;

namespace apfel {
  //_____________________________________________________________________________
  function<double(double const&)> TmdCrossSectionDY(double                                                                   const& Vs,
						    double                                                                   const& Q,
						    double                                                                   const& y,
						    function<Set<Distribution>(double const&, double const&, double const&)> const& EvolvedTMDPDFs,
						    function<double(double const&)>                                          const& Alphas,
						    function<vector<double>(double const&)>                                  const& fEWCharges,
						    int                                                                      const& PerturbativeOrder,
						    vector<double>                                                           const& Thresholds,
						    double                                                                   const& cmuf,
						    double                                                                   const& czetaf)
  {
    // Compute "muf" and "zetaf". They are assumed to be proportional
    // to Q and Q^2, respectively, through "cmuf" and "czetaf".
    const double muf   = cmuf * Q;
    const double zetaf = czetaf * Q * Q;

    // Compute values of x1 and x2.
    const double x1 = Q * exp(-y) / Vs;
    const double x2 = Q * exp(y) / Vs;

    // Get EW charges.
    const vector<double> Bq = fEWCharges(Q);

    // Tabulate input TMDs in the impact parameter to make the
    // integral faster.
    const auto TabFunc    = [] (double const& b) -> double{ return log(b); };
    const auto InvTabFunc = [] (double const& fb) -> double{ return exp(fb); };
    const TabulateObject<Set<Distribution>> TabulatedTMDs{[&] (double const& b) -> Set<Distribution>
	{ return EvolvedTMDPDFs(b, muf, zetaf); }, 50, 5e-5, 30, 3, {}, TabFunc, InvTabFunc};

    // Construct the TMD luminosisty in b scale to be fed to be
    // trasformed in qT space.
    const auto TMDLumib = [=] (double const& b) -> double
      {
	// Get map of the TMDs in "x1" and "x2" and rotate them into the
	// physical basis.
	map<int,double> TMD1 = QCDEvToPhys(TabulatedTMDs.EvaluateMapxQ(x1, b));
	map<int,double> TMD2 = QCDEvToPhys(TabulatedTMDs.EvaluateMapxQ(x2, b));

	// Construct the combination of TMDs weighted by the EW
	// charges.
	double lumi = 0;
	for (int i = 1; i <= 6; i++)
	  lumi += Bq[i-1] * ( TMD1[i] * TMD2[-i] + TMD1[-i] * TMD2[i] );

	// Divide by x1 and x2 because the "EvolvedTMDPDFs" function
	// returns x times the TMDs.
	lumi /= x1 * x2;

	return lumi;
      };

    // Compute hard coefficient, including the other kinematic
    // factors.
    const double hcs = ConvFact * FourPi * HardCrossSectionDY(PerturbativeOrder, Alphas(muf), NF(muf, Thresholds), muf/Q) / 9 / pow(Vs*Q,2);

    return [TMDLumib,hcs,Q] (double const& qT) -> double
      {
	const double cutPref = 1 + pow(qT/Q,2) / 2;
	return hcs * cutPref * OgataQuadrature(TMDLumib, qT);
	// For test purposes, one can comment out the line above, in
	// which the Ogata quadrature is used, and comment in the two
	// lines below to use the DGauss quadrature. The second will
	// be substantially slower but should return consistent
	// results.
	//const Integrator integrand{[=] (double const& bT) -> double{ return TMDLumib(bT) * j0(qT * bT) * bT / 2; }};
	//return hcs * cutPref * integrand.integrate(0.0000001,50,eps9);
      };
  }

  //_____________________________________________________________________________
  function<double(double const&)> TmdCrossSectionDY(double                                                                const& Vs,
						    double                                                                const& Qmin,
						    double                                                                const& Qmax,
						    double                                                                const& ymin,
						    double                                                                const& ymax,
						    function<Set<Distribution>(double const&)>                            const& InTMDPDFs,
						    function<vector<double>(double const&, double const&, double const&)> const& EvolFact,
						    function<double(double const&)>                                       const& Alphas,
						    function<vector<double>(double const&)>                               const& fEWCharges,
						    int                                                                   const& PerturbativeOrder,
						    vector<double>                                                        const& Thresholds,
						    double                                                                const& cmuf,
						    double                                                                const& czetaf,
						    double                                                                const& IntEps)
  {
    // Tabulation function and its inverse.
    const auto TabFunc    = [] (double const& b) -> double{ return log(b); };
    const auto InvTabFunc = [] (double const& fb) -> double{ return exp(fb); };

    // Tabulate initial TMD PDFs.
    const TabulateObject<Set<Distribution>> TabInTMDs{[&] (double const& b) -> Set<Distribution>
	{ return InTMDPDFs(b); }, 50, 5e-5, 30, 3, {}, TabFunc, InvTabFunc};

    // If integration over Q and y is required, it is advantageous to
    // perform these integrations before the integral over the impact
    // parameter b (Fourier transform).
    return [=] (double const& qT) -> double
      {
	const Integrator Qintegrand{
	  [&] (double const& lnQ) -> double
	    {
	      const double Q = exp(lnQ);
	      // Compute "muf" and "zetaf". They are assumed to be
	      // proportional to Q and Q^2, respectively, through
	      // "cmuf" and "czetaf".
	      const double muf   = cmuf * Q;
	      const double zetaf = czetaf * Q * Q;

	      // Get EW charges.
	      const vector<double> Bq = fEWCharges(Q);

	      // Tabulate input TMDs in the impact parameter to make the
	      // integral faster.
	      const TabulateObject<Set<Distribution>> TabEvTMDs{[&] (double const& b) -> Set<Distribution>
		  { return EvolFact(b, muf, zetaf) * TabInTMDs.Evaluate(b); }, 50, 5e-5, 30, 3, {}, TabFunc, InvTabFunc};

	      const auto bintegrand = [&] (double const& b) -> double
		{
		  // Define y integrand.
		  const Integrator yintegrand{
		    [&] (double const& ey) -> double
		    {
		      // Compute values of x1 and x2.
		      const double x1 = Q / ey / Vs;
		      const double x2 = Q * ey / Vs;

		      // Get map of the TMDs in "x1" and "x2" and
		      // rotate them into the physical basis.
		      map<int,double> TMD1 = QCDEvToPhys(TabEvTMDs.EvaluateMapxQ(x1, b));
		      map<int,double> TMD2 = QCDEvToPhys(TabEvTMDs.EvaluateMapxQ(x2, b));

		      // Construct the combination of TMDs weighted by
		      // the EW charges.
		      double lumi = 0;
		      for (int i = 1; i <= 6; i++)
			lumi += Bq[i-1] * ( TMD1[i] * TMD2[-i] + TMD1[-i] * TMD2[i] );

		      // Divide by x1 and x2 because "EvolvedTMDPDFs"
		      // returns x times the TMD PDFs.
		      lumi /= x1 * x2 * ey ;

		      return lumi;
		    }
		  };

		  // Integrate in y over [ymin:ymax] or return the
		  // value if ymin = ymax.
		  if (ymin == ymax)
		    return yintegrand.integrand(exp(ymin));
		  else
		    return yintegrand.integrate(exp(ymin), exp(ymax), IntEps);
		};

	      const double cutPref  = 1 + pow(qT/Q,2) / 2;
	      const double HardFact = HardCrossSectionDY(PerturbativeOrder, Alphas(muf), NF(muf, Thresholds), cmuf);
	      return cutPref * HardFact * OgataQuadrature(bintegrand, qT) / Q;
	    }
	};
	// Integrate in Q over [Qmin:Qmax] or return the
	// value if Qmin = Qmax.
	if (Qmin == Qmax)
	  return ConvFact * FourPi * Qintegrand.integrand(log(Qmin)) / 9 / Vs / Vs;
	else
	  return ConvFact * FourPi * Qintegrand.integrate(log(Qmin), log(Qmax), IntEps) / 9 / Vs / Vs;
      };
  }

  //_____________________________________________________________________________
  double HardCrossSectionDY(int const& pt, double const& Alphas, int const& nf, double const& kappa)
    {
      // Compute log anf its powers.
      const double lQ2  = 2 * log(kappa);
      const double lQ22 = lQ2 * lQ2;
      const double lQ23 = lQ22 * lQ2;
      const double lQ24 = lQ23 * lQ2;

      // Compute coupling and its powers.
      const double as  = Alphas / FourPi;
      const double as2 = as * as;

      // Now compute hard cross section according to the perturbative
      // order.
      double hxs = 1;
      if (pt > 0)
	hxs += 2 * as * CF * ( - lQ22 - 3 * lQ2 - 8 + 7 * Pi2 / 6 );
      if (pt > 1)
	hxs += 2 * as2 * CF *
	  ( CF * ( lQ24 + 6 * lQ23 + ( 25. - 7 * Pi2 / 3 ) * lQ22 + ( 93. / 2. - 5 * Pi2 - 24 * zeta3 ) * lQ2
		   + 511. / 8. - 83 * Pi2 / 6 - 30 * zeta3 + 67 * Pi2 * Pi2 / 60 ) +
	    CA * ( - 11 * lQ23 / 9 + ( - 233. / 18. + Pi2 / 3 ) * lQ22 + ( - 2545. / 54. + 22 * Pi2 / 9 + 26 * zeta3 ) * lQ2
		   - 51157. / 648. + 1061 * Pi2 / 108 + 313 * zeta3 / 9 - 4 * Pi2 * Pi2 / 45 ) +
	    TR * nf * ( 4 * lQ23 / 9 + 38 * lQ22 / 9 + ( 418. / 27. - 8 * Pi2 / 9 ) * lQ2
			+ 4085. / 162. - 91 * Pi2 / 27 + 4 * zeta3 / 9 ) );

      return hxs;
    };
}
