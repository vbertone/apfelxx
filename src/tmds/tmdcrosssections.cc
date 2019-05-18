//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/tmdcrosssections.h"
#include "apfel/hardfactors.h"
#include "apfel/tabulateobject.h"
#include "apfel/rotations.h"
#include "apfel/ogataquadrature.h"
#include "apfel/integrator.h"
#include "apfel/constants.h"
#include "apfel/tools.h"

namespace apfel
{
  //_____________________________________________________________________________
  std::function<double(double const&)> TmdCrossSectionDY(double                                                                        const& Vs,
                                                         double                                                                        const& Q,
                                                         double                                                                        const& y,
                                                         std::function<Set<Distribution>(double const&, double const&, double const&)> const& EvolvedTMDPDFs,
                                                         std::function<double(double const&)>                                          const& Alphas,
                                                         std::function<std::vector<double>(double const&)>                             const& fEWCharges,
                                                         int                                                                           const& PerturbativeOrder,
                                                         std::vector<double>                                                           const& Thresholds,
                                                         double                                                                        const& cmuf,
                                                         double                                                                        const& czetaf)
  {
    // Compute "muf" and "zetaf". They are assumed to be proportional
    // to Q and Q^2, respectively, through "cmuf" and "czetaf".
    const double muf   = cmuf * Q;
    const double zetaf = czetaf * Q * Q;

    // Compute values of x1 and x2.
    const double x1 = Q * exp(-y) / Vs;
    const double x2 = Q * exp(y) / Vs;

    // Get EW charges.
    const std::vector<double> Bq = fEWCharges(Q);

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
      std::map<int,double> TMD1 = QCDEvToPhys(TabulatedTMDs.EvaluateMapxQ(x1, b));
      std::map<int,double> TMD2 = QCDEvToPhys(TabulatedTMDs.EvaluateMapxQ(x2, b));

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
    const double hcs = ConvFact * FourPi * HardFactorDY(PerturbativeOrder, Alphas(muf), NF(muf, Thresholds), muf/Q) / 9 / pow(Vs*Q,2);

    return [TMDLumib,hcs,Q] (double const& qT) -> double
    {
      const double cutPref = 1 + pow(qT/Q,2) / 2;
      return 0;//hcs * cutPref * OgataQuadrature(TMDLumib, qT);
      // For test purposes, one can comment out the line above, in
      // which the Ogata quadrature is used, and comment in the two
      // lines below to use the DGauss quadrature. The second will
      // be substantially slower but should return consistent
      // results.
      //const Integrator integrand{[=] (double const& bT) -> double{ return TMDLumib(bT) * j0(qT * bT) * bT / 2; }};
      //return hcs * cutPref * integrand.integrate(0.00005, 30, eps9);
    };
  }

  //_____________________________________________________________________________
  std::function<double(double const&)> TmdCrossSectionSIDIS(double                                                                        const& Vs,
                                                            double                                                                        const& x,
                                                            double                                                                        const& y,
                                                            double                                                                        const& z,
                                                            std::function<Set<Distribution>(double const&, double const&, double const&)> const& EvolvedTMDPDFs,
                                                            std::function<Set<Distribution>(double const&, double const&, double const&)> const& EvolvedTMDFFs,
                                                            std::function<double(double const&)>                                          const& Alphas,
                                                            std::function<std::vector<double>(double const&)>                             const& fEWCharges,
                                                            int                                                                           const& PerturbativeOrder,
                                                            std::vector<double>                                                           const& Thresholds,
                                                            double                                                                        const& cmuf,
                                                            double                                                                        const& czetaf)
  {
    // Compute virtuality of the photon.
    const double Q = sqrt( x * y ) * Vs;

    // Compute "muf" and "zetaf". They are assumed to be proportional
    // to Q and Q^2, respectively, through "cmuf" and "czetaf".
    const double muf   = cmuf * Q;
    const double zetaf = czetaf * Q * Q;

    // Get EW charges.
    const std::vector<double> Bq = fEWCharges(Q);

    // Tabulate input TMD PDFs and FFs in the impact parameter to make
    // the integral faster.
    const auto TabFunc    = [] (double const& b) -> double{ return log(b); };
    const auto InvTabFunc = [] (double const& fb) -> double{ return exp(fb); };
    const TabulateObject<Set<Distribution>> TabulatedPDFs{[&] (double const& b) -> Set<Distribution>
      { return EvolvedTMDPDFs(b, muf, zetaf); }, 50, 5e-5, 30, 3, {}, TabFunc, InvTabFunc};
    const TabulateObject<Set<Distribution>> TabulatedFFs{[&] (double const& b) -> Set<Distribution>
      { return EvolvedTMDFFs(b, muf, zetaf); }, 50, 1e-2, 30, 3, {}, TabFunc, InvTabFunc};

    // Construct the TMD luminosisty in b scale to be fed to be
    // trasformed in qT space.
    const auto TMDLumib = [=] (double const& b) -> double
    {
      // Get map of the TMDs in "x1" and "x2" and rotate them into the
      // physical basis.
      std::map<int,double> PDFs = QCDEvToPhys(TabulatedPDFs.EvaluateMapxQ(x, b));
      std::map<int,double> FFs  = QCDEvToPhys(TabulatedFFs.EvaluateMapxQ(z, b));

      // Construct the combination of TMDs weighted by the EW
      // charges.
      double lumi = 0;
      for (int i = 1; i <= 6; i++)
        lumi += Bq[i-1] * ( PDFs[i] * FFs[i] + PDFs[-i] * FFs[-i] );

      // Divide by x1 and x2 because the "EvolvedTMDPDFs" function
      // returns x times the TMDs.
      lumi /= x * z;

      return lumi;
    };

    // Compute hard coefficient, including the other kinematic
    // factors.
    const double hcs = 2 * M_PI * HardFactorSIDIS(PerturbativeOrder, Alphas(muf), NF(muf, Thresholds), muf/Q) / x / y / pow(Q,2);
    //const double hcs = ConvFact * 2 * M_PI * HardFactorSIDIS(PerturbativeOrder, Alphas(muf), NF(muf, Thresholds), muf/Q) / x / y / pow(Q,2);

    return [TMDLumib,hcs] (double const& qT) -> double
    {
      return 0;//hcs * OgataQuadrature(TMDLumib, qT);
      // For test purposes, one can comment out the line above, in
      // which the Ogata quadrature is used, and comment in the two
      // lines below to use the DGauss quadrature. The second will
      // be substantially slower but should return consistent
      // results.
      //const Integrator integrand{[=] (double const& bT) -> double{ return TMDLumib(bT) * j0(qT * bT) * bT / 2; }};
      //return hcs * cutPref * integrand.integrate(0.00005,30,eps9);
    };
  }

  //_____________________________________________________________________________
  std::function<double(double const&)> TmdCrossSectionDY(double                                                                          const& Vs,
                                                         double                                                                          const& Qmin,
                                                         double                                                                          const& Qmax,
                                                         double                                                                          const& ymin,
                                                         double                                                                          const& ymax,
                                                         std::function<Set<Distribution>(double const&)>                                 const& InTMDPDFs,
                                                         std::function<std::vector<double>(double const&, double const&, double const&)> const& EvolFact,
                                                         std::function<double(double const&)>                                            const& Alphas,
                                                         std::function<std::vector<double>(double const&)>                               const& fEWCharges,
                                                         int                                                                             const& PerturbativeOrder,
                                                         std::vector<double>                                                             const& Thresholds,
                                                         double                                                                          const& cmuf,
                                                         double                                                                          const& czetaf,
                                                         double                                                                          const& IntEps)
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
          const std::vector<double> Bq = fEWCharges(Q);

          // Tabulate input TMDs in the impact parameter to make the
          // integral faster.
          const TabulateObject<Set<Distribution>> TabEvTMDs{
            [&] (double const& b) -> Set<Distribution>
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
                std::map<int,double> TMD1 = QCDEvToPhys(TabEvTMDs.EvaluateMapxQ(x1, b));
                std::map<int,double> TMD2 = QCDEvToPhys(TabEvTMDs.EvaluateMapxQ(x2, b));

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
              return yintegrand.integrand(exp(ymin)) * exp(ymin);
            else
              return yintegrand.integrate(exp(ymin), exp(ymax), IntEps);
          };

          const double cutPref  = 1 + pow(qT/Q,2) / 2;
          const double HardFact = HardFactorDY(PerturbativeOrder, Alphas(muf), NF(muf, Thresholds), cmuf);
          return 0;//2 * cutPref * HardFact * OgataQuadrature(bintegrand, qT);
        }
      };
      // Integrate in Q over [Qmin:Qmax] or return the
      // value if Qmin = Qmax.
      if (Qmin == Qmax)
        return ConvFact * FourPi * Qintegrand.integrand(log(Qmin)) / 9 / Vs / Vs / Qmin;
      else
        return ConvFact * FourPi * Qintegrand.integrate(log(Qmin), log(Qmax), IntEps) / 9 / Vs / Vs;
    };
  }
}
