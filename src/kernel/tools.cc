//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/tools.h"
#include "apfel/constants.h"
#include "apfel/distribution.h"
#include "apfel/set.h"
#include "apfel/doubleobject.h"
#include "apfel/alphaqcd.h"
#include "apfel/alphaqed.h"

#include <algorithm>
#include <math.h>
#include <numeric>

namespace apfel
{
  //_________________________________________________________________________
  int NF(double const& Q, std::vector<double> const& Thresholds)
  {
    // Compute number of active flavours the the PDF initial scale
    int nf = 0;
    for (auto const& v : Thresholds)
      if (Q > v)
        nf++;
      else
        break;
    return nf;
  }

  //_________________________________________________________________________
  double DeltaFun(double const& a, double const& b, double const& c)
  {
    return sqrt( a * a + b * b + c * c - 2 * ( a * b + b * c + c * a ) );
  }

  //_________________________________________________________________________
  std::vector<double> ElectroWeakCharges(double const& Q, bool const& virt, int const& sel)
  {
    // Relevant constants
    const double Q2    = Q * Q;
    const double MZ2   = apfel::ZMass * apfel::ZMass;
    const double GmZ2  = apfel::GammaZ * apfel::GammaZ;
    const double S2ThW = apfel::Sin2ThetaW;
    const double VD    = - 0.5 + 2 * S2ThW / 3;
    const double VU    = + 0.5 - 4 * S2ThW / 3;
    const double AD    = - 0.5;
    const double AU    = + 0.5;
    const double Ve    = - 0.5 + 2 * S2ThW;
    const double Ae    = - 0.5;
    const std::vector<double> Vq = {VD, VU, VD, VU, VD, VU};
    const std::vector<double> Aq = {AD, AU, AD, AU, AD, AU};

    // Propagator and its square
    double PZ;
    double PZ2;
    if (virt)
      {
        PZ  = Q2 * ( Q2 -  MZ2 ) / ( pow(Q2 - MZ2,2) + MZ2 * GmZ2 ) / ( 4 * S2ThW * ( 1 - S2ThW ) );
        PZ2 = pow(Q2,2) / ( pow(Q2 - MZ2,2) + MZ2 * GmZ2 ) / pow(4 * S2ThW * ( 1 - S2ThW ),2);
      }
    else
      {
        PZ  = Q2 / ( Q2 + MZ2 ) / ( 4 * S2ThW * ( 1 - S2ThW ) );
        PZ2 = PZ * PZ;
      }

    // Build electroweak charges
    std::vector<double> Charges(6, 0.);
    if (sel < 0 || sel > 5)
      for (auto i = 0; i < 6; i++)
        Charges[i] = apfel::QCh2[i]
                     - 2 * apfel::QCh[i] * Vq[i] * Ve * PZ
                     + ( Ve * Ve + Ae * Ae ) * ( Vq[i] * Vq[i] + Aq[i] * Aq[i] ) * PZ2;
    else
      Charges[sel] = apfel::QCh2[sel]
                     - 2 * apfel::QCh[sel] * Vq[sel] * Ve * PZ
                     + ( Ve * Ve + Ae * Ae ) * ( Vq[sel] * Vq[sel] + Aq[sel] * Aq[sel] ) * PZ2;

    return Charges;
  }

  //_________________________________________________________________________
  std::vector<double> ElectroWeakChargesNWA()
  {
    // Relevant constants
    const double S2ThW = apfel::Sin2ThetaW;
    const double VD    = - 0.5 + 2 * S2ThW / 3;
    const double VU    = + 0.5 - 4 * S2ThW / 3;
    const double AD    = - 0.5;
    const double AU    = + 0.5;
    const double Ve    = - 0.5 + 2 * S2ThW;
    const double Ae    = - 0.5;
    const std::vector<double> Vq = {VD, VU, VD, VU, VD, VU};
    const std::vector<double> Aq = {AD, AU, AD, AU, AD, AU};

    // Propagator and its square
    double PZ2 = pow(apfel::ZMass, 4) * M_PI / apfel::ZMass / apfel::ZMass / apfel::GammaZ / pow(4 * S2ThW * ( 1 - S2ThW ), 2) / 2;

    std::vector<double> Charges(6, 0.);
    for (auto i = 0; i < 6; i++)
      Charges[i] = ( Ve * Ve + Ae * Ae ) * ( Vq[i] * Vq[i] + Aq[i] * Aq[i] ) * PZ2;

    return Charges;
  }

  //_________________________________________________________________________________
  std::vector<double> ConcatenateAndSortVectors(std::vector<double> const& v1, std::vector<double> const& v2)
  {
    std::vector<double> v12 = v1;
    v12.insert(v12.end(), v2.begin(), v2.end());
    std::sort(v12.begin(), v12.end());
    return v12;
  }

  //_________________________________________________________________________________
  template<>
  double dabs(double const& d)
  {
    return std::abs(d);
  }

  //_________________________________________________________________________________
  template<>
  double dabs(Distribution const& d)
  {
    // The absolute value of a distribution is assumed to be equal to
    // its mean value over the definition interval, that is:
    //
    // 1/(b-a)\int_a^b dx f(x)
    //
    // where the integral is approximated through the rectangle rule.

    // Get joint grid and respective values
    const std::vector<double>& jg = d.GetGrid().GetJointGrid().GetGrid();
    const std::vector<double>& jv = d.GetDistributionJointGrid();

    // Now compute the integral
    double integ = 0;
    for (int i = 0; i < (int) jg.size() - 1; i++)
      integ += jv[i] * ( jg[i+1] - jg[i] );

    // Divide the integral by the interval and return the absolute
    // value of the result.
    return std::abs(integ / ( jg.back() - jg.front() ));
  }

  //_________________________________________________________________________________
  template<>
  double dabs(Set<Distribution> const& d)
  {
    // For a set of distributions, the absolute value is assumed to be
    // the average absolute value over the distributions of the set.
    double lgt = 0;
    const std::map<int, Distribution> objs = d.GetObjects();
    for (auto const& e : objs)
      lgt += dabs(e.second);
    return lgt / objs.size();
  }

  //_________________________________________________________________________________
  template<>
  double dabs(DoubleObject<Distribution> const& d)
  {
    // For a double object of distributions, the absolute value is
    // assumed to be the sum over all the terms, weighted by the
    // coefficients, of the product of the two distributions in each
    // term.
    double av = 0;
    for (auto const& e : d.GetTerms())
      av += e.coefficient * dabs(e.object1) * dabs(e.object2);

    return std::abs(av);
  }

  //_________________________________________________________________________________
  std::vector<double> ProductExpansion(std::vector<double> const& r)
  {
    const int k = r.size();
    std::vector<double> a(k + 1, 1.);
    std::vector<double> p(k + 1, 1.);
    for (int n = 0; n <= k; n++)
      {
        p[n] = a[0];
        std::vector<double> f(k + 1, 0.);
        for (int j = 0; j <= k; j++)
          for (int i = j + 1; i <= k; i++)
            f[j] += r[i-1] * a[i];
        a = f;
      }
    return p;
  }

  //_________________________________________________________________________________
  int factorial(int const& n)
  {
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
  }

  //_________________________________________________________________________________
  double GetSIATotalCrossSection(int                 const& PerturbativeOrder,
                                 double              const& Q,
                                 double              const& AlphaQCD,
                                 double              const& AlphaQED,
                                 std::vector<double> const& Thresholds,
                                 QuarkFlavour        const& Comp)
  {
    // Get time-like electroweak charges
    const std::vector<double> Bq = apfel::ElectroWeakCharges(Q, true, Comp);

    // Effective number of flavours
    double nf = NF(Q, Thresholds);

    // Sum the charges
    const double sumq = std::accumulate(Bq.begin(), Bq.begin() + nf, 0.);

    // Total born coss setion
    double sigma0tot = FourPi * pow(AlphaQED, 2) * NC * sumq / 3 / Q / Q;

    // QCD coupling
    const double as1 = AlphaQCD / FourPi;
    const double as2 = as1 * as1;

    // QCD correction factor
    double Ree = 1;
    if(PerturbativeOrder >= 1)
      Ree += as1 * CF * 3;
    if(PerturbativeOrder >= 2)
      Ree += as2 * ( CF * CF * ( - 3. / 2 )
                     + CA * CF * ( - 44 * zeta3 + 123. / 2 )
                     + nf * CF * TR * ( 16 * zeta3 - 22. ) );

    return 1e-3 * ConvFact * sigma0tot * Ree;
  }
}
