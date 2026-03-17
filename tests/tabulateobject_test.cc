//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

int main()
{
  // Mute APFEL++
  apfel::SetVerbosityLevel(0);

  // Vectors of masses and thresholds
  const std::vector<double> Thresholds = {0, 0, 0, sqrt(2), 4.5, 175};

  // Perturbative order
  const int PerturbativeOrder = 2;

  // Function to be tabulated in asref
  const std::function<apfel::TabulateObject<double>(double const&)> AlphasTT = [=] (double const& asref) -> apfel::TabulateObject<double>
  {
    // Tabulate the running coupling in mu
    apfel::AlphaQCD as{asref, 91.2, Thresholds, PerturbativeOrder};
    return apfel::TabulateObject<double>{as, 100, 0.9, 1001, 3};
  };

  // Tabulate in asref
  const apfel::TabulateObject AlphaTab{AlphasTT, 20, 0.100, 0.125, 3, Thresholds, [] (double const& x) -> double { return x; }, [] (double const& y) -> double { return y; }};

  // Print result of the double tabulation
  const int nmu = 50;
  const double mumin = 1;
  const double mumax = 91.2;
  const double mustp = exp(log(mumax / mumin) / ( nmu - 1 ));
  const std::vector<double> asrefv{0.115, 0.116, 0.117, 0.118, 0.119, 0.120};

  std::cout << std::scientific;
  std::cout << "\n     mu        "
            << " as(MZ)=0.115   "
            << " as(MZ)=0.116   "
            << " as(MZ)=0.117   "
            << " as(MZ)=0.118   "
            << " as(MZ)=0.119   "
            << " as(MZ)=0.120   "
            << std::endl;
  for (double mu = mumin; mu <= 1.000001 * mumax; mu *= mustp)
    {
      std::cout << mu << "\t";
      for (double asref : asrefv)
        std::cout << AlphaTab.Evaluate(asref).Evaluate(mu) << "\t";
      std::cout << std::endl;
    }

  return 0;
}
