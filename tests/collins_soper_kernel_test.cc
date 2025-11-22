//
// APFEL++ 2017
//
// Author: Chiara Bissolotti: cbissolotti@anl.gov
//

#include <apfel/apfelxx.h>

#include <iomanip>

int main()
{
  // Initializing grid
  const apfel::Grid g{{apfel::SubGrid{100, 1e-5, 3}, apfel::SubGrid{60, 1e-1, 3}, apfel::SubGrid{50, 6e-1, 3}}};

  // Vector of thresholds
  const std::vector<double> thresholds = {0, 0, 0, sqrt(2), 4.5, 175};

  // Initializing TMD objects
  const auto tmd_objects = apfel::InitializeTmdObjects(g, thresholds);

  // Strong coupling
  const double alphas_ref = 0.35;
  const double mu_ref = 2;
  const apfel::AlphaQCD alphas_obj{alphas_ref, mu_ref, thresholds, 3};
  const auto alphas = [&] (double const& mu) -> double{ return alphas_obj.Evaluate(mu); };

  // Computing Collins-Soper kernel
  const auto cs_kernel_LL    = apfel::CollinsSoperKernel(tmd_objects, alphas, 0);
  const auto cs_kernel_NLL   = apfel::CollinsSoperKernel(tmd_objects, alphas, 1);
  const auto cs_kernel_NNLL  = apfel::CollinsSoperKernel(tmd_objects, alphas, 2);
  const auto cs_kernel_NNNLL = apfel::CollinsSoperKernel(tmd_objects, alphas, 3);

  // Focus on the physically relevant range
  const std::vector<double> b_values = {0.2, 0.5, 0.8, 1.0};
  const std::vector<double> mu_values = {2.0, 5.0, 10.0, 91.0};

  std::cout << "\nCollins-Soper Kernel K(b, mu) values:" << std::endl;
  std::cout << "b [GeV^-1] |  mu [GeV]  |   K(b,mu) LL |  K(b,mu) NLL | K(b,mu) NNLL | K(b,mu) N3LL |" << std::endl;
  std::cout << "-------------------------------------------------------------------------------------" << std::endl;
  const auto is_valid_result = [] (double const& x) -> bool{ return !std::isnan(x) && !std::isinf(x) && std::abs(x) < 100.0; };
  for (const auto& b : b_values)
    {
      for (const auto& mu : mu_values)
        {
          double k_LL    = cs_kernel_LL(b, mu);
          double k_NLL   = cs_kernel_NLL(b, mu);
          double k_NNLL  = cs_kernel_NNLL(b, mu);
          double k_NNNLL = cs_kernel_NNNLL(b, mu);
          if (is_valid_result(k_LL) && is_valid_result(k_NLL) && is_valid_result(k_NNLL) && is_valid_result(k_NNNLL))
            {
              std::cout << std::fixed << std::setprecision(6)
                        << std::setw(10) << b << " | "
                        << std::setw(10) << mu << " | " << std::setprecision(8)
                        << std::setw(12) << k_LL << " | "
                        << std::setw(12) << k_NLL << " | "
                        << std::setw(12) << k_NNLL << " | "
                        << std::setw(12) << k_NNNLL << " | "
                        << std::endl;
            }
          else
            {
              std::cout << std::fixed << std::setprecision(6)
                        << std::setw(10) << b << " | "
                        << std::setw(10) << mu << " | " << std::setprecision(8)
                        << std::setw(12) << "N/A" << " | "
                        << std::setw(12) << "N/A" << " | "
                        << std::setw(12) << "N/A" << " | "
                        << " (numerical instability)" << std::endl;
            }
        }
      std::cout << "--------------------------------------------------------------------------------------" << std::endl;
    }

  const double fixed_mu = 10.0;
  std::cout << "\nCollins-Soper kernel as a function of b at mu = " << fixed_mu << " GeV:" << std::endl;
  std::cout << "b (GeV^-1) |  K(b,mu) NLL |" << std::endl;
  std::cout << "---------------------------" << std::endl;
  for (double b = 0.1; b <= 1.2; b += 0.1)
    {
      const double k = cs_kernel_NLL(b, fixed_mu);
      if (is_valid_result(k))
        {
          std::cout << std::fixed << std::setprecision(6)
                    << std::setw(10) << b << " | " << std::setprecision(8)
                    << std::setw(12) << k << " | "
                    << std::endl;
        }
      else
        {
          std::cout << std::fixed << std::setprecision(6)
                    << std::setw(10) << b << " | "
                    << std::setw(12) << "N/A" << " | "
                    << " (numerical instability)" << std::endl;
        }
    }
  std::cout << "\n";

  return 0;
}
