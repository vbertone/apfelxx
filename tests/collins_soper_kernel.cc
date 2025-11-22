#include "apfel/apfelxx.h"
#include "apfel/grid.h"
#include "apfel/tmdbuilder.h"
#include "apfel/alphaqcd.h"
#include "apfel/constants.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits>

int main()
{
  try {
    std::cout << "Initializing grid..." << std::endl;
    const std::vector<apfel::SubGrid> subgrid = {
      apfel::SubGrid{100, 1e-5, 3},
      apfel::SubGrid{60, 1e-1, 3},
      apfel::SubGrid{50, 6e-1, 3}
    };
    const apfel::Grid g{subgrid};

    const std::vector<double> thresholds = {0, 0, 0, sqrt(2), 4.5, 175};

    std::cout << "Initializing TMD objects..." << std::endl;
    const auto tmd_objects = apfel::InitializeTmdObjects(g, thresholds);
    
    const double alphas_ref = 0.35;
    const double mu_ref = 2.0;
    const apfel::AlphaQCD alphas_obj{alphas_ref, mu_ref, thresholds, 4};
    const auto alphas = [&] (double const& mu) -> double { 
      if (mu <= 0) throw std::invalid_argument("Scale mu must be positive");
      return alphas_obj.Evaluate(mu); 
    };
    
    std::cout << "Computing Collins-Soper kernel..." << std::endl;
    
    const auto cs_kernel_LL = apfel::CollinsSoperKernel(tmd_objects, alphas, 0);
    const auto cs_kernel_NLL = apfel::CollinsSoperKernel(tmd_objects, alphas, 1);
    const auto cs_kernel_NNLL = apfel::CollinsSoperKernel(tmd_objects, alphas, 2);
    
    // Focus on the physically relevant range
    const std::vector<double> b_values = {0.2, 0.5, 0.8, 1.0};
    const std::vector<double> mu_values = {2.0, 5.0, 10.0, 91.0};
    
    std::cout << "\nCollins-Soper Kernel K(b, mu) values:" << std::endl;
    std::cout << "b (GeV^-1) | mu (GeV)  | K(b,mu) LL  | K(b,mu) NLL | K(b,mu) NNLL" << std::endl;
    std::cout << "-------------------------------------------------------------------" << std::endl;
    
    auto is_valid_result = [](double x) {
        return !std::isnan(x) && !std::isinf(x) && 
               std::abs(x) < 100.0; // Physical values should be moderate
    };
    
    for (const auto& b : b_values) {
      if (b <= 0) {
        std::cerr << "Warning: Impact parameter b must be positive. Skipping b = " << b << std::endl;
        continue;
      }
      
      for (const auto& mu : mu_values) {
        if (mu <= 0) {
          std::cerr << "Warning: Scale mu must be positive. Skipping mu = " << mu << std::endl;
          continue;
        }
        
        try {
          double k_LL = cs_kernel_LL(b, mu);
          double k_NLL = cs_kernel_NLL(b, mu);
          double k_NNLL = cs_kernel_NNLL(b, mu);
          
          if (is_valid_result(k_LL) && is_valid_result(k_NLL) && is_valid_result(k_NNLL)) {
            std::cout << std::fixed << std::setprecision(2)
                      << std::setw(10) << b << " | "
                      << std::setw(10) << mu << " | "
                      << std::setw(12) << k_LL << " | "
                      << std::setw(12) << k_NLL << " | "
                      << std::setw(12) << k_NNLL
                      << std::endl;
          } else {
            std::cout << std::fixed << std::setprecision(2)
                      << std::setw(10) << b << " | "
                      << std::setw(10) << mu << " | "
                      << std::setw(12) << "N/A" << " | "
                      << std::setw(12) << "N/A" << " | "
                      << std::setw(12) << "N/A"
                      << " (numerical instability)" << std::endl;
          }
        } catch (const std::exception& e) {
          std::cout << std::fixed << std::setprecision(2)
                    << std::setw(10) << b << " | "
                    << std::setw(10) << mu << " | "
                    << std::setw(12) << "N/A" << " | "
                    << std::setw(12) << "N/A" << " | "
                    << std::setw(12) << "N/A"
                    << " (calculation failed)" << std::endl;
        }
      }
      std::cout << "-------------------------------------------------------------------" << std::endl;
    }
    
    const double fixed_mu = 10.0;
    std::cout << "\nCollins-Soper kernel as a function of b at mu = " << fixed_mu << " GeV:" << std::endl;
    std::cout << "b (GeV^-1) | K(b,mu) NLL" << std::endl;
    std::cout << "------------------------" << std::endl;
    
    for (double b = 0.1; b <= 1.2; b += 0.1) {
      try {
        double k = cs_kernel_NLL(b, fixed_mu);
        if (is_valid_result(k)) {
          std::cout << std::fixed << std::setprecision(2)
                    << std::setw(10) << b << " | "
                    << std::setw(12) << k
                    << std::endl;
        } else {
          std::cout << std::fixed << std::setprecision(2)
                    << std::setw(10) << b << " | "
                    << std::setw(12) << "N/A"
                    << " (numerical instability)" << std::endl;
        }
      } catch (const std::exception& e) {
        std::cout << std::fixed << std::setprecision(2)
                  << std::setw(10) << b << " | "
                  << std::setw(12) << "N/A"
                  << " (calculation failed)" << std::endl;
      }
    }
    
    std::cout << "\nComputation completed successfully!" << std::endl;
    
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  
  return 0;
}
