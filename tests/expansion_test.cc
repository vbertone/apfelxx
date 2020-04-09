//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

#include <iomanip>

int main()
{
  std::cout << std::setprecision(8) << std::scientific;

  // Vector of zero's
  const std::vector<double> zeros{0.1, 4, 1, 4, 12, 69, 76, 0.1};

  // Factorised function
  const std::function<double(double const&)> factF = [=] (double const& x) -> double
  {
    double res = 1;
    for (auto const& z : zeros)
      res *= x - z;
    return res;
  };

  // Expanded function
  const std::function<double(double const&)> expF = [=] (double const& x) -> double
  {
    const int k = zeros.size();
    std::vector<double> Pk(k + 1, 1.);
    double res = 0;
    for (int gamma = 0; gamma <= k; gamma++)
      {
        res += pow(-1, gamma) * Pk[0] * pow(x, k - gamma);
        Pk = apfel::VectorComposition(zeros, Pk);
      }
    return res;
  };

  const double x = 2;
  std::cout << "\nPolynomial of degree " << zeros.size() << " computed at x = " << x << std::endl;
  std::cout << "Zero's in: [  ";
  for (auto e : zeros)
    std::cout << e << "\t";
  std::cout << "]\n";
  std::cout << "\n   Factorised              Expanded" << std::endl;
  std::cout << factF(x) << "\t\t" << expF(x) << std::endl;
  std::cout << std::endl;

  return 0;
}
