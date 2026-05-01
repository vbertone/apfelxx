//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>
#include <apfel/splittingfunctionsunp_sl.h>

int main()
{
  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{500, 1e-11, 3}, apfel::SubGrid{400, 1e-1, 3}, apfel::SubGrid{400, 6e-1, 3}, apfel::SubGrid{400, 9e-1, 5}}};

  // Lower integration limit
  const double xmin = 1e-11;

  // Number of active flavours
  const int nf = 5;
  
  // Get N3LO splitting functions
  const apfel::Operator aO3nsp{g, apfel::aP3nsp{nf}, apfel::eps5};
  const apfel::Operator aO3nsm{g, apfel::aP3nsm{nf}, apfel::eps5};
  const apfel::Operator aO3nss{g, apfel::aP3nss{nf}, apfel::eps5};
  const apfel::Operator xO3nsp{g, apfel::P3nsp{nf}, apfel::eps5};
  const apfel::Operator xO3nsm{g, apfel::P3nsm{nf}, apfel::eps5};
  const apfel::Operator xO3nss{g, apfel::P3nss{nf}, apfel::eps5};

  // Functions to convolute with splitting and matching functions to
  // obtain the appropriate integrals
  std::cout << std::scientific;
  for (int N = 1; N <= 16; N++)
    {
      const apfel::Distribution Dist{g, [&] (double const& x) -> double { return pow(xmin / x, N - 1); }};
      std::cout << "N = " << N << "\t"
		<< "aPNSP = " << (aO3nsp * Dist).Evaluate(xmin) << "\t"
		<< "xPNSP = " << (xO3nsp * Dist).Evaluate(xmin) << "\t|\t"
		<< "aPNSM = " << (aO3nsm * Dist).Evaluate(xmin) << "\t"
		<< "xPNSM = " << (xO3nsm * Dist).Evaluate(xmin) << "\t|\t"
		<< "aPNSV = " << (aO3nss * Dist).Evaluate(xmin) << "\t"
		<< "xPNSV = " << (xO3nss * Dist).Evaluate(xmin) << "\t"
		<< std::endl;
    }

  const apfel::aP3nsp aE3nsp{nf};
  const apfel::aP3nsm aE3nsm{nf};
  const apfel::aP3nss aE3nss{nf};
  const apfel::P3nsp  xE3nsp{nf};
  const apfel::P3nsm  xE3nsm{nf};
  const apfel::P3nss  xE3nss{nf};
  const std::vector<double> xval = {0.0001, 0.0005, 0.001, 0.005, 0.01, 0.03, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.99};
  std::cout << "\n";
  for (double x : xval)
    std::cout << x << "\t\t"
	      << ( aE3nsm.Regular(x) + aE3nsm.Singular(x) ) / ( xE3nsm.Regular(x) + xE3nsm.Singular(x) ) << "\t\t"
	      << ( aE3nsp.Regular(x) + aE3nsp.Singular(x) ) / ( xE3nsp.Regular(x) + xE3nsp.Singular(x) ) << "\t\t"
	      << aE3nss.Regular(x) / xE3nss.Regular(x) << "\t"
	      << std::endl;
    
  return 0;
}
