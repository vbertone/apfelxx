//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>
#include <apfel/massivecoefficientfunctionsunp_sl.h>

#include <iomanip>

// LH Toy PDFs
double xglu(double const& x)
{
  return 1.7 * pow(x, -0.1) * pow(1 - x, 5);
}

int main()
{
  apfel::Timer t;

  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{100, 1e-5, 3}, apfel::SubGrid{100, 1e-1, 3}, apfel::SubGrid{100, 6e-1, 3}, apfel::SubGrid{80, 8.5e-1, 5}}};

  // Kinematics
  const double Q   = 2;
  const double mh  = sqrt(2);
  const double lxi = 2 * log(Q / mh);
  const double eta = 1 / ( 1 + pow(2 * mh / Q, 2) );

  // Integration accuracy
  const double IntEps = apfel::eps5;

  // Massive coeffient functions
  const apfel::Operator Om21g  {g, apfel::Cm21gNC{eta},     IntEps};
  const apfel::Operator Om22ns {g, apfel::Cm22nsNC{eta},    IntEps};
  const apfel::Operator Om22gc {g, apfel::Cm22gNC{eta},     IntEps};
  const apfel::Operator Om22gb {g, apfel::Cm22bargNC{eta},  IntEps};
  const apfel::Operator Om22psc{g, apfel::Cm22psNC{eta},    IntEps};
  const apfel::Operator Om22psb{g, apfel::Cm22barpsNC{eta}, IntEps};
  const apfel::Operator Om22ps = Om22psc + lxi * Om22psb;
  const apfel::Operator Om22g  = Om22gc  + lxi * Om22gb;

  // Construct extended operators
  const apfel::Operator Om21gExt {Om21g, eta};
  const apfel::Operator Om22nsExt{Om22ns, eta};
  const apfel::Operator Om22psExt{Om22ps, eta};
  const apfel::Operator Om22gExt {Om22g, eta};
  t.stop();

  // Define test distribution
  const apfel::Distribution TestDist{g, xglu};

  // Multiply the operators by the test distribution
  t.start();
  const apfel::Distribution Dm21g  = Om21g  * TestDist;
  const apfel::Distribution Dm22ns = Om22ns * TestDist;
  const apfel::Distribution Dm22ps = Om22ps * TestDist;
  const apfel::Distribution Dm22g  = Om22g  * TestDist;
  t.stop();

  t.start();
  const apfel::Distribution Dm21gExt  = Om21gExt  * TestDist;
  const apfel::Distribution Dm22nsExt = Om22nsExt * TestDist;
  const apfel::Distribution Dm22psExt = Om22psExt * TestDist;
  const apfel::Distribution Dm22gExt  = Om22gExt  * TestDist;
  t.stop();

  // Print results
  std::cout << std::scientific << std::endl;
  std::cout << "Q = " << Q << " GeV" << std::endl;
  std::cout << "mh = " << mh << " GeV" << std::endl;
  std::cout << "\nUnextended vs. exdented representation of massive operators" << std::endl;
  std::cout << "    x    \t"
            << "Cm21g    \t"
            << "Cm22ns   \t"
            << "Cm22ps   \t"
            << "Cm22g    \t"
            << std::endl;
  const int nx = 100;
  const double xmin = 1e-5;
  const double xmax = 0.999;
  const double xstp = exp( log(xmax / xmin) / ( nx - 1 ) );
  for (double x = xmin; x < xmax * 1.000001; x *= xstp)
    std::cout << std::setprecision(4) << x << "\t"
              << ( Dm21gExt.Evaluate(x)  - Dm21g.Evaluate(x / eta) )  / ( Dm21g.Evaluate(x / eta)  + 1 )<< "\t"
              << ( Dm22nsExt.Evaluate(x) - Dm22ns.Evaluate(x / eta) ) / ( Dm22ns.Evaluate(x / eta) + 1 ) << "\t"
              << ( Dm22psExt.Evaluate(x) - Dm22ps.Evaluate(x / eta) ) / ( Dm22ps.Evaluate(x / eta) + 1 ) << "\t"
              << ( Dm22gExt.Evaluate(x)  - Dm22g.Evaluate(x / eta) )  / ( Dm22g.Evaluate(x / eta)  + 1 ) << "\t"
              << std::endl;
  std::cout << std::endl;

  return 0;
}
