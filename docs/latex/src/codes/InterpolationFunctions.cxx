/*
  Author: Valerio Bertone
 */

// APFEL++ libs
#include <apfel/apfelxx.h>

int main() {

  // Initialize space- and time-like splitting functions.
  const apfel::Grid g{{{100, 1e-3, 0}}};
  const std::function<double(double const&)> f = [] (double const& x) -> double { return 2 + cos(M_PI * log(x)); };
  const apfel::Distribution d{g, f};
  const std::vector<double> xg = g.GetSubGrid(0).GetGrid();
  const int nnodes = 5;
/*
  for (int i = 0; i < nnodes + 3; i++)
    std::cout << xg[i] << "\t";
  std::cout << std::endl;
  exit(-1);
*/
  int nx = 10000;
  double xmin = 1e-3;
  double xmax = 1.413e-3;
  double xstp = exp( log( xmax / xmin ) / ( nx - 1 ) );
/*
  for (double x = xmin; x <= xmax; x *= xstp)
    {
      std::cout << std::scientific << x << "\t";
      for (int i = 0; i < nnodes; i++)
	std::cout << d.InterpolantLog(i, log(x), g.GetSubGrid(0)) << "\t";
      std::cout << std::endl;
    }
*/
  nx = 10000;
  xmin = 1e-3;
  xmax = 1;
  xstp = exp( log( xmax / xmin ) / ( nx - 1 ) );
  const apfel::Grid g0{{{100, 1e-3, 0}}};
  const apfel::Grid g1{{{100, 1e-3, 1}}};
  const apfel::Grid g2{{{100, 1e-3, 2}}};
  const apfel::Grid g3{{{100, 1e-3, 3}}};
  const apfel::Grid g4{{{100, 1e-3, 4}}};
  const apfel::Grid g5{{{100, 1e-3, 5}}};
  const apfel::Distribution d0{g0, f};
  const apfel::Distribution d1{g1, f};
  const apfel::Distribution d2{g2, f};
  const apfel::Distribution d3{g3, f};
  const apfel::Distribution d4{g4, f};
  const apfel::Distribution d5{g5, f};
  for (double x = xmin; x <= xmax; x *= xstp)
    {
      std::cout << std::scientific << x << "\t"
		<< f(x) << "\t"
		<< d1.Evaluate(x) << "\t"
		<< d2.Evaluate(x) << "\t"
		<< d3.Evaluate(x) << "\t"
		<< d4.Evaluate(x) << "\t"
		<< d5.Evaluate(x) << "\t"
		<< d0.Evaluate(x) << "\t"
		<< std::endl;
    }

/*
  // Derivative of the distribution
  const auto df = [&] (double const& x) -> double{ return - M_PI * sin(M_PI * log(x)) / x; };

  // Integral of the distribution
  const auto inf = [&] (double const& x) -> double{ return x * ( M_PI * sin(M_PI * log(x)) + cos(M_PI * log(x)) ) / ( 1 + M_PI * M_PI) + 2 * x; };

  // Grid
  const apfel::Grid gg{{{100, 9.9e-6, 5}, {100, 1e-1, 5}, {40, 8e-1, 5}}};

  // Interpolated distribution
  const apfel::Distribution xgluon{gg, f};
  nx = 10000;
  xmin = 1e-3;
  xmax = 1;
  xstp = exp( log( xmax / xmin ) / ( nx - 1 ) );
  for (double x = xmin; x <= xmax; x *= xstp)
    {
      const double original = f(x);
      const double derorig  = df(x);
      const double intorig  = inf(0.65) - inf(x);
      const double interpol = xgluon.Evaluate(x);
      const double derive   = xgluon.Derive(x);
      const double integr   = xgluon.Integrate(x, 0.65);
      std::cout << std::scientific << x << "\t"
                << original << "\t"
                << interpol << "\t"
                << interpol / original << "\t"
                << derorig  << "\t"
                << derive   << "\t"
                << derive / derorig << "\t"
                << intorig  << "\t"
                << integr   << "\t"
                << integr / intorig << "\t"
                << std::endl;
    }
*/
  return 0;
}
