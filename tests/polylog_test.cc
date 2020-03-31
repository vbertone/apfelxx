//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>
#include <apfel/hpolyweights.h>

int main()
{
  const double xx = 0.4;

  double *x = new double;
  int *nw = new int;
  int *wn1 = new int;
  int *wn2 = new int;
  *x = xx;
  *nw = 5;
  *wn1 = -1;
  *wn2 = 1;

  // HPLOG all at once up to weight 4
  double Hr1[3], Hr2[9], Hr3[27], Hr4[81], Hr5[243];
  apfel::hplog_(x, nw, Hr1, Hr2, Hr3, Hr4, Hr5, wn1, wn2);

  std::cout << " x = " << xx << std::endl;
  std::cout << "                        "
            << "   HPOLY        "
            << "   HPLOG_       "
            << std::endl;
  std::cout << std::scientific;
  std::cout << "HPL({-1}, x)          = " << apfel::hpoly({-1}, 0.4)          << ",\t" << Hr1[apfel::HPLogMap({-1})]          << std::endl;
  std::cout << "HPL({0,1}, x)         = " << apfel::hpoly({0,1}, 0.4)         << ",\t" << Hr2[apfel::HPLogMap({0,1})]         << std::endl;
  std::cout << "HPL({-1,0,0}, x)      = " << apfel::hpoly({-1,0,0}, 0.4)      << ",\t" << Hr3[apfel::HPLogMap({-1,0,0})]      << std::endl;
  std::cout << "HPL({-1,-1,1,0}, x)   = " << apfel::hpoly({-1,-1,1,0}, 0.4)   << ",\t" << Hr4[apfel::HPLogMap({-1,-1,1,0})]   << std::endl;
  std::cout << "HPL({-1,-1,1,0,1}, x) = " << apfel::hpoly({-1,-1,1,0,1}, 0.4) << ",\t" << Hr5[apfel::HPLogMap({-1,-1,1,0,1})] << std::endl;
  std::cout << "\n";

  // Performance test
  const double nc = 1000;
  apfel::Timer t;
  for (int i = 0; i < nc; i++)
    apfel::hpoly({-1,-1,1,0,1}, 0.4);
  std::cout << "Calling HPOLY  " << nc << " times... ";
  t.stop();

  t.start();
  for (int i = 0; i < nc; i++)
    apfel::hplog_(x, nw, Hr1, Hr2, Hr3, Hr4, Hr5, wn1, wn2);
  std::cout << "Calling HPLOG_ " << nc << " times... ";
  t.stop();
  std::cout << "\n";

  return 0;
}
