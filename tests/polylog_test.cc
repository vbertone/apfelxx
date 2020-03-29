//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

int main()
{
  const double xx = 0.4;
/*
  double *x = new double;
  int *nw = new int;
  int *wn1 = new int;
  int *wn2 = new int;
  *x = xx;
  *nw = 4;
  *wn1 = -1;
  *wn2 = 1;
*/
  // HPOLY all at once
  std::map<int, std::vector<double>> hpls = apfel::hpoly(xx);

  // HPLOG all at once
  double Hr1[3], Hr2[9], Hr3[27], Hr4[81];
  //apfel::hplog_(x, nw, Hr1, Hr2, Hr3, Hr4, wn1, wn2);

  std::cout << " x = " << xx << std::endl;
  std::cout << "                        "
            << "  HPOLY [1]     "
            << "  HPOLY [2]     "
            << "   HPLOG        "
            << std::endl;
  std::cout << std::scientific;
  std::cout << "HPL({-1}, x)          = "
            << apfel::hpoly(std::vector<int> {-1}, 0.4) << ",\t"
            << hpls.at(apfel::WeightAndIndex({-1}).first)[apfel::WeightAndIndex({-1}).second]
            << ",\t" << Hr1[0] << std::endl;
  std::cout << "HPL({0,1}, x)         = "
            << apfel::hpoly({0, 1}, 0.4) << ",\t"
            << hpls.at(apfel::WeightAndIndex({0,1}).first)[apfel::WeightAndIndex({0,1}).second]
            << ",\t" << Hr2[7] << std::endl;
  std::cout << "HPL({-1,0,0}, x)      = "
            << apfel::hpoly({-1,0,0}, 0.4) << ",\t"
            << hpls.at(apfel::WeightAndIndex({-1,0,0}).first)[apfel::WeightAndIndex({-1,0,0}).second]
            << ",\t" << Hr3[12] << std::endl;
  std::cout << "HPL({-1,-1,1,0}, x)   = "
            << apfel::hpoly({-1,-1,1,0}, 0.4) << ",\t"
            << hpls.at(apfel::WeightAndIndex({-1,-1,1,0}).first)[apfel::WeightAndIndex({-1,-1,1,0}).second]
            << ",\t" << Hr4[45] << std::endl;
  std::cout << "HPL({-1,-1,1,0,1}, x) = "
            << apfel::hpoly({-1,-1,1,0,1}, 0.4) << ",\t"
            << hpls.at(apfel::WeightAndIndex({-1,-1,1,0,1}).first)[apfel::WeightAndIndex({-1,-1,1,0,1}).second]
            << ",\t    ---" << std::endl;
  std::cout << "\n";

  // Performance test
  const double nc = 1000;
  apfel::Timer t;
  for (int i = 0; i < nc; i++)
    apfel::hpoly({-1,-1,1,0,1}, 0.4);
  std::cout << "Calling HPOLY [1] " << nc << " times... ";
  t.stop();

  t.start();
  for (int i = 0; i < nc; i++)
    apfel::hpoly(xx);
  std::cout << "Calling HPOLY [2] " << nc << " times... ";
  t.stop();
/*
  t.start();
  for (int i = 0; i < nc; i++)
    apfel::hplog_(x, nw, Hr1, Hr2, Hr3, Hr4, wn1, wn2);
  std::cout << "Calling HPLOG     " << nc << " times... ";
  t.stop();
  std::cout << "\n";

  delete x;
  delete nw;
  delete wn1;
  delete wn2;
*/
  return 0;
}
