#include "apfel/apfelxx.h"
#include "LHAPDF/LHAPDF.h"

int main()
{

  //constants
  double mc(1.3),mb(4.5),mt(175);
  double IntEps = 1e-5;
  int nQ = 100;
  double Qmin = 2;
  double Qmax = 225;
  int intdeg = 3;

  double x(0.1),Q(10);

  //define PDF function
  LHAPDF::PDF* dist = LHAPDF::mkPDF("CT18NNLO");
  auto my_PDF = [=] (double const& x, double const& Q) -> std::map<int,double> {return dist->xfxQ(x,Q);};

  //define alphaQCD
  const auto strongCoupling = [=] (double const& Q) -> double {return dist->alphasQ(Q);};

  // set perturbative order
  const int pto = 2;
  // define a x-grid through a two subgrids that have 100/50 interpolation nodes,
  // x_min=1e-5/1e-1 and are of degree 3/3
  const apfel::Grid g{{apfel::SubGrid{100,1e-5,3},apfel::SubGrid{50,1e-1,3}}};
  // define a vector of mass thresholds -> needed to not interpolate across a discontinuity
  const std::vector<double> Thresholds = {0,0,0,mc,mb,mt};
  // define a function that returns the effective electro-weak charges as a function of Q
  const auto fEW = [=] (double const& Q) -> std::vector<double>
  {
    return apfel::ElectroWeakCharges(Q,false);
  }; // use predefined function within APFEL++

  // rotate the PDFs into the QCD evolution basis
  const auto PDFrotated = [&] (double const& x, double const& Q) -> std::map<int,double>
  {
    return apfel::PhysToQCDEv(my_PDF(x,Q));
  };

  // the scaling variable chi(n)
  double n = 1;

  // Construct the operator grids for the F2 structure function
  const auto F2objects = apfel::InitializeF2NCObjectsASACOT(
                           g,Thresholds,IntEps,nQ,Qmin,Qmax,intdeg,n);

  // Build the structure function
  const std::map<int,apfel::Observable<>> F2 = apfel::BuildStructureFunctions(
                                                 F2objects,PDFrotated,pto,strongCoupling,fEW);

  // Evaluate the total F2 at x=0.1 and Q=10
  std::cout<<std::fixed<<std::setprecision(10);
  std::cout<<"\n\n";
  std::cout<<"Direct evaluation"<<std::endl;
  std::cout<<"F2(x=0.1,Q=10) = "<<F2.at(0).Evaluate(x,Q)<<std::endl;
  std::cout<<"\n\n";

  // Interpolate the Q-dependence of the total F2
  const apfel::TabulateObject<apfel::Distribution> F2total
  {
    [&] (double const& Q) -> apfel::Distribution{return F2.at(0).Evaluate(Q);},
    nQ,Qmin,Qmax,intdeg,Thresholds};

  // Evaluate the Q-interpolated total F2 at x=0.1 and Q=10
  std::cout<<"\n\n";
  std::cout<<"Tabulated evaluation"<<std::endl;
  std::cout<<"F2(x=0.1,Q=10) = "<<F2total.EvaluatexQ(x=0.1,Q=10)<<std::endl;
  std::cout<<"\n\n";

  return 0;
}