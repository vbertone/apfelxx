#include <fstream>
#include <iomanip>

#include "apfel/apfelxx.h"
#include "apfel/structurefunctionbuilder_ACOT.h"
#include "LHAPDF/LHAPDF.h"

const std::vector<double> x
{
  1.00000000e-06, 1.14975700e-06, 1.32194115e-06, 1.51991108e-06,
  1.74752840e-06, 2.00923300e-06, 2.31012970e-06, 2.65608778e-06,
  3.05385551e-06, 3.51119173e-06, 4.03701726e-06, 4.64158883e-06,
  5.33669923e-06, 6.13590727e-06, 7.05480231e-06, 8.11130831e-06,
  9.32603347e-06, 1.07226722e-05, 1.23284674e-05, 1.41747416e-05,
  1.62975083e-05, 1.87381742e-05, 2.15443469e-05, 2.47707636e-05,
  2.84803587e-05, 3.27454916e-05, 3.76493581e-05, 4.32876128e-05,
  4.97702356e-05, 5.72236766e-05, 6.57933225e-05, 7.56463328e-05,
  8.69749003e-05, 1.00000000e-04, 1.14975700e-04, 1.32194115e-04,
  1.51991108e-04, 1.74752840e-04, 2.00923300e-04, 2.31012970e-04,
  2.65608778e-04, 3.05385551e-04, 3.51119173e-04, 4.03701726e-04,
  4.64158883e-04, 5.33669923e-04, 6.13590727e-04, 7.05480231e-04,
  8.11130831e-04, 9.32603347e-04, 1.07226722e-03, 1.23284674e-03,
  1.41747416e-03, 1.62975083e-03, 1.87381742e-03, 2.15443469e-03,
  2.47707636e-03, 2.84803587e-03, 3.27454916e-03, 3.76493581e-03,
  4.32876128e-03, 4.97702356e-03, 5.72236766e-03, 6.57933225e-03,
  7.56463328e-03, 8.69749003e-03, 1.00000000e-02, 1.14975700e-02,
  1.32194115e-02, 1.51991108e-02, 1.74752840e-02, 2.00923300e-02,
  2.31012970e-02, 2.65608778e-02, 3.05385551e-02, 3.51119173e-02,
  4.03701726e-02, 4.64158883e-02, 5.33669923e-02, 6.13590727e-02,
  7.05480231e-02, 8.11130831e-02, 9.32603347e-02, 1.07226722e-01,
  1.23284674e-01, 1.41747416e-01, 1.62975083e-01, 1.87381742e-01,
  2.15443469e-01, 2.47707636e-01, 2.84803587e-01, 3.27454916e-01,
  3.76493581e-01, 4.32876128e-01, 4.97702356e-01, 5.72236766e-01,
  6.57933225e-01, 7.56463328e-01, 8.69749003e-01
};
const std::vector<double> Q{1.4,5,10,100};

int main()
{

  // path for output
  std::string path = "/home/peter/work/codes/apfelxx_ACOT/tests/ACOT_tests/APFELxx_results";

  // prepare PDFs and alphas
  LHAPDF::PDF* pdf = LHAPDF::mkPDF("CT18NLO_1_1",0);
  const auto PDFrotated = [&] (double const& x, double const& Q) -> std::map<int,double> {return apfel::PhysToQCDEv(pdf->xfxQ(x,Q));};
  const auto alphas = [&] (double const& Q) -> double{return pdf->alphasQ(Q);};
  LHAPDF::PDF* pdf_NNLO = LHAPDF::mkPDF("CT18NLO_1_1",0);
  const auto PDFrotated_NNLO = [&] (double const& x, double const& Q) -> std::map<int,double> {return apfel::PhysToQCDEv(pdf_NNLO->xfxQ(x,Q));};
  const auto alphas_NNLO = [&] (double const& Q) -> double{return pdf_NNLO->alphasQ(Q);};
  const std::vector<double> Thresholds = {0,0,0,pdf_NNLO->quarkMass(4),pdf_NNLO->quarkMass(5),pdf_NNLO->quarkMass(6)};

  // grid specifications
  double IntEps = 1e-5;
  int nQ = 60;
  double Qmin = 1.29;
  double Qmax = 101;
  int intdeg = 3;
  apfel::Grid g{{{25,1e-6,4},{20,1e-2,4},{10,1e-1,4},{5,5e-1,4}}};

  /////////////////////////////
  //////////// NLO ////////////
  /////////////////////////////
  int pto = 1;

  const auto fEW = [=] (double const& Q) -> std::vector<double> {return apfel::ElectroWeakCharges(Q,false);};
  const auto fPVEW = [=] (double const& Q) -> std::vector<double> {return apfel::ParityViolatingElectroWeakCharges(Q,false);};
  const auto fCKM2 = [=] (double const& ) -> std::vector<double> {return apfel::CKM2;};

  // Neutral current
  const auto F2objects = apfel::InitializeF2NCObjectsSACOT(g,Thresholds,IntEps,nQ,Qmin,Qmax,intdeg);
  const std::map<int,apfel::Observable<>> F2 = apfel::BuildStructureFunctions(F2objects,PDFrotated,pto,alphas,fEW);
  const apfel::TabulateObject<apfel::Distribution> F2total {[&] (double const& Q) -> apfel::Distribution{return F2.at(0).Evaluate(Q);},nQ,Qmin,Qmax,intdeg,Thresholds};

  const auto FLobjects = apfel::InitializeFLNCObjectsSACOT(g,Thresholds,IntEps,nQ,Qmin,Qmax,intdeg);
  const std::map<int,apfel::Observable<>> FL = apfel::BuildStructureFunctions(FLobjects,PDFrotated,pto,alphas,fEW);
  const apfel::TabulateObject<apfel::Distribution> FLtotal {[&] (double const& Q) -> apfel::Distribution{return FL.at(0).Evaluate(Q);},nQ,Qmin,Qmax,intdeg,Thresholds};

  const auto F3objects = apfel::InitializeF3NCObjectsSACOT(g,Thresholds,IntEps,nQ,Qmin,Qmax,intdeg);
  const std::map<int,apfel::Observable<>> F3 = apfel::BuildStructureFunctions(F3objects,PDFrotated,pto,alphas,fPVEW);
  const apfel::TabulateObject<apfel::Distribution> F3total {[&] (double const& Q) -> apfel::Distribution{return F3.at(0).Evaluate(Q);},nQ,Qmin,Qmax,intdeg,Thresholds};

  // charged current
  const auto F2plusObjects = apfel::InitializeF2CCPlusObjectsSACOT(g,Thresholds,IntEps,nQ,Qmin,Qmax,intdeg);
  const std::map<int,apfel::Observable<>> F2plus = apfel::BuildStructureFunctions(F2plusObjects,PDFrotated,pto,alphas,fCKM2);
  const auto F2minusObjects = apfel::InitializeF2CCMinusObjectsSACOT(g,Thresholds,IntEps,nQ,Qmin,Qmax,intdeg);
  const std::map<int,apfel::Observable<>> F2minus = apfel::BuildStructureFunctions(F2minusObjects,PDFrotated,pto,alphas,fCKM2);
  const apfel::TabulateObject<apfel::Distribution> F2WP {[&] (double const& Q) -> apfel::Distribution{return 2*(F2plus.at(0).Evaluate(Q)+F2minus.at(0).Evaluate(Q));},nQ,Qmin,Qmax,intdeg,Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F2WM {[&] (double const& Q) -> apfel::Distribution{return 2*(F2plus.at(0).Evaluate(Q)-F2minus.at(0).Evaluate(Q));},nQ,Qmin,Qmax,intdeg,Thresholds};

  const auto FLplusObjects = apfel::InitializeFLCCPlusObjectsSACOT(g,Thresholds,IntEps,nQ,Qmin,Qmax,intdeg);
  const std::map<int,apfel::Observable<>> FLplus = apfel::BuildStructureFunctions(FLplusObjects,PDFrotated,pto,alphas,fCKM2);
  const auto FLminusObjects = apfel::InitializeFLCCMinusObjectsSACOT(g,Thresholds,IntEps,nQ,Qmin,Qmax,intdeg);
  const std::map<int,apfel::Observable<>> FLminus = apfel::BuildStructureFunctions(FLminusObjects,PDFrotated,pto,alphas,fCKM2);
  const apfel::TabulateObject<apfel::Distribution> FLWP {[&] (double const& Q) -> apfel::Distribution{return 2*(FLplus.at(0).Evaluate(Q)+FLminus.at(0).Evaluate(Q));},nQ,Qmin,Qmax,intdeg,Thresholds};
  const apfel::TabulateObject<apfel::Distribution> FLWM {[&] (double const& Q) -> apfel::Distribution{return 2*(FLplus.at(0).Evaluate(Q)-FLminus.at(0).Evaluate(Q));},nQ,Qmin,Qmax,intdeg,Thresholds};

  const auto F3plusObjects = apfel::InitializeF3CCPlusObjectsSACOT(g,Thresholds,IntEps,nQ,Qmin,Qmax,intdeg);
  const std::map<int,apfel::Observable<>> F3plus = apfel::BuildStructureFunctions(F3plusObjects,PDFrotated,pto,alphas,fCKM2);
  const auto F3minusObjects = apfel::InitializeF3CCMinusObjectsSACOT(g,Thresholds,IntEps,nQ,Qmin,Qmax,intdeg);
  const std::map<int,apfel::Observable<>> F3minus = apfel::BuildStructureFunctions(F3minusObjects,PDFrotated,pto,alphas,fCKM2);
  const apfel::TabulateObject<apfel::Distribution> F3WP {[&] (double const& Q) -> apfel::Distribution{return 2*(F3plus.at(0).Evaluate(Q)+F3minus.at(0).Evaluate(Q));},nQ,Qmin,Qmax,intdeg,Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F3WM {[&] (double const& Q) -> apfel::Distribution{return -2*(F3plus.at(0).Evaluate(Q)-F3minus.at(0).Evaluate(Q));},nQ,Qmin,Qmax,intdeg,Thresholds};

  // dump to csv
  std::string filename = path + "/SACOT-chi_NLO.csv";

  std::ofstream file;
  file.open(filename);
  file<<std::fixed<<std::setprecision(10);
  file<<"x,Q2,NCF2,NCFL,NCF3,WMF2,WMFL,WMF3,WPF2,WPFL,WPF3\n";
  for(double xi: x)
    {
      for(double Qi: Q)
        {
          file<<xi<<","<<Qi*Qi;
          file<<","<<F2total.EvaluatexQ(xi,Qi);
          file<<","<<FLtotal.EvaluatexQ(xi,Qi);
          file<<","<<F3total.EvaluatexQ(xi,Qi)/xi;
          file<<","<<F2WM.EvaluatexQ(xi,Qi);
          file<<","<<FLWM.EvaluatexQ(xi,Qi);
          file<<","<<F3WM.EvaluatexQ(xi,Qi)/xi;
          file<<","<<F2WP.EvaluatexQ(xi,Qi);
          file<<","<<FLWP.EvaluatexQ(xi,Qi);
          file<<","<<F3WP.EvaluatexQ(xi,Qi)/xi;
          file<<"\n";
        }
    }
  file.close();

  //////////////////////////////
  //////////// NNLO ////////////
  //////////////////////////////

  pto = 2;
  double n = 1;

  // neutral current
  const auto F2objects_NNLO = apfel::InitializeF2NCObjectsASACOT(g,Thresholds,IntEps,nQ,Qmin,Qmax,intdeg,n);
  const std::map<int,apfel::Observable<>> F2_NNLO = apfel::BuildStructureFunctions(F2objects_NNLO,PDFrotated_NNLO,pto,alphas_NNLO,fEW);
  const apfel::TabulateObject<apfel::Distribution> F2total_NNLO {[&] (double const& Q) -> apfel::Distribution{return F2_NNLO.at(0).Evaluate(Q);},nQ,Qmin,Qmax,intdeg,Thresholds};

  const auto FLobjects_NNLO = apfel::InitializeFLNCObjectsASACOT(g,Thresholds,IntEps,nQ,Qmin,Qmax,intdeg,n);
  const std::map<int,apfel::Observable<>> FL_NNLO = apfel::BuildStructureFunctions(FLobjects_NNLO,PDFrotated_NNLO,pto,alphas_NNLO,fEW);
  const apfel::TabulateObject<apfel::Distribution> FLtotal_NNLO {[&] (double const& Q) -> apfel::Distribution{return FL_NNLO.at(0).Evaluate(Q);},nQ,Qmin,Qmax,intdeg,Thresholds};

  const auto F3objects_NNLO = apfel::InitializeF3NCObjectsASACOT(g,Thresholds,IntEps,nQ,Qmin,Qmax,intdeg);
  const std::map<int,apfel::Observable<>> F3_NNLO = apfel::BuildStructureFunctions(F3objects_NNLO,PDFrotated_NNLO,pto,alphas_NNLO,fPVEW);
  const apfel::TabulateObject<apfel::Distribution> F3total_NNLO {[&] (double const& Q) -> apfel::Distribution{return F3_NNLO.at(0).Evaluate(Q);},nQ,Qmin,Qmax,intdeg,Thresholds};

  // charged current
  const auto F2plusObjects_NNLO = apfel::InitializeF2CCPlusObjectsASACOT(g,Thresholds,IntEps,nQ,Qmin,Qmax,intdeg,n);
  const std::map<int,apfel::Observable<>> F2plus_NNLO = apfel::BuildStructureFunctions(F2plusObjects_NNLO,PDFrotated_NNLO,pto,alphas_NNLO,fCKM2);
  const auto F2minusObjects_NNLO = apfel::InitializeF2CCMinusObjectsASACOT(g,Thresholds,IntEps,nQ,Qmin,Qmax,intdeg,n);
  const std::map<int,apfel::Observable<>> F2minus_NNLO = apfel::BuildStructureFunctions(F2minusObjects_NNLO,PDFrotated_NNLO,pto,alphas_NNLO,fCKM2);
  const apfel::TabulateObject<apfel::Distribution> F2WP_NNLO {[&] (double const& Q) -> apfel::Distribution{return 2*(F2plus_NNLO.at(0).Evaluate(Q)+F2minus_NNLO.at(0).Evaluate(Q));},nQ,Qmin,Qmax,intdeg,Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F2WM_NNLO {[&] (double const& Q) -> apfel::Distribution{return 2*(F2plus_NNLO.at(0).Evaluate(Q)-F2minus_NNLO.at(0).Evaluate(Q));},nQ,Qmin,Qmax,intdeg,Thresholds};

  const auto FLplusObjects_NNLO = apfel::InitializeFLCCPlusObjectsASACOT(g,Thresholds,IntEps,nQ,Qmin,Qmax,intdeg,n);
  const std::map<int,apfel::Observable<>> FLplus_NNLO = apfel::BuildStructureFunctions(FLplusObjects_NNLO,PDFrotated_NNLO,pto,alphas_NNLO,fCKM2);
  const auto FLminusObjects_NNLO = apfel::InitializeFLCCMinusObjectsASACOT(g,Thresholds,IntEps,nQ,Qmin,Qmax,intdeg,n);
  const std::map<int,apfel::Observable<>> FLminus_NNLO = apfel::BuildStructureFunctions(FLminusObjects_NNLO,PDFrotated_NNLO,pto,alphas_NNLO,fCKM2);
  const apfel::TabulateObject<apfel::Distribution> FLWP_NNLO {[&] (double const& Q) -> apfel::Distribution{return 2*(FLplus_NNLO.at(0).Evaluate(Q)+FLminus_NNLO.at(0).Evaluate(Q));},nQ,Qmin,Qmax,intdeg,Thresholds};
  const apfel::TabulateObject<apfel::Distribution> FLWM_NNLO {[&] (double const& Q) -> apfel::Distribution{return 2*(FLplus_NNLO.at(0).Evaluate(Q)-FLminus_NNLO.at(0).Evaluate(Q));},nQ,Qmin,Qmax,intdeg,Thresholds};

  const auto F3plusObjects_NNLO = apfel::InitializeF3CCPlusObjectsASACOT(g,Thresholds,IntEps,nQ,Qmin,Qmax,intdeg);
  const std::map<int,apfel::Observable<>> F3plus_NNLO = apfel::BuildStructureFunctions(F3plusObjects_NNLO,PDFrotated_NNLO,pto,alphas_NNLO,fCKM2);
  const auto F3minusObjects_NNLO = apfel::InitializeF3CCMinusObjectsASACOT(g,Thresholds,IntEps,nQ,Qmin,Qmax,intdeg);
  const std::map<int,apfel::Observable<>> F3minus_NNLO = apfel::BuildStructureFunctions(F3minusObjects_NNLO,PDFrotated_NNLO,pto,alphas_NNLO,fCKM2);
  const apfel::TabulateObject<apfel::Distribution> F3WP_NNLO {[&] (double const& Q) -> apfel::Distribution{return 2*(F3plus_NNLO.at(0).Evaluate(Q)+F3minus_NNLO.at(0).Evaluate(Q));},nQ,Qmin,Qmax,intdeg,Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F3WM_NNLO {[&] (double const& Q) -> apfel::Distribution{return -2*(F3plus_NNLO.at(0).Evaluate(Q)-F3minus_NNLO.at(0).Evaluate(Q));},nQ,Qmin,Qmax,intdeg,Thresholds};

  // dump to csv
  filename = path + "/aSACOT-chi_NNLO_test.csv";

  file.open(filename);
  file<<std::fixed<<std::setprecision(10);
  file<<"x,Q2,NCF2,NCFL,NCF3,WMF2,WMFL,WMF3,WPF2,WPFL,WPF3\n";
  for(double xi: x)
    {
      for(double Qi: Q)
        {
          file<<xi<<","<<Qi*Qi;
          file<<","<<F2total_NNLO.EvaluatexQ(xi,Qi);
          file<<","<<FLtotal_NNLO.EvaluatexQ(xi,Qi);
          file<<","<<F3total_NNLO.EvaluatexQ(xi,Qi)/xi;
          file<<","<<F2WM_NNLO.EvaluatexQ(xi,Qi);
          file<<","<<FLWM_NNLO.EvaluatexQ(xi,Qi);
          file<<","<<F3WM_NNLO.EvaluatexQ(xi,Qi)/xi;
          file<<","<<F2WP_NNLO.EvaluatexQ(xi,Qi);
          file<<","<<FLWP_NNLO.EvaluatexQ(xi,Qi);
          file<<","<<F3WP_NNLO.EvaluatexQ(xi,Qi)/xi;
          file<<"\n";
        }
    }
  file.close();

  return 0;
}
