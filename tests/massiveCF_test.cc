//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>
#include <apfel/massivecoefficientfunctionsunp_sl.h>
#include <apfel/massivezerocoefficientfunctionsunp_sl.h>
#include <apfel/zeromasscoefficientfunctionsunp_sl.h>
#include <apfel/splittingfunctionspol_sl.h>
#include <apfel/matchingconditions_sl.h>
#include <apfel/betaqcd.h>

#include <iomanip>

// LH Toy PDFs
double xglu(double const& x)
{
  return 1.7 * pow(x,-0.1) * pow((1-x),5);
}

int main()
{
  /*
   * Test code to check the correctness and the accuracy of the
   * massive coefficient functions, both exact and asymptotic.
   */
  apfel::Timer t;
  t.start();

  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{100,1e-5,5}, apfel::SubGrid{60,1e-1,5}, apfel::SubGrid{50,6e-1,3}, apfel::SubGrid{50,8e-1,3}}};

  // Final scale
  const double Q2   = 200000;
  const double muF2 = Q2;
  const double M2   = 2;
  const double xi   = Q2 / M2;
  const double xiF  = M2 / muF2;
  const double lxi  = log(xi);
  const double lxi2 = lxi * lxi;
  const double lxiF = log(xiF);
  const double eta  = Q2 / ( Q2 + 4 * M2 );

  // Number of active flavours
  const int nf = 3;

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

  const apfel::Operator OmL1g  {g, apfel::CmL1gNC{eta},     IntEps};
  const apfel::Operator OmL2ns {g, apfel::CmL2nsNC{eta},    IntEps};
  const apfel::Operator OmL2gc {g, apfel::CmL2gNC{eta},     IntEps};
  const apfel::Operator OmL2gb {g, apfel::CmL2bargNC{eta},  IntEps};
  const apfel::Operator OmL2psc{g, apfel::CmL2psNC{eta},    IntEps};
  const apfel::Operator OmL2psb{g, apfel::CmL2barpsNC{eta}, IntEps};
  const apfel::Operator OmL2ps = OmL2psc + lxi * OmL2psb;
  const apfel::Operator OmL2g  = OmL2gc  + lxi * OmL2gb;

  // Massive zero coeffient functions
  const apfel::Operator Om021gc  {g, apfel::Cm021gNC_c{},   IntEps};
  const apfel::Operator Om021gl  {g, apfel::Cm021gNC_l{},   IntEps};
  const apfel::Operator Om022nsc {g, apfel::Cm022nsNC_c{},  IntEps};
  const apfel::Operator Om022nsl {g, apfel::Cm022nsNC_l{},  IntEps};
  const apfel::Operator Om022nsl2{g, apfel::Cm022nsNC_l2{}, IntEps};
  const apfel::Operator Om022psc {g, apfel::Cm022psNC_c{},  IntEps};
  const apfel::Operator Om022psl {g, apfel::Cm022psNC_l{},  IntEps};
  const apfel::Operator Om022psl2{g, apfel::Cm022psNC_l2{}, IntEps};
  const apfel::Operator Om022psf {g, apfel::Cm022psNC_f{},  IntEps};
  const apfel::Operator Om022pslf{g, apfel::Cm022psNC_lf{}, IntEps};
  const apfel::Operator Om022gc  {g, apfel::Cm022gNC_c{},   IntEps};
  const apfel::Operator Om022gl  {g, apfel::Cm022gNC_l{},   IntEps};
  const apfel::Operator Om022gl2 {g, apfel::Cm022gNC_l2{},  IntEps};
  const apfel::Operator Om022gf  {g, apfel::Cm022gNC_f{},   IntEps};
  const apfel::Operator Om022glf {g, apfel::Cm022gNC_lf{},  IntEps};
  const apfel::Operator Om021g  = Om021gc  + lxi * Om021gl;
  const apfel::Operator Om022ns = Om022nsc + lxi * Om022nsl + lxi2 * Om022nsl2;
  const apfel::Operator Om022ps = Om022psc + lxi * Om022psl + lxi2 * Om022psl2 + lxiF * Om022psf + lxi * lxiF * Om022pslf;
  const apfel::Operator Om022g  = Om022gc  + lxi * Om022gl  + lxi2 * Om022gl2  + lxiF * Om022gf  + lxi * lxiF * Om022glf;

  const apfel::Operator Om0L1g  {g, apfel::Cm0L1gNC_c{},  IntEps};
  const apfel::Operator Om0L2nsc{g, apfel::Cm0L2nsNC_c{}, IntEps};
  const apfel::Operator Om0L2nsl{g, apfel::Cm0L2nsNC_l{}, IntEps};
  const apfel::Operator Om0L2psc{g, apfel::Cm0L2psNC_c{}, IntEps};
  const apfel::Operator Om0L2psl{g, apfel::Cm0L2psNC_l{}, IntEps};
  const apfel::Operator Om0L2psf{g, apfel::Cm0L2psNC_f{}, IntEps};
  const apfel::Operator Om0L2gc {g, apfel::Cm0L2gNC_c{},  IntEps};
  const apfel::Operator Om0L2gl {g, apfel::Cm0L2gNC_l{},  IntEps};
  const apfel::Operator Om0L2gf {g, apfel::Cm0L2gNC_f{},  IntEps};
  const apfel::Operator Om0L2ns = Om0L2nsc + lxi * Om0L2nsl;
  const apfel::Operator Om0L2ps = Om0L2psc + lxi * Om0L2psl + lxiF * Om0L2psf;
  const apfel::Operator Om0L2g  = Om0L2gc  + lxi * Om0L2gl  + lxiF * Om0L2gf;

  /*
   * Massive0 coeffient functions. These coefficient functions are
   * explicitly reported in Appendix D of
   * https://arxiv.org/pdf/hep-ph/9601302.pdf. However, they can be
   * obtained by a proper combination of zero-mass coeffient
   * functions, splitting functions, matching conditions, and
   * beta-function coeffients. All these expressions are implemented
   * and they need to be combined using the following equations of the
   * same paper:
   *
   * - (4.6)  and (4.7)  for the O(as)   gluon        coefficient functions,
   * - (4.14) and (4.15) for the O(as^2) gluon        coefficient functions,
   * - (4.20) and (4.21) for the O(as^2) pure-singlet coefficient functions,
   * - (4.28) and (4.29) for the O(as^2) non-singlet  coefficient functions.
   *
   * The convolution will be performed numerically using the internal
   * functions. Interestingly, this way of computing the massive0
   * coeffient functions is particularly useful when implementing the
   * FONLL scheme because in the subtraction term (ZM-M0) the ZM
   * coefficient functions cancel explicitly.
   */

  // beta0 constants
  const double beta0  = apfel::beta0qcd(nf);
  const double beta0Q = beta0 - apfel::beta0qcd(nf-1);

  // Zero mass coefficient functions
  const apfel::Operator c21ns  {g, apfel::C21ns{},      IntEps};
  const apfel::Operator c21g   {g, apfel::C21g{},       IntEps};
  const apfel::Operator c22ps  {g, apfel::C22ps{},      IntEps};
  const apfel::Operator c22g   {g, apfel::C22g{},       IntEps};
  const apfel::Operator c22nsp1{g, apfel::C22nsp{nf},   IntEps};
  const apfel::Operator c22nsp2{g, apfel::C22nsp{nf+1}, IntEps};
  const apfel::Operator c22nsp = c22nsp2 - c22nsp1; // Select nf dependent part

  const apfel::Operator cL1ns  {g, apfel::CL1ns{},      IntEps};
  const apfel::Operator cL1g   {g, apfel::CL1g{},       IntEps};
  const apfel::Operator cL2ps  {g, apfel::CL2ps{},      IntEps};
  const apfel::Operator cL2g   {g, apfel::CL2g{},       IntEps};
  const apfel::Operator cL2nsp1{g, apfel::CL2nsp{nf},   IntEps};
  const apfel::Operator cL2nsp2{g, apfel::CL2nsp{nf+1}, IntEps};
  const apfel::Operator cL2nsp = cL2nsp2 - cL2nsp1; // Select nf dependent part

  // Splitting functions functions. Need to be multiplied by 2. This
  // is done in the expressions below.
  const apfel::Operator p0ns {g, apfel::P0ns{},      IntEps};
  const apfel::Operator p0qg {g, apfel::P0qg{1},     IntEps};
  const apfel::Operator p0gq {g, apfel::P0gq{},      IntEps};
  const apfel::Operator p0gg {g, apfel::P0gg{nf},    IntEps};
  const apfel::Operator p1ps1{g, apfel::P1ps{1},     IntEps};
  const apfel::Operator p1qg1{g, apfel::P1qg{1},     IntEps};
  const apfel::Operator p1ns1{g, apfel::P1nsp{nf},   IntEps};
  const apfel::Operator p1ns2{g, apfel::P1nsp{nf+1}, IntEps};
  const apfel::Operator p1ns = 2 * apfel::TR * ( p1ns2 - p1ns1 );  // Eq. (3.34)
  const apfel::Operator p1ps = 2 * apfel::TR * p1ps1;              // Eq. (3.30)
  const apfel::Operator p1qg = 2 * apfel::TR * p1qg1;              // Eq. (3.23)

  // Matching conditions (only non-log terms needed)
  const apfel::Operator a1Qg   {g, apfel::Null{},      IntEps};
  //const apfel::Operator a1Qgb = - zeta2 * p0qg / 4;
  const apfel::Operator aPS2Qq {g, apfel::APS2Hq_0{},  IntEps};
  const apfel::Operator aNS2qqQ{g, apfel::ANS2qqH_0{}, IntEps};
  const apfel::Operator a2Qg   {g, apfel::AS2Hg_0{},   IntEps};

  // Eqs. (4.6) and (4.7) for the O(as) gluon coefficient functions (!)
  const apfel::Operator HL1g = cL1g;
  const apfel::Operator H21g = p0qg * lxi + a1Qg + c21g;

  // Eqs. (4.14) and (4.15) for the O(as^2) gluon coefficient
  // functions.
  const apfel::Operator HL2g = ( - beta0 * cL1g + p0gg * cL1g + p0qg * cL1ns ) * lxi
                               + ( - beta0 * cL1g + p0gg * cL1g ) * lxiF + cL2g + a1Qg * cL1ns;
  const apfel::Operator H22g = ( p0qg * ( p0gg + p0ns ) / 2 - beta0 * p0qg / 2 ) * lxi2
                               + ( p1qg - beta0 * c21g + p0ns * a1Qg + p0gg * c21g + p0qg * c21ns ) * lxi
                               + ( p0qg * p0gg - beta0 * p0qg ) * lxi *lxiF
                               + ( - beta0 * ( c21g + a1Qg ) + p0gg * ( c21g + a1Qg ) ) * lxiF
                               + c22g + a2Qg;
  //+ 2 * beta0 * a1Qgb + c21ns * a1Qg + 2 * p0ns * a1Qgb -  2 * p0gg * a1Qgb; // These terms are already included in 'a2Qg'

  // Eqs. (4.20) and (4.21) for the O(as^2) pure-singlet coefficient
  // functions.
  const apfel::Operator HL2q = ( lxi + lxiF ) * p0gq * cL1g + cL2ps;
  const apfel::Operator H22q =  p0qg * p0gq * ( lxi2 / 2 + lxi * lxiF )
                                + ( p1ps + p0gq * c21g ) * lxi
                                + p0gq * ( a1Qg + c21g ) * lxiF + c22ps + aPS2Qq;
  //- 2 * p0gq * a1Qgb; // This term is already included in 'aPS2Qq'

  // Eqs. (4.28) and (4.29) for the O(as^2) non-singlet coefficient
  // functions (!).
  const apfel::Operator LL2NSq = - beta0Q * lxi * cL1ns + cL2nsp;
  const apfel::Operator L22NSq = - beta0Q * lxi2 * p0ns / 2 + ( p1ns - beta0Q * c21ns ) * lxi
                                 + c22nsp + aNS2qqQ;
  //+ beta0Q * zeta2 * p0ns / 2; // This term is already included in 'aNS2qqQ'

  // Define test distribution
  const apfel::Distribution TestDist{g, xglu};

  // Multiply the operators by the test distribution
  const apfel::Distribution Dm21g   = Om21g  * TestDist;
  const apfel::Distribution DmL1g   = OmL1g  * TestDist;
  const apfel::Distribution Dm22ns  = Om22ns * TestDist;
  const apfel::Distribution DmL2ns  = OmL2ns * TestDist;
  const apfel::Distribution Dm22ps  = Om22ps * TestDist;
  const apfel::Distribution DmL2ps  = OmL2ps * TestDist;
  const apfel::Distribution Dm22g   = Om22g  * TestDist;
  const apfel::Distribution DmL2g   = OmL2g  * TestDist;

  const apfel::Distribution Hm21g   = H21g   * TestDist;
  const apfel::Distribution HmL1g   = HL1g   * TestDist;
  const apfel::Distribution Hm22ns  = L22NSq * TestDist;
  const apfel::Distribution HmL2ns  = LL2NSq * TestDist;
  const apfel::Distribution Hm22ps  = H22q   * TestDist;
  const apfel::Distribution HmL2ps  = HL2q   * TestDist;
  const apfel::Distribution Hm22g   = H22g   * TestDist;
  const apfel::Distribution HmL2g   = HL2g   * TestDist;

  const apfel::Distribution Dm021g  = Om021g  * TestDist;
  const apfel::Distribution Dm0L1g  = Om0L1g  * TestDist;
  const apfel::Distribution Dm022ns = Om022ns * TestDist;
  const apfel::Distribution Dm0L2ns = Om0L2ns * TestDist;
  const apfel::Distribution Dm022ps = Om022ps * TestDist;
  const apfel::Distribution Dm0L2ps = Om0L2ps * TestDist;
  const apfel::Distribution Dm022g  = Om022g  * TestDist;
  const apfel::Distribution Dm0L2g  = Om0L2g  * TestDist;
  t.stop();

  const std::vector<double> xlha = {1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};

  std::cout << std::scientific << std::endl;
  std::cout << "xi = " << xi << std::endl;
  std::cout << "\nO(as) gluon coefficient functions" << std::endl;
  std::cout << "Asymptotic constructed vs. asymptotic explicit" << std::endl;
  std::cout << "    x    \t"
            << "  Cm21g  \t"
            << "  CmL1g  \t"
            << std::endl;
  for (double const& x : xlha)
    {
      std::cout << std::setprecision(4) << x << "\t"
                << Hm21g.Evaluate(x) / Dm021g.Evaluate(x) << "\t"
                << HmL1g.Evaluate(x) / Dm0L1g.Evaluate(x) << "\t"
                << std::endl;
    }
  std::cout << std::endl;
  std::cout << "Exact vs. asymptotic" << std::endl;
  std::cout << "    x    \t"
            << "  Cm21g  \t"
            << "  CmL1g  \t"
            << std::endl;
  for (double const& x : xlha)
    {
      std::cout << std::setprecision(4) << x << "\t"
                << Dm21g.Evaluate(x) / Dm021g.Evaluate(x) << "\t"
                << DmL1g.Evaluate(x) / Dm0L1g.Evaluate(x) << "\t"
                << std::endl;
    }
  std::cout << std::endl;

  std::cout << "O(as^2) non-singlet coefficient functions" << std::endl;
  std::cout << "Asymptotic constructed vs. asymptotic explicit" << std::endl;
  std::cout << "    x    \t"
            << "  Cm22ns \t"
            << "  CmL2ns \t"
            << std::endl;
  for (double const& x : xlha)
    {
      std::cout << std::setprecision(4) << x << "\t"
                << Hm22ns.Evaluate(x) / Dm022ns.Evaluate(x) << "\t"
                << HmL2ns.Evaluate(x) / Dm0L2ns.Evaluate(x) << "\t"
                << std::endl;
    }
  std::cout << std::endl;
  std::cout << "Exact vs. asymptotic" << std::endl;
  std::cout << "    x    \t"
            << "  Cm22ns \t"
            << "  CmL2ns \t"
            << std::endl;
  for (double const& x : xlha)
    {
      std::cout << std::setprecision(4) << x << "\t"
                << Dm22ns.Evaluate(x) / Dm022ns.Evaluate(x) << "\t"
                << DmL2ns.Evaluate(x) / Dm0L2ns.Evaluate(x) << "\t"
                << std::endl;
    }
  std::cout << std::endl;

  std::cout << "O(as^2) pure-singlet coefficient functions" << std::endl;
  std::cout << "Asymptotic constructed vs. asymptotic explicit" << std::endl;
  std::cout << "    x    \t"
            << "  Cm22ps \t"
            << "  CmL2ps \t"
            << std::endl;
  for (double const& x : xlha)
    {
      std::cout << std::setprecision(4) << x << "\t"
                << Hm22ps.Evaluate(x) / Dm022ps.Evaluate(x) << "\t"
                << HmL2ps.Evaluate(x) / Dm0L2ps.Evaluate(x) << "\t"
                << std::endl;
    }
  std::cout << std::endl;
  std::cout << "Exact vs. asymptotic" << std::endl;
  std::cout << "    x    \t"
            << "  Cm22ps \t"
            << "  CmL2ps \t"
            << std::endl;
  for (double const& x : xlha)
    {
      std::cout << std::setprecision(4) << x << "\t"
                << Dm22ps.Evaluate(x) / Dm022ps.Evaluate(x) << "\t"
                << DmL2ps.Evaluate(x) / Dm0L2ps.Evaluate(x) << "\t"
                << std::endl;
    }
  std::cout << std::endl;

  std::cout << "O(as^2) gluon coefficient functions" << std::endl;
  std::cout << "Asymptotic constructed vs. asymptotic explicit" << std::endl;
  std::cout << "    x    \t"
            << "  Cm22g \t"
            << "  CmL2g \t"
            << std::endl;
  for (double const& x : xlha)
    {
      std::cout << std::setprecision(4) << x << "\t"
                << Hm22g.Evaluate(x) / Dm022g.Evaluate(x) << "\t"
                << HmL2g.Evaluate(x) / Dm0L2g.Evaluate(x) << "\t"
                << std::endl;
    }
  std::cout << std::endl;
  std::cout << "Exact vs. asymptotic" << std::endl;
  std::cout << "    x    \t"
            << "  Cm22g \t"
            << "  CmL2g \t"
            << std::endl;
  for (double const& x : xlha)
    {
      std::cout << std::setprecision(4) << x << "\t"
                << Dm22g.Evaluate(x) / Dm022g.Evaluate(x) << "\t"
                << DmL2g.Evaluate(x) / Dm0L2g.Evaluate(x) << "\t"
                << std::endl;
    }
  std::cout << std::endl;

  return 0;
}
