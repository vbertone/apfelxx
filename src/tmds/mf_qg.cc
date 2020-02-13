//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/matchingfunctionspdf.h"
#include "apfel/constants.h"

#include <ginac/ginac.h>

namespace apfel
{
/////////////////////////////////////////////////////
  /*
    class C1qgpdf_new: public Expression
    {
    public:
      C1qgpdf_new(double const& Lp = 0);
      double Regular(double const& x) const;
    protected:
      double const _Lp;
    };

    class C2qgpdf_new: public Expression
    {
    public:
      C2qgpdf_new(double const& Lp = 0, double const& Lh = 0);
      double Regular(double const& x) const;
    protected:
      double const _Lp;
      double const _Lh;
    };

    class C3qgpdf: public Expression
    {
    public:
      C3qgpdf(int const& nf, double const& Lp = 0, double const& Lh = 0);
      double Regular(double const& x) const;
    protected:
      int    const _nf;
      double const _Lp;
      double const _Lh;
    };
  */
/////////////////////////////////////////////////////
  /*
    //_________________________________________________________________________________
    C1qgpdf::C1qgpdf(double const& Lp):
      Expression(),
      _Lp(Lp)
    {
    }
    double C1qgpdf::Regular(double const& z) const
    {
      return 2 * TR * ( 4 * ( 1 - z ) * z - 2 * _Lp * ( 1 - 2 * z + 2 * z * z ) );
    }
  */
  /*
    //_________________________________________________________________________________
    C2qgpdf::C2qgpdf(double const& Lp, double const& Lh):
      Expression(),
      _Lp(Lp),
      _Lh(Lh)
    {
    }
    double C2qgpdf::Regular(double const& z) const
    {
      // Symbol for GiNaC
      const GiNaC::symbol sz{"z"};

      // Useful definitions
      const double z2  = z * z;
      const double z3  = z * z2;
      const double omz = 1 - z;
      const double opz = 1 + z;

      // Logs of the scales
      const double Lp2 = _Lp * _Lp;
      const double Lp3 = _Lp * Lp2;

      // Harmonic polylogs:
      // Weight 1
      const double H0  = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{0},sz == z)).to_double();
      const double H1  = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1},sz == z)).to_double();
      const double Hm1 = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1},sz == z)).to_double();

      // Weight 2
      const double H00  = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{0,0},sz == z)).to_double();
      const double H10  = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,0},sz == z)).to_double();
      const double Hm10 = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,0},sz == z)).to_double();
      const double H11  = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,1},sz == z)).to_double();
      const double H2   = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{2},sz == z)).to_double();

      // Weight 3
      const double H000   = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{0,0,0},sz == z)).to_double();
      const double H100   = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,0,0},sz == z)).to_double();
      const double Hm100  = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,0,0},sz == z)).to_double();
      const double H20    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{2,0},sz == z)).to_double();
      const double Hm20   = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-2,0},sz == z)).to_double();
      const double H21    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{2,1},sz == z)).to_double();
      const double H12    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,2},sz == z)).to_double();
      const double H110   = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,1,0},sz == z)).to_double();
      const double H111   = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,1,1},sz == z)).to_double();
      const double Hm1m10 = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,-1,0},sz == z)).to_double();
      const double cf =
        CA*TR*((-2*(-172 + 315*z - 387*z2 + 298*z3))/(27.*z) +
  	     (4*(21 - 30*z + 68*z2)*H0)/9. + 2*z*(-3 + 4*z)*H1 +
  	     Lp2*((2*omz*(4 + 7*z + 31*z2))/(3.*z) + 4*(1 + 4*z)*H0 -
  		  4*(1 - 2*z + 2*z2)*H1) + 8*z*opz*Hm10 -
  	     (2*(3 - 12*z + 44*z2)*H00)/3. - (8*omz*(2 - z + 11*z2)*H10)/(3.*z) -
  	     8*omz*z*H11 + 4*(1 + 2*z)*H000 +
  	     4*(1 - 2*z + 2*z2)*(H12 + H110 - H111) +
  	     (8*(-2 + 3*z - 9*z2 + 11*z3)*zeta2)/(3.*z) +
  	     _Lp*((4*(-26 + 36*z - 135*z2 + 116*z3))/(9.*z) - (4*(3 + 44*z2)*H0)/3. +
  		 8*(1 + 2*z + 2*z2)*Hm10 + 8*(1 + 2*z)*H00 + 8*(1 - 2*z + 2*z2)*H11 +
  		 16*z*zeta2) + 4*(1 + 2*z + 2*z2)*
  	     (2*Hm20 - 2*Hm1m10 + Hm100 - Hm1*zeta2) - 8*z*(2*H20 - zeta3)) +
        CF*TR*(-13 + 75*z - 72*z2 + 2*Lp3*(1 - 2*z + 2*z2) +
  	     (8 + 15*z - 8*z2)*H0 - 2*z*(-3 + 4*z)*H1 +
  	     Lp2*(-7 + 12*z - 8*z2 + 4*_Lh*(1 - 2*z + 2*z2) - 2*(1 - 2*z + 4*z2)*H0 -
  		  4*(1 - 2*z + 2*z2)*H1) + (1 + 12*z - 8*z2)*H00 -
  	     2*(1 - 2*z + 4*z2)*H000 + 4*omz*z*(2*H2 + 2*H10 + 2*H11 - 3*zeta2) +
  	     _Lp*(-8*_Lh*omz*z - 2*(10 - 33*z + 28*z2) - 2*(1 - 8*z + 8*z2)*H0 +
  		 16*omz*z*H1 - 4*(1 - 2*z + 4*z2)*H00 +
  		 (4*(1 - 2*z + 2*z2)*(-18*H2 - 18*H10 - 18*H11 + (45*zeta2)/2.))/9.) +
  	     4*(1 - 2*z + 2*z2)*(H21 - H100 + H111 + 7*zeta3));
      return 2 * cf;
    }
  */
  /*
    //_________________________________________________________________________________
    C3qgpdf::C3qgpdf(int const& nf, double const& Lp, double const& Lh):
      Expression(),
      _nf(nf),
      _Lp(Lp),
      _Lh(Lh)
    {
    }
    double C3qgpdf::Regular(double const& z) const
    {
      // Symbol for GiNaC
      const GiNaC::symbol sz{"z"};

      // Useful definitions
      const double CF2 = CF * CF;
      const double CA2 = CA * CA;
      const double TR2 = TR * TR;

      const double z2   = z * z;
      const double z3   = z * z2;
      const double z4   = z * z3;
      const double z5   = z * z4;
      const double omz  = 1 - z;
      const double opz  = 1 + z;
      const double opz2 = opz * opz;

      // Logs of the scales
      const double Lp2 = _Lp * _Lp;
      const double Lp3 = _Lp * Lp2;
      const double Lp4 = _Lp * Lp3;
      const double Lp5 = _Lp * Lp4;
      const double Lh2 = _Lh * _Lh;

      // Harmonic polylogs:
      // Weight 1
      const double H0  = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{0},sz == z)).to_double();
      const double H1  = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1},sz == z)).to_double();
      const double Hm1 = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1},sz == z)).to_double();

      // Weight 2
      const double H00   = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{0,0},sz == z)).to_double();
      const double H10   = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,0},sz == z)).to_double();
      const double Hm10  = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,0},sz == z)).to_double();
      const double H11   = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,1},sz == z)).to_double();
      const double Hm1m1 = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,-1},sz == z)).to_double();
      const double H2    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{2},sz == z)).to_double();
      const double Hm2   = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-2},sz == z)).to_double();

      // Weight 3
      const double H000    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{0,0,0},sz == z)).to_double();
      const double H100    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,0,0},sz == z)).to_double();
      const double Hm100   = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,0,0},sz == z)).to_double();
      const double H20     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{2,0},sz == z)).to_double();
      const double Hm20    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-2,0},sz == z)).to_double();
      const double H3      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{3},sz == z)).to_double();
      const double Hm3     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-3},sz == z)).to_double();
      const double H21     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{2,1},sz == z)).to_double();
      const double H12     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,2},sz == z)).to_double();
      const double H110    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,1,0},sz == z)).to_double();
      const double H111    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,1,1},sz == z)).to_double();
      const double Hm1m10  = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,-1,0},sz == z)).to_double();
      const double Hm2m1   = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-2,-1},sz == z)).to_double();
      const double Hm1m2   = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,-2},sz == z)).to_double();
      const double Hm12    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,2},sz == z)).to_double();
      const double Hm1m1m1 = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,-1,-1},sz == z)).to_double();
      const double H1m2    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,-2},sz == z)).to_double();

      // Weight 4
      const double H0000    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{0,0,0,0},sz == z)).to_double();
      const double H200     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{2,0,0},sz == z)).to_double();
      const double Hm30     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-3,0},sz == z)).to_double();
      const double Hm2m10   = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-2,-1,0},sz == z)).to_double();
      const double Hm200    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-2,0,0},sz == z)).to_double();
      const double Hm1m20   = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,-2,0},sz == z)).to_double();
      const double Hm1m1m10 = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,-1,-1,0},sz == z)).to_double();
      const double Hm1m100  = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,-1,0,0},sz == z)).to_double();
      const double Hm1000   = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,0,0,0},sz == z)).to_double();
      const double H13      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,3},sz == z)).to_double();
      const double H112     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,1,2},sz == z)).to_double();
      const double H120     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,2,0},sz == z)).to_double();
      const double H121     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,2,1},sz == z)).to_double();
      const double H1000    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,0,0,0},sz == z)).to_double();
      const double H1100    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,1,0,0},sz == z)).to_double();
      const double H1110    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,1,1,0},sz == z)).to_double();
      const double H1111    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,1,1,1},sz == z)).to_double();
      const double H30      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{3,0},sz == z)).to_double();
      const double H210     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{2,1,0},sz == z)).to_double();
      const double H211     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{2,1,1},sz == z)).to_double();
      const double H4       = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{4},sz == z)).to_double();
      const double H22      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{2,2},sz == z)).to_double();
      const double H31      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{3,1},sz == z)).to_double();
      const double Hm120    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,2,0},sz == z)).to_double();
      const double H1m20    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,-2,0},sz == z)).to_double();
      const double Hm22     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-2,2},sz == z)).to_double();
      const double Hm13     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,3},sz == z)).to_double();
      const double Hm1m12   = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,-1,2},sz == z)).to_double();
      const double Hm121    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,2,1},sz == z)).to_double();

      // Weight 5
      const double H00000     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{0,0,0,0,0},sz == z)).to_double();
      const double Hm40       = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-4,0},sz == z)).to_double();
      const double H23        = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{2,3},sz == z)).to_double();
      const double H32        = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{3,2},sz == z)).to_double();
      const double Hm300      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-3,0,0},sz == z)).to_double();
      const double Hm2m20     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-2,-2,0},sz == z)).to_double();
      const double H2m20      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{2,-2,0},sz == z)).to_double();
      const double H212       = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{2,1,2},sz == z)).to_double();
      const double H220       = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{2,2,0},sz == z)).to_double();
      const double H221       = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{2,2,1},sz == z)).to_double();
      const double H300       = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{3,0,0},sz == z)).to_double();
      const double H310       = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{3,1,0},sz == z)).to_double();
      const double H311       = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{3,1,1},sz == z)).to_double();
      const double Hm2m100    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-2,-1,0,0},sz == z)).to_double();
      const double Hm2000     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-2,0,0,0},sz == z)).to_double();
      const double H2000      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{2,0,0,0},sz == z)).to_double();
      const double H2100      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{2,1,0,0},sz == z)).to_double();
      const double H2110      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{2,1,1,0},sz == z)).to_double();
      const double H2111      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{2,1,1,1},sz == z)).to_double();
      const double Hm3m10     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-3,-1,0},sz == z)).to_double();
      const double Hm2m1m10   = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-2,-1,-1,0},sz == z)).to_double();
      const double H5         = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{5},sz == z)).to_double();
      const double H40        = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{4,0},sz == z)).to_double();
      const double H41        = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{4,1},sz == z)).to_double();
      const double Hm1m30     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,-3,0},sz == z)).to_double();
      const double Hm130      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,3,0},sz == z)).to_double();
      const double Hm1m2m10   = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,-2,-1,0},sz == z)).to_double();
      const double Hm1m200    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,-2,0,0},sz == z)).to_double();
      const double Hm1m1m20   = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,-1,-2,0},sz == z)).to_double();
      const double Hm1200     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,2,0,0},sz == z)).to_double();
      const double Hm1210     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,2,1,0},sz == z)).to_double();
      const double Hm1m1m1m10 = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,-1,-1,-1,0},sz == z)).to_double();
      const double Hm1m1m100  = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,-1,-1,0,0},sz == z)).to_double();
      const double Hm1m1000   = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,-1,0,0,0},sz == z)).to_double();
      const double Hm10000    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,0,0,0,0},sz == z)).to_double();
      const double H14        = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,4},sz == z)).to_double();
      const double H1m30      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,-3,0},sz == z)).to_double();
      const double H113       = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,1,3},sz == z)).to_double();
      const double H122       = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,2,2},sz == z)).to_double();
      const double H130       = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,3,0},sz == z)).to_double();
      const double H131       = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,3,1},sz == z)).to_double();
      const double H1m200     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,-2,0,0},sz == z)).to_double();
      const double H11m20     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,1,-2,0},sz == z)).to_double();
      const double H1112      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,1,1,2},sz == z)).to_double();
      const double H1120      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,1,2,0},sz == z)).to_double();
      const double H1121      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,1,2,1},sz == z)).to_double();
      const double H1200      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,2,0,0},sz == z)).to_double();
      const double H1210      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,2,1,0},sz == z)).to_double();
      const double H1211      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,2,1,1},sz == z)).to_double();
      const double H10000     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,0,0,0,0},sz == z)).to_double();
      const double H11000     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,1,0,0,0},sz == z)).to_double();
      const double H11100     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,1,1,0,0},sz == z)).to_double();
      const double H11110     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,1,1,1,0},sz == z)).to_double();
      const double H11111     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,1,1,1,1},sz == z)).to_double();
      const double Hm32       = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-3,2},sz == z)).to_double();
      const double Hm23       = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-2,3},sz == z)).to_double();
      const double Hm2m12     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-2,-1,2},sz == z)).to_double();
      const double Hm220      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-2,2,0},sz == z)).to_double();
      const double Hm14       = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,4},sz == z)).to_double();
      const double Hm221      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-2,2,1},sz == z)).to_double();
      const double Hm1m22     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,-2,2},sz == z)).to_double();
      const double Hm1m13     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,-1,3},sz == z)).to_double();
      const double Hm122      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,2,2},sz == z)).to_double();
      const double Hm131      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,3,1},sz == z)).to_double();
      const double Hm1m1m12   = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,-1,-1,2},sz == z)).to_double();
      const double Hm1m121    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,-1,2,1},sz == z)).to_double();
      const double Hm1211     = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,2,1,1},sz == z)).to_double();
      const double H1m22      = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,-2,2},sz == z)).to_double();
      const double H1m2m10    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{1,-2,-1,0},sz == z)).to_double();
      const double Hm1m120    = GiNaC::ex_to<GiNaC::numeric>(GiNaC::H(GiNaC::lst{-1,-1,2,0},sz == z)).to_double();
      const double cf =
        CA*_nf*TR2*((4*(-1582 - 8697*z - 24198*z2 + 43324*z3))/(729.*z) +
                    (8*(-523 - 8119*z + 2280*z2)*H0)/243. +
                    (8*(-36 + 103*z - 1064*z2 + 668*z3)*H1)/(243.*z) +
                    Lp3*((-8*omz*(4 + 7*z + 31*z2))/(27.*z) - (16*(1 + 4*z)*H0)/9. +
                         (16*(1 - 2*z + 2*z2)*H1)/9.) + (16*(-3 - 23*z + 3*z2)*H2)/27. -
                    (16*(10 + 21*z + 20*z2)*Hm20)/9. - (8*(85 + 290*z + 299*z2)*Hm10)/81. -
                    (8*(32 + 1407*z + 391*z2)*H00)/81. - (8*(85 - 326*z + 290*z2)*H11)/81. +
                    (8*(12 + 65*z + 8*z2)*H20)/9. - (8*z*(5*H3 - 2*H21))/9. -
                    (16*(25 + 62*z + 62*z2)*Hm100)/27. + (32*(-4 - 74*z + 23*z2)*H000)/27. -
                    (32*(10 - 23*z + 23*z2)*(H12 - H111))/27. +
                    (16*(1 + 6*z + 2*z2)*H200)/3. - (32*(-1 + 12*z)*H0000)/9. +
                    (8*(12 + 27*z + 8*z2)*H0*zeta2)/9. +
                    Lp2*((4*(-40 + 6*z - 348*z2 + 409*z3))/(27.*z) +
                         (16*(-5 - 32*z + 9*z2)*H0)/9. + (80*(1 - 2*z + 2*z2)*H1)/9. -
                         (16*(1 + 2*z + 2*z2)*Hm10)/3. - (64*z*H00)/3. -
                         (16*(1 - 2*z + 2*z2)*H11)/3. - (32*z*zeta2)/3.) +
                    (-24*(-40 + 3*z - 123*z2 + 154*z3)*H10 -
                     8*(-120 - 18*z - 172*z2 + 435*z3)*zeta2)/(81.*z) +
                    (80*(1 + 2*z + 2*z2)*(2*Hm1m10 + Hm1*zeta2))/9. +
                    _Lp*((4*(48 - 23*z - 976*z2 + 1152*z3))/(27.*z) +
                        (16*(17 - 250*z + 83*z2)*H0)/27. + (8*(38 - 121*z + 112*z2)*H1)/27. -
                        (16*z*H2)/3. - (160*(1 + 2*z + 2*z2)*Hm10)/9. +
                        (8*(-13 - 114*z + 28*z2)*H00)/9. + (32*omz*(2 - z + 11*z2)*H10)/(9.*z) -
                        (160*(1 - 2*z + 2*z2)*H11)/9. + (64*z*H20)/3. - (16*(-1 + 10*z)*H000)/3. -
                        (16*(-1 + 2*z)*(4 + 2*z + 11*z2)*zeta2)/(9.*z) +
                        (4*(1 + 2*z + 2*z2)*(-24*Hm20 + 24*Hm1m10 - 36*Hm100 + 12*Hm1*zeta2))/
                        9. + (4*(1 - 2*z + 2*z2)*(-24*H12 + 12*H100 + 24*H111 + 12*H1*zeta2))/
                        9. + (64*z*zeta3)/3.) + (-16*(-24 + 31*z - 128*z2 + 116*z3)*H100 +
                                                 32*(-6 + 4*z - 29*z2 + 26*z3)*H110 +
                                                 16*(-12 + 23*z - 88*z2 + 82*z3)*H1*zeta2 +
                                                 48*(-4 + 18*z + 7*z2 + 30*z3)*zeta3)/(27.*z) -
                    (32*(1 + 2*z + 2*z2)*(3*Hm30 - 3*Hm2m10 + (9*Hm200)/2. - 3*Hm1m20 +
                                          3*Hm1m1m10 - (9*Hm1m100)/2. + 4*Hm1000 - (3*Hm2*zeta2)/2. +
                                          (3*Hm1m1*zeta2)/2. - (3*Hm10*zeta2)/2. - (3*Hm1*zeta3)/2.))/9. -
                    (16*(1 - 2*z + 2*z2)*(2*H13 - 8*H112 - 2*H120 - 5*H121 - 3*H1000 +
                                          4*H1100 - 4*H1110 + 5*H1111 + H10*zeta2 + 2*H11*zeta2 -
                                          10*H1*zeta3))/9. + (16*z*(12*H30 - 12*H210 - 12*H2*zeta2 - 19*zeta4))/9.) +
        CF*_nf*TR2*((448*_Lh*omz*z)/27. - (16*Lp4*(1 - 2*z + 2*z2))/9. -
                    (-160736 - 1294539*z + 1012650*z2 + 417704*z3)/(1458.*z) +
                    (8*(20752 + 7693*z + 17172*z2)*H0)/243. + (8*(-355 + 155*z + 268*z2)*H1)/243. +
                    (224*(4 - 11*z + 11*z2)*H2)/81. - (4*(-6023 + 6559*z + 7424*z2)*H00)/81. +
                    Lp3*((-8*_Lh*(1 - 2*z + 2*z2))/3. +
                         (4*(-16 + 165*z - 264*z2 + 94*z3))/(27.*z) - (8*(-11 + 4*z)*H0)/9. +
                         (32*(1 - 2*z + 2*z2)*H1)/9. - (16*(-1 + 2*z)*H00)/3.) +
                    (8*(112 - 353*z + 308*z2)*H11)/81. + (4*(1805 - 1282*z + 1344*z2)*H000)/27. +
                    Lp2*((-16*_Lh*(5 - 13*z + 13*z2))/9. +
                         (2*(256 + 3099*z - 3528*z2 + 272*z3))/(27.*z) +
                         (4*(385 + 58*z + 192*z2)*H0)/9. + (16*(5 - 16*z + 16*z2)*H1)/9. -
                         (32*(-5 + 7*z + 2*z2)*H00)/3. +
                         (4*(1 - 2*z + 2*z2)*(12*H2 + 24*H10 + 12*H11))/9. - 32*(-1 + 2*z)*H000) -
                    (32*(10 - 23*z + 23*z2)*(H21 + H111))/27. -
                    (16*(-47 + 94*z + 44*z2)*H0000)/9. +
                    (80*(1 - 2*z + 2*z2)*(H211 + 2*H1000 + H1111))/9. -
                    64*(-1 + 2*z)*H00000 + (32*(-96 + 131*z - 403*z2 + 382*z3)*H10 +
                                            8*(-384 + 356*z - 961*z2 + 877*z3)*zeta2)/(81.*z) -
                    (32*z*(-9 + 10*z)*(H20 + H0*zeta2))/9. -
                    (32*omz*(4 + 13*z + 22*z2)*(H110 + H1*zeta2))/(9.*z) +
                    _Lp*((-32*_Lh*(7 - 29*z + 29*z2))/27. +
                        (-5504 + 56933*z - 77614*z2 + 28218*z3)/(81.*z) -
                        (4*(-2733 - 630*z + 1208*z2)*H0)/27. + (16*z*(-17 + 20*z)*H1)/9. +
                        (8*(421 - 41*z + 232*z2)*H00)/9. -
                        (32*(-4 - 14*z + 7*z2 + 6*z3)*H10)/(9.*z) +
                        (4*(1 - 2*z + 2*z2)*(40*H2 + 40*H11))/9. + (32*(3 + 6*z + 4*z2)*H20)/3. -
                        (8*(-43 + 68*z + 32*z2)*H000)/3. +
                        (4*(1 - 2*z + 2*z2)*(-24*H21 + 48*H100 - 24*H111))/9. -
                        80*(-1 + 2*z)*H0000 - (4*(-32 - 37*z - 118*z2 + 222*z3)*zeta2)/(9.*z) +
                        (32*(3 + 6*z + 4*z2)*H0*zeta2)/3. - (8*(-23 - 242*z + 2*z2)*zeta3)/9.) +
                    (-16*(-48 - 133*z - 28*z2 + 184*z3)*H100 -
                     16*(24 + 154*z - 287*z2 + 221*z3)*zeta3)/(27.*z) +
                    (32*(3 + 6*z + 4*z2)*(H30 + 2*H200 - H210 - H2*zeta2 + H00*zeta2 +
                                          H0*zeta3))/3. - (8*(73 + 106*z + 104*z2)*zeta4)/9.) +
        CF2*TR*((2127 + 2406*z - 5020*z2)/24. + Lp5*(-1 + 2*z - 2*z2) +
                ((235 - 734*z - 708*z2)*H0)/3. - (2*(-505 + 191*z + 354*z2)*H1)/3. +
                Lp4*(7 - 14*z + 10*z2 - 4*_Lh*(1 - 2*z + 2*z2) + 2*(1 - 2*z + 4*z2)*H0 +
                     4*(1 - 2*z + 2*z2)*H1) - (4*(-70 - 343*z + 151*z2)*H2)/3. +
                (4*(49 - 11*z + 24*z2)*H3)/3. - (4*(5 - 136*z + 50*z2)*H4)/3. +
                (64*(7 + 16*z2)*Hm40)/3. + (16*(47 + 46*z + 53*z2)*Hm30)/3. +
                (8*(135 + 98*z + 63*z2)*Hm20)/3. + (4*(395 + 634*z + 251*z2)*Hm10)/3. +
                ((15 - 857*z - 1012*z2)*H00)/3. - (4*(77 - 26*z + 2*z2)*H10)/3. -
                (2*(110 - 353*z + 302*z2)*H11)/3. + (16*(2 - 9*z + 6*z2)*H12)/3. -
                (8*(-7 - 12*z + 25*z2)*H13)/3. - (4*(-7 + 108*z + 4*z2)*H20)/3. +
                (4*(35 - 47*z + 27*z2)*H21)/3. + (8*(2 + 13*z2)*H22)/3. -
                (32*(3 - 8*z + 3*z2)*H23)/3. + (8*(7 - 98*z + 33*z2)*H30)/3. -
                (8*(-3 - 8*z + 6*z2)*H31)/3. + (64*(1 - 2*z + 3*z2)*H32)/3. +
                (32*(11 + 18*z + 26*z2)*Hm300)/3. - (64*(3 + 4*z + 8*z2)*Hm2m20)/3. +
                (8*(55 + 138*z + 112*z2)*Hm200)/3. - (112*opz*(3 + 7*z)*Hm1m20)/3. +
                (4*(183 + 384*z + 209*z2)*Hm100)/3. + 64*z*opz*Hm120 +
                ((65 - 694*z - 792*z2)*H000)/3. + (32*(17 - 22*z + 2*z2)*H1m20)/3. +
                (4*opz*(-21 + 11*z)*H100)/3. - (8*(-6 - 5*z + 2*z2)*H110)/3. +
                (4*(15 - 38*z + 27*z2)*H111)/3. + (8*(8 - 12*z + 13*z2)*H112)/3. +
                (8*(23 - 68*z + 57*z2)*H120)/3. - (8*(-5 - 4*z + 6*z2)*H121)/3. +
                (64*(3 + 2*z2)*H2m20)/3. + (16*(2 - 35*z + 4*z2)*H200)/3. -
                (8*(-13 + 19*z2)*H210)/3. - (4*(-1 + 2*z)*(3 + 8*z)*H211)/3. +
                (16*(5 - 6*z + 2*z2)*H212)/3. + (32*(5 - 12*z + 17*z2)*H220)/3. +
                (16*(7 - 14*z + 16*z2)*H221)/3. - (32*(1 + 14*z + 16*z2)*H300)/3. +
                (16*(13 - 2*z + 48*z2)*H310)/3. + (16*(13 - 26*z + 28*z2)*H311)/3. -
                (32*(9 + 10*z + 24*z2)*Hm2m100)/3. + (64*(1 + z + 3*z2)*Hm2000)/3. -
                (64*opz*(9 + 14*z)*Hm1m100)/3. + (8*(27 + 82*z + 61*z2)*Hm1000)/3. -
                (8*(-4 + 21*z + 94*z2)*H0000)/3. - (8*(17 - 56*z + 33*z2)*H1000)/3. +
                (8*(45 - 50*z + 8*z2)*H1100)/3. - (8*(26 - 42*z + 19*z2)*H1110)/3. -
                (4*(1 - 20*z + 16*z2)*H1111)/3. - (16*(11 - 26*z + 40*z2)*H2000)/3. +
                (32*(1 + 2*z + 4*z2)*H2100)/3. - (16*(-5 + 14*z + 4*z2)*H2110)/3. +
                (8*(31 - 62*z + 60*z2)*H2111)/3. - (16*(3 + 4*z + 24*z2)*H00000)/3. +
                ((-241 - 81*z + 1768*z2)*zeta2)/3. + ((-220 - 177*z + 252*z2)*H0*zeta2)/3. +
                (2*(40 + 27*z)*H1*zeta2)/3. - (8*(-4 + 7*z + 7*z2)*H2*zeta2)/3. +
                (64*(1 + 2*z + 7*z2)*H3*zeta2)/3. - (64*z*(2 + 3*z)*Hm20*zeta2)/3. -
                (8*(-13 + 2*z + 21*z2)*Hm10*zeta2)/3. +
                ((17 - 708*z + 656*z2)*H00*zeta2)/3. + (8*(30 - 73*z + 61*z2)*H10*zeta2)/3. -
                (8*(17 - 15*z + 7*z2)*H11*zeta2)/3. + (32*(6 - 10*z + 9*z2)*H20*zeta2)/3. -
                (4*(-1 + 18*z + 30*z2)*H21*zeta2)/3. +
                Lp3*(-4*Lh2*(1 - 2*z + 2*z2) + (13 - 106*z + 114*z2)/3. -
                     4*(1 + 2*z2)*H0 - (4*(11 - 14*z + 6*z2)*H1)/3. +
                     _Lh*(2*(7 - 12*z + 8*z2) + 4*(1 - 2*z + 4*z2)*H0 +
                         8*(1 - 2*z + 2*z2)*H1) - (8*(-1 + 2*z + 2*z2)*H2)/3. +
                     (4*(1 - 2*z + 2*z2)*(6*H00 + 6*H10 - 6*H11))/9. +
                     (2*(-7 + 14*z + 2*z2)*zeta2)/3.) -
                (32*(5 + 6*z + 10*z2)*(2*Hm3m10 + Hm3*zeta2))/3. -
                (8*(21 + 38*z + 36*z2)*(2*Hm2m10 + Hm2*zeta2))/3. -
                (4*(41 + 96*z + 63*z2)*(2*Hm1m10 + Hm1*zeta2))/3. +
                32*opz*(1 + 3*z)*(2*Hm1m1m10 + Hm1m1*zeta2) +
                (4*(317 - 758*z + 757*z2)*zeta3)/3. - 32*opz*(1 + 2*z)*Hm1*zeta3 +
                (4*(57 - 344*z + 398*z2)*H0*zeta3)/3. - (8*omz*(-73 + 175*z)*H1*zeta3)/3. +
                (32*(20 - 34*z + 39*z2)*H2*zeta3)/3. -
                (4*(197 - 122*z + 410*z2)*zeta2*zeta3)/3. +
                Lp2*(8*Lh2*omz*z + (-99 + 150*z - 64*z2)/2. +
                     (-35 + 51*z - 148*z2)*H0 - 2*(26 - 91*z + 74*z2)*H1 -
                     8*(4 - 14*z + 13*z2)*H2 - 24*(1 - 2*z + 4*z2)*H3 +
                     16*(1 + 2*z + 4*z2)*Hm20 + 16*opz*(2 + 3*z)*Hm10 +
                     (-13 - 36*z - 120*z2)*H00 + 8*omz*(-4 + 9*z)*H10 -
                     4*(7 - 30*z + 26*z2)*H11 - 32*(1 - 2*z + 3*z2)*H20 -
                     4*(11 - 22*z + 26*z2)*H21 - 2*(3 + 10*z + 28*z2)*H000 +
                     _Lh*(4*(10 - 33*z + 28*z2) + 4*(1 - 8*z + 8*z2)*H0 - 32*omz*z*H1 +
                         8*(1 - 2*z + 4*z2)*H00 +
                         (4*(1 - 2*z + 2*z2)*(36*H2 + 36*H10 + 36*H11 - 45*zeta2))/9.) +
                     (63 - 92*z + 160*z2)*zeta2 + 2*(13 - 10*z + 52*z2)*H0*zeta2 +
                     (4*(1 + 2*z + 2*z2)*(-72*Hm1m10 + 36*Hm100 - 36*Hm1*zeta2))/9. +
                     (4*(1 - 2*z + 2*z2)*(-108*H12 - 27*H100 - 108*H110 - 117*H111 +
                                          81*H1*zeta2))/9. - 4*(13 - 42*z + 14*z2)*zeta3) +
                (32*(3 + 6*z + 8*z2)*(2*Hm2m1m10 + Hm2m1*zeta2 - Hm2*zeta3))/3. -
                (2*(1 - 2*z + 4*z2)*(52*H5 + 16*H40 - 16*H41 - 55*H000*zeta2 -
                                     172*H00*zeta3))/3. + ((688 + 907*z + 379*z2)*zeta4)/3. +
                (4*(109 - 118*z + 348*z2)*H0*zeta4)/3. +
                _Lp*((-5*(69 - 230*z + 186*z2))/4. + (-11 - 16*z - 372*z2)*H0 -
                    4*(-6 - 77*z + 93*z2)*H1 - 4*(11 - 94*z + 66*z2)*H2 -
                    4*(3 - 44*z + 32*z2)*H3 - 56*(1 - 2*z + 4*z2)*H4 +
                    32*(3 + 2*z + 8*z2)*Hm30 + 32*(6 + 7*z + 7*z2)*Hm20 +
                    16*(21 + 32*z + 13*z2)*Hm10 - 2*(14 + 123*z + 192*z2)*H00 -
                    16*(6 - 14*z + 11*z2)*H10 - 8*(9 - 35*z + 33*z2)*H11 +
                    128*omz*z*H12 - 8*(3 + 4*z2)*H20 - 8*z*(-11 + 9*z)*H21 -
                    64*(1 - z + 3*z2)*H30 - 8*(9 - 18*z + 20*z2)*H31 -
                    32*(3 + 2*z + 4*z2)*Hm2m10 + 16*(5 + 6*z + 12*z2)*Hm200 -
                    32*opz*(6 + 7*z)*Hm1m10 + 16*opz*(10 + 13*z)*Hm100 -
                    2*(1 + 60*z + 156*z2)*H000 + 8*omz*(2 + 13*z)*H100 +
                    32*omz*(2 + z)*H110 + 8*omz*(2 + 9*z)*H111 -
                    8*(3 - 14*z + 20*z2)*H200 + 16*(-1 + 2*z + 2*z2)*H210 -
                    8*(5 - 10*z + 8*z2)*H211 - 4*(5 + 6*z + 40*z2)*H0000 +
                    6*(16 - 23*z + 88*z2)*zeta2 - 16*(3 + 2*z + 4*z2)*Hm2*zeta2 -
                    16*opz*(6 + 7*z)*Hm1*zeta2 + 2*(7 - 64*z + 120*z2)*H0*zeta2 -
                    32*omz*(1 + 4*z)*H1*zeta2 + 8*(5 - 18*z + 22*z2)*H2*zeta2 +
                    60*(1 - 2*z + 4*z2)*H00*zeta2 +
                    _Lh*(2*(13 - 75*z + 72*z2) + 2*(-8 - 15*z + 8*z2)*H0 + 4*z*(-3 + 4*z)*H1 -
                        16*omz*z*H2 + 2*(-1 - 12*z + 8*z2)*H00 - 16*omz*z*H10 -
                        16*omz*z*H11 + 4*(1 - 2*z + 4*z2)*H000 + 24*omz*z*zeta2 +
                        (4*(1 - 2*z + 2*z2)*(-18*H21 + 18*H100 - 18*H111 - 126*zeta3))/9.) +
                    8*(19 - 30*z + 74*z2)*zeta3 + 16*(7 - 18*z + 28*z2)*H0*zeta3 +
                    (4*(1 + 2*z + 2*z2)*(-216*Hm1m20 + 72*Hm120 + 144*Hm1m1m10 -
                                         216*Hm1m100 + 72*Hm1000 + 72*Hm1m1*zeta2 - 36*Hm1*zeta3))/9. +
                    (4*(1 - 2*z + 2*z2)*(-252*H13 - 108*H22 + 72*H1m20 - 108*H112 -
                                         144*H120 - 180*H121 - 108*H1000 - 180*H1100 + 36*H1110 -
                                         72*H1111 + 270*H10*zeta2 + 198*H11*zeta2 + 468*H1*zeta3))/9. +
                    ((-153 + 370*z - 10*z2)*zeta4)/2.) -
                (16*(1 + 2*z + 2*z2)*(22*Hm1m30 - 6*Hm130 - 20*Hm1m2m10 + 30*Hm1m200 -
                                      16*Hm1m1m20 - 22*Hm1200 + 12*Hm1210 + 16*Hm1m1m1m10 -
                                      24*Hm1m1m100 + 6*Hm1m1000 - 5*Hm10000 - 10*Hm1m2*zeta2 +
                                      12*Hm12*zeta2 + 8*Hm1m1m1*zeta2 - 6*Hm1m10*zeta2 - 8*Hm1m1*zeta3 +
                                      4*Hm10*zeta3 + 22*Hm1*zeta4))/3. -
                4*(1 - 2*z + 2*z2)*(16*H14 - (40*H1m30)/3. + (8*H113)/3. - 8*H122 -
                                    (8*H130)/3. - 4*H131 + (16*H1m200)/3. - (16*H11m20)/3. - (4*H1112)/3. -
                                    (68*H1120)/3. - (28*H1121)/3. - (20*H1200)/3. - (52*H1210)/3. -
                                    (56*H1211)/3. + (28*H10000)/3. + (80*H11000)/3. - 4*H11100 +
                                    (4*H11110)/3. - 20*H11111 - (8*H12*zeta2)/3. - 17*H100*zeta2 -
                                    (32*H110*zeta2)/3. + 5*H111*zeta2 - 52*H10*zeta3 - 52*H11*zeta3 -
                                    (100*H1*zeta4)/3.) - (8*(63 - 454*z + 84*z2)*zeta5)/3.) +
        CA2*TR*((2*(-470494 + 450714*z - 2175249*z2 + 2153647*z3))/(729.*z) +
                (-2*(17152 + 167620*z + 62809*z2 + 747175*z3)*H0 -
                 2*(-6528 - 50692*z + 24401*z2 + 30151*z3)*H1)/(243.*z) +
                (4*(-321 - 483*z + 1460*z2)*H3)/27. - (8*(51 - 48*z + 137*z2)*H4)/9. -
                (4*(9 - 58*z + 72*z2)*H5)/3. + (8*(-41 + 74*z + 8*z2)*Hm40)/3. +
                (8*(216 + 180*z + 817*z2)*Hm30)/9. + (8*(29 + 46*z + 56*z2)*Hm32)/3. -
                (4*(-213 - 30*z + 2426*z2)*Hm20)/27. + (8*(31 + 37*z + 96*z2)*Hm22)/3. +
                (4*(49 + 110*z + 100*z2)*Hm23)/3. + (2*(20278 + 29235*z + 126201*z2)*H00)/81. -
                (4*(29 + 66*z + 28*z2)*H23)/3. + (4*(51 - 135*z + 317*z2)*H31)/9. -
                (8*(19 + 50*z + 20*z2)*H32)/3. - 192*z*H40 + 4*(-3 - 14*z + 8*z2)*H41 -
                (16*(3 + 22*z + 10*z2)*Hm3m10)/3. + (4*(-35 + 214*z + 52*z2)*Hm300)/3. +
                (8*(13 - 26*z + 8*z2)*Hm2m20)/3. - (8*(14 + 37*z + 78*z2)*Hm2m10)/3. -
                (8*(39 + 66*z + 76*z2)*Hm2m12)/3. + (4*(297 + 246*z + 1522*z2)*Hm200)/9. -
                (8*(3 + 18*z + 8*z2)*Hm220)/3. + (56*opz*Hm120)/3. +
                (opz*(12*(8 + 9*z + 60*z2)*Hm13 - 8*(16 + 53*z + 112*z2)*Hm1m12 -
                      32*(2 + 4*z + 17*z2)*Hm121))/(9.*z) + 4*(3 + 6*z + 4*z2)*Hm130 -
                (4*(1679 - 1712*z + 7385*z2)*H000)/27. + (16*(7 + 10*z + 8*z2)*H2m20)/3. -
                (4*(9 - 62*z + 16*z2)*H212)/3. + (8*(-1 - 62*z + 24*z2)*H220)/3. -
                (8*(1 - 26*z + 2*z2)*H221)/3. + (8*(17 - 100*z + 6*z2)*H300)/3. +
                (4*(-25 - 86*z + 4*z2)*H310)/3. + (16*(3 + 8*z)*H311)/3. -
                32*(1 + z2)*Hm2m1m10 - (4*(5 + 254*z + 80*z2)*Hm2m100)/3. +
                (8*(-7 + 78*z + 10*z2)*Hm2000)/3. + (8*(-1 - 2*z + 4*z2)*Hm1210)/3. -
                (8*(-55 - 345*z + 108*z2)*H0000)/9. - (4*(-11 - 118*z + 52*z2)*H2000)/3. -
                (4*(43 + 198*z + 20*z2)*H2100)/3. + (4*(27 + 94*z + 20*z2)*H2110)/3. +
                (16*(1 + 4*z)*H2111)/3. + 16*(-4 + 15*z)*H00000 -
                (16*(16 + 34*z + 33*z2)*Hm3*zeta2)/3. -
                (4*(76 + 111*z + 270*z2)*Hm2*zeta2)/3. + (4*(-13 + 18*z + 20*z2)*H3*zeta2)/3. +
                (8*(33 + 66*z + 70*z2)*Hm2m1*zeta2)/3. -
                (8*(17 + 90*z + 50*z2)*Hm20*zeta2)/3. + (4*(47 - 2*z + 88*z2)*H20*zeta2)/3. +
                (32*(6 - z + 8*z2)*H21*zeta2)/3. - 4*(7 + 14*z + 16*z2)*Hm100*zeta2 +
                12*(-1 + 2*z)*(-1 + 4*z)*H000*zeta2 +
                Lp3*((2*omz*(272 - 97*z + 1883*z2))/(27.*z) -
                     (4*(-8 - 11*z - 152*z2 + 18*z3)*H0)/(9.*z) -
                     (4*(-16 - z - 118*z2 + 146*z3)*H1)/(9.*z) + (32*(1 + 4*z)*H2)/3. +
                     (16*(-1 + 8*z)*H00)/3. + (4*(1 - 2*z + 2*z2)*(-12*H10 - 24*H11))/9. -
                     (32*(1 + 4*z)*zeta2)/3.) + (-4*(-1256 - 702*z - 12231*z2 + 98*z3)*H2 +
                                                 2*(5436 + 17579*z + 33484*z2 + 23245*z3)*Hm10 -
                                                 2*(-9490 + 9435*z - 43401*z2 + 43978*z3)*H10 +
                                                 2*(-918 + 7796*z - 12541*z2 + 6748*z3)*H11 -
                                                 2*(-12536 + 11376*z - 58657*z2 + 49619*z3)*zeta2)/(81.*z) +
                (-8*(116 + 96*z + 393*z2 + 422*z3)*Hm12 -
                 4*(-208 + 1463*z - 1678*z2 + 92*z3)*H12 -
                 4*(-116 - 111*z - 861*z2 + 445*z3)*H21 +
                 8*(40 + 54*z + 276*z2 + 7*z3)*Hm1m10 -
                 4*(484 + 124*z + 431*z2 + 393*z3)*Hm100 -
                 4*(40 + 175*z - 194*z2 + 337*z3)*H111 +
                 4*(272 + 246*z + 1062*z2 + 851*z3)*Hm1*zeta2)/(27.*z) +
                8*(5 + 10*z + 8*z2)*(Hm1200 - Hm12*zeta2) +
                (36*(-8 - 89*z - 138*z2 - 37*z3 + 19*z4)*H30 +
                 4*(84 + 393*z + 1024*z2 + 1133*z3 + 436*z4)*H210 +
                 4*(30 - 321*z - 323*z2 + 473*z3 + 436*z4)*H00*zeta2)/(9.*opz2)\
                - (16*(8 + 22*z + 19*z2)*Hm2*zeta3)/3. + 4*(27 + 30*z + 52*z2)*H2*zeta3 -
                (16*(14 + 28*z + 25*z2)*Hm10*zeta3)/3. - (16*(27 + 41*z)*H00*zeta3)/3. -
                (4*(93 + 242*z + 136*z2)*zeta2*zeta3)/3. +
                Lp2*((4*(-103 + 77*z - 534*z2 + 554*z3))/(3.*z) -
                     (4*(52 + 302*z + 146*z2 + 1407*z3)*H0)/(9.*z) +
                     (4*(-52 + 2*z - 178*z2 + 143*z3)*H1)/(9.*z) -
                     (8*(-4 - 3*z - 24*z2 + 44*z3)*H2)/(3.*z) + 32*(1 + 3*z)*H3 +
                     32*(-1 + 4*z)*Hm20 + (4*(16 - z + 118*z2 + 146*z3)*Hm10)/(3.*z) -
                     (16*(-3 - 23*z + 9*z2)*H00)/3. + (8*omz*(4 + 7*z + 31*z2)*H10)/(3.*z) +
                     (4*(-8 + 5*z - 70*z2 + 84*z3)*H11)/(3.*z) + 16*(1 + 4*z)*H20 -
                     16*(1 + 4*z)*H21 + 32*(-1 + 4*z)*H000 +
                     (8*(4 - 3*z + 35*z2 + 44*z3)*zeta2)/(3.*z) - 32*omz*H0*zeta2 +
                     (4*(1 + 2*z + 2*z2)*(36*Hm12 - 72*Hm1m10 + 54*Hm100 - 72*Hm1*zeta2))/9. +
                     (4*(1 - 2*z + 2*z2)*(-18*H100 + 108*H111 - 36*H1*zeta2))/9. +
                     16*(-2 + 17*z)*zeta3) + (6*(208 + 556*z + 1118*z2 + 2811*z3 + 2059*z4)*
                                              H20 + 4*(-1494 - 235*z - 6050*z2 + 123*z3 + 7486*z4)*H100 -
                                              2*(-1984 + 805*z - 6983*z2 - 1361*z3 + 8519*z4)*H110 +
                                              2*(624 + 2238*z + 5946*z2 + 7511*z3 + 3233*z4)*H0*zeta2 -
                                              2*(-1544 - 1467*z - 6989*z2 + 1403*z3 + 8577*z4)*H1*zeta2 -
                                              2*(-3288 + 2178*z - 8142*z2 - 13207*z3 + 509*z4)*zeta3)/(27.*z*opz) +
                (4*(-112 + 217*z - 761*z2 + 624*z3)*H13 +
                 12*(-16 - 59*z - 172*z2 + 97*z3)*H22 -
                 8*(16 + 39*z + 195*z2 + 205*z3)*Hm1m20 -
                 16*(-8 - 60*z + 69*z2 + 8*z3)*H1m20 -
                 4*(-64 + 220*z - 539*z2 + 471*z3)*H112 +
                 4*(-48 - 25*z - 352*z2 + 457*z3)*H120 -
                 4*(-48 + 124*z - 404*z2 + 383*z3)*H121 -
                 16*(-2 + 6*z2 + 53*z3)*H211 +
                 8*(8 + 51*z + 84*z2 + 74*z3)*Hm1m1m10 -
                 12*(24 + 16*z + 245*z2 + 286*z3)*Hm1m100 +
                 4*(48 + z + 494*z2 + 647*z3)*Hm1000 -
                 4*(-8 - 147*z - 63*z2 + 233*z3)*H1000 +
                 4*(-112 - 4*z - 685*z2 + 845*z3)*H1100 -
                 4*(-112 + 233*z - 811*z2 + 734*z3)*H1110 -
                 4*(-8 - 61*z + 14*z2)*H1111 +
                 4*(40 + 189*z + 414*z2 + 298*z3)*Hm1m1*zeta2 -
                 4*(56 - 30*z + 264*z2 + 401*z3)*Hm10*zeta2 -
                 4*(-96 + 97*z - 389*z2 + 287*z3)*H10*zeta2 -
                 8*(-20 - 5*z - 98*z2 + 112*z3)*H11*zeta2 -
                 12*(8 + 33*z + 72*z2 + 58*z3)*Hm1*zeta3 -
                 4*(-176 + 38*z - 307*z2 + 555*z3)*H1*zeta3)/(9.*z) +
                (8*(1 + 2*z + 2*z2)*(10*Hm14 - 16*Hm221 - 24*Hm1m30 - 24*Hm1m22 -
                                     11*Hm1m13 + 10*Hm122 - 10*Hm131 + 8*Hm1m2m10 - 16*Hm1m200 +
                                     2*Hm1m1m20 + 18*Hm1m1m12 + 16*Hm1m121 + 4*Hm1211 + 7*Hm1m1m100 -
                                     3*Hm1m1000 + 13*Hm10000 + 28*Hm1m2*zeta2 - 18*Hm1m1m1*zeta2 +
                                     18*Hm1m10*zeta2 + 4*Hm1m1*zeta3))/3. -
                (4*(83 + 166*z + 163*z2)*Hm1*zeta4)/3. + ((-425 - 38*z + 136*z2)*H0*zeta4)/3. +
                _Lp*((-2*(-28300 + 32664*z - 196185*z2 + 194638*z3))/(81.*z) +
                    (2*(1136 + 6994*z + 20623*z2 + 45240*z3)*H0)/(27.*z) -
                    (8*(1282 - 1691*z + 557*z2)*H1)/27. +
                    (4*(-104 - 291*z - 834*z2 + 197*z3)*H2)/(9.*z) -
                    4*(13 - 17*z + 44*z2)*H3 - 4*(-1 - 26*z + 24*z2)*H4 + 8*(-13 + 18*z)*Hm30 +
                    (8*(44 + 13*z + 234*z2)*Hm20)/3. -
                    (8*(122 - 31*z + 388*z2 + 447*z3)*Hm10)/(9.*z) +
                    (32*opz*(2 + 4*z + 11*z2)*Hm12)/(3.*z) -
                    (2*(1543 - 252*z + 6666*z2)*H00)/9. -
                    (4*omz*(96 + 31*z + 513*z2)*H10)/(3.*z) -
                    (4*(-78 + z - 218*z2 + 125*z3)*H11)/(9.*z) +
                    (8*(-12 + 20*z - 94*z2 + 97*z3)*H12)/(3.*z) +
                    (16*(-4 - 9*z - 56*z2 + 7*z3)*H20)/(3.*z) +
                    (8*(-4 - 18*z2 + 97*z3)*H21)/(3.*z) - 16*(1 + 8*z)*H22 - 128*z*H30 -
                    16*(3 + 8*z)*H31 - 16*(3 + 2*z + 4*z2)*Hm2m10 +
                    8*(-1 + 42*z + 12*z2)*Hm200 - (8*(47 + 76*z + 40*z2)*Hm1m10)/3. +
                    (4*(40 + 69*z + 372*z2 + 376*z3)*Hm100)/(3.*z) -
                    (4*(-37 - 176*z + 90*z2)*H000)/3. -
                    (4*(-8 - 100*z + 35*z2 + 84*z3)*H100)/(3.*z) -
                    (8*omz*(20 + 17*z + 137*z2)*H110)/(3.*z) -
                    (88*(1 - 2*z + 2*z2)*H111)/3. + 4*(21 + 34*z + 16*z2)*H200 -
                    16*(3 + 16*z)*H210 + 16*(-5 + 19*z)*H0000 +
                    (4*(-376 + 489*z - 1274*z2 + 1173*z3)*zeta2)/(9.*z) -
                    8*(11 + 18*z + 20*z2)*Hm2*zeta2 -
                    (4*(16 + 95*z + 196*z2 + 128*z3)*Hm1*zeta2)/(3.*z) +
                    (4*(-16 + 15*z - 261*z2 + 168*z3)*H0*zeta2)/(3.*z) +
                    (4*(-8 - 53*z + 28*z2 + 22*z3)*H1*zeta2)/(3.*z) -
                    8*(5 + 6*z + 4*z2)*H2*zeta2 + 4*(-1 - 22*z + 24*z2)*H00*zeta2 +
                    (8*(-44 - 21*z - 364*z2 + 408*z3)*zeta3)/(3.*z) - 16*(13 + 34*z)*H0*zeta3 +
                    (4*(1 + 2*z + 2*z2)*(144*Hm22 + 108*Hm13 - 36*Hm1m20 - 144*Hm1m12 -
                                         72*Hm121 - 72*Hm1m1m10 - 180*Hm1m100 + 90*Hm1000 +
                                         108*Hm1m1*zeta2 - 72*Hm10*zeta2 - 36*Hm1*zeta3))/9. +
                    (4*(1 - 2*z + 2*z2)*(-72*H13 + 108*H1m20 + 108*H112 + 144*H120 +
                                         36*H121 - 90*H1000 + 72*H1100 + 180*H1110 - 216*H1111 +
                                         180*H10*zeta2 + 72*H11*zeta2 + 540*H1*zeta3))/9. +
                    (-361 - 1294*z - 112*z2)*zeta4) +
                (-4*(144 + 426*z + 1830*z2 + 3167*z3 + 1756*z4 + 155*z5)*H200 +
                 12*(16 + 103*z + 398*z2 + 671*z3 + 456*z4 + 102*z5)*H2*zeta2 +
                 4*(48 - 42*z - 1470*z2 - 811*z3 + 2308*z4 + 1757*z5)*H0*zeta3 +
                 (1840 + 5012*z + 17187*z2 + 31855*z3 + 22925*z4 + 5103*z5)*zeta4)/
                (9.*z*opz2) - (2*(1 - 2*z + 2*z2)*
                               (4*H14 - 56*H1m30 + 16*H1m22 - 64*H113 - 28*H122 - 76*H130 -
                                20*H131 - 16*H1m2m10 + 140*H1112 - 8*H1120 + 80*H1121 -
                                116*H1200 - 4*H1210 + 40*H1211 + 44*H10000 + 76*H11000 -
                                52*H11100 + 172*H11110 - 120*H11111 - 24*H1m2*zeta2 -
                                4*H12*zeta2 - 44*H100*zeta2 + 32*H110*zeta2 + 16*H111*zeta2 -
                                120*H10*zeta3 + 172*H11*zeta3 + 57*H1*zeta4))/3. -
                (4*(77 - 1206*z + 136*z2)*zeta5)/3.) +
        CA*CF*TR*((-116208 - 834293*z + 59486*z2 + 803580*z3)/(1944.*z) +
                  ((-13165 - 32212*z - 26970*z2)*H0)/243. +
                  (2*(-2285 - 60229*z + 15479*z2 + 48075*z3)*H1)/(243.*z) +
                  Lp4*((2*(-12 + 13*z - 116*z2 + 137*z3))/(9.*z) - 4*(1 + 4*z)*H0 +
                       4*(1 - 2*z + 2*z2)*H1) - (8*(2576 + 3701*z + 2672*z2)*H2)/81. +
                  (4*(134 - 278*z + 179*z2)*H3)/9. - (4*(-6 + 78*z + 631*z2)*H4)/9. +
                  (16*(11 + 16*z + 8*z2)*H5)/3. - (16*(25 - 14*z + 28*z2)*Hm40)/3. -
                  (8*(139 + 84*z + 63*z2)*Hm30)/3. + (16*(23 + 10*z + 44*z2)*Hm32)/3. +
                  (4*(-339 + 194*z + 70*z2)*Hm20)/3. + 8*(15 + 28*z + 16*z2)*Hm22 +
                  (16*(23 + 30*z + 47*z2)*Hm23)/3. + (2*(-929 - 824*z + 125*z2)*Hm10)/3. +
                  (8*(55 + 82*z + 28*z2)*Hm12)/3. + 8*opz*(10 + 13*z)*Hm13 +
                  ((1745 + 26249*z - 47946*z2)*H00)/81. + (4*(-2 - 881*z + 694*z2)*H21)/27. -
                  (8*(-18 - 7*z + 14*z2)*H22)/3. + (16*(9 - 40*z + 3*z2)*H23)/3. +
                  (4*(3 + 186*z + 172*z2)*H31)/9. + (16*(1 - 14*z + 4*z2)*H32)/3. +
                  (64*omz*H40)/3. + (8*(-3 - 10*z + 4*z2)*H41)/3. -
                  (16*(1 - 10*z + 10*z2)*Hm3m10)/3. + (16*(4 + 8*z + 25*z2)*Hm300)/3. +
                  (16*(3 - 22*z + 8*z2)*Hm2m20)/3. + 8*(6 + 20*z + z2)*Hm2m10 -
                  (32*(17 + 14*z + 33*z2)*Hm2m12)/3. - (4*(79 + 45*z2)*Hm200)/3. +
                  (16*(3 + 22*z + 4*z2)*Hm220)/3. + (128*z*Hm221)/3. +
                  (8*(53 + 92*z + 27*z2)*Hm1m20)/3. - (8*(36 + 64*z + 35*z2)*Hm1m10)/3. -
                  16*opz*(11 + 8*z)*Hm1m12 + (4*(-133 - 98*z + 27*z2)*Hm100)/3. +
                  (8*(-14 - 18*z + 5*z2)*Hm120)/3. + (16*(-7 + z)*opz*Hm121)/3. -
                  16*(1 + 2*z + 3*z2)*Hm130 + ((299 + 11222*z + 1836*z2)*H000)/27. -
                  (32*(22 - 34*z + 9*z2)*H1m20)/3. +
                  (4*(-16 + 1056*z - 1290*z2 + 349*z3)*H12 +
                   8*(46 + 20*z + 74*z2 + 21*z3)*H111)/(27.*z) -
                  (32*(5 - 4*z + 2*z2)*H2m20)/3. + (4*(-37 + 128*z + 176*z2)*H211)/9. +
                  (8*(35 - 34*z + 76*z2)*H212)/3. - (32*(1 + 7*z + 2*z2)*H220)/3. +
                  (8*(11 - 46*z + 18*z2)*H221)/3. - (8*(17 - 22*z + 10*z2)*H300)/3. -
                  (16*(9 + 4*z + 31*z2)*H310)/3. - (32*(2 + 3*z + 2*z2)*H311)/3. -
                  (16*(5 + 2*z + 18*z2)*Hm2m1m10)/3. - (8*(27 + 70*z + 64*z2)*Hm2m100)/3. +
                  (8*(37 + 70*z + 90*z2)*Hm2000)/3. - 8*(10 + 16*z + z2)*Hm1m1m10 +
                  (4*(8 + 5*z)*(10 + 9*z)*Hm1m100)/3. - (4*(-9 - 26*z + 13*z2)*Hm1000)/3. -
                  (32*(2 + 4*z + 7*z2)*Hm1200)/3. + 32*(2 + 4*z + 5*z2)*Hm1210 -
                  (4*(121 - 230*z + 948*z2)*H0000)/9. + (16*(-2 - 24*z + 7*z2)*H2000)/3. -
                  (16*(-3 - 14*z + 7*z2)*H2100)/3. + (16*(10 - 14*z + 47*z2)*H2110)/3. -
                  8*(11 - 18*z + 20*z2)*H2111 + (32*(3 + 8*z)*H00000)/3. -
                  (8*(47 + 10*z + 98*z2)*Hm3*zeta2)/3. - 4*(24 + 36*z + 31*z2)*Hm2*zeta2 -
                  (4*(146 + 228*z + 91*z2)*Hm1*zeta2)/3. - 16*(5 - 10*z + 17*z2)*H3*zeta2 +
                  8*(21 + 18*z + 38*z2)*Hm2m1*zeta2 - (16*(31 + 30*z + 54*z2)*Hm20*zeta2)/3. +
                  4*(34 + 60*z + 31*z2)*Hm1m1*zeta2 - (4*(143 + 260*z + 99*z2)*Hm10*zeta2)/3. +
                  16*(3 + 6*z + 8*z2)*Hm12*zeta2 - (8*(25 - 92*z + 8*z2)*H20*zeta2)/3. -
                  (32*(6 - 11*z + z2)*H21*zeta2)/3. -
                  (4*(49 + 98*z + 110*z2)*Hm100*zeta2)/3. -
                  (4*(47 + 62*z + 32*z2)*H000*zeta2)/3. +
                  Lp3*((176 - 25*z + 910*z2 - 924*z3)/(9.*z) +
                       (4*(25 + 148*z + 26*z2)*H0)/9. -
                       (8*(-4 + 23*z - 79*z2 + 80*z3)*H1)/(9.*z) +
                       _Lh*((2*(-8 + 5*z - 70*z2 + 84*z3))/(3.*z) - 8*(1 + 4*z)*H0 +
                           8*(1 - 2*z + 2*z2)*H1) - (8*(-1 - 10*z + 4*z2)*H2)/3. -
                       8*(1 + 2*z + 2*z2)*Hm10 - (16*opz*H00)/3. +
                       (4*(1 - 2*z + 2*z2)*(-12*H10 - 42*H11))/9. + (4*(-5 - 26*z + 2*z2)*zeta2)/3.)
                  + (-2*(-404 - 1007*z - 5051*z2 + 6190*z3)*H10 +
                     4*(661 - 1894*z + 905*z2 + 564*z3)*H11 +
                     2*(-1812 + 13193*z - 10027*z2 + 13814*z3)*zeta2)/(81.*z) +
                  (2*(708 + 180*z + 1765*z2 + 2239*z3)*H20 -
                   6*(310 - 1272*z - 1097*z2 + 503*z3)*H0*zeta2)/(27.*opz) +
                  (omz*(-96*(2 + 17*z2)*H121 + 4*(16 + 49*z + 25*z2)*H11*zeta2))/(9.*z) -
                  (16*(23 + 18*z + 42*z2)*Hm2*zeta3)/3. - (8*(32 + 66*z + 37*z2)*Hm1*zeta3)/3. +
                  (8*(-31 - 38*z + 20*z2)*H2*zeta3)/3. - (8*(41 + 82*z + 70*z2)*Hm10*zeta3)/3. +
                  (8*(-33 - 66*z + 4*z2)*H00*zeta3)/3. + (8*(77 - 9*z + 128*z2)*zeta2*zeta3)/3. +
                  _Lh*((-1616*omz*z)/27. + 56*omz*z*zeta3) +
                  Lp2*((-1984 + 189*z + 7644*z2 - 7100*z3)/(54.*z) +
                       ((7 + 1894*z - 2088*z2)*H0)/9. -
                       (2*(104 + 110*z - 313*z2 + 368*z3)*H1)/(9.*z) -
                       (8*(13 - 14*z + 86*z2)*H2)/3. - 8*(-3 - 14*z + 4*z2)*H3 +
                       4*(3 + 4*z + 4*z2)*Hm10 - (2*(-5 - 176*z + 308*z2)*H00)/3. -
                       (8*(-6 + 11*z - 58*z2 + 64*z3)*H10)/(3.*z) -
                       (8*(-4 - 8*z - 14*z2 + 27*z3)*H11)/(3.*z) - 16*z*(-7 + 2*z)*H20 -
                       8*(1 - 14*z + 4*z2)*H21 + 4*(3 + 14*z)*H000 +
                       _Lh*((-4*(-52 + 5*z - 103*z2 + 65*z3))/(9.*z) + (8*(3 + 44*z2)*H0)/3. -
                           16*(1 + 2*z)*H00 - 16*(1 - 2*z + 2*z2)*H11 +
                           (4*(1 + 2*z + 2*z2)*(-36*Hm10 - 18*zeta2))/9.) +
                       (2*(4 + 15*z + 48*z2 + 275*z3)*zeta2)/(3.*z) +
                       4*(-7 - 32*z + 8*z2)*H0*zeta2 +
                       (4*(1 + 2*z + 2*z2)*(-18*Hm20 + 36*Hm12 + 18*Hm1m10 + 9*Hm100 -
                                            27*Hm1*zeta2))/9. + (4*(1 - 2*z + 2*z2)*
                                                                 (-45*H12 - 54*H100 - 45*H110 + 9*H111 + 45*H1*zeta2))/9. +
                       8*(1 - 9*z + 6*z2)*zeta3) +
                  (-4*(92 + 736*z - 575*z2 - 563*z3 + 710*z4)*H100 +
                   2*(-224 + 655*z - 1863*z2 - 763*z3 + 2087*z4)*H110 -
                   2*(-32 + 1999*z + 432*z2 - 1324*z3 + 167*z4)*H1*zeta2 -
                   2*(-1536 - 5458*z - 15548*z2 - 4303*z3 + 7215*z4)*zeta3)/(27.*z*opz) +
                  (-4*(84 - 546*z - 758*z2 + 494*z3 + 613*z4)*H30 -
                   12*(72 + 40*z + 45*z2 + 282*z3 + 199*z4)*H200 +
                   12*(-11 - 48*z - 24*z2 + 76*z3 + 57*z4)*H210 +
                   12*(-39 - 104*z - 72*z2 + 36*z3 + 37*z4)*H2*zeta2 +
                   2*(-3 + 438*z + 1847*z2 + 2296*z3 + 908*z4)*H00*zeta2 +
                   4*(255 + 2046*z + 4405*z2 + 3764*z3 + 1132*z4)*H0*zeta3)/
                  (9.*opz2) + (4*(-32 - 81*z - 72*z2 + 329*z3)*H13 +
                               8*(-8 + 69*z - 87*z2 + 17*z3)*H112 -
                               4*(16 + 114*z - 204*z2 + 173*z3)*H120 -
                               4*(8 + 239*z - 340*z2 + 203*z3)*H1000 -
                               4*(-16 + 174*z - 522*z2 + 373*z3)*H1100 +
                               4*(-16 + 273*z - 552*z2 + 277*z3)*H1110 +
                               8*(-4 - 29*z - 23*z2 + 24*z3)*H1111 -
                               4*(-44 + 177*z - 564*z2 + 683*z3)*H10*zeta2 -
                               24*(8 + 97*z - 187*z2 + 88*z3)*H1*zeta3)/(9.*z) +
                  (8*(1 + 2*z + 2*z2)*(20*Hm14 + 32*Hm1m30 - 60*Hm1m22 - 53*Hm1m13 +
                                       10*Hm122 - 6*Hm131 - 24*Hm1m2m10 - 22*Hm1m200 - 10*Hm1m1m20 +
                                       78*Hm1m1m12 - 4*Hm1m120 - 4*Hm1211 + 32*Hm1m1m1m10 +
                                       41*Hm1m1m100 - 51*Hm1m1000 + 15*Hm10000 + 48*Hm1m2*zeta2 -
                                       62*Hm1m1m1*zeta2 + 61*Hm1m10*zeta2 + 44*Hm1m1*zeta3))/3. +
                  ((312 - 2291*z - 10332*z2 - 15038*z3 - 9732*z4 - 2441*z5)*zeta4)/
                  (9.*z*opz2) + (8*(59 + 118*z + 121*z2)*Hm1*zeta4)/3. -
                  (2*(167 - 130*z + 148*z2)*H0*zeta4)/3. +
                  _Lp*((-864 - 89269*z + 35990*z2 + 24910*z3)/(324.*z) +
                      ((-1377 - 3324*z - 18776*z2)*H0)/27. -
                      (4*(-344 - 177*z - 294*z2 + 1022*z3)*H1)/(27.*z) +
                      (4*(244 + 124*z + 585*z2)*H2)/9. - (4*(27 + 126*z + 388*z2)*H3)/3. +
                      8*(11 + 18*z + 4*z2)*H4 - 64*(1 + z2)*Hm30 -
                      16*(13 - 3*z + 5*z2)*Hm20 + 32*(3 + 2*z + 6*z2)*Hm22 +
                      24*opz*(-14 + 5*z)*Hm10 + 16*opz*(7 + 5*z)*Hm12 -
                      (2*(-349 - 2704*z + 30*z2)*H00)/9. +
                      (4*(-140 + 61*z - 164*z2 + 91*z3)*H10)/(9.*z) +
                      (4*(-104 + 163*z - 848*z2 + 745*z3)*H11)/(9.*z) - 8*omz*(-2 + 7*z)*H12 -
                      8*(4 - 5*z + 56*z2)*H20 - (8*(-8 + 37*z + 54*z2)*H21)/3. +
                      8*(1 - 2*z + 4*z2)*H22 + 32*(2 + 7*z)*H30 + 16*(3 + 2*z + 4*z2)*H31 -
                      16*(1 + 2*z)*(1 + 4*z)*Hm2m10 + 8*(7 + 10*z + 20*z2)*Hm200 +
                      8*(19 + 32*z + 10*z2)*Hm1m10 - 4*(11 + 20*z + 12*z2)*Hm100 -
                      (4*(31 - 86*z + 404*z2)*H000)/3. - (8*(67 - 119*z + 65*z2)*H100)/3. -
                      (8*(-8 + 3*z - 63*z2 + 59*z3)*H110)/(3.*z) +
                      (8*(5 - 43*z + 49*z2)*H111)/3. - 24*(3 - 2*z + 4*z2)*H200 -
                      8*(-1 - 22*z + 12*z2)*H210 + 8*(11 - 22*z + 20*z2)*H211 +
                      8*(5 + 14*z)*H0000 - ((40 + 1007*z + 1862*z2 + 2914*z3)*zeta2)/(9.*z) -
                      8*(13 + 14*z + 32*z2)*Hm2*zeta2 - 4*(9 + 16*z + 10*z2)*Hm1*zeta2 +
                      8*(5 + 31*z + 66*z2)*H0*zeta2 - 4*(-13 + 6*z + 4*z2)*H1*zeta2 -
                      32*(1 - 3*z + 5*z2)*H2*zeta2 - 32*opz*(3 + z)*H00*zeta2 +
                      (2*(-96 + 539*z - 1294*z2 + 706*z3)*zeta3)/(9.*z) + 16*(-1 + 11*z)*H0*zeta3 +
                      _Lh*((4*(-172 + 517*z - 1193*z2 + 1104*z3))/(27.*z) -
                          (8*(21 - 30*z + 68*z2)*H0)/9. - 4*z*(-3 + 4*z)*H1 - 16*z*opz*Hm10 +
                          (4*(3 - 12*z + 44*z2)*H00)/3. + (16*omz*(2 - z + 11*z2)*H10)/(3.*z) +
                          16*omz*z*H11 + 32*z*H20 - 8*(1 + 2*z)*H000 +
                          (4*(1 - 2*z + 2*z2)*(-18*H12 - 18*H110 + 18*H111))/9. -
                          (16*(-2 + 3*z - 12*z2 + 14*z3)*zeta2)/(3.*z) +
                          (4*(1 + 2*z + 2*z2)*(-36*Hm20 + 36*Hm1m10 - 18*Hm100 + 18*Hm1*zeta2))/
                          9. - 4*(7 - 10*z + 14*z2)*zeta3) +
                      (4*(1 + 2*z + 2*z2)*(216*Hm13 - 36*Hm1m20 - 216*Hm1m12 + 36*Hm120 +
                                           72*Hm121 + 144*Hm1m1m10 - 180*Hm1m100 + 180*Hm1000 +
                                           288*Hm1m1*zeta2 - 342*Hm10*zeta2 - 288*Hm1*zeta3))/9. +
                      (4*(1 - 2*z + 2*z2)*(180*H13 - 108*H1m20 + 108*H112 - 36*H120 +
                                           252*H121 - 72*H1100 - 108*H1110 + 288*H1111 - 288*H10*zeta2 -
                                           270*H11*zeta2 - 396*H1*zeta3))/9. + 8*opz*(16 + 7*z)*zeta4) +
                  (2*(1 - 2*z + 2*z2)*(44*H14 - 8*H1m30 - 16*H1m22 + 24*H113 + 44*H122 -
                                       52*H130 + 52*H131 - 80*H1m2m10 + 48*H1m200 - 32*H11m20 +
                                       132*H1112 - 24*H1120 + 24*H1121 - 92*H1200 + 12*H1210 -
                                       72*H1211 - 36*H10000 + 68*H11000 - 28*H11100 + 180*H11110 -
                                       240*H11111 - 24*H1m2*zeta2 - 74*H12*zeta2 - 68*H100*zeta2 -
                                       38*H110*zeta2 + 46*H111*zeta2 - 80*H10*zeta3 - 20*H11*zeta3 +
                                       H1*zeta4))/3. + (4*(-49 - 974*z + 304*z2)*zeta5)/3.);
      return 2 * cf;
    }
  */
}
