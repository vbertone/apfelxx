//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/massivezerocoefficientfunctions.h"
#include "apfel/tools.h"
#include "apfel/specialfunctions.h"

#include <cmath>

namespace apfel
{
  //_________________________________________________________________________________
  Cm021gNC_c::Cm021gNC_c():
    Expression()
  {
  }
  double Cm021gNC_c::Regular(double const& x) const
  {
    const double x2 = x * x;
    const double cm021gc = TR * ( ( 4 - 8 * x + 8 * x2 ) * log( ( 1 - x ) / x ) - 4 + 32 * x - 32 * x2 );
    return cm021gc;
  }

  //_________________________________________________________________________________
  Cm021gNC_l::Cm021gNC_l():
    Expression()
  {
  }
  double Cm021gNC_l::Regular(double const& x) const
  {
    const double x2 = x * x;
    const double cm021gl =  TR * ( 4 - 8 * x + 8 * x2 );
    return cm021gl;
  }

  //_________________________________________________________________________________
  Cm0L1gNC_c::Cm0L1gNC_c():
    Expression()
  {
  }
  double Cm0L1gNC_c::Regular(double const& x) const
  {
    const double cm0l1gc = TR * 16 * x * ( 1 - x );
    return cm0l1gc;
  }

  //_________________________________________________________________________________
  Cm022nsNC_c::Cm022nsNC_c():
    Expression()
  {
  }
  double Cm022nsNC_c::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double dlx  = log(x);
    const double dlx2 = dlx * dlx;
    const double dlm  = log(1-x);
    const double dlm2 = dlm * dlm;
    const double spx  = wgplg(1,1,1-x);
    const double a2 = ( 1 + 13 * x ) * dlm - ( 3 + 23 * x ) * dlx + 29. / 6 - 295 * x / 6
      + ( 1 + x2 ) * ( - 4 * spx - 8 * dlx * dlm + 6 * dlx2 + 67 * dlx / 3) / ( 1 - x );
    const double b2 = ( - 1 - x ) * ( - 4 * zeta2 + 2 * dlm2 - 29 * dlm / 3 + 359. / 18 );
    const double cm022nscr = 2 * CF * TR * ( a2 + b2 ) / 3;
    return cm022nscr;
  }
  double Cm022nsNC_c::Singular(double const& x) const
  {
    const double dlm  = log( 1 - x );
    const double dlm2 = dlm * dlm;
    const double z    = 2 / ( 1 - x );
    const double cm022nscs = z * 2 * CF * TR * ( - 4 * zeta2 + 2 * dlm2 - 29 * dlm / 3 + 359. / 18 ) / 3;
    return cm022nscs;
  }
  double Cm022nsNC_c::Local(double const& x) const
  {
    const double dlm  = log( 1 - x );
    const double dlm2 = dlm * dlm;
    const double dlm3 = dlm2 * dlm;
    const double cm022nscl = 4 * CF * TR * ( - 4 * zeta2 * dlm + 2 * dlm3 / 3 - 29 * dlm2 / 6 + 359 * dlm / 18 ) / 3
      + CF * TR * ( 268 * zeta2 / 9 + 265. / 9 );
    return cm022nscl;
  }

  //_________________________________________________________________________________
  Cm022nsNC_l::Cm022nsNC_l():
    Expression()
  {
  }
  double Cm022nsNC_l::Regular(double const& x) const
  {
    const double x2  = x * x;
    const double dlx = log(x);
    const double dlm = log(1-x);
    const double a1  = ( -8 * ( 1 + x2 ) * dlx / ( 1 - x ) + 13 * x + 1 ) + ( - 1 - x ) * ( 4 * dlm - 29. / 3 );
    const double cm022nslr = 2 * CF * TR * a1 / 3;
    return cm022nslr;
  }
  double Cm022nsNC_l::Singular(double const& x) const
  {
    const double dlm = log( 1 - x );
    const double z   = 2 / ( 1 - x );
    const double cm022nsls = z * 2 * CF * TR * ( 4 * dlm - 29. / 3 ) / 3;
    return cm022nsls;
  }
  double Cm022nsNC_l::Local(double const& x) const
  {
    const double dlm  = log( 1 - x );
    const double dlm2 = dlm * dlm;
    const double cm022nsll = 4 * CF * TR * ( 2 * dlm2 - 29 * dlm / 3 ) / 3 - CF * TR * ( 32 * zeta2 / 3 + 38. / 3 );
    return cm022nsll;
  }

  //_________________________________________________________________________________
  Cm022nsNC_l2::Cm022nsNC_l2():
    Expression()
  {
  }
  double Cm022nsNC_l2::Regular(double const& x) const
  {
    const double b1 = 2 * ( - 1 - x );
    const double cm022nsl2r = 2 * CF * TR * b1 / 3;
    return cm022nsl2r;
  }
  double Cm022nsNC_l2::Singular(double const& x) const
  {
    const double z = 2 / ( 1 - x );
    const double cm022nsl2s = z * 4 * CF * TR / 3;
    return cm022nsl2s;
  }
  double Cm022nsNC_l2::Local(double const& x) const
  {
    const double dlm = log( 1 - x );
    const double cm022nsl2l = 8 * CF * TR * dlm / 3 + 2 * CF * TR;
    return cm022nsl2l;
  }

  //_________________________________________________________________________________
  Cm0L2nsNC_c::Cm0L2nsNC_c():
    Expression()
  {
  }
  double Cm0L2nsNC_c::Regular(double const& x) const
  {
    const double cm0l2nsc = 16 * CF * TR * ( x * log( 1 - x ) - 2 * x * log(x) - 25 * x / 6 + 1 ) / 3;
    return cm0l2nsc;
  }

  //_________________________________________________________________________________
  Cm0L2nsNC_l::Cm0L2nsNC_l():
    Expression()
  {
  }
  double Cm0L2nsNC_l::Regular(double const& x) const
  {
    const double cm0l2nsl = 16 * CF * TR * x / 3;
    return cm0l2nsl;
  }

  //_________________________________________________________________________________
  Cm022psNC_c::Cm022psNC_c():
    Expression()
  {
  }
  double Cm022psNC_c::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double dlx  = log(x);
    const double dlx2 = dlx*dlx;
    const double dlx3 = dlx*dlx2;
    const double dlm  = log(1-x);
    const double dlm2 = dlm*dlm;
    const double dlp  = log(1+x);
    const double s11  = wgplg(1,1,1-x);
    const double s12  = wgplg(1,2,1-x);
    const double s21  = wgplg(2,1,1-x);
    const double s11m = wgplg(1,1,-x);
    const double a3   = ( 1 + x ) * ( 16 * dlx3 / 3 - 16 * dlx2 * dlm + 8 * dlx * dlm2
				    - 32 * zeta2 * dlx + 16 * dlm * s11 + 32 * s12 - 16 * s21 )
      + ( 40 * x - 16 * x2 ) * dlx2 + 32 * x2 * dlx * dlm + ( 280. / 3 - 704 * x2 / 9 - 88 * x ) * dlx
      + ( 4 - 16 * x2 / 3 - 4 * x + 16 / x / 3 ) * dlm2
      + ( - 208. / 3 - 64 * x2 / 9 + 160 * x / 3 + 208 / x / 9 ) * dlm
      + ( - 16 + 64 * x2 / 3 - 16 * x - 32 / x ) * zeta2
      + ( 16 - 16 * x + 64 / x / 3 + 32 * x2 / 3 ) * s11
      + ( - 32 - 32 * x2 / 3 - 32 * x - 32 / x / 3 ) * ( s11m + dlx * dlp )
      + 304. / 9 + 832 * x2 / 9 - 1216 * x / 9 + 80 / x / 9;
    const double c2ps2am0_a0 = CF * TR * a3;
    return c2ps2am0_a0;
  }

  //_________________________________________________________________________________
  Cm022psNC_l::Cm022psNC_l():
    Expression()
  {
  }
  double Cm022psNC_l::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double dlx  = log(x);
    const double dlx2 = dlx*dlx;
    const double dlm  = log(1-x);
    const double s11  = wgplg(1,1,1-x);
    const double a2   = - 64 * x2 / 9 + 160 * x / 3 - 208. / 3 + 208 / x / 9 + 32 * x2 * dlx
      + 16 * ( 1 + x ) * ( - dlx2 + dlx * dlm + s11 )
      + ( - 32 * x2 / 3 - 8 * x + 8 + 32 / x / 3 ) * dlm;
    double const c2ps2am0_aq = CF * TR * a2;
    return c2ps2am0_aq;
  }

  //_________________________________________________________________________________
  Cm022psNC_l2::Cm022psNC_l2():
    Expression()
  {
  }
  double Cm022psNC_l2::Regular(double const& x) const
  {
    const double x2  = x * x;
    const double dlx = log(x);
    const double a1  = - 16 * x2 / 3 - 4 * x + 4 + 16 / x / 3 + 8 * ( 1 + x ) * dlx;
    const double c2ps2am0_aq2 = CF * TR * a1;
    return c2ps2am0_aq2;
  }

  //_________________________________________________________________________________
  Cm022psNC_f::Cm022psNC_f():
    Expression()
  {
  }
  double Cm022psNC_f::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double dlx  = log(x);
    const double dlx2 = dlx*dlx;
    const double dlm  = log(1-x);
    const double spx  = wgplg(1,1,1-x);
    const double c2   = 128 * x2 / 3 + 16 * x / 3 - 160. / 3 + 16 / x / 3
      + 8 * ( 1 + x ) * ( - dlx2 + 2 * dlx * dlm + 2 * spx )
      + ( - 32 * x2 / 3 - 8 * x + 8 + 32 / x / 3 ) * dlm
      + ( 32 * x2 / 3 - 40 * x - 8 ) * dlx;
    const double c2ps2am0_af = CF * TR * c2;
    return c2ps2am0_af;
  }

  //_________________________________________________________________________________
  Cm022psNC_lf::Cm022psNC_lf():
    Expression()
  {
  }
  double Cm022psNC_lf::Regular(double const& x) const
  {
    const double x2  = x * x;
    const double dlx = log(x);
    const double c1  = - 32 * x2 / 3 - 8 * x + 8 + 32 / x / 3 + 16 * ( 1 + x ) * dlx;
    const double c2ps2am0_aqf = CF * TR * c1;
    return c2ps2am0_aqf;
  }

  //_________________________________________________________________________________
  Cm0L2psNC_c::Cm0L2psNC_c():
    Expression()
  {
  }
  double Cm0L2psNC_c::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double dlx  = log(x);
    const double dlx2 = dlx*dlx;
    const double dlm  = log(1-x);
    const double spx  = wgplg(1,1,1-x);
    const double a2   = 32 * x * ( dlx2 - dlx * dlm - spx )
      + ( - 32 + 64 * x2 / 3 + 32 / x / 3 ) * dlm
      + ( 32 - 64 * x2 - 32 * x ) * dlx + 32. / 3
      + 320 * x2 / 9 - 128 * x / 3 - 32 / x / 9;
    const double clps2am0_a0 = CF * TR * a2;
    return clps2am0_a0;
  }

  //_________________________________________________________________________________
  Cm0L2psNC_l::Cm0L2psNC_l():
    Expression()
  {
  }
  double Cm0L2psNC_l::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double dlx = log(x);
    const double a1  = - 32 * x * dlx - 32 + 64 * x2 / 3 + 32 / x / 3;
    const double clps2am0_aq = CF * TR * a1;
    return clps2am0_aq;
  }

  //_________________________________________________________________________________
  Cm0L2psNC_f::Cm0L2psNC_f():
    Expression()
  {
  }
  double Cm0L2psNC_f::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double dlx = log(x);
    const double a1  = - 32 * x * dlx - 32 + 64 * x2 / 3 + 32 / x / 3;
    const double clps2am0_af = CF * TR * a1;
    return clps2am0_af;
  }

  //_________________________________________________________________________________
  Cm022gNC_c::Cm022gNC_c():
    Expression()
  {
  }
  double Cm022gNC_c::Regular(double const& x) const
  {
    const double x2     = x * x;
    const double x3     = x * x2;
    const double dlx    = log(x);
    const double dlx2   = dlx * dlx;
    const double dlx3   = dlx2 * dlx;
    const double dlm    = log(1-x);
    const double dlm2   = dlm * dlm;
    const double dlm3   = dlm2 * dlm;
    const double dlp    = log(1+x);
    const double dlp2   = dlp * dlp;
    const double s121mx = wgplg(1,2,1-x);
    const double s12mx  = wgplg(1,2,-x);
    const double s211mx = wgplg(2,1,1-x);
    const double s21mx  = wgplg(2,1,-x);
    const double s111mx = wgplg(1,1,1-x);
    const double s11mx  = wgplg(1,1,-x);
    const double z      = ( 1 - x ) / ( 1 + x );
    const double s21z   = wgplg(2,1,z);
    const double s21mz  = wgplg(2,1,-z);
    const double a31 = dlx3 * ( 16. / 3 + 16 * x ) +
      dlx2 * dlm * ( - 8 + 16 * x2 - 64 * x ) + 
      dlx2 * dlp * ( 12 + 32 * x2 + 24 * x ) + 
      dlx2 * ( - 114 * x2 + 184 * x ) + 
      dlx * dlm2 * ( - 16 * x2 + 48 * x ) + 
      dlx * dlm * dlp * ( - 16 - 32 * x2 - 32 * x ) + 
      dlx * dlm * ( 16 + 292 * x2 - 288 * x ) + 
      dlx * dlp2 * ( 8 + 16 * x );
    const double a32 = dlx * dlp * ( - 48 + 208 * x2 / 3 + 16 * x - 32 / x / 3 ) + 
      dlx * zeta2 * ( - 16 + 32 * x2 - 160 * x ) + 
      dlx * s111mx * (  + 32 * x ) + 
      dlx * s11mx * ( 24 + 32 * x2 + 48 * x ) + 
      dlx * ( 292. / 3 - 5780 * x2 / 9 + 332 * x ) + 
      dlm2 * ( - 6 - 214 * x2 / 3 + 64 * x + 16 / x / 3 ) + 
      zeta2 * dlm * ( - 40 - 64 * x2 + 48 * x );
    const double a33 = dlm * s111mx * ( 16 + 64 * x ) + 
      dlm * s11mx * ( - 16 - 32 * x2 - 32 * x ) + 
      dlm * ( - 112. / 3 + 2996 * x2 / 9 - 860 * x / 3 + 208 / x / 9 ) + 
      zeta2 * dlp * ( 8 + 16 * x ) + 
      dlp * s11mx * ( 16 + 32 * x ) + 
      s21mz * ( - 16 - 32 * x2 - 32 * x ) + 
      s21z * ( 16 + 32 * x2 + 32 * x );
    const double a34 = zeta2 * ( - 4 + 796 * x2 / 3 - 208 * x - 32 / x ) + 
      zeta3 * ( - 12 - 8 * x2 - 56 * x ) + 
      s111mx * ( 20 + 80 * x2 / 3 - 64 * x + 64 / x / 3 ) + 
      s211mx * ( - 16 - 128 * x ) + 
      s121mx * ( 40 + 144 * x ) + 
      s11mx * ( - 48 + 208 * x2 / 3 + 16 * x - 32 / x / 3 ) + 
      s21mx * ( - 24 - 48 * x ) + 
      s12mx * ( 16 + 32 * x ) + 80 / x / 9 + 
      466. / 9 - 878 * x2 / 9 + 260 * x / 9;
    const double a3 = a31 + a32 + a33 + a34;
    const double b31 = dlx3 * ( - 8. / 3 - 32 * x2 / 3 + 16 * x / 3 ) + 
      dlx2 * dlm * ( 16 + 48 * x2 - 32 * x ) + 
      dlx2 * dlp * ( 16 + 16 * x2 + 32 * x ) + 
      dlx2 * ( - 4 - 96 * x3 / 5 - 52 * x2 + 8 * x / 3 ) + 
      dlx * dlm2 * ( - 20 - 48 * x2 + 40 * x ) + 
      dlx * dlm * ( 24 + 168 * x2 - 160 * x ) + 
      dlx * dlp2 * ( - 32 - 32 * x2 - 64 * x ) + 
      dlx * dlp * ( 96 + 192 * x3 / 5 + 128 * x / 3 );
    const double b32 = 16 * dlx * dlp / x2 / 15 - 16 * dlx / x / 15 + 
      dlx * zeta2 * ( 32 + 64 * x2 - 64 * x ) + 
      dlx * s111mx * ( 32 * x2 ) + 
      dlx * s11mx * ( - 32 - 32 * x2 + 64 * x ) + 
      dlx * ( - 712. / 15 - 672 * x2 / 5 + 136 * x / 5 ) + 
      dlm3 * ( 8 + 16 * x2 - 16 * x ) + 
      dlm2 * ( - 22 - 84 * x2 + 88 * x ) + 
      zeta2 * dlm * ( - 32 * x2 );
    const double b33 = dlm * s111mx * ( 8 - 16 * x ) + 
      dlm * ( 28 + 96 * x2 - 132 * x ) + 
      zeta2 * dlp * ( - 32 - 32 * x2 - 64 * x ) + 
      dlp * s11mx * ( - 64 - 64 * x2 - 128 * x ) + 
      zeta2 * ( 48 + 192 * x3 / 5 + 104 * x2 - 208 * x / 3 ) + 
      zeta3 * ( 112 + 192 * x2 - 96 * x );
    const double b34 = s111mx * ( - 24 + 64 * x2 - 48 * x ) + 
      s211mx * ( - 24 - 32 * x2 + 48 * x ) + 
      s121mx * ( - 32 + 64 * x ) + 
      s11mx * ( 96 + 192 * x3 / 5 + 128 * x / 3 ) + 
      16 * s11mx / x2 / 15 + 16 / x / 15 + 
      s21mx * ( 96 + 96 * x2 - 64 * x ) + 
      s12mx * ( - 64 - 64 * x2 - 128 * x ) - 
      904. / 15 + 328 * x2 / 5 + 68 * x / 5;
    const double b3 = b31 + b32 + b33 + b34;
    const double c2g2am0_a0 = TR * ( CA * a3 + CF * b3 );
    return c2g2am0_a0;
  }

  //_________________________________________________________________________________
  Cm022gNC_l::Cm022gNC_l():
    Expression()
  {
  }
  double Cm022gNC_l::Regular(double const& x) const
  {
    const double x2     = x * x;
    const double dlx    = log(x);
    const double dlx2   = dlx * dlx;
    const double dlm    = log(1-x);
    const double dlm2   = dlm * dlm;
    const double dlp    = log(1+x);
    const double s111mx = wgplg(1,1,1-x);
    const double s11mx  = wgplg(1,1,-x);
    const double a2     =  - ( 16 + 32 * x2 ) * zeta2 + 1628 * x2 / 9
      - 368 * x / 3 - 220. / 3 + 208 / x / 9
      + ( 16 * x2 - 16 * x + 8 ) * dlm2 - ( 48 * x + 16 ) * dlx2
      + ( - 536 * x2 / 3 + 160 * x - 8 + 32 / x / 3 ) * dlm
      + ( 200 * x2 - 192 * x ) * dlx + ( 96 * x - 32 * x2 ) * dlx * dlm
      + ( 64 * x + 16 ) * s111mx - ( 32 * x2 + 32 * x + 16 ) * ( s11mx + dlx * dlp );
    const double b2     = ( - 64 * x2 + 64 * x - 32 ) * zeta2 + 16 * x2 - 68 * x + 36
      + ( 32 * x2 - 32 * x + 16 ) * dlm2 + ( 32 * x2 - 16 * x + 8 ) * dlx2
      + ( - 80 * x2 + 96 * x - 28 ) * dlm + ( 80 * x2 - 48 * x + 8 ) * dlx
      + ( - 64 * x2 + 48 * x - 24 ) * dlx * dlm + ( 8 - 16 * x ) * s111mx;
    const double c2g2am0_aq = TR * (  CA * a2 + CF * b2 );
    return c2g2am0_aq;
  }

  //_________________________________________________________________________________
  Cm022gNC_l2::Cm022gNC_l2():
    Expression()
  {
  }
  double Cm022gNC_l2::Regular(double const& x) const
  {
    const double x2  = x * x;
    const double dlx = log(x);
    const double dlm = log(1-x);
    const double a1  = 16 / x / 3 - 124 * x2 / 3 + 32 * x + 4
      + ( 16 * x2 - 16 * x + 8 ) * dlm + ( 32 * x + 8 ) * dlx;
    const double b1  = - 2 + 8 * x + ( 16 * x2 - 16 * x + 8 ) * dlm
      + ( - 16 * x2 + 8 * x - 4 ) * dlx;
    const double c2g2am0_aq2 = TR * ( CA * a1 + CF * b1 );
    return c2g2am0_aq2;
  }

  //_________________________________________________________________________________
  Cm022gNC_f::Cm022gNC_f():
    Expression()
  {
  }
  double Cm022gNC_f::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double dlx  = log(x);
    const double dlx2 = dlx * dlx;
    const double dlm  = log(1-x);
    const double dlm2 = dlm * dlm;
    const double s11  = wgplg(1,1,1-x);
    const double c2   = ( - 16 + 32 * x - 32 * x2 ) * zeta2
      + 1124 * x2 / 3 - 968 * x / 3 - 172. / 3 + 16 / x / 3
      + ( 16 - 32 * x + 32 * x2 ) * dlm2 - ( 8 + 32 * x ) * dlx2
      + ( 248 * x2 / 3 - 256 * x - 8 ) * dlx
      + ( 32 / x / 3 - 8 + 192 * x - 632 * x2 / 3 ) * dlm
      + ( 96 * x - 32 * x2 ) * dlx * dlm + ( 64 * x + 16 ) * s11;
    const double c2g2am0_af = TR * CA * c2;
    return c2g2am0_af;
  }

  //_________________________________________________________________________________
  Cm022gNC_lf::Cm022gNC_lf():
    Expression()
  {
  }
  double Cm022gNC_lf::Regular(double const& x) const
  {
    const double x2  = x * x;
    const double dlx = log(x);
    const double dlm = log(1-x);
    const double c1  = - 248 * x2 / 3 + 64 * x + 8 + 32 / x / 3
      + ( 32 * x2 - 32 * x + 16 ) * dlm + ( 64 * x + 16 ) * dlx;
    const double c2g2am0_aqf = TR * CA * c1;
    return c2g2am0_aqf;
  }

  //_________________________________________________________________________________
  Cm0L2gNC_c::Cm0L2gNC_c():
    Expression()
  {
  }
  double Cm0L2gNC_c::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double x3   = x * x2;
    const double dlx  = log(x);
    const double dlx2 = dlx * dlx;
    const double dlm  = log(1-x);
    const double dlm2 = dlm * dlm;
    const double dlp  = log(1+x);
    const double s11  = wgplg(1,1,1-x);
    const double s11m = wgplg(1,1,-x);
    const double a2   = 96 * x * dlx2 + ( 64 * x2 - 192 * x ) * dlx * dlm
      + ( 32 - 416 * x2 + 256 * x ) * dlx
      + 32 * x * ( 1 - x ) * dlm2 + ( - 32 + 928 * x2 / 3 - 288 * x + 32 / x / 3 ) * dlm
      + 64 * x2 * zeta2 - 128 * x * s11
      + 64 * x * ( 1 + x ) * ( s11m + dlx * dlp )
      + 32. / 3 - 1696 * x2 / 9 + 544 * x / 3 - 32 / x / 9;
    const double b2   = - ( 64 * x3 / 5 + 64 * x / 3 ) * dlx2
      + 32 * x * ( s11 + dlx * dlm )
      + ( - 208. / 15 + 192 * x2 / 5 - 416 * x / 5 - 64 / x / 15 ) * dlx
      + ( 16 - 64 * x2 + 48 * x ) * dlm + ( 128 * x3 / 5 - 64 * x / 3 ) * zeta2
      + ( 128 * x3 / 5 - 64 * x / 3 + 64 / x2 / 15 ) * ( s11m + dlx * dlp )
      - 256. / 15 + 672 * x2 / 5 - 608 * x / 5 + 64 / x / 15;
      const double clg2am0_a0 = TR * ( CA * a2 + CF * b2 );
      return clg2am0_a0;
  }

  //_________________________________________________________________________________
  Cm0L2gNC_l::Cm0L2gNC_l():
    Expression()
  {
  }
  double Cm0L2gNC_l::Regular(double const& x) const
  {
    const double x2  = x * x;
    const double dlx = log(x);
    const double dlm = log(1-x);
    const double a1  = 64 * x * ( 1 - x ) * dlm - 128 * x * dlx - 32 - 160 * x + 544 * x2 / 3 + 32 / x / 3;
    const double b1  = 32 * x * dlx + 16 * ( 1 - 2 * x2 + x );
    const double clg2am0_aq = TR * ( CA * a1 + CF * b1 );
    return clg2am0_aq;
  }

  //_________________________________________________________________________________
  Cm0L2gNC_f::Cm0L2gNC_f():
    Expression()
  {
  }
  double Cm0L2gNC_f::Regular(double const& x) const
  {
    const double x2  = x * x;
    const double dlx = log(x);
    const double dlm = log(1-x);
    const double c1  = 64 * x * ( 1 - x ) * dlm - 128 * x * dlx - 32 - 160 * x + 544 * x2 / 3 + 32 / x / 3;
    const double clg2am0_af = TR * CA * c1;
    return clg2am0_af;
  }
}
