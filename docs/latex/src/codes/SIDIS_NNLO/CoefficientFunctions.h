#include <apfel/apfelxx.h>

double Dn(int const& n, double const& y, bool const& integral)
{
  if (integral)
    return - pow(log(y), n + 1) / ( n + 1 );
  else
    return pow(log(y), n) / y;
}

double InvTanInt(double const& y)
{
  apfel::Integrator I{[] (double const& t) -> double { return atan(t) / t; }};
  return I.integrate(0, y, apfel::eps5);
}

/*
const double sqrtxz1  = sqrt(1 - 2*z + z*z + 4*x*z);
const double sqrtxz1i = 1. / sqrtxz1;
const double poly2    = 1 + 2*x + x*x - 4*x*z;
const double poly2i   = 1. / poly2;
const double poly2i2  = poly2i * poly2i;
const double sqrtxz2  = sqrt(poly2);
const double sqrtxz2i = 1. / sqrtxz2;
const double sqrtxz3  = sqrt(x/z);

const double NC   = apfel::NC;
const double NC2  = NC * NC;
const double NCi  = 1. / NC;
const double NCi2 = NCi * NCi;

const double l2  = log(2);
const double l22 = l2 * l2;
const double l23 = l2 * l22;

const double pi  = M_PI;
const double pi2 = pi * pi;
const double pi4 = pi2 * pi2;

const double x2 = x * x;
const double x3 = x * x2;
const double x4 = x * x3;
const double x5 = x * x4;
const double x6 = x * x5;
const double x7 = x * x6;

const double xi  = 1. / x;
const double xi2 = xi * xi;

const double z2 = z * z;
const double z3 = z * z2;

const double zi  = 1. / z;
const double zi2 = zi * zi;

const double omxi = 1. / ( 1 - x );
const double opxi = 1. / ( 1 + x );

const double omzi  = 1. / ( 1 - z );
const double opzi  = 1. / ( 1 + z );
const double mopzi = 1. / ( - 1 + z );

const double lx  = log(x);
const double lx2 = lx * lx;
const double lx3 = lx * lx2;

const double lz  = log(z);
const double lz2 = lz * lz;
const double lz3 = lz * lz2;

const double lomx  = log(1 - x);
const double lomx2 = lomx * lomx;
const double lomx3 = lomx * lomx2;

const double lomz  = log(1 - z);
const double lomz2 = lomz * lomz;
const double lomz3 = lomz * lomz2;

const double lopx  = log(1 + x);
const double lopx2 = lopx * lopx;
const double lopx3 = lopx * lopx2;

const double lopz  = log(1 + z);
const double lopz2 = lopz * lopz;

const double Li2x  = apfel::dilog(x);
const double Li2mx = apfel::dilog(-x);

const double Li2z  = apfel::dilog(z);
const double Li2mz = apfel::dilog(-z);

const double Li3x       = apfel::wgplg(2,1,x);
const double Li3omx     = apfel::wgplg(2,1,1 - x);
const double Li3omxh    = apfel::wgplg(2,1,0.5 - x/2.);
const double Li3opxh    = apfel::wgplg(2,1,(1 + x)/2.);
const double Li3omxopxi = apfel::wgplg(2,1,-((-1 + x)*opxi));

const double Li3z       = apfel::wgplg(2,1,z);
const double Li3mz      = apfel::wgplg(2,1,-z);
const double Li3omz     = apfel::wgplg(2,1,1 - z);
const double Li3omzh    = apfel::wgplg(2,1,0.5 - z/2.);
const double Li3opzh    = apfel::wgplg(2,1,(1 + z)/2.);
const double Li3opzi    = apfel::wgplg(2,1,opzi);
const double Li3omzopzi = apfel::wgplg(2,1,-((-1 + z)*opzi));

const double omxmzi  = 1. / ( 1 - x - z );
const double omxmzi2 = omxmzi * omxmzi;

const double xmzi  = 1. / ( x - z );
const double xmzi2 = xmzi * xmzi;

const double spec1i = 1. / ( 1 + sqrtxz1 - z );

const double lomxmz  = log(1 - x - z);
const double lmopxpz = log(-1 + x + z);
const double lxmz    = log(x - z);
const double lmxpz   = log(-x + z);
const double lxpz    = log(x + z);
const double lopxz   = log(1 + x*z);
const double lopxzi  = log(1 + x*zi);

const double li2omxzi = apfel::dilog(1 - x*zi);

const double lspec1   = log(1 + sqrtxz1 - z);
const double lspec1_2 = lspec1 * lspec1;
const double lspec2   = log(1 + sqrtxz1 + z);
const double lspec3   = log(sqrtxz3);
const double lspec4   = log(sqrtxz3*z);
const double lspec5   = log(1 - sqrtxz2 + x);
const double lspec6   = log(1 + sqrtxz2 + x);
const double lspec7   = log(1 - 2*z + 4*x*z + z2);
const double lspec7_2 = lspec7 * lspec7;

const double li2spec1  = apfel::dilog(0.5 - sqrtxz1/2. - z/2.);
const double li2spec2  = apfel::dilog(0.5 - sqrtxz1/2. + z/2.);
const double li2spec3  = apfel::dilog(0.5 - zi/2. - (sqrtxz1*zi)/2.);
const double li2spec4  = apfel::dilog(0.5 + zi/2. - (sqrtxz1*zi)/2.);
const double li2spec5  = apfel::dilog(0.5 - sqrtxz2/2. - x/2.);
const double li2spec6  = apfel::dilog(0.5 + sqrtxz2/2. - x/2.);
const double li2spec7  = apfel::dilog(0.5 - xi/2. - (sqrtxz2*xi)/2.);
const double li2spec8  = apfel::dilog(0.5 - xi/2. + (sqrtxz2*xi)/2.);
const double li2spec9  = apfel::dilog((1 - z)*omxi);
const double li2spec10 = apfel::dilog(x*(1 - z)*omxi*zi);
const double li2spec11 = apfel::dilog((1 - x)*omzi);
const double li2spec12 = apfel::dilog((1 - x)*z*xi*omzi);
const double li2spec13 = apfel::dilog(z*omxi);
const double li2spec14 = apfel::dilog((1 - z)*z*omxi*xi);
const double li2spec15 = apfel::dilog(x*z*omxi*omzi);
const double li2spec16 = apfel::dilog((1 - x)*zi);
const double li2spec17 = apfel::dilog((1 - x)*(1 - z)*xi*zi);
const double li2spec18 = apfel::dilog((1 - x)*x*omzi*zi);
const double li2spec19 = apfel::dilog(-(x*z));
const double li2spec20 = apfel::dilog(-(x*zi));
const double li2spec21 = apfel::dilog((1 - z)*sqrtxz1i);
const double li2spec22 = apfel::dilog(-spec1i + sqrtxz1*spec1i + z*spec1i);
const double li2spec23 = apfel::dilog(sqrtxz1*mopzi);

const double atanspec1 = atan(sqrtxz3);
const double atanspec2 = atan(sqrtxz3*z);

const double itani1 = InvTanInt(-sqrtxz3);
const double itani2 = InvTanInt(sqrtxz3);
const double itani3 = InvTanInt(sqrtxz3*z);

const double Tr1 = (z < 1 - x ? 1 : 0);
const double Tr2 = (z > 1 - x ? 1 : 0);

const double Tt1 = (z > x ? 1 : 0);
const double Tt2 = (z < x ? 1 : 0);

const double Tu1 = (z < 1 - x && z < x ? 1 : 0);
const double Tu2 = (z > 1 - x && z < x ? 1 : 0);
const double Tu3 = (z < 1 - x && z > x ? 1 : 0);
const double Tu4 = (z > 1 - x && z > x ? 1 : 0);
*/

//------------------------------------------------------
double C1LQ2Q_Lx_Lz(int const&)
{
  return 0;
};
double C1LQ2G_Lx_Lz(int const&)
{
  return 0;
};
double C1LG2Q_Lx_Lz(int const&)
{
  return 0;
};
double C1TQ2Q_Lx_Lz(int const&)
{
  const double NC  = apfel::NC;
  const double NCi = 1. / NC;
  return -4*NC + 4*NCi;
};
double C1TQ2G_Lx_Lz(int const&)
{
  return 0;
};
double C1TG2Q_Lx_Lz(int const&)
{
  return 0;
};
double C2LQ2G_Lx_Lz(int const&)
{
  return 0;
};
double C2LG2G_Lx_Lz(int const&)
{
  return 0;
};
double C2LG2Q_Lx_Lz(int const&)
{
  return 0;
};
double C2LQ2QNS_Lx_Lz(int const&)
{
  return 0;
};
double C2LQ2QPS_Lx_Lz(int const&)
{
  return 0;
};
double C2LQ2QB_Lx_Lz(int const&)
{
  return 0;
};
double C2LQ2QP1_Lx_Lz(int const&)
{
  return 0;
};
double C2LQ2QP2_Lx_Lz(int const&)
{
  return 0;
};
double C2LQ2QP3_Lx_Lz(int const&)
{
  return 0;
};
double C2TQ2G_Lx_Lz(int const&)
{
  return 0;
};
double C2TG2G_Lx_Lz(int const&)
{
  return 0;
};
double C2TG2Q_Lx_Lz(int const&)
{
  return 0;
};
double C2TQ2QNS_Lx_Lz(int const& nf)
{
  const double NC   = apfel::NC;
  const double NC2  = NC * NC;
  const double NCi  = 1. / NC;
  const double NCi2 = NCi * NCi;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double pi4 = pi2 * pi2;
  return (180 + 22860*NC*nf + 2880*apfel::zeta3 + 2880*NC*nf*apfel::zeta3 -  
     69165*NC2 + 29520*apfel::zeta3*NC2 + 320*pi2 +  
     1520*NC*nf*pi2 - 5540*NC2*pi2 -  
     20*nf*NCi*(1143 + 144*apfel::zeta3 + 76*pi2) +  
     NCi2*(68985 - 32400*apfel::zeta3 + 5220*pi2 - 168*pi4) -  
     72*pi4 + 240*NC2*pi4)/8640.;
};
double C2TQ2QPS_Lx_Lz(int const&)
{
  return 0;
};
double C2TQ2QB_Lx_Lz(int const&)
{
  return 0;
};
double C2TQ2QP1_Lx_Lz(int const&)
{
  return 0;
};
double C2TQ2QP2_Lx_Lz(int const&)
{
  return 0;
};
double C2TQ2QP3_Lx_Lz(int const&)
{
  return 0;
};

//------------------------------------------------------
double C1LQ2Q_Sx_Sz(int const&, double const&, double const&, bool const&, bool const&)
{
  return 0;
};
double C1LQ2G_Sx_Sz(int const&, double const&, double const&, bool const&, bool const&)
{
  return 0;
};
double C1LG2Q_Sx_Sz(int const&, double const&, double const&, bool const&, bool const&)
{
  return 0;
};

double C1TQ2Q_Sx_Sz(int const&, double const& x, double const& z, bool const& intx, bool const& intz)
{
  const double NC  = apfel::NC;
  const double NCi = 1. / NC;
  return Dn(0,1 - x, intx)*Dn(0,1 - z, intz)*(NC - NCi);
};
double C1TQ2G_Sx_Sz(int const&, double const&, double const&, bool const&, bool const&)
{
  return 0;
};
double C1TG2Q_Sx_Sz(int const&, double const&, double const&, bool const&, bool const&)
{
  return 0;
};
double C2LQ2G_Sx_Sz(int const&, double const&, double const&, bool const&, bool const&)
{
  return 0;
};
double C2LG2G_Sx_Sz(int const&, double const&, double const&, bool const&, bool const&)
{
  return 0;
};
double C2LG2Q_Sx_Sz(int const&, double const&, double const&, bool const&, bool const&)
{
  return 0;
};
double C2LQ2QNS_Sx_Sz(int const&, double const&, double const&, bool const&, bool const&)
{
  return 0;
};
double C2LQ2QPS_Sx_Sz(int const&, double const&, double const&, bool const&, bool const&)
{
  return 0;
};
double C2LQ2QB_Sx_Sz(int const&, double const&, double const&, bool const&, bool const&)
{
  return 0;
};
double C2LQ2QP1_Sx_Sz(int const&, double const&, double const&, bool const&, bool const&)
{
  return 0;
};
double C2LQ2QP2_Sx_Sz(int const&, double const&, double const&, bool const&, bool const&)
{
  return 0;
};
double C2LQ2QP3_Sx_Sz(int const&, double const&, double const&, bool const&, bool const&)
{
  return 0;
};
double C2TQ2G_Sx_Sz(int const&, double const&, double const&, bool const&, bool const&)
{
  return 0;
};
double C2TG2G_Sx_Sz(int const&, double const&, double const&, bool const&, bool const&)
{
  return 0;
};
double C2TG2Q_Sx_Sz(int const&, double const&, double const&, bool const&, bool const&)
{
  return 0;
};
double C2TQ2QNS_Sx_Sz(int const& nf, double const& x, double const& z, bool const& intx, bool const& intz)
{
  const double NC   = apfel::NC;
  const double NC2  = NC * NC;
  const double NCi  = 1. / NC;
  const double NCi2 = NCi * NCi;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  return (3*Dn(0,1 - z, intz)*Dn(1,1 - x, intx)* 
      (11 + 2*NC*nf - 2*nf*NCi - 11*NC2) +  
     3*Dn(0,1 - x, intx)*Dn(1,1 - z, intz)* 
      (11 + 2*NC*nf - 2*nf*NCi - 11*NC2) +  
     54*Dn(1,1 - x, intx)*Dn(1,1 - z, intz)*(-2 + NCi2 + NC2) +  
     27*Dn(0,1 - z, intz)*Dn(2,1 - x, intx)*(-2 + NCi2 + NC2) +  
     27*Dn(0,1 - x, intx)*Dn(2,1 - z, intz)*(-2 + NCi2 + NC2) -  
     Dn(0,1 - x, intx)*Dn(0,1 - z, intz)*(-77 + 10*NC*nf - 10*nf*NCi +  
        5*NC2 - 15*pi2 + 9*NC2*pi2 +  
        6*NCi2*(12 + pi2)))/18.;
};
double C2TQ2QPS_Sx_Sz(int const&, double const&, double const&, bool const&, bool const&)
{
  return 0;
};
double C2TQ2QB_Sx_Sz(int const&, double const&, double const&, bool const&, bool const&)
{
  return 0;
};
double C2TQ2QP1_Sx_Sz(int const&, double const&, double const&, bool const&, bool const&)
{
  return 0;
};
double C2TQ2QP2_Sx_Sz(int const&, double const&, double const&, bool const&, bool const&)
{
  return 0;
};
double C2TQ2QP3_Sx_Sz(int const&, double const&, double const&, bool const&, bool const&)
{
  return 0;
};
//------------------------------------------------------
double C1LQ2Q_Lx_Sz(int const&, double const&, bool const&)
{
  return 0;
};
double C1LQ2G_Lx_Sz(int const&, double const&, bool const&)
{
  return 0;
};
double C1LG2Q_Lx_Sz(int const&, double const&, bool const&)
{
  return 0;
};
double C1TQ2Q_Lx_Sz(int const&, double const& z, bool const& intz)
{
  const double NC  = apfel::NC;
  const double NCi = 1. / NC;
  return Dn(1,1 - z, intz)*(NC - NCi);
};
double C1TQ2G_Lx_Sz(int const&, double const&, bool const&)
{
  return 0;
};
double C1TG2Q_Lx_Sz(int const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2G_Lx_Sz(int const&, double const&, bool const&)
{
  return 0;
};
double C2LG2G_Lx_Sz(int const&, double const&, bool const&)
{
  return 0;
};
double C2LG2Q_Lx_Sz(int const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QNS_Lx_Sz(int const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QPS_Lx_Sz(int const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QB_Lx_Sz(int const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QP1_Lx_Sz(int const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QP2_Lx_Sz(int const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QP3_Lx_Sz(int const&, double const&, bool const&)
{
  return 0;
};
double C2TQ2G_Lx_Sz(int const&, double const&, bool const&)
{
  return 0;
};
double C2TG2G_Lx_Sz(int const&, double const&, bool const&)
{
  return 0;
};
double C2TG2Q_Lx_Sz(int const&, double const&, bool const&)
{
  return 0;
};
double C2TQ2QNS_Lx_Sz(int const& nf, double const& z, bool const& intz)
{
  const double NC   = apfel::NC;
  const double NC2  = NC * NC;
  const double NCi  = 1. / NC;
  const double NCi2 = NCi * NCi;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  return (9*(Dn(2,1 - z, intz)*(11 + 2*NC*nf - 2*nf*NCi - 11*NC2) +  
        6*Dn(3,1 - z, intz)*(-2 + NCi2 + NC2)) -  
     6*Dn(1,1 - z, intz)*(-77 + 10*NC*nf - 10*nf*NCi + 5*NC2 -  
        15*pi2 + 9*NC2*pi2 + 6*NCi2*(12 + pi2)) 
      + Dn(0,1 - z, intz)*(404 + 56*NC*nf - 810*apfel::zeta3 + 216*apfel::zeta3*NCi2 -  
        404*NC2 + 594*apfel::zeta3*NC2 - 33*pi2 -  
        6*NC*nf*pi2 + 33*NC2*pi2 +  
        2*nf*NCi*(-28 + 3*pi2)))/108.;
};
double C2TQ2QPS_Lx_Sz(int const&, double const&, bool const&)
{
  return 0;
};
double C2TQ2QB_Lx_Sz(int const&, double const&, bool const&)
{
  return 0;
};
double C2TQ2QP1_Lx_Sz(int const&, double const&, bool const&)
{
  return 0;
};
double C2TQ2QP2_Lx_Sz(int const&, double const&, bool const&)
{
  return 0;
};
double C2TQ2QP3_Lx_Sz(int const&, double const&, bool const&)
{
  return 0;
};
//------------------------------------------------------
double C1LQ2Q_Sx_Lz(int const&, double const&, bool const&)
{
  return 0;
};
double C1LQ2G_Sx_Lz(int const&, double const&, bool const&)
{
  return 0;
};
double C1LG2Q_Sx_Lz(int const&, double const&, bool const&)
{
  return 0;
};
double C1TQ2Q_Sx_Lz(int const&, double const& x, bool const& intx)
{
  const double NC  = apfel::NC;
  const double NCi = 1. / NC;
  return Dn(1,1 - x, intx)*(NC - NCi);
};
double C1TQ2G_Sx_Lz(int const&, double const&, bool const&)
{
  return 0;
};
double C1TG2Q_Sx_Lz(int const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2G_Sx_Lz(int const&, double const&, bool const&)
{
  return 0;
};
double C2LG2G_Sx_Lz(int const&, double const&, bool const&)
{
  return 0;
};
double C2LG2Q_Sx_Lz(int const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QNS_Sx_Lz(int const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QPS_Sx_Lz(int const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QB_Sx_Lz(int const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QP1_Sx_Lz(int const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QP2_Sx_Lz(int const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QP3_Sx_Lz(int const&, double const&, bool const&)
{
  return 0;
};
double C2TQ2G_Sx_Lz(int const&, double const&, bool const&)
{
  return 0;
};
double C2TG2G_Sx_Lz(int const&, double const&, bool const&)
{
  return 0;
};
double C2TG2Q_Sx_Lz(int const&, double const&, bool const&)
{
  return 0;
};
double C2TQ2QNS_Sx_Lz(int const& nf, double const& x, bool const& intx)
{
  const double NC   = apfel::NC;
  const double NC2  = NC * NC;
  const double NCi  = 1. / NC;
  const double NCi2 = NCi * NCi;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  return (9*(Dn(2,1 - x, intx)*(11 + 2*NC*nf - 2*nf*NCi - 11*NC2) +  
        6*Dn(3,1 - x, intx)*(-2 + NCi2 + NC2)) -  
     6*Dn(1,1 - x, intx)*(-77 + 10*NC*nf - 10*nf*NCi + 5*NC2 -  
        15*pi2 + 9*NC2*pi2 + 6*NCi2*(12 + pi2)) 
      + Dn(0,1 - x, intx)*(404 + 56*NC*nf - 810*apfel::zeta3 + 216*apfel::zeta3*NCi2 -  
        404*NC2 + 594*apfel::zeta3*NC2 - 33*pi2 -  
        6*NC*nf*pi2 + 33*NC2*pi2 +  
        2*nf*NCi*(-28 + 3*pi2)))/108.;
};
double C2TQ2QPS_Sx_Lz(int const&, double const&, bool const&)
{
  return 0;
};
double C2TQ2QB_Sx_Lz(int const&, double const&, bool const&)
{
  return 0;
};
double C2TQ2QP1_Sx_Lz(int const&, double const&, bool const&)
{
  return 0;
};
double C2TQ2QP2_Sx_Lz(int const&, double const&, bool const&)
{
  return 0;
};
double C2TQ2QP3_Sx_Lz(int const&, double const&, bool const&)
{
  return 0;
};
//------------------------------------------------------
double C1LQ2Q_Lx_Rz(int const&, double const&)
{
  return 0;
};
double C1LQ2G_Lx_Rz(int const&, double const&)
{
  return 0;
};
double C1LG2Q_Lx_Rz(int const&, double const&)
{
  return 0;
};
double C1TQ2Q_Lx_Rz(int const&, double const& z)
{
  const double NC  = apfel::NC;
  const double NCi = 1. / NC;
  const double lomz = log(1 - z);
  const double lz   = log(z);
  const double omzi = 1. / ( 1 - z );
  return -0.5*((NC - NCi)* 
     (-1 + z + (1 + z)*lomz + lz*(1 + z - 2*omzi)));
};
double C1TQ2G_Lx_Rz(int const&, double const& z)
{
  const double NC  = apfel::NC;
  const double NCi = 1. / NC;
  const double lomz = log(1 - z);
  const double zi  = 1. / z;
  const double lz   = log(z);
  return ((NC - NCi)*(z + lomz*(-2 + z + 2*zi) +  
       lz*(-2 + z + 2*zi)))/2.;
};
double C1TG2Q_Lx_Rz(int const&, double const&)
{
  return 0;
};
double C2LQ2G_Lx_Rz(int const&, double const&)
{
  return 0;
};
double C2LG2G_Lx_Rz(int const&, double const&)
{
  return 0;
};
double C2LG2Q_Lx_Rz(int const&, double const&)
{
  return 0;
};
double C2LQ2QNS_Lx_Rz(int const&, double const&)
{
  return 0;
};
double C2LQ2QPS_Lx_Rz(int const&, double const&)
{
  return 0;
};
double C2LQ2QB_Lx_Rz(int const&, double const&)
{
  return 0;
};
double C2LQ2QP1_Lx_Rz(int const&, double const&)
{
  return 0;
};
double C2LQ2QP2_Lx_Rz(int const&, double const&)
{
  return 0;
};
double C2LQ2QP3_Lx_Rz(int const&, double const&)
{
  return 0;
};
double C2TQ2G_Lx_Rz(int const&, double const& z)
{
  const double NC   = apfel::NC;
  const double NC2  = NC * NC;
  const double NCi  = 1. / NC;
  const double NCi2 = NCi * NCi;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double l2  = log(2);
  const double l22 = l2 * l2;
  const double l23 = l2 * l22;
  const double z2  = z * z;
  const double zi  = 1. / z;
  const double lz  = log(z);
  const double lz2 = lz * lz;
  const double lz3 = lz * lz2;
  const double lomz = log(1 - z);
  const double opzi = 1. / ( 1 + z );
  const double lomz2 = lomz * lomz;
  const double lomz3 = lomz * lomz2;
  const double lopz  = log(1 + z);
  const double lopz2 = lopz * lopz;
  const double Li2mz = apfel::dilog(-z);
  const double Li2z  = apfel::dilog(z);
  const double Li3opzh    = apfel::wgplg(2,1,(1 + z)/2.);
  const double Li3omz     = apfel::wgplg(2,1,1 - z);
  const double Li3omzh    = apfel::wgplg(2,1,0.5 - z/2.);
  const double Li3mz      = apfel::wgplg(2,1,-z);
  const double Li3z       = apfel::wgplg(2,1,z);
  const double Li3opzi    = apfel::wgplg(2,1,opzi);
  const double Li3omzopzi = apfel::wgplg(2,1,-((-1 + z)*opzi));
  return 9.694444444444445 + (253*z)/72. + (3*apfel::zeta3)/2. + (21*z*apfel::zeta3)/4. -  
   5*Li3omz - (5*z*Li3omz)/2. + 2*Li3omzh +  
   z*Li3omzh - Li3mz - (z*Li3mz)/2. + 5*Li3z -  
   (19*z*Li3z)/2. + 2*Li3opzh + z*Li3opzh -  
   2*Li3opzi - z*Li3opzi +  
   2*Li3omzopzi + z*Li3omzopzi -  
   (73*lomz)/12. + (23*z*lomz)/6. - (157*lz)/8. +  
   (7*z*lz)/8. - (3*z*lomz*lz)/2. -  
   2*l2*lomz*lopz - z*l2*lomz*lopz -  
   (z*lz*lopz)/2. - 2*lomz*lz*lopz -  
   z*lomz*lz*lopz - (3*NCi2)/2. -  
   (13*z*NCi2)/16. - (3*apfel::zeta3*NCi2)/2. +  
   (3*z*apfel::zeta3*NCi2)/4. + Li3omz*NCi2 -  
   (z*Li3omz*NCi2)/2. - 4*Li3z*NCi2 +  
   2*z*Li3z*NCi2 + (9*lomz*NCi2)/4. -  
   (5*z*lomz*NCi2)/4. + (67*lz*NCi2)/16. -  
   (45*z*lz*NCi2)/16. - (295*NC2)/36. -  
   (389*z*NC2)/144. - 6*z*apfel::zeta3*NC2 + 4*Li3omz*NC2 +  
   3*z*Li3omz*NC2 - 2*Li3omzh*NC2 -  
   z*Li3omzh*NC2 + Li3mz*NC2 +  
   (z*Li3mz*NC2)/2. - Li3z*NC2 + (15*z*Li3z*NC2)/2. -  
   2*Li3opzh*NC2 - z*Li3opzh*NC2 +  
   2*Li3opzi*NC2 + z*Li3opzi*NC2 -  
   2*Li3omzopzi*NC2 -  
   z*Li3omzopzi*NC2 +  
   (23*lomz*NC2)/6. - (31*z*lomz*NC2)/12. +  
   (247*lz*NC2)/16. + (31*z*lz*NC2)/16. +  
   (3*z*lomz*lz*NC2)/2. +  
   2*l2*lomz*lopz*NC2 +  
   z*l2*lomz*lopz*NC2 +  
   (z*lz*lopz*NC2)/2. +  
   2*lomz*lz*lopz*NC2 +  
   z*lomz*lz*lopz*NC2 + (7*pi2)/6. -  
   (5*z*pi2)/24. + (l2*pi2)/3. + (z*l2*pi2)/6. -  
   (lopz*pi2)/3. - (z*lopz*pi2)/6. -  
   (lomz*pi2)/2. + (7*z*lomz*pi2)/12. -  
   (lz*pi2)/6. + (7*z*lz*pi2)/12. -  
   (lopz*pi2)/6. - (z*lopz*pi2)/12. -  
   (NCi2*pi2)/3. + (z*NCi2*pi2)/16. +  
   (lomz*NCi2*pi2)/4. -  
   (z*lomz*NCi2*pi2)/8. +  
   (lz*NCi2*pi2)/4. - (z*lz*NCi2*pi2)/8. -  
   (5*NC2*pi2)/6. + (7*z*NC2*pi2)/48. -  
   (l2*NC2*pi2)/3. - (z*l2*NC2*pi2)/6. +  
   (lopz*NC2*pi2)/3. +  
   (z*lopz*NC2*pi2)/6. +  
   (lomz*NC2*pi2)/4. -  
   (11*z*lomz*NC2*pi2)/24. -  
   (lz*NC2*pi2)/12. - (11*z*lz*NC2*pi2)/24. +  
   (lopz*NC2*pi2)/6. +  
   (z*lopz*NC2*pi2)/12. - (1169*zi)/108. +  
   (9*apfel::zeta3*zi)/2. - 2*Li3omz*zi +  
   2*Li3omzh*zi - Li3mz*zi - 13*Li3z*zi +  
   2*Li3opzh*zi - 2*Li3opzi*zi +  
   2*Li3omzopzi*zi +  
   (34*lomz*zi)/9. + (5*lz*zi)/18. -  
   2*l2*lomz*lopz*zi -  
   2*lomz*lz*lopz*zi +  
   (5*apfel::zeta3*NCi2*zi)/2. -  
   (3*Li3omz*NCi2*zi)/2. + 3*Li3z*NCi2*zi -  
   (9*lomz*NCi2*zi)/4. -  
   (9*lz*NCi2*zi)/4. + (1169*NC2*zi)/108. -  
   7*apfel::zeta3*NC2*zi + (7*Li3omz*NC2*zi)/2. -  
   2*Li3omzh*NC2*zi + Li3mz*NC2*zi +  
   10*Li3z*NC2*zi - 2*Li3opzh*NC2*zi +  
   2*Li3opzi*NC2*zi -  
   2*Li3omzopzi*NC2*zi -  
   (55*lomz*NC2*zi)/36. +  
   (71*lz*NC2*zi)/36. +  
   2*l2*lomz*lopz*NC2*zi +  
   2*lomz*lz*lopz*NC2*zi +  
   (l2*pi2*zi)/3. - (lopz*pi2*zi)/3. +  
   (2*lomz*pi2*zi)/3. +  
   (2*lz*pi2*zi)/3. - (lopz*pi2*zi)/6. +  
   (NCi2*pi2*zi)/8. -  
   (lomz*NCi2*pi2*zi)/6. -  
   (lz*NCi2*pi2*zi)/6. -  
   (NC2*pi2*zi)/8. -  
   (l2*NC2*pi2*zi)/3. +  
   (lopz*NC2*pi2*zi)/3. -  
   (lomz*NC2*pi2*zi)/2. -  
   (lz*NC2*pi2*zi)/2. +  
   (lopz*NC2*pi2*zi)/6. +  
   (Li2mz*(-1 + NC2)*(z + 2*lomz*(2 + z + 2*zi) +  
        lz*(2 + z + 2*zi)))/2. + (269*z2)/108. +  
   (13*lomz*z2)/18. - (31*lz*z2)/18. -  
   (269*NC2*z2)/108. - (13*lomz*NC2*z2)/18. +  
   (31*lz*NC2*z2)/18. -  
   (Li2z*(-12 + 6*z*NCi2 + 12*NC2 - 6*z*NC2 +  
        98*zi - 9*NCi2*zi - 89*NC2*zi +  
        6*lz*(-5 + 3*NCi2 + 2*NC2)*(-2 + z + 2*zi) +  
        6*lomz*(3*(2 + z) + NCi2*(-4 + 2*z + 5*zi) -  
           NC2*(2 + 5*z + 5*zi)) - 8*z2 +  
        8*NC2*z2))/12. + lomz*l22 +  
   (z*lomz*l22)/2. + lopz*l22 +  
   (z*lopz*l22)/2. - lomz*NC2*l22 -  
   (z*lomz*NC2*l22)/2. -  
   lopz*NC2*l22 -  
   (z*lopz*NC2*l22)/2. +  
   lomz*zi*l22 + lopz*zi*l22 -  
   lomz*NC2*zi*l22 -  
   lopz*NC2*zi*l22 - (2*l23)/3. -  
   (z*l23)/3. + (2*NC2*l23)/3. +  
   (z*NC2*l23)/3. - (2*zi*l23)/3. +  
   (2*NC2*zi*l23)/3. - 4*lomz2 +  
   (z*lomz2)/8. - (3*z*lz*lomz2)/2. +  
   NCi2*lomz2 - (z*NCi2*lomz2)/16. +  
   (3*lz*NCi2*lomz2)/4. -  
   (3*z*lz*NCi2*lomz2)/8. +  
   3*NC2*lomz2 - (z*NC2*lomz2)/16. -  
   (3*lz*NC2*lomz2)/4. +  
   (15*z*lz*NC2*lomz2)/8. +  
   (49*zi*lomz2)/12. -  
   (3*lz*zi*lomz2)/2. -  
   (3*NCi2*zi*lomz2)/4. -  
   lz*NCi2*zi*lomz2 -  
   (10*NC2*zi*lomz2)/3. +  
   (5*lz*NC2*zi*lomz2)/2. -  
   (z2*lomz2)/3. +  
   (NC2*z2*lomz2)/3. + lomz3 -  
   (z*lomz3)/2. - (5*NCi2*lomz3)/12. +  
   (5*z*NCi2*lomz3)/24. -  
   (7*NC2*lomz3)/12. +  
   (7*z*NC2*lomz3)/24. - zi*lomz3 +  
   (5*NCi2*zi*lomz3)/12. +  
   (7*NC2*zi*lomz3)/12. + 3*lz2 -  
   (3*z*lz2)/16. - (lomz*lz2)/2. +  
   (z*lomz*lz2)/4. - (3*lopz*lz2)/2. -  
   (3*z*lopz*lz2)/4. - (NCi2*lz2)/2. +  
   (13*z*NCi2*lz2)/32. +  
   lomz*NCi2*lz2 -  
   (z*lomz*NCi2*lz2)/2. -  
   (5*NC2*lz2)/2. - (7*z*NC2*lz2)/32. -  
   (lomz*NC2*lz2)/2. +  
   (z*lomz*NC2*lz2)/4. +  
   (3*lopz*NC2*lz2)/2. +  
   (3*z*lopz*NC2*lz2)/4. +  
   (10*zi*lz2)/3. +  
   (lomz*zi*lz2)/2. -  
   (3*lopz*zi*lz2)/2. -  
   lomz*NCi2*zi*lz2 -  
   (10*NC2*zi*lz2)/3. +  
   (lomz*NC2*zi*lz2)/2. +  
   (3*lopz*NC2*zi*lz2)/2. +  
   (5*lz3)/12. + (35*z*lz3)/24. +  
   (5*NCi2*lz3)/24. - (5*z*NCi2*lz3)/48. -  
   (5*NC2*lz3)/8. - (65*z*NC2*lz3)/48. +  
   (5*zi*lz3)/3. -  
   (5*NC2*zi*lz3)/3. + lomz*lopz2 +  
   (z*lomz*lopz2)/2. -  
   lomz*NC2*lopz2 -  
   (z*lomz*NC2*lopz2)/2. +  
   lomz*zi*lopz2 -  
   lomz*NC2*zi*lopz2;
};
double C2TG2G_Lx_Rz(int const&, double const&)
{
  return 0;
};
double C2TG2Q_Lx_Rz(int const&, double const&)
{
  return 0;
};
double C2TQ2QNS_Lx_Rz(int const& nf, double const& z)
{
  const double NC   = apfel::NC;
  const double NC2  = NC * NC;
  const double NCi  = 1. / NC;
  const double NCi2 = NCi * NCi;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double z2 = z * z;
  const double zi  = 1. / z;
  const double omzi  = 1. / ( 1 - z );
  const double lz  = log(z);
  const double lz2 = lz * lz;
  const double lz3 = lz * lz2;
  const double lomz  = log(1 - z);
  const double lomz2 = lomz * lomz;
  const double lomz3 = lomz * lomz2;
  const double Li2z  = apfel::dilog(z);
  const double Li3z   = apfel::wgplg(2,1,z);
  const double Li3omz = apfel::wgplg(2,1,1 - z);
  return (-448 - 264*NC - 76*NC*nf - 1168*z - 96*NC*z - 148*NC*nf*z -  
     1404*apfel::zeta3 + 432*NC*apfel::zeta3 - 1404*z*apfel::zeta3 + 432*NC*z*apfel::zeta3 +  
     3024*Li3z - 432*NC*Li3z + 3024*z*Li3z - 432*NC*z*Li3z -  
     672*lomz - 720*NC*lomz + 48*NC*nf*lomz -  
     1176*z*lomz + 504*NC*z*lomz + 192*NC*nf*z*lomz -  
     2208*lz - 1080*NC*lz - 12*NC*nf*lz - 2172*z*lz -  
     972*NC*z*lz + 132*NC*nf*z*lz - 648*lomz*lz +  
     648*z*lomz*lz + 54*NCi2 - 54*z*NCi2 +  
     864*apfel::zeta3*NCi2 + 864*z*apfel::zeta3*NCi2 - 1296*Li3z*NCi2 -  
     1296*z*Li3z*NCi2 + 270*lomz*NCi2 +  
     1512*z*lomz*NCi2 + 810*lz*NCi2 +  
     1566*z*lz*NCi2 + 324*lomz*lz*NCi2 -  
     324*z*lomz*lz*NCi2 + 264*NCi + 76*nf*NCi +  
     96*z*NCi + 148*nf*z*NCi - 432*apfel::zeta3*NCi -  
     432*z*apfel::zeta3*NCi + 432*Li3z*NCi +  
     432*z*Li3z*NCi + 720*lomz*NCi -  
     48*nf*lomz*NCi - 504*z*lomz*NCi -  
     192*nf*z*lomz*NCi + 1080*lz*NCi +  
     12*nf*lz*NCi + 972*z*lz*NCi -  
     132*nf*z*lz*NCi + 394*NC2 + 1222*z*NC2 +  
     540*apfel::zeta3*NC2 + 540*z*apfel::zeta3*NC2 - 1728*Li3z*NC2 -  
     1728*z*Li3z*NC2 + 402*lomz*NC2 -  
     336*z*lomz*NC2 + 1398*lz*NC2 +  
     606*z*lz*NC2 + 324*lomz*lz*NC2 -  
     324*z*lomz*lz*NC2 + 300*pi2 - 36*NC*pi2 +  
     12*NC*nf*pi2 + 84*z*pi2 - 54*NC*z*pi2 +  
     12*NC*nf*z*pi2 - 216*lomz*pi2 +  
     36*NC*lomz*pi2 - 216*z*lomz*pi2 +  
     36*NC*z*lomz*pi2 - 180*lz*pi2 +  
     36*NC*lz*pi2 - 180*z*lz*pi2 +  
     36*NC*z*lz*pi2 - 126*NCi2*pi2 -  
     36*z*NCi2*pi2 + 72*lomz*NCi2*pi2 +  
     72*z*lomz*NCi2*pi2 + 72*lz*NCi2*pi2 +  
     72*z*lz*NCi2*pi2 + 36*NCi*pi2 -  
     12*nf*NCi*pi2 + 54*z*NCi*pi2 -  
     12*nf*z*NCi*pi2 - 36*lomz*NCi*pi2 -  
     36*z*lomz*NCi*pi2 - 36*lz*NCi*pi2 -  
     36*z*lz*NCi*pi2 - 174*NC2*pi2 -  
     48*z*NC2*pi2 + 144*lomz*NC2*pi2 +  
     144*z*lomz*NC2*pi2 + 108*lz*NC2*pi2 +  
     108*z*lz*NC2*pi2 + 5184*apfel::zeta3*omzi -  
     5184*Li3z*omzi + 3408*lz*omzi -  
     120*NC*nf*lz*omzi - 2160*apfel::zeta3*NCi2*omzi +  
     2160*Li3z*NCi2*omzi -  
     2052*lz*NCi2*omzi +  
     120*nf*lz*NCi*omzi -  
     3024*apfel::zeta3*NC2*omzi +  
     3024*Li3z*NC2*omzi -  
     1356*lz*NC2*omzi - 108*pi2*omzi +  
     288*lz*pi2*omzi +  
     54*NCi2*pi2*omzi +  
     36*lomz*NCi2*pi2*omzi -  
     108*lz*NCi2*pi2*omzi +  
     54*NC2*pi2*omzi -  
     36*lomz*NC2*pi2*omzi -  
     180*lz*NC2*pi2*omzi -  
     216*Li3omz*(-1 + NC - z + NC*z - (1 + z)*NCi + 2*NC2 +  
        2*z*NC2 - NCi2*(1 + z - 3*omzi) -  
        3*NC2*omzi) - 214*NC*zi -  
     84*NC*lomz*zi - 84*NC*lz*zi +  
     214*NCi*zi + 84*lomz*NCi*zi +  
     84*lz*NCi*zi + 574*NC*z2 +  
     300*NC*lomz*z2 - 372*NC*lz*z2 -  
     574*NCi*z2 - 300*lomz*NCi*z2 +  
     372*lz*NCi*z2 -  
     36*Li2z*(48 - 3*NC - 6*z - 12*NC*z - 27*NCi2 + 3*NCi +  
        12*z*NCi - 21*NC2 + 6*z*NC2 -  
        6*lz*(-5 + 2*NCi2 + 3*NC2)* 
         (1 + z - 2*omzi) - 18*omzi +  
        9*NCi2*omzi + 9*NC2*omzi +  
        6*lomz*(-1 + NC - z + NC*z - (1 + z)*NCi +  
           NC2*(1 + z - omzi) + NCi2*omzi) +  
        4*NC*zi - 4*NCi*zi - 4*NC*z2 +  
        4*NCi*z2) - 198*lomz2 +  
     54*NC*lomz2 - 36*NC*nf*lomz2 -  
     198*z*lomz2 - 54*NC*z*lomz2 -  
     36*NC*nf*z*lomz2 + 432*lz*lomz2 -  
     108*NC*lz*lomz2 + 432*z*lz*lomz2 -  
     108*NC*z*lz*lomz2 -  
     108*lz*NCi2*lomz2 -  
     108*z*lz*NCi2*lomz2 -  
     54*NCi*lomz2 + 36*nf*NCi*lomz2 +  
     54*z*NCi*lomz2 +  
     36*nf*z*NCi*lomz2 +  
     108*lz*NCi*lomz2 +  
     108*z*lz*NCi*lomz2 +  
     198*NC2*lomz2 + 198*z*NC2*lomz2 -  
     324*lz*NC2*lomz2 -  
     324*z*lz*NC2*lomz2 -  
     648*lz*omzi*lomz2 +  
     108*lz*NCi2*omzi*lomz2 +  
     540*lz*NC2*omzi*lomz2 +  
     72*NC*zi*lomz2 -  
     72*NCi*zi*lomz2 -  
     72*NC*z2*lomz2 +  
     72*NCi*z2*lomz2 + 216*lomz3 +  
     216*z*lomz3 - 108*NCi2*lomz3 -  
     108*z*NCi2*lomz3 - 108*NC2*lomz3 -  
     108*z*NC2*lomz3 + 1017*lz2 -  
     27*NC*lz2 + 18*NC*nf*lz2 + 369*z*lz2 -  
     135*NC*z*lz2 + 18*NC*nf*z*lz2 -  
     108*lomz*lz2 - 108*z*lomz*lz2 -  
     567*NCi2*lz2 - 189*z*NCi2*lz2 +  
     54*lomz*NCi2*lz2 +  
     54*z*lomz*NCi2*lz2 + 27*NCi*lz2 -  
     18*nf*NCi*lz2 + 135*z*NCi*lz2 -  
     18*nf*z*NCi*lz2 - 450*NC2*lz2 -  
     180*z*NC2*lz2 + 54*lomz*NC2*lz2 +  
     54*z*lomz*NC2*lz2 -  
     684*omzi*lz2 - 36*NC*nf*omzi*lz2 +  
     216*lomz*omzi*lz2 +  
     243*NCi2*omzi*lz2 -  
     108*lomz*NCi2*omzi*lz2 +  
     36*nf*NCi*omzi*lz2 +  
     441*NC2*omzi*lz2 -  
     108*lomz*NC2*omzi*lz2 +  
     72*NC*zi*lz2 - 72*NCi*zi*lz2 -  
     360*lz3 + 90*NC*lz3 - 360*z*lz3 +  
     90*NC*z*lz3 + 225*NCi2*lz3 +  
     225*z*NCi2*lz3 - 90*NCi*lz3 -  
     90*z*NCi*lz3 + 135*NC2*lz3 +  
     135*z*NC2*lz3 + 540*omzi*lz3 -  
     360*NCi2*omzi*lz3 -  
     180*NC2*omzi*lz3)/432.;
};
double C2TQ2QPS_Lx_Rz(int const&, double const&)
{
  return 0;
};
double C2TQ2QB_Lx_Rz(int const&, double const& z)
{
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double NCi2 = NCi * NCi;
  const double l2  = log(2);
  const double l22 = l2 * l2;
  const double l23 = l2 * l22;
  const double pi = M_PI;
  const double pi2 = pi * pi;
  const double z2 = z * z;
  const double zi = 1. / z;
  const double opzi = 1. / ( 1 + z );
  const double lz  = log(z);
  const double lz2 = lz * lz;
  const double lz3 = lz * lz2;
  const double lomz  = log(1 - z);
  const double lomz2 = lomz * lomz;
  const double lopz  = log(1 + z);
  const double lopz2 = lopz * lopz;
  const double Li2z  = apfel::dilog(z);
  const double Li2mz = apfel::dilog(-z);
  const double Li3z       = apfel::wgplg(2,1,z);
  const double Li3mz      = apfel::wgplg(2,1,-z);
  const double Li3omz     = apfel::wgplg(2,1,1 - z);
  const double Li3omzh    = apfel::wgplg(2,1,0.5 - z/2.);
  const double Li3opzh    = apfel::wgplg(2,1,(1 + z)/2.);
  const double Li3opzi    = apfel::wgplg(2,1,opzi);
  const double Li3omzopzi = apfel::wgplg(2,1,-((-1 + z)*opzi));
  return 0.125 - (11*NC)/18. - z/8. - (2*NC*z)/9. + NC*apfel::zeta3 + NC*z*apfel::zeta3 -  
   Li3omz - (NC*Li3omz)/2. + z*Li3omz - (NC*z*Li3omz)/2. +  
   Li3omzh - z*Li3omzh - Li3mz/2. + (z*Li3mz)/2. -  
   Li3z/2. - NC*Li3z + (z*Li3z)/2. - NC*z*Li3z + Li3opzh -  
   z*Li3opzh - Li3opzi + z*Li3opzi +  
   Li3omzopzi - z*Li3omzopzi -  
   lomz - (5*NC*lomz)/3. + z*lomz +  
   (7*NC*z*lomz)/6. - (7*lz)/8. - (5*NC*lz)/2. +  
   (z*lz)/8. - (9*NC*z*lz)/4. - l2*lomz*lopz +  
   z*l2*lomz*lopz + (lz*lopz)/2. +  
   (z*lz*lopz)/2. - lomz*lz*lopz +  
   z*lomz*lz*lopz + 2*Li3omz*opzi -  
   2*Li3omzh*opzi + Li3mz*opzi +  
   Li3z*opzi - 2*Li3opzh*opzi +  
   2*Li3opzi*opzi -  
   2*Li3omzopzi*opzi +  
   2*l2*lomz*lopz*opzi +  
   2*lomz*lz*lopz*opzi -  
   (Li2mz*(1 + z + 2*lomz*(-1 + z + 2*opzi) +  
        lz*(-1 + z + 2*opzi))*(-1 + NCi2))/2. -  
   NCi2/8. + (z*NCi2)/8. + Li3omz*NCi2 -  
   z*Li3omz*NCi2 - Li3omzh*NCi2 +  
   z*Li3omzh*NCi2 + (Li3mz*NCi2)/2. -  
   (z*Li3mz*NCi2)/2. + (Li3z*NCi2)/2. -  
   (z*Li3z*NCi2)/2. - Li3opzh*NCi2 +  
   z*Li3opzh*NCi2 + Li3opzi*NCi2 -  
   z*Li3opzi*NCi2 -  
   Li3omzopzi*NCi2 +  
   z*Li3omzopzi*NCi2 + lomz*NCi2 -  
   z*lomz*NCi2 + (7*lz*NCi2)/8. -  
   (z*lz*NCi2)/8. + l2*lomz*lopz*NCi2 -  
   z*l2*lomz*lopz*NCi2 -  
   (lz*lopz*NCi2)/2. - (z*lz*lopz*NCi2)/2. +  
   lomz*lz*lopz*NCi2 -  
   z*lomz*lz*lopz*NCi2 -  
   2*Li3omz*opzi*NCi2 +  
   2*Li3omzh*opzi*NCi2 -  
   Li3mz*opzi*NCi2 - Li3z*opzi*NCi2 +  
   2*Li3opzh*opzi*NCi2 -  
   2*Li3opzi*opzi*NCi2 +  
   2*Li3omzopzi*opzi*NCi2 -  
   2*l2*lomz*lopz*opzi*NCi2 -  
   2*lomz*lz*lopz*opzi*NCi2 +  
   (11*NCi)/18. + (2*z*NCi)/9. - apfel::zeta3*NCi -  
   z*apfel::zeta3*NCi + (Li3omz*NCi)/2. +  
   (z*Li3omz*NCi)/2. + Li3z*NCi + z*Li3z*NCi +  
   (5*lomz*NCi)/3. - (7*z*lomz*NCi)/6. +  
   (5*lz*NCi)/2. + (9*z*lz*NCi)/4. - pi2/24. -  
   (NC*pi2)/12. - (z*pi2)/24. - (NC*z*pi2)/8. +  
   (l2*pi2)/6. - (z*l2*pi2)/6. -  
   (lopz*pi2)/6. + (z*lopz*pi2)/6. -  
   (lomz*pi2)/12. + (NC*lomz*pi2)/12. +  
   (z*lomz*pi2)/12. + (NC*z*lomz*pi2)/12. +  
   (NC*lz*pi2)/12. + (NC*z*lz*pi2)/12. -  
   (lopz*pi2)/12. + (z*lopz*pi2)/12. -  
   (l2*opzi*pi2)/3. +  
   (lopz*opzi*pi2)/3. +  
   (lomz*opzi*pi2)/6. +  
   (lopz*opzi*pi2)/6. + (NCi2*pi2)/24. +  
   (z*NCi2*pi2)/24. - (l2*NCi2*pi2)/6. +  
   (z*l2*NCi2*pi2)/6. +  
   (lopz*NCi2*pi2)/6. -  
   (z*lopz*NCi2*pi2)/6. +  
   (lomz*NCi2*pi2)/12. -  
   (z*lomz*NCi2*pi2)/12. +  
   (lopz*NCi2*pi2)/12. -  
   (z*lopz*NCi2*pi2)/12. +  
   (l2*opzi*NCi2*pi2)/3. -  
   (lopz*opzi*NCi2*pi2)/3. -  
   (lomz*opzi*NCi2*pi2)/6. -  
   (lopz*opzi*NCi2*pi2)/6. +  
   (NCi*pi2)/12. + (z*NCi*pi2)/8. -  
   (lomz*NCi*pi2)/12. -  
   (z*lomz*NCi*pi2)/12. -  
   (lz*NCi*pi2)/12. - (z*lz*NCi*pi2)/12. -  
   (107*NC*zi)/216. - (7*NC*lomz*zi)/36. -  
   (7*NC*lz*zi)/36. + (107*NCi*zi)/216. +  
   (7*lomz*NCi*zi)/36. +  
   (7*lz*NCi*zi)/36. + (287*NC*z2)/216. +  
   (25*NC*lomz*z2)/36. - (31*NC*lz*z2)/36. -  
   (287*NCi*z2)/216. - (25*lomz*NCi*z2)/36. +  
   (31*lz*NCi*z2)/36. -  
   (Li2z*(-6 - 3*NC - 6*z - 12*NC*z + 6*(1 + z)*NCi2 +  
        6*(1 + z)*lomz*(NC - NCi) + 3*NCi +  
        12*z*NCi + 4*NC*zi - 4*NCi*zi -  
        4*NC*z2 + 4*NCi*z2))/12. +  
   (lomz*l22)/2. - (z*lomz*l22)/2. +  
   (lopz*l22)/2. - (z*lopz*l22)/2. -  
   lomz*opzi*l22 -  
   lopz*opzi*l22 -  
   (lomz*NCi2*l22)/2. +  
   (z*lomz*NCi2*l22)/2. -  
   (lopz*NCi2*l22)/2. +  
   (z*lopz*NCi2*l22)/2. +  
   lomz*opzi*NCi2*l22 +  
   lopz*opzi*NCi2*l22 - l23/3. +  
   (z*l23)/3. + (2*opzi*l23)/3. +  
   (NCi2*l23)/3. - (z*NCi2*l23)/3. -  
   (2*opzi*NCi2*l23)/3. +  
   (NC*lomz2)/8. - (NC*z*lomz2)/8. -  
   (NC*lz*lomz2)/4. - (NC*z*lz*lomz2)/4. -  
   (NCi*lomz2)/8. + (z*NCi*lomz2)/8. +  
   (lz*NCi*lomz2)/4. +  
   (z*lz*NCi*lomz2)/4. +  
   (NC*zi*lomz2)/6. -  
   (NCi*zi*lomz2)/6. -  
   (NC*z2*lomz2)/6. +  
   (NCi*z2*lomz2)/6. - lz2/2. -  
   (NC*lz2)/16. - (z*lz2)/2. -  
   (5*NC*z*lz2)/16. - (3*lopz*lz2)/4. +  
   (3*z*lopz*lz2)/4. +  
   (3*lopz*opzi*lz2)/2. +  
   (NCi2*lz2)/2. + (z*NCi2*lz2)/2. +  
   (3*lopz*NCi2*lz2)/4. -  
   (3*z*lopz*NCi2*lz2)/4. -  
   (3*lopz*opzi*NCi2*lz2)/2. +  
   (NCi*lz2)/16. + (5*z*NCi*lz2)/16. +  
   (NC*zi*lz2)/6. -  
   (NCi*zi*lz2)/6. + (5*lz3)/24. +  
   (5*NC*lz3)/24. - (5*z*lz3)/24. +  
   (5*NC*z*lz3)/24. - (5*opzi*lz3)/12. -  
   (5*NCi2*lz3)/24. + (5*z*NCi2*lz3)/24. +  
   (5*opzi*NCi2*lz3)/12. -  
   (5*NCi*lz3)/24. - (5*z*NCi*lz3)/24. +  
   (lomz*lopz2)/2. - (z*lomz*lopz2)/2. -  
   lomz*opzi*lopz2 -  
   (lomz*NCi2*lopz2)/2. +  
   (z*lomz*NCi2*lopz2)/2. +  
   lomz*opzi*NCi2*lopz2;
};
double C2TQ2QP1_Lx_Rz(int const&, double const& z)
{
  const double NC  = apfel::NC;
  const double NCi = 1. / NC;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double z2  = z * z;
  const double zi  = 1. / z;
  const double lz  = log(z);
  const double lz2 = lz * lz;
  const double lz3 = lz * lz2;
  const double lomz  = log(1 - z);
  const double lomz2 = lomz * lomz;
  const double Li2z  = apfel::dilog(z);
  const double Li3omz = apfel::wgplg(2,1,1 - z);
  const double Li3z   = apfel::wgplg(2,1,z);
  return ((NC - NCi)*(-264 - 96*z + 432*apfel::zeta3 + 432*z*apfel::zeta3 -  
       216*(1 + z)*Li3omz - 432*Li3z - 432*z*Li3z - 720*lomz +  
       504*z*lomz - 1080*lz - 972*z*lz - 36*pi2 -  
       54*z*pi2 + 36*lomz*pi2 + 36*z*lomz*pi2 +  
       36*lz*pi2 + 36*z*lz*pi2 - 214*zi -  
       84*lomz*zi - 84*lz*zi -  
       36*Li2z*(-3 - 12*z + 6*(1 + z)*lomz + 4*zi -  
          4*z2) + 574*z2 + 300*lomz*z2 -  
       372*lz*z2 + 54*lomz2 - 54*z*lomz2 -  
       108*lz*lomz2 - 108*z*lz*lomz2 +  
       72*zi*lomz2 - 72*z2*lomz2 -  
       27*lz2 - 135*z*lz2 + 72*zi*lz2 +  
       90*lz3 + 90*z*lz3))/432.;
};
double C2TQ2QP2_Lx_Rz(int const&, double const&)
{
  return 0;
};
double C2TQ2QP3_Lx_Rz(int const&, double const&)
{
  return 0;
};
//------------------------------------------------------
double C1LQ2Q_Rx_Lz(int const&, double const&)
{
  return 0;
};
double C1LQ2G_Rx_Lz(int const&, double const&)
{
  return 0;
};
double C1LG2Q_Rx_Lz(int const&, double const&)
{
  return 0;
};
double C1TQ2Q_Rx_Lz(int const&, double const& x)
{
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double lx  = log(x);
  const double lomx  = log(1 - x);
  const double omxi = 1. / ( 1 - x );
  return -0.5*((NC - NCi)* 
     (-1 + x + (1 + x)*lomx - lx*(1 + x - 2*omxi)));
};
double C1TQ2G_Rx_Lz(int const&, double const&)
{
  return 0;
};
double C1TG2Q_Rx_Lz(int const&, double const& x)
{
  const double x2 = x * x;
  const double lx  = log(x);
  const double lomx  = log(1 - x);
  return x + lx*(-0.5 + x - x2) - x2 +  
   lomx*(0.5 - x + x2);
};
double C2LQ2G_Rx_Lz(int const&, double const&)
{
  return 0;
};
double C2LG2G_Rx_Lz(int const&, double const&)
{
  return 0;
};
double C2LG2Q_Rx_Lz(int const&, double const&)
{
  return 0;
};
double C2LQ2QNS_Rx_Lz(int const&, double const&)
{
  return 0;
};
double C2LQ2QPS_Rx_Lz(int const&, double const&)
{
  return 0;
};
double C2LQ2QB_Rx_Lz(int const&, double const&)
{
  return 0;
};
double C2LQ2QP1_Rx_Lz(int const&, double const&)
{
  return 0;
};
double C2LQ2QP2_Rx_Lz(int const&, double const&)
{
  return 0;
};
double C2LQ2QP3_Rx_Lz(int const&, double const&)
{
  return 0;
};
double C2TQ2G_Rx_Lz(int const&, double const&)
{
  return 0;
};
double C2TG2G_Rx_Lz(int const&, double const&)
{
  return 0;
};
double C2TG2Q_Rx_Lz(int const& nf, double const& x)
{
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double l2  = log(2);
  const double l22 = l2 * l2;
  const double l23 = l2 * l22;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double x2 = x * x;
  const double xi  = 1. / x;
  const double opxi = 1. / ( 1 + x );
  const double lx  = log(x);
  const double lx2 = lx * lx;
  const double lx3 = lx * lx2;
  const double lomx  = log(1 - x);
  const double lomx2 = lomx * lomx;
  const double lomx3 = lomx * lomx2;
  const double lopx  = log(1 + x);
  const double lopx2 = lopx * lopx;
  const double lopx3 = lopx * lopx2;
  const double Li2x  = apfel::dilog(x);
  const double Li2mx = apfel::dilog(-x);
  const double Li3x       = apfel::wgplg(2,1,x);
  const double Li3omx     = apfel::wgplg(2,1,1 - x);
  const double Li3omxh    = apfel::wgplg(2,1,0.5 - x/2.);
  const double Li3opxh    = apfel::wgplg(2,1,(1 + x)/2.);
  const double Li3omxopxi = apfel::wgplg(2,1,-((-1 + x)*opxi));
  return (-43*NC)/18. + (1049*NC*x)/144. + (9*NC*apfel::zeta3)/4. -  
   (13*NC*x*apfel::zeta3)/2. - (3*NC*Li3omx)/2. - 7*NC*x*Li3omx +  
   NC*Li3omxh + 2*NC*x*Li3omxh - (9*NC*Li3x)/4. +  
   (5*NC*x*Li3x)/2. + NC*Li3opxh + 2*NC*x*Li3opxh +  
   NC*Li3omxopxi +  
   2*NC*x*Li3omxopxi - (59*NC*lomx)/24. +  
   (16*NC*x*lomx)/3. + (29*NC*lx)/6. + (143*NC*x*lx)/48. +  
   (3*NC*lomx*lx)/8. - (15*NC*x*lomx*lx)/2. -  
   NC*l2*lomx*lopx - 2*NC*x*l2*lomx*lopx +  
   NC*x*lx*lopx - NC*lomx*lx*lopx -  
   2*NC*x*lomx*lx*lopx + (17*NCi)/8. -  
   (81*x*NCi)/16. - 4*apfel::zeta3*NCi + 8*x*apfel::zeta3*NCi +  
   (3*Li3omx*NCi)/2. - 3*x*Li3omx*NCi +  
   (5*Li3x*NCi)/4. - (5*x*Li3x*NCi)/2. +  
   (3*lomx*NCi)/8. + (7*x*lomx*NCi)/4. -  
   (3*lx*NCi)/4. - (37*x*lx*NCi)/16. -  
   (7*lomx*lx*NCi)/8. +  
   (7*x*lomx*lx*NCi)/2. + (7*NC*pi2)/24. -  
   (NC*x*pi2)/4. + (NC*l2*pi2)/6. +  
   (NC*x*l2*pi2)/3. - (NC*lopx*pi2)/4. -  
   (NC*x*lopx*pi2)/2. - (NC*lomx*pi2)/6. +  
   NC*x*lomx*pi2 + (NC*lx*pi2)/12. -  
   (3*NC*x*lx*pi2)/2. + (NC*lopx*pi2)/12. +  
   (NC*x*lopx*pi2)/6. - (NCi*pi2)/8. +  
   (5*x*NCi*pi2)/12. + (lomx*NCi*pi2)/12. -  
   (x*lomx*NCi*pi2)/6. -  
   (lx*NCi*pi2)/4. + (x*lx*NCi*pi2)/2. +  
   (52*NC*xi)/27. + (13*NC*lomx*xi)/9. -  
   (2*NC*lomx*lx*xi)/3. - (3433*NC*x2)/432. +  
   (9*NC*apfel::zeta3*x2)/2. - (NC*Li3omx*x2)/2. +  
   2*NC*Li3omxh*x2 - (9*NC*Li3x*x2)/2. +  
   2*NC*Li3opxh*x2 +  
   2*NC*Li3omxopxi*x2 -  
   (383*NC*lomx*x2)/72. + (583*NC*lx*x2)/72. +  
   (49*NC*lomx*lx*x2)/6. -  
   2*NC*l2*lomx*lopx*x2 +  
   NC*lx*lopx*x2 -  
   2*NC*lomx*lx*lopx*x2 +  
   (57*NCi*x2)/16. - 8*apfel::zeta3*NCi*x2 +  
   (5*Li3omx*NCi*x2)/2. +  
   (5*Li3x*NCi*x2)/2. -  
   (7*lomx*NCi*x2)/8. +  
   (7*lx*NCi*x2)/8. -  
   3*lomx*lx*NCi*x2 + (29*NC*pi2*x2)/12. +  
   (NC*l2*pi2*x2)/3. -  
   (NC*lopx*pi2*x2)/2. -  
   (3*NC*lomx*pi2*x2)/4. + NC*lx*pi2*x2 +  
   (NC*lopx*pi2*x2)/6. -  
   (NCi*pi2*x2)/3. +  
   (lomx*NCi*pi2*x2)/4. -  
   (2*lx*NCi*pi2*x2)/3. +  
   NC*Li2mx*(x + x2 - lomx*(1 + 2*x + 2*x2) +  
      lx*(1 + 2*x + 2*x2)) -  
   (Li2x*(33*NC + 96*NC*x + 3*NCi + 16*NC*xi +  
        12*lomx*(NCi*(-1 + 2*x - x2) +  
           NC*(2 + 8*x - x2)) + 176*NC*x2 -  
        12*lx*(NCi*(-1 + 2*x - x2) +  
           NC*(5 + 2*x + 5*x2))))/24. +  
   (NC*lomx*l22)/2. + NC*x*lomx*l22 +  
   (NC*lopx*l22)/2. + NC*x*lopx*l22 +  
   NC*lomx*x2*l22 +  
   NC*lopx*x2*l22 - (NC*l23)/3. -  
   (2*NC*x*l23)/3. - (2*NC*x2*l23)/3. -  
   (3*NC*lomx2)/16. + (13*NC*x*lomx2)/4. -  
   (3*NC*lx*lomx2)/2. +  
   (7*NCi*lomx2)/16. -  
   (7*x*NCi*lomx2)/4. +  
   lx*NCi*lomx2 -  
   2*x*lx*NCi*lomx2 +  
   (NC*xi*lomx2)/3. -  
   (43*NC*x2*lomx2)/12. -  
   (7*NC*lx*x2*lomx2)/4. +  
   (3*NCi*x2*lomx2)/2. +  
   (7*lx*NCi*x2*lomx2)/4. +  
   (7*NC*lomx3)/24. - (7*NC*x*lomx3)/12. -  
   (5*NCi*lomx3)/24. +  
   (5*x*NCi*lomx3)/12. +  
   (7*NC*x2*lomx3)/12. -  
   (5*NCi*x2*lomx3)/12. -  
   (59*NC*lx2)/32. + (19*NC*x*lx2)/8. +  
   (11*NC*lomx*lx2)/8. -  
   (11*NC*x*lomx*lx2)/4. + NC*lopx*lx2 +  
   2*NC*x*lopx*lx2 + (7*NCi*lx2)/32. -  
   (15*x*NCi*lx2)/8. -  
   (3*lomx*NCi*lx2)/8. +  
   (3*x*lomx*NCi*lx2)/4. -  
   (137*NC*x2*lx2)/12. +  
   (11*NC*lomx*x2*lx2)/4. +  
   2*NC*lopx*x2*lx2 +  
   (3*NCi*x2*lx2)/2. -  
   (3*lomx*NCi*x2*lx2)/4. +  
   (19*NC*lx3)/48. + (29*NC*x*lx3)/24. +  
   (NCi*lx3)/48. - (x*NCi*lx3)/24. -  
   (NC*x2*lx3)/4. + (NCi*x2*lx3)/4. +  
   (NC*lomx*lopx2)/2. + NC*x*lomx*lopx2 +  
   NC*lomx*x2*lopx2 - (NC*lopx3)/6. -  
   (NC*x*lopx3)/3. - (NC*x2*lopx3)/3.;
};
double C2TQ2QNS_Rx_Lz(int const& nf, double const& x)
{
  const double NC   = apfel::NC;
  const double NC2  = NC * NC;
  const double NCi  = 1. / NC;
  const double NCi2 = NCi * NCi;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double x2 = x * x;
  const double xi  = 1. / x;
  const double omxi = 1. / ( 1 - x );
  const double lx  = log(x);
  const double lx2 = lx * lx;
  const double lx3 = lx * lx2;
  const double lomx  = log(1 - x);
  const double lomx2 = lomx * lomx;
  const double lomx3 = lomx * lomx2;
  const double Li2x  = apfel::dilog(x);
  const double Li3x       = apfel::wgplg(2,1,x);
  const double Li3omx     = apfel::wgplg(2,1,1 - x);
  return (1748 - 246*NC - 148*NC*nf - 3364*x + 102*NC*x - 76*NC*nf*x +  
     2916*apfel::zeta3 + 2916*x*apfel::zeta3 - 1296*Li3x - 1296*x*Li3x -  
     672*lomx - 396*NC*lomx + 48*NC*nf*lomx -  
     1176*x*lomx + 612*NC*x*lomx + 192*NC*nf*x*lomx +  
     1116*lx + 828*NC*lx - 180*NC*nf*lx + 1908*x*lx -  
     180*NC*x*lx - 252*NC*nf*x*lx + 612*lomx*lx -  
     108*NC*lomx*lx + 72*NC*nf*lomx*lx +  
     180*x*lomx*lx + 108*NC*x*lomx*lx +  
     72*NC*nf*x*lomx*lx - 1566*NCi2 + 1674*x*NCi2 -  
     1188*apfel::zeta3*NCi2 - 1188*x*apfel::zeta3*NCi2 +  
     756*Li3x*NCi2 + 756*x*Li3x*NCi2 +  
     216*lomx*NCi2 + 1458*x*lomx*NCi2 -  
     1134*lx*NCi2 - 2376*x*lx*NCi2 -  
     108*lomx*lx*NCi2 + 108*x*lomx*lx*NCi2 +  
     246*NCi + 148*nf*NCi - 102*x*NCi +  
     76*nf*x*NCi + 396*lomx*NCi -  
     48*nf*lomx*NCi - 612*x*lomx*NCi -  
     192*nf*x*lomx*NCi - 828*lx*NCi +  
     180*nf*lx*NCi + 180*x*lx*NCi +  
     252*nf*x*lx*NCi + 108*lomx*lx*NCi -  
     72*nf*lomx*lx*NCi -  
     108*x*lomx*lx*NCi -  
     72*nf*x*lomx*lx*NCi - 182*NC2 +  
     1690*x*NC2 - 1728*apfel::zeta3*NC2 - 1728*x*apfel::zeta3*NC2 +  
     540*Li3x*NC2 + 540*x*Li3x*NC2 +  
     456*lomx*NC2 - 282*x*lomx*NC2 +  
     18*lx*NC2 + 468*x*lx*NC2 -  
     504*lomx*lx*NC2 - 288*x*lomx*lx*NC2 +  
     42*pi2 + 36*NC*pi2 + 24*NC*nf*pi2 + 258*x*pi2 +  
     108*NC*x*pi2 + 24*NC*nf*x*pi2 - 216*lomx*pi2 +  
     36*NC*lomx*pi2 - 216*x*lomx*pi2 +  
     36*NC*x*lomx*pi2 + 360*lx*pi2 -  
     72*NC*lx*pi2 + 360*x*lx*pi2 -  
     72*NC*x*lx*pi2 + 36*NCi2*pi2 -  
     90*x*NCi2*pi2 + 108*lomx*NCi2*pi2 +  
     108*x*lomx*NCi2*pi2 -  
     144*lx*NCi2*pi2 - 144*x*lx*NCi2*pi2 -  
     36*NCi*pi2 - 24*nf*NCi*pi2 -  
     108*x*NCi*pi2 - 24*nf*x*NCi*pi2 -  
     36*lomx*NCi*pi2 -  
     36*x*lomx*NCi*pi2 + 72*lx*NCi*pi2 +  
     72*x*lx*NCi*pi2 - 78*NC2*pi2 -  
     168*x*NC2*pi2 + 108*lomx*NC2*pi2 +  
     108*x*lomx*NC2*pi2 - 216*lx*NC2*pi2 -  
     216*x*lx*NC2*pi2 - 2592*apfel::zeta3*omxi +  
     2592*Li3x*omxi - 1800*lx*omxi +  
     360*NC*nf*lx*omxi - 792*lomx*lx*omxi -  
     144*NC*nf*lomx*lx*omxi +  
     1512*apfel::zeta3*NCi2*omxi -  
     1512*Li3x*NCi2*omxi +  
     2160*lx*NCi2*omxi -  
     360*nf*lx*NCi*omxi +  
     144*nf*lomx*lx*NCi*omxi +  
     1080*apfel::zeta3*NC2*omxi -  
     1080*Li3x*NC2*omxi -  
     360*lx*NC2*omxi +  
     792*lomx*lx*NC2*omxi -  
     24*pi2*omxi - 24*NC*nf*pi2*omxi -  
     576*lx*pi2*omxi -  
     54*NCi2*pi2*omxi -  
     36*lomx*NCi2*pi2*omxi +  
     216*lx*NCi2*pi2*omxi +  
     24*nf*NCi*pi2*omxi +  
     78*NC2*pi2*omxi +  
     36*lomx*NC2*pi2*omxi +  
     360*lx*NC2*pi2*omxi -  
     216*Li3omx*(NC + NC*x - (1 + x)*NCi - 2*NC2 -  
        2*x*NC2 + NCi2*(2 + 2*x - 3*omxi) -  
        2*omxi + 5*NC2*omxi) + 416*NC*xi +  
     312*NC*lomx*xi - 144*NC*lomx*lx*xi -  
     416*NCi*xi - 312*lomx*NCi*xi +  
     144*lomx*lx*NCi*xi - 272*NC*x2 -  
     528*NC*lomx*x2 + 912*NC*lx*x2 +  
     144*NC*lomx*lx*x2 + 272*NCi*x2 +  
     528*lomx*NCi*x2 - 912*lx*NCi*x2 -  
     144*lomx*lx*NCi*x2 + 72*NC*pi2*x2 -  
     72*NCi*pi2*x2 -  
     36*Li2x*(5 + 9*NC + 2*NC*nf + 23*x + 15*NC*x + 2*NC*nf*x -  
        9*x*NCi2 - 9*NCi - 2*nf*NCi - 15*x*NCi -  
        2*nf*x*NCi - 5*NC2 - 14*x*NC2 - 4*omxi -  
        4*NC*nf*omxi - 9*NCi2*omxi +  
        4*nf*NCi*omxi + 13*NC2*omxi +  
        6*lomx*(-1 + NC - x + NC*x - (1 + x)*NCi +  
           NCi2*(1 + x - omxi) + NC2*omxi) -  
        6*lx*(2 + NC + 2*x + NC*x - (1 + x)*NCi -  
           6*omxi + NC2*omxi +  
           NCi2*(-2*(1 + x) + 5*omxi)) + 4*NC*xi -  
        4*NCi*xi + 8*NC*x2 - 8*NCi*x2) -  
     198*lomx2 + 54*NC*lomx2 -  
     36*NC*nf*lomx2 - 198*x*lomx2 -  
     54*NC*x*lomx2 - 36*NC*nf*x*lomx2 -  
     648*lx*lomx2 - 108*NC*lx*lomx2 -  
     648*x*lx*lomx2 - 108*NC*x*lx*lomx2 +  
     216*lx*NCi2*lomx2 +  
     216*x*lx*NCi2*lomx2 -  
     54*NCi*lomx2 + 36*nf*NCi*lomx2 +  
     54*x*NCi*lomx2 +  
     36*nf*x*NCi*lomx2 +  
     108*lx*NCi*lomx2 +  
     108*x*lx*NCi*lomx2 +  
     198*NC2*lomx2 + 198*x*NC2*lomx2 +  
     432*lx*NC2*lomx2 +  
     432*x*lx*NC2*lomx2 +  
     1512*lx*omxi*lomx2 -  
     540*lx*NCi2*omxi*lomx2 -  
     972*lx*NC2*omxi*lomx2 +  
     72*NC*xi*lomx2 -  
     72*NCi*xi*lomx2 -  
     72*NC*x2*lomx2 +  
     72*NCi*x2*lomx2 + 216*lomx3 +  
     216*x*lomx3 - 108*NCi2*lomx3 -  
     108*x*NCi2*lomx3 - 108*NC2*lomx3 -  
     108*x*NC2*lomx3 - 387*lx2 -  
     351*NC*lx2 - 90*NC*nf*lx2 - 819*x*lx2 -  
     351*NC*x*lx2 - 90*NC*nf*x*lx2 +  
     756*lomx*lx2 + 756*x*lomx*lx2 -  
     54*NCi2*lx2 + 324*x*NCi2*lx2 -  
     486*lomx*NCi2*lx2 -  
     486*x*lomx*NCi2*lx2 +  
     351*NCi*lx2 + 90*nf*NCi*lx2 +  
     351*x*NCi*lx2 + 90*nf*x*NCi*lx2 +  
     441*NC2*lx2 + 495*x*NC2*lx2 -  
     270*lomx*NC2*lx2 -  
     270*x*lomx*NC2*lx2 +  
     504*omxi*lx2 +  
     180*NC*nf*omxi*lx2 -  
     1512*lomx*omxi*lx2 +  
     243*NCi2*omxi*lx2 +  
     972*lomx*NCi2*omxi*lx2 -  
     180*nf*NCi*omxi*lx2 -  
     747*NC2*omxi*lx2 +  
     540*lomx*NC2*omxi*lx2 -  
     360*NC*x2*lx2 + 360*NCi*x2*lx2 -  
     144*lx3 + 90*NC*lx3 - 144*x*lx3 +  
     90*NC*x*lx3 + 45*NCi2*lx3 +  
     45*x*NCi2*lx3 - 90*NCi*lx3 -  
     90*x*NCi*lx3 + 99*NC2*lx3 +  
     99*x*NC2*lx3 + 108*omxi*lx3 -  
     108*NC2*omxi*lx3)/432.;
};
double C2TQ2QPS_Rx_Lz(int const&, double const&)
{
  return 0;
};
double C2TQ2QB_Rx_Lz(int const&, double const& x)
{
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double NCi2 = NCi * NCi;
  const double l2  = log(2);
  const double l22 = l2 * l2;
  const double l23 = l2 * l22;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double x2 = x * x;
  const double xi  = 1. / x;
  const double opxi = 1. / ( 1 + x );
  const double lx  = log(x);
  const double lx2 = lx * lx;
  const double lx3 = lx * lx2;
  const double lomx  = log(1 - x);
  const double lomx2 = lomx * lomx;
  const double lopx  = log(1 + x);
  const double lopx2 = lopx * lopx;
  const double lopx3 = lopx * lopx2;
  const double Li2x  = apfel::dilog(x);
  const double Li2mx = apfel::dilog(-x);
  const double Li3x       = apfel::wgplg(2,1,x);
  const double Li3omx     = apfel::wgplg(2,1,1 - x);
  const double Li3omxh    = apfel::wgplg(2,1,0.5 - x/2.);
  const double Li3opxh    = apfel::wgplg(2,1,(1 + x)/2.);
  const double Li3omxopxi = apfel::wgplg(2,1,-((-1 + x)*opxi));
  return -1.875 - (41*NC)/72. + (15*x)/8. + (17*NC*x)/72. - apfel::zeta3/2. +  
   (x*apfel::zeta3)/2. - Li3omx - (NC*Li3omx)/2. + x*Li3omx -  
   (NC*x*Li3omx)/2. + Li3omxh - x*Li3omxh - Li3x/2. +  
   (x*Li3x)/2. + Li3opxh - x*Li3opxh +  
   Li3omxopxi - x*Li3omxopxi -  
   lomx - (11*NC*lomx)/12. + x*lomx +  
   (17*NC*x*lomx)/12. - (3*lx)/8. + (23*NC*lx)/12. -  
   (19*x*lx)/8. - (5*NC*x*lx)/12. - (NC*lomx*lx)/4. +  
   (NC*x*lomx*lx)/4. - l2*lomx*lopx +  
   x*l2*lomx*lopx + (lx*lopx)/2. +  
   (x*lx*lopx)/2. - lomx*lx*lopx +  
   x*lomx*lx*lopx + apfel::zeta3*opxi +  
   2*Li3omx*opxi - 2*Li3omxh*opxi +  
   Li3x*opxi - 2*Li3opxh*opxi -  
   2*Li3omxopxi*opxi +  
   2*l2*lomx*lopx*opxi +  
   2*lomx*lx*lopx*opxi -  
   (Li2mx*(1 + x + 2*lomx*(-1 + x + 2*opxi) -  
        2*lx*(-1 + x + 2*opxi))*(-1 + NCi2))/2. +  
   (15*NCi2)/8. - (15*x*NCi2)/8. + (apfel::zeta3*NCi2)/2. -  
   (x*apfel::zeta3*NCi2)/2. + Li3omx*NCi2 -  
   x*Li3omx*NCi2 - Li3omxh*NCi2 +  
   x*Li3omxh*NCi2 + (Li3x*NCi2)/2. -  
   (x*Li3x*NCi2)/2. - Li3opxh*NCi2 +  
   x*Li3opxh*NCi2 - Li3omxopxi*NCi2 +  
   x*Li3omxopxi*NCi2 + lomx*NCi2 -  
   x*lomx*NCi2 + (3*lx*NCi2)/8. +  
   (19*x*lx*NCi2)/8. + l2*lomx*lopx*NCi2 -  
   x*l2*lomx*lopx*NCi2 -  
   (lx*lopx*NCi2)/2. - (x*lx*lopx*NCi2)/2. +  
   lomx*lx*lopx*NCi2 -  
   x*lomx*lx*lopx*NCi2 -  
   apfel::zeta3*opxi*NCi2 - 2*Li3omx*opxi*NCi2 +  
   2*Li3omxh*opxi*NCi2 -  
   Li3x*opxi*NCi2 +  
   2*Li3opxh*opxi*NCi2 +  
   2*Li3omxopxi*opxi*NCi2 -  
   2*l2*lomx*lopx*opxi*NCi2 -  
   2*lomx*lx*lopx*opxi*NCi2 +  
   (41*NCi)/72. - (17*x*NCi)/72. + (Li3omx*NCi)/2. +  
   (x*Li3omx*NCi)/2. + (11*lomx*NCi)/12. -  
   (17*x*lomx*NCi)/12. - (23*lx*NCi)/12. +  
   (5*x*lx*NCi)/12. + (lomx*lx*NCi)/4. -  
   (x*lomx*lx*NCi)/4. - pi2/24. + (NC*pi2)/12. -  
   (x*pi2)/24. + (NC*x*pi2)/4. + (l2*pi2)/6. -  
   (x*l2*pi2)/6. - (lopx*pi2)/4. +  
   (x*lopx*pi2)/4. - (lomx*pi2)/12. +  
   (NC*lomx*pi2)/12. + (x*lomx*pi2)/12. +  
   (NC*x*lomx*pi2)/12. + (lx*pi2)/6. -  
   (NC*lx*pi2)/6. - (x*lx*pi2)/6. -  
   (NC*x*lx*pi2)/6. + (lopx*pi2)/12. -  
   (x*lopx*pi2)/12. - (l2*opxi*pi2)/3. +  
   (lopx*opxi*pi2)/2. +  
   (lomx*opxi*pi2)/6. -  
   (lx*opxi*pi2)/3. -  
   (lopx*opxi*pi2)/6. + (NCi2*pi2)/24. +  
   (x*NCi2*pi2)/24. - (l2*NCi2*pi2)/6. +  
   (x*l2*NCi2*pi2)/6. +  
   (lopx*NCi2*pi2)/4. -  
   (x*lopx*NCi2*pi2)/4. +  
   (lomx*NCi2*pi2)/12. -  
   (x*lomx*NCi2*pi2)/12. -  
   (lx*NCi2*pi2)/6. + (x*lx*NCi2*pi2)/6. -  
   (lopx*NCi2*pi2)/12. +  
   (x*lopx*NCi2*pi2)/12. +  
   (l2*opxi*NCi2*pi2)/3. -  
   (lopx*opxi*NCi2*pi2)/2. -  
   (lomx*opxi*NCi2*pi2)/6. +  
   (lx*opxi*NCi2*pi2)/3. +  
   (lopx*opxi*NCi2*pi2)/6. -  
   (NCi*pi2)/12. - (x*NCi*pi2)/4. -  
   (lomx*NCi*pi2)/12. -  
   (x*lomx*NCi*pi2)/12. +  
   (lx*NCi*pi2)/6. + (x*lx*NCi*pi2)/6. +  
   (26*NC*xi)/27. + (13*NC*lomx*xi)/18. -  
   (NC*lomx*lx*xi)/3. - (26*NCi*xi)/27. -  
   (13*lomx*NCi*xi)/18. +  
   (lomx*lx*NCi*xi)/3. - (17*NC*x2)/27. -  
   (11*NC*lomx*x2)/9. + (19*NC*lx*x2)/9. +  
   (NC*lomx*lx*x2)/3. + (17*NCi*x2)/27. +  
   (11*lomx*NCi*x2)/9. -  
   (19*lx*NCi*x2)/9. -  
   (lomx*lx*NCi*x2)/3. + (NC*pi2*x2)/6. -  
   (NCi*pi2*x2)/6. +  
   (Li2x*(6 - 9*NC + 6*x - 15*NC*x - 6*NCi2 - 6*x*NCi2 -  
        6*(1 + x)*lomx*(NC - NCi) +  
        6*(1 + x)*lx*(NC - NCi) + 9*NCi + 15*x*NCi -  
        4*NC*xi + 4*NCi*xi - 8*NC*x2 +  
        8*NCi*x2))/12. + (lomx*l22)/2. -  
   (x*lomx*l22)/2. + (lopx*l22)/2. -  
   (x*lopx*l22)/2. - lomx*opxi*l22 -  
   lopx*opxi*l22 -  
   (lomx*NCi2*l22)/2. +  
   (x*lomx*NCi2*l22)/2. -  
   (lopx*NCi2*l22)/2. +  
   (x*lopx*NCi2*l22)/2. +  
   lomx*opxi*NCi2*l22 +  
   lopx*opxi*NCi2*l22 - l23/3. +  
   (x*l23)/3. + (2*opxi*l23)/3. +  
   (NCi2*l23)/3. - (x*NCi2*l23)/3. -  
   (2*opxi*NCi2*l23)/3. +  
   (NC*lomx2)/8. - (NC*x*lomx2)/8. -  
   (NC*lx*lomx2)/4. - (NC*x*lx*lomx2)/4. -  
   (NCi*lomx2)/8. + (x*NCi*lomx2)/8. +  
   (lx*NCi*lomx2)/4. +  
   (x*lx*NCi*lomx2)/4. +  
   (NC*xi*lomx2)/6. -  
   (NCi*xi*lomx2)/6. -  
   (NC*x2*lomx2)/6. +  
   (NCi*x2*lomx2)/6. - (13*NC*lx2)/16. +  
   (x*lx2)/2. - (13*NC*x*lx2)/16. +  
   lopx*lx2 - x*lopx*lx2 -  
   2*lopx*opxi*lx2 -  
   (x*NCi2*lx2)/2. - lopx*NCi2*lx2 +  
   x*lopx*NCi2*lx2 +  
   2*lopx*opxi*NCi2*lx2 +  
   (13*NCi*lx2)/16. + (13*x*NCi*lx2)/16. -  
   (5*NC*x2*lx2)/6. +  
   (5*NCi*x2*lx2)/6. - lx3/8. +  
   (5*NC*lx3)/24. + (x*lx3)/8. +  
   (5*NC*x*lx3)/24. + (opxi*lx3)/4. +  
   (NCi2*lx3)/8. - (x*NCi2*lx3)/8. -  
   (opxi*NCi2*lx3)/4. -  
   (5*NCi*lx3)/24. - (5*x*NCi*lx3)/24. +  
   (lomx*lopx2)/2. - (x*lomx*lopx2)/2. -  
   lomx*opxi*lopx2 -  
   (lomx*NCi2*lopx2)/2. +  
   (x*lomx*NCi2*lopx2)/2. +  
   lomx*opxi*NCi2*lopx2 -  
   lopx3/6. + (x*lopx3)/6. +  
   (opxi*lopx3)/3. + (NCi2*lopx3)/6. -  
   (x*NCi2*lopx3)/6. -  
   (opxi*NCi2*lopx3)/3.;
};
double C2TQ2QP1_Rx_Lz(int const&, double const&)
{
  return 0;
};
double C2TQ2QP2_Rx_Lz(int const&, double const& x)
{
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double x2 = x * x;
  const double xi  = 1. / x;
  const double lx  = log(x);
  const double lx2 = lx * lx;
  const double lx3 = lx * lx2;
  const double lomx  = log(1 - x);
  const double lomx2 = lomx * lomx;
  const double Li2x  = apfel::dilog(x);
  const double Li3omx     = apfel::wgplg(2,1,1 - x);
  return -0.0023148148148148147* 
   ((NC - NCi)*(246 - 102*x + 216*(1 + x)*Li3omx +  
       396*lomx - 612*x*lomx - 828*lx + 180*x*lx +  
       108*lomx*lx - 108*x*lomx*lx - 36*pi2 -  
       108*x*pi2 - 36*lomx*pi2 -  
       36*x*lomx*pi2 + 72*lx*pi2 +  
       72*x*lx*pi2 - 416*xi - 312*lomx*xi +  
       144*lomx*lx*xi + 272*x2 +  
       528*lomx*x2 - 912*lx*x2 -  
       144*lomx*lx*x2 - 72*pi2*x2 +  
       36*Li2x*(9 + 15*x + 6*(1 + x)*lomx - 6*(1 + x)*lx +  
          4*xi + 8*x2) - 54*lomx2 +  
       54*x*lomx2 + 108*lx*lomx2 +  
       108*x*lx*lomx2 - 72*xi*lomx2 +  
       72*x2*lomx2 + 351*lx2 +  
       351*x*lx2 + 360*x2*lx2 - 90*lx3 -  
       90*x*lx3));
};
double C2TQ2QP3_Rx_Lz(int const&, double const&)
{
  return 0;
};
//------------------------------------------------------
double C1LQ2Q_Sx_Rz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C1LQ2G_Sx_Rz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C1LG2Q_Sx_Rz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C1TQ2Q_Sx_Rz(int const&, double const& x, double const& z, bool const& intx)
{
  const double NC  = apfel::NC;
  const double NCi = 1. / NC;
  return -0.5*((1 + z)*Dn(0,1 - x, intx)*(NC - NCi));
};
double C1TQ2G_Sx_Rz(int const&, double const& x, double const& z, bool const& intx)
{
  const double NC  = apfel::NC;
  const double NCi = 1. / NC;
  const double zi  = 1. / z;
  return (Dn(0,1 - x, intx)*(NC - NCi)*(-2 + z + 2*zi))/2.;
};
double C1TG2Q_Sx_Rz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2G_Sx_Rz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2LG2G_Sx_Rz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2LG2Q_Sx_Rz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QNS_Sx_Rz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QPS_Sx_Rz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QB_Sx_Rz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QP1_Sx_Rz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QP2_Sx_Rz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QP3_Sx_Rz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2TQ2G_Sx_Rz(int const&, double const& x, double const& z, bool const& intx)
{
  const double NC   = apfel::NC;
  const double NC2  = NC * NC;
  const double NCi  = 1. / NC;
  const double NCi2 = NCi * NCi;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double z2 = z * z;
  const double zi  = 1. / z;
  const double lz  = log(z);
  const double lz2 = lz * lz;
  const double lomz  = log(1 - z);
  const double lomz2 = lomz * lomz;
  const double lopz  = log(1 + z);
  const double Li2z  = apfel::dilog(z);
  const double Li2mz = apfel::dilog(-z);
  return (54*Dn(2,1 - x, intx)*(-2 + NCi2 + NC2)* 
      (-2 + z + 2*zi) + 3*Dn(1,1 - x, intx)* 
      (-192 + 6*z + 48*NCi2 - 9*z*NCi2 + 144*NC2 +  
        3*z*NC2 + 196*zi - 36*NCi2*zi -  
        160*NC2*zi +  
        24*lomz*(-3 + NCi2 + 2*NC2)* 
         (-2 + z + 2*zi) +  
        6*lz*(6*(2 + z) + NCi2*(-2 + z + 4*zi) -  
           NC2*(10 + 7*z + 4*zi)) - 16*z2 +  
        16*NC2*z2) +  
     Dn(0,1 - x, intx)*(-438 + 276*z - 576*lomz + 18*z*lomz -  
        72*lz - 72*z*lz + 144*lomz*lz -  
        72*z*lomz*lz - 144*lz*lopz -  
        72*z*lz*lopz + 153*NCi2 - 90*z*NCi2 +  
        144*lomz*NCi2 - 27*z*lomz*NCi2 +  
        72*lomz*lz*NCi2 - 36*z*lomz*lz*NCi2 +  
        285*NC2 - 186*z*NC2 + 432*lomz*NC2 +  
        9*z*lomz*NC2 + 72*lz*NC2 +  
        72*z*lz*NC2 - 216*lomz*lz*NC2 +  
        108*z*lomz*lz*NC2 + 144*lz*lopz*NC2 +  
        72*z*lz*lopz*NC2 - 36*pi2 + 42*z*pi2 +  
        30*NCi2*pi2 - 15*z*NCi2*pi2 +  
        6*NC2*pi2 - 27*z*NC2*pi2 + 272*zi +  
        588*lomz*zi + 588*lz*zi -  
        144*lomz*lz*zi - 144*lz*lopz*zi -  
        162*NCi2*zi - 108*lomz*NCi2*zi -  
        54*lz*NCi2*zi -  
        72*lomz*lz*NCi2*zi - 110*NC2*zi -  
        480*lomz*NC2*zi - 534*lz*NC2*zi +  
        216*lomz*lz*NC2*zi +  
        144*lz*lopz*NC2*zi + 48*pi2*zi -  
        24*NCi2*pi2*zi -  
        24*NC2*pi2*zi +  
        72*Li2mz*(-1 + NC2)*(2 + z + 2*zi) -  
        18*Li2z*(6*(2 + z) + NCi2*(5*(-2 + z) + 12*zi) -  
           NC2*(2 + 11*z + 12*zi)) + 52*z2 -  
        48*lomz*z2 - 48*lz*z2 - 52*NC2*z2 +  
        48*lomz*NC2*z2 + 48*lz*NC2*z2 +  
        216*lomz2 - 108*z*lomz2 -  
        72*NCi2*lomz2 + 36*z*NCi2*lomz2 -  
        144*NC2*lomz2 + 72*z*NC2*lomz2 -  
        216*zi*lomz2 +  
        72*NCi2*zi*lomz2 +  
        144*NC2*zi*lomz2 + 72*lz2 +  
        216*z*lz2 + 36*NCi2*lz2 -  
        18*z*NCi2*lz2 - 108*NC2*lz2 -  
        198*z*NC2*lz2 + 216*zi*lz2 -  
        216*NC2*zi*lz2))/72.;
};
double C2TG2G_Sx_Rz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2TG2Q_Sx_Rz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2TQ2QNS_Sx_Rz(int const& nf, double const& x, double const& z, bool const& intx)
{
  const double NC   = apfel::NC;
  const double NC2  = NC * NC;
  const double NCi  = 1. / NC;
  const double NCi2 = NCi * NCi;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double z2 = z * z;
  const double zi  = 1. / z;
  const double omzi  = 1. / ( 1 - z );
  const double lz  = log(z);
  const double lz2 = lz * lz;
  const double lomz  = log(1 - z);
  const double lomz2 = lomz * lomz;
  const double Li2z  = apfel::dilog(z);
  return (-6*(9*(1 + z)*Dn(2,1 - x, intx)*(-2 + NCi2 + NC2) +  
        Dn(1,1 - x, intx)*(11 - 3*NC + 2*NC*nf + 11*z + 3*NC*z + 2*NC*nf*z +  
           3*NCi - 2*nf*NCi - 3*z*NCi -  
           2*nf*z*NCi - 11*NC2 - 11*z*NC2 +  
           18*(1 + z)*lomz*(-2 + NCi2 + NC2) -  
           3*(1 + z)*lz*(-2 + 2*NC + NCi2 - 2*NCi +  
              NC2) - 4*NC*zi + 4*NCi*zi +  
           4*NC*z2 - 4*NCi*z2)) -  
     Dn(0,1 - x, intx)*(112 + 120*NC - 8*NC*nf + 196*z - 84*NC*z - 32*NC*nf*z -  
        216*lz + 18*NC*lz - 36*z*lz + 72*NC*z*lz -  
        36*NCi2 - 252*z*NCi2 + 126*lz*NCi2 +  
        36*z*lz*NCi2 - 120*NCi + 8*nf*NCi +  
        84*z*NCi + 32*nf*z*NCi - 18*lz*NCi -  
        72*z*lz*NCi - 76*NC2 + 56*z*NC2 +  
        90*lz*NC2 + 18*(1 + z)*Li2z* 
         (-2 + 2*NC + NCi2 - 2*NCi + NC2) + 36*pi2 -  
        6*NC*pi2 + 36*z*pi2 - 6*NC*z*pi2 -  
        15*NCi2*pi2 - 15*z*NCi2*pi2 +  
        6*NCi*pi2 + 6*z*NCi*pi2 -  
        21*NC2*pi2 - 21*z*NC2*pi2 +  
        108*lz*omzi - 54*lz*NCi2*omzi -  
        54*lz*NC2*omzi + 14*NC*zi -  
        24*NC*lz*zi - 14*NCi*zi +  
        24*lz*NCi*zi - 50*NC*z2 +  
        24*NC*lz*z2 + 50*NCi*z2 -  
        24*lz*NCi*z2 +  
        6*lomz*(11 - 3*NC + 2*NC*nf + 11*z + 3*NC*z + 2*NC*nf*z -  
           11*NC2 - 11*z*NC2 +  
           6*lz*(-2 + NCi2 + NC2)*(1 + z - 2*omzi) -  
           4*NC*zi + 4*NC*z2 -  
           NCi*(-3 + 2*nf + 3*z + 2*nf*z - 4*zi + 4*z2)) 
         - 108*lomz2 - 108*z*lomz2 +  
        54*NCi2*lomz2 + 54*z*NCi2*lomz2 +  
        54*NC2*lomz2 + 54*z*NC2*lomz2 +  
        126*lz2 - 36*NC*lz2 + 126*z*lz2 -  
        36*NC*z*lz2 - 72*NCi2*lz2 -  
        72*z*NCi2*lz2 + 36*NCi*lz2 +  
        36*z*NCi*lz2 - 54*NC2*lz2 -  
        54*z*NC2*lz2 - 180*omzi*lz2 +  
        108*NCi2*omzi*lz2 +  
        72*NC2*omzi*lz2))/72.;
};
double C2TQ2QPS_Sx_Rz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2TQ2QB_Sx_Rz(int const& nf, double const& x, double const& z, bool const& intx)
{
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double NCi2 = NCi * NCi;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double z2 = z * z;
  const double zi  = 1. / z;
  const double opzi  = 1. / ( 1 + z );
  const double lz  = log(z);
  const double lz2 = lz * lz;
  const double lomz  = log(1 - z);
  const double lopz  = log(1 + z);
  const double Li2z  = apfel::dilog(z);
  const double Li2mz = apfel::dilog(-z);
  return (3*Dn(1,1 - x, intx)*(NC - NCi)* 
      (3 - 3*z + 6*(1 + z)*lz + 4*zi - 4*z2) -  
     Dn(0,1 - x, intx)*(36 + 60*NC - 36*z - 42*NC*z - 9*NC*lomz +  
        9*NC*z*lomz + 18*lz + 9*NC*lz + 18*z*lz +  
        36*NC*z*lz + 36*lz*lopz - 36*z*lz*lopz -  
        72*lz*lopz*opzi +  
        36*Li2mz*(-1 + z + 2*opzi)*(-1 + NCi2) -  
        36*NCi2 + 36*z*NCi2 - 18*lz*NCi2 -  
        18*z*lz*NCi2 - 36*lz*lopz*NCi2 +  
        36*z*lz*lopz*NCi2 +  
        72*lz*lopz*opzi*NCi2 +  
        18*(1 + z)*Li2z*(NC - NCi) - 60*NCi +  
        42*z*NCi + 9*lomz*NCi -  
        9*z*lomz*NCi - 9*lz*NCi -  
        36*z*lz*NCi + 3*pi2 - 3*NC*pi2 -  
        3*z*pi2 - 3*NC*z*pi2 - 6*opzi*pi2 -  
        3*NCi2*pi2 + 3*z*NCi2*pi2 +  
        6*opzi*NCi2*pi2 + 3*NCi*pi2 +  
        3*z*NCi*pi2 + 7*NC*zi -  
        12*NC*lomz*zi - 12*NC*lz*zi -  
        7*NCi*zi + 12*lomz*NCi*zi +  
        12*lz*NCi*zi - 25*NC*z2 +  
        12*NC*lomz*z2 + 12*NC*lz*z2 +  
        25*NCi*z2 - 12*lomz*NCi*z2 -  
        12*lz*NCi*z2 - 9*lz2 -  
        18*NC*lz2 + 9*z*lz2 - 18*NC*z*lz2 +  
        18*opzi*lz2 + 9*NCi2*lz2 -  
        9*z*NCi2*lz2 -  
        18*opzi*NCi2*lz2 +  
        18*NCi*lz2 + 18*z*NCi*lz2))/36.;
};
double C2TQ2QP1_Sx_Rz(int const&, double const& x, double const& z, bool const& intx)
{
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double z2 = z * z;
  const double zi  = 1. / z;
  const double lz  = log(z);
  const double lz2 = lz * lz;
  const double lomz  = log(1 - z);
  const double Li2z  = apfel::dilog(z);
  return -0.027777777777777776*((NC - NCi)* 
     (-3*Dn(1,1 - x, intx)*(3 - 3*z + 6*(1 + z)*lz + 4*zi -  
          4*z2) + Dn(0,1 - x, intx)* 
        (60 - 42*z + 18*(1 + z)*Li2z + 9*lz + 36*z*lz -  
          3*pi2 - 3*z*pi2 + 7*zi - 12*lz*zi -  
          25*z2 + 12*lz*z2 +  
          3*lomz*(-3 + 3*z - 4*zi + 4*z2) -  
          18*lz2 - 18*z*lz2)));
};
double C2TQ2QP2_Sx_Rz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2TQ2QP3_Sx_Rz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
//------------------------------------------------------
double C1LQ2Q_Rx_Sz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C1LQ2G_Rx_Sz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C1LG2Q_Rx_Sz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C1TQ2Q_Rx_Sz(int const&, double const& x, double const& z, bool const& intz)
{
  const double NC = apfel::NC;
  const double NCi = 1. / NC;
  return -0.5*((1 + x)*Dn(0,1 - z, intz)*(NC - NCi));
};
double C1TQ2G_Rx_Sz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C1TG2Q_Rx_Sz(int const&, double const& x, double const& z, bool const& intz)
{
  const double x2 = x * x;
  return Dn(0,1 - z, intz)*(0.5 - x + x2);
};
double C2LQ2G_Rx_Sz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2LG2G_Rx_Sz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2LG2Q_Rx_Sz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QNS_Rx_Sz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QPS_Rx_Sz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QB_Rx_Sz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QP1_Rx_Sz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QP2_Rx_Sz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2LQ2QP3_Rx_Sz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2TQ2G_Rx_Sz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2TG2G_Rx_Sz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2TG2Q_Rx_Sz(int const& nf, double const& x, double const& z, bool const& intz)
{
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double x2 = x * x;
  const double xi  = 1. / x;
  const double lx  = log(x);
  const double lx2 = lx * lx;
  const double lomx  = log(1 - x);
  const double lomx2 = lomx * lomx;
  const double lopx  = log(1 + x);
  const double Li2x  = apfel::dilog(x);
  const double Li2mx = apfel::dilog(-x);
  return (-54*Dn(2,1 - z, intz)*(NC - NCi)*(-1 + 2*x - 2*x2) -  
     3*Dn(1,1 - z, intz)*(9*NC - 168*NC*x - 21*NCi + 72*x*NCi -  
        16*NC*xi + 24*lomx*(2*NC - NCi)* 
         (-1 + 2*x - 2*x2) + 184*NC*x2 - 60*NCi*x2 -  
        6*lx*(NC*(1 + 22*x - 8*x2) +  
           NCi*(3 - 6*x + 8*x2))) -  
     Dn(0,1 - z, intz)*(177*NC - 357*NC*x + 27*NC*lomx -  
        504*NC*x*lomx - 126*NC*lx + 216*NC*x*lx +  
        144*NC*lomx*lx - 288*NC*x*lomx*lx +  
        72*NC*lx*lopx + 144*NC*x*lx*lopx -  
        27*NCi - 99*x*NCi - 63*lomx*NCi +  
        216*x*lomx*NCi + 54*lx*NCi -  
        216*x*lx*NCi - 72*lomx*lx*NCi +  
        144*x*lomx*lx*NCi +  
        18*Li2x*(NC*(5 + 14*x) + (-1 + 2*x)*NCi) + 9*NC*pi2 -  
        66*NC*x*pi2 - 9*NCi*pi2 +  
        18*x*NCi*pi2 - 104*NC*xi -  
        48*NC*lomx*xi + 347*NC*x2 +  
        552*NC*lomx*x2 - 1080*NC*lx*x2 +  
        288*NC*lomx*lx*x2 +  
        144*NC*lx*lopx*x2 + 27*NCi*x2 -  
        180*lomx*NCi*x2 + 180*lx*NCi*x2 -  
        144*lomx*lx*NCi*x2 +  
        48*NC*pi2*x2 - 24*NCi*pi2*x2 +  
        72*NC*Li2mx*(1 + 2*x + 2*x2) - 72*NC*lomx2 +  
        144*NC*x*lomx2 + 36*NCi*lomx2 -  
        72*x*NCi*lomx2 -  
        144*NC*x2*lomx2 +  
        72*NCi*x2*lomx2 + 54*NC*lx2 +  
        252*NC*x*lx2 + 18*NCi*lx2 -  
        36*x*NCi*lx2 - 72*NC*x2*lx2 +  
        72*NCi*x2*lx2))/72.;
};
double C2TQ2QNS_Rx_Sz(int const& nf, double const& x, double const& z, bool const& intz)
{
  const double NC   = apfel::NC;
  const double NC2  = NC * NC;
  const double NCi  = 1. / NC;
  const double NCi2 = NCi * NCi;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double x2 = x * x;
  const double xi  = 1. / x;
  const double omxi = 1. / ( 1 - x );
  const double lx  = log(x);
  const double lx2 = lx * lx;
  const double lomx  = log(1 - x);
  const double lomx2 = lomx * lomx;
  const double Li2x  = apfel::dilog(x);
  return (-6*(9*(1 + x)*Dn(2,1 - z, intz)*(-2 + NCi2 + NC2) +  
        Dn(1,1 - z, intz)*(11 - 3*NC + 2*NC*nf + 11*x + 3*NC*x + 2*NC*nf*x +  
           3*NCi - 2*nf*NCi - 3*x*NCi -  
           2*nf*x*NCi - 11*NC2 - 11*x*NC2 +  
           18*(1 + x)*lomx*(-2 + NCi2 + NC2) -  
           3*lx*(-10 + 2*NC - 10*x + 2*NC*x - 2*(1 + x)*NCi +  
              5*NC2 + 5*x*NC2 +  
              NCi2*(5 + 5*x - 8*omxi) + 16*omxi -  
              8*NC2*omxi) - 4*NC*xi +  
           4*NCi*xi + 4*NC*x2 - 4*NCi*x2)) -  
     Dn(0,1 - z, intz)*(112 + 66*NC - 8*NC*nf + 196*x - 102*NC*x - 32*NC*nf*x -  
        96*lx - 36*NC*lx - 24*NC*nf*lx - 204*x*lx -  
        108*NC*x*lx - 24*NC*nf*x*lx - 36*NCi2 -  
        252*x*NCi2 + 54*x*lx*NCi2 - 66*NCi +  
        8*nf*NCi + 102*x*NCi + 32*nf*x*NCi +  
        36*lx*NCi + 24*nf*lx*NCi +  
        108*x*lx*NCi + 24*nf*x*lx*NCi - 76*NC2 +  
        56*x*NC2 + 96*lx*NC2 + 150*x*lx*NC2 +  
        18*(1 + x)*Li2x*(-2 + 2*NC + NCi2 - 2*NCi +  
           NC2) + 36*pi2 - 6*NC*pi2 + 36*x*pi2 -  
        6*NC*x*pi2 - 15*NCi2*pi2 -  
        15*x*NCi2*pi2 + 6*NCi*pi2 +  
        6*x*NCi*pi2 - 21*NC2*pi2 -  
        21*x*NC2*pi2 + 156*lx*omxi +  
        48*NC*nf*lx*omxi + 54*lx*NCi2*omxi -  
        48*nf*lx*NCi*omxi -  
        210*lx*NC2*omxi - 52*NC*xi +  
        52*NCi*xi + 88*NC*x2 - 72*NC*lx*x2 -  
        88*NCi*x2 + 72*lx*NCi*x2 +  
        6*lomx*(11 - 3*NC + 2*NC*nf + 11*x + 3*NC*x + 2*NC*nf*x -  
           11*NC2 - 11*x*NC2 -  
           18*lx*(-2 + NCi2 + NC2)* 
            (1 + x - 2*omxi) - 4*NC*xi + 4*NC*x2 -  
           NCi*(-3 + 2*nf + 3*x + 2*nf*x - 4*xi + 4*x2)) 
         - 108*lomx2 - 108*x*lomx2 +  
        54*NCi2*lomx2 + 54*x*NCi2*lomx2 +  
        54*NC2*lomx2 + 54*x*NC2*lomx2 -  
        90*lx2 + 36*NC*lx2 - 90*x*lx2 +  
        36*NC*x*lx2 + 36*NCi2*lx2 +  
        36*x*NCi2*lx2 - 36*NCi*lx2 -  
        36*x*NCi*lx2 + 54*NC2*lx2 +  
        54*x*NC2*lx2 + 108*omxi*lx2 -  
        36*NCi2*omxi*lx2 -  
        72*NC2*omxi*lx2))/72.;
};
double C2TQ2QPS_Rx_Sz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2TQ2QB_Rx_Sz(int const&, double const& x, double const& z, bool const& intz)
{
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double NCi2 = NCi * NCi;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double x2 = x * x;
  const double xi  = 1. / x;
  const double opxi = 1. / ( 1 + x );
  const double lx  = log(x);
  const double lx2 = lx * lx;
  const double lomx  = log(1 - x);
  const double lopx  = log(1 + x);
  const double Li2x  = apfel::dilog(x);
  const double Li2mx = apfel::dilog(-x);
  return (3*Dn(1,1 - z, intz)*(NC - NCi)* 
      (3 - 3*x + 6*(1 + x)*lx + 4*xi - 4*x2) -  
     Dn(0,1 - z, intz)*(36 + 33*NC - 36*x - 51*NC*x - 9*NC*lomx +  
        9*NC*x*lomx + 18*lx - 18*NC*lx + 18*x*lx -  
        54*NC*x*lx + 36*lx*lopx - 36*x*lx*lopx -  
        72*lx*lopx*opxi +  
        36*Li2mx*(-1 + x + 2*opxi)*(-1 + NCi2) -  
        36*NCi2 + 36*x*NCi2 - 18*lx*NCi2 -  
        18*x*lx*NCi2 - 36*lx*lopx*NCi2 +  
        36*x*lx*lopx*NCi2 +  
        72*lx*lopx*opxi*NCi2 +  
        18*(1 + x)*Li2x*(NC - NCi) - 33*NCi +  
        51*x*NCi + 9*lomx*NCi -  
        9*x*lomx*NCi + 18*lx*NCi +  
        54*x*lx*NCi + 3*pi2 - 3*NC*pi2 -  
        3*x*pi2 - 3*NC*x*pi2 - 6*opxi*pi2 -  
        3*NCi2*pi2 + 3*x*NCi2*pi2 +  
        6*opxi*NCi2*pi2 + 3*NCi*pi2 +  
        3*x*NCi*pi2 - 26*NC*xi -  
        12*NC*lomx*xi + 26*NCi*xi +  
        12*lomx*NCi*xi + 44*NC*x2 +  
        12*NC*lomx*x2 - 36*NC*lx*x2 -  
        44*NCi*x2 - 12*lomx*NCi*x2 +  
        36*lx*NCi*x2 - 9*lx2 +  
        18*NC*lx2 + 9*x*lx2 + 18*NC*x*lx2 +  
        18*opxi*lx2 + 9*NCi2*lx2 -  
        9*x*NCi2*lx2 -  
        18*opxi*NCi2*lx2 -  
        18*NCi*lx2 - 18*x*NCi*lx2))/36.;
};
double C2TQ2QP1_Rx_Sz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
double C2TQ2QP2_Rx_Sz(int const&, double const& x, double const& z, bool const& intz)
{
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double x2 = x * x;
  const double xi  = 1. / x;
  const double lx  = log(x);
  const double lx2 = lx * lx;
  const double lomx  = log(1 - x);
  const double Li2x  = apfel::dilog(x);
  return -0.027777777777777776*((NC - NCi)* 
     (-3*Dn(1,1 - z, intz)*(3 - 3*x + 6*(1 + x)*lx + 4*xi -  
          4*x2) + Dn(0,1 - z, intz)* 
        (33 - 51*x + 18*(1 + x)*Li2x - 18*lx - 54*x*lx -  
          3*pi2 - 3*x*pi2 - 26*xi + 44*x2 -  
          36*lx*x2 + 3*lomx* 
           (-3 + 3*x - 4*xi + 4*x2) + 18*lx2 +  
          18*x*lx2)));
};
double C2TQ2QP3_Rx_Sz(int const&, double const&, double const&, bool const&)
{
  return 0;
};
//------------------------------------------------------
double C1LQ2Q_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double NC = apfel::NC;
  const double NCi = 1. / NC;
  return 2*x*z*(NC - NCi);
};
double C1LQ2G_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double NC = apfel::NC;
  const double NCi = 1. / NC;
  return -2*x*(-1 + z)*(NC - NCi);
};
double C1LG2Q_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double x2 = x * x;
  return 4*(x - x2);
};
double C1TQ2Q_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double NC = apfel::NC;
  const double NCi = 1. / NC;
  return (1 + x*z)*(NC - NCi);
};
double C1TQ2G_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double NC = apfel::NC;
  const double NCi = 1. / NC;
  const double zi  = 1. / z;
  return -0.5*((NC - NCi)*(-2 + 2*x*(-1 + z) + (1 + x)*zi));
};
double C1TG2Q_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double x2 = x * x;
  const double zi  = 1. / z;
  return -0.5*((-1 + 2*x - 2*x2)*(-2 + zi));
};
double C2LQ2G_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double sqrtxz1  = sqrt(1 - 2*z + z*z + 4*x*z);
  const double NC   = apfel::NC;
  const double NC2  = NC * NC;
  const double NCi  = 1. / NC;
  const double NCi2 = NCi * NCi;
  const double l2  = log(2);
  const double l22 = l2 * l2;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double x2 = x * x;
  const double x3 = x * x2;
  const double x4 = x * x3;
  const double xi  = 1. / x;
  const double z2 = z * z;
  const double zi  = 1. / z;
  const double omxi = 1. / ( 1 - x );
  const double omzi  = 1. / ( 1 - z );
  const double lx  = log(x);
  const double lx2 = lx * lx;
  const double lz  = log(z);
  const double lz2 = lz * lz;
  const double lomx  = log(1 - x);
  const double lomx2 = lomx * lomx;
  const double lomz  = log(1 - z);
  const double lomz2 = lomz * lomz;
  const double Li2x  = apfel::dilog(x);
  const double Li2z  = apfel::dilog(z);
  const double omxmzi  = 1. / ( 1 - x - z );
  const double omxmzi2 = omxmzi * omxmzi;
  const double lomxmz  = log(1 - x - z);
  const double lmopxpz = log(-1 + x + z);
  const double lxmz    = log(x - z);
  const double lmxpz   = log(-x + z);
  const double li2omxzi = apfel::dilog(1 - x*zi);
  const double lspec1   = log(1 + sqrtxz1 - z);
  const double lspec2   = log(1 + sqrtxz1 + z);
  const double li2spec1  = apfel::dilog(0.5 - sqrtxz1/2. - z/2.);
  const double li2spec2  = apfel::dilog(0.5 - sqrtxz1/2. + z/2.);
  const double li2spec3  = apfel::dilog(0.5 - zi/2. - (sqrtxz1*zi)/2.);
  const double li2spec4  = apfel::dilog(0.5 + zi/2. - (sqrtxz1*zi)/2.);
  const double li2spec9  = apfel::dilog((1 - z)*omxi);
  const double li2spec10 = apfel::dilog(x*(1 - z)*omxi*zi);
  const double li2spec11 = apfel::dilog((1 - x)*omzi);
  const double li2spec12 = apfel::dilog((1 - x)*z*xi*omzi);
  const double li2spec13 = apfel::dilog(z*omxi);
  const double li2spec14 = apfel::dilog((1 - z)*z*omxi*xi);
  const double li2spec15 = apfel::dilog(x*z*omxi*omzi);
  const double li2spec16 = apfel::dilog((1 - x)*zi);
  const double li2spec17 = apfel::dilog((1 - x)*(1 - z)*xi*zi);
  const double li2spec18 = apfel::dilog((1 - x)*x*omzi*zi);
  const double Tt1 = (z > x ? 1 : 0);
  const double Tt2 = (z < x ? 1 : 0);
  const double Tu1 = (z < 1 - x && z < x ? 1 : 0);
  const double Tu2 = (z > 1 - x && z < x ? 1 : 0);
  const double Tu3 = (z < 1 - x && z > x ? 1 : 0);
  const double Tu4 = (z > 1 - x && z > x ? 1 : 0);
  return -6 + 20*x + 5*z - 10*x*z + 2*x*li2spec1 +  
   2*x*z*li2spec1 - 2*x*li2spec2 -  
   2*x*z*li2spec2 + 2*x*Li2z + 10*x*z*Li2z -  
   2*x*li2spec3 -  
   2*x*z*li2spec3 +  
   2*x*li2spec4 +  
   2*x*z*li2spec4 +  
   4*x*li2omxzi - 4*x*z*li2omxzi -  
   2*sqrtxz1*x*l2 - 2*lomx - 3*x*lomx + 2*z*lomx +  
   9*x*z*lomx + 2*lx + 7*x*lx - sqrtxz1*x*lx +  
   4*(1 - x)*x*lx - 2*z*lx - 17*x*z*lx - 2*x*l2*lx -  
   2*x*z*l2*lx + 42*x*lomx*lx - 32*x*z*lomx*lx -  
   2*lomz + 3*x*lomz - 4*(1 - x)*x*lomz + 2*z*lomz +  
   5*x*z*lomz - 26*x*lomx*lomz +  
   20*x*z*lomx*lomz + 42*x*lx*lomz -  
   32*x*z*lx*lomz + 2*sqrtxz1*x*lspec1 +  
   4*x*l2*lspec1 + 4*x*z*l2*lspec1 -  
   6*lz - 8*x*lz - sqrtxz1*x*lz - 4*z*lz + 23*x*z*lz -  
   6*x*l2*lz - 6*x*z*l2*lz - 24*x*lomx*lz +  
   10*x*z*lomx*lz + 36*x*lx*lz - 16*x*z*lx*lz -  
   22*x*lomz*lz + 20*x*z*lomz*lz +  
   2*x*lspec1*lz + 2*x*z*lspec1*lz +  
   4*x*l2*lspec2 + 4*x*z*l2*lspec2 +  
   2*x*lx*lspec2 + 2*x*z*lx*lspec2 -  
   4*x*lspec1*lspec2 -  
   4*x*z*lspec1*lspec2 +  
   4*x*lz*lspec2 + 4*x*z*lz*lspec2 +  
   NCi2 - (17*x*NCi2)/2. - (3*z*NCi2)/2. +  
   7*x*z*NCi2 + x*Li2x*NCi2 - x*z*Li2x*NCi2 -  
   2*x*Li2z*NCi2 + x*z*Li2z*NCi2 -  
   2*x*li2omxzi*NCi2 +  
   2*x*z*li2omxzi*NCi2 - 2*sqrtxz1*x*l2*NCi2 +  
   lomx*NCi2 + (7*x*lomx*NCi2)/2. -  
   z*lomx*NCi2 - (9*x*z*lomx*NCi2)/2. -  
   lx*NCi2 - (7*x*lx*NCi2)/2. -  
   sqrtxz1*x*lx*NCi2 - 4*(1 - x)*x*lx*NCi2 +  
   z*lx*NCi2 + 8*x*z*lx*NCi2 -  
   28*x*lomx*lx*NCi2 + 18*x*z*lomx*lx*NCi2 +  
   lomz*NCi2 - (x*lomz*NCi2)/2. +  
   4*(1 - x)*x*lomz*NCi2 - z*lomz*NCi2 -  
   (9*x*z*lomz*NCi2)/2. + 18*x*lomx*lomz*NCi2 -  
   12*x*z*lomx*lomz*NCi2 -  
   31*x*lx*lomz*NCi2 + 21*x*z*lx*lomz*NCi2 +  
   2*sqrtxz1*x*lspec1*NCi2 + 2*lz*NCi2 +  
   6*x*lz*NCi2 - sqrtxz1*x*lz*NCi2 -  
   9*x*z*lz*NCi2 + 12*x*lomx*lz*NCi2 -  
   9*x*z*lomx*lz*NCi2 - 20*x*lx*lz*NCi2 +  
   17*x*z*lx*lz*NCi2 + 12*x*lomz*lz*NCi2 -  
   10*x*z*lomz*lz*NCi2 + 5*NC2 - (23*x*NC2)/2. -  
   (7*z*NC2)/2. + 3*x*z*NC2 - x*Li2x*NC2 +  
   x*z*Li2x*NC2 - 2*x*li2spec1*NC2 -  
   2*x*z*li2spec1*NC2 +  
   2*x*li2spec2*NC2 +  
   2*x*z*li2spec2*NC2 - 11*x*z*Li2z*NC2 +  
   2*x*li2spec3*NC2 +  
   2*x*z*li2spec3*NC2 -  
   2*x*li2spec4*NC2 -  
   2*x*z*li2spec4*NC2 -  
   2*x*li2omxzi*NC2 +  
   2*x*z*li2omxzi*NC2 + 4*sqrtxz1*x*l2*NC2 +  
   lomx*NC2 - (x*lomx*NC2)/2. -  
   z*lomx*NC2 - (9*x*z*lomx*NC2)/2. -  
   lx*NC2 - (7*x*lx*NC2)/2. +  
   2*sqrtxz1*x*lx*NC2 + z*lx*NC2 +  
   9*x*z*lx*NC2 + 2*x*l2*lx*NC2 +  
   2*x*z*l2*lx*NC2 - 14*x*lomx*lx*NC2 +  
   14*x*z*lomx*lx*NC2 + lomz*NC2 -  
   (5*x*lomz*NC2)/2. - z*lomz*NC2 -  
   (x*z*lomz*NC2)/2. + 8*x*lomx*lomz*NC2 -  
   8*x*z*lomx*lomz*NC2 - 11*x*lx*lomz*NC2 +  
   11*x*z*lx*lomz*NC2 -  
   4*sqrtxz1*x*lspec1*NC2 -  
   4*x*l2*lspec1*NC2 -  
   4*x*z*l2*lspec1*NC2 + 4*lz*NC2 +  
   2*x*lz*NC2 + 2*sqrtxz1*x*lz*NC2 +  
   4*z*lz*NC2 - 14*x*z*lz*NC2 +  
   6*x*l2*lz*NC2 + 6*x*z*l2*lz*NC2 +  
   12*x*lomx*lz*NC2 - x*z*lomx*lz*NC2 -  
   16*x*lx*lz*NC2 - x*z*lx*lz*NC2 +  
   10*x*lomz*lz*NC2 - 10*x*z*lomz*lz*NC2 -  
   2*x*lspec1*lz*NC2 -  
   2*x*z*lspec1*lz*NC2 -  
   4*x*l2*lspec2*NC2 -  
   4*x*z*l2*lspec2*NC2 -  
   2*x*lx*lspec2*NC2 -  
   2*x*z*lx*lspec2*NC2 +  
   4*x*lspec1*lspec2*NC2 +  
   4*x*z*lspec1*lspec2*NC2 -  
   4*x*lz*lspec2*NC2 -  
   4*x*z*lz*lspec2*NC2 + (13*x*pi2)/3. -  
   (14*x*z*pi2)/3. - (23*x*NCi2*pi2)/6. +  
   (7*x*z*NCi2*pi2)/3. - (x*NC2*pi2)/2. +  
   (7*x*z*NC2*pi2)/3. - 2*x2 + 2*lx*x2 -  
   2*lomz*x2 + 2*NCi2*x2 -  
   2*lx*NCi2*x2 + 2*lomz*NCi2*x2 +  
   2*lx*x3*omxmzi2 -  
   2*lomz*x3*omxmzi2 -  
   2*lx*NCi2*x3*omxmzi2 +  
   2*lomz*NCi2*x3*omxmzi2 -  
   2*lx*x4*omxmzi2 +  
   2*lomz*x4*omxmzi2 +  
   2*lx*NCi2*x4*omxmzi2 -  
   2*lomz*NCi2*x4*omxmzi2 +  
   2*x2*omxmzi + 2*lx*x2*omxmzi -  
   2*lomz*x2*omxmzi -  
   2*NCi2*x2*omxmzi -  
   2*lx*NCi2*x2*omxmzi +  
   2*lomz*NCi2*x2*omxmzi -  
   2*x3*omxmzi - 4*lx*x3*omxmzi +  
   4*lomz*x3*omxmzi +  
   2*NCi2*x3*omxmzi +  
   4*lx*NCi2*x3*omxmzi -  
   4*lomz*NCi2*x3*omxmzi - zi -  
   3*x*zi - 2*x*Li2z*zi + 2*sqrtxz1*x*l2*zi -  
   4*x*lomx*zi + 4*x*lx*zi +  
   sqrtxz1*x*lx*zi - 10*x*lomx*lx*zi -  
   4*x*lomz*zi + 6*x*lomx*lomz*zi -  
   10*x*lx*lomz*zi -  
   2*sqrtxz1*x*lspec1*zi - 7*x*lz*zi +  
   sqrtxz1*x*lz*zi + 4*x*lomx*lz*zi -  
   6*x*lx*lz*zi + 2*x*lomz*lz*zi +  
   (NCi2*zi)/2. - (x*NCi2*zi)/2. +  
   2*x*Li2z*NCi2*zi -  
   2*sqrtxz1*x*l2*NCi2*zi +  
   x*lomx*NCi2*zi - (x*lx*NCi2*zi)/2. -  
   sqrtxz1*x*lx*NCi2*zi +  
   10*x*lomx*lx*NCi2*zi +  
   x*lomz*NCi2*zi -  
   6*x*lomx*lomz*NCi2*zi +  
   10*x*lx*lomz*NCi2*zi +  
   2*sqrtxz1*x*lspec1*NCi2*zi +  
   2*x*lz*NCi2*zi - sqrtxz1*x*lz*NCi2*zi -  
   4*x*lomx*lz*NCi2*zi +  
   6*x*lx*lz*NCi2*zi -  
   2*x*lomz*lz*NCi2*zi + (NC2*zi)/2. +  
   (7*x*NC2*zi)/2. + 3*x*lomx*NC2*zi -  
   (7*x*lx*NC2*zi)/2. + 3*x*lomz*NC2*zi +  
   5*x*lz*NC2*zi - (4*x*pi2*zi)/3. +  
   (4*x*NCi2*pi2*zi)/3. + 2*z2 - 5*x*z2 +  
   2*x*lomx*z2 - 4*x*lx*z2 + 2*x*lomz*z2 +  
   2*x*lz*z2 - 2*NC2*z2 + 5*x*NC2*z2 -  
   2*x*lomx*NC2*z2 + 4*x*lx*NC2*z2 -  
   2*x*lomz*NC2*z2 - 2*x*lz*NC2*z2 -  
   4*x*l22 - 4*x*z*l22 + 4*x*NC2*l22 +  
   4*x*z*NC2*l22 - 15*x*lomx2 +  
   11*x*z*lomx2 + 10*x*NCi2*lomx2 -  
   6*x*z*NCi2*lomx2 + 5*x*NC2*lomx2 -  
   5*x*z*NC2*lomx2 + 4*x*zi*lomx2 -  
   4*x*NCi2*zi*lomx2 - 30*x*lx2 +  
   23*x*z*lx2 + 22*x*NCi2*lx2 -  
   15*x*z*NCi2*lx2 + 8*x*NC2*lx2 -  
   8*x*z*NC2*lx2 + 7*x*zi*lx2 -  
   7*x*NCi2*zi*lx2 - 15*x*lomz2 +  
   11*x*z*lomz2 + 12*x*NCi2*lomz2 -  
   8*x*z*NCi2*lomz2 + 3*x*NC2*lomz2 -  
   3*x*z*NC2*lomz2 + 4*x*zi*lomz2 -  
   4*x*NCi2*zi*lomz2 - 13*x*lz2 -  
   2*x*z*lz2 + 5*x*NCi2*lz2 -  
   5*x*z*NCi2*lz2 + 8*x*NC2*lz2 +  
   7*x*z*NC2*lz2 + x*zi*lz2 -  
   x*NCi2*zi*lz2 +  
   (-2*x*Li2z + 2*x*z*Li2z - 2*x*li2spec9 +  
      2*x*z*li2spec9 +  
      2*x*li2spec10 -  
      2*x*z*li2spec10 - 4*x*z*lomx +  
      6*x*z*lx - 12*x*lomx*lx + 12*x*z*lomx*lx -  
      2*x*z*lomz + 4*x*lomx*lomz -  
      4*x*z*lomx*lomz - 6*x*lx*lomz +  
      6*x*z*lx*lomz - 4*x*z*lz + 10*x*lomx*lz -  
      10*x*z*lomx*lz - 12*x*lx*lz + 12*x*z*lx*lz +  
      2*x*lomz*lz - 2*x*z*lomz*lz +  
      2*x*lx*lmxpz - 2*x*z*lx*lmxpz -  
      2*x*lz*lmxpz + 2*x*z*lz*lmxpz +  
      2*x*Li2z*NC2 - 2*x*z*Li2z*NC2 +  
      2*x*li2spec9*NC2 -  
      2*x*z*li2spec9*NC2 -  
      2*x*li2spec10*NC2 +  
      2*x*z*li2spec10*NC2 +  
      4*x*z*lomx*NC2 - 6*x*z*lx*NC2 +  
      12*x*lomx*lx*NC2 - 12*x*z*lomx*lx*NC2 +  
      2*x*z*lomz*NC2 - 4*x*lomx*lomz*NC2 +  
      4*x*z*lomx*lomz*NC2 +  
      6*x*lx*lomz*NC2 - 6*x*z*lx*lomz*NC2 +  
      4*x*z*lz*NC2 - 10*x*lomx*lz*NC2 +  
      10*x*z*lomx*lz*NC2 + 12*x*lx*lz*NC2 -  
      12*x*z*lx*lz*NC2 - 2*x*lomz*lz*NC2 +  
      2*x*z*lomz*lz*NC2 - 2*x*lx*lmxpz*NC2 +  
      2*x*z*lx*lmxpz*NC2 + 2*x*lz*lmxpz*NC2 -  
      2*x*z*lz*lmxpz*NC2 - x*pi2 + x*z*pi2 +  
      x*NC2*pi2 - x*z*NC2*pi2 +  
      4*x*lomx2 - 4*x*z*lomx2 -  
      4*x*NC2*lomx2 + 4*x*z*NC2*lomx2 +  
      7*x*lx2 - 7*x*z*lx2 - 7*x*NC2*lx2 +  
      7*x*z*NC2*lx2 + x*lomz2 -  
      x*z*lomz2 - x*NC2*lomz2 +  
      x*z*NC2*lomz2 + 5*x*lz2 -  
      5*x*z*lz2 - 5*x*NC2*lz2 +  
      5*x*z*NC2*lz2)*Tt1 +  
   (-2*x*Li2z + 2*x*z*Li2z + 2*x*li2spec11 -  
      2*x*z*li2spec11 -  
      2*x*li2spec12 +  
      2*x*z*li2spec12 - 4*x*z*lomx +  
      6*x*z*lx - 10*x*lomx*lx + 10*x*z*lomx*lx -  
      2*x*z*lomz + 4*x*lomx*lomz -  
      4*x*z*lomx*lomz - 8*x*lx*lomz +  
      8*x*z*lx*lomz + 2*x*lx*lxmz -  
      2*x*z*lx*lxmz - 4*x*z*lz + 8*x*lomx*lz -  
      8*x*z*lomx*lz - 10*x*lx*lz + 10*x*z*lx*lz +  
      4*x*lomz*lz - 4*x*z*lomz*lz -  
      2*x*lxmz*lz + 2*x*z*lxmz*lz +  
      2*x*Li2z*NC2 - 2*x*z*Li2z*NC2 -  
      2*x*li2spec11*NC2 +  
      2*x*z*li2spec11*NC2 +  
      2*x*li2spec12*NC2 -  
      2*x*z*li2spec12*NC2 +  
      4*x*z*lomx*NC2 - 6*x*z*lx*NC2 +  
      10*x*lomx*lx*NC2 - 10*x*z*lomx*lx*NC2 +  
      2*x*z*lomz*NC2 - 4*x*lomx*lomz*NC2 +  
      4*x*z*lomx*lomz*NC2 +  
      8*x*lx*lomz*NC2 - 8*x*z*lx*lomz*NC2 -  
      2*x*lx*lxmz*NC2 + 2*x*z*lx*lxmz*NC2 +  
      4*x*z*lz*NC2 - 8*x*lomx*lz*NC2 +  
      8*x*z*lomx*lz*NC2 + 10*x*lx*lz*NC2 -  
      10*x*z*lx*lz*NC2 - 4*x*lomz*lz*NC2 +  
      4*x*z*lomz*lz*NC2 + 2*x*lxmz*lz*NC2 -  
      2*x*z*lxmz*lz*NC2 - x*pi2 + x*z*pi2 +  
      x*NC2*pi2 - x*z*NC2*pi2 +  
      4*x*lomx2 - 4*x*z*lomx2 -  
      4*x*NC2*lomx2 + 4*x*z*NC2*lomx2 +  
      6*x*lx2 - 6*x*z*lx2 - 6*x*NC2*lx2 +  
      6*x*z*NC2*lx2 + x*lomz2 -  
      x*z*lomz2 - x*NC2*lomz2 +  
      x*z*NC2*lomz2 + 4*x*lz2 -  
      4*x*z*lz2 - 4*x*NC2*lz2 +  
      4*x*z*NC2*lz2)*Tt2 +  
   (4*x*Li2z - 2*x*z*Li2z - 2*x*li2spec13 -  
      2*x*li2spec14 +  
      2*x*z*li2spec14 -  
      2*x*li2spec11 + 2*x*z*li2spec11 +  
      4*x*li2spec15 -  
      2*x*z*li2spec15 + 6*x*lomx -  
      6*x*z*lomx - 12*x*lx + 12*x*z*lx -  
      32*x*lomx*lx + 20*x*z*lomx*lx + 8*x*lomz -  
      8*x*z*lomz + 26*x*lomx*lomz -  
      16*x*z*lomx*lomz - 34*x*lx*lomz +  
      22*x*z*lx*lomz + 6*x*lx*lomxmz -  
      4*x*z*lx*lomxmz - 6*x*lomz*lomxmz +  
      4*x*z*lomz*lomxmz + 2*x*lx*lxmz -  
      2*x*z*lx*lxmz + 6*x*lz - 6*x*z*lz +  
      12*x*lomx*lz - 8*x*z*lomx*lz - 22*x*lx*lz +  
      16*x*z*lx*lz + 16*x*lomz*lz -  
      12*x*z*lomz*lz - 2*x*lxmz*lz +  
      2*x*z*lxmz*lz - 4*x*Li2z*NCi2 +  
      2*x*z*Li2z*NCi2 + 2*x*li2spec13*NCi2 +  
      2*x*li2spec14*NCi2 -  
      2*x*z*li2spec14*NCi2 +  
      2*x*li2spec11*NCi2 -  
      2*x*z*li2spec11*NCi2 -  
      4*x*li2spec15*NCi2 +  
      2*x*z*li2spec15*NCi2 -  
      6*x*lomx*NCi2 + 6*x*z*lomx*NCi2 +  
      12*x*lx*NCi2 - 12*x*z*lx*NCi2 +  
      32*x*lomx*lx*NCi2 -  
      20*x*z*lomx*lx*NCi2 - 8*x*lomz*NCi2 +  
      8*x*z*lomz*NCi2 - 26*x*lomx*lomz*NCi2 +  
      16*x*z*lomx*lomz*NCi2 +  
      34*x*lx*lomz*NCi2 -  
      22*x*z*lx*lomz*NCi2 -  
      6*x*lx*lomxmz*NCi2 +  
      4*x*z*lx*lomxmz*NCi2 +  
      6*x*lomz*lomxmz*NCi2 -  
      4*x*z*lomz*lomxmz*NCi2 -  
      2*x*lx*lxmz*NCi2 + 2*x*z*lx*lxmz*NCi2 -  
      6*x*lz*NCi2 + 6*x*z*lz*NCi2 -  
      12*x*lomx*lz*NCi2 +  
      8*x*z*lomx*lz*NCi2 + 22*x*lx*lz*NCi2 -  
      16*x*z*lx*lz*NCi2 - 16*x*lomz*lz*NCi2 +  
      12*x*z*lomz*lz*NCi2 +  
      2*x*lxmz*lz*NCi2 - 2*x*z*lxmz*lz*NCi2 -  
      3*x*pi2 + (5*x*z*pi2)/3. + 3*x*NCi2*pi2 -  
      (5*x*z*NCi2*pi2)/3. - 2*x*Li2z*zi +  
      2*x*li2spec13*zi -  
      2*x*li2spec15*zi +  
      12*x*lomx*lx*zi -  
      10*x*lomx*lomz*zi +  
      12*x*lx*lomz*zi -  
      2*x*lx*lomxmz*zi +  
      2*x*lomz*lomxmz*zi -  
      4*x*lomx*lz*zi + 6*x*lx*lz*zi -  
      4*x*lomz*lz*zi + 2*x*Li2z*NCi2*zi -  
      2*x*li2spec13*NCi2*zi +  
      2*x*li2spec15*NCi2*zi -  
      12*x*lomx*lx*NCi2*zi +  
      10*x*lomx*lomz*NCi2*zi -  
      12*x*lx*lomz*NCi2*zi +  
      2*x*lx*lomxmz*NCi2*zi -  
      2*x*lomz*lomxmz*NCi2*zi +  
      4*x*lomx*lz*NCi2*zi -  
      6*x*lx*lz*NCi2*zi +  
      4*x*lomz*lz*NCi2*zi +  
      (4*x*pi2*zi)/3. -  
      (4*x*NCi2*pi2*zi)/3. + 9*x*lomx2 -  
      5*x*z*lomx2 - 9*x*NCi2*lomx2 +  
      5*x*z*NCi2*lomx2 - 4*x*zi*lomx2 +  
      4*x*NCi2*zi*lomx2 + 20*x*lx2 -  
      13*x*z*lx2 - 20*x*NCi2*lx2 +  
      13*x*z*NCi2*lx2 - 7*x*zi*lx2 +  
      7*x*NCi2*zi*lx2 + 13*x*lomz2 -  
      8*x*z*lomz2 - 13*x*NCi2*lomz2 +  
      8*x*z*NCi2*lomz2 - 5*x*zi*lomz2 +  
      5*x*NCi2*zi*lomz2 + 6*x*lz2 -  
      5*x*z*lz2 - 6*x*NCi2*lz2 +  
      5*x*z*NCi2*lz2 - x*zi*lz2 +  
      x*NCi2*zi*lz2)*Tu1 +  
   (4*x*Li2z - 2*x*z*Li2z - 2*x*li2spec11 +  
      2*x*z*li2spec11 + 2*x*li2spec16 -  
      4*x*li2spec17 +  
      2*x*z*li2spec17 +  
      2*x*li2spec18 -  
      2*x*z*li2spec18 + 6*x*lomx -  
      6*x*z*lomx - 12*x*lx + 12*x*z*lx -  
      26*x*lomx*lx + 16*x*z*lomx*lx + 8*x*lomz -  
      8*x*z*lomz + 20*x*lomx*lomz -  
      12*x*z*lomx*lomz - 32*x*lx*lomz +  
      22*x*z*lx*lomz + 2*x*lx*lxmz -  
      2*x*z*lx*lxmz + 6*x*lz - 6*x*z*lz +  
      12*x*lomx*lz - 8*x*z*lomx*lz - 28*x*lx*lz +  
      20*x*z*lx*lz + 22*x*lomz*lz -  
      16*x*z*lomz*lz - 2*x*lxmz*lz +  
      2*x*z*lxmz*lz + 6*x*lx*lmopxpz -  
      4*x*z*lx*lmopxpz - 6*x*lomz*lmopxpz +  
      4*x*z*lomz*lmopxpz - 4*x*Li2z*NCi2 +  
      2*x*z*Li2z*NCi2 + 2*x*li2spec11*NCi2 -  
      2*x*z*li2spec11*NCi2 -  
      2*x*li2spec16*NCi2 +  
      4*x*li2spec17*NCi2 -  
      2*x*z*li2spec17*NCi2 -  
      2*x*li2spec18*NCi2 +  
      2*x*z*li2spec18*NCi2 -  
      6*x*lomx*NCi2 + 6*x*z*lomx*NCi2 +  
      12*x*lx*NCi2 - 12*x*z*lx*NCi2 +  
      26*x*lomx*lx*NCi2 -  
      16*x*z*lomx*lx*NCi2 - 8*x*lomz*NCi2 +  
      8*x*z*lomz*NCi2 - 20*x*lomx*lomz*NCi2 +  
      12*x*z*lomx*lomz*NCi2 +  
      32*x*lx*lomz*NCi2 -  
      22*x*z*lx*lomz*NCi2 -  
      2*x*lx*lxmz*NCi2 + 2*x*z*lx*lxmz*NCi2 -  
      6*x*lz*NCi2 + 6*x*z*lz*NCi2 -  
      12*x*lomx*lz*NCi2 +  
      8*x*z*lomx*lz*NCi2 + 28*x*lx*lz*NCi2 -  
      20*x*z*lx*lz*NCi2 - 22*x*lomz*lz*NCi2 +  
      16*x*z*lomz*lz*NCi2 +  
      2*x*lxmz*lz*NCi2 - 2*x*z*lxmz*lz*NCi2 -  
      6*x*lx*lmopxpz*NCi2 +  
      4*x*z*lx*lmopxpz*NCi2 +  
      6*x*lomz*lmopxpz*NCi2 -  
      4*x*z*lomz*lmopxpz*NCi2 - 3*x*pi2 +  
      (5*x*z*pi2)/3. + 3*x*NCi2*pi2 -  
      (5*x*z*NCi2*pi2)/3. - 2*x*Li2z*zi -  
      2*x*li2spec16*zi +  
      2*x*li2spec17*zi +  
      10*x*lomx*lx*zi -  
      8*x*lomx*lomz*zi +  
      10*x*lx*lomz*zi - 4*x*lomx*lz*zi +  
      8*x*lx*lz*zi - 6*x*lomz*lz*zi -  
      2*x*lx*lmopxpz*zi +  
      2*x*lomz*lmopxpz*zi +  
      2*x*Li2z*NCi2*zi +  
      2*x*li2spec16*NCi2*zi -  
      2*x*li2spec17*NCi2*zi -  
      10*x*lomx*lx*NCi2*zi +  
      8*x*lomx*lomz*NCi2*zi -  
      10*x*lx*lomz*NCi2*zi +  
      4*x*lomx*lz*NCi2*zi -  
      8*x*lx*lz*NCi2*zi +  
      6*x*lomz*lz*NCi2*zi +  
      2*x*lx*lmopxpz*NCi2*zi -  
      2*x*lomz*lmopxpz*NCi2*zi +  
      (4*x*pi2*zi)/3. -  
      (4*x*NCi2*pi2*zi)/3. + 9*x*lomx2 -  
      5*x*z*lomx2 - 9*x*NCi2*lomx2 +  
      5*x*z*NCi2*lomx2 - 4*x*zi*lomx2 +  
      4*x*NCi2*zi*lomx2 + 19*x*lx2 -  
      13*x*z*lx2 - 19*x*NCi2*lx2 +  
      13*x*z*NCi2*lx2 - 6*x*zi*lx2 +  
      6*x*NCi2*zi*lx2 + 12*x*lomz2 -  
      8*x*z*lomz2 - 12*x*NCi2*lomz2 +  
      8*x*z*NCi2*lomz2 - 4*x*zi*lomz2 +  
      4*x*NCi2*zi*lomz2 + 6*x*lz2 -  
      5*x*z*lz2 - 6*x*NCi2*lz2 +  
      5*x*z*NCi2*lz2 - x*zi*lz2 +  
      x*NCi2*zi*lz2)*Tu2 +  
   (4*x*Li2z - 2*x*z*Li2z + 2*x*li2spec9 -  
      2*x*z*li2spec9 - 2*x*li2spec13 +  
      4*x*li2spec15 -  
      2*x*z*li2spec15 +  
      2*x*li2spec18 -  
      2*x*z*li2spec18 + 6*x*lomx -  
      6*x*z*lomx - 12*x*lx + 12*x*z*lx -  
      30*x*lomx*lx + 18*x*z*lomx*lx + 8*x*lomz -  
      8*x*z*lomz + 22*x*lomx*lomz -  
      12*x*z*lomx*lomz - 36*x*lx*lomz +  
      24*x*z*lx*lomz + 6*x*lx*lomxmz -  
      4*x*z*lx*lomxmz - 6*x*lomz*lomxmz +  
      4*x*z*lomz*lomxmz + 6*x*lz - 6*x*z*lz +  
      10*x*lomx*lz - 6*x*z*lomx*lz - 24*x*lx*lz +  
      18*x*z*lx*lz + 18*x*lomz*lz -  
      14*x*z*lomz*lz + 2*x*lx*lmxpz -  
      2*x*z*lx*lmxpz - 2*x*lz*lmxpz +  
      2*x*z*lz*lmxpz - 4*x*Li2z*NCi2 +  
      2*x*z*Li2z*NCi2 - 2*x*li2spec9*NCi2 +  
      2*x*z*li2spec9*NCi2 +  
      2*x*li2spec13*NCi2 -  
      4*x*li2spec15*NCi2 +  
      2*x*z*li2spec15*NCi2 -  
      2*x*li2spec18*NCi2 +  
      2*x*z*li2spec18*NCi2 -  
      6*x*lomx*NCi2 + 6*x*z*lomx*NCi2 +  
      12*x*lx*NCi2 - 12*x*z*lx*NCi2 +  
      30*x*lomx*lx*NCi2 -  
      18*x*z*lomx*lx*NCi2 - 8*x*lomz*NCi2 +  
      8*x*z*lomz*NCi2 - 22*x*lomx*lomz*NCi2 +  
      12*x*z*lomx*lomz*NCi2 +  
      36*x*lx*lomz*NCi2 -  
      24*x*z*lx*lomz*NCi2 -  
      6*x*lx*lomxmz*NCi2 +  
      4*x*z*lx*lomxmz*NCi2 +  
      6*x*lomz*lomxmz*NCi2 -  
      4*x*z*lomz*lomxmz*NCi2 - 6*x*lz*NCi2 +  
      6*x*z*lz*NCi2 - 10*x*lomx*lz*NCi2 +  
      6*x*z*lomx*lz*NCi2 + 24*x*lx*lz*NCi2 -  
      18*x*z*lx*lz*NCi2 - 18*x*lomz*lz*NCi2 +  
      14*x*z*lomz*lz*NCi2 -  
      2*x*lx*lmxpz*NCi2 +  
      2*x*z*lx*lmxpz*NCi2 +  
      2*x*lz*lmxpz*NCi2 -  
      2*x*z*lz*lmxpz*NCi2 - (13*x*pi2)/3. +  
      3*x*z*pi2 + (13*x*NCi2*pi2)/3. -  
      3*x*z*NCi2*pi2 - 2*x*Li2z*zi +  
      2*x*li2spec13*zi -  
      2*x*li2spec15*zi +  
      12*x*lomx*lx*zi -  
      10*x*lomx*lomz*zi +  
      12*x*lx*lomz*zi -  
      2*x*lx*lomxmz*zi +  
      2*x*lomz*lomxmz*zi -  
      4*x*lomx*lz*zi + 6*x*lx*lz*zi -  
      4*x*lomz*lz*zi + 2*x*Li2z*NCi2*zi -  
      2*x*li2spec13*NCi2*zi +  
      2*x*li2spec15*NCi2*zi -  
      12*x*lomx*lx*NCi2*zi +  
      10*x*lomx*lomz*NCi2*zi -  
      12*x*lx*lomz*NCi2*zi +  
      2*x*lx*lomxmz*NCi2*zi -  
      2*x*lomz*lomxmz*NCi2*zi +  
      4*x*lomx*lz*NCi2*zi -  
      6*x*lx*lz*NCi2*zi +  
      4*x*lomz*lz*NCi2*zi +  
      (4*x*pi2*zi)/3. -  
      (4*x*NCi2*pi2*zi)/3. + 11*x*lomx2 -  
      7*x*z*lomx2 - 11*x*NCi2*lomx2 +  
      7*x*z*NCi2*lomx2 - 4*x*zi*lomx2 +  
      4*x*NCi2*zi*lomx2 + 21*x*lx2 -  
      14*x*z*lx2 - 21*x*NCi2*lx2 +  
      14*x*z*NCi2*lx2 - 7*x*zi*lx2 +  
      7*x*NCi2*zi*lx2 + 15*x*lomz2 -  
      10*x*z*lomz2 - 15*x*NCi2*lomz2 +  
      10*x*z*NCi2*lomz2 - 5*x*zi*lomz2 +  
      5*x*NCi2*zi*lomz2 + 7*x*lz2 -  
      6*x*z*lz2 - 7*x*NCi2*lz2 +  
      6*x*z*NCi2*lz2 - x*zi*lz2 +  
      x*NCi2*zi*lz2)*Tu3 +  
   (4*x*Li2z - 2*x*z*Li2z + 2*x*li2spec9 -  
      2*x*z*li2spec9 -  
      2*x*li2spec14 +  
      2*x*z*li2spec14 +  
      2*x*li2spec16 -  
      4*x*li2spec17 +  
      2*x*z*li2spec17 + 6*x*lomx -  
      6*x*z*lomx - 12*x*lx + 12*x*z*lx -  
      28*x*lomx*lx + 18*x*z*lomx*lx + 8*x*lomz -  
      8*x*z*lomz + 20*x*lomx*lomz -  
      12*x*z*lomx*lomz - 30*x*lx*lomz +  
      20*x*z*lx*lomz + 6*x*lz - 6*x*z*lz +  
      14*x*lomx*lz - 10*x*z*lomx*lz -  
      26*x*lx*lz + 18*x*z*lx*lz + 20*x*lomz*lz -  
      14*x*z*lomz*lz + 2*x*lx*lmxpz -  
      2*x*z*lx*lmxpz - 2*x*lz*lmxpz +  
      2*x*z*lz*lmxpz + 6*x*lx*lmopxpz -  
      4*x*z*lx*lmopxpz - 6*x*lomz*lmopxpz +  
      4*x*z*lomz*lmopxpz - 4*x*Li2z*NCi2 +  
      2*x*z*Li2z*NCi2 - 2*x*li2spec9*NCi2 +  
      2*x*z*li2spec9*NCi2 +  
      2*x*li2spec14*NCi2 -  
      2*x*z*li2spec14*NCi2 -  
      2*x*li2spec16*NCi2 +  
      4*x*li2spec17*NCi2 -  
      2*x*z*li2spec17*NCi2 -  
      6*x*lomx*NCi2 + 6*x*z*lomx*NCi2 +  
      12*x*lx*NCi2 - 12*x*z*lx*NCi2 +  
      28*x*lomx*lx*NCi2 -  
      18*x*z*lomx*lx*NCi2 - 8*x*lomz*NCi2 +  
      8*x*z*lomz*NCi2 - 20*x*lomx*lomz*NCi2 +  
      12*x*z*lomx*lomz*NCi2 +  
      30*x*lx*lomz*NCi2 -  
      20*x*z*lx*lomz*NCi2 - 6*x*lz*NCi2 +  
      6*x*z*lz*NCi2 - 14*x*lomx*lz*NCi2 +  
      10*x*z*lomx*lz*NCi2 + 26*x*lx*lz*NCi2 -  
      18*x*z*lx*lz*NCi2 - 20*x*lomz*lz*NCi2 +  
      14*x*z*lomz*lz*NCi2 -  
      2*x*lx*lmxpz*NCi2 +  
      2*x*z*lx*lmxpz*NCi2 +  
      2*x*lz*lmxpz*NCi2 -  
      2*x*z*lz*lmxpz*NCi2 -  
      6*x*lx*lmopxpz*NCi2 +  
      4*x*z*lx*lmopxpz*NCi2 +  
      6*x*lomz*lmopxpz*NCi2 -  
      4*x*z*lomz*lmopxpz*NCi2 - 3*x*pi2 +  
      (5*x*z*pi2)/3. + 3*x*NCi2*pi2 -  
      (5*x*z*NCi2*pi2)/3. - 2*x*Li2z*zi -  
      2*x*li2spec16*zi +  
      2*x*li2spec17*zi +  
      10*x*lomx*lx*zi -  
      8*x*lomx*lomz*zi +  
      10*x*lx*lomz*zi - 4*x*lomx*lz*zi +  
      8*x*lx*lz*zi - 6*x*lomz*lz*zi -  
      2*x*lx*lmopxpz*zi +  
      2*x*lomz*lmopxpz*zi +  
      2*x*Li2z*NCi2*zi +  
      2*x*li2spec16*NCi2*zi -  
      2*x*li2spec17*NCi2*zi -  
      10*x*lomx*lx*NCi2*zi +  
      8*x*lomx*lomz*NCi2*zi -  
      10*x*lx*lomz*NCi2*zi +  
      4*x*lomx*lz*NCi2*zi -  
      8*x*lx*lz*NCi2*zi +  
      6*x*lomz*lz*NCi2*zi +  
      2*x*lx*lmopxpz*NCi2*zi -  
      2*x*lomz*lmopxpz*NCi2*zi +  
      (4*x*pi2*zi)/3. -  
      (4*x*NCi2*pi2*zi)/3. + 9*x*lomx2 -  
      5*x*z*lomx2 - 9*x*NCi2*lomx2 +  
      5*x*z*NCi2*lomx2 - 4*x*zi*lomx2 +  
      4*x*NCi2*zi*lomx2 + 18*x*lx2 -  
      12*x*z*lx2 - 18*x*NCi2*lx2 +  
      12*x*z*NCi2*lx2 - 6*x*zi*lx2 +  
      6*x*NCi2*zi*lx2 + 12*x*lomz2 -  
      8*x*z*lomz2 - 12*x*NCi2*lomz2 +  
      8*x*z*NCi2*lomz2 - 4*x*zi*lomz2 +  
      4*x*NCi2*zi*lomz2 + 5*x*lz2 -  
      4*x*z*lz2 - 5*x*NCi2*lz2 +  
      4*x*z*NCi2*lz2 - x*zi*lz2 +  
      x*NCi2*zi*lz2)*Tu4;
};
double C2LG2G_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double sqrtxz1  = sqrt(1 - 2*z + z*z + 4*x*z);
  const double sqrtxz3  = sqrt(x/z);
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double l2  = log(2);
  const double l22 = l2 * l2;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double x2 = x * x;
  const double xi  = 1. / x;
  const double xi2 = xi * xi;
  const double z2 = z * z;
  const double zi  = 1. / z;
  const double zi2 = zi * zi;
  const double opxi = 1. / ( 1 + x );
  const double lx  = log(x);
  const double lx2 = lx * lx;
  const double lz  = log(z);
  const double lz2 = lz * lz;
  const double lomx  = log(1 - x);
  const double lomz  = log(1 - z);
  const double lopx  = log(1 + x);
  const double Li2x  = apfel::dilog(x);
  const double Li2mx = apfel::dilog(-x);
  const double Li2z  = apfel::dilog(z);
  const double lopxz   = log(1 + x*z);
  const double lxpz    = log(x + z);
  const double lopxzi  = log(1 + x*zi);
  const double lspec1   = log(1 + sqrtxz1 - z);
  const double lspec1_2 = lspec1 * lspec1;
  const double lspec2   = log(1 + sqrtxz1 + z);
  const double lspec3   = log(sqrtxz3);
  const double lspec4   = log(sqrtxz3*z);
  const double li2spec1  = apfel::dilog(0.5 - sqrtxz1/2. - z/2.);
  const double li2spec2  = apfel::dilog(0.5 - sqrtxz1/2. + z/2.);
  const double li2spec3  = apfel::dilog(0.5 - zi/2. - (sqrtxz1*zi)/2.);
  const double li2spec4  = apfel::dilog(0.5 + zi/2. - (sqrtxz1*zi)/2.);
  const double li2spec19 = apfel::dilog(-(x*z));
  const double li2spec20 = apfel::dilog(-(x*zi));
  const double atanspec1 = atan(sqrtxz3);
  const double atanspec2 = atan(sqrtxz3*z);
  const double itani1 = InvTanInt(-sqrtxz3);
  const double itani2 = InvTanInt(sqrtxz3);
  const double itani3 = InvTanInt(sqrtxz3*z);
  return (-5*NC)/4. + (5*NC*x)/4. + (21*NC*z)/8. + (43*NC*x*z)/8. +  
   (7*NC*sqrtxz3*itani1)/4. -  
   (105*NC*sqrtxz3*x*z*itani1)/8. -  
   (7*NC*sqrtxz3*itani2)/4. +  
   (105*NC*sqrtxz3*x*z*itani2)/8. +  
   (7*NC*sqrtxz3*itani3)/2. -  
   (105*NC*sqrtxz3*x*z*itani3)/4. - 4*NC*Li2mx +  
   4*NC*x*Li2mx + 4*NC*z*Li2mx - 4*NC*x*z*Li2mx + 4*NC*Li2x -  
   4*NC*z*Li2x + 4*NC*x*z*li2spec1 -  
   4*NC*x*z*li2spec2 - 4*NC*x*Li2z + 2*NC*li2spec19 -  
   4*NC*x*li2spec19 - 2*NC*z*li2spec19 + 2*NC*li2spec20 -  
   2*NC*z*li2spec20 + 4*NC*x*z*li2spec20 -  
   4*NC*x*li2spec3 +  
   4*NC*x*li2spec4 +  
   (7*NC*sqrtxz3*atanspec1*lspec3)/2. -  
   (105*NC*sqrtxz3*x*z*atanspec1*lspec3)/4. + NC*lomx -  
   NC*x*lomx - NC*z*lomx - 3*NC*x*z*lomx -  
   (79*NC*lx)/8. - (7*NC*x*lx)/8. + (157*NC*z*lx)/16. +  
   (109*NC*x*z*lx)/16. + 4*NC*x*l2*lx - 8*NC*x*z*l2*lx +  
   4*NC*lomx*lx + 2*NC*x*lomx*lx -  
   4*NC*z*lomx*lx - 2*NC*x*z*lomx*lx -  
   4*NC*lx*lopx + 4*NC*x*lx*lopx +  
   4*NC*z*lx*lopx - 4*NC*x*z*lx*lopx + NC*lomz -  
   NC*x*lomz - NC*z*lomz - 3*NC*x*z*lomz +  
   2*NC*x*lx*lomz - 2*NC*x*z*lx*lomz -  
   4*NC*x*l2*lspec1 +  
   12*NC*x*z*l2*lspec1 -  
   4*NC*x*lx*lspec1 +  
   4*NC*x*z*lx*lspec1 - (17*NC*lz)/8. +  
   (NC*x*lz)/8. - (3*NC*z*lz)/16. + (3*NC*x*z*lz)/16. -  
   4*NC*x*l2*lz - 8*NC*x*z*l2*lz + 4*NC*x*lomx*lz -  
   2*NC*lx*lz - 8*NC*x*lx*lz + 2*NC*z*lx*lz -  
   4*NC*x*z*lx*lz + 4*NC*x*z*lspec1*lz -  
   (7*NC*sqrtxz3*atanspec2*lspec4)/2. +  
   (105*NC*sqrtxz3*x*z*atanspec2*lspec4)/4. +  
   4*NC*x*l2*lspec2 +  
   4*NC*x*z*l2*lspec2 +  
   4*NC*x*z*lx*lspec2 -  
   4*NC*x*lspec1*lspec2 -  
   4*NC*x*z*lspec1*lspec2 +  
   4*NC*x*lz*lspec2 +  
   4*NC*x*z*lz*lspec2 + 2*NC*x*lx*lxpz +  
   2*NC*x*z*lx*lxpz - 2*NC*x*lz*lxpz -  
   2*NC*x*z*lz*lxpz + 2*NC*lx*lopxz -  
   4*NC*x*lx*lopxz - 2*NC*z*lx*lopxz +  
   2*NC*lz*lopxz - 4*NC*x*lz*lopxz -  
   2*NC*z*lz*lopxz + 2*NC*lx*lopxzi -  
   2*NC*x*lx*lopxzi - 2*NC*z*lx*lopxzi +  
   2*NC*x*z*lx*lopxzi - 2*NC*lz*lopxzi +  
   2*NC*x*lz*lopxzi + 2*NC*z*lz*lopxzi -  
   2*NC*x*z*lz*lopxzi + 2*NCi + 14*x*NCi -  
   2*z*NCi - 6*x*z*NCi +  
   (sqrtxz3*itani1*NCi)/2. +  
   (15*sqrtxz3*x*z*itani1*NCi)/2. -  
   (sqrtxz3*itani2*NCi)/2. -  
   (15*sqrtxz3*x*z*itani2*NCi)/2. +  
   sqrtxz3*itani3*NCi +  
   15*sqrtxz3*x*z*itani3*NCi - 4*z*Li2mx*NCi +  
   4*x*z*Li2mx*NCi + 4*z*Li2x*NCi +  
   4*x*li2spec1*NCi -  
   4*x*z*li2spec1*NCi -  
   4*x*li2spec2*NCi +  
   4*x*z*li2spec2*NCi + 4*x*Li2z*NCi +  
   2*z*li2spec19*NCi + 2*z*li2spec20*NCi -  
   4*x*z*li2spec20*NCi +  
   sqrtxz3*atanspec1*lspec3*NCi +  
   15*sqrtxz3*x*z*atanspec1*lspec3*NCi -  
   lomx*NCi + x*lomx*NCi +  
   z*lomx*NCi + 3*x*z*lomx*NCi +  
   (3*lx*NCi)/2. + (5*x*lx*NCi)/2. -  
   (19*z*lx*NCi)/2. - (13*x*z*lx*NCi)/2. -  
   8*x*l2*lx*NCi + 8*x*z*l2*lx*NCi -  
   2*x*lomx*lx*NCi + 4*z*lomx*lx*NCi +  
   2*x*z*lomx*lx*NCi - 4*z*lx*lopx*NCi +  
   4*x*z*lx*lopx*NCi - lomz*NCi +  
   x*lomz*NCi + z*lomz*NCi +  
   3*x*z*lomz*NCi - 2*x*lx*lomz*NCi +  
   2*x*z*lx*lomz*NCi +  
   12*x*l2*lspec1*NCi -  
   12*x*z*l2*lspec1*NCi +  
   4*x*lx*lspec1*NCi -  
   4*x*z*lx*lspec1*NCi + (lz*NCi)/2. -  
   (9*x*lz*NCi)/2. + (z*lz*NCi)/2. -  
   (x*z*lz*NCi)/2. - 8*x*l2*lz*NCi +  
   8*x*z*l2*lz*NCi - 4*x*lomx*lz*NCi +  
   2*x*lx*lz*NCi - 2*z*lx*lz*NCi +  
   4*x*z*lx*lz*NCi +  
   4*x*lspec1*lz*NCi -  
   4*x*z*lspec1*lz*NCi -  
   sqrtxz3*atanspec2*lspec4*NCi -  
   15*sqrtxz3*x*z*atanspec2*lspec4*NCi +  
   4*x*l2*lspec2*NCi -  
   4*x*z*l2*lspec2*NCi +  
   4*x*lx*lspec2*NCi -  
   4*x*z*lx*lspec2*NCi -  
   4*x*lspec1*lspec2*NCi +  
   4*x*z*lspec1*lspec2*NCi +  
   4*x*lz*lspec2*NCi -  
   4*x*z*lz*lspec2*NCi -  
   2*x*z*lx*lxpz*NCi + 2*x*z*lz*lxpz*NCi +  
   2*z*lx*lopxz*NCi + 2*z*lz*lopxz*NCi +  
   2*z*lx*lopxzi*NCi -  
   2*x*z*lx*lopxzi*NCi -  
   2*z*lz*lopxzi*NCi +  
   2*x*z*lz*lopxzi*NCi - (2*NC*pi2)/3. +  
   (NC*x*pi2)/3. + (2*NC*z*pi2)/3. + (NC*x*z*pi2)/3. -  
   (x*NCi*pi2)/3. - (2*z*NCi*pi2)/3. -  
   (x*z*NCi*pi2)/3. +  
   (3*NC*sqrtxz3*itani1*xi2)/16. -  
   (3*NC*sqrtxz3*itani2*xi2)/16. +  
   (3*NC*sqrtxz3*itani3*xi2)/8. +  
   (3*NC*sqrtxz3*atanspec1*lspec3*xi2)/8. -  
   (3*NC*sqrtxz3*atanspec2*lspec4*xi2)/8. -  
   (3*NC*xi)/8. + (7*NC*sqrtxz3*z*itani1*xi)/8. -  
   (7*NC*sqrtxz3*z*itani2*xi)/8. +  
   (7*NC*sqrtxz3*z*itani3*xi)/4. -  
   4*NC*Li2mx*xi + 4*NC*z*Li2mx*xi +  
   4*NC*Li2x*xi - 4*NC*z*Li2x*xi +  
   2*NC*li2spec19*xi - 2*NC*z*li2spec19*xi +  
   2*NC*li2spec20*xi -  
   2*NC*z*li2spec20*xi +  
   (7*NC*sqrtxz3*z*atanspec1*lspec3*xi)/4. -  
   (125*NC*lx*xi)/16. + 8*NC*z*lx*xi +  
   4*NC*lomx*lx*xi - 4*NC*z*lomx*lx*xi -  
   4*NC*lx*lopx*xi + 4*NC*z*lx*lopx*xi +  
   (3*NC*lz*xi)/16. - 2*NC*lx*lz*xi +  
   2*NC*z*lx*lz*xi -  
   (7*NC*sqrtxz3*z*atanspec2*lspec4*xi)/4. +  
   2*NC*lx*lopxz*xi -  
   2*NC*z*lx*lopxz*xi +  
   2*NC*lz*lopxz*xi -  
   2*NC*z*lz*lopxz*xi +  
   2*NC*lx*lopxzi*xi -  
   2*NC*z*lx*lopxzi*xi -  
   2*NC*lz*lopxzi*xi +  
   2*NC*z*lz*lopxzi*xi -  
   (sqrtxz3*z*itani1*NCi*xi)/2. +  
   (sqrtxz3*z*itani2*NCi*xi)/2. -  
   sqrtxz3*z*itani3*NCi*xi -  
   4*z*Li2mx*NCi*xi + 4*z*Li2x*NCi*xi +  
   2*z*li2spec19*NCi*xi +  
   2*z*li2spec20*NCi*xi -  
   sqrtxz3*z*atanspec1*lspec3*NCi*xi -  
   8*z*lx*NCi*xi +  
   4*z*lomx*lx*NCi*xi -  
   4*z*lx*lopx*NCi*xi -  
   2*z*lx*lz*NCi*xi +  
   sqrtxz3*z*atanspec2*lspec4*NCi*xi +  
   2*z*lx*lopxz*NCi*xi +  
   2*z*lz*lopxz*NCi*xi +  
   2*z*lx*lopxzi*NCi*xi -  
   2*z*lz*lopxzi*NCi*xi -  
   (2*NC*pi2*xi)/3. + (2*NC*z*pi2*xi)/3. -  
   (2*z*NCi*pi2*xi)/3. + (3*NC*x2)/8. -  
   8*NC*z*x2 + (35*NC*sqrtxz3*itani1*x2)/16. -  
   (35*NC*sqrtxz3*itani2*x2)/16. +  
   (35*NC*sqrtxz3*itani3*x2)/8. + 4*NC*Li2mx*x2 -  
   4*NC*z*Li2mx*x2 - 2*NC*Li2x*x2 + 2*NC*z*Li2x*x2 -  
   4*NC*z*li2spec1*x2 +  
   4*NC*z*li2spec2*x2 + 4*NC*Li2z*x2 +  
   4*NC*z*li2spec19*x2 - 4*NC*li2spec20*x2 +  
   4*NC*li2spec3*x2 -  
   4*NC*li2spec4*x2 +  
   (35*NC*sqrtxz3*atanspec1*lspec3*x2)/8. +  
   4*NC*z*lomx*x2 + (35*NC*lx*x2)/16. -  
   4*NC*z*lx*x2 - 4*NC*l2*lx*x2 +  
   8*NC*z*l2*lx*x2 - 2*NC*lomx*lx*x2 +  
   2*NC*z*lomx*lx*x2 + 4*NC*lx*lopx*x2 -  
   4*NC*z*lx*lopx*x2 + 4*NC*z*lomz*x2 +  
   4*NC*l2*lspec1*x2 -  
   12*NC*z*l2*lspec1*x2 +  
   4*NC*lx*lspec1*x2 -  
   4*NC*z*lx*lspec1*x2 +  
   (29*NC*lz*x2)/16. + 4*NC*l2*lz*x2 +  
   8*NC*z*l2*lz*x2 - 4*NC*lomx*lz*x2 +  
   2*NC*lx*lz*x2 + 2*NC*z*lx*lz*x2 -  
   4*NC*z*lspec1*lz*x2 -  
   (35*NC*sqrtxz3*atanspec2*lspec4*x2)/8. -  
   4*NC*l2*lspec2*x2 -  
   4*NC*z*l2*lspec2*x2 -  
   4*NC*z*lx*lspec2*x2 +  
   4*NC*lspec1*lspec2*x2 +  
   4*NC*z*lspec1*lspec2*x2 -  
   4*NC*lz*lspec2*x2 -  
   4*NC*z*lz*lspec2*x2 -  
   2*NC*lx*lxpz*x2 - 2*NC*z*lx*lxpz*x2 +  
   2*NC*lz*lxpz*x2 + 2*NC*z*lz*lxpz*x2 +  
   4*NC*z*lx*lopxz*x2 +  
   4*NC*z*lz*lopxz*x2 -  
   2*NC*lx*lopxzi*x2 +  
   2*NC*z*lx*lopxzi*x2 +  
   2*NC*lz*lopxzi*x2 -  
   2*NC*z*lz*lopxzi*x2 - 16*NCi*x2 +  
   8*z*NCi*x2 + 4*z*Li2mx*NCi*x2 +  
   2*Li2x*NCi*x2 - 2*z*Li2x*NCi*x2 -  
   4*li2spec1*NCi*x2 +  
   4*z*li2spec1*NCi*x2 +  
   4*li2spec2*NCi*x2 -  
   4*z*li2spec2*NCi*x2 -  
   4*Li2z*NCi*x2 - 4*z*li2spec19*NCi*x2 -  
   4*z*lomx*NCi*x2 + 4*z*lx*NCi*x2 +  
   8*l2*lx*NCi*x2 -  
   8*z*l2*lx*NCi*x2 +  
   2*lomx*lx*NCi*x2 -  
   2*z*lomx*lx*NCi*x2 +  
   4*z*lx*lopx*NCi*x2 -  
   4*z*lomz*NCi*x2 -  
   12*l2*lspec1*NCi*x2 +  
   12*z*l2*lspec1*NCi*x2 -  
   4*lx*lspec1*NCi*x2 +  
   4*z*lx*lspec1*NCi*x2 +  
   4*lz*NCi*x2 + 8*l2*lz*NCi*x2 -  
   8*z*l2*lz*NCi*x2 +  
   4*lomx*lz*NCi*x2 -  
   2*z*lx*lz*NCi*x2 -  
   4*lspec1*lz*NCi*x2 +  
   4*z*lspec1*lz*NCi*x2 -  
   4*l2*lspec2*NCi*x2 +  
   4*z*l2*lspec2*NCi*x2 -  
   4*lx*lspec2*NCi*x2 +  
   4*z*lx*lspec2*NCi*x2 +  
   4*lspec1*lspec2*NCi*x2 -  
   4*z*lspec1*lspec2*NCi*x2 -  
   4*lz*lspec2*NCi*x2 +  
   4*z*lz*lspec2*NCi*x2 +  
   2*z*lx*lxpz*NCi*x2 -  
   2*z*lz*lxpz*NCi*x2 -  
   4*z*lx*lopxz*NCi*x2 -  
   4*z*lz*lopxz*NCi*x2 -  
   2*z*lx*lopxzi*NCi*x2 +  
   2*z*lz*lopxzi*NCi*x2 -  
   (2*NC*z*pi2*x2)/3. + (2*z*NCi*pi2*x2)/3. +  
   8*NC*Li2mx*opxi + 4*NC*x*Li2mx*opxi -  
   8*NC*z*Li2mx*opxi - 4*NC*x*z*Li2mx*opxi -  
   8*NC*Li2x*opxi - 4*NC*x*Li2x*opxi +  
   8*NC*z*Li2x*opxi + 4*NC*x*z*Li2x*opxi -  
   4*NC*li2spec19*opxi - 2*NC*x*li2spec19*opxi +  
   4*NC*z*li2spec19*opxi + 2*NC*x*z*li2spec19*opxi -  
   4*NC*li2spec20*opxi -  
   2*NC*x*li2spec20*opxi +  
   4*NC*z*li2spec20*opxi +  
   2*NC*x*z*li2spec20*opxi + 16*NC*lx*opxi +  
   8*NC*x*lx*opxi - 16*NC*z*lx*opxi -  
   8*NC*x*z*lx*opxi - 8*NC*lomx*lx*opxi -  
   4*NC*x*lomx*lx*opxi +  
   8*NC*z*lomx*lx*opxi +  
   4*NC*x*z*lomx*lx*opxi +  
   8*NC*lx*lopx*opxi +  
   4*NC*x*lx*lopx*opxi -  
   8*NC*z*lx*lopx*opxi -  
   4*NC*x*z*lx*lopx*opxi +  
   4*NC*lx*lz*opxi + 2*NC*x*lx*lz*opxi -  
   4*NC*z*lx*lz*opxi -  
   2*NC*x*z*lx*lz*opxi -  
   4*NC*lx*lopxz*opxi -  
   2*NC*x*lx*lopxz*opxi +  
   4*NC*z*lx*lopxz*opxi +  
   2*NC*x*z*lx*lopxz*opxi -  
   4*NC*lz*lopxz*opxi -  
   2*NC*x*lz*lopxz*opxi +  
   4*NC*z*lz*lopxz*opxi +  
   2*NC*x*z*lz*lopxz*opxi -  
   4*NC*lx*lopxzi*opxi -  
   2*NC*x*lx*lopxzi*opxi +  
   4*NC*z*lx*lopxzi*opxi +  
   2*NC*x*z*lx*lopxzi*opxi +  
   4*NC*lz*lopxzi*opxi +  
   2*NC*x*lz*lopxzi*opxi -  
   4*NC*z*lz*lopxzi*opxi -  
   2*NC*x*z*lz*lopxzi*opxi +  
   8*z*Li2mx*NCi*opxi +  
   4*x*z*Li2mx*NCi*opxi -  
   8*z*Li2x*NCi*opxi -  
   4*x*z*Li2x*NCi*opxi -  
   4*z*li2spec19*NCi*opxi -  
   2*x*z*li2spec19*NCi*opxi -  
   4*z*li2spec20*NCi*opxi -  
   2*x*z*li2spec20*NCi*opxi +  
   16*z*lx*NCi*opxi +  
   8*x*z*lx*NCi*opxi -  
   8*z*lomx*lx*NCi*opxi -  
   4*x*z*lomx*lx*NCi*opxi +  
   8*z*lx*lopx*NCi*opxi +  
   4*x*z*lx*lopx*NCi*opxi +  
   4*z*lx*lz*NCi*opxi +  
   2*x*z*lx*lz*NCi*opxi -  
   4*z*lx*lopxz*NCi*opxi -  
   2*x*z*lx*lopxz*NCi*opxi -  
   4*z*lz*lopxz*NCi*opxi -  
   2*x*z*lz*lopxz*NCi*opxi -  
   4*z*lx*lopxzi*NCi*opxi -  
   2*x*z*lx*lopxzi*NCi*opxi +  
   4*z*lz*lopxzi*NCi*opxi +  
   2*x*z*lz*lopxzi*NCi*opxi +  
   (4*NC*pi2*opxi)/3. + (2*NC*x*pi2*opxi)/3. -  
   (4*NC*z*pi2*opxi)/3. -  
   (2*NC*x*z*pi2*opxi)/3. +  
   (4*z*NCi*pi2*opxi)/3. +  
   (2*x*z*NCi*pi2*opxi)/3. +  
   4*NC*Li2mx*xi*opxi -  
   4*NC*z*Li2mx*xi*opxi -  
   4*NC*Li2x*xi*opxi +  
   4*NC*z*Li2x*xi*opxi -  
   2*NC*li2spec19*xi*opxi +  
   2*NC*z*li2spec19*xi*opxi -  
   2*NC*li2spec20*xi*opxi +  
   2*NC*z*li2spec20*xi*opxi +  
   8*NC*lx*xi*opxi -  
   8*NC*z*lx*xi*opxi -  
   4*NC*lomx*lx*xi*opxi +  
   4*NC*z*lomx*lx*xi*opxi +  
   4*NC*lx*lopx*xi*opxi -  
   4*NC*z*lx*lopx*xi*opxi +  
   2*NC*lx*lz*xi*opxi -  
   2*NC*z*lx*lz*xi*opxi -  
   2*NC*lx*lopxz*xi*opxi +  
   2*NC*z*lx*lopxz*xi*opxi -  
   2*NC*lz*lopxz*xi*opxi +  
   2*NC*z*lz*lopxz*xi*opxi -  
   2*NC*lx*lopxzi*xi*opxi +  
   2*NC*z*lx*lopxzi*xi*opxi +  
   2*NC*lz*lopxzi*xi*opxi -  
   2*NC*z*lz*lopxzi*xi*opxi +  
   4*z*Li2mx*NCi*xi*opxi -  
   4*z*Li2x*NCi*xi*opxi -  
   2*z*li2spec19*NCi*xi*opxi -  
   2*z*li2spec20*NCi*xi*opxi +  
   8*z*lx*NCi*xi*opxi -  
   4*z*lomx*lx*NCi*xi*opxi +  
   4*z*lx*lopx*NCi*xi*opxi +  
   2*z*lx*lz*NCi*xi*opxi -  
   2*z*lx*lopxz*NCi*xi*opxi -  
   2*z*lz*lopxz*NCi*xi*opxi -  
   2*z*lx*lopxzi*NCi*xi*opxi +  
   2*z*lz*lopxzi*NCi*xi*opxi +  
   (2*NC*pi2*xi*opxi)/3. -  
   (2*NC*z*pi2*xi*opxi)/3. +  
   (2*z*NCi*pi2*xi*opxi)/3. +  
   (3*NC*zi2)/8. - (3*NC*x*zi2)/8. +  
   (3*NC*sqrtxz3*itani1*zi2)/16. -  
   (3*NC*sqrtxz3*itani2*zi2)/16. +  
   (3*NC*sqrtxz3*itani3*zi2)/8. +  
   (3*NC*sqrtxz3*atanspec1*lspec3*zi2)/8. +  
   (3*NC*lx*zi2)/16. + (3*NC*x*lx*zi2)/16. -  
   (3*NC*lz*zi2)/16. + (3*NC*x*lz*zi2)/16. -  
   (3*NC*sqrtxz3*atanspec2*lspec4*zi2)/8. -  
   (7*NC*zi)/4. - (25*NC*x*zi)/4. +  
   (15*NC*sqrtxz3*x*itani1*zi)/8. -  
   (15*NC*sqrtxz3*x*itani2*zi)/8. +  
   (15*NC*sqrtxz3*x*itani3*zi)/4. -  
   8*NC*sqrtxz1*x*l2*zi +  
   (15*NC*sqrtxz3*x*atanspec1*lspec3*zi)/4. +  
   4*NC*x*lomx*zi - (NC*lx*zi)/8. -  
   (49*NC*x*lx*zi)/8. - 4*NC*sqrtxz1*x*lx*zi +  
   4*NC*x*lomz*zi +  
   8*NC*sqrtxz1*x*lspec1*zi - (NC*lz*zi)/8. +  
   (49*NC*x*lz*zi)/8. - 4*NC*sqrtxz1*x*lz*zi -  
   (15*NC*sqrtxz3*x*atanspec2*lspec4*zi)/4. -  
   8*x*NCi*zi - 4*x*lomx*NCi*zi +  
   4*x*lx*NCi*zi - 4*x*lomz*NCi*zi -  
   8*x*lz*NCi*zi + (3*NC*xi*zi)/8. -  
   (NC*sqrtxz3*itani1*xi*zi)/8. +  
   (NC*sqrtxz3*itani2*xi*zi)/8. -  
   (NC*sqrtxz3*itani3*xi*zi)/4. -  
   (NC*sqrtxz3*atanspec1*lspec3*xi*zi)/4. -  
   (3*NC*lx*xi*zi)/16. +  
   (3*NC*lz*xi*zi)/16. +  
   (NC*sqrtxz3*atanspec2*lspec4*xi*zi)/4. +  
   (61*NC*x2*zi)/8. + 8*NC*sqrtxz1*l2*x2*zi -  
   4*NC*lomx*x2*zi +  
   (29*NC*lx*x2*zi)/16. +  
   4*NC*sqrtxz1*lx*x2*zi -  
   4*NC*lomz*x2*zi -  
   8*NC*sqrtxz1*lspec1*x2*zi -  
   (99*NC*lz*x2*zi)/16. +  
   4*NC*sqrtxz1*lz*x2*zi + 8*NCi*x2*zi +  
   4*lomx*NCi*x2*zi -  
   4*lx*NCi*x2*zi +  
   4*lomz*NCi*x2*zi +  
   8*lz*NCi*x2*zi +  
   (19*NC*sqrtxz3*itani1*z2)/16. -  
   (19*NC*sqrtxz3*itani2*z2)/16. +  
   (19*NC*sqrtxz3*itani3*z2)/8. +  
   (19*NC*sqrtxz3*atanspec1*lspec3*z2)/8. -  
   (19*NC*sqrtxz3*atanspec2*lspec4*z2)/8. -  
   (3*sqrtxz3*itani1*NCi*z2)/2. +  
   (3*sqrtxz3*itani2*NCi*z2)/2. -  
   3*sqrtxz3*itani3*NCi*z2 -  
   3*sqrtxz3*atanspec1*lspec3*NCi*z2 +  
   3*sqrtxz3*atanspec2*lspec4*NCi*z2 -  
   8*NC*x*z*l22 - 8*x*NCi*l22 +  
   8*x*z*NCi*l22 + 8*NC*z*x2*l22 +  
   8*NCi*x2*l22 -  
   8*z*NCi*x2*l22 + 3*NC*lx2 -  
   NC*x*lx2 - 3*NC*z*lx2 + NC*x*z*lx2 +  
   x*NCi*lx2 + 3*z*NCi*lx2 -  
   x*z*NCi*lx2 + 3*NC*xi*lx2 -  
   3*NC*z*xi*lx2 + 3*z*NCi*xi*lx2 -  
   NC*x2*lx2 + NC*z*x2*lx2 +  
   NCi*x2*lx2 - z*NCi*x2*lx2 -  
   6*NC*opxi*lx2 - 3*NC*x*opxi*lx2 +  
   6*NC*z*opxi*lx2 +  
   3*NC*x*z*opxi*lx2 -  
   6*z*NCi*opxi*lx2 -  
   3*x*z*NCi*opxi*lx2 -  
   3*NC*xi*opxi*lx2 +  
   3*NC*z*xi*opxi*lx2 -  
   3*z*NCi*xi*opxi*lx2 +  
   4*NC*x*lspec1_2 -  
   4*NC*x*z*lspec1_2 -  
   4*x*NCi*lspec1_2 +  
   4*x*z*NCi*lspec1_2 -  
   4*NC*x2*lspec1_2 +  
   4*NC*z*x2*lspec1_2 +  
   4*NCi*x2*lspec1_2 -  
   4*z*NCi*x2*lspec1_2 - NC*lz2 +  
   4*NC*x*lz2 + NC*z*lz2 - 4*x*NCi*lz2 -  
   z*NCi*lz2 - NC*xi*lz2 +  
   NC*z*xi*lz2 - z*NCi*xi*lz2 -  
   2*NC*x2*lz2 - 2*NC*z*x2*lz2 +  
   4*NCi*x2*lz2 +  
   2*z*NCi*x2*lz2 + 2*NC*opxi*lz2 +  
   NC*x*opxi*lz2 - 2*NC*z*opxi*lz2 -  
   NC*x*z*opxi*lz2 +  
   2*z*NCi*opxi*lz2 +  
   x*z*NCi*opxi*lz2 +  
   NC*xi*opxi*lz2 -  
   NC*z*xi*opxi*lz2 +  
   z*NCi*xi*opxi*lz2;
};
double C2LG2Q_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double poly2    = 1 + 2*x + x*x - 4*x*z;
  const double poly2i   = 1. / poly2;
  const double poly2i2  = poly2i * poly2i;
  const double sqrtxz2  = sqrt(poly2);
  const double sqrtxz2i = 1. / sqrtxz2;
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double x2 = x * x;
  const double x3 = x * x2;
  const double x4 = x * x3;
  const double x5 = x * x4;
  const double x6 = x * x5;
  const double x7 = x * x6;
  const double xi  = 1. / x;
  const double z2 = z * z;
  const double z3 = z * z2;
  const double zi  = 1. / z;
  const double omxi = 1. / ( 1 - x );
  const double omzi  = 1. / ( 1 - z );
  const double lx  = log(x);
  const double lx2 = lx * lx;
  const double lz  = log(z);
  const double lz2 = lz * lz;
  const double lomx  = log(1 - x);
  const double lomx2 = lomx * lomx;
  const double lomz  = log(1 - z);
  const double lomz2 = lomz * lomz;
  const double lopx  = log(1 + x);
  const double Li2x  = apfel::dilog(x);
  const double Li2mx = apfel::dilog(-x);
  const double Li2z  = apfel::dilog(z);
  const double omxmzi  = 1. / ( 1 - x - z );
  const double omxmzi2 = omxmzi * omxmzi;
  const double xmzi  = 1. / ( x - z );
  const double xmzi2 = xmzi * xmzi;
  const double lomxmz  = log(1 - x - z);
  const double lmopxpz = log(-1 + x + z);
  const double lxmz    = log(x - z);
  const double lmxpz   = log(-x + z);
  const double li2omxzi = apfel::dilog(1 - x*zi);
  const double lspec5   = log(1 - sqrtxz2 + x);
  const double lspec6   = log(1 + sqrtxz2 + x);
  const double li2spec5  = apfel::dilog(0.5 - sqrtxz2/2. - x/2.);
  const double li2spec6  = apfel::dilog(0.5 + sqrtxz2/2. - x/2.);
  const double li2spec7  = apfel::dilog(0.5 - xi/2. - (sqrtxz2*xi)/2.);
  const double li2spec8  = apfel::dilog(0.5 - xi/2. + (sqrtxz2*xi)/2.);
  const double li2spec9  = apfel::dilog((1 - z)*omxi);
  const double li2spec10 = apfel::dilog(x*(1 - z)*omxi*zi);
  const double li2spec11 = apfel::dilog((1 - x)*omzi);
  const double li2spec12 = apfel::dilog((1 - x)*z*xi*omzi);
  const double li2spec13 = apfel::dilog(z*omxi);
  const double li2spec14 = apfel::dilog((1 - z)*z*omxi*xi);
  const double li2spec15 = apfel::dilog(x*z*omxi*omzi);
  const double li2spec16 = apfel::dilog((1 - x)*zi);
  const double li2spec17 = apfel::dilog((1 - x)*(1 - z)*xi*zi);
  const double li2spec18 = apfel::dilog((1 - x)*x*omzi*zi);
  const double Tu1 = (z < 1 - x && z < x ? 1 : 0);
  const double Tu2 = (z > 1 - x && z < x ? 1 : 0);
  const double Tu3 = (z < 1 - x && z > x ? 1 : 0);
  const double Tu4 = (z > 1 - x && z > x ? 1 : 0);
  return (-409*NC)/24. - (95*NC*x)/24. + (41*NC*z)/2. - (57*NC*x*z)/2. +  
   4*NC*x*z*Li2mx + 18*NC*x*Li2x - 4*NC*x*z*Li2x + 8*NC*x*Li2z -  
   6*NC*x*li2omxzi - 4*NC*lomx - 29*NC*x*lomx +  
   NC*z*lomx + 3*NC*x*z*lomx - 2*NC*lx +  
   (19*NC*x*lx)/2. + (39*NC*z*lx)/4. + (23*NC*x*z*lx)/4. -  
   18*NC*x*lomx*lx - 2*NC*x*z*lomx*lx +  
   4*NC*x*z*lx*lopx - 4*NC*lomz - 30*NC*x*lomz +  
   NC*z*lomz + 3*NC*x*z*lomz + 16*NC*x*lomx*lomz -  
   40*NC*x*lx*lomz + 2*NC*x*z*lx*lomz -  
   (27*NC*lz)/4. - (95*NC*x*lz)/4. + (47*NC*z*lz)/4. -  
   (51*NC*x*z*lz)/4. + 10*NC*x*lomx*lz -  
   22*NC*x*lx*lz - 2*NC*x*z*lx*lz +  
   22*NC*x*lomz*lz + (3*NCi)/8. + (37*x*NCi)/8. -  
   (z*NCi)/2. + (17*x*z*NCi)/2. + 4*x*Li2mx*NCi -  
   4*x*z*Li2mx*NCi + 2*x*Li2x*NCi +  
   4*x*z*Li2x*NCi - 8*x*Li2z*NCi +  
   2*x*li2omxzi*NCi + x*lomx*NCi -  
   z*lomx*NCi - 3*x*z*lomx*NCi +  
   (lx*NCi)/2. + x*lx*NCi + (z*lx*NCi)/4. +  
   (17*x*z*lx*NCi)/4. + 22*x*lomx*lx*NCi +  
   2*x*z*lomx*lx*NCi + 4*x*lx*lopx*NCi -  
   4*x*z*lx*lopx*NCi + 4*x*lomz*NCi -  
   z*lomz*NCi - 3*x*z*lomz*NCi -  
   14*x*lomx*lomz*NCi +  
   20*x*lx*lomz*NCi - 2*x*z*lx*lomz*NCi +  
   (lz*NCi)/4. - (15*x*lz*NCi)/4. -  
   (7*z*lz*NCi)/4. + (11*x*z*lz*NCi)/4. -  
   8*x*lomx*lz*NCi + 10*x*lx*lz*NCi +  
   2*x*z*lx*lz*NCi - 14*x*lomz*lz*NCi -  
   (17*NC*x*pi2)/3. + NC*x*z*pi2 + 3*x*NCi*pi2 -  
   x*z*NCi*pi2 + (3*NC*lx*poly2i2)/16. +  
   (9*NC*x*lx*poly2i2)/16. - (3*NC*lz*poly2i2)/16. +  
   (9*NC*x*lz*poly2i2)/16. -  
   (3*lx*NCi*poly2i2)/16. -  
   (9*x*lx*NCi*poly2i2)/16. +  
   (3*lz*NCi*poly2i2)/16. -  
   (9*x*lz*NCi*poly2i2)/16. + (NC*poly2i)/8. +  
   (NC*x*poly2i)/4. + (7*NC*lx*poly2i)/16. +  
   (3*NC*x*lx*poly2i)/8. - (7*NC*lz*poly2i)/16. +  
   (3*NC*x*lz*poly2i)/8. - (NCi*poly2i)/8. -  
   (x*NCi*poly2i)/4. + (lx*NCi*poly2i)/16. -  
   (3*x*lx*NCi*poly2i)/8. -  
   (lz*NCi*poly2i)/16. -  
   (3*x*lz*NCi*poly2i)/8. +  
   (5*NC*li2spec5*sqrtxz2i)/8. +  
   (151*NC*x*li2spec5*sqrtxz2i)/16. -  
   (5*NC*z*li2spec5*sqrtxz2i)/4. -  
   (51*NC*x*z*li2spec5*sqrtxz2i)/2. -  
   (5*NC*li2spec6*sqrtxz2i)/8. -  
   (151*NC*x*li2spec6*sqrtxz2i)/16. +  
   (5*NC*z*li2spec6*sqrtxz2i)/4. +  
   (51*NC*x*z*li2spec6*sqrtxz2i)/2. -  
   (5*NC*li2spec7*sqrtxz2i)/ 
    8. - (151*NC*x*li2spec7* 
      sqrtxz2i)/16. + (5*NC*z* 
      li2spec7*sqrtxz2i)/4. +  
   (51*NC*x*z*li2spec7* 
      sqrtxz2i)/2. + (5*NC* 
      li2spec8*sqrtxz2i)/8. +  
   (151*NC*x*li2spec8* 
      sqrtxz2i)/16. - (5*NC*z* 
      li2spec8*sqrtxz2i)/4. -  
   (51*NC*x*z*li2spec8* 
      sqrtxz2i)/2. - (5*NC*lx*lspec5* 
      sqrtxz2i)/8. - (151*NC*x*lx*lspec5* 
      sqrtxz2i)/16. + (5*NC*z*lx*lspec5* 
      sqrtxz2i)/4. + (51*NC*x*z*lx*lspec5* 
      sqrtxz2i)/2. + (5*NC*lx*lspec6* 
      sqrtxz2i)/8. + (151*NC*x*lx*lspec6* 
      sqrtxz2i)/16. - (5*NC*z*lx*lspec6* 
      sqrtxz2i)/4. - (51*NC*x*z*lx*lspec6* 
      sqrtxz2i)/2. - (li2spec5*NCi* 
      sqrtxz2i)/8. - (31*x*li2spec5*NCi* 
      sqrtxz2i)/16. + (z*li2spec5*NCi* 
      sqrtxz2i)/4. + (11*x*z*li2spec5*NCi* 
      sqrtxz2i)/2. + (li2spec6*NCi* 
      sqrtxz2i)/8. + (31*x*li2spec6*NCi* 
      sqrtxz2i)/16. - (z*li2spec6*NCi* 
      sqrtxz2i)/4. - (11*x*z*li2spec6*NCi* 
      sqrtxz2i)/2. + (li2spec7* 
      NCi*sqrtxz2i)/8. +  
   (31*x*li2spec7*NCi* 
      sqrtxz2i)/16. - (z* 
      li2spec7*NCi* 
      sqrtxz2i)/4. - (11*x*z* 
      li2spec7*NCi* 
      sqrtxz2i)/2. - (li2spec8* 
      NCi*sqrtxz2i)/8. -  
   (31*x*li2spec8*NCi* 
      sqrtxz2i)/16. + (z* 
      li2spec8*NCi* 
      sqrtxz2i)/4. + (11*x*z* 
      li2spec8*NCi* 
      sqrtxz2i)/2. + (lx*lspec5*NCi* 
      sqrtxz2i)/8. + (31*x*lx*lspec5*NCi* 
      sqrtxz2i)/16. - (z*lx*lspec5*NCi* 
      sqrtxz2i)/4. - (11*x*z*lx*lspec5*NCi* 
      sqrtxz2i)/2. - (lx*lspec6*NCi* 
      sqrtxz2i)/8. - (31*x*lx*lspec6*NCi* 
      sqrtxz2i)/16. + (z*lx*lspec6*NCi* 
      sqrtxz2i)/4. + (11*x*z*lx*lspec6*NCi* 
      sqrtxz2i)/2. - (3*NC*x*li2spec5* 
      poly2i2*sqrtxz2i)/8. +  
   (3*NC*x*li2spec6*poly2i2*sqrtxz2i)/8. +  
   (3*NC*x*li2spec7*poly2i2* 
      sqrtxz2i)/8. - (3*NC*x* 
      li2spec8*poly2i2* 
      sqrtxz2i)/8. + (3*NC*x*lx*lspec5*poly2i2* 
      sqrtxz2i)/8. - (3*NC*x*lx*lspec6*poly2i2* 
      sqrtxz2i)/8. + (3*x*li2spec5*NCi* 
      poly2i2*sqrtxz2i)/8. -  
   (3*x*li2spec6*NCi*poly2i2* 
      sqrtxz2i)/8. - (3*x* 
      li2spec7*NCi* 
      poly2i2*sqrtxz2i)/8. +  
   (3*x*li2spec8*NCi* 
      poly2i2*sqrtxz2i)/8. -  
   (3*x*lx*lspec5*NCi*poly2i2*sqrtxz2i)/ 
    8. + (3*x*lx*lspec6*NCi*poly2i2* 
      sqrtxz2i)/8. - (3*NC*x*li2spec5* 
      poly2i*sqrtxz2i)/8. +  
   (3*NC*x*li2spec6*poly2i*sqrtxz2i)/8. +  
   (3*NC*x*li2spec7*poly2i* 
      sqrtxz2i)/8. - (3*NC*x* 
      li2spec8*poly2i* 
      sqrtxz2i)/8. + (3*NC*x*lx*lspec5*poly2i* 
      sqrtxz2i)/8. - (3*NC*x*lx*lspec6*poly2i* 
      sqrtxz2i)/8. + (x*li2spec5*NCi* 
      poly2i*sqrtxz2i)/8. -  
   (x*li2spec6*NCi*poly2i*sqrtxz2i)/ 
    8. - (x*li2spec7*NCi* 
      poly2i*sqrtxz2i)/8. +  
   (x*li2spec8*NCi* 
      poly2i*sqrtxz2i)/8. -  
   (x*lx*lspec5*NCi*poly2i*sqrtxz2i)/ 
    8. + (x*lx*lspec6*NCi*poly2i* 
      sqrtxz2i)/8. + (169*NC*xi)/72. +  
   (4*NC*lomx*xi)/3. + (NC*lx*xi)/2. +  
   (4*NC*lomz*xi)/3. + (11*NC*lz*xi)/6. -  
   (NCi*xi)/8. - (3*NC*lx*poly2i2*xi)/16. -  
   (3*NC*lz*poly2i2*xi)/16. +  
   (3*lx*NCi*poly2i2*xi)/16. +  
   (3*lz*NCi*poly2i2*xi)/16. -  
   (NC*poly2i*xi)/8. -  
   (5*NC*lx*poly2i*xi)/16. -  
   (5*NC*lz*poly2i*xi)/16. +  
   (NCi*poly2i*xi)/8. -  
   (3*lx*NCi*poly2i*xi)/16. -  
   (3*lz*NCi*poly2i*xi)/16. -  
   (7*NC*li2spec5*sqrtxz2i*xi)/32. +  
   (7*NC*li2spec6*sqrtxz2i*xi)/32. +  
   (7*NC*li2spec7*sqrtxz2i* 
      xi)/32. - (7*NC*li2spec8* 
      sqrtxz2i*xi)/32. +  
   (7*NC*lx*lspec5*sqrtxz2i*xi)/32. -  
   (7*NC*lx*lspec6*sqrtxz2i*xi)/32. -  
   (li2spec5*NCi*sqrtxz2i*xi)/32. +  
   (li2spec6*NCi*sqrtxz2i*xi)/32. +  
   (li2spec7*NCi* 
      sqrtxz2i*xi)/32. -  
   (li2spec8*NCi* 
      sqrtxz2i*xi)/32. +  
   (lx*lspec5*NCi*sqrtxz2i*xi)/32. -  
   (lx*lspec6*NCi*sqrtxz2i*xi)/32. +  
   (3*NC*li2spec5*poly2i2*sqrtxz2i* 
      xi)/32. - (3*NC*li2spec6*poly2i2* 
      sqrtxz2i*xi)/32. -  
   (3*NC*li2spec7*poly2i2* 
      sqrtxz2i*xi)/32. +  
   (3*NC*li2spec8*poly2i2* 
      sqrtxz2i*xi)/32. -  
   (3*NC*lx*lspec5*poly2i2*sqrtxz2i*xi)/ 
    32. + (3*NC*lx*lspec6*poly2i2*sqrtxz2i* 
      xi)/32. - (3*li2spec5*NCi* 
      poly2i2*sqrtxz2i*xi)/32. +  
   (3*li2spec6*NCi*poly2i2*sqrtxz2i* 
      xi)/32. + (3*li2spec7* 
      NCi*poly2i2*sqrtxz2i*xi)/32. -  
   (3*li2spec8*NCi* 
      poly2i2*sqrtxz2i*xi)/32. +  
   (3*lx*lspec5*NCi*poly2i2*sqrtxz2i* 
      xi)/32. - (3*lx*lspec6*NCi* 
      poly2i2*sqrtxz2i*xi)/32. +  
   (NC*li2spec5*poly2i*sqrtxz2i*xi)/ 
    8. - (NC*li2spec6*poly2i*sqrtxz2i* 
      xi)/8. - (NC*li2spec7* 
      poly2i*sqrtxz2i*xi)/8. +  
   (NC*li2spec8*poly2i* 
      sqrtxz2i*xi)/8. -  
   (NC*lx*lspec5*poly2i*sqrtxz2i*xi)/ 
    8. + (NC*lx*lspec6*poly2i*sqrtxz2i* 
      xi)/8. + (li2spec5*NCi*poly2i* 
      sqrtxz2i*xi)/8. -  
   (li2spec6*NCi*poly2i*sqrtxz2i* 
      xi)/8. - (li2spec7* 
      NCi*poly2i*sqrtxz2i*xi)/8. +  
   (li2spec8*NCi*poly2i* 
      sqrtxz2i*xi)/8. -  
   (lx*lspec5*NCi*poly2i*sqrtxz2i* 
      xi)/8. + (lx*lspec6*NCi*poly2i* 
      sqrtxz2i*xi)/8. + (1343*NC*x2)/72. +  
   8*NC*z*x2 + 4*NC*z*Li2mx*x2 - 2*NC*Li2x*x2 +  
   2*NC*z*Li2x*x2 - 8*NC*Li2z*x2 +  
   6*NC*li2omxzi*x2 + (95*NC*lomx*x2)/3. -  
   4*NC*z*lomx*x2 - 64*NC*lx*x2 +  
   4*NC*z*lx*x2 + 18*NC*lomx*lx*x2 +  
   2*NC*z*lomx*lx*x2 + 4*NC*z*lx*lopx*x2 +  
   (98*NC*lomz*x2)/3. - 4*NC*z*lomz*x2 -  
   16*NC*lomx*lomz*x2 + 24*NC*lx*lomz*x2 +  
   (80*NC*lz*x2)/3. - 10*NC*lomx*lz*x2 +  
   18*NC*lx*lz*x2 - 22*NC*lomz*lz*x2 -  
   (39*NCi*x2)/8. - 8*z*NCi*x2 +  
   4*Li2mx*NCi*x2 - 4*z*Li2mx*NCi*x2 -  
   2*Li2x*NCi*x2 - 2*z*Li2x*NCi*x2 +  
   8*Li2z*NCi*x2 -  
   2*li2omxzi*NCi*x2 -  
   lomx*NCi*x2 + 4*z*lomx*NCi*x2 -  
   (7*lx*NCi*x2)/2. - 4*z*lx*NCi*x2 -  
   22*lomx*lx*NCi*x2 -  
   2*z*lomx*lx*NCi*x2 +  
   4*lx*lopx*NCi*x2 -  
   4*z*lx*lopx*NCi*x2 -  
   4*lomz*NCi*x2 + 4*z*lomz*NCi*x2 +  
   14*lomx*lomz*NCi*x2 -  
   20*lx*lomz*NCi*x2 +  
   (11*lz*NCi*x2)/2. +  
   8*lomx*lz*NCi*x2 -  
   14*lx*lz*NCi*x2 +  
   14*lomz*lz*NCi*x2 + 3*NC*pi2*x2 -  
   (7*NCi*pi2*x2)/3. -  
   (9*NC*lx*poly2i2*x2)/16. +  
   (9*NC*lz*poly2i2*x2)/16. +  
   (9*lx*NCi*poly2i2*x2)/16. -  
   (9*lz*NCi*poly2i2*x2)/16. -  
   (NC*poly2i*x2)/4. - (5*NC*lx*poly2i*x2)/8. +  
   (5*NC*lz*poly2i*x2)/8. +  
   (NCi*poly2i*x2)/4. +  
   (5*lx*NCi*poly2i*x2)/8. -  
   (5*lz*NCi*poly2i*x2)/8. +  
   (129*NC*li2spec5*sqrtxz2i*x2)/8. -  
   (129*NC*z*li2spec5*sqrtxz2i*x2)/4. -  
   (129*NC*li2spec6*sqrtxz2i*x2)/8. +  
   (129*NC*z*li2spec6*sqrtxz2i*x2)/4. -  
   (129*NC*li2spec7*sqrtxz2i* 
      x2)/8. + (129*NC*z*li2spec7*sqrtxz2i*x2)/4. +  
   (129*NC*li2spec8*sqrtxz2i* 
      x2)/8. - (129*NC*z*li2spec8*sqrtxz2i*x2)/4. -  
   (129*NC*lx*lspec5*sqrtxz2i*x2)/8. +  
   (129*NC*z*lx*lspec5*sqrtxz2i*x2)/4. +  
   (129*NC*lx*lspec6*sqrtxz2i*x2)/8. -  
   (129*NC*z*lx*lspec6*sqrtxz2i*x2)/4. -  
   (29*li2spec5*NCi*sqrtxz2i*x2)/8. +  
   (29*z*li2spec5*NCi*sqrtxz2i*x2)/ 
    4. + (29*li2spec6*NCi*sqrtxz2i* 
      x2)/8. - (29*z*li2spec6*NCi* 
      sqrtxz2i*x2)/4. +  
   (29*li2spec7*NCi* 
      sqrtxz2i*x2)/8. -  
   (29*z*li2spec7*NCi* 
      sqrtxz2i*x2)/4. -  
   (29*li2spec8*NCi* 
      sqrtxz2i*x2)/8. +  
   (29*z*li2spec8*NCi* 
      sqrtxz2i*x2)/4. +  
   (29*lx*lspec5*NCi*sqrtxz2i*x2)/8. -  
   (29*z*lx*lspec5*NCi*sqrtxz2i*x2)/ 
    4. - (29*lx*lspec6*NCi*sqrtxz2i*x2)/ 
    8. + (29*z*lx*lspec6*NCi*sqrtxz2i* 
      x2)/4. - (9*NC*lx*poly2i2*x3)/16. -  
   (9*NC*lz*poly2i2*x3)/16. +  
   (9*lx*NCi*poly2i2*x3)/16. +  
   (9*lz*NCi*poly2i2*x3)/16. -  
   (NC*poly2i*x3)/8. - (NC*lx*poly2i*x3)/16. -  
   (NC*lz*poly2i*x3)/16. +  
   (NCi*poly2i*x3)/8. +  
   (9*lx*NCi*poly2i*x3)/16. +  
   (9*lz*NCi*poly2i*x3)/16. +  
   (193*NC*li2spec5*sqrtxz2i*x3)/32. -  
   (193*NC*li2spec6*sqrtxz2i*x3)/32. -  
   (193*NC*li2spec7*sqrtxz2i* 
      x3)/32. + (193*NC*li2spec8*sqrtxz2i*x3)/32. -  
   (193*NC*lx*lspec5*sqrtxz2i*x3)/32. +  
   (193*NC*lx*lspec6*sqrtxz2i*x3)/32. -  
   (73*li2spec5*NCi*sqrtxz2i*x3)/ 
    32. + (73*li2spec6*NCi*sqrtxz2i* 
      x3)/32. + (73*li2spec7* 
      NCi*sqrtxz2i*x3)/32. -  
   (73*li2spec8*NCi* 
      sqrtxz2i*x3)/32. +  
   (73*lx*lspec5*NCi*sqrtxz2i*x3)/32. -  
   (73*lx*lspec6*NCi*sqrtxz2i*x3)/32. +  
   (9*NC*li2spec5*poly2i2*sqrtxz2i*x3)/ 
    16. - (9*NC*li2spec6*poly2i2*sqrtxz2i* 
      x3)/16. - (9*NC*li2spec7* 
      poly2i2*sqrtxz2i*x3)/16. +  
   (9*NC*li2spec8*poly2i2* 
      sqrtxz2i*x3)/16. -  
   (9*NC*lx*lspec5*poly2i2*sqrtxz2i*x3)/ 
    16. + (9*NC*lx*lspec6*poly2i2*sqrtxz2i* 
      x3)/16. - (9*li2spec5*NCi* 
      poly2i2*sqrtxz2i*x3)/16. +  
   (9*li2spec6*NCi*poly2i2*sqrtxz2i* 
      x3)/16. + (9*li2spec7* 
      NCi*poly2i2*sqrtxz2i*x3)/16. -  
   (9*li2spec8*NCi* 
      poly2i2*sqrtxz2i*x3)/16. +  
   (9*lx*lspec5*NCi*poly2i2*sqrtxz2i* 
      x3)/16. - (9*lx*lspec6*NCi*poly2i2* 
      sqrtxz2i*x3)/16. +  
   (3*NC*li2spec5*poly2i*sqrtxz2i*x3)/ 
    8. - (3*NC*li2spec6*poly2i*sqrtxz2i* 
      x3)/8. - (3*NC*li2spec7* 
      poly2i*sqrtxz2i*x3)/8. +  
   (3*NC*li2spec8*poly2i* 
      sqrtxz2i*x3)/8. -  
   (3*NC*lx*lspec5*poly2i*sqrtxz2i*x3)/ 
    8. + (3*NC*lx*lspec6*poly2i*sqrtxz2i* 
      x3)/8. - (5*li2spec5*NCi*poly2i* 
      sqrtxz2i*x3)/8. +  
   (5*li2spec6*NCi*poly2i*sqrtxz2i* 
      x3)/8. + (5*li2spec7* 
      NCi*poly2i*sqrtxz2i*x3)/8. -  
   (5*li2spec8*NCi* 
      poly2i*sqrtxz2i*x3)/8. +  
   (5*lx*lspec5*NCi*poly2i*sqrtxz2i* 
      x3)/8. - (5*lx*lspec6*NCi*poly2i* 
      sqrtxz2i*x3)/8. +  
   (9*NC*lx*poly2i2*x4)/16. -  
   (9*NC*lz*poly2i2*x4)/16. -  
   (9*lx*NCi*poly2i2*x4)/16. +  
   (9*lz*NCi*poly2i2*x4)/16. +  
   (NC*poly2i*x4)/8. +  
   (3*NC*lx*poly2i*x4)/16. -  
   (3*NC*lz*poly2i*x4)/16. -  
   (NCi*poly2i*x4)/8. -  
   (11*lx*NCi*poly2i*x4)/16. +  
   (11*lz*NCi*poly2i*x4)/16. +  
   (3*NC*lx*poly2i2*x5)/16. +  
   (3*NC*lz*poly2i2*x5)/16. -  
   (3*lx*NCi*poly2i2*x5)/16. -  
   (3*lz*NCi*poly2i2*x5)/16. -  
   (3*NC*li2spec5*poly2i2*sqrtxz2i*x5)/ 
    8. + (3*NC*li2spec6*poly2i2*sqrtxz2i* 
      x5)/8. + (3*NC*li2spec7* 
      poly2i2*sqrtxz2i*x5)/8. -  
   (3*NC*li2spec8*poly2i2* 
      sqrtxz2i*x5)/8. +  
   (3*NC*lx*lspec5*poly2i2*sqrtxz2i*x5)/ 
    8. - (3*NC*lx*lspec6*poly2i2*sqrtxz2i* 
      x5)/8. + (3*li2spec5*NCi*poly2i2* 
      sqrtxz2i*x5)/8. -  
   (3*li2spec6*NCi*poly2i2*sqrtxz2i* 
      x5)/8. - (3*li2spec7* 
      NCi*poly2i2*sqrtxz2i*x5)/8. +  
   (3*li2spec8*NCi* 
      poly2i2*sqrtxz2i*x5)/8. -  
   (3*lx*lspec5*NCi*poly2i2*sqrtxz2i* 
      x5)/8. + (3*lx*lspec6*NCi*poly2i2* 
      sqrtxz2i*x5)/8. -  
   (NC*li2spec5*poly2i*sqrtxz2i*x5)/ 
    8. + (NC*li2spec6*poly2i*sqrtxz2i* 
      x5)/8. + (NC*li2spec7* 
      poly2i*sqrtxz2i*x5)/8. -  
   (NC*li2spec8*poly2i* 
      sqrtxz2i*x5)/8. +  
   (NC*lx*lspec5*poly2i*sqrtxz2i*x5)/ 
    8. - (NC*lx*lspec6*poly2i*sqrtxz2i* 
      x5)/8. + (3*li2spec5*NCi*poly2i* 
      sqrtxz2i*x5)/8. -  
   (3*li2spec6*NCi*poly2i*sqrtxz2i* 
      x5)/8. - (3*li2spec7* 
      NCi*poly2i*sqrtxz2i*x5)/8. +  
   (3*li2spec8*NCi* 
      poly2i*sqrtxz2i*x5)/8. -  
   (3*lx*lspec5*NCi*poly2i*sqrtxz2i* 
      x5)/8. + (3*lx*lspec6*NCi*poly2i* 
      sqrtxz2i*x5)/8. -  
   (3*NC*lx*poly2i2*x6)/16. +  
   (3*NC*lz*poly2i2*x6)/16. +  
   (3*lx*NCi*poly2i2*x6)/16. -  
   (3*lz*NCi*poly2i2*x6)/16. +  
   (3*NC*li2spec5*poly2i2*sqrtxz2i*x7)/ 
    32. - (3*NC*li2spec6*poly2i2*sqrtxz2i* 
      x7)/32. - (3*NC*li2spec7* 
      poly2i2*sqrtxz2i*x7)/32. +  
   (3*NC*li2spec8*poly2i2* 
      sqrtxz2i*x7)/32. -  
   (3*NC*lx*lspec5*poly2i2*sqrtxz2i*x7)/ 
    32. + (3*NC*lx*lspec6*poly2i2*sqrtxz2i* 
      x7)/32. - (3*li2spec5*NCi* 
      poly2i2*sqrtxz2i*x7)/32. +  
   (3*li2spec6*NCi*poly2i2*sqrtxz2i* 
      x7)/32. + (3*li2spec7* 
      NCi*poly2i2*sqrtxz2i*x7)/32. -  
   (3*li2spec8*NCi* 
      poly2i2*sqrtxz2i*x7)/32. +  
   (3*lx*lspec5*NCi*poly2i2*sqrtxz2i* 
      x7)/32. - (3*lx*lspec6*NCi*poly2i2* 
      sqrtxz2i*x7)/32. - 2*x*NCi*omzi +  
   4*x*lomx*NCi*omzi -  
   6*x*lx*NCi*omzi +  
   2*x*lomz*NCi*omzi +  
   4*x*lz*NCi*omzi +  
   2*NCi*x2*omzi -  
   4*lomx*NCi*x2*omzi +  
   6*lx*NCi*x2*omzi -  
   2*lomz*NCi*x2*omzi -  
   4*lz*NCi*x2*omzi -  
   NC*lx*x2*omxmzi2 +  
   NC*lomz*x2*omxmzi2 +  
   lx*NCi*x2*omxmzi2 -  
   lomz*NCi*x2*omxmzi2 +  
   2*NC*lx*x3*omxmzi2 -  
   2*NC*lomz*x3*omxmzi2 -  
   2*lx*NCi*x3*omxmzi2 +  
   2*lomz*NCi*x3*omxmzi2 -  
   NC*lx*x4*omxmzi2 +  
   NC*lomz*x4*omxmzi2 +  
   lx*NCi*x4*omxmzi2 -  
   lomz*NCi*x4*omxmzi2 -  
   NC*x*omxmzi + NC*x*lx*omxmzi -  
   NC*x*lomz*omxmzi + x*NCi*omxmzi -  
   x*lx*NCi*omxmzi +  
   x*lomz*NCi*omxmzi +  
   2*NC*x2*omxmzi - 3*NC*lx*x2*omxmzi +  
   3*NC*lomz*x2*omxmzi -  
   2*NCi*x2*omxmzi +  
   lx*NCi*x2*omxmzi -  
   lomz*NCi*x2*omxmzi -  
   NC*x3*omxmzi + 2*NC*lx*x3*omxmzi -  
   2*NC*lomz*x3*omxmzi +  
   NCi*x3*omxmzi + NC*lx*x4*xmzi2 -  
   NC*lz*x4*xmzi2 -  
   lx*NCi*x4*xmzi2 +  
   lz*NCi*x4*xmzi2 -  
   2*lx*NCi*x2*xmzi +  
   2*lz*NCi*x2*xmzi -  
   3*NC*lx*x3*xmzi + 3*NC*lz*x3*xmzi +  
   5*lx*NCi*x3*xmzi -  
   5*lz*NCi*x3*xmzi + (NC*zi)/2. +  
   (NC*x*zi)/2. + NC*x*lx*zi - (NCi*zi)/2. -  
   (5*x*NCi*zi)/2. + 4*x*lomx*NCi*zi -  
   7*x*lx*NCi*zi + 4*x*lomz*NCi*zi +  
   2*x*lz*NCi*zi - NC*x2*zi +  
   3*NCi*x2*zi -  
   4*lomx*NCi*x2*zi +  
   6*lx*NCi*x2*zi -  
   4*lomz*NCi*x2*zi -  
   2*lz*NCi*x2*zi +  
   (51*NC*x*li2spec5*sqrtxz2i*z2)/2. -  
   (51*NC*x*li2spec6*sqrtxz2i*z2)/2. -  
   (51*NC*x*li2spec7*sqrtxz2i* 
      z2)/2. + (51*NC*x*li2spec8*sqrtxz2i*z2)/2. -  
   (51*NC*x*lx*lspec5*sqrtxz2i*z2)/2. +  
   (51*NC*x*lx*lspec6*sqrtxz2i*z2)/2. -  
   (11*x*li2spec5*NCi*sqrtxz2i*z2)/ 
    2. + (11*x*li2spec6*NCi*sqrtxz2i* 
      z2)/2. + (11*x*li2spec7* 
      NCi*sqrtxz2i*z2)/2. -  
   (11*x*li2spec8*NCi* 
      sqrtxz2i*z2)/2. +  
   (11*x*lx*lspec5*NCi*sqrtxz2i*z2)/ 
    2. - (11*x*lx*lspec6*NCi*sqrtxz2i* 
      z2)/2. - NC*x*lx*xmzi2*z3 +  
   NC*x*lz*xmzi2*z3 +  
   x*lx*NCi*xmzi2*z3 -  
   x*lz*NCi*xmzi2*z3 + 6*NC*x*lomx2 -  
   8*x*NCi*lomx2 - 6*NC*x2*lomx2 +  
   8*NCi*x2*lomx2 + 29*NC*x*lx2 -  
   NC*x*z*lx2 - 15*x*NCi*lx2 +  
   x*z*NCi*lx2 - 13*NC*x2*lx2 -  
   3*NC*z*x2*lx2 + 13*NCi*x2*lx2 +  
   3*z*NCi*x2*lx2 + 10*NC*x*lomz2 -  
   7*x*NCi*lomz2 - 10*NC*x2*lomz2 +  
   7*NCi*x2*lomz2 + 3*NC*x*lz2 -  
   2*x*NCi*lz2 - 3*NC*x2*lz2 +  
   2*NCi*x2*lz2 +  
   (-4*NC*x - 4*NC*x*li2spec13 +  
      4*NC*x*li2spec14 +  
      4*NC*x*li2spec11 + 4*NC*x*lomx -  
      12*NC*x*lx + 16*NC*x*lomx*lx + 8*NC*x*lomz -  
      12*NC*x*lomx*lomz + 20*NC*x*lx*lomz -  
      4*NC*x*lx*lomxmz + 4*NC*x*lomz*lomxmz -  
      4*NC*x*lx*lxmz + 8*NC*x*lz - 8*NC*x*lomx*lz +  
      20*NC*x*lx*lz - 16*NC*x*lomz*lz +  
      4*NC*x*lxmz*lz - 2*x*li2spec13*NCi +  
      2*x*li2spec11*NCi +  
      2*x*li2spec15*NCi -  
      2*x*li2spec12*NCi -  
      22*x*lomx*lx*NCi +  
      14*x*lomx*lomz*NCi -  
      20*x*lx*lomz*NCi +  
      2*x*lx*lomxmz*NCi -  
      2*x*lomz*lomxmz*NCi +  
      2*x*lx*lxmz*NCi + 12*x*lomx*lz*NCi -  
      16*x*lx*lz*NCi + 8*x*lomz*lz*NCi -  
      2*x*lxmz*lz*NCi + (2*NC*x*pi2)/3. -  
      (7*x*NCi*pi2)/3. + 4*NC*x2 +  
      4*NC*li2spec13*x2 -  
      4*NC*li2spec14*x2 -  
      4*NC*li2spec11*x2 - 4*NC*lomx*x2 +  
      12*NC*lx*x2 - 16*NC*lomx*lx*x2 -  
      8*NC*lomz*x2 + 12*NC*lomx*lomz*x2 -  
      20*NC*lx*lomz*x2 +  
      4*NC*lx*lomxmz*x2 -  
      4*NC*lomz*lomxmz*x2 +  
      4*NC*lx*lxmz*x2 - 8*NC*lz*x2 +  
      8*NC*lomx*lz*x2 - 20*NC*lx*lz*x2 +  
      16*NC*lomz*lz*x2 - 4*NC*lxmz*lz*x2 +  
      2*li2spec13*NCi*x2 -  
      2*li2spec11*NCi*x2 -  
      2*li2spec15*NCi*x2 +  
      2*li2spec12*NCi*x2 +  
      22*lomx*lx*NCi*x2 -  
      14*lomx*lomz*NCi*x2 +  
      20*lx*lomz*NCi*x2 -  
      2*lx*lomxmz*NCi*x2 +  
      2*lomz*lomxmz*NCi*x2 -  
      2*lx*lxmz*NCi*x2 -  
      12*lomx*lz*NCi*x2 +  
      16*lx*lz*NCi*x2 -  
      8*lomz*lz*NCi*x2 +  
      2*lxmz*lz*NCi*x2 -  
      (2*NC*pi2*x2)/3. + (7*NCi*pi2*x2)/3. +  
      2*x*NCi*omzi -  
      4*x*lomx*NCi*omzi +  
      6*x*lx*NCi*omzi -  
      2*x*lomz*NCi*omzi -  
      4*x*lz*NCi*omzi -  
      2*NCi*x2*omzi +  
      4*lomx*NCi*x2*omzi -  
      6*lx*NCi*x2*omzi +  
      2*lomz*NCi*x2*omzi +  
      4*lz*NCi*x2*omzi + 2*x*NCi*zi -  
      4*x*lomx*NCi*zi + 6*x*lx*NCi*zi -  
      4*x*lomz*NCi*zi - 2*x*lz*NCi*zi -  
      2*NCi*x2*zi +  
      4*lomx*NCi*x2*zi -  
      6*lx*NCi*x2*zi +  
      4*lomz*NCi*x2*zi +  
      2*lz*NCi*x2*zi - 2*NC*x*lomx2 +  
      8*x*NCi*lomx2 + 2*NC*x2*lomx2 -  
      8*NCi*x2*lomx2 - 12*NC*x*lx2 +  
      13*x*NCi*lx2 + 12*NC*x2*lx2 -  
      13*NCi*x2*lx2 - 6*NC*x*lomz2 +  
      6*x*NCi*lomz2 + 6*NC*x2*lomz2 -  
      6*NCi*x2*lomz2 - 8*NC*x*lz2 +  
      5*x*NCi*lz2 + 8*NC*x2*lz2 -  
      5*NCi*x2*lz2)*Tu1 +  
   (-4*NC*x + 4*NC*x*li2spec11 +  
      4*NC*x*li2spec16 -  
      4*NC*x*li2spec18 + 4*NC*x*lomx -  
      12*NC*x*lx + 12*NC*x*lomx*lx + 8*NC*x*lomz -  
      8*NC*x*lomx*lomz + 24*NC*x*lx*lomz -  
      4*NC*x*lx*lxmz + 8*NC*x*lz - 8*NC*x*lomx*lz +  
      24*NC*x*lx*lz - 20*NC*x*lomz*lz +  
      4*NC*x*lxmz*lz - 4*NC*x*lx*lmopxpz +  
      4*NC*x*lomz*lmopxpz +  
      2*x*li2spec11*NCi -  
      2*x*li2spec12*NCi +  
      2*x*li2spec16*NCi -  
      2*x*li2spec17*NCi -  
      20*x*lomx*lx*NCi +  
      12*x*lomx*lomz*NCi -  
      18*x*lx*lomz*NCi + 2*x*lx*lxmz*NCi +  
      12*x*lomx*lz*NCi - 18*x*lx*lz*NCi +  
      10*x*lomz*lz*NCi - 2*x*lxmz*lz*NCi +  
      2*x*lx*lmopxpz*NCi -  
      2*x*lomz*lmopxpz*NCi + (2*NC*x*pi2)/3. -  
      (7*x*NCi*pi2)/3. + 4*NC*x2 -  
      4*NC*li2spec11*x2 -  
      4*NC*li2spec16*x2 +  
      4*NC*li2spec18*x2 -  
      4*NC*lomx*x2 + 12*NC*lx*x2 -  
      12*NC*lomx*lx*x2 - 8*NC*lomz*x2 +  
      8*NC*lomx*lomz*x2 -  
      24*NC*lx*lomz*x2 + 4*NC*lx*lxmz*x2 -  
      8*NC*lz*x2 + 8*NC*lomx*lz*x2 -  
      24*NC*lx*lz*x2 + 20*NC*lomz*lz*x2 -  
      4*NC*lxmz*lz*x2 +  
      4*NC*lx*lmopxpz*x2 -  
      4*NC*lomz*lmopxpz*x2 -  
      2*li2spec11*NCi*x2 +  
      2*li2spec12*NCi*x2 -  
      2*li2spec16*NCi*x2 +  
      2*li2spec17*NCi*x2 +  
      20*lomx*lx*NCi*x2 -  
      12*lomx*lomz*NCi*x2 +  
      18*lx*lomz*NCi*x2 -  
      2*lx*lxmz*NCi*x2 -  
      12*lomx*lz*NCi*x2 +  
      18*lx*lz*NCi*x2 -  
      10*lomz*lz*NCi*x2 +  
      2*lxmz*lz*NCi*x2 -  
      2*lx*lmopxpz*NCi*x2 +  
      2*lomz*lmopxpz*NCi*x2 -  
      (2*NC*pi2*x2)/3. + (7*NCi*pi2*x2)/3. +  
      2*x*NCi*omzi -  
      4*x*lomx*NCi*omzi +  
      6*x*lx*NCi*omzi -  
      2*x*lomz*NCi*omzi -  
      4*x*lz*NCi*omzi -  
      2*NCi*x2*omzi +  
      4*lomx*NCi*x2*omzi -  
      6*lx*NCi*x2*omzi +  
      2*lomz*NCi*x2*omzi +  
      4*lz*NCi*x2*omzi + 2*x*NCi*zi -  
      4*x*lomx*NCi*zi + 6*x*lx*NCi*zi -  
      4*x*lomz*NCi*zi - 2*x*lz*NCi*zi -  
      2*NCi*x2*zi +  
      4*lomx*NCi*x2*zi -  
      6*lx*NCi*x2*zi +  
      4*lomz*NCi*x2*zi +  
      2*lz*NCi*x2*zi - 2*NC*x*lomx2 +  
      8*x*NCi*lomx2 + 2*NC*x2*lomx2 -  
      8*NCi*x2*lomx2 - 14*NC*x*lx2 +  
      12*x*NCi*lx2 + 14*NC*x2*lx2 -  
      12*NCi*x2*lx2 - 8*NC*x*lomz2 +  
      5*x*NCi*lomz2 + 8*NC*x2*lomz2 -  
      5*NCi*x2*lomz2 - 8*NC*x*lz2 +  
      5*x*NCi*lz2 + 8*NC*x2*lz2 -  
      5*NCi*x2*lz2)*Tu2 +  
   (-4*NC*x - 4*NC*x*li2spec9 -  
      4*NC*x*li2spec13 -  
      4*NC*x*li2spec18 + 4*NC*x*lomx -  
      12*NC*x*lx + 12*NC*x*lomx*lx + 8*NC*x*lomz -  
      4*NC*x*lomx*lomz + 24*NC*x*lx*lomz -  
      4*NC*x*lx*lomxmz + 4*NC*x*lomz*lomxmz +  
      8*NC*x*lz - 4*NC*x*lomx*lz + 24*NC*x*lx*lz -  
      20*NC*x*lomz*lz - 4*NC*x*lx*lmxpz +  
      4*NC*x*lz*lmxpz - 2*x*li2spec9*NCi -  
      2*x*li2spec13*NCi +  
      2*x*li2spec15*NCi +  
      2*x*li2spec10*NCi -  
      24*x*lomx*lx*NCi +  
      14*x*lomx*lomz*NCi -  
      18*x*lx*lomz*NCi +  
      2*x*lx*lomxmz*NCi -  
      2*x*lomz*lomxmz*NCi +  
      14*x*lomx*lz*NCi - 18*x*lx*lz*NCi +  
      6*x*lomz*lz*NCi + 2*x*lx*lmxpz*NCi -  
      2*x*lz*lmxpz*NCi + (10*NC*x*pi2)/3. -  
      (7*x*NCi*pi2)/3. + 4*NC*x2 +  
      4*NC*li2spec9*x2 +  
      4*NC*li2spec13*x2 +  
      4*NC*li2spec18*x2 -  
      4*NC*lomx*x2 + 12*NC*lx*x2 -  
      12*NC*lomx*lx*x2 - 8*NC*lomz*x2 +  
      4*NC*lomx*lomz*x2 -  
      24*NC*lx*lomz*x2 +  
      4*NC*lx*lomxmz*x2 -  
      4*NC*lomz*lomxmz*x2 - 8*NC*lz*x2 +  
      4*NC*lomx*lz*x2 - 24*NC*lx*lz*x2 +  
      20*NC*lomz*lz*x2 + 4*NC*lx*lmxpz*x2 -  
      4*NC*lz*lmxpz*x2 +  
      2*li2spec9*NCi*x2 +  
      2*li2spec13*NCi*x2 -  
      2*li2spec15*NCi*x2 -  
      2*li2spec10*NCi*x2 +  
      24*lomx*lx*NCi*x2 -  
      14*lomx*lomz*NCi*x2 +  
      18*lx*lomz*NCi*x2 -  
      2*lx*lomxmz*NCi*x2 +  
      2*lomz*lomxmz*NCi*x2 -  
      14*lomx*lz*NCi*x2 +  
      18*lx*lz*NCi*x2 -  
      6*lomz*lz*NCi*x2 -  
      2*lx*lmxpz*NCi*x2 +  
      2*lz*lmxpz*NCi*x2 -  
      (10*NC*pi2*x2)/3. + (7*NCi*pi2*x2)/3. +  
      2*x*NCi*omzi -  
      4*x*lomx*NCi*omzi +  
      6*x*lx*NCi*omzi -  
      2*x*lomz*NCi*omzi -  
      4*x*lz*NCi*omzi -  
      2*NCi*x2*omzi +  
      4*lomx*NCi*x2*omzi -  
      6*lx*NCi*x2*omzi +  
      2*lomz*NCi*x2*omzi +  
      4*lz*NCi*x2*omzi + 2*x*NCi*zi -  
      4*x*lomx*NCi*zi + 6*x*lx*NCi*zi -  
      4*x*lomz*NCi*zi - 2*x*lz*NCi*zi -  
      2*NCi*x2*zi +  
      4*lomx*NCi*x2*zi -  
      6*lx*NCi*x2*zi +  
      4*lomz*NCi*x2*zi +  
      2*lz*NCi*x2*zi - 6*NC*x*lomx2 +  
      8*x*NCi*lomx2 + 6*NC*x2*lomx2 -  
      8*NCi*x2*lomx2 - 14*NC*x*lx2 +  
      14*x*NCi*lx2 + 14*NC*x2*lx2 -  
      14*NCi*x2*lx2 - 10*NC*x*lomz2 +  
      6*x*NCi*lomz2 + 10*NC*x2*lomz2 -  
      6*NCi*x2*lomz2 - 10*NC*x*lz2 +  
      6*x*NCi*lz2 + 10*NC*x2*lz2 -  
      6*NCi*x2*lz2)*Tu3 +  
   (-4*NC*x - 4*NC*x*li2spec9 +  
      4*NC*x*li2spec14 +  
      4*NC*x*li2spec16 + 4*NC*x*lomx - 12*NC*x*lx +  
      16*NC*x*lomx*lx + 8*NC*x*lomz -  
      8*NC*x*lomx*lomz + 20*NC*x*lx*lomz +  
      8*NC*x*lz - 12*NC*x*lomx*lz + 20*NC*x*lx*lz -  
      16*NC*x*lomz*lz - 4*NC*x*lx*lmxpz +  
      4*NC*x*lz*lmxpz - 4*NC*x*lx*lmopxpz +  
      4*NC*x*lomz*lmopxpz -  
      2*x*li2spec9*NCi +  
      2*x*li2spec16*NCi +  
      2*x*li2spec10*NCi -  
      2*x*li2spec17*NCi -  
      22*x*lomx*lx*NCi +  
      12*x*lomx*lomz*NCi -  
      16*x*lx*lomz*NCi + 14*x*lomx*lz*NCi -  
      20*x*lx*lz*NCi + 8*x*lomz*lz*NCi +  
      2*x*lx*lmxpz*NCi - 2*x*lz*lmxpz*NCi +  
      2*x*lx*lmopxpz*NCi -  
      2*x*lomz*lmopxpz*NCi + (2*NC*x*pi2)/3. -  
      (7*x*NCi*pi2)/3. + 4*NC*x2 +  
      4*NC*li2spec9*x2 -  
      4*NC*li2spec14*x2 -  
      4*NC*li2spec16*x2 - 4*NC*lomx*x2 +  
      12*NC*lx*x2 - 16*NC*lomx*lx*x2 -  
      8*NC*lomz*x2 + 8*NC*lomx*lomz*x2 -  
      20*NC*lx*lomz*x2 - 8*NC*lz*x2 +  
      12*NC*lomx*lz*x2 - 20*NC*lx*lz*x2 +  
      16*NC*lomz*lz*x2 + 4*NC*lx*lmxpz*x2 -  
      4*NC*lz*lmxpz*x2 +  
      4*NC*lx*lmopxpz*x2 -  
      4*NC*lomz*lmopxpz*x2 +  
      2*li2spec9*NCi*x2 -  
      2*li2spec16*NCi*x2 -  
      2*li2spec10*NCi*x2 +  
      2*li2spec17*NCi*x2 +  
      22*lomx*lx*NCi*x2 -  
      12*lomx*lomz*NCi*x2 +  
      16*lx*lomz*NCi*x2 -  
      14*lomx*lz*NCi*x2 +  
      20*lx*lz*NCi*x2 -  
      8*lomz*lz*NCi*x2 -  
      2*lx*lmxpz*NCi*x2 +  
      2*lz*lmxpz*NCi*x2 -  
      2*lx*lmopxpz*NCi*x2 +  
      2*lomz*lmopxpz*NCi*x2 -  
      (2*NC*pi2*x2)/3. + (7*NCi*pi2*x2)/3. +  
      2*x*NCi*omzi -  
      4*x*lomx*NCi*omzi +  
      6*x*lx*NCi*omzi -  
      2*x*lomz*NCi*omzi -  
      4*x*lz*NCi*omzi -  
      2*NCi*x2*omzi +  
      4*lomx*NCi*x2*omzi -  
      6*lx*NCi*x2*omzi +  
      2*lomz*NCi*x2*omzi +  
      4*lz*NCi*x2*omzi + 2*x*NCi*zi -  
      4*x*lomx*NCi*zi + 6*x*lx*NCi*zi -  
      4*x*lomz*NCi*zi - 2*x*lz*NCi*zi -  
      2*NCi*x2*zi +  
      4*lomx*NCi*x2*zi -  
      6*lx*NCi*x2*zi +  
      4*lomz*NCi*x2*zi +  
      2*lz*NCi*x2*zi - 2*NC*x*lomx2 +  
      8*x*NCi*lomx2 + 2*NC*x2*lomx2 -  
      8*NCi*x2*lomx2 - 12*NC*x*lx2 +  
      13*x*NCi*lx2 + 12*NC*x2*lx2 -  
      13*NCi*x2*lx2 - 8*NC*x*lomz2 +  
      5*x*NCi*lomz2 + 8*NC*x2*lomz2 -  
      5*NCi*x2*lomz2 - 6*NC*x*lz2 +  
      6*x*NCi*lz2 + 6*NC*x2*lz2 -  
      6*NCi*x2*lz2)*Tu4;
};
double C2LQ2QNS_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double sqrtxz1  = sqrt(1 - 2*z + z*z + 4*x*z);
  const double poly2    = 1 + 2*x + x*x - 4*x*z;
  const double poly2i   = 1. / poly2;
  const double poly2i2  = poly2i * poly2i;
  const double sqrtxz2  = sqrt(poly2);
  const double sqrtxz2i = 1. / sqrtxz2;
  const double sqrtxz3  = sqrt(x/z);
  const double NC   = apfel::NC;
  const double NC2  = NC * NC;
  const double NCi  = 1. / NC;
  const double NCi2 = NCi * NCi;
  const double l2  = log(2);
  const double l22 = l2 * l2;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double x2 = x * x;
  const double x3 = x * x2;
  const double x4 = x * x3;
  const double x5 = x * x4;
  const double x6 = x * x5;
  const double x7 = x * x6;
  const double xi  = 1. / x;
  const double z2 = z * z;
  const double z3 = z * z2;
  const double zi  = 1. / z;
  const double omxi = 1. / ( 1 - x );
  const double opxi = 1. / ( 1 + x );
  const double omzi  = 1. / ( 1 - z );
  const double lx  = log(x);
  const double lx2 = lx * lx;
  const double lz  = log(z);
  const double lz2 = lz * lz;
  const double lomx  = log(1 - x);
  const double lomx2 = lomx * lomx;
  const double lomz  = log(1 - z);
  const double lomz2 = lomz * lomz;
  const double lopx  = log(1 + x);
  const double Li2x  = apfel::dilog(x);
  const double Li2mx = apfel::dilog(-x);
  const double Li2z  = apfel::dilog(z);
  const double xmzi  = 1. / ( x - z );
  const double xmzi2 = xmzi * xmzi;
  const double lomxmz  = log(1 - x - z);
  const double lmopxpz = log(-1 + x + z);
  const double lxmz    = log(x - z);
  const double lmxpz   = log(-x + z);
  const double lxpz    = log(x + z);
  const double lopxz   = log(1 + x*z);
  const double lopxzi  = log(1 + x*zi);
  const double li2omxzi = apfel::dilog(1 - x*zi);
  const double lspec1   = log(1 + sqrtxz1 - z);
  const double lspec1_2 = lspec1 * lspec1;
  const double lspec2   = log(1 + sqrtxz1 + z);
  const double lspec3   = log(sqrtxz3);
  const double lspec4   = log(sqrtxz3*z);
  const double lspec5   = log(1 - sqrtxz2 + x);
  const double lspec6   = log(1 + sqrtxz2 + x);
  const double li2spec1  = apfel::dilog(0.5 - sqrtxz1/2. - z/2.);
  const double li2spec2  = apfel::dilog(0.5 - sqrtxz1/2. + z/2.);
  const double li2spec3  = apfel::dilog(0.5 - zi/2. - (sqrtxz1*zi)/2.);
  const double li2spec4  = apfel::dilog(0.5 + zi/2. - (sqrtxz1*zi)/2.);
  const double li2spec5  = apfel::dilog(0.5 - sqrtxz2/2. - x/2.);
  const double li2spec6  = apfel::dilog(0.5 + sqrtxz2/2. - x/2.);
  const double li2spec7  = apfel::dilog(0.5 - xi/2. - (sqrtxz2*xi)/2.);
  const double li2spec8  = apfel::dilog(0.5 - xi/2. + (sqrtxz2*xi)/2.);
  const double li2spec9  = apfel::dilog((1 - z)*omxi);
  const double li2spec10 = apfel::dilog(x*(1 - z)*omxi*zi);
  const double li2spec11 = apfel::dilog((1 - x)*omzi);
  const double li2spec12 = apfel::dilog((1 - x)*z*xi*omzi);
  const double li2spec13 = apfel::dilog(z*omxi);
  const double li2spec14 = apfel::dilog((1 - z)*z*omxi*xi);
  const double li2spec15 = apfel::dilog(x*z*omxi*omzi);
  const double li2spec16 = apfel::dilog((1 - x)*zi);
  const double li2spec17 = apfel::dilog((1 - x)*(1 - z)*xi*zi);
  const double li2spec18 = apfel::dilog((1 - x)*x*omzi*zi);
  const double li2spec19 = apfel::dilog(-(x*z));
  const double li2spec20 = apfel::dilog(-(x*zi));
  const double atanspec1 = atan(sqrtxz3);
  const double atanspec2 = atan(sqrtxz3*z);
  const double itani1 = InvTanInt(-sqrtxz3);
  const double itani2 = InvTanInt(sqrtxz3);
  const double itani3 = InvTanInt(sqrtxz3*z);
  const double Tr1 = (z < 1 - x ? 1 : 0);
  const double Tr2 = (z > 1 - x ? 1 : 0);
  const double Tu1 = (z < 1 - x && z < x ? 1 : 0);
  const double Tu2 = (z > 1 - x && z < x ? 1 : 0);
  const double Tu3 = (z < 1 - x && z > x ? 1 : 0);
  const double Tu4 = (z > 1 - x && z > x ? 1 : 0);
  return 10.5 - (28*NC)/3. - (23*x)/2. + (29*NC*x)/6. - (52*z)/3. + 10*NC*z +  
   (2*NC*nf*z)/3. + (233*x*z)/9. - 9*NC*x*z - (16*NC*nf*x*z)/9. -  
   (5*sqrtxz3*itani1)/2. - (NC*sqrtxz3*itani1)/2. +  
   (9*sqrtxz3*x*z*itani1)/2. -  
   (3*NC*sqrtxz3*x*z*itani1)/2. +  
   (5*sqrtxz3*itani2)/2. + (NC*sqrtxz3*itani2)/2. -  
   (9*sqrtxz3*x*z*itani2)/2. +  
   (3*NC*sqrtxz3*x*z*itani2)/2. -  
   5*sqrtxz3*itani3 - NC*sqrtxz3*itani3 +  
   9*sqrtxz3*x*z*itani3 -  
   3*NC*sqrtxz3*x*z*itani3 + 4*x*Li2mx - 4*NC*x*Li2mx +  
   8*NC*x*z*Li2mx + 4*x*Li2x + 2*NC*x*Li2x + 2*x*z*Li2x +  
   4*NC*x*z*Li2x - 2*x*li2spec1 -  
   2*x*z*li2spec1 + 2*x*li2spec2 +  
   2*x*z*li2spec2 - 2*x*Li2z - NC*x*Li2z -  
   10*x*z*Li2z + 6*NC*x*z*Li2z + 4*x*li2spec19 -  
   4*x*z*li2spec20 + 4*NC*x*z*li2spec20 +  
   2*x*li2spec3 +  
   2*x*z*li2spec3 -  
   4*NC*x*z*li2spec3 -  
   2*x*li2spec4 -  
   2*x*z*li2spec4 +  
   4*NC*x*z*li2spec4 +  
   4*x*li2omxzi + 2*NC*x*li2omxzi +  
   8*x*z*li2omxzi - 4*NC*x*z*li2omxzi +  
   6*sqrtxz1*x*l2 - 4*NC*sqrtxz1*x*l2 -  
   5*sqrtxz3*atanspec1*lspec3 -  
   NC*sqrtxz3*atanspec1*lspec3 +  
   9*sqrtxz3*x*z*atanspec1*lspec3 -  
   3*NC*sqrtxz3*x*z*atanspec1*lspec3 - 2*NC*lomx +  
   2*x*lomx - 2*NC*x*lomx - 2*z*lomx -  
   (19*x*z*lomx)/3. + NC*x*z*lomx + (2*NC*nf*x*z*lomx)/3. +  
   4*lx - (NC*lx)/2. + (x*lx)/2. - (3*NC*x*lx)/2. +  
   3*sqrtxz1*x*lx - 2*NC*sqrtxz1*x*lx - (9*z*lx)/2. +  
   (9*NC*z*lx)/2. + (55*x*z*lx)/6. + (NC*x*z*lx)/2. -  
   (4*NC*nf*x*z*lx)/3. + 2*x*l2*lx + 2*x*z*l2*lx +  
   4*NC*x*z*l2*lx + 14*x*lomx*lx +  
   34*x*z*lomx*lx + 4*NC*x*z*lomx*lx +  
   4*x*lx*lopx - 4*NC*x*lx*lopx +  
   8*NC*x*z*lx*lopx - 2*NC*lomz + 2*x*lomz -  
   2*NC*x*lomz - 2*z*lomz - (7*x*z*lomz)/3. +  
   NC*x*z*lomz + (2*NC*nf*x*z*lomz)/3. -  
   4*x*lomx*lomz - 20*x*z*lomx*lomz +  
   6*x*lx*lomz - 2*NC*x*lx*lomz +  
   32*x*z*lx*lomz - 6*sqrtxz1*x*lspec1 +  
   4*NC*sqrtxz1*x*lspec1 - 4*x*l2*lspec1 -  
   4*x*z*l2*lspec1 - 4*NC*x*z*l2*lspec1 -  
   4*NC*x*z*lx*lspec1 + 5*lz - (7*NC*lz)/2. -  
   (21*x*lz)/2. - (NC*x*lz)/2. + 3*sqrtxz1*x*lz -  
   2*NC*sqrtxz1*x*lz - (9*z*lz)/2. + (5*NC*z*lz)/2. -  
   (37*x*z*lz)/2. + (3*NC*x*z*lz)/2. + 6*x*l2*lz +  
   6*x*z*l2*lz - 4*NC*x*z*l2*lz - 6*x*lomx*lz -  
   NC*x*lomx*lz - 14*x*z*lomx*lz -  
   2*NC*x*z*lomx*lz + 4*x*lx*lz + NC*x*lx*lz +  
   24*x*z*lx*lz + 6*NC*x*z*lx*lz - 6*x*lomz*lz -  
   2*NC*x*lomz*lz - 24*x*z*lomz*lz +  
   4*NC*x*z*lomz*lz - 2*x*lspec1*lz -  
   2*x*z*lspec1*lz +  
   5*sqrtxz3*atanspec2*lspec4 +  
   NC*sqrtxz3*atanspec2*lspec4 -  
   9*sqrtxz3*x*z*atanspec2*lspec4 +  
   3*NC*sqrtxz3*x*z*atanspec2*lspec4 -  
   4*x*l2*lspec2 - 4*x*z*l2*lspec2 +  
   4*NC*x*z*l2*lspec2 - 2*x*lx*lspec2 -  
   2*x*z*lx*lspec2 +  
   4*x*lspec1*lspec2 +  
   4*x*z*lspec1*lspec2 -  
   4*NC*x*z*lspec1*lspec2 -  
   4*x*lz*lspec2 - 4*x*z*lz*lspec2 +  
   4*NC*x*z*lz*lspec2 - 2*x*lx*lxpz -  
   2*x*z*lx*lxpz + 2*NC*x*z*lx*lxpz +  
   2*x*lz*lxpz + 2*x*z*lz*lxpz -  
   2*NC*x*z*lz*lxpz + 4*x*lx*lopxz +  
   4*x*lz*lopxz + 2*x*lx*lopxzi -  
   2*x*z*lx*lopxzi + 2*NC*x*z*lx*lopxzi -  
   2*x*lz*lopxzi + 2*x*z*lz*lopxzi -  
   2*NC*x*z*lz*lopxzi - (21*NCi2)/2. +  
   (23*x*NCi2)/2. + 18*z*NCi2 - 27*x*z*NCi2 +  
   (5*sqrtxz3*itani1*NCi2)/2. -  
   (9*sqrtxz3*x*z*itani1*NCi2)/2. -  
   (5*sqrtxz3*itani2*NCi2)/2. +  
   (9*sqrtxz3*x*z*itani2*NCi2)/2. +  
   5*sqrtxz3*itani3*NCi2 -  
   9*sqrtxz3*x*z*itani3*NCi2 - 4*x*Li2mx*NCi2 -  
   4*x*Li2x*NCi2 - 3*x*z*Li2x*NCi2 +  
   2*x*li2spec1*NCi2 +  
   2*x*z*li2spec1*NCi2 -  
   2*x*li2spec2*NCi2 -  
   2*x*z*li2spec2*NCi2 + 2*x*Li2z*NCi2 +  
   7*x*z*Li2z*NCi2 - 4*x*li2spec19*NCi2 +  
   4*x*z*li2spec20*NCi2 -  
   2*x*li2spec3*NCi2 -  
   2*x*z*li2spec3*NCi2 +  
   2*x*li2spec4*NCi2 +  
   2*x*z*li2spec4*NCi2 -  
   4*x*li2omxzi*NCi2 -  
   8*x*z*li2omxzi*NCi2 - 6*sqrtxz1*x*l2*NCi2 +  
   5*sqrtxz3*atanspec1*lspec3*NCi2 -  
   9*sqrtxz3*x*z*atanspec1*lspec3*NCi2 +  
   x*lomx*NCi2 + z*lomx*NCi2 +  
   6*x*z*lomx*NCi2 - 4*lx*NCi2 -  
   (9*x*lx*NCi2)/2. - 3*sqrtxz1*x*lx*NCi2 +  
   (11*z*lx*NCi2)/2. - (21*x*z*lx*NCi2)/2. -  
   2*x*l2*lx*NCi2 - 2*x*z*l2*lx*NCi2 -  
   14*x*lomx*lx*NCi2 - 22*x*z*lomx*lx*NCi2 -  
   4*x*lx*lopx*NCi2 + x*lomz*NCi2 +  
   z*lomz*NCi2 + 4*x*z*lomz*NCi2 +  
   4*x*lomx*lomz*NCi2 +  
   10*x*z*lomx*lomz*NCi2 -  
   6*x*lx*lomz*NCi2 - 19*x*z*lx*lomz*NCi2 +  
   6*sqrtxz1*x*lspec1*NCi2 +  
   4*x*l2*lspec1*NCi2 +  
   4*x*z*l2*lspec1*NCi2 - 5*lz*NCi2 +  
   (21*x*lz*NCi2)/2. - 3*sqrtxz1*x*lz*NCi2 +  
   (9*z*lz*NCi2)/2. + (25*x*z*lz*NCi2)/2. -  
   6*x*l2*lz*NCi2 - 6*x*z*l2*lz*NCi2 +  
   6*x*lomx*lz*NCi2 + 11*x*z*lomx*lz*NCi2 -  
   4*x*lx*lz*NCi2 - 19*x*z*lx*lz*NCi2 +  
   6*x*lomz*lz*NCi2 + 18*x*z*lomz*lz*NCi2 +  
   2*x*lspec1*lz*NCi2 +  
   2*x*z*lspec1*lz*NCi2 -  
   5*sqrtxz3*atanspec2*lspec4*NCi2 +  
   9*sqrtxz3*x*z*atanspec2*lspec4*NCi2 +  
   4*x*l2*lspec2*NCi2 +  
   4*x*z*l2*lspec2*NCi2 +  
   2*x*lx*lspec2*NCi2 +  
   2*x*z*lx*lspec2*NCi2 -  
   4*x*lspec1*lspec2*NCi2 -  
   4*x*z*lspec1*lspec2*NCi2 +  
   4*x*lz*lspec2*NCi2 +  
   4*x*z*lz*lspec2*NCi2 +  
   2*x*lx*lxpz*NCi2 + 2*x*z*lx*lxpz*NCi2 -  
   2*x*lz*lxpz*NCi2 - 2*x*z*lz*lxpz*NCi2 -  
   4*x*lx*lopxz*NCi2 - 4*x*lz*lopxz*NCi2 -  
   2*x*lx*lopxzi*NCi2 +  
   2*x*z*lx*lopxzi*NCi2 +  
   2*x*lz*lopxzi*NCi2 -  
   2*x*z*lz*lopxzi*NCi2 + (28*NCi)/3. -  
   (29*x*NCi)/6. - 10*z*NCi - (2*nf*z*NCi)/3. +  
   9*x*z*NCi + (16*nf*x*z*NCi)/9. +  
   (sqrtxz3*itani1*NCi)/2. +  
   (3*sqrtxz3*x*z*itani1*NCi)/2. -  
   (sqrtxz3*itani2*NCi)/2. -  
   (3*sqrtxz3*x*z*itani2*NCi)/2. +  
   sqrtxz3*itani3*NCi +  
   3*sqrtxz3*x*z*itani3*NCi + 4*x*Li2mx*NCi -  
   8*x*z*Li2mx*NCi - 2*x*Li2x*NCi -  
   4*x*z*Li2x*NCi + x*Li2z*NCi - 6*x*z*Li2z*NCi -  
   4*x*z*li2spec20*NCi +  
   4*x*z*li2spec3*NCi -  
   4*x*z*li2spec4*NCi -  
   2*x*li2omxzi*NCi +  
   4*x*z*li2omxzi*NCi + 4*sqrtxz1*x*l2*NCi +  
   sqrtxz3*atanspec1*lspec3*NCi +  
   3*sqrtxz3*x*z*atanspec1*lspec3*NCi +  
   2*lomx*NCi + 2*x*lomx*NCi -  
   x*z*lomx*NCi - (2*nf*x*z*lomx*NCi)/3. +  
   (lx*NCi)/2. + (3*x*lx*NCi)/2. +  
   2*sqrtxz1*x*lx*NCi - (9*z*lx*NCi)/2. -  
   (x*z*lx*NCi)/2. + (4*nf*x*z*lx*NCi)/3. -  
   4*x*z*l2*lx*NCi - 4*x*z*lomx*lx*NCi +  
   4*x*lx*lopx*NCi - 8*x*z*lx*lopx*NCi +  
   2*lomz*NCi + 2*x*lomz*NCi -  
   x*z*lomz*NCi - (2*nf*x*z*lomz*NCi)/3. +  
   2*x*lx*lomz*NCi -  
   4*sqrtxz1*x*lspec1*NCi +  
   4*x*z*l2*lspec1*NCi +  
   4*x*z*lx*lspec1*NCi + (7*lz*NCi)/2. +  
   (x*lz*NCi)/2. + 2*sqrtxz1*x*lz*NCi -  
   (5*z*lz*NCi)/2. - (3*x*z*lz*NCi)/2. +  
   4*x*z*l2*lz*NCi + x*lomx*lz*NCi +  
   2*x*z*lomx*lz*NCi - x*lx*lz*NCi -  
   6*x*z*lx*lz*NCi + 2*x*lomz*lz*NCi -  
   4*x*z*lomz*lz*NCi -  
   sqrtxz3*atanspec2*lspec4*NCi -  
   3*sqrtxz3*x*z*atanspec2*lspec4*NCi -  
   4*x*z*l2*lspec2*NCi +  
   4*x*z*lspec1*lspec2*NCi -  
   4*x*z*lz*lspec2*NCi -  
   2*x*z*lx*lxpz*NCi + 2*x*z*lz*lxpz*NCi -  
   2*x*z*lx*lopxzi*NCi +  
   2*x*z*lz*lopxzi*NCi - (2*z*NC2)/3. +  
   (10*x*z*NC2)/9. + x*z*Li2x*NC2 + 3*x*z*Li2z*NC2 -  
   3*x*lomx*NC2 + z*lomx*NC2 +  
   (x*z*lomx*NC2)/3. + 4*x*lx*NC2 -  
   z*lx*NC2 + (4*x*z*lx*NC2)/3. -  
   12*x*z*lomx*lx*NC2 - 3*x*lomz*NC2 +  
   z*lomz*NC2 - (5*x*z*lomz*NC2)/3. +  
   10*x*z*lomx*lomz*NC2 -  
   13*x*z*lx*lomz*NC2 + 6*x*z*lz*NC2 +  
   3*x*z*lomx*lz*NC2 - 5*x*z*lx*lz*NC2 +  
   6*x*z*lomz*lz*NC2 + (5*x*pi2)/3. -  
   (NC*x*pi2)/2. + 4*x*z*pi2 - NC*x*z*pi2 -  
   (5*x*NCi2*pi2)/3. - 2*x*z*NCi2*pi2 +  
   (x*NCi*pi2)/2. + x*z*NCi*pi2 -  
   2*x*z*NC2*pi2 - (3*lx*poly2i2)/4. -  
   (3*x*lx*poly2i2)/4. + (3*lz*poly2i2)/4. -  
   (3*x*lz*poly2i2)/4. + (3*lx*NCi2*poly2i2)/4. +  
   (3*x*lx*NCi2*poly2i2)/4. -  
   (3*lz*NCi2*poly2i2)/4. +  
   (3*x*lz*NCi2*poly2i2)/4. - poly2i/2. +  
   (3*lx*poly2i)/4. - (3*x*lx*poly2i)/2. -  
   (3*lz*poly2i)/4. - (3*x*lz*poly2i)/2. +  
   (NCi2*poly2i)/2. - (3*lx*NCi2*poly2i)/4. +  
   (3*x*lx*NCi2*poly2i)/2. +  
   (3*lz*NCi2*poly2i)/4. +  
   (3*x*lz*NCi2*poly2i)/2. -  
   (47*x*li2spec5*sqrtxz2i)/4. +  
   4*NC*x*li2spec5*sqrtxz2i +  
   18*x*z*li2spec5*sqrtxz2i -  
   16*NC*x*z*li2spec5*sqrtxz2i +  
   (47*x*li2spec6*sqrtxz2i)/4. -  
   4*NC*x*li2spec6*sqrtxz2i -  
   18*x*z*li2spec6*sqrtxz2i +  
   16*NC*x*z*li2spec6*sqrtxz2i +  
   (47*x*li2spec7*sqrtxz2i)/ 
    4. - 4*NC*x*li2spec7* 
    sqrtxz2i - 18*x*z*li2spec7* 
    sqrtxz2i + 16*NC*x*z* 
    li2spec7*sqrtxz2i -  
   (47*x*li2spec8*sqrtxz2i)/ 
    4. + 4*NC*x*li2spec8* 
    sqrtxz2i + 18*x*z*li2spec8* 
    sqrtxz2i - 16*NC*x*z* 
    li2spec8*sqrtxz2i +  
   (47*x*lx*lspec5*sqrtxz2i)/4. -  
   4*NC*x*lx*lspec5*sqrtxz2i -  
   18*x*z*lx*lspec5*sqrtxz2i +  
   16*NC*x*z*lx*lspec5*sqrtxz2i -  
   (47*x*lx*lspec6*sqrtxz2i)/4. +  
   4*NC*x*lx*lspec6*sqrtxz2i +  
   18*x*z*lx*lspec6*sqrtxz2i -  
   16*NC*x*z*lx*lspec6*sqrtxz2i +  
   (47*x*li2spec5*NCi2*sqrtxz2i)/4. -  
   18*x*z*li2spec5*NCi2*sqrtxz2i -  
   (47*x*li2spec6*NCi2*sqrtxz2i)/4. +  
   18*x*z*li2spec6*NCi2*sqrtxz2i -  
   (47*x*li2spec7*NCi2* 
      sqrtxz2i)/4. + 18*x*z* 
    li2spec7*NCi2* 
    sqrtxz2i + (47*x*li2spec8* 
      NCi2*sqrtxz2i)/4. -  
   18*x*z*li2spec8*NCi2* 
    sqrtxz2i - (47*x*lx*lspec5*NCi2* 
      sqrtxz2i)/4. + 18*x*z*lx*lspec5*NCi2* 
    sqrtxz2i + (47*x*lx*lspec6*NCi2* 
      sqrtxz2i)/4. - 18*x*z*lx*lspec6*NCi2* 
    sqrtxz2i - 4*x*li2spec5*NCi* 
    sqrtxz2i + 16*x*z*li2spec5*NCi* 
    sqrtxz2i + 4*x*li2spec6*NCi* 
    sqrtxz2i - 16*x*z*li2spec6*NCi* 
    sqrtxz2i + 4*x*li2spec7* 
    NCi*sqrtxz2i -  
   16*x*z*li2spec7*NCi* 
    sqrtxz2i - 4*x*li2spec8* 
    NCi*sqrtxz2i +  
   16*x*z*li2spec8*NCi* 
    sqrtxz2i + 4*x*lx*lspec5*NCi* 
    sqrtxz2i - 16*x*z*lx*lspec5*NCi* 
    sqrtxz2i - 4*x*lx*lspec6*NCi* 
    sqrtxz2i + 16*x*z*lx*lspec6*NCi* 
    sqrtxz2i + (3*x*li2spec5*poly2i2* 
      sqrtxz2i)/4. - (3*x*li2spec6*poly2i2* 
      sqrtxz2i)/4. - (3*x* 
      li2spec7*poly2i2* 
      sqrtxz2i)/4. + (3*x* 
      li2spec8*poly2i2* 
      sqrtxz2i)/4. - (3*x*lx*lspec5*poly2i2* 
      sqrtxz2i)/4. + (3*x*lx*lspec6*poly2i2* 
      sqrtxz2i)/4. - (3*x*li2spec5*NCi2* 
      poly2i2*sqrtxz2i)/4. +  
   (3*x*li2spec6*NCi2*poly2i2* 
      sqrtxz2i)/4. + (3*x* 
      li2spec7*NCi2* 
      poly2i2*sqrtxz2i)/4. -  
   (3*x*li2spec8*NCi2* 
      poly2i2*sqrtxz2i)/4. +  
   (3*x*lx*lspec5*NCi2*poly2i2*sqrtxz2i)/ 
    4. - (3*x*lx*lspec6*NCi2*poly2i2* 
      sqrtxz2i)/4. + (x*li2spec5*poly2i* 
      sqrtxz2i)/2. - (x*li2spec6*poly2i* 
      sqrtxz2i)/2. - (x*li2spec7*poly2i*sqrtxz2i)/2. +  
   (x*li2spec8*poly2i* 
      sqrtxz2i)/2. - (x*lx*lspec5*poly2i* 
      sqrtxz2i)/2. + (x*lx*lspec6*poly2i* 
      sqrtxz2i)/2. - (x*li2spec5*NCi2* 
      poly2i*sqrtxz2i)/2. +  
   (x*li2spec6*NCi2*poly2i*sqrtxz2i)/ 
    2. + (x*li2spec7*NCi2* 
      poly2i*sqrtxz2i)/2. -  
   (x*li2spec8*NCi2* 
      poly2i*sqrtxz2i)/2. +  
   (x*lx*lspec5*NCi2*poly2i*sqrtxz2i)/ 
    2. - (x*lx*lspec6*NCi2*poly2i* 
      sqrtxz2i)/2. - xi/2. + (10*NC*xi)/9. -  
   (3*sqrtxz3*z*itani1*xi)/2. +  
   (NC*sqrtxz3*z*itani1*xi)/2. +  
   (3*sqrtxz3*z*itani2*xi)/2. -  
   (NC*sqrtxz3*z*itani2*xi)/2. -  
   3*sqrtxz3*z*itani3*xi +  
   NC*sqrtxz3*z*itani3*xi + 4*Li2mx*xi -  
   4*z*Li2mx*xi + 4*NC*z*Li2mx*xi - 4*Li2x*xi +  
   4*z*Li2x*xi - 4*NC*z*Li2x*xi - 2*li2spec19*xi +  
   2*z*li2spec19*xi - 2*NC*z*li2spec19*xi -  
   2*li2spec20*xi + 2*z*li2spec20*xi -  
   2*NC*z*li2spec20*xi -  
   3*sqrtxz3*z*atanspec1*lspec3*xi +  
   NC*sqrtxz3*z*atanspec1*lspec3*xi +  
   (2*NC*lomx*xi)/3. + (17*lx*xi)/2. -  
   8*z*lx*xi + 8*NC*z*lx*xi -  
   4*lomx*lx*xi + 4*z*lomx*lx*xi -  
   4*NC*z*lomx*lx*xi + 4*lx*lopx*xi -  
   4*z*lx*lopx*xi + 4*NC*z*lx*lopx*xi +  
   (2*NC*lomz*xi)/3. + (lz*xi)/2. +  
   (2*NC*lz*xi)/3. + 2*lx*lz*xi -  
   2*z*lx*lz*xi + 2*NC*z*lx*lz*xi +  
   3*sqrtxz3*z*atanspec2*lspec4*xi -  
   NC*sqrtxz3*z*atanspec2*lspec4*xi -  
   2*lx*lopxz*xi + 2*z*lx*lopxz*xi -  
   2*NC*z*lx*lopxz*xi - 2*lz*lopxz*xi +  
   2*z*lz*lopxz*xi - 2*NC*z*lz*lopxz*xi -  
   2*lx*lopxzi*xi +  
   2*z*lx*lopxzi*xi -  
   2*NC*z*lx*lopxzi*xi +  
   2*lz*lopxzi*xi -  
   2*z*lz*lopxzi*xi +  
   2*NC*z*lz*lopxzi*xi + (NCi2*xi)/2. +  
   (3*sqrtxz3*z*itani1*NCi2*xi)/2. -  
   (3*sqrtxz3*z*itani2*NCi2*xi)/2. +  
   3*sqrtxz3*z*itani3*NCi2*xi -  
   4*Li2mx*NCi2*xi + 4*z*Li2mx*NCi2*xi +  
   4*Li2x*NCi2*xi - 4*z*Li2x*NCi2*xi +  
   2*li2spec19*NCi2*xi -  
   2*z*li2spec19*NCi2*xi +  
   2*li2spec20*NCi2*xi -  
   2*z*li2spec20*NCi2*xi +  
   3*sqrtxz3*z*atanspec1*lspec3*NCi2*xi -  
   (17*lx*NCi2*xi)/2. + 8*z*lx*NCi2*xi +  
   4*lomx*lx*NCi2*xi -  
   4*z*lomx*lx*NCi2*xi -  
   4*lx*lopx*NCi2*xi +  
   4*z*lx*lopx*NCi2*xi -  
   (lz*NCi2*xi)/2. - 2*lx*lz*NCi2*xi +  
   2*z*lx*lz*NCi2*xi -  
   3*sqrtxz3*z*atanspec2*lspec4*NCi2*xi +  
   2*lx*lopxz*NCi2*xi -  
   2*z*lx*lopxz*NCi2*xi +  
   2*lz*lopxz*NCi2*xi -  
   2*z*lz*lopxz*NCi2*xi +  
   2*lx*lopxzi*NCi2*xi -  
   2*z*lx*lopxzi*NCi2*xi -  
   2*lz*lopxzi*NCi2*xi +  
   2*z*lz*lopxzi*NCi2*xi -  
   (10*NCi*xi)/9. -  
   (sqrtxz3*z*itani1*NCi*xi)/2. +  
   (sqrtxz3*z*itani2*NCi*xi)/2. -  
   sqrtxz3*z*itani3*NCi*xi -  
   4*z*Li2mx*NCi*xi + 4*z*Li2x*NCi*xi +  
   2*z*li2spec19*NCi*xi +  
   2*z*li2spec20*NCi*xi -  
   sqrtxz3*z*atanspec1*lspec3*NCi*xi -  
   (2*lomx*NCi*xi)/3. - 8*z*lx*NCi*xi +  
   4*z*lomx*lx*NCi*xi -  
   4*z*lx*lopx*NCi*xi -  
   (2*lomz*NCi*xi)/3. -  
   (2*lz*NCi*xi)/3. -  
   2*z*lx*lz*NCi*xi +  
   sqrtxz3*z*atanspec2*lspec4*NCi*xi +  
   2*z*lx*lopxz*NCi*xi +  
   2*z*lz*lopxz*NCi*xi +  
   2*z*lx*lopxzi*NCi*xi -  
   2*z*lz*lopxzi*NCi*xi +  
   (2*pi2*xi)/3. - (2*z*pi2*xi)/3. +  
   (2*NC*z*pi2*xi)/3. - (2*NCi2*pi2*xi)/3. +  
   (2*z*NCi2*pi2*xi)/3. -  
   (2*z*NCi*pi2*xi)/3. +  
   (3*lx*poly2i2*xi)/4. +  
   (3*lz*poly2i2*xi)/4. -  
   (3*lx*NCi2*poly2i2*xi)/4. -  
   (3*lz*NCi2*poly2i2*xi)/4. +  
   (poly2i*xi)/2. - (5*lx*poly2i*xi)/4. -  
   (5*lz*poly2i*xi)/4. -  
   (NCi2*poly2i*xi)/2. +  
   (5*lx*NCi2*poly2i*xi)/4. +  
   (5*lz*NCi2*poly2i*xi)/4. -  
   (3*li2spec5*sqrtxz2i*xi)/8. +  
   (3*li2spec6*sqrtxz2i*xi)/8. +  
   (3*li2spec7*sqrtxz2i* 
      xi)/8. - (3*li2spec8* 
      sqrtxz2i*xi)/8. +  
   (3*lx*lspec5*sqrtxz2i*xi)/8. -  
   (3*lx*lspec6*sqrtxz2i*xi)/8. +  
   (3*li2spec5*NCi2*sqrtxz2i*xi)/8. -  
   (3*li2spec6*NCi2*sqrtxz2i*xi)/8. -  
   (3*li2spec7*NCi2* 
      sqrtxz2i*xi)/8. +  
   (3*li2spec8*NCi2* 
      sqrtxz2i*xi)/8. -  
   (3*lx*lspec5*NCi2*sqrtxz2i*xi)/8. +  
   (3*lx*lspec6*NCi2*sqrtxz2i*xi)/8. -  
   (3*li2spec5*poly2i2*sqrtxz2i*xi)/ 
    8. + (3*li2spec6*poly2i2*sqrtxz2i* 
      xi)/8. + (3*li2spec7* 
      poly2i2*sqrtxz2i*xi)/8. -  
   (3*li2spec8*poly2i2* 
      sqrtxz2i*xi)/8. +  
   (3*lx*lspec5*poly2i2*sqrtxz2i*xi)/ 
    8. - (3*lx*lspec6*poly2i2*sqrtxz2i* 
      xi)/8. + (3*li2spec5*NCi2* 
      poly2i2*sqrtxz2i*xi)/8. -  
   (3*li2spec6*NCi2*poly2i2*sqrtxz2i* 
      xi)/8. - (3*li2spec7* 
      NCi2*poly2i2*sqrtxz2i*xi)/8. +  
   (3*li2spec8*NCi2* 
      poly2i2*sqrtxz2i*xi)/8. -  
   (3*lx*lspec5*NCi2*poly2i2*sqrtxz2i* 
      xi)/8. + (3*lx*lspec6*NCi2*poly2i2* 
      sqrtxz2i*xi)/8. +  
   (3*li2spec5*poly2i*sqrtxz2i*xi)/ 
    4. - (3*li2spec6*poly2i*sqrtxz2i* 
      xi)/4. - (3*li2spec7* 
      poly2i*sqrtxz2i*xi)/4. +  
   (3*li2spec8*poly2i* 
      sqrtxz2i*xi)/4. -  
   (3*lx*lspec5*poly2i*sqrtxz2i*xi)/ 
    4. + (3*lx*lspec6*poly2i*sqrtxz2i* 
      xi)/4. - (3*li2spec5*NCi2* 
      poly2i*sqrtxz2i*xi)/4. +  
   (3*li2spec6*NCi2*poly2i*sqrtxz2i* 
      xi)/4. + (3*li2spec7* 
      NCi2*poly2i*sqrtxz2i*xi)/4. -  
   (3*li2spec8*NCi2* 
      poly2i*sqrtxz2i*xi)/4. +  
   (3*lx*lspec5*NCi2*poly2i*sqrtxz2i* 
      xi)/4. - (3*lx*lspec6*NCi2*poly2i* 
      sqrtxz2i*xi)/4. - x2/2. + (44*NC*x2)/9. +  
   (4*NC*lomx*x2)/3. + 5*lx*x2 - 6*NC*lx*x2 +  
   (4*NC*lomz*x2)/3. - 5*lz*x2 +  
   (10*NC*lz*x2)/3. - (NCi2*x2)/2. -  
   5*lx*NCi2*x2 + 5*lz*NCi2*x2 -  
   (44*NCi*x2)/9. - (4*lomx*NCi*x2)/3. +  
   6*lx*NCi*x2 - (4*lomz*NCi*x2)/3. -  
   (10*lz*NCi*x2)/3. + NC2*x2 +  
   (3*lx*poly2i2*x2)/4. -  
   (3*lz*poly2i2*x2)/4. -  
   (3*lx*NCi2*poly2i2*x2)/4. +  
   (3*lz*NCi2*poly2i2*x2)/4. +  
   (3*lx*poly2i*x2)/2. -  
   (3*lz*poly2i*x2)/2. -  
   (3*lx*NCi2*poly2i*x2)/2. +  
   (3*lz*NCi2*poly2i*x2)/2. -  
   (9*li2spec5*sqrtxz2i*x2)/2. +  
   4*NC*li2spec5*sqrtxz2i*x2 +  
   9*z*li2spec5*sqrtxz2i*x2 -  
   8*NC*z*li2spec5*sqrtxz2i*x2 +  
   (9*li2spec6*sqrtxz2i*x2)/2. -  
   4*NC*li2spec6*sqrtxz2i*x2 -  
   9*z*li2spec6*sqrtxz2i*x2 +  
   8*NC*z*li2spec6*sqrtxz2i*x2 +  
   (9*li2spec7*sqrtxz2i* 
      x2)/2. - 4*NC*li2spec7* 
    sqrtxz2i*x2 - 9*z* 
    li2spec7*sqrtxz2i*x2 
    + 8*NC*z*li2spec7*sqrtxz2i* 
    x2 - (9*li2spec8* 
      sqrtxz2i*x2)/2. +  
   4*NC*li2spec8*sqrtxz2i* 
    x2 + 9*z*li2spec8* 
    sqrtxz2i*x2 - 8*NC*z* 
    li2spec8*sqrtxz2i*x2 
    + (9*lx*lspec5*sqrtxz2i*x2)/2. -  
   4*NC*lx*lspec5*sqrtxz2i*x2 -  
   9*z*lx*lspec5*sqrtxz2i*x2 +  
   8*NC*z*lx*lspec5*sqrtxz2i*x2 -  
   (9*lx*lspec6*sqrtxz2i*x2)/2. +  
   4*NC*lx*lspec6*sqrtxz2i*x2 +  
   9*z*lx*lspec6*sqrtxz2i*x2 -  
   8*NC*z*lx*lspec6*sqrtxz2i*x2 +  
   (9*li2spec5*NCi2*sqrtxz2i*x2)/2. -  
   9*z*li2spec5*NCi2*sqrtxz2i*x2 -  
   (9*li2spec6*NCi2*sqrtxz2i*x2)/2. +  
   9*z*li2spec6*NCi2*sqrtxz2i*x2 -  
   (9*li2spec7*NCi2* 
      sqrtxz2i*x2)/2. +  
   9*z*li2spec7*NCi2* 
    sqrtxz2i*x2 + (9* 
      li2spec8*NCi2* 
      sqrtxz2i*x2)/2. -  
   9*z*li2spec8*NCi2* 
    sqrtxz2i*x2 - (9*lx*lspec5*NCi2* 
      sqrtxz2i*x2)/2. +  
   9*z*lx*lspec5*NCi2*sqrtxz2i*x2 +  
   (9*lx*lspec6*NCi2*sqrtxz2i*x2)/2. -  
   9*z*lx*lspec6*NCi2*sqrtxz2i*x2 -  
   4*li2spec5*NCi*sqrtxz2i*x2 +  
   8*z*li2spec5*NCi*sqrtxz2i*x2 +  
   4*li2spec6*NCi*sqrtxz2i*x2 -  
   8*z*li2spec6*NCi*sqrtxz2i*x2 +  
   4*li2spec7*NCi* 
    sqrtxz2i*x2 - 8*z* 
    li2spec7*NCi* 
    sqrtxz2i*x2 - 4* 
    li2spec8*NCi* 
    sqrtxz2i*x2 + 8*z* 
    li2spec8*NCi* 
    sqrtxz2i*x2 + 4*lx*lspec5*NCi* 
    sqrtxz2i*x2 - 8*z*lx*lspec5*NCi* 
    sqrtxz2i*x2 - 4*lx*lspec6*NCi* 
    sqrtxz2i*x2 + 8*z*lx*lspec6*NCi* 
    sqrtxz2i*x2 - (3*lx*poly2i2*x3)/4. -  
   (3*lz*poly2i2*x3)/4. +  
   (3*lx*NCi2*poly2i2*x3)/4. +  
   (3*lz*NCi2*poly2i2*x3)/4. -  
   (poly2i*x3)/2. - (5*lx*poly2i*x3)/4. -  
   (5*lz*poly2i*x3)/4. +  
   (NCi2*poly2i*x3)/2. +  
   (5*lx*NCi2*poly2i*x3)/4. +  
   (5*lz*NCi2*poly2i*x3)/4. +  
   (5*li2spec5*sqrtxz2i*x3)/8. -  
   (5*li2spec6*sqrtxz2i*x3)/8. -  
   (5*li2spec7*sqrtxz2i* 
      x3)/8. + (5*li2spec8* 
      sqrtxz2i*x3)/8. -  
   (5*lx*lspec5*sqrtxz2i*x3)/8. +  
   (5*lx*lspec6*sqrtxz2i*x3)/8. -  
   (5*li2spec5*NCi2*sqrtxz2i*x3)/8. +  
   (5*li2spec6*NCi2*sqrtxz2i*x3)/8. +  
   (5*li2spec7*NCi2* 
      sqrtxz2i*x3)/8. -  
   (5*li2spec8*NCi2* 
      sqrtxz2i*x3)/8. +  
   (5*lx*lspec5*NCi2*sqrtxz2i*x3)/8. -  
   (5*lx*lspec6*NCi2*sqrtxz2i*x3)/8. -  
   (li2spec5*poly2i*sqrtxz2i*x3)/4. +  
   (li2spec6*poly2i*sqrtxz2i*x3)/4. +  
   (li2spec7*poly2i* 
      sqrtxz2i*x3)/4. -  
   (li2spec8*poly2i* 
      sqrtxz2i*x3)/4. +  
   (lx*lspec5*poly2i*sqrtxz2i*x3)/4. -  
   (lx*lspec6*poly2i*sqrtxz2i*x3)/4. +  
   (li2spec5*NCi2*poly2i*sqrtxz2i* 
      x3)/4. - (li2spec6*NCi2*poly2i* 
      sqrtxz2i*x3)/4. -  
   (li2spec7*NCi2*poly2i* 
      sqrtxz2i*x3)/4. +  
   (li2spec8*NCi2*poly2i* 
      sqrtxz2i*x3)/4. -  
   (lx*lspec5*NCi2*poly2i*sqrtxz2i* 
      x3)/4. + (lx*lspec6*NCi2*poly2i* 
      sqrtxz2i*x3)/4. + (3*lx*poly2i2*x4)/4. -  
   (3*lz*poly2i2*x4)/4. -  
   (3*lx*NCi2*poly2i2*x4)/4. +  
   (3*lz*NCi2*poly2i2*x4)/4. +  
   (poly2i*x4)/2. + (7*lx*poly2i*x4)/4. -  
   (7*lz*poly2i*x4)/4. -  
   (NCi2*poly2i*x4)/2. -  
   (7*lx*NCi2*poly2i*x4)/4. +  
   (7*lz*NCi2*poly2i*x4)/4. +  
   (3*lx*poly2i2*x5)/4. +  
   (3*lz*poly2i2*x5)/4. -  
   (3*lx*NCi2*poly2i2*x5)/4. -  
   (3*lz*NCi2*poly2i2*x5)/4. -  
   (3*li2spec5*poly2i2*sqrtxz2i*x5)/ 
    4. + (3*li2spec6*poly2i2*sqrtxz2i* 
      x5)/4. + (3*li2spec7* 
      poly2i2*sqrtxz2i*x5)/4. -  
   (3*li2spec8*poly2i2* 
      sqrtxz2i*x5)/4. +  
   (3*lx*lspec5*poly2i2*sqrtxz2i*x5)/ 
    4. - (3*lx*lspec6*poly2i2*sqrtxz2i* 
      x5)/4. + (3*li2spec5*NCi2*poly2i2* 
      sqrtxz2i*x5)/4. -  
   (3*li2spec6*NCi2*poly2i2*sqrtxz2i* 
      x5)/4. - (3*li2spec7* 
      NCi2*poly2i2*sqrtxz2i*x5)/4. +  
   (3*li2spec8*NCi2* 
      poly2i2*sqrtxz2i*x5)/4. -  
   (3*lx*lspec5*NCi2*poly2i2*sqrtxz2i* 
      x5)/4. + (3*lx*lspec6*NCi2*poly2i2* 
      sqrtxz2i*x5)/4. -  
   li2spec5*poly2i*sqrtxz2i*x5 +  
   li2spec6*poly2i*sqrtxz2i*x5 +  
   li2spec7*poly2i* 
    sqrtxz2i*x5 - li2spec8*poly2i*sqrtxz2i*x5 +  
   lx*lspec5*poly2i*sqrtxz2i*x5 -  
   lx*lspec6*poly2i*sqrtxz2i*x5 +  
   li2spec5*NCi2*poly2i*sqrtxz2i* 
    x5 - li2spec6*NCi2*poly2i* 
    sqrtxz2i*x5 - li2spec7*NCi2*poly2i*sqrtxz2i* 
    x5 + li2spec8*NCi2* 
    poly2i*sqrtxz2i*x5 -  
   lx*lspec5*NCi2*poly2i*sqrtxz2i* 
    x5 + lx*lspec6*NCi2*poly2i* 
    sqrtxz2i*x5 - (3*lx*poly2i2*x6)/4. +  
   (3*lz*poly2i2*x6)/4. +  
   (3*lx*NCi2*poly2i2*x6)/4. -  
   (3*lz*NCi2*poly2i2*x6)/4. +  
   (3*li2spec5*poly2i2*sqrtxz2i*x7)/ 
    8. - (3*li2spec6*poly2i2*sqrtxz2i* 
      x7)/8. - (3*li2spec7* 
      poly2i2*sqrtxz2i*x7)/8. +  
   (3*li2spec8*poly2i2* 
      sqrtxz2i*x7)/8. -  
   (3*lx*lspec5*poly2i2*sqrtxz2i*x7)/ 
    8. + (3*lx*lspec6*poly2i2*sqrtxz2i* 
      x7)/8. - (3*li2spec5*NCi2*poly2i2* 
      sqrtxz2i*x7)/8. +  
   (3*li2spec6*NCi2*poly2i2*sqrtxz2i* 
      x7)/8. + (3*li2spec7* 
      NCi2*poly2i2*sqrtxz2i*x7)/8. -  
   (3*li2spec8*NCi2* 
      poly2i2*sqrtxz2i*x7)/8. +  
   (3*lx*lspec5*NCi2*poly2i2*sqrtxz2i* 
      x7)/8. - (3*lx*lspec6*NCi2*poly2i2* 
      sqrtxz2i*x7)/8. - 4*Li2mx*opxi +  
   4*z*Li2mx*opxi - 4*NC*z*Li2mx*opxi +  
   4*Li2x*opxi - 4*z*Li2x*opxi +  
   4*NC*z*Li2x*opxi + 2*li2spec19*opxi -  
   2*z*li2spec19*opxi + 2*NC*z*li2spec19*opxi +  
   2*li2spec20*opxi -  
   2*z*li2spec20*opxi +  
   2*NC*z*li2spec20*opxi - 8*lx*opxi +  
   8*z*lx*opxi - 8*NC*z*lx*opxi +  
   4*lomx*lx*opxi - 4*z*lomx*lx*opxi +  
   4*NC*z*lomx*lx*opxi -  
   4*lx*lopx*opxi + 4*z*lx*lopx*opxi -  
   4*NC*z*lx*lopx*opxi - 2*lx*lz*opxi +  
   2*z*lx*lz*opxi - 2*NC*z*lx*lz*opxi +  
   2*lx*lopxz*opxi -  
   2*z*lx*lopxz*opxi +  
   2*NC*z*lx*lopxz*opxi +  
   2*lz*lopxz*opxi -  
   2*z*lz*lopxz*opxi +  
   2*NC*z*lz*lopxz*opxi +  
   2*lx*lopxzi*opxi -  
   2*z*lx*lopxzi*opxi +  
   2*NC*z*lx*lopxzi*opxi -  
   2*lz*lopxzi*opxi +  
   2*z*lz*lopxzi*opxi -  
   2*NC*z*lz*lopxzi*opxi +  
   4*Li2mx*NCi2*opxi -  
   4*z*Li2mx*NCi2*opxi - 4*Li2x*NCi2*opxi +  
   4*z*Li2x*NCi2*opxi -  
   2*li2spec19*NCi2*opxi +  
   2*z*li2spec19*NCi2*opxi -  
   2*li2spec20*NCi2*opxi +  
   2*z*li2spec20*NCi2*opxi +  
   8*lx*NCi2*opxi - 8*z*lx*NCi2*opxi -  
   4*lomx*lx*NCi2*opxi +  
   4*z*lomx*lx*NCi2*opxi +  
   4*lx*lopx*NCi2*opxi -  
   4*z*lx*lopx*NCi2*opxi +  
   2*lx*lz*NCi2*opxi -  
   2*z*lx*lz*NCi2*opxi -  
   2*lx*lopxz*NCi2*opxi +  
   2*z*lx*lopxz*NCi2*opxi -  
   2*lz*lopxz*NCi2*opxi +  
   2*z*lz*lopxz*NCi2*opxi -  
   2*lx*lopxzi*NCi2*opxi +  
   2*z*lx*lopxzi*NCi2*opxi +  
   2*lz*lopxzi*NCi2*opxi -  
   2*z*lz*lopxzi*NCi2*opxi +  
   4*z*Li2mx*NCi*opxi -  
   4*z*Li2x*NCi*opxi -  
   2*z*li2spec19*NCi*opxi -  
   2*z*li2spec20*NCi*opxi +  
   8*z*lx*NCi*opxi -  
   4*z*lomx*lx*NCi*opxi +  
   4*z*lx*lopx*NCi*opxi +  
   2*z*lx*lz*NCi*opxi -  
   2*z*lx*lopxz*NCi*opxi -  
   2*z*lz*lopxz*NCi*opxi -  
   2*z*lx*lopxzi*NCi*opxi +  
   2*z*lz*lopxzi*NCi*opxi -  
   (2*pi2*opxi)/3. + (2*z*pi2*opxi)/3. -  
   (2*NC*z*pi2*opxi)/3. +  
   (2*NCi2*pi2*opxi)/3. -  
   (2*z*NCi2*pi2*opxi)/3. +  
   (2*z*NCi*pi2*opxi)/3. -  
   4*Li2mx*xi*opxi + 4*z*Li2mx*xi*opxi -  
   4*NC*z*Li2mx*xi*opxi +  
   4*Li2x*xi*opxi - 4*z*Li2x*xi*opxi +  
   4*NC*z*Li2x*xi*opxi +  
   2*li2spec19*xi*opxi -  
   2*z*li2spec19*xi*opxi +  
   2*NC*z*li2spec19*xi*opxi +  
   2*li2spec20*xi*opxi -  
   2*z*li2spec20*xi*opxi +  
   2*NC*z*li2spec20*xi*opxi -  
   8*lx*xi*opxi + 8*z*lx*xi*opxi -  
   8*NC*z*lx*xi*opxi +  
   4*lomx*lx*xi*opxi -  
   4*z*lomx*lx*xi*opxi +  
   4*NC*z*lomx*lx*xi*opxi -  
   4*lx*lopx*xi*opxi +  
   4*z*lx*lopx*xi*opxi -  
   4*NC*z*lx*lopx*xi*opxi -  
   2*lx*lz*xi*opxi +  
   2*z*lx*lz*xi*opxi -  
   2*NC*z*lx*lz*xi*opxi +  
   2*lx*lopxz*xi*opxi -  
   2*z*lx*lopxz*xi*opxi +  
   2*NC*z*lx*lopxz*xi*opxi +  
   2*lz*lopxz*xi*opxi -  
   2*z*lz*lopxz*xi*opxi +  
   2*NC*z*lz*lopxz*xi*opxi +  
   2*lx*lopxzi*xi*opxi -  
   2*z*lx*lopxzi*xi*opxi +  
   2*NC*z*lx*lopxzi*xi*opxi -  
   2*lz*lopxzi*xi*opxi +  
   2*z*lz*lopxzi*xi*opxi -  
   2*NC*z*lz*lopxzi*xi*opxi +  
   4*Li2mx*NCi2*xi*opxi -  
   4*z*Li2mx*NCi2*xi*opxi -  
   4*Li2x*NCi2*xi*opxi +  
   4*z*Li2x*NCi2*xi*opxi -  
   2*li2spec19*NCi2*xi*opxi +  
   2*z*li2spec19*NCi2*xi*opxi -  
   2*li2spec20*NCi2*xi*opxi +  
   2*z*li2spec20*NCi2*xi*opxi +  
   8*lx*NCi2*xi*opxi -  
   8*z*lx*NCi2*xi*opxi -  
   4*lomx*lx*NCi2*xi*opxi +  
   4*z*lomx*lx*NCi2*xi*opxi +  
   4*lx*lopx*NCi2*xi*opxi -  
   4*z*lx*lopx*NCi2*xi*opxi +  
   2*lx*lz*NCi2*xi*opxi -  
   2*z*lx*lz*NCi2*xi*opxi -  
   2*lx*lopxz*NCi2*xi*opxi +  
   2*z*lx*lopxz*NCi2*xi*opxi -  
   2*lz*lopxz*NCi2*xi*opxi +  
   2*z*lz*lopxz*NCi2*xi*opxi -  
   2*lx*lopxzi*NCi2*xi*opxi +  
   2*z*lx*lopxzi*NCi2*xi*opxi +  
   2*lz*lopxzi*NCi2*xi*opxi -  
   2*z*lz*lopxzi*NCi2*xi*opxi +  
   4*z*Li2mx*NCi*xi*opxi -  
   4*z*Li2x*NCi*xi*opxi -  
   2*z*li2spec19*NCi*xi*opxi -  
   2*z*li2spec20*NCi*xi*opxi +  
   8*z*lx*NCi*xi*opxi -  
   4*z*lomx*lx*NCi*xi*opxi +  
   4*z*lx*lopx*NCi*xi*opxi +  
   2*z*lx*lz*NCi*xi*opxi -  
   2*z*lx*lopxz*NCi*xi*opxi -  
   2*z*lz*lopxz*NCi*xi*opxi -  
   2*z*lx*lopxzi*NCi*xi*opxi +  
   2*z*lz*lopxzi*NCi*xi*opxi -  
   (2*pi2*xi*opxi)/3. +  
   (2*z*pi2*xi*opxi)/3. -  
   (2*NC*z*pi2*xi*opxi)/3. +  
   (2*NCi2*pi2*xi*opxi)/3. -  
   (2*z*NCi2*pi2*xi*opxi)/3. +  
   (2*z*NCi*pi2*xi*opxi)/3. -  
   4*x*omzi - 4*x*Li2x*omzi + 2*x*Li2z*omzi +  
   16*x*lx*omzi - 14*x*lomx*lx*omzi +  
   4*x*lomx*lomz*omzi -  
   6*x*lx*lomz*omzi - 8*x*lz*omzi +  
   6*x*lomx*lz*omzi - 4*x*lx*lz*omzi +  
   6*x*lomz*lz*omzi + 4*x*NCi2*omzi +  
   4*x*Li2x*NCi2*omzi -  
   2*x*Li2z*NCi2*omzi -  
   16*x*lx*NCi2*omzi +  
   14*x*lomx*lx*NCi2*omzi -  
   4*x*lomx*lomz*NCi2*omzi +  
   6*x*lx*lomz*NCi2*omzi +  
   8*x*lz*NCi2*omzi -  
   6*x*lomx*lz*NCi2*omzi +  
   4*x*lx*lz*NCi2*omzi -  
   6*x*lomz*lz*NCi2*omzi -  
   x*pi2*omzi + x*NCi2*pi2*omzi +  
   2*x*li2spec5*sqrtxz2i*omzi -  
   2*x*li2spec6*sqrtxz2i*omzi -  
   2*x*li2spec7*sqrtxz2i* 
    omzi + 2*x*li2spec8* 
    sqrtxz2i*omzi -  
   2*x*lx*lspec5*sqrtxz2i*omzi +  
   2*x*lx*lspec6*sqrtxz2i*omzi -  
   2*x*li2spec5*NCi2*sqrtxz2i* 
    omzi + 2*x*li2spec6*NCi2* 
    sqrtxz2i*omzi +  
   2*x*li2spec7*NCi2* 
    sqrtxz2i*omzi -  
   2*x*li2spec8*NCi2* 
    sqrtxz2i*omzi +  
   2*x*lx*lspec5*NCi2*sqrtxz2i*omzi -  
   2*x*lx*lspec6*NCi2*sqrtxz2i*omzi -  
   2*li2spec5*sqrtxz2i*x2*omzi +  
   2*li2spec6*sqrtxz2i*x2*omzi +  
   2*li2spec7*sqrtxz2i* 
    x2*omzi - 2*li2spec8*sqrtxz2i*x2*omzi +  
   2*lx*lspec5*sqrtxz2i*x2*omzi -  
   2*lx*lspec6*sqrtxz2i*x2*omzi +  
   2*li2spec5*NCi2*sqrtxz2i*x2* 
    omzi - 2*li2spec6*NCi2*sqrtxz2i* 
    x2*omzi - 2*li2spec7*NCi2*sqrtxz2i*x2* 
    omzi + 2*li2spec8* 
    NCi2*sqrtxz2i*x2*omzi -  
   2*lx*lspec5*NCi2*sqrtxz2i*x2* 
    omzi + 2*lx*lspec6*NCi2*sqrtxz2i* 
    x2*omzi + lx*NCi2*x3*xmzi2 -  
   lz*NCi2*x3*xmzi2 -  
   lx*NC2*x3*xmzi2 +  
   lz*NC2*x3*xmzi2 +  
   2*lx*x4*xmzi2 - 2*lz*x4*xmzi2 -  
   3*lx*NCi2*x4*xmzi2 +  
   3*lz*NCi2*x4*xmzi2 +  
   lx*NC2*x4*xmzi2 -  
   lz*NC2*x4*xmzi2 -  
   2*lx*x2*xmzi - 2*NC*lx*x2*xmzi +  
   2*lz*x2*xmzi + 2*NC*lz*x2*xmzi -  
   NCi2*x2*xmzi +  
   2*lx*NCi*x2*xmzi -  
   2*lz*NCi*x2*xmzi +  
   NC2*x2*xmzi +  
   2*lx*NC2*x2*xmzi -  
   2*lz*NC2*x2*xmzi -  
   8*lx*x3*xmzi + 2*NC*lx*x3*xmzi +  
   8*lz*x3*xmzi - 2*NC*lz*x3*xmzi +  
   NCi2*x3*xmzi +  
   9*lx*NCi2*x3*xmzi -  
   9*lz*NCi2*x3*xmzi -  
   2*lx*NCi*x3*xmzi +  
   2*lz*NCi*x3*xmzi -  
   NC2*x3*xmzi -  
   lx*NC2*x3*xmzi +  
   lz*NC2*x3*xmzi - 4*x*zi -  
   4*x*Li2mx*zi + 4*x*li2omxzi*zi +  
   2*sqrtxz1*x*l2*zi + 15*x*lx*zi +  
   sqrtxz1*x*lx*zi - 4*x*lx*lopx*zi -  
   2*sqrtxz1*x*lspec1*zi - 9*x*lz*zi +  
   sqrtxz1*x*lz*zi + 6*x*lx*lz*zi +  
   4*x*NCi2*zi + 4*x*Li2mx*NCi2*zi -  
   4*x*li2omxzi*NCi2*zi -  
   2*sqrtxz1*x*l2*NCi2*zi -  
   15*x*lx*NCi2*zi - sqrtxz1*x*lx*NCi2*zi +  
   4*x*lx*lopx*NCi2*zi +  
   2*sqrtxz1*x*lspec1*NCi2*zi +  
   9*x*lz*NCi2*zi - sqrtxz1*x*lz*NCi2*zi -  
   6*x*lx*lz*NCi2*zi - (x*pi2*zi)/3. +  
   (x*NCi2*pi2*zi)/3. +  
   2*x*li2spec5*sqrtxz2i*zi -  
   2*x*li2spec6*sqrtxz2i*zi -  
   2*x*li2spec7*sqrtxz2i* 
    zi + 2*x*li2spec8* 
    sqrtxz2i*zi -  
   2*x*lx*lspec5*sqrtxz2i*zi +  
   2*x*lx*lspec6*sqrtxz2i*zi -  
   2*x*li2spec5*NCi2*sqrtxz2i*zi +  
   2*x*li2spec6*NCi2*sqrtxz2i*zi +  
   2*x*li2spec7*NCi2* 
    sqrtxz2i*zi -  
   2*x*li2spec8*NCi2* 
    sqrtxz2i*zi +  
   2*x*lx*lspec5*NCi2*sqrtxz2i*zi -  
   2*x*lx*lspec6*NCi2*sqrtxz2i*zi +  
   2*li2spec5*sqrtxz2i*x2*zi -  
   2*li2spec6*sqrtxz2i*x2*zi -  
   2*li2spec7*sqrtxz2i* 
    x2*zi + 2*li2spec8* 
    sqrtxz2i*x2*zi -  
   2*lx*lspec5*sqrtxz2i*x2*zi +  
   2*lx*lspec6*sqrtxz2i*x2*zi -  
   2*li2spec5*NCi2*sqrtxz2i*x2* 
    zi + 2*li2spec6*NCi2*sqrtxz2i* 
    x2*zi + 2*li2spec7* 
    NCi2*sqrtxz2i*x2*zi -  
   2*li2spec8*NCi2* 
    sqrtxz2i*x2*zi +  
   2*lx*lspec5*NCi2*sqrtxz2i*x2* 
    zi - 2*lx*lspec6*NCi2*sqrtxz2i* 
    x2*zi + 4*x*omzi*zi -  
   4*x*li2omxzi*omzi*zi -  
   16*x*lx*omzi*zi + 8*x*lz*omzi*zi -  
   4*x*lx*lz*omzi*zi -  
   4*x*NCi2*omzi*zi +  
   4*x*li2omxzi*NCi2*omzi*zi +  
   16*x*lx*NCi2*omzi*zi -  
   8*x*lz*NCi2*omzi*zi +  
   4*x*lx*lz*NCi2*omzi*zi + NC*z2 -  
   (7*NC*x*z2)/2. - (5*sqrtxz3*itani1*z2)/2. +  
   (23*NC*sqrtxz3*itani1*z2)/2. +  
   (5*sqrtxz3*itani2*z2)/2. -  
   (23*NC*sqrtxz3*itani2*z2)/2. -  
   5*sqrtxz3*itani3*z2 +  
   23*NC*sqrtxz3*itani3*z2 - 8*NC*x*Li2x*z2 +  
   8*NC*x*li2spec1*z2 -  
   8*NC*x*li2spec2*z2 -  
   8*NC*x*li2spec19*z2 -  
   5*sqrtxz3*atanspec1*lspec3*z2 +  
   23*NC*sqrtxz3*atanspec1*lspec3*z2 +  
   NC*x*lomx*z2 - 2*NC*x*lx*z2 -  
   16*NC*x*l2*lx*z2 - 8*NC*x*lomx*lx*z2 +  
   NC*x*lomz*z2 + 24*NC*x*l2*lspec1*z2 +  
   8*NC*x*lx*lspec1*z2 + NC*x*lz*z2 -  
   16*NC*x*l2*lz*z2 - 4*NC*x*lx*lz*z2 +  
   8*NC*x*lspec1*lz*z2 +  
   5*sqrtxz3*atanspec2*lspec4*z2 -  
   23*NC*sqrtxz3*atanspec2*lspec4*z2 +  
   8*NC*x*l2*lspec2*z2 +  
   8*NC*x*lx*lspec2*z2 -  
   8*NC*x*lspec1*lspec2*z2 +  
   8*NC*x*lz*lspec2*z2 +  
   4*NC*x*lx*lxpz*z2 - 4*NC*x*lz*lxpz*z2 -  
   8*NC*x*lx*lopxz*z2 -  
   8*NC*x*lz*lopxz*z2 -  
   4*NC*x*lx*lopxzi*z2 +  
   4*NC*x*lz*lopxzi*z2 +  
   (5*sqrtxz3*itani1*NCi2*z2)/2. -  
   (5*sqrtxz3*itani2*NCi2*z2)/2. +  
   5*sqrtxz3*itani3*NCi2*z2 +  
   5*sqrtxz3*atanspec1*lspec3*NCi2*z2 -  
   5*sqrtxz3*atanspec2*lspec4*NCi2*z2 -  
   NCi*z2 + (7*x*NCi*z2)/2. -  
   (23*sqrtxz3*itani1*NCi*z2)/2. +  
   (23*sqrtxz3*itani2*NCi*z2)/2. -  
   23*sqrtxz3*itani3*NCi*z2 +  
   8*x*Li2x*NCi*z2 -  
   8*x*li2spec1*NCi*z2 +  
   8*x*li2spec2*NCi*z2 +  
   8*x*li2spec19*NCi*z2 -  
   23*sqrtxz3*atanspec1*lspec3*NCi*z2 -  
   x*lomx*NCi*z2 + 2*x*lx*NCi*z2 +  
   16*x*l2*lx*NCi*z2 +  
   8*x*lomx*lx*NCi*z2 -  
   x*lomz*NCi*z2 -  
   24*x*l2*lspec1*NCi*z2 -  
   8*x*lx*lspec1*NCi*z2 -  
   x*lz*NCi*z2 + 16*x*l2*lz*NCi*z2 +  
   4*x*lx*lz*NCi*z2 -  
   8*x*lspec1*lz*NCi*z2 +  
   23*sqrtxz3*atanspec2*lspec4*NCi*z2 -  
   8*x*l2*lspec2*NCi*z2 -  
   8*x*lx*lspec2*NCi*z2 +  
   8*x*lspec1*lspec2*NCi*z2 -  
   8*x*lz*lspec2*NCi*z2 -  
   4*x*lx*lxpz*NCi*z2 +  
   4*x*lz*lxpz*NCi*z2 +  
   8*x*lx*lopxz*NCi*z2 +  
   8*x*lz*lopxz*NCi*z2 +  
   4*x*lx*lopxzi*NCi*z2 -  
   4*x*lz*lopxzi*NCi*z2 +  
   (4*NC*x*pi2*z2)/3. - (4*x*NCi*pi2*z2)/3. -  
   18*x*li2spec5*sqrtxz2i*z2 +  
   16*NC*x*li2spec5*sqrtxz2i*z2 +  
   18*x*li2spec6*sqrtxz2i*z2 -  
   16*NC*x*li2spec6*sqrtxz2i*z2 +  
   18*x*li2spec7*sqrtxz2i* 
    z2 - 16*NC*x*li2spec7* 
    sqrtxz2i*z2 - 18*x* 
    li2spec8*sqrtxz2i*z2 
    + 16*NC*x*li2spec8* 
    sqrtxz2i*z2 + 18*x*lx*lspec5* 
    sqrtxz2i*z2 - 16*NC*x*lx*lspec5* 
    sqrtxz2i*z2 - 18*x*lx*lspec6* 
    sqrtxz2i*z2 + 16*NC*x*lx*lspec6* 
    sqrtxz2i*z2 + 18*x*li2spec5*NCi2* 
    sqrtxz2i*z2 - 18*x*li2spec6*NCi2* 
    sqrtxz2i*z2 - 18*x* 
    li2spec7*NCi2* 
    sqrtxz2i*z2 + 18*x* 
    li2spec8*NCi2* 
    sqrtxz2i*z2 - 18*x*lx*lspec5*NCi2* 
    sqrtxz2i*z2 + 18*x*lx*lspec6*NCi2* 
    sqrtxz2i*z2 - 16*x*li2spec5*NCi* 
    sqrtxz2i*z2 + 16*x*li2spec6*NCi* 
    sqrtxz2i*z2 + 16*x* 
    li2spec7*NCi* 
    sqrtxz2i*z2 - 16*x* 
    li2spec8*NCi* 
    sqrtxz2i*z2 + 16*x*lx*lspec5*NCi* 
    sqrtxz2i*z2 - 16*x*lx*lspec6*NCi* 
    sqrtxz2i*z2 - 8*NC*Li2mx*xi*z2 +  
   8*NC*Li2x*xi*z2 + 4*NC*li2spec19*xi*z2 +  
   4*NC*li2spec20*xi*z2 -  
   16*NC*lx*xi*z2 +  
   8*NC*lomx*lx*xi*z2 -  
   8*NC*lx*lopx*xi*z2 -  
   4*NC*lx*lz*xi*z2 +  
   4*NC*lx*lopxz*xi*z2 +  
   4*NC*lz*lopxz*xi*z2 +  
   4*NC*lx*lopxzi*xi*z2 -  
   4*NC*lz*lopxzi*xi*z2 +  
   8*Li2mx*NCi*xi*z2 -  
   8*Li2x*NCi*xi*z2 -  
   4*li2spec19*NCi*xi*z2 -  
   4*li2spec20*NCi*xi*z2 +  
   16*lx*NCi*xi*z2 -  
   8*lomx*lx*NCi*xi*z2 +  
   8*lx*lopx*NCi*xi*z2 +  
   4*lx*lz*NCi*xi*z2 -  
   4*lx*lopxz*NCi*xi*z2 -  
   4*lz*lopxz*NCi*xi*z2 -  
   4*lx*lopxzi*NCi*xi*z2 +  
   4*lz*lopxzi*NCi*xi*z2 -  
   (4*NC*pi2*xi*z2)/3. +  
   (4*NCi*pi2*xi*z2)/3. +  
   8*NC*Li2mx*opxi*z2 - 8*NC*Li2x*opxi*z2 -  
   4*NC*li2spec19*opxi*z2 -  
   4*NC*li2spec20*opxi*z2 +  
   16*NC*lx*opxi*z2 -  
   8*NC*lomx*lx*opxi*z2 +  
   8*NC*lx*lopx*opxi*z2 +  
   4*NC*lx*lz*opxi*z2 -  
   4*NC*lx*lopxz*opxi*z2 -  
   4*NC*lz*lopxz*opxi*z2 -  
   4*NC*lx*lopxzi*opxi*z2 +  
   4*NC*lz*lopxzi*opxi*z2 -  
   8*Li2mx*NCi*opxi*z2 +  
   8*Li2x*NCi*opxi*z2 +  
   4*li2spec19*NCi*opxi*z2 +  
   4*li2spec20*NCi*opxi*z2 -  
   16*lx*NCi*opxi*z2 +  
   8*lomx*lx*NCi*opxi*z2 -  
   8*lx*lopx*NCi*opxi*z2 -  
   4*lx*lz*NCi*opxi*z2 +  
   4*lx*lopxz*NCi*opxi*z2 +  
   4*lz*lopxz*NCi*opxi*z2 +  
   4*lx*lopxzi*NCi*opxi*z2 -  
   4*lz*lopxzi*NCi*opxi*z2 +  
   (4*NC*pi2*opxi*z2)/3. -  
   (4*NCi*pi2*opxi*z2)/3. +  
   8*NC*Li2mx*xi*opxi*z2 -  
   8*NC*Li2x*xi*opxi*z2 -  
   4*NC*li2spec19*xi*opxi*z2 -  
   4*NC*li2spec20*xi*opxi*z2 +  
   16*NC*lx*xi*opxi*z2 -  
   8*NC*lomx*lx*xi*opxi*z2 +  
   8*NC*lx*lopx*xi*opxi*z2 +  
   4*NC*lx*lz*xi*opxi*z2 -  
   4*NC*lx*lopxz*xi*opxi*z2 -  
   4*NC*lz*lopxz*xi*opxi*z2 -  
   4*NC*lx*lopxzi*xi*opxi*z2 +  
   4*NC*lz*lopxzi*xi*opxi*z2 -  
   8*Li2mx*NCi*xi*opxi*z2 +  
   8*Li2x*NCi*xi*opxi*z2 +  
   4*li2spec19*NCi*xi*opxi*z2 +  
   4*li2spec20*NCi*xi*opxi*z2 -  
   16*lx*NCi*xi*opxi*z2 +  
   8*lomx*lx*NCi*xi*opxi*z2 -  
   8*lx*lopx*NCi*xi*opxi*z2 -  
   4*lx*lz*NCi*xi*opxi*z2 +  
   4*lx*lopxz*NCi*xi*opxi*z2 +  
   4*lz*lopxz*NCi*xi*opxi*z2 +  
   4*lx*lopxzi*NCi*xi*opxi* 
    z2 - 4*lz*lopxzi*NCi*xi* 
    opxi*z2 + (4*NC*pi2*xi*opxi* 
      z2)/3. - (4*NCi*pi2*xi*opxi* 
      z2)/3. - 2*x*lx*xmzi2*z3 +  
   2*x*lz*xmzi2*z3 +  
   2*x*lx*NCi2*xmzi2*z3 -  
   2*x*lz*NCi2*xmzi2*z3 + 4*x*l22 +  
   4*x*z*l22 - 4*x*NCi2*l22 -  
   4*x*z*NCi2*l22 - 16*NC*x*z2*l22 +  
   16*x*NCi*z2*l22 - 4*x*lomx2 -  
   11*x*z*lomx2 + 4*x*NCi2*lomx2 +  
   6*x*z*NCi2*lomx2 + 5*x*z*NC2*lomx2 +  
   4*x*omzi*lomx2 -  
   4*x*NCi2*omzi*lomx2 - 10*x*lx2 +  
   4*NC*x*lx2 - 23*x*z*lx2 - 6*NC*x*z*lx2 +  
   10*x*NCi2*lx2 + 15*x*z*NCi2*lx2 -  
   4*x*NCi*lx2 + 6*x*z*NCi*lx2 +  
   8*x*z*NC2*lx2 - 3*xi*lx2 +  
   3*z*xi*lx2 - 3*NC*z*xi*lx2 +  
   3*NCi2*xi*lx2 -  
   3*z*NCi2*xi*lx2 +  
   3*z*NCi*xi*lx2 + 3*opxi*lx2 -  
   3*z*opxi*lx2 + 3*NC*z*opxi*lx2 -  
   3*NCi2*opxi*lx2 +  
   3*z*NCi2*opxi*lx2 -  
   3*z*NCi*opxi*lx2 +  
   3*xi*opxi*lx2 -  
   3*z*xi*opxi*lx2 +  
   3*NC*z*xi*opxi*lx2 -  
   3*NCi2*xi*opxi*lx2 +  
   3*z*NCi2*xi*opxi*lx2 -  
   3*z*NCi*xi*opxi*lx2 +  
   2*x*omzi*lx2 -  
   2*x*NCi2*omzi*lx2 - 5*x*zi*lx2 +  
   5*x*NCi2*zi*lx2 +  
   6*x*omzi*zi*lx2 -  
   6*x*NCi2*omzi*zi*lx2 +  
   4*NC*x*z2*lx2 - 4*x*NCi*z2*lx2 +  
   6*NC*xi*z2*lx2 -  
   6*NCi*xi*z2*lx2 -  
   6*NC*opxi*z2*lx2 +  
   6*NCi*opxi*z2*lx2 -  
   6*NC*xi*opxi*z2*lx2 +  
   6*NCi*xi*opxi*z2*lx2 -  
   x*lomz2 - 11*x*z*lomz2 +  
   x*NCi2*lomz2 + 6*x*z*NCi2*lomz2 +  
   5*x*z*NC2*lomz2 + x*omzi*lomz2 -  
   x*NCi2*omzi*lomz2 +  
   4*NC*x*z*lspec1_2 -  
   4*x*z*NCi*lspec1_2 -  
   8*NC*x*z2*lspec1_2 +  
   8*x*NCi*z2*lspec1_2 - 2*x*lz2 -  
   x*z*lz2 - 6*NC*x*z*lz2 +  
   2*x*NCi2*lz2 + x*z*NCi2*lz2 +  
   6*x*z*NCi*lz2 + xi*lz2 -  
   z*xi*lz2 + NC*z*xi*lz2 -  
   NCi2*xi*lz2 +  
   z*NCi2*xi*lz2 -  
   z*NCi*xi*lz2 - opxi*lz2 +  
   z*opxi*lz2 - NC*z*opxi*lz2 +  
   NCi2*opxi*lz2 -  
   z*NCi2*opxi*lz2 +  
   z*NCi*opxi*lz2 -  
   xi*opxi*lz2 +  
   z*xi*opxi*lz2 -  
   NC*z*xi*opxi*lz2 +  
   NCi2*xi*opxi*lz2 -  
   z*NCi2*xi*opxi*lz2 +  
   z*NCi*xi*opxi*lz2 +  
   x*omzi*lz2 - x*NCi2*omzi*lz2 +  
   4*NC*x*z2*lz2 - 4*x*NCi*z2*lz2 -  
   2*NC*xi*z2*lz2 +  
   2*NCi*xi*z2*lz2 +  
   2*NC*opxi*z2*lz2 -  
   2*NCi*opxi*z2*lz2 +  
   2*NC*xi*opxi*z2*lz2 -  
   2*NCi*xi*opxi*z2*lz2 +  
   (2*x*z*Li2z - 2*x*z*li2spec13 +  
      2*x*z*li2spec15 - 4*x*lomx +  
      4*x*z*lomx + 6*x*lx - 6*x*z*lx -  
      12*x*z*lomx*lx - 4*x*lomz + 4*x*z*lomz +  
      10*x*z*lomx*lomz - 12*x*z*lx*lomz +  
      2*x*z*lx*lomxmz - 2*x*z*lomz*lomxmz -  
      2*x*lz + 2*x*z*lz + 4*x*z*lomx*lz -  
      6*x*z*lx*lz + 4*x*z*lomz*lz -  
      2*x*z*Li2z*NC2 + 2*x*z*li2spec13*NC2 -  
      2*x*z*li2spec15*NC2 +  
      4*x*lomx*NC2 - 4*x*z*lomx*NC2 -  
      6*x*lx*NC2 + 6*x*z*lx*NC2 +  
      12*x*z*lomx*lx*NC2 + 4*x*lomz*NC2 -  
      4*x*z*lomz*NC2 - 10*x*z*lomx*lomz*NC2 +  
      12*x*z*lx*lomz*NC2 -  
      2*x*z*lx*lomxmz*NC2 +  
      2*x*z*lomz*lomxmz*NC2 + 2*x*lz*NC2 -  
      2*x*z*lz*NC2 - 4*x*z*lomx*lz*NC2 +  
      6*x*z*lx*lz*NC2 - 4*x*z*lomz*lz*NC2 -  
      (4*x*z*pi2)/3. + (4*x*z*NC2*pi2)/3. +  
      4*x*z*lomx2 - 4*x*z*NC2*lomx2 +  
      7*x*z*lx2 - 7*x*z*NC2*lx2 +  
      5*x*z*lomz2 - 5*x*z*NC2*lomz2 +  
      x*z*lz2 - x*z*NC2*lz2)*Tr1 +  
   (2*x*z*Li2z + 2*x*z*li2spec16 -  
      2*x*z*li2spec17 - 4*x*lomx +  
      4*x*z*lomx + 6*x*lx - 6*x*z*lx -  
      10*x*z*lomx*lx - 4*x*lomz + 4*x*z*lomz +  
      8*x*z*lomx*lomz - 10*x*z*lx*lomz - 2*x*lz +  
      2*x*z*lz + 4*x*z*lomx*lz - 8*x*z*lx*lz +  
      6*x*z*lomz*lz + 2*x*z*lx*lmopxpz -  
      2*x*z*lomz*lmopxpz - 2*x*z*Li2z*NC2 -  
      2*x*z*li2spec16*NC2 +  
      2*x*z*li2spec17*NC2 +  
      4*x*lomx*NC2 - 4*x*z*lomx*NC2 -  
      6*x*lx*NC2 + 6*x*z*lx*NC2 +  
      10*x*z*lomx*lx*NC2 + 4*x*lomz*NC2 -  
      4*x*z*lomz*NC2 - 8*x*z*lomx*lomz*NC2 +  
      10*x*z*lx*lomz*NC2 + 2*x*lz*NC2 -  
      2*x*z*lz*NC2 - 4*x*z*lomx*lz*NC2 +  
      8*x*z*lx*lz*NC2 - 6*x*z*lomz*lz*NC2 -  
      2*x*z*lx*lmopxpz*NC2 +  
      2*x*z*lomz*lmopxpz*NC2 - (4*x*z*pi2)/3. +  
      (4*x*z*NC2*pi2)/3. + 4*x*z*lomx2 -  
      4*x*z*NC2*lomx2 + 6*x*z*lx2 -  
      6*x*z*NC2*lx2 + 4*x*z*lomz2 -  
      4*x*z*NC2*lomz2 + x*z*lz2 -  
      x*z*NC2*lz2)*Tr2 +  
   (-2*x*Li2z - 2*x*z*Li2z + 2*x*z*li2spec13 -  
      2*x*z*li2spec14 +  
      2*x*li2spec11 -  
      2*x*li2spec12 -  
      2*x*z*li2spec12 + 6*x*z*lomx -  
      12*x*z*lx - 10*x*lomx*lx - 18*x*z*lomx*lx +  
      6*x*z*lomz + 4*x*lomx*lomz +  
      10*x*z*lomx*lomz - 8*x*lx*lomz -  
      18*x*z*lx*lomz + 2*x*z*lx*lomxmz -  
      2*x*z*lomz*lomxmz + 2*x*lx*lxmz +  
      4*x*z*lx*lxmz + 8*x*z*lz + 8*x*lomx*lz +  
      12*x*z*lomx*lz - 10*x*lx*lz - 20*x*z*lx*lz +  
      4*x*lomz*lz + 12*x*z*lomz*lz -  
      2*x*lxmz*lz - 4*x*z*lxmz*lz +  
      2*x*Li2z*NCi2 + 2*x*z*Li2z*NCi2 -  
      2*x*z*li2spec13*NCi2 +  
      2*x*z*li2spec14*NCi2 -  
      2*x*li2spec11*NCi2 +  
      2*x*li2spec12*NCi2 +  
      2*x*z*li2spec12*NCi2 -  
      6*x*z*lomx*NCi2 + 12*x*z*lx*NCi2 +  
      10*x*lomx*lx*NCi2 +  
      18*x*z*lomx*lx*NCi2 - 6*x*z*lomz*NCi2 -  
      4*x*lomx*lomz*NCi2 -  
      10*x*z*lomx*lomz*NCi2 +  
      8*x*lx*lomz*NCi2 +  
      18*x*z*lx*lomz*NCi2 -  
      2*x*z*lx*lomxmz*NCi2 +  
      2*x*z*lomz*lomxmz*NCi2 -  
      2*x*lx*lxmz*NCi2 - 4*x*z*lx*lxmz*NCi2 -  
      8*x*z*lz*NCi2 - 8*x*lomx*lz*NCi2 -  
      12*x*z*lomx*lz*NCi2 + 10*x*lx*lz*NCi2 +  
      20*x*z*lx*lz*NCi2 - 4*x*lomz*lz*NCi2 -  
      12*x*z*lomz*lz*NCi2 +  
      2*x*lxmz*lz*NCi2 + 4*x*z*lxmz*lz*NCi2 -  
      x*pi2 - (4*x*z*pi2)/3. + x*NCi2*pi2 +  
      (4*x*z*NCi2*pi2)/3. + 2*x*Li2z*omzi -  
      2*x*li2spec11*omzi +  
      2*x*li2spec12*omzi +  
      10*x*lomx*lx*omzi -  
      4*x*lomx*lomz*omzi +  
      8*x*lx*lomz*omzi -  
      2*x*lx*lxmz*omzi -  
      8*x*lomx*lz*omzi +  
      10*x*lx*lz*omzi -  
      4*x*lomz*lz*omzi +  
      2*x*lxmz*lz*omzi -  
      2*x*Li2z*NCi2*omzi +  
      2*x*li2spec11*NCi2*omzi -  
      2*x*li2spec12*NCi2*omzi -  
      10*x*lomx*lx*NCi2*omzi +  
      4*x*lomx*lomz*NCi2*omzi -  
      8*x*lx*lomz*NCi2*omzi +  
      2*x*lx*lxmz*NCi2*omzi +  
      8*x*lomx*lz*NCi2*omzi -  
      10*x*lx*lz*NCi2*omzi +  
      4*x*lomz*lz*NCi2*omzi -  
      2*x*lxmz*lz*NCi2*omzi +  
      x*pi2*omzi - x*NCi2*pi2*omzi +  
      4*x*lomx2 + 5*x*z*lomx2 -  
      4*x*NCi2*lomx2 - 5*x*z*NCi2*lomx2 -  
      4*x*omzi*lomx2 +  
      4*x*NCi2*omzi*lomx2 + 6*x*lx2 +  
      12*x*z*lx2 - 6*x*NCi2*lx2 -  
      12*x*z*NCi2*lx2 - 6*x*omzi*lx2 +  
      6*x*NCi2*omzi*lx2 + x*lomz2 +  
      4*x*z*lomz2 - x*NCi2*lomz2 -  
      4*x*z*NCi2*lomz2 -  
      x*omzi*lomz2 +  
      x*NCi2*omzi*lomz2 + 4*x*lz2 +  
      8*x*z*lz2 - 4*x*NCi2*lz2 -  
      8*x*z*NCi2*lz2 - 4*x*omzi*lz2 +  
      4*x*NCi2*omzi*lz2)*Tu1 +  
   (-2*x*Li2z - 2*x*z*Li2z + 2*x*li2spec11 -  
      2*x*li2spec12 -  
      2*x*z*li2spec12 -  
      2*x*z*li2spec16 +  
      2*x*z*li2spec18 + 6*x*z*lomx -  
      12*x*z*lx - 10*x*lomx*lx - 16*x*z*lomx*lx +  
      6*x*z*lomz + 4*x*lomx*lomz +  
      8*x*z*lomx*lomz - 8*x*lx*lomz -  
      20*x*z*lx*lomz + 2*x*lx*lxmz +  
      4*x*z*lx*lxmz + 8*x*z*lz + 8*x*lomx*lz +  
      12*x*z*lomx*lz - 10*x*lx*lz - 22*x*z*lx*lz +  
      4*x*lomz*lz + 14*x*z*lomz*lz -  
      2*x*lxmz*lz - 4*x*z*lxmz*lz +  
      2*x*z*lx*lmopxpz - 2*x*z*lomz*lmopxpz +  
      2*x*Li2z*NCi2 + 2*x*z*Li2z*NCi2 -  
      2*x*li2spec11*NCi2 +  
      2*x*li2spec12*NCi2 +  
      2*x*z*li2spec12*NCi2 +  
      2*x*z*li2spec16*NCi2 -  
      2*x*z*li2spec18*NCi2 -  
      6*x*z*lomx*NCi2 + 12*x*z*lx*NCi2 +  
      10*x*lomx*lx*NCi2 +  
      16*x*z*lomx*lx*NCi2 - 6*x*z*lomz*NCi2 -  
      4*x*lomx*lomz*NCi2 -  
      8*x*z*lomx*lomz*NCi2 +  
      8*x*lx*lomz*NCi2 +  
      20*x*z*lx*lomz*NCi2 -  
      2*x*lx*lxmz*NCi2 - 4*x*z*lx*lxmz*NCi2 -  
      8*x*z*lz*NCi2 - 8*x*lomx*lz*NCi2 -  
      12*x*z*lomx*lz*NCi2 + 10*x*lx*lz*NCi2 +  
      22*x*z*lx*lz*NCi2 - 4*x*lomz*lz*NCi2 -  
      14*x*z*lomz*lz*NCi2 +  
      2*x*lxmz*lz*NCi2 + 4*x*z*lxmz*lz*NCi2 -  
      2*x*z*lx*lmopxpz*NCi2 +  
      2*x*z*lomz*lmopxpz*NCi2 - x*pi2 -  
      (4*x*z*pi2)/3. + x*NCi2*pi2 +  
      (4*x*z*NCi2*pi2)/3. + 2*x*Li2z*omzi -  
      2*x*li2spec11*omzi +  
      2*x*li2spec12*omzi +  
      10*x*lomx*lx*omzi -  
      4*x*lomx*lomz*omzi +  
      8*x*lx*lomz*omzi -  
      2*x*lx*lxmz*omzi -  
      8*x*lomx*lz*omzi +  
      10*x*lx*lz*omzi -  
      4*x*lomz*lz*omzi +  
      2*x*lxmz*lz*omzi -  
      2*x*Li2z*NCi2*omzi +  
      2*x*li2spec11*NCi2*omzi -  
      2*x*li2spec12*NCi2*omzi -  
      10*x*lomx*lx*NCi2*omzi +  
      4*x*lomx*lomz*NCi2*omzi -  
      8*x*lx*lomz*NCi2*omzi +  
      2*x*lx*lxmz*NCi2*omzi +  
      8*x*lomx*lz*NCi2*omzi -  
      10*x*lx*lz*NCi2*omzi +  
      4*x*lomz*lz*NCi2*omzi -  
      2*x*lxmz*lz*NCi2*omzi +  
      x*pi2*omzi - x*NCi2*pi2*omzi +  
      4*x*lomx2 + 5*x*z*lomx2 -  
      4*x*NCi2*lomx2 - 5*x*z*NCi2*lomx2 -  
      4*x*omzi*lomx2 +  
      4*x*NCi2*omzi*lomx2 + 6*x*lx2 +  
      13*x*z*lx2 - 6*x*NCi2*lx2 -  
      13*x*z*NCi2*lx2 - 6*x*omzi*lx2 +  
      6*x*NCi2*omzi*lx2 + x*lomz2 +  
      5*x*z*lomz2 - x*NCi2*lomz2 -  
      5*x*z*NCi2*lomz2 -  
      x*omzi*lomz2 +  
      x*NCi2*omzi*lomz2 + 4*x*lz2 +  
      8*x*z*lz2 - 4*x*NCi2*lz2 -  
      8*x*z*NCi2*lz2 - 4*x*omzi*lz2 +  
      4*x*NCi2*omzi*lz2)*Tu2 +  
   (-2*x*Li2z - 2*x*z*Li2z - 2*x*li2spec9 +  
      2*x*z*li2spec13 +  
      2*x*li2spec10 +  
      2*x*z*li2spec10 +  
      2*x*z*li2spec18 + 6*x*z*lomx -  
      12*x*z*lx - 12*x*lomx*lx - 18*x*z*lomx*lx +  
      6*x*z*lomz + 4*x*lomx*lomz +  
      6*x*z*lomx*lomz - 6*x*lx*lomz -  
      18*x*z*lx*lomz + 2*x*z*lx*lomxmz -  
      2*x*z*lomz*lomxmz + 8*x*z*lz +  
      10*x*lomx*lz + 12*x*z*lomx*lz -  
      12*x*lx*lz - 24*x*z*lx*lz + 2*x*lomz*lz +  
      12*x*z*lomz*lz + 2*x*lx*lmxpz +  
      4*x*z*lx*lmxpz - 2*x*lz*lmxpz -  
      4*x*z*lz*lmxpz + 2*x*Li2z*NCi2 +  
      2*x*z*Li2z*NCi2 + 2*x*li2spec9*NCi2 -  
      2*x*z*li2spec13*NCi2 -  
      2*x*li2spec10*NCi2 -  
      2*x*z*li2spec10*NCi2 -  
      2*x*z*li2spec18*NCi2 -  
      6*x*z*lomx*NCi2 + 12*x*z*lx*NCi2 +  
      12*x*lomx*lx*NCi2 +  
      18*x*z*lomx*lx*NCi2 - 6*x*z*lomz*NCi2 -  
      4*x*lomx*lomz*NCi2 -  
      6*x*z*lomx*lomz*NCi2 +  
      6*x*lx*lomz*NCi2 +  
      18*x*z*lx*lomz*NCi2 -  
      2*x*z*lx*lomxmz*NCi2 +  
      2*x*z*lomz*lomxmz*NCi2 - 8*x*z*lz*NCi2 -  
      10*x*lomx*lz*NCi2 -  
      12*x*z*lomx*lz*NCi2 + 12*x*lx*lz*NCi2 +  
      24*x*z*lx*lz*NCi2 - 2*x*lomz*lz*NCi2 -  
      12*x*z*lomz*lz*NCi2 -  
      2*x*lx*lmxpz*NCi2 -  
      4*x*z*lx*lmxpz*NCi2 +  
      2*x*lz*lmxpz*NCi2 +  
      4*x*z*lz*lmxpz*NCi2 - x*pi2 -  
      (8*x*z*pi2)/3. + x*NCi2*pi2 +  
      (8*x*z*NCi2*pi2)/3. + 2*x*Li2z*omzi +  
      2*x*li2spec9*omzi -  
      2*x*li2spec10*omzi +  
      12*x*lomx*lx*omzi -  
      4*x*lomx*lomz*omzi +  
      6*x*lx*lomz*omzi -  
      10*x*lomx*lz*omzi +  
      12*x*lx*lz*omzi -  
      2*x*lomz*lz*omzi -  
      2*x*lx*lmxpz*omzi +  
      2*x*lz*lmxpz*omzi -  
      2*x*Li2z*NCi2*omzi -  
      2*x*li2spec9*NCi2*omzi +  
      2*x*li2spec10*NCi2*omzi -  
      12*x*lomx*lx*NCi2*omzi +  
      4*x*lomx*lomz*NCi2*omzi -  
      6*x*lx*lomz*NCi2*omzi +  
      10*x*lomx*lz*NCi2*omzi -  
      12*x*lx*lz*NCi2*omzi +  
      2*x*lomz*lz*NCi2*omzi +  
      2*x*lx*lmxpz*NCi2*omzi -  
      2*x*lz*lmxpz*NCi2*omzi +  
      x*pi2*omzi - x*NCi2*pi2*omzi +  
      4*x*lomx2 + 7*x*z*lomx2 -  
      4*x*NCi2*lomx2 - 7*x*z*NCi2*lomx2 -  
      4*x*omzi*lomx2 +  
      4*x*NCi2*omzi*lomx2 + 7*x*lx2 +  
      14*x*z*lx2 - 7*x*NCi2*lx2 -  
      14*x*z*NCi2*lx2 - 7*x*omzi*lx2 +  
      7*x*NCi2*omzi*lx2 + x*lomz2 +  
      6*x*z*lomz2 - x*NCi2*lomz2 -  
      6*x*z*NCi2*lomz2 -  
      x*omzi*lomz2 +  
      x*NCi2*omzi*lomz2 + 5*x*lz2 +  
      10*x*z*lz2 - 5*x*NCi2*lz2 -  
      10*x*z*NCi2*lz2 - 5*x*omzi*lz2 +  
      5*x*NCi2*omzi*lz2)*Tu3 +  
   (-2*x*Li2z - 2*x*z*Li2z - 2*x*li2spec9 -  
      2*x*z*li2spec14 -  
      2*x*z*li2spec16 +  
      2*x*li2spec10 +  
      2*x*z*li2spec10 + 6*x*z*lomx -  
      12*x*z*lx - 12*x*lomx*lx - 20*x*z*lomx*lx +  
      6*x*z*lomz + 4*x*lomx*lomz +  
      8*x*z*lomx*lomz - 6*x*lx*lomz -  
      16*x*z*lx*lomz + 8*x*z*lz + 10*x*lomx*lz +  
      16*x*z*lomx*lz - 12*x*lx*lz - 22*x*z*lx*lz +  
      2*x*lomz*lz + 10*x*z*lomz*lz +  
      2*x*lx*lmxpz + 4*x*z*lx*lmxpz -  
      2*x*lz*lmxpz - 4*x*z*lz*lmxpz +  
      2*x*z*lx*lmopxpz - 2*x*z*lomz*lmopxpz +  
      2*x*Li2z*NCi2 + 2*x*z*Li2z*NCi2 +  
      2*x*li2spec9*NCi2 +  
      2*x*z*li2spec14*NCi2 +  
      2*x*z*li2spec16*NCi2 -  
      2*x*li2spec10*NCi2 -  
      2*x*z*li2spec10*NCi2 -  
      6*x*z*lomx*NCi2 + 12*x*z*lx*NCi2 +  
      12*x*lomx*lx*NCi2 +  
      20*x*z*lomx*lx*NCi2 - 6*x*z*lomz*NCi2 -  
      4*x*lomx*lomz*NCi2 -  
      8*x*z*lomx*lomz*NCi2 +  
      6*x*lx*lomz*NCi2 +  
      16*x*z*lx*lomz*NCi2 - 8*x*z*lz*NCi2 -  
      10*x*lomx*lz*NCi2 -  
      16*x*z*lomx*lz*NCi2 + 12*x*lx*lz*NCi2 +  
      22*x*z*lx*lz*NCi2 - 2*x*lomz*lz*NCi2 -  
      10*x*z*lomz*lz*NCi2 -  
      2*x*lx*lmxpz*NCi2 -  
      4*x*z*lx*lmxpz*NCi2 +  
      2*x*lz*lmxpz*NCi2 +  
      4*x*z*lz*lmxpz*NCi2 -  
      2*x*z*lx*lmopxpz*NCi2 +  
      2*x*z*lomz*lmopxpz*NCi2 - x*pi2 -  
      (4*x*z*pi2)/3. + x*NCi2*pi2 +  
      (4*x*z*NCi2*pi2)/3. + 2*x*Li2z*omzi +  
      2*x*li2spec9*omzi -  
      2*x*li2spec10*omzi +  
      12*x*lomx*lx*omzi -  
      4*x*lomx*lomz*omzi +  
      6*x*lx*lomz*omzi -  
      10*x*lomx*lz*omzi +  
      12*x*lx*lz*omzi -  
      2*x*lomz*lz*omzi -  
      2*x*lx*lmxpz*omzi +  
      2*x*lz*lmxpz*omzi -  
      2*x*Li2z*NCi2*omzi -  
      2*x*li2spec9*NCi2*omzi +  
      2*x*li2spec10*NCi2*omzi -  
      12*x*lomx*lx*NCi2*omzi +  
      4*x*lomx*lomz*NCi2*omzi -  
      6*x*lx*lomz*NCi2*omzi +  
      10*x*lomx*lz*NCi2*omzi -  
      12*x*lx*lz*NCi2*omzi +  
      2*x*lomz*lz*NCi2*omzi +  
      2*x*lx*lmxpz*NCi2*omzi -  
      2*x*lz*lmxpz*NCi2*omzi +  
      x*pi2*omzi - x*NCi2*pi2*omzi +  
      4*x*lomx2 + 5*x*z*lomx2 -  
      4*x*NCi2*lomx2 - 5*x*z*NCi2*lomx2 -  
      4*x*omzi*lomx2 +  
      4*x*NCi2*omzi*lomx2 + 7*x*lx2 +  
      13*x*z*lx2 - 7*x*NCi2*lx2 -  
      13*x*z*NCi2*lx2 - 7*x*omzi*lx2 +  
      7*x*NCi2*omzi*lx2 + x*lomz2 +  
      5*x*z*lomz2 - x*NCi2*lomz2 -  
      5*x*z*NCi2*lomz2 -  
      x*omzi*lomz2 +  
      x*NCi2*omzi*lomz2 + 5*x*lz2 +  
      8*x*z*lz2 - 5*x*NCi2*lz2 -  
      8*x*z*NCi2*lz2 - 5*x*omzi*lz2 +  
      5*x*NCi2*omzi*lz2)*Tu4;
};
double C2LQ2QPS_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double sqrtxz3  = sqrt(x/z);
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double x2 = x * x;
  const double xi  = 1. / x;
  const double xi2 = xi * xi;
  const double z2 = z * z;
  const double zi  = 1. / z;
  const double zi2 = zi * zi;
  const double lx  = log(x);
  const double lz  = log(z);
  const double lspec3   = log(sqrtxz3);
  const double lspec4   = log(sqrtxz3*z);
  const double atanspec1 = atan(sqrtxz3);
  const double atanspec2 = atan(sqrtxz3*z);
  const double itani1 = InvTanInt(-sqrtxz3);
  const double itani2 = InvTanInt(sqrtxz3);
  const double itani3 = InvTanInt(sqrtxz3*z);
  return (3*NC)/4. + (NC*x)/4. + (5*NC*z)/8. - (5*NC*x*z)/8. -  
   (NC*sqrtxz3*itani1)/4. -  
   (9*NC*sqrtxz3*x*z*itani1)/8. +  
   (NC*sqrtxz3*itani2)/4. +  
   (9*NC*sqrtxz3*x*z*itani2)/8. -  
   (NC*sqrtxz3*itani3)/2. -  
   (9*NC*sqrtxz3*x*z*itani3)/4. -  
   (NC*sqrtxz3*atanspec1*lspec3)/2. -  
   (9*NC*sqrtxz3*x*z*atanspec1*lspec3)/4. - (3*NC*lx)/8. +  
   (9*NC*x*lx)/8. + (5*NC*z*lx)/16. + (5*NC*x*z*lx)/16. -  
   (13*NC*lz)/8. + (9*NC*x*lz)/8. + (5*NC*z*lz)/16. -  
   (5*NC*x*z*lz)/16. - 2*NC*x*lx*lz +  
   (NC*sqrtxz3*atanspec2*lspec4)/2. +  
   (9*NC*sqrtxz3*x*z*atanspec2*lspec4)/4. -  
   (3*NCi)/4. - (x*NCi)/4. - (5*z*NCi)/8. +  
   (5*x*z*NCi)/8. + (sqrtxz3*itani1*NCi)/4. +  
   (9*sqrtxz3*x*z*itani1*NCi)/8. -  
   (sqrtxz3*itani2*NCi)/4. -  
   (9*sqrtxz3*x*z*itani2*NCi)/8. +  
   (sqrtxz3*itani3*NCi)/2. +  
   (9*sqrtxz3*x*z*itani3*NCi)/4. +  
   (sqrtxz3*atanspec1*lspec3*NCi)/2. +  
   (9*sqrtxz3*x*z*atanspec1*lspec3*NCi)/4. +  
   (3*lx*NCi)/8. - (9*x*lx*NCi)/8. -  
   (5*z*lx*NCi)/16. - (5*x*z*lx*NCi)/16. +  
   (13*lz*NCi)/8. - (9*x*lz*NCi)/8. -  
   (5*z*lz*NCi)/16. + (5*x*z*lz*NCi)/16. +  
   2*x*lx*lz*NCi -  
   (sqrtxz3*atanspec2*lspec4*NCi)/2. -  
   (9*sqrtxz3*x*z*atanspec2*lspec4*NCi)/4. +  
   (3*NC*sqrtxz3*itani1*xi2)/16. -  
   (3*NC*sqrtxz3*itani2*xi2)/16. +  
   (3*NC*sqrtxz3*itani3*xi2)/8. +  
   (3*NC*sqrtxz3*atanspec1*lspec3*xi2)/8. -  
   (3*NC*sqrtxz3*atanspec2*lspec4*xi2)/8. -  
   (3*sqrtxz3*itani1*NCi*xi2)/16. +  
   (3*sqrtxz3*itani2*NCi*xi2)/16. -  
   (3*sqrtxz3*itani3*NCi*xi2)/8. -  
   (3*sqrtxz3*atanspec1*lspec3*NCi*xi2)/8. +  
   (3*sqrtxz3*atanspec2*lspec4*NCi*xi2)/8. -  
   (3*NC*xi)/8. + (3*NC*sqrtxz3*z*itani1*xi)/8. -  
   (3*NC*sqrtxz3*z*itani2*xi)/8. +  
   (3*NC*sqrtxz3*z*itani3*xi)/4. +  
   (3*NC*sqrtxz3*z*atanspec1*lspec3*xi)/4. +  
   (3*NC*lx*xi)/16. + (3*NC*lz*xi)/16. -  
   (3*NC*sqrtxz3*z*atanspec2*lspec4*xi)/4. +  
   (3*NCi*xi)/8. -  
   (3*sqrtxz3*z*itani1*NCi*xi)/8. +  
   (3*sqrtxz3*z*itani2*NCi*xi)/8. -  
   (3*sqrtxz3*z*itani3*NCi*xi)/4. -  
   (3*sqrtxz3*z*atanspec1*lspec3*NCi*xi)/4. -  
   (3*lx*NCi*xi)/16. -  
   (3*lz*NCi*xi)/16. +  
   (3*sqrtxz3*z*atanspec2*lspec4*NCi*xi)/4. -  
   (5*NC*x2)/8. - (5*NC*sqrtxz3*itani1*x2)/16. +  
   (5*NC*sqrtxz3*itani2*x2)/16. -  
   (5*NC*sqrtxz3*itani3*x2)/8. -  
   (5*NC*sqrtxz3*atanspec1*lspec3*x2)/8. -  
   (5*NC*lx*x2)/16. + (5*NC*lz*x2)/16. +  
   (5*NC*sqrtxz3*atanspec2*lspec4*x2)/8. +  
   (5*NCi*x2)/8. +  
   (5*sqrtxz3*itani1*NCi*x2)/16. -  
   (5*sqrtxz3*itani2*NCi*x2)/16. +  
   (5*sqrtxz3*itani3*NCi*x2)/8. +  
   (5*sqrtxz3*atanspec1*lspec3*NCi*x2)/8. +  
   (5*lx*NCi*x2)/16. - (5*lz*NCi*x2)/16. -  
   (5*sqrtxz3*atanspec2*lspec4*NCi*x2)/8. +  
   (3*NC*zi2)/8. - (3*NC*x*zi2)/8. +  
   (3*NC*sqrtxz3*itani1*zi2)/16. -  
   (3*NC*sqrtxz3*itani2*zi2)/16. +  
   (3*NC*sqrtxz3*itani3*zi2)/8. +  
   (3*NC*sqrtxz3*atanspec1*lspec3*zi2)/8. +  
   (3*NC*lx*zi2)/16. + (3*NC*x*lx*zi2)/16. -  
   (3*NC*lz*zi2)/16. + (3*NC*x*lz*zi2)/16. -  
   (3*NC*sqrtxz3*atanspec2*lspec4*zi2)/8. -  
   (3*NCi*zi2)/8. + (3*x*NCi*zi2)/8. -  
   (3*sqrtxz3*itani1*NCi*zi2)/16. +  
   (3*sqrtxz3*itani2*NCi*zi2)/16. -  
   (3*sqrtxz3*itani3*NCi*zi2)/8. -  
   (3*sqrtxz3*atanspec1*lspec3*NCi*zi2)/8. -  
   (3*lx*NCi*zi2)/16. -  
   (3*x*lx*NCi*zi2)/16. +  
   (3*lz*NCi*zi2)/16. -  
   (3*x*lz*NCi*zi2)/16. +  
   (3*sqrtxz3*atanspec2*lspec4*NCi*zi2)/8. -  
   (7*NC*zi)/4. + (3*NC*x*zi)/4. +  
   (3*NC*sqrtxz3*x*itani1*zi)/8. -  
   (3*NC*sqrtxz3*x*itani2*zi)/8. +  
   (3*NC*sqrtxz3*x*itani3*zi)/4. +  
   (3*NC*sqrtxz3*x*atanspec1*lspec3*zi)/4. -  
   (NC*lx*zi)/8. - (13*NC*x*lx*zi)/8. -  
   (NC*lz*zi)/8. - (3*NC*x*lz*zi)/8. -  
   (3*NC*sqrtxz3*x*atanspec2*lspec4*zi)/4. +  
   (7*NCi*zi)/4. - (3*x*NCi*zi)/4. -  
   (3*sqrtxz3*x*itani1*NCi*zi)/8. +  
   (3*sqrtxz3*x*itani2*NCi*zi)/8. -  
   (3*sqrtxz3*x*itani3*NCi*zi)/4. -  
   (3*sqrtxz3*x*atanspec1*lspec3*NCi*zi)/4. +  
   (lx*NCi*zi)/8. + (13*x*lx*NCi*zi)/8. +  
   (lz*NCi*zi)/8. + (3*x*lz*NCi*zi)/8. +  
   (3*sqrtxz3*x*atanspec2*lspec4*NCi*zi)/4. +  
   (3*NC*xi*zi)/8. -  
   (NC*sqrtxz3*itani1*xi*zi)/8. +  
   (NC*sqrtxz3*itani2*xi*zi)/8. -  
   (NC*sqrtxz3*itani3*xi*zi)/4. -  
   (NC*sqrtxz3*atanspec1*lspec3*xi*zi)/4. -  
   (3*NC*lx*xi*zi)/16. +  
   (3*NC*lz*xi*zi)/16. +  
   (NC*sqrtxz3*atanspec2*lspec4*xi*zi)/4. -  
   (3*NCi*xi*zi)/8. +  
   (sqrtxz3*itani1*NCi*xi*zi)/8. -  
   (sqrtxz3*itani2*NCi*xi*zi)/8. +  
   (sqrtxz3*itani3*NCi*xi*zi)/4. +  
   (sqrtxz3*atanspec1*lspec3*NCi*xi*zi)/4. +  
   (3*lx*NCi*xi*zi)/16. -  
   (3*lz*NCi*xi*zi)/16. -  
   (sqrtxz3*atanspec2*lspec4*NCi*xi*zi)/ 
    4. + (5*NC*x2*zi)/8. +  
   (5*NC*lx*x2*zi)/16. +  
   (5*NC*lz*x2*zi)/16. -  
   (5*NCi*x2*zi)/8. -  
   (5*lx*NCi*x2*zi)/16. -  
   (5*lz*NCi*x2*zi)/16. -  
   (5*NC*sqrtxz3*itani1*z2)/16. +  
   (5*NC*sqrtxz3*itani2*z2)/16. -  
   (5*NC*sqrtxz3*itani3*z2)/8. -  
   (5*NC*sqrtxz3*atanspec1*lspec3*z2)/8. +  
   (5*NC*sqrtxz3*atanspec2*lspec4*z2)/8. +  
   (5*sqrtxz3*itani1*NCi*z2)/16. -  
   (5*sqrtxz3*itani2*NCi*z2)/16. +  
   (5*sqrtxz3*itani3*NCi*z2)/8. +  
   (5*sqrtxz3*atanspec1*lspec3*NCi*z2)/8. -  
   (5*sqrtxz3*atanspec2*lspec4*NCi*z2)/8.;
};
double C2LQ2QB_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double sqrtxz1  = sqrt(1 - 2*z + z*z + 4*x*z);
  const double sqrtxz1i = 1. / sqrtxz1;
  const double poly2    = 1 + 2*x + x*x - 4*x*z;
  const double poly2i   = 1. / poly2;
  const double sqrtxz2  = sqrt(poly2);
  const double sqrtxz2i = 1. / sqrtxz2;
  const double sqrtxz3  = sqrt(x/z);
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double NCi2 = NCi * NCi;
  const double l2  = log(2);
  const double l22 = l2 * l2;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double x2 = x * x;
  const double x3 = x * x2;
  const double x4 = x * x3;
  const double x5 = x * x4;
  const double xi  = 1. / x;
  const double z2 = z * z;
  const double zi  = 1. / z;
  const double opxi = 1. / ( 1 + x );
  const double omzi  = 1. / ( 1 - z );
  const double opzi  = 1. / ( 1 + z );
  const double mopzi = 1. / ( - 1 + z );
  const double lx  = log(x);
  const double lx2 = lx * lx;
  const double lz  = log(z);
  const double lz2 = lz * lz;
  const double lomx  = log(1 - x);
  const double lomz  = log(1 - z);
  const double lomz2 = lomz * lomz;
  const double lopx  = log(1 + x);
  const double Li2x  = apfel::dilog(x);
  const double Li2mx = apfel::dilog(-x);
  const double Li2z  = apfel::dilog(z);
  const double xmzi  = 1. / ( x - z );
  const double spec1i = 1. / ( 1 + sqrtxz1 - z );
  const double lxpz    = log(x + z);
  const double lopxz   = log(1 + x*z);
  const double lopxzi  = log(1 + x*zi);
  const double li2omxzi = apfel::dilog(1 - x*zi);
  const double lspec1   = log(1 + sqrtxz1 - z);
  const double lspec1_2 = lspec1 * lspec1;
  const double lspec2   = log(1 + sqrtxz1 + z);
  const double lspec3   = log(sqrtxz3);
  const double lspec4   = log(sqrtxz3*z);
  const double lspec5   = log(1 - sqrtxz2 + x);
  const double lspec6   = log(1 + sqrtxz2 + x);
  const double lspec7   = log(1 - 2*z + 4*x*z + z2);
  const double lspec7_2 = lspec7 * lspec7;
  const double li2spec1  = apfel::dilog(0.5 - sqrtxz1/2. - z/2.);
  const double li2spec2  = apfel::dilog(0.5 - sqrtxz1/2. + z/2.);
  const double li2spec3  = apfel::dilog(0.5 - zi/2. - (sqrtxz1*zi)/2.);
  const double li2spec4  = apfel::dilog(0.5 + zi/2. - (sqrtxz1*zi)/2.);
  const double li2spec5  = apfel::dilog(0.5 - sqrtxz2/2. - x/2.);
  const double li2spec6  = apfel::dilog(0.5 + sqrtxz2/2. - x/2.);
  const double li2spec7  = apfel::dilog(0.5 - xi/2. - (sqrtxz2*xi)/2.);
  const double li2spec8  = apfel::dilog(0.5 - xi/2. + (sqrtxz2*xi)/2.);
  const double li2spec19 = apfel::dilog(-(x*z));
  const double li2spec20 = apfel::dilog(-(x*zi));
  const double li2spec21 = apfel::dilog((1 - z)*sqrtxz1i);
  const double li2spec22 = apfel::dilog(-spec1i + sqrtxz1*spec1i + z*spec1i);
  const double li2spec23 = apfel::dilog(sqrtxz1*mopzi);
  const double atanspec1 = atan(sqrtxz3);
  const double atanspec2 = atan(sqrtxz3*z);
  const double itani1 = InvTanInt(-sqrtxz3);
  const double itani2 = InvTanInt(sqrtxz3);
  const double itani3 = InvTanInt(sqrtxz3*z);
  return -6 - (34*NC)/3. + 6*x + (41*NC*x)/6. + 10*z + 12*NC*z - 10*x*z -  
   11*NC*x*z + (NC*sqrtxz3*itani1)/2. +  
   (3*NC*sqrtxz3*x*z*itani1)/2. -  
   (NC*sqrtxz3*itani2)/2. -  
   (3*NC*sqrtxz3*x*z*itani2)/2. +  
   NC*sqrtxz3*itani3 + 3*NC*sqrtxz3*x*z*itani3 +  
   4*x*Li2mx + 4*NC*x*Li2mx - 8*NC*x*z*Li2mx - 2*x*Li2x +  
   2*NC*x*Li2x + 2*x*z*Li2x - 4*NC*x*z*Li2x +  
   4*x*li2spec1 - 4*x*z*li2spec1 -  
   4*x*li2spec2 + 4*x*z*li2spec2 +  
   3*NC*x*Li2z - 2*NC*x*z*Li2z - 4*x*li2spec19 -  
   4*x*z*li2spec20 - 4*NC*x*z*li2spec20 -  
   4*x*li2spec3 +  
   4*x*z*li2spec3 +  
   4*NC*x*z*li2spec3 +  
   4*x*li2spec4 -  
   4*x*z*li2spec4 -  
   4*NC*x*z*li2spec4 -  
   2*NC*x*li2omxzi + 4*NC*x*z*li2omxzi +  
   4*NC*sqrtxz1*x*l2 + NC*sqrtxz3*atanspec1*lspec3 +  
   3*NC*sqrtxz3*x*z*atanspec1*lspec3 - 2*NC*lomx -  
   2*NC*x*lomx + NC*x*z*lomx - 2*lx - NC*lx -  
   2*x*lx - 7*NC*x*lx + 2*NC*sqrtxz1*x*lx + 4*z*lx +  
   (11*NC*z*lx)/2. + 4*x*z*lx + (11*NC*x*z*lx)/2. -  
   4*x*l2*lx + 4*x*z*l2*lx - 4*NC*x*z*l2*lx -  
   2*x*lomx*lx + 2*x*z*lomx*lx -  
   4*NC*x*z*lomx*lx + 4*x*lx*lopx +  
   4*NC*x*lx*lopx - 8*NC*x*z*lx*lopx - 2*NC*lomz -  
   2*NC*x*lomz + NC*x*z*lomz - 2*NC*x*lx*lomz -  
   4*NC*sqrtxz1*x*lspec1 + 8*x*l2*lspec1 -  
   8*x*z*l2*lspec1 + 4*NC*x*z*l2*lspec1 +  
   4*NC*x*z*lx*lspec1 - 5*NC*lz +  
   2*NC*sqrtxz1*x*lz + 2*z*lz + (7*NC*z*lz)/2. - 2*x*z*lz -  
   (7*NC*x*z*lz)/2. - 12*x*l2*lz + 12*x*z*l2*lz +  
   4*NC*x*z*l2*lz - NC*x*lomx*lz -  
   2*NC*x*z*lomx*lz - 2*x*lx*lz + NC*x*lx*lz +  
   4*x*z*lx*lz + 2*NC*x*z*lx*lz + 2*NC*x*lomz*lz -  
   4*NC*x*z*lomz*lz + 4*x*lspec1*lz -  
   4*x*z*lspec1*lz -  
   NC*sqrtxz3*atanspec2*lspec4 -  
   3*NC*sqrtxz3*x*z*atanspec2*lspec4 +  
   8*x*l2*lspec2 - 8*x*z*l2*lspec2 -  
   4*NC*x*z*l2*lspec2 + 4*x*lx*lspec2 -  
   4*x*z*lx*lspec2 -  
   8*x*lspec1*lspec2 +  
   8*x*z*lspec1*lspec2 +  
   4*NC*x*z*lspec1*lspec2 +  
   8*x*lz*lspec2 - 8*x*z*lz*lspec2 -  
   4*NC*x*z*lz*lspec2 + 2*x*lx*lxpz -  
   2*x*z*lx*lxpz - 2*NC*x*z*lx*lxpz -  
   2*x*lz*lxpz + 2*x*z*lz*lxpz +  
   2*NC*x*z*lz*lxpz - 4*x*lx*lopxz -  
   4*x*lz*lopxz - 2*x*lx*lopxzi -  
   2*x*z*lx*lopxzi - 2*NC*x*z*lx*lopxzi +  
   2*x*lz*lopxzi + 2*x*z*lz*lopxzi +  
   2*NC*x*z*lz*lopxzi + 6*NCi2 - 6*x*NCi2 -  
   10*z*NCi2 + 10*x*z*NCi2 - 4*x*Li2mx*NCi2 +  
   2*x*Li2x*NCi2 - 2*x*z*Li2x*NCi2 -  
   4*x*li2spec1*NCi2 +  
   4*x*z*li2spec1*NCi2 +  
   4*x*li2spec2*NCi2 -  
   4*x*z*li2spec2*NCi2 +  
   4*x*li2spec19*NCi2 + 4*x*z*li2spec20*NCi2 +  
   4*x*li2spec3*NCi2 -  
   4*x*z*li2spec3*NCi2 -  
   4*x*li2spec4*NCi2 +  
   4*x*z*li2spec4*NCi2 +  
   2*lx*NCi2 + 2*x*lx*NCi2 - 4*z*lx*NCi2 -  
   4*x*z*lx*NCi2 + 4*x*l2*lx*NCi2 -  
   4*x*z*l2*lx*NCi2 + 2*x*lomx*lx*NCi2 -  
   2*x*z*lomx*lx*NCi2 - 4*x*lx*lopx*NCi2 -  
   8*x*l2*lspec1*NCi2 +  
   8*x*z*l2*lspec1*NCi2 - 2*z*lz*NCi2 +  
   2*x*z*lz*NCi2 + 12*x*l2*lz*NCi2 -  
   12*x*z*l2*lz*NCi2 + 2*x*lx*lz*NCi2 -  
   4*x*z*lx*lz*NCi2 -  
   4*x*lspec1*lz*NCi2 +  
   4*x*z*lspec1*lz*NCi2 -  
   8*x*l2*lspec2*NCi2 +  
   8*x*z*l2*lspec2*NCi2 -  
   4*x*lx*lspec2*NCi2 +  
   4*x*z*lx*lspec2*NCi2 +  
   8*x*lspec1*lspec2*NCi2 -  
   8*x*z*lspec1*lspec2*NCi2 -  
   8*x*lz*lspec2*NCi2 +  
   8*x*z*lz*lspec2*NCi2 -  
   2*x*lx*lxpz*NCi2 + 2*x*z*lx*lxpz*NCi2 +  
   2*x*lz*lxpz*NCi2 - 2*x*z*lz*lxpz*NCi2 +  
   4*x*lx*lopxz*NCi2 + 4*x*lz*lopxz*NCi2 +  
   2*x*lx*lopxzi*NCi2 +  
   2*x*z*lx*lopxzi*NCi2 -  
   2*x*lz*lopxzi*NCi2 -  
   2*x*z*lz*lopxzi*NCi2 + (34*NCi)/3. -  
   (41*x*NCi)/6. - 12*z*NCi + 11*x*z*NCi -  
   (sqrtxz3*itani1*NCi)/2. -  
   (3*sqrtxz3*x*z*itani1*NCi)/2. +  
   (sqrtxz3*itani2*NCi)/2. +  
   (3*sqrtxz3*x*z*itani2*NCi)/2. -  
   sqrtxz3*itani3*NCi -  
   3*sqrtxz3*x*z*itani3*NCi - 4*x*Li2mx*NCi +  
   8*x*z*Li2mx*NCi - 2*x*Li2x*NCi +  
   4*x*z*Li2x*NCi - 3*x*Li2z*NCi +  
   2*x*z*Li2z*NCi + 4*x*z*li2spec20*NCi -  
   4*x*z*li2spec3*NCi +  
   4*x*z*li2spec4*NCi +  
   2*x*li2omxzi*NCi -  
   4*x*z*li2omxzi*NCi - 4*sqrtxz1*x*l2*NCi -  
   sqrtxz3*atanspec1*lspec3*NCi -  
   3*sqrtxz3*x*z*atanspec1*lspec3*NCi +  
   2*lomx*NCi + 2*x*lomx*NCi -  
   x*z*lomx*NCi + lx*NCi + 7*x*lx*NCi -  
   2*sqrtxz1*x*lx*NCi - (11*z*lx*NCi)/2. -  
   (11*x*z*lx*NCi)/2. + 4*x*z*l2*lx*NCi +  
   4*x*z*lomx*lx*NCi - 4*x*lx*lopx*NCi +  
   8*x*z*lx*lopx*NCi + 2*lomz*NCi +  
   2*x*lomz*NCi - x*z*lomz*NCi +  
   2*x*lx*lomz*NCi +  
   4*sqrtxz1*x*lspec1*NCi -  
   4*x*z*l2*lspec1*NCi -  
   4*x*z*lx*lspec1*NCi + 5*lz*NCi -  
   2*sqrtxz1*x*lz*NCi - (7*z*lz*NCi)/2. +  
   (7*x*z*lz*NCi)/2. - 4*x*z*l2*lz*NCi +  
   x*lomx*lz*NCi + 2*x*z*lomx*lz*NCi -  
   x*lx*lz*NCi - 2*x*z*lx*lz*NCi -  
   2*x*lomz*lz*NCi + 4*x*z*lomz*lz*NCi +  
   sqrtxz3*atanspec2*lspec4*NCi +  
   3*sqrtxz3*x*z*atanspec2*lspec4*NCi +  
   4*x*z*l2*lspec2*NCi -  
   4*x*z*lspec1*lspec2*NCi +  
   4*x*z*lz*lspec2*NCi +  
   2*x*z*lx*lxpz*NCi - 2*x*z*lz*lxpz*NCi +  
   2*x*z*lx*lopxzi*NCi -  
   2*x*z*lz*lopxzi*NCi + (x*pi2)/3. -  
   (NC*x*pi2)/2. - (2*x*z*pi2)/3. + (NC*x*z*pi2)/3. -  
   (x*NCi2*pi2)/3. + (2*x*z*NCi2*pi2)/3. +  
   (x*NCi*pi2)/2. - (x*z*NCi*pi2)/3. +  
   (NC*lx*poly2i)/2. + NC*x*lx*poly2i -  
   (NC*lz*poly2i)/2. + NC*x*lz*poly2i -  
   (lx*NCi*poly2i)/2. - x*lx*NCi*poly2i +  
   (lz*NCi*poly2i)/2. - x*lz*NCi*poly2i -  
   8*x*li2spec1*sqrtxz1i +  
   4*x*z*li2spec1*sqrtxz1i -  
   8*x*li2spec2*sqrtxz1i +  
   4*x*z*li2spec2*sqrtxz1i +  
   16*x*li2spec21*sqrtxz1i -  
   8*x*z*li2spec21*sqrtxz1i +  
   8*x*li2spec22*sqrtxz1i -  
   4*x*z*li2spec22*sqrtxz1i +  
   16*x*li2spec23*sqrtxz1i -  
   8*x*z*li2spec23*sqrtxz1i +  
   8*x*li2spec3*sqrtxz1i -  
   4*x*z*li2spec3*sqrtxz1i +  
   8*x*li2spec4*sqrtxz1i -  
   4*x*z*li2spec4*sqrtxz1i +  
   32*x*l2*sqrtxz1i - 16*x*z*l2*sqrtxz1i +  
   16*x*lx*sqrtxz1i - 8*x*z*lx*sqrtxz1i +  
   16*x*l2*lomz*sqrtxz1i -  
   8*x*z*l2*lomz*sqrtxz1i +  
   8*x*lx*lomz*sqrtxz1i -  
   4*x*z*lx*lomz*sqrtxz1i -  
   32*x*lspec1*sqrtxz1i +  
   16*x*z*lspec1*sqrtxz1i -  
   32*x*l2*lspec1*sqrtxz1i +  
   16*x*z*l2*lspec1*sqrtxz1i +  
   8*x*lx*lspec1*sqrtxz1i -  
   4*x*z*lx*lspec1*sqrtxz1i -  
   16*x*lomz*lspec1*sqrtxz1i +  
   8*x*z*lomz*lspec1*sqrtxz1i +  
   16*x*lz*sqrtxz1i - 8*x*z*lz*sqrtxz1i +  
   32*x*l2*lz*sqrtxz1i -  
   16*x*z*l2*lz*sqrtxz1i - 4*x*lx*lz*sqrtxz1i +  
   2*x*z*lx*lz*sqrtxz1i +  
   8*x*lomz*lz*sqrtxz1i -  
   4*x*z*lomz*lz*sqrtxz1i -  
   16*x*lspec1*lz*sqrtxz1i +  
   8*x*z*lspec1*lz*sqrtxz1i -  
   16*x*l2*lspec2*sqrtxz1i +  
   8*x*z*l2*lspec2*sqrtxz1i -  
   8*x*lx*lspec2*sqrtxz1i +  
   4*x*z*lx*lspec2*sqrtxz1i +  
   16*x*lspec1*lspec2*sqrtxz1i -  
   8*x*z*lspec1*lspec2*sqrtxz1i -  
   16*x*lz*lspec2*sqrtxz1i +  
   8*x*z*lz*lspec2*sqrtxz1i -  
   8*x*lomz*lspec7*sqrtxz1i +  
   4*x*z*lomz*lspec7*sqrtxz1i +  
   8*x*li2spec1*NCi2*sqrtxz1i -  
   4*x*z*li2spec1*NCi2*sqrtxz1i +  
   8*x*li2spec2*NCi2*sqrtxz1i -  
   4*x*z*li2spec2*NCi2*sqrtxz1i -  
   16*x*li2spec21*NCi2*sqrtxz1i +  
   8*x*z*li2spec21*NCi2*sqrtxz1i -  
   8*x*li2spec22*NCi2*sqrtxz1i +  
   4*x*z*li2spec22*NCi2*sqrtxz1i -  
   16*x*li2spec23*NCi2*sqrtxz1i +  
   8*x*z*li2spec23*NCi2*sqrtxz1i -  
   8*x*li2spec3*NCi2* 
    sqrtxz1i + 4*x*z*li2spec3* 
    NCi2*sqrtxz1i -  
   8*x*li2spec4*NCi2* 
    sqrtxz1i + 4*x*z*li2spec4* 
    NCi2*sqrtxz1i - 32*x*l2*NCi2*sqrtxz1i +  
   16*x*z*l2*NCi2*sqrtxz1i -  
   16*x*lx*NCi2*sqrtxz1i +  
   8*x*z*lx*NCi2*sqrtxz1i -  
   16*x*l2*lomz*NCi2*sqrtxz1i +  
   8*x*z*l2*lomz*NCi2*sqrtxz1i -  
   8*x*lx*lomz*NCi2*sqrtxz1i +  
   4*x*z*lx*lomz*NCi2*sqrtxz1i +  
   32*x*lspec1*NCi2*sqrtxz1i -  
   16*x*z*lspec1*NCi2*sqrtxz1i +  
   32*x*l2*lspec1*NCi2*sqrtxz1i -  
   16*x*z*l2*lspec1*NCi2*sqrtxz1i -  
   8*x*lx*lspec1*NCi2*sqrtxz1i +  
   4*x*z*lx*lspec1*NCi2*sqrtxz1i +  
   16*x*lomz*lspec1*NCi2*sqrtxz1i -  
   8*x*z*lomz*lspec1*NCi2*sqrtxz1i -  
   16*x*lz*NCi2*sqrtxz1i +  
   8*x*z*lz*NCi2*sqrtxz1i -  
   32*x*l2*lz*NCi2*sqrtxz1i +  
   16*x*z*l2*lz*NCi2*sqrtxz1i +  
   4*x*lx*lz*NCi2*sqrtxz1i -  
   2*x*z*lx*lz*NCi2*sqrtxz1i -  
   8*x*lomz*lz*NCi2*sqrtxz1i +  
   4*x*z*lomz*lz*NCi2*sqrtxz1i +  
   16*x*lspec1*lz*NCi2*sqrtxz1i -  
   8*x*z*lspec1*lz*NCi2*sqrtxz1i +  
   16*x*l2*lspec2*NCi2*sqrtxz1i -  
   8*x*z*l2*lspec2*NCi2*sqrtxz1i +  
   8*x*lx*lspec2*NCi2*sqrtxz1i -  
   4*x*z*lx*lspec2*NCi2*sqrtxz1i -  
   16*x*lspec1*lspec2*NCi2* 
    sqrtxz1i + 8*x*z*lspec1*lspec2* 
    NCi2*sqrtxz1i +  
   16*x*lz*lspec2*NCi2*sqrtxz1i -  
   8*x*z*lz*lspec2*NCi2*sqrtxz1i +  
   8*x*lomz*lspec7*NCi2* 
    sqrtxz1i - 4*x*z*lomz*lspec7* 
    NCi2*sqrtxz1i + (4*x*pi2*sqrtxz1i)/3. -  
   (2*x*z*pi2*sqrtxz1i)/3. -  
   (4*x*NCi2*pi2*sqrtxz1i)/3. +  
   (2*x*z*NCi2*pi2*sqrtxz1i)/3. +  
   (NC*li2spec5*sqrtxz2i)/2. -  
   NC*z*li2spec5*sqrtxz2i -  
   4*NC*x*z*li2spec5*sqrtxz2i -  
   (NC*li2spec6*sqrtxz2i)/2. +  
   NC*z*li2spec6*sqrtxz2i +  
   4*NC*x*z*li2spec6*sqrtxz2i -  
   (NC*li2spec7*sqrtxz2i)/2. +  
   NC*z*li2spec7*sqrtxz2i +  
   4*NC*x*z*li2spec7* 
    sqrtxz2i + (NC*li2spec8* 
      sqrtxz2i)/2. - NC*z* 
    li2spec8*sqrtxz2i -  
   4*NC*x*z*li2spec8* 
    sqrtxz2i - (NC*lx*lspec5*sqrtxz2i)/2. +  
   NC*z*lx*lspec5*sqrtxz2i +  
   4*NC*x*z*lx*lspec5*sqrtxz2i +  
   (NC*lx*lspec6*sqrtxz2i)/2. -  
   NC*z*lx*lspec6*sqrtxz2i -  
   4*NC*x*z*lx*lspec6*sqrtxz2i -  
   (li2spec5*NCi*sqrtxz2i)/2. +  
   z*li2spec5*NCi*sqrtxz2i +  
   4*x*z*li2spec5*NCi*sqrtxz2i +  
   (li2spec6*NCi*sqrtxz2i)/2. -  
   z*li2spec6*NCi*sqrtxz2i -  
   4*x*z*li2spec6*NCi*sqrtxz2i +  
   (li2spec7*NCi* 
      sqrtxz2i)/2. - z*li2spec7*NCi*sqrtxz2i -  
   4*x*z*li2spec7*NCi* 
    sqrtxz2i - (li2spec8* 
      NCi*sqrtxz2i)/2. +  
   z*li2spec8*NCi* 
    sqrtxz2i + 4*x*z*li2spec8* 
    NCi*sqrtxz2i +  
   (lx*lspec5*NCi*sqrtxz2i)/2. -  
   z*lx*lspec5*NCi*sqrtxz2i -  
   4*x*z*lx*lspec5*NCi*sqrtxz2i -  
   (lx*lspec6*NCi*sqrtxz2i)/2. +  
   z*lx*lspec6*NCi*sqrtxz2i +  
   4*x*z*lx*lspec6*NCi*sqrtxz2i -  
   (3*NC*x*li2spec5*poly2i*sqrtxz2i)/4. +  
   (3*NC*x*li2spec6*poly2i*sqrtxz2i)/4. +  
   (3*NC*x*li2spec7*poly2i* 
      sqrtxz2i)/4. - (3*NC*x* 
      li2spec8*poly2i* 
      sqrtxz2i)/4. + (3*NC*x*lx*lspec5*poly2i* 
      sqrtxz2i)/4. - (3*NC*x*lx*lspec6*poly2i* 
      sqrtxz2i)/4. + (3*x*li2spec5*NCi* 
      poly2i*sqrtxz2i)/4. -  
   (3*x*li2spec6*NCi*poly2i* 
      sqrtxz2i)/4. - (3*x* 
      li2spec7*NCi* 
      poly2i*sqrtxz2i)/4. +  
   (3*x*li2spec8*NCi* 
      poly2i*sqrtxz2i)/4. -  
   (3*x*lx*lspec5*NCi*poly2i*sqrtxz2i)/ 
    4. + (3*x*lx*lspec6*NCi*poly2i* 
      sqrtxz2i)/4. + (10*NC*xi)/9. -  
   (NC*sqrtxz3*z*itani1*xi)/2. +  
   (NC*sqrtxz3*z*itani2*xi)/2. -  
   NC*sqrtxz3*z*itani3*xi - 4*Li2mx*xi -  
   4*z*Li2mx*xi - 4*NC*z*Li2mx*xi + 4*Li2x*xi +  
   4*z*Li2x*xi + 4*NC*z*Li2x*xi + 2*li2spec19*xi +  
   2*z*li2spec19*xi + 2*NC*z*li2spec19*xi +  
   2*li2spec20*xi + 2*z*li2spec20*xi +  
   2*NC*z*li2spec20*xi -  
   NC*sqrtxz3*z*atanspec1*lspec3*xi +  
   (2*NC*lomx*xi)/3. - 8*lx*xi +  
   (NC*lx*xi)/2. - 8*z*lx*xi -  
   8*NC*z*lx*xi + 4*lomx*lx*xi +  
   4*z*lomx*lx*xi + 4*NC*z*lomx*lx*xi -  
   4*lx*lopx*xi - 4*z*lx*lopx*xi -  
   4*NC*z*lx*lopx*xi + (2*NC*lomz*xi)/3. +  
   (7*NC*lz*xi)/6. - 2*lx*lz*xi -  
   2*z*lx*lz*xi - 2*NC*z*lx*lz*xi +  
   NC*sqrtxz3*z*atanspec2*lspec4*xi +  
   2*lx*lopxz*xi + 2*z*lx*lopxz*xi +  
   2*NC*z*lx*lopxz*xi + 2*lz*lopxz*xi +  
   2*z*lz*lopxz*xi + 2*NC*z*lz*lopxz*xi +  
   2*lx*lopxzi*xi +  
   2*z*lx*lopxzi*xi +  
   2*NC*z*lx*lopxzi*xi -  
   2*lz*lopxzi*xi -  
   2*z*lz*lopxzi*xi -  
   2*NC*z*lz*lopxzi*xi +  
   4*Li2mx*NCi2*xi + 4*z*Li2mx*NCi2*xi -  
   4*Li2x*NCi2*xi - 4*z*Li2x*NCi2*xi -  
   2*li2spec19*NCi2*xi -  
   2*z*li2spec19*NCi2*xi -  
   2*li2spec20*NCi2*xi -  
   2*z*li2spec20*NCi2*xi +  
   8*lx*NCi2*xi + 8*z*lx*NCi2*xi -  
   4*lomx*lx*NCi2*xi -  
   4*z*lomx*lx*NCi2*xi +  
   4*lx*lopx*NCi2*xi +  
   4*z*lx*lopx*NCi2*xi +  
   2*lx*lz*NCi2*xi +  
   2*z*lx*lz*NCi2*xi -  
   2*lx*lopxz*NCi2*xi -  
   2*z*lx*lopxz*NCi2*xi -  
   2*lz*lopxz*NCi2*xi -  
   2*z*lz*lopxz*NCi2*xi -  
   2*lx*lopxzi*NCi2*xi -  
   2*z*lx*lopxzi*NCi2*xi +  
   2*lz*lopxzi*NCi2*xi +  
   2*z*lz*lopxzi*NCi2*xi -  
   (10*NCi*xi)/9. +  
   (sqrtxz3*z*itani1*NCi*xi)/2. -  
   (sqrtxz3*z*itani2*NCi*xi)/2. +  
   sqrtxz3*z*itani3*NCi*xi +  
   4*z*Li2mx*NCi*xi - 4*z*Li2x*NCi*xi -  
   2*z*li2spec19*NCi*xi -  
   2*z*li2spec20*NCi*xi +  
   sqrtxz3*z*atanspec1*lspec3*NCi*xi -  
   (2*lomx*NCi*xi)/3. -  
   (lx*NCi*xi)/2. + 8*z*lx*NCi*xi -  
   4*z*lomx*lx*NCi*xi +  
   4*z*lx*lopx*NCi*xi -  
   (2*lomz*NCi*xi)/3. -  
   (7*lz*NCi*xi)/6. +  
   2*z*lx*lz*NCi*xi -  
   sqrtxz3*z*atanspec2*lspec4*NCi*xi -  
   2*z*lx*lopxz*NCi*xi -  
   2*z*lz*lopxz*NCi*xi -  
   2*z*lx*lopxzi*NCi*xi +  
   2*z*lz*lopxzi*NCi*xi -  
   (2*pi2*xi)/3. - (2*z*pi2*xi)/3. -  
   (2*NC*z*pi2*xi)/3. + (2*NCi2*pi2*xi)/3. +  
   (2*z*NCi2*pi2*xi)/3. +  
   (2*z*NCi*pi2*xi)/3. -  
   (NC*lx*poly2i*xi)/2. -  
   (NC*lz*poly2i*xi)/2. +  
   (lx*NCi*poly2i*xi)/2. +  
   (lz*NCi*poly2i*xi)/2. -  
   (NC*li2spec5*sqrtxz2i*xi)/4. +  
   (NC*li2spec6*sqrtxz2i*xi)/4. +  
   (NC*li2spec7*sqrtxz2i* 
      xi)/4. - (NC*li2spec8* 
      sqrtxz2i*xi)/4. +  
   (NC*lx*lspec5*sqrtxz2i*xi)/4. -  
   (NC*lx*lspec6*sqrtxz2i*xi)/4. +  
   (li2spec5*NCi*sqrtxz2i*xi)/4. -  
   (li2spec6*NCi*sqrtxz2i*xi)/4. -  
   (li2spec7*NCi* 
      sqrtxz2i*xi)/4. +  
   (li2spec8*NCi* 
      sqrtxz2i*xi)/4. -  
   (lx*lspec5*NCi*sqrtxz2i*xi)/4. +  
   (lx*lspec6*NCi*sqrtxz2i*xi)/4. +  
   (NC*li2spec5*poly2i*sqrtxz2i*xi)/ 
    4. - (NC*li2spec6*poly2i*sqrtxz2i* 
      xi)/4. - (NC*li2spec7* 
      poly2i*sqrtxz2i*xi)/4. +  
   (NC*li2spec8*poly2i* 
      sqrtxz2i*xi)/4. -  
   (NC*lx*lspec5*poly2i*sqrtxz2i*xi)/ 
    4. + (NC*lx*lspec6*poly2i*sqrtxz2i* 
      xi)/4. - (li2spec5*NCi*poly2i* 
      sqrtxz2i*xi)/4. +  
   (li2spec6*NCi*poly2i*sqrtxz2i* 
      xi)/4. + (li2spec7* 
      NCi*poly2i*sqrtxz2i*xi)/4. -  
   (li2spec8*NCi*poly2i* 
      sqrtxz2i*xi)/4. +  
   (lx*lspec5*NCi*poly2i*sqrtxz2i* 
      xi)/4. - (lx*lspec6*NCi*poly2i* 
      sqrtxz2i*xi)/4. + (44*NC*x2)/9. +  
   (4*NC*lomx*x2)/3. - (5*NC*lx*x2)/2. +  
   (4*NC*lomz*x2)/3. - (NC*lz*x2)/6. -  
   (44*NCi*x2)/9. - (4*lomx*NCi*x2)/3. +  
   (5*lx*NCi*x2)/2. -  
   (4*lomz*NCi*x2)/3. + (lz*NCi*x2)/6. -  
   NC*lx*poly2i*x2 + NC*lz*poly2i*x2 +  
   lx*NCi*poly2i*x2 -  
   lz*NCi*poly2i*x2 +  
   16*li2spec1*sqrtxz1i*x2 +  
   16*li2spec2*sqrtxz1i*x2 -  
   32*li2spec21*sqrtxz1i*x2 -  
   16*li2spec22*sqrtxz1i*x2 -  
   32*li2spec23*sqrtxz1i*x2 -  
   16*li2spec3*sqrtxz1i* 
    x2 - 16*li2spec4* 
    sqrtxz1i*x2 - 64*l2*sqrtxz1i*x2 -  
   32*lx*sqrtxz1i*x2 -  
   32*l2*lomz*sqrtxz1i*x2 -  
   16*lx*lomz*sqrtxz1i*x2 +  
   64*lspec1*sqrtxz1i*x2 +  
   64*l2*lspec1*sqrtxz1i*x2 -  
   16*lx*lspec1*sqrtxz1i*x2 +  
   32*lomz*lspec1*sqrtxz1i*x2 -  
   32*lz*sqrtxz1i*x2 -  
   64*l2*lz*sqrtxz1i*x2 +  
   8*lx*lz*sqrtxz1i*x2 -  
   16*lomz*lz*sqrtxz1i*x2 +  
   32*lspec1*lz*sqrtxz1i*x2 +  
   32*l2*lspec2*sqrtxz1i*x2 +  
   16*lx*lspec2*sqrtxz1i*x2 -  
   32*lspec1*lspec2*sqrtxz1i*x2 +  
   32*lz*lspec2*sqrtxz1i*x2 +  
   16*lomz*lspec7*sqrtxz1i*x2 -  
   16*li2spec1*NCi2*sqrtxz1i*x2 -  
   16*li2spec2*NCi2*sqrtxz1i*x2 +  
   32*li2spec21*NCi2*sqrtxz1i*x2 +  
   16*li2spec22*NCi2*sqrtxz1i*x2 +  
   32*li2spec23*NCi2*sqrtxz1i*x2 +  
   16*li2spec3*NCi2* 
    sqrtxz1i*x2 + 16* 
    li2spec4*NCi2* 
    sqrtxz1i*x2 + 64*l2*NCi2*sqrtxz1i* 
    x2 + 32*lx*NCi2*sqrtxz1i*x2 +  
   32*l2*lomz*NCi2*sqrtxz1i*x2 +  
   16*lx*lomz*NCi2*sqrtxz1i*x2 -  
   64*lspec1*NCi2*sqrtxz1i*x2 -  
   64*l2*lspec1*NCi2*sqrtxz1i*x2 +  
   16*lx*lspec1*NCi2*sqrtxz1i*x2 -  
   32*lomz*lspec1*NCi2*sqrtxz1i*x2 +  
   32*lz*NCi2*sqrtxz1i*x2 +  
   64*l2*lz*NCi2*sqrtxz1i*x2 -  
   8*lx*lz*NCi2*sqrtxz1i*x2 +  
   16*lomz*lz*NCi2*sqrtxz1i*x2 -  
   32*lspec1*lz*NCi2*sqrtxz1i*x2 -  
   32*l2*lspec2*NCi2*sqrtxz1i*x2 -  
   16*lx*lspec2*NCi2*sqrtxz1i*x2 +  
   32*lspec1*lspec2*NCi2*sqrtxz1i* 
    x2 - 32*lz*lspec2*NCi2*sqrtxz1i* 
    x2 - 16*lomz*lspec7*NCi2* 
    sqrtxz1i*x2 - (8*pi2*sqrtxz1i*x2)/3. +  
   (8*NCi2*pi2*sqrtxz1i*x2)/3. -  
   (NC*li2spec5*sqrtxz2i*x2)/2. +  
   NC*z*li2spec5*sqrtxz2i*x2 +  
   (NC*li2spec6*sqrtxz2i*x2)/2. -  
   NC*z*li2spec6*sqrtxz2i*x2 +  
   (NC*li2spec7*sqrtxz2i* 
      x2)/2. - NC*z*li2spec7* 
    sqrtxz2i*x2 - (NC* 
      li2spec8*sqrtxz2i* 
      x2)/2. + NC*z*li2spec8* 
    sqrtxz2i*x2 + (NC*lx*lspec5* 
      sqrtxz2i*x2)/2. -  
   NC*z*lx*lspec5*sqrtxz2i*x2 -  
   (NC*lx*lspec6*sqrtxz2i*x2)/2. +  
   NC*z*lx*lspec6*sqrtxz2i*x2 +  
   (li2spec5*NCi*sqrtxz2i*x2)/2. -  
   z*li2spec5*NCi*sqrtxz2i*x2 -  
   (li2spec6*NCi*sqrtxz2i*x2)/2. +  
   z*li2spec6*NCi*sqrtxz2i*x2 -  
   (li2spec7*NCi* 
      sqrtxz2i*x2)/2. +  
   z*li2spec7*NCi* 
    sqrtxz2i*x2 + (li2spec8*NCi*sqrtxz2i*x2)/2. -  
   z*li2spec8*NCi* 
    sqrtxz2i*x2 - (lx*lspec5*NCi* 
      sqrtxz2i*x2)/2. +  
   z*lx*lspec5*NCi*sqrtxz2i*x2 +  
   (lx*lspec6*NCi*sqrtxz2i*x2)/2. -  
   z*lx*lspec6*NCi*sqrtxz2i*x2 -  
   (NC*lx*poly2i*x3)/2. -  
   (NC*lz*poly2i*x3)/2. +  
   (lx*NCi*poly2i*x3)/2. +  
   (lz*NCi*poly2i*x3)/2. +  
   (NC*li2spec5*sqrtxz2i*x3)/4. -  
   (NC*li2spec6*sqrtxz2i*x3)/4. -  
   (NC*li2spec7*sqrtxz2i* 
      x3)/4. + (NC*li2spec8* 
      sqrtxz2i*x3)/4. -  
   (NC*lx*lspec5*sqrtxz2i*x3)/4. +  
   (NC*lx*lspec6*sqrtxz2i*x3)/4. -  
   (li2spec5*NCi*sqrtxz2i*x3)/4. +  
   (li2spec6*NCi*sqrtxz2i*x3)/4. +  
   (li2spec7*NCi* 
      sqrtxz2i*x3)/4. -  
   (li2spec8*NCi* 
      sqrtxz2i*x3)/4. +  
   (lx*lspec5*NCi*sqrtxz2i*x3)/4. -  
   (lx*lspec6*NCi*sqrtxz2i*x3)/4. +  
   (3*NC*li2spec5*poly2i*sqrtxz2i*x3)/ 
    4. - (3*NC*li2spec6*poly2i*sqrtxz2i* 
      x3)/4. - (3*NC*li2spec7* 
      poly2i*sqrtxz2i*x3)/4. +  
   (3*NC*li2spec8*poly2i* 
      sqrtxz2i*x3)/4. -  
   (3*NC*lx*lspec5*poly2i*sqrtxz2i*x3)/ 
    4. + (3*NC*lx*lspec6*poly2i*sqrtxz2i* 
      x3)/4. - (3*li2spec5*NCi*poly2i* 
      sqrtxz2i*x3)/4. +  
   (3*li2spec6*NCi*poly2i*sqrtxz2i* 
      x3)/4. + (3*li2spec7* 
      NCi*poly2i*sqrtxz2i*x3)/4. -  
   (3*li2spec8*NCi* 
      poly2i*sqrtxz2i*x3)/4. +  
   (3*lx*lspec5*NCi*poly2i*sqrtxz2i* 
      x3)/4. - (3*lx*lspec6*NCi*poly2i* 
      sqrtxz2i*x3)/4. + (NC*lx*poly2i*x4)/2. -  
   (NC*lz*poly2i*x4)/2. -  
   (lx*NCi*poly2i*x4)/2. +  
   (lz*NCi*poly2i*x4)/2. -  
   (NC*li2spec5*poly2i*sqrtxz2i*x5)/ 
    4. + (NC*li2spec6*poly2i*sqrtxz2i* 
      x5)/4. + (NC*li2spec7* 
      poly2i*sqrtxz2i*x5)/4. -  
   (NC*li2spec8*poly2i* 
      sqrtxz2i*x5)/4. +  
   (NC*lx*lspec5*poly2i*sqrtxz2i*x5)/ 
    4. - (NC*lx*lspec6*poly2i*sqrtxz2i* 
      x5)/4. + (li2spec5*NCi*poly2i* 
      sqrtxz2i*x5)/4. -  
   (li2spec6*NCi*poly2i*sqrtxz2i* 
      x5)/4. - (li2spec7* 
      NCi*poly2i*sqrtxz2i*x5)/4. +  
   (li2spec8*NCi*poly2i* 
      sqrtxz2i*x5)/4. -  
   (lx*lspec5*NCi*poly2i*sqrtxz2i* 
      x5)/4. + (lx*lspec6*NCi*poly2i* 
      sqrtxz2i*x5)/4. + 4*Li2mx*opxi +  
   4*z*Li2mx*opxi + 4*NC*z*Li2mx*opxi -  
   4*Li2x*opxi - 4*z*Li2x*opxi -  
   4*NC*z*Li2x*opxi - 2*li2spec19*opxi -  
   2*z*li2spec19*opxi - 2*NC*z*li2spec19*opxi -  
   2*li2spec20*opxi -  
   2*z*li2spec20*opxi -  
   2*NC*z*li2spec20*opxi + 8*lx*opxi +  
   8*z*lx*opxi + 8*NC*z*lx*opxi -  
   4*lomx*lx*opxi - 4*z*lomx*lx*opxi -  
   4*NC*z*lomx*lx*opxi +  
   4*lx*lopx*opxi + 4*z*lx*lopx*opxi +  
   4*NC*z*lx*lopx*opxi + 2*lx*lz*opxi +  
   2*z*lx*lz*opxi + 2*NC*z*lx*lz*opxi -  
   2*lx*lopxz*opxi -  
   2*z*lx*lopxz*opxi -  
   2*NC*z*lx*lopxz*opxi -  
   2*lz*lopxz*opxi -  
   2*z*lz*lopxz*opxi -  
   2*NC*z*lz*lopxz*opxi -  
   2*lx*lopxzi*opxi -  
   2*z*lx*lopxzi*opxi -  
   2*NC*z*lx*lopxzi*opxi +  
   2*lz*lopxzi*opxi +  
   2*z*lz*lopxzi*opxi +  
   2*NC*z*lz*lopxzi*opxi -  
   4*Li2mx*NCi2*opxi -  
   4*z*Li2mx*NCi2*opxi + 4*Li2x*NCi2*opxi +  
   4*z*Li2x*NCi2*opxi +  
   2*li2spec19*NCi2*opxi +  
   2*z*li2spec19*NCi2*opxi +  
   2*li2spec20*NCi2*opxi +  
   2*z*li2spec20*NCi2*opxi -  
   8*lx*NCi2*opxi - 8*z*lx*NCi2*opxi +  
   4*lomx*lx*NCi2*opxi +  
   4*z*lomx*lx*NCi2*opxi -  
   4*lx*lopx*NCi2*opxi -  
   4*z*lx*lopx*NCi2*opxi -  
   2*lx*lz*NCi2*opxi -  
   2*z*lx*lz*NCi2*opxi +  
   2*lx*lopxz*NCi2*opxi +  
   2*z*lx*lopxz*NCi2*opxi +  
   2*lz*lopxz*NCi2*opxi +  
   2*z*lz*lopxz*NCi2*opxi +  
   2*lx*lopxzi*NCi2*opxi +  
   2*z*lx*lopxzi*NCi2*opxi -  
   2*lz*lopxzi*NCi2*opxi -  
   2*z*lz*lopxzi*NCi2*opxi -  
   4*z*Li2mx*NCi*opxi +  
   4*z*Li2x*NCi*opxi +  
   2*z*li2spec19*NCi*opxi +  
   2*z*li2spec20*NCi*opxi -  
   8*z*lx*NCi*opxi +  
   4*z*lomx*lx*NCi*opxi -  
   4*z*lx*lopx*NCi*opxi -  
   2*z*lx*lz*NCi*opxi +  
   2*z*lx*lopxz*NCi*opxi +  
   2*z*lz*lopxz*NCi*opxi +  
   2*z*lx*lopxzi*NCi*opxi -  
   2*z*lz*lopxzi*NCi*opxi +  
   (2*pi2*opxi)/3. + (2*z*pi2*opxi)/3. +  
   (2*NC*z*pi2*opxi)/3. -  
   (2*NCi2*pi2*opxi)/3. -  
   (2*z*NCi2*pi2*opxi)/3. -  
   (2*z*NCi*pi2*opxi)/3. +  
   4*Li2mx*xi*opxi + 4*z*Li2mx*xi*opxi +  
   4*NC*z*Li2mx*xi*opxi -  
   4*Li2x*xi*opxi - 4*z*Li2x*xi*opxi -  
   4*NC*z*Li2x*xi*opxi -  
   2*li2spec19*xi*opxi -  
   2*z*li2spec19*xi*opxi -  
   2*NC*z*li2spec19*xi*opxi -  
   2*li2spec20*xi*opxi -  
   2*z*li2spec20*xi*opxi -  
   2*NC*z*li2spec20*xi*opxi +  
   8*lx*xi*opxi + 8*z*lx*xi*opxi +  
   8*NC*z*lx*xi*opxi -  
   4*lomx*lx*xi*opxi -  
   4*z*lomx*lx*xi*opxi -  
   4*NC*z*lomx*lx*xi*opxi +  
   4*lx*lopx*xi*opxi +  
   4*z*lx*lopx*xi*opxi +  
   4*NC*z*lx*lopx*xi*opxi +  
   2*lx*lz*xi*opxi +  
   2*z*lx*lz*xi*opxi +  
   2*NC*z*lx*lz*xi*opxi -  
   2*lx*lopxz*xi*opxi -  
   2*z*lx*lopxz*xi*opxi -  
   2*NC*z*lx*lopxz*xi*opxi -  
   2*lz*lopxz*xi*opxi -  
   2*z*lz*lopxz*xi*opxi -  
   2*NC*z*lz*lopxz*xi*opxi -  
   2*lx*lopxzi*xi*opxi -  
   2*z*lx*lopxzi*xi*opxi -  
   2*NC*z*lx*lopxzi*xi*opxi +  
   2*lz*lopxzi*xi*opxi +  
   2*z*lz*lopxzi*xi*opxi +  
   2*NC*z*lz*lopxzi*xi*opxi -  
   4*Li2mx*NCi2*xi*opxi -  
   4*z*Li2mx*NCi2*xi*opxi +  
   4*Li2x*NCi2*xi*opxi +  
   4*z*Li2x*NCi2*xi*opxi +  
   2*li2spec19*NCi2*xi*opxi +  
   2*z*li2spec19*NCi2*xi*opxi +  
   2*li2spec20*NCi2*xi*opxi +  
   2*z*li2spec20*NCi2*xi*opxi -  
   8*lx*NCi2*xi*opxi -  
   8*z*lx*NCi2*xi*opxi +  
   4*lomx*lx*NCi2*xi*opxi +  
   4*z*lomx*lx*NCi2*xi*opxi -  
   4*lx*lopx*NCi2*xi*opxi -  
   4*z*lx*lopx*NCi2*xi*opxi -  
   2*lx*lz*NCi2*xi*opxi -  
   2*z*lx*lz*NCi2*xi*opxi +  
   2*lx*lopxz*NCi2*xi*opxi +  
   2*z*lx*lopxz*NCi2*xi*opxi +  
   2*lz*lopxz*NCi2*xi*opxi +  
   2*z*lz*lopxz*NCi2*xi*opxi +  
   2*lx*lopxzi*NCi2*xi*opxi +  
   2*z*lx*lopxzi*NCi2*xi*opxi -  
   2*lz*lopxzi*NCi2*xi*opxi -  
   2*z*lz*lopxzi*NCi2*xi*opxi -  
   4*z*Li2mx*NCi*xi*opxi +  
   4*z*Li2x*NCi*xi*opxi +  
   2*z*li2spec19*NCi*xi*opxi +  
   2*z*li2spec20*NCi*xi*opxi -  
   8*z*lx*NCi*xi*opxi +  
   4*z*lomx*lx*NCi*xi*opxi -  
   4*z*lx*lopx*NCi*xi*opxi -  
   2*z*lx*lz*NCi*xi*opxi +  
   2*z*lx*lopxz*NCi*xi*opxi +  
   2*z*lz*lopxz*NCi*xi*opxi +  
   2*z*lx*lopxzi*NCi*xi*opxi -  
   2*z*lz*lopxzi*NCi*xi*opxi +  
   (2*pi2*xi*opxi)/3. +  
   (2*z*pi2*xi*opxi)/3. +  
   (2*NC*z*pi2*xi*opxi)/3. -  
   (2*NCi2*pi2*xi*opxi)/3. -  
   (2*z*NCi2*pi2*xi*opxi)/3. -  
   (2*z*NCi*pi2*xi*opxi)/3. -  
   4*x*Li2x*omzi + 8*x*lx*omzi -  
   4*x*lomx*lx*omzi + 2*x*lx*lz*omzi +  
   4*x*Li2x*NCi2*omzi -  
   8*x*lx*NCi2*omzi +  
   4*x*lomx*lx*NCi2*omzi -  
   2*x*lx*lz*NCi2*omzi +  
   (2*x*pi2*omzi)/3. -  
   (2*x*NCi2*pi2*omzi)/3. +  
   2*NC*lx*x2*xmzi - 2*NC*lz*x2*xmzi -  
   2*lx*NCi*x2*xmzi +  
   2*lz*NCi*x2*xmzi -  
   2*NC*lx*x3*xmzi + 2*NC*lz*x3*xmzi +  
   2*lx*NCi*x3*xmzi -  
   2*lz*NCi*x3*xmzi + 4*x*Li2mx*zi -  
   4*x*Li2x*zi - 4*x*li2spec1*zi -  
   4*sqrtxz1*x*li2spec1*zi +  
   4*x*li2spec2*zi -  
   4*sqrtxz1*x*li2spec2*zi +  
   8*sqrtxz1*x*li2spec21*zi +  
   4*sqrtxz1*x*li2spec22*zi 
    + 8*sqrtxz1*x*li2spec23*zi -  
   4*x*li2spec20*zi +  
   4*x*li2spec3*zi +  
   4*sqrtxz1*x*li2spec3*zi -  
   4*x*li2spec4*zi +  
   4*sqrtxz1*x*li2spec4*zi +  
   16*sqrtxz1*x*l2*zi + 8*x*lx*zi +  
   8*sqrtxz1*x*lx*zi + 4*x*l2*lx*zi -  
   4*x*lomx*lx*zi + 4*x*lx*lopx*zi +  
   8*sqrtxz1*x*l2*lomz*zi +  
   4*sqrtxz1*x*lx*lomz*zi -  
   16*sqrtxz1*x*lspec1*zi -  
   8*x*l2*lspec1*zi -  
   16*sqrtxz1*x*l2*lspec1*zi +  
   4*sqrtxz1*x*lx*lspec1*zi -  
   8*sqrtxz1*x*lomz*lspec1*zi +  
   8*sqrtxz1*x*lz*zi + 12*x*l2*lz*zi +  
   16*sqrtxz1*x*l2*lz*zi + 4*x*lx*lz*zi -  
   2*sqrtxz1*x*lx*lz*zi +  
   4*sqrtxz1*x*lomz*lz*zi -  
   4*x*lspec1*lz*zi -  
   8*sqrtxz1*x*lspec1*lz*zi -  
   8*x*l2*lspec2*zi -  
   8*sqrtxz1*x*l2*lspec2*zi -  
   4*x*lx*lspec2*zi -  
   4*sqrtxz1*x*lx*lspec2*zi +  
   8*x*lspec1*lspec2*zi +  
   8*sqrtxz1*x*lspec1*lspec2*zi -  
   8*x*lz*lspec2*zi -  
   8*sqrtxz1*x*lz*lspec2*zi -  
   2*x*lx*lxpz*zi + 2*x*lz*lxpz*zi -  
   2*x*lx*lopxzi*zi +  
   2*x*lz*lopxzi*zi -  
   4*sqrtxz1*x*lomz*lspec7*zi -  
   4*x*Li2mx*NCi2*zi + 4*x*Li2x*NCi2*zi +  
   4*x*li2spec1*NCi2*zi +  
   4*sqrtxz1*x*li2spec1*NCi2*zi -  
   4*x*li2spec2*NCi2*zi +  
   4*sqrtxz1*x*li2spec2*NCi2*zi -  
   8*sqrtxz1*x*li2spec21*NCi2*zi -  
   4*sqrtxz1*x*li2spec22*NCi2* 
    zi - 8*sqrtxz1*x*li2spec23*NCi2* 
    zi + 4*x*li2spec20*NCi2*zi -  
   4*x*li2spec3*NCi2* 
    zi - 4*sqrtxz1*x*li2spec3* 
    NCi2*zi + 4*x*li2spec4*NCi2*zi -  
   4*sqrtxz1*x*li2spec4*NCi2* 
    zi - 16*sqrtxz1*x*l2*NCi2*zi -  
   8*x*lx*NCi2*zi -  
   8*sqrtxz1*x*lx*NCi2*zi -  
   4*x*l2*lx*NCi2*zi +  
   4*x*lomx*lx*NCi2*zi -  
   4*x*lx*lopx*NCi2*zi -  
   8*sqrtxz1*x*l2*lomz*NCi2*zi -  
   4*sqrtxz1*x*lx*lomz*NCi2*zi +  
   16*sqrtxz1*x*lspec1*NCi2*zi +  
   8*x*l2*lspec1*NCi2*zi +  
   16*sqrtxz1*x*l2*lspec1*NCi2*zi -  
   4*sqrtxz1*x*lx*lspec1*NCi2*zi +  
   8*sqrtxz1*x*lomz*lspec1*NCi2*zi -  
   8*sqrtxz1*x*lz*NCi2*zi -  
   12*x*l2*lz*NCi2*zi -  
   16*sqrtxz1*x*l2*lz*NCi2*zi -  
   4*x*lx*lz*NCi2*zi +  
   2*sqrtxz1*x*lx*lz*NCi2*zi -  
   4*sqrtxz1*x*lomz*lz*NCi2*zi +  
   4*x*lspec1*lz*NCi2*zi +  
   8*sqrtxz1*x*lspec1*lz*NCi2*zi +  
   8*x*l2*lspec2*NCi2*zi +  
   8*sqrtxz1*x*l2*lspec2*NCi2*zi +  
   4*x*lx*lspec2*NCi2*zi +  
   4*sqrtxz1*x*lx*lspec2*NCi2*zi -  
   8*x*lspec1*lspec2*NCi2*zi -  
   8*sqrtxz1*x*lspec1*lspec2*NCi2* 
    zi + 8*x*lz*lspec2*NCi2*zi +  
   8*sqrtxz1*x*lz*lspec2*NCi2*zi +  
   2*x*lx*lxpz*NCi2*zi -  
   2*x*lz*lxpz*NCi2*zi +  
   2*x*lx*lopxzi*NCi2*zi -  
   2*x*lz*lopxzi*NCi2*zi +  
   4*sqrtxz1*x*lomz*lspec7*NCi2* 
    zi + (2*x*pi2*zi)/3. +  
   (2*sqrtxz1*x*pi2*zi)/3. -  
   (2*x*NCi2*pi2*zi)/3. -  
   (2*sqrtxz1*x*NCi2*pi2*zi)/3. -  
   4*Li2mx*xi*zi + 4*Li2x*xi*zi +  
   2*li2spec19*xi*zi +  
   2*li2spec20*xi*zi - 8*lx*xi*zi +  
   4*lomx*lx*xi*zi -  
   4*lx*lopx*xi*zi -  
   2*lx*lz*xi*zi +  
   2*lx*lopxz*xi*zi +  
   2*lz*lopxz*xi*zi +  
   2*lx*lopxzi*xi*zi -  
   2*lz*lopxzi*xi*zi +  
   4*Li2mx*NCi2*xi*zi -  
   4*Li2x*NCi2*xi*zi -  
   2*li2spec19*NCi2*xi*zi -  
   2*li2spec20*NCi2*xi*zi +  
   8*lx*NCi2*xi*zi -  
   4*lomx*lx*NCi2*xi*zi +  
   4*lx*lopx*NCi2*xi*zi +  
   2*lx*lz*NCi2*xi*zi -  
   2*lx*lopxz*NCi2*xi*zi -  
   2*lz*lopxz*NCi2*xi*zi -  
   2*lx*lopxzi*NCi2*xi*zi +  
   2*lz*lopxzi*NCi2*xi*zi -  
   (2*pi2*xi*zi)/3. +  
   (2*NCi2*pi2*xi*zi)/3. +  
   4*Li2mx*opxi*zi - 4*Li2x*opxi*zi -  
   2*li2spec19*opxi*zi -  
   2*li2spec20*opxi*zi +  
   8*lx*opxi*zi -  
   4*lomx*lx*opxi*zi +  
   4*lx*lopx*opxi*zi +  
   2*lx*lz*opxi*zi -  
   2*lx*lopxz*opxi*zi -  
   2*lz*lopxz*opxi*zi -  
   2*lx*lopxzi*opxi*zi +  
   2*lz*lopxzi*opxi*zi -  
   4*Li2mx*NCi2*opxi*zi +  
   4*Li2x*NCi2*opxi*zi +  
   2*li2spec19*NCi2*opxi*zi +  
   2*li2spec20*NCi2*opxi*zi -  
   8*lx*NCi2*opxi*zi +  
   4*lomx*lx*NCi2*opxi*zi -  
   4*lx*lopx*NCi2*opxi*zi -  
   2*lx*lz*NCi2*opxi*zi +  
   2*lx*lopxz*NCi2*opxi*zi +  
   2*lz*lopxz*NCi2*opxi*zi +  
   2*lx*lopxzi*NCi2*opxi*zi -  
   2*lz*lopxzi*NCi2*opxi*zi +  
   (2*pi2*opxi*zi)/3. -  
   (2*NCi2*pi2*opxi*zi)/3. +  
   4*Li2mx*xi*opxi*zi -  
   4*Li2x*xi*opxi*zi -  
   2*li2spec19*xi*opxi*zi -  
   2*li2spec20*xi*opxi*zi +  
   8*lx*xi*opxi*zi -  
   4*lomx*lx*xi*opxi*zi +  
   4*lx*lopx*xi*opxi*zi +  
   2*lx*lz*xi*opxi*zi -  
   2*lx*lopxz*xi*opxi*zi -  
   2*lz*lopxz*xi*opxi*zi -  
   2*lx*lopxzi*xi*opxi*zi +  
   2*lz*lopxzi*xi*opxi*zi -  
   4*Li2mx*NCi2*xi*opxi*zi +  
   4*Li2x*NCi2*xi*opxi*zi +  
   2*li2spec19*NCi2*xi*opxi*zi +  
   2*li2spec20*NCi2*xi*opxi*zi -  
   8*lx*NCi2*xi*opxi*zi +  
   4*lomx*lx*NCi2*xi*opxi*zi -  
   4*lx*lopx*NCi2*xi*opxi*zi -  
   2*lx*lz*NCi2*xi*opxi*zi +  
   2*lx*lopxz*NCi2*xi*opxi*zi +  
   2*lz*lopxz*NCi2*xi*opxi*zi +  
   2*lx*lopxzi*NCi2*xi*opxi* 
    zi - 2*lz*lopxzi*NCi2*xi* 
    opxi*zi + (2*pi2*xi*opxi*zi)/ 
    3. - (2*NCi2*pi2*xi*opxi*zi)/3. -  
   4*x*Li2mx*omzi*zi + 4*x*Li2x*omzi*zi +  
   2*x*li2spec19*omzi*zi +  
   2*x*li2spec20*omzi*zi -  
   8*x*lx*omzi*zi +  
   4*x*lomx*lx*omzi*zi -  
   4*x*lx*lopx*omzi*zi -  
   2*x*lx*lz*omzi*zi +  
   2*x*lx*lopxz*omzi*zi +  
   2*x*lz*lopxz*omzi*zi +  
   2*x*lx*lopxzi*omzi*zi -  
   2*x*lz*lopxzi*omzi*zi +  
   4*x*Li2mx*NCi2*omzi*zi -  
   4*x*Li2x*NCi2*omzi*zi -  
   2*x*li2spec19*NCi2*omzi*zi -  
   2*x*li2spec20*NCi2*omzi*zi +  
   8*x*lx*NCi2*omzi*zi -  
   4*x*lomx*lx*NCi2*omzi*zi +  
   4*x*lx*lopx*NCi2*omzi*zi +  
   2*x*lx*lz*NCi2*omzi*zi -  
   2*x*lx*lopxz*NCi2*omzi*zi -  
   2*x*lz*lopxz*NCi2*omzi*zi -  
   2*x*lx*lopxzi*NCi2*omzi*zi +  
   2*x*lz*lopxzi*NCi2*omzi*zi -  
   (2*x*pi2*omzi*zi)/3. +  
   (2*x*NCi2*pi2*omzi*zi)/3. +  
   4*Li2mx*xi*omzi*zi -  
   4*Li2x*xi*omzi*zi -  
   2*li2spec19*xi*omzi*zi -  
   2*li2spec20*xi*omzi*zi +  
   8*lx*xi*omzi*zi -  
   4*lomx*lx*xi*omzi*zi +  
   4*lx*lopx*xi*omzi*zi +  
   2*lx*lz*xi*omzi*zi -  
   2*lx*lopxz*xi*omzi*zi -  
   2*lz*lopxz*xi*omzi*zi -  
   2*lx*lopxzi*xi*omzi*zi +  
   2*lz*lopxzi*xi*omzi*zi -  
   4*Li2mx*NCi2*xi*omzi*zi +  
   4*Li2x*NCi2*xi*omzi*zi +  
   2*li2spec19*NCi2*xi*omzi*zi +  
   2*li2spec20*NCi2*xi*omzi*zi -  
   8*lx*NCi2*xi*omzi*zi +  
   4*lomx*lx*NCi2*xi*omzi*zi -  
   4*lx*lopx*NCi2*xi*omzi*zi -  
   2*lx*lz*NCi2*xi*omzi*zi +  
   2*lx*lopxz*NCi2*xi*omzi*zi +  
   2*lz*lopxz*NCi2*xi*omzi*zi +  
   2*lx*lopxzi*NCi2*xi*omzi* 
    zi - 2*lz*lopxzi*NCi2*xi* 
    omzi*zi + (2*pi2*xi*omzi*zi)/ 
    3. - (2*NCi2*pi2*xi*omzi*zi)/3. -  
   4*Li2mx*opxi*omzi*zi +  
   4*Li2x*opxi*omzi*zi +  
   2*li2spec19*opxi*omzi*zi +  
   2*li2spec20*opxi*omzi*zi -  
   8*lx*opxi*omzi*zi +  
   4*lomx*lx*opxi*omzi*zi -  
   4*lx*lopx*opxi*omzi*zi -  
   2*lx*lz*opxi*omzi*zi +  
   2*lx*lopxz*opxi*omzi*zi +  
   2*lz*lopxz*opxi*omzi*zi +  
   2*lx*lopxzi*opxi*omzi*zi -  
   2*lz*lopxzi*opxi*omzi*zi +  
   4*Li2mx*NCi2*opxi*omzi*zi -  
   4*Li2x*NCi2*opxi*omzi*zi -  
   2*li2spec19*NCi2*opxi*omzi*zi -  
   2*li2spec20*NCi2*opxi*omzi*zi +  
   8*lx*NCi2*opxi*omzi*zi -  
   4*lomx*lx*NCi2*opxi*omzi*zi +  
   4*lx*lopx*NCi2*opxi*omzi*zi +  
   2*lx*lz*NCi2*opxi*omzi*zi -  
   2*lx*lopxz*NCi2*opxi*omzi*zi -  
   2*lz*lopxz*NCi2*opxi*omzi*zi -  
   2*lx*lopxzi*NCi2*opxi*omzi* 
    zi + 2*lz*lopxzi*NCi2*opxi* 
    omzi*zi - (2*pi2*opxi*omzi* 
      zi)/3. + (2*NCi2*pi2*opxi*omzi* 
      zi)/3. - 4*Li2mx*xi*opxi*omzi* 
    zi + 4*Li2x*xi*opxi*omzi*zi +  
   2*li2spec19*xi*opxi*omzi*zi +  
   2*li2spec20*xi*opxi*omzi*zi -  
   8*lx*xi*opxi*omzi*zi +  
   4*lomx*lx*xi*opxi*omzi*zi -  
   4*lx*lopx*xi*opxi*omzi*zi -  
   2*lx*lz*xi*opxi*omzi*zi +  
   2*lx*lopxz*xi*opxi*omzi*zi +  
   2*lz*lopxz*xi*opxi*omzi*zi +  
   2*lx*lopxzi*xi*opxi*omzi* 
    zi - 2*lz*lopxzi*xi*opxi* 
    omzi*zi + 4*Li2mx*NCi2*xi*opxi* 
    omzi*zi - 4*Li2x*NCi2*xi*opxi* 
    omzi*zi - 2*li2spec19*NCi2*xi* 
    opxi*omzi*zi -  
   2*li2spec20*NCi2*xi*opxi*omzi* 
    zi + 8*lx*NCi2*xi*opxi*omzi* 
    zi - 4*lomx*lx*NCi2*xi*opxi* 
    omzi*zi + 4*lx*lopx*NCi2*xi* 
    opxi*omzi*zi +  
   2*lx*lz*NCi2*xi*opxi*omzi* 
    zi - 2*lx*lopxz*NCi2*xi*opxi* 
    omzi*zi - 2*lz*lopxz*NCi2*xi* 
    opxi*omzi*zi -  
   2*lx*lopxzi*NCi2*xi*opxi* 
    omzi*zi + 2*lz*lopxzi*NCi2* 
    xi*opxi*omzi*zi -  
   (2*pi2*xi*opxi*omzi*zi)/3. +  
   (2*NCi2*pi2*xi*opxi*omzi*zi)/ 
    3. + NC*z2 - (7*NC*x*z2)/2. -  
   4*sqrtxz3*itani1*z2 -  
   (23*NC*sqrtxz3*itani1*z2)/2. +  
   4*sqrtxz3*itani2*z2 +  
   (23*NC*sqrtxz3*itani2*z2)/2. -  
   8*sqrtxz3*itani3*z2 -  
   23*NC*sqrtxz3*itani3*z2 + 8*NC*x*Li2x*z2 -  
   8*NC*x*li2spec1*z2 +  
   8*NC*x*li2spec2*z2 +  
   8*NC*x*li2spec19*z2 -  
   8*sqrtxz3*atanspec1*lspec3*z2 -  
   23*NC*sqrtxz3*atanspec1*lspec3*z2 +  
   NC*x*lomx*z2 - 2*NC*x*lx*z2 +  
   16*NC*x*l2*lx*z2 + 8*NC*x*lomx*lx*z2 +  
   NC*x*lomz*z2 - 24*NC*x*l2*lspec1*z2 -  
   8*NC*x*lx*lspec1*z2 + NC*x*lz*z2 +  
   16*NC*x*l2*lz*z2 + 4*NC*x*lx*lz*z2 -  
   8*NC*x*lspec1*lz*z2 +  
   8*sqrtxz3*atanspec2*lspec4*z2 +  
   23*NC*sqrtxz3*atanspec2*lspec4*z2 -  
   8*NC*x*l2*lspec2*z2 -  
   8*NC*x*lx*lspec2*z2 +  
   8*NC*x*lspec1*lspec2*z2 -  
   8*NC*x*lz*lspec2*z2 -  
   4*NC*x*lx*lxpz*z2 + 4*NC*x*lz*lxpz*z2 +  
   8*NC*x*lx*lopxz*z2 +  
   8*NC*x*lz*lopxz*z2 +  
   4*NC*x*lx*lopxzi*z2 -  
   4*NC*x*lz*lopxzi*z2 +  
   4*sqrtxz3*itani1*NCi2*z2 -  
   4*sqrtxz3*itani2*NCi2*z2 +  
   8*sqrtxz3*itani3*NCi2*z2 +  
   8*sqrtxz3*atanspec1*lspec3*NCi2*z2 -  
   8*sqrtxz3*atanspec2*lspec4*NCi2*z2 -  
   NCi*z2 + (7*x*NCi*z2)/2. +  
   (23*sqrtxz3*itani1*NCi*z2)/2. -  
   (23*sqrtxz3*itani2*NCi*z2)/2. +  
   23*sqrtxz3*itani3*NCi*z2 -  
   8*x*Li2x*NCi*z2 +  
   8*x*li2spec1*NCi*z2 -  
   8*x*li2spec2*NCi*z2 -  
   8*x*li2spec19*NCi*z2 +  
   23*sqrtxz3*atanspec1*lspec3*NCi*z2 -  
   x*lomx*NCi*z2 + 2*x*lx*NCi*z2 -  
   16*x*l2*lx*NCi*z2 -  
   8*x*lomx*lx*NCi*z2 -  
   x*lomz*NCi*z2 +  
   24*x*l2*lspec1*NCi*z2 +  
   8*x*lx*lspec1*NCi*z2 -  
   x*lz*NCi*z2 - 16*x*l2*lz*NCi*z2 -  
   4*x*lx*lz*NCi*z2 +  
   8*x*lspec1*lz*NCi*z2 -  
   23*sqrtxz3*atanspec2*lspec4*NCi*z2 +  
   8*x*l2*lspec2*NCi*z2 +  
   8*x*lx*lspec2*NCi*z2 -  
   8*x*lspec1*lspec2*NCi*z2 +  
   8*x*lz*lspec2*NCi*z2 +  
   4*x*lx*lxpz*NCi*z2 -  
   4*x*lz*lxpz*NCi*z2 -  
   8*x*lx*lopxz*NCi*z2 -  
   8*x*lz*lopxz*NCi*z2 -  
   4*x*lx*lopxzi*NCi*z2 +  
   4*x*lz*lopxzi*NCi*z2 -  
   (4*NC*x*pi2*z2)/3. + (4*x*NCi*pi2*z2)/3. +  
   4*NC*x*li2spec5*sqrtxz2i*z2 -  
   4*NC*x*li2spec6*sqrtxz2i*z2 -  
   4*NC*x*li2spec7*sqrtxz2i* 
    z2 + 4*NC*x*li2spec8* 
    sqrtxz2i*z2 - 4*NC*x*lx*lspec5* 
    sqrtxz2i*z2 + 4*NC*x*lx*lspec6* 
    sqrtxz2i*z2 - 4*x*li2spec5*NCi* 
    sqrtxz2i*z2 + 4*x*li2spec6*NCi* 
    sqrtxz2i*z2 + 4*x* 
    li2spec7*NCi* 
    sqrtxz2i*z2 - 4*x* 
    li2spec8*NCi* 
    sqrtxz2i*z2 + 4*x*lx*lspec5*NCi* 
    sqrtxz2i*z2 - 4*x*lx*lspec6*NCi* 
    sqrtxz2i*z2 + 8*NC*Li2mx*xi*z2 -  
   8*NC*Li2x*xi*z2 - 4*NC*li2spec19*xi*z2 -  
   4*NC*li2spec20*xi*z2 +  
   16*NC*lx*xi*z2 -  
   8*NC*lomx*lx*xi*z2 +  
   8*NC*lx*lopx*xi*z2 +  
   4*NC*lx*lz*xi*z2 -  
   4*NC*lx*lopxz*xi*z2 -  
   4*NC*lz*lopxz*xi*z2 -  
   4*NC*lx*lopxzi*xi*z2 +  
   4*NC*lz*lopxzi*xi*z2 -  
   8*Li2mx*NCi*xi*z2 +  
   8*Li2x*NCi*xi*z2 +  
   4*li2spec19*NCi*xi*z2 +  
   4*li2spec20*NCi*xi*z2 -  
   16*lx*NCi*xi*z2 +  
   8*lomx*lx*NCi*xi*z2 -  
   8*lx*lopx*NCi*xi*z2 -  
   4*lx*lz*NCi*xi*z2 +  
   4*lx*lopxz*NCi*xi*z2 +  
   4*lz*lopxz*NCi*xi*z2 +  
   4*lx*lopxzi*NCi*xi*z2 -  
   4*lz*lopxzi*NCi*xi*z2 +  
   (4*NC*pi2*xi*z2)/3. -  
   (4*NCi*pi2*xi*z2)/3. -  
   8*NC*Li2mx*opxi*z2 + 8*NC*Li2x*opxi*z2 +  
   4*NC*li2spec19*opxi*z2 +  
   4*NC*li2spec20*opxi*z2 -  
   16*NC*lx*opxi*z2 +  
   8*NC*lomx*lx*opxi*z2 -  
   8*NC*lx*lopx*opxi*z2 -  
   4*NC*lx*lz*opxi*z2 +  
   4*NC*lx*lopxz*opxi*z2 +  
   4*NC*lz*lopxz*opxi*z2 +  
   4*NC*lx*lopxzi*opxi*z2 -  
   4*NC*lz*lopxzi*opxi*z2 +  
   8*Li2mx*NCi*opxi*z2 -  
   8*Li2x*NCi*opxi*z2 -  
   4*li2spec19*NCi*opxi*z2 -  
   4*li2spec20*NCi*opxi*z2 +  
   16*lx*NCi*opxi*z2 -  
   8*lomx*lx*NCi*opxi*z2 +  
   8*lx*lopx*NCi*opxi*z2 +  
   4*lx*lz*NCi*opxi*z2 -  
   4*lx*lopxz*NCi*opxi*z2 -  
   4*lz*lopxz*NCi*opxi*z2 -  
   4*lx*lopxzi*NCi*opxi*z2 +  
   4*lz*lopxzi*NCi*opxi*z2 -  
   (4*NC*pi2*opxi*z2)/3. +  
   (4*NCi*pi2*opxi*z2)/3. -  
   8*NC*Li2mx*xi*opxi*z2 +  
   8*NC*Li2x*xi*opxi*z2 +  
   4*NC*li2spec19*xi*opxi*z2 +  
   4*NC*li2spec20*xi*opxi*z2 -  
   16*NC*lx*xi*opxi*z2 +  
   8*NC*lomx*lx*xi*opxi*z2 -  
   8*NC*lx*lopx*xi*opxi*z2 -  
   4*NC*lx*lz*xi*opxi*z2 +  
   4*NC*lx*lopxz*xi*opxi*z2 +  
   4*NC*lz*lopxz*xi*opxi*z2 +  
   4*NC*lx*lopxzi*xi*opxi*z2 -  
   4*NC*lz*lopxzi*xi*opxi*z2 +  
   8*Li2mx*NCi*xi*opxi*z2 -  
   8*Li2x*NCi*xi*opxi*z2 -  
   4*li2spec19*NCi*xi*opxi*z2 -  
   4*li2spec20*NCi*xi*opxi*z2 +  
   16*lx*NCi*xi*opxi*z2 -  
   8*lomx*lx*NCi*xi*opxi*z2 +  
   8*lx*lopx*NCi*xi*opxi*z2 +  
   4*lx*lz*NCi*xi*opxi*z2 -  
   4*lx*lopxz*NCi*xi*opxi*z2 -  
   4*lz*lopxz*NCi*xi*opxi*z2 -  
   4*lx*lopxzi*NCi*xi*opxi* 
    z2 + 4*lz*lopxzi*NCi*xi* 
    opxi*z2 - (4*NC*pi2*xi*opxi* 
      z2)/3. + (4*NCi*pi2*xi*opxi* 
      z2)/3. + 12*x*li2spec1*sqrtxz1i* 
    opzi - 4*x*z*li2spec1*sqrtxz1i* 
    opzi + 12*x*li2spec2*sqrtxz1i* 
    opzi - 4*x*z*li2spec2*sqrtxz1i* 
    opzi - 24*x*li2spec21*sqrtxz1i* 
    opzi + 8*x*z*li2spec21*sqrtxz1i* 
    opzi - 12*x*li2spec22* 
    sqrtxz1i*opzi +  
   4*x*z*li2spec22*sqrtxz1i*opzi -  
   24*x*li2spec23*sqrtxz1i*opzi +  
   8*x*z*li2spec23*sqrtxz1i*opzi -  
   12*x*li2spec3*sqrtxz1i* 
    opzi + 4*x*z*li2spec3* 
    sqrtxz1i*opzi -  
   12*x*li2spec4*sqrtxz1i* 
    opzi + 4*x*z*li2spec4* 
    sqrtxz1i*opzi -  
   48*x*l2*sqrtxz1i*opzi +  
   16*x*z*l2*sqrtxz1i*opzi -  
   24*x*lx*sqrtxz1i*opzi +  
   8*x*z*lx*sqrtxz1i*opzi -  
   24*x*l2*lomz*sqrtxz1i*opzi +  
   8*x*z*l2*lomz*sqrtxz1i*opzi -  
   12*x*lx*lomz*sqrtxz1i*opzi +  
   4*x*z*lx*lomz*sqrtxz1i*opzi +  
   48*x*lspec1*sqrtxz1i*opzi -  
   16*x*z*lspec1*sqrtxz1i*opzi +  
   48*x*l2*lspec1*sqrtxz1i*opzi -  
   16*x*z*l2*lspec1*sqrtxz1i*opzi -  
   12*x*lx*lspec1*sqrtxz1i*opzi +  
   4*x*z*lx*lspec1*sqrtxz1i*opzi +  
   24*x*lomz*lspec1*sqrtxz1i*opzi -  
   8*x*z*lomz*lspec1*sqrtxz1i*opzi -  
   24*x*lz*sqrtxz1i*opzi +  
   8*x*z*lz*sqrtxz1i*opzi -  
   48*x*l2*lz*sqrtxz1i*opzi +  
   16*x*z*l2*lz*sqrtxz1i*opzi +  
   6*x*lx*lz*sqrtxz1i*opzi -  
   2*x*z*lx*lz*sqrtxz1i*opzi -  
   12*x*lomz*lz*sqrtxz1i*opzi +  
   4*x*z*lomz*lz*sqrtxz1i*opzi +  
   24*x*lspec1*lz*sqrtxz1i*opzi -  
   8*x*z*lspec1*lz*sqrtxz1i*opzi +  
   24*x*l2*lspec2*sqrtxz1i*opzi -  
   8*x*z*l2*lspec2*sqrtxz1i*opzi +  
   12*x*lx*lspec2*sqrtxz1i*opzi -  
   4*x*z*lx*lspec2*sqrtxz1i*opzi -  
   24*x*lspec1*lspec2*sqrtxz1i* 
    opzi + 8*x*z*lspec1*lspec2* 
    sqrtxz1i*opzi +  
   24*x*lz*lspec2*sqrtxz1i*opzi -  
   8*x*z*lz*lspec2*sqrtxz1i*opzi +  
   12*x*lomz*lspec7*sqrtxz1i* 
    opzi - 4*x*z*lomz*lspec7* 
    sqrtxz1i*opzi -  
   12*x*li2spec1*NCi2*sqrtxz1i* 
    opzi + 4*x*z*li2spec1*NCi2* 
    sqrtxz1i*opzi -  
   12*x*li2spec2*NCi2*sqrtxz1i* 
    opzi + 4*x*z*li2spec2*NCi2* 
    sqrtxz1i*opzi +  
   24*x*li2spec21*NCi2*sqrtxz1i* 
    opzi - 8*x*z*li2spec21*NCi2* 
    sqrtxz1i*opzi +  
   12*x*li2spec22*NCi2*sqrtxz1i*opzi -  
   4*x*z*li2spec22*NCi2*sqrtxz1i*opzi +  
   24*x*li2spec23*NCi2*sqrtxz1i* 
    opzi - 8*x*z*li2spec23*NCi2* 
    sqrtxz1i*opzi +  
   12*x*li2spec3*NCi2* 
    sqrtxz1i*opzi -  
   4*x*z*li2spec3*NCi2* 
    sqrtxz1i*opzi +  
   12*x*li2spec4*NCi2* 
    sqrtxz1i*opzi -  
   4*x*z*li2spec4*NCi2* 
    sqrtxz1i*opzi +  
   48*x*l2*NCi2*sqrtxz1i*opzi -  
   16*x*z*l2*NCi2*sqrtxz1i*opzi +  
   24*x*lx*NCi2*sqrtxz1i*opzi -  
   8*x*z*lx*NCi2*sqrtxz1i*opzi +  
   24*x*l2*lomz*NCi2*sqrtxz1i*opzi -  
   8*x*z*l2*lomz*NCi2*sqrtxz1i*opzi +  
   12*x*lx*lomz*NCi2*sqrtxz1i*opzi -  
   4*x*z*lx*lomz*NCi2*sqrtxz1i*opzi -  
   48*x*lspec1*NCi2*sqrtxz1i*opzi +  
   16*x*z*lspec1*NCi2*sqrtxz1i*opzi -  
   48*x*l2*lspec1*NCi2*sqrtxz1i* 
    opzi + 16*x*z*l2*lspec1*NCi2* 
    sqrtxz1i*opzi +  
   12*x*lx*lspec1*NCi2*sqrtxz1i* 
    opzi - 4*x*z*lx*lspec1*NCi2* 
    sqrtxz1i*opzi -  
   24*x*lomz*lspec1*NCi2*sqrtxz1i* 
    opzi + 8*x*z*lomz*lspec1*NCi2* 
    sqrtxz1i*opzi +  
   24*x*lz*NCi2*sqrtxz1i*opzi -  
   8*x*z*lz*NCi2*sqrtxz1i*opzi +  
   48*x*l2*lz*NCi2*sqrtxz1i*opzi -  
   16*x*z*l2*lz*NCi2*sqrtxz1i*opzi -  
   6*x*lx*lz*NCi2*sqrtxz1i*opzi +  
   2*x*z*lx*lz*NCi2*sqrtxz1i*opzi +  
   12*x*lomz*lz*NCi2*sqrtxz1i*opzi -  
   4*x*z*lomz*lz*NCi2*sqrtxz1i*opzi -  
   24*x*lspec1*lz*NCi2*sqrtxz1i* 
    opzi + 8*x*z*lspec1*lz*NCi2* 
    sqrtxz1i*opzi -  
   24*x*l2*lspec2*NCi2*sqrtxz1i* 
    opzi + 8*x*z*l2*lspec2*NCi2* 
    sqrtxz1i*opzi -  
   12*x*lx*lspec2*NCi2*sqrtxz1i* 
    opzi + 4*x*z*lx*lspec2*NCi2* 
    sqrtxz1i*opzi +  
   24*x*lspec1*lspec2*NCi2*sqrtxz1i* 
    opzi - 8*x*z*lspec1*lspec2* 
    NCi2*sqrtxz1i*opzi -  
   24*x*lz*lspec2*NCi2*sqrtxz1i* 
    opzi + 8*x*z*lz*lspec2*NCi2* 
    sqrtxz1i*opzi -  
   12*x*lomz*lspec7*NCi2*sqrtxz1i* 
    opzi + 4*x*z*lomz*lspec7* 
    NCi2*sqrtxz1i*opzi -  
   2*x*pi2*sqrtxz1i*opzi +  
   (2*x*z*pi2*sqrtxz1i*opzi)/3. +  
   2*x*NCi2*pi2*sqrtxz1i*opzi -  
   (2*x*z*NCi2*pi2*sqrtxz1i*opzi)/3. -  
   16*li2spec1*sqrtxz1i*x2*opzi -  
   16*li2spec2*sqrtxz1i*x2*opzi +  
   32*li2spec21*sqrtxz1i*x2*opzi +  
   16*li2spec22*sqrtxz1i*x2*opzi +  
   32*li2spec23*sqrtxz1i*x2*opzi +  
   16*li2spec3*sqrtxz1i* 
    x2*opzi + 16* 
    li2spec4*sqrtxz1i*x2* 
    opzi + 64*l2*sqrtxz1i*x2*opzi +  
   32*lx*sqrtxz1i*x2*opzi +  
   32*l2*lomz*sqrtxz1i*x2*opzi +  
   16*lx*lomz*sqrtxz1i*x2*opzi -  
   64*lspec1*sqrtxz1i*x2*opzi -  
   64*l2*lspec1*sqrtxz1i*x2*opzi +  
   16*lx*lspec1*sqrtxz1i*x2*opzi -  
   32*lomz*lspec1*sqrtxz1i*x2* 
    opzi + 32*lz*sqrtxz1i*x2*opzi +  
   64*l2*lz*sqrtxz1i*x2*opzi -  
   8*lx*lz*sqrtxz1i*x2*opzi +  
   16*lomz*lz*sqrtxz1i*x2*opzi -  
   32*lspec1*lz*sqrtxz1i*x2*opzi -  
   32*l2*lspec2*sqrtxz1i*x2*opzi -  
   16*lx*lspec2*sqrtxz1i*x2*opzi +  
   32*lspec1*lspec2*sqrtxz1i*x2* 
    opzi - 32*lz*lspec2*sqrtxz1i*x2* 
    opzi - 16*lomz*lspec7* 
    sqrtxz1i*x2*opzi +  
   16*li2spec1*NCi2*sqrtxz1i*x2* 
    opzi + 16*li2spec2*NCi2* 
    sqrtxz1i*x2*opzi -  
   32*li2spec21*NCi2*sqrtxz1i*x2* 
    opzi - 16*li2spec22*NCi2* 
    sqrtxz1i*x2*opzi -  
   32*li2spec23*NCi2*sqrtxz1i*x2* 
    opzi - 16*li2spec3* 
    NCi2*sqrtxz1i*x2*opzi -  
   16*li2spec4*NCi2* 
    sqrtxz1i*x2*opzi -  
   64*l2*NCi2*sqrtxz1i*x2*opzi -  
   32*lx*NCi2*sqrtxz1i*x2*opzi -  
   32*l2*lomz*NCi2*sqrtxz1i*x2*opzi -  
   16*lx*lomz*NCi2*sqrtxz1i*x2*opzi +  
   64*lspec1*NCi2*sqrtxz1i*x2* 
    opzi + 64*l2*lspec1*NCi2*sqrtxz1i* 
    x2*opzi - 16*lx*lspec1*NCi2* 
    sqrtxz1i*x2*opzi +  
   32*lomz*lspec1*NCi2*sqrtxz1i*x2* 
    opzi - 32*lz*NCi2*sqrtxz1i*x2* 
    opzi - 64*l2*lz*NCi2*sqrtxz1i*x2* 
    opzi + 8*lx*lz*NCi2*sqrtxz1i*x2* 
    opzi - 16*lomz*lz*NCi2*sqrtxz1i*x2* 
    opzi + 32*lspec1*lz*NCi2*sqrtxz1i* 
    x2*opzi + 32*l2*lspec2*NCi2* 
    sqrtxz1i*x2*opzi +  
   16*lx*lspec2*NCi2*sqrtxz1i*x2* 
    opzi - 32*lspec1*lspec2*NCi2* 
    sqrtxz1i*x2*opzi +  
   32*lz*lspec2*NCi2*sqrtxz1i*x2* 
    opzi + 16*lomz*lspec7*NCi2* 
    sqrtxz1i*x2*opzi +  
   (8*pi2*sqrtxz1i*x2*opzi)/3. -  
   (8*NCi2*pi2*sqrtxz1i*x2*opzi)/3. +  
   4*x*li2spec1*zi*opzi +  
   4*sqrtxz1*x*li2spec1*zi*opzi -  
   4*x*li2spec2*zi*opzi +  
   4*sqrtxz1*x*li2spec2*zi*opzi -  
   2*x*li2spec19*zi*opzi -  
   8*sqrtxz1*x*li2spec21*zi*opzi -  
   4*sqrtxz1*x*li2spec22*zi* 
    opzi - 8*sqrtxz1*x*li2spec23*zi* 
    opzi + 2*x*li2spec20*zi*opzi -  
   4*x*li2spec3*zi* 
    opzi - 4*sqrtxz1*x* 
    li2spec3*zi*opzi 
    + 4*x*li2spec4*zi* 
    opzi - 4*sqrtxz1*x* 
    li2spec4*zi*opzi 
    - 16*sqrtxz1*x*l2*zi*opzi -  
   8*sqrtxz1*x*lx*zi*opzi -  
   4*x*l2*lx*zi*opzi -  
   8*sqrtxz1*x*l2*lomz*zi*opzi -  
   4*sqrtxz1*x*lx*lomz*zi*opzi +  
   16*sqrtxz1*x*lspec1*zi*opzi +  
   8*x*l2*lspec1*zi*opzi +  
   16*sqrtxz1*x*l2*lspec1*zi*opzi -  
   4*sqrtxz1*x*lx*lspec1*zi*opzi +  
   8*sqrtxz1*x*lomz*lspec1*zi*opzi -  
   8*sqrtxz1*x*lz*zi*opzi -  
   12*x*l2*lz*zi*opzi -  
   16*sqrtxz1*x*l2*lz*zi*opzi -  
   2*x*lx*lz*zi*opzi +  
   2*sqrtxz1*x*lx*lz*zi*opzi -  
   4*sqrtxz1*x*lomz*lz*zi*opzi +  
   4*x*lspec1*lz*zi*opzi +  
   8*sqrtxz1*x*lspec1*lz*zi*opzi +  
   8*x*l2*lspec2*zi*opzi +  
   8*sqrtxz1*x*l2*lspec2*zi*opzi +  
   4*x*lx*lspec2*zi*opzi +  
   4*sqrtxz1*x*lx*lspec2*zi*opzi -  
   8*x*lspec1*lspec2*zi*opzi -  
   8*sqrtxz1*x*lspec1*lspec2*zi* 
    opzi + 8*x*lz*lspec2*zi*opzi +  
   8*sqrtxz1*x*lz*lspec2*zi*opzi +  
   2*x*lx*lxpz*zi*opzi -  
   2*x*lz*lxpz*zi*opzi -  
   2*x*lx*lopxz*zi*opzi -  
   2*x*lz*lopxz*zi*opzi +  
   4*sqrtxz1*x*lomz*lspec7*zi* 
    opzi - 4*x*li2spec1*NCi2*zi* 
    opzi - 4*sqrtxz1*x*li2spec1*NCi2* 
    zi*opzi + 4*x*li2spec2*NCi2* 
    zi*opzi - 4*sqrtxz1*x*li2spec2* 
    NCi2*zi*opzi +  
   2*x*li2spec19*NCi2*zi*opzi +  
   8*sqrtxz1*x*li2spec21*NCi2*zi* 
    opzi + 4*sqrtxz1*x* 
    li2spec22*NCi2*zi*opzi +  
   8*sqrtxz1*x*li2spec23*NCi2*zi* 
    opzi - 2*x*li2spec20*NCi2*zi* 
    opzi + 4*x*li2spec3* 
    NCi2*zi*opzi +  
   4*sqrtxz1*x*li2spec3*NCi2* 
    zi*opzi - 4*x* 
    li2spec4*NCi2*zi* 
    opzi + 4*sqrtxz1*x* 
    li2spec4*NCi2*zi* 
    opzi + 16*sqrtxz1*x*l2*NCi2*zi*opzi +  
   8*sqrtxz1*x*lx*NCi2*zi*opzi +  
   4*x*l2*lx*NCi2*zi*opzi +  
   8*sqrtxz1*x*l2*lomz*NCi2*zi*opzi +  
   4*sqrtxz1*x*lx*lomz*NCi2*zi*opzi -  
   16*sqrtxz1*x*lspec1*NCi2*zi*opzi -  
   8*x*l2*lspec1*NCi2*zi*opzi -  
   16*sqrtxz1*x*l2*lspec1*NCi2*zi* 
    opzi + 4*sqrtxz1*x*lx*lspec1*NCi2* 
    zi*opzi - 8*sqrtxz1*x*lomz*lspec1* 
    NCi2*zi*opzi +  
   8*sqrtxz1*x*lz*NCi2*zi*opzi +  
   12*x*l2*lz*NCi2*zi*opzi +  
   16*sqrtxz1*x*l2*lz*NCi2*zi*opzi +  
   2*x*lx*lz*NCi2*zi*opzi -  
   2*sqrtxz1*x*lx*lz*NCi2*zi*opzi +  
   4*sqrtxz1*x*lomz*lz*NCi2*zi*opzi -  
   4*x*lspec1*lz*NCi2*zi*opzi -  
   8*sqrtxz1*x*lspec1*lz*NCi2*zi* 
    opzi - 8*x*l2*lspec2*NCi2*zi* 
    opzi - 8*sqrtxz1*x*l2*lspec2*NCi2* 
    zi*opzi - 4*x*lx*lspec2*NCi2* 
    zi*opzi - 4*sqrtxz1*x*lx*lspec2* 
    NCi2*zi*opzi +  
   8*x*lspec1*lspec2*NCi2*zi* 
    opzi + 8*sqrtxz1*x*lspec1*lspec2* 
    NCi2*zi*opzi -  
   8*x*lz*lspec2*NCi2*zi*opzi -  
   8*sqrtxz1*x*lz*lspec2*NCi2*zi* 
    opzi - 2*x*lx*lxpz*NCi2*zi* 
    opzi + 2*x*lz*lxpz*NCi2*zi* 
    opzi + 2*x*lx*lopxz*NCi2*zi* 
    opzi + 2*x*lz*lopxz*NCi2*zi* 
    opzi - 4*sqrtxz1*x*lomz*lspec7* 
    NCi2*zi*opzi -  
   (2*sqrtxz1*x*pi2*zi*opzi)/3. +  
   (2*sqrtxz1*x*NCi2*pi2*zi*opzi)/3. -  
   8*x*l22 + 8*x*z*l22 + 8*x*NCi2*l22 -  
   8*x*z*NCi2*l22 + 24*x*sqrtxz1i*l22 -  
   12*x*z*sqrtxz1i*l22 -  
   24*x*NCi2*sqrtxz1i*l22 +  
   12*x*z*NCi2*sqrtxz1i*l22 -  
   48*sqrtxz1i*x2*l22 +  
   48*NCi2*sqrtxz1i*x2*l22 +  
   8*x*zi*l22 + 12*sqrtxz1*x*zi*l22 -  
   8*x*NCi2*zi*l22 -  
   12*sqrtxz1*x*NCi2*zi*l22 +  
   16*NC*x*z2*l22 - 16*x*NCi*z2*l22 -  
   36*x*sqrtxz1i*opzi*l22 +  
   12*x*z*sqrtxz1i*opzi*l22 +  
   36*x*NCi2*sqrtxz1i*opzi*l22 -  
   12*x*z*NCi2*sqrtxz1i*opzi*l22 +  
   48*sqrtxz1i*x2*opzi*l22 -  
   48*NCi2*sqrtxz1i*x2*opzi*l22 -  
   8*x*zi*opzi*l22 -  
   12*sqrtxz1*x*zi*opzi*l22 +  
   8*x*NCi2*zi*opzi*l22 +  
   12*sqrtxz1*x*NCi2*zi*opzi*l22 +  
   x*lx2 + 6*NC*x*z*lx2 - x*NCi2*lx2 -  
   6*x*z*NCi*lx2 - 6*x*sqrtxz1i*lx2 +  
   3*x*z*sqrtxz1i*lx2 +  
   6*x*NCi2*sqrtxz1i*lx2 -  
   3*x*z*NCi2*sqrtxz1i*lx2 +  
   3*xi*lx2 + 3*z*xi*lx2 +  
   3*NC*z*xi*lx2 - 3*NCi2*xi*lx2 -  
   3*z*NCi2*xi*lx2 -  
   3*z*NCi*xi*lx2 +  
   12*sqrtxz1i*x2*lx2 -  
   12*NCi2*sqrtxz1i*x2*lx2 -  
   3*opxi*lx2 - 3*z*opxi*lx2 -  
   3*NC*z*opxi*lx2 +  
   3*NCi2*opxi*lx2 +  
   3*z*NCi2*opxi*lx2 +  
   3*z*NCi*opxi*lx2 -  
   3*xi*opxi*lx2 -  
   3*z*xi*opxi*lx2 -  
   3*NC*z*xi*opxi*lx2 +  
   3*NCi2*xi*opxi*lx2 +  
   3*z*NCi2*xi*opxi*lx2 +  
   3*z*NCi*xi*opxi*lx2 -  
   3*x*omzi*lx2 +  
   3*x*NCi2*omzi*lx2 - 3*x*zi*lx2 -  
   3*sqrtxz1*x*zi*lx2 +  
   3*x*NCi2*zi*lx2 +  
   3*sqrtxz1*x*NCi2*zi*lx2 +  
   3*xi*zi*lx2 -  
   3*NCi2*xi*zi*lx2 -  
   3*opxi*zi*lx2 +  
   3*NCi2*opxi*zi*lx2 -  
   3*xi*opxi*zi*lx2 +  
   3*NCi2*xi*opxi*zi*lx2 +  
   3*x*omzi*zi*lx2 -  
   3*x*NCi2*omzi*zi*lx2 -  
   3*xi*omzi*zi*lx2 +  
   3*NCi2*xi*omzi*zi*lx2 +  
   3*opxi*omzi*zi*lx2 -  
   3*NCi2*opxi*omzi*zi*lx2 +  
   3*xi*opxi*omzi*zi*lx2 -  
   3*NCi2*xi*opxi*omzi*zi* 
    lx2 - 4*NC*x*z2*lx2 +  
   4*x*NCi*z2*lx2 -  
   6*NC*xi*z2*lx2 +  
   6*NCi*xi*z2*lx2 +  
   6*NC*opxi*z2*lx2 -  
   6*NCi*opxi*z2*lx2 +  
   6*NC*xi*opxi*z2*lx2 -  
   6*NCi*xi*opxi*z2*lx2 +  
   9*x*sqrtxz1i*opzi*lx2 -  
   3*x*z*sqrtxz1i*opzi*lx2 -  
   9*x*NCi2*sqrtxz1i*opzi*lx2 +  
   3*x*z*NCi2*sqrtxz1i*opzi*lx2 -  
   12*sqrtxz1i*x2*opzi*lx2 +  
   12*NCi2*sqrtxz1i*x2*opzi*lx2 +  
   3*sqrtxz1*x*zi*opzi*lx2 -  
   3*sqrtxz1*x*NCi2*zi*opzi*lx2 +  
   8*x*sqrtxz1i*lomz2 -  
   4*x*z*sqrtxz1i*lomz2 -  
   8*x*NCi2*sqrtxz1i*lomz2 +  
   4*x*z*NCi2*sqrtxz1i*lomz2 -  
   16*sqrtxz1i*x2*lomz2 +  
   16*NCi2*sqrtxz1i*x2*lomz2 +  
   4*sqrtxz1*x*zi*lomz2 -  
   4*sqrtxz1*x*NCi2*zi*lomz2 -  
   12*x*sqrtxz1i*opzi*lomz2 +  
   4*x*z*sqrtxz1i*opzi*lomz2 +  
   12*x*NCi2*sqrtxz1i*opzi*lomz2 -  
   4*x*z*NCi2*sqrtxz1i*opzi*lomz2 +  
   16*sqrtxz1i*x2*opzi*lomz2 -  
   16*NCi2*sqrtxz1i*x2*opzi*lomz2 -  
   4*sqrtxz1*x*zi*opzi*lomz2 +  
   4*sqrtxz1*x*NCi2*zi*opzi*lomz2 -  
   4*NC*x*z*lspec1_2 +  
   4*x*z*NCi*lspec1_2 +  
   8*x*sqrtxz1i*lspec1_2 -  
   4*x*z*sqrtxz1i*lspec1_2 -  
   8*x*NCi2*sqrtxz1i*lspec1_2 +  
   4*x*z*NCi2*sqrtxz1i*lspec1_2 -  
   16*sqrtxz1i*x2*lspec1_2 +  
   16*NCi2*sqrtxz1i*x2*lspec1_2 +  
   4*sqrtxz1*x*zi*lspec1_2 -  
   4*sqrtxz1*x*NCi2*zi*lspec1_2 +  
   8*NC*x*z2*lspec1_2 -  
   8*x*NCi*z2*lspec1_2 -  
   12*x*sqrtxz1i*opzi*lspec1_2 +  
   4*x*z*sqrtxz1i*opzi*lspec1_2 +  
   12*x*NCi2*sqrtxz1i*opzi* 
    lspec1_2 -  
   4*x*z*NCi2*sqrtxz1i*opzi* 
    lspec1_2 +  
   16*sqrtxz1i*x2*opzi*lspec1_2 -  
   16*NCi2*sqrtxz1i*x2*opzi* 
    lspec1_2 -  
   4*sqrtxz1*x*zi*opzi*lspec1_2 +  
   4*sqrtxz1*x*NCi2*zi*opzi* 
    lspec1_2 - x*lz2 - 2*NC*x*lz2 +  
   3*x*z*lz2 + 2*NC*x*z*lz2 +  
   x*NCi2*lz2 - 3*x*z*NCi2*lz2 +  
   2*x*NCi*lz2 - 2*x*z*NCi*lz2 +  
   10*x*sqrtxz1i*lz2 - 5*x*z*sqrtxz1i*lz2 -  
   10*x*NCi2*sqrtxz1i*lz2 +  
   5*x*z*NCi2*sqrtxz1i*lz2 - xi*lz2 -  
   z*xi*lz2 - NC*z*xi*lz2 +  
   NCi2*xi*lz2 +  
   z*NCi2*xi*lz2 +  
   z*NCi*xi*lz2 -  
   20*sqrtxz1i*x2*lz2 +  
   20*NCi2*sqrtxz1i*x2*lz2 +  
   opxi*lz2 + z*opxi*lz2 +  
   NC*z*opxi*lz2 -  
   NCi2*opxi*lz2 -  
   z*NCi2*opxi*lz2 -  
   z*NCi*opxi*lz2 +  
   xi*opxi*lz2 +  
   z*xi*opxi*lz2 +  
   NC*z*xi*opxi*lz2 -  
   NCi2*xi*opxi*lz2 -  
   z*NCi2*xi*opxi*lz2 -  
   z*NCi*xi*opxi*lz2 +  
   3*x*zi*lz2 + 5*sqrtxz1*x*zi*lz2 -  
   3*x*NCi2*zi*lz2 -  
   5*sqrtxz1*x*NCi2*zi*lz2 -  
   xi*zi*lz2 +  
   NCi2*xi*zi*lz2 +  
   opxi*zi*lz2 -  
   NCi2*opxi*zi*lz2 +  
   xi*opxi*zi*lz2 -  
   NCi2*xi*opxi*zi*lz2 -  
   x*omzi*zi*lz2 +  
   x*NCi2*omzi*zi*lz2 +  
   xi*omzi*zi*lz2 -  
   NCi2*xi*omzi*zi*lz2 -  
   opxi*omzi*zi*lz2 +  
   NCi2*opxi*omzi*zi*lz2 -  
   xi*opxi*omzi*zi*lz2 +  
   NCi2*xi*opxi*omzi*zi*lz2 -  
   4*NC*x*z2*lz2 + 4*x*NCi*z2*lz2 +  
   2*NC*xi*z2*lz2 -  
   2*NCi*xi*z2*lz2 -  
   2*NC*opxi*z2*lz2 +  
   2*NCi*opxi*z2*lz2 -  
   2*NC*xi*opxi*z2*lz2 +  
   2*NCi*xi*opxi*z2*lz2 -  
   15*x*sqrtxz1i*opzi*lz2 +  
   5*x*z*sqrtxz1i*opzi*lz2 +  
   15*x*NCi2*sqrtxz1i*opzi*lz2 -  
   5*x*z*NCi2*sqrtxz1i*opzi*lz2 +  
   20*sqrtxz1i*x2*opzi*lz2 -  
   20*NCi2*sqrtxz1i*x2*opzi*lz2 -  
   2*x*zi*opzi*lz2 -  
   5*sqrtxz1*x*zi*opzi*lz2 +  
   2*x*NCi2*zi*opzi*lz2 +  
   5*sqrtxz1*x*NCi2*zi*opzi*lz2 +  
   2*x*sqrtxz1i*lspec7_2 -  
   x*z*sqrtxz1i*lspec7_2 -  
   2*x*NCi2*sqrtxz1i*lspec7_2 +  
   x*z*NCi2*sqrtxz1i*lspec7_2 -  
   4*sqrtxz1i*x2*lspec7_2 +  
   4*NCi2*sqrtxz1i*x2* 
    lspec7_2 +  
   sqrtxz1*x*zi*lspec7_2 -  
   sqrtxz1*x*NCi2*zi*lspec7_2 -  
   3*x*sqrtxz1i*opzi*lspec7_2 +  
   x*z*sqrtxz1i*opzi*lspec7_2 +  
   3*x*NCi2*sqrtxz1i*opzi* 
    lspec7_2 -  
   x*z*NCi2*sqrtxz1i*opzi* 
    lspec7_2 +  
   4*sqrtxz1i*x2*opzi* 
    lspec7_2 -  
   4*NCi2*sqrtxz1i*x2*opzi* 
    lspec7_2 -  
   sqrtxz1*x*zi*opzi*lspec7_2 +  
   sqrtxz1*x*NCi2*zi*opzi* 
    lspec7_2;
};
double C2LQ2QP1_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double z2 = z * z;
  const double lx  = log(x);
  const double lz  = log(z);
  const double lz2 = lz * lz;
  const double lomx  = log(1 - x);
  const double lomz  = log(1 - z);
  const double Li2z  = apfel::dilog(z);
  return -2*NC + (7*NC*x)/2. + NC*z + NC*x*Li2z + 2*NC*x*z*Li2z -  
   2*NC*x*lomx + NC*x*z*lomx + 4*NC*x*lx - 2*NC*x*z*lx -  
   2*NC*x*lomz + NC*x*z*lomz - NC*lz - NC*x*lz -  
   2*NC*z*lz + 4*NC*x*z*lz - NC*x*lomx*lz -  
   2*NC*x*z*lomx*lz + 2*NC*x*lx*lz +  
   4*NC*x*z*lx*lz + 2*NCi - (7*x*NCi)/2. -  
   z*NCi - x*Li2z*NCi - 2*x*z*Li2z*NCi +  
   2*x*lomx*NCi - x*z*lomx*NCi -  
   4*x*lx*NCi + 2*x*z*lx*NCi +  
   2*x*lomz*NCi - x*z*lomz*NCi +  
   lz*NCi + x*lz*NCi + 2*z*lz*NCi -  
   4*x*z*lz*NCi + x*lomx*lz*NCi +  
   2*x*z*lomx*lz*NCi - 2*x*lx*lz*NCi -  
   4*x*z*lx*lz*NCi - (NC*x*pi2)/6. -  
   (NC*x*z*pi2)/3. + (x*NCi*pi2)/6. +  
   (x*z*NCi*pi2)/3. + NC*z2 - (7*NC*x*z2)/2. +  
   NC*x*lomx*z2 - 2*NC*x*lx*z2 +  
   NC*x*lomz*z2 + NC*x*lz*z2 - NCi*z2 +  
   (7*x*NCi*z2)/2. - x*lomx*NCi*z2 +  
   2*x*lx*NCi*z2 - x*lomz*NCi*z2 -  
   x*lz*NCi*z2 - NC*x*lz2 -  
   2*NC*x*z*lz2 + x*NCi*lz2 +  
   2*x*z*NCi*lz2;
};
double C2LQ2QP2_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double poly2    = 1 + 2*x + x*x - 4*x*z;
  const double poly2i   = 1. / poly2;
  const double sqrtxz2  = sqrt(poly2);
  const double sqrtxz2i = 1. / sqrtxz2;
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double x2 = x * x;
  const double x3 = x * x2;
  const double x4 = x * x3;
  const double x5 = x * x4;
  const double xi  = 1. / x;
  const double z2 = z * z;
  const double lx  = log(x);
  const double lx2 = lx * lx;
  const double lz  = log(z);
  const double lomx  = log(1 - x);
  const double lomz  = log(1 - z);
  const double Li2x  = apfel::dilog(x);
  const double lspec5   = log(1 - sqrtxz2 + x);
  const double lspec6   = log(1 + sqrtxz2 + x);
  const double li2spec5  = apfel::dilog(0.5 - sqrtxz2/2. - x/2.);
  const double li2spec6  = apfel::dilog(0.5 + sqrtxz2/2. - x/2.);
  const double li2spec7  = apfel::dilog(0.5 - xi/2. - (sqrtxz2*xi)/2.);
  const double li2spec8  = apfel::dilog(0.5 - xi/2. + (sqrtxz2*xi)/2.);
  return (-25*NC)/3. + (7*NC*x)/3. + 10*NC*z - 10*NC*x*z + 2*NC*x*Li2x -  
   2*NC*lomx - (3*NC*lx)/4. - (33*NC*x*lx)/4. + 5*NC*z*lx +  
   5*NC*x*z*lx - 2*NC*lomz - 2*NC*x*lx*lomz -  
   (13*NC*lz)/4. + (3*NC*x*lz)/4. + 5*NC*z*lz - 5*NC*x*z*lz -  
   NC*x*lx*lz + (25*NCi)/3. - (7*x*NCi)/3. -  
   10*z*NCi + 10*x*z*NCi - 2*x*Li2x*NCi +  
   2*lomx*NCi + (3*lx*NCi)/4. +  
   (33*x*lx*NCi)/4. - 5*z*lx*NCi -  
   5*x*z*lx*NCi + 2*lomz*NCi +  
   2*x*lx*lomz*NCi + (13*lz*NCi)/4. -  
   (3*x*lz*NCi)/4. - 5*z*lz*NCi +  
   5*x*z*lz*NCi + x*lx*lz*NCi -  
   (NC*x*pi2)/3. + (x*NCi*pi2)/3. +  
   (NC*lx*poly2i)/4. + (NC*x*lx*poly2i)/2. -  
   (NC*lz*poly2i)/4. + (NC*x*lz*poly2i)/2. -  
   (lx*NCi*poly2i)/4. -  
   (x*lx*NCi*poly2i)/2. +  
   (lz*NCi*poly2i)/4. -  
   (x*lz*NCi*poly2i)/2. +  
   (NC*li2spec5*sqrtxz2i)/4. +  
   2*NC*x*li2spec5*sqrtxz2i -  
   (NC*z*li2spec5*sqrtxz2i)/2. -  
   10*NC*x*z*li2spec5*sqrtxz2i -  
   (NC*li2spec6*sqrtxz2i)/4. -  
   2*NC*x*li2spec6*sqrtxz2i +  
   (NC*z*li2spec6*sqrtxz2i)/2. +  
   10*NC*x*z*li2spec6*sqrtxz2i -  
   (NC*li2spec7*sqrtxz2i)/4. -  
   2*NC*x*li2spec7*sqrtxz2i +  
   (NC*z*li2spec7*sqrtxz2i)/ 
    2. + 10*NC*x*z*li2spec7* 
    sqrtxz2i + (NC*li2spec8* 
      sqrtxz2i)/4. + 2*NC*x* 
    li2spec8*sqrtxz2i -  
   (NC*z*li2spec8*sqrtxz2i)/ 
    2. - 10*NC*x*z*li2spec8* 
    sqrtxz2i - (NC*lx*lspec5*sqrtxz2i)/4. -  
   2*NC*x*lx*lspec5*sqrtxz2i +  
   (NC*z*lx*lspec5*sqrtxz2i)/2. +  
   10*NC*x*z*lx*lspec5*sqrtxz2i +  
   (NC*lx*lspec6*sqrtxz2i)/4. +  
   2*NC*x*lx*lspec6*sqrtxz2i -  
   (NC*z*lx*lspec6*sqrtxz2i)/2. -  
   10*NC*x*z*lx*lspec6*sqrtxz2i -  
   (li2spec5*NCi*sqrtxz2i)/4. -  
   2*x*li2spec5*NCi*sqrtxz2i +  
   (z*li2spec5*NCi*sqrtxz2i)/2. +  
   10*x*z*li2spec5*NCi*sqrtxz2i +  
   (li2spec6*NCi*sqrtxz2i)/4. +  
   2*x*li2spec6*NCi*sqrtxz2i -  
   (z*li2spec6*NCi*sqrtxz2i)/2. -  
   10*x*z*li2spec6*NCi*sqrtxz2i +  
   (li2spec7*NCi* 
      sqrtxz2i)/4. + 2*x* 
    li2spec7*NCi* 
    sqrtxz2i - (z*li2spec7* 
      NCi*sqrtxz2i)/2. -  
   10*x*z*li2spec7*NCi* 
    sqrtxz2i - (li2spec8* 
      NCi*sqrtxz2i)/4. -  
   2*x*li2spec8*NCi* 
    sqrtxz2i + (z*li2spec8* 
      NCi*sqrtxz2i)/2. +  
   10*x*z*li2spec8*NCi* 
    sqrtxz2i + (lx*lspec5*NCi* 
      sqrtxz2i)/4. + 2*x*lx*lspec5*NCi* 
    sqrtxz2i - (z*lx*lspec5*NCi* 
      sqrtxz2i)/2. - 10*x*z*lx*lspec5*NCi* 
    sqrtxz2i - (lx*lspec6*NCi* 
      sqrtxz2i)/4. - 2*x*lx*lspec6*NCi* 
    sqrtxz2i + (z*lx*lspec6*NCi* 
      sqrtxz2i)/2. + 10*x*z*lx*lspec6*NCi* 
    sqrtxz2i - (3*NC*x*li2spec5*poly2i* 
      sqrtxz2i)/8. + (3*NC*x*li2spec6* 
      poly2i*sqrtxz2i)/8. +  
   (3*NC*x*li2spec7*poly2i* 
      sqrtxz2i)/8. - (3*NC*x* 
      li2spec8*poly2i* 
      sqrtxz2i)/8. + (3*NC*x*lx*lspec5*poly2i* 
      sqrtxz2i)/8. - (3*NC*x*lx*lspec6*poly2i* 
      sqrtxz2i)/8. + (3*x*li2spec5*NCi* 
      poly2i*sqrtxz2i)/8. -  
   (3*x*li2spec6*NCi*poly2i* 
      sqrtxz2i)/8. - (3*x* 
      li2spec7*NCi* 
      poly2i*sqrtxz2i)/8. +  
   (3*x*li2spec8*NCi* 
      poly2i*sqrtxz2i)/8. -  
   (3*x*lx*lspec5*NCi*poly2i*sqrtxz2i)/ 
    8. + (3*x*lx*lspec6*NCi*poly2i* 
      sqrtxz2i)/8. + (10*NC*xi)/9. +  
   (2*NC*lomx*xi)/3. + (NC*lx*xi)/4. +  
   (2*NC*lomz*xi)/3. + (11*NC*lz*xi)/12. -  
   (10*NCi*xi)/9. - (2*lomx*NCi*xi)/3. -  
   (lx*NCi*xi)/4. -  
   (2*lomz*NCi*xi)/3. -  
   (11*lz*NCi*xi)/12. -  
   (NC*lx*poly2i*xi)/4. -  
   (NC*lz*poly2i*xi)/4. +  
   (lx*NCi*poly2i*xi)/4. +  
   (lz*NCi*poly2i*xi)/4. -  
   (NC*li2spec5*sqrtxz2i*xi)/8. +  
   (NC*li2spec6*sqrtxz2i*xi)/8. +  
   (NC*li2spec7*sqrtxz2i* 
      xi)/8. - (NC*li2spec8* 
      sqrtxz2i*xi)/8. +  
   (NC*lx*lspec5*sqrtxz2i*xi)/8. -  
   (NC*lx*lspec6*sqrtxz2i*xi)/8. +  
   (li2spec5*NCi*sqrtxz2i*xi)/8. -  
   (li2spec6*NCi*sqrtxz2i*xi)/8. -  
   (li2spec7*NCi* 
      sqrtxz2i*xi)/8. +  
   (li2spec8*NCi* 
      sqrtxz2i*xi)/8. -  
   (lx*lspec5*NCi*sqrtxz2i*xi)/8. +  
   (lx*lspec6*NCi*sqrtxz2i*xi)/8. +  
   (NC*li2spec5*poly2i*sqrtxz2i*xi)/ 
    8. - (NC*li2spec6*poly2i*sqrtxz2i* 
      xi)/8. - (NC*li2spec7* 
      poly2i*sqrtxz2i*xi)/8. +  
   (NC*li2spec8*poly2i* 
      sqrtxz2i*xi)/8. -  
   (NC*lx*lspec5*poly2i*sqrtxz2i*xi)/ 
    8. + (NC*lx*lspec6*poly2i*sqrtxz2i* 
      xi)/8. - (li2spec5*NCi*poly2i* 
      sqrtxz2i*xi)/8. +  
   (li2spec6*NCi*poly2i*sqrtxz2i* 
      xi)/8. + (li2spec7* 
      NCi*poly2i*sqrtxz2i*xi)/8. -  
   (li2spec8*NCi*poly2i* 
      sqrtxz2i*xi)/8. +  
   (lx*lspec5*NCi*poly2i*sqrtxz2i* 
      xi)/8. - (lx*lspec6*NCi*poly2i* 
      sqrtxz2i*xi)/8. + (44*NC*x2)/9. +  
   (4*NC*lomx*x2)/3. - (17*NC*lx*x2)/4. +  
   (4*NC*lomz*x2)/3. + (19*NC*lz*x2)/12. -  
   (44*NCi*x2)/9. - (4*lomx*NCi*x2)/3. +  
   (17*lx*NCi*x2)/4. -  
   (4*lomz*NCi*x2)/3. -  
   (19*lz*NCi*x2)/12. -  
   (NC*lx*poly2i*x2)/2. +  
   (NC*lz*poly2i*x2)/2. +  
   (lx*NCi*poly2i*x2)/2. -  
   (lz*NCi*poly2i*x2)/2. +  
   (7*NC*li2spec5*sqrtxz2i*x2)/4. -  
   (7*NC*z*li2spec5*sqrtxz2i*x2)/2. -  
   (7*NC*li2spec6*sqrtxz2i*x2)/4. +  
   (7*NC*z*li2spec6*sqrtxz2i*x2)/2. -  
   (7*NC*li2spec7*sqrtxz2i* 
      x2)/4. + (7*NC*z*li2spec7* 
      sqrtxz2i*x2)/2. +  
   (7*NC*li2spec8*sqrtxz2i* 
      x2)/4. - (7*NC*z*li2spec8* 
      sqrtxz2i*x2)/2. -  
   (7*NC*lx*lspec5*sqrtxz2i*x2)/4. +  
   (7*NC*z*lx*lspec5*sqrtxz2i*x2)/2. +  
   (7*NC*lx*lspec6*sqrtxz2i*x2)/4. -  
   (7*NC*z*lx*lspec6*sqrtxz2i*x2)/2. -  
   (7*li2spec5*NCi*sqrtxz2i*x2)/4. +  
   (7*z*li2spec5*NCi*sqrtxz2i*x2)/ 
    2. + (7*li2spec6*NCi*sqrtxz2i*x2)/ 
    4. - (7*z*li2spec6*NCi*sqrtxz2i* 
      x2)/2. + (7*li2spec7* 
      NCi*sqrtxz2i*x2)/4. -  
   (7*z*li2spec7*NCi* 
      sqrtxz2i*x2)/2. -  
   (7*li2spec8*NCi* 
      sqrtxz2i*x2)/4. +  
   (7*z*li2spec8*NCi* 
      sqrtxz2i*x2)/2. +  
   (7*lx*lspec5*NCi*sqrtxz2i*x2)/4. -  
   (7*z*lx*lspec5*NCi*sqrtxz2i*x2)/2. -  
   (7*lx*lspec6*NCi*sqrtxz2i*x2)/4. +  
   (7*z*lx*lspec6*NCi*sqrtxz2i*x2)/2. -  
   (NC*lx*poly2i*x3)/4. -  
   (NC*lz*poly2i*x3)/4. +  
   (lx*NCi*poly2i*x3)/4. +  
   (lz*NCi*poly2i*x3)/4. +  
   (NC*li2spec5*sqrtxz2i*x3)/8. -  
   (NC*li2spec6*sqrtxz2i*x3)/8. -  
   (NC*li2spec7*sqrtxz2i* 
      x3)/8. + (NC*li2spec8* 
      sqrtxz2i*x3)/8. -  
   (NC*lx*lspec5*sqrtxz2i*x3)/8. +  
   (NC*lx*lspec6*sqrtxz2i*x3)/8. -  
   (li2spec5*NCi*sqrtxz2i*x3)/8. +  
   (li2spec6*NCi*sqrtxz2i*x3)/8. +  
   (li2spec7*NCi* 
      sqrtxz2i*x3)/8. -  
   (li2spec8*NCi* 
      sqrtxz2i*x3)/8. +  
   (lx*lspec5*NCi*sqrtxz2i*x3)/8. -  
   (lx*lspec6*NCi*sqrtxz2i*x3)/8. +  
   (3*NC*li2spec5*poly2i*sqrtxz2i*x3)/ 
    8. - (3*NC*li2spec6*poly2i*sqrtxz2i* 
      x3)/8. - (3*NC*li2spec7* 
      poly2i*sqrtxz2i*x3)/8. +  
   (3*NC*li2spec8*poly2i* 
      sqrtxz2i*x3)/8. -  
   (3*NC*lx*lspec5*poly2i*sqrtxz2i*x3)/ 
    8. + (3*NC*lx*lspec6*poly2i*sqrtxz2i* 
      x3)/8. - (3*li2spec5*NCi*poly2i* 
      sqrtxz2i*x3)/8. +  
   (3*li2spec6*NCi*poly2i*sqrtxz2i* 
      x3)/8. + (3*li2spec7* 
      NCi*poly2i*sqrtxz2i*x3)/8. -  
   (3*li2spec8*NCi* 
      poly2i*sqrtxz2i*x3)/8. +  
   (3*lx*lspec5*NCi*poly2i*sqrtxz2i* 
      x3)/8. - (3*lx*lspec6*NCi*poly2i* 
      sqrtxz2i*x3)/8. + (NC*lx*poly2i*x4)/4. -  
   (NC*lz*poly2i*x4)/4. -  
   (lx*NCi*poly2i*x4)/4. +  
   (lz*NCi*poly2i*x4)/4. -  
   (NC*li2spec5*poly2i*sqrtxz2i*x5)/ 
    8. + (NC*li2spec6*poly2i*sqrtxz2i* 
      x5)/8. + (NC*li2spec7* 
      poly2i*sqrtxz2i*x5)/8. -  
   (NC*li2spec8*poly2i* 
      sqrtxz2i*x5)/8. +  
   (NC*lx*lspec5*poly2i*sqrtxz2i*x5)/ 
    8. - (NC*lx*lspec6*poly2i*sqrtxz2i* 
      x5)/8. + (li2spec5*NCi*poly2i* 
      sqrtxz2i*x5)/8. -  
   (li2spec6*NCi*poly2i*sqrtxz2i* 
      x5)/8. - (li2spec7* 
      NCi*poly2i*sqrtxz2i*x5)/8. +  
   (li2spec8*NCi*poly2i* 
      sqrtxz2i*x5)/8. -  
   (lx*lspec5*NCi*poly2i*sqrtxz2i* 
      x5)/8. + (lx*lspec6*NCi*poly2i* 
      sqrtxz2i*x5)/8. +  
   10*NC*x*li2spec5*sqrtxz2i*z2 -  
   10*NC*x*li2spec6*sqrtxz2i*z2 -  
   10*NC*x*li2spec7*sqrtxz2i* 
    z2 + 10*NC*x*li2spec8* 
    sqrtxz2i*z2 - 10*NC*x*lx*lspec5* 
    sqrtxz2i*z2 + 10*NC*x*lx*lspec6* 
    sqrtxz2i*z2 - 10*x*li2spec5*NCi* 
    sqrtxz2i*z2 + 10*x*li2spec6*NCi* 
    sqrtxz2i*z2 + 10*x* 
    li2spec7*NCi* 
    sqrtxz2i*z2 - 10*x* 
    li2spec8*NCi* 
    sqrtxz2i*z2 + 10*x*lx*lspec5*NCi* 
    sqrtxz2i*z2 - 10*x*lx*lspec6*NCi* 
    sqrtxz2i*z2 + 2*NC*x*lx2 -  
   2*x*NCi*lx2;
};
double C2LQ2QP3_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double sqrtxz1  = sqrt(1 - 2*z + z*z + 4*x*z);
  const double poly2    = 1 + 2*x + x*x - 4*x*z;
  const double poly2i   = 1. / poly2;
  const double sqrtxz2  = sqrt(poly2);
  const double sqrtxz2i = 1. / sqrtxz2;
  const double sqrtxz3  = sqrt(x/z);
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double l2  = log(2);
  const double l22 = l2 * l2;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double x2 = x * x;
  const double x3 = x * x2;
  const double x4 = x * x3;
  const double x5 = x * x4;
  const double xi  = 1. / x;
  const double z2 = z * z;
  const double zi  = 1. / z;
  const double opxi = 1. / ( 1 + x );
  const double lx  = log(x);
  const double lx2 = lx * lx;
  const double lz  = log(z);
  const double lz2 = lz * lz;
  const double lomx  = log(1 - x);
  const double lomz  = log(1 - z);
  const double lopx  = log(1 + x);
  const double Li2x  = apfel::dilog(x);
  const double Li2mx = apfel::dilog(-x);
  const double Li2z  = apfel::dilog(z);
  const double xmzi  = 1. / ( x - z );
  const double lxpz    = log(x + z);
  const double lopxz   = log(1 + x*z);
  const double lopxzi  = log(1 + x*zi);
  const double li2omxzi = apfel::dilog(1 - x*zi);
  const double lspec1   = log(1 + sqrtxz1 - z);
  const double lspec1_2 = lspec1 * lspec1;
  const double lspec2   = log(1 + sqrtxz1 + z);
  const double lspec3   = log(sqrtxz3);
  const double lspec4   = log(sqrtxz3*z);
  const double lspec5   = log(1 - sqrtxz2 + x);
  const double lspec6   = log(1 + sqrtxz2 + x);
  const double li2spec1  = apfel::dilog(0.5 - sqrtxz1/2. - z/2.);
  const double li2spec2  = apfel::dilog(0.5 - sqrtxz1/2. + z/2.);
  const double li2spec3  = apfel::dilog(0.5 - zi/2. - (sqrtxz1*zi)/2.);
  const double li2spec4  = apfel::dilog(0.5 + zi/2. - (sqrtxz1*zi)/2.);
  const double li2spec5  = apfel::dilog(0.5 - sqrtxz2/2. - x/2.);
  const double li2spec6  = apfel::dilog(0.5 + sqrtxz2/2. - x/2.);
  const double li2spec7  = apfel::dilog(0.5 - xi/2. - (sqrtxz2*xi)/2.);
  const double li2spec8  = apfel::dilog(0.5 - xi/2. + (sqrtxz2*xi)/2.);
  const double li2spec19 = apfel::dilog(-(x*z));
  const double li2spec20 = apfel::dilog(-(x*zi));
  const double atanspec1 = atan(sqrtxz3);
  const double atanspec2 = atan(sqrtxz3*z);
  const double itani1 = InvTanInt(-sqrtxz3);
  const double itani2 = InvTanInt(sqrtxz3);
  const double itani3 = InvTanInt(sqrtxz3*z);
  return NC - NC*x - NC*z + NC*x*z - (NC*sqrtxz3*itani1)/2. -  
   (3*NC*sqrtxz3*x*z*itani1)/2. +  
   (NC*sqrtxz3*itani2)/2. +  
   (3*NC*sqrtxz3*x*z*itani2)/2. -  
   NC*sqrtxz3*itani3 - 3*NC*sqrtxz3*x*z*itani3 -  
   4*NC*x*Li2mx + 8*NC*x*z*Li2mx + 4*NC*x*z*Li2x - 2*NC*x*Li2z +  
   4*NC*x*z*Li2z + 4*NC*x*z*li2spec20 -  
   4*NC*x*z*li2spec3 +  
   4*NC*x*z*li2spec4 +  
   2*NC*x*li2omxzi - 4*NC*x*z*li2omxzi -  
   4*NC*sqrtxz1*x*l2 - NC*sqrtxz3*atanspec1*lspec3 -  
   3*NC*sqrtxz3*x*z*atanspec1*lspec3 + (NC*lx)/4. +  
   (11*NC*x*lx)/4. - 2*NC*sqrtxz1*x*lx - (NC*z*lx)/2. -  
   (5*NC*x*z*lx)/2. + 4*NC*x*z*l2*lx +  
   4*NC*x*z*lomx*lx - 4*NC*x*lx*lopx +  
   8*NC*x*z*lx*lopx + 4*NC*sqrtxz1*x*lspec1 -  
   4*NC*x*z*l2*lspec1 -  
   4*NC*x*z*lx*lspec1 + (3*NC*lz)/4. -  
   (NC*x*lz)/4. - 2*NC*sqrtxz1*x*lz - (NC*z*lz)/2. +  
   (5*NC*x*z*lz)/2. - 4*NC*x*z*l2*lz + 2*NC*x*z*lx*lz -  
   2*NC*x*lomz*lz + 4*NC*x*z*lomz*lz +  
   NC*sqrtxz3*atanspec2*lspec4 +  
   3*NC*sqrtxz3*x*z*atanspec2*lspec4 +  
   4*NC*x*z*l2*lspec2 -  
   4*NC*x*z*lspec1*lspec2 +  
   4*NC*x*z*lz*lspec2 + 2*NC*x*z*lx*lxpz -  
   2*NC*x*z*lz*lxpz + 2*NC*x*z*lx*lopxzi -  
   2*NC*x*z*lz*lopxzi - NCi + x*NCi +  
   z*NCi - x*z*NCi +  
   (sqrtxz3*itani1*NCi)/2. +  
   (3*sqrtxz3*x*z*itani1*NCi)/2. -  
   (sqrtxz3*itani2*NCi)/2. -  
   (3*sqrtxz3*x*z*itani2*NCi)/2. +  
   sqrtxz3*itani3*NCi +  
   3*sqrtxz3*x*z*itani3*NCi + 4*x*Li2mx*NCi -  
   8*x*z*Li2mx*NCi - 4*x*z*Li2x*NCi +  
   2*x*Li2z*NCi - 4*x*z*Li2z*NCi -  
   4*x*z*li2spec20*NCi +  
   4*x*z*li2spec3*NCi -  
   4*x*z*li2spec4*NCi -  
   2*x*li2omxzi*NCi +  
   4*x*z*li2omxzi*NCi + 4*sqrtxz1*x*l2*NCi +  
   sqrtxz3*atanspec1*lspec3*NCi +  
   3*sqrtxz3*x*z*atanspec1*lspec3*NCi -  
   (lx*NCi)/4. - (11*x*lx*NCi)/4. +  
   2*sqrtxz1*x*lx*NCi + (z*lx*NCi)/2. +  
   (5*x*z*lx*NCi)/2. - 4*x*z*l2*lx*NCi -  
   4*x*z*lomx*lx*NCi + 4*x*lx*lopx*NCi -  
   8*x*z*lx*lopx*NCi -  
   4*sqrtxz1*x*lspec1*NCi +  
   4*x*z*l2*lspec1*NCi +  
   4*x*z*lx*lspec1*NCi - (3*lz*NCi)/4. +  
   (x*lz*NCi)/4. + 2*sqrtxz1*x*lz*NCi +  
   (z*lz*NCi)/2. - (5*x*z*lz*NCi)/2. +  
   4*x*z*l2*lz*NCi - 2*x*z*lx*lz*NCi +  
   2*x*lomz*lz*NCi - 4*x*z*lomz*lz*NCi -  
   sqrtxz3*atanspec2*lspec4*NCi -  
   3*sqrtxz3*x*z*atanspec2*lspec4*NCi -  
   4*x*z*l2*lspec2*NCi +  
   4*x*z*lspec1*lspec2*NCi -  
   4*x*z*lz*lspec2*NCi -  
   2*x*z*lx*lxpz*NCi + 2*x*z*lz*lxpz*NCi -  
   2*x*z*lx*lopxzi*NCi +  
   2*x*z*lz*lopxzi*NCi - (2*NC*x*z*pi2)/3. +  
   (2*x*z*NCi*pi2)/3. - (NC*lx*poly2i)/4. -  
   (NC*x*lx*poly2i)/2. + (NC*lz*poly2i)/4. -  
   (NC*x*lz*poly2i)/2. + (lx*NCi*poly2i)/4. +  
   (x*lx*NCi*poly2i)/2. -  
   (lz*NCi*poly2i)/4. +  
   (x*lz*NCi*poly2i)/2. -  
   (NC*li2spec5*sqrtxz2i)/4. +  
   2*NC*x*li2spec5*sqrtxz2i +  
   (NC*z*li2spec5*sqrtxz2i)/2. -  
   6*NC*x*z*li2spec5*sqrtxz2i +  
   (NC*li2spec6*sqrtxz2i)/4. -  
   2*NC*x*li2spec6*sqrtxz2i -  
   (NC*z*li2spec6*sqrtxz2i)/2. +  
   6*NC*x*z*li2spec6*sqrtxz2i +  
   (NC*li2spec7*sqrtxz2i)/4. -  
   2*NC*x*li2spec7*sqrtxz2i -  
   (NC*z*li2spec7*sqrtxz2i)/ 
    2. + 6*NC*x*z*li2spec7* 
    sqrtxz2i - (NC*li2spec8* 
      sqrtxz2i)/4. + 2*NC*x* 
    li2spec8*sqrtxz2i +  
   (NC*z*li2spec8*sqrtxz2i)/ 
    2. - 6*NC*x*z*li2spec8* 
    sqrtxz2i + (NC*lx*lspec5*sqrtxz2i)/4. -  
   2*NC*x*lx*lspec5*sqrtxz2i -  
   (NC*z*lx*lspec5*sqrtxz2i)/2. +  
   6*NC*x*z*lx*lspec5*sqrtxz2i -  
   (NC*lx*lspec6*sqrtxz2i)/4. +  
   2*NC*x*lx*lspec6*sqrtxz2i +  
   (NC*z*lx*lspec6*sqrtxz2i)/2. -  
   6*NC*x*z*lx*lspec6*sqrtxz2i +  
   (li2spec5*NCi*sqrtxz2i)/4. -  
   2*x*li2spec5*NCi*sqrtxz2i -  
   (z*li2spec5*NCi*sqrtxz2i)/2. +  
   6*x*z*li2spec5*NCi*sqrtxz2i -  
   (li2spec6*NCi*sqrtxz2i)/4. +  
   2*x*li2spec6*NCi*sqrtxz2i +  
   (z*li2spec6*NCi*sqrtxz2i)/2. -  
   6*x*z*li2spec6*NCi*sqrtxz2i -  
   (li2spec7*NCi* 
      sqrtxz2i)/4. + 2*x* 
    li2spec7*NCi* 
    sqrtxz2i + (z*li2spec7* 
      NCi*sqrtxz2i)/2. -  
   6*x*z*li2spec7*NCi* 
    sqrtxz2i + (li2spec8* 
      NCi*sqrtxz2i)/4. -  
   2*x*li2spec8*NCi* 
    sqrtxz2i - (z*li2spec8* 
      NCi*sqrtxz2i)/2. +  
   6*x*z*li2spec8*NCi* 
    sqrtxz2i - (lx*lspec5*NCi* 
      sqrtxz2i)/4. + 2*x*lx*lspec5*NCi* 
    sqrtxz2i + (z*lx*lspec5*NCi* 
      sqrtxz2i)/2. - 6*x*z*lx*lspec5*NCi* 
    sqrtxz2i + (lx*lspec6*NCi* 
      sqrtxz2i)/4. - 2*x*lx*lspec6*NCi* 
    sqrtxz2i - (z*lx*lspec6*NCi* 
      sqrtxz2i)/2. + 6*x*z*lx*lspec6*NCi* 
    sqrtxz2i + (3*NC*x*li2spec5*poly2i* 
      sqrtxz2i)/8. - (3*NC*x*li2spec6* 
      poly2i*sqrtxz2i)/8. -  
   (3*NC*x*li2spec7*poly2i* 
      sqrtxz2i)/8. + (3*NC*x* 
      li2spec8*poly2i* 
      sqrtxz2i)/8. - (3*NC*x*lx*lspec5*poly2i* 
      sqrtxz2i)/8. + (3*NC*x*lx*lspec6*poly2i* 
      sqrtxz2i)/8. - (3*x*li2spec5*NCi* 
      poly2i*sqrtxz2i)/8. +  
   (3*x*li2spec6*NCi*poly2i* 
      sqrtxz2i)/8. + (3*x* 
      li2spec7*NCi* 
      poly2i*sqrtxz2i)/8. -  
   (3*x*li2spec8*NCi* 
      poly2i*sqrtxz2i)/8. +  
   (3*x*lx*lspec5*NCi*poly2i*sqrtxz2i)/ 
    8. - (3*x*lx*lspec6*NCi*poly2i* 
      sqrtxz2i)/8. + (NC*sqrtxz3*z*itani1*xi)/2. -  
   (NC*sqrtxz3*z*itani2*xi)/2. +  
   NC*sqrtxz3*z*itani3*xi + 4*NC*z*Li2mx*xi -  
   4*NC*z*Li2x*xi - 2*NC*z*li2spec19*xi -  
   2*NC*z*li2spec20*xi +  
   NC*sqrtxz3*z*atanspec1*lspec3*xi -  
   (NC*lx*xi)/4. + 8*NC*z*lx*xi -  
   4*NC*z*lomx*lx*xi + 4*NC*z*lx*lopx*xi -  
   (NC*lz*xi)/4. + 2*NC*z*lx*lz*xi -  
   NC*sqrtxz3*z*atanspec2*lspec4*xi -  
   2*NC*z*lx*lopxz*xi -  
   2*NC*z*lz*lopxz*xi -  
   2*NC*z*lx*lopxzi*xi +  
   2*NC*z*lz*lopxzi*xi -  
   (sqrtxz3*z*itani1*NCi*xi)/2. +  
   (sqrtxz3*z*itani2*NCi*xi)/2. -  
   sqrtxz3*z*itani3*NCi*xi -  
   4*z*Li2mx*NCi*xi + 4*z*Li2x*NCi*xi +  
   2*z*li2spec19*NCi*xi +  
   2*z*li2spec20*NCi*xi -  
   sqrtxz3*z*atanspec1*lspec3*NCi*xi +  
   (lx*NCi*xi)/4. - 8*z*lx*NCi*xi +  
   4*z*lomx*lx*NCi*xi -  
   4*z*lx*lopx*NCi*xi +  
   (lz*NCi*xi)/4. -  
   2*z*lx*lz*NCi*xi +  
   sqrtxz3*z*atanspec2*lspec4*NCi*xi +  
   2*z*lx*lopxz*NCi*xi +  
   2*z*lz*lopxz*NCi*xi +  
   2*z*lx*lopxzi*NCi*xi -  
   2*z*lz*lopxzi*NCi*xi +  
   (2*NC*z*pi2*xi)/3. -  
   (2*z*NCi*pi2*xi)/3. +  
   (NC*lx*poly2i*xi)/4. +  
   (NC*lz*poly2i*xi)/4. -  
   (lx*NCi*poly2i*xi)/4. -  
   (lz*NCi*poly2i*xi)/4. +  
   (NC*li2spec5*sqrtxz2i*xi)/8. -  
   (NC*li2spec6*sqrtxz2i*xi)/8. -  
   (NC*li2spec7*sqrtxz2i* 
      xi)/8. + (NC*li2spec8* 
      sqrtxz2i*xi)/8. -  
   (NC*lx*lspec5*sqrtxz2i*xi)/8. +  
   (NC*lx*lspec6*sqrtxz2i*xi)/8. -  
   (li2spec5*NCi*sqrtxz2i*xi)/8. +  
   (li2spec6*NCi*sqrtxz2i*xi)/8. +  
   (li2spec7*NCi* 
      sqrtxz2i*xi)/8. -  
   (li2spec8*NCi* 
      sqrtxz2i*xi)/8. +  
   (lx*lspec5*NCi*sqrtxz2i*xi)/8. -  
   (lx*lspec6*NCi*sqrtxz2i*xi)/8. -  
   (NC*li2spec5*poly2i*sqrtxz2i*xi)/ 
    8. + (NC*li2spec6*poly2i*sqrtxz2i* 
      xi)/8. + (NC*li2spec7* 
      poly2i*sqrtxz2i*xi)/8. -  
   (NC*li2spec8*poly2i* 
      sqrtxz2i*xi)/8. +  
   (NC*lx*lspec5*poly2i*sqrtxz2i*xi)/ 
    8. - (NC*lx*lspec6*poly2i*sqrtxz2i* 
      xi)/8. + (li2spec5*NCi*poly2i* 
      sqrtxz2i*xi)/8. -  
   (li2spec6*NCi*poly2i*sqrtxz2i* 
      xi)/8. - (li2spec7* 
      NCi*poly2i*sqrtxz2i*xi)/8. +  
   (li2spec8*NCi*poly2i* 
      sqrtxz2i*xi)/8. -  
   (lx*lspec5*NCi*poly2i*sqrtxz2i* 
      xi)/8. + (lx*lspec6*NCi*poly2i* 
      sqrtxz2i*xi)/8. - (7*NC*lx*x2)/4. +  
   (7*NC*lz*x2)/4. + (7*lx*NCi*x2)/4. -  
   (7*lz*NCi*x2)/4. +  
   (NC*lx*poly2i*x2)/2. -  
   (NC*lz*poly2i*x2)/2. -  
   (lx*NCi*poly2i*x2)/2. +  
   (lz*NCi*poly2i*x2)/2. +  
   (9*NC*li2spec5*sqrtxz2i*x2)/4. -  
   (9*NC*z*li2spec5*sqrtxz2i*x2)/2. -  
   (9*NC*li2spec6*sqrtxz2i*x2)/4. +  
   (9*NC*z*li2spec6*sqrtxz2i*x2)/2. -  
   (9*NC*li2spec7*sqrtxz2i* 
      x2)/4. + (9*NC*z*li2spec7* 
      sqrtxz2i*x2)/2. +  
   (9*NC*li2spec8*sqrtxz2i* 
      x2)/4. - (9*NC*z*li2spec8* 
      sqrtxz2i*x2)/2. -  
   (9*NC*lx*lspec5*sqrtxz2i*x2)/4. +  
   (9*NC*z*lx*lspec5*sqrtxz2i*x2)/2. +  
   (9*NC*lx*lspec6*sqrtxz2i*x2)/4. -  
   (9*NC*z*lx*lspec6*sqrtxz2i*x2)/2. -  
   (9*li2spec5*NCi*sqrtxz2i*x2)/4. +  
   (9*z*li2spec5*NCi*sqrtxz2i*x2)/ 
    2. + (9*li2spec6*NCi*sqrtxz2i*x2)/ 
    4. - (9*z*li2spec6*NCi*sqrtxz2i* 
      x2)/2. + (9*li2spec7* 
      NCi*sqrtxz2i*x2)/4. -  
   (9*z*li2spec7*NCi* 
      sqrtxz2i*x2)/2. -  
   (9*li2spec8*NCi* 
      sqrtxz2i*x2)/4. +  
   (9*z*li2spec8*NCi* 
      sqrtxz2i*x2)/2. +  
   (9*lx*lspec5*NCi*sqrtxz2i*x2)/4. -  
   (9*z*lx*lspec5*NCi*sqrtxz2i*x2)/2. -  
   (9*lx*lspec6*NCi*sqrtxz2i*x2)/4. +  
   (9*z*lx*lspec6*NCi*sqrtxz2i*x2)/2. +  
   (NC*lx*poly2i*x3)/4. +  
   (NC*lz*poly2i*x3)/4. -  
   (lx*NCi*poly2i*x3)/4. -  
   (lz*NCi*poly2i*x3)/4. -  
   (NC*li2spec5*sqrtxz2i*x3)/8. +  
   (NC*li2spec6*sqrtxz2i*x3)/8. +  
   (NC*li2spec7*sqrtxz2i* 
      x3)/8. - (NC*li2spec8* 
      sqrtxz2i*x3)/8. +  
   (NC*lx*lspec5*sqrtxz2i*x3)/8. -  
   (NC*lx*lspec6*sqrtxz2i*x3)/8. +  
   (li2spec5*NCi*sqrtxz2i*x3)/8. -  
   (li2spec6*NCi*sqrtxz2i*x3)/8. -  
   (li2spec7*NCi* 
      sqrtxz2i*x3)/8. +  
   (li2spec8*NCi* 
      sqrtxz2i*x3)/8. -  
   (lx*lspec5*NCi*sqrtxz2i*x3)/8. +  
   (lx*lspec6*NCi*sqrtxz2i*x3)/8. -  
   (3*NC*li2spec5*poly2i*sqrtxz2i*x3)/ 
    8. + (3*NC*li2spec6*poly2i*sqrtxz2i* 
      x3)/8. + (3*NC*li2spec7* 
      poly2i*sqrtxz2i*x3)/8. -  
   (3*NC*li2spec8*poly2i* 
      sqrtxz2i*x3)/8. +  
   (3*NC*lx*lspec5*poly2i*sqrtxz2i*x3)/ 
    8. - (3*NC*lx*lspec6*poly2i*sqrtxz2i* 
      x3)/8. + (3*li2spec5*NCi*poly2i* 
      sqrtxz2i*x3)/8. -  
   (3*li2spec6*NCi*poly2i*sqrtxz2i* 
      x3)/8. - (3*li2spec7* 
      NCi*poly2i*sqrtxz2i*x3)/8. +  
   (3*li2spec8*NCi* 
      poly2i*sqrtxz2i*x3)/8. -  
   (3*lx*lspec5*NCi*poly2i*sqrtxz2i* 
      x3)/8. + (3*lx*lspec6*NCi*poly2i* 
      sqrtxz2i*x3)/8. - (NC*lx*poly2i*x4)/4. +  
   (NC*lz*poly2i*x4)/4. +  
   (lx*NCi*poly2i*x4)/4. -  
   (lz*NCi*poly2i*x4)/4. +  
   (NC*li2spec5*poly2i*sqrtxz2i*x5)/ 
    8. - (NC*li2spec6*poly2i*sqrtxz2i* 
      x5)/8. - (NC*li2spec7* 
      poly2i*sqrtxz2i*x5)/8. +  
   (NC*li2spec8*poly2i* 
      sqrtxz2i*x5)/8. -  
   (NC*lx*lspec5*poly2i*sqrtxz2i*x5)/ 
    8. + (NC*lx*lspec6*poly2i*sqrtxz2i* 
      x5)/8. - (li2spec5*NCi*poly2i* 
      sqrtxz2i*x5)/8. +  
   (li2spec6*NCi*poly2i*sqrtxz2i* 
      x5)/8. + (li2spec7* 
      NCi*poly2i*sqrtxz2i*x5)/8. -  
   (li2spec8*NCi*poly2i* 
      sqrtxz2i*x5)/8. +  
   (lx*lspec5*NCi*poly2i*sqrtxz2i* 
      x5)/8. - (lx*lspec6*NCi*poly2i* 
      sqrtxz2i*x5)/8. - 4*NC*z*Li2mx*opxi +  
   4*NC*z*Li2x*opxi + 2*NC*z*li2spec19*opxi +  
   2*NC*z*li2spec20*opxi - 8*NC*z*lx*opxi +  
   4*NC*z*lomx*lx*opxi -  
   4*NC*z*lx*lopx*opxi -  
   2*NC*z*lx*lz*opxi +  
   2*NC*z*lx*lopxz*opxi +  
   2*NC*z*lz*lopxz*opxi +  
   2*NC*z*lx*lopxzi*opxi -  
   2*NC*z*lz*lopxzi*opxi +  
   4*z*Li2mx*NCi*opxi -  
   4*z*Li2x*NCi*opxi -  
   2*z*li2spec19*NCi*opxi -  
   2*z*li2spec20*NCi*opxi +  
   8*z*lx*NCi*opxi -  
   4*z*lomx*lx*NCi*opxi +  
   4*z*lx*lopx*NCi*opxi +  
   2*z*lx*lz*NCi*opxi -  
   2*z*lx*lopxz*NCi*opxi -  
   2*z*lz*lopxz*NCi*opxi -  
   2*z*lx*lopxzi*NCi*opxi +  
   2*z*lz*lopxzi*NCi*opxi -  
   (2*NC*z*pi2*opxi)/3. +  
   (2*z*NCi*pi2*opxi)/3. -  
   4*NC*z*Li2mx*xi*opxi +  
   4*NC*z*Li2x*xi*opxi +  
   2*NC*z*li2spec19*xi*opxi +  
   2*NC*z*li2spec20*xi*opxi -  
   8*NC*z*lx*xi*opxi +  
   4*NC*z*lomx*lx*xi*opxi -  
   4*NC*z*lx*lopx*xi*opxi -  
   2*NC*z*lx*lz*xi*opxi +  
   2*NC*z*lx*lopxz*xi*opxi +  
   2*NC*z*lz*lopxz*xi*opxi +  
   2*NC*z*lx*lopxzi*xi*opxi -  
   2*NC*z*lz*lopxzi*xi*opxi +  
   4*z*Li2mx*NCi*xi*opxi -  
   4*z*Li2x*NCi*xi*opxi -  
   2*z*li2spec19*NCi*xi*opxi -  
   2*z*li2spec20*NCi*xi*opxi +  
   8*z*lx*NCi*xi*opxi -  
   4*z*lomx*lx*NCi*xi*opxi +  
   4*z*lx*lopx*NCi*xi*opxi +  
   2*z*lx*lz*NCi*xi*opxi -  
   2*z*lx*lopxz*NCi*xi*opxi -  
   2*z*lz*lopxz*NCi*xi*opxi -  
   2*z*lx*lopxzi*NCi*xi*opxi +  
   2*z*lz*lopxzi*NCi*xi*opxi -  
   (2*NC*z*pi2*xi*opxi)/3. +  
   (2*z*NCi*pi2*xi*opxi)/3. -  
   2*NC*lx*x2*xmzi + 2*NC*lz*x2*xmzi +  
   2*lx*NCi*x2*xmzi -  
   2*lz*NCi*x2*xmzi +  
   2*NC*lx*x3*xmzi - 2*NC*lz*x3*xmzi -  
   2*lx*NCi*x3*xmzi +  
   2*lz*NCi*x3*xmzi +  
   (23*NC*sqrtxz3*itani1*z2)/2. -  
   (23*NC*sqrtxz3*itani2*z2)/2. +  
   23*NC*sqrtxz3*itani3*z2 - 8*NC*x*Li2x*z2 +  
   8*NC*x*li2spec1*z2 -  
   8*NC*x*li2spec2*z2 -  
   8*NC*x*li2spec19*z2 +  
   23*NC*sqrtxz3*atanspec1*lspec3*z2 -  
   16*NC*x*l2*lx*z2 - 8*NC*x*lomx*lx*z2 +  
   24*NC*x*l2*lspec1*z2 +  
   8*NC*x*lx*lspec1*z2 -  
   16*NC*x*l2*lz*z2 - 4*NC*x*lx*lz*z2 +  
   8*NC*x*lspec1*lz*z2 -  
   23*NC*sqrtxz3*atanspec2*lspec4*z2 +  
   8*NC*x*l2*lspec2*z2 +  
   8*NC*x*lx*lspec2*z2 -  
   8*NC*x*lspec1*lspec2*z2 +  
   8*NC*x*lz*lspec2*z2 +  
   4*NC*x*lx*lxpz*z2 - 4*NC*x*lz*lxpz*z2 -  
   8*NC*x*lx*lopxz*z2 -  
   8*NC*x*lz*lopxz*z2 -  
   4*NC*x*lx*lopxzi*z2 +  
   4*NC*x*lz*lopxzi*z2 -  
   (23*sqrtxz3*itani1*NCi*z2)/2. +  
   (23*sqrtxz3*itani2*NCi*z2)/2. -  
   23*sqrtxz3*itani3*NCi*z2 +  
   8*x*Li2x*NCi*z2 -  
   8*x*li2spec1*NCi*z2 +  
   8*x*li2spec2*NCi*z2 +  
   8*x*li2spec19*NCi*z2 -  
   23*sqrtxz3*atanspec1*lspec3*NCi*z2 +  
   16*x*l2*lx*NCi*z2 +  
   8*x*lomx*lx*NCi*z2 -  
   24*x*l2*lspec1*NCi*z2 -  
   8*x*lx*lspec1*NCi*z2 +  
   16*x*l2*lz*NCi*z2 +  
   4*x*lx*lz*NCi*z2 -  
   8*x*lspec1*lz*NCi*z2 +  
   23*sqrtxz3*atanspec2*lspec4*NCi*z2 -  
   8*x*l2*lspec2*NCi*z2 -  
   8*x*lx*lspec2*NCi*z2 +  
   8*x*lspec1*lspec2*NCi*z2 -  
   8*x*lz*lspec2*NCi*z2 -  
   4*x*lx*lxpz*NCi*z2 +  
   4*x*lz*lxpz*NCi*z2 +  
   8*x*lx*lopxz*NCi*z2 +  
   8*x*lz*lopxz*NCi*z2 +  
   4*x*lx*lopxzi*NCi*z2 -  
   4*x*lz*lopxzi*NCi*z2 +  
   (4*NC*x*pi2*z2)/3. - (4*x*NCi*pi2*z2)/3. +  
   6*NC*x*li2spec5*sqrtxz2i*z2 -  
   6*NC*x*li2spec6*sqrtxz2i*z2 -  
   6*NC*x*li2spec7*sqrtxz2i* 
    z2 + 6*NC*x*li2spec8* 
    sqrtxz2i*z2 - 6*NC*x*lx*lspec5* 
    sqrtxz2i*z2 + 6*NC*x*lx*lspec6* 
    sqrtxz2i*z2 - 6*x*li2spec5*NCi* 
    sqrtxz2i*z2 + 6*x*li2spec6*NCi* 
    sqrtxz2i*z2 + 6*x* 
    li2spec7*NCi* 
    sqrtxz2i*z2 - 6*x* 
    li2spec8*NCi* 
    sqrtxz2i*z2 + 6*x*lx*lspec5*NCi* 
    sqrtxz2i*z2 - 6*x*lx*lspec6*NCi* 
    sqrtxz2i*z2 - 8*NC*Li2mx*xi*z2 +  
   8*NC*Li2x*xi*z2 + 4*NC*li2spec19*xi*z2 +  
   4*NC*li2spec20*xi*z2 -  
   16*NC*lx*xi*z2 +  
   8*NC*lomx*lx*xi*z2 -  
   8*NC*lx*lopx*xi*z2 -  
   4*NC*lx*lz*xi*z2 +  
   4*NC*lx*lopxz*xi*z2 +  
   4*NC*lz*lopxz*xi*z2 +  
   4*NC*lx*lopxzi*xi*z2 -  
   4*NC*lz*lopxzi*xi*z2 +  
   8*Li2mx*NCi*xi*z2 -  
   8*Li2x*NCi*xi*z2 -  
   4*li2spec19*NCi*xi*z2 -  
   4*li2spec20*NCi*xi*z2 +  
   16*lx*NCi*xi*z2 -  
   8*lomx*lx*NCi*xi*z2 +  
   8*lx*lopx*NCi*xi*z2 +  
   4*lx*lz*NCi*xi*z2 -  
   4*lx*lopxz*NCi*xi*z2 -  
   4*lz*lopxz*NCi*xi*z2 -  
   4*lx*lopxzi*NCi*xi*z2 +  
   4*lz*lopxzi*NCi*xi*z2 -  
   (4*NC*pi2*xi*z2)/3. +  
   (4*NCi*pi2*xi*z2)/3. +  
   8*NC*Li2mx*opxi*z2 - 8*NC*Li2x*opxi*z2 -  
   4*NC*li2spec19*opxi*z2 -  
   4*NC*li2spec20*opxi*z2 +  
   16*NC*lx*opxi*z2 -  
   8*NC*lomx*lx*opxi*z2 +  
   8*NC*lx*lopx*opxi*z2 +  
   4*NC*lx*lz*opxi*z2 -  
   4*NC*lx*lopxz*opxi*z2 -  
   4*NC*lz*lopxz*opxi*z2 -  
   4*NC*lx*lopxzi*opxi*z2 +  
   4*NC*lz*lopxzi*opxi*z2 -  
   8*Li2mx*NCi*opxi*z2 +  
   8*Li2x*NCi*opxi*z2 +  
   4*li2spec19*NCi*opxi*z2 +  
   4*li2spec20*NCi*opxi*z2 -  
   16*lx*NCi*opxi*z2 +  
   8*lomx*lx*NCi*opxi*z2 -  
   8*lx*lopx*NCi*opxi*z2 -  
   4*lx*lz*NCi*opxi*z2 +  
   4*lx*lopxz*NCi*opxi*z2 +  
   4*lz*lopxz*NCi*opxi*z2 +  
   4*lx*lopxzi*NCi*opxi*z2 -  
   4*lz*lopxzi*NCi*opxi*z2 +  
   (4*NC*pi2*opxi*z2)/3. -  
   (4*NCi*pi2*opxi*z2)/3. +  
   8*NC*Li2mx*xi*opxi*z2 -  
   8*NC*Li2x*xi*opxi*z2 -  
   4*NC*li2spec19*xi*opxi*z2 -  
   4*NC*li2spec20*xi*opxi*z2 +  
   16*NC*lx*xi*opxi*z2 -  
   8*NC*lomx*lx*xi*opxi*z2 +  
   8*NC*lx*lopx*xi*opxi*z2 +  
   4*NC*lx*lz*xi*opxi*z2 -  
   4*NC*lx*lopxz*xi*opxi*z2 -  
   4*NC*lz*lopxz*xi*opxi*z2 -  
   4*NC*lx*lopxzi*xi*opxi*z2 +  
   4*NC*lz*lopxzi*xi*opxi*z2 -  
   8*Li2mx*NCi*xi*opxi*z2 +  
   8*Li2x*NCi*xi*opxi*z2 +  
   4*li2spec19*NCi*xi*opxi*z2 +  
   4*li2spec20*NCi*xi*opxi*z2 -  
   16*lx*NCi*xi*opxi*z2 +  
   8*lomx*lx*NCi*xi*opxi*z2 -  
   8*lx*lopx*NCi*xi*opxi*z2 -  
   4*lx*lz*NCi*xi*opxi*z2 +  
   4*lx*lopxz*NCi*xi*opxi*z2 +  
   4*lz*lopxz*NCi*xi*opxi*z2 +  
   4*lx*lopxzi*NCi*xi*opxi* 
    z2 - 4*lz*lopxzi*NCi*xi* 
    opxi*z2 + (4*NC*pi2*xi*opxi* 
      z2)/3. - (4*NCi*pi2*xi*opxi* 
      z2)/3. - 16*NC*x*z2*l22 +  
   16*x*NCi*z2*l22 + 2*NC*x*lx2 -  
   6*NC*x*z*lx2 - 2*x*NCi*lx2 +  
   6*x*z*NCi*lx2 - 3*NC*z*xi*lx2 +  
   3*z*NCi*xi*lx2 +  
   3*NC*z*opxi*lx2 -  
   3*z*NCi*opxi*lx2 +  
   3*NC*z*xi*opxi*lx2 -  
   3*z*NCi*xi*opxi*lx2 +  
   4*NC*x*z2*lx2 - 4*x*NCi*z2*lx2 +  
   6*NC*xi*z2*lx2 -  
   6*NCi*xi*z2*lx2 -  
   6*NC*opxi*z2*lx2 +  
   6*NCi*opxi*z2*lx2 -  
   6*NC*xi*opxi*z2*lx2 +  
   6*NCi*xi*opxi*z2*lx2 +  
   4*NC*x*z*lspec1_2 -  
   4*x*z*NCi*lspec1_2 -  
   8*NC*x*z2*lspec1_2 +  
   8*x*NCi*z2*lspec1_2 + NC*x*lz2 -  
   4*NC*x*z*lz2 - x*NCi*lz2 +  
   4*x*z*NCi*lz2 + NC*z*xi*lz2 -  
   z*NCi*xi*lz2 - NC*z*opxi*lz2 +  
   z*NCi*opxi*lz2 -  
   NC*z*xi*opxi*lz2 +  
   z*NCi*xi*opxi*lz2 +  
   4*NC*x*z2*lz2 - 4*x*NCi*z2*lz2 -  
   2*NC*xi*z2*lz2 +  
   2*NCi*xi*z2*lz2 +  
   2*NC*opxi*z2*lz2 -  
   2*NCi*opxi*z2*lz2 +  
   2*NC*xi*opxi*z2*lz2 -  
   2*NCi*xi*opxi*z2*lz2;
};
double C2TQ2G_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double sqrtxz1  = sqrt(1 - 2*z + z*z + 4*x*z);
  const double NC   = apfel::NC;
  const double NC2  = NC * NC;
  const double NCi  = 1. / NC;
  const double NCi2 = NCi * NCi;
  const double l2  = log(2);
  const double l22 = l2 * l2;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double x2 = x * x;
  const double x3 = x * x2;
  const double x4 = x * x3;
  const double xi  = 1. / x;
  const double z2 = z * z;
  const double zi  = 1. / z;
  const double omxi = 1. / ( 1 - x );
  const double omzi  = 1. / ( 1 - z );
  const double lx  = log(x);
  const double lx2 = lx * lx;
  const double lz  = log(z);
  const double lz2 = lz * lz;
  const double lomx  = log(1 - x);
  const double lomx2 = lomx * lomx;
  const double lomz  = log(1 - z);
  const double lomz2 = lomz * lomz;
  const double lopz  = log(1 + z);
  const double Li2x  = apfel::dilog(x);
  const double Li2z  = apfel::dilog(z);
  const double Li2mz = apfel::dilog(-z);
  const double omxmzi  = 1. / ( 1 - x - z );
  const double omxmzi2 = omxmzi * omxmzi;
  const double lomxmz  = log(1 - x - z);
  const double lmopxpz = log(-1 + x + z);
  const double lxmz    = log(x - z);
  const double lmxpz   = log(-x + z);
  const double li2omxzi = apfel::dilog(1 - x*zi);
  const double lspec1   = log(1 + sqrtxz1 - z);
  const double lspec1_2 = lspec1 * lspec1;
  const double lspec2   = log(1 + sqrtxz1 + z);
  const double li2spec1  = apfel::dilog(0.5 - sqrtxz1/2. - z/2.);
  const double li2spec2  = apfel::dilog(0.5 - sqrtxz1/2. + z/2.);
  const double li2spec3  = apfel::dilog(0.5 - zi/2. - (sqrtxz1*zi)/2.);
  const double li2spec4  = apfel::dilog(0.5 + zi/2. - (sqrtxz1*zi)/2.);
  const double li2spec9  = apfel::dilog((1 - z)*omxi);
  const double li2spec10 = apfel::dilog(x*(1 - z)*omxi*zi);
  const double li2spec11 = apfel::dilog((1 - x)*omzi);
  const double li2spec12 = apfel::dilog((1 - x)*z*xi*omzi);
  const double li2spec13 = apfel::dilog(z*omxi);
  const double li2spec14 = apfel::dilog((1 - z)*z*omxi*xi);
  const double li2spec15 = apfel::dilog(x*z*omxi*omzi);
  const double li2spec16 = apfel::dilog((1 - x)*zi);
  const double li2spec17 = apfel::dilog((1 - x)*(1 - z)*xi*zi);
  const double li2spec18 = apfel::dilog((1 - x)*x*omzi*zi);
  const double Tt1 = (z > x ? 1 : 0);
  const double Tt2 = (z < x ? 1 : 0);
  const double Tu1 = (z < 1 - x && z < x ? 1 : 0);
  const double Tu2 = (z > 1 - x && z < x ? 1 : 0);
  const double Tu3 = (z < 1 - x && z > x ? 1 : 0);
  const double Tu4 = (z > 1 - x && z > x ? 1 : 0);
  return 11.916666666666666 + (53*x)/12. - (11*z)/12. - (35*x*z)/12. +  
   li2spec1 + x*li2spec1 +  
   x*z*li2spec1 - li2spec2 -  
   x*li2spec2 - x*z*li2spec2 +  
   Li2z + x*Li2z + 5*x*z*Li2z -  
   li2spec3 -  
   x*li2spec3 -  
   x*z*li2spec3 +  
   li2spec4 +  
   x*li2spec4 +  
   x*z*li2spec4 +  
   2*li2omxzi + 2*x*li2omxzi -  
   2*x*z*li2omxzi - sqrtxz1*x*l2 + (7*lomx)/2. +  
   (3*x*lomx)/2. + 2*z*lomx + (7*x*z*lomx)/2. -  
   (x*lx)/2. - (sqrtxz1*x*lx)/2. + 2*(1 - x)*x*lx -  
   (7*z*lx)/2. - 7*x*z*lx - l2*lx - x*l2*lx -  
   x*z*l2*lx + 13*lomx*lx + 23*x*lomx*lx -  
   z*lomx*lx - 17*x*z*lomx*lx + (7*lomz)/2. +  
   3*x*lomz - 2*(1 - x)*x*lomz + (7*z*lomz)/2. +  
   2*x*z*lomz - 9*lomx*lomz - 15*x*lomx*lomz +  
   z*lomx*lomz + 11*x*z*lomx*lomz +  
   13*lx*lomz + 23*x*lx*lomz - z*lx*lomz -  
   17*x*z*lx*lomz + sqrtxz1*x*lspec1 +  
   2*l2*lspec1 + 2*x*l2*lspec1 +  
   2*x*z*l2*lspec1 - 9*lz - (5*x*lz)/2. -  
   (sqrtxz1*x*lz)/2. - (3*z*lz)/2. + 11*x*z*lz -  
   3*l2*lz - 3*x*l2*lz - 3*x*z*l2*lz -  
   10*lomx*lz - 12*x*lomx*lz + 5*x*z*lomx*lz +  
   14*lx*lz + 18*x*lx*lz - 8*x*z*lx*lz -  
   9*lomz*lz - 11*x*lomz*lz + 10*x*z*lomz*lz +  
   lspec1*lz + x*lspec1*lz +  
   x*z*lspec1*lz + 2*l2*lspec2 +  
   2*x*l2*lspec2 + 2*x*z*l2*lspec2 +  
   lx*lspec2 + x*lx*lspec2 +  
   x*z*lx*lspec2 -  
   2*lspec1*lspec2 -  
   2*x*lspec1*lspec2 -  
   2*x*z*lspec1*lspec2 +  
   2*lz*lspec2 + 2*x*lz*lspec2 +  
   2*x*z*lz*lspec2 - (41*NCi2)/8. -  
   (23*x*NCi2)/8. + (z*NCi2)/4. + (11*x*z*NCi2)/4. +  
   (Li2x*NCi2)/2. + (x*Li2x*NCi2)/2. -  
   (x*z*Li2x*NCi2)/2. - x*Li2z*NCi2 +  
   (x*z*Li2z*NCi2)/2. - li2omxzi*NCi2 -  
   x*li2omxzi*NCi2 + x*z*li2omxzi*NCi2 -  
   sqrtxz1*x*l2*NCi2 - 2*lomx*NCi2 +  
   (3*x*lomx*NCi2)/2. - (z*lomx*NCi2)/4. -  
   2*x*z*lomx*NCi2 - (lx*NCi2)/2. -  
   2*x*lx*NCi2 - (sqrtxz1*x*lx*NCi2)/2. -  
   2*(1 - x)*x*lx*NCi2 + (z*lx*NCi2)/4. +  
   (15*x*z*lx*NCi2)/4. - 5*lomx*lx*NCi2 -  
   15*x*lomx*lx*NCi2 + (z*lomx*lx*NCi2)/2. +  
   (19*x*z*lomx*lx*NCi2)/2. - (lomz*NCi2)/2. +  
   2*(1 - x)*x*lomz*NCi2 - z*lomz*NCi2 -  
   (9*x*z*lomz*NCi2)/4. + 4*lomx*lomz*NCi2 +  
   10*x*lomx*lomz*NCi2 -  
   (z*lomx*lomz*NCi2)/2. -  
   (13*x*z*lomx*lomz*NCi2)/2. -  
   (13*lx*lomz*NCi2)/2. -  
   (33*x*lx*lomz*NCi2)/2. +  
   (z*lx*lomz*NCi2)/2. +  
   11*x*z*lx*lomz*NCi2 +  
   sqrtxz1*x*lspec1*NCi2 + 2*lz*NCi2 +  
   (7*x*lz*NCi2)/2. - (sqrtxz1*x*lz*NCi2)/2. +  
   (z*lz*NCi2)/2. - (9*x*z*lz*NCi2)/2. +  
   3*lomx*lz*NCi2 + 6*x*lomx*lz*NCi2 -  
   (9*x*z*lomx*lz*NCi2)/2. - 5*lx*lz*NCi2 -  
   10*x*lx*lz*NCi2 + (17*x*z*lx*lz*NCi2)/2. +  
   4*lomz*lz*NCi2 + 6*x*lomz*lz*NCi2 -  
   5*x*z*lomz*lz*NCi2 - (163*NC2)/24. -  
   (37*x*NC2)/24. + (2*z*NC2)/3. + (x*z*NC2)/6. -  
   (Li2x*NC2)/2. - (x*Li2x*NC2)/2. +  
   (x*z*Li2x*NC2)/2. - li2spec1*NC2 -  
   x*li2spec1*NC2 -  
   x*z*li2spec1*NC2 +  
   li2spec2*NC2 +  
   x*li2spec2*NC2 +  
   x*z*li2spec2*NC2 - Li2z*NC2 -  
   (11*x*z*Li2z*NC2)/2. +  
   li2spec3*NC2 +  
   x*li2spec3*NC2 +  
   x*z*li2spec3*NC2 -  
   li2spec4*NC2 -  
   x*li2spec4*NC2 -  
   x*z*li2spec4*NC2 -  
   li2omxzi*NC2 - x*li2omxzi*NC2 +  
   x*z*li2omxzi*NC2 + 2*sqrtxz1*x*l2*NC2 -  
   (3*lomx*NC2)/2. - 3*x*lomx*NC2 -  
   (7*z*lomx*NC2)/4. - (3*x*z*lomx*NC2)/2. +  
   (lx*NC2)/2. + (5*x*lx*NC2)/2. +  
   sqrtxz1*x*lx*NC2 + (13*z*lx*NC2)/4. +  
   (13*x*z*lx*NC2)/4. + l2*lx*NC2 +  
   x*l2*lx*NC2 + x*z*l2*lx*NC2 -  
   8*lomx*lx*NC2 - 8*x*lomx*lx*NC2 +  
   (z*lomx*lx*NC2)/2. +  
   (15*x*z*lomx*lx*NC2)/2. - 3*lomz*NC2 -  
   3*x*lomz*NC2 - (5*z*lomz*NC2)/2. +  
   (x*z*lomz*NC2)/4. + 5*lomx*lomz*NC2 +  
   5*x*lomx*lomz*NC2 -  
   (z*lomx*lomz*NC2)/2. -  
   (9*x*z*lomx*lomz*NC2)/2. -  
   (13*lx*lomz*NC2)/2. -  
   (13*x*lx*lomz*NC2)/2. +  
   (z*lx*lomz*NC2)/2. + 6*x*z*lx*lomz*NC2 -  
   2*sqrtxz1*x*lspec1*NC2 -  
   2*l2*lspec1*NC2 -  
   2*x*l2*lspec1*NC2 -  
   2*x*z*l2*lspec1*NC2 + 7*lz*NC2 -  
   x*lz*NC2 + sqrtxz1*x*lz*NC2 + z*lz*NC2 -  
   (13*x*z*lz*NC2)/2. + 3*l2*lz*NC2 +  
   3*x*l2*lz*NC2 + 3*x*z*l2*lz*NC2 +  
   7*lomx*lz*NC2 + 6*x*lomx*lz*NC2 -  
   (x*z*lomx*lz*NC2)/2. - 9*lx*lz*NC2 -  
   8*x*lx*lz*NC2 - (x*z*lx*lz*NC2)/2. +  
   5*lomz*lz*NC2 + 5*x*lomz*lz*NC2 -  
   5*x*z*lomz*lz*NC2 -  
   lspec1*lz*NC2 -  
   x*lspec1*lz*NC2 -  
   x*z*lspec1*lz*NC2 -  
   2*l2*lspec2*NC2 -  
   2*x*l2*lspec2*NC2 -  
   2*x*z*l2*lspec2*NC2 -  
   lx*lspec2*NC2 -  
   x*lx*lspec2*NC2 -  
   x*z*lx*lspec2*NC2 +  
   2*lspec1*lspec2*NC2 +  
   2*x*lspec1*lspec2*NC2 +  
   2*x*z*lspec1*lspec2*NC2 -  
   2*lz*lspec2*NC2 -  
   2*x*lz*lspec2*NC2 -  
   2*x*z*lz*lspec2*NC2 + (5*pi2)/6. +  
   (5*x*pi2)/2. - (z*pi2)/6. - (5*x*z*pi2)/2. -  
   (7*NCi2*pi2)/12. - (25*x*NCi2*pi2)/12. +  
   (z*NCi2*pi2)/12. + (5*x*z*NCi2*pi2)/4. -  
   (NC2*pi2)/4. - (5*x*NC2*pi2)/12. +  
   (z*NC2*pi2)/12. + (5*x*z*NC2*pi2)/4. -  
   (5*z*omxi)/2. + Li2x*omxi -  
   (z*Li2x*omxi)/2. -  
   3*li2spec1*omxi +  
   (z*li2spec1*omxi)/2. +  
   3*li2spec2*omxi -  
   (z*li2spec2*omxi)/2. +  
   2*Li2mz*omxi + z*Li2mz*omxi +  
   3*Li2z*omxi - (3*z*Li2z*omxi)/2. -  
   li2spec3*omxi +  
   (3*z*li2spec3*omxi)/2. +  
   li2spec4*omxi -  
   (3*z*li2spec4*omxi)/2. -  
   3*li2omxzi*omxi +  
   (3*z*li2omxzi*omxi)/2. +  
   3*sqrtxz1*l2*omxi + 5*lomx*omxi -  
   5*z*lomx*omxi - 4*lx*omxi +  
   (3*sqrtxz1*lx*omxi)/2. + (51*z*lx*omxi)/4. +  
   7*l2*lx*omxi - (5*z*l2*lx*omxi)/2. -  
   21*lomx*lx*omxi +  
   (21*z*lomx*lx*omxi)/2. + 5*lomz*omxi -  
   5*z*lomz*omxi + 10*lomx*lomz*omxi -  
   5*z*lomx*lomz*omxi -  
   24*lx*lomz*omxi +  
   12*z*lx*lomz*omxi -  
   3*sqrtxz1*lspec1*omxi -  
   10*l2*lspec1*omxi +  
   3*z*l2*lspec1*omxi -  
   4*lx*lspec1*omxi +  
   2*z*lx*lspec1*omxi + 5*lz*omxi +  
   (3*sqrtxz1*lz*omxi)/2. - (13*z*lz*omxi)/2. +  
   5*l2*lz*omxi + (z*l2*lz*omxi)/2. +  
   10*lomx*lz*omxi - 5*z*lomx*lz*omxi -  
   22*lx*lz*omxi + (9*z*lx*lz*omxi)/2. +  
   12*lomz*lz*omxi - 6*z*lomz*lz*omxi -  
   3*lspec1*lz*omxi +  
   (z*lspec1*lz*omxi)/2. +  
   2*lz*lopz*omxi + z*lz*lopz*omxi -  
   2*l2*lspec2*omxi -  
   z*l2*lspec2*omxi -  
   3*lx*lspec2*omxi +  
   (z*lx*lspec2*omxi)/2. +  
   2*lspec1*lspec2*omxi +  
   z*lspec1*lspec2*omxi -  
   2*lz*lspec2*omxi -  
   z*lz*lspec2*omxi +  
   2*z*NCi2*omxi - Li2x*NCi2*omxi +  
   (z*Li2x*NCi2*omxi)/2. +  
   li2spec1*NCi2*omxi -  
   (z*li2spec1*NCi2*omxi)/2. -  
   li2spec2*NCi2*omxi +  
   (z*li2spec2*NCi2*omxi)/2. -  
   3*Li2z*NCi2*omxi +  
   (3*z*Li2z*NCi2*omxi)/2. +  
   li2spec3*NCi2* 
    omxi - (z*li2spec3* 
      NCi2*omxi)/2. -  
   li2spec4*NCi2* 
    omxi + (z*li2spec4* 
      NCi2*omxi)/2. +  
   2*li2omxzi*NCi2*omxi -  
   z*li2omxzi*NCi2*omxi -  
   sqrtxz1*l2*NCi2*omxi -  
   3*lomx*NCi2*omxi +  
   3*z*lomx*NCi2*omxi +  
   6*lx*NCi2*omxi -  
   (sqrtxz1*lx*NCi2*omxi)/2. -  
   (51*z*lx*NCi2*omxi)/8. -  
   3*l2*lx*NCi2*omxi +  
   (3*z*l2*lx*NCi2*omxi)/2. +  
   13*lomx*lx*NCi2*omxi -  
   (13*z*lomx*lx*NCi2*omxi)/2. -  
   4*lomz*NCi2*omxi +  
   4*z*lomz*NCi2*omxi -  
   8*lomx*lomz*NCi2*omxi +  
   4*z*lomx*lomz*NCi2*omxi +  
   17*lx*lomz*NCi2*omxi -  
   (17*z*lx*lomz*NCi2*omxi)/2. +  
   sqrtxz1*lspec1*NCi2*omxi +  
   4*l2*lspec1*NCi2*omxi -  
   2*z*l2*lspec1*NCi2*omxi +  
   2*lx*lspec1*NCi2*omxi -  
   z*lx*lspec1*NCi2*omxi -  
   4*lz*NCi2*omxi -  
   (sqrtxz1*lz*NCi2*omxi)/2. +  
   (7*z*lz*NCi2*omxi)/2. -  
   l2*lz*NCi2*omxi +  
   (z*l2*lz*NCi2*omxi)/2. -  
   6*lomx*lz*NCi2*omxi +  
   3*z*lomx*lz*NCi2*omxi +  
   (25*lx*lz*NCi2*omxi)/2. -  
   (25*z*lx*lz*NCi2*omxi)/4. -  
   10*lomz*lz*NCi2*omxi +  
   5*z*lomz*lz*NCi2*omxi +  
   lspec1*lz*NCi2*omxi -  
   (z*lspec1*lz*NCi2*omxi)/2. +  
   lx*lspec2*NCi2*omxi -  
   (z*lx*lspec2*NCi2*omxi)/2. +  
   (z*NC2*omxi)/2. +  
   2*li2spec1*NC2*omxi -  
   2*li2spec2*NC2*omxi -  
   2*Li2mz*NC2*omxi - z*Li2mz*NC2*omxi -  
   z*li2spec3*NC2* 
    omxi + z*li2spec4* 
    NC2*omxi + li2omxzi*NC2*omxi -  
   (z*li2omxzi*NC2*omxi)/2. -  
   2*sqrtxz1*l2*NC2*omxi -  
   2*lomx*NC2*omxi +  
   2*z*lomx*NC2*omxi -  
   2*lx*NC2*omxi -  
   sqrtxz1*lx*NC2*omxi -  
   (51*z*lx*NC2*omxi)/8. -  
   4*l2*lx*NC2*omxi +  
   z*l2*lx*NC2*omxi +  
   8*lomx*lx*NC2*omxi -  
   4*z*lomx*lx*NC2*omxi -  
   lomz*NC2*omxi +  
   z*lomz*NC2*omxi -  
   2*lomx*lomz*NC2*omxi +  
   z*lomx*lomz*NC2*omxi +  
   7*lx*lomz*NC2*omxi -  
   (7*z*lx*lomz*NC2*omxi)/2. +  
   2*sqrtxz1*lspec1*NC2*omxi +  
   6*l2*lspec1*NC2*omxi -  
   z*l2*lspec1*NC2*omxi +  
   2*lx*lspec1*NC2*omxi -  
   z*lx*lspec1*NC2*omxi -  
   lz*NC2*omxi - sqrtxz1*lz*NC2*omxi +  
   3*z*lz*NC2*omxi -  
   4*l2*lz*NC2*omxi -  
   z*l2*lz*NC2*omxi -  
   4*lomx*lz*NC2*omxi +  
   2*z*lomx*lz*NC2*omxi +  
   (19*lx*lz*NC2*omxi)/2. +  
   (7*z*lx*lz*NC2*omxi)/4. -  
   2*lomz*lz*NC2*omxi +  
   z*lomz*lz*NC2*omxi +  
   2*lspec1*lz*NC2*omxi -  
   2*lz*lopz*NC2*omxi -  
   z*lz*lopz*NC2*omxi +  
   2*l2*lspec2*NC2*omxi +  
   z*l2*lspec2*NC2*omxi +  
   2*lx*lspec2*NC2*omxi -  
   2*lspec1*lspec2*NC2*omxi -  
   z*lspec1*lspec2*NC2*omxi +  
   2*lz*lspec2*NC2*omxi +  
   z*lz*lspec2*NC2*omxi -  
   (7*pi2*omxi)/3. + (4*z*pi2*omxi)/3. +  
   (5*NCi2*pi2*omxi)/3. -  
   (5*z*NCi2*pi2*omxi)/6. +  
   (2*NC2*pi2*omxi)/3. -  
   (z*NC2*pi2*omxi)/2. - x2 + lx*x2 -  
   lomz*x2 + NCi2*x2 - lx*NCi2*x2 +  
   lomz*NCi2*x2 - lx*x2*omxmzi2 +  
   lomz*x2*omxmzi2 +  
   lx*NCi2*x2*omxmzi2 -  
   lomz*NCi2*x2*omxmzi2 +  
   lx*x3*omxmzi2 -  
   lomz*x3*omxmzi2 -  
   lx*NCi2*x3*omxmzi2 +  
   lomz*NCi2*x3*omxmzi2 -  
   lx*x4*omxmzi2 +  
   lomz*x4*omxmzi2 +  
   lx*NCi2*x4*omxmzi2 -  
   lomz*NCi2*x4*omxmzi2 - x*omxmzi +  
   x*NCi2*omxmzi + x2*omxmzi +  
   lx*x2*omxmzi -  
   lomz*x2*omxmzi -  
   NCi2*x2*omxmzi -  
   lx*NCi2*x2*omxmzi +  
   lomz*NCi2*x2*omxmzi -  
   x3*omxmzi - 2*lx*x3*omxmzi +  
   2*lomz*x3*omxmzi +  
   NCi2*x3*omxmzi +  
   2*lx*NCi2*x3*omxmzi -  
   2*lomz*NCi2*x3*omxmzi - (377*zi)/36. +  
   (115*x*zi)/36. + Li2x*zi + x*Li2x*zi -  
   (li2spec1*zi)/2. -  
   (x*li2spec1*zi)/2. +  
   (li2spec2*zi)/2. +  
   (x*li2spec2*zi)/2. + (3*Li2z*zi)/2. +  
   (x*Li2z*zi)/2. - (3* 
      li2spec3*zi)/2. -  
   (3*x*li2spec3*zi)/2. +  
   (3*li2spec4*zi)/2. +  
   (3*x*li2spec4*zi)/2. -  
   li2omxzi*zi - x*li2omxzi*zi +  
   2*sqrtxz1*l2*zi + sqrtxz1*x*l2*zi -  
   (79*lomx*zi)/12. + (5*x*lomx*zi)/12. +  
   (83*lx*zi)/12. + sqrtxz1*lx*zi -  
   (31*x*lx*zi)/12. + (sqrtxz1*x*lx*zi)/2. +  
   (5*l2*lx*zi)/2. + (5*x*l2*lx*zi)/2. -  
   (13*lomx*lx*zi)/2. -  
   (23*x*lomx*lx*zi)/2. - (79*lomz*zi)/12. +  
   (5*x*lomz*zi)/12. + (11*lomx*lomz*zi)/2. +  
   (17*x*lomx*lomz*zi)/2. -  
   (15*lx*lomz*zi)/2. -  
   (25*x*lx*lomz*zi)/2. -  
   2*sqrtxz1*lspec1*zi -  
   sqrtxz1*x*lspec1*zi -  
   3*l2*lspec1*zi -  
   3*x*l2*lspec1*zi -  
   2*lx*lspec1*zi -  
   2*x*lx*lspec1*zi - (59*lz*zi)/6. +  
   sqrtxz1*lz*zi + (2*x*lz*zi)/3. +  
   (sqrtxz1*x*lz*zi)/2. - (l2*lz*zi)/2. -  
   (x*l2*lz*zi)/2. + 3*lomx*lz*zi +  
   5*x*lomx*lz*zi - 4*lx*lz*zi -  
   7*x*lx*lz*zi + (9*lomz*lz*zi)/2. +  
   (11*x*lomz*lz*zi)/2. -  
   (lspec1*lz*zi)/2. -  
   (x*lspec1*lz*zi)/2. +  
   l2*lspec2*zi +  
   x*l2*lspec2*zi -  
   (lx*lspec2*zi)/2. -  
   (x*lx*lspec2*zi)/2. -  
   lspec1*lspec2*zi -  
   x*lspec1*lspec2*zi +  
   lz*lspec2*zi +  
   x*lz*lspec2*zi + (9*NCi2*zi)/2. -  
   2*x*NCi2*zi - (3*Li2x*NCi2*zi)/4. -  
   (3*x*Li2x*NCi2*zi)/4. +  
   (li2spec1*NCi2*zi)/2. +  
   (x*li2spec1*NCi2*zi)/2. -  
   (li2spec2*NCi2*zi)/2. -  
   (x*li2spec2*NCi2*zi)/2. +  
   x*Li2z*NCi2*zi +  
   (li2spec3*NCi2*zi)/ 
    2. + (x*li2spec3*NCi2* 
      zi)/2. - (li2spec4* 
      NCi2*zi)/2. -  
   (x*li2spec4*NCi2*zi)/ 
    2. + (li2omxzi*NCi2*zi)/2. +  
   (x*li2omxzi*NCi2*zi)/2. -  
   2*sqrtxz1*l2*NCi2*zi -  
   sqrtxz1*x*l2*NCi2*zi +  
   (11*lomx*NCi2*zi)/8. -  
   (17*x*lomx*NCi2*zi)/8. -  
   (lx*NCi2*zi)/8. - sqrtxz1*lx*NCi2*zi +  
   (43*x*lx*NCi2*zi)/8. -  
   (sqrtxz1*x*lx*NCi2*zi)/2. -  
   (3*l2*lx*NCi2*zi)/2. -  
   (3*x*l2*lx*NCi2*zi)/2. +  
   (5*lomx*lx*NCi2*zi)/2. +  
   (15*x*lomx*lx*NCi2*zi)/2. +  
   (11*lomz*NCi2*zi)/8. -  
   (25*x*lomz*NCi2*zi)/8. -  
   (5*lomx*lomz*NCi2*zi)/2. -  
   (11*x*lomx*lomz*NCi2*zi)/2. +  
   (15*lx*lomz*NCi2*zi)/4. +  
   (35*x*lx*lomz*NCi2*zi)/4. +  
   2*sqrtxz1*lspec1*NCi2*zi +  
   sqrtxz1*x*lspec1*NCi2*zi +  
   2*l2*lspec1*NCi2*zi +  
   2*x*l2*lspec1*NCi2*zi +  
   lx*lspec1*NCi2*zi +  
   x*lx*lspec1*NCi2*zi +  
   (5*lz*NCi2*zi)/2. - sqrtxz1*lz*NCi2*zi -  
   (5*x*lz*NCi2*zi)/2. -  
   (sqrtxz1*x*lz*NCi2*zi)/2. -  
   (l2*lz*NCi2*zi)/2. -  
   (x*l2*lz*NCi2*zi)/2. -  
   (3*lomx*lz*NCi2*zi)/2. -  
   (7*x*lomx*lz*NCi2*zi)/2. +  
   (5*lx*lz*NCi2*zi)/2. +  
   (11*x*lx*lz*NCi2*zi)/2. -  
   2*lomz*lz*NCi2*zi -  
   3*x*lomz*lz*NCi2*zi +  
   (lspec1*lz*NCi2*zi)/2. +  
   (x*lspec1*lz*NCi2*zi)/2. +  
   (lx*lspec2*NCi2*zi)/2. +  
   (x*lx*lspec2*NCi2*zi)/2. +  
   (215*NC2*zi)/36. - (43*x*NC2*zi)/36. -  
   (Li2x*NC2*zi)/4. - (x*Li2x*NC2*zi)/4. -  
   (3*Li2z*NC2*zi)/2. - (3*x*Li2z*NC2*zi)/2. +  
   li2spec3*NC2*zi +  
   x*li2spec3*NC2*zi -  
   li2spec4*NC2*zi -  
   x*li2spec4*NC2*zi +  
   (li2omxzi*NC2*zi)/2. +  
   (x*li2omxzi*NC2*zi)/2. +  
   (125*lomx*NC2*zi)/24. +  
   (41*x*lomx*NC2*zi)/24. -  
   (163*lx*NC2*zi)/24. -  
   (67*x*lx*NC2*zi)/24. -  
   l2*lx*NC2*zi - x*l2*lx*NC2*zi +  
   4*lomx*lx*NC2*zi +  
   4*x*lomx*lx*NC2*zi +  
   (125*lomz*NC2*zi)/24. +  
   (65*x*lomz*NC2*zi)/24. -  
   3*lomx*lomz*NC2*zi -  
   3*x*lomx*lomz*NC2*zi +  
   (15*lx*lomz*NC2*zi)/4. +  
   (15*x*lx*lomz*NC2*zi)/4. +  
   l2*lspec1*NC2*zi +  
   x*l2*lspec1*NC2*zi +  
   lx*lspec1*NC2*zi +  
   x*lx*lspec1*NC2*zi +  
   (22*lz*NC2*zi)/3. +  
   (11*x*lz*NC2*zi)/6. + l2*lz*NC2*zi +  
   x*l2*lz*NC2*zi -  
   (3*lomx*lz*NC2*zi)/2. -  
   (3*x*lomx*lz*NC2*zi)/2. +  
   (3*lx*lz*NC2*zi)/2. +  
   (3*x*lx*lz*NC2*zi)/2. -  
   (5*lomz*lz*NC2*zi)/2. -  
   (5*x*lomz*lz*NC2*zi)/2. -  
   l2*lspec2*NC2*zi -  
   x*l2*lspec2*NC2*zi +  
   lspec1*lspec2*NC2*zi +  
   x*lspec1*lspec2*NC2*zi -  
   lz*lspec2*NC2*zi -  
   x*lz*lspec2*NC2*zi -  
   (5*pi2*zi)/4. - (23*x*pi2*zi)/12. +  
   (13*NCi2*pi2*zi)/24. +  
   (29*x*NCi2*pi2*zi)/24. +  
   (17*NC2*pi2*zi)/24. +  
   (17*x*NC2*pi2*zi)/24. -  
   Li2x*omxi*zi +  
   li2spec1*omxi*zi -  
   li2spec2*omxi*zi +  
   2*Li2mz*omxi*zi - 2*Li2z*omxi*zi +  
   3*li2spec3*omxi* 
    zi - 3*li2spec4* 
    omxi*zi + 2*li2omxzi*omxi* 
    zi - 3*sqrtxz1*l2*omxi*zi -  
   5*lomx*omxi*zi +  
   (lx*omxi*zi)/6. -  
   (3*sqrtxz1*lx*omxi*zi)/2. -  
   5*l2*lx*omxi*zi +  
   18*lomx*lx*omxi*zi -  
   5*lomz*omxi*zi -  
   8*lomx*lomz*omxi*zi +  
   19*lx*lomz*omxi*zi +  
   3*sqrtxz1*lspec1*omxi*zi +  
   6*l2*lspec1*omxi*zi +  
   4*lx*lspec1*omxi*zi -  
   (7*lz*omxi*zi)/2. -  
   (3*sqrtxz1*lz*omxi*zi)/2. +  
   l2*lz*omxi*zi -  
   8*lomx*lz*omxi*zi +  
   11*lx*lz*omxi*zi -  
   8*lomz*lz*omxi*zi +  
   lspec1*lz*omxi*zi +  
   2*lz*lopz*omxi*zi -  
   2*l2*lspec2*omxi*zi +  
   lx*lspec2*omxi*zi +  
   2*lspec1*lspec2*omxi*zi -  
   2*lz*lspec2*omxi*zi +  
   Li2x*NCi2*omxi*zi -  
   li2spec1*NCi2*omxi*zi +  
   li2spec2*NCi2*omxi*zi +  
   2*Li2z*NCi2*omxi*zi -  
   li2spec3*NCi2*omxi* 
    zi + li2spec4*NCi2* 
    omxi*zi - li2omxzi*NCi2*omxi* 
    zi + 3*sqrtxz1*l2*NCi2*omxi*zi +  
   3*lomx*NCi2*omxi*zi -  
   6*lx*NCi2*omxi*zi +  
   (3*sqrtxz1*lx*NCi2*omxi*zi)/2. +  
   3*l2*lx*NCi2*omxi*zi -  
   10*lomx*lx*NCi2*omxi*zi +  
   4*lomz*NCi2*omxi*zi +  
   6*lomx*lomz*NCi2*omxi*zi -  
   12*lx*lomz*NCi2*omxi*zi -  
   3*sqrtxz1*lspec1*NCi2*omxi*zi -  
   4*l2*lspec1*NCi2*omxi*zi -  
   2*lx*lspec1*NCi2*omxi*zi +  
   (3*lz*NCi2*omxi*zi)/2. +  
   (3*sqrtxz1*lz*NCi2*omxi*zi)/2. +  
   l2*lz*NCi2*omxi*zi +  
   4*lomx*lz*NCi2*omxi*zi -  
   8*lx*lz*NCi2*omxi*zi +  
   6*lomz*lz*NCi2*omxi*zi -  
   lspec1*lz*NCi2*omxi*zi -  
   lx*lspec2*NCi2*omxi*zi -  
   2*Li2mz*NC2*omxi*zi -  
   2*li2spec3*NC2*omxi* 
    zi + 2*li2spec4*NC2* 
    omxi*zi - li2omxzi*NC2*omxi* 
    zi + 2*lomx*NC2*omxi*zi +  
   (35*lx*NC2*omxi*zi)/6. +  
   2*l2*lx*NC2*omxi*zi -  
   8*lomx*lx*NC2*omxi*zi +  
   lomz*NC2*omxi*zi +  
   2*lomx*lomz*NC2*omxi*zi -  
   7*lx*lomz*NC2*omxi*zi -  
   2*l2*lspec1*NC2*omxi*zi -  
   2*lx*lspec1*NC2*omxi*zi +  
   2*lz*NC2*omxi*zi -  
   2*l2*lz*NC2*omxi*zi +  
   4*lomx*lz*NC2*omxi*zi -  
   3*lx*lz*NC2*omxi*zi +  
   2*lomz*lz*NC2*omxi*zi -  
   2*lz*lopz*NC2*omxi*zi +  
   2*l2*lspec2*NC2*omxi*zi -  
   2*lspec1*lspec2*NC2*omxi* 
    zi + 2*lz*lspec2*NC2*omxi* 
    zi + (7*pi2*omxi*zi)/3. -  
   (4*NCi2*pi2*omxi*zi)/3. -  
   NC2*pi2*omxi*zi + (31*z2)/18. -  
   (53*x*z2)/18. - (2*lomx*z2)/3. +  
   (4*x*lomx*z2)/3. + (4*lx*z2)/3. -  
   (8*x*lx*z2)/3. - (2*lomz*z2)/3. +  
   (4*x*lomz*z2)/3. - (2*lz*z2)/3. +  
   (4*x*lz*z2)/3. - (31*NC2*z2)/18. +  
   (53*x*NC2*z2)/18. + (2*lomx*NC2*z2)/3. -  
   (4*x*lomx*NC2*z2)/3. -  
   (4*lx*NC2*z2)/3. + (8*x*lx*NC2*z2)/3. +  
   (2*lomz*NC2*z2)/3. -  
   (4*x*lomz*NC2*z2)/3. +  
   (2*lz*NC2*z2)/3. - (4*x*lz*NC2*z2)/3. +  
   (4*lx*omxi*z2)/3. -  
   (4*lx*NC2*omxi*z2)/3. - 2*l22 -  
   2*x*l22 - 2*x*z*l22 + 2*NC2*l22 +  
   2*x*NC2*l22 + 2*x*z*NC2*l22 +  
   6*omxi*l22 - z*omxi*l22 -  
   2*NCi2*omxi*l22 +  
   z*NCi2*omxi*l22 -  
   4*NC2*omxi*l22 + zi*l22 +  
   x*zi*l22 - NCi2*zi*l22 -  
   x*NCi2*zi*l22 -  
   2*omxi*zi*l22 +  
   2*NCi2*omxi*zi*l22 -  
   (9*lomx2)/2. - (17*x*lomx2)/2. +  
   (z*lomx2)/2. + 6*x*z*lomx2 +  
   (3*NCi2*lomx2)/2. +  
   (11*x*NCi2*lomx2)/2. -  
   (z*NCi2*lomx2)/4. -  
   (13*x*z*NCi2*lomx2)/4. + 3*NC2*lomx2 +  
   3*x*NC2*lomx2 - (z*NC2*lomx2)/4. -  
   (11*x*z*NC2*lomx2)/4. +  
   5*omxi*lomx2 -  
   (5*z*omxi*lomx2)/2. -  
   3*NCi2*omxi*lomx2 +  
   (3*z*NCi2*omxi*lomx2)/2. -  
   2*NC2*omxi*lomx2 +  
   z*NC2*omxi*lomx2 +  
   (11*zi*lomx2)/4. +  
   (19*x*zi*lomx2)/4. -  
   NCi2*zi*lomx2 -  
   3*x*NCi2*zi*lomx2 -  
   (7*NC2*zi*lomx2)/4. -  
   (7*x*NC2*zi*lomx2)/4. -  
   (9*omxi*zi*lomx2)/2. +  
   (5*NCi2*omxi*zi*lomx2)/2. +  
   2*NC2*omxi*zi*lomx2 - 9*lx2 -  
   16*x*lx2 + (z*lx2)/2. + 12*x*z*lx2 +  
   (9*NCi2*lx2)/2. + (23*x*NCi2*lx2)/2. -  
   (z*NCi2*lx2)/4. - (31*x*z*NCi2*lx2)/4. +  
   (9*NC2*lx2)/2. + (9*x*NC2*lx2)/2. -  
   (z*NC2*lx2)/4. - (17*x*z*NC2*lx2)/4. +  
   17*omxi*lx2 - (17*z*omxi*lx2)/2. -  
   12*NCi2*omxi*lx2 +  
   6*z*NCi2*omxi*lx2 -  
   5*NC2*omxi*lx2 +  
   (5*z*NC2*omxi*lx2)/2. +  
   (11*zi*lx2)/2. + 9*x*zi*lx2 -  
   (11*NCi2*zi*lx2)/4. -  
   (25*x*NCi2*zi*lx2)/4. -  
   (11*NC2*zi*lx2)/4. -  
   (11*x*NC2*zi*lx2)/4. -  
   (27*omxi*zi*lx2)/2. +  
   (17*NCi2*omxi*zi*lx2)/2. +  
   5*NC2*omxi*zi*lx2 -  
   (9*lomz2)/2. - (17*x*lomz2)/2. +  
   (z*lomz2)/2. + 6*x*z*lomz2 +  
   (5*NCi2*lomz2)/2. +  
   (13*x*NCi2*lomz2)/2. -  
   (z*NCi2*lomz2)/4. -  
   (17*x*z*NCi2*lomz2)/4. + 2*NC2*lomz2 +  
   2*x*NC2*lomz2 - (z*NC2*lomz2)/4. -  
   (7*x*z*NC2*lomz2)/4. +  
   (13*omxi*lomz2)/2. -  
   (13*z*omxi*lomz2)/4. -  
   6*NCi2*omxi*lomz2 +  
   3*z*NCi2*omxi*lomz2 -  
   (NC2*omxi*lomz2)/2. +  
   (z*NC2*omxi*lomz2)/4. +  
   (11*zi*lomz2)/4. +  
   (19*x*zi*lomz2)/4. -  
   (3*NCi2*zi*lomz2)/2. -  
   (7*x*NCi2*zi*lomz2)/2. -  
   (5*NC2*zi*lomz2)/4. -  
   (5*x*NC2*zi*lomz2)/4. -  
   (9*omxi*zi*lomz2)/2. +  
   4*NCi2*omxi*zi*lomz2 +  
   (NC2*omxi*zi*lomz2)/2. +  
   4*omxi*lspec1_2 -  
   2*z*omxi*lspec1_2 -  
   2*NCi2*omxi*lspec1_2 +  
   z*NCi2*omxi*lspec1_2 -  
   2*NC2*omxi*lspec1_2 +  
   z*NC2*omxi*lspec1_2 +  
   2*zi*lspec1_2 +  
   2*x*zi*lspec1_2 -  
   NCi2*zi*lspec1_2 -  
   x*NCi2*zi*lspec1_2 -  
   NC2*zi*lspec1_2 -  
   x*NC2*zi*lspec1_2 -  
   4*omxi*zi*lspec1_2 +  
   2*NCi2*omxi*zi*lspec1_2 +  
   2*NC2*omxi*zi*lspec1_2 -  
   (13*lz2)/2. - (11*x*lz2)/2. - (z*lz2)/2. -  
   (3*x*z*lz2)/2. + 2*NCi2*lz2 +  
   2*x*NCi2*lz2 + (z*NCi2*lz2)/4. -  
   (9*x*z*NCi2*lz2)/4. + (9*NC2*lz2)/2. +  
   (7*x*NC2*lz2)/2. + (z*NC2*lz2)/4. +  
   (15*x*z*NC2*lz2)/4. + (9*omxi*lz2)/2. -  
   (7*z*omxi*lz2)/4. -  
   3*NCi2*omxi*lz2 +  
   (3*z*NCi2*omxi*lz2)/2. -  
   (3*NC2*omxi*lz2)/2. +  
   (z*NC2*omxi*lz2)/4. -  
   (3*zi*lz2)/4. - (x*zi*lz2)/4. -  
   (NCi2*zi*lz2)/2. -  
   x*NCi2*zi*lz2 +  
   (5*NC2*zi*lz2)/4. +  
   (5*x*NC2*zi*lz2)/4. -  
   2*omxi*zi*lz2 +  
   (3*NCi2*omxi*zi*lz2)/2. +  
   (NC2*omxi*zi*lz2)/2. +  
   (-1.5 - Li2z - x*Li2z + x*z*Li2z - li2spec9 -  
      x*li2spec9 + x*z*li2spec9 +  
      li2spec10 +  
      x*li2spec10 -  
      x*z*li2spec10 + 3*lomx -  
      2*x*z*lomx - (9*lx)/2. + 3*x*z*lx - 6*lomx*lx -  
      6*x*lomx*lx + 6*x*z*lomx*lx + (3*lomz)/2. -  
      x*z*lomz + 2*lomx*lomz + 2*x*lomx*lomz -  
      2*x*z*lomx*lomz - 3*lx*lomz -  
      3*x*lx*lomz + 3*x*z*lx*lomz + 3*lz -  
      2*x*z*lz + 5*lomx*lz + 5*x*lomx*lz -  
      5*x*z*lomx*lz - 6*lx*lz - 6*x*lx*lz +  
      6*x*z*lx*lz + lomz*lz + x*lomz*lz -  
      x*z*lomz*lz + lx*lmxpz + x*lx*lmxpz -  
      x*z*lx*lmxpz - lz*lmxpz - x*lz*lmxpz +  
      x*z*lz*lmxpz + (3*NC2)/2. + Li2z*NC2 +  
      x*Li2z*NC2 - x*z*Li2z*NC2 +  
      li2spec9*NC2 +  
      x*li2spec9*NC2 -  
      x*z*li2spec9*NC2 -  
      li2spec10*NC2 -  
      x*li2spec10*NC2 +  
      x*z*li2spec10*NC2 -  
      3*lomx*NC2 + 2*x*z*lomx*NC2 +  
      (9*lx*NC2)/2. - 3*x*z*lx*NC2 +  
      6*lomx*lx*NC2 + 6*x*lomx*lx*NC2 -  
      6*x*z*lomx*lx*NC2 - (3*lomz*NC2)/2. +  
      x*z*lomz*NC2 - 2*lomx*lomz*NC2 -  
      2*x*lomx*lomz*NC2 +  
      2*x*z*lomx*lomz*NC2 + 3*lx*lomz*NC2 +  
      3*x*lx*lomz*NC2 - 3*x*z*lx*lomz*NC2 -  
      3*lz*NC2 + 2*x*z*lz*NC2 -  
      5*lomx*lz*NC2 - 5*x*lomx*lz*NC2 +  
      5*x*z*lomx*lz*NC2 + 6*lx*lz*NC2 +  
      6*x*lx*lz*NC2 - 6*x*z*lx*lz*NC2 -  
      lomz*lz*NC2 - x*lomz*lz*NC2 +  
      x*z*lomz*lz*NC2 - lx*lmxpz*NC2 -  
      x*lx*lmxpz*NC2 + x*z*lx*lmxpz*NC2 +  
      lz*lmxpz*NC2 + x*lz*lmxpz*NC2 -  
      x*z*lz*lmxpz*NC2 - pi2/2. - (x*pi2)/2. +  
      (x*z*pi2)/2. + (NC2*pi2)/2. +  
      (x*NC2*pi2)/2. - (x*z*NC2*pi2)/2. +  
      (z*omxi)/2. + Li2z*omxi -  
      (z*Li2z*omxi)/2. +  
      li2spec9*omxi -  
      (z*li2spec9*omxi)/2. -  
      li2spec10*omxi +  
      (z*li2spec10*omxi)/2. -  
      2*lomx*omxi + 2*z*lomx*omxi +  
      3*lx*omxi - 3*z*lx*omxi +  
      6*lomx*lx*omxi -  
      3*z*lomx*lx*omxi - lomz*omxi +  
      z*lomz*omxi - 2*lomx*lomz*omxi +  
      z*lomx*lomz*omxi +  
      3*lx*lomz*omxi -  
      (3*z*lx*lomz*omxi)/2. - 2*lz*omxi +  
      2*z*lz*omxi - 5*lomx*lz*omxi +  
      (5*z*lomx*lz*omxi)/2. +  
      6*lx*lz*omxi - 3*z*lx*lz*omxi -  
      lomz*lz*omxi +  
      (z*lomz*lz*omxi)/2. -  
      lx*lmxpz*omxi +  
      (z*lx*lmxpz*omxi)/2. +  
      lz*lmxpz*omxi -  
      (z*lz*lmxpz*omxi)/2. -  
      (z*NC2*omxi)/2. - Li2z*NC2*omxi +  
      (z*Li2z*NC2*omxi)/2. -  
      li2spec9*NC2*omxi +  
      (z*li2spec9*NC2*omxi)/2. +  
      li2spec10*NC2*omxi -  
      (z*li2spec10*NC2*omxi)/2. +  
      2*lomx*NC2*omxi -  
      2*z*lomx*NC2*omxi -  
      3*lx*NC2*omxi + 3*z*lx*NC2*omxi -  
      6*lomx*lx*NC2*omxi +  
      3*z*lomx*lx*NC2*omxi +  
      lomz*NC2*omxi -  
      z*lomz*NC2*omxi +  
      2*lomx*lomz*NC2*omxi -  
      z*lomx*lomz*NC2*omxi -  
      3*lx*lomz*NC2*omxi +  
      (3*z*lx*lomz*NC2*omxi)/2. +  
      2*lz*NC2*omxi - 2*z*lz*NC2*omxi +  
      5*lomx*lz*NC2*omxi -  
      (5*z*lomx*lz*NC2*omxi)/2. -  
      6*lx*lz*NC2*omxi +  
      3*z*lx*lz*NC2*omxi +  
      lomz*lz*NC2*omxi -  
      (z*lomz*lz*NC2*omxi)/2. +  
      lx*lmxpz*NC2*omxi -  
      (z*lx*lmxpz*NC2*omxi)/2. -  
      lz*lmxpz*NC2*omxi +  
      (z*lz*lmxpz*NC2*omxi)/2. +  
      (pi2*omxi)/2. - (z*pi2*omxi)/4. -  
      (NC2*pi2*omxi)/2. +  
      (z*NC2*pi2*omxi)/4. + zi/2. -  
      (x*zi)/2. + (Li2z*zi)/2. + (x*Li2z*zi)/2. +  
      (li2spec9*zi)/2. +  
      (x*li2spec9*zi)/2. -  
      (li2spec10*zi)/2. -  
      (x*li2spec10*zi)/2. -  
      2*x*lomx*zi + 3*x*lx*zi +  
      3*lomx*lx*zi + 3*x*lomx*lx*zi -  
      x*lomz*zi - lomx*lomz*zi -  
      x*lomx*lomz*zi +  
      (3*lx*lomz*zi)/2. +  
      (3*x*lx*lomz*zi)/2. - 2*x*lz*zi -  
      (5*lomx*lz*zi)/2. -  
      (5*x*lomx*lz*zi)/2. + 3*lx*lz*zi +  
      3*x*lx*lz*zi - (lomz*lz*zi)/2. -  
      (x*lomz*lz*zi)/2. -  
      (lx*lmxpz*zi)/2. -  
      (x*lx*lmxpz*zi)/2. +  
      (lz*lmxpz*zi)/2. +  
      (x*lz*lmxpz*zi)/2. - (NC2*zi)/2. +  
      (x*NC2*zi)/2. - (Li2z*NC2*zi)/2. -  
      (x*Li2z*NC2*zi)/2. -  
      (li2spec9*NC2*zi)/2. -  
      (x*li2spec9*NC2*zi)/2. +  
      (li2spec10*NC2*zi)/2. +  
      (x*li2spec10*NC2*zi)/2. +  
      2*x*lomx*NC2*zi - 3*x*lx*NC2*zi -  
      3*lomx*lx*NC2*zi -  
      3*x*lomx*lx*NC2*zi +  
      x*lomz*NC2*zi +  
      lomx*lomz*NC2*zi +  
      x*lomx*lomz*NC2*zi -  
      (3*lx*lomz*NC2*zi)/2. -  
      (3*x*lx*lomz*NC2*zi)/2. +  
      2*x*lz*NC2*zi +  
      (5*lomx*lz*NC2*zi)/2. +  
      (5*x*lomx*lz*NC2*zi)/2. -  
      3*lx*lz*NC2*zi -  
      3*x*lx*lz*NC2*zi +  
      (lomz*lz*NC2*zi)/2. +  
      (x*lomz*lz*NC2*zi)/2. +  
      (lx*lmxpz*NC2*zi)/2. +  
      (x*lx*lmxpz*NC2*zi)/2. -  
      (lz*lmxpz*NC2*zi)/2. -  
      (x*lz*lmxpz*NC2*zi)/2. +  
      (pi2*zi)/4. + (x*pi2*zi)/4. -  
      (NC2*pi2*zi)/4. -  
      (x*NC2*pi2*zi)/4. - Li2z*omxi*zi -  
      li2spec9*omxi*zi +  
      li2spec10*omxi*zi +  
      2*lomx*omxi*zi -  
      3*lx*omxi*zi -  
      6*lomx*lx*omxi*zi +  
      lomz*omxi*zi +  
      2*lomx*lomz*omxi*zi -  
      3*lx*lomz*omxi*zi +  
      2*lz*omxi*zi +  
      5*lomx*lz*omxi*zi -  
      6*lx*lz*omxi*zi +  
      lomz*lz*omxi*zi +  
      lx*lmxpz*omxi*zi -  
      lz*lmxpz*omxi*zi +  
      Li2z*NC2*omxi*zi +  
      li2spec9*NC2*omxi*zi -  
      li2spec10*NC2*omxi* 
       zi - 2*lomx*NC2*omxi*zi +  
      3*lx*NC2*omxi*zi +  
      6*lomx*lx*NC2*omxi*zi -  
      lomz*NC2*omxi*zi -  
      2*lomx*lomz*NC2*omxi*zi +  
      3*lx*lomz*NC2*omxi*zi -  
      2*lz*NC2*omxi*zi -  
      5*lomx*lz*NC2*omxi*zi +  
      6*lx*lz*NC2*omxi*zi -  
      lomz*lz*NC2*omxi*zi -  
      lx*lmxpz*NC2*omxi*zi +  
      lz*lmxpz*NC2*omxi*zi -  
      (pi2*omxi*zi)/2. +  
      (NC2*pi2*omxi*zi)/2. +  
      2*lomx2 + 2*x*lomx2 - 2*x*z*lomx2 -  
      2*NC2*lomx2 - 2*x*NC2*lomx2 +  
      2*x*z*NC2*lomx2 - 2*omxi*lomx2 +  
      z*omxi*lomx2 +  
      2*NC2*omxi*lomx2 -  
      z*NC2*omxi*lomx2 -  
      zi*lomx2 - x*zi*lomx2 +  
      NC2*zi*lomx2 +  
      x*NC2*zi*lomx2 +  
      2*omxi*zi*lomx2 -  
      2*NC2*omxi*zi*lomx2 +  
      (7*lx2)/2. + (7*x*lx2)/2. -  
      (7*x*z*lx2)/2. - (7*NC2*lx2)/2. -  
      (7*x*NC2*lx2)/2. + (7*x*z*NC2*lx2)/2. -  
      (7*omxi*lx2)/2. +  
      (7*z*omxi*lx2)/4. +  
      (7*NC2*omxi*lx2)/2. -  
      (7*z*NC2*omxi*lx2)/4. -  
      (7*zi*lx2)/4. - (7*x*zi*lx2)/4. +  
      (7*NC2*zi*lx2)/4. +  
      (7*x*NC2*zi*lx2)/4. +  
      (7*omxi*zi*lx2)/2. -  
      (7*NC2*omxi*zi*lx2)/2. +  
      lomz2/2. + (x*lomz2)/2. -  
      (x*z*lomz2)/2. - (NC2*lomz2)/2. -  
      (x*NC2*lomz2)/2. +  
      (x*z*NC2*lomz2)/2. -  
      (omxi*lomz2)/2. +  
      (z*omxi*lomz2)/4. +  
      (NC2*omxi*lomz2)/2. -  
      (z*NC2*omxi*lomz2)/4. -  
      (zi*lomz2)/4. - (x*zi*lomz2)/4. +  
      (NC2*zi*lomz2)/4. +  
      (x*NC2*zi*lomz2)/4. +  
      (omxi*zi*lomz2)/2. -  
      (NC2*omxi*zi*lomz2)/2. +  
      (5*lz2)/2. + (5*x*lz2)/2. -  
      (5*x*z*lz2)/2. - (5*NC2*lz2)/2. -  
      (5*x*NC2*lz2)/2. + (5*x*z*NC2*lz2)/2. -  
      (5*omxi*lz2)/2. +  
      (5*z*omxi*lz2)/4. +  
      (5*NC2*omxi*lz2)/2. -  
      (5*z*NC2*omxi*lz2)/4. -  
      (5*zi*lz2)/4. - (5*x*zi*lz2)/4. +  
      (5*NC2*zi*lz2)/4. +  
      (5*x*NC2*zi*lz2)/4. +  
      (5*omxi*zi*lz2)/2. -  
      (5*NC2*omxi*zi*lz2)/2.)*Tt1 +  
   (-1.5 - Li2z - x*Li2z + x*z*Li2z + li2spec11 +  
      x*li2spec11 - x*z*li2spec11 -  
      li2spec12 -  
      x*li2spec12 +  
      x*z*li2spec12 + 3*lomx -  
      2*x*z*lomx - (9*lx)/2. + 3*x*z*lx - 5*lomx*lx -  
      5*x*lomx*lx + 5*x*z*lomx*lx + (3*lomz)/2. -  
      x*z*lomz + 2*lomx*lomz + 2*x*lomx*lomz -  
      2*x*z*lomx*lomz - 4*lx*lomz -  
      4*x*lx*lomz + 4*x*z*lx*lomz + lx*lxmz +  
      x*lx*lxmz - x*z*lx*lxmz + 3*lz - 2*x*z*lz +  
      4*lomx*lz + 4*x*lomx*lz - 4*x*z*lomx*lz -  
      5*lx*lz - 5*x*lx*lz + 5*x*z*lx*lz +  
      2*lomz*lz + 2*x*lomz*lz - 2*x*z*lomz*lz -  
      lxmz*lz - x*lxmz*lz + x*z*lxmz*lz +  
      (3*NC2)/2. + Li2z*NC2 + x*Li2z*NC2 -  
      x*z*Li2z*NC2 - li2spec11*NC2 -  
      x*li2spec11*NC2 +  
      x*z*li2spec11*NC2 +  
      li2spec12*NC2 +  
      x*li2spec12*NC2 -  
      x*z*li2spec12*NC2 -  
      3*lomx*NC2 + 2*x*z*lomx*NC2 +  
      (9*lx*NC2)/2. - 3*x*z*lx*NC2 +  
      5*lomx*lx*NC2 + 5*x*lomx*lx*NC2 -  
      5*x*z*lomx*lx*NC2 - (3*lomz*NC2)/2. +  
      x*z*lomz*NC2 - 2*lomx*lomz*NC2 -  
      2*x*lomx*lomz*NC2 +  
      2*x*z*lomx*lomz*NC2 + 4*lx*lomz*NC2 +  
      4*x*lx*lomz*NC2 - 4*x*z*lx*lomz*NC2 -  
      lx*lxmz*NC2 - x*lx*lxmz*NC2 +  
      x*z*lx*lxmz*NC2 - 3*lz*NC2 +  
      2*x*z*lz*NC2 - 4*lomx*lz*NC2 -  
      4*x*lomx*lz*NC2 + 4*x*z*lomx*lz*NC2 +  
      5*lx*lz*NC2 + 5*x*lx*lz*NC2 -  
      5*x*z*lx*lz*NC2 - 2*lomz*lz*NC2 -  
      2*x*lomz*lz*NC2 + 2*x*z*lomz*lz*NC2 +  
      lxmz*lz*NC2 + x*lxmz*lz*NC2 -  
      x*z*lxmz*lz*NC2 - pi2/2. - (x*pi2)/2. +  
      (x*z*pi2)/2. + (NC2*pi2)/2. +  
      (x*NC2*pi2)/2. - (x*z*NC2*pi2)/2. +  
      (z*omxi)/2. + Li2z*omxi -  
      (z*Li2z*omxi)/2. -  
      li2spec11*omxi +  
      (z*li2spec11*omxi)/2. +  
      li2spec12*omxi -  
      (z*li2spec12*omxi)/2. -  
      2*lomx*omxi + 2*z*lomx*omxi +  
      3*lx*omxi - 3*z*lx*omxi +  
      5*lomx*lx*omxi -  
      (5*z*lomx*lx*omxi)/2. - lomz*omxi +  
      z*lomz*omxi - 2*lomx*lomz*omxi +  
      z*lomx*lomz*omxi +  
      4*lx*lomz*omxi -  
      2*z*lx*lomz*omxi - lx*lxmz*omxi +  
      (z*lx*lxmz*omxi)/2. - 2*lz*omxi +  
      2*z*lz*omxi - 4*lomx*lz*omxi +  
      2*z*lomx*lz*omxi + 5*lx*lz*omxi -  
      (5*z*lx*lz*omxi)/2. -  
      2*lomz*lz*omxi + z*lomz*lz*omxi +  
      lxmz*lz*omxi -  
      (z*lxmz*lz*omxi)/2. -  
      (z*NC2*omxi)/2. - Li2z*NC2*omxi +  
      (z*Li2z*NC2*omxi)/2. +  
      li2spec11*NC2*omxi -  
      (z*li2spec11*NC2*omxi)/2. -  
      li2spec12*NC2*omxi +  
      (z*li2spec12*NC2*omxi)/2. +  
      2*lomx*NC2*omxi -  
      2*z*lomx*NC2*omxi -  
      3*lx*NC2*omxi + 3*z*lx*NC2*omxi -  
      5*lomx*lx*NC2*omxi +  
      (5*z*lomx*lx*NC2*omxi)/2. +  
      lomz*NC2*omxi -  
      z*lomz*NC2*omxi +  
      2*lomx*lomz*NC2*omxi -  
      z*lomx*lomz*NC2*omxi -  
      4*lx*lomz*NC2*omxi +  
      2*z*lx*lomz*NC2*omxi +  
      lx*lxmz*NC2*omxi -  
      (z*lx*lxmz*NC2*omxi)/2. +  
      2*lz*NC2*omxi - 2*z*lz*NC2*omxi +  
      4*lomx*lz*NC2*omxi -  
      2*z*lomx*lz*NC2*omxi -  
      5*lx*lz*NC2*omxi +  
      (5*z*lx*lz*NC2*omxi)/2. +  
      2*lomz*lz*NC2*omxi -  
      z*lomz*lz*NC2*omxi -  
      lxmz*lz*NC2*omxi +  
      (z*lxmz*lz*NC2*omxi)/2. +  
      (pi2*omxi)/2. - (z*pi2*omxi)/4. -  
      (NC2*pi2*omxi)/2. +  
      (z*NC2*pi2*omxi)/4. + zi/2. -  
      (x*zi)/2. + (Li2z*zi)/2. + (x*Li2z*zi)/2. -  
      (li2spec11*zi)/2. -  
      (x*li2spec11*zi)/2. +  
      (li2spec12*zi)/2. +  
      (x*li2spec12*zi)/2. -  
      2*x*lomx*zi + 3*x*lx*zi +  
      (5*lomx*lx*zi)/2. +  
      (5*x*lomx*lx*zi)/2. - x*lomz*zi -  
      lomx*lomz*zi - x*lomx*lomz*zi +  
      2*lx*lomz*zi + 2*x*lx*lomz*zi -  
      (lx*lxmz*zi)/2. - (x*lx*lxmz*zi)/2. -  
      2*x*lz*zi - 2*lomx*lz*zi -  
      2*x*lomx*lz*zi + (5*lx*lz*zi)/2. +  
      (5*x*lx*lz*zi)/2. - lomz*lz*zi -  
      x*lomz*lz*zi + (lxmz*lz*zi)/2. +  
      (x*lxmz*lz*zi)/2. - (NC2*zi)/2. +  
      (x*NC2*zi)/2. - (Li2z*NC2*zi)/2. -  
      (x*Li2z*NC2*zi)/2. +  
      (li2spec11*NC2*zi)/2. +  
      (x*li2spec11*NC2*zi)/2. -  
      (li2spec12*NC2*zi)/2. -  
      (x*li2spec12*NC2*zi)/2. +  
      2*x*lomx*NC2*zi - 3*x*lx*NC2*zi -  
      (5*lomx*lx*NC2*zi)/2. -  
      (5*x*lomx*lx*NC2*zi)/2. +  
      x*lomz*NC2*zi +  
      lomx*lomz*NC2*zi +  
      x*lomx*lomz*NC2*zi -  
      2*lx*lomz*NC2*zi -  
      2*x*lx*lomz*NC2*zi +  
      (lx*lxmz*NC2*zi)/2. +  
      (x*lx*lxmz*NC2*zi)/2. +  
      2*x*lz*NC2*zi +  
      2*lomx*lz*NC2*zi +  
      2*x*lomx*lz*NC2*zi -  
      (5*lx*lz*NC2*zi)/2. -  
      (5*x*lx*lz*NC2*zi)/2. +  
      lomz*lz*NC2*zi +  
      x*lomz*lz*NC2*zi -  
      (lxmz*lz*NC2*zi)/2. -  
      (x*lxmz*lz*NC2*zi)/2. +  
      (pi2*zi)/4. + (x*pi2*zi)/4. -  
      (NC2*pi2*zi)/4. -  
      (x*NC2*pi2*zi)/4. - Li2z*omxi*zi +  
      li2spec11*omxi*zi -  
      li2spec12*omxi*zi +  
      2*lomx*omxi*zi -  
      3*lx*omxi*zi -  
      5*lomx*lx*omxi*zi +  
      lomz*omxi*zi +  
      2*lomx*lomz*omxi*zi -  
      4*lx*lomz*omxi*zi +  
      lx*lxmz*omxi*zi +  
      2*lz*omxi*zi +  
      4*lomx*lz*omxi*zi -  
      5*lx*lz*omxi*zi +  
      2*lomz*lz*omxi*zi -  
      lxmz*lz*omxi*zi +  
      Li2z*NC2*omxi*zi -  
      li2spec11*NC2*omxi*zi +  
      li2spec12*NC2*omxi* 
       zi - 2*lomx*NC2*omxi*zi +  
      3*lx*NC2*omxi*zi +  
      5*lomx*lx*NC2*omxi*zi -  
      lomz*NC2*omxi*zi -  
      2*lomx*lomz*NC2*omxi*zi +  
      4*lx*lomz*NC2*omxi*zi -  
      lx*lxmz*NC2*omxi*zi -  
      2*lz*NC2*omxi*zi -  
      4*lomx*lz*NC2*omxi*zi +  
      5*lx*lz*NC2*omxi*zi -  
      2*lomz*lz*NC2*omxi*zi +  
      lxmz*lz*NC2*omxi*zi -  
      (pi2*omxi*zi)/2. +  
      (NC2*pi2*omxi*zi)/2. +  
      2*lomx2 + 2*x*lomx2 - 2*x*z*lomx2 -  
      2*NC2*lomx2 - 2*x*NC2*lomx2 +  
      2*x*z*NC2*lomx2 - 2*omxi*lomx2 +  
      z*omxi*lomx2 +  
      2*NC2*omxi*lomx2 -  
      z*NC2*omxi*lomx2 -  
      zi*lomx2 - x*zi*lomx2 +  
      NC2*zi*lomx2 +  
      x*NC2*zi*lomx2 +  
      2*omxi*zi*lomx2 -  
      2*NC2*omxi*zi*lomx2 +  
      3*lx2 + 3*x*lx2 - 3*x*z*lx2 -  
      3*NC2*lx2 - 3*x*NC2*lx2 +  
      3*x*z*NC2*lx2 - 3*omxi*lx2 +  
      (3*z*omxi*lx2)/2. +  
      3*NC2*omxi*lx2 -  
      (3*z*NC2*omxi*lx2)/2. -  
      (3*zi*lx2)/2. - (3*x*zi*lx2)/2. +  
      (3*NC2*zi*lx2)/2. +  
      (3*x*NC2*zi*lx2)/2. +  
      3*omxi*zi*lx2 -  
      3*NC2*omxi*zi*lx2 +  
      lomz2/2. + (x*lomz2)/2. -  
      (x*z*lomz2)/2. - (NC2*lomz2)/2. -  
      (x*NC2*lomz2)/2. +  
      (x*z*NC2*lomz2)/2. -  
      (omxi*lomz2)/2. +  
      (z*omxi*lomz2)/4. +  
      (NC2*omxi*lomz2)/2. -  
      (z*NC2*omxi*lomz2)/4. -  
      (zi*lomz2)/4. - (x*zi*lomz2)/4. +  
      (NC2*zi*lomz2)/4. +  
      (x*NC2*zi*lomz2)/4. +  
      (omxi*zi*lomz2)/2. -  
      (NC2*omxi*zi*lomz2)/2. +  
      2*lz2 + 2*x*lz2 - 2*x*z*lz2 -  
      2*NC2*lz2 - 2*x*NC2*lz2 +  
      2*x*z*NC2*lz2 - 2*omxi*lz2 +  
      z*omxi*lz2 +  
      2*NC2*omxi*lz2 -  
      z*NC2*omxi*lz2 - zi*lz2 -  
      x*zi*lz2 + NC2*zi*lz2 +  
      x*NC2*zi*lz2 +  
      2*omxi*zi*lz2 -  
      2*NC2*omxi*zi*lz2)*Tt2 +  
   (-3 + 2*x*Li2z - x*z*Li2z + li2spec13 -  
      x*li2spec13 - li2spec14 -  
      x*li2spec14 +  
      x*z*li2spec14 -  
      li2spec11 - x*li2spec11 +  
      x*z*li2spec11 +  
      2*x*li2spec15 -  
      x*z*li2spec15 + 3*x*lomx -  
      3*x*z*lomx - 6*x*lx + 6*x*z*lx - 4*lomx*lx -  
      16*x*lomx*lx + 10*x*z*lomx*lx + 4*x*lomz -  
      4*x*z*lomz + 3*lomx*lomz +  
      13*x*lomx*lomz - 8*x*z*lomx*lomz -  
      5*lx*lomz - 17*x*lx*lomz +  
      11*x*z*lx*lomz + lx*lomxmz +  
      3*x*lx*lomxmz - 2*x*z*lx*lomxmz -  
      lomz*lomxmz - 3*x*lomz*lomxmz +  
      2*x*z*lomz*lomxmz + lx*lxmz +  
      x*lx*lxmz - x*z*lx*lxmz + 3*x*lz -  
      3*x*z*lz + 2*lomx*lz + 6*x*lomx*lz -  
      4*x*z*lomx*lz - 5*lx*lz - 11*x*lx*lz +  
      8*x*z*lx*lz + 4*lomz*lz + 8*x*lomz*lz -  
      6*x*z*lomz*lz - lxmz*lz - x*lxmz*lz +  
      x*z*lxmz*lz + 3*NCi2 - 2*x*Li2z*NCi2 +  
      x*z*Li2z*NCi2 - li2spec13*NCi2 +  
      x*li2spec13*NCi2 +  
      li2spec14*NCi2 +  
      x*li2spec14*NCi2 -  
      x*z*li2spec14*NCi2 +  
      li2spec11*NCi2 +  
      x*li2spec11*NCi2 -  
      x*z*li2spec11*NCi2 -  
      2*x*li2spec15*NCi2 +  
      x*z*li2spec15*NCi2 -  
      3*x*lomx*NCi2 + 3*x*z*lomx*NCi2 +  
      6*x*lx*NCi2 - 6*x*z*lx*NCi2 +  
      4*lomx*lx*NCi2 + 16*x*lomx*lx*NCi2 -  
      10*x*z*lomx*lx*NCi2 - 4*x*lomz*NCi2 +  
      4*x*z*lomz*NCi2 - 3*lomx*lomz*NCi2 -  
      13*x*lomx*lomz*NCi2 +  
      8*x*z*lomx*lomz*NCi2 +  
      5*lx*lomz*NCi2 + 17*x*lx*lomz*NCi2 -  
      11*x*z*lx*lomz*NCi2 -  
      lx*lomxmz*NCi2 -  
      3*x*lx*lomxmz*NCi2 +  
      2*x*z*lx*lomxmz*NCi2 +  
      lomz*lomxmz*NCi2 +  
      3*x*lomz*lomxmz*NCi2 -  
      2*x*z*lomz*lomxmz*NCi2 -  
      lx*lxmz*NCi2 - x*lx*lxmz*NCi2 +  
      x*z*lx*lxmz*NCi2 - 3*x*lz*NCi2 +  
      3*x*z*lz*NCi2 - 2*lomx*lz*NCi2 -  
      6*x*lomx*lz*NCi2 + 4*x*z*lomx*lz*NCi2 +  
      5*lx*lz*NCi2 + 11*x*lx*lz*NCi2 -  
      8*x*z*lx*lz*NCi2 - 4*lomz*lz*NCi2 -  
      8*x*lomz*lz*NCi2 + 6*x*z*lomz*lz*NCi2 +  
      lxmz*lz*NCi2 + x*lxmz*lz*NCi2 -  
      x*z*lxmz*lz*NCi2 - pi2/6. - (3*x*pi2)/2. +  
      (5*x*z*pi2)/6. + (NCi2*pi2)/6. +  
      (3*x*NCi2*pi2)/2. - (5*x*z*NCi2*pi2)/6. +  
      2*z*omxi - Li2z*omxi +  
      (z*Li2z*omxi)/2. - li2spec13*omxi +  
      (z*li2spec13*omxi)/2. +  
      2*li2spec14*omxi -  
      z*li2spec14*omxi +  
      2*li2spec11*omxi -  
      z*li2spec11*omxi -  
      li2spec15*omxi +  
      (z*li2spec15*omxi)/2. -  
      3*lomx*omxi + 3*z*lomx*omxi +  
      6*lx*omxi - 6*z*lx*omxi +  
      14*lomx*lx*omxi -  
      7*z*lomx*lx*omxi - 4*lomz*omxi +  
      4*z*lomz*omxi - 11*lomx*lomz*omxi +  
      (11*z*lomx*lomz*omxi)/2. +  
      16*lx*lomz*omxi -  
      8*z*lx*lomz*omxi -  
      3*lx*lomxmz*omxi +  
      (3*z*lx*lomxmz*omxi)/2. +  
      3*lomz*lomxmz*omxi -  
      (3*z*lomz*lomxmz*omxi)/2. -  
      2*lx*lxmz*omxi + z*lx*lxmz*omxi -  
      3*lz*omxi + 3*z*lz*omxi -  
      6*lomx*lz*omxi +  
      3*z*lomx*lz*omxi + 13*lx*lz*omxi -  
      (13*z*lx*lz*omxi)/2. -  
      10*lomz*lz*omxi +  
      5*z*lomz*lz*omxi +  
      2*lxmz*lz*omxi - z*lxmz*lz*omxi -  
      2*z*NCi2*omxi + Li2z*NCi2*omxi -  
      (z*Li2z*NCi2*omxi)/2. +  
      li2spec13*NCi2*omxi -  
      (z*li2spec13*NCi2*omxi)/2. -  
      2*li2spec14*NCi2*omxi +  
      z*li2spec14*NCi2*omxi -  
      2*li2spec11*NCi2*omxi +  
      z*li2spec11*NCi2*omxi +  
      li2spec15*NCi2*omxi -  
      (z*li2spec15*NCi2*omxi)/2. +  
      3*lomx*NCi2*omxi -  
      3*z*lomx*NCi2*omxi -  
      6*lx*NCi2*omxi +  
      6*z*lx*NCi2*omxi -  
      14*lomx*lx*NCi2*omxi +  
      7*z*lomx*lx*NCi2*omxi +  
      4*lomz*NCi2*omxi -  
      4*z*lomz*NCi2*omxi +  
      11*lomx*lomz*NCi2*omxi -  
      (11*z*lomx*lomz*NCi2*omxi)/2. -  
      16*lx*lomz*NCi2*omxi +  
      8*z*lx*lomz*NCi2*omxi +  
      3*lx*lomxmz*NCi2*omxi -  
      (3*z*lx*lomxmz*NCi2*omxi)/2. -  
      3*lomz*lomxmz*NCi2*omxi +  
      (3*z*lomz*lomxmz*NCi2*omxi)/2. +  
      2*lx*lxmz*NCi2*omxi -  
      z*lx*lxmz*NCi2*omxi +  
      3*lz*NCi2*omxi -  
      3*z*lz*NCi2*omxi +  
      6*lomx*lz*NCi2*omxi -  
      3*z*lomx*lz*NCi2*omxi -  
      13*lx*lz*NCi2*omxi +  
      (13*z*lx*lz*NCi2*omxi)/2. +  
      10*lomz*lz*NCi2*omxi -  
      5*z*lomz*lz*NCi2*omxi -  
      2*lxmz*lz*NCi2*omxi +  
      z*lxmz*lz*NCi2*omxi +  
      pi2*omxi - (z*pi2*omxi)/2. -  
      NCi2*pi2*omxi +  
      (z*NCi2*pi2*omxi)/2. + 2*zi -  
      2*x*zi - x*Li2z*zi -  
      (li2spec13*zi)/2. +  
      (x*li2spec13*zi)/2. +  
      (li2spec14*zi)/2. +  
      (x*li2spec14*zi)/2. +  
      (li2spec11*zi)/2. +  
      (x*li2spec11*zi)/2. -  
      x*li2spec15*zi -  
      3*x*lomx*zi + 6*x*lx*zi +  
      2*lomx*lx*zi + 8*x*lomx*lx*zi -  
      4*x*lomz*zi - (3*lomx*lomz*zi)/2. -  
      (13*x*lomx*lomz*zi)/2. +  
      (5*lx*lomz*zi)/2. +  
      (17*x*lx*lomz*zi)/2. -  
      (lx*lomxmz*zi)/2. -  
      (3*x*lx*lomxmz*zi)/2. +  
      (lomz*lomxmz*zi)/2. +  
      (3*x*lomz*lomxmz*zi)/2. -  
      (lx*lxmz*zi)/2. - (x*lx*lxmz*zi)/2. -  
      3*x*lz*zi - lomx*lz*zi -  
      3*x*lomx*lz*zi + (5*lx*lz*zi)/2. +  
      (11*x*lx*lz*zi)/2. - 2*lomz*lz*zi -  
      4*x*lomz*lz*zi + (lxmz*lz*zi)/2. +  
      (x*lxmz*lz*zi)/2. - 2*NCi2*zi +  
      2*x*NCi2*zi + x*Li2z*NCi2*zi +  
      (li2spec13*NCi2*zi)/2. -  
      (x*li2spec13*NCi2*zi)/2. -  
      (li2spec14*NCi2*zi)/2. -  
      (x*li2spec14*NCi2*zi)/2. -  
      (li2spec11*NCi2*zi)/2. -  
      (x*li2spec11*NCi2*zi)/2. +  
      x*li2spec15*NCi2*zi +  
      3*x*lomx*NCi2*zi - 6*x*lx*NCi2*zi -  
      2*lomx*lx*NCi2*zi -  
      8*x*lomx*lx*NCi2*zi +  
      4*x*lomz*NCi2*zi +  
      (3*lomx*lomz*NCi2*zi)/2. +  
      (13*x*lomx*lomz*NCi2*zi)/2. -  
      (5*lx*lomz*NCi2*zi)/2. -  
      (17*x*lx*lomz*NCi2*zi)/2. +  
      (lx*lomxmz*NCi2*zi)/2. +  
      (3*x*lx*lomxmz*NCi2*zi)/2. -  
      (lomz*lomxmz*NCi2*zi)/2. -  
      (3*x*lomz*lomxmz*NCi2*zi)/2. +  
      (lx*lxmz*NCi2*zi)/2. +  
      (x*lx*lxmz*NCi2*zi)/2. +  
      3*x*lz*NCi2*zi +  
      lomx*lz*NCi2*zi +  
      3*x*lomx*lz*NCi2*zi -  
      (5*lx*lz*NCi2*zi)/2. -  
      (11*x*lx*lz*NCi2*zi)/2. +  
      2*lomz*lz*NCi2*zi +  
      4*x*lomz*lz*NCi2*zi -  
      (lxmz*lz*NCi2*zi)/2. -  
      (x*lxmz*lz*NCi2*zi)/2. +  
      (pi2*zi)/12. + (3*x*pi2*zi)/4. -  
      (NCi2*pi2*zi)/12. -  
      (3*x*NCi2*pi2*zi)/4. +  
      Li2z*omxi*zi -  
      li2spec14*omxi*zi -  
      li2spec11*omxi*zi +  
      li2spec15*omxi*zi +  
      3*lomx*omxi*zi -  
      6*lx*omxi*zi -  
      10*lomx*lx*omxi*zi +  
      4*lomz*omxi*zi +  
      8*lomx*lomz*omxi*zi -  
      11*lx*lomz*omxi*zi +  
      2*lx*lomxmz*omxi*zi -  
      2*lomz*lomxmz*omxi*zi +  
      lx*lxmz*omxi*zi +  
      3*lz*omxi*zi +  
      4*lomx*lz*omxi*zi -  
      8*lx*lz*omxi*zi +  
      6*lomz*lz*omxi*zi -  
      lxmz*lz*omxi*zi -  
      Li2z*NCi2*omxi*zi +  
      li2spec14*NCi2*omxi* 
       zi + li2spec11*NCi2*omxi* 
       zi - li2spec15*NCi2* 
       omxi*zi -  
      3*lomx*NCi2*omxi*zi +  
      6*lx*NCi2*omxi*zi +  
      10*lomx*lx*NCi2*omxi*zi -  
      4*lomz*NCi2*omxi*zi -  
      8*lomx*lomz*NCi2*omxi*zi +  
      11*lx*lomz*NCi2*omxi*zi -  
      2*lx*lomxmz*NCi2*omxi*zi +  
      2*lomz*lomxmz*NCi2*omxi*zi -  
      lx*lxmz*NCi2*omxi*zi -  
      3*lz*NCi2*omxi*zi -  
      4*lomx*lz*NCi2*omxi*zi +  
      8*lx*lz*NCi2*omxi*zi -  
      6*lomz*lz*NCi2*omxi*zi +  
      lxmz*lz*NCi2*omxi*zi -  
      (5*pi2*omxi*zi)/6. +  
      (5*NCi2*pi2*omxi*zi)/6. +  
      lomx2/2. + (9*x*lomx2)/2. -  
      (5*x*z*lomx2)/2. - (NCi2*lomx2)/2. -  
      (9*x*NCi2*lomx2)/2. +  
      (5*x*z*NCi2*lomx2)/2. -  
      3*omxi*lomx2 +  
      (3*z*omxi*lomx2)/2. +  
      3*NCi2*omxi*lomx2 -  
      (3*z*NCi2*omxi*lomx2)/2. -  
      (zi*lomx2)/4. -  
      (9*x*zi*lomx2)/4. +  
      (NCi2*zi*lomx2)/4. +  
      (9*x*NCi2*zi*lomx2)/4. +  
      (5*omxi*zi*lomx2)/2. -  
      (5*NCi2*omxi*zi*lomx2)/2. +  
      3*lx2 + 10*x*lx2 - (13*x*z*lx2)/2. -  
      3*NCi2*lx2 - 10*x*NCi2*lx2 +  
      (13*x*z*NCi2*lx2)/2. -  
      (19*omxi*lx2)/2. +  
      (19*z*omxi*lx2)/4. +  
      (19*NCi2*omxi*lx2)/2. -  
      (19*z*NCi2*omxi*lx2)/4. -  
      (3*zi*lx2)/2. - 5*x*zi*lx2 +  
      (3*NCi2*zi*lx2)/2. +  
      5*x*NCi2*zi*lx2 +  
      (13*omxi*zi*lx2)/2. -  
      (13*NCi2*omxi*zi*lx2)/2. +  
      (3*lomz2)/2. + (13*x*lomz2)/2. -  
      4*x*z*lomz2 - (3*NCi2*lomz2)/2. -  
      (13*x*NCi2*lomz2)/2. +  
      4*x*z*NCi2*lomz2 -  
      (11*omxi*lomz2)/2. +  
      (11*z*omxi*lomz2)/4. +  
      (11*NCi2*omxi*lomz2)/2. -  
      (11*z*NCi2*omxi*lomz2)/4. -  
      (3*zi*lomz2)/4. -  
      (13*x*zi*lomz2)/4. +  
      (3*NCi2*zi*lomz2)/4. +  
      (13*x*NCi2*zi*lomz2)/4. +  
      4*omxi*zi*lomz2 -  
      4*NCi2*omxi*zi*lomz2 +  
      2*lz2 + 3*x*lz2 - (5*x*z*lz2)/2. -  
      2*NCi2*lz2 - 3*x*NCi2*lz2 +  
      (5*x*z*NCi2*lz2)/2. -  
      (9*omxi*lz2)/2. +  
      (9*z*omxi*lz2)/4. +  
      (9*NCi2*omxi*lz2)/2. -  
      (9*z*NCi2*omxi*lz2)/4. -  
      zi*lz2 - (3*x*zi*lz2)/2. +  
      NCi2*zi*lz2 +  
      (3*x*NCi2*zi*lz2)/2. +  
      (5*omxi*zi*lz2)/2. -  
      (5*NCi2*omxi*zi*lz2)/2.)*Tu1 +  
   (-3 + 2*x*Li2z - x*z*Li2z - li2spec11 -  
      x*li2spec11 + x*z*li2spec11 -  
      li2spec16 + x*li2spec16 -  
      2*x*li2spec17 +  
      x*z*li2spec17 +  
      li2spec18 +  
      x*li2spec18 -  
      x*z*li2spec18 + 3*x*lomx -  
      3*x*z*lomx - 6*x*lx + 6*x*z*lx - 3*lomx*lx -  
      13*x*lomx*lx + 8*x*z*lomx*lx + 4*x*lomz -  
      4*x*z*lomz + 2*lomx*lomz +  
      10*x*lomx*lomz - 6*x*z*lomx*lomz -  
      6*lx*lomz - 16*x*lx*lomz +  
      11*x*z*lx*lomz + lx*lxmz + x*lx*lxmz -  
      x*z*lx*lxmz + 3*x*lz - 3*x*z*lz +  
      2*lomx*lz + 6*x*lomx*lz - 4*x*z*lomx*lz -  
      6*lx*lz - 14*x*lx*lz + 10*x*z*lx*lz +  
      5*lomz*lz + 11*x*lomz*lz -  
      8*x*z*lomz*lz - lxmz*lz - x*lxmz*lz +  
      x*z*lxmz*lz + lx*lmopxpz +  
      3*x*lx*lmopxpz - 2*x*z*lx*lmopxpz -  
      lomz*lmopxpz - 3*x*lomz*lmopxpz +  
      2*x*z*lomz*lmopxpz + 3*NCi2 -  
      2*x*Li2z*NCi2 + x*z*Li2z*NCi2 +  
      li2spec11*NCi2 +  
      x*li2spec11*NCi2 -  
      x*z*li2spec11*NCi2 +  
      li2spec16*NCi2 -  
      x*li2spec16*NCi2 +  
      2*x*li2spec17*NCi2 -  
      x*z*li2spec17*NCi2 -  
      li2spec18*NCi2 -  
      x*li2spec18*NCi2 +  
      x*z*li2spec18*NCi2 -  
      3*x*lomx*NCi2 + 3*x*z*lomx*NCi2 +  
      6*x*lx*NCi2 - 6*x*z*lx*NCi2 +  
      3*lomx*lx*NCi2 + 13*x*lomx*lx*NCi2 -  
      8*x*z*lomx*lx*NCi2 - 4*x*lomz*NCi2 +  
      4*x*z*lomz*NCi2 - 2*lomx*lomz*NCi2 -  
      10*x*lomx*lomz*NCi2 +  
      6*x*z*lomx*lomz*NCi2 +  
      6*lx*lomz*NCi2 + 16*x*lx*lomz*NCi2 -  
      11*x*z*lx*lomz*NCi2 - lx*lxmz*NCi2 -  
      x*lx*lxmz*NCi2 + x*z*lx*lxmz*NCi2 -  
      3*x*lz*NCi2 + 3*x*z*lz*NCi2 -  
      2*lomx*lz*NCi2 - 6*x*lomx*lz*NCi2 +  
      4*x*z*lomx*lz*NCi2 + 6*lx*lz*NCi2 +  
      14*x*lx*lz*NCi2 - 10*x*z*lx*lz*NCi2 -  
      5*lomz*lz*NCi2 - 11*x*lomz*lz*NCi2 +  
      8*x*z*lomz*lz*NCi2 + lxmz*lz*NCi2 +  
      x*lxmz*lz*NCi2 - x*z*lxmz*lz*NCi2 -  
      lx*lmopxpz*NCi2 -  
      3*x*lx*lmopxpz*NCi2 +  
      2*x*z*lx*lmopxpz*NCi2 +  
      lomz*lmopxpz*NCi2 +  
      3*x*lomz*lmopxpz*NCi2 -  
      2*x*z*lomz*lmopxpz*NCi2 - pi2/6. -  
      (3*x*pi2)/2. + (5*x*z*pi2)/6. + (NCi2*pi2)/6. +  
      (3*x*NCi2*pi2)/2. - (5*x*z*NCi2*pi2)/6. +  
      2*z*omxi - Li2z*omxi +  
      (z*Li2z*omxi)/2. +  
      2*li2spec11*omxi -  
      z*li2spec11*omxi +  
      li2spec16*omxi -  
      (z*li2spec16*omxi)/2. +  
      li2spec17*omxi -  
      (z*li2spec17*omxi)/2. -  
      2*li2spec18*omxi +  
      z*li2spec18*omxi -  
      3*lomx*omxi + 3*z*lomx*omxi +  
      6*lx*omxi - 6*z*lx*omxi +  
      11*lomx*lx*omxi -  
      (11*z*lomx*lx*omxi)/2. -  
      4*lomz*omxi + 4*z*lomz*omxi -  
      8*lomx*lomz*omxi +  
      4*z*lomx*lomz*omxi +  
      17*lx*lomz*omxi -  
      (17*z*lx*lomz*omxi)/2. -  
      2*lx*lxmz*omxi + z*lx*lxmz*omxi -  
      3*lz*omxi + 3*z*lz*omxi -  
      6*lomx*lz*omxi +  
      3*z*lomx*lz*omxi + 16*lx*lz*omxi -  
      8*z*lx*lz*omxi - 13*lomz*lz*omxi +  
      (13*z*lomz*lz*omxi)/2. +  
      2*lxmz*lz*omxi - z*lxmz*lz*omxi -  
      3*lx*lmopxpz*omxi +  
      (3*z*lx*lmopxpz*omxi)/2. +  
      3*lomz*lmopxpz*omxi -  
      (3*z*lomz*lmopxpz*omxi)/2. -  
      2*z*NCi2*omxi + Li2z*NCi2*omxi -  
      (z*Li2z*NCi2*omxi)/2. -  
      2*li2spec11*NCi2*omxi +  
      z*li2spec11*NCi2*omxi -  
      li2spec16*NCi2*omxi +  
      (z*li2spec16*NCi2*omxi)/2. -  
      li2spec17*NCi2*omxi +  
      (z*li2spec17*NCi2*omxi)/ 
       2. + 2*li2spec18*NCi2* 
       omxi - z*li2spec18*NCi2* 
       omxi + 3*lomx*NCi2*omxi -  
      3*z*lomx*NCi2*omxi -  
      6*lx*NCi2*omxi +  
      6*z*lx*NCi2*omxi -  
      11*lomx*lx*NCi2*omxi +  
      (11*z*lomx*lx*NCi2*omxi)/2. +  
      4*lomz*NCi2*omxi -  
      4*z*lomz*NCi2*omxi +  
      8*lomx*lomz*NCi2*omxi -  
      4*z*lomx*lomz*NCi2*omxi -  
      17*lx*lomz*NCi2*omxi +  
      (17*z*lx*lomz*NCi2*omxi)/2. +  
      2*lx*lxmz*NCi2*omxi -  
      z*lx*lxmz*NCi2*omxi +  
      3*lz*NCi2*omxi -  
      3*z*lz*NCi2*omxi +  
      6*lomx*lz*NCi2*omxi -  
      3*z*lomx*lz*NCi2*omxi -  
      16*lx*lz*NCi2*omxi +  
      8*z*lx*lz*NCi2*omxi +  
      13*lomz*lz*NCi2*omxi -  
      (13*z*lomz*lz*NCi2*omxi)/2. -  
      2*lxmz*lz*NCi2*omxi +  
      z*lxmz*lz*NCi2*omxi +  
      3*lx*lmopxpz*NCi2*omxi -  
      (3*z*lx*lmopxpz*NCi2*omxi)/2. -  
      3*lomz*lmopxpz*NCi2*omxi +  
      (3*z*lomz*lmopxpz*NCi2*omxi)/2. +  
      pi2*omxi - (z*pi2*omxi)/2. -  
      NCi2*pi2*omxi +  
      (z*NCi2*pi2*omxi)/2. + 2*zi -  
      2*x*zi - x*Li2z*zi +  
      (li2spec11*zi)/2. +  
      (x*li2spec11*zi)/2. +  
      (li2spec16*zi)/2. -  
      (x*li2spec16*zi)/2. +  
      x*li2spec17*zi -  
      (li2spec18*zi)/2. -  
      (x*li2spec18*zi)/2. -  
      3*x*lomx*zi + 6*x*lx*zi +  
      (3*lomx*lx*zi)/2. +  
      (13*x*lomx*lx*zi)/2. - 4*x*lomz*zi -  
      lomx*lomz*zi - 5*x*lomx*lomz*zi +  
      3*lx*lomz*zi + 8*x*lx*lomz*zi -  
      (lx*lxmz*zi)/2. - (x*lx*lxmz*zi)/2. -  
      3*x*lz*zi - lomx*lz*zi -  
      3*x*lomx*lz*zi + 3*lx*lz*zi +  
      7*x*lx*lz*zi - (5*lomz*lz*zi)/2. -  
      (11*x*lomz*lz*zi)/2. +  
      (lxmz*lz*zi)/2. + (x*lxmz*lz*zi)/2. -  
      (lx*lmopxpz*zi)/2. -  
      (3*x*lx*lmopxpz*zi)/2. +  
      (lomz*lmopxpz*zi)/2. +  
      (3*x*lomz*lmopxpz*zi)/2. -  
      2*NCi2*zi + 2*x*NCi2*zi +  
      x*Li2z*NCi2*zi -  
      (li2spec11*NCi2*zi)/2. -  
      (x*li2spec11*NCi2*zi)/2. -  
      (li2spec16*NCi2*zi)/2. +  
      (x*li2spec16*NCi2*zi)/2. -  
      x*li2spec17*NCi2*zi +  
      (li2spec18*NCi2*zi)/2. +  
      (x*li2spec18*NCi2*zi)/2. +  
      3*x*lomx*NCi2*zi - 6*x*lx*NCi2*zi -  
      (3*lomx*lx*NCi2*zi)/2. -  
      (13*x*lomx*lx*NCi2*zi)/2. +  
      4*x*lomz*NCi2*zi +  
      lomx*lomz*NCi2*zi +  
      5*x*lomx*lomz*NCi2*zi -  
      3*lx*lomz*NCi2*zi -  
      8*x*lx*lomz*NCi2*zi +  
      (lx*lxmz*NCi2*zi)/2. +  
      (x*lx*lxmz*NCi2*zi)/2. +  
      3*x*lz*NCi2*zi +  
      lomx*lz*NCi2*zi +  
      3*x*lomx*lz*NCi2*zi -  
      3*lx*lz*NCi2*zi -  
      7*x*lx*lz*NCi2*zi +  
      (5*lomz*lz*NCi2*zi)/2. +  
      (11*x*lomz*lz*NCi2*zi)/2. -  
      (lxmz*lz*NCi2*zi)/2. -  
      (x*lxmz*lz*NCi2*zi)/2. +  
      (lx*lmopxpz*NCi2*zi)/2. +  
      (3*x*lx*lmopxpz*NCi2*zi)/2. -  
      (lomz*lmopxpz*NCi2*zi)/2. -  
      (3*x*lomz*lmopxpz*NCi2*zi)/2. +  
      (pi2*zi)/12. + (3*x*pi2*zi)/4. -  
      (NCi2*pi2*zi)/12. -  
      (3*x*NCi2*pi2*zi)/4. +  
      Li2z*omxi*zi -  
      li2spec11*omxi*zi -  
      li2spec17*omxi*zi +  
      li2spec18*omxi*zi +  
      3*lomx*omxi*zi -  
      6*lx*omxi*zi -  
      8*lomx*lx*omxi*zi +  
      4*lomz*omxi*zi +  
      6*lomx*lomz*omxi*zi -  
      11*lx*lomz*omxi*zi +  
      lx*lxmz*omxi*zi +  
      3*lz*omxi*zi +  
      4*lomx*lz*omxi*zi -  
      10*lx*lz*omxi*zi +  
      8*lomz*lz*omxi*zi -  
      lxmz*lz*omxi*zi +  
      2*lx*lmopxpz*omxi*zi -  
      2*lomz*lmopxpz*omxi*zi -  
      Li2z*NCi2*omxi*zi +  
      li2spec11*NCi2*omxi*zi +  
      li2spec17*NCi2*omxi* 
       zi - li2spec18*NCi2* 
       omxi*zi -  
      3*lomx*NCi2*omxi*zi +  
      6*lx*NCi2*omxi*zi +  
      8*lomx*lx*NCi2*omxi*zi -  
      4*lomz*NCi2*omxi*zi -  
      6*lomx*lomz*NCi2*omxi*zi +  
      11*lx*lomz*NCi2*omxi*zi -  
      lx*lxmz*NCi2*omxi*zi -  
      3*lz*NCi2*omxi*zi -  
      4*lomx*lz*NCi2*omxi*zi +  
      10*lx*lz*NCi2*omxi*zi -  
      8*lomz*lz*NCi2*omxi*zi +  
      lxmz*lz*NCi2*omxi*zi -  
      2*lx*lmopxpz*NCi2*omxi*zi +  
      2*lomz*lmopxpz*NCi2*omxi*zi -  
      (5*pi2*omxi*zi)/6. +  
      (5*NCi2*pi2*omxi*zi)/6. +  
      lomx2/2. + (9*x*lomx2)/2. -  
      (5*x*z*lomx2)/2. - (NCi2*lomx2)/2. -  
      (9*x*NCi2*lomx2)/2. +  
      (5*x*z*NCi2*lomx2)/2. -  
      3*omxi*lomx2 +  
      (3*z*omxi*lomx2)/2. +  
      3*NCi2*omxi*lomx2 -  
      (3*z*NCi2*omxi*lomx2)/2. -  
      (zi*lomx2)/4. -  
      (9*x*zi*lomx2)/4. +  
      (NCi2*zi*lomx2)/4. +  
      (9*x*NCi2*zi*lomx2)/4. +  
      (5*omxi*zi*lomx2)/2. -  
      (5*NCi2*omxi*zi*lomx2)/2. +  
      (7*lx2)/2. + (19*x*lx2)/2. -  
      (13*x*z*lx2)/2. - (7*NCi2*lx2)/2. -  
      (19*x*NCi2*lx2)/2. +  
      (13*x*z*NCi2*lx2)/2. - 10*omxi*lx2 +  
      5*z*omxi*lx2 +  
      10*NCi2*omxi*lx2 -  
      5*z*NCi2*omxi*lx2 -  
      (7*zi*lx2)/4. - (19*x*zi*lx2)/4. +  
      (7*NCi2*zi*lx2)/4. +  
      (19*x*NCi2*zi*lx2)/4. +  
      (13*omxi*zi*lx2)/2. -  
      (13*NCi2*omxi*zi*lx2)/2. +  
      2*lomz2 + 6*x*lomz2 - 4*x*z*lomz2 -  
      2*NCi2*lomz2 - 6*x*NCi2*lomz2 +  
      4*x*z*NCi2*lomz2 -  
      6*omxi*lomz2 +  
      3*z*omxi*lomz2 +  
      6*NCi2*omxi*lomz2 -  
      3*z*NCi2*omxi*lomz2 -  
      zi*lomz2 - 3*x*zi*lomz2 +  
      NCi2*zi*lomz2 +  
      3*x*NCi2*zi*lomz2 +  
      4*omxi*zi*lomz2 -  
      4*NCi2*omxi*zi*lomz2 +  
      2*lz2 + 3*x*lz2 - (5*x*z*lz2)/2. -  
      2*NCi2*lz2 - 3*x*NCi2*lz2 +  
      (5*x*z*NCi2*lz2)/2. -  
      (9*omxi*lz2)/2. +  
      (9*z*omxi*lz2)/4. +  
      (9*NCi2*omxi*lz2)/2. -  
      (9*z*NCi2*omxi*lz2)/4. -  
      zi*lz2 - (3*x*zi*lz2)/2. +  
      NCi2*zi*lz2 +  
      (3*x*NCi2*zi*lz2)/2. +  
      (5*omxi*zi*lz2)/2. -  
      (5*NCi2*omxi*zi*lz2)/2.)*Tu2 +  
   (-3 + 2*x*Li2z - x*z*Li2z + li2spec9 +  
      x*li2spec9 - x*z*li2spec9 +  
      li2spec13 - x*li2spec13 +  
      2*x*li2spec15 -  
      x*z*li2spec15 +  
      li2spec18 +  
      x*li2spec18 -  
      x*z*li2spec18 + 3*x*lomx -  
      3*x*z*lomx - 6*x*lx + 6*x*z*lx - 3*lomx*lx -  
      15*x*lomx*lx + 9*x*z*lomx*lx + 4*x*lomz -  
      4*x*z*lomz + lomx*lomz + 11*x*lomx*lomz -  
      6*x*z*lomx*lomz - 6*lx*lomz -  
      18*x*lx*lomz + 12*x*z*lx*lomz +  
      lx*lomxmz + 3*x*lx*lomxmz -  
      2*x*z*lx*lomxmz - lomz*lomxmz -  
      3*x*lomz*lomxmz + 2*x*z*lomz*lomxmz +  
      3*x*lz - 3*x*z*lz + lomx*lz + 5*x*lomx*lz -  
      3*x*z*lomx*lz - 6*lx*lz - 12*x*lx*lz +  
      9*x*z*lx*lz + 5*lomz*lz + 9*x*lomz*lz -  
      7*x*z*lomz*lz + lx*lmxpz + x*lx*lmxpz -  
      x*z*lx*lmxpz - lz*lmxpz - x*lz*lmxpz +  
      x*z*lz*lmxpz + 3*NCi2 - 2*x*Li2z*NCi2 +  
      x*z*Li2z*NCi2 - li2spec9*NCi2 -  
      x*li2spec9*NCi2 +  
      x*z*li2spec9*NCi2 -  
      li2spec13*NCi2 + x*li2spec13*NCi2 -  
      2*x*li2spec15*NCi2 +  
      x*z*li2spec15*NCi2 -  
      li2spec18*NCi2 -  
      x*li2spec18*NCi2 +  
      x*z*li2spec18*NCi2 -  
      3*x*lomx*NCi2 + 3*x*z*lomx*NCi2 +  
      6*x*lx*NCi2 - 6*x*z*lx*NCi2 +  
      3*lomx*lx*NCi2 + 15*x*lomx*lx*NCi2 -  
      9*x*z*lomx*lx*NCi2 - 4*x*lomz*NCi2 +  
      4*x*z*lomz*NCi2 - lomx*lomz*NCi2 -  
      11*x*lomx*lomz*NCi2 +  
      6*x*z*lomx*lomz*NCi2 +  
      6*lx*lomz*NCi2 + 18*x*lx*lomz*NCi2 -  
      12*x*z*lx*lomz*NCi2 -  
      lx*lomxmz*NCi2 -  
      3*x*lx*lomxmz*NCi2 +  
      2*x*z*lx*lomxmz*NCi2 +  
      lomz*lomxmz*NCi2 +  
      3*x*lomz*lomxmz*NCi2 -  
      2*x*z*lomz*lomxmz*NCi2 - 3*x*lz*NCi2 +  
      3*x*z*lz*NCi2 - lomx*lz*NCi2 -  
      5*x*lomx*lz*NCi2 + 3*x*z*lomx*lz*NCi2 +  
      6*lx*lz*NCi2 + 12*x*lx*lz*NCi2 -  
      9*x*z*lx*lz*NCi2 - 5*lomz*lz*NCi2 -  
      9*x*lomz*lz*NCi2 + 7*x*z*lomz*lz*NCi2 -  
      lx*lmxpz*NCi2 - x*lx*lmxpz*NCi2 +  
      x*z*lx*lmxpz*NCi2 + lz*lmxpz*NCi2 +  
      x*lz*lmxpz*NCi2 - x*z*lz*lmxpz*NCi2 -  
      (5*pi2)/6. - (13*x*pi2)/6. + (3*x*z*pi2)/2. +  
      (5*NCi2*pi2)/6. + (13*x*NCi2*pi2)/6. -  
      (3*x*z*NCi2*pi2)/2. + 2*z*omxi -  
      Li2z*omxi + (z*Li2z*omxi)/2. -  
      2*li2spec9*omxi +  
      z*li2spec9*omxi -  
      li2spec13*omxi +  
      (z*li2spec13*omxi)/2. -  
      li2spec15*omxi +  
      (z*li2spec15*omxi)/2. -  
      2*li2spec18*omxi +  
      z*li2spec18*omxi -  
      3*lomx*omxi + 3*z*lomx*omxi +  
      6*lx*omxi - 6*z*lx*omxi +  
      12*lomx*lx*omxi -  
      6*z*lomx*lx*omxi - 4*lomz*omxi +  
      4*z*lomz*omxi - 7*lomx*lomz*omxi +  
      (7*z*lomx*lomz*omxi)/2. +  
      18*lx*lomz*omxi -  
      9*z*lx*lomz*omxi -  
      3*lx*lomxmz*omxi +  
      (3*z*lx*lomxmz*omxi)/2. +  
      3*lomz*lomxmz*omxi -  
      (3*z*lomz*lomxmz*omxi)/2. -  
      3*lz*omxi + 3*z*lz*omxi -  
      4*lomx*lz*omxi +  
      2*z*lomx*lz*omxi + 15*lx*lz*omxi -  
      (15*z*lx*lz*omxi)/2. -  
      12*lomz*lz*omxi +  
      6*z*lomz*lz*omxi -  
      2*lx*lmxpz*omxi +  
      z*lx*lmxpz*omxi +  
      2*lz*lmxpz*omxi -  
      z*lz*lmxpz*omxi - 2*z*NCi2*omxi +  
      Li2z*NCi2*omxi -  
      (z*Li2z*NCi2*omxi)/2. +  
      2*li2spec9*NCi2*omxi -  
      z*li2spec9*NCi2*omxi +  
      li2spec13*NCi2*omxi -  
      (z*li2spec13*NCi2*omxi)/2. +  
      li2spec15*NCi2*omxi -  
      (z*li2spec15*NCi2*omxi)/2. +  
      2*li2spec18*NCi2*omxi -  
      z*li2spec18*NCi2*omxi +  
      3*lomx*NCi2*omxi -  
      3*z*lomx*NCi2*omxi -  
      6*lx*NCi2*omxi +  
      6*z*lx*NCi2*omxi -  
      12*lomx*lx*NCi2*omxi +  
      6*z*lomx*lx*NCi2*omxi +  
      4*lomz*NCi2*omxi -  
      4*z*lomz*NCi2*omxi +  
      7*lomx*lomz*NCi2*omxi -  
      (7*z*lomx*lomz*NCi2*omxi)/2. -  
      18*lx*lomz*NCi2*omxi +  
      9*z*lx*lomz*NCi2*omxi +  
      3*lx*lomxmz*NCi2*omxi -  
      (3*z*lx*lomxmz*NCi2*omxi)/2. -  
      3*lomz*lomxmz*NCi2*omxi +  
      (3*z*lomz*lomxmz*NCi2*omxi)/2. +  
      3*lz*NCi2*omxi -  
      3*z*lz*NCi2*omxi +  
      4*lomx*lz*NCi2*omxi -  
      2*z*lomx*lz*NCi2*omxi -  
      15*lx*lz*NCi2*omxi +  
      (15*z*lx*lz*NCi2*omxi)/2. +  
      12*lomz*lz*NCi2*omxi -  
      6*z*lomz*lz*NCi2*omxi +  
      2*lx*lmxpz*NCi2*omxi -  
      z*lx*lmxpz*NCi2*omxi -  
      2*lz*lmxpz*NCi2*omxi +  
      z*lz*lmxpz*NCi2*omxi +  
      (7*pi2*omxi)/3. - (7*z*pi2*omxi)/6. -  
      (7*NCi2*pi2*omxi)/3. +  
      (7*z*NCi2*pi2*omxi)/6. + 2*zi -  
      2*x*zi - x*Li2z*zi -  
      (li2spec9*zi)/2. -  
      (x*li2spec9*zi)/2. -  
      (li2spec13*zi)/2. +  
      (x*li2spec13*zi)/2. -  
      x*li2spec15*zi -  
      (li2spec18*zi)/2. -  
      (x*li2spec18*zi)/2. -  
      3*x*lomx*zi + 6*x*lx*zi +  
      (3*lomx*lx*zi)/2. +  
      (15*x*lomx*lx*zi)/2. - 4*x*lomz*zi -  
      (lomx*lomz*zi)/2. -  
      (11*x*lomx*lomz*zi)/2. +  
      3*lx*lomz*zi + 9*x*lx*lomz*zi -  
      (lx*lomxmz*zi)/2. -  
      (3*x*lx*lomxmz*zi)/2. +  
      (lomz*lomxmz*zi)/2. +  
      (3*x*lomz*lomxmz*zi)/2. - 3*x*lz*zi -  
      (lomx*lz*zi)/2. -  
      (5*x*lomx*lz*zi)/2. + 3*lx*lz*zi +  
      6*x*lx*lz*zi - (5*lomz*lz*zi)/2. -  
      (9*x*lomz*lz*zi)/2. -  
      (lx*lmxpz*zi)/2. -  
      (x*lx*lmxpz*zi)/2. +  
      (lz*lmxpz*zi)/2. +  
      (x*lz*lmxpz*zi)/2. - 2*NCi2*zi +  
      2*x*NCi2*zi + x*Li2z*NCi2*zi +  
      (li2spec9*NCi2*zi)/2. +  
      (x*li2spec9*NCi2*zi)/2. +  
      (li2spec13*NCi2*zi)/2. -  
      (x*li2spec13*NCi2*zi)/2. +  
      x*li2spec15*NCi2*zi +  
      (li2spec18*NCi2*zi)/2. +  
      (x*li2spec18*NCi2*zi)/2. +  
      3*x*lomx*NCi2*zi - 6*x*lx*NCi2*zi -  
      (3*lomx*lx*NCi2*zi)/2. -  
      (15*x*lomx*lx*NCi2*zi)/2. +  
      4*x*lomz*NCi2*zi +  
      (lomx*lomz*NCi2*zi)/2. +  
      (11*x*lomx*lomz*NCi2*zi)/2. -  
      3*lx*lomz*NCi2*zi -  
      9*x*lx*lomz*NCi2*zi +  
      (lx*lomxmz*NCi2*zi)/2. +  
      (3*x*lx*lomxmz*NCi2*zi)/2. -  
      (lomz*lomxmz*NCi2*zi)/2. -  
      (3*x*lomz*lomxmz*NCi2*zi)/2. +  
      3*x*lz*NCi2*zi +  
      (lomx*lz*NCi2*zi)/2. +  
      (5*x*lomx*lz*NCi2*zi)/2. -  
      3*lx*lz*NCi2*zi -  
      6*x*lx*lz*NCi2*zi +  
      (5*lomz*lz*NCi2*zi)/2. +  
      (9*x*lomz*lz*NCi2*zi)/2. +  
      (lx*lmxpz*NCi2*zi)/2. +  
      (x*lx*lmxpz*NCi2*zi)/2. -  
      (lz*lmxpz*NCi2*zi)/2. -  
      (x*lz*lmxpz*NCi2*zi)/2. +  
      (5*pi2*zi)/12. + (13*x*pi2*zi)/12. -  
      (5*NCi2*pi2*zi)/12. -  
      (13*x*NCi2*pi2*zi)/12. +  
      Li2z*omxi*zi +  
      li2spec9*omxi*zi +  
      li2spec15*omxi*zi +  
      li2spec18*omxi*zi +  
      3*lomx*omxi*zi -  
      6*lx*omxi*zi -  
      9*lomx*lx*omxi*zi +  
      4*lomz*omxi*zi +  
      6*lomx*lomz*omxi*zi -  
      12*lx*lomz*omxi*zi +  
      2*lx*lomxmz*omxi*zi -  
      2*lomz*lomxmz*omxi*zi +  
      3*lz*omxi*zi +  
      3*lomx*lz*omxi*zi -  
      9*lx*lz*omxi*zi +  
      7*lomz*lz*omxi*zi +  
      lx*lmxpz*omxi*zi -  
      lz*lmxpz*omxi*zi -  
      Li2z*NCi2*omxi*zi -  
      li2spec9*NCi2*omxi*zi -  
      li2spec15*NCi2*omxi* 
       zi - li2spec18*NCi2* 
       omxi*zi -  
      3*lomx*NCi2*omxi*zi +  
      6*lx*NCi2*omxi*zi +  
      9*lomx*lx*NCi2*omxi*zi -  
      4*lomz*NCi2*omxi*zi -  
      6*lomx*lomz*NCi2*omxi*zi +  
      12*lx*lomz*NCi2*omxi*zi -  
      2*lx*lomxmz*NCi2*omxi*zi +  
      2*lomz*lomxmz*NCi2*omxi*zi -  
      3*lz*NCi2*omxi*zi -  
      3*lomx*lz*NCi2*omxi*zi +  
      9*lx*lz*NCi2*omxi*zi -  
      7*lomz*lz*NCi2*omxi*zi -  
      lx*lmxpz*NCi2*omxi*zi +  
      lz*lmxpz*NCi2*omxi*zi -  
      (3*pi2*omxi*zi)/2. +  
      (3*NCi2*pi2*omxi*zi)/2. +  
      (3*lomx2)/2. + (11*x*lomx2)/2. -  
      (7*x*z*lomx2)/2. - (3*NCi2*lomx2)/2. -  
      (11*x*NCi2*lomx2)/2. +  
      (7*x*z*NCi2*lomx2)/2. -  
      5*omxi*lomx2 +  
      (5*z*omxi*lomx2)/2. +  
      5*NCi2*omxi*lomx2 -  
      (5*z*NCi2*omxi*lomx2)/2. -  
      (3*zi*lomx2)/4. -  
      (11*x*zi*lomx2)/4. +  
      (3*NCi2*zi*lomx2)/4. +  
      (11*x*NCi2*zi*lomx2)/4. +  
      (7*omxi*zi*lomx2)/2. -  
      (7*NCi2*omxi*zi*lomx2)/2. +  
      (7*lx2)/2. + (21*x*lx2)/2. - 7*x*z*lx2 -  
      (7*NCi2*lx2)/2. - (21*x*NCi2*lx2)/2. +  
      7*x*z*NCi2*lx2 - (21*omxi*lx2)/2. +  
      (21*z*omxi*lx2)/4. +  
      (21*NCi2*omxi*lx2)/2. -  
      (21*z*NCi2*omxi*lx2)/4. -  
      (7*zi*lx2)/4. - (21*x*zi*lx2)/4. +  
      (7*NCi2*zi*lx2)/4. +  
      (21*x*NCi2*zi*lx2)/4. +  
      7*omxi*zi*lx2 -  
      7*NCi2*omxi*zi*lx2 +  
      (5*lomz2)/2. + (15*x*lomz2)/2. -  
      5*x*z*lomz2 - (5*NCi2*lomz2)/2. -  
      (15*x*NCi2*lomz2)/2. +  
      5*x*z*NCi2*lomz2 -  
      (15*omxi*lomz2)/2. +  
      (15*z*omxi*lomz2)/4. +  
      (15*NCi2*omxi*lomz2)/2. -  
      (15*z*NCi2*omxi*lomz2)/4. -  
      (5*zi*lomz2)/4. -  
      (15*x*zi*lomz2)/4. +  
      (5*NCi2*zi*lomz2)/4. +  
      (15*x*NCi2*zi*lomz2)/4. +  
      5*omxi*zi*lomz2 -  
      5*NCi2*omxi*zi*lomz2 +  
      (5*lz2)/2. + (7*x*lz2)/2. - 3*x*z*lz2 -  
      (5*NCi2*lz2)/2. - (7*x*NCi2*lz2)/2. +  
      3*x*z*NCi2*lz2 - (11*omxi*lz2)/2. +  
      (11*z*omxi*lz2)/4. +  
      (11*NCi2*omxi*lz2)/2. -  
      (11*z*NCi2*omxi*lz2)/4. -  
      (5*zi*lz2)/4. - (7*x*zi*lz2)/4. +  
      (5*NCi2*zi*lz2)/4. +  
      (7*x*NCi2*zi*lz2)/4. +  
      3*omxi*zi*lz2 -  
      3*NCi2*omxi*zi*lz2)*Tu3 +  
   (-3 + 2*x*Li2z - x*z*Li2z + li2spec9 +  
      x*li2spec9 - x*z*li2spec9 -  
      li2spec14 -  
      x*li2spec14 +  
      x*z*li2spec14 - li2spec16 +  
      x*li2spec16 -  
      2*x*li2spec17 +  
      x*z*li2spec17 + 3*x*lomx -  
      3*x*z*lomx - 6*x*lx + 6*x*z*lx - 4*lomx*lx -  
      14*x*lomx*lx + 9*x*z*lomx*lx + 4*x*lomz -  
      4*x*z*lomz + 2*lomx*lomz +  
      10*x*lomx*lomz - 6*x*z*lomx*lomz -  
      5*lx*lomz - 15*x*lx*lomz +  
      10*x*z*lx*lomz + 3*x*lz - 3*x*z*lz +  
      3*lomx*lz + 7*x*lomx*lz - 5*x*z*lomx*lz -  
      5*lx*lz - 13*x*lx*lz + 9*x*z*lx*lz +  
      4*lomz*lz + 10*x*lomz*lz -  
      7*x*z*lomz*lz + lx*lmxpz + x*lx*lmxpz -  
      x*z*lx*lmxpz - lz*lmxpz - x*lz*lmxpz +  
      x*z*lz*lmxpz + lx*lmopxpz +  
      3*x*lx*lmopxpz - 2*x*z*lx*lmopxpz -  
      lomz*lmopxpz - 3*x*lomz*lmopxpz +  
      2*x*z*lomz*lmopxpz + 3*NCi2 -  
      2*x*Li2z*NCi2 + x*z*Li2z*NCi2 -  
      li2spec9*NCi2 -  
      x*li2spec9*NCi2 +  
      x*z*li2spec9*NCi2 +  
      li2spec14*NCi2 +  
      x*li2spec14*NCi2 -  
      x*z*li2spec14*NCi2 +  
      li2spec16*NCi2 -  
      x*li2spec16*NCi2 +  
      2*x*li2spec17*NCi2 -  
      x*z*li2spec17*NCi2 -  
      3*x*lomx*NCi2 + 3*x*z*lomx*NCi2 +  
      6*x*lx*NCi2 - 6*x*z*lx*NCi2 +  
      4*lomx*lx*NCi2 + 14*x*lomx*lx*NCi2 -  
      9*x*z*lomx*lx*NCi2 - 4*x*lomz*NCi2 +  
      4*x*z*lomz*NCi2 - 2*lomx*lomz*NCi2 -  
      10*x*lomx*lomz*NCi2 +  
      6*x*z*lomx*lomz*NCi2 +  
      5*lx*lomz*NCi2 + 15*x*lx*lomz*NCi2 -  
      10*x*z*lx*lomz*NCi2 - 3*x*lz*NCi2 +  
      3*x*z*lz*NCi2 - 3*lomx*lz*NCi2 -  
      7*x*lomx*lz*NCi2 + 5*x*z*lomx*lz*NCi2 +  
      5*lx*lz*NCi2 + 13*x*lx*lz*NCi2 -  
      9*x*z*lx*lz*NCi2 - 4*lomz*lz*NCi2 -  
      10*x*lomz*lz*NCi2 +  
      7*x*z*lomz*lz*NCi2 - lx*lmxpz*NCi2 -  
      x*lx*lmxpz*NCi2 + x*z*lx*lmxpz*NCi2 +  
      lz*lmxpz*NCi2 + x*lz*lmxpz*NCi2 -  
      x*z*lz*lmxpz*NCi2 - lx*lmopxpz*NCi2 -  
      3*x*lx*lmopxpz*NCi2 +  
      2*x*z*lx*lmopxpz*NCi2 +  
      lomz*lmopxpz*NCi2 +  
      3*x*lomz*lmopxpz*NCi2 -  
      2*x*z*lomz*lmopxpz*NCi2 - pi2/6. -  
      (3*x*pi2)/2. + (5*x*z*pi2)/6. + (NCi2*pi2)/6. +  
      (3*x*NCi2*pi2)/2. - (5*x*z*NCi2*pi2)/6. +  
      2*z*omxi - Li2z*omxi +  
      (z*Li2z*omxi)/2. -  
      2*li2spec9*omxi +  
      z*li2spec9*omxi +  
      2*li2spec14*omxi -  
      z*li2spec14*omxi +  
      li2spec16*omxi -  
      (z*li2spec16*omxi)/2. +  
      li2spec17*omxi -  
      (z*li2spec17*omxi)/2. -  
      3*lomx*omxi + 3*z*lomx*omxi +  
      6*lx*omxi - 6*z*lx*omxi +  
      13*lomx*lx*omxi -  
      (13*z*lomx*lx*omxi)/2. -  
      4*lomz*omxi + 4*z*lomz*omxi -  
      8*lomx*lomz*omxi +  
      4*z*lomx*lomz*omxi +  
      15*lx*lomz*omxi -  
      (15*z*lx*lomz*omxi)/2. - 3*lz*omxi +  
      3*z*lz*omxi - 8*lomx*lz*omxi +  
      4*z*lomx*lz*omxi + 14*lx*lz*omxi -  
      7*z*lx*lz*omxi - 11*lomz*lz*omxi +  
      (11*z*lomz*lz*omxi)/2. -  
      2*lx*lmxpz*omxi +  
      z*lx*lmxpz*omxi +  
      2*lz*lmxpz*omxi -  
      z*lz*lmxpz*omxi -  
      3*lx*lmopxpz*omxi +  
      (3*z*lx*lmopxpz*omxi)/2. +  
      3*lomz*lmopxpz*omxi -  
      (3*z*lomz*lmopxpz*omxi)/2. -  
      2*z*NCi2*omxi + Li2z*NCi2*omxi -  
      (z*Li2z*NCi2*omxi)/2. +  
      2*li2spec9*NCi2*omxi -  
      z*li2spec9*NCi2*omxi -  
      2*li2spec14*NCi2*omxi +  
      z*li2spec14*NCi2*omxi -  
      li2spec16*NCi2*omxi +  
      (z*li2spec16*NCi2*omxi)/2. -  
      li2spec17*NCi2*omxi +  
      (z*li2spec17*NCi2*omxi)/ 
       2. + 3*lomx*NCi2*omxi -  
      3*z*lomx*NCi2*omxi -  
      6*lx*NCi2*omxi +  
      6*z*lx*NCi2*omxi -  
      13*lomx*lx*NCi2*omxi +  
      (13*z*lomx*lx*NCi2*omxi)/2. +  
      4*lomz*NCi2*omxi -  
      4*z*lomz*NCi2*omxi +  
      8*lomx*lomz*NCi2*omxi -  
      4*z*lomx*lomz*NCi2*omxi -  
      15*lx*lomz*NCi2*omxi +  
      (15*z*lx*lomz*NCi2*omxi)/2. +  
      3*lz*NCi2*omxi -  
      3*z*lz*NCi2*omxi +  
      8*lomx*lz*NCi2*omxi -  
      4*z*lomx*lz*NCi2*omxi -  
      14*lx*lz*NCi2*omxi +  
      7*z*lx*lz*NCi2*omxi +  
      11*lomz*lz*NCi2*omxi -  
      (11*z*lomz*lz*NCi2*omxi)/2. +  
      2*lx*lmxpz*NCi2*omxi -  
      z*lx*lmxpz*NCi2*omxi -  
      2*lz*lmxpz*NCi2*omxi +  
      z*lz*lmxpz*NCi2*omxi +  
      3*lx*lmopxpz*NCi2*omxi -  
      (3*z*lx*lmopxpz*NCi2*omxi)/2. -  
      3*lomz*lmopxpz*NCi2*omxi +  
      (3*z*lomz*lmopxpz*NCi2*omxi)/2. +  
      pi2*omxi - (z*pi2*omxi)/2. -  
      NCi2*pi2*omxi +  
      (z*NCi2*pi2*omxi)/2. + 2*zi -  
      2*x*zi - x*Li2z*zi -  
      (li2spec9*zi)/2. -  
      (x*li2spec9*zi)/2. +  
      (li2spec14*zi)/2. +  
      (x*li2spec14*zi)/2. +  
      (li2spec16*zi)/2. -  
      (x*li2spec16*zi)/2. +  
      x*li2spec17*zi -  
      3*x*lomx*zi + 6*x*lx*zi +  
      2*lomx*lx*zi + 7*x*lomx*lx*zi -  
      4*x*lomz*zi - lomx*lomz*zi -  
      5*x*lomx*lomz*zi +  
      (5*lx*lomz*zi)/2. +  
      (15*x*lx*lomz*zi)/2. - 3*x*lz*zi -  
      (3*lomx*lz*zi)/2. -  
      (7*x*lomx*lz*zi)/2. + (5*lx*lz*zi)/2. +  
      (13*x*lx*lz*zi)/2. - 2*lomz*lz*zi -  
      5*x*lomz*lz*zi - (lx*lmxpz*zi)/2. -  
      (x*lx*lmxpz*zi)/2. +  
      (lz*lmxpz*zi)/2. +  
      (x*lz*lmxpz*zi)/2. -  
      (lx*lmopxpz*zi)/2. -  
      (3*x*lx*lmopxpz*zi)/2. +  
      (lomz*lmopxpz*zi)/2. +  
      (3*x*lomz*lmopxpz*zi)/2. -  
      2*NCi2*zi + 2*x*NCi2*zi +  
      x*Li2z*NCi2*zi +  
      (li2spec9*NCi2*zi)/2. +  
      (x*li2spec9*NCi2*zi)/2. -  
      (li2spec14*NCi2*zi)/2. -  
      (x*li2spec14*NCi2*zi)/2. -  
      (li2spec16*NCi2*zi)/2. +  
      (x*li2spec16*NCi2*zi)/2. -  
      x*li2spec17*NCi2*zi +  
      3*x*lomx*NCi2*zi - 6*x*lx*NCi2*zi -  
      2*lomx*lx*NCi2*zi -  
      7*x*lomx*lx*NCi2*zi +  
      4*x*lomz*NCi2*zi +  
      lomx*lomz*NCi2*zi +  
      5*x*lomx*lomz*NCi2*zi -  
      (5*lx*lomz*NCi2*zi)/2. -  
      (15*x*lx*lomz*NCi2*zi)/2. +  
      3*x*lz*NCi2*zi +  
      (3*lomx*lz*NCi2*zi)/2. +  
      (7*x*lomx*lz*NCi2*zi)/2. -  
      (5*lx*lz*NCi2*zi)/2. -  
      (13*x*lx*lz*NCi2*zi)/2. +  
      2*lomz*lz*NCi2*zi +  
      5*x*lomz*lz*NCi2*zi +  
      (lx*lmxpz*NCi2*zi)/2. +  
      (x*lx*lmxpz*NCi2*zi)/2. -  
      (lz*lmxpz*NCi2*zi)/2. -  
      (x*lz*lmxpz*NCi2*zi)/2. +  
      (lx*lmopxpz*NCi2*zi)/2. +  
      (3*x*lx*lmopxpz*NCi2*zi)/2. -  
      (lomz*lmopxpz*NCi2*zi)/2. -  
      (3*x*lomz*lmopxpz*NCi2*zi)/2. +  
      (pi2*zi)/12. + (3*x*pi2*zi)/4. -  
      (NCi2*pi2*zi)/12. -  
      (3*x*NCi2*pi2*zi)/4. +  
      Li2z*omxi*zi +  
      li2spec9*omxi*zi -  
      li2spec14*omxi*zi -  
      li2spec17*omxi*zi +  
      3*lomx*omxi*zi -  
      6*lx*omxi*zi -  
      9*lomx*lx*omxi*zi +  
      4*lomz*omxi*zi +  
      6*lomx*lomz*omxi*zi -  
      10*lx*lomz*omxi*zi +  
      3*lz*omxi*zi +  
      5*lomx*lz*omxi*zi -  
      9*lx*lz*omxi*zi +  
      7*lomz*lz*omxi*zi +  
      lx*lmxpz*omxi*zi -  
      lz*lmxpz*omxi*zi +  
      2*lx*lmopxpz*omxi*zi -  
      2*lomz*lmopxpz*omxi*zi -  
      Li2z*NCi2*omxi*zi -  
      li2spec9*NCi2*omxi*zi +  
      li2spec14*NCi2*omxi* 
       zi + li2spec17*NCi2* 
       omxi*zi -  
      3*lomx*NCi2*omxi*zi +  
      6*lx*NCi2*omxi*zi +  
      9*lomx*lx*NCi2*omxi*zi -  
      4*lomz*NCi2*omxi*zi -  
      6*lomx*lomz*NCi2*omxi*zi +  
      10*lx*lomz*NCi2*omxi*zi -  
      3*lz*NCi2*omxi*zi -  
      5*lomx*lz*NCi2*omxi*zi +  
      9*lx*lz*NCi2*omxi*zi -  
      7*lomz*lz*NCi2*omxi*zi -  
      lx*lmxpz*NCi2*omxi*zi +  
      lz*lmxpz*NCi2*omxi*zi -  
      2*lx*lmopxpz*NCi2*omxi*zi +  
      2*lomz*lmopxpz*NCi2*omxi*zi -  
      (5*pi2*omxi*zi)/6. +  
      (5*NCi2*pi2*omxi*zi)/6. +  
      lomx2/2. + (9*x*lomx2)/2. -  
      (5*x*z*lomx2)/2. - (NCi2*lomx2)/2. -  
      (9*x*NCi2*lomx2)/2. +  
      (5*x*z*NCi2*lomx2)/2. -  
      3*omxi*lomx2 +  
      (3*z*omxi*lomx2)/2. +  
      3*NCi2*omxi*lomx2 -  
      (3*z*NCi2*omxi*lomx2)/2. -  
      (zi*lomx2)/4. -  
      (9*x*zi*lomx2)/4. +  
      (NCi2*zi*lomx2)/4. +  
      (9*x*NCi2*zi*lomx2)/4. +  
      (5*omxi*zi*lomx2)/2. -  
      (5*NCi2*omxi*zi*lomx2)/2. +  
      3*lx2 + 9*x*lx2 - 6*x*z*lx2 -  
      3*NCi2*lx2 - 9*x*NCi2*lx2 +  
      6*x*z*NCi2*lx2 - 9*omxi*lx2 +  
      (9*z*omxi*lx2)/2. +  
      9*NCi2*omxi*lx2 -  
      (9*z*NCi2*omxi*lx2)/2. -  
      (3*zi*lx2)/2. - (9*x*zi*lx2)/2. +  
      (3*NCi2*zi*lx2)/2. +  
      (9*x*NCi2*zi*lx2)/2. +  
      6*omxi*zi*lx2 -  
      6*NCi2*omxi*zi*lx2 +  
      2*lomz2 + 6*x*lomz2 - 4*x*z*lomz2 -  
      2*NCi2*lomz2 - 6*x*NCi2*lomz2 +  
      4*x*z*NCi2*lomz2 -  
      6*omxi*lomz2 +  
      3*z*omxi*lomz2 +  
      6*NCi2*omxi*lomz2 -  
      3*z*NCi2*omxi*lomz2 -  
      zi*lomz2 - 3*x*zi*lomz2 +  
      NCi2*zi*lomz2 +  
      3*x*NCi2*zi*lomz2 +  
      4*omxi*zi*lomz2 -  
      4*NCi2*omxi*zi*lomz2 +  
      (3*lz2)/2. + (5*x*lz2)/2. - 2*x*z*lz2 -  
      (3*NCi2*lz2)/2. - (5*x*NCi2*lz2)/2. +  
      2*x*z*NCi2*lz2 - (7*omxi*lz2)/2. +  
      (7*z*omxi*lz2)/4. +  
      (7*NCi2*omxi*lz2)/2. -  
      (7*z*NCi2*omxi*lz2)/4. -  
      (3*zi*lz2)/4. - (5*x*zi*lz2)/4. +  
      (3*NCi2*zi*lz2)/4. +  
      (5*x*NCi2*zi*lz2)/4. +  
      2*omxi*zi*lz2 -  
      2*NCi2*omxi*zi*lz2)*Tu4;
};
double C2TG2G_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double sqrtxz1  = sqrt(1 - 2*z + z*z + 4*x*z);
  const double sqrtxz1i = 1. / sqrtxz1;
  const double sqrtxz3  = sqrt(x/z);
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double l2  = log(2);
  const double l22 = l2 * l2;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double x2 = x * x;
  const double x3 = x * x2;
  const double xi  = 1. / x;
  const double xi2 = xi * xi;
  const double z2 = z * z;
  const double zi  = 1. / z;
  const double zi2 = zi * zi;
  const double opxi = 1. / ( 1 + x );
  const double omzi  = 1. / ( 1 - z );
  const double opzi  = 1. / ( 1 + z );
  const double mopzi = 1. / ( - 1 + z );
  const double lx  = log(x);
  const double lx2 = lx * lx;
  const double lz  = log(z);
  const double lz2 = lz * lz;
  const double lomx  = log(1 - x);
  const double lomx2 = lomx * lomx;
  const double lomz  = log(1 - z);
  const double lomz2 = lomz * lomz;
  const double lopx  = log(1 + x);
  const double Li2x  = apfel::dilog(x);
  const double Li2mx = apfel::dilog(-x);
  const double Li2z  = apfel::dilog(z);
  const double spec1i = 1. / ( 1 + sqrtxz1 - z );
  const double lxpz    = log(x + z);
  const double lopxz   = log(1 + x*z);
  const double lopxzi  = log(1 + x*zi);
  const double lspec1   = log(1 + sqrtxz1 - z);
  const double lspec1_2 = lspec1 * lspec1;
  const double lspec2   = log(1 + sqrtxz1 + z);
  const double lspec3   = log(sqrtxz3);
  const double lspec4   = log(sqrtxz3*z);
  const double lspec7   = log(1 - 2*z + 4*x*z + z2);
  const double lspec7_2 = lspec7 * lspec7;
  const double li2spec1  = apfel::dilog(0.5 - sqrtxz1/2. - z/2.);
  const double li2spec2  = apfel::dilog(0.5 - sqrtxz1/2. + z/2.);
  const double li2spec3  = apfel::dilog(0.5 - zi/2. - (sqrtxz1*zi)/2.);
  const double li2spec4  = apfel::dilog(0.5 + zi/2. - (sqrtxz1*zi)/2.);
  const double li2spec19 = apfel::dilog(-(x*z));
  const double li2spec20 = apfel::dilog(-(x*zi));
  const double li2spec21 = apfel::dilog((1 - z)*sqrtxz1i);
  const double li2spec22 = apfel::dilog(-spec1i + sqrtxz1*spec1i + z*spec1i);
  const double li2spec23 = apfel::dilog(sqrtxz1*mopzi);
  const double atanspec1 = atan(sqrtxz3);
  const double atanspec2 = atan(sqrtxz3*z);
  const double itani1 = InvTanInt(-sqrtxz3);
  const double itani2 = InvTanInt(sqrtxz3);
  const double itani3 = InvTanInt(sqrtxz3*z);
  return -0.25*NC + (41*NC*x)/8. + 2*NC*z - (9*NC*x*z)/16. -  
   (5*NC*sqrtxz3*itani1)/8. -  
   (105*NC*sqrtxz3*x*z*itani1)/16. +  
   (5*NC*sqrtxz3*itani2)/8. +  
   (105*NC*sqrtxz3*x*z*itani2)/16. -  
   (5*NC*sqrtxz3*itani3)/4. -  
   (105*NC*sqrtxz3*x*z*itani3)/8. - 2*NC*Li2mx +  
   NC*z*Li2mx - 2*NC*x*z*Li2mx + (3*NC*Li2x)/2. + NC*x*Li2x -  
   (3*NC*z*Li2x)/2. + NC*li2spec1 -  
   2*NC*x*li2spec1 - NC*z*li2spec1 +  
   2*NC*x*z*li2spec1 - NC*li2spec2 +  
   2*NC*x*li2spec2 + NC*z*li2spec2 -  
   2*NC*x*z*li2spec2 + NC*Li2z - 2*NC*x*Li2z +  
   NC*li2spec19 + NC*li2spec20 - NC*z*li2spec20 +  
   2*NC*x*z*li2spec20 +  
   NC*li2spec3 -  
   2*NC*x*li2spec3 -  
   NC*li2spec4 +  
   2*NC*x*li2spec4 -  
   (5*NC*sqrtxz3*atanspec1*lspec3)/4. -  
   (105*NC*sqrtxz3*x*z*atanspec1*lspec3)/8. +  
   (7*NC*lomx)/4. - 4*NC*x*lomx - (3*NC*z*lomx)/8. +  
   (NC*x*z*lomx)/4. - (67*NC*lx)/16. + (73*NC*x*lx)/16. +  
   (165*NC*z*lx)/32. + (45*NC*x*z*lx)/32. - 3*NC*l2*lx +  
   6*NC*x*l2*lx + 2*NC*z*l2*lx - 4*NC*x*z*l2*lx +  
   2*NC*lomx*lx - 2*NC*z*lomx*lx -  
   2*NC*lx*lopx + NC*z*lx*lopx -  
   2*NC*x*z*lx*lopx + (7*NC*lomz)/4. - 4*NC*x*lomz -  
   (3*NC*z*lomz)/8. + (NC*x*z*lomz)/4. -  
   NC*lomx*lomz + 2*NC*x*lomx*lomz +  
   (NC*z*lomx*lomz)/2. - NC*x*z*lomx*lomz +  
   (NC*lx*lomz)/2. - NC*x*lx*lomz -  
   (NC*z*lx*lomz)/2. + 4*NC*l2*lspec1 -  
   8*NC*x*l2*lspec1 - 3*NC*z*l2*lspec1 +  
   6*NC*x*z*l2*lspec1 + 2*NC*lx*lspec1 -  
   4*NC*x*lx*lspec1 - NC*z*lx*lspec1 +  
   2*NC*x*z*lx*lspec1 - (33*NC*lz)/16. -  
   (7*NC*x*lz)/16. + (NC*z*lz)/32. + (27*NC*x*z*lz)/32. -  
   NC*l2*lz + 2*NC*x*l2*lz + 2*NC*z*l2*lz -  
   4*NC*x*z*l2*lz - NC*lomx*lz + 2*NC*x*lomx*lz -  
   (5*NC*lx*lz)/2. - 3*NC*x*lx*lz +  
   (3*NC*z*lx*lz)/2. - 2*NC*x*z*lx*lz +  
   NC*lspec1*lz - 2*NC*x*lspec1*lz -  
   NC*z*lspec1*lz + 2*NC*x*z*lspec1*lz +  
   (5*NC*sqrtxz3*atanspec2*lspec4)/4. +  
   (105*NC*sqrtxz3*x*z*atanspec2*lspec4)/8. -  
   NC*z*l2*lspec2 + 2*NC*x*z*l2*lspec2 +  
   NC*lx*lspec2 - 2*NC*x*lx*lspec2 -  
   NC*z*lx*lspec2 + 2*NC*x*z*lx*lspec2 +  
   NC*z*lspec1*lspec2 -  
   2*NC*x*z*lspec1*lspec2 -  
   NC*z*lz*lspec2 + 2*NC*x*z*lz*lspec2 -  
   (NC*z*lx*lxpz)/2. + NC*x*z*lx*lxpz +  
   (NC*z*lz*lxpz)/2. - NC*x*z*lz*lxpz +  
   NC*lx*lopxz + NC*lz*lopxz +  
   NC*lx*lopxzi - (NC*z*lx*lopxzi)/2. +  
   NC*x*z*lx*lopxzi - NC*lz*lopxzi +  
   (NC*z*lz*lopxzi)/2. -  
   NC*x*z*lz*lopxzi + NCi/8. + (x*NCi)/2. -  
   (27*z*NCi)/16. + (x*z*NCi)/4. -  
   (3*sqrtxz3*itani1*NCi)/4. +  
   (15*sqrtxz3*x*z*itani1*NCi)/4. +  
   (3*sqrtxz3*itani2*NCi)/4. -  
   (15*sqrtxz3*x*z*itani2*NCi)/4. -  
   (3*sqrtxz3*itani3*NCi)/2. +  
   (15*sqrtxz3*x*z*itani3*NCi)/2. + Li2mx*NCi +  
   2*x*Li2mx*NCi - z*Li2mx*NCi + 2*x*z*Li2mx*NCi +  
   (Li2x*NCi)/2. - x*Li2x*NCi +  
   (3*z*Li2x*NCi)/2. - 2*li2spec1*NCi +  
   4*x*li2spec1*NCi +  
   z*li2spec1*NCi -  
   2*x*z*li2spec1*NCi +  
   2*li2spec2*NCi -  
   4*x*li2spec2*NCi -  
   z*li2spec2*NCi +  
   2*x*z*li2spec2*NCi - Li2z*NCi +  
   2*x*Li2z*NCi - 2*x*li2spec19*NCi -  
   li2spec20*NCi + z*li2spec20*NCi -  
   2*x*z*li2spec20*NCi -  
   (3*sqrtxz3*atanspec1*lspec3*NCi)/2. +  
   (15*sqrtxz3*x*z*atanspec1*lspec3*NCi)/2. -  
   (7*lomx*NCi)/4. + 4*x*lomx*NCi +  
   (3*z*lomx*NCi)/8. - (x*z*lomx*NCi)/4. +  
   (9*lx*NCi)/4. - (11*x*lx*NCi)/4. -  
   5*z*lx*NCi - (5*x*z*lx*NCi)/4. +  
   4*l2*lx*NCi - 8*x*l2*lx*NCi -  
   2*z*l2*lx*NCi + 4*x*z*l2*lx*NCi +  
   2*z*lomx*lx*NCi + lx*lopx*NCi +  
   2*x*lx*lopx*NCi - z*lx*lopx*NCi +  
   2*x*z*lx*lopx*NCi - (7*lomz*NCi)/4. +  
   4*x*lomz*NCi + (3*z*lomz*NCi)/8. -  
   (x*z*lomz*NCi)/4. + lomx*lomz*NCi -  
   2*x*lomx*lomz*NCi -  
   (z*lomx*lomz*NCi)/2. +  
   x*z*lomx*lomz*NCi - (lx*lomz*NCi)/2. +  
   x*lx*lomz*NCi + (z*lx*lomz*NCi)/2. -  
   6*l2*lspec1*NCi +  
   12*x*l2*lspec1*NCi +  
   3*z*l2*lspec1*NCi -  
   6*x*z*l2*lspec1*NCi -  
   2*lx*lspec1*NCi +  
   4*x*lx*lspec1*NCi +  
   z*lx*lspec1*NCi -  
   2*x*z*lx*lspec1*NCi + lz*NCi -  
   (3*x*lz*NCi)/4. + (z*lz*NCi)/8. -  
   x*z*lz*NCi + 4*l2*lz*NCi -  
   8*x*l2*lz*NCi - 2*z*l2*lz*NCi +  
   4*x*z*l2*lz*NCi + lomx*lz*NCi -  
   2*x*lomx*lz*NCi + lx*lz*NCi -  
   (3*z*lx*lz*NCi)/2. + 2*x*z*lx*lz*NCi -  
   2*lspec1*lz*NCi +  
   4*x*lspec1*lz*NCi +  
   z*lspec1*lz*NCi -  
   2*x*z*lspec1*lz*NCi +  
   (3*sqrtxz3*atanspec2*lspec4*NCi)/2. -  
   (15*sqrtxz3*x*z*atanspec2*lspec4*NCi)/2. -  
   2*l2*lspec2*NCi +  
   4*x*l2*lspec2*NCi +  
   z*l2*lspec2*NCi -  
   2*x*z*l2*lspec2*NCi -  
   2*lx*lspec2*NCi +  
   4*x*lx*lspec2*NCi +  
   z*lx*lspec2*NCi -  
   2*x*z*lx*lspec2*NCi +  
   2*lspec1*lspec2*NCi -  
   4*x*lspec1*lspec2*NCi -  
   z*lspec1*lspec2*NCi +  
   2*x*z*lspec1*lspec2*NCi -  
   2*lz*lspec2*NCi +  
   4*x*lz*lspec2*NCi +  
   z*lz*lspec2*NCi -  
   2*x*z*lz*lspec2*NCi -  
   (lx*lxpz*NCi)/2. + x*lx*lxpz*NCi +  
   (z*lx*lxpz*NCi)/2. - x*z*lx*lxpz*NCi +  
   (lz*lxpz*NCi)/2. - x*lz*lxpz*NCi -  
   (z*lz*lxpz*NCi)/2. + x*z*lz*lxpz*NCi -  
   2*x*lx*lopxz*NCi - 2*x*lz*lopxz*NCi -  
   (lx*lopxzi*NCi)/2. -  
   x*lx*lopxzi*NCi +  
   (z*lx*lopxzi*NCi)/2. -  
   x*z*lx*lopxzi*NCi +  
   (lz*lopxzi*NCi)/2. +  
   x*lz*lopxzi*NCi -  
   (z*lz*lopxzi*NCi)/2. +  
   x*z*lz*lopxzi*NCi - (NC*pi2)/12. -  
   (NC*x*pi2)/2. + (NC*z*pi2)/12. + (NC*x*z*pi2)/3. -  
   (NCi*pi2)/4. + (x*NCi*pi2)/2. -  
   (z*NCi*pi2)/12. - (x*z*NCi*pi2)/3. +  
   4*li2spec1*NCi*sqrtxz1i -  
   16*x*li2spec1*NCi*sqrtxz1i -  
   2*z*li2spec1*NCi*sqrtxz1i +  
   4*x*z*li2spec1*NCi*sqrtxz1i +  
   4*li2spec2*NCi*sqrtxz1i -  
   16*x*li2spec2*NCi*sqrtxz1i -  
   2*z*li2spec2*NCi*sqrtxz1i +  
   4*x*z*li2spec2*NCi*sqrtxz1i -  
   8*li2spec21*NCi*sqrtxz1i +  
   32*x*li2spec21*NCi*sqrtxz1i +  
   4*z*li2spec21*NCi*sqrtxz1i -  
   8*x*z*li2spec21*NCi*sqrtxz1i -  
   4*li2spec22*NCi*sqrtxz1i +  
   16*x*li2spec22*NCi*sqrtxz1i +  
   2*z*li2spec22*NCi*sqrtxz1i -  
   4*x*z*li2spec22*NCi*sqrtxz1i -  
   8*li2spec23*NCi*sqrtxz1i +  
   32*x*li2spec23*NCi*sqrtxz1i +  
   4*z*li2spec23*NCi*sqrtxz1i -  
   8*x*z*li2spec23*NCi*sqrtxz1i -  
   4*li2spec3*NCi* 
    sqrtxz1i + 16*x*li2spec3* 
    NCi*sqrtxz1i +  
   2*z*li2spec3*NCi* 
    sqrtxz1i - 4*x*z*li2spec3* 
    NCi*sqrtxz1i -  
   4*li2spec4*NCi* 
    sqrtxz1i + 16*x*li2spec4* 
    NCi*sqrtxz1i +  
   2*z*li2spec4*NCi* 
    sqrtxz1i - 4*x*z*li2spec4* 
    NCi*sqrtxz1i - 16*l2*NCi*sqrtxz1i +  
   64*x*l2*NCi*sqrtxz1i +  
   8*z*l2*NCi*sqrtxz1i -  
   16*x*z*l2*NCi*sqrtxz1i -  
   8*lx*NCi*sqrtxz1i +  
   32*x*lx*NCi*sqrtxz1i +  
   4*z*lx*NCi*sqrtxz1i -  
   8*x*z*lx*NCi*sqrtxz1i -  
   8*l2*lomz*NCi*sqrtxz1i +  
   32*x*l2*lomz*NCi*sqrtxz1i +  
   4*z*l2*lomz*NCi*sqrtxz1i -  
   8*x*z*l2*lomz*NCi*sqrtxz1i -  
   4*lx*lomz*NCi*sqrtxz1i +  
   16*x*lx*lomz*NCi*sqrtxz1i +  
   2*z*lx*lomz*NCi*sqrtxz1i -  
   4*x*z*lx*lomz*NCi*sqrtxz1i +  
   16*lspec1*NCi*sqrtxz1i -  
   64*x*lspec1*NCi*sqrtxz1i -  
   8*z*lspec1*NCi*sqrtxz1i +  
   16*x*z*lspec1*NCi*sqrtxz1i +  
   16*l2*lspec1*NCi*sqrtxz1i -  
   64*x*l2*lspec1*NCi*sqrtxz1i -  
   8*z*l2*lspec1*NCi*sqrtxz1i +  
   16*x*z*l2*lspec1*NCi*sqrtxz1i -  
   4*lx*lspec1*NCi*sqrtxz1i +  
   16*x*lx*lspec1*NCi*sqrtxz1i +  
   2*z*lx*lspec1*NCi*sqrtxz1i -  
   4*x*z*lx*lspec1*NCi*sqrtxz1i +  
   8*lomz*lspec1*NCi*sqrtxz1i -  
   32*x*lomz*lspec1*NCi*sqrtxz1i -  
   4*z*lomz*lspec1*NCi*sqrtxz1i +  
   8*x*z*lomz*lspec1*NCi*sqrtxz1i -  
   8*lz*NCi*sqrtxz1i +  
   32*x*lz*NCi*sqrtxz1i +  
   4*z*lz*NCi*sqrtxz1i -  
   8*x*z*lz*NCi*sqrtxz1i -  
   16*l2*lz*NCi*sqrtxz1i +  
   64*x*l2*lz*NCi*sqrtxz1i +  
   8*z*l2*lz*NCi*sqrtxz1i -  
   16*x*z*l2*lz*NCi*sqrtxz1i +  
   2*lx*lz*NCi*sqrtxz1i -  
   8*x*lx*lz*NCi*sqrtxz1i -  
   z*lx*lz*NCi*sqrtxz1i +  
   2*x*z*lx*lz*NCi*sqrtxz1i -  
   4*lomz*lz*NCi*sqrtxz1i +  
   16*x*lomz*lz*NCi*sqrtxz1i +  
   2*z*lomz*lz*NCi*sqrtxz1i -  
   4*x*z*lomz*lz*NCi*sqrtxz1i +  
   8*lspec1*lz*NCi*sqrtxz1i -  
   32*x*lspec1*lz*NCi*sqrtxz1i -  
   4*z*lspec1*lz*NCi*sqrtxz1i +  
   8*x*z*lspec1*lz*NCi*sqrtxz1i +  
   8*l2*lspec2*NCi*sqrtxz1i -  
   32*x*l2*lspec2*NCi*sqrtxz1i -  
   4*z*l2*lspec2*NCi*sqrtxz1i +  
   8*x*z*l2*lspec2*NCi*sqrtxz1i +  
   4*lx*lspec2*NCi*sqrtxz1i -  
   16*x*lx*lspec2*NCi*sqrtxz1i -  
   2*z*lx*lspec2*NCi*sqrtxz1i +  
   4*x*z*lx*lspec2*NCi*sqrtxz1i -  
   8*lspec1*lspec2*NCi*sqrtxz1i +  
   32*x*lspec1*lspec2*NCi* 
    sqrtxz1i + 4*z*lspec1*lspec2* 
    NCi*sqrtxz1i -  
   8*x*z*lspec1*lspec2*NCi* 
    sqrtxz1i + 8*lz*lspec2*NCi* 
    sqrtxz1i - 32*x*lz*lspec2*NCi* 
    sqrtxz1i - 4*z*lz*lspec2*NCi* 
    sqrtxz1i + 8*x*z*lz*lspec2*NCi* 
    sqrtxz1i + 4*lomz*lspec7*NCi* 
    sqrtxz1i - 16*x*lomz*lspec7* 
    NCi*sqrtxz1i -  
   2*z*lomz*lspec7*NCi* 
    sqrtxz1i + 4*x*z*lomz*lspec7* 
    NCi*sqrtxz1i -  
   (2*NCi*pi2*sqrtxz1i)/3. +  
   (8*x*NCi*pi2*sqrtxz1i)/3. +  
   (z*NCi*pi2*sqrtxz1i)/3. -  
   (2*x*z*NCi*pi2*sqrtxz1i)/3. -  
   (5*NC*sqrtxz3*itani1*xi2)/32. +  
   (5*NC*sqrtxz3*itani2*xi2)/32. -  
   (5*NC*sqrtxz3*itani3*xi2)/16. +  
   NC*z*Li2mx*xi2 - NC*z*Li2x*xi2 -  
   (NC*z*li2spec19*xi2)/2. -  
   (NC*z*li2spec20*xi2)/2. -  
   (5*NC*sqrtxz3*atanspec1*lspec3*xi2)/16. +  
   2*NC*z*lx*xi2 - NC*z*lomx*lx*xi2 +  
   NC*z*lx*lopx*xi2 + (NC*z*lx*lz*xi2)/2. +  
   (5*NC*sqrtxz3*atanspec2*lspec4*xi2)/16. -  
   (NC*z*lx*lopxz*xi2)/2. -  
   (NC*z*lz*lopxz*xi2)/2. -  
   (NC*z*lx*lopxzi*xi2)/2. +  
   (NC*z*lz*lopxzi*xi2)/2. -  
   Li2mx*NCi*xi2 - z*Li2mx*NCi*xi2 +  
   Li2x*NCi*xi2 + z*Li2x*NCi*xi2 +  
   (li2spec19*NCi*xi2)/2. +  
   (z*li2spec19*NCi*xi2)/2. +  
   (li2spec20*NCi*xi2)/2. +  
   (z*li2spec20*NCi*xi2)/2. -  
   2*lx*NCi*xi2 - 2*z*lx*NCi*xi2 +  
   lomx*lx*NCi*xi2 +  
   z*lomx*lx*NCi*xi2 -  
   lx*lopx*NCi*xi2 -  
   z*lx*lopx*NCi*xi2 -  
   (lx*lz*NCi*xi2)/2. -  
   (z*lx*lz*NCi*xi2)/2. +  
   (lx*lopxz*NCi*xi2)/2. +  
   (z*lx*lopxz*NCi*xi2)/2. +  
   (lz*lopxz*NCi*xi2)/2. +  
   (z*lz*lopxz*NCi*xi2)/2. +  
   (lx*lopxzi*NCi*xi2)/2. +  
   (z*lx*lopxzi*NCi*xi2)/2. -  
   (lz*lopxzi*NCi*xi2)/2. -  
   (z*lz*lopxzi*NCi*xi2)/2. +  
   (NC*z*pi2*xi2)/6. - (NCi*pi2*xi2)/6. -  
   (z*NCi*pi2*xi2)/6. + (5*NC*xi)/16. -  
   (21*NC*sqrtxz3*z*itani1*xi)/16. +  
   (21*NC*sqrtxz3*z*itani2*xi)/16. -  
   (21*NC*sqrtxz3*z*itani3*xi)/8. +  
   2*NC*z*Li2mx*xi - 2*NC*z*Li2x*xi -  
   NC*z*li2spec19*xi - NC*z*li2spec20*xi -  
   (21*NC*sqrtxz3*z*atanspec1*lspec3*xi)/8. -  
   (5*NC*lx*xi)/32. + 4*NC*z*lx*xi -  
   2*NC*z*lomx*lx*xi + 2*NC*z*lx*lopx*xi -  
   (5*NC*lz*xi)/32. + NC*z*lx*lz*xi +  
   (21*NC*sqrtxz3*z*atanspec2*lspec4*xi)/8. -  
   NC*z*lx*lopxz*xi - NC*z*lz*lopxz*xi -  
   NC*z*lx*lopxzi*xi +  
   NC*z*lz*lopxzi*xi +  
   (3*sqrtxz3*z*itani1*NCi*xi)/4. -  
   (3*sqrtxz3*z*itani2*NCi*xi)/4. +  
   (3*sqrtxz3*z*itani3*NCi*xi)/2. -  
   2*Li2mx*NCi*xi - 2*z*Li2mx*NCi*xi +  
   2*Li2x*NCi*xi + 2*z*Li2x*NCi*xi +  
   li2spec19*NCi*xi + z*li2spec19*NCi*xi +  
   li2spec20*NCi*xi +  
   z*li2spec20*NCi*xi +  
   (3*sqrtxz3*z*atanspec1*lspec3*NCi*xi)/2. -  
   4*lx*NCi*xi - 4*z*lx*NCi*xi +  
   2*lomx*lx*NCi*xi +  
   2*z*lomx*lx*NCi*xi -  
   2*lx*lopx*NCi*xi -  
   2*z*lx*lopx*NCi*xi -  
   lx*lz*NCi*xi -  
   z*lx*lz*NCi*xi -  
   (3*sqrtxz3*z*atanspec2*lspec4*NCi*xi)/2. +  
   lx*lopxz*NCi*xi +  
   z*lx*lopxz*NCi*xi +  
   lz*lopxz*NCi*xi +  
   z*lz*lopxz*NCi*xi +  
   lx*lopxzi*NCi*xi +  
   z*lx*lopxzi*NCi*xi -  
   lz*lopxzi*NCi*xi -  
   z*lz*lopxzi*NCi*xi +  
   (NC*z*pi2*xi)/3. - (NCi*pi2*xi)/3. -  
   (z*NCi*pi2*xi)/3. - (89*NC*x2)/16. -  
   (3*NC*z*x2)/2. + (35*NC*sqrtxz3*itani1*x2)/32. -  
   (35*NC*sqrtxz3*itani2*x2)/32. +  
   (35*NC*sqrtxz3*itani3*x2)/16. +  
   2*NC*Li2mx*x2 - 2*NC*z*Li2mx*x2 - NC*Li2x*x2 +  
   NC*z*Li2x*x2 - 2*NC*z*li2spec1*x2 +  
   2*NC*z*li2spec2*x2 + 2*NC*Li2z*x2 +  
   2*NC*z*li2spec19*x2 - 2*NC*li2spec20*x2 +  
   2*NC*li2spec3*x2 -  
   2*NC*li2spec4*x2 +  
   (35*NC*sqrtxz3*atanspec1*lspec3*x2)/16. +  
   (7*NC*lomx*x2)/2. + (NC*z*lomx*x2)/2. -  
   (77*NC*lx*x2)/32. - (NC*z*lx*x2)/2. -  
   2*NC*l2*lx*x2 + 4*NC*z*l2*lx*x2 +  
   NC*lomx*lx*x2 + 2*NC*lx*lopx*x2 -  
   2*NC*z*lx*lopx*x2 + (7*NC*lomz*x2)/2. +  
   (NC*z*lomz*x2)/2. - 2*NC*lomx*lomz*x2 +  
   NC*z*lomx*lomz*x2 + 2*NC*lx*lomz*x2 -  
   NC*z*lx*lomz*x2 +  
   2*NC*l2*lspec1*x2 -  
   6*NC*z*l2*lspec1*x2 +  
   2*NC*lx*lspec1*x2 -  
   2*NC*z*lx*lspec1*x2 +  
   (61*NC*lz*x2)/32. - NC*z*lz*x2 +  
   2*NC*l2*lz*x2 + 4*NC*z*l2*lz*x2 -  
   2*NC*lomx*lz*x2 + NC*lx*lz*x2 +  
   NC*z*lx*lz*x2 -  
   2*NC*z*lspec1*lz*x2 -  
   (35*NC*sqrtxz3*atanspec2*lspec4*x2)/16. -  
   2*NC*l2*lspec2*x2 -  
   2*NC*z*l2*lspec2*x2 -  
   2*NC*z*lx*lspec2*x2 +  
   2*NC*lspec1*lspec2*x2 +  
   2*NC*z*lspec1*lspec2*x2 -  
   2*NC*lz*lspec2*x2 -  
   2*NC*z*lz*lspec2*x2 -  
   NC*lx*lxpz*x2 - NC*z*lx*lxpz*x2 +  
   NC*lz*lxpz*x2 + NC*z*lz*lxpz*x2 +  
   2*NC*z*lx*lopxz*x2 +  
   2*NC*z*lz*lopxz*x2 -  
   NC*lx*lopxzi*x2 +  
   NC*z*lx*lopxzi*x2 +  
   NC*lz*lopxzi*x2 -  
   NC*z*lz*lopxzi*x2 - (9*NCi*x2)/4. +  
   (3*z*NCi*x2)/2. + 2*z*Li2mx*NCi*x2 +  
   Li2x*NCi*x2 - z*Li2x*NCi*x2 -  
   2*li2spec1*NCi*x2 +  
   2*z*li2spec1*NCi*x2 +  
   2*li2spec2*NCi*x2 -  
   2*z*li2spec2*NCi*x2 -  
   2*Li2z*NCi*x2 - 2*z*li2spec19*NCi*x2 -  
   (7*lomx*NCi*x2)/2. -  
   (z*lomx*NCi*x2)/2. +  
   (7*lx*NCi*x2)/2. + (z*lx*NCi*x2)/2. +  
   4*l2*lx*NCi*x2 -  
   4*z*l2*lx*NCi*x2 -  
   lomx*lx*NCi*x2 +  
   2*z*lx*lopx*NCi*x2 -  
   (7*lomz*NCi*x2)/2. -  
   (z*lomz*NCi*x2)/2. +  
   2*lomx*lomz*NCi*x2 -  
   z*lomx*lomz*NCi*x2 -  
   2*lx*lomz*NCi*x2 +  
   z*lx*lomz*NCi*x2 -  
   6*l2*lspec1*NCi*x2 +  
   6*z*l2*lspec1*NCi*x2 -  
   2*lx*lspec1*NCi*x2 +  
   2*z*lx*lspec1*NCi*x2 +  
   (lz*NCi*x2)/2. + z*lz*NCi*x2 +  
   4*l2*lz*NCi*x2 -  
   4*z*l2*lz*NCi*x2 +  
   2*lomx*lz*NCi*x2 -  
   z*lx*lz*NCi*x2 -  
   2*lspec1*lz*NCi*x2 +  
   2*z*lspec1*lz*NCi*x2 -  
   2*l2*lspec2*NCi*x2 +  
   2*z*l2*lspec2*NCi*x2 -  
   2*lx*lspec2*NCi*x2 +  
   2*z*lx*lspec2*NCi*x2 +  
   2*lspec1*lspec2*NCi*x2 -  
   2*z*lspec1*lspec2*NCi*x2 -  
   2*lz*lspec2*NCi*x2 +  
   2*z*lz*lspec2*NCi*x2 +  
   z*lx*lxpz*NCi*x2 -  
   z*lz*lxpz*NCi*x2 -  
   2*z*lx*lopxz*NCi*x2 -  
   2*z*lz*lopxz*NCi*x2 -  
   z*lx*lopxzi*NCi*x2 +  
   z*lz*lopxzi*NCi*x2 +  
   (NC*pi2*x2)/3. - (NC*z*pi2*x2)/2. -  
   (NCi*pi2*x2)/3. +  
   (z*NCi*pi2*x2)/2. +  
   20*li2spec1*NCi*sqrtxz1i*x2 -  
   2*z*li2spec1*NCi*sqrtxz1i*x2 +  
   20*li2spec2*NCi*sqrtxz1i*x2 -  
   2*z*li2spec2*NCi*sqrtxz1i*x2 -  
   40*li2spec21*NCi*sqrtxz1i*x2 +  
   4*z*li2spec21*NCi*sqrtxz1i*x2 -  
   20*li2spec22*NCi*sqrtxz1i*x2 +  
   2*z*li2spec22*NCi*sqrtxz1i*x2 -  
   40*li2spec23*NCi*sqrtxz1i*x2 +  
   4*z*li2spec23*NCi*sqrtxz1i*x2 -  
   20*li2spec3*NCi* 
    sqrtxz1i*x2 + 2*z* 
    li2spec3*NCi* 
    sqrtxz1i*x2 - 20* 
    li2spec4*NCi* 
    sqrtxz1i*x2 + 2*z* 
    li2spec4*NCi* 
    sqrtxz1i*x2 - 80*l2*NCi*sqrtxz1i* 
    x2 + 8*z*l2*NCi*sqrtxz1i*x2 -  
   40*lx*NCi*sqrtxz1i*x2 +  
   4*z*lx*NCi*sqrtxz1i*x2 -  
   40*l2*lomz*NCi*sqrtxz1i*x2 +  
   4*z*l2*lomz*NCi*sqrtxz1i*x2 -  
   20*lx*lomz*NCi*sqrtxz1i*x2 +  
   2*z*lx*lomz*NCi*sqrtxz1i*x2 +  
   80*lspec1*NCi*sqrtxz1i*x2 -  
   8*z*lspec1*NCi*sqrtxz1i*x2 +  
   80*l2*lspec1*NCi*sqrtxz1i*x2 -  
   8*z*l2*lspec1*NCi*sqrtxz1i*x2 -  
   20*lx*lspec1*NCi*sqrtxz1i*x2 +  
   2*z*lx*lspec1*NCi*sqrtxz1i*x2 +  
   40*lomz*lspec1*NCi*sqrtxz1i*x2 -  
   4*z*lomz*lspec1*NCi*sqrtxz1i*x2 -  
   40*lz*NCi*sqrtxz1i*x2 +  
   4*z*lz*NCi*sqrtxz1i*x2 -  
   80*l2*lz*NCi*sqrtxz1i*x2 +  
   8*z*l2*lz*NCi*sqrtxz1i*x2 +  
   10*lx*lz*NCi*sqrtxz1i*x2 -  
   z*lx*lz*NCi*sqrtxz1i*x2 -  
   20*lomz*lz*NCi*sqrtxz1i*x2 +  
   2*z*lomz*lz*NCi*sqrtxz1i*x2 +  
   40*lspec1*lz*NCi*sqrtxz1i*x2 -  
   4*z*lspec1*lz*NCi*sqrtxz1i*x2 +  
   40*l2*lspec2*NCi*sqrtxz1i*x2 -  
   4*z*l2*lspec2*NCi*sqrtxz1i*x2 +  
   20*lx*lspec2*NCi*sqrtxz1i*x2 -  
   2*z*lx*lspec2*NCi*sqrtxz1i*x2 -  
   40*lspec1*lspec2*NCi*sqrtxz1i* 
    x2 + 4*z*lspec1*lspec2*NCi* 
    sqrtxz1i*x2 + 40*lz*lspec2*NCi* 
    sqrtxz1i*x2 - 4*z*lz*lspec2*NCi* 
    sqrtxz1i*x2 + 20*lomz*lspec7* 
    NCi*sqrtxz1i*x2 -  
   2*z*lomz*lspec7*NCi*sqrtxz1i* 
    x2 - (10*NCi*pi2*sqrtxz1i*x2)/3. +  
   (z*NCi*pi2*sqrtxz1i*x2)/3. -  
   8*li2spec1*NCi*sqrtxz1i*x3 -  
   8*li2spec2*NCi*sqrtxz1i*x3 +  
   16*li2spec21*NCi*sqrtxz1i*x3 +  
   8*li2spec22*NCi*sqrtxz1i*x3 +  
   16*li2spec23*NCi*sqrtxz1i*x3 +  
   8*li2spec3*NCi* 
    sqrtxz1i*x3 + 8* 
    li2spec4*NCi* 
    sqrtxz1i*x3 + 32*l2*NCi*sqrtxz1i* 
    x3 + 16*lx*NCi*sqrtxz1i*x3 +  
   16*l2*lomz*NCi*sqrtxz1i*x3 +  
   8*lx*lomz*NCi*sqrtxz1i*x3 -  
   32*lspec1*NCi*sqrtxz1i*x3 -  
   32*l2*lspec1*NCi*sqrtxz1i*x3 +  
   8*lx*lspec1*NCi*sqrtxz1i*x3 -  
   16*lomz*lspec1*NCi*sqrtxz1i*x3 +  
   16*lz*NCi*sqrtxz1i*x3 +  
   32*l2*lz*NCi*sqrtxz1i*x3 -  
   4*lx*lz*NCi*sqrtxz1i*x3 +  
   8*lomz*lz*NCi*sqrtxz1i*x3 -  
   16*lspec1*lz*NCi*sqrtxz1i*x3 -  
   16*l2*lspec2*NCi*sqrtxz1i*x3 -  
   8*lx*lspec2*NCi*sqrtxz1i*x3 +  
   16*lspec1*lspec2*NCi*sqrtxz1i* 
    x3 - 16*lz*lspec2*NCi*sqrtxz1i* 
    x3 - 8*lomz*lspec7*NCi* 
    sqrtxz1i*x3 + (4*NCi*pi2*sqrtxz1i* 
      x3)/3. + 2*NC*Li2mx*opxi +  
   2*NC*x*Li2mx*opxi - 4*NC*z*Li2mx*opxi -  
   2*NC*x*z*Li2mx*opxi - 2*NC*Li2x*opxi -  
   2*NC*x*Li2x*opxi + 4*NC*z*Li2x*opxi +  
   2*NC*x*z*Li2x*opxi - NC*li2spec19*opxi -  
   NC*x*li2spec19*opxi + 2*NC*z*li2spec19*opxi +  
   NC*x*z*li2spec19*opxi - NC*li2spec20*opxi -  
   NC*x*li2spec20*opxi +  
   2*NC*z*li2spec20*opxi +  
   NC*x*z*li2spec20*opxi + 4*NC*lx*opxi +  
   4*NC*x*lx*opxi - 8*NC*z*lx*opxi -  
   4*NC*x*z*lx*opxi - 2*NC*lomx*lx*opxi -  
   2*NC*x*lomx*lx*opxi +  
   4*NC*z*lomx*lx*opxi +  
   2*NC*x*z*lomx*lx*opxi +  
   2*NC*lx*lopx*opxi +  
   2*NC*x*lx*lopx*opxi -  
   4*NC*z*lx*lopx*opxi -  
   2*NC*x*z*lx*lopx*opxi +  
   NC*lx*lz*opxi + NC*x*lx*lz*opxi -  
   2*NC*z*lx*lz*opxi - NC*x*z*lx*lz*opxi -  
   NC*lx*lopxz*opxi -  
   NC*x*lx*lopxz*opxi +  
   2*NC*z*lx*lopxz*opxi +  
   NC*x*z*lx*lopxz*opxi -  
   NC*lz*lopxz*opxi -  
   NC*x*lz*lopxz*opxi +  
   2*NC*z*lz*lopxz*opxi +  
   NC*x*z*lz*lopxz*opxi -  
   NC*lx*lopxzi*opxi -  
   NC*x*lx*lopxzi*opxi +  
   2*NC*z*lx*lopxzi*opxi +  
   NC*x*z*lx*lopxzi*opxi +  
   NC*lz*lopxzi*opxi +  
   NC*x*lz*lopxzi*opxi -  
   2*NC*z*lz*lopxzi*opxi -  
   NC*x*z*lz*lopxzi*opxi +  
   2*Li2mx*NCi*opxi +  
   4*z*Li2mx*NCi*opxi +  
   2*x*z*Li2mx*NCi*opxi -  
   2*Li2x*NCi*opxi - 4*z*Li2x*NCi*opxi -  
   2*x*z*Li2x*NCi*opxi -  
   li2spec19*NCi*opxi -  
   2*z*li2spec19*NCi*opxi -  
   x*z*li2spec19*NCi*opxi -  
   li2spec20*NCi*opxi -  
   2*z*li2spec20*NCi*opxi -  
   x*z*li2spec20*NCi*opxi +  
   4*lx*NCi*opxi + 8*z*lx*NCi*opxi +  
   4*x*z*lx*NCi*opxi -  
   2*lomx*lx*NCi*opxi -  
   4*z*lomx*lx*NCi*opxi -  
   2*x*z*lomx*lx*NCi*opxi +  
   2*lx*lopx*NCi*opxi +  
   4*z*lx*lopx*NCi*opxi +  
   2*x*z*lx*lopx*NCi*opxi +  
   lx*lz*NCi*opxi +  
   2*z*lx*lz*NCi*opxi +  
   x*z*lx*lz*NCi*opxi -  
   lx*lopxz*NCi*opxi -  
   2*z*lx*lopxz*NCi*opxi -  
   x*z*lx*lopxz*NCi*opxi -  
   lz*lopxz*NCi*opxi -  
   2*z*lz*lopxz*NCi*opxi -  
   x*z*lz*lopxz*NCi*opxi -  
   lx*lopxzi*NCi*opxi -  
   2*z*lx*lopxzi*NCi*opxi -  
   x*z*lx*lopxzi*NCi*opxi +  
   lz*lopxzi*NCi*opxi +  
   2*z*lz*lopxzi*NCi*opxi +  
   x*z*lz*lopxzi*NCi*opxi +  
   (NC*pi2*opxi)/3. + (NC*x*pi2*opxi)/3. -  
   (2*NC*z*pi2*opxi)/3. -  
   (NC*x*z*pi2*opxi)/3. +  
   (NCi*pi2*opxi)/3. +  
   (2*z*NCi*pi2*opxi)/3. +  
   (x*z*NCi*pi2*opxi)/3. -  
   NC*z*Li2mx*xi2*opxi +  
   NC*z*Li2x*xi2*opxi +  
   (NC*z*li2spec19*xi2*opxi)/2. +  
   (NC*z*li2spec20*xi2*opxi)/2. -  
   2*NC*z*lx*xi2*opxi +  
   NC*z*lomx*lx*xi2*opxi -  
   NC*z*lx*lopx*xi2*opxi -  
   (NC*z*lx*lz*xi2*opxi)/2. +  
   (NC*z*lx*lopxz*xi2*opxi)/2. +  
   (NC*z*lz*lopxz*xi2*opxi)/2. +  
   (NC*z*lx*lopxzi*xi2*opxi)/2. -  
   (NC*z*lz*lopxzi*xi2*opxi)/2. +  
   Li2mx*NCi*xi2*opxi +  
   z*Li2mx*NCi*xi2*opxi -  
   Li2x*NCi*xi2*opxi -  
   z*Li2x*NCi*xi2*opxi -  
   (li2spec19*NCi*xi2*opxi)/2. -  
   (z*li2spec19*NCi*xi2*opxi)/2. -  
   (li2spec20*NCi*xi2*opxi)/2. -  
   (z*li2spec20*NCi*xi2*opxi)/2. +  
   2*lx*NCi*xi2*opxi +  
   2*z*lx*NCi*xi2*opxi -  
   lomx*lx*NCi*xi2*opxi -  
   z*lomx*lx*NCi*xi2*opxi +  
   lx*lopx*NCi*xi2*opxi +  
   z*lx*lopx*NCi*xi2*opxi +  
   (lx*lz*NCi*xi2*opxi)/2. +  
   (z*lx*lz*NCi*xi2*opxi)/2. -  
   (lx*lopxz*NCi*xi2*opxi)/2. -  
   (z*lx*lopxz*NCi*xi2*opxi)/2. -  
   (lz*lopxz*NCi*xi2*opxi)/2. -  
   (z*lz*lopxz*NCi*xi2*opxi)/2. -  
   (lx*lopxzi*NCi*xi2*opxi)/2. -  
   (z*lx*lopxzi*NCi*xi2*opxi)/2. +  
   (lz*lopxzi*NCi*xi2*opxi)/2. +  
   (z*lz*lopxzi*NCi*xi2*opxi)/2. -  
   (NC*z*pi2*xi2*opxi)/6. +  
   (NCi*pi2*xi2*opxi)/6. +  
   (z*NCi*pi2*xi2*opxi)/6. -  
   3*NC*z*Li2mx*xi*opxi +  
   3*NC*z*Li2x*xi*opxi +  
   (3*NC*z*li2spec19*xi*opxi)/2. +  
   (3*NC*z*li2spec20*xi*opxi)/2. -  
   6*NC*z*lx*xi*opxi +  
   3*NC*z*lomx*lx*xi*opxi -  
   3*NC*z*lx*lopx*xi*opxi -  
   (3*NC*z*lx*lz*xi*opxi)/2. +  
   (3*NC*z*lx*lopxz*xi*opxi)/2. +  
   (3*NC*z*lz*lopxz*xi*opxi)/2. +  
   (3*NC*z*lx*lopxzi*xi*opxi)/2. -  
   (3*NC*z*lz*lopxzi*xi*opxi)/2. +  
   3*Li2mx*NCi*xi*opxi +  
   3*z*Li2mx*NCi*xi*opxi -  
   3*Li2x*NCi*xi*opxi -  
   3*z*Li2x*NCi*xi*opxi -  
   (3*li2spec19*NCi*xi*opxi)/2. -  
   (3*z*li2spec19*NCi*xi*opxi)/2. -  
   (3*li2spec20*NCi*xi*opxi)/2. -  
   (3*z*li2spec20*NCi*xi*opxi)/2. +  
   6*lx*NCi*xi*opxi +  
   6*z*lx*NCi*xi*opxi -  
   3*lomx*lx*NCi*xi*opxi -  
   3*z*lomx*lx*NCi*xi*opxi +  
   3*lx*lopx*NCi*xi*opxi +  
   3*z*lx*lopx*NCi*xi*opxi +  
   (3*lx*lz*NCi*xi*opxi)/2. +  
   (3*z*lx*lz*NCi*xi*opxi)/2. -  
   (3*lx*lopxz*NCi*xi*opxi)/2. -  
   (3*z*lx*lopxz*NCi*xi*opxi)/2. -  
   (3*lz*lopxz*NCi*xi*opxi)/2. -  
   (3*z*lz*lopxz*NCi*xi*opxi)/2. -  
   (3*lx*lopxzi*NCi*xi*opxi)/2. -  
   (3*z*lx*lopxzi*NCi*xi*opxi)/2. +  
   (3*lz*lopxzi*NCi*xi*opxi)/2. +  
   (3*z*lz*lopxzi*NCi*xi*opxi)/2. -  
   (NC*z*pi2*xi*opxi)/2. +  
   (NCi*pi2*xi*opxi)/2. +  
   (z*NCi*pi2*xi*opxi)/2. -  
   2*Li2x*NCi*omzi - 4*x*Li2x*NCi*omzi +  
   4*lx*NCi*omzi + 8*x*lx*NCi*omzi -  
   2*lomx*lx*NCi*omzi -  
   4*x*lomx*lx*NCi*omzi +  
   lx*lz*NCi*omzi +  
   2*x*lx*lz*NCi*omzi +  
   (NCi*pi2*omzi)/3. +  
   (2*x*NCi*pi2*omzi)/3. -  
   2*Li2x*NCi*x2*omzi +  
   4*lx*NCi*x2*omzi -  
   2*lomx*lx*NCi*x2*omzi +  
   lx*lz*NCi*x2*omzi +  
   (NCi*pi2*x2*omzi)/3. - (5*NC*zi2)/16. +  
   (5*NC*x*zi2)/16. - (5*NC*sqrtxz3*itani1*zi2)/ 
    32. + (5*NC*sqrtxz3*itani2*zi2)/32. -  
   (5*NC*sqrtxz3*itani3*zi2)/16. -  
   (5*NC*sqrtxz3*atanspec1*lspec3*zi2)/16. -  
   (5*NC*lx*zi2)/32. - (5*NC*x*lx*zi2)/32. +  
   (5*NC*lz*zi2)/32. - (5*NC*x*lz*zi2)/32. +  
   (5*NC*sqrtxz3*atanspec2*lspec4*zi2)/16. -  
   (11*NC*zi)/16. - (47*NC*x*zi)/8. -  
   (45*NC*sqrtxz3*x*itani1*zi)/16. +  
   (45*NC*sqrtxz3*x*itani2*zi)/16. -  
   (45*NC*sqrtxz3*x*itani3*zi)/8. +  
   NC*Li2mx*zi - 2*NC*x*Li2mx*zi -  
   (7*NC*Li2x*zi)/4. - (NC*x*Li2x*zi)/2. -  
   NC*li2spec1*zi +  
   2*NC*x*li2spec1*zi +  
   NC*li2spec2*zi -  
   2*NC*x*li2spec2*zi -  
   NC*li2spec20*zi + 2*NC*x*li2spec20*zi +  
   2*NC*sqrtxz1*l2*zi - 4*NC*sqrtxz1*x*l2*zi -  
   (45*NC*sqrtxz3*x*atanspec1*lspec3*zi)/8. -  
   (15*NC*lomx*zi)/8. + 5*NC*x*lomx*zi +  
   (59*NC*lx*zi)/16. + NC*sqrtxz1*lx*zi -  
   (109*NC*x*lx*zi)/16. - 2*NC*sqrtxz1*x*lx*zi +  
   2*NC*l2*lx*zi - 4*NC*x*l2*lx*zi -  
   (5*NC*lomx*lx*zi)/2. + NC*x*lomx*lx*zi +  
   NC*lx*lopx*zi - 2*NC*x*lx*lopx*zi -  
   (15*NC*lomz*zi)/8. + 5*NC*x*lomz*zi +  
   NC*lomx*lomz*zi -  
   2*NC*x*lomx*lomz*zi -  
   (3*NC*lx*lomz*zi)/4. +  
   (3*NC*x*lx*lomz*zi)/2. -  
   2*NC*sqrtxz1*lspec1*zi +  
   4*NC*sqrtxz1*x*lspec1*zi -  
   3*NC*l2*lspec1*zi +  
   6*NC*x*l2*lspec1*zi -  
   NC*lx*lspec1*zi +  
   2*NC*x*lx*lspec1*zi -  
   (63*NC*lz*zi)/16. + NC*sqrtxz1*lz*zi +  
   (93*NC*x*lz*zi)/16. - 2*NC*sqrtxz1*x*lz*zi +  
   2*NC*l2*lz*zi - 4*NC*x*l2*lz*zi +  
   (3*NC*lx*lz*zi)/4. - (7*NC*x*lx*lz*zi)/2. -  
   NC*lspec1*lz*zi +  
   2*NC*x*lspec1*lz*zi +  
   (45*NC*sqrtxz3*x*atanspec2*lspec4*zi)/8. -  
   NC*l2*lspec2*zi +  
   2*NC*x*l2*lspec2*zi -  
   NC*lx*lspec2*zi +  
   2*NC*x*lx*lspec2*zi +  
   NC*lspec1*lspec2*zi -  
   2*NC*x*lspec1*lspec2*zi -  
   NC*lz*lspec2*zi +  
   2*NC*x*lz*lspec2*zi -  
   (NC*lx*lxpz*zi)/2. + NC*x*lx*lxpz*zi +  
   (NC*lz*lxpz*zi)/2. - NC*x*lz*lxpz*zi -  
   (NC*lx*lopxzi*zi)/2. +  
   NC*x*lx*lopxzi*zi +  
   (NC*lz*lopxzi*zi)/2. -  
   NC*x*lz*lopxzi*zi +  
   (13*NCi*zi)/16. + (x*NCi*zi)/4. +  
   4*x*Li2mx*NCi*zi - (Li2x*NCi*zi)/4. -  
   (7*x*Li2x*NCi*zi)/2. +  
   2*li2spec1*NCi*zi +  
   2*sqrtxz1*li2spec1*NCi*zi -  
   4*x*li2spec1*NCi*zi -  
   4*sqrtxz1*x*li2spec1*NCi*zi -  
   2*li2spec2*NCi*zi +  
   2*sqrtxz1*li2spec2*NCi*zi +  
   4*x*li2spec2*NCi*zi -  
   4*sqrtxz1*x*li2spec2*NCi*zi -  
   li2spec19*NCi*zi -  
   4*sqrtxz1*li2spec21*NCi*zi +  
   8*sqrtxz1*x*li2spec21*NCi*zi -  
   2*sqrtxz1*li2spec22*NCi*zi +  
   4*sqrtxz1*x*li2spec22*NCi* 
    zi - 4*sqrtxz1*li2spec23*NCi*zi +  
   8*sqrtxz1*x*li2spec23*NCi*zi +  
   li2spec20*NCi*zi -  
   4*x*li2spec20*NCi*zi -  
   li2spec3*NCi*zi -  
   2*sqrtxz1*li2spec3*NCi* 
    zi + 2*x*li2spec3* 
    NCi*zi + 4*sqrtxz1*x* 
    li2spec3*NCi*zi +  
   li2spec4*NCi*zi -  
   2*sqrtxz1*li2spec4*NCi* 
    zi - 2*x*li2spec4* 
    NCi*zi + 4*sqrtxz1*x* 
    li2spec4*NCi*zi -  
   8*sqrtxz1*l2*NCi*zi +  
   16*sqrtxz1*x*l2*NCi*zi +  
   (15*lomx*NCi*zi)/8. -  
   5*x*lomx*NCi*zi - (7*lx*NCi*zi)/4. -  
   4*sqrtxz1*lx*NCi*zi + 13*x*lx*NCi*zi +  
   8*sqrtxz1*x*lx*NCi*zi -  
   3*l2*lx*NCi*zi +  
   6*x*l2*lx*NCi*zi +  
   (lomx*lx*NCi*zi)/2. -  
   5*x*lomx*lx*NCi*zi +  
   4*x*lx*lopx*NCi*zi +  
   (15*lomz*NCi*zi)/8. -  
   5*x*lomz*NCi*zi -  
   4*sqrtxz1*l2*lomz*NCi*zi +  
   8*sqrtxz1*x*l2*lomz*NCi*zi -  
   lomx*lomz*NCi*zi +  
   2*x*lomx*lomz*NCi*zi +  
   (3*lx*lomz*NCi*zi)/4. -  
   2*sqrtxz1*lx*lomz*NCi*zi -  
   (3*x*lx*lomz*NCi*zi)/2. +  
   4*sqrtxz1*x*lx*lomz*NCi*zi +  
   8*sqrtxz1*lspec1*NCi*zi -  
   16*sqrtxz1*x*lspec1*NCi*zi +  
   5*l2*lspec1*NCi*zi +  
   8*sqrtxz1*l2*lspec1*NCi*zi -  
   10*x*l2*lspec1*NCi*zi -  
   16*sqrtxz1*x*l2*lspec1*NCi*zi +  
   lx*lspec1*NCi*zi -  
   2*sqrtxz1*lx*lspec1*NCi*zi -  
   2*x*lx*lspec1*NCi*zi +  
   4*sqrtxz1*x*lx*lspec1*NCi*zi +  
   4*sqrtxz1*lomz*lspec1*NCi*zi -  
   8*sqrtxz1*x*lomz*lspec1*NCi*zi +  
   (23*lz*NCi*zi)/8. -  
   4*sqrtxz1*lz*NCi*zi - 7*x*lz*NCi*zi +  
   8*sqrtxz1*x*lz*NCi*zi -  
   5*l2*lz*NCi*zi -  
   8*sqrtxz1*l2*lz*NCi*zi +  
   10*x*l2*lz*NCi*zi +  
   16*sqrtxz1*x*l2*lz*NCi*zi -  
   (5*lx*lz*NCi*zi)/4. +  
   sqrtxz1*lx*lz*NCi*zi +  
   (9*x*lx*lz*NCi*zi)/2. -  
   2*sqrtxz1*x*lx*lz*NCi*zi -  
   2*sqrtxz1*lomz*lz*NCi*zi +  
   4*sqrtxz1*x*lomz*lz*NCi*zi +  
   2*lspec1*lz*NCi*zi +  
   4*sqrtxz1*lspec1*lz*NCi*zi -  
   4*x*lspec1*lz*NCi*zi -  
   8*sqrtxz1*x*lspec1*lz*NCi*zi +  
   3*l2*lspec2*NCi*zi +  
   4*sqrtxz1*l2*lspec2*NCi*zi -  
   6*x*l2*lspec2*NCi*zi -  
   8*sqrtxz1*x*l2*lspec2*NCi*zi +  
   2*lx*lspec2*NCi*zi +  
   2*sqrtxz1*lx*lspec2*NCi*zi -  
   4*x*lx*lspec2*NCi*zi -  
   4*sqrtxz1*x*lx*lspec2*NCi*zi -  
   3*lspec1*lspec2*NCi*zi -  
   4*sqrtxz1*lspec1*lspec2*NCi*zi +  
   6*x*lspec1*lspec2*NCi*zi +  
   8*sqrtxz1*x*lspec1*lspec2*NCi* 
    zi + 3*lz*lspec2*NCi*zi +  
   4*sqrtxz1*lz*lspec2*NCi*zi -  
   6*x*lz*lspec2*NCi*zi -  
   8*sqrtxz1*x*lz*lspec2*NCi*zi +  
   lx*lxpz*NCi*zi -  
   2*x*lx*lxpz*NCi*zi -  
   lz*lxpz*NCi*zi +  
   2*x*lz*lxpz*NCi*zi -  
   lx*lopxz*NCi*zi -  
   lz*lopxz*NCi*zi -  
   2*x*lx*lopxzi*NCi*zi +  
   2*x*lz*lopxzi*NCi*zi +  
   2*sqrtxz1*lomz*lspec7*NCi* 
    zi - 4*sqrtxz1*x*lomz*lspec7* 
    NCi*zi + (NC*pi2*zi)/24. +  
   (7*NC*x*pi2*zi)/12. +  
   (7*NCi*pi2*zi)/24. -  
   (sqrtxz1*NCi*pi2*zi)/3. +  
   (x*NCi*pi2*zi)/12. +  
   (2*sqrtxz1*x*NCi*pi2*zi)/3. +  
   NC*Li2mx*xi2*zi - NC*Li2x*xi2*zi -  
   (NC*li2spec19*xi2*zi)/2. -  
   (NC*li2spec20*xi2*zi)/2. +  
   2*NC*lx*xi2*zi -  
   NC*lomx*lx*xi2*zi +  
   NC*lx*lopx*xi2*zi +  
   (NC*lx*lz*xi2*zi)/2. -  
   (NC*lx*lopxz*xi2*zi)/2. -  
   (NC*lz*lopxz*xi2*zi)/2. -  
   (NC*lx*lopxzi*xi2*zi)/2. +  
   (NC*lz*lopxzi*xi2*zi)/2. -  
   2*Li2mx*NCi*xi2*zi +  
   2*Li2x*NCi*xi2*zi +  
   li2spec19*NCi*xi2*zi +  
   li2spec20*NCi*xi2*zi -  
   4*lx*NCi*xi2*zi +  
   2*lomx*lx*NCi*xi2*zi -  
   2*lx*lopx*NCi*xi2*zi -  
   lx*lz*NCi*xi2*zi +  
   lx*lopxz*NCi*xi2*zi +  
   lz*lopxz*NCi*xi2*zi +  
   lx*lopxzi*NCi*xi2*zi -  
   lz*lopxzi*NCi*xi2*zi +  
   (NC*pi2*xi2*zi)/6. -  
   (NCi*pi2*xi2*zi)/3. -  
   (5*NC*xi*zi)/16. -  
   (9*NC*sqrtxz3*itani1*xi*zi)/16. +  
   (9*NC*sqrtxz3*itani2*xi*zi)/16. -  
   (9*NC*sqrtxz3*itani3*xi*zi)/8. +  
   2*NC*Li2mx*xi*zi - 2*NC*Li2x*xi*zi -  
   NC*li2spec19*xi*zi -  
   NC*li2spec20*xi*zi -  
   (9*NC*sqrtxz3*atanspec1*lspec3*xi*zi)/8. +  
   (133*NC*lx*xi*zi)/32. -  
   2*NC*lomx*lx*xi*zi +  
   2*NC*lx*lopx*xi*zi -  
   (5*NC*lz*xi*zi)/32. +  
   NC*lx*lz*xi*zi +  
   (9*NC*sqrtxz3*atanspec2*lspec4*xi*zi)/8. -  
   NC*lx*lopxz*xi*zi -  
   NC*lz*lopxz*xi*zi -  
   NC*lx*lopxzi*xi*zi +  
   NC*lz*lopxzi*xi*zi -  
   4*Li2mx*NCi*xi*zi +  
   4*Li2x*NCi*xi*zi +  
   2*li2spec19*NCi*xi*zi +  
   2*li2spec20*NCi*xi*zi -  
   8*lx*NCi*xi*zi +  
   4*lomx*lx*NCi*xi*zi -  
   4*lx*lopx*NCi*xi*zi -  
   2*lx*lz*NCi*xi*zi +  
   2*lx*lopxz*NCi*xi*zi +  
   2*lz*lopxz*NCi*xi*zi +  
   2*lx*lopxzi*NCi*xi*zi -  
   2*lz*lopxzi*NCi*xi*zi +  
   (NC*pi2*xi*zi)/3. -  
   (2*NCi*pi2*xi*zi)/3. +  
   (117*NC*x2*zi)/16. - 2*NC*Li2mx*x2*zi +  
   NC*Li2x*x2*zi -  
   2*NC*li2spec1*x2*zi +  
   2*NC*li2spec2*x2*zi +  
   2*NC*li2spec19*x2*zi +  
   4*NC*sqrtxz1*l2*x2*zi -  
   (9*NC*lomx*x2*zi)/2. +  
   (109*NC*lx*x2*zi)/32. +  
   2*NC*sqrtxz1*lx*x2*zi +  
   4*NC*l2*lx*x2*zi -  
   NC*lomx*lx*x2*zi -  
   2*NC*lx*lopx*x2*zi -  
   (9*NC*lomz*x2*zi)/2. +  
   2*NC*lomx*lomz*x2*zi -  
   2*NC*lx*lomz*x2*zi -  
   4*NC*sqrtxz1*lspec1*x2*zi -  
   6*NC*l2*lspec1*x2*zi -  
   2*NC*lx*lspec1*x2*zi -  
   (131*NC*lz*x2*zi)/32. +  
   2*NC*sqrtxz1*lz*x2*zi +  
   4*NC*l2*lz*x2*zi +  
   NC*lx*lz*x2*zi -  
   2*NC*lspec1*lz*x2*zi -  
   2*NC*l2*lspec2*x2*zi -  
   2*NC*lx*lspec2*x2*zi +  
   2*NC*lspec1*lspec2*x2*zi -  
   2*NC*lz*lspec2*x2*zi -  
   NC*lx*lxpz*x2*zi +  
   NC*lz*lxpz*x2*zi +  
   2*NC*lx*lopxz*x2*zi +  
   2*NC*lz*lopxz*x2*zi +  
   NC*lx*lopxzi*x2*zi -  
   NC*lz*lopxzi*x2*zi +  
   (NCi*x2*zi)/2. +  
   2*Li2mx*NCi*x2*zi -  
   3*Li2x*NCi*x2*zi +  
   2*li2spec1*NCi*x2*zi +  
   2*sqrtxz1*li2spec1*NCi*x2*zi -  
   2*li2spec2*NCi*x2*zi +  
   2*sqrtxz1*li2spec2*NCi*x2*zi -  
   2*li2spec19*NCi*x2*zi -  
   4*sqrtxz1*li2spec21*NCi*x2*zi -  
   2*sqrtxz1*li2spec22*NCi*x2*zi -  
   4*sqrtxz1*li2spec23*NCi*x2*zi -  
   2*sqrtxz1*li2spec3*NCi* 
    x2*zi - 2*sqrtxz1* 
    li2spec4*NCi*x2* 
    zi - 8*sqrtxz1*l2*NCi*x2*zi +  
   (9*lomx*NCi*x2*zi)/2. -  
   (lx*NCi*x2*zi)/2. -  
   4*sqrtxz1*lx*NCi*x2*zi -  
   4*l2*lx*NCi*x2*zi -  
   lomx*lx*NCi*x2*zi +  
   2*lx*lopx*NCi*x2*zi +  
   (9*lomz*NCi*x2*zi)/2. -  
   4*sqrtxz1*l2*lomz*NCi*x2*zi -  
   2*lomx*lomz*NCi*x2*zi +  
   2*lx*lomz*NCi*x2*zi -  
   2*sqrtxz1*lx*lomz*NCi*x2*zi +  
   8*sqrtxz1*lspec1*NCi*x2*zi +  
   6*l2*lspec1*NCi*x2*zi +  
   8*sqrtxz1*l2*lspec1*NCi*x2*zi +  
   2*lx*lspec1*NCi*x2*zi -  
   2*sqrtxz1*lx*lspec1*NCi*x2*zi +  
   4*sqrtxz1*lomz*lspec1*NCi*x2*zi +  
   (13*lz*NCi*x2*zi)/2. -  
   4*sqrtxz1*lz*NCi*x2*zi -  
   4*l2*lz*NCi*x2*zi -  
   8*sqrtxz1*l2*lz*NCi*x2*zi +  
   sqrtxz1*lx*lz*NCi*x2*zi -  
   2*sqrtxz1*lomz*lz*NCi*x2*zi +  
   2*lspec1*lz*NCi*x2*zi +  
   4*sqrtxz1*lspec1*lz*NCi*x2*zi +  
   2*l2*lspec2*NCi*x2*zi +  
   4*sqrtxz1*l2*lspec2*NCi*x2*zi +  
   2*lx*lspec2*NCi*x2*zi +  
   2*sqrtxz1*lx*lspec2*NCi*x2*zi -  
   2*lspec1*lspec2*NCi*x2* 
    zi - 4*sqrtxz1*lspec1*lspec2* 
    NCi*x2*zi +  
   2*lz*lspec2*NCi*x2*zi +  
   4*sqrtxz1*lz*lspec2*NCi*x2*zi +  
   lx*lxpz*NCi*x2*zi -  
   lz*lxpz*NCi*x2*zi -  
   2*lx*lopxz*NCi*x2*zi -  
   2*lz*lopxz*NCi*x2*zi -  
   lx*lopxzi*NCi*x2*zi +  
   lz*lopxzi*NCi*x2*zi +  
   2*sqrtxz1*lomz*lspec7*NCi*x2* 
    zi - (2*NC*pi2*x2*zi)/3. +  
   NCi*pi2*x2*zi -  
   (sqrtxz1*NCi*pi2*x2*zi)/3. -  
   4*NC*Li2mx*opxi*zi -  
   2*NC*x*Li2mx*opxi*zi +  
   4*NC*Li2x*opxi*zi +  
   2*NC*x*Li2x*opxi*zi +  
   2*NC*li2spec19*opxi*zi +  
   NC*x*li2spec19*opxi*zi +  
   2*NC*li2spec20*opxi*zi +  
   NC*x*li2spec20*opxi*zi -  
   8*NC*lx*opxi*zi -  
   4*NC*x*lx*opxi*zi +  
   4*NC*lomx*lx*opxi*zi +  
   2*NC*x*lomx*lx*opxi*zi -  
   4*NC*lx*lopx*opxi*zi -  
   2*NC*x*lx*lopx*opxi*zi -  
   2*NC*lx*lz*opxi*zi -  
   NC*x*lx*lz*opxi*zi +  
   2*NC*lx*lopxz*opxi*zi +  
   NC*x*lx*lopxz*opxi*zi +  
   2*NC*lz*lopxz*opxi*zi +  
   NC*x*lz*lopxz*opxi*zi +  
   2*NC*lx*lopxzi*opxi*zi +  
   NC*x*lx*lopxzi*opxi*zi -  
   2*NC*lz*lopxzi*opxi*zi -  
   NC*x*lz*lopxzi*opxi*zi +  
   6*Li2mx*NCi*opxi*zi +  
   2*x*Li2mx*NCi*opxi*zi -  
   6*Li2x*NCi*opxi*zi -  
   2*x*Li2x*NCi*opxi*zi -  
   3*li2spec19*NCi*opxi*zi -  
   x*li2spec19*NCi*opxi*zi -  
   3*li2spec20*NCi*opxi*zi -  
   x*li2spec20*NCi*opxi*zi +  
   12*lx*NCi*opxi*zi +  
   4*x*lx*NCi*opxi*zi -  
   6*lomx*lx*NCi*opxi*zi -  
   2*x*lomx*lx*NCi*opxi*zi +  
   6*lx*lopx*NCi*opxi*zi +  
   2*x*lx*lopx*NCi*opxi*zi +  
   3*lx*lz*NCi*opxi*zi +  
   x*lx*lz*NCi*opxi*zi -  
   3*lx*lopxz*NCi*opxi*zi -  
   x*lx*lopxz*NCi*opxi*zi -  
   3*lz*lopxz*NCi*opxi*zi -  
   x*lz*lopxz*NCi*opxi*zi -  
   3*lx*lopxzi*NCi*opxi*zi -  
   x*lx*lopxzi*NCi*opxi*zi +  
   3*lz*lopxzi*NCi*opxi*zi +  
   x*lz*lopxzi*NCi*opxi*zi -  
   (2*NC*pi2*opxi*zi)/3. -  
   (NC*x*pi2*opxi*zi)/3. +  
   NCi*pi2*opxi*zi +  
   (x*NCi*pi2*opxi*zi)/3. -  
   NC*Li2mx*xi2*opxi*zi +  
   NC*Li2x*xi2*opxi*zi +  
   (NC*li2spec19*xi2*opxi*zi)/2. +  
   (NC*li2spec20*xi2*opxi*zi)/2. -  
   2*NC*lx*xi2*opxi*zi +  
   NC*lomx*lx*xi2*opxi*zi -  
   NC*lx*lopx*xi2*opxi*zi -  
   (NC*lx*lz*xi2*opxi*zi)/2. +  
   (NC*lx*lopxz*xi2*opxi*zi)/2. +  
   (NC*lz*lopxz*xi2*opxi*zi)/2. +  
   (NC*lx*lopxzi*xi2*opxi*zi)/2. -  
   (NC*lz*lopxzi*xi2*opxi*zi)/2. +  
   2*Li2mx*NCi*xi2*opxi*zi -  
   2*Li2x*NCi*xi2*opxi*zi -  
   li2spec19*NCi*xi2*opxi*zi -  
   li2spec20*NCi*xi2*opxi*zi +  
   4*lx*NCi*xi2*opxi*zi -  
   2*lomx*lx*NCi*xi2*opxi*zi +  
   2*lx*lopx*NCi*xi2*opxi*zi +  
   lx*lz*NCi*xi2*opxi*zi -  
   lx*lopxz*NCi*xi2*opxi*zi -  
   lz*lopxz*NCi*xi2*opxi*zi -  
   lx*lopxzi*NCi*xi2*opxi*zi +  
   lz*lopxzi*NCi*xi2*opxi*zi -  
   (NC*pi2*xi2*opxi*zi)/6. +  
   (NCi*pi2*xi2*opxi*zi)/3. -  
   3*NC*Li2mx*xi*opxi*zi +  
   3*NC*Li2x*xi*opxi*zi +  
   (3*NC*li2spec19*xi*opxi*zi)/2. +  
   (3*NC*li2spec20*xi*opxi*zi)/2. -  
   6*NC*lx*xi*opxi*zi +  
   3*NC*lomx*lx*xi*opxi*zi -  
   3*NC*lx*lopx*xi*opxi*zi -  
   (3*NC*lx*lz*xi*opxi*zi)/2. +  
   (3*NC*lx*lopxz*xi*opxi*zi)/2. +  
   (3*NC*lz*lopxz*xi*opxi*zi)/2. +  
   (3*NC*lx*lopxzi*xi*opxi*zi)/2. -  
   (3*NC*lz*lopxzi*xi*opxi*zi)/2. +  
   6*Li2mx*NCi*xi*opxi*zi -  
   6*Li2x*NCi*xi*opxi*zi -  
   3*li2spec19*NCi*xi*opxi*zi -  
   3*li2spec20*NCi*xi*opxi*zi +  
   12*lx*NCi*xi*opxi*zi -  
   6*lomx*lx*NCi*xi*opxi*zi +  
   6*lx*lopx*NCi*xi*opxi*zi +  
   3*lx*lz*NCi*xi*opxi*zi -  
   3*lx*lopxz*NCi*xi*opxi*zi -  
   3*lz*lopxz*NCi*xi*opxi*zi -  
   3*lx*lopxzi*NCi*xi*opxi* 
    zi + 3*lz*lopxzi*NCi*xi* 
    opxi*zi - (NC*pi2*xi*opxi* 
      zi)/2. + NCi*pi2*xi*opxi* 
    zi - 4*x*Li2mx*NCi*omzi*zi +  
   4*x*Li2x*NCi*omzi*zi +  
   2*x*li2spec19*NCi*omzi*zi +  
   2*x*li2spec20*NCi*omzi*zi -  
   8*x*lx*NCi*omzi*zi +  
   4*x*lomx*lx*NCi*omzi*zi -  
   4*x*lx*lopx*NCi*omzi*zi -  
   2*x*lx*lz*NCi*omzi*zi +  
   2*x*lx*lopxz*NCi*omzi*zi +  
   2*x*lz*lopxz*NCi*omzi*zi +  
   2*x*lx*lopxzi*NCi*omzi*zi -  
   2*x*lz*lopxzi*NCi*omzi*zi -  
   (2*x*NCi*pi2*omzi*zi)/3. +  
   2*Li2mx*NCi*xi2*omzi*zi -  
   2*Li2x*NCi*xi2*omzi*zi -  
   li2spec19*NCi*xi2*omzi*zi -  
   li2spec20*NCi*xi2*omzi*zi +  
   4*lx*NCi*xi2*omzi*zi -  
   2*lomx*lx*NCi*xi2*omzi*zi +  
   2*lx*lopx*NCi*xi2*omzi*zi +  
   lx*lz*NCi*xi2*omzi*zi -  
   lx*lopxz*NCi*xi2*omzi*zi -  
   lz*lopxz*NCi*xi2*omzi*zi -  
   lx*lopxzi*NCi*xi2*omzi*zi +  
   lz*lopxzi*NCi*xi2*omzi*zi +  
   (NCi*pi2*xi2*omzi*zi)/3. +  
   4*Li2mx*NCi*xi*omzi*zi -  
   4*Li2x*NCi*xi*omzi*zi -  
   2*li2spec19*NCi*xi*omzi*zi -  
   2*li2spec20*NCi*xi*omzi*zi +  
   8*lx*NCi*xi*omzi*zi -  
   4*lomx*lx*NCi*xi*omzi*zi +  
   4*lx*lopx*NCi*xi*omzi*zi +  
   2*lx*lz*NCi*xi*omzi*zi -  
   2*lx*lopxz*NCi*xi*omzi*zi -  
   2*lz*lopxz*NCi*xi*omzi*zi -  
   2*lx*lopxzi*NCi*xi*omzi* 
    zi + 2*lz*lopxzi*NCi*xi* 
    omzi*zi + (2*NCi*pi2*xi*omzi* 
      zi)/3. - 2*Li2mx*NCi*x2*omzi*zi +  
   2*Li2x*NCi*x2*omzi*zi +  
   li2spec19*NCi*x2*omzi*zi +  
   li2spec20*NCi*x2*omzi*zi -  
   4*lx*NCi*x2*omzi*zi +  
   2*lomx*lx*NCi*x2*omzi*zi -  
   2*lx*lopx*NCi*x2*omzi*zi -  
   lx*lz*NCi*x2*omzi*zi +  
   lx*lopxz*NCi*x2*omzi*zi +  
   lz*lopxz*NCi*x2*omzi*zi +  
   lx*lopxzi*NCi*x2*omzi*zi -  
   lz*lopxzi*NCi*x2*omzi*zi -  
   (NCi*pi2*x2*omzi*zi)/3. -  
   6*Li2mx*NCi*opxi*omzi*zi -  
   2*x*Li2mx*NCi*opxi*omzi*zi +  
   6*Li2x*NCi*opxi*omzi*zi +  
   2*x*Li2x*NCi*opxi*omzi*zi +  
   3*li2spec19*NCi*opxi*omzi*zi +  
   x*li2spec19*NCi*opxi*omzi*zi +  
   3*li2spec20*NCi*opxi*omzi*zi +  
   x*li2spec20*NCi*opxi*omzi*zi -  
   12*lx*NCi*opxi*omzi*zi -  
   4*x*lx*NCi*opxi*omzi*zi +  
   6*lomx*lx*NCi*opxi*omzi*zi +  
   2*x*lomx*lx*NCi*opxi*omzi*zi -  
   6*lx*lopx*NCi*opxi*omzi*zi -  
   2*x*lx*lopx*NCi*opxi*omzi*zi -  
   3*lx*lz*NCi*opxi*omzi*zi -  
   x*lx*lz*NCi*opxi*omzi*zi +  
   3*lx*lopxz*NCi*opxi*omzi*zi +  
   x*lx*lopxz*NCi*opxi*omzi*zi +  
   3*lz*lopxz*NCi*opxi*omzi*zi +  
   x*lz*lopxz*NCi*opxi*omzi*zi +  
   3*lx*lopxzi*NCi*opxi*omzi* 
    zi + x*lx*lopxzi*NCi*opxi* 
    omzi*zi - 3*lz*lopxzi*NCi* 
    opxi*omzi*zi -  
   x*lz*lopxzi*NCi*opxi*omzi* 
    zi - NCi*pi2*opxi*omzi*zi -  
   (x*NCi*pi2*opxi*omzi*zi)/3. -  
   2*Li2mx*NCi*xi2*opxi*omzi*zi +  
   2*Li2x*NCi*xi2*opxi*omzi*zi +  
   li2spec19*NCi*xi2*opxi*omzi*zi +  
   li2spec20*NCi*xi2*opxi*omzi* 
    zi - 4*lx*NCi*xi2*opxi*omzi* 
    zi + 2*lomx*lx*NCi*xi2*opxi* 
    omzi*zi - 2*lx*lopx*NCi*xi2* 
    opxi*omzi*zi -  
   lx*lz*NCi*xi2*opxi*omzi*zi +  
   lx*lopxz*NCi*xi2*opxi*omzi* 
    zi + lz*lopxz*NCi*xi2*opxi* 
    omzi*zi + lx*lopxzi*NCi* 
    xi2*opxi*omzi*zi -  
   lz*lopxzi*NCi*xi2*opxi* 
    omzi*zi - (NCi*pi2*xi2*opxi* 
      omzi*zi)/3. -  
   6*Li2mx*NCi*xi*opxi*omzi*zi +  
   6*Li2x*NCi*xi*opxi*omzi*zi +  
   3*li2spec19*NCi*xi*opxi*omzi*zi +  
   3*li2spec20*NCi*xi*opxi*omzi* 
    zi - 12*lx*NCi*xi*opxi*omzi* 
    zi + 6*lomx*lx*NCi*xi*opxi* 
    omzi*zi - 6*lx*lopx*NCi*xi* 
    opxi*omzi*zi -  
   3*lx*lz*NCi*xi*opxi*omzi* 
    zi + 3*lx*lopxz*NCi*xi*opxi* 
    omzi*zi + 3*lz*lopxz*NCi*xi* 
    opxi*omzi*zi +  
   3*lx*lopxzi*NCi*xi*opxi* 
    omzi*zi - 3*lz*lopxzi*NCi* 
    xi*opxi*omzi*zi -  
   NCi*pi2*xi*opxi*omzi*zi +  
   (19*NC*sqrtxz3*itani1*z2)/32. -  
   (19*NC*sqrtxz3*itani2*z2)/32. +  
   (19*NC*sqrtxz3*itani3*z2)/16. +  
   (19*NC*sqrtxz3*atanspec1*lspec3*z2)/16. -  
   (19*NC*sqrtxz3*atanspec2*lspec4*z2)/16. -  
   (3*sqrtxz3*itani1*NCi*z2)/4. +  
   (3*sqrtxz3*itani2*NCi*z2)/4. -  
   (3*sqrtxz3*itani3*NCi*z2)/2. -  
   (3*sqrtxz3*atanspec1*lspec3*NCi*z2)/2. +  
   (3*sqrtxz3*atanspec2*lspec4*NCi*z2)/2. -  
   6*li2spec1*NCi*sqrtxz1i*opzi +  
   20*x*li2spec1*NCi*sqrtxz1i* 
    opzi + 2*z*li2spec1*NCi* 
    sqrtxz1i*opzi -  
   4*x*z*li2spec1*NCi*sqrtxz1i* 
    opzi - 6*li2spec2*NCi*sqrtxz1i* 
    opzi + 20*x*li2spec2*NCi* 
    sqrtxz1i*opzi +  
   2*z*li2spec2*NCi*sqrtxz1i* 
    opzi - 4*x*z*li2spec2*NCi* 
    sqrtxz1i*opzi +  
   12*li2spec21*NCi*sqrtxz1i*opzi -  
   40*x*li2spec21*NCi*sqrtxz1i* 
    opzi - 4*z*li2spec21*NCi* 
    sqrtxz1i*opzi +  
   8*x*z*li2spec21*NCi*sqrtxz1i* 
    opzi + 6*li2spec22*NCi* 
    sqrtxz1i*opzi -  
   20*x*li2spec22*NCi*sqrtxz1i*opzi -  
   2*z*li2spec22*NCi*sqrtxz1i*opzi +  
   4*x*z*li2spec22*NCi*sqrtxz1i*opzi +  
   12*li2spec23*NCi*sqrtxz1i*opzi -  
   40*x*li2spec23*NCi*sqrtxz1i* 
    opzi - 4*z*li2spec23*NCi* 
    sqrtxz1i*opzi +  
   8*x*z*li2spec23*NCi*sqrtxz1i* 
    opzi + 6*li2spec3* 
    NCi*sqrtxz1i*opzi -  
   20*x*li2spec3*NCi* 
    sqrtxz1i*opzi -  
   2*z*li2spec3*NCi* 
    sqrtxz1i*opzi +  
   4*x*z*li2spec3*NCi* 
    sqrtxz1i*opzi +  
   6*li2spec4*NCi* 
    sqrtxz1i*opzi -  
   20*x*li2spec4*NCi* 
    sqrtxz1i*opzi -  
   2*z*li2spec4*NCi* 
    sqrtxz1i*opzi +  
   4*x*z*li2spec4*NCi* 
    sqrtxz1i*opzi +  
   24*l2*NCi*sqrtxz1i*opzi -  
   80*x*l2*NCi*sqrtxz1i*opzi -  
   8*z*l2*NCi*sqrtxz1i*opzi +  
   16*x*z*l2*NCi*sqrtxz1i*opzi +  
   12*lx*NCi*sqrtxz1i*opzi -  
   40*x*lx*NCi*sqrtxz1i*opzi -  
   4*z*lx*NCi*sqrtxz1i*opzi +  
   8*x*z*lx*NCi*sqrtxz1i*opzi +  
   12*l2*lomz*NCi*sqrtxz1i*opzi -  
   40*x*l2*lomz*NCi*sqrtxz1i*opzi -  
   4*z*l2*lomz*NCi*sqrtxz1i*opzi +  
   8*x*z*l2*lomz*NCi*sqrtxz1i*opzi +  
   6*lx*lomz*NCi*sqrtxz1i*opzi -  
   20*x*lx*lomz*NCi*sqrtxz1i*opzi -  
   2*z*lx*lomz*NCi*sqrtxz1i*opzi +  
   4*x*z*lx*lomz*NCi*sqrtxz1i*opzi -  
   24*lspec1*NCi*sqrtxz1i*opzi +  
   80*x*lspec1*NCi*sqrtxz1i*opzi +  
   8*z*lspec1*NCi*sqrtxz1i*opzi -  
   16*x*z*lspec1*NCi*sqrtxz1i*opzi -  
   24*l2*lspec1*NCi*sqrtxz1i*opzi +  
   80*x*l2*lspec1*NCi*sqrtxz1i* 
    opzi + 8*z*l2*lspec1*NCi* 
    sqrtxz1i*opzi -  
   16*x*z*l2*lspec1*NCi*sqrtxz1i* 
    opzi + 6*lx*lspec1*NCi*sqrtxz1i* 
    opzi - 20*x*lx*lspec1*NCi* 
    sqrtxz1i*opzi -  
   2*z*lx*lspec1*NCi*sqrtxz1i*opzi +  
   4*x*z*lx*lspec1*NCi*sqrtxz1i* 
    opzi - 12*lomz*lspec1*NCi* 
    sqrtxz1i*opzi +  
   40*x*lomz*lspec1*NCi*sqrtxz1i* 
    opzi + 4*z*lomz*lspec1*NCi* 
    sqrtxz1i*opzi -  
   8*x*z*lomz*lspec1*NCi*sqrtxz1i* 
    opzi + 12*lz*NCi*sqrtxz1i*opzi -  
   40*x*lz*NCi*sqrtxz1i*opzi -  
   4*z*lz*NCi*sqrtxz1i*opzi +  
   8*x*z*lz*NCi*sqrtxz1i*opzi +  
   24*l2*lz*NCi*sqrtxz1i*opzi -  
   80*x*l2*lz*NCi*sqrtxz1i*opzi -  
   8*z*l2*lz*NCi*sqrtxz1i*opzi +  
   16*x*z*l2*lz*NCi*sqrtxz1i*opzi -  
   3*lx*lz*NCi*sqrtxz1i*opzi +  
   10*x*lx*lz*NCi*sqrtxz1i*opzi +  
   z*lx*lz*NCi*sqrtxz1i*opzi -  
   2*x*z*lx*lz*NCi*sqrtxz1i*opzi +  
   6*lomz*lz*NCi*sqrtxz1i*opzi -  
   20*x*lomz*lz*NCi*sqrtxz1i*opzi -  
   2*z*lomz*lz*NCi*sqrtxz1i*opzi +  
   4*x*z*lomz*lz*NCi*sqrtxz1i*opzi -  
   12*lspec1*lz*NCi*sqrtxz1i*opzi +  
   40*x*lspec1*lz*NCi*sqrtxz1i* 
    opzi + 4*z*lspec1*lz*NCi* 
    sqrtxz1i*opzi -  
   8*x*z*lspec1*lz*NCi*sqrtxz1i* 
    opzi - 12*l2*lspec2*NCi*sqrtxz1i* 
    opzi + 40*x*l2*lspec2*NCi* 
    sqrtxz1i*opzi +  
   4*z*l2*lspec2*NCi*sqrtxz1i*opzi -  
   8*x*z*l2*lspec2*NCi*sqrtxz1i* 
    opzi - 6*lx*lspec2*NCi*sqrtxz1i* 
    opzi + 20*x*lx*lspec2*NCi* 
    sqrtxz1i*opzi +  
   2*z*lx*lspec2*NCi*sqrtxz1i*opzi -  
   4*x*z*lx*lspec2*NCi*sqrtxz1i* 
    opzi + 12*lspec1*lspec2*NCi* 
    sqrtxz1i*opzi -  
   40*x*lspec1*lspec2*NCi*sqrtxz1i* 
    opzi - 4*z*lspec1*lspec2*NCi* 
    sqrtxz1i*opzi +  
   8*x*z*lspec1*lspec2*NCi*sqrtxz1i* 
    opzi - 12*lz*lspec2*NCi*sqrtxz1i* 
    opzi + 40*x*lz*lspec2*NCi* 
    sqrtxz1i*opzi +  
   4*z*lz*lspec2*NCi*sqrtxz1i*opzi -  
   8*x*z*lz*lspec2*NCi*sqrtxz1i* 
    opzi - 6*lomz*lspec7*NCi* 
    sqrtxz1i*opzi +  
   20*x*lomz*lspec7*NCi*sqrtxz1i* 
    opzi + 2*z*lomz*lspec7*NCi* 
    sqrtxz1i*opzi -  
   4*x*z*lomz*lspec7*NCi* 
    sqrtxz1i*opzi +  
   NCi*pi2*sqrtxz1i*opzi -  
   (10*x*NCi*pi2*sqrtxz1i*opzi)/3. -  
   (z*NCi*pi2*sqrtxz1i*opzi)/3. +  
   (2*x*z*NCi*pi2*sqrtxz1i*opzi)/3. -  
   22*li2spec1*NCi*sqrtxz1i*x2* 
    opzi + 2*z*li2spec1*NCi* 
    sqrtxz1i*x2*opzi -  
   22*li2spec2*NCi*sqrtxz1i*x2* 
    opzi + 2*z*li2spec2*NCi* 
    sqrtxz1i*x2*opzi +  
   44*li2spec21*NCi*sqrtxz1i*x2* 
    opzi - 4*z*li2spec21*NCi* 
    sqrtxz1i*x2*opzi +  
   22*li2spec22*NCi*sqrtxz1i*x2* 
    opzi - 2*z*li2spec22*NCi* 
    sqrtxz1i*x2*opzi +  
   44*li2spec23*NCi*sqrtxz1i*x2* 
    opzi - 4*z*li2spec23*NCi* 
    sqrtxz1i*x2*opzi +  
   22*li2spec3*NCi* 
    sqrtxz1i*x2*opzi -  
   2*z*li2spec3*NCi* 
    sqrtxz1i*x2*opzi +  
   22*li2spec4*NCi* 
    sqrtxz1i*x2*opzi -  
   2*z*li2spec4*NCi* 
    sqrtxz1i*x2*opzi +  
   88*l2*NCi*sqrtxz1i*x2*opzi -  
   8*z*l2*NCi*sqrtxz1i*x2*opzi +  
   44*lx*NCi*sqrtxz1i*x2*opzi -  
   4*z*lx*NCi*sqrtxz1i*x2*opzi +  
   44*l2*lomz*NCi*sqrtxz1i*x2*opzi -  
   4*z*l2*lomz*NCi*sqrtxz1i*x2*opzi +  
   22*lx*lomz*NCi*sqrtxz1i*x2*opzi -  
   2*z*lx*lomz*NCi*sqrtxz1i*x2*opzi -  
   88*lspec1*NCi*sqrtxz1i*x2* 
    opzi + 8*z*lspec1*NCi*sqrtxz1i* 
    x2*opzi - 88*l2*lspec1*NCi* 
    sqrtxz1i*x2*opzi +  
   8*z*l2*lspec1*NCi*sqrtxz1i*x2* 
    opzi + 22*lx*lspec1*NCi*sqrtxz1i* 
    x2*opzi - 2*z*lx*lspec1*NCi* 
    sqrtxz1i*x2*opzi -  
   44*lomz*lspec1*NCi*sqrtxz1i*x2* 
    opzi + 4*z*lomz*lspec1*NCi* 
    sqrtxz1i*x2*opzi +  
   44*lz*NCi*sqrtxz1i*x2*opzi -  
   4*z*lz*NCi*sqrtxz1i*x2*opzi +  
   88*l2*lz*NCi*sqrtxz1i*x2*opzi -  
   8*z*l2*lz*NCi*sqrtxz1i*x2*opzi -  
   11*lx*lz*NCi*sqrtxz1i*x2*opzi +  
   z*lx*lz*NCi*sqrtxz1i*x2*opzi +  
   22*lomz*lz*NCi*sqrtxz1i*x2*opzi -  
   2*z*lomz*lz*NCi*sqrtxz1i*x2*opzi -  
   44*lspec1*lz*NCi*sqrtxz1i*x2* 
    opzi + 4*z*lspec1*lz*NCi* 
    sqrtxz1i*x2*opzi -  
   44*l2*lspec2*NCi*sqrtxz1i*x2* 
    opzi + 4*z*l2*lspec2*NCi* 
    sqrtxz1i*x2*opzi -  
   22*lx*lspec2*NCi*sqrtxz1i*x2* 
    opzi + 2*z*lx*lspec2*NCi* 
    sqrtxz1i*x2*opzi +  
   44*lspec1*lspec2*NCi*sqrtxz1i* 
    x2*opzi - 4*z*lspec1*lspec2* 
    NCi*sqrtxz1i*x2*opzi -  
   44*lz*lspec2*NCi*sqrtxz1i*x2* 
    opzi + 4*z*lz*lspec2*NCi* 
    sqrtxz1i*x2*opzi -  
   22*lomz*lspec7*NCi*sqrtxz1i* 
    x2*opzi + 2*z*lomz*lspec7* 
    NCi*sqrtxz1i*x2*opzi +  
   (11*NCi*pi2*sqrtxz1i*x2*opzi)/3. -  
   (z*NCi*pi2*sqrtxz1i*x2*opzi)/3. +  
   8*li2spec1*NCi*sqrtxz1i*x3* 
    opzi + 8*li2spec2*NCi*sqrtxz1i* 
    x3*opzi - 16*li2spec21*NCi* 
    sqrtxz1i*x3*opzi -  
   8*li2spec22*NCi*sqrtxz1i*x3* 
    opzi - 16*li2spec23*NCi*sqrtxz1i* 
    x3*opzi - 8*li2spec3*NCi*sqrtxz1i*x3* 
    opzi - 8*li2spec4* 
    NCi*sqrtxz1i*x3*opzi -  
   32*l2*NCi*sqrtxz1i*x3*opzi -  
   16*lx*NCi*sqrtxz1i*x3*opzi -  
   16*l2*lomz*NCi*sqrtxz1i*x3*opzi -  
   8*lx*lomz*NCi*sqrtxz1i*x3*opzi +  
   32*lspec1*NCi*sqrtxz1i*x3* 
    opzi + 32*l2*lspec1*NCi*sqrtxz1i* 
    x3*opzi - 8*lx*lspec1*NCi* 
    sqrtxz1i*x3*opzi +  
   16*lomz*lspec1*NCi*sqrtxz1i*x3* 
    opzi - 16*lz*NCi*sqrtxz1i*x3* 
    opzi - 32*l2*lz*NCi*sqrtxz1i*x3* 
    opzi + 4*lx*lz*NCi*sqrtxz1i*x3* 
    opzi - 8*lomz*lz*NCi*sqrtxz1i*x3* 
    opzi + 16*lspec1*lz*NCi*sqrtxz1i* 
    x3*opzi + 16*l2*lspec2*NCi* 
    sqrtxz1i*x3*opzi +  
   8*lx*lspec2*NCi*sqrtxz1i*x3* 
    opzi - 16*lspec1*lspec2*NCi* 
    sqrtxz1i*x3*opzi +  
   16*lz*lspec2*NCi*sqrtxz1i*x3* 
    opzi + 8*lomz*lspec7*NCi* 
    sqrtxz1i*x3*opzi -  
   (4*NCi*pi2*sqrtxz1i*x3*opzi)/3. -  
   2*li2spec1*NCi*zi*opzi -  
   2*sqrtxz1*li2spec1*NCi*zi* 
    opzi + 4*x*li2spec1*NCi*zi* 
    opzi + 4*sqrtxz1*x*li2spec1*NCi* 
    zi*opzi + 2*li2spec2*NCi* 
    zi*opzi - 2*sqrtxz1*li2spec2* 
    NCi*zi*opzi -  
   4*x*li2spec2*NCi*zi*opzi +  
   4*sqrtxz1*x*li2spec2*NCi*zi* 
    opzi + li2spec19*NCi*zi*opzi -  
   2*x*li2spec19*NCi*zi*opzi +  
   4*sqrtxz1*li2spec21*NCi*zi* 
    opzi - 8*sqrtxz1*x*li2spec21*NCi* 
    zi*opzi + 2*sqrtxz1* 
    li2spec22*NCi*zi*opzi -  
   4*sqrtxz1*x*li2spec22*NCi* 
    zi*opzi + 4*sqrtxz1*li2spec23* 
    NCi*zi*opzi -  
   8*sqrtxz1*x*li2spec23*NCi*zi* 
    opzi - li2spec20*NCi*zi*opzi +  
   2*x*li2spec20*NCi*zi*opzi +  
   2*li2spec3*NCi*zi* 
    opzi + 2*sqrtxz1*li2spec3*NCi*zi*opzi -  
   4*x*li2spec3*NCi*zi* 
    opzi - 4*sqrtxz1*x* 
    li2spec3*NCi*zi* 
    opzi - 2*li2spec4* 
    NCi*zi*opzi +  
   2*sqrtxz1*li2spec4*NCi* 
    zi*opzi + 4*x* 
    li2spec4*NCi*zi* 
    opzi - 4*sqrtxz1*x* 
    li2spec4*NCi*zi* 
    opzi + 8*sqrtxz1*l2*NCi*zi*opzi -  
   16*sqrtxz1*x*l2*NCi*zi*opzi +  
   4*sqrtxz1*lx*NCi*zi*opzi -  
   8*sqrtxz1*x*lx*NCi*zi*opzi +  
   2*l2*lx*NCi*zi*opzi -  
   4*x*l2*lx*NCi*zi*opzi +  
   4*sqrtxz1*l2*lomz*NCi*zi*opzi -  
   8*sqrtxz1*x*l2*lomz*NCi*zi*opzi +  
   2*sqrtxz1*lx*lomz*NCi*zi*opzi -  
   4*sqrtxz1*x*lx*lomz*NCi*zi*opzi -  
   8*sqrtxz1*lspec1*NCi*zi*opzi +  
   16*sqrtxz1*x*lspec1*NCi*zi*opzi -  
   4*l2*lspec1*NCi*zi*opzi -  
   8*sqrtxz1*l2*lspec1*NCi*zi*opzi +  
   8*x*l2*lspec1*NCi*zi*opzi +  
   16*sqrtxz1*x*l2*lspec1*NCi*zi* 
    opzi + 2*sqrtxz1*lx*lspec1*NCi* 
    zi*opzi - 4*sqrtxz1*x*lx*lspec1* 
    NCi*zi*opzi -  
   4*sqrtxz1*lomz*lspec1*NCi*zi* 
    opzi + 8*sqrtxz1*x*lomz*lspec1*NCi* 
    zi*opzi + 4*sqrtxz1*lz*NCi*zi* 
    opzi - 8*sqrtxz1*x*lz*NCi*zi*opzi +  
   6*l2*lz*NCi*zi*opzi +  
   8*sqrtxz1*l2*lz*NCi*zi*opzi -  
   12*x*l2*lz*NCi*zi*opzi -  
   16*sqrtxz1*x*l2*lz*NCi*zi*opzi +  
   lx*lz*NCi*zi*opzi -  
   sqrtxz1*lx*lz*NCi*zi*opzi -  
   2*x*lx*lz*NCi*zi*opzi +  
   2*sqrtxz1*x*lx*lz*NCi*zi*opzi +  
   2*sqrtxz1*lomz*lz*NCi*zi*opzi -  
   4*sqrtxz1*x*lomz*lz*NCi*zi*opzi -  
   2*lspec1*lz*NCi*zi*opzi -  
   4*sqrtxz1*lspec1*lz*NCi*zi*opzi +  
   4*x*lspec1*lz*NCi*zi*opzi +  
   8*sqrtxz1*x*lspec1*lz*NCi*zi* 
    opzi - 4*l2*lspec2*NCi*zi* 
    opzi - 4*sqrtxz1*l2*lspec2*NCi* 
    zi*opzi + 8*x*l2*lspec2*NCi* 
    zi*opzi + 8*sqrtxz1*x*l2*lspec2* 
    NCi*zi*opzi -  
   2*lx*lspec2*NCi*zi*opzi -  
   2*sqrtxz1*lx*lspec2*NCi*zi*opzi +  
   4*x*lx*lspec2*NCi*zi*opzi +  
   4*sqrtxz1*x*lx*lspec2*NCi*zi* 
    opzi + 4*lspec1*lspec2*NCi* 
    zi*opzi + 4*sqrtxz1*lspec1* 
    lspec2*NCi*zi*opzi -  
   8*x*lspec1*lspec2*NCi*zi* 
    opzi - 8*sqrtxz1*x*lspec1*lspec2* 
    NCi*zi*opzi -  
   4*lz*lspec2*NCi*zi*opzi -  
   4*sqrtxz1*lz*lspec2*NCi*zi*opzi +  
   8*x*lz*lspec2*NCi*zi*opzi +  
   8*sqrtxz1*x*lz*lspec2*NCi*zi* 
    opzi - lx*lxpz*NCi*zi*opzi +  
   2*x*lx*lxpz*NCi*zi*opzi +  
   lz*lxpz*NCi*zi*opzi -  
   2*x*lz*lxpz*NCi*zi*opzi +  
   lx*lopxz*NCi*zi*opzi -  
   2*x*lx*lopxz*NCi*zi*opzi +  
   lz*lopxz*NCi*zi*opzi -  
   2*x*lz*lopxz*NCi*zi*opzi -  
   2*sqrtxz1*lomz*lspec7*NCi*zi* 
    opzi + 4*sqrtxz1*x*lomz*lspec7* 
    NCi*zi*opzi +  
   (sqrtxz1*NCi*pi2*zi*opzi)/3. -  
   (2*sqrtxz1*x*NCi*pi2*zi*opzi)/3. -  
   2*li2spec1*NCi*x2*zi* 
    opzi - 2*sqrtxz1*li2spec1*NCi* 
    x2*zi*opzi +  
   2*li2spec2*NCi*x2*zi* 
    opzi - 2*sqrtxz1*li2spec2*NCi* 
    x2*zi*opzi +  
   li2spec19*NCi*x2*zi*opzi +  
   4*sqrtxz1*li2spec21*NCi*x2*zi* 
    opzi + 2*sqrtxz1*li2spec22*NCi* 
    x2*zi*opzi +  
   4*sqrtxz1*li2spec23*NCi*x2*zi* 
    opzi - li2spec20*NCi*x2*zi* 
    opzi + 2*li2spec3* 
    NCi*x2*zi*opzi +  
   2*sqrtxz1*li2spec3*NCi* 
    x2*zi*opzi -  
   2*li2spec4*NCi*x2* 
    zi*opzi + 2*sqrtxz1* 
    li2spec4*NCi*x2* 
    zi*opzi + 8*sqrtxz1*l2*NCi*x2*zi* 
    opzi + 4*sqrtxz1*lx*NCi*x2*zi* 
    opzi + 2*l2*lx*NCi*x2*zi* 
    opzi + 4*sqrtxz1*l2*lomz*NCi*x2*zi* 
    opzi + 2*sqrtxz1*lx*lomz*NCi*x2*zi* 
    opzi - 8*sqrtxz1*lspec1*NCi*x2* 
    zi*opzi - 4*l2*lspec1*NCi* 
    x2*zi*opzi -  
   8*sqrtxz1*l2*lspec1*NCi*x2*zi* 
    opzi + 2*sqrtxz1*lx*lspec1*NCi*x2* 
    zi*opzi - 4*sqrtxz1*lomz*lspec1* 
    NCi*x2*zi*opzi +  
   4*sqrtxz1*lz*NCi*x2*zi*opzi +  
   6*l2*lz*NCi*x2*zi*opzi +  
   8*sqrtxz1*l2*lz*NCi*x2*zi*opzi +  
   lx*lz*NCi*x2*zi*opzi -  
   sqrtxz1*lx*lz*NCi*x2*zi*opzi +  
   2*sqrtxz1*lomz*lz*NCi*x2*zi*opzi -  
   2*lspec1*lz*NCi*x2*zi* 
    opzi - 4*sqrtxz1*lspec1*lz*NCi*x2* 
    zi*opzi - 4*l2*lspec2*NCi* 
    x2*zi*opzi -  
   4*sqrtxz1*l2*lspec2*NCi*x2*zi* 
    opzi - 2*lx*lspec2*NCi*x2* 
    zi*opzi - 2*sqrtxz1*lx*lspec2* 
    NCi*x2*zi*opzi +  
   4*lspec1*lspec2*NCi*x2*zi* 
    opzi + 4*sqrtxz1*lspec1*lspec2* 
    NCi*x2*zi*opzi -  
   4*lz*lspec2*NCi*x2*zi* 
    opzi - 4*sqrtxz1*lz*lspec2*NCi*x2* 
    zi*opzi - lx*lxpz*NCi*x2*zi* 
    opzi + lz*lxpz*NCi*x2*zi* 
    opzi + lx*lopxz*NCi*x2*zi* 
    opzi + lz*lopxz*NCi*x2*zi* 
    opzi - 2*sqrtxz1*lomz*lspec7* 
    NCi*x2*zi*opzi +  
   (sqrtxz1*NCi*pi2*x2*zi*opzi)/3. -  
   2*NC*l22 + 4*NC*x*l22 + 2*NC*z*l22 -  
   4*NC*x*z*l22 + 4*NCi*l22 -  
   8*x*NCi*l22 - 2*z*NCi*l22 +  
   4*x*z*NCi*l22 -  
   12*NCi*sqrtxz1i*l22 +  
   48*x*NCi*sqrtxz1i*l22 +  
   6*z*NCi*sqrtxz1i*l22 -  
   12*x*z*NCi*sqrtxz1i*l22 +  
   4*NC*z*x2*l22 + 4*NCi*x2*l22 -  
   4*z*NCi*x2*l22 -  
   60*NCi*sqrtxz1i*x2*l22 +  
   6*z*NCi*sqrtxz1i*x2*l22 +  
   24*NCi*sqrtxz1i*x3*l22 +  
   2*NC*zi*l22 - 4*NC*x*zi*l22 -  
   4*NCi*zi*l22 -  
   6*sqrtxz1*NCi*zi*l22 +  
   8*x*NCi*zi*l22 +  
   12*sqrtxz1*x*NCi*zi*l22 +  
   4*NC*x2*zi*l22 -  
   4*NCi*x2*zi*l22 -  
   6*sqrtxz1*NCi*x2*zi*l22 +  
   18*NCi*sqrtxz1i*opzi*l22 -  
   60*x*NCi*sqrtxz1i*opzi*l22 -  
   6*z*NCi*sqrtxz1i*opzi*l22 +  
   12*x*z*NCi*sqrtxz1i*opzi*l22 +  
   66*NCi*sqrtxz1i*x2*opzi*l22 -  
   6*z*NCi*sqrtxz1i*x2*opzi*l22 -  
   24*NCi*sqrtxz1i*x3*opzi*l22 +  
   4*NCi*zi*opzi*l22 +  
   6*sqrtxz1*NCi*zi*opzi*l22 -  
   8*x*NCi*zi*opzi*l22 -  
   12*sqrtxz1*x*NCi*zi*opzi*l22 +  
   4*NCi*x2*zi*opzi*l22 +  
   6*sqrtxz1*NCi*x2*zi*opzi*l22 -  
   (NC*lomx2)/2. + NC*x*lomx2 +  
   (NC*z*lomx2)/4. - (NC*x*z*lomx2)/2. +  
   (NCi*lomx2)/2. - x*NCi*lomx2 -  
   (z*NCi*lomx2)/4. +  
   (x*z*NCi*lomx2)/2. - NC*x2*lomx2 +  
   (NC*z*x2*lomx2)/2. +  
   NCi*x2*lomx2 -  
   (z*NCi*x2*lomx2)/2. +  
   (NC*zi*lomx2)/2. - NC*x*zi*lomx2 -  
   (NCi*zi*lomx2)/2. +  
   x*NCi*zi*lomx2 +  
   NC*x2*zi*lomx2 -  
   NCi*x2*zi*lomx2 + NC*lx2 +  
   NC*x*lx2 - NC*z*lx2 + (NCi*lx2)/2. -  
   x*NCi*lx2 + z*NCi*lx2 +  
   3*NCi*sqrtxz1i*lx2 -  
   12*x*NCi*sqrtxz1i*lx2 -  
   (3*z*NCi*sqrtxz1i*lx2)/2. +  
   3*x*z*NCi*sqrtxz1i*lx2 -  
   (3*NC*z*xi2*lx2)/4. +  
   (3*NCi*xi2*lx2)/4. +  
   (3*z*NCi*xi2*lx2)/4. -  
   (3*NC*z*xi*lx2)/2. +  
   (3*NCi*xi*lx2)/2. +  
   (3*z*NCi*xi*lx2)/2. -  
   (3*NC*x2*lx2)/2. + NC*z*x2*lx2 +  
   (3*NCi*x2*lx2)/2. -  
   z*NCi*x2*lx2 +  
   15*NCi*sqrtxz1i*x2*lx2 -  
   (3*z*NCi*sqrtxz1i*x2*lx2)/2. -  
   6*NCi*sqrtxz1i*x3*lx2 -  
   (3*NC*opxi*lx2)/2. -  
   (3*NC*x*opxi*lx2)/2. +  
   3*NC*z*opxi*lx2 +  
   (3*NC*x*z*opxi*lx2)/2. -  
   (3*NCi*opxi*lx2)/2. -  
   3*z*NCi*opxi*lx2 -  
   (3*x*z*NCi*opxi*lx2)/2. +  
   (3*NC*z*xi2*opxi*lx2)/4. -  
   (3*NCi*xi2*opxi*lx2)/4. -  
   (3*z*NCi*xi2*opxi*lx2)/4. +  
   (9*NC*z*xi*opxi*lx2)/4. -  
   (9*NCi*xi*opxi*lx2)/4. -  
   (9*z*NCi*xi*opxi*lx2)/4. -  
   (3*NCi*omzi*lx2)/2. -  
   3*x*NCi*omzi*lx2 -  
   (3*NCi*x2*omzi*lx2)/2. -  
   NC*zi*lx2 - NC*x*zi*lx2 -  
   (NCi*zi*lx2)/2. +  
   (3*sqrtxz1*NCi*zi*lx2)/2. -  
   2*x*NCi*zi*lx2 -  
   3*sqrtxz1*x*NCi*zi*lx2 -  
   (3*NC*xi2*zi*lx2)/4. +  
   (3*NCi*xi2*zi*lx2)/2. -  
   (3*NC*xi*zi*lx2)/2. +  
   3*NCi*xi*zi*lx2 +  
   (3*NC*x2*zi*lx2)/2. -  
   3*NCi*x2*zi*lx2 +  
   (3*sqrtxz1*NCi*x2*zi*lx2)/2. +  
   3*NC*opxi*zi*lx2 +  
   (3*NC*x*opxi*zi*lx2)/2. -  
   (9*NCi*opxi*zi*lx2)/2. -  
   (3*x*NCi*opxi*zi*lx2)/2. +  
   (3*NC*xi2*opxi*zi*lx2)/4. -  
   (3*NCi*xi2*opxi*zi*lx2)/2. +  
   (9*NC*xi*opxi*zi*lx2)/4. -  
   (9*NCi*xi*opxi*zi*lx2)/2. +  
   3*x*NCi*omzi*zi*lx2 -  
   (3*NCi*xi2*omzi*zi*lx2)/2. -  
   3*NCi*xi*omzi*zi*lx2 +  
   (3*NCi*x2*omzi*zi*lx2)/2. +  
   (9*NCi*opxi*omzi*zi*lx2)/2. +  
   (3*x*NCi*opxi*omzi*zi*lx2)/2. +  
   (3*NCi*xi2*opxi*omzi*zi* 
      lx2)/2. + (9*NCi*xi*opxi*omzi* 
      zi*lx2)/2. -  
   (9*NCi*sqrtxz1i*opzi*lx2)/2. +  
   15*x*NCi*sqrtxz1i*opzi*lx2 +  
   (3*z*NCi*sqrtxz1i*opzi*lx2)/2. -  
   3*x*z*NCi*sqrtxz1i*opzi*lx2 -  
   (33*NCi*sqrtxz1i*x2*opzi*lx2)/2. +  
   (3*z*NCi*sqrtxz1i*x2*opzi*lx2)/2. +  
   6*NCi*sqrtxz1i*x3*opzi*lx2 -  
   (3*sqrtxz1*NCi*zi*opzi*lx2)/2. +  
   3*sqrtxz1*x*NCi*zi*opzi*lx2 -  
   (3*sqrtxz1*NCi*x2*zi*opzi*lx2)/2. -  
   (NC*lomz2)/2. + NC*x*lomz2 +  
   (NC*z*lomz2)/4. - (NC*x*z*lomz2)/2. +  
   (NCi*lomz2)/2. - x*NCi*lomz2 -  
   (z*NCi*lomz2)/4. +  
   (x*z*NCi*lomz2)/2. -  
   4*NCi*sqrtxz1i*lomz2 +  
   16*x*NCi*sqrtxz1i*lomz2 +  
   2*z*NCi*sqrtxz1i*lomz2 -  
   4*x*z*NCi*sqrtxz1i*lomz2 -  
   NC*x2*lomz2 + (NC*z*x2*lomz2)/2. +  
   NCi*x2*lomz2 -  
   (z*NCi*x2*lomz2)/2. -  
   20*NCi*sqrtxz1i*x2*lomz2 +  
   2*z*NCi*sqrtxz1i*x2*lomz2 +  
   8*NCi*sqrtxz1i*x3*lomz2 +  
   (NC*zi*lomz2)/2. - NC*x*zi*lomz2 -  
   (NCi*zi*lomz2)/2. -  
   2*sqrtxz1*NCi*zi*lomz2 +  
   x*NCi*zi*lomz2 +  
   4*sqrtxz1*x*NCi*zi*lomz2 +  
   NC*x2*zi*lomz2 -  
   NCi*x2*zi*lomz2 -  
   2*sqrtxz1*NCi*x2*zi*lomz2 +  
   6*NCi*sqrtxz1i*opzi*lomz2 -  
   20*x*NCi*sqrtxz1i*opzi*lomz2 -  
   2*z*NCi*sqrtxz1i*opzi*lomz2 +  
   4*x*z*NCi*sqrtxz1i*opzi*lomz2 +  
   22*NCi*sqrtxz1i*x2*opzi*lomz2 -  
   2*z*NCi*sqrtxz1i*x2*opzi*lomz2 -  
   8*NCi*sqrtxz1i*x3*opzi*lomz2 +  
   2*sqrtxz1*NCi*zi*opzi*lomz2 -  
   4*sqrtxz1*x*NCi*zi*opzi*lomz2 +  
   2*sqrtxz1*NCi*x2*zi*opzi*lomz2 -  
   2*NC*lspec1_2 + 4*NC*x*lspec1_2 +  
   NC*z*lspec1_2 - 2*NC*x*z*lspec1_2 +  
   2*NCi*lspec1_2 -  
   4*x*NCi*lspec1_2 -  
   z*NCi*lspec1_2 +  
   2*x*z*NCi*lspec1_2 -  
   4*NCi*sqrtxz1i*lspec1_2 +  
   16*x*NCi*sqrtxz1i*lspec1_2 +  
   2*z*NCi*sqrtxz1i*lspec1_2 -  
   4*x*z*NCi*sqrtxz1i*lspec1_2 -  
   2*NC*x2*lspec1_2 +  
   2*NC*z*x2*lspec1_2 +  
   2*NCi*x2*lspec1_2 -  
   2*z*NCi*x2*lspec1_2 -  
   20*NCi*sqrtxz1i*x2*lspec1_2 +  
   2*z*NCi*sqrtxz1i*x2*lspec1_2 +  
   8*NCi*sqrtxz1i*x3*lspec1_2 +  
   NC*zi*lspec1_2 -  
   2*NC*x*zi*lspec1_2 -  
   NCi*zi*lspec1_2 -  
   2*sqrtxz1*NCi*zi*lspec1_2 +  
   2*x*NCi*zi*lspec1_2 +  
   4*sqrtxz1*x*NCi*zi*lspec1_2 +  
   2*NC*x2*zi*lspec1_2 -  
   2*NCi*x2*zi*lspec1_2 -  
   2*sqrtxz1*NCi*x2*zi*lspec1_2 +  
   6*NCi*sqrtxz1i*opzi*lspec1_2 -  
   20*x*NCi*sqrtxz1i*opzi* 
    lspec1_2 -  
   2*z*NCi*sqrtxz1i*opzi*lspec1_2 +  
   4*x*z*NCi*sqrtxz1i*opzi* 
    lspec1_2 +  
   22*NCi*sqrtxz1i*x2*opzi* 
    lspec1_2 -  
   2*z*NCi*sqrtxz1i*x2*opzi* 
    lspec1_2 -  
   8*NCi*sqrtxz1i*x3*opzi* 
    lspec1_2 +  
   2*sqrtxz1*NCi*zi*opzi*lspec1_2 -  
   4*sqrtxz1*x*NCi*zi*opzi* 
    lspec1_2 +  
   2*sqrtxz1*NCi*x2*zi*opzi* 
    lspec1_2 - (NC*lz2)/2. -  
   (NC*z*lz2)/4. + (NC*x*z*lz2)/2. +  
   (NCi*lz2)/2. + (z*NCi*lz2)/4. -  
   (x*z*NCi*lz2)/2. -  
   5*NCi*sqrtxz1i*lz2 +  
   20*x*NCi*sqrtxz1i*lz2 +  
   (5*z*NCi*sqrtxz1i*lz2)/2. -  
   5*x*z*NCi*sqrtxz1i*lz2 +  
   (NC*z*xi2*lz2)/4. -  
   (NCi*xi2*lz2)/4. -  
   (z*NCi*xi2*lz2)/4. +  
   (NC*z*xi*lz2)/2. -  
   (NCi*xi*lz2)/2. -  
   (z*NCi*xi*lz2)/2. -  
   (3*NC*z*x2*lz2)/2. + NCi*x2*lz2 +  
   (3*z*NCi*x2*lz2)/2. -  
   25*NCi*sqrtxz1i*x2*lz2 +  
   (5*z*NCi*sqrtxz1i*x2*lz2)/2. +  
   10*NCi*sqrtxz1i*x3*lz2 +  
   (NC*opxi*lz2)/2. +  
   (NC*x*opxi*lz2)/2. - NC*z*opxi*lz2 -  
   (NC*x*z*opxi*lz2)/2. +  
   (NCi*opxi*lz2)/2. +  
   z*NCi*opxi*lz2 +  
   (x*z*NCi*opxi*lz2)/2. -  
   (NC*z*xi2*opxi*lz2)/4. +  
   (NCi*xi2*opxi*lz2)/4. +  
   (z*NCi*xi2*opxi*lz2)/4. -  
   (3*NC*z*xi*opxi*lz2)/4. +  
   (3*NCi*xi*opxi*lz2)/4. +  
   (3*z*NCi*xi*opxi*lz2)/4. +  
   (NCi*omzi*lz2)/2. -  
   x*NCi*omzi*lz2 +  
   (NCi*x2*omzi*lz2)/2. -  
   (NC*zi*lz2)/2. + NC*x*zi*lz2 +  
   (NCi*zi*lz2)/2. -  
   (5*sqrtxz1*NCi*zi*lz2)/2. +  
   5*sqrtxz1*x*NCi*zi*lz2 +  
   (NC*xi2*zi*lz2)/4. -  
   (NCi*xi2*zi*lz2)/2. +  
   (NC*xi*zi*lz2)/2. -  
   NCi*xi*zi*lz2 -  
   2*NC*x2*zi*lz2 +  
   (5*NCi*x2*zi*lz2)/2. -  
   (5*sqrtxz1*NCi*x2*zi*lz2)/2. -  
   NC*opxi*zi*lz2 -  
   (NC*x*opxi*zi*lz2)/2. +  
   (3*NCi*opxi*zi*lz2)/2. +  
   (x*NCi*opxi*zi*lz2)/2. -  
   (NC*xi2*opxi*zi*lz2)/4. +  
   (NCi*xi2*opxi*zi*lz2)/2. -  
   (3*NC*xi*opxi*zi*lz2)/4. +  
   (3*NCi*xi*opxi*zi*lz2)/2. -  
   x*NCi*omzi*zi*lz2 +  
   (NCi*xi2*omzi*zi*lz2)/2. +  
   NCi*xi*omzi*zi*lz2 -  
   (NCi*x2*omzi*zi*lz2)/2. -  
   (3*NCi*opxi*omzi*zi*lz2)/2. -  
   (x*NCi*opxi*omzi*zi*lz2)/2. -  
   (NCi*xi2*opxi*omzi*zi*lz2)/ 
    2. - (3*NCi*xi*opxi*omzi*zi* 
      lz2)/2. + (15*NCi*sqrtxz1i*opzi* 
      lz2)/2. - 25*x*NCi*sqrtxz1i*opzi* 
    lz2 - (5*z*NCi*sqrtxz1i*opzi* 
      lz2)/2. + 5*x*z*NCi*sqrtxz1i*opzi* 
    lz2 + (55*NCi*sqrtxz1i*x2*opzi* 
      lz2)/2. - (5*z*NCi*sqrtxz1i*x2* 
      opzi*lz2)/2. -  
   10*NCi*sqrtxz1i*x3*opzi*lz2 +  
   NCi*zi*opzi*lz2 +  
   (5*sqrtxz1*NCi*zi*opzi*lz2)/2. -  
   2*x*NCi*zi*opzi*lz2 -  
   5*sqrtxz1*x*NCi*zi*opzi*lz2 +  
   NCi*x2*zi*opzi*lz2 +  
   (5*sqrtxz1*NCi*x2*zi*opzi*lz2)/2. -  
   NCi*sqrtxz1i*lspec7_2 +  
   4*x*NCi*sqrtxz1i*lspec7_2 +  
   (z*NCi*sqrtxz1i*lspec7_2)/2. -  
   x*z*NCi*sqrtxz1i*lspec7_2 -  
   5*NCi*sqrtxz1i*x2* 
    lspec7_2 +  
   (z*NCi*sqrtxz1i*x2* 
      lspec7_2)/2. +  
   2*NCi*sqrtxz1i*x3* 
    lspec7_2 -  
   (sqrtxz1*NCi*zi*lspec7_2)/2. +  
   sqrtxz1*x*NCi*zi*lspec7_2 -  
   (sqrtxz1*NCi*x2*zi* 
      lspec7_2)/2. +  
   (3*NCi*sqrtxz1i*opzi* 
      lspec7_2)/2. -  
   5*x*NCi*sqrtxz1i*opzi* 
    lspec7_2 -  
   (z*NCi*sqrtxz1i*opzi* 
      lspec7_2)/2. +  
   x*z*NCi*sqrtxz1i*opzi* 
    lspec7_2 +  
   (11*NCi*sqrtxz1i*x2*opzi* 
      lspec7_2)/2. -  
   (z*NCi*sqrtxz1i*x2*opzi* 
      lspec7_2)/2. -  
   2*NCi*sqrtxz1i*x3*opzi* 
    lspec7_2 +  
   (sqrtxz1*NCi*zi*opzi* 
      lspec7_2)/2. -  
   sqrtxz1*x*NCi*zi*opzi* 
    lspec7_2 +  
   (sqrtxz1*NCi*x2*zi*opzi* 
      lspec7_2)/2.;
};
double C2TG2Q_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double poly2    = 1 + 2*x + x*x - 4*x*z;
  const double poly2i   = 1. / poly2;
  const double poly2i2  = poly2i * poly2i;
  const double sqrtxz2  = sqrt(poly2);
  const double sqrtxz2i = 1. / sqrtxz2;
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double x2 = x * x;
  const double x3 = x * x2;
  const double x4 = x * x3;
  const double x5 = x * x4;
  const double x6 = x * x5;
  const double x7 = x * x6;
  const double xi  = 1. / x;
  const double z2 = z * z;
  const double z3 = z * z2;
  const double zi  = 1. / z;
  const double omxi = 1. / ( 1 - x );
  const double omzi  = 1. / ( 1 - z );
  const double lx  = log(x);
  const double lx2 = lx * lx;
  const double lz  = log(z);
  const double lz2 = lz * lz;
  const double lomx  = log(1 - x);
  const double lomx2 = lomx * lomx;
  const double lomz  = log(1 - z);
  const double lomz2 = lomz * lomz;
  const double lopx  = log(1 + x);
  const double Li2x  = apfel::dilog(x);
  const double Li2mx = apfel::dilog(-x);
  const double Li2z  = apfel::dilog(z);
  const double omxmzi  = 1. / ( 1 - x - z );
  const double omxmzi2 = omxmzi * omxmzi;
  const double xmzi  = 1. / ( x - z );
  const double xmzi2 = xmzi * xmzi;
  const double lomxmz  = log(1 - x - z);
  const double lmopxpz = log(-1 + x + z);
  const double lxmz    = log(x - z);
  const double lmxpz   = log(-x + z);
  const double li2omxzi = apfel::dilog(1 - x*zi);
  const double lspec5   = log(1 - sqrtxz2 + x);
  const double lspec6   = log(1 + sqrtxz2 + x);
  const double li2spec5  = apfel::dilog(0.5 - sqrtxz2/2. - x/2.);
  const double li2spec6  = apfel::dilog(0.5 + sqrtxz2/2. - x/2.);
  const double li2spec7  = apfel::dilog(0.5 - xi/2. - (sqrtxz2*xi)/2.);
  const double li2spec8  = apfel::dilog(0.5 - xi/2. + (sqrtxz2*xi)/2.);
  const double li2spec9  = apfel::dilog((1 - z)*omxi);
  const double li2spec10 = apfel::dilog(x*(1 - z)*omxi*zi);
  const double li2spec11 = apfel::dilog((1 - x)*omzi);
  const double li2spec12 = apfel::dilog((1 - x)*z*xi*omzi);
  const double li2spec13 = apfel::dilog(z*omxi);
  const double li2spec14 = apfel::dilog((1 - z)*z*omxi*xi);
  const double li2spec15 = apfel::dilog(x*z*omxi*omzi);
  const double li2spec16 = apfel::dilog((1 - x)*zi);
  const double li2spec17 = apfel::dilog((1 - x)*(1 - z)*xi*zi);
  const double li2spec18 = apfel::dilog((1 - x)*x*omzi*zi);
  const double Tu1 = (z < 1 - x && z < x ? 1 : 0);
  const double Tu2 = (z > 1 - x && z < x ? 1 : 0);
  const double Tu3 = (z < 1 - x && z > x ? 1 : 0);
  const double Tu4 = (z > 1 - x && z > x ? 1 : 0);
  return (71*NC)/48. - (155*NC*x)/48. + (153*NC*z)/16. - 11*NC*x*z -  
   NC*Li2mx - 2*NC*x*Li2mx + NC*z*Li2mx + 2*NC*x*z*Li2mx +  
   (5*NC*Li2x)/2. + 8*NC*x*Li2x + (NC*z*Li2x)/2. - 2*NC*x*z*Li2x -  
   2*NC*Li2z + 4*NC*x*Li2z + (3*NC*li2omxzi)/2. -  
   3*NC*x*li2omxzi - (5*NC*lomx)/4. -  
   (27*NC*x*lomx)/2. + (3*NC*z*lomx)/8. -  
   (NC*x*z*lomx)/4. + (NC*lx)/16. + (19*NC*x*lx)/4. +  
   (37*NC*z*lx)/8. + (39*NC*x*z*lx)/8. +  
   (11*NC*lomx*lx)/2. - 11*NC*x*lomx*lx +  
   NC*z*lomx*lx - 2*NC*x*z*lomx*lx -  
   NC*lx*lopx - 2*NC*x*lx*lopx + NC*z*lx*lopx +  
   2*NC*x*z*lx*lopx - (7*NC*lomz)/2. - 14*NC*x*lomz +  
   (3*NC*z*lomz)/8. - (NC*x*z*lomz)/4. -  
   (9*NC*lomx*lomz)/2. + 9*NC*x*lomx*lomz -  
   (NC*z*lomx*lomz)/2. + NC*x*z*lomx*lomz +  
   4*NC*lx*lomz - 21*NC*x*lx*lomz +  
   (NC*z*lx*lomz)/2. - (63*NC*lz)/16. - (91*NC*x*lz)/8. +  
   (23*NC*z*lz)/4. - (57*NC*x*z*lz)/8. -  
   (5*NC*lomx*lz)/2. + 5*NC*x*lomx*lz +  
   (9*NC*lx*lz)/2. - 11*NC*x*lx*lz - NC*x*z*lx*lz -  
   (11*NC*lomz*lz)/2. + 11*NC*x*lomz*lz -  
   (21*NCi)/16. - (7*x*NCi)/16. + (7*z*NCi)/16. +  
   x*z*NCi + 2*Li2mx*NCi + 4*x*Li2mx*NCi -  
   z*Li2mx*NCi - 2*x*z*Li2mx*NCi -  
   (3*Li2x*NCi)/2. + 2*x*Li2x*NCi -  
   (z*Li2x*NCi)/2. + 2*x*z*Li2x*NCi + 2*Li2z*NCi -  
   4*x*Li2z*NCi - (li2omxzi*NCi)/2. +  
   x*li2omxzi*NCi + (7*lomx*NCi)/4. +  
   (3*x*lomx*NCi)/2. - (3*z*lomx*NCi)/8. +  
   (x*z*lomx*NCi)/4. - (53*lx*NCi)/16. -  
   (x*lx*NCi)/2. + (3*z*lx*NCi)/8. +  
   (x*z*lx*NCi)/8. - (13*lomx*lx*NCi)/2. +  
   13*x*lomx*lx*NCi - z*lomx*lx*NCi +  
   2*x*z*lomx*lx*NCi + 2*lx*lopx*NCi +  
   4*x*lx*lopx*NCi - z*lx*lopx*NCi -  
   2*x*z*lx*lopx*NCi + (lomz*NCi)/2. +  
   3*x*lomz*NCi - (3*z*lomz*NCi)/8. +  
   (x*z*lomz*NCi)/4. + 4*lomx*lomz*NCi -  
   8*x*lomx*lomz*NCi +  
   (z*lomx*lomz*NCi)/2. -  
   x*z*lomx*lomz*NCi - 5*lx*lomz*NCi +  
   11*x*lx*lomz*NCi - (z*lx*lomz*NCi)/2. +  
   (27*lz*NCi)/16. - (11*x*lz*NCi)/8. -  
   (3*z*lz*NCi)/4. + (17*x*z*lz*NCi)/8. +  
   2*lomx*lz*NCi - 4*x*lomx*lz*NCi -  
   (9*lx*lz*NCi)/2. + 5*x*lx*lz*NCi +  
   x*z*lx*lz*NCi + (7*lomz*lz*NCi)/2. -  
   7*x*lomz*lz*NCi + (NC*pi2)/4. - 3*NC*x*pi2 +  
   (NC*z*pi2)/12. + (NC*x*z*pi2)/3. -  
   (5*NCi*pi2)/12. + (5*x*NCi*pi2)/3. -  
   (z*NCi*pi2)/12. - (x*z*NCi*pi2)/3. -  
   (3*NC*lx*poly2i2)/32. - (3*NC*x*lx*poly2i2)/32. +  
   (3*NC*lz*poly2i2)/32. - (3*NC*x*lz*poly2i2)/32. +  
   (3*lx*NCi*poly2i2)/32. +  
   (3*x*lx*NCi*poly2i2)/32. -  
   (3*lz*NCi*poly2i2)/32. +  
   (3*x*lz*NCi*poly2i2)/32. - (NC*poly2i)/16. -  
   (5*NC*lx*poly2i)/32. + (23*NC*x*lx*poly2i)/16. +  
   (5*NC*lz*poly2i)/32. + (23*NC*x*lz*poly2i)/16. +  
   (NCi*poly2i)/16. - (3*lx*NCi*poly2i)/32. -  
   (15*x*lx*NCi*poly2i)/16. +  
   (3*lz*NCi*poly2i)/32. -  
   (15*x*lz*NCi*poly2i)/16. +  
   (9*NC*li2spec5*sqrtxz2i)/4. +  
   (211*NC*x*li2spec5*sqrtxz2i)/32. -  
   (9*NC*z*li2spec5*sqrtxz2i)/2. -  
   (51*NC*x*z*li2spec5*sqrtxz2i)/4. -  
   (9*NC*li2spec6*sqrtxz2i)/4. -  
   (211*NC*x*li2spec6*sqrtxz2i)/32. +  
   (9*NC*z*li2spec6*sqrtxz2i)/2. +  
   (51*NC*x*z*li2spec6*sqrtxz2i)/4. -  
   (9*NC*li2spec7*sqrtxz2i)/ 
    4. - (211*NC*x*li2spec7* 
      sqrtxz2i)/32. + (9*NC*z* 
      li2spec7*sqrtxz2i)/2. +  
   (51*NC*x*z*li2spec7* 
      sqrtxz2i)/4. + (9*NC* 
      li2spec8*sqrtxz2i)/4. +  
   (211*NC*x*li2spec8* 
      sqrtxz2i)/32. - (9*NC*z* 
      li2spec8*sqrtxz2i)/2. -  
   (51*NC*x*z*li2spec8* 
      sqrtxz2i)/4. - (9*NC*lx*lspec5* 
      sqrtxz2i)/4. - (211*NC*x*lx*lspec5* 
      sqrtxz2i)/32. + (9*NC*z*lx*lspec5* 
      sqrtxz2i)/2. + (51*NC*x*z*lx*lspec5* 
      sqrtxz2i)/4. + (9*NC*lx*lspec6* 
      sqrtxz2i)/4. + (211*NC*x*lx*lspec6* 
      sqrtxz2i)/32. - (9*NC*z*lx*lspec6* 
      sqrtxz2i)/2. - (51*NC*x*z*lx*lspec6* 
      sqrtxz2i)/4. - (li2spec5*NCi* 
      sqrtxz2i)/2. - (259*x*li2spec5*NCi* 
      sqrtxz2i)/32. + z*li2spec5*NCi* 
    sqrtxz2i + (11*x*z*li2spec5*NCi* 
      sqrtxz2i)/4. + (li2spec6*NCi* 
      sqrtxz2i)/2. + (259*x*li2spec6*NCi* 
      sqrtxz2i)/32. - z*li2spec6*NCi* 
    sqrtxz2i - (11*x*z*li2spec6*NCi* 
      sqrtxz2i)/4. + (li2spec7* 
      NCi*sqrtxz2i)/2. +  
   (259*x*li2spec7*NCi* 
      sqrtxz2i)/32. - z*li2spec7*NCi*sqrtxz2i -  
   (11*x*z*li2spec7*NCi* 
      sqrtxz2i)/4. - (li2spec8* 
      NCi*sqrtxz2i)/2. -  
   (259*x*li2spec8*NCi* 
      sqrtxz2i)/32. + z*li2spec8*NCi*sqrtxz2i +  
   (11*x*z*li2spec8*NCi* 
      sqrtxz2i)/4. + (lx*lspec5*NCi* 
      sqrtxz2i)/2. + (259*x*lx*lspec5*NCi* 
      sqrtxz2i)/32. - z*lx*lspec5*NCi* 
    sqrtxz2i - (11*x*z*lx*lspec5*NCi* 
      sqrtxz2i)/4. - (lx*lspec6*NCi* 
      sqrtxz2i)/2. - (259*x*lx*lspec6*NCi* 
      sqrtxz2i)/32. + z*lx*lspec6*NCi* 
    sqrtxz2i + (11*x*z*lx*lspec6*NCi* 
      sqrtxz2i)/4. + (3*NC*x*li2spec5* 
      poly2i2*sqrtxz2i)/32. -  
   (3*NC*x*li2spec6*poly2i2*sqrtxz2i)/32. -  
   (3*NC*x*li2spec7*poly2i2* 
      sqrtxz2i)/32. + (3*NC*x* 
      li2spec8*poly2i2* 
      sqrtxz2i)/32. - (3*NC*x*lx*lspec5* 
      poly2i2*sqrtxz2i)/32. +  
   (3*NC*x*lx*lspec6*poly2i2*sqrtxz2i)/32. -  
   (3*x*li2spec5*NCi*poly2i2* 
      sqrtxz2i)/32. + (3*x*li2spec6*NCi* 
      poly2i2*sqrtxz2i)/32. +  
   (3*x*li2spec7*NCi* 
      poly2i2*sqrtxz2i)/32. -  
   (3*x*li2spec8*NCi* 
      poly2i2*sqrtxz2i)/32. +  
   (3*x*lx*lspec5*NCi*poly2i2*sqrtxz2i)/ 
    32. - (3*x*lx*lspec6*NCi*poly2i2* 
      sqrtxz2i)/32. - (5*NC*x*li2spec5* 
      poly2i*sqrtxz2i)/8. +  
   (5*NC*x*li2spec6*poly2i*sqrtxz2i)/8. +  
   (5*NC*x*li2spec7*poly2i* 
      sqrtxz2i)/8. - (5*NC*x* 
      li2spec8*poly2i* 
      sqrtxz2i)/8. + (5*NC*x*lx*lspec5*poly2i* 
      sqrtxz2i)/8. - (5*NC*x*lx*lspec6*poly2i* 
      sqrtxz2i)/8. + (x*li2spec5*NCi* 
      poly2i*sqrtxz2i)/2. -  
   (x*li2spec6*NCi*poly2i*sqrtxz2i)/ 
    2. - (x*li2spec7*NCi* 
      poly2i*sqrtxz2i)/2. +  
   (x*li2spec8*NCi* 
      poly2i*sqrtxz2i)/2. -  
   (x*lx*lspec5*NCi*poly2i*sqrtxz2i)/ 
    2. + (x*lx*lspec6*NCi*poly2i* 
      sqrtxz2i)/2. - (425*NC*xi)/144. -  
   (4*NC*lomx*xi)/3. - (3*NC*lx*xi)/16. -  
   (4*NC*lomz*xi)/3. - (73*NC*lz*xi)/48. +  
   (NCi*xi)/16. - (lx*NCi*xi)/16. -  
   (lz*NCi*xi)/16. +  
   (3*NC*lx*poly2i2*xi)/32. +  
   (3*NC*lz*poly2i2*xi)/32. -  
   (3*lx*NCi*poly2i2*xi)/32. -  
   (3*lz*NCi*poly2i2*xi)/32. +  
   (NC*poly2i*xi)/16. +  
   (3*NC*lx*poly2i*xi)/32. +  
   (3*NC*lz*poly2i*xi)/32. -  
   (NCi*poly2i*xi)/16. +  
   (5*lx*NCi*poly2i*xi)/32. +  
   (5*lz*NCi*poly2i*xi)/32. +  
   (5*NC*li2spec5*sqrtxz2i*xi)/64. -  
   (5*NC*li2spec6*sqrtxz2i*xi)/64. -  
   (5*NC*li2spec7*sqrtxz2i* 
      xi)/64. + (5*NC*li2spec8* 
      sqrtxz2i*xi)/64. -  
   (5*NC*lx*lspec5*sqrtxz2i*xi)/64. +  
   (5*NC*lx*lspec6*sqrtxz2i*xi)/64. +  
   (3*li2spec5*NCi*sqrtxz2i*xi)/ 
    64. - (3*li2spec6*NCi*sqrtxz2i* 
      xi)/64. - (3*li2spec7* 
      NCi*sqrtxz2i*xi)/64. +  
   (3*li2spec8*NCi* 
      sqrtxz2i*xi)/64. -  
   (3*lx*lspec5*NCi*sqrtxz2i*xi)/64. +  
   (3*lx*lspec6*NCi*sqrtxz2i*xi)/64. -  
   (3*NC*li2spec5*poly2i2*sqrtxz2i* 
      xi)/64. + (3*NC*li2spec6*poly2i2* 
      sqrtxz2i*xi)/64. +  
   (3*NC*li2spec7*poly2i2* 
      sqrtxz2i*xi)/64. -  
   (3*NC*li2spec8*poly2i2* 
      sqrtxz2i*xi)/64. +  
   (3*NC*lx*lspec5*poly2i2*sqrtxz2i*xi)/ 
    64. - (3*NC*lx*lspec6*poly2i2*sqrtxz2i* 
      xi)/64. + (3*li2spec5*NCi* 
      poly2i2*sqrtxz2i*xi)/64. -  
   (3*li2spec6*NCi*poly2i2*sqrtxz2i* 
      xi)/64. - (3*li2spec7* 
      NCi*poly2i2*sqrtxz2i*xi)/64. +  
   (3*li2spec8*NCi* 
      poly2i2*sqrtxz2i*xi)/64. -  
   (3*lx*lspec5*NCi*poly2i2*sqrtxz2i* 
      xi)/64. + (3*lx*lspec6*NCi* 
      poly2i2*sqrtxz2i*xi)/64. -  
   (NC*li2spec5*poly2i*sqrtxz2i*xi)/ 
    32. + (NC*li2spec6*poly2i*sqrtxz2i* 
      xi)/32. + (NC*li2spec7* 
      poly2i*sqrtxz2i*xi)/32. -  
   (NC*li2spec8*poly2i* 
      sqrtxz2i*xi)/32. +  
   (NC*lx*lspec5*poly2i*sqrtxz2i*xi)/ 
    32. - (NC*lx*lspec6*poly2i*sqrtxz2i* 
      xi)/32. - (3*li2spec5*NCi* 
      poly2i*sqrtxz2i*xi)/32. +  
   (3*li2spec6*NCi*poly2i*sqrtxz2i* 
      xi)/32. + (3*li2spec7* 
      NCi*poly2i*sqrtxz2i*xi)/32. -  
   (3*li2spec8*NCi* 
      poly2i*sqrtxz2i*xi)/32. +  
   (3*lx*lspec5*NCi*poly2i*sqrtxz2i* 
      xi)/32. - (3*lx*lspec6*NCi* 
      poly2i*sqrtxz2i*xi)/32. +  
   (1127*NC*x2)/144. + (3*NC*z*x2)/2. + 2*NC*z*Li2mx*x2 -  
   NC*Li2x*x2 + NC*z*Li2x*x2 - 4*NC*Li2z*x2 +  
   3*NC*li2omxzi*x2 + (193*NC*lomx*x2)/12. -  
   (NC*z*lomx*x2)/2. - (129*NC*lx*x2)/4. +  
   (NC*z*lx*x2)/2. + 10*NC*lomx*lx*x2 +  
   2*NC*z*lomx*lx*x2 + 2*NC*z*lx*lopx*x2 +  
   (199*NC*lomz*x2)/12. - (NC*z*lomz*x2)/2. -  
   9*NC*lomx*lomz*x2 -  
   NC*z*lomx*lomz*x2 + 13*NC*lx*lomz*x2 +  
   NC*z*lx*lomz*x2 + (77*NC*lz*x2)/6. +  
   NC*z*lz*x2 - 5*NC*lomx*lz*x2 +  
   9*NC*lx*lz*x2 - 11*NC*lomz*lz*x2 -  
   (7*NCi*x2)/16. - (3*z*NCi*x2)/2. +  
   2*Li2mx*NCi*x2 - 2*z*Li2mx*NCi*x2 -  
   Li2x*NCi*x2 - z*Li2x*NCi*x2 +  
   4*Li2z*NCi*x2 - li2omxzi*NCi*x2 -  
   (11*lomx*NCi*x2)/4. +  
   (z*lomx*NCi*x2)/2. + (lx*NCi*x2)/2. -  
   (z*lx*NCi*x2)/2. -  
   12*lomx*lx*NCi*x2 -  
   2*z*lomx*lx*NCi*x2 +  
   2*lx*lopx*NCi*x2 -  
   2*z*lx*lopx*NCi*x2 -  
   (13*lomz*NCi*x2)/4. +  
   (z*lomz*NCi*x2)/2. +  
   8*lomx*lomz*NCi*x2 +  
   z*lomx*lomz*NCi*x2 -  
   11*lx*lomz*NCi*x2 -  
   z*lx*lomz*NCi*x2 +  
   (11*lz*NCi*x2)/4. - z*lz*NCi*x2 +  
   4*lomx*lz*NCi*x2 -  
   7*lx*lz*NCi*x2 +  
   7*lomz*lz*NCi*x2 + (5*NC*pi2*x2)/3. +  
   (NC*z*pi2*x2)/6. - (4*NCi*pi2*x2)/3. -  
   (z*NCi*pi2*x2)/6. +  
   (3*NC*lx*poly2i2*x2)/32. -  
   (3*NC*lz*poly2i2*x2)/32. -  
   (3*lx*NCi*poly2i2*x2)/32. +  
   (3*lz*NCi*poly2i2*x2)/32. -  
   (23*NC*lx*poly2i*x2)/16. +  
   (23*NC*lz*poly2i*x2)/16. +  
   (15*lx*NCi*poly2i*x2)/16. -  
   (15*lz*NCi*poly2i*x2)/16. +  
   (129*NC*li2spec5*sqrtxz2i*x2)/16. -  
   (129*NC*z*li2spec5*sqrtxz2i*x2)/8. -  
   (129*NC*li2spec6*sqrtxz2i*x2)/16. +  
   (129*NC*z*li2spec6*sqrtxz2i*x2)/8. -  
   (129*NC*li2spec7*sqrtxz2i* 
      x2)/16. + (129*NC*z* 
      li2spec7*sqrtxz2i* 
      x2)/8. + (129*NC*li2spec8* 
      sqrtxz2i*x2)/16. -  
   (129*NC*z*li2spec8*sqrtxz2i* 
      x2)/8. - (129*NC*lx*lspec5*sqrtxz2i* 
      x2)/16. + (129*NC*z*lx*lspec5*sqrtxz2i* 
      x2)/8. + (129*NC*lx*lspec6*sqrtxz2i* 
      x2)/16. - (129*NC*z*lx*lspec6*sqrtxz2i* 
      x2)/8. - (29*li2spec5*NCi* 
      sqrtxz2i*x2)/16. +  
   (29*z*li2spec5*NCi*sqrtxz2i*x2)/ 
    8. + (29*li2spec6*NCi*sqrtxz2i* 
      x2)/16. - (29*z*li2spec6*NCi* 
      sqrtxz2i*x2)/8. +  
   (29*li2spec7*NCi* 
      sqrtxz2i*x2)/16. -  
   (29*z*li2spec7*NCi* 
      sqrtxz2i*x2)/8. -  
   (29*li2spec8*NCi* 
      sqrtxz2i*x2)/16. +  
   (29*z*li2spec8*NCi* 
      sqrtxz2i*x2)/8. +  
   (29*lx*lspec5*NCi*sqrtxz2i*x2)/16. -  
   (29*z*lx*lspec5*NCi*sqrtxz2i*x2)/ 
    8. - (29*lx*lspec6*NCi*sqrtxz2i*x2)/ 
    16. + (29*z*lx*lspec6*NCi*sqrtxz2i* 
      x2)/8. - (3*NC*lx*poly2i2*x3)/32. -  
   (3*NC*lz*poly2i2*x3)/32. +  
   (3*lx*NCi*poly2i2*x3)/32. +  
   (3*lz*NCi*poly2i2*x3)/32. -  
   (NC*poly2i*x3)/16. - (NC*lx*poly2i*x3)/32. -  
   (NC*lz*poly2i*x3)/32. +  
   (NCi*poly2i*x3)/16. +  
   (9*lx*NCi*poly2i*x3)/32. +  
   (9*lz*NCi*poly2i*x3)/32. +  
   (193*NC*li2spec5*sqrtxz2i*x3)/64. -  
   (193*NC*li2spec6*sqrtxz2i*x3)/64. -  
   (193*NC*li2spec7*sqrtxz2i* 
      x3)/64. + (193*NC*li2spec8*sqrtxz2i*x3)/64. -  
   (193*NC*lx*lspec5*sqrtxz2i*x3)/64. +  
   (193*NC*lx*lspec6*sqrtxz2i*x3)/64. -  
   (73*li2spec5*NCi*sqrtxz2i*x3)/ 
    64. + (73*li2spec6*NCi*sqrtxz2i* 
      x3)/64. + (73*li2spec7* 
      NCi*sqrtxz2i*x3)/64. -  
   (73*li2spec8*NCi* 
      sqrtxz2i*x3)/64. +  
   (73*lx*lspec5*NCi*sqrtxz2i*x3)/64. -  
   (73*lx*lspec6*NCi*sqrtxz2i*x3)/64. +  
   (23*NC*li2spec5*poly2i*sqrtxz2i* 
      x3)/32. - (23*NC*li2spec6*poly2i* 
      sqrtxz2i*x3)/32. -  
   (23*NC*li2spec7*poly2i* 
      sqrtxz2i*x3)/32. +  
   (23*NC*li2spec8*poly2i* 
      sqrtxz2i*x3)/32. -  
   (23*NC*lx*lspec5*poly2i*sqrtxz2i*x3)/ 
    32. + (23*NC*lx*lspec6*poly2i*sqrtxz2i* 
      x3)/32. - (19*li2spec5*NCi* 
      poly2i*sqrtxz2i*x3)/32. +  
   (19*li2spec6*NCi*poly2i*sqrtxz2i* 
      x3)/32. + (19*li2spec7* 
      NCi*poly2i*sqrtxz2i*x3)/32. -  
   (19*li2spec8*NCi* 
      poly2i*sqrtxz2i*x3)/32. +  
   (19*lx*lspec5*NCi*poly2i*sqrtxz2i* 
      x3)/32. - (19*lx*lspec6*NCi* 
      poly2i*sqrtxz2i*x3)/32. +  
   (3*NC*lx*poly2i2*x4)/32. -  
   (3*NC*lz*poly2i2*x4)/32. -  
   (3*lx*NCi*poly2i2*x4)/32. +  
   (3*lz*NCi*poly2i2*x4)/32. +  
   (NC*poly2i*x4)/16. +  
   (3*NC*lx*poly2i*x4)/32. -  
   (3*NC*lz*poly2i*x4)/32. -  
   (NCi*poly2i*x4)/16. -  
   (11*lx*NCi*poly2i*x4)/32. +  
   (11*lz*NCi*poly2i*x4)/32. +  
   (3*NC*lx*poly2i2*x5)/32. +  
   (3*NC*lz*poly2i2*x5)/32. -  
   (3*lx*NCi*poly2i2*x5)/32. -  
   (3*lz*NCi*poly2i2*x5)/32. -  
   (3*NC*li2spec5*poly2i2*sqrtxz2i*x5)/ 
    32. + (3*NC*li2spec6*poly2i2*sqrtxz2i* 
      x5)/32. + (3*NC*li2spec7* 
      poly2i2*sqrtxz2i*x5)/32. -  
   (3*NC*li2spec8*poly2i2* 
      sqrtxz2i*x5)/32. +  
   (3*NC*lx*lspec5*poly2i2*sqrtxz2i*x5)/ 
    32. - (3*NC*lx*lspec6*poly2i2*sqrtxz2i* 
      x5)/32. + (3*li2spec5*NCi* 
      poly2i2*sqrtxz2i*x5)/32. -  
   (3*li2spec6*NCi*poly2i2*sqrtxz2i* 
      x5)/32. - (3*li2spec7* 
      NCi*poly2i2*sqrtxz2i*x5)/32. +  
   (3*li2spec8*NCi* 
      poly2i2*sqrtxz2i*x5)/32. -  
   (3*lx*lspec5*NCi*poly2i2*sqrtxz2i* 
      x5)/32. + (3*lx*lspec6*NCi*poly2i2* 
      sqrtxz2i*x5)/32. -  
   (NC*li2spec5*poly2i*sqrtxz2i*x5)/ 
    16. + (NC*li2spec6*poly2i*sqrtxz2i* 
      x5)/16. + (NC*li2spec7* 
      poly2i*sqrtxz2i*x5)/16. -  
   (NC*li2spec8*poly2i* 
      sqrtxz2i*x5)/16. +  
   (NC*lx*lspec5*poly2i*sqrtxz2i*x5)/ 
    16. - (NC*lx*lspec6*poly2i*sqrtxz2i* 
      x5)/16. + (3*li2spec5*NCi* 
      poly2i*sqrtxz2i*x5)/16. -  
   (3*li2spec6*NCi*poly2i*sqrtxz2i* 
      x5)/16. - (3*li2spec7* 
      NCi*poly2i*sqrtxz2i*x5)/16. +  
   (3*li2spec8*NCi* 
      poly2i*sqrtxz2i*x5)/16. -  
   (3*lx*lspec5*NCi*poly2i*sqrtxz2i* 
      x5)/16. + (3*lx*lspec6*NCi*poly2i* 
      sqrtxz2i*x5)/16. -  
   (3*NC*lx*poly2i2*x6)/32. +  
   (3*NC*lz*poly2i2*x6)/32. +  
   (3*lx*NCi*poly2i2*x6)/32. -  
   (3*lz*NCi*poly2i2*x6)/32. +  
   (3*NC*li2spec5*poly2i2*sqrtxz2i*x7)/ 
    64. - (3*NC*li2spec6*poly2i2*sqrtxz2i* 
      x7)/64. - (3*NC*li2spec7* 
      poly2i2*sqrtxz2i*x7)/64. +  
   (3*NC*li2spec8*poly2i2* 
      sqrtxz2i*x7)/64. -  
   (3*NC*lx*lspec5*poly2i2*sqrtxz2i*x7)/ 
    64. + (3*NC*lx*lspec6*poly2i2*sqrtxz2i* 
      x7)/64. - (3*li2spec5*NCi* 
      poly2i2*sqrtxz2i*x7)/64. +  
   (3*li2spec6*NCi*poly2i2*sqrtxz2i* 
      x7)/64. + (3*li2spec7* 
      NCi*poly2i2*sqrtxz2i*x7)/64. -  
   (3*li2spec8*NCi* 
      poly2i2*sqrtxz2i*x7)/64. +  
   (3*lx*lspec5*NCi*poly2i2*sqrtxz2i* 
      x7)/64. - (3*lx*lspec6*NCi*poly2i2* 
      sqrtxz2i*x7)/64. - 2*NC*omzi +  
   4*NC*x*omzi - 2*NC*Li2x*omzi +  
   4*NC*x*Li2x*omzi + (NC*Li2z*omzi)/2. -  
   NC*x*Li2z*omzi + (NC*li2omxzi*omzi)/2. -  
   NC*x*li2omxzi*omzi +  
   (NC*lomx*omzi)/2. + (13*NC*lx*omzi)/2. -  
   16*NC*x*lx*omzi - (7*NC*lomx*lx*omzi)/2. +  
   7*NC*x*lomx*lx*omzi + NC*lomz*omzi +  
   NC*lomx*lomz*omzi -  
   2*NC*x*lomx*lomz*omzi -  
   3*NC*lx*lomz*omzi +  
   6*NC*x*lx*lomz*omzi - 4*NC*lz*omzi +  
   13*NC*x*lz*omzi + (3*NC*lomx*lz*omzi)/2. -  
   3*NC*x*lomx*lz*omzi +  
   (NC*lx*lz*omzi)/2. + NC*x*lx*lz*omzi +  
   (5*NC*lomz*lz*omzi)/2. -  
   5*NC*x*lomz*lz*omzi + (NCi*omzi)/2. -  
   4*x*NCi*omzi + 2*Li2x*NCi*omzi -  
   4*x*Li2x*NCi*omzi - Li2z*NCi*omzi +  
   2*x*Li2z*NCi*omzi -  
   (li2omxzi*NCi*omzi)/2. +  
   x*li2omxzi*NCi*omzi -  
   lomx*NCi*omzi -  
   (13*lx*NCi*omzi)/2. +  
   16*x*lx*NCi*omzi +  
   (19*lomx*lx*NCi*omzi)/2. -  
   19*x*lomx*lx*NCi*omzi -  
   lomz*NCi*omzi -  
   4*lomx*lomz*NCi*omzi +  
   8*x*lomx*lomz*NCi*omzi +  
   6*lx*lomz*NCi*omzi -  
   12*x*lx*lomz*NCi*omzi +  
   (7*lz*NCi*omzi)/2. -  
   9*x*lz*NCi*omzi -  
   4*lomx*lz*NCi*omzi +  
   8*x*lomx*lz*NCi*omzi +  
   4*lx*lz*NCi*omzi -  
   8*x*lx*lz*NCi*omzi -  
   4*lomz*lz*NCi*omzi +  
   8*x*lomz*lz*NCi*omzi +  
   (NC*pi2*omzi)/6. - (NC*x*pi2*omzi)/3. +  
   (3*NCi*pi2*omzi)/4. -  
   (3*x*NCi*pi2*omzi)/2. +  
   (NC*li2spec5*sqrtxz2i*omzi)/2. -  
   (3*NC*x*li2spec5*sqrtxz2i*omzi)/2. -  
   (NC*li2spec6*sqrtxz2i*omzi)/2. +  
   (3*NC*x*li2spec6*sqrtxz2i*omzi)/2. -  
   (NC*li2spec7*sqrtxz2i* 
      omzi)/2. + (3*NC*x* 
      li2spec7*sqrtxz2i* 
      omzi)/2. + (NC*li2spec8*sqrtxz2i*omzi)/2. -  
   (3*NC*x*li2spec8*sqrtxz2i* 
      omzi)/2. - (NC*lx*lspec5*sqrtxz2i* 
      omzi)/2. + (3*NC*x*lx*lspec5*sqrtxz2i* 
      omzi)/2. + (NC*lx*lspec6*sqrtxz2i* 
      omzi)/2. - (3*NC*x*lx*lspec6*sqrtxz2i* 
      omzi)/2. - li2spec5*NCi* 
    sqrtxz2i*omzi +  
   3*x*li2spec5*NCi*sqrtxz2i* 
    omzi + li2spec6*NCi*sqrtxz2i* 
    omzi - 3*x*li2spec6*NCi* 
    sqrtxz2i*omzi +  
   li2spec7*NCi* 
    sqrtxz2i*omzi -  
   3*x*li2spec7*NCi* 
    sqrtxz2i*omzi -  
   li2spec8*NCi* 
    sqrtxz2i*omzi +  
   3*x*li2spec8*NCi* 
    sqrtxz2i*omzi +  
   lx*lspec5*NCi*sqrtxz2i*omzi -  
   3*x*lx*lspec5*NCi*sqrtxz2i*omzi -  
   lx*lspec6*NCi*sqrtxz2i*omzi +  
   3*x*lx*lspec6*NCi*sqrtxz2i*omzi +  
   (2*NC*lz*xi*omzi)/3. - 4*NC*x2*omzi -  
   4*NC*Li2x*x2*omzi + NC*Li2z*x2*omzi +  
   NC*li2omxzi*x2*omzi +  
   16*NC*lx*x2*omzi -  
   7*NC*lomx*lx*x2*omzi +  
   2*NC*lomx*lomz*x2*omzi -  
   6*NC*lx*lomz*x2*omzi -  
   (38*NC*lz*x2*omzi)/3. +  
   3*NC*lomx*lz*x2*omzi +  
   NC*lx*lz*x2*omzi +  
   5*NC*lomz*lz*x2*omzi +  
   3*NCi*x2*omzi +  
   2*Li2x*NCi*x2*omzi -  
   Li2z*NCi*x2*omzi -  
   li2omxzi*NCi*x2*omzi -  
   2*lomx*NCi*x2*omzi -  
   5*lx*NCi*x2*omzi +  
   12*lomx*lx*NCi*x2*omzi -  
   lomz*NCi*x2*omzi -  
   6*lomx*lomz*NCi*x2*omzi +  
   9*lx*lomz*NCi*x2*omzi +  
   3*lz*NCi*x2*omzi -  
   5*lomx*lz*NCi*x2*omzi +  
   6*lx*lz*NCi*x2*omzi -  
   5*lomz*lz*NCi*x2*omzi +  
   (NC*pi2*x2*omzi)/3. +  
   NCi*pi2*x2*omzi +  
   2*NC*li2spec5*sqrtxz2i*x2*omzi -  
   2*NC*li2spec6*sqrtxz2i*x2*omzi -  
   2*NC*li2spec7*sqrtxz2i* 
    x2*omzi + 2*NC* 
    li2spec8*sqrtxz2i*x2* 
    omzi - 2*NC*lx*lspec5*sqrtxz2i*x2* 
    omzi + 2*NC*lx*lspec6*sqrtxz2i*x2* 
    omzi - 3*li2spec5*NCi*sqrtxz2i* 
    x2*omzi + 3*li2spec6*NCi* 
    sqrtxz2i*x2*omzi +  
   3*li2spec7*NCi* 
    sqrtxz2i*x2*omzi -  
   3*li2spec8*NCi* 
    sqrtxz2i*x2*omzi +  
   3*lx*lspec5*NCi*sqrtxz2i*x2* 
    omzi - 3*lx*lspec6*NCi*sqrtxz2i* 
    x2*omzi - NC*li2spec5*sqrtxz2i* 
    x3*omzi + NC*li2spec6*sqrtxz2i* 
    x3*omzi + NC* 
    li2spec7*sqrtxz2i*x3* 
    omzi - NC*li2spec8* 
    sqrtxz2i*x3*omzi +  
   NC*lx*lspec5*sqrtxz2i*x3*omzi -  
   NC*lx*lspec6*sqrtxz2i*x3*omzi +  
   li2spec5*NCi*sqrtxz2i*x3* 
    omzi - li2spec6*NCi*sqrtxz2i* 
    x3*omzi - li2spec7* 
    NCi*sqrtxz2i*x3*omzi +  
   li2spec8*NCi* 
    sqrtxz2i*x3*omzi -  
   lx*lspec5*NCi*sqrtxz2i*x3* 
    omzi + lx*lspec6*NCi*sqrtxz2i* 
    x3*omzi + (NC*x*lx*omxmzi2)/4. -  
   (NC*x*lomz*omxmzi2)/4. -  
   (x*lx*NCi*omxmzi2)/4. +  
   (x*lomz*NCi*omxmzi2)/4. -  
   (3*NC*lx*x2*omxmzi2)/4. +  
   (3*NC*lomz*x2*omxmzi2)/4. +  
   (3*lx*NCi*x2*omxmzi2)/4. -  
   (3*lomz*NCi*x2*omxmzi2)/4. +  
   NC*lx*x3*omxmzi2 -  
   NC*lomz*x3*omxmzi2 -  
   lx*NCi*x3*omxmzi2 +  
   lomz*NCi*x3*omxmzi2 -  
   (NC*lx*x4*omxmzi2)/2. +  
   (NC*lomz*x4*omxmzi2)/2. +  
   (lx*NCi*x4*omxmzi2)/2. -  
   (lomz*NCi*x4*omxmzi2)/2. +  
   (NC*omxmzi)/4. - (3*NC*x*omxmzi)/4. -  
   (3*NC*lx*omxmzi)/4. + 2*NC*x*lx*omxmzi +  
   (3*NC*lomz*omxmzi)/4. -  
   2*NC*x*lomz*omxmzi - (NCi*omxmzi)/4. +  
   (3*x*NCi*omxmzi)/4. +  
   (3*lx*NCi*omxmzi)/4. -  
   (3*x*lx*NCi*omxmzi)/2. -  
   (3*lomz*NCi*omxmzi)/4. +  
   (3*x*lomz*NCi*omxmzi)/2. +  
   NC*x2*omxmzi -  
   (3*NC*lx*x2*omxmzi)/2. +  
   (3*NC*lomz*x2*omxmzi)/2. -  
   NCi*x2*omxmzi +  
   (lx*NCi*x2*omxmzi)/2. -  
   (lomz*NCi*x2*omxmzi)/2. -  
   (NC*x3*omxmzi)/2. +  
   NC*lx*x3*omxmzi -  
   NC*lomz*x3*omxmzi +  
   (NCi*x3*omxmzi)/2. +  
   (NC*lx*x4*xmzi2)/2. -  
   (NC*lz*x4*xmzi2)/2. -  
   (lx*NCi*x4*xmzi2)/2. +  
   (lz*NCi*x4*xmzi2)/2. +  
   (x*lx*NCi*xmzi)/2. -  
   (x*lz*NCi*xmzi)/2. -  
   lx*NCi*x2*xmzi +  
   lz*NCi*x2*xmzi -  
   (3*NC*lx*x3*xmzi)/2. +  
   (3*NC*lz*x3*xmzi)/2. +  
   (5*lx*NCi*x3*xmzi)/2. -  
   (5*lz*NCi*x3*xmzi)/2. - (199*NC*zi)/48. +  
   (28*NC*x*zi)/3. - (5*NC*Li2x*zi)/2. -  
   NC*x*Li2x*zi + (NC*Li2z*zi)/2. - NC*x*Li2z*zi +  
   (NC*li2omxzi*zi)/2. -  
   NC*x*li2omxzi*zi + (NC*lomx*zi)/4. +  
   (11*NC*x*lomx*zi)/2. + (61*NC*lx*zi)/8. -  
   (69*NC*x*lx*zi)/4. - 4*NC*lomx*lx*zi +  
   8*NC*x*lomx*lx*zi + (3*NC*lomz*zi)/4. +  
   (11*NC*x*lomz*zi)/2. + 2*NC*lomx*lomz*zi -  
   4*NC*x*lomx*lomz*zi -  
   2*NC*lx*lomz*zi + 10*NC*x*lx*lomz*zi -  
   (19*NC*lz*zi)/8. + (23*NC*x*lz*zi)/2. +  
   (3*NC*lomx*lz*zi)/2. -  
   3*NC*x*lomx*lz*zi + (5*NC*lx*lz*zi)/4. +  
   (3*NC*x*lx*lz*zi)/2. +  
   (5*NC*lomz*lz*zi)/2. -  
   5*NC*x*lomz*lz*zi + (9*NCi*zi)/16. -  
   3*x*NCi*zi - 2*Li2mx*NCi*zi -  
   4*x*Li2mx*NCi*zi + (Li2x*NCi*zi)/2. -  
   x*Li2x*NCi*zi -  
   (3*li2omxzi*NCi*zi)/2. +  
   3*x*li2omxzi*NCi*zi -  
   (lomx*NCi*zi)/4. -  
   (3*x*lomx*NCi*zi)/2. -  
   (53*lx*NCi*zi)/8. +  
   (69*x*lx*NCi*zi)/4. +  
   8*lomx*lx*NCi*zi -  
   16*x*lomx*lx*NCi*zi -  
   2*lx*lopx*NCi*zi -  
   4*x*lx*lopx*NCi*zi +  
   (lomz*NCi*zi)/4. -  
   (3*x*lomz*NCi*zi)/2. -  
   (9*lomx*lomz*NCi*zi)/2. +  
   9*x*lomx*lomz*NCi*zi +  
   7*lx*lomz*NCi*zi -  
   14*x*lx*lomz*NCi*zi +  
   (23*lz*NCi*zi)/8. -  
   (15*x*lz*NCi*zi)/2. -  
   (7*lomx*lz*NCi*zi)/2. +  
   7*x*lomx*lz*NCi*zi +  
   (17*lx*lz*NCi*zi)/4. -  
   (9*x*lx*lz*NCi*zi)/2. -  
   3*lomz*lz*NCi*zi +  
   6*x*lomz*lz*NCi*zi + (NC*pi2*zi)/6. +  
   (2*NC*x*pi2*zi)/3. + (3*NCi*pi2*zi)/4. -  
   (13*x*NCi*pi2*zi)/6. -  
   (NC*li2spec5*sqrtxz2i*zi)/2. -  
   (3*NC*x*li2spec5*sqrtxz2i*zi)/2. +  
   (NC*li2spec6*sqrtxz2i*zi)/2. +  
   (3*NC*x*li2spec6*sqrtxz2i*zi)/2. +  
   (NC*li2spec7*sqrtxz2i* 
      zi)/2. + (3*NC*x*li2spec7*sqrtxz2i*zi)/2. -  
   (NC*li2spec8*sqrtxz2i* 
      zi)/2. - (3*NC*x*li2spec8*sqrtxz2i*zi)/2. +  
   (NC*lx*lspec5*sqrtxz2i*zi)/2. +  
   (3*NC*x*lx*lspec5*sqrtxz2i*zi)/2. -  
   (NC*lx*lspec6*sqrtxz2i*zi)/2. -  
   (3*NC*x*lx*lspec6*sqrtxz2i*zi)/2. +  
   li2spec5*NCi*sqrtxz2i*zi +  
   3*x*li2spec5*NCi*sqrtxz2i*zi -  
   li2spec6*NCi*sqrtxz2i*zi -  
   3*x*li2spec6*NCi*sqrtxz2i*zi -  
   li2spec7*NCi* 
    sqrtxz2i*zi -  
   3*x*li2spec7*NCi* 
    sqrtxz2i*zi +  
   li2spec8*NCi* 
    sqrtxz2i*zi +  
   3*x*li2spec8*NCi* 
    sqrtxz2i*zi -  
   lx*lspec5*NCi*sqrtxz2i*zi -  
   3*x*lx*lspec5*NCi*sqrtxz2i*zi +  
   lx*lspec6*NCi*sqrtxz2i*zi +  
   3*x*lx*lspec6*NCi*sqrtxz2i*zi +  
   (13*NC*xi*zi)/9. +  
   (2*NC*lomx*xi*zi)/3. +  
   (2*NC*lomz*xi*zi)/3. +  
   (2*NC*lz*xi*zi)/3. - (161*NC*x2*zi)/18. -  
   3*NC*Li2x*x2*zi + NC*Li2z*x2*zi +  
   NC*li2omxzi*x2*zi -  
   (20*NC*lomx*x2*zi)/3. +  
   30*NC*lx*x2*zi -  
   8*NC*lomx*lx*x2*zi -  
   (20*NC*lomz*x2*zi)/3. +  
   4*NC*lomx*lomz*x2*zi -  
   6*NC*lx*lomz*x2*zi -  
   (79*NC*lz*x2*zi)/6. +  
   3*NC*lomx*lz*x2*zi +  
   NC*lx*lz*x2*zi +  
   5*NC*lomz*lz*x2*zi +  
   (11*NCi*x2*zi)/4. -  
   2*Li2mx*NCi*x2*zi +  
   Li2x*NCi*x2*zi -  
   Li2z*NCi*x2*zi -  
   li2omxzi*NCi*x2*zi -  
   (lomx*NCi*x2*zi)/2. -  
   (13*lx*NCi*x2*zi)/2. +  
   11*lomx*lx*NCi*x2*zi -  
   2*lx*lopx*NCi*x2*zi -  
   (lomz*NCi*x2*zi)/2. -  
   6*lomx*lomz*NCi*x2*zi +  
   9*lx*lomz*NCi*x2*zi +  
   3*lz*NCi*x2*zi -  
   5*lomx*lz*NCi*x2*zi +  
   6*lx*lz*NCi*x2*zi -  
   5*lomz*lz*NCi*x2*zi +  
   NCi*pi2*x2*zi -  
   2*NC*li2spec5*sqrtxz2i*x2*zi +  
   2*NC*li2spec6*sqrtxz2i*x2*zi +  
   2*NC*li2spec7*sqrtxz2i* 
    x2*zi - 2*NC*li2spec8*sqrtxz2i*x2*zi +  
   2*NC*lx*lspec5*sqrtxz2i*x2*zi -  
   2*NC*lx*lspec6*sqrtxz2i*x2*zi +  
   3*li2spec5*NCi*sqrtxz2i*x2* 
    zi - 3*li2spec6*NCi*sqrtxz2i* 
    x2*zi - 3*li2spec7* 
    NCi*sqrtxz2i*x2*zi +  
   3*li2spec8*NCi* 
    sqrtxz2i*x2*zi -  
   3*lx*lspec5*NCi*sqrtxz2i*x2* 
    zi + 3*lx*lspec6*NCi*sqrtxz2i* 
    x2*zi - NC*li2spec5*sqrtxz2i* 
    x3*zi + NC*li2spec6*sqrtxz2i* 
    x3*zi + NC*li2spec7* 
    sqrtxz2i*x3*zi -  
   NC*li2spec8*sqrtxz2i* 
    x3*zi + NC*lx*lspec5*sqrtxz2i* 
    x3*zi - NC*lx*lspec6*sqrtxz2i* 
    x3*zi + li2spec5*NCi* 
    sqrtxz2i*x3*zi -  
   li2spec6*NCi*sqrtxz2i*x3* 
    zi - li2spec7*NCi* 
    sqrtxz2i*x3*zi +  
   li2spec8*NCi* 
    sqrtxz2i*x3*zi -  
   lx*lspec5*NCi*sqrtxz2i*x3* 
    zi + lx*lspec6*NCi*sqrtxz2i* 
    x3*zi + 2*NC*omzi*zi -  
   4*NC*x*omzi*zi + 2*NC*Li2x*omzi*zi -  
   4*NC*x*Li2x*omzi*zi -  
   NC*li2omxzi*omzi*zi +  
   2*NC*x*li2omxzi*omzi*zi -  
   8*NC*lx*omzi*zi +  
   16*NC*x*lx*omzi*zi +  
   2*NC*lomx*lx*omzi*zi -  
   4*NC*x*lomx*lx*omzi*zi +  
   4*NC*lz*omzi*zi -  
   8*NC*x*lz*omzi*zi -  
   3*NC*lx*lz*omzi*zi +  
   6*NC*x*lx*lz*omzi*zi -  
   2*NCi*omzi*zi +  
   4*x*NCi*omzi*zi +  
   2*li2omxzi*NCi*omzi*zi -  
   4*x*li2omxzi*NCi*omzi*zi +  
   8*lx*NCi*omzi*zi -  
   16*x*lx*NCi*omzi*zi -  
   4*lz*NCi*omzi*zi +  
   8*x*lz*NCi*omzi*zi +  
   2*lx*lz*NCi*omzi*zi -  
   4*x*lx*lz*NCi*omzi*zi -  
   (NC*pi2*omzi*zi)/3. +  
   (2*NC*x*pi2*omzi*zi)/3. +  
   4*NC*x2*omzi*zi +  
   4*NC*Li2x*x2*omzi*zi -  
   2*NC*li2omxzi*x2*omzi*zi -  
   16*NC*lx*x2*omzi*zi +  
   4*NC*lomx*lx*x2*omzi*zi +  
   8*NC*lz*x2*omzi*zi -  
   6*NC*lx*lz*x2*omzi*zi -  
   2*NCi*x2*omzi*zi +  
   2*li2omxzi*NCi*x2*omzi*zi +  
   8*lx*NCi*x2*omzi*zi -  
   4*lz*NCi*x2*omzi*zi +  
   2*lx*lz*NCi*x2*omzi*zi -  
   (2*NC*pi2*x2*omzi*zi)/3. +  
   (51*NC*x*li2spec5*sqrtxz2i*z2)/4. -  
   (51*NC*x*li2spec6*sqrtxz2i*z2)/4. -  
   (51*NC*x*li2spec7*sqrtxz2i* 
      z2)/4. + (51*NC*x*li2spec8*sqrtxz2i*z2)/4. -  
   (51*NC*x*lx*lspec5*sqrtxz2i*z2)/4. +  
   (51*NC*x*lx*lspec6*sqrtxz2i*z2)/4. -  
   (11*x*li2spec5*NCi*sqrtxz2i*z2)/ 
    4. + (11*x*li2spec6*NCi*sqrtxz2i* 
      z2)/4. + (11*x*li2spec7* 
      NCi*sqrtxz2i*z2)/4. -  
   (11*x*li2spec8*NCi* 
      sqrtxz2i*z2)/4. +  
   (11*x*lx*lspec5*NCi*sqrtxz2i*z2)/ 
    4. - (11*x*lx*lspec6*NCi*sqrtxz2i* 
      z2)/4. - (NC*x*lx*xmzi2*z3)/2. +  
   (NC*x*lz*xmzi2*z3)/2. +  
   (x*lx*NCi*xmzi2*z3)/2. -  
   (x*lz*NCi*xmzi2*z3)/2. -  
   (7*NC*lomx2)/4. + (7*NC*x*lomx2)/2. -  
   (NC*z*lomx2)/4. + (NC*x*z*lomx2)/2. +  
   (9*NCi*lomx2)/4. -  
   (9*x*NCi*lomx2)/2. +  
   (z*NCi*lomx2)/4. -  
   (x*z*NCi*lomx2)/2. -  
   (7*NC*x2*lomx2)/2. -  
   (NC*z*x2*lomx2)/2. +  
   (9*NCi*x2*lomx2)/2. +  
   (z*NCi*x2*lomx2)/2. +  
   (NC*omzi*lomx2)/4. -  
   (NC*x*omzi*lomx2)/2. -  
   3*NCi*omzi*lomx2 +  
   6*x*NCi*omzi*lomx2 +  
   (NC*x2*omzi*lomx2)/2. -  
   4*NCi*x2*omzi*lomx2 +  
   (3*NC*zi*lomx2)/4. -  
   (3*NC*x*zi*lomx2)/2. -  
   3*NCi*zi*lomx2 +  
   6*x*NCi*zi*lomx2 +  
   (3*NC*x2*zi*lomx2)/2. -  
   4*NCi*x2*zi*lomx2 -  
   (5*NC*lx2)/4. + (33*NC*x*lx2)/2. -  
   NC*z*lx2 + (13*NCi*lx2)/4. -  
   (19*x*NCi*lx2)/2. + z*NCi*lx2 -  
   7*NC*x2*lx2 - 2*NC*z*x2*lx2 +  
   7*NCi*x2*lx2 +  
   2*z*NCi*x2*lx2 -  
   (NC*omzi*lx2)/4. +  
   (NC*x*omzi*lx2)/2. -  
   3*NCi*omzi*lx2 +  
   6*x*NCi*omzi*lx2 -  
   (NC*x2*omzi*lx2)/2. -  
   5*NCi*x2*omzi*lx2 -  
   2*NC*zi*lx2 - 2*NC*x*zi*lx2 -  
   (7*NCi*zi*lx2)/4. +  
   (11*x*NCi*zi*lx2)/2. -  
   2*NC*x2*zi*lx2 -  
   (7*NCi*x2*zi*lx2)/2. +  
   (5*NC*omzi*zi*lx2)/2. -  
   5*NC*x*omzi*zi*lx2 -  
   3*NCi*omzi*zi*lx2 +  
   6*x*NCi*omzi*zi*lx2 +  
   5*NC*x2*omzi*zi*lx2 -  
   3*NCi*x2*omzi*zi*lx2 -  
   (11*NC*lomz2)/4. + (11*NC*x*lomz2)/2. -  
   (NC*z*lomz2)/4. + (NC*x*z*lomz2)/2. +  
   2*NCi*lomz2 - 4*x*NCi*lomz2 +  
   (z*NCi*lomz2)/4. -  
   (x*z*NCi*lomz2)/2. -  
   (11*NC*x2*lomz2)/2. -  
   (NC*z*x2*lomz2)/2. +  
   4*NCi*x2*lomz2 +  
   (z*NCi*x2*lomz2)/2. +  
   NC*omzi*lomz2 -  
   2*NC*x*omzi*lomz2 -  
   (3*NCi*omzi*lomz2)/2. +  
   3*x*NCi*omzi*lomz2 +  
   2*NC*x2*omzi*lomz2 -  
   (5*NCi*x2*omzi*lomz2)/2. +  
   (5*NC*zi*lomz2)/4. -  
   (5*NC*x*zi*lomz2)/2. -  
   (5*NCi*zi*lomz2)/2. +  
   5*x*NCi*zi*lomz2 +  
   (5*NC*x2*zi*lomz2)/2. -  
   3*NCi*x2*zi*lomz2 -  
   (NC*lz2)/2. + NC*x*lz2 + (NC*z*lz2)/4. -  
   (NC*x*z*lz2)/2. + (NCi*lz2)/4. -  
   (x*NCi*lz2)/2. - (z*NCi*lz2)/4. +  
   (x*z*NCi*lz2)/2. - NC*x2*lz2 +  
   (NC*z*x2*lz2)/2. +  
   (NCi*x2*lz2)/2. -  
   (z*NCi*x2*lz2)/2. -  
   (NC*omzi*lz2)/4. +  
   (NC*x*omzi*lz2)/2. -  
   NCi*omzi*lz2 +  
   2*x*NCi*omzi*lz2 -  
   (NC*x2*omzi*lz2)/2. -  
   NCi*x2*omzi*lz2 +  
   (NC*zi*lz2)/4. - (NC*x*zi*lz2)/2. -  
   (5*NCi*zi*lz2)/4. +  
   (5*x*NCi*zi*lz2)/2. +  
   (NC*x2*zi*lz2)/2. -  
   2*NCi*x2*zi*lz2 +  
   (NC*omzi*zi*lz2)/2. -  
   NC*x*omzi*zi*lz2 +  
   NC*x2*omzi*zi*lz2 +  
   (-2*NC - 2*NC*x + NC*li2spec13 - 2*NC*x*li2spec13 -  
      NC*li2spec14 +  
      2*NC*x*li2spec14 -  
      NC*li2spec11 + 2*NC*x*li2spec11 +  
      (3*NC*lomx)/2. + 2*NC*x*lomx - (9*NC*lx)/2. -  
      6*NC*x*lx - 4*NC*lomx*lx + 8*NC*x*lomx*lx +  
      3*NC*lomz + 4*NC*x*lomz + 3*NC*lomx*lomz -  
      6*NC*x*lomx*lomz - 5*NC*lx*lomz +  
      10*NC*x*lx*lomz + NC*lx*lomxmz -  
      2*NC*x*lx*lomxmz - NC*lomz*lomxmz +  
      2*NC*x*lomz*lomxmz + NC*lx*lxmz -  
      2*NC*x*lx*lxmz + 3*NC*lz + 4*NC*x*lz +  
      2*NC*lomx*lz - 4*NC*x*lomx*lz -  
      5*NC*lx*lz + 10*NC*x*lx*lz + 4*NC*lomz*lz -  
      8*NC*x*lomz*lz - NC*lxmz*lz +  
      2*NC*x*lxmz*lz + (li2spec13*NCi)/2. -  
      x*li2spec13*NCi -  
      (li2spec11*NCi)/2. +  
      x*li2spec11*NCi -  
      (li2spec15*NCi)/2. +  
      x*li2spec15*NCi +  
      (li2spec12*NCi)/2. -  
      x*li2spec12*NCi -  
      2*lomx*NCi + 3*lx*NCi +  
      (11*lomx*lx*NCi)/2. -  
      11*x*lomx*lx*NCi - (3*lomz*NCi)/2. -  
      (7*lomx*lomz*NCi)/2. +  
      7*x*lomx*lomz*NCi + 5*lx*lomz*NCi -  
      10*x*lx*lomz*NCi -  
      (lx*lomxmz*NCi)/2. +  
      x*lx*lomxmz*NCi +  
      (lomz*lomxmz*NCi)/2. -  
      x*lomz*lomxmz*NCi -  
      (lx*lxmz*NCi)/2. + x*lx*lxmz*NCi -  
      (3*lz*NCi)/2. - 3*lomx*lz*NCi +  
      6*x*lomx*lz*NCi + 4*lx*lz*NCi -  
      8*x*lx*lz*NCi - 2*lomz*lz*NCi +  
      4*x*lomz*lz*NCi + (lxmz*lz*NCi)/2. -  
      x*lxmz*lz*NCi - (NC*pi2)/6. +  
      (NC*x*pi2)/3. + (7*NCi*pi2)/12. -  
      (7*x*NCi*pi2)/6. + 2*NC*x2 +  
      2*NC*li2spec13*x2 -  
      2*NC*li2spec14*x2 -  
      2*NC*li2spec11*x2 - 2*NC*lomx*x2 +  
      6*NC*lx*x2 - 8*NC*lomx*lx*x2 -  
      4*NC*lomz*x2 + 6*NC*lomx*lomz*x2 -  
      10*NC*lx*lomz*x2 +  
      2*NC*lx*lomxmz*x2 -  
      2*NC*lomz*lomxmz*x2 +  
      2*NC*lx*lxmz*x2 - 4*NC*lz*x2 +  
      4*NC*lomx*lz*x2 - 10*NC*lx*lz*x2 +  
      8*NC*lomz*lz*x2 - 2*NC*lxmz*lz*x2 +  
      li2spec13*NCi*x2 -  
      li2spec11*NCi*x2 -  
      li2spec15*NCi*x2 +  
      li2spec12*NCi*x2 +  
      11*lomx*lx*NCi*x2 -  
      7*lomx*lomz*NCi*x2 +  
      10*lx*lomz*NCi*x2 -  
      lx*lomxmz*NCi*x2 +  
      lomz*lomxmz*NCi*x2 -  
      lx*lxmz*NCi*x2 -  
      6*lomx*lz*NCi*x2 +  
      8*lx*lz*NCi*x2 -  
      4*lomz*lz*NCi*x2 +  
      lxmz*lz*NCi*x2 - (NC*pi2*x2)/3. +  
      (7*NCi*pi2*x2)/6. -  
      (NC*li2spec13*omzi)/2. +  
      NC*x*li2spec13*omzi +  
      (NC*li2spec14*omzi)/2. -  
      NC*x*li2spec14*omzi +  
      (NC*li2spec11*omzi)/2. -  
      NC*x*li2spec11*omzi -  
      (NC*lomx*omzi)/2. + (3*NC*lx*omzi)/2. +  
      2*NC*lomx*lx*omzi -  
      4*NC*x*lomx*lx*omzi - NC*lomz*omzi -  
      (3*NC*lomx*lomz*omzi)/2. +  
      3*NC*x*lomx*lomz*omzi +  
      (5*NC*lx*lomz*omzi)/2. -  
      5*NC*x*lx*lomz*omzi -  
      (NC*lx*lomxmz*omzi)/2. +  
      NC*x*lx*lomxmz*omzi +  
      (NC*lomz*lomxmz*omzi)/2. -  
      NC*x*lomz*lomxmz*omzi -  
      (NC*lx*lxmz*omzi)/2. +  
      NC*x*lx*lxmz*omzi - NC*lz*omzi -  
      NC*lomx*lz*omzi +  
      2*NC*x*lomx*lz*omzi +  
      (5*NC*lx*lz*omzi)/2. -  
      5*NC*x*lx*lz*omzi -  
      2*NC*lomz*lz*omzi +  
      4*NC*x*lomz*lz*omzi +  
      (NC*lxmz*lz*omzi)/2. -  
      NC*x*lxmz*lz*omzi +  
      (3*NCi*omzi)/2. -  
      (Li2z*NCi*omzi)/2. +  
      x*Li2z*NCi*omzi -  
      (li2spec13*NCi*omzi)/2. +  
      x*li2spec13*NCi*omzi +  
      li2spec11*NCi*omzi -  
      2*x*li2spec11*NCi*omzi +  
      (li2spec15*NCi*omzi)/2. -  
      x*li2spec15*NCi*omzi -  
      li2spec12*NCi*omzi +  
      2*x*li2spec12*NCi*omzi +  
      lomx*NCi*omzi -  
      (3*lx*NCi*omzi)/2. -  
      8*lomx*lx*NCi*omzi +  
      16*x*lomx*lx*NCi*omzi +  
      lomz*NCi*omzi +  
      (9*lomx*lomz*NCi*omzi)/2. -  
      9*x*lomx*lomz*NCi*omzi -  
      7*lx*lomz*NCi*omzi +  
      14*x*lx*lomz*NCi*omzi +  
      (lx*lomxmz*NCi*omzi)/2. -  
      x*lx*lomxmz*NCi*omzi -  
      (lomz*lomxmz*NCi*omzi)/2. +  
      x*lomz*lomxmz*NCi*omzi +  
      lx*lxmz*NCi*omzi -  
      2*x*lx*lxmz*NCi*omzi +  
      (lz*NCi*omzi)/2. +  
      5*lomx*lz*NCi*omzi -  
      10*x*lomx*lz*NCi*omzi -  
      (13*lx*lz*NCi*omzi)/2. +  
      13*x*lx*lz*NCi*omzi +  
      3*lomz*lz*NCi*omzi -  
      6*x*lomz*lz*NCi*omzi -  
      lxmz*lz*NCi*omzi +  
      2*x*lxmz*lz*NCi*omzi +  
      (NC*pi2*omzi)/12. - (NC*x*pi2*omzi)/6. -  
      (5*NCi*pi2*omzi)/6. +  
      (5*x*NCi*pi2*omzi)/3. -  
      NC*li2spec13*x2*omzi +  
      NC*li2spec14*x2*omzi +  
      NC*li2spec11*x2*omzi +  
      4*NC*lomx*lx*x2*omzi -  
      3*NC*lomx*lomz*x2*omzi +  
      5*NC*lx*lomz*x2*omzi -  
      NC*lx*lomxmz*x2*omzi +  
      NC*lomz*lomxmz*x2*omzi -  
      NC*lx*lxmz*x2*omzi -  
      2*NC*lomx*lz*x2*omzi +  
      5*NC*lx*lz*x2*omzi -  
      4*NC*lomz*lz*x2*omzi +  
      NC*lxmz*lz*x2*omzi -  
      NCi*x2*omzi -  
      li2spec13*NCi*x2*omzi +  
      li2spec11*NCi*x2*omzi +  
      li2spec15*NCi*x2* 
       omzi - li2spec12*NCi* 
       x2*omzi +  
      2*lomx*NCi*x2*omzi -  
      3*lx*NCi*x2*omzi -  
      11*lomx*lx*NCi*x2*omzi +  
      lomz*NCi*x2*omzi +  
      7*lomx*lomz*NCi*x2*omzi -  
      10*lx*lomz*NCi*x2*omzi +  
      lx*lomxmz*NCi*x2*omzi -  
      lomz*lomxmz*NCi*x2*omzi +  
      lx*lxmz*NCi*x2*omzi +  
      2*lz*NCi*x2*omzi +  
      6*lomx*lz*NCi*x2*omzi -  
      8*lx*lz*NCi*x2*omzi +  
      4*lomz*lz*NCi*x2*omzi -  
      lxmz*lz*NCi*x2*omzi +  
      (NC*pi2*x2*omzi)/6. -  
      (7*NCi*pi2*x2*omzi)/6. -  
      (NC*li2spec13*zi)/2. +  
      NC*x*li2spec13*zi +  
      (NC*li2spec14*zi)/2. -  
      NC*x*li2spec14*zi +  
      (NC*li2spec11*zi)/2. -  
      NC*x*li2spec11*zi -  
      (NC*lomx*zi)/2. + (3*NC*lx*zi)/2. +  
      2*NC*lomx*lx*zi - 4*NC*x*lomx*lx*zi -  
      NC*lomz*zi - (3*NC*lomx*lomz*zi)/2. +  
      3*NC*x*lomx*lomz*zi +  
      (5*NC*lx*lomz*zi)/2. -  
      5*NC*x*lx*lomz*zi -  
      (NC*lx*lomxmz*zi)/2. +  
      NC*x*lx*lomxmz*zi +  
      (NC*lomz*lomxmz*zi)/2. -  
      NC*x*lomz*lomxmz*zi -  
      (NC*lx*lxmz*zi)/2. +  
      NC*x*lx*lxmz*zi - NC*lz*zi -  
      NC*lomx*lz*zi + 2*NC*x*lomx*lz*zi +  
      (5*NC*lx*lz*zi)/2. - 5*NC*x*lx*lz*zi -  
      2*NC*lomz*lz*zi + 4*NC*x*lomz*lz*zi +  
      (NC*lxmz*lz*zi)/2. -  
      NC*x*lxmz*lz*zi + (3*NCi*zi)/2. +  
      (Li2z*NCi*zi)/2. - x*Li2z*NCi*zi -  
      li2spec13*NCi*zi +  
      2*x*li2spec13*NCi*zi +  
      (li2spec11*NCi*zi)/2. -  
      x*li2spec11*NCi*zi +  
      li2spec15*NCi*zi -  
      2*x*li2spec15*NCi*zi -  
      (li2spec12*NCi*zi)/2. +  
      x*li2spec12*NCi*zi +  
      lomx*NCi*zi - (3*lx*NCi*zi)/2. -  
      (17*lomx*lx*NCi*zi)/2. +  
      17*x*lomx*lx*NCi*zi +  
      (lomz*NCi*zi)/2. +  
      6*lomx*lomz*NCi*zi -  
      12*x*lomx*lomz*NCi*zi -  
      8*lx*lomz*NCi*zi +  
      16*x*lx*lomz*NCi*zi +  
      lx*lomxmz*NCi*zi -  
      2*x*lx*lomxmz*NCi*zi -  
      lomz*lomxmz*NCi*zi +  
      2*x*lomz*lomxmz*NCi*zi +  
      (lx*lxmz*NCi*zi)/2. -  
      x*lx*lxmz*NCi*zi +  
      lz*NCi*zi +  
      4*lomx*lz*NCi*zi -  
      8*x*lomx*lz*NCi*zi -  
      (11*lx*lz*NCi*zi)/2. +  
      11*x*lx*lz*NCi*zi +  
      3*lomz*lz*NCi*zi -  
      6*x*lomz*lz*NCi*zi -  
      (lxmz*lz*NCi*zi)/2. +  
      x*lxmz*lz*NCi*zi +  
      (NC*pi2*zi)/12. - (NC*x*pi2*zi)/6. -  
      (11*NCi*pi2*zi)/12. +  
      (11*x*NCi*pi2*zi)/6. -  
      NC*li2spec13*x2*zi +  
      NC*li2spec14*x2*zi +  
      NC*li2spec11*x2*zi +  
      4*NC*lomx*lx*x2*zi -  
      3*NC*lomx*lomz*x2*zi +  
      5*NC*lx*lomz*x2*zi -  
      NC*lx*lomxmz*x2*zi +  
      NC*lomz*lomxmz*x2*zi -  
      NC*lx*lxmz*x2*zi -  
      2*NC*lomx*lz*x2*zi +  
      5*NC*lx*lz*x2*zi -  
      4*NC*lomz*lz*x2*zi +  
      NC*lxmz*lz*x2*zi -  
      NCi*x2*zi -  
      li2spec13*NCi*x2*zi +  
      li2spec11*NCi*x2*zi +  
      li2spec15*NCi*x2*zi -  
      li2spec12*NCi*x2*zi +  
      2*lomx*NCi*x2*zi -  
      3*lx*NCi*x2*zi -  
      11*lomx*lx*NCi*x2*zi +  
      2*lomz*NCi*x2*zi +  
      7*lomx*lomz*NCi*x2*zi -  
      10*lx*lomz*NCi*x2*zi +  
      lx*lomxmz*NCi*x2*zi -  
      lomz*lomxmz*NCi*x2*zi +  
      lx*lxmz*NCi*x2*zi +  
      lz*NCi*x2*zi +  
      6*lomx*lz*NCi*x2*zi -  
      8*lx*lz*NCi*x2*zi +  
      4*lomz*lz*NCi*x2*zi -  
      lxmz*lz*NCi*x2*zi +  
      (NC*pi2*x2*zi)/6. -  
      (7*NCi*pi2*x2*zi)/6. +  
      (NC*lomx2)/2. - NC*x*lomx2 -  
      2*NCi*lomx2 + 4*x*NCi*lomx2 +  
      NC*x2*lomx2 -  
      4*NCi*x2*lomx2 -  
      (NC*omzi*lomx2)/4. +  
      (NC*x*omzi*lomx2)/2. +  
      3*NCi*omzi*lomx2 -  
      6*x*NCi*omzi*lomx2 -  
      (NC*x2*omzi*lomx2)/2. +  
      4*NCi*x2*omzi*lomx2 -  
      (NC*zi*lomx2)/4. +  
      (NC*x*zi*lomx2)/2. +  
      3*NCi*zi*lomx2 -  
      6*x*NCi*zi*lomx2 -  
      (NC*x2*zi*lomx2)/2. +  
      4*NCi*x2*zi*lomx2 +  
      3*NC*lx2 - 6*NC*x*lx2 -  
      (13*NCi*lx2)/4. + (13*x*NCi*lx2)/2. +  
      6*NC*x2*lx2 -  
      (13*NCi*x2*lx2)/2. -  
      (3*NC*omzi*lx2)/2. +  
      3*NC*x*omzi*lx2 +  
      (19*NCi*omzi*lx2)/4. -  
      (19*x*NCi*omzi*lx2)/2. -  
      3*NC*x2*omzi*lx2 +  
      (13*NCi*x2*omzi*lx2)/2. -  
      (3*NC*zi*lx2)/2. + 3*NC*x*zi*lx2 +  
      5*NCi*zi*lx2 -  
      10*x*NCi*zi*lx2 -  
      3*NC*x2*zi*lx2 +  
      (13*NCi*x2*zi*lx2)/2. +  
      (3*NC*lomz2)/2. - 3*NC*x*lomz2 -  
      (3*NCi*lomz2)/2. +  
      3*x*NCi*lomz2 + 3*NC*x2*lomz2 -  
      3*NCi*x2*lomz2 -  
      (3*NC*omzi*lomz2)/4. +  
      (3*NC*x*omzi*lomz2)/2. +  
      (7*NCi*omzi*lomz2)/4. -  
      (7*x*NCi*omzi*lomz2)/2. -  
      (3*NC*x2*omzi*lomz2)/2. +  
      3*NCi*x2*omzi*lomz2 -  
      (3*NC*zi*lomz2)/4. +  
      (3*NC*x*zi*lomz2)/2. +  
      (11*NCi*zi*lomz2)/4. -  
      (11*x*NCi*zi*lomz2)/2. -  
      (3*NC*x2*zi*lomz2)/2. +  
      3*NCi*x2*zi*lomz2 +  
      2*NC*lz2 - 4*NC*x*lz2 -  
      (5*NCi*lz2)/4. + (5*x*NCi*lz2)/2. +  
      4*NC*x2*lz2 -  
      (5*NCi*x2*lz2)/2. -  
      NC*omzi*lz2 + 2*NC*x*omzi*lz2 +  
      (9*NCi*omzi*lz2)/4. -  
      (9*x*NCi*omzi*lz2)/2. -  
      2*NC*x2*omzi*lz2 +  
      (5*NCi*x2*omzi*lz2)/2. -  
      NC*zi*lz2 + 2*NC*x*zi*lz2 +  
      (3*NCi*zi*lz2)/2. -  
      3*x*NCi*zi*lz2 -  
      2*NC*x2*zi*lz2 +  
      (5*NCi*x2*zi*lz2)/2.)*Tu1 +  
   (-2*NC - 2*NC*x - NC*li2spec11 +  
      2*NC*x*li2spec11 - NC*li2spec16 +  
      2*NC*x*li2spec16 +  
      NC*li2spec18 -  
      2*NC*x*li2spec18 + (3*NC*lomx)/2. +  
      2*NC*x*lomx - (9*NC*lx)/2. - 6*NC*x*lx -  
      3*NC*lomx*lx + 6*NC*x*lomx*lx + 3*NC*lomz +  
      4*NC*x*lomz + 2*NC*lomx*lomz -  
      4*NC*x*lomx*lomz - 6*NC*lx*lomz +  
      12*NC*x*lx*lomz + NC*lx*lxmz -  
      2*NC*x*lx*lxmz + 3*NC*lz + 4*NC*x*lz +  
      2*NC*lomx*lz - 4*NC*x*lomx*lz -  
      6*NC*lx*lz + 12*NC*x*lx*lz + 5*NC*lomz*lz -  
      10*NC*x*lomz*lz - NC*lxmz*lz +  
      2*NC*x*lxmz*lz + NC*lx*lmopxpz -  
      2*NC*x*lx*lmopxpz - NC*lomz*lmopxpz +  
      2*NC*x*lomz*lmopxpz -  
      (li2spec11*NCi)/2. +  
      x*li2spec11*NCi +  
      (li2spec12*NCi)/2. -  
      x*li2spec12*NCi -  
      (li2spec16*NCi)/2. +  
      x*li2spec16*NCi +  
      (li2spec17*NCi)/2. -  
      x*li2spec17*NCi -  
      2*lomx*NCi + 3*lx*NCi +  
      5*lomx*lx*NCi - 10*x*lomx*lx*NCi -  
      (3*lomz*NCi)/2. - 3*lomx*lomz*NCi +  
      6*x*lomx*lomz*NCi +  
      (9*lx*lomz*NCi)/2. -  
      9*x*lx*lomz*NCi - (lx*lxmz*NCi)/2. +  
      x*lx*lxmz*NCi - (3*lz*NCi)/2. -  
      3*lomx*lz*NCi + 6*x*lomx*lz*NCi +  
      (9*lx*lz*NCi)/2. - 9*x*lx*lz*NCi -  
      (5*lomz*lz*NCi)/2. +  
      5*x*lomz*lz*NCi + (lxmz*lz*NCi)/2. -  
      x*lxmz*lz*NCi -  
      (lx*lmopxpz*NCi)/2. +  
      x*lx*lmopxpz*NCi +  
      (lomz*lmopxpz*NCi)/2. -  
      x*lomz*lmopxpz*NCi - (NC*pi2)/6. +  
      (NC*x*pi2)/3. + (7*NCi*pi2)/12. -  
      (7*x*NCi*pi2)/6. + 2*NC*x2 -  
      2*NC*li2spec11*x2 -  
      2*NC*li2spec16*x2 +  
      2*NC*li2spec18*x2 -  
      2*NC*lomx*x2 + 6*NC*lx*x2 -  
      6*NC*lomx*lx*x2 - 4*NC*lomz*x2 +  
      4*NC*lomx*lomz*x2 -  
      12*NC*lx*lomz*x2 + 2*NC*lx*lxmz*x2 -  
      4*NC*lz*x2 + 4*NC*lomx*lz*x2 -  
      12*NC*lx*lz*x2 + 10*NC*lomz*lz*x2 -  
      2*NC*lxmz*lz*x2 +  
      2*NC*lx*lmopxpz*x2 -  
      2*NC*lomz*lmopxpz*x2 -  
      li2spec11*NCi*x2 +  
      li2spec12*NCi*x2 -  
      li2spec16*NCi*x2 +  
      li2spec17*NCi*x2 +  
      10*lomx*lx*NCi*x2 -  
      6*lomx*lomz*NCi*x2 +  
      9*lx*lomz*NCi*x2 -  
      lx*lxmz*NCi*x2 -  
      6*lomx*lz*NCi*x2 +  
      9*lx*lz*NCi*x2 -  
      5*lomz*lz*NCi*x2 +  
      lxmz*lz*NCi*x2 -  
      lx*lmopxpz*NCi*x2 +  
      lomz*lmopxpz*NCi*x2 -  
      (NC*pi2*x2)/3. + (7*NCi*pi2*x2)/6. +  
      (NC*li2spec11*omzi)/2. -  
      NC*x*li2spec11*omzi +  
      (NC*li2spec16*omzi)/2. -  
      NC*x*li2spec16*omzi -  
      (NC*li2spec18*omzi)/2. +  
      NC*x*li2spec18*omzi -  
      (NC*lomx*omzi)/2. + (3*NC*lx*omzi)/2. +  
      (3*NC*lomx*lx*omzi)/2. -  
      3*NC*x*lomx*lx*omzi - NC*lomz*omzi -  
      NC*lomx*lomz*omzi +  
      2*NC*x*lomx*lomz*omzi +  
      3*NC*lx*lomz*omzi -  
      6*NC*x*lx*lomz*omzi -  
      (NC*lx*lxmz*omzi)/2. +  
      NC*x*lx*lxmz*omzi - NC*lz*omzi -  
      NC*lomx*lz*omzi +  
      2*NC*x*lomx*lz*omzi +  
      3*NC*lx*lz*omzi - 6*NC*x*lx*lz*omzi -  
      (5*NC*lomz*lz*omzi)/2. +  
      5*NC*x*lomz*lz*omzi +  
      (NC*lxmz*lz*omzi)/2. -  
      NC*x*lxmz*lz*omzi -  
      (NC*lx*lmopxpz*omzi)/2. +  
      NC*x*lx*lmopxpz*omzi +  
      (NC*lomz*lmopxpz*omzi)/2. -  
      NC*x*lomz*lmopxpz*omzi +  
      (3*NCi*omzi)/2. -  
      (Li2z*NCi*omzi)/2. +  
      x*Li2z*NCi*omzi +  
      li2spec11*NCi*omzi -  
      2*x*li2spec11*NCi*omzi -  
      li2spec12*NCi*omzi +  
      2*x*li2spec12*NCi*omzi +  
      (li2spec16*NCi*omzi)/2. -  
      x*li2spec16*NCi*omzi -  
      (li2spec17*NCi*omzi)/ 
       2. + x*li2spec17*NCi* 
       omzi + lomx*NCi*omzi -  
      (3*lx*NCi*omzi)/2. -  
      (15*lomx*lx*NCi*omzi)/2. +  
      15*x*lomx*lx*NCi*omzi +  
      lomz*NCi*omzi +  
      4*lomx*lomz*NCi*omzi -  
      8*x*lomx*lomz*NCi*omzi -  
      (13*lx*lomz*NCi*omzi)/2. +  
      13*x*lx*lomz*NCi*omzi +  
      lx*lxmz*NCi*omzi -  
      2*x*lx*lxmz*NCi*omzi +  
      (lz*NCi*omzi)/2. +  
      5*lomx*lz*NCi*omzi -  
      10*x*lomx*lz*NCi*omzi -  
      7*lx*lz*NCi*omzi +  
      14*x*lx*lz*NCi*omzi +  
      (7*lomz*lz*NCi*omzi)/2. -  
      7*x*lomz*lz*NCi*omzi -  
      lxmz*lz*NCi*omzi +  
      2*x*lxmz*lz*NCi*omzi +  
      (lx*lmopxpz*NCi*omzi)/2. -  
      x*lx*lmopxpz*NCi*omzi -  
      (lomz*lmopxpz*NCi*omzi)/2. +  
      x*lomz*lmopxpz*NCi*omzi +  
      (NC*pi2*omzi)/12. - (NC*x*pi2*omzi)/6. -  
      (5*NCi*pi2*omzi)/6. +  
      (5*x*NCi*pi2*omzi)/3. +  
      NC*li2spec11*x2*omzi +  
      NC*li2spec16*x2*omzi -  
      NC*li2spec18*x2*omzi +  
      3*NC*lomx*lx*x2*omzi -  
      2*NC*lomx*lomz*x2*omzi +  
      6*NC*lx*lomz*x2*omzi -  
      NC*lx*lxmz*x2*omzi -  
      2*NC*lomx*lz*x2*omzi +  
      6*NC*lx*lz*x2*omzi -  
      5*NC*lomz*lz*x2*omzi +  
      NC*lxmz*lz*x2*omzi -  
      NC*lx*lmopxpz*x2*omzi +  
      NC*lomz*lmopxpz*x2*omzi -  
      NCi*x2*omzi +  
      li2spec11*NCi*x2*omzi -  
      li2spec12*NCi*x2* 
       omzi + li2spec16*NCi*x2* 
       omzi - li2spec17*NCi* 
       x2*omzi +  
      2*lomx*NCi*x2*omzi -  
      3*lx*NCi*x2*omzi -  
      10*lomx*lx*NCi*x2*omzi +  
      lomz*NCi*x2*omzi +  
      6*lomx*lomz*NCi*x2*omzi -  
      9*lx*lomz*NCi*x2*omzi +  
      lx*lxmz*NCi*x2*omzi +  
      2*lz*NCi*x2*omzi +  
      6*lomx*lz*NCi*x2*omzi -  
      9*lx*lz*NCi*x2*omzi +  
      5*lomz*lz*NCi*x2*omzi -  
      lxmz*lz*NCi*x2*omzi +  
      lx*lmopxpz*NCi*x2*omzi -  
      lomz*lmopxpz*NCi*x2*omzi +  
      (NC*pi2*x2*omzi)/6. -  
      (7*NCi*pi2*x2*omzi)/6. +  
      (NC*li2spec11*zi)/2. -  
      NC*x*li2spec11*zi +  
      (NC*li2spec16*zi)/2. -  
      NC*x*li2spec16*zi -  
      (NC*li2spec18*zi)/2. +  
      NC*x*li2spec18*zi -  
      (NC*lomx*zi)/2. + (3*NC*lx*zi)/2. +  
      (3*NC*lomx*lx*zi)/2. -  
      3*NC*x*lomx*lx*zi - NC*lomz*zi -  
      NC*lomx*lomz*zi +  
      2*NC*x*lomx*lomz*zi +  
      3*NC*lx*lomz*zi - 6*NC*x*lx*lomz*zi -  
      (NC*lx*lxmz*zi)/2. +  
      NC*x*lx*lxmz*zi - NC*lz*zi -  
      NC*lomx*lz*zi + 2*NC*x*lomx*lz*zi +  
      3*NC*lx*lz*zi - 6*NC*x*lx*lz*zi -  
      (5*NC*lomz*lz*zi)/2. +  
      5*NC*x*lomz*lz*zi +  
      (NC*lxmz*lz*zi)/2. -  
      NC*x*lxmz*lz*zi -  
      (NC*lx*lmopxpz*zi)/2. +  
      NC*x*lx*lmopxpz*zi +  
      (NC*lomz*lmopxpz*zi)/2. -  
      NC*x*lomz*lmopxpz*zi +  
      (3*NCi*zi)/2. + (Li2z*NCi*zi)/2. -  
      x*Li2z*NCi*zi +  
      (li2spec11*NCi*zi)/2. -  
      x*li2spec11*NCi*zi -  
      (li2spec12*NCi*zi)/2. +  
      x*li2spec12*NCi*zi +  
      li2spec16*NCi*zi -  
      2*x*li2spec16*NCi*zi -  
      li2spec17*NCi*zi +  
      2*x*li2spec17*NCi*zi +  
      lomx*NCi*zi - (3*lx*NCi*zi)/2. -  
      (15*lomx*lx*NCi*zi)/2. +  
      15*x*lomx*lx*NCi*zi +  
      (lomz*NCi*zi)/2. +  
      5*lomx*lomz*NCi*zi -  
      10*x*lomx*lomz*NCi*zi -  
      7*lx*lomz*NCi*zi +  
      14*x*lx*lomz*NCi*zi +  
      (lx*lxmz*NCi*zi)/2. -  
      x*lx*lxmz*NCi*zi +  
      lz*NCi*zi +  
      4*lomx*lz*NCi*zi -  
      8*x*lomx*lz*NCi*zi -  
      (13*lx*lz*NCi*zi)/2. +  
      13*x*lx*lz*NCi*zi +  
      4*lomz*lz*NCi*zi -  
      8*x*lomz*lz*NCi*zi -  
      (lxmz*lz*NCi*zi)/2. +  
      x*lxmz*lz*NCi*zi +  
      lx*lmopxpz*NCi*zi -  
      2*x*lx*lmopxpz*NCi*zi -  
      lomz*lmopxpz*NCi*zi +  
      2*x*lomz*lmopxpz*NCi*zi +  
      (NC*pi2*zi)/12. - (NC*x*pi2*zi)/6. -  
      (11*NCi*pi2*zi)/12. +  
      (11*x*NCi*pi2*zi)/6. +  
      NC*li2spec11*x2*zi +  
      NC*li2spec16*x2*zi -  
      NC*li2spec18*x2*zi +  
      3*NC*lomx*lx*x2*zi -  
      2*NC*lomx*lomz*x2*zi +  
      6*NC*lx*lomz*x2*zi -  
      NC*lx*lxmz*x2*zi -  
      2*NC*lomx*lz*x2*zi +  
      6*NC*lx*lz*x2*zi -  
      5*NC*lomz*lz*x2*zi +  
      NC*lxmz*lz*x2*zi -  
      NC*lx*lmopxpz*x2*zi +  
      NC*lomz*lmopxpz*x2*zi -  
      NCi*x2*zi +  
      li2spec11*NCi*x2*zi -  
      li2spec12*NCi*x2*zi +  
      li2spec16*NCi*x2*zi -  
      li2spec17*NCi*x2* 
       zi + 2*lomx*NCi*x2*zi -  
      3*lx*NCi*x2*zi -  
      10*lomx*lx*NCi*x2*zi +  
      2*lomz*NCi*x2*zi +  
      6*lomx*lomz*NCi*x2*zi -  
      9*lx*lomz*NCi*x2*zi +  
      lx*lxmz*NCi*x2*zi +  
      lz*NCi*x2*zi +  
      6*lomx*lz*NCi*x2*zi -  
      9*lx*lz*NCi*x2*zi +  
      5*lomz*lz*NCi*x2*zi -  
      lxmz*lz*NCi*x2*zi +  
      lx*lmopxpz*NCi*x2*zi -  
      lomz*lmopxpz*NCi*x2*zi +  
      (NC*pi2*x2*zi)/6. -  
      (7*NCi*pi2*x2*zi)/6. +  
      (NC*lomx2)/2. - NC*x*lomx2 -  
      2*NCi*lomx2 + 4*x*NCi*lomx2 +  
      NC*x2*lomx2 -  
      4*NCi*x2*lomx2 -  
      (NC*omzi*lomx2)/4. +  
      (NC*x*omzi*lomx2)/2. +  
      3*NCi*omzi*lomx2 -  
      6*x*NCi*omzi*lomx2 -  
      (NC*x2*omzi*lomx2)/2. +  
      4*NCi*x2*omzi*lomx2 -  
      (NC*zi*lomx2)/4. +  
      (NC*x*zi*lomx2)/2. +  
      3*NCi*zi*lomx2 -  
      6*x*NCi*zi*lomx2 -  
      (NC*x2*zi*lomx2)/2. +  
      4*NCi*x2*zi*lomx2 +  
      (7*NC*lx2)/2. - 7*NC*x*lx2 -  
      3*NCi*lx2 + 6*x*NCi*lx2 +  
      7*NC*x2*lx2 - 6*NCi*x2*lx2 -  
      (7*NC*omzi*lx2)/4. +  
      (7*NC*x*omzi*lx2)/2. +  
      (9*NCi*omzi*lx2)/2. -  
      9*x*NCi*omzi*lx2 -  
      (7*NC*x2*omzi*lx2)/2. +  
      6*NCi*x2*omzi*lx2 -  
      (7*NC*zi*lx2)/4. +  
      (7*NC*x*zi*lx2)/2. +  
      (9*NCi*zi*lx2)/2. -  
      9*x*NCi*zi*lx2 -  
      (7*NC*x2*zi*lx2)/2. +  
      6*NCi*x2*zi*lx2 +  
      2*NC*lomz2 - 4*NC*x*lomz2 -  
      (5*NCi*lomz2)/4. +  
      (5*x*NCi*lomz2)/2. +  
      4*NC*x2*lomz2 -  
      (5*NCi*x2*lomz2)/2. -  
      NC*omzi*lomz2 +  
      2*NC*x*omzi*lomz2 +  
      (3*NCi*omzi*lomz2)/2. -  
      3*x*NCi*omzi*lomz2 -  
      2*NC*x2*omzi*lomz2 +  
      (5*NCi*x2*omzi*lomz2)/2. -  
      NC*zi*lomz2 + 2*NC*x*zi*lomz2 +  
      (9*NCi*zi*lomz2)/4. -  
      (9*x*NCi*zi*lomz2)/2. -  
      2*NC*x2*zi*lomz2 +  
      (5*NCi*x2*zi*lomz2)/2. +  
      2*NC*lz2 - 4*NC*x*lz2 -  
      (5*NCi*lz2)/4. + (5*x*NCi*lz2)/2. +  
      4*NC*x2*lz2 -  
      (5*NCi*x2*lz2)/2. -  
      NC*omzi*lz2 + 2*NC*x*omzi*lz2 +  
      (9*NCi*omzi*lz2)/4. -  
      (9*x*NCi*omzi*lz2)/2. -  
      2*NC*x2*omzi*lz2 +  
      (5*NCi*x2*omzi*lz2)/2. -  
      NC*zi*lz2 + 2*NC*x*zi*lz2 +  
      (3*NCi*zi*lz2)/2. -  
      3*x*NCi*zi*lz2 -  
      2*NC*x2*zi*lz2 +  
      (5*NCi*x2*zi*lz2)/2.)*Tu2 +  
   (-2*NC - 2*NC*x + NC*li2spec9 -  
      2*NC*x*li2spec9 + NC*li2spec13 -  
      2*NC*x*li2spec13 +  
      NC*li2spec18 -  
      2*NC*x*li2spec18 + (3*NC*lomx)/2. +  
      2*NC*x*lomx - (9*NC*lx)/2. - 6*NC*x*lx -  
      3*NC*lomx*lx + 6*NC*x*lomx*lx + 3*NC*lomz +  
      4*NC*x*lomz + NC*lomx*lomz -  
      2*NC*x*lomx*lomz - 6*NC*lx*lomz +  
      12*NC*x*lx*lomz + NC*lx*lomxmz -  
      2*NC*x*lx*lomxmz - NC*lomz*lomxmz +  
      2*NC*x*lomz*lomxmz + 3*NC*lz + 4*NC*x*lz +  
      NC*lomx*lz - 2*NC*x*lomx*lz - 6*NC*lx*lz +  
      12*NC*x*lx*lz + 5*NC*lomz*lz -  
      10*NC*x*lomz*lz + NC*lx*lmxpz -  
      2*NC*x*lx*lmxpz - NC*lz*lmxpz +  
      2*NC*x*lz*lmxpz +  
      (li2spec9*NCi)/2. -  
      x*li2spec9*NCi +  
      (li2spec13*NCi)/2. -  
      x*li2spec13*NCi -  
      (li2spec15*NCi)/2. +  
      x*li2spec15*NCi -  
      (li2spec10*NCi)/2. +  
      x*li2spec10*NCi -  
      2*lomx*NCi + 3*lx*NCi +  
      6*lomx*lx*NCi - 12*x*lomx*lx*NCi -  
      (3*lomz*NCi)/2. -  
      (7*lomx*lomz*NCi)/2. +  
      7*x*lomx*lomz*NCi +  
      (9*lx*lomz*NCi)/2. -  
      9*x*lx*lomz*NCi -  
      (lx*lomxmz*NCi)/2. +  
      x*lx*lomxmz*NCi +  
      (lomz*lomxmz*NCi)/2. -  
      x*lomz*lomxmz*NCi - (3*lz*NCi)/2. -  
      (7*lomx*lz*NCi)/2. +  
      7*x*lomx*lz*NCi + (9*lx*lz*NCi)/2. -  
      9*x*lx*lz*NCi - (3*lomz*lz*NCi)/2. +  
      3*x*lomz*lz*NCi - (lx*lmxpz*NCi)/2. +  
      x*lx*lmxpz*NCi + (lz*lmxpz*NCi)/2. -  
      x*lz*lmxpz*NCi - (5*NC*pi2)/6. +  
      (5*NC*x*pi2)/3. + (7*NCi*pi2)/12. -  
      (7*x*NCi*pi2)/6. + 2*NC*x2 +  
      2*NC*li2spec9*x2 +  
      2*NC*li2spec13*x2 +  
      2*NC*li2spec18*x2 -  
      2*NC*lomx*x2 + 6*NC*lx*x2 -  
      6*NC*lomx*lx*x2 - 4*NC*lomz*x2 +  
      2*NC*lomx*lomz*x2 -  
      12*NC*lx*lomz*x2 +  
      2*NC*lx*lomxmz*x2 -  
      2*NC*lomz*lomxmz*x2 - 4*NC*lz*x2 +  
      2*NC*lomx*lz*x2 - 12*NC*lx*lz*x2 +  
      10*NC*lomz*lz*x2 + 2*NC*lx*lmxpz*x2 -  
      2*NC*lz*lmxpz*x2 +  
      li2spec9*NCi*x2 +  
      li2spec13*NCi*x2 -  
      li2spec15*NCi*x2 -  
      li2spec10*NCi*x2 +  
      12*lomx*lx*NCi*x2 -  
      7*lomx*lomz*NCi*x2 +  
      9*lx*lomz*NCi*x2 -  
      lx*lomxmz*NCi*x2 +  
      lomz*lomxmz*NCi*x2 -  
      7*lomx*lz*NCi*x2 +  
      9*lx*lz*NCi*x2 -  
      3*lomz*lz*NCi*x2 -  
      lx*lmxpz*NCi*x2 +  
      lz*lmxpz*NCi*x2 - (5*NC*pi2*x2)/3. +  
      (7*NCi*pi2*x2)/6. -  
      (NC*li2spec9*omzi)/2. +  
      NC*x*li2spec9*omzi -  
      (NC*li2spec13*omzi)/2. +  
      NC*x*li2spec13*omzi -  
      (NC*li2spec18*omzi)/2. +  
      NC*x*li2spec18*omzi -  
      (NC*lomx*omzi)/2. + (3*NC*lx*omzi)/2. +  
      (3*NC*lomx*lx*omzi)/2. -  
      3*NC*x*lomx*lx*omzi - NC*lomz*omzi -  
      (NC*lomx*lomz*omzi)/2. +  
      NC*x*lomx*lomz*omzi +  
      3*NC*lx*lomz*omzi -  
      6*NC*x*lx*lomz*omzi -  
      (NC*lx*lomxmz*omzi)/2. +  
      NC*x*lx*lomxmz*omzi +  
      (NC*lomz*lomxmz*omzi)/2. -  
      NC*x*lomz*lomxmz*omzi -  
      NC*lz*omzi - (NC*lomx*lz*omzi)/2. +  
      NC*x*lomx*lz*omzi +  
      3*NC*lx*lz*omzi - 6*NC*x*lx*lz*omzi -  
      (5*NC*lomz*lz*omzi)/2. +  
      5*NC*x*lomz*lz*omzi -  
      (NC*lx*lmxpz*omzi)/2. +  
      NC*x*lx*lmxpz*omzi +  
      (NC*lz*lmxpz*omzi)/2. -  
      NC*x*lz*lmxpz*omzi +  
      (3*NCi*omzi)/2. -  
      (Li2z*NCi*omzi)/2. +  
      x*Li2z*NCi*omzi -  
      li2spec9*NCi*omzi +  
      2*x*li2spec9*NCi*omzi -  
      (li2spec13*NCi*omzi)/2. +  
      x*li2spec13*NCi*omzi +  
      (li2spec15*NCi*omzi)/2. -  
      x*li2spec15*NCi*omzi +  
      li2spec10*NCi*omzi -  
      2*x*li2spec10*NCi*omzi +  
      lomx*NCi*omzi -  
      (3*lx*NCi*omzi)/2. -  
      9*lomx*lx*NCi*omzi +  
      18*x*lomx*lx*NCi*omzi +  
      lomz*NCi*omzi +  
      (9*lomx*lomz*NCi*omzi)/2. -  
      9*x*lomx*lomz*NCi*omzi -  
      6*lx*lomz*NCi*omzi +  
      12*x*lx*lomz*NCi*omzi +  
      (lx*lomxmz*NCi*omzi)/2. -  
      x*lx*lomxmz*NCi*omzi -  
      (lomz*lomxmz*NCi*omzi)/2. +  
      x*lomz*lomxmz*NCi*omzi +  
      (lz*NCi*omzi)/2. +  
      6*lomx*lz*NCi*omzi -  
      12*x*lomx*lz*NCi*omzi -  
      (15*lx*lz*NCi*omzi)/2. +  
      15*x*lx*lz*NCi*omzi +  
      2*lomz*lz*NCi*omzi -  
      4*x*lomz*lz*NCi*omzi +  
      lx*lmxpz*NCi*omzi -  
      2*x*lx*lmxpz*NCi*omzi -  
      lz*lmxpz*NCi*omzi +  
      2*x*lz*lmxpz*NCi*omzi +  
      (5*NC*pi2*omzi)/12. -  
      (5*NC*x*pi2*omzi)/6. -  
      (5*NCi*pi2*omzi)/6. +  
      (5*x*NCi*pi2*omzi)/3. -  
      NC*li2spec9*x2*omzi -  
      NC*li2spec13*x2*omzi -  
      NC*li2spec18*x2*omzi +  
      3*NC*lomx*lx*x2*omzi -  
      NC*lomx*lomz*x2*omzi +  
      6*NC*lx*lomz*x2*omzi -  
      NC*lx*lomxmz*x2*omzi +  
      NC*lomz*lomxmz*x2*omzi -  
      NC*lomx*lz*x2*omzi +  
      6*NC*lx*lz*x2*omzi -  
      5*NC*lomz*lz*x2*omzi -  
      NC*lx*lmxpz*x2*omzi +  
      NC*lz*lmxpz*x2*omzi -  
      NCi*x2*omzi -  
      li2spec9*NCi*x2*omzi -  
      li2spec13*NCi*x2*omzi +  
      li2spec15*NCi*x2* 
       omzi + li2spec10*NCi* 
       x2*omzi +  
      2*lomx*NCi*x2*omzi -  
      3*lx*NCi*x2*omzi -  
      12*lomx*lx*NCi*x2*omzi +  
      lomz*NCi*x2*omzi +  
      7*lomx*lomz*NCi*x2*omzi -  
      9*lx*lomz*NCi*x2*omzi +  
      lx*lomxmz*NCi*x2*omzi -  
      lomz*lomxmz*NCi*x2*omzi +  
      2*lz*NCi*x2*omzi +  
      7*lomx*lz*NCi*x2*omzi -  
      9*lx*lz*NCi*x2*omzi +  
      3*lomz*lz*NCi*x2*omzi +  
      lx*lmxpz*NCi*x2*omzi -  
      lz*lmxpz*NCi*x2*omzi +  
      (5*NC*pi2*x2*omzi)/6. -  
      (7*NCi*pi2*x2*omzi)/6. -  
      (NC*li2spec9*zi)/2. +  
      NC*x*li2spec9*zi -  
      (NC*li2spec13*zi)/2. +  
      NC*x*li2spec13*zi -  
      (NC*li2spec18*zi)/2. +  
      NC*x*li2spec18*zi -  
      (NC*lomx*zi)/2. + (3*NC*lx*zi)/2. +  
      (3*NC*lomx*lx*zi)/2. -  
      3*NC*x*lomx*lx*zi - NC*lomz*zi -  
      (NC*lomx*lomz*zi)/2. +  
      NC*x*lomx*lomz*zi +  
      3*NC*lx*lomz*zi - 6*NC*x*lx*lomz*zi -  
      (NC*lx*lomxmz*zi)/2. +  
      NC*x*lx*lomxmz*zi +  
      (NC*lomz*lomxmz*zi)/2. -  
      NC*x*lomz*lomxmz*zi - NC*lz*zi -  
      (NC*lomx*lz*zi)/2. +  
      NC*x*lomx*lz*zi + 3*NC*lx*lz*zi -  
      6*NC*x*lx*lz*zi -  
      (5*NC*lomz*lz*zi)/2. +  
      5*NC*x*lomz*lz*zi -  
      (NC*lx*lmxpz*zi)/2. +  
      NC*x*lx*lmxpz*zi +  
      (NC*lz*lmxpz*zi)/2. -  
      NC*x*lz*lmxpz*zi + (3*NCi*zi)/2. +  
      (Li2z*NCi*zi)/2. - x*Li2z*NCi*zi -  
      (li2spec9*NCi*zi)/2. +  
      x*li2spec9*NCi*zi -  
      li2spec13*NCi*zi +  
      2*x*li2spec13*NCi*zi +  
      li2spec15*NCi*zi -  
      2*x*li2spec15*NCi*zi +  
      (li2spec10*NCi*zi)/2. -  
      x*li2spec10*NCi*zi +  
      lomx*NCi*zi - (3*lx*NCi*zi)/2. -  
      9*lomx*lx*NCi*zi +  
      18*x*lomx*lx*NCi*zi +  
      (lomz*NCi*zi)/2. +  
      6*lomx*lomz*NCi*zi -  
      12*x*lomx*lomz*NCi*zi -  
      (15*lx*lomz*NCi*zi)/2. +  
      15*x*lx*lomz*NCi*zi +  
      lx*lomxmz*NCi*zi -  
      2*x*lx*lomxmz*NCi*zi -  
      lomz*lomxmz*NCi*zi +  
      2*x*lomz*lomxmz*NCi*zi +  
      lz*NCi*zi +  
      (9*lomx*lz*NCi*zi)/2. -  
      9*x*lomx*lz*NCi*zi -  
      6*lx*lz*NCi*zi +  
      12*x*lx*lz*NCi*zi +  
      (5*lomz*lz*NCi*zi)/2. -  
      5*x*lomz*lz*NCi*zi +  
      (lx*lmxpz*NCi*zi)/2. -  
      x*lx*lmxpz*NCi*zi -  
      (lz*lmxpz*NCi*zi)/2. +  
      x*lz*lmxpz*NCi*zi +  
      (5*NC*pi2*zi)/12. - (5*NC*x*pi2*zi)/6. -  
      (11*NCi*pi2*zi)/12. +  
      (11*x*NCi*pi2*zi)/6. -  
      NC*li2spec9*x2*zi -  
      NC*li2spec13*x2*zi -  
      NC*li2spec18*x2*zi +  
      3*NC*lomx*lx*x2*zi -  
      NC*lomx*lomz*x2*zi +  
      6*NC*lx*lomz*x2*zi -  
      NC*lx*lomxmz*x2*zi +  
      NC*lomz*lomxmz*x2*zi -  
      NC*lomx*lz*x2*zi +  
      6*NC*lx*lz*x2*zi -  
      5*NC*lomz*lz*x2*zi -  
      NC*lx*lmxpz*x2*zi +  
      NC*lz*lmxpz*x2*zi -  
      NCi*x2*zi -  
      li2spec9*NCi*x2*zi -  
      li2spec13*NCi*x2*zi +  
      li2spec15*NCi*x2*zi +  
      li2spec10*NCi*x2*zi +  
      2*lomx*NCi*x2*zi -  
      3*lx*NCi*x2*zi -  
      12*lomx*lx*NCi*x2*zi +  
      2*lomz*NCi*x2*zi +  
      7*lomx*lomz*NCi*x2*zi -  
      9*lx*lomz*NCi*x2*zi +  
      lx*lomxmz*NCi*x2*zi -  
      lomz*lomxmz*NCi*x2*zi +  
      lz*NCi*x2*zi +  
      7*lomx*lz*NCi*x2*zi -  
      9*lx*lz*NCi*x2*zi +  
      3*lomz*lz*NCi*x2*zi +  
      lx*lmxpz*NCi*x2*zi -  
      lz*lmxpz*NCi*x2*zi +  
      (5*NC*pi2*x2*zi)/6. -  
      (7*NCi*pi2*x2*zi)/6. +  
      (3*NC*lomx2)/2. - 3*NC*x*lomx2 -  
      2*NCi*lomx2 + 4*x*NCi*lomx2 +  
      3*NC*x2*lomx2 -  
      4*NCi*x2*lomx2 -  
      (3*NC*omzi*lomx2)/4. +  
      (3*NC*x*omzi*lomx2)/2. +  
      3*NCi*omzi*lomx2 -  
      6*x*NCi*omzi*lomx2 -  
      (3*NC*x2*omzi*lomx2)/2. +  
      4*NCi*x2*omzi*lomx2 -  
      (3*NC*zi*lomx2)/4. +  
      (3*NC*x*zi*lomx2)/2. +  
      3*NCi*zi*lomx2 -  
      6*x*NCi*zi*lomx2 -  
      (3*NC*x2*zi*lomx2)/2. +  
      4*NCi*x2*zi*lomx2 +  
      (7*NC*lx2)/2. - 7*NC*x*lx2 -  
      (7*NCi*lx2)/2. + 7*x*NCi*lx2 +  
      7*NC*x2*lx2 - 7*NCi*x2*lx2 -  
      (7*NC*omzi*lx2)/4. +  
      (7*NC*x*omzi*lx2)/2. +  
      (21*NCi*omzi*lx2)/4. -  
      (21*x*NCi*omzi*lx2)/2. -  
      (7*NC*x2*omzi*lx2)/2. +  
      7*NCi*x2*omzi*lx2 -  
      (7*NC*zi*lx2)/4. +  
      (7*NC*x*zi*lx2)/2. +  
      (21*NCi*zi*lx2)/4. -  
      (21*x*NCi*zi*lx2)/2. -  
      (7*NC*x2*zi*lx2)/2. +  
      7*NCi*x2*zi*lx2 +  
      (5*NC*lomz2)/2. - 5*NC*x*lomz2 -  
      (3*NCi*lomz2)/2. +  
      3*x*NCi*lomz2 + 5*NC*x2*lomz2 -  
      3*NCi*x2*lomz2 -  
      (5*NC*omzi*lomz2)/4. +  
      (5*NC*x*omzi*lomz2)/2. +  
      (7*NCi*omzi*lomz2)/4. -  
      (7*x*NCi*omzi*lomz2)/2. -  
      (5*NC*x2*omzi*lomz2)/2. +  
      3*NCi*x2*omzi*lomz2 -  
      (5*NC*zi*lomz2)/4. +  
      (5*NC*x*zi*lomz2)/2. +  
      (11*NCi*zi*lomz2)/4. -  
      (11*x*NCi*zi*lomz2)/2. -  
      (5*NC*x2*zi*lomz2)/2. +  
      3*NCi*x2*zi*lomz2 +  
      (5*NC*lz2)/2. - 5*NC*x*lz2 -  
      (3*NCi*lz2)/2. + 3*x*NCi*lz2 +  
      5*NC*x2*lz2 - 3*NCi*x2*lz2 -  
      (5*NC*omzi*lz2)/4. +  
      (5*NC*x*omzi*lz2)/2. +  
      (11*NCi*omzi*lz2)/4. -  
      (11*x*NCi*omzi*lz2)/2. -  
      (5*NC*x2*omzi*lz2)/2. +  
      3*NCi*x2*omzi*lz2 -  
      (5*NC*zi*lz2)/4. +  
      (5*NC*x*zi*lz2)/2. +  
      (7*NCi*zi*lz2)/4. -  
      (7*x*NCi*zi*lz2)/2. -  
      (5*NC*x2*zi*lz2)/2. +  
      3*NCi*x2*zi*lz2)*Tu3 +  
   (-2*NC - 2*NC*x + NC*li2spec9 -  
      2*NC*x*li2spec9 -  
      NC*li2spec14 +  
      2*NC*x*li2spec14 -  
      NC*li2spec16 + 2*NC*x*li2spec16 +  
      (3*NC*lomx)/2. + 2*NC*x*lomx - (9*NC*lx)/2. -  
      6*NC*x*lx - 4*NC*lomx*lx + 8*NC*x*lomx*lx +  
      3*NC*lomz + 4*NC*x*lomz + 2*NC*lomx*lomz -  
      4*NC*x*lomx*lomz - 5*NC*lx*lomz +  
      10*NC*x*lx*lomz + 3*NC*lz + 4*NC*x*lz +  
      3*NC*lomx*lz - 6*NC*x*lomx*lz -  
      5*NC*lx*lz + 10*NC*x*lx*lz + 4*NC*lomz*lz -  
      8*NC*x*lomz*lz + NC*lx*lmxpz -  
      2*NC*x*lx*lmxpz - NC*lz*lmxpz +  
      2*NC*x*lz*lmxpz + NC*lx*lmopxpz -  
      2*NC*x*lx*lmopxpz - NC*lomz*lmopxpz +  
      2*NC*x*lomz*lmopxpz +  
      (li2spec9*NCi)/2. -  
      x*li2spec9*NCi -  
      (li2spec16*NCi)/2. +  
      x*li2spec16*NCi -  
      (li2spec10*NCi)/2. +  
      x*li2spec10*NCi +  
      (li2spec17*NCi)/2. -  
      x*li2spec17*NCi -  
      2*lomx*NCi + 3*lx*NCi +  
      (11*lomx*lx*NCi)/2. -  
      11*x*lomx*lx*NCi - (3*lomz*NCi)/2. -  
      3*lomx*lomz*NCi +  
      6*x*lomx*lomz*NCi + 4*lx*lomz*NCi -  
      8*x*lx*lomz*NCi - (3*lz*NCi)/2. -  
      (7*lomx*lz*NCi)/2. +  
      7*x*lomx*lz*NCi + 5*lx*lz*NCi -  
      10*x*lx*lz*NCi - 2*lomz*lz*NCi +  
      4*x*lomz*lz*NCi - (lx*lmxpz*NCi)/2. +  
      x*lx*lmxpz*NCi + (lz*lmxpz*NCi)/2. -  
      x*lz*lmxpz*NCi -  
      (lx*lmopxpz*NCi)/2. +  
      x*lx*lmopxpz*NCi +  
      (lomz*lmopxpz*NCi)/2. -  
      x*lomz*lmopxpz*NCi - (NC*pi2)/6. +  
      (NC*x*pi2)/3. + (7*NCi*pi2)/12. -  
      (7*x*NCi*pi2)/6. + 2*NC*x2 +  
      2*NC*li2spec9*x2 -  
      2*NC*li2spec14*x2 -  
      2*NC*li2spec16*x2 - 2*NC*lomx*x2 +  
      6*NC*lx*x2 - 8*NC*lomx*lx*x2 -  
      4*NC*lomz*x2 + 4*NC*lomx*lomz*x2 -  
      10*NC*lx*lomz*x2 - 4*NC*lz*x2 +  
      6*NC*lomx*lz*x2 - 10*NC*lx*lz*x2 +  
      8*NC*lomz*lz*x2 + 2*NC*lx*lmxpz*x2 -  
      2*NC*lz*lmxpz*x2 +  
      2*NC*lx*lmopxpz*x2 -  
      2*NC*lomz*lmopxpz*x2 +  
      li2spec9*NCi*x2 -  
      li2spec16*NCi*x2 -  
      li2spec10*NCi*x2 +  
      li2spec17*NCi*x2 +  
      11*lomx*lx*NCi*x2 -  
      6*lomx*lomz*NCi*x2 +  
      8*lx*lomz*NCi*x2 -  
      7*lomx*lz*NCi*x2 +  
      10*lx*lz*NCi*x2 -  
      4*lomz*lz*NCi*x2 -  
      lx*lmxpz*NCi*x2 +  
      lz*lmxpz*NCi*x2 -  
      lx*lmopxpz*NCi*x2 +  
      lomz*lmopxpz*NCi*x2 -  
      (NC*pi2*x2)/3. + (7*NCi*pi2*x2)/6. -  
      (NC*li2spec9*omzi)/2. +  
      NC*x*li2spec9*omzi +  
      (NC*li2spec14*omzi)/2. -  
      NC*x*li2spec14*omzi +  
      (NC*li2spec16*omzi)/2. -  
      NC*x*li2spec16*omzi -  
      (NC*lomx*omzi)/2. + (3*NC*lx*omzi)/2. +  
      2*NC*lomx*lx*omzi -  
      4*NC*x*lomx*lx*omzi - NC*lomz*omzi -  
      NC*lomx*lomz*omzi +  
      2*NC*x*lomx*lomz*omzi +  
      (5*NC*lx*lomz*omzi)/2. -  
      5*NC*x*lx*lomz*omzi - NC*lz*omzi -  
      (3*NC*lomx*lz*omzi)/2. +  
      3*NC*x*lomx*lz*omzi +  
      (5*NC*lx*lz*omzi)/2. -  
      5*NC*x*lx*lz*omzi -  
      2*NC*lomz*lz*omzi +  
      4*NC*x*lomz*lz*omzi -  
      (NC*lx*lmxpz*omzi)/2. +  
      NC*x*lx*lmxpz*omzi +  
      (NC*lz*lmxpz*omzi)/2. -  
      NC*x*lz*lmxpz*omzi -  
      (NC*lx*lmopxpz*omzi)/2. +  
      NC*x*lx*lmopxpz*omzi +  
      (NC*lomz*lmopxpz*omzi)/2. -  
      NC*x*lomz*lmopxpz*omzi +  
      (3*NCi*omzi)/2. -  
      (Li2z*NCi*omzi)/2. +  
      x*Li2z*NCi*omzi -  
      li2spec9*NCi*omzi +  
      2*x*li2spec9*NCi*omzi +  
      (li2spec16*NCi*omzi)/2. -  
      x*li2spec16*NCi*omzi +  
      li2spec10*NCi*omzi -  
      2*x*li2spec10*NCi*omzi -  
      (li2spec17*NCi*omzi)/ 
       2. + x*li2spec17*NCi* 
       omzi + lomx*NCi*omzi -  
      (3*lx*NCi*omzi)/2. -  
      (17*lomx*lx*NCi*omzi)/2. +  
      17*x*lomx*lx*NCi*omzi +  
      lomz*NCi*omzi +  
      4*lomx*lomz*NCi*omzi -  
      8*x*lomx*lomz*NCi*omzi -  
      (11*lx*lomz*NCi*omzi)/2. +  
      11*x*lx*lomz*NCi*omzi +  
      (lz*NCi*omzi)/2. +  
      6*lomx*lz*NCi*omzi -  
      12*x*lomx*lz*NCi*omzi -  
      8*lx*lz*NCi*omzi +  
      16*x*lx*lz*NCi*omzi +  
      (5*lomz*lz*NCi*omzi)/2. -  
      5*x*lomz*lz*NCi*omzi +  
      lx*lmxpz*NCi*omzi -  
      2*x*lx*lmxpz*NCi*omzi -  
      lz*lmxpz*NCi*omzi +  
      2*x*lz*lmxpz*NCi*omzi +  
      (lx*lmopxpz*NCi*omzi)/2. -  
      x*lx*lmopxpz*NCi*omzi -  
      (lomz*lmopxpz*NCi*omzi)/2. +  
      x*lomz*lmopxpz*NCi*omzi +  
      (NC*pi2*omzi)/12. - (NC*x*pi2*omzi)/6. -  
      (5*NCi*pi2*omzi)/6. +  
      (5*x*NCi*pi2*omzi)/3. -  
      NC*li2spec9*x2*omzi +  
      NC*li2spec14*x2*omzi +  
      NC*li2spec16*x2*omzi +  
      4*NC*lomx*lx*x2*omzi -  
      2*NC*lomx*lomz*x2*omzi +  
      5*NC*lx*lomz*x2*omzi -  
      3*NC*lomx*lz*x2*omzi +  
      5*NC*lx*lz*x2*omzi -  
      4*NC*lomz*lz*x2*omzi -  
      NC*lx*lmxpz*x2*omzi +  
      NC*lz*lmxpz*x2*omzi -  
      NC*lx*lmopxpz*x2*omzi +  
      NC*lomz*lmopxpz*x2*omzi -  
      NCi*x2*omzi -  
      li2spec9*NCi*x2*omzi +  
      li2spec16*NCi*x2*omzi +  
      li2spec10*NCi*x2* 
       omzi - li2spec17*NCi* 
       x2*omzi +  
      2*lomx*NCi*x2*omzi -  
      3*lx*NCi*x2*omzi -  
      11*lomx*lx*NCi*x2*omzi +  
      lomz*NCi*x2*omzi +  
      6*lomx*lomz*NCi*x2*omzi -  
      8*lx*lomz*NCi*x2*omzi +  
      2*lz*NCi*x2*omzi +  
      7*lomx*lz*NCi*x2*omzi -  
      10*lx*lz*NCi*x2*omzi +  
      4*lomz*lz*NCi*x2*omzi +  
      lx*lmxpz*NCi*x2*omzi -  
      lz*lmxpz*NCi*x2*omzi +  
      lx*lmopxpz*NCi*x2*omzi -  
      lomz*lmopxpz*NCi*x2*omzi +  
      (NC*pi2*x2*omzi)/6. -  
      (7*NCi*pi2*x2*omzi)/6. -  
      (NC*li2spec9*zi)/2. +  
      NC*x*li2spec9*zi +  
      (NC*li2spec14*zi)/2. -  
      NC*x*li2spec14*zi +  
      (NC*li2spec16*zi)/2. -  
      NC*x*li2spec16*zi - (NC*lomx*zi)/2. +  
      (3*NC*lx*zi)/2. + 2*NC*lomx*lx*zi -  
      4*NC*x*lomx*lx*zi - NC*lomz*zi -  
      NC*lomx*lomz*zi +  
      2*NC*x*lomx*lomz*zi +  
      (5*NC*lx*lomz*zi)/2. -  
      5*NC*x*lx*lomz*zi - NC*lz*zi -  
      (3*NC*lomx*lz*zi)/2. +  
      3*NC*x*lomx*lz*zi +  
      (5*NC*lx*lz*zi)/2. - 5*NC*x*lx*lz*zi -  
      2*NC*lomz*lz*zi + 4*NC*x*lomz*lz*zi -  
      (NC*lx*lmxpz*zi)/2. +  
      NC*x*lx*lmxpz*zi +  
      (NC*lz*lmxpz*zi)/2. -  
      NC*x*lz*lmxpz*zi -  
      (NC*lx*lmopxpz*zi)/2. +  
      NC*x*lx*lmopxpz*zi +  
      (NC*lomz*lmopxpz*zi)/2. -  
      NC*x*lomz*lmopxpz*zi +  
      (3*NCi*zi)/2. + (Li2z*NCi*zi)/2. -  
      x*Li2z*NCi*zi -  
      (li2spec9*NCi*zi)/2. +  
      x*li2spec9*NCi*zi +  
      li2spec16*NCi*zi -  
      2*x*li2spec16*NCi*zi +  
      (li2spec10*NCi*zi)/2. -  
      x*li2spec10*NCi*zi -  
      li2spec17*NCi*zi +  
      2*x*li2spec17*NCi*zi +  
      lomx*NCi*zi - (3*lx*NCi*zi)/2. -  
      8*lomx*lx*NCi*zi +  
      16*x*lomx*lx*NCi*zi +  
      (lomz*NCi*zi)/2. +  
      5*lomx*lomz*NCi*zi -  
      10*x*lomx*lomz*NCi*zi -  
      (13*lx*lomz*NCi*zi)/2. +  
      13*x*lx*lomz*NCi*zi +  
      lz*NCi*zi +  
      (9*lomx*lz*NCi*zi)/2. -  
      9*x*lomx*lz*NCi*zi -  
      7*lx*lz*NCi*zi +  
      14*x*lx*lz*NCi*zi +  
      (7*lomz*lz*NCi*zi)/2. -  
      7*x*lomz*lz*NCi*zi +  
      (lx*lmxpz*NCi*zi)/2. -  
      x*lx*lmxpz*NCi*zi -  
      (lz*lmxpz*NCi*zi)/2. +  
      x*lz*lmxpz*NCi*zi +  
      lx*lmopxpz*NCi*zi -  
      2*x*lx*lmopxpz*NCi*zi -  
      lomz*lmopxpz*NCi*zi +  
      2*x*lomz*lmopxpz*NCi*zi +  
      (NC*pi2*zi)/12. - (NC*x*pi2*zi)/6. -  
      (11*NCi*pi2*zi)/12. +  
      (11*x*NCi*pi2*zi)/6. -  
      NC*li2spec9*x2*zi +  
      NC*li2spec14*x2*zi +  
      NC*li2spec16*x2*zi +  
      4*NC*lomx*lx*x2*zi -  
      2*NC*lomx*lomz*x2*zi +  
      5*NC*lx*lomz*x2*zi -  
      3*NC*lomx*lz*x2*zi +  
      5*NC*lx*lz*x2*zi -  
      4*NC*lomz*lz*x2*zi -  
      NC*lx*lmxpz*x2*zi +  
      NC*lz*lmxpz*x2*zi -  
      NC*lx*lmopxpz*x2*zi +  
      NC*lomz*lmopxpz*x2*zi -  
      NCi*x2*zi -  
      li2spec9*NCi*x2*zi +  
      li2spec16*NCi*x2*zi +  
      li2spec10*NCi*x2*zi -  
      li2spec17*NCi*x2* 
       zi + 2*lomx*NCi*x2*zi -  
      3*lx*NCi*x2*zi -  
      11*lomx*lx*NCi*x2*zi +  
      2*lomz*NCi*x2*zi +  
      6*lomx*lomz*NCi*x2*zi -  
      8*lx*lomz*NCi*x2*zi +  
      lz*NCi*x2*zi +  
      7*lomx*lz*NCi*x2*zi -  
      10*lx*lz*NCi*x2*zi +  
      4*lomz*lz*NCi*x2*zi +  
      lx*lmxpz*NCi*x2*zi -  
      lz*lmxpz*NCi*x2*zi +  
      lx*lmopxpz*NCi*x2*zi -  
      lomz*lmopxpz*NCi*x2*zi +  
      (NC*pi2*x2*zi)/6. -  
      (7*NCi*pi2*x2*zi)/6. +  
      (NC*lomx2)/2. - NC*x*lomx2 -  
      2*NCi*lomx2 + 4*x*NCi*lomx2 +  
      NC*x2*lomx2 -  
      4*NCi*x2*lomx2 -  
      (NC*omzi*lomx2)/4. +  
      (NC*x*omzi*lomx2)/2. +  
      3*NCi*omzi*lomx2 -  
      6*x*NCi*omzi*lomx2 -  
      (NC*x2*omzi*lomx2)/2. +  
      4*NCi*x2*omzi*lomx2 -  
      (NC*zi*lomx2)/4. +  
      (NC*x*zi*lomx2)/2. +  
      3*NCi*zi*lomx2 -  
      6*x*NCi*zi*lomx2 -  
      (NC*x2*zi*lomx2)/2. +  
      4*NCi*x2*zi*lomx2 +  
      3*NC*lx2 - 6*NC*x*lx2 -  
      (13*NCi*lx2)/4. + (13*x*NCi*lx2)/2. +  
      6*NC*x2*lx2 -  
      (13*NCi*x2*lx2)/2. -  
      (3*NC*omzi*lx2)/2. +  
      3*NC*x*omzi*lx2 +  
      5*NCi*omzi*lx2 -  
      10*x*NCi*omzi*lx2 -  
      3*NC*x2*omzi*lx2 +  
      (13*NCi*x2*omzi*lx2)/2. -  
      (3*NC*zi*lx2)/2. + 3*NC*x*zi*lx2 +  
      (19*NCi*zi*lx2)/4. -  
      (19*x*NCi*zi*lx2)/2. -  
      3*NC*x2*zi*lx2 +  
      (13*NCi*x2*zi*lx2)/2. +  
      2*NC*lomz2 - 4*NC*x*lomz2 -  
      (5*NCi*lomz2)/4. +  
      (5*x*NCi*lomz2)/2. +  
      4*NC*x2*lomz2 -  
      (5*NCi*x2*lomz2)/2. -  
      NC*omzi*lomz2 +  
      2*NC*x*omzi*lomz2 +  
      (3*NCi*omzi*lomz2)/2. -  
      3*x*NCi*omzi*lomz2 -  
      2*NC*x2*omzi*lomz2 +  
      (5*NCi*x2*omzi*lomz2)/2. -  
      NC*zi*lomz2 + 2*NC*x*zi*lomz2 +  
      (9*NCi*zi*lomz2)/4. -  
      (9*x*NCi*zi*lomz2)/2. -  
      2*NC*x2*zi*lomz2 +  
      (5*NCi*x2*zi*lomz2)/2. +  
      (3*NC*lz2)/2. - 3*NC*x*lz2 -  
      (3*NCi*lz2)/2. + 3*x*NCi*lz2 +  
      3*NC*x2*lz2 - 3*NCi*x2*lz2 -  
      (3*NC*omzi*lz2)/4. +  
      (3*NC*x*omzi*lz2)/2. +  
      (11*NCi*omzi*lz2)/4. -  
      (11*x*NCi*omzi*lz2)/2. -  
      (3*NC*x2*omzi*lz2)/2. +  
      3*NCi*x2*omzi*lz2 -  
      (3*NC*zi*lz2)/4. +  
      (3*NC*x*zi*lz2)/2. +  
      (7*NCi*zi*lz2)/4. -  
      (7*x*NCi*zi*lz2)/2. -  
      (3*NC*x2*zi*lz2)/2. +  
      3*NCi*x2*zi*lz2)*Tu4;
};
double C2TQ2QNS_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double sqrtxz1  = sqrt(1 - 2*z + z*z + 4*x*z);
  const double poly2    = 1 + 2*x + x*x - 4*x*z;
  const double poly2i   = 1. / poly2;
  const double poly2i2  = poly2i * poly2i;
  const double sqrtxz2  = sqrt(poly2);
  const double sqrtxz2i = 1. / sqrtxz2;
  const double sqrtxz3  = sqrt(x/z);
  const double NC   = apfel::NC;
  const double NC2  = NC * NC;
  const double NCi  = 1. / NC;
  const double NCi2 = NCi * NCi;
  const double l2  = log(2);
  const double l22 = l2 * l2;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double x2 = x * x;
  const double x3 = x * x2;
  const double x4 = x * x3;
  const double x5 = x * x4;
  const double x6 = x * x5;
  const double x7 = x * x6;
  const double xi  = 1. / x;
  const double xi2 = xi * xi;
  const double z2 = z * z;
  const double z3 = z * z2;
  const double zi  = 1. / z;
  const double omxi = 1. / ( 1 - x );
  const double opxi = 1. / ( 1 + x );
  const double omzi  = 1. / ( 1 - z );
  const double lx  = log(x);
  const double lx2 = lx * lx;
  const double lz  = log(z);
  const double lz2 = lz * lz;
  const double lomx  = log(1 - x);
  const double lomx2 = lomx * lomx;
  const double lomz  = log(1 - z);
  const double lomz2 = lomz * lomz;
  const double lopx  = log(1 + x);
  const double Li2x  = apfel::dilog(x);
  const double Li2mx = apfel::dilog(-x);
  const double Li2z  = apfel::dilog(z);
  const double omxmzi  = 1. / ( 1 - x - z );
  const double omxmzi2 = omxmzi * omxmzi;
  const double xmzi  = 1. / ( x - z );
  const double xmzi2 = xmzi * xmzi;
  const double lomxmz  = log(1 - x - z);
  const double lmopxpz = log(-1 + x + z);
  const double lxmz    = log(x - z);
  const double lmxpz   = log(-x + z);
  const double lxpz    = log(x + z);
  const double lopxz   = log(1 + x*z);
  const double lopxzi  = log(1 + x*zi);
  const double li2omxzi = apfel::dilog(1 - x*zi);
  const double lspec1   = log(1 + sqrtxz1 - z);
  const double lspec1_2 = lspec1 * lspec1;
  const double lspec2   = log(1 + sqrtxz1 + z);
  const double lspec3   = log(sqrtxz3);
  const double lspec4   = log(sqrtxz3*z);
  const double lspec5   = log(1 - sqrtxz2 + x);
  const double lspec6   = log(1 + sqrtxz2 + x);
  const double li2spec1  = apfel::dilog(0.5 - sqrtxz1/2. - z/2.);
  const double li2spec2  = apfel::dilog(0.5 - sqrtxz1/2. + z/2.);
  const double li2spec3  = apfel::dilog(0.5 - zi/2. - (sqrtxz1*zi)/2.);
  const double li2spec4  = apfel::dilog(0.5 + zi/2. - (sqrtxz1*zi)/2.);
  const double li2spec5  = apfel::dilog(0.5 - sqrtxz2/2. - x/2.);
  const double li2spec6  = apfel::dilog(0.5 + sqrtxz2/2. - x/2.);
  const double li2spec7  = apfel::dilog(0.5 - xi/2. - (sqrtxz2*xi)/2.);
  const double li2spec8  = apfel::dilog(0.5 - xi/2. + (sqrtxz2*xi)/2.);
  const double li2spec9  = apfel::dilog((1 - z)*omxi);
  const double li2spec10 = apfel::dilog(x*(1 - z)*omxi*zi);
  const double li2spec11 = apfel::dilog((1 - x)*omzi);
  const double li2spec12 = apfel::dilog((1 - x)*z*xi*omzi);
  const double li2spec13 = apfel::dilog(z*omxi);
  const double li2spec14 = apfel::dilog((1 - z)*z*omxi*xi);
  const double li2spec15 = apfel::dilog(x*z*omxi*omzi);
  const double li2spec16 = apfel::dilog((1 - x)*zi);
  const double li2spec17 = apfel::dilog((1 - x)*(1 - z)*xi*zi);
  const double li2spec18 = apfel::dilog((1 - x)*x*omzi*zi);
  const double li2spec19 = apfel::dilog(-(x*z));
  const double li2spec20 = apfel::dilog(-(x*zi));
  const double atanspec1 = atan(sqrtxz3);
  const double atanspec2 = atan(sqrtxz3*z);
  const double itani1 = InvTanInt(-sqrtxz3);
  const double itani2 = InvTanInt(sqrtxz3);
  const double itani3 = InvTanInt(sqrtxz3*z);
  const double Tr1 = (z < 1 - x ? 1 : 0);
  const double Tr2 = (z > 1 - x ? 1 : 0);
  const double Tu1 = (z < 1 - x && z < x ? 1 : 0);
  const double Tu2 = (z > 1 - x && z < x ? 1 : 0);
  const double Tu3 = (z < 1 - x && z > x ? 1 : 0);
  const double Tu4 = (z > 1 - x && z > x ? 1 : 0);
  return 12.527777777777779 + (5*NC)/3. - (5*NC*nf)/9. - (17*x)/4. + NC*x -  
   (20*z)/3. + (37*NC*z)/24. + (NC*nf*z)/3. + (403*x*z)/36. -  
   (95*NC*x*z)/24. - (8*NC*nf*x*z)/9. - (sqrtxz3*itani1)/4. +  
   (3*NC*sqrtxz3*itani1)/4. +  
   (9*sqrtxz3*x*z*itani1)/4. -  
   (3*NC*sqrtxz3*x*z*itani1)/4. +  
   (sqrtxz3*itani2)/4. - (3*NC*sqrtxz3*itani2)/4. -  
   (9*sqrtxz3*x*z*itani2)/4. +  
   (3*NC*sqrtxz3*x*z*itani2)/4. -  
   (sqrtxz3*itani3)/2. +  
   (3*NC*sqrtxz3*itani3)/2. +  
   (9*sqrtxz3*x*z*itani3)/2. -  
   (3*NC*sqrtxz3*x*z*itani3)/2. + 2*Li2mx - 2*NC*Li2mx +  
   2*x*Li2mx - 2*NC*x*Li2mx + 4*NC*z*Li2mx + 4*NC*x*z*Li2mx -  
   Li2x + NC*Li2x + 2*x*Li2x + x*z*Li2x + 2*NC*x*z*Li2x -  
   li2spec1 - x*li2spec1 +  
   NC*x*li2spec1 + 2*NC*z*li2spec1 -  
   x*z*li2spec1 + li2spec2 +  
   x*li2spec2 - NC*x*li2spec2 -  
   2*NC*z*li2spec2 + x*z*li2spec2 -  
   4*Li2z + (NC*Li2z)/2. - x*Li2z - (NC*x*Li2z)/2. - 5*x*z*Li2z +  
   3*NC*x*z*Li2z + 2*x*li2spec19 - NC*x*li2spec19 -  
   2*NC*z*li2spec19 - 2*li2spec20 + NC*li2spec20 -  
   2*x*z*li2spec20 + 2*NC*x*z*li2spec20 +  
   li2spec3 -  
   NC*li2spec3 +  
   x*li2spec3 +  
   x*z*li2spec3 -  
   2*NC*x*z*li2spec3 -  
   li2spec4 +  
   NC*li2spec4 -  
   x*li2spec4 -  
   x*z*li2spec4 +  
   2*NC*x*z*li2spec4 +  
   2*li2omxzi + 2*x*li2omxzi +  
   NC*x*li2omxzi + 4*x*z*li2omxzi -  
   2*NC*x*z*li2omxzi + 3*sqrtxz1*x*l2 -  
   2*NC*sqrtxz1*x*l2 - (sqrtxz3*atanspec1*lspec3)/2. +  
   (3*NC*sqrtxz3*atanspec1*lspec3)/2. +  
   (9*sqrtxz3*x*z*atanspec1*lspec3)/2. -  
   (3*NC*sqrtxz3*x*z*atanspec1*lspec3)/2. - lomx/6. -  
   (3*NC*lomx)/4. + (NC*nf*lomx)/3. + x*lomx -  
   (NC*x*lomx)/4. - (z*lomx)/2. + (3*NC*z*lomx)/4. -  
   (8*x*z*lomx)/3. + (NC*x*z*lomx)/4. +  
   (NC*nf*x*z*lomx)/3. + (85*lx)/12. - 2*NC*lx -  
   (2*NC*nf*lx)/3. + (x*lx)/4. - (3*NC*x*lx)/2. +  
   (3*sqrtxz1*x*lx)/2. - NC*sqrtxz1*x*lx - (11*z*lx)/4. +  
   (3*NC*z*lx)/4. + (49*x*z*lx)/12. + (3*NC*x*z*lx)/4. -  
   (2*NC*nf*x*z*lx)/3. + l2*lx + NC*l2*lx +  
   x*l2*lx - 2*NC*x*l2*lx - 4*NC*z*l2*lx +  
   x*z*l2*lx + 2*NC*x*z*l2*lx + 11*lomx*lx +  
   8*x*lomx*lx - NC*x*lomx*lx + z*lomx*lx +  
   18*x*z*lomx*lx + 2*NC*x*z*lomx*lx +  
   2*lx*lopx - 2*NC*lx*lopx + 2*x*lx*lopx -  
   2*NC*x*lx*lopx + 4*NC*z*lx*lopx +  
   4*NC*x*z*lx*lopx + (4*lomz)/3. - (3*NC*lomz)/4. +  
   (NC*nf*lomz)/3. - (NC*x*lomz)/4. - 2*z*lomz +  
   (3*NC*z*lomz)/4. - (7*x*z*lomz)/6. + (NC*x*z*lomz)/4. +  
   (NC*nf*x*z*lomz)/3. - 9*lomx*lomz -  
   3*x*lomx*lomz - z*lomx*lomz -  
   11*x*z*lomx*lomz + 14*lx*lomz -  
   NC*lx*lomz + 4*x*lx*lomz - NC*x*lx*lomz +  
   z*lx*lomz + 17*x*z*lx*lomz -  
   3*sqrtxz1*x*lspec1 + 2*NC*sqrtxz1*x*lspec1 -  
   2*l2*lspec1 - NC*l2*lspec1 -  
   2*x*l2*lspec1 + 3*NC*x*l2*lspec1 +  
   6*NC*z*l2*lspec1 - 2*x*z*l2*lspec1 -  
   2*NC*x*z*l2*lspec1 - NC*lx*lspec1 +  
   NC*x*lx*lspec1 + 2*NC*z*lx*lspec1 -  
   2*NC*x*z*lx*lspec1 + (3*lz)/4. + (3*NC*lz)/4. -  
   (25*x*lz)/4. + (NC*x*lz)/4. + (3*sqrtxz1*x*lz)/2. -  
   NC*sqrtxz1*x*lz - (5*z*lz)/4. + 2*NC*z*lz -  
   (37*x*z*lz)/4. + (NC*x*z*lz)/2. + 3*l2*lz -  
   NC*l2*lz + 3*x*l2*lz - 2*NC*x*l2*lz -  
   4*NC*z*l2*lz + 3*x*z*l2*lz - 2*NC*x*z*l2*lz -  
   4*lomx*lz - (NC*lomx*lz)/2. - 3*x*lomx*lz -  
   (NC*x*lomx*lz)/2. - 7*x*z*lomx*lz -  
   NC*x*z*lomx*lz + 7*lx*lz + NC*lx*lz +  
   2*x*lx*lz - NC*z*lx*lz + 12*x*z*lx*lz +  
   3*NC*x*z*lx*lz - 9*lomz*lz - 3*x*lomz*lz -  
   NC*x*lomz*lz - 12*x*z*lomz*lz +  
   2*NC*x*z*lomz*lz - lspec1*lz -  
   x*lspec1*lz + NC*x*lspec1*lz +  
   2*NC*z*lspec1*lz - x*z*lspec1*lz +  
   (sqrtxz3*atanspec2*lspec4)/2. -  
   (3*NC*sqrtxz3*atanspec2*lspec4)/2. -  
   (9*sqrtxz3*x*z*atanspec2*lspec4)/2. +  
   (3*NC*sqrtxz3*x*z*atanspec2*lspec4)/2. -  
   2*l2*lspec2 + NC*l2*lspec2 -  
   2*x*l2*lspec2 + NC*x*l2*lspec2 +  
   2*NC*z*l2*lspec2 - 2*x*z*l2*lspec2 +  
   2*NC*x*z*l2*lspec2 - lx*lspec2 -  
   x*lx*lspec2 + NC*x*lx*lspec2 +  
   2*NC*z*lx*lspec2 - x*z*lx*lspec2 +  
   2*lspec1*lspec2 -  
   NC*lspec1*lspec2 +  
   2*x*lspec1*lspec2 -  
   NC*x*lspec1*lspec2 -  
   2*NC*z*lspec1*lspec2 +  
   2*x*z*lspec1*lspec2 -  
   2*NC*x*z*lspec1*lspec2 -  
   2*lz*lspec2 + NC*lz*lspec2 -  
   2*x*lz*lspec2 + NC*x*lz*lspec2 +  
   2*NC*z*lz*lspec2 - 2*x*z*lz*lspec2 +  
   2*NC*x*z*lz*lspec2 - lx*lxpz +  
   (NC*lx*lxpz)/2. - x*lx*lxpz +  
   (NC*x*lx*lxpz)/2. + NC*z*lx*lxpz -  
   x*z*lx*lxpz + NC*x*z*lx*lxpz + lz*lxpz -  
   (NC*lz*lxpz)/2. + x*lz*lxpz -  
   (NC*x*lz*lxpz)/2. - NC*z*lz*lxpz +  
   x*z*lz*lxpz - NC*x*z*lz*lxpz +  
   2*x*lx*lopxz - NC*x*lx*lopxz -  
   2*NC*z*lx*lopxz + 2*x*lz*lopxz -  
   NC*x*lz*lopxz - 2*NC*z*lz*lopxz -  
   lx*lopxzi + (NC*lx*lopxzi)/2. +  
   x*lx*lopxzi - (NC*x*lx*lopxzi)/2. -  
   NC*z*lx*lopxzi - x*z*lx*lopxzi +  
   NC*x*z*lx*lopxzi + lz*lopxzi -  
   (NC*lz*lopxzi)/2. - x*lz*lopxzi +  
   (NC*x*lz*lopxzi)/2. + NC*z*lz*lopxzi +  
   x*z*lz*lopxzi - NC*x*z*lz*lopxzi -  
   (23*NCi2)/2. + (19*x*NCi2)/4. + 8*z*NCi2 -  
   (25*x*z*NCi2)/2. + (sqrtxz3*itani1*NCi2)/4. -  
   (9*sqrtxz3*x*z*itani1*NCi2)/4. -  
   (sqrtxz3*itani2*NCi2)/4. +  
   (9*sqrtxz3*x*z*itani2*NCi2)/4. +  
   (sqrtxz3*itani3*NCi2)/2. -  
   (9*sqrtxz3*x*z*itani3*NCi2)/2. -  
   2*Li2mx*NCi2 - 2*x*Li2mx*NCi2 + (Li2x*NCi2)/2. -  
   2*x*Li2x*NCi2 - (3*x*z*Li2x*NCi2)/2. +  
   li2spec1*NCi2 +  
   x*li2spec1*NCi2 +  
   x*z*li2spec1*NCi2 -  
   li2spec2*NCi2 -  
   x*li2spec2*NCi2 -  
   x*z*li2spec2*NCi2 + (5*Li2z*NCi2)/2. +  
   x*Li2z*NCi2 + (7*x*z*Li2z*NCi2)/2. -  
   2*x*li2spec19*NCi2 + 2*li2spec20*NCi2 +  
   2*x*z*li2spec20*NCi2 -  
   li2spec3*NCi2 -  
   x*li2spec3*NCi2 -  
   x*z*li2spec3*NCi2 +  
   li2spec4*NCi2 +  
   x*li2spec4*NCi2 +  
   x*z*li2spec4*NCi2 -  
   2*li2omxzi*NCi2 - 2*x*li2omxzi*NCi2 -  
   4*x*z*li2omxzi*NCi2 - 3*sqrtxz1*x*l2*NCi2 +  
   (sqrtxz3*atanspec1*lspec3*NCi2)/2. -  
   (9*sqrtxz3*x*z*atanspec1*lspec3*NCi2)/2. -  
   (3*lomx*NCi2)/4. + (3*x*lomx*NCi2)/4. +  
   (z*lomx*NCi2)/4. + (11*x*z*lomx*NCi2)/4. -  
   (9*lx*NCi2)/2. - (5*x*lx*NCi2)/2. -  
   (3*sqrtxz1*x*lx*NCi2)/2. + 3*z*lx*NCi2 -  
   5*x*z*lx*NCi2 - l2*lx*NCi2 -  
   x*l2*lx*NCi2 - x*z*l2*lx*NCi2 -  
   (9*lomx*lx*NCi2)/2. -  
   (15*x*lomx*lx*NCi2)/2. -  
   (z*lomx*lx*NCi2)/2. -  
   (23*x*z*lomx*lx*NCi2)/2. -  
   2*lx*lopx*NCi2 - 2*x*lx*lopx*NCi2 -  
   (3*lomz*NCi2)/2. + x*lomz*NCi2 +  
   z*lomz*NCi2 + 2*x*z*lomz*NCi2 +  
   (7*lomx*lomz*NCi2)/2. +  
   (5*x*lomx*lomz*NCi2)/2. +  
   (z*lomx*lomz*NCi2)/2. +  
   (11*x*z*lomx*lomz*NCi2)/2. -  
   7*lx*lomz*NCi2 - (7*x*lx*lomz*NCi2)/2. -  
   (z*lx*lomz*NCi2)/2. -  
   10*x*z*lx*lomz*NCi2 +  
   3*sqrtxz1*x*lspec1*NCi2 +  
   2*l2*lspec1*NCi2 +  
   2*x*l2*lspec1*NCi2 +  
   2*x*z*l2*lspec1*NCi2 - (11*lz*NCi2)/4. +  
   (23*x*lz*NCi2)/4. - (3*sqrtxz1*x*lz*NCi2)/2. +  
   (7*z*lz*NCi2)/4. + (25*x*z*lz*NCi2)/4. -  
   3*l2*lz*NCi2 - 3*x*l2*lz*NCi2 -  
   3*x*z*l2*lz*NCi2 + (5*lomx*lz*NCi2)/2. +  
   3*x*lomx*lz*NCi2 +  
   (11*x*z*lomx*lz*NCi2)/2. -  
   (9*lx*lz*NCi2)/2. - 2*x*lx*lz*NCi2 -  
   (19*x*z*lx*lz*NCi2)/2. + 6*lomz*lz*NCi2 +  
   3*x*lomz*lz*NCi2 + 9*x*z*lomz*lz*NCi2 +  
   lspec1*lz*NCi2 +  
   x*lspec1*lz*NCi2 +  
   x*z*lspec1*lz*NCi2 -  
   (sqrtxz3*atanspec2*lspec4*NCi2)/2. +  
   (9*sqrtxz3*x*z*atanspec2*lspec4*NCi2)/2. +  
   2*l2*lspec2*NCi2 +  
   2*x*l2*lspec2*NCi2 +  
   2*x*z*l2*lspec2*NCi2 +  
   lx*lspec2*NCi2 +  
   x*lx*lspec2*NCi2 +  
   x*z*lx*lspec2*NCi2 -  
   2*lspec1*lspec2*NCi2 -  
   2*x*lspec1*lspec2*NCi2 -  
   2*x*z*lspec1*lspec2*NCi2 +  
   2*lz*lspec2*NCi2 +  
   2*x*lz*lspec2*NCi2 +  
   2*x*z*lz*lspec2*NCi2 +  
   lx*lxpz*NCi2 + x*lx*lxpz*NCi2 +  
   x*z*lx*lxpz*NCi2 - lz*lxpz*NCi2 -  
   x*lz*lxpz*NCi2 - x*z*lz*lxpz*NCi2 -  
   2*x*lx*lopxz*NCi2 - 2*x*lz*lopxz*NCi2 +  
   lx*lopxzi*NCi2 -  
   x*lx*lopxzi*NCi2 +  
   x*z*lx*lopxzi*NCi2 -  
   lz*lopxzi*NCi2 +  
   x*lz*lopxzi*NCi2 -  
   x*z*lz*lopxzi*NCi2 - (5*NCi)/3. +  
   (5*nf*NCi)/9. - x*NCi - (37*z*NCi)/24. -  
   (nf*z*NCi)/3. + (95*x*z*NCi)/24. +  
   (8*nf*x*z*NCi)/9. - (3*sqrtxz3*itani1*NCi)/4. +  
   (3*sqrtxz3*x*z*itani1*NCi)/4. +  
   (3*sqrtxz3*itani2*NCi)/4. -  
   (3*sqrtxz3*x*z*itani2*NCi)/4. -  
   (3*sqrtxz3*itani3*NCi)/2. +  
   (3*sqrtxz3*x*z*itani3*NCi)/2. +  
   2*Li2mx*NCi + 2*x*Li2mx*NCi - 4*z*Li2mx*NCi -  
   4*x*z*Li2mx*NCi - Li2x*NCi - 2*x*z*Li2x*NCi -  
   x*li2spec1*NCi -  
   2*z*li2spec1*NCi +  
   x*li2spec2*NCi +  
   2*z*li2spec2*NCi - (Li2z*NCi)/2. +  
   (x*Li2z*NCi)/2. - 3*x*z*Li2z*NCi +  
   x*li2spec19*NCi + 2*z*li2spec19*NCi -  
   li2spec20*NCi - 2*x*z*li2spec20*NCi +  
   li2spec3*NCi +  
   2*x*z*li2spec3*NCi -  
   li2spec4*NCi -  
   2*x*z*li2spec4*NCi -  
   x*li2omxzi*NCi +  
   2*x*z*li2omxzi*NCi + 2*sqrtxz1*x*l2*NCi -  
   (3*sqrtxz3*atanspec1*lspec3*NCi)/2. +  
   (3*sqrtxz3*x*z*atanspec1*lspec3*NCi)/2. +  
   (3*lomx*NCi)/4. - (nf*lomx*NCi)/3. +  
   (x*lomx*NCi)/4. - (3*z*lomx*NCi)/4. -  
   (x*z*lomx*NCi)/4. - (nf*x*z*lomx*NCi)/3. +  
   2*lx*NCi + (2*nf*lx*NCi)/3. +  
   (3*x*lx*NCi)/2. + sqrtxz1*x*lx*NCi -  
   (3*z*lx*NCi)/4. - (3*x*z*lx*NCi)/4. +  
   (2*nf*x*z*lx*NCi)/3. - l2*lx*NCi +  
   2*x*l2*lx*NCi + 4*z*l2*lx*NCi -  
   2*x*z*l2*lx*NCi + x*lomx*lx*NCi -  
   2*x*z*lomx*lx*NCi + 2*lx*lopx*NCi +  
   2*x*lx*lopx*NCi - 4*z*lx*lopx*NCi -  
   4*x*z*lx*lopx*NCi + (3*lomz*NCi)/4. -  
   (nf*lomz*NCi)/3. + (x*lomz*NCi)/4. -  
   (3*z*lomz*NCi)/4. - (x*z*lomz*NCi)/4. -  
   (nf*x*z*lomz*NCi)/3. + lx*lomz*NCi +  
   x*lx*lomz*NCi -  
   2*sqrtxz1*x*lspec1*NCi +  
   l2*lspec1*NCi -  
   3*x*l2*lspec1*NCi -  
   6*z*l2*lspec1*NCi +  
   2*x*z*l2*lspec1*NCi +  
   lx*lspec1*NCi -  
   x*lx*lspec1*NCi -  
   2*z*lx*lspec1*NCi +  
   2*x*z*lx*lspec1*NCi - (3*lz*NCi)/4. -  
   (x*lz*NCi)/4. + sqrtxz1*x*lz*NCi -  
   2*z*lz*NCi - (x*z*lz*NCi)/2. +  
   l2*lz*NCi + 2*x*l2*lz*NCi +  
   4*z*l2*lz*NCi + 2*x*z*l2*lz*NCi +  
   (lomx*lz*NCi)/2. + (x*lomx*lz*NCi)/2. +  
   x*z*lomx*lz*NCi - lx*lz*NCi +  
   z*lx*lz*NCi - 3*x*z*lx*lz*NCi +  
   x*lomz*lz*NCi - 2*x*z*lomz*lz*NCi -  
   x*lspec1*lz*NCi -  
   2*z*lspec1*lz*NCi +  
   (3*sqrtxz3*atanspec2*lspec4*NCi)/2. -  
   (3*sqrtxz3*x*z*atanspec2*lspec4*NCi)/2. -  
   l2*lspec2*NCi -  
   x*l2*lspec2*NCi -  
   2*z*l2*lspec2*NCi -  
   2*x*z*l2*lspec2*NCi -  
   x*lx*lspec2*NCi -  
   2*z*lx*lspec2*NCi +  
   lspec1*lspec2*NCi +  
   x*lspec1*lspec2*NCi +  
   2*z*lspec1*lspec2*NCi +  
   2*x*z*lspec1*lspec2*NCi -  
   lz*lspec2*NCi -  
   x*lz*lspec2*NCi -  
   2*z*lz*lspec2*NCi -  
   2*x*z*lz*lspec2*NCi -  
   (lx*lxpz*NCi)/2. - (x*lx*lxpz*NCi)/2. -  
   z*lx*lxpz*NCi - x*z*lx*lxpz*NCi +  
   (lz*lxpz*NCi)/2. + (x*lz*lxpz*NCi)/2. +  
   z*lz*lxpz*NCi + x*z*lz*lxpz*NCi +  
   x*lx*lopxz*NCi + 2*z*lx*lopxz*NCi +  
   x*lz*lopxz*NCi + 2*z*lz*lopxz*NCi -  
   (lx*lopxzi*NCi)/2. +  
   (x*lx*lopxzi*NCi)/2. +  
   z*lx*lopxzi*NCi -  
   x*z*lx*lopxzi*NCi +  
   (lz*lopxzi*NCi)/2. -  
   (x*lz*lopxzi*NCi)/2. -  
   z*lz*lopxzi*NCi +  
   x*z*lz*lopxzi*NCi - (37*NC2)/36. -  
   (x*NC2)/2. - (4*z*NC2)/3. + (47*x*z*NC2)/36. +  
   (Li2x*NC2)/2. + (x*z*Li2x*NC2)/2. +  
   (3*Li2z*NC2)/2. + (3*x*z*Li2z*NC2)/2. +  
   (11*lomx*NC2)/12. - (7*x*lomx*NC2)/4. +  
   (z*lomx*NC2)/4. - (x*z*lomx*NC2)/12. -  
   (31*lx*NC2)/12. + (9*x*lx*NC2)/4. -  
   (z*lx*NC2)/4. + (11*x*z*lx*NC2)/12. -  
   (13*lomx*lx*NC2)/2. - (x*lomx*lx*NC2)/2. -  
   (z*lomx*lx*NC2)/2. -  
   (13*x*z*lomx*lx*NC2)/2. + (lomz*NC2)/6. -  
   x*lomz*NC2 + z*lomz*NC2 -  
   (5*x*z*lomz*NC2)/6. +  
   (11*lomx*lomz*NC2)/2. +  
   (x*lomx*lomz*NC2)/2. +  
   (z*lomx*lomz*NC2)/2. +  
   (11*x*z*lomx*lomz*NC2)/2. -  
   7*lx*lomz*NC2 - (x*lx*lomz*NC2)/2. -  
   (z*lx*lomz*NC2)/2. - 7*x*z*lx*lomz*NC2 +  
   2*lz*NC2 + (x*lz*NC2)/2. - (z*lz*NC2)/2. +  
   3*x*z*lz*NC2 + (3*lomx*lz*NC2)/2. +  
   (3*x*z*lomx*lz*NC2)/2. - (5*lx*lz*NC2)/2. -  
   (5*x*z*lx*lz*NC2)/2. + 3*lomz*lz*NC2 +  
   3*x*z*lomz*lz*NC2 + (11*pi2)/6. -  
   (5*NC*pi2)/12. + x*pi2 - (NC*x*pi2)/12. +  
   (z*pi2)/6. + (NC*z*pi2)/3. + (13*x*z*pi2)/6. -  
   (NC*x*z*pi2)/2. - (3*NCi2*pi2)/4. -  
   (11*x*NCi2*pi2)/12. - (z*NCi2*pi2)/12. -  
   (13*x*z*NCi2*pi2)/12. + (5*NCi*pi2)/12. +  
   (x*NCi*pi2)/12. - (z*NCi*pi2)/3. +  
   (x*z*NCi*pi2)/2. - (13*NC2*pi2)/12. -  
   (x*NC2*pi2)/12. - (z*NC2*pi2)/12. -  
   (13*x*z*NC2*pi2)/12. + (3*lx*poly2i2)/8. -  
   (3*x*lx*poly2i2)/8. - (3*lz*poly2i2)/8. -  
   (3*x*lz*poly2i2)/8. - (3*lx*NCi2*poly2i2)/8. +  
   (3*x*lx*NCi2*poly2i2)/8. +  
   (3*lz*NCi2*poly2i2)/8. +  
   (3*x*lz*NCi2*poly2i2)/8. + poly2i/4. -  
   (x*poly2i)/2. - (5*lx*poly2i)/8. +  
   (x*lx*poly2i)/4. + (5*lz*poly2i)/8. +  
   (x*lz*poly2i)/4. - (NCi2*poly2i)/4. +  
   (x*NCi2*poly2i)/2. + (5*lx*NCi2*poly2i)/8. -  
   (x*lx*NCi2*poly2i)/4. -  
   (5*lz*NCi2*poly2i)/8. -  
   (x*lz*NCi2*poly2i)/4. -  
   (9*li2spec5*sqrtxz2i)/4. +  
   2*NC*li2spec5*sqrtxz2i -  
   (51*x*li2spec5*sqrtxz2i)/8. +  
   NC*x*li2spec5*sqrtxz2i +  
   (9*z*li2spec5*sqrtxz2i)/2. -  
   4*NC*z*li2spec5*sqrtxz2i +  
   9*x*z*li2spec5*sqrtxz2i -  
   8*NC*x*z*li2spec5*sqrtxz2i +  
   (9*li2spec6*sqrtxz2i)/4. -  
   2*NC*li2spec6*sqrtxz2i +  
   (51*x*li2spec6*sqrtxz2i)/8. -  
   NC*x*li2spec6*sqrtxz2i -  
   (9*z*li2spec6*sqrtxz2i)/2. +  
   4*NC*z*li2spec6*sqrtxz2i -  
   9*x*z*li2spec6*sqrtxz2i +  
   8*NC*x*z*li2spec6*sqrtxz2i +  
   (9*li2spec7*sqrtxz2i)/4. -  
   2*NC*li2spec7*sqrtxz2i +  
   (51*x*li2spec7*sqrtxz2i)/ 
    8. - NC*x*li2spec7* 
    sqrtxz2i - (9*z*li2spec7* 
      sqrtxz2i)/2. + 4*NC*z* 
    li2spec7*sqrtxz2i -  
   9*x*z*li2spec7*sqrtxz2i +  
   8*NC*x*z*li2spec7* 
    sqrtxz2i - (9*li2spec8* 
      sqrtxz2i)/4. + 2*NC* 
    li2spec8*sqrtxz2i -  
   (51*x*li2spec8*sqrtxz2i)/ 
    8. + NC*x*li2spec8* 
    sqrtxz2i + (9*z*li2spec8* 
      sqrtxz2i)/2. - 4*NC*z* 
    li2spec8*sqrtxz2i +  
   9*x*z*li2spec8*sqrtxz2i -  
   8*NC*x*z*li2spec8* 
    sqrtxz2i + (9*lx*lspec5*sqrtxz2i)/4. -  
   2*NC*lx*lspec5*sqrtxz2i +  
   (51*x*lx*lspec5*sqrtxz2i)/8. -  
   NC*x*lx*lspec5*sqrtxz2i -  
   (9*z*lx*lspec5*sqrtxz2i)/2. +  
   4*NC*z*lx*lspec5*sqrtxz2i -  
   9*x*z*lx*lspec5*sqrtxz2i +  
   8*NC*x*z*lx*lspec5*sqrtxz2i -  
   (9*lx*lspec6*sqrtxz2i)/4. +  
   2*NC*lx*lspec6*sqrtxz2i -  
   (51*x*lx*lspec6*sqrtxz2i)/8. +  
   NC*x*lx*lspec6*sqrtxz2i +  
   (9*z*lx*lspec6*sqrtxz2i)/2. -  
   4*NC*z*lx*lspec6*sqrtxz2i +  
   9*x*z*lx*lspec6*sqrtxz2i -  
   8*NC*x*z*lx*lspec6*sqrtxz2i +  
   (9*li2spec5*NCi2*sqrtxz2i)/4. +  
   (51*x*li2spec5*NCi2*sqrtxz2i)/8. -  
   (9*z*li2spec5*NCi2*sqrtxz2i)/2. -  
   9*x*z*li2spec5*NCi2*sqrtxz2i -  
   (9*li2spec6*NCi2*sqrtxz2i)/4. -  
   (51*x*li2spec6*NCi2*sqrtxz2i)/8. +  
   (9*z*li2spec6*NCi2*sqrtxz2i)/2. +  
   9*x*z*li2spec6*NCi2*sqrtxz2i -  
   (9*li2spec7*NCi2* 
      sqrtxz2i)/4. - (51*x* 
      li2spec7*NCi2* 
      sqrtxz2i)/8. + (9*z* 
      li2spec7*NCi2* 
      sqrtxz2i)/2. + 9*x*z* 
    li2spec7*NCi2* 
    sqrtxz2i + (9*li2spec8* 
      NCi2*sqrtxz2i)/4. +  
   (51*x*li2spec8*NCi2* 
      sqrtxz2i)/8. - (9*z* 
      li2spec8*NCi2* 
      sqrtxz2i)/2. - 9*x*z* 
    li2spec8*NCi2* 
    sqrtxz2i - (9*lx*lspec5*NCi2* 
      sqrtxz2i)/4. - (51*x*lx*lspec5*NCi2* 
      sqrtxz2i)/8. + (9*z*lx*lspec5*NCi2* 
      sqrtxz2i)/2. + 9*x*z*lx*lspec5*NCi2* 
    sqrtxz2i + (9*lx*lspec6*NCi2* 
      sqrtxz2i)/4. + (51*x*lx*lspec6*NCi2* 
      sqrtxz2i)/8. - (9*z*lx*lspec6*NCi2* 
      sqrtxz2i)/2. - 9*x*z*lx*lspec6*NCi2* 
    sqrtxz2i - 2*li2spec5*NCi* 
    sqrtxz2i - x*li2spec5*NCi* 
    sqrtxz2i + 4*z*li2spec5*NCi* 
    sqrtxz2i + 8*x*z*li2spec5*NCi* 
    sqrtxz2i + 2*li2spec6*NCi* 
    sqrtxz2i + x*li2spec6*NCi* 
    sqrtxz2i - 4*z*li2spec6*NCi* 
    sqrtxz2i - 8*x*z*li2spec6*NCi* 
    sqrtxz2i + 2*li2spec7* 
    NCi*sqrtxz2i +  
   x*li2spec7*NCi* 
    sqrtxz2i - 4*z*li2spec7* 
    NCi*sqrtxz2i -  
   8*x*z*li2spec7*NCi* 
    sqrtxz2i - 2*li2spec8* 
    NCi*sqrtxz2i -  
   x*li2spec8*NCi* 
    sqrtxz2i + 4*z*li2spec8* 
    NCi*sqrtxz2i +  
   8*x*z*li2spec8*NCi* 
    sqrtxz2i + 2*lx*lspec5*NCi* 
    sqrtxz2i + x*lx*lspec5*NCi* 
    sqrtxz2i - 4*z*lx*lspec5*NCi* 
    sqrtxz2i - 8*x*z*lx*lspec5*NCi* 
    sqrtxz2i - 2*lx*lspec6*NCi* 
    sqrtxz2i - x*lx*lspec6*NCi* 
    sqrtxz2i + 4*z*lx*lspec6*NCi* 
    sqrtxz2i + 8*x*z*lx*lspec6*NCi* 
    sqrtxz2i - (5*omxi)/2. + (5*z*omxi)/2. -  
   (3*Li2x*omxi)/2. + (3*NC*Li2x*omxi)/2. -  
   (3*z*Li2x*omxi)/2. - 3*NC*z*Li2x*omxi +  
   2*li2spec1*omxi -  
   (3*NC*li2spec1*omxi)/2. +  
   z*li2spec1*omxi -  
   NC*z*li2spec1*omxi -  
   2*li2spec2*omxi +  
   (3*NC*li2spec2*omxi)/2. -  
   z*li2spec2*omxi +  
   NC*z*li2spec2*omxi +  
   (5*Li2z*omxi)/2. + (NC*Li2z*omxi)/2. +  
   (5*z*Li2z*omxi)/2. - NC*z*Li2z*omxi -  
   2*li2spec19*omxi + NC*li2spec19*omxi -  
   z*li2spec19*omxi + 2*NC*z*li2spec19*omxi +  
   2*li2spec20*omxi -  
   NC*li2spec20*omxi +  
   z*li2spec20*omxi -  
   2*NC*z*li2spec20*omxi -  
   2*li2spec3*omxi +  
   (NC*li2spec3*omxi)/2. -  
   z*li2spec3*omxi +  
   3*NC*z*li2spec3*omxi +  
   2*li2spec4*omxi -  
   (NC*li2spec4*omxi)/2. +  
   z*li2spec4*omxi -  
   3*NC*z*li2spec4*omxi -  
   (5*li2omxzi*omxi)/2. -  
   (NC*li2omxzi*omxi)/2. -  
   (5*z*li2omxzi*omxi)/2. +  
   NC*z*li2omxzi*omxi - 3*sqrtxz1*l2*omxi +  
   3*NC*sqrtxz1*l2*omxi + 5*z*lomx*omxi +  
   (lx*omxi)/4. - (9*NC*lx*omxi)/4. +  
   (NC*nf*lx*omxi)/2. - (3*sqrtxz1*lx*omxi)/2. +  
   (3*NC*sqrtxz1*lx*omxi)/2. - (21*z*lx*omxi)/2. +  
   (11*NC*z*lx*omxi)/4. + (NC*nf*z*lx*omxi)/2. -  
   2*l2*lx*omxi + (5*NC*l2*lx*omxi)/2. -  
   z*l2*lx*omxi - NC*z*l2*lx*omxi -  
   (25*lomx*lx*omxi)/2. +  
   (3*NC*lomx*lx*omxi)/2. -  
   (25*z*lomx*lx*omxi)/2. -  
   3*NC*z*lomx*lx*omxi + 5*z*lomz*omxi +  
   5*lomx*lomz*omxi +  
   5*z*lomx*lomz*omxi -  
   12*lx*lomz*omxi -  
   12*z*lx*lomz*omxi +  
   3*sqrtxz1*lspec1*omxi -  
   3*NC*sqrtxz1*lspec1*omxi +  
   4*l2*lspec1*omxi -  
   4*NC*l2*lspec1*omxi +  
   2*z*l2*lspec1*omxi -  
   NC*lx*lspec1*omxi +  
   2*NC*z*lx*lspec1*omxi +  
   3*lz*omxi - (3*NC*lz*omxi)/2. -  
   (3*sqrtxz1*lz*omxi)/2. +  
   (3*NC*sqrtxz1*lz*omxi)/2. + (13*z*lz*omxi)/2. -  
   (3*NC*z*lz*omxi)/2. - 6*l2*lz*omxi +  
   (7*NC*l2*lz*omxi)/2. - 3*z*l2*lz*omxi +  
   5*NC*z*l2*lz*omxi + 5*lomx*lz*omxi +  
   5*z*lomx*lz*omxi -  
   (15*lx*lz*omxi)/2. +  
   (NC*lx*lz*omxi)/2. -  
   (17*z*lx*lz*omxi)/2. -  
   (7*NC*z*lx*lz*omxi)/2. +  
   7*lomz*lz*omxi +  
   (NC*lomz*lz*omxi)/2. +  
   7*z*lomz*lz*omxi -  
   NC*z*lomz*lz*omxi +  
   2*lspec1*lz*omxi -  
   (3*NC*lspec1*lz*omxi)/2. +  
   z*lspec1*lz*omxi -  
   NC*z*lspec1*lz*omxi +  
   4*l2*lspec2*omxi -  
   2*NC*l2*lspec2*omxi +  
   2*z*l2*lspec2*omxi -  
   4*NC*z*l2*lspec2*omxi +  
   2*lx*lspec2*omxi -  
   (3*NC*lx*lspec2*omxi)/2. +  
   z*lx*lspec2*omxi -  
   NC*z*lx*lspec2*omxi -  
   4*lspec1*lspec2*omxi +  
   2*NC*lspec1*lspec2*omxi -  
   2*z*lspec1*lspec2*omxi +  
   4*NC*z*lspec1*lspec2*omxi +  
   4*lz*lspec2*omxi -  
   2*NC*lz*lspec2*omxi +  
   2*z*lz*lspec2*omxi -  
   4*NC*z*lz*lspec2*omxi +  
   2*lx*lxpz*omxi - NC*lx*lxpz*omxi +  
   z*lx*lxpz*omxi -  
   2*NC*z*lx*lxpz*omxi -  
   2*lz*lxpz*omxi + NC*lz*lxpz*omxi -  
   z*lz*lxpz*omxi +  
   2*NC*z*lz*lxpz*omxi -  
   2*lx*lopxz*omxi +  
   NC*lx*lopxz*omxi -  
   z*lx*lopxz*omxi +  
   2*NC*z*lx*lopxz*omxi -  
   2*lz*lopxz*omxi +  
   NC*lz*lopxz*omxi -  
   z*lz*lopxz*omxi +  
   2*NC*z*lz*lopxz*omxi + 2*NCi2*omxi -  
   2*z*NCi2*omxi + (3*Li2x*NCi2*omxi)/2. +  
   (3*z*Li2x*NCi2*omxi)/2. -  
   2*li2spec1*NCi2*omxi -  
   z*li2spec1*NCi2*omxi +  
   2*li2spec2*NCi2*omxi +  
   z*li2spec2*NCi2*omxi -  
   (5*Li2z*NCi2*omxi)/2. -  
   (5*z*Li2z*NCi2*omxi)/2. +  
   2*li2spec19*NCi2*omxi +  
   z*li2spec19*NCi2*omxi -  
   2*li2spec20*NCi2*omxi -  
   z*li2spec20*NCi2*omxi +  
   2*li2spec3*NCi2* 
    omxi + z*li2spec3* 
    NCi2*omxi - 2* 
    li2spec4*NCi2*omxi 
    - z*li2spec4*NCi2* 
    omxi + 3*li2omxzi*NCi2*omxi +  
   3*z*li2omxzi*NCi2*omxi +  
   3*sqrtxz1*l2*NCi2*omxi -  
   3*z*lomx*NCi2*omxi +  
   lx*NCi2*omxi +  
   (3*sqrtxz1*lx*NCi2*omxi)/2. +  
   (37*z*lx*NCi2*omxi)/4. +  
   2*l2*lx*NCi2*omxi +  
   z*l2*lx*NCi2*omxi +  
   (17*lomx*lx*NCi2*omxi)/2. +  
   (17*z*lomx*lx*NCi2*omxi)/2. -  
   3*z*lomz*NCi2*omxi -  
   3*lomx*lomz*NCi2*omxi -  
   3*z*lomx*lomz*NCi2*omxi +  
   8*lx*lomz*NCi2*omxi +  
   8*z*lx*lomz*NCi2*omxi -  
   3*sqrtxz1*lspec1*NCi2*omxi -  
   4*l2*lspec1*NCi2*omxi -  
   2*z*l2*lspec1*NCi2*omxi -  
   3*lz*NCi2*omxi +  
   (3*sqrtxz1*lz*NCi2*omxi)/2. -  
   (11*z*lz*NCi2*omxi)/2. +  
   6*l2*lz*NCi2*omxi +  
   3*z*l2*lz*NCi2*omxi -  
   4*lomx*lz*NCi2*omxi -  
   4*z*lomx*lz*NCi2*omxi +  
   (23*lx*lz*NCi2*omxi)/4. +  
   (27*z*lx*lz*NCi2*omxi)/4. -  
   (13*lomz*lz*NCi2*omxi)/2. -  
   (13*z*lomz*lz*NCi2*omxi)/2. -  
   2*lspec1*lz*NCi2*omxi -  
   z*lspec1*lz*NCi2*omxi -  
   4*l2*lspec2*NCi2*omxi -  
   2*z*l2*lspec2*NCi2*omxi -  
   2*lx*lspec2*NCi2*omxi -  
   z*lx*lspec2*NCi2*omxi +  
   4*lspec1*lspec2*NCi2*omxi +  
   2*z*lspec1*lspec2*NCi2*omxi -  
   4*lz*lspec2*NCi2*omxi -  
   2*z*lz*lspec2*NCi2*omxi -  
   2*lx*lxpz*NCi2*omxi -  
   z*lx*lxpz*NCi2*omxi +  
   2*lz*lxpz*NCi2*omxi +  
   z*lz*lxpz*NCi2*omxi +  
   2*lx*lopxz*NCi2*omxi +  
   z*lx*lopxz*NCi2*omxi +  
   2*lz*lopxz*NCi2*omxi +  
   z*lz*lopxz*NCi2*omxi -  
   (3*Li2x*NCi*omxi)/2. +  
   3*z*Li2x*NCi*omxi +  
   (3*li2spec1*NCi*omxi)/2. +  
   z*li2spec1*NCi*omxi -  
   (3*li2spec2*NCi*omxi)/2. -  
   z*li2spec2*NCi*omxi -  
   (Li2z*NCi*omxi)/2. + z*Li2z*NCi*omxi -  
   li2spec19*NCi*omxi -  
   2*z*li2spec19*NCi*omxi +  
   li2spec20*NCi*omxi +  
   2*z*li2spec20*NCi*omxi -  
   (li2spec3*NCi* 
      omxi)/2. - 3*z*li2spec3*NCi*omxi +  
   (li2spec4*NCi* 
      omxi)/2. + 3*z*li2spec4*NCi*omxi +  
   (li2omxzi*NCi*omxi)/2. -  
   z*li2omxzi*NCi*omxi -  
   3*sqrtxz1*l2*NCi*omxi +  
   (9*lx*NCi*omxi)/4. -  
   (nf*lx*NCi*omxi)/2. -  
   (3*sqrtxz1*lx*NCi*omxi)/2. -  
   (11*z*lx*NCi*omxi)/4. -  
   (nf*z*lx*NCi*omxi)/2. -  
   (5*l2*lx*NCi*omxi)/2. +  
   z*l2*lx*NCi*omxi -  
   (3*lomx*lx*NCi*omxi)/2. +  
   3*z*lomx*lx*NCi*omxi +  
   3*sqrtxz1*lspec1*NCi*omxi +  
   4*l2*lspec1*NCi*omxi +  
   lx*lspec1*NCi*omxi -  
   2*z*lx*lspec1*NCi*omxi +  
   (3*lz*NCi*omxi)/2. -  
   (3*sqrtxz1*lz*NCi*omxi)/2. +  
   (3*z*lz*NCi*omxi)/2. -  
   (7*l2*lz*NCi*omxi)/2. -  
   5*z*l2*lz*NCi*omxi -  
   (lx*lz*NCi*omxi)/2. +  
   (7*z*lx*lz*NCi*omxi)/2. -  
   (lomz*lz*NCi*omxi)/2. +  
   z*lomz*lz*NCi*omxi +  
   (3*lspec1*lz*NCi*omxi)/2. +  
   z*lspec1*lz*NCi*omxi +  
   2*l2*lspec2*NCi*omxi +  
   4*z*l2*lspec2*NCi*omxi +  
   (3*lx*lspec2*NCi*omxi)/2. +  
   z*lx*lspec2*NCi*omxi -  
   2*lspec1*lspec2*NCi*omxi -  
   4*z*lspec1*lspec2*NCi*omxi +  
   2*lz*lspec2*NCi*omxi +  
   4*z*lz*lspec2*NCi*omxi +  
   lx*lxpz*NCi*omxi +  
   2*z*lx*lxpz*NCi*omxi -  
   lz*lxpz*NCi*omxi -  
   2*z*lz*lxpz*NCi*omxi -  
   lx*lopxz*NCi*omxi -  
   2*z*lx*lopxz*NCi*omxi -  
   lz*lopxz*NCi*omxi -  
   2*z*lz*lopxz*NCi*omxi +  
   (NC2*omxi)/2. - (z*NC2*omxi)/2. -  
   (li2omxzi*NC2*omxi)/2. -  
   (z*li2omxzi*NC2*omxi)/2. -  
   2*z*lomx*NC2*omxi -  
   (5*lx*NC2*omxi)/4. +  
   (5*z*lx*NC2*omxi)/4. +  
   4*lomx*lx*NC2*omxi +  
   4*z*lomx*lx*NC2*omxi -  
   2*z*lomz*NC2*omxi -  
   2*lomx*lomz*NC2*omxi -  
   2*z*lomx*lomz*NC2*omxi +  
   4*lx*lomz*NC2*omxi +  
   4*z*lx*lomz*NC2*omxi -  
   z*lz*NC2*omxi -  
   lomx*lz*NC2*omxi -  
   z*lomx*lz*NC2*omxi +  
   (7*lx*lz*NC2*omxi)/4. +  
   (7*z*lx*lz*NC2*omxi)/4. -  
   (lomz*lz*NC2*omxi)/2. -  
   (z*lomz*lz*NC2*omxi)/2. -  
   (11*pi2*omxi)/12. - (5*NC*pi2*omxi)/12. -  
   (11*z*pi2*omxi)/12. + (5*NC*z*pi2*omxi)/6. +  
   (2*NCi2*pi2*omxi)/3. +  
   (2*z*NCi2*pi2*omxi)/3. +  
   (5*NCi*pi2*omxi)/12. -  
   (5*z*NCi*pi2*omxi)/6. +  
   (NC2*pi2*omxi)/4. +  
   (z*NC2*pi2*omxi)/4. + 2*Li2mx*xi2 -  
   NC*Li2mx*xi2 - 2*z*Li2mx*xi2 + 2*NC*z*Li2mx*xi2 -  
   2*Li2x*xi2 + NC*Li2x*xi2 + 2*z*Li2x*xi2 -  
   2*NC*z*Li2x*xi2 - li2spec19*xi2 +  
   (NC*li2spec19*xi2)/2. + z*li2spec19*xi2 -  
   NC*z*li2spec19*xi2 - li2spec20*xi2 +  
   (NC*li2spec20*xi2)/2. + z*li2spec20*xi2 -  
   NC*z*li2spec20*xi2 + 4*lx*xi2 -  
   2*NC*lx*xi2 - 4*z*lx*xi2 + 4*NC*z*lx*xi2 -  
   2*lomx*lx*xi2 + NC*lomx*lx*xi2 +  
   2*z*lomx*lx*xi2 - 2*NC*z*lomx*lx*xi2 +  
   2*lx*lopx*xi2 - NC*lx*lopx*xi2 -  
   2*z*lx*lopx*xi2 + 2*NC*z*lx*lopx*xi2 +  
   lx*lz*xi2 - (NC*lx*lz*xi2)/2. -  
   z*lx*lz*xi2 + NC*z*lx*lz*xi2 -  
   lx*lopxz*xi2 + (NC*lx*lopxz*xi2)/2. +  
   z*lx*lopxz*xi2 - NC*z*lx*lopxz*xi2 -  
   lz*lopxz*xi2 + (NC*lz*lopxz*xi2)/2. +  
   z*lz*lopxz*xi2 - NC*z*lz*lopxz*xi2 -  
   lx*lopxzi*xi2 +  
   (NC*lx*lopxzi*xi2)/2. +  
   z*lx*lopxzi*xi2 -  
   NC*z*lx*lopxzi*xi2 +  
   lz*lopxzi*xi2 -  
   (NC*lz*lopxzi*xi2)/2. -  
   z*lz*lopxzi*xi2 +  
   NC*z*lz*lopxzi*xi2 -  
   2*Li2mx*NCi2*xi2 + 2*z*Li2mx*NCi2*xi2 +  
   2*Li2x*NCi2*xi2 - 2*z*Li2x*NCi2*xi2 +  
   li2spec19*NCi2*xi2 - z*li2spec19*NCi2*xi2 +  
   li2spec20*NCi2*xi2 -  
   z*li2spec20*NCi2*xi2 -  
   4*lx*NCi2*xi2 + 4*z*lx*NCi2*xi2 +  
   2*lomx*lx*NCi2*xi2 -  
   2*z*lomx*lx*NCi2*xi2 -  
   2*lx*lopx*NCi2*xi2 +  
   2*z*lx*lopx*NCi2*xi2 -  
   lx*lz*NCi2*xi2 +  
   z*lx*lz*NCi2*xi2 +  
   lx*lopxz*NCi2*xi2 -  
   z*lx*lopxz*NCi2*xi2 +  
   lz*lopxz*NCi2*xi2 -  
   z*lz*lopxz*NCi2*xi2 +  
   lx*lopxzi*NCi2*xi2 -  
   z*lx*lopxzi*NCi2*xi2 -  
   lz*lopxzi*NCi2*xi2 +  
   z*lz*lopxzi*NCi2*xi2 +  
   Li2mx*NCi*xi2 - 2*z*Li2mx*NCi*xi2 -  
   Li2x*NCi*xi2 + 2*z*Li2x*NCi*xi2 -  
   (li2spec19*NCi*xi2)/2. +  
   z*li2spec19*NCi*xi2 -  
   (li2spec20*NCi*xi2)/2. +  
   z*li2spec20*NCi*xi2 +  
   2*lx*NCi*xi2 - 4*z*lx*NCi*xi2 -  
   lomx*lx*NCi*xi2 +  
   2*z*lomx*lx*NCi*xi2 +  
   lx*lopx*NCi*xi2 -  
   2*z*lx*lopx*NCi*xi2 +  
   (lx*lz*NCi*xi2)/2. -  
   z*lx*lz*NCi*xi2 -  
   (lx*lopxz*NCi*xi2)/2. +  
   z*lx*lopxz*NCi*xi2 -  
   (lz*lopxz*NCi*xi2)/2. +  
   z*lz*lopxz*NCi*xi2 -  
   (lx*lopxzi*NCi*xi2)/2. +  
   z*lx*lopxzi*NCi*xi2 +  
   (lz*lopxzi*NCi*xi2)/2. -  
   z*lz*lopxzi*NCi*xi2 +  
   (pi2*xi2)/3. - (NC*pi2*xi2)/6. -  
   (z*pi2*xi2)/3. + (NC*z*pi2*xi2)/3. -  
   (NCi2*pi2*xi2)/3. +  
   (z*NCi2*pi2*xi2)/3. +  
   (NCi*pi2*xi2)/6. -  
   (z*NCi*pi2*xi2)/3. + xi/4. -  
   (13*NC*xi)/9. + (9*sqrtxz3*z*itani1*xi)/4. -  
   (3*NC*sqrtxz3*z*itani1*xi)/4. -  
   (9*sqrtxz3*z*itani2*xi)/4. +  
   (3*NC*sqrtxz3*z*itani2*xi)/4. +  
   (9*sqrtxz3*z*itani3*xi)/2. -  
   (3*NC*sqrtxz3*z*itani3*xi)/2. - 2*Li2mx*xi +  
   NC*Li2mx*xi - 2*NC*z*Li2mx*xi + 2*Li2x*xi -  
   NC*Li2x*xi + 2*NC*z*Li2x*xi + li2spec19*xi -  
   (NC*li2spec19*xi)/2. + NC*z*li2spec19*xi +  
   li2spec20*xi - (NC*li2spec20*xi)/2. +  
   NC*z*li2spec20*xi +  
   (9*sqrtxz3*z*atanspec1*lspec3*xi)/2. -  
   (3*NC*sqrtxz3*z*atanspec1*lspec3*xi)/2. -  
   (2*NC*lomx*xi)/3. - (9*lx*xi)/2. +  
   2*NC*lx*xi - 4*NC*z*lx*xi +  
   2*lomx*lx*xi - NC*lomx*lx*xi +  
   2*NC*z*lomx*lx*xi - 2*lx*lopx*xi +  
   NC*lx*lopx*xi - 2*NC*z*lx*lopx*xi -  
   (2*NC*lomz*xi)/3. - (lz*xi)/2. -  
   (2*NC*lz*xi)/3. - lx*lz*xi +  
   (NC*lx*lz*xi)/2. - NC*z*lx*lz*xi -  
   (9*sqrtxz3*z*atanspec2*lspec4*xi)/2. +  
   (3*NC*sqrtxz3*z*atanspec2*lspec4*xi)/2. +  
   lx*lopxz*xi - (NC*lx*lopxz*xi)/2. +  
   NC*z*lx*lopxz*xi + lz*lopxz*xi -  
   (NC*lz*lopxz*xi)/2. +  
   NC*z*lz*lopxz*xi +  
   lx*lopxzi*xi -  
   (NC*lx*lopxzi*xi)/2. +  
   NC*z*lx*lopxzi*xi -  
   lz*lopxzi*xi +  
   (NC*lz*lopxzi*xi)/2. -  
   NC*z*lz*lopxzi*xi - (NCi2*xi)/4. -  
   (9*sqrtxz3*z*itani1*NCi2*xi)/4. +  
   (9*sqrtxz3*z*itani2*NCi2*xi)/4. -  
   (9*sqrtxz3*z*itani3*NCi2*xi)/2. +  
   2*Li2mx*NCi2*xi - 2*Li2x*NCi2*xi -  
   li2spec19*NCi2*xi -  
   li2spec20*NCi2*xi -  
   (9*sqrtxz3*z*atanspec1*lspec3*NCi2*xi)/2. +  
   (9*lx*NCi2*xi)/2. -  
   2*lomx*lx*NCi2*xi +  
   2*lx*lopx*NCi2*xi +  
   (lz*NCi2*xi)/2. + lx*lz*NCi2*xi +  
   (9*sqrtxz3*z*atanspec2*lspec4*NCi2*xi)/2. -  
   lx*lopxz*NCi2*xi -  
   lz*lopxz*NCi2*xi -  
   lx*lopxzi*NCi2*xi +  
   lz*lopxzi*NCi2*xi +  
   (13*NCi*xi)/9. +  
   (3*sqrtxz3*z*itani1*NCi*xi)/4. -  
   (3*sqrtxz3*z*itani2*NCi*xi)/4. +  
   (3*sqrtxz3*z*itani3*NCi*xi)/2. -  
   Li2mx*NCi*xi + 2*z*Li2mx*NCi*xi +  
   Li2x*NCi*xi - 2*z*Li2x*NCi*xi +  
   (li2spec19*NCi*xi)/2. -  
   z*li2spec19*NCi*xi +  
   (li2spec20*NCi*xi)/2. -  
   z*li2spec20*NCi*xi +  
   (3*sqrtxz3*z*atanspec1*lspec3*NCi*xi)/2. +  
   (2*lomx*NCi*xi)/3. - 2*lx*NCi*xi +  
   4*z*lx*NCi*xi + lomx*lx*NCi*xi -  
   2*z*lomx*lx*NCi*xi -  
   lx*lopx*NCi*xi +  
   2*z*lx*lopx*NCi*xi +  
   (2*lomz*NCi*xi)/3. +  
   (2*lz*NCi*xi)/3. -  
   (lx*lz*NCi*xi)/2. +  
   z*lx*lz*NCi*xi -  
   (3*sqrtxz3*z*atanspec2*lspec4*NCi*xi)/2. +  
   (lx*lopxz*NCi*xi)/2. -  
   z*lx*lopxz*NCi*xi +  
   (lz*lopxz*NCi*xi)/2. -  
   z*lz*lopxz*NCi*xi +  
   (lx*lopxzi*NCi*xi)/2. -  
   z*lx*lopxzi*NCi*xi -  
   (lz*lopxzi*NCi*xi)/2. +  
   z*lz*lopxzi*NCi*xi -  
   (pi2*xi)/3. + (NC*pi2*xi)/6. -  
   (NC*z*pi2*xi)/3. + (NCi2*pi2*xi)/3. -  
   (NCi*pi2*xi)/6. +  
   (z*NCi*pi2*xi)/3. -  
   (3*lx*poly2i2*xi)/8. -  
   (3*lz*poly2i2*xi)/8. +  
   (3*lx*NCi2*poly2i2*xi)/8. +  
   (3*lz*NCi2*poly2i2*xi)/8. -  
   (poly2i*xi)/4. + (7*lx*poly2i*xi)/8. +  
   (7*lz*poly2i*xi)/8. +  
   (NCi2*poly2i*xi)/4. -  
   (7*lx*NCi2*poly2i*xi)/8. -  
   (7*lz*NCi2*poly2i*xi)/8. +  
   (5*li2spec5*sqrtxz2i*xi)/16. -  
   (5*li2spec6*sqrtxz2i*xi)/16. -  
   (5*li2spec7*sqrtxz2i* 
      xi)/16. + (5*li2spec8* 
      sqrtxz2i*xi)/16. -  
   (5*lx*lspec5*sqrtxz2i*xi)/16. +  
   (5*lx*lspec6*sqrtxz2i*xi)/16. -  
   (5*li2spec5*NCi2*sqrtxz2i*xi)/ 
    16. + (5*li2spec6*NCi2*sqrtxz2i* 
      xi)/16. + (5*li2spec7* 
      NCi2*sqrtxz2i*xi)/16. -  
   (5*li2spec8*NCi2* 
      sqrtxz2i*xi)/16. +  
   (5*lx*lspec5*NCi2*sqrtxz2i*xi)/16. -  
   (5*lx*lspec6*NCi2*sqrtxz2i*xi)/16. +  
   (3*li2spec5*poly2i2*sqrtxz2i*xi)/ 
    16. - (3*li2spec6*poly2i2*sqrtxz2i* 
      xi)/16. - (3*li2spec7* 
      poly2i2*sqrtxz2i*xi)/16. +  
   (3*li2spec8*poly2i2* 
      sqrtxz2i*xi)/16. -  
   (3*lx*lspec5*poly2i2*sqrtxz2i*xi)/ 
    16. + (3*lx*lspec6*poly2i2*sqrtxz2i* 
      xi)/16. - (3*li2spec5*NCi2* 
      poly2i2*sqrtxz2i*xi)/16. +  
   (3*li2spec6*NCi2*poly2i2*sqrtxz2i* 
      xi)/16. + (3*li2spec7* 
      NCi2*poly2i2*sqrtxz2i*xi)/16. -  
   (3*li2spec8*NCi2* 
      poly2i2*sqrtxz2i*xi)/16. +  
   (3*lx*lspec5*NCi2*poly2i2*sqrtxz2i* 
      xi)/16. - (3*lx*lspec6*NCi2* 
      poly2i2*sqrtxz2i*xi)/16. -  
   (li2spec5*poly2i*sqrtxz2i*xi)/ 
    2. + (li2spec6*poly2i*sqrtxz2i* 
      xi)/2. + (li2spec7* 
      poly2i*sqrtxz2i*xi)/2. -  
   (li2spec8*poly2i* 
      sqrtxz2i*xi)/2. +  
   (lx*lspec5*poly2i*sqrtxz2i*xi)/2. -  
   (lx*lspec6*poly2i*sqrtxz2i*xi)/2. +  
   (li2spec5*NCi2*poly2i*sqrtxz2i* 
      xi)/2. - (li2spec6*NCi2*poly2i* 
      sqrtxz2i*xi)/2. -  
   (li2spec7*NCi2*poly2i* 
      sqrtxz2i*xi)/2. +  
   (li2spec8*NCi2*poly2i* 
      sqrtxz2i*xi)/2. -  
   (lx*lspec5*NCi2*poly2i*sqrtxz2i* 
      xi)/2. + (lx*lspec6*NCi2*poly2i* 
      sqrtxz2i*xi)/2. - x2/4. + (22*NC*x2)/9. +  
   (2*NC*lomx*x2)/3. + (5*lx*x2)/2. -  
   3*NC*lx*x2 + (2*NC*lomz*x2)/3. -  
   (5*lz*x2)/2. + (5*NC*lz*x2)/3. -  
   (NCi2*x2)/4. - (5*lx*NCi2*x2)/2. +  
   (5*lz*NCi2*x2)/2. - (22*NCi*x2)/9. -  
   (2*lomx*NCi*x2)/3. + 3*lx*NCi*x2 -  
   (2*lomz*NCi*x2)/3. -  
   (5*lz*NCi*x2)/3. + (NC2*x2)/2. +  
   (3*lx*poly2i2*x2)/8. -  
   (3*lz*poly2i2*x2)/8. -  
   (3*lx*NCi2*poly2i2*x2)/8. +  
   (3*lz*NCi2*poly2i2*x2)/8. +  
   (poly2i*x2)/2. + (lx*poly2i*x2)/4. -  
   (lz*poly2i*x2)/4. -  
   (NCi2*poly2i*x2)/2. -  
   (lx*NCi2*poly2i*x2)/4. +  
   (lz*NCi2*poly2i*x2)/4. -  
   (9*li2spec5*sqrtxz2i*x2)/4. +  
   2*NC*li2spec5*sqrtxz2i*x2 +  
   (9*z*li2spec5*sqrtxz2i*x2)/2. -  
   4*NC*z*li2spec5*sqrtxz2i*x2 +  
   (9*li2spec6*sqrtxz2i*x2)/4. -  
   2*NC*li2spec6*sqrtxz2i*x2 -  
   (9*z*li2spec6*sqrtxz2i*x2)/2. +  
   4*NC*z*li2spec6*sqrtxz2i*x2 +  
   (9*li2spec7*sqrtxz2i* 
      x2)/4. - 2*NC*li2spec7* 
    sqrtxz2i*x2 - (9*z* 
      li2spec7*sqrtxz2i* 
      x2)/2. + 4*NC*z*li2spec7* 
    sqrtxz2i*x2 - (9* 
      li2spec8*sqrtxz2i* 
      x2)/4. + 2*NC*li2spec8* 
    sqrtxz2i*x2 + (9*z* 
      li2spec8*sqrtxz2i* 
      x2)/2. - 4*NC*z*li2spec8* 
    sqrtxz2i*x2 + (9*lx*lspec5*sqrtxz2i* 
      x2)/4. - 2*NC*lx*lspec5*sqrtxz2i* 
    x2 - (9*z*lx*lspec5*sqrtxz2i*x2)/ 
    2. + 4*NC*z*lx*lspec5*sqrtxz2i*x2 -  
   (9*lx*lspec6*sqrtxz2i*x2)/4. +  
   2*NC*lx*lspec6*sqrtxz2i*x2 +  
   (9*z*lx*lspec6*sqrtxz2i*x2)/2. -  
   4*NC*z*lx*lspec6*sqrtxz2i*x2 +  
   (9*li2spec5*NCi2*sqrtxz2i*x2)/4. -  
   (9*z*li2spec5*NCi2*sqrtxz2i*x2)/ 
    2. - (9*li2spec6*NCi2*sqrtxz2i*x2)/ 
    4. + (9*z*li2spec6*NCi2*sqrtxz2i* 
      x2)/2. - (9*li2spec7* 
      NCi2*sqrtxz2i*x2)/4. +  
   (9*z*li2spec7*NCi2* 
      sqrtxz2i*x2)/2. +  
   (9*li2spec8*NCi2* 
      sqrtxz2i*x2)/4. -  
   (9*z*li2spec8*NCi2* 
      sqrtxz2i*x2)/2. -  
   (9*lx*lspec5*NCi2*sqrtxz2i*x2)/4. +  
   (9*z*lx*lspec5*NCi2*sqrtxz2i*x2)/2. +  
   (9*lx*lspec6*NCi2*sqrtxz2i*x2)/4. -  
   (9*z*lx*lspec6*NCi2*sqrtxz2i*x2)/2. -  
   2*li2spec5*NCi*sqrtxz2i*x2 +  
   4*z*li2spec5*NCi*sqrtxz2i*x2 +  
   2*li2spec6*NCi*sqrtxz2i*x2 -  
   4*z*li2spec6*NCi*sqrtxz2i*x2 +  
   2*li2spec7*NCi* 
    sqrtxz2i*x2 - 4*z* 
    li2spec7*NCi* 
    sqrtxz2i*x2 - 2* 
    li2spec8*NCi* 
    sqrtxz2i*x2 + 4*z* 
    li2spec8*NCi* 
    sqrtxz2i*x2 + 2*lx*lspec5*NCi* 
    sqrtxz2i*x2 - 4*z*lx*lspec5*NCi* 
    sqrtxz2i*x2 - 2*lx*lspec6*NCi* 
    sqrtxz2i*x2 + 4*z*lx*lspec6*NCi* 
    sqrtxz2i*x2 + (3*lx*poly2i2*x3)/8. +  
   (3*lz*poly2i2*x3)/8. -  
   (3*lx*NCi2*poly2i2*x3)/8. -  
   (3*lz*NCi2*poly2i2*x3)/8. -  
   (poly2i*x3)/4. - (5*lx*poly2i*x3)/8. -  
   (5*lz*poly2i*x3)/8. +  
   (NCi2*poly2i*x3)/4. +  
   (5*lx*NCi2*poly2i*x3)/8. +  
   (5*lz*NCi2*poly2i*x3)/8. +  
   (5*li2spec5*sqrtxz2i*x3)/16. -  
   (5*li2spec6*sqrtxz2i*x3)/16. -  
   (5*li2spec7*sqrtxz2i* 
      x3)/16. + (5*li2spec8* 
      sqrtxz2i*x3)/16. -  
   (5*lx*lspec5*sqrtxz2i*x3)/16. +  
   (5*lx*lspec6*sqrtxz2i*x3)/16. -  
   (5*li2spec5*NCi2*sqrtxz2i*x3)/16. +  
   (5*li2spec6*NCi2*sqrtxz2i*x3)/16. +  
   (5*li2spec7*NCi2* 
      sqrtxz2i*x3)/16. -  
   (5*li2spec8*NCi2* 
      sqrtxz2i*x3)/16. +  
   (5*lx*lspec5*NCi2*sqrtxz2i*x3)/16. -  
   (5*lx*lspec6*NCi2*sqrtxz2i*x3)/16. -  
   (3*li2spec5*poly2i2*sqrtxz2i*x3)/ 
    8. + (3*li2spec6*poly2i2*sqrtxz2i* 
      x3)/8. + (3*li2spec7* 
      poly2i2*sqrtxz2i*x3)/8. -  
   (3*li2spec8*poly2i2* 
      sqrtxz2i*x3)/8. +  
   (3*lx*lspec5*poly2i2*sqrtxz2i*x3)/ 
    8. - (3*lx*lspec6*poly2i2*sqrtxz2i* 
      x3)/8. + (3*li2spec5*NCi2*poly2i2* 
      sqrtxz2i*x3)/8. -  
   (3*li2spec6*NCi2*poly2i2*sqrtxz2i* 
      x3)/8. - (3*li2spec7* 
      NCi2*poly2i2*sqrtxz2i*x3)/8. +  
   (3*li2spec8*NCi2* 
      poly2i2*sqrtxz2i*x3)/8. -  
   (3*lx*lspec5*NCi2*poly2i2*sqrtxz2i* 
      x3)/8. + (3*lx*lspec6*NCi2*poly2i2* 
      sqrtxz2i*x3)/8. - (3*lx*poly2i2*x4)/8. +  
   (3*lz*poly2i2*x4)/8. +  
   (3*lx*NCi2*poly2i2*x4)/8. -  
   (3*lz*NCi2*poly2i2*x4)/8. +  
   (poly2i*x4)/4. + (7*lx*poly2i*x4)/8. -  
   (7*lz*poly2i*x4)/8. -  
   (NCi2*poly2i*x4)/4. -  
   (7*lx*NCi2*poly2i*x4)/8. +  
   (7*lz*NCi2*poly2i*x4)/8. +  
   (3*lx*poly2i2*x5)/8. +  
   (3*lz*poly2i2*x5)/8. -  
   (3*lx*NCi2*poly2i2*x5)/8. -  
   (3*lz*NCi2*poly2i2*x5)/8. -  
   (li2spec5*poly2i*sqrtxz2i*x5)/2. +  
   (li2spec6*poly2i*sqrtxz2i*x5)/2. +  
   (li2spec7*poly2i* 
      sqrtxz2i*x5)/2. -  
   (li2spec8*poly2i* 
      sqrtxz2i*x5)/2. +  
   (lx*lspec5*poly2i*sqrtxz2i*x5)/2. -  
   (lx*lspec6*poly2i*sqrtxz2i*x5)/2. +  
   (li2spec5*NCi2*poly2i*sqrtxz2i* 
      x5)/2. - (li2spec6*NCi2*poly2i* 
      sqrtxz2i*x5)/2. -  
   (li2spec7*NCi2*poly2i* 
      sqrtxz2i*x5)/2. +  
   (li2spec8*NCi2*poly2i* 
      sqrtxz2i*x5)/2. -  
   (lx*lspec5*NCi2*poly2i*sqrtxz2i* 
      x5)/2. + (lx*lspec6*NCi2*poly2i* 
      sqrtxz2i*x5)/2. - (3*lx*poly2i2*x6)/8. +  
   (3*lz*poly2i2*x6)/8. +  
   (3*lx*NCi2*poly2i2*x6)/8. -  
   (3*lz*NCi2*poly2i2*x6)/8. +  
   (3*li2spec5*poly2i2*sqrtxz2i*x7)/ 
    16. - (3*li2spec6*poly2i2*sqrtxz2i* 
      x7)/16. - (3*li2spec7* 
      poly2i2*sqrtxz2i*x7)/16. +  
   (3*li2spec8*poly2i2* 
      sqrtxz2i*x7)/16. -  
   (3*lx*lspec5*poly2i2*sqrtxz2i*x7)/ 
    16. + (3*lx*lspec6*poly2i2*sqrtxz2i* 
      x7)/16. - (3*li2spec5*NCi2* 
      poly2i2*sqrtxz2i*x7)/16. +  
   (3*li2spec6*NCi2*poly2i2*sqrtxz2i* 
      x7)/16. + (3*li2spec7* 
      NCi2*poly2i2*sqrtxz2i*x7)/16. -  
   (3*li2spec8*NCi2* 
      poly2i2*sqrtxz2i*x7)/16. +  
   (3*lx*lspec5*NCi2*poly2i2*sqrtxz2i* 
      x7)/16. - (3*lx*lspec6*NCi2*poly2i2* 
      sqrtxz2i*x7)/16. + 2*Li2mx*opxi -  
   NC*Li2mx*opxi + 2*NC*z*Li2mx*opxi -  
   2*Li2x*opxi + NC*Li2x*opxi -  
   2*NC*z*Li2x*opxi + li2spec19*opxi -  
   (NC*li2spec19*opxi)/2. - z*li2spec19*opxi +  
   NC*z*li2spec19*opxi + li2spec20*opxi -  
   (NC*li2spec20*opxi)/2. -  
   z*li2spec20*opxi +  
   NC*z*li2spec20*opxi + 4*lx*opxi -  
   2*NC*lx*opxi + 4*NC*z*lx*opxi -  
   2*lomx*lx*opxi + NC*lomx*lx*opxi -  
   2*NC*z*lomx*lx*opxi +  
   2*lx*lopx*opxi - NC*lx*lopx*opxi +  
   2*NC*z*lx*lopx*opxi - lx*lz*opxi +  
   (NC*lx*lz*opxi)/2. + z*lx*lz*opxi -  
   NC*z*lx*lz*opxi + lx*lopxz*opxi -  
   (NC*lx*lopxz*opxi)/2. -  
   z*lx*lopxz*opxi +  
   NC*z*lx*lopxz*opxi +  
   lz*lopxz*opxi -  
   (NC*lz*lopxz*opxi)/2. -  
   z*lz*lopxz*opxi +  
   NC*z*lz*lopxz*opxi +  
   lx*lopxzi*opxi -  
   (NC*lx*lopxzi*opxi)/2. -  
   z*lx*lopxzi*opxi +  
   NC*z*lx*lopxzi*opxi -  
   lz*lopxzi*opxi +  
   (NC*lz*lopxzi*opxi)/2. +  
   z*lz*lopxzi*opxi -  
   NC*z*lz*lopxzi*opxi -  
   2*Li2mx*NCi2*opxi + 2*Li2x*NCi2*opxi -  
   li2spec19*NCi2*opxi +  
   z*li2spec19*NCi2*opxi -  
   li2spec20*NCi2*opxi +  
   z*li2spec20*NCi2*opxi -  
   4*lx*NCi2*opxi +  
   2*lomx*lx*NCi2*opxi -  
   2*lx*lopx*NCi2*opxi +  
   lx*lz*NCi2*opxi -  
   z*lx*lz*NCi2*opxi -  
   lx*lopxz*NCi2*opxi +  
   z*lx*lopxz*NCi2*opxi -  
   lz*lopxz*NCi2*opxi +  
   z*lz*lopxz*NCi2*opxi -  
   lx*lopxzi*NCi2*opxi +  
   z*lx*lopxzi*NCi2*opxi +  
   lz*lopxzi*NCi2*opxi -  
   z*lz*lopxzi*NCi2*opxi +  
   Li2mx*NCi*opxi - 2*z*Li2mx*NCi*opxi -  
   Li2x*NCi*opxi + 2*z*Li2x*NCi*opxi +  
   (li2spec19*NCi*opxi)/2. -  
   z*li2spec19*NCi*opxi +  
   (li2spec20*NCi*opxi)/2. -  
   z*li2spec20*NCi*opxi +  
   2*lx*NCi*opxi - 4*z*lx*NCi*opxi -  
   lomx*lx*NCi*opxi +  
   2*z*lomx*lx*NCi*opxi +  
   lx*lopx*NCi*opxi -  
   2*z*lx*lopx*NCi*opxi -  
   (lx*lz*NCi*opxi)/2. +  
   z*lx*lz*NCi*opxi +  
   (lx*lopxz*NCi*opxi)/2. -  
   z*lx*lopxz*NCi*opxi +  
   (lz*lopxz*NCi*opxi)/2. -  
   z*lz*lopxz*NCi*opxi +  
   (lx*lopxzi*NCi*opxi)/2. -  
   z*lx*lopxzi*NCi*opxi -  
   (lz*lopxzi*NCi*opxi)/2. +  
   z*lz*lopxzi*NCi*opxi +  
   (2*pi2*opxi)/3. - (NC*pi2*opxi)/3. -  
   (z*pi2*opxi)/6. + (2*NC*z*pi2*opxi)/3. -  
   (2*NCi2*pi2*opxi)/3. +  
   (z*NCi2*pi2*opxi)/6. +  
   (NCi*pi2*opxi)/3. -  
   (2*z*NCi*pi2*opxi)/3. -  
   2*Li2mx*xi2*opxi + NC*Li2mx*xi2*opxi +  
   2*z*Li2mx*xi2*opxi -  
   2*NC*z*Li2mx*xi2*opxi +  
   2*Li2x*xi2*opxi - NC*Li2x*xi2*opxi -  
   2*z*Li2x*xi2*opxi +  
   2*NC*z*Li2x*xi2*opxi +  
   li2spec19*xi2*opxi -  
   (NC*li2spec19*xi2*opxi)/2. -  
   z*li2spec19*xi2*opxi +  
   NC*z*li2spec19*xi2*opxi +  
   li2spec20*xi2*opxi -  
   (NC*li2spec20*xi2*opxi)/2. -  
   z*li2spec20*xi2*opxi +  
   NC*z*li2spec20*xi2*opxi -  
   4*lx*xi2*opxi + 2*NC*lx*xi2*opxi +  
   4*z*lx*xi2*opxi -  
   4*NC*z*lx*xi2*opxi +  
   2*lomx*lx*xi2*opxi -  
   NC*lomx*lx*xi2*opxi -  
   2*z*lomx*lx*xi2*opxi +  
   2*NC*z*lomx*lx*xi2*opxi -  
   2*lx*lopx*xi2*opxi +  
   NC*lx*lopx*xi2*opxi +  
   2*z*lx*lopx*xi2*opxi -  
   2*NC*z*lx*lopx*xi2*opxi -  
   lx*lz*xi2*opxi +  
   (NC*lx*lz*xi2*opxi)/2. +  
   z*lx*lz*xi2*opxi -  
   NC*z*lx*lz*xi2*opxi +  
   lx*lopxz*xi2*opxi -  
   (NC*lx*lopxz*xi2*opxi)/2. -  
   z*lx*lopxz*xi2*opxi +  
   NC*z*lx*lopxz*xi2*opxi +  
   lz*lopxz*xi2*opxi -  
   (NC*lz*lopxz*xi2*opxi)/2. -  
   z*lz*lopxz*xi2*opxi +  
   NC*z*lz*lopxz*xi2*opxi +  
   lx*lopxzi*xi2*opxi -  
   (NC*lx*lopxzi*xi2*opxi)/2. -  
   z*lx*lopxzi*xi2*opxi +  
   NC*z*lx*lopxzi*xi2*opxi -  
   lz*lopxzi*xi2*opxi +  
   (NC*lz*lopxzi*xi2*opxi)/2. +  
   z*lz*lopxzi*xi2*opxi -  
   NC*z*lz*lopxzi*xi2*opxi +  
   2*Li2mx*NCi2*xi2*opxi -  
   2*z*Li2mx*NCi2*xi2*opxi -  
   2*Li2x*NCi2*xi2*opxi +  
   2*z*Li2x*NCi2*xi2*opxi -  
   li2spec19*NCi2*xi2*opxi +  
   z*li2spec19*NCi2*xi2*opxi -  
   li2spec20*NCi2*xi2*opxi +  
   z*li2spec20*NCi2*xi2*opxi +  
   4*lx*NCi2*xi2*opxi -  
   4*z*lx*NCi2*xi2*opxi -  
   2*lomx*lx*NCi2*xi2*opxi +  
   2*z*lomx*lx*NCi2*xi2*opxi +  
   2*lx*lopx*NCi2*xi2*opxi -  
   2*z*lx*lopx*NCi2*xi2*opxi +  
   lx*lz*NCi2*xi2*opxi -  
   z*lx*lz*NCi2*xi2*opxi -  
   lx*lopxz*NCi2*xi2*opxi +  
   z*lx*lopxz*NCi2*xi2*opxi -  
   lz*lopxz*NCi2*xi2*opxi +  
   z*lz*lopxz*NCi2*xi2*opxi -  
   lx*lopxzi*NCi2*xi2*opxi +  
   z*lx*lopxzi*NCi2*xi2*opxi +  
   lz*lopxzi*NCi2*xi2*opxi -  
   z*lz*lopxzi*NCi2*xi2*opxi -  
   Li2mx*NCi*xi2*opxi +  
   2*z*Li2mx*NCi*xi2*opxi +  
   Li2x*NCi*xi2*opxi -  
   2*z*Li2x*NCi*xi2*opxi +  
   (li2spec19*NCi*xi2*opxi)/2. -  
   z*li2spec19*NCi*xi2*opxi +  
   (li2spec20*NCi*xi2*opxi)/2. -  
   z*li2spec20*NCi*xi2*opxi -  
   2*lx*NCi*xi2*opxi +  
   4*z*lx*NCi*xi2*opxi +  
   lomx*lx*NCi*xi2*opxi -  
   2*z*lomx*lx*NCi*xi2*opxi -  
   lx*lopx*NCi*xi2*opxi +  
   2*z*lx*lopx*NCi*xi2*opxi -  
   (lx*lz*NCi*xi2*opxi)/2. +  
   z*lx*lz*NCi*xi2*opxi +  
   (lx*lopxz*NCi*xi2*opxi)/2. -  
   z*lx*lopxz*NCi*xi2*opxi +  
   (lz*lopxz*NCi*xi2*opxi)/2. -  
   z*lz*lopxz*NCi*xi2*opxi +  
   (lx*lopxzi*NCi*xi2*opxi)/2. -  
   z*lx*lopxzi*NCi*xi2*opxi -  
   (lz*lopxzi*NCi*xi2*opxi)/2. +  
   z*lz*lopxzi*NCi*xi2*opxi -  
   (pi2*xi2*opxi)/3. +  
   (NC*pi2*xi2*opxi)/6. +  
   (z*pi2*xi2*opxi)/3. -  
   (NC*z*pi2*xi2*opxi)/3. +  
   (NCi2*pi2*xi2*opxi)/3. -  
   (z*NCi2*pi2*xi2*opxi)/3. -  
   (NCi*pi2*xi2*opxi)/6. +  
   (z*NCi*pi2*xi2*opxi)/3. +  
   2*z*Li2mx*xi*opxi - 2*z*Li2x*xi*opxi -  
   z*li2spec19*xi*opxi -  
   z*li2spec20*xi*opxi +  
   4*z*lx*xi*opxi -  
   2*z*lomx*lx*xi*opxi +  
   2*z*lx*lopx*xi*opxi +  
   z*lx*lz*xi*opxi -  
   z*lx*lopxz*xi*opxi -  
   z*lz*lopxz*xi*opxi -  
   z*lx*lopxzi*xi*opxi +  
   z*lz*lopxzi*xi*opxi -  
   2*z*Li2mx*NCi2*xi*opxi +  
   2*z*Li2x*NCi2*xi*opxi +  
   z*li2spec19*NCi2*xi*opxi +  
   z*li2spec20*NCi2*xi*opxi -  
   4*z*lx*NCi2*xi*opxi +  
   2*z*lomx*lx*NCi2*xi*opxi -  
   2*z*lx*lopx*NCi2*xi*opxi -  
   z*lx*lz*NCi2*xi*opxi +  
   z*lx*lopxz*NCi2*xi*opxi +  
   z*lz*lopxz*NCi2*xi*opxi +  
   z*lx*lopxzi*NCi2*xi*opxi -  
   z*lz*lopxzi*NCi2*xi*opxi +  
   (z*pi2*xi*opxi)/3. -  
   (z*NCi2*pi2*xi*opxi)/3. -  
   (9*omzi)/2. - (3*x*omzi)/2. - 2*Li2x*omzi -  
   4*x*Li2x*omzi + (Li2z*omzi)/2. +  
   (3*x*Li2z*omzi)/2. + li2omxzi*omzi +  
   x*li2omxzi*omzi + 5*x*lomx*omzi +  
   8*lx*omzi + 7*x*lx*omzi -  
   6*lomx*lx*omzi - 13*x*lomx*lx*omzi +  
   5*x*lomz*omzi + 3*lomx*lomz*omzi +  
   5*x*lomx*lomz*omzi -  
   6*lx*lomz*omzi - 9*x*lx*lomz*omzi -  
   (9*lz*omzi)/2. - (NC*lz*omzi)/2. -  
   (5*x*lz*omzi)/2. + (NC*x*lz*omzi)/2. +  
   (5*lomx*lz*omzi)/2. +  
   (11*x*lomx*lz*omzi)/2. -  
   2*lx*lz*omzi - 4*x*lx*lz*omzi +  
   (7*lomz*lz*omzi)/2. +  
   (13*x*lomz*lz*omzi)/2. + 4*NCi2*omzi +  
   2*x*NCi2*omzi + 2*Li2x*NCi2*omzi +  
   4*x*Li2x*NCi2*omzi -  
   (Li2z*NCi2*omzi)/2. -  
   (3*x*Li2z*NCi2*omzi)/2. -  
   (li2omxzi*NCi2*omzi)/2. -  
   (x*li2omxzi*NCi2*omzi)/2. -  
   3*x*lomx*NCi2*omzi -  
   8*lx*NCi2*omzi - 10*x*lx*NCi2*omzi +  
   (7*lomx*lx*NCi2*omzi)/2. +  
   (21*x*lomx*lx*NCi2*omzi)/2. -  
   3*x*lomz*NCi2*omzi -  
   lomx*lomz*NCi2*omzi -  
   3*x*lomx*lomz*NCi2*omzi +  
   3*lx*lomz*NCi2*omzi +  
   6*x*lx*lomz*NCi2*omzi +  
   4*lz*NCi2*omzi + 4*x*lz*NCi2*omzi -  
   (3*lomx*lz*NCi2*omzi)/2. -  
   (9*x*lomx*lz*NCi2*omzi)/2. +  
   2*x*lx*lz*NCi2*omzi -  
   (5*lomz*lz*NCi2*omzi)/2. -  
   (11*x*lomz*lz*NCi2*omzi)/2. +  
   (lz*NCi*omzi)/2. -  
   (x*lz*NCi*omzi)/2. + (NC2*omzi)/2. -  
   (x*NC2*omzi)/2. -  
   (li2omxzi*NC2*omzi)/2. -  
   (x*li2omxzi*NC2*omzi)/2. -  
   2*x*lomx*NC2*omzi +  
   3*x*lx*NC2*omzi +  
   (5*lomx*lx*NC2*omzi)/2. +  
   (5*x*lomx*lx*NC2*omzi)/2. -  
   2*x*lomz*NC2*omzi -  
   2*lomx*lomz*NC2*omzi -  
   2*x*lomx*lomz*NC2*omzi +  
   3*lx*lomz*NC2*omzi +  
   3*x*lx*lomz*NC2*omzi +  
   (lz*NC2*omzi)/2. -  
   (3*x*lz*NC2*omzi)/2. -  
   lomx*lz*NC2*omzi -  
   x*lomx*lz*NC2*omzi +  
   2*lx*lz*NC2*omzi +  
   2*x*lx*lz*NC2*omzi -  
   lomz*lz*NC2*omzi -  
   x*lomz*lz*NC2*omzi -  
   (pi2*omzi)/12. - (7*x*pi2*omzi)/12. -  
   (NCi2*pi2*omzi)/6. +  
   (x*NCi2*pi2*omzi)/3. +  
   (NC2*pi2*omzi)/4. +  
   (x*NC2*pi2*omzi)/4. -  
   (3*li2spec5*sqrtxz2i*omzi)/2. +  
   x*li2spec5*sqrtxz2i*omzi +  
   (3*li2spec6*sqrtxz2i*omzi)/2. -  
   x*li2spec6*sqrtxz2i*omzi +  
   (3*li2spec7*sqrtxz2i* 
      omzi)/2. - x*li2spec7* 
    sqrtxz2i*omzi -  
   (3*li2spec8*sqrtxz2i* 
      omzi)/2. + x*li2spec8* 
    sqrtxz2i*omzi +  
   (3*lx*lspec5*sqrtxz2i*omzi)/2. -  
   x*lx*lspec5*sqrtxz2i*omzi -  
   (3*lx*lspec6*sqrtxz2i*omzi)/2. +  
   x*lx*lspec6*sqrtxz2i*omzi +  
   (3*li2spec5*NCi2*sqrtxz2i*omzi)/ 
    2. - x*li2spec5*NCi2*sqrtxz2i* 
    omzi - (3*li2spec6*NCi2* 
      sqrtxz2i*omzi)/2. +  
   x*li2spec6*NCi2*sqrtxz2i*omzi -  
   (3*li2spec7*NCi2* 
      sqrtxz2i*omzi)/2. +  
   x*li2spec7*NCi2* 
    sqrtxz2i*omzi +  
   (3*li2spec8*NCi2* 
      sqrtxz2i*omzi)/2. -  
   x*li2spec8*NCi2* 
    sqrtxz2i*omzi -  
   (3*lx*lspec5*NCi2*sqrtxz2i*omzi)/ 
    2. + x*lx*lspec5*NCi2*sqrtxz2i* 
    omzi + (3*lx*lspec6*NCi2*sqrtxz2i* 
      omzi)/2. - x*lx*lspec6*NCi2* 
    sqrtxz2i*omzi + 6*omxi*omzi +  
   6*Li2x*omxi*omzi -  
   2*Li2z*omxi*omzi -  
   2*li2omxzi*omxi*omzi -  
   5*lomx*omxi*omzi -  
   15*lx*omxi*omzi +  
   19*lomx*lx*omxi*omzi -  
   5*lomz*omxi*omzi -  
   8*lomx*lomz*omxi*omzi +  
   15*lx*lomz*omxi*omzi +  
   (17*lz*omxi*omzi)/2. -  
   8*lomx*lz*omxi*omzi +  
   6*lx*lz*omxi*omzi -  
   8*lomz*lz*omxi*omzi -  
   6*NCi2*omxi*omzi -  
   6*Li2x*NCi2*omxi*omzi +  
   2*Li2z*NCi2*omxi*omzi +  
   li2omxzi*NCi2*omxi*omzi +  
   3*lomx*NCi2*omxi*omzi +  
   18*lx*NCi2*omxi*omzi -  
   14*lomx*lx*NCi2*omxi*omzi +  
   3*lomz*NCi2*omxi*omzi +  
   4*lomx*lomz*NCi2*omxi*omzi -  
   9*lx*lomz*NCi2*omxi*omzi -  
   (35*lz*NCi2*omxi*omzi)/4. +  
   6*lomx*lz*NCi2*omxi*omzi -  
   2*lx*lz*NCi2*omxi*omzi +  
   7*lomz*lz*NCi2*omxi*omzi +  
   li2omxzi*NC2*omxi*omzi +  
   2*lomx*NC2*omxi*omzi -  
   3*lx*NC2*omxi*omzi -  
   5*lomx*lx*NC2*omxi*omzi +  
   2*lomz*NC2*omxi*omzi +  
   4*lomx*lomz*NC2*omxi*omzi -  
   6*lx*lomz*NC2*omxi*omzi +  
   (lz*NC2*omxi*omzi)/4. +  
   2*lomx*lz*NC2*omxi*omzi -  
   4*lx*lz*NC2*omxi*omzi +  
   lomz*lz*NC2*omxi*omzi +  
   (2*pi2*omxi*omzi)/3. -  
   (NCi2*pi2*omxi*omzi)/6. -  
   (NC2*pi2*omxi*omzi)/2. +  
   (NC*lz*xi*omzi)/3. -  
   (lz*NCi*xi*omzi)/3. -  
   (NC*lz*x2*omzi)/3. +  
   (lz*NCi*x2*omzi)/3. -  
   (3*li2spec5*sqrtxz2i*x2*omzi)/ 
    2. + (3*li2spec6*sqrtxz2i*x2* 
      omzi)/2. + (3*li2spec7* 
      sqrtxz2i*x2*omzi)/2. -  
   (3*li2spec8*sqrtxz2i* 
      x2*omzi)/2. +  
   (3*lx*lspec5*sqrtxz2i*x2*omzi)/ 
    2. - (3*lx*lspec6*sqrtxz2i*x2* 
      omzi)/2. + (3*li2spec5*NCi2* 
      sqrtxz2i*x2*omzi)/2. -  
   (3*li2spec6*NCi2*sqrtxz2i*x2* 
      omzi)/2. - (3*li2spec7* 
      NCi2*sqrtxz2i*x2*omzi)/2. +  
   (3*li2spec8*NCi2* 
      sqrtxz2i*x2*omzi)/2. -  
   (3*lx*lspec5*NCi2*sqrtxz2i*x2* 
      omzi)/2. + (3*lx*lspec6*NCi2* 
      sqrtxz2i*x2*omzi)/2. -  
   (x*lx*omxmzi2)/2. + (x*lomz*omxmzi2)/2. +  
   (x*lx*NCi2*omxmzi2)/4. -  
   (x*lomz*NCi2*omxmzi2)/4. +  
   (x*lx*NC2*omxmzi2)/4. -  
   (x*lomz*NC2*omxmzi2)/4. +  
   (lx*x2*omxmzi2)/2. -  
   (lomz*x2*omxmzi2)/2. -  
   (lx*NCi2*x2*omxmzi2)/4. +  
   (lomz*NCi2*x2*omxmzi2)/4. -  
   (lx*NC2*x2*omxmzi2)/4. +  
   (lomz*NC2*x2*omxmzi2)/4. -  
   omxmzi/2. + (x*omxmzi)/2. +  
   (3*lx*omxmzi)/2. - (3*lomz*omxmzi)/2. +  
   (NCi2*omxmzi)/4. - (x*NCi2*omxmzi)/4. -  
   (3*lx*NCi2*omxmzi)/4. +  
   (x*lx*NCi2*omxmzi)/4. +  
   (3*lomz*NCi2*omxmzi)/4. -  
   (x*lomz*NCi2*omxmzi)/4. +  
   (NC2*omxmzi)/4. - (x*NC2*omxmzi)/4. -  
   (3*lx*NC2*omxmzi)/4. -  
   (x*lx*NC2*omxmzi)/4. +  
   (3*lomz*NC2*omxmzi)/4. +  
   (x*lomz*NC2*omxmzi)/4. -  
   (lx*NCi2*x2*xmzi2)/4. +  
   (lz*NCi2*x2*xmzi2)/4. +  
   (lx*NC2*x2*xmzi2)/4. -  
   (lz*NC2*x2*xmzi2)/4. +  
   (lx*NCi2*x3*xmzi2)/2. -  
   (lz*NCi2*x3*xmzi2)/2. -  
   (lx*NC2*x3*xmzi2)/2. +  
   (lz*NC2*x3*xmzi2)/2. +  
   lx*x4*xmzi2 - lz*x4*xmzi2 -  
   (3*lx*NCi2*x4*xmzi2)/2. +  
   (3*lz*NCi2*x4*xmzi2)/2. +  
   (lx*NC2*x4*xmzi2)/2. -  
   (lz*NC2*x4*xmzi2)/2. -  
   (3*lx*xmzi)/2. - (3*x*lx*xmzi)/2. +  
   (NC*x*lx*xmzi)/2. + (3*lz*xmzi)/2. +  
   (3*x*lz*xmzi)/2. - (NC*x*lz*xmzi)/2. +  
   (x*NCi2*xmzi)/4. + (3*lx*NCi2*xmzi)/4. +  
   (3*x*lx*NCi2*xmzi)/4. -  
   (3*lz*NCi2*xmzi)/4. -  
   (3*x*lz*NCi2*xmzi)/4. -  
   (x*lx*NCi*xmzi)/2. +  
   (x*lz*NCi*xmzi)/2. - (x*NC2*xmzi)/4. +  
   (3*lx*NC2*xmzi)/4. +  
   (3*x*lx*NC2*xmzi)/4. -  
   (3*lz*NC2*xmzi)/4. -  
   (3*x*lz*NC2*xmzi)/4. +  
   (3*lx*omxi*xmzi)/2. -  
   (3*lz*omxi*xmzi)/2. -  
   (3*lx*NCi2*omxi*xmzi)/4. +  
   (3*lz*NCi2*omxi*xmzi)/4. -  
   (3*lx*NC2*omxi*xmzi)/4. +  
   (3*lz*NC2*omxi*xmzi)/4. -  
   lx*x2*xmzi - NC*lx*x2*xmzi +  
   lz*x2*xmzi + NC*lz*x2*xmzi -  
   (NCi2*x2*xmzi)/2. +  
   lx*NCi*x2*xmzi -  
   lz*NCi*x2*xmzi +  
   (NC2*x2*xmzi)/2. +  
   lx*NC2*x2*xmzi -  
   lz*NC2*x2*xmzi -  
   4*lx*x3*xmzi + NC*lx*x3*xmzi +  
   4*lz*x3*xmzi - NC*lz*x3*xmzi +  
   (NCi2*x3*xmzi)/2. +  
   (9*lx*NCi2*x3*xmzi)/2. -  
   (9*lz*NCi2*x3*xmzi)/2. -  
   lx*NCi*x3*xmzi +  
   lz*NCi*x3*xmzi -  
   (NC2*x3*xmzi)/2. -  
   (lx*NC2*x3*xmzi)/2. +  
   (lz*NC2*x3*xmzi)/2. - 2*zi -  
   (71*NC*zi)/72. - (15*x*zi)/4. + (121*NC*x*zi)/72. -  
   Li2mx*zi - x*Li2mx*zi - 2*Li2x*zi -  
   (NC*Li2x*zi)/2. - 2*x*Li2x*zi -  
   (NC*x*Li2x*zi)/2. - (li2spec1*zi)/2. -  
   (x*li2spec1*zi)/2. +  
   (li2spec2*zi)/2. +  
   (x*li2spec2*zi)/2. + li2spec19*zi -  
   x*li2spec20*zi +  
   (li2spec3*zi)/2. +  
   (x*li2spec3*zi)/2. -  
   (li2spec4*zi)/2. -  
   (x*li2spec4*zi)/2. +  
   li2omxzi*zi + 3*x*li2omxzi*zi +  
   2*sqrtxz1*l2*zi + sqrtxz1*x*l2*zi +  
   (NC*lomx*zi)/12. - (5*NC*x*lomx*zi)/12. +  
   (13*lx*zi)/2. + (5*NC*lx*zi)/6. +  
   sqrtxz1*lx*zi + 15*x*lx*zi +  
   (11*NC*x*lx*zi)/6. + (sqrtxz1*x*lx*zi)/2. +  
   (l2*lx*zi)/2. + (x*l2*lx*zi)/2. -  
   2*lomx*lx*zi - 2*x*lomx*lx*zi -  
   lx*lopx*zi - x*lx*lopx*zi +  
   (NC*lomz*zi)/12. - (5*NC*x*lomz*zi)/12. +  
   (NC*lx*lomz*zi)/2. +  
   (NC*x*lx*lomz*zi)/2. -  
   2*sqrtxz1*lspec1*zi -  
   sqrtxz1*x*lspec1*zi -  
   l2*lspec1*zi -  
   x*l2*lspec1*zi - 5*lz*zi +  
   (NC*lz*zi)/12. + sqrtxz1*lz*zi -  
   (17*x*lz*zi)/2. - (5*NC*x*lz*zi)/12. +  
   (sqrtxz1*x*lz*zi)/2. + (3*l2*lz*zi)/2. +  
   (3*x*l2*lz*zi)/2. + 3*lx*lz*zi +  
   (NC*lx*lz*zi)/2. + 6*x*lx*lz*zi +  
   (NC*x*lx*lz*zi)/2. -  
   (lspec1*lz*zi)/2. -  
   (x*lspec1*lz*zi)/2. -  
   l2*lspec2*zi -  
   x*l2*lspec2*zi -  
   (lx*lspec2*zi)/2. -  
   (x*lx*lspec2*zi)/2. +  
   lspec1*lspec2*zi +  
   x*lspec1*lspec2*zi -  
   lz*lspec2*zi -  
   x*lz*lspec2*zi -  
   (lx*lxpz*zi)/2. - (x*lx*lxpz*zi)/2. +  
   (lz*lxpz*zi)/2. + (x*lz*lxpz*zi)/2. +  
   lx*lopxz*zi + lz*lopxz*zi +  
   (lx*lopxzi*zi)/2. -  
   (x*lx*lopxzi*zi)/2. -  
   (lz*lopxzi*zi)/2. +  
   (x*lz*lopxzi*zi)/2. + 2*NCi2*zi +  
   (15*x*NCi2*zi)/4. + Li2mx*NCi2*zi +  
   x*Li2mx*NCi2*zi + 2*Li2x*NCi2*zi +  
   2*x*Li2x*NCi2*zi +  
   (li2spec1*NCi2*zi)/2. +  
   (x*li2spec1*NCi2*zi)/2. -  
   (li2spec2*NCi2*zi)/2. -  
   (x*li2spec2*NCi2*zi)/2. -  
   li2spec19*NCi2*zi +  
   x*li2spec20*NCi2*zi -  
   (li2spec3*NCi2*zi)/ 
    2. - (x*li2spec3*NCi2* 
      zi)/2. + (li2spec4* 
      NCi2*zi)/2. +  
   (x*li2spec4*NCi2*zi)/ 
    2. - li2omxzi*NCi2*zi -  
   3*x*li2omxzi*NCi2*zi -  
   2*sqrtxz1*l2*NCi2*zi -  
   sqrtxz1*x*l2*NCi2*zi -  
   (13*lx*NCi2*zi)/2. -  
   sqrtxz1*lx*NCi2*zi - 15*x*lx*NCi2*zi -  
   (sqrtxz1*x*lx*NCi2*zi)/2. -  
   (l2*lx*NCi2*zi)/2. -  
   (x*l2*lx*NCi2*zi)/2. +  
   2*lomx*lx*NCi2*zi +  
   2*x*lomx*lx*NCi2*zi +  
   lx*lopx*NCi2*zi +  
   x*lx*lopx*NCi2*zi +  
   2*sqrtxz1*lspec1*NCi2*zi +  
   sqrtxz1*x*lspec1*NCi2*zi +  
   l2*lspec1*NCi2*zi +  
   x*l2*lspec1*NCi2*zi +  
   5*lz*NCi2*zi - sqrtxz1*lz*NCi2*zi +  
   (17*x*lz*NCi2*zi)/2. -  
   (sqrtxz1*x*lz*NCi2*zi)/2. -  
   (3*l2*lz*NCi2*zi)/2. -  
   (3*x*l2*lz*NCi2*zi)/2. -  
   3*lx*lz*NCi2*zi -  
   6*x*lx*lz*NCi2*zi +  
   (lspec1*lz*NCi2*zi)/2. +  
   (x*lspec1*lz*NCi2*zi)/2. +  
   l2*lspec2*NCi2*zi +  
   x*l2*lspec2*NCi2*zi +  
   (lx*lspec2*NCi2*zi)/2. +  
   (x*lx*lspec2*NCi2*zi)/2. -  
   lspec1*lspec2*NCi2*zi -  
   x*lspec1*lspec2*NCi2*zi +  
   lz*lspec2*NCi2*zi +  
   x*lz*lspec2*NCi2*zi +  
   (lx*lxpz*NCi2*zi)/2. +  
   (x*lx*lxpz*NCi2*zi)/2. -  
   (lz*lxpz*NCi2*zi)/2. -  
   (x*lz*lxpz*NCi2*zi)/2. -  
   lx*lopxz*NCi2*zi -  
   lz*lopxz*NCi2*zi -  
   (lx*lopxzi*NCi2*zi)/2. +  
   (x*lx*lopxzi*NCi2*zi)/2. +  
   (lz*lopxzi*NCi2*zi)/2. -  
   (x*lz*lopxzi*NCi2*zi)/2. +  
   (71*NCi*zi)/72. - (121*x*NCi*zi)/72. +  
   (Li2x*NCi*zi)/2. + (x*Li2x*NCi*zi)/2. -  
   (lomx*NCi*zi)/12. +  
   (5*x*lomx*NCi*zi)/12. -  
   (5*lx*NCi*zi)/6. -  
   (11*x*lx*NCi*zi)/6. -  
   (lomz*NCi*zi)/12. +  
   (5*x*lomz*NCi*zi)/12. -  
   (lx*lomz*NCi*zi)/2. -  
   (x*lx*lomz*NCi*zi)/2. -  
   (lz*NCi*zi)/12. +  
   (5*x*lz*NCi*zi)/12. -  
   (lx*lz*NCi*zi)/2. -  
   (x*lx*lz*NCi*zi)/2. + (pi2*zi)/3. +  
   (NC*pi2*zi)/12. + (x*pi2*zi)/6. +  
   (NC*x*pi2*zi)/12. - (NCi2*pi2*zi)/3. -  
   (x*NCi2*pi2*zi)/6. -  
   (NCi*pi2*zi)/12. -  
   (x*NCi*pi2*zi)/12. +  
   (3*li2spec5*sqrtxz2i*zi)/2. +  
   x*li2spec5*sqrtxz2i*zi -  
   (3*li2spec6*sqrtxz2i*zi)/2. -  
   x*li2spec6*sqrtxz2i*zi -  
   (3*li2spec7*sqrtxz2i* 
      zi)/2. - x*li2spec7* 
    sqrtxz2i*zi +  
   (3*li2spec8*sqrtxz2i* 
      zi)/2. + x*li2spec8* 
    sqrtxz2i*zi -  
   (3*lx*lspec5*sqrtxz2i*zi)/2. -  
   x*lx*lspec5*sqrtxz2i*zi +  
   (3*lx*lspec6*sqrtxz2i*zi)/2. +  
   x*lx*lspec6*sqrtxz2i*zi -  
   (3*li2spec5*NCi2*sqrtxz2i*zi)/2. -  
   x*li2spec5*NCi2*sqrtxz2i*zi +  
   (3*li2spec6*NCi2*sqrtxz2i*zi)/2. +  
   x*li2spec6*NCi2*sqrtxz2i*zi +  
   (3*li2spec7*NCi2* 
      sqrtxz2i*zi)/2. +  
   x*li2spec7*NCi2* 
    sqrtxz2i*zi -  
   (3*li2spec8*NCi2* 
      sqrtxz2i*zi)/2. -  
   x*li2spec8*NCi2* 
    sqrtxz2i*zi +  
   (3*lx*lspec5*NCi2*sqrtxz2i*zi)/2. +  
   x*lx*lspec5*NCi2*sqrtxz2i*zi -  
   (3*lx*lspec6*NCi2*sqrtxz2i*zi)/2. -  
   x*lx*lspec6*NCi2*sqrtxz2i*zi +  
   6*omxi*zi + 4*Li2x*omxi*zi +  
   li2spec1*omxi*zi -  
   li2spec2*omxi*zi -  
   li2spec19*omxi*zi +  
   li2spec20*omxi*zi -  
   li2spec3*omxi*zi +  
   li2spec4*omxi*zi -  
   4*li2omxzi*omxi*zi -  
   3*sqrtxz1*l2*omxi*zi -  
   (87*lx*omxi*zi)/4. -  
   (2*NC*lx*omxi*zi)/3. -  
   (3*sqrtxz1*lx*omxi*zi)/2. -  
   l2*lx*omxi*zi +  
   4*lomx*lx*omxi*zi +  
   3*sqrtxz1*lspec1*omxi*zi +  
   2*l2*lspec1*omxi*zi +  
   (27*lz*omxi*zi)/2. -  
   (3*sqrtxz1*lz*omxi*zi)/2. -  
   3*l2*lz*omxi*zi -  
   8*lx*lz*omxi*zi +  
   lspec1*lz*omxi*zi +  
   2*l2*lspec2*omxi*zi +  
   lx*lspec2*omxi*zi -  
   2*lspec1*lspec2*omxi*zi +  
   2*lz*lspec2*omxi*zi +  
   lx*lxpz*omxi*zi -  
   lz*lxpz*omxi*zi -  
   lx*lopxz*omxi*zi -  
   lz*lopxz*omxi*zi -  
   6*NCi2*omxi*zi -  
   4*Li2x*NCi2*omxi*zi -  
   li2spec1*NCi2*omxi*zi +  
   li2spec2*NCi2*omxi*zi +  
   li2spec19*NCi2*omxi*zi -  
   li2spec20*NCi2*omxi*zi +  
   li2spec3*NCi2*omxi* 
    zi - li2spec4*NCi2* 
    omxi*zi + 4*li2omxzi*NCi2*omxi* 
    zi + 3*sqrtxz1*l2*NCi2*omxi*zi +  
   (87*lx*NCi2*omxi*zi)/4. +  
   (3*sqrtxz1*lx*NCi2*omxi*zi)/2. +  
   l2*lx*NCi2*omxi*zi -  
   4*lomx*lx*NCi2*omxi*zi -  
   3*sqrtxz1*lspec1*NCi2*omxi*zi -  
   2*l2*lspec1*NCi2*omxi*zi -  
   (27*lz*NCi2*omxi*zi)/2. +  
   (3*sqrtxz1*lz*NCi2*omxi*zi)/2. +  
   3*l2*lz*NCi2*omxi*zi +  
   8*lx*lz*NCi2*omxi*zi -  
   lspec1*lz*NCi2*omxi*zi -  
   2*l2*lspec2*NCi2*omxi*zi -  
   lx*lspec2*NCi2*omxi*zi +  
   2*lspec1*lspec2*NCi2*omxi* 
    zi - 2*lz*lspec2*NCi2*omxi* 
    zi - lx*lxpz*NCi2*omxi*zi +  
   lz*lxpz*NCi2*omxi*zi +  
   lx*lopxz*NCi2*omxi*zi +  
   lz*lopxz*NCi2*omxi*zi +  
   (2*lx*NCi*omxi*zi)/3. -  
   (2*pi2*omxi*zi)/3. +  
   (2*NCi2*pi2*omxi*zi)/3. -  
   Li2mx*xi2*zi + Li2x*xi2*zi +  
   (li2spec19*xi2*zi)/2. +  
   (li2spec20*xi2*zi)/2. -  
   2*lx*xi2*zi + lomx*lx*xi2*zi -  
   lx*lopx*xi2*zi -  
   (lx*lz*xi2*zi)/2. +  
   (lx*lopxz*xi2*zi)/2. +  
   (lz*lopxz*xi2*zi)/2. +  
   (lx*lopxzi*xi2*zi)/2. -  
   (lz*lopxzi*xi2*zi)/2. +  
   Li2mx*NCi2*xi2*zi -  
   Li2x*NCi2*xi2*zi -  
   (li2spec19*NCi2*xi2*zi)/2. -  
   (li2spec20*NCi2*xi2*zi)/2. +  
   2*lx*NCi2*xi2*zi -  
   lomx*lx*NCi2*xi2*zi +  
   lx*lopx*NCi2*xi2*zi +  
   (lx*lz*NCi2*xi2*zi)/2. -  
   (lx*lopxz*NCi2*xi2*zi)/2. -  
   (lz*lopxz*NCi2*xi2*zi)/2. -  
   (lx*lopxzi*NCi2*xi2*zi)/2. +  
   (lz*lopxzi*NCi2*xi2*zi)/2. -  
   (pi2*xi2*zi)/6. +  
   (NCi2*pi2*xi2*zi)/6. +  
   (13*NC*xi*zi)/18. + Li2mx*xi*zi -  
   Li2x*xi*zi - (li2spec19*xi*zi)/2. -  
   (li2spec20*xi*zi)/2. +  
   (NC*lomx*xi*zi)/3. + 2*lx*xi*zi -  
   lomx*lx*xi*zi +  
   lx*lopx*xi*zi +  
   (NC*lomz*xi*zi)/3. +  
   (NC*lz*xi*zi)/3. +  
   (lx*lz*xi*zi)/2. -  
   (lx*lopxz*xi*zi)/2. -  
   (lz*lopxz*xi*zi)/2. -  
   (lx*lopxzi*xi*zi)/2. +  
   (lz*lopxzi*xi*zi)/2. -  
   Li2mx*NCi2*xi*zi +  
   Li2x*NCi2*xi*zi +  
   (li2spec19*NCi2*xi*zi)/2. +  
   (li2spec20*NCi2*xi*zi)/2. -  
   2*lx*NCi2*xi*zi +  
   lomx*lx*NCi2*xi*zi -  
   lx*lopx*NCi2*xi*zi -  
   (lx*lz*NCi2*xi*zi)/2. +  
   (lx*lopxz*NCi2*xi*zi)/2. +  
   (lz*lopxz*NCi2*xi*zi)/2. +  
   (lx*lopxzi*NCi2*xi*zi)/2. -  
   (lz*lopxzi*NCi2*xi*zi)/2. -  
   (13*NCi*xi*zi)/18. -  
   (lomx*NCi*xi*zi)/3. -  
   (lomz*NCi*xi*zi)/3. -  
   (lz*NCi*xi*zi)/3. +  
   (pi2*xi*zi)/6. -  
   (NCi2*pi2*xi*zi)/6. -  
   (11*NC*x2*zi)/9. - (NC*lomx*x2*zi)/3. +  
   NC*lx*x2*zi - (NC*lomz*x2*zi)/3. -  
   (NC*lz*x2*zi)/3. +  
   (11*NCi*x2*zi)/9. +  
   (lomx*NCi*x2*zi)/3. -  
   lx*NCi*x2*zi +  
   (lomz*NCi*x2*zi)/3. +  
   (lz*NCi*x2*zi)/3. +  
   (3*li2spec5*sqrtxz2i*x2*zi)/2. -  
   (3*li2spec6*sqrtxz2i*x2*zi)/2. -  
   (3*li2spec7*sqrtxz2i* 
      x2*zi)/2. + (3* 
      li2spec8*sqrtxz2i* 
      x2*zi)/2. - (3*lx*lspec5*sqrtxz2i* 
      x2*zi)/2. + (3*lx*lspec6*sqrtxz2i* 
      x2*zi)/2. - (3*li2spec5*NCi2* 
      sqrtxz2i*x2*zi)/2. +  
   (3*li2spec6*NCi2*sqrtxz2i*x2* 
      zi)/2. + (3*li2spec7* 
      NCi2*sqrtxz2i*x2*zi)/2. -  
   (3*li2spec8*NCi2* 
      sqrtxz2i*x2*zi)/2. +  
   (3*lx*lspec5*NCi2*sqrtxz2i*x2* 
      zi)/2. - (3*lx*lspec6*NCi2* 
      sqrtxz2i*x2*zi)/2. -  
   Li2mx*opxi*zi + Li2x*opxi*zi -  
   (li2spec19*opxi*zi)/2. -  
   (li2spec20*opxi*zi)/2. -  
   2*lx*opxi*zi +  
   lomx*lx*opxi*zi -  
   lx*lopx*opxi*zi +  
   (lx*lz*opxi*zi)/2. -  
   (lx*lopxz*opxi*zi)/2. -  
   (lz*lopxz*opxi*zi)/2. -  
   (lx*lopxzi*opxi*zi)/2. +  
   (lz*lopxzi*opxi*zi)/2. +  
   Li2mx*NCi2*opxi*zi -  
   Li2x*NCi2*opxi*zi +  
   (li2spec19*NCi2*opxi*zi)/2. +  
   (li2spec20*NCi2*opxi*zi)/2. +  
   2*lx*NCi2*opxi*zi -  
   lomx*lx*NCi2*opxi*zi +  
   lx*lopx*NCi2*opxi*zi -  
   (lx*lz*NCi2*opxi*zi)/2. +  
   (lx*lopxz*NCi2*opxi*zi)/2. +  
   (lz*lopxz*NCi2*opxi*zi)/2. +  
   (lx*lopxzi*NCi2*opxi*zi)/2. -  
   (lz*lopxzi*NCi2*opxi*zi)/2. -  
   (pi2*opxi*zi)/3. +  
   (NCi2*pi2*opxi*zi)/3. +  
   Li2mx*xi2*opxi*zi -  
   Li2x*xi2*opxi*zi -  
   (li2spec19*xi2*opxi*zi)/2. -  
   (li2spec20*xi2*opxi*zi)/2. +  
   2*lx*xi2*opxi*zi -  
   lomx*lx*xi2*opxi*zi +  
   lx*lopx*xi2*opxi*zi +  
   (lx*lz*xi2*opxi*zi)/2. -  
   (lx*lopxz*xi2*opxi*zi)/2. -  
   (lz*lopxz*xi2*opxi*zi)/2. -  
   (lx*lopxzi*xi2*opxi*zi)/2. +  
   (lz*lopxzi*xi2*opxi*zi)/2. -  
   Li2mx*NCi2*xi2*opxi*zi +  
   Li2x*NCi2*xi2*opxi*zi +  
   (li2spec19*NCi2*xi2*opxi*zi)/2. +  
   (li2spec20*NCi2*xi2*opxi*zi)/2. -  
   2*lx*NCi2*xi2*opxi*zi +  
   lomx*lx*NCi2*xi2*opxi*zi -  
   lx*lopx*NCi2*xi2*opxi*zi -  
   (lx*lz*NCi2*xi2*opxi*zi)/2. +  
   (lx*lopxz*NCi2*xi2*opxi*zi)/2. +  
   (lz*lopxz*NCi2*xi2*opxi*zi)/2. +  
   (lx*lopxzi*NCi2*xi2*opxi*zi)/ 
    2. - (lz*lopxzi*NCi2*xi2*opxi* 
      zi)/2. + (pi2*xi2*opxi*zi)/6. -  
   (NCi2*pi2*xi2*opxi*zi)/6. +  
   2*omzi*zi + 4*x*omzi*zi +  
   2*Li2x*omzi*zi + 2*x*Li2x*omzi*zi -  
   li2omxzi*omzi*zi -  
   3*x*li2omxzi*omzi*zi -  
   8*lx*omzi*zi - 16*x*lx*omzi*zi +  
   2*lomx*lx*omzi*zi +  
   2*x*lomx*lx*omzi*zi +  
   4*lz*omzi*zi + 8*x*lz*omzi*zi -  
   3*lx*lz*omzi*zi -  
   5*x*lx*lz*omzi*zi -  
   2*NCi2*omzi*zi -  
   4*x*NCi2*omzi*zi -  
   2*Li2x*NCi2*omzi*zi -  
   2*x*Li2x*NCi2*omzi*zi +  
   li2omxzi*NCi2*omzi*zi +  
   3*x*li2omxzi*NCi2*omzi*zi +  
   8*lx*NCi2*omzi*zi +  
   16*x*lx*NCi2*omzi*zi -  
   2*lomx*lx*NCi2*omzi*zi -  
   2*x*lomx*lx*NCi2*omzi*zi -  
   4*lz*NCi2*omzi*zi -  
   8*x*lz*NCi2*omzi*zi +  
   3*lx*lz*NCi2*omzi*zi +  
   5*x*lx*lz*NCi2*omzi*zi -  
   (pi2*omzi*zi)/3. -  
   (x*pi2*omzi*zi)/3. +  
   (NCi2*pi2*omzi*zi)/3. +  
   (x*NCi2*pi2*omzi*zi)/3. -  
   6*omxi*omzi*zi -  
   4*Li2x*omxi*omzi*zi +  
   4*li2omxzi*omxi*omzi*zi +  
   24*lx*omxi*omzi*zi -  
   4*lomx*lx*omxi*omzi*zi -  
   12*lz*omxi*omzi*zi +  
   8*lx*lz*omxi*omzi*zi +  
   6*NCi2*omxi*omzi*zi +  
   4*Li2x*NCi2*omxi*omzi*zi -  
   4*li2omxzi*NCi2*omxi*omzi*zi -  
   24*lx*NCi2*omxi*omzi*zi +  
   4*lomx*lx*NCi2*omxi*omzi*zi +  
   12*lz*NCi2*omxi*omzi*zi -  
   8*lx*lz*NCi2*omxi*omzi*zi +  
   (2*pi2*omxi*omzi*zi)/3. -  
   (2*NCi2*pi2*omxi*omzi*zi)/3. +  
   (43*NC*z2)/36. - (77*NC*x*z2)/36. -  
   (5*sqrtxz3*itani1*z2)/4. +  
   (23*NC*sqrtxz3*itani1*z2)/4. +  
   (5*sqrtxz3*itani2*z2)/4. -  
   (23*NC*sqrtxz3*itani2*z2)/4. -  
   (5*sqrtxz3*itani3*z2)/2. +  
   (23*NC*sqrtxz3*itani3*z2)/2. -  
   4*NC*x*Li2x*z2 + 4*NC*x*li2spec1*z2 -  
   4*NC*x*li2spec2*z2 -  
   4*NC*x*li2spec19*z2 -  
   (5*sqrtxz3*atanspec1*lspec3*z2)/2. +  
   (23*NC*sqrtxz3*atanspec1*lspec3*z2)/2. -  
   (NC*lomx*z2)/3. + (2*NC*x*lomx*z2)/3. +  
   (2*NC*lx*z2)/3. - (4*NC*x*lx*z2)/3. -  
   8*NC*x*l2*lx*z2 - 4*NC*x*lomx*lx*z2 -  
   (NC*lomz*z2)/3. + (2*NC*x*lomz*z2)/3. +  
   12*NC*x*l2*lspec1*z2 +  
   4*NC*x*lx*lspec1*z2 - (NC*lz*z2)/3. +  
   (2*NC*x*lz*z2)/3. - 8*NC*x*l2*lz*z2 -  
   2*NC*x*lx*lz*z2 +  
   4*NC*x*lspec1*lz*z2 +  
   (5*sqrtxz3*atanspec2*lspec4*z2)/2. -  
   (23*NC*sqrtxz3*atanspec2*lspec4*z2)/2. +  
   4*NC*x*l2*lspec2*z2 +  
   4*NC*x*lx*lspec2*z2 -  
   4*NC*x*lspec1*lspec2*z2 +  
   4*NC*x*lz*lspec2*z2 +  
   2*NC*x*lx*lxpz*z2 - 2*NC*x*lz*lxpz*z2 -  
   4*NC*x*lx*lopxz*z2 -  
   4*NC*x*lz*lopxz*z2 -  
   2*NC*x*lx*lopxzi*z2 +  
   2*NC*x*lz*lopxzi*z2 +  
   (5*sqrtxz3*itani1*NCi2*z2)/4. -  
   (5*sqrtxz3*itani2*NCi2*z2)/4. +  
   (5*sqrtxz3*itani3*NCi2*z2)/2. +  
   (5*sqrtxz3*atanspec1*lspec3*NCi2*z2)/2. -  
   (5*sqrtxz3*atanspec2*lspec4*NCi2*z2)/2. -  
   (43*NCi*z2)/36. + (77*x*NCi*z2)/36. -  
   (23*sqrtxz3*itani1*NCi*z2)/4. +  
   (23*sqrtxz3*itani2*NCi*z2)/4. -  
   (23*sqrtxz3*itani3*NCi*z2)/2. +  
   4*x*Li2x*NCi*z2 -  
   4*x*li2spec1*NCi*z2 +  
   4*x*li2spec2*NCi*z2 +  
   4*x*li2spec19*NCi*z2 -  
   (23*sqrtxz3*atanspec1*lspec3*NCi*z2)/2. +  
   (lomx*NCi*z2)/3. -  
   (2*x*lomx*NCi*z2)/3. -  
   (2*lx*NCi*z2)/3. + (4*x*lx*NCi*z2)/3. +  
   8*x*l2*lx*NCi*z2 +  
   4*x*lomx*lx*NCi*z2 +  
   (lomz*NCi*z2)/3. -  
   (2*x*lomz*NCi*z2)/3. -  
   12*x*l2*lspec1*NCi*z2 -  
   4*x*lx*lspec1*NCi*z2 +  
   (lz*NCi*z2)/3. - (2*x*lz*NCi*z2)/3. +  
   8*x*l2*lz*NCi*z2 +  
   2*x*lx*lz*NCi*z2 -  
   4*x*lspec1*lz*NCi*z2 +  
   (23*sqrtxz3*atanspec2*lspec4*NCi*z2)/2. -  
   4*x*l2*lspec2*NCi*z2 -  
   4*x*lx*lspec2*NCi*z2 +  
   4*x*lspec1*lspec2*NCi*z2 -  
   4*x*lz*lspec2*NCi*z2 -  
   2*x*lx*lxpz*NCi*z2 +  
   2*x*lz*lxpz*NCi*z2 +  
   4*x*lx*lopxz*NCi*z2 +  
   4*x*lz*lopxz*NCi*z2 +  
   2*x*lx*lopxzi*NCi*z2 -  
   2*x*lz*lopxzi*NCi*z2 +  
   (2*NC*x*pi2*z2)/3. - (2*x*NCi*pi2*z2)/3. -  
   9*x*li2spec5*sqrtxz2i*z2 +  
   8*NC*x*li2spec5*sqrtxz2i*z2 +  
   9*x*li2spec6*sqrtxz2i*z2 -  
   8*NC*x*li2spec6*sqrtxz2i*z2 +  
   9*x*li2spec7*sqrtxz2i* 
    z2 - 8*NC*x*li2spec7* 
    sqrtxz2i*z2 - 9*x* 
    li2spec8*sqrtxz2i*z2 
    + 8*NC*x*li2spec8*sqrtxz2i* 
    z2 + 9*x*lx*lspec5*sqrtxz2i*z2 -  
   8*NC*x*lx*lspec5*sqrtxz2i*z2 -  
   9*x*lx*lspec6*sqrtxz2i*z2 +  
   8*NC*x*lx*lspec6*sqrtxz2i*z2 +  
   9*x*li2spec5*NCi2*sqrtxz2i*z2 -  
   9*x*li2spec6*NCi2*sqrtxz2i*z2 -  
   9*x*li2spec7*NCi2* 
    sqrtxz2i*z2 + 9*x* 
    li2spec8*NCi2* 
    sqrtxz2i*z2 - 9*x*lx*lspec5*NCi2* 
    sqrtxz2i*z2 + 9*x*lx*lspec6*NCi2* 
    sqrtxz2i*z2 - 8*x*li2spec5*NCi* 
    sqrtxz2i*z2 + 8*x*li2spec6*NCi* 
    sqrtxz2i*z2 + 8*x* 
    li2spec7*NCi* 
    sqrtxz2i*z2 - 8*x* 
    li2spec8*NCi* 
    sqrtxz2i*z2 + 8*x*lx*lspec5*NCi* 
    sqrtxz2i*z2 - 8*x*lx*lspec6*NCi* 
    sqrtxz2i*z2 + 4*NC*Li2x*omxi*z2 -  
   4*NC*li2spec1*omxi*z2 +  
   4*NC*li2spec2*omxi*z2 +  
   2*NC*li2spec19*omxi*z2 -  
   2*NC*li2spec20*omxi*z2 +  
   (2*NC*lx*omxi*z2)/3. +  
   8*NC*l2*lx*omxi*z2 +  
   4*NC*lomx*lx*omxi*z2 -  
   12*NC*l2*lspec1*omxi*z2 -  
   4*NC*lx*lspec1*omxi*z2 +  
   8*NC*l2*lz*omxi*z2 +  
   4*NC*lx*lz*omxi*z2 -  
   4*NC*lspec1*lz*omxi*z2 -  
   4*NC*l2*lspec2*omxi*z2 -  
   4*NC*lx*lspec2*omxi*z2 +  
   4*NC*lspec1*lspec2*omxi*z2 -  
   4*NC*lz*lspec2*omxi*z2 -  
   2*NC*lx*lxpz*omxi*z2 +  
   2*NC*lz*lxpz*omxi*z2 +  
   2*NC*lx*lopxz*omxi*z2 +  
   2*NC*lz*lopxz*omxi*z2 -  
   4*Li2x*NCi*omxi*z2 +  
   4*li2spec1*NCi*omxi*z2 -  
   4*li2spec2*NCi*omxi*z2 -  
   2*li2spec19*NCi*omxi*z2 +  
   2*li2spec20*NCi*omxi*z2 -  
   (2*lx*NCi*omxi*z2)/3. -  
   8*l2*lx*NCi*omxi*z2 -  
   4*lomx*lx*NCi*omxi*z2 +  
   12*l2*lspec1*NCi*omxi*z2 +  
   4*lx*lspec1*NCi*omxi*z2 -  
   8*l2*lz*NCi*omxi*z2 -  
   4*lx*lz*NCi*omxi*z2 +  
   4*lspec1*lz*NCi*omxi*z2 +  
   4*l2*lspec2*NCi*omxi*z2 +  
   4*lx*lspec2*NCi*omxi*z2 -  
   4*lspec1*lspec2*NCi*omxi* 
    z2 + 4*lz*lspec2*NCi*omxi* 
    z2 + 2*lx*lxpz*NCi*omxi*z2 -  
   2*lz*lxpz*NCi*omxi*z2 -  
   2*lx*lopxz*NCi*omxi*z2 -  
   2*lz*lopxz*NCi*omxi*z2 -  
   NC*pi2*omxi*z2 +  
   NCi*pi2*omxi*z2 -  
   4*NC*Li2mx*xi2*z2 + 4*NC*Li2x*xi2*z2 +  
   2*NC*li2spec19*xi2*z2 +  
   2*NC*li2spec20*xi2*z2 -  
   8*NC*lx*xi2*z2 +  
   4*NC*lomx*lx*xi2*z2 -  
   4*NC*lx*lopx*xi2*z2 -  
   2*NC*lx*lz*xi2*z2 +  
   2*NC*lx*lopxz*xi2*z2 +  
   2*NC*lz*lopxz*xi2*z2 +  
   2*NC*lx*lopxzi*xi2*z2 -  
   2*NC*lz*lopxzi*xi2*z2 +  
   4*Li2mx*NCi*xi2*z2 -  
   4*Li2x*NCi*xi2*z2 -  
   2*li2spec19*NCi*xi2*z2 -  
   2*li2spec20*NCi*xi2*z2 +  
   8*lx*NCi*xi2*z2 -  
   4*lomx*lx*NCi*xi2*z2 +  
   4*lx*lopx*NCi*xi2*z2 +  
   2*lx*lz*NCi*xi2*z2 -  
   2*lx*lopxz*NCi*xi2*z2 -  
   2*lz*lopxz*NCi*xi2*z2 -  
   2*lx*lopxzi*NCi*xi2*z2 +  
   2*lz*lopxzi*NCi*xi2*z2 -  
   (2*NC*pi2*xi2*z2)/3. +  
   (2*NCi*pi2*xi2*z2)/3. -  
   2*NC*li2spec19*opxi*z2 -  
   2*NC*li2spec20*opxi*z2 +  
   2*NC*lx*lz*opxi*z2 -  
   2*NC*lx*lopxz*opxi*z2 -  
   2*NC*lz*lopxz*opxi*z2 -  
   2*NC*lx*lopxzi*opxi*z2 +  
   2*NC*lz*lopxzi*opxi*z2 +  
   2*li2spec19*NCi*opxi*z2 +  
   2*li2spec20*NCi*opxi*z2 -  
   2*lx*lz*NCi*opxi*z2 +  
   2*lx*lopxz*NCi*opxi*z2 +  
   2*lz*lopxz*NCi*opxi*z2 +  
   2*lx*lopxzi*NCi*opxi*z2 -  
   2*lz*lopxzi*NCi*opxi*z2 -  
   (NC*pi2*opxi*z2)/3. +  
   (NCi*pi2*opxi*z2)/3. +  
   4*NC*Li2mx*xi2*opxi*z2 -  
   4*NC*Li2x*xi2*opxi*z2 -  
   2*NC*li2spec19*xi2*opxi*z2 -  
   2*NC*li2spec20*xi2*opxi*z2 +  
   8*NC*lx*xi2*opxi*z2 -  
   4*NC*lomx*lx*xi2*opxi*z2 +  
   4*NC*lx*lopx*xi2*opxi*z2 +  
   2*NC*lx*lz*xi2*opxi*z2 -  
   2*NC*lx*lopxz*xi2*opxi*z2 -  
   2*NC*lz*lopxz*xi2*opxi*z2 -  
   2*NC*lx*lopxzi*xi2*opxi*z2 +  
   2*NC*lz*lopxzi*xi2*opxi*z2 -  
   4*Li2mx*NCi*xi2*opxi*z2 +  
   4*Li2x*NCi*xi2*opxi*z2 +  
   2*li2spec19*NCi*xi2*opxi*z2 +  
   2*li2spec20*NCi*xi2*opxi*z2 -  
   8*lx*NCi*xi2*opxi*z2 +  
   4*lomx*lx*NCi*xi2*opxi*z2 -  
   4*lx*lopx*NCi*xi2*opxi*z2 -  
   2*lx*lz*NCi*xi2*opxi*z2 +  
   2*lx*lopxz*NCi*xi2*opxi*z2 +  
   2*lz*lopxz*NCi*xi2*opxi*z2 +  
   2*lx*lopxzi*NCi*xi2*opxi* 
    z2 - 2*lz*lopxzi*NCi*xi2* 
    opxi*z2 + (2*NC*pi2*xi2*opxi* 
      z2)/3. - (2*NCi*pi2*xi2*opxi* 
      z2)/3. + 4*NC*Li2mx*xi*opxi*z2 -  
   4*NC*Li2x*xi*opxi*z2 -  
   2*NC*li2spec19*xi*opxi*z2 -  
   2*NC*li2spec20*xi*opxi*z2 +  
   8*NC*lx*xi*opxi*z2 -  
   4*NC*lomx*lx*xi*opxi*z2 +  
   4*NC*lx*lopx*xi*opxi*z2 +  
   2*NC*lx*lz*xi*opxi*z2 -  
   2*NC*lx*lopxz*xi*opxi*z2 -  
   2*NC*lz*lopxz*xi*opxi*z2 -  
   2*NC*lx*lopxzi*xi*opxi*z2 +  
   2*NC*lz*lopxzi*xi*opxi*z2 -  
   4*Li2mx*NCi*xi*opxi*z2 +  
   4*Li2x*NCi*xi*opxi*z2 +  
   2*li2spec19*NCi*xi*opxi*z2 +  
   2*li2spec20*NCi*xi*opxi*z2 -  
   8*lx*NCi*xi*opxi*z2 +  
   4*lomx*lx*NCi*xi*opxi*z2 -  
   4*lx*lopx*NCi*xi*opxi*z2 -  
   2*lx*lz*NCi*xi*opxi*z2 +  
   2*lx*lopxz*NCi*xi*opxi*z2 +  
   2*lz*lopxz*NCi*xi*opxi*z2 +  
   2*lx*lopxzi*NCi*xi*opxi* 
    z2 - 2*lz*lopxzi*NCi*xi* 
    opxi*z2 + (2*NC*pi2*xi*opxi* 
      z2)/3. - (2*NCi*pi2*xi*opxi* 
      z2)/3. - x*lx*xmzi2*z3 +  
   x*lz*xmzi2*z3 +  
   x*lx*NCi2*xmzi2*z3 -  
   x*lz*NCi2*xmzi2*z3 + 2*l22 +  
   2*x*l22 - 2*NC*x*l22 - 4*NC*z*l22 +  
   2*x*z*l22 - 2*NCi2*l22 -  
   2*x*NCi2*l22 - 2*x*z*NCi2*l22 +  
   2*x*NCi*l22 + 4*z*NCi*l22 -  
   4*omxi*l22 + 3*NC*omxi*l22 -  
   2*z*omxi*l22 + 2*NC*z*omxi*l22 +  
   4*NCi2*omxi*l22 +  
   2*z*NCi2*omxi*l22 -  
   3*NCi*omxi*l22 -  
   2*z*NCi*omxi*l22 + zi*l22 +  
   x*zi*l22 - NCi2*zi*l22 -  
   x*NCi2*zi*l22 -  
   2*omxi*zi*l22 +  
   2*NCi2*omxi*zi*l22 -  
   8*NC*x*z2*l22 + 8*x*NCi*z2*l22 +  
   8*NC*omxi*z2*l22 -  
   8*NCi*omxi*z2*l22 - 4*lomx2 -  
   (5*x*lomx2)/2. - (z*lomx2)/2. -  
   6*x*z*lomx2 + (5*NCi2*lomx2)/4. +  
   (9*x*NCi2*lomx2)/4. +  
   (z*NCi2*lomx2)/4. +  
   (13*x*z*NCi2*lomx2)/4. +  
   (11*NC2*lomx2)/4. + (x*NC2*lomx2)/4. +  
   (z*NC2*lomx2)/4. +  
   (11*x*z*NC2*lomx2)/4. +  
   (5*omxi*lomx2)/2. +  
   (5*z*omxi*lomx2)/2. -  
   (3*NCi2*omxi*lomx2)/2. -  
   (3*z*NCi2*omxi*lomx2)/2. -  
   NC2*omxi*lomx2 -  
   z*NC2*omxi*lomx2 +  
   (5*omzi*lomx2)/4. +  
   (13*x*omzi*lomx2)/4. -  
   (NCi2*omzi*lomx2)/4. -  
   (9*x*NCi2*omzi*lomx2)/4. -  
   NC2*omzi*lomx2 -  
   x*NC2*omzi*lomx2 -  
   (9*omxi*omzi*lomx2)/2. +  
   (5*NCi2*omxi*omzi*lomx2)/2. +  
   2*NC2*omxi*omzi*lomx2 -  
   (17*lx2)/2. + 2*NC*lx2 - (11*x*lx2)/2. +  
   (5*NC*x*lx2)/2. - (z*lx2)/2. - 2*NC*z*lx2 -  
   12*x*z*lx2 - 3*NC*x*z*lx2 +  
   (17*NCi2*lx2)/4. + (21*x*NCi2*lx2)/4. +  
   (z*NCi2*lx2)/4. + (31*x*z*NCi2*lx2)/4. -  
   2*NCi*lx2 - (5*x*NCi*lx2)/2. +  
   2*z*NCi*lx2 + 3*x*z*NCi*lx2 +  
   (17*NC2*lx2)/4. + (x*NC2*lx2)/4. +  
   (z*NC2*lx2)/4. + (17*x*z*NC2*lx2)/4. +  
   9*omxi*lx2 - NC*omxi*lx2 +  
   9*z*omxi*lx2 + 2*NC*z*omxi*lx2 -  
   (13*NCi2*omxi*lx2)/2. -  
   (13*z*NCi2*omxi*lx2)/2. +  
   NCi*omxi*lx2 -  
   2*z*NCi*omxi*lx2 -  
   (5*NC2*omxi*lx2)/2. -  
   (5*z*NC2*omxi*lx2)/2. -  
   (3*xi2*lx2)/2. + (3*NC*xi2*lx2)/4. +  
   (3*z*xi2*lx2)/2. - (3*NC*z*xi2*lx2)/2. +  
   (3*NCi2*xi2*lx2)/2. -  
   (3*z*NCi2*xi2*lx2)/2. -  
   (3*NCi*xi2*lx2)/4. +  
   (3*z*NCi*xi2*lx2)/2. +  
   (3*xi*lx2)/2. - (3*NC*xi*lx2)/4. +  
   (3*NC*z*xi*lx2)/2. -  
   (3*NCi2*xi*lx2)/2. +  
   (3*NCi*xi*lx2)/4. -  
   (3*z*NCi*xi*lx2)/2. -  
   (5*opxi*lx2)/2. +  
   (5*NC*opxi*lx2)/4. +  
   (z*opxi*lx2)/2. -  
   (5*NC*z*opxi*lx2)/2. +  
   (5*NCi2*opxi*lx2)/2. -  
   (z*NCi2*opxi*lx2)/2. -  
   (5*NCi*opxi*lx2)/4. +  
   (5*z*NCi*opxi*lx2)/2. +  
   (3*xi2*opxi*lx2)/2. -  
   (3*NC*xi2*opxi*lx2)/4. -  
   (3*z*xi2*opxi*lx2)/2. +  
   (3*NC*z*xi2*opxi*lx2)/2. -  
   (3*NCi2*xi2*opxi*lx2)/2. +  
   (3*z*NCi2*xi2*opxi*lx2)/2. +  
   (3*NCi*xi2*opxi*lx2)/4. -  
   (3*z*NCi*xi2*opxi*lx2)/2. -  
   (3*z*xi*opxi*lx2)/2. +  
   (3*z*NCi2*xi*opxi*lx2)/2. +  
   (7*omzi*lx2)/4. +  
   (11*x*omzi*lx2)/4. +  
   (NCi2*omzi*lx2)/4. -  
   (3*x*NCi2*omzi*lx2)/4. -  
   2*NC2*omzi*lx2 -  
   2*x*NC2*omzi*lx2 -  
   (9*omxi*omzi*lx2)/2. +  
   (NCi2*omxi*omzi*lx2)/2. +  
   4*NC2*omxi*omzi*lx2 -  
   (5*zi*lx2)/2. - (NC*zi*lx2)/2. -  
   5*x*zi*lx2 - (NC*x*zi*lx2)/2. +  
   (5*NCi2*zi*lx2)/2. +  
   5*x*NCi2*zi*lx2 +  
   (NCi*zi*lx2)/2. +  
   (x*NCi*zi*lx2)/2. +  
   8*omxi*zi*lx2 -  
   8*NCi2*omxi*zi*lx2 +  
   (3*xi2*zi*lx2)/4. -  
   (3*NCi2*xi2*zi*lx2)/4. -  
   (3*xi*zi*lx2)/4. +  
   (3*NCi2*xi*zi*lx2)/4. +  
   (5*opxi*zi*lx2)/4. -  
   (5*NCi2*opxi*zi*lx2)/4. -  
   (3*xi2*opxi*zi*lx2)/4. +  
   (3*NCi2*xi2*opxi*zi*lx2)/4. +  
   (5*omzi*zi*lx2)/2. +  
   (11*x*omzi*zi*lx2)/2. -  
   (5*NCi2*omzi*zi*lx2)/2. -  
   (11*x*NCi2*omzi*zi*lx2)/2. -  
   8*omxi*omzi*zi*lx2 +  
   8*NCi2*omxi*omzi*zi*lx2 +  
   2*NC*x*z2*lx2 - 2*x*NCi*z2*lx2 -  
   NC*omxi*z2*lx2 +  
   NCi*omxi*z2*lx2 +  
   3*NC*xi2*z2*lx2 -  
   3*NCi*xi2*z2*lx2 +  
   NC*opxi*z2*lx2 -  
   NCi*opxi*z2*lx2 -  
   3*NC*xi2*opxi*z2*lx2 +  
   3*NCi*xi2*opxi*z2*lx2 -  
   3*NC*xi*opxi*z2*lx2 +  
   3*NCi*xi*opxi*z2*lx2 -  
   (11*lomz2)/2. - x*lomz2 -  
   (z*lomz2)/2. - 6*x*z*lomz2 +  
   (11*NCi2*lomz2)/4. +  
   (3*x*NCi2*lomz2)/4. +  
   (z*NCi2*lomz2)/4. +  
   (13*x*z*NCi2*lomz2)/4. +  
   (11*NC2*lomz2)/4. + (x*NC2*lomz2)/4. +  
   (z*NC2*lomz2)/4. +  
   (11*x*z*NC2*lomz2)/4. +  
   (13*omxi*lomz2)/4. +  
   (13*z*omxi*lomz2)/4. -  
   (9*NCi2*omxi*lomz2)/4. -  
   (9*z*NCi2*omxi*lomz2)/4. -  
   NC2*omxi*lomz2 -  
   z*NC2*omxi*lomz2 +  
   2*omzi*lomz2 +  
   (5*x*omzi*lomz2)/2. -  
   NCi2*omzi*lomz2 -  
   (3*x*NCi2*omzi*lomz2)/2. -  
   NC2*omzi*lomz2 -  
   x*NC2*omzi*lomz2 -  
   (9*omxi*omzi*lomz2)/2. +  
   (5*NCi2*omxi*omzi*lomz2)/2. +  
   2*NC2*omxi*omzi*lomz2 +  
   NC*lspec1_2 - NC*x*lspec1_2 -  
   2*NC*z*lspec1_2 +  
   2*NC*x*z*lspec1_2 -  
   NCi*lspec1_2 +  
   x*NCi*lspec1_2 +  
   2*z*NCi*lspec1_2 -  
   2*x*z*NCi*lspec1_2 +  
   NC*omxi*lspec1_2 -  
   2*NC*z*omxi*lspec1_2 -  
   NCi*omxi*lspec1_2 +  
   2*z*NCi*omxi*lspec1_2 -  
   4*NC*x*z2*lspec1_2 +  
   4*x*NCi*z2*lspec1_2 +  
   4*NC*omxi*z2*lspec1_2 -  
   4*NCi*omxi*z2*lspec1_2 +  
   lz2/2. - NC*lz2 - (x*lz2)/2. +  
   (NC*x*lz2)/2. + (z*lz2)/2. + NC*z*lz2 -  
   3*NC*x*z*lz2 - (NCi2*lz2)/4. +  
   (3*x*NCi2*lz2)/4. - (z*NCi2*lz2)/4. +  
   (x*z*NCi2*lz2)/4. + NCi*lz2 -  
   (x*NCi*lz2)/2. - z*NCi*lz2 +  
   3*x*z*NCi*lz2 - (NC2*lz2)/4. -  
   (x*NC2*lz2)/4. - (z*NC2*lz2)/4. -  
   (x*z*NC2*lz2)/4. + 2*omxi*lz2 -  
   (NC*omxi*lz2)/2. + 2*z*omxi*lz2 +  
   NC*z*omxi*lz2 -  
   (3*NCi2*omxi*lz2)/2. -  
   (3*z*NCi2*omxi*lz2)/2. +  
   (NCi*omxi*lz2)/2. -  
   z*NCi*omxi*lz2 -  
   (NC2*omxi*lz2)/2. -  
   (z*NC2*omxi*lz2)/2. +  
   (xi2*lz2)/2. - (NC*xi2*lz2)/4. -  
   (z*xi2*lz2)/2. + (NC*z*xi2*lz2)/2. -  
   (NCi2*xi2*lz2)/2. +  
   (z*NCi2*xi2*lz2)/2. +  
   (NCi*xi2*lz2)/4. -  
   (z*NCi*xi2*lz2)/2. - (xi*lz2)/2. +  
   (NC*xi*lz2)/4. - (NC*z*xi*lz2)/2. +  
   (NCi2*xi*lz2)/2. -  
   (NCi*xi*lz2)/4. +  
   (z*NCi*xi*lz2)/2. -  
   (opxi*lz2)/2. + (NC*opxi*lz2)/4. +  
   (z*opxi*lz2)/2. -  
   (NC*z*opxi*lz2)/2. +  
   (NCi2*opxi*lz2)/2. -  
   (z*NCi2*opxi*lz2)/2. -  
   (NCi*opxi*lz2)/4. +  
   (z*NCi*opxi*lz2)/2. -  
   (xi2*opxi*lz2)/2. +  
   (NC*xi2*opxi*lz2)/4. +  
   (z*xi2*opxi*lz2)/2. -  
   (NC*z*xi2*opxi*lz2)/2. +  
   (NCi2*xi2*opxi*lz2)/2. -  
   (z*NCi2*xi2*opxi*lz2)/2. -  
   (NCi*xi2*opxi*lz2)/4. +  
   (z*NCi*xi2*opxi*lz2)/2. +  
   (z*xi*opxi*lz2)/2. -  
   (z*NCi2*xi*opxi*lz2)/2. -  
   (omzi*lz2)/4. + (x*omzi*lz2)/4. +  
   (NCi2*omzi*lz2)/4. -  
   (x*NCi2*omzi*lz2)/4. -  
   (5*omxi*omzi*lz2)/2. +  
   (3*NCi2*omxi*omzi*lz2)/2. +  
   NC2*omxi*omzi*lz2 -  
   (3*zi*lz2)/4. - (x*zi*lz2)/4. +  
   (3*NCi2*zi*lz2)/4. +  
   (x*NCi2*zi*lz2)/4. +  
   omxi*zi*lz2 -  
   NCi2*omxi*zi*lz2 -  
   (xi2*zi*lz2)/4. +  
   (NCi2*xi2*zi*lz2)/4. +  
   (xi*zi*lz2)/4. -  
   (NCi2*xi*zi*lz2)/4. +  
   (opxi*zi*lz2)/4. -  
   (NCi2*opxi*zi*lz2)/4. +  
   (xi2*opxi*zi*lz2)/4. -  
   (NCi2*xi2*opxi*zi*lz2)/4. +  
   (omzi*zi*lz2)/2. +  
   (x*omzi*zi*lz2)/2. -  
   (NCi2*omzi*zi*lz2)/2. -  
   (x*NCi2*omzi*zi*lz2)/2. -  
   omxi*omzi*zi*lz2 +  
   NCi2*omxi*omzi*zi*lz2 +  
   2*NC*x*z2*lz2 - 2*x*NCi*z2*lz2 -  
   NC*omxi*z2*lz2 +  
   NCi*omxi*z2*lz2 -  
   NC*xi2*z2*lz2 +  
   NCi*xi2*z2*lz2 +  
   NC*opxi*z2*lz2 -  
   NCi*opxi*z2*lz2 +  
   NC*xi2*opxi*z2*lz2 -  
   NCi*xi2*opxi*z2*lz2 +  
   NC*xi*opxi*z2*lz2 -  
   NCi*xi*opxi*z2*lz2 +  
   (-1.5 + Li2z + x*z*Li2z - li2spec13 -  
      x*z*li2spec13 + li2spec15 +  
      x*z*li2spec15 + 3*lomx -  
      2*x*lomx + 2*x*z*lomx - (9*lx)/2. + 3*x*lx -  
      3*x*z*lx - 6*lomx*lx - 6*x*z*lomx*lx +  
      3*lomz - 2*x*lomz + 2*x*z*lomz +  
      5*lomx*lomz + 5*x*z*lomx*lomz -  
      6*lx*lomz - 6*x*z*lx*lomz + lx*lomxmz +  
      x*z*lx*lomxmz - lomz*lomxmz -  
      x*z*lomz*lomxmz + (3*lz)/2. - x*lz + x*z*lz +  
      2*lomx*lz + 2*x*z*lomx*lz - 3*lx*lz -  
      3*x*z*lx*lz + 2*lomz*lz + 2*x*z*lomz*lz +  
      (3*NC2)/2. - Li2z*NC2 - x*z*Li2z*NC2 +  
      li2spec13*NC2 + x*z*li2spec13*NC2 -  
      li2spec15*NC2 -  
      x*z*li2spec15*NC2 -  
      3*lomx*NC2 + 2*x*lomx*NC2 -  
      2*x*z*lomx*NC2 + (9*lx*NC2)/2. -  
      3*x*lx*NC2 + 3*x*z*lx*NC2 +  
      6*lomx*lx*NC2 + 6*x*z*lomx*lx*NC2 -  
      3*lomz*NC2 + 2*x*lomz*NC2 -  
      2*x*z*lomz*NC2 - 5*lomx*lomz*NC2 -  
      5*x*z*lomx*lomz*NC2 + 6*lx*lomz*NC2 +  
      6*x*z*lx*lomz*NC2 - lx*lomxmz*NC2 -  
      x*z*lx*lomxmz*NC2 +  
      lomz*lomxmz*NC2 +  
      x*z*lomz*lomxmz*NC2 - (3*lz*NC2)/2. +  
      x*lz*NC2 - x*z*lz*NC2 -  
      2*lomx*lz*NC2 - 2*x*z*lomx*lz*NC2 +  
      3*lx*lz*NC2 + 3*x*z*lx*lz*NC2 -  
      2*lomz*lz*NC2 - 2*x*z*lomz*lz*NC2 -  
      (2*pi2)/3. - (2*x*z*pi2)/3. + (2*NC2*pi2)/3. +  
      (2*x*z*NC2*pi2)/3. + omxi/2. -  
      (z*omxi)/2. - (Li2z*omxi)/2. -  
      (z*Li2z*omxi)/2. + (li2spec13*omxi)/2. +  
      (z*li2spec13*omxi)/2. -  
      (li2spec15*omxi)/2. -  
      (z*li2spec15*omxi)/2. -  
      2*z*lomx*omxi + 3*z*lx*omxi +  
      3*lomx*lx*omxi +  
      3*z*lomx*lx*omxi - 2*z*lomz*omxi -  
      (5*lomx*lomz*omxi)/2. -  
      (5*z*lomx*lomz*omxi)/2. +  
      3*lx*lomz*omxi +  
      3*z*lx*lomz*omxi -  
      (lx*lomxmz*omxi)/2. -  
      (z*lx*lomxmz*omxi)/2. +  
      (lomz*lomxmz*omxi)/2. +  
      (z*lomz*lomxmz*omxi)/2. -  
      z*lz*omxi - lomx*lz*omxi -  
      z*lomx*lz*omxi +  
      (3*lx*lz*omxi)/2. +  
      (3*z*lx*lz*omxi)/2. -  
      lomz*lz*omxi - z*lomz*lz*omxi -  
      (NC2*omxi)/2. + (z*NC2*omxi)/2. +  
      (Li2z*NC2*omxi)/2. +  
      (z*Li2z*NC2*omxi)/2. -  
      (li2spec13*NC2*omxi)/2. -  
      (z*li2spec13*NC2*omxi)/2. +  
      (li2spec15*NC2*omxi)/2. +  
      (z*li2spec15*NC2*omxi)/2. +  
      2*z*lomx*NC2*omxi -  
      3*z*lx*NC2*omxi -  
      3*lomx*lx*NC2*omxi -  
      3*z*lomx*lx*NC2*omxi +  
      2*z*lomz*NC2*omxi +  
      (5*lomx*lomz*NC2*omxi)/2. +  
      (5*z*lomx*lomz*NC2*omxi)/2. -  
      3*lx*lomz*NC2*omxi -  
      3*z*lx*lomz*NC2*omxi +  
      (lx*lomxmz*NC2*omxi)/2. +  
      (z*lx*lomxmz*NC2*omxi)/2. -  
      (lomz*lomxmz*NC2*omxi)/2. -  
      (z*lomz*lomxmz*NC2*omxi)/2. +  
      z*lz*NC2*omxi +  
      lomx*lz*NC2*omxi +  
      z*lomx*lz*NC2*omxi -  
      (3*lx*lz*NC2*omxi)/2. -  
      (3*z*lx*lz*NC2*omxi)/2. +  
      lomz*lz*NC2*omxi +  
      z*lomz*lz*NC2*omxi +  
      (pi2*omxi)/3. + (z*pi2*omxi)/3. -  
      (NC2*pi2*omxi)/3. -  
      (z*NC2*pi2*omxi)/3. + omzi/2. -  
      (x*omzi)/2. - (Li2z*omzi)/2. -  
      (x*Li2z*omzi)/2. + (li2spec13*omzi)/2. +  
      (x*li2spec13*omzi)/2. -  
      (li2spec15*omzi)/2. -  
      (x*li2spec15*omzi)/2. -  
      2*x*lomx*omzi + 3*x*lx*omzi +  
      3*lomx*lx*omzi +  
      3*x*lomx*lx*omzi - 2*x*lomz*omzi -  
      (5*lomx*lomz*omzi)/2. -  
      (5*x*lomx*lomz*omzi)/2. +  
      3*lx*lomz*omzi +  
      3*x*lx*lomz*omzi -  
      (lx*lomxmz*omzi)/2. -  
      (x*lx*lomxmz*omzi)/2. +  
      (lomz*lomxmz*omzi)/2. +  
      (x*lomz*lomxmz*omzi)/2. -  
      x*lz*omzi - lomx*lz*omzi -  
      x*lomx*lz*omzi +  
      (3*lx*lz*omzi)/2. +  
      (3*x*lx*lz*omzi)/2. -  
      lomz*lz*omzi - x*lomz*lz*omzi -  
      (NC2*omzi)/2. + (x*NC2*omzi)/2. +  
      (Li2z*NC2*omzi)/2. +  
      (x*Li2z*NC2*omzi)/2. -  
      (li2spec13*NC2*omzi)/2. -  
      (x*li2spec13*NC2*omzi)/2. +  
      (li2spec15*NC2*omzi)/2. +  
      (x*li2spec15*NC2*omzi)/2. +  
      2*x*lomx*NC2*omzi -  
      3*x*lx*NC2*omzi -  
      3*lomx*lx*NC2*omzi -  
      3*x*lomx*lx*NC2*omzi +  
      2*x*lomz*NC2*omzi +  
      (5*lomx*lomz*NC2*omzi)/2. +  
      (5*x*lomx*lomz*NC2*omzi)/2. -  
      3*lx*lomz*NC2*omzi -  
      3*x*lx*lomz*NC2*omzi +  
      (lx*lomxmz*NC2*omzi)/2. +  
      (x*lx*lomxmz*NC2*omzi)/2. -  
      (lomz*lomxmz*NC2*omzi)/2. -  
      (x*lomz*lomxmz*NC2*omzi)/2. +  
      x*lz*NC2*omzi +  
      lomx*lz*NC2*omzi +  
      x*lomx*lz*NC2*omzi -  
      (3*lx*lz*NC2*omzi)/2. -  
      (3*x*lx*lz*NC2*omzi)/2. +  
      lomz*lz*NC2*omzi +  
      x*lomz*lz*NC2*omzi +  
      (pi2*omzi)/3. + (x*pi2*omzi)/3. -  
      (NC2*pi2*omzi)/3. -  
      (x*NC2*pi2*omzi)/3. +  
      Li2z*omxi*omzi -  
      li2spec13*omxi*omzi +  
      li2spec15*omxi*omzi +  
      2*lomx*omxi*omzi -  
      3*lx*omxi*omzi -  
      6*lomx*lx*omxi*omzi +  
      2*lomz*omxi*omzi +  
      5*lomx*lomz*omxi*omzi -  
      6*lx*lomz*omxi*omzi +  
      lx*lomxmz*omxi*omzi -  
      lomz*lomxmz*omxi*omzi +  
      lz*omxi*omzi +  
      2*lomx*lz*omxi*omzi -  
      3*lx*lz*omxi*omzi +  
      2*lomz*lz*omxi*omzi -  
      Li2z*NC2*omxi*omzi +  
      li2spec13*NC2*omxi*omzi -  
      li2spec15*NC2*omxi* 
       omzi - 2*lomx*NC2*omxi*omzi +  
      3*lx*NC2*omxi*omzi +  
      6*lomx*lx*NC2*omxi*omzi -  
      2*lomz*NC2*omxi*omzi -  
      5*lomx*lomz*NC2*omxi*omzi +  
      6*lx*lomz*NC2*omxi*omzi -  
      lx*lomxmz*NC2*omxi*omzi +  
      lomz*lomxmz*NC2*omxi*omzi -  
      lz*NC2*omxi*omzi -  
      2*lomx*lz*NC2*omxi*omzi +  
      3*lx*lz*NC2*omxi*omzi -  
      2*lomz*lz*NC2*omxi*omzi -  
      (2*pi2*omxi*omzi)/3. +  
      (2*NC2*pi2*omxi*omzi)/3. +  
      2*lomx2 + 2*x*z*lomx2 -  
      2*NC2*lomx2 - 2*x*z*NC2*lomx2 -  
      omxi*lomx2 - z*omxi*lomx2 +  
      NC2*omxi*lomx2 +  
      z*NC2*omxi*lomx2 -  
      omzi*lomx2 - x*omzi*lomx2 +  
      NC2*omzi*lomx2 +  
      x*NC2*omzi*lomx2 +  
      2*omxi*omzi*lomx2 -  
      2*NC2*omxi*omzi*lomx2 +  
      (7*lx2)/2. + (7*x*z*lx2)/2. -  
      (7*NC2*lx2)/2. - (7*x*z*NC2*lx2)/2. -  
      (7*omxi*lx2)/4. -  
      (7*z*omxi*lx2)/4. +  
      (7*NC2*omxi*lx2)/4. +  
      (7*z*NC2*omxi*lx2)/4. -  
      (7*omzi*lx2)/4. -  
      (7*x*omzi*lx2)/4. +  
      (7*NC2*omzi*lx2)/4. +  
      (7*x*NC2*omzi*lx2)/4. +  
      (7*omxi*omzi*lx2)/2. -  
      (7*NC2*omxi*omzi*lx2)/2. +  
      (5*lomz2)/2. + (5*x*z*lomz2)/2. -  
      (5*NC2*lomz2)/2. -  
      (5*x*z*NC2*lomz2)/2. -  
      (5*omxi*lomz2)/4. -  
      (5*z*omxi*lomz2)/4. +  
      (5*NC2*omxi*lomz2)/4. +  
      (5*z*NC2*omxi*lomz2)/4. -  
      (5*omzi*lomz2)/4. -  
      (5*x*omzi*lomz2)/4. +  
      (5*NC2*omzi*lomz2)/4. +  
      (5*x*NC2*omzi*lomz2)/4. +  
      (5*omxi*omzi*lomz2)/2. -  
      (5*NC2*omxi*omzi*lomz2)/2. +  
      lz2/2. + (x*z*lz2)/2. -  
      (NC2*lz2)/2. - (x*z*NC2*lz2)/2. -  
      (omxi*lz2)/4. - (z*omxi*lz2)/4. +  
      (NC2*omxi*lz2)/4. +  
      (z*NC2*omxi*lz2)/4. -  
      (omzi*lz2)/4. - (x*omzi*lz2)/4. +  
      (NC2*omzi*lz2)/4. +  
      (x*NC2*omzi*lz2)/4. +  
      (omxi*omzi*lz2)/2. -  
      (NC2*omxi*omzi*lz2)/2.)*Tr1 +  
   (-1.5 + Li2z + x*z*Li2z + li2spec16 +  
      x*z*li2spec16 - li2spec17 -  
      x*z*li2spec17 + 3*lomx -  
      2*x*lomx + 2*x*z*lomx - (9*lx)/2. + 3*x*lx -  
      3*x*z*lx - 5*lomx*lx - 5*x*z*lomx*lx +  
      3*lomz - 2*x*lomz + 2*x*z*lomz +  
      4*lomx*lomz + 4*x*z*lomx*lomz -  
      5*lx*lomz - 5*x*z*lx*lomz + (3*lz)/2. -  
      x*lz + x*z*lz + 2*lomx*lz + 2*x*z*lomx*lz -  
      4*lx*lz - 4*x*z*lx*lz + 3*lomz*lz +  
      3*x*z*lomz*lz + lx*lmopxpz +  
      x*z*lx*lmopxpz - lomz*lmopxpz -  
      x*z*lomz*lmopxpz + (3*NC2)/2. - Li2z*NC2 -  
      x*z*Li2z*NC2 - li2spec16*NC2 -  
      x*z*li2spec16*NC2 +  
      li2spec17*NC2 +  
      x*z*li2spec17*NC2 -  
      3*lomx*NC2 + 2*x*lomx*NC2 -  
      2*x*z*lomx*NC2 + (9*lx*NC2)/2. -  
      3*x*lx*NC2 + 3*x*z*lx*NC2 +  
      5*lomx*lx*NC2 + 5*x*z*lomx*lx*NC2 -  
      3*lomz*NC2 + 2*x*lomz*NC2 -  
      2*x*z*lomz*NC2 - 4*lomx*lomz*NC2 -  
      4*x*z*lomx*lomz*NC2 + 5*lx*lomz*NC2 +  
      5*x*z*lx*lomz*NC2 - (3*lz*NC2)/2. +  
      x*lz*NC2 - x*z*lz*NC2 -  
      2*lomx*lz*NC2 - 2*x*z*lomx*lz*NC2 +  
      4*lx*lz*NC2 + 4*x*z*lx*lz*NC2 -  
      3*lomz*lz*NC2 - 3*x*z*lomz*lz*NC2 -  
      lx*lmopxpz*NC2 -  
      x*z*lx*lmopxpz*NC2 +  
      lomz*lmopxpz*NC2 +  
      x*z*lomz*lmopxpz*NC2 - (2*pi2)/3. -  
      (2*x*z*pi2)/3. + (2*NC2*pi2)/3. +  
      (2*x*z*NC2*pi2)/3. + omxi/2. -  
      (z*omxi)/2. - (Li2z*omxi)/2. -  
      (z*Li2z*omxi)/2. -  
      (li2spec16*omxi)/2. -  
      (z*li2spec16*omxi)/2. +  
      (li2spec17*omxi)/2. +  
      (z*li2spec17*omxi)/2. -  
      2*z*lomx*omxi + 3*z*lx*omxi +  
      (5*lomx*lx*omxi)/2. +  
      (5*z*lomx*lx*omxi)/2. -  
      2*z*lomz*omxi - 2*lomx*lomz*omxi -  
      2*z*lomx*lomz*omxi +  
      (5*lx*lomz*omxi)/2. +  
      (5*z*lx*lomz*omxi)/2. - z*lz*omxi -  
      lomx*lz*omxi - z*lomx*lz*omxi +  
      2*lx*lz*omxi + 2*z*lx*lz*omxi -  
      (3*lomz*lz*omxi)/2. -  
      (3*z*lomz*lz*omxi)/2. -  
      (lx*lmopxpz*omxi)/2. -  
      (z*lx*lmopxpz*omxi)/2. +  
      (lomz*lmopxpz*omxi)/2. +  
      (z*lomz*lmopxpz*omxi)/2. -  
      (NC2*omxi)/2. + (z*NC2*omxi)/2. +  
      (Li2z*NC2*omxi)/2. +  
      (z*Li2z*NC2*omxi)/2. +  
      (li2spec16*NC2*omxi)/2. +  
      (z*li2spec16*NC2*omxi)/2. -  
      (li2spec17*NC2*omxi)/2. -  
      (z*li2spec17*NC2*omxi)/ 
       2. + 2*z*lomx*NC2*omxi -  
      3*z*lx*NC2*omxi -  
      (5*lomx*lx*NC2*omxi)/2. -  
      (5*z*lomx*lx*NC2*omxi)/2. +  
      2*z*lomz*NC2*omxi +  
      2*lomx*lomz*NC2*omxi +  
      2*z*lomx*lomz*NC2*omxi -  
      (5*lx*lomz*NC2*omxi)/2. -  
      (5*z*lx*lomz*NC2*omxi)/2. +  
      z*lz*NC2*omxi +  
      lomx*lz*NC2*omxi +  
      z*lomx*lz*NC2*omxi -  
      2*lx*lz*NC2*omxi -  
      2*z*lx*lz*NC2*omxi +  
      (3*lomz*lz*NC2*omxi)/2. +  
      (3*z*lomz*lz*NC2*omxi)/2. +  
      (lx*lmopxpz*NC2*omxi)/2. +  
      (z*lx*lmopxpz*NC2*omxi)/2. -  
      (lomz*lmopxpz*NC2*omxi)/2. -  
      (z*lomz*lmopxpz*NC2*omxi)/2. +  
      (pi2*omxi)/3. + (z*pi2*omxi)/3. -  
      (NC2*pi2*omxi)/3. -  
      (z*NC2*pi2*omxi)/3. + omzi/2. -  
      (x*omzi)/2. - (Li2z*omzi)/2. -  
      (x*Li2z*omzi)/2. -  
      (li2spec16*omzi)/2. -  
      (x*li2spec16*omzi)/2. +  
      (li2spec17*omzi)/2. +  
      (x*li2spec17*omzi)/2. -  
      2*x*lomx*omzi + 3*x*lx*omzi +  
      (5*lomx*lx*omzi)/2. +  
      (5*x*lomx*lx*omzi)/2. -  
      2*x*lomz*omzi - 2*lomx*lomz*omzi -  
      2*x*lomx*lomz*omzi +  
      (5*lx*lomz*omzi)/2. +  
      (5*x*lx*lomz*omzi)/2. - x*lz*omzi -  
      lomx*lz*omzi - x*lomx*lz*omzi +  
      2*lx*lz*omzi + 2*x*lx*lz*omzi -  
      (3*lomz*lz*omzi)/2. -  
      (3*x*lomz*lz*omzi)/2. -  
      (lx*lmopxpz*omzi)/2. -  
      (x*lx*lmopxpz*omzi)/2. +  
      (lomz*lmopxpz*omzi)/2. +  
      (x*lomz*lmopxpz*omzi)/2. -  
      (NC2*omzi)/2. + (x*NC2*omzi)/2. +  
      (Li2z*NC2*omzi)/2. +  
      (x*Li2z*NC2*omzi)/2. +  
      (li2spec16*NC2*omzi)/2. +  
      (x*li2spec16*NC2*omzi)/2. -  
      (li2spec17*NC2*omzi)/2. -  
      (x*li2spec17*NC2*omzi)/ 
       2. + 2*x*lomx*NC2*omzi -  
      3*x*lx*NC2*omzi -  
      (5*lomx*lx*NC2*omzi)/2. -  
      (5*x*lomx*lx*NC2*omzi)/2. +  
      2*x*lomz*NC2*omzi +  
      2*lomx*lomz*NC2*omzi +  
      2*x*lomx*lomz*NC2*omzi -  
      (5*lx*lomz*NC2*omzi)/2. -  
      (5*x*lx*lomz*NC2*omzi)/2. +  
      x*lz*NC2*omzi +  
      lomx*lz*NC2*omzi +  
      x*lomx*lz*NC2*omzi -  
      2*lx*lz*NC2*omzi -  
      2*x*lx*lz*NC2*omzi +  
      (3*lomz*lz*NC2*omzi)/2. +  
      (3*x*lomz*lz*NC2*omzi)/2. +  
      (lx*lmopxpz*NC2*omzi)/2. +  
      (x*lx*lmopxpz*NC2*omzi)/2. -  
      (lomz*lmopxpz*NC2*omzi)/2. -  
      (x*lomz*lmopxpz*NC2*omzi)/2. +  
      (pi2*omzi)/3. + (x*pi2*omzi)/3. -  
      (NC2*pi2*omzi)/3. -  
      (x*NC2*pi2*omzi)/3. +  
      Li2z*omxi*omzi +  
      li2spec16*omxi*omzi -  
      li2spec17*omxi*omzi +  
      2*lomx*omxi*omzi -  
      3*lx*omxi*omzi -  
      5*lomx*lx*omxi*omzi +  
      2*lomz*omxi*omzi +  
      4*lomx*lomz*omxi*omzi -  
      5*lx*lomz*omxi*omzi +  
      lz*omxi*omzi +  
      2*lomx*lz*omxi*omzi -  
      4*lx*lz*omxi*omzi +  
      3*lomz*lz*omxi*omzi +  
      lx*lmopxpz*omxi*omzi -  
      lomz*lmopxpz*omxi*omzi -  
      Li2z*NC2*omxi*omzi -  
      li2spec16*NC2*omxi*omzi +  
      li2spec17*NC2*omxi* 
       omzi - 2*lomx*NC2*omxi*omzi +  
      3*lx*NC2*omxi*omzi +  
      5*lomx*lx*NC2*omxi*omzi -  
      2*lomz*NC2*omxi*omzi -  
      4*lomx*lomz*NC2*omxi*omzi +  
      5*lx*lomz*NC2*omxi*omzi -  
      lz*NC2*omxi*omzi -  
      2*lomx*lz*NC2*omxi*omzi +  
      4*lx*lz*NC2*omxi*omzi -  
      3*lomz*lz*NC2*omxi*omzi -  
      lx*lmopxpz*NC2*omxi*omzi +  
      lomz*lmopxpz*NC2*omxi*omzi -  
      (2*pi2*omxi*omzi)/3. +  
      (2*NC2*pi2*omxi*omzi)/3. +  
      2*lomx2 + 2*x*z*lomx2 -  
      2*NC2*lomx2 - 2*x*z*NC2*lomx2 -  
      omxi*lomx2 - z*omxi*lomx2 +  
      NC2*omxi*lomx2 +  
      z*NC2*omxi*lomx2 -  
      omzi*lomx2 - x*omzi*lomx2 +  
      NC2*omzi*lomx2 +  
      x*NC2*omzi*lomx2 +  
      2*omxi*omzi*lomx2 -  
      2*NC2*omxi*omzi*lomx2 +  
      3*lx2 + 3*x*z*lx2 - 3*NC2*lx2 -  
      3*x*z*NC2*lx2 - (3*omxi*lx2)/2. -  
      (3*z*omxi*lx2)/2. +  
      (3*NC2*omxi*lx2)/2. +  
      (3*z*NC2*omxi*lx2)/2. -  
      (3*omzi*lx2)/2. -  
      (3*x*omzi*lx2)/2. +  
      (3*NC2*omzi*lx2)/2. +  
      (3*x*NC2*omzi*lx2)/2. +  
      3*omxi*omzi*lx2 -  
      3*NC2*omxi*omzi*lx2 +  
      2*lomz2 + 2*x*z*lomz2 -  
      2*NC2*lomz2 - 2*x*z*NC2*lomz2 -  
      omxi*lomz2 - z*omxi*lomz2 +  
      NC2*omxi*lomz2 +  
      z*NC2*omxi*lomz2 -  
      omzi*lomz2 - x*omzi*lomz2 +  
      NC2*omzi*lomz2 +  
      x*NC2*omzi*lomz2 +  
      2*omxi*omzi*lomz2 -  
      2*NC2*omxi*omzi*lomz2 +  
      lz2/2. + (x*z*lz2)/2. -  
      (NC2*lz2)/2. - (x*z*NC2*lz2)/2. -  
      (omxi*lz2)/4. - (z*omxi*lz2)/4. +  
      (NC2*omxi*lz2)/4. +  
      (z*NC2*omxi*lz2)/4. -  
      (omzi*lz2)/4. - (x*omzi*lz2)/4. +  
      (NC2*omzi*lz2)/4. +  
      (x*NC2*omzi*lz2)/4. +  
      (omxi*omzi*lz2)/2. -  
      (NC2*omxi*omzi*lz2)/2.)*Tr2 +  
   (-3 - x*Li2z - x*z*Li2z + li2spec13 +  
      x*z*li2spec13 - li2spec14 -  
      x*z*li2spec14 -  
      li2spec11 + x*li2spec11 -  
      x*li2spec12 -  
      x*z*li2spec12 + 3*x*z*lomx -  
      6*x*z*lx - 4*lomx*lx - 5*x*lomx*lx -  
      9*x*z*lomx*lx + 3*x*z*lomz + 3*lomx*lomz +  
      2*x*lomx*lomz + 5*x*z*lomx*lomz -  
      5*lx*lomz - 4*x*lx*lomz - 9*x*z*lx*lomz +  
      lx*lomxmz + x*z*lx*lomxmz -  
      lomz*lomxmz - x*z*lomz*lomxmz +  
      lx*lxmz + x*lx*lxmz + 2*x*z*lx*lxmz +  
      4*x*z*lz + 2*lomx*lz + 4*x*lomx*lz +  
      6*x*z*lomx*lz - 5*lx*lz - 5*x*lx*lz -  
      10*x*z*lx*lz + 4*lomz*lz + 2*x*lomz*lz +  
      6*x*z*lomz*lz - lxmz*lz - x*lxmz*lz -  
      2*x*z*lxmz*lz + 3*NCi2 + x*Li2z*NCi2 +  
      x*z*Li2z*NCi2 - li2spec13*NCi2 -  
      x*z*li2spec13*NCi2 +  
      li2spec14*NCi2 +  
      x*z*li2spec14*NCi2 +  
      li2spec11*NCi2 -  
      x*li2spec11*NCi2 +  
      x*li2spec12*NCi2 +  
      x*z*li2spec12*NCi2 -  
      3*x*z*lomx*NCi2 + 6*x*z*lx*NCi2 +  
      4*lomx*lx*NCi2 + 5*x*lomx*lx*NCi2 +  
      9*x*z*lomx*lx*NCi2 - 3*x*z*lomz*NCi2 -  
      3*lomx*lomz*NCi2 -  
      2*x*lomx*lomz*NCi2 -  
      5*x*z*lomx*lomz*NCi2 +  
      5*lx*lomz*NCi2 + 4*x*lx*lomz*NCi2 +  
      9*x*z*lx*lomz*NCi2 - lx*lomxmz*NCi2 -  
      x*z*lx*lomxmz*NCi2 +  
      lomz*lomxmz*NCi2 +  
      x*z*lomz*lomxmz*NCi2 -  
      lx*lxmz*NCi2 - x*lx*lxmz*NCi2 -  
      2*x*z*lx*lxmz*NCi2 - 4*x*z*lz*NCi2 -  
      2*lomx*lz*NCi2 - 4*x*lomx*lz*NCi2 -  
      6*x*z*lomx*lz*NCi2 + 5*lx*lz*NCi2 +  
      5*x*lx*lz*NCi2 + 10*x*z*lx*lz*NCi2 -  
      4*lomz*lz*NCi2 - 2*x*lomz*lz*NCi2 -  
      6*x*z*lomz*lz*NCi2 + lxmz*lz*NCi2 +  
      x*lxmz*lz*NCi2 + 2*x*z*lxmz*lz*NCi2 -  
      pi2/6. - (x*pi2)/2. - (2*x*z*pi2)/3. +  
      (NCi2*pi2)/6. + (x*NCi2*pi2)/2. +  
      (2*x*z*NCi2*pi2)/3. + 2*omxi - 2*z*omxi +  
      (Li2z*omxi)/2. + (z*Li2z*omxi)/2. -  
      li2spec13*omxi -  
      z*li2spec13*omxi +  
      li2spec14*omxi +  
      z*li2spec14*omxi +  
      (li2spec11*omxi)/2. +  
      (z*li2spec11*omxi)/2. +  
      (li2spec12*omxi)/2. +  
      (z*li2spec12*omxi)/2. -  
      3*z*lomx*omxi + 6*z*lx*omxi +  
      (13*lomx*lx*omxi)/2. +  
      (13*z*lomx*lx*omxi)/2. -  
      3*z*lomz*omxi - 4*lomx*lomz*omxi -  
      4*z*lomx*lomz*omxi +  
      7*lx*lomz*omxi +  
      7*z*lx*lomz*omxi -  
      lx*lomxmz*omxi -  
      z*lx*lomxmz*omxi +  
      lomz*lomxmz*omxi +  
      z*lomz*lomxmz*omxi -  
      (3*lx*lxmz*omxi)/2. -  
      (3*z*lx*lxmz*omxi)/2. - 4*z*lz*omxi -  
      4*lomx*lz*omxi -  
      4*z*lomx*lz*omxi +  
      (15*lx*lz*omxi)/2. +  
      (15*z*lx*lz*omxi)/2. -  
      5*lomz*lz*omxi -  
      5*z*lomz*lz*omxi +  
      (3*lxmz*lz*omxi)/2. +  
      (3*z*lxmz*lz*omxi)/2. - 2*NCi2*omxi +  
      2*z*NCi2*omxi - (Li2z*NCi2*omxi)/2. -  
      (z*Li2z*NCi2*omxi)/2. +  
      li2spec13*NCi2*omxi +  
      z*li2spec13*NCi2*omxi -  
      li2spec14*NCi2*omxi -  
      z*li2spec14*NCi2*omxi -  
      (li2spec11*NCi2*omxi)/2. -  
      (z*li2spec11*NCi2*omxi)/2. -  
      (li2spec12*NCi2*omxi)/2. -  
      (z*li2spec12*NCi2*omxi)/ 
       2. + 3*z*lomx*NCi2*omxi -  
      6*z*lx*NCi2*omxi -  
      (13*lomx*lx*NCi2*omxi)/2. -  
      (13*z*lomx*lx*NCi2*omxi)/2. +  
      3*z*lomz*NCi2*omxi +  
      4*lomx*lomz*NCi2*omxi +  
      4*z*lomx*lomz*NCi2*omxi -  
      7*lx*lomz*NCi2*omxi -  
      7*z*lx*lomz*NCi2*omxi +  
      lx*lomxmz*NCi2*omxi +  
      z*lx*lomxmz*NCi2*omxi -  
      lomz*lomxmz*NCi2*omxi -  
      z*lomz*lomxmz*NCi2*omxi +  
      (3*lx*lxmz*NCi2*omxi)/2. +  
      (3*z*lx*lxmz*NCi2*omxi)/2. +  
      4*z*lz*NCi2*omxi +  
      4*lomx*lz*NCi2*omxi +  
      4*z*lomx*lz*NCi2*omxi -  
      (15*lx*lz*NCi2*omxi)/2. -  
      (15*z*lx*lz*NCi2*omxi)/2. +  
      5*lomz*lz*NCi2*omxi +  
      5*z*lomz*lz*NCi2*omxi -  
      (3*lxmz*lz*NCi2*omxi)/2. -  
      (3*z*lxmz*lz*NCi2*omxi)/2. +  
      (5*pi2*omxi)/12. + (5*z*pi2*omxi)/12. -  
      (5*NCi2*pi2*omxi)/12. -  
      (5*z*NCi2*pi2*omxi)/12. + 2*omzi -  
      2*x*omzi + x*Li2z*omzi -  
      (li2spec13*omzi)/2. -  
      (x*li2spec13*omzi)/2. +  
      (li2spec14*omzi)/2. +  
      (x*li2spec14*omzi)/2. +  
      (li2spec11*omzi)/2. -  
      (x*li2spec11*omzi)/2. +  
      x*li2spec12*omzi -  
      3*x*lomx*omzi + 6*x*lx*omzi +  
      2*lomx*lx*omzi +  
      7*x*lomx*lx*omzi - 3*x*lomz*omzi -  
      (3*lomx*lomz*omzi)/2. -  
      (7*x*lomx*lomz*omzi)/2. +  
      (5*lx*lomz*omzi)/2. +  
      (13*x*lx*lomz*omzi)/2. -  
      (lx*lomxmz*omzi)/2. -  
      (x*lx*lomxmz*omzi)/2. +  
      (lomz*lomxmz*omzi)/2. +  
      (x*lomz*lomxmz*omzi)/2. -  
      (lx*lxmz*omzi)/2. -  
      (3*x*lx*lxmz*omzi)/2. - 4*x*lz*omzi -  
      lomx*lz*omzi - 5*x*lomx*lz*omzi +  
      (5*lx*lz*omzi)/2. +  
      (15*x*lx*lz*omzi)/2. -  
      2*lomz*lz*omzi -  
      4*x*lomz*lz*omzi +  
      (lxmz*lz*omzi)/2. +  
      (3*x*lxmz*lz*omzi)/2. - 2*NCi2*omzi +  
      2*x*NCi2*omzi - x*Li2z*NCi2*omzi +  
      (li2spec13*NCi2*omzi)/2. +  
      (x*li2spec13*NCi2*omzi)/2. -  
      (li2spec14*NCi2*omzi)/2. -  
      (x*li2spec14*NCi2*omzi)/ 
       2. - (li2spec11*NCi2*omzi)/2. +  
      (x*li2spec11*NCi2*omzi)/2. -  
      x*li2spec12*NCi2*omzi +  
      3*x*lomx*NCi2*omzi -  
      6*x*lx*NCi2*omzi -  
      2*lomx*lx*NCi2*omzi -  
      7*x*lomx*lx*NCi2*omzi +  
      3*x*lomz*NCi2*omzi +  
      (3*lomx*lomz*NCi2*omzi)/2. +  
      (7*x*lomx*lomz*NCi2*omzi)/2. -  
      (5*lx*lomz*NCi2*omzi)/2. -  
      (13*x*lx*lomz*NCi2*omzi)/2. +  
      (lx*lomxmz*NCi2*omzi)/2. +  
      (x*lx*lomxmz*NCi2*omzi)/2. -  
      (lomz*lomxmz*NCi2*omzi)/2. -  
      (x*lomz*lomxmz*NCi2*omzi)/2. +  
      (lx*lxmz*NCi2*omzi)/2. +  
      (3*x*lx*lxmz*NCi2*omzi)/2. +  
      4*x*lz*NCi2*omzi +  
      lomx*lz*NCi2*omzi +  
      5*x*lomx*lz*NCi2*omzi -  
      (5*lx*lz*NCi2*omzi)/2. -  
      (15*x*lx*lz*NCi2*omzi)/2. +  
      2*lomz*lz*NCi2*omzi +  
      4*x*lomz*lz*NCi2*omzi -  
      (lxmz*lz*NCi2*omzi)/2. -  
      (3*x*lxmz*lz*NCi2*omzi)/2. +  
      (pi2*omzi)/12. + (7*x*pi2*omzi)/12. -  
      (NCi2*pi2*omzi)/12. -  
      (7*x*NCi2*pi2*omzi)/12. -  
      Li2z*omxi*omzi +  
      li2spec13*omxi*omzi -  
      li2spec14*omxi*omzi -  
      li2spec12*omxi*omzi +  
      3*lomx*omxi*omzi -  
      6*lx*omxi*omzi -  
      9*lomx*lx*omxi*omzi +  
      3*lomz*omxi*omzi +  
      5*lomx*lomz*omxi*omzi -  
      9*lx*lomz*omxi*omzi +  
      lx*lomxmz*omxi*omzi -  
      lomz*lomxmz*omxi*omzi +  
      2*lx*lxmz*omxi*omzi +  
      4*lz*omxi*omzi +  
      6*lomx*lz*omxi*omzi -  
      10*lx*lz*omxi*omzi +  
      6*lomz*lz*omxi*omzi -  
      2*lxmz*lz*omxi*omzi +  
      Li2z*NCi2*omxi*omzi -  
      li2spec13*NCi2*omxi*omzi +  
      li2spec14*NCi2*omxi* 
       omzi + li2spec12*NCi2* 
       omxi*omzi -  
      3*lomx*NCi2*omxi*omzi +  
      6*lx*NCi2*omxi*omzi +  
      9*lomx*lx*NCi2*omxi*omzi -  
      3*lomz*NCi2*omxi*omzi -  
      5*lomx*lomz*NCi2*omxi*omzi +  
      9*lx*lomz*NCi2*omxi*omzi -  
      lx*lomxmz*NCi2*omxi*omzi +  
      lomz*lomxmz*NCi2*omxi*omzi -  
      2*lx*lxmz*NCi2*omxi*omzi -  
      4*lz*NCi2*omxi*omzi -  
      6*lomx*lz*NCi2*omxi*omzi +  
      10*lx*lz*NCi2*omxi*omzi -  
      6*lomz*lz*NCi2*omxi*omzi +  
      2*lxmz*lz*NCi2*omxi*omzi -  
      (2*pi2*omxi*omzi)/3. +  
      (2*NCi2*pi2*omxi*omzi)/3. +  
      lomx2/2. + 2*x*lomx2 +  
      (5*x*z*lomx2)/2. - (NCi2*lomx2)/2. -  
      2*x*NCi2*lomx2 -  
      (5*x*z*NCi2*lomx2)/2. -  
      (3*omxi*lomx2)/2. -  
      (3*z*omxi*lomx2)/2. +  
      (3*NCi2*omxi*lomx2)/2. +  
      (3*z*NCi2*omxi*lomx2)/2. -  
      (omzi*lomx2)/4. -  
      (9*x*omzi*lomx2)/4. +  
      (NCi2*omzi*lomx2)/4. +  
      (9*x*NCi2*omzi*lomx2)/4. +  
      (5*omxi*omzi*lomx2)/2. -  
      (5*NCi2*omxi*omzi*lomx2)/2. +  
      3*lx2 + 3*x*lx2 + 6*x*z*lx2 -  
      3*NCi2*lx2 - 3*x*NCi2*lx2 -  
      6*x*z*NCi2*lx2 - (9*omxi*lx2)/2. -  
      (9*z*omxi*lx2)/2. +  
      (9*NCi2*omxi*lx2)/2. +  
      (9*z*NCi2*omxi*lx2)/2. -  
      (3*omzi*lx2)/2. -  
      (9*x*omzi*lx2)/2. +  
      (3*NCi2*omzi*lx2)/2. +  
      (9*x*NCi2*omzi*lx2)/2. +  
      6*omxi*omzi*lx2 -  
      6*NCi2*omxi*omzi*lx2 +  
      (3*lomz2)/2. + (x*lomz2)/2. +  
      2*x*z*lomz2 - (3*NCi2*lomz2)/2. -  
      (x*NCi2*lomz2)/2. -  
      2*x*z*NCi2*lomz2 -  
      (7*omxi*lomz2)/4. -  
      (7*z*omxi*lomz2)/4. +  
      (7*NCi2*omxi*lomz2)/4. +  
      (7*z*NCi2*omxi*lomz2)/4. -  
      (3*omzi*lomz2)/4. -  
      (5*x*omzi*lomz2)/4. +  
      (3*NCi2*omzi*lomz2)/4. +  
      (5*x*NCi2*omzi*lomz2)/4. +  
      2*omxi*omzi*lomz2 -  
      2*NCi2*omxi*omzi*lomz2 +  
      2*lz2 + 2*x*lz2 + 4*x*z*lz2 -  
      2*NCi2*lz2 - 2*x*NCi2*lz2 -  
      4*x*z*NCi2*lz2 - 3*omxi*lz2 -  
      3*z*omxi*lz2 +  
      3*NCi2*omxi*lz2 +  
      3*z*NCi2*omxi*lz2 -  
      omzi*lz2 - 3*x*omzi*lz2 +  
      NCi2*omzi*lz2 +  
      3*x*NCi2*omzi*lz2 +  
      4*omxi*omzi*lz2 -  
      4*NCi2*omxi*omzi*lz2)*Tu1 +  
   (-3 - x*Li2z - x*z*Li2z - li2spec11 +  
      x*li2spec11 -  
      x*li2spec12 -  
      x*z*li2spec12 - li2spec16 -  
      x*z*li2spec16 + li2spec18 +  
      x*z*li2spec18 + 3*x*z*lomx -  
      6*x*z*lx - 3*lomx*lx - 5*x*lomx*lx -  
      8*x*z*lomx*lx + 3*x*z*lomz + 2*lomx*lomz +  
      2*x*lomx*lomz + 4*x*z*lomx*lomz -  
      6*lx*lomz - 4*x*lx*lomz -  
      10*x*z*lx*lomz + lx*lxmz + x*lx*lxmz +  
      2*x*z*lx*lxmz + 4*x*z*lz + 2*lomx*lz +  
      4*x*lomx*lz + 6*x*z*lomx*lz - 6*lx*lz -  
      5*x*lx*lz - 11*x*z*lx*lz + 5*lomz*lz +  
      2*x*lomz*lz + 7*x*z*lomz*lz - lxmz*lz -  
      x*lxmz*lz - 2*x*z*lxmz*lz +  
      lx*lmopxpz + x*z*lx*lmopxpz -  
      lomz*lmopxpz - x*z*lomz*lmopxpz +  
      3*NCi2 + x*Li2z*NCi2 + x*z*Li2z*NCi2 +  
      li2spec11*NCi2 -  
      x*li2spec11*NCi2 +  
      x*li2spec12*NCi2 +  
      x*z*li2spec12*NCi2 +  
      li2spec16*NCi2 +  
      x*z*li2spec16*NCi2 -  
      li2spec18*NCi2 -  
      x*z*li2spec18*NCi2 -  
      3*x*z*lomx*NCi2 + 6*x*z*lx*NCi2 +  
      3*lomx*lx*NCi2 + 5*x*lomx*lx*NCi2 +  
      8*x*z*lomx*lx*NCi2 - 3*x*z*lomz*NCi2 -  
      2*lomx*lomz*NCi2 -  
      2*x*lomx*lomz*NCi2 -  
      4*x*z*lomx*lomz*NCi2 +  
      6*lx*lomz*NCi2 + 4*x*lx*lomz*NCi2 +  
      10*x*z*lx*lomz*NCi2 - lx*lxmz*NCi2 -  
      x*lx*lxmz*NCi2 - 2*x*z*lx*lxmz*NCi2 -  
      4*x*z*lz*NCi2 - 2*lomx*lz*NCi2 -  
      4*x*lomx*lz*NCi2 - 6*x*z*lomx*lz*NCi2 +  
      6*lx*lz*NCi2 + 5*x*lx*lz*NCi2 +  
      11*x*z*lx*lz*NCi2 - 5*lomz*lz*NCi2 -  
      2*x*lomz*lz*NCi2 - 7*x*z*lomz*lz*NCi2 +  
      lxmz*lz*NCi2 + x*lxmz*lz*NCi2 +  
      2*x*z*lxmz*lz*NCi2 -  
      lx*lmopxpz*NCi2 -  
      x*z*lx*lmopxpz*NCi2 +  
      lomz*lmopxpz*NCi2 +  
      x*z*lomz*lmopxpz*NCi2 - pi2/6. -  
      (x*pi2)/2. - (2*x*z*pi2)/3. + (NCi2*pi2)/6. +  
      (x*NCi2*pi2)/2. + (2*x*z*NCi2*pi2)/3. +  
      2*omxi - 2*z*omxi + (Li2z*omxi)/2. +  
      (z*Li2z*omxi)/2. +  
      (li2spec11*omxi)/2. +  
      (z*li2spec11*omxi)/2. +  
      (li2spec12*omxi)/2. +  
      (z*li2spec12*omxi)/2. +  
      li2spec16*omxi +  
      z*li2spec16*omxi -  
      li2spec18*omxi -  
      z*li2spec18*omxi -  
      3*z*lomx*omxi + 6*z*lx*omxi +  
      (11*lomx*lx*omxi)/2. +  
      (11*z*lomx*lx*omxi)/2. -  
      3*z*lomz*omxi - 3*lomx*lomz*omxi -  
      3*z*lomx*lomz*omxi +  
      8*lx*lomz*omxi +  
      8*z*lx*lomz*omxi -  
      (3*lx*lxmz*omxi)/2. -  
      (3*z*lx*lxmz*omxi)/2. - 4*z*lz*omxi -  
      4*lomx*lz*omxi -  
      4*z*lomx*lz*omxi +  
      (17*lx*lz*omxi)/2. +  
      (17*z*lx*lz*omxi)/2. -  
      6*lomz*lz*omxi -  
      6*z*lomz*lz*omxi +  
      (3*lxmz*lz*omxi)/2. +  
      (3*z*lxmz*lz*omxi)/2. -  
      lx*lmopxpz*omxi -  
      z*lx*lmopxpz*omxi +  
      lomz*lmopxpz*omxi +  
      z*lomz*lmopxpz*omxi -  
      2*NCi2*omxi + 2*z*NCi2*omxi -  
      (Li2z*NCi2*omxi)/2. -  
      (z*Li2z*NCi2*omxi)/2. -  
      (li2spec11*NCi2*omxi)/2. -  
      (z*li2spec11*NCi2*omxi)/2. -  
      (li2spec12*NCi2*omxi)/2. -  
      (z*li2spec12*NCi2*omxi)/ 
       2. - li2spec16*NCi2*omxi -  
      z*li2spec16*NCi2*omxi +  
      li2spec18*NCi2*omxi +  
      z*li2spec18*NCi2*omxi +  
      3*z*lomx*NCi2*omxi -  
      6*z*lx*NCi2*omxi -  
      (11*lomx*lx*NCi2*omxi)/2. -  
      (11*z*lomx*lx*NCi2*omxi)/2. +  
      3*z*lomz*NCi2*omxi +  
      3*lomx*lomz*NCi2*omxi +  
      3*z*lomx*lomz*NCi2*omxi -  
      8*lx*lomz*NCi2*omxi -  
      8*z*lx*lomz*NCi2*omxi +  
      (3*lx*lxmz*NCi2*omxi)/2. +  
      (3*z*lx*lxmz*NCi2*omxi)/2. +  
      4*z*lz*NCi2*omxi +  
      4*lomx*lz*NCi2*omxi +  
      4*z*lomx*lz*NCi2*omxi -  
      (17*lx*lz*NCi2*omxi)/2. -  
      (17*z*lx*lz*NCi2*omxi)/2. +  
      6*lomz*lz*NCi2*omxi +  
      6*z*lomz*lz*NCi2*omxi -  
      (3*lxmz*lz*NCi2*omxi)/2. -  
      (3*z*lxmz*lz*NCi2*omxi)/2. +  
      lx*lmopxpz*NCi2*omxi +  
      z*lx*lmopxpz*NCi2*omxi -  
      lomz*lmopxpz*NCi2*omxi -  
      z*lomz*lmopxpz*NCi2*omxi +  
      (5*pi2*omxi)/12. + (5*z*pi2*omxi)/12. -  
      (5*NCi2*pi2*omxi)/12. -  
      (5*z*NCi2*pi2*omxi)/12. + 2*omzi -  
      2*x*omzi + x*Li2z*omzi +  
      (li2spec11*omzi)/2. -  
      (x*li2spec11*omzi)/2. +  
      x*li2spec12*omzi +  
      (li2spec16*omzi)/2. +  
      (x*li2spec16*omzi)/2. -  
      (li2spec18*omzi)/2. -  
      (x*li2spec18*omzi)/2. -  
      3*x*lomx*omzi + 6*x*lx*omzi +  
      (3*lomx*lx*omzi)/2. +  
      (13*x*lomx*lx*omzi)/2. -  
      3*x*lomz*omzi - lomx*lomz*omzi -  
      3*x*lomx*lomz*omzi +  
      3*lx*lomz*omzi +  
      7*x*lx*lomz*omzi -  
      (lx*lxmz*omzi)/2. -  
      (3*x*lx*lxmz*omzi)/2. - 4*x*lz*omzi -  
      lomx*lz*omzi - 5*x*lomx*lz*omzi +  
      3*lx*lz*omzi + 8*x*lx*lz*omzi -  
      (5*lomz*lz*omzi)/2. -  
      (9*x*lomz*lz*omzi)/2. +  
      (lxmz*lz*omzi)/2. +  
      (3*x*lxmz*lz*omzi)/2. -  
      (lx*lmopxpz*omzi)/2. -  
      (x*lx*lmopxpz*omzi)/2. +  
      (lomz*lmopxpz*omzi)/2. +  
      (x*lomz*lmopxpz*omzi)/2. -  
      2*NCi2*omzi + 2*x*NCi2*omzi -  
      x*Li2z*NCi2*omzi -  
      (li2spec11*NCi2*omzi)/2. +  
      (x*li2spec11*NCi2*omzi)/2. -  
      x*li2spec12*NCi2*omzi -  
      (li2spec16*NCi2*omzi)/2. -  
      (x*li2spec16*NCi2*omzi)/2. +  
      (li2spec18*NCi2*omzi)/2. +  
      (x*li2spec18*NCi2*omzi)/ 
       2. + 3*x*lomx*NCi2*omzi -  
      6*x*lx*NCi2*omzi -  
      (3*lomx*lx*NCi2*omzi)/2. -  
      (13*x*lomx*lx*NCi2*omzi)/2. +  
      3*x*lomz*NCi2*omzi +  
      lomx*lomz*NCi2*omzi +  
      3*x*lomx*lomz*NCi2*omzi -  
      3*lx*lomz*NCi2*omzi -  
      7*x*lx*lomz*NCi2*omzi +  
      (lx*lxmz*NCi2*omzi)/2. +  
      (3*x*lx*lxmz*NCi2*omzi)/2. +  
      4*x*lz*NCi2*omzi +  
      lomx*lz*NCi2*omzi +  
      5*x*lomx*lz*NCi2*omzi -  
      3*lx*lz*NCi2*omzi -  
      8*x*lx*lz*NCi2*omzi +  
      (5*lomz*lz*NCi2*omzi)/2. +  
      (9*x*lomz*lz*NCi2*omzi)/2. -  
      (lxmz*lz*NCi2*omzi)/2. -  
      (3*x*lxmz*lz*NCi2*omzi)/2. +  
      (lx*lmopxpz*NCi2*omzi)/2. +  
      (x*lx*lmopxpz*NCi2*omzi)/2. -  
      (lomz*lmopxpz*NCi2*omzi)/2. -  
      (x*lomz*lmopxpz*NCi2*omzi)/2. +  
      (pi2*omzi)/12. + (7*x*pi2*omzi)/12. -  
      (NCi2*pi2*omzi)/12. -  
      (7*x*NCi2*pi2*omzi)/12. -  
      Li2z*omxi*omzi -  
      li2spec12*omxi*omzi -  
      li2spec16*omxi*omzi +  
      li2spec18*omxi*omzi +  
      3*lomx*omxi*omzi -  
      6*lx*omxi*omzi -  
      8*lomx*lx*omxi*omzi +  
      3*lomz*omxi*omzi +  
      4*lomx*lomz*omxi*omzi -  
      10*lx*lomz*omxi*omzi +  
      2*lx*lxmz*omxi*omzi +  
      4*lz*omxi*omzi +  
      6*lomx*lz*omxi*omzi -  
      11*lx*lz*omxi*omzi +  
      7*lomz*lz*omxi*omzi -  
      2*lxmz*lz*omxi*omzi +  
      lx*lmopxpz*omxi*omzi -  
      lomz*lmopxpz*omxi*omzi +  
      Li2z*NCi2*omxi*omzi +  
      li2spec12*NCi2*omxi* 
       omzi + li2spec16*NCi2*omxi* 
       omzi - li2spec18*NCi2* 
       omxi*omzi -  
      3*lomx*NCi2*omxi*omzi +  
      6*lx*NCi2*omxi*omzi +  
      8*lomx*lx*NCi2*omxi*omzi -  
      3*lomz*NCi2*omxi*omzi -  
      4*lomx*lomz*NCi2*omxi*omzi +  
      10*lx*lomz*NCi2*omxi*omzi -  
      2*lx*lxmz*NCi2*omxi*omzi -  
      4*lz*NCi2*omxi*omzi -  
      6*lomx*lz*NCi2*omxi*omzi +  
      11*lx*lz*NCi2*omxi*omzi -  
      7*lomz*lz*NCi2*omxi*omzi +  
      2*lxmz*lz*NCi2*omxi*omzi -  
      lx*lmopxpz*NCi2*omxi*omzi +  
      lomz*lmopxpz*NCi2*omxi*omzi -  
      (2*pi2*omxi*omzi)/3. +  
      (2*NCi2*pi2*omxi*omzi)/3. +  
      lomx2/2. + 2*x*lomx2 +  
      (5*x*z*lomx2)/2. - (NCi2*lomx2)/2. -  
      2*x*NCi2*lomx2 -  
      (5*x*z*NCi2*lomx2)/2. -  
      (3*omxi*lomx2)/2. -  
      (3*z*omxi*lomx2)/2. +  
      (3*NCi2*omxi*lomx2)/2. +  
      (3*z*NCi2*omxi*lomx2)/2. -  
      (omzi*lomx2)/4. -  
      (9*x*omzi*lomx2)/4. +  
      (NCi2*omzi*lomx2)/4. +  
      (9*x*NCi2*omzi*lomx2)/4. +  
      (5*omxi*omzi*lomx2)/2. -  
      (5*NCi2*omxi*omzi*lomx2)/2. +  
      (7*lx2)/2. + 3*x*lx2 + (13*x*z*lx2)/2. -  
      (7*NCi2*lx2)/2. - 3*x*NCi2*lx2 -  
      (13*x*z*NCi2*lx2)/2. - 5*omxi*lx2 -  
      5*z*omxi*lx2 +  
      5*NCi2*omxi*lx2 +  
      5*z*NCi2*omxi*lx2 -  
      (7*omzi*lx2)/4. -  
      (19*x*omzi*lx2)/4. +  
      (7*NCi2*omzi*lx2)/4. +  
      (19*x*NCi2*omzi*lx2)/4. +  
      (13*omxi*omzi*lx2)/2. -  
      (13*NCi2*omxi*omzi*lx2)/2. +  
      2*lomz2 + (x*lomz2)/2. +  
      (5*x*z*lomz2)/2. - 2*NCi2*lomz2 -  
      (x*NCi2*lomz2)/2. -  
      (5*x*z*NCi2*lomz2)/2. -  
      (9*omxi*lomz2)/4. -  
      (9*z*omxi*lomz2)/4. +  
      (9*NCi2*omxi*lomz2)/4. +  
      (9*z*NCi2*omxi*lomz2)/4. -  
      omzi*lomz2 -  
      (3*x*omzi*lomz2)/2. +  
      NCi2*omzi*lomz2 +  
      (3*x*NCi2*omzi*lomz2)/2. +  
      (5*omxi*omzi*lomz2)/2. -  
      (5*NCi2*omxi*omzi*lomz2)/2. +  
      2*lz2 + 2*x*lz2 + 4*x*z*lz2 -  
      2*NCi2*lz2 - 2*x*NCi2*lz2 -  
      4*x*z*NCi2*lz2 - 3*omxi*lz2 -  
      3*z*omxi*lz2 +  
      3*NCi2*omxi*lz2 +  
      3*z*NCi2*omxi*lz2 -  
      omzi*lz2 - 3*x*omzi*lz2 +  
      NCi2*omzi*lz2 +  
      3*x*NCi2*omzi*lz2 +  
      4*omxi*omzi*lz2 -  
      4*NCi2*omxi*omzi*lz2)*Tu2 +  
   (-3 - x*Li2z - x*z*Li2z + li2spec9 -  
      x*li2spec9 + li2spec13 +  
      x*z*li2spec13 + x*li2spec10 +  
      x*z*li2spec10 +  
      li2spec18 +  
      x*z*li2spec18 + 3*x*z*lomx -  
      6*x*z*lx - 3*lomx*lx - 6*x*lomx*lx -  
      9*x*z*lomx*lx + 3*x*z*lomz + lomx*lomz +  
      2*x*lomx*lomz + 3*x*z*lomx*lomz -  
      6*lx*lomz - 3*x*lx*lomz - 9*x*z*lx*lomz +  
      lx*lomxmz + x*z*lx*lomxmz -  
      lomz*lomxmz - x*z*lomz*lomxmz +  
      4*x*z*lz + lomx*lz + 5*x*lomx*lz +  
      6*x*z*lomx*lz - 6*lx*lz - 6*x*lx*lz -  
      12*x*z*lx*lz + 5*lomz*lz + x*lomz*lz +  
      6*x*z*lomz*lz + lx*lmxpz + x*lx*lmxpz +  
      2*x*z*lx*lmxpz - lz*lmxpz - x*lz*lmxpz -  
      2*x*z*lz*lmxpz + 3*NCi2 + x*Li2z*NCi2 +  
      x*z*Li2z*NCi2 - li2spec9*NCi2 +  
      x*li2spec9*NCi2 -  
      li2spec13*NCi2 - x*z*li2spec13*NCi2 -  
      x*li2spec10*NCi2 -  
      x*z*li2spec10*NCi2 -  
      li2spec18*NCi2 -  
      x*z*li2spec18*NCi2 -  
      3*x*z*lomx*NCi2 + 6*x*z*lx*NCi2 +  
      3*lomx*lx*NCi2 + 6*x*lomx*lx*NCi2 +  
      9*x*z*lomx*lx*NCi2 - 3*x*z*lomz*NCi2 -  
      lomx*lomz*NCi2 -  
      2*x*lomx*lomz*NCi2 -  
      3*x*z*lomx*lomz*NCi2 +  
      6*lx*lomz*NCi2 + 3*x*lx*lomz*NCi2 +  
      9*x*z*lx*lomz*NCi2 - lx*lomxmz*NCi2 -  
      x*z*lx*lomxmz*NCi2 +  
      lomz*lomxmz*NCi2 +  
      x*z*lomz*lomxmz*NCi2 - 4*x*z*lz*NCi2 -  
      lomx*lz*NCi2 - 5*x*lomx*lz*NCi2 -  
      6*x*z*lomx*lz*NCi2 + 6*lx*lz*NCi2 +  
      6*x*lx*lz*NCi2 + 12*x*z*lx*lz*NCi2 -  
      5*lomz*lz*NCi2 - x*lomz*lz*NCi2 -  
      6*x*z*lomz*lz*NCi2 - lx*lmxpz*NCi2 -  
      x*lx*lmxpz*NCi2 - 2*x*z*lx*lmxpz*NCi2 +  
      lz*lmxpz*NCi2 + x*lz*lmxpz*NCi2 +  
      2*x*z*lz*lmxpz*NCi2 - (5*pi2)/6. -  
      (x*pi2)/2. - (4*x*z*pi2)/3. + (5*NCi2*pi2)/6. +  
      (x*NCi2*pi2)/2. + (4*x*z*NCi2*pi2)/3. +  
      2*omxi - 2*z*omxi + (Li2z*omxi)/2. +  
      (z*Li2z*omxi)/2. -  
      (li2spec9*omxi)/2. -  
      (z*li2spec9*omxi)/2. -  
      li2spec13*omxi -  
      z*li2spec13*omxi -  
      (li2spec10*omxi)/2. -  
      (z*li2spec10*omxi)/2. -  
      li2spec18*omxi -  
      z*li2spec18*omxi -  
      3*z*lomx*omxi + 6*z*lx*omxi +  
      6*lomx*lx*omxi +  
      6*z*lomx*lx*omxi - 3*z*lomz*omxi -  
      2*lomx*lomz*omxi -  
      2*z*lomx*lomz*omxi +  
      (15*lx*lomz*omxi)/2. +  
      (15*z*lx*lomz*omxi)/2. -  
      lx*lomxmz*omxi -  
      z*lx*lomxmz*omxi +  
      lomz*lomxmz*omxi +  
      z*lomz*lomxmz*omxi - 4*z*lz*omxi -  
      (7*lomx*lz*omxi)/2. -  
      (7*z*lomx*lz*omxi)/2. +  
      9*lx*lz*omxi + 9*z*lx*lz*omxi -  
      (11*lomz*lz*omxi)/2. -  
      (11*z*lomz*lz*omxi)/2. -  
      (3*lx*lmxpz*omxi)/2. -  
      (3*z*lx*lmxpz*omxi)/2. +  
      (3*lz*lmxpz*omxi)/2. +  
      (3*z*lz*lmxpz*omxi)/2. -  
      2*NCi2*omxi + 2*z*NCi2*omxi -  
      (Li2z*NCi2*omxi)/2. -  
      (z*Li2z*NCi2*omxi)/2. +  
      (li2spec9*NCi2*omxi)/2. +  
      (z*li2spec9*NCi2*omxi)/2. +  
      li2spec13*NCi2*omxi +  
      z*li2spec13*NCi2*omxi +  
      (li2spec10*NCi2*omxi)/2. +  
      (z*li2spec10*NCi2*omxi)/ 
       2. + li2spec18*NCi2*omxi +  
      z*li2spec18*NCi2*omxi +  
      3*z*lomx*NCi2*omxi -  
      6*z*lx*NCi2*omxi -  
      6*lomx*lx*NCi2*omxi -  
      6*z*lomx*lx*NCi2*omxi +  
      3*z*lomz*NCi2*omxi +  
      2*lomx*lomz*NCi2*omxi +  
      2*z*lomx*lomz*NCi2*omxi -  
      (15*lx*lomz*NCi2*omxi)/2. -  
      (15*z*lx*lomz*NCi2*omxi)/2. +  
      lx*lomxmz*NCi2*omxi +  
      z*lx*lomxmz*NCi2*omxi -  
      lomz*lomxmz*NCi2*omxi -  
      z*lomz*lomxmz*NCi2*omxi +  
      4*z*lz*NCi2*omxi +  
      (7*lomx*lz*NCi2*omxi)/2. +  
      (7*z*lomx*lz*NCi2*omxi)/2. -  
      9*lx*lz*NCi2*omxi -  
      9*z*lx*lz*NCi2*omxi +  
      (11*lomz*lz*NCi2*omxi)/2. +  
      (11*z*lomz*lz*NCi2*omxi)/2. +  
      (3*lx*lmxpz*NCi2*omxi)/2. +  
      (3*z*lx*lmxpz*NCi2*omxi)/2. -  
      (3*lz*lmxpz*NCi2*omxi)/2. -  
      (3*z*lz*lmxpz*NCi2*omxi)/2. +  
      (13*pi2*omxi)/12. + (13*z*pi2*omxi)/12. -  
      (13*NCi2*pi2*omxi)/12. -  
      (13*z*NCi2*pi2*omxi)/12. + 2*omzi -  
      2*x*omzi + x*Li2z*omzi -  
      (li2spec9*omzi)/2. +  
      (x*li2spec9*omzi)/2. -  
      (li2spec13*omzi)/2. -  
      (x*li2spec13*omzi)/2. -  
      x*li2spec10*omzi -  
      (li2spec18*omzi)/2. -  
      (x*li2spec18*omzi)/2. -  
      3*x*lomx*omzi + 6*x*lx*omzi +  
      (3*lomx*lx*omzi)/2. +  
      (15*x*lomx*lx*omzi)/2. -  
      3*x*lomz*omzi -  
      (lomx*lomz*omzi)/2. -  
      (5*x*lomx*lomz*omzi)/2. +  
      3*lx*lomz*omzi +  
      6*x*lx*lomz*omzi -  
      (lx*lomxmz*omzi)/2. -  
      (x*lx*lomxmz*omzi)/2. +  
      (lomz*lomxmz*omzi)/2. +  
      (x*lomz*lomxmz*omzi)/2. -  
      4*x*lz*omzi - (lomx*lz*omzi)/2. -  
      (11*x*lomx*lz*omzi)/2. +  
      3*lx*lz*omzi + 9*x*lx*lz*omzi -  
      (5*lomz*lz*omzi)/2. -  
      (7*x*lomz*lz*omzi)/2. -  
      (lx*lmxpz*omzi)/2. -  
      (3*x*lx*lmxpz*omzi)/2. +  
      (lz*lmxpz*omzi)/2. +  
      (3*x*lz*lmxpz*omzi)/2. -  
      2*NCi2*omzi + 2*x*NCi2*omzi -  
      x*Li2z*NCi2*omzi +  
      (li2spec9*NCi2*omzi)/2. -  
      (x*li2spec9*NCi2*omzi)/2. +  
      (li2spec13*NCi2*omzi)/2. +  
      (x*li2spec13*NCi2*omzi)/2. +  
      x*li2spec10*NCi2*omzi +  
      (li2spec18*NCi2*omzi)/2. +  
      (x*li2spec18*NCi2*omzi)/ 
       2. + 3*x*lomx*NCi2*omzi -  
      6*x*lx*NCi2*omzi -  
      (3*lomx*lx*NCi2*omzi)/2. -  
      (15*x*lomx*lx*NCi2*omzi)/2. +  
      3*x*lomz*NCi2*omzi +  
      (lomx*lomz*NCi2*omzi)/2. +  
      (5*x*lomx*lomz*NCi2*omzi)/2. -  
      3*lx*lomz*NCi2*omzi -  
      6*x*lx*lomz*NCi2*omzi +  
      (lx*lomxmz*NCi2*omzi)/2. +  
      (x*lx*lomxmz*NCi2*omzi)/2. -  
      (lomz*lomxmz*NCi2*omzi)/2. -  
      (x*lomz*lomxmz*NCi2*omzi)/2. +  
      4*x*lz*NCi2*omzi +  
      (lomx*lz*NCi2*omzi)/2. +  
      (11*x*lomx*lz*NCi2*omzi)/2. -  
      3*lx*lz*NCi2*omzi -  
      9*x*lx*lz*NCi2*omzi +  
      (5*lomz*lz*NCi2*omzi)/2. +  
      (7*x*lomz*lz*NCi2*omzi)/2. +  
      (lx*lmxpz*NCi2*omzi)/2. +  
      (3*x*lx*lmxpz*NCi2*omzi)/2. -  
      (lz*lmxpz*NCi2*omzi)/2. -  
      (3*x*lz*lmxpz*NCi2*omzi)/2. +  
      (5*pi2*omzi)/12. + (11*x*pi2*omzi)/12. -  
      (5*NCi2*pi2*omzi)/12. -  
      (11*x*NCi2*pi2*omzi)/12. -  
      Li2z*omxi*omzi +  
      li2spec13*omxi*omzi +  
      li2spec10*omxi*omzi +  
      li2spec18*omxi*omzi +  
      3*lomx*omxi*omzi -  
      6*lx*omxi*omzi -  
      9*lomx*lx*omxi*omzi +  
      3*lomz*omxi*omzi +  
      3*lomx*lomz*omxi*omzi -  
      9*lx*lomz*omxi*omzi +  
      lx*lomxmz*omxi*omzi -  
      lomz*lomxmz*omxi*omzi +  
      4*lz*omxi*omzi +  
      6*lomx*lz*omxi*omzi -  
      12*lx*lz*omxi*omzi +  
      6*lomz*lz*omxi*omzi +  
      2*lx*lmxpz*omxi*omzi -  
      2*lz*lmxpz*omxi*omzi +  
      Li2z*NCi2*omxi*omzi -  
      li2spec13*NCi2*omxi*omzi -  
      li2spec10*NCi2*omxi* 
       omzi - li2spec18*NCi2* 
       omxi*omzi -  
      3*lomx*NCi2*omxi*omzi +  
      6*lx*NCi2*omxi*omzi +  
      9*lomx*lx*NCi2*omxi*omzi -  
      3*lomz*NCi2*omxi*omzi -  
      3*lomx*lomz*NCi2*omxi*omzi +  
      9*lx*lomz*NCi2*omxi*omzi -  
      lx*lomxmz*NCi2*omxi*omzi +  
      lomz*lomxmz*NCi2*omxi*omzi -  
      4*lz*NCi2*omxi*omzi -  
      6*lomx*lz*NCi2*omxi*omzi +  
      12*lx*lz*NCi2*omxi*omzi -  
      6*lomz*lz*NCi2*omxi*omzi -  
      2*lx*lmxpz*NCi2*omxi*omzi +  
      2*lz*lmxpz*NCi2*omxi*omzi -  
      (4*pi2*omxi*omzi)/3. +  
      (4*NCi2*pi2*omxi*omzi)/3. +  
      (3*lomx2)/2. + 2*x*lomx2 +  
      (7*x*z*lomx2)/2. - (3*NCi2*lomx2)/2. -  
      2*x*NCi2*lomx2 -  
      (7*x*z*NCi2*lomx2)/2. -  
      (5*omxi*lomx2)/2. -  
      (5*z*omxi*lomx2)/2. +  
      (5*NCi2*omxi*lomx2)/2. +  
      (5*z*NCi2*omxi*lomx2)/2. -  
      (3*omzi*lomx2)/4. -  
      (11*x*omzi*lomx2)/4. +  
      (3*NCi2*omzi*lomx2)/4. +  
      (11*x*NCi2*omzi*lomx2)/4. +  
      (7*omxi*omzi*lomx2)/2. -  
      (7*NCi2*omxi*omzi*lomx2)/2. +  
      (7*lx2)/2. + (7*x*lx2)/2. + 7*x*z*lx2 -  
      (7*NCi2*lx2)/2. - (7*x*NCi2*lx2)/2. -  
      7*x*z*NCi2*lx2 - (21*omxi*lx2)/4. -  
      (21*z*omxi*lx2)/4. +  
      (21*NCi2*omxi*lx2)/4. +  
      (21*z*NCi2*omxi*lx2)/4. -  
      (7*omzi*lx2)/4. -  
      (21*x*omzi*lx2)/4. +  
      (7*NCi2*omzi*lx2)/4. +  
      (21*x*NCi2*omzi*lx2)/4. +  
      7*omxi*omzi*lx2 -  
      7*NCi2*omxi*omzi*lx2 +  
      (5*lomz2)/2. + (x*lomz2)/2. +  
      3*x*z*lomz2 - (5*NCi2*lomz2)/2. -  
      (x*NCi2*lomz2)/2. -  
      3*x*z*NCi2*lomz2 -  
      (11*omxi*lomz2)/4. -  
      (11*z*omxi*lomz2)/4. +  
      (11*NCi2*omxi*lomz2)/4. +  
      (11*z*NCi2*omxi*lomz2)/4. -  
      (5*omzi*lomz2)/4. -  
      (7*x*omzi*lomz2)/4. +  
      (5*NCi2*omzi*lomz2)/4. +  
      (7*x*NCi2*omzi*lomz2)/4. +  
      3*omxi*omzi*lomz2 -  
      3*NCi2*omxi*omzi*lomz2 +  
      (5*lz2)/2. + (5*x*lz2)/2. + 5*x*z*lz2 -  
      (5*NCi2*lz2)/2. - (5*x*NCi2*lz2)/2. -  
      5*x*z*NCi2*lz2 - (15*omxi*lz2)/4. -  
      (15*z*omxi*lz2)/4. +  
      (15*NCi2*omxi*lz2)/4. +  
      (15*z*NCi2*omxi*lz2)/4. -  
      (5*omzi*lz2)/4. -  
      (15*x*omzi*lz2)/4. +  
      (5*NCi2*omzi*lz2)/4. +  
      (15*x*NCi2*omzi*lz2)/4. +  
      5*omxi*omzi*lz2 -  
      5*NCi2*omxi*omzi*lz2)*Tu3 +  
   (-3 - x*Li2z - x*z*Li2z + li2spec9 -  
      x*li2spec9 - li2spec14 -  
      x*z*li2spec14 - li2spec16 -  
      x*z*li2spec16 + x*li2spec10 +  
      x*z*li2spec10 + 3*x*z*lomx -  
      6*x*z*lx - 4*lomx*lx - 6*x*lomx*lx -  
      10*x*z*lomx*lx + 3*x*z*lomz + 2*lomx*lomz +  
      2*x*lomx*lomz + 4*x*z*lomx*lomz -  
      5*lx*lomz - 3*x*lx*lomz - 8*x*z*lx*lomz +  
      4*x*z*lz + 3*lomx*lz + 5*x*lomx*lz +  
      8*x*z*lomx*lz - 5*lx*lz - 6*x*lx*lz -  
      11*x*z*lx*lz + 4*lomz*lz + x*lomz*lz +  
      5*x*z*lomz*lz + lx*lmxpz + x*lx*lmxpz +  
      2*x*z*lx*lmxpz - lz*lmxpz - x*lz*lmxpz -  
      2*x*z*lz*lmxpz + lx*lmopxpz +  
      x*z*lx*lmopxpz - lomz*lmopxpz -  
      x*z*lomz*lmopxpz + 3*NCi2 + x*Li2z*NCi2 +  
      x*z*Li2z*NCi2 - li2spec9*NCi2 +  
      x*li2spec9*NCi2 +  
      li2spec14*NCi2 +  
      x*z*li2spec14*NCi2 +  
      li2spec16*NCi2 +  
      x*z*li2spec16*NCi2 -  
      x*li2spec10*NCi2 -  
      x*z*li2spec10*NCi2 -  
      3*x*z*lomx*NCi2 + 6*x*z*lx*NCi2 +  
      4*lomx*lx*NCi2 + 6*x*lomx*lx*NCi2 +  
      10*x*z*lomx*lx*NCi2 - 3*x*z*lomz*NCi2 -  
      2*lomx*lomz*NCi2 -  
      2*x*lomx*lomz*NCi2 -  
      4*x*z*lomx*lomz*NCi2 +  
      5*lx*lomz*NCi2 + 3*x*lx*lomz*NCi2 +  
      8*x*z*lx*lomz*NCi2 - 4*x*z*lz*NCi2 -  
      3*lomx*lz*NCi2 - 5*x*lomx*lz*NCi2 -  
      8*x*z*lomx*lz*NCi2 + 5*lx*lz*NCi2 +  
      6*x*lx*lz*NCi2 + 11*x*z*lx*lz*NCi2 -  
      4*lomz*lz*NCi2 - x*lomz*lz*NCi2 -  
      5*x*z*lomz*lz*NCi2 - lx*lmxpz*NCi2 -  
      x*lx*lmxpz*NCi2 - 2*x*z*lx*lmxpz*NCi2 +  
      lz*lmxpz*NCi2 + x*lz*lmxpz*NCi2 +  
      2*x*z*lz*lmxpz*NCi2 -  
      lx*lmopxpz*NCi2 -  
      x*z*lx*lmopxpz*NCi2 +  
      lomz*lmopxpz*NCi2 +  
      x*z*lomz*lmopxpz*NCi2 - pi2/6. -  
      (x*pi2)/2. - (2*x*z*pi2)/3. + (NCi2*pi2)/6. +  
      (x*NCi2*pi2)/2. + (2*x*z*NCi2*pi2)/3. +  
      2*omxi - 2*z*omxi + (Li2z*omxi)/2. +  
      (z*Li2z*omxi)/2. -  
      (li2spec9*omxi)/2. -  
      (z*li2spec9*omxi)/2. +  
      li2spec14*omxi +  
      z*li2spec14*omxi +  
      li2spec16*omxi +  
      z*li2spec16*omxi -  
      (li2spec10*omxi)/2. -  
      (z*li2spec10*omxi)/2. -  
      3*z*lomx*omxi + 6*z*lx*omxi +  
      7*lomx*lx*omxi +  
      7*z*lomx*lx*omxi - 3*z*lomz*omxi -  
      3*lomx*lomz*omxi -  
      3*z*lomx*lomz*omxi +  
      (13*lx*lomz*omxi)/2. +  
      (13*z*lx*lomz*omxi)/2. - 4*z*lz*omxi -  
      (11*lomx*lz*omxi)/2. -  
      (11*z*lomx*lz*omxi)/2. +  
      8*lx*lz*omxi + 8*z*lx*lz*omxi -  
      (9*lomz*lz*omxi)/2. -  
      (9*z*lomz*lz*omxi)/2. -  
      (3*lx*lmxpz*omxi)/2. -  
      (3*z*lx*lmxpz*omxi)/2. +  
      (3*lz*lmxpz*omxi)/2. +  
      (3*z*lz*lmxpz*omxi)/2. -  
      lx*lmopxpz*omxi -  
      z*lx*lmopxpz*omxi +  
      lomz*lmopxpz*omxi +  
      z*lomz*lmopxpz*omxi -  
      2*NCi2*omxi + 2*z*NCi2*omxi -  
      (Li2z*NCi2*omxi)/2. -  
      (z*Li2z*NCi2*omxi)/2. +  
      (li2spec9*NCi2*omxi)/2. +  
      (z*li2spec9*NCi2*omxi)/2. -  
      li2spec14*NCi2*omxi -  
      z*li2spec14*NCi2*omxi -  
      li2spec16*NCi2*omxi -  
      z*li2spec16*NCi2*omxi +  
      (li2spec10*NCi2*omxi)/2. +  
      (z*li2spec10*NCi2*omxi)/ 
       2. + 3*z*lomx*NCi2*omxi -  
      6*z*lx*NCi2*omxi -  
      7*lomx*lx*NCi2*omxi -  
      7*z*lomx*lx*NCi2*omxi +  
      3*z*lomz*NCi2*omxi +  
      3*lomx*lomz*NCi2*omxi +  
      3*z*lomx*lomz*NCi2*omxi -  
      (13*lx*lomz*NCi2*omxi)/2. -  
      (13*z*lx*lomz*NCi2*omxi)/2. +  
      4*z*lz*NCi2*omxi +  
      (11*lomx*lz*NCi2*omxi)/2. +  
      (11*z*lomx*lz*NCi2*omxi)/2. -  
      8*lx*lz*NCi2*omxi -  
      8*z*lx*lz*NCi2*omxi +  
      (9*lomz*lz*NCi2*omxi)/2. +  
      (9*z*lomz*lz*NCi2*omxi)/2. +  
      (3*lx*lmxpz*NCi2*omxi)/2. +  
      (3*z*lx*lmxpz*NCi2*omxi)/2. -  
      (3*lz*lmxpz*NCi2*omxi)/2. -  
      (3*z*lz*lmxpz*NCi2*omxi)/2. +  
      lx*lmopxpz*NCi2*omxi +  
      z*lx*lmopxpz*NCi2*omxi -  
      lomz*lmopxpz*NCi2*omxi -  
      z*lomz*lmopxpz*NCi2*omxi +  
      (5*pi2*omxi)/12. + (5*z*pi2*omxi)/12. -  
      (5*NCi2*pi2*omxi)/12. -  
      (5*z*NCi2*pi2*omxi)/12. + 2*omzi -  
      2*x*omzi + x*Li2z*omzi -  
      (li2spec9*omzi)/2. +  
      (x*li2spec9*omzi)/2. +  
      (li2spec14*omzi)/2. +  
      (x*li2spec14*omzi)/2. +  
      (li2spec16*omzi)/2. +  
      (x*li2spec16*omzi)/2. -  
      x*li2spec10*omzi -  
      3*x*lomx*omzi + 6*x*lx*omzi +  
      2*lomx*lx*omzi +  
      8*x*lomx*lx*omzi - 3*x*lomz*omzi -  
      lomx*lomz*omzi -  
      3*x*lomx*lomz*omzi +  
      (5*lx*lomz*omzi)/2. +  
      (11*x*lx*lomz*omzi)/2. - 4*x*lz*omzi -  
      (3*lomx*lz*omzi)/2. -  
      (13*x*lomx*lz*omzi)/2. +  
      (5*lx*lz*omzi)/2. +  
      (17*x*lx*lz*omzi)/2. -  
      2*lomz*lz*omzi -  
      3*x*lomz*lz*omzi -  
      (lx*lmxpz*omzi)/2. -  
      (3*x*lx*lmxpz*omzi)/2. +  
      (lz*lmxpz*omzi)/2. +  
      (3*x*lz*lmxpz*omzi)/2. -  
      (lx*lmopxpz*omzi)/2. -  
      (x*lx*lmopxpz*omzi)/2. +  
      (lomz*lmopxpz*omzi)/2. +  
      (x*lomz*lmopxpz*omzi)/2. -  
      2*NCi2*omzi + 2*x*NCi2*omzi -  
      x*Li2z*NCi2*omzi +  
      (li2spec9*NCi2*omzi)/2. -  
      (x*li2spec9*NCi2*omzi)/2. -  
      (li2spec14*NCi2*omzi)/2. -  
      (x*li2spec14*NCi2*omzi)/ 
       2. - (li2spec16*NCi2*omzi)/2. -  
      (x*li2spec16*NCi2*omzi)/2. +  
      x*li2spec10*NCi2*omzi +  
      3*x*lomx*NCi2*omzi -  
      6*x*lx*NCi2*omzi -  
      2*lomx*lx*NCi2*omzi -  
      8*x*lomx*lx*NCi2*omzi +  
      3*x*lomz*NCi2*omzi +  
      lomx*lomz*NCi2*omzi +  
      3*x*lomx*lomz*NCi2*omzi -  
      (5*lx*lomz*NCi2*omzi)/2. -  
      (11*x*lx*lomz*NCi2*omzi)/2. +  
      4*x*lz*NCi2*omzi +  
      (3*lomx*lz*NCi2*omzi)/2. +  
      (13*x*lomx*lz*NCi2*omzi)/2. -  
      (5*lx*lz*NCi2*omzi)/2. -  
      (17*x*lx*lz*NCi2*omzi)/2. +  
      2*lomz*lz*NCi2*omzi +  
      3*x*lomz*lz*NCi2*omzi +  
      (lx*lmxpz*NCi2*omzi)/2. +  
      (3*x*lx*lmxpz*NCi2*omzi)/2. -  
      (lz*lmxpz*NCi2*omzi)/2. -  
      (3*x*lz*lmxpz*NCi2*omzi)/2. +  
      (lx*lmopxpz*NCi2*omzi)/2. +  
      (x*lx*lmopxpz*NCi2*omzi)/2. -  
      (lomz*lmopxpz*NCi2*omzi)/2. -  
      (x*lomz*lmopxpz*NCi2*omzi)/2. +  
      (pi2*omzi)/12. + (7*x*pi2*omzi)/12. -  
      (NCi2*pi2*omzi)/12. -  
      (7*x*NCi2*pi2*omzi)/12. -  
      Li2z*omxi*omzi -  
      li2spec14*omxi*omzi -  
      li2spec16*omxi*omzi +  
      li2spec10*omxi*omzi +  
      3*lomx*omxi*omzi -  
      6*lx*omxi*omzi -  
      10*lomx*lx*omxi*omzi +  
      3*lomz*omxi*omzi +  
      4*lomx*lomz*omxi*omzi -  
      8*lx*lomz*omxi*omzi +  
      4*lz*omxi*omzi +  
      8*lomx*lz*omxi*omzi -  
      11*lx*lz*omxi*omzi +  
      5*lomz*lz*omxi*omzi +  
      2*lx*lmxpz*omxi*omzi -  
      2*lz*lmxpz*omxi*omzi +  
      lx*lmopxpz*omxi*omzi -  
      lomz*lmopxpz*omxi*omzi +  
      Li2z*NCi2*omxi*omzi +  
      li2spec14*NCi2*omxi* 
       omzi + li2spec16*NCi2*omxi* 
       omzi - li2spec10*NCi2* 
       omxi*omzi -  
      3*lomx*NCi2*omxi*omzi +  
      6*lx*NCi2*omxi*omzi +  
      10*lomx*lx*NCi2*omxi*omzi -  
      3*lomz*NCi2*omxi*omzi -  
      4*lomx*lomz*NCi2*omxi*omzi +  
      8*lx*lomz*NCi2*omxi*omzi -  
      4*lz*NCi2*omxi*omzi -  
      8*lomx*lz*NCi2*omxi*omzi +  
      11*lx*lz*NCi2*omxi*omzi -  
      5*lomz*lz*NCi2*omxi*omzi -  
      2*lx*lmxpz*NCi2*omxi*omzi +  
      2*lz*lmxpz*NCi2*omxi*omzi -  
      lx*lmopxpz*NCi2*omxi*omzi +  
      lomz*lmopxpz*NCi2*omxi*omzi -  
      (2*pi2*omxi*omzi)/3. +  
      (2*NCi2*pi2*omxi*omzi)/3. +  
      lomx2/2. + 2*x*lomx2 +  
      (5*x*z*lomx2)/2. - (NCi2*lomx2)/2. -  
      2*x*NCi2*lomx2 -  
      (5*x*z*NCi2*lomx2)/2. -  
      (3*omxi*lomx2)/2. -  
      (3*z*omxi*lomx2)/2. +  
      (3*NCi2*omxi*lomx2)/2. +  
      (3*z*NCi2*omxi*lomx2)/2. -  
      (omzi*lomx2)/4. -  
      (9*x*omzi*lomx2)/4. +  
      (NCi2*omzi*lomx2)/4. +  
      (9*x*NCi2*omzi*lomx2)/4. +  
      (5*omxi*omzi*lomx2)/2. -  
      (5*NCi2*omxi*omzi*lomx2)/2. +  
      3*lx2 + (7*x*lx2)/2. + (13*x*z*lx2)/2. -  
      3*NCi2*lx2 - (7*x*NCi2*lx2)/2. -  
      (13*x*z*NCi2*lx2)/2. -  
      (19*omxi*lx2)/4. -  
      (19*z*omxi*lx2)/4. +  
      (19*NCi2*omxi*lx2)/4. +  
      (19*z*NCi2*omxi*lx2)/4. -  
      (3*omzi*lx2)/2. - 5*x*omzi*lx2 +  
      (3*NCi2*omzi*lx2)/2. +  
      5*x*NCi2*omzi*lx2 +  
      (13*omxi*omzi*lx2)/2. -  
      (13*NCi2*omxi*omzi*lx2)/2. +  
      2*lomz2 + (x*lomz2)/2. +  
      (5*x*z*lomz2)/2. - 2*NCi2*lomz2 -  
      (x*NCi2*lomz2)/2. -  
      (5*x*z*NCi2*lomz2)/2. -  
      (9*omxi*lomz2)/4. -  
      (9*z*omxi*lomz2)/4. +  
      (9*NCi2*omxi*lomz2)/4. +  
      (9*z*NCi2*omxi*lomz2)/4. -  
      omzi*lomz2 -  
      (3*x*omzi*lomz2)/2. +  
      NCi2*omzi*lomz2 +  
      (3*x*NCi2*omzi*lomz2)/2. +  
      (5*omxi*omzi*lomz2)/2. -  
      (5*NCi2*omxi*omzi*lomz2)/2. +  
      (3*lz2)/2. + (5*x*lz2)/2. + 4*x*z*lz2 -  
      (3*NCi2*lz2)/2. - (5*x*NCi2*lz2)/2. -  
      4*x*z*NCi2*lz2 - (11*omxi*lz2)/4. -  
      (11*z*omxi*lz2)/4. +  
      (11*NCi2*omxi*lz2)/4. +  
      (11*z*NCi2*omxi*lz2)/4. -  
      (3*omzi*lz2)/4. -  
      (13*x*omzi*lz2)/4. +  
      (3*NCi2*omzi*lz2)/4. +  
      (13*x*NCi2*omzi*lz2)/4. +  
      4*omxi*omzi*lz2 -  
      4*NCi2*omxi*omzi*lz2)*Tu4;
};
double C2TQ2QPS_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double sqrtxz3  = sqrt(x/z);
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double x2 = x * x;
  const double xi  = 1. / x;
  const double xi2 = xi * xi;
  const double z2 = z * z;
  const double zi  = 1. / z;
  const double zi2 = zi * zi;
  const double lx  = log(x);
  const double lz  = log(z);
  const double lspec3   = log(sqrtxz3);
  const double lspec4   = log(sqrtxz3*z);
  const double atanspec1 = atan(sqrtxz3);
  const double atanspec2 = atan(sqrtxz3*z);
  const double itani1 = InvTanInt(-sqrtxz3);
  const double itani2 = InvTanInt(sqrtxz3);
  const double itani3 = InvTanInt(sqrtxz3*z);
  return (15*NC)/8. - (15*NC*x)/8. + (5*NC*z)/16. - (5*NC*x*z)/16. -  
   (13*NC*sqrtxz3*itani1)/8. -  
   (9*NC*sqrtxz3*x*z*itani1)/16. +  
   (13*NC*sqrtxz3*itani2)/8. +  
   (9*NC*sqrtxz3*x*z*itani2)/16. -  
   (13*NC*sqrtxz3*itani3)/4. -  
   (9*NC*sqrtxz3*x*z*itani3)/8. -  
   (13*NC*sqrtxz3*atanspec1*lspec3)/4. -  
   (9*NC*sqrtxz3*x*z*atanspec1*lspec3)/8. + (17*NC*lx)/16. +  
   (17*NC*x*lx)/16. + (5*NC*z*lx)/32. + (5*NC*x*z*lx)/32. -  
   (17*NC*lz)/16. + (17*NC*x*lz)/16. + (5*NC*z*lz)/32. -  
   (5*NC*x*z*lz)/32. - NC*lx*lz - NC*x*lx*lz +  
   (13*NC*sqrtxz3*atanspec2*lspec4)/4. +  
   (9*NC*sqrtxz3*x*z*atanspec2*lspec4)/8. -  
   (15*NCi)/8. + (15*x*NCi)/8. - (5*z*NCi)/16. +  
   (5*x*z*NCi)/16. + (13*sqrtxz3*itani1*NCi)/8. +  
   (9*sqrtxz3*x*z*itani1*NCi)/16. -  
   (13*sqrtxz3*itani2*NCi)/8. -  
   (9*sqrtxz3*x*z*itani2*NCi)/16. +  
   (13*sqrtxz3*itani3*NCi)/4. +  
   (9*sqrtxz3*x*z*itani3*NCi)/8. +  
   (13*sqrtxz3*atanspec1*lspec3*NCi)/4. +  
   (9*sqrtxz3*x*z*atanspec1*lspec3*NCi)/8. -  
   (17*lx*NCi)/16. - (17*x*lx*NCi)/16. -  
   (5*z*lx*NCi)/32. - (5*x*z*lx*NCi)/32. +  
   (17*lz*NCi)/16. - (17*x*lz*NCi)/16. -  
   (5*z*lz*NCi)/32. + (5*x*z*lz*NCi)/32. +  
   lx*lz*NCi + x*lx*lz*NCi -  
   (13*sqrtxz3*atanspec2*lspec4*NCi)/4. -  
   (9*sqrtxz3*x*z*atanspec2*lspec4*NCi)/8. -  
   (5*NC*sqrtxz3*itani1*xi2)/32. +  
   (5*NC*sqrtxz3*itani2*xi2)/32. -  
   (5*NC*sqrtxz3*itani3*xi2)/16. -  
   (5*NC*sqrtxz3*atanspec1*lspec3*xi2)/16. +  
   (5*NC*sqrtxz3*atanspec2*lspec4*xi2)/16. +  
   (5*sqrtxz3*itani1*NCi*xi2)/32. -  
   (5*sqrtxz3*itani2*NCi*xi2)/32. +  
   (5*sqrtxz3*itani3*NCi*xi2)/16. +  
   (5*sqrtxz3*atanspec1*lspec3*NCi*xi2)/16. -  
   (5*sqrtxz3*atanspec2*lspec4*NCi*xi2)/16. +  
   (5*NC*xi)/16. - (9*NC*sqrtxz3*z*itani1*xi)/ 
    16. + (9*NC*sqrtxz3*z*itani2*xi)/16. -  
   (9*NC*sqrtxz3*z*itani3*xi)/8. -  
   (9*NC*sqrtxz3*z*atanspec1*lspec3*xi)/8. -  
   (5*NC*lx*xi)/32. - (5*NC*lz*xi)/32. +  
   (9*NC*sqrtxz3*z*atanspec2*lspec4*xi)/8. -  
   (5*NCi*xi)/16. +  
   (9*sqrtxz3*z*itani1*NCi*xi)/16. -  
   (9*sqrtxz3*z*itani2*NCi*xi)/16. +  
   (9*sqrtxz3*z*itani3*NCi*xi)/8. +  
   (9*sqrtxz3*z*atanspec1*lspec3*NCi*xi)/8. +  
   (5*lx*NCi*xi)/32. +  
   (5*lz*NCi*xi)/32. -  
   (9*sqrtxz3*z*atanspec2*lspec4*NCi*xi)/8. -  
   (5*NC*x2)/16. - (5*NC*sqrtxz3*itani1*x2)/32. +  
   (5*NC*sqrtxz3*itani2*x2)/32. -  
   (5*NC*sqrtxz3*itani3*x2)/16. -  
   (5*NC*sqrtxz3*atanspec1*lspec3*x2)/16. -  
   (5*NC*lx*x2)/32. + (5*NC*lz*x2)/32. +  
   (5*NC*sqrtxz3*atanspec2*lspec4*x2)/16. +  
   (5*NCi*x2)/16. +  
   (5*sqrtxz3*itani1*NCi*x2)/32. -  
   (5*sqrtxz3*itani2*NCi*x2)/32. +  
   (5*sqrtxz3*itani3*NCi*x2)/16. +  
   (5*sqrtxz3*atanspec1*lspec3*NCi*x2)/16. +  
   (5*lx*NCi*x2)/32. - (5*lz*NCi*x2)/32. -  
   (5*sqrtxz3*atanspec2*lspec4*NCi*x2)/16. -  
   (5*NC*zi2)/16. + (5*NC*x*zi2)/16. -  
   (5*NC*sqrtxz3*itani1*zi2)/32. +  
   (5*NC*sqrtxz3*itani2*zi2)/32. -  
   (5*NC*sqrtxz3*itani3*zi2)/16. -  
   (5*NC*sqrtxz3*atanspec1*lspec3*zi2)/16. -  
   (5*NC*lx*zi2)/32. - (5*NC*x*lx*zi2)/32. +  
   (5*NC*lz*zi2)/32. - (5*NC*x*lz*zi2)/32. +  
   (5*NC*sqrtxz3*atanspec2*lspec4*zi2)/16. +  
   (5*NCi*zi2)/16. - (5*x*NCi*zi2)/16. +  
   (5*sqrtxz3*itani1*NCi*zi2)/32. -  
   (5*sqrtxz3*itani2*NCi*zi2)/32. +  
   (5*sqrtxz3*itani3*NCi*zi2)/16. +  
   (5*sqrtxz3*atanspec1*lspec3*NCi*zi2)/16. +  
   (5*lx*NCi*zi2)/32. +  
   (5*x*lx*NCi*zi2)/32. -  
   (5*lz*NCi*zi2)/32. +  
   (5*x*lz*NCi*zi2)/32. -  
   (5*sqrtxz3*atanspec2*lspec4*NCi*zi2)/16. -  
   (15*NC*zi)/8. + (15*NC*x*zi)/8. -  
   (9*NC*sqrtxz3*x*itani1*zi)/16. +  
   (9*NC*sqrtxz3*x*itani2*zi)/16. -  
   (9*NC*sqrtxz3*x*itani3*zi)/8. -  
   (9*NC*sqrtxz3*x*atanspec1*lspec3*zi)/8. -  
   (17*NC*lx*zi)/16. - (17*NC*x*lx*zi)/16. -  
   (17*NC*lz*zi)/16. + (17*NC*x*lz*zi)/16. -  
   NC*lx*lz*zi - NC*x*lx*lz*zi +  
   (9*NC*sqrtxz3*x*atanspec2*lspec4*zi)/8. +  
   (15*NCi*zi)/8. - (15*x*NCi*zi)/8. +  
   (9*sqrtxz3*x*itani1*NCi*zi)/16. -  
   (9*sqrtxz3*x*itani2*NCi*zi)/16. +  
   (9*sqrtxz3*x*itani3*NCi*zi)/8. +  
   (9*sqrtxz3*x*atanspec1*lspec3*NCi*zi)/8. +  
   (17*lx*NCi*zi)/16. +  
   (17*x*lx*NCi*zi)/16. +  
   (17*lz*NCi*zi)/16. -  
   (17*x*lz*NCi*zi)/16. +  
   lx*lz*NCi*zi +  
   x*lx*lz*NCi*zi -  
   (9*sqrtxz3*x*atanspec2*lspec4*NCi*zi)/8. -  
   (5*NC*xi*zi)/16. -  
   (9*NC*sqrtxz3*itani1*xi*zi)/16. +  
   (9*NC*sqrtxz3*itani2*xi*zi)/16. -  
   (9*NC*sqrtxz3*itani3*xi*zi)/8. -  
   (9*NC*sqrtxz3*atanspec1*lspec3*xi*zi)/8. +  
   (5*NC*lx*xi*zi)/32. -  
   (5*NC*lz*xi*zi)/32. +  
   (9*NC*sqrtxz3*atanspec2*lspec4*xi*zi)/8. +  
   (5*NCi*xi*zi)/16. +  
   (9*sqrtxz3*itani1*NCi*xi*zi)/16. -  
   (9*sqrtxz3*itani2*NCi*xi*zi)/16. +  
   (9*sqrtxz3*itani3*NCi*xi*zi)/8. +  
   (9*sqrtxz3*atanspec1*lspec3*NCi*xi*zi)/ 
    8. - (5*lx*NCi*xi*zi)/32. +  
   (5*lz*NCi*xi*zi)/32. -  
   (9*sqrtxz3*atanspec2*lspec4*NCi*xi* 
      zi)/8. + (5*NC*x2*zi)/16. +  
   (5*NC*lx*x2*zi)/32. +  
   (5*NC*lz*x2*zi)/32. -  
   (5*NCi*x2*zi)/16. -  
   (5*lx*NCi*x2*zi)/32. -  
   (5*lz*NCi*x2*zi)/32. -  
   (5*NC*sqrtxz3*itani1*z2)/32. +  
   (5*NC*sqrtxz3*itani2*z2)/32. -  
   (5*NC*sqrtxz3*itani3*z2)/16. -  
   (5*NC*sqrtxz3*atanspec1*lspec3*z2)/16. +  
   (5*NC*sqrtxz3*atanspec2*lspec4*z2)/16. +  
   (5*sqrtxz3*itani1*NCi*z2)/32. -  
   (5*sqrtxz3*itani2*NCi*z2)/32. +  
   (5*sqrtxz3*itani3*NCi*z2)/16. +  
   (5*sqrtxz3*atanspec1*lspec3*NCi*z2)/16. -  
   (5*sqrtxz3*atanspec2*lspec4*NCi*z2)/16.;
};
double C2TQ2QB_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double sqrtxz1  = sqrt(1 - 2*z + z*z + 4*x*z);
  const double sqrtxz1i = 1. / sqrtxz1;
  const double poly2    = 1 + 2*x + x*x - 4*x*z;
  const double poly2i   = 1. / poly2;
  const double sqrtxz2  = sqrt(poly2);
  const double sqrtxz2i = 1. / sqrtxz2;
  const double sqrtxz3  = sqrt(x/z);
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double NCi2 = NCi * NCi;
  const double l2  = log(2);
  const double l22 = l2 * l2;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double x2 = x * x;
  const double x3 = x * x2;
  const double x4 = x * x3;
  const double x5 = x * x4;
  const double xi  = 1. / x;
  const double xi2 = xi * xi;
  const double z2 = z * z;
  const double zi  = 1. / z;
  const double omxi = 1. / ( 1 - x );
  const double opxi = 1. / ( 1 + x );
  const double omzi  = 1. / ( 1 - z );
  const double opzi  = 1. / ( 1 + z );
  const double mopzi = 1. / ( - 1 + z );
  const double lx  = log(x);
  const double lx2 = lx * lx;
  const double lz  = log(z);
  const double lz2 = lz * lz;
  const double lomx  = log(1 - x);
  const double lomz  = log(1 - z);
  const double lomz2 = lomz * lomz;
  const double lopx  = log(1 + x);
  const double lopz  = log(1 + z);
  const double Li2x  = apfel::dilog(x);
  const double Li2mx = apfel::dilog(-x);
  const double Li2z  = apfel::dilog(z);
  const double Li2mz = apfel::dilog(-z);
  const double xmzi  = 1. / ( x - z );
  const double spec1i = 1. / ( 1 + sqrtxz1 - z );
  const double lxpz    = log(x + z);
  const double lopxz   = log(1 + x*z);
  const double lopxzi  = log(1 + x*zi);
  const double li2omxzi = apfel::dilog(1 - x*zi);
  const double lspec1   = log(1 + sqrtxz1 - z);
  const double lspec1_2 = lspec1 * lspec1;
  const double lspec2   = log(1 + sqrtxz1 + z);
  const double lspec3   = log(sqrtxz3);
  const double lspec4   = log(sqrtxz3*z);
  const double lspec5   = log(1 - sqrtxz2 + x);
  const double lspec6   = log(1 + sqrtxz2 + x);
  const double lspec7   = log(1 - 2*z + 4*x*z + z2);
  const double lspec7_2 = lspec7 * lspec7;
  const double li2spec1  = apfel::dilog(0.5 - sqrtxz1/2. - z/2.);
  const double li2spec2  = apfel::dilog(0.5 - sqrtxz1/2. + z/2.);
  const double li2spec3  = apfel::dilog(0.5 - zi/2. - (sqrtxz1*zi)/2.);
  const double li2spec4  = apfel::dilog(0.5 + zi/2. - (sqrtxz1*zi)/2.);
  const double li2spec5  = apfel::dilog(0.5 - sqrtxz2/2. - x/2.);
  const double li2spec6  = apfel::dilog(0.5 + sqrtxz2/2. - x/2.);
  const double li2spec7  = apfel::dilog(0.5 - xi/2. - (sqrtxz2*xi)/2.);
  const double li2spec8  = apfel::dilog(0.5 - xi/2. + (sqrtxz2*xi)/2.);
  const double li2spec19 = apfel::dilog(-(x*z));
  const double li2spec20 = apfel::dilog(-(x*zi));
  const double li2spec21 = apfel::dilog((1 - z)*sqrtxz1i);
  const double li2spec22 = apfel::dilog(-spec1i + sqrtxz1*spec1i + z*spec1i);
  const double li2spec23 = apfel::dilog(sqrtxz1*mopzi);
  const double atanspec1 = atan(sqrtxz3);
  const double atanspec2 = atan(sqrtxz3*z);
  const double itani1 = InvTanInt(-sqrtxz3);
  const double itani2 = InvTanInt(sqrtxz3);
  const double itani3 = InvTanInt(sqrtxz3*z);
  return -4 - (13*NC)/3. + (7*x)/2. + 2*NC*x + 5*z + (157*NC*z)/24. -  
   (19*x*z)/4. - (119*NC*x*z)/24. - (3*NC*sqrtxz3*itani1)/4. +  
   (3*NC*sqrtxz3*x*z*itani1)/4. +  
   (3*NC*sqrtxz3*itani2)/4. -  
   (3*NC*sqrtxz3*x*z*itani2)/4. -  
   (3*NC*sqrtxz3*itani3)/2. +  
   (3*NC*sqrtxz3*x*z*itani3)/2. + 2*Li2mx + 2*NC*Li2mx +  
   2*x*Li2mx + 2*NC*x*Li2mx - 4*NC*z*Li2mx - 4*NC*x*z*Li2mx -  
   Li2x + NC*Li2x - x*Li2x + 2*NC*x*Li2x + x*z*Li2x -  
   2*NC*x*z*Li2x + 2*x*li2spec1 -  
   NC*x*li2spec1 - 2*NC*z*li2spec1 -  
   2*x*z*li2spec1 - 2*x*li2spec2 +  
   NC*x*li2spec2 + 2*NC*z*li2spec2 +  
   2*x*z*li2spec2 + (NC*Li2z)/2. + (3*NC*x*Li2z)/2. -  
   NC*x*z*Li2z - 2*x*li2spec19 + NC*x*li2spec19 + 2*NC*z*li2spec19 -  
   NC*li2spec20 - 2*x*z*li2spec20 -  
   2*NC*x*z*li2spec20 +  
   NC*li2spec3 -  
   2*x*li2spec3 +  
   2*x*z*li2spec3 +  
   2*NC*x*z*li2spec3 -  
   NC*li2spec4 +  
   2*x*li2spec4 -  
   2*x*z*li2spec4 -  
   2*NC*x*z*li2spec4 -  
   NC*x*li2omxzi + 2*NC*x*z*li2omxzi +  
   2*NC*sqrtxz1*x*l2 - (3*NC*sqrtxz3*atanspec1*lspec3)/2. +  
   (3*NC*sqrtxz3*x*z*atanspec1*lspec3)/2. -  
   (3*NC*lomx)/4. - (NC*x*lomx)/4. + (3*NC*z*lomx)/4. +  
   (NC*x*z*lomx)/4. - lx - (7*NC*lx)/4. - x*lx -  
   (17*NC*x*lx)/4. + NC*sqrtxz1*x*lx + 2*z*lx +  
   (5*NC*z*lx)/4. + 2*x*z*lx + (13*NC*x*z*lx)/4. -  
   NC*l2*lx - 2*x*l2*lx + 2*NC*x*l2*lx +  
   4*NC*z*l2*lx + 2*x*z*l2*lx - 2*NC*x*z*l2*lx -  
   lomx*lx - x*lomx*lx + NC*x*lomx*lx +  
   x*z*lomx*lx - 2*NC*x*z*lomx*lx + 2*lx*lopx +  
   2*NC*lx*lopx + 2*x*lx*lopx +  
   2*NC*x*lx*lopx - 4*NC*z*lx*lopx -  
   4*NC*x*z*lx*lopx - (3*NC*lomz)/4. - (NC*x*lomz)/4. +  
   (3*NC*z*lomz)/4. + (NC*x*z*lomz)/4. - NC*lx*lomz -  
   NC*x*lx*lomz - 2*NC*sqrtxz1*x*lspec1 +  
   NC*l2*lspec1 + 4*x*l2*lspec1 -  
   3*NC*x*l2*lspec1 - 6*NC*z*l2*lspec1 -  
   4*x*z*l2*lspec1 + 2*NC*x*z*l2*lspec1 +  
   NC*lx*lspec1 - NC*x*lx*lspec1 -  
   2*NC*z*lx*lspec1 +  
   2*NC*x*z*lx*lspec1 - lz/2. - (7*NC*lz)/2. +  
   (x*lz)/2. + (NC*x*lz)/2. + NC*sqrtxz1*x*lz + z*lz +  
   (5*NC*z*lz)/2. - x*z*lz - 2*NC*x*z*lz + NC*l2*lz -  
   6*x*l2*lz + 2*NC*x*l2*lz + 4*NC*z*l2*lz +  
   6*x*z*l2*lz + 2*NC*x*z*l2*lz - (NC*lomx*lz)/2. -  
   (NC*x*lomx*lz)/2. - NC*x*z*lomx*lz - lx*lz -  
   x*lx*lz + NC*x*lx*lz + NC*z*lx*lz +  
   2*x*z*lx*lz + NC*x*z*lx*lz + NC*x*lomz*lz -  
   2*NC*x*z*lomz*lz + 2*x*lspec1*lz -  
   NC*x*lspec1*lz - 2*NC*z*lspec1*lz -  
   2*x*z*lspec1*lz +  
   (3*NC*sqrtxz3*atanspec2*lspec4)/2. -  
   (3*NC*sqrtxz3*x*z*atanspec2*lspec4)/2. -  
   NC*l2*lspec2 + 4*x*l2*lspec2 -  
   NC*x*l2*lspec2 - 2*NC*z*l2*lspec2 -  
   4*x*z*l2*lspec2 - 2*NC*x*z*l2*lspec2 +  
   2*x*lx*lspec2 - NC*x*lx*lspec2 -  
   2*NC*z*lx*lspec2 - 2*x*z*lx*lspec2 +  
   NC*lspec1*lspec2 -  
   4*x*lspec1*lspec2 +  
   NC*x*lspec1*lspec2 +  
   2*NC*z*lspec1*lspec2 +  
   4*x*z*lspec1*lspec2 +  
   2*NC*x*z*lspec1*lspec2 -  
   NC*lz*lspec2 + 4*x*lz*lspec2 -  
   NC*x*lz*lspec2 - 2*NC*z*lz*lspec2 -  
   4*x*z*lz*lspec2 - 2*NC*x*z*lz*lspec2 -  
   (NC*lx*lxpz)/2. + x*lx*lxpz -  
   (NC*x*lx*lxpz)/2. - NC*z*lx*lxpz -  
   x*z*lx*lxpz - NC*x*z*lx*lxpz +  
   (NC*lz*lxpz)/2. - x*lz*lxpz +  
   (NC*x*lz*lxpz)/2. + NC*z*lz*lxpz +  
   x*z*lz*lxpz + NC*x*z*lz*lxpz -  
   2*x*lx*lopxz + NC*x*lx*lopxz +  
   2*NC*z*lx*lopxz - 2*x*lz*lopxz +  
   NC*x*lz*lopxz + 2*NC*z*lz*lopxz -  
   (NC*lx*lopxzi)/2. - x*lx*lopxzi +  
   (NC*x*lx*lopxzi)/2. + NC*z*lx*lopxzi -  
   x*z*lx*lopxzi - NC*x*z*lx*lopxzi +  
   (NC*lz*lopxzi)/2. + x*lz*lopxzi -  
   (NC*x*lz*lopxzi)/2. - NC*z*lz*lopxzi +  
   x*z*lz*lopxzi + NC*x*z*lz*lopxzi +  
   4*NCi2 - (7*x*NCi2)/2. - 5*z*NCi2 +  
   (19*x*z*NCi2)/4. - 2*Li2mx*NCi2 - 2*x*Li2mx*NCi2 +  
   Li2x*NCi2 + x*Li2x*NCi2 - x*z*Li2x*NCi2 -  
   2*x*li2spec1*NCi2 +  
   2*x*z*li2spec1*NCi2 +  
   2*x*li2spec2*NCi2 -  
   2*x*z*li2spec2*NCi2 +  
   2*x*li2spec19*NCi2 + 2*x*z*li2spec20*NCi2 +  
   2*x*li2spec3*NCi2 -  
   2*x*z*li2spec3*NCi2 -  
   2*x*li2spec4*NCi2 +  
   2*x*z*li2spec4*NCi2 +  
   lx*NCi2 + x*lx*NCi2 - 2*z*lx*NCi2 -  
   2*x*z*lx*NCi2 + 2*x*l2*lx*NCi2 -  
   2*x*z*l2*lx*NCi2 + lomx*lx*NCi2 +  
   x*lomx*lx*NCi2 - x*z*lomx*lx*NCi2 -  
   2*lx*lopx*NCi2 - 2*x*lx*lopx*NCi2 -  
   4*x*l2*lspec1*NCi2 +  
   4*x*z*l2*lspec1*NCi2 + (lz*NCi2)/2. -  
   (x*lz*NCi2)/2. - z*lz*NCi2 + x*z*lz*NCi2 +  
   6*x*l2*lz*NCi2 - 6*x*z*l2*lz*NCi2 +  
   lx*lz*NCi2 + x*lx*lz*NCi2 -  
   2*x*z*lx*lz*NCi2 -  
   2*x*lspec1*lz*NCi2 +  
   2*x*z*lspec1*lz*NCi2 -  
   4*x*l2*lspec2*NCi2 +  
   4*x*z*l2*lspec2*NCi2 -  
   2*x*lx*lspec2*NCi2 +  
   2*x*z*lx*lspec2*NCi2 +  
   4*x*lspec1*lspec2*NCi2 -  
   4*x*z*lspec1*lspec2*NCi2 -  
   4*x*lz*lspec2*NCi2 +  
   4*x*z*lz*lspec2*NCi2 -  
   x*lx*lxpz*NCi2 + x*z*lx*lxpz*NCi2 +  
   x*lz*lxpz*NCi2 - x*z*lz*lxpz*NCi2 +  
   2*x*lx*lopxz*NCi2 + 2*x*lz*lopxz*NCi2 +  
   x*lx*lopxzi*NCi2 +  
   x*z*lx*lopxzi*NCi2 -  
   x*lz*lopxzi*NCi2 -  
   x*z*lz*lopxzi*NCi2 + (13*NCi)/3. -  
   2*x*NCi - (157*z*NCi)/24. + (119*x*z*NCi)/24. +  
   (3*sqrtxz3*itani1*NCi)/4. -  
   (3*sqrtxz3*x*z*itani1*NCi)/4. -  
   (3*sqrtxz3*itani2*NCi)/4. +  
   (3*sqrtxz3*x*z*itani2*NCi)/4. +  
   (3*sqrtxz3*itani3*NCi)/2. -  
   (3*sqrtxz3*x*z*itani3*NCi)/2. -  
   2*Li2mx*NCi - 2*x*Li2mx*NCi + 4*z*Li2mx*NCi +  
   4*x*z*Li2mx*NCi - Li2x*NCi - 2*x*Li2x*NCi +  
   2*x*z*Li2x*NCi + x*li2spec1*NCi +  
   2*z*li2spec1*NCi -  
   x*li2spec2*NCi -  
   2*z*li2spec2*NCi - (Li2z*NCi)/2. -  
   (3*x*Li2z*NCi)/2. + x*z*Li2z*NCi -  
   x*li2spec19*NCi - 2*z*li2spec19*NCi +  
   li2spec20*NCi + 2*x*z*li2spec20*NCi -  
   li2spec3*NCi -  
   2*x*z*li2spec3*NCi +  
   li2spec4*NCi +  
   2*x*z*li2spec4*NCi +  
   x*li2omxzi*NCi -  
   2*x*z*li2omxzi*NCi - 2*sqrtxz1*x*l2*NCi +  
   (3*sqrtxz3*atanspec1*lspec3*NCi)/2. -  
   (3*sqrtxz3*x*z*atanspec1*lspec3*NCi)/2. +  
   (3*lomx*NCi)/4. + (x*lomx*NCi)/4. -  
   (3*z*lomx*NCi)/4. - (x*z*lomx*NCi)/4. +  
   (7*lx*NCi)/4. + (17*x*lx*NCi)/4. -  
   sqrtxz1*x*lx*NCi - (5*z*lx*NCi)/4. -  
   (13*x*z*lx*NCi)/4. + l2*lx*NCi -  
   2*x*l2*lx*NCi - 4*z*l2*lx*NCi +  
   2*x*z*l2*lx*NCi - x*lomx*lx*NCi +  
   2*x*z*lomx*lx*NCi - 2*lx*lopx*NCi -  
   2*x*lx*lopx*NCi + 4*z*lx*lopx*NCi +  
   4*x*z*lx*lopx*NCi + (3*lomz*NCi)/4. +  
   (x*lomz*NCi)/4. - (3*z*lomz*NCi)/4. -  
   (x*z*lomz*NCi)/4. + lx*lomz*NCi +  
   x*lx*lomz*NCi +  
   2*sqrtxz1*x*lspec1*NCi -  
   l2*lspec1*NCi +  
   3*x*l2*lspec1*NCi +  
   6*z*l2*lspec1*NCi -  
   2*x*z*l2*lspec1*NCi -  
   lx*lspec1*NCi +  
   x*lx*lspec1*NCi +  
   2*z*lx*lspec1*NCi -  
   2*x*z*lx*lspec1*NCi + (7*lz*NCi)/2. -  
   (x*lz*NCi)/2. - sqrtxz1*x*lz*NCi -  
   (5*z*lz*NCi)/2. + 2*x*z*lz*NCi -  
   l2*lz*NCi - 2*x*l2*lz*NCi -  
   4*z*l2*lz*NCi - 2*x*z*l2*lz*NCi +  
   (lomx*lz*NCi)/2. + (x*lomx*lz*NCi)/2. +  
   x*z*lomx*lz*NCi - x*lx*lz*NCi -  
   z*lx*lz*NCi - x*z*lx*lz*NCi -  
   x*lomz*lz*NCi + 2*x*z*lomz*lz*NCi +  
   x*lspec1*lz*NCi +  
   2*z*lspec1*lz*NCi -  
   (3*sqrtxz3*atanspec2*lspec4*NCi)/2. +  
   (3*sqrtxz3*x*z*atanspec2*lspec4*NCi)/2. +  
   l2*lspec2*NCi +  
   x*l2*lspec2*NCi +  
   2*z*l2*lspec2*NCi +  
   2*x*z*l2*lspec2*NCi +  
   x*lx*lspec2*NCi +  
   2*z*lx*lspec2*NCi -  
   lspec1*lspec2*NCi -  
   x*lspec1*lspec2*NCi -  
   2*z*lspec1*lspec2*NCi -  
   2*x*z*lspec1*lspec2*NCi +  
   lz*lspec2*NCi +  
   x*lz*lspec2*NCi +  
   2*z*lz*lspec2*NCi +  
   2*x*z*lz*lspec2*NCi +  
   (lx*lxpz*NCi)/2. + (x*lx*lxpz*NCi)/2. +  
   z*lx*lxpz*NCi + x*z*lx*lxpz*NCi -  
   (lz*lxpz*NCi)/2. - (x*lz*lxpz*NCi)/2. -  
   z*lz*lxpz*NCi - x*z*lz*lxpz*NCi -  
   x*lx*lopxz*NCi - 2*z*lx*lopxz*NCi -  
   x*lz*lopxz*NCi - 2*z*lz*lopxz*NCi +  
   (lx*lopxzi*NCi)/2. -  
   (x*lx*lopxzi*NCi)/2. -  
   z*lx*lopxzi*NCi +  
   x*z*lx*lopxzi*NCi -  
   (lz*lopxzi*NCi)/2. +  
   (x*lz*lopxzi*NCi)/2. +  
   z*lz*lopxzi*NCi -  
   x*z*lz*lopxzi*NCi + pi2/3. -  
   (NC*pi2)/12. + (x*pi2)/6. - (5*NC*x*pi2)/12. -  
   (NC*z*pi2)/3. - (x*z*pi2)/3. + (NC*x*z*pi2)/6. -  
   (NCi2*pi2)/3. - (x*NCi2*pi2)/6. +  
   (x*z*NCi2*pi2)/3. + (NCi*pi2)/12. +  
   (5*x*NCi*pi2)/12. + (z*NCi*pi2)/3. -  
   (x*z*NCi*pi2)/6. - (NC*lx*poly2i)/4. +  
   (NC*lz*poly2i)/4. + (lx*NCi*poly2i)/4. -  
   (lz*NCi*poly2i)/4. +  
   8*li2spec1*sqrtxz1i -  
   4*x*li2spec1*sqrtxz1i +  
   2*x*z*li2spec1*sqrtxz1i +  
   8*li2spec2*sqrtxz1i -  
   4*x*li2spec2*sqrtxz1i +  
   2*x*z*li2spec2*sqrtxz1i -  
   16*li2spec21*sqrtxz1i +  
   8*x*li2spec21*sqrtxz1i -  
   4*x*z*li2spec21*sqrtxz1i -  
   8*li2spec22*sqrtxz1i +  
   4*x*li2spec22*sqrtxz1i -  
   2*x*z*li2spec22*sqrtxz1i -  
   16*li2spec23*sqrtxz1i +  
   8*x*li2spec23*sqrtxz1i -  
   4*x*z*li2spec23*sqrtxz1i -  
   8*li2spec3*sqrtxz1i +  
   4*x*li2spec3*sqrtxz1i -  
   2*x*z*li2spec3*sqrtxz1i -  
   8*li2spec4*sqrtxz1i +  
   4*x*li2spec4*sqrtxz1i -  
   2*x*z*li2spec4*sqrtxz1i -  
   32*l2*sqrtxz1i + 16*x*l2*sqrtxz1i -  
   8*x*z*l2*sqrtxz1i - 16*lx*sqrtxz1i +  
   8*x*lx*sqrtxz1i - 4*x*z*lx*sqrtxz1i -  
   16*l2*lomz*sqrtxz1i +  
   8*x*l2*lomz*sqrtxz1i -  
   4*x*z*l2*lomz*sqrtxz1i -  
   8*lx*lomz*sqrtxz1i +  
   4*x*lx*lomz*sqrtxz1i -  
   2*x*z*lx*lomz*sqrtxz1i +  
   32*lspec1*sqrtxz1i -  
   16*x*lspec1*sqrtxz1i +  
   8*x*z*lspec1*sqrtxz1i +  
   32*l2*lspec1*sqrtxz1i -  
   16*x*l2*lspec1*sqrtxz1i +  
   8*x*z*l2*lspec1*sqrtxz1i -  
   8*lx*lspec1*sqrtxz1i +  
   4*x*lx*lspec1*sqrtxz1i -  
   2*x*z*lx*lspec1*sqrtxz1i +  
   16*lomz*lspec1*sqrtxz1i -  
   8*x*lomz*lspec1*sqrtxz1i +  
   4*x*z*lomz*lspec1*sqrtxz1i -  
   16*lz*sqrtxz1i + 8*x*lz*sqrtxz1i -  
   4*x*z*lz*sqrtxz1i - 32*l2*lz*sqrtxz1i +  
   16*x*l2*lz*sqrtxz1i - 8*x*z*l2*lz*sqrtxz1i +  
   4*lx*lz*sqrtxz1i - 2*x*lx*lz*sqrtxz1i +  
   x*z*lx*lz*sqrtxz1i - 8*lomz*lz*sqrtxz1i +  
   4*x*lomz*lz*sqrtxz1i -  
   2*x*z*lomz*lz*sqrtxz1i +  
   16*lspec1*lz*sqrtxz1i -  
   8*x*lspec1*lz*sqrtxz1i +  
   4*x*z*lspec1*lz*sqrtxz1i +  
   16*l2*lspec2*sqrtxz1i -  
   8*x*l2*lspec2*sqrtxz1i +  
   4*x*z*l2*lspec2*sqrtxz1i +  
   8*lx*lspec2*sqrtxz1i -  
   4*x*lx*lspec2*sqrtxz1i +  
   2*x*z*lx*lspec2*sqrtxz1i -  
   16*lspec1*lspec2*sqrtxz1i +  
   8*x*lspec1*lspec2*sqrtxz1i -  
   4*x*z*lspec1*lspec2*sqrtxz1i +  
   16*lz*lspec2*sqrtxz1i -  
   8*x*lz*lspec2*sqrtxz1i +  
   4*x*z*lz*lspec2*sqrtxz1i +  
   8*lomz*lspec7*sqrtxz1i -  
   4*x*lomz*lspec7*sqrtxz1i +  
   2*x*z*lomz*lspec7*sqrtxz1i -  
   8*li2spec1*NCi2*sqrtxz1i +  
   4*x*li2spec1*NCi2*sqrtxz1i -  
   2*x*z*li2spec1*NCi2*sqrtxz1i -  
   8*li2spec2*NCi2*sqrtxz1i +  
   4*x*li2spec2*NCi2*sqrtxz1i -  
   2*x*z*li2spec2*NCi2*sqrtxz1i +  
   16*li2spec21*NCi2*sqrtxz1i -  
   8*x*li2spec21*NCi2*sqrtxz1i +  
   4*x*z*li2spec21*NCi2*sqrtxz1i +  
   8*li2spec22*NCi2*sqrtxz1i -  
   4*x*li2spec22*NCi2*sqrtxz1i +  
   2*x*z*li2spec22*NCi2*sqrtxz1i +  
   16*li2spec23*NCi2*sqrtxz1i -  
   8*x*li2spec23*NCi2*sqrtxz1i +  
   4*x*z*li2spec23*NCi2*sqrtxz1i +  
   8*li2spec3*NCi2* 
    sqrtxz1i - 4*x*li2spec3* 
    NCi2*sqrtxz1i +  
   2*x*z*li2spec3*NCi2* 
    sqrtxz1i + 8*li2spec4* 
    NCi2*sqrtxz1i -  
   4*x*li2spec4*NCi2* 
    sqrtxz1i + 2*x*z*li2spec4* 
    NCi2*sqrtxz1i + 32*l2*NCi2*sqrtxz1i -  
   16*x*l2*NCi2*sqrtxz1i +  
   8*x*z*l2*NCi2*sqrtxz1i +  
   16*lx*NCi2*sqrtxz1i -  
   8*x*lx*NCi2*sqrtxz1i +  
   4*x*z*lx*NCi2*sqrtxz1i +  
   16*l2*lomz*NCi2*sqrtxz1i -  
   8*x*l2*lomz*NCi2*sqrtxz1i +  
   4*x*z*l2*lomz*NCi2*sqrtxz1i +  
   8*lx*lomz*NCi2*sqrtxz1i -  
   4*x*lx*lomz*NCi2*sqrtxz1i +  
   2*x*z*lx*lomz*NCi2*sqrtxz1i -  
   32*lspec1*NCi2*sqrtxz1i +  
   16*x*lspec1*NCi2*sqrtxz1i -  
   8*x*z*lspec1*NCi2*sqrtxz1i -  
   32*l2*lspec1*NCi2*sqrtxz1i +  
   16*x*l2*lspec1*NCi2*sqrtxz1i -  
   8*x*z*l2*lspec1*NCi2*sqrtxz1i +  
   8*lx*lspec1*NCi2*sqrtxz1i -  
   4*x*lx*lspec1*NCi2*sqrtxz1i +  
   2*x*z*lx*lspec1*NCi2*sqrtxz1i -  
   16*lomz*lspec1*NCi2*sqrtxz1i +  
   8*x*lomz*lspec1*NCi2*sqrtxz1i -  
   4*x*z*lomz*lspec1*NCi2*sqrtxz1i +  
   16*lz*NCi2*sqrtxz1i -  
   8*x*lz*NCi2*sqrtxz1i +  
   4*x*z*lz*NCi2*sqrtxz1i +  
   32*l2*lz*NCi2*sqrtxz1i -  
   16*x*l2*lz*NCi2*sqrtxz1i +  
   8*x*z*l2*lz*NCi2*sqrtxz1i -  
   4*lx*lz*NCi2*sqrtxz1i +  
   2*x*lx*lz*NCi2*sqrtxz1i -  
   x*z*lx*lz*NCi2*sqrtxz1i +  
   8*lomz*lz*NCi2*sqrtxz1i -  
   4*x*lomz*lz*NCi2*sqrtxz1i +  
   2*x*z*lomz*lz*NCi2*sqrtxz1i -  
   16*lspec1*lz*NCi2*sqrtxz1i +  
   8*x*lspec1*lz*NCi2*sqrtxz1i -  
   4*x*z*lspec1*lz*NCi2*sqrtxz1i -  
   16*l2*lspec2*NCi2*sqrtxz1i +  
   8*x*l2*lspec2*NCi2*sqrtxz1i -  
   4*x*z*l2*lspec2*NCi2*sqrtxz1i -  
   8*lx*lspec2*NCi2*sqrtxz1i +  
   4*x*lx*lspec2*NCi2*sqrtxz1i -  
   2*x*z*lx*lspec2*NCi2*sqrtxz1i +  
   16*lspec1*lspec2*NCi2*sqrtxz1i -  
   8*x*lspec1*lspec2*NCi2*sqrtxz1i +  
   4*x*z*lspec1*lspec2*NCi2* 
    sqrtxz1i - 16*lz*lspec2*NCi2* 
    sqrtxz1i + 8*x*lz*lspec2*NCi2* 
    sqrtxz1i - 4*x*z*lz*lspec2*NCi2* 
    sqrtxz1i - 8*lomz*lspec7*NCi2* 
    sqrtxz1i + 4*x*lomz*lspec7* 
    NCi2*sqrtxz1i -  
   2*x*z*lomz*lspec7*NCi2* 
    sqrtxz1i - (4*pi2*sqrtxz1i)/3. +  
   (2*x*pi2*sqrtxz1i)/3. - (x*z*pi2*sqrtxz1i)/3. +  
   (4*NCi2*pi2*sqrtxz1i)/3. -  
   (2*x*NCi2*pi2*sqrtxz1i)/3. +  
   (x*z*NCi2*pi2*sqrtxz1i)/3. -  
   (NC*li2spec5*sqrtxz2i)/4. -  
   (3*NC*x*li2spec5*sqrtxz2i)/4. +  
   (NC*z*li2spec5*sqrtxz2i)/2. -  
   2*NC*x*z*li2spec5*sqrtxz2i +  
   (NC*li2spec6*sqrtxz2i)/4. +  
   (3*NC*x*li2spec6*sqrtxz2i)/4. -  
   (NC*z*li2spec6*sqrtxz2i)/2. +  
   2*NC*x*z*li2spec6*sqrtxz2i +  
   (NC*li2spec7*sqrtxz2i)/4. +  
   (3*NC*x*li2spec7*sqrtxz2i)/ 
    4. - (NC*z*li2spec7* 
      sqrtxz2i)/2. + 2*NC*x*z* 
    li2spec7*sqrtxz2i -  
   (NC*li2spec8*sqrtxz2i)/4. -  
   (3*NC*x*li2spec8*sqrtxz2i)/ 
    4. + (NC*z*li2spec8* 
      sqrtxz2i)/2. - 2*NC*x*z* 
    li2spec8*sqrtxz2i +  
   (NC*lx*lspec5*sqrtxz2i)/4. +  
   (3*NC*x*lx*lspec5*sqrtxz2i)/4. -  
   (NC*z*lx*lspec5*sqrtxz2i)/2. +  
   2*NC*x*z*lx*lspec5*sqrtxz2i -  
   (NC*lx*lspec6*sqrtxz2i)/4. -  
   (3*NC*x*lx*lspec6*sqrtxz2i)/4. +  
   (NC*z*lx*lspec6*sqrtxz2i)/2. -  
   2*NC*x*z*lx*lspec6*sqrtxz2i +  
   (li2spec5*NCi*sqrtxz2i)/4. +  
   (3*x*li2spec5*NCi*sqrtxz2i)/4. -  
   (z*li2spec5*NCi*sqrtxz2i)/2. +  
   2*x*z*li2spec5*NCi*sqrtxz2i -  
   (li2spec6*NCi*sqrtxz2i)/4. -  
   (3*x*li2spec6*NCi*sqrtxz2i)/4. +  
   (z*li2spec6*NCi*sqrtxz2i)/2. -  
   2*x*z*li2spec6*NCi*sqrtxz2i -  
   (li2spec7*NCi* 
      sqrtxz2i)/4. - (3*x* 
      li2spec7*NCi* 
      sqrtxz2i)/4. + (z*li2spec7*NCi*sqrtxz2i)/2. -  
   2*x*z*li2spec7*NCi* 
    sqrtxz2i + (li2spec8* 
      NCi*sqrtxz2i)/4. +  
   (3*x*li2spec8*NCi* 
      sqrtxz2i)/4. - (z*li2spec8*NCi*sqrtxz2i)/2. +  
   2*x*z*li2spec8*NCi* 
    sqrtxz2i - (lx*lspec5*NCi* 
      sqrtxz2i)/4. - (3*x*lx*lspec5*NCi* 
      sqrtxz2i)/4. + (z*lx*lspec5*NCi* 
      sqrtxz2i)/2. - 2*x*z*lx*lspec5*NCi* 
    sqrtxz2i + (lx*lspec6*NCi* 
      sqrtxz2i)/4. + (3*x*lx*lspec6*NCi* 
      sqrtxz2i)/4. - (z*lx*lspec6*NCi* 
      sqrtxz2i)/2. + 2*x*z*lx*lspec6*NCi* 
    sqrtxz2i + (NC*x*li2spec5*poly2i* 
      sqrtxz2i)/8. - (NC*x*li2spec6*poly2i* 
      sqrtxz2i)/8. - (NC*x* 
      li2spec7*poly2i* 
      sqrtxz2i)/8. + (NC*x* 
      li2spec8*poly2i* 
      sqrtxz2i)/8. - (NC*x*lx*lspec5*poly2i* 
      sqrtxz2i)/8. + (NC*x*lx*lspec6*poly2i* 
      sqrtxz2i)/8. - (x*li2spec5*NCi* 
      poly2i*sqrtxz2i)/8. +  
   (x*li2spec6*NCi*poly2i*sqrtxz2i)/ 
    8. + (x*li2spec7*NCi* 
      poly2i*sqrtxz2i)/8. -  
   (x*li2spec8*NCi* 
      poly2i*sqrtxz2i)/8. +  
   (x*lx*lspec5*NCi*poly2i*sqrtxz2i)/ 
    8. - (x*lx*lspec6*NCi*poly2i* 
      sqrtxz2i)/8. + 2*Li2x*omxi -  
   (3*NC*Li2x*omxi)/2. - z*Li2x*omxi +  
   3*NC*z*Li2x*omxi -  
   (3*li2spec1*omxi)/2. +  
   (3*NC*li2spec1*omxi)/2. +  
   (3*z*li2spec1*omxi)/2. +  
   NC*z*li2spec1*omxi +  
   (3*li2spec2*omxi)/2. -  
   (3*NC*li2spec2*omxi)/2. -  
   (3*z*li2spec2*omxi)/2. -  
   NC*z*li2spec2*omxi + Li2mz*omxi -  
   z*Li2mz*omxi - (NC*Li2z*omxi)/2. +  
   NC*z*Li2z*omxi + li2spec19*omxi -  
   NC*li2spec19*omxi - z*li2spec19*omxi -  
   2*NC*z*li2spec19*omxi - li2spec20*omxi +  
   NC*li2spec20*omxi +  
   z*li2spec20*omxi +  
   2*NC*z*li2spec20*omxi +  
   (3*li2spec3*omxi)/2. -  
   (NC*li2spec3*omxi)/2. -  
   (3*z*li2spec3*omxi)/2. -  
   3*NC*z*li2spec3*omxi -  
   (3*li2spec4*omxi)/2. +  
   (NC*li2spec4*omxi)/2. +  
   (3*z*li2spec4*omxi)/2. +  
   3*NC*z*li2spec4*omxi +  
   (NC*li2omxzi*omxi)/2. -  
   NC*z*li2omxzi*omxi - sqrtxz1*l2*omxi -  
   3*NC*sqrtxz1*l2*omxi - lx*omxi -  
   (NC*lx*omxi)/4. - (sqrtxz1*lx*omxi)/2. -  
   (3*NC*sqrtxz1*lx*omxi)/2. + (z*lx*omxi)/4. -  
   (NC*z*lx*omxi)/4. + (3*l2*lx*omxi)/2. -  
   (5*NC*l2*lx*omxi)/2. -  
   (3*z*l2*lx*omxi)/2. + NC*z*l2*lx*omxi +  
   2*lomx*lx*omxi -  
   (3*NC*lomx*lx*omxi)/2. -  
   z*lomx*lx*omxi +  
   3*NC*z*lomx*lx*omxi +  
   sqrtxz1*lspec1*omxi +  
   3*NC*sqrtxz1*lspec1*omxi -  
   3*l2*lspec1*omxi +  
   4*NC*l2*lspec1*omxi +  
   3*z*l2*lspec1*omxi +  
   NC*lx*lspec1*omxi -  
   2*NC*z*lx*lspec1*omxi +  
   (lz*omxi)/2. + (3*NC*lz*omxi)/2. -  
   (sqrtxz1*lz*omxi)/2. -  
   (3*NC*sqrtxz1*lz*omxi)/2. + (z*lz*omxi)/2. +  
   (3*NC*z*lz*omxi)/2. + (9*l2*lz*omxi)/2. -  
   (7*NC*l2*lz*omxi)/2. -  
   (9*z*l2*lz*omxi)/2. -  
   5*NC*z*l2*lz*omxi + lx*lz*omxi -  
   (7*NC*lx*lz*omxi)/2. - z*lx*lz*omxi +  
   (NC*z*lx*lz*omxi)/2. -  
   (NC*lomz*lz*omxi)/2. +  
   NC*z*lomz*lz*omxi -  
   (3*lspec1*lz*omxi)/2. +  
   (3*NC*lspec1*lz*omxi)/2. +  
   (3*z*lspec1*lz*omxi)/2. +  
   NC*z*lspec1*lz*omxi +  
   lz*lopz*omxi - z*lz*lopz*omxi -  
   3*l2*lspec2*omxi +  
   2*NC*l2*lspec2*omxi +  
   3*z*l2*lspec2*omxi +  
   4*NC*z*l2*lspec2*omxi -  
   (3*lx*lspec2*omxi)/2. +  
   (3*NC*lx*lspec2*omxi)/2. +  
   (3*z*lx*lspec2*omxi)/2. +  
   NC*z*lx*lspec2*omxi +  
   3*lspec1*lspec2*omxi -  
   2*NC*lspec1*lspec2*omxi -  
   3*z*lspec1*lspec2*omxi -  
   4*NC*z*lspec1*lspec2*omxi -  
   3*lz*lspec2*omxi +  
   2*NC*lz*lspec2*omxi +  
   3*z*lz*lspec2*omxi +  
   4*NC*z*lz*lspec2*omxi -  
   lx*lxpz*omxi + NC*lx*lxpz*omxi +  
   z*lx*lxpz*omxi +  
   2*NC*z*lx*lxpz*omxi + lz*lxpz*omxi -  
   NC*lz*lxpz*omxi - z*lz*lxpz*omxi -  
   2*NC*z*lz*lxpz*omxi +  
   lx*lopxz*omxi - NC*lx*lopxz*omxi -  
   z*lx*lopxz*omxi -  
   2*NC*z*lx*lopxz*omxi +  
   lz*lopxz*omxi - NC*lz*lopxz*omxi -  
   z*lz*lopxz*omxi -  
   2*NC*z*lz*lopxz*omxi -  
   2*Li2x*NCi2*omxi + z*Li2x*NCi2*omxi +  
   (3*li2spec1*NCi2*omxi)/2. -  
   (3*z*li2spec1*NCi2*omxi)/2. -  
   (3*li2spec2*NCi2*omxi)/2. +  
   (3*z*li2spec2*NCi2*omxi)/2. -  
   Li2mz*NCi2*omxi + z*Li2mz*NCi2*omxi -  
   li2spec19*NCi2*omxi +  
   z*li2spec19*NCi2*omxi +  
   li2spec20*NCi2*omxi -  
   z*li2spec20*NCi2*omxi -  
   (3*li2spec3*NCi2* 
      omxi)/2. + (3*z*li2spec3*NCi2*omxi)/2. +  
   (3*li2spec4*NCi2* 
      omxi)/2. - (3*z*li2spec4*NCi2*omxi)/2. +  
   sqrtxz1*l2*NCi2*omxi +  
   lx*NCi2*omxi +  
   (sqrtxz1*lx*NCi2*omxi)/2. -  
   (z*lx*NCi2*omxi)/4. -  
   (3*l2*lx*NCi2*omxi)/2. +  
   (3*z*l2*lx*NCi2*omxi)/2. -  
   2*lomx*lx*NCi2*omxi +  
   z*lomx*lx*NCi2*omxi -  
   sqrtxz1*lspec1*NCi2*omxi +  
   3*l2*lspec1*NCi2*omxi -  
   3*z*l2*lspec1*NCi2*omxi -  
   (lz*NCi2*omxi)/2. +  
   (sqrtxz1*lz*NCi2*omxi)/2. -  
   (z*lz*NCi2*omxi)/2. -  
   (9*l2*lz*NCi2*omxi)/2. +  
   (9*z*l2*lz*NCi2*omxi)/2. -  
   lx*lz*NCi2*omxi +  
   z*lx*lz*NCi2*omxi +  
   (3*lspec1*lz*NCi2*omxi)/2. -  
   (3*z*lspec1*lz*NCi2*omxi)/2. -  
   lz*lopz*NCi2*omxi +  
   z*lz*lopz*NCi2*omxi +  
   3*l2*lspec2*NCi2*omxi -  
   3*z*l2*lspec2*NCi2*omxi +  
   (3*lx*lspec2*NCi2*omxi)/2. -  
   (3*z*lx*lspec2*NCi2*omxi)/2. -  
   3*lspec1*lspec2*NCi2*omxi +  
   3*z*lspec1*lspec2*NCi2*omxi +  
   3*lz*lspec2*NCi2*omxi -  
   3*z*lz*lspec2*NCi2*omxi +  
   lx*lxpz*NCi2*omxi -  
   z*lx*lxpz*NCi2*omxi -  
   lz*lxpz*NCi2*omxi +  
   z*lz*lxpz*NCi2*omxi -  
   lx*lopxz*NCi2*omxi +  
   z*lx*lopxz*NCi2*omxi -  
   lz*lopxz*NCi2*omxi +  
   z*lz*lopxz*NCi2*omxi +  
   (3*Li2x*NCi*omxi)/2. -  
   3*z*Li2x*NCi*omxi -  
   (3*li2spec1*NCi*omxi)/2. -  
   z*li2spec1*NCi*omxi +  
   (3*li2spec2*NCi*omxi)/2. +  
   z*li2spec2*NCi*omxi +  
   (Li2z*NCi*omxi)/2. - z*Li2z*NCi*omxi +  
   li2spec19*NCi*omxi +  
   2*z*li2spec19*NCi*omxi -  
   li2spec20*NCi*omxi -  
   2*z*li2spec20*NCi*omxi +  
   (li2spec3*NCi* 
      omxi)/2. + 3*z*li2spec3*NCi*omxi -  
   (li2spec4*NCi* 
      omxi)/2. - 3*z*li2spec4*NCi*omxi -  
   (li2omxzi*NCi*omxi)/2. +  
   z*li2omxzi*NCi*omxi +  
   3*sqrtxz1*l2*NCi*omxi +  
   (lx*NCi*omxi)/4. +  
   (3*sqrtxz1*lx*NCi*omxi)/2. +  
   (z*lx*NCi*omxi)/4. +  
   (5*l2*lx*NCi*omxi)/2. -  
   z*l2*lx*NCi*omxi +  
   (3*lomx*lx*NCi*omxi)/2. -  
   3*z*lomx*lx*NCi*omxi -  
   3*sqrtxz1*lspec1*NCi*omxi -  
   4*l2*lspec1*NCi*omxi -  
   lx*lspec1*NCi*omxi +  
   2*z*lx*lspec1*NCi*omxi -  
   (3*lz*NCi*omxi)/2. +  
   (3*sqrtxz1*lz*NCi*omxi)/2. -  
   (3*z*lz*NCi*omxi)/2. +  
   (7*l2*lz*NCi*omxi)/2. +  
   5*z*l2*lz*NCi*omxi +  
   (7*lx*lz*NCi*omxi)/2. -  
   (z*lx*lz*NCi*omxi)/2. +  
   (lomz*lz*NCi*omxi)/2. -  
   z*lomz*lz*NCi*omxi -  
   (3*lspec1*lz*NCi*omxi)/2. -  
   z*lspec1*lz*NCi*omxi -  
   2*l2*lspec2*NCi*omxi -  
   4*z*l2*lspec2*NCi*omxi -  
   (3*lx*lspec2*NCi*omxi)/2. -  
   z*lx*lspec2*NCi*omxi +  
   2*lspec1*lspec2*NCi*omxi +  
   4*z*lspec1*lspec2*NCi*omxi -  
   2*lz*lspec2*NCi*omxi -  
   4*z*lz*lspec2*NCi*omxi -  
   lx*lxpz*NCi*omxi -  
   2*z*lx*lxpz*NCi*omxi +  
   lz*lxpz*NCi*omxi +  
   2*z*lz*lxpz*NCi*omxi +  
   lx*lopxz*NCi*omxi +  
   2*z*lx*lopxz*NCi*omxi +  
   lz*lopxz*NCi*omxi +  
   2*z*lz*lopxz*NCi*omxi -  
   (pi2*omxi)/4. + (5*NC*pi2*omxi)/12. +  
   (z*pi2*omxi)/12. - (5*NC*z*pi2*omxi)/6. +  
   (NCi2*pi2*omxi)/4. -  
   (z*NCi2*pi2*omxi)/12. -  
   (5*NCi*pi2*omxi)/12. +  
   (5*z*NCi*pi2*omxi)/6. -  
   4*li2spec1*sqrtxz1i*omxi -  
   2*z*li2spec1*sqrtxz1i*omxi -  
   4*li2spec2*sqrtxz1i*omxi -  
   2*z*li2spec2*sqrtxz1i*omxi +  
   8*li2spec21*sqrtxz1i*omxi +  
   4*z*li2spec21*sqrtxz1i*omxi +  
   4*li2spec22*sqrtxz1i*omxi +  
   2*z*li2spec22*sqrtxz1i*omxi +  
   8*li2spec23*sqrtxz1i*omxi +  
   4*z*li2spec23*sqrtxz1i*omxi +  
   4*li2spec3*sqrtxz1i* 
    omxi + 2*z*li2spec3* 
    sqrtxz1i*omxi +  
   4*li2spec4*sqrtxz1i* 
    omxi + 2*z*li2spec4* 
    sqrtxz1i*omxi + 16*l2*sqrtxz1i*omxi +  
   8*z*l2*sqrtxz1i*omxi +  
   8*lx*sqrtxz1i*omxi +  
   4*z*lx*sqrtxz1i*omxi +  
   8*l2*lomz*sqrtxz1i*omxi +  
   4*z*l2*lomz*sqrtxz1i*omxi +  
   4*lx*lomz*sqrtxz1i*omxi +  
   2*z*lx*lomz*sqrtxz1i*omxi -  
   16*lspec1*sqrtxz1i*omxi -  
   8*z*lspec1*sqrtxz1i*omxi -  
   16*l2*lspec1*sqrtxz1i*omxi -  
   8*z*l2*lspec1*sqrtxz1i*omxi +  
   4*lx*lspec1*sqrtxz1i*omxi +  
   2*z*lx*lspec1*sqrtxz1i*omxi -  
   8*lomz*lspec1*sqrtxz1i*omxi -  
   4*z*lomz*lspec1*sqrtxz1i*omxi +  
   8*lz*sqrtxz1i*omxi +  
   4*z*lz*sqrtxz1i*omxi +  
   16*l2*lz*sqrtxz1i*omxi +  
   8*z*l2*lz*sqrtxz1i*omxi -  
   2*lx*lz*sqrtxz1i*omxi -  
   z*lx*lz*sqrtxz1i*omxi +  
   4*lomz*lz*sqrtxz1i*omxi +  
   2*z*lomz*lz*sqrtxz1i*omxi -  
   8*lspec1*lz*sqrtxz1i*omxi -  
   4*z*lspec1*lz*sqrtxz1i*omxi -  
   8*l2*lspec2*sqrtxz1i*omxi -  
   4*z*l2*lspec2*sqrtxz1i*omxi -  
   4*lx*lspec2*sqrtxz1i*omxi -  
   2*z*lx*lspec2*sqrtxz1i*omxi +  
   8*lspec1*lspec2*sqrtxz1i* 
    omxi + 4*z*lspec1*lspec2* 
    sqrtxz1i*omxi -  
   8*lz*lspec2*sqrtxz1i*omxi -  
   4*z*lz*lspec2*sqrtxz1i*omxi -  
   4*lomz*lspec7*sqrtxz1i* 
    omxi - 2*z*lomz*lspec7* 
    sqrtxz1i*omxi +  
   4*li2spec1*NCi2*sqrtxz1i*omxi +  
   2*z*li2spec1*NCi2*sqrtxz1i* 
    omxi + 4*li2spec2*NCi2*sqrtxz1i* 
    omxi + 2*z*li2spec2*NCi2* 
    sqrtxz1i*omxi -  
   8*li2spec21*NCi2*sqrtxz1i*omxi -  
   4*z*li2spec21*NCi2*sqrtxz1i* 
    omxi - 4*li2spec22*NCi2* 
    sqrtxz1i*omxi -  
   2*z*li2spec22*NCi2*sqrtxz1i*omxi -  
   8*li2spec23*NCi2*sqrtxz1i*omxi -  
   4*z*li2spec23*NCi2*sqrtxz1i*omxi -  
   4*li2spec3*NCi2* 
    sqrtxz1i*omxi -  
   2*z*li2spec3*NCi2* 
    sqrtxz1i*omxi -  
   4*li2spec4*NCi2* 
    sqrtxz1i*omxi -  
   2*z*li2spec4*NCi2* 
    sqrtxz1i*omxi -  
   16*l2*NCi2*sqrtxz1i*omxi -  
   8*z*l2*NCi2*sqrtxz1i*omxi -  
   8*lx*NCi2*sqrtxz1i*omxi -  
   4*z*lx*NCi2*sqrtxz1i*omxi -  
   8*l2*lomz*NCi2*sqrtxz1i*omxi -  
   4*z*l2*lomz*NCi2*sqrtxz1i*omxi -  
   4*lx*lomz*NCi2*sqrtxz1i*omxi -  
   2*z*lx*lomz*NCi2*sqrtxz1i*omxi +  
   16*lspec1*NCi2*sqrtxz1i*omxi +  
   8*z*lspec1*NCi2*sqrtxz1i*omxi +  
   16*l2*lspec1*NCi2*sqrtxz1i*omxi +  
   8*z*l2*lspec1*NCi2*sqrtxz1i*omxi -  
   4*lx*lspec1*NCi2*sqrtxz1i*omxi -  
   2*z*lx*lspec1*NCi2*sqrtxz1i*omxi +  
   8*lomz*lspec1*NCi2*sqrtxz1i* 
    omxi + 4*z*lomz*lspec1*NCi2* 
    sqrtxz1i*omxi -  
   8*lz*NCi2*sqrtxz1i*omxi -  
   4*z*lz*NCi2*sqrtxz1i*omxi -  
   16*l2*lz*NCi2*sqrtxz1i*omxi -  
   8*z*l2*lz*NCi2*sqrtxz1i*omxi +  
   2*lx*lz*NCi2*sqrtxz1i*omxi +  
   z*lx*lz*NCi2*sqrtxz1i*omxi -  
   4*lomz*lz*NCi2*sqrtxz1i*omxi -  
   2*z*lomz*lz*NCi2*sqrtxz1i*omxi +  
   8*lspec1*lz*NCi2*sqrtxz1i*omxi +  
   4*z*lspec1*lz*NCi2*sqrtxz1i*omxi +  
   8*l2*lspec2*NCi2*sqrtxz1i*omxi +  
   4*z*l2*lspec2*NCi2*sqrtxz1i*omxi +  
   4*lx*lspec2*NCi2*sqrtxz1i*omxi +  
   2*z*lx*lspec2*NCi2*sqrtxz1i*omxi -  
   8*lspec1*lspec2*NCi2*sqrtxz1i* 
    omxi - 4*z*lspec1*lspec2*NCi2* 
    sqrtxz1i*omxi +  
   8*lz*lspec2*NCi2*sqrtxz1i*omxi +  
   4*z*lz*lspec2*NCi2*sqrtxz1i*omxi +  
   4*lomz*lspec7*NCi2*sqrtxz1i* 
    omxi + 2*z*lomz*lspec7*NCi2* 
    sqrtxz1i*omxi +  
   (2*pi2*sqrtxz1i*omxi)/3. +  
   (z*pi2*sqrtxz1i*omxi)/3. -  
   (2*NCi2*pi2*sqrtxz1i*omxi)/3. -  
   (z*NCi2*pi2*sqrtxz1i*omxi)/3. -  
   2*Li2mx*xi2 + NC*Li2mx*xi2 - 2*z*Li2mx*xi2 -  
   2*NC*z*Li2mx*xi2 + 2*Li2x*xi2 - NC*Li2x*xi2 +  
   2*z*Li2x*xi2 + 2*NC*z*Li2x*xi2 + li2spec19*xi2 -  
   (NC*li2spec19*xi2)/2. + z*li2spec19*xi2 +  
   NC*z*li2spec19*xi2 + li2spec20*xi2 -  
   (NC*li2spec20*xi2)/2. + z*li2spec20*xi2 +  
   NC*z*li2spec20*xi2 - 4*lx*xi2 +  
   2*NC*lx*xi2 - 4*z*lx*xi2 - 4*NC*z*lx*xi2 +  
   2*lomx*lx*xi2 - NC*lomx*lx*xi2 +  
   2*z*lomx*lx*xi2 + 2*NC*z*lomx*lx*xi2 -  
   2*lx*lopx*xi2 + NC*lx*lopx*xi2 -  
   2*z*lx*lopx*xi2 - 2*NC*z*lx*lopx*xi2 -  
   lx*lz*xi2 + (NC*lx*lz*xi2)/2. -  
   z*lx*lz*xi2 - NC*z*lx*lz*xi2 +  
   lx*lopxz*xi2 - (NC*lx*lopxz*xi2)/2. +  
   z*lx*lopxz*xi2 + NC*z*lx*lopxz*xi2 +  
   lz*lopxz*xi2 - (NC*lz*lopxz*xi2)/2. +  
   z*lz*lopxz*xi2 + NC*z*lz*lopxz*xi2 +  
   lx*lopxzi*xi2 -  
   (NC*lx*lopxzi*xi2)/2. +  
   z*lx*lopxzi*xi2 +  
   NC*z*lx*lopxzi*xi2 -  
   lz*lopxzi*xi2 +  
   (NC*lz*lopxzi*xi2)/2. -  
   z*lz*lopxzi*xi2 -  
   NC*z*lz*lopxzi*xi2 +  
   2*Li2mx*NCi2*xi2 + 2*z*Li2mx*NCi2*xi2 -  
   2*Li2x*NCi2*xi2 - 2*z*Li2x*NCi2*xi2 -  
   li2spec19*NCi2*xi2 - z*li2spec19*NCi2*xi2 -  
   li2spec20*NCi2*xi2 -  
   z*li2spec20*NCi2*xi2 +  
   4*lx*NCi2*xi2 + 4*z*lx*NCi2*xi2 -  
   2*lomx*lx*NCi2*xi2 -  
   2*z*lomx*lx*NCi2*xi2 +  
   2*lx*lopx*NCi2*xi2 +  
   2*z*lx*lopx*NCi2*xi2 +  
   lx*lz*NCi2*xi2 +  
   z*lx*lz*NCi2*xi2 -  
   lx*lopxz*NCi2*xi2 -  
   z*lx*lopxz*NCi2*xi2 -  
   lz*lopxz*NCi2*xi2 -  
   z*lz*lopxz*NCi2*xi2 -  
   lx*lopxzi*NCi2*xi2 -  
   z*lx*lopxzi*NCi2*xi2 +  
   lz*lopxzi*NCi2*xi2 +  
   z*lz*lopxzi*NCi2*xi2 -  
   Li2mx*NCi*xi2 + 2*z*Li2mx*NCi*xi2 +  
   Li2x*NCi*xi2 - 2*z*Li2x*NCi*xi2 +  
   (li2spec19*NCi*xi2)/2. -  
   z*li2spec19*NCi*xi2 +  
   (li2spec20*NCi*xi2)/2. -  
   z*li2spec20*NCi*xi2 -  
   2*lx*NCi*xi2 + 4*z*lx*NCi*xi2 +  
   lomx*lx*NCi*xi2 -  
   2*z*lomx*lx*NCi*xi2 -  
   lx*lopx*NCi*xi2 +  
   2*z*lx*lopx*NCi*xi2 -  
   (lx*lz*NCi*xi2)/2. +  
   z*lx*lz*NCi*xi2 +  
   (lx*lopxz*NCi*xi2)/2. -  
   z*lx*lopxz*NCi*xi2 +  
   (lz*lopxz*NCi*xi2)/2. -  
   z*lz*lopxz*NCi*xi2 +  
   (lx*lopxzi*NCi*xi2)/2. -  
   z*lx*lopxzi*NCi*xi2 -  
   (lz*lopxzi*NCi*xi2)/2. +  
   z*lz*lopxzi*NCi*xi2 -  
   (pi2*xi2)/3. + (NC*pi2*xi2)/6. -  
   (z*pi2*xi2)/3. - (NC*z*pi2*xi2)/3. +  
   (NCi2*pi2*xi2)/3. +  
   (z*NCi2*pi2*xi2)/3. -  
   (NCi*pi2*xi2)/6. +  
   (z*NCi*pi2*xi2)/3. - (13*NC*xi)/9. +  
   (3*NC*sqrtxz3*z*itani1*xi)/4. -  
   (3*NC*sqrtxz3*z*itani2*xi)/4. +  
   (3*NC*sqrtxz3*z*itani3*xi)/2. -  
   NC*Li2mx*xi + 2*NC*z*Li2mx*xi + NC*Li2x*xi -  
   2*NC*z*Li2x*xi + (NC*li2spec19*xi)/2. -  
   NC*z*li2spec19*xi + (NC*li2spec20*xi)/2. -  
   NC*z*li2spec20*xi +  
   (3*NC*sqrtxz3*z*atanspec1*lspec3*xi)/2. -  
   (2*NC*lomx*xi)/3. - (9*NC*lx*xi)/4. +  
   4*NC*z*lx*xi + NC*lomx*lx*xi -  
   2*NC*z*lomx*lx*xi - NC*lx*lopx*xi +  
   2*NC*z*lx*lopx*xi - (2*NC*lomz*xi)/3. -  
   (11*NC*lz*xi)/12. - (NC*lx*lz*xi)/2. +  
   NC*z*lx*lz*xi -  
   (3*NC*sqrtxz3*z*atanspec2*lspec4*xi)/2. +  
   (NC*lx*lopxz*xi)/2. -  
   NC*z*lx*lopxz*xi +  
   (NC*lz*lopxz*xi)/2. -  
   NC*z*lz*lopxz*xi +  
   (NC*lx*lopxzi*xi)/2. -  
   NC*z*lx*lopxzi*xi -  
   (NC*lz*lopxzi*xi)/2. +  
   NC*z*lz*lopxzi*xi +  
   (13*NCi*xi)/9. -  
   (3*sqrtxz3*z*itani1*NCi*xi)/4. +  
   (3*sqrtxz3*z*itani2*NCi*xi)/4. -  
   (3*sqrtxz3*z*itani3*NCi*xi)/2. +  
   Li2mx*NCi*xi - 2*z*Li2mx*NCi*xi -  
   Li2x*NCi*xi + 2*z*Li2x*NCi*xi -  
   (li2spec19*NCi*xi)/2. +  
   z*li2spec19*NCi*xi -  
   (li2spec20*NCi*xi)/2. +  
   z*li2spec20*NCi*xi -  
   (3*sqrtxz3*z*atanspec1*lspec3*NCi*xi)/2. +  
   (2*lomx*NCi*xi)/3. +  
   (9*lx*NCi*xi)/4. - 4*z*lx*NCi*xi -  
   lomx*lx*NCi*xi +  
   2*z*lomx*lx*NCi*xi +  
   lx*lopx*NCi*xi -  
   2*z*lx*lopx*NCi*xi +  
   (2*lomz*NCi*xi)/3. +  
   (11*lz*NCi*xi)/12. +  
   (lx*lz*NCi*xi)/2. -  
   z*lx*lz*NCi*xi +  
   (3*sqrtxz3*z*atanspec2*lspec4*NCi*xi)/2. -  
   (lx*lopxz*NCi*xi)/2. +  
   z*lx*lopxz*NCi*xi -  
   (lz*lopxz*NCi*xi)/2. +  
   z*lz*lopxz*NCi*xi -  
   (lx*lopxzi*NCi*xi)/2. +  
   z*lx*lopxzi*NCi*xi +  
   (lz*lopxzi*NCi*xi)/2. -  
   z*lz*lopxzi*NCi*xi -  
   (NC*pi2*xi)/6. + (NC*z*pi2*xi)/3. +  
   (NCi*pi2*xi)/6. -  
   (z*NCi*pi2*xi)/3. +  
   (NC*lx*poly2i*xi)/4. +  
   (NC*lz*poly2i*xi)/4. -  
   (lx*NCi*poly2i*xi)/4. -  
   (lz*NCi*poly2i*xi)/4. +  
   (NC*li2spec5*sqrtxz2i*xi)/8. -  
   (NC*li2spec6*sqrtxz2i*xi)/8. -  
   (NC*li2spec7*sqrtxz2i* 
      xi)/8. + (NC*li2spec8* 
      sqrtxz2i*xi)/8. -  
   (NC*lx*lspec5*sqrtxz2i*xi)/8. +  
   (NC*lx*lspec6*sqrtxz2i*xi)/8. -  
   (li2spec5*NCi*sqrtxz2i*xi)/8. +  
   (li2spec6*NCi*sqrtxz2i*xi)/8. +  
   (li2spec7*NCi* 
      sqrtxz2i*xi)/8. -  
   (li2spec8*NCi* 
      sqrtxz2i*xi)/8. +  
   (lx*lspec5*NCi*sqrtxz2i*xi)/8. -  
   (lx*lspec6*NCi*sqrtxz2i*xi)/8. -  
   (NC*li2spec5*poly2i*sqrtxz2i*xi)/ 
    8. + (NC*li2spec6*poly2i*sqrtxz2i* 
      xi)/8. + (NC*li2spec7* 
      poly2i*sqrtxz2i*xi)/8. -  
   (NC*li2spec8*poly2i* 
      sqrtxz2i*xi)/8. +  
   (NC*lx*lspec5*poly2i*sqrtxz2i*xi)/ 
    8. - (NC*lx*lspec6*poly2i*sqrtxz2i* 
      xi)/8. + (li2spec5*NCi*poly2i* 
      sqrtxz2i*xi)/8. -  
   (li2spec6*NCi*poly2i*sqrtxz2i* 
      xi)/8. - (li2spec7* 
      NCi*poly2i*sqrtxz2i*xi)/8. +  
   (li2spec8*NCi*poly2i* 
      sqrtxz2i*xi)/8. -  
   (lx*lspec5*NCi*poly2i*sqrtxz2i* 
      xi)/8. + (lx*lspec6*NCi*poly2i* 
      sqrtxz2i*xi)/8. + (22*NC*x2)/9. +  
   (2*NC*lomx*x2)/3. - (5*NC*lx*x2)/4. +  
   (2*NC*lomz*x2)/3. - (NC*lz*x2)/12. -  
   (22*NCi*x2)/9. - (2*lomx*NCi*x2)/3. +  
   (5*lx*NCi*x2)/4. -  
   (2*lomz*NCi*x2)/3. + (lz*NCi*x2)/12. +  
   8*li2spec1*sqrtxz1i*x2 +  
   8*li2spec2*sqrtxz1i*x2 -  
   16*li2spec21*sqrtxz1i*x2 -  
   8*li2spec22*sqrtxz1i*x2 -  
   16*li2spec23*sqrtxz1i*x2 -  
   8*li2spec3*sqrtxz1i* 
    x2 - 8*li2spec4* 
    sqrtxz1i*x2 - 32*l2*sqrtxz1i*x2 -  
   16*lx*sqrtxz1i*x2 -  
   16*l2*lomz*sqrtxz1i*x2 -  
   8*lx*lomz*sqrtxz1i*x2 +  
   32*lspec1*sqrtxz1i*x2 +  
   32*l2*lspec1*sqrtxz1i*x2 -  
   8*lx*lspec1*sqrtxz1i*x2 +  
   16*lomz*lspec1*sqrtxz1i*x2 -  
   16*lz*sqrtxz1i*x2 -  
   32*l2*lz*sqrtxz1i*x2 +  
   4*lx*lz*sqrtxz1i*x2 -  
   8*lomz*lz*sqrtxz1i*x2 +  
   16*lspec1*lz*sqrtxz1i*x2 +  
   16*l2*lspec2*sqrtxz1i*x2 +  
   8*lx*lspec2*sqrtxz1i*x2 -  
   16*lspec1*lspec2*sqrtxz1i*x2 +  
   16*lz*lspec2*sqrtxz1i*x2 +  
   8*lomz*lspec7*sqrtxz1i*x2 -  
   8*li2spec1*NCi2*sqrtxz1i*x2 -  
   8*li2spec2*NCi2*sqrtxz1i*x2 +  
   16*li2spec21*NCi2*sqrtxz1i*x2 +  
   8*li2spec22*NCi2*sqrtxz1i*x2 +  
   16*li2spec23*NCi2*sqrtxz1i*x2 +  
   8*li2spec3*NCi2* 
    sqrtxz1i*x2 + 8* 
    li2spec4*NCi2* 
    sqrtxz1i*x2 + 32*l2*NCi2*sqrtxz1i* 
    x2 + 16*lx*NCi2*sqrtxz1i*x2 +  
   16*l2*lomz*NCi2*sqrtxz1i*x2 +  
   8*lx*lomz*NCi2*sqrtxz1i*x2 -  
   32*lspec1*NCi2*sqrtxz1i*x2 -  
   32*l2*lspec1*NCi2*sqrtxz1i*x2 +  
   8*lx*lspec1*NCi2*sqrtxz1i*x2 -  
   16*lomz*lspec1*NCi2*sqrtxz1i*x2 +  
   16*lz*NCi2*sqrtxz1i*x2 +  
   32*l2*lz*NCi2*sqrtxz1i*x2 -  
   4*lx*lz*NCi2*sqrtxz1i*x2 +  
   8*lomz*lz*NCi2*sqrtxz1i*x2 -  
   16*lspec1*lz*NCi2*sqrtxz1i*x2 -  
   16*l2*lspec2*NCi2*sqrtxz1i*x2 -  
   8*lx*lspec2*NCi2*sqrtxz1i*x2 +  
   16*lspec1*lspec2*NCi2*sqrtxz1i* 
    x2 - 16*lz*lspec2*NCi2*sqrtxz1i* 
    x2 - 8*lomz*lspec7*NCi2* 
    sqrtxz1i*x2 - (4*pi2*sqrtxz1i*x2)/3. +  
   (4*NCi2*pi2*sqrtxz1i*x2)/3. -  
   (NC*li2spec5*sqrtxz2i*x2)/4. +  
   (NC*z*li2spec5*sqrtxz2i*x2)/2. +  
   (NC*li2spec6*sqrtxz2i*x2)/4. -  
   (NC*z*li2spec6*sqrtxz2i*x2)/2. +  
   (NC*li2spec7*sqrtxz2i* 
      x2)/4. - (NC*z*li2spec7* 
      sqrtxz2i*x2)/2. -  
   (NC*li2spec8*sqrtxz2i* 
      x2)/4. + (NC*z*li2spec8* 
      sqrtxz2i*x2)/2. +  
   (NC*lx*lspec5*sqrtxz2i*x2)/4. -  
   (NC*z*lx*lspec5*sqrtxz2i*x2)/2. -  
   (NC*lx*lspec6*sqrtxz2i*x2)/4. +  
   (NC*z*lx*lspec6*sqrtxz2i*x2)/2. +  
   (li2spec5*NCi*sqrtxz2i*x2)/4. -  
   (z*li2spec5*NCi*sqrtxz2i*x2)/2. -  
   (li2spec6*NCi*sqrtxz2i*x2)/4. +  
   (z*li2spec6*NCi*sqrtxz2i*x2)/2. -  
   (li2spec7*NCi* 
      sqrtxz2i*x2)/4. +  
   (z*li2spec7*NCi* 
      sqrtxz2i*x2)/2. +  
   (li2spec8*NCi* 
      sqrtxz2i*x2)/4. -  
   (z*li2spec8*NCi* 
      sqrtxz2i*x2)/2. -  
   (lx*lspec5*NCi*sqrtxz2i*x2)/4. +  
   (z*lx*lspec5*NCi*sqrtxz2i*x2)/2. +  
   (lx*lspec6*NCi*sqrtxz2i*x2)/4. -  
   (z*lx*lspec6*NCi*sqrtxz2i*x2)/2. -  
   (NC*lx*poly2i*x3)/4. -  
   (NC*lz*poly2i*x3)/4. +  
   (lx*NCi*poly2i*x3)/4. +  
   (lz*NCi*poly2i*x3)/4. +  
   (NC*li2spec5*sqrtxz2i*x3)/8. -  
   (NC*li2spec6*sqrtxz2i*x3)/8. -  
   (NC*li2spec7*sqrtxz2i* 
      x3)/8. + (NC*li2spec8* 
      sqrtxz2i*x3)/8. -  
   (NC*lx*lspec5*sqrtxz2i*x3)/8. +  
   (NC*lx*lspec6*sqrtxz2i*x3)/8. -  
   (li2spec5*NCi*sqrtxz2i*x3)/8. +  
   (li2spec6*NCi*sqrtxz2i*x3)/8. +  
   (li2spec7*NCi* 
      sqrtxz2i*x3)/8. -  
   (li2spec8*NCi* 
      sqrtxz2i*x3)/8. +  
   (lx*lspec5*NCi*sqrtxz2i*x3)/8. -  
   (lx*lspec6*NCi*sqrtxz2i*x3)/8. +  
   (NC*li2spec5*poly2i*sqrtxz2i*x3)/ 
    8. - (NC*li2spec6*poly2i*sqrtxz2i* 
      x3)/8. - (NC*li2spec7* 
      poly2i*sqrtxz2i*x3)/8. +  
   (NC*li2spec8*poly2i* 
      sqrtxz2i*x3)/8. -  
   (NC*lx*lspec5*poly2i*sqrtxz2i*x3)/ 
    8. + (NC*lx*lspec6*poly2i*sqrtxz2i* 
      x3)/8. - (li2spec5*NCi*poly2i* 
      sqrtxz2i*x3)/8. +  
   (li2spec6*NCi*poly2i*sqrtxz2i* 
      x3)/8. + (li2spec7* 
      NCi*poly2i*sqrtxz2i*x3)/8. -  
   (li2spec8*NCi*poly2i* 
      sqrtxz2i*x3)/8. +  
   (lx*lspec5*NCi*poly2i*sqrtxz2i* 
      x3)/8. - (lx*lspec6*NCi*poly2i* 
      sqrtxz2i*x3)/8. + (NC*lx*poly2i*x4)/4. -  
   (NC*lz*poly2i*x4)/4. -  
   (lx*NCi*poly2i*x4)/4. +  
   (lz*NCi*poly2i*x4)/4. -  
   (NC*li2spec5*poly2i*sqrtxz2i*x5)/ 
    8. + (NC*li2spec6*poly2i*sqrtxz2i* 
      x5)/8. + (NC*li2spec7* 
      poly2i*sqrtxz2i*x5)/8. -  
   (NC*li2spec8*poly2i* 
      sqrtxz2i*x5)/8. +  
   (NC*lx*lspec5*poly2i*sqrtxz2i*x5)/ 
    8. - (NC*lx*lspec6*poly2i*sqrtxz2i* 
      x5)/8. + (li2spec5*NCi*poly2i* 
      sqrtxz2i*x5)/8. -  
   (li2spec6*NCi*poly2i*sqrtxz2i* 
      x5)/8. - (li2spec7* 
      NCi*poly2i*sqrtxz2i*x5)/8. +  
   (li2spec8*NCi*poly2i* 
      sqrtxz2i*x5)/8. -  
   (lx*lspec5*NCi*poly2i*sqrtxz2i* 
      x5)/8. + (lx*lspec6*NCi*poly2i* 
      sqrtxz2i*x5)/8. + NC*Li2mx*opxi -  
   2*NC*z*Li2mx*opxi - NC*Li2x*opxi +  
   2*NC*z*Li2x*opxi - li2spec19*opxi +  
   (NC*li2spec19*opxi)/2. - z*li2spec19*opxi -  
   NC*z*li2spec19*opxi - li2spec20*opxi +  
   (NC*li2spec20*opxi)/2. -  
   z*li2spec20*opxi -  
   NC*z*li2spec20*opxi + 2*NC*lx*opxi -  
   4*NC*z*lx*opxi - NC*lomx*lx*opxi +  
   2*NC*z*lomx*lx*opxi +  
   NC*lx*lopx*opxi -  
   2*NC*z*lx*lopx*opxi + lx*lz*opxi -  
   (NC*lx*lz*opxi)/2. + z*lx*lz*opxi +  
   NC*z*lx*lz*opxi - lx*lopxz*opxi +  
   (NC*lx*lopxz*opxi)/2. -  
   z*lx*lopxz*opxi -  
   NC*z*lx*lopxz*opxi -  
   lz*lopxz*opxi +  
   (NC*lz*lopxz*opxi)/2. -  
   z*lz*lopxz*opxi -  
   NC*z*lz*lopxz*opxi -  
   lx*lopxzi*opxi +  
   (NC*lx*lopxzi*opxi)/2. -  
   z*lx*lopxzi*opxi -  
   NC*z*lx*lopxzi*opxi +  
   lz*lopxzi*opxi -  
   (NC*lz*lopxzi*opxi)/2. +  
   z*lz*lopxzi*opxi +  
   NC*z*lz*lopxzi*opxi +  
   li2spec19*NCi2*opxi +  
   z*li2spec19*NCi2*opxi +  
   li2spec20*NCi2*opxi +  
   z*li2spec20*NCi2*opxi -  
   lx*lz*NCi2*opxi -  
   z*lx*lz*NCi2*opxi +  
   lx*lopxz*NCi2*opxi +  
   z*lx*lopxz*NCi2*opxi +  
   lz*lopxz*NCi2*opxi +  
   z*lz*lopxz*NCi2*opxi +  
   lx*lopxzi*NCi2*opxi +  
   z*lx*lopxzi*NCi2*opxi -  
   lz*lopxzi*NCi2*opxi -  
   z*lz*lopxzi*NCi2*opxi -  
   Li2mx*NCi*opxi + 2*z*Li2mx*NCi*opxi +  
   Li2x*NCi*opxi - 2*z*Li2x*NCi*opxi -  
   (li2spec19*NCi*opxi)/2. +  
   z*li2spec19*NCi*opxi -  
   (li2spec20*NCi*opxi)/2. +  
   z*li2spec20*NCi*opxi -  
   2*lx*NCi*opxi + 4*z*lx*NCi*opxi +  
   lomx*lx*NCi*opxi -  
   2*z*lomx*lx*NCi*opxi -  
   lx*lopx*NCi*opxi +  
   2*z*lx*lopx*NCi*opxi +  
   (lx*lz*NCi*opxi)/2. -  
   z*lx*lz*NCi*opxi -  
   (lx*lopxz*NCi*opxi)/2. +  
   z*lx*lopxz*NCi*opxi -  
   (lz*lopxz*NCi*opxi)/2. +  
   z*lz*lopxz*NCi*opxi -  
   (lx*lopxzi*NCi*opxi)/2. +  
   z*lx*lopxzi*NCi*opxi +  
   (lz*lopxzi*NCi*opxi)/2. -  
   z*lz*lopxzi*NCi*opxi -  
   (pi2*opxi)/6. + (NC*pi2*opxi)/3. -  
   (z*pi2*opxi)/6. - (2*NC*z*pi2*opxi)/3. +  
   (NCi2*pi2*opxi)/6. +  
   (z*NCi2*pi2*opxi)/6. -  
   (NCi*pi2*opxi)/3. +  
   (2*z*NCi*pi2*opxi)/3. +  
   2*Li2mx*xi2*opxi - NC*Li2mx*xi2*opxi +  
   2*z*Li2mx*xi2*opxi +  
   2*NC*z*Li2mx*xi2*opxi -  
   2*Li2x*xi2*opxi + NC*Li2x*xi2*opxi -  
   2*z*Li2x*xi2*opxi -  
   2*NC*z*Li2x*xi2*opxi -  
   li2spec19*xi2*opxi +  
   (NC*li2spec19*xi2*opxi)/2. -  
   z*li2spec19*xi2*opxi -  
   NC*z*li2spec19*xi2*opxi -  
   li2spec20*xi2*opxi +  
   (NC*li2spec20*xi2*opxi)/2. -  
   z*li2spec20*xi2*opxi -  
   NC*z*li2spec20*xi2*opxi +  
   4*lx*xi2*opxi - 2*NC*lx*xi2*opxi +  
   4*z*lx*xi2*opxi +  
   4*NC*z*lx*xi2*opxi -  
   2*lomx*lx*xi2*opxi +  
   NC*lomx*lx*xi2*opxi -  
   2*z*lomx*lx*xi2*opxi -  
   2*NC*z*lomx*lx*xi2*opxi +  
   2*lx*lopx*xi2*opxi -  
   NC*lx*lopx*xi2*opxi +  
   2*z*lx*lopx*xi2*opxi +  
   2*NC*z*lx*lopx*xi2*opxi +  
   lx*lz*xi2*opxi -  
   (NC*lx*lz*xi2*opxi)/2. +  
   z*lx*lz*xi2*opxi +  
   NC*z*lx*lz*xi2*opxi -  
   lx*lopxz*xi2*opxi +  
   (NC*lx*lopxz*xi2*opxi)/2. -  
   z*lx*lopxz*xi2*opxi -  
   NC*z*lx*lopxz*xi2*opxi -  
   lz*lopxz*xi2*opxi +  
   (NC*lz*lopxz*xi2*opxi)/2. -  
   z*lz*lopxz*xi2*opxi -  
   NC*z*lz*lopxz*xi2*opxi -  
   lx*lopxzi*xi2*opxi +  
   (NC*lx*lopxzi*xi2*opxi)/2. -  
   z*lx*lopxzi*xi2*opxi -  
   NC*z*lx*lopxzi*xi2*opxi +  
   lz*lopxzi*xi2*opxi -  
   (NC*lz*lopxzi*xi2*opxi)/2. +  
   z*lz*lopxzi*xi2*opxi +  
   NC*z*lz*lopxzi*xi2*opxi -  
   2*Li2mx*NCi2*xi2*opxi -  
   2*z*Li2mx*NCi2*xi2*opxi +  
   2*Li2x*NCi2*xi2*opxi +  
   2*z*Li2x*NCi2*xi2*opxi +  
   li2spec19*NCi2*xi2*opxi +  
   z*li2spec19*NCi2*xi2*opxi +  
   li2spec20*NCi2*xi2*opxi +  
   z*li2spec20*NCi2*xi2*opxi -  
   4*lx*NCi2*xi2*opxi -  
   4*z*lx*NCi2*xi2*opxi +  
   2*lomx*lx*NCi2*xi2*opxi +  
   2*z*lomx*lx*NCi2*xi2*opxi -  
   2*lx*lopx*NCi2*xi2*opxi -  
   2*z*lx*lopx*NCi2*xi2*opxi -  
   lx*lz*NCi2*xi2*opxi -  
   z*lx*lz*NCi2*xi2*opxi +  
   lx*lopxz*NCi2*xi2*opxi +  
   z*lx*lopxz*NCi2*xi2*opxi +  
   lz*lopxz*NCi2*xi2*opxi +  
   z*lz*lopxz*NCi2*xi2*opxi +  
   lx*lopxzi*NCi2*xi2*opxi +  
   z*lx*lopxzi*NCi2*xi2*opxi -  
   lz*lopxzi*NCi2*xi2*opxi -  
   z*lz*lopxzi*NCi2*xi2*opxi +  
   Li2mx*NCi*xi2*opxi -  
   2*z*Li2mx*NCi*xi2*opxi -  
   Li2x*NCi*xi2*opxi +  
   2*z*Li2x*NCi*xi2*opxi -  
   (li2spec19*NCi*xi2*opxi)/2. +  
   z*li2spec19*NCi*xi2*opxi -  
   (li2spec20*NCi*xi2*opxi)/2. +  
   z*li2spec20*NCi*xi2*opxi +  
   2*lx*NCi*xi2*opxi -  
   4*z*lx*NCi*xi2*opxi -  
   lomx*lx*NCi*xi2*opxi +  
   2*z*lomx*lx*NCi*xi2*opxi +  
   lx*lopx*NCi*xi2*opxi -  
   2*z*lx*lopx*NCi*xi2*opxi +  
   (lx*lz*NCi*xi2*opxi)/2. -  
   z*lx*lz*NCi*xi2*opxi -  
   (lx*lopxz*NCi*xi2*opxi)/2. +  
   z*lx*lopxz*NCi*xi2*opxi -  
   (lz*lopxz*NCi*xi2*opxi)/2. +  
   z*lz*lopxz*NCi*xi2*opxi -  
   (lx*lopxzi*NCi*xi2*opxi)/2. +  
   z*lx*lopxzi*NCi*xi2*opxi +  
   (lz*lopxzi*NCi*xi2*opxi)/2. -  
   z*lz*lopxzi*NCi*xi2*opxi +  
   (pi2*xi2*opxi)/3. -  
   (NC*pi2*xi2*opxi)/6. +  
   (z*pi2*xi2*opxi)/3. +  
   (NC*z*pi2*xi2*opxi)/3. -  
   (NCi2*pi2*xi2*opxi)/3. -  
   (z*NCi2*pi2*xi2*opxi)/3. +  
   (NCi*pi2*xi2*opxi)/6. -  
   (z*NCi*pi2*xi2*opxi)/3. +  
   2*Li2mx*xi*opxi + 2*z*Li2mx*xi*opxi -  
   2*Li2x*xi*opxi - 2*z*Li2x*xi*opxi -  
   li2spec19*xi*opxi -  
   z*li2spec19*xi*opxi -  
   li2spec20*xi*opxi -  
   z*li2spec20*xi*opxi +  
   4*lx*xi*opxi + 4*z*lx*xi*opxi -  
   2*lomx*lx*xi*opxi -  
   2*z*lomx*lx*xi*opxi +  
   2*lx*lopx*xi*opxi +  
   2*z*lx*lopx*xi*opxi +  
   lx*lz*xi*opxi +  
   z*lx*lz*xi*opxi -  
   lx*lopxz*xi*opxi -  
   z*lx*lopxz*xi*opxi -  
   lz*lopxz*xi*opxi -  
   z*lz*lopxz*xi*opxi -  
   lx*lopxzi*xi*opxi -  
   z*lx*lopxzi*xi*opxi +  
   lz*lopxzi*xi*opxi +  
   z*lz*lopxzi*xi*opxi -  
   2*Li2mx*NCi2*xi*opxi -  
   2*z*Li2mx*NCi2*xi*opxi +  
   2*Li2x*NCi2*xi*opxi +  
   2*z*Li2x*NCi2*xi*opxi +  
   li2spec19*NCi2*xi*opxi +  
   z*li2spec19*NCi2*xi*opxi +  
   li2spec20*NCi2*xi*opxi +  
   z*li2spec20*NCi2*xi*opxi -  
   4*lx*NCi2*xi*opxi -  
   4*z*lx*NCi2*xi*opxi +  
   2*lomx*lx*NCi2*xi*opxi +  
   2*z*lomx*lx*NCi2*xi*opxi -  
   2*lx*lopx*NCi2*xi*opxi -  
   2*z*lx*lopx*NCi2*xi*opxi -  
   lx*lz*NCi2*xi*opxi -  
   z*lx*lz*NCi2*xi*opxi +  
   lx*lopxz*NCi2*xi*opxi +  
   z*lx*lopxz*NCi2*xi*opxi +  
   lz*lopxz*NCi2*xi*opxi +  
   z*lz*lopxz*NCi2*xi*opxi +  
   lx*lopxzi*NCi2*xi*opxi +  
   z*lx*lopxzi*NCi2*xi*opxi -  
   lz*lopxzi*NCi2*xi*opxi -  
   z*lz*lopxzi*NCi2*xi*opxi +  
   (pi2*xi*opxi)/3. +  
   (z*pi2*xi*opxi)/3. -  
   (NCi2*pi2*xi*opxi)/3. -  
   (z*NCi2*pi2*xi*opxi)/3. -  
   2*x*Li2x*omzi + 4*x*lx*omzi -  
   2*x*lomx*lx*omzi - (lz*omzi)/2. -  
   (NC*lz*omzi)/2. + (x*lz*omzi)/2. +  
   (NC*x*lz*omzi)/2. + (lx*lz*omzi)/2. +  
   (x*lx*lz*omzi)/2. + 2*x*Li2x*NCi2*omzi -  
   4*x*lx*NCi2*omzi +  
   2*x*lomx*lx*NCi2*omzi +  
   (lz*NCi2*omzi)/2. -  
   (x*lz*NCi2*omzi)/2. -  
   (lx*lz*NCi2*omzi)/2. -  
   (x*lx*lz*NCi2*omzi)/2. +  
   (lz*NCi*omzi)/2. -  
   (x*lz*NCi*omzi)/2. + (x*pi2*omzi)/3. -  
   (x*NCi2*pi2*omzi)/3. +  
   (NC*lz*xi*omzi)/3. -  
   (lz*NCi*xi*omzi)/3. -  
   (NC*lz*x2*omzi)/3. +  
   (lz*NCi*x2*omzi)/3. -  
   2*Li2x*opxi*omzi +  
   4*lx*opxi*omzi -  
   2*lomx*lx*opxi*omzi +  
   2*Li2x*NCi2*opxi*omzi -  
   4*lx*NCi2*opxi*omzi +  
   2*lomx*lx*NCi2*opxi*omzi +  
   (pi2*opxi*omzi)/3. -  
   (NCi2*pi2*opxi*omzi)/3. -  
   (NC*x*lx*xmzi)/2. + (NC*x*lz*xmzi)/2. +  
   (x*lx*NCi*xmzi)/2. -  
   (x*lz*NCi*xmzi)/2. +  
   NC*lx*x2*xmzi - NC*lz*x2*xmzi -  
   lx*NCi*x2*xmzi +  
   lz*NCi*x2*xmzi -  
   NC*lx*x3*xmzi + NC*lz*x3*xmzi +  
   lx*NCi*x3*xmzi -  
   lz*NCi*x3*xmzi + 2*zi -  
   (71*NC*zi)/72. - (7*x*zi)/4. + (121*NC*x*zi)/72. +  
   2*x*Li2mx*zi + (Li2x*zi)/2. - (NC*Li2x*zi)/2. -  
   (3*x*Li2x*zi)/2. - (NC*x*Li2x*zi)/2. -  
   2*x*li2spec1*zi -  
   2*sqrtxz1*x*li2spec1*zi +  
   2*x*li2spec2*zi -  
   2*sqrtxz1*x*li2spec2*zi +  
   4*sqrtxz1*x*li2spec21*zi +  
   2*sqrtxz1*x*li2spec22*zi 
    + 4*sqrtxz1*x*li2spec23*zi -  
   2*x*li2spec20*zi +  
   2*x*li2spec3*zi +  
   2*sqrtxz1*x*li2spec3*zi -  
   2*x*li2spec4*zi +  
   2*sqrtxz1*x*li2spec4*zi +  
   8*sqrtxz1*x*l2*zi + (NC*lomx*zi)/12. -  
   (5*NC*x*lomx*zi)/12. + (lx*zi)/2. +  
   (5*NC*lx*zi)/6. + (9*x*lx*zi)/2. +  
   (11*NC*x*lx*zi)/6. + 4*sqrtxz1*x*lx*zi +  
   2*x*l2*lx*zi + (lomx*lx*zi)/2. -  
   (3*x*lomx*lx*zi)/2. + 2*x*lx*lopx*zi +  
   (NC*lomz*zi)/12. - (5*NC*x*lomz*zi)/12. +  
   4*sqrtxz1*x*l2*lomz*zi +  
   (NC*lx*lomz*zi)/2. +  
   (NC*x*lx*lomz*zi)/2. +  
   2*sqrtxz1*x*lx*lomz*zi -  
   8*sqrtxz1*x*lspec1*zi -  
   4*x*l2*lspec1*zi -  
   8*sqrtxz1*x*l2*lspec1*zi +  
   2*sqrtxz1*x*lx*lspec1*zi -  
   4*sqrtxz1*x*lomz*lspec1*zi +  
   (NC*lz*zi)/12. - (5*NC*x*lz*zi)/12. +  
   4*sqrtxz1*x*lz*zi + 6*x*l2*lz*zi +  
   8*sqrtxz1*x*l2*lz*zi + (NC*lx*lz*zi)/2. +  
   2*x*lx*lz*zi + (NC*x*lx*lz*zi)/2. -  
   sqrtxz1*x*lx*lz*zi +  
   2*sqrtxz1*x*lomz*lz*zi -  
   2*x*lspec1*lz*zi -  
   4*sqrtxz1*x*lspec1*lz*zi -  
   4*x*l2*lspec2*zi -  
   4*sqrtxz1*x*l2*lspec2*zi -  
   2*x*lx*lspec2*zi -  
   2*sqrtxz1*x*lx*lspec2*zi +  
   4*x*lspec1*lspec2*zi +  
   4*sqrtxz1*x*lspec1*lspec2*zi -  
   4*x*lz*lspec2*zi -  
   4*sqrtxz1*x*lz*lspec2*zi -  
   x*lx*lxpz*zi + x*lz*lxpz*zi -  
   x*lx*lopxzi*zi +  
   x*lz*lopxzi*zi -  
   2*sqrtxz1*x*lomz*lspec7*zi -  
   2*NCi2*zi + (7*x*NCi2*zi)/4. -  
   2*x*Li2mx*NCi2*zi - (Li2x*NCi2*zi)/2. +  
   (3*x*Li2x*NCi2*zi)/2. +  
   2*x*li2spec1*NCi2*zi +  
   2*sqrtxz1*x*li2spec1*NCi2*zi -  
   2*x*li2spec2*NCi2*zi +  
   2*sqrtxz1*x*li2spec2*NCi2*zi -  
   4*sqrtxz1*x*li2spec21*NCi2*zi -  
   2*sqrtxz1*x*li2spec22*NCi2* 
    zi - 4*sqrtxz1*x*li2spec23*NCi2* 
    zi + 2*x*li2spec20*NCi2*zi -  
   2*x*li2spec3*NCi2* 
    zi - 2*sqrtxz1*x*li2spec3* 
    NCi2*zi + 2*x*li2spec4*NCi2*zi -  
   2*sqrtxz1*x*li2spec4*NCi2* 
    zi - 8*sqrtxz1*x*l2*NCi2*zi -  
   (lx*NCi2*zi)/2. - (9*x*lx*NCi2*zi)/2. -  
   4*sqrtxz1*x*lx*NCi2*zi -  
   2*x*l2*lx*NCi2*zi -  
   (lomx*lx*NCi2*zi)/2. +  
   (3*x*lomx*lx*NCi2*zi)/2. -  
   2*x*lx*lopx*NCi2*zi -  
   4*sqrtxz1*x*l2*lomz*NCi2*zi -  
   2*sqrtxz1*x*lx*lomz*NCi2*zi +  
   8*sqrtxz1*x*lspec1*NCi2*zi +  
   4*x*l2*lspec1*NCi2*zi +  
   8*sqrtxz1*x*l2*lspec1*NCi2*zi -  
   2*sqrtxz1*x*lx*lspec1*NCi2*zi +  
   4*sqrtxz1*x*lomz*lspec1*NCi2*zi -  
   4*sqrtxz1*x*lz*NCi2*zi -  
   6*x*l2*lz*NCi2*zi -  
   8*sqrtxz1*x*l2*lz*NCi2*zi -  
   2*x*lx*lz*NCi2*zi +  
   sqrtxz1*x*lx*lz*NCi2*zi -  
   2*sqrtxz1*x*lomz*lz*NCi2*zi +  
   2*x*lspec1*lz*NCi2*zi +  
   4*sqrtxz1*x*lspec1*lz*NCi2*zi +  
   4*x*l2*lspec2*NCi2*zi +  
   4*sqrtxz1*x*l2*lspec2*NCi2*zi +  
   2*x*lx*lspec2*NCi2*zi +  
   2*sqrtxz1*x*lx*lspec2*NCi2*zi -  
   4*x*lspec1*lspec2*NCi2*zi -  
   4*sqrtxz1*x*lspec1*lspec2*NCi2* 
    zi + 4*x*lz*lspec2*NCi2*zi +  
   4*sqrtxz1*x*lz*lspec2*NCi2*zi +  
   x*lx*lxpz*NCi2*zi -  
   x*lz*lxpz*NCi2*zi +  
   x*lx*lopxzi*NCi2*zi -  
   x*lz*lopxzi*NCi2*zi +  
   2*sqrtxz1*x*lomz*lspec7*NCi2* 
    zi + (71*NCi*zi)/72. -  
   (121*x*NCi*zi)/72. + (Li2x*NCi*zi)/2. +  
   (x*Li2x*NCi*zi)/2. -  
   (lomx*NCi*zi)/12. +  
   (5*x*lomx*NCi*zi)/12. -  
   (5*lx*NCi*zi)/6. -  
   (11*x*lx*NCi*zi)/6. -  
   (lomz*NCi*zi)/12. +  
   (5*x*lomz*NCi*zi)/12. -  
   (lx*lomz*NCi*zi)/2. -  
   (x*lx*lomz*NCi*zi)/2. -  
   (lz*NCi*zi)/12. +  
   (5*x*lz*NCi*zi)/12. -  
   (lx*lz*NCi*zi)/2. -  
   (x*lx*lz*NCi*zi)/2. - (pi2*zi)/12. +  
   (NC*pi2*zi)/12. + (x*pi2*zi)/4. +  
   (NC*x*pi2*zi)/12. + (sqrtxz1*x*pi2*zi)/3. +  
   (NCi2*pi2*zi)/12. -  
   (x*NCi2*pi2*zi)/4. -  
   (sqrtxz1*x*NCi2*pi2*zi)/3. -  
   (NCi*pi2*zi)/12. -  
   (x*NCi*pi2*zi)/12. - Li2x*omxi*zi +  
   2*li2spec1*omxi*zi +  
   2*sqrtxz1*li2spec1*omxi*zi -  
   2*li2spec2*omxi*zi +  
   2*sqrtxz1*li2spec2*omxi*zi -  
   li2spec19*omxi*zi -  
   4*sqrtxz1*li2spec21*omxi*zi -  
   2*sqrtxz1*li2spec22*omxi*zi -  
   4*sqrtxz1*li2spec23*omxi*zi +  
   li2spec20*omxi*zi -  
   2*li2spec3*omxi* 
    zi - 2*sqrtxz1*li2spec3* 
    omxi*zi + 2* 
    li2spec4*omxi*zi 
    - 2*sqrtxz1*li2spec4* 
    omxi*zi - 8*sqrtxz1*l2*omxi*zi +  
   (3*lx*omxi*zi)/4. -  
   (2*NC*lx*omxi*zi)/3. -  
   4*sqrtxz1*lx*omxi*zi -  
   2*l2*lx*omxi*zi -  
   lomx*lx*omxi*zi -  
   4*sqrtxz1*l2*lomz*omxi*zi -  
   2*sqrtxz1*lx*lomz*omxi*zi +  
   8*sqrtxz1*lspec1*omxi*zi +  
   4*l2*lspec1*omxi*zi +  
   8*sqrtxz1*l2*lspec1*omxi*zi -  
   2*sqrtxz1*lx*lspec1*omxi*zi +  
   4*sqrtxz1*lomz*lspec1*omxi*zi -  
   4*sqrtxz1*lz*omxi*zi -  
   6*l2*lz*omxi*zi -  
   8*sqrtxz1*l2*lz*omxi*zi -  
   lx*lz*omxi*zi +  
   sqrtxz1*lx*lz*omxi*zi -  
   2*sqrtxz1*lomz*lz*omxi*zi +  
   2*lspec1*lz*omxi*zi +  
   4*sqrtxz1*lspec1*lz*omxi*zi +  
   4*l2*lspec2*omxi*zi +  
   4*sqrtxz1*l2*lspec2*omxi*zi +  
   2*lx*lspec2*omxi*zi +  
   2*sqrtxz1*lx*lspec2*omxi*zi -  
   4*lspec1*lspec2*omxi*zi -  
   4*sqrtxz1*lspec1*lspec2*omxi* 
    zi + 4*lz*lspec2*omxi*zi +  
   4*sqrtxz1*lz*lspec2*omxi*zi +  
   lx*lxpz*omxi*zi -  
   lz*lxpz*omxi*zi -  
   lx*lopxz*omxi*zi -  
   lz*lopxz*omxi*zi +  
   2*sqrtxz1*lomz*lspec7*omxi* 
    zi + Li2x*NCi2*omxi*zi -  
   2*li2spec1*NCi2*omxi*zi -  
   2*sqrtxz1*li2spec1*NCi2*omxi* 
    zi + 2*li2spec2*NCi2*omxi* 
    zi - 2*sqrtxz1*li2spec2*NCi2* 
    omxi*zi + li2spec19*NCi2*omxi* 
    zi + 4*sqrtxz1*li2spec21*NCi2* 
    omxi*zi + 2*sqrtxz1* 
    li2spec22*NCi2*omxi*zi +  
   4*sqrtxz1*li2spec23*NCi2*omxi*zi -  
   li2spec20*NCi2*omxi*zi +  
   2*li2spec3*NCi2* 
    omxi*zi + 2*sqrtxz1* 
    li2spec3*NCi2*omxi* 
    zi - 2*li2spec4*NCi2* 
    omxi*zi + 2*sqrtxz1* 
    li2spec4*NCi2*omxi* 
    zi + 8*sqrtxz1*l2*NCi2*omxi*zi -  
   (3*lx*NCi2*omxi*zi)/4. +  
   4*sqrtxz1*lx*NCi2*omxi*zi +  
   2*l2*lx*NCi2*omxi*zi +  
   lomx*lx*NCi2*omxi*zi +  
   4*sqrtxz1*l2*lomz*NCi2*omxi*zi +  
   2*sqrtxz1*lx*lomz*NCi2*omxi*zi -  
   8*sqrtxz1*lspec1*NCi2*omxi*zi -  
   4*l2*lspec1*NCi2*omxi*zi -  
   8*sqrtxz1*l2*lspec1*NCi2*omxi*zi +  
   2*sqrtxz1*lx*lspec1*NCi2*omxi*zi -  
   4*sqrtxz1*lomz*lspec1*NCi2*omxi* 
    zi + 4*sqrtxz1*lz*NCi2*omxi*zi +  
   6*l2*lz*NCi2*omxi*zi +  
   8*sqrtxz1*l2*lz*NCi2*omxi*zi +  
   lx*lz*NCi2*omxi*zi -  
   sqrtxz1*lx*lz*NCi2*omxi*zi +  
   2*sqrtxz1*lomz*lz*NCi2*omxi*zi -  
   2*lspec1*lz*NCi2*omxi*zi -  
   4*sqrtxz1*lspec1*lz*NCi2*omxi*zi -  
   4*l2*lspec2*NCi2*omxi*zi -  
   4*sqrtxz1*l2*lspec2*NCi2*omxi*zi -  
   2*lx*lspec2*NCi2*omxi*zi -  
   2*sqrtxz1*lx*lspec2*NCi2*omxi*zi +  
   4*lspec1*lspec2*NCi2*omxi* 
    zi + 4*sqrtxz1*lspec1*lspec2* 
    NCi2*omxi*zi -  
   4*lz*lspec2*NCi2*omxi*zi -  
   4*sqrtxz1*lz*lspec2*NCi2*omxi*zi -  
   lx*lxpz*NCi2*omxi*zi +  
   lz*lxpz*NCi2*omxi*zi +  
   lx*lopxz*NCi2*omxi*zi +  
   lz*lopxz*NCi2*omxi*zi -  
   2*sqrtxz1*lomz*lspec7*NCi2* 
    omxi*zi + (2*lx*NCi*omxi*zi)/ 
    3. + (pi2*omxi*zi)/6. -  
   (sqrtxz1*pi2*omxi*zi)/3. -  
   (NCi2*pi2*omxi*zi)/6. +  
   (sqrtxz1*NCi2*pi2*omxi*zi)/3. -  
   2*Li2mx*xi2*zi + 2*Li2x*xi2*zi +  
   li2spec19*xi2*zi +  
   li2spec20*xi2*zi - 4*lx*xi2*zi +  
   2*lomx*lx*xi2*zi -  
   2*lx*lopx*xi2*zi -  
   lx*lz*xi2*zi +  
   lx*lopxz*xi2*zi +  
   lz*lopxz*xi2*zi +  
   lx*lopxzi*xi2*zi -  
   lz*lopxzi*xi2*zi +  
   2*Li2mx*NCi2*xi2*zi -  
   2*Li2x*NCi2*xi2*zi -  
   li2spec19*NCi2*xi2*zi -  
   li2spec20*NCi2*xi2*zi +  
   4*lx*NCi2*xi2*zi -  
   2*lomx*lx*NCi2*xi2*zi +  
   2*lx*lopx*NCi2*xi2*zi +  
   lx*lz*NCi2*xi2*zi -  
   lx*lopxz*NCi2*xi2*zi -  
   lz*lopxz*NCi2*xi2*zi -  
   lx*lopxzi*NCi2*xi2*zi +  
   lz*lopxzi*NCi2*xi2*zi -  
   (pi2*xi2*zi)/3. +  
   (NCi2*pi2*xi2*zi)/3. +  
   (13*NC*xi*zi)/18. + (NC*lomx*xi*zi)/3. +  
   (NC*lomz*xi*zi)/3. +  
   (NC*lz*xi*zi)/3. -  
   (13*NCi*xi*zi)/18. -  
   (lomx*NCi*xi*zi)/3. -  
   (lomz*NCi*xi*zi)/3. -  
   (lz*NCi*xi*zi)/3. -  
   (11*NC*x2*zi)/9. - (NC*lomx*x2*zi)/3. +  
   NC*lx*x2*zi - (NC*lomz*x2*zi)/3. -  
   (NC*lz*x2*zi)/3. +  
   (11*NCi*x2*zi)/9. +  
   (lomx*NCi*x2*zi)/3. -  
   lx*NCi*x2*zi +  
   (lomz*NCi*x2*zi)/3. +  
   (lz*NCi*x2*zi)/3. +  
   2*Li2mx*opxi*zi - 2*Li2x*opxi*zi -  
   li2spec19*opxi*zi -  
   li2spec20*opxi*zi +  
   4*lx*opxi*zi -  
   2*lomx*lx*opxi*zi +  
   2*lx*lopx*opxi*zi +  
   lx*lz*opxi*zi -  
   lx*lopxz*opxi*zi -  
   lz*lopxz*opxi*zi -  
   lx*lopxzi*opxi*zi +  
   lz*lopxzi*opxi*zi -  
   2*Li2mx*NCi2*opxi*zi +  
   2*Li2x*NCi2*opxi*zi +  
   li2spec19*NCi2*opxi*zi +  
   li2spec20*NCi2*opxi*zi -  
   4*lx*NCi2*opxi*zi +  
   2*lomx*lx*NCi2*opxi*zi -  
   2*lx*lopx*NCi2*opxi*zi -  
   lx*lz*NCi2*opxi*zi +  
   lx*lopxz*NCi2*opxi*zi +  
   lz*lopxz*NCi2*opxi*zi +  
   lx*lopxzi*NCi2*opxi*zi -  
   lz*lopxzi*NCi2*opxi*zi +  
   (pi2*opxi*zi)/3. -  
   (NCi2*pi2*opxi*zi)/3. +  
   2*Li2mx*xi2*opxi*zi -  
   2*Li2x*xi2*opxi*zi -  
   li2spec19*xi2*opxi*zi -  
   li2spec20*xi2*opxi*zi +  
   4*lx*xi2*opxi*zi -  
   2*lomx*lx*xi2*opxi*zi +  
   2*lx*lopx*xi2*opxi*zi +  
   lx*lz*xi2*opxi*zi -  
   lx*lopxz*xi2*opxi*zi -  
   lz*lopxz*xi2*opxi*zi -  
   lx*lopxzi*xi2*opxi*zi +  
   lz*lopxzi*xi2*opxi*zi -  
   2*Li2mx*NCi2*xi2*opxi*zi +  
   2*Li2x*NCi2*xi2*opxi*zi +  
   li2spec19*NCi2*xi2*opxi*zi +  
   li2spec20*NCi2*xi2*opxi*zi -  
   4*lx*NCi2*xi2*opxi*zi +  
   2*lomx*lx*NCi2*xi2*opxi*zi -  
   2*lx*lopx*NCi2*xi2*opxi*zi -  
   lx*lz*NCi2*xi2*opxi*zi +  
   lx*lopxz*NCi2*xi2*opxi*zi +  
   lz*lopxz*NCi2*xi2*opxi*zi +  
   lx*lopxzi*NCi2*xi2*opxi*zi -  
   lz*lopxzi*NCi2*xi2*opxi*zi +  
   (pi2*xi2*opxi*zi)/3. -  
   (NCi2*pi2*xi2*opxi*zi)/3. +  
   2*Li2mx*xi*opxi*zi -  
   2*Li2x*xi*opxi*zi -  
   li2spec19*xi*opxi*zi -  
   li2spec20*xi*opxi*zi +  
   4*lx*xi*opxi*zi -  
   2*lomx*lx*xi*opxi*zi +  
   2*lx*lopx*xi*opxi*zi +  
   lx*lz*xi*opxi*zi -  
   lx*lopxz*xi*opxi*zi -  
   lz*lopxz*xi*opxi*zi -  
   lx*lopxzi*xi*opxi*zi +  
   lz*lopxzi*xi*opxi*zi -  
   2*Li2mx*NCi2*xi*opxi*zi +  
   2*Li2x*NCi2*xi*opxi*zi +  
   li2spec19*NCi2*xi*opxi*zi +  
   li2spec20*NCi2*xi*opxi*zi -  
   4*lx*NCi2*xi*opxi*zi +  
   2*lomx*lx*NCi2*xi*opxi*zi -  
   2*lx*lopx*NCi2*xi*opxi*zi -  
   lx*lz*NCi2*xi*opxi*zi +  
   lx*lopxz*NCi2*xi*opxi*zi +  
   lz*lopxz*NCi2*xi*opxi*zi +  
   lx*lopxzi*NCi2*xi*opxi*zi -  
   lz*lopxzi*NCi2*xi*opxi*zi +  
   (pi2*xi*opxi*zi)/3. -  
   (NCi2*pi2*xi*opxi*zi)/3. -  
   2*x*Li2mx*omzi*zi + 2*x*Li2x*omzi*zi +  
   x*li2spec19*omzi*zi +  
   x*li2spec20*omzi*zi -  
   4*x*lx*omzi*zi +  
   2*x*lomx*lx*omzi*zi -  
   2*x*lx*lopx*omzi*zi -  
   x*lx*lz*omzi*zi +  
   x*lx*lopxz*omzi*zi +  
   x*lz*lopxz*omzi*zi +  
   x*lx*lopxzi*omzi*zi -  
   x*lz*lopxzi*omzi*zi +  
   2*x*Li2mx*NCi2*omzi*zi -  
   2*x*Li2x*NCi2*omzi*zi -  
   x*li2spec19*NCi2*omzi*zi -  
   x*li2spec20*NCi2*omzi*zi +  
   4*x*lx*NCi2*omzi*zi -  
   2*x*lomx*lx*NCi2*omzi*zi +  
   2*x*lx*lopx*NCi2*omzi*zi +  
   x*lx*lz*NCi2*omzi*zi -  
   x*lx*lopxz*NCi2*omzi*zi -  
   x*lz*lopxz*NCi2*omzi*zi -  
   x*lx*lopxzi*NCi2*omzi*zi +  
   x*lz*lopxzi*NCi2*omzi*zi -  
   (x*pi2*omzi*zi)/3. +  
   (x*NCi2*pi2*omzi*zi)/3. +  
   2*Li2mx*xi2*omzi*zi -  
   2*Li2x*xi2*omzi*zi -  
   li2spec19*xi2*omzi*zi -  
   li2spec20*xi2*omzi*zi +  
   4*lx*xi2*omzi*zi -  
   2*lomx*lx*xi2*omzi*zi +  
   2*lx*lopx*xi2*omzi*zi +  
   lx*lz*xi2*omzi*zi -  
   lx*lopxz*xi2*omzi*zi -  
   lz*lopxz*xi2*omzi*zi -  
   lx*lopxzi*xi2*omzi*zi +  
   lz*lopxzi*xi2*omzi*zi -  
   2*Li2mx*NCi2*xi2*omzi*zi +  
   2*Li2x*NCi2*xi2*omzi*zi +  
   li2spec19*NCi2*xi2*omzi*zi +  
   li2spec20*NCi2*xi2*omzi*zi -  
   4*lx*NCi2*xi2*omzi*zi +  
   2*lomx*lx*NCi2*xi2*omzi*zi -  
   2*lx*lopx*NCi2*xi2*omzi*zi -  
   lx*lz*NCi2*xi2*omzi*zi +  
   lx*lopxz*NCi2*xi2*omzi*zi +  
   lz*lopxz*NCi2*xi2*omzi*zi +  
   lx*lopxzi*NCi2*xi2*omzi*zi -  
   lz*lopxzi*NCi2*xi2*omzi*zi +  
   (pi2*xi2*omzi*zi)/3. -  
   (NCi2*pi2*xi2*omzi*zi)/3. -  
   2*Li2mx*opxi*omzi*zi +  
   2*Li2x*opxi*omzi*zi +  
   li2spec19*opxi*omzi*zi +  
   li2spec20*opxi*omzi*zi -  
   4*lx*opxi*omzi*zi +  
   2*lomx*lx*opxi*omzi*zi -  
   2*lx*lopx*opxi*omzi*zi -  
   lx*lz*opxi*omzi*zi +  
   lx*lopxz*opxi*omzi*zi +  
   lz*lopxz*opxi*omzi*zi +  
   lx*lopxzi*opxi*omzi*zi -  
   lz*lopxzi*opxi*omzi*zi +  
   2*Li2mx*NCi2*opxi*omzi*zi -  
   2*Li2x*NCi2*opxi*omzi*zi -  
   li2spec19*NCi2*opxi*omzi*zi -  
   li2spec20*NCi2*opxi*omzi*zi +  
   4*lx*NCi2*opxi*omzi*zi -  
   2*lomx*lx*NCi2*opxi*omzi*zi +  
   2*lx*lopx*NCi2*opxi*omzi*zi +  
   lx*lz*NCi2*opxi*omzi*zi -  
   lx*lopxz*NCi2*opxi*omzi*zi -  
   lz*lopxz*NCi2*opxi*omzi*zi -  
   lx*lopxzi*NCi2*opxi*omzi* 
    zi + lz*lopxzi*NCi2*opxi* 
    omzi*zi - (pi2*opxi*omzi* 
      zi)/3. + (NCi2*pi2*opxi*omzi* 
      zi)/3. - 2*Li2mx*xi2*opxi*omzi* 
    zi + 2*Li2x*xi2*opxi*omzi*zi +  
   li2spec19*xi2*opxi*omzi*zi +  
   li2spec20*xi2*opxi*omzi*zi -  
   4*lx*xi2*opxi*omzi*zi +  
   2*lomx*lx*xi2*opxi*omzi*zi -  
   2*lx*lopx*xi2*opxi*omzi*zi -  
   lx*lz*xi2*opxi*omzi*zi +  
   lx*lopxz*xi2*opxi*omzi*zi +  
   lz*lopxz*xi2*opxi*omzi*zi +  
   lx*lopxzi*xi2*opxi*omzi* 
    zi - lz*lopxzi*xi2*opxi* 
    omzi*zi + 2*Li2mx*NCi2*xi2*opxi* 
    omzi*zi - 2*Li2x*NCi2*xi2*opxi* 
    omzi*zi - li2spec19*NCi2*xi2*opxi* 
    omzi*zi - li2spec20*NCi2*xi2* 
    opxi*omzi*zi +  
   4*lx*NCi2*xi2*opxi*omzi*zi -  
   2*lomx*lx*NCi2*xi2*opxi*omzi* 
    zi + 2*lx*lopx*NCi2*xi2*opxi* 
    omzi*zi + lx*lz*NCi2*xi2* 
    opxi*omzi*zi -  
   lx*lopxz*NCi2*xi2*opxi*omzi* 
    zi - lz*lopxz*NCi2*xi2*opxi* 
    omzi*zi - lx*lopxzi*NCi2* 
    xi2*opxi*omzi*zi +  
   lz*lopxzi*NCi2*xi2*opxi* 
    omzi*zi - (pi2*xi2*opxi* 
      omzi*zi)/3. +  
   (NCi2*pi2*xi2*opxi*omzi*zi)/ 
    3. - 2*Li2mx*xi*opxi*omzi*zi +  
   2*Li2x*xi*opxi*omzi*zi +  
   li2spec19*xi*opxi*omzi*zi +  
   li2spec20*xi*opxi*omzi*zi -  
   4*lx*xi*opxi*omzi*zi +  
   2*lomx*lx*xi*opxi*omzi*zi -  
   2*lx*lopx*xi*opxi*omzi*zi -  
   lx*lz*xi*opxi*omzi*zi +  
   lx*lopxz*xi*opxi*omzi*zi +  
   lz*lopxz*xi*opxi*omzi*zi +  
   lx*lopxzi*xi*opxi*omzi* 
    zi - lz*lopxzi*xi*opxi* 
    omzi*zi + 2*Li2mx*NCi2*xi*opxi* 
    omzi*zi - 2*Li2x*NCi2*xi*opxi* 
    omzi*zi - li2spec19*NCi2*xi*opxi* 
    omzi*zi - li2spec20*NCi2*xi* 
    opxi*omzi*zi +  
   4*lx*NCi2*xi*opxi*omzi*zi -  
   2*lomx*lx*NCi2*xi*opxi*omzi* 
    zi + 2*lx*lopx*NCi2*xi*opxi* 
    omzi*zi + lx*lz*NCi2*xi* 
    opxi*omzi*zi -  
   lx*lopxz*NCi2*xi*opxi*omzi* 
    zi - lz*lopxz*NCi2*xi*opxi* 
    omzi*zi - lx*lopxzi*NCi2* 
    xi*opxi*omzi*zi +  
   lz*lopxzi*NCi2*xi*opxi* 
    omzi*zi - (pi2*xi*opxi* 
      omzi*zi)/3. +  
   (NCi2*pi2*xi*opxi*omzi*zi)/ 
    3. + (43*NC*z2)/36. - (77*NC*x*z2)/36. -  
   2*sqrtxz3*itani1*z2 -  
   (23*NC*sqrtxz3*itani1*z2)/4. +  
   2*sqrtxz3*itani2*z2 +  
   (23*NC*sqrtxz3*itani2*z2)/4. -  
   4*sqrtxz3*itani3*z2 -  
   (23*NC*sqrtxz3*itani3*z2)/2. +  
   4*NC*x*Li2x*z2 - 4*NC*x*li2spec1*z2 +  
   4*NC*x*li2spec2*z2 +  
   4*NC*x*li2spec19*z2 -  
   4*sqrtxz3*atanspec1*lspec3*z2 -  
   (23*NC*sqrtxz3*atanspec1*lspec3*z2)/2. -  
   (NC*lomx*z2)/3. + (2*NC*x*lomx*z2)/3. +  
   (2*NC*lx*z2)/3. - (4*NC*x*lx*z2)/3. +  
   8*NC*x*l2*lx*z2 + 4*NC*x*lomx*lx*z2 -  
   (NC*lomz*z2)/3. + (2*NC*x*lomz*z2)/3. -  
   12*NC*x*l2*lspec1*z2 -  
   4*NC*x*lx*lspec1*z2 - (NC*lz*z2)/3. +  
   (2*NC*x*lz*z2)/3. + 8*NC*x*l2*lz*z2 +  
   2*NC*x*lx*lz*z2 -  
   4*NC*x*lspec1*lz*z2 +  
   4*sqrtxz3*atanspec2*lspec4*z2 +  
   (23*NC*sqrtxz3*atanspec2*lspec4*z2)/2. -  
   4*NC*x*l2*lspec2*z2 -  
   4*NC*x*lx*lspec2*z2 +  
   4*NC*x*lspec1*lspec2*z2 -  
   4*NC*x*lz*lspec2*z2 -  
   2*NC*x*lx*lxpz*z2 + 2*NC*x*lz*lxpz*z2 +  
   4*NC*x*lx*lopxz*z2 +  
   4*NC*x*lz*lopxz*z2 +  
   2*NC*x*lx*lopxzi*z2 -  
   2*NC*x*lz*lopxzi*z2 +  
   2*sqrtxz3*itani1*NCi2*z2 -  
   2*sqrtxz3*itani2*NCi2*z2 +  
   4*sqrtxz3*itani3*NCi2*z2 +  
   4*sqrtxz3*atanspec1*lspec3*NCi2*z2 -  
   4*sqrtxz3*atanspec2*lspec4*NCi2*z2 -  
   (43*NCi*z2)/36. + (77*x*NCi*z2)/36. +  
   (23*sqrtxz3*itani1*NCi*z2)/4. -  
   (23*sqrtxz3*itani2*NCi*z2)/4. +  
   (23*sqrtxz3*itani3*NCi*z2)/2. -  
   4*x*Li2x*NCi*z2 +  
   4*x*li2spec1*NCi*z2 -  
   4*x*li2spec2*NCi*z2 -  
   4*x*li2spec19*NCi*z2 +  
   (23*sqrtxz3*atanspec1*lspec3*NCi*z2)/2. +  
   (lomx*NCi*z2)/3. -  
   (2*x*lomx*NCi*z2)/3. -  
   (2*lx*NCi*z2)/3. + (4*x*lx*NCi*z2)/3. -  
   8*x*l2*lx*NCi*z2 -  
   4*x*lomx*lx*NCi*z2 +  
   (lomz*NCi*z2)/3. -  
   (2*x*lomz*NCi*z2)/3. +  
   12*x*l2*lspec1*NCi*z2 +  
   4*x*lx*lspec1*NCi*z2 +  
   (lz*NCi*z2)/3. - (2*x*lz*NCi*z2)/3. -  
   8*x*l2*lz*NCi*z2 -  
   2*x*lx*lz*NCi*z2 +  
   4*x*lspec1*lz*NCi*z2 -  
   (23*sqrtxz3*atanspec2*lspec4*NCi*z2)/2. +  
   4*x*l2*lspec2*NCi*z2 +  
   4*x*lx*lspec2*NCi*z2 -  
   4*x*lspec1*lspec2*NCi*z2 +  
   4*x*lz*lspec2*NCi*z2 +  
   2*x*lx*lxpz*NCi*z2 -  
   2*x*lz*lxpz*NCi*z2 -  
   4*x*lx*lopxz*NCi*z2 -  
   4*x*lz*lopxz*NCi*z2 -  
   2*x*lx*lopxzi*NCi*z2 +  
   2*x*lz*lopxzi*NCi*z2 -  
   (2*NC*x*pi2*z2)/3. + (2*x*NCi*pi2*z2)/3. +  
   2*NC*x*li2spec5*sqrtxz2i*z2 -  
   2*NC*x*li2spec6*sqrtxz2i*z2 -  
   2*NC*x*li2spec7*sqrtxz2i* 
    z2 + 2*NC*x*li2spec8* 
    sqrtxz2i*z2 - 2*NC*x*lx*lspec5* 
    sqrtxz2i*z2 + 2*NC*x*lx*lspec6* 
    sqrtxz2i*z2 - 2*x*li2spec5*NCi* 
    sqrtxz2i*z2 + 2*x*li2spec6*NCi* 
    sqrtxz2i*z2 + 2*x* 
    li2spec7*NCi* 
    sqrtxz2i*z2 - 2*x* 
    li2spec8*NCi* 
    sqrtxz2i*z2 + 2*x*lx*lspec5*NCi* 
    sqrtxz2i*z2 - 2*x*lx*lspec6*NCi* 
    sqrtxz2i*z2 - 4*NC*Li2x*omxi*z2 +  
   4*NC*li2spec1*omxi*z2 -  
   4*NC*li2spec2*omxi*z2 -  
   2*NC*li2spec19*omxi*z2 +  
   2*NC*li2spec20*omxi*z2 +  
   (2*NC*lx*omxi*z2)/3. -  
   8*NC*l2*lx*omxi*z2 -  
   4*NC*lomx*lx*omxi*z2 +  
   12*NC*l2*lspec1*omxi*z2 +  
   4*NC*lx*lspec1*omxi*z2 -  
   8*NC*l2*lz*omxi*z2 -  
   4*NC*lx*lz*omxi*z2 +  
   4*NC*lspec1*lz*omxi*z2 +  
   4*NC*l2*lspec2*omxi*z2 +  
   4*NC*lx*lspec2*omxi*z2 -  
   4*NC*lspec1*lspec2*omxi*z2 +  
   4*NC*lz*lspec2*omxi*z2 +  
   2*NC*lx*lxpz*omxi*z2 -  
   2*NC*lz*lxpz*omxi*z2 -  
   2*NC*lx*lopxz*omxi*z2 -  
   2*NC*lz*lopxz*omxi*z2 +  
   4*Li2x*NCi*omxi*z2 -  
   4*li2spec1*NCi*omxi*z2 +  
   4*li2spec2*NCi*omxi*z2 +  
   2*li2spec19*NCi*omxi*z2 -  
   2*li2spec20*NCi*omxi*z2 -  
   (2*lx*NCi*omxi*z2)/3. +  
   8*l2*lx*NCi*omxi*z2 +  
   4*lomx*lx*NCi*omxi*z2 -  
   12*l2*lspec1*NCi*omxi*z2 -  
   4*lx*lspec1*NCi*omxi*z2 +  
   8*l2*lz*NCi*omxi*z2 +  
   4*lx*lz*NCi*omxi*z2 -  
   4*lspec1*lz*NCi*omxi*z2 -  
   4*l2*lspec2*NCi*omxi*z2 -  
   4*lx*lspec2*NCi*omxi*z2 +  
   4*lspec1*lspec2*NCi*omxi* 
    z2 - 4*lz*lspec2*NCi*omxi* 
    z2 - 2*lx*lxpz*NCi*omxi*z2 +  
   2*lz*lxpz*NCi*omxi*z2 +  
   2*lx*lopxz*NCi*omxi*z2 +  
   2*lz*lopxz*NCi*omxi*z2 +  
   NC*pi2*omxi*z2 -  
   NCi*pi2*omxi*z2 +  
   4*NC*Li2mx*xi2*z2 - 4*NC*Li2x*xi2*z2 -  
   2*NC*li2spec19*xi2*z2 -  
   2*NC*li2spec20*xi2*z2 +  
   8*NC*lx*xi2*z2 -  
   4*NC*lomx*lx*xi2*z2 +  
   4*NC*lx*lopx*xi2*z2 +  
   2*NC*lx*lz*xi2*z2 -  
   2*NC*lx*lopxz*xi2*z2 -  
   2*NC*lz*lopxz*xi2*z2 -  
   2*NC*lx*lopxzi*xi2*z2 +  
   2*NC*lz*lopxzi*xi2*z2 -  
   4*Li2mx*NCi*xi2*z2 +  
   4*Li2x*NCi*xi2*z2 +  
   2*li2spec19*NCi*xi2*z2 +  
   2*li2spec20*NCi*xi2*z2 -  
   8*lx*NCi*xi2*z2 +  
   4*lomx*lx*NCi*xi2*z2 -  
   4*lx*lopx*NCi*xi2*z2 -  
   2*lx*lz*NCi*xi2*z2 +  
   2*lx*lopxz*NCi*xi2*z2 +  
   2*lz*lopxz*NCi*xi2*z2 +  
   2*lx*lopxzi*NCi*xi2*z2 -  
   2*lz*lopxzi*NCi*xi2*z2 +  
   (2*NC*pi2*xi2*z2)/3. -  
   (2*NCi*pi2*xi2*z2)/3. +  
   2*NC*li2spec19*opxi*z2 +  
   2*NC*li2spec20*opxi*z2 -  
   2*NC*lx*lz*opxi*z2 +  
   2*NC*lx*lopxz*opxi*z2 +  
   2*NC*lz*lopxz*opxi*z2 +  
   2*NC*lx*lopxzi*opxi*z2 -  
   2*NC*lz*lopxzi*opxi*z2 -  
   2*li2spec19*NCi*opxi*z2 -  
   2*li2spec20*NCi*opxi*z2 +  
   2*lx*lz*NCi*opxi*z2 -  
   2*lx*lopxz*NCi*opxi*z2 -  
   2*lz*lopxz*NCi*opxi*z2 -  
   2*lx*lopxzi*NCi*opxi*z2 +  
   2*lz*lopxzi*NCi*opxi*z2 +  
   (NC*pi2*opxi*z2)/3. -  
   (NCi*pi2*opxi*z2)/3. -  
   4*NC*Li2mx*xi2*opxi*z2 +  
   4*NC*Li2x*xi2*opxi*z2 +  
   2*NC*li2spec19*xi2*opxi*z2 +  
   2*NC*li2spec20*xi2*opxi*z2 -  
   8*NC*lx*xi2*opxi*z2 +  
   4*NC*lomx*lx*xi2*opxi*z2 -  
   4*NC*lx*lopx*xi2*opxi*z2 -  
   2*NC*lx*lz*xi2*opxi*z2 +  
   2*NC*lx*lopxz*xi2*opxi*z2 +  
   2*NC*lz*lopxz*xi2*opxi*z2 +  
   2*NC*lx*lopxzi*xi2*opxi*z2 -  
   2*NC*lz*lopxzi*xi2*opxi*z2 +  
   4*Li2mx*NCi*xi2*opxi*z2 -  
   4*Li2x*NCi*xi2*opxi*z2 -  
   2*li2spec19*NCi*xi2*opxi*z2 -  
   2*li2spec20*NCi*xi2*opxi*z2 +  
   8*lx*NCi*xi2*opxi*z2 -  
   4*lomx*lx*NCi*xi2*opxi*z2 +  
   4*lx*lopx*NCi*xi2*opxi*z2 +  
   2*lx*lz*NCi*xi2*opxi*z2 -  
   2*lx*lopxz*NCi*xi2*opxi*z2 -  
   2*lz*lopxz*NCi*xi2*opxi*z2 -  
   2*lx*lopxzi*NCi*xi2*opxi* 
    z2 + 2*lz*lopxzi*NCi*xi2* 
    opxi*z2 - (2*NC*pi2*xi2*opxi* 
      z2)/3. + (2*NCi*pi2*xi2*opxi* 
      z2)/3. - 4*NC*Li2mx*xi*opxi*z2 +  
   4*NC*Li2x*xi*opxi*z2 +  
   2*NC*li2spec19*xi*opxi*z2 +  
   2*NC*li2spec20*xi*opxi*z2 -  
   8*NC*lx*xi*opxi*z2 +  
   4*NC*lomx*lx*xi*opxi*z2 -  
   4*NC*lx*lopx*xi*opxi*z2 -  
   2*NC*lx*lz*xi*opxi*z2 +  
   2*NC*lx*lopxz*xi*opxi*z2 +  
   2*NC*lz*lopxz*xi*opxi*z2 +  
   2*NC*lx*lopxzi*xi*opxi*z2 -  
   2*NC*lz*lopxzi*xi*opxi*z2 +  
   4*Li2mx*NCi*xi*opxi*z2 -  
   4*Li2x*NCi*xi*opxi*z2 -  
   2*li2spec19*NCi*xi*opxi*z2 -  
   2*li2spec20*NCi*xi*opxi*z2 +  
   8*lx*NCi*xi*opxi*z2 -  
   4*lomx*lx*NCi*xi*opxi*z2 +  
   4*lx*lopx*NCi*xi*opxi*z2 +  
   2*lx*lz*NCi*xi*opxi*z2 -  
   2*lx*lopxz*NCi*xi*opxi*z2 -  
   2*lz*lopxz*NCi*xi*opxi*z2 -  
   2*lx*lopxzi*NCi*xi*opxi* 
    z2 + 2*lz*lopxzi*NCi*xi* 
    opxi*z2 - (2*NC*pi2*xi*opxi* 
      z2)/3. + (2*NCi*pi2*xi*opxi* 
      z2)/3. - 8*li2spec1*sqrtxz1i* 
    opzi + 6*x*li2spec1*sqrtxz1i* 
    opzi - 2*x*z*li2spec1*sqrtxz1i* 
    opzi - 8*li2spec2*sqrtxz1i* 
    opzi + 6*x*li2spec2*sqrtxz1i* 
    opzi - 2*x*z*li2spec2*sqrtxz1i* 
    opzi + 16*li2spec21*sqrtxz1i* 
    opzi - 12*x*li2spec21*sqrtxz1i* 
    opzi + 4*x*z*li2spec21*sqrtxz1i* 
    opzi + 8*li2spec22* 
    sqrtxz1i*opzi -  
   6*x*li2spec22*sqrtxz1i*opzi +  
   2*x*z*li2spec22*sqrtxz1i*opzi +  
   16*li2spec23*sqrtxz1i*opzi -  
   12*x*li2spec23*sqrtxz1i*opzi +  
   4*x*z*li2spec23*sqrtxz1i*opzi +  
   8*li2spec3*sqrtxz1i* 
    opzi - 6*x*li2spec3* 
    sqrtxz1i*opzi +  
   2*x*z*li2spec3*sqrtxz1i* 
    opzi + 8*li2spec4* 
    sqrtxz1i*opzi -  
   6*x*li2spec4*sqrtxz1i* 
    opzi + 2*x*z*li2spec4* 
    sqrtxz1i*opzi + 32*l2*sqrtxz1i*opzi -  
   24*x*l2*sqrtxz1i*opzi +  
   8*x*z*l2*sqrtxz1i*opzi +  
   16*lx*sqrtxz1i*opzi -  
   12*x*lx*sqrtxz1i*opzi +  
   4*x*z*lx*sqrtxz1i*opzi +  
   16*l2*lomz*sqrtxz1i*opzi -  
   12*x*l2*lomz*sqrtxz1i*opzi +  
   4*x*z*l2*lomz*sqrtxz1i*opzi +  
   8*lx*lomz*sqrtxz1i*opzi -  
   6*x*lx*lomz*sqrtxz1i*opzi +  
   2*x*z*lx*lomz*sqrtxz1i*opzi -  
   32*lspec1*sqrtxz1i*opzi +  
   24*x*lspec1*sqrtxz1i*opzi -  
   8*x*z*lspec1*sqrtxz1i*opzi -  
   32*l2*lspec1*sqrtxz1i*opzi +  
   24*x*l2*lspec1*sqrtxz1i*opzi -  
   8*x*z*l2*lspec1*sqrtxz1i*opzi +  
   8*lx*lspec1*sqrtxz1i*opzi -  
   6*x*lx*lspec1*sqrtxz1i*opzi +  
   2*x*z*lx*lspec1*sqrtxz1i*opzi -  
   16*lomz*lspec1*sqrtxz1i*opzi +  
   12*x*lomz*lspec1*sqrtxz1i*opzi -  
   4*x*z*lomz*lspec1*sqrtxz1i*opzi +  
   16*lz*sqrtxz1i*opzi -  
   12*x*lz*sqrtxz1i*opzi +  
   4*x*z*lz*sqrtxz1i*opzi +  
   32*l2*lz*sqrtxz1i*opzi -  
   24*x*l2*lz*sqrtxz1i*opzi +  
   8*x*z*l2*lz*sqrtxz1i*opzi -  
   4*lx*lz*sqrtxz1i*opzi +  
   3*x*lx*lz*sqrtxz1i*opzi -  
   x*z*lx*lz*sqrtxz1i*opzi +  
   8*lomz*lz*sqrtxz1i*opzi -  
   6*x*lomz*lz*sqrtxz1i*opzi +  
   2*x*z*lomz*lz*sqrtxz1i*opzi -  
   16*lspec1*lz*sqrtxz1i*opzi +  
   12*x*lspec1*lz*sqrtxz1i*opzi -  
   4*x*z*lspec1*lz*sqrtxz1i*opzi -  
   16*l2*lspec2*sqrtxz1i*opzi +  
   12*x*l2*lspec2*sqrtxz1i*opzi -  
   4*x*z*l2*lspec2*sqrtxz1i*opzi -  
   8*lx*lspec2*sqrtxz1i*opzi +  
   6*x*lx*lspec2*sqrtxz1i*opzi -  
   2*x*z*lx*lspec2*sqrtxz1i*opzi +  
   16*lspec1*lspec2*sqrtxz1i* 
    opzi - 12*x*lspec1*lspec2* 
    sqrtxz1i*opzi +  
   4*x*z*lspec1*lspec2*sqrtxz1i* 
    opzi - 16*lz*lspec2*sqrtxz1i* 
    opzi + 12*x*lz*lspec2*sqrtxz1i* 
    opzi - 4*x*z*lz*lspec2*sqrtxz1i* 
    opzi - 8*lomz*lspec7* 
    sqrtxz1i*opzi +  
   6*x*lomz*lspec7*sqrtxz1i* 
    opzi - 2*x*z*lomz*lspec7* 
    sqrtxz1i*opzi +  
   8*li2spec1*NCi2*sqrtxz1i*opzi -  
   6*x*li2spec1*NCi2*sqrtxz1i* 
    opzi + 2*x*z*li2spec1*NCi2* 
    sqrtxz1i*opzi +  
   8*li2spec2*NCi2*sqrtxz1i*opzi -  
   6*x*li2spec2*NCi2*sqrtxz1i* 
    opzi + 2*x*z*li2spec2*NCi2* 
    sqrtxz1i*opzi -  
   16*li2spec21*NCi2*sqrtxz1i*opzi +  
   12*x*li2spec21*NCi2*sqrtxz1i* 
    opzi - 4*x*z*li2spec21*NCi2* 
    sqrtxz1i*opzi -  
   8*li2spec22*NCi2*sqrtxz1i*opzi +  
   6*x*li2spec22*NCi2*sqrtxz1i*opzi -  
   2*x*z*li2spec22*NCi2*sqrtxz1i*opzi -  
   16*li2spec23*NCi2*sqrtxz1i*opzi +  
   12*x*li2spec23*NCi2*sqrtxz1i* 
    opzi - 4*x*z*li2spec23*NCi2* 
    sqrtxz1i*opzi -  
   8*li2spec3*NCi2* 
    sqrtxz1i*opzi +  
   6*x*li2spec3*NCi2* 
    sqrtxz1i*opzi -  
   2*x*z*li2spec3*NCi2* 
    sqrtxz1i*opzi -  
   8*li2spec4*NCi2* 
    sqrtxz1i*opzi +  
   6*x*li2spec4*NCi2* 
    sqrtxz1i*opzi -  
   2*x*z*li2spec4*NCi2* 
    sqrtxz1i*opzi -  
   32*l2*NCi2*sqrtxz1i*opzi +  
   24*x*l2*NCi2*sqrtxz1i*opzi -  
   8*x*z*l2*NCi2*sqrtxz1i*opzi -  
   16*lx*NCi2*sqrtxz1i*opzi +  
   12*x*lx*NCi2*sqrtxz1i*opzi -  
   4*x*z*lx*NCi2*sqrtxz1i*opzi -  
   16*l2*lomz*NCi2*sqrtxz1i*opzi +  
   12*x*l2*lomz*NCi2*sqrtxz1i*opzi -  
   4*x*z*l2*lomz*NCi2*sqrtxz1i*opzi -  
   8*lx*lomz*NCi2*sqrtxz1i*opzi +  
   6*x*lx*lomz*NCi2*sqrtxz1i*opzi -  
   2*x*z*lx*lomz*NCi2*sqrtxz1i*opzi +  
   32*lspec1*NCi2*sqrtxz1i*opzi -  
   24*x*lspec1*NCi2*sqrtxz1i*opzi +  
   8*x*z*lspec1*NCi2*sqrtxz1i*opzi +  
   32*l2*lspec1*NCi2*sqrtxz1i*opzi -  
   24*x*l2*lspec1*NCi2*sqrtxz1i* 
    opzi + 8*x*z*l2*lspec1*NCi2* 
    sqrtxz1i*opzi -  
   8*lx*lspec1*NCi2*sqrtxz1i*opzi +  
   6*x*lx*lspec1*NCi2*sqrtxz1i*opzi -  
   2*x*z*lx*lspec1*NCi2*sqrtxz1i* 
    opzi + 16*lomz*lspec1*NCi2* 
    sqrtxz1i*opzi -  
   12*x*lomz*lspec1*NCi2*sqrtxz1i* 
    opzi + 4*x*z*lomz*lspec1*NCi2* 
    sqrtxz1i*opzi -  
   16*lz*NCi2*sqrtxz1i*opzi +  
   12*x*lz*NCi2*sqrtxz1i*opzi -  
   4*x*z*lz*NCi2*sqrtxz1i*opzi -  
   32*l2*lz*NCi2*sqrtxz1i*opzi +  
   24*x*l2*lz*NCi2*sqrtxz1i*opzi -  
   8*x*z*l2*lz*NCi2*sqrtxz1i*opzi +  
   4*lx*lz*NCi2*sqrtxz1i*opzi -  
   3*x*lx*lz*NCi2*sqrtxz1i*opzi +  
   x*z*lx*lz*NCi2*sqrtxz1i*opzi -  
   8*lomz*lz*NCi2*sqrtxz1i*opzi +  
   6*x*lomz*lz*NCi2*sqrtxz1i*opzi -  
   2*x*z*lomz*lz*NCi2*sqrtxz1i*opzi +  
   16*lspec1*lz*NCi2*sqrtxz1i*opzi -  
   12*x*lspec1*lz*NCi2*sqrtxz1i* 
    opzi + 4*x*z*lspec1*lz*NCi2* 
    sqrtxz1i*opzi +  
   16*l2*lspec2*NCi2*sqrtxz1i*opzi -  
   12*x*l2*lspec2*NCi2*sqrtxz1i* 
    opzi + 4*x*z*l2*lspec2*NCi2* 
    sqrtxz1i*opzi +  
   8*lx*lspec2*NCi2*sqrtxz1i*opzi -  
   6*x*lx*lspec2*NCi2*sqrtxz1i*opzi +  
   2*x*z*lx*lspec2*NCi2*sqrtxz1i* 
    opzi - 16*lspec1*lspec2*NCi2* 
    sqrtxz1i*opzi +  
   12*x*lspec1*lspec2*NCi2*sqrtxz1i* 
    opzi - 4*x*z*lspec1*lspec2* 
    NCi2*sqrtxz1i*opzi +  
   16*lz*lspec2*NCi2*sqrtxz1i*opzi -  
   12*x*lz*lspec2*NCi2*sqrtxz1i* 
    opzi + 4*x*z*lz*lspec2*NCi2* 
    sqrtxz1i*opzi +  
   8*lomz*lspec7*NCi2*sqrtxz1i* 
    opzi - 6*x*lomz*lspec7*NCi2* 
    sqrtxz1i*opzi +  
   2*x*z*lomz*lspec7*NCi2* 
    sqrtxz1i*opzi +  
   (4*pi2*sqrtxz1i*opzi)/3. -  
   x*pi2*sqrtxz1i*opzi +  
   (x*z*pi2*sqrtxz1i*opzi)/3. -  
   (4*NCi2*pi2*sqrtxz1i*opzi)/3. +  
   x*NCi2*pi2*sqrtxz1i*opzi -  
   (x*z*NCi2*pi2*sqrtxz1i*opzi)/3. -  
   2*Li2mz*omxi*opzi -  
   2*lz*lopz*omxi*opzi +  
   2*Li2mz*NCi2*omxi*opzi +  
   2*lz*lopz*NCi2*omxi*opzi -  
   (pi2*omxi*opzi)/6. +  
   (NCi2*pi2*omxi*opzi)/6. +  
   2*li2spec1*sqrtxz1i*omxi* 
    opzi + 2*z*li2spec1*sqrtxz1i* 
    omxi*opzi +  
   2*li2spec2*sqrtxz1i*omxi* 
    opzi + 2*z*li2spec2*sqrtxz1i* 
    omxi*opzi -  
   4*li2spec21*sqrtxz1i*omxi* 
    opzi - 4*z*li2spec21*sqrtxz1i* 
    omxi*opzi -  
   2*li2spec22*sqrtxz1i*omxi*opzi 
    - 2*z*li2spec22*sqrtxz1i*omxi*opzi 
    - 4*li2spec23*sqrtxz1i*omxi* 
    opzi - 4*z*li2spec23*sqrtxz1i* 
    omxi*opzi -  
   2*li2spec3*sqrtxz1i* 
    omxi*opzi -  
   2*z*li2spec3*sqrtxz1i* 
    omxi*opzi -  
   2*li2spec4*sqrtxz1i* 
    omxi*opzi -  
   2*z*li2spec4*sqrtxz1i* 
    omxi*opzi -  
   8*l2*sqrtxz1i*omxi*opzi -  
   8*z*l2*sqrtxz1i*omxi*opzi -  
   4*lx*sqrtxz1i*omxi*opzi -  
   4*z*lx*sqrtxz1i*omxi*opzi -  
   4*l2*lomz*sqrtxz1i*omxi*opzi -  
   4*z*l2*lomz*sqrtxz1i*omxi*opzi -  
   2*lx*lomz*sqrtxz1i*omxi*opzi -  
   2*z*lx*lomz*sqrtxz1i*omxi*opzi +  
   8*lspec1*sqrtxz1i*omxi*opzi +  
   8*z*lspec1*sqrtxz1i*omxi*opzi +  
   8*l2*lspec1*sqrtxz1i*omxi* 
    opzi + 8*z*l2*lspec1*sqrtxz1i* 
    omxi*opzi -  
   2*lx*lspec1*sqrtxz1i*omxi* 
    opzi - 2*z*lx*lspec1*sqrtxz1i* 
    omxi*opzi +  
   4*lomz*lspec1*sqrtxz1i*omxi* 
    opzi + 4*z*lomz*lspec1*sqrtxz1i* 
    omxi*opzi -  
   4*lz*sqrtxz1i*omxi*opzi -  
   4*z*lz*sqrtxz1i*omxi*opzi -  
   8*l2*lz*sqrtxz1i*omxi*opzi -  
   8*z*l2*lz*sqrtxz1i*omxi*opzi +  
   lx*lz*sqrtxz1i*omxi*opzi +  
   z*lx*lz*sqrtxz1i*omxi*opzi -  
   2*lomz*lz*sqrtxz1i*omxi*opzi -  
   2*z*lomz*lz*sqrtxz1i*omxi*opzi +  
   4*lspec1*lz*sqrtxz1i*omxi* 
    opzi + 4*z*lspec1*lz*sqrtxz1i* 
    omxi*opzi +  
   4*l2*lspec2*sqrtxz1i*omxi* 
    opzi + 4*z*l2*lspec2*sqrtxz1i* 
    omxi*opzi +  
   2*lx*lspec2*sqrtxz1i*omxi* 
    opzi + 2*z*lx*lspec2*sqrtxz1i* 
    omxi*opzi -  
   4*lspec1*lspec2*sqrtxz1i*omxi* 
    opzi - 4*z*lspec1*lspec2* 
    sqrtxz1i*omxi*opzi +  
   4*lz*lspec2*sqrtxz1i*omxi* 
    opzi + 4*z*lz*lspec2*sqrtxz1i* 
    omxi*opzi +  
   2*lomz*lspec7*sqrtxz1i*omxi* 
    opzi + 2*z*lomz*lspec7* 
    sqrtxz1i*omxi*opzi -  
   2*li2spec1*NCi2*sqrtxz1i*omxi* 
    opzi - 2*z*li2spec1*NCi2* 
    sqrtxz1i*omxi*opzi -  
   2*li2spec2*NCi2*sqrtxz1i*omxi* 
    opzi - 2*z*li2spec2*NCi2* 
    sqrtxz1i*omxi*opzi +  
   4*li2spec21*NCi2*sqrtxz1i*omxi* 
    opzi + 4*z*li2spec21*NCi2* 
    sqrtxz1i*omxi*opzi +  
   2*li2spec22*NCi2*sqrtxz1i*omxi* 
    opzi + 2*z*li2spec22*NCi2* 
    sqrtxz1i*omxi*opzi +  
   4*li2spec23*NCi2*sqrtxz1i*omxi* 
    opzi + 4*z*li2spec23*NCi2* 
    sqrtxz1i*omxi*opzi +  
   2*li2spec3*NCi2* 
    sqrtxz1i*omxi*opzi +  
   2*z*li2spec3*NCi2* 
    sqrtxz1i*omxi*opzi +  
   2*li2spec4*NCi2* 
    sqrtxz1i*omxi*opzi +  
   2*z*li2spec4*NCi2* 
    sqrtxz1i*omxi*opzi +  
   8*l2*NCi2*sqrtxz1i*omxi*opzi +  
   8*z*l2*NCi2*sqrtxz1i*omxi*opzi +  
   4*lx*NCi2*sqrtxz1i*omxi*opzi +  
   4*z*lx*NCi2*sqrtxz1i*omxi*opzi +  
   4*l2*lomz*NCi2*sqrtxz1i*omxi* 
    opzi + 4*z*l2*lomz*NCi2*sqrtxz1i* 
    omxi*opzi +  
   2*lx*lomz*NCi2*sqrtxz1i*omxi* 
    opzi + 2*z*lx*lomz*NCi2*sqrtxz1i* 
    omxi*opzi -  
   8*lspec1*NCi2*sqrtxz1i*omxi* 
    opzi - 8*z*lspec1*NCi2*sqrtxz1i* 
    omxi*opzi -  
   8*l2*lspec1*NCi2*sqrtxz1i*omxi* 
    opzi - 8*z*l2*lspec1*NCi2* 
    sqrtxz1i*omxi*opzi +  
   2*lx*lspec1*NCi2*sqrtxz1i*omxi* 
    opzi + 2*z*lx*lspec1*NCi2* 
    sqrtxz1i*omxi*opzi -  
   4*lomz*lspec1*NCi2*sqrtxz1i*omxi* 
    opzi - 4*z*lomz*lspec1*NCi2* 
    sqrtxz1i*omxi*opzi +  
   4*lz*NCi2*sqrtxz1i*omxi*opzi +  
   4*z*lz*NCi2*sqrtxz1i*omxi*opzi +  
   8*l2*lz*NCi2*sqrtxz1i*omxi*opzi +  
   8*z*l2*lz*NCi2*sqrtxz1i*omxi*opzi -  
   lx*lz*NCi2*sqrtxz1i*omxi*opzi -  
   z*lx*lz*NCi2*sqrtxz1i*omxi*opzi +  
   2*lomz*lz*NCi2*sqrtxz1i*omxi* 
    opzi + 2*z*lomz*lz*NCi2*sqrtxz1i* 
    omxi*opzi -  
   4*lspec1*lz*NCi2*sqrtxz1i*omxi* 
    opzi - 4*z*lspec1*lz*NCi2* 
    sqrtxz1i*omxi*opzi -  
   4*l2*lspec2*NCi2*sqrtxz1i*omxi* 
    opzi - 4*z*l2*lspec2*NCi2* 
    sqrtxz1i*omxi*opzi -  
   2*lx*lspec2*NCi2*sqrtxz1i*omxi* 
    opzi - 2*z*lx*lspec2*NCi2* 
    sqrtxz1i*omxi*opzi +  
   4*lspec1*lspec2*NCi2*sqrtxz1i* 
    omxi*opzi +  
   4*z*lspec1*lspec2*NCi2*sqrtxz1i* 
    omxi*opzi -  
   4*lz*lspec2*NCi2*sqrtxz1i*omxi* 
    opzi - 4*z*lz*lspec2*NCi2* 
    sqrtxz1i*omxi*opzi -  
   2*lomz*lspec7*NCi2*sqrtxz1i* 
    omxi*opzi -  
   2*z*lomz*lspec7*NCi2*sqrtxz1i* 
    omxi*opzi -  
   (pi2*sqrtxz1i*omxi*opzi)/3. -  
   (z*pi2*sqrtxz1i*omxi*opzi)/3. +  
   (NCi2*pi2*sqrtxz1i*omxi*opzi)/3. +  
   (z*NCi2*pi2*sqrtxz1i*omxi*opzi)/3. -  
   8*li2spec1*sqrtxz1i*x2*opzi -  
   8*li2spec2*sqrtxz1i*x2*opzi +  
   16*li2spec21*sqrtxz1i*x2*opzi +  
   8*li2spec22*sqrtxz1i*x2*opzi +  
   16*li2spec23*sqrtxz1i*x2*opzi +  
   8*li2spec3*sqrtxz1i* 
    x2*opzi + 8*li2spec4*sqrtxz1i*x2*opzi +  
   32*l2*sqrtxz1i*x2*opzi +  
   16*lx*sqrtxz1i*x2*opzi +  
   16*l2*lomz*sqrtxz1i*x2*opzi +  
   8*lx*lomz*sqrtxz1i*x2*opzi -  
   32*lspec1*sqrtxz1i*x2*opzi -  
   32*l2*lspec1*sqrtxz1i*x2*opzi +  
   8*lx*lspec1*sqrtxz1i*x2*opzi -  
   16*lomz*lspec1*sqrtxz1i*x2* 
    opzi + 16*lz*sqrtxz1i*x2*opzi +  
   32*l2*lz*sqrtxz1i*x2*opzi -  
   4*lx*lz*sqrtxz1i*x2*opzi +  
   8*lomz*lz*sqrtxz1i*x2*opzi -  
   16*lspec1*lz*sqrtxz1i*x2*opzi -  
   16*l2*lspec2*sqrtxz1i*x2*opzi -  
   8*lx*lspec2*sqrtxz1i*x2*opzi +  
   16*lspec1*lspec2*sqrtxz1i*x2* 
    opzi - 16*lz*lspec2*sqrtxz1i*x2* 
    opzi - 8*lomz*lspec7* 
    sqrtxz1i*x2*opzi +  
   8*li2spec1*NCi2*sqrtxz1i*x2* 
    opzi + 8*li2spec2*NCi2*sqrtxz1i* 
    x2*opzi - 16*li2spec21*NCi2* 
    sqrtxz1i*x2*opzi -  
   8*li2spec22*NCi2*sqrtxz1i*x2* 
    opzi - 16*li2spec23*NCi2*sqrtxz1i* 
    x2*opzi - 8*li2spec3*NCi2*sqrtxz1i*x2* 
    opzi - 8*li2spec4* 
    NCi2*sqrtxz1i*x2*opzi -  
   32*l2*NCi2*sqrtxz1i*x2*opzi -  
   16*lx*NCi2*sqrtxz1i*x2*opzi -  
   16*l2*lomz*NCi2*sqrtxz1i*x2*opzi -  
   8*lx*lomz*NCi2*sqrtxz1i*x2*opzi +  
   32*lspec1*NCi2*sqrtxz1i*x2* 
    opzi + 32*l2*lspec1*NCi2*sqrtxz1i* 
    x2*opzi - 8*lx*lspec1*NCi2* 
    sqrtxz1i*x2*opzi +  
   16*lomz*lspec1*NCi2*sqrtxz1i*x2* 
    opzi - 16*lz*NCi2*sqrtxz1i*x2* 
    opzi - 32*l2*lz*NCi2*sqrtxz1i*x2* 
    opzi + 4*lx*lz*NCi2*sqrtxz1i*x2* 
    opzi - 8*lomz*lz*NCi2*sqrtxz1i*x2* 
    opzi + 16*lspec1*lz*NCi2*sqrtxz1i* 
    x2*opzi + 16*l2*lspec2*NCi2* 
    sqrtxz1i*x2*opzi +  
   8*lx*lspec2*NCi2*sqrtxz1i*x2* 
    opzi - 16*lspec1*lspec2*NCi2* 
    sqrtxz1i*x2*opzi +  
   16*lz*lspec2*NCi2*sqrtxz1i*x2* 
    opzi + 8*lomz*lspec7*NCi2* 
    sqrtxz1i*x2*opzi +  
   (4*pi2*sqrtxz1i*x2*opzi)/3. -  
   (4*NCi2*pi2*sqrtxz1i*x2*opzi)/3. +  
   2*x*li2spec1*zi*opzi +  
   2*sqrtxz1*x*li2spec1*zi*opzi -  
   2*x*li2spec2*zi*opzi +  
   2*sqrtxz1*x*li2spec2*zi*opzi -  
   x*li2spec19*zi*opzi -  
   4*sqrtxz1*x*li2spec21*zi*opzi -  
   2*sqrtxz1*x*li2spec22*zi* 
    opzi - 4*sqrtxz1*x*li2spec23*zi* 
    opzi + x*li2spec20*zi*opzi -  
   2*x*li2spec3*zi* 
    opzi - 2*sqrtxz1*x* 
    li2spec3*zi*opzi 
    + 2*x*li2spec4*zi* 
    opzi - 2*sqrtxz1*x* 
    li2spec4*zi*opzi 
    - 8*sqrtxz1*x*l2*zi*opzi -  
   4*sqrtxz1*x*lx*zi*opzi -  
   2*x*l2*lx*zi*opzi -  
   4*sqrtxz1*x*l2*lomz*zi*opzi -  
   2*sqrtxz1*x*lx*lomz*zi*opzi +  
   8*sqrtxz1*x*lspec1*zi*opzi +  
   4*x*l2*lspec1*zi*opzi +  
   8*sqrtxz1*x*l2*lspec1*zi*opzi -  
   2*sqrtxz1*x*lx*lspec1*zi*opzi +  
   4*sqrtxz1*x*lomz*lspec1*zi*opzi -  
   4*sqrtxz1*x*lz*zi*opzi -  
   6*x*l2*lz*zi*opzi -  
   8*sqrtxz1*x*l2*lz*zi*opzi -  
   x*lx*lz*zi*opzi +  
   sqrtxz1*x*lx*lz*zi*opzi -  
   2*sqrtxz1*x*lomz*lz*zi*opzi +  
   2*x*lspec1*lz*zi*opzi +  
   4*sqrtxz1*x*lspec1*lz*zi*opzi +  
   4*x*l2*lspec2*zi*opzi +  
   4*sqrtxz1*x*l2*lspec2*zi*opzi +  
   2*x*lx*lspec2*zi*opzi +  
   2*sqrtxz1*x*lx*lspec2*zi*opzi -  
   4*x*lspec1*lspec2*zi*opzi -  
   4*sqrtxz1*x*lspec1*lspec2*zi* 
    opzi + 4*x*lz*lspec2*zi*opzi +  
   4*sqrtxz1*x*lz*lspec2*zi*opzi +  
   x*lx*lxpz*zi*opzi -  
   x*lz*lxpz*zi*opzi -  
   x*lx*lopxz*zi*opzi -  
   x*lz*lopxz*zi*opzi +  
   2*sqrtxz1*x*lomz*lspec7*zi* 
    opzi - 2*x*li2spec1*NCi2*zi* 
    opzi - 2*sqrtxz1*x*li2spec1*NCi2* 
    zi*opzi + 2*x*li2spec2*NCi2* 
    zi*opzi - 2*sqrtxz1*x*li2spec2* 
    NCi2*zi*opzi +  
   x*li2spec19*NCi2*zi*opzi +  
   4*sqrtxz1*x*li2spec21*NCi2*zi* 
    opzi + 2*sqrtxz1*x* 
    li2spec22*NCi2*zi*opzi +  
   4*sqrtxz1*x*li2spec23*NCi2*zi* 
    opzi - x*li2spec20*NCi2*zi* 
    opzi + 2*x*li2spec3* 
    NCi2*zi*opzi +  
   2*sqrtxz1*x*li2spec3*NCi2* 
    zi*opzi - 2*x* 
    li2spec4*NCi2*zi* 
    opzi + 2*sqrtxz1*x* 
    li2spec4*NCi2*zi* 
    opzi + 8*sqrtxz1*x*l2*NCi2*zi*opzi +  
   4*sqrtxz1*x*lx*NCi2*zi*opzi +  
   2*x*l2*lx*NCi2*zi*opzi +  
   4*sqrtxz1*x*l2*lomz*NCi2*zi*opzi +  
   2*sqrtxz1*x*lx*lomz*NCi2*zi*opzi -  
   8*sqrtxz1*x*lspec1*NCi2*zi*opzi -  
   4*x*l2*lspec1*NCi2*zi*opzi -  
   8*sqrtxz1*x*l2*lspec1*NCi2*zi* 
    opzi + 2*sqrtxz1*x*lx*lspec1*NCi2* 
    zi*opzi - 4*sqrtxz1*x*lomz*lspec1* 
    NCi2*zi*opzi +  
   4*sqrtxz1*x*lz*NCi2*zi*opzi +  
   6*x*l2*lz*NCi2*zi*opzi +  
   8*sqrtxz1*x*l2*lz*NCi2*zi*opzi +  
   x*lx*lz*NCi2*zi*opzi -  
   sqrtxz1*x*lx*lz*NCi2*zi*opzi +  
   2*sqrtxz1*x*lomz*lz*NCi2*zi*opzi -  
   2*x*lspec1*lz*NCi2*zi*opzi -  
   4*sqrtxz1*x*lspec1*lz*NCi2*zi* 
    opzi - 4*x*l2*lspec2*NCi2*zi* 
    opzi - 4*sqrtxz1*x*l2*lspec2*NCi2* 
    zi*opzi - 2*x*lx*lspec2*NCi2* 
    zi*opzi - 2*sqrtxz1*x*lx*lspec2* 
    NCi2*zi*opzi +  
   4*x*lspec1*lspec2*NCi2*zi* 
    opzi + 4*sqrtxz1*x*lspec1*lspec2* 
    NCi2*zi*opzi -  
   4*x*lz*lspec2*NCi2*zi*opzi -  
   4*sqrtxz1*x*lz*lspec2*NCi2*zi* 
    opzi - x*lx*lxpz*NCi2*zi*opzi +  
   x*lz*lxpz*NCi2*zi*opzi +  
   x*lx*lopxz*NCi2*zi*opzi +  
   x*lz*lopxz*NCi2*zi*opzi -  
   2*sqrtxz1*x*lomz*lspec7*NCi2* 
    zi*opzi - (sqrtxz1*x*pi2*zi*opzi)/ 
    3. + (sqrtxz1*x*NCi2*pi2*zi*opzi)/3. -  
   2*li2spec1*omxi*zi*opzi -  
   2*sqrtxz1*li2spec1*omxi*zi* 
    opzi + 2*li2spec2*omxi*zi* 
    opzi - 2*sqrtxz1*li2spec2*omxi* 
    zi*opzi + li2spec19*omxi*zi* 
    opzi + 4*sqrtxz1*li2spec21*omxi* 
    zi*opzi + 2*sqrtxz1* 
    li2spec22*omxi*zi*opzi +  
   4*sqrtxz1*li2spec23*omxi*zi* 
    opzi - li2spec20*omxi*zi* 
    opzi + 2*li2spec3* 
    omxi*zi*opzi +  
   2*sqrtxz1*li2spec3*omxi* 
    zi*opzi - 2* 
    li2spec4*omxi*zi* 
    opzi + 2*sqrtxz1*li2spec4*omxi*zi*opzi +  
   8*sqrtxz1*l2*omxi*zi*opzi +  
   4*sqrtxz1*lx*omxi*zi*opzi +  
   2*l2*lx*omxi*zi*opzi +  
   4*sqrtxz1*l2*lomz*omxi*zi*opzi +  
   2*sqrtxz1*lx*lomz*omxi*zi*opzi -  
   8*sqrtxz1*lspec1*omxi*zi*opzi -  
   4*l2*lspec1*omxi*zi*opzi -  
   8*sqrtxz1*l2*lspec1*omxi*zi* 
    opzi + 2*sqrtxz1*lx*lspec1*omxi* 
    zi*opzi - 4*sqrtxz1*lomz*lspec1* 
    omxi*zi*opzi +  
   4*sqrtxz1*lz*omxi*zi*opzi +  
   6*l2*lz*omxi*zi*opzi +  
   8*sqrtxz1*l2*lz*omxi*zi*opzi +  
   lx*lz*omxi*zi*opzi -  
   sqrtxz1*lx*lz*omxi*zi*opzi +  
   2*sqrtxz1*lomz*lz*omxi*zi*opzi -  
   2*lspec1*lz*omxi*zi*opzi -  
   4*sqrtxz1*lspec1*lz*omxi*zi* 
    opzi - 4*l2*lspec2*omxi*zi* 
    opzi - 4*sqrtxz1*l2*lspec2*omxi* 
    zi*opzi - 2*lx*lspec2*omxi* 
    zi*opzi - 2*sqrtxz1*lx*lspec2* 
    omxi*zi*opzi +  
   4*lspec1*lspec2*omxi*zi* 
    opzi + 4*sqrtxz1*lspec1*lspec2* 
    omxi*zi*opzi -  
   4*lz*lspec2*omxi*zi*opzi -  
   4*sqrtxz1*lz*lspec2*omxi*zi* 
    opzi - lx*lxpz*omxi*zi*opzi +  
   lz*lxpz*omxi*zi*opzi +  
   lx*lopxz*omxi*zi*opzi +  
   lz*lopxz*omxi*zi*opzi -  
   2*sqrtxz1*lomz*lspec7*omxi* 
    zi*opzi + 2*li2spec1*NCi2* 
    omxi*zi*opzi +  
   2*sqrtxz1*li2spec1*NCi2*omxi*zi* 
    opzi - 2*li2spec2*NCi2*omxi* 
    zi*opzi + 2*sqrtxz1*li2spec2* 
    NCi2*omxi*zi*opzi -  
   li2spec19*NCi2*omxi*zi*opzi -  
   4*sqrtxz1*li2spec21*NCi2*omxi*zi* 
    opzi - 2*sqrtxz1*li2spec22*NCi2* 
    omxi*zi*opzi -  
   4*sqrtxz1*li2spec23*NCi2*omxi*zi* 
    opzi + li2spec20*NCi2*omxi*zi* 
    opzi - 2*li2spec3* 
    NCi2*omxi*zi*opzi -  
   2*sqrtxz1*li2spec3*NCi2* 
    omxi*zi*opzi +  
   2*li2spec4*NCi2* 
    omxi*zi*opzi -  
   2*sqrtxz1*li2spec4*NCi2* 
    omxi*zi*opzi -  
   8*sqrtxz1*l2*NCi2*omxi*zi*opzi -  
   4*sqrtxz1*lx*NCi2*omxi*zi*opzi -  
   2*l2*lx*NCi2*omxi*zi*opzi -  
   4*sqrtxz1*l2*lomz*NCi2*omxi*zi* 
    opzi - 2*sqrtxz1*lx*lomz*NCi2*omxi* 
    zi*opzi + 8*sqrtxz1*lspec1*NCi2* 
    omxi*zi*opzi +  
   4*l2*lspec1*NCi2*omxi*zi* 
    opzi + 8*sqrtxz1*l2*lspec1*NCi2* 
    omxi*zi*opzi -  
   2*sqrtxz1*lx*lspec1*NCi2*omxi*zi* 
    opzi + 4*sqrtxz1*lomz*lspec1*NCi2* 
    omxi*zi*opzi -  
   4*sqrtxz1*lz*NCi2*omxi*zi*opzi -  
   6*l2*lz*NCi2*omxi*zi*opzi -  
   8*sqrtxz1*l2*lz*NCi2*omxi*zi*opzi -  
   lx*lz*NCi2*omxi*zi*opzi +  
   sqrtxz1*lx*lz*NCi2*omxi*zi*opzi -  
   2*sqrtxz1*lomz*lz*NCi2*omxi*zi* 
    opzi + 2*lspec1*lz*NCi2*omxi* 
    zi*opzi + 4*sqrtxz1*lspec1*lz* 
    NCi2*omxi*zi*opzi +  
   4*l2*lspec2*NCi2*omxi*zi* 
    opzi + 4*sqrtxz1*l2*lspec2*NCi2* 
    omxi*zi*opzi +  
   2*lx*lspec2*NCi2*omxi*zi* 
    opzi + 2*sqrtxz1*lx*lspec2*NCi2* 
    omxi*zi*opzi -  
   4*lspec1*lspec2*NCi2*omxi* 
    zi*opzi - 4*sqrtxz1*lspec1* 
    lspec2*NCi2*omxi*zi*opzi +  
   4*lz*lspec2*NCi2*omxi*zi* 
    opzi + 4*sqrtxz1*lz*lspec2*NCi2* 
    omxi*zi*opzi +  
   lx*lxpz*NCi2*omxi*zi*opzi -  
   lz*lxpz*NCi2*omxi*zi*opzi -  
   lx*lopxz*NCi2*omxi*zi*opzi -  
   lz*lopxz*NCi2*omxi*zi*opzi +  
   2*sqrtxz1*lomz*lspec7*NCi2* 
    omxi*zi*opzi +  
   (sqrtxz1*pi2*omxi*zi*opzi)/3. -  
   (sqrtxz1*NCi2*pi2*omxi*zi*opzi)/3. -  
   4*x*l22 + 2*NC*x*l22 + 4*NC*z*l22 +  
   4*x*z*l22 + 4*x*NCi2*l22 -  
   4*x*z*NCi2*l22 - 2*x*NCi*l22 -  
   4*z*NCi*l22 - 24*sqrtxz1i*l22 +  
   12*x*sqrtxz1i*l22 - 6*x*z*sqrtxz1i*l22 +  
   24*NCi2*sqrtxz1i*l22 -  
   12*x*NCi2*sqrtxz1i*l22 +  
   6*x*z*NCi2*sqrtxz1i*l22 +  
   3*omxi*l22 - 3*NC*omxi*l22 -  
   3*z*omxi*l22 - 2*NC*z*omxi*l22 -  
   3*NCi2*omxi*l22 +  
   3*z*NCi2*omxi*l22 +  
   3*NCi*omxi*l22 +  
   2*z*NCi*omxi*l22 +  
   12*sqrtxz1i*omxi*l22 +  
   6*z*sqrtxz1i*omxi*l22 -  
   12*NCi2*sqrtxz1i*omxi*l22 -  
   6*z*NCi2*sqrtxz1i*omxi*l22 -  
   24*sqrtxz1i*x2*l22 +  
   24*NCi2*sqrtxz1i*x2*l22 +  
   4*x*zi*l22 + 6*sqrtxz1*x*zi*l22 -  
   4*x*NCi2*zi*l22 -  
   6*sqrtxz1*x*NCi2*zi*l22 -  
   4*omxi*zi*l22 -  
   6*sqrtxz1*omxi*zi*l22 +  
   4*NCi2*omxi*zi*l22 +  
   6*sqrtxz1*NCi2*omxi*zi*l22 +  
   8*NC*x*z2*l22 - 8*x*NCi*z2*l22 -  
   8*NC*omxi*z2*l22 +  
   8*NCi*omxi*z2*l22 +  
   24*sqrtxz1i*opzi*l22 -  
   18*x*sqrtxz1i*opzi*l22 +  
   6*x*z*sqrtxz1i*opzi*l22 -  
   24*NCi2*sqrtxz1i*opzi*l22 +  
   18*x*NCi2*sqrtxz1i*opzi*l22 -  
   6*x*z*NCi2*sqrtxz1i*opzi*l22 -  
   6*sqrtxz1i*omxi*opzi*l22 -  
   6*z*sqrtxz1i*omxi*opzi*l22 +  
   6*NCi2*sqrtxz1i*omxi*opzi*l22 +  
   6*z*NCi2*sqrtxz1i*omxi*opzi*l22 +  
   24*sqrtxz1i*x2*opzi*l22 -  
   24*NCi2*sqrtxz1i*x2*opzi*l22 -  
   4*x*zi*opzi*l22 -  
   6*sqrtxz1*x*zi*opzi*l22 +  
   4*x*NCi2*zi*opzi*l22 +  
   6*sqrtxz1*x*NCi2*zi*opzi*l22 +  
   4*omxi*zi*opzi*l22 +  
   6*sqrtxz1*omxi*zi*opzi*l22 -  
   4*NCi2*omxi*zi*opzi*l22 -  
   6*sqrtxz1*NCi2*omxi*zi*opzi*l22 +  
   (x*lx2)/2. - (NC*x*lx2)/2. + 2*NC*z*lx2 +  
   3*NC*x*z*lx2 - (x*NCi2*lx2)/2. +  
   (x*NCi*lx2)/2. - 2*z*NCi*lx2 -  
   3*x*z*NCi*lx2 + 6*sqrtxz1i*lx2 -  
   3*x*sqrtxz1i*lx2 +  
   (3*x*z*sqrtxz1i*lx2)/2. -  
   6*NCi2*sqrtxz1i*lx2 +  
   3*x*NCi2*sqrtxz1i*lx2 -  
   (3*x*z*NCi2*sqrtxz1i*lx2)/2. -  
   omxi*lx2 + NC*omxi*lx2 +  
   (z*omxi*lx2)/2. - 2*NC*z*omxi*lx2 +  
   NCi2*omxi*lx2 -  
   (z*NCi2*omxi*lx2)/2. -  
   NCi*omxi*lx2 +  
   2*z*NCi*omxi*lx2 -  
   3*sqrtxz1i*omxi*lx2 -  
   (3*z*sqrtxz1i*omxi*lx2)/2. +  
   3*NCi2*sqrtxz1i*omxi*lx2 +  
   (3*z*NCi2*sqrtxz1i*omxi*lx2)/2. +  
   (3*xi2*lx2)/2. - (3*NC*xi2*lx2)/4. +  
   (3*z*xi2*lx2)/2. + (3*NC*z*xi2*lx2)/2. -  
   (3*NCi2*xi2*lx2)/2. -  
   (3*z*NCi2*xi2*lx2)/2. +  
   (3*NCi*xi2*lx2)/4. -  
   (3*z*NCi*xi2*lx2)/2. +  
   (3*NC*xi*lx2)/4. - (3*NC*z*xi*lx2)/2. -  
   (3*NCi*xi*lx2)/4. +  
   (3*z*NCi*xi*lx2)/2. +  
   6*sqrtxz1i*x2*lx2 -  
   6*NCi2*sqrtxz1i*x2*lx2 +  
   (opxi*lx2)/2. - (5*NC*opxi*lx2)/4. +  
   (z*opxi*lx2)/2. +  
   (5*NC*z*opxi*lx2)/2. -  
   (NCi2*opxi*lx2)/2. -  
   (z*NCi2*opxi*lx2)/2. +  
   (5*NCi*opxi*lx2)/4. -  
   (5*z*NCi*opxi*lx2)/2. -  
   (3*xi2*opxi*lx2)/2. +  
   (3*NC*xi2*opxi*lx2)/4. -  
   (3*z*xi2*opxi*lx2)/2. -  
   (3*NC*z*xi2*opxi*lx2)/2. +  
   (3*NCi2*xi2*opxi*lx2)/2. +  
   (3*z*NCi2*xi2*opxi*lx2)/2. -  
   (3*NCi*xi2*opxi*lx2)/4. +  
   (3*z*NCi*xi2*opxi*lx2)/2. -  
   (3*xi*opxi*lx2)/2. -  
   (3*z*xi*opxi*lx2)/2. +  
   (3*NCi2*xi*opxi*lx2)/2. +  
   (3*z*NCi2*xi*opxi*lx2)/2. -  
   (3*x*omzi*lx2)/2. +  
   (3*x*NCi2*omzi*lx2)/2. -  
   (3*opxi*omzi*lx2)/2. +  
   (3*NCi2*opxi*omzi*lx2)/2. -  
   (zi*lx2)/4. - (NC*zi*lx2)/2. -  
   (7*x*zi*lx2)/4. - (NC*x*zi*lx2)/2. -  
   (3*sqrtxz1*x*zi*lx2)/2. +  
   (NCi2*zi*lx2)/4. +  
   (7*x*NCi2*zi*lx2)/4. +  
   (3*sqrtxz1*x*NCi2*zi*lx2)/2. +  
   (NCi*zi*lx2)/2. +  
   (x*NCi*zi*lx2)/2. +  
   (omxi*zi*lx2)/2. +  
   (3*sqrtxz1*omxi*zi*lx2)/2. -  
   (NCi2*omxi*zi*lx2)/2. -  
   (3*sqrtxz1*NCi2*omxi*zi*lx2)/2. +  
   (3*xi2*zi*lx2)/2. -  
   (3*NCi2*xi2*zi*lx2)/2. -  
   (3*opxi*zi*lx2)/2. +  
   (3*NCi2*opxi*zi*lx2)/2. -  
   (3*xi2*opxi*zi*lx2)/2. +  
   (3*NCi2*xi2*opxi*zi*lx2)/2. -  
   (3*xi*opxi*zi*lx2)/2. +  
   (3*NCi2*xi*opxi*zi*lx2)/2. +  
   (3*x*omzi*zi*lx2)/2. -  
   (3*x*NCi2*omzi*zi*lx2)/2. -  
   (3*xi2*omzi*zi*lx2)/2. +  
   (3*NCi2*xi2*omzi*zi*lx2)/2. +  
   (3*opxi*omzi*zi*lx2)/2. -  
   (3*NCi2*opxi*omzi*zi*lx2)/2. +  
   (3*xi2*opxi*omzi*zi*lx2)/2. -  
   (3*NCi2*xi2*opxi*omzi*zi* 
      lx2)/2. + (3*xi*opxi*omzi*zi* 
      lx2)/2. - (3*NCi2*xi*opxi*omzi* 
      zi*lx2)/2. - 2*NC*x*z2*lx2 +  
   2*x*NCi*z2*lx2 +  
   NC*omxi*z2*lx2 -  
   NCi*omxi*z2*lx2 -  
   3*NC*xi2*z2*lx2 +  
   3*NCi*xi2*z2*lx2 -  
   NC*opxi*z2*lx2 +  
   NCi*opxi*z2*lx2 +  
   3*NC*xi2*opxi*z2*lx2 -  
   3*NCi*xi2*opxi*z2*lx2 +  
   3*NC*xi*opxi*z2*lx2 -  
   3*NCi*xi*opxi*z2*lx2 -  
   6*sqrtxz1i*opzi*lx2 +  
   (9*x*sqrtxz1i*opzi*lx2)/2. -  
   (3*x*z*sqrtxz1i*opzi*lx2)/2. +  
   6*NCi2*sqrtxz1i*opzi*lx2 -  
   (9*x*NCi2*sqrtxz1i*opzi*lx2)/2. +  
   (3*x*z*NCi2*sqrtxz1i*opzi*lx2)/2. +  
   (3*sqrtxz1i*omxi*opzi*lx2)/2. +  
   (3*z*sqrtxz1i*omxi*opzi*lx2)/2. -  
   (3*NCi2*sqrtxz1i*omxi*opzi*lx2)/ 
    2. - (3*z*NCi2*sqrtxz1i*omxi*opzi* 
      lx2)/2. - 6*sqrtxz1i*x2*opzi* 
    lx2 + 6*NCi2*sqrtxz1i*x2*opzi* 
    lx2 + (3*sqrtxz1*x*zi*opzi*lx2)/2. -  
   (3*sqrtxz1*x*NCi2*zi*opzi*lx2)/2. -  
   (3*sqrtxz1*omxi*zi*opzi*lx2)/2. +  
   (3*sqrtxz1*NCi2*omxi*zi*opzi*lx2)/ 
    2. - 8*sqrtxz1i*lomz2 +  
   4*x*sqrtxz1i*lomz2 -  
   2*x*z*sqrtxz1i*lomz2 +  
   8*NCi2*sqrtxz1i*lomz2 -  
   4*x*NCi2*sqrtxz1i*lomz2 +  
   2*x*z*NCi2*sqrtxz1i*lomz2 +  
   4*sqrtxz1i*omxi*lomz2 +  
   2*z*sqrtxz1i*omxi*lomz2 -  
   4*NCi2*sqrtxz1i*omxi*lomz2 -  
   2*z*NCi2*sqrtxz1i*omxi*lomz2 -  
   8*sqrtxz1i*x2*lomz2 +  
   8*NCi2*sqrtxz1i*x2*lomz2 +  
   2*sqrtxz1*x*zi*lomz2 -  
   2*sqrtxz1*x*NCi2*zi*lomz2 -  
   2*sqrtxz1*omxi*zi*lomz2 +  
   2*sqrtxz1*NCi2*omxi*zi*lomz2 +  
   8*sqrtxz1i*opzi*lomz2 -  
   6*x*sqrtxz1i*opzi*lomz2 +  
   2*x*z*sqrtxz1i*opzi*lomz2 -  
   8*NCi2*sqrtxz1i*opzi*lomz2 +  
   6*x*NCi2*sqrtxz1i*opzi*lomz2 -  
   2*x*z*NCi2*sqrtxz1i*opzi*lomz2 -  
   2*sqrtxz1i*omxi*opzi*lomz2 -  
   2*z*sqrtxz1i*omxi*opzi*lomz2 +  
   2*NCi2*sqrtxz1i*omxi*opzi* 
    lomz2 + 2*z*NCi2*sqrtxz1i*omxi* 
    opzi*lomz2 +  
   8*sqrtxz1i*x2*opzi*lomz2 -  
   8*NCi2*sqrtxz1i*x2*opzi*lomz2 -  
   2*sqrtxz1*x*zi*opzi*lomz2 +  
   2*sqrtxz1*x*NCi2*zi*opzi*lomz2 +  
   2*sqrtxz1*omxi*zi*opzi*lomz2 -  
   2*sqrtxz1*NCi2*omxi*zi*opzi* 
    lomz2 - NC*lspec1_2 +  
   NC*x*lspec1_2 + 2*NC*z*lspec1_2 -  
   2*NC*x*z*lspec1_2 +  
   NCi*lspec1_2 -  
   x*NCi*lspec1_2 -  
   2*z*NCi*lspec1_2 +  
   2*x*z*NCi*lspec1_2 -  
   8*sqrtxz1i*lspec1_2 +  
   4*x*sqrtxz1i*lspec1_2 -  
   2*x*z*sqrtxz1i*lspec1_2 +  
   8*NCi2*sqrtxz1i*lspec1_2 -  
   4*x*NCi2*sqrtxz1i*lspec1_2 +  
   2*x*z*NCi2*sqrtxz1i*lspec1_2 -  
   NC*omxi*lspec1_2 +  
   2*NC*z*omxi*lspec1_2 +  
   NCi*omxi*lspec1_2 -  
   2*z*NCi*omxi*lspec1_2 +  
   4*sqrtxz1i*omxi*lspec1_2 +  
   2*z*sqrtxz1i*omxi*lspec1_2 -  
   4*NCi2*sqrtxz1i*omxi*lspec1_2 -  
   2*z*NCi2*sqrtxz1i*omxi*lspec1_2 -  
   8*sqrtxz1i*x2*lspec1_2 +  
   8*NCi2*sqrtxz1i*x2*lspec1_2 +  
   2*sqrtxz1*x*zi*lspec1_2 -  
   2*sqrtxz1*x*NCi2*zi*lspec1_2 -  
   2*sqrtxz1*omxi*zi*lspec1_2 +  
   2*sqrtxz1*NCi2*omxi*zi*lspec1_2 +  
   4*NC*x*z2*lspec1_2 -  
   4*x*NCi*z2*lspec1_2 -  
   4*NC*omxi*z2*lspec1_2 +  
   4*NCi*omxi*z2*lspec1_2 +  
   8*sqrtxz1i*opzi*lspec1_2 -  
   6*x*sqrtxz1i*opzi*lspec1_2 +  
   2*x*z*sqrtxz1i*opzi*lspec1_2 -  
   8*NCi2*sqrtxz1i*opzi*lspec1_2 +  
   6*x*NCi2*sqrtxz1i*opzi*lspec1_2 -  
   2*x*z*NCi2*sqrtxz1i*opzi* 
    lspec1_2 -  
   2*sqrtxz1i*omxi*opzi* 
    lspec1_2 -  
   2*z*sqrtxz1i*omxi*opzi* 
    lspec1_2 +  
   2*NCi2*sqrtxz1i*omxi*opzi* 
    lspec1_2 +  
   2*z*NCi2*sqrtxz1i*omxi*opzi* 
    lspec1_2 +  
   8*sqrtxz1i*x2*opzi*lspec1_2 -  
   8*NCi2*sqrtxz1i*x2*opzi* 
    lspec1_2 -  
   2*sqrtxz1*x*zi*opzi*lspec1_2 +  
   2*sqrtxz1*x*NCi2*zi*opzi* 
    lspec1_2 +  
   2*sqrtxz1*omxi*zi*opzi* 
    lspec1_2 -  
   2*sqrtxz1*NCi2*omxi*zi*opzi* 
    lspec1_2 - (x*lz2)/2. -  
   (3*NC*x*lz2)/2. - NC*z*lz2 +  
   (3*x*z*lz2)/2. + NC*x*z*lz2 +  
   (x*NCi2*lz2)/2. - (3*x*z*NCi2*lz2)/2. +  
   (3*x*NCi*lz2)/2. + z*NCi*lz2 -  
   x*z*NCi*lz2 - 10*sqrtxz1i*lz2 +  
   5*x*sqrtxz1i*lz2 -  
   (5*x*z*sqrtxz1i*lz2)/2. +  
   10*NCi2*sqrtxz1i*lz2 -  
   5*x*NCi2*sqrtxz1i*lz2 +  
   (5*x*z*NCi2*sqrtxz1i*lz2)/2. +  
   (omxi*lz2)/4. + (NC*omxi*lz2)/2. -  
   (z*omxi*lz2)/4. - NC*z*omxi*lz2 -  
   (NCi2*omxi*lz2)/4. +  
   (z*NCi2*omxi*lz2)/4. -  
   (NCi*omxi*lz2)/2. +  
   z*NCi*omxi*lz2 +  
   5*sqrtxz1i*omxi*lz2 +  
   (5*z*sqrtxz1i*omxi*lz2)/2. -  
   5*NCi2*sqrtxz1i*omxi*lz2 -  
   (5*z*NCi2*sqrtxz1i*omxi*lz2)/2. -  
   (xi2*lz2)/2. + (NC*xi2*lz2)/4. -  
   (z*xi2*lz2)/2. - (NC*z*xi2*lz2)/2. +  
   (NCi2*xi2*lz2)/2. +  
   (z*NCi2*xi2*lz2)/2. -  
   (NCi*xi2*lz2)/4. +  
   (z*NCi*xi2*lz2)/2. -  
   (NC*xi*lz2)/4. + (NC*z*xi*lz2)/2. +  
   (NCi*xi*lz2)/4. -  
   (z*NCi*xi*lz2)/2. -  
   10*sqrtxz1i*x2*lz2 +  
   10*NCi2*sqrtxz1i*x2*lz2 +  
   (opxi*lz2)/2. - (NC*opxi*lz2)/4. +  
   (z*opxi*lz2)/2. +  
   (NC*z*opxi*lz2)/2. -  
   (NCi2*opxi*lz2)/2. -  
   (z*NCi2*opxi*lz2)/2. +  
   (NCi*opxi*lz2)/4. -  
   (z*NCi*opxi*lz2)/2. +  
   (xi2*opxi*lz2)/2. -  
   (NC*xi2*opxi*lz2)/4. +  
   (z*xi2*opxi*lz2)/2. +  
   (NC*z*xi2*opxi*lz2)/2. -  
   (NCi2*xi2*opxi*lz2)/2. -  
   (z*NCi2*xi2*opxi*lz2)/2. +  
   (NCi*xi2*opxi*lz2)/4. -  
   (z*NCi*xi2*opxi*lz2)/2. +  
   (xi*opxi*lz2)/2. +  
   (z*xi*opxi*lz2)/2. -  
   (NCi2*xi*opxi*lz2)/2. -  
   (z*NCi2*xi*opxi*lz2)/2. +  
   (3*x*zi*lz2)/2. +  
   (5*sqrtxz1*x*zi*lz2)/2. -  
   (3*x*NCi2*zi*lz2)/2. -  
   (5*sqrtxz1*x*NCi2*zi*lz2)/2. -  
   omxi*zi*lz2 -  
   (5*sqrtxz1*omxi*zi*lz2)/2. +  
   NCi2*omxi*zi*lz2 +  
   (5*sqrtxz1*NCi2*omxi*zi*lz2)/2. -  
   (xi2*zi*lz2)/2. +  
   (NCi2*xi2*zi*lz2)/2. +  
   (opxi*zi*lz2)/2. -  
   (NCi2*opxi*zi*lz2)/2. +  
   (xi2*opxi*zi*lz2)/2. -  
   (NCi2*xi2*opxi*zi*lz2)/2. +  
   (xi*opxi*zi*lz2)/2. -  
   (NCi2*xi*opxi*zi*lz2)/2. -  
   (x*omzi*zi*lz2)/2. +  
   (x*NCi2*omzi*zi*lz2)/2. +  
   (xi2*omzi*zi*lz2)/2. -  
   (NCi2*xi2*omzi*zi*lz2)/2. -  
   (opxi*omzi*zi*lz2)/2. +  
   (NCi2*opxi*omzi*zi*lz2)/2. -  
   (xi2*opxi*omzi*zi*lz2)/2. +  
   (NCi2*xi2*opxi*omzi*zi*lz2)/ 
    2. - (xi*opxi*omzi*zi*lz2)/2. +  
   (NCi2*xi*opxi*omzi*zi*lz2)/ 
    2. - 2*NC*x*z2*lz2 +  
   2*x*NCi*z2*lz2 +  
   NC*omxi*z2*lz2 -  
   NCi*omxi*z2*lz2 +  
   NC*xi2*z2*lz2 -  
   NCi*xi2*z2*lz2 -  
   NC*opxi*z2*lz2 +  
   NCi*opxi*z2*lz2 -  
   NC*xi2*opxi*z2*lz2 +  
   NCi*xi2*opxi*z2*lz2 -  
   NC*xi*opxi*z2*lz2 +  
   NCi*xi*opxi*z2*lz2 +  
   10*sqrtxz1i*opzi*lz2 -  
   (15*x*sqrtxz1i*opzi*lz2)/2. +  
   (5*x*z*sqrtxz1i*opzi*lz2)/2. -  
   10*NCi2*sqrtxz1i*opzi*lz2 +  
   (15*x*NCi2*sqrtxz1i*opzi*lz2)/2. -  
   (5*x*z*NCi2*sqrtxz1i*opzi*lz2)/2. +  
   (omxi*opzi*lz2)/2. -  
   (NCi2*omxi*opzi*lz2)/2. -  
   (5*sqrtxz1i*omxi*opzi*lz2)/2. -  
   (5*z*sqrtxz1i*omxi*opzi*lz2)/2. +  
   (5*NCi2*sqrtxz1i*omxi*opzi*lz2)/ 
    2. + (5*z*NCi2*sqrtxz1i*omxi*opzi* 
      lz2)/2. + 10*sqrtxz1i*x2*opzi* 
    lz2 - 10*NCi2*sqrtxz1i*x2*opzi* 
    lz2 - x*zi*opzi*lz2 -  
   (5*sqrtxz1*x*zi*opzi*lz2)/2. +  
   x*NCi2*zi*opzi*lz2 +  
   (5*sqrtxz1*x*NCi2*zi*opzi*lz2)/2. +  
   omxi*zi*opzi*lz2 +  
   (5*sqrtxz1*omxi*zi*opzi*lz2)/2. -  
   NCi2*omxi*zi*opzi*lz2 -  
   (5*sqrtxz1*NCi2*omxi*zi*opzi*lz2)/ 
    2. - 2*sqrtxz1i*lspec7_2 +  
   x*sqrtxz1i*lspec7_2 -  
   (x*z*sqrtxz1i*lspec7_2)/2. +  
   2*NCi2*sqrtxz1i*lspec7_2 -  
   x*NCi2*sqrtxz1i*lspec7_2 +  
   (x*z*NCi2*sqrtxz1i*lspec7_2)/ 
    2. + sqrtxz1i*omxi* 
    lspec7_2 +  
   (z*sqrtxz1i*omxi*lspec7_2)/ 
    2. - NCi2*sqrtxz1i*omxi* 
    lspec7_2 -  
   (z*NCi2*sqrtxz1i*omxi* 
      lspec7_2)/2. -  
   2*sqrtxz1i*x2*lspec7_2 +  
   2*NCi2*sqrtxz1i*x2* 
    lspec7_2 +  
   (sqrtxz1*x*zi*lspec7_2)/2. -  
   (sqrtxz1*x*NCi2*zi*lspec7_2)/ 
    2. - (sqrtxz1*omxi*zi* 
      lspec7_2)/2. +  
   (sqrtxz1*NCi2*omxi*zi* 
      lspec7_2)/2. +  
   2*sqrtxz1i*opzi*lspec7_2 -  
   (3*x*sqrtxz1i*opzi*lspec7_2)/ 
    2. + (x*z*sqrtxz1i*opzi* 
      lspec7_2)/2. -  
   2*NCi2*sqrtxz1i*opzi* 
    lspec7_2 +  
   (3*x*NCi2*sqrtxz1i*opzi* 
      lspec7_2)/2. -  
   (x*z*NCi2*sqrtxz1i*opzi* 
      lspec7_2)/2. -  
   (sqrtxz1i*omxi*opzi* 
      lspec7_2)/2. -  
   (z*sqrtxz1i*omxi*opzi* 
      lspec7_2)/2. +  
   (NCi2*sqrtxz1i*omxi*opzi* 
      lspec7_2)/2. +  
   (z*NCi2*sqrtxz1i*omxi*opzi* 
      lspec7_2)/2. +  
   2*sqrtxz1i*x2*opzi* 
    lspec7_2 -  
   2*NCi2*sqrtxz1i*x2*opzi* 
    lspec7_2 -  
   (sqrtxz1*x*zi*opzi*lspec7_2)/ 
    2. + (sqrtxz1*x*NCi2*zi*opzi* 
      lspec7_2)/2. +  
   (sqrtxz1*omxi*zi*opzi* 
      lspec7_2)/2. -  
   (sqrtxz1*NCi2*omxi*zi*opzi* 
      lspec7_2)/2.;
};
double C2TQ2QP1_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double z2 = z * z;
  const double zi  = 1. / z;
  const double omxi = 1. / ( 1 - x );
  const double lx  = log(x);
  const double lz  = log(z);
  const double lz2 = lz * lz;
  const double lomx  = log(1 - x);
  const double lomz  = log(1 - z);
  const double Li2z  = apfel::dilog(z);
  return -0.16666666666666666*NC + (4*NC*x)/3. - (23*NC*z)/24. +  
   (13*NC*x*z)/24. + (NC*Li2z)/2. + (NC*x*Li2z)/2. + NC*x*z*Li2z -  
   (NC*lomx)/4. - (3*NC*x*lomx)/4. + (3*NC*z*lomx)/4. +  
   (NC*x*z*lomx)/4. + (NC*lx)/2. + (3*NC*x*lx)/2. -  
   (3*NC*z*lx)/2. - (NC*x*z*lx)/2. - (NC*lomz)/4. -  
   (3*NC*x*lomz)/4. + (3*NC*z*lomz)/4. + (NC*x*z*lomz)/4. -  
   (3*NC*lz)/4. - (NC*x*lz)/4. - (NC*z*lz)/4. +  
   (7*NC*x*z*lz)/4. - (NC*lomx*lz)/2. -  
   (NC*x*lomx*lz)/2. - NC*x*z*lomx*lz +  
   NC*lx*lz + NC*x*lx*lz + 2*NC*x*z*lx*lz +  
   NCi/6. - (4*x*NCi)/3. + (23*z*NCi)/24. -  
   (13*x*z*NCi)/24. - (Li2z*NCi)/2. -  
   (x*Li2z*NCi)/2. - x*z*Li2z*NCi +  
   (lomx*NCi)/4. + (3*x*lomx*NCi)/4. -  
   (3*z*lomx*NCi)/4. - (x*z*lomx*NCi)/4. -  
   (lx*NCi)/2. - (3*x*lx*NCi)/2. +  
   (3*z*lx*NCi)/2. + (x*z*lx*NCi)/2. +  
   (lomz*NCi)/4. + (3*x*lomz*NCi)/4. -  
   (3*z*lomz*NCi)/4. - (x*z*lomz*NCi)/4. +  
   (3*lz*NCi)/4. + (x*lz*NCi)/4. +  
   (z*lz*NCi)/4. - (7*x*z*lz*NCi)/4. +  
   (lomx*lz*NCi)/2. + (x*lomx*lz*NCi)/2. +  
   x*z*lomx*lz*NCi - lx*lz*NCi -  
   x*lx*lz*NCi - 2*x*z*lx*lz*NCi -  
   (NC*pi2)/12. - (NC*x*pi2)/12. - (NC*x*z*pi2)/6. +  
   (NCi*pi2)/12. + (x*NCi*pi2)/12. +  
   (x*z*NCi*pi2)/6. - (5*NC*lx*omxi)/4. +  
   (5*NC*z*lx*omxi)/4. - (3*NC*lx*lz*omxi)/2. -  
   (3*NC*z*lx*lz*omxi)/2. +  
   (5*lx*NCi*omxi)/4. -  
   (5*z*lx*NCi*omxi)/4. +  
   (3*lx*lz*NCi*omxi)/2. +  
   (3*z*lx*lz*NCi*omxi)/2. - (5*NC*zi)/72. +  
   (19*NC*x*zi)/72. - (NC*lomx*zi)/6. -  
   (NC*x*lomx*zi)/6. + (NC*lx*zi)/3. +  
   (NC*x*lx*zi)/3. - (NC*lomz*zi)/6. -  
   (NC*x*lomz*zi)/6. - (NC*lz*zi)/6. -  
   (NC*x*lz*zi)/6. + (5*NCi*zi)/72. -  
   (19*x*NCi*zi)/72. + (lomx*NCi*zi)/6. +  
   (x*lomx*NCi*zi)/6. -  
   (lx*NCi*zi)/3. - (x*lx*NCi*zi)/3. +  
   (lomz*NCi*zi)/6. +  
   (x*lomz*NCi*zi)/6. +  
   (lz*NCi*zi)/6. + (x*lz*NCi*zi)/6. -  
   (2*NC*lx*omxi*zi)/3. +  
   (2*lx*NCi*omxi*zi)/3. + (43*NC*z2)/36. -  
   (77*NC*x*z2)/36. - (NC*lomx*z2)/3. +  
   (2*NC*x*lomx*z2)/3. + (2*NC*lx*z2)/3. -  
   (4*NC*x*lx*z2)/3. - (NC*lomz*z2)/3. +  
   (2*NC*x*lomz*z2)/3. - (NC*lz*z2)/3. +  
   (2*NC*x*lz*z2)/3. - (43*NCi*z2)/36. +  
   (77*x*NCi*z2)/36. + (lomx*NCi*z2)/3. -  
   (2*x*lomx*NCi*z2)/3. -  
   (2*lx*NCi*z2)/3. + (4*x*lx*NCi*z2)/3. +  
   (lomz*NCi*z2)/3. -  
   (2*x*lomz*NCi*z2)/3. +  
   (lz*NCi*z2)/3. - (2*x*lz*NCi*z2)/3. +  
   (2*NC*lx*omxi*z2)/3. -  
   (2*lx*NCi*omxi*z2)/3. - (NC*lz2)/2. -  
   (NC*x*lz2)/2. - NC*x*z*lz2 +  
   (NCi*lz2)/2. + (x*NCi*lz2)/2. +  
   x*z*NCi*lz2;
};
double C2TQ2QP2_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double poly2    = 1 + 2*x + x*x - 4*x*z;
  const double poly2i   = 1. / poly2;
  const double sqrtxz2  = sqrt(poly2);
  const double sqrtxz2i = 1. / sqrtxz2;
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double x2 = x * x;
  const double x3 = x * x2;
  const double x4 = x * x3;
  const double x5 = x * x4;
  const double xi  = 1. / x;
  const double z2 = z * z;
  const double zi  = 1. / z;
  const double omzi  = 1. / ( 1 - z );
  const double lx  = log(x);
  const double lx2 = lx * lx;
  const double lz  = log(z);
  const double lomx  = log(1 - x);
  const double lomz  = log(1 - z);
  const double Li2x  = apfel::dilog(x);
  const double lspec5   = log(1 - sqrtxz2 + x);
  const double lspec6   = log(1 + sqrtxz2 + x);
  const double li2spec5  = apfel::dilog(0.5 - sqrtxz2/2. - x/2.);
  const double li2spec6  = apfel::dilog(0.5 + sqrtxz2/2. - x/2.);
  const double li2spec7  = apfel::dilog(0.5 - xi/2. - (sqrtxz2*xi)/2.);
  const double li2spec8  = apfel::dilog(0.5 - xi/2. + (sqrtxz2*xi)/2.);
  return (-7*NC)/6. + (NC*x)/6. + 5*NC*z - 5*NC*x*z + NC*Li2x +  
   NC*x*Li2x - (NC*lomx)/2. + (NC*x*lomx)/2. -  
   (19*NC*lx)/8. - (35*NC*x*lx)/8. + (5*NC*z*lx)/2. +  
   (5*NC*x*z*lx)/2. - (NC*lomz)/2. + (NC*x*lomz)/2. -  
   NC*lx*lomz - NC*x*lx*lomz - (5*NC*lz)/8. +  
   (5*NC*x*lz)/8. + (5*NC*z*lz)/2. - (5*NC*x*z*lz)/2. -  
   (NC*lx*lz)/2. - (NC*x*lx*lz)/2. + (7*NCi)/6. -  
   (x*NCi)/6. - 5*z*NCi + 5*x*z*NCi -  
   Li2x*NCi - x*Li2x*NCi + (lomx*NCi)/2. -  
   (x*lomx*NCi)/2. + (19*lx*NCi)/8. +  
   (35*x*lx*NCi)/8. - (5*z*lx*NCi)/2. -  
   (5*x*z*lx*NCi)/2. + (lomz*NCi)/2. -  
   (x*lomz*NCi)/2. + lx*lomz*NCi +  
   x*lx*lomz*NCi + (5*lz*NCi)/8. -  
   (5*x*lz*NCi)/8. - (5*z*lz*NCi)/2. +  
   (5*x*z*lz*NCi)/2. + (lx*lz*NCi)/2. +  
   (x*lx*lz*NCi)/2. - (NC*pi2)/6. -  
   (NC*x*pi2)/6. + (NCi*pi2)/6. +  
   (x*NCi*pi2)/6. - (NC*lx*poly2i)/8. +  
   (NC*lz*poly2i)/8. + (lx*NCi*poly2i)/8. -  
   (lz*NCi*poly2i)/8. +  
   (7*NC*li2spec5*sqrtxz2i)/8. +  
   (NC*x*li2spec5*sqrtxz2i)/8. -  
   (7*NC*z*li2spec5*sqrtxz2i)/4. -  
   5*NC*x*z*li2spec5*sqrtxz2i -  
   (7*NC*li2spec6*sqrtxz2i)/8. -  
   (NC*x*li2spec6*sqrtxz2i)/8. +  
   (7*NC*z*li2spec6*sqrtxz2i)/4. +  
   5*NC*x*z*li2spec6*sqrtxz2i -  
   (7*NC*li2spec7*sqrtxz2i)/ 
    8. - (NC*x*li2spec7* 
      sqrtxz2i)/8. + (7*NC*z* 
      li2spec7*sqrtxz2i)/4. +  
   5*NC*x*z*li2spec7* 
    sqrtxz2i + (7*NC*li2spec8* 
      sqrtxz2i)/8. + (NC*x* 
      li2spec8*sqrtxz2i)/8. -  
   (7*NC*z*li2spec8*sqrtxz2i)/ 
    4. - 5*NC*x*z*li2spec8* 
    sqrtxz2i - (7*NC*lx*lspec5*sqrtxz2i)/8. -  
   (NC*x*lx*lspec5*sqrtxz2i)/8. +  
   (7*NC*z*lx*lspec5*sqrtxz2i)/4. +  
   5*NC*x*z*lx*lspec5*sqrtxz2i +  
   (7*NC*lx*lspec6*sqrtxz2i)/8. +  
   (NC*x*lx*lspec6*sqrtxz2i)/8. -  
   (7*NC*z*lx*lspec6*sqrtxz2i)/4. -  
   5*NC*x*z*lx*lspec6*sqrtxz2i -  
   (7*li2spec5*NCi*sqrtxz2i)/8. -  
   (x*li2spec5*NCi*sqrtxz2i)/8. +  
   (7*z*li2spec5*NCi*sqrtxz2i)/4. +  
   5*x*z*li2spec5*NCi*sqrtxz2i +  
   (7*li2spec6*NCi*sqrtxz2i)/8. +  
   (x*li2spec6*NCi*sqrtxz2i)/8. -  
   (7*z*li2spec6*NCi*sqrtxz2i)/4. -  
   5*x*z*li2spec6*NCi*sqrtxz2i +  
   (7*li2spec7*NCi* 
      sqrtxz2i)/8. + (x*li2spec7*NCi*sqrtxz2i)/8. -  
   (7*z*li2spec7*NCi* 
      sqrtxz2i)/4. - 5*x*z* 
    li2spec7*NCi* 
    sqrtxz2i - (7*li2spec8* 
      NCi*sqrtxz2i)/8. -  
   (x*li2spec8*NCi* 
      sqrtxz2i)/8. + (7*z* 
      li2spec8*NCi* 
      sqrtxz2i)/4. + 5*x*z* 
    li2spec8*NCi* 
    sqrtxz2i + (7*lx*lspec5*NCi* 
      sqrtxz2i)/8. + (x*lx*lspec5*NCi* 
      sqrtxz2i)/8. - (7*z*lx*lspec5*NCi* 
      sqrtxz2i)/4. - 5*x*z*lx*lspec5*NCi* 
    sqrtxz2i - (7*lx*lspec6*NCi* 
      sqrtxz2i)/8. - (x*lx*lspec6*NCi* 
      sqrtxz2i)/8. + (7*z*lx*lspec6*NCi* 
      sqrtxz2i)/4. + 5*x*z*lx*lspec6*NCi* 
    sqrtxz2i + (NC*x*li2spec5*poly2i* 
      sqrtxz2i)/16. - (NC*x*li2spec6*poly2i* 
      sqrtxz2i)/16. - (NC*x* 
      li2spec7*poly2i* 
      sqrtxz2i)/16. + (NC*x* 
      li2spec8*poly2i* 
      sqrtxz2i)/16. - (NC*x*lx*lspec5*poly2i* 
      sqrtxz2i)/16. + (NC*x*lx*lspec6*poly2i* 
      sqrtxz2i)/16. - (x*li2spec5*NCi* 
      poly2i*sqrtxz2i)/16. +  
   (x*li2spec6*NCi*poly2i*sqrtxz2i)/ 
    16. + (x*li2spec7*NCi* 
      poly2i*sqrtxz2i)/16. -  
   (x*li2spec8*NCi* 
      poly2i*sqrtxz2i)/16. +  
   (x*lx*lspec5*NCi*poly2i*sqrtxz2i)/ 
    16. - (x*lx*lspec6*NCi*poly2i* 
      sqrtxz2i)/16. - (13*NC*xi)/9. -  
   (2*NC*lomx*xi)/3. - (NC*lx*xi)/8. -  
   (2*NC*lomz*xi)/3. - (19*NC*lz*xi)/24. +  
   (13*NCi*xi)/9. + (2*lomx*NCi*xi)/3. +  
   (lx*NCi*xi)/8. +  
   (2*lomz*NCi*xi)/3. +  
   (19*lz*NCi*xi)/24. +  
   (NC*lx*poly2i*xi)/8. +  
   (NC*lz*poly2i*xi)/8. -  
   (lx*NCi*poly2i*xi)/8. -  
   (lz*NCi*poly2i*xi)/8. +  
   (NC*li2spec5*sqrtxz2i*xi)/16. -  
   (NC*li2spec6*sqrtxz2i*xi)/16. -  
   (NC*li2spec7*sqrtxz2i* 
      xi)/16. + (NC*li2spec8* 
      sqrtxz2i*xi)/16. -  
   (NC*lx*lspec5*sqrtxz2i*xi)/16. +  
   (NC*lx*lspec6*sqrtxz2i*xi)/16. -  
   (li2spec5*NCi*sqrtxz2i*xi)/16. +  
   (li2spec6*NCi*sqrtxz2i*xi)/16. +  
   (li2spec7*NCi* 
      sqrtxz2i*xi)/16. -  
   (li2spec8*NCi* 
      sqrtxz2i*xi)/16. +  
   (lx*lspec5*NCi*sqrtxz2i*xi)/16. -  
   (lx*lspec6*NCi*sqrtxz2i*xi)/16. -  
   (NC*li2spec5*poly2i*sqrtxz2i*xi)/ 
    16. + (NC*li2spec6*poly2i*sqrtxz2i* 
      xi)/16. + (NC*li2spec7* 
      poly2i*sqrtxz2i*xi)/16. -  
   (NC*li2spec8*poly2i* 
      sqrtxz2i*xi)/16. +  
   (NC*lx*lspec5*poly2i*sqrtxz2i*xi)/ 
    16. - (NC*lx*lspec6*poly2i*sqrtxz2i* 
      xi)/16. + (li2spec5*NCi*poly2i* 
      sqrtxz2i*xi)/16. -  
   (li2spec6*NCi*poly2i*sqrtxz2i* 
      xi)/16. - (li2spec7* 
      NCi*poly2i*sqrtxz2i*xi)/16. +  
   (li2spec8*NCi*poly2i* 
      sqrtxz2i*xi)/16. -  
   (lx*lspec5*NCi*poly2i*sqrtxz2i* 
      xi)/16. + (lx*lspec6*NCi*poly2i* 
      sqrtxz2i*xi)/16. + (22*NC*x2)/9. +  
   (2*NC*lomx*x2)/3. - (17*NC*lx*x2)/8. +  
   (2*NC*lomz*x2)/3. + (19*NC*lz*x2)/24. -  
   (22*NCi*x2)/9. - (2*lomx*NCi*x2)/3. +  
   (17*lx*NCi*x2)/8. -  
   (2*lomz*NCi*x2)/3. -  
   (19*lz*NCi*x2)/24. +  
   (7*NC*li2spec5*sqrtxz2i*x2)/8. -  
   (7*NC*z*li2spec5*sqrtxz2i*x2)/4. -  
   (7*NC*li2spec6*sqrtxz2i*x2)/8. +  
   (7*NC*z*li2spec6*sqrtxz2i*x2)/4. -  
   (7*NC*li2spec7*sqrtxz2i* 
      x2)/8. + (7*NC*z*li2spec7* 
      sqrtxz2i*x2)/4. +  
   (7*NC*li2spec8*sqrtxz2i* 
      x2)/8. - (7*NC*z*li2spec8* 
      sqrtxz2i*x2)/4. -  
   (7*NC*lx*lspec5*sqrtxz2i*x2)/8. +  
   (7*NC*z*lx*lspec5*sqrtxz2i*x2)/4. +  
   (7*NC*lx*lspec6*sqrtxz2i*x2)/8. -  
   (7*NC*z*lx*lspec6*sqrtxz2i*x2)/4. -  
   (7*li2spec5*NCi*sqrtxz2i*x2)/8. +  
   (7*z*li2spec5*NCi*sqrtxz2i*x2)/ 
    4. + (7*li2spec6*NCi*sqrtxz2i*x2)/ 
    8. - (7*z*li2spec6*NCi*sqrtxz2i* 
      x2)/4. + (7*li2spec7* 
      NCi*sqrtxz2i*x2)/8. -  
   (7*z*li2spec7*NCi* 
      sqrtxz2i*x2)/4. -  
   (7*li2spec8*NCi* 
      sqrtxz2i*x2)/8. +  
   (7*z*li2spec8*NCi* 
      sqrtxz2i*x2)/4. +  
   (7*lx*lspec5*NCi*sqrtxz2i*x2)/8. -  
   (7*z*lx*lspec5*NCi*sqrtxz2i*x2)/4. -  
   (7*lx*lspec6*NCi*sqrtxz2i*x2)/8. +  
   (7*z*lx*lspec6*NCi*sqrtxz2i*x2)/4. -  
   (NC*lx*poly2i*x3)/8. -  
   (NC*lz*poly2i*x3)/8. +  
   (lx*NCi*poly2i*x3)/8. +  
   (lz*NCi*poly2i*x3)/8. +  
   (NC*li2spec5*sqrtxz2i*x3)/16. -  
   (NC*li2spec6*sqrtxz2i*x3)/16. -  
   (NC*li2spec7*sqrtxz2i* 
      x3)/16. + (NC*li2spec8* 
      sqrtxz2i*x3)/16. -  
   (NC*lx*lspec5*sqrtxz2i*x3)/16. +  
   (NC*lx*lspec6*sqrtxz2i*x3)/16. -  
   (li2spec5*NCi*sqrtxz2i*x3)/16. +  
   (li2spec6*NCi*sqrtxz2i*x3)/16. +  
   (li2spec7*NCi* 
      sqrtxz2i*x3)/16. -  
   (li2spec8*NCi* 
      sqrtxz2i*x3)/16. +  
   (lx*lspec5*NCi*sqrtxz2i*x3)/16. -  
   (lx*lspec6*NCi*sqrtxz2i*x3)/16. +  
   (NC*li2spec5*poly2i*sqrtxz2i*x3)/ 
    16. - (NC*li2spec6*poly2i*sqrtxz2i* 
      x3)/16. - (NC*li2spec7* 
      poly2i*sqrtxz2i*x3)/16. +  
   (NC*li2spec8*poly2i* 
      sqrtxz2i*x3)/16. -  
   (NC*lx*lspec5*poly2i*sqrtxz2i*x3)/ 
    16. + (NC*lx*lspec6*poly2i*sqrtxz2i* 
      x3)/16. - (li2spec5*NCi*poly2i* 
      sqrtxz2i*x3)/16. +  
   (li2spec6*NCi*poly2i*sqrtxz2i* 
      x3)/16. + (li2spec7* 
      NCi*poly2i*sqrtxz2i*x3)/16. -  
   (li2spec8*NCi*poly2i* 
      sqrtxz2i*x3)/16. +  
   (lx*lspec5*NCi*poly2i*sqrtxz2i* 
      x3)/16. - (lx*lspec6*NCi*poly2i* 
      sqrtxz2i*x3)/16. + (NC*lx*poly2i*x4)/8. -  
   (NC*lz*poly2i*x4)/8. -  
   (lx*NCi*poly2i*x4)/8. +  
   (lz*NCi*poly2i*x4)/8. -  
   (NC*li2spec5*poly2i*sqrtxz2i*x5)/ 
    16. + (NC*li2spec6*poly2i*sqrtxz2i* 
      x5)/16. + (NC*li2spec7* 
      poly2i*sqrtxz2i*x5)/16. -  
   (NC*li2spec8*poly2i* 
      sqrtxz2i*x5)/16. +  
   (NC*lx*lspec5*poly2i*sqrtxz2i*x5)/ 
    16. - (NC*lx*lspec6*poly2i*sqrtxz2i* 
      x5)/16. + (li2spec5*NCi*poly2i* 
      sqrtxz2i*x5)/16. -  
   (li2spec6*NCi*poly2i*sqrtxz2i* 
      x5)/16. - (li2spec7* 
      NCi*poly2i*sqrtxz2i*x5)/16. +  
   (li2spec8*NCi*poly2i* 
      sqrtxz2i*x5)/16. -  
   (lx*lspec5*NCi*poly2i*sqrtxz2i* 
      x5)/16. + (lx*lspec6*NCi*poly2i* 
      sqrtxz2i*x5)/16. - (NC*lz*omzi)/2. +  
   (NC*x*lz*omzi)/2. + (lz*NCi*omzi)/2. -  
   (x*lz*NCi*omzi)/2. +  
   (NC*lz*xi*omzi)/3. -  
   (lz*NCi*xi*omzi)/3. -  
   (NC*lz*x2*omzi)/3. +  
   (lz*NCi*x2*omzi)/3. - (11*NC*zi)/12. +  
   (17*NC*x*zi)/12. - (NC*Li2x*zi)/2. -  
   (NC*x*Li2x*zi)/2. + (NC*lomx*zi)/4. -  
   (NC*x*lomx*zi)/4. + (NC*lx*zi)/2. +  
   (3*NC*x*lx*zi)/2. + (NC*lomz*zi)/4. -  
   (NC*x*lomz*zi)/4. + (NC*lx*lomz*zi)/2. +  
   (NC*x*lx*lomz*zi)/2. + (NC*lz*zi)/4. -  
   (NC*x*lz*zi)/4. + (NC*lx*lz*zi)/2. +  
   (NC*x*lx*lz*zi)/2. + (11*NCi*zi)/12. -  
   (17*x*NCi*zi)/12. + (Li2x*NCi*zi)/2. +  
   (x*Li2x*NCi*zi)/2. -  
   (lomx*NCi*zi)/4. +  
   (x*lomx*NCi*zi)/4. -  
   (lx*NCi*zi)/2. - (3*x*lx*NCi*zi)/2. -  
   (lomz*NCi*zi)/4. +  
   (x*lomz*NCi*zi)/4. -  
   (lx*lomz*NCi*zi)/2. -  
   (x*lx*lomz*NCi*zi)/2. -  
   (lz*NCi*zi)/4. + (x*lz*NCi*zi)/4. -  
   (lx*lz*NCi*zi)/2. -  
   (x*lx*lz*NCi*zi)/2. + (NC*pi2*zi)/12. +  
   (NC*x*pi2*zi)/12. - (NCi*pi2*zi)/12. -  
   (x*NCi*pi2*zi)/12. + (13*NC*xi*zi)/18. +  
   (NC*lomx*xi*zi)/3. +  
   (NC*lomz*xi*zi)/3. +  
   (NC*lz*xi*zi)/3. -  
   (13*NCi*xi*zi)/18. -  
   (lomx*NCi*xi*zi)/3. -  
   (lomz*NCi*xi*zi)/3. -  
   (lz*NCi*xi*zi)/3. -  
   (11*NC*x2*zi)/9. - (NC*lomx*x2*zi)/3. +  
   NC*lx*x2*zi - (NC*lomz*x2*zi)/3. -  
   (NC*lz*x2*zi)/3. +  
   (11*NCi*x2*zi)/9. +  
   (lomx*NCi*x2*zi)/3. -  
   lx*NCi*x2*zi +  
   (lomz*NCi*x2*zi)/3. +  
   (lz*NCi*x2*zi)/3. +  
   5*NC*x*li2spec5*sqrtxz2i*z2 -  
   5*NC*x*li2spec6*sqrtxz2i*z2 -  
   5*NC*x*li2spec7*sqrtxz2i* 
    z2 + 5*NC*x*li2spec8* 
    sqrtxz2i*z2 - 5*NC*x*lx*lspec5* 
    sqrtxz2i*z2 + 5*NC*x*lx*lspec6* 
    sqrtxz2i*z2 - 5*x*li2spec5*NCi* 
    sqrtxz2i*z2 + 5*x*li2spec6*NCi* 
    sqrtxz2i*z2 + 5*x* 
    li2spec7*NCi* 
    sqrtxz2i*z2 - 5*x* 
    li2spec8*NCi* 
    sqrtxz2i*z2 + 5*x*lx*lspec5*NCi* 
    sqrtxz2i*z2 - 5*x*lx*lspec6*NCi* 
    sqrtxz2i*z2 + NC*lx2 + NC*x*lx2 -  
   NCi*lx2 - x*NCi*lx2 -  
   (NC*zi*lx2)/2. - (NC*x*zi*lx2)/2. +  
   (NCi*zi*lx2)/2. +  
   (x*NCi*zi*lx2)/2.;
};
double C2TQ2QP3_Rx_Rz(int const& nf, double const& x, double const& z)
{
  const double sqrtxz1  = sqrt(1 - 2*z + z*z + 4*x*z);
  const double poly2    = 1 + 2*x + x*x - 4*x*z;
  const double poly2i   = 1. / poly2;
  const double sqrtxz2  = sqrt(poly2);
  const double sqrtxz2i = 1. / sqrtxz2;
  const double sqrtxz3  = sqrt(x/z);
  const double NC   = apfel::NC;
  const double NCi  = 1. / NC;
  const double l2  = log(2);
  const double l22 = l2 * l2;
  const double pi  = M_PI;
  const double pi2 = pi * pi;
  const double x2 = x * x;
  const double x3 = x * x2;
  const double x4 = x * x3;
  const double x5 = x * x4;
  const double xi  = 1. / x;
  const double xi2 = xi * xi;
  const double z2 = z * z;
  const double zi  = 1. / z;
  const double omxi = 1. / ( 1 - x );
  const double opxi = 1. / ( 1 + x );
  const double lx  = log(x);
  const double lx2 = lx * lx;
  const double lz  = log(z);
  const double lz2 = lz * lz;
  const double lomx  = log(1 - x);
  const double lomz  = log(1 - z);
  const double lopx  = log(1 + x);
  const double Li2x  = apfel::dilog(x);
  const double Li2mx = apfel::dilog(-x);
  const double Li2z  = apfel::dilog(z);
  const double xmzi  = 1. / ( x - z );
  const double lxpz    = log(x + z);
  const double lopxz   = log(1 + x*z);
  const double lopxzi  = log(1 + x*zi);
  const double li2omxzi = apfel::dilog(1 - x*zi);
  const double lspec1   = log(1 + sqrtxz1 - z);
  const double lspec1_2 = lspec1 * lspec1;
  const double lspec2   = log(1 + sqrtxz1 + z);
  const double lspec3   = log(sqrtxz3);
  const double lspec4   = log(sqrtxz3*z);
  const double lspec5   = log(1 - sqrtxz2 + x);
  const double lspec6   = log(1 + sqrtxz2 + x);
  const double li2spec1  = apfel::dilog(0.5 - sqrtxz1/2. - z/2.);
  const double li2spec2  = apfel::dilog(0.5 - sqrtxz1/2. + z/2.);
  const double li2spec3  = apfel::dilog(0.5 - zi/2. - (sqrtxz1*zi)/2.);
  const double li2spec4  = apfel::dilog(0.5 + zi/2. - (sqrtxz1*zi)/2.);
  const double li2spec5  = apfel::dilog(0.5 - sqrtxz2/2. - x/2.);
  const double li2spec6  = apfel::dilog(0.5 + sqrtxz2/2. - x/2.);
  const double li2spec7  = apfel::dilog(0.5 - xi/2. - (sqrtxz2*xi)/2.);
  const double li2spec8  = apfel::dilog(0.5 - xi/2. + (sqrtxz2*xi)/2.);
  const double li2spec19 = apfel::dilog(-(x*z));
  const double li2spec20 = apfel::dilog(-(x*zi));
  const double atanspec1 = atan(sqrtxz3);
  const double atanspec2 = atan(sqrtxz3*z);
  const double itani1 = InvTanInt(-sqrtxz3);
  const double itani2 = InvTanInt(sqrtxz3);
  const double itani3 = InvTanInt(sqrtxz3*z);
  return 3*NC - (NC*x)/2. - (5*NC*z)/2. + (NC*x*z)/2. +  
   (3*NC*sqrtxz3*itani1)/4. -  
   (3*NC*sqrtxz3*x*z*itani1)/4. -  
   (3*NC*sqrtxz3*itani2)/4. +  
   (3*NC*sqrtxz3*x*z*itani2)/4. +  
   (3*NC*sqrtxz3*itani3)/2. -  
   (3*NC*sqrtxz3*x*z*itani3)/2. - 2*NC*Li2mx -  
   2*NC*x*Li2mx + 4*NC*z*Li2mx + 4*NC*x*z*Li2mx - NC*x*Li2x +  
   2*NC*x*z*Li2x + NC*x*li2spec1 +  
   2*NC*z*li2spec1 - NC*x*li2spec2 -  
   2*NC*z*li2spec2 - NC*x*Li2z + 2*NC*x*z*Li2z -  
   NC*x*li2spec19 - 2*NC*z*li2spec19 + NC*li2spec20 +  
   2*NC*x*z*li2spec20 -  
   NC*li2spec3 -  
   2*NC*x*z*li2spec3 +  
   NC*li2spec4 +  
   2*NC*x*z*li2spec4 +  
   NC*x*li2omxzi - 2*NC*x*z*li2omxzi -  
   2*NC*sqrtxz1*x*l2 + (3*NC*sqrtxz3*atanspec1*lspec3)/2. -  
   (3*NC*sqrtxz3*x*z*atanspec1*lspec3)/2. - (NC*lx)/8. +  
   (11*NC*x*lx)/8. - NC*sqrtxz1*x*lx - (NC*z*lx)/4. -  
   (5*NC*x*z*lx)/4. + NC*l2*lx - 2*NC*x*l2*lx -  
   4*NC*z*l2*lx + 2*NC*x*z*l2*lx - NC*x*lomx*lx +  
   2*NC*x*z*lomx*lx - 2*NC*lx*lopx -  
   2*NC*x*lx*lopx + 4*NC*z*lx*lopx +  
   4*NC*x*z*lx*lopx + 2*NC*sqrtxz1*x*lspec1 -  
   NC*l2*lspec1 + 3*NC*x*l2*lspec1 +  
   6*NC*z*l2*lspec1 -  
   2*NC*x*z*l2*lspec1 - NC*lx*lspec1 +  
   NC*x*lx*lspec1 + 2*NC*z*lx*lspec1 -  
   2*NC*x*z*lx*lspec1 + (17*NC*lz)/8. -  
   (NC*x*lz)/8. - NC*sqrtxz1*x*lz - (NC*z*lz)/4. +  
   (5*NC*x*z*lz)/4. - NC*l2*lz - 2*NC*x*l2*lz -  
   4*NC*z*l2*lz - 2*NC*x*z*l2*lz + (NC*lx*lz)/2. -  
   (NC*x*lx*lz)/2. - NC*z*lx*lz + NC*x*z*lx*lz -  
   NC*x*lomz*lz + 2*NC*x*z*lomz*lz +  
   NC*x*lspec1*lz + 2*NC*z*lspec1*lz -  
   (3*NC*sqrtxz3*atanspec2*lspec4)/2. +  
   (3*NC*sqrtxz3*x*z*atanspec2*lspec4)/2. +  
   NC*l2*lspec2 + NC*x*l2*lspec2 +  
   2*NC*z*l2*lspec2 +  
   2*NC*x*z*l2*lspec2 + NC*x*lx*lspec2 +  
   2*NC*z*lx*lspec2 -  
   NC*lspec1*lspec2 -  
   NC*x*lspec1*lspec2 -  
   2*NC*z*lspec1*lspec2 -  
   2*NC*x*z*lspec1*lspec2 +  
   NC*lz*lspec2 + NC*x*lz*lspec2 +  
   2*NC*z*lz*lspec2 +  
   2*NC*x*z*lz*lspec2 + (NC*lx*lxpz)/2. +  
   (NC*x*lx*lxpz)/2. + NC*z*lx*lxpz +  
   NC*x*z*lx*lxpz - (NC*lz*lxpz)/2. -  
   (NC*x*lz*lxpz)/2. - NC*z*lz*lxpz -  
   NC*x*z*lz*lxpz - NC*x*lx*lopxz -  
   2*NC*z*lx*lopxz - NC*x*lz*lopxz -  
   2*NC*z*lz*lopxz + (NC*lx*lopxzi)/2. -  
   (NC*x*lx*lopxzi)/2. - NC*z*lx*lopxzi +  
   NC*x*z*lx*lopxzi - (NC*lz*lopxzi)/2. +  
   (NC*x*lz*lopxzi)/2. + NC*z*lz*lopxzi -  
   NC*x*z*lz*lopxzi - 3*NCi + (x*NCi)/2. +  
   (5*z*NCi)/2. - (x*z*NCi)/2. -  
   (3*sqrtxz3*itani1*NCi)/4. +  
   (3*sqrtxz3*x*z*itani1*NCi)/4. +  
   (3*sqrtxz3*itani2*NCi)/4. -  
   (3*sqrtxz3*x*z*itani2*NCi)/4. -  
   (3*sqrtxz3*itani3*NCi)/2. +  
   (3*sqrtxz3*x*z*itani3*NCi)/2. +  
   2*Li2mx*NCi + 2*x*Li2mx*NCi - 4*z*Li2mx*NCi -  
   4*x*z*Li2mx*NCi + x*Li2x*NCi - 2*x*z*Li2x*NCi -  
   x*li2spec1*NCi -  
   2*z*li2spec1*NCi +  
   x*li2spec2*NCi +  
   2*z*li2spec2*NCi + x*Li2z*NCi -  
   2*x*z*Li2z*NCi + x*li2spec19*NCi +  
   2*z*li2spec19*NCi - li2spec20*NCi -  
   2*x*z*li2spec20*NCi +  
   li2spec3*NCi +  
   2*x*z*li2spec3*NCi -  
   li2spec4*NCi -  
   2*x*z*li2spec4*NCi -  
   x*li2omxzi*NCi +  
   2*x*z*li2omxzi*NCi + 2*sqrtxz1*x*l2*NCi -  
   (3*sqrtxz3*atanspec1*lspec3*NCi)/2. +  
   (3*sqrtxz3*x*z*atanspec1*lspec3*NCi)/2. +  
   (lx*NCi)/8. - (11*x*lx*NCi)/8. +  
   sqrtxz1*x*lx*NCi + (z*lx*NCi)/4. +  
   (5*x*z*lx*NCi)/4. - l2*lx*NCi +  
   2*x*l2*lx*NCi + 4*z*l2*lx*NCi -  
   2*x*z*l2*lx*NCi + x*lomx*lx*NCi -  
   2*x*z*lomx*lx*NCi + 2*lx*lopx*NCi +  
   2*x*lx*lopx*NCi - 4*z*lx*lopx*NCi -  
   4*x*z*lx*lopx*NCi -  
   2*sqrtxz1*x*lspec1*NCi +  
   l2*lspec1*NCi -  
   3*x*l2*lspec1*NCi -  
   6*z*l2*lspec1*NCi +  
   2*x*z*l2*lspec1*NCi +  
   lx*lspec1*NCi -  
   x*lx*lspec1*NCi -  
   2*z*lx*lspec1*NCi +  
   2*x*z*lx*lspec1*NCi - (17*lz*NCi)/8. +  
   (x*lz*NCi)/8. + sqrtxz1*x*lz*NCi +  
   (z*lz*NCi)/4. - (5*x*z*lz*NCi)/4. +  
   l2*lz*NCi + 2*x*l2*lz*NCi +  
   4*z*l2*lz*NCi + 2*x*z*l2*lz*NCi -  
   (lx*lz*NCi)/2. + (x*lx*lz*NCi)/2. +  
   z*lx*lz*NCi - x*z*lx*lz*NCi +  
   x*lomz*lz*NCi - 2*x*z*lomz*lz*NCi -  
   x*lspec1*lz*NCi -  
   2*z*lspec1*lz*NCi +  
   (3*sqrtxz3*atanspec2*lspec4*NCi)/2. -  
   (3*sqrtxz3*x*z*atanspec2*lspec4*NCi)/2. -  
   l2*lspec2*NCi -  
   x*l2*lspec2*NCi -  
   2*z*l2*lspec2*NCi -  
   2*x*z*l2*lspec2*NCi -  
   x*lx*lspec2*NCi -  
   2*z*lx*lspec2*NCi +  
   lspec1*lspec2*NCi +  
   x*lspec1*lspec2*NCi +  
   2*z*lspec1*lspec2*NCi +  
   2*x*z*lspec1*lspec2*NCi -  
   lz*lspec2*NCi -  
   x*lz*lspec2*NCi -  
   2*z*lz*lspec2*NCi -  
   2*x*z*lz*lspec2*NCi -  
   (lx*lxpz*NCi)/2. - (x*lx*lxpz*NCi)/2. -  
   z*lx*lxpz*NCi - x*z*lx*lxpz*NCi +  
   (lz*lxpz*NCi)/2. + (x*lz*lxpz*NCi)/2. +  
   z*lz*lxpz*NCi + x*z*lz*lxpz*NCi +  
   x*lx*lopxz*NCi + 2*z*lx*lopxz*NCi +  
   x*lz*lopxz*NCi + 2*z*lz*lopxz*NCi -  
   (lx*lopxzi*NCi)/2. +  
   (x*lx*lopxzi*NCi)/2. +  
   z*lx*lopxzi*NCi -  
   x*z*lx*lopxzi*NCi +  
   (lz*lopxzi*NCi)/2. -  
   (x*lz*lopxzi*NCi)/2. -  
   z*lz*lopxzi*NCi +  
   x*z*lz*lopxzi*NCi - (NC*pi2)/6. +  
   (NC*x*pi2)/6. + (NC*z*pi2)/3. - (NC*x*z*pi2)/3. +  
   (NCi*pi2)/6. - (x*NCi*pi2)/6. -  
   (z*NCi*pi2)/3. + (x*z*NCi*pi2)/3. +  
   (NC*lx*poly2i)/8. - (NC*lz*poly2i)/8. -  
   (lx*NCi*poly2i)/8. +  
   (lz*NCi*poly2i)/8. +  
   (9*NC*li2spec5*sqrtxz2i)/8. +  
   (7*NC*x*li2spec5*sqrtxz2i)/8. -  
   (9*NC*z*li2spec5*sqrtxz2i)/4. -  
   3*NC*x*z*li2spec5*sqrtxz2i -  
   (9*NC*li2spec6*sqrtxz2i)/8. -  
   (7*NC*x*li2spec6*sqrtxz2i)/8. +  
   (9*NC*z*li2spec6*sqrtxz2i)/4. +  
   3*NC*x*z*li2spec6*sqrtxz2i -  
   (9*NC*li2spec7*sqrtxz2i)/ 
    8. - (7*NC*x*li2spec7* 
      sqrtxz2i)/8. + (9*NC*z* 
      li2spec7*sqrtxz2i)/4. +  
   3*NC*x*z*li2spec7* 
    sqrtxz2i + (9*NC*li2spec8* 
      sqrtxz2i)/8. + (7*NC*x* 
      li2spec8*sqrtxz2i)/8. -  
   (9*NC*z*li2spec8*sqrtxz2i)/ 
    4. - 3*NC*x*z*li2spec8* 
    sqrtxz2i - (9*NC*lx*lspec5*sqrtxz2i)/8. -  
   (7*NC*x*lx*lspec5*sqrtxz2i)/8. +  
   (9*NC*z*lx*lspec5*sqrtxz2i)/4. +  
   3*NC*x*z*lx*lspec5*sqrtxz2i +  
   (9*NC*lx*lspec6*sqrtxz2i)/8. +  
   (7*NC*x*lx*lspec6*sqrtxz2i)/8. -  
   (9*NC*z*lx*lspec6*sqrtxz2i)/4. -  
   3*NC*x*z*lx*lspec6*sqrtxz2i -  
   (9*li2spec5*NCi*sqrtxz2i)/8. -  
   (7*x*li2spec5*NCi*sqrtxz2i)/8. +  
   (9*z*li2spec5*NCi*sqrtxz2i)/4. +  
   3*x*z*li2spec5*NCi*sqrtxz2i +  
   (9*li2spec6*NCi*sqrtxz2i)/8. +  
   (7*x*li2spec6*NCi*sqrtxz2i)/8. -  
   (9*z*li2spec6*NCi*sqrtxz2i)/4. -  
   3*x*z*li2spec6*NCi*sqrtxz2i +  
   (9*li2spec7*NCi* 
      sqrtxz2i)/8. + (7*x* 
      li2spec7*NCi* 
      sqrtxz2i)/8. - (9*z* 
      li2spec7*NCi* 
      sqrtxz2i)/4. - 3*x*z* 
    li2spec7*NCi* 
    sqrtxz2i - (9*li2spec8* 
      NCi*sqrtxz2i)/8. -  
   (7*x*li2spec8*NCi* 
      sqrtxz2i)/8. + (9*z* 
      li2spec8*NCi* 
      sqrtxz2i)/4. + 3*x*z* 
    li2spec8*NCi* 
    sqrtxz2i + (9*lx*lspec5*NCi* 
      sqrtxz2i)/8. + (7*x*lx*lspec5*NCi* 
      sqrtxz2i)/8. - (9*z*lx*lspec5*NCi* 
      sqrtxz2i)/4. - 3*x*z*lx*lspec5*NCi* 
    sqrtxz2i - (9*lx*lspec6*NCi* 
      sqrtxz2i)/8. - (7*x*lx*lspec6*NCi* 
      sqrtxz2i)/8. + (9*z*lx*lspec6*NCi* 
      sqrtxz2i)/4. + 3*x*z*lx*lspec6*NCi* 
    sqrtxz2i - (NC*x*li2spec5*poly2i* 
      sqrtxz2i)/16. + (NC*x*li2spec6*poly2i* 
      sqrtxz2i)/16. + (NC*x* 
      li2spec7*poly2i* 
      sqrtxz2i)/16. - (NC*x* 
      li2spec8*poly2i* 
      sqrtxz2i)/16. + (NC*x*lx*lspec5*poly2i* 
      sqrtxz2i)/16. - (NC*x*lx*lspec6*poly2i* 
      sqrtxz2i)/16. + (x*li2spec5*NCi* 
      poly2i*sqrtxz2i)/16. -  
   (x*li2spec6*NCi*poly2i*sqrtxz2i)/ 
    16. - (x*li2spec7*NCi* 
      poly2i*sqrtxz2i)/16. +  
   (x*li2spec8*NCi* 
      poly2i*sqrtxz2i)/16. -  
   (x*lx*lspec5*NCi*poly2i*sqrtxz2i)/ 
    16. + (x*lx*lspec6*NCi*poly2i* 
      sqrtxz2i)/16. + (3*NC*Li2x*omxi)/2. -  
   3*NC*z*Li2x*omxi -  
   (3*NC*li2spec1*omxi)/2. -  
   NC*z*li2spec1*omxi +  
   (3*NC*li2spec2*omxi)/2. +  
   NC*z*li2spec2*omxi +  
   (NC*Li2z*omxi)/2. - NC*z*Li2z*omxi +  
   NC*li2spec19*omxi + 2*NC*z*li2spec19*omxi -  
   NC*li2spec20*omxi -  
   2*NC*z*li2spec20*omxi +  
   (NC*li2spec3*omxi)/2. +  
   3*NC*z*li2spec3*omxi -  
   (NC*li2spec4*omxi)/2. -  
   3*NC*z*li2spec4*omxi -  
   (NC*li2omxzi*omxi)/2. +  
   NC*z*li2omxzi*omxi +  
   3*NC*sqrtxz1*l2*omxi - NC*lx*omxi +  
   (3*NC*sqrtxz1*lx*omxi)/2. +  
   (3*NC*z*lx*omxi)/2. + (5*NC*l2*lx*omxi)/2. -  
   NC*z*l2*lx*omxi +  
   (3*NC*lomx*lx*omxi)/2. -  
   3*NC*z*lomx*lx*omxi -  
   3*NC*sqrtxz1*lspec1*omxi -  
   4*NC*l2*lspec1*omxi -  
   NC*lx*lspec1*omxi +  
   2*NC*z*lx*lspec1*omxi -  
   (3*NC*lz*omxi)/2. + (3*NC*sqrtxz1*lz*omxi)/2. -  
   (3*NC*z*lz*omxi)/2. + (7*NC*l2*lz*omxi)/2. +  
   5*NC*z*l2*lz*omxi + 2*NC*lx*lz*omxi -  
   2*NC*z*lx*lz*omxi +  
   (NC*lomz*lz*omxi)/2. -  
   NC*z*lomz*lz*omxi -  
   (3*NC*lspec1*lz*omxi)/2. -  
   NC*z*lspec1*lz*omxi -  
   2*NC*l2*lspec2*omxi -  
   4*NC*z*l2*lspec2*omxi -  
   (3*NC*lx*lspec2*omxi)/2. -  
   NC*z*lx*lspec2*omxi +  
   2*NC*lspec1*lspec2*omxi +  
   4*NC*z*lspec1*lspec2*omxi -  
   2*NC*lz*lspec2*omxi -  
   4*NC*z*lz*lspec2*omxi -  
   NC*lx*lxpz*omxi -  
   2*NC*z*lx*lxpz*omxi +  
   NC*lz*lxpz*omxi +  
   2*NC*z*lz*lxpz*omxi +  
   NC*lx*lopxz*omxi +  
   2*NC*z*lx*lopxz*omxi +  
   NC*lz*lopxz*omxi +  
   2*NC*z*lz*lopxz*omxi -  
   (3*Li2x*NCi*omxi)/2. +  
   3*z*Li2x*NCi*omxi +  
   (3*li2spec1*NCi*omxi)/2. +  
   z*li2spec1*NCi*omxi -  
   (3*li2spec2*NCi*omxi)/2. -  
   z*li2spec2*NCi*omxi -  
   (Li2z*NCi*omxi)/2. + z*Li2z*NCi*omxi -  
   li2spec19*NCi*omxi -  
   2*z*li2spec19*NCi*omxi +  
   li2spec20*NCi*omxi +  
   2*z*li2spec20*NCi*omxi -  
   (li2spec3*NCi* 
      omxi)/2. - 3*z*li2spec3*NCi*omxi +  
   (li2spec4*NCi* 
      omxi)/2. + 3*z*li2spec4*NCi*omxi +  
   (li2omxzi*NCi*omxi)/2. -  
   z*li2omxzi*NCi*omxi -  
   3*sqrtxz1*l2*NCi*omxi +  
   lx*NCi*omxi -  
   (3*sqrtxz1*lx*NCi*omxi)/2. -  
   (3*z*lx*NCi*omxi)/2. -  
   (5*l2*lx*NCi*omxi)/2. +  
   z*l2*lx*NCi*omxi -  
   (3*lomx*lx*NCi*omxi)/2. +  
   3*z*lomx*lx*NCi*omxi +  
   3*sqrtxz1*lspec1*NCi*omxi +  
   4*l2*lspec1*NCi*omxi +  
   lx*lspec1*NCi*omxi -  
   2*z*lx*lspec1*NCi*omxi +  
   (3*lz*NCi*omxi)/2. -  
   (3*sqrtxz1*lz*NCi*omxi)/2. +  
   (3*z*lz*NCi*omxi)/2. -  
   (7*l2*lz*NCi*omxi)/2. -  
   5*z*l2*lz*NCi*omxi -  
   2*lx*lz*NCi*omxi +  
   2*z*lx*lz*NCi*omxi -  
   (lomz*lz*NCi*omxi)/2. +  
   z*lomz*lz*NCi*omxi +  
   (3*lspec1*lz*NCi*omxi)/2. +  
   z*lspec1*lz*NCi*omxi +  
   2*l2*lspec2*NCi*omxi +  
   4*z*l2*lspec2*NCi*omxi +  
   (3*lx*lspec2*NCi*omxi)/2. +  
   z*lx*lspec2*NCi*omxi -  
   2*lspec1*lspec2*NCi*omxi -  
   4*z*lspec1*lspec2*NCi*omxi +  
   2*lz*lspec2*NCi*omxi +  
   4*z*lz*lspec2*NCi*omxi +  
   lx*lxpz*NCi*omxi +  
   2*z*lx*lxpz*NCi*omxi -  
   lz*lxpz*NCi*omxi -  
   2*z*lz*lxpz*NCi*omxi -  
   lx*lopxz*NCi*omxi -  
   2*z*lx*lopxz*NCi*omxi -  
   lz*lopxz*NCi*omxi -  
   2*z*lz*lopxz*NCi*omxi -  
   (5*NC*pi2*omxi)/12. + (5*NC*z*pi2*omxi)/6. +  
   (5*NCi*pi2*omxi)/12. -  
   (5*z*NCi*pi2*omxi)/6. - NC*Li2mx*xi2 +  
   2*NC*z*Li2mx*xi2 + NC*Li2x*xi2 - 2*NC*z*Li2x*xi2 +  
   (NC*li2spec19*xi2)/2. - NC*z*li2spec19*xi2 +  
   (NC*li2spec20*xi2)/2. -  
   NC*z*li2spec20*xi2 - 2*NC*lx*xi2 +  
   4*NC*z*lx*xi2 + NC*lomx*lx*xi2 -  
   2*NC*z*lomx*lx*xi2 - NC*lx*lopx*xi2 +  
   2*NC*z*lx*lopx*xi2 - (NC*lx*lz*xi2)/2. +  
   NC*z*lx*lz*xi2 + (NC*lx*lopxz*xi2)/2. -  
   NC*z*lx*lopxz*xi2 +  
   (NC*lz*lopxz*xi2)/2. -  
   NC*z*lz*lopxz*xi2 +  
   (NC*lx*lopxzi*xi2)/2. -  
   NC*z*lx*lopxzi*xi2 -  
   (NC*lz*lopxzi*xi2)/2. +  
   NC*z*lz*lopxzi*xi2 +  
   Li2mx*NCi*xi2 - 2*z*Li2mx*NCi*xi2 -  
   Li2x*NCi*xi2 + 2*z*Li2x*NCi*xi2 -  
   (li2spec19*NCi*xi2)/2. +  
   z*li2spec19*NCi*xi2 -  
   (li2spec20*NCi*xi2)/2. +  
   z*li2spec20*NCi*xi2 +  
   2*lx*NCi*xi2 - 4*z*lx*NCi*xi2 -  
   lomx*lx*NCi*xi2 +  
   2*z*lomx*lx*NCi*xi2 +  
   lx*lopx*NCi*xi2 -  
   2*z*lx*lopx*NCi*xi2 +  
   (lx*lz*NCi*xi2)/2. -  
   z*lx*lz*NCi*xi2 -  
   (lx*lopxz*NCi*xi2)/2. +  
   z*lx*lopxz*NCi*xi2 -  
   (lz*lopxz*NCi*xi2)/2. +  
   z*lz*lopxz*NCi*xi2 -  
   (lx*lopxzi*NCi*xi2)/2. +  
   z*lx*lopxzi*NCi*xi2 +  
   (lz*lopxzi*NCi*xi2)/2. -  
   z*lz*lopxzi*NCi*xi2 -  
   (NC*pi2*xi2)/6. + (NC*z*pi2*xi2)/3. +  
   (NCi*pi2*xi2)/6. -  
   (z*NCi*pi2*xi2)/3. -  
   (3*NC*sqrtxz3*z*itani1*xi)/4. +  
   (3*NC*sqrtxz3*z*itani2*xi)/4. -  
   (3*NC*sqrtxz3*z*itani3*xi)/2. +  
   NC*Li2mx*xi - 2*NC*z*Li2mx*xi - NC*Li2x*xi +  
   2*NC*z*Li2x*xi - (NC*li2spec19*xi)/2. +  
   NC*z*li2spec19*xi - (NC*li2spec20*xi)/2. +  
   NC*z*li2spec20*xi -  
   (3*NC*sqrtxz3*z*atanspec1*lspec3*xi)/2. +  
   (17*NC*lx*xi)/8. - 4*NC*z*lx*xi -  
   NC*lomx*lx*xi + 2*NC*z*lomx*lx*xi +  
   NC*lx*lopx*xi - 2*NC*z*lx*lopx*xi +  
   (NC*lz*xi)/8. + (NC*lx*lz*xi)/2. -  
   NC*z*lx*lz*xi +  
   (3*NC*sqrtxz3*z*atanspec2*lspec4*xi)/2. -  
   (NC*lx*lopxz*xi)/2. +  
   NC*z*lx*lopxz*xi -  
   (NC*lz*lopxz*xi)/2. +  
   NC*z*lz*lopxz*xi -  
   (NC*lx*lopxzi*xi)/2. +  
   NC*z*lx*lopxzi*xi +  
   (NC*lz*lopxzi*xi)/2. -  
   NC*z*lz*lopxzi*xi +  
   (3*sqrtxz3*z*itani1*NCi*xi)/4. -  
   (3*sqrtxz3*z*itani2*NCi*xi)/4. +  
   (3*sqrtxz3*z*itani3*NCi*xi)/2. -  
   Li2mx*NCi*xi + 2*z*Li2mx*NCi*xi +  
   Li2x*NCi*xi - 2*z*Li2x*NCi*xi +  
   (li2spec19*NCi*xi)/2. -  
   z*li2spec19*NCi*xi +  
   (li2spec20*NCi*xi)/2. -  
   z*li2spec20*NCi*xi +  
   (3*sqrtxz3*z*atanspec1*lspec3*NCi*xi)/2. -  
   (17*lx*NCi*xi)/8. + 4*z*lx*NCi*xi +  
   lomx*lx*NCi*xi -  
   2*z*lomx*lx*NCi*xi -  
   lx*lopx*NCi*xi +  
   2*z*lx*lopx*NCi*xi -  
   (lz*NCi*xi)/8. -  
   (lx*lz*NCi*xi)/2. +  
   z*lx*lz*NCi*xi -  
   (3*sqrtxz3*z*atanspec2*lspec4*NCi*xi)/2. +  
   (lx*lopxz*NCi*xi)/2. -  
   z*lx*lopxz*NCi*xi +  
   (lz*lopxz*NCi*xi)/2. -  
   z*lz*lopxz*NCi*xi +  
   (lx*lopxzi*NCi*xi)/2. -  
   z*lx*lopxzi*NCi*xi -  
   (lz*lopxzi*NCi*xi)/2. +  
   z*lz*lopxzi*NCi*xi +  
   (NC*pi2*xi)/6. - (NC*z*pi2*xi)/3. -  
   (NCi*pi2*xi)/6. +  
   (z*NCi*pi2*xi)/3. -  
   (NC*lx*poly2i*xi)/8. -  
   (NC*lz*poly2i*xi)/8. +  
   (lx*NCi*poly2i*xi)/8. +  
   (lz*NCi*poly2i*xi)/8. -  
   (NC*li2spec5*sqrtxz2i*xi)/16. +  
   (NC*li2spec6*sqrtxz2i*xi)/16. +  
   (NC*li2spec7*sqrtxz2i* 
      xi)/16. - (NC*li2spec8* 
      sqrtxz2i*xi)/16. +  
   (NC*lx*lspec5*sqrtxz2i*xi)/16. -  
   (NC*lx*lspec6*sqrtxz2i*xi)/16. +  
   (li2spec5*NCi*sqrtxz2i*xi)/16. -  
   (li2spec6*NCi*sqrtxz2i*xi)/16. -  
   (li2spec7*NCi* 
      sqrtxz2i*xi)/16. +  
   (li2spec8*NCi* 
      sqrtxz2i*xi)/16. -  
   (lx*lspec5*NCi*sqrtxz2i*xi)/16. +  
   (lx*lspec6*NCi*sqrtxz2i*xi)/16. +  
   (NC*li2spec5*poly2i*sqrtxz2i*xi)/ 
    16. - (NC*li2spec6*poly2i*sqrtxz2i* 
      xi)/16. - (NC*li2spec7* 
      poly2i*sqrtxz2i*xi)/16. +  
   (NC*li2spec8*poly2i* 
      sqrtxz2i*xi)/16. -  
   (NC*lx*lspec5*poly2i*sqrtxz2i*xi)/ 
    16. + (NC*lx*lspec6*poly2i*sqrtxz2i* 
      xi)/16. - (li2spec5*NCi*poly2i* 
      sqrtxz2i*xi)/16. +  
   (li2spec6*NCi*poly2i*sqrtxz2i* 
      xi)/16. + (li2spec7* 
      NCi*poly2i*sqrtxz2i*xi)/16. -  
   (li2spec8*NCi*poly2i* 
      sqrtxz2i*xi)/16. +  
   (lx*lspec5*NCi*poly2i*sqrtxz2i* 
      xi)/16. - (lx*lspec6*NCi*poly2i* 
      sqrtxz2i*xi)/16. - (7*NC*lx*x2)/8. +  
   (7*NC*lz*x2)/8. + (7*lx*NCi*x2)/8. -  
   (7*lz*NCi*x2)/8. +  
   (9*NC*li2spec5*sqrtxz2i*x2)/8. -  
   (9*NC*z*li2spec5*sqrtxz2i*x2)/4. -  
   (9*NC*li2spec6*sqrtxz2i*x2)/8. +  
   (9*NC*z*li2spec6*sqrtxz2i*x2)/4. -  
   (9*NC*li2spec7*sqrtxz2i* 
      x2)/8. + (9*NC*z*li2spec7* 
      sqrtxz2i*x2)/4. +  
   (9*NC*li2spec8*sqrtxz2i* 
      x2)/8. - (9*NC*z*li2spec8* 
      sqrtxz2i*x2)/4. -  
   (9*NC*lx*lspec5*sqrtxz2i*x2)/8. +  
   (9*NC*z*lx*lspec5*sqrtxz2i*x2)/4. +  
   (9*NC*lx*lspec6*sqrtxz2i*x2)/8. -  
   (9*NC*z*lx*lspec6*sqrtxz2i*x2)/4. -  
   (9*li2spec5*NCi*sqrtxz2i*x2)/8. +  
   (9*z*li2spec5*NCi*sqrtxz2i*x2)/ 
    4. + (9*li2spec6*NCi*sqrtxz2i*x2)/ 
    8. - (9*z*li2spec6*NCi*sqrtxz2i* 
      x2)/4. + (9*li2spec7* 
      NCi*sqrtxz2i*x2)/8. -  
   (9*z*li2spec7*NCi* 
      sqrtxz2i*x2)/4. -  
   (9*li2spec8*NCi* 
      sqrtxz2i*x2)/8. +  
   (9*z*li2spec8*NCi* 
      sqrtxz2i*x2)/4. +  
   (9*lx*lspec5*NCi*sqrtxz2i*x2)/8. -  
   (9*z*lx*lspec5*NCi*sqrtxz2i*x2)/4. -  
   (9*lx*lspec6*NCi*sqrtxz2i*x2)/8. +  
   (9*z*lx*lspec6*NCi*sqrtxz2i*x2)/4. +  
   (NC*lx*poly2i*x3)/8. +  
   (NC*lz*poly2i*x3)/8. -  
   (lx*NCi*poly2i*x3)/8. -  
   (lz*NCi*poly2i*x3)/8. -  
   (NC*li2spec5*sqrtxz2i*x3)/16. +  
   (NC*li2spec6*sqrtxz2i*x3)/16. +  
   (NC*li2spec7*sqrtxz2i* 
      x3)/16. - (NC*li2spec8* 
      sqrtxz2i*x3)/16. +  
   (NC*lx*lspec5*sqrtxz2i*x3)/16. -  
   (NC*lx*lspec6*sqrtxz2i*x3)/16. +  
   (li2spec5*NCi*sqrtxz2i*x3)/16. -  
   (li2spec6*NCi*sqrtxz2i*x3)/16. -  
   (li2spec7*NCi* 
      sqrtxz2i*x3)/16. +  
   (li2spec8*NCi* 
      sqrtxz2i*x3)/16. -  
   (lx*lspec5*NCi*sqrtxz2i*x3)/16. +  
   (lx*lspec6*NCi*sqrtxz2i*x3)/16. -  
   (NC*li2spec5*poly2i*sqrtxz2i*x3)/ 
    16. + (NC*li2spec6*poly2i*sqrtxz2i* 
      x3)/16. + (NC*li2spec7* 
      poly2i*sqrtxz2i*x3)/16. -  
   (NC*li2spec8*poly2i* 
      sqrtxz2i*x3)/16. +  
   (NC*lx*lspec5*poly2i*sqrtxz2i*x3)/ 
    16. - (NC*lx*lspec6*poly2i*sqrtxz2i* 
      x3)/16. + (li2spec5*NCi*poly2i* 
      sqrtxz2i*x3)/16. -  
   (li2spec6*NCi*poly2i*sqrtxz2i* 
      x3)/16. - (li2spec7* 
      NCi*poly2i*sqrtxz2i*x3)/16. +  
   (li2spec8*NCi*poly2i* 
      sqrtxz2i*x3)/16. -  
   (lx*lspec5*NCi*poly2i*sqrtxz2i* 
      x3)/16. + (lx*lspec6*NCi*poly2i* 
      sqrtxz2i*x3)/16. - (NC*lx*poly2i*x4)/8. +  
   (NC*lz*poly2i*x4)/8. +  
   (lx*NCi*poly2i*x4)/8. -  
   (lz*NCi*poly2i*x4)/8. +  
   (NC*li2spec5*poly2i*sqrtxz2i*x5)/ 
    16. - (NC*li2spec6*poly2i*sqrtxz2i* 
      x5)/16. - (NC*li2spec7* 
      poly2i*sqrtxz2i*x5)/16. +  
   (NC*li2spec8*poly2i* 
      sqrtxz2i*x5)/16. -  
   (NC*lx*lspec5*poly2i*sqrtxz2i*x5)/ 
    16. + (NC*lx*lspec6*poly2i*sqrtxz2i* 
      x5)/16. - (li2spec5*NCi*poly2i* 
      sqrtxz2i*x5)/16. +  
   (li2spec6*NCi*poly2i*sqrtxz2i* 
      x5)/16. + (li2spec7* 
      NCi*poly2i*sqrtxz2i*x5)/16. -  
   (li2spec8*NCi*poly2i* 
      sqrtxz2i*x5)/16. +  
   (lx*lspec5*NCi*poly2i*sqrtxz2i* 
      x5)/16. - (lx*lspec6*NCi*poly2i* 
      sqrtxz2i*x5)/16. - NC*Li2mx*opxi +  
   2*NC*z*Li2mx*opxi + NC*Li2x*opxi -  
   2*NC*z*Li2x*opxi - (NC*li2spec19*opxi)/2. +  
   NC*z*li2spec19*opxi -  
   (NC*li2spec20*opxi)/2. +  
   NC*z*li2spec20*opxi - 2*NC*lx*opxi +  
   4*NC*z*lx*opxi + NC*lomx*lx*opxi -  
   2*NC*z*lomx*lx*opxi -  
   NC*lx*lopx*opxi +  
   2*NC*z*lx*lopx*opxi +  
   (NC*lx*lz*opxi)/2. - NC*z*lx*lz*opxi -  
   (NC*lx*lopxz*opxi)/2. +  
   NC*z*lx*lopxz*opxi -  
   (NC*lz*lopxz*opxi)/2. +  
   NC*z*lz*lopxz*opxi -  
   (NC*lx*lopxzi*opxi)/2. +  
   NC*z*lx*lopxzi*opxi +  
   (NC*lz*lopxzi*opxi)/2. -  
   NC*z*lz*lopxzi*opxi +  
   Li2mx*NCi*opxi - 2*z*Li2mx*NCi*opxi -  
   Li2x*NCi*opxi + 2*z*Li2x*NCi*opxi +  
   (li2spec19*NCi*opxi)/2. -  
   z*li2spec19*NCi*opxi +  
   (li2spec20*NCi*opxi)/2. -  
   z*li2spec20*NCi*opxi +  
   2*lx*NCi*opxi - 4*z*lx*NCi*opxi -  
   lomx*lx*NCi*opxi +  
   2*z*lomx*lx*NCi*opxi +  
   lx*lopx*NCi*opxi -  
   2*z*lx*lopx*NCi*opxi -  
   (lx*lz*NCi*opxi)/2. +  
   z*lx*lz*NCi*opxi +  
   (lx*lopxz*NCi*opxi)/2. -  
   z*lx*lopxz*NCi*opxi +  
   (lz*lopxz*NCi*opxi)/2. -  
   z*lz*lopxz*NCi*opxi +  
   (lx*lopxzi*NCi*opxi)/2. -  
   z*lx*lopxzi*NCi*opxi -  
   (lz*lopxzi*NCi*opxi)/2. +  
   z*lz*lopxzi*NCi*opxi -  
   (NC*pi2*opxi)/3. + (2*NC*z*pi2*opxi)/3. +  
   (NCi*pi2*opxi)/3. -  
   (2*z*NCi*pi2*opxi)/3. +  
   NC*Li2mx*xi2*opxi -  
   2*NC*z*Li2mx*xi2*opxi -  
   NC*Li2x*xi2*opxi +  
   2*NC*z*Li2x*xi2*opxi -  
   (NC*li2spec19*xi2*opxi)/2. +  
   NC*z*li2spec19*xi2*opxi -  
   (NC*li2spec20*xi2*opxi)/2. +  
   NC*z*li2spec20*xi2*opxi +  
   2*NC*lx*xi2*opxi -  
   4*NC*z*lx*xi2*opxi -  
   NC*lomx*lx*xi2*opxi +  
   2*NC*z*lomx*lx*xi2*opxi +  
   NC*lx*lopx*xi2*opxi -  
   2*NC*z*lx*lopx*xi2*opxi +  
   (NC*lx*lz*xi2*opxi)/2. -  
   NC*z*lx*lz*xi2*opxi -  
   (NC*lx*lopxz*xi2*opxi)/2. +  
   NC*z*lx*lopxz*xi2*opxi -  
   (NC*lz*lopxz*xi2*opxi)/2. +  
   NC*z*lz*lopxz*xi2*opxi -  
   (NC*lx*lopxzi*xi2*opxi)/2. +  
   NC*z*lx*lopxzi*xi2*opxi +  
   (NC*lz*lopxzi*xi2*opxi)/2. -  
   NC*z*lz*lopxzi*xi2*opxi -  
   Li2mx*NCi*xi2*opxi +  
   2*z*Li2mx*NCi*xi2*opxi +  
   Li2x*NCi*xi2*opxi -  
   2*z*Li2x*NCi*xi2*opxi +  
   (li2spec19*NCi*xi2*opxi)/2. -  
   z*li2spec19*NCi*xi2*opxi +  
   (li2spec20*NCi*xi2*opxi)/2. -  
   z*li2spec20*NCi*xi2*opxi -  
   2*lx*NCi*xi2*opxi +  
   4*z*lx*NCi*xi2*opxi +  
   lomx*lx*NCi*xi2*opxi -  
   2*z*lomx*lx*NCi*xi2*opxi -  
   lx*lopx*NCi*xi2*opxi +  
   2*z*lx*lopx*NCi*xi2*opxi -  
   (lx*lz*NCi*xi2*opxi)/2. +  
   z*lx*lz*NCi*xi2*opxi +  
   (lx*lopxz*NCi*xi2*opxi)/2. -  
   z*lx*lopxz*NCi*xi2*opxi +  
   (lz*lopxz*NCi*xi2*opxi)/2. -  
   z*lz*lopxz*NCi*xi2*opxi +  
   (lx*lopxzi*NCi*xi2*opxi)/2. -  
   z*lx*lopxzi*NCi*xi2*opxi -  
   (lz*lopxzi*NCi*xi2*opxi)/2. +  
   z*lz*lopxzi*NCi*xi2*opxi +  
   (NC*pi2*xi2*opxi)/6. -  
   (NC*z*pi2*xi2*opxi)/3. -  
   (NCi*pi2*xi2*opxi)/6. +  
   (z*NCi*pi2*xi2*opxi)/3. +  
   (NC*x*lx*xmzi)/2. - (NC*x*lz*xmzi)/2. -  
   (x*lx*NCi*xmzi)/2. +  
   (x*lz*NCi*xmzi)/2. -  
   NC*lx*x2*xmzi + NC*lz*x2*xmzi +  
   lx*NCi*x2*xmzi -  
   lz*NCi*x2*xmzi +  
   NC*lx*x3*xmzi - NC*lz*x3*xmzi -  
   lx*NCi*x3*xmzi +  
   lz*NCi*x3*xmzi +  
   (23*NC*sqrtxz3*itani1*z2)/4. -  
   (23*NC*sqrtxz3*itani2*z2)/4. +  
   (23*NC*sqrtxz3*itani3*z2)/2. -  
   4*NC*x*Li2x*z2 + 4*NC*x*li2spec1*z2 -  
   4*NC*x*li2spec2*z2 -  
   4*NC*x*li2spec19*z2 +  
   (23*NC*sqrtxz3*atanspec1*lspec3*z2)/2. -  
   8*NC*x*l2*lx*z2 - 4*NC*x*lomx*lx*z2 +  
   12*NC*x*l2*lspec1*z2 +  
   4*NC*x*lx*lspec1*z2 -  
   8*NC*x*l2*lz*z2 - 2*NC*x*lx*lz*z2 +  
   4*NC*x*lspec1*lz*z2 -  
   (23*NC*sqrtxz3*atanspec2*lspec4*z2)/2. +  
   4*NC*x*l2*lspec2*z2 +  
   4*NC*x*lx*lspec2*z2 -  
   4*NC*x*lspec1*lspec2*z2 +  
   4*NC*x*lz*lspec2*z2 +  
   2*NC*x*lx*lxpz*z2 - 2*NC*x*lz*lxpz*z2 -  
   4*NC*x*lx*lopxz*z2 -  
   4*NC*x*lz*lopxz*z2 -  
   2*NC*x*lx*lopxzi*z2 +  
   2*NC*x*lz*lopxzi*z2 -  
   (23*sqrtxz3*itani1*NCi*z2)/4. +  
   (23*sqrtxz3*itani2*NCi*z2)/4. -  
   (23*sqrtxz3*itani3*NCi*z2)/2. +  
   4*x*Li2x*NCi*z2 -  
   4*x*li2spec1*NCi*z2 +  
   4*x*li2spec2*NCi*z2 +  
   4*x*li2spec19*NCi*z2 -  
   (23*sqrtxz3*atanspec1*lspec3*NCi*z2)/2. +  
   8*x*l2*lx*NCi*z2 +  
   4*x*lomx*lx*NCi*z2 -  
   12*x*l2*lspec1*NCi*z2 -  
   4*x*lx*lspec1*NCi*z2 +  
   8*x*l2*lz*NCi*z2 +  
   2*x*lx*lz*NCi*z2 -  
   4*x*lspec1*lz*NCi*z2 +  
   (23*sqrtxz3*atanspec2*lspec4*NCi*z2)/2. -  
   4*x*l2*lspec2*NCi*z2 -  
   4*x*lx*lspec2*NCi*z2 +  
   4*x*lspec1*lspec2*NCi*z2 -  
   4*x*lz*lspec2*NCi*z2 -  
   2*x*lx*lxpz*NCi*z2 +  
   2*x*lz*lxpz*NCi*z2 +  
   4*x*lx*lopxz*NCi*z2 +  
   4*x*lz*lopxz*NCi*z2 +  
   2*x*lx*lopxzi*NCi*z2 -  
   2*x*lz*lopxzi*NCi*z2 +  
   (2*NC*x*pi2*z2)/3. - (2*x*NCi*pi2*z2)/3. +  
   3*NC*x*li2spec5*sqrtxz2i*z2 -  
   3*NC*x*li2spec6*sqrtxz2i*z2 -  
   3*NC*x*li2spec7*sqrtxz2i* 
    z2 + 3*NC*x*li2spec8* 
    sqrtxz2i*z2 - 3*NC*x*lx*lspec5* 
    sqrtxz2i*z2 + 3*NC*x*lx*lspec6* 
    sqrtxz2i*z2 - 3*x*li2spec5*NCi* 
    sqrtxz2i*z2 + 3*x*li2spec6*NCi* 
    sqrtxz2i*z2 + 3*x* 
    li2spec7*NCi* 
    sqrtxz2i*z2 - 3*x* 
    li2spec8*NCi* 
    sqrtxz2i*z2 + 3*x*lx*lspec5*NCi* 
    sqrtxz2i*z2 - 3*x*lx*lspec6*NCi* 
    sqrtxz2i*z2 + 4*NC*Li2x*omxi*z2 -  
   4*NC*li2spec1*omxi*z2 +  
   4*NC*li2spec2*omxi*z2 +  
   2*NC*li2spec19*omxi*z2 -  
   2*NC*li2spec20*omxi*z2 +  
   8*NC*l2*lx*omxi*z2 +  
   4*NC*lomx*lx*omxi*z2 -  
   12*NC*l2*lspec1*omxi*z2 -  
   4*NC*lx*lspec1*omxi*z2 +  
   8*NC*l2*lz*omxi*z2 +  
   4*NC*lx*lz*omxi*z2 -  
   4*NC*lspec1*lz*omxi*z2 -  
   4*NC*l2*lspec2*omxi*z2 -  
   4*NC*lx*lspec2*omxi*z2 +  
   4*NC*lspec1*lspec2*omxi*z2 -  
   4*NC*lz*lspec2*omxi*z2 -  
   2*NC*lx*lxpz*omxi*z2 +  
   2*NC*lz*lxpz*omxi*z2 +  
   2*NC*lx*lopxz*omxi*z2 +  
   2*NC*lz*lopxz*omxi*z2 -  
   4*Li2x*NCi*omxi*z2 +  
   4*li2spec1*NCi*omxi*z2 -  
   4*li2spec2*NCi*omxi*z2 -  
   2*li2spec19*NCi*omxi*z2 +  
   2*li2spec20*NCi*omxi*z2 -  
   8*l2*lx*NCi*omxi*z2 -  
   4*lomx*lx*NCi*omxi*z2 +  
   12*l2*lspec1*NCi*omxi*z2 +  
   4*lx*lspec1*NCi*omxi*z2 -  
   8*l2*lz*NCi*omxi*z2 -  
   4*lx*lz*NCi*omxi*z2 +  
   4*lspec1*lz*NCi*omxi*z2 +  
   4*l2*lspec2*NCi*omxi*z2 +  
   4*lx*lspec2*NCi*omxi*z2 -  
   4*lspec1*lspec2*NCi*omxi* 
    z2 + 4*lz*lspec2*NCi*omxi* 
    z2 + 2*lx*lxpz*NCi*omxi*z2 -  
   2*lz*lxpz*NCi*omxi*z2 -  
   2*lx*lopxz*NCi*omxi*z2 -  
   2*lz*lopxz*NCi*omxi*z2 -  
   NC*pi2*omxi*z2 +  
   NCi*pi2*omxi*z2 -  
   4*NC*Li2mx*xi2*z2 + 4*NC*Li2x*xi2*z2 +  
   2*NC*li2spec19*xi2*z2 +  
   2*NC*li2spec20*xi2*z2 -  
   8*NC*lx*xi2*z2 +  
   4*NC*lomx*lx*xi2*z2 -  
   4*NC*lx*lopx*xi2*z2 -  
   2*NC*lx*lz*xi2*z2 +  
   2*NC*lx*lopxz*xi2*z2 +  
   2*NC*lz*lopxz*xi2*z2 +  
   2*NC*lx*lopxzi*xi2*z2 -  
   2*NC*lz*lopxzi*xi2*z2 +  
   4*Li2mx*NCi*xi2*z2 -  
   4*Li2x*NCi*xi2*z2 -  
   2*li2spec19*NCi*xi2*z2 -  
   2*li2spec20*NCi*xi2*z2 +  
   8*lx*NCi*xi2*z2 -  
   4*lomx*lx*NCi*xi2*z2 +  
   4*lx*lopx*NCi*xi2*z2 +  
   2*lx*lz*NCi*xi2*z2 -  
   2*lx*lopxz*NCi*xi2*z2 -  
   2*lz*lopxz*NCi*xi2*z2 -  
   2*lx*lopxzi*NCi*xi2*z2 +  
   2*lz*lopxzi*NCi*xi2*z2 -  
   (2*NC*pi2*xi2*z2)/3. +  
   (2*NCi*pi2*xi2*z2)/3. -  
   2*NC*li2spec19*opxi*z2 -  
   2*NC*li2spec20*opxi*z2 +  
   2*NC*lx*lz*opxi*z2 -  
   2*NC*lx*lopxz*opxi*z2 -  
   2*NC*lz*lopxz*opxi*z2 -  
   2*NC*lx*lopxzi*opxi*z2 +  
   2*NC*lz*lopxzi*opxi*z2 +  
   2*li2spec19*NCi*opxi*z2 +  
   2*li2spec20*NCi*opxi*z2 -  
   2*lx*lz*NCi*opxi*z2 +  
   2*lx*lopxz*NCi*opxi*z2 +  
   2*lz*lopxz*NCi*opxi*z2 +  
   2*lx*lopxzi*NCi*opxi*z2 -  
   2*lz*lopxzi*NCi*opxi*z2 -  
   (NC*pi2*opxi*z2)/3. +  
   (NCi*pi2*opxi*z2)/3. +  
   4*NC*Li2mx*xi2*opxi*z2 -  
   4*NC*Li2x*xi2*opxi*z2 -  
   2*NC*li2spec19*xi2*opxi*z2 -  
   2*NC*li2spec20*xi2*opxi*z2 +  
   8*NC*lx*xi2*opxi*z2 -  
   4*NC*lomx*lx*xi2*opxi*z2 +  
   4*NC*lx*lopx*xi2*opxi*z2 +  
   2*NC*lx*lz*xi2*opxi*z2 -  
   2*NC*lx*lopxz*xi2*opxi*z2 -  
   2*NC*lz*lopxz*xi2*opxi*z2 -  
   2*NC*lx*lopxzi*xi2*opxi*z2 +  
   2*NC*lz*lopxzi*xi2*opxi*z2 -  
   4*Li2mx*NCi*xi2*opxi*z2 +  
   4*Li2x*NCi*xi2*opxi*z2 +  
   2*li2spec19*NCi*xi2*opxi*z2 +  
   2*li2spec20*NCi*xi2*opxi*z2 -  
   8*lx*NCi*xi2*opxi*z2 +  
   4*lomx*lx*NCi*xi2*opxi*z2 -  
   4*lx*lopx*NCi*xi2*opxi*z2 -  
   2*lx*lz*NCi*xi2*opxi*z2 +  
   2*lx*lopxz*NCi*xi2*opxi*z2 +  
   2*lz*lopxz*NCi*xi2*opxi*z2 +  
   2*lx*lopxzi*NCi*xi2*opxi* 
    z2 - 2*lz*lopxzi*NCi*xi2* 
    opxi*z2 + (2*NC*pi2*xi2*opxi* 
      z2)/3. - (2*NCi*pi2*xi2*opxi* 
      z2)/3. + 4*NC*Li2mx*xi*opxi*z2 -  
   4*NC*Li2x*xi*opxi*z2 -  
   2*NC*li2spec19*xi*opxi*z2 -  
   2*NC*li2spec20*xi*opxi*z2 +  
   8*NC*lx*xi*opxi*z2 -  
   4*NC*lomx*lx*xi*opxi*z2 +  
   4*NC*lx*lopx*xi*opxi*z2 +  
   2*NC*lx*lz*xi*opxi*z2 -  
   2*NC*lx*lopxz*xi*opxi*z2 -  
   2*NC*lz*lopxz*xi*opxi*z2 -  
   2*NC*lx*lopxzi*xi*opxi*z2 +  
   2*NC*lz*lopxzi*xi*opxi*z2 -  
   4*Li2mx*NCi*xi*opxi*z2 +  
   4*Li2x*NCi*xi*opxi*z2 +  
   2*li2spec19*NCi*xi*opxi*z2 +  
   2*li2spec20*NCi*xi*opxi*z2 -  
   8*lx*NCi*xi*opxi*z2 +  
   4*lomx*lx*NCi*xi*opxi*z2 -  
   4*lx*lopx*NCi*xi*opxi*z2 -  
   2*lx*lz*NCi*xi*opxi*z2 +  
   2*lx*lopxz*NCi*xi*opxi*z2 +  
   2*lz*lopxz*NCi*xi*opxi*z2 +  
   2*lx*lopxzi*NCi*xi*opxi* 
    z2 - 2*lz*lopxzi*NCi*xi* 
    opxi*z2 + (2*NC*pi2*xi*opxi* 
      z2)/3. - (2*NCi*pi2*xi*opxi* 
      z2)/3. - 2*NC*x*l22 - 4*NC*z*l22 +  
   2*x*NCi*l22 + 4*z*NCi*l22 +  
   3*NC*omxi*l22 + 2*NC*z*omxi*l22 -  
   3*NCi*omxi*l22 -  
   2*z*NCi*omxi*l22 -  
   8*NC*x*z2*l22 + 8*x*NCi*z2*l22 +  
   8*NC*omxi*z2*l22 -  
   8*NCi*omxi*z2*l22 + NC*lx2 +  
   (3*NC*x*lx2)/2. - 2*NC*z*lx2 -  
   3*NC*x*z*lx2 - NCi*lx2 -  
   (3*x*NCi*lx2)/2. + 2*z*NCi*lx2 +  
   3*x*z*NCi*lx2 - NC*omxi*lx2 +  
   2*NC*z*omxi*lx2 +  
   NCi*omxi*lx2 -  
   2*z*NCi*omxi*lx2 +  
   (3*NC*xi2*lx2)/4. - (3*NC*z*xi2*lx2)/2. -  
   (3*NCi*xi2*lx2)/4. +  
   (3*z*NCi*xi2*lx2)/2. -  
   (3*NC*xi*lx2)/4. + (3*NC*z*xi*lx2)/2. +  
   (3*NCi*xi*lx2)/4. -  
   (3*z*NCi*xi*lx2)/2. +  
   (5*NC*opxi*lx2)/4. -  
   (5*NC*z*opxi*lx2)/2. -  
   (5*NCi*opxi*lx2)/4. +  
   (5*z*NCi*opxi*lx2)/2. -  
   (3*NC*xi2*opxi*lx2)/4. +  
   (3*NC*z*xi2*opxi*lx2)/2. +  
   (3*NCi*xi2*opxi*lx2)/4. -  
   (3*z*NCi*xi2*opxi*lx2)/2. +  
   2*NC*x*z2*lx2 - 2*x*NCi*z2*lx2 -  
   NC*omxi*z2*lx2 +  
   NCi*omxi*z2*lx2 +  
   3*NC*xi2*z2*lx2 -  
   3*NCi*xi2*z2*lx2 +  
   NC*opxi*z2*lx2 -  
   NCi*opxi*z2*lx2 -  
   3*NC*xi2*opxi*z2*lx2 +  
   3*NCi*xi2*opxi*z2*lx2 -  
   3*NC*xi*opxi*z2*lx2 +  
   3*NCi*xi*opxi*z2*lx2 +  
   NC*lspec1_2 - NC*x*lspec1_2 -  
   2*NC*z*lspec1_2 +  
   2*NC*x*z*lspec1_2 -  
   NCi*lspec1_2 +  
   x*NCi*lspec1_2 +  
   2*z*NCi*lspec1_2 -  
   2*x*z*NCi*lspec1_2 +  
   NC*omxi*lspec1_2 -  
   2*NC*z*omxi*lspec1_2 -  
   NCi*omxi*lspec1_2 +  
   2*z*NCi*omxi*lspec1_2 -  
   4*NC*x*z2*lspec1_2 +  
   4*x*NCi*z2*lspec1_2 +  
   4*NC*omxi*z2*lspec1_2 -  
   4*NCi*omxi*z2*lspec1_2 -  
   (NC*lz2)/2. + NC*x*lz2 + NC*z*lz2 -  
   2*NC*x*z*lz2 + (NCi*lz2)/2. -  
   x*NCi*lz2 - z*NCi*lz2 +  
   2*x*z*NCi*lz2 - (NC*omxi*lz2)/2. +  
   NC*z*omxi*lz2 +  
   (NCi*omxi*lz2)/2. -  
   z*NCi*omxi*lz2 -  
   (NC*xi2*lz2)/4. + (NC*z*xi2*lz2)/2. +  
   (NCi*xi2*lz2)/4. -  
   (z*NCi*xi2*lz2)/2. +  
   (NC*xi*lz2)/4. - (NC*z*xi*lz2)/2. -  
   (NCi*xi*lz2)/4. +  
   (z*NCi*xi*lz2)/2. +  
   (NC*opxi*lz2)/4. -  
   (NC*z*opxi*lz2)/2. -  
   (NCi*opxi*lz2)/4. +  
   (z*NCi*opxi*lz2)/2. +  
   (NC*xi2*opxi*lz2)/4. -  
   (NC*z*xi2*opxi*lz2)/2. -  
   (NCi*xi2*opxi*lz2)/4. +  
   (z*NCi*xi2*opxi*lz2)/2. +  
   2*NC*x*z2*lz2 - 2*x*NCi*z2*lz2 -  
   NC*omxi*z2*lz2 +  
   NCi*omxi*z2*lz2 -  
   NC*xi2*z2*lz2 +  
   NCi*xi2*z2*lz2 +  
   NC*opxi*z2*lz2 -  
   NCi*opxi*z2*lz2 +  
   NC*xi2*opxi*z2*lz2 -  
   NCi*xi2*opxi*z2*lz2 +  
   NC*xi*opxi*z2*lz2 -  
   NCi*xi*opxi*z2*lz2;
};
