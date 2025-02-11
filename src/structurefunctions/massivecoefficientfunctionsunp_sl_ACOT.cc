#include "apfel/massivecoefficientfunctionsunp_sl_ACOT.h"
#include "apfel/constants.h"
#include "apfel/specialfunctions.h"
#include "apfel/integrator.h"

#include <iostream>

namespace apfel
{
  //_________________________________________________________________________________
  Cm20qNC_ACOT::Cm20qNC_ACOT(double const& eta): Expression(eta,true)
  {
    _xi = 4./(std::pow((2./eta-1),2)-1);
  }

  double Cm20qNC_ACOT::Local(double const& x) const
  {
    if(x>=1)
      {
        return 0;
      }
    return std::sqrt(1+4/_xi);
  }

  //_________________________________________________________________________________
  Cm21gNC_ACOT::Cm21gNC_ACOT(double const& eta):
    Expression(eta,true)
  {
  }
  double Cm21gNC_ACOT::Regular(double const& x) const
  {
    const double eta   = this->_eta;
    const double z     = eta * x;
    const double z2    = z * z;
    const double epsi  = ( 1 - eta ) / eta / 4;
    const double epsi2 = epsi * epsi;
    const double v     = sqrt( 1 - 4 * epsi * z / ( 1 - z ));
    return 4 * TR * ( ( z2 + ( 1 - z ) * ( 1 - z ) + 4 * epsi * z * ( 1 - 3 * z )
                        - 8 * epsi2 * z2 ) * log( ( 1 + v ) / ( 1 - v ) )
                      + ( 8 * z * ( 1 - z ) - 1 - 4 * epsi * z * ( 1 - z ) ) * v );
  }

  //_________________________________________________________________________________
  Cm21gNC_sub_ACOT::Cm21gNC_sub_ACOT(double const& eta): Expression(eta,true)
  {
    _xi = 4/(std::pow((2./eta-1),2)-1);
  }

  double Cm21gNC_sub_ACOT::Regular(double const& z) const
  {
    return 4 * TR *std::sqrt(1+4/_xi)* ( 1 - 2 * z + 2 * z * z );
  }


  //_________________________________________________________________________________
  Cm21qNC_ACOT::Cm21qNC_ACOT(double const& eta): Expression(eta,true)
  {
    _xi = 4/(std::pow((2./eta-1),2)-1);
    double D = std::sqrt(1+4/_xi);
    double D2 = D*D;
    double Spp = 2/_xi +1;
    double Spm = 1;
    double I1 = std::log((Spp+D)/(Spp-D))/D;
    double Li2deltaOdiff = dilog(2*D/(D-Spp));
    double Li2deltaOsum = dilog(2*D/(D+Spp));
    double log_1 = std::log(D2*_xi);
    double log2_diff = std::pow(std::log(0.5*std::fabs(D-Spm)),2);
    double log2_sum =  std::pow(std::log(0.5*(D+Spm)),2);
    double Li2_diff = dilog(0.5*(D-Spm)/D);
    double Li2_sum = dilog(0.5*(D+Spm)/D);

    _S2 = 2 + Spp*(I1*D+Li2deltaOdiff-Li2deltaOsum)/D
          + log_1*(Spp*I1-2);
    _V2 = (0.5*D2 + Spp*(1+std::log(1/D)))*I1 + 2*std::log(1/_xi)
          - 4 + Spp*(log2_diff-log2_sum-2*Li2_diff + 2*Li2_sum)/D;
    _sud = 2*(Spp*std::log((Spp+D)/(Spp-D))-2*D)/D;
  }

  double Cm21qNC_ACOT::Local(double const& z) const
  {
    if(z>=1-eps4)
      {
        return 0;
      }
    const Integrator Integrand{[&] (double const& z)->double{
        return (g(z)-_sud)/(1-z);
      }};
    double extra_int = Integrand.integrate(z,1,eps5);
    return 2*apfel::CF*std::sqrt(1+4/_xi)*(_S2+_V2+_sud*std::log(1-z)+extra_int)/(1.+1./_xi);
  }

  double Cm21qNC_ACOT::g(double const& x) const
  {
    double epsi = 1/_xi;
    double z = x;
    double Spp = 1+2*epsi;
    double Spm = 1;
    double D = std::sqrt( 1 + 4 * epsi);
    double D2 = D*D;
    double D4 = D2*D2;
    double s1 = 0.5*(1-z)*((D-Spm)*z+D+Spm)/z;
    double s12 = s1*s1;
    double Dd = std::sqrt(s1*(s1+2) +4*epsi +1);
    double Dd2 = Dd*Dd;
    double Dd4 = Dd2*Dd2;
    double Lxi = std::log((Spp+s1-Dd)/(Spp+s1+Dd));
    double Ixi = (s1+2*epsi)/s12 + (s1+epsi)*Spp*Lxi/Dd/s12;

    double pref = s1*Dd2/(s1+epsi)/D/16;
    double term1 = -2*D4*Ixi;
    double term2 = 2*epsi*((s1+epsi)*(Dd2-6*epsi)*Lxi/Dd - 0.5*Dd2*(s1+Spp)/(s1+epsi) + 2*Dd2 -3*(s1+Spp));
    double term3 = -2*(D2-6*epsi)*(s1+epsi) - 4*epsi*s12 - 9*epsi*Spm*Spm + D2*(2*Spp-epsi) + 2*s1*(2*D2-4*epsi*Spm) + 0.5*(Dd2-6*(epsi+s1))*Spp*(s1+Spp)/(s1+epsi) - 2*D2*(D2+2*(2*epsi+s1)*Spm)/s1;
    double term4 = (s1+epsi)*(-2*D2*(D2+2*Spm*Spp)/s1 - 2*s1*(D2-6*epsi) - (Dd2-18*epsi)*Spp - 2*D2*(Spp+2*Spm))*Lxi/Dd;
    double fQ2 = 16*(term1+term2+term3+term4)/Dd4;

    return (1-x)*pref*fQ2;
  }

  double Cm21qNC_ACOT::Singular(double const& z) const
  {
    if(z>=1)
      {
        return 0;
      }
    return 2*apfel::CF*std::sqrt(1+4/_xi)*g(z)/(1-z)/(1.+1./_xi);
  }

  //_________________________________________________________________________________
  Cm21qNC_sub_ACOT::Cm21qNC_sub_ACOT(double const& eta): Expression(eta,true)
  {
    _xi = 4/(std::pow((2./eta-1),2)-1);
  }

  double Cm21qNC_sub_ACOT::Local(double const& z) const
  {
    if(z>=1)
      {
        return 0;
      }
    double lnxi = std::log(_xi);
    double ln1mz = std::log(1-z);
    double term1 = -ln1mz*(z*(z+2)+2*ln1mz-1) + 2*z;
    double term2 = lnxi*(2*ln1mz+0.5*z*(z+2));
    return 2 * apfel::CF * std::sqrt(1+4/_xi)* (term1+term2)/(1.+1./_xi);
  }

  double Cm21qNC_sub_ACOT::Singular(double const& z) const
  {
    if(z>=1)
      {
        return 0;
      }
    double lnxi = std::log(_xi);
    double ln1mz = std::log(1-z);
    double term1 = (1+z*z)*(lnxi-1-2*ln1mz)/(1-z);
    return 2 * apfel::CF * std::sqrt(1+4/_xi) * term1/(1.+1./_xi);
  }

  //_________________________________________________________________________________
  CmL1gNC_ACOT::CmL1gNC_ACOT(double const& eta):
    Expression(eta,true)
  {
  }
  double CmL1gNC_ACOT::Regular(double const& x) const
  {
    if (x >= 1)
      {
        return 0;
      }
    const double eta   = this->_eta;
    const double z     = eta * x;
    const double z2    = z * z;
    const double epsi  = ( 1 - eta ) / eta / 4;
    const double v     = sqrt( 1 - 4 * z / ( 1 - z ) * epsi );
    return 4 * TR * ( - 8 * epsi * z2 * log( ( 1 + v ) / ( 1 - v ) )
                      + 4 * v * z * ( 1 - z ) );
  }

  //_________________________________________________________________________________
  Cm20qNC_ACOT_chi::Cm20qNC_ACOT_chi(double const& eta): Expression(eta,true) {}

  double Cm20qNC_ACOT_chi::Local(double const& x) const
  {
    if(x>=1)
      {
        return 0;
      }
    return 1;
  }

  //_________________________________________________________________________________
  Cm21gNC_sub_ACOT_chi::Cm21gNC_sub_ACOT_chi(double const& eta): Expression(eta,true)
  {

  }

  double Cm21gNC_sub_ACOT_chi::Regular(double const& z) const
  {
    return 4 * TR * ( 1 - 2 * z + 2 * z * z );
  }

  //_________________________________________________________________________________
  Cm21qNC_ACOT_chi::Cm21qNC_ACOT_chi(double const& eta):
    Expression(eta,true)
  {
  }
  double Cm21qNC_ACOT_chi::Singular(double const& z) const
  {
    if(z>1)
      {
        return 0;
      }
    return 2 * apfel::CF *( (1+z*z)*(log((1-z)/z)-3./4.)/(1-z) + (9.+5*z)/4.);
  }
  double Cm21qNC_ACOT_chi::Local(double const& z) const
  {
    double ln1mz = log(1-z);
    double lnz = log(z);
    double z2 = z*z;
    double term1 = apfel::Pi2/3. + 7.*z/2. + z2;
    double term2 = ln1mz*(3. - z - z2/2. - ln1mz);
    double term3 = (z+z2/2.)*lnz-2*apfel::dilog(1-z);
    return -2 * apfel::CF * (term1 + term2 +term3 );
  }

  //_________________________________________________________________________________
  Cm31qNC_ACOT_chi::Cm31qNC_ACOT_chi(double const& eta):
    Expression(eta,true)
  {
  }
  double Cm31qNC_ACOT_chi::Regular(double const& x) const
  {
    return 2 * CF * ( - ( 1 + x ) * log( 1 - x ) - ( 1 + pow(x,2) ) * log(x) / ( 1 - x ) + 2 + x );
  }
  double Cm31qNC_ACOT_chi::Singular(double const& x) const
  {
    return 2 * CF * ( 2 * log( 1 - x ) - 3 / 2. ) / ( 1 - x );
  }
  double Cm31qNC_ACOT_chi::Local(double const& x) const
  {
    return 2 * CF * ( pow(log(1-x),2) - 3 * log( 1 - x ) / 2 - ( 2 * zeta2 + 9 / 2. ) );
  }

  //_________________________________________________________________________________
  CmL1qNC_ACOT_chi::CmL1qNC_ACOT_chi(double const& eta):
    Expression(eta,true)
  {
  }
  double CmL1qNC_ACOT_chi::Regular(double const& x) const
  {
    return 4 * CF * x;
  }

  //_________________________________________________________________________________
  Cm21gCC_sub_ACOT::Cm21gCC_sub_ACOT(double const& eta, double const& xi):
    Expression(eta,true),
    _xi(xi)
  {
  }

  double Cm21gCC_sub_ACOT::Regular(double const& z) const
  {
    return 4 * TR *  ( 1 - 2 * z + 2 * z * z );
  }

  //_________________________________________________________________________________
  Cm21gCC_general_mass::Cm21gCC_general_mass(double const& eta, double const& xi, double const& m1, double const& m2):
    Expression(eta, true),
    _xi(xi),
    _m1(m1 == 0.0 ? 0.01 : m1),
    _m2(m2 == 0.0 ? 0.01 : m2)
  {
  }

  double Cm21gCC_general_mass::Regular(double const& x) const
  {
    const double Q = sqrt(_xi)*(_m1+_m2);

    const double m1 = _m1/Q;
    const double m12 = m1*m1;
    const double m14 = m12*m12;
    const double m16 = m14*m12;
    const double m2 = _m2/Q;
    const double m22 = m2*m2;
    const double m24 = m22*m22;
    const double m26 = m24*m22;

    const double z = _eta*x;
    const double z2 = z*z;
    const double omz = 1-z;
    const double z1mz = z*omz;
    const double zo1mz = z/omz;

    const double nu =  sqrt(1-pow(m1+m2,2)*zo1mz);
    const double nub = sqrt(1-pow(m1-m2,2)*zo1mz);

    const double L =  log((1+(m12-m22)*zo1mz + nu*nub)/(1+(m12-m22)*zo1mz - nu*nub));
    const double Lt = log((1+(m22-m12)*zo1mz + nu*nub)/(1+(m22-m12)*zo1mz - nu*nub));

    const double term1 = -1 + z1mz*( 8 - 2*(m12+m22 - pow(m12-m22,2)));

    const double term2 = 1-2*z1mz + m12*(1+8*z-18*z2) + m22*(1-4*z+6*z2);
    const double term3 = -(m14+m24)*2*z*(1-3*z) + m12*m22*4*z*(1-5*z);
    const double term4 = (m16-m14*m22-m12*m24+m26)*2*z2;
    // term2,term3 and term4 with m1 <-> m2
    const double term6 = 1-2*z1mz + m22*(1+8*z-18*z2) + m12*(1-4*z+6*z2);
    const double term7 = term3; // term3 is symmetric
    const double term8 = (m26-m24*m12-m22*m14+m16)*2*z2;

    return 4*TR*(nu*nub*term1 + 0.5*L*(term2+term3+term4) + 0.5*Lt*(term6+term7+ term8));
  }

  //_________________________________________________________________________________
  CmL1gCC_general_mass::CmL1gCC_general_mass(double const& eta, double const& xi, double const& m1, double const& m2):
    Expression(eta, true),
    _xi(xi),
    _m1(m1 == 0.0 ? 0.01 : m1),
    _m2(m2 == 0.0 ? 0.01 : m2)
  {
  }

  double CmL1gCC_general_mass::Regular(double const& x) const
  {
    //TODO
    //can be calculated analytically!
    const double Q = sqrt(_xi)*(_m1+_m2);

    const double m1 = _m1/Q;
    const double m12 = m1*m1;
    const double m14 = m12*m12;
    const double m16 = m14*m12;
    const double m2 = _m2/Q;
    const double m22 = m2*m2;
    const double m24 = m22*m22;
    const double m26 = m24*m22;

    const double z = _eta*x;
    const double z2 = z*z;
    const double omz = 1-z;
    const double z1mz = z*omz;
    const double zo1mz = z/omz;

    const double nu =  sqrt(1-pow(m1+m2,2)*zo1mz);
    const double nub = sqrt(1-pow(m1-m2,2)*zo1mz);

    const double L =  log((1+(m12-m22)*zo1mz + nu*nub)/(1+(m12-m22)*zo1mz - nu*nub));
    const double Lt = log((1+(m22-m12)*zo1mz + nu*nub)/(1+(m22-m12)*zo1mz - nu*nub));

    const double termF2_1 = -1 + z1mz*( 8 - 2*(m12+m22 - pow(m12-m22,2)));

    const double termF2_2 = 1-2*z1mz + m12*(1+8*z-18*z2) + m22*(1-4*z+6*z2);
    const double termF2_3 = -(m14+m24)*2*z*(1-3*z) + m12*m22*4*z*(1-5*z);
    const double termF2_4 = (m16-m14*m22-m12*m24+m26)*2*z2;
    // termF2_2,termF2_3 and termF2_4 with m1 <-> m2
    const double termF2_6 = 1-2*z1mz + m22*(1+8*z-18*z2) + m12*(1-4*z+6*z2);
    const double termF2_7 = termF2_3; // termF2_3 is symmetric
    const double termF2_8 = (m26-m24*m12-m22*m14+m16)*2*z2;

    const double termF1_1 = -(1-4*z1mz);

    const double termF1_2 = (1-2*z1mz+(m12-m22)*2*z*(1-2*z));
    const double termF1_3 = pow(m12-m22,2)*2*z2;

    const double termF1_4 = (1-2*z1mz+(m22-m12)*2*z*(1-2*z));
    const double termF1_5 = termF1_3; //termF1_3 is symmetric

    const double result1 = nu*nub*(termF2_1-termF1_1);
    const double result2 = termF2_2+termF2_3+termF2_4 - termF1_2 - termF1_3;
    const double result3 = termF2_6+termF2_7+termF2_8 - termF1_4 - termF1_5;



    return 4*TR*(result1 + 0.5*L*result2 + 0.5*Lt*result3);
  }

  //_________________________________________________________________________________
  Cm31gCC_general_mass::Cm31gCC_general_mass(double const& eta, double const& xi, double const& m1, double const& m2):
    Expression(eta, true),
    _xi(xi),
    _m1(m1 == 0.0 ? 0.01 : m1),
    _m2(m2 == 0.0 ? 0.01 : m2)
  {
  }

  double Cm31gCC_general_mass::Regular(double const& x) const
  {
    const double Q = sqrt(_xi)*(_m1+_m2);

    const double m1 = _m1/Q;
    const double m12 = m1*m1;
    const double m14 = m12*m12;
    const double m2 = _m2/Q;
    const double m22 = m2*m2;
    const double m24 = m22*m22;

    const double z = _eta*x;
    const double z2 = z*z;
    const double omz = 1-z;
    const double z1mz = z*omz;
    const double zo1mz = z/omz;

    const double nu =  sqrt(1-pow(m1+m2,2)*zo1mz);
    const double nub = sqrt(1-pow(m1-m2,2)*zo1mz);

    const double L =  log((1+(m12-m22)*zo1mz + nu*nub)/(1+(m12-m22)*zo1mz - nu*nub));
    const double Lt = log((1+(m22-m12)*zo1mz + nu*nub)/(1+(m22-m12)*zo1mz - nu*nub));

    const double term1 = (m12-m22)*4*z1mz; // maybe put this to ZERO!

    const double term2 = 1-2*z1mz+(m12-m22)*2*z*(1-2*z) - (m14-m24)*2*z2;

    const double term3 = 1-2*z1mz+(m22-m12)*2*z*(1-2*z) - (m24-m14)*2*z2;

    return 2*TR*(nu*nub*term1 - L*term2 + Lt*term3);
  }

  ////////////////////////////////////////////
  /// F2-minus
  //_________________________________________________________________________________
  C22nsm_aSACOT_chi_0::C22nsm_aSACOT_chi_0(double const& eta, bool const& xdependent):
    Expression((xdependent ? eta : 1.),xdependent)
  {}

  double C22nsm_aSACOT_chi_0::Regular(double const& x) const
  {
    double const dl      = log(x);
    double const dl_2    = dl * dl;
    double const dl_3    = dl_2 * dl;
    double const dl1     = log(1-x);
    double const dl1_2   = dl1 * dl1;
    double const dl1_3   = dl1_2 * dl1;
    return
      - 84.18 - 1010 * x - 3.748 * dl_3 - 19.56 * dl_2 - 1.235 * dl
      - 17.19 * dl1_3 + 71.08 * dl1_2 - 663.0 * dl1 - 192.4 * dl * dl1_2 + 80.41 * dl_2 * dl1;
  }

  double C22nsm_aSACOT_chi_0::Singular(double const& x) const
  {
    double const dl1    = log(1-x);
    double const dl1_2  = dl1 * dl1;
    double const dl1_3  = dl1_2 * dl1;
    double const c2ns2b =
      + 14.2222 * dl1_3 - 61.3333 * dl1_2 - 31.105 * dl1 + 188.64;
    return c2ns2b / ( 1 - x );
  }

  double C22nsm_aSACOT_chi_0::Local(double const& x) const
  {
    double const dl1     = log(1-x);
    double const dl1_2   = dl1 * dl1;
    double const dl1_3   = dl1_2 * dl1;
    double const dl1_4   = dl1_3 * dl1;
    return + 3.55555 * dl1_4 - 20.4444 * dl1_3 - 15.5525 * dl1_2 + 188.64 * dl1 - 338.531 + 0.537;
  }

  //_________________________________________________________________________________
  C22nsm_aSACOT_chi_nf::C22nsm_aSACOT_chi_nf(double const& eta, bool const& xdependent):
    Expression((xdependent ? eta : 1.),xdependent)
  {}

  double C22nsm_aSACOT_chi_nf::Regular(double const& x) const
  {
    double const dl      = log(x);
    double const dl_2    = dl * dl;
    double const dl1     = log(1-x);
    double const dl1_2   = dl1 * dl1;
    return
      - 5.691 - 37.91 * x + 2.244 * dl_2 + 5.770 * dl - 1.707 * dl1_2  + 22.95 * dl1
      + 3.036 * dl_2 * dl1 + 17.97 * dl * dl1;
  }

  double C22nsm_aSACOT_chi_nf::Singular(double const& x) const
  {
    double const dl1    = log(1-x);
    double const dl1_2  = dl1 * dl1;
    double const c2ns2b =
      1.77778 * dl1_2 - 8.5926 * dl1 + 6.3489;
    return c2ns2b / ( 1 - x );
  }

  double C22nsm_aSACOT_chi_nf::Local(double const& x) const
  {
    double const dl1     = log(1-x);
    double const dl1_2   = dl1 * dl1;
    double const dl1_3   = dl1_2 * dl1;
    return
      0.592593 * dl1_3 - 4.2963 * dl1_2 + 6.3489 * dl1 + 46.844 - 0.0035;
  }

  ////////////////////////////////////////////
  /// F2-plus
  //_________________________________________________________________________________
  C22nsp_aSACOT_chi_0::C22nsp_aSACOT_chi_0(double const& eta, bool const& xdependent):
    Expression((xdependent ? eta : 1.),xdependent)
  {}

  double C22nsp_aSACOT_chi_0::Regular(double const& x) const
  {
    double const dl      = log(x);
    double const dl_2    = dl * dl;
    double const dl_3    = dl_2 * dl;
    double const dl1     = log(1-x);
    double const dl1_2   = dl1 * dl1;
    double const dl1_3   = dl1_2 * dl1;
    return
      - 69.59 - 1008 * x - 2.835 * dl_3 - 17.08 * dl_2 + 5.986 * dl
      - 17.19 * dl1_3 + 71.08 * dl1_2 - 660.7 * dl1 - 174.8 * dl * dl1_2 + 95.09 * dl_2 * dl1;
  }

  double C22nsp_aSACOT_chi_0::Singular(double const& x) const
  {
    double const dl1    = log(1-x);
    double const dl1_2  = dl1 * dl1;
    double const dl1_3  = dl1_2 * dl1;
    double const c2ns2b =
      + 14.2222 * dl1_3 - 61.3333 * dl1_2 - 31.105 * dl1 + 188.64;
    return c2ns2b / ( 1 - x );
  }

  double C22nsp_aSACOT_chi_0::Local(double const& x) const
  {
    double const dl1     = log(1-x);
    double const dl1_2   = dl1 * dl1;
    double const dl1_3   = dl1_2 * dl1;
    double const dl1_4   = dl1_3 * dl1;
    return + 3.55555 * dl1_4 - 20.4444 * dl1_3 - 15.5525 * dl1_2 + 188.64 * dl1 - 338.531 + 0.485;
  }

  //_________________________________________________________________________________
  C22nsp_aSACOT_chi_nf::C22nsp_aSACOT_chi_nf(double const& eta, bool const& xdependent):
    Expression((xdependent ? eta : 1.),xdependent)
  {}

  double C22nsp_aSACOT_chi_nf::Regular(double const& x) const
  {
    double const dl      = log(x);
    double const dl_2    = dl * dl;
    double const dl1     = log(1-x);
    double const dl1_2   = dl1 * dl1;
    return
      - 5.691 - 37.91 * x + 2.244 * dl_2 + 5.770 * dl - 1.707 * dl1_2
      + 22.95 * dl1 + 3.036 * dl_2 * dl1 + 17.97 * dl * dl1;
  }

  double C22nsp_aSACOT_chi_nf::Singular(double const& x) const
  {
    double const dl1    = log(1-x);
    double const dl1_2  = dl1 * dl1;
    double const c2ns2b =
      1.77778 * dl1_2 - 8.5926 * dl1 + 6.3489;
    return c2ns2b / ( 1 - x );
  }

  double C22nsp_aSACOT_chi_nf::Local(double const& x) const
  {
    double const dl1     = log(1-x);
    double const dl1_2   = dl1 * dl1;
    double const dl1_3   = dl1_2 * dl1;
    return 0.592593 * dl1_3 - 4.2963 * dl1_2 + 6.3489 * dl1 + 46.844 - 0.0035;
  }

  //_________________________________________________________________________________
  C22ps_aSACOT_chi::C22ps_aSACOT_chi(double const& eta, bool const& xdependent):
    Expression(eta,xdependent)
  {
  }
  double C22ps_aSACOT_chi::Regular(double const& x) const
  {
    double const dl     = log(x);
    double const dl_2   = dl * dl;
    double const dl_3   = dl_2 * dl;
    double const dl1    = log(1-x);
    double const dl1_2  = dl1 * dl1;
    double const dl1_3  = dl1_2 * dl1;
    return
      5.290 * ( 1 / x - 1 ) + 4.310 * dl_3 - 2.086 * dl_2 + 39.78 * dl - 0.101 * ( 1 - x ) * dl1_3 - ( 24.75 - 13.80 * x ) * dl_2 * dl1 + 30.23 * dl * dl1;
  }

  //_________________________________________________________________________________
  C22g_aSACOT_chi::C22g_aSACOT_chi(double const& eta, bool const& xdependent):
    Expression(eta,xdependent)
  {
  }
  double C22g_aSACOT_chi::Regular(double const& x) const
  {
    double const dl    = log(x);
    double const dl_2  = dl * dl;
    double const dl_3  = dl_2 * dl;
    double const dl1   = log(1-x);
    double const dl1_2 = dl1 * dl1;
    double const dl1_3 = dl1_2 * dl1;
    return
      1 / x * ( 11.90 + 1494.* dl1 ) + 5.319 * dl_3 - 59.48 * dl_2 - 284.8 * dl + 392.4 - 1483 * dl1
      + ( 6.445 + 209.4 * ( 1 - x ) ) * dl1_3 - 24.00 * dl1_2 - 724.1 * dl_2 * dl1 - 871.8 * dl * dl1_2;
  }
  double C22g_aSACOT_chi::Local(double const&) const
  {
    return - 0.28;
  }

  ////////////////////////////////////////////
  /// FL-minus
  //_________________________________________________________________________________
  CL2nsm_aSACOT_chi_0::CL2nsm_aSACOT_chi_0(double const& eta, bool const& xdependent):
    Expression((xdependent ? eta : 1.),xdependent)
  {}

  double CL2nsm_aSACOT_chi_0::Regular(double const& x) const
  {
    double const dl      = log(x);
    double const dl_2    = dl * dl;
    double const dl1     = log(1-x);
    double const dl1_2   = dl1 * dl1;
    return
      - 52.27 + 100.8 * x + ( 23.29 * x - 0.043 ) * dl_2 - 22.21 * dl + 13.30 * dl1_2
      - 59.12 * dl1 - 141.7 * dl * dl1;
  }

  double CL2nsm_aSACOT_chi_0::Local(double const&) const
  {
    return - 0.150;
  }

  //_________________________________________________________________________________
  CL2nsm_aSACOT_chi_nf::CL2nsm_aSACOT_chi_nf(double const& eta, bool const& xdependent):
    Expression((xdependent ? eta : 1.),xdependent)
  {}

  double CL2nsm_aSACOT_chi_nf::Regular(double const& x) const
  {
    double const dl      = log(x);
    double const dl1     = log(1-x);
    return
      16 / 27. * ( 6 * x * dl1 - 12 * x * dl - 25 * x + 6 );
  }

  ////////////////////////////////////////////
  /// FL-plus
  //_________________________________________________________________________________
  CL2nsp_aSACOT_chi_0::CL2nsp_aSACOT_chi_0(double const& eta, bool const& xdependent):
    Expression((xdependent ? eta : 1.),xdependent)
  {}

  double CL2nsp_aSACOT_chi_0::Regular(double const& x) const
  {
    double const dl      = log(x);
    double const dl_2    = dl * dl;
    double const dl1     = log(1-x);
    double const dl1_2   = dl1 * dl1;
    return
      - 40.41 + 97.48 * x + ( 26.56 * x - 0.031 ) * dl_2 - 14.85 * dl + 13.62 * dl1_2
      - 55.79 * dl1 - 150.5 * dl * dl1;
  }

  double CL2nsp_aSACOT_chi_0::Local(double const&) const
  {
    return - 0.164;
  }

  //_________________________________________________________________________________
  CL2nsp_aSACOT_chi_nf::CL2nsp_aSACOT_chi_nf(double const& eta, bool const& xdependent):
    Expression((xdependent ? eta : 1.),xdependent)
  {}

  double CL2nsp_aSACOT_chi_nf::Regular(double const& x) const
  {
    double const dl      = log(x);
    double const dl1     = log(1-x);
    return
      16 / 27. * ( 6 * x * dl1 - 12 * x * dl - 25 * x + 6 );
  }

  //_________________________________________________________________________________
  CL2g_aSACOT_chi::CL2g_aSACOT_chi(double const& eta, bool const& xdependent):
    Expression((xdependent ? eta : 1.),xdependent)
  {}

  double CL2g_aSACOT_chi::Regular(double const& x) const
  {
    double const dl    = log(x);
    double const dl_2  = dl * dl;
    double const dl1   = log(1-x);
    double const dl1_2 = dl1 * dl1;
    double const omx   = 1 - x;
    return
      ( 94.74 - 49.20 * x ) * omx * dl1_2 + 864.8 * omx * dl1 + 1161 * x * dl * dl1
      + 60.06 * x * dl_2 + 39.66 * omx * dl - 5.333 * ( 1 / x - 1 );
  }

  //_________________________________________________________________________________
  CL2ps_aSACOT_chi::CL2ps_aSACOT_chi(double const& eta, bool const& xdependent):
    Expression((xdependent ? eta : 1.),xdependent)
  {}

  double CL2ps_aSACOT_chi::Regular(double const& x) const
  {
    double const dl     = log(x);
    double const dl_2   = dl * dl;
    double const dl1    = log(1-x);
    double const omx    = 1 - x;
    double const omx2   = omx * omx;
    double const omx3   = omx2 * omx;
    return
      ( 15.94 - 5.212 * x ) * omx2 * dl1 + ( 0.421 + 1.520 * x ) * dl_2 + 28.09 * omx * dl
      - ( 2.370 / x - 19.27 ) * omx3;
  }

  ////////////////////////////////////////////
  /// F3-plus

  // note here the change in notation from the files from A. Vogt!
  // minus will be F(nu)-F(anti-nu)
  // plus  will be F(nu)+F(nu)
  // there are the other way around in A. Vogts files
  //_________________________________________________________________________________
  C32nsm_aSACOT_chi_0::C32nsm_aSACOT_chi_0(double const& eta, bool const& xdependent):
    Expression((xdependent ? eta : 1.),xdependent)
  {}

  double C32nsm_aSACOT_chi_0::Regular(double const& x) const
  {
    double const dl    = log(x);
    double const dl_2  = dl*dl;
    double const dl_3  = dl_2*dl;
    double const dl1   = log(1.-x);
    double const dl1_2 = dl1*dl1;
    double const dl1_3 = dl1_2*dl1;
    return
      - 206.1 - 576.8 * x - 3.922 * dl_3 - 33.31 * dl_2 - 67.60 * dl
      - 15.20 * dl1_3 + 94.61 * dl1_2 - 409.6 * dl1
      - 147.9 * dl * dl1_2;
  }

  double C32nsm_aSACOT_chi_0::Singular(double const& x) const
  {
    double const dl1    = log(1-x);
    double const dl1_2  = dl1 * dl1;
    double const dl1_3  = dl1_2 * dl1;
    double const c2ns2b =
      + 14.2222 * dl1_3 - 61.3333 * dl1_2 - 31.105 * dl1 + 188.64;
    return c2ns2b / ( 1 - x );
  }

  double C32nsm_aSACOT_chi_0::Local(double const& x) const
  {
    double const dl1   = log(1-x);
    double const dl1_2 = dl1 * dl1;
    double const dl1_3 = dl1_2 * dl1;
    double const dl1_4 = dl1_3 * dl1;
    return
      + 3.55555 * dl1_4 - 20.4444 * dl1_3 - 15.5525 * dl1_2
      + 188.64 * dl1 - 338.531 - 0.104;
  }

  //_________________________________________________________________________________
  C32nsm_aSACOT_chi_nf::C32nsm_aSACOT_chi_nf(double const& eta, bool const& xdependent):
    Expression((xdependent ? eta : 1.),xdependent)
  {}

  double C32nsm_aSACOT_chi_nf::Regular(double const& x) const
  {
    double const dl    = log(x);
    double const dl_2  = dl*dl;
    double const dl1   = log(1.-x);
    double const dl1_2 = dl1*dl1;
    double const dl1_3 = dl1_2*dl1;
    return
      - 6.337 - 14.97 * x + 2.207 * dl_2 + 8.683 * dl + 0.042 * dl1_3 - 0.808 * dl1_2 + 25.00 * dl1
      + 9.684 * dl * dl1 ;
  }

  double C32nsm_aSACOT_chi_nf::Singular(double const& x) const
  {
    double const dl1    = log(1-x);
    double const dl1_2  = dl1 * dl1;
    double const c2ns2b =
      + 1.77778 * dl1_2 - 8.5926 * dl1 + 6.3489;
    return c2ns2b / ( 1 - x );
  }

  double C32nsm_aSACOT_chi_nf::Local(double const& x) const
  {
    double const dl1   = log(1-x);
    double const dl1_2 = dl1 * dl1;
    double const dl1_3 = dl1_2 * dl1;
    return
      0.592593 * dl1_3 - 4.2963 * dl1_2 + 6.3489 * dl1 + 46.844 + 0.013;
  }

  //again note the change in notation!
  //_________________________________________________________________________________
  C32nsp_aSACOT_chi_0::C32nsp_aSACOT_chi_0(double const& eta, bool const& xdependent):
    Expression((xdependent ? eta : 1.),xdependent)
  {}

  double C32nsp_aSACOT_chi_0::Regular(double const& x) const
  {
    double const dl      = log(x);
    double const dl_2    = dl * dl;
    double const dl_3    = dl_2 * dl;
    double const dl1     = log(1-x);
    double const dl1_2   = dl1 * dl1;
    double const dl1_3   = dl1_2 * dl1;
    return
      - 242.9 - 467.2 * x - 3.049 * dl_3 - 30.14 * dl_2 - 79.14 * dl - 15.20 * dl1_3
      + 94.61 * dl1_2 - 396.1 * dl1 - 92.43 * dl * dl1_2;

  }

  double C32nsp_aSACOT_chi_0::Singular(double const& x) const
  {
    double const dl1    = log(1-x);
    double const dl1_2  = dl1 * dl1;
    double const dl1_3  = dl1_2 * dl1;
    double const c2ns2b =
      + 14.2222 * dl1_3 - 61.3333 * dl1_2 - 31.105 * dl1 + 188.64;
    return c2ns2b / ( 1 - x );
  }

  double C32nsp_aSACOT_chi_0::Local(double const& x) const
  {
    double const dl1   = log(1-x);
    double const dl1_2 = dl1 * dl1;
    double const dl1_3 = dl1_2 * dl1;
    double const dl1_4 = dl1_3 * dl1;
    return
      + 3.55555 * dl1_4 - 20.4444 * dl1_3 - 15.5525 * dl1_2 + 188.64 * dl1 - 338.531 - 0.152;
  }

  //_________________________________________________________________________________
  C32nsp_aSACOT_chi_nf::C32nsp_aSACOT_chi_nf(double const& eta, bool const& xdependent):
    Expression((xdependent ? eta : 1.),xdependent)
  {}

  double C32nsp_aSACOT_chi_nf::Regular(double const& x) const
  {
    double const dl      = log(x);
    double const dl_2    = dl * dl;
    double const dl1     = log(1-x);
    double const dl1_2   = dl1 * dl1;
    double const dl1_3   = dl1_2 * dl1;
    return
      - 6.337 - 14.97 * x + 2.207 * dl_2 + 8.683 * dl + 0.042 * dl1_3 - 0.808 * dl1_2
      + 25.00 * dl1 + 9.684 * dl * dl1;
  }

  double C32nsp_aSACOT_chi_nf::Singular(double const& x) const
  {
    double const dl1    = log(1-x);
    double const dl1_2  = dl1 * dl1;
    double const c2ns2b =
      1.77778 * dl1_2 - 8.5926 * dl1 + 6.3489;
    return c2ns2b / ( 1 - x );
  }

  double C32nsp_aSACOT_chi_nf::Local(double const& x) const
  {
    double const dl1   = log(1-x);
    double const dl1_2 = dl1 * dl1;
    double const dl1_3 = dl1_2 * dl1;
    return
      0.592593 * dl1_3 - 4.2963 * dl1_2 + 6.3489 * dl1 + 46.844 + 0.013;
  }
}