//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <cmath>

using namespace std;

namespace apfel
{
  /**
   * @brief Real dilogarithm \f$\mathrm{Li}_2(z)\f$
   * @param x real argument
   * @note Implementation translated by R.Brun from CERNLIB DILOG function C332
   * @return \f$\mathrm{Li}_2(z)\f$
   */
  double dilog(double const& x)
  {
     const double PI  = M_PI;
     const double HF  = 0.5;
     const double PI2 = PI * PI;
     const double PI3 = PI2 / 3;
     const double PI6 = PI2 / 6;
     const double PI12 = PI2 / 12;
     const double C[20] = { 0.42996693560813697,
			    0.40975987533077105,
			   -0.01858843665014592,
			    0.00145751084062268,
			   -0.00014304184442340,
			    0.00001588415541880,
			   -0.00000190784959387,
			    0.00000024195180854,
			   -0.00000003193341274,
			    0.00000000434545063,
		           -0.00000000060578480,
			    0.00000000008612098,
			   -0.00000000001244332,
			    0.00000000000182256,
			   -0.00000000000027007,
			    0.00000000000004042,
			   -0.00000000000000610,
			    0.00000000000000093,
			   -0.00000000000000014,
			    0.00000000000000002};

     double T, H, Y, S, A, ALFA, B1, B2, B0;

     if (x == 1) {
       H = PI6;
     } else if (x == -1) {
       H = -PI12;
     } else {
       T = -x;
       if (T <= -2) {
	 Y = -1/(1+T);
	 S = 1;
	 B1= log(-T);
	 B2= log(1+1/T);
	 A = -PI3+HF*(B1*B1-B2*B2);
       } else if (T < -1) {
	 Y = -1-T;
	 S = -1;
	 A = log(-T);
	 A = -PI6+A*(A+log(1+1/T));
       } else if (T <= -0.5) {
	 Y = -(1+T)/T;
	 S = 1;
	 A = log(-T);
	 A = -PI6+A*(-HF*A+log(1+T));
       } else if (T < 0) {
	 Y = -T/(1+T);
	 S = -1;
	 B1= log(1+T);
	 A = HF*B1*B1;
       } else if (T <= 1) {
	 Y = T;
	 S = 1;
	 A = 0;
       } else {
	 Y = 1/T;
	 S = -1;
	 B1= log(T);
	 A = PI6+HF*B1*B1;
       }
       H    = Y+Y-1;
       ALFA = H+H;
       B1   = 0;
       B2   = 0;
       for (int i=19;i>=0;i--){
	 B0 = C[i] + ALFA*B1-B2;
	 B2 = B1;
	 B1 = B0;
       }
       H = -(S*(B0-H*B2)+A);
     }
     return H;
  };

}

