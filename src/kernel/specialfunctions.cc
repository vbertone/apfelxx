//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/specialfunctions.h"
#include "apfel/messages.h"
#include "apfel/hpolyweights.h"

#include <cmath>
#include <cstdio>

namespace apfel
{
  //_________________________________________________________________________________
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
                           0.00000000000000002
                         };

    double T, H, Y, S, A, ALFA, B1, B2, B0;

    if (x == 1)
      {
        H = PI6;
      }
    else if (x == -1)
      {
        H = -PI12;
      }
    else
      {
        T = -x;
        if (T <= -2)
          {
            Y = -1/(1+T);
            S = 1;
            B1= log(-T);
            B2= log(1+1/T);
            A = -PI3+HF*(B1*B1-B2*B2);
          }
        else if (T < -1)
          {
            Y = -1-T;
            S = -1;
            A = log(-T);
            A = -PI6+A*(A+log(1+1/T));
          }
        else if (T <= -0.5)
          {
            Y = -(1+T)/T;
            S = 1;
            A = log(-T);
            A = -PI6+A*(-HF*A+log(1+T));
          }
        else if (T < 0)
          {
            Y = -T/(1+T);
            S = -1;
            B1= log(1+T);
            A = HF*B1*B1;
          }
        else if (T <= 1)
          {
            Y = T;
            S = 1;
            A = 0;
          }
        else
          {
            Y = 1/T;
            S = -1;
            B1= log(T);
            A = PI6+HF*B1*B1;
          }
        H    = Y+Y-1;
        ALFA = H+H;
        B1   = 0;
        B2   = 0;
        for (int i=19; i>=0; i--)
          {
            B0 = C[i] + ALFA*B1-B2;
            B2 = B1;
            B1 = B0;
          }
        H = -(S*(B0-H*B2)+A);
      }
    return H;
  }

  //_________________________________________________________________________________
  double wgplg(int const& n, int const& p, double const& x)
  {
    int p1;
    int i,l,k,m,n1;
    double u[4],s1[5][5],c[5][5];
    double a[31][11];
    double x1,h,alfa,r,q,c1,c2,b0=0,b1,b2;
    double result;
    double v[5],sk,sm;

    double  fct[]= {1.0,1.0,2.0,6.0,24.0};
    double  sgn[]= {1.0,-1.0,1.0,-1.0,1.0 };

    // DATA INDEX /1,2,3,4,6*0,5,6,7,7*0,8,9,8*0,10/
    int index[]= {0,1,2,3,4,
                  0,0,0,0,0,0,
                  5,6,7,
                  0,0,0,0,0,0,0,
                  8,9,
                  0,0,0,0,0,0,0,0,
                  10
                 };

    int nc[]= {0,24,26,28,30,22,24,26,19,22,17};

    c1=1.3333333333333;
    c2=0.3333333333333;

    s1[1][1]=1.6449340668482;
    s1[1][2]=1.2020569031596;
    s1[1][3]=1.0823232337111;
    s1[1][4]=1.0369277551434;
    s1[2][1]=1.2020569031596;
    s1[2][2]=2.7058080842778e-1;
    s1[2][3]=9.6551159989444e-2;
    s1[3][1]=1.0823232337111;
    s1[3][2]=9.6551159989444e-2;
    s1[4][1]=1.0369277551434;

    c[1][1]=1.6449340668482;
    c[1][2]=1.2020569031596;
    c[1][3]=1.0823232337111;
    c[1][4]=1.0369277551434;
    c[2][1]=0.0000000000000;
    c[2][2]=-1.8940656589945;
    c[2][3]=-3.0142321054407;
    c[3][1]=1.8940656589945;
    c[3][2]=3.0142321054407;
    c[4][1]=0.0000000000000;

    a[0][1]=.96753215043498;
    a[1][1]=.16607303292785;
    a[2][1]=.02487932292423;
    a[3][1]=.00468636195945;
    a[4][1]=.00100162749616;
    a[5][1]=.00023200219609;
    a[6][1]=.00005681782272;
    a[7][1]=.00001449630056;
    a[8][1]=.00000381632946;
    a[9][1]=.00000102990426;
    a[10][1]=.00000028357538;
    a[11][1]=.00000007938705;
    a[12][1]=.00000002253670;
    a[13][1]=.00000000647434;
    a[14][1]=.00000000187912;
    a[15][1]=.00000000055029;
    a[16][1]=.00000000016242;
    a[17][1]=.00000000004827;
    a[18][1]=.00000000001444;
    a[19][1]=.00000000000434;
    a[20][1]=.00000000000131;
    a[21][1]=.00000000000040;
    a[22][1]=.00000000000012;
    a[23][1]=.00000000000004;
    a[24][1]=.00000000000001;

    a[0][2]=.95180889127832;
    a[1][2]=.43131131846532;
    a[2][2]=.10002250714905;
    a[3][2]=.02442415595220;
    a[4][2]=.00622512463724;
    a[5][2]=.00164078831235;
    a[6][2]=.00044407920265;
    a[7][2]=.00012277494168;
    a[8][2]=.00003453981284;
    a[9][2]=.00000985869565;
    a[10][2]=.00000284856995;
    a[11][2]=.00000083170847;
    a[12][2]=.00000024503950;
    a[13][2]=.00000007276496;
    a[14][2]=.00000002175802;
    a[15][2]=.00000000654616;
    a[16][2]=.00000000198033;
    a[17][2]=.00000000060204;
    a[18][2]=.00000000018385;
    a[19][2]=.00000000005637;
    a[20][2]=.00000000001735;
    a[21][2]=.00000000000536;
    a[22][2]=.00000000000166;
    a[23][2]=.00000000000052;
    a[24][2]=.00000000000016;
    a[25][2]=.00000000000005;
    a[26][2]=.00000000000002;

    a[0][3]=.98161027991365;
    a[1][3]=.72926806320726;
    a[2][3]=.22774714909321;
    a[3][3]=.06809083296197;
    a[4][3]=.02013701183064;
    a[5][3]=.00595478480197;
    a[6][3]=.00176769013959;
    a[7][3]=.00052748218502;
    a[8][3]=.00015827461460;
    a[9][3]=.00004774922076;
    a[10][3]=.00001447920408;
    a[11][3]=.00000441154886;
    a[12][3]=.00000135003870;
    a[13][3]=.00000041481779;
    a[14][3]=.00000012793307;
    a[15][3]=.00000003959070;
    a[16][3]=.00000001229055;
    a[17][3]=.00000000382658;
    a[18][3]=.00000000119459;
    a[19][3]=.00000000037386;
    a[20][3]=.00000000011727;
    a[21][3]=.00000000003687;
    a[22][3]=.00000000001161;
    a[23][3]=.00000000000366;
    a[24][3]=.00000000000116;
    a[25][3]=.00000000000037;
    a[26][3]=.00000000000012;
    a[27][3]=.00000000000004;
    a[28][3]=.00000000000001;

    a[0][4]=1.0640521184614;
    a[1][4]=1.0691720744981;
    a[2][4]=.41527193251768;
    a[3][4]=.14610332936222;
    a[4][4]=.04904732648784;
    a[5][4]=.01606340860396;
    a[6][4]=.00518889350790;
    a[7][4]=.00166298717324;
    a[8][4]=.00053058279969;
    a[9][4]=.00016887029251;
    a[10][4]=.00005368328059;
    a[11][4]=.00001705923313;
    a[12][4]=.00000542174374;
    a[13][4]=.00000172394082;
    a[14][4]=.00000054853275;
    a[15][4]=.00000017467795;
    a[16][4]=.00000005567550;
    a[17][4]=.00000001776234;
    a[18][4]=.00000000567224;
    a[19][4]=.00000000181313;
    a[20][4]=.00000000058012;
    a[21][4]=.00000000018579;
    a[22][4]=.00000000005955;
    a[23][4]=.00000000001911;
    a[24][4]=.00000000000614;
    a[25][4]=.00000000000197;
    a[26][4]=.00000000000063;
    a[27][4]=.00000000000020;
    a[28][4]=.00000000000007;
    a[29][4]=.00000000000002;
    a[30][4]=.00000000000001;

    a[0][5]=.97920860669175;
    a[1][5]=.08518813148683;
    a[2][5]=.00855985222013;
    a[3][5]=.00121177214413;
    a[4][5]=.00020722768531;
    a[5][5]=.00003996958691;
    a[6][5]=.00000838064065;
    a[7][5]=.00000186848945;
    a[8][5]=.00000043666087;
    a[9][5]=.00000010591733;
    a[10][5]=.00000002647892;
    a[11][5]=.00000000678700;
    a[12][5]=.00000000177654;
    a[13][5]=.00000000047342;
    a[14][5]=.00000000012812;
    a[15][5]=.00000000003514;
    a[16][5]=.00000000000975;
    a[17][5]=.00000000000274;
    a[18][5]=.00000000000077;
    a[19][5]=.00000000000022;
    a[20][5]=.00000000000006;
    a[21][5]=.00000000000002;
    a[22][5]=.00000000000001;

    a[0][6]=.95021851963952;
    a[1][6]=.29052529161433;
    a[2][6]=.05081774061716;
    a[3][6]=.00995543767280;
    a[4][6]=.00211733895031;
    a[5][6]=.00047859470550;
    a[6][6]=.00011334321308;
    a[7][6]=.00002784733104;
    a[8][6]=.00000704788108;
    a[9][6]=.00000182788740;
    a[10][6]=.00000048387492;
    a[11][6]=.00000013033842;
    a[12][6]=.00000003563769;
    a[13][6]=.00000000987174;
    a[14][6]=.00000000276586;
    a[15][6]=.00000000078279;
    a[16][6]=.00000000022354;
    a[17][6]=.00000000006435;
    a[18][6]=.00000000001866;
    a[19][6]=.00000000000545;
    a[20][6]=.00000000000160;
    a[21][6]=.00000000000047;
    a[22][6]=.00000000000014;
    a[23][6]=.00000000000004;
    a[24][6]=.00000000000001;

    a[0][7]=.95064032186777;
    a[1][7]=.54138285465171;
    a[2][7]=.13649979590321;
    a[3][7]=.03417942328207;
    a[4][7]=.00869027883583;
    a[5][7]=.00225284084155;
    a[6][7]=.00059516089806;
    a[7][7]=.00015995617766;
    a[8][7]=.00004365213096;
    a[9][7]=.00001207474688;
    a[10][7]=.00000338018176;
    a[11][7]=.00000095632476;
    a[12][7]=.00000027313129;
    a[13][7]=.00000007866968;
    a[14][7]=.00000002283195;
    a[15][7]=.00000000667205;
    a[16][7]=.00000000196191;
    a[17][7]=.00000000058018;
    a[18][7]=.00000000017246;
    a[19][7]=.00000000005151;
    a[20][7]=.00000000001545;
    a[21][7]=.00000000000465;
    a[22][7]=.00000000000141;
    a[23][7]=.00000000000043;
    a[24][7]=.00000000000013;
    a[25][7]=.00000000000004;
    a[26][7]=.00000000000001;

    a[0][8]=.98800011672229;
    a[1][8]=.04364067609601;
    a[2][8]=.00295091178278;
    a[3][8]=.00031477809720;
    a[4][8]=.00004314846029;
    a[5][8]=.00000693818230;
    a[6][8]=.00000124640350;
    a[7][8]=.00000024293628;
    a[8][8]=.00000005040827;
    a[9][8]=.00000001099075;
    a[10][8]=.00000000249467;
    a[11][8]=.00000000058540;
    a[12][8]=.00000000014127;
    a[13][8]=.00000000003492;
    a[14][8]=.00000000000881;
    a[15][8]=.00000000000226;
    a[16][8]=.00000000000059;
    a[17][8]=.00000000000016;
    a[18][8]=.00000000000004;
    a[19][8]=.00000000000001;

    a[0][9]=.95768506546350;
    a[1][9]=.19725249679534;
    a[2][9]=.02603370313918;
    a[3][9]=.00409382168261;
    a[4][9]=.00072681707110;
    a[5][9]=.00014091879261;
    a[6][9]=.00002920458914;
    a[7][9]=.00000637631144;
    a[8][9]=.00000145167850;
    a[9][9]=.00000034205281;
    a[10][9]=.00000008294302;
    a[11][9]=.00000002060784;
    a[12][9]=.00000000522823;
    a[13][9]=.00000000135066;
    a[14][9]=.00000000035451;
    a[15][9]=.00000000009436;
    a[16][9]=.00000000002543;
    a[17][9]=.00000000000693;
    a[18][9]=.00000000000191;
    a[19][9]=.00000000000053;
    a[20][9]=.00000000000015;
    a[21][9]=.00000000000004;
    a[22][9]=.00000000000001;

    a[0][10]=.99343651671347;
    a[1][10]=.02225770126826;
    a[2][10]=.00101475574703;
    a[3][10]=.00008175156250;
    a[4][10]=.00000899973547;
    a[5][10]=.00000120823987;
    a[6][10]=.00000018616913;
    a[7][10]=.00000003174723;
    a[8][10]=.00000000585215;
    a[9][10]=.00000000114739;
    a[10][10]=.00000000023652;
    a[11][10]=.00000000005082;
    a[12][10]=.00000000001131;
    a[13][10]=.00000000000259;
    a[14][10]=.00000000000061;
    a[15][10]=.00000000000015;
    a[16][10]=.00000000000004;
    a[17][10]=.00000000000001;

    if((n<1) || (n>4) || (p<1) || (p>4) || (n+p>5))
      {
        result=0.0;
        printf("****CERN SUBROUTINE RESULT...ILLEGAL VALUES n=%d p=%d",n,p);
        return(result);
      }

    if(x==sgn[0])
      {
        result=s1[n][p];
        return(result);
      }

    if((x>fct[2]) || (x<sgn[1]))
      {
        x1=sgn[0]/x;
        h=c1*x1+c2;
        alfa=h+h;
        v[0]=sgn[0];
        v[1]=log(-x);
        for(l=2; l<=n+p; l++)
          {
            v[l]=v[1]*v[l-1]/l;
          }
        sk=0.0;
        for(k=0; k<=p-1; k++)
          {
            p1=p-k;
            r=pow(x1,p1)/(fct[p1]*fct[n-1]);
            sm=0.0;
            for(m=0; m<=k; m++)
              {
                n1=n+k-m;
                l=index[10*n1+p1-10];
                b1=0.0;
                b2=0.0;
                for(i=nc[l]; i>=0; i--)
                  {
                    b0=a[i][l]+alfa*b1-b2;
                    b2=b1;
                    b1=b0;
                  }
                q=(fct[n1-1]/fct[k-m])*(b0-h*b2)*r/pow(p1,n1);
                sm=sm+v[m]*q;
              }
            sk=sk+sgn[k]*sm;
          }
        sm=0.0;
        for(m=0; m<=n-1; m++)
          {
            sm=sm+v[m]*c[n-m][p];
          }

        result=sgn[n]*sk+sgn[p]*(sm+v[n+p]);
        return(result);
      }

    if(x>0.5)
      {
        x1=sgn[0]-x;
        h=c1*x1+c2;
        alfa=h+h;
        v[0]=sgn[0];
        u[0]=sgn[0];
        v[1]=log(x1);
        u[1]=log(x);
        for(l=2; l<=p; l++)
          {
            v[l]=v[1]*v[l-1]/l;
          }
        for(l=2; l<=n; l++)
          {
            u[l]=u[1]*u[l-1]/l;
          }
        sk=0.0;
        for(k=0; k<=n-1; k++)
          {
            p1=n-k;
            r=pow(x1,p1)/fct[p1];
            sm=0.0;
            for(m=0; m<=p-1; m++)
              {
                n1=p-m;
                l=index[10*n1+p1-10];
                b1=0.0;
                b2=0.0;
                for(i=nc[l]; i>=0; i--)
                  {
                    b0=a[i][l]+alfa*b1-b2;
                    b2=b1;
                    b1=b0;
                  }
                q=sgn[m]*(b0-h*b2)*r/pow(p1,n1);
                sm=sm+v[m]*q;
              }
            sk=sk+u[k]*(s1[p1][p]-sm);
          }
        result=sk+sgn[p]*u[n]*v[p];
        return(result);
      }

    l=index[10*n+p-10];
    h=c1*x+c2;
    alfa=h+h;
    b1=0.0;
    b2=0.0;
    for(i=nc[l]; i>=0; i--)
      {
        b0=a[i][l]+alfa*b1-b2;
        b2=b1;
        b1=b0;
      }
    result=(b0-h*b2)*pow(x,p)/(fct[p]*pow(p,n));
    return(result);
  }

  //_________________________________________________________________________________
  std::pair<int, int> WeightAndIndex(std::vector<int> const& w)
  {
    // Unpack vector of weights
    std::vector<int> uw;
    for (auto const n : w)
      {
        const int d = std::max(std::abs(n), 1);
        uw.resize(uw.size() + d);
        uw.back() = n / d;
      }

    // Now return the weight size and the corresponding index in the
    // basis.
    return {uw.size(), WeightIndex.at(uw)};
  }

  //_________________________________________________________________________________
  double hpoly(std::vector<int> const& w, double const& x)
  {
    // Make sure that the argument is inside the validity range.
    if (x <= 0 || x > sqrt(2) - 1)
      throw std::runtime_error(error("hpoly", "Argument out of range."));

    // Get weight size and index in the basis
    const std::pair<int, int> wi = WeightAndIndex(w);

    // Make sure that the weight does not exceed 5
    if (wi.first < 1 || wi.first > 5)
      throw std::runtime_error(error("hpoly", "Weight out of range."));

    // If the weight is one use the explicit expressions
    if (wi.first == 1)
      if (w[0] == - 1)
        return log(1 + x);
      else if (w[0] == 0)
        return log(x);
      else
        return - log(1 - x);
    else
      {
        // Find position on the coefficient vectors
        const int step = bv[wi.first-1] * bs;
        std::vector<double>::const_iterator in;
        if (wi.first == 2)
          in = da2.begin() + bs * wi.second;
        else if (wi.first == 3)
          in = da3.begin() + bs * wi.second;
        else if (wi.first == 4)
          in = da4.begin() + bs * wi.second;
        else
          in = da5.begin() + bs * wi.second;

        // Compute HPL
        const double u   =   log(1 + x);
        const double v   = - log(1 - x);
        const double dlu =   log(u);
        const double dlv =   log(v);
        double uk = u;
        double vk = v;
        std::vector<double> tu(wi.first, 0);
        std::vector<double> tv(wi.first - 1, 0);
        for (std::vector<double>::const_iterator it = in; it < in + bs; it++)
          {
            int i = 0;
            for (int k = 0; k < wi.first; k++)
              tu[k] += uk * *(it + i++ * step);
            for (int k = 0; k < wi.first - 1; k++)
              tv[k] += vk * *(it + i++ * step);
            uk *= u;
            vk *= v;
          }
        double hpl = tu[0];
        for (int k = 1; k < wi.first; k++)
          hpl += tu[k] * pow(dlu, k) + tv[k-1] * pow(dlv, k - 1);

        return hpl;
      }
  }

  //_________________________________________________________________________________
  std::map<int, std::vector<double>> hpoly(double const& x, int const& wmax)
  {
    // Make sure that the argument is inside the validity range.
    if (x <= 0 || x > sqrt(2) - 1)
      throw std::runtime_error(error("hpoly", "Argument out of range."));

    // Make sure that wmax is within the validity range.
    if (wmax < 1 || wmax > 5)
      throw std::runtime_error(error("hpoly", "Max weight out of range."));

    // Relevant variables
    const double u =   log(1 + x);
    const double v = - log(1 - x);
    std::vector<double> dlu{1, log(u)};
    std::vector<double> dlv{1, log(v)};

    // Starting points of the iterator
    std::vector<std::vector<double>::const_iterator> vin{da2.begin(), da3.begin(), da4.begin(), da5.begin()};

    // Initialise output
    std::map<int, std::vector<double>> hpls;

    // Weight 1
    hpls.insert({1, {u, log(x), v}});

    // Loop over remaining weights
    for (int iw = 2; iw <= wmax; iw++)
      {
        const int step = bv[iw-1] * bs;
        std::vector<double> hpb(bv[iw-1], 0);
        for (int ib = 0; ib < bv[iw-1]; ib++)
          {
            // Compute HPL
            double uk = u;
            double vk = v;
            std::vector<double> tu(iw, 0);
            std::vector<double> tv(iw - 1, 0);
            for (std::vector<double>::const_iterator it = vin[iw-2]; it < vin[iw-2] + bs; it++)
              {
                int i = 0;
                for (int k = 0; k < iw; k++)
                  tu[k] += uk * *(it + i++ * step);
                for (int k = 0; k < iw - 1; k++)
                  tv[k] += vk * *(it + i++ * step);
                uk *= u;
                vk *= v;
              }
            hpb[ib] = tu[0];
            for (int k = 1; k < iw; k++)
              hpb[ib] += tu[k] * dlu[k] + tv[k-1] * dlv[k-1];
            vin[iw-2] += bs;
          }
        hpls.insert({iw, hpb});
        dlu.push_back(dlu[1] * dlu.back());
        dlu.push_back(dlv[1] * dlv.back());
      }
    return hpls;
  }
}
