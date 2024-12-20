//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/disbasis.h"
#include "apfel/messages.h"
#include "apfel/operator.h"

#include <numeric>

namespace apfel
{
  //_________________________________________________________________________________
  DISNCBasis::DISNCBasis(int const& k, double const& fact):
    ConvolutionMap{"DISNCBasis_" + std::to_string(k)}
  {
    _rules[GLUON]   = { {CG, GLUON, fact} };
    _rules[SIGMA]   = { {CS, SIGMA, fact / 6} };
    _rules[VALENCE] = { {CS, VALENCE, fact / 6} };
    for (int i = 2; i <= 6; i++)
      {
        double coef = 0;
        if (i == k)
          coef = - 1. / i;
        else if (i >= k + 1)
          coef = 1. / i / ( i - 1 );

        // Change sign to T3 and V3
        if (i == 2)
          coef *= - 1;

        // Multiply the coefficient by the overall factor
        coef *= fact;

        _rules[2 * i - 1] = { {CNS, 2 * i - 1, coef} };
        _rules[2 * i]     = { {CNS, 2 * i,     coef} };
      }
  }

  //_________________________________________________________________________________
  DISNCBasis::DISNCBasis(std::vector<double> const& Ch):
    ConvolutionMap{"DISNCBasis_tot"}
  {
    if (Ch.size() != 6)
      throw std::runtime_error(error("DISNCBasis", "The charge vector must have 6 entries."));

    // Sum of the charges
    const double SumCh = accumulate(Ch.begin(), Ch.end(), 0.);

    _rules[GLUON]   = { {CG, GLUON, SumCh} };
    _rules[SIGMA]   = { {CS, SIGMA, SumCh / 6} };
    _rules[VALENCE] = { {CS, VALENCE, SumCh / 6} };
    for (int j = 2; j <= 6; j++)
      {
        double coef = 0;
        for (int i = 1; i <= j; i++)
          if (i < j)
            coef += Ch[i - 1];
          else
            coef += Ch[i - 1] * ( 1 - j );
        coef /= j * ( j - 1 );

        // Change sign to T3 and V3
        if (j == 2)
          coef *= - 1;

        _rules[2 * j - 1] = { {CNS, 2 * j - 1, coef} };
        _rules[2 * j]     = { {CNS, 2 * j,     coef} };
      }
  }

  //_________________________________________________________________________________
  DISCCBasis::DISCCBasis(int const& l, bool const& Is3, double const& fact):
    ConvolutionMap{"DISCCBasis_" + std::to_string(l) + "_" + std::to_string(Is3)}
  {
    // Retrieve CKM matrix element
    const int i = Vij.at(l).first;
    const int j = Vij.at(l).second;

    _rules[GLUON]   = { {CG, GLUON, fact} };
    _rules[SIGMA]   = { {CS, SIGMA, fact / 6} };
    _rules[VALENCE] = { {CS, VALENCE, fact / 6} };
    for (int k = 2; k <= 6; k++)
      {
        double coefp = 0;
        double coefm = 0;
        if (k >= 2 * j - 1)
          {
            double a = 1;
            if (k == 2 * j - 1)
              a -= k;

            coefp += a;
            coefm += a;
          }
        if (k >= 2 * i)
          {
            double b = 1;
            if (k == 2 * i)
              b -= k;

            coefp += b;
            coefm -= b;
          }

        coefp /= 2 * k * ( k - 1 );
        coefm /= 2 * k * ( k - 1 );

        // Change sign to T3 and V3
        if (k == 2)
          {
            coefp *= - 1;
            coefm *= - 1;
          }

        // Multiply the coefficient by the overall factor
        coefp *= fact;
        coefm *= fact;

        if (Is3)
          {
            _rules[2 * k - 1] = { {CNS, 2 * k - 1, coefm} };
            _rules[2 * k]     = { {CNS, 2 * k,     coefp} };
          }
        else
          {
            _rules[2 * k - 1] = { {CNS, 2 * k - 1, coefp} };
            _rules[2 * k]     = { {CNS, 2 * k,     coefm} };
          }
      }
  }

  //_________________________________________________________________________________
  DISCCBasis::DISCCBasis(std::vector<double> const& CKM, bool const& Is3):
    ConvolutionMap{"DISCCBasis_tot"}
  {
    if (CKM.size() != 9)
      throw std::runtime_error(error("DISCCBasis", "The CKM vector must have 9 entries."));

    _rules[GLUON]   = { {CG, GLUON, 0} };
    _rules[SIGMA]   = { {CS, SIGMA, 0} };
    _rules[VALENCE] = { {CS, VALENCE, 0} };
    for (int k = 3; k <= 12; k++)
      _rules[k] = { {CNS, k, 0} };

    // Fill in rules according to the input CKM matrix
    for (auto i = 1; i <= (int) CKM.size(); i++)
      {
        if (CKM[i - 1] == 0)
          continue;

        // Get basis of the i-th component
        const DISCCBasis basis{i, Is3};

        // Get the rules
        const auto rules = basis.GetRules();

        // Update rules
        for (int k = 0; k <= 12; k++)
          if (rules.count(k) != 0)
            _rules[k][0].coefficient += CKM[i - 1] * rules.at(k)[0].coefficient;
      }
  }
  
  //_________________________________________________________________________________
  DISNCBasis_ACOT::DISNCBasis_ACOT(std::vector<double> const& Ch):
    ConvolutionMap{"DISNCBasis_ACOT_tot"},
    _Ch(Ch)
  {
    if (Ch.size() != 6)
      throw std::runtime_error(error("DISNCBasis_ACOT", "The charge vector must have 6 entries."));

    _rules[0] = {{0,0,1}};
    for(int i=1;i<=6;i++){
      _rules[2*i-1] = {{2*i-1,2*i-1,1}};
      _rules[2*i] = {{2*i,2*i,1}};
    }
    _avCh.resize(6);
    for(int i=1; i<=6; i++){
      for(int j=1; j<=i; j++){
        _avCh[i-1] += Ch[j-1];
      }
      _avCh[i-1] /= (double)i;
    }
  }
  //_________________________________________________________________________________
  std::vector<std::map<int,Operator>> DISNCBasis_ACOT::get_light_operators(bool isPV, std::vector<Operator> gluon, std::vector<std::map<int,Operator>> ns, std::vector<Operator> ps){
    std::map<int,Operator> ClLO;
    std::map<int,Operator> ClNLO;
    std::map<int,Operator> ClNNLO;

    //switches between convolution basis for PV charges (pv=1) or not (pv=0)
    int pv = isPV ? 1 : 0;

    //calc the change of basis explicitly
    //PS only enters from O(alpha_S^2) on

    ClLO.insert({      OGLUON , 3.*_avCh[2]*gluon.at(0)});
    ClNLO.insert({     OGLUON , 3.*_avCh[2]*gluon.at(1)});
    ClNNLO.insert({    OGLUON , 3.*_avCh[2]*gluon.at(2)});
    ClLO.insert({  SIGMA + pv , 0.5*_avCh[2]*ns.at(0).at(3)});
    ClNLO.insert({ SIGMA + pv , 0.5*_avCh[2]*ns.at(1).at(3)});
    ClNNLO.insert({SIGMA + pv , 0.5*_avCh[2]*(ns.at(2).at(3)+3.*ps.at(2))});
    ClLO.insert({     T3 + pv , 0.5*(_Ch[0]-_Ch[1])*ns.at(0).at(3)});
    ClNLO.insert({    T3 + pv , 0.5*(_Ch[0]-_Ch[1])*ns.at(1).at(3)});
    ClNNLO.insert({   T3 + pv , 0.5*(_Ch[0]-_Ch[1])*ns.at(2).at(3)});
    ClLO.insert({     T8 + pv , (_Ch[0]+_Ch[1]-2*_Ch[2])/6.*ns.at(0).at(3)});
    ClNLO.insert({    T8 + pv , (_Ch[0]+_Ch[1]-2*_Ch[2])/6.*ns.at(1).at(3)});
    ClNNLO.insert({   T8 + pv , (_Ch[0]+_Ch[1]-2*_Ch[2])/6.*ns.at(2).at(3)});
    ClLO.insert({    T15 + pv , _avCh[2]/4.*ns.at(0).at(3)});
    ClNLO.insert({   T15 + pv , _avCh[2]/4.*ns.at(1).at(3)});
    ClNNLO.insert({  T15 + pv , _avCh[2]/4.*(ns.at(2).at(3)+3.*ps.at(2))});
    ClLO.insert({    T24 + pv , 3.*_avCh[2]/20.*ns.at(0).at(3)});
    ClNLO.insert({   T24 + pv , 3.*_avCh[2]/20.*ns.at(1).at(3)});
    ClNNLO.insert({  T24 + pv , 3.*_avCh[2]/20.*(ns.at(2).at(3)+3.*ps.at(2))});
    ClLO.insert({    T35 + pv , _avCh[2]/10.*ns.at(0).at(3)});
    ClNLO.insert({   T35 + pv , _avCh[2]/10.*ns.at(1).at(3)});
    ClNNLO.insert({  T35 + pv , _avCh[2]/10.*(ns.at(2).at(3)+3.*ps.at(2))});
    return {ClLO,ClNLO,ClNNLO};
  }

  //_________________________________________________________________________________
  std::vector<std::map<int,Operator>> DISNCBasis_ACOT::get_charm_operators(bool isPV, std::vector<Operator> gluon, std::vector<std::map<int,Operator>> ns, std::vector<Operator> ps){
    std::map<int,Operator> CcLO;
    std::map<int,Operator> CcNLO;
    std::map<int,Operator> CcNNLO;

    //switches between convolution basis for PV charges (pv=1) or not (pv=0)
    int pv = isPV ? 1 : 0;


    //calc the change of basis explicitly
    //PS only enters from O(alpha_S^2) on

    CcLO.insert({       OGLUON , _Ch[3]*gluon.at(0)});
    CcNLO.insert({      OGLUON , _Ch[3]*gluon.at(1)});
    CcNNLO.insert({     OGLUON , _Ch[3]*gluon.at(2)});
    CcLO.insert({   SIGMA + pv , 2.*_avCh[3]/3.*ns.at(0).at(4)-0.5*_avCh[2]*ns.at(0).at(3)});
    CcNLO.insert({  SIGMA + pv , 2.*_avCh[3]/3.*ns.at(1).at(4)-0.5*_avCh[2]*ns.at(1).at(3)});
    CcNNLO.insert({ SIGMA + pv , 2.*_avCh[3]/3.*ns.at(2).at(4)-0.5*_avCh[2]*ns.at(2).at(3)+0.5*(_Ch[3]+4./3.*_avCh[3])*ps.at(2)});
    CcLO.insert({      T3 + pv , (_avCh[0]-_Ch[1])/2.*(ns.at(0).at(4)-ns.at(0).at(3))});
    CcNLO.insert({     T3 + pv , (_avCh[0]-_Ch[1])/2.*(ns.at(1).at(4)-ns.at(1).at(3))});
    CcNNLO.insert({    T3 + pv , (_avCh[0]-_Ch[1])/2.*(ns.at(2).at(4)-ns.at(2).at(3))});
    CcLO.insert({      T8 + pv , (_avCh[1]-_Ch[2])/3.*(ns.at(0).at(4)-ns.at(0).at(3))});
    CcNLO.insert({     T8 + pv , (_avCh[1]-_Ch[2])/3.*(ns.at(1).at(4)-ns.at(1).at(3))});
    CcNNLO.insert({    T8 + pv , (_avCh[1]-_Ch[2])/3.*(ns.at(2).at(4)-ns.at(2).at(3))});
    CcLO.insert({     T15 + pv , (_avCh[2]-_Ch[3])/4.*ns.at(0).at(4)-_avCh[2]/4.*ns.at(0).at(3)});
    CcNLO.insert({    T15 + pv , (_avCh[2]-_Ch[3])/4.*ns.at(1).at(4)-_avCh[2]/4.*ns.at(1).at(3)});
    CcNNLO.insert({   T15 + pv , (_avCh[2]-_Ch[3])/4.*ns.at(2).at(4)-_avCh[2]/4.*ns.at(2).at(3)+(_Ch[3]/4.-_avCh[3])*ps.at(2)});
    CcLO.insert({     T24 + pv , _avCh[3]/5.*ns.at(0).at(4)-3.*_avCh[2]/20.*ns.at(0).at(3)});
    CcNLO.insert({    T24 + pv , _avCh[3]/5.*ns.at(1).at(4)-3.*_avCh[2]/20.*ns.at(1).at(3)});
    CcNNLO.insert({   T24 + pv , _avCh[3]/5.*ns.at(2).at(4)-3.*_avCh[2]/20.*ns.at(2).at(3)+3.*(_Ch[3]+4./3.*_avCh[3])/20.*ps.at(2)});
    CcLO.insert({     T35 + pv , 2.*_avCh[3]/15.*ns.at(0).at(4)-_avCh[2]/10.*ns.at(0).at(3)});
    CcNLO.insert({    T35 + pv , 2.*_avCh[3]/15.*ns.at(1).at(4)-_avCh[2]/10.*ns.at(1).at(3)});
    CcNNLO.insert({   T35 + pv , 2.*_avCh[3]/15.*ns.at(2).at(4)-_avCh[2]/10.*ns.at(2).at(3)+(_Ch[3]+4./3.*_avCh[3])/10.*ps.at(2)});
    return {CcLO,CcNLO,CcNNLO};
  }

  //_________________________________________________________________________________
  std::vector<std::map<int,Operator>> DISNCBasis_ACOT::get_bottom_operators(bool isPV, std::vector<Operator> gluon, std::vector<std::map<int,Operator>> ns, std::vector<Operator> ps){
    std::map<int,Operator> CbLO;
    std::map<int,Operator> CbNLO;
    std::map<int,Operator> CbNNLO;

    //switches between convolution basis for PV charges (pv=1) or not (pv=0)
    int pv = isPV ? 1 : 0;

    //calc the change of basis explicitly
    //PS only enters from O(alpha_S^2) on

    CbLO.insert({       OGLUON , _Ch[4]*gluon.at(0)});
    CbNLO.insert({      OGLUON , _Ch[4]*gluon.at(1)});
    CbNNLO.insert({     OGLUON , _Ch[4]*gluon.at(2)});
    CbLO.insert({   SIGMA + pv , 5.*_avCh[4]/6.*ns.at(0).at(5)-2.*_avCh[3]/3.*ns.at(0).at(4)});
    CbNLO.insert({  SIGMA + pv , 5.*_avCh[4]/6.*ns.at(1).at(5)-2.*_avCh[3]/3.*ns.at(1).at(4)});
    CbNNLO.insert({ SIGMA + pv , 5.*_avCh[4]/6.*ns.at(2).at(5)-2.*_avCh[3]/3.*ns.at(2).at(4)+2.*(_Ch[4]+5./4.*_avCh[4])/3.*ps.at(2)});
    CbLO.insert({      T3 + pv , (_avCh[0]-_Ch[1])/2.*(ns.at(0).at(5)-ns.at(0).at(4))});
    CbNLO.insert({     T3 + pv , (_avCh[0]-_Ch[1])/2.*(ns.at(1).at(5)-ns.at(1).at(4))});
    CbNNLO.insert({    T3 + pv , (_avCh[0]-_Ch[1])/2.*(ns.at(2).at(5)-ns.at(2).at(4))});
    CbLO.insert({      T8 + pv , (_avCh[1]-_Ch[2])/3.*(ns.at(0).at(5)-ns.at(0).at(4))});
    CbNLO.insert({     T8 + pv , (_avCh[1]-_Ch[2])/3.*(ns.at(1).at(5)-ns.at(1).at(4))});
    CbNNLO.insert({    T8 + pv , (_avCh[1]-_Ch[2])/3.*(ns.at(2).at(5)-ns.at(2).at(4))});
    CbLO.insert({     T15 + pv , (_avCh[2]-_Ch[3])/4.*(ns.at(0).at(5)-ns.at(0).at(4))});
    CbNLO.insert({    T15 + pv , (_avCh[2]-_Ch[3])/4.*(ns.at(1).at(5)-ns.at(1).at(4))});
    CbNNLO.insert({   T15 + pv , (_avCh[2]-_Ch[3])/4.*(ns.at(2).at(5)-ns.at(2).at(4))});
    CbLO.insert({     T24 + pv , (_avCh[3]-_Ch[4])/5.*ns.at(0).at(5)-_avCh[3]/5.*ns.at(0).at(4)});
    CbNLO.insert({    T24 + pv , (_avCh[3]-_Ch[4])/5.*ns.at(1).at(5)-_avCh[3]/5.*ns.at(1).at(4)});
    CbNNLO.insert({   T24 + pv , (_avCh[3]-_Ch[4])/5.*ns.at(2).at(5)-_avCh[3]/5.*ns.at(2).at(4)+(_Ch[4]-5.*_avCh[4])/5.*ps.at(2)});
    CbLO.insert({     T35 + pv , _avCh[4]/6.*ns.at(0).at(5)-2.*_avCh[3]/15.*ns.at(0).at(4)});
    CbNLO.insert({    T35 + pv , _avCh[4]/6.*ns.at(1).at(5)-2.*_avCh[3]/15.*ns.at(1).at(4)});
    CbNNLO.insert({   T35 + pv , _avCh[4]/6.*ns.at(2).at(5)-2.*_avCh[3]/15.*ns.at(2).at(4)+2.*(_Ch[4]+5./4.*_avCh[4])/15.*ps.at(2)});
    return {CbLO,CbNLO,CbNNLO};
  }

  //_________________________________________________________________________________
  std::vector<std::map<int,Operator>> DISNCBasis_ACOT::get_top_operators(bool isPV, std::vector<Operator> gluon, std::vector<std::map<int,Operator>> ns, std::vector<Operator> ps){
    std::map<int,Operator> CtLO;
    std::map<int,Operator> CtNLO;
    std::map<int,Operator> CtNNLO;

    //switches between convolution basis for PV charges (pv=1) or not (pv=0)
    int pv = isPV ? 1 : 0;

    //calc the change of basis explicitly
    //PS only enters from O(alpha_S^2) on

    CtLO.insert({       OGLUON , _Ch[5]*gluon.at(0)});
    CtNLO.insert({      OGLUON , _Ch[5]*gluon.at(1)});
    CtNNLO.insert({     OGLUON , _Ch[5]*gluon.at(2)});
    CtLO.insert({   SIGMA + pv , _avCh[5]*ns.at(0).at(6)-5.*_avCh[4]/6.*ns.at(0).at(5)});
    CtNLO.insert({  SIGMA + pv , _avCh[5]*ns.at(1).at(6)-5.*_avCh[4]/6.*ns.at(1).at(5)});
    CtNNLO.insert({ SIGMA + pv , _avCh[5]*ns.at(2).at(6)-5.*_avCh[4]/6.*ns.at(2).at(5)+5.*(_Ch[5]+6./5.*_avCh[5])/6.*ps.at(2)});
    CtLO.insert({      T3 + pv , (_avCh[0]-_Ch[1])/2.*(ns.at(0).at(6)-ns.at(0).at(5))});
    CtNLO.insert({     T3 + pv , (_avCh[0]-_Ch[1])/2.*(ns.at(1).at(6)-ns.at(1).at(5))});
    CtNNLO.insert({    T3 + pv , (_avCh[0]-_Ch[1])/2.*(ns.at(2).at(6)-ns.at(2).at(5))});
    CtLO.insert({      T8 + pv , (_avCh[1]-_Ch[2])/3.*(ns.at(0).at(6)-ns.at(0).at(5))});
    CtNLO.insert({     T8 + pv , (_avCh[1]-_Ch[2])/3.*(ns.at(1).at(6)-ns.at(1).at(5))});
    CtNNLO.insert({    T8 + pv , (_avCh[1]-_Ch[2])/3.*(ns.at(2).at(6)-ns.at(2).at(5))});
    CtLO.insert({     T15 + pv , (_avCh[2]-_Ch[3])/4.*(ns.at(0).at(6)-ns.at(0).at(5))});
    CtNLO.insert({    T15 + pv , (_avCh[2]-_Ch[3])/4.*(ns.at(1).at(6)-ns.at(1).at(5))});
    CtNNLO.insert({   T15 + pv , (_avCh[2]-_Ch[3])/4.*(ns.at(2).at(6)-ns.at(2).at(5))});
    CtLO.insert({     T24 + pv , (_avCh[3]-_Ch[4])/5.*(ns.at(0).at(6)-ns.at(0).at(5))});
    CtNLO.insert({    T24 + pv , (_avCh[3]-_Ch[4])/5.*(ns.at(1).at(6)-ns.at(1).at(5))});
    CtNNLO.insert({   T24 + pv , (_avCh[3]-_Ch[4])/5.*(ns.at(2).at(6)-ns.at(2).at(5))});
    CtLO.insert({     T35 + pv , (_avCh[4]-_Ch[5])/6.*ns.at(0).at(6)-_avCh[4]/6.*ns.at(0).at(5)});
    CtNLO.insert({    T35 + pv , (_avCh[4]-_Ch[5])/6.*ns.at(1).at(6)-_avCh[4]/6.*ns.at(1).at(5)});
    CtNNLO.insert({   T35 + pv , (_avCh[4]-_Ch[5])/6.*ns.at(2).at(6)-_avCh[4]/6.*ns.at(2).at(5)+(_Ch[5]-6.*_avCh[5])/6.*ps.at(2)});
    return {CtLO,CtNLO,CtNNLO};
  }

  //_________________________________________________________________________________
  std::vector<std::map<int,Operator>> DISNCBasis_ACOT::get_tot_operators(bool isPV, std::vector<std::vector<std::map<int,Operator>>> coeff){
    std::map<int,Operator> CLO;
    std::map<int,Operator> CNLO;
    std::map<int,Operator> CNNLO;

    //switches between convolution basis for PV charges (pv=1) or not (pv=0)
    int pv = isPV ? 1 : 0;

    //simply add all coefficients from all contributions for each distribution

    CLO.insert({       OGLUON , coeff[0][0].at(OGLUON)+coeff[1][0].at(OGLUON)+coeff[2][0].at(OGLUON)+coeff[3][0].at(OGLUON)});
    CNLO.insert({      OGLUON , coeff[0][1].at(OGLUON)+coeff[1][1].at(OGLUON)+coeff[2][1].at(OGLUON)+coeff[3][1].at(OGLUON)});
    CNNLO.insert({     OGLUON , coeff[0][2].at(OGLUON)+coeff[1][2].at(OGLUON)+coeff[2][2].at(OGLUON)+coeff[3][2].at(OGLUON)});
    CLO.insert({   SIGMA + pv , coeff[0][0].at(SIGMA + pv)+coeff[1][0].at(SIGMA + pv)+coeff[2][0].at(SIGMA + pv)+coeff[3][0].at(SIGMA + pv)});
    CNLO.insert({  SIGMA + pv , coeff[0][1].at(SIGMA + pv)+coeff[1][1].at(SIGMA + pv)+coeff[2][1].at(SIGMA + pv)+coeff[3][1].at(SIGMA + pv)});
    CNNLO.insert({ SIGMA + pv , coeff[0][2].at(SIGMA + pv)+coeff[1][2].at(SIGMA + pv)+coeff[2][2].at(SIGMA + pv)+coeff[3][2].at(SIGMA + pv)});
    CLO.insert({      T3 + pv , coeff[0][0].at(T3 + pv)+coeff[1][0].at(T3 + pv)+coeff[2][0].at(T3 + pv)+coeff[3][0].at(T3 + pv)});
    CNLO.insert({     T3 + pv , coeff[0][1].at(T3 + pv)+coeff[1][1].at(T3 + pv)+coeff[2][1].at(T3 + pv)+coeff[3][1].at(T3 + pv)});
    CNNLO.insert({    T3 + pv , coeff[0][2].at(T3 + pv)+coeff[1][2].at(T3 + pv)+coeff[2][2].at(T3 + pv)+coeff[3][2].at(T3 + pv)});
    CLO.insert({      T8 + pv , coeff[0][0].at(T8 + pv)+coeff[1][0].at(T8 + pv)+coeff[2][0].at(T8 + pv)+coeff[3][0].at(T8 + pv)});
    CNLO.insert({     T8 + pv , coeff[0][1].at(T8 + pv)+coeff[1][1].at(T8 + pv)+coeff[2][1].at(T8 + pv)+coeff[3][1].at(T8 + pv)});
    CNNLO.insert({    T8 + pv , coeff[0][2].at(T8 + pv)+coeff[1][2].at(T8 + pv)+coeff[2][2].at(T8 + pv)+coeff[3][2].at(T8 + pv)});
    CLO.insert({     T15 + pv , coeff[0][0].at(T15 + pv)+coeff[1][0].at(T15 + pv)+coeff[2][0].at(T15 + pv)+coeff[3][0].at(T15 + pv)});
    CNLO.insert({    T15 + pv , coeff[0][1].at(T15 + pv)+coeff[1][1].at(T15 + pv)+coeff[2][1].at(T15 + pv)+coeff[3][1].at(T15 + pv)});
    CNNLO.insert({   T15 + pv , coeff[0][2].at(T15 + pv)+coeff[1][2].at(T15 + pv)+coeff[2][2].at(T15 + pv)+coeff[3][2].at(T15 + pv)});
    CLO.insert({     T24 + pv , coeff[0][0].at(T24 + pv)+coeff[1][0].at(T24 + pv)+coeff[2][0].at(T24 + pv)+coeff[3][0].at(T24 + pv)});
    CNLO.insert({    T24 + pv , coeff[0][1].at(T24 + pv)+coeff[1][1].at(T24 + pv)+coeff[2][1].at(T24 + pv)+coeff[3][1].at(T24 + pv)});
    CNNLO.insert({   T24 + pv , coeff[0][2].at(T24 + pv)+coeff[1][2].at(T24 + pv)+coeff[2][2].at(T24 + pv)+coeff[3][2].at(T24 + pv)});
    CLO.insert({     T35 + pv , coeff[0][0].at(T35 + pv)+coeff[1][0].at(T35 + pv)+coeff[2][0].at(T35 + pv)+coeff[3][0].at(T35 + pv)});
    CNLO.insert({    T35 + pv , coeff[0][1].at(T35 + pv)+coeff[1][1].at(T35 + pv)+coeff[2][1].at(T35 + pv)+coeff[3][1].at(T35 + pv)});
    CNNLO.insert({   T35 + pv , coeff[0][2].at(T35 + pv)+coeff[1][2].at(T35 + pv)+coeff[2][2].at(T35 + pv)+coeff[3][2].at(T35 + pv)});

    return {CLO,CNLO,CNNLO};
  }
}
