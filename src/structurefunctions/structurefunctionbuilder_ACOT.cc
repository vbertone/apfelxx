/**
 * @file structurefunctionbuilderACOT.cc
 * @author Peter Risse
 * @brief Same as structurefunctionbuilder but for the ACOT prescription
 * @version 1.0
 * @date 2024-07-31 
 */

#include "apfel/structurefunctionbuilder.h"
#include "apfel/structurefunctionbuilder_ACOT.h"
#include "apfel/timer.h"
#include "apfel/tools.h"
#include "apfel/constants.h"
#include "apfel/zeromasscoefficientfunctionsunp_sl.h"
#include "apfel/zeromasscoefficientfunctionspol_sl.h"
#include "apfel/massivecoefficientfunctionsunp_sl.h"
#include "apfel/massivecoefficientfunctionsunp_sl_ACOT.h"
#include "apfel/massivezerocoefficientfunctionsunp_sl.h"
#include "apfel/zeromasscoefficientfunctionsunp_tl.h"
#include "apfel/tabulateobject.h"

#include "apfel/splittingfunctionsunp_sl.h"

#include <numeric>

namespace apfel
{ 

  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF2NCObjectsACOT(Grid                const& g,
                                                                                                                 std::vector<double> const& Masses,
                                                                                                                 double              const& IntEps,
                                                                                                                 int                 const& nQ,
                                                                                                                 double              const& Qmin,
                                                                                                                 double              const& Qmax,
                                                                                                                 int                 const& intdeg)
  {  
    Timer t;
    // Make sure that the vector of masses contains all the 6 masses.
    if (Masses.size() != 6)
      throw std::runtime_error(error("InitializeF2NCObjectsACOT", "The vector of masses has to contain exactly 6 ordered masses."));

    report("Initializing StructureFunctionObjects for F2 NC full ACOT up to NLO.\n");


    // calc gridding range
    const double mc(Masses[3]),mt(Masses[5]);
    const double Qgridmin = 0.9*Qmin;
    const double Qgridmax = 1.1*Qmax;
    const double sximin = Qgridmin/mt;
    const double sximax = Qgridmax/mc;
    const double lambda = 0.99*sximin;

    // Determine number of active flavours
    int actnf = 0;
    for (auto const& m : Masses){
      if (m < eps8){
        actnf++;
      }
    }

    // Initalise DGLAP objects needed for scale variations
    const auto PDFObj = apfel::InitializeDglapObjectsQCD(g, std::vector<double>(actnf, 0.));

    // Zero Mass coefficient functions
    const Operator Zero {g, Null{}, IntEps};
    const Operator Id {g,Identity{}, IntEps};
    const Operator O21q_L{g, C21ns{}, IntEps};
    const Operator O21g_L {g, C21g{},  IntEps};
  // prepare light components
    std::vector<Operator> L_gluon;
    std::vector<Operator> L_ps;
    std::vector<std::map<int,Operator>> L_ns;

    L_gluon.push_back(Zero); 
    L_gluon.push_back(O21g_L);
    L_gluon.push_back(Zero);

    L_ps.push_back(Zero); 
    L_ps.push_back(Zero); 
    L_ps.push_back(Zero);

    L_ns.push_back({{3,Id}});
    L_ns.push_back({{3,O21q_L}});
    L_ns.push_back({{3,Zero}}); 

    // Massive coefficient functions
    const auto fO20ns = [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 2./(1+std::sqrt( 1 + 4 / xi ));
      const Operator O20{g,Cm20qNC_ACOT{eta},IntEps}; 
      return O20;
    };
    const TabulateObject<Operator> TabO20q_H{fO20ns, nQ, sximin, sximax, intdeg, {1.}, lambda};

    // NLO 
    const auto fO21g = [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta_g = 1/( 1 + 4 / xi );
      const double eta_sub = 2./(1+std::sqrt( 1 + 4 / xi ));
      const Operator Om21gNC{g,Cm21gNC_ACOT{eta_g},IntEps};
      const Operator Om21gNC_sub{g,Cm21gNC_sub_ACOT{eta_sub},IntEps};
      return Om21gNC + (xi>=1 ?-log(xi)*Om21gNC_sub : Zero);
    };
    const TabulateObject<Operator> TabO21g_H{fO21g, nQ, sximin, sximax, intdeg, {1.}, lambda};
    
    const auto fO21ns =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 2./(1+std::sqrt( 1 + 4 / xi ));
      const Operator O21nsNC{g,Cm21qNC_ACOT{eta},IntEps};
      const Operator Om21qNC_sub{g,Cm21qNC_sub_ACOT{eta},IntEps};
      return O21nsNC - Om21qNC_sub;
    };
    const TabulateObject<Operator> TabO21q_H{fO21ns, nQ, sximin, sximax, intdeg, {1.}, lambda};


    // Vector of distributions to skip
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

    const auto F2Obj = [=,&g] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects{
      if(Q<Qgridmin || Q>Qgridmax){
        throw std::runtime_error(error("F2NCfullACOT NLO", "Q out of range ["+std::to_string(Qmin)+","+std::to_string(Qmax)+"]. Q=" + std::to_string(Q)));
      }
      StructureFunctionObjects FObj;
      FObj.skip = skip;
      FObj.nf = actnf;
      FObj.P = PDFObj.at(actnf);

      // change charges ->  change the order of up and down quark to match the transformation into the QCD-evolution basis
      std::vector<double> myCh = Ch;
      myCh[0] = Ch[1]; 
      myCh[1] = Ch[0];

      std::vector<std::vector<Operator>> H_gluon(3);
      std::vector<std::vector<Operator>> H_ps(3);
      std::vector<std::vector<std::map<int,Operator>>> H_ns(3);

      // insert heavy quark coefficients
      for(int j = 4; j<=6; j++){
        const double sxi = Q/Masses[j-1];
        H_gluon.at(j-4).push_back(Zero); 
        H_gluon.at(j-4).push_back(Zero+TabO21g_H.Evaluate(sxi));
        H_gluon.at(j-4).push_back(Zero);

        H_ps.at(j-4).push_back(Zero);
        H_ps.at(j-4).push_back(Zero);
        H_ps.at(j-4).push_back(Zero);

        H_ns.at(j-4).push_back({{j-1,Zero+TabO20q_H.Evaluate(sxi)},{j,Zero+TabO20q_H.Evaluate(sxi)}});
        H_ns.at(j-4).push_back({{j-1,Zero+TabO21q_H.Evaluate(sxi)},{j,Zero+TabO21q_H.Evaluate(sxi)}});
        H_ns.at(j-4).push_back({{j-1,Zero},{j,Zero}});
      }

      //construction of the total structure function by calculating light, charm, bottom and top part individually 
      //such that we can access them as well
      DISNCBasis_ACOT DISbasis_L(myCh);
      auto C2_light_coeff = DISbasis_L.get_light_operators(false,L_gluon,L_ns,L_ps);
      FObj.ConvBasis.insert({1,DISbasis_L});
      FObj.C0.insert({1,Set<Operator>{FObj.ConvBasis.at(1), C2_light_coeff[0]}});
      FObj.C1.insert({1,Set<Operator>{FObj.ConvBasis.at(1), C2_light_coeff[1]}});
      FObj.C2.insert({1,Set<Operator>{FObj.ConvBasis.at(1), C2_light_coeff[2]}});

      DISNCBasis_ACOT DISbasis_C(myCh);
      auto C2_charm_coeff = DISbasis_C.get_charm_operators(false,H_gluon.at(0),H_ns.at(0),H_ps.at(0));
      FObj.ConvBasis.insert({2,DISbasis_C});
      FObj.C0.insert({2,Set<Operator>{FObj.ConvBasis.at(2), C2_charm_coeff.at(0)}});
      FObj.C1.insert({2,Set<Operator>{FObj.ConvBasis.at(2), C2_charm_coeff.at(1)}});
      FObj.C2.insert({2,Set<Operator>{FObj.ConvBasis.at(2), C2_charm_coeff.at(2)}});

      DISNCBasis_ACOT DISbasis_B(myCh);
      auto C2_bottom_coeff = DISbasis_B.get_bottom_operators(false,H_gluon.at(1),H_ns.at(1),H_ps.at(1));
      FObj.ConvBasis.insert({3,DISbasis_B});
      FObj.C0.insert({3,Set<Operator>{FObj.ConvBasis.at(3), C2_bottom_coeff.at(0)}});
      FObj.C1.insert({3,Set<Operator>{FObj.ConvBasis.at(3), C2_bottom_coeff.at(1)}});
      FObj.C2.insert({3,Set<Operator>{FObj.ConvBasis.at(3), C2_bottom_coeff.at(2)}});

      DISNCBasis_ACOT DISbasis_T(myCh);
      auto C2_top_coeff = DISbasis_T.get_top_operators(false,H_gluon.at(2),H_ns.at(2),H_ps.at(2));
      FObj.ConvBasis.insert({4,DISbasis_T});
      FObj.C0.insert({4,Set<Operator>{FObj.ConvBasis.at(4), C2_top_coeff.at(0)}});
      FObj.C1.insert({4,Set<Operator>{FObj.ConvBasis.at(4), C2_top_coeff.at(1)}});
      FObj.C2.insert({4,Set<Operator>{FObj.ConvBasis.at(4), C2_top_coeff.at(2)}});

      DISNCBasis_ACOT DISbasis_all(myCh);
      auto C2_tot_coeff = DISbasis_all.get_tot_operators(false,{C2_light_coeff,C2_charm_coeff,C2_bottom_coeff,C2_top_coeff});
      FObj.ConvBasis.insert({0, DISbasis_all});
      FObj.C0.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C2_tot_coeff.at(0)}});
      FObj.C1.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C2_tot_coeff.at(1)}});
      FObj.C2.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C2_tot_coeff.at(2)}});
   
      return FObj;
    };

    t.stop();
    return F2Obj;
  }

std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF2NCObjectsSACOT(Grid                const& g,
                                                                                                                std::vector<double> const& Masses,
                                                                                                                double              const& IntEps,
                                                                                                                int                 const& nQ,
                                                                                                                double              const& Qmin,
                                                                                                                double              const& Qmax,
                                                                                                                int                 const& intdeg)
  {  
    Timer t;
    // Make sure that the vector of masses contains all the 6 masses.
    if (Masses.size() != 6)
      throw std::runtime_error(error("InitializeF2NCObjectsSACOT", "The vector of masses has to contain exactly 6 ordered masses."));
    
    report("Initializing StructureFunctionObjects for F2 NC simplified ACOT up to NLO.\n");


    // calc gridding range
    const double mc(Masses[3]),mt(Masses[5]);
    const double Qgridmin = 0.9*Qmin;
    const double Qgridmax = 1.1*Qmax;
    const double sximin = Qgridmin/mt;
    const double sximax = Qgridmax/mc;
    const double lambda = 0.99*sximin;

    // Determine number of active flavours
    int actnf = 0;
    for (auto const& m : Masses){
      if (m < eps8){
        actnf++;
      }
    }

    // Initalise DGLAP objects needed for scale variations
    const auto PDFObj = apfel::InitializeDglapObjectsQCD(g, std::vector<double>(actnf, 0.));

    // Zero Mass coefficient functions
    const Operator Zero {g, Null{}, IntEps};
    const Operator Id {g,Identity{}, IntEps};
    const Operator O21q_L{g, C21ns{}, IntEps};
    const Operator O21g_L {g, C21g{},  IntEps};
    // prepare light components
    std::vector<Operator> L_gluon;
    std::vector<Operator> L_ps;
    std::vector<std::map<int,Operator>> L_ns;

    L_gluon.push_back(Zero); 
    L_gluon.push_back(O21g_L);
    L_gluon.push_back(Zero);

    L_ps.push_back(Zero); 
    L_ps.push_back(Zero); 
    L_ps.push_back(Zero);

    L_ns.push_back({{3,Id}});
    L_ns.push_back({{3,O21q_L}});
    L_ns.push_back({{3,Zero}});

    // Massive coefficient functions
    const auto fO20ns = [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 4 / xi );
      const Operator O20{g,Cm20qNC_ACOT_chi{eta},IntEps}; 
      return O20;
    };
    const TabulateObject<Operator> TabO20q_H{fO20ns, nQ, sximin, sximax, intdeg, {1.}, lambda};
    
    // NLO 
    const auto fO21g = [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 4 / xi );
      const Operator Om21gNC{g,Cm21gNC_ACOT{eta},IntEps};
      const Operator Om21gNC_sub{g,Cm21gNC_sub_ACOT_chi{eta},IntEps};
      return Om21gNC + (xi>=1 ?-log(xi)*Om21gNC_sub : Zero);
    };
    const TabulateObject<Operator> TabO21g_H{fO21g, nQ, sximin, sximax, intdeg, {1.}, lambda};
    
    const auto fO21ns =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 4 / xi );
      const Operator O21nsNC{g,Cm21qNC_ACOT_chi{eta},IntEps};
      return O21nsNC;
    };
    const TabulateObject<Operator> TabO21q_H{fO21ns, nQ, sximin, sximax, intdeg, {1.}, lambda};

    // Vector of distributions to skip
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

    const auto F2Obj = [=] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects{
      if(Q<Qgridmin || Q>Qgridmax){
        throw std::runtime_error(error("F2NCsimACOT NLO", "Q out of range ["+std::to_string(Qmin)+","+std::to_string(Qmax)+"]. Q=" + std::to_string(Q)));
      }
      StructureFunctionObjects FObj;
      FObj.skip = skip;
      FObj.nf = actnf;
      FObj.P = PDFObj.at(actnf);

      // change charges ->  change the order of up and down quark to match the transformation into the QCD-evolution basis
      std::vector<double> myCh = Ch;
      myCh[0] = Ch[1]; 
      myCh[1] = Ch[0];

      std::vector<std::vector<Operator>> H_gluon(3);
      std::vector<std::vector<Operator>> H_ps(3);
      std::vector<std::vector<std::map<int,Operator>>> H_ns(3);

      // insert heavy quark coefficients
      for(int j = 4; j<=6; j++){
        const double sxi = Q/Masses[j-1];
        H_gluon.at(j-4).push_back(Zero); 
        H_gluon.at(j-4).push_back(Zero+TabO21g_H.Evaluate(sxi));
        H_gluon.at(j-4).push_back(Zero);

        H_ps.at(j-4).push_back(Zero);
        H_ps.at(j-4).push_back(Zero);
        H_ps.at(j-4).push_back(Zero);

        H_ns.at(j-4).push_back({{j-1,Zero+TabO20q_H.Evaluate(sxi)},{j,Zero+TabO20q_H.Evaluate(sxi)}});
        H_ns.at(j-4).push_back({{j-1,Zero+TabO21q_H.Evaluate(sxi)},{j,Zero+TabO21q_H.Evaluate(sxi)}});
        H_ns.at(j-4).push_back({{j-1,Zero},{j,Zero}});
      }

      //construction of the total structure function by calculating light, charm, bottom and top part individually 
      //such that we can access them as well
      DISNCBasis_ACOT DISbasis_L(myCh);
      auto C2_light_coeff = DISbasis_L.get_light_operators(false,L_gluon,L_ns,L_ps);
      FObj.ConvBasis.insert({1,DISbasis_L});
      FObj.C0.insert({1,Set<Operator>{FObj.ConvBasis.at(1), C2_light_coeff[0]}});
      FObj.C1.insert({1,Set<Operator>{FObj.ConvBasis.at(1), C2_light_coeff[1]}});
      FObj.C2.insert({1,Set<Operator>{FObj.ConvBasis.at(1), C2_light_coeff[2]}});

      DISNCBasis_ACOT DISbasis_C(myCh);
      auto C2_charm_coeff = DISbasis_C.get_charm_operators(false,H_gluon.at(0),H_ns.at(0),H_ps.at(0));
      FObj.ConvBasis.insert({2,DISbasis_C});
      FObj.C0.insert({2,Set<Operator>{FObj.ConvBasis.at(2), C2_charm_coeff.at(0)}});
      FObj.C1.insert({2,Set<Operator>{FObj.ConvBasis.at(2), C2_charm_coeff.at(1)}});
      FObj.C2.insert({2,Set<Operator>{FObj.ConvBasis.at(2), C2_charm_coeff.at(2)}});

      DISNCBasis_ACOT DISbasis_B(myCh);
      auto C2_bottom_coeff = DISbasis_B.get_bottom_operators(false,H_gluon.at(1),H_ns.at(1),H_ps.at(1));
      FObj.ConvBasis.insert({3,DISbasis_B});
      FObj.C0.insert({3,Set<Operator>{FObj.ConvBasis.at(3), C2_bottom_coeff.at(0)}});
      FObj.C1.insert({3,Set<Operator>{FObj.ConvBasis.at(3), C2_bottom_coeff.at(1)}});
      FObj.C2.insert({3,Set<Operator>{FObj.ConvBasis.at(3), C2_bottom_coeff.at(2)}});

      DISNCBasis_ACOT DISbasis_T(myCh);
      auto C2_top_coeff = DISbasis_T.get_top_operators(false,H_gluon.at(2),H_ns.at(2),H_ps.at(2));
      FObj.ConvBasis.insert({4,DISbasis_T});
      FObj.C0.insert({4,Set<Operator>{FObj.ConvBasis.at(4), C2_top_coeff.at(0)}});
      FObj.C1.insert({4,Set<Operator>{FObj.ConvBasis.at(4), C2_top_coeff.at(1)}});
      FObj.C2.insert({4,Set<Operator>{FObj.ConvBasis.at(4), C2_top_coeff.at(2)}});

      DISNCBasis_ACOT DISbasis_all(myCh);
      auto C2_tot_coeff = DISbasis_all.get_tot_operators(false,{C2_light_coeff,C2_charm_coeff,C2_bottom_coeff,C2_top_coeff});
      FObj.ConvBasis.insert({0, DISbasis_all});
      FObj.C0.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C2_tot_coeff.at(0)}});
      FObj.C1.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C2_tot_coeff.at(1)}});
      FObj.C2.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C2_tot_coeff.at(2)}});
   
      return FObj;
    };
    t.stop();
    return F2Obj;
  }

  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF3NCObjectsSACOT(Grid                const& g,
                                                                                                                std::vector<double> const& Masses,
                                                                                                                double              const& IntEps,
                                                                                                                int                 const& nQ,
                                                                                                                double              const& Qmin,
                                                                                                                double              const& Qmax,
                                                                                                                int                 const& intdeg)
  {  
    Timer t;

    const int num_flavours = Masses.size();
    if (num_flavours != 6)
      throw std::runtime_error(error("InitializeF3NCObjectsSACOT", "The vector of masses has to contain exactly 6 ordered masses."));

    report("Initializing StructureFunctionObjects for F3 NC simplified ACOT up to NLO\n");

    // calc gridding range
    const double mc(Masses[3]),mt(Masses[5]);
    const double Qgridmin = 0.9*Qmin;
    const double Qgridmax = 1.1*Qmax;
    const double sximin = Qgridmin/mt;
    const double sximax = Qgridmax/mc;
    const double lambda = 0.99*sximin;

    // Determine number of active flavours
    int actnf = 0;
    for (auto const& m : Masses){
      if (m < eps8){
        actnf++;
      }
    }

    // Zero Mass coefficient functions
    const Operator Zero {g, Null{}, IntEps};
    const Operator Id {g,Identity{}, IntEps};
    const Operator O31_ns{g, C31ns{}, IntEps};

    // Massive coefficient functions
    const auto fO30ns = [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 4 / xi );
      const Operator O30{g,Cm20qNC_ACOT_chi{eta},IntEps}; 
      return O30;
    };
    const TabulateObject<Operator> TabO30ns_H{fO30ns, nQ, sximin, sximax, intdeg, {}, lambda};

    const auto fO31ns =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 4 / xi );
      const Operator O31nsNC{g,Cm31qNC_ACOT_chi{eta},IntEps};
      return O31nsNC;
    };
    const TabulateObject<Operator> TabO31ns_H{fO31ns, nQ, sximin, sximax, intdeg, {}, lambda};

    // Vector of distributions to skip
    const std::vector<int> skip = {1, 3, 5, 7, 9, 11};

    const auto F3Obj = [=,&g] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects{
      if(Q<Qgridmin || Q>Qgridmax){
        throw std::runtime_error(error("F3NCsimACOT NLO", "Q out of range ["+std::to_string(Qmin)+","+std::to_string(Qmax)+"]. Q=" + std::to_string(Q)));
      }
      StructureFunctionObjects FObj;
      FObj.skip = skip;

      // change charges ->  change the order of up and down quark to match the transformation into the QCD-evolution basis
      std::vector<double> myCh = Ch;
      myCh[0] = Ch[1]; 
      myCh[1] = Ch[0];

      std::vector<Operator> L_gluon;
      std::vector<Operator> L_ps;
      std::vector<std::map<int,Operator>> L_ns;
      std::vector<std::vector<Operator>> H_gluon(3);
      std::vector<std::vector<Operator>> H_ps(3);
      std::vector<std::vector<std::map<int,Operator>>> H_ns(3);
      // insert light components
      L_gluon.push_back(Zero); //Zero for F3
      L_gluon.push_back(Zero);
      L_gluon.push_back(Zero);

      L_ps.push_back(Zero); //Zero for F3
      L_ps.push_back(Zero); 
      L_ps.push_back(Zero);

      L_ns.push_back({{3,Id}});
      L_ns.push_back({{3,O31_ns}});
      L_ns.push_back({{3,Zero}});
      // insert heavy quark coefficients
      for(int j = 4; j<=6; j++){
        const double sxi = Q/Masses[j-1];
        H_gluon.at(j-4).push_back(Zero); //Zero for F3
        H_gluon.at(j-4).push_back(Zero);
        H_gluon.at(j-4).push_back(Zero);

        H_ps.at(j-4).push_back(Zero); //Zero for F3
        H_ps.at(j-4).push_back(Zero);
        H_ps.at(j-4).push_back(Zero);

        H_ns.at(j-4).push_back({{j-1,TabO30ns_H.Evaluate(sxi)},{j,TabO30ns_H.Evaluate(sxi)}});
        H_ns.at(j-4).push_back({{j-1,TabO31ns_H.Evaluate(sxi)},{j,TabO31ns_H.Evaluate(sxi)}});
        H_ns.at(j-4).push_back({{j-1,Zero},{j,Zero}});
      }
      //construction of the total structure function by calculating light, charm, bottom and top part individually 
      //such that we can access them as well
      DISNCBasis_ACOT DISbasis_L(myCh);
      auto C3_light_coeff = DISbasis_L.get_light_operators(true,L_gluon,L_ns,L_ps);
      FObj.ConvBasis.insert({1,DISbasis_L});
      FObj.C0.insert({1,Set<Operator>{FObj.ConvBasis.at(1), C3_light_coeff[0]}});
      FObj.C1.insert({1,Set<Operator>{FObj.ConvBasis.at(1), C3_light_coeff[1]}});
      FObj.C2.insert({1,Set<Operator>{FObj.ConvBasis.at(1), C3_light_coeff[2]}});

      DISNCBasis_ACOT DISbasis_C(myCh);
      auto C3_charm_coeff = DISbasis_C.get_charm_operators(true,H_gluon.at(0),H_ns.at(0),H_ps.at(0));
      FObj.ConvBasis.insert({2,DISbasis_C});
      FObj.C0.insert({2,Set<Operator>{FObj.ConvBasis.at(2), C3_charm_coeff.at(0)}});
      FObj.C1.insert({2,Set<Operator>{FObj.ConvBasis.at(2), C3_charm_coeff.at(1)}});
      FObj.C2.insert({2,Set<Operator>{FObj.ConvBasis.at(2), C3_charm_coeff.at(2)}});

      DISNCBasis_ACOT DISbasis_B(myCh);
      auto C3_bottom_coeff = DISbasis_B.get_bottom_operators(true,H_gluon.at(1),H_ns.at(1),H_ps.at(1));
      FObj.ConvBasis.insert({3,DISbasis_B});
      FObj.C0.insert({3,Set<Operator>{FObj.ConvBasis.at(3), C3_bottom_coeff.at(0)}});
      FObj.C1.insert({3,Set<Operator>{FObj.ConvBasis.at(3), C3_bottom_coeff.at(1)}});
      FObj.C2.insert({3,Set<Operator>{FObj.ConvBasis.at(3), C3_bottom_coeff.at(2)}});

      DISNCBasis_ACOT DISbasis_T(myCh);
      auto C3_top_coeff = DISbasis_T.get_top_operators(true,H_gluon.at(2),H_ns.at(2),H_ps.at(2));
      FObj.ConvBasis.insert({4,DISbasis_T});
      FObj.C0.insert({4,Set<Operator>{FObj.ConvBasis.at(4), C3_top_coeff.at(0)}});
      FObj.C1.insert({4,Set<Operator>{FObj.ConvBasis.at(4), C3_top_coeff.at(1)}});
      FObj.C2.insert({4,Set<Operator>{FObj.ConvBasis.at(4), C3_top_coeff.at(2)}});

      DISNCBasis_ACOT DISbasis_all(myCh);
      auto C3_tot_coeff = DISbasis_all.get_tot_operators(true,{C3_light_coeff,C3_charm_coeff,C3_bottom_coeff,C3_top_coeff});
      FObj.ConvBasis.insert({0, DISbasis_all});
      FObj.C0.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C3_tot_coeff.at(0)}});
      FObj.C1.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C3_tot_coeff.at(1)}});
      FObj.C2.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C3_tot_coeff.at(2)}});
      
      return FObj;
    };

    t.stop();
    return F3Obj;
  }

  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeFLNCObjectsSACOT(Grid                const& g,
                                                                                                                std::vector<double> const& Masses,
                                                                                                                double              const& IntEps,
                                                                                                                int                 const& nQ,
                                                                                                                double              const& Qmin,
                                                                                                                double              const& Qmax,
                                                                                                                int                 const& intdeg)
  {  
    Timer t;

    const int num_flavours = Masses.size();
    if (num_flavours != 6)
      throw std::runtime_error(error("InitializeFLNCObjectsSACOT", "The vector of masses has to contain exactly 6 ordered masses."));

    report("Initializing StructureFunctionObjects for FL NC simplified ACOT up to NLO\n");

    // calc gridding range
    const double mc(Masses[3]),mt(Masses[5]);
    const double Qgridmin = 0.9*Qmin;
    const double Qgridmax = 1.1*Qmax;
    const double sximin = Qgridmin/mt;
    const double sximax = Qgridmax/mc;
    const double lambda = 0.99*sximin;

    // Determine number of active flavours
    int actnf = 0;
    for (auto const& m : Masses){
      if (m < eps8){
        actnf++;
      }
    }

    // Zero Mass coefficient functions
    const Operator Zero {g, Null{}, IntEps};
    const Operator Id {g,Identity{}, IntEps};
    const Operator OL1_ns{g, CL1ns{}, IntEps};
    const Operator OL1_g{g, CL1g{}, IntEps};

    // Massive coefficient functions
    const auto fOL1ns =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 4 / xi );
      const Operator OL1nsNC{g,CmL1qNC_ACOT_chi{eta},IntEps};
      return OL1nsNC;
    };
    const TabulateObject<Operator> TabOL1ns_H{fOL1ns, nQ, sximin, sximax, intdeg, {}, lambda};

    const auto fOL1g =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 4 / xi );
      const Operator OL1gNC{g,CmL1gNC_ACOT{eta},IntEps};
      return OL1gNC;
    };
    const TabulateObject<Operator> TabOL1g_H{fOL1g, nQ, sximin, sximax, intdeg, {}, lambda};

    // Vector of distributions to skip
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

    const auto FLObj = [=,&g] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects{
      if(Q<Qgridmin || Q>Qgridmax){
        throw std::runtime_error(error("FLNCsimACOT NLO", "Q out of range ["+std::to_string(Qmin)+","+std::to_string(Qmax)+"]. Q=" + std::to_string(Q)));
      }
      StructureFunctionObjects FObj;
      FObj.skip = skip;

      // change charges ->  change the order of up and down quark to match the transformation into the QCD-evolution basis
      std::vector<double> myCh = Ch;
      myCh[0] = Ch[1]; 
      myCh[1] = Ch[0];

      std::vector<Operator> L_gluon;
      std::vector<Operator> L_ps;
      std::vector<std::map<int,Operator>> L_ns;
      std::vector<std::vector<Operator>> H_gluon(3);
      std::vector<std::vector<Operator>> H_ps(3);
      std::vector<std::vector<std::map<int,Operator>>> H_ns(3);
      // insert light components
      L_gluon.push_back(Zero); 
      L_gluon.push_back(OL1_g);
      L_gluon.push_back(Zero);

      L_ps.push_back(Zero); 
      L_ps.push_back(Zero); 
      L_ps.push_back(Zero);

      L_ns.push_back({{3,Zero}}); //Zero at LO
      L_ns.push_back({{3,OL1_ns}});
      L_ns.push_back({{3,Zero}});
      // insert heavy quark coefficients
      for(int j = 4; j<=6; j++){
        const double sxi = Q/Masses[j-1];
        H_gluon.at(j-4).push_back(Zero); 
        H_gluon.at(j-4).push_back(TabOL1g_H.Evaluate(sxi));
        H_gluon.at(j-4).push_back(Zero);

        H_ps.at(j-4).push_back(Zero); 
        H_ps.at(j-4).push_back(Zero);
        H_ps.at(j-4).push_back(Zero);

        H_ns.at(j-4).push_back({{j-1,Zero},{j,Zero}});
        H_ns.at(j-4).push_back({{j-1,TabOL1ns_H.Evaluate(sxi)},{j,TabOL1ns_H.Evaluate(sxi)}});
        H_ns.at(j-4).push_back({{j-1,Zero},{j,Zero}});
      }
      //construction of the total structure function by calculating light, charm, bottom and top part individually 
      //such that we can access them as well
      DISNCBasis_ACOT DISbasis_L(myCh);
      auto CL_light_coeff = DISbasis_L.get_light_operators(false,L_gluon,L_ns,L_ps);
      FObj.ConvBasis.insert({1,DISbasis_L});
      FObj.C0.insert({1,Set<Operator>{FObj.ConvBasis.at(1), CL_light_coeff[0]}});
      FObj.C1.insert({1,Set<Operator>{FObj.ConvBasis.at(1), CL_light_coeff[1]}});
      FObj.C2.insert({1,Set<Operator>{FObj.ConvBasis.at(1), CL_light_coeff[2]}});

      DISNCBasis_ACOT DISbasis_C(myCh);
      auto CL_charm_coeff = DISbasis_C.get_charm_operators(false,H_gluon.at(0),H_ns.at(0),H_ps.at(0));
      FObj.ConvBasis.insert({2,DISbasis_C});
      FObj.C0.insert({2,Set<Operator>{FObj.ConvBasis.at(2), CL_charm_coeff.at(0)}});
      FObj.C1.insert({2,Set<Operator>{FObj.ConvBasis.at(2), CL_charm_coeff.at(1)}});
      FObj.C2.insert({2,Set<Operator>{FObj.ConvBasis.at(2), CL_charm_coeff.at(2)}});

      DISNCBasis_ACOT DISbasis_B(myCh);
      auto CL_bottom_coeff = DISbasis_B.get_bottom_operators(false,H_gluon.at(1),H_ns.at(1),H_ps.at(1));
      FObj.ConvBasis.insert({3,DISbasis_B});
      FObj.C0.insert({3,Set<Operator>{FObj.ConvBasis.at(3), CL_bottom_coeff.at(0)}});
      FObj.C1.insert({3,Set<Operator>{FObj.ConvBasis.at(3), CL_bottom_coeff.at(1)}});
      FObj.C2.insert({3,Set<Operator>{FObj.ConvBasis.at(3), CL_bottom_coeff.at(2)}});

      DISNCBasis_ACOT DISbasis_T(myCh);
      auto CL_top_coeff = DISbasis_T.get_top_operators(false,H_gluon.at(2),H_ns.at(2),H_ps.at(2));
      FObj.ConvBasis.insert({4,DISbasis_T});
      FObj.C0.insert({4,Set<Operator>{FObj.ConvBasis.at(4), CL_top_coeff.at(0)}});
      FObj.C1.insert({4,Set<Operator>{FObj.ConvBasis.at(4), CL_top_coeff.at(1)}});
      FObj.C2.insert({4,Set<Operator>{FObj.ConvBasis.at(4), CL_top_coeff.at(2)}});

      DISNCBasis_ACOT DISbasis_all(myCh);
      auto CL_tot_coeff = DISbasis_all.get_tot_operators(false,{CL_light_coeff,CL_charm_coeff,CL_bottom_coeff,CL_top_coeff});
      FObj.ConvBasis.insert({0, DISbasis_all});
      FObj.C0.insert({0, Set<Operator>{FObj.ConvBasis.at(0), CL_tot_coeff.at(0)}});
      FObj.C1.insert({0, Set<Operator>{FObj.ConvBasis.at(0), CL_tot_coeff.at(1)}});
      FObj.C2.insert({0, Set<Operator>{FObj.ConvBasis.at(0), CL_tot_coeff.at(2)}});
      
      return FObj;
    };

    t.stop();
    return FLObj;
  }

  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF2CCPlusObjectsSACOT(Grid                const& g,
                                                                                                                    std::vector<double> const& Masses,
                                                                                                                    double              const& IntEps,
                                                                                                                    int                 const& nQ,
                                                                                                                    double              const& Qmin,
                                                                                                                    double              const& Qmax,
                                                                                                                    int                 const& intdeg)
  {  
    Timer t;

    const int num_flavours = Masses.size();
    if (num_flavours != 6)
      throw std::runtime_error(error("InitializeF2CCPlusObjectsSACOT", "The vector of masses has to contain exactly 6 ordered masses."));
    
    report("Initializing StructureFunctionObjects for F2 CC plus simplified ACOT up to NLO\n");

    const Operator Zero {g, Null{}, IntEps};

    const double ml = Masses[0] == 0 ? 0.01 : Masses[0];
    const double mc = Masses[3] == 0 ? 0.01 : Masses[3];
    const double mb = Masses[4] == 0 ? 0.01 : Masses[4];
    const double mt = Masses[5] == 0 ? 0.01 : Masses[5];
    const double mll = ml+ml;
    const double ml2 = ml*ml;
    const double mc2 = mc*mc;
    const double mlc = ml+mc;
    const double mcb = mc+mb;
    const double mb2 = mb*mb;
    const double mlb = ml+mb;
    const double mt2 = mt*mt;
    const double mlt = mt+ml;
    const double mtb = mt+mb;
    const std::vector<double> m({ml,ml,ml,mc,mb,mt});

    // calc gridding range
    const double Qgridmin = 0.95*Qmin;
    const double Qgridmax = 1.05*Qmax;
    const double sximin = Qgridmin/(mt+mb);
    const double sximax = Qgridmax/ml;
    const double lambda = 0.99*sximin;

    const auto fO20ns = [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1 / xi );
      const Operator O20{g,Cm20qNC_ACOT_chi{eta},IntEps}; 
      return O20;
    };
    const TabulateObject<Operator> TabO20ns_H{fO20ns, nQ, sximin, sximax, intdeg, {}, lambda};

    const auto fO21ns =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1 / xi );
      const Operator O21nsNC{g,Cm21qNC_ACOT_chi{eta},IntEps};
      return O21nsNC;
    };
    const TabulateObject<Operator> TabO21ns_H{fO21ns, nQ, sximin, sximax, intdeg, {}, lambda};

    const auto fO21g_light_light =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1./ xi );
      const Operator O21CCg_gen_mass{g,Cm21gCC_general_mass{eta,xi,ml,ml},IntEps};
      return O21CCg_gen_mass;
    };
    const TabulateObject<Operator> TabO21g_ll{fO21g_light_light, nQ, Qgridmin/mll, Qgridmax/mll, intdeg, {}, lambda};

    const auto fO21g_light_charm =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1./ xi );
      const Operator O21CCg_gen_mass{g,Cm21gCC_general_mass{eta,xi,ml,mc},IntEps};
      return O21CCg_gen_mass;
    };
    const TabulateObject<Operator> TabO21g_lc{fO21g_light_charm, nQ, Qgridmin/mlc, Qgridmax/mlc, intdeg, {}, lambda};

    const auto fO21g_charm_bottom =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1./ xi );
      const Operator O21CCg_gen_mass{g,Cm21gCC_general_mass{eta,xi,mb,mc},IntEps};
      return O21CCg_gen_mass;
    };
    const TabulateObject<Operator> TabO21g_cb{fO21g_charm_bottom, nQ, Qgridmin/mcb, Qgridmax/mcb, intdeg, {}, lambda};

    const auto fO21g_light_bottom =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1./ xi );
      const Operator O21CCg_gen_mass{g,Cm21gCC_general_mass{eta,xi,ml,mb},IntEps};
      return O21CCg_gen_mass;
    };
    const TabulateObject<Operator> TabO21g_lb{fO21g_light_bottom, nQ, Qgridmin/mlb, Qgridmax/mlb, intdeg, {}, lambda};

    const auto fO21g_top_bottom =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1./ xi );
      const Operator O21CCg_gen_mass{g,Cm21gCC_general_mass{eta,xi,mb,mt},IntEps};
      return O21CCg_gen_mass;
    };
    const TabulateObject<Operator> TabO21g_tb{fO21g_top_bottom, nQ, Qgridmin/mtb, Qgridmax/mtb, intdeg, {}, lambda};

    const auto fO21g_light_top =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1./ xi );
      const Operator O21CCg_gen_mass{g,Cm21gCC_general_mass{eta,xi,ml,mt},IntEps};
      return O21CCg_gen_mass;
    };
    const TabulateObject<Operator> TabO21g_lt{fO21g_light_top, nQ, Qgridmin/mlt, Qgridmax/mlt, intdeg, {}, lambda};

    const auto fO21g_sub =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1./ xi );
      const Operator Om21CCg_sub{g,Cm21gNC_sub_ACOT_chi{eta},IntEps};
      return 0.5*Om21CCg_sub;
    };
    const TabulateObject<Operator> TabO21g_sub{fO21g_sub, nQ, sximin, sximax, intdeg, {}, lambda};

    // Vector of distributions to skip
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

    const auto F2Obj = [=,&g] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects{
      if(Q<Qgridmin || Q>Qgridmax){
        throw std::runtime_error(error("F2CCsimACOT plus", "Q out of range ["+std::to_string(Qmin)+","+std::to_string(Qmax)+"]. Q=" + std::to_string(Q)));
      }
      StructureFunctionObjects FObj;
      FObj.skip = skip;
      const double Q2 = Q*Q;
      const double log_ml = Q>=ml ? log(Q2/ml2) : 0;
      const double log_mc = Q>=mc ? log(Q2/mc2) : 0;
      const double log_mb = Q>=mb ? log(Q2/mb2) : 0;
      const double log_mt = Q>=mt ? log(Q2/mt2) : 0;  

      std::vector<std::map<int,Operator>> coef_tot(3);
      std::vector<std::map<int,Operator>> coef_light(3);
      std::vector<std::map<int,Operator>> coef_charm(3);
      std::vector<std::map<int,Operator>> coef_bottom(3);
      std::vector<std::map<int,Operator>> coef_top(3);

      std::vector<Operator> OqLO({ TabO20ns_H.Evaluate(Q/mll),
                                   TabO20ns_H.Evaluate(Q/mlc),
                                   TabO20ns_H.Evaluate(Q/mcb),
                                   TabO20ns_H.Evaluate(Q/mlb),
                                   TabO20ns_H.Evaluate(Q/mtb),
                                   TabO20ns_H.Evaluate(Q/mlt)});
      std::vector<Operator> OqNLO({TabO21ns_H.Evaluate(Q/mll),   //0
                                   TabO21ns_H.Evaluate(Q/mlc),   //1
                                   TabO21ns_H.Evaluate(Q/mcb),   //2
                                   TabO21ns_H.Evaluate(Q/mlb),   //3
                                   TabO21ns_H.Evaluate(Q/mtb),   //4
                                   TabO21ns_H.Evaluate(Q/mlt)}); //5
      std::vector<std::vector<Operator>> Oq({OqLO,OqNLO});

      std::vector<Operator> OgLO({Zero,Zero,Zero,Zero,Zero,Zero});
      std::vector<Operator> OgNLO({TabO21g_ll.Evaluate(Q/mll) - 2*log_ml       *TabO21g_sub.Evaluate(Q/mll),
                                   TabO21g_lc.Evaluate(Q/mlc) - (log_ml+log_mc)*TabO21g_sub.Evaluate(Q/mlc),
                                   TabO21g_cb.Evaluate(Q/mcb) - (log_mc+log_mb)*TabO21g_sub.Evaluate(Q/mcb),
                                   TabO21g_lb.Evaluate(Q/mlb) - (log_ml+log_mb)*TabO21g_sub.Evaluate(Q/mlb),
                                   TabO21g_tb.Evaluate(Q/mtb) - (log_mt+log_mb)*TabO21g_sub.Evaluate(Q/mtb),
                                   TabO21g_lt.Evaluate(Q/mlt) - (log_ml+log_mt)*TabO21g_sub.Evaluate(Q/mlt)});
      std::vector<std::vector<Operator>> Og({OgLO,OgNLO});

      //total
      for(int l=0; l<2; l++){ // gluon pieces
        coef_tot.at(l).insert({0,(Ch.at(0)+Ch.at(1))*Og.at(l).at(0)}); 
        coef_tot.at(l).at(0) += (Ch.at(3)+Ch.at(4))*Og.at(l).at(1);
        coef_tot.at(l).at(0) += Ch.at(5)*Og.at(l).at(2);
        coef_tot.at(l).at(0) += Ch.at(2)*Og.at(l).at(3); 
        coef_tot.at(l).at(0) += Ch.at(8)*Og.at(l).at(4);
        coef_tot.at(l).at(0) += (Ch.at(6)+Ch.at(7))*Og.at(l).at(5);
      }
      for(int l=0; l<2;l++){ // quark pieces
        coef_tot.at(l).insert({1,Ch.at(0)*Oq.at(l).at(0)+Ch.at(3)*Oq.at(l).at(1)+Ch.at(6)*Oq.at(l).at(5)});
        coef_tot.at(l).insert({2,(Ch.at(0)+Ch.at(1))*Oq.at(l).at(0)+Ch.at(2)*Oq.at(l).at(3)});
        coef_tot.at(l).insert({3,Ch.at(1)*Oq.at(l).at(0)+Ch.at(4)*Oq.at(l).at(1)+Ch.at(7)*Oq.at(l).at(5)});
        coef_tot.at(l).insert({4,(Ch.at(3)+Ch.at(4))*Oq.at(l).at(1)+Ch.at(5)*Oq.at(l).at(2)});
        coef_tot.at(l).insert({5,Ch.at(2)*Oq.at(l).at(3)+Ch.at(5)*Oq.at(l).at(2)+Ch.at(8)*Oq.at(l).at(4)});
        coef_tot.at(l).insert({6,(Ch.at(6)+Ch.at(7))*Oq.at(l).at(5)+Ch.at(8)*Oq.at(l).at(4)});
      }
      for(int i=0;i<=6;i++){
        coef_tot.at(2).insert({i,Zero}); 
      }
      //light
      for(int k=0; k<2; k++){ //gluon pieces
        coef_light.at(k).insert({0,(Ch.at(0)+Ch.at(1))*Og.at(k).at(0)}); 
        coef_light.at(k).at(0) += (Ch.at(3)+Ch.at(4))*Og.at(k).at(1);
        coef_light.at(k).at(0) += Ch.at(2)*Og.at(k).at(3); 
        coef_light.at(k).at(0) += (Ch.at(6)+Ch.at(7))*Og.at(k).at(5);
      }
      for(int k=0; k<2;k++){ //quark pieces
        coef_light.at(k).insert({1,Ch.at(0)*Oq.at(k).at(0)+Ch.at(3)*Oq.at(k).at(1)+Ch.at(6)*Oq.at(k).at(5)});
        coef_light.at(k).insert({2,(Ch.at(0)+Ch.at(1))*Oq.at(k).at(0)+Ch.at(2)*Oq.at(k).at(3)});
        coef_light.at(k).insert({3,Ch.at(1)*Oq.at(k).at(0)+Ch.at(4)*Oq.at(k).at(1)+Ch.at(7)*Oq.at(k).at(5)});
        coef_light.at(k).insert({4,(Ch.at(3)+Ch.at(4))*Oq.at(k).at(1)});
        coef_light.at(k).insert({5,Ch.at(2)*Oq.at(k).at(3)});
        coef_light.at(k).insert({6,(Ch.at(6)+Ch.at(7))*Oq.at(k).at(5)});
      }
      for(int i=0;i<=6;i++){
        coef_light.at(2).insert({i,Zero}); 
      }
      //charm
      for(int k=0; k<2; k++){
        coef_charm.at(k).insert({0,(Ch.at(3)+Ch.at(4))*Og.at(k).at(1)});
        coef_charm.at(k).at(0) += Ch.at(5)*Og.at(k).at(2);
      }
      for(int k=0; k<2; k++){
        coef_charm.at(k).insert({1,Ch.at(3)*Oq.at(k).at(1)});
        coef_charm.at(k).insert({2,Zero});
        coef_charm.at(k).insert({3,Ch.at(4)*Oq.at(k).at(1)});
        coef_charm.at(k).insert({4,(Ch.at(4)+Ch.at(3))*Oq.at(k).at(1)+Ch.at(5)*Oq.at(k).at(2)});
        coef_charm.at(k).insert({5,Ch.at(5)*Oq.at(k).at(2)});
        coef_charm.at(k).insert({6,Zero});
      }
      for(int i=0;i<=6;i++){
        coef_charm.at(2).insert({i,Zero}); 
      }
      //bottom
      for(int k=0; k<2; k++){ // gluon pieces
        coef_bottom.at(k).insert({0,Ch.at(5)*Og.at(k).at(2)});
        coef_bottom.at(k).at(0) += Ch.at(2)*Og.at(k).at(3); 
        coef_bottom.at(k).at(0) += Ch.at(8)*Og.at(k).at(4);
      }
      for(int k=0; k<2;k++){ // quark pieces
        coef_bottom.at(k).insert({1,Zero});
        coef_bottom.at(k).insert({2,Ch.at(2)*Oq.at(k).at(3)});
        coef_bottom.at(k).insert({3,Zero});
        coef_bottom.at(k).insert({4,Ch.at(5)*Oq.at(k).at(2)});
        coef_bottom.at(k).insert({5,Ch.at(2)*Oq.at(k).at(3)+Ch.at(5)*Oq.at(k).at(2)+Ch.at(8)*Oq.at(k).at(4)});
        coef_bottom.at(k).insert({6,Ch.at(8)*Oq.at(k).at(4)});
      }
      for(int i=0;i<=6;i++){
        coef_bottom.at(2).insert({i,Zero}); 
      }
      //top
      for(int k=0; k<2; k++){ // gluon pieces
        coef_top.at(k).insert({0,Ch.at(8)*Og.at(k).at(4)});
        coef_top.at(k).at(0) += (Ch.at(6)+Ch.at(7))*Og.at(k).at(5);
      }
      for(int k=0; k<2;k++){ // quark pieces
        coef_top.at(k).insert({1,Ch.at(6)*Oq.at(k).at(5)});
        coef_top.at(k).insert({2,Zero});
        coef_top.at(k).insert({3,Ch.at(7)*Oq.at(k).at(5)});
        coef_top.at(k).insert({4,Zero});
        coef_top.at(k).insert({5,Ch.at(8)*Oq.at(k).at(4)});
        coef_top.at(k).insert({6,(Ch.at(6)+Ch.at(7))*Oq.at(k).at(5)+Ch.at(8)*Oq.at(k).at(4)});
      }
      for(int i=0;i<=6;i++){
        coef_top.at(2).insert({i,Zero}); 
      }
      
      DISCCBasis_ACOT DISbasis_all(Ch,Zero);
      auto C_tot = DISbasis_all.get_operators_plus(coef_tot);
      FObj.ConvBasis.insert({0, DISbasis_all});
      FObj.C0.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C_tot.at(0)}});
      FObj.C1.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C_tot.at(1)}});
      FObj.C2.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C_tot.at(2)}});

      DISCCBasis_ACOT DISbasis_light(Ch,Zero);
      auto C_light = DISbasis_light.get_operators_plus(coef_light);
      FObj.ConvBasis.insert({1, DISbasis_light});
      FObj.C0.insert({1, Set<Operator>{FObj.ConvBasis.at(1), C_light.at(0)}});
      FObj.C1.insert({1, Set<Operator>{FObj.ConvBasis.at(1), C_light.at(1)}});
      FObj.C2.insert({1, Set<Operator>{FObj.ConvBasis.at(1), C_light.at(2)}});

      DISCCBasis_ACOT DISbasis_charm(Ch,Zero);
      auto C_charm = DISbasis_charm.get_operators_plus(coef_charm);
      FObj.ConvBasis.insert({2, DISbasis_charm});
      FObj.C0.insert({2, Set<Operator>{FObj.ConvBasis.at(2), C_charm.at(0)}});
      FObj.C1.insert({2, Set<Operator>{FObj.ConvBasis.at(2), C_charm.at(1)}});
      FObj.C2.insert({2, Set<Operator>{FObj.ConvBasis.at(2), C_charm.at(2)}});

      DISCCBasis_ACOT DISbasis_bottom(Ch,Zero);
      auto C_bottom = DISbasis_bottom.get_operators_plus(coef_bottom);
      FObj.ConvBasis.insert({3, DISbasis_bottom});
      FObj.C0.insert({3, Set<Operator>{FObj.ConvBasis.at(3), C_bottom.at(0)}});
      FObj.C1.insert({3, Set<Operator>{FObj.ConvBasis.at(3), C_bottom.at(1)}});
      FObj.C2.insert({3, Set<Operator>{FObj.ConvBasis.at(3), C_bottom.at(2)}});

      DISCCBasis_ACOT DISbasis_top(Ch,Zero);
      auto C_top = DISbasis_top.get_operators_plus(coef_top);
      FObj.ConvBasis.insert({4, DISbasis_top});
      FObj.C0.insert({4, Set<Operator>{FObj.ConvBasis.at(4), C_top.at(0)}});
      FObj.C1.insert({4, Set<Operator>{FObj.ConvBasis.at(4), C_top.at(1)}});
      FObj.C2.insert({4, Set<Operator>{FObj.ConvBasis.at(4), C_top.at(2)}});

      
      return FObj;
    };

    t.stop();
    return F2Obj;
  }

  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF2CCMinusObjectsSACOT(Grid                const& g,
                                                                                                                     std::vector<double> const& Masses,
                                                                                                                     double              const& IntEps,
                                                                                                                     int                 const& nQ,
                                                                                                                     double              const& Qmin,
                                                                                                                     double              const& Qmax,
                                                                                                                     int                 const& intdeg)
  {  
    Timer t;

    const int num_flavours = Masses.size();
    if (num_flavours != 6)
      throw std::runtime_error(error("InitializeF2CCMinusObjectsSACOT", "The vector of masses has to contain exactly 6 ordered masses."));
    
    report("Initializing StructureFunctionObjects for F2 CC minus simplified ACOT up to NLO\n");

    const Operator Zero {g, Null{}, IntEps};

    const double ml = Masses[0] == 0 ? 0.01 : Masses[0];
    const double mc = Masses[3] == 0 ? 0.01 : Masses[3];
    const double mb = Masses[4] == 0 ? 0.01 : Masses[4];
    const double mt = Masses[5] == 0 ? 0.01 : Masses[5];
    const double mll = ml+ml;
    const double mlc = ml+mc;
    const double mcb = mc+mb;
    const double mlb = ml+mb;
    const double mlt = mt+ml;
    const double mtb = mt+mb;
    const std::vector<double> m({ml,ml,ml,mc,mb,mt});

    // calc gridding range
    const double Qgridmin = 0.95*Qmin;
    const double Qgridmax = 1.05*Qmax;
    const double sximin = Qgridmin/(mt+mb);
    const double sximax = Qgridmax/ml;
    const double lambda = 0.99*sximin;


    // Massive coefficient functions
    const auto fO20ns = [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1 / xi );
      const Operator O20{g,Cm20qNC_ACOT_chi{eta},IntEps}; 
      return O20;
      // return Zero;
    };
    const TabulateObject<Operator> TabO20ns_H{fO20ns, nQ, sximin, sximax, intdeg, {}, lambda};

    const auto fO21ns =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1 / xi );
      const Operator O21nsNC{g,Cm21qNC_ACOT_chi{eta},IntEps};
      return O21nsNC;
      // return Zero;
    };
    const TabulateObject<Operator> TabO21ns_H{fO21ns, nQ, sximin, sximax, intdeg, {}, lambda};

    // Vector of distributions to skip
    const std::vector<int> skip = {0, 1, 2, 3, 5, 7, 9, 11};

    const auto F2Obj = [=,&g] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects{
      if(Q<Qgridmin || Q>Qgridmax){
        throw std::runtime_error(error("F2CCsimACOT NLO minus", "Q out of range ["+std::to_string(Qmin)+","+std::to_string(Qmax)+"]. Q=" + std::to_string(Q)));
      }
      StructureFunctionObjects FObj;
      FObj.skip = skip;
      std::vector<Operator> OqLO({ TabO20ns_H.Evaluate(Q/mll),
                                   TabO20ns_H.Evaluate(Q/mlc),
                                   TabO20ns_H.Evaluate(Q/mcb),
                                   TabO20ns_H.Evaluate(Q/mlb),
                                   TabO20ns_H.Evaluate(Q/mtb),
                                   TabO20ns_H.Evaluate(Q/mlt)});
      std::vector<Operator> OqNLO({TabO21ns_H.Evaluate(Q/mll),   //0
                                   TabO21ns_H.Evaluate(Q/mlc),   //1
                                   TabO21ns_H.Evaluate(Q/mcb),   //2
                                   TabO21ns_H.Evaluate(Q/mlb),   //3
                                   TabO21ns_H.Evaluate(Q/mtb),   //4
                                   TabO21ns_H.Evaluate(Q/mlt)}); //5
      std::vector<std::vector<Operator>> Oq({OqLO,OqNLO});

      std::vector<std::map<int,Operator>> coef_tot(3);
      std::vector<std::map<int,Operator>> coef_light(3);
      std::vector<std::map<int,Operator>> coef_charm(3);
      std::vector<std::map<int,Operator>> coef_bottom(3);
      std::vector<std::map<int,Operator>> coef_top(3);

      //tot
      for(int k=0; k<3; k++){//gluon pieces
        coef_tot.at(k).insert({0,Zero});
      }
      for(int k=0; k<2;k++){// quark pieces
        coef_tot.at(k).insert({1,Ch.at(0)*Oq.at(k).at(0)+Ch.at(3)*Oq.at(k).at(1)+Ch.at(6)*Oq.at(k).at(5)});
        coef_tot.at(k).insert({2,-(Ch.at(0)+Ch.at(1))*Oq.at(k).at(0)-Ch.at(2)*Oq.at(k).at(3)});
        coef_tot.at(k).insert({3,Ch.at(1)*Oq.at(k).at(0)+Ch.at(4)*Oq.at(k).at(1)+Ch.at(7)*Oq.at(k).at(5)});
        coef_tot.at(k).insert({4,-(Ch.at(4)+Ch.at(3))*Oq.at(k).at(1)-Ch.at(5)*Oq.at(k).at(2)});
        coef_tot.at(k).insert({5,Ch.at(2)*Oq.at(k).at(3)+Ch.at(5)*Oq.at(k).at(2)+Ch.at(8)*Oq.at(k).at(4)});
        coef_tot.at(k).insert({6,-(Ch.at(6)+Ch.at(7))*Oq.at(k).at(5)-Ch.at(8)*Oq.at(k).at(4)});
      }
      for(int i=1;i<=6;i++){
        coef_tot.at(2).insert({i,Zero}); 
      }
      //light
      for(int k=0; k<3; k++){ //gluon pieces
        coef_light.at(k).insert({0,Zero});
      }
      for(int k=0; k<2;k++){ //quark pieces
        coef_light.at(k).insert({1,Ch.at(0)*Oq.at(k).at(0)+Ch.at(3)*Oq.at(k).at(1)+Ch.at(6)*Oq.at(k).at(5)});
        coef_light.at(k).insert({2,-(Ch.at(0)+Ch.at(1))*Oq.at(k).at(0)-Ch.at(2)*Oq.at(k).at(3)});
        coef_light.at(k).insert({3,Ch.at(1)*Oq.at(k).at(0)+Ch.at(4)*Oq.at(k).at(1)+Ch.at(7)*Oq.at(k).at(5)});
        coef_light.at(k).insert({4,-(Ch.at(3)+Ch.at(4))*Oq.at(k).at(1)});
        coef_light.at(k).insert({5,Ch.at(2)*Oq.at(k).at(3)});
        coef_light.at(k).insert({6,-(Ch.at(6)+Ch.at(7))*Oq.at(k).at(5)});
      }
      for(int i=1;i<=6;i++){
        coef_light.at(2).insert({i,Zero}); 
      }
      //charm
      for(int k=0; k<3; k++){ //gluon pieces
        coef_charm.at(k).insert({0,Zero});
      }
      //down,up
      for(int k=0; k<2; k++){
        coef_charm.at(k).insert({1,Ch.at(3)*Oq.at(k).at(1)});
        coef_charm.at(k).insert({2,Zero});
        coef_charm.at(k).insert({3,Ch.at(4)*Oq.at(k).at(1)});
        coef_charm.at(k).insert({4,-(Ch.at(4)+Ch.at(3))*Oq.at(k).at(1)-Ch.at(5)*Oq.at(k).at(2)});
        coef_charm.at(k).insert({5,Ch.at(5)*Oq.at(k).at(2)});
        coef_charm.at(k).insert({6,Zero});
      }
      for(int i=1;i<=6;i++){
        coef_charm.at(2).insert({i,Zero}); 
      }
      //bottom
      for(int k=0; k<3; k++){ // gluon pieces
        coef_bottom.at(k).insert({0,Zero});
      }
      for(int k=0; k<2;k++){ // quark pieces
        coef_bottom.at(k).insert({1,Zero});
        coef_bottom.at(k).insert({2,-Ch.at(2)*Oq.at(k).at(3)});
        coef_bottom.at(k).insert({3,Zero});
        coef_bottom.at(k).insert({4,-Ch.at(5)*Oq.at(k).at(2)});
        coef_bottom.at(k).insert({5,Ch.at(2)*Oq.at(k).at(3)+Ch.at(5)*Oq.at(k).at(2)+Ch.at(8)*Oq.at(k).at(4)});
        coef_bottom.at(k).insert({6,-Ch.at(8)*Oq.at(k).at(4)});
      }
      for(int i=1;i<=6;i++){
        coef_bottom.at(2).insert({i,Zero}); 
      }
      //top
      for(int k=0; k<3; k++){ // gluon pieces
        coef_top.at(k).insert({0,Zero});
      }
      for(int k=0; k<2;k++){ // quark pieces
        coef_top.at(k).insert({1,Ch.at(6)*Oq.at(k).at(5)});
        coef_top.at(k).insert({2,Zero});
        coef_top.at(k).insert({3,Ch.at(7)*Oq.at(k).at(5)});
        coef_top.at(k).insert({4,Zero});
        coef_top.at(k).insert({5,Ch.at(8)*Oq.at(k).at(4)});
        coef_top.at(k).insert({6,-(Ch.at(6)+Ch.at(7))*Oq.at(k).at(5)-Ch.at(8)*Oq.at(k).at(4)});
      }
      for(int i=1;i<=6;i++){
        coef_top.at(2).insert({i,Zero}); 
      }

      
      DISCCBasis_ACOT DISbasis_all(Ch,Zero);
      auto C_tot = DISbasis_all.get_operators_minus(coef_tot);
      FObj.ConvBasis.insert({0, DISbasis_all});
      FObj.C0.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C_tot.at(0)}});
      FObj.C1.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C_tot.at(1)}});
      FObj.C2.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C_tot.at(2)}});

      DISCCBasis_ACOT DISbasis_light(Ch,Zero);
      auto C_light = DISbasis_light.get_operators_minus(coef_light);
      FObj.ConvBasis.insert({1, DISbasis_light});
      FObj.C0.insert({1, Set<Operator>{FObj.ConvBasis.at(1), C_light.at(0)}});
      FObj.C1.insert({1, Set<Operator>{FObj.ConvBasis.at(1), C_light.at(1)}});
      FObj.C2.insert({1, Set<Operator>{FObj.ConvBasis.at(1), C_light.at(2)}});

      DISCCBasis_ACOT DISbasis_charm(Ch,Zero);
      auto C_charm = DISbasis_charm.get_operators_minus(coef_charm);
      FObj.ConvBasis.insert({2, DISbasis_charm});
      FObj.C0.insert({2, Set<Operator>{FObj.ConvBasis.at(2), C_charm.at(0)}});
      FObj.C1.insert({2, Set<Operator>{FObj.ConvBasis.at(2), C_charm.at(1)}});
      FObj.C2.insert({2, Set<Operator>{FObj.ConvBasis.at(2), C_charm.at(2)}});

      DISCCBasis_ACOT DISbasis_bottom(Ch,Zero);
      auto C_bottom = DISbasis_bottom.get_operators_minus(coef_bottom);
      FObj.ConvBasis.insert({3, DISbasis_bottom});
      FObj.C0.insert({3, Set<Operator>{FObj.ConvBasis.at(3), C_bottom.at(0)}});
      FObj.C1.insert({3, Set<Operator>{FObj.ConvBasis.at(3), C_bottom.at(1)}});
      FObj.C2.insert({3, Set<Operator>{FObj.ConvBasis.at(3), C_bottom.at(2)}});

      DISCCBasis_ACOT DISbasis_top(Ch,Zero);
      auto C_top = DISbasis_top.get_operators_minus(coef_top);
      FObj.ConvBasis.insert({4, DISbasis_top});
      FObj.C0.insert({4, Set<Operator>{FObj.ConvBasis.at(4), C_top.at(0)}});
      FObj.C1.insert({4, Set<Operator>{FObj.ConvBasis.at(4), C_top.at(1)}});
      FObj.C2.insert({4, Set<Operator>{FObj.ConvBasis.at(4), C_top.at(2)}});

      
      return FObj;
    };

    t.stop();
    return F2Obj;
  }

  std::function<StructureFunctionObjects(double const &, std::vector<double> const &)> InitializeFLCCPlusObjectsSACOT(Grid const &g,
                                                                                                                      std::vector<double> const &Masses,
                                                                                                                      double const &IntEps,
                                                                                                                      int const& nQ,
                                                                                                                      double const& Qmin,
                                                                                                                      double const& Qmax,
                                                                                                                      int    const& intdeg)
  {
    Timer t;

    const int num_flavours = Masses.size();
    if (num_flavours != 6)
      throw std::runtime_error(error("InitializeFLCCPlusObjectsSACOT", "The vector of masses has to contain exactly 6 ordered masses."));
    
    report("Initializing StructureFunctionObjects for FL CC plus simplified ACOT up to NLO\n");

    const Operator Zero {g, Null{}, IntEps};

    const double ml = Masses[0] == 0 ? 0.01 : Masses[0];
    const double mc = Masses[3] == 0 ? 0.01 : Masses[3];
    const double mb = Masses[4] == 0 ? 0.01 : Masses[4];
    const double mt = Masses[5] == 0 ? 0.01 : Masses[5];
    const double mll = ml+ml;
    const double mlc = ml+mc;
    const double mcb = mc+mb;
    const double mlb = ml+mb;
    const double mlt = mt+ml;
    const double mtb = mt+mb;
    const std::vector<double> m({ml,ml,ml,mc,mb,mt});

    // calc grid range
    const double Qgridmin = 0.95*Qmin;
    const double Qgridmax = 1.05*Qmax;
    const double sximin = Qgridmin/(mt+mb);
    const double sximax = Qgridmax/ml;
    const double lambda = 0.99*sximin;

    const auto fOL1ns =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1 / xi );
      const Operator OL1nsNC{g,CmL1qNC_ACOT_chi{eta},IntEps};
      return OL1nsNC;
    };
    const TabulateObject<Operator> TabOL1ns_H{fOL1ns, nQ, sximin, sximax, intdeg, {}, lambda};

    const auto fOL1g_light_light =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1./ xi );
      const Operator OL1CCg_gen_mass{g,CmL1gCC_general_mass{eta,xi,ml,ml},IntEps};
      return OL1CCg_gen_mass;
    };
    const TabulateObject<Operator> TabOL1g_ll{fOL1g_light_light, nQ, Qgridmin/mll, Qgridmax/mll, intdeg, {}, lambda};

    const auto fOL1g_light_charm =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1./ xi );
      const Operator OL1CCg_gen_mass{g,CmL1gCC_general_mass{eta,xi,ml,mc},IntEps};
      return OL1CCg_gen_mass;
    };
    const TabulateObject<Operator> TabOL1g_lc{fOL1g_light_charm, nQ, Qgridmin/mlc, Qgridmax/mlc, intdeg, {}, lambda};

    const auto fOL1g_charm_bottom =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1./ xi );
      const Operator OL1CCg_gen_mass{g,CmL1gCC_general_mass{eta,xi,mb,mc},IntEps};
      return OL1CCg_gen_mass;
    };
    const TabulateObject<Operator> TabOL1g_cb{fOL1g_charm_bottom, nQ, Qgridmin/mcb, Qgridmax/mcb, intdeg, {}, lambda};

    const auto fOL1g_light_bottom =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1./ xi );
      const Operator OL1CCg_gen_mass{g,CmL1gCC_general_mass{eta,xi,ml,mb},IntEps};
      return OL1CCg_gen_mass;
    };
    const TabulateObject<Operator> TabOL1g_lb{fOL1g_light_bottom, nQ, Qgridmin/mlb, Qgridmax/mlb, intdeg, {}, lambda};

    const auto fOL1g_top_bottom =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1./ xi );
      const Operator OL1CCg_gen_mass{g,CmL1gCC_general_mass{eta,xi,mb,mt},IntEps};
      return OL1CCg_gen_mass;
    };
    const TabulateObject<Operator> TabOL1g_tb{fOL1g_top_bottom, nQ, Qgridmin/mtb, Qgridmax/mtb, intdeg, {}, lambda};

    const auto fOL1g_light_top =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1./ xi );
      const Operator OL1CCg_gen_mass{g,CmL1gCC_general_mass{eta,xi,ml,mt},IntEps};
      return OL1CCg_gen_mass;
    };
    const TabulateObject<Operator> TabOL1g_lt{fOL1g_light_top, nQ, Qgridmin/mlt, Qgridmax/mlt, intdeg, {}, lambda};

    // Vector of distributions to skip
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

    const auto FLObj = [=,&g] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects{
      if(Q<Qgridmin || Q>Qgridmax){
      throw std::runtime_error(error("FLCCsimACOT plus", "Q out of range ["+std::to_string(Qmin)+","+std::to_string(Qmax)+"]. Q=" + std::to_string(Q)));
    }
    StructureFunctionObjects FObj;
    FObj.skip = skip;

    std::vector<Operator> OqLO( {Zero,Zero,Zero,Zero,Zero,Zero});
    std::vector<Operator> OqNLO({TabOL1ns_H.Evaluate(Q/mll),   //0
                                 TabOL1ns_H.Evaluate(Q/mlc),   //1
                                 TabOL1ns_H.Evaluate(Q/mcb),   //2
                                 TabOL1ns_H.Evaluate(Q/mlb),   //3
                                 TabOL1ns_H.Evaluate(Q/mtb),   //4
                                 TabOL1ns_H.Evaluate(Q/mlt)}); //5
    std::vector<std::vector<Operator>> Oq({OqLO,OqNLO});

    std::vector<Operator> OgLO( {Zero,Zero,Zero,Zero,Zero,Zero});
    std::vector<Operator> OgNLO({TabOL1g_ll.Evaluate(Q/mll),
                                 TabOL1g_lc.Evaluate(Q/mlc),
                                 TabOL1g_cb.Evaluate(Q/mcb),
                                 TabOL1g_lb.Evaluate(Q/mlb),
                                 TabOL1g_tb.Evaluate(Q/mtb),
                                 TabOL1g_lt.Evaluate(Q/mlt)});
    std::vector<std::vector<Operator>> Og({OgLO,OgNLO});

      std::vector<std::map<int,Operator>> coef_tot(3);
      std::vector<std::map<int,Operator>> coef_light(3);
      std::vector<std::map<int,Operator>> coef_charm(3);
      std::vector<std::map<int,Operator>> coef_bottom(3);
      std::vector<std::map<int,Operator>> coef_top(3);

      //total
      for(int k=0; k<2; k++){ // gluon pieces
        coef_tot.at(k).insert({0,(Ch.at(0)+Ch.at(1))*Og.at(k).at(0)}); 
        coef_tot.at(k).at(0) += (Ch.at(3)+Ch.at(4))*Og.at(k).at(1);
        coef_tot.at(k).at(0) += Ch.at(5)*Og.at(k).at(2);
        coef_tot.at(k).at(0) += Ch.at(2)*Og.at(k).at(3); 
        coef_tot.at(k).at(0) += Ch.at(8)*Og.at(k).at(4);
        coef_tot.at(k).at(0) += (Ch.at(6)+Ch.at(7))*Og.at(k).at(5);
      }
      for(int k=0; k<2;k++){ // quark pieces
        coef_tot.at(k).insert({1,Ch.at(0)*Oq.at(k).at(0)+Ch.at(3)*Oq.at(k).at(1)+Ch.at(6)*Oq.at(k).at(5)});
        coef_tot.at(k).insert({2,(Ch.at(0)+Ch.at(1))*Oq.at(k).at(0)+Ch.at(2)*Oq.at(k).at(3)});
        coef_tot.at(k).insert({3,Ch.at(1)*Oq.at(k).at(0)+Ch.at(4)*Oq.at(k).at(1)+Ch.at(7)*Oq.at(k).at(5)});
        coef_tot.at(k).insert({4,(Ch.at(3)+Ch.at(4))*Oq.at(k).at(1)+Ch.at(5)*Oq.at(k).at(2)});
        coef_tot.at(k).insert({5,Ch.at(2)*Oq.at(k).at(3)+Ch.at(5)*Oq.at(k).at(2)+Ch.at(8)*Oq.at(k).at(4)});
        coef_tot.at(k).insert({6,(Ch.at(6)+Ch.at(7))*Oq.at(k).at(5)+Ch.at(8)*Oq.at(k).at(4)});
      }
      for(int i=0;i<=6;i++){
        coef_tot.at(2).insert({i,Zero}); 
      }
      //light
      for(int k=0; k<2; k++){ //gluon pieces
        coef_light.at(k).insert({0,(Ch.at(0)+Ch.at(1))*Og.at(k).at(0)}); 
        coef_light.at(k).at(0) += (Ch.at(3)+Ch.at(4))*Og.at(k).at(1);
        coef_light.at(k).at(0) += Ch.at(2)*Og.at(k).at(3); 
        coef_light.at(k).at(0) += (Ch.at(6)+Ch.at(7))*Og.at(k).at(5);
      }
      for(int k=0; k<2;k++){ //quark pieces
        coef_light.at(k).insert({1,Ch.at(0)*Oq.at(k).at(0)+Ch.at(3)*Oq.at(k).at(1)+Ch.at(6)*Oq.at(k).at(5)});
        coef_light.at(k).insert({2,(Ch.at(0)+Ch.at(1))*Oq.at(k).at(0)+Ch.at(2)*Oq.at(k).at(3)});
        coef_light.at(k).insert({3,Ch.at(1)*Oq.at(k).at(0)+Ch.at(4)*Oq.at(k).at(1)+Ch.at(7)*Oq.at(k).at(5)});
        coef_light.at(k).insert({4,(Ch.at(3)+Ch.at(4))*Oq.at(k).at(1)});
        coef_light.at(k).insert({5,Ch.at(2)*Oq.at(k).at(3)});
        coef_light.at(k).insert({6,(Ch.at(6)+Ch.at(7))*Oq.at(k).at(5)});
      }
      for(int i=0;i<=6;i++){
        coef_light.at(2).insert({i,Zero}); 
      }
      //charm
      for(int k=0; k<2; k++){
        coef_charm.at(k).insert({0,(Ch.at(3)+Ch.at(4))*Og.at(k).at(1)});
        coef_charm.at(k).at(0) += Ch.at(5)*Og.at(k).at(2);
      }
      for(int k=0; k<2; k++){
        coef_charm.at(k).insert({1,Ch.at(3)*Oq.at(k).at(1)});
        coef_charm.at(k).insert({2,Zero});
        coef_charm.at(k).insert({3,Ch.at(4)*Oq.at(k).at(1)});
        coef_charm.at(k).insert({4,(Ch.at(4)+Ch.at(3))*Oq.at(k).at(1)+Ch.at(5)*Oq.at(k).at(2)});
        coef_charm.at(k).insert({5,Ch.at(5)*Oq.at(k).at(2)});
        coef_charm.at(k).insert({6,Zero});
      }
      for(int i=0;i<=6;i++){
        coef_charm.at(2).insert({i,Zero}); 
      }
      //bottom
      for(int k=0; k<2; k++){ // gluon pieces
        coef_bottom.at(k).insert({0,Ch.at(5)*Og.at(k).at(2)});
        coef_bottom.at(k).at(0) += Ch.at(2)*Og.at(k).at(3); 
        coef_bottom.at(k).at(0) += Ch.at(8)*Og.at(k).at(4);
      }
      for(int k=0; k<2;k++){ // quark pieces
        coef_bottom.at(k).insert({1,Zero});
        coef_bottom.at(k).insert({2,Ch.at(2)*Oq.at(k).at(3)});
        coef_bottom.at(k).insert({3,Zero});
        coef_bottom.at(k).insert({4,Ch.at(5)*Oq.at(k).at(2)});
        coef_bottom.at(k).insert({5,Ch.at(2)*Oq.at(k).at(3)+Ch.at(5)*Oq.at(k).at(2)+Ch.at(8)*Oq.at(k).at(4)});
        coef_bottom.at(k).insert({6,Ch.at(8)*Oq.at(k).at(4)});
      }
      for(int i=0;i<=6;i++){
        coef_bottom.at(2).insert({i,Zero}); 
      }
      //top
      for(int k=0; k<2; k++){ // gluon pieces
        coef_top.at(k).insert({0,Ch.at(8)*Og.at(k).at(4)});
        coef_top.at(k).at(0) += (Ch.at(6)+Ch.at(7))*Og.at(k).at(5);
      }
      for(int k=0; k<2;k++){ // quark pieces
        coef_top.at(k).insert({1,Ch.at(6)*Oq.at(k).at(5)});
        coef_top.at(k).insert({2,Zero});
        coef_top.at(k).insert({3,Ch.at(7)*Oq.at(k).at(5)});
        coef_top.at(k).insert({4,Zero});
        coef_top.at(k).insert({5,Ch.at(8)*Oq.at(k).at(4)});
        coef_top.at(k).insert({6,(Ch.at(6)+Ch.at(7))*Oq.at(k).at(5)+Ch.at(8)*Oq.at(k).at(4)});
      }
      for(int i=0;i<=6;i++){
        coef_top.at(2).insert({i,Zero}); 
      }
      
      DISCCBasis_ACOT DISbasis_all(Ch,Zero);
      auto C_tot = DISbasis_all.get_operators_plus(coef_tot);
      FObj.ConvBasis.insert({0, DISbasis_all});
      FObj.C0.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C_tot.at(0)}});
      FObj.C1.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C_tot.at(1)}});
      FObj.C2.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C_tot.at(2)}});

      DISCCBasis_ACOT DISbasis_light(Ch,Zero);
      auto C_light = DISbasis_light.get_operators_plus(coef_light);
      FObj.ConvBasis.insert({1, DISbasis_light});
      FObj.C0.insert({1, Set<Operator>{FObj.ConvBasis.at(1), C_light.at(0)}});
      FObj.C1.insert({1, Set<Operator>{FObj.ConvBasis.at(1), C_light.at(1)}});
      FObj.C2.insert({1, Set<Operator>{FObj.ConvBasis.at(1), C_light.at(2)}});

      DISCCBasis_ACOT DISbasis_charm(Ch,Zero);
      auto C_charm = DISbasis_charm.get_operators_plus(coef_charm);
      FObj.ConvBasis.insert({2, DISbasis_charm});
      FObj.C0.insert({2, Set<Operator>{FObj.ConvBasis.at(2), C_charm.at(0)}});
      FObj.C1.insert({2, Set<Operator>{FObj.ConvBasis.at(2), C_charm.at(1)}});
      FObj.C2.insert({2, Set<Operator>{FObj.ConvBasis.at(2), C_charm.at(2)}});

      DISCCBasis_ACOT DISbasis_bottom(Ch,Zero);
      auto C_bottom = DISbasis_bottom.get_operators_plus(coef_bottom);
      FObj.ConvBasis.insert({3, DISbasis_bottom});
      FObj.C0.insert({3, Set<Operator>{FObj.ConvBasis.at(3), C_bottom.at(0)}});
      FObj.C1.insert({3, Set<Operator>{FObj.ConvBasis.at(3), C_bottom.at(1)}});
      FObj.C2.insert({3, Set<Operator>{FObj.ConvBasis.at(3), C_bottom.at(2)}});

      DISCCBasis_ACOT DISbasis_top(Ch,Zero);
      auto C_top = DISbasis_top.get_operators_plus(coef_top);
      FObj.ConvBasis.insert({4, DISbasis_top});
      FObj.C0.insert({4, Set<Operator>{FObj.ConvBasis.at(4), C_top.at(0)}});
      FObj.C1.insert({4, Set<Operator>{FObj.ConvBasis.at(4), C_top.at(1)}});
      FObj.C2.insert({4, Set<Operator>{FObj.ConvBasis.at(4), C_top.at(2)}});

      
      return FObj;
    };

    t.stop();
    return FLObj;
  }

  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeFLCCMinusObjectsSACOT(Grid                const& g,
                                                                                                                     std::vector<double> const& Masses,
                                                                                                                     double              const& IntEps,
                                                                                                                     int                 const& nQ,
                                                                                                                     double              const& Qmin,
                                                                                                                     double              const& Qmax,
                                                                                                                     int                 const& intdeg)
  {  
    Timer t;

    const int num_flavours = Masses.size();
    if (num_flavours != 6)
      throw std::runtime_error(error("InitializeFLCCMinusObjectsSACOT", "The vector of masses has to contain exactly 6 ordered masses."));
    
    report("Initializing StructureFunctionObjects for FL CC minus simplified ACOT up to NLO\n");

    const Operator Zero {g, Null{}, IntEps};

    const double ml = Masses[0] == 0 ? 0.01 : Masses[0];
    const double mc = Masses[3] == 0 ? 0.01 : Masses[3];
    const double mb = Masses[4] == 0 ? 0.01 : Masses[4];
    const double mt = Masses[5] == 0 ? 0.01 : Masses[5];
    const double mll = ml+ml;
    const double mlc = ml+mc;
    const double mcb = mc+mb;
    const double mlb = ml+mb;
    const double mlt = mt+ml;
    const double mtb = mt+mb;
    const std::vector<double> m({ml,ml,ml,mc,mb,mt});

    // calc gridding range
    const double Qgridmin = 0.95*Qmin;
    const double Qgridmax = 1.05*Qmax;
    const double sximin = Qgridmin/(mt+mb);
    const double sximax = Qgridmax/ml;
    const double lambda = 0.99*sximin;

    const auto fOL1ns =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1 / xi );
      const Operator OL1nsNC{g,CmL1qNC_ACOT_chi{eta},IntEps};
      return OL1nsNC;
    };
    const TabulateObject<Operator> TabOL1ns_H{fOL1ns, nQ, sximin, sximax, intdeg, {}, lambda};

    // Vector of distributions to skip
    const std::vector<int> skip = {0, 1, 2, 3, 5, 7, 9, 11};

    const auto FLObj = [=,&g] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects{
      if(Q<Qgridmin || Q>Qgridmax){
        throw std::runtime_error(error("FLCCsimACOT minus", "Q out of range ["+std::to_string(Qmin)+","+std::to_string(Qmax)+"]. Q=" + std::to_string(Q)));
      }
      StructureFunctionObjects FObj;
      FObj.skip = skip;

      std::vector<Operator> OqLO( {Zero,Zero,Zero,Zero,Zero,Zero});
      std::vector<Operator> OqNLO({TabOL1ns_H.Evaluate(Q/mll),   //0
                                   TabOL1ns_H.Evaluate(Q/mlc),   //1
                                   TabOL1ns_H.Evaluate(Q/mcb),   //2
                                   TabOL1ns_H.Evaluate(Q/mlb),   //3
                                   TabOL1ns_H.Evaluate(Q/mtb),   //4
                                   TabOL1ns_H.Evaluate(Q/mlt)}); //5
      std::vector<std::vector<Operator>> Oq({OqLO,OqNLO});

      std::vector<std::map<int,Operator>> coef_tot(3);
      std::vector<std::map<int,Operator>> coef_light(3);
      std::vector<std::map<int,Operator>> coef_charm(3);
      std::vector<std::map<int,Operator>> coef_bottom(3);
      std::vector<std::map<int,Operator>> coef_top(3);

      //tot
      for(int k=0; k<3; k++){//gluon pieces
        coef_tot.at(k).insert({0,Zero});
      }
      for(int k=0; k<2;k++){// quark pieces
        coef_tot.at(k).insert({1,Ch.at(0)*Oq.at(k).at(0)+Ch.at(3)*Oq.at(k).at(1)+Ch.at(6)*Oq.at(k).at(5)});
        coef_tot.at(k).insert({2,-(Ch.at(0)+Ch.at(1))*Oq.at(k).at(0)-Ch.at(2)*Oq.at(k).at(3)});
        coef_tot.at(k).insert({3,Ch.at(1)*Oq.at(k).at(0)+Ch.at(4)*Oq.at(k).at(1)+Ch.at(7)*Oq.at(k).at(5)});
        coef_tot.at(k).insert({4,-(Ch.at(4)+Ch.at(3))*Oq.at(k).at(1)-Ch.at(5)*Oq.at(k).at(2)});
        coef_tot.at(k).insert({5,Ch.at(2)*Oq.at(k).at(3)+Ch.at(5)*Oq.at(k).at(2)+Ch.at(8)*Oq.at(k).at(4)});
        coef_tot.at(k).insert({6,-(Ch.at(6)+Ch.at(7))*Oq.at(k).at(5)-Ch.at(8)*Oq.at(k).at(4)});
      }
      for(int i=1;i<=6;i++){
        coef_tot.at(2).insert({i,Zero}); 
      }
      //light
      for(int k=0; k<3; k++){ //gluon pieces
        coef_light.at(k).insert({0,Zero});
      }
      for(int k=0; k<2;k++){ //quark pieces
        coef_light.at(k).insert({1,Ch.at(0)*Oq.at(k).at(0)+Ch.at(3)*Oq.at(k).at(1)+Ch.at(6)*Oq.at(k).at(5)});
        coef_light.at(k).insert({2,-(Ch.at(0)+Ch.at(1))*Oq.at(k).at(0)-Ch.at(2)*Oq.at(k).at(3)});
        coef_light.at(k).insert({3,Ch.at(1)*Oq.at(k).at(0)+Ch.at(4)*Oq.at(k).at(1)+Ch.at(7)*Oq.at(k).at(5)});
        coef_light.at(k).insert({4,-(Ch.at(3)+Ch.at(4))*Oq.at(k).at(1)});
        coef_light.at(k).insert({5,Ch.at(2)*Oq.at(k).at(3)});
        coef_light.at(k).insert({6,-(Ch.at(6)+Ch.at(7))*Oq.at(k).at(5)});
      }
      for(int i=1;i<=6;i++){
        coef_light.at(2).insert({i,Zero}); 
      }
      //charm
      for(int k=0; k<3; k++){ //gluon pieces
        coef_charm.at(k).insert({0,Zero});
      }
      //down,up
      for(int k=0; k<2; k++){
        coef_charm.at(k).insert({1,Ch.at(3)*Oq.at(k).at(1)});
        coef_charm.at(k).insert({2,Zero});
        coef_charm.at(k).insert({3,Ch.at(4)*Oq.at(k).at(1)});
        coef_charm.at(k).insert({4,-(Ch.at(4)+Ch.at(3))*Oq.at(k).at(1)-Ch.at(5)*Oq.at(k).at(2)});
        coef_charm.at(k).insert({5,Ch.at(5)*Oq.at(k).at(2)});
        coef_charm.at(k).insert({6,Zero});
      }
      for(int i=1;i<=6;i++){
        coef_charm.at(2).insert({i,Zero}); 
      }
      //bottom
      for(int k=0; k<3; k++){ // gluon pieces
        coef_bottom.at(k).insert({0,Zero});
      }
      for(int k=0; k<2;k++){ // quark pieces
        coef_bottom.at(k).insert({1,Zero});
        coef_bottom.at(k).insert({2,-Ch.at(2)*Oq.at(k).at(3)});
        coef_bottom.at(k).insert({3,Zero});
        coef_bottom.at(k).insert({4,-Ch.at(5)*Oq.at(k).at(2)});
        coef_bottom.at(k).insert({5,Ch.at(2)*Oq.at(k).at(3)+Ch.at(5)*Oq.at(k).at(2)+Ch.at(8)*Oq.at(k).at(4)});
        coef_bottom.at(k).insert({6,-Ch.at(8)*Oq.at(k).at(4)});
      }
      for(int i=1;i<=6;i++){
        coef_bottom.at(2).insert({i,Zero}); 
      }
      //top
      for(int k=0; k<3; k++){ // gluon pieces
        coef_top.at(k).insert({0,Zero});
      }
      for(int k=0; k<2;k++){ // quark pieces
        coef_top.at(k).insert({1,Ch.at(6)*Oq.at(k).at(5)});
        coef_top.at(k).insert({2,Zero});
        coef_top.at(k).insert({3,Ch.at(7)*Oq.at(k).at(5)});
        coef_top.at(k).insert({4,Zero});
        coef_top.at(k).insert({5,Ch.at(8)*Oq.at(k).at(4)});
        coef_top.at(k).insert({6,-(Ch.at(6)+Ch.at(7))*Oq.at(k).at(5)-Ch.at(8)*Oq.at(k).at(4)});
      }
      for(int i=1;i<=6;i++){
        coef_top.at(2).insert({i,Zero}); 
      }

      
      DISCCBasis_ACOT DISbasis_all(Ch,Zero);
      auto C_tot = DISbasis_all.get_operators_minus(coef_tot);
      FObj.ConvBasis.insert({0, DISbasis_all});
      FObj.C0.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C_tot.at(0)}});
      FObj.C1.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C_tot.at(1)}});
      FObj.C2.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C_tot.at(2)}});

      DISCCBasis_ACOT DISbasis_light(Ch,Zero);
      auto C_light = DISbasis_light.get_operators_minus(coef_light);
      FObj.ConvBasis.insert({1, DISbasis_light});
      FObj.C0.insert({1, Set<Operator>{FObj.ConvBasis.at(1), C_light.at(0)}});
      FObj.C1.insert({1, Set<Operator>{FObj.ConvBasis.at(1), C_light.at(1)}});
      FObj.C2.insert({1, Set<Operator>{FObj.ConvBasis.at(1), C_light.at(2)}});

      DISCCBasis_ACOT DISbasis_charm(Ch,Zero);
      auto C_charm = DISbasis_charm.get_operators_minus(coef_charm);
      FObj.ConvBasis.insert({2, DISbasis_charm});
      FObj.C0.insert({2, Set<Operator>{FObj.ConvBasis.at(2), C_charm.at(0)}});
      FObj.C1.insert({2, Set<Operator>{FObj.ConvBasis.at(2), C_charm.at(1)}});
      FObj.C2.insert({2, Set<Operator>{FObj.ConvBasis.at(2), C_charm.at(2)}});

      DISCCBasis_ACOT DISbasis_bottom(Ch,Zero);
      auto C_bottom = DISbasis_bottom.get_operators_minus(coef_bottom);
      FObj.ConvBasis.insert({3, DISbasis_bottom});
      FObj.C0.insert({3, Set<Operator>{FObj.ConvBasis.at(3), C_bottom.at(0)}});
      FObj.C1.insert({3, Set<Operator>{FObj.ConvBasis.at(3), C_bottom.at(1)}});
      FObj.C2.insert({3, Set<Operator>{FObj.ConvBasis.at(3), C_bottom.at(2)}});

      DISCCBasis_ACOT DISbasis_top(Ch,Zero);
      auto C_top = DISbasis_top.get_operators_minus(coef_top);
      FObj.ConvBasis.insert({4, DISbasis_top});
      FObj.C0.insert({4, Set<Operator>{FObj.ConvBasis.at(4), C_top.at(0)}});
      FObj.C1.insert({4, Set<Operator>{FObj.ConvBasis.at(4), C_top.at(1)}});
      FObj.C2.insert({4, Set<Operator>{FObj.ConvBasis.at(4), C_top.at(2)}});
      
      return FObj;
    };

    t.stop();
    return FLObj;
  }

  std::function<StructureFunctionObjects(double const &, std::vector<double> const &)> InitializeF3CCPlusObjectsSACOT(Grid const &g,
                                                                                                                      std::vector<double> const &Masses,
                                                                                                                      double const &IntEps,
                                                                                                                      int const& nQ,
                                                                                                                      double const& Qmin,
                                                                                                                      double const& Qmax,
                                                                                                                      int const& intdeg)
  {
    Timer t;

    const int num_flavours = Masses.size();
    if (num_flavours != 6)
      throw std::runtime_error(error("InitializeF3CCPlusObjectsSACOT", "The vector of masses has to contain exactly 6 ordered masses."));
    
    report("Initializing StructureFunctionObjects for F3 CC plus simplified ACOT up to NLO\n");

    const Operator Zero {g, Null{}, IntEps};

    const double ml = Masses[0] == 0 ? 0.01 : Masses[0];
    const double mc = Masses[3] == 0 ? 0.01 : Masses[3];
    const double mb = Masses[4] == 0 ? 0.01 : Masses[4];
    const double mt = Masses[5] == 0 ? 0.01 : Masses[5];
    const double mll = ml+ml;
    const double ml2 = ml*ml;
    const double mc2 = mc*mc;
    const double mlc = ml+mc;
    const double mcb = mc+mb;
    const double mb2 = mb*mb;
    const double mlb = ml+mb;
    const double mt2 = mt*mt;
    const double mlt = mt+ml;
    const double mtb = mt+mb;
    const std::vector<double> m({ml,ml,ml,mc,mb,mt});

    // calc gridding range
    const double Qgridmin = 0.95*Qmin;
    const double Qgridmax = 1.05*Qmax;
    const double sximin = Qgridmin/(mt+mb);
    const double sximax = Qgridmax/ml;
    const double lambda = 0.99*sximin;

    const auto fO30ns =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1 / xi );
      const Operator O30nsNC{g,Cm20qNC_ACOT_chi{eta},IntEps};
      return O30nsNC;
    };
    const TabulateObject<Operator> TabO30ns_H{fO30ns, nQ, sximin, sximax, intdeg, {}, lambda};

    const auto fO31ns =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1 / xi );
      const Operator O31nsNC{g,Cm31qNC_ACOT_chi{eta},IntEps};
      return O31nsNC;
    };
    const TabulateObject<Operator> TabO31ns_H{fO31ns, nQ, sximin, sximax, intdeg, {}, lambda};

    const auto fO31g_light_charm =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1./ xi );
      const Operator O31CCg_gen_mass{g,Cm31gCC_general_mass{eta,xi,mc,ml},IntEps};
      return O31CCg_gen_mass;
    };
    const TabulateObject<Operator> TabO31g_lc{fO31g_light_charm, nQ, Qgridmin/mlc, Qgridmax/mlc, intdeg, {}, 9e-6};

    const auto fO31g_charm_bottom =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1./ xi );
      const Operator O31CCg_gen_mass{g,Cm31gCC_general_mass{eta,xi,mc,mb},IntEps};
      return O31CCg_gen_mass;
    };
    const TabulateObject<Operator> TabO31g_cb{fO31g_charm_bottom, nQ, Qgridmin/mcb, Qgridmax/mcb, intdeg, {}, lambda};

    const auto fO31g_light_bottom =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1./ xi );
      const Operator O31CCg_gen_mass{g,Cm31gCC_general_mass{eta,xi,ml,mb},IntEps};
      return O31CCg_gen_mass;
    };
    const TabulateObject<Operator> TabO31g_lb{fO31g_light_bottom, nQ, Qgridmin/mlb, Qgridmax/mlb, intdeg, {}, lambda};

    const auto fO31g_top_bottom =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1./ xi );
      const Operator O31CCg_gen_mass{g,Cm31gCC_general_mass{eta,xi,mt,mb},IntEps};
      return O31CCg_gen_mass;
    };
    const TabulateObject<Operator> TabO31g_tb{fO31g_top_bottom, nQ, Qgridmin/mtb, Qgridmax/mtb, intdeg, {}, lambda};

    const auto fO31g_light_top =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1./ xi );
      const Operator O31CCg_gen_mass{g,Cm31gCC_general_mass{eta,xi,mt,ml},IntEps};
      return O31CCg_gen_mass;
    };
    const TabulateObject<Operator> TabO31g_lt{fO31g_light_top, nQ, Qgridmin/mlt, Qgridmax/mlt, intdeg, {}, lambda};

    const auto fO31g_sub =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1./ xi );
      const Operator Om31CCg_sub{g,Cm21gNC_sub_ACOT_chi{eta},IntEps};
      return 0.5*Om31CCg_sub;
    };
    const TabulateObject<Operator> TabO31g_sub{fO31g_sub, nQ, sximin, sximax, intdeg, {1}, lambda};

    // Vector of distributions to skip
    const std::vector<int> skip = {2, 4, 6, 8, 10, 12};

    const auto F3Obj = [=,&g] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects{
      if(Q<Qgridmin || Q>Qgridmax){
        throw std::runtime_error(error("F3CCsimACOT plus", "Q out of range ["+std::to_string(Qmin)+","+std::to_string(Qmax)+"]. Q=" + std::to_string(Q)));
      }
      StructureFunctionObjects FObj;
      FObj.skip = skip;
      const double Q2 = Q*Q;
      const double log_ml = Q>=ml ? log(Q2/ml2) : 0;
      const double log_mc = Q>=mc ? log(Q2/mc2) : 0;
      const double log_mb = Q>=mb ? log(Q2/mb2) : 0;
      const double log_mt = Q>=mt ? log(Q2/mt2) : 0;  

      std::vector<Operator> OqLO({  TabO30ns_H.Evaluate(Q/mll),   //0
                                    TabO30ns_H.Evaluate(Q/mlc),   //1
                                    TabO30ns_H.Evaluate(Q/mcb),   //2
                                    TabO30ns_H.Evaluate(Q/mlb),   //3
                                    TabO30ns_H.Evaluate(Q/mtb),   //4
                                    TabO30ns_H.Evaluate(Q/mlt)}); //5
      std::vector<Operator> OqNLO({ TabO31ns_H.Evaluate(Q/mll),   //0
                                    TabO31ns_H.Evaluate(Q/mlc),   //1
                                    TabO31ns_H.Evaluate(Q/mcb),   //2
                                    TabO31ns_H.Evaluate(Q/mlb),   //3
                                    TabO31ns_H.Evaluate(Q/mtb),   //4
                                    TabO31ns_H.Evaluate(Q/mlt)}); //5
      std::vector<std::vector<Operator>> Oq({OqLO,OqNLO});

      std::vector<Operator> OgLO( { Zero,Zero,Zero,Zero,Zero,Zero});
      std::vector<Operator> OgNLO({ Zero, //light light is zero due to equal masses
                                    TabO31g_lc.Evaluate(Q/mlc) - (log_ml-log_mc)*TabO31g_sub.Evaluate(Q/mlc),
                                    TabO31g_cb.Evaluate(Q/mcb) - (log_mb-log_mc)*TabO31g_sub.Evaluate(Q/mcb),
                                    TabO31g_lb.Evaluate(Q/mlb) - (log_mb-log_ml)*TabO31g_sub.Evaluate(Q/mlb),
                                    TabO31g_tb.Evaluate(Q/mtb) - (log_mb-log_mt)*TabO31g_sub.Evaluate(Q/mtb),
                                    TabO31g_lt.Evaluate(Q/mlt) - (log_ml-log_mt)*TabO31g_sub.Evaluate(Q/mlt)});
      std::vector<std::vector<Operator>> Og({OgLO,OgNLO});

      std::vector<std::map<int,Operator>> coef_tot(3);
      std::vector<std::map<int,Operator>> coef_light(3);
      std::vector<std::map<int,Operator>> coef_charm(3);
      std::vector<std::map<int,Operator>> coef_bottom(3);
      std::vector<std::map<int,Operator>> coef_top(3);

      //total
      for(int k=0; k<2; k++){ // gluon pieces
        coef_tot.at(k).insert({0,(Ch.at(0)+Ch.at(1))*Og.at(k).at(0)}); 
        coef_tot.at(k).at(0) += (Ch.at(3)+Ch.at(4))*Og.at(k).at(1);
        coef_tot.at(k).at(0) += Ch.at(5)*Og.at(k).at(2);
        coef_tot.at(k).at(0) += Ch.at(2)*Og.at(k).at(3); 
        coef_tot.at(k).at(0) += Ch.at(8)*Og.at(k).at(4);
        coef_tot.at(k).at(0) += (Ch.at(6)+Ch.at(7))*Og.at(k).at(5);
      }
      for(int k=0; k<2;k++){ // quark pieces
        coef_tot.at(k).insert({1,Ch.at(0)*Oq.at(k).at(0)+Ch.at(3)*Oq.at(k).at(1)+Ch.at(6)*Oq.at(k).at(5)});
        coef_tot.at(k).insert({2,-(Ch.at(0)+Ch.at(1))*Oq.at(k).at(0)-Ch.at(2)*Oq.at(k).at(3)});
        coef_tot.at(k).insert({3,Ch.at(1)*Oq.at(k).at(0)+Ch.at(4)*Oq.at(k).at(1)+Ch.at(7)*Oq.at(k).at(5)});
        coef_tot.at(k).insert({4,-(Ch.at(3)+Ch.at(4))*Oq.at(k).at(1)-Ch.at(5)*Oq.at(k).at(2)});
        coef_tot.at(k).insert({5,Ch.at(2)*Oq.at(k).at(3)+Ch.at(5)*Oq.at(k).at(2)+Ch.at(8)*Oq.at(k).at(4)});
        coef_tot.at(k).insert({6,-(Ch.at(6)+Ch.at(7))*Oq.at(k).at(5)-Ch.at(8)*Oq.at(k).at(4)});
      }
      for(int i=0;i<=6;i++){
        coef_tot.at(2).insert({i,Zero});
      }
      //light
      for(int k=0; k<2; k++){ //gluon pieces
        coef_light.at(k).insert({0,(Ch.at(0)+Ch.at(1))*Og.at(k).at(0)}); 
        coef_light.at(k).at(0) += (Ch.at(3)+Ch.at(4))*Og.at(k).at(1);
        coef_light.at(k).at(0) += Ch.at(2)*Og.at(k).at(3); 
        coef_light.at(k).at(0) += (Ch.at(6)+Ch.at(7))*Og.at(k).at(5);
      }
      for(int k=0; k<2;k++){ //quark pieces
        coef_light.at(k).insert({1,Ch.at(0)*Oq.at(k).at(0)+Ch.at(3)*Oq.at(k).at(1)+Ch.at(6)*Oq.at(k).at(5)});
        coef_light.at(k).insert({2,-(Ch.at(0)+Ch.at(1))*Oq.at(k).at(0)-Ch.at(2)*Oq.at(k).at(3)});
        coef_light.at(k).insert({3,Ch.at(1)*Oq.at(k).at(0)+Ch.at(4)*Oq.at(k).at(1)+Ch.at(7)*Oq.at(k).at(5)});
        coef_light.at(k).insert({4,-(Ch.at(3)+Ch.at(4))*Oq.at(k).at(1)});
        coef_light.at(k).insert({5,Ch.at(2)*Oq.at(k).at(3)});
        coef_light.at(k).insert({6,-(Ch.at(6)+Ch.at(7))*Oq.at(k).at(5)});
      }
      for(int i=0;i<=6;i++){
        coef_light.at(2).insert({i,Zero}); 
      }
      //charm
      for(int k=0; k<2; k++){
        coef_charm.at(k).insert({0,(Ch.at(3)+Ch.at(4))*Og.at(k).at(1)});
        coef_charm.at(k).at(0) += Ch.at(5)*Og.at(k).at(2);
      }
      for(int k=0; k<2; k++){
        coef_charm.at(k).insert({1,Ch.at(3)*Oq.at(k).at(1)});
        coef_charm.at(k).insert({2,Zero});
        coef_charm.at(k).insert({3,Ch.at(4)*Oq.at(k).at(1)});
        coef_charm.at(k).insert({4,-(Ch.at(4)+Ch.at(3))*Oq.at(k).at(1)-Ch.at(5)*Oq.at(k).at(2)});
        coef_charm.at(k).insert({5,Ch.at(5)*Oq.at(k).at(2)});
        coef_charm.at(k).insert({6,Zero});
      }
      for(int i=0;i<=6;i++){
        coef_charm.at(2).insert({i,Zero}); 
      }
      //bottom
      for(int k=0; k<2; k++){ // gluon pieces
        coef_bottom.at(k).insert({0,Ch.at(5)*Og.at(k).at(2)});
        coef_bottom.at(k).at(0) += Ch.at(2)*Og.at(k).at(3); 
        coef_bottom.at(k).at(0) += Ch.at(8)*Og.at(k).at(4);
      }
      for(int k=0; k<2;k++){ // quark pieces
        coef_bottom.at(k).insert({1,Zero});
        coef_bottom.at(k).insert({2,-Ch.at(2)*Oq.at(k).at(3)});
        coef_bottom.at(k).insert({3,Zero});
        coef_bottom.at(k).insert({4,-Ch.at(5)*Oq.at(k).at(2)});
        coef_bottom.at(k).insert({5,Ch.at(2)*Oq.at(k).at(3)+Ch.at(5)*Oq.at(k).at(2)+Ch.at(8)*Oq.at(k).at(4)});
        coef_bottom.at(k).insert({6,-Ch.at(8)*Oq.at(k).at(4)});
      }
      for(int i=0;i<=6;i++){
        coef_bottom.at(2).insert({i,Zero}); 
      }
      //top
      for(int k=0; k<2; k++){ // gluon pieces
        coef_top.at(k).insert({0,Ch.at(8)*Og.at(k).at(4)});
        coef_top.at(k).at(0) += (Ch.at(6)+Ch.at(7))*Og.at(k).at(5);
      }
      for(int k=0; k<2;k++){ // quark pieces
        coef_top.at(k).insert({1,Ch.at(6)*Oq.at(k).at(5)});
        coef_top.at(k).insert({2,Zero});
        coef_top.at(k).insert({3,Ch.at(7)*Oq.at(k).at(5)});
        coef_top.at(k).insert({4,Zero});
        coef_top.at(k).insert({5,Ch.at(8)*Oq.at(k).at(4)});
        coef_top.at(k).insert({6,-(Ch.at(6)+Ch.at(7))*Oq.at(k).at(5)-Ch.at(8)*Oq.at(k).at(4)});
      }
      for(int i=0;i<=6;i++){
        coef_top.at(2).insert({i,Zero}); 
      }

      DISCCBasis_ACOT DISbasis_all(Ch,Zero);
      auto C_tot = DISbasis_all.get_operators_plus(coef_tot);
      FObj.ConvBasis.insert({0, DISbasis_all});
      FObj.C0.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C_tot.at(0)}});
      FObj.C1.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C_tot.at(1)}});
      FObj.C2.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C_tot.at(2)}});

      DISCCBasis_ACOT DISbasis_light(Ch,Zero);
      auto C_light = DISbasis_light.get_operators_plus(coef_light);
      FObj.ConvBasis.insert({1, DISbasis_light});
      FObj.C0.insert({1, Set<Operator>{FObj.ConvBasis.at(1), C_light.at(0)}});
      FObj.C1.insert({1, Set<Operator>{FObj.ConvBasis.at(1), C_light.at(1)}});
      FObj.C2.insert({1, Set<Operator>{FObj.ConvBasis.at(1), C_light.at(2)}});

      DISCCBasis_ACOT DISbasis_charm(Ch,Zero);
      auto C_charm = DISbasis_charm.get_operators_plus(coef_charm);
      FObj.ConvBasis.insert({2, DISbasis_charm});
      FObj.C0.insert({2, Set<Operator>{FObj.ConvBasis.at(2), C_charm.at(0)}});
      FObj.C1.insert({2, Set<Operator>{FObj.ConvBasis.at(2), C_charm.at(1)}});
      FObj.C2.insert({2, Set<Operator>{FObj.ConvBasis.at(2), C_charm.at(2)}});

      DISCCBasis_ACOT DISbasis_bottom(Ch,Zero);
      auto C_bottom = DISbasis_bottom.get_operators_plus(coef_bottom);
      FObj.ConvBasis.insert({3, DISbasis_bottom});
      FObj.C0.insert({3, Set<Operator>{FObj.ConvBasis.at(3), C_bottom.at(0)}});
      FObj.C1.insert({3, Set<Operator>{FObj.ConvBasis.at(3), C_bottom.at(1)}});
      FObj.C2.insert({3, Set<Operator>{FObj.ConvBasis.at(3), C_bottom.at(2)}});

      DISCCBasis_ACOT DISbasis_top(Ch,Zero);
      auto C_top = DISbasis_top.get_operators_plus(coef_top);
      FObj.ConvBasis.insert({4, DISbasis_top});
      FObj.C0.insert({4, Set<Operator>{FObj.ConvBasis.at(4), C_top.at(0)}});
      FObj.C1.insert({4, Set<Operator>{FObj.ConvBasis.at(4), C_top.at(1)}});
      FObj.C2.insert({4, Set<Operator>{FObj.ConvBasis.at(4), C_top.at(2)}});
      
      return FObj;
    };

    t.stop();
    return F3Obj;
  }

  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF3CCMinusObjectsSACOT(Grid                const& g,
                                                                                                                     std::vector<double> const& Masses,
                                                                                                                     double              const& IntEps,
                                                                                                                     int                 const& nQ,
                                                                                                                     double              const& Qmin,
                                                                                                                     double              const& Qmax,
                                                                                                                     int                 const& intdeg)
  {  
    Timer t;

    const int num_flavours = Masses.size();
    if (num_flavours != 6)
      throw std::runtime_error(error("InitializeF3CCMinusObjectsSACOT", "The vector of masses has to contain exactly 6 ordered masses."));
    
    report("Initializing StructureFunctionObjects for F3 CC minus simplified ACOT up to NLO\n");

    const Operator Zero {g, Null{}, IntEps};

    const double ml = Masses[0] == 0 ? 0.01 : Masses[0];
    const double mc = Masses[3] == 0 ? 0.01 : Masses[3];
    const double mb = Masses[4] == 0 ? 0.01 : Masses[4];
    const double mt = Masses[5] == 0 ? 0.01 : Masses[5];
    const double mll = ml+ml;
    const double mlc = ml+mc;
    const double mcb = mc+mb;
    const double mlb = ml+mb;
    const double mlt = mt+ml;
    const double mtb = mt+mb;
    const std::vector<double> m({ml,ml,ml,mc,mb,mt});

    // calc gridding range
    const double Qgridmin = 0.95*Qmin;
    const double Qgridmax = 1.05*Qmax;
    const double sximin = Qgridmin/(mt+mb);
    const double sximax = Qgridmax/ml;
    const double lambda = 0.99*sximin;

    const auto fO30ns =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1 / xi );
      const Operator O30nsNC{g,Cm20qNC_ACOT_chi{eta},IntEps};
      return O30nsNC;
    };
    const TabulateObject<Operator> TabO30ns_H{fO30ns, nQ, sximin, sximax, intdeg, {}, lambda};

    const auto fO31ns =  [=,&g] (double const& sxi) -> Operator{
      const double xi = sxi*sxi;
      const double eta = 1/( 1 + 1 / xi );
      const Operator O31nsNC{g,Cm31qNC_ACOT_chi{eta},IntEps};
      return O31nsNC;
    };
    const TabulateObject<Operator> TabO31ns_H{fO31ns, nQ, sximin, sximax, intdeg, {}, lambda};

    // Vector of distributions to skip
    const std::vector<int> skip = {0, 1, 3, 5, 7, 9, 11};

    const auto F3Obj = [=,&g] (double const& Q, std::vector<double> const& Ch) -> StructureFunctionObjects{
      if(Q<Qgridmin || Q>Qgridmax){
        throw std::runtime_error(error("F3CCsimACOT NLO minus", "Q out of range ["+std::to_string(Qmin)+","+std::to_string(Qmax)+"]. Q=" + std::to_string(Q)));
      }
      StructureFunctionObjects FObj;
      FObj.skip = skip;
      const std::vector<double> m({ml,ml,ml,mc,mb,mt});

      std::vector<Operator> OqLO({  TabO30ns_H.Evaluate(Q/mll),   //0
                                    TabO30ns_H.Evaluate(Q/mlc),   //1
                                    TabO30ns_H.Evaluate(Q/mcb),   //2
                                    TabO30ns_H.Evaluate(Q/mlb),   //3
                                    TabO30ns_H.Evaluate(Q/mtb),   //4
                                    TabO30ns_H.Evaluate(Q/mlt)}); //5
      std::vector<Operator> OqNLO({ TabO31ns_H.Evaluate(Q/mll),   //0
                                    TabO31ns_H.Evaluate(Q/mlc),   //1
                                    TabO31ns_H.Evaluate(Q/mcb),   //2
                                    TabO31ns_H.Evaluate(Q/mlb),   //3
                                    TabO31ns_H.Evaluate(Q/mtb),   //4
                                    TabO31ns_H.Evaluate(Q/mlt)}); //5
      std::vector<std::vector<Operator>> Oq({OqLO,OqNLO});

      std::vector<std::map<int,Operator>> coef_tot(3);
      std::vector<std::map<int,Operator>> coef_light(3);
      std::vector<std::map<int,Operator>> coef_charm(3);
      std::vector<std::map<int,Operator>> coef_bottom(3);
      std::vector<std::map<int,Operator>> coef_top(3);

      //total
      for(int k=0; k<3; k++){ //gluon pieces
        coef_tot.at(k).insert({0,Zero});
      }
      for(int k=0; k<2;k++){// quark pieces
        coef_tot.at(k).insert({1,Ch.at(0)*Oq.at(k).at(0)+Ch.at(3)*Oq.at(k).at(1)+Ch.at(6)*Oq.at(k).at(5)});
        coef_tot.at(k).insert({2,(Ch.at(0)+Ch.at(1))*Oq.at(k).at(0)+Ch.at(2)*Oq.at(k).at(3)});
        coef_tot.at(k).insert({3,Ch.at(1)*Oq.at(k).at(0)+Ch.at(4)*Oq.at(k).at(1)+Ch.at(7)*Oq.at(k).at(5)});
        coef_tot.at(k).insert({4,(Ch.at(4)+Ch.at(3))*Oq.at(k).at(1)+Ch.at(5)*Oq.at(k).at(2)});
        coef_tot.at(k).insert({5,Ch.at(2)*Oq.at(k).at(3)+Ch.at(5)*Oq.at(k).at(2)+Ch.at(8)*Oq.at(k).at(4)});
        coef_tot.at(k).insert({6,(Ch.at(6)+Ch.at(7))*Oq.at(k).at(5)+Ch.at(8)*Oq.at(k).at(4)});
      }
      for(int i=1;i<=6;i++){
        coef_tot.at(2).insert({i,Zero}); 
      }
      //light
      for(int k=0; k<3; k++){ //gluon pieces
        coef_light.at(k).insert({0,Zero});
      }
      for(int k=0; k<2;k++){ //quark pieces
        coef_light.at(k).insert({1,Ch.at(0)*Oq.at(k).at(0)+Ch.at(3)*Oq.at(k).at(1)+Ch.at(6)*Oq.at(k).at(5)});
        coef_light.at(k).insert({2,(Ch.at(0)+Ch.at(1))*Oq.at(k).at(0)+Ch.at(2)*Oq.at(k).at(3)});
        coef_light.at(k).insert({3,Ch.at(1)*Oq.at(k).at(0)+Ch.at(4)*Oq.at(k).at(1)+Ch.at(7)*Oq.at(k).at(5)});
        coef_light.at(k).insert({4,(Ch.at(3)+Ch.at(4))*Oq.at(k).at(1)});
        coef_light.at(k).insert({5,Ch.at(2)*Oq.at(k).at(3)});
        coef_light.at(k).insert({6,(Ch.at(6)+Ch.at(7))*Oq.at(k).at(5)});
      }
      for(int i=1;i<=6;i++){
        coef_light.at(2).insert({i,Zero}); 
      }
      //charm
      for(int k=0; k<3; k++){ //gluon pieces
        coef_charm.at(k).insert({0,Zero});
      }
      //down,up
      for(int k=0; k<2; k++){
        coef_charm.at(k).insert({1,Ch.at(3)*Oq.at(k).at(1)});
        coef_charm.at(k).insert({2,Zero});
        coef_charm.at(k).insert({3,Ch.at(4)*Oq.at(k).at(1)});
        coef_charm.at(k).insert({4,(Ch.at(4)+Ch.at(3))*Oq.at(k).at(1)+Ch.at(5)*Oq.at(k).at(2)});
        coef_charm.at(k).insert({5,Ch.at(5)*Oq.at(k).at(2)});
        coef_charm.at(k).insert({6,Zero});
      }
      for(int i=1;i<=6;i++){
        coef_charm.at(2).insert({i,Zero}); 
      }
      //bottom
      for(int k=0; k<3; k++){ // gluon pieces
        coef_bottom.at(k).insert({0,Zero});
      }
      for(int k=0; k<2;k++){ // quark pieces
        coef_bottom.at(k).insert({1,Zero});
        coef_bottom.at(k).insert({2,Ch.at(2)*Oq.at(k).at(3)});
        coef_bottom.at(k).insert({3,Zero});
        coef_bottom.at(k).insert({4,Ch.at(5)*Oq.at(k).at(2)});
        coef_bottom.at(k).insert({5,Ch.at(2)*Oq.at(k).at(3)+Ch.at(5)*Oq.at(k).at(2)+Ch.at(8)*Oq.at(k).at(4)});
        coef_bottom.at(k).insert({6,Ch.at(8)*Oq.at(k).at(4)});
      }
      for(int i=1;i<=6;i++){
        coef_bottom.at(2).insert({i,Zero}); 
      }
      //top
      for(int k=0; k<3; k++){ // gluon pieces
        coef_top.at(k).insert({0,Zero});
      }
      for(int k=0; k<2;k++){ // quark pieces
        coef_top.at(k).insert({1,Ch.at(6)*Oq.at(k).at(5)});
        coef_top.at(k).insert({2,Zero});
        coef_top.at(k).insert({3,Ch.at(7)*Oq.at(k).at(5)});
        coef_top.at(k).insert({4,Zero});
        coef_top.at(k).insert({5,Ch.at(8)*Oq.at(k).at(4)});
        coef_top.at(k).insert({6,(Ch.at(6)+Ch.at(7))*Oq.at(k).at(5)+Ch.at(8)*Oq.at(k).at(4)});
      }
      for(int i=1;i<=6;i++){
        coef_top.at(2).insert({i,Zero}); 
      }

      DISCCBasis_ACOT DISbasis_all(Ch,Zero);
      auto C_tot = DISbasis_all.get_operators_minus(coef_tot);
      FObj.ConvBasis.insert({0, DISbasis_all});
      FObj.C0.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C_tot.at(0)}});
      FObj.C1.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C_tot.at(1)}});
      FObj.C2.insert({0, Set<Operator>{FObj.ConvBasis.at(0), C_tot.at(2)}});

      DISCCBasis_ACOT DISbasis_light(Ch,Zero);
      auto C_light = DISbasis_light.get_operators_minus(coef_light);
      FObj.ConvBasis.insert({1, DISbasis_light});
      FObj.C0.insert({1, Set<Operator>{FObj.ConvBasis.at(1), C_light.at(0)}});
      FObj.C1.insert({1, Set<Operator>{FObj.ConvBasis.at(1), C_light.at(1)}});
      FObj.C2.insert({1, Set<Operator>{FObj.ConvBasis.at(1), C_light.at(2)}});

      DISCCBasis_ACOT DISbasis_charm(Ch,Zero);
      auto C_charm = DISbasis_charm.get_operators_minus(coef_charm);
      FObj.ConvBasis.insert({2, DISbasis_charm});
      FObj.C0.insert({2, Set<Operator>{FObj.ConvBasis.at(2), C_charm.at(0)}});
      FObj.C1.insert({2, Set<Operator>{FObj.ConvBasis.at(2), C_charm.at(1)}});
      FObj.C2.insert({2, Set<Operator>{FObj.ConvBasis.at(2), C_charm.at(2)}});

      DISCCBasis_ACOT DISbasis_bottom(Ch,Zero);
      auto C_bottom = DISbasis_bottom.get_operators_minus(coef_bottom);
      FObj.ConvBasis.insert({3, DISbasis_bottom});
      FObj.C0.insert({3, Set<Operator>{FObj.ConvBasis.at(3), C_bottom.at(0)}});
      FObj.C1.insert({3, Set<Operator>{FObj.ConvBasis.at(3), C_bottom.at(1)}});
      FObj.C2.insert({3, Set<Operator>{FObj.ConvBasis.at(3), C_bottom.at(2)}});

      DISCCBasis_ACOT DISbasis_top(Ch,Zero);
      auto C_top = DISbasis_top.get_operators_minus(coef_top);
      FObj.ConvBasis.insert({4, DISbasis_top});
      FObj.C0.insert({4, Set<Operator>{FObj.ConvBasis.at(4), C_top.at(0)}});
      FObj.C1.insert({4, Set<Operator>{FObj.ConvBasis.at(4), C_top.at(1)}});
      FObj.C2.insert({4, Set<Operator>{FObj.ConvBasis.at(4), C_top.at(2)}});
      
      return FObj;
    };

    t.stop();
    return F3Obj;
  }

}