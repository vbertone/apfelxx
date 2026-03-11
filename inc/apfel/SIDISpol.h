//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/apfelxx.h"
#include "apfel/SIDIS.h"
#include "apfel/sidiscoefficientfunctionspol.h"

namespace apfel
{
  /**
   * @brief Structure that contains the precomputed SIDIS polarised
   * hard cross sections.
   */
  struct SidisPolObjects
  {
    DoubleObject<Operator> G10qq;
    DoubleObject<Operator> G11qq;
    DoubleObject<Operator> G11gq;
    DoubleObject<Operator> G11qg;
    std::map<int, DoubleOperator>         G12ns;
    std::map<int, DoubleOperator>         G12gq;
    std::map<int, DoubleOperator>         G12qg;
    std::map<std::string, DoubleOperator> G12;
  };

  // Expressions needed for the computation of the SIDIS cross sections
  // F2
  class lrqqpol: public Expression
  {
  public:
    lrqqpol(): Expression() {}
    double Regular(double const& x) const
    {
      const double expr = ( 1 + x * x ) * log(x) / ( 1 - x ) + 1 - x - ( 1 + x ) * log( 1 - x );
      return 2 * CF * expr;
    }
  };

  class srqqpol: public Expression
  {
  public:
    srqqpol(): Expression() {}
    double Regular(double const& x) const { return - 2 * CF * ( 1 + x ); }
  };

  class rlqqpol: public Expression
  {
  public:
    rlqqpol(): Expression() {}
    double Regular(double const& x) const
    {
      const double expr = - ( 1 + x * x ) * log(x) / ( 1 - x ) + 1 - x - ( 1 + x ) * log( 1 - x );
      return 2 * CF * expr;
    }
  };

  class rsqqpol: public Expression
  {
  public:
    rsqqpol(): Expression() {}
    double Regular(double const& x) const { return - 2 * CF * ( 1 + x ); }
  };

  class r11qqpol: public Expression
  {
  public:
    r11qqpol(): Expression() {}
    double Regular(double const& x) const { return x; }
  };

  class r12qqpol: public Expression
  {
  public:
    r12qqpol(): Expression() {}
    double Regular(double const&) const { return 1; }
  };

  class r21qqpol: public Expression
  {
  public:
    r21qqpol(): Expression() {}
    double Regular(double const&) const { return 1; }
  };

  class r22qqpol: public Expression
  {
  public:
    r22qqpol(): Expression() {}
    double Regular(double const& x) const { return x; }
  };

  class lrgqpol: public Expression
  {
  public:
    lrgqpol(): Expression() {}
    double Regular(double const& x) const
    {
      const double omx = ( 1 - x );
      const double expr = ( 1 + omx * omx ) * log( x * omx ) / x + x;
      return 2 * CF * expr;
    }
  };

  class srgqpol: public Expression
  {
  public:
    srgqpol(): Expression() {}
    double Regular(double const& x) const
    {
      const double omx = ( 1 - x );
      const double expr = ( 1 + omx * omx ) / x;
      return 2 * CF * expr;
    }
  };

  class r11gqpol: public Expression
  {
  public:
    r11gqpol(): Expression() {}
    double Regular(double const& x) const { return 1 + x; }
  };

  class r12gqpol: public Expression
  {
  public:
    r12gqpol(): Expression() {}
    double Regular(double const& x) const { return 2 - 1 / x; }
  };

  class r21gqpol: public Expression
  {
  public:
    r21gqpol(): Expression() {}
    double Regular(double const&) const { return 1; }
  };

  class r22gqpol: public Expression
  {
  public:
    r22gqpol(): Expression() {}
    double Regular(double const& x) const { return x; }
  };

  class rlqgpol: public Expression
  {
  public:
    rlqgpol(): Expression() {}
    double Regular(double const& x) const
    {
      const double omx = ( 1 - x );
      const double expr = ( x * x - omx * omx ) * log( omx / x ) + 2 * omx;
      return expr;
    }
  };

  class rsqgpol: public Expression
  {
  public:
    rsqgpol(): Expression() {}
    double Regular(double const& x) const
    {
      const double omx = ( 1 - x );
      const double expr = ( x * x - omx * omx );
      return expr;
    }
  };

  class r11qgpol: public Expression
  {
  public:
    r11qgpol(): Expression() {}
    double Regular(double const& x) const { return x * x - ( 1 - x ) * ( 1 - x ); }
  };

  class r12qgpol: public Expression
  {
  public:
    r12qgpol(): Expression() {}
    double Regular(double const& x) const { return 1 / x - 2; }
  };

  // Functions that fills in the SIDIS hard cross sections on two
  // different grids.
  SidisPolObjects InitializeSIDISpol(Grid const& gx, Grid const& gz, std::vector<double> const& Thresholds, std::vector<int> exclude = {}, double const& IntEps = 1e-3)
  {
    report("Initializing SIDIS longitudinally polarised hard cross sections... ");
    Timer t;

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the threshold
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
        nfi++;

    // Define object of the structure containing the DglapObjects.
    SidisPolObjects SidisObj;

    // ====================================================
    // Compute SIDIS partonic cross sections.
    // Expressions taken from Appendix C of hep-ph/9711387.
    // ====================================================
    // LO contribution
    const Operator odeltax{gx, delta{}};
    const Operator odeltaz{gz, delta{}};

    if (std::find(exclude.begin(), exclude.end(), 0) == exclude.end()) SidisObj.G10qq.AddTerm({1, odeltax, odeltaz});

    // NLO contributions
    // G1
    const Operator oD0x{gx, DDn{0}};
    const Operator oD1x{gx, DDn{1}};
    const Operator oD0z{gz, DDn{0}};
    const Operator oD1z{gz, DDn{1}};

    const double LLqq = - 16 * CF;
    const double LSqq = 4 * CF;
    const double SLqq = 4 * CF;
    const double SSqq = 4 * CF;
    const double K1qq = 4 * CF;
    const double K2qq = 4 * CF;

    const Operator orlqqx{gx,  rlqqpol{}};
    const Operator orsqqx{gx,  rsqqpol{}};
    const Operator or11qqx{gx, r11qqpol{}};
    const Operator or21qqx{gx, r21qqpol{}};
    const Operator olrqqz{gz,  lrqqpol{}};
    const Operator osrqqz{gz,  srqqpol{}};
    const Operator or12qqz{gz, r12qqpol{}};
    const Operator or22qqz{gz, r22qqpol{}};

    if (std::find(exclude.begin(), exclude.end(), 1)  == exclude.end()) SidisObj.G11qq.AddTerm({LLqq, odeltax, odeltaz});
    if (std::find(exclude.begin(), exclude.end(), 2)  == exclude.end()) SidisObj.G11qq.AddTerm({LSqq, odeltax, oD1z   });
    if (std::find(exclude.begin(), exclude.end(), 3)  == exclude.end()) SidisObj.G11qq.AddTerm({1,    odeltax, olrqqz });
    if (std::find(exclude.begin(), exclude.end(), 4)  == exclude.end()) SidisObj.G11qq.AddTerm({SLqq, oD1x,    odeltaz});
    if (std::find(exclude.begin(), exclude.end(), 5)  == exclude.end()) SidisObj.G11qq.AddTerm({SSqq, oD0x,    oD0z   });
    if (std::find(exclude.begin(), exclude.end(), 6)  == exclude.end()) SidisObj.G11qq.AddTerm({1,    oD0x,    osrqqz });
    if (std::find(exclude.begin(), exclude.end(), 7)  == exclude.end()) SidisObj.G11qq.AddTerm({1,    orlqqx,  odeltaz});
    if (std::find(exclude.begin(), exclude.end(), 8)  == exclude.end()) SidisObj.G11qq.AddTerm({1,    orsqqx,  oD0z   });
    if (std::find(exclude.begin(), exclude.end(), 9)  == exclude.end()) SidisObj.G11qq.AddTerm({K1qq, or11qqx, or12qqz});
    if (std::find(exclude.begin(), exclude.end(), 10) == exclude.end()) SidisObj.G11qq.AddTerm({K2qq, or21qqx, or22qqz});

    const double K1gq = 2 * CF;
    const double K2gq = - 4 * CF;

    const Operator olrgqz{gz,  lrgqpol{}};
    const Operator osrgqz{gz,  srgqpol{}};
    const Operator or11gqx{gx, r11gqpol{}};
    const Operator or12gqz{gz, r12gqpol{}};
    const Operator or21gqx{gx, r21gqpol{}};
    const Operator or22gqz{gz, r22gqpol{}};

    if (std::find(exclude.begin(), exclude.end(), 11) == exclude.end()) SidisObj.G11gq.AddTerm({1,    odeltax, olrgqz });
    if (std::find(exclude.begin(), exclude.end(), 12) == exclude.end()) SidisObj.G11gq.AddTerm({1,    oD0x,    osrgqz });
    if (std::find(exclude.begin(), exclude.end(), 13) == exclude.end()) SidisObj.G11gq.AddTerm({K1gq, or11gqx, or12gqz});
    if (std::find(exclude.begin(), exclude.end(), 14) == exclude.end()) SidisObj.G11gq.AddTerm({K2gq, or21gqx, or22gqz});

    const double K1qg = 1;

    const Operator orlqgx{gx,  rlqgpol{}};
    const Operator orsqgx{gx,  rsqgpol{}};
    const Operator or11qgx{gx, r11qgpol{}};
    const Operator or12qgz{gz, r12qgpol{}};

    if (std::find(exclude.begin(), exclude.end(), 15) == exclude.end()) SidisObj.G11qg.AddTerm({1,    orlqgx,  odeltaz});
    if (std::find(exclude.begin(), exclude.end(), 16) == exclude.end()) SidisObj.G11qg.AddTerm({1,    orsqgx,  oD0z   });
    if (std::find(exclude.begin(), exclude.end(), 17) == exclude.end()) SidisObj.G11qg.AddTerm({K1qg, or11qgx, or12qgz});

    // ====================================================
    // Exact NNLO corrections.
    // Polarised G1 expressions from arXiv:2404.08597.
    // ====================================================
    // NNLO G1: nf-independent channels
    SidisObj.G12.emplace("gg",   DoubleOperator{gx, gz, DC2G2G{},   IntEps});
    SidisObj.G12.emplace("ps",   DoubleOperator{gx, gz, DC2Q2QPS{}, IntEps});
    SidisObj.G12.emplace("qbq",  DoubleOperator{gx, gz, DC2Q2QB{},  IntEps});
    SidisObj.G12.emplace("qpq1", DoubleOperator{gx, gz, DC2Q2QP1{}, IntEps});
    SidisObj.G12.emplace("qpq2", DoubleOperator{gx, gz, DC2Q2QP2{}, IntEps});
    SidisObj.G12.emplace("qpq3", DoubleOperator{gx, gz, DC2Q2QP3{}, IntEps});

    // NNLO G1: nf-dependent channels (all three are nf-dependent for polarised)
    for (int nf = nfi; nf <= nff; nf++)
      {
        SidisObj.G12ns.emplace(nf, DoubleOperator{gx, gz, DC2Q2QNS{nf}, IntEps});
        SidisObj.G12gq.emplace(nf, DoubleOperator{gx, gz, DC2Q2G{nf},   IntEps});
        SidisObj.G12qg.emplace(nf, DoubleOperator{gx, gz, DC2G2Q{nf},   IntEps});
      }
    t.stop();

    return SidisObj;
  }

  // Functions that fills in the SIDIS hard cross sections on one
  // single grid and exchanges the last defaulted arguments.
  SidisPolObjects InitializeSIDISpol(Grid const& gx, std::vector<double> const& Thresholds, std::vector<int> exclude = {}, double const& IntEps = 1e-3)
  {
    return InitializeSIDISpol(gx, gx, Thresholds, exclude, IntEps);
  }
}
