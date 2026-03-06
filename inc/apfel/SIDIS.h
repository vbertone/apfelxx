//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/apfelxx.h"
#include "apfel/sidiscoefficientfunctionsunp.h"

#include <algorithm>

namespace apfel
{
  /**
   * @brief Structure that contains the precomputed SIDIS hard cross
   * sections.
   */
  struct SidisObjects
  {
    DoubleObject<Operator> C20qq;
    DoubleObject<Operator> C21qq;
    DoubleObject<Operator> C21gq;
    DoubleObject<Operator> C21qg;
    DoubleObject<Operator> CL1qq;
    DoubleObject<Operator> CL1gq;
    DoubleObject<Operator> CL1qg;
    std::map<int, DoubleOperator>         C2Tns;
    std::map<std::string, DoubleOperator> C2T;
    std::map<int, DoubleOperator>         C2Lns;
    std::map<std::string, DoubleOperator> C2L;
  };

  // Recurring expressions
  class one: public Expression
  {
  public:
    one(): Expression() {}
    double Regular(double const&) const { return 1; }
  };

  class delta: public Expression
  {
  public:
    delta(): Expression() {}
    double Local(double const&) const { return 1; }
  };

  class DDn: public Expression
  {
  public:
    DDn(int const& n): Expression(), _n(n) {}
    double Singular(double const& x) const { return pow(log( 1 - x ), _n) / ( 1 - x ); }
    double Local(double const& x) const { return pow(log( 1 - x ), _n + 1) / ( _n + 1 ); }
  private:
    int const _n;
  };

  class ln: public Expression
  {
  public:
    ln(int const& n): Expression(), _n(n) {}
    double Regular(double const& x) const { return pow(log( 1 - x ), _n); }
  private:
    int const _n;
  };

  // Expressions needed for the computation of the SIDIS cross sections
  // F2
  class lrqq: public Expression
  {
  public:
    lrqq(): Expression() {}
    double Regular(double const& x) const
    {
      const double expr = ( 1 + x * x ) * log(x) / ( 1 - x ) + 1 - x - ( 1 + x ) * log( 1 - x );
      return 2 * CF * expr;
    }
  };

  class srqq: public Expression
  {
  public:
    srqq(): Expression() {}
    double Regular(double const& x) const { return - 2 * CF * ( 1 + x ); }
  };

  class rlqq: public Expression
  {
  public:
    rlqq(): Expression() {}
    double Regular(double const& x) const
    {
      const double expr = - ( 1 + x * x ) * log(x) / ( 1 - x ) + 1 - x - ( 1 + x ) * log( 1 - x );
      return 2 * CF * expr;
    }
  };

  class rsqq: public Expression
  {
  public:
    rsqq(): Expression() {}
    double Regular(double const& x) const { return - 2 * CF * ( 1 + x ); }
  };

  class r11qq: public Expression
  {
  public:
    r11qq(): Expression() {}
    double Regular(double const&) const { return 1; }
  };

  class r12qq: public Expression
  {
  public:
    r12qq(): Expression() {}
    double Regular(double const&) const { return 1; }
  };

  class r21qq: public Expression
  {
  public:
    r21qq(): Expression() {}
    double Regular(double const& x) const { return x; }
  };

  class r22qq: public Expression
  {
  public:
    r22qq(): Expression() {}
    double Regular(double const& x) const { return x; }
  };

  class lrgq: public Expression
  {
  public:
    lrgq(): Expression() {}
    double Regular(double const& x) const
    {
      const double omx = ( 1 - x );
      const double expr = ( 1 + omx * omx ) * log( x * omx ) / x + x;
      return 2 * CF * expr;
    }
  };

  class srgq: public Expression
  {
  public:
    srgq(): Expression() {}
    double Regular(double const& x) const
    {
      const double omx = ( 1 - x );
      const double expr = ( 1 + omx * omx ) / x;
      return 2 * CF * expr;
    }
  };

  class r11gq: public Expression
  {
  public:
    r11gq(): Expression() {}
    double Regular(double const& x) const { return 1 + 3 * x; }
  };

  class r12gq: public Expression
  {
  public:
    r12gq(): Expression() {}
    double Regular(double const&) const { return 1; }
  };

  class r21gq: public Expression
  {
  public:
    r21gq(): Expression() {}
    double Regular(double const& x) const { return x; }
  };

  class r22gq: public Expression
  {
  public:
    r22gq(): Expression() {}
    double Regular(double const& x) const { return x; }
  };

  class r31gq: public Expression
  {
  public:
    r31gq(): Expression() {}
    double Regular(double const& x) const { return 1 + x; }
  };

  class r32gq: public Expression
  {
  public:
    r32gq(): Expression() {}
    double Regular(double const& x) const { return 1 / x; }
  };

  class rlqg: public Expression
  {
  public:
    rlqg(): Expression() {}
    double Regular(double const& x) const
    {
      const double omx = ( 1 - x );
      const double expr = ( x * x + omx * omx ) * log( omx / x ) + 2 * x * omx;
      return expr;
    }
  };

  class rsqg: public Expression
  {
  public:
    rsqg(): Expression() {}
    double Regular(double const& x) const
    {
      const double omx = ( 1 - x );
      const double expr = ( x * x + omx * omx );
      return expr;
    }
  };

  class r11qg: public Expression
  {
  public:
    r11qg(): Expression() {}
    double Regular(double const& x) const { return - 1 + 6 * x - 6 * x * x; }
  };

  class r12qg: public Expression
  {
  public:
    r12qg(): Expression() {}
    double Regular(double const&) const { return 1; }
  };

  class r21qg: public Expression
  {
  public:
    r21qg(): Expression() {}
    double Regular(double const& x) const { return x * x + ( 1 - x ) * ( 1 - x ); }
  };

  class r22qg: public Expression
  {
  public:
    r22qg(): Expression() {}
    double Regular(double const& x) const { return 1 / x; }
  };

  // FL
  class r11Lqq: public Expression
  {
  public:
    r11Lqq(): Expression() {}
    double Regular(double const& x) const { return x; }
  };

  class r12Lqq: public Expression
  {
  public:
    r12Lqq(): Expression() {}
    double Regular(double const& x) const { return x; }
  };

  class r11Lgq: public Expression
  {
  public:
    r11Lgq(): Expression() {}
    double Regular(double const& x) const { return x; }
  };

  class r12Lgq: public Expression
  {
  public:
    r12Lgq(): Expression() {}
    double Regular(double const& x) const { return 1 - x; }
  };

  class r11Lqg: public Expression
  {
  public:
    r11Lqg(): Expression() {}
    double Regular(double const& x) const { return x * ( 1 - x ); }
  };

  class r12Lqg: public Expression
  {
  public:
    r12Lqg(): Expression() {}
    double Regular(double const&) const { return 1; }
  };

  // Functions that fills in the SIDIS hard cross sections on two
  // different grids.
  SidisObjects InitializeSIDIS(Grid const& gx, Grid const& gz, std::vector<double> const& Thresholds, std::vector<int> exclude = {}, double const& IntEps = 1e-3)
  {
    report("Initializing SIDIS hard cross sections... ");
    Timer t;

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the threshold
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
        nfi++;

    // Define object of the structure containing the DglapObjects
    SidisObjects SidisObj;

    // ====================================================
    // Compute SIDIS partonic cross sections.
    // Expressions taken from Appendix C of hep-ph/9711387.
    // ====================================================
    // LO contribution
    const Operator odeltax{gx, delta{}};
    const Operator odeltaz{gz, delta{}};

    if (std::find(exclude.begin(), exclude.end(), 0) == exclude.end()) SidisObj.C20qq.AddTerm({1, odeltax, odeltaz});

    // NLO contributions
    // F2
    const Operator oD0x{gx, DDn{0}};
    const Operator oD1x{gx, DDn{1}};
    const Operator oD0z{gz, DDn{0}};
    const Operator oD1z{gz, DDn{1}};

    const double LLqq = - 16 * CF;
    const double LSqq = 4 * CF;
    const double SLqq = 4 * CF;
    const double SSqq = 4 * CF;
    const double K1qq = 4 * CF;
    const double K2qq = 12 * CF;

    const Operator orlqqx{gx,  rlqq{}};
    const Operator orsqqx{gx,  rsqq{}};
    const Operator or11qqx{gx, r11qq{}};
    const Operator or21qqx{gx, r21qq{}};
    const Operator olrqqz{gz,  lrqq{}};
    const Operator osrqqz{gz,  srqq{}};
    const Operator or12qqz{gz, r12qq{}};
    const Operator or22qqz{gz, r22qq{}};

    if (std::find(exclude.begin(), exclude.end(), 1)  == exclude.end()) SidisObj.C21qq.AddTerm({LLqq, odeltax, odeltaz});
    if (std::find(exclude.begin(), exclude.end(), 2)  == exclude.end()) SidisObj.C21qq.AddTerm({LSqq, odeltax, oD1z   });
    if (std::find(exclude.begin(), exclude.end(), 3)  == exclude.end()) SidisObj.C21qq.AddTerm({1,    odeltax, olrqqz });
    if (std::find(exclude.begin(), exclude.end(), 4)  == exclude.end()) SidisObj.C21qq.AddTerm({SLqq, oD1x,    odeltaz});
    if (std::find(exclude.begin(), exclude.end(), 5)  == exclude.end()) SidisObj.C21qq.AddTerm({SSqq, oD0x,    oD0z   });
    if (std::find(exclude.begin(), exclude.end(), 6)  == exclude.end()) SidisObj.C21qq.AddTerm({1,    oD0x,    osrqqz });
    if (std::find(exclude.begin(), exclude.end(), 7)  == exclude.end()) SidisObj.C21qq.AddTerm({1,    orlqqx,  odeltaz});
    if (std::find(exclude.begin(), exclude.end(), 8)  == exclude.end()) SidisObj.C21qq.AddTerm({1,    orsqqx,  oD0z   });
    if (std::find(exclude.begin(), exclude.end(), 9)  == exclude.end()) SidisObj.C21qq.AddTerm({K1qq, or11qqx, or12qqz});
    if (std::find(exclude.begin(), exclude.end(), 10) == exclude.end()) SidisObj.C21qq.AddTerm({K2qq, or21qqx, or22qqz});

    const double K1gq = 4 * CF;
    const double K2gq = - 12 * CF;
    const double K3gq = - 2 * CF;

    const Operator olrgqz{gz,  lrgq{}};
    const Operator osrgqz{gz,  srgq{}};
    const Operator or11gqx{gx, r11gq{}};
    const Operator or12gqz{gz, r12gq{}};
    const Operator or21gqx{gx, r21gq{}};
    const Operator or22gqz{gz, r22gq{}};
    const Operator or31gqx{gx, r31gq{}};
    const Operator or32gqz{gz, r32gq{}};

    if (std::find(exclude.begin(), exclude.end(), 11) == exclude.end()) SidisObj.C21gq.AddTerm({1,    odeltax, olrgqz });
    if (std::find(exclude.begin(), exclude.end(), 12) == exclude.end()) SidisObj.C21gq.AddTerm({1,    oD0x,    osrgqz });
    if (std::find(exclude.begin(), exclude.end(), 13) == exclude.end()) SidisObj.C21gq.AddTerm({K1gq, or11gqx, or12gqz});
    if (std::find(exclude.begin(), exclude.end(), 14) == exclude.end()) SidisObj.C21gq.AddTerm({K2gq, or21gqx, or22gqz});
    if (std::find(exclude.begin(), exclude.end(), 15) == exclude.end()) SidisObj.C21gq.AddTerm({K3gq, or31gqx, or32gqz});

    const double K1qg = 2;
    const double K2qg = 1;

    const Operator orlqgx{gx,  rlqg{}};
    const Operator orsqgx{gx,  rsqg{}};
    const Operator or11qgx{gx, r11qg{}};
    const Operator or12qgz{gz, r12qg{}};
    const Operator or21qgx{gx, r21qg{}};
    const Operator or22qgz{gz, r22qg{}};

    if (std::find(exclude.begin(), exclude.end(), 16) == exclude.end()) SidisObj.C21qg.AddTerm({1,    orlqgx,  odeltaz});
    if (std::find(exclude.begin(), exclude.end(), 17) == exclude.end()) SidisObj.C21qg.AddTerm({1,    orsqgx,  oD0z   });
    if (std::find(exclude.begin(), exclude.end(), 18) == exclude.end()) SidisObj.C21qg.AddTerm({K1qg, or11qgx, or12qgz});
    if (std::find(exclude.begin(), exclude.end(), 19) == exclude.end()) SidisObj.C21qg.AddTerm({K2qg, or21qgx, or22qgz});

    // FL
    const double K1Lqq = 8 * CF;

    const Operator or11Lqqx{gx, r11Lqq{}};
    const Operator or12Lqqz{gz, r12Lqq{}};

    if (std::find(exclude.begin(), exclude.end(), 20) == exclude.end()) SidisObj.CL1qq.AddTerm({K1Lqq, or11Lqqx, or12Lqqz});

    const double K1Lgq = 8 * CF;

    const Operator or11Lgqx{gx, r11Lgq{}};
    const Operator or12Lgqz{gz, r12Lgq{}};

    if (std::find(exclude.begin(), exclude.end(), 21) == exclude.end()) SidisObj.CL1gq.AddTerm({K1Lgq, or11Lgqx, or12Lgqz});

    const double K1Lqg = 8;

    const Operator or11Lqgx{gx, r11Lqg{}};
    const Operator or12Lqgz{gz, r12Lqg{}};

    if (std::find(exclude.begin(), exclude.end(), 22) == exclude.end()) SidisObj.CL1qg.AddTerm({K1Lqg, or11Lqgx, or12Lqgz});

    // ====================================================
    // Exact NNLO corrections.
    // Unpolarised FT and FL expressions from arXiv:2401.16281.
    // ====================================================
    // NNLO FT: nf-independent channels
    SidisObj.C2T.emplace("gq",   DoubleOperator{gx, gz, C2TQ2G{},   IntEps});
    SidisObj.C2T.emplace("qg",   DoubleOperator{gx, gz, C2TG2Q{},   IntEps});
    SidisObj.C2T.emplace("gg",   DoubleOperator{gx, gz, C2TG2G{},   IntEps});
    SidisObj.C2T.emplace("ps",   DoubleOperator{gx, gz, C2TQ2QPS{}, IntEps});
    SidisObj.C2T.emplace("qbq",  DoubleOperator{gx, gz, C2TQ2QB{},  IntEps});
    SidisObj.C2T.emplace("qpq1", DoubleOperator{gx, gz, C2TQ2QP1{}, IntEps});
    SidisObj.C2T.emplace("qpq2", DoubleOperator{gx, gz, C2TQ2QP2{}, IntEps});
    SidisObj.C2T.emplace("qpq3", DoubleOperator{gx, gz, C2TQ2QP3{}, IntEps});

    // NNLO FL: nf-independent channels
    SidisObj.C2L.emplace("gq",   DoubleOperator{gx, gz, C2LQ2G{},   IntEps});
    SidisObj.C2L.emplace("qg",   DoubleOperator{gx, gz, C2LG2Q{},   IntEps});
    SidisObj.C2L.emplace("gg",   DoubleOperator{gx, gz, C2LG2G{},   IntEps});
    SidisObj.C2L.emplace("ps",   DoubleOperator{gx, gz, C2LQ2QPS{}, IntEps});
    SidisObj.C2L.emplace("qbq",  DoubleOperator{gx, gz, C2LQ2QB{},  IntEps});
    SidisObj.C2L.emplace("qpq1", DoubleOperator{gx, gz, C2LQ2QP1{}, IntEps});
    SidisObj.C2L.emplace("qpq2", DoubleOperator{gx, gz, C2LQ2QP2{}, IntEps});
    SidisObj.C2L.emplace("qpq3", DoubleOperator{gx, gz, C2LQ2QP3{}, IntEps});

    // NNLO FT and FL: nf-dependent non-singlet
    for (int nf = nfi; nf <= nff; nf++)
      {
        SidisObj.C2Tns.emplace(nf, DoubleOperator{gx, gz, C2TQ2QNS{nf}, IntEps});
        SidisObj.C2Lns.emplace(nf, DoubleOperator{gx, gz, C2LQ2QNS{nf}, IntEps});
      }
    t.stop();

    return SidisObj;
  }

  // Functions that fills in the SIDIS hard cross sections on one
  // single grid and exchanges the last defaulted arguments.
  SidisObjects InitializeSIDIS(Grid const& gx, std::vector<double> const& Thresholds, std::vector<int> exclude = {}, double const& IntEps = 1e-3)
  {
    return InitializeSIDIS(gx, gx, Thresholds, exclude, IntEps);
  }
}
