//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/apfelxx.h"

namespace apfel
{
  /**
   * @brief Structure that contains the precomputed SIDIS hard cross
   * sections.
   */
  struct SidisObjects
  {
    apfel::DoubleObject<apfel::Operator> C20qq;
    apfel::DoubleObject<apfel::Operator> C21qq;
    apfel::DoubleObject<apfel::Operator> C21gq;
    apfel::DoubleObject<apfel::Operator> C21qg;
    apfel::DoubleObject<apfel::Operator> CL1qq;
    apfel::DoubleObject<apfel::Operator> CL1gq;
    apfel::DoubleObject<apfel::Operator> CL1qg;
  };

  // Expressions needed for the computation of the SIDIS cross sections.
  // F2
  class delta: public apfel::Expression
  {
  public:
    delta(): Expression() {}
    double Local(double const&) const { return 1; }
  };

  class s0: public apfel::Expression
  {
  public:
    s0(): Expression() {}
    double Singular(double const& x) const { return 1 / ( 1 - x ); }
    double Local(double const& x) const { return log( 1 - x ); }
  };

  class s1: public apfel::Expression
  {
  public:
    s1(): Expression() {}
    double Singular(double const& x) const { return log( 1 - x ) / ( 1 - x ); }
    double Local(double const& x) const { return pow(log( 1 - x ), 2) / 2; }
  };

  class lrqq: public apfel::Expression
  {
  public:
    lrqq(): Expression() {}
    double Regular(double const& x) const
    {
      const double expr = ( 1 + x * x ) * log(x) / ( 1 - x ) + 1 - x - ( 1 + x ) * log( 1 - x );
      return 2 * apfel::CF * expr;
    }
  };

  class srqq: public apfel::Expression
  {
  public:
    srqq(): Expression() {}
    double Regular(double const& x) const { return - 2 * apfel::CF * ( 1 + x ); }
  };

  class rlqq: public apfel::Expression
  {
  public:
    rlqq(): Expression() {}
    double Regular(double const& x) const
    {
      const double expr = - ( 1 + x * x ) * log(x) / ( 1 - x ) + 1 - x - ( 1 + x ) * log( 1 - x );
      return 2 * apfel::CF * expr;
    }
  };

  class rsqq: public apfel::Expression
  {
  public:
    rsqq(): Expression() {}
    double Regular(double const& x) const { return - 2 * apfel::CF * ( 1 + x ); }
  };

  class r11qq: public apfel::Expression
  {
  public:
    r11qq(): Expression() {}
    double Regular(double const&) const { return 1; }
  };

  class r12qq: public apfel::Expression
  {
  public:
    r12qq(): Expression() {}
    double Regular(double const&) const { return 1; }
  };

  class r21qq: public apfel::Expression
  {
  public:
    r21qq(): Expression() {}
    double Regular(double const& x) const { return x; }
  };

  class r22qq: public apfel::Expression
  {
  public:
    r22qq(): Expression() {}
    double Regular(double const& x) const { return x; }
  };

  class lrgq: public apfel::Expression
  {
  public:
    lrgq(): Expression() {}
    double Regular(double const& x) const
    {
      const double omx = ( 1 - x );
      const double expr = ( 1 + omx * omx ) * log( x * omx ) / x + x;
      return 2 * apfel::CF * expr;
    }
  };

  class srgq: public apfel::Expression
  {
  public:
    srgq(): Expression() {}
    double Regular(double const& x) const
    {
      const double omx = ( 1 - x );
      const double expr = ( 1 + omx * omx ) / x;
      return 2 * apfel::CF * expr;
    }
  };

  class r11gq: public apfel::Expression
  {
  public:
    r11gq(): Expression() {}
    double Regular(double const& x) const { return 1 + 3 * x; }
  };

  class r12gq: public apfel::Expression
  {
  public:
    r12gq(): Expression() {}
    double Regular(double const&) const { return 1; }
  };

  class r21gq: public apfel::Expression
  {
  public:
    r21gq(): Expression() {}
    double Regular(double const& x) const { return x; }
  };

  class r22gq: public apfel::Expression
  {
  public:
    r22gq(): Expression() {}
    double Regular(double const& x) const { return x; }
  };

  class r31gq: public apfel::Expression
  {
  public:
    r31gq(): Expression() {}
    double Regular(double const& x) const { return 1 + x; }
  };

  class r32gq: public apfel::Expression
  {
  public:
    r32gq(): Expression() {}
    double Regular(double const& x) const { return 1 / x; }
  };

  class rlqg: public apfel::Expression
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

  class rsqg: public apfel::Expression
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

  class r11qg: public apfel::Expression
  {
  public:
    r11qg(): Expression() {}
    double Regular(double const& x) const { return - 1 + 6 * x - 6 * x * x; }
  };

  class r12qg: public apfel::Expression
  {
  public:
    r12qg(): Expression() {}
    double Regular(double const&) const { return 1; }
  };

  class r21qg: public apfel::Expression
  {
  public:
    r21qg(): Expression() {}
    double Regular(double const& x) const { return x * x + ( 1 - x ) * ( 1 - x ); }
  };

  class r22qg: public apfel::Expression
  {
  public:
    r22qg(): Expression() {}
    double Regular(double const& x) const { return 1 / x; }
  };

  // FL
  class r11Lqq: public apfel::Expression
  {
  public:
    r11Lqq(): Expression() {}
    double Regular(double const& x) const { return x; }
  };

  class r12Lqq: public apfel::Expression
  {
  public:
    r12Lqq(): Expression() {}
    double Regular(double const& x) const { return x; }
  };

  class r11Lgq: public apfel::Expression
  {
  public:
    r11Lgq(): Expression() {}
    double Regular(double const& x) const { return x; }
  };

  class r12Lgq: public apfel::Expression
  {
  public:
    r12Lgq(): Expression() {}
    double Regular(double const& x) const { return 1 - x; }
  };

  class r11Lqg: public apfel::Expression
  {
  public:
    r11Lqg(): Expression() {}
    double Regular(double const& x) const { return x * ( 1 - x ); }
  };

  class r12Lqg: public apfel::Expression
  {
  public:
    r12Lqg(): Expression() {}
    double Regular(double const&) const { return 1; }
  };

  // Functions that fills in the SIDIS hard cross sections on two
  // different grids.
  SidisObjects InitializeSIDIS(apfel::Grid const& gx, apfel::Grid const& gz, std::vector<int> exclude = {})
  {
    apfel::report("Initializing SIDIS hard cross sections... ");
    apfel::Timer t;

    // Define object of the structure containing the DglapObjects.
    SidisObjects SidisObj;

    // ====================================================
    // Compute SIDIS partonic cross sections.
    // Expressions taken from Appendix C of hep-ph/9711387.
    // ====================================================
    // LO contribution.
    const apfel::Operator odeltax{gx, delta{}};
    const apfel::Operator odeltaz{gz, delta{}};

    if (std::find(exclude.begin(), exclude.end(), 0) == exclude.end()) SidisObj.C20qq.AddTerm({1, odeltax, odeltaz});  //0

    // NLO contributions
    // F2
    const apfel::Operator os0x{gx, s0{}};
    const apfel::Operator os1x{gx, s1{}};
    const apfel::Operator os0z{gz, s0{}};
    const apfel::Operator os1z{gz, s1{}};

    const double LLqq = - 16 * apfel::CF;
    const double LSqq = 4 * apfel::CF;
    const double SLqq = 4 * apfel::CF;
    const double SSqq = 4 * apfel::CF;
    const double K1qq = 4 * apfel::CF;
    const double K2qq = 12 * apfel::CF;

    const apfel::Operator orlqqx{gx,  rlqq{}};
    const apfel::Operator orsqqx{gx,  rsqq{}};
    const apfel::Operator or11qqx{gx, r11qq{}};
    const apfel::Operator or21qqx{gx, r21qq{}};
    const apfel::Operator olrqqz{gz,  lrqq{}};
    const apfel::Operator osrqqz{gz,  srqq{}};
    const apfel::Operator or12qqz{gz, r12qq{}};
    const apfel::Operator or22qqz{gz, r22qq{}};

    if (std::find(exclude.begin(), exclude.end(), 1)  == exclude.end()) SidisObj.C21qq.AddTerm({LLqq, odeltax, odeltaz}); //1
    if (std::find(exclude.begin(), exclude.end(), 2)  == exclude.end()) SidisObj.C21qq.AddTerm({LSqq, odeltax, os1z   }); //2
    if (std::find(exclude.begin(), exclude.end(), 3)  == exclude.end()) SidisObj.C21qq.AddTerm({1,    odeltax, olrqqz});  //3
    if (std::find(exclude.begin(), exclude.end(), 4)  == exclude.end()) SidisObj.C21qq.AddTerm({SLqq, os1x,    odeltaz}); //4
    if (std::find(exclude.begin(), exclude.end(), 5)  == exclude.end()) SidisObj.C21qq.AddTerm({SSqq, os0x,    os0z   }); //5
    if (std::find(exclude.begin(), exclude.end(), 6)  == exclude.end()) SidisObj.C21qq.AddTerm({1,    os0x,    osrqqz});  //6
    if (std::find(exclude.begin(), exclude.end(), 7)  == exclude.end()) SidisObj.C21qq.AddTerm({1,    orlqqx,  odeltaz}); //7
    if (std::find(exclude.begin(), exclude.end(), 8)  == exclude.end()) SidisObj.C21qq.AddTerm({1,    orsqqx,  os0z   }); //8
    if (std::find(exclude.begin(), exclude.end(), 9)  == exclude.end()) SidisObj.C21qq.AddTerm({K1qq, or11qqx, or12qqz}); //9
    if (std::find(exclude.begin(), exclude.end(), 10) == exclude.end()) SidisObj.C21qq.AddTerm({K2qq, or21qqx, or22qqz}); //10

    const double K1gq = 4 * apfel::CF;
    const double K2gq = - 12 * apfel::CF;
    const double K3gq = - 2 * apfel::CF;

    const apfel::Operator olrgqz{gz,  lrgq{}};
    const apfel::Operator osrgqz{gz,  srgq{}};
    const apfel::Operator or11gqx{gx, r11gq{}};
    const apfel::Operator or12gqz{gz, r12gq{}};
    const apfel::Operator or21gqx{gx, r21gq{}};
    const apfel::Operator or22gqz{gz, r22gq{}};
    const apfel::Operator or31gqx{gx, r31gq{}};
    const apfel::Operator or32gqz{gz, r32gq{}};

    if (std::find(exclude.begin(), exclude.end(), 11) == exclude.end()) SidisObj.C21gq.AddTerm({1,    odeltax, olrgqz});  //11
    if (std::find(exclude.begin(), exclude.end(), 12) == exclude.end()) SidisObj.C21gq.AddTerm({1,    os0x,    osrgqz});  //12
    if (std::find(exclude.begin(), exclude.end(), 13) == exclude.end()) SidisObj.C21gq.AddTerm({K1gq, or11gqx, or12gqz}); //13
    if (std::find(exclude.begin(), exclude.end(), 14) == exclude.end()) SidisObj.C21gq.AddTerm({K2gq, or21gqx, or22gqz}); //14
    if (std::find(exclude.begin(), exclude.end(), 15) == exclude.end()) SidisObj.C21gq.AddTerm({K3gq, or31gqx, or32gqz}); //15

    const double K1qg = 2;
    const double K2qg = 1;

    const apfel::Operator orlqgx{gx,  rlqg{}};
    const apfel::Operator orsqgx{gx,  rsqg{}};
    const apfel::Operator or11qgx{gx, r11qg{}};
    const apfel::Operator or12qgz{gz, r12qg{}};
    const apfel::Operator or21qgx{gx, r21qg{}};
    const apfel::Operator or22qgz{gz, r22qg{}};

    if (std::find(exclude.begin(), exclude.end(), 16) == exclude.end()) SidisObj.C21qg.AddTerm({1,    orlqgx,  odeltaz}); //16
    if (std::find(exclude.begin(), exclude.end(), 17) == exclude.end()) SidisObj.C21qg.AddTerm({1,    orsqgx,  os0z   }); //17
    if (std::find(exclude.begin(), exclude.end(), 18) == exclude.end()) SidisObj.C21qg.AddTerm({K1qg, or11qgx, or12qgz}); //18
    if (std::find(exclude.begin(), exclude.end(), 19) == exclude.end()) SidisObj.C21qg.AddTerm({K2qg, or21qgx, or22qgz}); //19

    // FL
    const double K1Lqq = 8 * apfel::CF;

    const apfel::Operator or11Lqqx{gx, r11Lqq{}};
    const apfel::Operator or12Lqqz{gz, r12Lqq{}};

    if (std::find(exclude.begin(), exclude.end(), 20) == exclude.end()) SidisObj.CL1qq.AddTerm({K1Lqq, or11Lqqx, or12Lqqz});  //20

    const double K1Lgq = 8 * apfel::CF;

    const apfel::Operator or11Lgqx{gx, r11Lgq{}};
    const apfel::Operator or12Lgqz{gz, r12Lgq{}};

    if (std::find(exclude.begin(), exclude.end(), 21) == exclude.end()) SidisObj.CL1gq.AddTerm({K1Lgq, or11Lgqx, or12Lgqz});  //21

    const double K1Lqg = 8;

    const apfel::Operator or11Lqgx{gx, r11Lqg{}};
    const apfel::Operator or12Lqgz{gz, r12Lqg{}};

    if (std::find(exclude.begin(), exclude.end(), 22) == exclude.end()) SidisObj.CL1qg.AddTerm({K1Lqg, or11Lqgx, or12Lqgz});  //22
    t.stop();

    return SidisObj;
  }

  // Functions that fills in the SIDIS hard cross sections on one
  // single grid.
  SidisObjects InitializeSIDIS(apfel::Grid const& gx, std::vector<int> exclude = {})
  {
    return InitializeSIDIS(gx, gx, exclude);
  }
}
