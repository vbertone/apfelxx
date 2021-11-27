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
    DoubleObject<Operator> C20qq;
    DoubleObject<Operator> C21qq;
    DoubleObject<Operator> C21gq;
    DoubleObject<Operator> C21qg;
    DoubleObject<Operator> CL1qq;
    DoubleObject<Operator> CL1gq;
    DoubleObject<Operator> CL1qg;
    std::map<int, DoubleObject<Operator>> C22qq;
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

  class Dn: public Expression
  {
  public:
    Dn(int const& n): Expression(), _n(n) {}
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

  // Expressions needed for the computation of the SIDIS cross sections.
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
  SidisObjects InitializeSIDIS(Grid const& gx, Grid const& gz, std::vector<double> const& Thresholds, std::vector<int> exclude = {})
  {
    report("Initializing SIDIS hard cross sections... ");
    Timer t;

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the thresholds
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
        nfi++;

    // Define object of the structure containing the DglapObjects.
    SidisObjects SidisObj;

    // ====================================================
    // Compute SIDIS partonic cross sections.
    // Expressions taken from Appendix C of hep-ph/9711387.
    // ====================================================
    // LO contribution.
    const Operator odeltax{gx, delta{}};
    const Operator odeltaz{gz, delta{}};

    if (std::find(exclude.begin(), exclude.end(), 0) == exclude.end()) SidisObj.C20qq.AddTerm({1, odeltax, odeltaz});  //0

    // NLO contributions
    // F2
    const Operator oD0x{gx, Dn{0}};
    const Operator oD1x{gx, Dn{1}};
    const Operator oD0z{gz, Dn{0}};
    const Operator oD1z{gz, Dn{1}};

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

    if (std::find(exclude.begin(), exclude.end(), 1)  == exclude.end()) SidisObj.C21qq.AddTerm({LLqq, odeltax, odeltaz}); //1
    if (std::find(exclude.begin(), exclude.end(), 2)  == exclude.end()) SidisObj.C21qq.AddTerm({LSqq, odeltax, oD1z   }); //2
    if (std::find(exclude.begin(), exclude.end(), 3)  == exclude.end()) SidisObj.C21qq.AddTerm({1,    odeltax, olrqqz});  //3
    if (std::find(exclude.begin(), exclude.end(), 4)  == exclude.end()) SidisObj.C21qq.AddTerm({SLqq, oD1x,    odeltaz}); //4
    if (std::find(exclude.begin(), exclude.end(), 5)  == exclude.end()) SidisObj.C21qq.AddTerm({SSqq, oD0x,    oD0z   }); //5
    if (std::find(exclude.begin(), exclude.end(), 6)  == exclude.end()) SidisObj.C21qq.AddTerm({1,    oD0x,    osrqqz});  //6
    if (std::find(exclude.begin(), exclude.end(), 7)  == exclude.end()) SidisObj.C21qq.AddTerm({1,    orlqqx,  odeltaz}); //7
    if (std::find(exclude.begin(), exclude.end(), 8)  == exclude.end()) SidisObj.C21qq.AddTerm({1,    orsqqx,  oD0z   }); //8
    if (std::find(exclude.begin(), exclude.end(), 9)  == exclude.end()) SidisObj.C21qq.AddTerm({K1qq, or11qqx, or12qqz}); //9
    if (std::find(exclude.begin(), exclude.end(), 10) == exclude.end()) SidisObj.C21qq.AddTerm({K2qq, or21qqx, or22qqz}); //10

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

    if (std::find(exclude.begin(), exclude.end(), 11) == exclude.end()) SidisObj.C21gq.AddTerm({1,    odeltax, olrgqz});  //11
    if (std::find(exclude.begin(), exclude.end(), 12) == exclude.end()) SidisObj.C21gq.AddTerm({1,    oD0x,    osrgqz});  //12
    if (std::find(exclude.begin(), exclude.end(), 13) == exclude.end()) SidisObj.C21gq.AddTerm({K1gq, or11gqx, or12gqz}); //13
    if (std::find(exclude.begin(), exclude.end(), 14) == exclude.end()) SidisObj.C21gq.AddTerm({K2gq, or21gqx, or22gqz}); //14
    if (std::find(exclude.begin(), exclude.end(), 15) == exclude.end()) SidisObj.C21gq.AddTerm({K3gq, or31gqx, or32gqz}); //15

    const double K1qg = 2;
    const double K2qg = 1;

    const Operator orlqgx{gx,  rlqg{}};
    const Operator orsqgx{gx,  rsqg{}};
    const Operator or11qgx{gx, r11qg{}};
    const Operator or12qgz{gz, r12qg{}};
    const Operator or21qgx{gx, r21qg{}};
    const Operator or22qgz{gz, r22qg{}};

    if (std::find(exclude.begin(), exclude.end(), 16) == exclude.end()) SidisObj.C21qg.AddTerm({1,    orlqgx,  odeltaz}); //16
    if (std::find(exclude.begin(), exclude.end(), 17) == exclude.end()) SidisObj.C21qg.AddTerm({1,    orsqgx,  oD0z   }); //17
    if (std::find(exclude.begin(), exclude.end(), 18) == exclude.end()) SidisObj.C21qg.AddTerm({K1qg, or11qgx, or12qgz}); //18
    if (std::find(exclude.begin(), exclude.end(), 19) == exclude.end()) SidisObj.C21qg.AddTerm({K2qg, or21qgx, or22qgz}); //19

    // FL
    const double K1Lqq = 8 * CF;

    const Operator or11Lqqx{gx, r11Lqq{}};
    const Operator or12Lqqz{gz, r12Lqq{}};

    if (std::find(exclude.begin(), exclude.end(), 20) == exclude.end()) SidisObj.CL1qq.AddTerm({K1Lqq, or11Lqqx, or12Lqqz});  //20

    const double K1Lgq = 8 * CF;

    const Operator or11Lgqx{gx, r11Lgq{}};
    const Operator or12Lgqz{gz, r12Lgq{}};

    if (std::find(exclude.begin(), exclude.end(), 21) == exclude.end()) SidisObj.CL1gq.AddTerm({K1Lgq, or11Lgqx, or12Lgqz});  //21

    const double K1Lqg = 8;

    const Operator or11Lqgx{gx, r11Lqg{}};
    const Operator or12Lqgz{gz, r12Lqg{}};

    if (std::find(exclude.begin(), exclude.end(), 22) == exclude.end()) SidisObj.CL1qg.AddTerm({K1Lqg, or11Lqgx, or12Lqgz});  //22

    // ====================================================
    // Approximated NNLO corrections derived from threshold
    // resummation. They only contribute to the qq channel of F2.
    // Expressions taken from Appendix B of arXiv:2109.00847.
    // ====================================================
    // Additional singular terms
    const Operator oD2x{gx, Dn{2}};
    const Operator oD3x{gx, Dn{3}};
    const Operator oD2z{gz, Dn{2}};
    const Operator oD3z{gz, Dn{3}};

    // Non-singular (next-to-leading power) terms
    const Operator ol1x{gx, ln{1}};
    const Operator ol2x{gx, ln{2}};
    const Operator ol3x{gx, ln{3}};
    const Operator ol1z{gz, ln{1}};
    const Operator ol2z{gz, ln{2}};
    const Operator ol3z{gz, ln{3}};

    // Constant function
    const Operator oonex{gx, one{}};
    const Operator oonez{gz, one{}};

    // Overall constant (16 = 4^2 due to the different normalisation
    // of the expansion parameter)
    const double ovc = 16. * CF;

    // Loop over nf
    for (int nf = nfi; nf <= nff; nf++)
      {
        DoubleObject<Operator> cnf{};

        // CF and leading-power terms (Eq. 62)
        const double c1 = ovc * CF / 2.;
        if (std::find(exclude.begin(), exclude.end(), 23) == exclude.end()) cnf.AddTerm({c1, odeltax, oD3z});
        if (std::find(exclude.begin(), exclude.end(), 24) == exclude.end()) cnf.AddTerm({c1, oD3x, odeltaz});

        const double c2 = ovc * CF * 3. / 2.;
        if (std::find(exclude.begin(), exclude.end(), 25) == exclude.end()) cnf.AddTerm({c2, oD0x, oD2z});
        if (std::find(exclude.begin(), exclude.end(), 26) == exclude.end()) cnf.AddTerm({c2, oD2x, oD0z});
        if (std::find(exclude.begin(), exclude.end(), 27) == exclude.end()) cnf.AddTerm({2 * c2, oD1x, oD1z});

        const double c3 = - ovc * CF * ( 4. + Pi2 / 3. );
        if (std::find(exclude.begin(), exclude.end(), 28) == exclude.end()) cnf.AddTerm({c3, oD0x, oD0z});
        if (std::find(exclude.begin(), exclude.end(), 29) == exclude.end()) cnf.AddTerm({c3, odeltax, oD1z});
        if (std::find(exclude.begin(), exclude.end(), 30) == exclude.end()) cnf.AddTerm({c3, oD1x, odeltaz});

        const double c4 = ovc * CF * 2. * zeta3;
        if (std::find(exclude.begin(), exclude.end(), 31) == exclude.end()) cnf.AddTerm({c4, odeltax, oD0z});
        if (std::find(exclude.begin(), exclude.end(), 32) == exclude.end()) cnf.AddTerm({c4, oD0x, odeltaz});

        const double c5 = ovc * CF * ( 511. / 64. - 15. * zeta3 / 4. + 29. * Pi2 / 48. - 7. * Pi2 * Pi2 / 360. );
        if (std::find(exclude.begin(), exclude.end(), 33) == exclude.end()) cnf.AddTerm({c5, odeltax, odeltaz});

        // CF and next-to-leading-power terms (Eq. 63)
        const double c6 = - ovc * CF * 3. / 2.;
        if (std::find(exclude.begin(), exclude.end(), 34) == exclude.end()) cnf.AddTerm({c6, oD2x, oonez});
        if (std::find(exclude.begin(), exclude.end(), 35) == exclude.end()) cnf.AddTerm({c6, oonex, oD2z});
        if (std::find(exclude.begin(), exclude.end(), 36) == exclude.end()) cnf.AddTerm({2 * c6, oD1x, ol1z});
        if (std::find(exclude.begin(), exclude.end(), 37) == exclude.end()) cnf.AddTerm({2 * c6, ol1x, oD1z});
        if (std::find(exclude.begin(), exclude.end(), 38) == exclude.end()) cnf.AddTerm({c6, oD0x, ol2z});
        if (std::find(exclude.begin(), exclude.end(), 39) == exclude.end()) cnf.AddTerm({c6, ol2x, oD0z});

        const double c7 = - ovc * CF / 2.;
        if (std::find(exclude.begin(), exclude.end(), 40) == exclude.end()) cnf.AddTerm({c7, odeltax, ol3z});
        if (std::find(exclude.begin(), exclude.end(), 41) == exclude.end()) cnf.AddTerm({c7, ol3x, odeltaz});

        // CA terms (Eq. 64)
        const double c8 = - ovc * CA * 11. / 24.;
        if (std::find(exclude.begin(), exclude.end(), 42) == exclude.end()) cnf.AddTerm({c8, odeltax, oD2z});
        if (std::find(exclude.begin(), exclude.end(), 43) == exclude.end()) cnf.AddTerm({c8, oD2x, odeltaz});
        if (std::find(exclude.begin(), exclude.end(), 44) == exclude.end()) cnf.AddTerm({2 * c8, oD0x, oD1z});
        if (std::find(exclude.begin(), exclude.end(), 45) == exclude.end()) cnf.AddTerm({2 * c8, oD1x, oD0z});

        const double c9 = ovc * CA * ( 67. / 36. - Pi2 / 12. );
        if (std::find(exclude.begin(), exclude.end(), 46) == exclude.end()) cnf.AddTerm({c9, oD0x, oD0z});
        if (std::find(exclude.begin(), exclude.end(), 47) == exclude.end()) cnf.AddTerm({c9, odeltax, oD1z});
        if (std::find(exclude.begin(), exclude.end(), 48) == exclude.end()) cnf.AddTerm({c9, oD1x, odeltaz});

        const double c10 = ovc * CA * ( 7. * zeta3 / 4. + 11. * Pi2 / 72. - 101. / 54. );
        if (std::find(exclude.begin(), exclude.end(), 49) == exclude.end()) cnf.AddTerm({c10, odeltax, oD0z});
        if (std::find(exclude.begin(), exclude.end(), 50) == exclude.end()) cnf.AddTerm({c10, oD0x, odeltaz});

        const double c11 = ovc * CA * ( 43. * zeta3 / 12. + 17. * Pi2 * Pi2 / 720. - 1535. / 192. - 269. * Pi2 / 432. );
        if (std::find(exclude.begin(), exclude.end(), 51) == exclude.end()) cnf.AddTerm({c11, odeltax, odeltaz});

        // nf terms (Eq. 65)
        const double c12 = ovc * nf / 12.;
        if (std::find(exclude.begin(), exclude.end(), 52) == exclude.end()) cnf.AddTerm({c12, odeltax, oD2z});
        if (std::find(exclude.begin(), exclude.end(), 53) == exclude.end()) cnf.AddTerm({c12, oD2z, odeltaz});
        if (std::find(exclude.begin(), exclude.end(), 54) == exclude.end()) cnf.AddTerm({2 * c12, oD0x, oD1z});
        if (std::find(exclude.begin(), exclude.end(), 55) == exclude.end()) cnf.AddTerm({2 * c12, oD1x, oD0z});

        const double c13 = - ovc * nf * 5. / 18.;
        if (std::find(exclude.begin(), exclude.end(), 56) == exclude.end()) cnf.AddTerm({c13, oD0x, oD0z});
        if (std::find(exclude.begin(), exclude.end(), 57) == exclude.end()) cnf.AddTerm({c13, odeltax, oD1z});
        if (std::find(exclude.begin(), exclude.end(), 58) == exclude.end()) cnf.AddTerm({c13, oD1x, odeltaz});

        const double c14 = ovc * nf * ( 7. / 27. - Pi2 / 36. );
        if (std::find(exclude.begin(), exclude.end(), 59) == exclude.end()) cnf.AddTerm({c14, odeltax, oD0z});
        if (std::find(exclude.begin(), exclude.end(), 60) == exclude.end()) cnf.AddTerm({c14, oD0x, odeltaz});

        const double c15 = ovc * nf * ( zeta3 / 6. + 19. * Pi2 / 216. + 127. / 96. );
        if (std::find(exclude.begin(), exclude.end(), 61) == exclude.end()) cnf.AddTerm({c15, odeltax, odeltaz});

        SidisObj.C22qq.insert({nf, cnf});
      }
    t.stop();

    return SidisObj;
  }

  // Functions that fills in the SIDIS hard cross sections on one
  // single grid and exchanges the last defaulted arguments.
  SidisObjects InitializeSIDIS(Grid const& gx, std::vector<double> const& Thresholds, std::vector<int> exclude = {})
  {
    return InitializeSIDIS(gx, gx, Thresholds, exclude);
  }
}
