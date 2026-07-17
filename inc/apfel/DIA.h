//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/apfelxx.h"
#include "apfel/SIDIS.h"

namespace apfel
{
  /**
   * @brief Structure that contains the precomputed DIA hard cross
   * sections.
   */
  struct DiaObjects
  {
    DoubleObject<Operator> C0qq;
    DoubleObject<Operator> C1qq;
    DoubleObject<Operator> C1gq;
    DoubleObject<Operator> C1qg;
  };

  // Expressions needed for the computation of the DIA cross sections
  class lrqqdia: public Expression
  {
  public:
    lrqqdia(): Expression() {}
    double Regular(double const& u) const
    {
      const double expr = ( 1 + u * u ) * log(u) / ( 1 - u ) + 1 - u - ( 1 + u ) * log(1 - u);
      return 2 * CF * expr;
    }
  };

  class rlqqdia: public Expression
  {
  public:
    rlqqdia(): Expression() {}
    double Regular(double const& z) const
    {
      const double expr = 2 * ( 1 + z * z ) * log(z) / ( 1 - z ) + 1 - z - ( 1 + z ) * log(1 - z);
      return 2 * CF * expr;
    }
  };

  class srqqdia: public Expression
  {
  public:
    srqqdia(): Expression() {}
    double Regular(double const& u) const { return - 2 * CF * ( 1 + u ); }
  };

  class rsqqdia: public Expression
  {
  public:
    rsqqdia(): Expression() {}
    double Regular(double const& z) const { return - 2 * CF * ( 1 + z ); }
  };

  class r11qqdia: public Expression
  {
  public:
    r11qqdia(): Expression() {}
    double Regular(double const& z) const { return z; }
  };

  class r12qqdia: public Expression
  {
  public:
    r12qqdia(): Expression() {}
    double Regular(double const&) const { return 1; }
  };

  class r21qqdia: public Expression
  {
  public:
    r21qqdia(): Expression() {}
    double Regular(double const&) const { return 1; }
  };

  class r22qqdia: public Expression
  {
  public:
    r22qqdia(): Expression() {}
    double Regular(double const& u) const { return u; }
  };

  class r31qqdia: public Expression
  {
  public:
    r31qqdia(): Expression() {}
    double Regular(double const& z) const { return z; }
  };

  class r32qqdia: public Expression
  {
  public:
    r32qqdia(): Expression() {}
    double Regular(double const& u) const { return u; }
  };

  class r41qqdia: public Expression
  {
  public:
    r41qqdia(): Expression() {}
    double Regular(double const&) const { return 1; }
  };

  class r42qqdia: public Expression
  {
  public:
    r42qqdia(): Expression() {}
    double Regular(double const&) const { return 1; }
  };

  class lrqgdia: public Expression
  {
  public:
    lrqgdia(): Expression() {}
    double Regular(double const& u) const
    {
      const double omu = ( 1 - u );
      const double expr = ( 1 + omu * omu ) * ( log(u * omu) + 1 ) / u - 2 * omu / u;
      return 2 * CF * expr;
    }
  };

  class srqgdia: public Expression
  {
  public:
    srqgdia(): Expression() {}
    double Regular(double const& u) const
    {
      const double omu = ( 1 - u );
      const double expr = ( 1 + omu * omu ) / u;
      return 2 * CF * expr;
    }
  };

  class r11qgdia: public Expression
  {
  public:
    r11qgdia(): Expression() {}
    double Regular(double const& z) const { return 1; }
  };

  class r12qgdia: public Expression
  {
  public:
    r12qgdia(): Expression() {}
    double Regular(double const& u) const { return 2 - u - 1 / u; }
  };

  class r21qgdia: public Expression
  {
  public:
    r21qgdia(): Expression() {}
    double Regular(double const& z) const { return z; }
  };

  class r22qgdia: public Expression
  {
  public:
    r22qgdia(): Expression() {}
    double Regular(double const& u) const { return - u - 1 / u; }
  };

  class rlgqdia: public Expression
  {
  public:
    rlgqdia(): Expression() {}
    double Regular(double const& z) const
    {
      const double omz = ( 1 - z );
      const double expr = ( 1 + omz * omz ) * ( log(pow(z, 2) * omz) + 1 ) / z - 2 * omz / z;
      return 2 * CF * expr;
    }
  };

  class rsgqdia: public Expression
  {
  public:
    rsgqdia(): Expression() {}
    double Regular(double const& z) const
    {
      const double omz = ( 1 - z );
      const double expr = ( 1 + omz * omz ) / z;
      return 2 * CF * expr;
    }
  };

  class r11gqdia: public Expression
  {
  public:
    r11gqdia(): Expression() {}
    double Regular(double const& z) const { return ( 1 + ( 1 - z ) * ( 1 - z ) ) / z; }
  };

  class r12gqdia: public Expression
  {
  public:
    r12gqdia(): Expression() {}
    double Regular(double const& u) const { return 1 / u; }
  };

  class r21gqdia: public Expression
  {
  public:
    r21gqdia(): Expression() {}
    double Regular(double const& z) const { return z; }
  };

  class r22gqdia: public Expression
  {
  public:
    r22gqdia(): Expression() {}
    double Regular(double const& u) const { return 1; }
  };

  // Functions that fills in the unpolarised DIA hard cross sections
  // on two different grids.
  DiaObjects InitializeDIA(Grid const& gz, Grid const& gu, std::vector<double> const& Thresholds)
  {
    report("Initialising DIA hard cross sections... ");
    Timer t;

    // Define object of the structure containing the DiaObjects
    DiaObjects DiaObj;

    // ====================================================
    // Compute DIA partonic cross sections
    // ====================================================
    // LO contribution
    const Operator odeltaz{gz, delta{}};
    const Operator odeltau{gu, delta{}};

    DiaObj.C0qq.AddTerm({1, odeltaz, odeltau});

    // NLO contributions
    const Operator oD0z{gz, DDn{0}};
    const Operator oD1z{gz, DDn{1}};
    const Operator oD0u{gu, DDn{0}};
    const Operator oD1u{gu, DDn{1}};

    const double LLqq = 2 * CF * ( Pi2 - 8 );
    const double LSqq = 4 * CF;
    const double SLqq = 4 * CF;
    const double SSqq = 4 * CF;

    const double K1qq = - 2 * CF;
    const double K2qq = 2 * CF;
    const double K3qq = 2 * CF;
    const double K4qq = 2 * CF;

    const Operator orlqqz{gz,  rlqqdia{}};
    const Operator orsqqz{gz,  rsqqdia{}};
    const Operator or11qqz{gz, r11qqdia{}};
    const Operator or21qqz{gz, r21qqdia{}};
    const Operator or31qqz{gz, r31qqdia{}};
    const Operator or41qqz{gz, r41qqdia{}};
    const Operator olrqqu{gu,  lrqqdia{}};
    const Operator osrqqu{gu,  srqqdia{}};
    const Operator or12qqu{gu, r12qqdia{}};
    const Operator or22qqu{gu, r22qqdia{}};
    const Operator or32qqu{gu, r32qqdia{}};
    const Operator or42qqu{gu, r42qqdia{}};

    DiaObj.C1qq.AddTerm({LLqq, odeltaz, odeltau});
    DiaObj.C1qq.AddTerm({LSqq, odeltaz, oD1u   });
    DiaObj.C1qq.AddTerm({1,    odeltaz, olrqqu });
    DiaObj.C1qq.AddTerm({SLqq, oD1z,    odeltau});
    DiaObj.C1qq.AddTerm({SSqq, oD0z,    oD0u   });
    DiaObj.C1qq.AddTerm({1,    oD0z,    osrqqu });
    DiaObj.C1qq.AddTerm({1,    orlqqz,  odeltau});
    DiaObj.C1qq.AddTerm({1,    orsqqz,  oD0u   });
    DiaObj.C1qq.AddTerm({K1qq, or11qqz, or12qqu});
    DiaObj.C1qq.AddTerm({K2qq, or21qqz, or22qqu});
    DiaObj.C1qq.AddTerm({K3qq, or31qqz, or32qqu});
    DiaObj.C1qq.AddTerm({K4qq, or41qqz, or42qqu});

    const double K1qg = 2 * CF;
    const double K2qg = 2 * CF;

    const Operator olrqgu {gu,  lrqgdia{}};
    const Operator osrqgu {gu,  srqgdia{}};
    const Operator or11qgz{gz, r11qgdia{}};
    const Operator or12qgu{gu, r12qgdia{}};
    const Operator or21qgz{gz, r21qgdia{}};
    const Operator or22qgu{gu, r22qgdia{}};

    DiaObj.C1qg.AddTerm({1,    odeltaz, olrqgu });
    DiaObj.C1qg.AddTerm({1,    oD0z,    osrqgu });
    DiaObj.C1qg.AddTerm({K1qg, or11qgz, or12qgu});
    DiaObj.C1qg.AddTerm({K2qg, or21qgz, or22qgu});

    const double K1gq = 2 * CF;
    const double K2gq = - 4 * CF;

    const Operator orlgqz{gz,  rlgqdia{}};
    const Operator orsgqz{gz,  rsgqdia{}};
    const Operator or11gqz{gz, r11gqdia{}};
    const Operator or12gqu{gu, r12gqdia{}};
    const Operator or21gqz{gz, r21gqdia{}};
    const Operator or22gqu{gu, r22gqdia{}};

    DiaObj.C1gq.AddTerm({1,    orlgqz,  odeltau});
    DiaObj.C1gq.AddTerm({1,    orsgqz,  oD0u   });
    DiaObj.C1gq.AddTerm({K1gq, or11gqz, or12gqu});
    DiaObj.C1gq.AddTerm({K2gq, or21gqz, or22gqu});
    t.stop();

    return DiaObj;
  }

  // Functions that fills in the DIA unpolarised hard cross sections
  // on one single grid and exchanges the last defaulted arguments.
  DiaObjects InitializeDIA(Grid const& gz, std::vector<double> const& Thresholds)
  {
    return InitializeDIA(gz, gz, Thresholds);
  }
}
