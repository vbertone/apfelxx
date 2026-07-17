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
  struct DiaTransObjects
  {
    DoubleObject<Operator> dC0qq;
    DoubleObject<Operator> dC1qq;
  };

  // Expressions needed for the computation of the DIA cross sections
  class lrqqdiatrans: public Expression
  {
  public:
    lrqqdiatrans(): Expression() {}
    double Regular(double const& u) const
    {
      const double expr = 2 * u * log(u) / ( 1 - u ) + 1 - u - 2 * log(1 - u);
      return 2 * CF * expr;
    }
  };

  class rlqqdiatrans: public Expression
  {
  public:
    rlqqdiatrans(): Expression() {}
    double Regular(double const& z) const
    {
      const double expr = 4 * z * log(z) / ( 1 - z ) + 1 - z - 2 * log(1 - z);
      return 2 * CF * expr;
    }
  };

  class srqqdiatrans: public Expression
  {
  public:
    srqqdiatrans(): Expression() {}
    double Regular(double const&) const { return - 4 * CF; }
  };

  class rsqqdiatrans: public Expression
  {
  public:
    rsqqdiatrans(): Expression() {}
    double Regular(double const&) const { return - 4 * CF; }
  };

  class r11qqdiatrans: public Expression
  {
  public:
    r11qqdiatrans(): Expression() {}
    double Regular(double const&) const { return 1; }
  };

  class r12qqdiatrans: public Expression
  {
  public:
    r12qqdiatrans(): Expression() {}
    double Regular(double const&) const { return 1; }
  };

  // Functions that fills in the DIS transversely polarised hard cross
  // sections on two different grids.
  DiaTransObjects InitializeDIAtrans(Grid const& gz, Grid const& gu, std::vector<double> const& Thresholds)
  {
    report("Initialising DIA transversely polarised hard cross sections...");
    Timer t;

    // Define object of the structure containing the DiaTransObjects
    DiaTransObjects DiaObj;

    // ====================================================
    // Compute DIA partonic cross sections
    // ====================================================
    // LO contribution
    const Operator odeltaz{gz, delta{}};
    const Operator odeltau{gu, delta{}};

    DiaObj.dC0qq.AddTerm({1, odeltaz, odeltau});

    // NLO contributions
    const Operator oD0z{gz, DDn{0}};
    const Operator oD1z{gz, DDn{1}};
    const Operator oD0u{gu, DDn{0}};
    const Operator oD1u{gu, DDn{1}};

    const double LLqq = 2 * CF * ( Pi2 - 8 );
    const double LSqq = 4 * CF;
    const double SLqq = 4 * CF;
    const double SSqq = 4 * CF;

    const double K1qq = 4 * CF;

    const Operator orlqqz{gz,  rlqqdiatrans{}};
    const Operator orsqqz{gz,  rsqqdiatrans{}};
    const Operator or11qqz{gz, r11qqdiatrans{}};
    const Operator olrqqu{gu,  lrqqdiatrans{}};
    const Operator osrqqu{gu,  srqqdiatrans{}};
    const Operator or12qqu{gu, r12qqdiatrans{}};

    DiaObj.dC1qq.AddTerm({LLqq, odeltaz, odeltau});
    DiaObj.dC1qq.AddTerm({LSqq, odeltaz, oD1u   });
    DiaObj.dC1qq.AddTerm({1,    odeltaz, olrqqu });
    DiaObj.dC1qq.AddTerm({SLqq, oD1z,    odeltau});
    DiaObj.dC1qq.AddTerm({SSqq, oD0z,    oD0u   });
    DiaObj.dC1qq.AddTerm({1,    oD0z,    osrqqu });
    DiaObj.dC1qq.AddTerm({1,    orlqqz,  odeltau});
    DiaObj.dC1qq.AddTerm({1,    orsqqz,  oD0u   });
    DiaObj.dC1qq.AddTerm({K1qq, or11qqz, or12qqu});

    return DiaObj;
  }

  // Functions that fills in the DIA transversely polarised hard cross
  // sections on one single grid and exchanges the last defaulted
  // arguments.
  DiaTransObjects InitializeDIAtrans(Grid const& gz, std::vector<double> const& Thresholds)
  {
    return InitializeDIAtrans(gz, gz, Thresholds);
  }
}
