//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/expression.h"
#include "apfel/doubleobject.h"
#include "apfel/operator.h"
#include "apfel/timer.h"

// SIDIS hard cross sections.
apfel::DoubleObject<apfel::Operator> C20qq;
apfel::DoubleObject<apfel::Operator> C21qq;
apfel::DoubleObject<apfel::Operator> C21gq;
apfel::DoubleObject<apfel::Operator> C21qg;

apfel::DoubleObject<apfel::Operator> CL1qq;
apfel::DoubleObject<apfel::Operator> CL1gq;
apfel::DoubleObject<apfel::Operator> CL1qg;

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
  double Local(double const& x) const { double l = log( 1 - x ); return l * l / 2; }
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

// Functions that fills in the SIDIS hard cross sections.
void InitializeSIDIS(apfel::Grid const& g)
{
  std::cout << "Initializing SIDIS hard cross sections... ";
  apfel::Timer t;
  t.start();

  // ====================================================
  // Compute SIDIS partonic cross sections.
  // Expressions taken from Appendix C of hep-ph/9711387.
  // ====================================================
  // LO contribution.
  const apfel::Operator odelta{g, delta{}};

  C20qq.AddTerm({1, odelta, odelta});

  // NLO contributions
  // F2
  const apfel::Operator os0{g, s0{}};
  const apfel::Operator os1{g, s1{}};

  const double LLqq = - 16 * apfel::CF;
  const double LSqq = 4 * apfel::CF;
  const double SLqq = 4 * apfel::CF;
  const double SSqq = 4 * apfel::CF;
  const double K1qq = 4 * apfel::CF;
  const double K2qq = 12 * apfel::CF;

  const apfel::Operator olrqq{g, lrqq{}};
  const apfel::Operator osrqq{g, srqq{}};
  const apfel::Operator orlqq{g, rlqq{}};
  const apfel::Operator orsqq{g, rsqq{}};
  const apfel::Operator or11qq{g, r11qq{}};
  const apfel::Operator or12qq{g, r12qq{}};
  const apfel::Operator or21qq{g, r21qq{}};
  const apfel::Operator or22qq{g, r22qq{}};

  C21qq.AddTerm({LLqq, odelta, odelta}); //1
  C21qq.AddTerm({LSqq, odelta, os1   }); //2
  C21qq.AddTerm({1, odelta, olrqq });    //3
  C21qq.AddTerm({SLqq, os1, odelta});    //4
  C21qq.AddTerm({SSqq, os0, os0   });    //5
  C21qq.AddTerm({1, os0, osrqq });       //6
  C21qq.AddTerm({1, orlqq, odelta});     //7
  C21qq.AddTerm({1, orsqq, os0   });     //8
  C21qq.AddTerm({K1qq, or11qq, or12qq}); //9
  C21qq.AddTerm({K2qq, or21qq, or22qq}); //10

  const double K1gq = 4 * apfel::CF;
  const double K2gq = - 12 * apfel::CF;
  const double K3gq = - 2 * apfel::CF;

  const apfel::Operator olrgq{g, lrgq{}};
  const apfel::Operator osrgq{g, srgq{}};
  const apfel::Operator or11gq{g, r11gq{}};
  const apfel::Operator or12gq{g, r12gq{}};
  const apfel::Operator or21gq{g, r21gq{}};
  const apfel::Operator or22gq{g, r22gq{}};
  const apfel::Operator or31gq{g, r31gq{}};
  const apfel::Operator or32gq{g, r32gq{}};

  C21gq.AddTerm({1, odelta, olrgq });    //1
  C21gq.AddTerm({1, os0, osrgq });       //2
  C21gq.AddTerm({K1gq, or11gq, or12gq}); //3
  C21gq.AddTerm({K2gq, or21gq, or22gq}); //4
  C21gq.AddTerm({K3gq, or31gq, or32gq}); //5

  const double K1qg = 2;
  const double K2qg = 1;

  const apfel::Operator orlqg{g, rlqg{}};
  const apfel::Operator orsqg{g, rsqg{}};
  const apfel::Operator or11qg{g, r11qg{}};
  const apfel::Operator or12qg{g, r12qg{}};
  const apfel::Operator or21qg{g, r21qg{}};
  const apfel::Operator or22qg{g, r22qg{}};

  C21qg.AddTerm({1, orlqg, odelta});     //1
  C21qg.AddTerm({1, orsqg, os0   });     //2
  C21qg.AddTerm({K1qg, or11qg, or12qg}); //3
  C21qg.AddTerm({K2qg, or21qg, or22qg}); //4

  // FL
  const double K1Lqq = 8 * apfel::CF;

  const apfel::Operator or11Lqq{g, r11Lqq{}};
  const apfel::Operator or12Lqq{g, r12Lqq{}};

  CL1qq.AddTerm({K1Lqq, or11Lqq, or12Lqq});

  const double K1Lgq = 8 * apfel::CF;

  const apfel::Operator or11Lgq{g, r11Lgq{}};
  const apfel::Operator or12Lgq{g, r12Lgq{}};

  CL1gq.AddTerm({K1Lgq, or11Lgq, or12Lgq});

  const double K1Lqg = 8;

  const apfel::Operator or11Lqg{g, r11Lqg{}};
  const apfel::Operator or12Lqg{g, r12Lqg{}};

  CL1qg.AddTerm({K1Lqg, or11Lqg, or12Lqg});
  t.stop();
}
