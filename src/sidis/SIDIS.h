/*
  Author: Valerio Bertone
 */

#include "apfel/expression.h"
#include "apfel/doubleobject.h"
#include "apfel/operator.h"
#include "apfel/timer.h"

using namespace std;
using namespace apfel;

// SIDIS hard cross sections.
DoubleObject<Operator> C20qq;
DoubleObject<Operator> C21qq;
DoubleObject<Operator> C21gq;
DoubleObject<Operator> C21qg;

DoubleObject<Operator> CL1qq;
DoubleObject<Operator> CL1gq;
DoubleObject<Operator> CL1qg;

// Expressions needed for the computation of the SIDIS cross sections.
// F2
class delta: public Expression
{
public:
  delta(): Expression() {}
  double Local(double const&) const { return 1; }
};

class s0: public Expression
{
public:
  s0(): Expression() {}
  double Singular(double const& x) const { return 1 / ( 1 - x ); }
  double Local(double const& x) const { return log( 1 - x ); }
};

class s1: public Expression
{
public:
  s1(): Expression() {}
  double Singular(double const& x) const { return log( 1 - x ) / ( 1 - x ); }
  double Local(double const& x) const { double l = log( 1 - x ); return l * l / 2; }
};

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

// Functions that fills in the SIDIS hard cross sections.
void InitializeSIDIS(Grid const& g)
{
  cout << "Initializing SIDIS hard cross sections... ";
  Timer t;
  t.start();

  // ====================================================
  // Compute SIDIS partonic cross sections.
  // Expressions taken from Appendix C of hep-ph/9711387.
  // ====================================================
  // LO contribution.
  const Operator odelta{g, delta{}};

  C20qq.AddTerm({1, odelta, odelta});

  // NLO contributions
  // F2
  const Operator os0{g, s0{}};
  const Operator os1{g, s1{}};

  const double LLqq = - 16 * CF;
  const double LSqq = 4 * CF;
  const double SLqq = 4 * CF;
  const double SSqq = 4 * CF;
  const double K1qq = 4 * CF;
  const double K2qq = 12 * CF;

  const Operator olrqq{g, lrqq{}};
  const Operator osrqq{g, srqq{}};
  const Operator orlqq{g, rlqq{}};
  const Operator orsqq{g, rsqq{}};
  const Operator or11qq{g, r11qq{}};
  const Operator or12qq{g, r12qq{}};
  const Operator or21qq{g, r21qq{}};
  const Operator or22qq{g, r22qq{}};

  C21qq.AddTerm({LLqq, odelta, odelta}); //1
  C21qq.AddTerm({LSqq, odelta, os1   }); //2
  C21qq.AddTerm({1   , odelta, olrqq }); //3
  C21qq.AddTerm({SLqq, os1   , odelta}); //4
  C21qq.AddTerm({SSqq, os0   , os0   }); //5
  C21qq.AddTerm({1   , os0   , osrqq }); //6
  C21qq.AddTerm({1   , orlqq , odelta}); //7
  C21qq.AddTerm({1   , orsqq , os0   }); //8
  C21qq.AddTerm({K1qq, or11qq, or12qq}); //9
  C21qq.AddTerm({K2qq, or21qq, or22qq}); //10

  const double K1gq = 4 * CF;
  const double K2gq = - 12 * CF;
  const double K3gq = - 2 * CF;

  const Operator olrgq{g, lrgq{}};
  const Operator osrgq{g, srgq{}};
  const Operator or11gq{g, r11gq{}};
  const Operator or12gq{g, r12gq{}};
  const Operator or21gq{g, r21gq{}};
  const Operator or22gq{g, r22gq{}};
  const Operator or31gq{g, r31gq{}};
  const Operator or32gq{g, r32gq{}};

  C21gq.AddTerm({1   , odelta, olrgq }); //1
  C21gq.AddTerm({1   , os0   , osrgq }); //2
  C21gq.AddTerm({K1gq, or11gq, or12gq}); //3
  C21gq.AddTerm({K2gq, or21gq, or22gq}); //4
  C21gq.AddTerm({K3gq, or31gq, or32gq}); //5

  const double K1qg = 2;
  const double K2qg = 1;

  const Operator orlqg{g, rlqg{}};
  const Operator orsqg{g, rsqg{}};
  const Operator or11qg{g, r11qg{}};
  const Operator or12qg{g, r12qg{}};
  const Operator or21qg{g, r21qg{}};
  const Operator or22qg{g, r22qg{}};

  C21qg.AddTerm({1   , orlqg , odelta}); //1
  C21qg.AddTerm({1   , orsqg , os0   }); //2
  C21qg.AddTerm({K1qg, or11qg, or12qg}); //3
  C21qg.AddTerm({K2qg, or21qg, or22qg}); //4

  // FL
  const double K1Lqq = 8 * CF;

  const Operator or11Lqq{g, r11Lqq{}};
  const Operator or12Lqq{g, r12Lqq{}};

  CL1qq.AddTerm({K1Lqq, or11Lqq, or12Lqq});

  const double K1Lgq = 8 * CF;

  const Operator or11Lgq{g, r11Lgq{}};
  const Operator or12Lgq{g, r12Lgq{}};

  CL1gq.AddTerm({K1Lgq, or11Lgq, or12Lgq});

  const double K1Lqg = 8;

  const Operator or11Lqg{g, r11Lqg{}};
  const Operator or12Lqg{g, r12Lqg{}};

  CL1qg.AddTerm({K1Lqg, or11Lqg, or12Lqg});
  t.stop();
}
