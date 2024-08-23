//
// APFEL++ 2017
//
// Author:Valerio Bertone:valerio.bertone@cern.ch
//

#pragma once

#include "apfel/doubleexpression.h"

namespace apfel
{

  /**
   * @brief Description of the class
   */
  class C1LQ2Q: public DoubleExpression
  {
  public:
    C1LQ2Q();
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class C1LQ2G: public DoubleExpression
  {
  public:
    C1LQ2G();
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class C1LG2Q: public DoubleExpression
  {
  public:
    C1LG2Q();
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class C1TQ2Q: public DoubleExpression
  {
  public:
    C1TQ2Q();
    double LocalLocal(double const& x, double const& z) const override;
    double LocalSingular(double const& x, double const& z) const override;
    double LocalRegular(double const& x, double const& z) const override;
    double SingularLocal(double const& x, double const& z) const override;
    double SingularSingular(double const& x, double const& z) const override;
    double SingularRegular(double const& x, double const& z) const override;
    double RegularLocal(double const& x, double const& z) const override;
    double RegularSingular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class C1TQ2G: public DoubleExpression
  {
  public:
    C1TQ2G();
    double LocalRegular(double const& x, double const& z) const override;
    double SingularRegular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class C1TG2Q: public DoubleExpression
  {
  public:
    C1TG2Q();
    double RegularLocal(double const& x, double const& z) const override;
    double RegularSingular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class C2LQ2G: public DoubleExpression
  {
  public:
    C2LQ2G();
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class C2LG2G: public DoubleExpression
  {
  public:
    C2LG2G();
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class C2LG2Q: public DoubleExpression
  {
  public:
    C2LG2Q();
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class C2LQ2QNS: public DoubleExpression
  {
  public:
    C2LQ2QNS(int const& nf);
    double RegularRegular(double const& x, double const& z) const override;
  private:
    int const _nf;
  };

  /**
   * @brief Description of the class
   */
  class C2LQ2QPS: public DoubleExpression
  {
  public:
    C2LQ2QPS();
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class C2LQ2QB: public DoubleExpression
  {
  public:
    C2LQ2QB();
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class C2LQ2QP1: public DoubleExpression
  {
  public:
    C2LQ2QP1();
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class C2LQ2QP2: public DoubleExpression
  {
  public:
    C2LQ2QP2();
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class C2LQ2QP3: public DoubleExpression
  {
  public:
    C2LQ2QP3();
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class C2TQ2G: public DoubleExpression
  {
  public:
    C2TQ2G();
    double LocalRegular(double const& x, double const& z) const override;
    double SingularRegular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class C2TG2G: public DoubleExpression
  {
  public:
    C2TG2G();
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class C2TG2Q: public DoubleExpression
  {
  public:
    C2TG2Q();
    double RegularLocal(double const& x, double const& z) const override;
    double RegularSingular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class C2TQ2QNS: public DoubleExpression
  {
  public:
    C2TQ2QNS(int const& nf);
    double LocalLocal(double const& x, double const& z) const override;
    double LocalSingular(double const& x, double const& z) const override;
    double LocalRegular(double const& x, double const& z) const override;
    double SingularLocal(double const& x, double const& z) const override;
    double SingularSingular(double const& x, double const& z) const override;
    double SingularRegular(double const& x, double const& z) const override;
    double RegularLocal(double const& x, double const& z) const override;
    double RegularSingular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  private:
    int const _nf;
  };

  /**
   * @brief Description of the class
   */
  class C2TQ2QPS: public DoubleExpression
  {
  public:
    C2TQ2QPS();
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class C2TQ2QB: public DoubleExpression
  {
  public:
    C2TQ2QB();
    double LocalRegular(double const& x, double const& z) const override;
    double SingularRegular(double const& x, double const& z) const override;
    double RegularLocal(double const& x, double const& z) const override;
    double RegularSingular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class C2TQ2QP1: public DoubleExpression
  {
  public:
    C2TQ2QP1();
    double LocalRegular(double const& x, double const& z) const override;
    double SingularRegular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class C2TQ2QP2: public DoubleExpression
  {
  public:
    C2TQ2QP2();
    double RegularLocal(double const& x, double const& z) const override;
    double RegularSingular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class C2TQ2QP3: public DoubleExpression
  {
  public:
    C2TQ2QP3();
    double RegularRegular(double const& x, double const& z) const override;
  };
}
