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
  class DC1Q2Q: public DoubleExpression
  {
  public:
    DC1Q2Q();
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
  class DC1Q2G: public DoubleExpression
  {
  public:
    DC1Q2G();
    double LocalRegular(double const& x, double const& z) const override;
    double SingularRegular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class DC1G2Q: public DoubleExpression
  {
  public:
    DC1G2Q();
    double RegularLocal(double const& x, double const& z) const override;
    double RegularSingular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class DC2Q2G: public DoubleExpression
  {
  public:
    DC2Q2G(int const& nf);
    double LocalRegular(double const& x, double const& z) const override;
    double SingularRegular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  private:
    int const _nf;
  };

  /**
   * @brief Description of the class
   */
  class DC2G2G: public DoubleExpression
  {
  public:
    DC2G2G();
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class DC2G2Q: public DoubleExpression
  {
  public:
    DC2G2Q(int const& nf);
    double RegularLocal(double const& x, double const& z) const override;
    double RegularSingular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  private:
    int const _nf;
  };

  /**
   * @brief Description of the class
   */
  class DC2Q2QNS: public DoubleExpression
  {
  public:
    DC2Q2QNS(int const& nf);
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
  class DC2Q2QPS: public DoubleExpression
  {
  public:
    DC2Q2QPS();
    double RegularLocal(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class DC2Q2QB: public DoubleExpression
  {
  public:
    DC2Q2QB();
    double LocalRegular(double const& x, double const& z) const override;
    double SingularRegular(double const& x, double const& z) const override;
    double RegularLocal(double const& x, double const& z) const override;
    double RegularSingular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class DC2Q2QP1: public DoubleExpression
  {
  public:
    DC2Q2QP1();
    double LocalRegular(double const& x, double const& z) const override;
    double SingularRegular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class DC2Q2QP2: public DoubleExpression
  {
  public:
    DC2Q2QP2();
    double RegularLocal(double const& x, double const& z) const override;
    double RegularSingular(double const& x, double const& z) const override;
    double RegularRegular(double const& x, double const& z) const override;
  };

  /**
   * @brief Description of the class
   */
  class DC2Q2QP3: public DoubleExpression
  {
  public:
    DC2Q2QP3();
    double RegularRegular(double const& x, double const& z) const override;
  };
}
