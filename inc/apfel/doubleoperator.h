//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/doubleexpression.h"
#include "apfel/grid.h"
#include "apfel/matrix.h"
#include "apfel/operator.h"
#include "apfel/doubledistribution.h"
#include "apfel/distributionoperator.h"
#include "apfel/operatordistribution.h"
#include "apfel/config.h"

#include <iosfwd>
#include <memory>

namespace apfel
{
  /**
   * @brief This class defines the basic object DoubleOperator which
   * essentially contains the convolution between a "DoubleExpression"
   * object and two sets of interpolanting functions defined on two
   * independent grids. This is the two-dimensional generalisation of
   * the object "Operator".
   * @note For the moment, only the inclusive massless case can be
   * handled within this class.
   */
  class DoubleOperator
  {
  public:
    DoubleOperator() = delete;
    DoubleOperator(DoubleOperator const&) = default;
    DoubleOperator(DoubleOperator&&)      = default;

    /**
     * @brief The DoubleOperator constructor.
     * @param gr1: the Grid object for the first variable
     * @param gr2: the Grid object for the second variable
     * @param dexpr: the double expression to be convoluted
     * @param eps: relative accuracy of the numerical integrations (default: 10<SUP>-3</SUP>)
     */
    DoubleOperator(Grid const& gr1, Grid const& gr2, DoubleExpression const& dexpr, double const& eps = 1e-3);

    /**
     * @brief The DoubleOperator constructor.
     * @param O1: the first single operator
     * @param O2: the second single operator
     */
    DoubleOperator(Operator const& O1, Operator const& O2);

    /**
     * @brief The DoubleOperator constructor.
     * @param DObj: double object of operators
     */
    DoubleOperator(DoubleObject<Operator> const& DObj);

    /**
     * @brief Function that dumps the DoubleOperator object onto an
     * output stream in cereal portable-binary format. Only the
     * integration accuracy, the DoubleExpression name, the grid
     * descriptors (for a consistency check on reload) and the operator
     * container are stored; the grids themselves are supplied again at
     * load time.
     * @param os: the output stream the object is written to
     */
    void EmitDoubleOperatorBinary(std::ostream& os) const;

    /**
     * @brief Convenience overload that dumps the DoubleOperator object
     * to a file in cereal portable-binary format.
     * @param filename: the path of the file the object is written to
     */
    void EmitDoubleOperatorBinary(std::string const& filename) const;

    /**
     * @brief Function that reconstructs a DoubleOperator object from an
     * input stream produced by EmitDoubleOperatorBinary.
     * @param is: the input stream the object is read from
     * @param gr1: the Grid object for the first variable
     * @param gr2: the Grid object for the second variable
     * @param dexpr: the double expression, needed to compare names
     * @return the reconstructed DoubleOperator object
     */
    static DoubleOperator ReadBinary(std::istream& is, Grid const& gr1, Grid const& gr2, DoubleExpression const& dexpr);

    /**
     * @brief Convenience overload that reconstructs a DoubleOperator
     * object from a file produced by EmitDoubleOperatorBinary.
     * @param filename: the path of the file the object is read from
     * @param gr1: the Grid object for the first variable
     * @param gr2: the Grid object for the second variable
     * @param dexpr: the double expression, needed to compare names
     * @return the reconstructed DoubleOperator object
     */
    static DoubleOperator ReadBinary(std::string const& filename, Grid const& gr1, Grid const& gr2, DoubleExpression const& dexpr);

    /**
     * @brief Self-contained reader that reconstructs a DoubleOperator
     * from an input stream produced by EmitDoubleOperatorBinary without
     * any external information: the grids are rebuilt from the stored
     * descriptors and owned by the returned object. Use this when the
     * original grids are not available (e.g. reading the operator in a
     * different program). The grids can then be retrieved through
     * GetFirstGrid()/GetSecondGrid().
     * @param is: the input stream the object is read from
     * @return the reconstructed, self-contained DoubleOperator object
     */
    static DoubleOperator ReadBinary(std::istream& is);

    /**
     * @brief Convenience overload of the self-contained reader taking a
     * file path.
     * @param filename: the path of the file the object is read from
     * @return the reconstructed, self-contained DoubleOperator object
     */
    static DoubleOperator ReadBinary(std::string const& filename);

    /**
     * @brief The Operator virtual destructor.
     */
    virtual ~DoubleOperator() {}

    /**
     * @name Binary operators
     */
    ///@{
    DistributionOperator MultiplyFirstBy(Distribution const& d) const;                      //!< this times a distribution convoluted with the first variable
    OperatorDistribution MultiplySecondBy(Distribution const& d) const;                     //!< this times a distribution convoluted with the second variable
    DoubleDistribution operator *= (DoubleDistribution const& d) const;                     //!< this *= DoubleDistribution
    DoubleOperator&    operator *= (double const& s);                                       //!< this *= Scalar
    DoubleOperator&    operator /= (double const& s);                                       //!< this /= Scalar
    DoubleOperator&    operator += (DoubleOperator const& o);                               //!< this += DoubleOperator
    DoubleOperator&    operator -= (DoubleOperator const& o);                               //!< this -= DoubleOperator
    DoubleOperator&    operator *= (DoubleOperator const& o);                               //!< this *= DoubleOperator
    DoubleOperator&    operator  = (DoubleOperator const& o);                               //!< this  = DoubleOperator
    DoubleOperator&    operator *= (std::function<double(double const&, double const&)> f); //!< this *= 2D-function
    DoubleOperator&    operator *= (std::function<double(double const&)> f);                //!< this *= 1D-Function
    ///@}

    /**
     * @brief Function that returns the first Grid object associated
     * with the double operator.
     */
    Grid const& GetFirstGrid() const { return _grid1; }

    /**
     * @brief Function that returns the second Grid object associated
     * with the double operator.
     */
    Grid const& GetSecondGrid() const { return _grid2; }

    /**
     * @brief Function that returns the integration accuracy.
     */
    double const& GetIntegrationAccuracy() const { return _eps; }

    /**
     * @brief Function that returns the name of the DoubleExpression
     * object associated with the operator.
     */
    std::string const& GetDoubleExpressionName() const { return _dexprName; }

    /**
     * @brief Function that returns the DoubleOperator container.
     */
    std::vector<std::vector<matrix<matrix<double>>>> GetDoubleOperator() const { return _dOperator; }

    /**
     * @brief Function that prints the DoubleOperator object.
     */
    void Print() const { std::cout << *this << std::endl; }

  private:
    /**
     * @brief Private all-arguments constructor used by the external-grid
     * ReadBinary to initialise the const and reference members directly
     * from the deserialised data, sidestepping the lack of a default
     * constructor. The grids are external (not owned).
     */
    DoubleOperator(Grid const& gr1, Grid const& gr2, double const& eps, std::string const& dexprName,
                   std::vector<std::vector<matrix<matrix<double>>>> dOperator);

    /**
     * @brief Private all-arguments constructor used by the
     * self-contained ReadBinary: the object takes shared ownership of
     * the rebuilt grids and binds its grid references to them.
     */
    DoubleOperator(std::shared_ptr<const Grid> grid1, std::shared_ptr<const Grid> grid2, double const& eps,
                   std::string const& dexprName, std::vector<std::vector<matrix<matrix<double>>>> dOperator);

  protected:
    std::shared_ptr<const Grid>                      _grid1Owned; //!< Owns _grid1 for self-contained operators (null otherwise)
    std::shared_ptr<const Grid>                      _grid2Owned; //!< Owns _grid2 for self-contained operators (null otherwise)
    Grid                                      const& _grid1;      //!< First grid on which to compute the operator
    Grid                                      const& _grid2;      //!< Second grid on which to compute the operator
    double                                    const  _eps;        //!< Integration accuracy
    std::string                               const  _dexprName;  //!< Name of the double expression
    std::vector<std::vector<matrix<matrix<double>>>> _dOperator;  //!< DoubleOperator container

    friend std::ostream& operator << (std::ostream& os, DoubleOperator const& dop);
  };

  /**
   * @name Ternary operators
   */
  ///@{
  DoubleDistribution   operator * (DoubleOperator const& lhs, DoubleDistribution const& rhs);                 //!< DoubleOperator*DoubleDistribution
  DistributionOperator operator * (Distribution const& lhs, DoubleOperator const& rhs);                       //!< Distribution*DoubleOperator = DistributionOperator
  OperatorDistribution operator * (DoubleOperator const& lhs, Distribution const& rhs);                       //!< DoubleOperator*Distribution = OperatorDistribution
  DoubleOperator       operator * (DoubleOperator lhs, DoubleOperator const& rhs);                            //!< DoubleOperator*DoubleOperator
  DoubleOperator       operator * (double const& s, DoubleOperator rhs);                                      //!< Scalar*DoubleOperator
  DoubleOperator       operator * (DoubleOperator lhs, double const& s);                                      //!< DoubleOperator*Scalar
  DoubleOperator       operator * (std::function<double(double const&, double const)> f, DoubleOperator rhs); //!< 2D-function*DoubleOperator
  DoubleOperator       operator * (DoubleOperator lhs, std::function<double(double const&, double const)> f); //!< DoubleOperator*2D-function
  DoubleOperator       operator * (std::function<double(double const&)> f, DoubleOperator rhs);               //!< 1D-function*DoubleOperator
  DoubleOperator       operator * (DoubleOperator lhs, std::function<double(double const&)> f);               //!< DoubleOperator*1D-function
  DoubleOperator       operator / (DoubleOperator lhs, double const& s);                                      //!< DoubleOperator/Scalar
  DoubleOperator       operator + (DoubleOperator lhs, DoubleOperator const& rhs);                            //!< DoubleOperator+DoubleOperator
  DoubleOperator       operator - (DoubleOperator lhs, DoubleOperator const& rhs);                            //!< DoubleOperator-DoubleOperator
  ///@}

  /**
   * @brief Method which prints DoubleOperator with std::cout <<.
   */
  std::ostream& operator << (std::ostream& os, DoubleOperator const& dop);
}

