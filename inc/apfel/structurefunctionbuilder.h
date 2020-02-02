//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/grid.h"
#include "apfel/dglap.h"
#include "apfel/observable.h"
#include "apfel/disbasis.h"

namespace apfel
{
  /**
   * @brief Structure that contains all the precomputed quantities
   * needed to compute the DIS structure functions, i.e. the
   * perturbative coefficients of the coefficient functions for
   * F<SUB>2</SUB>, F<SUB>L</SUB>, and xF<SUB>3</SUB>.
   */
  struct StructureFunctionObjects
  {
    std::vector<int>             skip;
    std::map<int,ConvolutionMap> ConvBasis;
    std::map<int,Set<Operator>>  C0;
    std::map<int,Set<Operator>>  C1;
    std::map<int,Set<Operator>>  C2;
  };

  /**
   * @name DIS structure function object initializers
   * Collection of functions that initialise StructureFunctionObjects
   * structure for the different kinds of structure functions
   * available.
   */
  ///@{
  /**
   * @brief The InitializeF2NCObjectsZM precomputes the perturbative
   * coefficients of coefficient functions for NC F2 in the ZM scheme
   * and store them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF2NCObjectsZM(Grid                const& g,
                                                                                                             std::vector<double> const& Thresholds,
                                                                                                             double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeFLNCObjectsZM precomputes the perturbative
   * coefficients of coefficient functions for NC FL in the ZM scheme
   * and store them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeFLNCObjectsZM(Grid                const& g,
                                                                                                             std::vector<double> const& Thresholds,
                                                                                                             double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeF3NCObjectsZM precomputes the perturbative
   * coefficients of coefficient functions for NC xF3 in the ZM scheme
   * and store them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF3NCObjectsZM(Grid                const& g,
                                                                                                             std::vector<double> const& Thresholds,
                                                                                                             double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeF2CCPlusObjectsZM precomputes the perturbative
   * coefficients of coefficient functions for combination ( F2(nu) +
   * F2(nubar) ) / 2 in the ZM scheme and store them in the
   * 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF2CCPlusObjectsZM(Grid                const& g,
                                                                                                                 std::vector<double> const& Thresholds,
                                                                                                                 double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeF2CCMinusObjectsZM precomputes the perturbative
   * coefficients of coefficient functions for combination ( F2(nu) -
   * F2(nubar) ) / 2 in the ZM scheme and store them in the
   * 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF2CCMinusObjectsZM(Grid                const& g,
                                                                                                                  std::vector<double> const& Thresholds,
                                                                                                                  double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeFLCCPlusObjectsZM precomputes the perturbative
   * coefficients of coefficient functions for combination ( FL(nu) +
   * FL(nubar) ) / 2 in the ZM scheme and store them in the
   * 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeFLCCPlusObjectsZM(Grid                const& g,
                                                                                                                 std::vector<double> const& Thresholds,
                                                                                                                 double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeFLCCMinusObjectsZM precomputes the perturbative
   * coefficients of coefficient functions for combination ( FL(nu) -
   * FL(nubar) ) / 2 in the ZM scheme and store them in the
   * 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeFLCCMinusObjectsZM(Grid                const& g,
                                                                                                                  std::vector<double> const& Thresholds,
                                                                                                                  double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeF3CCPlusObjectsZM precomputes the perturbative
   * coefficients of coefficient functions for combination ( F3(nu) +
   * F3(nubar) ) / 2 in the ZM scheme and store them in the
   * 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF3CCPlusObjectsZM(Grid                const& g,
                                                                                                                 std::vector<double> const& Thresholds,
                                                                                                                 double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeF3CCMinusObjectsZM precomputes the perturbative
   * coefficients of coefficient functions for combination ( F3(nu) -
   * F3(nubar) ) / 2 in the ZM scheme and store them in the
   * 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF3CCMinusObjectsZM(Grid                const& g,
                                                                                                                  std::vector<double> const& Thresholds,
                                                                                                                  double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeF2NCObjectsMassive precomputes the
   * perturbative coefficients of coefficient functions for
   * combination NC F2 in the massive scheme and store them in the
   * 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy quark masses
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @param nxi: the number of nodes of the grid in &xi; = Q<SUP>2</SUP>/M<SUP>2</SUP> (default: 150)
   * @param ximin: the lower bound of the grid in &xi; (default: 0.001)
   * @param ximax: the upper bound of the grid in &xi; (default: 100000)
   * @param intdeg: the interpolation degree on the grid in &xi; (default: 3)
   * @param lambda: the value of the parameter in the function ln(ln(&xi;/&Lambda;<SUP>2</SUP>)) used for the tabulation (default: 0.0005)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF2NCObjectsMassive(Grid                const& g,
                                                                                                                  std::vector<double> const& Masses,
                                                                                                                  double              const& IntEps = 1e-5,
                                                                                                                  int                 const& nxi    = 150,
                                                                                                                  double              const& ximin  = 0.001,
                                                                                                                  double              const& ximax  = 100000,
                                                                                                                  int                 const& intdeg = 3,
                                                                                                                  double              const& lambda = 0.0005);

  /**
   * @brief The InitializeFLNCObjectsMassive precomputes the
   * perturbative coefficients of coefficient functions for
   * combination NC FL in the massive scheme and store them in the
   * 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy quark masses
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @param nxi: the number of nodes of the grid in &xi; = Q<SUP>2</SUP>/M<SUP>2</SUP> (default: 150)
   * @param ximin: the lower bound of the grid in &xi; (default: 0.001)
   * @param ximax: the upper bound of the grid in &xi; (default: 100000)
   * @param intdeg: the interpolation degree on the grid in &xi; (default: 3)
   * @param lambda: the value of the parameter in the function ln(ln(&xi;/&Lambda;<SUP>2</SUP>)) used for the tabulation (default: 0.0005)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeFLNCObjectsMassive(Grid                const& g,
                                                                                                                  std::vector<double> const& Masses,
                                                                                                                  double              const& IntEps = 1e-5,
                                                                                                                  int                 const& nxi    = 150,
                                                                                                                  double              const& ximin  = 0.001,
                                                                                                                  double              const& ximax  = 100000,
                                                                                                                  int                 const& intdeg = 3,
                                                                                                                  double              const& lambda = 0.0005);

  /**
   * @brief The InitializeF2NCObjectsMassiveZero precomputes the
   * perturbative coefficients of coefficient functions for
   * combination NC F2 in the massless limit of the massive scheme and
   * store them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy quark masses
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @param nxi: the number of nodes of the grid in &xi; = Q<SUP>2</SUP>/M<SUP>2</SUP> (default: 150)
   * @param ximin: the lower bound of the grid in &xi; (default: 0.001)
   * @param ximax: the upper bound of the grid in &xi; (default: 100000)
   * @param intdeg: the interpolation degree on the grid in &xi; (default: 3)
   * @param lambda: the value of the parameter in the function ln(ln(&xi;/&Lambda;<SUP>2</SUP>)) used for the tabulation (default: 0.0005)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF2NCObjectsMassiveZero(Grid                const& g,
                                                                                                                      std::vector<double> const& Masses,
                                                                                                                      double              const& IntEps = 1e-5,
                                                                                                                      int                 const& nxi    = 150,
                                                                                                                      double              const& ximin  = 0.001,
                                                                                                                      double              const& ximax  = 100000,
                                                                                                                      int                 const& intdeg = 3,
                                                                                                                      double              const& lambda = 0.0005);

  /**
   * @brief The InitializeFLNCObjectsMassiveZero precomputes the
   * perturbative coefficients of coefficient functions for
   * combination NC FL in the massless limit of the massive scheme and
   * store them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy quark masses
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @param nxi: the number of nodes of the grid in &xi; = Q<SUP>2</SUP>/M<SUP>2</SUP> (default: 150)
   * @param ximin: the lower bound of the grid in &xi; (default: 0.001)
   * @param ximax: the upper bound of the grid in &xi; (default: 100000)
   * @param intdeg: the interpolation degree on the grid in &xi; (default: 3)
   * @param lambda: the value of the parameter in the function ln(ln(&xi;/&Lambda;<SUP>2</SUP>)) used for the tabulation (default: 0.0005)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeFLNCObjectsMassiveZero(Grid                const& g,
                                                                                                                      std::vector<double> const& Masses,
                                                                                                                      double              const& IntEps = 1e-5,
                                                                                                                      int                 const& nxi    = 150,
                                                                                                                      double              const& ximin  = 0.001,
                                                                                                                      double              const& ximax  = 100000,
                                                                                                                      int                 const& intdeg = 3,
                                                                                                                      double              const& lambda = 0.0005);
  ///@}

  /**
   * @name SIA structure function object initializers
   * Collection of functions that initialise StructureFunctionObjects
   * structure for the different kinds of structure functions
   * available.
   * @note For now only Zero-Mass structure functions up to
   * O(&alpha;<SUB>s</SUB>) are implemented.
   */
  ///@{
  /**
   * @brief The InitializeF2NCObjectsZMT precomputes the perturbative
   * coefficients of coefficient functions for NC F2 for SIA in the ZM
   * scheme and store them in the 'StructureFunctionObjects'
   * structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF2NCObjectsZMT(Grid                const& g,
                                                                                                              std::vector<double> const& Thresholds,
                                                                                                              double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeFLNCObjectsZMT precomputes the perturbative
   * coefficients of coefficient functions for NC FL for SIA in the ZM
   * scheme and store them in the 'StructureFunctionObjects'
   * structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeFLNCObjectsZMT(Grid                const& g,
                                                                                                              std::vector<double> const& Thresholds,
                                                                                                              double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeF3NCObjectsZMT precomputes the perturbative
   * coefficients of coefficient functions for NC xF3 for SIA in the
   * ZM scheme and store them in the 'StructureFunctionObjects'
   * structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeF3NCObjectsZMT(Grid                const& g,
                                                                                                              std::vector<double> const& Thresholds,
                                                                                                              double              const& IntEps = 1e-5);
  ///@}

  /**
   * @name Structure function builders
   * Collection of functions that build a map of Observable objects
   * corresponding to the different component of the structure
   * functions.
   */
  ///@{
  /**
   * @brief The BuildStructureFunctions function constructs a map of
   * "Observable" objects.
   * @param FObj: the StructureFunctionObjects-valued for the structure function objects
   * @param InDistFunc: the distribution to be convoluted with as a map<int,double>-valued function of x and Q
   * @param PerturbativeOrder: the perturbative order
   * @param Alphas: the strong coupling function
   * @param Couplings: the vector-valued function of (non-QCD) couplings
   * @return A map of "Observable" objects, one for number of active flavours
   */
  std::map<int,Observable<>> BuildStructureFunctions(std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> const& FObj,
                                                     std::function<std::map<int,double>(double const&, double const&)>                  const& InDistFunc,
                                                     int                                                                                const& PerturbativeOrder,
                                                     std::function<double(double const&)>                                               const& Alphas,
                                                     std::function<std::vector<double>(double const&)>                                  const& Couplings);

  /**
   * @brief The BuildStructureFunctions function constructs a map of
   * "Observable" objects.
   * @param FObj: the StructureFunctionObjects-valued for the structure function objects
   * @param InDistFunc: the distribution to be convoluted with as a double-valued function of i, x, and Q
   * @param PerturbativeOrder: the perturbative order
   * @param Alphas: the strong coupling function
   * @param Couplings: the vector-valued function of (non-QCD) couplings
   * @return A map of "Observable" objects, one for number of active flavours
   */
  std::map<int,Observable<>> BuildStructureFunctions(std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> const& FObj,
                                                     std::function<double(int const&, double const&, double const&)>                    const& InDistFunc,
                                                     int                                                                                const& PerturbativeOrder,
                                                     std::function<double(double const&)>                                               const& Alphas,
                                                     std::function<std::vector<double>(double const&)>                                  const& Couplings);

  /**
   * @brief The BuildStructureFunctions function constructs an
   * "Observable" object.

   * @param FObjQ: the StructureFunctionObjects at the scale Q
   * @param InDistFuncQ: the distribution to be convoluted with at the scale Q as a map<int,Distribution>
   * @param PerturbativeOrder: the perturbative order
   * @param AlphasQ: the strong coupling at the scale Q
   * @param k: the observable index
   * @return A "Distribution" object
   */
  Distribution BuildStructureFunctions(StructureFunctionObjects   const& FObjQ,
                                       std::map<int,Distribution> const& InDistFuncQ,
                                       int                        const& PerturbativeOrder,
                                       double                     const& AlphasQ,
                                       int                        const& k);

  /**
   * @brief The BuildStructureFunctions function constructs a map of
   * "Observable" objects.
   * @param FObjQ: the StructureFunctionObjects at the scale Q
   * @param InDistFuncQ: the distribution to be convoluted with at the scale Q as a map<int,Distribution>
   * @param PerturbativeOrder: the perturbative order
   * @param AlphasQ: the strong coupling at the scale Q
   * @return A map of "Distribution" objects, one for number of active flavours
   */
  std::map<int,Distribution> BuildStructureFunctions(StructureFunctionObjects   const& FObjQ,
                                                     std::map<int,Distribution> const& InDistFuncQ,
                                                     int                        const& PerturbativeOrder,
                                                     double                     const& AlphasQ);
  ///@}
}
