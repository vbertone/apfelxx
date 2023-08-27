//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include <vector>
#include <string>

namespace apfel
{
  /**
   * @name Tools
   * Collection of useful tools.
   */
  ///@{
  /// Quark enumerator
  enum QuarkFlavour: int {TOTAL, DOWN, UP, STRANGE, CHARM, BOTTOM, TOP};

  /**
   * @brief Return the number of active flavours at the scale Q given
   * the (ordered) vector of thresholds.
   * @param Q: the scale
   * @param Thresholds: the vector of thresholds
   * @return number of active flavours at Q
   */
  int NF(double const& Q, std::vector<double> const& Thresholds);

  /**
   * @brief Utility function used by the heavy-quark initiated massive
   * coefficient functions.
   * @param a: first parameter
   * @param b: second parameter
   * @param c: third parameter
   * @return Triangular function
   */
  double DeltaFun(double const& a, double const& b, double const& c);

  /**
   * @brief Utility function for the computation of the electroweak
   * charges, for both time-like and space-like virtualities
   * (Reference: https://arxiv.org/pdf/hep-ph/9711387.pdf).
   * @param Q: absolute value the virtuality of the vector boson
   * @param virt: virtuality (true: time-like, false: space-like)
   * @param Comp: the flavour selector (default: TOTAL, i.e. all flavours are computed)
   * @return the std::vector of the electroweak charges
   */
  std::vector<double> ElectroWeakCharges(double const& Q, bool const& virt, int const& Comp = TOTAL);

  /**
   * @brief Utility function for the computation of the
   * parity-violating electroweak charges, for both time-like and
   * space-like virtualities.
   * @param Q: absolute value the virtuality of the vector boson
   * @param virt: virtuality (true: time-like, false: space-like)
   * @param Comp: the flavour selector (default: TOTAL, i.e. all flavours are computed)
   * @return the std::vector of the electroweak charges
   */
  std::vector<double> ParityViolatingElectroWeakCharges(double const& Q, bool const& virt, int const& Comp = TOTAL);

  /**
   * @brief Utility function for the computation of the electroweak
   * charges for Drell-Yan in narrow-width appriximation
   * @return the std::vector of the electroweak charges
   */
  std::vector<double> ElectroWeakChargesNWA();

  /**
   * @brief Utility function that concatenates and sort the input
   * vectors.
   * @param v1: first vector
   * @param v2: second vector
   * @return a std::vector containing the sorted entries of 'v1' and 'v2'
   */
  std::vector<double> ConcatenateAndSortVectors(std::vector<double> const& v1, std::vector<double> const& v2);

  /**
   * @brief Template utility function that returns the even entries of
   * the input vector.
   * @param v: input vector
   * @return a std::vector containing the even entries of the input vector.
   */
  template<typename T>
  std::vector<T> EvenVector(std::vector<T> const& v);

  /**
   * @brief Template utility function that returns the odd entries of
   * the input vector.
   * @param v: input vector
   * @return a std::vector containing the odd entries of the input vector.
   */
  template<typename T>
  std::vector<T> OddVector(std::vector<T> const& v);

  /**
   * @brief Absolute value of the object T. In the case of a
   * Distribution, this is computed like the squared mean average of
   * the entries of the joint grid. In the case of a set of
   * distributions, the minimum dabs over the distributions is
   * returned.
   * @param d: input object
   * @return the absolute value
   */
  template<typename T>
  double dabs(T const& d);

  /**
   * @brief Function that computes the coefficients of the expansion
   * of a product of n binomials with zero's in r.
   * @param r: input vector of zero's
   */
  std::vector<double> ProductExpansion(std::vector<double> const& r);

  /**
   * @brief Factorial of an integer
   * @param n: input integer
   */
  int factorial(int const& n);

  /**
   * @brief Function that computes the total cross section in a
   * electron-positron annihilation process.
   * @param PerturbativeOrder: perturbative order of the computation
   * @param Q: vector-boson invariant mass
   * @param AlphaQCD: value of the strong coupling at Q
   * @param AlphaQED: value of the electromagnetic coupling at Q
   * @param Thresholds: heavy-quark thresholds
   * @param Comp: component of the cross section, e.g. charm, bottom, etc. (default = TOTAL)
   * @param NoCharges: whether to exclude the sum over the charge of the active flavours (default = false)
   * @return the total cross section (in nbarn)
   * @note The QCD corrections to the total cross section in a
   * electron-positron annihilation process are taken from
   * Phys.Lett. B259 (1991) 144–150.
   */
  double GetSIATotalCrossSection(int                 const& PerturbativeOrder,
                                 double              const& Q,
                                 double              const& AlphaQCD,
                                 double              const& AlphaQED,
                                 std::vector<double> const& Thresholds,
                                 QuarkFlavour        const& Comp = TOTAL,
                                 bool                const& NoCharges = false);
  ///@}
}
