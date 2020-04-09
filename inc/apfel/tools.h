//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include <vector>

namespace apfel
{
  /**
   * @name Tools
   * Collection of useful tools.
   */
  ///@{
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
   * @param sel: the flavour selector (default: -1, i.e. all flavours are computed)
   * @return the std::vector of the electroweak charges
   */
  std::vector<double> ElectroWeakCharges(double const& Q, bool const& virt, int const& sel = -1);

  /**
   * @brief Utility function for the computation of the electroweak
   * charges for Drell-Yan in narrow-width appriximatiob
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
   * @brief Function used for the recursive expansion \prod_{n=1}^k
   * (x-r_n) = \sum_{m=0}^k (-1)^m p(m) x^m used in turn for the
   * integration on the grid.
   * @param r: first input vector
   * @param a: second input vector
   */
  std::vector<double> VectorComposition(std::vector<double> const& r, std::vector<double> const& a);
  ///@}
}
