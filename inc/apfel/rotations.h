//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/distribution.h"
#include "apfel/operator.h"

namespace apfel
{
  /**
   * @name Basis rotations
   * Collection of functions to rotate distributions from the QCD
   * evolution basis to the phsycal basis and back. Specifically, the
   * QCD evolution basis is:
   *
   * {g, &Sigma;, V, T<SUB>3</SUB>, V<SUB>3</SUB>, T<SUB>8</SUB>,
   * V<SUB>8</SUB>, T<SUB>15</SUB>, V<SUB>15</SUB>, T<SUB>24</SUB>,
   * V<SUB>24</SUB>, T<SUB>35</SUB>, T<SUB>34</SUB>}
   *
   * while the physical basis is:
   *
   * { <span style="text-decoration: overline">t</span>, <span
   * style="text-decoration: overline">b</span>, <span
   * style="text-decoration: overline">c</span>, <span
   * style="text-decoration: overline">s</span>, <span
   * style="text-decoration: overline">u</span>, <span
   * style="text-decoration: overline">d</span>, g, d, u, s, c, b, t}
   */
  ///@{
  /**
   * @brief Rotation matrix from the QCD evolution to the physical
   * basis.
   */
  const double RotQCDEvToPhys[6][6] =
  {
    {1./12., -1./4., 1./12., 1./24.,  1./40.,  1./60.},    // d
    {1./12.,  1./4., 1./12., 1./24.,  1./40.,  1./60.},    // u
    {1./12.,     0., -1./6., 1./24.,  1./40.,  1./60.},    // s
    {1./12.,     0.,     0., -1./8.,  1./40.,  1./60.},    // c
    {1./12.,     0.,     0.,     0., -1./10.,  1./60.},    // b
    {1./12.,     0.,     0.,     0.,      0., -1./12.}     // t
  };

  /**
   * @brief Rotation from the physical to the QCD evolution basis.
   * @param InPhysMap: the map in the physical basis
   * @return The map in the QCD evolution basis
   */
  std::map<int, double> PhysToQCDEv(std::map<int, double> const& InPhysMap);

  /**
   * @brief Rotation from the QCD evolution to the physical basis.
   * @param QCDEvMap: The map in the QCD evolution basis
   * @return the map in the physical basis
   */
  std::map<int, double> QCDEvToPhys(std::map<int, double> const& QCDEvMap);

  /**
   * @brief Rotation from the QCD evolution to the physical basis of a
   * map of distributions.
   * @param QCDEvMap: The map of distributions in the QCD evolution basis
   * @return the map of distributions in the physical basis
   */
  std::map<int, Distribution> QCDEvToPhys(std::map<int, Distribution> const& QCDEvMap);

  /**
   * @brief Rotation from the QCD evolution to the physical basis of a
   * map of operators.
   * @param QCDEvMap: The map of operatoirs in the QCD evolution basis
   * @return the map of operators in the physical basis
   */
  std::map<int, Operator> QCDEvToPhys(std::map<int, Operator> const& QCDEvMap);
  ///@}
}
