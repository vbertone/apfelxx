//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
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
   * basis. Inverse of RotPhysToQCDEv.
   */
  const double RotQCDEvToPhys[6][6] =
  {
    {1./6., -1./2.,  1./6., 1./12., 1./20., 1./30.},    // T/V --> d^{+-}
    {1./6.,  1./2.,  1./6., 1./12., 1./20., 1./30.},    //         u^{+-}
    {1./6.,     0., -1./3., 1./12., 1./20., 1./30.},    //         s^{+-}
    {1./6.,     0.,     0., -1./4., 1./20., 1./30.},    //         c^{+-}
    {1./6.,     0.,     0.,     0., -1./5., 1./30.},    //         b^{+-}
    {1./6.,     0.,     0.,     0.,     0., -1./6.}     //         t^{+-}
  };

  /**
   * @brief Rotation matrix from the QCD evolution to the physical
   * basis for the full bases. Inverse of RotPhysToQCDEvFull.
   */
  const double RotQCDEvToPhysFull[13][13] =
  {
    { 0.,  1./12., -1./12.,     0.,     0.,      0.,      0.,      0.,      0.,     0.,       0., -1./12.,  1./12.},
    { 0.,  1./12., -1./12.,     0.,     0.,      0.,      0.,      0.,      0., -1./10.,  1./10.,  1./60., -1./60.},
    { 0.,  1./12., -1./12.,     0.,     0.,      0.,      0.,  -1./8.,   1./8.,  1./40., -1./40.,  1./60., -1./60.},
    { 0.,  1./12., -1./12.,     0.,     0.,  -1./6.,   1./6.,  1./24., -1./24.,  1./40., -1./40.,  1./60., -1./60.},
    { 0.,  1./12., -1./12.,  1./4., -1./4.,  1./12., -1./12.,  1./24., -1./24.,  1./40., -1./40.,  1./60., -1./60.},
    { 0.,  1./12., -1./12., -1./4.,  1./4.,  1./12., -1./12.,  1./24., -1./24.,  1./40., -1./40.,  1./60., -1./60.},
    { 1.,      0.,      0.,     0.,     0.,      0.,      0.,      0.,      0.,     0.,       0.,      0.,      0.},
    { 0.,  1./12.,  1./12., -1./4., -1./4.,  1./12.,  1./12.,  1./24.,  1./24.,  1./40.,  1./40.,  1./60.,  1./60.},
    { 0.,  1./12.,  1./12.,  1./4.,  1./4.,  1./12.,  1./12.,  1./24.,  1./24.,  1./40.,  1./40.,  1./60.,  1./60.},
    { 0.,  1./12.,  1./12.,     0.,     0.,  -1./6.,  -1./6.,  1./24.,  1./24.,  1./40.,  1./40.,  1./60.,  1./60.},
    { 0.,  1./12.,  1./12.,     0.,     0.,      0.,      0.,  -1./8.,  -1./8.,  1./40.,  1./40.,  1./60.,  1./60.},
    { 0.,  1./12.,  1./12.,     0.,     0.,      0.,      0.,      0.,      0., -1./10., -1./10.,  1./60.,  1./60.},
    { 0.,  1./12.,  1./12.,     0.,     0.,      0.,      0.,      0.,      0.,     0.,       0., -1./12., -1./12.}
  };

  /**
   * @brief Rotation matrix from the physical to the QCD evolution
   * basis. Inverse of RotQCDEvToPhys.
   */
  const double RotPhysToQCDEv[6][6] =
  {
    {  1.,  1.,  1.,  1.,  1.,  1.},    // q^{+-} --> Sigma/V
    { -1.,  1.,  0.,  0.,  0.,  0.},    //            T3/V3
    {  1.,  1., -2.,  0.,  0.,  0.},    //            T8/V8
    {  1.,  1.,  1., -3.,  0.,  0.},    //            T15/V15
    {  1.,  1.,  1.,  1., -4.,  0.},    //            T24/V24
    {  1.,  1.,  1.,  1.,  1., -5.}     //            T35/V35
  };

  /**
   * @brief Rotation matrix from the physical to the QCD evolution
   * basis for the full bases. Inverse of RotQCDEvToPhysFull.
   */
  const double RotPhysToQCDEvFull[13][13] =
  {
    // tb   bb   cb   sb   ub   db   g    d    u    s    c    b    t
    {  0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.},   // g
    {  1.,  1.,  1.,  1.,  1.,  1.,  0.,  1.,  1.,  1.,  1.,  1.,  1.},   // Sigma
    { -1., -1., -1., -1., -1., -1.,  0.,  1.,  1.,  1.,  1.,  1.,  1.},   // V
    {  0.,  0.,  0.,  0.,  1., -1.,  0., -1.,  1.,  0.,  0.,  0.,  0.},   // T3
    {  0.,  0.,  0.,  0., -1.,  1.,  0., -1.,  1.,  0.,  0.,  0.,  0.},   // V3
    {  0.,  0.,  0., -2.,  1.,  1.,  0.,  1.,  1., -2.,  0.,  0.,  0.},   // T8
    {  0.,  0.,  0.,  2., -1., -1.,  0.,  1.,  1., -2.,  0.,  0.,  0.},   // V8
    {  0.,  0., -3.,  1.,  1.,  1.,  0.,  1.,  1.,  1., -3.,  0.,  0.},   // T15
    {  0.,  0.,  3., -1., -1., -1.,  0.,  1.,  1.,  1., -3.,  0.,  0.},   // V15
    {  0., -4.,  1.,  1.,  1.,  1.,  0.,  1.,  1.,  1.,  1., -4.,  0.},   // T24
    {  0.,  4., -1., -1., -1., -1.,  0.,  1.,  1.,  1.,  1., -4.,  0.},   // V24
    { -5.,  1.,  1.,  1.,  1.,  1.,  0.,  1.,  1.,  1.,  1.,  1., -5.},   // T35
    {  5., -1., -1., -1., -1., -1.,  0.,  1.,  1.,  1.,  1.,  1., -5.}    // V35
  };

  /**
   * @brief Rotation matrix from the physical basis to the basis of
   * q_{+/-} = q +/- qbar.
   */
  const double RotPhysToPlusMinus[13][13] =
  {
    // tb   bb   cb   sb   ub   db   g    d    u    s    c    b    t
    {  0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.},   // g
    {  0.,  0.,  0.,  0.,  0.,  1.,  0.,  1.,  0.,  0.,  0.,  0.,  0.},   // d+
    {  0.,  0.,  0.,  0.,  0., -1.,  0.,  1.,  0.,  0.,  0.,  0.,  0.},   // d-
    {  0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.},   // u+
    {  0.,  0.,  0.,  0., -1.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.},   // u-
    {  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.},   // s+
    {  0.,  0.,  0., -1.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.},   // s-
    {  0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,  0.},   // c+
    {  0.,  0., -1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,  0.},   // c-
    {  0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.},   // b+
    {  0., -1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.},   // b-
    {  1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.},   // t+
    { -1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.}    // t-
  };

  /**
   * @brief Rotation matrix from the basis of q_{+/-} = q +/- qbar to
   * the physical basis. Inverse of RotPhysToPlusMinus.
   */
  const double RotPlusMinusToPhys[13][13] =
  {
    //  g     d+    d-    u+    u-    s+    s-    c+    c-    b+    b-    t+    t-
    {   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  0.5, -0.5},  // tb
    {   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  0.5, -0.5,   0.,   0.},  // bb
    {   0.,   0.,   0.,   0.,   0.,   0.,   0.,  0.5, -0.5,   0.,   0.,   0.,   0.},  // cb
    {   0.,   0.,   0.,   0.,   0.,  0.5, -0.5,   0.,   0.,   0.,   0.,   0.,   0.},  // sb
    {   0.,   0.,   0.,  0.5, -0.5,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.},  // ub
    {   0.,  0.5, -0.5,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.},  // db
    {   1.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.},  // g
    {   0.,  0.5,  0.5,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.},  // d
    {   0.,   0.,   0.,  0.5, 0.5,    0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.},  // u
    {   0.,   0.,   0.,   0.,   0.,  0.5,  0.5,   0.,   0.,   0.,   0.,   0.,   0.},  // s
    {   0.,   0.,   0.,   0.,   0.,   0.,   0.,  0.5,  0.5,   0.,   0.,   0.,   0.},  // c
    {   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  0.5,  0.5,   0.,   0.},  // b
    {   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  0.5,  0.5}   // t
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

  /**
   * @brief Rotation from the physical to the PlusMinus basis.
   * @param InPhysMap: the map in the physical basis
   * @return The map in the PlusMinus basis
   */
  std::map<int, double> PhysToPlusMinus(std::map<int, double> const& InPhysMap);

  /**
   * @brief Rotation from to the PlusMinus basis to the physical basis.
   * @param PlusMinusMap: the map in the PlusMinus basis
   * @return The map in the physical basis
   */
  std::map<int, double> PlusMinusToPhys(std::map<int, double> const& PlusMinusMap);

  /**
   * @brief Rotation from to the PlusMinus basis to the physical basis
   * of a map of distributions.
   * @param PlusMinusMap: the map in the PlusMinus basis
   * @return The map of distributions in the physical basis
   */
  std::map<int, Distribution> PlusMinusToPhys(std::map<int, Distribution> const& PlusMinusMap);
  ///@}
}
