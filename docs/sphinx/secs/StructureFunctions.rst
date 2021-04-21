===================
Structure functions
===================

.. contents::
   :depth: 3
..

In this part the technology developed to numerically compute Melling
convolutions exploiting interpolation functions will be applied to the
computation of totally structure functions. On top of the numerical
aspects of the computation, more “theoretical” aspects will be touched
such as the computation of scale variations and of the target-mass
corrections.

Zero-mass structure functions
=============================

A generic structure function in the zero-mass (ZM) scheme, in which the
finitness of the quark mass is neglected, computed at the value
:math:`x` of the Bjorken variable and at the (absolute) value :math:`Q`
of the vector-boson virtuality is given by the following Mellin
convolution:

.. math::
   F(x,Q) = \sum_{i=g,q}x\int_x^1\frac{dy}y
     C_i\left(y,\alpha_s(Q)\right)q_i\left(\frac{x}{y},Q\right)=\sum_{i=g,q}\int_x^1dy\,
     C_i\left(y,\alpha_s(Q)\right)\frac{x}{y}q_i\left(\frac{x}{y},Q\right)\,.
   :label: eq:structfunc

\ The functions :math:`C_i`, often dubbed coefficient functions, are
perturbatively computable and thus, to :math:`N`-th order, admit the
expansion:

.. math:: C_i\left(y,\alpha_s(Q)\right) = \sum_{n=0}^{N}\left(\frac{\alpha_s(Q)}{4\pi}\right)^nC_i^{(n)}(x)\,.

The integral in Eq. :eq:`eq:structfunc` can be
computed using the numerical techniques developed in the part of this
documentation devoted to the structure of the Mellin-convolution
integrals. Specifically, introducing a grid :math:`\{x_\alpha\}` having
a set of interpolating functions of degree :math:`k`,
:math:`\{w_\alpha^{(k)}(y)\}`, associated, the structure function on the
grid point is computed as:

.. math::
   F(x_\beta,Q) =
     \sum_{i=g,q}\sum^{N_{x}}_{\alpha=0}\sum_{n=0}^{N}\left(\frac{\alpha_s(Q)}{4\pi}\right)^n\underbrace{\int_{{\rm
             max}(x_\beta,x_\beta/x_{\alpha+1})}^{{\rm
             min}(1,x_\beta/x_{\alpha-k})}dy\,
         C_i^{(n)}\left(y\right)
         w_{\alpha}^{(k)}\left(\frac{x_\beta}{y}\right)}_{\Gamma_{i,\beta\alpha}^{(n)}}\overline{q}_{\alpha}\,.
   :label: eq:stffast

with :math:`\overline{q}_\alpha=x_\alpha q(x_\alpha,Q)`. The value of
the structure function for any generic value of :math:`x` can then be
obtained through interpolation. The integrals
:math:`\Gamma_{i,\beta\alpha}^{(n)}` are efficiently computed according
to the procedure discussed in the part on the integration. Notice, that
this integrals are totally scale independent in that the :math:`Q`
dependence of the structure functions is totally driven by
:math:`\alpha_s` and the PDFs :math:`q`. This allows for a
pre-computation of these integrals that can be used to compute structure
functions at any scale :math:`Q`.

Before moving to considering the massive structure functions, it is
useful to point out that in the ZM scheme there is no need to
distinguish between charged- and neutral-current coefficient functions.
In other words, when considering the exchange of a neutral vector boson
:math:`Z/\gamma^*` or of a charged vector boson :math:`W^\pm`, the
difference in functional form of the coefficient functions :math:`C_i`
only amounts to an overall factor associated to the electroweak charge.

Massive structure functions
===========================

When computing structure functions retaining the mass of the quarks, the
convolution integrals become more complicated. In addition, contrary to
the massless case, the single precomputed intergrals
:math:`\Gamma_{i,\beta\alpha}^{(n)}`, introduced in
Eq. :eq:`eq:stffast`, will depend on the scale
:math:`Q`. This prevents a scale independent pre-computation of the
coefficient functions on :math:`x`-space grid. A possible solution to
this problem relies on interpolation. Specifically, the integrals
:math:`\Gamma_{i,\beta\alpha}^{(n)}` are tabulated over a grid in the
variable:

.. math::
   \xi = \frac{Q^2}{m_H^2}\,,
   :label: eq:scalingvar

\ where :math:`m_H` is the mass of the heavy quark under consideration,
and subsequently interpolated to obtain the structure function for any
value of :math:`Q`.

.. container::
   :name: neutral-current-coefficient-functions

   .. rubric:: Neutral Current Coefficient Functions
      :name: neutral-current-coefficient-functions

As far as the neutral-current coefficient functions are concerned,
explicit expressions beyond :math:`\mathcal{O}(\alpha_s)` are
complicated and thus not suitable for fast numerical computations.
Fortunately, the authors of Ref. (Riemersma, Smith, and Neerven 1995)
have tabulated these functions and published the corresponding
interpolating routines. It is these routines that ``APFEL++`` uses for
the :math:`\mathcal{O}(\alpha_s^2)` corrections to the massive neutral
current coefficient functions. The only exception is the pure-singlet
:math:`\mathcal{O}(\alpha_s^2)` coefficient functions (sometimes called
gluon-radiation terms) in which case the analytical expressions given in
Appendix A of Ref. (Buza et al. 1996) are used. At
:math:`\mathcal{O}(\alpha_s)` the analytic expressions are instead
used. [1]_

In the implementation of a so-called General-Mass (GM) scheme, it is
usually required to know the massless limit of the massive coefficient
functions analytically. To be more precise, this limit prescribes to set
to zero all mass suppressed contributions and retain only the
mass-independent and the logarithmically enhanced contributions. In this
case, exact expressions up to :math:`\mathcal{O}(\alpha_s^2)` have been
computed in Ref. (Buza et al. 1996) and reported in Appendix D of this
paper. These expressions are implemented in ``APFEL++``.

As a final remark, the massive coefficient functions for the neutral
current structure functions are presently known only for :math:`F_2` and
:math:`F_L`. For the parity-violating structure function :math:`F_3` the
massless coefficient functions are instead used. This is usually
acceptable because the neutral-current :math:`F_3` in the typical
kinematics covered by modern experiments is sizeable at large :math:`Q`
only where mass effects are negligible.

.. container::
   :name: charged-current-coefficient-functions

   .. rubric:: Charged Current Coefficient Functions
      :name: charged-current-coefficient-functions

We can now consider the charged-current sector in which massive
coefficient functions are know up to :math:`\mathcal{O}(\alpha_s^2)`.
However, the :math:`\mathcal{O}(\alpha_s^2)` corrections, recently
computed in Ref. (Gao 2018), are not publicly available. Therefore a
computation of charged-current structure functions at NNLO in
``APFEL++`` is currently impossible.

The structure functions associated to the heavy quark :math:`H` in the
approximation of diagonal CKM matrix [2]_ are given in terms of the
following convolutions:

.. math::
   F_1^H(x,Q,m_H)=\frac12\int_{\chi}^{1}\frac{dy}{y}\left[C_{1,q}(y,Q)s\left(\frac{\chi}{y},Q\right)+C_{1,g}(y,Q)g\left(\frac{\chi}{y},Q\right)\right]
   :label: F1

.. math::
   F_2^H(x,Q,m_H)=\chi\int_{\chi}^{1}\frac{dy}{y}\left[C_{2,q}(y,Q)s\left(\frac{\chi}{y},Q\right)+C_{2,g}(y,Q)g\left(\frac{\chi}{y},Q\right)\right]
   :label: F2

.. math::
   F_3^H(x,Q,m_H)=\int_{\chi}^{1}\frac{dy}{y}\left[C_{3,q}(y,Q)s\left(\frac{\chi}{y},Q\right)+C_{3,g}(y,Q)g\left(\frac{\chi}{y},Q\right)\right]
   :label: F3

with:

.. math::
   \chi = x\left(1+\frac{m_H^2}{Q^2}\right) =
     \frac{x}{\lambda}\,,

where:

.. math::
   \lambda = \frac{Q^2}{Q^2+m_H^2} =
     \frac{\xi}{1+\xi}\,,

with :math:`xi` given in Eq. :eq:`eq:scalingvar`.
Defining:

.. math:: F_L^H(x,Q,m_H) = F_2^H(x,Q,m_H) - 2xF_1^H(x,Q,m_H)\,,

one has that:

.. math::
   F_L^H(x,Q,m_H)=\chi\int_{\chi}^{1}\frac{dy}{y}\left[C_{L,q}(y,Q)s\left(\frac{\chi}{y},Q\right)+C_{L,g}(y,Q)g\left(\frac{\chi}{y},Q\right)\right]\,,
   :label: FL

with:

.. math::
   C_{L,q(g)}(y,Q) =
     C_{2,q(g)}(y,Q)-\lambda C_{1,q(g)}(y,Q)
   :label: clll

As usual, the coefficient functions admit a perturbative expansion that
at N\ :math:`^N`\ LO reads:

.. math::
   C_{k,q(g)}(y,Q) = \sum_{n=0}^N \left(\frac{\alpha_s(Q)}{4\pi}\right)^n
     C_{k,q(g)}^{(n)}(y,\xi)\,,\quad k = 1,2,3,L\,.

\ In the following we will truncate the expansion at :math:`N=1`, *i.e.*
at NLO.

At LO the coefficient functions read:

.. math::
   \begin{array}{l}
     \displaystyle C^{(0)}_{1,q}(x,\xi) = \delta(1-x)\,,\\
     \\ \displaystyle C^{(0)}_{2,q}(x,\xi) = \delta(1-x) \,,\\ \\
     \displaystyle C^{(0)}_{3,q}(x,\xi) = \delta(1-x) \,,\\ \\
     \displaystyle C^{(0)}_{L,q}(x,\xi) = (1-\lambda)\delta(1-x) \,,\\ \\
     \displaystyle  C_{k,g}^{(0)}(y,\xi)=0\,, \quad k=1,2,3,L\,.
   \end{array}

The :math:`\mathcal{O}(\alpha_s)` (NLO) charged-current massive
coefficient have been computed and reported in Appendix A of
Ref. (Gluck, Kretzer, and Reya 1996). However, their implementation in
`` APFEL++`` requires some manipulations. We start by defining:

.. math:: K_A=\frac{1}{\lambda}(1-\lambda)\ln(1-\lambda)\quad \mbox{and}\quad K_F =\frac{Q}{\mu_F}\,.

The explicit expressions of the NLO quark coefficient functions then
read:

.. math::
   \begin{array}{rcl}
     C^{(1)}_{1,q}&=&\displaystyle 2C_F \bigg\{
                      \bigg(-4-\frac{1}{2\lambda}-2\zeta_2-\frac{1+3\lambda}{2\lambda}K_A+\frac32
                      \ln\frac{K_F^2}{\lambda}\bigg)\delta(1-z)\\ \\ &-&\displaystyle
                                                                         \frac{(1+z^2)\ln z}{1-z} +
                                                                         \left(-\ln\frac{K_F^2}{\lambda}-2\ln(1-z)+\ln(1-\lambda
                                                                         z)\right)(1+z)+(3-z)+\frac{1}{\lambda^2}+\frac{z-1}{\lambda}\\ \\
                  &+&\displaystyle 2 \left[\frac{2\ln(1-z)-\ln(1-\lambda
                      z)}{1-z}\right]_++
                      2\left(-1+\ln\frac{K_F^2}{\lambda}\right)\left[\frac{1}{1-z}\right]_+\\
     \\ &+& \displaystyle
            \frac{\lambda-1}{\lambda^2}\left[\frac{1}{1-\lambda z}\right]_+
            +\frac{1}{2}\left[\frac{1-z}{(1-\lambda z)^2}\right]_+\bigg\}\,,
   \end{array}

.. math::
   \begin{array}{rcl}
     C^{(1)}_{2,q}&=&\displaystyle 2C_F \bigg\{
                      \bigg(-4-\frac{1}{2\lambda}-2\zeta_2-\frac{1+\lambda}{2\lambda}K_A+\frac32
                      \ln\frac{K_F^2}{\lambda}\bigg)\delta(1-z)\\ \\ &-&\displaystyle
                                                                         \frac{(1+z^2)\ln z}{1-z} +
                                                                         \left(2-\ln\frac{K_F^2}{\lambda}-2\ln(1-z)+\ln(1-\lambda
                                                                         z)\right)(1+z)+\frac{1}{\lambda}\\ \\ &+&\displaystyle 2
                                                                                                                   \left[\frac{2\ln(1-z)-\ln(1-\lambda z)}{1-z}\right]_++
                                                                                                                   2\left(-1+\ln\frac{K_F^2}{\lambda}\right)\left[\frac{1}{1-z}\right]_+\\
     \\ &+& \displaystyle
            \frac{2\lambda^2-\lambda-1}{\lambda}\left[\frac{1}{1-\lambda
            z}\right]_+ +\frac{1}{2}\left[\frac{1-z}{(1-\lambda
            z)^2}\right]_+\bigg\}\,,
   \end{array}

.. math::
   \begin{array}{rcl}
     C^{(1)}_{3,q}&=&\displaystyle 2C_F \bigg\{
                      \bigg(-4-\frac{1}{2\lambda}-2\zeta_2-\frac{1+3\lambda}{2\lambda}K_A+\frac32
                      \ln\frac{K_F^2}{\lambda}\bigg)\delta(1-z)\\ \\ &-&\displaystyle
                                                                         \frac{(1+z^2)\ln z}{1-z} +
                                                                         \left(1-\ln\frac{K_F^2}{\lambda}-2\ln(1-z)+\ln(1-\lambda
                                                                         z)\right)(1+z)+\frac{1}{\lambda}\\ \\ &+&\displaystyle 2
                                                                                                                   \left[\frac{2\ln(1-z)-\ln(1-\lambda z)}{1-z}\right]_++
                                                                                                                   2\left(-1+\ln\frac{K_F^2}{\lambda}\right)\left[\frac{1}{1-z}\right]_+\\
     \\ &+& \displaystyle \frac{\lambda-1}{\lambda}\left[\frac{1}{1-\lambda
            z}\right]_+ +\frac{1}{2}\left[\frac{1-z}{(1-\lambda
            z)^2}\right]_+\bigg\}\,,
   \end{array}

.. math::
   \begin{array}{rcl}
     C^{(1)}_{L,q} &=&\displaystyle 2C_F
                       (1-\lambda)\bigg\{
                       \bigg(-4-\frac{1}{2\lambda}-2\zeta_2-\frac{1+\lambda}{2\lambda}K_A+\frac32
                       \ln\frac{K_F^2}{\lambda}\bigg)\delta(1-z)\\ \\ &-&\displaystyle
                                                                          \frac{(1+z^2)\ln z}{1-z} +
                                                                          \left(-\ln\frac{K_F^2}{\lambda}-2\ln(1-z)+\ln(1-\lambda
                                                                          z)\right)(1+z)+3\\ \\ &+&\displaystyle 2
                                                                                                    \left[\frac{2\ln(1-z)-\ln(1-\lambda z)}{1-z}\right]_++ 2
                                                                                                    \left(-1+\ln\frac{K_F^2}{\lambda}\right)\left[\frac{1}{1-z}\right]_+\\
     \\ &-& \displaystyle 2\left[\frac{1}{1-\lambda z}\right]_+
            +\frac{1}{2}\left[\frac{1-z}{(1-\lambda z)^2}\right]_+\bigg\} + 2C_F
            \left[\lambda K_A \delta(1-z) + (1+\lambda)z\right]\,.
   \end{array}

\ In order to proceed, it is useful to work out the effect of the
:math:`+`-prescription in the presence of an incomplete integration:

.. math::
   \displaystyle \int_x^1 dz\left[f(z)\right]_+g(z) =
     \int_x^1
     dz\,f(z)\left[g(z)-g(1)\right]-g(1)\underbrace{\int_0^xdz\,f(z)}_{-R_f(x)}=\displaystyle \int_x^1
     dz\left\{\left[f(z)\right]_{x}+R_f(x)\delta(1-z)\right\}g(z)\,.

where the :math:`x`-prescription in the r.h.s. of the equation above
should be understood in as a usual :math:`+`-prescription regardless of
the integration bounds.

Despite most of the time the residual :math:`R_f(x)` function can be
evaluated analytically, sometimes it needs to be evaluated numerically.
The :math:`+`-prescribed functions that present in the expressions above
give raise to the following residual functions that can be computed
analytically:

.. math:: -\int_0^x\frac{dz}{1-z} = \ln(1-x)\,,

.. math::
   -\int_0^xdz\frac{\ln(1-z)}{1-z} = \frac12
     \ln^2(1-x)\,,

.. math::
   -\int_0^x\frac{dz}{1-\lambda z} =
     \frac{1}{\lambda}\ln(1-\lambda x)\,,

.. math::
   -\int_0^xdz\frac{1-z}{(1-\lambda z)^2} =
     \frac{1}{\lambda^2}\ln(1-\lambda
     x)+\frac{1-\lambda}{\lambda}\frac{x}{1-\lambda x}\,,

while the following integral: [3]_

.. math:: R(x)=-\int_0^xdz\frac{\ln(1-\lambda z)}{1-z}

is computed numerically. The advantage of the :math:`x`-prescription is
that, when convoluting the coefficient functions above with PDFs at the
point :math:`x`, one can treat the :math:`+`-prescribed functions using
the standard definition at the price of adding to the local terms the
following functions:

.. math::
   \begin{array}{rcl}
     C^{(1)}_{1,q}&\rightarrow& \displaystyle
                                C^{(1)}_{1,q} +
                                2C_F\bigg[2\ln^2(1-x)-2R(x)+2\left(-1+\ln\frac{K_F^2}{\lambda}\right)\ln(1-x)\\
     \\ &&\displaystyle +\frac{\lambda-1}{\lambda^3}\ln(1-\lambda
           x)+\frac{1}{2\lambda^2}\ln(1-\lambda
           x)+\frac{1-\lambda}{2\lambda}\frac{x}{1-\lambda x}\bigg]\delta(1-z) \,,
   \end{array}

.. math::
   \begin{array}{rcl}
     C^{(1)}_{2,q}&\rightarrow& \displaystyle
                                C^{(1)}_{2,q} +
                                2C_F\bigg[2\ln^2(1-x)-2R(x)+2\left(-1+\ln\frac{K_F^2}{\lambda}\right)\ln(1-x)\\
     \\ &&\displaystyle
           +\frac{2\lambda^2-\lambda-1}{\lambda^2}\ln(1-\lambda
           x)+\frac{1}{2\lambda^2}\ln(1-\lambda
           x)+\frac{1-\lambda}{2\lambda}\frac{x}{1-\lambda x}\bigg]\delta(1-z)
   \end{array}\,,

.. math::
   \begin{array}{rcl}
     C^{(1)}_{3,q}&\rightarrow& \displaystyle
                                C^{(1)}_{3,q} +
                                2C_F\bigg[2\ln^2(1-x)-2R(x)+2\left(-1+\ln\frac{K_F^2}{\lambda}\right)\ln(1-x)\\
     \\ &&\displaystyle +\frac{\lambda-1}{\lambda^2}\ln(1-\lambda
           x)+\frac{1}{2\lambda^2}\ln(1-\lambda
           x)+\frac{1-\lambda}{2\lambda}\frac{x}{1-\lambda x}\bigg]\delta(1-z) \,,
   \end{array}

.. math::
   \begin{array}{rcl}
     C^{(1)}_{L,q}&\rightarrow& \displaystyle
                                C^{(1)}_{L,q} +
                                2C_F(1-\lambda)\bigg[2\ln^2(1-x)-2R(x)+2\left(-1+\ln\frac{K_F^2}{\lambda}\right)\ln(1-x)\\
     \\ &&\displaystyle -\frac{2}{\lambda}\ln(1-\lambda
           x)+\frac{1}{2\lambda^2}\ln(1-\lambda
           x)+\frac{1-\lambda}{2\lambda}\frac{x}{1-\lambda x}\bigg]\delta(1-z) \,.
   \end{array}

Now let us consider the gluon coefficient functions. They read:

.. math::
   \begin{array}{rcl}
     C^{(1)}_{1,g}&=&\displaystyle
                      2T_R\bigg\{[z^2+(1-z)^2]\left[\ln\left(\frac{1-z}{z}\right)
                      -\frac12\ln(1-\lambda) +\frac12\ln\frac{K_F^2}{\lambda} \right]+\\ \\
                  &&\displaystyle 4z(1-z) - 1+\\ \\ &&\displaystyle
                                                       (1-\lambda)\left[-4z(1-z) + \frac{z}{1-\lambda z} +2z(1-2\lambda
                                                       z)\ln\frac{1-\lambda z}{(1-\lambda)z}\right]\bigg\}\,,
   \end{array}

\ 

.. math::
   \begin{array}{rcl}
     C^{(1)}_{2,g}&=&\displaystyle 2T_R
                      \bigg\{[z^2+(1-z)^2]\left[\ln\left(\frac{1-z}{z}\right)
                      -\frac12\ln(1-\lambda) +\frac12\ln\frac{K_F^2}{\lambda} \right]+\\ \\
                  &&\displaystyle 8z(1-z)- 1+\\ \\ &&\displaystyle
                                                      (1-\lambda)\left[-6(1+2\lambda) z(1-z)+\frac{1}{1-\lambda z} +
                                                      6\lambda z(1-2\lambda z)\ln\frac{1-\lambda
                                                      z}{(1-\lambda)z}\right]\bigg\}\,,
   \end{array}

.. math::
   \begin{array}{rcl}
     C^{(1)}_{3,g}&=&\displaystyle 2T_R
                      \bigg\{[z^2+(1-z)^2]\left[ 2\ln\left(\frac{1-z}{1-\lambda
                      z}\right)+\frac12\ln(1-\lambda)
                      +\frac12\ln\frac{K_F^2}{\lambda}\right]+\\ \\ &&\displaystyle
                                                                       (1-\lambda)\left[2 z(1-z) - 2z[1-(1+\lambda )z]\ln\frac{1-\lambda
                                                                       z}{(1-\lambda)z}\right]\bigg\}\,,
   \end{array}

.. math::
   \begin{array}{rcl}
     C^{(1)}_{L,g}&=&\displaystyle 2T_R
                      \bigg\{(1-\lambda)[z^2+(1-z)^2]\left[\ln\left(\frac{1-z}{z}\right)
                      -\frac12\ln(1-\lambda) +\frac12\ln\frac{K_F^2}{\lambda} \right]+\\ \\
                  &&\displaystyle 4(2-\lambda)z(1-z)+\\ \\ &&\displaystyle
                                                              (1-\lambda)\left[-2(3+4\lambda) z(1-z)+ 4\lambda z(1-2\lambda
                                                              z)\ln\frac{1-\lambda z}{(1-\lambda)z}\right]\bigg\}\,.
   \end{array}

Since these functions do not contain any :math:`+`-prescribed functions,
they can be implemented as they are.

We now consider the massless limit of the above massive coefficient
functions. In the limit :math:`m_H\rightarrow 0`, one has:

.. math::
   \begin{array}{l} 
     \lambda \rightarrow 1\quad\mbox{and}\quad K_A \rightarrow 0
   \end{array}\,.

\ For the quark coefficient functions this gives:

.. math::
   \begin{array}{rcl}
     C^{(1)}_{1,q}
     \displaystyle\mathop{\longrightarrow}_{m_H\rightarrow 0}
     C^{0,(1)}_{1,q} &=&\displaystyle 2C_F \bigg\{
                         -\left(\frac{9}{2}+2\zeta_2-\frac32 \ln K_F^2\right)\delta(1-z)\\ \\
                     &-&\displaystyle \frac{(1+z^2)\ln z}{1-z} -\left(\ln(1-z)+\ln
                         K_F^2\right)(1+z)+3\\ \\ &+&\displaystyle 2
                                                      \left[\frac{\ln(1-z)}{1-z}\right]_{x} -\left(\frac{3}{2}-2\ln
                                                      K_F^2\right)\left[\frac{1}{1-z}\right]_{x}\bigg\}\,,
   \end{array}

.. math::
   \begin{array}{rcl}
     C^{(1)}_{2,q}\displaystyle\mathop{\longrightarrow}_{m_H\rightarrow
     0}C^{0,(1)}_{2,q}&=&\displaystyle 2C_F \bigg\{
                          -\left(\frac{9}{2}+2\zeta_2-\frac32 \ln K_F^2\right)\delta(1-z)\\ \\
                      &-&\displaystyle \frac{(1+z^2)\ln z}{1-z} -\left(\ln(1-z)+\ln
                          K_F^2\right)(1+z)+2z+3\\ \\ &+&\displaystyle 2
                                                          \left[\frac{\ln(1-z)}{1-z}\right]_{x} -\left(\frac{3}{2}-2\ln
                                                          K_F^2\right)\left[\frac{1}{1-z}\right]_{x}\bigg\}\,,
   \end{array}

.. math::
   \begin{array}{rcl}
     C^{(1)}_{3,q}\displaystyle\mathop{\longrightarrow}_{m_H\rightarrow 0}
     C^{0,(1)}_{3,q}&=&\displaystyle 2C_F \bigg\{
                         -\left(\frac{9}{2}+2\zeta_2-\frac32 \ln K_F^2\right)\delta(1-z)\\ \\
                     &-&\displaystyle \frac{(1+z^2)\ln z}{1-z} -\left(\ln(1-z)+\ln
                         K_F^2\right)(1+z)+z+2\\ \\ &+&\displaystyle 2
                                                        \left[\frac{\ln(1-z)}{1-z}\right]_{x} -\left(\frac{3}{2}-2\ln
                                                        K_F^2\right)\left[\frac{1}{1-z}\right]_{x}\bigg\}\,,
   \end{array}

.. math::
   C^{(1)}_{L,q}\mathop{\longrightarrow}_{m_H\rightarrow
   0} C^{0,(1)}_{L,q} = 4C_F z \,.

Considering that:

.. math::
   R(x) \mathop{\longrightarrow}_{m_H\rightarrow 0}
     \frac12\ln(1-x)^2\,,

\ the local terms to be added to the quark coefficient functions are:

.. math::
   C^{0,(1)}_{1,q}\rightarrow C^{0, (1)}_{1,q} +
     2C_F\left[\ln^2(1-x)-\left(\frac32-2\ln
         K_F^2\right)\ln(1-x)\right]\delta(1-z)\,,

.. math::
   C^{0, (1)}_{2,q}\rightarrow C^{0, (1)}_{2,q} +
     2C_F\left[\ln^2(1-x)-\left(\frac32-2\ln
         K_F^2\right)\ln(1-x)\right]\delta(1-z)\,,

.. math::
   C^{0, (1)}_{3,q}\rightarrow C^{0, (1)}_{3,q} +
     2C_F\left[\ln^2(1-x)-\left(\frac32-2\ln
         K_F^2\right)\ln(1-x)\right]\delta(1-z)\,,

while no local term needs to be added to :math:`C^{0, (1)}_{L,q}`.

Now we turn to consider the zero-mass limit of the gluon coefficient
functions for which we need to know that:

.. math::
   \ln(1-\lambda)
     \mathop{\longrightarrow}_{m_H\rightarrow 0}
     -\ln\left(\frac{Q^2}{m_H^2}\right)\,.

\ This term is clearly divergent in the zero-mass limit and embeds a
collinear divergence typical of any massless calculations. However, we
retain all such terms so that:

.. math::
   C^{(1)}_{1,g}\mathop{\longrightarrow}_{m_H\rightarrow
       0} C^{0, (1)}_{1,g} =
     2T_R\left\{[z^2+(1-z)^2]\left[\ln\left(\frac{1-z}{z}\right) +\frac12
         \ln\left(\frac{Q^2}{m_H^2}\right) +\frac12\ln K_F^2\right]+ 4z(1-z) -
       1\right\}\,,

.. math::
   C^{(1)}_{2,g}\mathop{\longrightarrow}_{m_H\rightarrow
       0} C^{0, (1)}_{2,g} =
     2T_R\left\{[z^2+(1-z)^2]\left[\ln\left(\frac{1-z}{z}\right) +\frac12
         \ln\left(\frac{Q^2}{m_H^2}\right) +\frac12\ln K_F^2\right]+ 8z(1-z) -
       1\right\}\,,

.. math::
   C^{(1)}_{3,g}\mathop{\longrightarrow}_{m_H\rightarrow
       0} C^{0, (1)}_{3,g} = 2T_R[z^2+(1-z)^2]\left[-\frac12
       \ln\left(\frac{Q^2}{m_H^2}\right) +\frac12\ln K_F^2\right]\,,

.. math::
   C^{(1)}_{L,g}\mathop{\longrightarrow}_{m_H\rightarrow
       0} C^{0,(1)}_{L,g} = 2T_R\left[4z(1-z)\right]\,.

We also note that in the limit :math:`m_H\rightarrow 0`, the convolution
integrals in Eqs. :eq:`F1`, :eq:`F2`, :eq:`F3`
and :eq:`FL` will extend from :math:`x` to 1 rather than from
:math:`\chi` to 1.

The massive structure functions will eventually need to be combined to
the massless and the massive-zero ones to compute the GM structure
functions. Since the latter are computed through Mellin convolutions
whose integral extends from :math:`x` to 1, it is convenient to rewrite
Eqs. :eq:`F1`, :eq:`F2`, :eq:`F3`
and :eq:`FL` in such a way that the lower integration bound is
also :math:`x` rather than :math:`\chi`. To this end, let us consider
the integral:

.. math::
   I=\int_\chi^1\frac{dy}{y}
     C(y)f\left(\frac{\chi}{y}\right)\,,

\ where :math:`\chi=x/\lambda`. By changing of integration variable
:math:`z = \lambda y`, integral above becomes:

.. math::
   I=\int_x^\lambda\frac{dz}{z}
     C\left(\frac{z}{\lambda}\right)f\left(\frac{x}{y}\right) =
     \int_x^1\frac{dz}{z}
     \widetilde{C}(z,\lambda)f\left(\frac{x}{y}\right)\,,

where:

.. math:: \widetilde{C}(z,\lambda)=\theta(\lambda-z)C\left(\frac{z}{\lambda}\right)\,.

In this way we have achieved the goal of expressing the “reduced”
convolution in Eqs. :eq:`F1`, :eq:`F2`, :eq:`F3`
and :eq:`FL` as a “standard” convolution between :math:`x` and
1. As already mentioned above, this does not need to be done in the
massive-zero case as the convolution already extends between :math:`x`
and :math:`1`.

Target Mass Corrections
=======================

Kinematic corrections due to the finite mass of the hadron :math:`M_p`
which recoils against the vector boson might be relevant in the
small-:math:`Q` region. The leading contributions to these corrections
have been computed long time ago in Ref. (Georgi and Politzer 1976) and,
denoting the target-mass corrected structure functions with the symbol
:math:`\widetilde{\quad}`, they take the form:

.. math::
   \begin{array}{rcl}
     \displaystyle \widetilde{F}_2(x,Q) &=&
                                            \displaystyle \frac{x^2}{\xi^2 \tau^{3/2}} F_2(\xi,Q) + \frac{6\rho
                                            x^3}{\tau^2} I_2(\xi,Q)\,,\\ \\ \displaystyle \widetilde{F}_L(x,Q) &=&
                                                                                                                   \displaystyle F_L(\xi,Q)+\frac{x^2(1-\tau)}{\xi^2 \tau^{3/2}}
                                                                                                                   F_2(\xi,Q) + \frac{\rho x^3(6-2\tau)}{\tau^2} I_2(\xi,Q)\,,\\ \\
     \displaystyle x\widetilde{F}_3(x,Q) &=& \displaystyle \frac{x^2}{\xi^2
                                             \tau} \xi F_3(\xi,Q) + \frac{2\rho x^3}{\tau^{3/2}} I_3(\xi,Q)\,,
   \end{array}

\ where:

.. math::
   \rho = \frac{M_p^2}{Q^2}\,,\qquad \tau = 1 + 4\rho
     x^2\,,\qquad \xi = \frac{2x}{1+\sqrt{\tau}}\,,

and:

.. math::
   I_2(\xi,Q) = \int_\xi^1dy \frac{F_2(y,Q)}{y^2}\,,
     \qquad I_3(\xi,Q) = \int_\xi^1dy \frac{yF_3(y,Q)}{y^2}\,.

These integrals can be efficiently computed relying on the procedure
based on the interpolating functions discussed in the part of the
documentation devoted to the interpolation. From the equations above, it
is clear that in the limit :math:`M_p \rightarrow 0`, that implies
:math:`\rho \rightarrow 0`, :math:`\tau \rightarrow 1`, and
:math:`\xi \rightarrow x`, all structure functions reduce to the usual
formulas.

Renormalisation and factorisation scale variations
==================================================

In the previous sections, when discussing the implementation of the
structure functions in ``APFEL++``, we implicitly assumed that the
renormalisation scale :math:`\mu_R` and the factorisation scale
:math:`\mu_F` were identified with the scale :math:`Q`. The purpose of
this section is to relax this assumption. To do so, it is necessary to
consider the expansion of the DGLAP and of the renormalisation-group
(RG) equation for :math:`\alpha_s` up to :math:`\mathcal{O}(\alpha_s^2)`
that is the maximum order considered in the structure functions
discussed above. These equations read:

.. math::
   \frac{\partial f_{i}}{\partial\ln\mu_F^2} =
     \frac{\alpha_s(\mu_F)}{4\pi}\left[ P_{ij}^{(0)}(x) +
       \frac{\alpha_s(\mu_F)}{4\pi} P_{ij}^{(1)}(x) + \dots\right]\otimes
     f_j(x,\mu_F)\,,

\ where a summation over repeated indices is understood, and:

.. math::
   \frac{\partial
     }{\partial\ln\mu_R^2}\left(\frac{\alpha_s}{4\pi}\right) =
     -\left(\frac{\alpha_s(\mu_R)}{4\pi}\right)^2\left[ \beta_0 +
       \frac{\alpha_s(\mu_R)}{4\pi}\beta_1 + \dots\right]\,.

Defining:

.. math::
   \xi_R\equiv\frac{\mu_R}{Q}\,,\quad\xi_F\equiv\frac{\mu_F}{Q}\quad\mbox{and}\quad
     a_s=\frac{\alpha_s}{4\pi}

and:

.. math::
   t_R \equiv \ln\xi_R^2\quad\mbox{and}\quad t_F \equiv
     \ln\xi_F^2\,,

they can be written more compactly as:

.. math::
   \frac{\partial f_{i}}{\partial t_F}
     = a_s(t_F)\left[ P_{ij}^{(0)} + a_s(t_F) P_{ij}^{(1)} +
       \dots\right]\otimes f_j(t_F)\,,
   :label: DGLAPsimp

and:

.. math::
   \frac{\partial a_s}{\partial t_R} =
     -a_s^2(t_R)\left[ \beta_0 + a_s(t_R)\beta_1 + \dots\right]\,.
   :label: BETAsimp

Now, the Taylor expansion of :math:`f_i(t)` around :math:`t=t_F`, up to
second order is:

.. math::
   f_i(t) = f_i(t_F)+\frac{\partial
       f_{i}}{\partial t}\bigg|_{t=t_F} (t-t_F) + \frac12 \frac{\partial^2
       f_{i}}{\partial t^2}\bigg|_{t=t_F} (t-t_F)^2 + \dots
   :label: expDGLAP

Using Eqs. :eq:`DGLAPsimp`
and :eq:`BETAsimp`, we have that:

.. math::
   \begin{array}{l} 
     \displaystyle \frac{\partial f_{i}}{\partial
     t}\bigg|_{t=t_F} = \left[ a_s(t_F) P_{ij}^{(0)} + a_s^2(t_F)
     P_{ij}^{(1)}\right]\otimes f_j(t_F) + \mathcal{O}(a_s^3)\\ \\
     \displaystyle \frac{\partial^2 f_{i}}{\partial t^2}\bigg|_{t=t_F} =
     a_s^2(t_F)\left[ P_{il}^{(0)}\otimes P_{lj}^{(0)} - \beta_0
     P_{ij}^{(0)} \right]\otimes f_j(t_F) + \mathcal{O}(a_s^3)
   \end{array}

Choosing :math:`t=0` in Eq. :eq:`expDGLAP`, that is
equivalent to setting :math:`\mu_F=Q`, gives:

.. math::
   f_i(0) = \left\{1-a_s(t_F) t_F
       P_{ij}^{(0)} + a_s^2(t_F)\left[-t_F P_{ij}^{(1)}+ t_F^2 \frac12\left(
           P_{il}^{(0)}\otimes P_{lj}^{(0)} - \beta_0 P_{ij}^{(0)}
         \right)\right]\right\}\otimes f_j(t_F) + \mathcal{O}(a_s^3)\,.
   :label: expDGLAP1

In addition, using Eq. :eq:`BETAsimp`, one easily finds:

.. math::
   a_s(t_F)=
     a_s(t_R)\left[1+a_s(t_R)\beta_{0}(t_R-t_F)+\mathcal{O}(a_s^2)\right]\,,
   :label: BETAexp

which can be plugged into Eq. :eq:`expDGLAP1` to give:

.. math::
   f_i(0) = \left\{1-a_s(t_R) t_F
       P_{ij}^{(0)} + a_s^2(t_R)\left[-t_F P_{ij}^{(1)}+ t_F^2 \frac12\left(
           P_{il}^{(0)}\otimes P_{lj}^{(0)} + \beta_0 P_{ij}^{(0)} \right)-
         t_Ft_R\beta_{0}P_{ij}^{(0)}\right]\right\}\otimes f_j(t_F) +
     \mathcal{O}(a_s^3)\,.
   :label: expDGLAP2

It is also useful to consider :math:`t_F=0` in
Eq. :eq:`BETAexp`, that is equivalent to set
:math:`\mu_R=Q`, which gives:

.. math::
   a_s(0)=
     a_s(t_R)\left[1+a_s(t_R)\beta_{0}t_R+\mathcal{O}(a_s^2)\right]\,.
   :label: BETAexp1

We are now ready to use these equations to derive the scale variation
terms to be included in ZM structure functions. Truncating the
perturbative series to :math:`\mathcal{O}(\alpha_s^2)`, they are written
in terms of PDFs and coefficient functions as:

.. math::
   F(t_R,t_F) / x = \left[\sum_{k=0}^{2} a_s^k(t_R)
       \widetilde{\mathcal{C}}_i^{(k)}(t_R,t_F)\right]\otimes f_i(t_F)\,,
   :label: NonZeroScales

\ where the symbol :math:`\otimes` represents the convolution that is
not necessary to write explicitly here. Since structure functions are
physically observable, they must be renormalisation and factorisation
scale invariant, that is:

.. math::
   F(t_R,t_F) = F(0,0)\,,
   :label: invariance

up to subleading, *i.e.* :math:`\mathcal{O}(\alpha_s^3)`, corrections.
Since:

.. math::
   F(0,0) / x = \left[\sum_{k=0}^{2}
       a_s^k(0) \widetilde{C}_i^{(k)}\right]\otimes
     f_i(0)\,,
   :label: ZeroScales

where :math:`\widetilde{C}_i^{(k)}` are the usual perturbative
contributions to the ZM coefficient functions, one can plug
Eqs. :eq:`expDGLAP2` and :eq:`BETAexp1`
into Eq. :eq:`ZeroScales` and impose the identity in
Eq. :eq:`invariance`. By doing so, one finds the
explicit expression of the “generalised” coefficient functions
:math:`\widetilde{\mathcal{C}}_i^{(k)}(t_R,t_F)` that include also the
scale variation terms. The result is:

.. math::
   \begin{array}{rcl}
       F(0,0)/x &=&\bigg\{ \widetilde{C}_j^{(0)} \\ \\
                &+&\displaystyle a_s(t_R)\left[\widetilde{C}_j^{(1)}- t_F
                    \widetilde{C}_i^{(0)} \otimes P_{ij}^{(0)}\right]\\ \\
                &+&\displaystyle a_s^2(t_R)\bigg[\widetilde{C}_j^{(2)} + t_R\beta_0
                    \widetilde{C}_j^{(1)} -t_F \left(\widetilde{C}_i^{(0)} \otimes
                    P_{ij}^{(1)}+\widetilde{C}_i^{(1)} \otimes P_{ij}^{(0)}\right)\\ \\
                &+&\displaystyle \frac{t_F^2}2 \widetilde{C}_i^{(0)} \otimes \left(
                    P_{il}^{(0)}\otimes P_{lj}^{(0)} + \beta_0 P_{ij}^{(0)} \right)-
                    t_Ft_R\beta_{0}\widetilde{C}_i^{(0)} \otimes
                    P_{ij}^{(0)}\bigg]\bigg\}\otimes f_j(t_F)+\mathcal{O}(a_s^3)\,.
   \end{array}

Finally, using the identity in Eq. :eq:`invariance`, it
is easy to find that:

.. math::
   \begin{array}{rcl}
     \displaystyle
     \widetilde{\mathcal{C}}_j^{(0)}(t_R,t_F) &=& \displaystyle
                                                  \widetilde{C}_j^{(0)} \\ \\ \displaystyle
     \widetilde{\mathcal{C}}_j^{(1)}(t_R,t_F) &=& \displaystyle
                                                  \widetilde{C}_j^{(1)}-t_F \widetilde{C}_i^{(0)} \otimes P_{ij}^{(0)}
     \\ \\ \displaystyle \widetilde{\mathcal{C}}_j^{(2)}(t_R,t_F) &=&
                                                                      \displaystyle \widetilde{C}_j^{(2)} + t_R\beta_0 \widetilde{C}_j^{(1)}
                                                                      -t_F \left(\widetilde{C}_i^{(0)} \otimes
                                                                      P_{ij}^{(1)}+\widetilde{C}_i^{(1)} \otimes P_{ij}^{(0)}\right)\\ \\
                                              &+&\displaystyle\frac{t_F^2}2 \widetilde{C}_i^{(0)} \otimes \left(
                                                  P_{il}^{(0)}\otimes P_{lj}^{(0)} + \beta_0 P_{ij}^{(0)} \right)-
                                                  t_Ft_R\beta_{0}\widetilde{C}_i^{(0)} \otimes P_{ij}^{(0)}\,.
   \end{array}
   :label: generalizedCF

Unsurprisingly, setting :math:`\mu_F=\mu_R=Q`, that results in
:math:`t_F=t_R=0`, one finds
:math:`\widetilde{\mathcal{C}}_j^{(k)}(0,0) =\widetilde{C}_j^{(k)}` as
required by construction.

In order to provide an operative formulation of scale variations, it is
necessary to specify the basis in which PDFs are expressed. The
preferred choice in ``APFEL`` is the so-called QCD evolution basis:

.. math:: \{g,\Sigma,V,T_{3},V_{3},T_{8},V_{8},T_{15},V_{15},T_{24},V_{24},T_{35},V_{35}\}\,.

The distributions in the QCD evolution basis can be written in terms
distributions in the more familiar “physical” basis, *i.e.*
:math:`\{\overline{t},\overline{b},\overline{c},\overline{s},\overline{u},\overline{d},d,u,s,c,b,t\}`,
as follows:

.. math::
   \begin{array}{rcl}
   \Sigma&=& \sum_{q}q^+\,,\\
   V &=& \sum_{q}q^-\,,\\
   T_3&=& u^+-d^+\,,\\
   V_3&=& u^--d^-\,,\\
   T_8&=& u^++d^+-2s^+\,, \\
   V_8&=& u^-+d^- -2s^-\,, \\
   T_{15}&=& u^++d^++s^+-3c^{+}\,, \\
   V_{15}&=& u^-+d^- +s^--3c^{-}\,, \\
   T_{24}&=& u^++d^++s^++c^{+}-4b^+\,, \\
   V_{24}&=& u^-+d^- +s^-+c^{-}-4b^-\,, \\
   T_{35}&=& u^++d^++s^++c^{+}+b^+-5t^{+}\,, \\
   V_{35}&=& u^-+d^- +s^-+c^{-}+b^--5t^{-}\,.
   \end{array}

where the notation :math:`q^{\pm}\equiv q\pm \overline{q}` is used. As
the name suggests, the QCD evolution basis is particularly useful when
evolving PDFs because in this basis the DGLAP evolution equations take a
maximally diagonalised form. Adopting the QCD evolution basis implies
that the indices :math:`i`, :math:`j`, and :math:`l` in
Eq. :eq:`generalizedCF` run between 0 and 12 over
this basis with 0 corresponding to the gluon component.

In the following, we will consider a ZM neutral current structure
function for :math:`Q\ll M_Z` in such a way that only the photon
contributes. We will extend the treatment to the general neutral-current
case and to the charge-current one below. Omitting for simplicity the
convolution sign and an overall factor :math:`x`, the starting point is
the definition of the structure function in terms of the PDFs in the
physical basis that reads:

.. math::
   F=\langle e_q^2 \rangle \left\{C_gg +\sum_{i=u}^{t}\underbrace{\theta(Q^2-m_i^2)\left[C_{\rm PS}+\frac{e_i^2}{\langle e_q^2
         \rangle}C_+\right]}_{\hat{C}_i}q_i^+\right\}\,,
   :label: StructFuncDef

\ where :math:`C_+` and :math:`C_{\rm PS}` correspond to the non-singlet
and pure-singlet coefficient functions, that are usually the quantities
computed in perturbation theory, and where:

.. math:: \langle e_q^2 \rangle = \sum_{i=u}^t e_i^2\theta(Q^2-m_i^2)\,.

with :math:`m_i` is the mass of the :math:`i`-th quark flavour. Now, in
order to express the structure function in
Eq. :eq:`StructFuncDef` in the evolution basis, we
need to find the tranformation :math:`T` such that:

.. math::
   q_i^+ = \sum_{j=1}^6T_{ij}f_j\,,
   :label: Rotation

where :math:`f_j` belongs to the evolution basis, that is:
:math:`f_1=\Sigma`, :math:`f_2=T_3`, :math:`f_3=T_8` and so on. One can
show that the trasformation matrix :math:`T` takes the form:

.. math::
   \begin{array}{l}
   \displaystyle T_{ij}=\theta_{ji}\frac{1-\delta_{ij}j}{j(j-1)}\quad j\geq 2\,,\\
   \\
   \displaystyle T_{i1} = \frac{1}{6}\,,
   \end{array}
   :label: TransDef

with :math:`\theta_{ji}=1` for :math:`j\geq i` and zero otherwise. In
addition, one can show that:

.. math:: \sum_{j=1}^6T_{ij} = 0\,,\quad\mbox{and}\quad \sum_{i=1}^6T_{ij} = \delta_{1j}\,.

Now, we can plug Eq. :eq:`Rotation` into
Eq. :eq:`StructFuncDef` and, using
Eq. :eq:`TransDef`, we get:

.. math::
   F=\langle e_q^2 \rangle \left\{C_gg +\frac16\left(C_++n_f C_{\rm
         PS}\right)\Sigma+\sum_{j=2}^{6}\frac{1}{j(j-1)}\left[\sum_{i=1}^j\hat{C}_i-j\hat{C}_j\right]f_j\right\}\,,
   :label: StructFuncDefEvol

where we have transmuted the sum over :math:`u`, :math:`d` and so on
into a sum between 1 and 6 and we have defined the number of active
flavours :math:`n_f` as:

.. math:: n_f = \sum_{i=1}^6\theta(Q^2-m_i^2)\,.

We can now express the term in square brackets in terms :math:`C_+` and
:math:`C_{\rm PS}`. In particular:

.. math::
   \sum_{i=1}^j\hat{C}_i-j\hat{C}_j = \sum_{i=1}^j \theta(Q^2-m_i^2)\left(C_{\rm PS}+\frac{e_i^2}{\langle e_q^2
         \rangle}C_+\right) - j \theta(Q^2-m_j^2)\left(C_{\rm PS}+\frac{e_j^2}{\langle e_q^2
         \rangle}C_+\right)\,.

We can distinguish two cases. The first is :math:`Q^2<m_j^2` and under
this assumption we have:

.. math:: \sum_{i=1}^j\hat{C}_i-j\hat{C}_j = C_++n_fC_{\rm PS}\,.

If instead :math:`Q^2 \geq m_j^2`, then:

.. math:: \sum_{i=1}^j\hat{C}_i-j\hat{C}_j = K_j C_+\,,

with:

.. math:: K_j=\frac{1}{\langle e_q^2\rangle}\left(\sum_{i=1}^{j}e_i^2-je_j^2\right)=\frac{1}{\langle e_q^2\rangle}\left(\sum_{i=1}^{j-1}e_i^2-(j-1)e_j^2\right)\,.

Both cases can be gathered in one single formula as:

.. math::
   \sum_{i=1}^j\hat{C}_i-j\hat{C}_j = \theta(m_j^2-Q^2-\epsilon)\left[ C_++n_f C_{\rm PS}\right]
   +\theta(Q^2-m_j^2)\left[K_jC_+\right]\,.
   :label: SumCoef

where :math:`\epsilon` is an infinitesimal positive number that ensures
that the case :math:`Q^2=m_j^2` is included in the second term of the
r.h.s. of Eq. :eq:`SumCoef`. Eq. :eq:`SumCoef`
can be rewritten as:

.. math:: \sum_{i=1}^j\hat{C}_i-j\hat{C}_j = \theta_{n_fj}\left[K_jC_+\right]+\theta_{j,n_f-1}\left[ C_++n_f C_{\rm PS}\right]\,.

In addition, neglecting intrinsic heavy-quark contributions, one can
easily see that:

.. math::
   f_j = \theta_{n_fj}f_j+\theta_{j,n_f+1}\Sigma\,,
   :label: SingletReduction

and thus:

.. math::
   \sum_{j=2}^{6}\frac{1}{j(j-1)}\left[\sum_{i=1}^j\hat{C}_i-j\hat{C}_j\right]f_j
   = C_+\left[\sum_{j=2}^{n_f}\frac{K_j}{j(j-1)}f_j\right]+\left[\sum_{j=n_f+1}^{6}\frac{1}{j(j-1)}\right] \left[ C_++n_f C_{\rm PS}\right]\Sigma\,.

But:

.. math:: \sum_{j=n_f+1}^{6}\frac{1}{j(j-1)}=\frac{1}{n_f}-\frac{1}{6}\,,

and moreover:

.. math:: \frac{K_j}{j(j-1)} = \frac{1}{\langle e_q^2\rangle}\frac{1}{j(j-1)}\left(\sum_{i=1}^je_i^2-je_j^2\right)=\frac{1}{\langle e_q^2\rangle}\underbrace{\frac{1}{j(j-1)}\sum_{i=1}^6e_i^2\left[\theta_{ji}-j\delta_{ij}\right]}_{d_j}\,.

Finally, putting all pieces together, gives:

.. math::
   \begin{array}{rcl}
   F&=&\displaystyle \langle e_q^2 \rangle \left[C_g g +\left(C_{\rm PS} +\frac1{n_f}C_+
     \right) \Sigma\right]+C_+\sum_{j=2}^{n_f} d_j f_j\\
   \\
   &=&\displaystyle \langle e_q^2
   \rangle \left[C_g g +\frac{1}{6}\left(C_+ + 6C_{\rm PS}
     \right) \Sigma\right]+C_+\sum_{j=2}^{6} d_j f_j\,.
   \end{array}
   :label: StructureFunctionEvol

The second version of the equation above is particularly useful for the
implementation because the flavours structure does not have any explicit
dependence on :math:`n_f` and does not rely on the absence of intrinsic
heavy-quark contributions
(Eq. :eq:`SingletReduction`). It is useful to
separate the contributions deriving from the coupling of the vector
boson to the different quark flavours. To do so, one just needs to
select the contributions proportional to, say, the :math:`k`-th charge
:math:`e_k^2`. This is easily done with the following replacement:

.. math:: e_i^2\rightarrow \delta_{ik}e_i^2\,.

This gives:

.. math:: \langle e_q^2 \rangle \rightarrow \theta(Q^2-m_k^2) e_k^2\,,

and:

.. math:: d_j \rightarrow \frac{e_k^2\left[\theta_{jk}-j\delta_{kj}\right]}{j(j-1)}=\theta(Q^2-m_k^2) e_k^2\frac{\left[\theta_{jk}-j\delta_{kj}\right]}{j(j-1)}\,,

so that the component of the structure function :math:`F` associated to
the :math:`k`-th quark flavour is:

.. math::
   \begin{array}{rcl}
   F^{(k)} &=&\displaystyle  \theta(Q^2-m_k^2) e_k^2\left\{\left[C_g g +\frac{1}{6}\left(C_++6 C_{\rm PS} \right) \Sigma\right]+C_+\sum_{j=2}^{6}
               \frac{\left[\theta_{jk}-j\delta_{kj}\right]}{j(j-1)} f_j
               \right\}\\
   \\
   &=& \displaystyle \theta(Q^2-m_k^2) e_k^2\left\{\left[C_g g +\frac{1}{6}\left(C_++6 C_{\rm PS} \right) \Sigma\right]-\frac{1}{k}C_+f_k+C_+\sum_{j=k+1}^{6}
               \frac{1}{j(j-1)} f_j
               \right\}
   \end{array}
   :label: eq:StructureFunctionSingle

and it is such that:

.. math::
   F = \sum_{k=1}^6 F^{(k)}\,.
   :label: eq:StructureFunctionSingleSum

Phenomenologically relevant combinations are the light component, that
includes the contribution of down, up, and strange, and the three heavy
quark components (even though the top component is typically not
relevant). The light component is defined as:

.. math:: F^l = \sum_{k=1}^3 F^{(k)}= \langle e_l^2 \rangle \left[C_g g +\frac{1}{6}\left(C_++6 C_{\rm PS} \right) \Sigma\right]+C_+\sum_{j=2}^{6} d_j^{(l)} f_j\,.

where:

.. math:: \langle e_l^2\rangle = \sum_{i=1}^3e_i^2\,,

and:

.. math::
   d_j^{(l)}=
   \frac{1}{j(j-1)}\sum_{i=1}^3e_i^2\left[\theta_{ji}-j\delta_{ij}\right]=
   \left\{
   \begin{array}{ll}
   \frac{1}{2}(e_u^2-e_d^2)\,,\quad& j= 2 \\
   \\
   \frac{1}{6}(e_u^2+e_d^2-2e_s^2)\,,\quad& j= 3 \\
   \\
   \frac{\langle e_l^2\rangle}{j(j-1)}\,,\quad & j\geq 4
   \end{array}
   \right.\,,

no need of the :math:`\theta`-functions as the scale :math:`Q` will
always be above the strange threshold. Therefore, the explicit form of
:math:`F^l` is:

.. math::
   F^l = \langle e_l^2 \rangle \left[C_g g +\frac{1}{6}\left(C_++6 C_{\rm
         PS} \right) \Sigma\right]+\frac{1}{2}(e_u^2-e_d^2)C_+T_3+\frac{1}{6}(e_u^2+e_d^2-2e_s^2)C_+T_8+\langle
   e_l^2 \rangle C_+\sum_{j=4}^{6} \frac{1}{j(j-1)} f_j\,.
   :label: LightSF

The heavy-quark components are instead defined as:

.. math::
   \begin{array}{rcl}
   F^c &=& \displaystyle \theta(Q^2-m_c^2) e_c^2\left\{\left[C_g g +\frac{1}{6}\left(C_++6 C_{\rm PS} \right) \Sigma\right]-\frac{1}{4}C_+T_{15}+C_+\sum_{j=5}^{6}
               \frac{1}{j(j-1)} f_j
               \right\}\,,\\
   \\
   F^b &=& \displaystyle \theta(Q^2-m_b^2) e_b^2\left\{\left[C_g g +\frac{1}{6}\left(C_++6 C_{\rm PS} \right) \Sigma\right]-\frac{1}{5}C_+T_{24}+C_+\sum_{j=6}^{6}
               \frac{1}{j(j-1)} f_j
               \right\}\,,\\
   \\
   F^t &=& \displaystyle \theta(Q^2-m_t^2) e_t^2\left\{\left[C_g g +\frac{1}{6}\left(C_++6 C_{\rm PS} \right) \Sigma\right]-\frac{1}{6}C_+T_{35}\right\}\,.
   \end{array}
   :label: HeavySF

To conclude the treatment of all the structure functions, it should be
mentioned that
Eq. :eq:`StructureFunctionEvol` is valid only
for :math:`F_2` and :math:`F_L`. However, :math:`F_3`, that appears when
weak contributions are included, can be easily derived following the
same steps with the only differences being that:

-  the distributions :math:`\{\Sigma,T_3,T_8,T_{15},T_{24},T_{35}\}`
   must be replaced with :math:`\{V,V_3,V_8,V_{15},V_{24},V_{35}\}`,

-  :math:`C_+` must be replaced with :math:`C_-`,

-  the gluon and the pure-singlet coefficient functions are identically
   zero,

-  the squared electric charges must be replaced with the appropriate
   electroweak charges :math:`c_i`.

Following this recipe, one finds:

.. math::
   F_3=\langle c_q^2 \rangle \frac1{6}C_- V+C_-\sum_{j=2}^{6} d_j
   g_j\,,
   :label: StructureFunctionEvolF3

\ where :math:`g_j` runs over
:math:`\{V,V_3,V_8,V_{15},V_{24},V_{35}\}`.

Despite
Eq. :eq:`eq:StructureFunctionSingle` has
been derived in the ZM scheme, it can be generalised to the massive
scheme. The complication is that, due to the explicit dependence of the
coefficient functions on the :math:`k`-th quark mass, [4]_ the
:math:`k`-th structure function will read:

.. math::
   F^{{\rm M},(k)} =  e_k^2\left\{\left[C_g^{(k)} g
       +\frac{1}{6}\left(C_+^{(k)} + 6C_{\rm PS}^{(k)}
     \right) \Sigma\right]+C_+^{(k)}\sum_{j=2}^{6}
               \frac{\left[\theta_{jk}-j\delta_{kj}\right]}{j(j-1)} f_j
               \right\}\,,
   :label: eq:StructureFunctionSingleHeavy

\ where an explicit dependence of the index :math:`k` has been
introduced in the coefficient functions. This dependence invalidates the
equality in
Eq. :eq:`eq:StructureFunctionSingleSum`
where :math:`F` is given in
Eq. :eq:`StructureFunctionEvol`. Therefore,
in order to compute inclusive structure functions in the massive case a
different combination is to be taken. Specifically one defines:

.. math::
   F^{{\rm M}} = \sum_{k=1}^6F^{{\rm M},(k)} =  \left[\langle C_g\rangle g +\frac{1}{6}\langle C_{q}\rangle \Sigma\right]+\sum_{j=2}^{6}
               \frac{\langle C_{\rm +}\rangle_j}{j(j-1)} f_j\,,
   :label: eq:StructureFunctionHeavy

with the definitions:

.. math::
   \begin{array}{rcl}
     \langle C_g\rangle &=&\displaystyle \sum_{k=1}^6e_k^2 C_g^{(k)}\,,\\
     \\
     \langle C_{q}\rangle &=&\displaystyle \sum_{k=1}^6e_k^2
                                  \left(C_+^{(k)} +6C_{\rm PS}^{(k)} \right)\,,\\
     \\
     \langle C_{\rm +}\rangle_j &=&\displaystyle \sum_{k=1}^6e_k^2\left[\theta_{jk}-j\delta_{kj}\right]C_+^{(k)}=\sum_{k=1}^{j-1}e_k^2C_+^{(k)}+e_j^2(1-j)C_+^{(j)}\,.
   \end{array}

Eq. :eq:`StructureFunctionEvol` is the result
that allows us to implement the scale variation formulae given in
Eq. :eq:`generalizedCF` in ``APFEL++``. An important
aspect of Eq. :eq:`StructureFunctionEvol` is
that it is written in terms of the fundamental coefficient functions
:math:`C_g`, :math:`C_+` and :math:`C_{\rm PS}` and PDFs appear in the
evolution basis in which the splitting-function matrix diagonalizes. In
particular, up to :math:`\mathcal{O}(\alpha_s^2)`, one has that:

.. math::
   \begin{array}{ll}
     P_{ij}^{(k)} \rightarrow P_{ij}^{(k)} &\quad i,j=g,q(\Sigma)\\
     P_{ij}^{(k)} \rightarrow \delta_{ij}P_+^{(k)} &\quad i,j=T_{3},T_{8},V_{15},T_{24},T_{35}\\
     P_{ij}^{(k)} \rightarrow \delta_{ij}P_-^{(k)} &\quad i,j=V,V_{3},V_{8},V_{15},V_{24},V_{35}
   \end{array}

\ Also, defining:

.. math:: C_q=C_{\rm PS} +\frac{1}{n_f}C_+\,,

one can connect Eq. :eq:`ZeroScales` and
Eq. :eq:`StructureFunctionEvol` by observing
that:

.. math::
   \begin{array}{ll} 
     \widetilde{C}_j^{(k)} \rightarrow \langle e_q^2\rangle C_j^{(k)}
     &\quad j=g,q(\Sigma)\\
     \widetilde{C}_j^{(k)}
     \rightarrow d_jC_+^{(k)} &\quad j=T_{3},T_{8},T_{15},T_{24},T_{35}\\
     \widetilde{C}_j^{(k)}\rightarrow d_jC_-^{(k)} &\quad j=V_{3},V_{8},V_{15},V_{24},V_{35}
   \end{array}

where we have also considered the “minus” distributions that appear in
the :math:`F_3` structure function. Of course, the same relations must
hold also for Eq. :eq:`NonZeroScales`:

.. math::
   \begin{array}{ll}
     \widetilde{\mathcal{C}}_j^{(k)} \rightarrow
     \langle e_q^2\rangle \mathcal{C}_j^{(k)} &\quad j=g,q(\Sigma)\\
     \widetilde{\mathcal{C}}_j^{(k)} \rightarrow d_j\mathcal{C}_+^{(k)}
                            &\quad j=T_{3},T_{8},T_{15},T_{24},T_{35}\\
     \widetilde{\mathcal{C}}_j^{(k)} \rightarrow d_j\mathcal{C}_-^{(k)}
                            &\quad j=V_{3},V_{8},V_{15},V_{24},V_{35}
   \end{array}

with

.. math:: \mathcal{C}_q=\mathcal{C}_{\rm PS} +\frac{1}{n_f}\mathcal{C}_+\,.

In addition, in the following we will make use of the following
identity:

.. math:: P_-^{(0)}=P_+^{(0)}=P_{qq}^{(0)}\,.

Now, considering that :math:`C_{\rm PS}` starts at
:math:`\mathcal{O}(\alpha_s^2)`, we can write:

.. math::
   \begin{array}{l}
     \displaystyle C_-^{(0)}(x) = C_+^{(0)}(x) =
     \Delta_{\rm SF}\delta(1-x)\\
     \displaystyle C_j^{(0)}(x) = \left(\Delta_{\rm
     SF}/n_f\right)\delta_{qj}\delta(1-x) \quad \mbox{for } j=q,g
   \end{array}

where :math:`\Delta_{\rm SF}=1` for :math:`F_2` and :math:`F_3` and
:math:`\Delta_{\rm SF} = 0` for :math:`F_L`. From
Eq. :eq:`generalizedCF` it follows that:

.. math::
   \begin{array}{rcl}
     \displaystyle \mathcal{C}_\pm^{(0)}(t_R,t_F) &=&
                                                      \displaystyle \Delta_{\rm SF}\delta(1-x) \\ \\ \displaystyle
     \mathcal{C}_\pm^{(1)}(t_R,t_F) &=& \displaystyle
                                        C_\pm^{(1)}-\Delta_{\rm SF} t_F P_{qq}^{(0)} \\ \\ \displaystyle
     \mathcal{C}_\pm^{(2)}(t_R,t_F) &=& \displaystyle C_\pm^{(2)} +
                                        t_R\beta_0 C_\pm^{(1)} -t_F \left(\Delta_{\rm SF}
                                        P_\pm^{(1)}+C_\pm^{(1)} \otimes P_{qq}^{(0)}\right)\\ \\
                                                  &+&\displaystyle\Delta_{\rm SF} \frac{t_F^2}2 \left(
                                                      P_{qq}^{(0)}\otimes P_{qq}^{(0)} + \beta_0 P_{qq}^{(0)} \right)-
                                                      \Delta_{\rm SF} t_Ft_R\beta_{0}P_{qq}^{(0)}\,,
   \end{array}
   :label: NonSingletCF

that can be rearranged as:

.. math::
   \begin{array}{rcl}
     \displaystyle \mathcal{C}_\pm^{(0)}(t_R,t_F) &=&
                                                      \displaystyle \Delta_{\rm SF}\delta(1-x) \\ \\ \displaystyle
     \mathcal{C}_\pm^{(1)}(t_R,t_F) &=& \displaystyle
                                        C_\pm^{(1)}-\Delta_{\rm SF} t_F P_{qq}^{(0)} \\ \\ \displaystyle
     \mathcal{C}_\pm^{(2)}(t_R,t_F) &=& \displaystyle C_\pm^{(2)} +
                                        t_R\beta_0 C_\pm^{(1)} -t_F C_\pm^{(1)} \otimes P_{qq}^{(0)}\\ \\
                                                  &+&\displaystyle\Delta_{\rm SF}\frac{t_F^2}2 \left(
                                                      P_{qq}^{(0)}\otimes P_{qq}^{(0)} - \beta_0 P_{qq}^{(0)} \right)-
                                                      \Delta_{\rm SF} t_F\left[P_\pm^{(1)} - (t_F-t_R)\beta_{0}
                                                      P_{qq}^{(0)}\right]\,.
   \end{array}
   :label: NonSingletCF1

The term in square brackets in the r.h.s. of the third line corresponds
to what we define :math:`\widetilde{\mathcal{P}}_\pm^{(1)}(t_R,t_F)`,
that is the NLO contribution to the non-singlet splitting functions in
the presence of scale variations (:math:`\mu_R\neq\mu_F`).

Let us now consider the singlet sector that becomes:

.. math::
   \begin{array}{rcl}
     \displaystyle {\mathcal{C}}_j^{(0)}(t_R,t_F) &=&
                                                      \displaystyle \frac{\Delta_{\rm SF}}{n_f}\delta_{qj}\delta(1-x) \\ \\
     \displaystyle {\mathcal{C}}_j^{(1)}(t_R,t_F) &=& \displaystyle
                                                      {C}_j^{(1)}- \frac{\Delta_{\rm SF}}{n_f} t_F P_{qj}^{(0)} \\ \\ \displaystyle
     {\mathcal{C}}_j^{(2)}(t_R,t_F) &=& \displaystyle {C}_j^{(2)} +
                                        t_R\beta_0 {C}_j^{(1)} -t_F {C}_i^{(1)} \otimes P_{ij}^{(0)}\\ \\
                                                  &+&\displaystyle\frac{\Delta_{\rm SF}}{n_f} \frac{t_F^2}2 \left(
                                                      P_{qi}^{(0)}\otimes P_{ij}^{(0)} - \beta_0 P_{qj}^{(0)} \right)-
                                                      \frac{\Delta_{\rm SF}}{n_f} t_F\widetilde{P}_{qj}^{(1)}\,,
   \end{array}
   :label: SingletCF

\ for :math:`j=g,q`. Taking into account
Eq. :eq:`NonSingletCF1` and considering also that
:math:`C_{\rm PS}^{(0)}=C_{\rm PS}^{(1)}=0`, it is easy to see that:

.. math::
   \begin{array}{rcl}
     \displaystyle {\mathcal{C}}_{\rm PS}^{(0)}(t_R,t_F)
     &=& \displaystyle 0 \\ \\ \displaystyle {\mathcal{C}}_{\rm
     PS}^{(1)}(t_R,t_F) &=& \displaystyle 0 \\ \\ \displaystyle
     {\mathcal{C}}_{\rm PS}^{(2)}(t_R,t_F) &=& \displaystyle {C}_{\rm
                                               PS}^{(2)}-t_FC_g^{(1)}\otimes
     P_{gq}^{(0)}+\frac{\Delta_{\rm SF}}{n_f}\frac{t_F^2}{2}P_{qg}^{(0)}\otimes P_{gq}^{(0)}-\frac{\Delta_{\rm SF}}{n_f}t_F\left[\widetilde{P}_{qq}^{(1)}-\widetilde{P}_{+}^{(1)}\right]\,.
   \end{array}
   :label: PureSingletCF

Now we consider the massive case. In the neutral-current sector the
leading-order coefficient functions :math:`C_i^{(0)}` are identically
zero which substantially simplifies the structure of the coefficient
functions in the presence of scale variations:

.. math::
   \begin{array}{rcl}
     \displaystyle \mathcal{C}_j^{(0)}(t_R,t_F) &=&
                                                    \displaystyle 0 \\ \\ \displaystyle \mathcal{C}_j^{(1)}(t_R,t_F) &=&
                                                                                                                         \displaystyle C_j^{(1)} \\ \\ \displaystyle
     \mathcal{C}_j^{(2)}(t_R,t_F) &=& \displaystyle C_j^{(2)} + t_R\beta_0
                                      C_j^{(1)} -t_F C_i^{(1)} \otimes P_{ij}^{(0)}\,.
   \end{array}

\ In addition, the factorisation scale variation terms are already
present in the implementation of the massive coefficient functions in
``APFEL++``. As a consequence, only the renormalisation variation terms
need to be implemented.

As far as the massive charged-current sector is concerned, no
:math:`\mathcal{O}(\alpha_s^2)` are presently available and thus only
the first two lines of Eq. :eq:`generalizedCF` are
actually required. Also in this case the factorisation scale variation
terms are already present in the implementation of the coefficient
functions and again this avoids the pre-computation of additional terms.

Single-inclusive :math:`e^+e^-` annihilation structure functions
================================================================

The implementation of the Single-Inclusive :math:`e^+e^-` Annihilation
(SIA) structure functions in ``APFEL++`` is not very complicated. The
reason is that SIA is structurally identical to the case of
deep-inelastic scattering (DIS) that was (implicitly) discussed above.
As a matter fact, one can regard SIA as the time-like counterpart of DIS
and the differences are only at the level of coefficient functions and
splitting functions. Presently, the coefficient functions for SIA are
known up to :math:`\mathcal{O}(\alpha_s^2)` (NNLO) in the zero-mass
scheme and they have been computed in Ref. (Mitov and Moch 2006) and the
:math:`x`-space expressions reported in Appendix C of that paper.

The way in which the SIA expressions are reported is slightly different
w.r.t. DIS.It is then useful to reduce the SIA expressions to the same
form of DIS in such a way to use the DIS-based structure of ``APFEL++``
also for SIA. In particular, the SIA cross section in Ref. (Mitov and
Moch 2006) is expressed in terms of the three structure functions:
:math:`F_T`, :math:`F_L` and :math:`F_A`. However, comparing the SIA
cross section with the DIS one it is easy to realize that defining:

.. math::
   \begin{array}{l}
   F_2(x,Q) = F_T(x,Q) + F_L(x,Q)\,,\\
    F_L(x,Q) =
   F_L(x,Q)\,,\\
   F_A(x,Q) = xF_3(x,Q)\,,
   \end{array}
   :label: SIAtoDIS

\ the SIA cross section reduces to the same structure of DIS. Upon this
identification, the usual factorised form for the structure functions
applies: [5]_

.. math::
   F_k(x,Q) = \sum_{j=q,g} x\int_x^1\frac{dy}{y}
     c_{k,j}(\alpha_s(Q),x)\mathcal{D}_j\left(\frac{x}{y},Q\right)\,,\quad\mbox{with}\quad
     k = 2,L,3\,,

where :math:`\mathcal{D}_j` is the fragmentation function of the flavour
:math:`j` and the coefficient functions :math:`c_{k,j}` admit the
perturbative expansion:

.. math::
   c_{k,j}(\alpha_s(Q),x) = \sum_{n=0}^N \left(\frac{\alpha_s(Q)}{4\pi}\right)^n
     c_{k,j}^{(n)}(x)\,.

The leading-order coefficient functions are trivially:

.. math::
   \begin{array}{l}
     c_{k,g}^{(0)}(x) = 0 \,,\quad k = 2,L,3\,,\\ \\
     c_{L,q}^{(0)}(x) = 0\,, \\ \\ c_{2,q}^{(0)}(x) = c_{3,q}^{(0)}(x) =
     \delta(1-x)\,.
   \end{array}

Now we consider the NLO coefficient functions. Their explicit
expressions are give in Eqs. (C.13)-(C.17) of Ref. (Mitov and Moch 2006)
but, in order to write them in a form suitable for the implementation in
``APFEL++``, we need to isolate regular, singular, and local terms and
finally combine them according to Eq. :eq:`SIAtoDIS`.

.. math::
   \begin{array}{lcl}
     \displaystyle c_{L,q}^{(1)}(x) = 2C_F\,, & &\\ \\
     \displaystyle c_{L,g}^{(1)}(x) = 2C_F\frac{4(1-x)}{x}\,, & &\\ \\
     \displaystyle c_{2,q}^{(1)}(x) = c_{T,q}^{(1)}(x) + c_{L,q}^{(1)}(x)
                                           &=& \displaystyle
                                               2C_F\bigg[2\left(\frac{\ln(1-x)}{1-x}\right)_+-\frac{3}{2}\left(\frac{1}{1-x}\right)_+
                                               - (1+x)\ln(1-x)\\ \\ & &\displaystyle +2\frac{1+x^2}{1-x}\ln x
                                                                        +\frac{5}{2} -\frac{3}{2} x
                                                                        +\left(4\zeta_2-\frac{9}{2}\right)\delta(1-x)\bigg]\,,\\ \\ \displaystyle
     c_{2,g}^{(1)}(x) = c_{T,g}^{(1)}(x) + c_{L,g}^{(1)}(x) &=&
                                                                \displaystyle 4C_F\frac{1+(1-x)^2}{x} \ln[x^2(1-x)]\,,\\ \\ \displaystyle
     c_{3,q}^{(1)}(x) &=& \displaystyle
                          2C_F\bigg[2\left(\frac{\ln(1-x)}{1-x}\right)_+-\frac{3}{2}\left(\frac{1}{1-x}\right)_+
                          - (1+x)\ln(1-x)\\ \\ & &\displaystyle +2\frac{1+x^2}{1-x}\ln x
                                                   +\frac{1}{2} -\frac{1}{2} x
                                                   +\left(4\zeta_2-\frac{9}{2}\right)\delta(1-x)\bigg]\,,\\ \\ \displaystyle
     c_{3,g}^{(1)}(x) = 0\,.
   \end{array}

\ The NNLO coefficient functions, despite implemented in ``APFEL++``,
are not reported here because too lengthy.

Longitudinally polarised structure functions
============================================

Let us now consider the differential cross sections for unpolarised and
polarised DIS (see *e.g.* Eq. (19.16) of Sec. 19 in Ref. (Olive and
others 2014)):

.. math::
   \begin{array}{lcl}
   \displaystyle \frac{d^2\sigma^i}{dxdy}
   & = &
   \displaystyle \frac{2\pi\alpha^2}{xyQ^2}\eta^i
   \left[
   +Y_+ F_2^i \mp Y_- x F_3^i - y^2 F_L^i
   \right]
   \\ \\
   \displaystyle \frac{d^2\Delta\sigma^i}{dxdy}
   & = &
   \displaystyle \frac{2\pi\alpha^2}{xyQ^2}\eta^i
   \left[
   -Y_+ g_4^i \mp Y_- 2x g_1^i + y^2 g_L^i
   \right]
   \mbox{\,,}
   \end{array}

\ where :math:`i={\rm NC, CC}`, :math:`Y_{\pm}=1 \pm (1-y)^2`,
:math:`\eta^{\rm NC}=1`, :math:`\eta^{\rm CC}=(1\pm \lambda)^2\eta_W`
(with :math:`\lambda=\pm 1` is the helicity of the incoming lepton and
:math:`\eta_W=\frac{1}{2}
\left(\frac{G_FM_W}{4\pi\alpha}\frac{Q^2}{Q^2+M_W^2} \right)^2`), and

.. math::
   \begin{array}{lcl}
   \displaystyle F_L^i & = & \displaystyle F_2^i - 2xF_1^i\\ \\
   \displaystyle F_L^i & = & \displaystyle g_4^i - 2xg_5^i
   \mbox{\,.}
   \end{array}

Because the same tensor structure occurs in the spin-dependent and
spin-independent parts of the DIS hadronic tensor (in the limit
:math:`M^2/Q^2\to 0`), the polarised cross section can be obtained from
the unpolarised cross section with the following replacement

.. math::
   \displaystyle F_2^i \rightarrow -2g_4^i
   \ \ \ \ \ \ \ \ \ 
   \displaystyle F_3^i \rightarrow +4g_1^i
   \ \ \ \ \ \ \ \ \
   \displaystyle F_L^i \rightarrow -2g_L^i 
   \mbox{\,.}

Note that the extra factor two is due to the fact that the total cross
section is an average over initial-state polarisations.

The *polarised* structure functions :math:`g_4`, :math:`g_1` and
:math:`g_L` are expressed as a convolution of coefficient functions,
:math:`\Delta c_{k,j}`, and polarised PDFs, :math:`\Delta f_j`, (summed
over all flavors :math:`j`)

.. math::
   g_k(x,Q) 
   =
   \sum_{j=q,g}x \int_x^1 \frac{dy}{y} \Delta c_{k,j}(\alpha_s(Q),x)
   \Delta f_j\left(\frac{x}{y},Q\right)
   \mbox{\,,}
   \ \ \ \ \ \ \ \ \ \ \
   {\rm with} 
   \ \ k=4,1,L
   \mbox{\,.}

\ The coefficient functions :math:`\Delta c_{k,j}` allow for the usual
perturbative expansion

.. math::
   \Delta c_{k,j}(\alpha_s(Q),x) 
   =
   \sum_{n=0}^N\left(\frac{\alpha_s(Q)}{4\pi}\right)^n\Delta c_{k,j}^{(n)}(x)
   \,\mbox{,}

where the coefficients :math:`\Delta c_{k,j}^{(n)}(x)` are known up to
NLO, *i.e.* :math:`n=1` (see *e.g.* (Florian and Rotstein Habarnau 2013)
and references therein). At LO they are

.. math::
   \begin{array}{lcl}
   \Delta c_{4,q}^{(0)}(x) = \Delta c_{1,q}^{(0)}(x) &=& \delta(1-x)
   \\ \\
   \Delta c_{L,q}^{(0)}(x) &=& 0
   \,\mbox{,}
   \\ \\
   \Delta c_{k,g}^{(0)}(x) &=& 0
   \ \ \ \ \ \ \ \ \ \
   {\rm with} \ k=4,1,L
   \,\mbox{.}
   \end{array}

At NLO they read

.. math::
   \begin{array}{lcl}
     %C4q
     \displaystyle \Delta c_{4,q}^{(1)}(x) 
     &=& 
     \displaystyle 2C_F\,
     \bigg\{
     2\left[\frac{\ln(1-x)}{1-x}\right]_+ 
     - \frac{3}{2}\left[\frac{1}{1-x}\right]_+
     - (1+x)\ln(1-x)
     \\ \\ 
     & &\displaystyle - \frac{1+x^2}{1-x}\ln x
     + 3 + 2x
     -\left(\frac{9}{2} + 2\zeta_2\right)\delta(1-x)\bigg\}
     \,\mbox{,}
     \\ \\ 
     %C4g
     \displaystyle \Delta c_{4,g}^{(1)}(x) 
     &=& 
     0
     \,\mbox{,} 
     \\ \\
     %C1q
     \displaystyle \Delta c_{1,q}^{(1)}(x) 
     &=& 
     \displaystyle 2C_F\,
     \bigg\{
     2\left[\frac{\ln(1-x)}{1-x}\right]_+ 
     - \frac{3}{2}\left[\frac{1}{1-x}\right]_+
     - (1+x)\ln(1-x)
     \\ \\ 
     & &\displaystyle - \frac{1+x^2}{1-x}\ln x
     + 2 + x
     -\left(\frac{9}{2} + 2\zeta_2\right)\delta(1-x)\bigg\}
     \,\mbox{,}
     \\ \\ 
     %C1g
     \displaystyle \Delta c_{1,g}^{(1)}(x) 
     &=&
     \displaystyle 4T_R
     \bigg\{(2x - 1)\ln \frac{1-x}{x} - 4x +3 \bigg\} 
     \,\mbox{,}
     \\ \\ 
     %CLq
     \displaystyle \Delta c_{L,q}^{(1)}(x) 
     &=& 
     2C_F\, 2x 
     \,\mbox{,}
     \\ \\
     %CLg
     \displaystyle \Delta c_{L,g}^{(1)}(x) 
     &=& 
     0 
     \,\mbox{.}
   \end{array}

In the NC case the couplings can be written as:

.. math::
   \begin{array}{l}
   \displaystyle B_q(Q^2) = -e_qA_q(V_e\pm \lambda A_e)P_Z+V_qA_q(V_e^2+A_e^2\pm2\lambda V_eA_e)P_Z^2 \,\mbox{,}\\
   \\
   \displaystyle D_q(Q^2) = \pm\frac12 \lambda e_q^2 - e_qV_q(A_e\pm\lambda V_e)P_Z +\frac12(V_q^2+A_q^2)\left[2V_eA_e\pm\lambda (V_e^2+A_e^2)\right]P_Z^2\,\mbox{.}
   \end{array}

\ where :math:`\lambda` corresponds to the polarisation of the incoming
lepton. It should be stressed that :math:`B_q` multiplies :math:`g_4`
and :math:`g_L` while :math:`D_q` multiplies :math:`g_1`.

**References**

**References**

.. container:: references csl-bib-body hanging-indent
   :name: refs

   .. container:: csl-entry
      :name: ref-Buza:1995ie

      Buza, M., Y. Matiounine, J. Smith, R. Migneron, and W. L. van
      Neerven. 1996. “Heavy quark coefficient functions at asymptotic
      values Q**2 :math:`>`\ :math:`>` m**2.” *Nucl. Phys. B* 472:
      611–58. https://doi.org/10.1016/0550-3213(96)00228-3.

   .. container:: csl-entry
      :name: ref-deFlorian:2012wk

      Florian, Daniel de, and Yamila Rotstein Habarnau. 2013. “Polarized
      semi-inclusive electroweak structure functions at
      next-to-leading-order.” *Eur. Phys. J. C* 73 (3): 2356.
      https://doi.org/10.1140/epjc/s10052-013-2356-3.

   .. container:: csl-entry
      :name: ref-Forte:2010ta

      Forte, Stefano, Eric Laenen, Paolo Nason, and Juan Rojo. 2010.
      “Heavy quarks in deep-inelastic scattering.” *Nucl. Phys. B* 834:
      116–62. https://doi.org/10.1016/j.nuclphysb.2010.03.014.

   .. container:: csl-entry
      :name: ref-Gao:2017kkx

      Gao, Jun. 2018. “Massive charged-current coefficient functions in
      deep-inelastic scattering at NNLO and impact on strange-quark
      distributions.” *JHEP* 02: 026.
      https://doi.org/10.1007/JHEP02(2018)026.

   .. container:: csl-entry
      :name: ref-Georgi:1976ve

      Georgi, Howard, and H.David Politzer. 1976. “Freedom at Moderate
      Energies: Masses in Color Dynamics.” *Phys. Rev. D* 14: 1829.
      https://doi.org/10.1103/PhysRevD.14.1829.

   .. container:: csl-entry
      :name: ref-Gluck:1996ve

      Gluck, M., S. Kretzer, and E. Reya. 1996. “The Strange sea density
      and charm production in deep inelastic charged current processes.”
      *Phys. Lett. B* 380: 171–76.
      https://doi.org/10.1016/0370-2693(96)00456-X.

   .. container:: csl-entry
      :name: ref-Mitov:2006wy

      Mitov, Alexander, and Sven-Olaf Moch. 2006. “QCD Corrections to
      Semi-Inclusive Hadron Production in Electron-Positron Annihilation
      at Two Loops.” *Nucl. Phys. B* 751: 18–52.
      https://doi.org/10.1016/j.nuclphysb.2006.05.018.

   .. container:: csl-entry
      :name: ref-Agashe:2014kda

      Olive, K. A., and others. 2014. “Review of Particle Physics.”
      *Chin. Phys. C* 38: 090001.
      https://doi.org/10.1088/1674-1137/38/9/090001.

   .. container:: csl-entry
      :name: ref-Riemersma:1994hv

      Riemersma, S., J. Smith, and W. L. van Neerven. 1995. “Rates for
      inclusive deep inelastic electroproduction of charm quarks at
      HERA.” *Phys. Lett. B* 347: 143–51.
      https://doi.org/10.1016/0370-2693(95)00036-K.

.. [1]
   Actually, if neglecting intrinsic heavy-quark components, at
   :math:`\mathcal{O}(\alpha_s)` there is one single coefficient
   function for each structure functions, available for example in
   available in Ref. (Forte et al. 2010), that is convoluted with the
   gluon distribution.

.. [2]
   This approximation will be released later.

.. [3]
   Actually, also this integral can be computed analytically but it
   involves the di-log and I will leave it like this for now.

.. [4]
   In general, the :math:`k`-th structure function will not only depend
   on the :math:`k`-th mass but also on the other masses. However, in
   the neutral-current case and neglecting intrinsic heavy-quark
   contributions, the dependence on masses other than the :math:`k`-th
   one in the coefficient functions only enters at
   :math:`\mathcal{O}(\alpha_s^3)` and therefore will be neglected here.

.. [5]
   Notice that, to uniform the notation, we understood the factor
   :math:`x` in front of :math:`F_3`.
