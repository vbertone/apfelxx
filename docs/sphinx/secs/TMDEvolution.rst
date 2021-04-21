===========================================================
Transverse momentum resummation and matching to fixed order
===========================================================

:Author: Valerio Bertone

.. contents::
   :depth: 3
..

TMD evolution and its perturbative expansion
============================================

.. container::
   :name: evolution-of-the-tmds

   .. rubric:: Evolution of the TMDs
      :name: evolution-of-the-tmds

In this section I will show how to solve the renormalisation-group
equation (RGE) and the rapidity-evolution equation (often referred to as
Collins-Soper (CS) equation) of a TMD distribution :math:`F`. The
distribution :math:`F` can be either a PDF or a FF and can be associated
to either to a quark (anti)flavour or to the gluon: the structure of the
solution of the evolution equations is exactly the same. In fact, the
solution of the equations only depends on whether we are evolving a
quark or a gluon, while it does not distinguish between PDFs and FFs. In
the impact-parameter space, :math:`F` is a function of the
transverse-momentum fraction :math:`x`, of the bidimensional
impact-parameter vector :math:`\mathbf{b}_T`, of the renormalisation
scale :math:`\mu`, and of the rapidity scale :math:`\zeta`, *i.e.*
:math:`F\equiv F(x, \mathbf{b}_T,\mu,\zeta)`. Since the evolution
equations govern the behaviour of :math:`F` w.r.t. the scale :math:`\mu`
and :math:`\zeta`, in order to simplify the notation, in this section I
will drop the dependence on :math:`x` and :math:`\mathbf{b}_T`, *i.e.*
:math:`F\equiv F(\mu,\zeta)`.( [1]_)

The solution of the evolution equations allows one to express the
distribution :math:`F` at some final scales :math:`(\mu,\zeta)` in terms
of the same distribution at the initial scales :math:`(\mu_0,\zeta_0)`.
It will turn out that this is accomplished by computing the evolution
kernel :math:`R\left[(\mu,\zeta)\leftarrow (\mu_0,\zeta_0)\right]`, such
that:

.. math::
   :label: eq:evkernel
     F(\mu,\zeta) = R\left[(\mu,\zeta)\leftarrow
       (\mu_0,\zeta_0)\right]F(\mu_0,\zeta_0)\,.

\ The evolution kernel :math:`R` can be expressed in terms of
perturbatively computable quantities. A collateral aspect that will be
discussed below is the independence from the path that connects the
initial and final scales :math:`(\mu_0,\zeta_0)` and
:math:`(\mu,\zeta)`. This in turn concerns the resummation of large
logarithms that is required to ensure that the perturbative convergence
is not spoiled.

The RGE and the CS equations read:

.. math::
   :label: eq:eveqs
   \begin{array}{l}
   \displaystyle \frac{\partial \ln F}{\partial \ln \sqrt{\zeta}} =
     K(\mu)\,,\\
   \\
   \displaystyle \frac{\partial \ln F}{\partial \ln \mu} = \gamma(\mu,\zeta)\,,
   \end{array}

\ where :math:`\gamma` and :math:`K` are the anomalous dimensions of the
evolution in :math:`\mu` and :math:`\sqrt{\zeta}`, respectively. The
equations above can be solved as follows. The first equation gives:

.. math::
   :label: eq:zetadir
   F(\mu,\zeta) = \exp\left[  K(\mu)\ln\frac{\sqrt{\zeta}}{\sqrt{\zeta_0}}\right]F(\mu,\zeta_0)\,.

The factor :math:`F(\mu,\zeta_0)` can then be evolved in :math:`\mu`
using the second equation:

.. math::
   :label: eq:mudir
   F(\mu,\zeta_0) = \exp\left[ \int_{\mu_0}^{\mu}\frac{d\mu'}{\mu'}\gamma(\mu',\zeta_0)\right]F(\mu_0,\zeta_0)\,,

such that:

.. math::
   :label: eq:solution1
     F(\mu,\zeta) = \exp\left[ K(\mu)\ln\frac{\sqrt{\zeta}}{\sqrt{\zeta_0}}+\int_{\mu_0}^{\mu}\frac{d\mu'}{\mu'}\gamma(\mu',\zeta_0)\right]F(\mu_0,\zeta_0)\,.

This equation has exactly the structure of
Eq. :eq:`eq:evkernel`. We now need to express the
argument of the exponent in terms of perturbatively computable
quantities.

In order to do so, we use the fact that the rapidity anomalous dimension
:math:`K` needs to be renormalised and thus it obeys its own RGE, that
reads:

.. math::
   :label: eq:CuspRGE
   \frac{\partial K}{\partial \ln \mu} = - \gamma_K(\alpha_s(\mu))\,.

:math:`\gamma_K` is said cusp anomalous dimension and obeys the
perturbative expansion:

.. math::
   :label: eq:GKexp
     \gamma_K(\alpha_s(\mu)) = \sum_{n=0}^{\infty}\left(\frac{\alpha_s(\mu)}{4\pi}\right)^{n+1}\gamma_K^{(n)}\,,

where :math:`\gamma_K^{(n)}` are numerical coefficients. The value of
the coefficients up to :math:`n=3` can be read from Eq. (59) of
Ref. (John Collins and Rogers 2017). They coincide with those reported
in Eq. (D.6) of Ref. (Echevarria, Scimemi, and Vladimirov 2016) up to a
factor two due to a different normalisation of the rapidity anomalous
dimension :math:`K` whose derivative w.r.t. :math:`\zeta` is exactly
:math:`\gamma_K`. In addition, the cusp anomalous dimensions for quarks
and gluon are equal up to a factor :math:`C_F` in the quark case and
:math:`C_A` in the gluon case.

Eq. :eq:`eq:CuspRGE` can be easily solved obtaining:

.. math::
   :label: eq:CuspEv
   K(\mu) = K(\mu_0) - \int_{\mu_0}^{\mu}\frac{d\mu'}{\mu'}\gamma_K(\alpha_s(\mu'))\,.

We anticipate that in the :math:`\overline{\mbox{MS}}` renormalisation
scheme, there exists a particular scale,
:math:`\mu_b = 2e^{-\gamma_E}/b_T`, such that :math:`K` computed at
:math:`\mu_b` is free of logaritms and admits the following perturbative
expansion:

.. math::
   :label: eq:Kexp
   K(\mu_b) = \sum_{n=0}^{\infty}\left(\frac{\alpha_s(\mu_b)}{4\pi}\right)^{n+1}K^{(n,0)}\,,

where :math:`K^{(n,0)}` are numerical coefficients. Therefore, if
:math:`\mu_0\simeq \mu_b` the first term in the r.h.s. of
Eq. :eq:`eq:CuspEv` is free of large logs and thus its
perturbative expansion, that reads:

.. math::
   :label: eq:rapandimlog
   K(\mu_0) = \sum_{n=0}^{\infty}\left(\frac{\alpha_s(\mu_0)}{4\pi}\right)^{n+1}\sum_{m=0}^{n+1}K^{(n,m)}\ln^m\left(\frac{\mu_0}{\mu_b}\right)\,,

is reliable. The logarithmic terms in
Eq. :eq:`eq:rapandimlog` can be computed by using
Eq. :eq:`eq:CuspRGE` or
Eq. :eq:`eq:CuspEv` and expanding :math:`\alpha_s(\mu_0)`
around :math:`\alpha_s(\mu_b)` using the appropriate RGE. Again, such an
expansion is reliable only if :math:`\mu_0` and :math:`\mu_b` are
comparable such not to generate large logarithms. The second term in the
r.h.s. of Eq. :eq:`eq:CuspEv` instead takes care, through
the evolution of :math:`\alpha_s`, of resumming large logarithms in the
case in which :math:`\mu\gg \mu_0`. The coefficients :math:`K^{(n,m)}`
up to :math:`n=2` are reported in Eq. (D.9) of Ref. (Echevarria,
Scimemi, and Vladimirov 2016) and up to :math:`n=1` in Eq. (69) of
Ref. (John Collins and Rogers 2017). They differ by a factor :math:`-2`
due to a different definition of :math:`K`. In addition, the logarithmic
expansion is done in terms of :math:`\ln(\mu_0/\mu_b)` in Ref. (John
Collins and Rogers 2017) and in terms of :math:`\ln(\mu_0^2/\mu_b^2)` in
Ref. (Echevarria, Scimemi, and Vladimirov 2016). Therefore, each
coefficient differs by an additional factor :math:`2^m`, where :math:`m`
is the power of the logarithm that multiplies the coefficient itself.

For the sake of completeness, we compute the coefficients
:math:`K^{(n,m)}` up to :math:`n=3`. To do so, we need to know the
coefficients of the RGE of the coupling :math:`\alpha_s`, that we write:

.. math:: \frac{d}{d\ln\mu}\left(\frac{\alpha_s(\mu)}{4\pi}\right) = -2 \left(\frac{\alpha_s(\mu)}{4\pi}\right)^{2}\sum_{n=0}^{\infty}\left(\frac{\alpha_s(\mu)}{4\pi}\right)^{n}\beta_n \,,

whose solution for the evolution between :math:`\mu_b` and :math:`\mu'`
expanded to two loops reads:

.. math::
   \left(\frac{\alpha_s(\mu')}{4\pi}\right) =
   \left(\frac{\alpha_s(\mu_b)}{4\pi}\right)+\left(\frac{\alpha_s(\mu_b)}{4\pi}\right)^2\left(-2\beta_0\ln\frac{\mu'}{\mu_b}\right) +\left(\frac{\alpha_s(\mu_b)}{4\pi}\right)^3\left(4\beta_0^2 \ln^2\frac{\mu'}{\mu_b}-2\beta_1 \ln\frac{\mu'}{\mu_b}\right) +\mathcal{O}(\alpha_s^4)\,,

such that:

.. math::
   \left(\frac{\alpha_s(\mu')}{4\pi}\right)^2 =
   \left(\frac{\alpha_s(\mu_b)}{4\pi}\right)^2+\left(\frac{\alpha_s(\mu_b)}{4\pi}\right)^3\left(-4\beta_0\ln\frac{\mu'}{\mu_b}\right) +\mathcal{O}(\alpha_s^4)\,,

and:

.. math::
   \left(\frac{\alpha_s(\mu')}{4\pi}\right)^3 =
   \left(\frac{\alpha_s(\mu_b)}{4\pi}\right)^3+\mathcal{O}(\alpha_s^4)\,.

Therefore, using Eq. :eq:`eq:GKexp` up to three loops, one
finds:

.. math::
   :label: eq:GKexpForInt
   \begin{array}{rcl}
     \gamma_K(\alpha_s(\mu')) &=&\displaystyle
     \left(\frac{\alpha_s(\mu_b)}{4\pi}\right) \gamma_K^{(0)}+\left(\frac{\alpha_s(\mu_b)}{4\pi}\right)^2\left[-2\beta_0 \gamma_K^{(0)}\ln\frac{\mu'}{\mu_b}+\gamma_K^{(1)}\right]
     \\
   \\
   &+&\displaystyle \left(\frac{\alpha_s(\mu_b)}{4\pi}\right)^3\left[4\beta_0^2 \gamma_K^{(0)}
       \ln^2\frac{\mu'}{\mu_b}+\left(-2\beta_1 \gamma_K^{(0)}-4\beta_0 \gamma_K^{(1)}\right)\ln\frac{\mu'}{\mu_b}+\gamma_K^{(2)}\right]
     +\mathcal{O}(\alpha_s^4)
     \,.
   \end{array}

Using Eq. :eq:`eq:GKexpForInt` in the integral in
the r.h.s. of Eq. :eq:`eq:CuspEv` computed between
:math:`\mu_b` and :math:`\mu_0`, one finds:

.. math::
   \begin{array}{rcl}
   \displaystyle \int_{\mu_b}^{\mu_0}\frac{d\mu'}{\mu'}\gamma_K(\alpha_s(\mu')) &=&
     \displaystyle \int_{0}^{\ln\frac{\mu_0}{\mu_b}}d\ln\frac{\mu'}{\mu_b}\gamma_K(\alpha_s(\mu'))\\
   \\
   &=&\displaystyle 
     \left(\frac{\alpha_s(\mu_b)}{4\pi}\right) \gamma_K^{(0)}\ln\frac{\mu_0}{\mu_b}+\left(\frac{\alpha_s(\mu_b)}{4\pi}\right)^2\left[-\beta_0 \gamma_K^{(0)}\ln^2\frac{\mu_0}{\mu_b}+\gamma_K^{(1)}\ln\frac{\mu_0}{\mu_b}\right]
     \\
   \\
   &+&\displaystyle \left(\frac{\alpha_s(\mu_b)}{4\pi}\right)^3\left[\frac{4}{3}\beta_0^2 \gamma_K^{(0)}
       \ln^3\frac{\mu_0}{\mu_b}+\left(-\beta_1 \gamma_K^{(0)}-2\beta_0 \gamma_K^{(1)}\right)\ln^2\frac{\mu_0}{\mu_b}+\gamma_K^{(2)}\ln\frac{\mu_0}{\mu_b}\right]
     +\mathcal{O}(\alpha_s^4)
     \,.
   \end{array}

Replacing :math:`\alpha_s(\mu_b)` with :math:`\alpha_s(\mu_0)` through
the expansion:

.. math::
   :label: eq:asinvexp
   \left(\frac{\alpha_s(\mu_b)}{4\pi}\right)=\left(\frac{\alpha_s(\mu_0)}{4\pi}\right)+2\beta_0\ln\frac{\mu_0}{\mu_b}\left(\frac{\alpha_s(\mu_0)}{4\pi}\right)^2+\left(2\beta_1\ln\frac{\mu_0}{\mu_b}+4\beta_0^2\ln^2\frac{\mu_0}{\mu_b}\right)\left(\frac{\alpha_s(\mu_0)}{4\pi}\right)^3+\mathcal{O}(\alpha_s^4)\,,

the final result is:

.. math::
   :label: eq:CSkernelScales
   K(\mu_0) =\sum_{n=0}^{2}\left(\frac{\alpha_s(\mu_0)}{4\pi}\right)^{n+1}\sum_{m=0}^{n+1}K^{(n,m)}\ln^m\left(\frac{\mu_0}{\mu_b}\right)\,,

where we read off:

.. math::
   \begin{array}{l}
   \displaystyle K^{(0,1)}=-\gamma_K^{(0)}\,,\\
   \\
   \displaystyle K^{(1,1)}=-\gamma_K^{(1)}+2\beta_0 K^{(0,0)}\,,\quad K^{(1,2)}=-\beta_0 \gamma_K^{(0)}\\
   \\
   \displaystyle K^{(2,1)}=-\gamma_K^{(2)}+2\beta_1 K^{(0,0)}+4\beta_0 K^{(1,0)}\,,\quad K^{(2,2)}=-\beta_1 \gamma_K^{(0)}-2\beta_0 \gamma_K^{(1)}+4\beta_0^2K^{(0,0)}\,,\quad K^{(2,3)}=-\frac{4}{3}\beta_0^2 \gamma_K^{(0)}\\
   \end{array}

A further important property of the anomalous dimensions can be derived
by considering the fact that the crossed double derivarives of :math:`F`
must be equal:

.. math::
   \frac{\partial}{\partial \ln \mu} \frac{\partial \ln F}{\partial \ln
     \sqrt{\zeta}} = \frac{\partial}{\partial \ln \sqrt{\zeta}} \frac{\partial \ln F}{\partial \ln
     \mu}\,.

\ Using Eqs. :eq:`eq:eveqs` and
:eq:`eq:CuspRGE` leads to the following additional
differential equation:

.. math::
   \frac{\partial \gamma }{\partial \ln
     \sqrt{\zeta}} = - \gamma_K(\alpha_s(\mu))\,.

Using the point :math:`\zeta = \mu^2` as a boundary condition, the
solution of the equation above is:

.. math:: \gamma(\mu,\zeta) = \gamma(\mu,\mu^2) - \gamma_K(\alpha_s(\mu))\ln\frac{\sqrt{\zeta}}{\mu}\,.

It turns out that :math:`\gamma(\mu,\mu^2)` has a purely perturbative
expansion:

.. math::
   :label: eq:GFexp
     \gamma(\mu,\mu^2) \equiv \gamma_F(\alpha_s(\mu)) = \sum_{n=0}^{\infty}\left(\frac{\alpha_s(\mu)}{4\pi}\right)^{n+1}\gamma_F^{(n)}\,,

where :math:`\gamma_F^{(n)}` are again numerical coefficients that are
given in Eq. (58) of Ref. (John Collins and Rogers 2017) and Eq. (D.7)
of Ref. (Echevarria, Scimemi, and Vladimirov 2016) (Eq. (D.8) reports
the coefficients for the gluon anomalous dimension). The two sets of
coefficients differ by a minus sign due to the different definition of
the constant (non-logarithmic) term of the RGE anomalous dimension.
Therefore:

.. math::
   :label: eq:GFEv
   \gamma(\mu,\zeta) = \gamma_F(\alpha_s(\mu)) - \gamma_K(\alpha_s(\mu))\ln\frac{\sqrt{\zeta}}{\mu}\,.

Finally, plugging Eqs. :eq:`eq:CuspEv`
and :eq:`eq:GFEv` into
Eq. :eq:`eq:solution1`, one finds:

.. math::
   :label: eq:solution2
     F(\mu,\zeta) = \exp\left\{ K(\mu_0)\ln\frac{\sqrt{\zeta}}{\sqrt{\zeta_0}}+\int_{\mu_0}^{\mu}\frac{d\mu'}{\mu'}\left[\gamma_F(\alpha_s(\mu')) - \gamma_K(\alpha_s(\mu'))\ln\frac{\sqrt{\zeta}}{\mu'}\right]\right\}F(\mu_0,\zeta_0)\,.

Comparing Eq. :eq:`eq:solution2` to
Eq. :eq:`eq:evkernel` allows one to give an explicit
expression to the evolution kernel:

.. math::
   :label: eq:evkernelexp
     R\left[(\mu,\zeta)\leftarrow
       (\mu_0,\zeta_0)\right] = \exp\left\{ K(\mu_0)\ln\frac{\sqrt{\zeta}}{\sqrt{\zeta_0}}+\int_{\mu_0}^{\mu}\frac{d\mu'}{\mu'}\left[\gamma_F(\alpha_s(\mu')) - \gamma_K(\alpha_s(\mu'))\ln\frac{\sqrt{\zeta}}{\mu'}\right]\right\}\,.

Eq. :eq:`eq:solution2` has been obtained evolving the
TMD distribution :math:`F` first in the :math:`\zeta` direction
(Eq. :eq:`eq:zetadir`) and then in the :math:`\mu`
direction (Eq. :eq:`eq:mudir`). However, it is easy to
verify that exchanging the order of the evolutions leads to the exact
same result, Eq. :eq:`eq:solution2`. In particular,
the following relation holds:

.. math::
   R\left[(\mu,\zeta)\leftarrow
       (\mu_0,\zeta)\right] R\left[(\mu_0,\zeta)\leftarrow
       (\mu_0,\zeta_0)\right] = R\left[(\mu,\zeta)\leftarrow
       (\mu,\zeta_0)\right] R\left[(\mu,\zeta_0)\leftarrow
       (\mu_0,\zeta_0)\right]=R\left[(\mu,\zeta)\leftarrow
       (\mu_0,\zeta_0)\right]\,.

\ This is a direct consequence of the independence of evolution kernel
:math:`R` in Eq. :eq:`eq:evkernelexp` from the path
:math:`\mathcal{P}` followed to connect the points
:math:`(\mu_0,\zeta_0)` to the point :math:`(\mu,\zeta)`:

.. math::
   R\left[(\mu,\zeta)\mathop{\leftarrow}_{\mathcal{P}}
       (\mu_0,\zeta_0)\right] \equiv R\left[(\mu,\zeta)\leftarrow
       (\mu_0,\zeta_0)\right]\,.

Another important piece of information comes from the fact that, for
small values of :math:`b_T`, the TMD :math:`F` can be matched onto the
respective collinear distribution :math:`f` (a PDF or a FF) through the
perturbative coefficients :math:`C`\ ( [2]_):

.. math::
   :label: eq:matching
   F(\mu,\zeta) = C(\mu,\zeta) \otimes f(\mu)\,,

\ so that:

.. math::
   :label: eq:solution3
     F(\mu,\zeta) = \exp\left\{
       K(\mu_0)\ln\frac{\sqrt{\zeta}}{\sqrt{\zeta_0}}+\int_{\mu_0}^{\mu}\frac{d\mu'}{\mu'}\left[\gamma_F(\alpha_s(\mu'))
         -
         \gamma_K(\alpha_s(\mu'))\ln\frac{\sqrt{\zeta}}{\mu'}\right]\right\}C(\mu_0,\zeta_0)\otimes
     f(\mu_0)\,.

Exactly as in the case of :math:`K`, for
:math:`\mu_0=\sqrt{\zeta_0}=\mu_b` the matching function admits the
expansion:

.. math::
   :label: eq:MatchCoeffExp
   C(\mu_b,\mu_b^2) =\sum_{n=0}^{\infty}\left(\frac{\alpha_s(\mu_b)}{4\pi}\right)^{n}C^{(n,0)}\,,

where the coefficients :math:`C^{(n,0)}` are functions of :math:`x`
only. In order to be able to compute the function :math:`C` for generic
values of the scales :math:`\mu` and :math:`\zeta`, evolution equations
can be derived. Deriving Eq. :eq:`eq:matching` with
respect to :math:`\mu` and :math:`\zeta` gives:

.. math::
   :label: eq:derivs1
   \begin{array}{l}
   \displaystyle \frac{\partial F}{\partial \ln\sqrt{\zeta}} =
     \frac{\partial C}{\partial \ln\sqrt{\zeta}}\otimes f(\mu)\,,\\
   \\
   \displaystyle \frac{\partial F}{\partial \ln\mu} = \frac{\partial
     C}{\partial \ln\mu}\otimes f(\mu) + C(\mu,\zeta) \otimes
     \frac{\partial f}{\partial \ln\mu} = \left[\frac{\partial
     C}{\partial \ln\mu}+C(\mu,\zeta) \otimes 2P(\alpha_s(\mu)) \right]\otimes f(\mu)\,.
   \end{array}

In the r.h.s. of the second equation I have used the DGLAP equation:

.. math::
   :label: eq:DGLAPeq
   \frac{\partial f}{\partial \ln\mu} = 2P(\alpha_s(\mu)) \otimes f(\mu)\,.

One can also take the derivative of
Eq. :eq:`eq:solution2` and, using
Eq. :eq:`eq:matching`, the result is:

.. math::
   :label: eq:derivs2
   \begin{array}{l}
   \displaystyle \frac{\partial F}{\partial \ln\sqrt{\zeta}} =\left[
     K(\mu_0)-\int_{\mu_0}^{\mu}\frac{d\mu'}{\mu'}\gamma_K(\alpha_s(\mu'))\right]F(\mu,\zeta)= \left[
     K(\mu_0)-\int_{\mu_0}^{\mu}\frac{d\mu'}{\mu'}\gamma_K(\alpha_s(\mu'))\right]C(\mu,\zeta)\otimes f(\mu)
   \,,\\
   \\
   \displaystyle \frac{\partial F}{\partial \ln\mu} = \left[\gamma_F(\alpha_s(\mu))
         -
         \gamma_K(\alpha_s(\mu))\ln\frac{\sqrt{\zeta}}{\mu}\right]
     F(\mu,\zeta) = \left[\gamma_F(\alpha_s(\mu))
         -
         \gamma_K(\alpha_s(\mu))\ln\frac{\sqrt{\zeta}}{\mu}\right]
     C(\mu,\zeta)\otimes f(\mu)\,.
   \end{array}

Equating Eq. :eq:`eq:derivs1` and
Eq. :eq:`eq:derivs2` and dropping the distribution
:math:`f`, the evolution equations for :math:`C` are obtained:

.. math::
   :label: eq:CEvEq
   \begin{array}{l}
     \displaystyle \frac{\partial C}{\partial \ln\sqrt{\zeta}} = \left[
     K(\mu_0)-\int_{\mu_0}^{\mu}\frac{d\mu'}{\mu'}\gamma_K(\alpha_s(\mu'))\right]C(\mu,\zeta)
     \,,\\
     \\
     \displaystyle \frac{\partial
     C}{\partial \ln\mu} = \left\{\left[\gamma_F(\alpha_s(\mu))
     -
     \gamma_K(\alpha_s(\mu))\ln\frac{\sqrt{\zeta}}{\mu}\right]\Delta-2P(\alpha_s(\mu))\right\}\otimes
     C(\mu,\zeta)\,,
   \end{array}

where the :math:`\Delta` is a matrix in flavour space whose components
are defined as:

.. math:: \Delta_{ij}(x) = \delta_{ij}\delta(1-x)\,,

being :math:`i` and :math:`j` flavour indices. The equations above can
be solved to determine the evolution of the matching function :math:`C`.
If one assumes :math:`\mu_0=\sqrt{\zeta_0}\simeq \mu_b` in
Eq. :eq:`eq:solution2`, the matching function
:math:`C` can be reliably expanded as:

.. math::
   :label: eq:MatchExp
     C(\mu_0,\mu_0^2) = \sum_{n=0}^{\infty}\left(\frac{\alpha_s(\mu_0)}{4\pi}\right)^{n}\sum_{m=0}^{2n}C^{(n,m)}\ln^m\left(\frac{\mu_0}{\mu_b}\right)\,.

The coefficient functions :math:`C^{(n,m)}` have been computed for both
PDFs and FFs in Ref. (Echevarria, Scimemi, and Vladimirov 2016). The
same functions have been computed also in Ref. (Catani et al. 2012) and
reported in Ref. (John Collins and Rogers 2017). The authors of the
latter paper have verified the equality of the two sets of functions.

For the sake of completeness, we rederive here the functions
:math:`C^{(n,m)}`, for :math:`m\neq0`, up to :math:`n=2`. This is done
by solving the equations in Eq. :eq:`eq:CEvEq` to evolve
the coefficient :math:`C` between :math:`(\mu_b,\mu_b^2)` and
:math:`(\mu_0,\zeta_0)` and, under the assumption that
:math:`\mu_0\simeq\sqrt{\zeta_0}\simeq \mu_b`, expand the solution to
:math:`\mathcal{O}(\alpha_s^2)`. Contrary to the usual assumption, we
take :math:`\mu_0\neq\sqrt{\zeta_0}` but in order to simplify the
notation we introduce the following definitions:

.. math::
   a_b = \frac{\alpha_s(\mu_b)}{4\pi}\,,\quad a_0 = \frac{\alpha_s(\mu_0)}{4\pi}\,,\quad
   L_\mu=\ln\frac{\mu_0}{\mu_b}\,,\quad L_\zeta=\ln\frac{\sqrt{\zeta_0}}{\mu_b}\,,

and we set :math:`C^{(0,0)}=\Delta\equiv 1` keeping in mind that any
scalar carries the flavour and :math:`x` structure of :math:`\Delta`.
Due to the fact the evolution kernel of the first equation in
Eq. :eq:`eq:CEvEq` is in dependent of :math:`\zeta`, the
evolution of :math:`C` in :math:`\zeta` reads:

.. math::
   :label: eq:Cexp1
   \begin{array}{rcl}
   C(\mu_b,\zeta_0) &=&\displaystyle
   \exp\left[K(\mu_b)L_{\zeta}\right]
   C(\mu_b,\mu_b^2) \\
   \\
   &=&\displaystyle
       1+a_b\left[C^{(1,0)}+K^{(0,0)}L_{\zeta}\right]+
       a_b^2\left[C^{(2,0)}+\left(K^{(0,0)}
       C^{(1,0)}+K^{(1,0)}
       \right)L_{\zeta} + \frac12\left(K^{(0,0)}\right)^2 L_{\zeta}^2\right]+\mathcal{O}(\alpha_s^3)
   \end{array}

Now we need to solve also the second equation in
Eq. :eq:`eq:CEvEq` to evolve :math:`C` from
:math:`(\mu_b,\zeta_0)` to :math:`(\mu_0,\zeta_0)`. The solution is:

.. math::
   :label: eq:Cexp2
     C(\mu_0,\zeta_0) =
                          \exp\left[\sum_{n=0}^{\infty}\int_{\mu_b}^{\mu_0}\frac{d\mu'}{\mu'}\left(\frac{\alpha_s(\mu')}{4\pi}\right)^{n+1}\left\{\left[\gamma_F^{(n)}
                          -
                          \gamma_K^{(n)}\ln\frac{\sqrt{\zeta_0}}{\mu'}\right]\Delta-2P^{(n)}\right\}\right]\otimes
                          C(\mu_b,\zeta_0) \\

\ Finally, we need to plug Eq. :eq:`eq:Cexp1` into
Eq. :eq:`eq:Cexp2`, replace :math:`\alpha_s(\mu_b)` with
:math:`\alpha_s(\mu_0)` using Eq. :eq:`eq:asinvexp`,
and truncate the result to :math:`\mathcal{O}(\alpha_s^2)`. We refrain
from writing the full expression because too lengthy that however can be
written in a compact form as:

.. math::
   :label: eq:MatchExp2
     C(\mu_0,\zeta_0) = \sum_{n=0}^{2}a_0^{n}\sum_{m=0}^{2n}\sum_{k=0}^{m}\bar{C}^{(n,m-k,k)}L_\mu^{m-k}L_\zeta^{k} +\mathcal{O}(\alpha_s^3)\,,

with:

.. math:: C^{(n,m)}=\sum_{k=0}^{m}\bar{C}^{(n,m-k,k)}\,.

In order to simplify the structure, one may take:

.. math:: L_\mu=L_\zeta=L_0\,,

and accounting for the fact that :math:`K^{(0,0)}=0`, we have:

.. math::
   \begin{array}{rcl}
     C(\mu_0,\mu_0^2) &=&\displaystyle 1+a_0\left[
                          C^{(1,0)}+\left(\gamma_F^{(0)}
                          -2P^{(0)}\right)L_0-\frac12\gamma_K^{(0)}L_0^2\right]\\
     \\
                      &+&\displaystyle a_0^2\bigg[C^{(2,0)}+\left(\gamma_F^{(1)} +K^{(1,0)}-2P^{(1)}+\left(\gamma_F^{(0)}+2\beta_0\right)C^{(1,0)}
                          -2C^{(1,0)}\otimes P^{(0)}\right)L_0
                          \\
     \\
                      &+&\displaystyle \left(\beta_0 \gamma_F^{(0)}+\frac12\left(\gamma_F^{(0)}\right)^2 -2 \left(\beta_0
                          +
                          \gamma_F^{(0)}\right) P^{(0)} +2P^{(0)}\otimes P^{(0)}-\frac12\gamma_K^{(1)}-\frac12\gamma_K^{(0)}C^{(1,0)}\right)L_0^2\\
     \\
                      &+&\displaystyle 
                          \left(-\frac12\gamma_K^{(0)}\gamma_F^{(0)}
                          +\gamma_K^{(0)} P^{(0)}-\frac23 \beta_0 \gamma_K^{(0)}\right) L_0^3+
                          \frac18\left(\gamma_K^{(0)}\right)^2L_0^4\bigg]+\mathcal{O}(\alpha_s^3)\,.
   \end{array}

In order to use Eq. :eq:`eq:solution3` in
phenomenological applications, one needs to define the values of both
the initial and final pairs of scales, :math:`(\mu_0,\zeta_0)` and
:math:`(\mu,\zeta)`. As discussed above, in the
:math:`\overline{\mbox{MS}}` renormalisation scheme the natural value
for the initial scales is
:math:`\mu_0=\sqrt{\zeta_0}=\mu_b = 2e^{-\gamma_E}/b_T`. This particular
choice is such that all logarithms of the initial scales nullify.
However, in order to estimate higer-order corrections, one can displace
the initial scales around the central value by a modest factors
:math:`C_i^\mu` and :math:`C_i^\zeta` such that
:math:`(\mu_0,\zeta_0) = (C_i^\mu\mu_b,(C_i^\zeta\mu_b)^2)`. To
implement these variations, one needs to used
Eq. :eq:`eq:CSkernelScales` for the Collins-Soper
evolution kernel :math:`K` and Eq. :eq:`eq:MatchExp2`
for the matching functions :math:`C`.

Now we discuss the final scales :math:`\mu` and :math:`\zeta`. First it
is important to stress that the choice of :math:`\zeta` in the
computation of physical observables is in fact *immaterial*. The reason
is that a physical observable that takes place at the typical hard scale
:math:`Q` is proportional to *two* TMD distributions:

.. math:: \sigma(Q)\propto H\left(Q,\mu\right) F_1(\mu,\zeta_1) F_2(\mu,\zeta_2)

where :math:`H` is the appropriate hard function and with the constraint
:math:`\zeta_1\zeta_2=Q^4`. Using the evolution kernel :math:`R`, we
find:

.. math::
   :label: eq:physobsev
   \sigma(Q)\propto H\left(Q,\mu\right) R\left[(\mu,\zeta_1)\leftarrow(\mu_0,\zeta_0)\right]R\left[(\mu,\zeta_2)\leftarrow(\mu_0,\zeta_0)\right]F_1(\mu_0,\zeta_0) F_2(\mu_0,\zeta_0)\,.

But using Eq. :eq:`eq:evkernelexp`:

.. math::
   \begin{array}{l}
     R\left[(\mu,\zeta_1)\leftarrow(\mu_0,\zeta_0)\right]R\left[(\mu,\zeta_2)\leftarrow(\mu_0,\zeta_0)\right]\\
   \\
     =\displaystyle \exp\left\{K(\mu_0)\ln\frac{\sqrt{\zeta_1
         \zeta_2}}{\zeta_0}+\int_{\mu_0}^{\mu}\frac{d\mu'}{\mu'}\left[2\gamma_F(\alpha_s(\mu'))
         -
         \gamma_K(\alpha_s(\mu'))\ln\frac{\sqrt{\zeta_1\zeta_2}}{\mu'^2}\right]\right\}\\
     \\
     =\displaystyle \exp\left\{2K(\mu_0)\ln\frac{Q}{\sqrt{\zeta_0}}+2\int_{\mu_0}^{\mu}\frac{d\mu'}{\mu'}\left[\gamma_F(\alpha_s(\mu'))
         -
         \gamma_K(\alpha_s(\mu'))\ln\frac{Q}{\mu'}\right]\right\}\,,
   \end{array}

where in the second line we have used the constraint
:math:`\zeta_1\zeta_2=Q^4`. It thus turns out that the combination of
the evolution kernels does not depend on the scales :math:`\zeta_1` and
:math:`\zeta_2`. Therefore, any choice for these scales that obeys the
constraint above leads to *exactly* the same results. This being clear,
the somewhat obvious choice is :math:`\zeta_1=\zeta_2=Q^2` and one does
not need to consider any variations.

Finally, we consider the scale :math:`\mu`. The natural choice is
:math:`\mu=Q`, however, in order to estimate higer-order corrections,
one may consider variations around this value by a modest factor
:math:`C_f`, :math:`\mu=C_fQ`. This variations generate logarithms of
:math:`\mu/Q` in the perturbative expansion of the hard factor :math:`H`
that can be computed by requiring that the physical observable
:math:`\sigma` is independent of :math:`\mu`. Specifically, one
requires:

.. math:: \frac{d\sigma}{d\ln\mu}=0\,,

\ that, using Eq. :eq:`eq:physobsev`, translates into
the following evolution equation for :math:`H`:

.. math::
   :label: eq:HEv
   \frac{d \ln H}{d \ln\mu} =-2\gamma_F(\alpha_s(\mu))
         -
         2\gamma_K(\alpha_s(\mu))\ln\frac{\mu}{Q}\,.

whose solution for the evolution between :math:`Q` and :math:`\mu` is:

.. math::
   :label: eq:HEvSol
   H(Q,\mu) = \exp\left\{-2\int_Q^\mu\frac{d\mu'}{\mu'}\left[\gamma_F(\alpha_s(\mu'))
         +
         \gamma_K(\alpha_s(\mu'))\ln\frac {\mu'}{Q}\right]\right\}H(Q,Q)\,.

Knowing that:

.. math::
   \begin{array}{rcl}
     H(Q,Q) &=&\displaystyle
                1+\sum_{n=1}^\infty\left(\frac{\alpha_s(Q)}{4\pi}\right)^nH^{(n)}\,,
   \end{array}

and assuming that :math:`Q` and :math:`\mu` are not too far apart,
Eq. :eq:`eq:HEvSol` can be reliably expanded and
truncated to :math:`\mathcal{O}(\alpha_s^2)`. Defining:

.. math:: a_\mu = \frac{\alpha_s(\mu)}{4\pi}\quad\mbox{and}\quad L_{Q} = \ln\frac{\mu}{Q}\,,

we find:

.. math::
   :label: eq:Hexp2
   \begin{array}{rcl}
     H(Q,\mu) &=&\displaystyle 1+a_\mu\left[H^{(1)}-2\gamma_F^{(0)}L_{Q}
                          -
                          \gamma_K^{(0)}L_{Q}^2\right]
                          \\
   \\
              &+&\displaystyle a_\mu^{2}\bigg[H^{(2)} + 
                  \left(-2\gamma_F^{(1)}+2\beta_0H^{(1)}-2\gamma_F^{(0)}H^{(1)}\right)L_{Q}\\
   \\
   &+&\displaystyle \left(-\gamma_K^{(1)}-2\beta_0 \gamma_F^{(0)}+2\left(\gamma_F^{(0)}\right)^2 -
                          \gamma_K^{(0)}H^{(1)}\right)L_{Q}^2\\
     \\
                      &+&\displaystyle\left(-\frac23
                          \beta_0 \gamma_K^{(0)}+ 2\gamma_F^{(0)}
                          \gamma_K^{(0)}\right)L_{Q}^3+
                          \frac12\left(\gamma_K^{(0)}\right)^2L_{Q}^4
                          \bigg]
                          +\mathcal{O}(\alpha_s^3)\,.
   \end{array}

The formula above allows us to vary the final renormalisation scale
:math:`\mu` around the natural value :math:`Q`. To summarise, there
three possible variations that allow one to estimate higher-order
corrections. The first two are associated to the variation of the
initial scales :math:`\mu_0` and :math:`\sqrt{\zeta_0}` around the
natural scale :math:`\mu_b` and they are parametrised by the modest
factors :math:`C_i^\mu` and :math:`C_i^\zeta`, such that
:math:`L_\mu=\ln C_i^\mu` and :math:`L_\zeta=\ln C_i^\zeta`. The third
variation is associated to the final renormalisation scale :math:`\mu`
around the hard scale :math:`Q` and is parameterised by the factor
:math:`C_f`, such that :math:`L_Q=\ln C_f`. In conclusion, the cross
section formula with generic scales is proportional to:

.. math::
   :label: eq:physobsev2
   \sigma(Q)\propto H\left(Q,C_fQ\right)
     \left(R\left[(C_fQ,Q^2)\leftarrow(C_i^\mu\mu_b,(C_i^\zeta\mu_b)^2)\right]\right)^2 F_1(C_i^\mu\mu_b,(C_i^\zeta\mu_b)^2)
     F_2(C_i^\mu\mu_b,(C_i^\zeta\mu_b)^2)\,.

.. container::
   :name: sec:nonpert

   .. rubric:: Non-perturbative component
      :name: sec:nonpert

In the previous section, I have considered the evolution of TMDs and
thus I concentrated on their dependence on the renormalisation and
rapidity scales :math:`\mu` and :math:`\zeta`, leaving aside the
dependence on :math:`b_T`. The computation of the rapidity evolution
kernel in the :math:`\overline{\mbox{MS}}` scheme has led to the
introduction of the scale:

.. math:: \mu_b = \frac{2e^{-\gamma_E}}{b_T}\,,

\ This scale, within a modest factor, provides a natural choice for the
evolution initial scales :math:`\mu_0` and :math:`\sqrt{\zeta_0}` that
prevents the appearance of large logarithms. Crucially, the strong
coupling :math:`\alpha_s` has to be computed in the vicinity of the
scale :math:`\mu_b`. Therefore, if the impact parameters :math:`b_T`
becomes large, :math:`\alpha_s(\mu_b)` may potentially become very large
invalidating any perturbative calculation. Since the computation of
:math:`q_T` dependent observables requires Fourier transforming TMDs,
they need to be accessed also at large values of :math:`b_T` where the
perturbative computation of the previous section breaks down. To
overcome this limitation, it is customary to introduce an *arbitrary*
scale :math:`b_{\rm max}` that denotes the maximum value of :math:`b_T`
at which one trusts perturbation theory. The value has to be such that:

.. math:: \alpha_s\left(\frac{2e^{-\gamma_E}}{b_{\rm max}}\right) \ll 1\,.

Then one introduces a monotonic function :math:`b_*` of :math:`b_T` with
the following behaviour:

.. math::
   :label: eq:bstardef
   \begin{array}{lll}
   b_*(b_T)\simeq b_T &\mbox{ for } & b_T\rightarrow 0\,,\\
   b_*(b_T) \rightarrow b_{\rm max} &\mbox{ for } & b_T\rightarrow \infty\,.
   \end{array}

A common choice is:

.. math::
   :label: eq:bstardefCSS
   b_*(b_T) = \frac{b_T}{\sqrt{1+ b_T^2/b_{\rm max}^2}}\,,

but more complicated functions can be employed to have better control on
the low-:math:`b_T` region. Now, one writes:

.. math::
   :label: eq:separatation
     F(x, {b}_T,\mu,\zeta) = \left[\frac{F(x,
         {b}_T,\mu,\zeta)}{F(x,
         b_*({b}_T),\mu,\zeta)}\right]F(x,
     b_*({b}_T),\mu,\zeta) \equiv f_{\rm NP}(x,
     {b}_T,\zeta) F(x, b_*({b}_T),\mu,\zeta)\,.

This separation, often referred to as CSS prescription after Collins,
Soper, and Sterman, is advantageous because, due to the behaviour of
:math:`b_*(b_T)`, :math:`F(x, b_*({b}_T),\mu,\zeta)` can be computed in
perturbation theory for any value of :math:`b_T`, while
:math:`f_{\rm NP}(x, {b}_T,\zeta)` embodies the non-perturbative
dependence. It is important to stress that this separation is arbitrary
and dependes of the particular choice of :math:`b_*` and
:math:`b_{\rm max}`. Therefore, for any particular choice, only the
combination in Eq. :eq:`eq:separatation` is
meaningful and it is misleading to refer to :math:`f_{\rm NP}` as to the
non-perturbative part of TMDs in a universal sense. The reason why the
function :math:`f_{\rm NP}` in
Eq. :eq:`eq:separatation` does not depend on the
renormalisation scale :math:`\mu` is that this dependence cancels out in
the ratio. To be more specific, if one chooses
:math:`\mu_0 = \mu_b = 2e^{-\gamma_E}/{b_T}` and uses
Eq. :eq:`eq:solution2`, one finds:

.. math::
   \begin{array}{l}
   \displaystyle  f_{\rm NP}(x,
     {b}_T,\zeta) = \frac{F(x,
         {b}_T,\mu,\zeta)}{F(x,
         b_*({b}_T),\mu,\zeta)} = \\
   \\
   \displaystyle \exp\left\{
       K(\mu_b)\ln\frac{\sqrt{\zeta}}{\sqrt{\zeta_b}}-
     K(\mu_{b_*})\ln\frac{\sqrt{\zeta}}{\sqrt{\zeta_{b_*}}}+\int_{\mu_b}^{\mu_{b_*}}\frac{d\mu'}{\mu'}\left[\gamma_F(\alpha_s(\mu'))
     -
     \gamma_K(\alpha_s(\mu'))\ln\frac{\sqrt{\zeta}}{\mu'}\right]\right\}\frac{F(\mu_b,\zeta_b)}{F(\mu_{b_*},\zeta_{b_*})}\,,
   \end{array}

with :math:`\mu_{b_*} = \sqrt{\zeta_{b_*}}=2e^{-\gamma_E}/{b_*(b_T)}`.
It is thus apparent that the dependence on :math:`\mu` cancels. Due to
the required behaviour of :math:`b_*`, if :math:`b_T` becomes small,
:math:`b_T` and :math:`b_*` get closer and so do :math:`\mu_{b_*}` and
:math:`\mu_b` as well as :math:`\zeta_{b_*}` and :math:`\zeta_b`.
Therefore, from the equation above, we see that
:math:`f_{\rm NP}\rightarrow 1` for :math:`b_T\rightarrow 0`.
Conversely, if :math:`b_T` becomes large :math:`b_*` saturates to
:math:`b_{\rm max}`. In this limit, we can write:

.. math::
   :label: eq:fNPexplicit
   \begin{array}{l}
   \displaystyle  f_{\rm NP}(x,
     {b}_T,\zeta) \mathop{\longrightarrow}_{b_T\rightarrow \infty}\\
   \\
   \displaystyle \exp\left\{
       K(\mu_b)\ln\frac{\sqrt{\zeta}}{\sqrt{\zeta_b}}-
     K(\mu_{\rm min})\ln\frac{\sqrt{\zeta}}{\sqrt{\zeta_{\rm min}}}+\int_{\mu _{b}}^{\mu_{\rm min}}\frac{d\mu'}{\mu'}\left[\gamma_F(\alpha_s(\mu'))
     -
     \gamma_K(\alpha_s(\mu'))\ln\frac{\sqrt{\zeta}}{\mu'}\right]
     \right\}\frac{F(\mu_b,\zeta_b)}{F(\mu_{\rm min},\zeta_{\rm min})}\,,
   \end{array}

where I have defined
:math:`\mu_{\rm min} = \sqrt{\zeta_{\rm min}} = 2e^{-\gamma_E}/b_{\rm max}`.
As :math:`b_T` becomes large, :math:`\mu_b` and :math:`\zeta_b` become
increasingly small and the exponential takes the form of an evolution
kernel between two sets of scales far apart from each other. As well
known, the evolution kernel suppresses the distribution for large final
scales. Extrapolating, one would then expect :math:`f_{\rm NP}` to be
exponentially suppressed as :math:`b_T` becomes large.

The definition of :math:`f_{\rm NP}` in
Eq. :eq:`eq:fNPexplicit` allows for a more
systematic study of this quantity that may help adopting an optimal
parameterisation when fitting it to data. The idea is to plot the ratio
in the r.h.s. of the first line of
Eq. :eq:`eq:fNPexplicit` as a function of
:math:`b_T` for fixed values of :math:`x` and :math:`\zeta`. We know
that for :math:`b_T` larger than some value, the numerator of this ratio
will become unreliable. To identify this value, one can change the
definition of :math:`\mu_b` (as well as that for :math:`\mu_{b_*}`) by
introducing a factor :math:`C_2`:

.. math:: \mu_b = C_2\frac{2e^{-\gamma_E}}{b_T}\,.

\ and vary :math:`C_2` around one by, say, a factor two up and down.
When the effect of varying :math:`C_2` on the ratio becomes large, one
can say that non-perturbative effects are large. On the other hand, in
the region of :math:`b_T` where variations are small perturbation theory
is still reliable. Therefore this region can be used to put some
constraint on :math:`f_{\rm NP}`.

.. container::
   :name: alternative-to-the-css-prescription

   .. rubric:: Alternative to the CSS prescription
      :name: alternative-to-the-css-prescription

As discussed above, the CSS prescription has the purpose to avoid
integrating the Landau pole both in the strong coupling :math:`\alpha_s`
and in the collinear distributions. This is implemented by *globally*
replacing :math:`b_T` with :math:`b_*` in
Eq. :eq:`eq:bstardefCSS` in the TMDs that in turn
naturally leads to the introduction of the non-perturbative function
:math:`f_{\rm NP}` through
Eq. :eq:`eq:separatation`. Despite very simple and
natural, this approach implies using :math:`b_*` even where not strictly
required. Specifically in the logarithms appearing explicitly in the
exponent of the Sudakov form factor that will now depend on :math:`b_*`
rather than on :math:`b_T` even though they are not problematic. Of
course, the presence of :math:`f_{\rm NP}` will reabsorb these effects
as clear from Eq. :eq:`eq:separatation`.

However, one can decide take a less invasive perspective to leave the
explicit logarithms in the Sudakov unchanged while only regularising
only :math:`\alpha_s` and the collinear distributions. This can be done
by computing them at some modified scale :math:`\mu_*(mu)` su that:

.. math::
   \begin{array}{lll}
   \mu_*(\mu)\simeq \mu &\mbox{ for } & \mu\rightarrow \infty\,,\\
   \mu_*(\mu) \rightarrow \mu_{\rm min} &\mbox{ for } & \mu\rightarrow 0\,.
   \end{array}

\ This can be easily connected to the :math:`b_*` prescription as:

.. math:: \mu_*(\mu)=\frac{b_0}{b_*\left(\frac{b_0}{\mu}\right)}\,,

with :math:`b_0=2e^{-\gamma_E}`. Adopting the particular :math:`b_*` in
Eq. :eq:`eq:bstardefCSS` and choosing
:math:`b_{\rm max}=b_0` gives the surprisingly simple relation:

.. math:: \mu_*(\mu)=\sqrt{\mu^2+1}\,.

Therefore, the alternative recipe to the CSS procedure amounts to
computing :math:`\alpha_s` and collinear distributions in
:math:`\sqrt{\mu^2+1}` rather than in :math:`\mu`. This suffices to
avoid integrating over the Landau pole, regularising only what actually
needs to be regularised and leaving the rest untouched. The one
unpleasant drawback of this approach is the fact that it does not make
as transparent as in CSS the introduction of the function
:math:`f_{\rm NP}` that in any case has to be present to encapsulate
non-perturbative effects.

The difference between the standard CSS prescription and its alternative
discussed here is displayed in Fig. `1 <#fig:NewPrescription_qT>`__ for
the :math:`q_T` distribution of the Drell-Yan cross section at
:math:`\sqrt{s} = 13` TeV, :math:`Q=M_Z`, and :math:`y=0`.

.. container::
   :name: fig:NewPrescription_qT

   .. figure:: ../../latex/src/../../latex/src/plots/NewPrescription_qT.pdf
   .. :alt: .[fig:NewPrescription_qT]
      :name: fig:NewPrescription_qT
      :width: 60.0%

      .[fig:NewPrescription_qT]

.. container::
   :name: extraction-of-the-singular-behaviour-for-q_trightarrow-0

   .. rubric:: Extraction of the singular behaviour for
      :math:`q_T\rightarrow 0`
      :name: extraction-of-the-singular-behaviour-for-q_trightarrow-0

In this section we derive the expansion of the different ingredients of
TMDs in order to finally extract the singular behaviour of a cross
section involving TMDs for :math:`q_T\rightarrow 0`. This will
eventually serve to match the low-:math:`q_T` resummed calculation to
the fixed-order one valid at :math:`q_T\simeq Q`.

.. container::
   :name: expansion-of-the-evolution-kernel

   .. rubric:: Expansion of the evolution kernel
      :name: expansion-of-the-evolution-kernel

In this subsection, we work out the perturbative expansion of the
evolution kernel :math:`R` up to :math:`\mathcal{O}(\alpha_s^2)`. This
expansion is useful to extract the (singular) asymptotic behaviour of a
cross section computed at fixed order in pQCD as :math:`q_T` tends to
zero. Since TMDs always appear in pairs in the computation of a cross
section and that the evolution kernel is universal, it is convenient to
compute the expansion of :math:`R^2`. To do so, we start from
Eq. :eq:`eq:evkernelexp` where we set
:math:`\mu_0=\sqrt{\zeta_0}=\mu_b=2e^{-\gamma_E}/{b_T}` and
:math:`\mu=\sqrt{\zeta}=Q`, possible scale variations can be reinstated
at a later stage:

.. math:: R^2 = \exp\left\{ K(\mu_b)\ln\frac{Q^2}{\mu_b^2}+\int_{\mu_b^2}^{Q^2}\frac{d\mu'^2}{\mu'^2}\left[\gamma_F(\alpha_s(\mu')) - \frac12\gamma_K(\alpha_s(\mu'))\ln\frac{Q^2}{\mu'^2}\right]\right\}\,.

Now we use the expansions in Eqs. :eq:`eq:GKexp`,
:eq:`eq:Kexp` and :eq:`eq:GFexp`:

.. math::
   :label: eq:Sudexp
     R^2 = \exp\left\{\sum_{n=0}^\infty a_s^{n+1}(\mu_b)K^{(n,0)}\ln\frac{Q^2}{\mu_b^2}+\int_{\mu_b^2}^{Q^2}\frac{d\mu'^2}{\mu'^2}a_s^{n+1}(\mu')\left(\gamma_F^{(n)} - \frac12\gamma_K^{(n)}\ln\frac{Q^2}{\mu'^2}\right)\right\}\,,

where we have defined:

.. math:: a_s = \frac{\alpha_s}{4\pi}\,.

In order to carry out the perturbative expansion we need to write the
argument of the exponential in terms of a common value of
:math:`\alpha_s` computed at some hard scale: the natural choice in view
of the matching is :math:`\alpha_s(Q)`. This can be achieved by using
the RGE for :math:`a_s` at leading order:

.. math::
   :label: eq:RGEas
   \mu^2\frac{da_s}{d\mu^2}=-\beta_0a_s^2(\mu)\,,

whose solution is:

.. math::
   :label: eq:asexp
   a_s(\mu) = \frac{a_s(Q)}{1+a_s(Q)\beta_0\ln(\mu^2/Q^2)}\simeq a_s(Q)\left[1+a_s(Q)\beta_0\ln(Q^2/\mu^2)+\mathcal{O}(a_s^2)\right]\,.

Plugging the expansion in the r.h.s. of Eq. :eq:`eq:asexp`
into Eq. :eq:`eq:Sudexp` and retaining only terms up to
:math:`\mathcal{O}(a_s^2)`, one finds:

.. math::
   :label: eq:SudExp1
   \begin{array}{rcl}
     R^2 &=&\displaystyle \exp\Bigg\{ a_s(Q)\left[
             K^{(0,0)}\ln\frac{Q^2}{\mu_b^2}+\int_{\mu_b^2}^{Q^2}\frac{d\mu'^2}{\mu'^2}\left(\gamma_F^{(0)}
             - \frac12\gamma_K^{(0)}\ln\frac{Q^2}{\mu'^2}\right)\right]\\
   \\
   &+& \displaystyle a_s^2(Q)\beta_0\left[K^{(0,0)}\ln^2\frac{Q^2}{\mu_b^2}+\int_{\mu_b^2}^{Q^2}\frac{d\mu'^2}{\mu'^2}\left(\gamma_F^{(0)}\ln\frac{Q^2}{\mu'^2}
             - \frac12\gamma_K^{(0)}\ln^2\frac{Q^2}{\mu'^2}\right)\right]\\
   \\
   &+&\displaystyle a_s^2(Q)\left[K^{(1,0)}\ln\frac{Q^2}{\mu_b^2}+\int_{\mu_b^2}^{Q^2}\frac{d\mu'^2}{\mu'^2}\left(\gamma_F^{(1)}
             - \frac12\gamma_K^{(1)}\ln\frac{Q^2}{\mu'^2}\right)\right]+\mathcal{O}(a_s^3)\Bigg\}\,.
   \end{array}

The final step before carrying out the expansion is computing the
integrals:

.. math::
   \int_{\mu_b^2}^{Q^2}\frac{d\mu'^2}{\mu'^2}\ln^k\left(\frac{Q^2}{\mu'^2}\right)
   = \int_{\ln(\mu_b^2)}^{\ln
     Q^2}d\ln\mu'^2\ln^k\left(\frac{Q^2}{\mu'^2}\right) = \int^{\ln(Q^2/\mu_b^2)}_{0}dx\,x^k=\frac{1}{k+1}\ln^{k+1}\left(\frac{Q^2}{\mu_b^2}\right)\,.

so that:

.. math::
   :label: eq:SudExp0
   \begin{array}{rcl}
     R^2 &\simeq&\displaystyle \exp\Bigg\{ a_s(Q)\left[
             \left(K^{(0,0)}+\gamma_F^{(0)}\right)L - \frac14\gamma_K^{(0)}L^2
             \right]\\
   \\
   &+& \displaystyle a_s^2(Q)\bigg[\left(K^{(1,0)}+\gamma_F^{(1)}\right)L +\left(\beta_0 K^{(0,0)}+\frac12\beta_0\gamma_F^{(0)} - \frac14\gamma_K^{(1)}\right)L^2- \frac16 \beta_0\gamma_K^{(0)}L^3\bigg]\Bigg\}\,,
   \end{array}

where we have defined:

.. math::
   :label: eq:reslogs
   L\equiv \ln\frac{Q^2}{\mu_b^2}\,,

Eq. :eq:`eq:SudExp1` can be conveniently written as:

.. math::
   :label: eq:SudCompact
     R^2=\exp\left\{\sum_{n=1}^2a_s^n(Q)\sum_{k=1}^{n+1}S^{(n,k)}L^k\right\}\,,

with:

.. math::
   \begin{array}{l}
   \displaystyle \displaystyle S^{(1,1)} =K^{(0,0)}+\gamma_F^{(0)}\,,\quad\displaystyle
                                                      S^{(1,2)} =- \frac14\gamma_K^{(0)}\,,\\
   \\
   \displaystyle S^{(2,1)} =K^{(1,0)}+\gamma_F^{(1)}\,,\quad\displaystyle S^{(2,2)} = \beta_0 K^{(0,0)}+\frac12\beta_0\gamma_F^{(0)} - \frac14\gamma_K^{(1)}\,,\quad\displaystyle S^{(2,3)} = - \frac16 \beta_0\gamma_K^{(0)}\,.
   \end{array}

Eq. :eq:`eq:SudCompact` can be easily expanded up to
order :math:`a_s^2` as:

.. math::
   :label: eq:exp1
   \begin{array}{rcl}
     \displaystyle R^2&=&\displaystyle
                          1+a_s(Q)\sum_{k=1}^{2}S^{(1,k)}L^k+a_s^2(Q)\left[\sum_{k=1}^{3}S^{(2,k)}L^k+\frac12\left(\sum_{k=1}^{2}S^{(1,k)}L^k\right)^2\right]+\mathcal{O}(a_s^3)\\
     \\
                      &=&\displaystyle
                          1+a_s(Q)\sum_{k=1}^{2}S^{(1,k)}L^k+a_s^2(Q)
                          \sum_{k=1}^{4}\widetilde{S}^{(2,k)}L^k
                          +\mathcal{O}(a_s^3)\,,
   \end{array}

with:

.. math::
   \begin{array}{ll}
   \displaystyle    \widetilde{S}^{(2,1)}={S}^{(2,1)}\,,&\displaystyle \quad \widetilde{S}^{(2,2)}={S}^{(2,2)}+\frac12\left[S^{(1,1)}\right]^2\,,\\
   \\
   \displaystyle
     \widetilde{S}^{(2,3)}={S}^{(2,3)}+{S}^{(1,1)}{S}^{(1,2)}\,,&\displaystyle
                                                                   \quad \widetilde{S}^{(2,4)}=\frac12\left[S^{(1,2)}\right]^2\,.
   \end{array}

.. container::
   :name: expansion-of-the-dglap-evolution

   .. rubric:: Expansion of the DGLAP evolution
      :name: expansion-of-the-dglap-evolution

A further step towards the extraction of the singular behaviour of the
resummed cross section in the limit :math:`q_T\rightarrow 0` requires
expanding the solution of the DGLAP equation that also resums logarithms
as that in Eq. :eq:`eq:reslogs`. The solution of the
DGLAP equation in Eq. :eq:`eq:DGLAPeq` that evolves PDFs
from the scale :math:`\mu_0=\mu_b` to the scale :math:`Q` can be written
as( [3]_):

.. math:: f(Q) = \Gamma(Q,\mu_b)\otimes f(\mu_b)\,,

\ where the evolution operator :math:`\Gamma` obeys the equation:

.. math::
   :label: eq:DGLAPeqOp
   \frac{\partial \Gamma(Q,\mu_b)}{\partial \ln Q^2} = P(Q) \otimes \Gamma(Q,\mu_b)\,,

with the splitting functions :math:`P` having the perturbative
expansion:

.. math:: P(Q) =\sum_{n=0}^{\infty}a_s^{n+1}(Q) P^{(n)}\,.

The coefficients :math:`P^{(n)}` are functions of :math:`x` (or
:math:`z`) only and thus (in :math:`\overline{\mbox{MS}}`) :math:`P`
depends on :math:`Q` only through the coupling :math:`a_s`. We now need
to expand :math:`\Gamma` in powers of :math:`a_s` up to
:math:`\mathcal{O}(a_s^2)`. To this end, we assume that :math:`\Gamma`
obeys the expansion:

.. math:: \Gamma(Q,\mu_b) = \Delta + \sum_{n=1}^{\infty} a_s^{n}(Q)\sum_{k=1}^{n}L^k\Gamma^{(n,k)}\,,

where the coefficients :math:`\Gamma^{(n,k)}` are functions of :math:`x`
only that we need to determine up to :math:`n=2`.
Eq. :eq:`eq:DGLAPeqOp` can be solved iteratively. At
:math:`\mathcal{O}(a_s)`, considering that:

.. math:: \frac{da_s}{d \ln Q^2}=\mathcal{O}(a_s^2)\,,

this gives:

.. math:: \Gamma^{(1,1)}=P^{(0)}\,.

At :math:`\mathcal{O}(a_s^2)`, this gives:

.. math::
   -L\beta_0 P^{(0)} +\Gamma^{(2,1)}+2L\Gamma^{(2,2)}= P^{(1)}+L
   P^{(0)}\otimes P^{(0)}\,,

so that:

.. math::
   \begin{array}{l}
   \displaystyle \Gamma^{(2,1)} = P^{(1)}\,,\\
   \\
   \displaystyle \Gamma^{(2,2)} = \frac12\beta_0 P^{(0)} + \frac12 P^{(0)}\otimes P^{(0)}\,.
   \end{array}

Putting all pieces together, we find:

.. math::
   \Gamma(Q,\mu_b) = \Delta + a_s(Q)L P^{(0)}
   +a_s^2(Q)\left[LP^{(1)}+ \frac{L^2}2 \left(\beta_0 P^{(0)} +P^{(0)}\otimes P^{(0)}\right)\right]+\mathcal{O}(a_s^3)\,.

However, the quantity we are interested in is the inverse of the
operator :math:`\Gamma`, that evolves the PDFs from the scale :math:`Q`
to :math:`\mu_b`:

.. math::
   :label: eq:DGLAPeqOp1
   \begin{array}{l}
   \displaystyle \Gamma(\mu_b,Q)=\Gamma^{-1}(Q,\mu_b) \\
   \\\displaystyle
    = \Delta - a_s(Q)L P^{(0)}
   - a_s^2(Q)\left[LP^{(1)}+ \frac{L^2}2 \left(\beta_0 P^{(0)}-P^{(0)}\otimes P^{(0)}\right)\right]+\mathcal{O}(a_s^3)\,.
   \end{array}

.. container::
   :name: expansion-of-the-matching-coefficients-at-the-hard-scale

   .. rubric:: Expansion of the matching coefficients at the hard scale
      :name: expansion-of-the-matching-coefficients-at-the-hard-scale

The final step is the expansion of the matching coefficients :math:`C`
defined in Eq. :eq:`eq:matching` for
:math:`\mu_0=\mu_b` around the scale :math:`Q` up to order
:math:`a_s^2`. To do this, we use the expansion in
Eq. :eq:`eq:MatchCoeffExp`, in which :math:`a_s`
is computed at :math:`\mu_b`, and combine it with
Eq. :eq:`eq:asexp`. This yields:

.. math::
   :label: eq:MatchCoeffExp1
     C(\mu_b,\mu_b^2) = \Delta + a_s(Q) C^{(1,0)} + a_s^2(Q) \left[C^{(2,0)}+L\beta_0C^{(1,0)}\right]+\mathcal{O}(a_s^3)\,.

.. container::
   :name: expansion-of-the-single-low-scale-tmd

   .. rubric:: Expansion of the single low-scale TMD
      :name: expansion-of-the-single-low-scale-tmd

Before going into the more convoluted (but more relevant) case of a
cross section, it is useful to compute the expansion of a single
low-scale TMD, *i.e.* the convolution between matching coefficients and
collinear distributions at the scale :math:`\mu_0=\mu_b`. This can be
achieved combining Eqs. :eq:`eq:DGLAPeqOp1`
and :eq:`eq:MatchCoeffExp1`:

.. math::
   :label: eq:LowScaleTMD
   \begin{array}{rcl}
   C(\mu_b,\mu_b^2)\otimes f(\mu_b) &=& C(\mu_b,\mu_b^2)\otimes
   \Gamma(\mu_b,Q) \otimes f(Q)\\
   \\
   &=& \bigg\{\Delta + a_s(Q)\left[C^{(1,0)}-L P^{(0)}\right]\\
   \\
   &+&\displaystyle a_s^2(Q)\bigg[C^{(2,0)}+L\left(-P^{(1)}-C^{(1,0)}\otimes P^{(0)}+\beta_0C^{(1,0)}\right)\\
   \\
   &+&\displaystyle \frac{L^2}2 \left(P^{(0)}\otimes P^{(0)}-\beta_0 P^{(0)}\right)\bigg]\bigg\}\otimes f(Q)+\mathcal{O}(a_s^3)\,.
   \end{array}

.. container::
   :name: expansion-of-the-cross-section

   .. rubric:: Expansion of the cross section
      :name: expansion-of-the-cross-section

We are finally in the position to gather all pieces and write down the
perturbative expansion of a cross section involving TMDs, relevant in
the limit :math:`q_T\rightarrow 0`. The only additional information
required is the hard factor :math:`H` that also admits a perturbative
expansion. It is now opportune to reinstate the flavour indices. In
general, the hard factor :math:`H` has a flavour structure, meaning that
it carries a pair of flavour indices. Choosing :math:`\mu = Q`, the
expansion of the hard factor reads:

.. math::
   :label: eq:HardFactExp
   H_{ij}(Q) = \sum_{n=0}^\infty a_s^n(Q) H_{ij}^{(n)}=H_{ij}^{(0)}+a_s(Q) H_{ij}^{(1)}+a_s^2(Q) H_{ij}^{(2)}+\mathcal{O}(a_s^3)\,,

where :math:`H_{ij}^{(n)}` are numerical factors. In relevant processes,
such as Drell-Yan and SIDIS, it turns out that up to
:math:`\mathcal{O}(\alpha_s^2)` the coefficients of the expansion of
hard factor :math:`H` are diagonal in flavour. Specifically, for scales
sufficiently below the :math:`Z`-mass scale, we have:

.. math:: H_{ij}^{(n)}\propto e_i^2\delta_{ij}\,,\quad\mbox{for}\quad n=0,1,2\,,

where :math:`e_i` is the electric charge of the :math:`i`-th quark
flavour. For higher scales, one just needs to replace the electric
charges with the effective electro-weak charges. This simplification
does not hold beyond :math:`\mathcal{O}(\alpha_s^2)`.

A physical cross section differential in :math:`q_T` can finally be
computed as:

.. math:: \displaystyle \frac{d\sigma}{dq_T^2} \propto \sum_{ij}H_{ij} (Q) \int d^2\mathbf{b}_T e^{i \mathbf{b}_T\cdot \mathbf{q}_T}F_i(x, \mathbf{b}_T,Q,Q^2) D_j(z, \mathbf{b}_T,Q,Q^2)

where :math:`F_i` and :math:`D_j` are two different sets of TMD
distributions (*e.g.* a set of PDFs and a set of FFs) indexed by a
flavour index. Assuming that the flavour indices :math:`i` and :math:`j`
run *either* on quarks *or* on the gluon and writing explicitly all
perturbative factors, one has:

.. math::
   :label: eq:xsecexp
   \begin{array}{rcl}
   \displaystyle \frac{d\sigma}{dq_T^2} &\propto&\displaystyle  \sum_{ij}H_{ij} (Q) \int \frac{d^2\mathbf{b}_T}{4\pi} e^{i \mathbf{b}_T\cdot
     \mathbf{q}_T} R^2\left[(\mu_b,\mu_b^2)\rightarrow
     (Q,Q^2)\right] \\
   \\
   &\times&\displaystyle \sum_{k}\left[ \mathcal{C}_{ik}(x;\mu_b,\mu_b^2)\otimes
     f_k(x,\mu_b)\right]\sum_{l}\left[
     \mathbb{C}_{jl}(z;\mu_b,\mu_b^2)\otimes
     d_l(z,\mu_b)\right]\\
   \\
   &=&\displaystyle \sum_{n=0}^\infty a_s^n(Q)\sum_{p=0}^{2n} I_p(q_T,Q)
       \sum_{ij,kl} B_{ij,kl}^{(n,p)}(x,z) \mathop{\otimes}_x f_k(x,Q) \mathop{\otimes}_z d_l(z,Q)\,.
   \end{array}

In the equation above we used different symbols to distinguish the
perturbative components of the distributions :math:`F` and :math:`D`
because in general they are different. In addition, we have defined:

.. math::
   :label: eq:IkPert
   I_p(q_T,Q)\equiv \int \frac{d^2\mathbf{b}_T}{4\pi} e^{i \mathbf{b}_T\cdot
     \mathbf{q}_T} \ln^p\left(\frac{b_T^2Q^2}{4e^{-2\gamma_E}}\right) = \frac12\int_0^\infty db_T\,b_T J_0(b_Tq_T) \ln^p\left(\frac{b_T^2Q^2}{4e^{-2\gamma_E}}\right)\,.

Results for :math:`I_p` have been computed up to :math:`p=4` in
Eq. (136) of Appendix B of Ref. (Bozzi et al. 2006). Specifically, and
including the trivial tranform with :math:`p=0`, they read:

.. math::
   :label: eq:PureLogs
   \begin{array}{l}
   \displaystyle I_0(q_T,Q) = \delta(q_T)\,,\\
   \\
   \displaystyle I_1(q_T,Q) = - \frac{1}{q_T^2}\,,\\
   \\
   \displaystyle I_2(q_T,Q) = - \frac{2}{q_T^2}\ln\left(\frac{Q^2}{q_T^2}\right) \,,\\
   \\
   \displaystyle I_3(q_T,Q) = - \frac{3}{q_T^2}\ln^2\left(\frac{Q^2}{q_T^2}\right) \,,\\
   \\
   \displaystyle I_4(q_T,Q) = - \frac{4}{q_T^2}\left[\ln^3\left(\frac{Q^2}{q_T^2}\right)-4\zeta_3\right]\,.
   \end{array}

The functions above are compute under the assumption that there is no
non-perturbative component. Upon this assumption, that we will relax
below, all terms with :math:`p=0` will be proportional to
:math:`\delta(q_T)`. Analogous terms are *not* included in the
fixed-order calculation. Therefore, when matching the resummed
calculation to the fixed-order one, one should remove from
Eq. :eq:`eq:xsecexp` all terms with :math:`p=0`.
Importantly, this means removing the full :math:`\mathcal{O}(1)` terms
such that the leading-order term for :math:`q_T>0` is
:math:`\mathcal{O}(a_s)`.

The coefficients :math:`B^{(n,k)}` in
Eq. :eq:`eq:xsecexp` can be determined by using the
expansions worked out above in Eqs. :eq:`eq:exp1`,
:eq:`eq:LowScaleTMD`,
and :eq:`eq:HardFactExp`. The leading order is
obtained by just retaining the :math:`\mathcal{O}(a_s^0)` terms in all
the expansions. This gives:

.. math:: B_{ij,kl}^{(0,0)}(x,z) = H_{ij}^{(0)}\delta_{ik}\delta_{jl}\delta(1-x)\delta(1-z)\,.

The :math:`\mathcal{O}(a_s)` terms are obtained by combining the
leading-order ones of all expansions but one. Finally organising the
terms in powers of :math:`L`, this yields:

.. math::
   \begin{array}{rcl}
   B_{ij,kl}^{(1,0)}(x,z) &=& \displaystyle H_{ij}^{(1)}\delta_{ik}\delta(1-x) \delta_{jl}\delta(1-z)+H_{ij}^{(0)}\delta_{ik}\delta(1-x)\mathbb{C}_{jl}^{(1,0)}(z)+H_{ij}^{(0)}\mathcal{C}_{ik}^{(1,0)}(x)\delta_{jl}\delta(1-z)\,,\\
   \\
   B_{ij,kl}^{(1,1)}(x,z) &=& \displaystyle H_{ij}^{(0)} S^{(1,1)}\delta_{ik}\delta(1-x) \delta_{jl}\delta(1-z)-H_{ij}^{(0)}\delta_{ik}\delta(1-x)\mathbb{P}_{jl}^{(0)}(z)-H_{ij}^{(0)}\mathcal{P}_{ik}^{(0)}(x)\delta_{jl}\delta(1-z)\,,\\
   \\
   B_{ij,kl}^{(1,2)}(x,z) &=& \displaystyle H_{ij}^{(0)} S^{(1,2)}\delta_{ik}\delta(1-x) \delta_{jl}\delta(1-z)\,.
   \end{array}

Finally, the :math:`\mathcal{O}(a_s^2)` terms are more convoluted and
read:

.. math::
   \begin{array}{rcl}
     {B}_{ij,kl}^{(2,0)}(x,z) &=&\displaystyle H_{ij}^{(2)}
                                  \delta_{ik}\delta(1-x) \delta_{jl}\delta(1-z) \,,\\
     \\
                              &+&\displaystyle
                                  H_{ij}^{(1)}\left[\mathcal{C}_{ik}^{(1,0)}(x)\delta_{jl}\delta(1-z)
                                  +
                                  \delta_{ik}\delta(1-x)\mathbb{C}_{jl}^{(1,0)}(z)\right]\\
     \\
                              &+&\displaystyle H_{ij}^{(0)}\mathcal{C}_{ik}^{(2,0)}(x)\delta_{jl}\delta(1-z) + 
                                  H_{ij}^{(0)}\mathcal{C}_{ik}^{(1,0)}(x)\mathbb{C}_{jl}^{(1,0)}(z)
                                  +H_{ij}^{(0)}\delta_{ik}\delta(1-x)\mathbb{C}_{jl}^{(2,0)}(z)\,,\\
     \\
     {B}_{ij,kl}^{(2,1)}(x,z) &=&\displaystyle \left(H_{ij}^{(0)}\widetilde{S}^{(2,1)}+H_{ij}^{(1)}S^{(1,1)}\right)
                                 \delta_{ik}\delta(1-x) \delta_{jl}\delta(1-z)\\
     \\
                             &-& \displaystyle
                                 H_{ij}^{(1)}\left[\mathcal{P}_{ik}^{(0)}(x)\delta_{jl}\delta(1-z)
                                 +
                                 \delta_{ik}\delta(1-x)\mathbb{P}_{jl}^{(0)}(z)\right]\\
     \\
                             &+&
                                 H_{ij}^{(0)} S^{(1,1)}\left[\mathcal{C}_{ik}^{(1,0)}(x)\delta_{jl}\delta(1-z)
                                 +
                                 \delta_{ik}\delta(1-x)\mathbb{C}_{jl}^{(1,0)}(z)\right]\\
     \\
                             &+&\displaystyle 
                                 H_{ij}^{(0)}\left(-\mathcal{P}^{(1)}(x)-\sum_\alpha
                                 \mathcal{C}_{i\alpha}^{(1,0)}(x)\otimes
                                 \mathcal{P}_{\alpha k}^{(0)}(x)+\beta_0\mathcal{C}_{ik}^{(1,0)}(x)\right)\delta_{jl}\delta(1-z)\\
   \\
    &-&\displaystyle H_{ij}^{(0)}\mathcal{P}_{ik}^{(0)}(x)\mathbb{C}_{jl}^{(1,0)}(z)
        -H_{ij}^{(0)}\mathcal{C}_{ik}^{(1,0)}(x)\mathbb{P}_{jl}^{(0)}(z)\\
   \\
                                 &+&\displaystyle H_{ij}^{(0)}\delta_{jl}\delta(1-z) \left(-\mathbb{P}^{(1)}(z)-\sum_\beta
                                 \mathbb{C}_{j\beta}^{(1,0)}(z)\otimes
                                 \mathbb{P}_{\beta l}^{(0)}(z)+\beta_0\mathbb{C}_{jl}^{(1,0)}(z)\right) \,,\\
     \\
     {B}_{ij,kl}^{(2,2)}(x,z) &=&\displaystyle \left( H_{ij}^{(0)}\widetilde{S}^{(2,2)}+H_{ij}^{(1)}S^{(1,2)}\right)
                                 \delta_{ik}\delta(1-x) \delta_{jl}\delta(1-z)\\
     \\
                             &+&\displaystyle
                                 H_{ij}^{(0)}{S}^{(1,2)}\left[\mathcal{C}_{ik}^{(1,0)}(x)\delta_{jl}\delta(1-z)+\delta_{ik}\delta(1-x)\mathbb{C}_{jl}^{(1,0)}(z)\right]\\
     \\
                             &-&\displaystyle
                                 H_{ij}^{(0)} S^{(1,1)}\left[\mathcal{P}_{ik}^{(0)}(x)\delta_{jl}\delta(1-z)+\delta_{ik}\delta(1-x)\mathbb{P}_{jl}^{(0)}(z)\right]\\
     \\
    &+& \displaystyle
        H_{ij}^{(0)}\mathcal{P}_{ik}^{(0)}(x)\mathbb{P}_{jl}^{(0)}(z)\\
   \\
                             &+& \displaystyle H_{ij}^{(0)}\frac12\left(\sum_{\alpha}\mathcal{P}^{(0)}_{i\alpha}(x)\otimes
       \mathcal{P}_{\alpha k}^{(0)}(x)-\beta_0 \mathcal{P}_{ik}^{(0)}(x)\right)\delta_{jl}\delta(1-z)\\
   \\
                                 &+&\displaystyle H_{ij}^{(0)}\delta_{ik}\delta(1-x) \frac12 \left(\sum_{\beta}\mathbb{P}^{(0)}_{j\beta}(z)\otimes
       \mathbb{P}_{\beta l}^{(0)}(z) -\beta_0 \mathbb{P}_{jl}^{(0)}(z)\right)\,,\\
     \\
     {B}_{ij,kl}^{(2,3)}(x,z) &=&\displaystyle
                                 H_{ij}^{(0)}\widetilde{S}^{(2,3)}\delta_{ik}\delta(1-x) \delta_{jl}\delta(1-z)
                                 - H_{ij}^{(0)} S^{(1,2)}\left[\mathcal{P}_{ik}^{(0)}(x)\delta_{jl}\delta(1-z)
                                 +
                                 \delta_{ik}\delta(1-x)\mathbb{P}_{jl}^{(0)}(z)\right]\,,\\
     \\
     {B}_{ij,kl}^{(2,4)}(x,z) &=&\displaystyle H_{ij}^{(0)}\widetilde{S}^{(2,4)}\delta_{ik}\delta(1-x) \delta_{jl}\delta(1-z) \,.
   \end{array}

.. container::
   :name: a-note-on-the-logarithmic-ordering

   .. rubric:: A note on the logarithmic ordering
      :name: a-note-on-the-logarithmic-ordering

The logarithmic accuracy is driven by the evolution kernel :math:`R`
defined in Eq. :eq:`eq:evkernel`. As it should be clear
from Eq. :eq:`eq:exp1`, the evolution kernel squared
:math:`R^2`, also referred to as Sudakov form factor, admits the
following expansion:

.. math:: R^2 = \sum_{n=0}^\infty a_s^{n}\sum_{k=1}^{2n}\widetilde{S}^{(n,k)}L^k\,.

This expansion naturally allows us to define a logarithmic ordering such
that:

.. math:: R^2 = \sum_{m=0}^{\infty}  R_{\rm N^mLL}^2\,,

with:

.. math:: R_{\rm N^mLL}^2 = \sum_{n=[m/2]}^\infty \widetilde{S}^{(n,2n-m)}\,a_s^{n}L^{2n-m}\,,

where :math:`[m/2]` stands for “integer part of :math:`m/2`.” Now, I
multiply :math:`R_{\rm N^mLL}^2` by :math:`a_s` to a generic power
:math:`p` and do some manipulation defining :math:`n+p=j`. This way, I
get:

.. math::
   :label: eq:LogAcc
   a_s^p R_{\rm N^mLL}^2 = \sum_{j=[(m+2p)/2]}^\infty
   \widetilde{S}^{(j-p,2j-(m+2p))}\,a_s^{j}L^{2j-(m+2p)}\sim R_{\rm N^{m+2p}LL}^2\,,

where the symbol :math:`\sim` means that the terms at two sides have the
same logarithmic accuracy. This finding is relevant in that, in order to
compute a cross section, the Sudakov from factor has to be multiplied by
the hard factor :math:`H` and a pair of matching functions :math:`C`.
Both these components admit an expansion in powers of :math:`a_s`
without enhancing logaritms :math:`L` (see
Eqs. :eq:`eq:MatchCoeffExp`
and :eq:`eq:HardFactExp`). Therefore :math:`R^2`
will eventually be multiplied by some power of :math:`a_s`.

Suppose we compute the Sudakov form factor to leading-logarithmic (LL)
accuracy (:math:`m=0`): Eq. :eq:`eq:LogAcc` tells us
that, in terms of logarithmic accuracy, including
:math:`\mathcal{O}(a_s)` corrections in :math:`H` and :math:`C` implies
introducing next-to-next-to-leading-logarithmic (NNLL) corrections. In
practice, this means that to LL and NLL the functions :math:`H` and
:math:`C` can be taken at leading order, at NNLL we only need to
introduce :math:`\mathcal{O}(a_s)` corrections in :math:`H` and
:math:`C`, while at N\ :math:`^3`\ LL we have to introduce
:math:`\mathcal{O}(a_s^2)` corrections in :math:`H` and :math:`C`.

As a final remark, we notice that, as opposed to the Sudakov form
factor, the evolution of collinear distributions and :math:`\alpha_s`
provide resummation of single logarithms. Therefore, the inclusion of
perturbative corrections to the :math:`\beta`-functions and splitting
functions should be accordingly. Specifically, the perturbative accuracy
at which they are computed should be the same as that of the non-cusp
anomalous dimension :math:`\gamma_F`.

.. container::
   :name: non-perturbative-effect-on-the-asymptotic-cross-section

   .. rubric:: Non-perturbative effect on the asymptotic cross section
      :name: non-perturbative-effect-on-the-asymptotic-cross-section

As discussed in Sect. `2 <#sec:nonpert>`__, TMDs have a non-perturbative
component that can be parameterised as in
Eq. :eq:`eq:separatation`. This implies two main
ingredients: the introduction of a function :math:`b_*(b_T)` that
prevents the perturbative component to enter the non-perturbative regime
and a function :math:`f_{\rm NP}` that parametrises the actual
non-perturbative component. These changes leave almost unchanged the
derivation of the asymptotic behavior of a cross section for
:math:`q_T\rightarrow 0`. The only change concerns the functions
:math:`I_p` defined in Eq. :eq:`eq:IkPert`, that have to
be replaced with:

.. math::
   :label: eq:IkNonPert
     \widetilde{I}_p(x,z,q_T,Q)\equiv \frac12\int_0^\infty db_T\,b_T J_0(b_Tq_T) f_{\rm NP}^{(1)}(x,
     {b}_T,Q^2) f_{\rm NP}^{(2)}(z,
     {b}_T,Q^2)\ln^p\left(\frac{b_*^2(b_T)Q^2}{4e^{-2\gamma_E}}\right)\,,

where :math:`f_{\rm NP}^{(1)}` and :math:`f_{\rm NP}^{(2)}` are the
non-perturbative functions associated to the two TMD distributions
:math:`F` and :math:`D` involved in the cross sections. Here we are
assuming :math:`\mu=\sqrt{\zeta}=Q` and
:math:`\mu_0=\sqrt{\zeta_0}=\mu_b`. In general, the integral in
Eq. :eq:`eq:IkNonPert` cannot be evaluated
analytically but one can use the Ogata quadrature method to compute it
numerically.( [4]_) In the case :math:`k=0`
Eq. :eq:`eq:IkNonPert` has an intersting consequence.
Specifically, contrary to the case in which no non-perturbative
contribution is introduced, :math:`\widetilde{I}_0` is different from
zero also for :math:`q_T>0`. Therefore, also the terms in the expansion
in Eq. :eq:`eq:xsecexp` that are independent of
:math:`L` contribute to the cross section differential in :math:`q_T`.
Clearly, any non-perturbative contribution is expected to give a
substantial contribution only in the vicinity of :math:`q_T=0`.

In this respect, one should check that the introduction of the
non-perturbative components does not affect the large-:math:`q_T`
region. This is amounts to show that:

.. math::
   :label: eq:LogAsympt
   \widetilde{I}_p(x,z,q_T,Q)\mathop{\sim}_{q_T\gtrsim Q} I_p(q_T,Q)\,.

Since for small values of :math:`b_T` we (must) have that
:math:`b_*(b_T)\simeq b_T` and that
:math:`f_{\rm NP}(x,b_T,\zeta)\simeq 1`, this implies that for large
values of :math:`q_T` (that is the conjugate variable of :math:`b_T`),
:math:`\widetilde{I}_p` has to tend to :math:`I_p`. This can be shown
numerically on a case-by-case base. In particular, we can compute the
integral in Eq. :eq:`eq:IkNonPert` using the Ogata
quadrature method and compare the results to those in
Eq. :eq:`eq:PureLogs`. For the comparison we use
Eq. :eq:`eq:bstardefCSS` with
:math:`b_{\rm max} = 2e^{-\gamma_E}` and choose as non-perturbative
functions :math:`f_{\rm NP}` the following :math:`x`-independent form:

.. math:: f_{\rm NP}^{(1)}({b}_T,\zeta) = f_{\rm NP}^{(2)}({b}_T,\zeta) = \exp\left[\left( - g_1 - g_2 \ln\left( \frac{\sqrt{\zeta}}{2 Q_0}\right) \right) \frac{b_T^2}{2}\right]\,,

with :math:`g_1 = 0.02`, :math:`g_1=0.5`, and :math:`Q_0 = 1.6` GeV.

.. container::
   :name: fig:logs

   .. figure:: ../../latex/src/../../latex/src/plots/logs.pdf
   .. :alt: Comparison of the integral in
      Eq. :eq:`eq:IkNonPert` for :math:`p=0,1,2,3`
      (solid lines) to the expressions in
      Eq. :eq:`eq:PureLogs` (dashed lines) as functions
      of :math:`q_T` at :math:`Q = 10` GeV. The black dashed line is not
      present because it corresponds to :math:`\delta(q_T)`.[fig:logs]
      :name: fig:logs
      :width: 45.0%

      Comparison of the integral in
      Eq. :eq:`eq:IkNonPert` for :math:`p=0,1,2,3`
      (solid lines) to the expressions in
      Eq. :eq:`eq:PureLogs` (dashed lines) as functions
      of :math:`q_T` at :math:`Q = 10` GeV. The black dashed line is not
      present because it corresponds to :math:`\delta(q_T)`.[fig:logs]

In Fig. `2 <#fig:logs>`__, I compare the integral in
Eq. :eq:`eq:IkNonPert` for :math:`p=0,1,2,3` (solid
lines) to the expressions in Eq. :eq:`eq:PureLogs`
(dashed lines) at :math:`Q = 10` GeV as functions of :math:`q_T`. Notice
that the black dashed line is not present because it corresponds to
:math:`\delta(q_T)`. It is thus true that :math:`\widetilde{I}_p` tends
to :math:`I_p` as :math:`q_T` increases. However, this is not enough to
satisfy Eq. :eq:`eq:LogAsympt`. As a matter of fact,
if one zooms in around :math:`q_T\lesssim Q`, one finds a substantial
difference between :math:`\widetilde{I}_p` and :math:`I_p`.

.. container::
   :name: fig:logsRatio

   .. figure:: ../../latex/src/../../latex/src/plots/logsRatio.pdf
   .. :alt: Ratio between :math:`\widetilde{I}_1` and :math:`I_1` in the
      region of :math:`q_T` close to :math:`Q`.[fig:logsRatio]
      :name: fig:logsRatio
      :width: 45.0%

      Ratio between :math:`\widetilde{I}_1` and :math:`I_1` in the
      region of :math:`q_T` close to :math:`Q`.[fig:logsRatio]

In Fig `3 <#fig:logsRatio>`__ the zoom is showed for :math:`p=1`. Even
if :math:`q_T` is abundantly in the fixed-order region, the difference
between :math:`\widetilde{I}_1` and :math:`I_1` is still large.
Therefore, it appears that the relation in
Eq. :eq:`eq:LogAsympt` is not plainly fulfilled.
Moreover, the discrepancy tends to become larger as :math:`p` increases.
Therefore, when performing the matching between resummed and fixed-order
cross sections, it is very important to include non-perturbative effects
in the expanded calculation in order to ensure a proper cancellation
with the resummed calculation in the large-:math:`q_T` region.

On the other hand, by design, the inclusion of the non-perturbative
effects modifies significantly also the low-:math:`q_T` region. The
modification is such to prevent the cancellation at small :math:`q_T`
between the expansion of the resummed calculation and the fixed-order
one. A possible solution is to use yet another definition of the
functions :math:`I_p` in the expanded calculation. These functions have
to be such to tend to those in Eq. :eq:`eq:PureLogs`
for small values of :math:`q_T` but converge to the definition in
Eq. :eq:`eq:IkNonPert` more rapidly as :math:`q_T`
increases. To do this, one may exploit the fact that power-suppressed
contributions do not affect the logarithmic expansion of the resummed
calculation to combine the two definitions as follows:

.. math::
   :label: eq:logsTrans
   \hat{I}_p(x,z,q_T,Q) = \left[1-\left(\frac{q_T}{Q}\right)^S\right]I_p(q_T,Q) + \left(\frac{q_T}{Q}\right)^S\widetilde{I}_p(x,z,q_T,Q)\,,

where the exponent :math:`S` can be adjusted to make the transition from
one regime to the other more or less strong. Clearly, this definition
should not be pushed to values of :math:`q_T` much above :math:`Q`. As
an example, Fig. `4 <#fig:logsTrans>`__ shows the shape of
:math:`\hat{I}_1` defined in Eq. :eq:`eq:logsTrans`
for :math:`S = 0.2` along with its components given in
Eqs. :eq:`eq:PureLogs`
and :eq:`eq:IkNonPert`.

.. container::
   :name: fig:logsTrans

   .. figure:: ../../latex/src/../../latex/src/plots/logsTrans.pdf
   .. :alt: Behaviour of :math:`\hat{I}_1` defined in
      Eq. :eq:`eq:logsTrans` for :math:`S = 0.2`
      compared to the single components defined in Comparison of the
      integral in Eqs. :eq:`eq:PureLogs`
      and :eq:`eq:IkNonPert`.[fig:logsTrans]
      :name: fig:logsTrans
      :width: 45.0%

      Behaviour of :math:`\hat{I}_1` defined in
      Eq. :eq:`eq:logsTrans` for :math:`S = 0.2`
      compared to the single components defined in Comparison of the
      integral in Eqs. :eq:`eq:PureLogs`
      and :eq:`eq:IkNonPert`.[fig:logsTrans]

Semi-inclusive deep-inelastic-scattering
========================================

.. container::
   :name: fixed-order-and-asymptotic-limit

   .. rubric:: Fixed order and asymptotic limit
      :name: fixed-order-and-asymptotic-limit

In order to validate the results above, it is opportune to compare the
:math:`\mathcal{O}(a_s)` expressions to those present in the literature.
To this end, we write explicitly the expression for the SIDIS cross
section differential in :math:`q_T` for :math:`q_T>0`, *i.e.* without
the :math:`\delta(q_T)` terms, and with no non-perturbative effect.
Considering that, :math:`S^{(1,1)} = 6C_F`, :math:`S^{(1,2)} = -2C_F`,
and :math:`H_{ij}^{(0)} = e_i^2\delta_{ij}`, this yields:

.. math::
   :label: eq:asymptfromres
   \begin{array}{rcl}
     \displaystyle  \frac{d\sigma}{dx dy dz dq_T^2} &\propto &\displaystyle a_s(Q)\sum_{ij,kl}\left[
                                                      -{B}_{ij,kl}^{(1,1)}(x,z)\frac{1}{q_T^2}
                                                      -{B}_{ij,kl}^{(1,2)}(x,z)\frac{2}{q_T^2}\ln\left(\frac{Q^2}{q_T^2}\right)\right]\mathop{\otimes}_x f_k(x,Q) \mathop{\otimes}_z d_l(z,Q)\\
     \\
                                           &=&\displaystyle
                                               \frac{a_s(Q)}{q_T^2}\sum_{ij,kl}e_i^2\delta_{ij}\bigg[4C_F\left(\ln\left(\frac{Q^2}{q_T^2}\right)-\frac{3}{2}\right)\delta_{ik}\delta_{jl}\delta(1-x)\delta(1-z)\\
     \\
                                           &+&\displaystyle
                                               \mathcal{P}_{ik}^{(0)}(x)\delta_{jl}\delta(1-z)
                                               +
                                               \delta_{ik}\delta(1-x)\mathbb{P}_{jl}^{(0)}(z)\bigg]\mathop{\otimes}_x
                                               f_k(x,Q) \mathop{\otimes}_z
                                               d_l(z,Q)\\
     \\
                                           &=&\displaystyle \frac{a_s(Q)}{q_T^2}\sum_{i}e_i^2\bigg[4C_F\left(\ln\left(\frac{Q^2}{q_T^2}\right)-\frac{3}{2}\right) f_i(x,Q)d_i(z,Q)\\
     \\
                                           &+&\displaystyle \left(\sum_{k}\mathcal{P}_{ik}^{(0)}(x) \mathop{\otimes}_x
                                               f_k(x,Q)\right)d_i(z,Q) + f_i(x,Q)\left(\sum_l\mathbb{P}_{il}^{(0)}(z) \mathop{\otimes}_zd_l(z,Q)\right)\bigg] +\mathcal{O}(a_s^2)\,.
   \end{array}

\ This result, up to pre-factors that will be made explicit below,
nicely agrees with that of, *e.g.*, Refs. (Meng, Olness, and Soper 1996;
J. Collins et al. 2016; Nadolsky, Stump, and Yuan 2000).

In order to check that the matching is actually removing the double
counting terms, it is instructive to derive
Eq. :eq:`eq:asymptfromres` extracting the
asymptote from the fixed-order computation at :math:`\mathcal{O}(a_s)`.
We take the expressions for the coefficient functions from
Eqs. (106)-(109) of Appendix B of Ref. (Nadolsky, Stump, and Yuan 2000)
or from Eqs. (4.6)-(4.20) of Ref. (Bacchetta et al. 2008). Referring to
the second reference, some simplifications apply. In particular, we
consider cross sections with unpolarised projectiles
(:math:`\lambda_e=0`) on unpolarised targets (:math:`S_\perp^\mu = 0`)
and integrated over the azimuthal angles :math:`\phi_H` and
:math:`\phi_S`. By doing so and after a simple manipulation, the cross
section simplifies greatly and can be written in terms of structure
functions as:

.. math::
   :label: eq:xsexinsf
     \frac{d\sigma}{dx dy dz dq_T^2} = \frac{2\pi\alpha^2}{xyQ^2}\left[Y_+ F_{UU,T}+2(1-y)F_{UU,L}\right]=\frac{2\pi\alpha^2}{xyQ^2}Y_+\left[F_{UU,2}-\frac{y^2}{Y_+}F_{UU,L}\right]\,,

with:

.. math:: Y_+\equiv 1+(1-y)^2\,,

and where we have defined the structure function:

.. math:: F_{UU,2}\equiv F_{UU,T} + F_{UU,L}\,.

Notice that, as compared to Ref. (Bacchetta et al. 2008), we have
factored out from the structure functions a factor
:math:`1/(\pi z^2)`\ ( [5]_) so that they factorize as:

.. math::
   :label: eq:FOxsec
   \begin{array}{rcl}
     F_{UU,S}&=&\displaystyle 
     a_s\frac{x}{Q^2}\sum_{i}e_i^2\int_x^1\frac{d\bar{x}}{\bar{x}}\int_z^1\frac{d\bar{z}}{\bar{z}}\delta\left(\frac{q_T^2}{Q^2}-\frac{(1-\bar{x})(1-\bar{z})}{\bar{x}\bar{z}}\right)\bigg[\hat{B}_{qq}^{S,
       \rm FO}(\bar{x},\bar{z},q_T)f_i\left(\frac{x}{\bar{x}}\right)
                 d_i\left(\frac{z}{\bar{z}}\right)\\
   \\
   &+&\displaystyle \hat{B}_{qg}^{S,
       \rm FO}(\bar{x},\bar{z},q_T)f_g\left(\frac{x}{\bar{x}}\right)
                 d_i\left(\frac{z}{\bar{z}}\right)+ \hat{B}_{gq}^{S,
       \rm FO}(\bar{x},\bar{z},q_T)f_i\left(\frac{x}{\bar{x}}\right)
                 d_g\left(\frac{z}{\bar{z}}\right)\bigg]
   +\mathcal{O}(a_s^2)\,.
   \end{array}

with :math:`S=2,L` and where the sum over :math:`i` runs over the active
quark and antiquark flavours. The explicit expressions for the
coefficient functions are:

.. math::
   :label: eq:Bacchettaetal
   \begin{array}{l}
   \displaystyle \hat{B}_{qq}^{2,\rm FO}(x,z,q_T) = 2C_F\left[(1-x)(1-z)+4xz+\frac{1+x^2z^2}{xz}\frac{Q^2}{q_T^2}\right]\,,\\
   \\
   \displaystyle \hat{B}_{qq}^{L,\rm FO}(x,z,q_T) = 8C_Fxz\,,\\
   \\
   \displaystyle \hat{B}_{qg}^{2,\rm FO}(x,z,q_T) = 2T_R\left[[x^2+(1-x)^2][z^2+(1-z)^2]\frac{1-x}{xz^2}\frac{Q^2}{q_T^2} +8x(1-x)\right]\,,\\
   \\
   \displaystyle \hat{B}_{qg}^{L,\rm FO}(x,z,q_T) = 16T_Rx(1-x)\,.\\
   \\
   \displaystyle \hat{B}_{gq}^{2,\rm FO}(x,z,q_T) = 2C_F\left[(1-x)z+4x(1-z)+\frac{1+x^2(1-z)^2}{xz}\frac{1-z}{z}\frac{Q^2}{q_T^2}\right]\,,\\
   \\
   \displaystyle \hat{B}_{gq}^{L,\rm FO}(x,z,q_T) = 8C_Fx(1-z)\,,
   \end{array}

These expressions are enough to compute the SIDIS cross section at
:math:`\mathcal{O}(a_s)` in the region :math:`q_T \lesssim Q`. In order
to match Eq. :eq:`eq:asymptfromres`, one has to
take the limit :math:`q_T/Q\rightarrow 0` and retain in the coefficient
functions only the terms enhanced as :math:`\ln(Q^2/q_T^2)`. This
automatically means that :math:`F_{UU,L}` does not contribute in this
limit because it contains no logarithmic enhancements such that:

.. math:: F_{UU,L}\mathop{\ll}_{q_T/Q\rightarrow 0} F_{UU,2}\,.

Another crucial observation is that the :math:`\delta`-function in
Eq. :eq:`eq:FOxsec` can be expanded as follows( [6]_):

.. math::
   :label: eq:deltaexpansion
     \delta\left(\frac{q_T^2}{Q^2}-\frac{(1-x)(1-z)}{xz}\right)
     \mathop{\longrightarrow}_{q_T^2/Q^2\rightarrow0}\ln\left(\frac{Q^2}{q_T^2}\right)\delta(1-x) \delta(1-z) +
     \frac{x\delta(1-z)}{(1-x)_+}+ \frac{z\delta(1-x)}{(1-z)_+}\,,

so that:

.. math::
   :label: eq:FOxsecasy
   \begin{array}{rcl}
     F_{UU,2}&\displaystyle \mathop{\longrightarrow}_{q_T/Q\rightarrow 0}&\displaystyle 
     a_s\frac{x}{q_T^2}\sum_{i}e_i^2\int_x^1\frac{d\bar{x}}{\bar{x}}\int_z^1\frac{d\bar{z}}{\bar{z}}\bigg[\hat{B}_{qq}^{2,
       \rm asy}(\bar{x},\bar{z},q_T)f_i\left(\frac{x}{\bar{x}}\right)
                 d_i\left(\frac{z}{\bar{z}}\right)\\
   \\
   &+&\displaystyle \hat{B}_{qg}^{2,
       \rm asy}(\bar{x},\bar{z},q_T)f_g\left(\frac{x}{\bar{x}}\right)
                 d_i\left(\frac{z}{\bar{z}}\right)+ \hat{B}_{gq}^{2,
       \rm asy}(\bar{x},\bar{z},q_T)f_i\left(\frac{x}{\bar{x}}\right)
                 d_g\left(\frac{z}{\bar{z}}\right)\bigg]
   +\mathcal{O}(a_s^2)\,.
   \end{array}

with:

.. math::
   \begin{array}{rcl}
     \hat{B}_{qq}^{2,\rm asy}(x,z,q_T) &=& \displaystyle  2C_F\left[2
                                           \ln\left(\frac{Q^2}{q_T^2}\right)+\frac{1+x^2}{(1-x)_+}\delta(1-z) +\delta(1-x)\frac{1+z^2}{(1-z)_+} \right]\\
     \\
                                       &=&\displaystyle 2C_F\left[2\ln\left(\frac{Q^2}{q_T^2}\right) -3\right]\delta(1-x) \delta(1-z)+\mathcal{P}_{qq}^{(0)}(x)\delta(1-z)
                                           +\delta(1-x) \mathbb{P}_{qq}^{(0)}(z)\,,\\
     \\
     \hat{B}_{qg}^{2,\rm asy}(x,z,q_T) &=& \displaystyle
                                           2T_R\left[x^2+(1-x)^2\right]\delta(1-z)=\mathcal{P}_{qg}^{(0)}(x)
                                           \delta(1-z)\,,\\
   \\
     \hat{B}_{gq}^{2,\rm asy}(x,z,q_T) &=& \displaystyle  \delta(1-x)2C_F\left[\frac{1+(1-z)^2}{z}\right]=\delta(1-x)\mathbb{P}_{qg}^{(0)}(z)\,.
   \end{array}

It is thus easy to see that we can rewrite
Eq. :eq:`eq:FOxsecasy` as:

.. math::
   :label: eq:FOxsecasy2
   \begin{array}{rcl}
     F_{UU,2}&\displaystyle \mathop{\longrightarrow}_{q_T/Q\rightarrow 0}&\displaystyle 
     a_s\frac{x}{q_T^2}\sum_{i}e_i^2 \Bigg[4C_F\left(\ln\left(\frac{Q^2}{q_T^2}\right) -\frac{3}{2}\right) f_i\left(x\right)
                 d_i\left(z\right)\\
   \\
   &+&\displaystyle
       \left(\sum_{k=q,g}\mathcal{P}^{(0)}_{qk}(x)\otimes f_k\left(x\right) \right)
                 d_i\left(z\right)+ f_i\left(x\right)
                 \left(\sum_{k=q,g}\mathbb{P}^{(0)}_{qk}(z) \otimes d_k\left(z\right)\right)\Bigg]
   +\mathcal{O}(a_s^2)\,.
   \end{array}

Therefore, up to factors, Eq. :eq:`eq:FOxsecasy2`
agrees with Eq. :eq:`eq:asymptfromres`. This
confirms that the expansion of the resummed calculation, as well as the
asymptotic limit of the fixed order, removes the double-counting terms
when doing the matching.

In order to provide a version of Eq. :eq:`eq:FOxsec` that
can be readily implemented, we need to perform one of the integrals
making use of the :math:`\delta`-function. We integrate over
:math:`\bar{x}` so that we write:

.. math::
   \delta\left(\frac{q_T^2}{Q^2}-\frac{(1-\bar{x})(1-\bar{z})}{\bar{x}\bar{z}}\right)
   = \frac{\bar{z}\bar{x}_0^2}{1-\bar{z}}\delta(\bar{x} - \bar{x}_0)\,,

with:

.. math::
   :label: eq:x0bar
     \bar{x}_0 = \frac{1-\bar{z}}{1-\bar{z}\left(1-\frac{q_T^2}{Q^2}\right)}\,.

This allows us to write:

.. math::
   :label: eq:FOxsecFin
   \begin{array}{rcl}
     F_{UU,S}&=&\displaystyle 
     a_s\frac{x }{Q^2}\sum_{i}e_i^2\int_z^{z_{\rm max}} \frac{d\bar{z}}{1-\bar{z}}\bar{x}_0\bigg[\hat{B}_{qq}^{S,
       \rm FO}(\bar{x}_0,\bar{z},q_T)f_i\left(\frac{x}{\bar{x}_0}\right)
                 d_i\left(\frac{z}{\bar{z}}\right)\\
   \\
   &+&\displaystyle \hat{B}_{qg}^{S,
       \rm FO}(\bar{x}_0,\bar{z},q_T)f_g\left(\frac{x}{\bar{x}_0}\right)
                 d_i\left(\frac{z}{\bar{z}}\right)+ \hat{B}_{gq}^{S,
       \rm FO}(\bar{x}_0,\bar{z},q_T)f_i\left(\frac{x}{\bar{x}_0}\right)
                 d_g\left(\frac{z}{\bar{z}}\right)\bigg]
   +\mathcal{O}(a_s^2)\,,
   \end{array}

with:

.. math:: z_{\rm max} = \frac{1-x}{1-x\left(1-\frac{q_T^2}{Q^2}\right)}\,.

Now we can rewrite the cross section above in such a way that it matches
that at :math:`\mathcal{O}(a_s)` of Ref. (Daleo, Florian, and Sassot
2005). That would allow us to confidently use the
:math:`\mathcal{O}(a_s^2)` calculation presented in that reference for
the matching to the resummed calculation. This is made tricky by the
different notation used in Ref. (Daleo, Florian, and Sassot 2005) and
from the fact that in that paper the cross section is differential in a
different set of variables. Specifically, we would like it to be
differential in :math:`x`, :math:`y`, :math:`z`, and :math:`q_T^2` while
in Ref. (Daleo, Florian, and Sassot 2005) it is differential in
:math:`x`, :math:`Q^2`, :math:`\eta`, and :math:`p_T^2`, where the last
two are the rapidity and the transverse momentum of the outgoing hadron,
respectively. Eq. (13) of Ref. (Daleo, Florian, and Sassot 2005), can be
translated into our notation by noticing that:

.. math::
   :label: eq:Daleoetal
   \frac{d\sigma}{dxdQ^2d p_T^2 d\eta} = \frac{x}{z Q^2}\sum_{i,j}\int_z^{z_{\rm
       max}} \frac{d\bar{z}}{1-\bar{z}} f_i\left(\frac{x}{\bar{x}_0}\right)d_j\left(\frac{z}{\bar{z}}\right) \frac{d\sigma_{ij}^{(1)}}{dxdQ^2d p_T^2 d\eta}+\mathcal{O}(a_s^2)\,,

where we have exploited the :math:`\delta`-functions in Eqs. (18)-(20)
to get rid of the integral over :math:`z`.( [7]_)

The :math:`\mathcal{O}(a_s)` partonic cross sections in Eqs. (18)-(20)
of Ref. (Daleo, Florian, and Sassot 2005), setting
:math:`\varepsilon=0`, can be written as:

.. math::
   \frac{d\sigma_{ij}^{(1)}}{dxdQ^2d p_T^2 d\eta}=\frac{2\pi \alpha^2 a_s
     e_q^2\bar{x}_0}{x Q^4}Y_+ \left[ \underbrace{\left( F_{UU,M}^{ij}(\bar{x}_0,\bar{z})+\frac{3}{2}F_{UU,L}^{ij}(\bar{x}_0,\bar{z})\right)}_{F_{UU,2}^{ij}} - \frac{y^2}{Y_+}F_{UU,L}^{ij}(\bar{x}_0,\bar{z})\right]\,.

One can verify that :math:`F_{UU,2}^{qq}(\bar{x}_0,\bar{z})`,
:math:`F_{UU,L}^{qq}(\bar{x}_0,\bar{z})`,
:math:`F_{UU,M}^{qg}(\bar{x}_0,\bar{z})`,
:math:`F_{UU,L}^{qg}(\bar{x}_0,\bar{z})`,
:math:`F_{UU,M}^{gq}(\bar{x}_0,\bar{z})`, and
:math:`F_{UU,L}^{gq}(\bar{x}_0,\bar{z})` correspond exactly to
:math:`\hat{B}^{2, \rm FO}_{qq}(\bar{x}_0,\bar{z},q_T)`,
:math:`\hat{B}^{L, \rm FO}_{qq}(\bar{x}_0,\bar{z},q_T)`,
:math:`\hat{B}^{M, \rm FO}_{qg}(\bar{x}_0,\bar{z},q_T)`,
:math:`\hat{B}^{L, \rm FO}_{qg}(\bar{x}_0,\bar{z},q_T)`,
:math:`\hat{B}^{M, \rm FO}_{gq}(\bar{x}_0,\bar{z},q_T)`, and
:math:`\hat{B}^{L, \rm FO}_{gq}(\bar{x}_0,\bar{z},q_T)` in
Eq. :eq:`eq:Bacchettaetal` of Ref. (Bacchetta et
al. 2008). It is crucial to notice that the correspondence holds only if
:math:`\bar{x}_0` defined in Eq. :eq:`eq:x0bar` as a
function of :math:`\bar{z}` is used. As an example, taking into account
the different factors (a factor 2 for :math:`F_{UU,M}` and a factor 4
for :math:`F_{UU,L}`) and using our notation for the integration
variables (:math:`y\rightarrow \bar{z}`,
:math:`\rho(z=0)\rightarrow \bar{x}_0`), we read off from Eq. (18) of
Ref. (Daleo, Florian, and Sassot 2005):

.. math::
   \begin{array}{rcl}
   F_{UU,M}^{qq}(\bar{x}_0,\bar{z}) &=& \displaystyle 2C_F\left[\frac{(\bar{x}_0+\bar{z})^2+2(1-\bar{x}_0-\bar{z})}{(1-\bar{x}_0)(1-\bar{z})}\right]\,,\\
   \\
   F_{UU,L}^{qq}(\bar{x}_0,\bar{z}) &=& \displaystyle 8C_F\bar{x_0}\bar{z}\,,
   \end{array}

so that:

.. math::
   \begin{array}{rcl}
   F_{UU,2}^{qq}(\bar{x}_0,\bar{z}) &=& \displaystyle F_{UU,M}^{qq}(\bar{x}_0,\bar{z})+\frac{3}{2}F_{UU,L}^{qq}(\bar{x}_0,\bar{z})\\
   \\
    &=& \displaystyle 2
        C_F\left[\frac{(\bar{x}_0+\bar{z})^2+2(1-\bar{x}_0-\bar{z})}{(1-\bar{x}_0)(1-\bar{z})}
        + 3\bar{x_0}\bar{z}\right]\\
   \\
   &=& \displaystyle 2
        C_F\left[(1-\bar{x}_0)(1-\bar{z})+4\bar{x}_0\bar{z} +\frac{1+\bar{x}_0^2\bar{z}^2}{\bar{x}_0\bar{z}}\left(\frac{\bar{x}_0\bar{z}}{(1-\bar{x}_0)(1-\bar{z})}\right)\right]\,.
   \end{array}

Using Eq. :eq:`eq:x0bar`, it is easy to see that the
factor in round brackets in the last line of the equation above is equal
to :math:`Q^2/q_T^2`. Therefore, it reduces exactly to the first
relation in Eq. :eq:`eq:Bacchettaetal` of
Ref. (Bacchetta et al. 2008). The same holds also for the two remaining
partonic channels.

Putting all pieces together, Eq. :eq:`eq:Daleoetal`
can be recast as:

.. math::
   \frac{d\sigma}{dxdQ^2d p_T^2 d\eta} = \frac{2\pi \alpha^2
     }{z x Q^4}Y_+ \left[F_{UU,2} - \frac{y^2}{Y_+}F_{UU,L}\right]\,,

with:

.. math::
   \begin{array}{rcl}
     F_{UU,S}&=&\displaystyle 
     a_s\frac{x }{Q^2}\sum_{i}e_i^2\int_z^{z_{\rm max}}
                 \frac{d\bar{z}}{1-\bar{z}}\bar{x}_0\bigg[F^{qq}_{UU,L}(\bar{x}_0,\bar{z})f_i\left(\frac{x}{\bar{x}_0}\right)
                 d_i\left(\frac{z}{\bar{z}}\right)\\
   \\
   &+&\displaystyle F^{qg}_{UU, L}(\bar{x}_0,\bar{z})f_g\left(\frac{x}{\bar{x}_0}\right)
                 d_i\left(\frac{z}{\bar{z}}\right)+ F^{gq}_{UU, L}(\bar{x}_0,\bar{z})f_i\left(\frac{x}{\bar{x}_0}\right)
                 d_g\left(\frac{z}{\bar{z}}\right)\bigg]
   +\mathcal{O}(a_s^2)\,.
   \end{array}

Therefore, the structure of the observables is exactly the same. What is
left to work out is the Jacobian to express the cross section
differential in the same variables as in
Eq. :eq:`eq:xsexinsf`. What we need is to know is how
the variables :math:`Q^2`, :math:`p_T`, and :math:`\eta` are related to
:math:`y`, :math:`z`, and :math:`q_T`. The relevant relations are:

.. math::
   :label: eq:variablechange
   \left\{\begin{array}{rcl}
   Q^2 &=& xyS_H\\
   \\
   p_T^2 &=& z^2 q_T^2\\
   \\
   \eta &=& \displaystyle \frac{1}{2}\ln\left(\frac{y(1-x)S_H}{q_T^2}\right)
   \end{array}\right.\quad\Longrightarrow\quad dQ^2dp_T^2d\eta = \frac{zQ^2}{y} dydzdq_T^2\,.

where :math:`S_H` is the squared center-of-mass energy of the collision,
so that:

.. math::
   :label: eq:differentialversion2
     \frac{d\sigma}{dx dy dz dq_T^2} = \frac{zQ^2}{y}\frac{d\sigma}{dxdQ^2d p_T^2
       d\eta} = \frac{2\pi\alpha^2}{xyQ^2}Y_+\left[F_{UU,2}-\frac{y^2}{Y_+}F_{UU,L}\right]\,,

exactly like in Eq. :eq:`eq:xsexinsf`. This confirms
that the computation of Ref. (Daleo, Florian, and Sassot 2005) can be
matched to the resummed calculation provided that the correct Jacobian
is taken into account.

In Fig. `5 <#fig:APFELvsTIMBA>`__, a comparison between the
:math:`\mathcal{O}(\alpha_s)` computation of Ref. (Daleo, Florian, and
Sassot 2005) for the differential cross section in
Eq. :eq:`eq:Daleoetal`, implemented in the computer
code ``TIMBA``, and the implementation in the ``APFEL++``
framework (Bertone 2018) of the expressions in
Eq. :eq:`eq:FOxsecFin` is presented. The cross section
as a function of the hadron transverse momentum :math:`p_T` obtained
with ``APFEL++`` is displayed as ratios to ``TIMBA`` for representative
values of the kinematic variables: the agreement between the two codes
is excellent. Since we will be using ``TIMBA`` for the NLO (*i.e.*
:math:`\mathcal{O}(a_s^2)`) cross section calculation to be matched to
the NNLL resummed computation, this provides a solid ground to start
from.

.. container::
   :name: fig:APFELvsTIMBA

   .. figure:: ../../latex/src/../../latex/src/plots/APFELvsTIMBA.pdf
   .. :alt: Comparison at :math:`\mathcal{O}(\alpha_s)` between the
      `` TIMBA`` code, implementation of the results of Ref. (Daleo,
      Florian, and Sassot 2005), and the implementation of the
      expressions in Eq. :eq:`eq:FOxsecFin` in the
      ``APFEL++`` code (Bertone 2018).[fig:APFELvsTIMBA]
      :name: fig:APFELvsTIMBA
      :width: 50.0%

      Comparison at :math:`\mathcal{O}(\alpha_s)` between the `` TIMBA``
      code, implementation of the results of Ref. (Daleo, Florian, and
      Sassot 2005), and the implementation of the expressions in
      Eq. :eq:`eq:FOxsecFin` in the ``APFEL++``
      code (Bertone 2018).[fig:APFELvsTIMBA]

In order to be able to compare the theoretic predictions to the data, it
is useful to work out the for of the cross section differential in the
variables :math:`x`, :math:`Q^2`, :math:`z`, and :math:`p_T^2`. This is
easily done starting from
Eq. :eq:`eq:differentialversion2` and
considering that:

.. math::
   \eta = \ln\left(\frac{zW}{p_T}\right)\,,\quad\mbox{with}\quad W
   =\sqrt{\frac{Q^2(1-x)}{x}}\,,\quad\Rightarrow\quad d\eta = \frac{dz}{z}\,.

It follows that:

.. math::
   :label: eq:differentialversion3
    \frac{d\sigma}{dxdQ^2dz dp_T^2} = \frac{2\pi\alpha^2}{xQ^4}Y_+\left[F_{UU,2}-\frac{y^2}{Y_+}F_{UU,L}\right]\,,

with :math:`y=Q^2/xS_H` (see
Eq. :eq:`eq:variablechange`).

.. container::
   :name: integrating-over-q_t

   .. rubric:: Integrating over :math:`q_T`
      :name: integrating-over-q_t

It is interesting to carry out the integration over :math:`q_T` of
Eq. :eq:`eq:FOxsec` analytically. The result can then be
compared to those presented in Appendix C of Ref. (Florian, Stratmann,
and Vogelsang 1998). Exploiting the :math:`\delta`-function yields:

.. math::
   \begin{array}{rcl}
     \displaystyle \int_0^{\infty} dq_T^2\,F_{UU,S}&=&\displaystyle 
                                                       a_sx\sum_{i}e_i^2\int_x^1\frac{d\bar{x}}{\bar{x}}\int_z^1\frac{d\bar{z}}{\bar{z}}\bigg[\hat{B}_{qq}^{S,
                                                       \rm FO}\left(\bar{x},\bar{z},\frac{(1-\bar{x})(1-\bar{z})}{\bar{x}\bar{z}}\right)f_i\left(\frac{x}{\bar{x}}\right)
                                                       d_i\left(\frac{z}{\bar{z}}\right)\\
     \\
                                                   &+&\displaystyle \hat{B}_{qg}^{S,
                                                       \rm FO}\left(\bar{x},\bar{z},\frac{(1-\bar{x})(1-\bar{z})}{\bar{x}\bar{z}}\right)f_g\left(\frac{x}{\bar{x}}\right)
                                                       d_i\left(\frac{z}{\bar{z}}\right)\\
     \\
                                                   &+& \displaystyle \hat{B}_{gq}^{S,
                                                       \rm FO}\left(\bar{x},\bar{z},\frac{(1-\bar{x})(1-\bar{z})}{\bar{x}\bar{z}}\right)f_i\left(\frac{x}{\bar{x}}\right)
                                                       d_g\left(\frac{z}{\bar{z}}\right)\bigg]
                                                       +\mathcal{O}(a_s^2)\,.
   \end{array}

\ replacing :math:`q_T^2/Q^2` in
Eq. :eq:`eq:Bacchettaetal` has no effect on the
coefficient functions for :math:`F_{UU,L}`. The reason is that
:math:`F_{UU,L}` is independent of :math:`q_T` and therefore integrating
over :math:`q_T` has essentially the effect of removing the
:math:`\delta`-function. On the other hand, integrating :math:`F_{UU,2}`
gives:

.. math::
   \begin{array}{l}
   \displaystyle \hat{B}_{qq}^{2,\rm FO}\left(x,z,\frac{(1-{x})(1-{z})}{{x}{z}}\right) = 2C_F\left[(1-x)(1-z)+4xz+\frac{1+x^2z^2}{(1-x)(1-z)}\right]\,,\\
   \\
   \displaystyle \hat{B}_{qg}^{2,\rm FO}\left
     (x,z,\frac{(1-{x})(1-{z})}{{x}{z}}\right) = 2T_R\left[[x^2+(1-x)^2]\left(\frac{1}{1-z}+\frac1{z}-2\right)+8x(1-x)\right]\,,\\
   \\
   \displaystyle \hat{B}_{gq}^{2,\rm FO}\left (x,z,\frac{(1-{x})(1-{z})}{{x}{z}}\right) =
     2C_F\left[(1-x)z+4x(1-z)+\frac{1+x^2(1-z)^2}{z(1-x)}\right]\,,
   \end{array}

Instating the plus prescription in all terms that behave as
:math:`(1-x)^{-1}` and :math:`(1-z)^{-1}`\ ( [8]_), one can see that up
to “local” terms (*i.e.* those terms proportional to the
:math:`\delta`-function) the expressions of Appendix C of Ref. (Florian,
Stratmann, and Vogelsang 1998) are reproduced( [9]_). The local terms
correspond to virtual corrections that thus give contributions at
:math:`q_T=0`. As such, they should coincide with the terms proportional
to :math:`\delta(q_T)` generated in the expansion of the resummed
calculation. Specifically, they should be equal to the coefficients
:math:`{B}_{ij,kl}^{(n,0)}` in the expansion in
Eq. :eq:`eq:xsecexp`. This is evidently the case for
:math:`{B}_{ij,kl}^{(0,0)}` that gives the leading-order contribution to
the integrated cross section. Unfortunately, this is not true at
next-to-leading order because the expressions for the local terms in
Appendix C of Ref. (Florian, Stratmann, and Vogelsang 1998)
(Eqs. (C.2)-(C.4)) do not match :math:`{B}_{ij,kl}^{(1,0)}`. However,
when checking the expressions explicitly an interesting pattern emerges:
if all the occurrences of :math:`\ln(1-x)` (or :math:`\ln(1-z)`) in
Eqs. (C.2)-(C.4) of Ref. (Florian, Stratmann, and Vogelsang 1998) are
replaced with :math:`\ln(x)` (or :math:`\ln(z)`) the two sets of
expressions coincide. In fact, the expression to
:math:`\mathcal{O}(\alpha_s)` that we derived above are free of large
threshold logarithms associated to soft-gluon emission.

This should not be surprising because the CS equation

.. container::
   :name: integrating-over-z

   .. rubric:: Integrating over :math:`z`
      :name: integrating-over-z

In order to understand which of the two is the correct set of
expressions, whether those computed in Ref. (Furmanski and Petronzio
1982) or those derived here, we may try to integrate over :math:`z`
summing over all possible hadrons produced in the final state. By doing
so, we should obtain the expressions of the hard cross sections for
inclusive deep-inelastic scattering. This last set of expressions has
also been computed in Ref. (Furmanski and Petronzio 1982) and is very
likely to be correct. The presence of terms proportional to
:math:`\ln(1-x)` in the inclusive DIS hard cross sections suggests that
the expressions reported in Ref. (Florian, Stratmann, and Vogelsang
1998) and computed in Ref. (Furmanski and Petronzio 1982) are the
correct ones. To integrate the expressions in Ref. (Florian, Stratmann,
and Vogelsang 1998) we need to use the fact that the SIDIS cross section
integrate in :math:`q_T` for the production of a hadron :math:`h` has
the following structure:

.. math::
   :label: eq:intxsec
   \frac{d\sigma^h}{dx dy dz} = \sum_{ij} \int_x^1\frac{d\xi}{\xi}
   \int_z^1\frac{d\zeta}{\zeta}
   C_{ij}\left(\frac{x}{\xi},\frac{z}{\zeta}\right) d_{i/h}(\zeta) f_j(\xi)\,,

where :math:`d_{i/h}` is the fragmentation function for the parton
:math:`i` fragmenting into the hadron species :math:`h`. What we need to
compute is:

.. math:: \frac{d\sigma}{dx dy} = \sum_{h}\int_0^1 dz\,z\,\frac{d\sigma^h}{dx dy dz}\,,

where the sum over :math:`h` extends over all possible hadronic species.
In order to compute the integral above we need to use the following
property of the Mellin convolution:

.. math::
   C(z)\otimes d(z) = \int_z^1\frac{d\zeta}{\zeta}
   C\left(\frac{z}{\zeta}\right)d(\zeta) = \int_0^1d\zeta \int_0^1d\eta\,C(\zeta)f(\eta)\delta(z-\zeta\eta)\,.

By using this property we have:

.. math::
   \begin{array}{rcl}
     \displaystyle \frac{d\sigma}{dx dy} &=& \displaystyle \sum_{ij} \int_x^1\frac{d\xi}{\xi}
                                             \int_0^1 d\zeta \,
                                             C_{ij}\left(\frac{x}{\xi},\zeta\right)f_j(\xi) \sum_h\int_0^1
                                             d\eta\,d_{i/h}(\eta)\,\int_0^1dz\,z\delta(z-\zeta\eta)\\
     \\
                                         &=& \displaystyle \sum_{j} \int_x^1\frac{d\xi}{\xi}
                                             \left[\int_0^1 d\zeta \,\zeta
                                             \sum_i  C_{ij}\left(\frac{x}{\xi},\zeta\right)\right]f_j(\xi) = \sum_{j} \int_x^1\frac{d\xi}{\xi}\,
                                             \widetilde{C}_{j}\left(\frac{x}{\xi}\right)f_j(\xi)\,.
   \end{array}

where we have defined the inclusive DIS coefficient functions as:

.. math::
   :label: eq:incDIS
   \widetilde{C}_{j}\left(\xi\right) \equiv \int_0^1 d\zeta \,\zeta
                                             \sum_i  C_{ij}\left(\xi,\zeta\right)\,,

and we have used the momentum sum rule:

.. math:: \sum_h\int_0^1 d\eta\,\eta d_{i/h}(\eta) = 1\,,

that essentially encodes the fact that the total probability for a
parton to fragment into any hadron in one. The correctness of
Eq. :eq:`eq:incDIS` can be verified explicitly by
integrating over :math:`z` the expressions reported in Appendix C of
Ref. (Florian, Stratmann, and Vogelsang 1998) for both structure
functions :math:`F_2` and :math:`F_L`. The leading order is trivial:

.. math::
   \widetilde{C}_{2,q}^{(0)}\left(\xi\right) \equiv \int_0^1 d\zeta \,\zeta
                                            \delta(1-\xi)\delta(1-\zeta) = \delta(1-\xi)\,,

as well know, while the gluon as well as the longitudinal coefficient
functions are zero. The next-to-leading order coefficient functions are
computed as:

.. math::
   \begin{array}{l}
   \displaystyle \widetilde{C}_{q}^{(1)}\left(\xi\right) \equiv \int_0^1
     d\zeta\,\zeta
     \left[C_{qq}^{(1)}\left(\xi,\zeta\right)+C_{gq}^{(1)}\left(\xi,\zeta\right)\right]\,,\\
   \\
   \displaystyle \widetilde{C}_{g}^{(1)}\left(\xi\right) \equiv \int_0^1
     d\zeta\,\zeta C_{qg}^{(1)}\left(\xi,\zeta\right)\,.
   \end{array}

For :math:`F_L` they read:

.. math::
   \begin{array}{l}
   \displaystyle \widetilde{C}_{L,q}^{(1)}\left(\xi\right) = 8C_F\xi\,,\\
   \\
   \displaystyle \widetilde{C}_{L,g}^{(1)}\left(\xi\right) = 8T_R\xi(1-\xi)\,,
   \end{array}

while for :math:`F_2`:

.. math::
   \begin{array}{rcl}
   \displaystyle \widetilde{C}_{2,q}^{(1)}\left(\xi\right) &=&
                                                               \displaystyle 2C_F\bigg[2\left(\frac{\ln(1-\xi)}{1-\xi}\right)_+-\frac32\left(\frac{1}{1-\xi}\right)_+-(1+\xi)\ln(1-\xi)\\
   \\
   &-&\displaystyle \frac{1+\xi^2}{1-\xi}\ln\xi+3+2\xi-\left(\frac92+2\zeta_2\right)\delta(1-\xi)\bigg]\,,\\
   \\
   \displaystyle \widetilde{C}_{2,g}^{(1)}\left(\xi\right)
   &=&\displaystyle 2T_R\left[(\xi^2+(1-\xi)^2)\ln\left(\frac{1-\xi}{\xi}\right)+6\xi(1-\xi)-1\right]\,.
   \end{array}

These expressions coincide with those reported, *e.g.*, in Chapter 4 of
Ref. (Ellis, Stirling, and Webber 2011) with the only exception of
:math:`\widetilde{C}_{2,g}^{(1)}` where the factor 6 in front of the
term :math:`\xi(1-\xi)`, according to Ref. (Ellis, Stirling, and Webber
2011), should be a 8. This should be investigated further.

Drell-Yan production
====================

Appendices
==========

.. container::
   :name: sec:deltaexpansion

   .. rubric:: Expansion of the kinematic :math:`\delta`-function
      :name: sec:deltaexpansion

In order to prove the equality in
Eq. :eq:`eq:deltaexpansion`, we consider the
following integral:

.. math::
   :label: eq:testint
   I(\varepsilon) =\int_0^1 dx \int_0^1 dy\,\delta(xy-\varepsilon) f(x,y)\,,

where :math:`f(x,y)` is a test function “well-behaved” over the
integration region. Now I split the integral as follows:

.. math:: I(\varepsilon) =\left(\int_0^1 dx \int_x^1 dy+\int_0^1 dy \int_y^1 dx\right)\,\delta(xy-\varepsilon) f(x,y)\,,

where the first term in the r.h.s. corresponds to the integral over the
grey region *above* the blue line while the second over the grey region
*below* the blue line in Fig. `6 <#fig:deltaexpansion>`__.

.. container::
   :name: fig:deltaexpansion

   .. figure:: ../../latex/src/../../latex/src/plots/DeltaExpansion.pdf
   .. :alt: Integration region of the integral in
      Eq. :eq:`eq:testint`. The integral is along the
      red curve defined by the
      :math:`\delta`-function.[fig:deltaexpansion]
      :name: fig:deltaexpansion
      :width: 50.0%

      Integration region of the integral in
      Eq. :eq:`eq:testint`. The integral is along the
      red curve defined by the
      :math:`\delta`-function.[fig:deltaexpansion]

Now I use the following equalities:

.. math::
   \delta(xy-\varepsilon) = \left\{
   \begin{array}{ll}
   \displaystyle \frac{1}{x}\delta\left(y-\frac{\varepsilon}{x}\right)\theta(y-\sqrt{\varepsilon}) &
                                                                         \quad\mbox{integral
                                                                         over
                                                                         $y$,}\\
   \\
   \displaystyle \frac{1}{y}\delta\left(x-\frac{\varepsilon}{y}\right) \theta(x-\sqrt{\varepsilon}) &
                                                                         \quad\mbox{integral
                                                                         over
                                                                         $x$.}
   \end{array}
   \right.

\ The :math:`\theta`-functions arise from the fact that the first
integral has to be done along the upper branch on the red curve while
the second along the lower branch. The two branches are joint at the
point :math:`x=y=\sqrt{\varepsilon}` and thus the integration ranges are
bounded from below by this point. Therefore, I find:

.. math::
   I(\varepsilon) =\int_{\sqrt{\varepsilon}}^1
   \frac{dx}{x}f\left(x,\frac{\varepsilon}{x}\right) + \int_{\sqrt{\varepsilon}}^1 \frac{dy}{y}f\left(\frac{\varepsilon}{y},y\right)\,.

It is crucial to realise that in the first and the second integral the
following conditions hold:
:math:`\varepsilon < \sqrt{\varepsilon}\leq x` and
:math:`\varepsilon < \sqrt{\varepsilon}\leq y`, respectively. Therefore,
in the limit :math:`\varepsilon\rightarrow 0`, the arguments
:math:`\varepsilon/x` and :math:`\varepsilon/y` of the function
:math:`f` will tend to zero. I now add and subtract a term proportional
to :math:`f(0,0)` to both integrals, so that:

.. math::
   I(\varepsilon) =\int_{\sqrt{\varepsilon}}^1
   \frac{dx}{x}\left[f\left(x,\frac{\varepsilon}{x}\right)-f(0,0)\right] + \int_{\sqrt{\varepsilon}}^1 \frac{dy}{y}\left[f\left(\frac{\varepsilon}{y},y\right)-f(0,0)\right]+2f(0,0) \underbrace{\int_{\sqrt{\varepsilon}}^1 \frac{d\xi}{\xi}}_{-\ln\sqrt{\varepsilon}}\,.

Finally, I take the limit for :math:`\varepsilon\rightarrow 0`:

.. math::
   \lim_{\varepsilon\rightarrow 0}I(\varepsilon) =\int_0^1
   \frac{dx}{x}\left[f\left(x,0\right)-f(0,0)\right] + \int_{0}^1 \frac{dy}{y}\left[f\left(0,y\right)-f(0,0)\right]-\ln\left(\varepsilon\right) f(0,0)\,,

that I rewrite as:

.. math:: \lim_{\varepsilon\rightarrow 0}I(\varepsilon) =\int_0^1dx \int_0^1dy\left\{\frac{\delta(y)}{[x]_+}+\frac{\delta(x)}{[y]_+}-\ln\left(\varepsilon\right)\delta(x)\delta(y)\right\}f(x,y)\,.

Comparing the equation above with Eq. :eq:`eq:testint`,
one deduces that:

.. math::
   \delta(xy-\varepsilon)\mathop{\longrightarrow}_{\varepsilon\rightarrow 0}
   \frac{\delta(y)}{[x]_+}+\frac{\delta(x)}{[y]_+}-\ln\left(\varepsilon\right)\delta(x)\delta(y)\,.

Finally, substituting:

.. math:: x\rightarrow \frac{1-x}{x}\,,\quad y\rightarrow \frac{1-z}{z}\,,\quad\mbox{and}\quad\varepsilon\rightarrow\frac{q_T^2}{Q^2}\,,

it is easy to recover
Eq. :eq:`eq:deltaexpansion`. In particular, one
needs to use the fact that:

.. math:: \delta\left(\frac{1-x}{x}\right)=\delta(1-x)\,.

.. container::
   :name: solution-of-an-rge-in-perturbation-theory

   .. rubric:: Solution of an RGE in perturbation theory
      :name: solution-of-an-rge-in-perturbation-theory

In this section we address the question of how a renormalisation group
equation (RGE) can be solved within the framework of perturbation
theory. The relevance of this issue concernes the estimation of
theoretical uncertainties stemming from higher-order corrections. More
specifically, we will show that the so-called *resummation scale*
originates from a perturbative solution of the TMD evolution equations.

In order to define the problem, we consider a generic RGE in QCD:

.. math::
   :label: eq:RGEproto
   \frac{d\ln R}{d\ln \mu} = \gamma(\alpha_s(\mu))\,,

\ where :math:`R` is a generic renormalisation-scheme dependent (and
thus unobservable) quantity, :math:`\mu` is the scale resulting from the
renormalisation procedure, and :math:`\gamma` is the log-free anomalous
dimension that is computable in perturbation theory as a power series in
the strong coupling :math:`\alpha_s` truncated at order :math:`k`:

.. math:: \gamma(\alpha_s(\mu)) = \frac{\alpha_s(\mu)}{4\pi}\sum_{n=0}^k\left(\frac{\alpha_s(\mu)}{4\pi}\right)^n\gamma^{(n)}\,.

Notice that Eq. :eq:`eq:RGEproto` applies to different
relevant quantities in QCD such as the strong coupling itself and to
collinear and TMD distributions. Of course,
Eq. :eq:`eq:RGEproto` is an ordinary differential
equation that defines a family of solutions. In order to identify a
particular solution, we need a boundary condition given by the knowledge
of :math:`R` at some particular scale :math:`\mu_0`:
:math:`R(\mu_0)=R_0`.

.. container::
   :name: sec:strongcoupling

   .. rubric:: The strong coupling
      :name: sec:strongcoupling

The question that we want to address here is how to solve
Eq. :eq:`eq:RGEproto`, along with its boundary
condition :math:`R(\mu_0)=R_0`, within the framework of perturbation
theory. To rephrase it: how do we compute :math:`R(\mu)` to some
perturbative accuracy for a generic scale :math:`\mu` given
Eq. :eq:`eq:RGEproto` and its boundary condition
:math:`R(\mu_0)=R_0`? The first case to be considered, because it
underlies all others, is the strong coupling, *i.e.* we take
:math:`R=a_s=\alpha_s/4\pi` that obeys the RGE [10]_:

.. math::
   :label: eq:RGEalphas
   \frac{d\ln a_s}{d\ln \mu} = \overline{\beta}(a_s(\mu)) =a_s(\mu)\sum_{n=0}^ka_s^n(\mu)\beta^{(n)}\,.

The possibly more natural perturbative approach to this problem is that
to assume that :math:`\alpha_s(\mu)` admits a perturbative expansion
around :math:`a_s(\mu_0)=0` truncated to the same order :math:`k` of the
anomalous dimension:

.. math::
   :label: eq:pertanstz
   a_s(\mu) = a_s(\mu_0) \sum_{n=0}^{k}c_n(\mu) a_s^n(\mu_0)\,.

Plugging Eq. :eq:`eq:pertanstz` into both left- and
right-hand side of Eq. :eq:`eq:RGEalphas`, one finds
the following recursive equations for the first few coefficients
:math:`c_n`:

.. math::
   \begin{array}{lcl}
   \displaystyle \frac{d c_0}{d \ln \mu} &=& 0 \\
   \\
   \displaystyle \frac{d c_1}{d \ln \mu} &=& \beta^{(0)}c_0^2 \\
   \\
   \displaystyle \frac{d c_2}{d \ln \mu} &=& 2 \beta^{(0)} c_0 c_1 +
                                             \beta^{(1)} c_0^3\\
   \\
   \displaystyle \frac{d c_3}{d \ln \mu} &=& \beta^{(0)} (c_1^2 + 2 c_0
                                             c_2) + 3 \beta^{(1)} c_0^2
                                             c_1 + \beta^{(2)} c_0^4 \\
   \\
   \vdots
   \end{array}

Integrating these differential equations between :math:`\mu_0` and
:math:`\mu` and using the boundary condition gives:

.. math::
   :label: eq:pertcoefs
   \begin{array}{lcl}
   c_0 &=& 1 \\
   \\
   c_1 &=& \displaystyle \beta^{(0)}\ln\left(\frac{\mu}{\mu_0}\right) \\
   \\
   c_2 &=& \displaystyle \left(\beta^{(0)}\right)^2\ln^2\left(\frac{\mu}{\mu_0}\right) + \beta^{(1)}\ln\left(\frac{\mu}{\mu_0}\right)\\
   \\
   c_3 &=& \displaystyle \left(\beta^{(0)}\right)^3\ln^3\left(\frac{\mu}{\mu_0}\right) + \frac{5}{2}\beta^{(0)}\beta^{(1)}\ln^2\left(\frac{\mu}{\mu_0}\right) + \beta^{(2)}\ln\left(\frac{\mu}{\mu_0}\right)  \\
   \\
   \vdots
   \end{array}

Computing these coefficients up to order :math:`k` allows one to obtain
the perturbative expansion of :math:`a_s(\mu)` in terms of
:math:`a_s(\mu_0)` up the fixed-order perturbative accuracy
N\ :math:`^{k}`\ LO. However, it should be clear from
Eq. :eq:`eq:pertcoefs` that, for :math:`\mu\gg \mu_0`
or :math:`\mu\ll \mu_0`, :math:`c_n\sim \ln^{n}(\mu/\mu_0)\gg 1`.
Therefore, despite the smallness of :math:`\alpha_s(\mu_0)`, the
presence of such large logarithms is such that the truncated series in
Eq. :eq:`eq:pertanstz`, *i.e.* a fixed-order
computation, may not provide an accurate approximation of
:math:`\alpha_s(\mu)`.

However, nothing prevents us from extending the series in
Eq. :eq:`eq:pertanstz` from :math:`k+1` terms to and
infinite number of terms even though the anomalous dimension in
Eq. :eq:`eq:RGEalphas` is truncated to order
:math:`k`:

.. math::
   :label: eq:pertanstzinf
   a_s(\mu) = a_s(\mu_0) \sum_{n=0}^{\infty}c_n(\mu) a_s^n(\mu_0)\,.

\ The coefficients :math:`c_n` can be computed iteratively exactly as
shown above for any value of :math:`n`. The calculation above can be
used to infer the structure of the coefficients :math:`c_n`:

.. math:: c_n(\mu) = \sum_{i=n-k}^n c_{n,i}\ln^i\left(\frac{\mu}{\mu_0}\right)\,,

where :math:`c_{n,i}` are numerical coefficients given by combinations
of the coefficients :math:`\beta_j`, :math:`j=0,\dots,k`. Using this
equation and rearranging the series finally gives:

.. math::
   :label: eq:pertanstzinfrearr
     a_s(\mu) = a_s(\mu_0) \sum_{l=0}^k a_s^{l}(\mu_0) \sum_{n=0}^{\infty} c_{n+l,n}\left[a_s(\mu_0)\ln\left(\frac{\mu}{\mu_0}\right) \right]^{n}\,.

This rearrangment exposes the resummation procedure: for each power of
:math:`\alpha_s` at which the :math:`\beta`-function is known, all
powers of the kind :math:`(a_s\ln(\mu/\mu_0))^n` are included. In
jargon, one usually says that the knowledge of the anomalous dimension
to order :math:`k` allows N\ :math:`^k`\ LL resummation. The question
that remains to be answered is how do we calculate all the
:math:`c_{n+l,n}` coefficients to effectively enable resummation? This
is usually done by observing that
Eq. :eq:`eq:RGEalphas` with :math:`k=0` can be easily
solved in a closed form and gives:

.. math::
   :label: eq:llsol
   a_s^{\rm LL}(\mu) = \frac{a_s(\mu_0)}{1-\beta^{(0)}a_s(\mu_0)\ln\left(\frac{\mu}{\mu_0}\right)}\,,

which provides the LL resummation of the running coupling. Comparing
this result with
Eq. :eq:`eq:pertanstzinfrearr` with
:math:`k=0` tells us that:

.. math:: c_{n,n} = \left(\beta^{(0)}\right)^n\,,

which is also consistent with Eq. :eq:`eq:pertcoefs`.
In order to go further, one needs to include the term :math:`k=1` in
Eq. :eq:`eq:RGEalphas` and integrate the differential
equation between :math:`\mu_0` and :math:`\mu`. The result is a
transcendental equation in :math:`a_s`:

.. math::
   :label: eq:transcend
   \frac1{\beta^{(0)}}\left(-\frac{1}{a_s(\mu)}+\frac{1}{a_s(\mu_0)}\right) + \frac{b_1}{\beta^{(0)}}\ln\left(\frac{a_s(\mu_0) (1+b_1a_s(\mu))}{a_s(\mu)  (1+b_1a_s(\mu_0))}\right)=\ln\left(\frac{\mu}{\mu_0}\right)

with :math:`b_n=\beta^{(n)}/\beta^{(0)}`. Of course, the solution of
this transcendental equation solves exactly the differential equation in
Eq. :eq:`eq:RGEalphas`. However, in order to obtain an
analytic equation for :math:`a_s(\mu)` we need to resort on perturbation
theory. This is done by first rearranging the above equation as:

.. math:: \frac{a_s(\mu_0)}{a_s(\mu)}=1-a_s(\mu_0)\beta^{(0)}\ln\left(\frac{\mu}{\mu_0}\right)  + a_s(\mu_0) b_1\ln\left(\frac{a_s(\mu_0) (1+b_1a_s(\mu))}{a_s(\mu)  (1+b_1a_s(\mu_0))}\right)\,.

This equation tells us that the second logarithm appearing in the
r.h.s., being proportional to :math:`\alpha_s` can be expanded and
truncated to leading order. Notice that this operation is allowed
because no potentially large logarithms of the scales are present.
Therefore, we can use Eq. :eq:`eq:llsol` inside this term,
obtaining:

.. math::
   :label: eq:firstexp
   \ln\left(\frac{a_s(\mu_0) (1+b_1a_s(\mu))}{a_s(\mu)  (1+b_1a_s(\mu_0))}\right)=\ln\left(\frac{a_s(\mu_0)}{a_s^{\rm LL}(\mu)}\right)+\mathcal{O}(a_s)\,,

where we have also used the fact that:

.. math::
   :label: eq:approxLL
   a_s(\mu) = a_s^{\rm LL}(\mu)\left[1+\mathcal{O}(a_s) \right]\,,

that descends directly from
Eq. :eq:`eq:pertanstzinfrearr`. This leads us
to:

.. math:: a_s(\mu)=\frac{a_s^{\rm LL}(\mu)}{ 1  + a_s^{\rm LL}(\mu) b_1 \ln\left(\frac{a_s(\mu_0)}{a_s^{\rm LL}(\mu)}\right)}+\mathcal{O}(a_s^3)\,

that can be further expanded, finally obtaining:

.. math::
   :label: eq:ansolut
   a_s^{\rm NLL}(\mu)=a_s^{\rm LL}(\mu) \left[1 - b_1  a_s^{\rm LL}(\mu) \ln\left(\frac{a_s(\mu_0)}{a_s^{\rm LL}(\mu)}\right)\right]\,,

which gives the NLL solution for the evolution of the strong coupling
(see also Ref. (Del Debbio et al. 2007)). The same kind of procedure can
be applied iteratively to obtain perturbatively more accurate solutions.
In the following we will derive the NNLL-accurate solution for the
running of the strong coupling that we will need below. This can be
achieved by including :math:`k=2` term in
Eq. :eq:`eq:RGEalphas` and integrating the resulting
differential equation between :math:`\mu_0` and :math:`\mu`. This
results in the following transcendental equation:

.. math::
   :label: eq:transc2
   \begin{array}{rcl}
   \displaystyle \beta^{(0)}\ln\left(\frac{\mu}{\mu_0}\right) &=&\displaystyle \left. -\frac{1}{a_s}\
   -b_1\ln(a_s)
                                                    +\mbox{const.}+a_s(b_1^2-b_2)+\mathcal{O}(a_s^2)\right|_{a_s(\mu_0)}^{a_s(\mu)}\\
   \\
   &=&\displaystyle
       -\frac{1}{a_s(\mu)}+\frac{1}{a_s(\mu_0)}+b_1\ln\left(\frac{a_s(\mu_0)}{a_s^{\rm
       NLL}(\mu)}\right)
       + (b_1^2-b_2)(a_s^{\rm LL}(\mu)-a_s(\mu_0)) +\mathcal{O}(a_s^2)\,.
   \end{array}

where we have replaced :math:`a_s(\mu)` with :math:`a_s^{\rm LL}` and
:math:`a_s^{\rm
  NLL}` where appropriate using Eq. :eq:`eq:approxLL`
and:

.. math::
   :label: eq:approxNLL
   a_s(\mu) = a_s^{\rm NLL}(\mu)\left[1+\mathcal{O}(a_s^2) \right]\,.

Eq. :eq:`eq:transc2` can finally be solved for
:math:`a_s(\mu)` obtaining:

.. math::
   :label: eq:ansolutnnll
   \begin{array}{rcl}
     a_s^{\rm NNLL}(\mu)&=&\displaystyle a_s^{\rm LL}(\mu)\bigg[1+b_1 a_s^{\rm LL}(\mu)\ln\left(\frac{a_s^{\rm
           LL}(\mu)}{a_s(\mu_0)}\right)+\displaystyle (b_2-b_1^2) a_s^{\rm LL}(\mu) \left[a_s^{\rm
         LL}(\mu)-a_s(\mu_0)\right]\\
   \\
   &+&\displaystyle b_1^2 \left(a_s^{\rm LL}(\mu)\right)^2 \ln\left(\frac{a_s^{\rm
           LL}(\mu)}{a_s(\mu_0)}\right)\left(1+\ln\left(\frac{a_s^{\rm
           LL}(\mu)}{a_s(\mu_0)}\right)\right) \bigg]\\
   \\
   &=&\displaystyle a_s^{\rm LL}(\mu)\left[1+b_1 \left[a_s^{\rm NLL}(\mu)+b_1 \left(a_s^{\rm LL}(\mu)\right)^2\right]\ln\left(\frac{a_s^{\rm
           LL}(\mu)}{a_s(\mu_0)}\right)+\displaystyle (b_2-b_1^2) a_s^{\rm LL}(\mu) \left[a_s^{\rm
         LL}(\mu)-a_s(\mu_0)\right]\right]\,,
   \end{array}

*i.e.* the NNLL-accurate running of the strong coupling.

It is now important to stress an important aspect that is central to
this section, *i.e.* the difference between, *e.g.*
Eqs. :eq:`eq:transcend`
and :eq:`eq:ansolut`.
Eq. :eq:`eq:ansolut` has been obtained from
Eq. :eq:`eq:transcend` by neglecting genuinely
subleading corrections. This makes the two solutions equivalent from the
point of view of accuracy. However, while
Eq. :eq:`eq:transcend` satisfies the evolution
equation Eq. :eq:`eq:RGEalphas` with :math:`k=1`
*exactly*, Eq. :eq:`eq:ansolut` violates it by
subleading terms. In other words, one finds:

.. math::
   :label: eq:ansolNLL
   \frac{d \ln a_s^{\rm NLL}(\mu)}{d\ln\mu}=\beta^{(0)}a_s^{\rm NLL}(\mu) +\beta^{(1)} (a_s^{\rm NLL}(\mu))^2 +\mathcal{O}(a_s^3)\,.

Therefore, the neglect of subleading terms required to go from
Eq. :eq:`eq:transcend` to
Eq. :eq:`eq:ansolut` has the effect of introducing
subleading terms in the evolution equation that is no longer exactly
satisfied. We point out that these spurious subleading terms originates
from insisting on obtaining a closed-form analytic expression for the
evolution of the coupling.

Eq. :eq:`eq:ansolut` has the additional (perhaps
undesired) feature that it makes the evolution of the coupling
*non-conservative*. In practice, this means that evolving
:math:`\alpha_s` from :math:`\mu_0` to :math:`\mu` using
Eq. :eq:`eq:ansolut` and then back from :math:`\mu` to
:math:`\mu_0` one does not recover the initial value. This is clearly
not the case if one uses the exact solution deriving from
Eq. :eq:`eq:transcend`. To show this,
Fig. `7 <#fig:Alphas>`__ displays the behaviour of the NLL analytical
solution, Eq. :eq:`eq:ansolNLL`, with :math:`n_f=5`
active flavours evolved from :math:`\alpha_s(M_Z)=0.118` down to 1 GeV
and then up back to :math:`M_Z`. As clear, the backward and forward
evolution curves do not coincide leading to a mismatch of the value of
:math:`\alpha_s(M_Z)`.

.. container::
   :name: fig:Alphas

   .. figure:: ../../latex/src/../../latex/src/plots/Alphas.pdf
   .. :alt: Non-conservativity of the NLL analytic evolution of the
      strong coupling.[fig:Alphas]
      :name: fig:Alphas
      :width: 60.0%

      Non-conservativity of the NLL analytic evolution of the strong
      coupling.[fig:Alphas]

Of course, one can argue that this effect is subleading and thus
allowed. On the other hand, one may wonder whether a violation of the
evolution equation introduced with the sole purpose of obtaining an
analytic expression for the evolution is justified. In fact, it appears
to be a better practice to employ an evolution that obeys the evolution
equation exactly order by order in perturbation theory. The reason is
simply a better control on the theoretical uncertainty related to
missing higher-order corrections. To be more specific, the subleading
terms introduced to obtain Eq. :eq:`eq:ansolut` would
affect the computation of an observable at fixed-order in perturbation
theory in the exact same way as a renormalisation scale variation would
do. Therefore, it appears more transparent to use an exact solution to
the evolution equation and leave the estimate of theoretical
uncertainties to renormalisation-scale variations.

.. container::
   :name: dglap-evolution

   .. rubric:: DGLAP evolution
      :name: dglap-evolution

As a second application of Eq. :eq:`eq:RGEproto`, we
consider the DGLAP evolution in which the quantity :math:`R` is
identified with a collinear distribution :math:`f`. In order to simplify
the discussion, we consider the Mellin moments of a non-singlet
distribution whose evolution equation can be written as exactly as in
Eq. :eq:`eq:RGEproto` with :math:`R=f` and
:math:`\gamma` the appropriate anomalous dimension. Our aim is to find a
solution to this equation in perturbation theory by computing :math:`f`
at some final scale :math:`\mu` when knowing :math:`f` at some initial
scale :math:`\mu_0`. Of course, the solution to the evolution equation
is:

.. math::
   :label: eq:dglapsol
   \begin{array}{rcl}
   f(\mu)&=&\displaystyle
   \exp\left[\int_{\mu_0}^{\mu}d\ln\mu'\,\gamma(a_s(\mu'))\right]f(\mu_0)\\
   \\
   &=&\displaystyle
       \exp\left[\sum_{n=0}^k\gamma^{(n)}\int_{\mu_0}^{\mu}d\ln\mu'\,a_s^{n+1}(\mu')\right]f(\mu_0)\\
   \\
   &=&\displaystyle f(\mu_0) \prod_{n=0}^k\exp\left[\gamma^{(n)}I_n\right]f(\mu_0)\,,
   \end{array}

\ with:

.. math:: I_n=\int_{a_s(\mu_0)}^{a_s(\mu)}da_s\left(\frac{a_s^{n}}{\overline{\beta}(a_s)}\right)\,.

In the second line of Eq. :eq:`eq:dglapsol` we have
used the perturbative expansion of the anomalous dimension
:math:`\gamma` and in the third line made a change of variable in the
integral using Eq. :eq:`eq:RGEalphas`. It is now
opportune to analyse the integrals :math:`I_n`. The first question that
we ask is whether the perturbative expansion of :math:`\gamma` and
:math:`\overline{\beta}` are related. More specifically, given a
:math:`k+1` terms in the expansion of :math:`\overline{\beta}`, is there
any bound on the number of terms in the expansion of :math:`\gamma` in
terms of logarithmic accuracy? In order to answer this question, it is
necessary to first understand what kind of logarithms
Eq. :eq:`eq:dglapsol` is actually resumming. Indeed,
one may wonder whether the evolution of :math:`\alpha_s` alone is enough
to rusum all the large logarithms. If so, we would be enabled to expand
and truncate the exponential to the desired perturbative accuracy
because higher orders would be suppressed by powers of :math:`\alpha_s`.
The answer to this question boils down to determining whether
:math:`I_n` is logarithmically enhanced for some value of :math:`n`. We
assume that :math:`\overline\beta` is truncated at some power
:math:`k'+1`, with :math:`k'` generally different from :math:`k` in
Eq. :eq:`eq:dglapsol`. To give a quantitative argument
let us take :math:`k'=2`, such that

.. math:: \overline{\beta}(a_s) =\beta^{(0)}a_s\left[1+b_1a_s\right]\,.

We then compute :math:`I_0`, :math:`I_1`, and :math:`I_2`:

.. math::
   :label: eq:intsI
   \begin{array}{rcl}
     I_0&=&\displaystyle\frac1{\beta^{(0)}}\int_{a_s(\mu_0)}^{a_s(\mu)}da_s\left(\frac{1}{a_s+b_1a_s^2}\right)=\frac1{\beta^{(0)}}\left[\ln
            a_s-\ln(1+b_1a_s)\right]\bigg|_{a_s(\mu_0)}^{a_s(\mu)}=\frac{1}{\beta^{(0)}}\ln\left(\frac{a_s(\mu)}{a_s(\mu_0)}\right)+\mathcal{O}(a_s)\,,\\
   \\
     I_1&=&\displaystyle \frac1{\beta^{(0)}}\int_{a_s(\mu_0)}^{a_s(\mu)}da_s\left(\frac{1}{1+b_1a_s}\right)=\frac1{\beta^{(0)}}\left[\frac{\ln(1+b_1a_s)}{b_1}\right]\bigg|_{a_s(\mu_0)}^{a_s(\mu)}=\mathcal{O}(a_s) \,,\\
   \\
     I_2&=&\displaystyle\frac1{\beta^{(0)}}\int_{a_s(\mu_0)}^{a_s(\mu)}da_s\left(\frac{a_s}{1+b_1a_s}\right)=\frac1{\beta^{(0)}}\left[\frac{a_s}{b_1}-\frac{\ln(1+b_1a_s)}{b_1^2}\right]\bigg|_{a_s(\mu_0)}^{a_s(\mu)}=\mathcal{O}(a_s^2)\,.
   \end{array}

In general one finds that, for :math:`n\geq1`,
:math:`I_n=\mathcal{O}(a_s^n)`. Therefore, the exponential in
Eq. :eq:`eq:dglapsol` can safely be expanded:

.. math:: \exp\left[\gamma^{(n)}I_n\right]=1+\gamma^{((n))}\times\mathcal{O}(a_s^n)\,,\quad n\geq 1\,.

Somewhat surprisingly, this reveals that the inclusion of perturbative
corrections beyond leading order to the anomalous dimension
:math:`\gamma` does not directly improve the logarithmic accuracy of the
resummation provided by the evolution. The exception to this pattern is
:math:`I_0` that cannot be (fully) expanded due to the presence of a
logarithm of :math:`a_s(\mu)/a_s(\mu_0)` that, as can be see from
Eq. :eq:`eq:llsol`, is potentially large:

.. math::
   :label: eq:expres
   \exp\left[\gamma^{(0)}I_0\right]=\exp \left[\frac{\gamma^{(0)}}{\beta^{(0)}}\ln\left(\frac{a_s(\mu)}{a_s(\mu_0)}\right)\right]+\mathcal{O}(a_s)=\left(\frac{a_s(\mu)}{a_s(\mu_0)}\right)^{\frac{\gamma^{(0)}}{\beta^{(0)}}}+\mathcal{O}(a_s)\,.

This factor alone provides the resummation up to some logarithmic
accuracy determined uniquely by the evolution of :math:`a_s(\mu)`.
Finally one has:

.. math::
   :label: eq:dglapanal
   f(\mu) = \displaystyle \left[1+\mathcal{O}(a_s)\right]\left(\frac{a_s(\mu)}{a_s(\mu_0)}\right)^{\frac{\gamma^{(0)}}{\beta^{(0)}}}f(\mu_0) \,.

However, the inclusion of correction terms in the square brackets in the
r.h.s. of Eq. :eq:`eq:dglapanal` is still required to
attain some given logarithmic accuracy. Since the exponential in
Eq. :eq:`eq:expres` resums single logarithms of the kind
:math:`a_s^m\ln ^m(\mu/\mu_0)`, a term proportional to :math:`a_s^n`
(without any logarithm) belongs to the tower
:math:`a_s^{m+n}\ln ^m(\mu/\mu_0)` and thus contributes to
N\ :math:`^n`\ LL. To make an example, suppose we want to achieve NLL
resummation. To do so, we need to include the :math:`\mathcal{O}(a_s)`
in the squared bracket of Eq. :eq:`eq:dglapanal`
because that contributes to NLL. That term can be easily computed by
retaining the :math:`\mathcal{O}(a_s)` terms in
Eq. :eq:`eq:intsI`. Moreover, :math:`a_s(\mu)` in the last
factor has to be computed using the NLL evolution given in
Eq. :eq:`eq:ansolNLL`. This finally gives:

.. math::
   f^{\rm NLL}(\mu) = \displaystyle\left[1+\frac1{\beta^{(0)}}\left(\gamma^{(1)}-b_1
       \gamma^{(0)}\right) \left(a_s^{\rm LL}(\mu)-a_s(\mu_0)\right)\right]
   \left(\frac{a_s^{\rm NLL}(\mu)}{a_s(\mu_0)}\right)^{\frac{\gamma^{(0)}}{\beta^{(0)}}} f(\mu_0) \,.

Also note that in the squared bracket we are entitled to replaced
:math:`a_s` with :math:`a_s^{\rm LL}` becase the difference is
:math:`\mathcal{O}(a_s^2)` and thus NNLL.

In order to make a step towards the notation typically used in
:math:`q_T` resummation, we introduce the quantity:

.. math::
   :label: eq:lambdadef
   \lambda \equiv \lambda(\mu) = a_s(\mu_0)\beta^{(0)}\ln\left(\frac{\mu}{\mu_0}\right)\,,

that allows us to write the LL and NLL RGE solution for :math:`a_s` as:

.. math:: a_s^{\rm LL}(\mu) = a_s(\mu_0)\frac{1}{1-\lambda}\,,

and:

.. math:: a_s^{\rm NLL}(\mu)=a_s(\mu_0)\frac{1}{1-\lambda}\left[1 - a_s(\mu_0) \frac{b_1 \ln\left(1-\lambda\right)}{1-\lambda}\right]\,,

such that:

.. math::
   \begin{array}{rcl}
     f^{\rm NLL}(\mu) &=&\displaystyle \displaystyle \left[1+a_s(\mu_0)\frac1{\beta^{(0)}}\left(\gamma^{(1)}-b_1
                               \gamma^{(0)}\right)
                               \frac{\lambda}{1-\lambda}\right]\\
     \\
                      &\times&\displaystyle                        \exp\left[-\frac{\gamma^{(0)}}{\beta^{(0)}}\ln\left(1-\lambda\right)-a_s(\mu_0) \frac{\gamma^{(0)}b_1}{\beta^{(0)}}\frac{
                          \ln\left(1-\lambda\right)}{1-\lambda}\right]f(\mu_0)\,.
   \end{array}

We now introduce the functions:

.. math::
   \begin{array}{rcl}
   g^{(0)}(\lambda) &=& \displaystyle  1+a_s(\mu_0)\frac1{\beta^{(0)}}\left(\gamma^{(1)}-b_1
                               \gamma^{(0)}\right)
                               \frac{\lambda}{1-\lambda}\,,\\
   \\
   g^{(1)}(\lambda) &=& \displaystyle
                        -\frac{\gamma^{(0)}}{\beta^{(0)}}\ln\left(1-\lambda\right)\,,\\
   \\
   g^{(2)}(\lambda) &=& \displaystyle -\frac{\gamma^{(0)}b_1}{\beta^{(0)}}\frac{
                          \ln\left(1-\lambda\right)}{1-\lambda}\,,
   \end{array}

such that:

.. math::
   :label: eq:DGLAPalaSoftGlue
     f^{\rm NLL}(\mu) = g^{(0)}(\lambda)\exp\left[g^{(1)}(\lambda)+a_s(\mu_0) g^{(2)}(\lambda)\right]f(\mu_0) \,.

This procedure can be extended to N\ :math:`^k`\ LL accuracy just by
introducing the appropriate :math:`g^{(n)}` functions, with
:math:`n\leq k+1`, in the exponential and including corrections up to
:math:`\mathcal{O}(a_s^k)` in :math:`g^{(0)}`. The notation is purposely
reminiscent of soft-gluon and :math:`q_T` resummation.

We can use the parallel with soft-gluon resummation to introduce the
analogous of a “resummation” scale. This is done by introducing an
*arbitrary* scale :math:`\mu_{\rm R}\sim \mu_0` along with the trivial
equality:

.. math::
   \ln\left(\frac{\mu}{\mu_0}\right) = \ln\left(\frac{\mu}{\mu_{\rm R}}\right)
   + \ln\left(\frac{\mu_{\rm R}}{\mu_0}\right)\,,

\ such that :math:`\ln\left(\mu_{\rm R}/\mu_0\right)\sim 1`. Using
Eq. :eq:`eq:lambdadef`, one has:

.. math::
   \lambda =
     a_s(\mu_0)\beta^{(0)}\ln\left(\frac{\mu}{\mu_{\rm R}}\right)+a_s(\mu_0)\beta^{(0)}\ln\left(\frac{\mu_{\rm
           R}}{\mu_0}\right)=\overline{\lambda}+a_s(\mu_0)\beta^{(0)}\ln\left(\frac{\mu_{\rm
           R}}{\mu_0}\right)\,,

where the second term in the r.h.s. is not enhanced by a large logarithm
and can thus be treated perturbatively. Writing the functions
:math:`g^{(n)}` in terms of :math:`\overline{\lambda}` and neglecting
subleading terms gives the new functions:

.. math::
   \begin{array}{rcl}
   \overline{g}^{(0)}(\overline{\lambda}) &=& \displaystyle  1+a_s(\mu_0)\frac1{\beta^{(0)}}\left(\gamma^{(1)}-b_1
                               \gamma^{(0)}\right)
                               \frac{\overline{\lambda}}{1-\overline{\lambda}}\,,\\
   \\
   \overline{g}^{(1)}(\overline{\lambda}) &=& \displaystyle
                        -\frac{\gamma^{(0)}}{\beta^{(0)}}\ln\left(1-\overline{\lambda}\right)\,,\\
   \\
   \overline{g}^{(2)}(\overline{\lambda}) &=& \displaystyle -\frac{\gamma^{(0)}}{1-\overline{\lambda}}\left[\frac{b_1}{\beta^{(0)}}\ln\left(1-\overline{\lambda}\right)-\ln\left(\frac{\mu_{\rm
           R}}{\mu_0}\right)\right]\,,
   \end{array}

so that:

.. math::
   :label: eq:DGLAPalaSoftGlueMuR
     f^{\rm NLL}(\mu) = \overline{g}^{(0)}(\overline{\lambda})\exp\left[\overline{g}^{(1)}(\overline{\lambda})+a_s(\mu_0) \overline{g}^{(2)}(\overline{\lambda})\right]f(\mu_0) \,.

Eqs. :eq:`eq:DGLAPalaSoftGlueMuR`
and :eq:`eq:DGLAPalaSoftGlue` are both NLL
accurate and only differ by subleading terms. In fact,
Eq. :eq:`eq:DGLAPalaSoftGlueMuR` generalises
Eq. :eq:`eq:DGLAPalaSoftGlue` because the
latter reduces to the former for :math:`\mu_{\rm R}=\mu_0`. Therefore,
one should be able to use modest variations of :math:`\mu_{\rm R}`
around :math:`\mu_0` to estimate the possible impact of higher-order
corrections to the anomalous dimension.

Eq. :eq:`eq:DGLAPalaSoftGlueMuR` provides a
closed-form, fully analytical solution to the RGE,
Eq. :eq:`eq:RGEproto`, for collinear distributions.
However, exactly like in the case of the analytical solution for
:math:`a_s`, Eq. :eq:`eq:dglapanal` will satisfy
Eq. :eq:`eq:RGEproto` only up to subleading terms with
a consequent non-conservative pattern of the evolution. Conversely, a
numerical evaluation of the solution in
Eq. :eq:`eq:dglapsol` will obey exactly
Eq. :eq:`eq:RGEproto` and the evolution will be
conservative.

.. container::
   :name: the-sudakov-form-factor

   .. rubric:: The Sudakov form factor
      :name: the-sudakov-form-factor

Now, we turn to address the question of solving the TMD evolution
equations in Eq. :eq:`eq:eveqs` in perturbation theory. We
have already derived the solution,
Eq. :eq:`eq:solution2`, in terms of the relevant
anomalous dimensions. However, this solution contains an integral. A
possibility is that of computing that integral numerically. This ensures
that the results obeys *exactly* the evolution equations in
Eq. :eq:`eq:eveqs` and the evolution is conservative.
Another possibility is that of attempting to solve the integral
analytically. This turns out to be possible provided that one uses the
appropriate analytic solution for the running of :math:`a_s`. As
discussed above, this comes at the price that the result violates the
equations in Eq. :eq:`eq:eveqs` by subleading terms. As we
will see below, this very violation gives origin to the resummation
scale.

In order to derive a fully analytic expression for the expression in
Eq. :eq:`eq:solution2` we need to consider the
perturbative expansion of the relevant anomalous dimensions. This is
done in Eq. :eq:`eq:Sudexp` for the square of the
evolution operator that is relevant to Drell-Yan production. In order to
match the notation of this section. we identify :math:`\mu_b=\mu` and
:math:`Q=\mu_0` in Eq. :eq:`eq:Sudexp`. In addition, we
limit ourselves to consider the NNLL accuracy that means truncating the
expansion at :math:`n=2` for :math:`\gamma_K` and at :math:`n=1` for
:math:`K` and :math:`\gamma_F`. This finally gives:

.. math::
   \begin{array}{rcl}
     R^{\rm NNLL} &=&\displaystyle  \exp\bigg\{-\sum_{n=0}^1 K^{(n,0)}
           \left[a_s^{{\rm N}^{1-n}{\rm
           LL}}(\mu)\right]^{n+1}\ln\frac{\mu}{\mu_0}-\sum_{n=0}^1
           \gamma_F^{(n)}\int_{\mu_0}^{\mu}\frac{d\mu'}{\mu'}\left[a_s^{{\rm
           N}^{1-n}{\rm LL}}(\mu')\right]^{n+1}\\
   \\
    &-&\displaystyle \sum_{n=0}^2
        \gamma_K^{(n)}\int_{\mu_0}^{\mu}\frac{d\mu'}{\mu'}\left[a_s^{{\rm
        N}^{2-n}{\rm
        LL}}(\mu')\right]^{n+1}\ln\frac{\mu'}{\mu_0}\bigg\}\,,
   \end{array}

\ where we have also used the fact that:

.. math:: a_s = a_s^{{\rm N}^{k}{\rm LL}}\left[1+\mathcal{O}(a_s^k)\right]\,.

It is convenient to express the above equation in terms of the variable
:math:`\lambda` defined in Eq. :eq:`eq:lambdadef`.
This leads to:

.. math::
   \begin{array}{rcl}
     R^{\rm NNLL} &=&\displaystyle  \exp\Bigg\{-\frac{1}{\beta^{(0)}a_s(\mu_0)}\sum_{n=0}^1\left( K^{(n,0)}
           \lambda\left[a_s^{{\rm N}^{1-n}{\rm
           LL}}(\lambda)\right]^{n+1}+
           \gamma_F^{(n)}\int_{0}^{\lambda}d\lambda'\left[a_s^{{\rm
           N}^{1-n}{\rm LL}}(\lambda')\right]^{n+1}\right)\\
   \\
    &-&\displaystyle \left(\frac{1}{\beta^{(0)}a_s(\mu_0)}\right)^2\sum_{n=0}^2
        \gamma_K^{(n)}\int_{0}^{\lambda}d\lambda'\,\lambda'\left[a_s^{{\rm
        N}^{2-n}{\rm
        LL}}(\lambda')\right]^{n+1}\Bigg\}\,,
   \end{array}

where :math:`\lambda'=\lambda(\mu')` and
:math:`a_s(\lambda)\equiv a_s(\mu)`. We write it more explicitly as:

.. math::
   \begin{array}{rcl}
     R^{\rm NNLL} &=&\displaystyle  \exp\Bigg\{-\frac{1}{\beta^{(0)}a_s(\mu_0)}\left( K^{(0,0)}
           \lambda a_s^{{\rm NLL}}(\lambda)+
           \gamma_F^{(0)}\int_{0}^{\lambda}d\lambda' a_s^{{\rm NLL}}(\lambda')+ K^{(1,0)}
           \lambda\left[a_s^{{\rm
           LL}}(\lambda)\right]^{2}+
           \gamma_F^{(1)}\int_{0}^{\lambda}d\lambda'\left[a_s^{{\rm LL}}(\lambda')\right]^{2}\right)\\
   \\
    &-&\displaystyle \left(\frac{1}{\beta^{(0)}a_s(\mu_0)}\right)^2\left(
        \gamma_K^{(0)}\int_{0}^{\lambda}d\lambda'\,\lambda' a_s^{{\rm
        NNLL}}(\lambda')+\gamma_K^{(1)}\int_{0}^{\lambda}d\lambda'\,\lambda'\left[a_s^{{\rm
        NLL}}(\lambda')\right]^{2}+\gamma_K^{(2)}\int_{0}^{\lambda}d\lambda'\,\lambda'\left[a_s^{{\rm
        LL}}(\lambda')\right]^{3}\right)\Bigg\}\,.
   \end{array}

We now need to compute the integrals involving the strong coupling
:math:`a_s`. This can be done using the analytic expressions derived in
Sect. `2.1 <#sec:strongcoupling>`__.The task is simplified by the fact
that :math:`a_s^{{\rm N}^k{\rm LL}}`, for :math:`k>0`, can be written in
terms of :math:`a_s^{\rm LL}` (see Eqs. :eq:`eq:ansolut`
and :eq:`eq:ansolutnnll`) so that we can exploit the
relations:

.. math::
   \int_0^\lambda d\lambda'\,f\left(a_s^{\rm LL}(\lambda')\right) =
   a_s(\mu_0)\int_{a_s(\mu_0)}^{a_s^{\rm LL}(\lambda)}\frac{da}{a^2}\,f(a)\,,

.. math::
   \int_0^\lambda d\lambda'\,\lambda'\,f\left(a_s^{\rm LL}(\lambda')\right) =
   \int_0^\lambda d\lambda'\,f\left(a_s^{\rm LL}(\lambda')\right)-a_s^2(\mu_0)\int_{a_s(\mu_0)}^{a_s^{\rm LL}(\lambda)}\frac{da}{a^3}\,f(a)\,.

.. container::
   :name: iterative-solution-of-the-strong-coupling-running

   .. rubric:: Iterative solution of the strong coupling running
      :name: iterative-solution-of-the-strong-coupling-running

In this section, I will derive a general formula to compute iteratively
the running of the strong coupling. In order to simplify the notation, I
will refer to the strong coupling computed at some generic scale as to
:math:`a`, while I will call :math:`a_0` the value of coupling at the
reference scale. The corresponding RGE will then be:

.. math::
   :label: eq:allorders
   \frac{d\ln a}{d\ln \mu} = \beta(a)\,,

\ with:

.. math:: \beta(a) = a\sum_{n=0}^\infty \beta_n a^n\,.

If one truncates the :math:`\beta` function at the first order, the
resulting evolution equation reads:

.. math:: \frac{d\ln a_{\rm LL}}{d\ln \mu} = a_{\rm LL}\beta_0\,,

that admits a simple closed-form exact solution:

.. math::
   :label: eq:leadingorder
   a_{\rm LL}(\mu) = \frac{a_0}{1-a_0\beta_0\ln\left(\frac{\mu}{\mu_0}\right) }\,.

Since we know how to write :math:`a_{\rm LL}` in a closed form, we can
combine Eqs. :eq:`eq:allorders`
and :eq:`eq:leadingorder` as follows:

.. math:: \frac{d\ln a}{d\ln a_{\rm LL}} =\frac{a}{a_{\rm LL}}\sum_{n=0}^\infty b_na^n\,,\quad\mbox{with}\quad b_n=\frac{\beta_n}{\beta_0}\,.

Now we can attempt to solve this equation whose general solution is:

.. math::
   :label: eq:intrun
   \int_{a_0}^{a_{\rm LL}}\frac{dy}{y^2}=\int_{a_0}^{a}\frac{dx}{x^2}\left(\sum_{n=0}^\infty b_nx^n\right)^{-1}\,.

Assuming that :math:`x\ll 1` (as it should be in the perturbative
regime), one can write:

.. math:: \left(\sum_{n=0}^\infty b_nx^n\right)^{-1}=\sum_{n=0}^\infty \overline{b}_nx^n\,,

where the coefficients :math:`\overline{b}_n` can be iteratively
determined in terms of the coefficients :math:`{b}_n`. The terms up to
:math:`x^3` read:

.. math::
   \begin{array}{rcl}
   \overline{b}_1  &=& - b_1\,,\\
   \overline{b}_2  &=& - b_2+b_1^2\,,\\
   \overline{b}_3  &=& - b_3+2b_2b_1-b_1^3\,.
   \end{array}

Plugging this expansion into Eq. :eq:`eq:intrun` and
solving the integrals gives:

.. math:: \frac{1}{a_{\rm LL}} = \frac{1}{a}-\overline{b}_1\ln\left(\frac{a}{a_0}\right)-\sum_{k=1}^\infty\frac{\overline{b}_{k+1}}{k}\left(a^k-a_0^k\right)\,.

This is a transcendental equation in :math:`a` that can be solved
iteratively assuming that :math:`a,a_0\ll 1`. Upon this assumption, the
terms in the r.h.s. of the equation above are increasingly smaller.
Therefore, at the first order one has:

.. math::
   \frac{1}{a_{\rm LL}} = \frac{1}{a}\quad\Longrightarrow\quad a = a_{\rm
   LL}\,,

as expected. At second order instead:

.. math:: \frac{1}{a_{\rm LL}} = \frac{1}{a}-\overline{b}_1\ln\left(\frac{a}{a_0}\right)\,.

This is a transcendental equation that can be solved by observing that
the second term in the r.h.s. is subleading w.r.t. the first and thus
one can plug into that term the first order solution, that is:

.. math::
   \frac{1}{a_{\rm LL}} = \frac{1}{a}-\overline{b}_1\ln\left(\frac{a_{\rm
       LL}}{a_0}\right)\,.

This equation can now be resolved for :math:`a` obtaining:

.. math::
   a = a_{\rm NLL}=\frac{a_{\rm LL}}{\displaystyle 1+a_{\rm LL}\overline{b}_1\ln\left(\frac{a_{\rm
             LL}}{a_0}\right)}\,.

One can keep iterating in the same way finally obtaining the general
recursive formula for the analytic N\ :math:`^{n}`\ LL solution:

.. math::
   :label: eq:masterf
   a_{{\rm N}^{n}{\rm LL}}=a_{\rm LL}\left\{1+a_{\rm LL}\left[\overline{b}_1\ln\left(\frac{a_{{\rm
             N}^{n-1}{\rm LL}}}{a_0}\right)+\sum_{k=1}^{n-1}\frac{\overline{b}_{k+1}}{k}\left(a_{{\rm
             N}^{n-k-1}{\rm
             LL}}^k-a_0^k\right)\right]\right\}^{-1}+\mathcal{O}(a^{n+2})\,,\quad n\geq 1\,.

Note also that the curly bracket, that appears in the denominator, can
be further expanded and truncated to the relevant order because the
single terms appearing there are increasingly subleading.

At this stage it is important to highlight a feature of the analytic
solution of running of the strong coupling. Due to its iterative nature,
the N\ :math:`^{n}`\ LL solution will only depend on :math:`a_{\rm LL}`.
This allows us to write the solution to all other RGEs, such as the
DGLAP evolution equations and the Collins-Soper equation, in terms of
:math:`a_{\rm LL}` only.

Additionally, it is interesting to observe that there are *two*
resummation structures. The first is contained in :math:`a_{\rm LL}` and
given by the geometric summation in
Eq. :eq:`eq:leadingorder`. The second is
:math:`\ln\left(\frac{a_{{\rm N}^{n-1}{\rm
        LL}}}{a_0}\right)`. This structure first appears at NLL as:

.. math:: a_{{\rm NLL}}=a_{\rm LL}\left[1+a_{\rm LL}\overline{b}_1\ln\left(\frac{a_{\rm LL}}{a_0}\right)\right]^{-1}+\mathcal{O}(a^3)=a_{\rm LL}\left[1-a_{\rm LL}\overline{b}_1\ln\left(\frac{a_{\rm LL}}{a_0}\right)\right]+\mathcal{O}(a^3)\,.

When going further, the content of the logarithm will be the coupling
computed at higher logarithmic orders. However, beyond LL, the logarithm
can be expanded around :math:`a_{\rm LL}=0` finally leading to
:math:`\ln\left(\frac{a_{\rm LL}}{a_0}\right)`. As an example, we
consider :math:`a_{\rm NNLL}` that, using
Eq. :eq:`eq:masterf`, reads:

.. math::
   \begin{array}{rcl}
     a_{\rm NNLL}&=&\displaystyle a_{\rm LL}\left\{1+a_{\rm
                     LL}\left[\overline{b}_1\ln\left(\frac{a_{\rm NLL}}{a_0}\right)+\overline{b}_{2}\left(a_{\rm
                     LL}-a_0\right)\right]\right\}^{-1}+\mathcal{O}(a^4)\\
     \\
                 &=&\displaystyle a_{\rm LL}\left\{1-a_{\rm LL}\left[\overline{b}_1 \ln\left(\frac{a_{\rm
                     LL}}{a_0}\right)+\overline{b}_{2}\left(a_{\rm
                     LL}-a_0\right) \right]+a_{\rm
                     LL}^2 \overline{b}_1^2\ln\left(\frac{a_{\rm LL}}{a_0}\right)\left(1+\ln\left(\frac{a_{\rm LL}}{a_0}\right)\right)\right\}+\mathcal{O}(a^4)\,.
   \end{array}

.. container::
   :name: analytic-perturbation-theory

   .. rubric:: Analytic perturbation theory
      :name: analytic-perturbation-theory

The fact highlighted above that the solution to all other RGEs relevant
to resum large logarithms of :math:`q_T` can be written in terms of
:math:`a_{\rm LL}` only is particularly relevant. As a matter of fact,
this will eventually allow us to eliminate the Landau pole simply by
requiring that :math:`a_{\rm LL}`, or equivalently the one-loop
:math:`\beta`-function, has the correct analytic structure in the
complex :math:`\mu` plane. This approach is often referred to analytic
perturbation theory (APT) and was originally suggested by Bogolyubov,
Logiunov, and Shirkov back in 1959 (Bogolyubov, Logunov, and Shirkov
1960) to obtain the correct causality properties of the QED Green’s
functions also when working at fixed-order in perturbation theory. At
one-loop order, this boils down to requiring that the strong coupling
:math:`a_{\rm LL}` admits a spectral representation allowing to restore
the analytic properties of the two-point Green’s function required for
causality to be fulfilled. The direct consequence of the procedure is
the natural appearance of power suppressed terms of the form
:math:`\Lambda/\mu` that automatically eliminate the Landau pole making
the strong coupling :math:`a_{\rm LL}` finite for all values of
:math:`\mu`.

However, :math:`a_{\rm LL}` does not only appears as integer powers but
also within logs. It is therefore crucial to be able to treat these
terms correctly preserving the correct analytical structure on the
complex :math:`\mu`. Here fractional analytical perturbation theory
(FAPT) comes to rescue (Bakulev, Mikhailov, and Stefanis 2005).

.. container::
   :name: derivation-of-analytic-expression-for-g-functions-in-q_t-resummation

   .. rubric:: Derivation of analytic expression for g-functions in
      :math:`q_{T}`-resummation
      :name: derivation-of-analytic-expression-for-g-functions-in-q_t-resummation

.. container::
   :name: the-sudakov-integral

   .. rubric:: The Sudakov integral
      :name: the-sudakov-integral

The Sudakov exponential has the following explicit form:

.. math:: S=-\int_{b_{0}^{2}/{b^{2}}}^{Q^{2}} \frac{dq^{2}}{q^{2}}\left[A\left(\alpha_{s}\left(q^{2}\right)\right)\ln\frac{M^{2}}{q^{2}}+\tilde B \left(\alpha_{s}\left(q^{2}\right)\right)\right],

where :math:`M` is the invariant mass of the final state, :math:`Q` is
the resummation scale, and :math:`A, \tilde B` are power series in
:math:`\alpha_{s}\left(Q^{2}\right)`.

The :math:`A` function is completely determined by the cusp anomalous
dimensions.

The :math:`\tilde B` function contains the usual :math:`B` coefficients
and two additional terms coming from PDF evolution and coefficient
functions.

From the PDF evolution equations in Mellin space

.. math::
   \begin{aligned}
   f_{a,N}\left(\mu_{F}^{2}\right)&=&\sum _{b} U_{ab,N}\left(\frac{b_{0}^{2}}{b^{2}},\mu _{F}^{2}\right) f_{b,N}\mu _{F}^{2}\\
   \frac{dU_{ab,N}\left(\mu ^{2},\mu _{F}^{2}\right)}{d\ln\mu^{2}}&=&\sum _{c}\gamma _{ac,N}\left(\alpha _{s}\left(\mu^{2}\right)\right) U_{cb,N}\left(\mu ^{2},\mu _{F}^{2}\right)\end{aligned}

one obtains (for a single species of partons) the solution

.. math:: U_{N}\left(b_{0}^{2}/b^{2},\mu _{F}^{2}\right)=\exp\left\{-\int_{b_{0}^{2}/{b^{2}}}^{\mu_{F}^{2}} \frac{dq^{2}}{q^{2}} \gamma_{N}\left(\alpha_{s}\left(q^{2}\right)\right)\right\}

Thus we have the first additional contribution to :math:`B`
coefficients:
:math:`B\longrightarrow B + 2 \, \gamma_{N} \, \alpha_{s}`.

As for the coefficient functions (again, in Mellin space),we can use the
following renormalisation-group identity:

.. math:: C_{N}\left(\alpha_{s}\left(b_{0}^{2}/b^{2}\right)\right)=C_{N}\left(\alpha_{s}\left(M^{2}\right)\right)\exp\left\{-\int_{b_{0}^{2}/{b^{2}}}^{M^{2}} \frac{dq^{2}}{q^{2}} \beta\left(\alpha_{s}\left(q^{2}\right)\right)\frac{d\ln C_{N}\left(\alpha_{s}\left(q^{2}\right)\right)}{d\ln\alpha_{s}\left(q^{2}\right)}\right\}

to derive the second additional contribution to :math:`B` coefficients:
:math:`B\longrightarrow B + 2 \, \beta\left(\alpha_{s}\right) \frac{d\ln C_{N}\left(\alpha_{s}\right)}{d\ln\alpha_{s}}`.

The :math:`\tilde B` function is thus given by

.. math:: \tilde B = B + 2 \, \gamma_{N} \, \alpha_{s} + 2 \, \beta\left(\alpha_{s}\right) \frac{d\ln C_{N}\left(\alpha_{s}\right)}{d\ln\alpha_{s}}.

With simple math and a change of the integration variable, we have

.. math::
   \begin{aligned}
   S&=&-\int_{b_{0}^{2}/{b^{2}}}^{Q^{2}} \frac{dq^{2}}{q^{2}}\left[A\left(\alpha_{s}\left(q^{2}\right)\right)\left(\ln\frac{M^{2}}{Q^{2}}+\ln\frac{Q^{2}}{q^{2}}\right)+\tilde B \left(\alpha_{s}\left(q^{2}\right)\right)\right]\\
   &=&-\int_{b_{0}^{2}/{b^{2}}}^{Q^{2}} \frac{dq^{2}}{q^{2}}\left[A\left(\alpha_{s}\left(q^{2}\right)\right)\ln\frac{Q^{2}}{q^{2}}+\bar B \left(\alpha_{s}\left(q^{2}\right)\right)\right]\\
   &=&\int_{\lambda}^{0} \frac{dt}{\beta_{0}\alpha_{s}}\left[A\left(\alpha_{s}\left(t\right)\right)\frac{t}{\beta_{0}\alpha_{s}}+\bar B \left(\alpha_{s}\left(t\right)\right)\right],\end{aligned}

where
:math:`\lambda=\beta_{0}\alpha_{s}\ln\left(Q^{2}b^{2}/b_{0}^{2}\right)=\beta_{0}\alpha_{s}L`,
and the :math:`A,\bar B` functions are given by

.. math::
   \begin{aligned}
   A&=&\sum _{n=1}^{\infty}A_{n}\alpha _{s}^{n},\\
   \bar{B}&=&\sum _{n=1}^{\infty}\left(\tilde B_{n}+A_{n}\ln\frac{M^{2}}{Q^{2}}\right)\alpha _{s}^{n} \end{aligned}

i.e., the :math:`\ln\left(M^{2}/Q^{2}\right)` dependence is entirely
embodied in the shift :math:`\tilde B\longrightarrow \bar B`.

We will now explicitly compute the :math:`g_{n}` functions, including
the :math:`\mu_{R}`-dependent terms coming from the replacement
:math:`Q\longrightarrow\mu_{R}`.

.. container::
   :name: alpha_s-evolution

   .. rubric:: :math:`\alpha_{s}` evolution
      :name: alpha_s-evolution

The evolution of :math:`\alpha_{s}` from :math:`Q^{2}` to :math:`q^{2}`
up to NNLL order is given by:

.. math:: \alpha_{s}\left(q^{2}\right)=\frac{\alpha_{s}}{1-t}-\left(\frac{\alpha_{s}}{1-t}\right)^{2}\frac{\beta_{1}}{\beta_{0}}\ln\left(1-t\right)+\left(\frac{\alpha_{s}}{1-t}\right)^{3}\left(\frac{\beta_{1}^{2}}{\beta_{0}^{2}}\left(\ln^{2}\left(1-t\right)-\ln\left(1-t\right)-t\right)+\frac{\beta_{2}}{\beta_{0}}t\right)+{\cal O}\left(\alpha_{s}^{4}\right)

with :math:`\alpha_{s}=\alpha_{s}\left(Q^{2}\right)` and
:math:`t=\beta_{0}\alpha_{s}\ln\frac{Q^{2}}{q^{2}}`.

From now on, I will call the first three terms of this expansion
LL\ :math:`\alpha_{s}`, NLL\ :math:`\alpha_{s}` and
NNLL\ :math:`\alpha_{s}`, respectively.

For future reference, using the expansions
:math:`\frac{1}{1-t}=1+t+t^{2}+\dots` and
:math:`\ln\left(1-t\right)=-t-\frac{1}{2}t^{2}+\dots`, we get the
following expansion for :math:`\alpha_{s}`:

.. math::
   :label: alphasPT
   \alpha_{s}\left(q^{2}\right)=\alpha_{s}-\alpha_{s}^{2}\beta_{0}\ln\frac{q^{2}}{Q^{2}}-\alpha_{s}^{3}\left(\beta_{1}\ln\frac{q^{2}}{Q^{2}}-\beta_{0}^{2}\ln^{2}\frac{q^{2}}{Q^{2}}\right)

.. container::
   :name: g_1

   .. rubric:: :math:`g_{1}`
      :name: g_1

The LL function is easily obtained by integrating the :math:`A_{1}`
term:

.. math:: \int _{\lambda }^{0}\frac{t\text{dt}}{\beta _{0}^{2}\alpha _{s}^{2}}\left[A_{1}\alpha _{s}^{\text{LL}}\left(t\right)\right]

| The only possibility here is to use the LL expression for
  :math:`\alpha_{s}`, since any other choice would be subleading. The
  integrand is of :math:`{\cal O}(\alpha_{s}^{-1})`.
| The result (multiplied by :math:`\beta_{0}\alpha_{s}/\lambda` in order
  to write the integral as L\ :math:`g_{1}`) is:

  .. math:: g_1(\lambda)=\frac{A_1 (\lambda +\ln(1-\lambda))}{\beta _{0}\lambda }

.. container::
   :name: g_2

   .. rubric:: :math:`g_{2}`
      :name: g_2

The NLL function has three contributions:

.. math:: \int _{\lambda }^{0}\frac{t\text{dt}}{\beta _{0}^{2}\alpha _{s}^{2}}\left[A_{1}\alpha _{s}^{\text{NLL}}\left(t\right)+A_{2}\left(\alpha _{s}^{\text{LL}}\left(t\right)\right)^{2}\right]+\int _{\lambda }^{0}\frac{\text{dt}}{\beta _{0}\alpha _{s}}\left[{\bar B}_{1}\alpha _{s}^{\text{LL}}\left(t\right)\right]

| The combinations are chosen in such a way that all the integrands are
  :math:`{\cal O}(\alpha_{s}^{0})`.
| We obtain

.. math:: g_2(\lambda)=\frac{\bar{B}_1}{\beta _0}\ln(1-\lambda )+\frac{A_1 \beta _1}{\beta _0^3} \left(\frac{\lambda }{1-\lambda }+\frac{1}{2} \log ^2(1-\lambda )+\frac{\log (1-\lambda )}{1-\lambda }\right)-\frac{A_2}{\beta _0^2}\left(\frac{\lambda }{1-\lambda }+\log (1-\lambda )\right).

.. container::
   :name: g_3

   .. rubric:: :math:`g_{3}`
      :name: g_3

The NNLL function has five contributions:

.. math:: \int _{\lambda }^{0}\frac{t\text{dt}}{\beta _{0}^{2}\alpha _{s}^{2}}\left[A_{1}\alpha _{s}^{\text{NNLL}}\left(t\right)+A_{2}\left(2\alpha _{s}^{\text{LL}}\left(t\right)\alpha_{s}^{\text{NLL}}\left(t\right)\right)+A_{3}\left(\alpha _{s}^{\text{LL}}\left(t\right)\right)^{3}\right]+\int _{\lambda }^{0}\frac{\text{dt}}{\beta _{0}\alpha _{s}}\left[{\bar B}_{1}\alpha _{s}^{\text{NLL}}\left(t\right)+{\bar B}_{2}\left(\alpha _{s}^{\text{LL}}\left(t\right)\right)^{2}\right]

| The combinations are chosen in such a way that all the integrands are
  :math:`{\cal O}(\alpha_{s})`.
| Summing up the five terms, we obtain

.. math::
   \begin{aligned}
   g_3(\lambda)&=&-\frac{\bar{B}_2}{\beta _0}\frac{\lambda}{(1-\lambda)}+\frac{ \bar{B}_1\beta _1}{\beta _0^2 }\frac{\lambda+\ln(1-\lambda)}{1-\lambda} -\frac{A_3}{2 \beta _0^2}\frac{\lambda ^2}{(1-\lambda)^2}\nonumber\\&+&A_1 \left(\frac{\lambda\left(\beta _1^2 \lambda +\beta _0 \beta _2 (2-3 \lambda)\right)}{2 \beta _0^4 (1-\lambda)^2}+\frac{\beta _1^2 (1-2 \lambda) \ln^2(1-\lambda)}{\left(2 \beta _0^4\right) (1-\lambda)^2}+\left(\frac{\beta _1^2}{\beta _0^4 (1-\lambda)}+\frac{\beta _0 \beta _2-\beta _1^2}{\beta _0^4}\right) \ln(1-\lambda)\right)\nonumber\\&+&\frac{A_2 \beta _1}{\beta _0^3}\left(\frac{\lambda  (3 \lambda -2)}{2 (1-\lambda)^2}-\frac{(1-2 \lambda) \ln(1-\lambda)}{(1-\lambda )^2}\right).\end{aligned}

.. container::
   :name: introduction-of-mu_r-dependence

   .. rubric:: Introduction of :math:`\mu_{R}`-dependence
      :name: introduction-of-mu_r-dependence

We now allow the resummation scale :math:`Q` to be different from the
renormalisation scale :math:`\mu_{R}`, and consequently compute the
variation of :math:`g_{1}` with respect to this change of scale.

Considering the expansion of :math:`\alpha_{s}` up to
:math:`{\cal O}(\alpha_{s}^{2})` in Eq.\ `[alphasPT] <#alphasPT>`__, we
derive the following formulae for :math:`g_{1}` variation:

.. math::
   \begin{aligned}
   \lambda&\longrightarrow&\lambda-\beta_{0}^{2}\alpha_{s}^{2}(\mu_{R}^{2})L\ln\left(\frac{Q^{2}}{\mu_{R}^{2}}\right)\\
   Lg_{1}&\longrightarrow&Lg_{1}+Lg_{1}^{'}(\lambda)\Delta\lambda=Lg_{1}-\lambda^{2}g_{1}^{'}(\lambda)\ln\left(\frac{Q^{2}}{\mu_{R}^{2}}\right)\end{aligned}

The result is a next-to-leading order contribution that has to be added
to :math:`g_{2}`:

.. math:: \Delta g_{1}^{LL}(\lambda)=\frac{A_1}{\beta _0}\left(\frac{\lambda}{1-\lambda}+\ln(1-\lambda)\right)\ln\left(\frac{Q^2}{\mu _R^2}\right)

Considering the expansion of :math:`\alpha_{s}` up to
:math:`{\cal O}(\alpha_{s}^{3})` in Eq.\ `[alphasPT] <#alphasPT>`__, we
get two shifts:

.. math::
   \begin{aligned}
   \lambda&\longrightarrow&\lambda-\beta_{0}\beta_{1}\alpha_{s}^{3}(\mu_{R}^{2})L\ln\left(\frac{Q^{2}}{\mu_{R}^{2}}\right)\\
   Lg_{1}&\longrightarrow&Lg_{1}+Lg_{1}^{'}(\lambda)\Delta\lambda=Lg_{1}-\alpha_{s}\frac{\beta_{1}}{\beta_{0}}\lambda^{2}g_{1}^{'}(\lambda)\ln\left(\frac{Q^{2}}{\mu_{R}^{2}}\right)\end{aligned}

and

.. math::
   \begin{aligned}
   \lambda&\longrightarrow&\lambda+\beta_{0}^{3}\alpha_{s}^{3}(\mu_{R}^{2})L\ln\left(\frac{Q^{2}}{\mu_{R}^{2}}\right)\\
   Lg_{1}&\longrightarrow&Lg_{1}+Lg_{1}^{'}(\lambda)\Delta\lambda=Lg_{1}+\beta_{0}\alpha_{s}\lambda^{2}g_{1}^{'}(\lambda)\ln\left(\frac{Q^{2}}{\mu_{R}^{2}}\right)\end{aligned}

These are both next-to-next-to-leading order contributions and will be
added to :math:`g_{3}`:

.. math:: \Delta g_{1}^{NLL1}(\lambda)=\alpha_s\frac{A_1}{\beta _0^2}\left(\frac{\lambda}{1-\lambda}+\ln(1-\lambda)\right)\left(\beta_1-\beta _0^2\ln\left(\frac{Q^2}{\mu _R^2}\right)\right)\ln\left(\frac{Q^2}{\mu _R^2}\right)

A further term generating a contribution of the same order is the third
term of the Taylor expansion (computed with first order squared
variation of :math:`\alpha_{s}`:

.. math::
   \begin{aligned}
   \lambda&\longrightarrow&\lambda-\beta_{0}^{2}\alpha_{s}^{2}(\mu_{R}^{2})L\ln\left(\frac{Q^{2}}{\mu_{R}^{2}}\right)\\
   Lg_{1}&\longrightarrow&Lg_{1}+\frac{1}{2}g_{1}^{''}(\lambda)(\Delta\lambda)^{2}=Lg_{1}+\frac{1}{2}\beta_{0}\alpha_{s}\lambda^{3}g_{1}^{''}(\lambda)\ln^{2}(\frac{Q^{2}}{\mu_{R}^{2}})\end{aligned}

which gives

.. math:: \Delta g_{1}^{NLL2}(\lambda)=\alpha _s A_1\left(\frac{\lambda(2-3\lambda)}{2(\lambda-1)^2}+\ln(1-\lambda)\right)\ln ^2\left(\frac{Q^2}{\mu _R^2}\right)

Applying the same variation at :math:`g_{2}` we get

.. math::
   \begin{aligned}
   \lambda&\longrightarrow&\lambda-\beta_{0}^{2}\alpha_{s}^{2}(\mu_{R}^{2})L\ln\left(\frac{Q^{2}}{\mu_{R}^{2}}\right)\\
   g_{2}&\longrightarrow&g_{2}+g_{2}^{'}(\lambda)\Delta\lambda=g_{2}-\beta_{0}\alpha_{s}(\mu_{R}^{2})\lambda g_{2}^{'}(\lambda)\ln\left(\frac{Q^{2}}{\mu_{R}^{2}}\right)\end{aligned}

The resulting next-to-next-to-leading contribution to add to
:math:`g_{3}` is:

.. math:: \Delta g_{2}^{LL}(\lambda)=\alpha _s\left(\bar{B}_1\frac{\lambda}{1-\lambda }+\frac{A_2}{\beta _0} \frac{\lambda ^2}{(1-\lambda )^2}-\frac{A_1\beta_1}{\beta _0^2}\frac{\lambda ^2}{(1-\lambda )^2}\ln(1-\lambda)\right)\ln\left(\frac{Q^2}{\mu _R^2}\right).

.. container::
   :name: final-expressions

   .. rubric:: Final expressions
      :name: final-expressions

The final expression for :math:`g`-functions including
renormalisation-scale dependence is:

.. math::
   \begin{aligned}
   g_1(\lambda)&=&\frac{A_1 (\lambda +\ln(1-\lambda))}{\beta _{0}\lambda}\\
   g_2(\lambda)&=&\frac{\bar{B}_1}{\beta _0}\ln(1-\lambda )+\frac{A_1 \beta _1}{\beta _0^3} \left(\frac{\lambda }{1-\lambda }+\frac{1}{2} \log ^2(1-\lambda )+\frac{\log (1-\lambda )}{1-\lambda }\right)-\frac{A_2}{\beta _0^2}\left(\frac{\lambda }{1-\lambda }+\log (1-\lambda )\right)\nonumber\\&+&\frac{A_1}{\beta _0}\left(\frac{\lambda}{1-\lambda}+\ln(1-\lambda)\right)\ln\left(\frac{Q^2}{\mu _R^2}\right)\\
   g_3(\lambda)&=&-\frac{\bar{B}_2}{\beta _0}\frac{\lambda}{(1-\lambda)}+\frac{ \bar{B}_1\beta _1}{\beta _0^2 }\frac{\lambda+\ln(1-\lambda)}{1-\lambda} -\frac{A_3}{2 \beta _0^2}\frac{\lambda ^2}{(1-\lambda)^2}\nonumber\\&+&A_1 \left(\frac{\lambda\left(\beta _1^2 \lambda +\beta _0 \beta _2 (2-3 \lambda)\right)}{2 \beta _0^4 (1-\lambda)^2}+\frac{\beta _1^2 (1-2 \lambda) \ln^2(1-\lambda)}{\left(2 \beta _0^4\right) (1-\lambda)^2}+\left(\frac{\beta _1^2}{\beta _0^4 (1-\lambda)}+\frac{\beta _0 \beta _2-\beta _1^2}{\beta _0^4}\right) \ln(1-\lambda)\right)\nonumber\\&+&\frac{A_2 \beta _1}{\beta _0^3}\left(\frac{\lambda  (3 \lambda -2)}{2 (1-\lambda)^2}-\frac{(1-2 \lambda) \ln(1-\lambda)}{(1-\lambda )^2}\right)-\frac{A_{1}}{2}\frac{\lambda^{2}}{(1-\lambda)^{2}}\ln^{2}\left(\frac{Q^2}{\mu _R^2}\right)\nonumber\\&+&\left(\bar{B}_1\frac{\lambda}{1-\lambda }+\frac{A_2}{\beta _0}\frac{\lambda ^2}{(1-\lambda )^2}+\frac{A_1\beta_1}{\beta _0^2}\left(\frac{\lambda}{1-\lambda}+\frac{(1-2\lambda)\ln(1-\lambda)}{(1-\lambda)^{2}}\right)\right)\ln\left(\frac{Q^2}{\mu _R^2}\right)\end{aligned}

**References**

**References**

.. container:: references csl-bib-body hanging-indent
   :name: refs

   .. container:: csl-entry
      :name: ref-Bacchetta:2008xw

      Bacchetta, Alessandro, Daniel Boer, Markus Diehl, and Piet J.
      Mulders. 2008. “Matches and mismatches in the descriptions of
      semi-inclusive processes at low and high transverse momentum.”
      *JHEP* 08: 023. https://doi.org/10.1088/1126-6708/2008/08/023.

   .. container:: csl-entry
      :name: ref-Bakulev:2005gw

      Bakulev, A. P., S. V. Mikhailov, and N. G. Stefanis. 2005. “QCD
      analytic perturbation theory: From integer powers to any power of
      the running coupling.” *Phys. Rev. D* 72: 074014.
      https://doi.org/10.1103/PhysRevD.72.074014.

   .. container:: csl-entry
      :name: ref-Bertone:2017gds

      Bertone, Valerio. 2018. “APFEL++: A new PDF evolution library in
      C++.” Edited by Uta Klein. *PoS* DIS2017: 201.
      https://doi.org/10.22323/1.297.0201.

   .. container:: csl-entry
      :name: ref-Bogolyubov:1960cha

      Bogolyubov, N. N., A. A. Logunov, and D. V. Shirkov. 1960. “The
      method of dispersion relations and perturbation theory.” *Sov.
      Phys. JETP* 10 (3): 574–81.

   .. container:: csl-entry
      :name: ref-Bozzi:2005wk

      Bozzi, Giuseppe, Stefano Catani, Daniel de Florian, and
      Massimiliano Grazzini. 2006. “Transverse-momentum resummation and
      the spectrum of the Higgs boson at the LHC.” *Nucl. Phys. B* 737:
      73–120. https://doi.org/10.1016/j.nuclphysb.2005.12.022.

   .. container:: csl-entry
      :name: ref-Catani:2012qa

      Catani, Stefano, Leandro Cieri, Daniel de Florian, Giancarlo
      Ferrera, and Massimiliano Grazzini. 2012. “Vector boson production
      at hadron colliders: hard-collinear coefficients at the NNLO.”
      *Eur. Phys. J. C* 72: 2195.
      https://doi.org/10.1140/epjc/s10052-012-2195-7.

   .. container:: csl-entry
      :name: ref-Collins:2016hqq

      Collins, J., L. Gamberg, A. Prokudin, T. C. Rogers, N. Sato, and
      B. Wang. 2016. “Relating Transverse Momentum Dependent and
      Collinear Factorization Theorems in a Generalized Formalism.”
      *Phys. Rev. D* 94 (3): 034014.
      https://doi.org/10.1103/PhysRevD.94.034014.

   .. container:: csl-entry
      :name: ref-Collins:2017oxh

      Collins, John, and Ted C. Rogers. 2017. “Connecting Different TMD
      Factorization Formalisms in QCD.” *Phys. Rev. D* 96 (5): 054011.
      https://doi.org/10.1103/PhysRevD.96.054011.

   .. container:: csl-entry
      :name: ref-Daleo:2004pn

      Daleo, A., D. de Florian, and R. Sassot. 2005. “O(alpha**2(s)) QCD
      corrections to the electroproduction of hadrons with high
      transverse momentum.” *Phys. Rev. D* 71: 034013.
      https://doi.org/10.1103/PhysRevD.71.034013.

   .. container:: csl-entry
      :name: ref-DelDebbio:2007ee

      Del Debbio, Luigi, Stefano Forte, Jose I. Latorre, Andrea
      Piccione, and Joan Rojo. 2007. “Neural network determination of
      parton distributions: The Nonsinglet case.” *JHEP* 03: 039.
      https://doi.org/10.1088/1126-6708/2007/03/039.

   .. container:: csl-entry
      :name: ref-Echevarria:2016scs

      Echevarria, Miguel G., Ignazio Scimemi, and Alexey Vladimirov.
      2016. “Unpolarized Transverse Momentum Dependent Parton
      Distribution and Fragmentation Functions at
      next-to-next-to-leading order.” *JHEP* 09: 004.
      https://doi.org/10.1007/JHEP09(2016)004.

   .. container:: csl-entry
      :name: ref-Ellis:1991qj

      Ellis, R.Keith, W.James Stirling, and B. R. Webber. 2011. *QCD and
      collider physics*. Vol. 8. Cambridge University Press.

   .. container:: csl-entry
      :name: ref-deFlorian:1997zj

      Florian, D. de, M. Stratmann, and W. Vogelsang. 1998. “QCD
      analysis of unpolarized and polarized Lambda baryon production in
      leading and next-to-leading order.” *Phys. Rev. D* 57: 5811–24.
      https://doi.org/10.1103/PhysRevD.57.5811.

   .. container:: csl-entry
      :name: ref-Furmanski:1981cw

      Furmanski, W., and R. Petronzio. 1982. “Lepton - Hadron Processes
      Beyond Leading Order in Quantum Chromodynamics.” *Z. Phys. C* 11:
      293. https://doi.org/10.1007/BF01578280.

   .. container:: csl-entry
      :name: ref-Meng:1995yn

      Meng, Ruibin, Fredrick I. Olness, and Davison E. Soper. 1996.
      “Semiinclusive deeply inelastic scattering at small q(T).” *Phys.
      Rev. D* 54: 1919–35. https://doi.org/10.1103/PhysRevD.54.1919.

   .. container:: csl-entry
      :name: ref-Nadolsky:1999kb

      Nadolsky, Pavel M., D. R. Stump, and C. P. Yuan. 2000.
      “Semiinclusive hadron production at HERA: The Effect of QCD gluon
      resummation.” *Phys. Rev. D* 61: 014003.
      https://doi.org/10.1103/PhysRevD.64.059903.

.. [1]
   Notice that, despite the variables :math:`x` and :math:`\mathbf{b}_T`
   will not appear explicitly, the symbol :math:`\otimes` indicates the
   Mellin convolution integral w.r.t. :math:`x` while :math:`b_T`
   indicates the length of the vector :math:`\mathbf{b}_T`.

.. [2]
   A sum over flavours is understood. As a matter of fact, the matching
   function :math:`C` has to be regarded as a matrix in flavour space
   multipling a vector of collinear PDFs/FFs.

.. [3]
   A summation over the flavour indices is understood

.. [4]
   The convergence of the Ogata quadrature method crucially depends on
   the function :math:`f_{\rm NP}` being exponentially suppressed for
   large values of :math:`b_T`, such that full integrand is dumped in
   this region.

.. [5]
   The factor :math:`z^{2}` is the consequence of the fact that we are
   writing the cross section differential in :math:`q_T^2` that is the
   transverse momentum of the exchanged photon while in Ref. (Bacchetta
   et al. 2008) the cross section is differentila in :math:`p_T^2` that
   is the transverse momentum of the of the outgoing hadrons. Since
   :math:`p_T=z q_T`, the factor :math:`z^2` cancels.

.. [6]
   The proof of this relation is given in
   Appendix `1 <#sec:deltaexpansion>`__.

.. [7]
   Notice that the :math:`z` variable of Ref. (Daleo, Florian, and
   Sassot 2005) does not coincide with our definition.

.. [8]
   This step has the effect of removing the uncanceled soft
   divergencies.

.. [9]
   The expressions reported in Ref. (Florian, Stratmann, and Vogelsang
   1998) are those computed long time ago in Ref. (Furmanski and
   Petronzio 1982).

.. [10]
   Notice that :math:`\overline{\beta}` in
   Eq. :eq:`eq:RGEalphas` is related to the usual QCD
   :math:`\beta`-function by :math:`\overline{\beta}=-\beta/a_s`.
