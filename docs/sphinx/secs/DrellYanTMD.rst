============================================
Drell-Yan cross section in TMD factorisation
============================================

:Author: Valerio Bertone

.. role:: raw-latex(raw)
   :format: latex
..

.. contents::
   :depth: 3
..

Structure of the observables
============================

Let us start from Eq. (2.6) of Ref. (Scimemi and Vladimirov 2018), that
is the fully differential cross section for lepton-pair production in
the region in which the TMD factorisation applies, :math:`i.e.`
:math:`q_T \ll
Q`. After some minor manipulations, it reads:

.. math::
   :label: eq:crosssection
     \frac{d\sigma}{dQ dy dq_T} =
     \frac{16\pi\alpha^2q_T}{9 Q^3} H(Q,\mu) \sum_q C_q(Q)
     \int\frac{d^2\mathbf{b}}{4\pi} e^{i \mathbf{b}\cdot \mathbf{q}_T} \overline{F}_q(x_1,\mathbf{b};\mu,\zeta) \overline{F}_{\bar{q}}(x_2,\mathbf{b};\mu,\zeta)\,,

where :math:`Q`, :math:`y`, and :math:`q_T` are the invariant mass, the
rapidity, and the transverse momentum of the lepton pair, respectively,
while :math:`\alpha` is the electromagnetic coupling, :math:`H` is the
appropriate QCD hard factor that can be perturbatively computed, and
:math:`C_q` are the effective electroweak charges. In addition, the
variables :math:`x_1` and :math:`x_2` are functions of :math:`Q` and
:math:`y` and are given by:

.. math::
   :label: eq:Bjorkenx12
     x_{1,2} = \frac{Q}{\sqrt{s}}e^{\pm y}\,,

being :math:`\sqrt{s}` the centre-of-mass energy of the collision. In
Eq. :eq:`eq:crosssection` we are using the
short-hand notation:

.. math:: \overline{F}_q(x,\mathbf{b};\mu,\zeta) \equiv xF_q(x,\mathbf{b};\mu,\zeta)\,,

that is convenient for the implementation. The scales :math:`\mu` and
:math:`\zeta` are introduced as a consequence of the removal of UV and
rapidity divergences in the definition of the TMDs. Despite these scales
are arbitrary scales, they are typically chosen
:math:`\mu=\sqrt{\zeta}=Q`. Therefore, for all practical purposes their
presence is fictitious.

The computation-intensive part of
Eq.:eq:`eq:crosssection` has the form of the
integral:

.. math::
   :label: eq:integral
   I_{ij}(x_1,x_2,q_T;\mu,\zeta)=\int\frac{d^2\mathbf{b}}{4\pi} e^{i \mathbf{b}\cdot \mathbf{q}_T} \overline{F}_i(x_1,\mathbf{b};\mu,\zeta) \overline{F}_{j}(x_2,\mathbf{b};\mu,\zeta)\,.

where :math:`\overline{F}_{i(j)}` are combinations of evolved TMD PDFs.
At this stage, for convenience, :math:`i` and :math:`j` do not coincide
with :math:`q` and :math:`\bar{q}` but they are linked through a simple
linear transformation. The integral over the bidimensional impact
parameter **b** has to be taken. However, :math:`\overline{F}_{i(j)}`
only depend on the absolute value of **b**, therefore
Eq. :eq:`eq:integral` can be written as:

.. math::
   :label: eq:integral2
   I_{ij}(x_1,x_2,q_T;\mu,\zeta)=\frac12\int_0^\infty db\,b J_0(bq_T)  \overline{F}_i(x_1,b;\mu,\zeta) \overline{F}_{j}(x_2,b;\mu,\zeta)\,.

where :math:`J_0` is the zero-th order Bessel function of the first kind
whose integral representation is:

.. math:: J_0(x) = \frac1{2\pi}\int_0^{2\pi} d\theta e^{ix\cos(\theta)}\,.

The evolved quark TMD PDF :math:`\overline{F}_i` at the final scales
:math:`\mu` and :math:`\zeta` is obtained by multiplying the same
distribution at the initial scales :math:`\mu_0` and :math:`\zeta_0` by
a single evolution factor :math:`R_q`\ ( [1]_). that is:

.. math::
   :label: eq:evolution
     \overline{F}_i(x,b;\mu,\zeta) = R_q(\mu_0,\zeta_0\rightarrow \mu,\zeta;b)
     \overline{F}_i(x,b;\mu_0,\zeta_0)\,.

The initial scale TMD PDFs at small values :math:`b` can be written as:

.. math::
   :label: eq:LOconv
   \overline{F}_i(x,b;\mu_0,\zeta_0) = \sum_{j=g,q(\bar{q})}x\int_x^1\frac{dy}{y}C_{ij}(y;\mu_0,\zeta_0)f_j\left(\frac{x}{y},\mu_0\right)\,,

where :math:`f_j` are the collinear PDFs (including the gluon) and
:math:`C_{ij}` are the so-called matching functions that are
perturbatively computable and are currently known to NNLO, :math:`i.e.`
:math:`\mathcal{O}(\alpha_s^2)`. If we define:

.. math:: \overline{f}_i\left(x,\mu_0\right) = xf_i\left(x,\mu_0\right)\,,

Eq. :eq:`eq:LOconv` can be written as:

.. math::
   :label: eq:LOconvMod
   \overline{F}_i(x,b;\mu_0,\zeta_0) =
   \sum_{j=g,q(\bar{q})}\int_x^1dy\,C_{ij}(y;\mu_0,\zeta_0)  \overline{f}_i\left(\frac{x}{y},\mu_0\right)\,.

Putting Eqs. :eq:`eq:evolution`
and :eq:`eq:LOconvMod`, one finds:

.. math::
   :label: eq:pertTMD
     \overline{F}_i(x,b;\mu,\zeta) = R_q(\mu_0,\zeta_0\rightarrow \mu,\zeta;b)
     \sum_{j=g,q(\bar{q})}\int_x^1dy\,C_{ij}(y;\mu_0,\zeta_0)  \overline{f}_i\left(\frac{x}{y},\mu_0\right)\,.

Matching and evolution are affected by non-perturbative effects that
become relevant at large :math:`b`. In order to account for such
effects, one usually introduces a phenomenological function
:math:`f_{\rm NP}`. In the traditional approach (CSS (Collins 2013)),
the :math:`b`-space TMDs get a multiplicative correction that does not
depend on the flavour. In addition, the perturbative content of the TMDs
is smoothly damped away at large :math:`b` by introducing the so-called
:math:`b_*`-prescription:

.. math::
   :label: eq:LOconvNP1
     \overline{F}_i(x,b;\mu,\zeta) \rightarrow \overline{F}_i(x,b_*(b);\mu,\zeta) f_{\rm NP}(x,b,\zeta)\,,

where :math:`b_*\equiv b_*(b)` is a monotonic function of the impact
parameter :math:`b` such that:

.. math::
   \lim_{b\rightarrow 0}
     b_*(b) = b_{\rm min}\quad\mbox{and}\quad\lim_{b\rightarrow \infty}
     b_*(b) = b_{\rm max}\,,

being :math:`b_{\rm min}` and :math:`b_{\rm max}` constant values both
in the perturbative region. Including the non-perturbative function,
Eq. :eq:`eq:integral2` becomes:

.. math::
   :label: eq:integral3
   \begin{array}{l}
   \displaystyle I_{ij}(x_1,x_2,q_T;\mu,\zeta) = \displaystyle \int_0^\infty db\,J_0(bq_T)\left[\frac{b}2 
   \overline{F}_i(x_1,b_*(b);\mu,\zeta) \overline{F}_{j}(x_2,b_*(b);\mu,\zeta) f_{\rm NP}(x_1,b,\zeta)
     f_{\rm NP}(x_2,b,\zeta) \right]\\
   \\
   \displaystyle =\frac{1}{q_T}\int_0^\infty d\bar{b}\,J_0(\bar{b})\left[\frac{\bar{b}}{2q_T} 
   \overline{F}_i(x_1,b_*\left(\frac{\bar{b}}{q_T}\right);\mu,\zeta) \overline{F}_{j}(x_2,,b_*\left(\frac{\bar{b}}{q_T}\right);\mu,\zeta) f_{\rm NP}\left(x_1,\frac{\bar{b}}{q_T},\zeta\right)
     f_{\rm NP}\left(x_2,\frac{\bar{b}}{q_T},\zeta\right) \right]
   \,.
   \end{array}

Eq. :eq:`eq:integral3` is a Hankel tranform and can be
efficiently computed using the so-called Ogata quadrature (Ogata, n.d.).
Effectively, the computation of the integral in
Eq. :eq:`eq:integral` is achieved through a weighted
sum:

.. math::
   :label: eq:ogataquadrature
   \begin{array}{rcl}
   I_{ij}(x_1,x_2,q_T;\mu,\zeta) &\simeq& \displaystyle
                                                        \frac1{q_T}\sum_{n=1}^N
                                                        \frac{w_n^{(0)}z_n^{(0)}}{2q_T} 
   \overline{F}_i\left(x_1,b_*\left
                                          (\frac{z_n^{(0)}}{q_T}\right);\mu,\zeta\right) \overline{F}_j\left(x_2,b_*\left (\frac{z_n^{(0)}}{q_T}\right);\mu,\zeta\right)\\
   \\
   &\times&\displaystyle f_{\rm NP}\left(x_1,\frac{z_n^{(0)}}{q_T},\zeta\right)
     f_{\rm NP}\left(x_2,\frac{z_n^{(0)}}{q_T},\zeta\right)\,,
   \end{array}

where the unscaled coordinates :math:`z_n^{(0)}` and the weights
:math:`w_n^{(0)}` can be precomputed in terms of the zero’s of the
Bessel function :math:`J_0` and one single parameter (see Ref. (Ogata,
n.d.) for more details, specifically Eqs. (5.1) and (5.2) or
Appendix `3 <#app:OgataQuadrature>`__ for the relevant formula to
compute the unscaled coordinates and the weights) [2]_. Based on the
(empirically verified) assumption that the absolute value of each term
in the sum in the r.h.s. of
Eq. :eq:`eq:ogataquadrature` is smaller than
that of the preceding one, the truncation number :math:`N` is chosen
dynamically in such a way that the :math:`(N+1)`-th term is smaller in
absolute value than a user-defined cutoff relatively to the sum of the
preceding :math:`N` terms.

Eq. :eq:`eq:ogataquadrature` factors out the
non-perturbative part of the calculation represented by
:math:`f_{\rm NP}` from the perturbative content. This is done on
purpose to devise a method in which the perturbative content is
precomputed and numerically convoluted with the non-perturbative
functions *a posteriori*. This is convenient in view of a fit of the
function :math:`f_{\rm NP}`.

As customary in QCD, the most convenient basis for the matching in
Eq. :eq:`eq:LOconv` is the so-called “evolution” basis
(*i.e.* :math:`\Sigma`, :math:`V`, :math:`T_3`, :math:`V_3`, etc.). In
fact, in this basis the operator matrix :math:`C_{ij}` is almost
diagonal with the only exception of crossing terms that couple the gluon
and the singlet :math:`\Sigma` distributions. As a consequence, this is
the most convenient basis for the computation of :math:`I_{ij}`. On the
other hand, TMDs in Eq. :eq:`eq:crosssection`
appear in the so-called “physical” basis (*i.e.* :math:`d`,
:math:`\bar{d}`, :math:`u`, :math:`\bar{u}`, etc.). Therefore, we need
to rotate :math:`F_{i(j)}` from the evolution basis, over which the
indices :math:`i` and :math:`j` run, to the physical basis. This is done
by means of an appropriate constant matrix :math:`T`, so that:

.. math::
   :label: eq:lumiInterRot
   \overline{F}_{q}(x_1,b;\mu,\zeta)= \sum_{i}T_{qi}F_{i}(x_1,b;\mu,\zeta)\,,

and similarly for :math:`\overline{F}_{\bar{q}}`. Putting all pieces
together, one can conveniently write the cross section in
Eq. :eq:`eq:crosssection` as:

.. math::
   :label: eq:fastdiffxsec
     \frac{d\sigma}{dQ dy dq_T} \simeq\sum_{n=1}^N w_n^{(0)} \frac{z_n^{(0)}}{q_T}S\left(x_1,x_2,\frac{z_n^{(0)}}{q_T};\mu,\zeta\right) f_{\rm NP}\left(x_1,\frac{z_n^{(0)}}{q_T},\zeta\right) f_{\rm NP}\left(x_2,\frac{z_n^{(0)}}{q_T},\zeta\right)\,,

with:

.. math::
   :label: eq:pertfact
   S(x_1,x_2,b;\mu,\zeta) =\frac{8\pi\alpha^2}{9 Q^3}
       H(Q,\mu) \sum_q C_q(Q) \left[
   \overline{F}_q\left(x_1,b_*(b);\mu,\zeta \right)\right] \left[
   \overline{F}_{\bar{q}}\left(x_2,b_*(b);\mu,\zeta \right)\right]\,.

Eq. :eq:`eq:fastdiffxsec` allows one to precompute
the weights :math:`S` in such a way that the differential cross section
in Eq. :eq:`eq:crosssection` can be computed as a
simple weighted sum of the non-perturbative contribution. A misleading
aspect of Eq. :eq:`eq:pertfact` is the fact that
:math:`S` has five arguments. In actual facts, :math:`S` only depends on
three independent variables. The reason is that :math:`\mu` and
:math:`\zeta` are usually taken to be proportional to :math:`Q` by a
constant factor. In addition :math:`x_1` and :math:`x_2` depend on
:math:`Q` and :math:`y` through
Eq. :eq:`eq:Bjorkenx12`. Therefore, the full
dependence on the kinematics of the final state of
Eq. :eq:`eq:crosssection` can be specified by
:math:`Q`, :math:`y` and :math:`q_T`.

Integrating over the final-state kinematic variables
====================================================

Despite Eq. :eq:`eq:fastdiffxsec` provides a
powerful tool for a fast computation of cross sections, it is often not
sufficient to allow for a direct comparison to experimental data. The
reason is that experimental measurements of differential distributions
are usually delivered as integrated over finite regions of the
final-state kinematic phase space. In other words, experiments measure
quantities like:

.. math::
   :label: eq:Intcrosssection
   \widetilde{\sigma}=\int_{Q_{\rm min}}^{Q_{\rm max}}dQ \int_{y_{\rm min}}^{y_{\rm max}}dy \int_{q_{T,\rm min}}^{q_{T,\rm max}}dq_T\left[\frac{d\sigma}{dQ dy dq_T} \right]\,.

As a consequence, in order to guarantee performance, we need to include
the integrations above in the precomputed factors.

.. container::
   :name: integrating-over-q_t

   .. rubric:: Integrating over :math:`q_T`
      :name: integrating-over-q_t

The integration over bins in :math:`q_T` can be carried out analytically
exploiting the following property of Bessel’s function:

.. math:: \frac{d}{dx}\left[x^mJ_m(x)\right]=x^mJ_{m-1}(x)\,,

\ that leads to:

.. math::
   :label: eq:besselproperty
   \int dx\,x J_0(x) = xJ_1(x)\quad\Rightarrow\quad \int_{x_1}^{x_2}
   dx\,x J_0(x) = x_2J_1(x_2) - x_1J_1(x_1)\,.

To see it, we observe that the differential cross section in
Eq. :eq:`eq:crosssection` has the following
structure:

.. math:: \frac{d\sigma}{dQ dy dq_T} \propto \int_0^\infty db\, q_T  J_0(bq_T)\dots

where the ellipses indicate terms that do not depend on :math:`q_T`.
Therefore, using Eq. :eq:`eq:besselproperty` we
find:

.. math::
   \begin{array}{l}
   \displaystyle \int_{q_{T,\rm min}}^{q_{T,\rm
     max}}dq_T\left[\frac{d\sigma}{dQ dy dq_T} \right] \propto \int_0^\infty db\,
     \int_{q_{T,\rm min}}^{q_{T,\rm
     max}} dq_T\,   q_T J_0(bq_T)\dots= \\
   \\
   \displaystyle \int_0^\infty \frac{db}{b^2}\,
     \int_{bq_{T,\rm min}}^{bq_{T,\rm
     max}} dx\,   x J_0(x)\dots=\int_0^\infty \frac{db}{b}\left[q_{T,\rm
     max}J_1(bq_{T,\rm max}) - q_{T,\rm
     min}J_1(bq_{T,\rm min})\right]\dots\,.
   \end{array}

Therefore, defining:

.. math:: K(q_T) \equiv \int dq_T\left[\frac{d\sigma}{dQ dy dq_T} \right]

as the indefinite integral over :math:`q_T` of the cross section in
Eq. :eq:`eq:crosssection`, we have that:

.. math::
   :label: eq:primitive
   \int_{q_{T,\rm min}}^{q_{T,\rm
     max}}dq_T\left[\frac{d\sigma}{dQ dy dq_T} \right] = K(Q,y,q_{T,\rm max})
   - K(Q,y,q_{T,\rm min})\,,

with:

.. math::
   :label: eq:Kexplicit
   \begin{array}{l}
    \displaystyle K(Q,y,q_T) =
     \frac{8\pi\alpha^2q_T}{9 Q^3} H(Q,\mu) \\
   \\
   \displaystyle \times
     \int_0^\infty db\, J_1(bq_T) \sum_q C_q(Q)\overline{F}_q(x_1,b;\mu,\zeta) \overline{F}_{\bar{q}}(x_2,b;\mu,\zeta) f_{\rm NP}(x_1,b,\zeta)
     f_{\rm NP}(x_2,b,\zeta)\,,
   \end{array}

that can be computed using the Ogata quadrature as:

.. math::
   :label: eq:OgataPrimitive
     K(Q,y,q_T) \simeq \sum_{n=1}^N w_n^{(1)} S\left(x_1,x_2,\frac{z_n^{(1)}}{q_T};\mu,\zeta\right) f_{\rm NP}\left(x_1,\frac{z_n^{(1)}}{q_T},\zeta\right) f_{\rm NP}\left(x_2,\frac{z_n^{(1)}}{q_T},\zeta\right)\,,

with :math:`S` defined in Eq. :eq:`eq:pertfact`. The
unscaled coordinates :math:`z_n^{(1)}` and the weights :math:`w_n^{(1)}`
can again be precomputed and stored in terms of the zero’s of the Bessel
function :math:`J_1`. Eq. :eq:`eq:primitive` reduces
the integration in :math:`q_T` to a calculation completely analogous to
the unintegrated cross section. This is particularly convenient because
it avoids the computation a numerical integration.

.. container::
   :name: sec:kincuts

   .. rubric:: Kinematic cuts
      :name: sec:kincuts

In the presence of kinematic cuts, such as those on the final-state
leptons in Drell-Yan, the analytic integration over :math:`q_T`
discussed above cannot be performed. The reason is that the
implementation of these cuts effectively introduces a
:math:`q_T`-dependent function :math:`\mathcal{P}`\ ( [3]_) in the
integral:

.. math:: \frac{d\sigma}{dQ dy dq_T} \propto \int_0^\infty db\, q_T  J_0(bq_T)\mathcal{P}(q_T)\dots\,,

that prevents the direct use of
Eq. :eq:`eq:besselproperty`. Since
:math:`\mathcal{P}` is a slowly-varying function of :math:`q_T` over the
typical bin size, we can approximate the integral over the bins in
:math:`q_T` as:

.. math::
   :label: eq:intbyparts
   \begin{array}{rcl}
   \displaystyle  \int_{q_{T,\rm min}}^{q_{T,\rm
         max}} dq_T\, q_T J_0(bq_T) \mathcal{P}(q_T)&\simeq&\displaystyle
     \mathcal{P}\left(\frac{q_{T,\rm max}+q_{T,\rm min}}2\right)\int_{q_{T,\rm min}}^{q_{T,\rm
         max}} dq_T\, q_T J_0(bq_T) \\
   \\
   &=& \displaystyle\mathcal{P}\left(\frac{q_{T,\rm max}+q_{T,\rm
       min}}2\right) \frac{1}{b}\left[q_{T,\rm max} J_1(bq_{T,\rm max}) - q_{T,\rm min} J_1(bq_{T,\rm min}) \right]\,.
   \end{array}

Unfortunately, this structure is inconvenient because it mixes different
bin bounds and prevents a recursive computation. However, we can try to
go further and, assuming that the bin width is small enough, we can
expand :math:`\mathcal{P}` is the following ways:

.. math::
   \begin{array}{l}
   \displaystyle\mathcal{P}\left(\frac{q_{T,\rm max}+q_{T,\rm
       min}}2\right)= \mathcal{P}\left(q_{T,\rm min}+\Delta q_T\right) =
     \mathcal{P}\left(q_{T,\rm min}\right) + \mathcal{P}'\left(q_{T,\rm
     min}\right)\Delta q_T +\mathcal{O}\left(\Delta q_T^2\right)\,,\\
   \\
   \displaystyle\mathcal{P}\left(\frac{q_{T,\rm max}+q_{T,\rm
       min}}2\right)= \mathcal{P}\left(q_{T,\rm max}-\Delta q_T\right) =
     \mathcal{P}\left(q_{T,\rm max}\right) - \mathcal{P}'\left(q_{T,\rm
     max}\right)\Delta q_T +\mathcal{O}\left(\Delta q_T^2\right)\,,
   \end{array}

with:

.. math::
   :label: eq:halfqTinterval
   \Delta q_T = \frac{q_{T,\rm max}- q_{T,\rm min}}2\,.

Therefore:

.. math::
   :label: eq:lastexpP
   \begin{array}{rcl}
     \displaystyle  b\int_{q_{T,\rm min}}^{q_{T,\rm
     max}} dq_T\, q_T J_0(bq_T) \mathcal{P}(q_T)&\simeq&\displaystyle
                                                         q_{T,\rm
                                                         max}
                                                         J_1(bq_{T,\rm
                                                         max})\left[
                                                         \mathcal{P}\left(q_{T,\rm
                                                         max}\right)
                                                         - 
                                                         \mathcal{P}'\left(q_{T,\rm
                                                         max}\right)\Delta
                                                         q_T\right]\\
   \\
   &-&\displaystyle q_{T,\rm
                                                         min}
                                                         J_1(bq_{T,\rm
                                                         min})\left[
                                                         \mathcal{P}\left(q_{T,\rm
                                                         min}\right)
                                                         +
                                                         \mathcal{P}'\left(q_{T,\rm
                                                         min}\right)\Delta
                                                         q_T\right]\,.
   \end{array}

The advantage of this formula as compared to
Eq. :eq:`eq:intbyparts` is that each single term
depends on one single bin-bound in :math:`q_T` rather than on a
combination of two consecutive bounds. Therefore, in the presence of
kinematic cuts, the actual form of the primitive function :math:`K`
defined in Eq. :eq:`eq:primitive` and given explicitly
in Eq. :eq:`eq:Kexplicit` is:

.. math::
   :label: eq:KexplicitCuts
   \begin{array}{l}
    \displaystyle K(Q,y,q_T) =
     \frac{8\pi\alpha^2q_T}{9 Q^3} H(Q,\mu)
     \left[\mathcal{P}\left(Q,y,q_{T}\right) \pm
     \mathcal{P}'\left(Q,y,q_{T}\right)\Delta q_T\right]\\
   \\
   \displaystyle \times
     \int_0^\infty db\, J_1(bq_T) \sum_q C_q(Q)\overline{F}_q(x_1,b;\mu,\zeta) \overline{F}_{\bar{q}}(x_2,b;\mu,\zeta) f_{\rm NP}(x_1,b,\zeta)
     f_{\rm NP}(x_2,b,\zeta)\,,
   \end{array}

where I have explicitly reinstated the dependence of the function
:math:`\mathcal{P}` and its derivative with respect to :math:`q_T`,
:math:`\mathcal{P}'`, on :math:`Q` and :math:`y`. In the square bracket
in Eq. :eq:`eq:KexplicitCuts`, the minus sign
applies when :math:`q_T` is the upper bound of the bin and the plus sign
when it is the lower bound (see Eq. :eq:`eq:lastexpP`).
As discussed below, when integrating over bins in :math:`Q` and
:math:`y`, one should also integrate the functions :math:`\mathcal{P}`
and :math:`\mathcal{P}'`. However, we will argue that, in the
interpolation procedure discussed below, these functions can be
extracted from the integrals in :math:`Q` and :math:`y` in a proper
manner in such a way to avoid computing the expensive function
:math:`\mathcal{P}` many times and, moreover, simplify enormously the
structure of the resulting interpolation tables.

.. container::
   :name: on-the-position-of-the-peak-of-the-q_t-distribution

   .. rubric:: On the position of the peak of the :math:`q_T`
      distribution
      :name: on-the-position-of-the-peak-of-the-q_t-distribution

It is interesting at this point to take a short detour to discuss the
position of the peak on the distribution in :math:`q_T` of the cross
section in Eq. :eq:`eq:crosssection`. The peak can
be located by setting the derivative in :math:`q_T` of the cross section
equal to zero. To do so, we use another property of Bessel’s functions:

.. math:: \frac{dJ_0(x)}{dx} = -J_1(x)\,.

\ Using this relation, it is easy to see that:

.. math::
   \begin{array}{l}
   \displaystyle 0 = \frac{d}{dq_T}  \left[\frac{d\sigma}{dQ dy dq_T}\right]
     =\\
   \\
   \displaystyle  \frac{8\pi\alpha^2}{9 Q^3} H(Q,\mu) 
     \int_0^\infty db\,b \left[J_0(bq_T) -bq_TJ_1(bq_T)\right] 
   \sum_q C_q(Q)\overline{F}_q(x_1,b_*(b);\mu,\zeta)
     \overline{F}_{\bar{q}}(x_2,b_*(b);\mu,\zeta)\\
   \\
   \displaystyle \times f_{\rm NP}(x_1,b,\zeta)
     f_{\rm NP}(x_2,b,\zeta)\,,
   \end{array}

that is equivalent to require that:

.. math::
   \int_0^\infty db\,b \left[J_0(bq_T) -bq_TJ_1(bq_T)\right] 
     \sum_q C_q(Q)\overline{F}_q(x_1,b_*(b);\mu,\zeta) \overline{F}_{\bar{q}}(x_2,b_*(b);\mu,\zeta) f_{\rm NP}(x_1,b,\zeta)
     f_{\rm NP}(x_2,b,\zeta) = 0\,.

The integral above can be solved numerically using the technique
discussed above and the value of :math:`q_T` that satisfies this
equation represents the position of the peak of the :math:`q_T`
distribution.

.. container::
   :name: sec:QyInt

   .. rubric:: Integrating over :math:`Q` and :math:`y`
      :name: sec:QyInt

As a final step, we need to perform the integrals over :math:`Q` and
:math:`y` defined in
Eq. :eq:`eq:Intcrosssection`. To compute these
integrals we can only rely on numerical methods. Having reduced the
integration in :math:`q_T` to the difference of the two terms in the
r.h.s. of Eq. :eq:`eq:primitive`( [4]_), we can
concentrate on integrating the function :math:`K` over :math:`Q` and
:math:`y` for a fixed value of :math:`q_T`:

.. math::
   :label: eq:integralK
   \widetilde{K}(q_T)=\int_{Q_{\rm min}}^{Q_{\rm max}}dQ \int_{y_{\rm
       min}}^{y_{\rm max}}dy\,K(Q,y,q_T)\,,

\ such that:

.. math:: \widetilde{\sigma} = \widetilde{K} (q_{T,\rm max})- \widetilde{K} (q_{T,\rm min})\,.

To this purpose, it is convenient to make explicit the dependence of
:math:`x_1` and :math:`x_2` on :math:`Q` and :math:`y` using
Eq. :eq:`eq:Bjorkenx12`. In addition, for the sake of
simplicity we will identify the scales :math:`\mu` and
:math:`\sqrt{\zeta}` with :math:`Q` (possible scale variations can be
easily reinstated at a later stage) and thus drop one of the arguments
from the TMD distributions :math:`\overline{F}` and from the hard factor
:math:`H`. This yields:

.. math::
   :label: eq:finalintegral
   \begin{array}{rcl}
   \displaystyle  \widetilde{K}(q_T) &=& \displaystyle \frac{8\pi q_T}{9} \int_0^\infty db\, J_1(bq_T)
     \int_{Q_{\rm min}}^{Q_{\rm max}}
     dQ \int_{e^{y_{\rm
       min}}}^{e^{y_{\rm max}}}\frac{d\xi}{\xi}\\
   \\
   &\times& \displaystyle 
                            \frac{1}{Q^3} \alpha^2(Q) H(Q)\sum_q C_q(Q)\overline{F}_q\left(\frac{Q}{\sqrt{s}}\xi,b_*(b);Q\right)
                            \overline{F}_{\bar{q}}\left(\frac{Q}{\sqrt{s}}\frac1{\xi},b_*(b);Q\right) \\
   \\
   &\times& \displaystyle f_{\rm NP}\left(\frac{Q}{\sqrt{s}}\xi,b;Q\right)
     f_{\rm NP}\left(\frac{Q}{\sqrt{s}}\frac1{\xi},b;Q\right)\,,
   \end{array}

where we have performed the change of variable :math:`e^{y} = \xi`. Now
we define one grid in :math:`\xi`, :math:`\{\xi_\alpha\}` with
:math:`\alpha=0,\dots,N_\xi`, and one grid in :math:`Q`,
:math:`\{Q_\tau\}` with :math:`\tau=0,\dots,N_Q`, each of which with a
set of interpolating functions :math:`\mathcal{I}` associated. In
addition, the grids are such that: :math:`\xi_0 = e^{y_{\rm min}}` and
:math:`\xi_{N_\xi} = e^{y_{\rm max}}`, and :math:`Q_0 = Q_{\rm min}` and
:math:`Q_{N_Q} = Q_{\rm max}`. This allows us to interpolate the pair of
functions :math:`f_{\rm NP}` in
Eq. :eq:`eq:finalintegral` for generic values of
:math:`\xi` and :math:`Q` as:

.. math::
   :label: eq:interpolation
   f_{\rm NP}\left(\frac{Q}{\sqrt{s}}\xi,b;Q\right) f_{\rm NP}\left(\frac{Q}{\sqrt{s}}\frac1{\xi},b;Q\right) \simeq \sum_{\alpha=0}^{N_\xi}\sum_{\tau=0}^{N_Q}\mathcal{I}_\alpha(\xi)\mathcal{I}_\tau(Q) f_{\rm NP}\left(\frac{Q_\tau}{\sqrt{s}}\xi_\alpha,b;Q_\tau\right) f_{\rm NP}\left(\frac{Q_\tau}{\sqrt{s}}\frac1{\xi_\alpha},b;Q_\tau\right)\,.

Plugging the equation above into
Eq. :eq:`eq:finalintegral` we obtain:

.. math::
   \begin{array}{rcl}
   \displaystyle  \widetilde{K}(q_T) &\simeq& \displaystyle \frac{8\pi q_T}{9} \int_0^\infty db\, J_1(bq_T)
     \sum_{\tau=0}^{N_Q}\sum_{\alpha=0}^{N_\xi}\Bigg[\int_{Q_{\rm min}}^{Q_{\rm max}}dQ\,\mathcal{I}_\tau(Q)\, 
     \frac{1}{Q^3} \alpha^2(Q) H(Q) 
     \\
   \\
   &\times& \displaystyle 
                            \int_{e^{y_{\rm
       min}}}^{e^{y_{\rm max}}}d\xi\,\mathcal{I}_\alpha(\xi)\,\frac{1}{\xi} \sum_q C_q(Q)\overline{F}_q\left(\frac{Q}{\sqrt{s}}\xi,b_*(b);Q\right)
                            \overline{F}_{\bar{q}}\left(\frac{Q}{\sqrt{s}}\frac1{\xi},b_*(b);Q\right)\Bigg] \\
   \\
   &\times& \displaystyle f_{\rm NP}\left(\frac{Q_\tau}{\sqrt{s}}\xi_\alpha,b;Q_\tau\right) f_{\rm NP}\left(\frac{Q_\tau}{\sqrt{s}}\frac1{\xi_\alpha},b;Q_\tau\right)\,.
   \end{array}

Finally, the integration over :math:`b` can be performed using the Ogata
quadrature as discussed above, so that:

.. math::
   \begin{array}{rcl}
   \displaystyle  \widetilde{K}(q_T) &\simeq& \displaystyle \sum_{n=1}^N
     \sum_{\tau=0}^{N_Q}\sum_{\alpha=0}^{N_\xi}\Bigg[\frac{8\pi}{9} w_n^{(1)}\int_{Q_{\rm min}}^{Q_{\rm max}}dQ\,\mathcal{I}_\tau(Q)\, 
     \frac{1}{Q^3} \alpha^2(Q) H(Q) 
   \\
   \\
   &\times& \displaystyle 
                            \int_{e^{y_{\rm
       min}}}^{e^{y_{\rm max}}}d\xi\,\mathcal{I}_\alpha(\xi)\,\frac{1}{\xi} \sum_q C_q(Q)\overline{F}_q\left(\frac{Q}{\sqrt{s}}\xi,b_*\left(\frac{z_n}{q_T}\right);Q\right)
                            \overline{F}_{\bar{q}}\left(\frac{Q}{\sqrt{s}}\frac1{\xi},b_*\left(\frac{z_n}{q_T}\right);Q\right)\Bigg] \\
   \\
   &\times& \displaystyle f_{\rm NP}\left(\frac{Q_\tau}{\sqrt{s}}\xi_\alpha,\frac{z_n}{q_T};Q_\tau\right) f_{\rm NP}\left(\frac{Q_\tau}{\sqrt{s}}\frac1{\xi_\alpha},\frac{z_n}{q_T};Q_\tau\right)\,.
   \end{array}

In conclusion, if we define:

.. math::
   :label: eq:weights
   \begin{array}{rcl}
     \displaystyle  W_{n\tau\alpha}(q_T) & \equiv & \displaystyle w_n^{(1)}\frac{8\pi}{9} \int_{Q_{\rm min}}^{Q_{\rm max}}dQ\,\mathcal{I}_\tau(Q)\, 
                                               \frac{\alpha^2(Q)}{Q^3} H(Q) 
                                               \\
     \\
                                    &\times& \displaystyle 
                                             \int_{e^{y_{\rm
                                             min}}}^{e^{y_{\rm max}}}d\xi\,\mathcal{I}_\alpha(\xi)\,\frac{1}{\xi} \sum_q C_q(Q)\overline{F}_q\left(\frac{Q}{\sqrt{s}}\xi,b_*\left(\frac{z_n}{q_T}\right);Q\right)
                                             \overline{F}_{\bar{q}}\left(\frac{Q}{\sqrt{s}}\frac1{\xi},b_*\left(\frac{z_n}{q_T}\right);Q\right)\,,
   \end{array}

the quantity :math:`\widetilde{K}(q_T)` can be computed as:

.. math::
   :label: eq:finalinterpolated
   \widetilde{K}(q_T) \simeq \sum_{n=1}^N
     \sum_{\tau=0}^{N_Q}\sum_{\alpha=0}^{N_\xi} W_{n\tau\alpha}(q_T) f_{\rm NP}\left(\frac{Q_\tau}{\sqrt{s}}\xi_\alpha,\frac{z_n}{q_T};Q_\tau\right) f_{\rm NP}\left(\frac{Q_\tau}{\sqrt{s}}\frac1{\xi_\alpha},\frac{z_n}{q_T};Q_\tau\right)\,.

The advantage of
Eq. :eq:`eq:finalinterpolated` is that the
weights :math:`W_{n\alpha\tau}`, that clearly depend on :math:`q_T` but
also on the intervals :math:`[Q_{\rm min}:Q_{\rm max}]` and
:math:`[y_{\rm min}:y_{\rm max}]`, can be precomputed once and for all
for each of the experimental points included in a fit and used to
determine the function :math:`f_{\rm NP}`. This provides a fast tool for
the computation of predictions that makes the extraction of the
non-perturbative part of the TMDs much easier.

It is now time to discuss how the weights defined in
Eq. :eq:`eq:weights` are affected by the presence of
cuts as discussed in Sect. `2.1.1 <#sec:kincuts>`__. In principle, the
function between square brackets in
Eq. :eq:`eq:KexplicitCuts` should be inside the
integrals in Eq. :eq:`eq:weights` and integrated over
the variable :math:`Q` and :math:`\xi=e^y`. However, this turns out to
be numerically problematic because the phase-space-reduction function
:math:`\mathcal{P}` is expensive to compute. On top of this, the fact
that the factor between square brackets in
Eq. :eq:`eq:KexplicitCuts` depends on whether
:math:`q_T` is a lower or an upper integration bound would lead to a
duplication of the weights to compute. In order to simplify the
computation, we assume that the function :math:`\mathcal{P}` and its
derivative :math:`\mathcal{P}'` are slowly varying functions of
:math:`Q` and :math:`y` over the typical grid interval of the grids in
:math:`Q` and :math:`\xi`. In addition, the interpolating functions
:math:`\mathcal{I}_\tau(Q)` and :math:`\mathcal{I}_\alpha(\xi)` are
strongly peaked at :math:`Q_\tau` and :math:`\xi_\alpha`, respectively.
These considerations allow us to avoid integrating explicitly
:math:`\mathcal{P}` and :math:`\mathcal{P}'` over :math:`Q` and
:math:`\xi` and to replace the weights in
Eq. :eq:`eq:weights` with:

.. math::
   :label: eq:weightswithcuts
     \displaystyle  W_{n\tau\alpha}(q_T) \rightarrow \left[\mathcal{P}\left(Q_\tau,\ln(\xi_\alpha),q_{T}\right) \pm
     \mathcal{P}'\left(Q_\tau,\ln(\xi_\alpha),q_{T}\right)\Delta q_T\right] W_{n\tau\alpha}(q_T)\,.

At the end of the day, the only additional information required to
implement cuts on the final state is the value of the
phase-space-reduction function :math:`\mathcal{P}` and its derivative
:math:`\mathcal{P}'` on all points of the bidimensional grid in
:math:`Q` and :math:`\xi` for all :math:`q_T` bin bounds.
Eq. :eq:`eq:weightswithcuts` will then allow one
to use the weights computed over the full phase space. We will check the
accuracy of this procedure by comparing it to the explicit integration.

.. container::
   :name: cross-section-differential-in-x_f

   .. rubric:: Cross section differential in :math:`x_F`
      :name: cross-section-differential-in-x_f

In some cases, the Drell-Yan differential cross section may be presented
as differential in the invariant mass of the lepton pair :math:`Q` and,
instead of the rapidity :math:`y`, of the Feynman variable :math:`x_F`
defined as:

.. math::
   x_F = 
     \frac{Q}{\sqrt{s}}\left(e^{y} - e^{-y}\right) =
     \frac{2Q}{\sqrt{s}}\sinh y = x_1-x_2\,,

\ so that:

.. math:: \frac{dx_F}{dy} = \frac{2Q}{\sqrt{s}}\cosh y=x_1+x_2\,.

Therefore:

.. math::
   \frac{d\sigma}{dQ dx_F dq_T} =
     \frac{dy}{dx_F}\frac{d\sigma}{dQ dy dq_T}=
   \frac{\sqrt{s}}{2Q\cosh y}\frac{d\sigma}{dQ dy dq_T}=\frac1{x_1+x_2}\frac{d\sigma}{dQ dy dq_T}

with:

.. math::
   y(x_F,Q) =
     \sinh^{-1}\left(\frac{x_F\sqrt{s}}{2Q}\right) =
     \ln\left[\frac{\sqrt{s}}{2Q}\left(x_F+\sqrt{x_F^2 + \frac{4Q^2}{s}}\right)\right]\,,

so that:

.. math::
   :label: eq:x12ofxFQ
   x_1 = \frac12\left(x_F+\sqrt{x_F^2 +
       \frac{4Q^2}{s}}\right)\quad\mbox{and}\quad x_2 = \frac{Q^2}{sx_1}\,.

Therefore, we can compute the integral:

.. math::
   \widetilde{I}(q_T)=\int_{Q_{\rm min}}^{Q_{\rm max}}dQ \int_{x_{F,\rm
       min}}^{x_{F,\rm max}}dx_F\,I(Q,x_F,q_T)\,,

where :math:`I` is the primitive in :math:`q_T` of the cross section
differential in :math:`x_F`:

.. math:: I(Q,x_F,q_T) = \int dq_T\left[\frac{d\sigma}{dQ dx_F dq_T}\right]\,,

following the same steps of Sect. `2.3 <#sec:QyInt>`__. This leads to:

.. math::
   \begin{array}{rcl}
   \displaystyle  \widetilde{I}(q_T) &\simeq& \displaystyle \sum_{n=1}^N
     \sum_{\tau=0}^{N_Q}\sum_{\alpha=0}^{N_x}\overline{W}_{n\tau\alpha}(q_T) f_{\rm NP}\left(x_{1,\alpha\tau},\frac{z_n}{q_T};Q_\tau\right) f_{\rm NP}\left(x_{2,\alpha\tau},\frac{z_n}{q_T};Q_\tau\right)\,,
   \end{array}

with:

.. math::
   :label: eq:xFtens
   \begin{array}{rcl}
   \displaystyle  \overline{W}_{n\tau\alpha}(q_T) &\equiv& \displaystyle w_n^{(1)}\frac{8\pi}{9} \int_{Q_{\rm min}}^{Q_{\rm max}}dQ\,\mathcal{I}_\tau(Q)\, 
     \frac{1}{Q^3} \alpha^2(Q) H(Q) 
   \\
   \\
   &\times& \displaystyle 
                            \int_{x_{F,\rm
       min}}^{x_{F,\rm max}}dx_F\,\mathcal{I}_\alpha(x_F)\,\frac{1}{x_1+x_2} \sum_q C_q(Q)\overline{F}_q\left(x_1,b_*\left(\frac{z_n}{q_T}\right);Q\right)
                            \overline{F}_{\bar{q}}\left(x_2,b_*\left(\frac{z_n}{q_T}\right);Q\right)\,,
   \end{array}

where :math:`x_1` and :math:`x_2` are functions of :math:`x_F` and
:math:`Q` through Eq. :eq:`eq:x12ofxFQ`. In addition,
we have defined a grid in :math:`x_F`, :math:`\{x_{F,\alpha}\}` with
:math:`\alpha = 0,\dots,N_x`, that allowed us to define
:math:`x_{1(2),\alpha\tau}\equiv x_{1(2)}(x_{F,\alpha},Q_\tau)`.

.. container::
   :name: flavour-dependence

   .. rubric:: Flavour dependence
      :name: flavour-dependence

It may be advantageous to introduce a flavour dependence of the
non-perturbative contributions to TMDs. This can be easily done by
observing that the tensor :math:`W_{n\tau\alpha}` defined in
Eq. :eq:`eq:weights` can be decomposed as [5]_:

.. math:: W_{n\tau\alpha}(q_T) = \sum_q W_{n\tau\alpha}^{(q)} (q_T)\,,

\ with:

.. math::
   :label: eq:weightsfl
   \begin{array}{rcl}
     \displaystyle  W_{n\tau\alpha}^{(q)}(q_T) & \equiv & \displaystyle w_n^{(1)}\frac{8\pi}{9} \int_{Q_{\rm min}}^{Q_{\rm max}}dQ\,\mathcal{I}_\tau(Q)\, 
                                               \frac{\alpha^2(Q)}{Q^3} H(Q) C_q(Q)
                                               \\
     \\
                                    &\times& \displaystyle 
                                             \int_{e^{y_{\rm
                                             min}}}^{e^{y_{\rm max}}}d\xi\,\mathcal{I}_\alpha(\xi)\,\frac{1}{\xi} \overline{F}_q\left(\frac{Q}{\sqrt{s}}\xi,b_*\left(\frac{z_n}{q_T}\right);Q\right)
                                             \overline{F}_{\bar{q}}\left(\frac{Q}{\sqrt{s}}\frac1{\xi},b_*\left(\frac{z_n}{q_T}\right);Q\right)\,.
   \end{array}

This allows for an independent parameterisation of the non-perturbative
contribution such that
Eq. :eq:`eq:finalinterpolated` can be written
as:

.. math::
   :label: eq:finalinterpolatedfl
   \widetilde{K}(q_T) \simeq \sum_q\sum_{n=1}^N
     \sum_{\tau=0}^{N_Q}\sum_{\alpha=0}^{N_\xi} W_{n\tau\alpha}^{(q)}(q_T) f_{\rm NP}^{(q)}\left(\frac{Q_\tau}{\sqrt{s}}\xi_\alpha,\frac{z_n}{q_T};Q_\tau\right) f_{\rm NP}^{(q)}\left(\frac{Q_\tau}{\sqrt{s}}\frac1{\xi_\alpha},\frac{z_n}{q_T};Q_\tau\right)\,,

where :math:`f_{\rm NP}^{(q)}` parametrises the non-perturbative
component of the TMD with flavour :math:`q`.

.. container::
   :name: gradient-with-respect-to-the-free-parameters

   .. rubric:: Gradient with respect to the free parameters
      :name: gradient-with-respect-to-the-free-parameters

A very appealing implication of the computation of cross section in
terms of precomputed table as in
Eqs. :eq:`eq:finalinterpolated`
and :eq:`eq:finalinterpolatedfl` is the fact
that it exposes the free parameters of the non-perturbative functions.
To be more specific, the non-perturbative function :math:`f_{\rm NP}`,
on top of being a function of :math:`x`, :math:`b`, and :math:`\zeta`,
depends parameterically on a set of :math:`N_p` parameters
:math:`\{\theta_k\}`, :math:`k=1,\dots,N_p`, that are typically
determined by fits to data, in other words:

.. math:: f_{\rm NP}\equiv f_{\rm NP}\left(x,b,\zeta;\{\theta_k\}\right)\,.

Now, when performing a fit, it is very useful to be able to compute the
derivative of the figure of merit (usually the :math:`\chi^2`) with
respect to the parameters to be determined. In turn, this immediately
implies being able to compute the derivative of the observables.
Referring to Eq. :eq:`eq:finalinterpolated`,
the relevant quantity is:

.. math::
   \frac{d\widetilde{K}}{d\theta_k}=
   \sum_{n=1}^N
     \sum_{\tau=0}^{N_Q}\sum_{\alpha=0}^{N_\xi} W_{n\tau\alpha}(q_T) \left[\frac{d f_{\rm NP}^{(1)}}{d\theta_k} f_{\rm NP}^{(2)}+f_{\rm NP}^{(1)}\frac{d f_{\rm NP}^{(2)}}{d\theta_k} \right]\,,

where :math:`f_{\rm NP}^{(1)}` and :math:`f_{\rm NP}^{(2)}` refer to the
non-perturbative function :math:`f_{\rm NP}` computed in :math:`x_1` and
:math:`x_2`, respectively. It is thus clear that the derivatives w.r.t.
the free parameters penetrates the observable. Since in most cases the
derivative of :math:`f_{\rm NP}` can be computed analytically, this
allows one to compute the gradient of the figure of merit analytically.
This potentially makes any fit much simpler.

.. container::
   :name: narrow-width-approximation

   .. rubric:: Narrow-width approximation
      :name: narrow-width-approximation

A possible alternative to the numerical integration in :math:`Q` when
the integration region includes the :math:`Z`-peak region is the
so-called narrow-width approximation (NWA). In the NWA one assumes that
the width of the :math:`Z` boson, :math:`\Gamma_Z`, is much smaller than
its mass, :math:`M_Z`. This way one can approximate the peaked behaviour
of the couplings :math:`C_q(Q)` around :math:`Q=M_Z` with a
:math:`\delta`-function, *i.e.* :math:`C_q(Q)\sim \delta(Q^2-M_Z^2)`.
Therefore, the integration over :math:`Q` can be done analytically. The
exact structure of the electroweak couplings is the following:

.. math::
   :label: eq:fullcoup
   C_q(Q) = e_q^2 - 2 e_q V_q V_e \chi_1(Q) + (V_e^2 + A_e^2)(V_q^2 + A_q^2)\chi_2(Q)\,,

with:

.. math::
   \begin{array}{l}
   \displaystyle \chi_1(Q) = \frac{1}{4 \sin^2\theta_W \cos^2\theta_W } \frac{Q^2 ( Q^2 -  M_Z^2 )}{ (Q^2 - M_Z^2)^2 + M_Z^2 \Gamma_Z^2} \,,\\
   \displaystyle \chi_2(Q) = \frac{1}{16 \sin^4\theta_W\cos^4\theta_W} \frac{Q^4}{ (Q^2 - M_Z^2)^2 + M_Z^2 \Gamma_Z^2} \,.
   \end{array}

In the limit :math:`\Gamma_Z/M_Z\rightarrow 0`, the leading contribution
to the coupling in Eq. :eq:`eq:fullcoup` comes from the
region :math:`Q\simeq M_Z` and is that proportional to :math:`\chi_2`:

.. math::
   :label: eq:partlead
   C_q(Q) \simeq (V_e^2 + A_e^2)(V_q^2 + A_q^2)\chi_2(Q)\,,\quad Q\simeq M_Z\,.

In addition, in this limit one can show that:

.. math::
   :label: eq:breitwigner
   \frac{1}{ (Q^2 - M_Z^2)^2 + M_Z^2 \Gamma_Z^2}\rightarrow
   \frac{\pi}{M_Z\Gamma_Z}\delta(Q^2-M_Z^2) = \frac{\pi}{2M_Z^2\Gamma_Z}\delta(Q-M_Z)\,.

Therefore, considering that:

.. math:: \Gamma_Z = \frac{\alpha M_Z}{\sin^2\theta_W \cos^2\theta_W}\,,

the electroweak couplings in the NWA have the following form:

.. math::
   C_q(Q) \simeq \frac{\pi M_Z (V_e^2 + A_e^2)(V_q^2 + A_q^2) }{32 \alpha
       \sin^2\theta_W\cos^2\theta_W} \delta(Q-M_Z)=\widetilde{C}_q(Q) \delta(Q-M_Z)\,.

Therefore, using Eq. :eq:`eq:partlead` the integral of
the cross section over :math:`Q` under the condition that
:math:`Q_{\rm min}<M_Z<Q_{\rm max}` has the consequence of adjusting the
couplings and of setting :math:`Q=M_Z` in the computation. This yields:

.. math::
   \int_{Q_{\rm min}}^{Q_{\rm max}}dQ\,\frac{d\sigma}{dQ dy dq_T} =
    \frac{16\pi\alpha^2q_T}{9 M_Z^3} H(M_Z,M_Z) \sum_q \widetilde{C}_q(M_Z)
     I_{q\bar{q}}(x_1,x_2,q_T;M_Z,M_Z^2)\,,

where we are also assuming that :math:`\mu=\sqrt{\zeta}=M_Z`. As a final
step, one may want to let the :math:`Z` boson decay into leptons. At
leading order in the EW sector and assuming an equal decay rate for
electrons, muons, and tauons, this can be done by multiplying the cross
section above by three times the branching ratio for the :math:`Z`
decaying into any pair of leptons,
:math:`3\mbox{Br}(Z\rightarrow \ell^+\ell^-)`.

.. _app:OgataQuadrature:

Ogata quadrature
================

In this section we limit ourselves to write the formulas for the
computation of the unscaled coordinates :math:`z_n^{(\nu)}` and weights
:math:`w_n^{(\nu)}` required to compute the following integral:

.. math::
   :label: eq:OgataQuadMast
   I_\nu(q_T)=\int_0^\infty db J_\nu(bq_T) f\left(b\right) =
   \frac1{q_T}\int_0^\infty d\bar{b} J_\nu(\bar{b})
   f\left(\frac{\bar{b}}{q_T}\right) \simeq
   \frac{1}{q_T}\sum_{n=1}^\infty
   w_n^{(\nu)}f\left(\frac{z_n^{(\nu)}}{q_T}\right)\quad \nu =0,1,\dots\,,

using the Ogata-quadrature algorithm. More details can be found in
Ref. (Ogata, n.d.). There relevant formulas are:

.. math::
   \begin{array}{l}
   \displaystyle z_n^{(\nu)} = \frac{\pi}{h}  \psi\left(\frac{h\xi_{\nu
     n}}{\pi}\right)\,,\\
   \\
   \displaystyle w_n^{(\nu)}  = \pi\frac{Y_\nu(\xi_{\nu
     n})}{J_{\nu+1}(\xi_{\nu n})}  J_\nu(z_n^{(\nu)})  \psi'\left(\frac{h\xi_{\nu
     n}}{\pi}\right)\,.
   \end{array}

where:

-  :math:`h` is a free parameter of the algorithm that has to be
   typically small (we choose :math:`h = 10^{-3}`),

-  :math:`\xi_{\nu n}` are the zero’s of :math:`J_\nu`, *i.e.*
   :math:`J_\nu(\xi_{\nu n}) = 0` :math:`\forall\, n`,

-  :math:`J_\nu` and :math:`Y_\nu` are the Bessel functions of first and
   second kind, respectively, of degree :math:`\nu`,

-  :math:`\psi` is the following function:

   .. math:: \psi(t) = t\tanh\left(\frac{\pi}{2}\sinh t\right)

   \ and its derivative:

   .. math::
      \psi'(t) =  \frac{\pi t \cosh t + \sinh( \pi \sinh t ) }{1 +
          \cosh( \pi \sinh t  ) }\,.

Cuts on the final-state leptons
===============================

In this section we derive explicitly the phase-space reduction factor
:math:`\mathcal{P}` introduced in Sect. `2.1.1 <#sec:kincuts>`__. This
factor is defined as:

.. math::
   :label: eq:PSredDef
   \mathcal{P}(Q,y,q_T) = \mathcal{P}(q) = \frac{\displaystyle \int_{\mbox{\footnotesize fid.
       reg.}}d^4p_1 d^4p_2 \,\delta(p_1^2) \delta(p_2^2)\theta(p_{1,0}) \theta(p_{2,0})\delta^{(4)}(p_1+p_2-q) g_{\mu\nu}L^{\mu\nu}(p_1,p_2)}{\displaystyle \int d^4p_1 d^4p_2\, \delta(p_1^2) \delta(p_2^2) \theta(p_{1,0}) \theta(p_{2,0})\delta^{(4)}(p_1+p_2-q) g_{\mu\nu}L^{\mu\nu}(p_1,p_2)}\,,

where :math:`p_1` and :math:`p_2` are the four-momenta of the outgoing
leptons and :math:`L^{\mu\nu}` is the leptonic tensor that, assuming
massless leptons, reads:

.. math::
   :label: eq:lepttens
   L^{\mu\nu}(p_1,p_2) = 4(p_1^{\mu}p_2^{\nu}+p_2^{\mu}p_1^{\nu}-g^{\mu\nu}p_1p_2)\,,

so that:

.. math:: g_{\mu\nu}L^{\mu\nu}(p_1,p_2) = -8(p_1p_2) = -4(p_1+p_2)^2\,.

In the last step we have used the on-shell-ness of the leptons
(:math:`p_1^2=p_2^2=0`). The integral in the denominator of
Eq. :eq:`eq:PSredDef` is restricted to some *fiducial
region*. Finally, we find:

.. math::
   :label: eq:PSredDef2
   \mathcal{P}(q) = \frac{\displaystyle \int_{\mbox{\footnotesize fid.
       reg.}}d^4p_1 d^4p_2 \,\delta(p_1^2) \delta(p_2^2) \theta(p_{1,0}) \theta(p_{2,0})\delta^{(4)}(p_1+p_2-q) (p_1+p_2)^2}{\displaystyle \int d^4p_1 d^4p_2\, \delta(p_1^2) \delta(p_2^2) \theta(p_{1,0}) \theta(p_{2,0})\delta^{(4)}(p_1+p_2-q) (p_1+p_2)^2}\,.

The effect of integrating over the fiducial region can be implemented by
defining a generalised :math:`\theta`-function, :math:`\Phi(p_1,p_2)`,
that is equal to one inside the fiducial region and zero outside. This
allows one to integrate also the numerator of
Eq. :eq:`eq:PSredDef2` over the full phase-space of
the two outgoing leptons:

.. math::
   :label: eq:PSredDef3
   \mathcal{P}(q) = \frac{\displaystyle \int d^4p_1 d^4p_2 \,\delta(p_1^2) \delta(p_2^2) \theta(p_{1,0}) \theta(p_{2,0})\delta^{(4)}(p_1+p_2-q) \Phi(p_1,p_2) (p_1+p_2)^2}{\displaystyle \int d^4p_1 d^4p_2\, \delta(p_1^2) \delta(p_2^2) \theta(p_{1,0}) \theta(p_{2,0})\delta^{(4)}(p_1+p_2-q) (p_1+p_2)^2}\,.

Now we can integrate over one of the outgoing momenta, say :math:`p_2`,
exploiting the momentum-conservation :math:`\delta`-function both in the
numerator and in the denominator. Specifically, the numerator of
Eq. :eq:`eq:PSredDef3` gives:

.. math::
   \begin{array}{c}
   \displaystyle \int d^4p_1 d^4p_2\, \delta(p_1^2)
   \delta(p_2^2) \theta(p_{1,0}) \theta(p_{2,0})\delta^{(4)}(p_1+p_2-q)
     \Phi(p_1,p_2) (p_1+p_2)^2 = \\
   \\
   \displaystyle Q^2 \int d^4p_1 \delta(p_1^2)
   \delta((q-p_1)^2) \theta(p_{1,0})
     \theta(q_0-p_{1,0})\Phi(p_1,q-p_1)\,,
   \end{array}

and likewise in the denominator setting :math:`\Phi(p_1,p_2)=1`.
Finally, renaming :math:`p_1=p`, the phase-space reduction factor reads:

.. math::
   :label: eq:PSredDef4
     \mathcal{P}(q) = \frac{\displaystyle \int d^4p \delta(p^2) \delta((q-p)^2) \theta(p_{0})
       \theta(q_0-p_{0})  \Phi(p,q-p)}{\displaystyle \int d^4p \delta(p^2) \delta((q-p)^2) \theta(p_{0})
       \theta(q_0-p_{0})  }\,.

The :math:`\delta`-functions can now be used to constrain two of the
four components of the momentum :math:`p`. The first,
:math:`\delta(p_2^2)`, is usually used to set the first component of
:math:`p`, the energy, to the on-shell value. Since the leptons are
assumed to be massless, this produces:

.. math::
   :label: eq:phasespacemeasure
   \int d^4p\delta(p^2)\theta(p_0) = \int d^4p\delta(E^2-|\mathbf{p}|^2)\theta(E)=\int\frac{dEd^3\mathbf{p}}{2|\mathbf{p}|}\delta(E-|\mathbf{p}|)=\int\frac{d^3\mathbf{p}}{2|\mathbf{p}|}\,.

Of course, the four-momentum :math:`p` appearing in the rest of the
integrand has to be set on shell (:math:`E=|\mathbf{p}|`). Now we
express the three-dimensional measure :math:`d^3\mathbf{p}` in spherical
coordinates as:

.. math:: d^3\mathbf{p} = |\mathbf{p}|^2d|\mathbf{p}|d(\cos\theta) d\phi\,.

Then we make a change of variable from :math:`(|\mathbf{p}|,\cos\theta)`
to :math:`(|\mathbf{p}_T|,\eta)`: the second set of variables are
exactly those on which kinematic cuts are imposed. We do so by knowing
that:

.. math::
   \left\{
   \begin{array}{l}
   |\mathbf{p}| = |\mathbf{p}_T|\cosh\eta\,,\\
   \cos\theta =\tanh\eta\,.
   \end{array}
   \right.

This leads to:

.. math:: \int\frac{d^3\mathbf{p}}{2 |\mathbf{p}|} = \frac12\int|\mathbf{p}|d|\mathbf{p}|d(\cos\theta) d\phi=\frac12\int|\mathbf{p}_T|d|\mathbf{p}_T|d\eta d\phi=\frac12\int d^2\mathbf{p}_T d\eta\,.

Now we consider the second :math:`\delta`-function:

.. math::
   :label: eq:integralyeah!
   \frac12\int d^2\mathbf{p}_T d\eta\,\delta((q-p)^2)\theta(q_0-p_0)=\frac12\int_{-\infty}^\infty d\eta
   \int_0^{2\pi} d\phi \int_0^\infty|\mathbf{p}_T|d|\mathbf{p}_T|\,\delta(Q^2-2p\cdot q) \theta(q_0-p_0)\,,

being :math:`q^2=Q^2` and :math:`p^2=0`. It is convenient to express the
four-vector :math:`q` in terms of :math:`Q`, :math:`y`, and
:math:`\mathbf{q}_T`:

.. math::
   :label: eq:qexplicit
   q=\left(M\cosh y,\mathbf{q}_T,M\sinh y\right)\,.

with :math:`M=\sqrt{Q^2+|\mathbf{q}_T|^2}`. While:

.. math::
   :label: eq:pexplicit
   p=\left(|\mathbf{p}_T|\cosh\eta,\mathbf{p}_T,|\mathbf{p}_T|\sinh\eta\right)\,,

so that:

.. math::
   p\cdot q=|\mathbf{p}_T|M\left(\cosh\eta \cosh y-\sinh\eta\sinh
     y\right)-\mathbf{p}_T\cdot
   \mathbf{q}_T=|\mathbf{p}_T|M\cosh\left(\eta - y\right)-\mathbf{p}_T\cdot \mathbf{q}_T\,.

We can now assume that the two-dimensional vector :math:`\mathbf{q}_T`
is aligned with the :math:`x` axis so that
:math:`\mathbf{p}_T\cdot \mathbf{q}_T =
|\mathbf{p}_T||\mathbf{q}_T|\cos\phi`\ ( [6]_). Therefore, the argument
of the :math:`\delta`-function in
Eq. (`[eq:integralyeah!] <#eq:integralyeahux21>`__) becomes:

.. math::
   :label: eq:deltaargument
   f(|\mathbf{p}_T|,\eta,\phi) = Q^2-2 |\mathbf{p}_T|\left[M\cosh\left(\eta - y\right)-|\mathbf{q}_T|\cos\phi\right]\,.

and that of the :math:`\vartheta`-function
:math:`M\cosh y-|\mathbf{p}_T|\cosh\eta`. It thus appears convenient to
integrate Eq. (`[eq:integralyeah!] <#eq:integralyeahux21>`__) over
:math:`|\mathbf{p}_T|` first:

.. math::
   :label: eq:firstintegral
   \frac12\int_0^\infty|\mathbf{p}_T|d|\mathbf{p}_T|\,\delta(Q^2-2p\cdot q) \theta(q_0-p_0)=\frac{\overline{p}_T^2}{2Q^2}\vartheta(M\cosh y-\overline{p}_T\cosh\eta)=\frac{\overline{p}_T^2}{2Q^2}\,,

with( [7]_):

.. math::
   :label: eq:overpT
   \overline{p}_T(\cos\phi) = \frac{Q^2}{2 \left[M\cosh\left(\eta - y\right)-|\mathbf{q}_T|\cos\phi\right]}=\frac{Q^2}{2 |\mathbf{q}_T|}\frac1{\left[\frac{M\cosh\left(\eta - y\right)}{|\mathbf{q}_T|}-\cos\phi\right]}\,.

Now we turn to consider the integral in :math:`d\phi`. To this end, the
following relations are useful:

.. math::
   :label: eq:intoverphi
   \int_0^{2\pi}d\phi\, f(\cos\phi) = \int_{-1}^1\frac{dx}{\sqrt{1-x^2}}\left[f(x)+f(-x)\right]\,.

and:

.. math::
   :label: eq:complicatedintegral
   \int \frac{dx}{(a\pm
     x)^2\sqrt{1-x^2}}=\frac{\sqrt{1-x^2}}{(a^2-1)(x\pm
     a)}\pm\frac{a}{(a^2-1)^{3/2}}\tan^{-1}\left(\frac{1\pm ax}{\sqrt{a^2-1}\sqrt{1-x^2}}\right)\,.

The last integral is such that:

.. math::
   :label: eq:defintoverphi
   \int_{-1}^{1} \frac{dx}{(a\pm x)^2\sqrt{1-x^2}}=\frac{\pi a}{(a^2-1)^{3/2}}\,.

We now use
Eqs. :eq:`eq:intoverphi`-:eq:`eq:defintoverphi`
to compute:

.. math::
   \begin{array}{rcl}
   &&\displaystyle
     \frac{1}{2Q^2}\int_0^{2\pi}d\phi\,[\overline{p}_T(\cos\phi)]^2 =
                                                                       \displaystyle
                                                                       \frac{Q^2}{4
                                                                       |\mathbf{q}_T|^2}\int_0^{2\pi}\frac{d\phi}{\left[\frac{M\cosh\left(\eta
                                                                       -
                                                                       y\right)}{|\mathbf{q}_T|}-\cos\phi\right]^2}\\
   \\
   &=&\displaystyle \frac{Q^2}{8|\mathbf{q}_T|^2}\int_{-1}^{1}\frac{dx}{\sqrt{1-x^2}}\left[\frac{1}{\left(\frac{M\cosh\left(\eta- y\right)}{|\mathbf{q}_T|}-x\right)^2}+\frac{1}{\left(\frac{M\cosh\left(\eta- y\right)}{|\mathbf{q}_T|}+x\right)^2}\right]\\
   \\
   &=&\displaystyle \frac{Q^2}{8}\Bigg\{\frac{|\mathbf{q}_T|^2 x\sqrt{1-x^2}}{(M ^2\cosh ^2\left(\eta- y\right)-|\mathbf{q}_T|^2)(x^2 |\mathbf{q}_T|^2-
     M ^2\cosh ^2\left(\eta- y\right))}\\
   \\
   &-&\displaystyle \frac{M\cosh\left(\eta-
       y\right)}{(M^2\cosh^2\left(\eta-
       y\right)-|\mathbf{q}_T|^2)^{3/2}}\Bigg[\tan^{-1}\left(\frac{|\mathbf{q}_T|-
       xM\cosh\left(\eta-y\right)}{\sqrt{(M^2\cosh^2\left(\eta-y\right)-|\mathbf{q}_T|^2)}\sqrt{1-x^2}}\right)\\
   \\
   &-&\displaystyle\tan^{-1}\left(\frac{|\mathbf{q}_T|+
       xM\cosh\left(\eta-y\right)}{\sqrt{(M^2\cosh^2\left(\eta-y\right)-|\mathbf{q}_T|^2)}\sqrt{1-x^2}}\right)\Bigg]\Bigg\}_{-1}^{1}\\
   \\
   &=&\displaystyle \frac{\pi
     Q^2M\cosh\left(\eta -
       y\right)}{4(M^2\cosh^2\left(\eta -
       y\right)-|\mathbf{q}_T|^2)^{3/2}}
   \end{array}

We can go further and solve also the integral in :math:`\eta`:

.. math::
   :label: eq:remarkableintegral
   \begin{array}{l}
     \displaystyle \int d^4p \delta(p^2) \delta((q-p)^2) \theta(p_{0})
     \theta(q_0-p_{0})=\int_{-\infty}^\infty
     d\eta\frac{\pi Q^2M\cosh(\eta-y)}{4(M^2\cosh^2(\eta-y)-|\mathbf{q}_T|^2)^{3/2}}= \\
     \\
     \displaystyle\frac{\pi }{4}\frac{Q^2}{M^2}\int_{-\infty}^\infty
     \frac{d(\sinh\eta)}{\left(\sinh^2(\eta-y)+\frac{Q^2}{M^2}\right)^{3/2}}= \frac{\pi}{4} \frac{Q^2}{M^2}\left[\frac{M^2}{Q^2}\frac{\sinh\eta}{\sqrt{\sinh^2\eta+\frac{Q^2}{M^2}}}\right]_{-\infty}^{\infty}= \frac{\pi}{2}\,.
   \end{array}

Remarkably, this result gives us the denominator of
Eq. :eq:`eq:PSredDef4`. We now need to compute the
numerator by inserting the appropriate function :math:`\Phi`. In our
case, the kinematic cuts are identical for the outgoing leptons and
read:

.. math::
   \eta_{\rm
     min} < \eta_{1(2)} < \eta_{\rm max}\quad\mbox{and}\quad |\mathbf{p}_{T,1(2)}| > p_{T,\rm min}\,.

Therefore, the function :math:`\Phi` factorises into two identical
functions as:

.. math:: \Phi(p_1,p_2) = \Theta(p_1)\Theta(p_2)\,,

with:

.. math:: \Theta(p) = \vartheta(\eta - \eta_{\rm min})\vartheta(\eta_{\rm max}-\eta) \vartheta(|\mathbf{p}_{T}| - p_{T,\rm min}) \,.

Regerring to Eq. :eq:`eq:PSredDef4`, and considering
that:

.. math:: q-p=\left(M\cosh y-|\mathbf{p}_T|\cosh\eta,\mathbf{q}_T-\mathbf{p}_T,M\sinh y-|\mathbf{p}_T|\sinh\eta\right)\,.

we thus have:

.. math::
   :label: eq:intdomain
   \begin{array}{ll}
   &\Phi(p,q-p) = \Theta(p) \Theta(q-p)=\\
   \\
   &\displaystyle 
   \vartheta(\eta - \eta_{\rm min}) \vartheta(\eta_{\rm max}-\eta)
     \times\\
   \\
   &\displaystyle \vartheta(|\mathbf{p}_{T}| - p_{T,\rm min})\times\\
   \\
   &\displaystyle \vartheta\left(\frac12\ln\left(\frac{M\cosh y-|\mathbf{p}_T|\cosh\eta+M\sinh y-|\mathbf{p}_T|\sinh\eta}{M\cosh y-|\mathbf{p}_T|\cosh\eta-M\sinh y+|\mathbf{p}_T|\sinh\eta}\right)-\eta_{\rm min}\right)\times\\
   \\
   &\displaystyle \vartheta\left(\eta_{\rm
     max}-\frac12\ln\left(\frac{M\cosh y-|\mathbf{p}_T|\cosh\eta+M\sinh
     y-|\mathbf{p}_T|\sinh\eta}{M\cosh y-|\mathbf{p}_T|\cosh\eta-M\sinh
     y+|\mathbf{p}_T|\sinh\eta}\right)\right)\times\\
   \\
   &\vartheta(|\mathbf{q}_{T}-\mathbf{p}_{T}| - p_{T,\rm min})=\\
   \\
   1):\quad&\displaystyle \vartheta(\eta-\eta_{\rm min}) \times\vartheta(\eta_{\rm max}-\eta) \times\\
   \\
   2):\quad&\displaystyle \vartheta(\overline{p}_T - p_{T,\rm min})\times\\
   \\
   3):\quad&\displaystyle 
     \vartheta\left(\frac12\ln\left(\frac{Me^y-\overline{p}_Te^\eta}{Me^{-y}-\overline{p}_Te^{-
   \eta}}\right)-\eta_{\rm min}\right)\times\vartheta\left(\eta_{\rm max}-\frac12\ln\left(\frac{Me^y-\overline{p}_Te^\eta}{Me^{-y}-\overline{p}_Te^{-
   \eta}}\right)\right)\times\\
   \\
   4):\quad&\displaystyle \vartheta(\sqrt{|\mathbf{q}_T|^2+\overline{p}_T^2-2 |\mathbf{q}_T|\overline{p}_T\cos\phi} - p_{T,\rm min})\,,
   \end{array}

where in the last step we have replaced :math:`|\mathbf{p}_T|` with
:math:`\overline{p}_T` defined Eq. :eq:`eq:overpT`. Now
the question is identifying the integration domain defined by
:math:`\Phi(p,q-p)` on the :math:`(\eta,\cos\phi)`-plane. Since the
:math:`\theta`-functions in Eq. :eq:`eq:PSredDef4`
will be used inside a double nested integral over :math:`x=\cos\phi`
first and :math:`\eta` second, it is convenient to rewrite the function
:math:`\Phi(p,q-p)` in Eq. :eq:`eq:intdomain` as
follows:

.. math::
   :label: eq:almostfinal
   \begin{array}{rcl}
   \Phi(p,q-p) &=&\displaystyle \vartheta(\eta-\eta_{\rm min}) \times
   \vartheta(\eta_{\rm max}-\eta) \\
   \\
   &\times& \vartheta(x - f^{(2)}(\eta,
                   p_{T,\rm min})) \\
   \\
   &\times&\displaystyle
            \vartheta(f^{(3)}(\eta,\eta_{\rm min})-x) \times \vartheta(f^{(3)}(\eta,\eta_{\rm max})-x)\\
   \\
   &\times&\vartheta(f^{(4)}(\eta,
            p_{T,\rm min})-x)\,,
   \end{array}

with:

.. math::
   :label: eq:relevantfuncs
   \begin{array}{rcl}
   f^{(2)}(\eta, p_{T,\rm cut}) & = &\displaystyle \frac{2M p_{T,\rm min}\cosh(\eta-y) 
                       - Q^{2}}{2p_{T,\rm cut}
                       |\mathbf{q}_T|}\,, \\
   \\
   f^{(3)}(\eta,\eta_{\rm cut}) & = &\displaystyle \frac{M \cosh(\eta-y)}{|\mathbf{q}_T|
                       }-\frac{Q^{2} \left(\sinh(\eta
                                      -y)\coth(y-\eta_{\rm cut})+\cosh(\eta-y)\right)}{2|\mathbf{q}_T|  M}\,,\\
   \\
   f^{(4)}(\eta, p_{T,\rm cut}) & = &\displaystyle \frac{M \cosh(\eta-y)(Q^2 - 2
                       p_{T,\rm cut}^{2} + 2 |\mathbf{q}_T|^2)- Q^{2} \sqrt{M^{2} \sinh ^{2} (\eta-y) + p_{T,\rm min}^{2} }}{2 |\mathbf{q}_T| \left(M^{2} - p_{T,\rm min}^{2}\right)}\,.
   \end{array}

Considering that :math:`1\leq\cos\phi\leq 1`, the integration domain is
limited to this region. Therefore,
Eq. :eq:`eq:almostfinal` can be written in an even
more convenient way as:

.. math::
   :label: eq:final
   \begin{array}{rcl}
   \Phi(p,q-p) &=& \vartheta(\eta-\eta_{\rm min})\vartheta(\eta_{\rm
     max}-\eta)  \\
   \\
   &\times&\vartheta(x -
     \mbox{max}[f^{(2)}(\eta,p_{T,\rm min}),-1])\\
   \\
   &\times&\vartheta(\mbox{min}[f^{(3)}(\eta,\eta_{\rm min}),f^{(3)}(\eta,\eta_{\rm
            max}), f^{(4)}(\eta,p_{T,\rm
            min}),1]-x)
   \end{array}

such that a double integral over :math:`\eta` and :math:`x` would read:

.. math::
   \int_{-\infty}^{\infty}d\eta\int_{-1}^{1}dx\,\Phi(p,q-p)\dots =
   \int_{\eta_{\rm min}}^{\eta_{\rm
       max}}d\eta\,\vartheta(x_2(\eta)-x_1(\eta))\int_{x_1(\eta)}^{x_2(\eta)}dx\dots\,.

with:

.. math:: x_1(\eta) = \mbox{max}[f^{(2)}(\eta,p_{T,\rm min}),-1]

and:

.. math::
   x_2(\eta) = \mbox{min}[f^{(3)}(\eta,\eta_{\rm min}),f^{(3)}(\eta,\eta_{\rm
            max}),f^{(4)}(\eta,p_{T,\rm
            min}),1]\,.

.. container::
   :name: fig:IntDomain

   .. figure:: ../../latex/src/../../latex/src/plots/IntDomain.pdf
   .. :alt: The red area indicates the integration domain of the
      numerator in of the phase-space reduction factor
      Eq. :eq:`eq:PSredDef4` for
      :math:`p_{T,\rm min}=20` GeV and
      :math:`-\eta_{\rm min}=\eta_{\rm max}=2.4` at :math:`Q=91` GeV,
      :math:`|\mathbf{q}_T|=10` GeV and :math:`y=1`.[fig:IntDomain]
      :name: fig:IntDomain
      :width: 80.0%

      The red area indicates the integration domain of the numerator in
      of the phase-space reduction factor
      Eq. :eq:`eq:PSredDef4` for
      :math:`p_{T,\rm min}=20` GeV and
      :math:`-\eta_{\rm min}=\eta_{\rm max}=2.4` at :math:`Q=91` GeV,
      :math:`|\mathbf{q}_T|=10` GeV and :math:`y=1`.[fig:IntDomain]

As an example, Fig. `1 <#fig:IntDomain>`__ shows the integration domain
of the numerator in of the phase-space reduction factor
Eq. :eq:`eq:PSredDef4` for :math:`p_{T,\rm min}=20`
GeV and :math:`-\eta_{\rm min}=\eta_{\rm max}=2.4` at :math:`Q=91` GeV,
:math:`|\mathbf{q}_T|=10` GeV and :math:`y=1`. The gray band corresponds
to the region :math:`1\leq\cos\phi\leq 1`. The :math:`\theta`-function
1) in Eq. :eq:`eq:intdomain` limits the region to the
vertical stip defined by :math:`\eta_{\rm min} < \eta < \eta_{\rm max}`
(black vertical lines), the :math:`\theta`-function 2) gives the red
lines, the :math:`\theta`-functions 3) the blue lines , and the
:math:`\theta`-function 4) the green lines.

Gathering all pieces, the final expression for the phase-space reduction
factor reads:

.. math::
   :label: eq:finalformula
     \mathcal{P}(Q,y,q_T)=\displaystyle \int_{\eta_{\rm
         min}}^{\eta_{\rm
         max}}d\eta\,\vartheta(x_2(\eta)-x_1(\eta))\left[F(x_2(\eta),\eta)-F(x_1(\eta) ,\eta)\right]

with:

.. math::
   :label: eq:integrandF
   \begin{array}{rcl}
   \displaystyle F(x ,\eta)&=&\displaystyle \frac{1}{4\pi}\frac{Q^2}{E_q^2-q_T^2}\Bigg\{\frac{q_T^2 x\sqrt{1-x^2}}{x^2 q_T^2-
     E_q^2}\\
   \\
   &-&\displaystyle \frac{E_q}{\sqrt{E_q^2-q_T^2}}\left[\tan^{-1}\left(\frac{q_T-
       xE_q}{\sqrt{E_q^2-q_T^2}\sqrt{1-x^2}}\right)-\displaystyle\tan^{-1}\left(\frac{q_T+
       xE_q}{\sqrt{E_q^2-q_T^2}\sqrt{1-x^2}}\right)\right]\Bigg\}\
   \end{array}

where we have defined :math:`E_q = M\cosh(\eta-y)` and
:math:`q_T=|\mathbf{q}_T|`.

Let us consider the case :math:`y=q_T=0`. For simplicity, we also take
:math:`\eta_{\rm min} = -\eta_{\rm max}`. In these conditions,
Eq. :eq:`eq:integrandF` reduces to:

.. math::
   :label: eq:F00
   F(x,\eta)=\displaystyle \frac{1}{4\pi}\frac{1}{\cosh^2(\eta)}\left[\tan^{-1}\left(\frac{
       x}{\sqrt{1-x^2}}\right)-\displaystyle\tan^{-1}\left(-\frac{
       x}{\sqrt{1-x^2}}\right)\right]\,.

\ As evident from Eq. :eq:`eq:relevantfuncs`, for
:math:`q_T=0` all functions :math:`f^{(2)}`, :math:`f^{(3)}`, and
:math:`f^{(4)}` diverge. The relevant question, though, is whether they
go to plus or minus infinity depending on the value of :math:`Q`. Of
course, :math:`q_T` will tend to zero positively so we find:

.. math::
   \begin{array}{l}
   \displaystyle f^{(2)}(\eta,p_{T,\rm min})\rightarrow \infty\times\mbox{sign}\left[2 p_{T,\rm min}\cosh(\eta)-Q\right]\,,\\
   \\
   \displaystyle f^{(3)}(\eta,\eta_{\rm min}= -\eta_{\rm max}) \rightarrow +\infty\\
   \\
   \displaystyle f^{(3)}(\eta,\eta_{\rm max}) \rightarrow -\infty\,,\\
   \\
   \displaystyle f^{(4)}(\eta,p_{T,\rm min}) \rightarrow \infty\times\mbox{sign}\left[\frac{\cosh(\eta)(Q^2 - 2
                       p_{T,\rm min}^{2})- Q\sqrt{Q^{2} \sinh^{2}(\eta) + p_{T,\rm min}^{2} }}{Q^{2} - p_{T,\rm min}^{2}}\right]\,, 
   \end{array}

Therefore, :math:`f^{(3)}` never actually contributes. In addition, for
the :math:`\theta`-function in
Eq. :eq:`eq:finalformula` to be different from
zero, we need :math:`f^{(2)}(\eta)\rightarrow -\infty` and
:math:`f^{(3)}(\eta)\rightarrow \infty`. These both translate into
:math:`Q\geq 2 p_{T,\rm min}\cosh(\eta)`. This inequality is satisfied
only if :math:`Q\geq 2p_{T,\rm min}` for:

.. math::
   :label: eq:etabardef
   -\overline{\eta}\leq\eta\leq \overline{\eta}\quad\mbox{with}\quad \overline{\eta} =\cosh^{-1}\left(\frac{Q}{2p_{T,\rm min}}\right)\,.

Therefore, the phase-space reduction factor eventually becomes:

.. math::
   \begin{array}{rcl}
     \mathcal{P}(Q,0,0)&=&\displaystyle\frac12\vartheta(Q- 2p_{T,\rm min})\displaystyle \int_{-\eta_{\rm
         max}}^{\eta_{\rm
         max}}\frac{d\eta}{\cosh^2\eta}\,\vartheta(\eta+\overline{\eta})
     \vartheta(\overline{\eta}-\eta)\\
   \\
   &=&\displaystyle\vartheta(Q- 2p_{T,\rm
       min})\tanh(\mbox{max}[\eta_{\rm max},\overline{\eta}])\,.
   \end{array}

This result can be written more explicitly as:

.. math::
   :label: eq:partcase
   \mathcal{P}(Q,0,0) = 
   \left\{
   \begin{array}{ll}
   0 & \quad Q< 2p_{T,\rm min}\,,\\
   \displaystyle \tanh(\overline{\eta})=\left(1+\frac{2p_{T,\rm min}}{Q}\right)\sqrt{1-\frac{4 p_{T,\rm min}}{Q+2p_{T,\rm min}}}& \quad 2p_{T,\rm min} \leq Q < 2p_{T,\rm min}\cosh\eta_{\rm max}\,,\\
   \tanh(\eta_{\rm max}) & \quad Q \geq 2p_{T,\rm min}\cosh\eta_{\rm max}\,.
   \end{array}
   \right.

This differs from Eq. (24) of Ref. (Scimemi and Vladimirov 2018).
Despite the three different regions coincide, the behavior of the
phase-space reduction factor for all regions but for
:math:`Q< 2p_{T,\rm min}` is different. In favour of our result there is
the fact that :math:`\mathcal{P}(Q,0,0)` in
Eq. :eq:`eq:partcase` is continuous at
:math:`Q = 2p_{T,\rm min}\cosh\eta_{\rm max}` while that of
Ref. (Scimemi and Vladimirov 2018) is not. In addition, when setting
:math:`p_{T,\rm min} = 0` and :math:`\eta_{\rm max}=\infty`, *i.e.* no
cuts, our result tends to :math:`\mathcal{P}(Q,0,0)=\tanh(\infty)=1`, as
it should. While the result in Eq. (24) of Ref. (Scimemi and Vladimirov
2018) actually diverges in this limit.

The integrand of Eq. :eq:`eq:finalformula`, due to
the behaviour of :math:`x_1` and :math:`x_2` as functions of
:math:`\eta`, a piecewise function. Therefore, its numerical integration
is problematic in that quadrature algorithms assume the integrand be
continuos over the integration range. The solution is to identify the
discontinuity points and integrate the function separately over the
resulting ranges. However, the complexity of the integration region
(*e.g.* see Fig. `1 <#fig:IntDomain>`__) makes the analytical
identification of the discontinuity points very hard to achieve.

.. container::
   :name: contracting-the-leptonic-tensor-with-g_perpmunu

   .. rubric:: Contracting the leptonic tensor with
      :math:`g_\perp^{\mu\nu}`
      :name: contracting-the-leptonic-tensor-with-g_perpmunu

The calculation done in the previous section holds when contracting the
leptonic tensor :math:`L_{\mu\nu}` with the metric tensor
:math:`g^{\mu\nu}` associated with the Lorentz structure of the hadronic
tensor. However, at small values of :math:`|\mathbf{q}_T|`, the
leading-power Lorentz structure that one needs to multiply the leptonic
tensor for is:

.. math:: g_\perp^{\mu\nu} = g^{\mu\nu}+z^\mu z^\nu-t^\mu t^\nu

\ where the vectors :math:`z^\mu` and :math:`t^\mu` in the Collins-Soper
(CS) frame are defined as:

.. math::
   :label: eq:auxvects
   \begin{array}{l}
   \displaystyle z^\mu = (\sinh y,\mathbf{0},\cosh y)\,,\\
   \\
   \displaystyle t^\mu = \frac{q^\mu}{Q}\,,
   \end{array}

and they are such that :math:`z^2=-1`, :math:`t^2=1` and :math:`zq = 0`.
if we use the on-shell-ness of :math:`p_1` and :math:`p_2`
(:math:`p_1^2=p_2^2=0`) and the momentum conservation
(:math:`p\equiv p_1`, :math:`p_2=q-p`), we find that
:math:`t^\mu t^\nu L_{\mu\nu}=0` and the quantity above reduces to:

.. math::
   :label: eq:LT
     L_\perp = g_\perp^{\mu\nu}L_{\mu\nu} = 4\left[\frac12q^2 + 2(zp)^2\right] = 2Q^2\left[1+4 \sinh^2(y-\eta)\frac{|\mathbf{p}_T|^2}{Q^2}\right]\,.

Therefore, we need to introduce this factor in both the numerator and
the denominator of Eq. :eq:`eq:PSredDef4`. Following
the same steps of the previous section, up to a factor :math:`2Q^2`,
this leads us to replace the integral in
Eq. :eq:`eq:firstintegral` with( [8]_):

.. math::
   :label: eq:firstintegralNew
     \frac12 \int_0^\infty|\mathbf{p}_T|\left[1+4\sinh^2(y-\eta)
     \frac{|\mathbf{p}_T|^2}{Q^2}\right]d|\mathbf{p}_T|\,\delta(Q^2-2p\cdot
     q) =\frac{2\overline{p}_T^2}{Q^2}+
     2\sinh^2(y-\eta)\frac{\overline{p}_T^4} {Q^4}\,.

We can still use
Eq. :eq:`eq:complicatedintegral` for the
first term in the r.h.s. of the equation above. For the second, instead
we need to use:

.. math::
   :label: eq:complicatedintegral2
   \begin{array}{rcl}
   \displaystyle\int \frac{dx}{(a\pm
     x)^4\sqrt{1-x^2}}&=&\displaystyle\frac{\sqrt{1-x^2}\left[(11a^2+4)x^2\pm
                          3 a(9a^2+1)x + (18a^4-5a^2+2)\right]}{6(a^2-1)^3(x\pm
     a)^3}\\
   \\
   &\pm&\displaystyle\frac{a(2a^2+3)}{2(a^2-1)^{7/2}}\tan^{-1}\left(\frac{1\pm
         ax}{\sqrt{a^2-1}\sqrt{1-x^2}}\right)\,,
   \end{array}

that is such that:

.. math::
   \int_{-1}^{1} \frac{dx}{(a\pm
     x)^4\sqrt{1-x^2}}=\frac{\pi a(2a^2+3)}{2(a^2-1)^{7/2}}\,.

In our particular case, the integrand we are considering is the second
term in the r.h.s. term of
Eq. :eq:`eq:firstintegralNew`:

.. math::
   :label: eq:integralmmm
   \begin{array}{l}
   \displaystyle
     \frac{2\sinh^2(y-\eta)}{Q^4}\int_{-1}^{1}d(\cos\phi)\overline{p}_T^4(\cos\phi)=\\
   \\
   \displaystyle\frac{Q^4}{8q_T^4} \sinh^2(y-\eta) \int_{-1}^{1} \frac{dx}{\sqrt{1-x^2}}\left[\frac{1}{(a+
         x)^4}+\frac1{(a-
         x)^4}\right]= \frac{\pi Q^4}{8q_T^4}\sinh^2(y-\eta) \frac{
     a(2a^2+3)}{(a^2-1)^{7/2}}\,.
   \end{array}

with:

.. math:: a =\frac{M}{q_T}\cosh(y-\eta)\,.

Now we need to integrate Eq. :eq:`eq:integralmmm`
over :math:`\eta`:

.. math::
   \frac{\pi Q^4}{8q_T^4} \int_{-\infty}^{\infty}d\eta\sinh^2(y-\eta)\frac{
       a(2a^2+3)}{(a^2-1)^{7/2}}\,.

If we make the following change of variable in the integral above:

.. math:: z=\frac{M}{q_T}\sinh(y-\eta)

such that:

.. math:: a^2 = z^2+\frac{M^2}{q_T^2}\quad\mbox{and}\quad dz = -ad\eta\,,

the integral above becomes:

.. math::
   \frac{\pi Q^4}{4M^2 q_T^2} \int_{-\infty}^{\infty}dz\frac{
     z^2\left(z^2+\frac{M^2}{q_T^2}+\frac32\right)}{\left(z^2+\frac{M^2}{q_T^2}-1\right)^{7/2}}=\frac{\pi}{6}\,.

Putting this result together with
Eq. :eq:`eq:remarkableintegral` and taking
into account the factor :math:`2Q^2` in Eq. :eq:`eq:LT`, we
find that:

.. math::
   :label: eq:normalisation
   \begin{array}{l}
     \displaystyle \int d^4p \delta(p^2) \delta((q-p)^2) \theta(p_{0})
     \theta(q_0-p_{0})L_\perp = \frac{4\pi}{3}Q^2\,.
   \end{array}

This result agrees with Eq. (2.38) of Ref. (Scimemi and Vladimirov
2018), up to a factor four. This provides the denominator of the
phase-space-reduction factor :math:`\mathcal{P}`. The structure of
:math:`\mathcal{P}` will be exactly like that in
Eq. :eq:`eq:finalformula`, the only thing we need
to do is to identify the correct function :math:`F(x,\eta)`. To this
end, we need to make the following replacement for the function
:math:`F` given in Eq. :eq:`eq:integrandF` with:

.. math::
   :label: eq:FGcombination
   F(x,\eta)\rightarrow \overline{F}(x,\eta) = \frac34 F(x,\eta)+ \frac14 G(x,\eta)\,,

where:

.. math::
   :label: eq:integrandG
   \begin{array}{rcl}
   \displaystyle G(x ,\eta)&=&\displaystyle 
                               \frac{1}{16\pi }\sinh^2(y-\eta)\frac{Q^4}{(E_q^2-q_T^2)^3}
                               \Bigg\{\sqrt{1-x^2}q_T\\
   \\
   &\times&\displaystyle\Bigg[\frac{(11E_q^2q_T^2+4q_T^4)x^2+
                          3 E_qq_T(9E_q^2+q_T^2)x + (18E_q^4-5E_q^2q_T^2+2q_T^4)}{(xq_T+
     E_q)^3}\\
   \\
   &+&\displaystyle\frac{(11E_q^2q_T^2+4q_T^4)x^2-
                          3 E_qq_T(9E_q^2+q_T^2)x + (18E_q^4-5E_q^2q_T^2+2q_T^4)}{(xq_T-
     E_q)^3}\Bigg]\\
   \\
   &-&\displaystyle\frac{6E_q (2E_q^2+3q_T^2)}{\sqrt{E_q^2-q_T^2}}\left[\tan^{-1}\left(\frac{q_T-
         xE_q}{\sqrt{E_q^2-q_T^2}\sqrt{1-x^2}}\right)-\tan^{-1}\left(\frac{q_T+
         xE_q}{\sqrt{E_q^2-q_T^2}\sqrt{1-x^2}}\right)\right]
   \Bigg\}\,.
   \end{array}

Finally, combining the functions :math:`F` and :math:`G` given in
Eqs. :eq:`eq:integrandF`
and :eq:`eq:integrandG`, respectively, according to
Eq. :eq:`eq:FGcombination` to obtain
:math:`\overline{F}`, and replacing :math:`F` with :math:`\overline{F}`
in Eq. :eq:`eq:finalformula` gives the
phase-space-reduction factor :math:`\mathcal{P}(Q,y,q_T)` when the
leptonic tensor :math:`L_{\mu\nu}` is contracted with the transverse
metrics :math:`g_\perp^{\mu\nu}`.

As a final check, it is interesting to compute :math:`\mathcal{P}` when
:math:`y=q_T=0`, again assuming
:math:`\eta_{\rm min} = -\eta_{\rm max}`. In this limit :math:`F`
reduces to the expression in Eq. :eq:`eq:F00`, while
:math:`G` becomes:

.. math::
   \begin{array}{rcl}
   \displaystyle G(x ,\eta)&=&\displaystyle 
                               \frac{3}{4\pi }\frac{\sinh^2(\eta)}{\cosh^4(\eta)}
                               \Bigg\{\left[\tan^{-1}\left(\frac{
         x}{\sqrt{1-x^2}}\right)-\tan^{-1}\left(-\frac{
         x}{\sqrt{1-x^2}}\right)\right]
   \Bigg\}\,,
   \end{array}

\ such that:

.. math::
   \begin{array}{rcl}
     \mathcal{P}(Q,0,0)&=&\displaystyle\frac38\vartheta(Q- 2p_{T,\rm min})\displaystyle \int_{-\eta_{\rm
         max}}^{\eta_{\rm
         max}}d\eta\left[\frac{\cosh^2(\eta)+\sinh^2(\eta)}{\cosh^4(\eta)}\right]\,\vartheta(\eta+\overline{\eta})
     \vartheta(\overline{\eta}-\eta)\\
   \\
   &=&\displaystyle\vartheta(Q- 2p_{T,\rm
       min})\tanh(\mbox{max}[\eta_{\rm max},\overline{\eta}])\left[1-\frac{1}{4\cosh^2(\mbox{max}[\eta_{\rm max},\overline{\eta}])}\right]\,,
   \end{array}

with :math:`\overline{\eta}` defined in
Eq. :eq:`eq:etabardef`. The relation above can be
written more explicitly as:

.. math::
   :label: eq:partcase2
   {\footnotesize
   \mathcal{P}(Q,0,0) = 
   \left\{
   \begin{array}{ll}
   0 & \quad Q< 2p_{T,\rm min}\,,\\
    \tanh(\overline{\eta})\left[1-\frac{1}{4\cosh^2(\overline{\eta})}\right]=\left(1-\frac{p_{T,\rm min}^2}{Q^2}\right)\sqrt{1-\frac{4 p_{T,\rm min}^2}{Q^2}}& \quad 2p_{T,\rm min} \leq Q < 2p_{T,\rm min}\cosh\eta_{\rm max}\,,\\
    \tanh(\eta_{\rm max})\left[1-\frac{1}{4\cosh^2(\eta_{\rm max})}\right] & \quad Q \geq 2p_{T,\rm min}\cosh\eta_{\rm max}\,.
   \end{array}
   \right.}

The reason why we kept :math:`\eta_{\rm min} \neq -\eta_{\rm max}` is
that in some cases it may be required to implement an asymmetric cut,
:math:`\eta_{\rm min}<\eta<\eta_{\rm max}`. This is the case, for
example, of the LHCb experiment that delivers data only in the forward
region (:math:`2 < \eta < 4.5`).

.. container::
   :name: fig:IntDomainAsy

   .. figure:: ../../latex/src/../../latex/src/plots/IntDomainAsy.pdf
   .. :alt: Same as Fig. `1 <#fig:IntDomain>`__ for the asymmetric
      rapidity cut :math:`2<\eta < 4.5` at
      :math:`y=3`.[fig:IntDomainAsy]
      :name: fig:IntDomainAsy
      :width: 80.0%

      Same as Fig. `1 <#fig:IntDomain>`__ for the asymmetric rapidity
      cut :math:`2<\eta < 4.5` at :math:`y=3`.[fig:IntDomainAsy]

As an example, Fig. `2 <#fig:IntDomainAsy>`__ shows the integration
domain of the phase-space reduction factor
Eq. :eq:`eq:PSredDef4` for :math:`p_{T,\rm min}=20`
GeV and :math:`2<\eta < 4.5` at :math:`Q=91` GeV,
:math:`|\mathbf{q}_T|=10` GeV and :math:`y=3`.

.. container::
   :name: parity-violating-contribution

   .. rubric:: Parity-violating contribution
      :name: parity-violating-contribution

In the presence of cuts on the final-state leptons and for invariant
masses around the :math:`Z` mass, parity-violating effects arise. As we
will show below, these effects integrate to zero when removing the
leptonic cuts. This contribution stems from interference of the
antisymmetric contributions to the lepton tensor, proportional to
:math:`p_1^{\mu}
p_2^{\nu}\epsilon_{\mu\nu\rho\sigma}`, and the hadronic tensor,
proportional to :math:`\epsilon_{\perp}^{\mu\nu}` defined as:

.. math:: \epsilon_{\perp}^{\mu\nu}\equiv \epsilon^{\mu\nu\rho\sigma}t_\rho z_\sigma\,,

where :math:`t^\mu` and :math:`z^\mu` are given in
Eq. :eq:`eq:auxvects`. Therefore, the contribution we
are after results from the contraction of the following Lorentz
structures:

.. math::
   L_{\rm PV}\equiv p_1^{\mu}
   p_2^{\nu}\epsilon_{\mu\nu\rho\sigma}\epsilon_{\perp}^{\rho\sigma}\,.

After some manipulation, one finds:

.. math::
   :label: eq:pvphasespace
   \begin{array}{rcl}
   \displaystyle L_{\rm PV}&=& \displaystyle
   p_1^{\mu}
   p_2^{\nu}\epsilon_{\mu\nu\rho\sigma}\epsilon^{\rho\sigma\alpha\beta}t_\alpha
   z_\beta= -2 p_1^{\mu}
   p_2^{\nu}\delta_\mu^\alpha \delta_\nu^\beta t_\alpha \\
   \\
   &=&\displaystyle -2 (p_1 t)(p_2 z) = -2 (pt)\left[(qz)-(pz)\right] = 2(pt)(pz)\\
   \\
   &=&\displaystyle \frac{2|\mathbf{p}_T|^2}{Q}\sinh(y-\eta)\left[M\cosh(y-\eta)-|\mathbf{q}_T|\cos\phi\right]\,,
   \end{array}

where I have defined :math:`p_1\equiv p` and used the equalities
:math:`p_2 = q - p`, and :math:`zq=0`. To compute the third line I have
used the explicit parameterisation of :math:`q` and :math:`p` given in
Eqs. :eq:`eq:qexplicit`
and :eq:`eq:pexplicit`, respectively. The presence of
:math:`\sinh(y-\eta)` in Eq. :eq:`eq:pvphasespace`
is such that integrating over the full range in the lepton rapidity
:math:`\eta` nullifies this contribution:

.. math::
   :label: eq:noPVcontr
   \int_{-\infty}^{\infty} d\eta\,L_{\rm PV} = 0\,.

Therefore, it turns out that, for observables inclusive in the lepton
phase space, the parity violating term does not give any contribution.
Conversely, the presence of cuts on the final-state leptons may prevent
Eq. :eq:`eq:noPVcontr` from being satisfied leaving a
residual contribution. In order to quantify this effect, we take the
same steps performed in the previous sections to integrate
:math:`L_{\rm PV}` over the fiducial region. As above, we start
integrating over the full range in :math:`|\mathbf{p}_T|` using the
on-shell-ness :math:`\delta`-function:

.. math::
   \begin{array}{c}
   \displaystyle \int_0^\infty d|\mathbf{p}_T||\mathbf{p}_T|\,L_{\rm PV}=\frac{\sinh(y-\eta)}{Q}\left[M\cosh(y-\eta)-|\mathbf{q}_T|\cos\phi\right]
   \int_0^\infty d|\mathbf{p}_T||\mathbf{p}_T|^3\delta(Q^2-2pq)\\
   \\
   \displaystyle = \frac{\overline{p}_T^4}{Q^3}\sinh(y-\eta)\left[M\cosh(y-\eta)-|\mathbf{q}_T|\cos\phi\right]\,,
   \end{array}

with :math:`\overline{p}_T` defined in
Eq. :eq:`eq:overpT`. Now we compute the indefinte
integral over :math:`\cos\phi`. To do so, we need to use
Eq. :eq:`eq:intoverphi` along with the equality:

.. math::
   :label: eq:complicatedintegral3
   \int \frac{dx}{(a\pm
     x)^3\sqrt{1-x^2}}=\frac{\sqrt{1-x^2}\left[3ax\pm(4a^2-1) \right]}{2(a^2-1)^2(x\pm
     a)^2}\pm\displaystyle\frac{(2a^2+1)}{2(a^2-1)^{5/2}}\tan^{-1}\left(\frac{1\pm
         ax}{\sqrt{a^2-1}\sqrt{1-x^2}}\right)\,.

This allows us to compute the integral( [9]_):

.. math::
   :label: eq:Hdefinition
   \begin{array}{rcl}
   \displaystyle H(x,\eta)&=&\displaystyle\left(\frac {4\pi Q^2}{3}\right)^{-1}\frac{\sinh(y-\eta)}{Q^3}\left[M\cosh(y-\eta) \int
     d(\cos\phi)\overline{p}_T^4(\cos\phi)- |\mathbf{q}_T|
     \int d(\cos\phi) \cos\phi\,\overline{p}_T^4(\cos\phi)\right]\\
   \\
   &=&\displaystyle\frac{3Q^3 \sinh(y-\eta)}{64\pi |\mathbf{q}_T|^4}\left[M\cosh(y-\eta) \int
     \frac{d(\cos\phi)}{\left(a-\cos\phi\right)^4}-|\mathbf{q}_T|\int
     \frac{\cos\phi \,d(\cos\phi)}{\left(a-\cos\phi\right)^4}\right]\\
   \\
   &=&\displaystyle\frac{3Q^3 \sinh(y-\eta)}{64\pi q_T^3} \int
     \frac{dx}{\sqrt{1-x^2}}\left[\frac{1}{(a-x)^3}+\frac{1}{(a+x)^3}\right]\\
   \\
   &=&\displaystyle\frac{3Q^3 \sinh(y-\eta)}{128\pi q_T^3}\Bigg\{\frac{\sqrt{1-x^2}}{(a^2-1)^2}\left[\frac{3ax-(4a^2-1) }{(x-
     a)^2}+\frac{3ax+(4a^2-1) }{(x+
     a)^2}\right]\\
   \\
   &-&\displaystyle\frac{(2a^2+1)}{(a^2-1)^{5/2}}\left[\tan^{-1}\left(\frac{1-
         ax}{\sqrt{a^2-1}\sqrt{1-x^2}}\right) -\tan^{-1}\left(\frac{1+
         ax}{\sqrt{a^2-1}\sqrt{1-x^2}}\right)\right]\Bigg\}
   \end{array}

with :math:`E_q=M\cosh\left(\eta - y\right)`,
:math:`q_T=|\mathbf{q}_T|`, and:

.. math:: a = \frac{E_q}{q_T}\,,

so that:

.. math::
   :label: eq:Hfinal
   \begin{array}{rcl}
   \displaystyle H(x,\eta)&=&\displaystyle\frac{3Q^3 \sinh(y-\eta)}{128\pi (E_q^2-q_T^2)^{2}}\Bigg\{\sqrt{1-x^2}q_T\left[\frac{3E_qq_Tx-(4E_q^2-q_T^2) }{(xq_T-
     E_q)^2}+\frac{3E_qq_Tx+(4E_q^2-q_T^2) }{(q_Tx+
     E_q)^2}\right]\\
   \\
   &-&\displaystyle\frac{(2E_q^2+q_T^2)}{\sqrt{E_q^2-q_T^2}}\left[\tan^{-1}\left(\frac{q_T-
         xE_q}{\sqrt{E_q^2-q_T^2}\sqrt{1-x^2}}\right) -\tan^{-1}\left(\frac{q_T+
         xE_q}{\sqrt{E_q^2-q_T^2}\sqrt{1-x^2}}\right)\right]\Bigg\}\,.
   \end{array}

Finally, using the definition of :math:`H` in
Eq. :eq:`eq:Hfinal`, one can perform the integral over
the fiducial phase space as discussed above in
Eq. :eq:`eq:finalformula`:

.. math::
   :label: eq:finalformulaH
     \mathcal{P}_{\rm PV}(Q,y,q_T)=\displaystyle \int_{\eta_{\rm
         min}}^{\eta_{\rm
         max}}d\eta\,\vartheta(x_2(\eta)-x_1(\eta))\left[H(x_2(\eta),\eta)-H(x_1(\eta) ,\eta)\right]\,.

This allows one to estimate the impact of the parity-violating
contribution to the phase-space reduction factor.

In order to quantify numerically the impact of
Eq. :eq:`eq:finalformulaH`,
Fig. `3 <#fig:PhaseSpaceRedFactor>`__ displays the size
:math:`\mathcal{P}_{\rm PV}` relative to the parity-conserving
phase-space reduction factor as a function of :math:`y` for three
different values of :math:`q_T` at :math:`Q=M_Z` and for the following
lepton cuts: :math:`p_{T,\ell}>20` GeV and
:math:`-2.4 < \eta_\ell < 2.4`.

.. container::
   :name: fig:PhaseSpaceRedFactor

   .. figure:: ../../latex/src/../../latex/src/plots/PhaseSpaceRedFactor.pdf
   .. :alt: Ratio between the parity-violating phase-space reduction
      factor :math:`\mathcal{P}_{\rm PV}` in
      Eq. :eq:`eq:finalformulaH` and the
      respective parity-conserving factor as a function of the :math:`Z`
      rapidity :math:`y` at :math:`Q=M_Z` and for three different values
      of :math:`q_T`, with lepton cuts equal to
      :math:`p_{T,\ell}>20` GeV and
      :math:`-2.4 < \eta_\ell < 2.4`.[fig:PhaseSpaceRedFactor]
      :name: fig:PhaseSpaceRedFactor
      :width: 80.0%

      Ratio between the parity-violating phase-space reduction factor
      :math:`\mathcal{P}_{\rm PV}` in
      Eq. :eq:`eq:finalformulaH` and the
      respective parity-conserving factor as a function of the :math:`Z`
      rapidity :math:`y` at :math:`Q=M_Z` and for three different values
      of :math:`q_T`, with lepton cuts equal to
      :math:`p_{T,\ell}>20` GeV and
      :math:`-2.4 < \eta_\ell < 2.4`.[fig:PhaseSpaceRedFactor]

It turns out that the size of :math:`\mathcal{P}_{\rm PV}` relative to
:math:`\mathcal{P}` is never larger than :math:`2\times10^{-6}`. In
addition, the rapid oscillations with :math:`y` contribute to suppress
even more the integral over realistic bins in :math:`y`. One can thus
conclude that, for realistic kinematic configurations, the impact of
parity violating effects is completely negligible.

We finally notice that also following contraction enters the game in the
presence of cuts on the leptonic final state:

.. math:: L_\phi=(z^\mu t^\nu+z^\nu t^\mu)L_{\mu\nu}\,.

\ We find that:

.. math:: L_\phi=8(zp)\left[Q-2(tp)\right]\,.

Using Eqs. :eq:`eq:pexplicit`
and :eq:`eq:auxvects`, we have that:

.. math:: (zp) = p_T(\cosh\eta\sinh y-\sinh\eta\cosh y) = p_T\sinh(y-\eta)\,,

and:

.. math::
   (tp) = \frac1{Q}\left[Mp_T\cosh(y-\eta)-\mathbf{q}_T\cdot
     \mathbf{p}_T\right] = \frac{p_T}{Q}\left[M\cosh(y-\eta)-q_T\cos\phi\right]\,.

Therefore:

.. math::
   L_\phi=16 \frac{p_T^2}{Q}\sinh(y-\eta)\left[\frac{Q^2}{2 p_T
     }-M\cosh(y-\eta)+q_T\cos\phi\right]=16\left[\frac12p_TQ\sinh(y-\eta)-L_{\rm PV}\right]\,.

Due to the presence of the overall factor :math:`\sinh(y-\eta)`, also
this contribution is expected to be subdominant.

Differential cross section in the leptonic variables
====================================================

The calculation of the phase-space reduction factor carried out in the
previous section can be used to express the Drell-Yan cross section in
Eq. :eq:`eq:crosssection` as differential in the
kinematic variables of the single leptons. Loosely speaking, this
amounts to removing the integral sign in the numerator in
Eq. :eq:`eq:PSredDef` but taking into account kinematic
constraints. Using the transverse metric tensor
:math:`g_\perp^{\mu\nu}`, one finds:

.. math:: d\mathcal{P} = \frac{\displaystyle d^4p_1 d^4p_2 \,\delta(p_1^2) \delta(p_2^2)\theta(p_{1,0}) \theta(p_{2,0})\delta^{(4)}(p_1+p_2-q) L_\perp(p_1,p_2)}{\displaystyle \int d^4p_1 d^4p_2\, \delta(p_1^2) \delta(p_2^2) \theta(p_{1,0}) \theta(p_{2,0})\delta^{(4)}(p_1+p_2-q) L_\perp(p_1,p_2)}\,,

From Eq. :eq:`eq:normalisation`, we know the value
of the denominator. In the numerator, we can make use of the
momentum-conservation and one of the on-shell-ness
:math:`\delta`-functions. Using the r.h.s. of
Eq. :eq:`eq:firstintegralNew`, but integrating
over :math:`\phi` rather than :math:`|\mathbf{p}_T|`, leads to:

.. math::
   \frac{d\mathcal{P}}{d|\mathbf{p}_T|d\eta} = \frac {3 |\mathbf{p}_T|}{4\pi}\left[1+4 \sinh^2(y-\eta)\frac{|\mathbf{p}_T|^2}{Q^2}\right]
    \int_0^{2\pi} d\phi\,\delta(Q^2-2 |\mathbf{p}_T|\left[M\cosh\left(\eta - y\right)-|\mathbf{q}_T|\cos\phi\right])\,,

where we have used
Eqs. :eq:`eq:phasespacemeasure`-(`[eq:integralyeah!] <#eq:integralyeahux21>`__)
and Eq. :eq:`eq:deltaargument`. Finally, to
perform the integral over :math:`\phi` we use
Eq. :eq:`eq:intoverphi` to get:

.. math::
   :label: eq:difflepphasespace
   \frac{d\mathcal{P}}{d|\mathbf{p}_T|d\eta} = \frac {3|\mathbf{p}_T|}{2\pi Q^2}
    \frac{Q^2+4 |\mathbf{p}_T|^2 \sinh^2(\eta-y)}{\sqrt{4|\mathbf{p}_T|^2|\mathbf{q}_T|^2-(2|\mathbf{p}_T|M\cosh(\eta-y)-Q^2)^2}}\,.

Getting rid of the absolute value of the transverse vectors, this allows
one to get the Drell-Yan cross section differential in the leptonic
variables :math:`|\mathbf{p}_T|` and :math:`\eta`:

.. math::
   \frac{d\sigma}{dQ dy dq_T d\eta dp_T} =
                                                             \left[\frac
                                                             {3p_T^2}{\pi
                                                             Q^2 M}\frac{\left(\frac{Q^2}{4 p_T^2}-1\right)+\cosh^2(\eta-y)}{\sqrt{\frac{q_T^2}{M^2}-\left(\cosh(\eta-y)-\frac{Q^2}{2p_TM}\right)^2}}\right]\frac{d\sigma}{dQ dy dq_T}\,.

Due to the square root in the denominator, for fixed values of
:math:`Q`, :math:`q_T`, and :math:`y`, the expression above is defined
for values of :math:`p_T` and :math:`\eta` such that:

.. math::
   :label: eq:kinconstr
   \frac{Q^2}{2p_TM}-\frac{q_T }{M}< \cosh(\eta-y) < \frac{Q^2}{2p_TM}+\frac{q_T}{M}\,.

In order to focus on the :math:`p_T` dependence of the cross section,
one may want to integrate of the lepton rapidity :math:`\eta`. Using the
analogous of Eq. :eq:`eq:intoverphi` for
:math:`\cosh(\eta)`:

.. math::
   :label: eq:intovereta
   \int_{-\infty}^{\infty}d\eta\, f(\cos\eta) = \int_{1}^{\infty}\frac{dx}{\sqrt{x^2-1}}\left[f(x)+f(-x)\right]\,,

and taking into account the constraint in
Eq. :eq:`eq:kinconstr`, one finds:

.. math::
   :label: eq:LeptonicModulation
   \begin{array}{l}
   \displaystyle \frac{d\sigma}{dQ dy dq_T dp_T} =\frac{d\sigma}{dQ dy
     dq_T} \times\\
   \\
   \displaystyle\frac
                                                             {3p_T}{2\pi
                                                             Q^2}\int_{-\frac{q_T}{M}}^{\frac{q_T}{M}}
     \frac{dx}{\sqrt{\frac{q_T^2}{M^2}-x^2}}\left(\frac{\left(\frac{Q^2}{4
     p_T^2}-1\right)+\left(x+\frac{Q^2}{2p_TM}\right)^2}{\sqrt{\left(x+\frac{Q^2}{2p_TM}\right)^2-1}}+\frac{\left(\frac{Q^2}{4
     p_T^2}-1\right)+\left(x-\frac{Q^2}{2p_TM}\right)^2}{\sqrt{\left(x-\frac{Q^2}{2p_TM}\right)^2-1}}\right)\,.
   \end{array}

For fixed values of :math:`Q`, :math:`q_T`, and :math:`y`, the integral
above can be solved numerically and plotted as a function of
:math:`p_T`. The result is shown in
Fig. `4 <#fig:LeptonicModulation>`__.

.. container::
   :name: fig:LeptonicModulation

   .. figure:: ../../latex/src/../../latex/src/plots/LeptonicModulation.pdf
   .. :alt: Behaviour of the second line of
      Eq. :eq:`eq:LeptonicModulation` at
      :math:`Q=M_Z` and :math:`y=0` as a function of the lepton
      transverse momentum :math:`p_{T,\ell}`.[fig:LeptonicModulation]
      :name: fig:LeptonicModulation
      :width: 80.0%

      Behaviour of the second line of
      Eq. :eq:`eq:LeptonicModulation` at
      :math:`Q=M_Z` and :math:`y=0` as a function of the lepton
      transverse momentum :math:`p_{T,\ell}`.[fig:LeptonicModulation]

Convolution in transverse-momentum space
========================================

The “transverse-momentum” version of
Eq. :eq:`eq:crosssection` reads:

.. math::
   :label: eq:crosssectionkt
   \begin{array}{rcl}
     \displaystyle \frac{d\sigma}{dQ dy dq_T} &=&\displaystyle 
     \frac{16\pi\alpha^2q_T}{9 Q^3} H(Q,\mu) \sum_q C_q(Q)\\
   \\
   &\times&\displaystyle 
     \int d^2\mathbf{k}_{T,1}d^2\mathbf{k}_{T,2}
     \overline{F}_q(x_1,\mathbf{k}_{T,1};\mu,\zeta)
     \overline{F}_{\bar{q}}(x_2,\mathbf{k}_{T,2};\mu,\zeta)\delta^{(2)}(\mathbf{q}_{T}-\mathbf{k}_{T,1}-\mathbf{k}_{T,2})\,,
   \end{array}

\ where, with abuse of notation, we used the same symbol for the
distributions :math:`F` both in :math:`k_T` and :math:`b_T` space. We
can make use of the :math:`\delta`-function to reduce the expression
above to:

.. math::
   :label: eq:crosssectionkt2
   \frac{d\sigma}{dQ dy dq_T} =\frac{16\pi\alpha^2q_T}{9 Q^3} H(Q,\mu) \sum_q C_q(Q)\int d^2\mathbf{k}_{T}
     \overline{F}_q(x_1,\mathbf{k}_{T};\mu,\zeta)
     \overline{F}_{\bar{q}}(x_2,\mathbf{q}_{T}-\mathbf{k}_{T};\mu,\zeta)\,.

Since the distributions :math:`F` only depend on the absolute value of
their argument (be it :math:`k_T` or :math:`b_T`), the expression above
can be further reduced to:

.. math::
   :label: eq:crosssectionkt3
   \begin{array}{rcl}
   \displaystyle \frac{d\sigma}{dQ dy dq_T} &=&\displaystyle
                                                \frac{16\pi\alpha^2q_T}{9
                                                Q^3} H(Q,\mu) \sum_q
                                                C_q(Q)\\
   \\
   &\times&\displaystyle \int_0^\infty dk_{T}\,k_{T}\int_0^{2\pi}d\theta\,
     \overline{F}_q(x_1,k_T;\mu,\zeta)
     \overline{F}_{\bar{q}}(x_2,\sqrt{q_{T}^2+k_{T}^2
     -2q_{T}k_{T}\cos\theta};\mu,\zeta)\,,
   \end{array}

where, without loss of generality, we are assuming that the vector
:math:`\mathbf{q}_T` is directed along the :math:`x` axis. Knowing the
functions :math:`F` in :math:`k_T` space, the formula above can be
implemented numerically.

The :math:`A_0` coefficient
===========================

Using Eq. (5) of Ref. (Richter-Was and Was 2016), the azimuthal angle
:math:`\theta` of the (negative) lepton in the Collins-Soper frame can
be related to the momenta on the vector boson and of one of the leptons
in the laboratory frame as:

.. math::
   :label: eq:RitterSport
   \cos\theta(\eta,p_T) = \mbox{sign}(y)\frac{2p_T}{Q}\sinh(y-\eta).

\ This can be related to the coefficient :math:`A_0` as (Gauld et al.
2017):

.. math:: A_0= 4 - 10\,\langle\cos^2\theta\rangle\,,

where:

.. math::
   \begin{array}{rcl}
   \displaystyle \langle\cos^2\theta\rangle &\equiv& \displaystyle  \left[\int dp_Td\eta
     \frac{d\sigma}{dQdydq_T d\eta dp_T}\right]^{-1}\int
     dp_Td\eta\left[\cos\theta(\eta, p_T)\right]^2 \frac{d\sigma}{dQdydq_T d\eta dp_T}\,.
   \end{array}

If the integration in both numerator and denominator in the r.h.s. of
the quantity above is performed over the full phase space, the result at
small :math:`q_T` is:

.. math::
   \langle\cos^2\theta\rangle =\int
     dp_Td\eta\left[\cos\theta(\eta, p_T)\right]^2 \frac{d\mathcal{P}}{dp_Td\eta}\,

with the element of phase space :math:`d\mathcal{P}` is given in
Eq. :eq:`eq:difflepphasespace`. Therefore, the
“QCD dynamics” totally cancels in the ration and one is left with a
relatively interesting kinematic factor.

Implementing the :math:`\phi*` distribution
===========================================

A physical observable that is particularly useful to probe the intrinsic
transverse dynamics of hadrons is the so-called :math:`\phi^*`
distribution in Drell-Yan production, defined as (Banfi et al. 2011):

.. math:: \phi^*=\tan\left(\frac{\pi-\Delta\phi}{2}\right)\sin\theta^*=\tan\left(\frac{\pi-\Delta\phi}{2}\right)\sqrt{1-\tanh^2\left(\frac{\Delta\eta}{2}\right)}\,,

where :math:`\Delta\phi = \phi_1-\phi_2` and
:math:`\Delta\eta = \eta_1-\eta_2` are the separation in azimuthal angle
and pseudo-rapidity of the two outgoing leptons. The reason why
:math:`\phi^*` is interesting is that experimental measurements of this
quantity reach an extremely high accuracy. It is therefore useful to
express the Drell-Yan cross section differential in :math:`q_T` into a
cross section differential in :math:`\phi^*`. This is achieved through
the chain rule:

.. math:: \frac{d\sigma}{dQdydq_T d\eta dp_T} = \frac{d\phi^*}{dq_T}\frac{d\sigma}{dQdyd\phi^*d\eta dp_T}\,.

Of course, this requires the knowledge of :math:`\phi^*` as a function
of :math:`q_T` which is what we will work out in the following. The
starting point is the four-momentum conservation :math:`q=p_1+p_2`.
Taking the absolute value of this relation and using the
parameterisations in Eq. :eq:`eq:qexplicit` for
:math:`q` and in Eq. :eq:`eq:pexplicit` for both
:math:`p_1` and :math:`p_2` immediately gives:

.. math::
   :label: eq:coshdeta
   \cosh\Delta\eta = 1+\frac{Q^2}{2p_{T1}p_{T2}}\,.

In addition, equating the absolute values of the transverse component
(:math:`|\mathbf{q}_T|^2=|\mathbf{p}_{T1}+\mathbf{p}_{T2}|^2`) leads to:

.. math::
   :label: eq:cosdphi
   \cos\Delta\phi =\frac{p_{T1}^2+p_{T2}^2-q_T^2}{2p_{T1}p_{T2}}\,.

Using Eq. :eq:`eq:coshdeta` one has:

.. math::
   \tanh\left(\frac{\Delta\eta}{2}\right) =
   \sqrt{\frac{\cosh\Delta\eta-1}{\cosh\Delta\eta+1}} = \sqrt{\frac{Q^2}{Q^2+4p_{T1}p_{T2}}}\,,

so that:

.. math::
   \sin\theta^*=\sqrt{1-\tanh^2\left(\frac{\Delta\eta}{2}\right)}
   =\sqrt{\frac{4p_{T1}p_{T2}}{Q^2 + 4p_{T1}p_{T2}}}\,.

On the other hand, using Eq. :eq:`eq:cosdphi` one finds:

.. math:: \tan\left(\frac{\pi-\Delta\phi}{2}\right) = \sqrt{\frac{1-\cos(\pi-\Delta\phi)}{1+\cos(\pi-\Delta\phi)}}=\sqrt{\frac{1+\cos\Delta\phi}{1-\cos\Delta\phi}}=\sqrt{\frac{p_{T1}^2+p_{T2}^2-q_T^2+2p_{T1}p_{T2}}{p_{T1}^2+p_{T2}^2-q_T^2-2p_{T1}p_{T2}}}\,.

Now let us define :math:`p_{T1}=p_T` and :math:`\eta=\eta_1` and trade
:math:`p_{T2}` for :math:`\cosh(y-\eta)` using the relation:

.. math:: p_{T2}^2=M^2+p_T^2-2p_TM\cosh(y-\eta)

with :math:`M=\sqrt{Q^2+q_T^2}`. This gives:

.. math:: \sin\theta^*=\sqrt{\frac{4p_{T}\sqrt{M^2+p_T^2-2p_TM\cosh(y-\eta)}}{Q^2 +4p_{T}\sqrt{M^2+p_T^2-2p_TM\cosh(y-\eta)}}}\,.

and:

.. math:: \tan\left(\frac{\pi-\Delta\phi}{2}\right) = \sqrt{\frac{Q^2+2p_T\left[p_T-M\cosh(y-\eta)+\sqrt{M^2+p_T^2-2p_TM\cosh(y-\eta)}\right]}{Q^2+2p_T\left[p_T-M\cosh(y-\eta)-\sqrt{M^2+p_T^2-2p_TM\cosh(y-\eta)}\right]}}\,.

Gathering all pieces together finally gives:

.. math::
   :label: eq:phistarmaster
   \begin{array}{rcl}
   \phi^*&=&\displaystyle
             \sqrt{\frac{Q^2+2p_T\left[p_T-M\cosh(y-\eta)+\sqrt{M^2+p_T^2-2p_TM\cosh(y-\eta)}\right]}{Q^2+2p_T\left[p_T-M\cosh(y-\eta)-\sqrt{M^2+p_T^2-2p_TM\cosh(y-\eta)}\right]}}\\
   \\
   &\times& \displaystyle\sqrt{\frac{4p_{T}\sqrt{M^2+p_T^2-2p_TM\cosh(y-\eta)}}{Q^2 +4p_{T}\sqrt{M^2+p_T^2-2p_TM\cosh(y-\eta)}}}\,.
   \end{array}

One can show that the boost required to go from the laboratory frame, in
which the vectors :math:`q` and :math:`p` are given by
Eqs. :eq:`eq:qexplicit`
and :eq:`eq:pexplicit`, respectively, with :math:`p_T`
given by Eq. :eq:`eq:overpT`, to the Collins-Soper frame,
in which:

.. math::
   :label: eq:CSdecomp
   q = Q(1,0,0,0)\quad\mbox{and}\quad p = \frac{Q}2(1,\sin\theta\cos\phi, \sin\theta\sin\phi,\cos\theta)\,,

is given by:

.. math::
   :label: eq:boostLabToCS
   \Lambda_{{\rm Lab}\rightarrow {\rm CS}} = 
   \begin{pmatrix}
   \frac{M}{Q}\cosh y & -\frac{q_T}{Q} & 0 & -\frac{M}{Q}\sinh y\\
   -\frac{q_T}{Q}\cosh y & \frac{M}{Q} & 0 & \frac{q_T}{Q}\sinh y\\
   0 & 0 & 1 & 0 \\
   -\sinh y & 0 & 0 & \cosh y
   \end{pmatrix}\,,

with :math:`M=\sqrt{Q^2+q_T^2}`. Likewise, the boost to go from the
Collins-Soper to the laboratory frame is the inverse of the above
transformation and reads:

.. math::
   :label: eq:boostCSToLab
   \Lambda_{ {\rm CS}\rightarrow{\rm Lab}} = \Lambda_{{\rm Lab}\rightarrow {\rm CS}}^{-1}=
   \begin{pmatrix}
   \frac{M}{Q}\cosh y & \frac{q_T}{Q}\cosh y & 0 & \sinh y\\
   \frac{q_T}{Q} & \frac{M}{Q} & 0 & 0\\
   0 & 0 & 1 & 0 \\
   \frac{M}{Q}\sinh y & \frac{q_T}{Q}\sinh y & 0 & \cosh y
   \end{pmatrix}\,.

Notice that rotating the vector in
Eq. :eq:`eq:pexplicit` using the transformation in
Eq. :eq:`eq:boostLabToCS` and equating the
:math:`z` component to that of the vector :math:`p` in
Eq. :eq:`eq:CSdecomp` immediately gives
Eq. :eq:`eq:RitterSport`. This immediately allows
one to write:

.. math:: \sin\theta = \frac{\sqrt{Q^2-4p_T\sinh^2(\eta-y)}}{Q}

**References**

**References**

.. container:: references csl-bib-body hanging-indent
   :name: refs

   .. container:: csl-entry
      :name: ref-Banfi:2010cf

      Banfi, A., S. Redford, M. Vesterinen, P. Waller, and T. R. Wyatt.
      2011. “Optimisation of variables for studying dilepton transverse
      momentum distributions at hadron colliders.” *Eur. Phys. J. C* 71:
      1600. https://doi.org/10.1140/epjc/s10052-011-1600-y.

   .. container:: csl-entry
      :name: ref-Collins:2011zzd

      Collins, John. 2013. *Foundations of perturbative QCD*. Vol. 32.
      Cambridge University Press.

   .. container:: csl-entry
      :name: ref-Gauld:2017tww

      Gauld, R., A. Gehrmann-De Ridder, T. Gehrmann, E. W. N. Glover,
      and A. Huss. 2017. “Precise predictions for the angular
      coefficients in Z-boson production at the LHC.” *JHEP* 11: 003.
      https://doi.org/10.1007/JHEP11(2017)003.

   .. container:: csl-entry
      :name: ref-Ogata:quadrature

      Ogata, H. n.d. “A Numerical Integration Formula Based on the
      Bessel Functions.”
      `http://www.kurims.kyoto-u.ac.jp/$\sim$okamoto/paper/Publ\_RIMS\_DE/41-4-40.pdf <http://www.kurims.kyoto-u.ac.jp/$\sim$okamoto/paper/Publ\_RIMS\_DE/41-4-40.pdf>`__.

   .. container:: csl-entry
      :name: ref-Richter-Was:2016mal

      Richter-Was, E., and Z. Was. 2016. “Separating electroweak and
      strong interactions in Drell:raw-latex:`\textendash{}`Yan
      processes at LHC: leptons angular distributions and reference
      frames.” *Eur. Phys. J. C* 76 (8): 473.
      https://doi.org/10.1140/epjc/s10052-016-4319-y.

   .. container:: csl-entry
      :name: ref-Scimemi:2017etj

      Scimemi, Ignazio, and Alexey Vladimirov. 2018. “Analysis of vector
      boson production within TMD factorization.” *Eur. Phys. J. C* 78
      (2): 89. https://doi.org/10.1140/epjc/s10052-018-5557-y.

.. [1]
   Note that in Eq. :eq:`eq:crosssection` the gluon
   TMD PDF :math:`\overline{F}_g` is not involved. If also the gluon TMD
   PDF was involved, it would evolve by means of a different evolution
   factor :math:`R_g`.

.. [2]
   The superscript 0 in :math:`z_n^{(0)}` and :math:`w_n^{(0)}`
   indicates that here we are performing a Hankel tranform that involves
   the Bessel function of degree zero :math:`J_0`. This is useful in
   view of the next section in which the integration over :math:`q_T`
   will give rise to a similar Hankel transform with :math:`J_0`
   replaced by :math:`J_1`. Also in that case the Ogata quadrature
   algorithm can be applied but coordinates and weights will be
   different.

.. [3]
   In fact, :math:`\mathcal{P}` also depends on the invariant mass
   :math:`Q` and the rapidity :math:`y` of the lepton pair that also
   need to be integrated over.

.. [4]
   For the moment we ignore the complication introduced by the presence
   of cuts on the final state discussed in
   Sect. `2.1.1 <#sec:kincuts>`__. We will come back on this issue at
   the end of the section.

.. [5]
   The same procedure applies to the tensor
   :math:`\overline{W}_{n\tau\alpha}` defined in
   Eq. :eq:`eq:xFtens`.

.. [6]
   In the general case in which :math:`\mathbf{q}_T` forms an angle
   :math:`\beta` with the :math:`x` axis, the scalar product would
   result in :math:`|\mathbf{p}_T||\mathbf{q}_T|\cos(\phi-\beta)`.
   However, the angle :math:`\beta` could always be reabsorbed in a
   redefinition of the integration angle :math:`\phi` in
   Eq. (`[eq:integralyeah!] <#eq:integralyeahux21>`__).

.. [7]
   Notice that the :math:`\vartheta`-function has no effect. I have
   verified it numerically but I cannot see it analytically.

.. [8]
   We removed the :math:`\theta`-function as we know it does not have
   any effect.

.. [9]
   The factor :math:`\left(\frac {4\pi Q^2}{3}\right)^{-1}` in
   Eq. :eq:`eq:Hdefinition` corresponds to the full
   phase-space in integral of :math:`L_\perp` in
   Eq. :eq:`eq:normalisation` that provides the
   natural normalisation.
