========================================
SIDIS cross section in TMD factorisation
========================================

:Author: Valerio Bertone

.. contents::
   :depth: 3
..

Structure of the observable
===========================

In this document we report the relevant formulas for the computation of
semi-inclusive deep-inelastic scattering (SIDIS) multiplicities under
the assumption that the (negative) virtuality of the :math:`Q^2` of the
exchanged vector boson is much smaller than the :math:`Z` mass. This
allows us to neglect weak contributions and write the cross section in
TMD factorisation as:

.. math::
   :label: eq:sidisxsec
     \frac{d\sigma}{dxdQdz d q_T} = \frac{4\pi \alpha^2q_T}{z x Q^3}Y_+ H(Q,\mu) \sum_q e_q^2
     \int_0^\infty db \,b J_0\left(bq_T\right)\overline{F}_q(x,b;\mu,\zeta_1) \overline{D}_{q}(z,b;\mu,\zeta_2)\,,

with :math:`\zeta_1\zeta_2=Q^4` and:

.. math:: Y_+=1+(1-y)^2=1+\left(1-\frac{Q^2}{xs}\right)^2\,,

where :math:`s` is the squared center of mass energy. The single TMDs
are evolved and matched onto the respective collinear functions as
usual:

.. math:: \overline{F}_q(x,b;\mu,\zeta) =xF_q(x,b;\mu,\zeta) = R_q(\mu_0,\zeta_0\rightarrow \mu,\zeta;b) \sum_{j}\int_x^1dy\,\mathcal{C}_{qj}(y;\mu_0,\zeta_0)\left[\frac{x}{y}f_j\left(\frac{x}{y},\mu_0\right)\right]\,,

and:

.. math:: \overline{D}_{i}(z,b;\mu,\zeta) =z^3D_{q}(z,b;\mu,\zeta) = R_q(\mu_0,\zeta_0\rightarrow \mu,\zeta;b) \sum_{j}\int_z^1dy\,\left[y^2\mathbb{C}_{qj}(y;\mu_0,\zeta_0)\right]\left[\frac{z}{y}d_j\left(\frac{z}{y},\mu_0\right)\right]\,.

Notice that here we limit to the case :math:`Q\ll M_Z` such that we can
neglect the contribution of the :math:`Z` boson and thus the electroweak
couplings are given by the squared electric charges.

As usual, low-:math:`q_T` non-perturbative corrections are taken into
account by introducing the monotonic function :math:`b_*(b)` that
behaves as:

.. math::
   \lim_{b\rightarrow 0}
     b_*(b) = b_{\rm min}\quad\mbox{and}\quad\lim_{b\rightarrow \infty}
     b_*(b) = b_{\rm max}\,.

\ This allows us to replace the TMDs in
Eq. :eq:`eq:sidisxsec` with their “regularised”
version:

.. math::
   :label: eq:fandDNP
   \begin{array}{rcl}
     \overline{F}_q(x,b;\mu,\zeta) &\rightarrow&
     \overline{F}_q(x,b_*(b);\mu,\zeta) f_{\rm NP}(x,b,\zeta)\,,\\
   \\
     \overline{D}_q(z,b;\mu,\zeta) &\rightarrow&
     \overline{D}_q(z,b_*(b);\mu,\zeta) D_{\rm NP}(z,b,\zeta)\,,
   \end{array}

where we have introduced the non-perturbative functions
:math:`f_{\rm NP}` and :math:`D_{\rm NP}`. It is important to stress
that these functions further factorise as follows:

.. math::
   :label: NPfuncts
   \begin{array}{rcl}
   f_{\rm NP}(x,b,\zeta)&=&\displaystyle\widetilde{f}_{\rm NP}(x,b)\exp\left[g_K(b)\ln\left(\frac{\zeta}{Q_0^2}\right)\right]\,,\\
   \\
   D_{\rm NP}(z,b,\zeta) &=&\displaystyle\widetilde{D}_{\rm NP}(x,b)\exp\left[g_K(b)\ln\left(\frac{\zeta}{Q_0^2}\right)\right]\,.
   \end{array}

The common exponential function represents the non-perturbative
corrections to TMD evolution and the specific functional form is driven
by the solution of the Collins-Soper equation where :math:`Q_0` is some
initial scale. Finally the set of non-perturbative functions to be
determined from fits to data are :math:`\widetilde{f}_{\rm NP}`,
:math:`\widetilde{D}_{\rm NP}`, and :math:`g_K(b)`. It is worth noticing
that by definition

.. math:: f_{\rm NP}(x,b,\zeta) = \frac{\overline{F}_q(x,b;\mu,\zeta)}{\overline{F}_q(x,b_*(b);\mu,\zeta)} \,,\\

and similarly for :math:`D_{\rm NP}`. Therefore, one has a partial
handle on the :math:`b`-dependence of these functions from the region in
which :math:`b` is small enough to make both numerator and denominator
perturbatively computable. Making use of
Eq. :eq:`NPfuncts` and setting :math:`\zeta_1=\zeta_2=Q^2`
allows us to rewrite Eq. :eq:`eq:sidisxsec` as:

.. math::
   :label: eq:sidisxsec3
   \begin{array}{rcl}
   \displaystyle   \frac{d\sigma}{dxdQdz d q_{T}} &=&\displaystyle
   \frac{4\pi \alpha^2q_{T}}{xzQ^3}Y_+ H(Q,\mu) \sum_q e_q^2\\
   \\
   &\times&\displaystyle 
     \int_0^\infty db\, J_0\left(bq_T\right) \,b\overline{F}_q(x,b_*(b);\mu,Q^2) \overline{D}_q(z,b_*(b);\mu,Q^2) f_{\rm NP}(x,b,Q^2) D_{\rm NP}(z,b,Q^2)\,.
   \end{array}

The integral in the r.h.s. can be numerically computed using the Ogata
quadrature of zero-th degree (because :math:`J_0` enters the integral):

.. math::
   \begin{array}{rcl}
   \displaystyle   \frac{d\sigma}{dxdQdz d q_{T}} &\simeq&\displaystyle
   \frac{4\pi \alpha^2}{xzQ^3}Y_+ H(Q,\mu) \sum_q e_q^2\\
   \\
   &\times&\displaystyle
     \sum_{n=1}^{N}w_n^{(0)}\,\frac{\xi_n^{(0)}}{q_{T}}\overline{F}_q\left(x,b_*\left(\frac{\xi_n^{(0)}}{q_{T}}\right);\mu,Q^2\right)
            \overline{D}_q\left(z,b_*\left(\frac{\xi_n^{(0)}}{q_{T}}\right);\mu,Q^2\right)\\
   \\
   &\times& \displaystyle f_{\rm NP}\left(x, \frac{\xi_n^{(0)}}{q_{T}},Q^2\right) D_{\rm NP}\left(z, \frac{\xi_n^{(0)}}{q_{T}},Q^2\right)\,,
   \end{array}

where :math:`w_n^{(0)}` and :math:`\xi_n^{(0)}` are the Ogata weights
and coordinates, respectively, and the sum over :math:`n` is truncated
to the :math:`N`-th term that should be chosen in such a way to
guarantee a given target accuracy. The equation above can be
conveniently recasted as follows:

.. math::
   \displaystyle \frac{d\sigma}{dxdQdz d q_{T}} \simeq
     \sum_{n=1}^{N}w_n^{(0)}\,\frac{\xi_n^{(0)}}{q_{T}}S\left(x,z,\frac{\xi_n^{(0)}}{q_{T}};\mu,Q^2\right) f_{\rm NP}\left(x, \frac{\xi_n^{(0)}}{q_{T}},Q^2\right) D_{\rm NP}\left(z, \frac{\xi_n^{(0)}}{q_{T}},Q^2\right)\,,

where:

.. math::
   :label: eq:Sall
   \displaystyle S\left(x,z,b;\mu,Q^2\right)=
   \frac{4\pi \alpha^2}{xzQ^3}Y_+ H(Q,\mu) \sum_q e_q^2 \left[\overline{F}_q\left(x,b_*(b);\mu,Q^2\right)\right]\left[\overline{D}_q\left(z,b_*(b);\mu,Q^2\right)\right]\,.

Integrating over the final-state kinematic variables
====================================================

Experimental measurements of differential distributions for SIDIS
production are often delivered as integrated over finite regions of the
final-state kinematic phase space.

More specifically, the cross section is not integrated of the transverse
momentum of the vector boson, :math:`q_T`, but over the transverse
momentum of the outgoing hadron, :math:`p_{Th}`, that is connected to
the former through:

.. math:: p_{Th} = zq_T\,.

\ The integrated cross section then reads:

.. math::
   :label: eq:IntcrosssectionSIDIS
     \widetilde{\sigma}=\int_{Q_{\rm min}}^{Q_{\rm max}}dQ \int_{x_{\rm min}}^{x_{\rm max}}dx
     \int_{z_{\rm
         min}}^{z_{\rm max}}dz \int_{p_{Th,\rm min}/z}^{p_{Th,\rm max}/z}dq_{T}\left[\frac{d\sigma}{dxdQ
         dz dq_{T}} \right]\,.

One can exploit a property of the Bessel functions to compute the
indefinite integral in :math:`q_{T}` of the cross section. Specifically,
we now compute:

.. math:: K(x,z,Q,q_{T}) = \int dq_{T}\left[\frac{d\sigma}{dxdQdz dq_{T}}\right]\,.

This is easily done by using the following property of the Bessel
functions:

.. math:: \int dx\,x J_0(x) = xJ_1(x)\,,

that is equivalent to:

.. math:: \int dq_{T}\,q_{T} J_0\left(bq_T\right) = \frac{q_T}{b}J_1\left(bq_T\right)\,.

Therefore:

.. math::
   \begin{array}{rcl}
   \displaystyle   K(x,z,Q,q_{T}) &=&\displaystyle
   \frac{4\pi \alpha^2q_{T}}{xzQ^3}Y_+ H(Q,\mu) \sum_q e_q^2\\
   \\
   &\times&\displaystyle 
     \int_0^\infty db\, J_1\left(bq_T\right)
            \,\overline{F}_q(x,b_*(b);\mu,Q^2)
            \overline{D}_q(z,b_*(b);\mu,Q^2) f_{\rm NP}(x,b,Q^2) D_{\rm NP}(z,b,Q^2)\,.
   \end{array}

The integral can again be computed using the Ogata quadrature as:

.. math::
   :label: eq:primitivepTh
   \displaystyle   K(x,z,Q,q_{T}) \simeq \displaystyle
     \sum_{n=1}^N w_n^{(1)}
            \,S\left(x,z,\frac{\xi_n^{(1)}}{q_{T}};\mu,Q^2\right) f_{\rm NP}\left(x, \frac{\xi_n^{(1)}}{q_{T}},Q^2\right) D_{\rm NP}\left(z, \frac{\xi_n^{(1)}}{q_{T}},Q^2\right)\,,

with :math:`S` given in Eq. :eq:`eq:Sall`. Once :math:`K`
is known, the integral of the cross section over the bin
:math:`q_{T}\in [p_{Th,\rm min}/z: p_{Th,\rm max}/z]` is computed as:

.. math::
   \int_{p_{Th,\rm min}/z}^{p_{Th,\rm max}/z} dq_{T}\left[\frac{d\sigma}{dxdQdz d
         q_{T}}\right]=K(x,z,Q, p_{Th,\rm max}/z)-K(x,z,Q, p_{Th,\rm min}/z)\,.

This allows one to compute analytically one of the integrals that are
often required to compare predictions to data.

.. container::
   :name: integrating-over-x-z-and-q

   .. rubric:: Integrating over :math:`x`, :math:`z`, and :math:`Q`
      :name: integrating-over-x-z-and-q

We now move to considering the integral of the cross section over
:math:`x`, :math:`z`, and :math:`Q`. Since these integrals usually come
together with an integration in :math:`q_{T}`, in the following we will
consider the primitive function :math:`K` in
Eq. :eq:`eq:primitivepTh` rather than the cross
section itself, that is:

.. math::
   :label: eq:IntcrosssectionSIDISPrim
     \widetilde{K}(p_{Th})=\int_{Q_{\rm min}}^{Q_{\rm max}}dQ \int_{z_{\rm
         min}}^{z_{\rm max}}dz \int_{x_{\rm min}}^{x_{\rm max}}dx\,
     K(x,z,Q,p_{Th}/z)\,,

\ so that:

.. math:: \widetilde{\sigma}=\widetilde{K}(p_{Th,\rm max})-\widetilde{K}(p_{Th,\rm min})\,.

The amount of numerical computation required to carry out the
integration of a single bin is very large. Indicatively, it amounts to
computing a three-dimensional integral for each of the terms of the
Ogata quadrature that usually range from a few tens to hundreds.
Therefore, in order to be able to do the integrations in a reasonable
amount of time and yet obtain accurate results, it is necessary to put
in place an efficient integration strategy. This goal can be achieved by
exploiting a numerical integration based on interpolation techniques to
precompute the relevant quantities. To this purpose, we first define one
grid in :math:`x`, :math:`\{x_\alpha\}` with :math:`\alpha=0,\dots,N_x`,
one grid in :math:`z`, :math:`\{z_\beta\}` with
:math:`\beta=0,\dots,N_z`, and one grid in :math:`Q`, :math:`\{Q_\tau\}`
with :math:`\tau=0,\dots,N_Q`, each of which with a set of interpolating
functions :math:`\mathcal{I}` associated. The grids should be such to
span the full kinematic range covered by given data set. Then the value
of :math:`K` in Eq. :eq:`eq:primitivepTh` for any
kinematics can be obtained through interpolation as:

.. math::
   \begin{array}{rcl}
   \displaystyle   K(x,z,Q,q_T) &\simeq& \displaystyle 
     \sum_{n=1}^N w_n^{(1)} S\left(x,z,\frac{\xi_n^{(1)}}{q_T};\mu,Q^2\right)
            \sum_{\alpha=1}^{N_x}\sum_{\beta=1}^{N_z}\sum_{\tau=1}^{N_Q}
                                            
                                            \mathcal{I}_\alpha(x)
                                            \mathcal{I}_\beta(z)
                                            \mathcal{I}_\tau(Q) \\
   \\
   &\times&\displaystyle f_{\rm NP}\left(x_\alpha, \frac{\xi_n^{(1)}}{ q_{T}},Q_\tau^2\right) D_{\rm NP}\left(z_\beta, \frac{\xi_n^{(1)}}{q_{T}},Q_\tau^2\right)\,.
   \end{array}

\ Once we have :math:`K` in this form, the integration over :math:`x`,
:math:`z`, and :math:`Q` in
Eq. :eq:`eq:IntcrosssectionSIDISPrim`
does not involve the non-perturbative functions :math:`f_{\rm NP}` and
:math:`D_{\rm NP}` and can be written as:

.. math::
   :label: eq:numintDIS
     \widetilde{K}(p_{Th}) = \sum_{n=1}^N 
                                                    \sum_{\alpha=1}^{N_x}\sum_{\beta=1}^{N_z}\sum_{\tau=1}^{N_Q}
                                                    W_{n\alpha\beta\tau}(p_{Th}) f_{\rm NP}\left(x_\alpha, \frac{z_\beta\xi_n^{(1)}}{p_{Th}},Q_\tau^2\right) D_{\rm NP}\left(z_\beta, \frac{z_\beta\xi_n^{(1)}}{p_{Th}},Q_\tau^2\right)\,,

with:

.. math::
   :label: eq:WDIS
     W_{n\alpha\beta\tau}(p_{Th}) = 
     w_n^{(1)} \int_{Q_{\rm min}}^{Q_{\rm
         max}}dQ\,\mathcal{I}_\tau(Q) \int_{z_{\rm min}}^{z_{\rm max}}dz\,
     \mathcal{I}_\beta(z) \int_{x_{\rm min}}^{x_{\rm
         max}}dx\, \mathcal{I}_\alpha(x) S\left(x,z,\frac{z\xi_n^{(1)}}{p_{Th}};\mu,Q^2\right)\,.

Since the aim is to fit the functions :math:`f_{\rm NP}` and
:math:`D_{\rm NP}` to data, one can precompute and store the
coefficients :math:`W` defined in Eq. :eq:`eq:WDIS` and
compute the cross sections in a fast way making use of
Eq. :eq:`eq:numintDIS`.

It is often the case that the integrated cross section,
Eq. :eq:`eq:IntcrosssectionSIDIS`, is given
within a certain acceptance region which is typically defined as:

.. math::
   W=\sqrt{\frac{(1-x)Q^2}{x}}\geq W_{\rm min}\,,\quad y_{\rm min}\leq y
   \left(=\frac{Q^2}{sx}\right) \leq y_{\rm max}\,.

\ These constraints can be expressed as constraints on the variable
:math:`x` for a fixed value of :math:`Q`:

.. math::
   x\leq \frac{Q^2}{W_{\rm min}^2+Q^2}\,,\quad x\geq\frac{Q^2}{s y_{\rm
         max}}\,,\quad x\leq 
     \frac{Q^2}{sy_{\rm min}}\,.

Therefore, in order to implement the acceptance cuts in the computation
of the integrated cross sections, it is enough to replace the
integration bounds of the integral in :math:`x` in
Eq. :eq:`eq:IntcrosssectionSIDIS` as
follows:

.. math::
   x_{\rm min}\rightarrow \overline{x}_{\rm min}(Q)=\mbox{max}\left[x_{\rm min},\frac{Q^2}{s y_{\rm
           max}}\right]\,,\quad   x_{\rm max}\rightarrow \overline{x}_{\rm max}(Q)=\mbox{min}\left[x_{\rm max},\frac{Q^2}{s y_{\rm
           min}}, \frac{Q^2}{W_{\rm min}^2+Q^2}\right]\,.

Given the number and complexity of integrals as that in
Eq. :eq:`eq:WDIS` to be computed, a numerical
implementation requires devising a strategy that maximises the
efficiency. In order to do so, we need to “unpac” the function :math:`S`
and perform the integrations in a way that unnecessary computations are
avoided as much as possible. Taking :math:`\mu = Q` and
:math:`b = z\xi_n^{(1)} /p_{Th}`, using the fact that the Sudakov
evolution factor for quarks :math:`R_q` is universal, and dropping the
unnecessary dependencies, Eq. :eq:`eq:Sall` can be
conveniently reorganised as:

.. math::
   \begin{array}{rcl}
   \displaystyle
     S\left(x,z,\frac{z\xi_n^{(1)}}{p_{Th}},Q\right)&=&\displaystyle
                                                        4\pi \frac{\alpha^2(Q)H(Q)}{Q^3}R_q^2\left(b_*\left(\frac{z\xi_n^{(1)}}{p_{Th}}\right),Q\right)\\
   \\
   &\times&\displaystyle \frac1{z}\sum_qe_q^2\overline{D}_q\left(z,b_*\left(\frac{z\xi_n^{(1)}}{p_{Th}}\right)\right)\frac{Y_+(x,Q)}{x}\overline{F}_q\left(x,b_*\left(\frac{z\xi_n^{(1)}}{p_{Th}}\right)\right)\,.
   \end{array}

\ with:

.. math:: \overline{F}_q(x,b) = \sum_{j}\int_x^1dy\,\mathcal{C}_{qj}(y)\left[\frac{x}{y}f_j\left(\frac{x}{y},\frac{b_0}{b}\right)\right]\,,

and:

.. math:: \overline{D}_{i}(z,b) = \sum_{j}\int_z^1dy\,\left[y^2\mathbb{C}_{qj}(y)\right]\left[\frac{z}{y}d_j\left(\frac{z}{y}, \frac{b_0}{b}\right)\right]\,,

being the initial-scale TMD distributions before Sudakov evolution. The
weights in Eq. :eq:`eq:WDIS` can then be computed as:

.. math::
   \begin{array}{rcl}
   \displaystyle
     W_{n\alpha\beta\tau}(p_{Th})&=&\displaystyle w_n^{(1)}
                                     4\pi
            \int_{Q_{\rm min}}^{Q_{\rm
         max}}dQ\,\mathcal{I}_\tau(Q) \frac{\alpha^2(Q)H(Q)}{Q^3}\\
   \\
   &\times&\displaystyle\int_{z_{\rm min}}^{z_{\rm max}}dz\,
     \mathcal{I}_\beta(z) 
                                                        R_q^2\left(b_*\left(\frac{z\xi_n^{(1)}}{p_{Th}},Q\right),Q\right)\frac1{z}\sum_qe_q^2 \overline{D}_q\left(z,b_*\left(\frac{z\xi_n^{(1)}}{p_{Th}},Q\right)\right)\\
   \\
   &\times&\displaystyle \int_{\overline{x}_{\rm min}(Q)}^{\overline{x}_{\rm
         max}(Q)}dx\, \mathcal{I}_\alpha(x)\frac{Y_+(x,Q)}{x}\overline{F}_q\left(x,b_*\left(\frac{z\xi_n^{(1)}}{p_{Th}},Q\right)\right)\,.
   \end{array}

This nesting should optimise the computation of the integrals.
