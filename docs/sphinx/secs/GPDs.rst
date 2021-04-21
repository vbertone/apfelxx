================================
Generalised parton distributions
================================

:Author: Valerio Bertone

.. contents::
   :depth: 3
..

Introduction
============

In this set of notes I collect the technical aspects concerning
generalised parton distributions (GPDs). Since the computation GPDs
introduces new kinds of convolution integrals, a strategy aimed at
optimising the numerics needs to be devised.

Evolution equation
==================

The evolution equation for GPDs reads: [1]_

.. math::
   :label: eq:eveq
   \mu^2\frac{d}{d\mu^2}f(x,\xi) = \int_{-1}^{1}\frac{dx'}{\left|2\xi\right|}\mathbb{V}\left(\frac{x}{\xi},\frac{x'}{\xi}\right)f(x',\xi)\,.

In general, the GPD :math:`f` and the evolution kernel
:math:`\mathbb{V}` should be respectively interpreted as a vector and a
matrix in flavour space. However, for now, we will just be concerned
with the integral in the r.h.s. of Eq. :eq:`eq:eveq`
regardless of the flavour structure.

The support of the evolution kernel
:math:`\mathbb{V}\left(\frac{x}{\xi},\frac{x'}{\xi}\right)` is shown in
Fig. `1 <#fig:GPDIntDomain>`__.

.. container::
   :name: fig:GPDIntDomain

   .. figure:: ../../latex/src/../../latex/src/plots/GPDIntDomain.pdf
   .. :alt: Support domain of the evolution kernel.[fig:GPDIntDomain]
      :name: fig:GPDIntDomain
      :width: 70.0%

      Support domain of the evolution kernel.[fig:GPDIntDomain]

The Knowledge of the support region of the evolution kernel allows us to
rearrange Eq. :eq:`eq:eveq` as follows:

.. math:: \displaystyle\mu^2\frac{d}{d\mu^2}f(\pm x,\xi) =\int_{b(x)}^{1}\frac{dx'}{x'}\left[\frac{x'}{\left|2\xi\right|}\mathbb{V}\left(\pm \frac{x}{\xi},\frac{x'}{\xi}\right)f(x',\xi)+\frac{x'}{\left|2\xi\right|}\mathbb{V}\left(\mp \frac{x}{\xi},\frac{x'}{\xi}\right)f(-x',\xi)\right]\,,

with:

.. math::
   :label: eq:lowintb
   b(x) = |x|\theta\left(\left|\frac{x}{\xi}\right|-1\right)\,,

and where we have used the symmetry property of the evolution kernels:
:math:`\mathbb{V}(y,y')=\mathbb{V}(-y,-y')`. In the unpolarised case, it
is useful to define: [2]_

.. math::
   :label: eq:pmdef
   \begin{array}{rcl}
   \displaystyle f^{\pm}(x,\xi) &=&\displaystyle  f(x,\xi) \mp
                          f(-x,\xi)\,,\\
   \\
   \displaystyle \mathbb{V}^{\pm}(y,y') &=&\displaystyle  \mathbb{V}(y,y') \mp \mathbb{V}(-y,y')\,,
   \end{array}

so that the evolution equation for :math:`f^{\pm}` reads:

.. math::
   :label: eq:eveq2
   \displaystyle\mu^2\frac{d}{d\mu^2}f^{\pm}(x,\xi) = \int_{b(x)}^{1}\frac{dx'}{x'}\frac{x'}{\left|2\xi\right|}
                                                            \mathbb{V}^{\pm}\left(\frac{x}{\xi},\frac{x'}{\xi}\right)f^{\pm}(x',\xi)\,.

The :math:`f^{\pm}` distributions can be regarded as the GPD analogous
of the :math:`\pm` forward distributions that can then be used to
construct the usual singlet and non-singlet distributions in the QCD
evolution basis. This uniquely determines the flavour structure of the
evolution kernels :math:`\mathbb{V}^{\pm}`.

It is relevant to observe that the presence of the
:math:`\theta`-function in the lower integration bound :math:`b`,
Eq. :eq:`eq:lowintb`, is such that for :math:`|x|>|\xi|`
the evolution equation has the exact form of the DGLAP evolution
equation which corresponds to integrating over the blue regions in
Fig. `1 <#fig:GPDIntDomain>`__ (DGLAP region, henceforth). Conversely,
for :math:`|x|\leq|\xi|` the lower integration bound becomes zero and
the evolution equation assumes the form of the so-called ERBL equation
that describes the evolution of meson distribution amplitudes (DAs).
This corresponds to integrating over the red region (ERBL region,
henceforth). Crucially, in the limits :math:`\xi\rightarrow 0` and
:math:`\xi\rightarrow \pm1` Eq. :eq:`eq:eveq2` needs to
recover the DGLAP and ERBL equations, respectively.

.. container::
   :name: sec:endpoint

   .. rubric:: End-point contributions
      :name: sec:endpoint

Some of the expressions for the anomalous dimensions discussed below
contain :math:`+`-prescribed terms. It is thus important to treat these
terms properly. We are generally dealing with objects defined as:

.. math::
   :label: eq:pludistributionn
     \left[\mathbb{V}\left(x,x'\right)\right]_+=
     \mathbb{V}\left(x,x'\right)- \delta(x-x')\int_{-1}^{1}dx\,\mathbb{V}\left(x,x'\right)\,.

where the function :math:`\mathbb{V}` has a pole at :math:`x'=x`.

Let us take as an example the one-loop non-singlet anomalous dimension.
For definiteness, we will refer for the precise expression to Eq. (101)
of Ref. (Diehl 2003) and report it here for convenience (up to a factor
:math:`\alpha_s/4\pi`):

.. math::
   :label: eq:diehlexpr
   V_{\rm NS}^{(0)}(x,x') = 2C_F\left[\rho(x,x')\left\{\frac{1+x}{1+x'}\left(1+\frac{2}{x'-x}\right)\right\}+(x\rightarrow -x, x'\rightarrow -x')\right]_+\,,

with: [3]_

.. math::
   :label: eq:supportdiehl
   \rho(x,x')=\theta(-x + x')\theta(1 + x) - \theta(x - x')\theta(1 - x)

In order for Eq. :eq:`eq:diehlexpr` to be consistent
with the forward evolution, one should find:

.. math::
   :label: eq:forwardlimit
     \lim_{\xi\rightarrow 0}\frac{1}{|2\xi|} V_{\rm NS}^{(0)}\left(\frac{x}{\xi},\frac{x'}{\xi}\right) \mathop{=}^?
     \frac{1}{x'} P_{\rm NS}\left(\frac{x}{x'}\right) = \frac{1}{x'}2C_F\left[\theta\left(\frac{x}{x'}\right)\theta\left(1-\frac{x}{x'}\right)\frac{1+\left(\frac{x}{x'}\right)^2}{1-\left(\frac{x}{x'}\right)}\right]_+\,,

such that Eq. :eq:`eq:eveq` exactly reduces to the DGLAP
equation. However, if one takes the explicit limit for
:math:`\xi\rightarrow 0` of Eq. :eq:`eq:diehlexpr` one
finds: [4]_

.. math::
   :label: eq:forwardlimit2
   \lim_{\xi\rightarrow 0}\frac{1}{|2\xi|} V_{\rm NS}^{(0)}\left(\frac{x}{\xi},\frac{x'}{\xi}\right) = 2C_F\left[\frac{1}{x'}\theta\left(\frac{x}{x'}\right)\left(1-\frac{x}{x'}\right)\frac{1+\left(\frac{x}{x'}\right)^2}{1-\left(\frac{x}{x'}\right)}\right]_+\,.

Therefore, as compared to
Eq. :eq:`eq:forwardlimit`, the factor :math:`1/x'`
in Eq. :eq:`eq:forwardlimit2` appears *inside* the
+-prescription sign rather than outside which makes the two expressions
effectively different under integration. The difference amounts to a
local term that can be quantified by knowing that:

.. math:: \left[yg(y) \right]_+=y\left[g(y)\right]_+ + \delta(1-y)\int_0^1dz\,(1-z)g(z)\,.

Notice that, thanks to the factor :math:`(1-z)`, the integral in the
r.h.s. of the above equation converges despite the singularity of
:math:`g`. For example:

.. math::
   \left[\frac{y}{1-y}\right]_+=y\left[\frac1{1-y}\right]_+
     + \delta(1-y)\,.

Finally, one finds that the forward limit of
Eq. :eq:`eq:diehlexpr` gives:

.. math::
   \lim_{\xi\rightarrow 0}\frac{1}{|2\xi|} V_{\rm
       NS}^{(0)}\left(\frac{x}{\xi},\frac{x'}{\xi}\right) =
     \frac{1}{x'} \left[P_{\rm NS}\left(\frac{x}{x'}\right)+\frac{4}{3}C_F\delta\left(1-\frac{x}{x'}\right)\right]\,,

which does *not* reproduce the DGLAP equation due to the presence of an
additional local term.

.. container::
   :name: on-vinnikovs-code

   .. rubric:: On Vinnikov’s code
      :name: on-vinnikovs-code

The purpose of this section is to draw the attention on a possible
incongruence of the GPD evolution code developed by Vinnikov and
presented in Ref. (Vinnikov 2006). For definiteness, we concentrate on
the non-singlet :math:`H_{\rm NS}` GPD in the DGLAP region
:math:`x>\xi`, whose evolution equation is given in Eq. (29). For
convenience, we report that equation here:

.. math::
   :label: eq:Vinnikov
   \begin{array}{rcl}
   \displaystyle \frac{d H_{\rm NS}(x,\xi,Q^2)}{d\ln
     Q^2}&=&\displaystyle\frac{2\alpha_s(Q^2)}{3\pi}\Bigg[\int_x^1 dy
             \frac{x^2+y^2-2\xi^2}{(y-x)(y^2-\xi^2)}\left(H_{\rm
             NS}(y,\xi,Q^2)-H_{\rm NS}(x,\xi,Q^2)\right)\\
   \\
   &+&\displaystyle H_{\rm
       NS}(x,\xi,Q^2)\bigg(\frac32+2\ln(1-x)+\frac{x-\xi}{2\xi}\ln((x-\xi)(1+\xi))\\
   \\
   &-&\displaystyle \frac{x+\xi}{2\xi}\ln((x+\xi)(1-\xi))\bigg)\Bigg]\,,
   \end{array}

\ and take the forward limit :math:`\xi\rightarrow 0`, obtaining:

.. math::
   \begin{array}{rcl}
   \displaystyle \frac{d H_{\rm NS}(x,0,Q^2)}{d\ln
     Q^2}&=&\displaystyle\frac{2\alpha_s(Q^2)}{3\pi}\Bigg[\int_x^1 dy
             \frac{x^2+y^2}{y^2 (y-x)}\left(H_{\rm
             NS}(y,0,Q^2)-H_{\rm NS}(x,0,Q^2)\right)\\
   \\
   &+&\displaystyle H_{\rm  NS}(x,0,Q^2)\left(\frac32+2\ln(1-x)\right)\Bigg]\,,
   \end{array}

The limit for :math:`\xi\rightarrow 0` of the equation above should
reproduce the usual DGLAP evolution equation:

.. math::
   \displaystyle \frac{d H_{\rm 
   NS}(x,0,Q^2)}{d\ln
     Q^2}=\frac{\alpha_s(Q^2)}{4\pi}\int_x^1 \frac{dy}{y}\left[\hat{P}_{\rm
     NS}\left(\frac{x}{y}\right)\right]_+H_{\rm NS}\left(y,0,Q^2\right)\,,

where:

.. math:: \hat{P}_{\rm NS}\left(z\right)=2C_F\frac{1+z^2}{1-z}\,,

with :math:`C_F=4/3`. Written explicitly and accounting for the
additional local term arising from the incompleteness of the convolution
integral, one finds:

.. math::
   :label: eq:DGLAPexpl
   \begin{array}{rcl}
     \displaystyle \frac{d H_{\rm NS}(x,0,Q^2)}{d\ln
       Q^2}&=&\displaystyle \frac{2\alpha_s(Q^2)}{3\pi}\Bigg[\int_x^1
               dy\,\frac{x^2+y^2}{y^3(y-x)}\left(y
               H_{\rm NS}\left(y,0,Q^2\right)-xH_{\rm
         NS}(x,0,Q^2)\right)\\
   \\
    &+&\displaystyle H_{\rm NS}(x,0,Q^2)\left(\frac{x(x+2)}{2}+2\ln(1-x)\right)\Bigg]\,,
   \end{array}

which evidently differs from Eq. :eq:`eq:Vinnikov`. By
inspection, one observes that the difference can be partially traced
back to the issue discussed in Sect. (`2.1 <#sec:endpoint>`__). An
interesting observation is that, for :math:`x\rightarrow 1`, the two
expressions tend to coincide. This means that the difference is larger
at small values of :math:`x`. This fact may have concurred to cause the
oversight of this discrepancy in past numerical comparisons.

.. container::
   :name: on-jis-evolution-equation

   .. rubric:: On Ji’s evolution equation
      :name: on-jis-evolution-equation

In this section we discuss the evolution equations derived by Ji in
Ref. (Ji 1997). This form of the evolution equation is dubbed
“near-forward” in Ref. (Blumlein, Geyer, and Robaschik 1999) because it
closely resembles the DGLAP equation. However, in Ref. (Ji 1997) two
different equations apply to the regions :math:`x<\xi` and
:math:`x>\xi`. In this section, we will unify them showing that the
resulting one-loop non-singlet off-forward anomalous dimension cannot be
written as a fully :math:`+`-prescribed distribution.

We start by considering Eqs. (15)-(17) of Ref. (Ji 1997). The first step
is to replace :math:`\xi/2` with :math:`\xi` to match our notation. Then
we consider the subtraction integrals in Eq. (16) keeping in mind that
they apply to both regions :math:`x<\xi` and :math:`x>\xi`: [5]_

.. math::
   \int_{\pm \xi}^x\frac{dy}{y-x} = -\int_{\pm \kappa}^1\frac{dz}{1-z}
     = -\int_{0}^1\frac{dz}{1-z}+\int_{1\mp \kappa}^1\frac{dt}{t}=-\int_{0}^1\frac{dz}{1-z}-\ln(|1\mp\kappa|)\,,

with:

.. math::
   :label: eq:kappadef
     \kappa = \frac{\xi}{x}\,,

such that the full local term in Eq. (16) becomes:

.. math::
   \frac{3}{2}+ \int_{\xi}^x\frac{dy}{y-x}+\int_{-\xi}^x\frac{dy}{y-x}
     = \frac{3}{2}-2\int_{0}^1\frac{dz}{1-z}-\ln\left(|1-\kappa^2|\right)\,,

Considering the symmetry for :math:`\xi\leftrightarrow -\xi` of the
evolution kernel in Eq. (17) of Ref. (Ji 1997), we can write Eq. (15)
valid for :math:`\kappa<1` in a more compact way as:

.. math::
   :label: eq:ji1
   \mu^2\frac{d}{d\mu^2}f^-(x,\xi) = \frac{\alpha_s(\mu)}{4\pi}\int_x^1\frac{dy}{y}\mathcal{P}_1^{-,(0)}(y,\kappa)f^-\left(\frac{x}{y},\xi\right)\,,

with:

.. math::
   :label: eq:p1minus0
   \begin{array}{rcl}
   \displaystyle \mathcal{P}_1^{-,(0)}(y,\kappa) &=& \displaystyle
                                                     2P_{\rm NS}(y,2\kappa y) + \delta(1-y)2C_F\left(\frac{3}{2}-2\int_{0}^{1}\frac{dz}{1-z}-\ln(|1-\kappa^2|)\right) \\
   \\
   &=& \displaystyle  2C_F\left\{\left(\frac{2}{1-y}\right)_+-\frac{1
     +y}{1-\kappa^2y^2}+\delta(1-y)\left[\frac32-
     \ln\left(|1-\kappa^2|\right)\right]\right\}\\
   \\
   &=& \displaystyle  2C_F\left\{\left[\frac{1
     +(1-2\kappa^2)y^2}{(1-y)(1-\kappa^2y^2)}\right]_++\delta(1-y)\left[\frac32+\left(\frac{1}{2\kappa^2}-1\right)
     \ln\left(|1-\kappa^2|\right)+\frac{1}{2\kappa}\ln\left(\left|\frac{1-\kappa}{1+\kappa}\right|\right)\right]\right\}\,,
   \end{array}

where :math:`P_{\rm NS}` is given in Eq. (17) of Ref. (Ji 1997). The
splitting function :math:`\mathcal{P}_1^{-,(0)}` is such that:

.. math::
   :label: eq:nonzeroplusGPD
   \int_0^1dy\,\mathcal{P}_1^{-,(0)}(y,\kappa) =2C_F\left[\frac32+\left(\frac{1}{2k^2}-1\right) \ln\left(|1-\kappa^2|\right)+\frac{1}{2\kappa}\ln\left(\left|\frac{1-\kappa}{1+\kappa}\right|\right)\right]\,,

which means that it cannot be written as a fully :math:`+`-prescribed
distribution. However, the integral above correctly tends to zero as
:math:`\kappa\rightarrow 0` allowing one to recover the usual DGLAP
splitting function in the forward limit:

.. math::
   :label: eq:dglapsplitting
     \lim_{\kappa\rightarrow 0}\mathcal{P}_1^{-,(0)}(y,\kappa) = 2C_F\left[\frac{1+y^2}{1-y}\right]_+\,.

It should also be pointed out that also the limit for
:math:`\kappa\rightarrow 1` of Eq. :eq:`eq:p1minus0` is
well-behaved:

.. math:: \lim_{\kappa\rightarrow 1}\mathcal{P}_1^{-,(0)}(y,\kappa) = 2C_F\left\{\left[\frac{1}{1-y}\right]_++\delta(1-y)\left[\frac32-\ln(2)\right]\right\}\,.

which is necessary to have a smooth transition of the GPDs from the
DGLAP (:math:`x>\xi`) to the ERBL (:math:`x<\xi`) region.

We now consider Eqs. (18) and (19) of Ref. (Ji 1997) valid for
:math:`\kappa>1`. Interestingly, after some algebra, we find:

.. math::
   :label: eq:ji2
   \mu^2\frac{d}{d\mu^2}f^-(x,\xi) = \frac{\alpha_s(\mu)}{4\pi}\left[\int_x^1\frac{dy}{y}\mathcal{P}_1^{-,(0)}(y,\kappa)f^-\left(\frac{x}{y},\xi\right)+\int_x^\infty\frac{dy}{y}\mathcal{P}_2^{-,(0)}(y,\kappa)f^-\left(\frac{x}{y},\xi\right)\right]\,,

with :math:`\mathcal{P}_1^{-,(0)}` given by:

.. math::
   :label: eq:p1minus02
   \displaystyle \mathcal{P}_1^{-,(0)}(y,\kappa)=2P_{\rm NS}'(y,2\kappa y) +2P_{\rm NS}'(y,-2\kappa y)+ \delta(1-x)2C_F\left(\frac{3}{2}-2\int_{0}^{1}\frac{dy}{1-y}-\ln(|1-\kappa^2|)\right)\,,

with :math:`P_{\rm NS}'` is given in Eq. (19) of Ref. (Ji 1997) and
remarkably equal to the expression in
Eq. :eq:`eq:p1minus0` signifying that:

.. math::
   :label: eq:primevsnonprime
     P_{\rm NS}(y,2\kappa y) =  P_{\rm NS}'(y,2\kappa y) +P_{\rm NS}'(y,-2\kappa y)\,.

While:

.. math::
   :label: eq:p2minus0
     \mathcal{P}_2^{-,(0)}(y,\kappa) = -2P_{\rm NS}'(y,-2\kappa y) +2P_{\rm
       NS}'(-y,2\kappa y) = 2C_F(\kappa-1)\frac{y+(1+2\kappa)y^3}{(1-y^2)(1-\kappa^2y^2)}\,.

It is very interesting to notice that :math:`\mathcal{P}_2^{-,(0)}` is
proportional to :math:`(\kappa-1)` that finally guarantees the
continuity of GPDs at :math:`\kappa=1`.

We observe that, within the integration interval, the splitting function
:math:`\mathcal{P}_2^{-,(0)}` is singular at :math:`y=1`. [6]_ However,
as pointed out above, the second integral on the r.h.s. of
Eq. :eq:`eq:ji2` has to be regarded as principal-valued
therefore it is well-defined. In order to treat this integral
numerically we consider the specific case:

.. math:: I=\int_x^\infty dy\,\frac{f(y)}{1-y}\,,

\ where :math:`f` is a test function well-behaved over the full
integration range. If one subtracts and adds back the divergence at
:math:`y=1`, *i.e.*:

.. math:: f(1)\int_0^1\frac{dy}{1-y}\,,

one can rearrange the integral as follows:

.. math::
   :label: eq:plusplusdist
   I=\int_x^\infty\frac{dy}{1-y}\left[f(y)-f(1)\left(1+\theta(y-1)\frac{1-y}{y}\right)\right]+f(1)\ln(1-x)\equiv
   \int_x^\infty dy\left(\frac{1}{1-y}\right)_{++}f(y)\,,

which effectively defines the :math:`++`-distribution. It should be
noticed that this definition is specific to the function
:math:`1/(1-y)`. In case of a different singular function the function
that multiplies :math:`\theta(y-1)` would be different. The advantage of
this rearrangement is that the integrand is free of the divergence at
:math:`y=1` and is thus amenable to numerical integration. Also, the
:math:`++`-distribution reduces to the standard :math:`+`-distribution
when the upper integration bound is one rather than infinity. In this
sense the :math:`++`-distribution generalises the :math:`+`-distribution
to ERBL-like integrals.

In view of the use of Eq. :eq:`eq:plusplusdist`, it
is convenient to rewrite Eq. :eq:`eq:p2minus0` as
follows:

.. math::
   :label: eq:p2minus01
     \mathcal{P}_2^{-,(0)}(y,\kappa) =
     2C_F\left[\frac{1+(1+\kappa)y+(1+\kappa -\kappa^2)y^2}{(1+y)(1-\kappa^2y^2)}-\left(\frac{1}{1-y}\right)_{++}\right]\,,

where the first term in the squared bracket is regular at :math:`y=1`.

Finally, Eqs. :eq:`eq:ji1` and Eq. :eq:`eq:ji2`
can be combined as follows:

.. math::
   :label: eq:jitot
     \mu^2\frac{d}{d\mu^2}f^-(x,\xi) = \frac{\alpha_s(\mu)}{4\pi}\left[\int_x^1\frac{dy}{y}\mathcal{P}_1^{-,(0)}(y,\kappa)f^-\left(\frac{x}{y},\xi\right)+\theta(\kappa-1)\int_x^\infty\frac{dy}{y}\mathcal{P}_2^{-,(0)}(y,\kappa)f^-\left(\frac{x}{y},\xi\right)\right]\,,

or even more compactly as:

.. math::
   :label: eq:jitotcomp
     \mu^2\frac{d}{d\mu^2}f^-(x,\xi) = \frac{\alpha_s(\mu)}{4\pi}\int_x^\infty\frac{dy}{y}\mathcal{P}^{-,(0)}(y,\kappa)f^-\left(\frac{x}{y},\xi\right)\,,

with:

.. math:: \mathcal{P}^{-,(0)}(y,\kappa) = \theta(1-y)\mathcal{P}_1^{-,(0)}(y,\kappa)+\theta(\kappa-1)\mathcal{P}_2^{-,(0)}(y,\kappa)\,,

to obtain a single DGLAP-like evolution equation valid for all values of
:math:`\kappa`. In fact, it should be pointed out that, when performing
the integrals numerically, the form in
Eq. :eq:`eq:jitotcomp` has to be adopted. The reason
is that both functions :math:`\mathcal{P}_1^{-,(0)}` and
:math:`\mathcal{P}_2^{-,(0)}`, due to the factor :math:`1-\kappa^2y^2`,
are affected by a pole at :math:`y = |\kappa|^{-1}` that, for
:math:`|\kappa|>1` or equivalently :math:`|x|<|\xi|` (*i.e.* in the ERBL
region) have to cancel to give a finite result. Using the explicit
expressions for :math:`\mathcal{P}_1^{-,(0)}` and
:math:`\mathcal{P}_2^{-,(0)}`, we find:

.. math::
   \lim_{y\rightarrow \kappa^{-1}} (1-\kappa^2y^2)
   \mathcal{P}_1^{-,(0)}(y,\kappa) = -2 C_F \frac{1+\kappa}{\kappa}\,,

and:

.. math::
   \lim_{y\rightarrow \kappa^{-1}} (1-\kappa^2y^2)
   \mathcal{P}_2^{-,(0)}(y,\kappa) = 2 C_F \frac{1+\kappa}{\kappa}\,.

Since the coefficient of the pole are equal in absolute value and
opposite in sign they cancel in the integral. Below, we will explicitly
verify this property also for the anomalous dimensions of the singlet
sector.

Importantly, in the limit :math:`\kappa\rightarrow 0`, the second
integral in the r.h.s. of Eq. :eq:`eq:jitot` drops and the
splitting function :math:`\mathcal{P}_1^{-,(0)}` reduces to the one-loop
non-singlet DGLAP splitting function (see
Eq. :eq:`eq:dglapsplitting`) so that, as
expected, Eq. :eq:`eq:jitot` becomes the DGLAP equation.

.. container::
   :name: the-erbl-equation

   .. rubric:: The ERBL equation
      :name: the-erbl-equation

It is also interesting to verify that also the ERBL equation is
recovered in the limit :math:`\xi\rightarrow 1`. Given the definition of
:math:`\kappa`, Eq. :eq:`eq:kappadef`, this limit is
attained by taking :math:`\kappa \rightarrow 1/x`. However, the limit
procedure is more subtle than in the DGLAP case due to the presence of
:math:`+`-prescriptions and explicit local terms that need to cooperate
to give the right result.

We make use of Eqs. :eq:`eq:p1minus02`
and :eq:`eq:p2minus0` to write the evolution equation
in terms of the function :math:`P_{\rm NS}'` in a form similar to that
originally given in Ref. (Ji 1997) but more compactly as:

.. math::
   :label: eq:erbl2
     \mu^2\frac{d}{d\mu^2}f^-(x,1) =
     \frac{\alpha_s(\mu)}{4\pi}\left[\int_{-1}^1 dy\,V_{\rm NS}^{(0)}(x,y)f^-\left(y,1\right)\right] \,.

with:

.. math::
   :label: eq:ERBLkernelxy
   \begin{array}{rcl}
   V_{\rm NS}^{(0)}(x,y) &=&\displaystyle \theta(y-x)\left[\frac{2}{y}P_{\rm
       NS}'\left(\frac{x}{y},\frac{2}{y}\right)
                             \right]-2C_F\delta\left(y-x\right)\int_{-1}^1 dz \,\frac{\theta(z-x)}{z-x}\\
   \\
   &+&\displaystyle \theta(x-y)\left[-\frac{2}{y}P_{\rm
       NS}'\left(\frac{x}{y},-\frac{2}{y}\right)
                             \right]+2C_F\delta\left(x-y\right)\int_{-1}^{1}dz\,\frac{\theta(x-z)}{z-x}\\
   \\
   &+&\displaystyle 3C_F\delta\left(y-x\right)\,,
   \end{array}

where:

.. math:: \frac{2}{y}P_{\rm NS}'\left(\frac{x}{y},\frac{2}{y}\right)=2C_F\frac{1+x}{1+y}\left(\frac{1}{2}+\frac{1}{y-x}\right)\,.

In order to make a step towards the ERBL equation, we change the
variables :math:`x` and :math:`y` with:

.. math::
   :label: eq:ERBLvariables
   \begin{array}{l}
   \displaystyle t= \frac12\left(x + 1\right)\,,\\
   \\
   \displaystyle u = \frac12\left(y + 1\right)\,,\\
   \end{array}

such that the evolution variable becomes:

.. math:: \mu^2\frac{d}{d\mu^2}\Phi^-(t) = \frac{\alpha_s(\mu)}{4\pi}\left[\int_{0}^1 du\,\overline{V}_{\rm NS}^{(0)}(t,u) \Phi^-\left(u\right)\right] \,.

with :math:`\Phi^-(t) = f^-(x,1)` and:

.. math::
   :label: eq:anomdim
   \begin{array}{rcl}
   \overline{V}_{\rm NS}^{(0)}(t,u) &=&\displaystyle
   C_F\bigg[\theta(u-t)\left(\frac{t-1}{u}+\frac1{u-t}-\delta(u-t)\int_0^1 du'\frac{\theta(u'-t)}{u'-t}\right)\\
   \\
   &-&\displaystyle \theta(t-u)\left(\frac{t}{1-u}+\frac1{u-t}-\delta(t-u)\int_0^1 du'\frac{\theta(t-u')}{u'-t}\right) +\frac32\delta\left(u-t\right)\bigg]\,.
   \end{array}

Now we define:

.. math:: \left[f(t,u)\right]_+\equiv f(t,u)-\delta(u-t)\int_0^1du'\,f(t,u') \,,

where :math:`f` has a single pole at :math:`u=t`, so that we can write
Eq. :eq:`eq:anomdim` more compactly as:

.. math::
   :label: eq:anomdim1
     \overline{V}_{\rm NS}^{(0)}(t,u) =
     C_F\left\{\left[\theta(u-t)\frac{t-1}{u}+\left(\frac{\theta(u-t)}{u-t}\right)_+\right]-\left[\theta(t-u)\frac{t}{1-u}+\left(\frac{\theta(t-u)}{u-t}\right)_+\right]+\frac32\delta\left(u-t\right)\right\}\,.

This confirms the result of Ref. (Blumlein, Geyer, and Robaschik 1999)
modulo the fact that, for achieving a correct cancellation of the
divergencies, the :math:`\theta`-function for the :math:`+`-prescribed
terms needs to be inside the :math:`+`-prescription sign itself rather
than outside. In fact, it is the very presence of the
:math:`\theta`-function that generates the necessity of a
:math:`+`-prescription. Without the :math:`\theta`-function, the
integral could be interpreted as principal-valued integral that requires
no regularisation. The presence of the :math:`\theta`-function
interdicts the cancellation between left and right sides of the pole
resulting in a singular integral. The :math:`+`-prescription reabsorbes
this singularity producing a well-behaved integral.

One can check that integrating :math:`\overline{V}_{\rm NS}^{(0)}` over
:math:`t` gives zero: [7]_

.. math::
   :label: eq:nonzeroplusERBL
     \int_0^1dt\,\overline{V}_{\rm NS}^{(0)}(t,u) = 0\,.

\ This finally confirms that :math:`\overline{V}_{\rm NS}^{(0)}` as
derived from Ref. (Ji 1997) admits a fully :math:`+`-prescribed form,
that is:

.. math::
   \overline{V}_{\rm NS}^{(0)}(t,u) =
     C_F\left\{\theta(u-t)\left[\frac{t-1}{u}+\frac{1}{u-t}\right]-\theta(t-u)\left[\frac{t}{1-u}+\frac{1}{u-t}\right]\right\}_+\,.

This was also explicitly derived in Ref. (Mikhailov and Radyushkin 1985)
and argued that this property must hold for symmetry reasons.

It is now interesting to derive an ERBL-like evolution equation for GPDs
at one loop. This equation can be written as:

.. math::
   :label: eq:ERBLlikeeveq
   \frac{d}{d\ln\mu^2} f^{-}(x,\xi)=
   \frac{\alpha_s(\mu)}{4\pi}\int_{-1}^1\frac{dy}{\xi} \mathbb{V}_{\rm NS}^{(0)}\left(\frac{x}{\xi},\frac{y}{\xi}\right) f^{-}(y,\xi)\,,

with:

.. math::
   :label: eq:evkernelERBLgenxi
   \begin{array}{rcl}
   \displaystyle \frac{1}{\xi}\mathbb{V}_{\rm NS}^{(0)}\left(\frac{x}{\xi}, \frac{y}{\xi}\right) &=& \displaystyle
   \frac{2}{y}\bigg\{\theta(x-\xi) \theta(y-x)P_{\rm
                                                 NS}\left(\frac{x}{y},\frac{2\xi}{y}\right)-\theta(-x-\xi)
                                                 \theta(x-y)P_{\rm
                                                 NS}\left(\frac{x}{y},\frac{2\xi}{y}\right)\\
   \\
   &+&\displaystyle \theta(\xi-x)\theta(x+\xi) \left[\theta(y-x)P_{\rm
       NS}'\left(\frac{x}{y},\frac{2\xi}{y}\right)-\theta(x-y)P_{\rm
       NS}'\left(\frac{x}{y},-\frac{2\xi}{y}\right)\right]\\
   \\
   &+&\displaystyle \delta\left(1-\frac{x}{y}\right) C_F\left[\frac{3}{2}+\int_\xi^{x}\frac{dz}{z-x}+\int_{-\xi}^{x}\frac{dz}{z-x}\right]\bigg\}\,.
   \end{array}

We can now use Eq. :eq:`eq:primevsnonprime` to
replace :math:`P_{\rm NS}` with :math:`P_{\rm NS}'` which gives:

.. math::
   \begin{array}{rcl}
   \displaystyle \frac{1}{\xi}\mathbb{V}_{\rm NS}^{(0)}\left(\frac{x}{\xi}, \frac{y}{\xi}\right) &=& \displaystyle
   \frac{2}{y}\bigg\{\left[\theta(x-\xi) \theta(y-x)+\theta(\xi-x)\theta(x+\xi)\theta(y-x)-\theta(-x-\xi)
                                                 \theta(x-y)\right]P_{\rm NS}'\left(\frac{x}{y},\frac{2\xi}{y}\right)\\
   \\
   &+&\displaystyle \left[\theta(x-\xi) \theta(y-x)-\theta(\xi-x)\theta(x+\xi)\theta(x-y)-\theta(-x-\xi)
                                                 \theta(x-y)\right]P_{\rm NS}'\left(\frac{x}{y},-\frac{2\xi}{y}\right)\\
   \\
   &+&\displaystyle \delta\left(1-\frac{x}{y}\right) C_F\left[\frac{3}{2}+\int_\xi^{x}\frac{dz}{z-x}+\int_{-\xi}^{x}\frac{dz}{z-x}\right]\bigg\}\,.
   \end{array}

.. container::
   :name: conformal-moments

   .. rubric:: Conformal moments
      :name: conformal-moments

A relevant question is whether
Eq. :eq:`eq:evkernelERBLgenxi` is such that
the so-called conformal moments of GPDs (see Eq. (111) of Ref. (Diehl
2003)):

.. math:: \mathcal{C}_n(\xi) = \xi^{n}\int_{-1}^{1}dx\,C_n^{3/2}\left(\frac{x}{\xi}\right)f^{-}(x,\xi)\,,

where :math:`C_n^{3/2}` are Gegenbauer polynomials, “diagonalise” the
leading-order evolution equation in
Eq. :eq:`eq:ERBLlikeeveq`. To state it differently,
we would like to check whether conformal moments evolve
multiplicatively. Multiplying
Eq. :eq:`eq:ERBLlikeeveq` by
:math:`\xi^{n}C_n^{3/2}\left(x/\xi\right)` and integrating over
:math:`x` between :math:`-1` and 1, yields:

.. math::
   :label: eq:confmomprel
   \frac{d \mathcal{C}_n(\xi)}{d\ln\mu^2} =
   \frac{\alpha_s(\mu)}{4\pi}\xi^{n}\int_{-1}^1 dy f^{-}(y,\xi)\int_{-1}^{1}\frac{dx}{\xi}\,C_n^{3/2}\left(\frac{x}{\xi}\right) \mathbb{V}_{\rm NS}^{(0)}\left(\frac{x}{\xi},\frac{y}{\xi}\right) \,.

Therefore, the aim is to check whether the following equality holds:

.. math::
   :label: eq:confmomgpd
     \int_{-1}^1\frac{dx}{\xi}\,C_n^{3/2}\left(\frac{x}{\xi}\right) \mathbb{V}_{\rm NS}^{(0)}\left(\frac{x}{\xi},\frac{y}{\xi}\right) =\mathcal{V}_{n}^{(0)}(\xi)C_n^{3/2}\left(\frac{y}{\xi}\right)\,,

where :math:`\mathcal{V}_{n}^{(0)}` is a number that may generally
depend on :math:`\xi`. If Eq. :eq:`eq:confmomgpd`
held, Eq. :eq:`eq:confmomprel` would become:

.. math::
   :label: eq:confmom
     \frac{d \mathcal{C}_n(\xi)}{d\ln\mu^2} =\frac{\alpha_s(\mu)}{4\pi}\mathcal{V}_{n}^{(0)}(\xi)\mathcal{C}_n(\xi)\,,

that is to say that conformal moments evolve multiplicative exactly like
Mellin moments for forward distributions. As a matter of fact, up to
numerical factor, conformal moments do coincide with Mellin moments in
the limit :math:`\xi\rightarrow 0`. This can be seen by observing that
the Gegenbauer polynomials admit the following expansion:

.. math::
   \xi^{n}C_n^{3/2}\left(\frac{x}{\xi}\right)=\frac{n+1}{2^{n+1}}\sum_{k=0}^n{n
   \choose k}{n+2 \choose k+1}(x+\xi)^k (x-\xi)^{n-k}\,,

that is such that:

.. math::
   :label: eq:gegenlimit
   \lim_{\xi\rightarrow 0}
   \xi^{n}C_n^{3/2}\left(\frac{x}{\xi}\right)= \frac{(2n+1)!}{2^n(n!)^2}x^{n}\,.

Therefore, the conformal moments of the non-singlet distribution in the
forward limit become:

.. math:: \lim_{\xi\rightarrow 0}\mathcal{C}_n(\xi) = \frac{(2n+1)!}{2^n(n!)^2}[1+(-1)^n]f^{-}(n+1)\,,

where, with abuse of notation, we have defined the Mellin moments of the
forward distribution as:

.. math:: f^{-}(n) = \lim_{\xi\rightarrow 0} \int_0^1dx\,x^{n-1}f^{-}(x,\xi)\,,

that are known to diagonalise the DGLAP equation to all orders. Using
Eq. :eq:`eq:gegenlimit` and the fact that:

.. math::
   \lim_{\xi\rightarrow 0} \frac{1}{\xi}\mathbb{V}_{\rm NS}^{(0)}\left(\frac{x}{\xi},
     \frac{y}{\xi}\right) = \left[\theta(x)\theta(y-x)-\theta(-x)\theta(x-y)\right]\frac{1}{y} P_{\rm NS}^{(0)}\left(\frac{x}{y}\right)\,,

where :math:`P_{\rm NS}^{(0)}` is the leading-order non-singlet DGLAP
splitting function, one finds:

.. math::
   \lim_{\xi\rightarrow 0}
    \int_{-1}^1\frac{dx}{\xi}\,C_n^{3/2}\left(\frac{x}{\xi}\right)
    \mathbb{V}_{\rm NS}^{(0)}\left(\frac{x}{\xi},\frac{y}{\xi}\right) = \frac{(2n+1)!}{2^n(n!)^2}y^{n} P_{\rm NS}^{(0)}(n+1) =\left[\lim_{\xi\rightarrow 0}
   \xi^{n}C_n^{3/2}\left(\frac{y}{\xi}\right)\right]P_{\rm NS}^{(0)}(n+1)\,,

where, again with abuse of notation, the Mellin moments of the splitting
functions are given by (see, *e.g.*, Eq. (4.129) of Ref. (Ellis,
Stirling, and Webber 2011)):

.. math:: P_{\rm NS}^{(0)}(n) = \int_{0}^{1}dx\,x^{n-1}P_{\rm NS}^{(0)}(x) = 2C_F\left[-\frac12+\frac{1}{n(n+1)}-2\sum_{k=2}^{n}\frac{1}{k}\right]\,.

Finally, Eq. :eq:`eq:confmom` is fulfilled and the
evolution kernel can be directly read off:

.. math::
   :label: eq:forwardlimitconfmom
    \lim_{\xi\rightarrow 0} \mathcal{V}_{n}^{(0)}(\xi) = P_{\rm
      NS}^{(0)}(n+1)= 2C_F\left[-\frac12+\frac{1}{(n+1)(n+2)}-2\sum_{k=2}^{n+1}\frac{1}{k}\right]

On the other hand, it is well known that the one-loop ERBL kernel
:math:`V_{\rm NS}^{(0)}` obeys the relation:

.. math::
   :label: eq:confmomerbl
     \int_{-1}^{1}dx\,C_n^{3/2}\left(x\right) V_{\rm
       NS}^{(0)}\left(x,y\right) =V_{n}^{(0)} C_n^{3/2}\left(y\right)\,.

Evidently, :math:`V_{n}^{(0)}` is the evolution kernel in the
:math:`\xi\rightarrow 1` limit and its expression can be read off from
Eq. (2.19) of Ref. (Lepage and Brodsky 1980) where, to match our
notation, we have to flip the sign and remove the
:math:`\beta`-function:

.. math:: \lim_{\xi\rightarrow 1}\mathcal{V}_{n}(\xi)= V_{n}^{(0)}=2C_F\left[-\frac12+\frac{1}{(n+1)(n+2)}-2\sum_{k=2}^{n+1}\frac{1}{k}\right]\,.

Remarkably, this is identical to
Eq. :eq:`eq:forwardlimitconfmom`. This
finally means that Eq. :eq:`eq:confmom` is fulfilled in
both the :math:`\xi\rightarrow 0` and :math:`\xi\rightarrow 1` limits
and that in these limits the evolution kernels are equal. In other
words, conformal moments of forward distributions and distribution
amplitudes at leading order evolve multiplicatively with the same
anomalous dimension. This suggests that the kernels
:math:`\mathcal{V}_{n}` may be constant in :math:`\xi` and equal to
Eq. :eq:`eq:forwardlimitconfmom`.

To prove this we need to compute
Eq. :eq:`eq:confmomgpd` for a generic value of
:math:`\xi`. To do so, we use
Eq. :eq:`eq:evkernelERBLgenxi` for
:math:`\mathbb{V}_{\rm NS}^{(0)}` that allows us to decompose
Eq. :eq:`eq:confmomgpd` as follows:

.. math::
   :label: eq:GPDevkernelalaERBL
   \begin{array}{rcl}
     \displaystyle
     \int_{-1}^1\frac{dx}{\xi}\,C_n^{3/2}\left(\frac{x}{\xi}\right)
     \mathbb{V}_{\rm NS}^{(0)}\left(\frac{x}{\xi},\frac{y}{\xi}\right)&=&\displaystyle 3C_F C_n^{3/2}\left(\frac{y}{\xi}\right) \\
   \\
   &-&\displaystyle C_F\int_{\xi}^{y}dx\left[\frac{x+\xi}{\xi(y-\xi)}C_n^{3/2}\left(\frac{x}{\xi}\right) -2\frac{C_n^{3/2}\left(x/\xi\right)-C_n^{3/2}\left(y/\xi\right)}{y-x}\right]\\
   \\
   &+&\displaystyle C_F\int_{-\xi}^{y}dx\left[\frac{x-\xi}{\xi(y+\xi)}C_n^{3/2}\left(\frac{x}{\xi}\right)+2\frac{C_n^{3/2}\left(x/\xi\right)-C_n^{3/2}\left(y/\xi\right)}{y-x}\right]\,.
   \end{array}

\ Of course, in the limits :math:`\xi\rightarrow 0` and
:math:`\xi\rightarrow 1` this expression produces the Gegenbauer
polynomials in the respective limits times :math:`V_n^{0}`. The
calculation for *all* :math:`n` and for a generic value of :math:`\xi`
is not straight forward. To start with, we compute the integral above
for :math:`n=0` where we expected to obtain zero because
:math:`V^{(0)}=0`. Since :math:`C_0^{\alpha}(z)=1`, the zero-th
conformal moment evaluates to:

.. math::
   \begin{array}{rcl}
   \displaystyle
     \int_{-1}^1\frac{dx}{\xi}\,C_n^{3/2}\left(\frac{x}{\xi}\right)
     \mathbb{V}_{\rm NS}^{(0)}\left(\frac{x}{\xi},\frac{y}{\xi}\right)  &=& \displaystyle
     C_F\left[3 - \int_{\xi}^{y}dx \frac{x+\xi}{\xi(y-\xi)}
   +\int_{-\xi}^{y}dx\frac{x-\xi}{\xi(y+\xi)}\right]\\
   \\
   &=& \displaystyle C_F\left[3 -\frac{y+\xi}{2\xi}-1+\frac{y-\xi}{2\xi}-1\right]=0\,,
   \end{array}

providing further evidence that indeed Gegenbauer polynomials do
diagonalise GPD leading-order evolution with evolution kernels
:math:`\mathcal{V}_n^{(0)}(\xi)=V_n^{(0)}` independent of :math:`\xi`.

Now let us consider the case :math:`n=1`. In the case the Gegenbauer
polynomial is:

.. math:: C_1^{3/2}\left(\frac{x}{\xi}\right)=\frac{3x}{\xi}\,.

Since:

.. math:: V_1^{(0)} = -\frac{8C_F}{3}\,,

we expect:

.. math::
   \int_{-1}^1\frac{dx}{\xi}\,C_1^{3/2}\left(\frac{x}{\xi}\right)
     \mathbb{V}_{\rm NS}^{(0)}\left(\frac{x}{\xi},\frac{y}{\xi}\right) \mathop{=}^{?} -\frac{8C_F y}{\xi}\,.

Using Eq. :eq:`eq:GPDevkernelalaERBL` with
:math:`n=1`, we find:

.. math::
   \begin{array}{rcl}
     \displaystyle
     \int_{-1}^1\frac{dx}{\xi}\,C_1^{3/2}\left(\frac{x}{\xi}\right)
     \mathbb{V}_{\rm
     NS}^{(0)}\left(\frac{x}{\xi},\frac{y}{\xi}\right)&=&\displaystyle
                                                          \frac{3C_F}{\xi}\left[3y
                                                          -\int_{\xi}^{y}dx\left[\frac{x(x+\xi)}{\xi(y-\xi)}
                                                          +2\right]+\int_{-\xi}^{y}dx\left[\frac{x(x-\xi)}{\xi(y+\xi)}-2\right]\right]\\
   \\
   &=&\displaystyle \frac{3C_F}{\xi}\left[3y-4y-\frac{5}{3}y\right] = -\frac{8C_Fy}{\xi}\,,
   \end{array}

as expected. Despite an all-:math:`n` proof would desirable, these
calculations provide compelling evidence that indeed conformal moments
evolve multiplicatively as:

.. math:: \frac{d \mathcal{C}_n(\xi)}{d\ln\mu^2} =\frac{\alpha_s(\mu)}{4\pi}{V}_{n}^{(0)}\mathcal{C}_n(\xi)\,.

The solution to this equation is simply:

.. math::
   \mathcal{C}_n(\xi,\mu)
   =\exp\left[\frac{V_{n}^{(0)}}{4\pi}\int_{\mu_0}^\mu d\ln\mu^2\alpha_s(\mu)\right]\mathcal{C}_n(\xi,\mu_0)\,,

in which remarkably the evolution kernel does not depend on :math:`\xi`.
This in turns means that the :math:`\xi` dependence of conformal moments
is entirely encoded in the initial scale moments and that evolution
leave it unchanged.

The question that remains to be answered is whether the knowledge of the
conformal moments of a GPD is enough to reconstruct its :math:`x`
dependence.

.. container::
   :name: eq:sumruleNS

   .. rubric:: Sum rules
      :name: eq:sumruleNS

A question arises: does the fact that the GPD anomalous dimension
(*cfr.* Eq. :eq:`eq:nonzeroplusGPD`) does not
admit a fully :math:`+`-prescribed form violate any conservation law? To
answer this question, we notice that the fact that the non-singlet DGLAP
anomalous dimension integrates to zero (see
Eq. :eq:`eq:dglapsplitting`), and thus admits a
:math:`+`-prescribed form, derives from the conservation of the total
number of quarks minus anti-quarks (valence sum rule):

.. math::
   :label: eq:valsumrulef
   \int_0^1dx\, f^-(x,0) = \mbox{constant}\,.

\ Taking the derivative of this equation w.r.t. :math:`\ln\mu^2` and
using the DGLAP equation gives:

.. math::
   :label: eq:valsumrule
   \begin{array}{rcl}
   0&=&\displaystyle\int_0^1dx\int_x^1\frac{dy}{y}
        \mathcal{P}_1^-(y,0)f^-\left(\frac{x}{y},0\right) =
        \int_0^1dy\,
        \mathcal{P}_1^-(y,0)\int_0^y\frac{dx}{y}\,f^-\left(\frac{x}{y},0\right)\\
   \\
   &=&\displaystyle \int_0^1dy\,
        \mathcal{P}_1^-(y,0)\int_0^1dz\,f^-\left(z,0\right) =
       \mbox{constant}\times\int_0^1dy\,\mathcal{P}_1^-(y,0)\quad\Leftrightarrow
       \quad \int_0^1dy\,\mathcal{P}_1^-(y,0) = 0\,.
   \end{array}

This clearly justifies the requirement for :math:`\mathcal{P}_1^-(y,0)`
to be fully :math:`+`-prescribed.

One may try to apply the same argument to GPDs. In this case the valence
sum rule generalises in:

.. math::
   :label: eq:sumrule
   \int_0^1 dx\,f^-(x,\xi) = F\,,

\ where :math:`F` is independent of :math:`\mu` and :math:`\xi` but may
(and does) depend on the momentum transfer :math:`t`. :math:`F` is
usually referred to as form factor. One should now take the derivative
w.r.t. :math:`\ln\mu^2` and use Eq. :eq:`eq:jitot` but in
doing this one needs to take into account that :math:`\kappa=\xi/x`:

.. math::
   \begin{array}{rcl}
     0&=&\displaystyle \int_0^1
          dx\int_0^1\frac{dy}{y}\left[\theta(y-x)\mathcal{P}_1^{-,(0)}\left(\frac{x}{y},\frac{\xi}{x}\right)+\theta(\xi-x)
          \mathcal{P}_2^{-,(0)}\left(\frac{x}{y},\frac{\xi}{x}\right)\right]f^-\left(y,\xi\right)\\
   \\
   &=&\displaystyle \int_0^1 dy\,f^-\left(y,\xi\right) \left[\int_0^y\frac{dx}{y}\,\mathcal{P}_1^{-,(0)}\left(\frac{x}{y},\frac{\xi}{x}\right)+
          \int_0^\xi \frac{dx}{y}\, \mathcal{P}_2^{-,(0)}\left(\frac{x}{y},\frac{\xi}{x}\right)\right]\\
   \\
   &=&\displaystyle \int_0^1 dy\,f^-\left(y,\xi\right) \left[\int_0^1dz\,\mathcal{P}_1^{-,(0)}\left(z,\frac{\xi}{yz}\right)+
          \int_0^{\xi/y} dz\, \mathcal{P}_2^{-,(0)}\left(z,\frac{\xi}{yz}\right)\right]\,,
   \end{array}

In order for this relation to be identically true, it is necessary that:

.. math::
   :label: eq:valsumruleGPD
     \int_0^1dz\,\mathcal{P}_1^{-,(0)}\left(z,\frac{\xi}{yz}\right)+
     \int_0^{\xi/y} dz\, \mathcal{P}_2^{-,(0)}\left(z,\frac{\xi}{yz}\right)=0\,.

Notice that for :math:`\xi\rightarrow 0`, the equality above reduces to
Eq. :eq:`eq:valsumrule`. It is interesting to verify
Eq. :eq:`eq:valsumruleGPD` plugging in the
explicit expressions for :math:`\mathcal{P}_1^{-,(0)}`,
Eq. :eq:`eq:p1minus0`, and
:math:`\mathcal{P}_2^{-,(0)}`, Eq. :eq:`eq:p2minus0`.
One finds:

.. math:: \int_0^1dz\,\mathcal{P}_1^{-,(0)}\left(z,\frac{\xi}{yz}\right)=2C_F\left[-\frac{3}{2}\frac{\xi^2}{\xi^2-y^2}-\ln\left(\left|1-\frac{\xi^2}{y^2}\right|\right)\right]\,,

that correctly tends to zero as :math:`\xi\rightarrow 0`, and:

.. math:: \int_0^{\xi/y} dz\, \mathcal{P}_2^{-,(0)}\left(z,\frac{\xi}{yz}\right)=2C_F\left[\frac{3}{2}\frac{\xi^2}{\xi^2-y^2}+\ln\left(\left|1-\frac{\xi^2}{y^2}\right|\right)\right]\,,

such that Eq. :eq:`eq:valsumruleGPD` is fulfilled.
Despite Eq. :eq:`eq:valsumruleGPD` has been
explicitly proved at one-loop, the same relation must hold order by
order in perturbation theory.

It is important to notice that the constraint on the non-singlet GPD
anomalous dimensions deriving from the valence sum rule, and resulting
in Eq. :eq:`eq:valsumruleGPD`, does not take the
form of a :math:`+`-prescription,
Eq. :eq:`eq:pludistributionn`. A further proof
can be given by considering the non-singlet GPD evolution equation given
in Eq. :eq:`eq:eveq` (see also Eq. (99) of Ref. (Diehl
2003)). The independence of the form factor from :math:`\mu` immediately
leads to:

.. math::
   \int_{-1}^1dx'
   f(x',\xi)\int_{-1}^1dx\left[\hat{V}_{\rm NS}\left(\frac{x}{\xi},\frac{x'}{\xi}\right)\right]_+ = 0\,,

where :math:`\hat{V}_{\rm NS}` is nothing but :math:`V _{\rm NS}`
stripped of the supposedly global :math:`+`-prescription. Assuming
:math:`\xi` positive and different from zero (:math:`\xi>0`), a simple
change of variables gives:

.. math::
   \int_{-1/\xi}^{1/\xi}dy'
     f(\xi y',\xi)\int_{-1/\xi}^{1/\xi}dy\left[\hat{V}_{\rm NS}\left(y,y'\right)\right]_+ = 0\,.

Given the fact that the bounds of the inner integral are not :math:`-1`
and :math:`1`, the effect of the :math:`+`-prescription as given in
Eq. :eq:`eq:pludistributionn` cannot give zero.
This prevents the above equation to be identically fulfilled violating
polynomiality of GPDs. We can thus conclude that :math:`V _{\rm NS}`
*cannot* be written as a fully :math:`+`-prescribed function.

.. container::
   :name: the-singlet-sector

   .. rubric:: The singlet sector
      :name: the-singlet-sector

Having ascertained that the evolution equation from Ref. (Ji 1997) is
well-behaved for the non-singlet distribution :math:`f^-`, we move on to
consider the singlet :math:`f_{\rm S}` and gluon :math:`f_{\rm G}`
distributions. As in the standard DGLAP evolution equation, singlet and
gluon GPDs couple under evolution. Defining :math:`f^+` as the column
vector of singlet and gluon GPDs, the corresponding anomalous dimension
:math:`\mathcal{P}^+` is a matrix in flavour space:

.. math::
   :label: eq:singmat
     \mathcal{P}^+=\begin{pmatrix}
       \mathcal{P}_{\rm SS}& \mathcal{P}_{\rm SG}\\
       \mathcal{P}_{\rm GS}& \mathcal{P}_{\rm GG}
   \end{pmatrix}\,.

\ Following the same procedure discussed above for the non-singlet
distribution :math:`f^-`, the one-loop evolution equation for
:math:`f^+` reads:

.. math::
   :label: eq:jitotsing
     \mu^2\frac{d}{d\mu^2}f^+(x,\xi) = \frac{\alpha_s(\mu)}{4\pi}\int_x^\infty\frac{dy}{y}\mathcal{P}^{+,(0)}(y,\kappa)f^+\left(\frac{x}{y},\xi\right)\,,

with:

.. math:: \mathcal{P}^{+,(0)}(y,\kappa)= \theta(1-y)\mathcal{P}_1^{+,(0)}(y,\kappa)+\theta(\kappa-1)\mathcal{P}_2^{+,(0)}(y,\kappa)\,.

The single splitting function matrices :math:`\mathcal{P}_{1}` and
:math:`\mathcal{P}_2` are derived from the expression Ref. (Ji 1997) as:

.. math::
   \begin{array}{rcl}
   \mathcal{P}_{1,\rm IJ}(y,\kappa) &=& 2P_{\rm IJ}(y,2\kappa y)= 2P_{\rm IJ}'(y,2\kappa y) +
                                   2P_{\rm IJ}'(y,-2\kappa y)\,,\\
   \\
   \mathcal{P}_{2,\rm IJ}(y,\kappa) &=& - 2P_{\rm IJ}'(y,-2\kappa y) - 2P_{\rm IJ}'(-y,2\kappa y)\,,
   \end{array}

with :math:`\rm I,J=S,G`. This leads to:

.. math::
   \left\{\begin{array}{rcl}
   \mathcal{P}_{1,\rm SS}^{(0)}(y,\kappa) &=& \mathcal{P}_{1}^{-,(0)}(y,\kappa)\,,\\
   \\
   \mathcal{P}_{2,\rm SS}^{(0)}(y,\kappa) &=& \displaystyle
                                              2C_F\left[\frac{1+y+\kappa y+\kappa^3y^2}{\kappa(1+y)(1-\kappa^2y^2)}-\left(\frac{1}{1-y}\right)_{++}\right]\,,
   \end{array}\right.

.. math::
   \left\{\begin{array}{rcl}
   \mathcal{P}_{1,\rm SG}^{(0)}(y,\kappa) &=& \displaystyle 4n_f T_R\left[\frac{y^2+(1-y)^2-\kappa^2y^2}{(1-\kappa^2y^2)^2}\right]\,,\\
   \\
   \mathcal{P}_{2,\rm SG}^{(0)}(y,\kappa) &=& \displaystyle 4n_f T_R(1-\kappa)\left[\frac{1-\kappa(\kappa+2)y^2}{\kappa(1-\kappa^2y^2)^2}\right]\,,
   \end{array}\right.

.. math::
   \left\{\begin{array}{rcl}
     \mathcal{P}_{1,\rm GS}^{(0)}(y,\kappa) &=& \displaystyle 2C_F\left[\frac{1+(1-y)^2-\kappa^2y^2}{y(1-\kappa^2y^2)}\right]\,,\\
     \\
     \mathcal{P}_{2,\rm GS}^{(0)}(y,\kappa) &=& \displaystyle-2C_F\frac{(1-\kappa)^2}{\kappa(1-\kappa^2 y^2)}\,,
   \end{array}\right.

.. math::
   \left\{\begin{array}{rcl}
              \mathcal{P}_{1,\rm GG}^{(0)}(y,\kappa) &=& \displaystyle
                                                         4C_A\left[\left(\frac{1}{1-y}\right)_+-\frac{1+\kappa^2y}{1-\kappa^2y^2}+\frac{1}{(1-\kappa^2y^2)^2}\left(\frac{1-y}{y}+y(1-y)\right)\right]\\
              \\
                                                     &+&\displaystyle\delta(1-y)\left[\left(\frac{11C_A-4 n_f
                                                         T_R}{3}\right) -2C_A\ln(|1-\kappa^2|)\right]\,,\\
              \\
              % &=& \displaystyle
              %                                            4C_A\bigg[\left(\frac{1}{1-y}-\frac{\kappa^2y}{1-\kappa^2y^2}\right)_++\frac{1}{(1-\kappa^2y^2)^2}\left(\frac{1-y}{y}+y(1-y)-1+\kappa^2y^2\right)\\
              % \\
              %                                        &+&\displaystyle\delta(1-y)\left(\frac{11C_A-4 n_f T_R}{3}\right)\,,\\
              % \\
              \mathcal{P}_{2,\rm GG}^{(0)}(y,\kappa) &=& \displaystyle
                                                         2C_A\left[\frac{2(1-\kappa)(1+y^2)}{(1-\kappa^2y^2)^2}+\frac{\kappa^2(1+y)}{1-\kappa^2y^2}+\frac{1-\kappa^2}{1-\kappa^2y^2}\left(2-\frac{1}{\kappa}-\frac{1}{1+y}\right)
                                                         -\left(\frac{1}{1-y}\right)_{++}\right]
                                                         \,.
   \end{array}\right.

In all cases, the limit for :math:`\kappa\rightarrow 0` of
:math:`\mathcal{P}_1` reproduces the one-loop DGLAP splitting functions.
In addition, we also notice that all :math:`\mathcal{P}_2`\ ’s are
proportional to :math:`\kappa-1`. Along with the fact that all
:math:`\mathcal{P}_1` are well-behaved at :math:`\kappa=1`:

.. math::
   \begin{array}{rcl}
   \mathcal{P}_{1,\rm SS}^{(0)}(y,1) &=& \displaystyle 2C_F\left\{\left[\frac{1}{1-y}\right]_++\delta(1-y)\left[\frac32-\ln(2)\right]\right\}\,,\\
   \\
   \mathcal{P}_{1,\rm SG}^{(0)}(y,1) &=& \displaystyle \frac{4n_f T_R}{(1+y)^2}\,,\\
   \\
     \mathcal{P}_{1,\rm GS}^{(0)}(y,1) &=& \displaystyle \frac{4C_F}{y(1+y)}\,,\\
     \\
     \mathcal{P}_{1,\rm GG}^{(0)}(y,1) &=& \displaystyle 4C_A\left[\left(\frac{1}{1-y^2}\right)_++\frac{1}{y(1+y)^2}\right]+\delta(1-y)\left(\frac{11C_A-4 n_f T_R }{3}\right)\,,
   \end{array}

\ this guarantees the continuity of GPDs across the point :math:`x=\xi`.

Now, we explicitly verify that the pole at :math:`y=|\kappa|^{-1}` that
affects all splitting functions above cancels between
:math:`\mathcal{P}_{1,\rm IJ}^{(0)}` and
:math:`\mathcal{P}_{2,\rm IJ}^{(0)}`. We find:

.. math::
   \lim_{y\rightarrow \kappa^{-1}} (1-\kappa^2y^2)
   \mathcal{P}_{1,\rm SS}^{(0)}(y,\kappa) = - \lim_{y\rightarrow \kappa^{-1}} (1-\kappa^2y^2)
   \mathcal{P}_{2,\rm SS}^{(0)}(y,\kappa) = -2C_F\frac{1+\kappa}{\kappa}\,,

.. math::
   \lim_{y\rightarrow \kappa^{-1}} (1-\kappa^2y^2)^2
   \mathcal{P}_{1,\rm SG}^{(0)}(y,\kappa) = - \lim_{y\rightarrow \kappa^{-1}} (1-\kappa^2y^2)^2
   \mathcal{P}_{2,\rm SG}^{(0)}(y,\kappa) = \frac{8 n_f T_R(1-\kappa)}{\kappa}\,,

.. math::
   \lim_{y\rightarrow \kappa^{-1}} (1-\kappa^2y^2)
   \mathcal{P}_{1,\rm GS}^{(0)}(y,\kappa) = - \lim_{y\rightarrow \kappa^{-1}} (1-\kappa^2y^2)
   \mathcal{P}_{2,\rm GS}^{(0)}(y,\kappa) = 2C_F\frac{(1-\kappa)^2}{\kappa}\,,

.. math::
   \lim_{y\rightarrow \kappa^{-1}} (1-\kappa^2y^2)^2
   \mathcal{P}_{1,\rm GG}^{(0)}(y,\kappa) = - \lim_{y\rightarrow \kappa^{-1}} (1-\kappa^2y^2)^2
   \mathcal{P}_{2,\rm GG}^{(0)}(y,\kappa) = 4C_A\frac{(\kappa-1)(\kappa^2+1)}{\kappa^2}\,.

These results confirm the cancellation of the pole at
:math:`y=|\kappa|^{-1}` in the integral in the r.h.s. of
Eq. :eq:`eq:jitotsing`.

We now compute the ERBL limit by taking :math:`\kappa\rightarrow 1/x`.
To do so, we use the ERBL-compliant form of the evolution equation:

.. math::
   :label: eq:erblsing
     \mu^2\frac{d}{d\mu^2}\Phi^+(t) =
     \frac{\alpha_s(\mu)}{4\pi}\left[\int_{0}^1 du\,\overline{V}_{\rm S}^{(0)}(t,u) \Phi^+\left(u\right)\right] \,,

with :math:`\Phi^+(t) = f^+(x,1)` and where :math:`u` and :math:`t` are
defined in Eq. :eq:`eq:ERBLvariables`.
:math:`\overline{V}_{\rm S}^{(0)}` is a matrix in flavour space with the
same structure of :math:`\mathcal{P}^+` in
Eq. :eq:`eq:singmat`, whose components can be written in
terms of the :math:`P_{\rm IJ}'` functions as follows:

.. math::
   \begin{array}{rcl}
   V_{\rm SS}^{(0)}(t,u) &=&\displaystyle \theta(u-t)\left[\frac{2}{2u-1}P_{\rm
       SS}'\left(\frac{2t-1}{2u-1},\frac{2}{2u-1}\right)
                             \right]-C_F\delta\left(u-t\right)\int_{0}^1 du' \,\frac{\theta(u'-t)}{u'-t}\\
   \\
   &+&\displaystyle \theta(t-u)\left[-\frac{2}{2u-1}P_{\rm
       SS}'\left(\frac{2t-1}{2u-1},-\frac{2}{2u-1}\right)
                             \right]+C_F\delta\left(t-u\right)\int_{0}^{1}du'\,\frac{\theta(t-u')}{u'-t}\\
   \\
   &+&\displaystyle \frac32C_F\delta\left(u-t\right)\,, \\
   \\
   V_{\rm SG,GS}^{(0)}(t,u) &=&\displaystyle \theta(u-t)\left[\frac{2}{2u-1}P_{\rm
       SG,GS}'\left(\frac{2t-1}{2u-1},\frac{2}{2u-1}\right)
                             \right]\\
   \\
   &+&\displaystyle \theta(t-u)\left[-\frac{2}{2u-1}P_{\rm
       SG,GS}'\left(\frac{2t-1}{2u-1},-\frac{2}{2u-1}\right)
                             \right]\,,\\
   \\
   V_{\rm GG}^{(0)}(t,u) &=&\displaystyle \theta(u-t)\left[\frac{2}{2u-1}P_{\rm
       GG}'\left(\frac{2t-1}{2u-1},\frac{2}{2u-1}\right)
                             \right]-C_A\delta\left(u-t\right)\int_{0}^1 du' \,\frac{\theta(u'-t)}{u'-t}\\
   \\
   &+&\displaystyle \theta(t-u)\left[-\frac{2}{2u-1}P_{\rm
       GG}'\left(\frac{2t-1}{2u-1},-\frac{2}{2u-1}\right)
                             \right]+C_A\delta\left(t-u\right)\int_{0}^{1}du'\,\frac{\theta(t-u')}{u'-t}\\
   \\
   &+&\displaystyle \left(\frac{11C_A-4n_fT_R}{6}\right)\delta\left(u-t\right)\,.
   \end{array}

The explicit expressions read:

.. math::
   \begin{array}{rcl}
     V_{\rm SS}^{(0)}(t,u) &=&\displaystyle V_{\rm NS}^{(0)}(t,u)\,,\\
     \\
     V_{\rm SG}^{(0)}(t,u) &=&\displaystyle 2n_fT_R\left(\frac{2u-1}{2}\right) \left[\theta(u-t)\frac{t}{u}\left(\frac{2t-1}{u}-2\frac{1-t}{1-u}\right)-\theta(t-u)\left(\frac{1-t}{1-u}\right)\left(\frac{1-2t}{1-u}-2\frac{t}{u}\right)\right]\,,\\
     \\
     V_{\rm GS}^{(0)}(t,u) &=&\displaystyle C_F \left(\frac{2}{2t-1}\right)
                               \left[\theta(u-t)\left(2t-\frac{t^2}{u}\right)
                               -\theta(t-u) \left(2(1-t)-\frac{(1-t)^2}{1-u}\right)\right]\,,\\
     \\
     V_{\rm GG}^{(0)}(t,u) &=&\displaystyle \mbox{This is a little
                               convoluted\dots leave it for when I feel
                               like doing the calculation}\,.
   \end{array}

.. container::
   :name: eq:sumruleSG

   .. rubric:: Sum rules
      :name: eq:sumruleSG

As done above for the non-singlet splitting function, here we exploit
polynomiality and the GPDs and the momentum sum rule (MSR) of PDFs to
check the correctness of the singlet splitting functions. Before
considering GPDs, let us derive the constrain that the MSR has on the
forward splitting functions. The MSR states that the integral over the
longitudinal momentum fractions :math:`x` of the sum of all PDFs
weighted by :math:`x` represents the momentum fraction of all partons in
the hadron and thus is to be equal to one at all scales. The respective
formula reads:

.. math:: \int_0^1dx\,x\left[f_{\rm S}(x,0)+f_{\rm G}(x,0)\right] = 1\,.

\ We can now that the derivative w.r.t. :math:`\ln\mu^2` of both side of
the the equation above and, using the evolution equation we find:

.. math::
   \int_0^1dy\,y\left[\mathcal{P}_{1,\rm SS}(y,0)+\mathcal{P}_{1,\rm
       GS}(y,0)\right]\left[\int_0^1dz\,z f_{\rm S}(z,0)\right] + \int_0^1dy\,y\left[\mathcal{P}_{1,\rm SG}(y,0)+\mathcal{P}_{1,\rm
       GG}(y,0)\right]\left[\int_0^1dz\,z f_{\rm G}(z,0)\right] = 0\,.

In order for this equation to be identically fulfilled for any pair of
distributions :math:`f_{\rm S}` and :math:`f_{\rm G}`, one needs to
have:

.. math::
   :label: eq:MSRforward
   \begin{array}{rcl}
   \displaystyle \int_0^1dy\,y\left[\mathcal{P}_{1,\rm SS}(y,0)+\mathcal{P}_{1,\rm
       GS}(y,0)\right] &=& 0 \,,\\
   \\
    \displaystyle \int_0^1dy\,y\left[\mathcal{P}_{1,\rm SG}(y,0)+\mathcal{P}_{1,\rm
       GG}(y,0)\right] &=& 0\,.
   \end{array}

Effectively, these equations are fulfilled by the DGLAP splitting
functions. We now wish to extend this argument to GPDs to find similar
relations for the splitting functions that apply to :math:`\xi\neq 0`.
To do so, we use polynomiality of GPDs that for both singlet and gluon
reads:

.. math:: \int_0^1dx\,xf_{\rm S(G)}(x,\xi) = A_{q(g)}+\xi^2 D_{q(g)}\,,

where the coefficients :math:`A` and :math:`D` generally depend on the
scale :math:`\mu`. However, it is well known that helicity-conserving
(:math:`H`) and helicity-flip (:math:`E`) GPDs have the same D-term but
with opposite sign such that, if we define :math:`f_{\rm S(G)}` as the
sum of :math:`H` and :math:`E` the :math:`\xi^2` term cancels and one is
left with:

.. math:: \int_0^1dx\,xf_{\rm S(G)}(x,\xi) = A_{q(g)}\,.

Therefore the sum of :math:`f_{\rm S}` and :math:`f_{\rm G}` gives:

.. math:: \int_0^1dx\,x\left[f_{\rm S}(x,\xi) + f_{\rm G}(x,\xi)\right]= A\,,

with :math:`A=A_q+A_g`. Importantly :math:`A` is an observable and thus
independent of the scale :math:`\mu`. In addition, :math:`H` and
:math:`E` obey the same evolution equations. This allows us to take the
derivative with respect to :math:`\ln\mu^2` of both sides of the
equation above and use the evolution equations. Following the same steps
as in Sect. `2.3.3 <#eq:sumruleNS>`__ leads to the following equalities:

.. math::
   :label: eq:MSRoffforward
   \begin{array}{rcl}
   \displaystyle 
     \int_0^1dz\,z\left[\mathcal{P}_{1,\rm SS}\left(z,\frac{\xi}{yz}\right)+\mathcal{P}_{1,\rm GS}\left(z,\frac{\xi}{yz}\right)\right]+
     \int_0^{\xi/y} dz\, z\left[\mathcal{P}_{2,\rm SS}\left(z,\frac{\xi}{yz}\right)+\mathcal{P}_{2,\rm GS}\left(z,\frac{\xi}{yz}\right)\right]&=&0\,,\\
   \\
    \displaystyle \int_0^1dz\,z\left[\mathcal{P}_{1,\rm SG}\left(z,\frac{\xi}{yz}\right)+\mathcal{P}_{1,\rm GG}\left(z,\frac{\xi}{yz}\right)\right]+
     \int_0^{\xi/y} dz\, z\left[\mathcal{P}_{2,\rm SG}\left(z,\frac{\xi}{yz}\right)+\mathcal{P}_{2,\rm GG}\left(z,\frac{\xi}{yz}\right)\right]&=&0\,,
   \end{array}

that hold at each single order in :math:`\alpha_s` and automatically
reduce to Eq. :eq:`eq:MSRforward` for
:math:`\xi\rightarrow 0`. It is now important to verify that the
one-loop splitting functions written above obey these relations. We
find:

.. math::
   \int_0^1dz\,z\left[\mathcal{P}^{(0)}_{1,\rm
       SS}\left(z,\frac{\xi}{yz}\right)+\mathcal{P}^{(0)}_{1,\rm
       GS}\left(z,\frac{\xi}{yz}\right)\right] = 2C_F\left[-\frac12\frac{\xi^2}{y^2-\xi^2}-\ln\left(\left|1-\frac{\xi^2}{y^2}\right|\right)\right]\,,

and:

.. math:: \int_0^{\xi/y} dz\, z\left[\mathcal{P}^{(0)}_{2,\rm SS}\left(z,\frac{\xi}{yz}\right)+\mathcal{P}^{(0)}_{2,\rm GS}\left(z,\frac{\xi}{yz}\right)\right]= 2C_F\left[\frac12\frac{\xi^2}{y^2-\xi^2}+\ln\left(\left|1-\frac{\xi^2}{y^2}\right|\right)\right]\,,

that evidently cancel so that the first of
Eq. :eq:`eq:MSRoffforward` is satisfied. In
addition, they both tend to zero as :math:`\xi\rightarrow 0` as required
by the first of Eq. :eq:`eq:MSRforward`. Then we
find:

.. math::
   \begin{array}{rcl}
   \displaystyle\int_0^1dz\,z\left[\mathcal{P}^{(0)}_{1,\rm
       SG}\left(z,\frac{\xi}{yz}\right)+\mathcal{P}^{(0)}_{1,\rm
       GG}\left(z,\frac{\xi}{yz}\right)\right] &=&\displaystyle
     \frac{y^2\xi^2}{3(y^2-\xi^2)^2}\left[C_A\left(\frac{11\xi^2}{y^2}-4\right)+2n_fT_R\left(1-\frac{2\xi^2}{y^2}\right)\right]\\
   \\
   &-&\displaystyle 2C_A\ln\left(\left|\frac{y^2-\xi^2}{y^2}\right|\right)\,,
   \end{array}

and:

.. math::
   \begin{array}{rcl}
   \displaystyle \int_0^{\xi/y} dz\, z\left[\mathcal{P}^{(0)}_{2,\rm
     SG}\left(z,\frac{\xi}{yz}\right)+\mathcal{P}^{(0)}_{2,\rm
     GG}\left(z,\frac{\xi}{yz}\right)\right]&=& \displaystyle
                                                -\frac{y^2\xi^2}{3(y^2-\xi^2)^2}\left[C_A\left(\frac{11\xi^2}{y^2}-4\right)+2n_fT_R\left(1-\frac{2\xi^2}{y^2}\right)\right]\\
   \\
   &+&\displaystyle 2C_A\ln\left(\left|\frac{y^2-\xi^2}{y^2}\right|\right)\,,
   \end{array}

so that also the second of
Eq. :eq:`eq:MSRoffforward` is satisfied. Again,
their limit for :math:`\xi\rightarrow 0` tends to zero that is
consistent with the second of
Eq. :eq:`eq:MSRforward`.

Compton form factor
===================

In this section, we discuss the numerical implementation of the Compton
form factors (CFFs) for double deeply-virtual Compton scattering
(DDVCS). The CFFs can be regarded as the GPD counterpart of the
deep-inelastic scattering (DIS) structure function. The relation is
easily established by observing that the diagrams necessary to compute
the DIS structure functions at a given order are exactly the same
diagrams of DDVCS in the forward limit with a cut on the final state. It
is well-known that cut and uncut digrams are related through dispersion
relation. Specifically, (:math:`2\pi i` times) a given cut diagram with
an external on-shell momentum :math:`p` is equal to the discontinuity of
the uncut diagram in the forward limit along the branch cut on the real
axis determined by :math:`p^2\geq m^2`, where :math:`m` is the rest mass
of the external particle. For the CCF :math:`\mathcal{H}`, this
translates into the following relation with the DIS structure function
:math:`F_1`:( [8]_)

.. math::
   \lim_{\xi\rightarrow 0}
   \left[\mathcal{H}(x+i\varepsilon,\xi,t,Q)-\mathcal{H}(x-i\varepsilon,\xi,t,Q)\right]
   = 4\pi i F_1(x,Q)\,.

\ Analyticity of the amplitude implies that
:math:`\mathcal{H}(x-i\varepsilon,\xi,t,Q) =
\mathcal{H}^*(x+i\varepsilon,\xi,t,Q)` which finally gives:

.. math::
   :label: eq:opticalth
     \frac{1}{2\pi} \lim_{\xi\rightarrow 0} \mbox{Im}\left[\mathcal{H}(x,\xi,t,Q)\right] = F_1(x,Q)\,.

An analogous relation exists between the CFF
:math:`\widetilde{\mathcal{H}}` and the polarised structure function
:math:`g_5` as well as between the longitudinal projections
:math:`\mathcal{H}_L` and :math:`\widetilde{\mathcal{H}}_L`, and
:math:`F_L` and :math:`g_L`, respectively. These relation can be used to
check that coefficient functions for the CCFs reproduce the well-know
DIS coefficient functions in the forward limit.

Factorisation for DDVCS allows one to write CFFs as convolution between
short-distance partonic cross sections (or coefficient functions) and
nucleon GPDs. For :math:`\mathcal{H}` factorisation reads:

.. math::
   :label: eq:CCFfact
     \mathcal{H}(x,\xi,t,Q) = \sum_{q}e_q^2\int_{-1}^{1}\frac{dy}{y}\left[C_q\left(\frac{x}{y},\frac{\xi}{x},Q\right)F_q(y,\xi,t,Q)+C_g\left(\frac{x}{y},\frac{\xi}{x},Q\right)F_g(y,\xi,t,Q)\right]

where :math:`e_q^2` is the electric charge of the flavour
:math:`q`\ ( [9]_) with the sum running over the *active* flavours (and
not over the anti-flavours). :math:`F_q` and :math:`F_g` are the quark
:math:`q` and gluon GPDs, respectively, and :math:`C_q` and :math:`C_g`
the respective coefficient functions that admit a perturbative expansion
in powers of :math:`\alpha_s`:

.. math::
   :label: eq:cffalstructfuncERBL
   C_q\left(y,\kappa,Q\right) = \sum_{n=0}^{\infty}\left(\frac{\alpha_s(Q)}{4\pi}\right)^nC_q^{(n)}\left(y,\kappa\right)\quad\mbox{and}\quad C_g\left(y,\kappa,Q\right) = \sum_{n=1}^{\infty}\left(\frac{\alpha_s(Q)}{4\pi}\right)^nC_g^{(n)}\left(y,\kappa\right)\,.

Currently, corrections up to :math:`n=2` are known (Mankiewicz et al.
1998; Ji and Osborne 1998; Pire, Szymanowski, and Wagner 2011; Braun et
al. 2020). For implementation purposes, it is convenient to rewrite
Eq. :eq:`eq:CCFfact` as follows:

.. math::
   :label: eq:cffalstructfunc
   \begin{array}{rcl}
     x\mathcal{H}(x,\xi,t,Q) &=&\displaystyle
                                \sum_{q}e_q^2\int_{x}^{\infty}dy\bigg[C_q\left(y,\kappa,Q\right) \frac{x}{y} F_q\left(\frac{x}{y},\xi,t,Q\right)
                                + C_q\left(y,-\kappa,Q\right) \frac{x}{y} F_{\overline{q}}\left(\frac{x}{y},\xi,t,Q\right)\\
     \\
                            &+& \displaystyle  \left[C_g\left(y,\kappa,Q\right)+C_g\left(y,-\kappa,Q\right)\right]\frac{x}{y} F_g\left(\frac{x}{y},\xi,t,Q\right)\bigg]\,,
   \end{array}

with :math:`\kappa=\xi/x` and where we have exploited the symmetries
:math:`C_{q,g}\left(y,\kappa,Q\right)=C_{q,g}\left(-y,-\kappa,Q\right)`
and :math:`F_g(y,\xi,t,Q)=-F_g(-y,\xi,t,Q)`, and defined
:math:`F_{\overline{q}}(y,\xi,t,Q)=-F_q(-y,\xi,t,Q)`.

At first order (:math:`n=0`), only the quark channel contributes and the
respective coefficient function is:

.. math:: C_q^{(0)}\left(y,\kappa\right) = -\frac{1}{1-y+i\varepsilon}-\frac{1}{1+y-i\varepsilon}\,.

In order to get rid of the :math:`i\varepsilon` that complicates the
implementation, we use the following identity:

.. math::
   :label: eq:princval
   \frac{1}{X \mp i\varepsilon} = {\rm PV}\frac{1}{X} \pm \pi i\delta(X)\,

where PV stands for principal value, so that:

.. math:: C_q^{(0)}\left(y,\kappa\right) = {\rm PV}\left[-\frac{1}{1- y}-\frac{1}{1+ y}\right] + \pi i\left[\delta(1-y) - \delta(1+y)\right]\,.

Plugging this into Eq. :eq:`eq:cffalstructfunc`,
defining :math:`\overline{F}_q^{+}(y,\xi,t,Q) =
y(F_q(y,\xi,t,Q)+F_{\overline{q}}(y,\xi,t,Q))`, and understanding the
principal value gives:

.. math::
   :label: eq:cffalstructfunc
   \begin{array}{rcl}
     x\mathcal{H}^{\rm LO}(x,\xi,t,Q) &=&\displaystyle
                                          \sum_{q}e_q^2\int_{x}^{\infty}dy\left[-\left(\frac{1}{1-
                                          y}\right)_{++}-\frac{1}{1+ y}+\pi i\delta(1-y)\right]
                                          \overline{F}_q^{+}\left(\frac{x}{y},\xi,t,Q\right)\,,
   \end{array}

where we have also introduced the :math:`++`-prescription defined in
Eq. :eq:`eq:plusplusdist` that makes the integrand
finite over the integration range making it numerically computable. It
is intersting to observe that applying
Eq. :eq:`eq:opticalth` to the leading-order CFF gives
the expected result:

.. math::
   \frac{1}{2\pi} \lim_{\xi\rightarrow 0} \mbox{Im}\left[x\mathcal{H}^{\rm
         LO}(x,\xi,t,Q)\right] = \frac12\sum_{q}e_q^2xf_q^{+}(x,Q) =
     xF_1^{\rm LO}(x,Q)\,,

where :math:`f_q^+` are the forward PDFs. The particular structure of
Eq. :eq:`eq:cffalstructfunc` allows one to
relate real and imaginary part of the CFF by observing that:

.. math::
   \begin{array}{rcl}
   \displaystyle \mbox{Im}\left[x\mathcal{H}^{\rm
         LO}(x,\xi,t,Q)\right] &=& \displaystyle  \pi \sum_{q}e_q^2
     \overline{F}_q^{+}\left(x,\xi,t,Q\right)\,,\\
   \\
   \displaystyle \mbox{Re}\left[x\mathcal{H}^{\rm
         LO}(x,\xi,t,Q)\right] &=& \displaystyle \int_{x}^{\infty}dy\left[-\left(\frac{1}{1-
                                          y}\right)_{++}-\frac{1}{1+ y}\right]
                                          \sum_{q}e_q^2\overline{F}_q^{+}\left(\frac{x}{y},\xi,t,Q\right)\,.
   \end{array}

Therefore, one can plug the first into the second obtaining:

.. math::
   \mbox{Re}\left[x\mathcal{H}^{\rm
         LO}(x,\xi,t,Q)\right] =\frac{1}{\pi}\int_{x}^{\infty}dy\left[-\left(\frac{1}{1-
                                          y}\right)_{++}-\frac{1}{1+ y}\right]\mbox{Im}\left[x\mathcal{H}^{\rm
         LO}\left(\frac{x}{y},\xi,t,Q\right)\right]\,.

This relation, often referred to as *dispersion* relation, is strictly
valid only at leading order and is a useful tool when extracting GPDs
(or CFFs) from experimental data.

We can now move to consider the :math:`\mathcal{O}(\alpha_s)`
corrections to the CCF :math:`\mathcal{H}`. At this order, both the
quark and the gluon channels contribute. We take the analytic
expressions from Ref. (Pire, Szymanowski, and Wagner 2011). For the
quark channel we use the finite part of Eq. (46) and for the gluon the
finite part of Eq. (47). Since we are assuming :math:`x_{\rm B}>0`
(:math:`=x` in our notation), we need to take :math:`Q^2>0`. In
addition, our definition of coefficient functions in
Eq. :eq:`eq:cffalstructfuncERBL` is such
that the expressions of Ref. (Pire, Szymanowski, and Wagner 2011) have
to be identified with
:math:`(1/x)C_{q,g}^{(1)}(x_{\rm B}/x,\xi/x_{\rm B})`. This gives:

.. math::
   :label: eq:DDVCSceofffuncsNLO
   \begin{array}{rcl}
   C_{q}^{(1)}(y,\kappa) &=& \displaystyle -\frac{C_F}{1-y+i\varepsilon}\bigg\{ -
                             9 +
                             3\frac{1-\kappa^2y^2-2y}{1-\kappa^2y^2}\ln\left(-\frac{1-y}{y}-i\varepsilon\right)+\frac{1+y^2 -2\kappa^2y^2}{1-\kappa^2y^2} \ln^2\left(-\frac{1-y}{y}-i\varepsilon\right)\\
   \\
   &+&\displaystyle 6\frac{\kappa(1-\kappa)y^2}{1-\kappa^2
       y^2}\ln(1-\kappa-i\varepsilon)+
   \frac{(1-\kappa)(1-y-2\kappa^2 y^2)}{\kappa(1-\kappa^2
       y^2)}\ln^2(1-\kappa-i\varepsilon)\bigg\}\\
   \\
   &+& (y\rightarrow -y, \kappa\rightarrow -\kappa) \,,\\
   \\
   C_{g}^{(1)}(y,\kappa) &=& \displaystyle T_R \left(-\frac{1}{1-y+i\varepsilon}-\frac{1}{1+y-i\varepsilon}\right)\\
   \\
   &\times&\displaystyle \bigg\{ 
                             2\ln\left(-\frac{1-y}{y}-i\varepsilon\right)
                             +\frac{1-2y+2y^2-\kappa^2 y^2}{1-\kappa^2
                             y^2}\ln^2\left(-\frac{1-y}{y}-i\varepsilon\right)\\
   \\
   &+&\displaystyle
       \frac{2(1-\kappa)}{\kappa}\ln(1-\kappa-i\varepsilon)+\frac{(1-\kappa)(1-2\kappa
       y^2-\kappa^2 y^2)}{\kappa(1-\kappa^2
       y^2)}\ln^2(1-\kappa-i\varepsilon)
   \bigg\}\\
   \\
   &+& (y\rightarrow -y, \kappa\rightarrow -\kappa) \,.
   \end{array}

\ Now, we need to split the expressions above into real and imaginary
parts. To do so, we use Eq. :eq:`eq:princval` and:

.. math::
   \ln\left(-\frac{1-y}{y}-i\varepsilon\right) =
     \ln\left(\left|\frac{1-y}{y}\right|\right)-i\pi\theta(1-y)\,.

The presence of the :math:`\theta`-function in front of :math:`i\pi` is
a consequence of the fact that that term only arises if the argument of
the logarithm is negative that only happens if :math:`y<1`. It follows
that:

.. math::
   \ln^2\left(-\frac{1-y}{y}-i\varepsilon\right) =
     \ln^2\left(\left|\frac{1-y}{y}\right|\right)-\pi^2\theta(1-y)-2\pi i\theta(1-y) \ln\left(\frac{1-y}{y}\right)\,,

Similarly:

.. math:: \ln(1-\kappa-i\varepsilon) = \ln(|1-\kappa|)-i\pi \theta(\kappa-1)\,.

Before computing the full general expressions for real and imaginary
parts, it is instructive to compute the imaginary part in the forward
limit to check that Eq. :eq:`eq:opticalth` is
fulfilled also at :math:`\mathcal{O}(\alpha_s)`. Let us start with the
quark coefficient function. The forward limit is obtained by taking
:math:`\kappa\rightarrow 0`. By doing so, we observe that the terms
:math:`(y\rightarrow -y, \kappa\rightarrow -\kappa)` in
Eq. :eq:`eq:DDVCSceofffuncsNLO` do not
contribute to the imaginary part. The reason is that the integral in
Eq. :eq:`eq:cffalstructfunc` only covers
positive values of :math:`y` where those terms are regular, so that the
:math:`i\varepsilon` prescription has no effect and no imaginary terms
are generated. The limit produces:( [10]_)

.. math::
   \begin{array}{rcl}
     \displaystyle \frac{1}{2\pi}\lim_{\kappa\rightarrow 0}{\rm
     Im}\left[C_{q}^{(1)}(y,\kappa)\right] &=& \displaystyle C_F\theta(1-y)\left\{
                                               2\frac{\ln\left(1-y\right)}{1-y} -\frac{3}{2}\frac{1}{1-y} -(1+y) \ln\left(1-y\right)-\frac{1+y^2}{1-y}\ln\left(y\right)+3\right\}\\
     \\
                                           &-& \displaystyle C_F\delta(1-y)\left\{
                                               \frac{9}{2} +\pi^2+
                                               \frac{3}{2}\ln\left(1-y\right)-\ln^2\left(1-y\right)\right\}\,.
   \end{array}

In is important to notice the :math:`\theta`-function in the first term
of the r.h.s.. Despite it presence exposes the singularity of
:math:`1/(1-y)`, this is exactly what is needed to cancel the
singularities generated by :math:`\delta(1-y)\ln(1-y)` and
:math:`\delta(1-y)\ln^2(1-y)`. To see this, we write the logarithms in
the second line as:

.. math::
   \delta(1-y)\ln(1-y) = -\delta(1-y)\int_0^{1}\frac{dz}{1-z}\quad\mbox{and}\quad \delta(1-y)\ln^2(1-y) =
   -2 \delta(1-y)\int_0^{1}dz \frac{\ln(1-z)}{1-z}\,,

that combined with the :math:`1/(1-y)` and :math:`\ln(1-y)/(1-y)` terms
in the first line gives rise to the :math:`+`-prescription so that the
final result is:

.. math::
   :label: eq:coefff1
   \begin{array}{rcl}
     \displaystyle \frac{1}{2\pi}\lim_{\kappa\rightarrow 0}{\rm
     Im}\left[C_{q}^{(1)}(y,\kappa)\right] &=& \displaystyle C_F\theta(1-y)\bigg[
                                               2\left(\frac{\ln\left(1-y\right)}{1-y}\right)_+
                                               -\frac{3}{2}\left(\frac{1}{1-y}\right)_+
                                               -(1+y)
                                               \ln\left(1-y\right)\\
     \\
                                           &-&\displaystyle \frac{1+y^2}{1-y}\ln\left(y\right)+3- \delta(1-y)\left(\pi^2+\frac{9}{2}\right)\bigg]\,.
   \end{array}

This result *almost* perfectly agrees with the expected result. The
expression for the :math:`\mathcal{O}(\alpha_s)` correction to the quark
channel of :math:`F_1` can be read off from Eq. (4.85) of Ref. (Ellis,
Stirling, and Webber 2011). First we notice that the expansion parameter
in that reference is :math:`\alpha_s/(2\pi)` while we use
:math:`\alpha_s/(4\pi)`. Therefore, Eq. (4.85) has to be multiplied by a
factor 2 in the first place to match our notation. In addition,
Eq. (4.85) corresponds to :math:`F_2`. But since :math:`2xF_1=F_2-F_L`
and the corresponding coefficient function for :math:`F_L` is
:math:`C_{L,q}^{(1)}(x)=4C_Fx` (see Eq. (4.86) of Ref. (Ellis, Stirling,
and Webber 2011)), it is sufficient to subtract this term and finally
divide by 2 to obtain the coefficient function for :math:`F_1`. The net
result should be exactly Eq. (4.85) of Ref. (Ellis, Stirling, and Webber
2011) with the term :math:`2z` removed. Our
Eq. :eq:`eq:coefff1` matches up to the :math:`\pi^2`
term in the coefficient of :math:`\delta(1-y)` that multiplies a
different factor. This discrepancy needs to be investigated further.

We now consider the gluon coefficient function. In this case, due to the
overall factor that is singular in both :math:`y=1` and :math:`y=-1`, we
also need to consider the
:math:`(y\rightarrow -y, \kappa\rightarrow -\kappa)` term. This limit
gives:

.. math::
   \begin{array}{rcl}
   \lim_{\kappa\rightarrow 0}{\rm
    Im}\left[C_{g}^{(1)}(y,\kappa)\right] &=& \displaystyle 2\pi T_R \left(\frac{1}{1-y}+\frac{1}{1+y}\right) \left\{ 1+(1-2y+2y^2)\ln\left(\frac{1-y}{y}\right)\right\}\\
   \\
   &+& (y\rightarrow -y)\\
   \\
   &+& \displaystyle \pi T_R \left(\delta(1-y)-\delta(1+y)\right)\\
   \\
   &\times&\displaystyle \bigg\{ 
                            2\ln\left(\frac{1-y}{y}\right)
                            +(1-2y+2y^2)\left[\ln^2\left(\frac{1-y}{y}\right)-\pi^2\right]\bigg\}\\
   \\
   &+& (y\rightarrow -y) \,.
   \end{array}

**References**

**References**

.. container:: references csl-bib-body hanging-indent
   :name: refs

   .. container:: csl-entry
      :name: ref-Blumlein:1999sc

      Blumlein, Johannes, Bodo Geyer, and Dieter Robaschik. 1999. “The
      Virtual Compton amplitude in the generalized Bjorken region:
      twist-2 contributions.” *Nucl. Phys. B* 560: 283–344.
      https://doi.org/10.1016/S0550-3213(99)00418-6.

   .. container:: csl-entry
      :name: ref-Braun:2020yib

      Braun, V. M., A. N. Manashov, S. Moch, and J. Schoenleber. 2020.
      “Two-loop coefficient function for DVCS: vector contributions.”
      *JHEP* 09: 117. https://doi.org/10.1007/JHEP09(2020)117.

   .. container:: csl-entry
      :name: ref-Diehl:2003ny

      Diehl, M. 2003. “Generalized parton distributions.” *Phys. Rept.*
      388: 41–277. https://doi.org/10.1016/j.physrep.2003.08.002.

   .. container:: csl-entry
      :name: ref-Ellis:1991qj

      Ellis, R.Keith, W.James Stirling, and B. R. Webber. 2011. *QCD and
      collider physics*. Vol. 8. Cambridge University Press.

   .. container:: csl-entry
      :name: ref-Ji:1996nm

      Ji, Xiang-Dong. 1997. “Deeply virtual Compton scattering.” *Phys.
      Rev. D* 55: 7114–25. https://doi.org/10.1103/PhysRevD.55.7114.

   .. container:: csl-entry
      :name: ref-Ji:1998xh

      Ji, Xiang-Dong, and Jonathan Osborne. 1998. “One loop corrections
      and all order factorization in deeply virtual Compton scattering.”
      *Phys. Rev. D* 58: 094018.
      https://doi.org/10.1103/PhysRevD.58.094018.

   .. container:: csl-entry
      :name: ref-Lepage:1980fj

      Lepage, G. Peter, and Stanley J. Brodsky. 1980. “Exclusive
      Processes in Perturbative Quantum Chromodynamics.” *Phys. Rev. D*
      22: 2157. https://doi.org/10.1103/PhysRevD.22.2157.

   .. container:: csl-entry
      :name: ref-Mankiewicz:1997bk

      Mankiewicz, L., G. Piller, E. Stein, M. Vanttinen, and T. Weigl.
      1998. “NLO corrections to deeply virtual Compton scattering.”
      *Phys. Lett. B* 425: 186–92.
      https://doi.org/10.1016/S0370-2693(98)00190-7.

   .. container:: csl-entry
      :name: ref-Mikhailov:1984ii

      Mikhailov, S. V., and A. V. Radyushkin. 1985. “Evolution Kernels
      in {QCD}: Two Loop Calculation in Feynman Gauge.” *Nucl. Phys. B*
      254: 89–126. https://doi.org/10.1016/0550-3213(85)90213-5.

   .. container:: csl-entry
      :name: ref-Pire:2011st

      Pire, B., L. Szymanowski, and J. Wagner. 2011. “NLO corrections to
      timelike, spacelike and double deeply virtual Compton scattering.”
      *Phys. Rev. D* 83: 034009.
      https://doi.org/10.1103/PhysRevD.83.034009.

   .. container:: csl-entry
      :name: ref-Vinnikov:2006xw

      Vinnikov, A. V. 2006. “Code for prompt numerical computation of
      the leading order GPD evolution,” April.
      https://arxiv.org/abs/hep-ph/0604248.

.. [1]
   It should be noticed that the integration bounds of the integration
   in Eq. :eq:`eq:eveq` are dictated by the operator
   defintion of the distribution :math:`f` on the light cone and not by
   the kernel :math:`\mathbb{V}`.

.. [2]
   Notice the seemingly unusual fact that :math:`f^{+}` is defined as
   difference and :math:`f^{-}` as sum of GPDs computed at opposite
   values of :math:`x`. This can be understood from the fact that, in
   the forward limit, :math:`f(-x)= -\overline{f}(x)`, *i.e.* the PDF of
   a quark computed at :math:`-x` equals the PDF of the corresponding
   antiquark computed at :math:`x` with opposite sign. The opposite sign
   is absent in the longitudinally polarised case.

.. [3]
   There is probably a typo in Eq. (102) of Ref. (Diehl 2003) as the
   second :math:`-1` should actually be a :math:`+1`.

.. [4]
   The factor :math:`\theta\left(\frac{x}{x'}\right)` comes from the
   factor :math:`\theta(-x+x')` in
   Eq. :eq:`eq:supportdiehl` that can be rewritten
   as
   :math:`\theta\left(\frac{x}{x'}\right)\theta\left(1-\frac{x}{x'}\right)`.

.. [5]
   Note that all divergent integrals considered here are implicitly
   assumed to be principal-valued integrals such that:

   .. math:: \int_{-1}^{1}\frac{dt}{t}=0\,.

   \ \ This allows us to omit the :math:`\pm i\epsilon` terms.

.. [6]
   While the singularity at :math:`y=-1` is placed below the lower bound
   :math:`y=x` and thus does not cause any problem, there is an
   additional singularity at :math:`y=\pm1/\kappa` that needs to be
   considered. As discussed below, this singularity in
   :math:`\mathcal{P}_2^{-,(0)}` exactly cancels against an opposite
   singularity in :math:`\mathcal{P}_1^{-,(0)}`.

.. [7]
   Note that the two :math:`+`-prescribed terms when integrated over
   :math:`t` do not individually give zero but their combination does.

.. [8]
   An additional factor 2 is included by convention in the definition of
   :math:`F_1`.

.. [9]
   We are assuming that :math:`Q\ll M_Z` so that weak corrections to the
   electroweak coupling can be neglected

.. [10]
   Despite not contributing to the imaginary part, notice that:

   .. math:: \lim_{\kappa\rightarrow 0}\frac{\ln(1-k)}{\kappa}=-1\quad\mbox{and}\quad  \lim_{\kappa\rightarrow 0}\frac{\ln^2(1-k)}{\kappa}=0\,.
