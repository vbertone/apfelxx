.. contents::
   :depth: 3
..

| **:math:`\overline{\mbox{MS}}` Mass Implementation**
| VB

Introduction
============

Until recently, heavy quark masses in PDF fits were mostly treated in
the pole-mass (or on-shell) renormalization scheme (REFERENCE). However,
in a recent paper (Alekhin and Moch 2011) it has been shown that the
implementation of the :math:`\overline{\mbox{MS}}` masses results in an
improvement of the perturbative stability and in the consequent
reduction of the theoretical uncertainty due to variations of
renormalization and factorization scales. As a consequence, this allows
one to a more reliabe determination of the numerical value of the heavy
quark masses themselves.

In this set of notes we will describe the implemntation of the
:math:`\overline{\mbox{MS}}` heavy-quark running masses in ``APFEL``.
Strarting from the more used defition of structure functions and in
terms of pole masses, our goal is to consistently replace them with the
:math:`\overline{\mbox{MS}}` mass definition.

:math:`\overline{\mbox{MS}}` Running Mass vs. Pole Mass
=======================================================

The (scale independent) pole mass :math:`M` and the (scale dependent)
:math:`\overline{\mbox{MS}}` mass :math:`m(\mu)` arise from two
different renormalization procedures and in perturbation theory they can
be expressed on in terms of the other. The relation connecting pole and
:math:`\overline{\mbox{MS}}` mass definitions has been computed in
Ref. (Chetyrkin and Steinhauser 2000) up to four loops. , and in
particular from eq. (8) of that paper we read the ratio between
:math:`M` and :math:`m(\mu)`:

.. math::
   :label: straight
   \frac{m(\mu)}{M} = 1 + z^{(1)}a_s+\underbrace{C_F\left[C_Fz_1^{(2)}+C_Az_2^{(2)}+T_R\left(N_Lz_3^{(2)}+z_4^{(2)}\right)\right]}_{z^{(2)}}a_s^2+\mathcal{O}(a_s^3)\,,

where:

.. math::
   :label: coeff
   \begin{array}{l}
   \displaystyle z^{(1)}(\mu,M) = C_F\left(- 4 - 3L_{\mu M}\right)\\
   \\
   \displaystyle z_{1}^{(2)}(\mu,M) = \frac{7}{8} - 30\zeta_2 - 12\zeta_3 + 48\zeta_2\ln(2) + \frac{21}2L_{\mu M} + \frac92 L_{\mu M}^2\\
   \\
   \displaystyle z_{2}^{(2)}(\mu,M) = -\frac{1111}{24} + 8\zeta_2 +6\zeta_3 - 24\zeta_2\ln(2) - \frac{185}6L_{\mu M} - \frac{11}2 L_{\mu M}^2\\
   \\
   \displaystyle z_{3}^{(2)}(\mu,M) = \frac{71}{6} + 8\zeta_2 + \frac{26}3L_{\mu M} + 2 L_{\mu M}^2\\
   \\
   \displaystyle z_{4}^{(2)}(\mu,M) = \frac{143}{6} -16\zeta_2 + \frac{26}3L_{\mu M} + 2 L_{\mu M}^2\,,\\
   \end{array}

where we have defined:

.. math:: L_{\mu M} = \ln\frac{\mu^2}{M^2}

and where :math:`N_L` is the number of light (massless) quarks (i.e.
:math:`N_L=3` for the charm and :math:`N_L=4` for the bottom). Moreover,
in our notation we use( [1]_):

.. math:: a_s\equiv a_s(\mu) = \frac{\alpha_s(\mu)}{4\pi}\,.

But what we really need is the inverse of eq.
:eq:`straight` in terms of :math:`\mu` and :math:`m(\mu)`.
We could derive it inverting eq. :eq:`straight` and then
writing :math:`M` in terms of :math:`m(\mu)`, but more easily we can
read it from eq. (31) of (Chetyrkin and Steinhauser 2000) and write it
as:

.. math::
   :label: inverse
   \frac{M}{m(\mu)} = 1 + h^{(1)}a_s+\underbrace{C_F\left[C_Fh_1^{(2)}+C_Ah_2^{(2)}+T_R\left(N_Lh_3^{(2)}+h_4^{(2)}\right)\right]}_{h^{(2)}}a_s^2+\mathcal{O}(a_s^3)\,,

with:

.. math::
   :label: coeffinv
   \begin{array}{l}
   \displaystyle h^{(1)}(\mu,m(\mu)) = C_F\left(4 + 3L_{\mu m}\right)\\
   \\
   \displaystyle h_{1}^{(2)}(\mu,m(\mu)) = -\frac{7}{8} + 30\zeta_2 + 12\zeta_3 - 48\zeta_2\ln(2) - \frac{9}2L_{\mu m} + \frac92 L_{\mu m}^2\\
   \\
   \displaystyle h_{2}^{(2)}(\mu,m(\mu)) = \frac{1111}{24} - 8\zeta_2 - 6\zeta_3 + 24\zeta_2\ln(2) + \frac{185}6L_{\mu m} + \frac{11}2 L_{\mu m}^2\\
   \\
   \displaystyle h_{3}^{(2)}(\mu,m(\mu)) = - \frac{71}{6} - 8\zeta_2 - \frac{26}3L_{\mu m} - 2 L_{\mu m}^2\\
   \\
   \displaystyle h_{4}^{(2)}(\mu,m(\mu)) = -\frac{143}{6} + 16\zeta_2 - \frac{26}3L_{\mu m} - 2 L_{\mu m}^2\,,\\
   \end{array}

where now we have defined:

.. math:: L_{\mu m} = \ln\frac{\mu^2}{m^2(\mu)}\,.

In the following we will use eq. :eq:`inverse` to replace
the pole mass :math:`M` with the :math:`\overline{\mbox{MS}}` mass
:math:`m(\mu)`.

RGE Solution for the :math:`\overline{\mbox{MS}}` Running Mass
==============================================================

Actally, what we also need is to know how :math:`m(\mu)` runs with the
renormalization scale :math:`\mu`. To this end, we use the RGE to obtain
the following differential equation for the running:

.. math::
   :label: andim
   \mu^2\frac{dm}{d\mu^2} = m(\mu)\gamma_m(a_s) = -m(\mu)\sum_{n=0}^{\infty}\gamma_m^{(n)}a_s^{n+1}

and from Eqs. (46), (47) and (48) of (Chetyrkin and Retey 2000) for
:math:`SU(3)` and taking into account a factor 4 difference in the
definition of :math:`a_s`, we read:

.. math:: \gamma_m^{(0)} = 4

.. math:: \gamma_m^{(1)} = \frac{202}3 - \frac{20}{9}N_f

.. math:: \gamma_m^{(2)} = 1249 - \left(\frac{2216}{27}+\frac{160}{3}\zeta_3\right)N_f-\frac{140}{81}N_f^2\,,

where :math:`N_f` is the number of active flavours. But we also know how
:math:`a_s` runs, that is:

.. math::
   :label: betaf
   \mu^2\frac{da_s}{d\mu^2} = \beta(a_s) = -\sum_{n=0}^{\infty}\beta_n a_s^{n+2}\,,

with:

.. math:: \beta_0 = 11-\frac23 N_f

\ 

.. math:: \beta_1 = 102 - \frac{38}3 N_f

.. math:: \beta_2 = \frac{2857}{2} - \frac{5033}{18}N_f + \frac{325}{54}N_f^2

and from eq. :eq:`betaf` it follows that:

.. math:: \mu^2\frac{dm}{d\mu^2} = \beta(a_s)\frac{d m}{da_s}\,,

\ so that the differential equation in eq. :eq:`andim` can be
written as:

.. math::
   :label: runmass
   \frac{dm}{da_s} = \frac{\gamma_m(a_s)}{\beta(a_s)}m(a_s)\,.

The formal solution of eq. :eq:`runmass` reads:

.. math::
   :label: numsol
   m(\mu) = m(\mu_0)\exp\left[\int_{a_s(\mu_0)}^{a_s(\mu)}\frac{\gamma_m(a_s)}{\beta(a_s)}da_s\right]\,,

then we expand the integrand in eq. :eq:`numsol` using the
perturbative expansion of the :math:`\gamma_m(a_s)` and
:math:`\beta(a_s)` functions given in eqs. :eq:`andim` and
:eq:`betaf` obtaining the following polynomial:

.. math:: \frac{\gamma_m(a)}{\beta(a)} = \frac1{a}\left[c_0 + (c_1-b_1c_0)a+(c_2-c_1b_1-b_2c_0+b_1^2c_0)a^2+\mathcal{O}(a^3)\right]

where we have defined:

.. math::
   :label: jhgkgfkgf
   \left\{\begin{array}{l}
   \displaystyle b_i = \frac{\beta_i}{\beta_0}\\
   \\
   \displaystyle c_i = \frac{\gamma_m^{(i)}}{\beta_0}
   \end{array}\right.\,.

We integrate eq. :eq:`numsol` getting:

.. math::
   :label: integral
   \int_{a_0}^a\frac{\gamma_m(a)}{\beta(a)}da = c_0\ln\frac{a}{a_0} + (c_1-b_1c_0)(a-a_0)+\frac12(c_2-c_1b_1-b_2c_0+b_1^2c_0)(a^2-a_0^2)\,,

where :math:`a\equiv a_s(\mu)` and :math:`a_0\equiv a_s(\mu_0)`. After
that, we put it in the exponential function and expand again, finally
obtaining:

.. math::
   :label: ansol
   m(\mu)=m(\mu_0)\left(\frac{a}{a_0}\right)^{c_0}\frac{1+(c_1-b_1c_0)a+\frac12[c_2-c_1b_1-b_2c_0+b_1^2c_0+(c_1-b_1c_0)^2]a^2}{1+(c_1-b_1c_0)a_0+\frac12[c_2-c_1b_1-b_2c_0+b_1^2c_0+(c_1-b_1c_0)^2]a_0^2}\,,

which gives the NNLO running of :math:`m(\mu)`. Of course, to obatin the
NLO running one has just to disregard the terms proportional to
:math:`a^2` and :math:`a_0^2` in the ration, while at LO also the terms
proportional to :math:`a` and :math:`a_0` should be omitted( [2]_).

Matching Conditions
===================

As one can see from eq. :eq:`ansol`, the running of the
:math:`\overline{\mbox{MS}}` mass :math:`m(\mu)` requires the value of
:math:`\alpha_s` at the scales :math:`\mu` and :math:`\mu_0` (i.e.
:math:`a` and :math:`a_0`). But in turn the running of :math:`\alpha_s`
itself depends on the values of the heavy quark mass thresholds by means
of the so called matching conditions, which tell essentially us how to
perform the switching of the running from :math:`N_f` to :math:`N_f+1`
active flavours. So, there seems to be a circular problem. But in the
following we will see how to get out of this using the fact that the
scale where to perform the matching is arbitrary. We will use this
arbitrariness to see how to compute the running of :math:`\alpha_s`
without knowing the running of :math:`m(\mu)`.

In general the matching conditions give rise to a discontinuity of
:math:`\alpha_s` at the matching scale and in the present code they are
written in terms of the pole masses. These masses are scale independent
and are given as input parameters, therefore they don’t give any
problem. Moreover, one can show that if the matching point :math:`\mu`
is chosen to be equal to the pole mass :math:`M`, the discontinuity
appears only at NNLO.

Now, the first step to replace the pole mass :math:`M` with the
:math:`\overline{\mbox{MS}}` mass :math:`m(\mu)` is to rewrite the
matching conditions for :math:`\alpha_s` in terms of the
:math:`\overline{\mbox{MS}}` mass rather than the pole mass. This is
exactly what we are going to do in the following. Then we will find that
choosing this time :math:`\mu=m(\mu)` the discontinuity appears again
only at NNLO but with a different coefficient.

As known, the same problem holds for PDFs. In fact, also PDFs need to be
matched and in the following we will discuss also how to write the
matching conditions for PDFs, which originally are given for the pole
mass :math:`M`, in terms of the :math:`\overline{\mbox{MS}}` mass
:math:`m(\mu)`.

.. container::
   :name: matching-of-alpha_smu

   .. rubric:: Matching of :math:`\alpha_s(\mu)`
      :name: matching-of-alpha_smu

In this section we will show how to express in terms of the
:math:`\overline{\mbox{MS}}` mass :math:`m(\mu)` the matching condition
for :math:`\alpha_s`. We took the matching condition for
:math:`\alpha_s` from eq. (2.41) of (Vogt 2005), which in turn was taken
from eq. (9) of (Chetyrkin, Kniehl, and Steinhauser 1997). Here we write
this equation (up to NNLO and taking into account a factor 4 coming from
the different definition of :math:`a`) as follows:

.. math::
   :label: alphaspole
   \frac{a^{(n-1)}(\mu)}{a^{(n)}(\mu)}=1-\frac23 L_{\mu M}a^{(n)}(\mu)+\left(\frac49L_{\mu M}^2-\frac{38}3L_{\mu M}-\frac{14}3\right)[a^{(n)}(\mu)]^2\,.

being :math:`M` the pole mass of the :math:`n`-th flavour. But from eq.
:eq:`inverse` we read:

.. math:: \ln M^2 = \ln m^2(\mu) + 2\ln[1+h^{(1)}(\mu)a^{(n)}(\mu)]\quad\mbox{with}\quad h^{(1)}(\mu) = \frac{16}3+4L_{\mu m}

that, using the expansion:

.. math:: \ln(1+x)=\sum_{k=1}^{\infty}\frac{(-1)^{k+1}}{k}x^k\,,

can be written as:

.. math::
   :label: cacchiocacchio
   \ln M^2 = \ln m^2(\mu) + 2h^{(1)}(\mu)a^{(n)}(\mu)+\mathcal{O}([a^{(n)}]^2)\,.

Therefore it is straightforward to see that:

.. math:: L_{\mu M}  = L_{\mu m} - 2h^{(1)}a^{(n)}=L_{\mu m}-\left(\frac{32}3+8L_{\mu m}\right)a^{(n)}\,,

so that:

.. math::
   :label: alphasmsbar
   \frac{a^{(n-1)}(\mu)}{a^{(n)}(\mu)}=1-\frac23 L_{\mu m}a^{(n)}(\mu)+\left(\frac49L_{\mu m}^2-\frac{22}3L_{\mu m}+\frac{22}9\right)[a^{(n)}(\mu)]^2\,.

Now, in order to get rid of the logarithmic terms, we choose to match
:math:`a^{(n-1)}` and :math:`a^{(n)}` at :math:`\mu=m(\mu)=m(m)` so that
we get:

.. math:: a^{(n-1)}(m)=a^{(n)}(m)\left(1+\frac{22}9[a^{(n)}(m)]^2\right)\,,

which can be easily inverted obtaining:

.. math::
   :label: pollopollo
   a^{(n)}(m)=a^{(n-1)}(m)\left(1-\frac{22}9[a^{(n-1)}(m)]^2\right)\,.

So, exactly as it happened in the case of the pole mass, also in the
:math:`\overline{\mbox{MS}}` mass case we can make the matching
condition for :math:`\alpha_s` start to play a role only at NNLO. But
the difference is that now the coefficient of the matching is
:math:`-22/9` rather than :math:`14/3`. It is interesting to observe
that, in order to perform the matching as described above, we just need
to know the value of :math:`m(m)`. This is the so called RG-invariant
:math:`\overline{\mbox{MS}}` mass and this will be given as input
parameter in place of the pole mass :math:`M`. It is not by chance that
the PDG provides exactly the values for :math:`m_c(m_c)` and
:math:`m_b(m_b)`. This finally allows to evaluate the values of
:math:`a` and :math:`a_0` that enter in eq. :eq:`ansol`
needed to compute the running of :math:`m(\mu)`.

.. container::
   :name: matching-of-mmu

   .. rubric:: Matching of :math:`m(\mu)`
      :name: matching-of-mmu

On the other hand, also the running of :math:`m(\mu)` needs to be
matched. In particular, we might need to match the :math:`(n-1)` with
:math:`(n)` scheme of the mass :math:`m_q(\mu)`, with :math:`q=c,b,t`,
being :math:`m_h(\mu)`, with again :math:`h=c,b,t` the :math:`n`-th
threshold. From eqs. (26) and (27) of (Chetyrkin, Kuhn, and Steinhauser
2000) one reads that:

.. math::
   :label: mqmc
   \frac{m_q^{(n-1)}(\mu)}{m_q^{(n)}(\mu)}=1+\left(\frac43L_{\mu m}^{(h)2}-\frac{20}9L_{\mu m}^{(h)}+\frac{89}{27}\right)[a^{(n)}(\mu)]^2\,,

where:

.. math:: L_{\mu m}^{(h)} =\ln\frac{\mu^2}{m_h^2(\mu)}\,.

Exactly as before, if we choose to match the two schemes at the scale
:math:`\mu=m_h(\mu)=m_h(m_h)`, the logarithmic terms vanish and we are
left with:

.. math:: m_q^{(n-1)}(m_h)=\left(1+\frac{89}{27}[a^{(n)}(m_h)]^2\right)m_q^{(n)}(m_h)=\zeta_n^{(\mbox{\tiny down})}(m_h)m_q^{(n)}(m_h)

whose inverse is:

.. math:: m_q^{(n)}(m_h)=\left(1-\frac{89}{27}[a^{(n-1)}(m_h)]^2\right)m_q^{(n-1)}(m_h)=\zeta_n^{(\mbox{\tiny up})}(m_h)m_q^{(n-1)}(m_h)

.. container::
   :name: implementation-of-the-running-in-the-vfn-scheme

   .. rubric:: Implementation of the Running in the VFN Scheme
      :name: implementation-of-the-running-in-the-vfn-scheme

In this section we will see how the running for :math:`m_q(\mu)` in VFN
scheme can be implemented in terms of the :math:`\overline{\mbox{MS}}`
masses. Let us start with an example: we want to show how to compute
:math:`m_c^{(5)}(\mu)`, assuming that :math:`m_b<\mu<M_Z`. Our input
parameters are :math:`a_s^{(5)}(M_Z)`, :math:`m_c^{(4)}(m_c)` and
:math:`m_b^{(5)}(m_b)`. First of all, starting from
:math:`a_s^{(5)}(M_Z)`, using the equation for the running and the
matching condition in eq. :eq:`pollopollo`, we evaluate
in order: :math:`a_s^{(5)}(\mu)`, :math:`a_s^{(5)}(m_b)`,
:math:`a_s^{(4)}(m_b)` and :math:`a_s^{(4)}(m_c)`. Once we have this
values, we use eq. :eq:`numsol` to write the evolution as
follows:

.. math:: m_c^{(5)}(\mu) = m_c^{(4)}(m_c)\underbrace{\exp\left[\int_{a_s^{(4)}(m_c)}^{a_s^{(4)}(m_b)}\frac{\gamma_m^{(4)}(a)}{\beta^{(4)}(a)}da\right]}_{E_m^{(4)}(m_c,m_b)}\zeta_5^{(\mbox{\tiny up})}(m_b)\underbrace{\exp\left[\int_{a_s^{(5)}(m_b)}^{a_s^{(5)}(\mu)}\frac{\gamma_m^{(5)}(a)}{\beta^{(5)}(a)}da\right]}_{E_m^{(5)}(m_b,\mu)}\,,

where, using eq. :eq:`ansol` up to NNLO, we have that:

.. math::
   \begin{array}{l}
   \displaystyle E_m^{(4)}(m_c,m_b)=\left(\frac{a_s^{(4)}(m_b)}{a_s^{(4)}(m_c)}\right)^{c_0}\frac{1+(c_1-b_1c_0)a_s^{(4)}(m_b)+\frac12[c_2-c_1b_1-b_2c_0+b_1^2c_0+(c_1-b_1c_0)^2][a_s^{(4)}(m_b)]^2}{1+(c_1-b_1c_0)a_s^{(4)}(m_c)+\frac12[c_2-c_1b_1-b_2c_0+b_1^2c_0+(c_1-b_1c_0)^2][a_s^{(4)}(m_c)]^2}\\
   \\
   \displaystyle E_m^{(5)}(m_b,\mu)=\left(\frac{a_s^{(5)}(\mu)}{a_s^{(5)}(m_b)}\right)^{c_0}\frac{1+(c_1-b_1c_0)a_s^{(5)}(\mu)+\frac12[c_2-c_1b_1-b_2c_0+b_1^2c_0+(c_1-b_1c_0)^2][a_s^{(5)}(\mu)]^2}{1+(c_1-b_1c_0)a_s^{(5)}(m_b)+\frac12[c_2-c_1b_1-b_2c_0+b_1^2c_0+(c_1-b_1c_0)^2][a_s^{(5)}(m_b)]^2}\,,
   \end{array}

where we are also assuming that the coefficients :math:`b_i` and
:math:`c_i`, given in eq. :eq:`jhgkgfkgf`, in :math:`E_1`
are compute with :math:`N_f=4` while in :math:`E_2` with :math:`N_f=5`.

Looking a the above example and noticing that:

.. math:: E_m^{(i)}(\mu_1,\mu_2) = \frac1{E_m^{(i)}(\mu_2,\mu_1) }\,,

\ one can write down the general pattern for the mass evolution of the
charm:

.. math::
   m_c^{(n)}(\mu)=\left\{
   \begin{array}{ll}
   \displaystyle \frac{1}{E_m^{(3)}(\mu,m_c)}\zeta_4^{(\mbox{\tiny down})}(m_c)m_c^{(4)}(m_c)&\quad \mu<m_c\\
   \\
   \displaystyle m_c^{(4)}(m_c)E_m^{(4)}(m_c,\mu)&\quad m_c\leq\mu<m_b\\
   \\
   \displaystyle m_c^{(4)}(m_c)E_m^{(4)}(m_c,m_b)\zeta_5^{(\mbox{\tiny up})}(m_b)E_m^{(5)}(m_b,\mu)&\quad m_b\leq\mu<m_t\\
   \\
   \displaystyle m_c^{(4)}(m_c)E_m^{(4)}(m_c,m_b)\zeta_5^{(\mbox{\tiny up})}(m_b)E_m^{(5)}(m_b,m_t)\zeta_6^{(\mbox{\tiny up})}(m_t)E^{(6)}(m_t,\mu)&\quad\mu \geq m_t
   \end{array}\right.

and guess also the pattern for bottom:

.. math::
   m_b^{(n)}(\mu)=\left\{
   \begin{array}{ll}
   \displaystyle \frac{1}{E_m^{(3)}(\mu,m_c)}\zeta_4^{(\mbox{\tiny down})}(m_c)\frac{1}{E_m^{(4)}(m_c,m_b)}\zeta_5^{(\mbox{\tiny down})}(m_b)m_b^{(5)}(m_b)&\quad \mu<m_c\\
   \\
   \displaystyle \frac{1}{E_m^{(4)}(\mu,m_b)}\zeta_5^{(\mbox{\tiny down})}(m_b)m_b^{(5)}(m_b)&\quad m_c\leq\mu<m_b\\
   \\
   \displaystyle m_b^{(5)}(m_b)E_m^{(5)}(m_b,\mu)&\quad m_b\leq\mu<m_t\\
   \\
   \displaystyle m_b^{(5)}(m_b)E_m^{(5)}(m_b,m_t)\zeta_6^{(\mbox{\tiny up})}(m_t)E^{(6)}(m_t,\mu)&\quad\mu \geq m_t
   \end{array}\right.

and top:

.. math::
   m_t^{(n)}(\mu)=\left\{
   \begin{array}{ll}
   \displaystyle \frac{1}{E_m^{(3)}(\mu,m_c)}\zeta_4^{(\mbox{\tiny down})}(m_c)\frac{1}{E_m^{(4)}(m_c,m_b)}\zeta_5^{(\mbox{\tiny down})}(m_b)\frac{1}{E_m^{(5)}(m_b,m_t)}\zeta_6^{(\mbox{\tiny down})}(m_t)m_t^{(6)}(m_t)&\quad \mu<m_c\\
   \\
   \displaystyle \frac{1}{E_m^{(4)}(\mu,m_b)}\zeta_5^{(\mbox{\tiny down})}(m_b)\frac{1}{E_m^{(5)}(m_b,m_t)}\zeta_6^{(\mbox{\tiny down})}(m_t)m_t^{(6)}(m_t)&\quad m_c\leq\mu<m_b\\
   \\
   \displaystyle \frac{1}{E_m^{(5)}(\mu,m_t)}\zeta_6^{(\mbox{\tiny down})}(m_t)m_t^{(6)}(m_t)&\quad m_b\leq\mu<m_t\\
   \\
   \displaystyle m_t^{(6)}(m_t)E^{(6)}(m_t,\mu)&\quad\mu \geq m_t
   \end{array}\right.\,.

As one can see from the above equations, having the evolution equation
for a given number of active flavours, having the matching conditions
and having the input parameters :math:`m_i^{(n)}(m_i)` we can evaluate
the value of :math:`m_i(\mu)` in the VFN scheme for any value of
:math:`\mu`. Notice that we are assuming that the input parameters
:math:`m_i(m_i)` are always given just beyond the threshold, i.e. we
assume them to be equal to :math:`m_c^{(4)}(m_c)`,
:math:`m_b^{(5)}(m_b)` and :math:`m_t^{(6)}(m_t)`.

In the FFN scheme with :math:`N_f` flavours, instead, the situation is
much easier. In fact, the evolution reduces to( [3]_):

.. math:: m_i^{(N_f)}(\mu) = m_i(m_i)E_m^{(N_f)}(m_i,\mu)\,.

.. container::
   :name: pdfs

   .. rubric:: PDFs
      :name: pdfs

Now we see how to replace the pole mass with the
:math:`\overline{\mbox{MS}}` mass in the matching conditions for the
PDFs. One can write the singlet and the gluon in the :math:`N_f+1`
scheme in terms of singlet and gluon in :math:`N_f` scheme at any scale
:math:`\mu` as follows( [4]_ is undersato to be in the :math:`N_f`
scheme.):

.. math::
   :label: couple1
   \begin{array}{c}
   \displaystyle {\Sigma^{(N_f+1)}(\mu) \choose g^{(N_f+1)}(\mu)}=\begin{pmatrix}1+a_s^2[A_{qq,h}^{N\!S,(2)}+\tilde{A}^{S,(2)}_{hq}] & a_s\tilde{A}^{S,(1)}_{hg}+a_s^2\tilde{A}^{S,(2)}_{hg}\\
   a_s^2A^{S,(2)}_{gq,h} & 1+a_sA_{gg,h}^{S,(1)}+a_s^2A_{gg,h}^{S,(2)}\end{pmatrix}\left(x,\frac{\mu^2}{M^2}\right){\Sigma^{(N_f)}(\mu) \choose g^{(N_f)}(\mu)}\,,
   \end{array}

\ where the form of the (:math:`x`-space) functions entering the
transformation matrix above are given in Appendix B of (Buza et al.
1998) and where the "1"s are indeed delta functions. We notice here
that, at the generic scale :math:`\mu`, also the
:math:`\mathcal{O}(a_s)` functions :math:`\tilde{A}^{S,(1)}_{hg}` and
:math:`A_{gg,h}^{S,(1)}` are there. But given that they are both
proportional to :math:`\ln(\mu^2/M^2)`, they disappear if one chooses to
perform the matching at the threshold :math:`\mu=M`. We omit the
matching conditions for the non-singlet PDFs because such
:math:`\mathcal{O}(a_s)` functions are not present there and they are
the point of the discussion.

Now, if we want to replace the pole mass :math:`M` with the
:math:`\overline{\mbox{MS}}` mass :math:`m(\mu)`, we just have to plug
eq. :eq:`cacchiocacchio` into eq.
:eq:`couple1`. In the :math:`\mathcal{O}(a_s^2)` functions
the second term in the l.h.s. of eq.
:eq:`cacchiocacchio` would give rise to subleading
terms. Therefore in those functions it is enough to replace :math:`M`
with :math:`m(\mu)`. On the other hand, the second term in the l.h.s. of
eq. :eq:`cacchiocacchio` is important in the
:math:`\mathcal{O}(a_s)` functions. Since both the functions
:math:`\tilde{A}^{S,(1)}_{hg}` and :math:`A_{gg,h}^{S,(1)}` are
proportional to :math:`\ln(\mu^2/M^2)`, they can be written as:

.. math::
   :label: ahciccio!
   \begin{array}{l}
   \displaystyle \tilde{A}^{S,(1)}_{hg}\left(x,\frac{\mu^2}{M^2}\right) = f_1(x)\ln\frac{\mu^2}{M^2}\\
   \\
   \displaystyle A^{S,(1)}_{gg,h}\left(x,\frac{\mu^2}{M^2}\right) = f_2(x)\ln\frac{\mu^2}{M^2}
   \end{array}\,,

\ where:

.. math::
   :label: ahbelli!
   \begin{array}{l}
   \displaystyle f_1(x)= 4 T_R[x^2+(1-x)^2]\\
   \\
   \displaystyle f_2(x)= -\frac43 T_R \delta(1-x)
   \end{array}\,.

Replacing :math:`M` with :math:`m` in eqs.
(`[ahciccio!] <#ahciccioux21>`__) using eq.
:eq:`cacchiocacchio`, we get:

.. math::
   \begin{array}{l}
   \displaystyle \tilde{A}^{S,(1)}_{hg}\left(x,\frac{\mu^2}{m^2}\right) = f_1(x)\ln\frac{\mu^2}{m^2}-2h^{(1)}(\mu)f_1(x)a_s(\mu)\\
   \\
   \displaystyle A^{S,(1)}_{gg,h}\left(x,\frac{\mu^2}{m^2}\right) = f_2(x)\ln\frac{\mu^2}{m^2}-2h^{(1)}(\mu)f_2(x)a_s(\mu)
   \end{array}\,.

Therefore eq. :eq:`couple1` in terms of :math:`m` becomes:

.. math::
   \begin{array}{c}
   \displaystyle {\Sigma^{(N_f+1)}\choose g^{(N_f+1)}}=\begin{pmatrix}1+a_s^2[A_{qq,h}^{N\!S,(2)}+\tilde{A}^{S,(2)}_{hq}] & a_s\tilde{A}^{S,(1)}_{hg}+a_s^2[\tilde{A}^{S,(2)}_{hg}-2h^{(1)}f_1]\\
   a_s^2A^{S,(2)}_{gq,h} & 1+a_sA_{gg,h}^{S,(1)}+a_s^2[A_{gg,h}^{S,(2)}-2h^{(1)}f_2]\end{pmatrix}\left(x,\frac{\mu^2}{m^2}\right){\Sigma^{(N_f)} \choose g^{(N_f)}}
   \end{array}\,.

Given that our code works in the Mellin space, the functions in eq.
(`[ahbelli!] <#ahbelliux21>`__) need to be Mellin-transformed obtaining:

.. math::
   :label: ahbellissimi!
   \begin{array}{l}
   \displaystyle f_1(N) = \mathbf{M}[f_1(x)])(N)= 4 T_R\left[\frac2{N+2}-\frac2{N+1}+\frac1N\right]\\
   \\
   \displaystyle f_2(N) = \mathbf{M}[f_2(x)])(N)= -\frac43 T_R
   \end{array}

Now, we choose to match :math:`N_f+1` and :math:`N_f` schemes at
:math:`\mu = m(\mu) = m(m)` so that all the logarithmic terms vanish
(including the functions :math:`\tilde{A}^{S,(1)}_{hg}` and
:math:`A_{gg,h}^{S,(1)}`) obtaining:

.. math::
   \begin{array}{c}
   \displaystyle {\Sigma^{(N_f+1)}\choose g^{(N_f+1)}}=\begin{pmatrix}1+a_s^2[A_{qq,h}^{N\!S,(2)}+\tilde{A}^{S,(2)}_{hq}] & a_s^2[\tilde{A}^{S,(2)}_{hg}-2h^{(1)}f_1]\\
   a_s^2A^{S,(2)}_{gq,h} & 1+a_s^2[A_{gg,h}^{S,(2)}-2h^{(1)}f_2]\end{pmatrix}(x){\Sigma^{(N_f)} \choose g^{(N_f)}}\,.
   \end{array}

Again, as it happened in the the pole mass case, if one chooses cleverly
the matching point also in the :math:`\overline{\mbox{MS}}` mass case
the matching conditions for PDFs start to play a role only at NNLO, but
the difference now is that some coefficients of the matching matrix for
gluon and singlet result modified by some simple functions.

One might want to invert the matching conditions for PDFs and considerig
that eq :eq:`couple1` can be written as:

.. math::
   :label: couple1exp
   \begin{array}{c}
   \displaystyle {\Sigma^{(N_f+1)}(\mu) \choose g^{(N_f+1)}(\mu)}=\left[
   \begin{pmatrix}
   1 & 0 \\
   0 &1
   \end{pmatrix}
   +a_s
   \begin{pmatrix}
   0 & \tilde{A}^{S,(1)}_{hg}\\
   0 & A_{gg,h}^{S,(1)}
   \end{pmatrix}
   +a_s^2
   \begin{pmatrix}
   [A_{qq,h}^{N\!S,(2)}+\tilde{A}^{S,(2)}_{hq}] & \tilde{A}^{S,(2)}_{hg}\\
   A^{S,(2)}_{gq,h} & A_{gg,h}^{S,(2)}
   \end{pmatrix}
   \right]
   {\Sigma^{(N_f)}(\mu) \choose g^{(N_f)}(\mu)}\,,
   \end{array}

\ and using the expansion:

.. math:: \frac{1}{1+b_1x+b_2x^2} =  1-b_1x-(b_2-b_1^2)x^2 + \mathcal{O}(x^3)\,,

we can immediately write:

.. math::
   :label: couple1expInv
   \begin{array}{c}
   \displaystyle {\Sigma^{(N_f)}(\mu) \choose g^{(N_f)}(\mu)}=\left[
   \begin{pmatrix}
   1 & 0 \\
   0 &1
   \end{pmatrix}
   -a_s
   \begin{pmatrix}
   0 & \tilde{A}^{S,(1)}_{hg}\\
   0 & A_{gg,h}^{S,(1)}
   \end{pmatrix}
   -a_s^2
   \begin{pmatrix}
   [A_{qq,h}^{N\!S,(2)}+\tilde{A}^{S,(2)}_{hq}] & \tilde{A}^{S,(2)}_{hg}
   - \tilde{A}^{S,(1)}_{hg} \otimes A_{gg,h}^{S,(1)} \\
   A^{S,(2)}_{gq,h} & A_{gg,h}^{S,(2)} - A_{gg,h}^{S,(1)} \otimes A_{gg,h}^{S,(1)}
   \end{pmatrix}
   \right]
   {\Sigma^{(N_f+1)}(\mu) \choose g^{(N_f+1)}(\mu)}\,,
   \end{array}

and using eqs. (`[ahciccio!] <#ahciccioux21>`__)
and (`[ahbelli!] <#ahbelliux21>`__), we get:

.. math::
   \begin{array}{l}
   \displaystyle   - \tilde{A}^{S,(1)}_{hg} \otimes A_{gg,h}^{S,(1)} =
      \frac{16}{3}T_R^2\left[x^2+(1-x)^2\right]\ln^2\frac{\mu^2}{M^2} = \frac{4}{3}T_R \ln\frac{\mu^2}{M^2}\tilde{A}^{S,(1)}_{hg}\\
   \\
   \displaystyle - A_{gg,h}^{S,(1)} \otimes A_{gg,h}^{S,(1)} =
      - \frac{16}{9}T_R^2\delta(1-x) \ln^2\frac{\mu^2}{M^2} = \frac{4}{3}T_R \ln\frac{\mu^2}{M^2}A_{gg,h}^{S,(1)}
      \end{array}

.. container::
   :name: renormalization-scale-variation

   .. rubric:: Renormalization Scale Variation
      :name: renormalization-scale-variation

The scale :math:`\mu` that appears in the running of :math:`a_s` and
:math:`m_q` is the *renormalization* scale, which now we will indicate
as :math:`\mu_R`, while the scale that explicitly appears in the PDF
evolution is the *factorization* scale, which now we will call
:math:`\mu_F`, and in principle they could be taken different and in
general one can write :math:`\mu_R = \kappa \mu_F`, where :math:`\kappa`
can be any real number( [5]_).

What one usually does for the choice of the matching points is to set
:math:`\mu_F` to heavy quark thresholds (:math:`M_c`, :math:`M_b` and
:math:`M_t` in the Pole Mass scheme and :math:`m_c(m_c)`,
:math:`m_b(m_b)` and :math:`m_t(m_t)` in the
:math:`\overline{\mbox{MS}}` scheme). In this way the logarithmic terms
in the PDF matching conditions are assured to vanish guaranteeing the
same matching pattern for PDFs. But if :math:`\kappa` is different from
one, the logarithmic terms in the matching conditions for
:math:`a_s(\mu_R)` and :math:`m_q(\mu_R)` don’t vanish anymore. This is
exactly the case when one wants to perform the renormalization scale
variation and in the following we will just show how to implement the
matching condition for :math:`a_s` and :math:`m_q` in terms of
:math:`\kappa`.

.. container::
   :name: alpha_s

   .. rubric:: :math:`\alpha_s`
      :name: alpha_s

Let us start with :math:`\alpha_s`. Using the expansion:

.. math:: x=\frac{y}{1+b_1x+b_2x^2} = y\left[1-b_1x-(b_2-b_1^2)x^2 + \mathcal{O}(x^3)\right] = y \left[1-b_1y-(b_2-2b_1^2)y^2 + \mathcal{O}(y^3)\right]

we can invert eqs. :eq:`alphaspole` and
:eq:`alphasmsbar` obtaining:

.. math:: \frac{a^{(n)}(\mu_R)}{a^{(n-1)}(\mu_R)} = 1 + c_1a^{(n-1)}(\mu_R) + c_2 [a^{(n-1)}(\mu_R)]^2

where:

.. math::
   c_1 = \left\{
   \begin{array}{ll}
   \displaystyle \frac23 L_{\mu M}&\quad\mbox{Pole Mass}\\
   \\
   \displaystyle \frac23 L_{\mu m}&\quad\overline{\mbox{MS}}
   \end{array}\right.\quad\mbox{and}\quad
   c_2 = \left\{
   \begin{array}{ll}
   \displaystyle \frac49L_{\mu M}^2+\frac{38}3L_{\mu M}+\frac{14}3 &\quad\mbox{Pole Mass}\\
   \\
   \displaystyle \frac49L_{\mu m}^2+\frac{22}3L_{\mu m}-\frac{22}9&\quad\overline{\mbox{MS}}
   \end{array}\right.

Now, setting :math:`\mu_F=\kappa \mu_F` we have that:

.. math:: L_{\mu M} = \ln\frac{\mu_R}{M}=\ln\frac{\kappa\mu_F}{M}\quad\mbox{and}\quad L_{\mu m} = \ln\frac{\mu_R}{m(\mu_R)}=\ln\frac{\kappa\mu_F}{m(\kappa \mu_F)}\,.

In the case of the pole mass scheme, choosing :math:`\mu_F = M`, we have
directly that :math:`L_{\mu M}\rightarrow \ln\kappa` so that the full
matching condition reads:

.. math:: a^{(n-1)}(\kappa M)=a^{(n)}(\kappa M)\left\{1-\frac23\ln\kappa\,a^{(n)}(\kappa M)+\left(\frac49\ln^2\kappa-\frac{38}3\ln\kappa-\frac{14}3\right)[a^{(n)}(\kappa M)]^2\right\}

and:

.. math:: a^{(n)}(\kappa M)=a^{(n-1)}(\kappa M)\left\{1+\frac23\ln\kappa\,a^{(n-1)}(\kappa M)+\left(\frac49\ln^2\kappa+\frac{38}3\ln\kappa+\frac{14}3\right)[a^{(n-1)}(\kappa M)]^2\right\}\,.

In the case of the :math:`\overline{\mbox{MS}}` mass scheme, instead,
one chooses :math:`\mu_F = m(m)`, so that:

.. math:: L_{\mu m} \rightarrow \ln\kappa + \ln\frac{m(m)}{m(\kappa m)}\,.

From eq. :eq:`numsol` we see that:

.. math::
   :label: chezebedei
   \ln\frac{m(m)}{m(\kappa m)}=\int_{a_s(\kappa m)}^{a_s(m)}\frac{\gamma_m(a_s)}{\beta(a_s)}da_s\,,

and integrating, from eq. :eq:`integral`, we have that:

.. math::
   :label: hdsfgajfh
   \int_{a_s(\kappa m)}^{a_s(m)}\frac{\gamma_m(a_s)}{\beta(a_s)}da_s = c_0\ln\frac{a_s(m)}{a_s(\kappa m)} + (c_1-b_1c_0)[a_s(m)-a_s(\kappa m)]+\dots\,.

But from the perturbative expansion of the running of :math:`a_s` we
have:

.. math:: a_s(m) = a_s(\kappa m)\left[1-a_s(\kappa m)\beta_0\ln\kappa\right]\quad\Rightarrow\quad a_s(m)-a_s(\kappa m) = \mathcal{O}[a_s^2(\kappa m)]\,,

therefore the second term in the braket of the r.h.s. of eq.
:eq:`hdsfgajfh`, being of order :math:`a_s^2(\kappa m)`,
can be neglected because it would contribute only to the term
proportional to :math:`a_s^3(\kappa m)`. On the other hand:

.. math:: \ln\frac{a_s(m)}{a_s(\kappa m)} = \ln\left[1-a_s(\kappa m)\beta_0\ln\kappa \right]=a_s(\kappa m)\beta_0\ln\kappa+\mathcal{O}[a_s^2(\kappa m)]\,.

At the end of the day we find:

.. math:: \ln\frac{m(m)}{m(\kappa m)}=a_s(\kappa m)\gamma_m^{(0)}\ln\kappa+\mathcal{O}[a_s^2(\kappa m)]

so that:

.. math::
   :label: bugni
   L_{\mu m} \rightarrow \ln\kappa[1+\gamma_m^{(0)}a_s(\kappa m)]\,.

In the above equation, since
:math:`a_s^{(n-1)}=a_s^{(n)}+\mathcal{O}([a_s^{(n)}]^2)`, it doesn’t
matter whether one puts :math:`a_s^{(n)}(\kappa m)` or
:math:`a_s^{(n-1)}(\kappa m)` because in any case the difference would
be subleading.

Therefore, setting :math:`\mu=\mu_R=\kappa m(m) = \kappa m` into eq.
:eq:`alphasmsbar` and plugging eq.
:eq:`bugni`, one gets:

.. math:: a^{(n-1)}(\kappa m)=a^{(n)}(\kappa m)\left\{1-\frac23 \ln\kappa\,a^{(n)}(\kappa m)+\left[\frac49\ln^2\kappa-\frac{2}3\left(\gamma_m^{(0)}+11\right)\ln\kappa+\frac{22}9\right][a^{(n)}(\kappa m)]^2\right\}\,,

whose inverse is:

.. math:: a^{(n)}(\kappa m)=a^{(n-1)}(\kappa m)\left\{1+\frac23 \ln\kappa\,a^{(n-1)}(\kappa m)+\left[\frac49\ln^2\kappa+\frac{2}3\left(\gamma_m^{(0)}+11\right)\ln\kappa-\frac{22}9\right][a^{(n-1)}(\kappa m)]^2\right\}\,.

.. container::
   :name: m_q

   .. rubric:: :math:`m_q`
      :name: m_q

Now let us turn to :math:`m_q`. In this case everything is much easier.
First of all, we work only in the :math:`\overline{\mbox{MS}}` scheme,
secondly, given that also for an arbitary matching point the matching
condition for the running of the :math:`\overline{\mbox{MS}}` mass
starts at :math:`\mathcal{O}(\alpha_s^2)` (cfr. eq.
:eq:`mqmc`), writing :math:`L_{\mu m}` in terms of
:math:`\ln\kappa` would give rise to subleading terms (see eq.
:eq:`bugni`). It turns out that the matching condition for
the running of the :math:`\overline{\mbox{MS}}` mass in terms of
:math:`\ln\kappa` looks like this:

.. math:: m_q^{(n-1)}(\kappa m_h)=\left[1+\left(\frac43\ln^2\kappa-\frac{20}9\ln\kappa+\frac{89}{27}\right)[a^{(n)}(\kappa m_h)]^2\right]m_q^{(n)}(\kappa m_h)=\zeta_n^{(\mbox{\tiny down})}(\kappa m_h)m_q^{(n)}(\kappa m_h)

and the inverse is:

.. math:: m_q^{(n)}(\kappa m_h)=\left[1-\left(\frac43\ln^2\kappa-\frac{20}9\ln\kappa+\frac{89}{27}\right)[a^{(n-1)}(\kappa m_h)]^2\right]m_q^{(n-1)}(\kappa m_h)=\zeta_n^{(\mbox{\tiny up})}(\kappa m_h)m_q^{(n-1)}(\kappa m_h)\,.

Structure Functions
===================

.. container::
   :name: neutral-current

   .. rubric:: Neutral Current
      :name: neutral-current

In this section we discuss the explicit substitution of the
:math:`\overline{\mbox{MS}}` mass in the NC massive structure functions
(:math:`F_2` and :math:`F_L`). In our notation we define:

.. math:: M =\;\mbox{pole mass},\quad m\equiv m(\mu) =\;\overline{\mbox{MS}}\mbox{ mass},\quad a_s\equiv a_s(\mu),\quad h^{(l)}\equiv h^{(l)}(\mu,m(\mu))\,.

Dropping all the unnecessary dependences, the NC massive structure
function up to :math:`\mathcal{O}(a_s^2)` has the form:

.. math:: F(M) = a_sF^{(0)}(M) + a_s^2F^{(1)}(M) + \mathcal{O}(a_s^3)\,.

\ Now we want to explicitly replace the pole mass :math:`M` with the
:math:`\overline{\mbox{MS}}` mass :math:`m` using eq.
:eq:`inverse`, which in short reads:

.. math:: M = m(1 + a_sh^{(1)}) + \mathcal{O}(a_s^2)

To this end we expand :math:`F^{(0)}(M)` and :math:`F^{(1)}(M)` around
:math:`M=m` using the Taylor series in thy way:

.. math:: F^{(l)}(M) = \sum_{n=0}^{\infty}\frac1{n!}\frac{d^n F^{(l)}}{dM^n}\bigg|_{M=m}(M-m)^n\,,

so that, up to :math:`\mathcal{O}(a_s^2)`, what we need is:

.. math::
   \begin{array}{l}
   \displaystyle F^{(0)}(m) = F^{(0)}(m) + \frac{dF^{(0)}}{dM}\bigg|_{M=m}\underbrace{(M-m)}_{a_smh^{(1)}} = F^{(0)}(m) + a_smh^{(1)}\frac{dF^{(0)}}{dM}\bigg|_{M=m}\\
   \\
   \displaystyle F^{(1)}(M) = F^{(1)}(m)
   \end{array}\,.

Finally we have that:

.. math::
   :label: changescheme
   F(m) = a_sF^{(0)}(m) + a_s^2\left[F^{(1)}(m)+mh^{(1)}\frac{dF^{(0)}}{dM}\bigg|_{M=m}\right]\,.

In order to implement this structure function, we need to evaluate
explicitly the derivative in eq. :eq:`changescheme`.
First of all we observe that:

.. math::
   :label: pippo
   F^{(0)}(M) = x\int_x^{x_{\mbox{\tiny max}}(M)} \frac{dz}{z}g\left(\frac{x}{z}\right)C_g^{(0)}(\eta(z,M),\xi(M),\chi(M))\,,

where we have defined:

.. math:: x_{\mbox{\tiny max}}(M)=\frac1{1+\frac{4M^2}{Q^2}},\quad\eta(z,M) = \frac{Q^2}{4M^2}\left(\frac1z-1\right)-1,\quad \xi(M) =\frac{Q^2}{M^2},\quad \chi(M) =\frac{\mu^2}{M^2}\,.

But defining:

.. math:: G(z,M)=\frac{x}{z}g\left(\frac{x}{z}\right)C_g^{(0)}(\eta(z,M),\xi(M),\chi(M))\,,

eq. :eq:`pippo` can be written as:

.. math:: F^{(0)}(M) = \int_x^{x_{\mbox{\tiny max}}(M)} dz\,G(z,M)\,.

Therefore:

.. math::
   \begin{array}{c}
   \displaystyle \frac{dF^{(0)}}{dM} = \frac{d}{dM}\int_x^{x_{\mbox{\tiny max}}(M)} dzG(z,M) = \frac{d}{dM}\left[\widetilde{G}(x_{\mbox{\tiny max}}(M),M)-\widetilde{G}(x,M)\right]=\\
   \\
   \displaystyle \frac{d\widetilde{G}(x_{\mbox{\tiny max}}(M),M)}{dM}-\frac{d\widetilde{G}(x,M)}{dM}\,,
   \end{array}

where :math:`\widetilde{G}(z,M)` is the primitive of :math:`G(z,M)` with
respect of :math:`z` (i.e.
:math:`\partial\widetilde{G}/\partial z = G`). But:

.. math:: \frac{d\widetilde{G}(x_{\mbox{\tiny max}}(M),M)}{dM} = \frac{\partial \widetilde{G}(x_{\mbox{\tiny max}},M)}{\partial M}+\frac{dx_{\mbox{\tiny max}}}{dM}\underbrace{\frac{\partial \widetilde{G}(x_{\mbox{\tiny max}},M)}{\partial x_{\mbox{\tiny max}}}}_{G(x_{\mbox{\tiny max}},M)}

and:

.. math:: \frac{d \widetilde{G}(x,M)}{d M}=\frac{\partial\widetilde{G}(x,M)}{\partial M}\,,

thus:

.. math::
   :label: ignazio
   \begin{array}{c}
   \displaystyle \frac{dF^{(0)}}{dM} = \frac{\partial \widetilde{G}(x_{\mbox{\tiny max}},M)}{\partial M}-\frac{\partial\widetilde{G}(x,M)}{\partial M}+\frac{dx_{\mbox{\tiny max}}}{dM}G(x_{\mbox{\tiny max}},M) =\\
   \\
   \displaystyle \int_x^{x_{\mbox{\tiny max}}(M)} dz\frac{\partial G(z,M)}{\partial M}+\frac{dx_{\mbox{\tiny max}}}{dM}G(x_{\mbox{\tiny max}},M)\,.
   \end{array}

But in (Alekhin and Moch 2011) has been shown that the bounduary term in
eq. :eq:`ignazio` vanishes, thus it can be omitted.

Finally, since:

.. math:: \frac{\partial G(z,M)}{\partial M} = \frac{x}{z}g\left(\frac{x}{z}\right)\frac{\partial C_g^{(0)}}{\partial M}\,,

we have that:

.. math::
   :label: pollo
   \frac{dF^{(0)}}{dM}\bigg|_{M=m}=\left[x\int_x^{x_{\mbox{\tiny max}}(M)}\frac{dz}{z}g\left(\frac{x}{z}\right)\frac{\partial C_g^{(0)}}{\partial M}\right]\Bigg|_{M=m}=x\int_x^{x_{\mbox{\tiny max}}(m)}\frac{dz}{z}g\left(\frac{x}{z}\right)\left[\frac{\partial C_g^{(0)}}{\partial M}\right]\Bigg|_{M=m}\,.

Now, taking into account that:

.. math:: F^{(1)}(M) = \sum_{i=q,\overline{q},g}x\int_x^{x_{\mbox{\tiny max}}(M)} \frac{dz}{z}q_i\left(\frac{x}{z}\right)C_i^{(1)}(z,M)

and using eqs. :eq:`changescheme` and
:eq:`pollo`, one can explicitly write down the entire NNLO
massive structure function in terms of the :math:`\overline{\mbox{MS}}`
mass as follows:

.. math::
   :label: master
   \begin{array}{c}
   \displaystyle F(m) = x\int_x^{x_{\mbox{\tiny max}}(m)} \frac{dz}{z}g\left(\frac{x}{z}\right)\left[a_sC_g^{(0)}(z,m)+a_s^2\left(C_g^{(1)}(z,m)+mh^{(1)}\left[\frac{\partial C_g^{(0)}}{\partial M}\right]\Bigg|_{M=m}\right)\right]+\\
   \\
   \displaystyle \sum_{i=q,\overline{q}}x\int_x^{x_{\mbox{\tiny max}}(M)} \frac{dz}{z}q_i\left(\frac{x}{z}\right)a_s^2C_i^{(1)}(z,M)\,.
   \end{array}

We can now read the recipe for the implementation: in order to implement
the :math:`\mathcal{O}(a_s^2)` massive structure function (:math:`F_2`
or :math:`F_L`) in terms of the :math:`\overline{\mbox{MS}}` mass
:math:`m`, one has just to replace everywhere the pole mass :math:`M`
with :math:`m` and add to the :math:`\mathcal{O}(a_s^2)` gluon
coefficient function :math:`C_g^{(1)}(z,m)` the term:

.. math:: m(\mu)h^{(1)}(\mu,m(\mu))\left[\frac{\partial C_g^{(0)}}{\partial M}\right]\Bigg|_{M=m(\mu)}\,.

Of course, for the massless limit of the massive structure function
(massive0) the same recipe holds, with the only obvious difference that
one has to replace the massive coefficient functions with the massive0
ones.

Now we need to evaluate explicitly the derivative of :math:`C_g^{(0)}`
in eq. :eq:`master` and this must be done separately for
:math:`F_2` and :math:`F_L`.

.. container::
   :name: f_2

   .. rubric:: :math:`F_2`
      :name: f_2

We consider first :math:`F_2`. Since the NNPDF code works in the Mellin
space, it is better to calculate directly the derivative of the Mellin
transform of :math:`C_{2,g}^{(0)}`, which in the massive case is:

.. math::
   \begin{array}{rl}
   \displaystyle C_{2,g}^{(0)}(N,Q^2,M^2)=&\displaystyle T_R\Big\{2(1-6\epsilon-4\epsilon^2)I_2(a,N)-2(1-2\epsilon)I_1(a,N)+I_0(a,N)+\\
   \\
   & \displaystyle -4(2-\epsilon)J_2(a,N)+4(2-\epsilon)J_1(a,N)-J_0(a,N)\Big\}\,,
   \end{array}

\ where:

.. math::
   :label: iq
   I_q(a,N) = \frac{a^{N+q}}{N+q}\frac{\Gamma(N+q)\Gamma(\frac12)}{\Gamma(N+q+\frac12)} {_2F_1}\left(\frac12,N+q,N+q+\frac12;a\right)

.. math::
   :label: jq
   \begin{array}{rl}
   J_q(a,N) =&\displaystyle a^{N+q}\frac{\Gamma(N+q)\Gamma(\frac12)}{\Gamma(N+q+\frac12)}\bigg\{{_2F_1}\left(\frac12,N+q,N+q+\frac12;a\right)\\
   \\
   &\displaystyle -\frac{N+q}{N+q+\frac12}{_2F_1}\left(\frac12,N+q+1,N+q+\frac32;a\right)\bigg\}
   \end{array}\,,

with:

.. math::
   :label: definitions1
   \epsilon = \frac{M^2}{Q^2}\quad\mbox{and}\quad a=\frac1{1+4\epsilon}\,.

From the definitions in eq. :eq:`definitions1` we
obtain:

.. math::
   :label: gngngngngng
   \begin{array}{rl}
   \displaystyle \frac{\partial}{\partial M} &\displaystyle = \frac{\partial \epsilon}{\partial M} \frac{\partial}{\partial \epsilon} = \frac{2\epsilon}{M}\frac{\partial}{\partial \epsilon}\\
   \\
   \displaystyle \frac{\partial}{\partial M} &\displaystyle = \frac{\partial \epsilon}{\partial M} \frac{\partial a}{\partial \epsilon} \frac{\partial}{\partial a} = -\frac{8a^2\epsilon}{M}\frac{\partial}{\partial a}
   \end{array}\,.

Therefore:

.. math::
   \begin{array}{rl}
   \displaystyle \frac{\partial C_{2,g}^{(0)}}{\partial M}=&\displaystyle T_R\Bigg\{\frac{2\epsilon}M\Big[2(-6-8\epsilon)I_2+4I_1+4J_2-4J_1\Big]\\
   \\
   &\displaystyle -\frac{8a^2\epsilon}{M}\Bigg[2(1-6\epsilon-4\epsilon^2)\frac{\partial I_2}{\partial a}-2(1-2\epsilon)\frac{\partial I_1}{\partial a}+\frac{\partial I_0}{\partial a}\\
   \\
   & \displaystyle -4(2-\epsilon)\frac{\partial J_2}{\partial a}+4(2-\epsilon)\frac{\partial J_1}{\partial a}-\frac{\partial J_0}{\partial a}\Bigg]\Bigg\}\,.
   \end{array}

Now, starting from eqs. :eq:`iq` and :eq:`jq`, we need
to evaluate the derivative of :math:`I_q` and :math:`J_q` and we do it
using the relation valid for the derivative of the hypergeometric
function:

.. math:: \frac{\partial}{\partial x} {_2F_1}(a,b,c;x) = \frac{b}{x}\left[ {_2F_1}(a,b+1,c;x) - {_2F_1}(a,b,c;x)\right]\,,

we have that:

.. math::
   \begin{array}{ll}
   \displaystyle \frac{d}{da} a^{N+q} {_2F_1}\left(\frac12,N+q,N+q+\frac12;a\right) = &\displaystyle a^{N+q-1}(N+q) {_2F_1}\left(\frac12,N+q+1,N+q+\frac12;a\right)\\
   \\
   \displaystyle \frac{d}{da} a^{N+q} {_2F_1}\left(\frac12,N+q+1,N+q+\frac32;a\right) = &\displaystyle a^{N+q-1}(N+q+1) {_2F_1}\left(\frac12,N+q+2,N+q+\frac32;a\right)\\
   &\displaystyle -a^{N+q-1}{_2F_1}\left(\frac12,N+q+1,N+q+\frac32;a\right)\,,
   \end{array}

so that we get:

.. math::
   :label: diq
   \begin{array}{rl}
   \displaystyle \frac{\partial I_q}{\partial a} &\displaystyle = a^{N+q-1}\frac{\Gamma(N+q)\Gamma(\frac12)}{\Gamma(N+q+\frac12)}{_2F_1}\left(\frac12,N+q+1,N+q+\frac12;a\right)
   \end{array}

and:

.. math::
   :label: djq
   \begin{array}{rl}
   \displaystyle \frac{\partial J_q}{\partial a} &\displaystyle =  a^{N+q-1}\frac{\Gamma(N+q+1)\Gamma(\frac12)}{\Gamma(N+q+\frac12)}\bigg\{{_2F_1}\left(\frac12,N+q+1,N+q+\frac12;a\right)\\
   \\
   &\displaystyle - \frac{N+q+1}{N+q+\frac12}{_2F_1}\left(\frac12,N+q+2,N+q+\frac32;a\right)\\
   \\
   &\displaystyle +\frac1{N+q+\frac12}{_2F_1}\left(\frac12,N+q+1,N+q+\frac32;a\right)\bigg\}\,.
   \end{array}

Looking at these expressions, one can see that in these derivatives, a
part from hypergeometric functions of the form
:math:`{_2F_1(a,b,a+b;x)}` which were already present in
:math:`C_{2,g}^{(0)}` itself, also hypergeometric functions of the form
:math:`{_2F_1(a,b,a+b-1;x)}` appear. This raises a technical problem
because the NNPDF code uses a fast routine for the hypergeometric
function which is accurate both around :math:`x=0` and :math:`x=1`, but
with the limitation :math:`c=a+b`. Now, instead, we need also the case
:math:`{c=a+b-1}`, therefore we need to extend our routine including
this possibility. We can do this using the expansion around :math:`x=1`
reported in eq. (15.3.12) of (Abramowitz and Stegun 1964).

Now we consider the NC massive0 structure function :math:`F_2^0`. In
this limit the gluon coefficient function takes the form:

.. math::
   \begin{array}{c}
   \displaystyle C_{2,g}^{0,(0)}(N,Q^2,M^2)=\\
   \\
   \displaystyle T_R\Bigg[2\left(\ln\frac{Q^2}{M^2}-4\right)\frac1{N+2}-2\left(\ln\frac{Q^2}{M^2}-4\right)\frac1{N+1}+\left(\ln\frac{Q^2}{M^2}-1\right)\frac1N\\
   \\
   \displaystyle -2\frac{S_1(N+2)}{N+2}+2\frac{S_1(N+1)}{N+1}-\frac{S_1(N)}{N}+\frac2{(N+2)^2}-\frac2{(N+1)^2}+\frac1{N^2}\Bigg]\,.
   \end{array}

\ Therefore, considering that:

.. math:: \frac{\partial}{\partial M} \ln\frac{Q^2}{M^2} = - \frac2M\,,

the derivative of :math:`C_{2,g}^{0,(0)}` is given by:

.. math:: \frac{\partial C_{2,g}^{0,(0)}}{\partial M}= -T_R\frac{2}{M}\left[\frac2{N+2}-\frac2{N+1}+\frac1N\right]

.. container::
   :name: f_l

   .. rubric:: :math:`F_L`
      :name: f_l

Now we consider :math:`F_L`. In this case the Mellin transform of the
gluon coefficient function takes the simpler form:

.. math::
   :label: cgnf1L
   C_{L,g}^{(0)}\left(N,Q^2,M^2\right)= T_R\left[-8\epsilon I_2(a,N)-4J_2(a,N)+4J_1(a,N)\right]\,.

where :math:`I_q` and :math:`J_q` are given in eqs. :eq:`iq` and
:eq:`jq`, respectively. Therefore, using eq.
:eq:`gngngngngng`, we get:

.. math:: \frac{\partial C_{L,g}^{(0)}}{\partial M}= T_R\left\{-\frac{16\epsilon}{M} I_2 -\frac{8a^2\epsilon}{M}\left[-8\epsilon \frac{\partial I_2}{\partial a}-4\frac{\partial J_2}{\partial a}+4\frac{\partial J_1}{\partial a}\right]\right\}

where the derivatives of :math:`I_q` and :math:`J_q` with respect of
:math:`a` are given in eqs. :eq:`diq` and :eq:`djq`,
respectively.

The massive0 gluon coefficient function :math:`C_{L,g}^{0,(0)}`,
instead, turns out to be independent from :math:`M`. This means that:

.. math:: \frac{\partial C_{L,g}^{0,(0)}}{\partial M}= 0

Finally, having the derivative with respect of :math:`M` of the
:math:`\mathcal{O}(a_s)` gluon coefficient function for both :math:`F_2`
and :math:`F_L` in both the massive and massive0 schemes, we can plug it
into eq. :eq:`master` and obtain the neutral current
structure function in terms of the :math:`\overline{\mbox{MS}}` mass
:math:`m`.

.. container::
   :name: charged-current

   .. rubric:: Charged Current
      :name: charged-current

In this section we consider the generic CC massive structure function.
The treatment is exactly the same of the NC structure functions, with
the only difference that the CC case they start at
:math:`\mathcal{O}(a_s^0)` and they are presently known up to
:math:`\mathcal{O}(a_s)`. This means that their perturbative expansion
in terms of the pole mass :math:`M` looks like this:

.. math:: F_k(M) = F_k^{(0)}(M) + a_sF_k^{(1)}(M) + \mathcal{O}(a_s^2)\,,

\ with :math:`k=2,3,L`. Therefore, expanding :math:`F^{(0)}` and
:math:`F^{(1)}` around :math:`M=m` and keeping only the terms up to
:math:`\mathcal{O}(a_s)`, one obtains:

.. math::
   :label: gigi
   F_k(m) = F_k^{(0)}(m) + a_s\left[F_k^{(1)}(m)+mh^{(1)}\frac{dF_k^{(0)}}{dM}\bigg|_{M=m}\right]\,.

One can show that:

.. math:: F^{(0)}_k(M) = b_k(M)s'(\xi(M))\,,

\ where:

.. math::
   :label: definitions
   \xi = x\underbrace{\left(1+\frac{M^2}{Q^2}\right)}_{\frac1\lambda}=\frac{x}\lambda\quad\mbox{and}\quad
   \left\{\begin{array}{l}
   b_2 = \xi\\
   b_3 = 1\\
   b_L = (1-\lambda)\xi
   \end{array}
   \right.

and where we have also defined:

.. math:: s'\equiv 2|V_{cs}|^2s+2|V_{cd}|^2[f\,d+(1-f)u]\quad\mbox{with}\quad f=\frac{N_p}{N_p+N_n}\,.

Therefore:

.. math::
   :label: pincopanco
   mh^{(1)}\frac{dF^{(0)}_k}{dM}\bigg|_{M=m} = mh^{(1)}\frac{d\xi}{dM}\frac{dF^{(0)}_k}{d\xi}\bigg|_{M=m} = 2h^{(1)}(1-\lambda)\xi\left[\frac{db_k}{d\xi}s'(\xi)+b_k(\xi)\frac{ds'}{d\xi}\right]\bigg|_{M=m}\,,

that can be conveniently rewritten as:

.. math:: mh^{(1)}\frac{dF^{(0)}_k}{dM}\bigg|_{M=m} = 2h^{(1)}(1-\lambda)\left[\left(\frac{db_k}{d\xi}-\frac{b_k}{\xi}\right)+b_k(\xi)\frac{d}{d\xi}\right]\xi s'(\xi)\bigg|_{M=m}\,,

so that, using eq. :eq:`definitions`, we have that:

.. math::
   \begin{array}{rcl}
   \displaystyle mh^{(1)}\frac{dF^{(0)}_2}{dM}\bigg|_{M=m} &=& \displaystyle
   2h^{(1)}(1-\lambda)\xi\frac{d}{d\xi} \xi s'(\xi)\bigg|_{M=m}\\
   \\
   \displaystyle mh^{(1)}\frac{dF^{(0)}_3}{dM}\bigg|_{M=m} &=& \displaystyle
   2h^{(1)}(1-\lambda)\frac{1}{\xi}\left[ \xi\frac{d}{d\xi}-1\right]\xi
   s'(\xi)\bigg|_{M=m}\\
   \\
   \displaystyle mh^{(1)}\frac{dF^{(0)}_L}{dM}\bigg|_{M=m} &=& \displaystyle
   2h^{(1)}(1-\lambda)^2\xi \frac{d}{d\xi}\xi s'(\xi)\bigg|_{M=m}
   \end{array}\,.

Eqs. :eq:`pincopanco`, though apparently very easy,
involve the derivative of the PDF :math:`s'` and this makes the
implementation a little bit more troblesome.

Using the same arguments of eq. :eq:`ignazio`, one can show
that:

.. math:: \frac{ds'}{d\xi} = \frac{d}{d\xi}\int_\xi^1\frac{dy}{y}\delta(1-y)s'\left(\frac{\xi}{y}\right) = \int_\xi^1\frac{dy}{y}\delta(1-y)\frac{d}{d\xi}s'\left(\frac{\xi}{y}\right)\,.

therefore eq. :eq:`pincopanco` can be written as:

.. math::
   :label: pancopinco
   mh^{(1)}\frac{dF^{(0)}_k}{dM}\bigg|_{M=m} = 2h^{(1)}(1-\lambda)\xi\int_\xi^1\frac{dy}{y}\delta(1-y)\left[\frac{db_k}{d\xi}+b_k(\xi)\frac{d}{d\xi}\right]s'\left(\frac{\xi}{y}\right)\,,

where in the r.h.s. we are understanding that the pole mass :math:`M`,
which appears only through :math:`\xi`, must be replaced everywhere with
the :math:`\overline{\mbox{MS}}` mass :math:`m`. But since:

.. math:: \frac{d}{dx}f\left(\frac{x}{y}\right) = \frac1{y}\frac{d}{d\left(\frac{x}{y}\right)}f\left(\frac{x}{y}\right) = \frac1{xy}\frac{d}{d\left(\frac{1}{y}\right)}f\left(\frac{x}{y}\right)

and:

.. math:: d\left(\frac1y\right) = -\frac1{y^2}dy\quad\Rightarrow \quad \frac{d}{dx}f\left(\frac{x}{y}\right) = -\frac{y}{x}\frac{d}{dy}f\left(\frac{x}{y}\right)\,,

it follows that:

.. math:: \int_\xi^1\frac{dy}y\delta(1-y)\frac{d}{d\xi}s'\left(\frac{\xi}{y}\right) = -\int_\xi^1\frac{dy}\xi\delta(1-y)\frac{d}{dy}s'\left(\frac{\xi}{y}\right)\,.

Now, integrating by parts the r.h.s. of the equation above, one gets:

.. math:: \int_\xi^1\frac{dy}y\delta(1-y)\frac{d}{d\xi}s'\left(\frac{\xi}{y}\right)=\int_\xi^1\frac{dy}y\left\{\frac{y}{\xi}\left[\frac{d}{dy}\delta(1-y)\right]\right\}s'\left(\frac{\xi}{y}\right)\,.

Therefore eq. :eq:`pancopinco` can we rewritten as:

.. math:: mh^{(1)}\frac{dF^{(0)}_k}{dM}\bigg|_{M=m} = 2h^{(1)}(1-\lambda)\xi\int_\xi^1\frac{dy}{y}\left\{\frac{db_k}{d\xi}\delta(1-y)+\frac{b_k(\xi)}{\xi}\left[y\frac{d}{dy}\delta(1-y)\right]\right\}s'\left(\frac{\xi}{y}\right)\,.

In the above equation a coefficient function can be isolated and,
considering the form of :math:`b_k` given in eq.
:eq:`definitions`, we write:

.. math::
   :label: RoughCF
   \begin{array}{l}
   \displaystyle \widetilde{C}_{2,q}(y)=2h^{(1)}(1-\lambda)\left\{\delta(1-y)+\left[y\frac{d}{dy}\delta(1-y)\right]\right\}\\
   \\
   \displaystyle \widetilde{C}_{3,q}(y)=2h^{(1)}(1-\lambda)\left[y\frac{d}{dy}\delta(1-y)\right]\\
   \\
   \displaystyle \widetilde{C}_{L,q}(y)=2h^{(1)}(1-\lambda)^2\left\{\delta(1-y)+\left[y\frac{d}{dy}\delta(1-y)\right]\right\}
   \end{array}\,,

in such a way that:

.. math::
   \begin{array}{l}
   \displaystyle mh^{(1)}\frac{dF^{(0)}_k}{dM}\bigg|_{M=m} = \xi \int_\xi^1\frac{dy}{y}\widetilde{C}_{k,q}(y)s'\left(\frac{\xi}{y}\right)\\
   \\
   \displaystyle mh^{(1)}\frac{dF^{(0)}_3}{dM}\bigg|_{M=m} = \int_\xi^1\frac{dy}{y}\widetilde{C}_{3,q}(y)s'\left(\frac{\xi}{y}\right)
   \end{array}\,,

where now :math:`k=2,L` and whose Mellin transforms, taking into account
that:

.. math:: \mathbf{M}\left[y\frac{d}{dy}\delta(1-y)\right](N) = - N\,,

can be easly evaluated obtaining( [6]_):

.. math::
   :label: greatresult
   \begin{array}{l}
   \displaystyle \widetilde{C}_{2,q}(N)=\mathbf{M}[\widetilde{C}_{2,q}(y)](N)=2h^{(1)}(1-\lambda)(1-N)\\
   \\
   \displaystyle \widetilde{C}_{3,q}(N)=\mathbf{M}[\widetilde{C}_{3,q}(y)](N)=-2h^{(1)}(1-\lambda)N\\
   \\
   \displaystyle \widetilde{C}_{L,q}(N)=\mathbf{M}[\widetilde{C}_{L,q}(y)](N)=2h^{(1)}(1-\lambda)^2(1-N)
   \end{array}\,.

In order to carry out the :math:`x`-space implementation, one can show
that:

.. math::
   :label: DeltaDerivative
   \frac{d}{dy}\delta(1-y) = \left[\frac{\delta(1-y)}{1-y}\right]_+\,,

which is a pretty formal expression that however helps in manipulating
the coefficient functions in the presence of
:math:`\overline{\mbox{MS}}` masses. In fact, using
eq. :eq:`DeltaDerivative`, one can easily show
that:

.. math:: y\frac{d}{dy}\delta(1-y) = \left[\frac{\delta(1-y)}{1-y}\right]_+ - \delta(1-y)\,,

so that eqs. :eq:`RoughCF` become:

.. math::
   \begin{array}{l}
   \displaystyle \widetilde{C}_{2,q}(y)=2h^{(1)}(1-\lambda) \left[\frac{\delta(1-y)}{1-y}\right]_+\\
   \\
   \displaystyle \widetilde{C}_{3,q}(y)=2h^{(1)}(1-\lambda)\left\{\left[\frac{\delta(1-y)}{1-y}\right]_+ - \delta(1-y)\right\}\\
   \\
   \displaystyle \widetilde{C}_{L,q}(y)=2h^{(1)}(1-\lambda)^2 \left[\frac{\delta(1-y)}{1-y}\right]_+
   \end{array}\,.

Now, since in :eq:`gigi` we have that:

.. math:: F_k^{(1)}(m)=\xi\int_{\xi}^{1}\frac{dy}{y}\left\{C_{k,q}(y)s'\left(\frac{\xi}{y}\right)+C_{k,g}(y)g\left(\frac{\xi}{y}\right)\right\}

.. math:: F_3^{(1)}(m)=\int_{\xi}^{1}\frac{dy}{y}\left\{C_{3,q}(y)s'\left(\frac{\xi}{y}\right)+C_{3,g}(y)g\left(\frac{\xi}{y}\right)\right\}\,.

This means that the whole :math:`\mathcal{O}(a_s)` in eq.
:eq:`gigi` can be written as:

.. math:: F_k^{(1)}(m)+mh^{(1)}\frac{dF_k^{(0)}}{dM}\bigg|_{M=m}=\xi\int_{\xi}^{1}\frac{dy}{y}\left\{\left[C_{k,q}(y)+\widetilde{C}_{k,q}(y)\right]s'\left(\frac{\xi}{y}\right)+C_{k,g}(y)g\left(\frac{\xi}{y}\right)\right\}

.. math:: F_3^{(1)}(m)+mh^{(1)}\frac{dF_3^{(0)}}{dM}\bigg|_{M=m}=\int_{\xi}^{1}\frac{dy}{y}\left\{\left[C_{3,q}(y)+\widetilde{C}_{3,q}(y)\right]s'\left(\frac{\xi}{y}\right)+C_{3,g}(y)g\left(\frac{\xi}{y}\right)\right\}\,.

Therefore, in order to consistently replace the pole mass :math:`M` with
the :math:`\overline{\mbox{MS}}` mass :math:`m` in the charge current
massive coefficient functions, one has just to naively replace :math:`M`
with :math:`m` and then correct the :math:`\mathcal{O}(a_s)` quark
coefficient functions adding (in the Mellin space) the contributions
given in eq. :eq:`greatresult`.

It is interesting to observe that in massless limit, where
:math:`\lambda\rightarrow 1`, all the coefficient functions in eq.
:eq:`greatresult` vanish, with the consequence that the
CC massive0 structure functions up to :math:`\mathcal{O}(a_s)` in terms
of :math:`M` or :math:`m` look exactly the same.

Thresholds
==========

The evolution schemes that are involved in the FONLL scheme, meaning ZM
(Zero-Mass-Variable-Flavour-Number) and FFN (massive) schemes, require
the presence of a mass threshold for each heavy flavour. These
thresholds are basically the points from where the respective heavy
quark structure functions start contributing to the total structure
function.

In the ZM scheme, these thresholds don’t have a uniquely defined
physical meaning, but rather they just represent a convenient choice of
:math:`Q^2` where to perform the matching between the :math:`N_f` and
the :math:`N_f+1` scheme. If one chooses to write the observables in
terms of the pole masses :math:`M_c`, :math:`M_b` and :math:`M_t`, the
most natural choice for the thresholds are the pole masses themselves,
so that ZM structure function can be written as:

.. math:: F^{(zm)}(x,Q^2) = F^{(zm),l}(x,Q^2)+\sum_{i=c,b,t}\theta(Q^2-M_i^2)F^{(zm),i}(x,Q^2,M_i)\,.

This is also justified by the fact that :math:`\alpha_s` as well as PDFs
are conveniently matched in correspondence of these thresholds.

When one instead chooses to write the observables in terms of the
:math:`\overline{\mbox{MS}}` masses, as we have seen, the most covenient
choice for the matching thresholds of the PDFs, :math:`\alpha_s` and
masses evolution are the RG-invariant masses :math:`m_c(m_c)`,
:math:`m_b(m_b)` and :math:`m_t(m_t)`. For this reason, in the
:math:`\overline{\mbox{MS}}` framework, it looks more natural to choose
the same thresholds also for the structure functions, so that in term of
the :math:`\overline{\mbox{MS}}` masses, the generic ZM structure
function looks like this:

.. math:: F^{(zm)}(x,Q^2) = F^{(zm),l}(x,Q^2)+\sum_{i=c,b,t}\theta(Q^2-m_i^2(m_i))F^{(zm),i}(x,Q^2,m_i(\mu))\,.

This means that now the :math:`i`-th heavy quark structure function
switchs on at the scale :math:`Q^2=m_i(m_i)` rather than
:math:`Q^2=M_i`.

In the FFN scheme, instead, the heavy quark mass thresholds assume a
precise physical meaning. In fact, they tell us whether the invariant
mass of the incoming particles (the photon and the parton)
:math:`W=\sqrt{Q^2(1-x)/x}` is big enough for producing (up to NNLO)
one, in the CC case, or two, in the NC, heavy quarks. In terms of the
pole masses :math:`M_c`, :math:`M_b` and :math:`M_t`, the kinematical
threshold for producing :math:`i`-th species of heavy quarks is given
by:

.. math::
   W^2\geq \kappa M_i^2\quad\mbox{with}\quad
   \left\{\begin{array}{l}
   \kappa = 4\quad\mbox{for NC}\\
   \kappa = 1\quad\mbox{for CC}
   \end{array}\right.\,,

\ so that FFN structure function can be written as:

.. math:: F^{(m)}(x,Q^2) = F^{(m),l}(x,Q^2)+\sum_{i=c,b,t}\theta(W^2-\kappa M_i^2)F^{(m),i}(x,Q^2,M_i)\,.

Now the question is: what are the right thresholds if we write this
structure function in terms of the :math:`\overline{\mbox{MS}}` masses
:math:`m_c(\mu)`, :math:`m_b(\mu)` and :math:`m_t(\mu)`? One more time,
the most natural choice seems to be the RG-invariant masses
:math:`m_c(m_c)`, :math:`m_b(m_b)` and :math:`m_t(m_t)`. The reason is
the following.

Given that the :math:`\overline{\mbox{MS}}` masses run, we are
interested in knowing the value of the heavy quark mass :math:`m_i(\mu)`
when the scale of the process :math:`Q^2\simeq m_i^2`, but since in any
case :math:`\mu^2\simeq Q^2`, it seems to be the most reasonable choice
to take as a threshold the value :math:`m_i(m_i)`. This means that in
terms of :math:`\overline{\mbox{MS}}` masses, the FFN structure function
takes the form:

.. math:: F^{(m)}(x,Q^2) = F^{(m),l}(x,Q^2)+\sum_{i=c,b,t}\theta(W^2-\kappa m_i^2(m_i))F^{(m),i}(x,Q^2,m_i)\,.

As a conclusion, also in the massive case it turns out to be convenient
to replace at the thresholds the pole masses :math:`M_c`, :math:`M_b`
and :math:`M_t` with the :math:`\overline{\mbox{MS}}` RG-invariant
masses :math:`m_c(m_c)`, :math:`m_b(m_b)` and :math:`m_t(m_t)`.

**References**

**References**

.. container:: references csl-bib-body hanging-indent
   :name: refs

   .. container:: csl-entry
      :name: ref-AbramowitzStegun

      Abramowitz, Milton, and Irene A. Stegun. 1964. *Handbook of
      Mathematical Functions with Formulas, Graphs, and Mathematical
      Tables*. Ninth Dover printing, tenth GPO printing. New York:
      Dover.

   .. container:: csl-entry
      :name: ref-Alekhin:2010sv

      Alekhin, S., and S. Moch. 2011. “Heavy-quark deep-inelastic
      scattering with a running mass.” *Phys. Lett. B* 699: 345–53.
      https://doi.org/10.1016/j.physletb.2011.04.026.

   .. container:: csl-entry
      :name: ref-Buza:1996wv

      Buza, M., Y. Matiounine, J. Smith, and W. L. van Neerven. 1998.
      “Charm electroproduction viewed in the variable flavor number
      scheme versus fixed order perturbation theory.” *Eur. Phys. J. C*
      1: 301–20. https://doi.org/10.1007/BF01245820.

   .. container:: csl-entry
      :name: ref-Chetyrkin:1997sg

      Chetyrkin, K. G., Bernd A. Kniehl, and M. Steinhauser. 1997.
      “Strong coupling constant with flavor thresholds at four loops in
      the MS scheme.” *Phys. Rev. Lett.* 79: 2184–87.
      https://doi.org/10.1103/PhysRevLett.79.2184.

   .. container:: csl-entry
      :name: ref-Chetyrkin:2000yt

      Chetyrkin, K. G., Johann H. Kuhn, and M. Steinhauser. 2000.
      “RunDec: A Mathematica package for running and decoupling of the
      strong coupling and quark masses.” *Comput. Phys. Commun.* 133:
      43–65. https://doi.org/10.1016/S0010-4655(00)00155-7.

   .. container:: csl-entry
      :name: ref-Chetyrkin:1999pq

      Chetyrkin, K. G., and A. Retey. 2000. “Renormalization and running
      of quark mass and field in the regularization invariant and MS-bar
      schemes at three loops and four loops.” *Nucl. Phys. B* 583: 3–34.
      https://doi.org/10.1016/S0550-3213(00)00331-X.

   .. container:: csl-entry
      :name: ref-Chetyrkin:1999qi

      Chetyrkin, K. G., and M. Steinhauser. 2000. “The Relation between
      the MS-bar and the on-shell quark mass at order alpha(s)**3.”
      *Nucl. Phys. B* 573: 617–51.
      https://doi.org/10.1016/S0550-3213(99)00784-1.

   .. container:: csl-entry
      :name: ref-Melnikov:2000qh

      Melnikov, Kirill, and Timo van Ritbergen. 2000. “The Three loop
      relation between the MS-bar and the pole quark masses.” *Phys.
      Lett. B* 482: 99–108.
      https://doi.org/10.1016/S0370-2693(00)00507-4.

   .. container:: csl-entry
      :name: ref-Vogt:2004ns

      Vogt, A. 2005. “Efficient evolution of unpolarized and polarized
      parton distributions with QCD-PEGASUS.” *Comput. Phys. Commun.*
      170: 65–92. https://doi.org/10.1016/j.cpc.2005.03.103.

.. [1]
   As a consistency check, note that setting :math:`\mu^2=M^2` and
   taking into account the fact that :math:`\zeta_2 = \pi^2/6`, the
   coefficients in eq. :eq:`coeff` reduce to the first five
   coefficients in the equation between eq. (10) and (11) of (Melnikov
   and Ritbergen 2000).

.. [2]
   Note that, to be consistent, the evaluation of :math:`a` and
   :math:`a_0` must be done at the same perturbative order of
   :math:`m(\mu)`. So, for instance, if we want to evaluate the NNLO
   running of :math:`m(\mu)` also the value of :math:`a` and :math:`a_0`
   must be computed using the NNLO running.

.. [3]
   It should be noticed that when considering the FFNS the number of
   active flavours stays the same for all scales. In particular, given
   that in the approach discussed above the number of active flavours
   for each of the input masses :math:`m_i(m_i)` is assumed to be equal
   to the number of flavours right above the tresholds (:math:`e.g.`
   :math:`m_c^{(4)}(m_c)`), this is not the same parameter as
   :math:`m_c^{(3)}(m_c)` that instead would be used in the
   :math:`N_f=3` FFNS. In fact, beyond NLO, due to the presence of the
   matching conditions, they differ by :math:`\mathcal{O}(\alpha_s)`
   terms.

.. [4]
   It should be noticed that :math:`a_s` enetring
   eq. :eq:`couple1`

.. [5]
   It should be noticed that in the case :math:`\kappa\neq1` PDFs aquire
   an implicit dependence on :math:`\mu_R` that essentially comes from
   the redefinition of the splitting functions that in turn derives from
   the expansion of :math:`\alpha_s(\mu_R)` around :math:`\mu_R=\mu_F`
   that appears in the DGLAP equation.

.. [6]

   .. container::
      :name: alternative-calculation

      .. rubric:: Alternative Calculation
         :name: alternative-calculation
         :class: unnumbered

   We sketch here an alternative calculation that, under some point of
   view, looks more transparent ant confirms the result found in eq.
   :eq:`greatresult`. We start directly calculating the
   Mellin transform of the term proportional to the derivative of eq.
   :eq:`pincopanco`, that is, disregarding the overall
   constant:

   .. math:: I_k(N)=\int_0^1d\xi \xi^{N-1}\left[\xi b_k(\xi)\frac{ds'}{d\xi}\right]\,.

   Bu since:

   .. math:: \frac{d}{d\xi}\xi^{N}b_k(\xi)s'(\xi)=\left[\frac{d}{d\xi}\xi^{N}b_k(\xi)\right]s'(\xi)+\xi^{N}b_k(\xi)\frac{ds'}{d\xi}\,,

   it follows that:

   .. math:: I_k(N)=\underbrace{\xi^{N}b_k(\xi)s'(\xi)\Big|_0^1}_{=0}-\int_0^1d\xi\left[\frac{d}{d\xi}\xi^{N}b_k(\xi)\right]s'(\xi)\,.

   Now, using the definition of :math:`b_k(\xi)` given in eq.
   :eq:`definitions`, we can easily find that:

   .. math::
      \begin{array}{l}
        \displaystyle I_2(N) = -(N+1)s'(N+1)\\
        \displaystyle I_L(N) = -(1-\lambda)(N+1)s'(N+1)\\
        \displaystyle I_3(N) = -Ns'(N)
        \end{array}

   from which one can extract the coefficient functions in the Mellin
   space.
