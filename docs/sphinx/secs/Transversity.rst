==========================
Transversity Distributions
==========================

.. contents::
   :depth: 3
..

Perturbative evolution
======================

In this section we discuss the structure of the DGLAP evolution
equations for the transversity distributions. The same structure holds
for both PDFs and FFs. Therefore we first discuss the structure of the
evolution equations in terms of distributions in the so-called
“evolution” basis and then report the splitting functions up to
:math:`\mathcal{O}(\alpha_s^2)`, *i.e.* next-to-leading order (NLO),
separately for PDFs and FFs. Contrary to unpolarised and longitudinally
polarised collinear distributions, no transversely polarised gluon
distribution exists. This simplifies the structure of the evolution
equations that, when written in the evolution basis, are completely
decoupled.

As a first step we define the evolution basis. Given a set of quark
distributions in the more familiar “physical” basis, *i.e.*
:math:`\{\overline{t},\overline{b},\overline{c},\overline{s},\overline{u},\overline{d},d,u,s,c,b,t\}`,
and defining :math:`q^{\pm}\equiv q\pm \overline{q}`, the “evolution”
basis is defined as follows:

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

\ It is possible to show that in this basis the general form of the
DGLAP evolution equations for the transversity distribution reduces to
the following set of *decoupled* integro-differential equation:

.. math::
   :label: eq:evoleqs
   \begin{array}{rcl}
   \displaystyle \mu^2\frac{d\Sigma}{d\mu^2} &=& \displaystyle P_{qq}\otimes \Sigma\,,\\
   \\
   \displaystyle \mu^2\frac{dV}{d\mu^2} &=& \displaystyle P^{V}\otimes V\,,\\
   \\
   \displaystyle \mu^2\frac{dT_i}{d\mu^2} &=& \displaystyle P^+\otimes T_i\,,\\
   \\
   \displaystyle \mu^2\frac{dV_i}{d\mu^2} &=& \displaystyle P^-\otimes V_i\,,\\
   \\
   \end{array}

with :math:`i=3,8,15,24,35` and where the Mellin convolution symbol
:math:`\otimes` is defined as:

.. math:: f(x)\otimes g(x)\equiv \int_x^1\frac{dy}{y}f(y)g\left(\frac{x}{y}\right)=\int_x^1\frac{dy}{y}f\left(\frac{x}{y}\right)g(y)\,.

The splitting functions :math:`P_{qq}`, :math:`P^V`, :math:`P^+`, and
:math:`P^-` are usually decomposed as follows:

.. math::
   \begin{array}{l}
   \displaystyle P^\pm \equiv P_{qq}^V \pm P_{q\overline{q}}^V\,, \\
   \\
   \displaystyle P_{qq} \equiv P^+ + n_f (P_{qq}^S + P_{q\overline{q}}^S)\,,\\
   \\
   \displaystyle P^V \equiv P^- + n_f (P_{qq}^S - P_{q\overline{q}}^S)\,,
   \end{array}\,,

where :math:`n_f` is the number of active flavours at a given scale
:math:`\mu` and the splitting functions :math:`P_{qq}^V`,
:math:`P_{q\overline{q}}^V`, :math:`P_{qq}^S`,
:math:`P_{q\overline{q}}^S` have the usual perturbative expansion:

.. math:: P(x,\mu)=\sum_{n=0}\left(\frac{\alpha_s(\mu)}{4\pi}\right)^{n+1}P^{(n)}(x)\,.

Given the expansion above, one can show that at
:math:`\mathcal{O}(\alpha_s)`, *i.e.* leading order, all coefficients
but :math:`P_{qq}^{V,(0)}` vanish. It is then easy to see that:

.. math:: P_{qq}^{(0)} = P^{V,(0)} = P ^{+,(0)}= P^{-,(0)}=P_{qq}^{V,(0)}\,.

This means that the evolution equations in
Eq. :eq:`eq:evoleqs` have all the same evolution kernel.

If one wants to include NLO corrections, one finds that the
:math:`\mathcal{O}(\alpha_s^2)` coefficients :math:`P_{qq}^{V,(1)}` and
:math:`P_{q\overline{q}}^{V,(1)}` are different from zero while
:math:`P_{qq}^{S,(1)}` and :math:`P_{q\overline{q}}^{S,(1)}` vanish.
This immediately implies that :math:`P_{qq} = P^+` and
:math:`P^V = P^-`. Therefore, at NLO, the evolution equations are fully
determined by the functions :math:`P_{qq}^{V,(0)}`,
:math:`P_{qq}^{V,(1)}`, and :math:`P_{q\overline{q}}^{V,(1)}`. We are
now in the position to discuss the specific expressions of these
functions for both PDFs and FFs. A further simplification is given by
the fact that the function :math:`P_{qq}^{V,(0)}` is the same for both
PDFs and FFs. However, this is no longer the case for
:math:`P_{qq}^{V,(1)}` and :math:`P_{q\overline{q}}^{V,(1)}` whose form
differs between PDFs and FFs.

In order to carry out the implementation of the splitting functions in
``APFEL``, it is necessary to make sure that the expressions of the
single coefficients of the perturbative expansions have the following
structure:

.. math::
   :label: eq:splittingfuncs
   P(y) = R(y) + S \left(\frac1{1-x}\right)_++ L\delta(1-y)\,,

\ where :math:`R` is a regular function in :math:`x=1`, and :math:`S`
and :math:`L` are numerical coefficients. The plus prescription used
above has the following definition upon integration with a test function
:math:`f`:

.. math:: \int_0^1dy\left(\frac1{1-y}\right)_+f(y)\equiv \int_0^1dy\frac{f(y)-f(1)}{1-y}\,.

An important detail to notice is that the definition above is strictly
true only if the lower integration bound is equal to zero. In actual
facts, this is never the case because Mellin convolutions involving
plus-prescripted functions have the following structure:

.. math:: \int_x^1dy\left(\frac1{1-y}\right)_+f(y)\,.

with :math:`0<x<1`. This integral can be manipulated as follows:

.. math::
   \begin{array}{rcl}
   \displaystyle \int_x^1dy\left(\frac1{1-y}\right)_+f(y) &=&
                                                              \displaystyle
                                                              \int_0^1dy\left(\frac1{1-y}\right)_+f(y)
                                                              -
                                                              \displaystyle
                                                              \int_0^xdy\left(\frac1{1-y}\right)_+f(y)\\
   \\
   &=& 
                                                              \displaystyle
                                                              \int_0^1dy \frac{f(y)-f(1)}{1-y}
                                                              -
                                                              \displaystyle
                                                              \int_0^xdy\frac{f(y)}{1-y}\\
   \\
   &=& 
                                                              \displaystyle
                                                              \int_x^1dy \frac{f(y)-f(1)}{1-y}
                                                              -
                                                              \displaystyle
                                                              f(0)\int_0^x\frac{dy}{1-y}\\
   \\
   &=& 
                                                              \displaystyle
                                                              \int_x^1dy \left(\frac{1}{1-y}\right)_{\oplus}f(y)
                                                              +f(0)\ln(1-x)\\
   \\
   &=& 
                                                              \displaystyle
                                                              \int_x^1dy \left[\left(\frac{1}{1-y}\right)_{\oplus}+\ln(1-x)\delta(1-y)\right]f(y)\,,
   \end{array}

where we have defined a “generalised” plus prescription that holds in
its form regardless of the lower integration bound:

.. math:: \int_x^1dy\left(\frac1{1-y}\right)_\oplus f(y)\equiv \int_x^1dy\frac{f(y)-f(1)}{1-y}\,.

Therefore, a Mellin-like convolution of the splitting function in
Eq. :eq:`eq:splittingfuncs` with the test
function :math:`f` will take the form:

.. math::
   \int_x^1 dy P(y) f(y) = \int_x^1 dy \left[ R(x)
     +S\left(\frac1{1-y}\right)_\oplus + \left(L+S\ln(1-x)\right)\delta(1-y)\right] f(y)\,.

This provides a suitable expression for the implementation in
`` APFEL``. Therefore, one just needs to manipulate the expressions
given in the literature to reduce them to the form of
Eq. :eq:`eq:splittingfuncs`. This is typically an
easy task.

Let us start with :math:`P_{qq}^{V,(0)}` that we take from Eq. (38) of
Ref. (Vogelsang 1998). [1]_ After a simple manipulation, it takes the
form:

.. math::
   :label: eq:LOsplitting
   P_{qq}^{V,(0)}(y) = 2C_F\left[-2+2\left(\frac{1}{1-y}\right)_++\frac32\delta(1-y)\right]\,.

In this form it is easy to identify the elements introduced in
Eq. :eq:`eq:splittingfuncs`. In particular, we
find that :math:`R(x)=-4C_F`, :math:`S=4C_F`, and :math:`L=3C_F`. As
mentioned above, :math:`P_{qq}^{V,(0)}` is the same for PDFs and FFs,
therefore the expression in Eq. :eq:`eq:LOsplitting`
is all one needs to implement the LO evolution of both transversity PDFs
and FFs

We now consider the NLO corrections. In order to distinguish between
PDFs and FF we will use the symbols :math:`\mathcal{P}` and
:math:`\mathbb{P}`, respectively, for the splitting functions. We first
consider the PDF splitting functions that we again take from
Ref. (Vogelsang 1998). We observe that
:math:`\mathcal{P}_{q\overline{q}}^{V,(1)}`, taken from Eq. (44) of this
paper, is a purely regular functions with no plus-prescripted and
:math:`\delta`-function terms. Therefore, it needs no manipulation:

.. math::
   :label: eq:PDFLOsplittingqqb
   \mathcal{P}_{q\overline{q}}^{V,(1)} (y) =
   4C_F\left(C_F-\frac12C_A\right)\left[ - 1 + y -\frac{4S_2(y)}{1+y}\right]\,,

with:

.. math:: S_2(y) = -2\mbox{Li}_2(-y) - 2\ln y\ln(1+y)+\frac12\ln^2y-\frac{\pi^2}{6}\,.

The function :math:`\mathcal{P}_{qq}^{V,(1)}` from Eq. (43) of
Ref. (Vogelsang 1998) is instead more complicated but can be recasted in
the form of Eq. :eq:`eq:LOsplitting` as:

.. math::
   :label: eq:PDFLOsplittingqq
   \begin{array}{rcl} 
     \mathcal{P}_{qq}^{V,(1)} (y) &=& \Bigg\{\displaystyle 4C_F^2 \left[ 1-y -\left( \frac{3}{2} + 
                                      2 \ln (1-y) \right) \frac{2y\ln
                                      y}{1-y}\right]\\
     \\
                                  &+& \displaystyle 2C_F C_A\left[
                                      - \frac{143}{9} +
                                      \frac{2\pi^2}{3} +y + \left( \frac{11}{3} 
                                      + \ln y  \right) \frac{2y\ln
                                      y}{1-y}\right] + \displaystyle \frac{8}{3} n_f C_F
                                      T_R  \left[ - \frac{2y\ln y}{1-y}+ \frac{10}{3} \right]\Bigg\}\\
     \\
                                  &+& \displaystyle \Bigg\{2C_F C_A\left(\frac{134}{9} -
                                      \frac{2\pi^2}{3} \right) - \frac{80}{9} n_f C_F
                                      T_R  \Bigg\}\left(\frac{1}{1-y}\right)_+\\
     \\
                                  &+& \displaystyle\Bigg\{ 4C_F^2 \left( \frac{3}{8} -\frac{\pi^2}{2} + 6\zeta (3) 
                                      \right) + 2C_FC_A \left( \frac{17}{12} + \frac{11 \pi^2}{9} -
                                      6 \zeta (3) \right) - \frac{8}{3} n_f C_F
                                      T_R\left( \frac{1}{4} + \frac{\pi^2}{3} \right)\Bigg\} \delta(1-y)\,,
   \end{array}

\ where the regular, plus-prescripted, and :math:`\delta`-function terms
are enclosed between curly brackets.

We can now turn to consider the splitting functions for FFs. In this
case we take the expressions from Ref. (Stratmann and Vogelsang 2002).
From Eq. (17) of this paper we immediately see that:

.. math:: \mathbb{P}_{q\overline{q}}^{V,(1)} (y)=\mathcal{P}_{q\overline{q}}^{V,(1)} (y)\,,

where :math:`\mathcal{P}_{q\overline{q}}^{V,(1)}` is given in
Eq. :eq:`eq:PDFLOsplittingqqb`. For
:math:`\mathbb{P}_{qq}^{V,(1)}` we instead find:

.. math::
   \begin{array}{rcl} 
     \mathbb{P}_{qq}^{V,(1)} (y) &=& \Bigg\{\displaystyle 4C_F^2 \left[ 1-y +\left( \frac{3}{2} + 
                                      2 \ln (1-y) -2\ln y\right) \frac{2y\ln
                                      y}{1-y}\right]\\
     \\
                                  &+& \displaystyle 2C_F C_A\left[
                                      - \frac{143}{9} +
                                      \frac{2\pi^2}{3} +y + \left( \frac{11}{3} 
                                      + \ln y  \right) \frac{2y\ln
                                      y}{1-y}\right] + \displaystyle \frac{8}{3} n_f C_F
                                      T_R  \left[ - \frac{2y\ln y}{1-y}+ \frac{10}{3} \right]\Bigg\}\\
     \\
                                  &+& \displaystyle \Bigg\{2C_F C_A\left(\frac{134}{9} -
                                      \frac{2\pi^2}{3} \right) - \frac{80}{9} n_f C_F
                                      T_R  \Bigg\}\left(\frac{1}{1-y}\right)_+\\
     \\
                                  &+& \displaystyle\Bigg\{ 4C_F^2 \left( \frac{3}{8} -\frac{\pi^2}{2} + 6\zeta (3) 
                                      \right) + 2C_FC_A \left( \frac{17}{12} + \frac{11 \pi^2}{9} -
                                      6 \zeta (3) \right) - \frac{8}{3} n_f C_F
                                      T_R\left( \frac{1}{4} + \frac{\pi^2}{3} \right)\Bigg\} \delta(1-y)\,,
   \end{array}

that is just a small difference in the regular term as compared to
:math:`\mathcal{P}_{qq}^{V,(1)}` in
Eq. :eq:`eq:PDFLOsplittingqq`.

**References**

**References**

.. container:: references csl-bib-body hanging-indent
   :name: refs

   .. container:: csl-entry
      :name: ref-Stratmann:2001pt

      Stratmann, Marco, and Werner Vogelsang. 2002. “Next-to-leading
      order QCD evolution of transversity fragmentation functions.”
      *Phys. Rev. D* 65: 057502.
      https://doi.org/10.1103/PhysRevD.65.057502.

   .. container:: csl-entry
      :name: ref-Vogelsang:1997ak

      Vogelsang, Werner. 1998. “Next-to-leading order evolution of
      transversity distributions and Soffer’s inequality.” *Phys. Rev.
      D* 57: 1886–94. https://doi.org/10.1103/PhysRevD.57.1886.

.. [1]
   A factor 2 is introduced to account for the different expansion
   parameter, here :math:`\alpha_s/4\pi` rather than
   :math:`\alpha_s/2\pi`
