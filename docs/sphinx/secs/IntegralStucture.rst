=====================
Convolution integrals
=====================

.. contents::
   :depth: 3
..

This section discusses the approach implemented in ``APFEL++`` to
compute Mellin-convolution integrals. This kind of convolutions appears
very often in the presence of collinearly factorised quantities as
integrals over the longitudinal momentum fraction of partons. This is
for example the case of all cross sections in collinear factorisation
with hadrons in the initial and/or final state. As a consequence,
Mellin-convolution integrals also appear in the DGLAP evolution
equations. Also in transverse-momentum-dependent (TMD) factorisation
they play an important role when matching TMD distributions onto the
respective collinear ones. It is therefore of central importance to put
in place an efficient and accurate strategy to numerically compute these
convolutions. In doing so, it is also important to take into account as
much as possible all possible flavours in which they may appear. For
example, different variants apply in the presence of heavy-quark masses
or when considering “non-forward” integrals in the framework of
exclusive physics.

Structure of the integrals
==========================

The general structure of the integrals considered in ``APFEL++`` has the
form of a Mellin convolution between an *operator* :math:`O`, this is a
potentially complicated known function such as a splitting function or a
partonic cross section, [1]_ and a *distribution* function :math:`d`, as
for example a parton distribution function (PDF) or a fragmentation
function (FF). The explicit form of these integrals reads:

.. math::
   I(x) = x\int_0^1dz\int_0^1dy\, O(z)d(y)\delta(x-yz) =
   x\int_x^1\frac{dz}{z} O \left(\frac{x}{z}\right) d(z) =
   x\int_x^1\frac{dz}{z} O(z)d\left(\frac{x}{z}\right)\,.
   :label: eq:ZMconv

\ Typically in the presence of mass effects, the integration phase space
may be modified and the convolution in
Eq. :eq:`eq:ZMconv` is generalised as:

.. math::
   I(x,\eta) = x\int_{x/\eta}^1\frac{dz}{z} O(z,\eta)d\left(\frac{x}{\eta
       z}\right)\,,
   :label: eq:Genconv

where :math:`\eta\leq1`, with :math:`\eta = 1` reproducing
Eq. :eq:`eq:ZMconv`. However, it is convenient to rewrite
the integral in Eq. :eq:`eq:Genconv` in a form in which
the lower bound of the integral is equal to :math:`x`. This is easily
done by performing the change of variable :math:`y =\eta z`, so that:

.. math::
   I(x,\eta)= \int_x^\eta dy\,
     O\left(\frac{y}{\eta},\eta\right)\,\frac{x}{y}d\left(\frac{x}{y}\right)\,.
   :label: eq:Genconv1

In order to enable the precomputation the expensive part of the integral
in Eq. :eq:`eq:Genconv1`, we use the standard
interpolation formula on the distribution :math:`d`:

.. math::
   \frac{x}{y} d\left(\frac{x}{y}\right) = \sum_{\alpha=0}^{N_x} x_\alpha
   d(x_\alpha) w_{\alpha}^{(k)}\left(\frac{x}{y}\right)\,,
   :label: eq:Interpolation

\ where :math:`\alpha` runs over the nodes of a grid in :math:`x` and
the weights :math:`w_{\alpha}` are assumed to be Langrange polynomials
of degree :math:`k`. Now we plug
Eq. :eq:`eq:Interpolation` into
Eq. :eq:`eq:Genconv1` and specialise the computation of
the integral :math:`I` to the girid point :math:`x_\beta`. This gives:

.. math::
   I(x_\beta,\eta)=  \sum_{\alpha=0}^{N_x} \overline{d}_\alpha\int_{x_\beta}^\eta dy\,
     O\left(\frac{y}{\eta},\eta\right)\,w_{\alpha}^{(k)}\left(\frac{x_\beta}{y}\right)\,.

where we have defined
:math:`\overline{d}_\alpha = x_\alpha d(x_\alpha)`. As shown in the
section devoted to the interpolation procedure, under some specific
conditions the interpolating functions are such that:

.. math:: w_{\alpha}^{(k)}\left(\frac{x_\beta}{y}\right)\neq 0\quad\mbox{for}\quad c < y < d\,,

with:

.. math::
   c =
     \mbox{max}(x_\beta,x_\beta/x_{\alpha+1}) \quad\mbox{and}\quad d =
     \mbox{min}(\eta,x_\beta/x_{\alpha-k}) \,.
   :label: eq:intlims

Therefore, Eq. :eq:`eq:Genconv1` can be adjusted as:

.. math::
   I(x_\beta,\eta)=  \sum_{\alpha=0}^{N_x} \overline{d}_\alpha\int_c^d dy\,
     O\left(\frac{y}{\eta},\eta\right)\,w_{\alpha}^{(k)}\left(\frac{x_\beta}{y}\right)\,.
   :label: eq:Genconv2

Finally, we change the integration variable back to :math:`z=y/\eta` so
that Eq. :eq:`eq:Genconv2` becomes:

.. math::
   I(x_\beta,\eta)=  \sum_{\alpha=0}^{N_x} \overline{d}_\alpha \underbrace{\left[\eta\int_{c/\eta}^{d/\eta} dz\,
     O\left(z,\eta\right)\,w_{\alpha}^{(k)}\left(\frac{x_\beta}{\eta
         z}\right)\right]}_{\Gamma_{\alpha\beta}}\,.
   :label: eq:Genconv3

To summarise, Eq. :eq:`eq:Genconv3` allows one to
compute the integral :math:`I` on the nodes of a given grid as a
weighted sum of the values of distribution :math:`d` on the nodes
themselves. The weights are given by the appropriate integral of the
operator :math:`O` and the interpolating functions
:math:`w_{\alpha}^{(k)}`. The value of :math:`I` for a generic value of
:math:`x` can finally be obtained by interpolation using a formula
similar to Eq. :eq:`eq:Interpolation`.
Importantly, Eq. :eq:`eq:Genconv3` provides a fast way
of computing convolution integrals. Indeed, the weights
:math:`\Gamma_{\alpha\beta}` can be precomputed and stored once and for
all and the integral :math:`I` can be quickly computed with different
distributions :math:`d` by simply taking the weighted average.

We now need an operational way to compute the weights
:math:`\Gamma_{\alpha\beta}`. To do so, it is necessary to know the
general structure of the operator :math:`O`. Quite generally, the
operator :math:`O` splits as follows:

.. math::
   O(z,\eta) = R(z,\eta)+\sum_{i}\left[P^{(i)}(z,\eta)\right]_+
     S^{(i)}(z,\eta)+ L(\eta)\delta(1-z)\,,
   :label: eq:CoeffFuncs

\ where :math:`R` and :math:`S^{(i)}` are a regular functions at
:math:`z=1`, that is:

.. math:: R(1,\eta) =  \lim_{z\rightarrow 1} R(z,\eta) = K(\eta)\quad\mbox{and}\quad S^{(i)}(1,\eta) =  \lim_{z\rightarrow 1} S^{(i)}(z,\eta) = J^{(i)}(\eta)\,,

being :math:`K`, :math:`J^{(i)}`, and :math:`L` a finite function of
:math:`\eta`. The functions :math:`P^{(i)}` are instead singular at
:math:`z=1` and non-integrable in the limit :math:`\eta\rightarrow1` and
thus are regularised through the :math:`+`-prescription. Plugging
Eq. :eq:`eq:CoeffFuncs` into the definition of
:math:`\Gamma_{\alpha\beta}` in Eq. :eq:`eq:Genconv3`
and making use of the definition of :math:`+`-prescription gives:

.. math::
   \begin{array}{rcl}
   \Gamma_{\beta\alpha} &=& \eta\displaystyle \int_{c/\eta}^{d/\eta} dz\left\{R\left(z,\eta\right)w_{\alpha}^{(k)}\left(\frac{x_\beta}{\eta z}\right)+\sum_{i}P^{(i)}(z,\eta)\left[S^{(i)}\left(z,\eta\right)w_{\alpha}^{(k)}\left(\frac{x_\beta}{\eta z}\right)-S^{(i)}(1,\eta)w_{\alpha}^{(k)}\left(\frac{x_\beta}{\eta}\right)\right]\right\}\\
   \\
   &+&\displaystyle \eta\left[L(\eta)-\sum_{i}S^{(i)}(1,\eta)\int_0^{c/\eta}dz
       P^{(i)}(z,\eta)\right]w_{\alpha}^{(k)}\left(\frac{x_\beta}{\eta}\right)\,,
   \end{array}
   :label: eq:FinalExpression

that can be further manipulated changing the integration variable into
:math:`y=\eta z`:

.. math::
   \begin{array}{rcl}
     \Gamma_{\beta\alpha} &=& \displaystyle \int_{c}^{d} dy\left\{R\left(\frac{y}{\eta},\eta\right)w_{\alpha}^{(k)}\left(\frac{x_\beta}{y}\right)+\sum_{i}P^{(i)}\left(\frac{y}{\eta},\eta\right)\left[S^{(i)}\left(\frac{y}{\eta},\eta\right)w_{\alpha}^{(k)}\left(\frac{x_\beta}{y}\right)-S^{(i)}(1,\eta)w_{\alpha}^{(k)}\left(\frac{x_\beta}{\eta}\right)\right]\right\}\\
     \\
                          &+&\displaystyle \eta \left[L(\eta)+\sum_{i}S^{(i)}(1,\eta)Q^{(i)}\left(\frac{c}{\eta},\eta\right)\right]w_{\alpha}^{(k)}\left(\frac{x_\beta}{\eta}\right)\,.
   \end{array}
   :label: eq:FinalExpression2

where we have defined:

.. math:: Q^{(i)}(a,\eta)\equiv-\int_0^{a}dz P^{(i)}\left(z,\eta\right).

These integrals can, most of the times, be computed analytically.

Eqs. :eq:`eq:FinalExpression`-:eq:`eq:FinalExpression2`
expose the full potential complexity of the task of computing the
weights :math:`\Gamma_{\alpha\beta}`. However, in many practical
applications some simplifications apply. For example, in the case of the
perturbative coefficients of the (:math:`\overline{\mbox{MS}}`)
splitting functions there are two simplications: the first is that
:math:`\eta=1`, the second is that there is one single term in the sum
over :math:`i` (:math:`i=0`) such that and the general form of the
function :math:`P` is:

.. math:: P^{(0)}(z,\eta)\rightarrow \frac1{1-z}\,,

\ so that:

.. math:: Q^{(0)}(a,\eta)=\ln(1-a)\,.

Considering that:

.. math:: w_\alpha^{(k)}(x_\beta) = \delta_{\alpha\beta}\,,

and that the expressions can always be manipulated in such a way that
the coefficient of the :math:`+`-prescripted term :math:`S` is a
constant, one finds that:

.. math::
   \Gamma_{\beta\alpha} = \int_{c}^{d}
   dz\left\{R\left(z\right)w_{\alpha}^{(k)}\left(\frac{x_\beta}{z}\right)+\frac{S}{1-z}\left[w_{\alpha}^{(k)}\left(\frac{x_\beta}{z}\right)-\delta_{\alpha\beta}\right]\right\}+\displaystyle
   \left[S\ln\left(1-c\right)+L \right]\delta_{\alpha\beta}\,.
   :label: eq:SplittingFunctions

The same kind of simplifications apply to the case of the perturbative
Zero-Mass (ZM) coefficient functions with the only exception that the
sum over :math:`i` extends to more terms depending on the perturbative
order. For example at :math:`\mathcal{O}(\alpha_s)`, :math:`i.e.` at
NLO, the general form of the coefficient functions reads:

.. math::
   \begin{array}{rcl}
     \Gamma_{\beta\alpha} &=&\displaystyle \int_{c}^{d}
                              dz\left\{R\left(z\right)w_{\alpha}^{(k)}\left(\frac{x_\beta}{z}\right)
     + \frac{S^{(0)}+S^{(1)}\ln(1-z)}{1-z}\left[w_{\alpha}^{(k)}\left(\frac{x_\beta}{z}\right)-\delta_{\alpha\beta}\right]\right\}\\
   \\
   &+&\displaystyle
       \left[S^{(0)}\ln\left(1-c\right)+\frac12
       S^{(1)}\ln^2\left(1-c\right)+L \right]\delta_{\alpha\beta}\,.
   \end{array}
   :label: eq:ZMCoeffFunctions

Thing get more complicated when considering massive coefficient
functions. In this case, the role of :math:`\eta` is played by the
ration :math:`m/Q`, where :math:`m` is the mass of the quark and
:math:`Q` the hard scale. In addition, the :math:`+`-prescripted
contributions can become more convoluted.

As an aside, it is interesting to notice that we have silently decided
to use one of the two possible ways to compute the convolution integral
:math:`I(x)` in Eq. :eq:`eq:ZMconv`. Specifically, we
have chosen the rightmost equality. However, the equality:

.. math:: I(x) = x\int_x^1\frac{dz}{z} O \left(\frac{x}{z}\right) d(z)\,,

\ is equally valid. If one choses this particular version of the
convolution integral and follows the interpolation procedure discussed
above, the result for :math:`I` computed in :math:`x_\beta` would be
identical to Eq. :eq:`eq:Genconv3` (with
:math:`\eta=1`) except for the fact that the weights
:math:`\Gamma_{\beta \alpha}` are computed as:

.. math::
   \Gamma_{\beta\alpha} = \int_{a}^{b}dz\,
   O\left(\frac{x_\beta}{z}\right)w_{\alpha}^{(k)}(z)
   :label: eq:oldway

with:

.. math:: a \equiv \mbox{max}(x_\beta,x_{\alpha-k})\quad\mbox{and}\quad b \equiv \mbox{min}(1,x_{\alpha+1})\,.

This immediately implies that the integral in
Eq. :eq:`eq:Genconv3` and that in
Eq. :eq:`eq:oldway` have to coincide. This is indeed the
case and the only reason to choose
Eq. :eq:`eq:Genconv3` rather than
Eq. :eq:`eq:oldway` is convenience. Specifically,
Eq. :eq:`eq:oldway` implies computing the operator
:math:`O` at :math:`x_\beta/z` rather than at :math:`z` and this makes
the treatment of :math:`\delta`-functions and :math:`+`-distributions
possibly present in :math:`O` a little more delicate. However,
Eq. :eq:`eq:oldway` remains a valid alternative that we
just decided not to pursue.

.. container::
   :name: double-convolutions

   .. rubric:: Double convolutions
      :name: double-convolutions

The technique discussed above to compute convolution integrals as a
weighted sums of a given input distributions :math:`d` can be, under
certain circumstances, generalised to convolution integrals involving
*two* distributions. Relevant examples are the Drell-Yan (DY) and the
semi-inclusive deep-inelastic scattering (SIDIS) cross sections
integrated over the transverse momentum of the virtual boson. In these
cases the cross section is proportional to a double convolutions between
partonic cross sections and a pair of non-perturbative distributions.
The general form of this double convolution is:

.. math::
   J(x_1,x_2) =
   x_1 x_2\int_{x_1}^1\frac{dy_1}{y_1}\int_{x_2}^1\frac{dy_2}{y_2}
   O(y_1,y_2) d^{(1)}\left(\frac{x_1}{y_1}\right)
     d^{(2)}\left(\frac{x_2}{y_2}\right)\,.

\ Applying the same interpolation procedure as in the single-convolution
case gives:

.. math::
   \begin{array}{rcl}
   J(x_\delta,x_\gamma) &=&\displaystyle
   \int_{x_\delta}^1dy_1\int_{x_\gamma}^1dy_2\,
   O(y_1,y_2) \left[\frac{x_\delta}{y_1}d^{(1)}\left(\frac{x_\delta}{y_1}\right)\right]
     \left[\frac{x_\gamma}{y_2}d^{(2)}\left(\frac{x_\gamma}{y_2}\right)\right]\\
   \\
   &=&\displaystyle
       \sum_{\alpha=0}^{N_x}\sum_{\beta=0}^{N_x} \overline{d}^{(1)}_\beta
       \overline{d}^{(2)}_\alpha \underbrace{\left[\int_{x_\delta}^1 dy_1 \int_{x_\gamma}^1 dy_2\,
       O(y_1,y_2)\, w_{\beta}^{(k)}\left(\frac{x_\delta}{y_1}\right) 
       w_{\alpha}^{(k)}\left(\frac{x_\gamma}{y_2}\right)\right]}_{\Theta^{\beta\delta,\alpha\gamma}}\,.
   \end{array}
   :label: eq:DoubleZMconv

In Eq. :eq:`eq:DoubleZMconv` we have assumes that
there are no mass corrections and thus the convolutions take the
simplest form. In the case of double convolutions, the operator
:math:`O` is a function of two variables, :math:`y_1` and :math:`y_2`,
and, as in the case of the single convolutions, it receives three kinds
of contribution in both these variables: local terms proportional to
:math:`\delta`-functions, singular terms proportional to
:math:`+`-prescripted functions, and regular terms. The complication
here is that these contributions from the two variables :math:`y_1` and
:math:`y_2` mix and thus, for example, terms local in :math:`y_1` and
singular in :math:`y_2` may also appear. It is thus necessary to
identify the general structure of the function :math:`O` to see whether
it is possible to decompose the double operator
:math:`\Theta^{\beta\delta,\alpha\gamma}` into products of single
operators.

As we will explicitly show below, in the case of SIDIS up to
:math:`\mathcal{O}(\alpha_s)` (NLO) the general structure of the
function :math:`O` can be inferred looking at Eqs. (C.2)-(C.7) of
Ref. (Florian, Stratmann, and Vogelsang 1998):

.. math::
   \begin{array}{rclclcl}
     O(y_1,y_2) &=&\displaystyle {\rm LL}\,\delta(1-y_1)\delta(1-y_2) &+&\displaystyle {\rm LS}\,
                                                                          \delta(1-y_1)\left[\frac{\ln(1-y_2)}{1-y_2}\right]_+ &+&\displaystyle \delta(1-y_1)\,{\rm LR}(y_2) \\
     \\
                &+&\displaystyle {\rm SL}\,\left[\frac{\ln(1-y_1)}{1-y_1}\right]_+\delta(1-y_2)
                                                                      &+&\displaystyle {\rm SS}\,\left[\frac1{1-y_1}\right]_+
                                                                          \left[\frac1{1-y_2}\right]_+ &+&\displaystyle \left[\frac1{1-y_1}\right]_+\,{\rm
                                                                                                           SR}(y_2)\\
     \\
                &+&\displaystyle {\rm RL}(y_1)\,\delta(1-y_2)   &+&\displaystyle
                                                                    {\rm RS}(y_1)
                                                                    \,\left[\frac1{1-y_2}\right]_+&+&\displaystyle
                                                                                                      \sum_iK_i{\rm
                                                                                                      R}_i^{(1)}(y_1){\rm
                                                                                                      R}_i^{(2)}(y_2)\,.
   \end{array}
   :label: eq:DoubleFuncStruct

\ It is clear that in
Eq. :eq:`eq:DoubleFuncStruct` all terms
factorise into a part that only depends on :math:`y_1` and a part that
only depends on :math:`y_2`. This is the crucial feature that enables
use the technology developed above for the single convolutions. Plugging
Eq. :eq:`eq:DoubleFuncStruct` into
Eq. :eq:`eq:DoubleZMconv`, one finds that:

.. math::
   \begin{array}{rclclcl}
     \Theta^{\beta\delta,\alpha\gamma} &=& {\rm LL} \,
                                           \Gamma^{\rm L}_{\beta\delta}\Gamma^{\rm L}_{\alpha\gamma}
     &+& \displaystyle {\rm LS}\,\Gamma^{\rm L}_{\beta\delta}\Gamma_{\alpha\gamma}^{\rm S1}
     &+&\displaystyle \Gamma^{\rm L}_{\beta\delta} \Gamma_{\alpha\gamma}^{\rm LR}\\
     \\
                                       &+&{\rm SL}\,\displaystyle\Gamma_{\beta\delta}^{\rm S1}\Gamma^{\rm L}_{\alpha\gamma} 
     &+&
        {\rm SS}\,\Gamma_{\beta\delta}^{\rm
         S0}\Gamma_{\alpha\gamma}^{\rm S0}&+& \Gamma_{\beta\delta}^{\rm S0}\Gamma_{\alpha\gamma}^{\rm SR}\\
   \\
   &+& \displaystyle \Gamma_{\beta\delta}^{\rm RL}\Gamma^{\rm L}_{\alpha\gamma} &+& \Gamma_{\beta\delta}^{\rm RS}\Gamma_{\alpha\gamma}^{\rm
         S0} &+& \displaystyle \sum_i K_i \Gamma_{\beta\delta}^{R_i^{(1)}}\Gamma_{\alpha\gamma}^{R_i^{(2)}}
   \end{array}
   :label: eq:MasterFormula

with:

.. math::
   \begin{array}{rcl}
     \Gamma_{\alpha\beta}^{\rm L} &=& \displaystyle  \int_{c_{\alpha\beta}}^{d_{\alpha\beta}}
                                          dz\,\delta(1-z)w_{\beta}^{(k)}\left(\frac{x_\alpha}{z}\right)
                                      =
                                      w_{\beta}^{(k)}\left(x_\alpha\right)
                                      =  \delta_{\alpha\beta}\\
   \\
     \Gamma_{\alpha\beta}^{{\rm S}n} &=& \displaystyle  \int_{c_{\alpha\beta}}^{d_{\alpha\beta}}
                                          dz\frac{\ln^n(1-z)}{1-z}\left[w_{\beta}^{(k)}\left(\frac{x_\alpha}{z}\right)-\delta_{\alpha\beta}\right]+\frac1{(n+1)!}\ln^{n+1}\left(1-c_{\alpha\beta}\right)\delta_{\alpha\beta}\\
     \\
     \Gamma_{\alpha\beta}^{f} &=&\displaystyle \int_{c_{\alpha\beta}}^{d_{\alpha\beta}}
                                        dz\,f(z)\,w_{\beta}^{(k)}\left(\frac{x_\alpha}{z}\right)
   \end{array}
   :label: eq:PrecompOp

where :math:`f` is a regular function and the integration bounds are
defined as:

.. math::
   c_{\alpha\beta} =
     \mbox{max}(x_\alpha,x_\alpha/x_{\beta+1}) \quad\mbox{and}\quad d_{\alpha\beta} =
     \mbox{min}(1,x_\alpha/x_{\beta-k}) \,.

In general terms, we assume that it is always be possible to write an
object of the kind of :math:`\Theta^{\beta\delta,\alpha\gamma}` as
series of bilinear terms:

.. math::
   \Theta^{\beta\delta,\alpha\gamma} = \sum_j
   C_j\Gamma_j^{(1),\beta\delta}\Gamma_j^{(2),\alpha\gamma}
   :label: eq:OpSeries

where :math:`C_j` are scalar coefficients, and the weights
:math:`\Gamma_j^{(1),\beta\delta}` and
:math:`\Gamma_j^{(2),\alpha\gamma}` can be computed using the technology
discussed in the previous section. Plugging
Eq. :eq:`eq:OpSeries` into
Eq. :eq:`eq:DoubleZMconv`, one finds that:

.. math::
   J(x_\delta,x_\gamma) =
     \sum_j C_j  f_j^{(1),\delta} f_j^{(2),\gamma}\,,
   :label: eq:DoubleZMconv2

where we have defined:

.. math:: f_j^{(1),\delta}\equiv \sum_{\beta=0}^{N_x} \overline{d}^{(1)}_\beta\Gamma_j^{(1),\beta\delta}\quad\mbox{and}\quad f_j^{(2),\gamma}\equiv \sum_{\alpha=0}^{N_x}\overline{d}^{(2)}_\alpha \Gamma_j^{(2),\alpha\gamma}\,.

Eq. :eq:`eq:DoubleZMconv2` shows that, under the
hypothesis that the operator :math:`O(y_1,y_2)` can be expressed as a
series of terms whose dependence on :math:`y_1` and :math:`y_2`
factorizes, [2]_ the double convolution in
Eq. :eq:`eq:DoubleZMconv` is given by a series of
bilinear distributions (:math:`f_j^{(1),\delta}` and
:math:`f_j^{(2),\gamma}`) singularly obtained as convolutions of single
operators with the input distributions :math:`d^{(1)}` and
:math:`d^{(2)}`. This is a particularly useful achievement that allows
us to compute double convolutions without the need of extending the
integration and the interpolation procedures to two dimensions with an
obvious gain in accuracy and performance. As a matter of fact, the same
argument can be extended to a multiple convolution of the function
:math:`O(\{y_i\})`, which again can be expressed as a series of
:math:`n`-linear terms, with :math:`i=1,\dots,n`, with :math:`n`
distributions:

.. math:: J(\{x_{\alpha_i}\}) = \sum_j C_j  \prod_{i=1}^nf_j^{(i),\alpha_i}\,,

with:

.. math:: f_j^{(i),\alpha_i} \equiv \sum_{\beta=0}^{N_x} \overline{d}^{(i)}_\beta\Gamma_j^{(i),\beta\alpha_i}\,.

This technology could be useful for more complicated observables, like
cross sections in :math:`pp` collisions with an identified hadron in the
final state, that requires for example three convolutions.

The challenging part of the procedure just presented resides in the
“pre-processing” of the function :math:`O(y_1,y_2)` that has to be
analytically manipulated to disentangle the single terms. This step,
however, has to be taken only once.

Before going into a concrete application, it is useful to connect
Eq. :eq:`eq:MasterFormula` to
Eq. :eq:`eq:OpSeries` by identifying number and form of
the coefficients and weights involved. Specifically, assuming that the
series in the last term in the r.h.s. of
Eq. :eq:`eq:MasterFormula` has :math:`r` terms,
the series in Eq. :eq:`eq:OpSeries` will have
:math:`8+r` terms, that is:

.. math::
   \Theta^{\beta\delta,\alpha\gamma} = \sum_{j=1}^{8+r}
   C_j\Gamma_j^{(1),\beta\delta}\Gamma_j^{(2),\alpha\gamma}
   :label: eq:OpSeriesCon

\ with:

.. math::
   \begin{array}{lclll}
   j = 1 &:& C_1 = {\rm LL}, & \Gamma_1^{(1),\beta\delta} = \Gamma^{\rm L}_{\beta\delta}, & \Gamma_1^{(2),\alpha\gamma} = \Gamma^{\rm L}_{\alpha\gamma},\\
   j = 2 &:& C_2 = {\rm LS}, & \Gamma_2^{(1),\beta\delta} = \Gamma^{\rm L}_{\beta\delta}, & \Gamma_2^{(2),\alpha\gamma} = \Gamma_{\alpha\gamma}^{\rm S1},\\
   j = 3 &:& C_3 = 1, & \Gamma_3^{(1),\beta\delta} = \Gamma^{\rm L}_{\beta\delta}, & \Gamma_3^{(2),\alpha\gamma} = \Gamma_{\alpha\gamma}^{\rm LR},\\
   j = 4 &:& C_4 = {\rm SL}, & \Gamma_4^{(1),\beta\delta} = \Gamma_{\beta\delta}^{\rm S1}, & \Gamma_4^{(2),\alpha\gamma} = \Gamma^{\rm L}_{\alpha\gamma},\\
   j = 5 &:& C_5 = {\rm SS},& \Gamma_5^{(1),\beta\delta} = \Gamma_{\beta\delta}^{\rm S0}, & \Gamma_5^{(2),\alpha\gamma} = \Gamma_{\alpha\gamma}^{\rm S0},\\
   j = 6 &:& C_6 = 1, & \Gamma_6^{(1),\beta\delta} = \Gamma_{\beta\delta}^{\rm S0}, & \Gamma_6^{(2),\alpha\gamma} = \Gamma_{\alpha\gamma}^{\rm SR},\\
   j = 7 &:& C_7 = 1, & \Gamma_7^{(1),\beta\delta} = \Gamma_{\beta\delta}^{\rm RL}, & \Gamma_7^{(2),\alpha\gamma} = \Gamma^{\rm L}_{\alpha\gamma},\\
   j = 8 &:& C_8 = 1, & \Gamma_8^{(1),\beta\delta} = \Gamma_{\beta\delta}^{\rm RS}, & \Gamma_8^{(2),\alpha\gamma} = \Gamma_{\alpha\gamma}^{\rm S0},\\
   j = 9 &:& C_9 = K_1 & \Gamma_9^{(1),\beta\delta} = \Gamma_{\beta\delta}^{R_1^{(1)}}, & \Gamma_9^{(2),\alpha\gamma} = \Gamma_{\alpha\gamma}^{R_1^{(2)}},\\
   \vdots & & & & \\
   j = 8 + r &:& C_{8+r} = K_r & \Gamma_{8+r}^{(1),\beta\delta} = \Gamma_{\beta\delta}^{R_r^{(1)}}, & \Gamma_{8+r}^{(2),\alpha\gamma} = \Gamma_{\alpha\gamma}^{R_r^{(2)}}
   \end{array}

It should be noted that, despite the large number of terms in the series
in Eq. :eq:`eq:OpSeriesCon`, the number of weights
to be precomputed is usually pretty limited. In addition, in many cases
many of the terms of the series are zero so that the number of
contributions is further reduced. We can now apply this procedure up to
NLO in QCD to two specific cases: SIDIS first and DY second
(incomplete).

.. container::
   :name: semi-inclusive-deep-inelastic-scattering-sidis

   .. rubric:: Semi-inclusive deep inelastic scattering (SIDIS)
      :name: semi-inclusive-deep-inelastic-scattering-sidis

the structure of the (:math:`p_T`-integrated) SIDIS cross section and
the expressions of the respective hard coefficient functions can be
found in Ref. (Florian, Stratmann, and Vogelsang 1998). Following this
paper, the SIDIS differential cross section for the exchange of a
virtual photon can be written as:

.. math::
   \frac{d^3\sigma}{dxdydz} = 
     \frac{2\, \pi\alpha^2}{xyQ^2} 
     \left[ (1+(1-y)^2) 2xF_1(x,z,Q^2) + 
       2 (1-y) F_L(x,z,Q^2) \right]\,,
   :label: eq:sidis

\ with :math:`Q^2 = - q^2` the (negative) virtuality of the exchanged
photon, :math:`x` and :math:`z` the momentum fractions that enter the
PDFs and the FFs, and :math:`y = Q^2/xs` the inelasticity given in terms
of :math:`Q`, :math:`x`, and the collision energy in the center of mass
:math:`s`. Notice that, as compared to Ref. (Florian, Stratmann, and
Vogelsang 1998), we have absorbed a factor :math:`x` into the definition
of :math:`F_L` as customary for the longitudinal structure function in
inclusive DIS.

We now use the Callan-Gross relation:

.. math::
   F_2 = 2xF_1 + F_L
   :label: eq:CallanGross

\ to replace :math:`2xF_1` with :math:`F_2` in
Eq. :eq:`eq:sidis`:

.. math::
   \frac{d^3\sigma}{dxdydz} = 
     \frac{2\, \pi\alpha^2}{xyQ^2} 
     \left[ Y_+ F_2(x,z,Q^2)
       -y^2 F_L(x,z,Q^2) \right]\,,
   :label: eq:sidis2

where we have defined:

.. math:: Y_+ = 1+(1-y)^2\,.

It is also useful to write Eq. :eq:`eq:sidis2` as
differential in :math:`x`, :math:`Q^2`, and :math:`z`:

.. math::
   \frac{d^3\sigma}{dx\, dQ^2\, dz} = 
     \frac{2\, \pi\alpha^2}{xQ^4} 
     \left[ Y_+ F_2(x,z,Q)
       -y^2 F_L(x,z,Q) \right]\,.
   :label: eq:sidis3

The structure functions :math:`F_2` and :math:`F_L` are given at NLO by:

.. math::
   \begin{array}{rcl}
   \displaystyle  F_{2,L}(x,z,Q) &=& \displaystyle x\sum_{q,\overline{q}} e_q^2 \bigg[ q(x,Q)
       \otimes  C^{2,L}_{qq}(x,z) \otimes D_q(z,Q)  \\
   &+&\displaystyle q(x,Q)  \otimes  C^{2,L}_{gq}(x,z) \otimes D_g(z,Q)+  g(x,Q)  \otimes
     C^{2,L}_{qg}(x,z) \otimes D_q(z,Q) \bigg]\,,
   \end{array}
   :label: eq:f1sidis

where :math:`\{q,g\}` are the quark and gluon PDFs and
:math:`\{D_q,D_g\}` are the quark and gluon FFs, :math:`e_q` is the
electric charge of the quark :math:`q` and
:math:`\{C^{2,L}_{qq},C^{2,L}_{qg},C^{2,L}_{gq}\}` are the relevant
partonic cross sections. The partonic cross sections allow for a
perturbative expansion in power of :math:`\alpha_s`:

.. math:: C = \sum_{n=0} \left(\frac{\alpha_s}{4\pi}\right)^nC^{(n)}

that we truncate to NLO, :math:`i.e.` to :math:`n=1`. At LO
(:math:`n=0`) we have the simple expression:

.. math::
   \begin{array}{rcl}
   C^{2,(0)}_{qq}(x,z) &=& \delta(1-x)\delta(1-z)\,,\\
   \\
   C^{2,(0)}_{qg}(x,z) &=&C^{2,(0)}_{gq}(x,z) = 0\,.
   \end{array}

At NLO (:math:`n=1`) we take the expressions from Appendix C of
Ref. (Florian, Stratmann, and Vogelsang 1998) being careful to take into
account an additional factor two due to the difference in the expansion
parameter (:math:`\alpha_s/4\pi` rather than :math:`\alpha_s/2\pi`). We
also need to combine the expressions for :math:`F_1` and :math:`F_L`
using Eq. :eq:`eq:CallanGross` to obtain the
partonic cross sections for :math:`F_2`. We start with the partonic
cross sections for :math:`F_L` that read:

.. math::
   \begin{array}{rcl}
   C_{qq}^{L,(1)} &=& 8 C_F x z\,, \\
   C_{gq}^{L,(1)} &=& 8 C_F x (1-z)\,, \\
   C_{qg}^{L,(1)} &=& 8 x(1-x)\,, 
   \end{array}
   :label: eq:cfFL

\ while those for :math:`F_2` read:

.. math::
   \begin{array}{rcl}
     \displaystyle \frac{C_{qq}^{2,(1)}}{2C_F} &=& \displaystyle -8\delta(1-x)\delta(1-z)+2\delta(1-x)
                                     \left(\frac{\ln
                                     (1-z)}{1-z}\right)_++ 
                                     \delta(1-x) \left[\frac{1+z^2}{1-z}\ln z+(1-z)-(1+z)\ln(1-z)\right] \\
     \\
                                 &+&\displaystyle 
                                     2\left(\frac{\ln
                                     (1-x)}{1-x}\right)_+\delta(1-z)+2\left(\frac{1}{1-x}\right)_+\left(\frac{1}{1-z}\right)_+- \left(\frac{1}{1-x}\right)_+(1+z)\\
     \\
                                 &+& \displaystyle \left[ 
                                     -\frac{1+x^2}{1-x}\ln x+(1-x)
                                     -(1+x)\ln(1-x)\right]\delta(1-z) -
                                     (1+x)\left(\frac{1}{1-z}\right)_++(2+6xz)\,,\\
     \\
     \displaystyle \frac{C_{gq}^{2,(1)}}{2C_F} &=& \displaystyle \delta (1-x) \left[\frac{1+(1-z)^2}{z} \ln\left[
                                                  z(1-z)\right]+z\right]\\
     \\
                                               &+&\displaystyle \left(\frac{1}{1-x}\right)_+\frac{1+(1-z)^2}{z}\\
     \\
                                 &+&  \displaystyle 2(1+3x)-6xz-(1+x)\frac{1}{z}\,,\\
     \\
     C_{qg}^{2,(1)} &=& \displaystyle \left[(x^2+(1-x)^2)
                        \ln\left(\frac{1-x}{x}\right)
                        +2x(1-x)\right]\delta (1-z) +(x^2+(1-x)^2)
                        \left(\frac{1}{1-z}\right)_+\\
   \\
   &+&\displaystyle 2(-1+6x-6x^2) + (x^2+(1-x)^2) \frac{1}{z} \,.
   \end{array}
   :label: eq:cfF2

By inspection of Eqs. :eq:`eq:cfFL`
and :eq:`eq:cfF2` we can deduce the various coefficients of
Eq. :eq:`eq:DoubleFuncStruct`. :math:`F_L`
involves only regular functions so that all contributions are zero but
the fully regular ones:

.. math::
   \begin{array}{rcl}
   C_{qq}^{L,(1)}(x,z) &:& K_1 = 8 C_F\,,\quad R_1^{(1)}(x) = x\,,\quad
                           R_1^{(2)}(z) = z\,,\\
   C_{gq}^{L,(1)}(x,z) &:& K_1 = 8 C_F\,,\quad R_1^{(1)}(x) = x\,,\quad
                           R_1^{(2)}(z) = 1-z\,,\\
   C_{qg}^{L,(1)}(x,z) &:& K_1 = 8\,,\quad R_1^{(1)}(x) = x(1-x)\,,\quad
                           R_1^{(2)}(z) = 1\,.\\
   \end{array}

The situation is more complicated for :math:`F_2` but we can still
identify the different contributions:

.. math::
   \begin{array}{rcl}
   C_{qq}^{2,(1)} &:& \displaystyle {\rm LL} = -16 C_F\,,\quad {\rm LS} =
                           4C_F\,,\quad {\rm LR}(z) = 2C_F
                           \left[\frac{1+z^2}{1-z}\ln
                           z+(1-z)-(1+z)\ln(1-z)\right]\\
   \\
   && \displaystyle {\rm SL} = 4C_F\,,\quad {\rm SS} = 4C_F\,,\quad {\rm
      SR}(z) = -2C_F(1+z)\\
   \\
   && \displaystyle {\rm RL}(x) = 2C_F \left[ 
                                     -\frac{1+x^2}{1-x}\ln x+(1-x)
                                     -(1+x)\ln(1-x)\right]\,,\quad {\rm
      RS}(x) = -2C_F(1+x)\,,\\
   \\
   && \displaystyle \left\{K_1 =  4C_F,\, R_1^{(1)}(x) = 1,\,
      R_1^{(2)}(z) = 1\right\}\,,\\
   \\
   && \displaystyle \left\{K_2 =  12C_F,\, R_2^{(1)}(x) = x,\,
      R_2^{(2)}(z) = z\right\}\,,\\
   \\
   \end{array}

\ 

.. math::
   \begin{array}{rcl}
   C_{gq}^{2,(1)} &:& \displaystyle {\rm LR}(z) = 2C_F\left[\frac{1+(1-z)^2}{z} \ln\left[
                                                  z(1-z)\right]+z\right]\,,\\
     \\
   && \displaystyle {\rm SR}(z) = 2C_F\left[\frac{1+(1-z)^2}{z}\right]\,,\\
   \\
   && \displaystyle \left\{K_1 =  4C_F,\, R_1^{(1)}(x) = 1 + 3x,\,
      R_1^{(2)}(z) = 1\right\}\,,\\
   \\
   && \displaystyle \left\{K_2 =  -12C_F,\, R_2^{(1)}(x) = x,\,
      R_2^{(2)}(z) = z\right\}\,,\\
   \\
   && \displaystyle \left\{K_3 =  -2C_F,\, R_3^{(1)}(x) = 1+x,\,
      R_3^{(2)}(z) = \frac1{z}\right\}\,,\\
   \\
   \end{array}

.. math::
   \begin{array}{rcl}
   C_{qg}^{2,(1)} &:&\displaystyle {\rm RL}(x) = \left[x^2+(1-x)^2\right]
                        \ln\left(\frac{1-x}{x}\right)
                        +2x(1-x)\,,\quad  {\rm RS}(x) = x^2+(1-x)^2\,,\\
   \\
   && \displaystyle \left\{K_1 = 2,\, R_1^{(1)}(x) = - 1 + 6x-6x^2,\,
      R_1^{(2)}(z) = 1\right\}\,,\\
   \\
   && \displaystyle \left\{K_2 =  1,\, R_2^{(1)}(x) = x^2+(1-x)^2,\,
      R_2^{(2)}(z) = \frac1{z}\right\}\,.
   \end{array}

Analogously, for the only LO partonic cross sections we find that:

.. math::
   \begin{array}{rcl}
   C_{qq}^{2,(0)} &:& \displaystyle {\rm LL} =1 \,.
   \end{array}

All the coefficients that are not mentioned are equal to zero. We can
now explicitly implement
Eq. :eq:`eq:MasterFormula`. The one thing that is
left to sort out is the structure of :math:`F_2` and :math:`F_L` in
terms of the appropriate PDF and FF combinations. Looking
Eq. :eq:`eq:f1sidis`, we observe that none of the
coefficient functions depends on the particular quark flavour (this is a
feature of the ZM scheme). Therefore, simplifying the notation, we can
rewrite Eq. :eq:`eq:f1sidis` as:

.. math::
   \begin{array}{rcl}
     \displaystyle  F &=& \displaystyle C_{qq}  \sum_{q} e_q^2 \left[q  D_q +\overline{q}  D_{\overline{q}}\right]  + C_{gq}  \sum_{q} e_q^2
                          \left[q+\overline{q}\right]\,D_g  +  C_{qg}  g \sum_{q} e_q^2
                          \left [D_q +D_{\overline{q}}\right]\,,
   \end{array}
   :label: eq:structF2L

where now the sums run only over the quark flavours and not over the
antiflavours.

.. container::
   :name: drell-yan-dy

   .. rubric:: Drell Yan (DY)
      :name: drell-yan-dy

In this section we apply to the Drell-Yan (DY) process the same
procedure followed above for SIDIS. As a matter of fact, SIDIS and DY
are strictly connected in that DY can be regarded as the time-like
counterpart of SIDIS. As a consequence, the structure of the relevant
observables as well as the form of the expressions involved are very
similar. Therefore, the application of the method described above is
straightforward.

.. container::
   :name: advantage-of-a-logarithmic-grid

   .. rubric:: Advantage of a logarithmic grid
      :name: advantage-of-a-logarithmic-grid

Given the particular structure of the integral :math:`I` in
Eq. :eq:`eq:ZMconv`, it turns out to be very convenient
to use a logarithmically distributed grid along with Lagrange
interpolating functions. Let us specifically consider the (massless)
integrals:

.. math::
   I_{\beta\alpha} =
   \int_{x_\beta}^1dy\,O(y)w_{\alpha}\left(\frac{x_\beta}{y}\right)\,.
   :label: eq:integralex

\ A logarithmically-spaced grid is defined such that
:math:`\ln x_{n+1}=\ln x_{n}+\delta x` with :math:`\delta x` a positive
constant. In addition, we consider a set of Lagrange interpolating
functions of degree :math:`\kappa` polynomial in :math:`\ln z`,
:math:`\{w_{\alpha}(z)\}`, that thus have the following form (see
section on the interpolation):

.. math::
   w_{\alpha}^{(k)}(z) = \sum_{j=0\atop j \leq
       \alpha}^{k}\theta(z-x_{\alpha-j})\theta(x_{\alpha-j+1}-z)\prod^{k}_{\delta=0\atop \delta\ne
       j}\left[\frac{\ln(z)-\ln(x_{\alpha-j+\delta})}{\ln(x_{\alpha})-\ln(x_{\alpha-j+\delta})}\right]\,.
   :label: eq:LagrangeFormulaLog

Due to the fact that the grid is logarithmically distributed with step
:math:`\delta x`, the function above can be rearranged as follows:

.. math:: w_{\alpha}^{(k)}(z) = \sum_{j=0\atop j \leq \alpha}^{k}\theta(z-x_{\alpha-j})\theta(x_{\alpha-j+1}-z)\prod^{k}_{\delta=0\atop \delta\ne j}\left[\frac{1}{\delta x} \ln\left(\frac{z}{x_\alpha}\right)\frac{1}{j-\delta}+1\right].

It is finally easy to see that
Eq. :eq:`eq:LagrangeFormulaLog` is such that:

.. math::
   w_{\alpha}(z) =
   \widetilde{w}\left(\ln\frac{z}{x_\alpha}\right)\quad\Rightarrow\quad
   w_{\alpha}\left(\frac{x_\beta}{y}\right) =
   \widetilde{w}\left(\ln\frac{x_\beta}{x_\alpha}-\ln y\right)
   =\widetilde{w}\left((\beta-\alpha)\delta x-\ln y\right)\,.

Therefore, the integrand of the integral in
Eq. :eq:`eq:integralex` only depends on the
difference :math:`\beta-\alpha` and not on :math:`\beta` and
:math:`\alpha` separately. Since the lower bound is :math:`x_\beta`,
this symmetry seems to broken at the level of the integral. However, the
symmetry is preserved thanks to the support properties of the
interpolating functions :math:`w_\alpha` and the fact that the relevant
functions (PDFs or FFs) are zero at :math:`x=1`. To see this, we
consider the integration limits in Eq. :eq:`eq:intlims`
with :math:`\eta=1`. They can be written as:

.. math::
   c =
     \mbox{max}(x_\beta,e^{(\beta-\alpha-1)\delta x}) \quad\mbox{and}\quad d =
     \mbox{min}(1,e^{(\beta-\alpha+\kappa)\delta x}) \,.
   :label: eq:intlims1

While the limit :math:`d` is manifestly only dependent on the difference
:math:`\beta-\alpha`, the limit :math:`c` is not. However, :math:`c`
does not have this symmetry only when :math:`x_\beta` is selected in
place of :math:`e^{(\beta-\alpha-1)\delta x}` and this can only happen
when:

.. math::
   x_\beta > e^{(\beta-\alpha-1)\delta x}\,.
   :label: eq:ineq

Since the last point of the grid is :math:`x_{N_x}=1`, being :math:`N_x`
the number of grid intervals, one can write:

.. math:: x_\beta = \frac{x_\beta}{x_{N_x}} = e^{(\beta-N_x)\delta x}.

Finally, relying on the monotonicity of the exponential function, the
inequality in Eq. :eq:`eq:ineq` becomes:

.. math:: \beta-N_x > \beta-\alpha-1\quad\Leftrightarrow\quad \alpha > N_x-1\quad\Leftrightarrow\quad \alpha = N_x\,.

Therefore, the integrals :math:`I_{\beta N_x}` do not respect the
“:math:`\beta-\alpha`” symmetry. However, as mentioned above,
:math:`I_{\beta N_x}` will always multiply a function computed in
:math:`x_{N_x}=1`. In all cases of interest this function is identically
zero at :math:`x_{N_x}=1` and thus the symmetry is effectively
preserved. In addition, :math:`c` in
Eq. :eq:`eq:intlims1` is such that if
:math:`\beta>\alpha` one has :math:`c\geq1`. But being :math:`c` the
lower integration bound of and since in
Eq. :eq:`eq:integralex` the upper bound is 1, one
immediately has that :math:`I_{\beta\alpha}=0` for :math:`\beta>\alpha`.
The consequence of these observations is that computing the integrals
:math:`a_\alpha=I_{0\alpha}` for :math:`\alpha=0,\dots,N_x` is enough to
reconstruct the full set of integrals :math:`I_{\beta\alpha}` because,
in matricial representation, :math:`I` will look like this:

.. math::
   \displaystyle I_{\beta\alpha} = 
   \begin{pmatrix}
   a_0 &  a_1 & a_2 & \cdots & a_{N_x} \\
    0  & a_0 & a_1 & \cdots & a_{N_x-1} \\
    0  & 0   &  a_0 & \cdots & a_{N_x-2} \\
   \vdots & \vdots & \vdots & \ddots & \vdots \\
    0  &   0  &   0 & \cdots & a_0 
   \end{pmatrix}\,.
   :label: eq:MatrixRep

In conclusion, adopting a logarithmically-spaced grid allows one to
compute :math:`N_x+1` integrals rather than :math:`(N_x+1)(N_x+2)/2`
integrals.

There is another aspect that matters in terms of numerical efficiency of
the computation of the integrals :math:`I_{\beta\alpha}`. Given the
support region of the interpolating functions :math:`w_\alpha`, the
integral in Eq. :eq:`eq:integralex` effectively
reads:

.. math:: I_{\beta\alpha} = \int_c^ddy\,O(y)w_{\alpha}\left(\frac{x_\beta}{y}\right)\,,

with the integration limits given in
Eq. :eq:`eq:intlims`. These limits can be rearranged as
follows:

.. math::
   c=\frac{x_\beta}{x_{{\rm min}[N_x,\alpha+1]}}\quad\mbox{and}\quad
     {d} = \frac{x_\beta}{x_{{\rm max}[\beta,\alpha-\kappa]}}\,,
   :label: eq:intlimsind

which makes it manifest the index range covered by the integration
range. The basic observation is that the functions :math:`w_\alpha` are
piecewise in correspondence of the grid nodes. This feature makes a
numerical integration over the full range defined in
Eq. :eq:`eq:intlimsind` hard to converge due to the
cusps at the grid nodes. However, the functions :math:`w_\alpha` are
smooth between two consecutive nodes. Therefore, it turns out to be
convenient to compute the integrals :math:`I_{\beta\alpha}` by breaking
the integration range as follows:

.. math::
   I_{\beta\alpha} = \sum_{j={\rm max}[0, \alpha + 1 - N_x]}^{{\rm
       min}[\kappa, \alpha - \beta]}
   \int_{x_\beta/x_{\alpha-j+1}}^{x_\beta/x_{\alpha-j}}
   dy\,O(y)w_{\alpha}\left(\frac{x_\beta}{y}\right)\,,
   :label: eq:finalformula

in such a way that the integrand of each single integral is a smooth
function and thus easier to integrate. Despite the number of intergrals
to be computed increases, this procedure makes the computation faster
and more accurate. Finally, if the grid is logarithmically distributed,
and one defines:

.. math:: s = \exp\left[\delta x\right]\,,

Eq. :eq:`eq:finalformula` can also be written as:

.. math::
   I_{\beta\alpha} = \sum_{j={\rm max}[0, \alpha + 1 - N_x]}^{{\rm
       min}[\kappa, \alpha - \beta]}  \int_{s^{\beta -\alpha + j-1}}^{s^{\beta -\alpha + j}} dy\,O(y)w_{\alpha}\left(\frac{x_\beta}{y}\right)\,.

These are the basic objects computed by ``APFEL++`` to define a
DGLAP-like operator.

GPD-related integrals
=====================

When considering computations involving GPDs, another kind of integral
structure comes into play, that is:

.. math::
   J(x) = x\int_0^1\frac{dz}{z} O \left(\frac{x}{z},x\right) d(z)\,.
   :label: eq:ZMconvERBL

\ This integral differs from that in Eq. :eq:`eq:ZMconv`
in two respects: the lower integration bound is zero rather than
:math:`x` and the operator :math:`O` may also depend explicitly on the
external variable :math:`x`. These differences make the strategy devised
for the numerical computation of convolution integrals as in
Eq. :eq:`eq:ZMconv` partially invalid. Let us discuss it
in detail to motivate the particular strategy adopted to compute
Eq. :eq:`eq:ZMconvERBL` on a grid.

As discussed above, logarithmically distributed grids are particularly
advantageous for integrals like those in
Eq. :eq:`eq:ZMconv` because they allow for a substantial
reduction of the number of integrals to be computed. However,
logarithmic grids have two main drawbacks. First, logarithmic grids that
start from a relatively low value of :math:`x` tend to be relatively
sparse at large values of :math:`x`. This is a problem because all
integrals that we are considering are such that the function being
interpolated is integrated up to :math:`x=1` and the interpolation there
can thus potentially degrade in accuracy. A possible solution to this
problem is to increase the density in a stepwise fashion as :math:`x`
gets closer to one. This produces locally logarithmically distributed
grids that allow one to exploit the symmetry discussed above while
making the grid denser, and thus more accurate, at large :math:`x`. The
implementation of this procedure in `` APFEL++`` is achieved by means of
*locked* subgrids. In practice, one starts with a logarithmic grid with
a given lower bound, *e.g.* :math:`x_{\rm min}^{(0)}=10^{-5}`. Starting
from a given node, :math:`x_{\rm min}^{(1)}`, the density of the grid is
increased by some integer factor. This procedure can be repeated an
arbitrary number of times as one moves towards large :math:`x`,
effectively defining logarithmic subgrids that are increasingly denser
and thus guarantee a better interpolation accuracy. When dealing with
integrals such as that in Eq. :eq:`eq:integralex`,
the simplest way to exploit the subgrid structure is to switch to a
denser grid at the level of the *integral*, essentially using one grid
when :math:`x_\beta` is below the transition node and the other when it
is above. The advantage of this approach is that the integration
procedure discussed above applies verbatim with the only difference
that, according to the position of :math:`x_\beta`, one grid is used
rather than another. The disadvantage of this procedure is that
integrals with low values of :math:`x_\beta` do not take advantage of
the denser grids at large :math:`x`. When dealing with functions like
PDFs that vanish rapidly as :math:`x` tends to one, this is typically
fine because the bulk of the integrals is typically due to the region at
:math:`x\gtrsim x_\beta`. This usually makes possible interpolation
inaccuracies at large :math:`x` negligible. Unfortunately, when applied
to integrals of the kind of Eq. :eq:`eq:ZMconvERBL`,
this strategy may lead to severe inaccuracies at large values of
:math:`x_\beta`. This is due to the fact that the integral extends down
to zero and this procedure does not make use of the low-:math:`x` grids,
effectively truncating the integral to increasingly larger values of
:math:`x` as :math:`x_\beta` approaches one.

The second problem is that logarithmic grids get down to :math:`x = 0`.
This is not an issue for the computation in
Eq. :eq:`eq:ZMconv` because the lower integration bound
is :math:`x` and thus, as long as the integral in not computed below the
lower bound of the grid, this does not introduce any inaccuracy.
Conversely, this is a potential problem for integrals like those in
Eq. :eq:`eq:ZMconvERBL` that extend down to zero for
any value of :math:`x`. This may suggest that logarithmic grids are not
suitable in the first place for integrals as in
Eq. :eq:`eq:ZMconvERBL`. However, we observe that
GPD-related integrands are usually well-behaved around :math:`y=0`
implying that a logarithmic grid with lower bound :math:`x_0` close
enough to zero (*e.g.* :math:`x_0=10^{-5}`) is expected to be enough to
make the contribution to the integral due to the region
:math:`y\in[0,x_0]` negligibly small. Therefore, also for the
computation of Eq. :eq:`eq:ZMconvERBL`, we stick to
logarithmically distributed grids which allows us to use the same class
for constructing and managing grids developed to compute
Eq. :eq:`eq:ZMconv`.

In addition, GPD-related integrals, due to the explicit dependence on
the external variable of the operator (see
Eq. :eq:`eq:ZMconvERBL`), do not enjoy
:math:`\beta-\alpha` symmetry discussed above. This is a major
limitation in that it does not allow us to reduce the number of
integrals to be computed.

These observation, along with the discussion on the interpolation grid,
leads us to take an approach in which the operator is computed for all
possible pairs :math:`(\alpha,\beta)` over a grid that extends as low as
possible in :math:`x` and is dense at large :math:`x`: *i.e.* the joint
grid.

Using the usual interpolation formula gives us that the integral in
Eq. :eq:`eq:ZMconvERBL` is computed on the grid as
follows:

.. math:: J(x_\beta) = \sum_{\alpha=}^{N_x-1}J_{\beta\alpha}\overline{d}_\alpha\,.

with :math:`\overline{d}_\alpha = x_\alpha d(x_\alpha)` and the
assumption :math:`d(1)=0`, and with:

.. math::
   J_{\beta\alpha} =
   \sum_{j=0}^{{\rm min}[\alpha,\kappa]}\int_{x_\beta/x_{\alpha-j+1}}^{x_\beta/x_{\alpha-j}} dy\,O(y,x_\beta)w_\alpha^{(\kappa)}\left(\frac{x_\beta}{y}\right)\,,

where the indices :math:`\alpha` and :math:`\beta` are intended to run
over the joint grid. We point out again that in this case it is
necessary to compute all the :math:`N_x^2` integrals making up
:math:`J_{\beta\alpha}`.

.. container::
   :name: principal-valued-integrals

   .. rubric:: Principal-valued integrals
      :name: principal-valued-integrals

Expressions for the operator :math:`O` in
Eq. :eq:`eq:ZMconvERBL` often involve singular terms
that are to interpreted as Cauchy principal-valued distributions. We are
specifically concerned with integrals of the following kind:

.. math:: I={\rm PV}\int_x^{\infty}dy\,\frac{f(y)}{1-y}\,

\ where :math:`f(y)` is a smooth function over the integration range. In
order to numerically treat this kind of integrals, in the document
devoted to GPDs the so-called :math:`++`-distribution has been
introduced such that:

.. math:: I=\int_x^{\infty}dy\,\left(\frac{1}{1-y}\right)_{++}f(y)\,

Similarly to the popular :math:`+`-distribution, the
:math:`++`-distribution acts upon integration in the following way:

.. math::
   \int_x^{\infty}dy\,\left(\frac{1}{1-y}\right)_{++}f(y) =
   \int_x^{\infty}\frac{dy}{1-y}\left[f(y)-f(1)\left(1+\theta(y-1)\frac{1-y}{y}\right)\right]+f(1)\ln(1-x)
   \,.
   :label: eq:pluplusprescription

However, differently from the :math:`+`-distribution, the
:math:`++`-distribution does not regularise an otherwise divergen
integral. It rather rearranges the integrand in such a way that it is
finite over the integration range (the divergence at :math:`y=1` is no
longer there) making it numerically treatable. We now need to devise a
way to compute this integral within the interpolation framework. Indeed,
the interpolation procedure produces integrals as this:

.. math:: I_{\beta\alpha}=\int_{x_{\beta}}^{\infty}dy\,\left(\frac{1}{1-y}\right)_{++}w_{\alpha}\left(\frac{x_{\beta}}{y}\right) \,,

that need to be treated with care due to the support region of the
interpolation functions :math:`w_\alpha`. In particular, we know that
the function :math:`w_\alpha(x_\beta/y)` has support
:math:`y\in [x_\beta/x_{\alpha+1},x_\beta/x_{\alpha-\kappa}]`, so that:

.. math::
   \begin{array}{rcl}
   I_{\beta\alpha}&=&\displaystyle
                      \int_{c}^{d}dy\,\left(\frac{1}{1-y}\right)_{++}w_{\alpha}\left(\frac{x_{\beta}}{y}\right)=\int_{c}^{\infty}dy\,\left(\frac{1}{1-y}\right)_{++}w_{\alpha}\left(\frac{x_{\beta}}{y}\right)\\
   \\
   &=&\displaystyle
       \int_c^{\infty}\frac{dy}{1-y}\left[w_{\alpha}\left(\frac{x_{\beta}}{y}\right)-\delta_{\beta\alpha}\left(1+\theta(y-1)\frac{1-y}{y}\right)\right]+\delta_{\beta\alpha}\ln(1-c)\\
   \\
   &=&\displaystyle
       \int_c^{d}\frac{dy}{1-y}\left[w_{\alpha}\left(\frac{x_{\beta}}{y}\right)-\delta_{\beta\alpha}\left(1+\theta(y-1)\frac{1-y}{y}\right)\right]+\delta_{\beta\alpha}\left[\ln(1-c)-\int_d^{\infty}\frac{dy}{1-y}\left(1+\theta(y-1)\frac{1-y}{y}\right)\right]\\
   \\
   &=&\displaystyle
       \int_c^{d}\frac{dy}{1-y}\left[w_{\alpha}\left(\frac{x_{\beta}}{y}\right)-\delta_{\beta\alpha}\left(1+\theta(y-1)\frac{1-y}{y}\right)\right]+\delta_{\beta\alpha}\left[\ln(1-c)-\int_d^{\infty}dy\left(\frac{1}{1-y}+\frac{1}{y}\right)\right]\\
   \\
   &=&\displaystyle
       \int_c^{d}\frac{dy}{1-y}\left[w_{\alpha}\left(\frac{x_{\beta}}{y}\right)-\delta_{\beta\alpha}\left(1+\theta(y-1)\frac{1-y}{y}\right)\right]+\delta_{\beta\alpha}\left[\ln(1-c)-\int^{1-d}_{-\infty}\frac{dz}{z}-\int_d^{\infty}\frac{dy}{y}\right] \\
   \\
   &=&\displaystyle
       \int_c^{d}\frac{dy}{1-y}\left[w_{\alpha}\left(\frac{x_{\beta}}{y}\right)-\delta_{\beta\alpha}\left(1+\theta(y-1)\frac{1-y}{y}\right)\right]+\delta_{\beta\alpha}\left[\ln(1-c)+\int_{1-d}^{d}\frac{dy}{y}\right] \\
   \\
   &=&\displaystyle
       \int_c^{d}\frac{dy}{1-y}\left[w_{\alpha}\left(\frac{x_{\beta}}{y}\right)-\delta_{\beta\alpha}\left(1+\theta(y-1)\frac{1-y}{y}\right)\right]+\delta_{\beta\alpha}\left[\ln(1-c)+\int_{d-1}^{d}\frac{dy}{y}\right] \\
   \\
   &=&\displaystyle
       \int_c^{d}\frac{dy}{1-y}\left[w_{\alpha}\left(\frac{x_{\beta}}{y}\right)-\delta_{\beta\alpha}\left(1+\theta(y-1)\frac{1-y}{y}\right)\right]+\delta_{\beta\alpha}\left[\ln(1-c)-\ln\left(1-\frac{1}{d}\right)\right] \,,
   \end{array}

with:

.. math:: c = \frac{x_\beta}{x_{\alpha+1}}\quad\mbox{and}\quad d = \frac{x_\beta}{x_{\alpha-\rm{min}[\alpha,\kappa]}}\,,

and where we have used the fact that :math:`d\geq 1` for
:math:`\beta=\alpha` and, according to the Cauchy principal-value
prescription:

.. math:: \int_{-a}^{a}\frac{dy}{y} = 0\,,

where :math:`a` could also be :math:`\infty`.

**References**

**References**

.. container:: references csl-bib-body hanging-indent
   :name: refs

   .. container:: csl-entry
      :name: ref-deFlorian:1997zj

      Florian, D. de, M. Stratmann, and W. Vogelsang. 1998. “QCD
      analysis of unpolarized and polarized Lambda baryon production in
      leading and next-to-leading order.” *Phys. Rev. D* 57: 5811–24.
      https://doi.org/10.1103/PhysRevD.57.5811.

.. [1]
   More precisely, :math:`O` is in general a distribution that may
   contain :math:`\delta` functions and :math:`+`-distributions.

.. [2]
   This is the case for SIDIS and DY up to NLO. However, one may expect
   that this feature holds also beyond.
