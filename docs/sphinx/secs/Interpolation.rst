=============
Interpolation
=============

.. contents::
   :depth: 3
..

In this part of the documentation we address the question of
interpolation. Interpolation is an extremely powerful tool that allows
one to reconstruct, within some accuracy, continuos functions when
knowing them only in a finite number of points. Unsurprisingly,
interpolation techniques are massively used in numerical codes that deal
with physical quantities that are continuos functions of space-time
variables, such as position and momentum. In this respect, ``APFEL++``
makes no exception.

In order to make the best of interpolation, ``APFEL++`` has been
designed around a specific choice of the interpolation strategy.
Specifically, the very computational core of ``APFEL++`` relies on
*Langrange* polynomials and their properties. Despite Langrange
polynomials do not enjoy the smoothness of splines, they enjoy a number
of very useful properties that allow extending the use of Lagrange
polynomials from just interpolation to derivation and integration. In
the following, the interpolation strategy is detailed along with its
implementation and applications ``APFEL++``.

.. _sec:LagrangeInterpolation:

Lagrange interpolation
======================

In this section we will derive a general expression for the Lagrange
interpolating functions :math:`w`. These functions are ubiquitous in
`` APFEL++`` in that they are used, not only for simple interpolations,
but also for computing convolutions, integrals in general, and also
derivatives. It is then important to give them a thorough derivation in
terms of Lagrange polynomial and to understand how they behave upon
derivation and integration.

Suppose one wants to interpolate the test function :math:`g` in the
point :math:`x` using a set of Lagrange polynomials of degree :math:`k`.
This requires a subset of :math:`k+1` consecutive points on an
interpolation grid, say :math:`\{x_{\alpha},\dots,x_{\alpha+k}\}`. The
relative position between the point :math:`x` and the subset of points
used for the interpolation is arbitrary. It is convenient to choose the
subset of points such that :math:`x_\alpha < x \leq x_{\alpha+k}`. [1]_
However, the ambiguity remains because there are :math:`k` possible
choices according to whether :math:`x_\alpha < x \leq x_{\alpha+1}`, or
:math:`x_{\alpha+1} < x \leq x_{\alpha+2}`, and so on. For now we assume
that:

.. math::
   x_{\alpha} < x \leq x_{\alpha+1}\,,
   :label: eq:assumption1

\ but we will release this assumption below. Using the standard Lagrange
interpolation procedure, one can approximate the function :math:`g` in
:math:`x` as:

.. math::
   g(x) = \sum_{i=0}^k\ell_i^{(k)}(x)g(x_{\alpha+i})\,,
   :label: particularCase

where :math:`\ell_i^{(k)}` is the :math:`i`-th Lagrange polynomial of
degree :math:`k` which can be written as:

.. math::
   \ell_i^{(k)}(x) = \prod^{k}_{m=0\atop m\ne
   i}\frac{x-x_{\alpha+m}}{x_{\alpha+i}-x_{\alpha+m}}\,.
   :label: eq:LagPoly

Since we have assumed that :math:`x_\alpha < x \leq x_{\alpha+1}` (see
Eq. :eq:`eq:assumption1`),
Eq. :eq:`particularCase` can be written as:

.. math::
   g(x) =
     \theta(x-x_{\alpha})\theta(x_{\alpha+1}-x)\sum_{i=0}^k
     g(x_{\alpha+i})\prod^{k}_{m=0\atop m\ne
       i}\frac{x-x_{\alpha+m}}{x_{\alpha+i}-x_{\alpha+m}}\,.
   :label: particularCaseTheta

In order to make Eq. :eq:`particularCaseTheta`
valid for all values of :math:`\alpha`, one just has to sum over all
:math:`N_x+1` nodes of the *global* interpolation grid
:math:`\{x_0,\dots,x_{N_x}\}`, that is:

.. math::
   g(x) =
     \sum_{\alpha=0}^{N_x-1}\theta(x-x_{\alpha})\theta(x_{\alpha+1}-x)\sum_{i=0}^k
     g(x_{\alpha+i})\prod^{k}_{m=0\atop m\ne i}\frac{x-x_{\alpha+m}}{x_{\alpha+i}-x_{\alpha+m}}\,.
   :label: generalCase

\ Defining :math:`\beta=\alpha+i`, one can rearrange the equation above
as:

.. math::
   g(x) =
     \sum_{\beta=0}^{N_x+k-1}w_\beta^{(k)}(x) g(x_{\beta})\,,
   :label: generalCase2

that leads to the definition of the interpolating functions:

.. math::
   w_\beta^{(k)}(x) = \sum_{i=0\atop i\leq\beta}^k
     \theta(x-x_{\beta-i})\theta(x_{\beta-i+1}-x) \prod^{k}_{m=0\atop m\ne i}\frac{x-x_{\beta-i+m}}{x_{\beta}-x_{\beta-i+m}}\,.
   :label: eq:intfunc

Notice that the condition :math:`i\leq\beta` comes from the condition
:math:`\alpha\geq 0`. It is important to observe that the sum in
Eq. :eq:`generalCase2` extends up to the
:math:`(N_x+k-1)`-th node. Therefore, the original grid needs to be
extended by :math:`k-1` nodes. However, the range of validity of the
interpolation remains that defined by the original grid, *i.e.*
:math:`x_0 \leq x \leq x_{N_x}`.

When implementing this interpolation procedure, it is important to
realise that, typically, only a small number of terms in the sum in
Eq. :eq:`generalCase2` is different from zero. For any
given value of :math:`x`, it is possible to determine the values of the
index :math:`\beta` for which the interpolating functions
:math:`w_\beta^{(k)}` are different from zero, reducing (often
dramatically) the amount of sums required to carry out an interpolation.
The range of :math:`\beta` is easily determined by observing that in
Eq. :eq:`particularCaseTheta` the summation
extends on the nodes between :math:`x_\alpha` and :math:`x_{\alpha+k}`.
But since :math:`\beta` is defined like :math:`\alpha+i` this exactly
defines the range in :math:`\beta`:

.. math:: \alpha(x) \leq \beta \leq \alpha(x) + k\,,

\ where the function :math:`\alpha(x)` is implicitly defined through
Eq. :eq:`eq:assumption1`. Therefore,
Eq. :eq:`generalCase2` becomes:

.. math::
   g(x) =
     \sum_{\beta=\alpha(x)}^{\alpha(x)+k}w_\beta^{(k)}(x) g(x_{\beta})\,,
   :label: eq:limitedsum

Since the interpolation functions :math:`w_\beta^{(k)}(x)` often appear
inside integrals, it is very useful to use the fact that they are
different from zero only over a limited interval, specifically:

.. math::
   w_\beta^{(k)}(x) \neq 0\quad \Leftrightarrow\quad
   x_{\beta-k}<x < x_{\beta+1}\,.
   :label: eq:limits

\ This allows one to optimise the integration restricting the
integration region only to where the interpolating functions are
different from zero.

.. container::
   :name: generalised-interpolation

   .. rubric:: Generalised interpolation
      :name: generalised-interpolation

In most of the applications within ``APFEL++`` the assumption in
Eq. :eq:`eq:assumption1` is used. However, sometimes
it may be useful to release this assumption. A situation in which this
is advantageous is in the presence of non-smooth or discontinuos
functions (such as PDFs and FFs as function of the factorisation scale
:math:`\mu` in correspondence of the heavy-quark thresholds). When
interpolating these functions one should not interpolate over the
discontinuities. To do so and yet retain a given interpolation degree,
one can release the assumption in
Eq. :eq:`eq:assumption1`. Specifically, we
generalise it to:

.. math::
   x_{\alpha+t} < x \leq
     x_{\alpha+t+1}\quad\mbox{with}\quad t = 0,\dots,k-1\,,
   :label: IntAssumptionGen

\ such that the interpolation formula becomes:

.. math::
   g(x) =
     \sum_{\alpha=-t}^{N_x-t-1}\theta(x-x_{\alpha+t})\theta(x_{\alpha+t+1}-x)\sum_{i=0}^k
     g(x_{\alpha+i})\prod^{k}_{m=0\atop m\ne i}\frac{x-x_{\alpha+m}}{x_{\alpha+i}-x_{\alpha+m}}\,,
   :label: MoreGeneralCase

that can be rearranged as:

.. math::
   g(x) =
   \sum_{\beta=-t}^{N_x+k-t-1}w_{\beta,t}^{(k)}(x) g(x_{\beta})\,,
   :label: generalCase3

with:

.. math::
   w_{\beta,t}^{(k)}(x) = \sum_{i=0,i\leq\beta}^k
   \theta(x-x_{\beta-i+t})\theta(x_{\beta-i+t+1}-x) \prod^{k}_{m=0,m\ne
   i}\frac{x-x_{\beta-i+m}}{x_{\beta}-x_{\beta-i+m}}\,,
   :label: eq:generalisedintfuncs

being the “generalised” interpolation functions. We observe that the
support region of :math:`w_{\beta,t}^{(k)}` is:

.. math:: w_{\beta,t}^{(k)}(x)\neq 0 \quad\Leftrightarrow\quad x_{\beta+t-k} < x < x_{\beta+t+1}\,,

that generalises Eq. :eq:`eq:limits`. The generalised
interpolation functions can be used to avoid interpolating over some
particular grid nodes. In order to avoid interpolation over a specific
node of the grid, one can choose :math:`t` dynamically in such a way
that :math:`\beta+t` in
Eq. :eq:`eq:generalisedintfuncs` never
corresponds to that particular node. This mechanism is implemented in
``APFEL++`` as follows. The interpolation grid is chosen to have two
nodes in correspondence of the *threshold* :math:`x_T`, but slightly
displaced up and down by an “infinitesimal” amount :math:`\epsilon` to
keep them separate, that is:

.. math:: \{x_0,\dots,x_{\alpha_t-1},x_{\alpha_t},\dots,x_{N_x}\}\quad\mbox{with}\quad x_{\alpha_t-1}=x_T-\epsilon\quad\mbox{and}\quad x_{\alpha_t}=x_T+\epsilon\,.

The aim is then to avoid interpolating over the nodes
:math:`x_{\alpha_t-1}` and :math:`x_{\alpha_t}`. By default we assume
:math:`t=0` in
Eq. :eq:`eq:generalisedintfuncs` so that we
automatically reduce to Eq. :eq:`eq:intfunc`. In this
situation, we are assuming Eq. :eq:`eq:assumption1`
where effectively the index :math:`\alpha` is determined dynamically
according to the values of :math:`x`. Therefore, we can effectively
write:

.. math:: x_{\alpha(x)} < x \leq x_{\alpha(x)+1}\,,

which implicitly defines the function :math:`\alpha(x)`.
Eq. :eq:`particularCaseTheta` then requires
summing over the :math:`k+1` nodes of the grid
:math:`\{x_{\alpha(x)},\dots,x_{\alpha(x)+k}\}`. However, when the point
:math:`x` approaches :math:`x_T` from below, the range
:math:`\{x_{\alpha(x)},\dots,x_{\alpha(x)+k}\}` may end up enclosing
both nodes :math:`x_{\alpha_t-1}` and :math:`x_{\alpha_t}`. To avoid
this, we promote the index :math:`t` in
Eq. :eq:`IntAssumptionGen` to a function of
:math:`x` defined through the inequalities:

.. math::
   \left\{\begin{array}{l}
   x< x_T\,,\\
   \\
   x_{\alpha_t-2} < x_{\alpha(x)-t(x)+k}\leq x_{\alpha_t-1}\,,
   \end{array}\right.

that translate into:

.. math::
   \left\{\begin{array}{l}
   \alpha(x) \leq \alpha_t -2\,,\\
   \\
   \alpha_t-2 <  \alpha(x)-t(x)+k \leq \alpha_t-1\,.
   \end{array}\right.

Imposing the unnecessary but convenient constraint :math:`t(x)\geq 0`,
finally gives:

.. math::
   t(x) = \mbox{max}\left[\mbox{min}\left[\alpha(x),\alpha_t -2\right]
     -\alpha_t+k + 1, 0\right] \,,
   :label: eq:tdef

that also obeys:

.. math:: 0\leq t(x) \leq k - 1\,,

as required. In addition, as in
Eq. :eq:`eq:limitedsum`, the summation over
:math:`\beta` in Eq. :eq:`generalCase3` can be
restricted to a range of :math:`k+1` nodes as:

.. math::
   g(x) =
   \sum_{\beta=\alpha(x)-t(x)}^{\alpha(x)-t(x)+k}w_{\beta,t}^{(k)}(x)
   g(x_{\beta})\,.
   :label: generalCase4

.. container::
   :name: bi-dimensional-interpolation

   .. rubric:: Bi-dimensional interpolation
      :name: bi-dimensional-interpolation

As discussed in the section devoted to the computation of the
convolution integrals, the interpolation functions can be used to
compute integrals of the following kind:

.. math:: I_1 = \int_{x_0}^{x_{N_x}}dx\,g(x)f(x)\,,

\ where :math:`f` is a smooth function. Using
Eqs. :eq:`generalCase2`
and :eq:`eq:limits` we have that:

.. math:: I_1 = \sum_{\beta=0}^{N_x+k-1} W_\beta g(x_{\beta})\,,

with:

.. math::
   W_\beta = \int_{x_{{\rm max}(0,\beta-k)}}^{x_{{\rm
         min}(N_x,\beta+1)}}dx \,w_\beta^{(k)}(x)f(x)\,.
   :label: eq:monodim

The equation above can be easily generalised to a bidimensional integral
as:

.. math:: I_2 = \int_{x_0}^{x_{N_x}}dx \int_{y_0}^{y_{N_y}}dy\,g(x,y)f(x,y) = \sum_{\alpha=0}^{N_x+k-1} \sum_{\beta=0}^{N_y+l-1} W_{\alpha\beta} g(x_{\alpha},y_{\beta})\,,

with:

.. math::
   W_{\alpha\beta} = \int_{x_{{\rm max}(0,\alpha-k)}}^{x_{{\rm
         min}(N_x,\alpha+1)}}dx \int_{y_{{\rm max}(0,\beta-k)}}^{y_{{\rm
         min}(N_y,\beta+1)}}dy
   \,w_\alpha^{(k)}(x)\,w_\beta^{(l)}(y)\,f(x,y)\,.
   :label: eq:bidim

As discussed below in much detail, we mention that the functions
:math:`w` are piecewise. In particular, while they are continuous in
correspondence of the nodes of the grid, their first derivative is not.
As a consequence, the result of the numerical integrals in
Eqs. :eq:`eq:monodim` and :eq:`eq:bidim`
may be inaccurate. To overcome this problem, it is sufficient to split
the integrals in sub-integrals over the intervals delimited by two
consecutive nodes. Using Eq. :eq:`eq:limits`, it is easy
to see that, for an interpolation of degree :math:`k`, one needs to
compute :math:`k+1` integrals over the intervals included between the
:math:`(\beta-k)`-th and the :math:`(\beta+1)`-th node.

.. container::
   :name: derivative-and-integral-of-an-interpolated-function

   .. rubric:: Derivative and integral of an interpolated function
      :name: derivative-and-integral-of-an-interpolated-function

The simple form of the Lagrange polynomials allows, in some cases, for
the analytic handling of operations such has derivation and integration
of the interpolated functions. In this section we discuss how derivation
and integration can be carried out analytically exploiting some specific
properties of the Lagrange polynomials.

Referring Eq. :eq:`particularCase`, the main
observation is that the interpolating functions :math:`\ell_i^{(k)}` are
solutions of the following differential-equation system:

.. math::
   \left\{
   \begin{array}{l}
   \displaystyle \frac{d\ell_i^{(k)}}{dx} = \left(\sum_{n=0 \atop n\neq
     i}^k\frac{1}{x-x_{\alpha+n}}\right) \ell_i^{(k)}(x)\\
   \\
   \ell_i^{(k)}(x_{\alpha+i}) = 1
   \end{array}
   \right.\,,
   :label: eq:diffeqlagpol

\ whose solution is Eq. :eq:`eq:LagPoly`. This allows us
to compute the derivative of the test function :math:`g` in
Eq. :eq:`eq:limitedsum` only knowing its values on
the grid points as:

.. math::
   \frac{dg}{dx} =
     \sum_{\beta=\alpha(x)}^{\alpha(x)+k}\mathcal{D}_\beta^{(k)}(x)
     g(x_{\beta})\,,
   :label: DerGeneralCase2

where:

.. math::
   \mathcal{D}_\beta^{(k)}(x) = \sum_{i=0}^{{\rm min}(k,\beta)}
     \theta(x-x_{\beta-i})\theta(x_{\beta-i+1}-x) \left(\sum_{n=0 \atop n\neq
     i}^k\frac1{x_\beta-x_{\beta-i+n}}\prod^{k}_{m=0\atop m\ne
       i,n}\frac{x-x_{\beta-i+m}}{x_{\beta}-x_{\beta-i+m}}\right)\,.
   :label: eq:derintfunc

It should be pointed out that, due to the discontinuity of the
derivative of the interpolation functions on the grid nodes, numerical
derivation through Eq. :eq:`DerGeneralCase2` does
not work when :math:`x` coincides with a grid node. This problem can be
overcome by simply slightly displacing the value of :math:`x` in case in
falls on a grid node.

Eq. :eq:`eq:diffeqlagpol` can also be used to
compute integrals of the function :math:`g`. In particular, suppose we
want to compute the integral:

.. math::
   I(a,b) = \int_a^b dy\,g(y)\,.
   :label: eq:IntegralInt

\ Using Eq. :eq:`eq:limitedsum`, we find that:

.. math::
   I(a,b)=
     \int_a^b dy
     \sum_{\beta=\alpha(y)}^{\alpha(y)+k}w_\beta^{(k)}(y)
     g(x_{\beta}) = \int_a^b dy\left[w_{\alpha(y)}^{(k)}(y)
     g(x_{\alpha(y)})+\dots+w_{\alpha(y)+k}^{(k)}(y)
     g(x_{\alpha(y)+k})\right]\,.

The dependence on the integration variable :math:`y` of the index
:math:`\alpha` complicates the solution of this integral. However, it is
easy to derive a close form for the function :math:`\alpha(y)`:

.. math:: \alpha(y) = -1 + \sum_{\beta = 0}^{N_x}\theta(y-x_\beta)\,,

that (obviously) means that :math:`\alpha(y)` is constant on the
separate intervals :math:`[x_0:x_1]`, :math:`[x_1:x_2]`, and so on. This
allows us to break the integral above as follows:

.. math::
   \int_a^b dy \sum_{\beta=\alpha(y)}^{\alpha(y)+k}w_\beta^{(k)}(y)
     g(x_{\beta})= \left[\int_a^{x_{\alpha(a)+1}}+ \sum_{\gamma=\alpha(a)+1}^{\alpha(b)-1}\int_{x_\gamma}^{x_{\gamma+1}}+\int_{x_{\alpha(b)}}^b\right]dy \sum_{\beta=\alpha(y)}^{\alpha(y)+k}w_\beta^{(k)}(y)
     g(x_{\beta})\,,

such that in each single integral the integrand has a constant value of
:math:`\alpha`. The three terms above can be separately rearranged as:

.. math:: \sum_{\beta=\alpha(a)}^{\alpha(a)+k}g(x_{\beta})\int_a^{x_{\alpha(a)+1}}dy\,w_\beta^{(k)}(y)\,,

.. math:: \sum_{\gamma=\alpha(a)+1}^{\alpha(b)-1}\sum_{\beta=\gamma}^{\gamma+k}g(x_{\beta})\int_{x_\gamma}^{x_{\gamma+1}}dy\,w_\beta^{(k)}(y)\,,

.. math::
   \sum_{\beta=\alpha(b)}^{\alpha(b)+k}g(x_{\beta})\int_{x_{\alpha(b)}}^bdy\,
   w_\beta^{(k)}(y)\,,

where we have used the fact that :math:`\alpha(x_\gamma)=\gamma`. Notice
that with this we managed to obtain direct integrations of the
interpolating functions :math:`w_\beta^{(k)}`. Now we note that the
three terms above span the range of nodes
:math:`\beta\in[\alpha(a):\alpha(b)+k]`. Therefore, the integral in
Eq. :eq:`eq:IntegralInt` can be written as:

.. math::
   I(a,b) =
     \sum_{\beta=\alpha(a)}^{\alpha(b)+k}\mathcal{G}_\beta^{(k)}(a,b)
     g(x_{\beta})\,,
   :label: eq:integralT0

with:

.. math::
   \begin{array}{rcl}
   \displaystyle \mathcal{G}_\beta^{(k)}(a,b) &=& \displaystyle
                                                  \int_{{\rm
                                                  max}(a,x_{\beta-k})}^{{\rm
                                                  min}(b,x_{\beta+1})}
                                                  dy\,w_\beta^{(k)}(y)\\
   \\
   &=&\displaystyle 
   \sum_{i=0}^{{\rm min}(k,\beta)}
     \int_{{\rm max}(a,x_{\beta-k})}^{{\rm min}(b,x_{\beta+1})}dy\,\theta(y-x_{\beta-i})\theta(x_{\beta-i+1}-y)
     \prod^{k}_{m=0\atop m\ne
       i}\frac{y-x_{\beta-i+m}}{x_{\beta}-x_{\beta-i+m}}\,.
   \end{array}

\ Now, assuming that :math:`a\leq b`, we can write the integral above
as:

.. math::
   \int_{{\rm max}(a,x_{\beta-k})}^{{\rm min}(b,x_{\beta+1})} dy\,\theta(y-x_{\beta-i})\theta(x_{\beta-i+1}-y)\dots
   =\theta(b-x_{\beta-i})\theta(x_{\beta-i+1}-a)\int_{{\rm max}(a,x_{\beta-i})}^{{\rm min}(b,x_{\beta-i+1})}dy\dots\,,

so that:

.. math::
   \mathcal{G}_\beta^{(k)}(a,b) = \sum_{i=0}^{{\rm min}(k,\beta)}
     \theta(b-x_{\beta-i})\theta(x_{\beta-i+1}-a)\prod^{k}_{m=0\atop m\ne
       i}\frac{1}{x_{\beta}-x_{\beta-i+m}}\int_{{\rm max}(a,x_{\beta-i})}^{{\rm min}(b,x_{\beta-i+1})}dy\prod^{k}_{m=0\atop m\ne
       i}(y-x_{\beta-i+m})\,.

It is easy to see that:

.. math::
   \prod^{k}_{m=0\atop m\ne
   i}(y-x_{\beta-i+m})=\sum_{n=0}^k (-1)^n p_{\beta}^{(k)}(n)y^{k-n}\,,
   :label: eq:expansion

such that:

.. math::
   \begin{array}{rcl}
   \displaystyle \mathcal{G}_\beta^{(k)}(a,b) &=&\displaystyle  \sum_{i=0}^{{\rm min}(k,\beta)}
     \theta(b-x_{\beta-i})\theta(x_{\beta-i+1}-a)\\
   \\
   &\times&\displaystyle \left[\prod^{k}_{m=0\atop m\ne
       i}\frac{1}{x_{\beta}-x_{\beta-i+m}}\right]\sum_{n=0}^k
            \frac{(-1)^n p_{\beta}^{(k)}(n)}{k-n+1}
            \left(\overline{b}^{k-n+1}-\overline{a}^{k-n+1}\right)\,.
   \end{array}
   :label: eq:linIntegral

with :math:`\overline{a} = {\rm max}(a,x_{\beta-i})` and
:math:`\overline{b} = {\rm min}(b,x_{\beta-i+1})`. What is left to do is
to determine the coefficients :math:`p_{\beta}^{(k)}`. For convenience,
let us define the set
:math:`\{r_1,\dots,r_k\}=\{x_{\beta-i},\dots,x_{\beta-1},x_{\beta+1},\dots,x_{\beta-i+k}\}`. [2]_
The coefficients :math:`p_{\beta}^{(k)}` defined in
Eq. :eq:`eq:expansion` can be expressed as:

.. math::
   p_{\beta}^{(k)}(n)=\sum_{i_1=1}^k r_{i_1}\sum_{i_2=i_1+1}^k
     r_{i_2}\dots\sum_{i_n=i_{n-1}+1}^k r_{i_n}\,.
   :label: eq:pncoef

In order to obtain a convenient algorithm to compute the coefficients
:math:`p_{\beta}^{(k)}`, we define the vectorial function:

.. math:: \mathbf{f}^{(k)}(\mathbf{r},\mathbf{a})\quad\mbox{with components}\quad  f_j^{(k)}\left(\mathbf{r},\mathbf{a}\right)=\sum_{i=j+1}^k r_{i}a_{i}\,.

It is important to notice that while :math:`\mathbf{r}` is a
:math:`k`-dimensional vector with index running between 1 and :math:`k`,
:math:`\mathbf{a}` and :math:`\mathbf{f}^{(k)}` are in principle
infinite-dimensional vectors with index running between :math:`-\infty`
and :math:`+\infty`. However, given the definition of its components in
Eq. :eq:`eq:pncoef`, it turns out that
:math:`f_j^{(k)}=0` for :math:`j\geq k`. The same applies to
:math:`\mathbf{a}` because, as we will see below, it has to be
identified with some :math:`\mathbf{f}^{(k)}`. The function
:math:`\mathbf{f}^{(k)}` can now be used to define the vector
:math:`\mathbf{P}^{(k)}(n)` recursively. The relevant recursive relation
is:

.. math:: \mathbf{P}^{(k)}(n+1) = \mathbf{f}^{(k)}(\mathbf{r},\mathbf{P}^{(k)}(n))\quad\mbox{with}\quad\mathbf{P}^{(k)}(0) = \mathbf{1}\,.

Finally, the coefficient :math:`p_{\beta}^{(k)}(n)` is the zero-th
component of the vector :math:`\mathbf{P}^{(k)}(n)`, *i.e.*
:math:`p_{\beta}^{(k)}(n)\equiv P_0^{(k)}(n)`.

Now, we turn to derive derivation and integration formulas for the
generalised interpolation formula in
Eq. :eq:`generalCase4`. The derivative is as simple as
in the :math:`t=0` case, that is:

.. math::
   \frac{dg}{dx} =
   \sum_{\beta=\alpha(x)-t(x)}^{\alpha(x)-t(x)+k}\mathcal{D}_{\beta,t}^{(k)}(x) g(x_{\beta})\,,

with:

.. math::
   \mathcal{D}_{\beta,t}^{(k)}(x) = \sum_{i=0}^{{\rm min}(k,\beta)}
   \theta(x-x_{\beta-i+t(x)})\theta(x_{\beta-i+t(x)+1}-x) \left(\sum^{k}_{n=0\atop m\ne
   i}\frac{1}{x_{\beta}-x_{\beta-i+n}}\prod^{k}_{m=0\atop m\ne
   i,n}\frac{x-x_{\beta-i+m}}{x_{\beta}-x_{\beta-i+m}}\right)\,.

Due to the presence of the function :math:`t(x)`, defined in
Eq. :eq:`eq:tdef`, the integration procedure is more
involved. We want to compute:

.. math::
   \begin{array}{rcl}
     I(a,b)&=&\displaystyle
     \int_a^b dy
     \sum_{\beta=\alpha(y)-t(y)}^{\alpha(y)-t(y)+k}w_{\beta,t}^{(k)}(y)
     g(x_{\beta}) \\
   \\
   &=&\displaystyle \left[\int_a^{x_{\alpha(a)+1}}+ \sum_{\gamma=\alpha(a)+1}^{\alpha(b)}\int_{x_\gamma}^{x_{\gamma+1}}-\int_b^{x_{\alpha(b)+1}}\right]dy \sum_{\beta=\alpha(y) -t(y)}^{\alpha(y)-t(y)+k}w_{\beta,t}^{(k)}(y)
     g(x_{\beta})\,.
   \end{array}

\ As in the case of :math:`t=0`, each of the integrals above has a
constant value of :math:`\alpha(y)` and :math:`t(y)` so that sum over
:math:`\beta` and the integral sign can be exchanged:

.. math::
   \begin{array}{rcl}
   I(a,b) &=&\displaystyle  
   \sum_{\beta=\alpha(a)
              -t(a)}^{\alpha(a)-t(a)+k}g(x_{\beta})\int_a^{x_{\alpha(a)+1}}dy\,w_{\beta,t}^{(k)}(y)\\
   \\
   &+&\displaystyle
       \sum_{\gamma=\alpha(a)+1}^{\alpha(b)}\sum_{\beta=\gamma-t_\gamma}^{\gamma-t_\gamma+k}g(x_{\beta})\int_{x_\gamma}^{x_{\gamma+1}}dy\,w_{\beta,t}^{(k)}(y)\\
   \\
   &-&\displaystyle \sum_{\beta=\alpha(b)-t(b)}^{\alpha(b)
       -t(b)+k}g(x_{\beta})\int_b^{x_{\alpha(b)+1}}dy\,w_{\beta,t}^{(k)}(y)\,,
   \end{array}
   :label: eq:integralT

with:

.. math:: t_\gamma = t(x_\gamma) = \mbox{max}\left[\mbox{min}\left[\gamma,\alpha_t -2\right] -\alpha_t+k + 1, 0\right] \,.

It turns out to be very complicated to write the expression in
Eq. :eq:`eq:integralT` in a more compact form, such as
that in Eq. :eq:`eq:integralT0`. Therefore, we
implement Eq. :eq:`eq:integralT` as it is. The one
thing left to do is to solve the integrals. This is done by using the
general formula:

.. math::
   \begin{array}{rcl}
     \displaystyle\int_{\ell_1}^{\ell_2}dy\,w_{\beta,t}^{(k)}(y)
     &=&\displaystyle \sum_{i=0}^{{\rm min}(k,\beta)}
                                                                         \theta(\ell_2-x_{\beta-i+t_\ell})\theta(x_{\beta-i+t_\ell+1}-\ell_1)\\
     \\
                                                                     &\times&\displaystyle \left[\prod^{k}_{m=0\atop m\ne
                                                                              i}\frac{1}{x_{\beta}-x_{\beta-i+m}}\right]
                                                                              \sum_{n=0}^k \frac{(-1)^n p_{\beta}^{(k)}(n)}{k-n+1} \left(\overline{\ell}_2^{k-n+1}-\overline{\ell}_1^{k-n+1}\right)\,,
   \end{array}

with :math:`t(\ell_1)=t(\ell_2)=t_\ell`, and
:math:`\overline{\ell}_1={\rm max}(\ell_1,x_{\beta-i+t_\ell})` and
:math:`\overline{\ell}_2={\rm min}(\ell_2,x_{\beta-i+t_\ell+1})`.

.. [1]
   Actually, it is not necessary to impose the constraint
   :math:`x_\alpha < x \leq x_{\alpha+k}`. In case this relation is not
   fulfilled one usually speaks about *extrapolation* rather than
   *interpolation*. If not necessary, this option is typically not
   convenient because it may lead to a substantial deterioration in the
   accuracy with which :math:`g(x)` is determined.

.. [2]
   From the definition of the set :math:`\{r_1,\dots,r_k\}`, it becomes
   clear that the subscript :math:`\beta` of the coefficients
   :math:`p_{\beta}^{(k)}` refers to the fact that this coefficients are
   computed on a vector of :math:`k` nodes where the node
   :math:`x_{\beta}` has been removed.
