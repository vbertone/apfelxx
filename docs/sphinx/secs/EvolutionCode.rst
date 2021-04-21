.. contents::
   :depth: 3
..

QCD Evolution Implementation
============================

The QCD DGLAP equation looks like this:

.. math::
   :label: dglap
   \mu^{2}\frac{\partial q_{i}(x,\mu^{2})}{\partial \mu^{2}}=\int^{1}_{x}\frac{dy}y P_{ij}\left(\frac{x}{y},\alpha_{s}(\mu^{2})\right)q_{j}(y,\mu^{2})

Now let us make the following definitions :math:`t=\ln(\mu^{2})`,
:math:`\tilde{q}(x,t)=xq(x,\mu^{2})`,
:math:`\tilde{P}_{ij}(x,t)=xP_{ij}(x,\alpha_{s}(\mu^{2}))` so that eq
:eq:`dglap` becomes:

.. math::
   :label: dglap2
   \frac{\partial \tilde{q}_{i}(x,t)}{\partial t}=\int^{1}_{x}\frac{dy}y \tilde{P}_{ij}\left(\frac{x}{y},t\right)\tilde{q}_{j}(y,t)

Interpolation
=============

In order to numerically solve the above equation, we need to write PDFs
as interpolated functions over an :math:`x` grid. In particular we want
to have something like this:

.. math:: \tilde{q}(y,t)=\sum^{N_{x}}_{\alpha=0}w_{\alpha}^{(k)}(y)\tilde{q}(x_{\alpha},t)\,,

where :math:`w_{\alpha}^{(k)}` are the interpolation functions of degree
:math:`k` we are looking for.

Using the Lagrange formula, we find that:

.. math::
   :label: LagrangeFormula
   w_{\alpha}^{(k)}(x) = \sum_{j=0,j \leq \alpha}^{k}\theta(x-x_{\alpha-j})\theta(x_{\alpha-j+1}-x)\prod^{k}_{\delta=0,\delta\ne j}\left[\frac{x-x_{\alpha-j+\delta}}{x_{\alpha}-x_{\alpha-j+\delta}}\right]\,.

This automatically means that:

.. math::
   :label: nonzero
   w_{\alpha}^{(k)}(x) \neq 0 \quad\mbox{for}\quad x_{\alpha-k} < x < x_{\alpha+1}.

This way we have that eq. :eq:`dglap2` becomes:

.. math::
   :label: dglap3
   \frac{\partial \tilde{q}_{i}(x,t)}{\partial t}=\sum_{\alpha}\int^{1}_{x}\frac{dy}y \tilde{P}_{ij}\left(\frac{x}{y},t\right)w_{\alpha}^{(k)}(y)\tilde{q}_{j}(x_{\alpha},t)

more in particular, if :math:`x` is one of the grid points
:math:`x_\beta`, simplifying a bit the notation we have that:

.. math::
   :label: dglap4
   \frac{\partial \tilde{q}_{i}(x_\beta,t)}{\partial t}=\sum_{\alpha}\underbrace{\left[\int^{1}_{x_\beta}\frac{dy}y \tilde{P}_{ij}\left(\frac{x_\beta}{y},t\right)w_{\alpha}^{(k)}(y)\right]}_{\Pi_{ij,\beta\alpha}(t)}\tilde{q}_{j}(x_{\alpha},t)\,.

Given eq. :eq:`nonzero`, it follows the condition:

.. math::
   :label: nonzero2
   \Pi_{ij,\beta\alpha}(t) \neq 0 \quad\mbox{for}\quad \beta \leq \alpha.

In addition, the integral in eq. :eq:`dglap4` which gives
:math:`\Pi_{ij,\beta\alpha}` can be optimized again using
eq. :eq:`nonzero` and it can be written as:

.. math::
   :label: optimization
   \Pi_{ij,\beta\alpha}(t) = \int^{b}_{a}\frac{dy}y \tilde{P}_{ij}\left(\frac{x_\beta}{y},t\right)w_{\alpha}^{(k)}(y)

where:

.. math:: a =  \mbox{max}(x_\beta,x_{\alpha-k})\quad\mbox{and}\quad b = \mbox{min}(1,x_{\alpha+1})\,.

However, performing a change of variables, :math:`\Pi_{ij,\beta\alpha}`
can also be written as:

.. math::
   :label: optimization2
   \Pi_{ij,\beta\alpha}(t) = \int^{d}_{c}\frac{dy}y \tilde{P}_{ij}(y,t)w_{\alpha}\left(\frac{x_\beta}{y}\right)

where this time:

.. math::
   :label: bounds2
   c =  \mbox{max}(x_\beta,x_\beta/x_{\alpha+1}) \quad\mbox{and}\quad d = \mbox{min}(1,x_\beta/x_{\alpha-k}) \,.

Verifying that eqs. :eq:`optimization`
and :eq:`optimization2` give the same numerical
result provides a cross-check of the correctness of the procedure.

Now, if rather than eq. :eq:`LagrangeFormula`, one
uses a logarithmic interpolation of the form:

.. math::
   :label: LagrangeFormulaLog
   w_{\alpha}^{(k)}(x) = \sum_{j=0,j \leq \alpha}^{k}\theta(x-x_{\alpha-j})\theta(x_{\alpha-j+1}-x)\prod^{k}_{\delta=0,\delta\ne j}\left[\frac{\ln(x)-\ln(x_{\alpha-j+\delta})}{\ln(x_{\alpha})-\ln(x_{\alpha-j+\delta})}\right]

over a logarithmically distributed grid, i.e. such that
:math:`\ln(x_{\beta})-\ln(x_{\alpha})=(\beta-\alpha)\delta x`, where the
step :math:`\delta x` is a constant, one has that:

.. math:: w_{\alpha}^{(k)}(x) = \sum_{j=0,j \leq \alpha}^{k}\theta(x-x_{\alpha-j})\theta(x_{\alpha-j+1}-x)\prod^{k}_{\delta=0,\delta\ne j}\left[\frac{1}{\delta x} \ln\left(\frac{x}{x_\alpha}\right)\frac{1}{j-\delta}+1\right]

that in general means that
:math:`w_\alpha^{(k)}(x)\equiv w_\alpha^{(k)}[\ln(x)-\ln(x_\alpha)]`.
Therefore in eq. :eq:`optimization2` we have that:

.. math::
   \begin{array}{rcl}
   \displaystyle w_{\alpha}^{(k)}\left(\frac{x_\beta}{y}\right) &\equiv& w_\alpha^{(k)}[\ln(x_{\beta})-\ln(x_\alpha)-\ln(y)]\\
   \\
                                         &=&w_\alpha^{(k)}[(\beta-\alpha)\delta x-\ln(y)]
   \end{array}

which means that :math:`w_{\alpha}^{(k)}(x_\beta/y)` only depends on the
difference :math:`(\beta-\alpha)` with the consequence that also
:math:`\Pi_{ij,\beta\alpha}` only depends on :math:`(\beta-\alpha)`.
Now, one can use this information with eq. :eq:`nonzero2`
to represent :math:`\Pi_{ij,\beta\alpha}(t)` as a matrix, where
:math:`\beta` is the row index and :math:`\alpha` the column index. In
this way :math:`\Pi_{ij,\beta\alpha}(t)` would look like this:

.. math::
   :label: MatrixRep
   \displaystyle \Pi_{ij}(t) = 
   \begin{pmatrix}
   a_0 &  a_1 & a_2 & \cdots & a_{N_x} \\
    0  & a_0 & a_1 & \cdots & a_{N_x-1} \\
    0  & 0   &  a_0 & \cdots & a_{N_x-2} \\
   \vdots & \vdots & \vdots & \ddots & \vdots \\
    0  &   0  &   0 & \cdots & a_0 
   \end{pmatrix}

therefore, if one knows the first raw of the matrix above, i.e.
:math:`\Pi_{ij,0\alpha}(t)`, it is possible to reconstruct the whole
matrix. Of course, this feature must be numerically verified but it
makes the computation of the evolution operators much faster because it
reduces the number of integrals to be computed by a factor :math:`N_x`.

Splitting Functions Treatment
=============================

The general form DGLAP splitting functions is the following:

.. math:: \tilde{P}_{ij}(x,t) = xP_{ij}^{R}(x,t) + \frac{xP_{ij}^{S}(t)}{(1-x)_+} + P_{ij}^{L}(t)x\delta(1-x)

where :math:`P_{ij}^{R}(x,t)` is the regular term that can be integrated
without any problem over any range, :math:`P_{ij}^{S}(x,t)` is instead
the function that multiplies the singular term which is regularized by
means of the plus prescription whose definition, referring to
eq. :eq:`optimization2`, is:

.. math::
   \begin{array}{c}
   \displaystyle \int_c^d dy \frac{f(y)}{(1-y)_+} = \int_c^d dy \frac{f(y) - f(1)\theta(d-1)}{1-y} - f(1)\theta(d-1)\int_0^c \frac{dy}{1-y} = \\
   \\
   \displaystyle \int_c^d dy \frac{f(y) - f(1)\theta(d-1)}{1-y} + f(1)\ln(1-c)\theta(d-1)\,.
   \end{array}

Finally :math:`P_{ij}^{L}(t)` is the coefficient of the local term, i.e.
the term proportional to :math:`\delta(1-x)`. Each of these terms has a
perturbative expansion that at N\ :math:`^k`\ LO looks like this:

.. math:: P_{ij}^{J}(x,t) = \sum_{n=0}^{k}a_s^{n+1}(t)P_{ij}^{J,(n)}(x)\qquad\mbox{with}\qquad J=R,S,L

Therefore one has that:

.. math::
   \begin{array}{c}
   \displaystyle \Pi_{ij,\beta\alpha}(t) = \\
   \\
   \displaystyle \sum_{n=0}^{k} a_s^{n+1}(t) \bigg\{\int^{d}_{c}dy\left[{P}_{ij}^{R,(n)}(y)w_{\alpha}\left(\frac{x_\beta}{y}\right)+\frac{{P}_{ij}^{S,(n)}}{1-y}\left(w_{\alpha}\left(\frac{x_\beta}{y}\right)-w_{\alpha}^{(k)}(x_\beta)\theta(d-1)\right)\right]\\
   \\
   \displaystyle +\left[{P}_{ij}^{S,(n)}\ln(1-c)\theta(d-1)+{P}_{ij}^{L,(n)}\right]w_{\alpha}^{(k)}(x_\beta)\bigg\}\,.
   \end{array}

Moreover it is easy to see that
:math:`w_{\alpha}^{(k)}(x_\beta)=\delta_{\beta\alpha}`, so that:

.. math::
   \begin{array}{c}
   \displaystyle \Pi_{ij,\beta\alpha}(t) = \\
   \\
   \displaystyle \sum_{n=0}^{k} a_s^{n+1}(t) \bigg\{\int^{d}_{c}dy\left[{P}_{ij}^{R,(n)}(y)w_{\alpha}\left(\frac{x_\beta}{y}\right)+\frac{{P}_{ij}^{S,(n)}}{1-y}\left(w_{\alpha}\left(\frac{x_\beta}{y}\right)-\delta_{\beta\alpha}\theta(d-1)\right)\right]\\
   \\
   \displaystyle +\left[{P}_{ij}^{S,(n)}\ln(1-c)\theta(d-1)+{P}_{ij}^{L,(n)}\right]\delta_{\beta\alpha}\bigg\}\,.
   \end{array}

Calling:

.. math::
   :label: pertExp
   \begin{array}{c}
   \displaystyle \Pi_{ij,\beta\alpha}^{(n)} = \\
   \\
   \displaystyle \int^{d}_{c}dy\left[{P}_{ij}^{R,(n)}(y)w_{\alpha}\left(\frac{x_\beta}{y}\right)+\frac{{P}_{ij}^{S,(n)}}{1-y}\left(w_{\alpha}\left(\frac{x_\beta}{y}\right)-\delta_{\beta\alpha}\theta(d-1)\right)\right]\\
   \\
   \displaystyle +\left[{P}_{ij}^{S,(n)}\ln(1-c)\theta(d-1)+{P}_{ij}^{L,(n)}\right]\delta_{\beta\alpha}\,,
   \end{array}

\ we have that:

.. math::
   :label: splittingexp
   \Pi_{ij,\beta\alpha}(t) = \sum_{n=0}^{k} a_s^{n+1}(t) \Pi_{ij,\beta\alpha}^{(n)}\,,

and the integrals :math:`\Pi_{ij,\beta\alpha}^{(n)}` do not depend on
the energy therefore, once the grid (and the number of active flavours)
has been fixed, they can be evaluate once and for all at the beginning
and used for the evolution to any scale.

It is not very easy to see that eq. :eq:`pertExp` respects
the symmetry described in eq. :eq:`MatrixRep`. To show
this, we distinguish two case: 1) :math:`d < 1` and 2) :math:`d = 1`. In
the case 1), due to the presence of the Heaviside’s functions
:math:`\theta(d-1)`, eq. :eq:`pertExp` reduces to:

.. math::
   :label: pertExp2
   \Pi_{ij,\beta\alpha}^{(n)} = \int^{d}_{c}dy\left[{P}_{ij}^{R,(n)}(y)+\frac{{P}_{ij}^{S,(n)}}{1-y}\right]w_{\alpha}\left(\frac{x_\beta}{y}\right) + {P}_{ij}^{L,(n)}\delta_{\beta\alpha}\,,

which evidently obeys eq. :eq:`MatrixRep`. In the case
2), instead, we have:

.. math::
   :label: pertExp3
   \begin{array}{c}
   \displaystyle \Pi_{ij,\beta\alpha}^{(n)} = \int^{1}_{c}dy\left[{P}_{ij}^{R,(n)}(y)w_{\alpha}\left(\frac{x_\beta}{y}\right)+\frac{{P}_{ij}^{S,(n)}}{1-y}\left(w_{\alpha}\left(\frac{x_\beta}{y}\right)-\delta_{\beta\alpha}\right)\right]\\
   \\
   \displaystyle +\left[{P}_{ij}^{S,(n)}\ln(1-c)+{P}_{ij}^{L,(n)}\right]\delta_{\beta\alpha}\,,
   \end{array}

and apparently, if :math:`\alpha=\beta`, the term :math:`\ln(1-c)` seems
to break the symmetry. However, this is not the case. In fact, from
eq. :eq:`bounds2`, we know that in this particular case:

.. math:: c = \mbox{max}(x_\beta,x_\beta/x_{\beta+1}) = \frac{x_\beta}{x_{\beta+1}}

because :math:`x_{\beta+1}<1`. In addition, on a logarithmically
distributed grid, :math:`x_{\beta+1}=x_{\beta}\exp(\delta x)`, where
:math:`\delta x` is the constant step. Therefore, it turns out that:

.. math:: \ln(1-c) = \ln\left(1-\frac{x_{\beta}}{x_{\beta+1}}\right) = \ln[1 - \exp(-\delta x)]\,,

that is a constant which does not depend on the indices :math:`\alpha`
and :math:`\beta` and therefore does not break the symmetry given in
eq. :eq:`MatrixRep`.

Solution of the DGLAP Equation
==============================

As a consequence of the DGLAP equation form, we can assume that
:math:`\tilde{q}_{i}(x_\beta,t)\equiv q_{i,\beta}(t)` evolves between
the energies :math:`t` and :math:`t_0` according to the following
(discretized) evolution equation:

.. math:: \tilde{q}_{i,\beta}(t) = \sum_{\gamma,k} M_{ik,\beta\gamma}(t,t_0)\tilde{q}_{k,\gamma}(t_0)

with the boundary condition
:math:`M_{ik,\beta\gamma}(t_0,t_0)=\delta_{ik}\delta_{\beta\gamma}`. It
follows that eq. :eq:`dglap4` takes the form:

.. math::
   :label: tosolve
   \left\{\begin{array}{l}
   \displaystyle \frac{\partial  M_{ij,\alpha\beta}(t,t_0)}{\partial t}=\sum_{\gamma,k} \Pi_{ik,\alpha\gamma}(t)M_{kj,\gamma\beta}(t,t_0)\\
   \\
   \displaystyle M_{ij,\alpha\beta}(t_0,t_0)=\delta_{ij}\delta_{\alpha\beta}
   \end{array}\right.

which is a first order linear differential equation in the quantity
:math:`M_{ij,\alpha\beta}(t,t_0)` that, as we actually do, can be
numerically solved using the fourth order Adaptive Step-size Control
Runge-Kutta algorithm. Using the arguments we discussed above, we do not
need to compute all the entries of :math:`\Pi_{ik,\gamma\alpha}(t)`. In
addition, as we have already shown, the perturbative contributions to
:math:`\Pi_{ik,\gamma\alpha}(t)` can be precomputed before solving the
differential equation in eq. :eq:`tosolve`.

.. container::
   :name: the-non-singlet

   .. rubric:: The Non Singlet
      :name: the-non-singlet

The non-singlet case is the easiest one because the differential
equations in eq. :eq:`tosolve` decouple in the flavour pair
:math:`(i,j)`, and we can write them as:

.. math::
   :label: tosolveNS
   \left\{\begin{array}{l}
   \displaystyle \frac{\partial  \mathcal{M}_{\alpha\beta}^{(i)}(t,t_0)}{\partial t}=\sum_{\gamma=0}^{N_x} \mathcal{P}_{\alpha\gamma}^{(i)}(t)\mathcal{M}_{\gamma\beta}^{(i)}(t,t_0)\\
   \\
   \displaystyle \mathcal{M}_{\alpha\beta}^{(i)}(t_0,t_0)=\delta_{\alpha\beta}
   \end{array}\right.\quad\mbox{with } i=+,-,V\,.

\ As a further simplification, we can use the fact that at LO :math:`+`,
:math:`-` and :math:`V` all the evolution operators are equal, therefore
solving only one of them is enough. at NLO instead only :math:`-` and
:math:`V` are equal while at NNLO they are all different.

Now, given the symmetries carried by :math:`\Pi_{ij,\alpha\beta}`, we
can write:

.. math:: \mathcal{P}_{\alpha\gamma}^{(i)} = \mathcal{P}_{0(\gamma-\alpha)}^{(i)}\theta(\gamma-\alpha)\,,

so that eq. :eq:`tosolveNS` becomes:

.. math::
   :label: tosolveNS1
   \left\{\begin{array}{l}
   \displaystyle \frac{\partial  \mathcal{M}_{\alpha\beta}^{(i)}(t,t_0)}{\partial t}= \sum_{\delta=0}^{N_x-\alpha} \mathcal{P}_{0\delta}^{(i)}(t)\mathcal{M}_{(\alpha+\delta)\beta}^{(i)}(t,t_0)\\
   \\
   \displaystyle \mathcal{M}_{\alpha\beta}^{(i)}(t_0,t_0)=\delta_{\alpha\beta}
   \end{array}\right.\quad\mbox{with } i=+,-,V\,.

.. container::
   :name: the-singlet

   .. rubric:: The Singlet
      :name: the-singlet

The singlet sector is totally analogous to the non-singlet one, the only
difference is that there is one additional summation over the flavours.
In practice we have that:

.. math::
   :label: tosolveSG
   \left\{\begin{array}{l}
   \displaystyle \frac{\partial  \mathcal{M}_{ij,\alpha\beta}^{\rm SG}(t,t_0)}{\partial t}= \sum_k \sum_{\delta=0}^{N_x-\alpha} \mathcal{P}_{ik,0\delta}^{\rm SG}(t)\mathcal{M}_{kj,(\alpha+\delta)\beta}^{\rm SG}(t,t_0)\\
   \\
   \displaystyle \mathcal{M}_{ij,\alpha\beta}^{\rm SG}(t_0,t_0)=\delta_{ij}\delta_{\alpha\beta}
   \end{array}\right.\quad\mbox{with }i,j,k=q,g\,.

Alternative Solutions of the DGLAP Equation
===========================================

In this section we show how it is possible to solve the DGLAP equation
in an alternative way with respect to that shown in the previous
sections exploiting the RGE of the running coupling :math:`\alpha_s`.
This will lead to a different equation that admits two preturbatively
equivalent solutions: the first, that we will refer to as "exact"
solution and that reproduces the solution seen above, and the second,
the so-called "expanded" solution, that reproduces the solution usually
adopted in the N-space code (like NNPDF).

The starting point is the RGE:

.. math::
   :label: RGEalpha
   \mu^2\frac{\partial a_s}{\partial \mu^2} = \frac{\partial
     a_s}{\partial t} = \beta(a_s)\,,

\ where:

.. math:: a_s \equiv \frac{\alpha_s}{4\pi}

and:

.. math:: \beta(a_s) = -a_s^2 \sum_{n=0}^{N} a_s^n \beta_n \,.

where :math:`N` represents the desired preturbative order. Using eq.
:eq:`RGEalpha` and eq.
:eq:`splittingexp`, we can rewrite eq.
:eq:`tosolve` as:

.. math::
   :label: tosolvealpha
   \left\{\begin{array}{l}
   \displaystyle \frac{\partial  M_{ij,\alpha\beta}(t,t_0)}{\partial
     a_s}= -\frac{1}{a_s}\sum_{\gamma,k} \left[\frac{\displaystyle \sum_{n=0}^N a_s^n
       \Pi_{ik,\alpha\gamma}^{(n)}}{ \displaystyle \sum_{n=0}^{N} a_s^n \beta_n}\right]M_{kj,\gamma\beta}(t,t_0)\\
   \\
   \displaystyle M_{ij,\alpha\beta}(t_0,t_0)=\delta_{ij}\delta_{\alpha\beta}
   \end{array}\right.

Now there are two possibile way to solve eq.
:eq:`tosolvealpha`: either we solve directly it
numerically as it is or we first expand the term in the square brackets
as a series of :math:`a_s` keeping only the terms up to order
:math:`a_s^N` and then we solve the equation. It is obvious that the
first way to solve eq. :eq:`tosolvealpha` must be
numerical equal to the solution of eq. :eq:`tosolve`. The
second way instead is not numerically equal but is perturbatively
equivalent. This second solution is referred to as N-space solution as
it is usually used in the N-space approach because it permits to solve
analytically the DGLAP equation.

To expand the term in the square brackets, we notice that up to NNLO
(:math:`N=2`) we have that:

.. math::
   \frac{1}{\displaystyle \sum_{n=0}^{2} a_s^n \beta_n} =
   \frac1{\beta_0}\left[ 1 - \frac{\beta_1}{\beta_0} a_s + \left(\frac{\beta_1^2}{\beta_0^2}- \frac{\beta_2}{\beta_0}\right)a_s^2\right] + \mathcal{O}(a_s^{3})\,,

so that:

.. math::
   \begin{array}{rcl}
   \displaystyle\frac{\displaystyle \sum_{n=0}^2 a_s^n
       \Pi_{ik,\alpha\gamma}^{(n)}}{ \displaystyle \sum_{n=0}^{2} a_s^n
       \beta_n} &=&\displaystyle \frac{1}{\beta_0}\left\{\Pi_{ik,\alpha\gamma}^{(0)} +
   a_s\left[\Pi_{ik,\alpha\gamma}^{(1)} - b_1 \Pi_{ik,\alpha\gamma}^{(0)}\right]\right.\\
   \\
   &+& \displaystyle
   \left. a_s^2\left[\Pi_{ik,\alpha\gamma}^{(2)} - b_1
   \Pi_{ik,\alpha\gamma}^{(1)} +\left(b_1^2-b_2\right)
   \Pi_{ik,\alpha\gamma}^{(0)} \right]\right\}  + \mathcal{O}(a_s^{3})
   \end{array}

where we have defined:

.. math:: b_n \equiv \frac{\beta_n}{\beta_0}\,.

Finally, defining:

.. math::
   \begin{array}{l}
   \displaystyle\widetilde{\Pi}_{ik,\alpha\gamma}^{(0)} \equiv
   \Pi_{ik,\alpha\gamma}^{(0)} \,,\\
   \displaystyle\widetilde{\Pi}_{ik,\alpha\gamma}^{(1)} \equiv
   \Pi_{ik,\alpha\gamma}^{(1)} - b_1 \Pi_{ik,\alpha\gamma}^{(0)} \,,\\
   \displaystyle\widetilde{\Pi}_{ik,\alpha\gamma}^{(2)} \equiv \Pi_{ik,\alpha\gamma}^{(2)} - b_1
   \Pi_{ik,\alpha\gamma}^{(1)} +\left(b_1^2-b_2\right)
   \Pi_{ik,\alpha\gamma}^{(0)}\,,
   \end{array}

we can write eq. :eq:`tosolvealpha` up to NNLO as:

.. math::
   :label: tosolvealphatrn
   \left\{\begin{array}{l}
   \displaystyle \frac{\partial  M_{ij,\alpha\beta}(t,t_0)}{\partial
     a_s}= -\frac{1}{a_s\beta_0}\sum_{\gamma,k} \left[\displaystyle \sum_{n=0}^2 a_s^n
       \widetilde{\Pi}_{ik,\alpha\gamma}^{(n)}\right]M_{kj,\gamma\beta}(t,t_0)\\
   \\
   \displaystyle M_{ij,\alpha\beta}(t_0,t_0)=\delta_{ij}\delta_{\alpha\beta}
   \end{array}\right.\,.

Solving eq. :eq:`tosolvealphatrn` provides the
so-called "expanded" solution.

Small-:math:`x` Resummation using ``Hell``
==========================================

The implemenatation of the small-:math:`x` resummation in ``APFEL`` is
done by interfacing it to the code ``Hell`` by Marco Bonvini. ``Hell``
provides, amongst other things, the small-:math:`x` resummed singlet
splitting functions (times :math:`x`) up to NLL to be matched to the
unresummed LO, NLO and NNLO splitting functions. In practice, the user
can specify the logarithmic accuracy (LL or NLL) and, for each
logarithmic accuracy, the perturbative accuracy to which the rusummed
splitting functions are to be matched.

With the inclusion of the small-:math:`x` resummation, the splitting
functions :math:`{P}_{ij}` in eq. :eq:`dglap` should be
interpreted as:

.. math::
   :label: ResummedSplittings
   {P}_{ij}(x,\alpha_s) = \underbrace{{P}_{ij}^{\rm
       FO}(x,\alpha_s)}_{\mbox{\tiny Already present in {\tt APFEL}}} +
   \underbrace{{P}_{ij}^{\rm Res-FO}(x,\alpha_s)}_{\mbox{\tiny Provided
       by {\tt Hell}}}

\ where "FO" stands for Fixed Order and "Res" for Resummed. The
essential difference between :math:`{P}_{ij}^{\rm FO}` and
:math:`{P}_{ij}^{\rm Res-FO}` stems from the fact that the former admits
the usual pertubative expansion:

.. math:: {P}_{ij}^{\rm FO}(x,\alpha_s) = \sum_{n=0}^N a_s^{n+1} P_{ij}^{(n)}(x)\,,

while the latter, by definition, does not. This feature of
:math:`{P}_{ij}^{\rm Res-FO}` forbids to precompute the perturbative
coefficients in the r.h.s. of the DGLAP equation before solving it. In
principle then, one should recompute the integrals of the splitting
functions on the :math:`x`-space interpolation grid at every step of the
algorithm that numerically solves the discretized DGLAP equation. This
is clearly very unefficient and would enormously inflate the computation
time. In order to keep the computation time under control, we use an
interpolation grid also in :math:`\alpha_s`. In practice, we precompute
the integrals of the resummed splitting function over the
:math:`x`-space interpolation grid for several values of
:math:`\alpha_s` (logarithmically?) distributed over a reasonable range,
so that the values of the same integrals for any value of
:math:`\alpha_s` needed during the numerical solution of the DGLAP
equation would be obtained by interpolation.

Using the notation of eq. :eq:`pertExp`, we have that the
integral of the resummed splitting functions on the :math:`x`-space grid
would take the form:

.. math:: \Pi_{ij,\beta\alpha}^{\rm Res}(\alpha_s) = \int^{d}_{c}dy{P}_{ij}^{\rm Res-FO}(y,\alpha_s)w_{\alpha}\left(\frac{x_\beta}{y}\right)\,,

so that eq. :eq:`tosolve` would become:

.. math::
   :label: tosolveRes
   \left\{\begin{array}{l}
   \displaystyle \frac{\partial  M_{ij,\alpha\beta}(t,t_0)}{\partial t}=\sum_{\gamma,k} \left[\Pi_{ik,\alpha\gamma}(t)+\Pi_{ij,\beta\alpha}^{\rm Res}(\alpha_s(t))\right]M_{kj,\gamma\beta}(t,t_0)\\
   \\
   \displaystyle M_{ij,\alpha\beta}(t_0,t_0)=\delta_{ij}\delta_{\alpha\beta}
   \end{array}\right.\,.

As we mentioned above, :math:`\Pi_{ij,\beta\alpha}^{\rm Res}` cannot be
expandend as a truncated series of :math:`\alpha_s` and thus we cannot
precompute the perturbative coefficients making the numerical solution
of the DGLAP equation efficient. A possible way out is to precompute
:math:`\Pi_{ij,\beta\alpha}^{\rm Res}` over a grid in :math:`\alpha_s`,
say :math:`\alpha_s^{(\tau)}` with :math:`\tau = 0,\dots,m`, so that,
after the initialization step, we have the set of integrals:

.. math::
   \Pi_{ij,\beta\alpha,\tau}^{\rm Res} = \int^{d}_{c}dy{P}_{ij}^{\rm
     Res-FO}(y,\alpha_s^{(\tau)})w_{\alpha}\left(\frac{x_\beta}{y}\right)\quad
   \tau = 0,\dots,m\,,

\ and for obtaining the value of :math:`\Pi_{ij,\beta\alpha}^{\rm Res}`
for a generic value of :math:`\alpha_s` we use the linear interpolation.
Supposing that
:math:`\alpha_s^{(\tau)}\leq\alpha_s<\alpha_s^{(\tau+1)}`, we have that:

.. math::
   \Pi_{ij,\beta\alpha}^{\rm Res}(\alpha_s) =
   \left(\frac{\alpha_s^{(\tau+1)}-\alpha_s}{\alpha_s^{(\tau+1)}-\alpha_s^{(\tau)}}\right)\Pi_{ij,\beta\alpha,\tau}^{\rm Res} +\left(\frac{\alpha_s-\alpha_s^{(\tau)}}{\alpha_s^{(\tau+1)}-\alpha_s^{(\tau)}}\right)\Pi_{ij,\beta\alpha,\tau+1}^{\rm Res}

There is a further complication that we need to deal with that is the
number of active flavours. In fact, any integral must be computed with
the correct number :math:`n_f` of active flavours. Assuming to be
working in the VFNS (the FFNS is instead trivial), what we need then is
a grid in :math:`\alpha_s` that has a node for each crossing point, that
is there must be one value of the index :math:`\tau` such that
:math:`\alpha_s^{(\tau)}=\alpha_s(m_h)` where :math:`m_h` is the mass of
any heavy flavour and such that for :math:`\alpha_s<\alpha_s^{(\tau)}`
there are, say, :math:`n_f` active flavours, while for
:math:`\alpha_s\geq\alpha_s^{(\tau)}` there are :math:`n_f+1` active
flavours. This in practice means that the grid in :math:`\alpha_s` has
as many fixed points as potentially active flavours. Unfortunately this
is not enough because, when considering the NNLO evolution in the VFNS,
the evolution of :math:`\alpha_s`, as well as that of PDFs, has a
discontinuity in correspondence of the heavy quark thresholds( [1]_). To
overcome this problem we may assign to the point of the grid
corresponding to the heavy threshold :math:`m_h` two value, :math:`i.e.`
the values :math:`\alpha_s^{(\tau)}=\alpha_s(m_h-\varepsilon)` and
:math:`\alpha_s^{(\tau+1)}=\alpha_s(m_h)`. This trick, in conjuction
with the linear interpolation, does to job without any further
assumption.

We will now try to derive the form of the “expanded” solution in the
presence of small-:math:`x` resummation. In order to do so, we need to
recognise that the resummed spliiting functions in
eq. :eq:`ResummedSplittings` admit an expansion
in :math:`\alpha_s` for fixed

A Remark on the Interpolation Functions
=======================================

Just for the record, it is useful to derive the expression for the
interpolation functions given in
eq. :eq:`LagrangeFormula` and show how this is not
the only possible choice.

Suppose we want to perform an interpolation of degree :math:`k` of the
test function :math:`g` in the point :math:`x`. As is well known, we
will need a subset of :math:`k+1` consecutive points on the total
interpolation grids :math:`\{x_{\alpha},\dots,x_{\alpha+k}\}`. However,
the relative position between the point :math:`x` and the subset of
points used for the interpolation is arbitrary. In principle, it is not
even required that :math:`x` is somewhere between :math:`x_\alpha` and
:math:`x_{\alpha+k}`. However, in this case one would talk about
extrapolation rather than interpolation and this is clearly not a
convenient option because it would lead to a substantial deteriotation
in the accuracy with which :math:`g(x)` is determined. As a consequence,
it is convenient to choose the subset of points in such a way that
:math:`x_\alpha < x \leq x_{\alpha+k}`. However, the ambiguity remains
because there are :math:`k` possible choices of the subset of points
accordind to which :math:`x_\alpha < x \leq x_{\alpha+1}`, or
:math:`x_{\alpha+1} < x \leq x_{\alpha+2}`, and so on.

In particular, to derive eq. :eq:`LagrangeFormula`
we have assumed that:

.. math::
   :label: IntAssumption1
   x_\alpha < x \leq x_{\alpha+1}\,.

Let’s see how eq. :eq:`LagrangeFormula` comes out.
Using the standard Lagrange interpolation procedure, we can approximate
the function :math:`g` in :math:`x` as:

.. math::
   :label: particularCase
   g(x) = \sum_{i=0}^k\ell_i^{(k)}(x)g(x_{\alpha+i})

\ where :math:`\ell_i^{(k)}` is the :math:`i`-th Lagrange polynomial of
degree :math:`k` which can be written as:

.. math:: \ell_i^{(k)}(x) = \prod^{k}_{m=0,m\ne i}\frac{x-x_{\alpha+m}}{x_{\alpha+i}-x_{\alpha+m}}\,.

However, as we said, we impose that
eq. :eq:`particularCase` applies only for the
assumption in eq. :eq:`IntAssumption1` is fulfilled.
We can then generalize it by writing:

.. math::
   :label: particularCaseTheta
   g(x) = \theta(x-x_{\alpha})\theta(x_{\alpha+1}-x)\sum_{i=0}^k g(x_{\alpha+i})\prod^{k}_{m=0,m\ne i}\frac{x-x_{\alpha+m}}{x_{\alpha+i}-x_{\alpha+m}}\,.

Now, if we want to relax the restriction in
eq. :eq:`IntAssumption1`, we just have to sum over
all nodes of the global interpolation grid, that is:

.. math::
   :label: generalCase
   g(x) = \sum_{\alpha=0}^{N_x}\theta(x-x_{\alpha})\theta(x_{\alpha+1}-x)\sum_{i=0}^k g(x_{\alpha+i})\prod^{k}_{m=0,m\ne i}\frac{x-x_{\alpha+m}}{x_{\alpha+i}-x_{\alpha+m}}\,.

Defining :math:`\beta=\alpha+i`, we can rewrite the equation above as:

.. math::
   :label: generalCase2
   g(x) = \sum_{\beta=0}^{N_x}g(x_{\beta}) \sum_{i=0,i\leq\beta}^k \theta(x-x_{\beta-i})\theta(x_{\beta-i+1}-x) \prod^{k}_{m=0,m\ne i}\frac{x-x_{\beta-i+m}}{x_{\beta}-x_{\beta-i+m}}\,,

where the additional condition :math:`i\leq\beta` comes from the
condition :math:`\alpha\geq 0`. Eq. :eq:`generalCase2`
is clearly equivalent to eq. :eq:`LagrangeFormula`,
assuming that:

.. math:: w_\beta^{(k)}(x) = \sum_{i=0,i\leq\beta}^k \theta(x-x_{\beta-i})\theta(x_{\beta-i+1}-x) \prod^{k}_{m=0,m\ne i}\frac{x-x_{\beta-i+m}}{x_{\beta}-x_{\beta-i+m}}\,.

Now, instead of starting from the assumption in
eq. :eq:`IntAssumption1`, we start from the more
general condition:

.. math::
   :label: IntAssumptionGen
   x_{\alpha+t} < x \leq x_{\alpha+t+1}\quad\mbox{with}\quad t = 0,\dots,k-1\,,

we have that the interpolation formula would look like this:

.. math::
   :label: MoreGeneralCase
   g(x) = \sum_{\alpha=0}^{N_x}\theta(x-x_{\alpha+t})\theta(x_{\alpha+t+1}-x)\sum_{i=0}^k g(x_{\alpha+i})\prod^{k}_{m=0,m\ne i}\frac{x-x_{\alpha+m}}{x_{\alpha+i}-x_{\alpha+m}}\,,

that can be rearranged as:

.. math:: g(x) = \sum_{\beta=0}^{N_x}g(x_{\beta}) \sum_{i=0,i\leq\beta}^k \theta(x-x_{\beta-i+t})\theta(x_{\beta-i+t+1}-x) \prod^{k}_{m=0,m\ne i}\frac{x-x_{\beta-i+m}}{x_{\beta}-x_{\beta-i+m}}\,.

Therefore the “generalized” interpolation functions are:

.. math:: w_{\beta,t}^{(k)}(x) = \sum_{i=0,i\leq\beta}^k \theta(x-x_{\beta-i+t})\theta(x_{\beta-i+t+1}-x) \prod^{k}_{m=0,m\ne i}\frac{x-x_{\beta-i+m}}{x_{\beta}-x_{\beta-i+m}}\,,

and they assume that:

.. math:: x_{\alpha+t} < x \leq x_{\alpha+t+1}\quad\mbox{with}\quad t = 0,\dots,k-1\,.

.. [1]
   A discontinuity appears also at NLO if factorization and
   renormalization scales are not equal.
