==========
Statistics
==========

.. contents::
   :depth: 3
..

``APFEL++`` is often used to fit non-perturbative functions, such as
parton distribution functions (PDFs) and fragmentation functions (FFs),
to experimental data. As a consequence, it is necessary to make full use
of the appropriate statistical tools that allow for a correct treatment
of the experimental information when comparing predictions to data. This
section is devoted to a detailed discussion of the concepts and tools
that are often used when fitting PDFs and FFs to data.

Most of these tools are not directly implemented in ``APFEL++`` but can
be found in other codes, as for instance ``NangaParbat`` (see
https://github.com/vbertone/NangaParbat).

The :math:`\chi^2` in the presence of correlations
==================================================

Suppose to have an ensamble of :math:`n` measurements having the
following structure:

.. math::
   m_i\pm \sigma_{i,\rm stat} \pm \sigma_{i,\rm unc} \pm \sigma_{i,\rm
       corr}^{(1)}\pm\dots \pm \sigma_{i,\rm
       corr}^{(k)}\,,

\ where :math:`m_i`, with :math:`i=1,\dots, n`, is the central value of
the :math:`i`-th measurement, :math:`\sigma_{i,\rm stat}` its
(uncorrelated) statistical uncertainty, :math:`\sigma_{i,\rm unc}` its
uncorrelated systematic uncertainty [1]_, and
:math:`\sigma_{i,\rm corr}^{(l)}`, with :math:`l=1,\dots,k`, its
correlated systematic uncertainties. With this information at hand, one
can construct the full covariance matrix :math:`V_{ij}` as follows (see
for example Ref. (Ball and others 2013)):

.. math::
   V_{ij}=\left(\sigma_{i,\rm stat}^2 +\sigma_{i,\rm unc}^2\right)\delta_{ij} + \sum_{l=1}^{k}\sigma_{i,\rm
       corr}^{(l)}\sigma_{j,\rm
       corr}^{(l)}\,.
   :label: eq:covmat

This is a clearly symmetric matrix. Given a set of predictions
:math:`t_i` corresponding to the :math:`n` measurements of the ensamble,
the :math:`\chi^2` takes the form:

.. math::
   \chi^2=
     \sum_{i,j=1}^{n}\left(m_i-t_i\right)V_{ij}^{-1}\left(m_j-t_j\right) =
     \mathbf{y}^{T} \cdot \mathbf{V}^{-1} \cdot \mathbf{y}\,,
   :label: eq:chi2cov

where in the second equality I have used the matricial notation and
defined :math:`y_i = m_i-t_i`. A convenient way to compute the
:math:`\chi^2` relies on the Cholesky decomposition of the covariance
matrix :math:`\mathbf{V}`. In particular, it can be proven that any
symmetric and positive definite matrix :math:`\mathbf{V}` can be
decomposed as:

.. math::
   \mathbf{V} = \mathbf{L}\cdot\mathbf{L}^{T}\,,
   :label: eq:choleskydec

where :math:`\mathbf{L}` is a lower triangular matrix whose entries are
related recursively to those of :math:`\mathbf{V}` as follows:

.. math::
   \begin{array}{rcl}
     L_{kk} &=&\displaystyle \sqrt{V_{kk}-\sum_{j=1}^{k-1}L_{kj}^2}\,,\\
     \\
     L_{ik} &=&\displaystyle
                \frac{1}{L_{kk}}\left(V_{ik}-\sum_{j=1}^{k-1}L_{ij}L_{kj}\right)\,,\quad
                k < i\,,\\
   \\
     L_{ik} &=&\displaystyle 0\,,\quad
                k > i\,.\\
   \end{array}
   :label: eq:cholalg

It is then easy to see that the :math:`\chi^2` can be written as:

.. math:: \chi^2 = \left|\mathbf{L}^{-1}\cdot \mathbf{y}\right|^2\,.

But the vector :math:`{\vec \chi} \equiv \mathbf{L}^{-1}\cdot \mathbf{y}`
is the solution of the linear system:

.. math::
   \mathbf{L} \cdot {\vec \chi} = \mathbf{y}\,,
   :label: eq:chivdef

that can be efficiently solved by forward substitution, so that:

.. math:: \chi^2 = \left|{\vec \chi}\right|^2\,.

Following this procedure, one does not need to compute explicitly the
inverse of the covariance matrix :math:`\mathbf{V}`, simplifying
significantly the computation of the :math:`\chi^2`.

Additive and multiplicative uncertainties
=========================================

The correlated systematic uncertainties
:math:`\sigma_{i,\rm corr}^{(l)}` may be either *additive* or
*multiplicative*. The nature of the single uncertainties is typically
provided by the experiments that release the measurements. A typical
example of multiplicative uncertainty is the luminosity uncertainty but
there can be others.

Now let us express all the correlated systematic uncertainties
:math:`\sigma_{i,\rm corr}^{(l)}` as relative to the associate central
value :math:`m_i`, so that one defines [2]_:

.. math:: \sigma_{i,\rm corr}^{(l)}\equiv  \delta_{i,\rm corr}^{(l)} m_i

\ and let us also define
:math:`s_i^2\equiv \sigma_{i,\rm stat}^2 +\sigma_{i,\rm unc}^2` so that
Eq. :eq:`eq:covmat` can be rewritten as:

.. math::
   V_{ij}=s_i^2\delta_{ij} + \left(\sum_{l=1}^{k}\delta_{i,\rm
       corr}^{(l)}\delta_{j,\rm
       corr}^{(l)}\right)m_im_j\,.
   :label: eq:covmat2

Now I split the correlated systematic uncertainties into :math:`k_a`
additive uncertainties and :math:`k_m` multiplicative uncertainties,
such that :math:`k_a+k_m=k`. This way
Eq. :eq:`eq:covmat2` takes the form:

.. math::
   V_{ij}=s_i^2\delta_{ij} + \left(\sum_{l=1}^{k_a}\delta_{i,\rm
       add}^{(l)}\delta_{j,\rm
       add}^{(l)}+\sum_{l=1}^{k_m}\delta_{i,\rm
       mult}^{(l)}\delta_{j,\rm
       mult}^{(l)}\right)m_im_j\,.
   :label: eq:covmat3

It is believed that this definition of the covariance matrix is
problematic in that it results in the so-called D’Agostini bias of the
multiplicative uncertainties (D’Agostini 1994). A possible solution to
this problem is the so-called :math:`t_0`-prescription (Ball et al.
2010), where the experimental central value :math:`m_i` in the
multiplicative term is replaced by a fixed theoretical predictions
:math:`t_i^{(0)}`, typically computed in a previous fit in which the
“standard” definition of the covariance matrix in
Eq. :eq:`eq:covmat` (often referred to as *experimental*
definition) is used. Applying the :math:`t_0` prescription, the
covariance matrix takes the form:

.. math::
   V_{ij}=s_i^2\delta_{ij} + \sum_{l=1}^{k_a}\delta_{i,\rm
       add}^{(l)}\delta_{j,\rm
       add}^{(l)}m_im_j+\sum_{l=1}^{k_m}\delta_{i,\rm
       mult}^{(l)}\delta_{j,\rm
       mult}^{(l)}t_i^{(0)}t_j^{(0)}\,.
   :label: eq:covmat4

.. _eq:sysshifts:

Determining the systematic shifts
=================================

In order to visualise the effect of systematic uncertainties, it is
instructive to compute the *systematic shift* generated by the
systematic uncertainties. To do so, one needs to write the
:math:`\chi^2` in terms of the so-called “nuisance parameters”
:math:`\lambda_\alpha`. One can show that the definition of the
:math:`\chi^2` in Eq. :eq:`eq:chi2cov` is equivalent
to (Ball and others 2013):

.. math::
   \chi^2 = \sum_{i=1}^n\frac{1}{s_i^2}\left(m_i -t_i
     -\sum_{\alpha=1}^k\lambda_\alpha \sigma_{i,\rm corr}^{(\alpha)}
   \right)^2 + \sum_{\alpha=1}^k\lambda_\alpha^2\,.
   :label: eq:chi2nuis

\ The optimal value of the nuisance parameters can be computed by
minimising the :math:`\chi^2` with respect to them, that is imposing:

.. math:: \frac{\partial \chi^2}{\partial \lambda_\beta} = 0\,.

This yields the system:

.. math::
   \sum_{\beta=1}^kA_{\alpha\beta}\lambda_\beta =\rho_\alpha\,,
   :label: eq:nuissys

with:

.. math::
   A_{\alpha\beta}= \delta_{\alpha\beta}+\sum_{i=1}^n\frac{\sigma_{i,\rm
       corr}^{(\alpha)}\sigma_{i,\rm
       corr}^{(\beta)}}{s_i^2}\quad\mbox{and}\quad
   \rho_\alpha=\sum_{i=1}^n\frac{m_i-t_i}{s_i^2}\sigma_{i,\rm
     corr}^{(\alpha)}\,,
   :label: eq:sysing

that determines the values of :math:`\lambda_\beta`. The quantity:

.. math:: d_i =\sum_{\alpha=1}^k\lambda_\alpha \sigma_{i,\rm corr}^{(\alpha)}

in Eq. :eq:`eq:chi2nuis` can be interpreted as a shift
caused by the correlated systematic uncertainties. Defining the shifted
predictions as:

.. math:: \overline{t}_i =t_i+d_i\,,

the :math:`\chi^2` reads:

.. math::
   \chi^2 = \sum_{i=1}^n\left(\frac{m_i -\overline{t}_i}{s_i}\right)^2
     + \sum_{\alpha=1}^k\lambda_\alpha^2\,.
   :label: eq:chi2nuisshift

Therefore, up to a penalty term given by the sum of the square of the
nuisance parameters, the :math:`\chi^2` takes the form of the
uncorrelated definition. In order to achieve a visual assessment of the
agreement between data and theory, it appears natural to compare the
central experimental values :math:`m_i` to the shifted theoretical
predictions :math:`\overline{t}_i` in units of the uncorrelated
uncertainty :math:`s_i`.

Effect of cuts on the :math:`\chi^2`
====================================

A relevant question to ask is what is the effect of possible cuts on the
input data set on the form of the :math:`\chi^2` in the presence of
correlations. In order to address this question, I consider a simple
data set made of two datapoints equipped with a single uncorrelated
uncertainty and a single fully-correlated uncertainty: [3]_

.. math::
   S:\quad\{m_i\pm \sigma_i^u\pm \sigma_i^c\}\,,\quad i  =1,2\,.
   :label: eq:fullset

\ The resulting covariance matrix reads:

.. math::
   \mathbf{V}=\begin{pmatrix}
   (\sigma_1^u)^2+ (\sigma_1^c)^2 & \sigma_1^c \sigma_2^c \\
   \sigma_1^c \sigma_2^c & (\sigma_2^u)^2+ (\sigma_2^c)^2
   \end{pmatrix}\,,

with inverse:

.. math::
   \mathbf{V}^{-1}=\frac{1}{(\sigma_1^u \sigma_2^u)^2+ (\sigma_1^c \sigma_2^u)^2+ (\sigma_1^u \sigma_2^c)^2}\begin{pmatrix}
   (\sigma_2^u)^2+ (\sigma_2^c)^2 & -\sigma_1^c \sigma_2^c \\
   -\sigma_1^c \sigma_2^c & (\sigma_1^u)^2+ (\sigma_1^c)^2
   \end{pmatrix}\,.

In addition, the set of two theoretical predictions :math:`\{t_i\}` that
allows one to define a column vector of residuals
:math:`\mathbf{r}^{\rm T}=(r_1,\;r_2)` with :math:`r_i=m_i-t_i`. The
:math:`\chi^2` then takes the following explicit expression:

.. math::
   \chi^2= \mathbf{r}^{\rm T}\cdot \mathbf{V}^{-1}\cdot \mathbf{r}=
     \frac{r_1^2\left[(\sigma_2^u)^2+
         (\sigma_2^c)^2\right]+r_2^2\left[(\sigma_1^u)^2+
         (\sigma_1^c)^2\right]-2 r_1 r_2\sigma_1^c
       \sigma_2^c}{(\sigma_1^u \sigma_2^u)^2+ (\sigma_1^c \sigma_2^u)^2+
       (\sigma_1^u \sigma_2^c)^2}\,.
   :label: eq:fullchi2

Now, suppose one wants to exclude the datapoint with :math:`i=2` from
the definition of the :math:`\chi^2`. The possibly most natural way to
proceed is to exclude the point, along with its uncertainties, directly
from the set in Eq. :eq:`eq:fullset`. This
straightforwardly leads to the following :math:`\chi^2`:

.. math::
   {\chi}_{\rm cut\;1}^2 = \frac{r_1^2}{(\sigma_1^u)^2+
     (\sigma_1^c)^2}\,.
   :label: eq:redchi21

However, this procedure might be questioned in that it artificially
modifies the original data set by effectively defining a *new* data set
where the point :math:`i=2` is not present. In the absence of
correlations this is a sound procedure because each single datapoint can
be regarder an independent subset of the original set. However, this is
not the case when correlations are present.

Therefore, it is necessary to devise a method that avoids any
modification of the data set but yet allows one to prevent specific
datapoints to contribute to the :math:`\chi^2`. In general, the purpose
of excluding datapoints from a fit is to avoid that the model used to
compute the predictions is forced outside its application region. One
can therefore *enforce* that the model exactly agrees with the
datapoints to be excluded. This effectively amounts to set
:math:`t_2=m_2`, or equivalently :math:`r_2=0`, reducing
Eq. :eq:`eq:fullchi2` to:

.. math::
   \chi_{\rm cut\;2}^2= 
     \frac{r_1^2\left[(\sigma_2^u)^2+ (\sigma_2^c)^2\right]}{(\sigma_1^u
       \sigma_2^u)^2+ (\sigma_1^c \sigma_2^u)^2+ (\sigma_1^u
       \sigma_2^c)^2}\,.
   :label: eq:fullchi2trunc

\ This has the clear advantage that the experimental information,
including correlations, remains totally intact. In order to show that
Eq. :eq:`eq:fullchi2trunc` has the right
properties, I first take the limit for :math:`\sigma_2^c\rightarrow 0`.
This gives:

.. math::
   \lim_{\sigma_2^c\rightarrow 0}\chi_{\rm cut}^2= 
     \frac{r_1^2}{(\sigma_1^u)^2+ (\sigma_1^c)^2}\,,

which reproduces the result of Eq. :eq:`eq:redchi21`.
This was to be expected because no reference to the point :math:`i=2`
can remain after the removal of its correlation to the point
:math:`i=1`. In addition, it is also very instructive to take the limit
:math:`\sigma_2^u\rightarrow 0`
Eq. :eq:`eq:fullchi2trunc`, which gives:

.. math::
   \lim_{\sigma_2^u\rightarrow 0}\chi_{\rm cut}^2= 
   \frac{r_1^2}{ (\sigma_1^u)^2}\,.

In this case, not only any reference to the point :math:`i=2` drops out,
but also the dependence of the :math:`\chi^2` on the correlation
uncertainty of the point :math:`i=1` disappears. Despite at first sight
this looks counterintuitive, it seems to be the correct behaviour. To
see this, observe that, if the point :math:`i=2` has no uncorrelated
uncertainty, its fluctuations are totally driven by the fluctuations of
the point :math:`i=1`. Therefore, the point :math:`i=2` is completely
dependent on the point :math:`i=1` and thus cannot give any contribution
to the :math:`\chi^2`. In addition, the correlation uncertainty of the
point :math:`i=1` must also be irrelevant because, being :math:`i=2`
totally correlated to :math:`i=1`, any correlation of :math:`i=1` to
:math:`i=2` translates into a correlation to itself and this cannot
contribute to the :math:`\chi^2` either. Importantly, if one takes the
limit :math:`\sigma_1^u\rightarrow 0`, the :math:`\chi^2` diverges no
matter the value of :math:`\sigma_1^c`. This is consistent with setting
:math:`\sigma_1^u=\sigma_2^u` since the beginning in
Eq. :eq:`eq:fullset` which indeed gives a singular
:math:`\chi^2` because the covariance matrix is singular.
Eq. :eq:`eq:redchi21` does not produce these features
hanging in favour of Eq. :eq:`eq:fullchi2trunc`.

As a further argument in favour of
Eq. :eq:`eq:fullchi2trunc` and against
Eq. :eq:`eq:redchi21`, I consider a larger data set
with, say, :math:`n` data points and split it into two subsets with
:math:`n_1` and :math:`n_2` data points, respectively, such that
:math:`n_1+n_2=n`. As I will explicitly show below,
Eq. :eq:`eq:fullchi2trunc` guarantees that the
:math:`\chi^2`\ ’s of the two subsets, :math:`\chi_1^2` and
:math:`\chi_2^2`, combine sensibly to give the total :math:`\chi^2`. I
will also show that this is not the case for
Eq. :eq:`eq:redchi21`. The fulfilment of this condition
is crucial to ensure that the experimental information carried by the
data set as a whole is preserved also when one is either not able or
does not want to make predictions for subset of data points. This
property is particularly relevant in fits that adopt the
cross-validation method as stopping criterion. In this context, the full
data set is indeed split into a training subset and a validation subset
and only the :math:`\chi^2` of the training set is minimised while the
validation :math:`\chi^2` is monitored. Since the cross-validation
method relies on the exploitation of the full experimental information
to stop the fit, it is necessary that training and validation
:math:`\chi^2`\ ’s combine meaningfully.

Now, I move to explicitly show that the generalisation of
Eq. :eq:`eq:fullchi2trunc` to an :math:`n`-point
data set where :math:`n_2` are cut away and :math:`n_1` remain fulfils
the requirement:

.. math::
   \chi^2 = \left|{\vec \chi}_1+{\vec \chi}_2\right|^2\,,
   :label: eq:compisition1

\ where the vectors :math:`{\vec \chi}_1` and :math:`{\vec \chi}_2` are
defined according to Eq. :eq:`eq:chivdef`. To define
them more precisely, without loss of generality, one can order the data
points in such a way that those labelled with :math:`i=1,\dots,n_1`
belong to the first subset and the points labelled with
:math:`i=n_1+1,\dots,n_1+n_2(=n)` belong to the second subset. The
entries of the :math:`{\vec \chi}_1` and :math:`{\vec \chi}_2` are then
given by:

.. math::
   \vec \chi_{1,k} = \sum_{j=1}^{n_1} L_{kj}^{-1}y_j\,,\quad\mbox{and}\quad
   \vec \chi_{2,k} = \sum_{j=n_1+1}^{n} L_{kj}^{-1}y_j\,,

where :math:`L_{ij}^{-1}` are the entries of the inverse of the Cholesky
decomposition :math:`\mathbf{L}` of the covariance :math:`\mathbf{V}`.
It should be clear that, according to the generalisation of the
prescription of Eq. :eq:`eq:fullchi2trunc`, one
has:

.. math:: \chi_1^2= \left|\vec \chi_{1}\right|^2\,,\quad\mbox{and}\quad \chi_2^2= \left|\vec \chi_{2}\right|^2\,.

The proof of Eq. :eq:`eq:compisition1` runs as
follows:

.. math::
   \begin{array}{rcl}
     \chi^2 &=&\displaystyle 
                \sum_{i=1}^n\sum_{j=1}^nr_iV_{ij}^{-1}r_j=\sum_{i=1}^{n_1}\sum_{j=1}^{n_1}r_iV_{ij}^{-1}r_j
                + \sum_{i=n_1+1}^{n}\sum_{j=n_1+1}^{n}r_iV_{ij}^{-1}r_j+2
                \sum_{i=1}^{n_1}\sum_{j=n_1+1}^{n}r_iV_{ij}^{-1}r_j\\
     \\
            &=& \displaystyle \sum_{k=1}^{n}\left(\sum_{i=1}^{n_1}L^{-1}_{ki}r_i\right)\left(\sum_{j=1}^{n_1}L^{-1}_{kj}r_j\right)
                +
                \sum_{k=1}^{n}\left(\sum_{i=n_1+1}^{n}L^{-1}_{ki}r_i\right)\left(\sum_{j=n_1+1}^{n}L^{-1}_{kj}r_j\right)\\
   \\
    &+&\displaystyle 2
                \sum_{k=1}^n\left(\sum_{i=1}^{n_1}L_{ki}^{-1}r_i\right)\left(\sum_{j=n_1+1}^{n}L_{kj}^{-1}r_j\right)\\
     \\
            &=&\displaystyle {\vec \chi}_1^{\rm T}\cdot{\vec \chi}_1+{\vec \chi}_2^{\rm T}\cdot{\vec \chi}_2 +2
                {\vec \chi}_1^{\rm T}\cdot{\vec \chi}_2 =\left|{\vec \chi}_1+{\vec \chi}_2\right|^2\,,
   \end{array}
   :label: eq:combinechi2s

where in the last step of the first line I have used the fact that the
covariance matrix (and thus its inverse) is symmetric upon exchange of
the indices while in the second and third lines I have used the Cholesky
decomposition of the covariance matrix. It straightforward to see that
Eq. :eq:`eq:redchi21` does not fulfils the same
relation. The fact that Eq. :eq:`eq:compisition1`
is consistent with how the :math:`\chi^2`\ ’s should combine can be
better seen by giving it a geometrical interpretation.

.. container::
   :name: a-geometrical-point-of-view

   .. rubric:: A geometrical point of view
      :name: a-geometrical-point-of-view

It is instructive to consider the question of cuts from a geometrical
point of view. In fact, the whole bunch of questions related to the
computation of the :math:`\chi^2` can be nicely visualised as a set of
geometrical operations in a multidimensional space :math:`\mathcal{M}`
(manifold) equipped with the appropriate metric tensor. Based on the
definition in Eq. :eq:`eq:chi2cov`, the :math:`\chi^2`
can be interpreted as the squared absolute value of the
:math:`n`-dimensional vector :math:`\mathbf{y}` given by the difference
between the vector of experimental central values :math:`\mathbf{m}` and
the vector of theoretical predictions :math:`\mathbf{t}`. Importantly,
the scalar product required to compute this absolute value has to be
computed using the inverse of the covariance matrix
:math:`\mathbf{V}^{-1}` as a covariant positive-definite metric
tensor. [4]_ Borrowing the Einstein notation to represent scalar
products, *i.e.* defining contravariant vectors with upper index as:

.. math:: y^i = V^{-1}_{ij}y_j\,,

where the sum over :math:`j` between 1 and :math:`n` is understood, the
:math:`\chi^2` is the *norm* squared of the vector :math:`\mathbf{y}` in
the space :math:`\mathcal{M}`.

.. math::
   \chi^2 = \mathbf{y}^{\rm T}\cdot \mathbf{V}^{-1}\cdot \mathbf{y} =
     y^iy_i\equiv \left\|\mathbf{y}\right\|^2\,.

It is interesting to notice that the minimisation of the :math:`\chi^2`
can then be regarded as the attempt to align the vector
:math:`\mathbf{t}` to the vector :math:`\mathbf{m}` in the space
:math:`\mathcal{M}` in such a way that their difference is as small as
possible.

The enforcement of cuts can be accommodated in this framework by
requiring that the :math:`\chi^2` is the squared absolute value (or
length) of the projection of :math:`\mathbf{y}` in the subspace of
:math:`\mathcal{M}`, :math:`\mathcal{M}_{1}`, defined by the components
of the vector :math:`\mathbf{m}` corresponding to the points that pass
the cut. Introducing the operator :math:`\mathbf{P}` that projects any
vector in :math:`\mathcal{M}` onto :math:`\mathcal{M}_{1}`, the
corresponding :math:`\chi^2` will read:

.. math:: \chi^2_{1} = \left(\mathbf{P}\cdot \mathbf{y}\right)^{\rm T}\cdot \mathbf{V}^{-1}\cdot \left(\mathbf{P}\cdot \mathbf{y}\right)=\left(\mathbf{P}\cdot \mathbf{y}\right)^i \left(\mathbf{P}\cdot \mathbf{y}\right)_i=\left\|\mathbf{P}\cdot\mathbf{y}\right\|^2\,.

The matrix corresponding to the projector :math:`\mathbf{P}` is simply
the diagonal matrix with zero’s and one’s on the diagonal depending on
whether the point is excluded or included, respectively. [5]_ The net
effect of the projector :math:`\mathbf{P}` is to set to zero the
components of the vector of residuals :math:`\mathbf{y}` corresponding
to the points excluded by the cut. This is in agreement with
Eq. :eq:`eq:fullchi2trunc` in the case of a
two-dimensional space.

It is interesting to define the projector
:math:`\mathbf{Q}=\mathbb{I}-\mathbf{P}` orthogonal to
:math:`\mathbf{P}`\  [6]_ such that the corresponding :math:`\chi^2`
reads:

.. math:: \chi^2_{2} = \left(\mathbf{Q}\cdot \mathbf{y}\right)^{\rm T}\cdot \mathbf{V}^{-1}\cdot \left(\mathbf{Q}\cdot \mathbf{y}\right) =\left(\mathbf{Q}\cdot \mathbf{y}\right)^i \left(\mathbf{Q}\cdot \mathbf{y}\right)_i =\left\|\mathbf{Q}\cdot\mathbf{y}\right\|^2\,.

With these definitions at hands, one can write the total :math:`\chi^2`
as:

.. math:: \chi^2= \left\|\mathbf{y}\right\|^2 = \left\|\mathbf{P}\cdot\mathbf{y}+\mathbf{Q}\cdot\mathbf{y}\right\|^2 = \chi^2_{1} + \chi^2_{2} + 2 \left(\mathbf{P}\cdot \mathbf{y}\right)^i \left(\mathbf{Q}\cdot \mathbf{y}\right)_i\,.

that agrees with Eq. :eq:`eq:combinechi2s`. In
conclusion, the introduction of a cut on the full data set has a nice
geometrical interpretation: the resulting :math:`\chi^2` corresponds to
the norm of the projection of the vector or residual :math:`\mathbf{y}`,
defined on the manifold :math:`\mathcal{M}` with covariant metric
:math:`\mathbf{V}^{-1}`, on to the subspace :math:`\mathcal{M}_1`
defined by the projector :math:`\mathbf{P}`. This is consistent with
Eq. :eq:`eq:compisition1` and provides it with a
nice geometrical interpretation. [7]_

Monte Carlo replica generation
==============================

In this section, I consider the Monte Carlo generation of artificial
replicas. The aim is to obtain a formula for the generation of replicas
that gives rise to a figure of merit distributed correctly, *i.e.* that
follows a :math:`\chi^2` distribution with as many degrees of freedom as
number of data points. Following Ref. (contributors 2020), one requires
that the :math:`\chi^2`, that in matricial form can be expressed as:

.. math::
   \chi^2 = \mathbf{y}^{T}\cdot \mathbf{V}^{-1} \cdot \mathbf{y}\,,
   :label: eq:chi2matrix

\ where :math:`\mathbf{V}` is the :math:`n\times n` covariance matrix,
being :math:`n` the number of datapoints, and :math:`\mathbf{y}` is a
column vector with :math:`n` entries that can be expressed as:

.. math:: y_j = f_j-m_j\,,

with :math:`m_j` a given value (to be thought as the central value) and
:math:`f_j` its fluctuated value, is distributed like:

.. math::
   \chi^2\sim \sum_{i=1}^{n}z_i^2=|\mathbf{z}|^2\,,
   :label: eq:wikipedia

where :math:`z_1,\dots , z_n` are independent, standard normal random
variables, *i.e.* such that:

.. math::
   \langle z_i\rangle=0\quad\mbox{and}\quad \langle z_iz_j\rangle =
   \delta_{ij}\,,
   :label: eq:averages

where the symbol :math:`\langle\,\dots\rangle` indicates the
average. [8]_ This immediately tells that:

.. math:: \langle\chi^2\rangle=n\,.

As already discussed above, Eq. :eq:`eq:chi2matrix`
can alternatively be written as:

.. math::
   \chi^2 = \left|\mathbf{L}^{-1}\cdot \mathbf{y}\right|^2\,,
   :label: explicitchi2chol

where :math:`\mathbf{L}` is the Cholesky decomposition of the matrix
:math:`\mathbf{V}`, such that:

.. math::
   \mathbf{V} = \mathbf{L}\cdot\mathbf{L}^{T}\,.
   :label: eq:choleskydecagain

Comparing Eqs. :eq:`explicitchi2chol`
and :eq:`eq:wikipedia`, one can infer that the
following relation should hold:

.. math:: \mathbf{L}^{-1}\cdot \mathbf{y}\sim\mathbf{z}\,,

which immediately implies:

.. math:: \mathbf{y}\sim\mathbf{L}\cdot\mathbf{z}\,,

or equivalently:

.. math::
   f_j\sim m_j+\sum_{i=1}^nL_{ji}z_i\,.
   :label: eq:deviations

In conclusion, in order for the :math:`\chi^2` to have the correct
distribution, the fluctuations :math:`f_j` around the experimental
central values :math:`m_j` of a given Monte Carlo replica have to be
computed according to Eq. :eq:`eq:deviations`.

I now show that Eq. :eq:`eq:deviations` reproduces
the statistical properties of the original ensamble of data. To do so, I
compute average and correlations of
Eq. :eq:`eq:deviations` using
Eq. :eq:`eq:averages`. The average gives:

.. math:: \langle f_j\rangle = m_j+\sum_{i=1}^nL_{ji}\langle z_i \rangle = m_j\,,

as expected. The correlation instead gives:

.. math::
   \begin{array}{rcl}
   \langle f_j f_k\rangle &=& \displaystyle \left\langle
     \left(m_j+\sum_{i=1}^nL_{ji}z_i\right)
     \left(m_k+\sum_{i'=1}^nL_{ki'}z_{i'}\right)\right\rangle\\
   \\
   &=& \displaystyle m_jm_k+\sum_{i,i'=1}^nL_{ji'}L_{ki'}\langle z_{i'}
       z_{i'}\rangle=\langle f_j \rangle\langle f_k\rangle+V_{jk}\,,
   \end{array}

again as expected. These are fundamental requirements of a sound Monte
Carlo ensamble (see also Eq. (10) of Ref. (Ball et al. 2010)).

Eq. :eq:`eq:deviations` should be compared to the
more “standard” Monte-Carlo replica generation formulas such as those
given in Refs. (Ball et al. 2009; Ball and others 2015). I now
specifically consider Eq. (20) of Ref. (Ball and others 2015) that I
report here in terms of relative uncertainties :math:`\delta`\ ’s
matching our notation:

.. math::
   \widetilde{f}_j\sim m_j\prod_{m=1}^{k_m}\left(1+\delta_{{\rm mult},j}^{(m)} x_m\right)\left(1+\delta_{{\rm unc},j} z_j+\sum_{a=1}^{k_a} \delta_{{\rm add},j}^{(a)} y_a\right)\,,
   :label: eq:nnpdfmc

\ where :math:`z_j`, :math:`y_a`, and :math:`x_m` are independent random
variables. I now expand the product in the front up to second order in
:math:`x`

.. math::
   \prod_{m=1}^{k_m}\left(1+\delta_{{\rm mult},j}^{(m)}
       x_m\right)=1+\sum_{m=1}^{k_m}\delta_{{\rm mult},j}^{(m)}
     x_m+\sum_{m=1}^{k_m}\sum_{m'=1\atop m'\neq m}^{k_m}\delta_{{\rm
         mult},j}^{(m)}\delta_{{\rm mult},j}^{(m')}x_m x_{m'}+\mbox{trilinear}\,,

where the “trilinear” corrections correspond to combinations of the kind
:math:`x_ix_jx_k`, with :math:`i\neq j\neq k\neq i`, such that their
average vanishes:
:math:`\langle x_ix_jx_k\rangle = \langle x_i\rangle \langle x_j
\rangle\langle x_k\rangle = 0`. This allows one to expand
Eq. :eq:`eq:nnpdfmc` as:

.. math::
   \begin{array}{rcl}
   \widetilde{f}_j&\sim&\displaystyle m_j\Bigg(1+\delta_{{\rm unc},j} z_j+\sum_{a=1}^{k_a} \delta_{{\rm add},j}^{(a)} y_a +\sum_{m=1}^{k_m}\delta_{{\rm mult},j}^{(m)}
     x_m\\
   \\
   &+&\displaystyle\delta_{{\rm unc},j} \sum_{m=1}^{k_m} \delta_{{\rm mult},j}^{(m)}
       z_jx_m+\sum_{a=1}^{k_a} \sum_{m=1}^{k_m}\delta_{{\rm add},j}^{(a)} \delta_{{\rm mult},j}^{(m)} y_a 
     x_m\\
   \\
   &+&\displaystyle \sum_{m=1}^{k_m}\sum_{m'=1\atop m'\neq m}^{k_m}\delta_{{\rm
         mult},j}^{(m)}\delta_{{\rm mult},j}^{(m')}x_m x_{m'}\Bigg)+\mbox{trilinear}\,.
   \end{array}
   :label: eq:nnpdfmcexp

Since Eq. :eq:`eq:nnpdfmcexp` involves only terms
that have at most one power of the same random variable, the average of
:math:`\widetilde{f}_j` automatically yields:

.. math::
   \langle \widetilde {f}_j\rangle = m_j\,,
   :label: eq:centralvaluenn

in agreement with the expectation. I now compute the correlations making
use of the combination rules of random variables. This gives:

.. math::
   \begin{array}{rcl}
     \langle \widetilde{f}_j \widetilde{f}_k\rangle&\sim&\displaystyle m_jm_k\Bigg(1+\delta_{{\rm unc},j}^2\delta_{jk}+\sum_{a=1}^{k_a} \delta_{{\rm add},j}^{(a)}\delta_{{\rm add},k}^{(a)} +\sum_{m=1}^{k_m}\delta_{{\rm mult},j}^{(m)}\delta_{{\rm mult},k}^{(m)}\\
   \\
   &+&\displaystyle\delta_{{\rm unc},j}^2 \sum_{m=1}^{k_m} \delta_{{\rm mult},j}^{(m)}\delta_{{\rm  mult},k}^{(m)}
     +\sum_{a=1}^{k_a} \sum_{m=1}^{k_m}\delta_{{\rm add},j}^{(a)}\delta_{{\rm add},k}^{(a)} \delta_{{\rm mult},j}^{(m)} \delta_{{\rm mult},k}^{(m)} \\
   \\
   &+&\displaystyle \sum_{m=1}^{k_m}\sum_{m'=1\atop m'\neq m}^{k_m}\delta_{{\rm
         mult},j}^{(m)}\delta_{{\rm
         mult},k}^{(m)}\delta_{{\rm mult},j}^{(m')}\delta_{{\rm
       mult},k}^{(m')}+\mathcal{O}(\delta^6)\Bigg)\\
   \\
   &=&\displaystyle \langle \widetilde{f}_j\rangle \langle
       \widetilde{f}_k\rangle +V_{jk}+\mathcal{O}(\sigma_{\rm mult}^2 s^2, \sigma_{\rm mult}^2 \sigma_{\rm add}^2)\,.
   \end{array}
   :label: eq:corrnn

In conclusion, Eq. :eq:`eq:nnpdfmc` reproduces exactly
the central values of the original ensemble (see
Eq. :eq:`eq:centralvaluenn`) while it reproduces
correlations only up to corrections proportional to the fourth power of
the typical uncertainty. Interestingly, these corrections are all
proportional to :math:`\sigma_{\rm mult}^2`. This means that, in the
absence of multiplicative uncertainties, somewhat expectedly
Eq. :eq:`eq:nnpdfmc` *exactly* reproduces the
correlations. On the other hand, if multiplicative uncertainties are
large these correction become potentially sizeable preventing
Eq. :eq:`eq:nnpdfmc` from correctly reproducing the
correlations of the original ensemble.

The case of Eq. (13) of Ref. (Ball et al. 2009) is more involved because
of the presence of the square-root sign for the relative multiplicative
uncertainties. This case can be handled again be expanding to second
order the factor depending on multiplicative uncertainties with the only
difference that, due to the expansion:

.. math:: \sqrt{1+x}=1+\frac12x-\frac18x^2+\mathcal{O}(x^3)\,,

\ the square root produces a factor :math:`1/2` in front of the linear
terms and a factor :math:`-1/8` in front of the quadratic terms. This is
expected to reduce the impact of the corrections terms in
Eq. :eq:`eq:corrnn`.

It is instructive to work out whether
Eq. :eq:`eq:deviations` can be directly related to
the Monte-Carlo replicas generation formulas given in Refs. (Ball et al.
2009; Ball and others 2015). The aim is to find a closed-form expression
for the Cholesky-decomposition matrix :math:`\mathbf{L}` of the
covariance matrix :math:`\mathbf{V}` in terms of the single
uncertainties. Through Eq. :eq:`eq:deviations`, this
will finally allow one to derive a formula of the same kind of those
given in Refs. (Ball et al. 2009; Ball and others 2015). The basic
observation is that the covariance matrix can be written as:

.. math::
   \mathbf{V} = \mathbf{D}+ \sum_{l=1}^k 
   \mathbf{V}^{(l)} = \sum_{l=0}^k \mathbf{L}^{(l)}\cdot \mathbf{L}^{(l)T}\,,

with:

.. math::
   L_{ij}^{(l)}=\left\{
   \begin{array}{ll}
   s_i\delta_{ij}&\quad l = 0\,,\\
   \sigma_i^{(l)}\delta_{j1}&\quad 1\leq l\leq k\,.
   \end{array}
   \right.

On the other hand, one also knows that:

.. math:: \mathbf{V} = \mathbf{L}\cdot \mathbf{L}^T\,,

therefore:

.. math:: \mathbf{L}\cdot \mathbf{L}^T=\sum_{l=0}^k \mathbf{L}^{(l)}\cdot \mathbf{L}^{(l)T}\,.

The question is to express :math:`\mathbf{L}` in terms of all the
:math:`\mathbf{L}^{(l)}`. Unfortunately though this equation has no
general solution. However, a step forward can be done be exploiting the
fact that the matrices :math:`\mathbf{V}^{(l)}` are positive
semidefinite (*i.e.* they have zero determinant) and thus their Cholesky
decomposition is not unique. In particular, provided that
:math:`k\leq n`, one can choose:

.. math::
   \widetilde{L}_{ij}^{(l)}=
   \sigma_i^{(l)}\delta_{jl}\,,\quad 1\leq l\leq k\leq n\,,

such that:

.. math::
   \left(\sum_{l=1}^k\widetilde{\mathbf{L}}^{(l)}\right)\cdot\left(
       \sum_{l'=1}^k\widetilde{\mathbf{L}}^{(l')}\right)=\sum_{l=1}^k 
     \mathbf{V}^{(l)}\,.
   :label: eq:choldecomp2

However, choosing:

.. math::
   \widetilde{\mathbf{L}} = \mathbf{L}^{(0)} +
   \sum_{l=1}^k\widetilde{\mathbf{L}}^{(l)}\,,
   :label: eq:CholVapprox

and multiplying it for its transpose would give:

.. math:: \widetilde{\mathbf{L}}\cdot \widetilde{\mathbf{L}}^T = \mathbf{V} + \sum_{l=1}^k\left(\mathbf{L}^{(0)}\cdot \widetilde{\mathbf{L}}^{(l)T}+\widetilde{\mathbf{L}}^{(l)}\cdot \mathbf{L}^{(0)T}\right)\,,

that gives back the covariance matrix plus a spurious term. Writing the
expression above in terms of single entries, one finds:

.. math:: \sum_{m=1}^n\widetilde{L}_{im}\widetilde{L}_{jm} = V_{ij} + s_i \sigma_{j}^{(i)}+s_j\sigma_{i}^{(j)}\,.

The peculiarity here is that the upper index of :math:`\sigma_{j}^{(i)}`
and :math:`\sigma_{i}^{(j)}`, since :math:`k\leq n` (and in fact usually
:math:`k\ll n`, *i.e.* the number of correlated uncertainties is
typically much smaller than the number of data points), exceeds the
allowed range. To fix this problem, I just assume
:math:`\sigma_{j}^{(i)}=0` for :math:`i>k`. This finally gives:

.. math::
   \sum_{m=1}^n\widetilde{L}_{im}\widetilde{L}_{jm} = V_{ij}\,,\quad
   \mbox{for } i,j>k\,,

that means that, by adopting the Cholesky decomposition in
Eq. :eq:`eq:choldecomp2`, part of the covariance
matrix (usually the largest part) is exactly recovered. Therefore,
:math:`\widetilde{\mathbf{L}}`
Eq. :eq:`eq:CholVapprox` seems to be a generally
good approximator to the real Cholesky decomposition of the covariance
matrix :math:`\mathbf{V}`. Plugging this equation into
Eq. :eq:`eq:deviations` immediately gives:

.. math::
   f_j\sim m_j+z_js_j + \sum_{l=1}^{k}z_l \sigma_{j,\rm
       corr}^{(l)}\,.
   :label: eq:deviations2

It is interesting to focus on similarities and differences with the
formulas present in the literature. I will specifically consider
Eq. (13) of Ref. (Ball et al. 2009) and Eq. (20) of Ref. (Ball and
others 2015).

Let us start by considering the case in which there are only
uncorrelated and correlated additive uncertainties. In this case,
Refs. (Ball et al. 2009; Ball and others 2015) are in mutual agreement
and only in partial agreement with our
Eq. :eq:`eq:deviations2`, the difference being that
the random numbers used to generate the fluctuation associated to the
:math:`k` correlated uncertainties are equal to the first :math:`k` of
the :math:`n` used for the uncorrelated uncertainties (remember the
requirement :math:`k\leq n`). Of course, if :math:`k\ll n`, as it often
happens, the agreement between
Eq. :eq:`eq:deviations2` and Refs. (Ball et al.
2009; Ball and others 2015) is asymptotically good.

When considering also correlated multiplicative uncertainties there are
differences across formulas of Refs. (Ball et al. 2009; Ball and others
2015) and Eq. :eq:`eq:deviations2`. While
Eq. :eq:`eq:deviations2` treats multiplicative
uncertainties on the same footing as additive uncertainties
(consistently with the construction of the covariance matrix in
Eq. :eq:`eq:covmat`), Refs. (Ball et al. 2009; Ball and
others 2015) introduce a multiplicative factor. In addition, Ref. (Ball
et al. 2009) distinguishes between relative and absolute multiplicative
uncertainties.

|image| |image1|

The D’Agostini “bias”
=====================

The main topic of this section is the so-called D’Agostini “bias” (DAB)
originally presented in Ref. (D’Agostini 1994). Quoting from the
abstract of that paper: “Best fits to data which are affected by
systematic uncertainties on the normalisation factor have the tendency
to produce curves lower than expected [...]. This bias can become
unacceptable if the normalisation error is large, or a large number of
data points are used.” In the following, I will summarise the arguments
behind this conclusion. I will then move to giving my own view on this
question arguing that the DAB is in fact not a bias but rather a feature
of correlation uncertainties in general and not of multiplicative
uncertainties in particular. This implies that the distinction between
additive and normalisation uncertainties is effectively immaterial.
Specifically, what really triggers the DAB is *not* the fact the
uncertainty is on the multiplication factor but rather that it typically
gives raise to correlated uncertainties that are different in *absolute*
size. To show this, below I will show that the DAB takes place also in
the presence of additive uncertainties that are different in absolute
size. Conversely, I will show that the DAB does not takes place when
considering multiplicative uncertainties and the central values are
identical. The bottom line is that any difference in absolute size
between correlated uncertainties causes a deviation from the expectation
according to which the best fit value should be given by the weighted
average of the experimental central values weighted by the inverse of
the squared uncorrelated uncertainties. This deviation from the
expectation is the origin of the DAB.

In order to present the argument, Ref. (D’Agostini 1994) considers a
data set made of two data points with the respective uncorrelated
uncertainty, :math:`x_i\pm\sigma_i` with :math:`i=1,2`. In addition, it
is assumed that :math:`x_i` measure the same quantity and thus they can
be fitted with a single constant :math:`\mu`. Two cases are treated. In
the first case the two data points are also affected by a single offset
(additive) correlated uncertainty :math:`\sigma_c`. In the second case,
instead, the data points are affected by a normalisation correlated
uncertainty whose size is proportional to the central values through the
factor :math:`\delta_f`. Given the simplicity of the problem, it is
possible in both cases to analytically compute the best value of
:math:`\mu` and its uncertainty :math:`\delta\mu` by minimising the
:math:`\chi^2` and computing its second derivative at the minimum. As
usual the :math:`\chi^2` is defined as:

.. math::
   \chi^2 = \sum_{i,j=1}^2y_iV^{-1}_{ij}y_j = \mathbf{y}^T\cdot
   \mathbf{V}^{-1}\cdot \mathbf{y}\,,

\ where :math:`y_i=\mu-x_i` and :math:`V^{-1}_{ij}` are the entries of
the inverse of the covariance matrix.

Let us start with the offset case in which the covariance matrix takes
the form:

.. math::
   \mathbf{V}=
   \begin{pmatrix}
   \sigma_1^2+\sigma_c^2& \sigma_c^2\\
   \sigma_c^2& \sigma_2^2+\sigma_c^2
   \end{pmatrix}\,,

\ whose inverse is:

.. math::
   \mathbf{V}^{-1}=\frac{1}{D}
   \begin{pmatrix}
   \sigma_2^2+\sigma_c^2& -\sigma_c^2\\
   -\sigma_c^2& \sigma_1^2+\sigma_c^2
   \end{pmatrix}\,,

with :math:`D=\sigma_1^2 \sigma_2^2+(\sigma_1^2+\sigma_2^2)
\sigma_c^2`. With this I can straightforwardly compute the
:math:`\chi^2`:

.. math:: \chi^2=\frac{\sigma_2^2y_1^2+\sigma_1^2y_2^2+(y_1-y_2)^2\sigma_c^2}{D}\,,

that can be minimised w.r.t. to :math:`\mu` by requiring that:

.. math::
   \left.\frac{d\chi^2}{d\mu}\right|_{\mu=\overline{\mu}} =
     \left.2\frac{\sigma_2^2y_1+\sigma_1^2y_2}{D}\right|_{\mu=\overline{\mu}}=0\,,
   :label: eq:minchi2

where :math:`\overline{\mu}` is the best fit value of :math:`\mu`.
Notice the cancellation in the numerator of the correlated uncertainty
:math:`\sigma_c`. This “accidental” cancellation in the offset case,
caused by the fact that :math:`\sigma_c` is the same for both
:math:`x_1` and :math:`x_2`, can be regarded as the (misleading)
evidence that the DAB only takes place in the normalisation case. The
solution to Eq. :eq:`eq:minchi2` is:

.. math::
   \overline{\mu} =
   \frac{\sigma_2^2x_1+\sigma_1^2x_2}{\sigma_1^2+\sigma_2^2}=\frac{\frac{x_1}{\sigma_1^2}+\frac{x_2}{\sigma_2^2}}{\frac{1}{\sigma_1^2}+\frac{1}{\sigma_2^2}}\,,
   :label: eq:offcv

that fulfils the intuitive expectation that the best value of
:math:`\mu` should be a weighted average of the central values weighted
by the inverse of the uncorrelated uncertainties squared. One can also
compute the uncertainty of :math:`\overline{\mu}`, :math:`\delta\mu`, by
simply computing the inverse of the second derivative of the
:math:`\chi^2` at :math:`\mu=\overline{\mu}`: [9]_

.. math::
   \delta\mu^2 =
     \left[\frac12\frac{d^2\chi^2}{d\mu^2}\right]_{\mu=\overline{\mu}}^{-1}
     = \frac {\sigma_1^2 \sigma_2^2+(\sigma_1^2+\sigma_2^2)
       \sigma_c^2}{\sigma_2^2+\sigma_1^2}=\frac
     {1+\frac{\sigma_c^2}{\sigma_1^2}+\frac{\sigma_c^2}{\sigma_2^2}
     }{\frac{1}{\sigma_1^2}+\frac{1}{\sigma_2^2}}\,.
   :label: eq:offcvunc

I now move to considering the normalisation case in which the covariance
matrix reads:

.. math::
   \mathbf{V}=
   \begin{pmatrix}
   \sigma_1^2+\delta_f^2x_1^2& \delta_f^2x_1x_2\\
   \delta_f^2x_1x_2& \sigma_2^2+\delta_f^2x_2^2
   \end{pmatrix}\,,

\ whose inverse is:

.. math::
   \mathbf{V}^{-1}=\frac{1}{D}
   \begin{pmatrix}
   \sigma_2^2+\delta_f^2x_2^2& -\delta_f^2x_1x_2\\
   -\delta_f^2x_1x_2& \sigma_1^2+\delta_f^2x_1^2
   \end{pmatrix}\,,

with
:math:`D=\sigma_1^2\sigma_2^2+\delta_f^2(\sigma_1^2x_2^2+\sigma_2^2x_1^2)`.
The resulting :math:`\chi^2` is:

.. math:: \chi^2=\frac{1}{D}\left[\sigma_2^2y_1^2+ \sigma_1^2y_2^2+\delta_f^2(x_2y_1 -x_1y_2)^2\right]\,.

The minimum is computed as:

.. math:: \left.\frac{d\chi^2}{d\mu}\right|_{\mu=\overline{\mu}} = \left.2\frac{\sigma_2^2y_1+ \sigma_1^2y_2+\delta_f^2(x_2 -x_1)(x_2y_1 -x_1y_2)}{D}\right|_{\mu=\overline{\mu}}=0\,.

This time, differently from the offset case, the correlated
uncertainties in the denominator do not cancel and the result of the
equation above is:

.. math::
   \overline{\mu} =
   \frac{\sigma_2^2x_1+\sigma_1^2x_2}{\sigma_1^2+\sigma_2^2+(\delta_f
     x_1-\delta_f
     x_2)^2}=\frac{\frac{x_1}{\sigma_1^2}+\frac{x_2}{\sigma_2^2}}{\frac{1}{\sigma_1^2}+\frac{1}{\sigma_2^2}+\left(\frac{\delta_f
         x_1-\delta_f x_2}{\sigma_1\sigma_2}\right)^2}\,,
   :label: eq:normcv

while:

.. math::
   \delta\mu^2 =
     \left[\frac12\frac{d^2\chi^2}{d\mu^2}\right]_{\mu=\overline{\mu}}^{-1}
     =
     \frac{\sigma_1^2\sigma_2^2+\delta_f^2(\sigma_1^2x_2^2+\sigma_2^2x_1^2)}{\sigma_2^2+
       \sigma_1^2+\delta_f^2(x_1 -x_2)^2}=\frac{1+\frac{\delta_f^2
         x_1^2}{\sigma_1^2}+\frac{\delta_f^2
         x_2^2}{\sigma_2^2}}{\frac{1}{\sigma_1^2}+
       \frac{1}{\sigma_2^2}+\left(\frac{\delta_f x_1 -\delta_f
           x_2}{\sigma_1\sigma_2}\right)^2}\,.
   :label: eq:normcvunc

According to Ref. (D’Agostini 1994), the term proportional to
:math:`\delta_f` in the denominator of the r.h.s. of
Eq. :eq:`eq:normcv` is a source of bias in that it
deviates from the intuitive expectation of
Eq. :eq:`eq:offcv`. Being this additional term positive,
it tends to diminish the best fit value.

I will now show that the additional term in the denominator of
Eq. :eq:`eq:normcv` is not a specific feature of
normalisation uncertainties. To show this, I consider the same data set
of two data points but where the data point :math:`x_1` is affected by
the additive uncertainty :math:`\sigma_{c,1}` and the point :math:`x_2`
by the additive uncertainty :math:`\sigma_{c,2}`. According to the same
procedure followed above the final result for best-fit value and
uncertainty of :math:`\mu` is:

.. math::
   \overline{\mu} =
   \frac{\frac{x_1}{\sigma_1^2}+\frac{x_2}{\sigma_2^2}}{\frac{1}{\sigma_1^2}+\frac{1}{\sigma_2^2}+\left(\frac{\sigma_{c,1}-\sigma_{c,2}}{\sigma_1\sigma_2}\right)^2}\,,
   :label: eq:add2cv

\ and:

.. math:: \delta\mu^2=\frac{1+\frac{\sigma_{c,1}^2}{\sigma_1^2}+\frac{\sigma_{c,2}^2}{\sigma_2^2}}{\frac{1}{\sigma_1^2}+ \frac{1}{\sigma_2^2}+\left(\frac{\sigma_{c,1} -\sigma_{c,2}}{\sigma_1\sigma_2}\right)^2}\,.

These equations are (somewhat expectedly) equal to
Eqs. :eq:`eq:normcv`
and :eq:`eq:normcvunc` with the only difference that
:math:`\delta_fx_i` is replaced by :math:`\sigma_{c,i}`. Therefore, the
DAB shows up also when additive uncertainties are considered as long as
they are different in absolute value. Therefore, the DAB is not a
prerogative of normalisation uncertainties only. As expected, the
equations above reduce to Eqs. :eq:`eq:offcv`
and :eq:`eq:offcvunc` when
:math:`\sigma_{c,1}=\sigma_{c,2}=\sigma_{c}`.

Another important aspect of Eqs. :eq:`eq:normcv`
and :eq:`eq:normcvunc` is that the term that causes
the “bias” is proportional to the difference :math:`x_1-x_2`. Therefore,
if the central values happen to be the same that term vanishes and thus
no bias is present. This observation, along with the fact that the DAB
is potentially present also when additive uncertainties are considered,
leads to a simple conclusion: the DAB is present in all cases whenever
there are correlated uncertainties different in absolute value. This is
trivially the case when one introduces different correlated additive
uncertainties, :math:`\sigma_{c,1}` and :math:`\sigma_{c,2}` with
:math:`\sigma_{c,1}\neq \sigma_{c,2}`. But this also the case in the
presence of a normalisation uncertainties when central values are
different (:math:`x_1\neq x_2`). On the contrary, if
:math:`\sigma_{c,1}= \sigma_{c,2}` in the additive case or
:math:`x_1=x_2` in the multiplicative case, no “bias” is present. The
conclusion is that the very distinction between additive and
multiplicative uncertainties is immaterial.

Now that I have established that the DAB is not a prerogative of
normalisation uncertainties only, I will argue that this is in fact not
a bias but a manifestation of large correlated uncertainties in general
that indicates an intrinsic inconsistency of the data set being
considered. I start by observing that, in the presence of correlations,
a visual comparison between data and best fit value is meaningful only
when including the systematic shifts introduces in
Sect. `3 <#eq:sysshifts>`__. The task is particularly simple due to the
small number of points and the fact that there one single correlated
uncertainty. In particular, there is a single nuisance parameter that in
the normalisation case is given by:

.. math:: \lambda=\frac{\frac{\delta_f x_1^2}{\sigma_1^2}+\frac{\delta_f x_2^2}{\sigma_2^2}}{1+\frac{\delta_f^2 x_1^2}{\sigma_1^2}+\frac{\delta_f^2 x_2^2}{\sigma_2^2}}-\overline{\mu}\frac{\frac{\delta_f x_1}{\sigma_1^2}+\frac{\delta_f x_2}{\sigma_2^2}}{1+\frac{\delta_f^2 x_1^2}{\sigma_1^2}+\frac{\delta_f^2 x_2^2}{\sigma_2^2}}\,.

It is the quantity:

.. math::
   \hat{\mu}_i =\overline{\mu}+\lambda \delta_fx_i\,,
   :label: eq:shiftedvals

that can be directly compared to the data-point central values.
Unfortunately, the analytic form of this quantity is not very
informative but to see the effect of the systematic shifts I consider
the very same example given in the introduction of Ref. (D’Agostini
1994), *i.e.* a data set of two measurements: :math:`8.0\pm 2\%` and
:math:`8.5\pm 2\%` having a common 10% normalisation uncertainty.

.. container::
   :name: fig:DAgostiniExample

   .. figure:: ../../latex/src/../../latex/src/plots/DAgostiniExample.pdf
   .. :alt: Numerical example taken from Ref. (D’Agostini 1994) of a
      data set of two points with central value 8.0 and 8.5 (black
      points) each affected by a 2% uncorrelated uncertainty (solid
      error bars) and a common 10% normalisation uncertainty (the dashed
      error bars indicate the sum in quadrature of correlated and
      uncorrelated uncertainties). The horizontal red line corresponds
      to the best fit value given in Eq. :eq:`eq:offcv`
      while the horizontal blue lines correspond to the shifted values
      according to
      Eq. :eq:`eq:shiftedvals`.[fig:DAgostiniExample]
      :name: fig:DAgostiniExample
      :width: 70.0%

      Numerical example taken from Ref. (D’Agostini 1994) of a data set
      of two points with central value 8.0 and 8.5 (black points) each
      affected by a 2% uncorrelated uncertainty (solid error bars) and a
      common 10% normalisation uncertainty (the dashed error bars
      indicate the sum in quadrature of correlated and uncorrelated
      uncertainties). The horizontal red line corresponds to the best
      fit value given in Eq. :eq:`eq:offcv` while the
      horizontal blue lines correspond to the shifted values according
      to
      Eq. :eq:`eq:shiftedvals`.[fig:DAgostiniExample]

Fig. `3 <#fig:DAgostiniExample>`__ shows the result of minimising the
:math:`\chi^2` to determine the best fit value of the parameter
:math:`\mu`. The solid error bars display the uncorrelated uncertainties
while the dashed error bars display the sum in quadrature of correlated
and uncorrelated uncertainties. The best fit value
:math:`\overline{\mu}`, Eq. :eq:`eq:offcv`, is represented
by the horizontal red line. As argued in
Fig. `3 <#fig:DAgostiniExample>`__, this line lies below both the
central experimental values (:math:`\overline{\mu}\simeq 7.87`)
contradicting the common expectation that best fit value should be
somewhere between the experimental central values. However, as explained
above, a visual comparison in the presence of (large) correlations is
meaningful only when including the systematic shifts. The inclusion of
the shifts yields the two horizontal blue lines in
Fig. `3 <#fig:DAgostiniExample>`__. The shifted values are
:math:`\hat{\mu}_1 \simeq 8.22` and :math:`\hat{\mu}_2 \simeq 8.25` that
indeed fulfil the common expectation. It is worth noticing that the best
fit :math:`\chi^2` per degree of freedom is around 4.4, signifying that
the measurements, that are assumed to refer to the same quantity, are in
strong tension. This is a first suggestion that that there is no problem
with the standard inclusion of correlated normalisation uncertainties in
the covariance matrix while the problem seems to be in the data set.

One may ask what is the meaning of a value of :math:`\overline{\mu}`
that can be smaller than both :math:`x_1` and :math:`x_2` or in general
far from both of them. To answer this question I first observe that,
according to Eq. :eq:`eq:offcv`, :math:`\overline{\mu}`
tends to zero as the uncorrelated uncertainties :math:`\sigma_1` and
:math:`\sigma_2` tend to zero. This can be seen more explicitly by
setting :math:`\sigma_1=\sigma_2=\varepsilon` and expanding the
:math:`\chi^2` around :math:`\varepsilon = 0`:

.. math::
   \chi^2\rightarrow
   \frac{\mu^2(x_1-x_2)^2}{\varepsilon^2(x_1^2+x_2^2)}+\frac{\left[\mu(x_1+x_2)-(x_1^2+x_2^2)\right]^2}{\delta_f^2
     (x_1^2+ x_2^2)^2}+\mathcal{O}(\varepsilon^2)\,.
   :label: eq:chi2exp

\ Therefore, if :math:`x_1\neq x_2`, for
:math:`\varepsilon\rightarrow 0` the dominant term is the first one and
the minimum of the :math:`\chi^2` in this limit is clearly at
:math:`\mu=0`. If :math:`x_1 = x_2=x`, instead, the first term drops and
the best fit value of the parameter :math:`\mu` is determined by the
second term that has its minimum at :math:`\mu=x`.

This can be interpreted as follows: as :math:`\sigma_1` and
:math:`\sigma_2` tend to zero, the covariance matrix tends to become
singular signifying that :math:`x_1` and :math:`x_2` are no longer
independent. However, if :math:`x_1=x_2=x` the singularity cancels. This
is due to the fact that, without uncorrelated uncertainties, the number
of points of the data set effectively reduces to one. As a matter of
fact, in this case :math:`\varepsilon` can safely be set to zero and the
:math:`\chi^2` becomes:

.. math:: \chi^2= \left(\frac{\mu - x}{\delta_f  x}\right)^2\,,

\ that is evidently the :math:`\chi^2` associated to one single data
point with central value :math:`x` and uncertainty
:math:`\delta_f x`. [10]_ As a consequence, the best fit value of the
parameters :math:`\mu` is simply given by :math:`x`, as one would
expect.

If :math:`x_1\neq x_2`, instead, the singularity of the covariance
matrix does not cancel and the first term in
Eq. :eq:`eq:chi2exp` dominates in the limit
:math:`\varepsilon\rightarrow0`. This is a situation in which two
totally dependent measurements of the same quantity have different
central values. It is thus to be expected that the :math:`\chi^2`
becomes arbitrarily large because the two data points are inconsistent
and thus impossible to be simultaneously fitted. In this situation, the
parameter :math:`\mu` attempts to minimise the :math:`\chi^2` by tending
to zero. More specifically, when correlations dominate over uncorrelated
uncertainties, the :math:`\chi^2` is proportional to (the square of) the
cross difference :math:`\delta_fx_1y_2-\delta_fx_2y_1`, [11]_ with
:math:`y_i = \mu - x_i`. Therefore, if :math:`x_1\neq x_2`, this
difference is minimal when :math:`y_i\simeq x_i` which is equivalent to
:math:`\mu\simeq 0`. It is also interesting to observe that the first
term in Eq. :eq:`eq:chi2exp` does not depend of
:math:`\delta_f`. Therefore, the size of the :math:`\chi^2` is entirely
driven by the difference :math:`x_1-x_2`, no matter the precise value of
normalisation uncertainty.

In conclusion, it appears that data sets dominated by normalisation (or
equivalently widely spread correlated) uncertainties carry an intrinsic
degree of inconsistency that, when attempting to fit them, manifests
into the DAB. Therefore, from this point of view, the DAB is not a flaw
of the fitting procedure but rather an indication that the data set is
intrinsically inconsistent. It thus appears counterproductive, if not
wrong, to adjust the definition of the :math:`\chi^2` to accommodate an
inconsistency of the data set.

.. container::
   :name: the-effect-of-the-t_0-prescription

   .. rubric:: The effect of the :math:`t_0` prescription
      :name: the-effect-of-the-t_0-prescription

In this section, I consider the effect of the iterative :math:`t_0`
prescription on the best fit value :math:`\overline{\mu}`. The
:math:`t_0` prescription essentially replaces the normalisation
uncertainties :math:`\delta_fx_1` and :math:`\delta_fx_2` with one
single correlated uncertainties :math:`\delta_f t_0`, where :math:`t_0`
is a “guess” for :math:`\mu`. The net effect of this replacement is that
the terms :math:`\delta_fx_1` and :math:`\delta_fx_2` in the r.h.s. of
both Eqs. :eq:`eq:normcv`
and :eq:`eq:normcvunc` are replaced by
:math:`\delta_ft_0` and thus cancel, reducing to
Eqs.:eq:`eq:offcv` and :eq:`eq:offcvunc`.
But this should not be surprising because the :math:`t_0` prescription
is essentially replacing two correlated uncertainties that are different
in value with two correlated uncertainties that are instead *equal*.
This makes the normalisation case artificially equal to the offset case
in which the DAB is not present.

Now, what happens if one iterates the :math:`t_0`? Let us indicate with
:math:`\overline{\mu}_i` and :math:`\delta\mu_i^2` best fit value and
variance at the :math:`i`-th iteration. As explained above, the zero-th
iteration gives:

.. math:: \overline{\mu}_0 =\frac{\frac{x_1}{\sigma_1^2}+\frac{x_2}{\sigma_2^2}}{\frac{1}{\sigma_1^2}+\frac{1}{\sigma_2^2}}\,,

and:

.. math:: \delta\mu_0^2 = \frac{1+\frac{\delta_f^2 t_0^2}{\sigma_1^2}+\frac{\delta_f^2 t_0^2}{\sigma_2^2}}{\frac{1}{\sigma_1^2}+ \frac{1}{\sigma_2^2}}\,.

Since :math:`\overline{\mu}_0` does not depend on :math:`t_0`, all
following iterations will give the exact same result, *i.e.*
:math:`\overline{\mu}_i=\overline{\mu}_0` :math:`\forall i>0`. The
stabilisation happens only one step further for the uncertainty for
which one finds :math:`\delta\mu_i^2=\delta\mu_1^2` :math:`\forall i>1`,
with:

.. math:: \delta\mu_1^2 = \frac{1+\frac{\delta_f^2 \overline{\mu}_0^2}{\sigma_1^2}+\frac{\delta_f^2 \overline{\mu}_0^2}{\sigma_2^2}}{\frac{1}{\sigma_1^2}+ \frac{1}{\sigma_2^2}}=\frac{1}{\frac{1}{\sigma_1^2}+ \frac{1}{\sigma_2^2}}+\delta_f^2\overline{\mu}_0^2\,.

As above, I set :math:`\sigma_1=\sigma_2=\varepsilon` getting:

.. math:: \overline{\mu}_0 =\frac{x_1+x_2}{2}\,,

and:

.. math:: \delta\mu_1^2 = \frac{\varepsilon^2}{2}+\frac{\delta_f^2  (x_1+x_2)^2}{4}\,.

Now, I plug :math:`\overline{\mu}_0` into the original definition of the
:math:`\chi^2`. The result is:

.. math:: \chi^2=\left(\frac{x_1- x_2}{2}\right)^2\left[\frac{1}{\varepsilon^2+\delta_f^2(x_1^2+x_2^2)}+\frac{1}{\varepsilon^2}\right]\,.

**References**

**References**

.. container:: references csl-bib-body hanging-indent
   :name: refs

   .. container:: csl-entry
      :name: ref-Ball:2008by

      Ball, Richard D., Luigi Del Debbio, Stefano Forte, Alberto
      Guffanti, Jose I. Latorre, Andrea Piccione, Juan Rojo, and Maria
      Ubiali. 2009. “A Determination of parton distributions with
      faithful uncertainty estimation.” *Nucl. Phys. B* 809: 1–63.
      https://doi.org/10.1016/j.nuclphysb.2008.09.037.

   .. container:: csl-entry
      :name: ref-Ball:2009qv

      Ball, Richard D., Luigi Del Debbio, Stefano Forte, Alberto
      Guffanti, Jose I. Latorre, Juan Rojo, and Maria Ubiali. 2010.
      “Fitting Parton Distribution Data with Multiplicative
      Normalization Uncertainties.” *JHEP* 05: 075.
      https://doi.org/10.1007/JHEP05(2010)075.

   .. container:: csl-entry
      :name: ref-Ball:2012wy

      Ball, Richard D., and others. 2013. “Parton Distribution
      Benchmarking with LHC Data.” *JHEP* 04: 125.
      https://doi.org/10.1007/JHEP04(2013)125.

   .. container:: csl-entry
      :name: ref-Ball:2014uwa

      ———. 2015. “Parton distributions for the LHC Run II.” *JHEP* 04:
      040. https://doi.org/10.1007/JHEP04(2015)040.

   .. container:: csl-entry
      :name: ref-wiki:xxx

      contributors, Wikipedia. 2020. “Chi-square distribution —
      Wikipedia, The Free Encyclopedia.”
      `url{https://en.wikipedia.org/w/index.php?title=Chi-square_distribution&oldid=979585636} <url{https://en.wikipedia.org/w/index.php?title=Chi-square_distribution&oldid=979585636}>`__.

   .. container:: csl-entry
      :name: ref-DAgostini:1993arp

      D’Agostini, G. 1994. “On the use of the covariance matrix to fit
      correlated data.” *Nucl. Instrum. Meth. A* 346: 306–11.
      https://doi.org/10.1016/0168-9002(94)90719-6.

.. [1]
   There could be more than one uncorrelated systematic uncertainty. In
   this case, :math:`\sigma_{i,\rm unc}` is just the square root of the
   sum in quadrature of all the uncorrelated systematic uncertainties.

.. [2]
   Note that this redefinition does not change the nature of the
   uncertainties, additive uncertainties remain additive as well as
   multiplicative uncertainties remain multiplicative.

.. [3]
   Note that one needs to include *also* an uncorrelated uncertainty
   because otherwise the covariance matrix would be singular. This a
   symptom of the fact that if no uncorrelated uncertainties are
   present, the two points would be totally dependent on each other.
   This observation will be relevant in the following.

.. [4]
   The contravariant counterpart of :math:`\mathbf{V}^{-1}` is clearly
   :math:`\mathbf{V}`, given that
   :math:`\mathbf{V}^{-1}\cdot\mathbf{V}=\mathbb{I}`.

.. [5]
   Notice that :math:`\mathbf{P}` is indeed idempotent, *i.e.*
   :math:`\mathbf{P}^2=\mathbf{P}`, as required.

.. [6]
   :math:`\mathbf{Q}` is such that :math:`\mathbf{Q}^2=\mathbf{Q}`,
   :math:`\mathbf{P}\cdot\mathbf{Q}=0`, and
   :math:`\mathbf{P}+\mathbf{Q}=\mathbb{I}`.

.. [7]
   As an interesting aside, the difference between the following
   equalities:

   .. math::
      \chi^2 = \left|\mathbf{L}^{-1}\cdot
          \mathbf{y}\right|^2\quad\mbox{and}\quad \chi^2= \left\|\mathbf{y}\right\|^2\,,

   can be paralleled with the norm of a vector in special relativity in
   the Euclidean space in the first case (where the time component is
   multiplied by the imaginary unit :math:`i` and the space components
   are left unchanged) and in Minkowski space in the second case.
   Indeed, the Cholesky decomposition of the metric tensor
   :math:`g_{\mu\nu}=\mbox{diag}(-1,1,1,1)` is
   :math:`L_{\mu\nu}=\mbox{diag}(i,1,1,1)`.

.. [8]
   The relevant identities are:

   .. math::
      \langle z_i^{2n-1}\rangle = 0\quad\mbox{and}\quad \langle
          z_i^{2n}\rangle= (2n-1)!!\quad n=1,2,3,\dots\,,

   \ \ where :math:`(2n-1)!!=1\cdot3\cdot\dots\cdot(2n-1)` is the odd
   semi-factorial. When considering combinations of two random variables
   the result is:

   .. math::
      \langle z_i^n z_j^m\rangle = \left\{
        \begin{array}{ll}
        \langle z_i^{n+m}\rangle&\quad i = j\\
        \langle z_i^n\rangle\langle z_j^m\rangle&\quad i \neq j
        \end{array}
        \right.

.. [9]
   Actually, the second derivative of the :math:`\chi^2` does not depend
   on :math:`\mu` and thus it is not necessary to set
   :math:`\mu=\overline{\mu}`.

.. [10]
   Note that in the case of one-point data set there is no meaning in
   talking about correlated uncertainties. Therefore, :math:`\delta_f x`
   plays the role of an uncorrelated uncertainty.

.. [11]
   In general, this can be seen as
   :math:`\sigma_{c,1}y_2-\sigma_{c,2}y_1`.

.. |image| image:: ../../latex/src/../../latex/src/plots/Chi2DistBABAR.pdf
   :width: 49.0%
.. |image1| image:: ../../latex/src/../../latex/src/plots/Chi2DistALEPH.pdf
   :width: 49.0%
