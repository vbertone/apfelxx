======================================================
Computation of the charged current structure functions
======================================================

.. contents::
   :depth: 3
..

The structure of the observables
================================

The structure in terms of PDFs of the charged current (CC) structure
functions is complicated by the mixing between down- and up-type quarks
provided by the CKM matrix. As a first step, we write the
:math:`\mathcal{O}(\alpha_s)` contribution to :math:`F=F_2,F_L` (we will
consider :math:`F_3` later) in a convenient way as:

.. math::
   :label: compactNu
   F^{\nu} = \sum_{U=u,c,t}\sum_{D=d,s,b}|V_{UD}|^2\left[C_\pm\left(D +\overline{U}\right) +2 C_gg\right]

and:

.. math::
   :label: compactNub
   F^{\overline{\nu}} = \sum_{U=u,c,t}\sum_{D=d,s,b}|V_{UD}|^2\left[C_\pm\left(\overline{D} +U\right) +2 C_gg\right]

where we have omitted the convolution symbol and an overall factor
:math:`2x`. At this order we don’t have to worry about whether
:math:`C_+` or :math:`C_-` has to be used because they coincide.
However, in the following it will appear naturally which one has be used
and where. One can combine the expressions above conveniently as:

.. math::
   F^{\pm} \equiv \frac{F^{\nu} \pm F^{\overline{\nu}}}2=\frac12
   \sum_{U=u,c,t}\sum_{D=d,s,b}|V_{UD}|^2\left[C_\pm\left(D^\pm
       \pm U^\pm\right) + P^{\pm}4C_gg\right]

where we have used the usual definition
:math:`q^{\pm} = q\pm\overline{q}` and defined the projector:

.. math:: P^\pm=\frac{1\pm1}2\,.

It should be noted that the subscript :math:`\pm` to the quark
coefficient function :math:`C_\pm` because is now associated to each of
:math:`F^\pm`.

Now we need to express these observables in terms of PDFs in the
evolution basis. The starting point is the relation:

.. math::
   :label: TranformationBella
   q_i^\pm = \sum_{j=1}^6M_{ij}d^\pm_j\,,

\ where :math:`d^\pm_j` belong to the QCD evolution basis, that is:
:math:`d^+_1=\Sigma`, :math:`d^+_2=-T_3`, :math:`d^+_3=T_8`,
:math:`d^+_4=T_{15}`, :math:`d^+_5=T_{24}`, and :math:`d^+_6=T_{35}` and
:math:`d^-_1=V`, :math:`d^-_2=-V_3`, :math:`d^-_3=V_8`,
:math:`d^-_4=V_{15}`, :math:`d^-_5=V_{24}`, and :math:`d^-_6=V_{35}`.
Note that here we are using the more “natural” ordering for the
distibutions :math:`q_i=\{d,u,s,c,b,t\}` rather than that where
:math:`u` comes before :math:`d`; this is the reason of the minus sign
in front of :math:`T_3` and :math:`V_3`. The trasformation matrix
:math:`M_{ij}` can be written as:

.. math::
   :label: TransDef
   \begin{array}{l}
   \displaystyle M_{ij}=\theta_{ji}\frac{1-\delta_{ij}j}{j(j-1)}\quad j\geq 2\,,\\
   \\
   \displaystyle M_{i1} = \frac{1}{6}\,,
   \end{array}

with :math:`\theta_{ji}=1` for :math:`j\geq i` and zero otherwise. In
addition, one can show that :math:`M_{ij}` is such that:

.. math::
   :label: eq:properties
   \sum_{j=1}^6M_{ij} = 0\,,\quad\mbox{and}\quad \sum_{i=1}^6M_{ij} = \delta_{1j}\,.

Using eq. :eq:`TranformationBella` we can make
the following identifications:

.. math::
   D^{\pm} = q_{2j-1}^\pm\quad\mbox{and}\quad U^{\pm} =
   q_{2j}^\pm\,,\quad j=1,2,3\,,

\ so that we can write:

.. math::
   F^\pm=
   \frac12\sum_{i=1}^3\sum_{j=1}^3|V_{2i,(2j-1)}|^2\left[C_\pm\left(q_{2j-1}^\pm
       \pm q_{2i}^\pm\right) + 4P^{\pm} C_g g\right]\,.

Using the definition of :math:`M_{ij}` in
eq. :eq:`TransDef`, we can rewrite :math:`F^{\pm}` in
terms of PDFs in the evolution basis as:

.. math::
   :label: eq:decompF2L
   F^\pm=
   \sum_{i=1}^3\sum_{j=1}^3|V_{2i,(2j-1)}|^2 F_{ij}^\pm\,,

with:

.. math::
   :label: F2Ldef
   F_{ij}^\pm=
   C_g 2P^\pm g
   +
   C_\pm^{\rm S} P^\pm \frac16 d_1^\pm
   + C_\pm\sum_{k=2}^6\frac{\theta_{k,2j-1}(1-\delta_{2j-1,k}k)\pm \theta_{k,2i}(1-\delta_{2i,k}k) }{2k(k-1)}d_k^\pm\,.

Eq. :eq:`F2Ldef` is valid only for :math:`F_2` and
:math:`F_3`. In order to obtain a similar equation also for :math:`F_3`,
one needs to change sign to the antiquark distributions, :math:`i.e.`
:math:`\overline{q}_i\rightarrow - \overline{q}_i`. In the QCD evolution
basis, this has the consequence of exchanging the :math:`T`-like
distributions with the :math:`V`-like ones, that is to say
:math:`d_k^+\leftrightarrow d_k^-`. It is the easy to see that:

.. math::
   F_3^\pm=
   \sum_{i=1}^3\sum_{j=1}^3|V_{2i,(2j-1)}|^2 F_{3,ij}^\pm\,,

\ with:

.. math::
   :label: F3def
   F_{3,ij}^\pm=
   C_g 2P^\pm g
   +
   C_\pm^{\rm S} P^\mp \frac16 d_1^\pm
   + C_\pm\sum_{k=2}^6\frac{\theta_{k,2j-1}(1-\delta_{2j-1,k}k)\mp \theta_{k,2i}(1-\delta_{2i,k}k) }{2k(k-1)}d_k^\pm\,.

It is now useful to consider the inclusive structure functions and
exploit the unitarity of the CKM matrix elements :math:`V_{UD}`:

.. math:: \sum_{i=1}^3|V_{2i,(2j-1)}|^2 = \sum_{j=1}^3|V_{2i,(2j-1)}|^2 = 1\quad\Rightarrow\quad \sum_{i=1}^3\sum_{j=1}^3|V_{2i,(2j-1)}|^2 = 3\,.

Summing over :math:`i` and :math:`j` in
eq. :eq:`eq:decompF2L` and using
eq. :eq:`F2Ldef`, one obtains:

.. math::
   F^\pm=
   C_g 6P^\pm g
   +
   C_\pm^{\rm S} P^\pm \frac12 d_1^\pm
   + \frac12 C_\pm\sum_{k=2}^6 d_k^\pm\sum_{l=1}^6(\pm 1)^{l+1}M_{lk}\,.

Considering separately :math:`F^+` and :math:`F^-` and using
eq. :eq:`eq:properties`, one finds:

.. math:: F^+= C_g 6g + C_+^{\rm S} \frac12 d_1^+

and:

.. math:: F^-= \frac12 C_-\sum_{k=2}^6\left[\frac{P^+}{k-1}-\frac{P^-}{k}\right]d_k^-\,,

with the even/odd projectors defined as:

.. math:: P_k^{\pm} = \frac{1\pm(-1)^k}{2}\,.

It should be pointed out that such simple expressions (independent of
the CMK matrix elements) is achievable only if it is possible to
factorize the non-singlet coefficient functions as implicitly done in
eqs. :eq:`F2Ldef` and :eq:`F3def`. In fact, this
is possible only in the ZM case in which the coefficient functions of
each PDF combination is the same.

For :math:`F_3` we find:

.. math:: F_3^+= C_g 6g + C_-^{\rm S} \frac12 d_1^-

\ and:

.. math:: F_3^-= \frac12 C_+\sum_{k=2}^6\left[\frac{P_k^+}{k-1}-\frac{P_k^-}{k}\right]d_k^+\,.
