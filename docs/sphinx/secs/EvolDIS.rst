=====================================
Combining evolution and DIS operators
=====================================

.. contents::
   :depth: 3
..

The structure of the observables
================================

In all cases the inclusive DIS structure functions are conveniently
expressed in terms of PDF combinations in the so-called physical basis
:math:`\{q_i^\pm\}`, with :math:`q_i^\mp = q_i\pm\overline{q}_i`, where
:math:`q_i` and :math:`\overline{q}_i` are the PDFs of the :math:`i-th`
quark-flavour, with :math:`i=u,d,s,c,b,t`, and its antiflavour,
respectively. Schematically, a DIS structure function can be written as:

.. math:: F = C_g g + \sum_{i}\left(C_i^+q_i^++C_i^-q_i^-\right)\,,

\ being :math:`C` the appropriate coefficient functions. Conversely, the
evolution of PDFs is usually compute in the so-called QCD evolution
basis :math:`\{d_i^\pm\}`, with :math:`d^+_1=\Sigma`,
:math:`d^+_2=-T_3`, :math:`d^+_3=T_8`, :math:`d^+_4=T_{15}`,
:math:`d^+_5=T_{24}`, and :math:`d^+_6=T_{35}` and :math:`d^-_1=V`,
:math:`d^-_2=-V_3`, :math:`d^-_3=V_8`, :math:`d^-_4=V_{15}`,
:math:`d^-_5=V_{24}`, and :math:`d^-_6=V_{35}`. The gluon remains
unchanged.

.. math:: g = \Gamma_{gg}g_0 + \Gamma_{gq}d_{1,0}^+

while:

.. math::
   \begin{array}{rcl}
   d_i^\pm&=&\displaystyle \theta_{i2}\theta(Q-m_i)\Gamma^{\pm}d_{i,0}^\pm+\theta(m_i-Q)
   \left\{\begin{array}{ll}
   \Gamma_{qq}d_{1,0}^++\Gamma_{qg}g_0&\quad\mbox{for }+\\
   \Gamma^vd_{1,0}^-&\quad\mbox{for }-
   \end{array}\right.\\
   \\
   &=&\left\{\begin{array}{l}
   \Gamma_{qq}d_{1,0}^++\Gamma_{qg}g_0\\
   \Gamma^vd_{1,0}^-
   \end{array}\right.+\theta(Q-m_i) \left[\theta_{i2}\Gamma^{\pm}d_{i,0}^\pm-
   \left\{\begin{array}{l}
   \Gamma_{qq}d_{1,0}^++\Gamma_{qg}g_0\\
   \Gamma^vd_{1,0}^-
   \end{array}\right.\right]\,.
   \end{array}

so that:

.. math::
   \begin{array}{rcl}
   \displaystyle  q_i^- &=& \displaystyle \delta_{i1}
     \Gamma^vd_{1,0}^-+\Gamma^{\pm}\sum_{j=2}^6
     \theta_{ji}\frac{1-\delta_{ij}j}{j(j-1)} 
     \theta(Q-m_j)d_{j,0}^\pm
   \end{array}

Therefore, we need to relate these two bases. This is done through the
linear transformation:

.. math::
   :label: TranformationBella
   q_i^\pm = \sum_{j=1}^6M_{ij}d^\pm_j\,,

\ where the trasformation matrix :math:`M_{ij}` can be written as:

.. math::
   :label: TransDef
   \begin{array}{l}
   \displaystyle M_{ij}=\theta_{ji}\frac{1-\delta_{ij}j}{j(j-1)}\quad j\geq 2\,,\\
   \\
   \displaystyle M_{i1} = \frac{1}{6}\,,
   \end{array}

with :math:`\theta_{ji}=1` for :math:`j\geq i` and zero otherwise.

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
eq. (`[eq:properties] <#eq:properties>`__), one finds:

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
