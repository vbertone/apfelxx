Ancillary files for "The four-loop non-singlet splitting functions in QCD"

Authors: Thomas Gehrmann, Andreas von Manteuffel, Vasily Sotnikov, and Tong-Zhi Yang

------------------------------------------------------------------------

File descriptions:

- GammaR.m provide the result for four-loop rapidity anomalous dimension,
i.e. eq. (15) in the paper.

- Gqq3NSp.m, Gqq3NSm.m, and Gqq3NSs.m provide the four-loop non-singlet
  anomalous dimensions gamma_ns^(3)+, gamma_ns^(3)-, and gamma_ns^(3)s,
  respectively.

- Pqq3NSp.m, Pqq3NSm.m, and Pqq3NSs.m provide the four-loop non-singlet
  splitting functions P_ns^(3)+, P_ns^(3)-, and P_ns^(3)s, respectively.

- Pqq3NSp_Fit.m, Pqq3NSm_Fit.m, and Pqq3NSs_Fit.m provide high-quality
  parametrizations of the four-loop non-singlet splitting functions
  P_ns^(3)+, P_ns^(3)-, and P_ns^(3)s, respectively. In these files, two
  labels separate the endpoint contributions from the regular part:
  flag[plusDelta] denotes the contributions from delta(1-x) and the plus
  distributions, while flag[NLP] denotes the next-to-leading power (NLP)
  corrections in the limits x->1 and x->0. The coefficients of both labels
  are given in analytic form.

- Pqq3NSpmXto1.m provides the result for P_ns^(3)+ and P_ns^(3)- in the
  x->1 limit up to next-to-leading power, i.e. eq. (12) in the paper. 
  The four-loop virtual anomalous dimension B_4 in eq. (13) can be extracted 
  by taking the coefficient of delta[1-x] of the above result.

- Pqq3NSpXto0.m, Pqq3NSmXto0.m, and Pqq3NSsXto0.m provide the results for
  P_ns^(3)+, P_ns^(3)-, and P_ns^(3)s, respectively, in the x->0 limit up
  to next-to-leading power, i.e. eqs. (16), (17), (18) in the paper.

- PTqq3NSp.m, PTqq3NSm.m and PTqq3NSs.m provide the four-loop non-singlet
 time-like splitting functions PT_ns^(3)+, PT_ns^(3)-, and PT_ns^(3)s, respectively.

------------------------------------------------------------------------

Color structures:

The files above contain results decomposed into 15 color structures:

  color[ca^3*cf], color[ca^2*cf^2], color[ca*cf^3], color[cf^4],
  color[d4FA/nc], color[ca^2*cf*nf], color[ca*cf^2*nf], color[cf^3*nf],
  color[ca*d33c*nf/nc], color[cf*d33c*nf/nc], color[d4FF*nf/nc],
  color[ca*cf*nf^2], color[cf^2*nf^2], color[d33c*nf^2/nc],
  color[cf*nf^3],

where for SU(3):

  cf      = (nc^2-1)/(2*nc)                            = 4/3
  ca      = nc                                         = 3
  d33c/nc = (nc^2-1)*(nc^2-4)/(16*nc^2)                = 5/18
  d4FF/nc = (nc^2-1)*(nc^4-6*nc^2+18)/(96*nc^3)        = 5/36
  d4FA/nc = (nc^2-1)*(nc^2+6)/(48*nc)                  = 5/2

and nf denotes the number of massless quark flavors.

------------------------------------------------------------------------

Notation:

- delta(1-x)      : the Dirac delta function;
- plusD[n, 1-x]   = (log(1-x)^n / (1-x))_+ : the plus distribution;
- zeta[n]         : Riemann zeta values;
- HPL[{...}, x]     : harmonic polylogarithms (HPLs);
- S[..., n]       : harmonic sums.
