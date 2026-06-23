//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/kcs.h"
#include "apfel/betaqcd.h"
#include "apfel/gammak.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________
  double KCS00()
  {
    return 0;
  }

  //_________________________________________________________________________
  double KCS01()
  {
    return - gammaK0();
  }

  //_________________________________________________________________________
  double KCS10(int const& nf)
  {
    return CA * ( 28 * zeta3 - 808. / 27. ) + 224 * TR * nf / 27;
  }

  //_________________________________________________________________________
  double KCS11(int const& nf)
  {
    return - gammaK1(nf);
  }

  //_________________________________________________________________________
  double KCS12(int const& nf)
  {
    return - ( - 2 ) * beta0qcd(nf) * gammaK0() / 2;
  }

  //_________________________________________________________________________
  double KCS20(int const& nf)
  {
    return 2 * ( CA * CA * ( - 176 * zeta3 * zeta2 / 3 + 6392 * zeta2 / 81 + 12328 * zeta3 / 27
                             + 154 * zeta4 / 3 - 192 * zeta5 - 297029. / 729. ) / 2
                 + CA * TR * nf * ( - 824 * zeta2 / 81 - 904 * zeta3 / 27 + 20 * zeta4 / 3 + 62626. / 729. )
                 + 2 * TR * TR * nf * nf * ( - 32 * zeta3 / 9 - 1856. / 729. )
                 + CF * TR * nf * ( - 304 * zeta3 / 9 - 16 * zeta4 + 1711. / 27. ) );
  }

  //_________________________________________________________________________
  double KCS21(int const& nf)
  {
    return - gammaK2(nf);
  }

  //_________________________________________________________________________
  double KCS22(int const& nf)
  {
    const double bt0 = - 2 * beta0qcd(nf);
    const double b1 = beta1qcd(nf) / beta0qcd(nf);
    return - bt0 * ( 2 * gammaK1(nf) + b1 * gammaK0() ) / 2;
  }

  //_________________________________________________________________________
  double KCS23(int const& nf)
  {
    const double bt02 = pow(- 2 * beta0qcd(nf), 2);
    return - bt02 * gammaK0() / 3;
  }

  //_________________________________________________________________________
  double KCS30(int const& nf)
  {
    const int nf2       = nf * nf;
    const int nf3       = nf2 * nf;
    const double CA2    = CA * CA;
    const double CA3    = CA * CA2;
    const double CF2    = CF * CF;
    const double zeta32 = zeta3 * zeta3;
    // Quartic Casimirs normalised to the dimension of the
    // fundamental representation, d_F^{abcd} X^{abcd} / N_R, for SU(3)
    // (same values as in gammaK3).
    const double dFA4nc = 5. / 2.;
    const double dFF4nc = 5. / 36.;
    // The four genuinely-new four-loop constants, known only
    // numerically. Values from eq. (10) of
    // https://arxiv.org/pdf/2205.02242.
    const double b4C4AF    = - 998.0;    // b_{4,q,C_4^{AF}}
    const double b4C4FF    = - 143.6;    // b_{4,q,C_4^{FF}}
    const double b4nfCACF2 = - 455.247;  // b_{4,q,nf CF^2 CA}
    const double b4nfCF3   =   80.780;   // b_{4,q,nf CF^3}
    return
      CA3 * ( - b4C4AF / 12 - 13052. / 9. * zeta2 * zeta3 + 1376. / 3. * zeta2 * zeta5 + 389083. / 243. * zeta2
              - 10582. / 9. * zeta32 + 2114. / 3. * zeta3 * zeta4 + 600872. / 81. * zeta3 + 4144. / 9. * zeta4
              - 90962. / 27. * zeta5 - 31790. / 27. * zeta6 + 11071. / 6. * zeta7 - 28290079. / 4374. )
      + CA2 * nf * ( - b4C4FF / 24 - b4nfCACF2 - b4nfCF3 / 2 + 1040. / 3. * zeta2 * zeta3 - 91067. / 243. * zeta2
                     - 4292. / 9. * zeta32 - 123826. / 81. * zeta3 + 21812. / 27. * zeta4 - 8968. / 27. * zeta5
                     + 791. / 27. * zeta6 + 10761379. / 5832. )
      + CA * CF * nf * ( 2 * b4nfCACF2 - 4432. / 9. * zeta2 * zeta3 + 2561. / 27. * zeta2 + 3400. / 3. * zeta32
                         - 946. / 9. * zeta3 - 61108. / 27. * zeta4 + 10952. / 9. * zeta5 - 718 * zeta6 + 2149049. / 972. )
      + CF2 * nf * ( 2 * b4nfCF3 + 512. / 3. * zeta2 * zeta3 - 324 * zeta2 - 368 * zeta32 + 1120. / 9. * zeta3
                     + 334 * zeta4 - 3872. / 3. * zeta5 + 14668. / 9. * zeta6 - 27949. / 108. )
      + CA * nf2 * ( 112. / 9. * zeta2 * zeta3 + 3376. / 243. * zeta2 - 11128. / 81. * zeta3 + 80. / 9. * zeta4
                     + 736. / 9. * zeta5 - 898033. / 5832. )
      + CF * nf2 * ( 3464. / 27. * zeta3 + 80. / 3. * zeta4 + 16 * zeta5 - 110059. / 486. )
      + nf3 * ( 80. / 9. * zeta3 - 8. / 9. * zeta4 + 5216. / 2187. )
      + dFA4nc / CF * ( 2 * b4C4AF + 1792 * zeta2 * zeta3 - 1024 * zeta2 * zeta5 + 2176. / 3. * zeta2 + 3344. / 3. * zeta32
                        + 368 * zeta3 * zeta4 + 7808. / 9. * zeta3 - 112. / 3. * zeta4 + 1840. / 9. * zeta5
                        - 3476. / 9. * zeta6 - 3484 * zeta7 - 192 )
      + dFF4nc * nf / CF * ( 2 * b4C4FF - 128 * zeta2 * zeta3 - 4544. / 3. * zeta2 - 1216. / 3. * zeta32 + 5312. / 9. * zeta3
                             + 800. / 3. * zeta4 + 21760. / 9. * zeta5 - 1184. / 9. * zeta6 + 384 );
  }

  //_________________________________________________________________________
  double KCS31(int const& nf)
  {
    return 0. * nf;
  }

  //_________________________________________________________________________
  double KCS32(int const& nf)
  {
    const double bt0 = - 2 * beta0qcd(nf);
    const double b1 = beta1qcd(nf) / beta0qcd(nf);
    const double b2 = beta2qcd(nf) / beta0qcd(nf);
    return - bt0 * ( 3 * gammaK2(nf) + 2 * gammaK1(nf) * b1 + gammaK0() *b2 ) / 2;
  }

  //_________________________________________________________________________
  double KCS33(int const& nf)
  {
    const double bt02 = pow(- 2 * beta0qcd(nf), 2);
    const double b1 = beta1qcd(nf) / beta0qcd(nf);
    return - bt02 * ( 6 * gammaK1(nf) + 5 * gammaK0() * b1 ) / 6;
  }

  //_________________________________________________________________________
  double KCS34(int const& nf)
  {
    const double bt03 = pow(- 2 * beta0qcd(nf), 3);
    return - bt03 * gammaK0() / 4;
  }

  //_________________________________________________________________________
  double KCS30g(int const& nf)
  {
    const int nf2       = nf * nf;
    const int nf3       = nf2 * nf;
    const double CA2    = CA * CA;
    const double CA3    = CA * CA2;
    const double CF2    = CF * CF;
    const double zeta32 = zeta3 * zeta3;
    // Quartic Casimirs normalised to the dimension of the adjoint
    // representation, d_X^{abcd} d_A^{abcd} / N_A, for SU(3).
    const double dAANA = 135. / 8.;  // d_A^{abcd} d_A^{abcd} / N_A
    const double dFANA = 15. / 16.;  // d_F^{abcd} d_A^{abcd} / N_A
    // The four genuinely-new four-loop constants, known only
    // numerically. Values from eq. (10) of
    // https://arxiv.org/pdf/2205.02242. By generalised Casimir scaling
    // the gluon result shares the same constants and bracket
    // coefficients as the quark one (KCS30); only the colour
    // representation (C_F -> C_A) and the quartic Casimirs differ.
    const double b4C4AF    = - 998.0;    // b_{4,g,C_4^{AA}}
    const double b4C4FF    = - 143.6;    // b_{4,g,C_4^{FA}}
    const double b4nfCACF2 = - 455.247;  // b_{4,g,nf CF^2 CA}
    const double b4nfCF3   =   80.780;   // b_{4,g,nf CF^3}
    return
      CA3 * ( - b4C4AF / 12 - 13052. / 9. * zeta2 * zeta3 + 1376. / 3. * zeta2 * zeta5 + 389083. / 243. * zeta2
              - 10582. / 9. * zeta32 + 2114. / 3. * zeta3 * zeta4 + 600872. / 81. * zeta3 + 4144. / 9. * zeta4
              - 90962. / 27. * zeta5 - 31790. / 27. * zeta6 + 11071. / 6. * zeta7 - 28290079. / 4374. )
      + CA2 * nf * ( - b4C4FF / 24 - b4nfCACF2 - b4nfCF3 / 2 + 1040. / 3. * zeta2 * zeta3 - 91067. / 243. * zeta2
                     - 4292. / 9. * zeta32 - 123826. / 81. * zeta3 + 21812. / 27. * zeta4 - 8968. / 27. * zeta5
                     + 791. / 27. * zeta6 + 10761379. / 5832. )
      + CA * CF * nf * ( 2 * b4nfCACF2 - 4432. / 9. * zeta2 * zeta3 + 2561. / 27. * zeta2 + 3400. / 3. * zeta32
                         - 946. / 9. * zeta3 - 61108. / 27. * zeta4 + 10952. / 9. * zeta5 - 718 * zeta6 + 2149049. / 972. )
      + CF2 * nf * ( 2 * b4nfCF3 + 512. / 3. * zeta2 * zeta3 - 324 * zeta2 - 368 * zeta32 + 1120. / 9. * zeta3
                     + 334 * zeta4 - 3872. / 3. * zeta5 + 14668. / 9. * zeta6 - 27949. / 108. )
      + CA * nf2 * ( 112. / 9. * zeta2 * zeta3 + 3376. / 243. * zeta2 - 11128. / 81. * zeta3 + 80. / 9. * zeta4
                     + 736. / 9. * zeta5 - 898033. / 5832. )
      + CF * nf2 * ( 3464. / 27. * zeta3 + 80. / 3. * zeta4 + 16 * zeta5 - 110059. / 486. )
      + nf3 * ( 80. / 9. * zeta3 - 8. / 9. * zeta4 + 5216. / 2187. )
      + dAANA / CA * ( 2 * b4C4AF + 1792 * zeta2 * zeta3 - 1024 * zeta2 * zeta5 + 2176. / 3. * zeta2 + 3344. / 3. * zeta32
                       + 368 * zeta3 * zeta4 + 7808. / 9. * zeta3 - 112. / 3. * zeta4 + 1840. / 9. * zeta5
                       - 3476. / 9. * zeta6 - 3484 * zeta7 - 192 )
      + dFANA * nf / CA * ( 2 * b4C4FF - 128 * zeta2 * zeta3 - 4544. / 3. * zeta2 - 1216. / 3. * zeta32 + 5312. / 9. * zeta3
                            + 800. / 3. * zeta4 + 21760. / 9. * zeta5 - 1184. / 9. * zeta6 + 384 );
  }
}
