//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/zeromasscoefficientfunctionspol_sl.h"
#include "apfel/specialfunctions.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________________
  G41ns::G41ns():
    Expression()
  {
  }
  double G41ns::Regular(double const& x) const
  {
    return 2 * CF * ( - ( 1 + x ) * log(1 - x) - ( 1 + pow(x, 2) ) * log(x) / ( 1 - x ) + 3 + 2 * x );
  }
  double G41ns::Singular(double const& x) const
  {
    return 2 * CF * ( 2 * log(1 - x) - 3 / 2. ) / ( 1 - x );
  }
  double G41ns::Local(double const& x) const
  {
    return 2 * CF * ( pow(log(1 - x), 2) - 3 * log(1 - x) / 2 - ( 2 * zeta2 + 9 / 2. ) );
  }

  //_________________________________________________________________________________
  GL1ns::GL1ns():
    Expression()
  {
  }
  double GL1ns::Regular(double const& x) const
  {
    return 4 * CF * x;
  }

  //_________________________________________________________________________________
  G11ns::G11ns():
    Expression()
  {
  }
  double G11ns::Regular(double const& x) const
  {
    return 2 * CF * ( - ( 1 + x ) * log(1 - x) - ( 1 + pow(x, 2) ) * log(x) / ( 1 - x ) + 2 + x );
  }
  double G11ns::Singular(double const& x) const
  {
    return 2 * CF * ( 2 * log(1 - x) - 3 / 2. ) / ( 1 - x );
  }
  double G11ns::Local(double const& x) const
  {
    return 2 * CF * ( pow(log(1 - x), 2) - 3 * log(1 - x) / 2 - ( 2 * zeta2 + 9 / 2. ) );
  }

  //_________________________________________________________________________________
  G11g::G11g():
    Expression()
  {
  }
  double G11g::Regular(double const& x) const
  {
    return 4 * TR * ( ( 2 * x - 1 ) * log( ( 1 - x ) / x ) - 4 * x + 3 );
  }

  //_________________________________________________________________________________
  G12nsp::G12nsp(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double G12nsp::Regular(double const& z) const
  {
    const double z2     = z * z;
    const double S121mz = apfel::wgplg(1, 2, 1 - z);
    const double S12mz  = apfel::wgplg(1, 2, - z);
    const double Li31mz = apfel::wgplg(2, 1, 1 - z);
    const double Li3mz  = apfel::wgplg(2, 1, - z);
    const double Li3r   = apfel::wgplg(2, 1, ( 1 - z ) / ( 1 + z ));
    const double Li3mr  = apfel::wgplg(2, 1, - ( 1 - z ) / ( 1 + z ));
    const double Li21mz = dilog(1 - z);
    const double Li2mz  = dilog(- z);
    const double ln1mz  = log(1 - z);
    const double ln1mz2 = ln1mz * ln1mz;
    const double ln1mz3 = ln1mz * ln1mz2;
    const double lnz    = log(z);
    const double lnz2   = lnz * lnz;
    const double lnz3   = lnz * lnz2;
    const double ln1pz  = log(1 + z);
    const double ln1pz2 = ln1pz * ln1pz;
    const double c2NSp = CF * CF * ( ( 1 + z2 ) / ( 1 - z ) * ( - 12 * S121mz + 12 * Li31mz + 48 * Li3mz + 36 * zeta3
                                                                - 6 * Li21mz - 24 * lnz * Li2mz + 24 * zeta2 * lnz - 4 * ln1mz * Li21mz
                                                                + 12 * lnz2 * ln1mz - 14 * lnz * ln1mz2 - 4 * lnz3 / 3. - 3 * lnz2 / 2.
                                                                + 12 * lnz * ln1mz + 61 * lnz / 2. )
                                     + ( 1 + z ) * ( - 4 * Li31mz + 4 * ln1mz * Li21mz - 4 * ln1mz3 - 4 * lnz * Li21mz - 4 * zeta2 * lnz
                                                     + 2 * lnz * ln1mz2 - 4 * lnz2 * ln1mz + 5 * lnz3 / 3. + 4 * zeta3 )
                                     + ( 1 - z ) * ( 8 * S121mz - 16 * Li3mz + 8 * lnz * Li2mz )
                                     - 8 * ( 1 + z + z2 + 1. / z ) * ( Li2mz + lnz * ln1pz ) - 2 * ( 5 + 13 * z ) * Li21mz
                                     - 4 * ( 7 + 2 * z2 + 7 * z ) * zeta2 + 2 * ( 5 + 7 * z ) * ln1mz2 + 8 * ( 1 + 3 * z ) * zeta2 * ln1mz
                                     + 2 * ( 5 + 3 * z ) * ln1mz + ( 29 / 2. + 4 * z2 + 41 * z / 2. ) * lnz2 - 16 * ( 1 + 2 * z ) * lnz * ln1mz
                                     + 3 * ( 3 - z ) * lnz / 2. - 41 - 10 * z )
                         + CA * CF * ( ( 1 + z2 ) / ( 1 - z ) * ( 12 * S121mz - 12 * Li31mz - 24 * Li3mz - 18 * zeta3 + 22 * Li21mz / 3. + 12 * lnz * Li2mz
                                                                  + 4 * ln1mz * Li21mz + 4 * lnz * Li21mz + 2 * lnz2 * ln1mz
                                                                  + 44 * lnz * ln1mz / 3. - lnz3 - 55 * lnz2 / 6. - 239 * lnz / 6. )
                                       + ( 1 + z ) * ( 4 * Li21mz + 11 * ln1mz2 / 3. + 4 * lnz * ln1mz - 20 * zeta3 )
                                       + ( 1 - z ) * ( - 4 * S121mz + 8 * Li3mz - 4 * lnz * Li2mz )
                                       + ( 8 * zeta2 - 134 / 9. - 314 * z / 9. ) * ln1mz + 4 * ( 1 + z + z2 + 1. / z ) * ( Li2mz + lnz * ln1pz )
                                       + 4 * zeta2 * ( 3 * z2 - 4 * z - 4 ) / 3. - 2 * ( 2 + z2 + 2 * z ) * lnz2 + ( 157 * z - 47 ) * lnz / 6.
                                       + 5 * ( 145 + 364 * z ) / 27. )
                         + _nf * CF * ( ( 1 + z2 ) / ( 1 - z ) * ( - 4 * Li21mz / 3. - 8 * lnz * ln1mz / 3. + 5 * lnz2 / 3. + 19 * lnz / 3. )
                                        + ( 1 + z ) * ( 4 * zeta2 / 3. - 2 * ln1mz2 / 3. ) + 4 * ( 5 + 14 * z ) * ln1mz / 9.
                                        + ( 1 + 11 * z ) * lnz / 3. - 2 * ( 58 + 151 * z ) / 27. );
    const double c2NSm = ( pow(CF, 2) - CA * CF / 2. )
                         * ( ( 1 + z2 ) / ( 1 + z ) * ( 16 * Li3r - 16 * Li3mr + 8 * S121mz - 16 * Li31mz - 16 * S12mz + 8 * Li3mz
                                                        - 16 * ln1mz * Li2mz - 16 * ln1pz * Li2mz + 8 * lnz * Li21mz
                                                        + 16 * lnz * Li2mz - 16 * lnz * ln1pz * ln1mz + 20 * lnz2 * ln1pz
                                                        + 4 * lnz2 * ln1mz - 8 * lnz * ln1pz2 - 8 * zeta2 * ln1mz
                                                        - 8 * zeta2 * ln1pz - 2 * lnz3 - 8 * lnz + 8 * zeta3 )
                             + ( 1 + z ) * ( 16 * S12mz - 8 * Li3mz + 16 * ln1pz * Li2mz
                                             + 8 * zeta2 * ln1pz + 8 * lnz * ln1pz2 - 4 * lnz2 * ln1pz + 8 * Li21mz
                                             + 8 * lnz * ln1mz - 8 * zeta3 )
                             + ( 1 - z ) * ( 16 * ln1mz + 30 ) + 8 * ( z2 + 1. / z ) * ( Li2mz + lnz * ln1pz )
                             - 4 * ( 2 + z2 + z ) * lnz2 + 4 * ( 1 + 2 * z2 - z ) * zeta2 + ( 6 + 38 * z ) * lnz );
    return c2NSp + c2NSm;
  }
  double G12nsp::Singular(double const& z) const
  {
    const double ln1mz  = log(1 - z);
    const double ln1mz2 = ln1mz * ln1mz;
    const double ln1mz3 = ln1mz * ln1mz2;
    return CF * CF * ( 8 * ln1mz3 / ( 1 - z ) - 18 * ln1mz2 / ( 1 - z )
                       - ( 32 * zeta2 + 27 ) * ln1mz / ( 1 - z )
                       + ( - 8 * zeta3 + 36 * zeta2 + 51 / 2. ) / ( 1 - z ) )
           + CA * CF * ( - 22 * ln1mz2 / ( 1 - z ) / 3. + ( - 8 * zeta2 + 367 / 9. ) * ln1mz / ( 1 - z )
                         + ( 40 * zeta3 + 44 * zeta2 / 3. - 3155 / 54. ) / ( 1 - z ) )
           + _nf * CF * ( 4 * ln1mz2 / ( 1 - z ) / 3. - 58 * ln1mz / ( 1 - z ) / 9.
                          + ( - 8 * zeta2 / 3. + 247 / 27. ) / ( 1 - z ) );
  }
  double G12nsp::Local(double const& z) const
  {
    const double ln1mz  = log(1 - z);
    const double ln1mz2 = ln1mz * ln1mz;
    const double ln1mz3 = ln1mz * ln1mz2;
    const double ln1mz4 = ln1mz * ln1mz3;
    return CF * CF * ( 2 * ln1mz4 - 6 * ln1mz3 - ( 32 * zeta2 + 27 ) * ln1mz2 / 2. + ( - 8 * zeta3 + 36 * zeta2 + 51 / 2. ) * ln1mz
                       + 6 * pow(zeta2, 2) - 78 * zeta3 + 69 * zeta2 + 331 / 8. )
           + CA * CF * ( - 22 * ln1mz3 / 9. + ( - 4 * zeta2 + 367 / 18. ) * ln1mz2 + ( 40 * zeta3 + 44 * zeta2 / 3. - 3155 / 54. ) * ln1mz
                         + 71 * pow(zeta2, 2) / 5. + 140 * zeta3 / 3. - 251 * zeta2 / 3 - 5465 / 72. )
           + _nf * CF * ( 4 * ln1mz3 / 9. - 29 * ln1mz2 / 9. + ( - 8 * zeta2 / 3. + 247 / 27. ) * ln1mz
                          + 4 * zeta3 / 3. + 38 * zeta2 / 3. + 457 / 36. );
  }

  //_________________________________________________________________________________
  G12ps::G12ps():
    Expression()
  {
  }
  double G12ps::Regular(double const& z) const
  {
    const double z2     = z * z;
    const double Li31mz = apfel::wgplg(2, 1, 1 - z);
    const double Li21mz = dilog(1 - z);
    const double Li2mz  = dilog(- z);
    const double ln1mz  = log(1 - z);
    const double ln1mz2 = ln1mz * ln1mz;
    const double lnz    = log(z);
    const double lnz2   = lnz * lnz;
    const double lnz3   = lnz * lnz2;
    const double ln1pz  = log(1 + z);
    return CF * TR * ( ( 1 + z ) * ( - 16 * Li31mz + 16 * ln1mz * Li21mz - 16 * lnz * Li21mz - 16 * zeta2 * lnz
                                     + 8 * lnz * ln1mz2 - 16 * lnz2 * ln1mz * 20 * lnz3 / 3. )
                       // Erratum of 2007
                       + ( 1 - z ) * ( 20 * ln1mz2 - 88 * ln1mz + 760 / 3. )
                       - 32 * ( 1 + z2 / 3. + z + 1 / z / 3. ) * ( Li2mz + lnz * ln1pz )
                       + ( 50 + 16 * z2 / 3. - 10 * z ) * lnz2 - 32 * ( 2 - z ) * lnz * ln1mz
                       + 4 * ( 119 - 13 * z ) * lnz / 3. - ( 72 + 32 * z2 / 3. - 40 * z ) * zeta2
                       - 8 * ( 3 + z ) * Li21mz
                       //
                       /* Amended part
                       + ( 1 - z ) * ( 20 * ln1mz2 - 88 * ln1mz + 808 / 3. )
                       - 32 * ( 1 + z2 / 3. + z + 1 / z / 3. ) * ( Li2mz + lnz * ln1pz )
                       + ( 58 + 16 * z2 / 3. - 6 * z ) * lnz2 - 32 * ( 2 - z ) * lnz * ln1mz
                       + 4 * ( 137 - 19 * z ) * lnz / 3. - ( 72 + 32 * z2 / 3. - 40 * z ) * zeta2
                       - 8 * ( 3 + z ) * Li21mz
                       */
                     );
  }

  //_________________________________________________________________________________
  G12g::G12g():
    Expression()
  {
  }
  double G12g::Regular(double const& z) const
  {
    const double z2     = z * z;
    const double S121mz = apfel::wgplg(1, 2, 1 - z);
    const double S12mz  = apfel::wgplg(1, 2, - z);
    const double Li31mz = apfel::wgplg(2, 1, 1 - z);
    const double Li3mz  = apfel::wgplg(2, 1, - z);
    const double Li3r   = apfel::wgplg(2, 1, ( 1 - z ) / ( 1 + z ));
    const double Li3mr  = apfel::wgplg(2, 1, - ( 1 - z ) / ( 1 + z ));
    const double Li21mz = dilog(1 - z);
    const double Li2mz  = dilog(- z);
    const double ln1mz  = log(1 - z);
    const double ln1mz2 = ln1mz * ln1mz;
    const double ln1mz3 = ln1mz * ln1mz2;
    const double lnz    = log(z);
    const double lnz2   = lnz * lnz;
    const double lnz3   = lnz * lnz2;
    const double ln1pz  = log(1 + z);
    const double ln1pz2 = ln1pz * ln1pz;
    return CF * TR * ( ( 1 - 2 * z ) * ( 32 * Li31mz - 16 * ln1mz * Li21mz - 8 * lnz * Li21mz - 24 * zeta2 * lnz
                                         - 20 * ln1mz3 / 3. + 16 * lnz * ln1mz2 - 16 * lnz2 * ln1mz + 10 * lnz3 / 3. )
                       - 16 * ( 1 + z2 + 2 * z ) * ( 4 * S12mz + 4 * ln1pz * Li2mz + 2 * lnz * ln1pz2 - lnz2 * ln1pz + 2 * zeta2 * ln1pz )
                       - 32 * ( 1 + z2 - 6 * z ) * Li3mz + 8 * ( 1 + 4 * z2 - 2 * z ) * S121mz
                       + 16 * ( 13 * z2 + 12 * z + 4. / z ) * ( Li2mz + lnz * ln1pz ) / 3.
                       // Erratum of 2007
                       + 4 * ( 5 - 12 * z ) * Li21mz + 32 * ( 1 + z2 - 2 * z ) * lnz * Li2mz
                       + ( 123 - 104 * z2 - 48 * z ) * lnz2 / 3. - ( 88 - 96 * z ) * lnz * ln1mz
                       + 6 * ( 9 - 12 * z ) * ln1mz2 - 32 * zeta2 * z2 * ln1mz - 4 * ( 31 - 4 * z2 - 26 * z ) * ln1mz
                       + ( 416 - 48 * z2 - 274 * z ) * lnz / 3. - 8 * ( 5 - 4 * z2 - 26 * z ) * zeta3
                       - 4 * ( 81 - 52 * z2 - 108 * z ) * zeta2 / 3. + 2 * ( 233 - 239 * z ) / 3.
                       //
                       /* Amended part
                       - 4 * ( 3 + 28 * z ) * Li21mz + 32 * ( 1 + z2 - 2 * z ) * lnz * Li2mz
                       + ( 171 - 104 * z2 + 48 * z ) * lnz2 / 3. - ( 120 - 32 * z ) * lnz * ln1mz
                       + 6 * ( 9 - 12 * z ) * ln1mz2 - 32 * zeta2 * z2 * ln1mz - 4 * ( 55 - 4 * z2 - 50 * z ) * ln1mz
                       + ( 800 - 48 * z2 - 82 * z ) * lnz / 3. - 8 * ( 5 - 4 * z2 - 26 * z ) * zeta3
                       - 4 * ( 81 - 52 * z2 - 108 * z ) * zeta2 / 3. + 2 * ( 641 - 647 * z ) / 3.
                       */
                     )
           + CA * TR * ( 16 * ( 1 + 2 * z ) * ( Li3r - Li3mr - ln1mz * Li2mz - lnz * Li21mz - lnz * ln1mz * ln1pz )
                         + 16 * ( 1 + z2 + 2 * z ) * ( 2 * S12mz + lnz * ln1pz2 + 2 * ln1pz * Li2mz + zeta2 * ln1pz )
                         + 8 * ( 1 - 2 * z2 + 2 * z ) * S121mz - 8 * ( 9 + 2 * z ) * Li31mz
                         - 8 * ( 1 - z2 + 2 * z ) * ( 2 * Li3mz - lnz2 * ln1pz - 2 * lnz * Li2mz )
                         - 16 * ( 2 + z ) * lnz2 * ln1mz + 24 * ( lnz * ln1mz2 - 2 * zeta2 * lnz - Li21mz )
                         - 16 * ( 6 + 11 * z2 + 12 * z + 2. / z ) * ( Li2mz + lnz * ln1pz ) / 3.
                         - 4 * ( 1 - 2 * z ) * ln1mz3 / 3. + 8 * ( 6 - 7 * z ) * ln1mz2
                         + 8 * ( 3 + 2 * z2 - 10 * z ) * zeta2 * ln1mz + 4 * ( 7 + 10 * z ) * lnz3 / 3.
                         + 2 * ( 135 + 44 * z2 - 48 * z ) * lnz2 / 3. - 8 * ( 17 - 16 * z ) * lnz * ln1mz
                         + 8 * ( 118 + 3 * z2 - 26 * z ) * lnz / 3. - 4 * ( 44  + 2 * z2 - 53 * z ) * ln1mz
                         - 4 * ( 3 + 4 * z2 + 10 * z ) * zeta3 - 16 * ( 27 + 11 * z2 - 24 * z ) * zeta2 / 3.
                         + 8 * ( 5 + 2 * z ) * ln1mz * Li21mz + 4 * ( 355 - 367 * z ) / 3. );
           // Erraturm of 1997
           + CF * TR * ( 16 * ( 1 + 2 * z ) * ( 2 * Li21mz + 2 * lnz * ln1mz - lnz2 )
                         + 96 * ( 1 - z ) * ln1mz - ( 144 + 64 * z ) * lnz - 304 * ( 1 - z ) );
    //
  }
}
