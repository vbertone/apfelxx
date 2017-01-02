//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/distribution.h"

namespace apfel
{
  //_________________________________________________________________________
  Distribution::Distribution()
  {
  }

  //_________________________________________________________________________
  void Distribution::BASIS2EVLN(array<double,14> const& basis, array<double,14> &evln) const
  {
    LHA2EVLN(basis, evln);
  }

  //_________________________________________________________________________
  void Distribution::EVLN2BASIS(array<double,14> const& evln, array<double,14> &basis) const
  {
    EVLN2LHA(evln, basis);
  }

  //_________________________________________________________________________
  void Distribution::LHA2EVLN(array<double,14> const& LHA, array<double,14> &EVLN)
  {
    const double uplus = LHA[U] + LHA[UBAR];
    const double uminus = LHA[U] - LHA[UBAR];

    const double dplus = LHA[D] + LHA[DBAR];
    const double dminus = LHA[D] - LHA[DBAR];

    const double cplus = LHA[C] + LHA[CBAR];
    const double cminus = LHA[C] - LHA[CBAR];

    const double splus = LHA[S] + LHA[SBAR];
    const double sminus = LHA[S] - LHA[SBAR];

    const double tplus = LHA[T] + LHA[TBAR];
    const double tminus = LHA[T] - LHA[TBAR];

    const double bplus = LHA[B] + LHA[BBAR];
    const double bminus = LHA[B] - LHA[BBAR];

    EVLN[0]= LHA[PHT]; // photon
    EVLN[1]=(uplus + dplus + cplus + splus + tplus + bplus); //Singlet
    EVLN[2]=(LHA[GLUON]); // Gluon

    EVLN[3]=( uminus + dminus + sminus + cminus + bminus + tminus ); //V
    EVLN[4]=( uminus - dminus ); // V3
    EVLN[5]=( uminus + dminus - 2*sminus); // V8
    EVLN[6]=( uminus + dminus + sminus - 3*cminus); //V15
    EVLN[7]=( uminus + dminus + sminus + cminus - 4*bminus ); //V24
    EVLN[8]=( uminus + dminus + sminus + cminus + bminus - 5*tminus); // V35

    EVLN[9]=(  uplus - dplus ); // T3
    EVLN[10]=( uplus + dplus - 2*splus ); // T8
    EVLN[11]=( uplus + dplus + splus - 3*cplus ); //T15
    EVLN[12]=( uplus + dplus + splus + cplus - 4*bplus ); //T24
    EVLN[13]=( uplus + dplus + splus + cplus + bplus - 5*tplus ); // T35

  }

  //_________________________________________________________________________
  void Distribution::EVLN2LHA(array<double,14> const& EVLN, array<double, 14> &LHA)
  {
    // Basis {"PHT","SNG","GLU","VAL","V03","V08","V15","V24","V35","T03","T08","T15","T24","T35"};
    LHA[PHT] = EVLN[0];

    LHA[GLUON] = EVLN[2];

    LHA[U] = ( 10*EVLN[1]
               + 30*EVLN[9] + 10*EVLN[10] + 5*EVLN[11] + 3*EVLN[12] + 2*EVLN[13]
               + 10*EVLN[3] + 30*EVLN[4] + 10*EVLN[5] + 5*EVLN[6] + 3*EVLN[7] + 2*EVLN[8] )
      / 120;

    LHA[UBAR] = ( 10*EVLN[1]
                  + 30*EVLN[9] + 10*EVLN[10] + 5*EVLN[11] + 3*EVLN[12] + 2*EVLN[13]
                  - 10*EVLN[3] - 30*EVLN[4] - 10*EVLN[5] - 5*EVLN[6] - 3*EVLN[7] - 2*EVLN[8] )
      / 120;


    LHA[D] = ( 10*EVLN[1]
               - 30*EVLN[9] + 10*EVLN[10] + 5*EVLN[11] + 3*EVLN[12] + 2*EVLN[13]
               + 10*EVLN[3] - 30*EVLN[4] + 10*EVLN[5] + 5*EVLN[6] + 3*EVLN[7] + 2*EVLN[8] )
      / 120;

    LHA[DBAR] = ( 10*EVLN[1]
                  - 30*EVLN[9] + 10*EVLN[10] + 5*EVLN[11] + 3*EVLN[12] + 2*EVLN[13]
                  - 10*EVLN[3] + 30*EVLN[4] - 10*EVLN[5] - 5*EVLN[6] - 3*EVLN[7] - 2*EVLN[8] )
      / 120;

    LHA[S] = ( 10*EVLN[1]
               - 20*EVLN[10] + 5*EVLN[11] + 3*EVLN[12] + 2*EVLN[13]
               + 10*EVLN[3] - 20*EVLN[5] + 5*EVLN[6] + 3*EVLN[7] + 2*EVLN[8] )
      / 120;

    LHA[SBAR] = ( 10*EVLN[1]
                  - 20*EVLN[10] + 5*EVLN[11] + 3*EVLN[12] + 2*EVLN[13]
                  - 10*EVLN[3] + 20*EVLN[5] - 5*EVLN[6] - 3*EVLN[7] - 2*EVLN[8] )
      / 120;

    LHA[C] = ( 10*EVLN[1]
               - 15*EVLN[11] + 3*EVLN[12] + 2*EVLN[13]
               + 10*EVLN[3] - 15*EVLN[6] + 3*EVLN[7] + 2*EVLN[8] )
      / 120;

    LHA[CBAR] = ( 10*EVLN[1]
                  - 15*EVLN[11] + 3*EVLN[12] + 2*EVLN[13]
                  - 10*EVLN[3] + 15*EVLN[6] - 3*EVLN[7] - 2*EVLN[8] )
      / 120;

    LHA[B] = ( 5*EVLN[1]
               - 6*EVLN[12] + EVLN[13]
               + 5*EVLN[3] - 6*EVLN[7] + EVLN[8] )
      / 60;

    LHA[BBAR] = ( 5*EVLN[1]
                  - 6*EVLN[12] + EVLN[13]
                  - 5*EVLN[3] + 6*EVLN[7] - EVLN[8] )
      / 60;

    LHA[T] = ( EVLN[1]
               - EVLN[13]
               + EVLN[3] - EVLN[8] )
      / 12;

    LHA[TBAR] = ( EVLN[1]
                  - EVLN[13]
                  - EVLN[3] + EVLN[8] )
      / 12;
  }
}
