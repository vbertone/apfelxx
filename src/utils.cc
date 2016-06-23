/*
  utils.cc:

  Author: Valerio Bertone
*/

using namespace std;

namespace apfel {

  // ================================================================================
  // Rotate physical basis into evolution basis
  enum { TBAR, BBAR, CBAR, SBAR, UBAR, DBAR, GLUON, D, U, S, C, B, T, PHT };
  void Phys2Evol(double const *Phys, double *Evol)
  {
    const double uplus  = Phys[U] + Phys[UBAR];
    const double uminus = Phys[U] - Phys[UBAR];
    
    const double dplus  = Phys[D] + Phys[DBAR];
    const double dminus = Phys[D] - Phys[DBAR];
    
    const double cplus  = Phys[C] + Phys[CBAR];
    const double cminus = Phys[C] - Phys[CBAR];
    
    const double splus  = Phys[S] + Phys[SBAR];
    const double sminus = Phys[S] - Phys[SBAR];
    
    const double tplus  = Phys[T] + Phys[TBAR];
    const double tminus = Phys[T] - Phys[TBAR];
    
    const double bplus  = Phys[B] + Phys[BBAR];
    const double bminus = Phys[B] - Phys[BBAR];

    Evol[0]  = Phys[PHT];                                                   // Photon

    Evol[1]  = ( uplus + dplus + cplus + splus + tplus + bplus );           // Singlet
    Evol[2]  = Phys[GLUON];                                                 // Gluon
        
    Evol[3]  = ( uplus - dplus );                                           // T3
    Evol[4]  = ( uplus + dplus - 2 * splus );                               // T8
    Evol[5]  = ( uplus + dplus + splus - 3 * cplus );                       // T15
    Evol[6]  = ( uplus + dplus + splus + cplus - 4 * bplus );               // T24
    Evol[7]  = ( uplus + dplus + splus + cplus + bplus - 5 * tplus );       // T35

    Evol[8]  = ( uminus + dminus + sminus + cminus + bminus + tminus );     // Total valence

    Evol[9]  = ( uminus - dminus );                                         // V3
    Evol[10] = ( uminus + dminus - 2 * sminus);                             // V8
    Evol[11] = ( uminus + dminus + sminus - 3 * cminus);                    // V15
    Evol[12] = ( uminus + dminus + sminus + cminus - 4 * bminus );          // V24
    Evol[13] = ( uminus + dminus + sminus + cminus + bminus - 5 * tminus ); // V35
  }

}
