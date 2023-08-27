//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/evolutionbasisqcdqed.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________
  EvolutionBasisQCDQED::EvolutionBasisQCDQED(int const& nu, int const& nd, int const& nl):
    ConvolutionMap{"EvolutionBasisQCDQED_nu" + std::to_string(nu) + "_nd" + std::to_string(nd) + "_nl" + std::to_string(nl)}
  {
    // Define relevant constants
    const double nf      = nu + nd;
    const double dnf     = nu - nd;
    const double dnfrel  = dnf / nf;
    const double eSigma2 = NC * ( nu * eu2 + nd * ed2 );
    const double dSigma2 = NC * ( nu * eu2 - nd * ed2 );;
    const double etap    = ( eu2 + ed2 ) / 2;
    const double etam    = ( eu2 - ed2 ) / 2;

    // Numbering
    // 0      1      2      3       4       5        6         7         8    9    10   11   12   13   14   15   16   17   18   19
    // GLUON, GAMMA, SIGMA, DSIGMA, SIGMAL, VALENCE, DVALENCE, VALENCEL, T1U, V1U, T2U, V2U, T1D, V1D, T2D, V2D, T1L, V1L, T2L, V2L

    // Singlet
    _rules[GLUON]  = { {PGG, GLUON, 1}, {PGQ, SIGMA, 1},
      {PGGQED, GLUON, eSigma2}, {PGGMQED, GAMMA, eSigma2}, {PGQQED, SIGMA, etap}, {PGQQED, DSIGMA, etap}
    };
    _rules[GAMMA]  = { {PGMGQED, GLUON, eSigma2}, {PGMGMQED, GAMMA, eSigma2}, {PGMQQED, SIGMA, etap}, {PGMQQED, DSIGMA, etap},
      {PGMLQED, SIGMAL, 1}
    };
    _rules[SIGMA]  = { {PQG, GLUON, 1}, {PQQ, SIGMA, 1},
      {PQGQED, GLUON, 2. * eSigma2}, {PQGMQED, GAMMA, 2. * eSigma2},
      {PQQQED, SIGMA,  etap * eSigma2 / nf}, {PNSPQED, SIGMA,  etap - etap * eSigma2 / nf},
      {PQQQED, DSIGMA, etam * eSigma2 / nf}, {PNSPQED, DSIGMA, etam - etam * eSigma2 / nf},
      {PQLQED, SIGMAL, 2. * eSigma2}
    };
    _rules[DSIGMA] = { {PQG, GLUON, dnfrel}, {PQQ, SIGMA, dnfrel}, {PNSP, SIGMA, - dnfrel}, {PNSP, DSIGMA, 1},
      {PQGQED, GLUON, 2. * dSigma2}, {PQGMQED, GAMMA, 2. * dSigma2},
      {PQQQED, SIGMA,  etap * dSigma2 / nf}, {PNSPQED, SIGMA,  etam - etap * dSigma2 / nf},
      {PQQQED, DSIGMA, etam * dSigma2 / nf}, {PNSPQED, DSIGMA, etap - etam * dSigma2 / nf},
      {PQLQED, SIGMAL, 2. * dSigma2}
    };
    _rules[SIGMAL] = { {PLGMQED, GAMMA, 2. * nl}, {PQLQED, SIGMA,  2. * nl * etap}, {PQLQED, DSIGMA, 2. * nl * etam}, {PLLQED, SIGMAL, 1} };

    // Coupled total valences
    _rules[VALENCE]  = { {PNSV, VALENCE, 1},
      {PNSMQED, VALENCE, etap}, {PNSMQED, DVALENCE, etam}
    };
    _rules[DVALENCE] = { {PNSV, VALENCE, dnfrel}, {PNSM, VALENCE, - dnfrel}, {PNSM, DVALENCE, 1},
      {PNSMQED, VALENCE, etam}, {PNSMQED, DVALENCE, etap}
    };

    // Lepton total valence
    _rules[VALENCEL] = { {PNSMQED, VALENCEL, 1} };

    // Non-singlet distributions
    for (int i = 0; i < 2; i++)
      {
        // Up-type
        if (nu > i + 1)
          {
            _rules[8 + 2 * i]     = { {PNSP, 8 + 2 * i,     1}, {PNSPQED, 8 + 2 * i,     eu2} };
            _rules[8 + 2 * i + 1] = { {PNSM, 8 + 2 * i + 1, 1}, {PNSMQED, 8 + 2 * i + 1, eu2} };
          }
        else
          {
            _rules[8 + 2 * i]     = ( _rules[SIGMA]   + _rules[DSIGMA]   ) / 2.;
            _rules[8 + 2 * i + 1] = ( _rules[VALENCE] + _rules[DVALENCE] ) / 2.;
          }
        // Down-type
        if (nd > i + 1)
          {
            _rules[12 + 2 * i]     = { {PNSP, 12 + 2 * i,     1}, {PNSPQED, 12 + 2 * i,     eu2} };
            _rules[12 + 2 * i + 1] = { {PNSM, 12 + 2 * i + 1, 1}, {PNSMQED, 12 + 2 * i + 1, eu2} };
          }
        else
          {
            _rules[12 + 2 * i]     = ( _rules[SIGMA]   - _rules[DSIGMA]   ) / 2.;
            _rules[12 + 2 * i + 1] = ( _rules[VALENCE] - _rules[DVALENCE] ) / 2.;
          }
        // Leptons
        if (nl > i + 1)
          {
            _rules[16 + 2 * i]     = { {PNSPQED, 16 + 2 * i,     1} };
            _rules[16 + 2 * i + 1] = { {PNSMQED, 16 + 2 * i + 1, 1} };
          }
        else
          {
            _rules[16 + 2 * i]     = _rules[SIGMAL];
            _rules[16 + 2 * i + 1] = _rules[VALENCEL];
          }
      }
  }
}
