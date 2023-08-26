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

    // Singlet
    _rules[GLUON]  = { {PGG, GLUON, 1}, {PGQ, SIGMA, 1},
      {PGGQED, GLUON, eSigma2}, {PGGMQED, GAMMA, eSigma2}, {PGQQED, SIGMA, etap}, {PGQQED, DSIGMA, etap}
    };
    _rules[GAMMA]  = { {PGMGQED, GLUON, eSigma2}, {PGMGMQED, GAMMA, eSigma2}, {PGMQQED, SIGMA, etap}, {PGMQQED, DSIGMA, etap},
      {PGMLQED, SIGMAL, 1}
    };
    _rules[SIGMA]  = { {PQG, GLUON, 1}, {PQQ, SIGMA, 1},
      {PQGQED, GLUON, 2. * eSigma2}, {PQGMQED, GAMMA, 2. * eSigma2},
      {PQQQED, SIGMA,  etap * eSigma2 / nf}, {PNSPQED, SIGMA,  etap * ( 1 - eSigma2 / nf )},
      {PQQQED, DSIGMA, etam * eSigma2 / nf}, {PNSPQED, DSIGMA, etam * ( 1 - eSigma2 / nf )},
      {PQLQED, SIGMAL, 2. * eSigma2}
    };
    _rules[DSIGMA] = { {PQG, GLUON, dnfrel}, {PQQ, SIGMA, dnfrel}, {PNSP, SIGMA, - dnfrel}, {PQQ, DSIGMA, 1},
      {PQGQED, GLUON, 2. * dSigma2}, {PQGMQED, GAMMA, 2. * dSigma2},
      {PQQQED, SIGMA,  etap * dSigma2 / nf}, {PNSPQED, SIGMA,  etam - etap * dSigma2 / nf},
      {PQQQED, DSIGMA, etam * dSigma2 / nf}, {PNSPQED, DSIGMA, etap - etam * dSigma2 / nf},
      {PQLQED, SIGMAL, 2. * dSigma2}
    };
    _rules[SIGMAL] = { {PLGMQED, GAMMA, 2. * nl}, {PQLQED, SIGMA,  2. * nl * etap}, {PQLQED, DSIGMA, 2. * nl * etam}, {PLLQED, SIGMAL, 1} };

    // Coupled non-singlet
    _rules[VALENCE]  = { {PNSV, VALENCE, 1},
      {PNSMQED, VALENCE, etap}, {PNSMQED, DVALENCE, etam}
    };
    _rules[DVALENCE] = { {PNSV, VALENCE, dnfrel}, {PNSM, VALENCE, - dnfrel}, {PNSM, DVALENCE, 1},
      {PNSMQED, VALENCE, etam}, {PNSMQED, DVALENCE, etap}
    };

    // Lepton non-singlet
    _rules[VALENCEL] =  { {PNSMQED, VALENCEL, 1} };

    /*
        for (int j = nf + 1; j <= 6; j++)
          _rules[GLUON].push_back({PGQ, 2 * j - 1, 6. / j / ( j - 1 )});

        _rules[SIGMA]   = { {PQG, GLUON, 1}, {PQQ, SIGMA, 1} };
        for (int j = nf + 1; j <= 6; j++)
          _rules[SIGMA].push_back({PQQ, 2 * j - 1, 6. / j / ( j - 1 )});

        _rules[VALENCE] = { {PNSV, VALENCE, 1} };

        for (int i = 2; i <= nf; i++)
          {
            _rules[2 * i - 1] = { {PNSP, 2 * i - 1, 1} };
            _rules[2 * i]     = { {PNSM, 2 * i,     1} };
          }

        for (int i = nf + 1; i <= 6; i++)
          {
            _rules[2 * i - 1] = { {PQG, GLUON, 1}, {PQQ, SIGMA, 1} };
            for (int j = nf + 1; j <= 6; j++)
              _rules[2 * i - 1].push_back({PQQ, 2 * j - 1, 6. / j / ( j - 1 )});
            _rules[2 * i] = { {PNSV, 2 * i, 1} };
          }
    */
  }
}
