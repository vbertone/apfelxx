//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/tmdbuilder.h"
#include "apfel/timer.h"
#include "apfel/matchingfunctionspdf.h"
#include "apfel/matchingfunctionsff.h"
#include "apfel/evolutionbasisqcd.h"
#include "apfel/betaqcd.h"
#include "apfel/gammak.h"
#include "apfel/gammaf.h"
#include "apfel/kcs.h"
#include "apfel/djet.h"
#include "apfel/hardfactors.h"
#include "apfel/tools.h"
#include "apfel/constants.h"
#include "apfel/integrator.h"

namespace apfel
{
  //_____________________________________________________________________________
  std::map<int, TmdObjects> InitializeTmdObjects(Grid                const& g,
                                                 std::vector<double> const& Thresholds,
                                                 double              const& IntEps)
  {
    // Initialise space-like and time-like splitting functions on the
    // grid required to compute the log terms of the matching
    // functions.
    const std::map<int, DglapObjects> DglapObjpdf = InitializeDglapObjectsQCD(g, Thresholds, false, IntEps);
    const std::map<int, DglapObjects> DglapObjff  = InitializeDglapObjectsQCDT(g, Thresholds, false, IntEps);

    report("Initializing TMD objects for matching and evolution... ");
    Timer t;

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the threshold
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
        nfi++;

    // ===============================================================
    // LO matching functions operators
    std::map<int, std::map<int, Operator>> C00;
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, Id});
        OM.insert({EvolutionBasisQCD::PNSM, Id});
        OM.insert({EvolutionBasisQCD::PNSV, Id});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * Id});
        OM.insert({EvolutionBasisQCD::PQG,  Zero});
        OM.insert({EvolutionBasisQCD::PGQ,  Zero});
        OM.insert({EvolutionBasisQCD::PGG,  Id});
        C00.insert({nf, OM});
      }

    // ===============================================================
    // NLO matching functions operators
    // PDFs
    std::map<int, std::map<int, Operator>> C10pdf;
    const Operator O1nspdf{g, C1nspdf{}, IntEps};
    const Operator O1qgpdf{g, C1qgpdf{}, IntEps};
    const Operator O1gqpdf{g, C1gqpdf{}, IntEps};
    const Operator O1ggpdf{g, C1ggpdf{}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O1nspdf});
        OM.insert({EvolutionBasisQCD::PNSM, O1nspdf});
        OM.insert({EvolutionBasisQCD::PNSV, O1nspdf});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O1nspdf});
        OM.insert({EvolutionBasisQCD::PQG,           nf * O1qgpdf});
        OM.insert({EvolutionBasisQCD::PGQ,  ( nf / 6. ) * O1gqpdf});
        OM.insert({EvolutionBasisQCD::PGG,                O1ggpdf});
        C10pdf.insert({nf, OM});
      }

    // Terms proportion to one power of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C11pdf;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O11gmVq = gammaFq0() * Id;
        const Operator O11gmVg = gammaFg0(nf) * Id;
        const auto P0 = DglapObjpdf.at(nf).SplittingFunctions.at(0);
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O11gmVq - 2 * P0.at(0)});
        OM.insert({EvolutionBasisQCD::PNSM, O11gmVq - 2 * P0.at(1)});
        OM.insert({EvolutionBasisQCD::PNSV, O11gmVq - 2 * P0.at(2)});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * ( O11gmVq - 2 * P0.at(3) )});
        OM.insert({EvolutionBasisQCD::PQG,                          - 2 * P0.at(4)});
        OM.insert({EvolutionBasisQCD::PGQ,                - ( nf / 3. ) * P0.at(5)});
        OM.insert({EvolutionBasisQCD::PGG,                  O11gmVg - 2 * P0.at(6)});
        C11pdf.insert({nf, OM});
      }

    // Terms proportion to two powers of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C12;
    const Operator O12gmKq = - CF * gammaK0() / 2 * Id;
    const Operator O12gmKg = - CA * gammaK0() / 2 * Id;
    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O12gmKq});
        OM.insert({EvolutionBasisQCD::PNSM, O12gmKq});
        OM.insert({EvolutionBasisQCD::PNSV, O12gmKq});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O12gmKq});
        OM.insert({EvolutionBasisQCD::PQG,                Zero});
        OM.insert({EvolutionBasisQCD::PGQ,                Zero});
        OM.insert({EvolutionBasisQCD::PGG,                O12gmKg});
        C12.insert({nf, OM});
      }

    // FFs
    std::map<int, std::map<int, Operator>> C10ff;
    const Operator O1nsff{g, C1nsff{}, IntEps};
    const Operator O1qgff{g, C1qgff{}, IntEps};
    const Operator O1gqff{g, C1gqff{}, IntEps};
    const Operator O1ggff{g, C1ggff{}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O1nsff});
        OM.insert({EvolutionBasisQCD::PNSM, O1nsff});
        OM.insert({EvolutionBasisQCD::PNSV, O1nsff});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O1nsff});
        OM.insert({EvolutionBasisQCD::PQG,           nf * O1qgff});
        OM.insert({EvolutionBasisQCD::PGQ,  ( nf / 6. ) * O1gqff});
        OM.insert({EvolutionBasisQCD::PGG,                O1ggff});
        C10ff.insert({nf, OM});
      }

    // Terms proportion to one power of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C11ff;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O11gmVq = gammaFq0() * Id;
        const Operator O11gmVg = gammaFg0(nf) * Id;
        const auto P0 = DglapObjff.at(nf).SplittingFunctions.at(0);
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O11gmVq - 2 * P0.at(0)});
        OM.insert({EvolutionBasisQCD::PNSM, O11gmVq - 2 * P0.at(1)});
        OM.insert({EvolutionBasisQCD::PNSV, O11gmVq - 2 * P0.at(2)});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * ( O11gmVq - 2 * P0.at(3) )});
        OM.insert({EvolutionBasisQCD::PQG,                          - 2 * P0.at(4)});
        OM.insert({EvolutionBasisQCD::PGQ,                - ( nf / 3. ) * P0.at(5)});
        OM.insert({EvolutionBasisQCD::PGG,                  O11gmVg - 2 * P0.at(6)});
        C11ff.insert({nf, OM});
      }

    // Terms proportion to two powers of log(mu0/mub) equal to that of
    // PDFs.

    // ===============================================================
    // NNLO matching functions operators
    // PDFs
    std::map<int, std::map<int, Operator>> C20pdf;
    const Operator O2pspdf{g, C2pspdf{}, IntEps};
    const Operator O2qgpdf{g, C2qgpdf{}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O2gqpdf {g, C2gqpdf{nf},  IntEps};
        const Operator O2ggpdf {g, C2ggpdf{nf},  IntEps};
        const Operator O2nsppdf{g, C2nsppdf{nf}, IntEps};
        const Operator O2nsmpdf{g, C2nsmpdf{nf}, IntEps};
        const Operator O2qqpdf = O2nsppdf + nf * O2pspdf;
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O2nsppdf});
        OM.insert({EvolutionBasisQCD::PNSM, O2nsmpdf});
        OM.insert({EvolutionBasisQCD::PNSV, O2nsmpdf});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O2qqpdf});
        OM.insert({EvolutionBasisQCD::PQG,           nf * O2qgpdf});
        OM.insert({EvolutionBasisQCD::PGQ,  ( nf / 6. ) * O2gqpdf});
        OM.insert({EvolutionBasisQCD::PGG,                O2ggpdf});
        C20pdf.insert({nf, OM});
      }

    // Terms proportion to one power of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C21pdf;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const double b0   = beta0qcd(nf);
        const double gFq0 = gammaFq0();
        const double gFg0 = gammaFg0(nf);
        const Operator O21gmVq = ( gammaFq1(nf) + CF * KCS10(nf) ) * Id;
        const Operator O21gmVg = ( gammaFg1(nf) + CA * KCS10(nf) ) * Id;
        const auto P0 = DglapObjpdf.at(nf).SplittingFunctions.at(0);
        const auto P1 = DglapObjpdf.at(nf).SplittingFunctions.at(1);
        const auto C1 = C10pdf.at(nf);
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O21gmVq - 2 * P1.at(0) + ( gFq0 + 2 * b0 ) * C1.at(0) - 2 * C1.at(0) * P0.at(0)});
        OM.insert({EvolutionBasisQCD::PNSM, O21gmVq - 2 * P1.at(1) + ( gFq0 + 2 * b0 ) * C1.at(1) - 2 * C1.at(1) * P0.at(1)});
        OM.insert({EvolutionBasisQCD::PNSV, O21gmVq - 2 * P1.at(2) + ( gFq0 + 2 * b0 ) * C1.at(2) - 2 * C1.at(2) * P0.at(2)});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * ( O21gmVq - 2 * P1.at(3) + ( gFq0 + 2 * b0 ) * C1.at(3) - 2 * ( C1.at(3) * P0.at(3) + C1.at(4) * P0.at(5) ) )});
        OM.insert({EvolutionBasisQCD::PQG,                          - 2 * P1.at(4) + ( gFq0 + 2 * b0 ) * C1.at(4) - 2 * ( C1.at(3) * P0.at(4) + C1.at(4) * P0.at(6) )});
        OM.insert({EvolutionBasisQCD::PGQ,          ( nf / 6. ) * ( - 2 * P1.at(5) + ( gFg0 + 2 * b0 ) * C1.at(5) - 2 * ( C1.at(5) * P0.at(3) + C1.at(6) * P0.at(5) ) )});
        OM.insert({EvolutionBasisQCD::PGG,                  O21gmVg - 2 * P1.at(6) + ( gFg0 + 2 * b0 ) * C1.at(6) - 2 * ( C1.at(5) * P0.at(4) + C1.at(6) * P0.at(6) )});
        C21pdf.insert({nf, OM});
      }

    // Terms proportion to two powers of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C22pdf;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const double b0   = beta0qcd(nf);
        const double gK0  = gammaK0();
        const double gK1  = gammaK1(nf);
        const double gFq0 = gammaFq0();
        const double gFg0 = gammaFg0(nf);
        const Operator O22gmVq = ( b0 * gFq0 + pow(gFq0, 2) / 2 - CF * gK1 / 2 ) * Id;
        const Operator O22gmVg = ( b0 * gFg0 + pow(gFg0, 2) / 2 - CA * gK1 / 2 ) * Id;
        const auto P0 = DglapObjpdf.at(nf).SplittingFunctions.at(0);
        const auto C1 = C10pdf.at(nf);
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O22gmVq - 2 * ( b0 + gFq0 ) * P0.at(0) - CF * gK0 / 2 * C1.at(0) + 2 * P0.at(0) * P0.at(0)});
        OM.insert({EvolutionBasisQCD::PNSM, O22gmVq - 2 * ( b0 + gFq0 ) * P0.at(1) - CF * gK0 / 2 * C1.at(1) + 2 * P0.at(1) * P0.at(1)});
        OM.insert({EvolutionBasisQCD::PNSV, O22gmVq - 2 * ( b0 + gFq0 ) * P0.at(2) - CF * gK0 / 2 * C1.at(2) + 2 * P0.at(2) * P0.at(2)});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * (O22gmVq - 2 * ( b0 + gFq0 ) * P0.at(3) - CF * gK0 / 2 * C1.at(3) + 2 * ( P0.at(3) * P0.at(3) + P0.at(4) * P0.at(5) ) )});
        OM.insert({EvolutionBasisQCD::PQG,                         - 2 * ( b0 + gFq0 ) * P0.at(4) - CF * gK0 / 2 * C1.at(4) + 2 * ( P0.at(3) * P0.at(4) + P0.at(4) * P0.at(6) )});
        OM.insert({EvolutionBasisQCD::PGQ,          ( nf / 6. ) * (- 2 * ( b0 + gFg0 ) * P0.at(5) - CA * gK0 / 2 * C1.at(5) + 2 * ( P0.at(5) * P0.at(3) + P0.at(6) * P0.at(5) ) )});
        OM.insert({EvolutionBasisQCD::PGG,                 O22gmVg - 2 * ( b0 + gFg0 ) * P0.at(6) - CA * gK0 / 2 * C1.at(6) + 2 * ( P0.at(5) * P0.at(4) + P0.at(6) * P0.at(6) )});
        C22pdf.insert({nf, OM});
      }

    // Terms proportion to three powers of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C23pdf;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const double b0   = beta0qcd(nf);
        const double gK0  = gammaK0();
        const double gFq0 = gammaFq0();
        const double gFg0 = gammaFg0(nf);
        const Operator O23gmVq = CF * gK0 * ( - gFq0 / 2 - 2 * b0 / 3 ) * Id;
        const Operator O23gmVg = CA * gK0 * ( - gFg0 / 2 - 2 * b0 / 3 ) * Id;
        const auto P0 = DglapObjpdf.at(nf).SplittingFunctions.at(0);
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O23gmVq + CF * gK0 * P0.at(0)});
        OM.insert({EvolutionBasisQCD::PNSM, O23gmVq + CF * gK0 * P0.at(1)});
        OM.insert({EvolutionBasisQCD::PNSV, O23gmVq + CF * gK0 * P0.at(2)});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * ( O23gmVq + CF * gK0 * P0.at(3) )});
        OM.insert({EvolutionBasisQCD::PQG,                          + CF * gK0 * P0.at(4)});
        OM.insert({EvolutionBasisQCD::PGQ,           + ( nf / 6. ) *  CA * gK0 * P0.at(5)});
        OM.insert({EvolutionBasisQCD::PGG,                  O23gmVg + CA * gK0 * P0.at(6)});
        C23pdf.insert({nf, OM});
      }

    // Terms proportion to four powers of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C24;
    const double gK0 = gammaK0();
    const Operator O24gmVq = pow(CF * gK0, 2) / 8 * Id;
    const Operator O24gmVg = pow(CA * gK0, 2) / 8 * Id;
    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O24gmVq});
        OM.insert({EvolutionBasisQCD::PNSM, O24gmVq});
        OM.insert({EvolutionBasisQCD::PNSV, O24gmVq});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O24gmVq});
        OM.insert({EvolutionBasisQCD::PQG,                Zero});
        OM.insert({EvolutionBasisQCD::PGQ,                Zero});
        OM.insert({EvolutionBasisQCD::PGG,                O24gmVg});
        C24.insert({nf, OM});
      }

    // FFs
    std::map<int, std::map<int, Operator>> C20ff;
    const Operator O2psff{g, C2psff{}, IntEps};
    const Operator O2qgff{g, C2qgff{}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O2gqff {g, C2gqff{nf},  IntEps};
        const Operator O2ggff {g, C2ggff{nf},  IntEps};
        const Operator O2nspff{g, C2nspff{nf}, IntEps};
        const Operator O2nsmff{g, C2nsmff{nf}, IntEps};
        const Operator O2qqff = O2nspff + nf * O2psff;
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O2nspff});
        OM.insert({EvolutionBasisQCD::PNSM, O2nsmff});
        OM.insert({EvolutionBasisQCD::PNSV, O2nsmff});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O2qqff});
        OM.insert({EvolutionBasisQCD::PQG,           nf * O2qgff});
        OM.insert({EvolutionBasisQCD::PGQ,  ( nf / 6. ) * O2gqff});
        OM.insert({EvolutionBasisQCD::PGG,                O2ggff});
        C20ff.insert({nf, OM});
      }

    // Terms proportion to one power of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C21ff;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const double b0   = beta0qcd(nf);
        const double gFq0 = gammaFq0();
        const double gFg0 = gammaFg0(nf);
        const Operator O21gmVq = ( gammaFq1(nf) + CF * KCS10(nf) ) * Id;
        const Operator O21gmVg = ( gammaFg1(nf) + CA * KCS10(nf) ) * Id;
        const auto P0 = DglapObjff.at(nf).SplittingFunctions.at(0);
        const auto P1 = DglapObjff.at(nf).SplittingFunctions.at(1);
        const auto C1 = C10ff.at(nf);
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O21gmVq - 2 * P1.at(0) + ( gFq0 + 2 * b0 ) * C1.at(0) - 2 * C1.at(0) * P0.at(0)});
        OM.insert({EvolutionBasisQCD::PNSM, O21gmVq - 2 * P1.at(1) + ( gFq0 + 2 * b0 ) * C1.at(1) - 2 * C1.at(1) * P0.at(1)});
        OM.insert({EvolutionBasisQCD::PNSV, O21gmVq - 2 * P1.at(2) + ( gFq0 + 2 * b0 ) * C1.at(2) - 2 * C1.at(2) * P0.at(2)});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * ( O21gmVq - 2 * P1.at(3) + ( gFq0 + 2 * b0 ) * C1.at(3) - 2 * ( C1.at(3) * P0.at(3) + C1.at(4) * P0.at(5) ) )});
        OM.insert({EvolutionBasisQCD::PQG,                          - 2 * P1.at(4) + ( gFq0 + 2 * b0 ) * C1.at(4) - 2 * ( C1.at(3) * P0.at(4) + C1.at(4) * P0.at(6) )});
        OM.insert({EvolutionBasisQCD::PGQ,          ( nf / 6. ) * ( - 2 * P1.at(5) + ( gFg0 + 2 * b0 ) * C1.at(5) - 2 * ( C1.at(5) * P0.at(3) + C1.at(6) * P0.at(5) ) )});
        OM.insert({EvolutionBasisQCD::PGG,                  O21gmVg - 2 * P1.at(6) + ( gFg0 + 2 * b0 ) * C1.at(6) - 2 * ( C1.at(5) * P0.at(4) + C1.at(6) * P0.at(6) )});
        C21ff.insert({nf, OM});
      }

    // Terms proportion to two powers of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C22ff;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const double b0    = beta0qcd(nf);
        const double gv0q  = gammaFq0();
        const double gv0g  = gammaFg0(nf);
        const double gcp   = gammaK1(nf);
        const double factq = - ( gv0q + b0 ) / 2;
        const double factg = - ( gv0g + b0 ) / 2;
        const Operator O22gmVq = ( CF * gcp / 8 + gv0q * ( 2 * b0 + gv0q ) / 8 ) * Id;
        const Operator O22gmVg = ( CA * gcp / 8 + gv0g * ( 2 * b0 + gv0g ) / 8 ) * Id;
        const auto P0 = DglapObjff.at(nf).SplittingFunctions.at(0);
        const auto C1 = C10ff.at(nf);
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O22gmVq + factq * P0.at(0) + CF * C1.at(0) + P0.at(0) * P0.at(0) / 2});
        OM.insert({EvolutionBasisQCD::PNSM, O22gmVq + factq * P0.at(1) + CF * C1.at(1) + P0.at(1) * P0.at(1) / 2});
        OM.insert({EvolutionBasisQCD::PNSV, O22gmVq + factq * P0.at(2) + CF * C1.at(2) + P0.at(2) * P0.at(2) / 2});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * ( O22gmVq + factq * P0.at(3) + CF * C1.at(3) + ( P0.at(3) * P0.at(3) + P0.at(4) * P0.at(5) ) / 2 )});
        OM.insert({EvolutionBasisQCD::PQG,                            factq * P0.at(4) + CF * C1.at(4) + ( P0.at(3) * P0.at(4) + P0.at(4) * P0.at(6) ) / 2});
        OM.insert({EvolutionBasisQCD::PGQ,            ( nf / 6. ) * ( factg * P0.at(5) + CA * C1.at(5) + ( P0.at(5) * P0.at(3) + P0.at(6) * P0.at(5) ) / 2 )});
        OM.insert({EvolutionBasisQCD::PGG,                  O22gmVg + factg * P0.at(6) + CA * C1.at(6) + ( P0.at(5) * P0.at(4) + P0.at(6) * P0.at(6) ) / 2});
        C22ff.insert({nf, OM});
      }

    // Terms proportion to three powers of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C23ff;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const double b0   = beta0qcd(nf);
        const double gK0  = gammaK0();
        const double gFq0 = gammaFq0();
        const double gFg0 = gammaFg0(nf);
        const Operator O23gmVq = CF * gK0 * ( - gFq0 / 2 - 2 * b0 / 3 ) * Id;
        const Operator O23gmVg = CA * gK0 * ( - gFg0 / 2 - 2 * b0 / 3 ) * Id;
        const auto P0 = DglapObjff.at(nf).SplittingFunctions.at(0);
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O23gmVq + CF * gK0 * P0.at(0)});
        OM.insert({EvolutionBasisQCD::PNSM, O23gmVq + CF * gK0 * P0.at(1)});
        OM.insert({EvolutionBasisQCD::PNSV, O23gmVq + CF * gK0 * P0.at(2)});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * ( O23gmVq + CF * gK0 * P0.at(3) )});
        OM.insert({EvolutionBasisQCD::PQG,                          + CF * gK0 * P0.at(4)});
        OM.insert({EvolutionBasisQCD::PGQ,            + ( nf / 6. ) * CA * gK0 * P0.at(5)});
        OM.insert({EvolutionBasisQCD::PGG,                  O23gmVg + CA * gK0 * P0.at(6)});
        C23ff.insert({nf, OM});
      }

    // Terms proportion to four powers of log(mu0/mub) equal to that
    // of PDFs.

    // ===============================================================
    // NNNLO matching functions operators
    // PDFs
    std::map<int, std::map<int, Operator>> C30pdf;
    const Operator O3pvpdf{g, C3pvpdf{}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O3pspdf {g, C3pspdf{nf},  IntEps};
        const Operator O3qgpdf {g, C3qgpdf{nf},  IntEps};
        const Operator O3gqpdf {g, C3gqpdf{nf},  IntEps};
        const Operator O3ggpdf {g, C3ggpdf{nf},  IntEps};
        const Operator O3nsppdf{g, C3nsppdf{nf}, IntEps};
        const Operator O3nsmpdf{g, C3nsmpdf{nf}, IntEps};
        const Operator O3qqpdf  = O3nsppdf + nf * O3pspdf;
        const Operator O3nsvpdf = O3nsmpdf + nf * O3pvpdf;
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O3nsppdf});
        OM.insert({EvolutionBasisQCD::PNSM, O3nsmpdf});
        OM.insert({EvolutionBasisQCD::PNSV, O3nsvpdf});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O3qqpdf});
        OM.insert({EvolutionBasisQCD::PQG,           nf * O3qgpdf});
        OM.insert({EvolutionBasisQCD::PGQ,  ( nf / 6. ) * O3gqpdf});
        OM.insert({EvolutionBasisQCD::PGG,                O3ggpdf});
        C30pdf.insert({nf, OM});
      }

    // FFs
    std::map<int, std::map<int, Operator>> C30ff;
    const Operator O3pvff{g, C3pvff{}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O3psff {g, C3psff{nf},  IntEps};
        const Operator O3qgff {g, C3qgff{nf},  IntEps};
        const Operator O3gqff {g, C3gqff{nf},  IntEps};
        const Operator O3ggff {g, C3ggff{nf},  IntEps};
        const Operator O3nspff{g, C3nspff{nf}, IntEps};
        const Operator O3nsmff{g, C3nsmff{nf}, IntEps};
        const Operator O3qqff  = O3nspff + nf * O3psff;
        const Operator O3nsvff = O3nsmff + nf * O3pvff;
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O3nspff});
        OM.insert({EvolutionBasisQCD::PNSM, O3nsmff});
        OM.insert({EvolutionBasisQCD::PNSV, O3nsvff});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O3qqff});
        OM.insert({EvolutionBasisQCD::PQG,           nf * O3qgff});
        OM.insert({EvolutionBasisQCD::PGQ,  ( nf / 6. ) * O3gqff});
        OM.insert({EvolutionBasisQCD::PGG,                O3ggff});
        C30ff.insert({nf, OM});
      }

    // Define map containing the TmdObjects for each nf
    std::map<int, TmdObjects> TmdObj;

    // Construct a set of operators for each perturbative order for
    // the matching functions. Initialize also coefficients of: beta
    // function, gammaK, gammaF, and Collins-Soper anomalous
    // dimensions.
    for (int nf = nfi; nf <= nff; nf++)
      {
        TmdObjects obj;

        // Threshold
        obj.Threshold = Thresholds[nf-1];

        // Beta function
        obj.Beta.insert({0, beta0qcd(nf)});
        obj.Beta.insert({1, beta1qcd(nf)});
        obj.Beta.insert({2, beta2qcd(nf)});

        // GammaF quark
        obj.GammaFq.insert({0, gammaFq0()});
        obj.GammaFq.insert({1, gammaFq1(nf)});
        obj.GammaFq.insert({2, gammaFq2(nf)});

        // GammaF gluon
        obj.GammaFg.insert({0, gammaFg0(nf)});
        obj.GammaFg.insert({1, gammaFg1(nf)});
        obj.GammaFg.insert({2, gammaFg2(nf)});

        // gammaK (multiply by CF for quarks and by CA for gluons)
        obj.GammaK.insert({0, gammaK0()});
        obj.GammaK.insert({1, gammaK1(nf)});
        obj.GammaK.insert({2, gammaK2(nf)});
        obj.GammaK.insert({3, gammaK3(nf)});

        // Collins-Soper anomalous dimensions (multiply by CF for
        // quarks and by CA for gluons).
        obj.KCS.insert({0, {KCS00(),   KCS01()}});
        obj.KCS.insert({1, {KCS10(nf), KCS11(nf), KCS12(nf)}});
        obj.KCS.insert({2, {KCS20(nf), KCS21(nf), KCS22(nf), KCS23(nf)}});

        // Matching functions
        const EvolutionBasisQCD evb{nf};

        // PDFs
        obj.MatchingFunctionsPDFs.insert({0, {{evb, C00.at(nf)}}});
        obj.MatchingFunctionsPDFs.insert({1, {{evb, C10pdf.at(nf)}, {evb, C11pdf.at(nf)}, {evb, C12.at(nf)}}});
        obj.MatchingFunctionsPDFs.insert({2, {{evb, C20pdf.at(nf)}, {evb, C21pdf.at(nf)}, {evb, C22pdf.at(nf)}, {evb, C23pdf.at(nf)}, {evb, C24.at(nf)}}});
        obj.MatchingFunctionsPDFs.insert({3, {{evb, C30pdf.at(nf)}}});

        // FFs
        obj.MatchingFunctionsFFs.insert({0, {{evb, C00.at(nf)}}});
        obj.MatchingFunctionsFFs.insert({1, {{evb, C10ff.at(nf)}, {evb, C11ff.at(nf)}, {evb, C12.at(nf)}}});
        obj.MatchingFunctionsFFs.insert({2, {{evb, C20ff.at(nf)}, {evb, C21ff.at(nf)}, {evb, C22ff.at(nf)}, {evb, C23ff.at(nf)}, {evb, C24.at(nf)}}});
        obj.MatchingFunctionsFFs.insert({3, {{evb, C30ff.at(nf)}}});

        // Hard factors (set to zero when unknown). In addition,
        // H3Ch() should be multiplied by N_{nf,j} =
        // (\sum_{i=1}^{nf}e_q) / e_j channel by channel. Here we set
        // this factor to the constant 1 / 4 because it results from
        // the average ( N_{5,u} + N_{5,d} ) / 2. The numerical
        // difference for different choices is however very small.
        obj.HardFactors.insert({"DY",    {{1, H1DY()},    {2, H2DY(nf)},    {3, H3DY(nf)    + H3Ch() / 4}}});
        obj.HardFactors.insert({"SIDIS", {{1, H1SIDIS()}, {2, H2SIDIS(nf)}, {3, H3SIDIS(nf) + H3Ch() / 4}}});
        obj.HardFactors.insert({"ggH",   {{1, H1ggH()},   {2, H2ggH(nf)},   {3, 0}}});

        // Insert full object
        TmdObj.insert({nf, obj});
      }
    t.stop();

    return TmdObj;
  }

  //_____________________________________________________________________________
  std::map<int, TmdObjects> InitializeTmdObjectsDYResScheme(Grid                const& g,
                                                            std::vector<double> const& Thresholds,
                                                            double              const& IntEps)
  {
    // Get TMD objects
    std::map<int, TmdObjects> TmdObjs = InitializeTmdObjects(g, Thresholds, IntEps);

    // Run of all active flavours
    for (auto& to : TmdObjs)
      {
        // Transform non-cusp anomalous dimension of the quarks
        to.second.GammaFq.at(2) += - to.second.Beta.at(1) * to.second.HardFactors.at("DY").at(1)
                                   - 2 * to.second.Beta.at(0) * ( to.second.HardFactors.at("DY").at(2) - pow(to.second.HardFactors.at("DY").at(1), 2) / 2 );
        to.second.GammaFq.at(1) += - to.second.Beta.at(0) * to.second.HardFactors.at("DY").at(1);

        // Transform PDF matching functions. Only non log terms need
        // to me modified.
        to.second.MatchingFunctionsPDFs.at(2)[0] += ( to.second.HardFactors.at("DY").at(1) / 2 ) * to.second.MatchingFunctionsPDFs.at(1)[0]
                                                    + ( ( to.second.HardFactors.at("DY").at(2) - pow(to.second.HardFactors.at("DY").at(1), 2) / 4 ) / 2 ) * to.second.MatchingFunctionsPDFs.at(0)[0];
        to.second.MatchingFunctionsPDFs.at(1)[0] += ( to.second.HardFactors.at("DY").at(1) / 2 ) * to.second.MatchingFunctionsPDFs.at(0)[0];

        // Transform hard function for DY, set all coefficients to
        // zero.
        for (auto& h : to.second.HardFactors.at("DY"))
          h.second = 0;
      }

    return TmdObjs;
  }

  //_____________________________________________________________________________
  std::map<int, TmdObjects> InitializeTmdObjectsBM(Grid                const& g,
                                                   std::vector<double> const& Thresholds,
                                                   double              const& IntEps)
  {
    // Initialise space-like splitting functions on the grid required
    // to compute the log terms of the matching functions.
    const std::map<int, DglapObjects> DglapObjpdf = InitializeDglapObjectsQCD(g, Thresholds, false, IntEps);

    report("Initializing TMD objects for matching and evolution of the Boer-Mulders gluon TMD... ");
    Timer t;

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the threshold
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
        nfi++;

    // ===============================================================
    // Identity and zero operators
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};

    // ===============================================================
    // NLO matching functions operators
    // PDFs
    std::map<int, std::map<int, Operator>> C10pdf;
    const Operator O1gqpdf{g, C1gqpdfBM{}, IntEps};
    const Operator O1ggpdf{g, C1ggpdfBM{}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        for (int iOp = 0; iOp < 5; iOp++)
          OM.insert({iOp, Zero});
        OM.insert({EvolutionBasisQCD::PGQ,  ( nf / 6. ) * O1gqpdf});
        OM.insert({EvolutionBasisQCD::PGG,                O1ggpdf});
        C10pdf.insert({nf, OM});
      }

    // Terms proportion to one power of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C11pdf;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O11gmVg = gammaFg0(nf) * Id;
        const auto P0 = DglapObjpdf.at(nf).SplittingFunctions.at(0);
        std::map<int, Operator> OM;
        for (int iOp = 0; iOp < 5; iOp++)
          OM.insert({iOp, Zero});
        OM.insert({EvolutionBasisQCD::PGQ,  - ( nf / 3. ) * P0.at(5)});
        OM.insert({EvolutionBasisQCD::PGG,    O11gmVg - 2 * P0.at(6)});
        C11pdf.insert({nf, OM});
      }

    // Terms proportion to two powers of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C12;
    const Operator O12gmKg = - CA * gammaK0() / 2 * Id;
    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        for (int iOp = 0; iOp < 5; iOp++)
          OM.insert({iOp, Zero});
        OM.insert({EvolutionBasisQCD::PGQ, Zero});
        OM.insert({EvolutionBasisQCD::PGG, O12gmKg});
        C12.insert({nf, OM});
      }

    // Terms proportion to two powers of log(mu0/mub) equal to that of
    // PDFs.

    // ===============================================================
    // NNLO matching functions operators
    // PDFs
    std::map<int, std::map<int, Operator>> C20pdf;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O2gqpdf{g, C2gqpdf{nf}, IntEps};
        const Operator O2ggpdf{g, C2ggpdf{nf}, IntEps};
        std::map<int, Operator> OM;
        for (int iOp = 0; iOp < 5; iOp++)
          OM.insert({iOp, Zero});
        OM.insert({EvolutionBasisQCD::PGQ,  ( nf / 6. ) * O2gqpdf});
        OM.insert({EvolutionBasisQCD::PGG,                O2ggpdf});
        C20pdf.insert({nf, OM});
      }

    // Terms proportion to one power of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C21pdf;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const double b0   = beta0qcd(nf);
        const double gFg0 = gammaFg0(nf);
        const Operator O21gmVg = ( gammaFg1(nf) + CA * KCS10(nf) ) * Id;
        const auto P0 = DglapObjpdf.at(nf).SplittingFunctions.at(0);
        const auto P1 = DglapObjpdf.at(nf).SplittingFunctions.at(1);
        const auto C1 = C10pdf.at(nf);
        std::map<int, Operator> OM;
        for (int iOp = 0; iOp < 5; iOp++)
          OM.insert({iOp, Zero});
        OM.insert({EvolutionBasisQCD::PGQ, ( nf / 6. ) * ( - 2 * P1.at(5) + ( gFg0 + 2 * b0 ) * C1.at(5) - 2 * ( C1.at(5) * P0.at(3) + C1.at(6) * P0.at(5) ) )});
        OM.insert({EvolutionBasisQCD::PGG,         O21gmVg - 2 * P1.at(6) + ( gFg0 + 2 * b0 ) * C1.at(6) - 2 * ( C1.at(5) * P0.at(4) + C1.at(6) * P0.at(6) )});
        C21pdf.insert({nf, OM});
      }

    // Terms proportion to two powers of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C22pdf;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const double b0   = beta0qcd(nf);
        const double gK0  = gammaK0();
        const double gK1  = gammaK1(nf);
        const double gFg0 = gammaFg0(nf);
        const Operator O22gmVg = ( b0 * gFg0 + pow(gFg0, 2) / 2 - CA * gK1 / 2 ) * Id;
        const auto P0 = DglapObjpdf.at(nf).SplittingFunctions.at(0);
        const auto C1 = C10pdf.at(nf);
        std::map<int, Operator> OM;
        for (int iOp = 0; iOp < 5; iOp++)
          OM.insert({iOp, Zero});
        OM.insert({EvolutionBasisQCD::PGQ, ( nf / 6. ) * (- 2 * ( b0 + gFg0 ) * P0.at(5) - CA * gK0 / 2 * C1.at(5) + 2 * ( P0.at(5) * P0.at(3) + P0.at(6) * P0.at(5) ) )});
        OM.insert({EvolutionBasisQCD::PGG,        O22gmVg - 2 * ( b0 + gFg0 ) * P0.at(6) - CA * gK0 / 2 * C1.at(6) + 2 * ( P0.at(5) * P0.at(4) + P0.at(6) * P0.at(6) )});
        C22pdf.insert({nf, OM});
      }

    // Terms proportion to three powers of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C23pdf;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const double b0   = beta0qcd(nf);
        const double gK0  = gammaK0();
        const double gFg0 = gammaFg0(nf);
        const Operator O23gmVg = CA * gK0 * ( - gFg0 / 2 - 2 * b0 / 3 ) * Id;
        const auto P0 = DglapObjpdf.at(nf).SplittingFunctions.at(0);
        std::map<int, Operator> OM;
        for (int iOp = 0; iOp < 5; iOp++)
          OM.insert({iOp, Zero});
        OM.insert({EvolutionBasisQCD::PGQ, + ( nf / 6. ) *  CA * gK0 * P0.at(5)});
        OM.insert({EvolutionBasisQCD::PGG,        O23gmVg + CA * gK0 * P0.at(6)});
        C23pdf.insert({nf, OM});
      }

    // Terms proportion to four powers of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C24;
    const double gK0 = gammaK0();
    const Operator O24gmVg = pow(CA * gK0, 2) / 8 * Id;
    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        for (int iOp = 0; iOp < 5; iOp++)
          OM.insert({iOp, Zero});
        OM.insert({EvolutionBasisQCD::PGQ, Zero});
        OM.insert({EvolutionBasisQCD::PGG, O24gmVg});
        C24.insert({nf, OM});
      }

    // FFs (unknown, thus set to zero)
    std::map<int, Operator> ZeroOp;
    for (int iOp = 0; iOp < 7; iOp++)
      ZeroOp.insert({iOp, Zero});

    // Define map containing the TmdObjects for each nf
    std::map<int, TmdObjects> TmdObj;

    // Construct a set of operators for each perturbative order for
    // the matching functions. Initialize also coefficients of: beta
    // function, gammaK, gammaF, and Collins-Soper anomalous
    // dimensions.
    for (int nf = nfi; nf <= nff; nf++)
      {
        TmdObjects obj;

        // Threshold
        obj.Threshold = Thresholds[nf-1];

        // Beta function
        obj.Beta.insert({0, beta0qcd(nf)});
        obj.Beta.insert({1, beta1qcd(nf)});
        obj.Beta.insert({2, beta2qcd(nf)});

        // GammaF quark
        obj.GammaFq.insert({0, 0});
        obj.GammaFq.insert({1, 0});
        obj.GammaFq.insert({2, 0});

        // GammaF gluon
        obj.GammaFg.insert({0, gammaFg0(nf)});
        obj.GammaFg.insert({1, gammaFg1(nf)});
        obj.GammaFg.insert({2, gammaFg2(nf)});

        // gammaK (multiply by CF for quarks and by CA for gluons)
        obj.GammaK.insert({0, gammaK0()});
        obj.GammaK.insert({1, gammaK1(nf)});
        obj.GammaK.insert({2, gammaK2(nf)});
        obj.GammaK.insert({3, gammaK3(nf)});

        // Collins-Soper anomalous dimensions (multiply by CF for
        // quarks and by CA for gluons).
        obj.KCS.insert({0, {KCS00(),   KCS01()}});
        obj.KCS.insert({1, {KCS10(nf), KCS11(nf), KCS12(nf)}});
        obj.KCS.insert({2, {KCS20(nf), KCS21(nf), KCS22(nf), KCS23(nf)}});

        // Matching functions
        const EvolutionBasisQCD evb{nf};

        // PDFs
        obj.MatchingFunctionsPDFs.insert({0, {{evb, ZeroOp}}});
        obj.MatchingFunctionsPDFs.insert({1, {{evb, C10pdf.at(nf)}, {evb, C11pdf.at(nf)}, {evb, C12.at(nf)}}});
        obj.MatchingFunctionsPDFs.insert({2, {{evb, C20pdf.at(nf)}, {evb, C21pdf.at(nf)}, {evb, C22pdf.at(nf)}, {evb, C23pdf.at(nf)}, {evb, C24.at(nf)}}});
        obj.MatchingFunctionsPDFs.insert({3, {{evb, ZeroOp}}});

        // FFs
        obj.MatchingFunctionsFFs.insert({0, {{evb, ZeroOp}}});
        obj.MatchingFunctionsFFs.insert({1, {{evb, ZeroOp}, {evb, ZeroOp}, {evb, ZeroOp}}});
        obj.MatchingFunctionsFFs.insert({2, {{evb, ZeroOp}, {evb, ZeroOp}, {evb, ZeroOp}, {evb, ZeroOp}, {evb, ZeroOp}}});
        obj.MatchingFunctionsFFs.insert({3, {{evb, ZeroOp}}});

        // Hard factors different from zero only for Higgs production.
        obj.HardFactors.insert({"DY",    {{1, 0}, {2, 0}, {3, 0}}});
        obj.HardFactors.insert({"SIDIS", {{1, 0}, {2, 0}, {3, 0}}});
        obj.HardFactors.insert({"ggH",   {{1, H1ggH()}, {2, H2ggH(nf)}, {3, 0}}});

        // Insert full object
        TmdObj.insert({nf, obj});
      }
    t.stop();

    return TmdObj;
  }

  //_____________________________________________________________________________
  std::map<int, TmdObjects> InitializeTmdObjectsSivers(Grid                const& g,
                                                       std::vector<double> const& Thresholds,
                                                       double              const& IntEps)
  {
    // Do not compute the log terms to the matching function because
    // the collinear evolution of the Qiu-Sterman is not exactly
    // known.

    report("Initializing TMD objects for matching and evolution of the Sivers quark TMD... ");
    Timer t;

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the threshold
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
        nfi++;

    // ===============================================================
    // LO matching functions operators
    std::map<int, std::map<int, Operator>> C00;
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, Id});
        OM.insert({EvolutionBasisQCD::PNSM, Id});
        OM.insert({EvolutionBasisQCD::PNSV, Id});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * Id});
        OM.insert({EvolutionBasisQCD::PQG,  Zero});
        OM.insert({EvolutionBasisQCD::PGQ,  Zero});
        OM.insert({EvolutionBasisQCD::PGG,  Zero});
        C00.insert({nf, OM});
      }

    // ===============================================================
    // NLO matching functions operators
    // PDFs
    std::map<int, std::map<int, Operator>> C10pdf;
    const Operator O1nspdf{g, C1nspdfSivers{}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O1nspdf});
        OM.insert({EvolutionBasisQCD::PNSM, O1nspdf});
        OM.insert({EvolutionBasisQCD::PNSV, O1nspdf});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O1nspdf});
        OM.insert({EvolutionBasisQCD::PQG,  Zero});
        OM.insert({EvolutionBasisQCD::PGQ,  Zero});
        OM.insert({EvolutionBasisQCD::PGG,  Zero});
        C10pdf.insert({nf, OM});
      }

    // FFs (not implemented, set to zero)
    std::map<int, Operator> ZeroOp;
    for (int iOp = 0; iOp < 7; iOp++)
      ZeroOp.insert({iOp, Zero});

    // Define map containing the TmdObjects for each nf
    std::map<int, TmdObjects> TmdObj;

    // Construct a set of operators for each perturbative order for
    // the matching functions. Initialize also coefficients of: beta
    // function, gammaK, gammaF, and Collins-Soper anomalous
    // dimensions.
    for (int nf = nfi; nf <= nff; nf++)
      {
        TmdObjects obj;

        // Threshold
        obj.Threshold = Thresholds[nf-1];

        // Beta function
        obj.Beta.insert({0, beta0qcd(nf)});
        obj.Beta.insert({1, beta1qcd(nf)});
        obj.Beta.insert({2, beta2qcd(nf)});

        // GammaF quark
        obj.GammaFq.insert({0, gammaFq0()});
        obj.GammaFq.insert({1, gammaFq1(nf)});
        obj.GammaFq.insert({2, gammaFq2(nf)});

        // GammaF gluon
        obj.GammaFg.insert({0, gammaFg0(nf)});
        obj.GammaFg.insert({1, gammaFg1(nf)});
        obj.GammaFg.insert({2, gammaFg2(nf)});

        // gammaK (multiply by CF for quarks and by CA for gluons)
        obj.GammaK.insert({0, gammaK0()});
        obj.GammaK.insert({1, gammaK1(nf)});
        obj.GammaK.insert({2, gammaK2(nf)});
        obj.GammaK.insert({3, gammaK3(nf)});

        // Collins-Soper anomalous dimensions (multiply by CF for
        // quarks and by CA for gluons).
        obj.KCS.insert({0, {KCS00(),   KCS01()}});
        obj.KCS.insert({1, {KCS10(nf), KCS11(nf), KCS12(nf)}});
        obj.KCS.insert({2, {KCS20(nf), KCS21(nf), KCS22(nf), KCS23(nf)}});

        // Matching functions
        const EvolutionBasisQCD evb{nf};

        // PDFs
        obj.MatchingFunctionsPDFs.insert({0, {{evb, C00.at(nf)}}});
        obj.MatchingFunctionsPDFs.insert({1, {{evb, C10pdf.at(nf)}, {evb, ZeroOp}, {evb, ZeroOp}}});
        obj.MatchingFunctionsPDFs.insert({2, {{evb, ZeroOp}, {evb, ZeroOp}, {evb, ZeroOp}, {evb, ZeroOp}, {evb, ZeroOp}}});
        obj.MatchingFunctionsPDFs.insert({3, {{evb, ZeroOp}}});

        // FFs
        obj.MatchingFunctionsFFs.insert({0, {{evb, ZeroOp}}});
        obj.MatchingFunctionsFFs.insert({1, {{evb, ZeroOp}, {evb, ZeroOp}, {evb, ZeroOp}}});
        obj.MatchingFunctionsFFs.insert({2, {{evb, ZeroOp}, {evb, ZeroOp}, {evb, ZeroOp}, {evb, ZeroOp}, {evb, ZeroOp}}});
        obj.MatchingFunctionsFFs.insert({3, {{evb, ZeroOp}}});

        // Hard factors (set to zero when unknown). In addition,
        // H3Ch() should be multiplied by N_{nf,j} =
        // (\sum_{i=1}^{nf}e_q) / e_j channel by channel. Here we set
        // this factor to the constant 1 / 4 because it results from
        // the average ( N_{5,u} + N_{5,d} ) / 2. The numerical
        // difference for different choices is however very small.
        obj.HardFactors.insert({"DY",    {{1, H1DY()},    {2, H2DY(nf)},    {3, H3DY(nf)    + H3Ch() / 4}}});
        obj.HardFactors.insert({"SIDIS", {{1, H1SIDIS()}, {2, H2SIDIS(nf)}, {3, H3SIDIS(nf) + H3Ch() / 4}}});
        obj.HardFactors.insert({"ggH",   {{1, H1ggH()},   {2, H2ggH(nf)},   {3, 0}}});

        // Insert full object
        TmdObj.insert({nf, obj});
      }
    t.stop();

    return TmdObj;
  }

  //_____________________________________________________________________________
  std::map<int, TmdObjects> InitializeTmdObjectsg1(Grid                const& g,
                                                   std::vector<double> const& Thresholds,
                                                   double              const& IntEps)
  {
    // Initialise space-like longitudinally polarised splitting
    // functions on the grid required to compute the log terms of the
    // matching functions.
    const std::map<int, DglapObjects> DglapObjpdf = InitializeDglapObjectsQCDpol(g, Thresholds, false, IntEps);

    report("Initializing TMD objects for matching and evolution of g1... ");
    Timer t;

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the threshold
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
        nfi++;

    // ===============================================================
    // LO matching functions operators
    std::map<int, std::map<int, Operator>> C00;
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, Id});
        OM.insert({EvolutionBasisQCD::PNSM, Id});
        OM.insert({EvolutionBasisQCD::PNSV, Id});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * Id});
        OM.insert({EvolutionBasisQCD::PQG,  Zero});
        OM.insert({EvolutionBasisQCD::PGQ,  Zero});
        OM.insert({EvolutionBasisQCD::PGG,  Id});
        C00.insert({nf, OM});
      }

    // ===============================================================
    // NLO matching functions operators
    // PDFs
    std::map<int, std::map<int, Operator>> C10pdf;
    const Operator O1nspdf{g, C1nspdfg1{}, IntEps};
    const Operator O1qgpdf{g, C1qgpdfg1{}, IntEps};
    const Operator O1gqpdf{g, C1gqpdfg1{}, IntEps};
    const Operator O1ggpdf{g, C1ggpdfg1{}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O1nspdf});
        OM.insert({EvolutionBasisQCD::PNSM, O1nspdf});
        OM.insert({EvolutionBasisQCD::PNSV, O1nspdf});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O1nspdf});
        OM.insert({EvolutionBasisQCD::PQG,           nf * O1qgpdf});
        OM.insert({EvolutionBasisQCD::PGQ,  ( nf / 6. ) * O1gqpdf});
        OM.insert({EvolutionBasisQCD::PGG,                O1ggpdf});
        C10pdf.insert({nf, OM});
      }

    // Terms proportion to one power of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C11pdf;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O11gmVq = gammaFq0() * Id;
        const Operator O11gmVg = gammaFg0(nf) * Id;
        const auto P0 = DglapObjpdf.at(nf).SplittingFunctions.at(0);
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O11gmVq - 2 * P0.at(0)});
        OM.insert({EvolutionBasisQCD::PNSM, O11gmVq - 2 * P0.at(1)});
        OM.insert({EvolutionBasisQCD::PNSV, O11gmVq - 2 * P0.at(2)});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * ( O11gmVq - 2 * P0.at(3) )});
        OM.insert({EvolutionBasisQCD::PQG,                          - 2 * P0.at(4)});
        OM.insert({EvolutionBasisQCD::PGQ,                - ( nf / 3. ) * P0.at(5)});
        OM.insert({EvolutionBasisQCD::PGG,                  O11gmVg - 2 * P0.at(6)});
        C11pdf.insert({nf, OM});
      }

    // Terms proportion to two powers of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C12;
    const Operator O12gmKq = - CF * gammaK0() / 2 * Id;
    const Operator O12gmKg = - CA * gammaK0() / 2 * Id;
    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O12gmKq});
        OM.insert({EvolutionBasisQCD::PNSM, O12gmKq});
        OM.insert({EvolutionBasisQCD::PNSV, O12gmKq});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O12gmKq});
        OM.insert({EvolutionBasisQCD::PQG,                Zero});
        OM.insert({EvolutionBasisQCD::PGQ,                Zero});
        OM.insert({EvolutionBasisQCD::PGG,                O12gmKg});
        C12.insert({nf, OM});
      }

    // Construct zero set of operators for the unknown bits
    std::map<int, Operator> ZeroOp;
    for (int iOp = 0; iOp < 7; iOp++)
      ZeroOp.insert({iOp, Zero});

    // Define map containing the TmdObjects for each nf
    std::map<int, TmdObjects> TmdObj;

    // Construct a set of operators for each perturbative order for
    // the matching functions. Initialize also coefficients of: beta
    // function, gammaK, gammaF, and Collins-Soper anomalous
    // dimensions.
    for (int nf = nfi; nf <= nff; nf++)
      {
        TmdObjects obj;

        // Threshold
        obj.Threshold = Thresholds[nf-1];

        // Beta function
        obj.Beta.insert({0, beta0qcd(nf)});
        obj.Beta.insert({1, beta1qcd(nf)});
        obj.Beta.insert({2, beta2qcd(nf)});

        // GammaF quark
        obj.GammaFq.insert({0, gammaFq0()});
        obj.GammaFq.insert({1, gammaFq1(nf)});
        obj.GammaFq.insert({2, gammaFq2(nf)});

        // GammaF gluon
        obj.GammaFg.insert({0, gammaFg0(nf)});
        obj.GammaFg.insert({1, gammaFg1(nf)});
        obj.GammaFg.insert({2, gammaFg2(nf)});

        // gammaK (multiply by CF for quarks and by CA for gluons)
        obj.GammaK.insert({0, gammaK0()});
        obj.GammaK.insert({1, gammaK1(nf)});
        obj.GammaK.insert({2, gammaK2(nf)});
        obj.GammaK.insert({3, gammaK3(nf)});

        // Collins-Soper anomalous dimensions (multiply by CF for
        // quarks and by CA for gluons).
        obj.KCS.insert({0, {KCS00(),   KCS01()}});
        obj.KCS.insert({1, {KCS10(nf), KCS11(nf), KCS12(nf)}});
        obj.KCS.insert({2, {KCS20(nf), KCS21(nf), KCS22(nf), KCS23(nf)}});

        // Matching functions
        const EvolutionBasisQCD evb{nf};

        // PDFs
        obj.MatchingFunctionsPDFs.insert({0, {{evb, C00.at(nf)}}});
        obj.MatchingFunctionsPDFs.insert({1, {{evb, C10pdf.at(nf)}, {evb, C11pdf.at(nf)}, {evb, C12.at(nf)}}});
        obj.MatchingFunctionsPDFs.insert({2, {{evb, ZeroOp}, {evb, ZeroOp}, {evb, ZeroOp}, {evb, ZeroOp}, {evb, ZeroOp}}});
        obj.MatchingFunctionsPDFs.insert({3, {{evb, ZeroOp}}});

        // FFs
        obj.MatchingFunctionsFFs.insert({0, {{evb, C00.at(nf)}}});
        obj.MatchingFunctionsFFs.insert({1, {{evb, ZeroOp}, {evb, ZeroOp}, {evb, ZeroOp}}});
        obj.MatchingFunctionsFFs.insert({2, {{evb, ZeroOp}, {evb, ZeroOp}, {evb, ZeroOp}, {evb, ZeroOp}, {evb, ZeroOp}}});
        obj.MatchingFunctionsFFs.insert({3, {{evb, ZeroOp}}});

        // Hard factors (set to zero when unknown). In addition,
        // H3Ch() should be multiplied by N_{nf,j} =
        // (\sum_{i=1}^{nf}e_q) / e_j channel by channel. Here we set
        // this factor to the constant 1 / 4 because it results from
        // the average ( N_{5,u} + N_{5,d} ) / 2. The numerical
        // difference for different choices is however very small.
        obj.HardFactors.insert({"DY",    {{1, H1DY()},    {2, H2DY(nf)},    {3, H3DY(nf)    + H3Ch() / 4}}});
        obj.HardFactors.insert({"SIDIS", {{1, H1SIDIS()}, {2, H2SIDIS(nf)}, {3, H3SIDIS(nf) + H3Ch() / 4}}});
        obj.HardFactors.insert({"ggH",   {{1, H1ggH()},   {2, H2ggH(nf)},   {3, 0}}});

        // Insert full object
        TmdObj.insert({nf, obj});
      }
    t.stop();

    return TmdObj;
  }

  //_____________________________________________________________________________
  std::function<Set<Distribution>(double const&, double const&, double const&)> BuildTmdPDFs(std::map<int, TmdObjects>                       const& TmdObj,
                                                                                             std::function<Set<Distribution>(double const&)> const& CollPDFs,
                                                                                             std::function<double(double const&)>            const& Alphas,
                                                                                             int                                             const& PerturbativeOrder,
                                                                                             double                                          const& Ci,
                                                                                             double                                          const& IntEps)
  {
    // Match TMDs onto collinear PDFs
    const std::function<Set<Distribution>(double const&)> MatchedTmdPDFs = MatchTmdPDFs(TmdObj, CollPDFs, Alphas, PerturbativeOrder, Ci);

    // Compute TMD evolution factors
    const std::function<std::vector<double>(double const&, double const&, double const&)> EvolFactors = EvolutionFactors(TmdObj, Alphas, PerturbativeOrder, Ci, IntEps);

    // Computed TMDPDFs at the final scale by multiplying the initial
    // scale TMDPDFs by the evolution factor.
    const auto EvolvedTMDs = [=] (double const& b, double const& muf, double const& zetaf) -> Set<Distribution>
    {
      return EvolFactors(b, muf, zetaf) * MatchedTmdPDFs(b);
    };

    return EvolvedTMDs;
  }

  //_____________________________________________________________________________
  std::function<Set<Distribution>(double const&, double const&, double const&)> BuildTmdFFs(std::map<int, TmdObjects>                       const& TmdObj,
                                                                                            std::function<Set<Distribution>(double const&)> const& CollFFs,
                                                                                            std::function<double(double const&)>            const& Alphas,
                                                                                            int                                             const& PerturbativeOrder,
                                                                                            double                                          const& Ci,
                                                                                            double                                          const& IntEps)
  {
    // Match TMDs onto collinear FFs
    const std::function<Set<Distribution>(double const&)> MatchedTmdFFs = MatchTmdFFs(TmdObj, CollFFs, Alphas, PerturbativeOrder, Ci);

    // Compute TMD evolution factors
    const std::function<std::vector<double>(double const&, double const&, double const&)> EvolFactors = EvolutionFactors(TmdObj, Alphas, PerturbativeOrder, Ci, IntEps);

    // Computed TMDFFs at the final scale by multiplying the initial
    // scale TMDFFs by the evolution factor.
    const auto EvolvedTMDs = [=] (double const& b, double const& muf, double const& zetaf) -> Set<Distribution>
    {
      return EvolFactors(b, muf, zetaf) * MatchedTmdFFs(b);
    };

    return EvolvedTMDs;
  }

  //_____________________________________________________________________________
  std::function<double(double const&, double const&, double const&)> BuildTmdJet(std::map<int, TmdObjects>            const& TmdObj,
                                                                                 JetAlgorithm                         const& JetAlgo,
                                                                                 double                               const& JetR,
                                                                                 std::function<double(double const&)> const& Alphas,
                                                                                 int                                  const& PerturbativeOrder,
                                                                                 double                               const& CJ,
                                                                                 double                               const& Ci,
                                                                                 double                               const& IntEps)
  {
    // Stop the code if Ci is different from one. Ci variations not
    // available yet.
    if (Ci != 1)
      throw std::runtime_error(error("BuildTmdJet", "Ci variations not available yet."));

    // Get initial-scale jet TMD
    const std::function<double(double const&, double const&)> MatchedTmdJet = MatchTmdJet(TmdObj, JetAlgo, tan(JetR/2), Alphas, PerturbativeOrder, CJ, Ci, IntEps);

    // Compute TMD evolution factors
    const std::function<double(double const&, double const&, double const&)> EvolFactors = QuarkEvolutionFactor(TmdObj, Alphas, PerturbativeOrder, Ci, IntEps);

    // Computed jet TMD at the final scale by multiplying the initial
    // scale jet TMD by the evolution factor.
    const auto EvolvedTMDs = [=] (double const& b, double const& muf, double const& zetaf) -> double
    {
      return EvolFactors(b, muf, zetaf) * MatchedTmdJet(b, muf);
    };

    return EvolvedTMDs;
  }

  //_____________________________________________________________________________
  std::function<Set<Distribution>(double const&)> MatchTmdPDFs(std::map<int, TmdObjects>                       const& TmdObj,
                                                               std::function<Set<Distribution>(double const&)> const& CollPDFs,
                                                               std::function<double(double const&)>            const& Alphas,
                                                               int                                             const& PerturbativeOrder,
                                                               double                                          const& Ci)
  {
    // Get matching functions
    const std::function<Set<Operator>(double const&)> MatchFunc = MatchingFunctionsPDFs(TmdObj, Alphas, PerturbativeOrder, Ci);

    // Construct function that returns the product of matching
    // functions and collinear PDFs.
    const auto MatchedTMDs = [=] (double const& b) -> Set<Distribution>
    {
      // Define lower scales
      const double mu0 = Ci * 2 * exp(- emc) / b;

      // Convolute matching functions with the collinear PDFs and
      // return.
      return MatchFunc(mu0) * CollPDFs(mu0);
    };

    return MatchedTMDs;
  }

  //_____________________________________________________________________________
  std::function<Set<Distribution>(double const&)> MatchTmdFFs(std::map<int, TmdObjects>                       const& TmdObj,
                                                              std::function<Set<Distribution>(double const&)> const& CollFFs,
                                                              std::function<double(double const&)>            const& Alphas,
                                                              int                                             const& PerturbativeOrder,
                                                              double                                          const& Ci)
  {
    // Get matching functions
    const std::function<Set<Operator>(double const&)> MatchFunc = MatchingFunctionsFFs(TmdObj, Alphas, PerturbativeOrder, Ci);

    // Construct function that returns the product of matching
    // functions and collinear FFs. Includes a factor z^2 typical of
    // FFs.
    const auto MatchedTMDs = [=] (double const& b) -> Set<Distribution>
    {
      // Define lower scales
      const double mu0 = Ci * 2 * exp(- emc) / b;

      // Convolute matching functions with the collinear FFs and
      // return.
      return MatchFunc(mu0) * CollFFs(mu0);
    };

    return MatchedTMDs;
  }

  //_____________________________________________________________________________
  std::function<double(double const&, double const&)> MatchTmdJet(std::map<int, TmdObjects>            const& TmdObj,
                                                                  JetAlgorithm                         const& JetAlgo,
                                                                  double                               const& tR,
                                                                  std::function<double(double const&)> const& Alphas,
                                                                  int                                  const& PerturbativeOrder,
                                                                  double                               const& CJ,
                                                                  double                               const& Ci,
                                                                  double                               const& IntEps)
  {
    // Retrieve thresholds from "TmdObj"
    std::vector<double> thrs;
    for(auto const& obj : TmdObj)
      {
        const int    nf  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) thrs.size() < nf)
          thrs.resize(nf);
        thrs[nf-1] = thr;
      }

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the threshold
    // vector entries are ordered).
    int nfi = 0;
    int nff = thrs.size();
    for (auto const& v : thrs)
      if (v <= 0)
        nfi++;

    // Compute log and its powers
    const double lQ  = log(CJ);
    const double lQ2 = lQ * lQ;

    // Select coefficients according to the jet algorithm
    std::map<int, double> gK0;
    std::map<int, double> gK1;
    std::map<int, double> gK2;
    std::map<int, double> gF0;
    std::map<int, double> gF1;
    for (int nf = nfi; nf <= nff; nf++)
      {
        gK0.insert({nf, CF * TmdObj.at(nf).GammaK.at(0)});
        gK1.insert({nf, CF * TmdObj.at(nf).GammaK.at(1)});
        gK2.insert({nf, CF * TmdObj.at(nf).GammaK.at(2)});
        gF0.insert({nf, TmdObj.at(nf).GammaFq.at(0)});
        gF1.insert({nf, TmdObj.at(nf).GammaFq.at(1)});
      }

    double djet1 = 0;
    if (JetAlgo == CONE)
      djet1 = dJetqCone1();
    else if (JetAlgo == KT)
      djet1 = dJetqkT1();
    else
      throw std::runtime_error(error("MatchTmdJet", "Unknown jet algorithm."));

    // Construct function that returns the product of matching
    // functions and collinear FFs. Includes a factor z^2 typical of
    // FFs.
    const auto MatchedTMDs = [=] (double const& b, double const& mu) -> double
    {
      // Define jet scale
      const double muJ = CJ * tR * mu;

      // The strong coupling at muJ
      const double coup = Alphas(muJ) / FourPi;

      // Number of active flavours
      const int nf = NF(mu, thrs);

      // Compute low-scale jet function
      double J = 1;
      if (PerturbativeOrder > 1 || PerturbativeOrder < 0)
        J += coup * ( djet1 + gF0.at(nf) * lQ + gK0.at(nf) * lQ2 / 2 );

      // Define integrand of the evolution factor
      const Integrator gammaJ{
        [=] (double const& mup) -> double
        {
          // The strong coupling at mup
          const double coup = Alphas(mup) / FourPi;

          // Log of the scales
          const double L = log(muJ / mup);

          // Number of active flavours
          const int nfp = NF(mup, thrs);

          // LL
          double gJ = - coup * gK0.at(nfp) * L;

          // NLL
          if (PerturbativeOrder != 0)
            gJ += coup * ( gF0.at(nfp) - coup * gK1.at(nfp) * L );
          // NNLL
          if (PerturbativeOrder != 1)
            gJ += coup * coup * ( gF1.at(nfp)  - coup * gK2.at(nfp)  * L );

          return gJ / mup;
        }
      };

      // Define lower scales
      const double mu0 = Ci * 2 * exp(- emc) / b;

      // Convolute matching functions with the collinear FFs and
      // return.
      return J * exp(- gammaJ.integrate(mu0, muJ, IntEps) );
    };

    return MatchedTMDs;
  }

  //_____________________________________________________________________________
  std::function<Set<Operator>(double const&)> MatchingFunctionsPDFs(std::map<int, TmdObjects>            const& TmdObj,
                                                                    std::function<double(double const&)> const& Alphas,
                                                                    int                                  const& PerturbativeOrder,
                                                                    double                               const& Ci)
  {
    // Retrieve thresholds from "TmdObj"
    std::vector<double> thrs;
    for(auto const& obj : TmdObj)
      {
        const int    nf  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) thrs.size() < nf)
          thrs.resize(nf);
        thrs[nf-1] = thr;
      }

    // Define the log(Ci) to assess scale variations
    const double Lmu = log(Ci);

    // Matching functions as functions of the absolute value of the
    // impact parameter b.
    std::function<Set<Operator>(double const&)> MatchFunc;
    if (PerturbativeOrder == LL || PerturbativeOrder == NLL)
      MatchFunc = [=] (double const& mu) -> Set<Operator>
      {
        return TmdObj.at(NF(mu, thrs)).MatchingFunctionsPDFs.at(0)[0];
      };
    else if (PerturbativeOrder == NNLL || PerturbativeOrder == NLLp)
      MatchFunc = [=] (double const& mu) -> Set<Operator>
      {
        const double coup = Alphas(mu) / FourPi;
        const auto& mf = TmdObj.at(NF(mu, thrs)).MatchingFunctionsPDFs;
        const auto c0  = mf.at(0);
        const auto c1  = mf.at(1);
        const auto lo  = c0[0];
        const auto nlo = c1[0] + Lmu * ( c1[1] + Lmu * c1[2] );
        return lo + coup * nlo;
      };
    else if (PerturbativeOrder == NNNLL || PerturbativeOrder == NNLLp)
      MatchFunc = [=] (double const& mu) -> Set<Operator>
      {
        const double coup = Alphas(mu) / FourPi;
        const auto& mf  = TmdObj.at(NF(mu, thrs)).MatchingFunctionsPDFs;
        const auto c0   = mf.at(0);
        const auto c1   = mf.at(1);
        const auto c2   = mf.at(2);
        const auto lo   = c0[0];
        const auto nlo  = c1[0] + Lmu * ( c1[1] + Lmu * c1[2] );
        const auto nnlo = c2[0] + Lmu * ( c2[1] + Lmu * ( c2[2] + Lmu * ( c2[3] + Lmu * c2[4] ) ) );
        return lo + coup * ( nlo + coup * nnlo );
      };
    else if (PerturbativeOrder == NNNLLp)
      MatchFunc = [=] (double const& mu) -> Set<Operator>
      {
        const double coup = Alphas(mu) / FourPi;
        const auto& mf   = TmdObj.at(NF(mu, thrs)).MatchingFunctionsPDFs;
        const auto c0    = mf.at(0);
        const auto c1    = mf.at(1);
        const auto c2    = mf.at(2);
        const auto lo    = c0[0];
        const auto nlo   = c1[0] + Lmu * ( c1[1] + Lmu * c1[2] );
        const auto nnlo  = c2[0] + Lmu * ( c2[1] + Lmu * ( c2[2] + Lmu * ( c2[3] + Lmu * c2[4] ) ) );
        const auto nnnlo = mf.at(3)[0];
        return lo + coup * ( nlo + coup * ( nnlo + coup * nnnlo ) );
      };

    // Return matching functions
    return MatchFunc;
  }

  //_____________________________________________________________________________
  std::function<Set<Operator>(double const&)> MatchingFunctionsFFs(std::map<int, TmdObjects>             const& TmdObj,
                                                                   std::function<double(double const&)>  const& Alphas,
                                                                   int                                   const& PerturbativeOrder,
                                                                   double                                const& Ci)
  {
    // Retrieve thresholds from "TmdObj"
    std::vector<double> thrs;
    for(auto const& obj : TmdObj)
      {
        const int    nf  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) thrs.size() < nf)
          thrs.resize(nf);
        thrs[nf-1] = thr;
      }

    // Define the log(Ci) to assess scale variations
    const double Lmu = log(Ci);

    // Matching functions as functions of the absolute value of the
    // impact parameter b.
    std::function<Set<Operator>(double const&)> MatchFunc;
    if (PerturbativeOrder == LL || PerturbativeOrder == NLL)
      MatchFunc = [=] (double const& mu) -> Set<Operator>
      {
        return TmdObj.at(NF(mu, thrs)).MatchingFunctionsFFs.at(0)[0];
      };
    else if (PerturbativeOrder == NNLL || PerturbativeOrder == NLLp)
      MatchFunc = [=] (double const& mu) -> Set<Operator>
      {
        const double coup = Alphas(mu) / FourPi;
        const auto& mf = TmdObj.at(NF(mu, thrs)).MatchingFunctionsFFs;
        const auto c0  = mf.at(0);
        const auto c1  = mf.at(1);
        const auto lo  = c0[0];
        const auto nlo = c1[0] + Lmu * ( c1[1] + Lmu * c1[2] );
        return lo + coup * nlo;
      };
    else if (PerturbativeOrder == NNNLL || PerturbativeOrder == NNLLp)
      MatchFunc = [=] (double const& mu) -> Set<Operator>
      {
        const double coup = Alphas(mu) / FourPi;
        const auto& mf  = TmdObj.at(NF(mu, thrs)).MatchingFunctionsFFs;
        const auto c0   = mf.at(0);
        const auto c1   = mf.at(1);
        const auto c2   = mf.at(2);
        const auto lo   = c0[0];
        const auto nlo  = c1[0] + Lmu * ( c1[1] + Lmu * c1[2] );
        const auto nnlo = c2[0] + Lmu * ( c2[1] + Lmu * ( c2[2] + Lmu * ( c2[3] + Lmu * c2[4] ) ) );
        return lo + coup * ( nlo + coup * nnlo );
      };
    else if (PerturbativeOrder == NNNLLp)
      MatchFunc = [=] (double const& mu) -> Set<Operator>
      {
        const double coup = Alphas(mu) / FourPi;
        const auto& mf   = TmdObj.at(NF(mu, thrs)).MatchingFunctionsFFs;
        const auto c0    = mf.at(0);
        const auto c1    = mf.at(1);
        const auto c2    = mf.at(2);
        const auto lo    = c0[0];
        const auto nlo   = c1[0] + Lmu * ( c1[1] + Lmu * c1[2] );
        const auto nnlo  = c2[0] + Lmu * ( c2[1] + Lmu * ( c2[2] + Lmu * ( c2[3] + Lmu * c2[4] ) ) );
        const auto nnnlo = mf.at(3)[0];
        return lo + coup * ( nlo + coup * ( nnlo + coup * nnnlo ) );
      };

    // Return matching functions
    return MatchFunc;
  }

  //_____________________________________________________________________________
  std::function<std::vector<double>(double const&, double const&, double const&)> EvolutionFactors(std::map<int, TmdObjects>            const& TmdObj,
                                                                                                   std::function<double(double const&)> const& Alphas,
                                                                                                   int                                  const& PerturbativeOrder,
                                                                                                   double                               const& Ci,
                                                                                                   double                               const& IntEps)
  {
    // Retrieve thresholds from "TmdObj"
    std::vector<double> thrs;
    for(auto const& obj : TmdObj)
      {
        const int    nf  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) thrs.size() < nf)
          thrs.resize(nf);
        thrs[nf-1] = thr;
      }

    // Define the log(Ci) to assess scale variations
    const double Lmu = log(Ci);

    // Create functions needed for the TMD evolution
    std::function<double(double const&)> gammaFq;
    std::function<double(double const&)> gammaFg;
    std::function<double(double const&)> gammaK;
    std::function<double(double const&)> K;
    // LL
    if (PerturbativeOrder == LL)
      {
        gammaFq = [=] (double const&) -> double{ return 0; };
        gammaFg = [=] (double const&) -> double{ return 0; };
        gammaK  = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * TmdObj.at(NF(mu, thrs)).GammaK.at(0);
        };
        K = [=] (double const&) -> double{ return 0; };
      }
    // NLL
    else if (PerturbativeOrder == NLL || PerturbativeOrder == NLLp)
      {
        gammaFq = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * TmdObj.at(NF(mu, thrs)).GammaFq.at(0);
        };
        gammaFg = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * TmdObj.at(NF(mu, thrs)).GammaFg.at(0);
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * gc.at(1) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          return coup * lo;
        };
      }
    // NNLL
    else if (PerturbativeOrder == NNLL || PerturbativeOrder == NNLLp)
      {
        gammaFq = [=] (double const& mu) -> double
        {
          const auto& gv    = TmdObj.at(NF(mu, thrs)).GammaFq;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * gv.at(1) );
        };
        gammaFg = [=] (double const& mu) -> double
        {
          const auto& gv    = TmdObj.at(NF(mu, thrs)).GammaFg;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * gv.at(1) );
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * gc.at(2) ) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const std::vector<double> d1 = d.at(1);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          const double nlo  = d1[0] + Lmu * ( d1[1] + Lmu * d1[2] );
          return coup * ( lo + coup * nlo );
        };
      }
    // N3LL
    else if (PerturbativeOrder == NNNLL || PerturbativeOrder == NNNLLp)
      {
        gammaFq = [=] (double const& mu) -> double
        {
          const auto& gv    = TmdObj.at(NF(mu, thrs)).GammaFq;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * ( gv.at(1) + coup * gv.at(2) ) );
        };
        gammaFg = [=] (double const& mu) -> double
        {
          const auto& gv    = TmdObj.at(NF(mu, thrs)).GammaFg;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * ( gv.at(1) + coup * gv.at(2) ) );
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * ( gc.at(2) + coup * gc.at(3) ) ) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const std::vector<double> d1 = d.at(1);
          const std::vector<double> d2 = d.at(2);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          const double nlo  = d1[0] + Lmu * ( d1[1] + Lmu * d1[2] );
          const double nnlo = d2[0] + Lmu * ( d2[1] + Lmu * ( d2[2] + Lmu * d2[3] ) );
          return coup * ( lo + coup * ( nlo + coup * nnlo ) );
        };
      }

    // Define the integrands
    const Integrator I1q{[=] (double const& mu) -> double{ return gammaFq(mu) / mu; }};
    const Integrator I1g{[=] (double const& mu) -> double{ return gammaFg(mu) / mu; }};
    const Integrator I2 {[=] (double const& mu) -> double{ return gammaK(mu) / mu; }};
    const Integrator I3 {[=] (double const& mu) -> double{ return gammaK(mu) * log(mu) / mu; }};

    // Construct function that returns the perturbative evolution
    // kernel.
    const auto EvolFactors = [=] (double const& b, double const& muf, double const& zetaf) -> std::vector<double>
    {
      // Define lower scales
      const double mu0   = Ci * 2 * exp(- emc) / b;
      const double zeta0 = mu0 * mu0;

      // Compute argument of the exponent of the evolution factors
      const double IntI1q = I1q.integrate(mu0, muf, thrs, IntEps);
      const double IntI1g = I1g.integrate(mu0, muf, thrs, IntEps);
      const double IntI2  = I2.integrate(mu0, muf, thrs, IntEps) * log(zetaf);
      const double IntI3  = I3.integrate(mu0, muf, thrs, IntEps);

      // Compute the evolution factors
      const double Klz = ( K(mu0) * log( zetaf / zeta0 ) - IntI2 ) / 2 + IntI3;
      const double Rq  = exp( CF * Klz + IntI1q );
      const double Rg  = exp( CA * Klz + IntI1g );

      // Return vector of evolution factors
      return std::vector<double>{Rg, Rq, Rq, Rq, Rq, Rq, Rq, Rq, Rq, Rq, Rq, Rq, Rq};
    };

    return EvolFactors;
  }

  //_____________________________________________________________________________
  std::function<std::vector<double>(double const&, double const&, double const&)> EvolutionFactorsK(std::map<int, TmdObjects>            const& TmdObj,
                                                                                                    std::function<double(double const&)> const& Alphas,
                                                                                                    int                                  const& PerturbativeOrder,
                                                                                                    double                               const& Ci,
                                                                                                    double                               const& IntEps)
  {
    // Retrieve thresholds from "TmdObj"
    std::vector<double> thrs;
    for(auto const& obj : TmdObj)
      {
        const int    nf  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) thrs.size() < nf)
          thrs.resize(nf);
        thrs[nf-1] = thr;
      }

    // Define the log(Ci) to assess scale variations
    const double Lmu = log(Ci);

    // Create functions needed for the TMD evolution
    std::function<double(double const&)> gammaFq;
    std::function<double(double const&)> gammaFg;
    std::function<double(double const&)> gammaK1;
    std::function<double(double const&)> gammaK2;
    std::function<double(double const&)> K;
    // LL
    if (PerturbativeOrder == LL)
      {
        gammaFq = [=] (double const&) -> double{ return 0; };
        gammaFg = [=] (double const&) -> double{ return 0; };
        gammaK1 = [=] (double const&) -> double{ return 0; };
        gammaK2 = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * TmdObj.at(NF(mu, thrs)).GammaK.at(0);
        };
        K = [=] (double const&) -> double{ return 0; };
      }
    // NLL
    else if (PerturbativeOrder == NLL || PerturbativeOrder == NLLp)
      {
        gammaFq = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * TmdObj.at(NF(mu, thrs)).GammaFq.at(0);
        };
        gammaFg = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * TmdObj.at(NF(mu, thrs)).GammaFg.at(0);
        };
        gammaK1 = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * TmdObj.at(NF(mu, thrs)).GammaK.at(0);
        };
        gammaK2 = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * gc.at(1) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          return coup * lo;
        };
      }
    // NNLL
    else if (PerturbativeOrder == NNLL || PerturbativeOrder == NNLLp)
      {
        gammaFq = [=] (double const& mu) -> double
        {
          const auto& gv    = TmdObj.at(NF(mu, thrs)).GammaFq;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * gv.at(1) );
        };
        gammaFg = [=] (double const& mu) -> double
        {
          const auto& gv    = TmdObj.at(NF(mu, thrs)).GammaFg;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * gv.at(1) );
        };
        gammaK1 = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * gc.at(1) );
        };
        gammaK2 = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * gc.at(2) ) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const std::vector<double> d1 = d.at(1);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          const double nlo  = d1[0] + Lmu * ( d1[1] + Lmu * d1[2] );
          return coup * ( lo + coup * nlo );
        };
      }
    // N3LL
    else if (PerturbativeOrder == NNNLL || PerturbativeOrder == NNNLLp)
      {
        gammaFq = [=] (double const& mu) -> double
        {
          const auto& gv    = TmdObj.at(NF(mu, thrs)).GammaFq;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * ( gv.at(1) + coup * gv.at(2) ) );
        };
        gammaFg = [=] (double const& mu) -> double
        {
          const auto& gv    = TmdObj.at(NF(mu, thrs)).GammaFg;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * ( gv.at(1) + coup * gv.at(2) ) );
        };
        gammaK1 = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * gc.at(2) ) );
        };
        gammaK2 = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * ( gc.at(2) + coup * gc.at(3) ) ) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const std::vector<double> d1 = d.at(1);
          const std::vector<double> d2 = d.at(2);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          const double nlo  = d1[0] + Lmu * ( d1[1] + Lmu * d1[2] );
          const double nnlo = d2[0] + Lmu * ( d2[1] + Lmu * ( d2[2] + Lmu * d2[3] ) );
          return coup * ( lo + coup * ( nlo + coup * nnlo ) );
        };
      }

    // Define the integrands
    const Integrator I1q{[=] (double const& mu) -> double{ return gammaFq(mu) / mu; }};
    const Integrator I1g{[=] (double const& mu) -> double{ return gammaFg(mu) / mu; }};
    const Integrator I2 {[=] (double const& mu) -> double{ return gammaK2(mu) / mu; }};
    const Integrator I3 {[=] (double const& mu) -> double{ return gammaK2(mu) * log(mu) / mu; }};
    const Integrator I4 {[=] (double const& mu) -> double{ return gammaK1(mu) / mu; }};

    // Construct function that returns the perturbative evolution
    // kernel.
    const auto EvolFactors = [=] (double const& b, double const& muf, double const& zetaf) -> std::vector<double>
    {
      // Define lower scales
      const double mu0   = Ci * 2 * exp(- emc) / b;
      const double zeta0 = mu0 * mu0;

      // Compute argument of the exponent of the evolution factors
      const double IntI1q = I1q.integrate(mu0, muf, thrs, IntEps);
      const double IntI1g = I1g.integrate(mu0, muf, thrs, IntEps);
      const double IntI2  = I2.integrate(mu0, muf, thrs, IntEps) * log(muf);
      const double IntI4  = I4.integrate(mu0, muf, thrs, IntEps) * log(zetaf / pow(muf, 2));
      const double IntI3  = I3.integrate(mu0, muf, thrs, IntEps);

      // Compute the evolution factors
      const double Klz = ( K(mu0) * log( zetaf / zeta0 ) - IntI4 ) / 2 - IntI2 + IntI3;
      const double Rq  = exp( CF * Klz + IntI1q );
      const double Rg  = exp( CA * Klz + IntI1g );

      // Return vector of evolution factors
      return std::vector<double>{Rg, Rq, Rq, Rq, Rq, Rq, Rq, Rq, Rq, Rq, Rq, Rq, Rq};
    };

    return EvolFactors;
  }

  //_____________________________________________________________________________
  std::function<double(double const&, double const&, double const&)> QuarkEvolutionFactor(std::map<int, TmdObjects>            const& TmdObj,
                                                                                          std::function<double(double const&)> const& Alphas,
                                                                                          int                                  const& PerturbativeOrder,
                                                                                          double                               const& Ci,
                                                                                          double                               const& IntEps)
  {
    // Retrieve thresholds from "TmdObj"
    std::vector<double> thrs;
    for(auto const& obj : TmdObj)
      {
        const int    nf  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) thrs.size() < nf)
          thrs.resize(nf);
        thrs[nf-1] = thr;
      }

    // Define the log(Ci) to assess scale variations
    const double Lmu = log(Ci);

    // Create functions needed for the TMD evolution
    std::function<double(double const&)> gammaFq;
    std::function<double(double const&)> gammaK;
    std::function<double(double const&)> K;
    // LL
    if (PerturbativeOrder == LL)
      {
        gammaFq = [=] (double const&) -> double{ return 0; };
        gammaK  = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * TmdObj.at(NF(mu, thrs)).GammaK.at(0);
        };
        K = [=] (double const&) -> double{ return 0; };
      }
    // NLL
    else if (PerturbativeOrder == NLL || PerturbativeOrder == NLLp)
      {
        gammaFq = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * TmdObj.at(NF(mu, thrs)).GammaFq.at(0);
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * gc.at(1) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          return coup * lo;
        };
      }
    // NNLL
    else if (PerturbativeOrder == NNLL || PerturbativeOrder == NNLLp)
      {
        gammaFq = [=] (double const& mu) -> double
        {
          const auto& gv    = TmdObj.at(NF(mu, thrs)).GammaFq;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * gv.at(1) );
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * gc.at(2) ) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const std::vector<double> d1 = d.at(1);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          const double nlo  = d1[0] + Lmu * ( d1[1] + Lmu * d1[2] );
          return coup * ( lo + coup * nlo );
        };
      }
    // N3LL
    else if (PerturbativeOrder == NNNLL || PerturbativeOrder == NNNLLp)
      {
        gammaFq = [=] (double const& mu) -> double
        {
          const auto& gv    = TmdObj.at(NF(mu, thrs)).GammaFq;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * ( gv.at(1) + coup * gv.at(2) ) );
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * ( gc.at(2) + coup * gc.at(3) ) ) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const std::vector<double> d1 = d.at(1);
          const std::vector<double> d2 = d.at(2);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          const double nlo  = d1[0] + Lmu * ( d1[1] + Lmu * d1[2] );
          const double nnlo = d2[0] + Lmu * ( d2[1] + Lmu * ( d2[2] + Lmu * d2[3] ) );
          return coup * ( lo + coup * ( nlo + coup * nnlo ) );
        };
      }

    // Define the integrands
    const Integrator I1{[=] (double const& mu) -> double{ return gammaFq(mu) / mu; }};
    const Integrator I2{[=] (double const& mu) -> double{ return gammaK(mu) / mu; }};
    const Integrator I3{[=] (double const& mu) -> double{ return gammaK(mu) * log(mu) / mu; }};

    // Construct function that returns the perturbative evolution
    // kernel.
    const auto EvolFactor = [=] (double const& b, double const& muf, double const& zetaf) -> double
    {
      // Define lower scales
      const double mu0   = Ci * 2 * exp(- emc) / b;
      const double zeta0 = mu0 * mu0;

      // Compute argument of the exponent of the evolution factors
      const double IntI1 = I1.integrate(mu0, muf, thrs, IntEps);
      const double IntI2 = I2.integrate(mu0, muf, thrs, IntEps) * log(zetaf);
      const double IntI3 = I3.integrate(mu0, muf, thrs, IntEps);

      // Compute the evolution factors
      const double Klz = ( K(mu0) * log( zetaf / zeta0 ) - IntI2 ) / 2 + IntI3;
      const double Rq  = exp( CF * Klz + IntI1 );

      // Return the evolution factor
      return Rq;
    };

    return EvolFactor;
  }

  //_____________________________________________________________________________
  std::function<double(double const&, double const&, double const&)> GluonEvolutionFactor(std::map<int, TmdObjects>            const& TmdObj,
                                                                                          std::function<double(double const&)> const& Alphas,
                                                                                          int                                  const& PerturbativeOrder,
                                                                                          double                               const& Ci,
                                                                                          double                               const& IntEps)
  {
    // Retrieve thresholds from "TmdObj"
    std::vector<double> thrs;
    for(auto const& obj : TmdObj)
      {
        const int    nf  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) thrs.size() < nf)
          thrs.resize(nf);
        thrs[nf-1] = thr;
      }

    // Define the log(Ci) to assess scale variations
    const double Lmu = log(Ci);

    // Create functions needed for the TMD evolution
    std::function<double(double const&)> gammaFg;
    std::function<double(double const&)> gammaK;
    std::function<double(double const&)> K;
    // LL
    if (PerturbativeOrder == LL)
      {
        gammaFg = [=] (double const&) -> double{ return 0; };
        gammaK  = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * TmdObj.at(NF(mu, thrs)).GammaK.at(0);
        };
        K = [=] (double const&) -> double{ return 0; };
      }
    // NLL
    else if (PerturbativeOrder == NLL || PerturbativeOrder == NLLp)
      {
        gammaFg = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * TmdObj.at(NF(mu, thrs)).GammaFg.at(0);
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * gc.at(1) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          return coup * lo;
        };
      }
    // NNLL
    else if (PerturbativeOrder == NNLL || PerturbativeOrder == NNLLp)
      {
        gammaFg = [=] (double const& mu) -> double
        {
          const auto& gv    = TmdObj.at(NF(mu, thrs)).GammaFg;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * gv.at(1) );
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * gc.at(2) ) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const std::vector<double> d1 = d.at(1);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          const double nlo  = d1[0] + Lmu * ( d1[1] + Lmu * d1[2] );
          return coup * ( lo + coup * nlo );
        };
      }
    // N3LL
    else if (PerturbativeOrder == NNNLL || PerturbativeOrder == NNNLLp)
      {
        gammaFg = [=] (double const& mu) -> double
        {
          const auto& gv    = TmdObj.at(NF(mu, thrs)).GammaFg;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * ( gv.at(1) + coup * gv.at(2) ) );
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * ( gc.at(2) + coup * gc.at(3) ) ) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const std::vector<double> d1 = d.at(1);
          const std::vector<double> d2 = d.at(2);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          const double nlo  = d1[0] + Lmu * ( d1[1] + Lmu * d1[2] );
          const double nnlo = d2[0] + Lmu * ( d2[1] + Lmu * ( d2[2] + Lmu * d2[3] ) );
          return coup * ( lo + coup * ( nlo + coup * nnlo ) );
        };
      }

    // Define the integrands
    const Integrator I1{[=] (double const& mu) -> double{ return gammaFg(mu) / mu; }};
    const Integrator I2{[=] (double const& mu) -> double{ return gammaK(mu) / mu; }};
    const Integrator I3{[=] (double const& mu) -> double{ return gammaK(mu) * log(mu) / mu; }};

    // Construct function that returns the perturbative evolution
    // kernel.
    const auto EvolFactor = [=] (double const& b, double const& muf, double const& zetaf) -> double
    {
      // Define lower scales
      const double mu0   = Ci * 2 * exp(- emc) / b;
      const double zeta0 = mu0 * mu0;

      // Compute argument of the exponent of the evolution factors
      const double IntI1 = I1.integrate(mu0, muf, thrs, IntEps);
      const double IntI2 = I2.integrate(mu0, muf, thrs, IntEps) * log(zetaf);
      const double IntI3 = I3.integrate(mu0, muf, thrs, IntEps);

      // Compute the evolution factors
      const double Klz = ( K(mu0) * log( zetaf / zeta0 ) - IntI2 ) / 2 + IntI3;
      const double Rg  = exp( CA * Klz + IntI1 );

      // Return the evolution factor
      return Rg;
    };

    return EvolFactor;
  }

  //_____________________________________________________________________________
  std::function<double(double const&, double const&)> CollinsSoperKernel(std::map<int, TmdObjects>            const& TmdObj,
                                                                         std::function<double(double const&)> const& Alphas,
                                                                         int                                  const& PerturbativeOrder,
                                                                         double                               const& Ci,
                                                                         double                               const& IntEps)
  {
    // Retrieve thresholds from "TmdObj"
    std::vector<double> thrs;
    for(auto const& obj : TmdObj)
      {
        const int    nf  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) thrs.size() < nf)
          thrs.resize(nf);
        thrs[nf-1] = thr;
      }

    // Define the log(Ci) to assess scale variations
    const double Lmu = log(Ci);

    // Create functions needed for the TMD evolution
    std::function<double(double const&)> gammaK;
    std::function<double(double const&)> K;
    // LL
    if (PerturbativeOrder == LL)
      {
        gammaK  = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * TmdObj.at(NF(mu, thrs)).GammaK.at(0);
        };
        K = [=] (double const&) -> double{ return 0; };
      }
    // NLL
    else if (PerturbativeOrder == NLL || PerturbativeOrder == NLLp)
      {
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * gc.at(1) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          return coup * lo;
        };
      }
    // NNLL
    else if (PerturbativeOrder == NNLL || PerturbativeOrder == NNLLp)
      {
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * gc.at(2) ) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const std::vector<double> d1 = d.at(1);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          const double nlo  = d1[0] + Lmu * ( d1[1] + Lmu * d1[2] );
          return coup * ( lo + coup * nlo );
        };
      }
    // N3LL
    else if (PerturbativeOrder == NNNLL || PerturbativeOrder == NNNLLp)
      {
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu, thrs)).GammaK;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * ( gc.at(2) + coup * gc.at(3) ) ) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu, thrs)).KCS;
          const std::vector<double> d0 = d.at(0);
          const std::vector<double> d1 = d.at(1);
          const std::vector<double> d2 = d.at(2);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          const double nlo  = d1[0] + Lmu * ( d1[1] + Lmu * d1[2] );
          const double nnlo = d2[0] + Lmu * ( d2[1] + Lmu * ( d2[2] + Lmu * d2[3] ) );
          return coup * ( lo + coup * ( nlo + coup * nnlo ) );
        };
      }

    // Define the integrands
    const Integrator I2{[=] (double const& mu) -> double{ return gammaK(mu) / mu; }};

    // Construct function that returns the perturbative evolution
    // kernel.
    const auto CSKernel = [=] (double const& b, double const& muf) -> double
    {
      // Define lower scales
      const double mu0 = Ci * 2 * exp(- emc) / b;

      // Compute argument of the exponent of the evolution factors
      const double IntI2 = I2.integrate(mu0, muf, thrs, IntEps);

      // Compute the evolution factors
      const double Klz = CF * ( K(mu0) - IntI2 );

      // Return the evolution factor
      return Klz;
    };

    return CSKernel;
  }

  //_____________________________________________________________________________
  std::function<double(double const&)> HardFactor(std::string                          const& Process,
                                                  std::map<int, TmdObjects>            const& TmdObj,
                                                  std::function<double(double const&)> const& Alphas,
                                                  int                                  const& PerturbativeOrder,
                                                  double                               const& Cf)
  {
    // Retrieve thresholds from "TmdObj"
    std::vector<double> thrs;
    for(auto const& obj : TmdObj)
      {
        const int    nf  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) thrs.size() < nf)
          thrs.resize(nf);
        thrs[nf-1] = thr;
      }

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the threshold
    // vector entries are ordered).
    int nfi = 0;
    int nff = thrs.size();
    for (auto const& v : thrs)
      if (v <= 0)
        nfi++;

    // Compute log and its powers
    const double lQ  = log(Cf);
    const double lQ2 = lQ * lQ;
    const double lQ3 = lQ * lQ2;
    const double lQ4 = lQ * lQ3;
    const double lQ5 = lQ * lQ4;
    const double lQ6 = lQ * lQ5;

    // Select coefficients according to the process
    std::map<int, double> b0;
    std::map<int, double> b1;
    std::map<int, double> gK0;
    std::map<int, double> gK1;
    std::map<int, double> gK2;
    std::map<int, double> gF0;
    std::map<int, double> gF1;
    std::map<int, double> gF2;
    std::map<int, double> H1;
    std::map<int, double> H2;
    std::map<int, double> H3;
    if (Process == "DY" || Process == "SIDIS")
      for (int nf = nfi; nf <= nff; nf++)
        {
          b0.insert({nf, TmdObj.at(nf).Beta.at(0)});
          b1.insert({nf, TmdObj.at(nf).Beta.at(1)});
          gK0.insert({nf, CF * TmdObj.at(nf).GammaK.at(0)});
          gK1.insert({nf, CF * TmdObj.at(nf).GammaK.at(1)});
          gK2.insert({nf, CF * TmdObj.at(nf).GammaK.at(2)});
          gF0.insert({nf, TmdObj.at(nf).GammaFq.at(0)});
          gF1.insert({nf, TmdObj.at(nf).GammaFq.at(1)});
          gF2.insert({nf, TmdObj.at(nf).GammaFq.at(2)});
          H1.insert({nf, TmdObj.at(nf).HardFactors.at(Process).at(1)});
          H2.insert({nf, TmdObj.at(nf).HardFactors.at(Process).at(2)});
          H3.insert({nf, TmdObj.at(nf).HardFactors.at(Process).at(3)});
        }
    else if (Process == "ggH")
      for (int nf = nfi; nf <= nff; nf++)
        {
          b0.insert({nf, TmdObj.at(nf).Beta.at(0)});
          b1.insert({nf, TmdObj.at(nf).Beta.at(1)});
          gK0.insert({nf, CA * TmdObj.at(nf).GammaK.at(0)});
          gK1.insert({nf, CA * TmdObj.at(nf).GammaK.at(1)});
          gK2.insert({nf, CA * TmdObj.at(nf).GammaK.at(2)});
          gF0.insert({nf, TmdObj.at(nf).GammaFg.at(0)});
          gF1.insert({nf, TmdObj.at(nf).GammaFg.at(1)});
          gF2.insert({nf, TmdObj.at(nf).GammaFg.at(2)});
          H1.insert({nf, TmdObj.at(nf).HardFactors.at(Process).at(1)});
          H2.insert({nf, TmdObj.at(nf).HardFactors.at(Process).at(2)});
          H3.insert({nf, TmdObj.at(nf).HardFactors.at(Process).at(3)});
        }
    else
      throw std::runtime_error(error("HardFactor", "Process not available."));

    // Construct function that returns the hard function
    const auto HardFactor = [=] (double const& Q) -> double
    {
      // The strong coupling
      const double coup = Alphas(Q) / FourPi;

      // Number of active flavours
      const int nf = NF(Q, thrs);

      // Compute hard function
      double H = 1;
      if (PerturbativeOrder > 1 || PerturbativeOrder < 0)
        H += coup * ( H1.at(nf)
                      - 2 * gF0.at(nf) * lQ - gK0.at(nf) * lQ2 );
      if (PerturbativeOrder > 2 || PerturbativeOrder < -1)
        H += pow(coup, 2) * ( H2.at(nf)
                              + ( - 2 * H1.at(nf) * gF0.at(nf) - 2 * gF1.at(nf) ) * lQ
                              + ( - gK0.at(nf) * H1.at(nf) + 2 * b0.at(nf) * gF0.at(nf) + 2 * pow(gF0.at(nf), 2) - gK1.at(nf) ) * lQ2
                              + ( 4 * b0.at(nf) * gK0.at(nf) / 3 + 2 * gF0.at(nf) * gK0.at(nf) ) * lQ3
                              + pow(gK0.at(nf), 2) / 2 * lQ4 );
      if (PerturbativeOrder > 3 || PerturbativeOrder < -2)
        H += pow(coup, 3) * ( H3.at(nf)
                              + ( - 2 * gF0.at(nf) * H2.at(nf) - 2 * gF1.at(nf) * H1.at(nf) - 2 * gF2.at(nf) ) * lQ
                              + ( 2 * b0.at(nf) * gF0.at(nf) * H1.at(nf) + 2 * pow(gF0.at(nf), 2) * H1.at(nf)
                                  + 4 * gF0.at(nf) * gF1.at(nf) - gK0.at(nf) * H2.at(nf) + 2 * b1.at(nf) * gF0.at(nf)
                                  + 4 * b0.at(nf) * gF1.at(nf) - H1.at(nf) * gK1.at(nf) - gK2.at(nf) ) * lQ2
                              + ( 2 * gF0.at(nf) * gK0.at(nf) * H1.at(nf) + 2 * gF0.at(nf) * gK1.at(nf)
                                  + 2 * gK0.at(nf) * gF1.at(nf) + 4 * b0.at(nf) * gK0.at(nf) * H1.at(nf) / 3
                                  - 4 * b0.at(nf) * pow(gF0.at(nf), 2) - 8 * pow(b0.at(nf), 2) * gF0.at(nf) / 3
                                  - 4 * pow(gF0.at(nf), 3) / 3 + 4 * b1.at(nf) * gK0.at(nf) / 3 + 8 * b0.at(nf) * gK1.at(nf) ) * lQ3
                              + ( pow(gK0.at(nf), 2) * H1.at(nf) / 2 + gK0.at(nf) * gK1.at(nf)
                                  - 14 * b0.at(nf) * gF0.at(nf) * gK0.at(nf) / 3 - 2 * pow(gF0.at(nf), 2) * gK0.at(nf)
                                  - 2 * pow(b0.at(nf), 2) * gK0.at(nf) ) * lQ4
                              + ( 4 * b0.at(nf) * pow(gK0.at(nf), 2) / 3 - gF0.at(nf) * pow(gK0.at(nf), 2) ) * lQ5
                              - pow(gK0.at(nf), 3) / 6 * lQ6 );

      // Return hard factor
      return H;
    };

    return HardFactor;
  }
}
