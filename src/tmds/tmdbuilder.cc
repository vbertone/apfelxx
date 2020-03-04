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
#include "apfel/gammacusp.h"
#include "apfel/gammav.h"
#include "apfel/gammacs.h"
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
    const std::map<int, DglapObjects> DglapObjpdf = InitializeDglapObjectsQCD(g, Thresholds, IntEps);
    const std::map<int, DglapObjects> DglapObjff  = InitializeDglapObjectsQCDT(g, Thresholds, IntEps);

    report("Initializing TMD objects for matching and evolution... ");
    Timer t;

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the thresholds
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
        nfi++;

    // ===============================================================
    // LO matching functions operators.
    std::map<int, Operator> C00;
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};
    C00.insert({EvolutionBasisQCD::PNSP, Id});
    C00.insert({EvolutionBasisQCD::PNSM, Id});
    C00.insert({EvolutionBasisQCD::PNSV, Id});
    C00.insert({EvolutionBasisQCD::PQQ,  Id});
    C00.insert({EvolutionBasisQCD::PQG,  Zero});
    C00.insert({EvolutionBasisQCD::PGQ,  Zero});
    C00.insert({EvolutionBasisQCD::PGG,  Id});

    // ===============================================================
    // NLO matching functions operators.
    // PDFs
    std::map<int, std::map<int, Operator>> C10pdf;
    const Operator O1nspdf{g, C1nspdf{}, IntEps};
    const Operator O1qgpdf{g, C1qgpdf{}, IntEps};
    const Operator O1gqpdf{g, C1gqpdf{}, IntEps};
    const Operator O1ggpdf{g, C1ggpdf{}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O1qgpdfnf = nf * O1qgpdf;
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O1nspdf});
        OM.insert({EvolutionBasisQCD::PNSM, O1nspdf});
        OM.insert({EvolutionBasisQCD::PNSV, O1nspdf});
        OM.insert({EvolutionBasisQCD::PQQ,  O1nspdf});
        OM.insert({EvolutionBasisQCD::PQG,  O1qgpdfnf});
        OM.insert({EvolutionBasisQCD::PGQ,  O1gqpdf});
        OM.insert({EvolutionBasisQCD::PGG,  O1ggpdf});
        C10pdf.insert({nf, OM});
      }

    // Terms proportion to one power of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C11pdf;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O11gmVq = gammaVq0() * Id / 2;
        const Operator O11gmVg = gammaVg0(nf) * Id / 2;
        const auto P0 = DglapObjpdf.at(nf).SplittingFunctions.at(0);
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O11gmVq - P0.at(0)});
        OM.insert({EvolutionBasisQCD::PNSM, O11gmVq - P0.at(1)});
        OM.insert({EvolutionBasisQCD::PNSV, O11gmVq - P0.at(2)});
        OM.insert({EvolutionBasisQCD::PQQ,  O11gmVq - P0.at(3)});
        OM.insert({EvolutionBasisQCD::PQG,     (-1) * P0.at(4)});
        OM.insert({EvolutionBasisQCD::PGQ,     (-1) * P0.at(5)});
        OM.insert({EvolutionBasisQCD::PGG,  O11gmVg - P0.at(6)});
        C11pdf.insert({nf, OM});
      }

    // Terms proportion to two powers of log(mu0/mub)
    std::map<int, Operator> C12;
    C12.insert({EvolutionBasisQCD::PNSP, CF * Id});
    C12.insert({EvolutionBasisQCD::PNSM, CF * Id});
    C12.insert({EvolutionBasisQCD::PNSV, CF * Id});
    C12.insert({EvolutionBasisQCD::PQQ,  CF * Id});
    C12.insert({EvolutionBasisQCD::PQG,  Zero});
    C12.insert({EvolutionBasisQCD::PGQ,  Zero});
    C12.insert({EvolutionBasisQCD::PGG,  CA * Id});

    // FFs
    std::map<int, std::map<int, Operator>> C10ff;
    const Operator O1nsff{g, C1nsff{}, IntEps};
    const Operator O1qgff{g, C1qgff{}, IntEps};
    const Operator O1gqff{g, C1gqff{}, IntEps};
    const Operator O1ggff{g, C1ggff{}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O1qgffnf = nf * O1qgff;
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O1nsff});
        OM.insert({EvolutionBasisQCD::PNSM, O1nsff});
        OM.insert({EvolutionBasisQCD::PNSV, O1nsff});
        OM.insert({EvolutionBasisQCD::PQQ,  O1nsff});
        OM.insert({EvolutionBasisQCD::PQG,  O1qgffnf});
        OM.insert({EvolutionBasisQCD::PGQ,  O1gqff});
        OM.insert({EvolutionBasisQCD::PGG,  O1ggff});
        C10ff.insert({nf, OM});
      }

    // Terms proportion to one power of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C11ff;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O11gmVq = gammaVq0() * Id / 2;
        const Operator O11gmVg = gammaVg0(nf) * Id / 2;
        const auto P0 = DglapObjff.at(nf).SplittingFunctions.at(0);
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O11gmVq - P0.at(0)});
        OM.insert({EvolutionBasisQCD::PNSM, O11gmVq - P0.at(1)});
        OM.insert({EvolutionBasisQCD::PNSV, O11gmVq - P0.at(2)});
        OM.insert({EvolutionBasisQCD::PQQ,  O11gmVq - P0.at(3)});
        OM.insert({EvolutionBasisQCD::PQG,     (-1) * P0.at(4)});
        OM.insert({EvolutionBasisQCD::PGQ,     (-1) * P0.at(5)});
        OM.insert({EvolutionBasisQCD::PGG,  O11gmVg - P0.at(6)});
        C11ff.insert({nf, OM});
      }

    // Terms proportion to two powers of log(mu0/mub) equal to that of
    // PDFs.

    // ===============================================================
    // NNLO matching functions operators.
    // PDFs
    std::map<int, std::map<int, Operator>> C20pdf;
    const Operator O2Vqqbpdf{g, C2Vqqbpdf{}, IntEps};
    const Operator O2pspdf{g, C2pspdf{}, IntEps};
    const Operator O2qgpdf{g, C2qgpdf{}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O2Vqqpdf{g, C2Vqqpdf{nf}, IntEps};
        const Operator O2qgpdfnf = nf * O2qgpdf;
        const Operator O2gqpdf{g, C2gqpdf{nf}, IntEps};
        const Operator O2ggpdf{g, C2ggpdf{nf}, IntEps};
        const Operator O2nsppdf = O2Vqqpdf + O2Vqqbpdf;
        const Operator O2nsmpdf = O2Vqqpdf - O2Vqqbpdf;
        const Operator O2qqpdf  = O2nsppdf + nf * O2pspdf;
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O2nsppdf});
        OM.insert({EvolutionBasisQCD::PNSM, O2nsmpdf});
        OM.insert({EvolutionBasisQCD::PNSV, O2nsmpdf});
        OM.insert({EvolutionBasisQCD::PQQ,  O2qqpdf});
        OM.insert({EvolutionBasisQCD::PQG,  O2qgpdfnf});
        OM.insert({EvolutionBasisQCD::PGQ,  O2gqpdf});
        OM.insert({EvolutionBasisQCD::PGG,  O2ggpdf});
        C20pdf.insert({nf, OM});
      }

    // Terms proportion to one power of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C21pdf;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const double factq = beta0qcd(nf) + gammaVq0() / 2;
        const double factg = beta0qcd(nf) + gammaVg0(nf) / 2;
        const Operator O21gmVq = gammaVq1(nf) * Id / 2;
        const Operator O21gmVg = gammaVg1(nf) * Id / 2;
        const auto P0 = DglapObjpdf.at(nf).SplittingFunctions.at(0);
        const auto P1 = DglapObjpdf.at(nf).SplittingFunctions.at(1);
        const auto C1 = C10pdf.at(nf);
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O21gmVq - P1.at(0) + factq * C1.at(0) - C1.at(0) * P0.at(0)});
        OM.insert({EvolutionBasisQCD::PNSM, O21gmVq - P1.at(1) + factq * C1.at(1) - C1.at(1) * P0.at(1)});
        OM.insert({EvolutionBasisQCD::PNSV, O21gmVq - P1.at(2) + factq * C1.at(2) - C1.at(2) * P0.at(2)});
        OM.insert({EvolutionBasisQCD::PQQ,  O21gmVq - P1.at(3) + factq * C1.at(3) - C1.at(3) * P0.at(3) - C1.at(4) * P0.at(5)});
        OM.insert({EvolutionBasisQCD::PQG,     (-1) * P1.at(4) + factq * C1.at(4) - C1.at(3) * P0.at(4) - C1.at(4) * P0.at(6)});
        OM.insert({EvolutionBasisQCD::PGQ,     (-1) * P1.at(5) + factg * C1.at(5) - C1.at(5) * P0.at(3) - C1.at(6) * P0.at(5)});
        OM.insert({EvolutionBasisQCD::PGG,  O21gmVg - P1.at(6) + factg * C1.at(6) - C1.at(5) * P0.at(4) - C1.at(6) * P0.at(6)});
        C21pdf.insert({nf, OM});
      }

    // Terms proportion to two powers of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C22pdf;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const double b0    = beta0qcd(nf);
        const double gv0q  = gammaVq0();
        const double gv0g  = gammaVg0(nf);
        const double gcp   = GammaCusp1(nf);
        const double factq = - ( gv0q + b0 ) / 2;
        const double factg = - ( gv0g + b0 ) / 2;
        const Operator O22gmVq = ( CF * gcp / 8 + gv0q * ( 2 * b0 + gv0q ) / 8 ) * Id;
        const Operator O22gmVg = ( CA * gcp / 8 + gv0g * ( 2 * b0 + gv0g ) / 8 ) * Id;
        const auto P0 = DglapObjpdf.at(nf).SplittingFunctions.at(0);
        const auto C1 = C10pdf.at(nf);
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O22gmVq + factq * P0.at(0) + CF * C1.at(0) + P0.at(0) * P0.at(0) / 2});
        OM.insert({EvolutionBasisQCD::PNSM, O22gmVq + factq * P0.at(1) + CF * C1.at(1) + P0.at(1) * P0.at(1) / 2});
        OM.insert({EvolutionBasisQCD::PNSV, O22gmVq + factq * P0.at(2) + CF * C1.at(2) + P0.at(2) * P0.at(2) / 2});
        OM.insert({EvolutionBasisQCD::PQQ,  O22gmVq + factq * P0.at(3) + CF * C1.at(3) + ( P0.at(3) * P0.at(3) + P0.at(4) * P0.at(5) ) / 2});
        OM.insert({EvolutionBasisQCD::PQG,            factq * P0.at(4) + CF * C1.at(4) + ( P0.at(3) * P0.at(4) + P0.at(4) * P0.at(6) ) / 2});
        OM.insert({EvolutionBasisQCD::PGQ,            factg * P0.at(5) + CA * C1.at(5) + ( P0.at(5) * P0.at(3) + P0.at(6) * P0.at(5) ) / 2});
        OM.insert({EvolutionBasisQCD::PGG,  O22gmVg + factg * P0.at(6) + CA * C1.at(6) + ( P0.at(5) * P0.at(4) + P0.at(6) * P0.at(6) ) / 2});
        C22pdf.insert({nf, OM});
      }

    // Terms proportion to three powers of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C23pdf;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O23gmVq = CF * ( beta0qcd(nf) / 3 + gammaVq0() / 2 ) * Id;
        const Operator O23gmVg = CA * ( beta0qcd(nf) / 3 + gammaVg0(nf) / 2 ) * Id;;
        const auto P0 = DglapObjpdf.at(nf).SplittingFunctions.at(0);
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O23gmVq - CF * P0.at(0)});
        OM.insert({EvolutionBasisQCD::PNSM, O23gmVq - CF * P0.at(1)});
        OM.insert({EvolutionBasisQCD::PNSV, O23gmVq - CF * P0.at(2)});
        OM.insert({EvolutionBasisQCD::PQQ,  O23gmVq - CF * P0.at(3)});
        OM.insert({EvolutionBasisQCD::PQG,          - CF * P0.at(4)});
        OM.insert({EvolutionBasisQCD::PGQ,          - CA * P0.at(5)});
        OM.insert({EvolutionBasisQCD::PGG,  O23gmVg - CA * P0.at(6)});
        C23pdf.insert({nf, OM});
      }

    // Terms proportion to four powers of log(mu0/mub)
    std::map<int, Operator> C24;
    C24.insert({EvolutionBasisQCD::PNSP, CF * CF * Id / 2});
    C24.insert({EvolutionBasisQCD::PNSM, CF * CF * Id / 2});
    C24.insert({EvolutionBasisQCD::PNSV, CF * CF * Id / 2});
    C24.insert({EvolutionBasisQCD::PQQ,  CF * CF * Id / 2});
    C24.insert({EvolutionBasisQCD::PQG,  Zero});
    C24.insert({EvolutionBasisQCD::PGQ,  Zero});
    C24.insert({EvolutionBasisQCD::PGG,  CA * CA * Id / 2});

    // FFs
    std::map<int, std::map<int, Operator>> C20ff;
    const Operator O2Vqqbff{g, C2Vqqbff{}, IntEps};
    const Operator O2psff{g, C2psff{}, IntEps};
    const Operator O2qgff{g, C2qgff{}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O2Vqqff{g, C2Vqqff{nf}, IntEps};
        const Operator O2qgffnf = nf * O2qgff;
        const Operator O2gqff{g, C2gqff{nf}, IntEps};
        const Operator O2ggff{g, C2ggff{nf}, IntEps};
        const Operator O2nspff = O2Vqqff + O2Vqqbff;
        const Operator O2nsmff = O2Vqqff - O2Vqqbff;
        const Operator O2qqff  = O2nspff + nf * O2psff;
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O2nspff});
        OM.insert({EvolutionBasisQCD::PNSM, O2nsmff});
        OM.insert({EvolutionBasisQCD::PNSV, O2nsmff});
        OM.insert({EvolutionBasisQCD::PQQ,  O2qqff});
        OM.insert({EvolutionBasisQCD::PQG,  O2qgffnf});
        OM.insert({EvolutionBasisQCD::PGQ,  O2gqff});
        OM.insert({EvolutionBasisQCD::PGG,  O2ggff});
        C20ff.insert({nf, OM});
      }

    // Terms proportion to one power of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C21ff;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const double factq = beta0qcd(nf) + gammaVq0() / 2;
        const double factg = beta0qcd(nf) + gammaVg0(nf) / 2;
        const Operator O21gmVq = gammaVq1(nf) * Id / 2;
        const Operator O21gmVg = gammaVg1(nf) * Id / 2;
        const auto P0 = DglapObjff.at(nf).SplittingFunctions.at(0);
        const auto P1 = DglapObjff.at(nf).SplittingFunctions.at(1);
        const auto C1 = C10ff.at(nf);
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O21gmVq - P1.at(0) + factq * C1.at(0) - C1.at(0) * P0.at(0)});
        OM.insert({EvolutionBasisQCD::PNSM, O21gmVq - P1.at(1) + factq * C1.at(1) - C1.at(1) * P0.at(1)});
        OM.insert({EvolutionBasisQCD::PNSV, O21gmVq - P1.at(2) + factq * C1.at(2) - C1.at(2) * P0.at(2)});
        OM.insert({EvolutionBasisQCD::PQQ,  O21gmVq - P1.at(3) + factq * C1.at(3) - C1.at(3) * P0.at(3) - C1.at(4) * P0.at(5)});
        OM.insert({EvolutionBasisQCD::PQG,     (-1) * P1.at(4) + factq * C1.at(4) - C1.at(3) * P0.at(4) - C1.at(4) * P0.at(6)});
        OM.insert({EvolutionBasisQCD::PGQ,     (-1) * P1.at(5) + factg * C1.at(5) - C1.at(5) * P0.at(3) - C1.at(6) * P0.at(5)});
        OM.insert({EvolutionBasisQCD::PGG,  O21gmVg - P1.at(6) + factg * C1.at(6) - C1.at(5) * P0.at(4) - C1.at(6) * P0.at(6)});
        C21ff.insert({nf, OM});
      }

    // Terms proportion to two powers of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C22ff;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const double b0    = beta0qcd(nf);
        const double gv0q  = gammaVq0();
        const double gv0g  = gammaVg0(nf);
        const double gcp   = GammaCusp1(nf);
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
        OM.insert({EvolutionBasisQCD::PQQ,  O22gmVq + factq * P0.at(3) + CF * C1.at(3) + ( P0.at(3) * P0.at(3) + P0.at(4) * P0.at(5) ) / 2});
        OM.insert({EvolutionBasisQCD::PQG,            factq * P0.at(4) + CF * C1.at(4) + ( P0.at(3) * P0.at(4) + P0.at(4) * P0.at(6) ) / 2});
        OM.insert({EvolutionBasisQCD::PGQ,            factg * P0.at(5) + CA * C1.at(5) + ( P0.at(5) * P0.at(3) + P0.at(6) * P0.at(5) ) / 2});
        OM.insert({EvolutionBasisQCD::PGG,  O22gmVg + factg * P0.at(6) + CA * C1.at(6) + ( P0.at(5) * P0.at(4) + P0.at(6) * P0.at(6) ) / 2});
        C22ff.insert({nf, OM});
      }

    // Terms proportion to three powers of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C23ff;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O23gmVq = CF * ( beta0qcd(nf) / 3 + gammaVq0() / 2 ) * Id;
        const Operator O23gmVg = CA * ( beta0qcd(nf) / 3 + gammaVg0(nf) / 2 ) * Id;;
        const auto P0 = DglapObjff.at(nf).SplittingFunctions.at(0);
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O23gmVq - CF * P0.at(0)});
        OM.insert({EvolutionBasisQCD::PNSM, O23gmVq - CF * P0.at(1)});
        OM.insert({EvolutionBasisQCD::PNSV, O23gmVq - CF * P0.at(2)});
        OM.insert({EvolutionBasisQCD::PQQ,  O23gmVq - CF * P0.at(3)});
        OM.insert({EvolutionBasisQCD::PQG,          - CF * P0.at(4)});
        OM.insert({EvolutionBasisQCD::PGQ,          - CA * P0.at(5)});
        OM.insert({EvolutionBasisQCD::PGG,  O23gmVg - CA * P0.at(6)});
        C23ff.insert({nf, OM});
      }

    // Terms proportion to four powers of log(mu0/mub) equal to that
    // of PDFs.

    // Define map containing the TmdObjects for each nf.
    std::map<int, TmdObjects> TmdObj;

    // Construct sets of operators for each perturbative order for the
    // matching functions. Initialize also coefficients of: beta
    // function, GammaCusp, gammaV, and Collins-Soper anomalous
    // dimensions.
    std::map<int, std::map<int, Set<Operator>>> MatchingFunctionsPDFs;
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
        obj.GammaFq.insert({0, gammaVq0()});
        obj.GammaFq.insert({1, gammaVq1(nf)});
        obj.GammaFq.insert({2, gammaVq2(nf)});

        // GammaF gluon
        obj.GammaFg.insert({0, gammaVg0(nf)});
        obj.GammaFg.insert({1, gammaVg1(nf)});
        obj.GammaFg.insert({2, gammaVg2(nf)});

        // GammaCusp (multiply by CF for quarks and by CA for gluons)
        obj.GammaCusp.insert({0, GammaCusp0()});
        obj.GammaCusp.insert({1, GammaCusp1(nf)});
        obj.GammaCusp.insert({2, GammaCusp2(nf)});
        obj.GammaCusp.insert({3, GammaCusp3(nf)});

        // Collins-Soper anomalous dimensions (multiply by CF for
        // quarks and by CA for gluons).
        obj.GammaCS.insert({0, {CSd10(), CSd11()}});
        obj.GammaCS.insert({1, {CSd20(nf), CSd21(nf), CSd22(nf)}});
        obj.GammaCS.insert({2, {CSd30(nf), CSd31(nf), CSd32(nf), CSd33(nf)}});

        // Matching functions.
        const EvolutionBasisQCD evb{nf};

        // PDFs
        obj.MatchingFunctionsPDFs.insert({0, {{evb, C00}}});
        obj.MatchingFunctionsPDFs.insert({1, {{evb, C10pdf.at(nf)}, {evb, C11pdf.at(nf)}, {evb, C12}}});
        obj.MatchingFunctionsPDFs.insert({2, {{evb, C20pdf.at(nf)}, {evb, C21pdf.at(nf)}, {evb, C22pdf.at(nf)}, {evb, C23pdf.at(nf)}, {evb, C24}}});

        // FFs
        obj.MatchingFunctionsFFs.insert({0, {{evb, C00}}});
        obj.MatchingFunctionsFFs.insert({1, {{evb, C10ff.at(nf)}, {evb, C11ff.at(nf)}, {evb, C12}}});
        obj.MatchingFunctionsFFs.insert({2, {{evb, C20ff.at(nf)}, {evb, C21ff.at(nf)}, {evb, C22ff.at(nf)}, {evb, C23ff.at(nf)}, {evb, C24}}});

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
    // Match TMDs on collinear PDFs.
    const std::function<Set<Distribution>(double const&)> MatchedTmdPDFs = MatchTmdPDFs(TmdObj, CollPDFs, Alphas, PerturbativeOrder, Ci);

    // Compute TMD evolution factors.
    const std::function<std::vector<double>(double const&, double const&, double const&)> EvolFactors = EvolutionFactors(TmdObj, Alphas, PerturbativeOrder, Ci, IntEps);

    // Computed TMDPDFs at the final scale by multiplying initial
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
    // Match TMDs on collinear FFs.
    const std::function<Set<Distribution>(double const&)> MatchedTmdFFs = MatchTmdFFs(TmdObj, CollFFs, Alphas, PerturbativeOrder, Ci);

    // Compute TMD evolution factors.
    const std::function<std::vector<double>(double const&, double const&, double const&)> EvolFactors = EvolutionFactors(TmdObj, Alphas, PerturbativeOrder, Ci, IntEps);

    // Computed TMDFFs at the final scale by multiplying initial scale
    // TMDFFs by the evolution factor.
    const auto EvolvedTMDs = [=] (double const& b, double const& muf, double const& zetaf) -> Set<Distribution>
    {
      return EvolFactors(b, muf, zetaf) * MatchedTmdFFs(b);
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
    // Retrieve thresholds from "TmdObj".
    std::vector<double> thrs;
    for(auto const& obj : TmdObj)
      {
        const int    nf  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) thrs.size() < nf)
          thrs.resize(nf);
        thrs[nf-1] = thr;
      }

    // Define the log(Ci) to assess scale variations.
    const double Lmu = 2 * log(Ci);

    // Matching functions as functions of the absolute value of the
    // impact parameter b.
    std::function<Set<Operator>(double const&)> MatchFunc;
    if (PerturbativeOrder == LL || PerturbativeOrder == NLL)
      MatchFunc = [=] (double const& mu) -> Set<Operator>
      {
        return TmdObj.at(NF(mu,thrs)).MatchingFunctionsPDFs.at(0)[0];
      };
    else if (PerturbativeOrder == NNLL || PerturbativeOrder == NLLp)
      MatchFunc = [=] (double const& mu) -> Set<Operator>
      {
        const double coup = Alphas(mu) / FourPi;
        const auto& mf = TmdObj.at(NF(mu,thrs)).MatchingFunctionsPDFs;
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
        const auto& mf  = TmdObj.at(NF(mu,thrs)).MatchingFunctionsPDFs;
        const auto c0   = mf.at(0);
        const auto c1   = mf.at(1);
        const auto c2   = mf.at(2);
        const auto lo   = c0[0];
        const auto nlo  = c1[0] + Lmu * ( c1[1] + Lmu * c1[2] );
        const auto nnlo = c2[0] + Lmu * ( c2[1] + Lmu * ( c2[2] + Lmu * ( c2[3] + Lmu * c2[4] ) ) );
        return lo + coup * ( nlo + coup * nnlo );
      };

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
    // Retrieve thresholds from "TmdObj".
    std::vector<double> thrs;
    for(auto const& obj : TmdObj)
      {
        const int    nf  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) thrs.size() < nf)
          thrs.resize(nf);
        thrs[nf-1] = thr;
      }

    // Define the log(Ci) to assess scale variations.
    const double Lmu = 2 * log(Ci);

    // Matching functions as functions of the absolute value of the
    // impact parameter b.
    std::function<Set<Operator>(double const&)> MatchFunc;
    if (PerturbativeOrder == LL || PerturbativeOrder == NLL)
      MatchFunc = [=] (double const& mu) -> Set<Operator>
      {
        return TmdObj.at(NF(mu,thrs)).MatchingFunctionsFFs.at(0)[0];
      };
    else if (PerturbativeOrder == NNLL || PerturbativeOrder == NLLp)
      MatchFunc = [=] (double const& mu) -> Set<Operator>
      {
        const double coup = Alphas(mu) / FourPi;
        const auto& mf = TmdObj.at(NF(mu,thrs)).MatchingFunctionsFFs;
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
        const auto& mf  = TmdObj.at(NF(mu,thrs)).MatchingFunctionsFFs;
        const auto c0   = mf.at(0);
        const auto c1   = mf.at(1);
        const auto c2   = mf.at(2);
        const auto lo   = c0[0];
        const auto nlo  = c1[0] + Lmu * ( c1[1] + Lmu * c1[2] );
        const auto nnlo = c2[0] + Lmu * ( c2[1] + Lmu * ( c2[2] + Lmu * ( c2[3] + Lmu * c2[4] ) ) );
        return lo + coup * ( nlo + coup * nnlo );
      };

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
  std::function<std::vector<double>(double const&, double const&, double const&)> EvolutionFactors(std::map<int, TmdObjects>            const& TmdObj,
                                                                                                   std::function<double(double const&)> const& Alphas,
                                                                                                   int                                  const& PerturbativeOrder,
                                                                                                   double                               const& Ci,
                                                                                                   double                               const& IntEps)
  {
    // Retrieve thresholds from "TmdObj".
    std::vector<double> thrs;
    for(auto const& obj : TmdObj)
      {
        const int    nf  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) thrs.size() < nf)
          thrs.resize(nf);
        thrs[nf-1] = thr;
      }

    // Define the log(Ci) to assess scale variations.
    const double Lmu = 2 * log(Ci);

    // Create functions needed for the TMD evolution.
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
          return coup * TmdObj.at(NF(mu,thrs)).GammaCusp.at(0);
        };
        K = [=] (double const&) -> double{ return 0; };
      }
    // NLL
    else if (PerturbativeOrder == NLL || PerturbativeOrder == NLLp)
      {
        gammaFq = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * TmdObj.at(NF(mu,thrs)).GammaFq.at(0);
        };
        gammaFg = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * TmdObj.at(NF(mu,thrs)).GammaFg.at(0);
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu,thrs)).GammaCusp;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * gc.at(1) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu,thrs)).GammaCS;
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
          const auto& gv    = TmdObj.at(NF(mu,thrs)).GammaFq;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * gv.at(1) );
        };
        gammaFg = [=] (double const& mu) -> double
        {
          const auto& gv    = TmdObj.at(NF(mu,thrs)).GammaFg;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * gv.at(1) );
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu,thrs)).GammaCusp;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * gc.at(2) ) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu,thrs)).GammaCS;
          const std::vector<double> d0 = d.at(0);
          const std::vector<double> d1 = d.at(1);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          const double nlo  = d1[0] + Lmu * ( d1[1] + Lmu * d1[2] );
          return coup * ( lo + coup * nlo );
        };
      }
    // N3LL
    else if (PerturbativeOrder == NNNLL)
      {
        gammaFq = [=] (double const& mu) -> double
        {
          const auto& gv    = TmdObj.at(NF(mu,thrs)).GammaFq;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * ( gv.at(1) + coup * gv.at(2) ) );
        };
        gammaFg = [=] (double const& mu) -> double
        {
          const auto& gv    = TmdObj.at(NF(mu,thrs)).GammaFg;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * ( gv.at(1) + coup * gv.at(2) ) );
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu,thrs)).GammaCusp;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * ( gc.at(2) + coup * gc.at(3) ) ) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu,thrs)).GammaCS;
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

    // Define the integrands.
    Integrator I1q{[=] (double const& mu) -> double{ return gammaFq(mu) / mu; }};
    Integrator I1g{[=] (double const& mu) -> double{ return gammaFg(mu) / mu; }};
    Integrator I3 {[=] (double const& mu) -> double{ return gammaK(mu) * log(mu) / mu; }};
    Integrator I2 {[=] (double const& mu) -> double{ return gammaK(mu) / mu; }};

    // Construct function that returns the perturbative evolution
    // kernel.
    const auto EvolFactors = [=] (double const& b, double const& muf, double const& zetaf) -> std::vector<double>
    {
      // Define lower scales
      const double mu0   = Ci * 2 * exp(- emc) / b;
      const double zeta0 = mu0 * mu0;

      // Compute argument of the exponent of the evolution factors.
      const double IntI1q = I1q.integrate(mu0, muf, thrs, IntEps);
      const double IntI1g = I1g.integrate(mu0, muf, thrs, IntEps);
      const double IntI2  = I2.integrate(mu0, muf, thrs, IntEps) * log(zetaf);
      const double IntI3  = I3.integrate(mu0, muf, thrs, IntEps);

      // Compute the evolution factors.
      const double Klz = ( K(mu0) * log( zetaf / zeta0 ) - IntI2 ) / 2 + IntI3;
      const double Rq  = exp( CF * Klz + IntI1q );
      const double Rg  = exp( CA * Klz + IntI1g );

      // Return vector of evolution factors.
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
    // Retrieve thresholds from "TmdObj".
    std::vector<double> thrs;
    for(auto const& obj : TmdObj)
      {
        const int    nf  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) thrs.size() < nf)
          thrs.resize(nf);
        thrs[nf-1] = thr;
      }

    // Define the log(Ci) to assess scale variations.
    const double Lmu = 2 * log(Ci);

    // Create functions needed for the TMD evolution.
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
          return coup * TmdObj.at(NF(mu,thrs)).GammaCusp.at(0);
        };
        K = [=] (double const&) -> double{ return 0; };
      }
    // NLL
    else if (PerturbativeOrder == NLL || PerturbativeOrder == NLLp)
      {
        gammaFq = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * TmdObj.at(NF(mu,thrs)).GammaFq.at(0);
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu,thrs)).GammaCusp;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * gc.at(1) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu,thrs)).GammaCS;
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
          const auto& gv    = TmdObj.at(NF(mu,thrs)).GammaFq;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * gv.at(1) );
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu,thrs)).GammaCusp;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * gc.at(2) ) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu,thrs)).GammaCS;
          const std::vector<double> d0 = d.at(0);
          const std::vector<double> d1 = d.at(1);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          const double nlo  = d1[0] + Lmu * ( d1[1] + Lmu * d1[2] );
          return coup * ( lo + coup * nlo );
        };
      }
    // N3LL
    else if (PerturbativeOrder == NNNLL)
      {
        gammaFq = [=] (double const& mu) -> double
        {
          const auto& gv    = TmdObj.at(NF(mu,thrs)).GammaFq;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * ( gv.at(1) + coup * gv.at(2) ) );
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu,thrs)).GammaCusp;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * ( gc.at(2) + coup * gc.at(3) ) ) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu,thrs)).GammaCS;
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

    // Define the integrands.
    Integrator I1{[=] (double const& mu) -> double{ return gammaFq(mu) / mu; }};
    Integrator I3{[=] (double const& mu) -> double{ return gammaK(mu) * log(mu) / mu; }};
    Integrator I2{[=] (double const& mu) -> double{ return gammaK(mu) / mu; }};

    // Construct function that returns the perturbative evolution
    // kernel.
    const auto EvolFactor = [=] (double const& b, double const& muf, double const& zetaf) -> double
    {
      // Define lower scales
      const double mu0   = Ci * 2 * exp(- emc) / b;
      const double zeta0 = mu0 * mu0;

      // Compute argument of the exponent of the evolution factors.
      const double IntI1 = I1.integrate(mu0, muf, thrs, IntEps);
      const double IntI2 = I2.integrate(mu0, muf, thrs, IntEps) * log(zetaf);
      const double IntI3 = I3.integrate(mu0, muf, thrs, IntEps);

      // Compute the evolution factors.
      const double Klz = ( K(mu0) * log( zetaf / zeta0 ) - IntI2 ) / 2 + IntI3;
      const double Rq  = exp( CF * Klz + IntI1 );

      // Return vector of evolution factors.
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
    // Retrieve thresholds from "TmdObj".
    std::vector<double> thrs;
    for(auto const& obj : TmdObj)
      {
        const int    nf  = obj.first;
        const double thr = obj.second.Threshold;
        if ((int) thrs.size() < nf)
          thrs.resize(nf);
        thrs[nf-1] = thr;
      }

    // Define the log(Ci) to assess scale variations.
    const double Lmu = 2 * log(Ci);

    // Create functions needed for the TMD evolution.
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
          return coup * TmdObj.at(NF(mu,thrs)).GammaCusp.at(0);
        };
        K = [=] (double const&) -> double{ return 0; };
      }
    // NLL
    else if (PerturbativeOrder == NLL || PerturbativeOrder == NLLp)
      {
        gammaFg = [=] (double const& mu) -> double
        {
          const double coup = Alphas(mu) / FourPi;
          return coup * TmdObj.at(NF(mu,thrs)).GammaFg.at(0);
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu,thrs)).GammaCusp;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * gc.at(1) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu,thrs)).GammaCS;
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
          const auto& gv    = TmdObj.at(NF(mu,thrs)).GammaFg;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * gv.at(1) );
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu,thrs)).GammaCusp;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * gc.at(2) ) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu,thrs)).GammaCS;
          const std::vector<double> d0 = d.at(0);
          const std::vector<double> d1 = d.at(1);
          const double coup = Alphas(mu) / FourPi;
          const double lo   = d0[0] + Lmu * d0[1];
          const double nlo  = d1[0] + Lmu * ( d1[1] + Lmu * d1[2] );
          return coup * ( lo + coup * nlo );
        };
      }
    // N3LL
    else if (PerturbativeOrder == NNNLL)
      {
        gammaFg = [=] (double const& mu) -> double
        {
          const auto& gv    = TmdObj.at(NF(mu,thrs)).GammaFg;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gv.at(0) + coup * ( gv.at(1) + coup * gv.at(2) ) );
        };
        gammaK = [=] (double const& mu) -> double
        {
          const auto& gc    = TmdObj.at(NF(mu,thrs)).GammaCusp;
          const double coup = Alphas(mu) / FourPi;
          return coup * ( gc.at(0) + coup * ( gc.at(1) + coup * ( gc.at(2) + coup * gc.at(3) ) ) );
        };
        K = [=] (double const& mu) -> double
        {
          const auto& d = TmdObj.at(NF(mu,thrs)).GammaCS;
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

    // Define the integrands.
    Integrator I1{[=] (double const& mu) -> double{ return gammaFg(mu) / mu; }};
    Integrator I3{[=] (double const& mu) -> double{ return gammaK(mu) * log(mu) / mu; }};
    Integrator I2{[=] (double const& mu) -> double{ return gammaK(mu) / mu; }};

    // Construct function that returns the perturbative evolution
    // kernel.
    const auto EvolFactor = [=] (double const& b, double const& muf, double const& zetaf) -> double
    {
      // Define lower scales
      const double mu0   = Ci * 2 * exp(- emc) / b;
      const double zeta0 = mu0 * mu0;

      // Compute argument of the exponent of the evolution factors.
      const double IntI1 = I1.integrate(mu0, muf, thrs, IntEps);
      const double IntI2 = I2.integrate(mu0, muf, thrs, IntEps) * log(zetaf);
      const double IntI3 = I3.integrate(mu0, muf, thrs, IntEps);

      // Compute the evolution factors.
      const double Klz = ( K(mu0) * log( zetaf / zeta0 ) - IntI2 ) / 2 + IntI3;
      const double Rg  = exp( CA * Klz + IntI1 );

      // Return vector of evolution factors.
      return Rg;
    };

    return EvolFactor;
  }

  //_____________________________________________________________________________
  double HardFactorDY(int const& PerturbativeOrder, double const& Alphas, int const& nf, double const& kappa)
  {
    // Compute log and its powers.
    const double lQ2  = 2 * log(kappa);
    const double lQ22 = lQ2 * lQ2;
    const double lQ23 = lQ22 * lQ2;
    const double lQ24 = lQ23 * lQ2;

    // Compute coupling and its powers.
    const double as  = Alphas / FourPi;
    const double as2 = as * as;

    // Now compute hard factor according to the perturbative order.
    double hfct = 1;
    if (PerturbativeOrder > 1 || PerturbativeOrder < 0)
      hfct += 2 * as * CF * ( - lQ22 - 3 * lQ2 - 8 + 7 * Pi2 / 6 );
    if (PerturbativeOrder > 2 || PerturbativeOrder < -1)
      hfct += 2 * as2 * CF *
              ( CF * ( lQ24 + 6 * lQ23 + ( 25. - 7 * Pi2 / 3 ) * lQ22 + ( 93. / 2. - 5 * Pi2 - 24 * zeta3 ) * lQ2
                       + 511. / 8. - 83 * Pi2 / 6 - 30 * zeta3 + 67 * Pi2 * Pi2 / 60 ) +
                CA * ( - 11 * lQ23 / 9 + ( - 233. / 18. + Pi2 / 3 ) * lQ22 + ( - 2545. / 54. + 22 * Pi2 / 9 + 26 * zeta3 ) * lQ2
                       - 51157. / 648. + 1061 * Pi2 / 108 + 313 * zeta3 / 9 - 4 * Pi2 * Pi2 / 45 ) +
                TR * nf * ( 4 * lQ23 / 9 + 38 * lQ22 / 9 + ( 418. / 27. - 8 * Pi2 / 9 ) * lQ2
                            + 4085. / 162. - 91 * Pi2 / 27 + 4 * zeta3 / 9 ) );

    return hfct;
  }

  //_____________________________________________________________________________
  double HardFactorSIDIS(int const& PerturbativeOrder, double const& Alphas, int const& nf, double const& kappa)
  {
    // Compute log and its powers.
    const double lQ2  = 2 * log(kappa);
    const double lQ22 = lQ2 * lQ2;
    const double lQ23 = lQ22 * lQ2;
    const double lQ24 = lQ23 * lQ2;

    // Compute coupling and its powers.
    const double as  = Alphas / FourPi;
    const double as2 = as * as;

    // Now compute hard factor according to the perturbative order.
    double hfct = 1;
    if (PerturbativeOrder > 1 || PerturbativeOrder < 0)
      hfct += 2 * as * CF * ( - lQ22 - 3 * lQ2 - 8 + zeta2 );
    if (PerturbativeOrder > 2 || PerturbativeOrder < -1)
      hfct += 2 * as2 * CF *
              ( CF * ( lQ24 + 6 * lQ23 + ( 25. - 2 * zeta2 ) * lQ22 + ( 93. / 2. + 6 * zeta2 - 24 * zeta3 ) * lQ2
                       + 511. / 8. + 13 * zeta2 - 30 * zeta3 + 39 * zeta4 / 2 ) +
                CA * ( - 11 * lQ23 / 9 + ( - 233. / 18. + 2 * zeta2 ) * lQ22 + ( - 2545. / 54. - 22 * zeta2 / 3 + 26 * zeta3 ) * lQ2
                       - 51157. / 648. - 337 * zeta2 / 18 + 313 * zeta3 / 9 + 22 * zeta4 ) +
                TR * nf * ( 4 * lQ23 / 9 + 38 * lQ22 / 9 + ( 418. / 27. + 8 * zeta2 / 3 ) * lQ2
                            + 4085. / 162. + 46 * zeta2 / 9 + 4 * zeta3 / 9 ) );

    return hfct;
  }

  //_____________________________________________________________________________
  double HardFactorggHiggs(int const& PerturbativeOrder, double const& Alphas, int const& nf, double const& kappa)
  {
    // Scale variations are not implemented yet. The relevant formulas
    // can be found in the appendix of
    // https://arxiv.org/pdf/1805.05916.pdf.
    if (kappa != 1)
      throw std::runtime_error(error("HardFactorggHiggs", "Scale variations not implemented yet."));
    // Compute log and its powers.
    //const double lQ2  = 2 * log(kappa);
    //const double lQ22 = lQ2 * lQ2;
    //const double lQ23 = lQ22 * lQ2;
    //const double lQ24 = lQ23 * lQ2;

    // Compute coupling and its powers.
    const double as  = Alphas / FourPi;
    const double as2 = as * as;

    // PDG values of the Higgs mass
    // (http://pdg.lbl.gov/2019/listings/rpp2019-list-higgs-boson.pdf)
    // and of the top pole mass
    // (http://pdg.lbl.gov/2019/listings/rpp2019-list-t-quark.pdf) in
    // GeV as of July 2019
    const double mH = 125.1;
    const double mt = 172.9;
    const double Lt = 2 * log( mH / mt );

    // Now compute hard factor according to the perturbative order.
    double hfct = 1;
    if (PerturbativeOrder > 1 || PerturbativeOrder < 0)
      hfct += 2 * as * ( CA * ( 5 + 7 * Pi2 / 6 ) - 3 * CF );
    if (PerturbativeOrder > 2 || PerturbativeOrder < -1)
      hfct += as2 * ( - 135 * CA + 23827 * CA * CA - 216 * CF - 15660 * CA * CF + 5832 * CF * CF
                      + 2268 * CA * CA * Lt - 3564 * CA * CF * Lt - 4510 * CA * nf - 4428 * CF * nf
                      + 1296 * CF * nf * Lt + 6795 * CA * CA * Pi2 - 2268 * CA * CF * Pi2
                      - 450 * Pi2 * CA * nf + 333 * Pi2 * Pi2 * CA * CA - 5148 * CA * CA * zeta3
                      - 1656 * zeta3 * CA * nf + 2592 * zeta3 * CF * nf ) / 162;

    return hfct;
  }
  /*
    //_____________________________________________________________________________
    void test(std::map<int,TmdObjects>             const& TmdObj,
              std::function<double(double const&)> const& Alphas,
              int                                  const& PerturbativeOrder,
              double                               const& Ci)
    {
      const double Q   = 1;
      const double muR = 1;
      const double xQ  = 1;

      const double Q2   = Q * Q;
      const double muR2 = muR * muR;
      const double xQ2  = xQ * xQ;

      //const auto& gc    = TmdObj.at(NF(muR, thrs)).GammaCusp;
      //const double coup = Alphas(muR) / FourPi;

      // Coefficients A's in terms of GammaCusp, GammaCS (K) and BetaQCD
      // (see Eq. (74) of https://arxiv.org/pdf/1007.4005.pdf)
      const double A1 = 0;
      const double A2 = 0;
      const double A3 = 0;
      const double A4 = 0;

      // Coefficients B's in terms of GammaF, GammaCS (K) and BetaQCD
      // (see Eq. (74) of https://arxiv.org/pdf/1007.4005.pdf)
      const double B1 = 0;
      const double B2 = 0;
      const double B3 = 0;

      //lambda = ( as(muR) / FourPi ) * beta0qcd(nf) * L;

      const auto g1 = [=] (int const& nf, double const& lambda) -> double
      {
        const double twl    = 2 * lambda;
        const double omtwl  = 1 - twl;
        const double lomtwl = log(omtwl);
        return 4 * A1 / beta0qcd(nf) * ( twl + lomtwl ) / twl;
      };

      const auto g2 = [=] (int const& nf, double const& lambda) -> double
      {
        const double twl    = 2 * lambda;
        const double omtwl  = 1 - twl;
        const double lomtwl = log( omtwl);
        return
        2 / beta0qcd(nf) * lomtwl * ( A1 * log( 1 / xQ2 ) + B1 )
        - A2 / pow(beta0qcd(nf) / 2, 2) * ( twl + omtwl * lomtwl ) / omtwl
        + A1 * ( - ( ( beta1qcd(nf) / pow(FourPi, 2) ) ) / ( 4 * M_PI * pow(beta0qcd(nf) / FourPi, 3) )
                 * ( lomtwl * ( ( twl - 1 ) * lomtwl - 2 ) - 4 * lambda ) / omtwl
                 - 1 / ( beta0qcd(nf) / 2 ) * ( ( twl * ( 1 - lomtwl ) + lomtwl ) ) / omtwl * log( muR2 / xQ2 / Q2 ) );
      };

      const auto g3 = [=] (int const& nf,  double const& lambda) -> double
      {
        const double twl    = 2 * lambda;
        const double twl2   = twl * twl;
        const double omtwl  = 1 - twl;
        const double lomtwl = log(omtwl);
        return
        ( A1 * log( 1 / xQ2 ) + B1 ) * ( - lambda / omtwl * log( muR2 / xQ2 / Q2 ) + ( ( beta1qcd(nf) / pow(FourPi, 2) ) )
                                         / ( 2 * pow(beta0qcd(nf) / FourPi, 2) ) * ( twl + lomtwl ) / omtwl )
        - 1 / ( beta0qcd(nf) / 2 ) * lambda / omtwl * ( A2 * log( 1 / xQ2 ) + B2 )
        - A3 / ( pow(beta0qcd(nf) / 2, 2) ) * twl2 / pow(omtwl, 2)
        + A2 * ( ( ( beta1qcd(nf) / pow(FourPi, 2) ) ) / ( 4 * M_PI * pow(beta0qcd(nf) / FourPi, 3) )
                 * ( twl * ( 3 * lambda - 1 ) + ( 4 * lambda - 1 ) * lomtwl ) / pow(omtwl, 2)
                 - 1 / ( beta0qcd(nf) / 4 ) * twl2 / pow(omtwl, 2) * log( muR2 / xQ2 / Q2 ) )
        + A1 * ( ( lambda * ( beta0qcd(nf) / FourPi * ( beta2qcd(nf) / pow(FourPi, 3) ) * ( 1 - 3 * lambda )
                              + pow(beta1qcd(nf) / pow(FourPi, 2), 2) * lambda ) )
                 / ( pow(beta0qcd(nf) / FourPi, 4) * pow(omtwl, 2) )
                 + ( omtwl * lomtwl * ( beta0qcd(nf) / FourPi * ( beta2qcd(nf) / pow(FourPi, 3) ) * omtwl
                                        + 2 * pow(beta1qcd(nf) / pow(FourPi, 2), 2) * lambda ) )
                 / ( 2 * pow(beta0qcd(nf) / FourPi, 4) * pow(omtwl, 2) )
                 + ( pow(beta1qcd(nf) / pow(FourPi, 2), 2) ) / ( 4 * pow(beta0qcd(nf) / FourPi, 4) )
                 * ( ( 1 - 4 * lambda ) * pow(lomtwl, 2) ) / pow(omtwl, 2)
                 - twl2 / pow(omtwl, 2) * pow(log( muR2 / xQ2 / Q2 ), 2)
                 - ( ( beta1qcd(nf) / pow(FourPi, 2) ) ) / ( 2 * pow(beta0qcd(nf) / FourPi, 2) )
                 * ( ( twl * omtwl + ( 1 - 4 * lambda ) * lomtwl ) ) / pow(omtwl, 2) * log( muR2 / xQ2 / Q2 ) );
      };

      const auto g4 = [=] (int const& nf,  double const& lambda) -> double
      {
        const double twl    = 2 * lambda;
        const double twl2   = twl * twl;
        const double twl3   = twl2 * twl;
        const double omtwl  = 1 - twl;
        const double lomtwl = log(omtwl);
        return
        ( A4 * ( 3 - twl ) * twl2 ) / ( 6 * pow(beta0qcd(nf) / 2, 2) * pow(twl - 1, 3) )
        + A3 / ( 48 * M_PI * pow(beta0qcd(nf) / FourPi, 3) * pow(twl - 1, 3) )
        * ( 3 * ( beta1qcd(nf) / pow(FourPi, 2) ) * ( 1 - 6 * lambda ) * lomtwl
            + twl * ( ( beta1qcd(nf) / pow(FourPi, 2) ) * ( 5 * lambda * ( twl - 3 ) + 3 )
                      + 6 * pow(beta0qcd(nf) / FourPi, 2) * ( 3 - twl ) * lambda * log( muR2 / xQ2 / Q2 ) )
            + 12 * pow(beta0qcd(nf) / FourPi, 2) * ( lambda - 1 ) * lambda * ( twl - 1 ) * log( 1 / xQ2 ) )
        + A2 / ( 24 * pow(beta0qcd(nf) / FourPi, 4) * pow(twl - 1, 3) )
        * ( 32 * beta0qcd(nf) / FourPi * ( beta2qcd(nf) / pow(FourPi, 3) ) * twl3
            - 2 * pow(beta1qcd(nf) / pow(FourPi, 2), 2) * lambda * ( lambda * ( 11 * twl - 9 ) + 3 )
            + 12 * pow(beta0qcd(nf) / FourPi, 4) * ( 3 - twl ) * twl2 * pow(log( muR2 / xQ2 / Q2 ), 2)
            + 6 * pow(beta0qcd(nf) / FourPi, 2) * log( muR2 / xQ2 / Q2 )
            * ( ( beta1qcd(nf) / pow(FourPi, 2) ) * ( 1 - 6 * lambda ) * lomtwl
                + 2 * ( lambda - 1 ) * lambda * ( twl - 1 ) * ( ( beta1qcd(nf) / pow(FourPi, 2) )
                                                                + 2 * pow(beta0qcd(nf) / FourPi, 2) * log( 1 / xQ2 ) ) )
            + 3 * ( beta1qcd(nf) / pow(FourPi, 2) ) * ( ( beta1qcd(nf) / pow(FourPi, 2) ) * lomtwl
                                                        * ( twl + ( 6 * lambda - 1 ) * lomtwl - 1 )
                                                        - 2 * pow(beta0qcd(nf) / FourPi, 2) * ( twl - 1 )
                                                        * ( 2 * ( lambda - 1 ) * lambda - lomtwl ) * log( 1 / xQ2 ) ) )
        + ( M_PI * A1 ) / ( 12 * pow(beta0qcd(nf) / FourPi, 5) * pow(twl - 1, 3) )
        * ( pow(beta1qcd(nf) / pow(FourPi, 2), 3) * ( 1 - 6 * lambda ) * pow(lomtwl, 3)
            + 3 * lomtwl * ( pow(beta0qcd(nf) / FourPi, 2) * ( beta3qcd(nf) / pow(FourPi, 4) ) * pow(twl - 1, 3)
                             + beta0qcd(nf) / FourPi * ( beta1qcd(nf) / pow(FourPi, 2) ) * ( beta2qcd(nf) / pow(FourPi, 3) )
                             * ( omtwl * ( 8 * twl2 - 4 * lambda + 3 ) )
                             + 4 * pow(beta1qcd(nf) / pow(FourPi, 2), 3) * twl2 * ( twl + 1 )
                             + pow(beta0qcd(nf) / FourPi, 2) * ( beta1qcd(nf) / pow(FourPi, 2) ) * log( muR2 / xQ2 / Q2 )
                             * ( pow(beta0qcd(nf) / FourPi, 2) * ( 1 - 6 * lambda ) * log( muR2 / xQ2 / Q2 )
                                 - 4 * ( beta1qcd(nf) / pow(FourPi, 2) ) * lambda ) )
            + 3 * pow(beta1qcd(nf) / pow(FourPi, 2), 2) * pow(lomtwl, 2)
            * ( 2 * ( beta1qcd(nf) / pow(FourPi, 2) ) * lambda + pow(beta0qcd(nf) / FourPi, 2)
                * ( 6 * lambda - 1 ) * log( muR2 / xQ2 / Q2 ) )
            + 3 * pow(beta0qcd(nf) / FourPi, 2) * ( twl - 1 ) * log( 1 / xQ2 )
            * ( - pow(beta1qcd(nf) / pow(FourPi, 2), 2) * pow(lomtwl, 2) + 2 * pow(beta0qcd(nf) / FourPi, 2)
                * ( beta1qcd(nf) / pow(FourPi, 2) ) * lomtwl * log( muR2 / xQ2 / Q2 )
                + 4 * lambda * ( lambda * ( pow(beta1qcd(nf) / pow(FourPi, 2), 2)
                                            - beta0qcd(nf) / FourPi * ( beta2qcd(nf) / pow(FourPi, 3) ) )
                                 + pow(beta0qcd(nf) / FourPi, 4) * ( lambda - 1 ) * pow(log( muR2 / xQ2 / Q2 ), 2) ) )
            + twl * ( pow(beta0qcd(nf) / FourPi, 2) * ( beta3qcd(nf) / pow(FourPi, 4) ) * ( ( 15 - 14 * lambda ) * lambda - 3 )
                      + beta0qcd(nf) / FourPi * ( beta1qcd(nf) / pow(FourPi, 2) )
                      * ( beta2qcd(nf) / pow(FourPi, 3) ) * ( 5 * lambda * ( twl - 3 ) + 3 )
                      + 4 * pow(beta1qcd(nf) / pow(FourPi, 2), 3) * twl2
                      + 2 * pow(beta0qcd(nf) / FourPi, 6) * ( 3 - twl ) * lambda * pow(log( muR2 / xQ2 / Q2 ), 3)
                      + 3 * pow(beta0qcd(nf) / FourPi, 4) * ( beta1qcd(nf) / pow(FourPi, 2) ) * pow(log( muR2 / xQ2 / Q2 ), 2)
                      + 6 * pow(beta0qcd(nf) / FourPi, 2) * lambda * ( twl + 1 )
                      * ( beta0qcd(nf) / FourPi * ( beta2qcd(nf) / pow(FourPi, 3) )
                          - pow(beta1qcd(nf) / pow(FourPi, 2), 2) ) * log( muR2 / xQ2 / Q2 )
                      - 8 * pow(beta0qcd(nf) / FourPi, 6) * ( 4 * twl2 - 6 * lambda + 3 ) * zeta3 ) )
        + ( B3 * ( lambda - 1 ) * lambda ) / ( beta0qcd(nf) * pow(omtwl, 2) )
        + ( B2 * ( ( beta1qcd(nf) / pow(FourPi, 2) ) * lomtwl
                   - 2 * ( lambda - 1 ) * lambda * ( ( beta1qcd(nf) / pow(FourPi, 2) )
                                                     - 2 * pow(beta0qcd(nf) / FourPi, 2) * log( muR2 / xQ2 / Q2 ) ) ) )
        / ( 4 * pow(beta0qcd(nf) / FourPi, 2) * pow(omtwl, 2) )
        + ( M_PI * B1 ) / ( 4 * pow(beta0qcd(nf) / FourPi, 3) * pow(omtwl, 2) )
        * ( 4 * lambda * ( lambda * ( pow(beta1qcd(nf) / pow(FourPi, 2), 2)
                                      - beta0qcd(nf) / FourPi * ( beta2qcd(nf) / pow(FourPi, 3) ) )
                           + pow(beta0qcd(nf) / FourPi, 4) * ( lambda - 1 ) * pow(log( muR2 / xQ2 / Q2 ), 2) )
            - pow(beta1qcd(nf) / pow(FourPi, 2), 2) * pow(lomtwl, 2)
            + 2 * pow(beta0qcd(nf) / FourPi, 2) * ( beta1qcd(nf) / pow(FourPi, 2) ) * lomtwl * log( muR2 / xQ2 / Q2 ) );
      };
    }
  */
}
