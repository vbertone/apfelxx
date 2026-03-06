// APFEL++ 2017

#include "apfel/sidisbuilder.h"
#include "apfel/sidiscoefficientfunctionsunp.h"
#include "apfel/sidiscoefficientfunctionspol.h"
#include "apfel/doubleexpression.h"
#include "apfel/distribution.h"
#include "apfel/constants.h"
#include "apfel/tools.h"
#include "apfel/messages.h"
#include "apfel/timer.h"

namespace apfel
{
  SidisNNLOObjects InitializeSidisObjects(Grid                const& gx,
                                          Grid                const& gz,
                                          std::vector<double> const& Thresholds,
                                          double              const& IntEps)
  {
    report("Initializing exact NNLO SIDIS hard cross sections... ");
    Timer t;

    // Determine the range of active flavour numbers from the threshold vector.
    // Thresholds with value <= 0 are treated as massless (already active at Q = 0).
    int nfi = 0;
    int nff = static_cast<int>(Thresholds.size());
    for (auto const& v : Thresholds)
      if (v <= 0)
        nfi++;

    SidisNNLOObjects obj;

    // =========================================================================
    // LO: delta(1-x) * delta(1-z)
    // Used by FT and G1 (FL has no LO contribution).
    // =========================================================================
    const DoubleOperator OLO{gx, gz, DoubleIdentity{}, IntEps};
    obj.FT.emplace("DoubleIdentity", OLO);
    obj.G1.emplace("DoubleIdentity", OLO);

    // =========================================================================
    // NLO coefficient functions
    // =========================================================================

    // --- NLO FT ---
    obj.FT.emplace("C1TQ2Q", DoubleOperator{gx, gz, C1TQ2Q{}, IntEps});
    obj.FT.emplace("C1TQ2G", DoubleOperator{gx, gz, C1TQ2G{}, IntEps});
    obj.FT.emplace("C1TG2Q", DoubleOperator{gx, gz, C1TG2Q{}, IntEps});

    // --- NLO FL ---
    obj.FL.emplace("C1LQ2Q", DoubleOperator{gx, gz, C1LQ2Q{}, IntEps});
    obj.FL.emplace("C1LQ2G", DoubleOperator{gx, gz, C1LQ2G{}, IntEps});
    obj.FL.emplace("C1LG2Q", DoubleOperator{gx, gz, C1LG2Q{}, IntEps});

    // --- NLO G1 ---
    obj.G1.emplace("DC1Q2Q", DoubleOperator{gx, gz, DC1Q2Q{}, IntEps});
    obj.G1.emplace("DC1Q2G", DoubleOperator{gx, gz, DC1Q2G{}, IntEps});
    obj.G1.emplace("DC1G2Q", DoubleOperator{gx, gz, DC1G2Q{}, IntEps});

    // =========================================================================
    // NNLO nf-independent coefficient functions
    // =========================================================================

    // --- NNLO FT (nf-independent) ---
    obj.FT.emplace("C2TQ2G",   DoubleOperator{gx, gz, C2TQ2G{},   IntEps});
    obj.FT.emplace("C2TG2Q",   DoubleOperator{gx, gz, C2TG2Q{},   IntEps});
    obj.FT.emplace("C2TG2G",   DoubleOperator{gx, gz, C2TG2G{},   IntEps});
    obj.FT.emplace("C2TQ2QPS", DoubleOperator{gx, gz, C2TQ2QPS{}, IntEps});
    obj.FT.emplace("C2TQ2QB",  DoubleOperator{gx, gz, C2TQ2QB{},  IntEps});
    obj.FT.emplace("C2TQ2QP1", DoubleOperator{gx, gz, C2TQ2QP1{}, IntEps});
    obj.FT.emplace("C2TQ2QP2", DoubleOperator{gx, gz, C2TQ2QP2{}, IntEps});
    obj.FT.emplace("C2TQ2QP3", DoubleOperator{gx, gz, C2TQ2QP3{}, IntEps});

    // --- NNLO FL (nf-independent) ---
    obj.FL.emplace("C2LQ2G",   DoubleOperator{gx, gz, C2LQ2G{},   IntEps});
    obj.FL.emplace("C2LG2Q",   DoubleOperator{gx, gz, C2LG2Q{},   IntEps});
    obj.FL.emplace("C2LG2G",   DoubleOperator{gx, gz, C2LG2G{},   IntEps});
    obj.FL.emplace("C2LQ2QPS", DoubleOperator{gx, gz, C2LQ2QPS{}, IntEps});
    obj.FL.emplace("C2LQ2QB",  DoubleOperator{gx, gz, C2LQ2QB{},  IntEps});
    obj.FL.emplace("C2LQ2QP1", DoubleOperator{gx, gz, C2LQ2QP1{}, IntEps});
    obj.FL.emplace("C2LQ2QP2", DoubleOperator{gx, gz, C2LQ2QP2{}, IntEps});
    obj.FL.emplace("C2LQ2QP3", DoubleOperator{gx, gz, C2LQ2QP3{}, IntEps});

    // --- NNLO G1 (nf-independent) ---
    obj.G1.emplace("DC2G2G",   DoubleOperator{gx, gz, DC2G2G{},   IntEps});
    obj.G1.emplace("DC2Q2QPS", DoubleOperator{gx, gz, DC2Q2QPS{}, IntEps});
    obj.G1.emplace("DC2Q2QB",  DoubleOperator{gx, gz, DC2Q2QB{},  IntEps});
    obj.G1.emplace("DC2Q2QP1", DoubleOperator{gx, gz, DC2Q2QP1{}, IntEps});
    obj.G1.emplace("DC2Q2QP2", DoubleOperator{gx, gz, DC2Q2QP2{}, IntEps});
    obj.G1.emplace("DC2Q2QP3", DoubleOperator{gx, gz, DC2Q2QP3{}, IntEps});

    // =========================================================================
    // NNLO nf-dependent coefficient functions
    // Loop over all active-flavour numbers in the range [nfi, nff].
    // =========================================================================
    for (int nf = nfi; nf <= nff; nf++)
      {
        const std::string snf = "_nf" + std::to_string(nf);

        // --- NNLO FT non-singlet (nf-dependent) ---
        obj.FT.emplace("C2TQ2QNS" + snf, DoubleOperator{gx, gz, C2TQ2QNS{nf}, IntEps});

        // --- NNLO FL non-singlet (nf-dependent) ---
        obj.FL.emplace("C2LQ2QNS" + snf, DoubleOperator{gx, gz, C2LQ2QNS{nf}, IntEps});

        // --- NNLO G1 non-singlet (nf-dependent) ---
        obj.G1.emplace("DC2Q2QNS" + snf, DoubleOperator{gx, gz, DC2Q2QNS{nf}, IntEps});

        // --- NNLO G1 gluon-in-quark (nf-dependent for polarised) ---
        obj.G1.emplace("DC2Q2G" + snf, DoubleOperator{gx, gz, DC2Q2G{nf}, IntEps});

        // --- NNLO G1 quark-in-gluon (nf-dependent for polarised) ---
        obj.G1.emplace("DC2G2Q" + snf, DoubleOperator{gx, gz, DC2G2Q{nf}, IntEps});
      }

    t.stop();

    return obj;
  }

  // ===========================================================================
  // Helper: build the 9 partonic DoubleDistribution channels from PDF/FF maps.
  // Shared by all three Build functions.
  // ===========================================================================
  namespace
  {
    struct SidisChannels
    {
      DoubleDistribution distns;
      DoubleDistribution distgq;
      DoubleDistribution distqg;
      DoubleDistribution distgg;
      DoubleDistribution distps;
      DoubleDistribution distqbq;
      DoubleDistribution distqpq1;
      DoubleDistribution distqpq2;
      DoubleDistribution distqpq3;
    };

    SidisChannels BuildChannels(Grid                             const& gx,
                                Grid                             const& gz,
                                std::map<int, Distribution>      const& dPDF,
                                std::map<int, Distribution>      const& dFF,
                                int                              const  nf)
    {
      const double sumch2 = SumCh2[nf];

      Distribution eq2PDFs{gx};
      Distribution eq2FFs{gz};
      for (int j = -nf; j <= nf; j++)
        {
          if (j == 0) continue;
          eq2PDFs += QCh2[std::abs(j) - 1] * dPDF.at(j);
          eq2FFs  += QCh2[std::abs(j) - 1] * dFF.at(j);
        }

      SidisChannels ch{
        {gx, gz}, {eq2PDFs, dFF.at(21)}, {dPDF.at(21), eq2FFs},
        {sumch2 * dPDF.at(21), dFF.at(21)},
        {gx, gz}, {gx, gz}, {gx, gz}, {gx, gz}, {gx, gz}
      };

      for (int j = -nf; j <= nf; j++)
        {
          if (j == 0) continue;
          ch.distns  += DoubleDistribution{QCh2[std::abs(j) - 1] * dPDF.at(j), dFF.at(j)};
          ch.distps  += DoubleDistribution{sumch2 * dPDF.at(j), dFF.at(j)};
          ch.distqbq += DoubleDistribution{QCh2[std::abs(j) - 1] * dPDF.at(j), dFF.at(-j)};
          for (int i = -nf; i <= nf; i++)
            {
              if (i == 0 || i == j || i == -j) continue;
              ch.distqpq1 += DoubleDistribution{QCh2[std::abs(j) - 1] * dPDF.at(j), dFF.at(i)};
              ch.distqpq2 += DoubleDistribution{dPDF.at(j), QCh2[std::abs(i) - 1] * dFF.at(i)};
              ch.distqpq3 += DoubleDistribution{((j > 0 ? 1 : -1) * QCh[std::abs(j) - 1]) * dPDF.at(j),
                                                ((i > 0 ? 1 : -1) * QCh[std::abs(i) - 1]) * dFF.at(i)};
            }
        }
      return ch;
    }
  }

  // ===========================================================================
  // BuildSidisUnpFT
  // ===========================================================================
  std::function<DoubleDistribution(double const&)> BuildSidisUnpFT(
      SidisNNLOObjects                                                          const& obj,
      std::function<std::map<int, double>(double const&, double const&)> const& InPDFs,
      std::function<std::map<int, double>(double const&, double const&)> const& InFFs,
      std::function<double(double const&)>                               const& Alphas,
      std::vector<double>                                                const& Thresholds,
      int                                                                const& PerturbativeOrder)
  {
    return [=, &obj] (double const& Q) -> DoubleDistribution
      {
        Grid const& gx = obj.FT.begin()->second.GetFirstGrid();
        Grid const& gz = obj.FT.begin()->second.GetSecondGrid();

        const std::function<double(double const&, double const&)> ixz =
          [] (double const& x, double const& z) -> double { return 1. / x / z; };

        const int    nf    = NF(Q, Thresholds);
        const double coup  = Alphas(Q) / FourPi;
        const double coup2 = coup * coup;

        const std::map<int, Distribution> dPDF = DistributionMap(gx, InPDFs, Q);
        const std::map<int, Distribution> dFF  = DistributionMap(gz, InFFs,  Q);

        const SidisChannels ch = BuildChannels(gx, gz, dPDF, dFF, nf);

        const DoubleOperator OZero{gx, gz, DoubleNull{}};
        const std::string    snf = "_nf" + std::to_string(nf);

        DoubleOperator Ons   = OZero;
        DoubleOperator Ogq   = OZero;
        DoubleOperator Oqg   = OZero;
        DoubleOperator Ogg   = OZero;
        DoubleOperator Ops   = OZero;
        DoubleOperator Oqbq  = OZero;
        DoubleOperator Oqpq1 = OZero;
        DoubleOperator Oqpq2 = OZero;
        DoubleOperator Oqpq3 = OZero;

        if (PerturbativeOrder >= 0)
          Ons += obj.FT.at("DoubleIdentity");

        if (PerturbativeOrder >= 1)
          {
            Ons += coup * obj.FT.at("C1TQ2Q");
            Ogq += coup * obj.FT.at("C1TQ2G");
            Oqg += coup * obj.FT.at("C1TG2Q");
          }

        if (PerturbativeOrder >= 2)
          {
            Ons   += coup2 * obj.FT.at("C2TQ2QNS" + snf);
            Ogq   += coup2 * obj.FT.at("C2TQ2G");
            Oqg   += coup2 * obj.FT.at("C2TG2Q");
            Ogg   += coup2 * obj.FT.at("C2TG2G");
            Ops   += coup2 * obj.FT.at("C2TQ2QPS");
            Oqbq  += coup2 * obj.FT.at("C2TQ2QB");
            Oqpq1 += coup2 * obj.FT.at("C2TQ2QP1");
            Oqpq2 += coup2 * obj.FT.at("C2TQ2QP2");
            Oqpq3 += coup2 * obj.FT.at("C2TQ2QP3");
          }

        return ixz * (Ons   * ch.distns   + Ogq  * ch.distgq  + Oqg  * ch.distqg +
                      Ogg   * ch.distgg   + Ops  * ch.distps  + Oqbq * ch.distqbq +
                      Oqpq1 * ch.distqpq1 + Oqpq2 * ch.distqpq2 + Oqpq3 * ch.distqpq3);
      };
  }

  // ===========================================================================
  // BuildSidisUnpFL
  // ===========================================================================
  std::function<DoubleDistribution(double const&)> BuildSidisUnpFL(
      SidisNNLOObjects                                                          const& obj,
      std::function<std::map<int, double>(double const&, double const&)> const& InPDFs,
      std::function<std::map<int, double>(double const&, double const&)> const& InFFs,
      std::function<double(double const&)>                               const& Alphas,
      std::vector<double>                                                const& Thresholds,
      int                                                                const& PerturbativeOrder)
  {
    return [=, &obj] (double const& Q) -> DoubleDistribution
      {
        Grid const& gx = obj.FL.begin()->second.GetFirstGrid();
        Grid const& gz = obj.FL.begin()->second.GetSecondGrid();

        const std::function<double(double const&, double const&)> ixz =
          [] (double const& x, double const& z) -> double { return 1. / x / z; };

        const int    nf    = NF(Q, Thresholds);
        const double coup  = Alphas(Q) / FourPi;
        const double coup2 = coup * coup;

        const std::map<int, Distribution> dPDF = DistributionMap(gx, InPDFs, Q);
        const std::map<int, Distribution> dFF  = DistributionMap(gz, InFFs,  Q);

        const SidisChannels ch = BuildChannels(gx, gz, dPDF, dFF, nf);

        const DoubleOperator OZero{gx, gz, DoubleNull{}};
        const std::string    snf = "_nf" + std::to_string(nf);

        DoubleOperator Ons   = OZero;
        DoubleOperator Ogq   = OZero;
        DoubleOperator Oqg   = OZero;
        DoubleOperator Ogg   = OZero;
        DoubleOperator Ops   = OZero;
        DoubleOperator Oqbq  = OZero;
        DoubleOperator Oqpq1 = OZero;
        DoubleOperator Oqpq2 = OZero;
        DoubleOperator Oqpq3 = OZero;

        // FL has no LO; NLO starts at O(αs)
        if (PerturbativeOrder >= 1)
          {
            Ons += coup * obj.FL.at("C1LQ2Q");
            Ogq += coup * obj.FL.at("C1LQ2G");
            Oqg += coup * obj.FL.at("C1LG2Q");
          }

        if (PerturbativeOrder >= 2)
          {
            Ons   += coup2 * obj.FL.at("C2LQ2QNS" + snf);
            Ogq   += coup2 * obj.FL.at("C2LQ2G");
            Oqg   += coup2 * obj.FL.at("C2LG2Q");
            Ogg   += coup2 * obj.FL.at("C2LG2G");
            Ops   += coup2 * obj.FL.at("C2LQ2QPS");
            Oqbq  += coup2 * obj.FL.at("C2LQ2QB");
            Oqpq1 += coup2 * obj.FL.at("C2LQ2QP1");
            Oqpq2 += coup2 * obj.FL.at("C2LQ2QP2");
            Oqpq3 += coup2 * obj.FL.at("C2LQ2QP3");
          }

        return ixz * (Ons   * ch.distns   + Ogq  * ch.distgq  + Oqg  * ch.distqg +
                      Ogg   * ch.distgg   + Ops  * ch.distps  + Oqbq * ch.distqbq +
                      Oqpq1 * ch.distqpq1 + Oqpq2 * ch.distqpq2 + Oqpq3 * ch.distqpq3);
      };
  }

  // ===========================================================================
  // BuildSidisG1
  // ===========================================================================
  std::function<DoubleDistribution(double const&)> BuildSidisG1(
      SidisNNLOObjects                                                          const& obj,
      std::function<std::map<int, double>(double const&, double const&)> const& InPDFs,
      std::function<std::map<int, double>(double const&, double const&)> const& InFFs,
      std::function<double(double const&)>                               const& Alphas,
      std::vector<double>                                                const& Thresholds,
      int                                                                const& PerturbativeOrder)
  {
    return [=, &obj] (double const& Q) -> DoubleDistribution
      {
        Grid const& gx = obj.G1.begin()->second.GetFirstGrid();
        Grid const& gz = obj.G1.begin()->second.GetSecondGrid();

        const std::function<double(double const&, double const&)> ixz =
          [] (double const& x, double const& z) -> double { return 1. / x / z; };

        const int    nf    = NF(Q, Thresholds);
        const double coup  = Alphas(Q) / FourPi;
        const double coup2 = coup * coup;

        const std::map<int, Distribution> dPDF = DistributionMap(gx, InPDFs, Q);
        const std::map<int, Distribution> dFF  = DistributionMap(gz, InFFs,  Q);

        const SidisChannels ch = BuildChannels(gx, gz, dPDF, dFF, nf);

        const DoubleOperator OZero{gx, gz, DoubleNull{}};
        const std::string    snf = "_nf" + std::to_string(nf);

        DoubleOperator Ons   = OZero;
        DoubleOperator Ogq   = OZero;
        DoubleOperator Oqg   = OZero;
        DoubleOperator Ogg   = OZero;
        DoubleOperator Ops   = OZero;
        DoubleOperator Oqbq  = OZero;
        DoubleOperator Oqpq1 = OZero;
        DoubleOperator Oqpq2 = OZero;
        DoubleOperator Oqpq3 = OZero;

        if (PerturbativeOrder >= 0)
          Ons += obj.G1.at("DoubleIdentity");

        if (PerturbativeOrder >= 1)
          {
            Ons += coup * obj.G1.at("DC1Q2Q");
            Ogq += coup * obj.G1.at("DC1Q2G");
            Oqg += coup * obj.G1.at("DC1G2Q");
          }

        if (PerturbativeOrder >= 2)
          {
            // For G1: ns, gq, and qg are all nf-dependent
            Ons   += coup2 * obj.G1.at("DC2Q2QNS" + snf);
            Ogq   += coup2 * obj.G1.at("DC2Q2G"   + snf);
            Oqg   += coup2 * obj.G1.at("DC2G2Q"   + snf);
            Ogg   += coup2 * obj.G1.at("DC2G2G");
            Ops   += coup2 * obj.G1.at("DC2Q2QPS");
            Oqbq  += coup2 * obj.G1.at("DC2Q2QB");
            Oqpq1 += coup2 * obj.G1.at("DC2Q2QP1");
            Oqpq2 += coup2 * obj.G1.at("DC2Q2QP2");
            Oqpq3 += coup2 * obj.G1.at("DC2Q2QP3");
          }

        return ixz * (Ons   * ch.distns   + Ogq  * ch.distgq  + Oqg  * ch.distqg +
                      Ogg   * ch.distgg   + Ops  * ch.distps  + Oqbq * ch.distqbq +
                      Oqpq1 * ch.distqpq1 + Oqpq2 * ch.distqpq2 + Oqpq3 * ch.distqpq3);
      };
  }
}
