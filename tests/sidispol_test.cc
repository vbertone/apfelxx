//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/apfelxx.h"
#include "apfel/SIDISpol.h"

// Double expression for the leading-order SIDIS G1
class G10nsSidis: public apfel::DoubleExpression
{
public:
  G10nsSidis(): DoubleExpression() {}
  double LocalLocal(double const&, double const&) const override
  {
    return 1;
  }
};

// Double expression for the next-to-leading-order SIDIS G1 qq
class G11nsSidis: public apfel::DoubleExpression
{
public:
  G11nsSidis(): DoubleExpression() {}
  double LocalLocal(double const& x, double const& z) const override
  {
    return 2 * apfel::CF * ( - 8 + pow(log(1 - z), 2) + pow(log(1 - x), 2) + 2 * log(1 - x) * log(1 - z) );
  }
  double LocalSingular(double const& x, double const& z) const override
  {
    return 2 * apfel::CF * ( 2 * log(1 - z) / ( 1 - z ) + 2 * log( 1 - x ) / ( 1 - z ) );
  }
  double LocalRegular(double const& x, double const& z) const override
  {
    return 2 * apfel::CF * ( ( 1 + z * z ) * log(z) / ( 1 - z ) + 1 - z - ( 1 + z ) * log(1 - z) - log(1 - x) * ( 1 + z )  );
  }
  double SingularLocal(double const& x, double const& z) const override
  {
    return 2 * apfel::CF * ( 2 * log(1 - x) / ( 1 - x ) + 2 * log(1 - z) / ( 1 - x ) );
  }
  double SingularSingular(double const& x, double const& z) const override
  {
    return 4 * apfel::CF / ( 1 - x ) / ( 1 - z );
  }
  double SingularRegular(double const& x, double const& z) const override
  {
    return - 2 * apfel::CF * ( 1 + z ) / ( 1 - x );
  }
  double RegularLocal(double const& x, double const& z) const override
  {
    return 2 * apfel::CF * ( - ( 1 + x * x ) * log(x) / ( 1 - x ) + 1 - x - ( 1 + x ) * log(1 - x) - log(1 - z) * ( 1 + x ) );
  }
  double RegularSingular(double const& x, double const& z) const override
  {
    return - 2 * apfel::CF * ( 1 + x ) / ( 1 - z );
  }
  double RegularRegular(double const& x, double const& z) const override
  {
    return 2 * apfel::CF * ( 2 * x + 2 * z );
  }
};

// Double expression for the next-to-leading-order SIDIS F2 gq
class G11gqSidis: public apfel::DoubleExpression
{
public:
  G11gqSidis(): DoubleExpression() {}
  double LocalRegular(double const& x, double const& z) const override
  {
    return 2 * apfel::CF * ( ( 1 + ( 1 - z ) * ( 1 - z ) ) * log(z * ( 1 - z )) / z + z + log(1 - x) * ( 1 + ( 1 - z ) * ( 1 - z ) ) / z );
  }
  double SingularRegular(double const& x, double const& z) const override
  {
    return 2 * apfel::CF * ( 1 + ( 1 - z ) * ( 1 - z ) ) / z / ( 1 - x );
  }
  double RegularRegular(double const& x, double const& z) const override
  {
    return 2 * apfel::CF * ( ( 1 + x ) * ( 2 - 1 / z ) - 2 * z );
  }
};

// Double expression for the next-to-leading-order SIDIS F2 qg
class G11qgSidis: public apfel::DoubleExpression
{
public:
  G11qgSidis(): DoubleExpression() {}
  double RegularLocal(double const& x, double const& z) const override
  {
    return ( x * x - ( 1 - x ) * ( 1 - x ) ) * log(( 1 - x ) / x) + 2 * ( 1 - x ) + ( x * x - ( 1 - x ) * ( 1 - x ) ) * log(1 - z);
  }
  double RegularSingular(double const& x, double const& z) const override
  {
    return ( x * x - ( 1 - x ) * ( 1 - x ) ) / ( 1 - z );
  }
  double RegularRegular(double const& x, double const& z) const override
  {
    return ( x * x - ( 1 - x ) * ( 1 - x ) ) * ( 1 / z - 2 );
  }
};

int main()
{
  // Initial scale
  const double mu0 = sqrt(2);

  // Vectors of masses and thresholds
  const std::vector<double> Thresholds = {0, 0, 0, sqrt(2), 4.5, 175};

  // Perturbative order
  const int PerturbativeOrder = 1;

  // Running coupling
  apfel::AlphaQCD a{0.35, sqrt(2), Thresholds, PerturbativeOrder};
  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // x-space and z-space grids
  const apfel::Grid gx{{apfel::SubGrid{100, 1e-3, 3}, apfel::SubGrid{60, 1e-1, 3}, apfel::SubGrid{50, 8e-1, 3}}};
  const apfel::Grid gz{{apfel::SubGrid{100, 1e-2, 3}, apfel::SubGrid{60, 2e-1, 3}, apfel::SubGrid{50, 8e-1, 3}}};

  // Construct the DGLAP objects
  const auto EvolvedpPDFs = BuildDglap(InitializeDglapObjectsQCDpol(gx, Thresholds), apfel::LHToyPDFsPol, mu0, PerturbativeOrder, as);
  const auto EvolvedFFs   = BuildDglap(InitializeDglapObjectsQCDT(gz, Thresholds), apfel::LHToyFFs, mu0, PerturbativeOrder, as);

  // Tabulate PDFs and FFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedpPDFs{*EvolvedpPDFs, 50, 1, 100, 3};
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedFFs{*EvolvedFFs, 50, 1, 100, 3};

  // Initialize SIDIS objects
  const apfel::SidisPolObjects spo = InitializeSIDISpol(gx, gz, Thresholds);

  // Define function to compute polarised SIDIS cross section at O(as) using the DoubleObject class
  const std::function<apfel::DoubleObject<apfel::Distribution>(double const&)> CrossSectionDOpol = [=] (double const& Q) -> apfel::DoubleObject<apfel::Distribution>
  {
    // Compute number of active flavours
    const int nf = apfel::NF(Q, Thresholds);

    // Coupling from PDF set
    const double coup = as(Q) / apfel::FourPi;

    // Compute distributions for PDFs and FFs
    const std::map<int, apfel::Distribution> dpPDF = apfel::QCDEvToPhys(TabulatedpPDFs.Evaluate(Q).GetObjects());
    const std::map<int, apfel::Distribution> dFF   = apfel::QCDEvToPhys(TabulatedFFs.Evaluate(Q).GetObjects());

    apfel::DoubleObject<apfel::Distribution> distqq;
    apfel::DoubleObject<apfel::Distribution> distgq;
    apfel::DoubleObject<apfel::Distribution> distqg;
    for (auto j = - nf; j <= nf; j++)
      {
        // Skip the gluon
        if (j == 0)
          continue;

        distqq.AddTerm({apfel::QCh2[abs(j)-1], dpPDF.at(j), dFF.at(j)});
        distgq.AddTerm({apfel::QCh2[abs(j)-1], dpPDF.at(j), dFF.at(0)});
        distqg.AddTerm({apfel::QCh2[abs(j)-1], dpPDF.at(0), dFF.at(j)});
      }

    // Assemble double distribution for G1
    return ( spo.G10qq + coup * spo.G11qq ) * distqq + coup * ( spo.G11gq * distgq + spo.G11qg * distqg );
  };
  const apfel::TabulateObject<apfel::DoubleObject<apfel::Distribution>> TabCrossSectionDOpol{CrossSectionDOpol, 50, 1, 10, 3, Thresholds};

  // Compute relevant double polarised operators
  apfel::Timer t;
  const apfel::DoubleOperator Q20qq{gx, gz, G10nsSidis{}};
  const apfel::DoubleOperator Q21qq{gx, gz, G11nsSidis{}};
  const apfel::DoubleOperator Q21gq{gx, gz, G11gqSidis{}};
  const apfel::DoubleOperator Q21qg{gx, gz, G11qgSidis{}};
  t.stop();

  // Define function to compute polarised SIDIS cross section at O(as) using the DoubleDistribution (and DoubleOperator) class
  const std::function<apfel::DoubleDistribution(double const&)> CrossSectionDDpol = [=, &gx, &gz] (double const& Q) -> apfel::DoubleDistribution
  {
    // Compute number of active flavours
    const int nf = apfel::NF(Q, Thresholds);

    // Coupling from PDF set
    const double coup = as(Q) / apfel::FourPi;

    // Compute distributions for PDFs and FFs
    const std::map<int, apfel::Distribution> dpPDF = apfel::QCDEvToPhys(TabulatedpPDFs.Evaluate(Q).GetObjects());
    const std::map<int, apfel::Distribution> dFF   = apfel::QCDEvToPhys(TabulatedFFs.Evaluate(Q).GetObjects());

    // Initialize double distributions
    apfel::DoubleDistribution distqq{gx, gz};
    apfel::DoubleDistribution distgq{gx, gz};
    apfel::DoubleDistribution distqg{gx, gz};
    for (auto j = - nf; j <= nf; j++)
      {
        // Skip the gluon
        if (j == 0)
          continue;

        distqq += apfel::QCh2[abs(j)-1] * apfel::DoubleDistribution{dpPDF.at(j), dFF.at(j)};
        distgq += apfel::QCh2[abs(j)-1] * apfel::DoubleDistribution{dpPDF.at(j), dFF.at(0)};
        distqg += apfel::QCh2[abs(j)-1] * apfel::DoubleDistribution{dpPDF.at(0), dFF.at(j)};
      }

    // Assemble double distribution for G1
    return ( Q20qq + coup * Q21qq ) * distqq + coup * ( Q21gq * distgq + Q21qg * distqg );
  };
  const apfel::TabulateObject<apfel::DoubleDistribution> TabCrossSectionDDpol{CrossSectionDDpol, 50, 1, 10, 3, Thresholds};

  // Define kinematics
  const double x  = 0.01;
  const double z  = 0.4;
  const double Q  = 2;
  const double op = 1 + apfel::eps3;
  const double om = 1 - apfel::eps3;
  const double dn =  pow(op - om, 3) * Q * x * z;

  // Print value
  std::cout << std::scientific << std::endl;
  t.start();
  std::cout << "Reduced pol. SIDIS cross section (Q = " << Q << " GeV, x = " << x << ", z = " << z << "): " << CrossSectionDOpol(Q).Evaluate(x, z) << std::endl;
  std::cout << "Reduced pol. SIDIS cross section (Q = " << Q << " GeV, x = " << x << ", z = " << z << "): " << TabCrossSectionDOpol.Evaluate(Q).Evaluate(x, z) << std::endl;
  std::cout << "Reduced pol. SIDIS cross section (Q = " << Q << " GeV, x = " << x << ", z = " << z << "): " << TabCrossSectionDOpol.Integrate(Q*om, Q*op).Integrate(x*om, x*op, z*om, z*op) / dn << std::endl;
  t.stop();
  t.start();
  std::cout << "Reduced pol. SIDIS cross section (Q = " << Q << " GeV, x = " << x << ", z = " << z << "): " << CrossSectionDDpol(Q).Evaluate(x, z) << std::endl;
  std::cout << "Reduced pol. SIDIS cross section (Q = " << Q << " GeV, x = " << x << ", z = " << z << "): " << CrossSectionDDpol(Q).Evaluate1(x).Evaluate(z) << std::endl;
  std::cout << "Reduced pol. SIDIS cross section (Q = " << Q << " GeV, x = " << x << ", z = " << z << "): " << CrossSectionDDpol(Q).Evaluate2(z).Evaluate(x) << std::endl;
  std::cout << "Reduced pol. SIDIS cross section (Q = " << Q << " GeV, x = " << x << ", z = " << z << "): " << TabCrossSectionDDpol.Evaluate(Q).Evaluate(x, z) << std::endl;
  std::cout << "Reduced pol. SIDIS cross section (Q = " << Q << " GeV, x = " << x << ", z = " << z << "): " << TabCrossSectionDDpol.Integrate(Q*om, Q*op).Integrate(x*om, x*op, z*om, z*op) / dn << std::endl;
  std::cout << "Reduced pol. SIDIS cross section (Q = " << Q << " GeV, x = " << x << ", z = " << z << "): " << TabCrossSectionDDpol.Integrate(Q*om, Q*op).Integrate1(x*om, x*op).Integrate(z*om, z*op) / dn << std::endl;
  std::cout << "Reduced pol. SIDIS cross section (Q = " << Q << " GeV, x = " << x << ", z = " << z << "): " << TabCrossSectionDDpol.Integrate(Q*om, Q*op).Integrate2(z*om, z*op).Integrate(x*om, x*op) / dn << std::endl;
  t.stop();
  std::cout << std::endl;

  return 0;
}
