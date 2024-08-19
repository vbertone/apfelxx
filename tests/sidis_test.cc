//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/apfelxx.h"
#include "apfel/SIDIS.h"

// Double expression for the leading-order SIDIS F2
class C20nsSidis: public apfel::DoubleExpression
{
public:
  C20nsSidis(): DoubleExpression() {}
  double LocalLocal(double const&, double const&) const override
  {
    return 1;
  }
};

// Double expression for the next-to-leading-order SIDIS FL qq
class CL1nsSidis: public apfel::DoubleExpression
{
public:
  CL1nsSidis(): DoubleExpression() {}
  double RegularRegular(double const& x, double const& z) const override
  {
    return 8 * apfel::CF * x * z;
  }
};

// Double expression for the next-to-leading-order SIDIS FL gq
class CL1gqSidis: public apfel::DoubleExpression
{
public:
  CL1gqSidis(): DoubleExpression() {}
  double RegularRegular(double const& x, double const& z) const override
  {
    return 8 * apfel::CF * x * ( 1 - z );
  }
};

// Double expression for the next-to-leading-order SIDIS FL qg
class CL1qgSidis: public apfel::DoubleExpression
{
public:
  CL1qgSidis(): DoubleExpression() {}
  double RegularRegular(double const& x, double const&) const override
  {
    return 8 * x * ( 1 - x );
  }
};

// Double expression for the next-to-leading-order SIDIS F2 qq
class C21nsSidis: public apfel::DoubleExpression
{
public:
  C21nsSidis(): DoubleExpression() {}

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
    return 2 * apfel::CF * ( 2 + 6 * x * z );
  }
};

// Double expression for the next-to-leading-order SIDIS F2 gq
class C21gqSidis: public apfel::DoubleExpression
{
public:
  C21gqSidis(): DoubleExpression() {}
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
    return 2 * apfel::CF * ( 2 * ( 1 + 3 * x ) - 6 * x * z - ( 1 + x ) / z );
  }
};

// Double expression for the next-to-leading-order SIDIS F2 qg
class C21qgSidis: public apfel::DoubleExpression
{
public:
  C21qgSidis(): DoubleExpression() {}
  double RegularLocal(double const& x, double const& z) const override
  {
    return ( x * x + ( 1 - x ) * ( 1 - x ) ) * log(( 1 - x ) / x) + 2 * x * ( 1 - x ) + ( x * x + ( 1 - x ) * ( 1 - x ) ) * log(1 - z);
  }
  double RegularSingular(double const& x, double const& z) const override
  {
    return ( x * x + ( 1 - x ) * ( 1 - x ) ) / ( 1 - z );
  }
  double RegularRegular(double const& x, double const& z) const override
  {
    return 2 * ( - 1 + 6 * x - 6 * x * x ) + ( x * x + ( 1 - x ) * ( 1 - x ) ) / z;
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
  const auto EvolvedPDFs = BuildDglap(InitializeDglapObjectsQCD(gx, Thresholds), apfel::LHToyPDFs, mu0, PerturbativeOrder, as);
  const auto EvolvedFFs  = BuildDglap(InitializeDglapObjectsQCDT(gz, Thresholds), apfel::LHToyFFs, mu0, PerturbativeOrder, as);

  // Tabulate PDFs and FFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 100, 3};
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedFFs{*EvolvedFFs, 50, 1, 100, 3};

  // Initialize SIDIS objects
  const apfel::SidisObjects so = InitializeSIDIS(gx, gz, Thresholds);

  // Center-of-mass energy
  const double Vs = 17.35;

  // Fine structure constant squared
  const double alpha2 = pow(apfel::alphaem, 2);

  // Define function to compute SIDIS cross section at O(as) using the DoubleObject class
  const std::function<apfel::DoubleObject<apfel::Distribution>(double const&)> CrossSectionDO = [=] (double const& Q) -> apfel::DoubleObject<apfel::Distribution>
  {
    const std::function<double(double const&)> y2 = [=] (double const& x) -> double{ return pow(pow(Q / Vs, 2) / x, 2) / x; };
    const std::function<double(double const&)> yp = [=] (double const& x) -> double{ return ( 1 + pow(1 - pow(Q / Vs, 2) / x, 2) ) / x; };
    const std::function<double(double const&)> iz = [=] (double const& z) -> double{ return 1 / z; };

    // Compute number of active flavours
    const int nf = apfel::NF(Q, Thresholds);

    // Coupling from PDF set
    const double coup = as(Q) / apfel::FourPi;

    // Compute distributions for PDFs and FFs
    const std::map<int, apfel::Distribution> dPDF = apfel::QCDEvToPhys(TabulatedPDFs.Evaluate(Q).GetObjects());
    const std::map<int, apfel::Distribution> dFF  = apfel::QCDEvToPhys(TabulatedFFs.Evaluate(Q).GetObjects());

    apfel::DoubleObject<apfel::Distribution> distqq;
    apfel::DoubleObject<apfel::Distribution> distgq;
    apfel::DoubleObject<apfel::Distribution> distqg;
    for (auto j = - nf; j <= nf; j++)
      {
        // Skip the gluon
        if (j == 0)
          continue;

        distqq.AddTerm({apfel::QCh2[abs(j)-1], dPDF.at(j), dFF.at(j)});
        distgq.AddTerm({apfel::QCh2[abs(j)-1], dPDF.at(j), dFF.at(0)});
        distqg.AddTerm({apfel::QCh2[abs(j)-1], dPDF.at(0), dFF.at(j)});
      }

    // Assemble double distribution for the reduced cross section as Y^+ F2 - y^2 FL
    return ( 4 * M_PI * alpha2 / pow(Q, 3) ) * ( ( ( so.C20qq + coup * so.C21qq ) * distqq + coup * ( so.C21gq * distgq + so.C21qg * distqg ) ).MultiplyBy(yp, iz)
                                                 - ( coup * ( so.CL1qq * distqq + so.CL1gq * distgq + so.CL1qg * distqg ) ).MultiplyBy(y2, iz) );
  };
  const apfel::TabulateObject<apfel::DoubleObject<apfel::Distribution>> TabCrossSectionDO{CrossSectionDO, 50, 1, 10, 3, Thresholds};

  // Compute relevant double operators
  apfel::Timer t;
  const apfel::DoubleOperator O20qq{gx, gz, C20nsSidis{}};
  const apfel::DoubleOperator OL1qq{gx, gz, CL1nsSidis{}};
  const apfel::DoubleOperator OL1gq{gx, gz, CL1gqSidis{}};
  const apfel::DoubleOperator OL1qg{gx, gz, CL1qgSidis{}};
  const apfel::DoubleOperator O21qq{gx, gz, C21nsSidis{}};
  const apfel::DoubleOperator O21gq{gx, gz, C21gqSidis{}};
  const apfel::DoubleOperator O21qg{gx, gz, C21qgSidis{}};
  t.stop();

  // Define function to compute SIDIS cross section at O(as) using the DoubleDistribution (and DoubleOperator) class
  const std::function<apfel::DoubleDistribution(double const&)> CrossSectionDD = [=, &gx, &gz] (double const& Q) -> apfel::DoubleDistribution
  {
    const std::function<double(double const&, double const&)> y2 = [=] (double const& x, double const& z) -> double{ return pow(pow(Q / Vs, 2) / x, 2) / x / z; };
    const std::function<double(double const&, double const&)> yp = [=] (double const& x, double const& z) -> double{ return ( 1 + pow(1 - pow(Q / Vs, 2) / x, 2) ) / x / z; };

    // Compute number of active flavours
    const int nf = apfel::NF(Q, Thresholds);

    // Coupling from PDF set
    const double coup = as(Q) / apfel::FourPi;

    // Compute distributions for PDFs and FFs
    const std::map<int, apfel::Distribution> dPDF = apfel::QCDEvToPhys(TabulatedPDFs.Evaluate(Q).GetObjects());
    const std::map<int, apfel::Distribution> dFF  = apfel::QCDEvToPhys(TabulatedFFs.Evaluate(Q).GetObjects());

    // Initialize double distributions
    apfel::DoubleDistribution distqq{gx, gz};
    apfel::DoubleDistribution distgq{gx, gz};
    apfel::DoubleDistribution distqg{gx, gz};
    for (auto j = - nf; j <= nf; j++)
      {
        // Skip the gluon
        if (j == 0)
          continue;

        distqq += apfel::QCh2[abs(j)-1] * apfel::DoubleDistribution{dPDF.at(j), dFF.at(j)};
        distgq += apfel::QCh2[abs(j)-1] * apfel::DoubleDistribution{dPDF.at(j), dFF.at(0)};
        distqg += apfel::QCh2[abs(j)-1] * apfel::DoubleDistribution{dPDF.at(0), dFF.at(j)};
      }

    // Assemble double distribution for the reduced cross section as Y^+ F2 - y^2 FL
    return ( 4 * M_PI * alpha2 / pow(Q, 3) ) * ( yp * ( ( O20qq + coup * O21qq ) * distqq + coup * ( O21gq * distgq + O21qg * distqg ) )
                                                 - coup * ( y2 * ( OL1qq * distqq + OL1gq * distgq + OL1qg * distqg ) ) );
  };
  const apfel::TabulateObject<apfel::DoubleDistribution> TabCrossSectionDD{CrossSectionDD, 50, 1, 10, 3, Thresholds};

  // Define kinematics
  const double x  = 0.01;
  const double z  = 0.4;
  const double Q  = 2;
  const double op = 1 + apfel::eps3;
  const double om = 1 - apfel::eps3;

  // Print value
  std::cout << std::scientific << std::endl;
  t.start();
  std::cout << "Reduced SIDIS cross section (Q = " << Q << " GeV, x = " << x << ", z = " << z << "): " << CrossSectionDO(Q).Evaluate(x, z) << std::endl;
  std::cout << "Reduced SIDIS cross section (Q = " << Q << " GeV, x = " << x << ", z = " << z << "): " << TabCrossSectionDO.Evaluate(Q).Evaluate(x, z) << std::endl;
  std::cout << "Reduced SIDIS cross section (Q = " << Q << " GeV, x = " << x << ", z = " << z << "): " << TabCrossSectionDO.Integrate(Q*om, Q*op).Integrate(x*om, x*op, z*om, z*op) / ( ( op - om ) * Q ) / ( ( op - om ) * x ) / ( ( op - om ) * z ) << std::endl;
  t.stop();
  t.start();
  std::cout << "Reduced SIDIS cross section (Q = " << Q << " GeV, x = " << x << ", z = " << z << "): " << CrossSectionDD(Q).Evaluate(x, z) << std::endl;
  std::cout << "Reduced SIDIS cross section (Q = " << Q << " GeV, x = " << x << ", z = " << z << "): " << CrossSectionDD(Q).Evaluate1(x).Evaluate(z) << std::endl;
  std::cout << "Reduced SIDIS cross section (Q = " << Q << " GeV, x = " << x << ", z = " << z << "): " << CrossSectionDD(Q).Evaluate2(z).Evaluate(x) << std::endl;
  std::cout << "Reduced SIDIS cross section (Q = " << Q << " GeV, x = " << x << ", z = " << z << "): " << TabCrossSectionDD.Evaluate(Q).Evaluate(x, z) << std::endl;
  std::cout << "Reduced SIDIS cross section (Q = " << Q << " GeV, x = " << x << ", z = " << z << "): " << TabCrossSectionDD.Integrate(Q*om, Q*op).Integrate(x*om, x*op, z*om, z*op) / ( ( op - om ) * Q ) / ( ( op - om ) * x ) / ( ( op - om ) * z ) << std::endl;
  t.stop();
  std::cout << std::endl;

  return 0;
}
