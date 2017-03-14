//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/dglap.h"
#include "apfel/grid.h"
#include "apfel/distribution.h"

#include <functional>
#include <memory>

using std::function;
using std::unique_ptr;

namespace apfel
{
  /**
   * @brief The DglapQCD class.
   */
  class DglapQCD
  {
  public:

    DglapQCD() = delete;

    /**
     * @brief AlphaQCD default constructor.
     * @param AlphaRef the reference value of the coupling.
     * @param MuRef the reference value of the scale.
     * @param Masses vector of masses.
     * @param Thresholds vector of thresholds.
     * @param nsteps number of steps of the ODE solver.
     */
    DglapQCD(Grid                                       const& g,
	     function<double(int const&,double const&)> const& InPDFsFunc,
	     double                                     const& MuRef,
	     vector<double>                             const& Masses,
	     vector<double>                             const& Thresholds,
	     int                                        const& PerturbativeOrder,
	     function<double(double const&)>            const& Alphas,
	     double                                     const& IntEps = 1e-5,
	     int                                        const& nsteps = 10);

    DglapQCD(Grid                                       const& g,
	     function<double(int const&,double const&)> const& InPDFsFunc,
	     double                                     const& MuRef,
	     vector<double>                             const& Masses,
	     int                                        const& PerturbativeOrder,
	     function<double(double const&)>            const& Alphas,
	     double                                     const& IntEps = 1e-5,
	     int                                        const& nsteps = 10);

    // The PDF class
    class PDF: public Distribution
    {
      public:
      PDF(Grid                                        const& g,
	  function<double(int const&, double const&)> const& InPDFsFunc,
	  int                                         const& ipdf);
    };

    // Get methods
    Grid                                       const& GetGrid()                 const { return _g; }
    function<double(int const&,double const&)> const& GetInitialDistributions() const { return _InPDFsFunc; }
    double                                     const& GetMuRef()                const { return _MuRef; }
    vector<double>                             const& GetMasses()               const { return _Masses; }
    vector<double>                             const& GetThresholds()           const { return _Thresholds; }
    int                                        const& GetPerturbativeOrder()    const { return _PerturbativeOrder; }
    function<double(double const&)>            const& GetAlphas()               const { return _Alphas; }
    double                                     const& GetIntEps()               const { return _IntEps; }
    int                                        const& GetNumberOfSteps()        const { return _nsteps; }
    Dglap                                      const& GetDglapObject()          const { return *_DglapObj; }

  private:
    Grid                                       const& _g;
    function<double(int const&,double const&)> const& _InPDFsFunc;
    double                                     const& _MuRef;
    vector<double>                             const& _Masses;
    vector<double>                             const& _Thresholds;
    int                                        const& _PerturbativeOrder;
    function<double(double const&)>            const& _Alphas;
    double                                     const& _IntEps;
    int                                               _nsteps;
    shared_ptr<Dglap>                                 _DglapObj;
  };

}
