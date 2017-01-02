//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <array>
using std::array;

namespace apfel
{
  /**
   * @brief The Distribution class for PDFs.
   *
   * This class provides methods to inherit a custom PDF distribution
   * in a generic basis.
   */
  class Distribution
  {
  public:
    /**
     * LHA-style flavour basis
     */
    enum lhaBasis {TBAR,BBAR,CBAR,SBAR,UBAR,DBAR,GLUON,D,U,S,C,B,T,PHT};

    /**
     * EVLN basis
     */
    enum evolBasis {EVLN_GAM, EVLN_SNG, EVLN_GLU, EVLN_VAL, EVLN_V3,
                    EVLN_V8, EVLN_V15, EVLN_V24, EVLN_V35,  EVLN_T3,
                    EVLN_T8, EVLN_T15, EVLN_T24, EVLN_T35 };
  protected:
    /**
     * @brief Distribution constructors.
     */
    Distribution();

  public:
    /**
     * @brief GetPDF return the PDF value.
     *
     * This method must be implemented in the inherited class.
     *
     * @param id the flavor id number following the pdg-id convention.
     * @param x the momentum fraction
     * @param q the energy scale
     * @param n the replica number
     * @return the evaluated PDF
     */
    virtual double GetPDF(int const& id, double const& x, double const& q, int const& n = 0) const = 0;

    /**
     * @brief Virtual method which rotates from custom basis to evolution basis.
     *
     * @param[in] basis custom basis 14 entries, suggested to list the basis in a enumerator.
     * @param[out] evln array with 14 entries following the \c evolBasis enumerator.
     */
    virtual void BASIS2EVLN(array<double,14> const& basis, array<double,14> &evln) const;

    /**
     * @brief Virtual method which rotates from evolution basis to custom basis.
     *
     * @param[in] evln array with 14 entries following the \c evolBasis enumerator.
     * @param[out] basis custom basis 14 entries, suggested to list the basis in a enumerator.
     */
    virtual void EVLN2BASIS(array<double,14> const& evln, array<double,14> &basis) const;

    /**
     * @brief Static method for rotating from the LHA basis to the Evolution basis.
     *
     * @param[in] LHA array with 14 entries using the enumerator order defined in \c lhaBasis
     * @param[out] EVLN array with 14 entries following the \c evolBasis enumerator.
     */
    static void LHA2EVLN(array<double,14> const& LHA, array<double,14> &EVLN);

    /**
     * @brief Static method for rotating from the Evolution basis to the LHA basis.
     *
     * @param[in] EVLN array with 14 entries following the \c evolBasis enumerator.
     * @param[out] LHA array with 14 entries using the enumerator order defined in \c lhaBasis
     */
    static void EVLN2LHA(array<double,14> const& EVLN, array<double,14> &LHA);
  };

}
