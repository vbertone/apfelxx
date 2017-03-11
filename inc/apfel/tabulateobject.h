//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/qgrid.h"
#include "apfel/matchedevolution.h"

namespace apfel
{
  /**
   * @brief
   */
  template<class T>
  class TabulateObject: public QGrid<T>
  {
  public:
    /**
     * @brief TabulateObject default constructor.
     */
    TabulateObject(MatchedEvolution<T> *Object,
		   int                 const& nQ,
		   double              const& QMin,
		   double              const& QMax,
		   int                 const& InterDegree);

  };

}
