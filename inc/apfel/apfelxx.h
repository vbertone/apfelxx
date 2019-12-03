//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

/**
 * @name Main include of APFEL++
 */
///@{
/**
 * @brief Utilities
 */
#include "apfel/constants.h"
#include "apfel/lhtoypdfs.h"
#include "apfel/messages.h"
#include "apfel/rotations.h"
#include "apfel/specialfunctions.h"
#include "apfel/timer.h"
#include "apfel/tools.h"
#include "apfel/version.h"

/**
 * @brief Basic objects
 */
#include "apfel/distribution.h"
#include "apfel/doubleobject.h"
#include "apfel/expression.h"
#include "apfel/grid.h"
#include "apfel/integrator.h"
#include "apfel/matchedevolution.h"
#include "apfel/observable.h"
#include "apfel/ogataquadrature.h"
#include "apfel/doubleexponentialquadrature.h"
#include "apfel/operator.h"
#include "apfel/qgrid.h"
#include "apfel/set.h"
#include "apfel/tabulateobject.h"
#include "apfel/evolutionbasisqcd.h"

/**
 * @brief Convolution maps
 */
#include "apfel/disbasis.h"
#include "apfel/evolutionbasisqcd.h"
#include "apfel/matchingbasisqcd.h"

/**
 * @brief High-level interface
 */
#include "apfel/evolutionsetup.h"
#include "apfel/initialiseevolution.h"

/**
 * @brief Builders and relevant physical quantities
 */
#include "apfel/dglapbuilder.h"
#include "apfel/structurefunctionbuilder.h"
#include "apfel/tmdbuilder.h"
#include "apfel/alphaqcd.h"
#include "apfel/alphaqed.h"
#include "apfel/twobodyphasespace.h"
///@}
