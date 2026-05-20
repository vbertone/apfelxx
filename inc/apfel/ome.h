/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

/**
 * \file
 * \brief Main header for libome
 *
 * \details
 * This header includes all necessary further headers to use the library.
 * It works for both C and C++. If it is included from C sources, only the
 * C interface will be exposed.
 */

#ifndef LIBOME_OME_H
#define LIBOME_OME_H

#ifdef __cplusplus
#include "apfel/ome/ome_type_aliases.h"
#endif

#include "apfel/ome/AqqQNSEven.h"
#include "apfel/ome/AqqQNSOdd.h"
#include "apfel/ome/AQqPS.h"
#include "apfel/ome/AQqPSs.h"
#include "apfel/ome/AqqQPS.h"
#include "apfel/ome/AqgQ.h"
#include "apfel/ome/AgqQ.h"
#include "apfel/ome/AggQ.h"
#include "apfel/ome/AQg.h"

#include "apfel/ome/polAqqQNSEven.h"
#include "apfel/ome/polAqqQNSOdd.h"
#include "apfel/ome/polAQqPS.h"
#include "apfel/ome/polAqqQPS.h"
#include "apfel/ome/polAqgQ.h"
#include "apfel/ome/polAgqQ.h"
#include "apfel/ome/polAggQ.h"
#include "apfel/ome/polAQg.h"

#endif
