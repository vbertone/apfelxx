//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <sstream>
#include <iostream>

using namespace std;

namespace apfel
{
  /**
   * @brief Set Verbosity level
   */
  void SetVerbosityLevel(int const& vl);

  /**
   * @brief Get Verbosity level
   */
  int GetVerbosityLevel();

  /**
   * @brief info
   * @param tag
   * @param what
   */
  void report(string const& what);

  /**
   * @brief info
   * @param tag
   * @param what
   */
  void info(string const& tag, string const& what);

  /**
   * @brief warning
   * @param tag
   * @param what
   */
  void warning(string const& tag, string const& what);

  /**
   * @brief error
   * @param tag
   * @param what
   * @return
   */
  string error(string const& tag, string const& what);

  /**
   * @brief Banner
   */
  void Banner();
}
