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

  // Enumerator for code warning
  enum code {red = 31, green = 32, yellow = 33, blue = 34, normal = 39};

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
