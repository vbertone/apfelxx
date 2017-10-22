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
   * @brief Verbosity enumerator
   *
   * LOW    = No informative and warning messages are displayed. Error messages are printed anyway.
   * MEDIUM = No informative messages are displayed. Warning and error messages are printed.
   * HIGH   = All messages are displayed.
   */
  enum verbosity: int {LOW, MEDIUM, HIGH};

  /**
   * @brief Current Verbosity level
   */
  static int VerbosityLevel = MEDIUM;

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
}
