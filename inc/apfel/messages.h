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

/**
 * @brief Namespace for all APFEL++ functions and classes
 */
namespace apfel
{
  /**
   * @name Message functions
   * Collection of functions related to the verbosity of the code.
   */
  ///@{
  /**
   * @brief Set Verbosity level
   * @param vl: verbosity level
   * @note possible values of vl:
   * LOW    = No informative and warning messages are displayed. Error messages are printed anyway
   * MEDIUM = No informative messages are displayed. Warning and error messages are printed
   * HIGH   = All messages are displayed
   */
  void SetVerbosityLevel(int const& vl);

  /**
   * @brief Get Verbosity level
   * @return the current verbosity level
   */
  int GetVerbosityLevel();

  /**
   * @brief Function that prints information on screen. Effective
   * according to the verbosity level.
   * @param what: the message to report
   */
  void report(string const& what);

  /**
   * @brief Function that prints information on screen. Effective
   * according to the verbosity level.
   * @param what: the message to report
   */
  void info(string const& tag, string const& what);

  /**
   * @brief Function that prints warnings on screen. Effective
   * according to the verbosity level.
   * @param what: the warning to report
   */
  void warning(string const& tag, string const& what);

  /**
   * @brief Function that prints information on screen. Always
   * effective.
   * @param tag: the function that generates the error
   * @param what: the error to report
   * @return the error message
   */
  string error(string const& tag, string const& what);

  /**
   * @brief Function that prints the APFEL++ banner on
   * screen. Effective according to the verbosity level.
   */
  void Banner();
  ///@}
}
