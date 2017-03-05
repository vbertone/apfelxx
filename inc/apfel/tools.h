//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <iostream>
#include <exception>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

namespace apfel
{
  /**
    * Some general tools for the output and exception traitment.
    */
  const double eps2  = 1e-2;
  const double eps3  = 1e-3;
  const double eps4  = 1e-4;
  const double eps5  = 1e-5;
  const double eps6  = 1e-6;
  const double eps7  = 1e-7;
  const double eps8  = 1e-8;
  const double eps9  = 1e-9;
  const double eps10 = 1e-10;
  const double eps11 = 1e-11;
  const double eps12 = 1e-12;
  const double eps13 = 1e-13;
  const double eps14 = 1e-14;
  const double eps15 = 1e-15;
  const double eps25 = 1e-25;

  // Color factors
  const double TR = 0.5;
  const double CF = 4./3.;
  const double CA = 3.;

  // Useful constants
  const double FourPi = 4 * M_PI;
  const double Pi2    = M_PI * M_PI;
  const double emc    = 0.5772156649;
  const double zeta2  = 1.644934067;
  const double zeta3  = 1.2020569031;
  const double zeta4  = 1.0823232337;
  const double zeta5  = 1.0369277551;

  //! enumerator for code warning
  enum code { red = 31, green = 32, yellow = 33, blue = 34, normal = 39};

  /**
   * @brief info
   * @param tag
   * @param what
   */
  void info(std::string const& tag, std::string const &what);

  /**
   * @brief warning
   * @param tag
   * @param what
   */
  void warning(std::string const& tag, std::string const &what);

  /**
   * @brief success
   * @param tag
   * @param what
   */
  void success(std::string const& tag, std::string const &what);

  /**
   * @brief error
   * @param tag
   * @param what
   * @return
   */
  std::string error(std::string const& tag, std::string const &what);

  /**
   * @brief The runtime_exception class
   */
  class runtime_exception: public std::runtime_error
  {
  public:
    runtime_exception(const std::string &tag,
                      const std::string &what):
      std::runtime_error(error(tag,what)) {}
  };

  /**
   * @brief The logic_exception class
   */
  class logic_exception: public std::logic_error
  {
  public:
    logic_exception(const std::string &tag,
                    const std::string &what):
      std::logic_error(error(tag,what)) {}
  };
}

namespace std {
  /**
   * @brief implementation of operator<< for apfel::code
   */
  std::ostream& operator<<(std::ostream& os, apfel::code code);
}
