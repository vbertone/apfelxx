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

  //! enumerator for code warning
  enum code { red = 31, green = 32, yellow = 33, blue = 34, normal = 39};

  /**
   * @brief info
   * @param tag
   * @param what
   */
  inline void info(std::string const& tag, std::string const &what)
  {
    std::cout << code::blue << "[" << tag << "] info: " << what << code::normal << "\n";
  }

  /**
   * @brief warning
   * @param tag
   * @param what
   */
  inline void warning(std::string const& tag, std::string const &what)
  {
    std::cout << code::yellow << "[" << tag << "] warning: " << what << code::normal << "\n";
  }

  /**
   * @brief success
   * @param tag
   * @param what
   */
  inline void success(std::string const& tag, std::string const &what)
  {
    std::cout << code::green << "[" << tag << "] success: " << what << code::normal << "\n";
  }

  /**
   * @brief error
   * @param tag
   * @param what
   * @return
   */
  inline std::string error(std::string const& tag, std::string const &what)
  {
    std::stringstream ss("");
    ss << code::red << "[" << tag << "] error: " << what << code::normal << "\n";
    return ss.str();
  }

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
