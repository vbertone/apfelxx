/*
 * APFEL++ 2017
 *
 * Authors: Valerio Bertone: valerio.bertone@vu.nl
 *          Stefano Carrazza: stefano.carrazza@cern.ch
 */

#pragma once

#include <iostream>
#include <exception>
#include <stdexcept>
#include <string>
#include <sstream>

namespace apfel {

  /**
    * Some general tools for the output and exception traitment.
    */

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
