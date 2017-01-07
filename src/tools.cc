//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/tools.h"

namespace apfel
{
  //_________________________________________________________________________
  void info(std::string const& tag, std::string const &what)
  {
    std::cout << code::blue << "[" << tag << "] info: " << what << code::normal << "\n";
  }

  //_________________________________________________________________________
  void warning(std::string const& tag, std::string const &what)
  {
    std::cout << code::yellow << "[" << tag << "] warning: " << what << code::normal << "\n";
  }

  //_________________________________________________________________________
  void success(std::string const& tag, std::string const &what)
  {
    std::cout << code::green << "[" << tag << "] success: " << what << code::normal << "\n";
  }

  //_________________________________________________________________________
  std::string error(std::string const& tag, std::string const &what)
  {
    std::stringstream ss("");
    ss << code::red << "[" << tag << "] error: " << what << code::normal << "\n";
    return ss.str();
  }
}

namespace std {

  ostream& operator<<(std::ostream& os, apfel::code code)
  {
    return os << "\033[" << static_cast<int>(code) << "m";
  }

}
