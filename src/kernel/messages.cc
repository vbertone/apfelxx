//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/messages.h"
#include "apfel/version.h"

namespace apfel
{
  //_________________________________________________________________________
  void SetVerbosityLevel(int const& vl)
  {
    VerbosityLevel = vl;
  }

  //_________________________________________________________________________
  int GetVerbosityLevel()
  {
    return VerbosityLevel;
  }

  //_________________________________________________________________________
  void report(string const& what)
  {
    if (VerbosityLevel > MEDIUM)
      cout << what;
  }

  //_________________________________________________________________________
  void info(string const& tag, string const& what)
  {
    if (VerbosityLevel > MEDIUM)
      cout << "\033[" << code::blue << "m[" << tag << "] Info: " << what << "\033[" << code::normal << "m\n";
  }

  //_________________________________________________________________________
  void warning(string const& tag, string const& what)
  {
    if (VerbosityLevel > LOW)
      cout << "\033[" << code::yellow << "m[" << tag << "] Warning: " << what << "\033[" << code::normal << "m\n";
  }

  //_________________________________________________________________________
  string error(string const& tag, string const& what)
  {
    stringstream ss;
    ss << "\033[" << code::red << "m[" << tag << "] Error: " << what << "\033[" << code::normal << "m\n";
    return ss.str();
  }

  //_________________________________________________________________________
  void Banner()
  {
    if (VerbosityLevel > MEDIUM)
      {
	cout << "\033[" << code::green << "m\n";
	cout << "     _/_/_/   _/_/_/_/  _/_/_/_/  _/_/_/_/  _/\n";
	cout << "   _/    _/  _/    _/  _/        _/        _/        _/     _/\n";
	cout << "  _/_/_/_/  _/_/_/_/  _/_/_/    _/_/_/    _/      _/_/_/ _/_/_/\n";
	cout << " _/    _/  _/        _/        _/        _/        _/     _/\n";
	cout << "_/    _/  _/        _/        _/_/_/_/  _/_/_/_/\n\n";
	cout << "_____v" << VERSION << ": A new PDF evolution library in C++, arXiv:1708.00911\n";
	cout << "     Authors: V. Bertone, S. Carrazza\n";
	cout << "\033[" << code::normal << "m\n";
      }
  }
}
