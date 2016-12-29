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

  /**
    * Parameters of the dgauss integration
    */
  const double cst = 1e-25;
  const vector<vector<double> > gq_x = {
    {0.3399810435848562,
     0.8611363115940525},
    {0.1834346424956498,
     0.5255324099163289,
     0.7966664774136267,
     0.9602898564975362},
    {0.0950125098376374,
     0.2816035507792589,
     0.4580167776572273,
     0.6178762444026437,
     0.7554044083550030,
     0.8656312023878317,
     0.9445750230732325,
     0.9894009349916499},
    {0.0483076656877383,
     0.1444719615827964,
     0.2392873622521370,
     0.3318686022821276,
     0.4213512761306353,
     0.5068999089322293,
     0.5877157572407623,
     0.6630442669302152,
     0.7321821187402896,
     0.7944837959679424,
     0.8493676137325699,
     0.8963211557660521,
     0.9349060759377396,
     0.9647622555875064,
     0.9856115115452683,
     0.9972638618494815},
    {0.0243502926634244,
     0.0729931217877990,
     0.1214628192961205,
     0.1696444204239928,
     0.2174236437400070,
     0.2646871622087674,
     0.3113228719902109,
     0.3572201583376681,
     0.4022701579639916,
     0.4463660172534640,
     0.4894031457070529,
     0.5312794640198945,
     0.5718956462026340,
     0.6111553551723932,
     0.6489654712546573,
     0.6852363130542332,
     0.7198818501716108,
     0.7528199072605318,
     0.7839723589433414,
     0.8132653151227975,
     0.8406292962525803,
     0.8659993981540928,
     0.8893154459951141,
     0.9105221370785028,
     0.9295691721319395,
     0.9464113748584028,
     0.9610087996520537,
     0.9733268277899109,
     0.9833362538846259,
     0.9910133714767443,
     0.9963401167719552,
     0.9993050417357721}};
  const vector<vector<double> > gq_w = {
  // 4-point integration
    {0.6521451548625461,
     0.3478548451374538},
  // 8-point integration
    {0.3626837833783619,
     0.3137066458778872,
     0.2223810344533744,
     0.1012285362903762},
  // 16-point integration
    {0.1894506104550684,
     0.1826034150449235,
     0.1691565193950025,
     0.1495959888165767,
     0.1246289712555338,
     0.0951585116824927,
     0.0622535239386478,
     0.0271524594117540},
  // 32-point integration
    {0.0965400885147278,
     0.0956387200792748,
     0.0938443990808045,
     0.0911738786957638,
     0.0876520930044038,
     0.0833119242269467,
     0.0781938957870703,
     0.0723457941088485,
     0.0658222227763618,
     0.0586840934785355,
     0.0509980592623761,
     0.0428358980222266,
     0.0342738629130214,
     0.0253920653092620,
     0.0162743947309056,
     0.0070186100094700},
  // 64-point integration
    {0.0486909570091397,
     0.0485754674415034,
     0.0483447622348029,
     0.0479993885964583,
     0.0475401657148303,
     0.0469681828162100,
     0.0462847965813144,
     0.0454916279274181,
     0.0445905581637565,
     0.0435837245293234,
     0.0424735151236535,
     0.0412625632426235,
     0.0399537411327203,
     0.0385501531786156,
     0.0370551285402400,
     0.0354722132568823,
     0.0338051618371416,
     0.0320579283548515,
     0.0302346570724024,
     0.0283396726142594,
     0.0263774697150546,
     0.0243527025687108,
     0.0222701738083832,
     0.0201348231535302,
     0.0179517157756973,
     0.0157260304760247,
     0.0134630478967186,
     0.0111681394601311,
     0.0088467598263639,
     0.0065044579689783,
     0.0041470332605624,
     0.0017832807216964}};

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
