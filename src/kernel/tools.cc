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
  int NF(double const& Q, std::vector<double> const& Thresholds)
  {
    // Compute number of active flavours the the PDF initial scale
    int nf = 0;
    for (auto const& v : Thresholds)
      if (Q > v)
	nf++;
      else
	break;
    return nf;
  }

  //_________________________________________________________________________
  double beta0(int const& nf)
  {
    const double coeff = ( 11 * CA - 4 * TR * nf ) / 3;
    return coeff;
  }

  //_________________________________________________________________________
  double beta1(int const& nf)
  {
    const double coeff = 34 * CA * CA / 3 - 20 * CA * TR * nf / 3 - 4 * CF * TR * nf;
    return coeff;
  }

  //_________________________________________________________________________
  double beta2(int const& nf)
  {
    const double coeff = 2857 * pow(CA,3) / 54
      + ( 2 * CF * CF - 205 * CF * CA / 9 - 1415 * CA * CA / 27 ) * TR * nf
      + ( 44 * CF / 9 + 158 * CA / 27 ) * TR * TR * nf * nf;
    return coeff;
  }

  //_________________________________________________________________________
  double GammaCusp0()
  {
    return 1;
  }

  //_________________________________________________________________________
  double GammaCusp1(int const& nf)
  {
    const double coeff = ( 67. / 9. - Pi2 / 3 ) * CA - 20 * TR * nf / 3;
    return coeff;
  }

  //_________________________________________________________________________
  double GammaCusp2(int const& nf)
  {
    const double coeff = ( 245. / 6. - 134 * Pi2 / 27 + 11 * Pi2 * Pi2 / 45 + 22 * zeta3 / 3 ) * CA * CA
      + ( - 418. / 27. + 40 * Pi2 / 27 - 56 * zeta3 / 3 ) * CA * TR * nf
      + ( - 55. / 3. + 16 * zeta3 ) * CF * TR * nf
      - 16 * TR * TR * nf * nf / 27;
    return coeff;
  }

  //_________________________________________________________________________
  double gammaVq0()
  {
    return - 6 * CF;
  }

  //_________________________________________________________________________
  double gammaVq1(int const& nf)
  {
    const double coeff = CF * CF * ( - 3 + 4 * Pi2 - 48 * zeta3 )
      + CF * CA * ( - 961. / 27. - 11 * Pi2 / 3 + 52 * zeta3 )
      + CF * TR * nf * (  260. / 27. + 4 * Pi2 / 3 );
    return coeff;
  }

  //_________________________________________________________________________
  double gammaVq2(int const& nf)
  {
    const double coeff = pow(CF,3) * (  - 29 - 6 * Pi2 -  16 * Pi2 * Pi2 / 5 - 136 * zeta3 +  32 * Pi2 * zeta3 / 3  + 480 * zeta5 ) 
      + CF * CF * CA * ( - 151. / 2. +  410 * Pi2 / 9 + 494 * Pi2 * Pi2 / 135 - 1688 * zeta3 / 3 - 16 * Pi2 * zeta3 / 3 - 240 * zeta5 ) 
      + CF * CA * CA * ( - 139345. / 1458. - 7163 * Pi2 / 243 - 83 * Pi2 * Pi2 / 45 + 7052 * zeta3 / 9 - 88 * Pi2 * zeta3 / 9 - 272 * zeta5 ) 
      + CF * CF * TR * nf * ( 5906. / 27. - 52 * Pi2 / 9 - 56 * Pi2 * Pi2 / 27 + 1024 * zeta3 / 9 ) 
      + CF * CA * TR * nf * ( - 34636. / 729. + 5188 * Pi2 / 243 + 44 * Pi2 * Pi2 / 45 - 3856 * zeta3 / 27 ) 
      + CF * TR * TR * nf * nf * ( 19336. / 729. - 80 * Pi2 / 27 - 64 * zeta3 / 27 );
    return coeff;
  }

  //_________________________________________________________________________
  double gammaVg0(int const& nf)
  {
    const double coeff = - 22 * CA / 3 + 8 * TR * nf / 3;
    return coeff;
  }

  //_________________________________________________________________________
  double gammaVg1(int const& nf)
  {
    const double coeff = CA * CA * ( - 1384. / 27. + 11 * Pi2 / 9 + 4 * zeta3 )
      + CA * TR * nf * ( 512. / 27. - 4 * Pi2 / 9 ) + 8 * CF * TR * nf;
    return coeff;
  }

  //_________________________________________________________________________
  double gammaVg2(int const& nf)
  {
    const double coeff = 2 * pow(CA,3) * (  - 97186. / 729. + 6109 * Pi2 / 486 - 319  * Pi2 * Pi2 / 270  + 122 * zeta3 / 3 - 20 * Pi2 * zeta3 / 9 - 16 * zeta5 )
      + 2 * CA * CA * TR * nf * ( 30715. / 729. - 1198 * Pi2 / 243 + 82 * Pi2 * Pi2 / 135 + 712 * zeta3 / 27 )
      + 2 * CA * CF * TR * nf * ( 2434. / 27. - 2 * Pi2 / 3 - 8 * Pi2 * Pi2 / 45  - 304 * zeta3 / 9 ) 
      - 4 * CF * CF * TR * nf
      + 2 * CA * TR * TR * nf * nf * ( - 538. / 729. + 40 * Pi2 / 81 - 224 * zeta3 / 27 ) 
      - 88 * CF *  TR * TR * nf * nf / 9;
    return coeff;
  };

  //_________________________________________________________________________
  double CSd10()
  {
    return 0;
  }

  //_________________________________________________________________________
  double CSd11()
  {
    return 2;
  }

  //_________________________________________________________________________
  double CSd20(int const& nf)
  {
    const double coeff = CA * ( 404. / 27. - 14 * zeta3 ) - 112 * TR * nf / 27;
    return coeff;
  }

  //_________________________________________________________________________
  double CSd21(int const& nf)
  {
    return 2 * GammaCusp1(nf);
  }

  //_________________________________________________________________________
  double CSd22(int const& nf)
  {
    return beta0(nf);
  }

  //_________________________________________________________________________
  double CSd30(int const& nf)
  {
    const double coeff = - CA * CA * ( - 176 * zeta3 * zeta2 / 3 + 6392 * zeta2 / 81 + 12328 * zeta3 / 27 + 154 * zeta4 / 3 - 192 * zeta5 - 297029. / 729. ) / 2
      - CA * TR * nf * ( - 824 * zeta2 / 81 - 904 * zeta3 / 27 + 20 * zeta4 / 3 + 62626. / 729. )
      - 2 * TR * TR * nf * nf * ( - 32 * zeta3 / 9 - 1856. / 729. )
      - CF * TR * nf * ( - 304 * zeta3 / 9 - 16 * zeta4 + 1711. / 27. );
    return coeff;
  }

  //_________________________________________________________________________
  double CSd31(int const& nf)
  {
    return 2 * beta0(nf) * CSd20(nf) + 2 * GammaCusp2(nf);
  }

  //_________________________________________________________________________
  double CSd32(int const& nf)
  {
    return 2 * GammaCusp1(nf) * beta0(nf) + beta1(nf);
  }

  //_________________________________________________________________________
  double CSd33(int const& nf)
  {
    const double b0 = beta0(nf);
    return 2 * b0 * b0 / 3;
  }

  //_________________________________________________________________________
  double Lzetaq10()
  {
    return - 3. / 2.;
  }

  //_________________________________________________________________________
  double Lzetaq11()
  {
    return 1. / 2.;
  }

  //_________________________________________________________________________
  double Lzetaq20(int const& nf)
  {
    const double coeff = CF * ( - 3. / 4. + Pi2 - 12 * zeta3 )
      + CA * ( 649. / 108. - 17 * Pi2 / 12 + 19 * zeta3 / 2 )
      + TR * nf * ( - 53. / 27. + Pi2 / 3 );
    return coeff;
  }

  //_________________________________________________________________________
  double Lzetaq22(int const& nf)
  {
    const double coeff = ( 11 * CA - 4 * TR * nf ) / 36;
    return coeff;
  }

  //_________________________________________________________________________
  double Lzetag10(int const& nf)
  {
    const double coeff = - 11. / 6. + 2 * TR * nf / 3 / CA;
    return coeff;
  }

  //_________________________________________________________________________
  double Lzetag11()
  {
    return 1. / 2.;
  }

  //_________________________________________________________________________
  double Lzetag20(int const& nf)
  {
    const double coeff = CA * ( 247. / 54. - 11 * Pi2 / 36 - 5 * zeta3 / 2 )
      + TR * nf * ( - 16. / 3. + Pi2 / 9 )
      + TR * nf * ( 2 * CF + 40 * TR * nf / 27 ) / CA;
    return coeff;
  }

  //_________________________________________________________________________
  double Lzetag22(int const& nf)
  {
    const double coeff = ( 11 * CA - 4 * TR * nf ) / 36;
    return coeff;
  }

  //_________________________________________________________________________
  void info(std::string const& tag, std::string const& what)
  {
    std::cout << code::blue << "[" << tag << "] info: " << what << code::normal << "\n";
  }

  //_________________________________________________________________________
  void warning(std::string const& tag, std::string const& what)
  {
    std::cout << code::yellow << "[" << tag << "] warning: " << what << code::normal << "\n";
  }

  //_________________________________________________________________________
  void success(std::string const& tag, std::string const& what)
  {
    std::cout << code::green << "[" << tag << "] success: " << what << code::normal << "\n";
  }

  //_________________________________________________________________________
  std::string error(std::string const& tag, std::string const& what)
  {
    std::stringstream ss("");
    ss << code::red << "[" << tag << "] error: " << what << code::normal << "\n";
    return ss.str();
  }
}

namespace std {
  //_________________________________________________________________________
  ostream& operator<<(std::ostream& os, apfel::code code)
  {
    return os << "\033[" << static_cast<int>(code) << "m";
  }
}
