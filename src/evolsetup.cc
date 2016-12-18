/*
  evolsetup.cc:

  Author: Valerio Bertone
*/

#include <stdexcept>

#include "apfel/evolsetup.h"

using namespace std;

namespace apfel {

  // ================================================================================
  // Setters
  void evolsetup::EnableWelcomeMessage(bool const& WelcomeMessage_)
  {
    _WelcomeMessage = WelcomeMessage_;
    return;
  }

  void evolsetup::SetQLimits(double const& Qmin_, double const& Qmax_)
  {
    _Qmin = Qmin_;
    _Qmax = Qmax_;
    return;
  }

  void evolsetup::SetPerturbativeOrder(int const& PerturbativeOrder_)
  {
    _PerturbativeOrder = PerturbativeOrder_;
    return;
  }

  void evolsetup::SetVFNS()
  {
    _FlavourScheme = VFNS;
    return;
  }

  void evolsetup::SetFFNS(int const& Nf_FF_)
  {
    _FlavourScheme = VFNS;
    _Nf_FF = Nf_FF_;
    return;
  }

  void evolsetup::SetTheory(theory const& Theory_)
  {
    _Theory = Theory_;
    return;
  }

  void evolsetup::SetFastEvolution(bool const& FastEvolution_)
  {
    _FastEvolution = FastEvolution_;
    return;
  }

  void evolsetup::SetTimeLikeEvolution(bool const& TimeLikeEvolution_)
  {
    _TimeLikeEvolution = TimeLikeEvolution_;
    return;
  }

  void evolsetup::SetPolarizedEvolution(bool const& PolarizedEvolution_)
  {
    _PolarizedEvolution = PolarizedEvolution_;
    return;
  }

  void evolsetup::SetSmallxResummation(bool const& SmallxResummation_, logaccuracy const& LogAccuracy_)
  {
    _SmallxResummation = SmallxResummation_;
    _LogAccuracy = LogAccuracy_;
    return;
  }

  void evolsetup::SetAlphaQCDRef(double const& AlphaQCDRef_, double const& QQCDRef_)
  {
    _AlphaQCDRef = AlphaQCDRef_;
    _QQCDRef = QQCDRef_;
    return;
  }

  void evolsetup::SetAlphaQEDRef(double const& AlphaQEDRef_, double const& QQEDRef_)
  {
    _AlphaQEDRef = AlphaQEDRef_;
    _QQEDRef = QQEDRef_;
    return;
  }

  void evolsetup::SetLambdaQCDRef(double const& LambdaQCDRef_, int const& nLambdaQCDRef_)
  {
    _LambdaQCDRef = LambdaQCDRef_;
    _nLambdaQCDRef = nLambdaQCDRef_;
    return;
  }

  void evolsetup::SetEpsilonTruncation(double const& EpsilonTruncation_)
  {
    _EpsilonTruncation = EpsilonTruncation_;
    return;
  }

  void evolsetup::SetAlphaEvolution(alphaevolution const& AlphaEvolution_)
  {
    _AlphaEvolution = AlphaEvolution_;
    return;
  }

  void evolsetup::SetPDFEvolution(pdfevolution const& PDFEvolution_)
  {
    _PDFEvolution = PDFEvolution_;
    return;
  }

  void evolsetup::SetRenFacRatio(double const& RenFacRatio_)
  {
    _RenFacRatio = RenFacRatio_;
    return;
  }

  void evolsetup::SetPoleMasses(double const& MCharm_, double const& MBottom_, double const& MTop_)
  {
    _MCharm = MCharm_;
    _MBottom = MBottom_;
    _MTop = MTop_;
    return;
  }

  void evolsetup::SetMSBarMasses(double const& MCharm_, double const& MBottom_, double const& MTop_)
  {
    _MCharm = MCharm_;
    _MBottom = MBottom_;
    _MTop = MTop_;
    return;
  }

  void evolsetup::SetMassMatchingScales(double const& kThCharm_, double const& kThBottom_, double const& kThTop_)
  {
    _kThCharm = kThCharm_;
    _kThBottom = kThBottom_;
    _kThTop = kThTop_;
    return;
  }

  void evolsetup::SetMassScaleReference(double const& QRefMCharm_, double const& QRefMBottom_, double const& QRefMTop_)
  {
    _QRefMCharm = QRefMCharm_;
    _QRefMBottom = QRefMBottom_;
    _QRefMTop = QRefMTop_;
    return;
  }

  void evolsetup::SetTauMass(double const& TauMass_)
  {
    _TauMass = TauMass_;
    return;
  }

  void evolsetup::EnableMassRunning(bool const& MassRunning_)
  {
    _MassRunning = MassRunning_;
    return;
  }

  void evolsetup::SetMaxFlavourPDFs(int const& MaxFlavourPDFs_)
  {
    _MaxFlavourPDFs = MaxFlavourPDFs_;
    return;
  }

  void evolsetup::SetMaxFlavourAlpha(int const& MaxFlavourAlpha_)
  {
    _MaxFlavourAlpha = MaxFlavourAlpha_;
    return;
  }

  void evolsetup::SetPDFSet(string const& PDFSet_)
  {
    _PDFSet = PDFSet_;
    return;
  }

  void evolsetup::SetReplica(int const& Replica_)
  {
    _Replica = Replica_;
    return;
  }

  void evolsetup::EnableEvolutionOperator(bool const& EvolutionOperator_)
  {
    _EvolutionOperator = EvolutionOperator_;
    return;
  }

  void evolsetup::EnableLeptonEvolution(bool const& LeptonEvolution_)
  {
    _LeptonEvolution = LeptonEvolution_;
    return;
  }

  void evolsetup::SetLHgridParameters(int const& nxLHA_, int const& nxmLHA_, double const& xminLHA_, double const& xmLHA_, double const& xmaxLHA_,
				      int const& nQ2LHA_, double const& Q2minLHA_, double const& Q2maxLHA_)
  {
    _nxLHA    = nxLHA_;
    _nxmLHA   = nxmLHA_;
    _xminLHA  = xminLHA_;
    _xmLHA    = xmLHA_;
    _xmaxLHA  = xmaxLHA_;
    _nQ2LHA   = nQ2LHA_;
    _Q2minLHA = Q2minLHA_;
    _Q2maxLHA = Q2maxLHA_;
    return;
  }

  void evolsetup::SetQGridParameters(int const& nQ2g_, int const& InterDegreeQ_)
  {
    _nQ2g = nQ2g_;
    _InterDegreeQ = InterDegreeQ_;
    return;
  }

  void evolsetup::LockGrids(bool const& Locked_)
  {
    _Locked = Locked_;
    return;
  }

  void evolsetup::SetGridParameters(int const& nx_, int const& id_, double const& xmin_)
  {
    // Assign grid parameters
    gridparams gp = {nx_, xmin_, id_, NULL};
    _GridParams.push_back(gp);
    return;
  }

  void evolsetup::SetGridParameters(int const& nx_, int const& id_, double *xgext_)
  {
    // Assign grid parameters
    gridparams gp = {nx_, xgext_[0], id_, xgext_};
    _GridParams.push_back(gp);
    return;
  }

  void evolsetup::SetGaussPoints(int const& GaussPoints_)
  {
    _GaussPoints = GaussPoints_;
    return;
  }

  // ================================================================================
  // Getters
  int evolsetup::nx(int const& ig_) const
  {
    if(ig_ < 0 || ig_ > _GridParams.size() ) throw runtime_error("Index out of range.");
    return _GridParams[ig_].nx;
  }

  int evolsetup::InterDegree(int const& ig_) const
  {
    if(ig_ < 0 || ig_ > _GridParams.size() ) throw runtime_error("Index out of range.");
    return _GridParams[ig_].id;
  }

  double evolsetup::xMin(int const& ig_) const
  {
    if(ig_ < 0 || ig_ > _GridParams.size() ) throw runtime_error("Index out of range.");
    return _GridParams[ig_].xmin;
  }

  int evolsetup::NfFinPDF() const
  {
    int nf;
    if(_FlavourScheme == FFNS)               nf = _Nf_FF;
    else if(_FlavourScheme == VFNS) {
      if     (_Qmax > _kThTop * _MTop)       nf = 6;
      else if(_Qmax > _kThBottom * _MBottom) nf = 5;
      else if(_Qmax > _kThCharm * _MCharm)   nf = 4;
      else                                   nf = 3;
      // Limit to maximum number of flavours allowed
      if(nf > _MaxFlavourPDFs)               nf = _MaxFlavourPDFs;
    }
    return nf;
  }

  int evolsetup::NfIniPDF() const
  {
    int nf;
    if(_FlavourScheme == FFNS)               nf = _Nf_FF;
    else if(_FlavourScheme == VFNS) {
      if     (_Qmin > _kThTop * _MTop)       nf = 6;
      else if(_Qmin > _kThBottom * _MBottom) nf = 5;
      else if(_Qmin > _kThCharm * _MCharm)   nf = 4;
      else                                   nf = 3;
      // Limit to maximum number of flavours allowed
      if(nf > _MaxFlavourPDFs)               nf = _MaxFlavourPDFs;
    }
    return nf;
  }

  int evolsetup::NfFinAlpha() const
  {
    int nf;
    if(_FlavourScheme == FFNS)                              nf = _Nf_FF;
    else if(_FlavourScheme == VFNS) {
      if     (_Qmax > _RenFacRatio * _kThTop * _MTop)       nf = 6;
      else if(_Qmax > _RenFacRatio * _kThBottom * _MBottom) nf = 5;
      else if(_Qmax > _RenFacRatio * _kThCharm * _MCharm)   nf = 4;
      else                                                  nf = 3;
      // Limit to maximum number of flavours allowed
      if(nf > _MaxFlavourPDFs)                              nf = _MaxFlavourAlpha;
    }
    return nf;
  }

  int evolsetup::NfIniAlpha() const
  {
    int nf;
    if(_FlavourScheme == FFNS)                              nf = _Nf_FF;
    else if(_FlavourScheme == VFNS) {
      if     (_Qmin > _RenFacRatio * _kThTop * _MTop)       nf = 6;
      else if(_Qmin > _RenFacRatio * _kThBottom * _MBottom) nf = 5;
      else if(_Qmin > _RenFacRatio * _kThCharm * _MCharm)   nf = 4;
      else                                                  nf = 3;
      // Limit to maximum number of flavours allowed
      if(nf > _MaxFlavourPDFs)                              nf = _MaxFlavourAlpha;
    }
    return nf;
  }

}
