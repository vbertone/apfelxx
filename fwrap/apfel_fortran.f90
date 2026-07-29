!
! apfelxx Fortran wrapper
!
! Public, old-APFEL-style procedure names on top of bind(C) interfaces
! into the apfelfortran.cc shim. Keeping the raw bind(C) interfaces
! private and re-exporting them under the classic names (CheckAPFEL,
! InitializeAPFEL, EvolveAPFEL, xPDF, ...) lets the C symbol names and
! the user-facing Fortran names differ, which avoids clashing with the
! genuine Fortran symbols of the old APFEL library if both are ever
! linked into the same program.
!
module apfel_fortran
  use iso_c_binding
  implicit none

  private
  public :: CheckAPFEL
  public :: SetPerturbativeOrder
  public :: SetPoleMasses
  public :: SetAlphaQCDRef
  public :: SetQLimits
  public :: SetNumberOfGrids
  public :: SetGridParameters
  public :: SetQGridParameters
  public :: InitializeAPFEL
  public :: alphaQCD
  public :: GetPerturbativeOrder
  public :: SetPDFSet
  public :: EvolveAPFEL
  public :: CachePDFsAPFEL
  public :: xPDF
  public :: xPDFxQ
  public :: SetRenFacRatio
  public :: SetReplica

  interface
    function apfelxxf_checkapfel() bind(C, name="apfelxxf_checkapfel") result(res)
      import :: c_bool
      logical(c_bool) :: res
    end function apfelxxf_checkapfel

    subroutine apfelxxf_setperturbativeorder(pto) bind(C, name="apfelxxf_setperturbativeorder")
      import :: c_int
      integer(c_int), value :: pto
    end subroutine apfelxxf_setperturbativeorder

    subroutine apfelxxf_setpolemasses(mc, mb, mt) bind(C, name="apfelxxf_setpolemasses")
      import :: c_double
      real(c_double), value :: mc, mb, mt
    end subroutine apfelxxf_setpolemasses

    subroutine apfelxxf_setalphaqcdref(alpharef, qref) bind(C, name="apfelxxf_setalphaqcdref")
      import :: c_double
      real(c_double), value :: alpharef, qref
    end subroutine apfelxxf_setalphaqcdref

    subroutine apfelxxf_setqlimits(qmin, qmax) bind(C, name="apfelxxf_setqlimits")
      import :: c_double
      real(c_double), value :: qmin, qmax
    end subroutine apfelxxf_setqlimits

    subroutine apfelxxf_setnumberofgrids(n) bind(C, name="apfelxxf_setnumberofgrids")
      import :: c_int
      integer(c_int), value :: n
    end subroutine apfelxxf_setnumberofgrids

    subroutine apfelxxf_setgridparameters(i, np, deg, x) bind(C, name="apfelxxf_setgridparameters")
      import :: c_int, c_double
      integer(c_int), value :: i, np, deg
      real(c_double), value :: x
    end subroutine apfelxxf_setgridparameters

    subroutine apfelxxf_setqgridparameters(npq, degq) bind(C, name="apfelxxf_setqgridparameters")
      import :: c_int
      integer(c_int), value :: npq, degq
    end subroutine apfelxxf_setqgridparameters

    subroutine apfelxxf_initializeapfel() bind(C, name="apfelxxf_initializeapfel")
    end subroutine apfelxxf_initializeapfel

    function apfelxxf_alphaqcd(q) bind(C, name="apfelxxf_alphaqcd") result(res)
      import :: c_double
      real(c_double), value :: q
      real(c_double) :: res
    end function apfelxxf_alphaqcd

    function apfelxxf_getperturbativeorder() bind(C, name="apfelxxf_getperturbativeorder") result(res)
      import :: c_int
      integer(c_int) :: res
    end function apfelxxf_getperturbativeorder

    subroutine apfelxxf_setpdfset(name, length) bind(C, name="apfelxxf_setpdfset")
      import :: c_char, c_int
      character(kind=c_char), dimension(*), intent(in) :: name
      integer(c_int), value :: length
    end subroutine apfelxxf_setpdfset

    subroutine apfelxxf_evolveapfel(q0, q) bind(C, name="apfelxxf_evolveapfel")
      import :: c_double
      real(c_double), value :: q0, q
    end subroutine apfelxxf_evolveapfel

    subroutine apfelxxf_cachepdfsapfel(q0) bind(C, name="apfelxxf_cachepdfsapfel")
      import :: c_double
      real(c_double), value :: q0
    end subroutine apfelxxf_cachepdfsapfel

    function apfelxxf_xpdf(i, x) bind(C, name="apfelxxf_xpdf") result(res)
      import :: c_int, c_double
      integer(c_int), value :: i
      real(c_double), value :: x
      real(c_double) :: res
    end function apfelxxf_xpdf

    function apfelxxf_xpdfxq(i, x, q) bind(C, name="apfelxxf_xpdfxq") result(res)
      import :: c_int, c_double
      integer(c_int), value :: i
      real(c_double), value :: x, q
      real(c_double) :: res
    end function apfelxxf_xpdfxq

    subroutine apfelxxf_setrenfacratio(xi) bind(C, name="apfelxxf_setrenfacratio")
      import :: c_double
      real(c_double), value :: xi
    end subroutine apfelxxf_setrenfacratio

    subroutine apfelxxf_setreplica(irep) bind(C, name="apfelxxf_setreplica")
      import :: c_int
      integer(c_int), value :: irep
    end subroutine apfelxxf_setreplica
  end interface

contains

  logical function CheckAPFEL()
    CheckAPFEL = logical(apfelxxf_checkapfel())
  end function CheckAPFEL

  subroutine SetPerturbativeOrder(pto)
    integer, intent(in) :: pto
    call apfelxxf_setperturbativeorder(int(pto, c_int))
  end subroutine SetPerturbativeOrder

  subroutine SetPoleMasses(mc, mb, mt)
    double precision, intent(in) :: mc, mb, mt
    call apfelxxf_setpolemasses(real(mc, c_double), real(mb, c_double), real(mt, c_double))
  end subroutine SetPoleMasses

  subroutine SetAlphaQCDRef(alpharef, qref)
    double precision, intent(in) :: alpharef, qref
    call apfelxxf_setalphaqcdref(real(alpharef, c_double), real(qref, c_double))
  end subroutine SetAlphaQCDRef

  subroutine SetQLimits(qmin, qmax)
    double precision, intent(in) :: qmin, qmax
    call apfelxxf_setqlimits(real(qmin, c_double), real(qmax, c_double))
  end subroutine SetQLimits

  subroutine SetNumberOfGrids(n)
    integer, intent(in) :: n
    call apfelxxf_setnumberofgrids(int(n, c_int))
  end subroutine SetNumberOfGrids

  subroutine SetGridParameters(i, np, deg, x)
    integer, intent(in) :: i, np, deg
    double precision, intent(in) :: x
    call apfelxxf_setgridparameters(int(i, c_int), int(np, c_int), int(deg, c_int), real(x, c_double))
  end subroutine SetGridParameters

  subroutine SetQGridParameters(npq, degq)
    integer, intent(in) :: npq, degq
    call apfelxxf_setqgridparameters(int(npq, c_int), int(degq, c_int))
  end subroutine SetQGridParameters

  subroutine InitializeAPFEL()
    call apfelxxf_initializeapfel()
  end subroutine InitializeAPFEL

  double precision function alphaQCD(q)
    double precision, intent(in) :: q
    alphaQCD = apfelxxf_alphaqcd(real(q, c_double))
  end function alphaQCD

  integer function GetPerturbativeOrder()
    GetPerturbativeOrder = apfelxxf_getperturbativeorder()
  end function GetPerturbativeOrder

  subroutine SetPDFSet(name)
    character(len=*), intent(in) :: name
    call apfelxxf_setpdfset(name, len(name, kind=c_int))
  end subroutine SetPDFSet

  subroutine EvolveAPFEL(q0, q)
    double precision, intent(in) :: q0, q
    call apfelxxf_evolveapfel(real(q0, c_double), real(q, c_double))
  end subroutine EvolveAPFEL

  subroutine CachePDFsAPFEL(q0)
    double precision, intent(in) :: q0
    call apfelxxf_cachepdfsapfel(real(q0, c_double))
  end subroutine CachePDFsAPFEL

  double precision function xPDF(i, x)
    integer, intent(in) :: i
    double precision, intent(in) :: x
    xPDF = apfelxxf_xpdf(int(i, c_int), real(x, c_double))
  end function xPDF

  double precision function xPDFxQ(i, x, q)
    integer, intent(in) :: i
    double precision, intent(in) :: x, q
    xPDFxQ = apfelxxf_xpdfxq(int(i, c_int), real(x, c_double), real(q, c_double))
  end function xPDFxQ

  subroutine SetRenFacRatio(xi)
    double precision, intent(in) :: xi
    call apfelxxf_setrenfacratio(real(xi, c_double))
  end subroutine SetRenFacRatio

  subroutine SetReplica(irep)
    integer, intent(in) :: irep
    call apfelxxf_setreplica(int(irep, c_int))
  end subroutine SetReplica

end module apfel_fortran
