!
! Validates SetRenFacRatio: sets xi = muR/muF = 2 (BuildDglap's scale-
! variation parameter, hardcoded to its default of 1 before this
! addition) and checks that EvolveAPFEL/xPDF picks it up, against
! reference values from an independent, standalone C++ program built
! directly against apfelxx with the same BuildDglap(..., xi = 2) call.
!
program renfacratio_test
  use apfel_fortran
  implicit none

  double precision, parameter :: sqrt2 = 1.4142135623730951d0
  double precision, parameter :: Q0 = sqrt2, Q = 100d0, x = 0.1d0
  double precision, parameter :: tol = 1d-9

  double precision, parameter :: ref_gluon = 0.8795946987924586d0
  double precision, parameter :: ref_up    = 0.6578866471090744d0
  double precision, parameter :: ref_ubar  = 0.0952085953462092d0

  call SetPerturbativeOrder(1)
  call SetPoleMasses(sqrt2, 4.5d0, 175d0)
  call SetAlphaQCDRef(0.118d0, 91.1876d0)
  call SetQLimits(1d0, 1000d0)

  call SetNumberOfGrids(4)
  call SetGridParameters(1, 100, 3, 1d-5)
  call SetGridParameters(2, 100, 3, 1d-1)
  call SetGridParameters(3, 100, 3, 6d-1)
  call SetGridParameters(4, 80,  5, 8.5d-1)
  call SetQGridParameters(50, 3)

  call InitializeAPFEL()
  call SetPDFSet("toyLH")

  call SetRenFacRatio(2d0)
  call EvolveAPFEL(Q0, Q)

  call check("xPDF(0,x)  [gluon]", xPDF(0, x),  ref_gluon)
  call check("xPDF(2,x)  [up]   ", xPDF(2, x),  ref_up)
  call check("xPDF(-2,x) [ubar] ", xPDF(-2, x), ref_ubar)

  print *, "PASS"

contains

  subroutine check(label, computed, expected)
    character(len=*), intent(in) :: label
    double precision, intent(in) :: computed, expected
    if (abs(computed - expected) > tol) then
      print *, "FAIL: ", label, " = ", computed, " (expected ", expected, ")"
      stop 1
    end if
  end subroutine check

end program renfacratio_test
