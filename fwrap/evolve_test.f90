!
! Validates SetPDFSet/EvolveAPFEL/xPDF and, independently,
! CachePDFsAPFEL/xPDFxQ: evolves the built-in toy PDF set from Q0 =
! sqrt(2) GeV to Q = 100 GeV at NLO and checks both against reference
! values obtained from an independent, standalone C++ program built
! directly against apfelxx (same Grid, AlphaQCD,
! InitializeDglapObjectsQCDOme, BuildDglap, LHToyPDFs, QCDEvToPhys
! construction).
!
! xPDF/EvolveAPFEL (direct dglap->Evaluate) and xPDFxQ/CachePDFsAPFEL
! (TabulateObject::EvaluateMapxQ) are deliberately independent of each
! other here: CachePDFsAPFEL is called without ever calling
! EvolveAPFEL, and xPDFxQ is checked without touching xPDF.
!
program evolve_test
  use apfel_evolution
  implicit none

  double precision, parameter :: sqrt2 = 1.4142135623730951d0
  double precision, parameter :: Q0 = sqrt2, Q = 100d0, x = 0.1d0
  double precision, parameter :: tol = 1d-9

  ! Reference values (EvolveAPFEL/xPDF: direct dglap->Evaluate)
  double precision, parameter :: ref_gluon = 0.8524539070427345d0
  double precision, parameter :: ref_up    = 0.6450437853482118d0
  double precision, parameter :: ref_ubar  = 0.0925297723004884d0

  ! Reference values (CachePDFsAPFEL/xPDFxQ: TabulateObject::EvaluateMapxQ)
  double precision, parameter :: refc_gluon = 0.8524542968574440d0
  double precision, parameter :: refc_up    = 0.6450438064814029d0
  double precision, parameter :: refc_ubar  = 0.0925297722520588d0

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

  ! EvolveAPFEL/xPDF
  call EvolveAPFEL(Q0, Q)
  call check("xPDF(0,x)  [gluon]", xPDF(0, x),  ref_gluon)
  call check("xPDF(2,x)  [up]   ", xPDF(2, x),  ref_up)
  call check("xPDF(-2,x) [ubar] ", xPDF(-2, x), ref_ubar)

  ! CachePDFsAPFEL/xPDFxQ, independent of EvolveAPFEL/xPDF above
  call CachePDFsAPFEL(Q0)
  call check("xPDFxQ(0,x,Q)  [gluon]", xPDFxQ(0, x, Q),  refc_gluon)
  call check("xPDFxQ(2,x,Q)  [up]   ", xPDFxQ(2, x, Q),  refc_up)
  call check("xPDFxQ(-2,x,Q) [ubar] ", xPDFxQ(-2, x, Q), refc_ubar)

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

end program evolve_test
