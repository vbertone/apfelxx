!
! Validates InitializeAPFEL: stages perturbative order, pole masses,
! alpha_s reference value, Q limits and the default apfelxx x/Q grids,
! then checks alphaQCD against the reference value obtained from the
! existing apfelpy wrapper with the same construction (AlphaQCDRef =
! 0.118, MuQCDRef = 91.1876, NLO, same grid parameters).
!
program initialize_test
  use apfel_fortran
  implicit none

  double precision, parameter :: sqrt2 = 1.4142135623730951d0
  double precision, parameter :: expected = 0.2975093430684206d0
  double precision, parameter :: tol = 1d-10
  double precision :: computed

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

  if (GetPerturbativeOrder() /= 1) then
    print *, "FAIL: GetPerturbativeOrder() = ", GetPerturbativeOrder(), " (expected 1)"
    stop 1
  end if

  computed = alphaQCD(2d0)
  if (abs(computed - expected) > tol) then
    print *, "FAIL: alphaQCD(2) = ", computed, " (expected ", expected, ")"
    stop 1
  end if

  print *, "PASS: alphaQCD(2) = ", computed
end program initialize_test
