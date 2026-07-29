!
! Smoke test for the apfelxx Fortran wrapper build wiring: confirms that
! a real Fortran program can link against libapfelxxfortran and call
! through to the C++ shim.
!
program checkapfel_test
  use apfel_evolution
  implicit none

  if (.not. CheckAPFEL()) then
    print *, "FAIL: CheckAPFEL() returned .false."
    stop 1
  end if

  print *, "PASS: CheckAPFEL() returned .true."
end program checkapfel_test
