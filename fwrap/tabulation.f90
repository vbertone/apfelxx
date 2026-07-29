!
! Fortran driver reproducing the structure of the old APFEL example
! program examples/Tabulation.f on top of the new apfelxx Fortran
! wrapper (fwrap/apfel_evolution.f90), evolving from Q0 = sqrt(2) GeV
! to Q = 100 GeV (non-interactive, unlike the original).
!
! No SetXxx calls are made at all: InitializeAPFEL falls back to its
! built-in defaults for the perturbative order, quark masses/thresholds,
! Q range, x/Q grids, and the "toyLH" input PDF set, all matching
! apfel::EvolutionSetup's own defaults (inc/apfel/evolutionsetup.h).
!
program tabulation
  use apfel_evolution
  implicit none

  integer, parameter :: nlha = 9
  double precision :: xlha(nlha)
  double precision, parameter :: Q0 = 1.4142135623730951d0, Q = 100d0
  integer :: ilha

  data xlha / 1d-5, 1d-4, 1d-3, 1d-2, 1d-1, 3d-1, 5d-1, 7d-1, 9d-1 /

  ! Initializes integrals on the grids
  call InitializeAPFEL()

  ! Evolve PDFs on the grids from Q0 = sqrt(2) GeV to Q = 100 GeV
  call EvolveAPFEL(Q0, Q)

  ! Tabulate PDFs for the LHA x values
  write(6,*) "alpha_QCD(mu2F) =", alphaQCD(Q)
  write(6,*) " "

  write(6,*) "Standard evolution:"
  write(6,'(a7,5a12)') "x", "u-ubar", "d-dbar", "2(ubr+dbr)", "c+cbar", "gluon"
  do ilha = 1, nlha
    write(6,'(es7.1,5es12.4)') xlha(ilha), &
      xPDF(2, xlha(ilha)) - xPDF(-2, xlha(ilha)), &
      xPDF(1, xlha(ilha)) - xPDF(-1, xlha(ilha)), &
      2d0 * (xPDF(-1, xlha(ilha)) + xPDF(-2, xlha(ilha))), &
      xPDF(4, xlha(ilha)) + xPDF(-4, xlha(ilha)), &
      xPDF(0, xlha(ilha))
  end do
  write(6,*) " "

  ! Cached PDFs
  call CachePDFsAPFEL(Q0)

  write(6,*) "Cached evolution:"
  write(6,'(a7,5a12)') "x", "u-ubar", "d-dbar", "2(ubr+dbr)", "c+cbar", "gluon"
  do ilha = 1, nlha
    write(6,'(es7.1,5es12.4)') xlha(ilha), &
      xPDFxQ(2, xlha(ilha), Q) - xPDFxQ(-2, xlha(ilha), Q), &
      xPDFxQ(1, xlha(ilha), Q) - xPDFxQ(-1, xlha(ilha), Q), &
      2d0 * (xPDFxQ(-1, xlha(ilha), Q) + xPDFxQ(-2, xlha(ilha), Q)), &
      xPDFxQ(4, xlha(ilha), Q) + xPDFxQ(-4, xlha(ilha), Q), &
      xPDFxQ(0, xlha(ilha), Q)
  end do
  write(6,*) " "

end program tabulation
