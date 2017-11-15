program main

#include "FortranTestMacros.h"
#include "ForTrilinosBelos_config.hpp"

  use ISO_FORTRAN_ENV
  implicit none

  call test_belos_types()
contains

  subroutine test_belos_types()
    use ISO_FORTRAN_ENV
    use, intrinsic :: ISO_C_BINDING
    use forbelos
    implicit none

    integer :: status, output_detail

    status = Passed
    output_detail = Errors + Warnings + FinalSummary

  end subroutine

end program
