! Copyright 2017, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
program main

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

    status = BelosPassed
    output_detail = BelosErrors + BelosWarnings + BelosFinalSummary

  end subroutine

end program
