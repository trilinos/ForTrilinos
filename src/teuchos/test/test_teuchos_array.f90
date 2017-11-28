!Copyright 2017, UT-Battelle, LLC
!
!SPDX-License-Identifier: BSD-3-Clause
!License-Filename: LICENSE
program main

#include "FortranTestMacros.h"
#include "ForTrilinosTeuchos_config.hpp"

  use ISO_FORTRAN_ENV
  implicit none

  call test_array()
contains

  subroutine test_array()
    use ISO_FORTRAN_ENV
    use, intrinsic :: ISO_C_BINDING
    use forteuchos
    implicit none

    type(TeuchosArrayInt) :: arr_int

    call arr_int%create()
    EXPECT_EQ(0, ierr)

    EXPECT_EQ(arr_int%size(), 0)

    call arr_int%resize(10)
    EXPECT_EQ(0, ierr)

    ! arr_int(0) = 2
    ! EXPECT_EQ(arr_int(0), 2)

    call arr_int%release()

    call arr_int%create(10)
    EXPECT_EQ(ierr, 0)
    EXPECT_EQ(arr_int%size(), 10)

    call arr_int%release()
  end subroutine

end program
