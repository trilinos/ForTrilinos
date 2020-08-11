! vi: ft=fortran
! Copyright 2017-2018, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
#ifndef FORTRILINOS_H
#define FORTRILINOS_H

use forerror

#define FORTRILINOS_CHECK_IERR() \
 IF(FORTRILINOS_IERR/=0) THEN; \
   WRITE(error_unit, '(A)') '*** ForTrilinos caught exception!'//NEW_LINE('A')//FORTRILINOS_GET_SERR(); \
   STOP 1; \
 ENDIF

#endif
