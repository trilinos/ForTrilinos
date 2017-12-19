! vi: ft=fortran
! Copyright 2017, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
#ifndef FORTILINOS_H
#define FORTILINOS_H

use forerror

#define FORTRILINOS_CHECK_IERR() \
 IF(FORTRILINOS_IERR/=0) THEN; \
   WRITE(error_unit, '(A)') '*** ForTrilinos caught exception!'//NEW_LINE('A')//TRIM(FORTRILINOS_GET_SERR()); \
   STOP 1; \
 ENDIF

#endif
