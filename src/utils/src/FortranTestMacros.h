! vi: ft=fortran
! Copyright 2017, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
#ifndef FOTRANTESTMACROS_H
#define FOTRANTESTMACROS_H

! This header contains several `use` statements and, therefore, *must* be `use`d
! *inside* a `program` or `module`.  Otherwise, the compile will fail

#include "DBCF.h"
use DBCF_M
use fortest

#define ASSERT_C_ASSOC(X) \
 Insist(C_ASSOCIATED(X%instance_ptr), "X not C_ASSOCIATED")
#define ASSERT_NOT_C_ASSOC(X) \
 Insist(.NOT.C_ASSOCIATED(X%instance_ptr), "X C_ASSOCIATED")
#define ASSERT_C_ASSOC2(X, Y) \
 Insist(C_ASSOCIATED(X%instance_ptr, Y%instance_ptr), \
        "X not C_ASSOCIATED with Y")
#define ASSERT_NOT_C_ASSOC2(X, Y) \
 Insist(.NOT.C_ASSOCIATED(X%instance_ptr, Y%instance_ptr), \
        "X C_ASSOCIATED with Y")

! NOTE: GFortran doesn't support the '#' token to convert code to a string, so
! the given arguments MUST NOT have quotes in them
#define EXPECT_EQ(REF, TEST) \
 IF(.NOT. REF == TEST) THEN; \
 WRITE(0, *) "Expected: ", REF; \
 WRITE(0, *) "Actual: ", TEST; \
 Insist(.FALSE., "EXPECT_EQ(REF,TEST)"); \
 ENDIF

#define EXPECT_TRUE(TEST) \
 IF(.NOT. TEST) THEN; \
 WRITE(0, '(A)') "Expected: TRUE"; \
 WRITE(0, '(A)') "Actual: FALSE"; \
 Insist(.FALSE., "EXPECT_TRUE(TEST)"); \
 ENDIF

#define EXPECT_FALSE(TEST) \
 IF(TEST) THEN; \
 WRITE(0, '(A)') "Expected: FALSE"; \
 WRITE(0, '(A)') "Actual: TRUE"; \
 Insist(.FALSE., "EXPECT_FALSE(TEST)"); \
 ENDIF

#define CHECK_IERR() \
 IF(IERR/=0) THEN; \
 Insist(.FALSE.,'*** ForTrilinos caught exception!'//NEW_LINE('A')//TRIM(SERR)); \
 ENDIF

#define SETUP_TEST() \
 LOGICAL :: SUCCESS; \
 CALL SETUP_TEST2()

#define ADD_SUBTEST_AND_RUN(NAME) \
 CALL SETUP_SUBTEST(SUCCESS); \
 CALL NAME(SUCCESS); \
 CALL TEARDOWN_SUBTEST("NAME", SUCCESS)

#define TEARDOWN_TEST() \
 CHECK_IERR(); \
 CALL TEARDOWN_TEST2()

#define FORTRILINOS_UNIT_TEST(NAME) \
 SUBROUTINE NAME(SUCCESS); \
 IMPLICIT NONE; \
 LOGICAL :: SUCCESS

#define END_FORTRILINOS_UNIT_TEST(NAME) \
 SUCCESS = SUCCESS .AND. (IERR == 0); \
 IERR = 0; \
 RETURN; \
 END SUBROUTINE NAME

! The TEST_* macros to follow are intended to be called from within
! a FORTRILINOS_UNIT_TEST.  Instead of stopping calculations at the first sign of
! error, they set the variable "SUCCESS" to .FALSE. and return.  This informs the
! ADD_SUBTEST_AND_RUN macro that the subtest failed.

#define TEST_IERR() \
 CALL FORTEST_IERR(SUCCESS, __FILE__, __LINE__, IERR, SERR); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

#define TEST_EQUALITY(A, B) \
 CALL FORTEST_EQUALITY(SUCCESS, __FILE__, __LINE__, "A", int(A), "B", int(B)); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

#define TEST_INEQUALITY(A, B) \
 CALL FORTEST_INEQUALITY(SUCCESS, __FILE__, __LINE__, "A", int(A), "B", int(B)); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

#define TEST_ARRAY_EQUALITY(A, B) \
 CALL FORTEST_ARRAY_EQUALITY(SUCCESS, __FILE__, __LINE__, "A", A, "B", B); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

#define TEST_ARRAY_INEQUALITY(A, B) \
 CALL FORTEST_ARRAY_INEQUALITY(SUCCESS, __FILE__, __LINE__, "A", A, "B", B); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

#define TEST_FLOATING_EQUALITY(A, B, T) \
 CALL FORTEST_EQUALITY(SUCCESS, __FILE__, __LINE__, "A", A, "B", B, T); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

#define TEST_FLOATING_INEQUALITY(A, B, T) \
 CALL FORTEST_INEQUALITY(SUCCESS, __FILE__, __LINE__, "A", A, "V", B, T); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

#define TEST_FLOATING_ARRAY_EQUALITY(A, B, T) \
 CALL FORTEST_ARRAY_EQUALITY(SUCCESS, __FILE__, __LINE__, "A", A, "B", B, T); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

#define TEST_FLOATING_ARRAY_INEQUALITY(A, B, T) \
 CALL FORTEST_ARRAY_INEQUALITY(SUCCESS, __FILE__, __LINE__, "A", A, "B", B, T); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

#define TEST_ASSERT(C) \
 CALL FORTEST_ASSERT(SUCCESS, __FILE__, __LINE__, "C", C); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

#define TEST_THROW(CODE) \
 CODE; \
 CALL FORTEST_THROW(SUCCESS, __FILE__, __LINE__, "CODE", IERR); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

#define TEST_NOTHROW(CODE) \
 CODE; \
 CALL FORTEST_NOTHROW(SUCCESS, __FILE__, __LINE__, "CODE", IERR); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

#define OUT0(S) \
 CALL WRITE_TO_PROC0(S)

#endif /* FOTRANTESTMACROS_H */
