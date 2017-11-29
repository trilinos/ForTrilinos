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

! Setup the test.  This procedure should be called *once* after program
! declarations are made.
#define SETUP_TEST() \
 LOGICAL :: SUCCESS; \
 CALL SETUP_TEST2()

! Add a subtest to the test and run it, checking for errors at its completion
#define ADD_SUBTEST_AND_RUN(NAME) \
 CALL SETUP_SUBTEST(SUCCESS); \
 CALL NAME(SUCCESS); \
 CALL TEARDOWN_SUBTEST("NAME", SUCCESS)

! Teardown the test.  No subtests should be added and ran after teardown
#define TEARDOWN_TEST() \
 CHECK_IERR(); \
 CALL TEARDOWN_TEST2()

! Macro defining a ForTrilinos unit test.  Defines subroutine signature and
! declares local variable `SUCCESS` for monitoring test status
#define FORTRILINOS_UNIT_TEST(NAME) \
 SUBROUTINE NAME(SUCCESS); \
 IMPLICIT NONE; \
 LOGICAL :: SUCCESS

! Macro defining the end of a ForTrilinos unit test.  Sets the value of SUCCESS
! and resets IERR so that downstream tests can run.
#define END_FORTRILINOS_UNIT_TEST(NAME) \
 SUCCESS = SUCCESS .AND. (IERR == 0); \
 IERR = 0; \
 RETURN; \
 END SUBROUTINE NAME

! The TEST_* macros to follow are intended to be called from within
! a FORTRILINOS_UNIT_TEST.  Instead of stopping calculations at the first sign of
! error, they set the variable "SUCCESS" to .FALSE. and return.  This informs the
! ADD_SUBTEST_AND_RUN macro that the subtest failed.

! Checks for the value of IERR and toggles SUCCESS to .FALSE. if IERR /= 0.
! If the SUCCESS flag is toggled, the test is exited
#define TEST_IERR() \
 CALL FORTEST_IERR(SUCCESS, __FILE__, __LINE__, IERR, SERR); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

! Checks integers A and B are the same and toggles SUCCESS to .FALSE. if they
! are not. If the SUCCESS flag is toggled, the test is exited
#define TEST_EQUALITY(A, B) \
 CALL FORTEST_EQUALITY(SUCCESS, __FILE__, __LINE__, "A", A, "B", B); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

! Checks integers A and B are not the same and toggles SUCCESS to .FALSE. if they
! are. If the SUCCESS flag is toggled, the test is exited
#define TEST_INEQUALITY(A, B) \
 CALL FORTEST_INEQUALITY(SUCCESS, __FILE__, __LINE__, "A", A, "B", B); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

! Checks integer arrays A and B are the same and toggles SUCCESS to .FALSE.
! if they are not. If the SUCCESS flag is toggled, the test is exited
#define TEST_ARRAY_EQUALITY(A, B) \
 CALL FORTEST_ARRAY_EQUALITY(SUCCESS, __FILE__, __LINE__, "A", A, "B", B); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

! Checks integer arrays A and B are not the same and toggles SUCCESS to .FALSE.
! if they are. If the SUCCESS flag is toggled, the test is exited
#define TEST_ARRAY_INEQUALITY(A, B) \
 CALL FORTEST_ARRAY_INEQUALITY(SUCCESS, __FILE__, __LINE__, "A", A, "B", B); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

! Checks real numbers A and B are within tolerance and toggles SUCCESS to
! .FALSE. if they are not. If the SUCCESS flag is toggled, the test is exited
#define TEST_FLOATING_EQUALITY(A, B, T) \
 CALL FORTEST_EQUALITY(SUCCESS, __FILE__, __LINE__, "A", A, "B", B, T); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

! Checks real numbers A and B are not within tolerance and toggles SUCCESS to
! .FALSE. if they are. If the SUCCESS flag is toggled, the test is exited
#define TEST_FLOATING_INEQUALITY(A, B, T) \
 CALL FORTEST_INEQUALITY(SUCCESS, __FILE__, __LINE__, "A", A, "V", B, T); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

! Checks array of real numbers A and B are within tolerance and toggles SUCCESS
! to .FALSE. if they are not. If the SUCCESS flag is toggled, the test is exited
#define TEST_FLOATING_ARRAY_EQUALITY(A, B, T) \
 CALL FORTEST_ARRAY_EQUALITY(SUCCESS, __FILE__, __LINE__, "A", A, "B", B, T); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

! Checks array of real numbers A and B are not within tolerance and toggles
! SUCCESS to .FALSE. if they are. If the SUCCESS flag is toggled, the test
! is exited
#define TEST_FLOATING_ARRAY_INEQUALITY(A, B, T) \
 CALL FORTEST_ARRAY_INEQUALITY(SUCCESS, __FILE__, __LINE__, "A", A, "B", B, T); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

! Checks that condition C is .TRUE. and toggles success to .FALSE. if not.  If
! the SUCCESS flag is toggled, the test is exited
#define TEST_ASSERT(C) \
 CALL FORTEST_ASSERT(SUCCESS, __FILE__, __LINE__, "C", C); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

! Checks that instruction CODE throws an exception toggles success to .FALSE. if
! not.  If the SUCCESS flag is toggled, the test is exited
#define TEST_THROW(CODE) \
 CODE; \
 CALL FORTEST_THROW(SUCCESS, __FILE__, __LINE__, "CODE", IERR); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

! Checks that instruction CODE does not throw an exception toggles success to
! .FALSE. if not.  If the SUCCESS flag is toggled, the test is exited
#define TEST_NOTHROW(CODE) \
 CODE; \
 CALL FORTEST_NOTHROW(SUCCESS, __FILE__, __LINE__, "CODE", IERR); \
 IF (.NOT.SUCCESS) THEN; IERR = 0; RETURN; ENDIF

! Writes string S to stderr (file unit 0) on processor 0
#define OUT0(S) \
 CALL WRITE_TO_PROC0(S)

#endif /* FOTRANTESTMACROS_H */
