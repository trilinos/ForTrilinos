! vi: ft=fortran
! Copyright 2017-2018, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
#ifndef FOTRANTESTUTILITIES_H
#define FOTRANTESTUTILITIES_H

! This header contains several `use` statements and, therefore, *must* be `use`d
! *inside* a `program` or `module`.  Otherwise, the compile will fail

! Several macros are defined to ease writing tests in ForTrilinos.  Macros exist
! for setting up a test, defining unit tests, tearing down tests, checking
! variable equality, etc.  Creating a unit test involves several steps:
!
! - Declaring the main program
! - Including this header file inside of the program
! - Declaring program-wide variables used in the main program's program units
! - Setting up the test globally.  Use the macro SETUP_TEST:
!
!      SETUP_TEST()
!
!   The SETUP_TEST macro initializes module variables in the fortest module and
!   initializes MPI for MPI builds of ForTrilinos (meaning that MPI should not
!   be initialized elsewhere in the test!)
!
! - Adding and running subtests.  Use the macro ADD_SUBTEST_AND_RUN:
!
!      ADD_SUBTEST_AND_RUN(TestName)
!
!   where TestName is the name of a test defined in the program unit.  The test
!   should be defined using the FORTRILINOS_UNIT_TEST and
!   END_FORTRILINOS_UNIT_TEST macros (after the program's "contains" directive):
!
!      FORTRILINOS_UNIT_TEST(TestName)
!        ! Test definition
!      END_FORTRILINOS_UNIT_TEST(TestName)
!
!   These macros define a subroutine named "TestName" and declare an in/out
!   logical "SUCCESS" that can be set to .FALSE. if the test fails.
!   FORTRILINOS_UNIT_TEST works directly with ADD_SUBTEST_AND_RUN, which calls
!   the FORTRILINOS_UNIT_TEST, setting SUCCESS=.TRUE. on entry, and incrementing
!   the count of tests run.
!
!   The macro END_FORTRILINOS_UNIT_TEST ends the subroutine and does a global
!   reduction (for MPI builds of ForTrilinos) on SUCCESS to determine the test
!   pass/fail result over all processors.  A test is considered to have failed
!   if SUCCESS is .FALSE. on any processor.  Error states at the completion of
!   a subtest are reset so that other subtests can be run.
!
! - Tearing down the test.  Use the macro TEARDOWN_TEST
!
!      TEARDOWN_TEST()
!
!   This macro prints any error diagnostics to the console, shuts down MPI (on
!   MPI builds of ForTrilinos), and stops the program with an error flag if any
!   of the subtests fail.  This, in turn, is communicated to CTest to get test
!   failure/pass metrics.
!
! - Ending the test program
!
! Additionally, several macros for testing equality of objects, asserting
! required conditionals, testing code, etc. are defined.  Look further down in
! this file for macros starting "TEST_".
!
! See tests in src/tpetra/test for examples of creating unit tests and using the
! TEST_* macros.

#include "ForTrilinos.h"
use fortest

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
 FORTRILINOS_CHECK_IERR(); \
 CALL TEARDOWN_TEST2()

! Macro defining a ForTrilinos unit test.  Defines subroutine signature and
! declares local variable `SUCCESS` for monitoring test status
#define FORTRILINOS_UNIT_TEST(NAME) \
 SUBROUTINE NAME(SUCCESS); \
 IMPLICIT NONE; \
 LOGICAL :: SUCCESS

! Macro defining the end of a ForTrilinos unit test.  Sets the value of SUCCESS
! and resets FORTRILINOS_IERR so that downstream tests can run.
#define END_FORTRILINOS_UNIT_TEST(NAME) \
 SUCCESS = SUCCESS .AND. (FORTRILINOS_IERR == 0); \
 FORTRILINOS_IERR = 0; \
 RETURN; \
 END SUBROUTINE NAME

! The TEST_* macros to follow are intended to be called from within
! a FORTRILINOS_UNIT_TEST.  Instead of stopping calculations at the first sign of
! error, they set the variable "SUCCESS" to .FALSE. and return.  This informs the
! ADD_SUBTEST_AND_RUN macro that the subtest failed.

! Checks for the value of FORTRILINOS_IERR and toggles SUCCESS to .FALSE. if FORTRILINOS_IERR /= 0.
! If the SUCCESS flag is toggled, the test is exited
#define TEST_IERR() \
 CALL FORTEST_IERR(SUCCESS, FILENAME, __LINE__); \
 IF (.NOT.SUCCESS) THEN; FORTRILINOS_IERR = 0; RETURN; ENDIF

! Checks integers A and B are the same and toggles SUCCESS to .FALSE. if they
! are not. If the SUCCESS flag is toggled, the test is exited
#define TEST_EQUALITY(A, B) \
 CALL FORTEST_EQUALITY(SUCCESS, FILENAME, __LINE__, "A", A, "B", B); \
 IF (.NOT.SUCCESS) THEN; FORTRILINOS_IERR = 0; RETURN; ENDIF

! Checks integers A and B are not the same and toggles SUCCESS to .FALSE. if they
! are. If the SUCCESS flag is toggled, the test is exited
#define TEST_INEQUALITY(A, B) \
 CALL FORTEST_INEQUALITY(SUCCESS, FILENAME, __LINE__, "A", A, "B", B); \
 IF (.NOT.SUCCESS) THEN; FORTRILINOS_IERR = 0; RETURN; ENDIF

! Checks integer arrays A and B are the same and toggles SUCCESS to .FALSE.
! if they are not. If the SUCCESS flag is toggled, the test is exited
#define TEST_ARRAY_EQUALITY(A, B) \
 CALL FORTEST_ARRAY_EQUALITY(SUCCESS, FILENAME, __LINE__, "A", A, "B", B); \
 IF (.NOT.SUCCESS) THEN; FORTRILINOS_IERR = 0; RETURN; ENDIF

! Checks integer arrays A and B are not the same and toggles SUCCESS to .FALSE.
! if they are. If the SUCCESS flag is toggled, the test is exited
#define TEST_ARRAY_INEQUALITY(A, B) \
 CALL FORTEST_ARRAY_INEQUALITY(SUCCESS, FILENAME, __LINE__, "A", A, "B", B); \
 IF (.NOT.SUCCESS) THEN; FORTRILINOS_IERR = 0; RETURN; ENDIF

! Checks real numbers A and B are within tolerance and toggles SUCCESS to
! .FALSE. if they are not. If the SUCCESS flag is toggled, the test is exited
#define TEST_FLOATING_EQUALITY(A, B, T) \
 CALL FORTEST_EQUALITY(SUCCESS, FILENAME, __LINE__, "A", A, "B", B, T); \
 IF (.NOT.SUCCESS) THEN; FORTRILINOS_IERR = 0; RETURN; ENDIF

! Checks real numbers A and B are not within tolerance and toggles SUCCESS to
! .FALSE. if they are. If the SUCCESS flag is toggled, the test is exited
#define TEST_FLOATING_INEQUALITY(A, B, T) \
 CALL FORTEST_INEQUALITY(SUCCESS, FILENAME, __LINE__, "A", A, "V", B, T); \
 IF (.NOT.SUCCESS) THEN; FORTRILINOS_IERR = 0; RETURN; ENDIF

! Checks array of real numbers A and B are within tolerance and toggles SUCCESS
! to .FALSE. if they are not. If the SUCCESS flag is toggled, the test is exited
#define TEST_FLOATING_ARRAY_EQUALITY(A, B, T) \
 CALL FORTEST_ARRAY_EQUALITY(SUCCESS, FILENAME, __LINE__, "A", A, "B", B, T); \
 IF (.NOT.SUCCESS) THEN; FORTRILINOS_IERR = 0; RETURN; ENDIF

! Checks array of real numbers A and B are not within tolerance and toggles
! SUCCESS to .FALSE. if they are. If the SUCCESS flag is toggled, the test
! is exited
#define TEST_FLOATING_ARRAY_INEQUALITY(A, B, T) \
 CALL FORTEST_ARRAY_INEQUALITY(SUCCESS, FILENAME, __LINE__, "A", A, "B", B, T); \
 IF (.NOT.SUCCESS) THEN; FORTRILINOS_IERR = 0; RETURN; ENDIF

! Checks that condition C is .TRUE. and toggles success to .FALSE. if not.  If
! the SUCCESS flag is toggled, the test is exited
#define TEST_ASSERT(C) \
 CALL FORTEST_ASSERT(SUCCESS, FILENAME, __LINE__, "C", C); \
 IF (.NOT.SUCCESS) THEN; FORTRILINOS_IERR = 0; RETURN; ENDIF

! Checks that instruction CODE throws an exception toggles success to .FALSE. if
! not.  If the SUCCESS flag is toggled, the test is exited
#define TEST_THROW(CODE) \
 CODE; \
 CALL FORTEST_THROW(SUCCESS, FILENAME, __LINE__, "CODE"); \
 IF (.NOT.SUCCESS) THEN; FORTRILINOS_IERR = 0; RETURN; ENDIF

! Checks that instruction CODE does not throw an exception toggles success to
! .FALSE. if not.  If the SUCCESS flag is toggled, the test is exited
#define TEST_NOTHROW(CODE) \
 CODE; \
 CALL FORTEST_NOTHROW(SUCCESS, FILENAME, __LINE__, "CODE"); \
 IF (.NOT.SUCCESS) THEN; FORTRILINOS_IERR = 0; RETURN; ENDIF

! Writes string S to stderr (file unit 0) on processor 0
#define OUT0(S) \
 CALL WRITE_TO_PROC0(S)

#endif /* FOTRANTESTUTILITIES_H */
