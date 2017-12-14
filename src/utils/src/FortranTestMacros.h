! vi: ft=fortran
! Copyright 2017, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
#ifndef FOTRANTESTMACROS_H
#define FOTRANTESTMACROS_H

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

#include "DBCF.h"
use DBCF_M

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

#endif /* FOTRANTESTMACROS_H */
