! vi: ft=fortran
#ifndef FOTRANTESTMACROS_H
#define FOTRANTESTMACROS_H

#include "ForTrilinosUtils_config.hpp"
/* Options */

#include "DBCF.h"
use DBCF_M

#define ASSERT_C_ASSOC( X ) \
    Insist( C_ASSOCIATED( X % instance_ptr ), "X not C_ASSOCIATED" )
#define ASSERT_NOT_C_ASSOC( X ) \
    Insist(.NOT.C_ASSOCIATED( X % instance_ptr ), "X C_ASSOCIATED" )
#define ASSERT_C_ASSOC2( X, Y )                                 \
    Insist( C_ASSOCIATED( X % instance_ptr, Y % instance_ptr ), \
            "X not C_ASSOCIATED with Y" )
#define ASSERT_NOT_C_ASSOC2( X, Y )                                 \
    Insist(.NOT.C_ASSOCIATED( X % instance_ptr, Y % instance_ptr ), \
           "X C_ASSOCIATED with Y" )

! NOTE: GFortran doesn't support the '#' token to convert code to a string, so
! the given arguments MUST NOT have quotes in them
#define EXPECT_EQ( REF, TEST )               \
    if(.NOT. REF == TEST ) then;              \
    WRITE( 0, * ) "Expected: ", REF;         \
    WRITE( 0, * ) "Actual: ", TEST;          \
    Insist(.false., "EXPECT_EQ(REF,TEST)" ); \
    endif

#define EXPECT_TRUE( TEST )                \
    if(.NOT. TEST ) then;                   \
    WRITE( 0, * ) "Expected: TRUE";        \
    WRITE( 0, * ) "Actual: FALSE";         \
    Insist(.false., "EXPECT_TRUE(TEST)" ); \
    endif

#define EXPECT_FALSE( TEST )                \
    if( TEST ) then;                        \
    WRITE( 0, * ) "Expected: FALSE";        \
    WRITE( 0, * ) "Actual: TRUE";           \
    Insist(.false., "EXPECT_FALSE(TEST)" ); \
    endif

#define CHECK_IERR( ) \
    if(.NOT. ierr == 0 ) then;                                    \
    WRITE( 0, * ) "Expected ierr = 0, but got ", ierr;            \
    WRITE( *, * ) "With associated serr = ", trim(serr);          \
    Insist(.false., "Expected ierr = 0" );                        \
    end if

#define FORTRILINOS_UNIT_TEST(NAME)                               \
    subroutine NAME(success);                                     \
    implicit none;                                                \
    logical :: success

#define END_FORTRILINOS_UNIT_TEST(NAME)                           \
    if (success) success = (ierr == 0);                           \
    ierr = 0;                                                     \
    return;                                                       \
    end subroutine NAME

! Two versions of several macros are provided for with/without MPI.
#ifdef HAVE_MPI

#define DECLARE_TEST_VARIABLES()                                  \
  use mpi;                                                        \
  implicit none;                                                  \
  logical LOCAL_SUCCESS;                                          \
  character(len=256) JNKSTR;                                      \
  integer ERR1(1), ERR2(1), IERROR;                               \
  integer COMM_RANK, FAILED_TESTS, TOTAL_TESTS

#define INITIALIZE_TEST()                                         \
  ERR1(1) = 0; ERR2(1) = 0; FAILED_TESTS = 0; TOTAL_TESTS=0;      \
  call MPI_INIT(ierr);                                            \
  call MPI_COMM_RANK(MPI_COMM_WORLD, COMM_RANK, IERROR)

#define ADD_SUBTEST_AND_RUN(NAME)                                 \
    LOCAL_SUCCESS = .true.;                                       \
    TOTAL_TESTS = TOTAL_TESTS + 1;                                \
    call NAME(LOCAL_SUCCESS);                                     \
    if (LOCAL_SUCCESS) then;                                      \
    ERR1(1) = 0;                                                  \
    else;                                                         \
    ERR1(1) = 1;                                                  \
    end if;                                                       \
    call MPI_ALLREDUCE(ERR1, ERR2, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, IERROR); \
    if ( ERR2(1) /= 0 ) then;                                     \
    if (COMM_RANK == 0) write(0, *) "Test FAILED!";               \
    FAILED_TESTS = FAILED_TESTS + 1;                              \
    end if

#define SHUTDOWN_TEST()                                           \
  if (FAILED_TESTS == 0) then;                                    \
  if (COMM_RANK == 0) then;                                       \
  write(0,'(A,I3,A)') "100% of ", TOTAL_TESTS, " tests PASSED";   \
  end if;                                                         \
  else;                                                           \
  if (COMM_RANK == 0) then;                                       \
  write(*,'(I3,A,I3,A)') FAILED_TESTS, " of ", TOTAL_TESTS, " tests FAILED";   \
  end if;                                                         \
  Insist(.false., "FAILED TESTS ENCOUNTERED" );                   \
  end if;                                                         \
  call MPI_FINALIZE(ierr)

#else

#define DECLARE_TEST_VARIABLES()                                  \
  implicit none;                                                  \
  logical LOCAL_SUCCESS;                                          \
  character(len=256) JNKSTR;                                      \
  integer COMM_RANK, TOTAL_TESTS, FAILED_TESTS

#define INITIALIZE_TEST()                                         \
  FAILED_TESTS = 0; COMM_RANK = 0;

#define ADD_SUBTEST_AND_RUN(NAME)                                 \
    LOCAL_SUCCESS = .true.;                                       \
    TOTAL_TESTS = TOTAL_TESTS + 1;                                \
    call NAME(LOCAL_SUCCESS);                                     \
    if (.not. LOCAL_SUCCESS) then;                                \
    if (COMM_RANK == 0) write(0, *) "Test FAILED!";               \
    FAILED_TESTS = FAILED_TESTS + 1;                              \
    end if

#define SHUTDOWN_TEST()                                           \
  if (FAILED_TESTS == 0) then;                                    \
  write(0,'(A,I3,A)') "100% of ", TOTAL_TESTS, " tests PASSED";   \
  else;                                                           \
  write(*,'(I3,A,I3,A)') FAILED_TESTS, " of ", TOTAL_TESTS, " tests FAILED";   \
  Insist(.false., "FAILED TESTS ENCOUNTERED" );                   \
  end if

#endif

! The TEST_* macros to follow are intended to be called from within
! a FORTRILINOS_UNIT_TEST.  Insetad of stopping calculations at the first sign of
! error, they set the variable "success" to .false. and return.  This informs the
! ADD_SUBTEST_AND_RUN macro that the subtest failed.
#define TEST_IERR()                                               \
  if (ierr /= 0) then;                                            \
  if (COMM_RANK == 0) then;                                       \
  WRITE(0, '(A,A,A,I6)') "File '", __FILE__, "', line ", __LINE__; \
  write(0, '(A)') "  TEST_IERR()";                                \
  write(0, '(A)') "Error:", trim(serr);                           \
  end if;                                                         \
  success = .false.; ierr = 0;                                    \
  return;                                                         \
  end if

#define TEST_COMPARE_FLOATING_ARRAYS(ARR1, ARR2, TOL)             \
  IF(ABS(MAXVAL(ARR1 - ARR2)) > TOL) THEN;                        \
  if (COMM_RANK == 0) then;                                       \
  WRITE(0, '(A,A,A,I6)') "File '", __FILE__, "', line ", __LINE__; \
  write(0, '(A)') "  TEST_COMPARE_FLOATING_ARRAYS(ARR1,ARR2,TOL)";\
  write(0, '(A)') "Error: Expected ARR1 == ARR2";                 \
  end if;                                                         \
  success = .false.; ierr = 0;                                    \
  return;                                                         \
  end if

#define TEST_ARRAY_EQUALITY(ARR, VAL, TOL)                        \
  if(abs(maxval(ARR - VAL)) > TOL) then;                          \
  if (COMM_RANK == 0) then;                                       \
  WRITE(0, '(A,A,A,I6)') "File '", __FILE__, "', line ", __LINE__; \
  write(0, '(A)') "  TEST_ARRAY_EQUALITY(ARR1,ARR2,TOL)";         \
  write(0, '(A)') "Error: Expected ARR1 == ARR2";                 \
  end if;                                                         \
  success = .false.; ierr = 0;                                    \
  return;                                                         \
  end if

#define TEST_ARRAY_INEQUALITY(ARR, VAL, TOL)                      \
  if(abs(maxval(ARR - VAL)) <= TOL) then;                         \
  if (COMM_RANK == 0) then;                                       \
  WRITE(0, '(A,A,A,I6)') "File '", __FILE__, "', line ", __LINE__; \
  write(0, '(A)') "  TEST_ARRAY_INEQUALITY(ARR1,ARR2,TOL)";       \
  write(0, '(A)') "Error: Expected ARR1 /= ARR2";                 \
  end if;                                                         \
  success = .false.; ierr = 0;                                    \
  return;                                                         \
  end if

#define TEST_FLOATING_EQUALITY(VAL1, VAL2, TOL)                   \
  if(abs(VAL1 - VAL2) > TOL) then;                                \
  if (COMM_RANK == 0) then;                                       \
  WRITE(0, '(A,A,A,I6)') "File '", __FILE__, "', line ", __LINE__; \
  WRITE(0, '(A)') "  TEST_FLOATING_EQUALITY(VAL1,VAL2,TOL)";      \
  write(0, *) "Error: Expected ", VAL1, "==", VAL2;           \
  end if;                                                         \
  success = .false.; ierr = 0;                                    \
  return;                                                         \
  end if

#define TEST_EQUALITY(VAL1, VAL2)                                 \
  if(VAL1 /= VAL2) then;                                          \
  if (COMM_RANK == 0) then;                                       \
  WRITE(0, '(A,A,A,I6)') "File '", __FILE__, "', line ", __LINE__; \
  WRITE(0, '(A)') "  TEST_EQUALITY(VAL1,VAL2)";                   \
  write(0, *) "Error: Expected ", VAL1, "==", VAL2;           \
  end if;                                                         \
  success = .false.; ierr = 0;                                    \
  return;                                                         \
  end if

#define TEST_ASSERT(COND)                                         \
  if(.not. COND) then;                                            \
  if (COMM_RANK == 0) then;                                       \
  WRITE(0, '(A,A,A,I6)') "File '", __FILE__, "', line ", __LINE__; \
  WRITE(0, '(A)') "  TEST_ASSERT(COND)";                          \
  WRITE(0, '(A)') "Error: Expected COND to evaluate to .TRUE.";   \
  end if;                                                         \
  success = .false.; ierr = 0;                                    \
  return;                                                         \
  end if

#define TEST_THROW(CODE)                                          \
  CODE;                                                           \
  if (ierr == 0) then;                                            \
  if (COMM_RANK == 0) then;                                       \
  WRITE(0, '(A,A,A,I6)') "File '", __FILE__, "', line ", __LINE__; \
  WRITE(0, '(A)') "  TEST_THROW(CODE)";                           \
  WRITE(0, '(A)') "Error: Expected CODE to set ierr /= 0";        \
  end if;                                                         \
  success = .false.; ierr = 0;                                    \
  return;                                                         \
  end if

#define TEST_NOTHROW(CODE)                                        \
  CODE;                                                           \
  if (ierr /= 0) then;                                            \
  if (COMM_RANK == 0) then;                                       \
  WRITE(0, '(A,A,A,I6)') "File '", __FILE__, "', line ", __LINE__; \
  WRITE(0, '(A)') "  TEST_NOTHROW(CODE)";                         \
  WRITE(0, '(A)') "Error: Expected CODE to not set ierr /= 0";    \
  WRITE(0, '(A,A)') "serr: ", trim(serr);                         \
  end if;                                                         \
  success = .false.; ierr = 0;                                    \
  return;                                                         \
  end if

#define OUT0(STRING)                                              \
    if (COMM_RANK == 0) WRITE(0, *) STRING


#endif /* FOTRANTESTMACROS_H */
