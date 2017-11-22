#ifndef FOTRANTESTMACROS_H
#define FOTRANTESTMACROS_H

#include "ForTrilinosUtils_config.hpp"
/* Options */

#include "DBCF.h"

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
    if(.NOT. ierr == 0 ) then;                         \
    WRITE( 0, * ) "Expected ierr = 0, but got ", ierr; \
    Insist(.false., "Expected ierr = 0" );             \
    end if

#define TEST_IERR() \
    if (ierr /= 0) then; \
      success = .false.; \
      ierr = 0; \
      return; \
    end if

#define TEST_COMPARE_FLOATING_ARRAYS(ARR1, ARR2, TOL)          \
    if(abs(maxval(ARR1 - ARR2)) > TOL) then;                   \
    write( 0, * ) "Floating point arrays are not the same";    \
    success = .false.; ierr = 0;                               \
    return;                                                    \
    end if

#define TEST_ARRAY_EQUALITY(ARR, VAL, TOL)                     \
    if(abs(maxval(ARR - VAL)) > TOL) then;                     \
    write( 0, * ) "Elements of array are NOT equal to ", VAL;  \
    success = .false.; ierr = 0;                               \
    return;                                                    \
    end if

#define TEST_ARRAY_INEQUALITY(ARR, VAL, TOL)                   \
    if(abs(maxval(ARR - VAL)) <= TOL) then;                    \
    write( 0, * ) "Elements of array ARE equal to ", VAL;      \
    success = .false.; ierr = 0;                               \
    return;                                                    \
    end if

#define TEST_FLOATING_EQUALITY(VAL1, VAL2, TOL)              \
    if(abs(VAL1 - VAL2) > TOL) then;                         \
    WRITE( 0, * ) "Floating point numbers not the same" ;    \
    WRITE( 0, * ) VAL1, "/=", VAL2;                          \
    success = .false.; ierr = 0;                             \
    return;                                                  \
    end if

#define TEST_EQUALITY(VAL1, VAL2)                            \
    if(VAL1 /= VAL2) then;                                   \
    WRITE( 0, * ) "Values are not the same";                 \
    WRITE( 0, * ) VAL1, "/=", VAL2;                          \
    success = .false.; ierr = 0;                             \
    return;                                                  \
    end if

#define TEST_LOGICAL_EQUALITY(VAL1, VAL2)                    \
    if(VAL1 .neqv. VAL2) then;                               \
    WRITE( 0, * ) "Values are not the same";                 \
    WRITE( 0, * ) VAL1, "/=", VAL2;                          \
    success = .false.; ierr = 0;                             \
    return;                                                  \
    end if

#define FORTRILINOS_UNIT_TEST(NAME) \
    subroutine NAME(success);       \
    implicit none;                  \
    logical :: success

#define END_FORTRILINOS_UNIT_TEST(NAME) \
    if (success) success = (ierr == 0); \
    ierr = 0;                           \
    return;                             \
    end subroutine NAME

#define OUT0(STRING) \
    if (COMM_RANK == 0) WRITE(0, *) STRING

use DBCF_M

#ifdef HAVE_MPI

#define DECLARE_TEST_VARIABLES()                                \
  use mpi;                                                      \
  implicit none;                                                \
  logical LOCAL_SUCCESS;                                        \
  integer ERR1(1), ERR2(1), IERROR, COMM_RANK, ERROR_COUNTER

#define INITIALIZE_TEST()                                       \
  ERR1(1) = 0; ERR2(1) = 0; ERROR_COUNTER = 0;                  \
  call MPI_INIT(ierr);                                          \
  call MPI_COMM_RANK(MPI_COMM_WORLD, COMM_RANK, IERROR)

#define ADD_SUBTEST_AND_RUN(NAME)                               \
    LOCAL_SUCCESS = .true.;                                     \
    call NAME(LOCAL_SUCCESS);                                   \
    if (LOCAL_SUCCESS) then;                                    \
    ERR1(1) = 0;                                                \
    else;                                                       \
    ERR1(1) = 1;                                                \
    end if;                                                     \
    call MPI_ALLREDUCE(ERR1, ERR2, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, IERROR); \
    if ( ERR2(1) /= 0 ) then;                                   \
      if (COMM_RANK == 0) write(0, *) "Test FAILED!";           \
      ERROR_COUNTER = ERROR_COUNTER + 1;                        \
    endif

#define SHUTDOWN_TEST()                                         \
  if (ERROR_COUNTER == 0) then;                                 \
    if (COMM_RANK == 0) then;                                   \
      write(*,*) "Test PASSED";                                 \
    end if;                                                     \
  else;                                                         \
    if (COMM_RANK == 0) then;                                   \
      write(*,*) "A total of ", ERROR_COUNTER, " tests FAILED"; \
    end if;                                                     \
    Insist(.false., "FAILED TESTS ENCOUNTERED" );               \
  end if;                                                       \
  call MPI_FINALIZE(ierr)

#else

#define DECLARE_TEST_VARIABLES() \
  implicit none; \
  logical LOCAL_SUCCESS; \
  integer ERROR_COUNTER, COMM_RANK

#define INITIALIZE_TEST() \
  ERROR_COUNTER = 0; COMM_RANK = 0;

#define ADD_SUBTEST_AND_RUN(NAME) \
    LOCAL_SUCCESS = .true.
    call NAME(LOCAL_SUCCESS); \
    if ( .not. LOCAL_SUCCESS ) then; \
      write(0, *) "Test FAILED!"; \
      ERROR_COUNTER = ERROR_COUNTER + 1; \
    endif

#define SHUTDOWN_TEST() \
  if (ERROR_COUNTER == 0) then; \
    write(*,*) "Test PASSED"; \
  else; \
    write(*,*) "A total of ", ERROR_COUNTER, " tests FAILED"; \
    Insist(.false., "FAILED TESTS ENCOUNTERED" ); \
  end if

#endif

#endif /* FOTRANTESTMACROS_H */
