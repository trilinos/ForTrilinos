#ifndef ScaleSTL_FortranTestMacros_H
#define ScaleSTL_FortranTestMacros_H

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

! NOTE: Gfortran and ifort do not support variadic macros, so we need several
! macros, instead of just 1
#define CALL_AND_CHECK_IERR( NAME, PROCEDURE ) \
    call PROCEDURE(); \
    if (ierr /= 0) then; \
      NAME = ierr; \
      ierr = 0; \
      return; \
    endif

#define CALL_AND_CHECK_IERR_1( NAME, PROCEDURE, X1 ) \
    call PROCEDURE( X1 ); \
    if (ierr /= 0) then; \
      NAME = ierr; \
      ierr = 0; \
      return; \
    endif

#define CALL_AND_CHECK_IERR_2( NAME, PROCEDURE, X1, X2 ) \
    call PROCEDURE( X1, X2 ); \
    if (ierr /= 0) then; \
      NAME = ierr; \
      ierr = 0; \
      return; \
    endif

#define CALL_AND_CHECK_IERR_3( NAME, PROCEDURE, X1, X2, X3 ) \
    call PROCEDURE( X1, X2, X3 ); \
    if (ierr /= 0) then; \
      NAME = ierr; \
      ierr = 0; \
      return; \
    endif

#define CALL_AND_CHECK_IERR_4( NAME, PROCEDURE, X1, X2, X3, X4 ) \
    call PROCEDURE( X1, X2, X3, X4 ); \
    if (ierr /= 0) then; \
      NAME = ierr; \
      ierr = 0; \
      return; \
    endif

#define CALL_AND_CHECK_IERR_5( NAME, PROCEDURE, X1, X2, X3, X4, X5 ) \
    call PROCEDURE( X1, X2, X3, X4, X5 ); \
    if (ierr /= 0) then; \
      NAME = ierr; \
      ierr = 0; \
      return; \
    endif

#define CALL_FCN_AND_CHECK_IERR( NAME, PROCEDURE, RES) \
    RES = PROCEDURE( ); \
    if (ierr /= 0) then; \
      NAME = ierr; \
      ierr = 0; \
      return; \
    endif

#define ADD_TEST(TEST_NAME, ERROR_COUNTER) \
    if ( TEST_NAME() /= 0 ) then; \
      write(0, *) "Test FAILED!"; \
      ERROR_COUNTER = ERROR_COUNTER + 1; \
    endif

#define SET_ERROR_COUNT_AND_RESET_IERR(TEST_NAME, ERROR_COUNT) \
    TEST_NAME = ERROR_COUNT + ierr; \
    ierr = 0


use DBCF_M

#endif
