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

#define CHECK_IERR( ) \
    if(.NOT. ierr == 0 ) then;                         \
    WRITE( 0, * ) "Expected ierr = 0, but got ", ierr; \
    Insist(.false., "Expected ierr = 0" );             \
    endif

#define TEST_FOR_IERR( NAME ) \
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

#define SET_ERROR_COUNT_AND_RETURN(NAME, ERROR_COUNT) \
    NAME = ERROR_COUNT + ierr; \
    ierr = 0; \
    return

use DBCF_M

#endif
