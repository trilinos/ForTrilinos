/* DO NOT INCLUDE ANY CONFIG-DEPENDENT MACROS HERE! */

#include "breaklines.h"

#define PRINTIT(x) print *, #x

#define DOIT(x) x

#define ECHO(x) \
  PRINTIT(x) BREAKLINE \
  DOIT(x)

#define TEST_EQUALITY(x, y) \
  print *, "TEST: ", #x, " = ", (x), " ?==? ", #y, " = ", (y) BREAKLINE \
  if ((x) .NE. (y)) then BREAKLINE \
    success = .FALSE. BREAKLINE \
    print *, "Assertion failed on line ", STRINGIFY(__LINE__) BREAKLINE \
  endif

#define TEST_INEQUALITY(x, y) \
  print *, "TEST: ", #x, " = ", (x), " ?!=? ", #y, " = ", (y) BREAKLINE \
  if ((x) .EQ. (y)) then BREAKLINE \
    success = .FALSE. BREAKLINE \
    print *, "Assertion failed on line ", STRINGIFY(__LINE__) BREAKLINE \
  endif

