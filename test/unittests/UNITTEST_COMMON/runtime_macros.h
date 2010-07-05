/* DO NOT INCLUDE ANY CONFIG-DEPENDENT MACROS HERE! */

#define BREAKLINE ;

/* Add double quotes around anything */

#define STRINGIFY0(WHATEVER) #WHATEVER
#define STRINGIFY(WHATEVER) STRINGIFY0(WHATEVER)

/* Macros that can be called from within unittests */

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

#define TEST_EQUIV(x, y) \
  print *, "TEST: ", #x, " = ", (x), " ?.eqv.? ", #y, " = ", (y) BREAKLINE \
  if ((x) .NEQV. (y)) then BREAKLINE \
    success = .FALSE. BREAKLINE \
    print *, "Assertion failed on line ", STRINGIFY(__LINE__) BREAKLINE \
  endif


#define TEST_LESSEQUAL(x, y) \
  print *, "TEST: ", #x, " = ", (x), " ?<=? ", #y, " = ", (y) BREAKLINE \
    if (.not.((x) .LE. (y))) then BREAKLINE				\
    success = .FALSE. BREAKLINE \
    print *, "Assertion failed on line ", STRINGIFY(__LINE__) BREAKLINE \
  endif
