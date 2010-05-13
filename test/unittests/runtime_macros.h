#define PRINTIT(x) print *, #x

#define DOIT(x) x

#define ECHO(x) \
  PRINTIT(x) @@ \
  DOIT(x)

#define TEST_EQUALITY(x, y) \
  print *, "TEST: ", #x, " = ", (x), " ?==? ", #y, " = ", (y) @@ \
  if ((x) .NE. (y)) then @@ \
    success = .FALSE. @@ \
    print *, "Assertion failed on line ", STRINGIFY(__LINE__) @@ \
  endif

#define TEST_INEQUALITY(x, y) \
  print *, "TEST: ", #x, " = ", (x), " ?!=? ", #y, " = ", (y) @@ \
  if ((x) .EQ. (y)) then @@ \
    success = .FALSE. @@ \
    print *, "Assertion failed on line ", STRINGIFY(__LINE__) @@ \
  endif

