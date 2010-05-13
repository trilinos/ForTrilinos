/**********************************************************
 ** Fortran driver code
 **********************************************************/

#if defined(MAIN_UNITTEST_DRIVER)

/* Module level */

#  define FORTRILINOS_UNITTEST_MODULE_DEF0(CLASSNAME) \
    module TEST_CALLS_FILE(CLASSNAME) @@ \
      use TEST_FILE(CLASSNAME) @@ \
      implicit none @@ \
      public @@ \
      contains

#  define FORTRILINOS_UNITTEST_MODULE_DEF(CLASSNAME) \
    FORTRILINOS_UNITTEST_MODULE_DEF0(CLASSNAME)

#  define FORTRILINOS_UNITTEST_MODULE_BEGIN(CLASSNAME) \
    logical function select_test(which_test) result(success) @@ \
      character(len=50),intent(in) :: which_test @@ \
      success=.FALSE.

#  define FORTRILINOS_UNITTEST_MODULE_END0(CLASSNAME) \
    end module TEST_CALLS_FILE(CLASSNAME)

#  define FORTRILINOS_UNITTEST_MODULE_END(CLASSNAME) \
    ; stop "Missing test. TEST FAILED" @@ endif @@ end function @@ \
    FORTRILINOS_UNITTEST_MODULE_END0(CLASSNAME)

/* Unittest level */

#  define FORTRILINOS_UNITTEST_DEF(CLASSNAME, METHODNAME) \
    if (which_test==STRINGIFY(METHODNAME)) then @@ \
      success=TEST_NAME(CLASSNAME, METHODNAME)() @@ else &

#  define FORTRILINOS_UNITTEST_BEGIN #if 0

#  define FORTRILINOS_UNITTEST_END #endif

#endif
