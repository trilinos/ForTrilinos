/* Add double quotes around anything */

#define STRINGIFY0(WHATEVER) #WHATEVER
#define STRINGIFY(WHATEVER) STRINGIFY0(WHATEVER)

/* Stem of file that will contain the test implementations */

#define TEST_FILE0(CLASSNAME) CLASSNAME##_test_impls
#define TEST_FILE(CLASSNAME) TEST_FILE0(CLASSNAME)
#define TEST_FILE_STR(CLASSNAME) STRINGIFY(TEST_FILE(CLASSNAME))

/* Stem of file that will contain the test calls */

#define TEST_CALLS_FILE0(CLASSNAME) CLASSNAME##_test_calls
#define TEST_CALLS_FILE(CLASSNAME) TEST_CALLS_FILE0(CLASSNAME)
#define TEST_CALLS_FILE_STR(CLASSNAME) STRINGIFY(TEST_CALLS_FILE(CLASSNAME))

/* Executable used to call one unittest for class */

#define TEST_EXEC_FILE0(CLASSNAME) CLASSNAME##_run_one_test.exe
#define TEST_EXEC_FILE(CLASSNAME) TEST_EXEC_FILE0(CLASSNAME)

/* Name of an individual method unittest */

#define TEST_SUFFIX_STR() STRINGIFY(_UnitTest)
#define TEST_NAME0(CLASSNAME, METHODNAME) CLASSNAME##_test_##METHODNAME##_UnitTest
#define TEST_NAME(CLASSNAME, METHODNAME) TEST_NAME0(CLASSNAME, METHODNAME)
#define TEST_NAME_STR(CLASSNAME, METHODNAME) STRINGIFY(TEST_NAME(CLASSNAME, METHODNAME))

/* Short name of an individual method unittest */

#define SHORT_TEST_NAME0(METHODNAME) METHODNAME##_UnitTest
#define SHORT_TEST_NAME(METHODNAME) SHORT_TEST_NAME0(METHODNAME)
#define SHORT_TEST_NAME_STR(METHODNAME) STRINGIFY(SHORT_TEST_NAME(METHODNAME))

/* More macros, depending on what we're building */

#ifdef MAIN_UNITTEST_DRIVER
#  include "fortran_driver_macros.h"
#else
#  ifdef MAIN_UNITTEST_DRIVER_CXX
#    include "cxx_driver_macros.h"
#  else
#    include "fortran_impl_macros.h"
#  endif
#endif

