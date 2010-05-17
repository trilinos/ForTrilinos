#include "breaklines.h"

/**********************************************************
 ** Unittest implementations
 **********************************************************/

#if !defined(MAIN_UNITTEST_DRIVER) && !defined(MAIN_UNITTEST_DRIVER_CXX)

/* Module level */

#  define FORTRILINOS_UNITTEST_MODULE_DEF0(CLASSNAME) \
      module TEST_FILE(CLASSNAME) BREAKLINE \
        implicit none BREAKLINE \
        public BREAKLINE \
        contains

#  define FORTRILINOS_UNITTEST_MODULE_DEF(CLASSNAME) \
      FORTRILINOS_UNITTEST_MODULE_DEF0(CLASSNAME)

#  define FORTRILINOS_UNITTEST_MODULE_BEGIN(CLASSNAME)

#  define FORTRILINOS_UNITTEST_MODULE_END0(CLASSNAME) \
      end module TEST_FILE(CLASSNAME)

#  define FORTRILINOS_UNITTEST_MODULE_END(CLASSNAME) \
      FORTRILINOS_UNITTEST_MODULE_END0(CLASSNAME)

/* Unittest level */

#  define FORTRILINOS_UNITTEST_DEF(CLASSNAME, METHODNAME) \
      logical function TEST_NAME(CLASSNAME, METHODNAME)() result(success)

#  define FORTRILINOS_UNITTEST_BEGIN

#  define FORTRILINOS_UNITTEST_END \
      end function

#endif
