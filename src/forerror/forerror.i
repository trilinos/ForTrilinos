%module forerror

%include "fortrilinos_copyright.i"

/* Rename exception variables and set up exception handling */
#define SWIG_FORTRAN_ERROR_INT fortrilinos_ierr
#define SWIG_FORTRAN_ERROR_STR fortrilinos_get_serr
%include <exception.i>
