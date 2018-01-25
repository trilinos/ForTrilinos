%module forerror

%include "copyright.i"

#define SWIG_FORTRAN_ERROR_INT fortrilinos_ierr
#define SWIG_FORTRAN_ERROR_STR fortrilinos_get_serr

%include <exception.i>
