%module forerror

%include "fortrilinos_copyright.i"

/* Rename exception variables and set up exception handling */
#define SWIG_FORTRAN_ERROR_INT fortrilinos_ierr
#define SWIG_FORTRAN_ERROR_STR fortrilinos_get_serr
%include <exception.i>

/* -------------------------------------------------------------------------
 * Version information
 * ------------------------------------------------------------------------- */

%apply char* { const char fortrilinos_version[] };
%fortranbindc fortrilinos_version_major;
%fortranbindc fortrilinos_version_minor;
%fortranbindc fortrilinos_version_patch;

// These symbols are defined in the CMake-generated `fortrilinos_version.cpp`
%inline %{
extern "C" {
extern const char fortrilinos_version[];
extern const int fortrilinos_version_major;
extern const int fortrilinos_version_minor;
extern const int fortrilinos_version_patch;
}
%}
