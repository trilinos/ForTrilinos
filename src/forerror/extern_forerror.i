/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

/* Include this file directly instead of importing forerror to generate the
 * correct external linkage for the error variables:
      %include <extern_forerror.i>
 * rather than
      %import <forerror.i>
 *
 * This must be included before any other exception handling code, so it's best
 * placed at the very top of a module.
 */
#define SWIG_FORTRAN_ERROR_INT fortrilinos_ierr
#define SWIG_FORTRAN_ERROR_STR fortrilinos_get_serr

%include <extern_exception.i>
%import "forerror/forerror.i"
