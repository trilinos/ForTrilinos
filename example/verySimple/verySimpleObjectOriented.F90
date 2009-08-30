program main

  ! This file is the object-oriented equivalent of verySimple.F90.  In Trilinos 10.0,
  ! this is a snapshot of an unstable (evolving) file expected to become stable in a
  ! subsequent release.  This file exercises the derived types defined in 
  ! ForTrilinos/src/epetra/Epetra*.F90, which wrap the interface bodies in 
  ! ForTrilinos/src/epetra/forepetra.F90.   
    
  ! This file represents the preferred style for using ForTrilinos and is recommended for 
  ! Fortran users whose compilers support the object-oriented features of Fortran 2003.
  ! As of the Trilinos 10.0 release date, the latest versions of the IBM and Cray compilers 
  ! nominally support the required features.  The Numerical Algorithms Group (NAG) and Intel 
  ! compilers support all features but one: final subroutines.  (In each case, the support is
  ! somewhat immature and buggy.)  ForTrilinos/src/ForTrilinos_hermetic.F90 contains utilities 
  ! that help users work around the lack of final subroutines.

  ! Parts of this file that will not yet compile (in the absence of unimplemented modules) 
  ! are commented but included to exemplify the style we expect to support in a subsequent release.
 
  use ,intrinsic :: iso_c_binding ,only : c_int,c_double
! use epetra_map_module           ,only : epetra_map
! use epetra_vector_module        ,only : epetra_vector
  use interoperability_check      ,only : valid_kind_parameters
  implicit none

  if (.not. valid_kind_parameters()) stop 'C interoperability not supported on this platform.'

  ! Data declarations 

! type(epetra_map)    :: map
! type(epetra_vector) :: x, b
  integer(c_int) :: numGlobalElements_local, numGlobalElements_return
  real(c_double) :: bnorm, xnorm,err_tol,expected_bnorm,expected_xnorm,bnorm_err,xnorm_err 
  real(c_double) :: two = 2.0, zero = 0.0
  logical        :: success = .true.
  
  ! Executable code
   
! Create a map 
  numGlobalElements_local = 4;
! map = epetra_map(numGlobalElements_local);

! numGlobalElements_return = map%NumGlobalElements()
  print *,'NumGlobalElements = ', numGlobalElements_return
  if ( numGlobalElements /= numGlobalElements_return ) &
    stop 'In ForTrilinos (verySimpleObjectOriented.F90: return mismatch'
   
  ! Create vectors
! x = epetra_vector(map)
! b = epetra_vector(map)
 
  ! Do some vector operations
! call b%PutScalar(two)
! call x%Update(two, b, zero) ! /* x = 2*b */
 
!  bnorm = b%Norm2()
!  xnorm = x%Norm2()
 
   print *, "2 norm of x = ", xnorm 
   print *, "2 norm of b = ", bnorm 

! Test the expected value 
 
   err_tol = 1.0e-14;
   expected_bnorm = sqrt( 2.0 * 2.0 * numGlobalElements_return );
   expected_xnorm = sqrt( 4.0 * 4.0 * numGlobalElements_return );
   bnorm_err = abs( expected_bnorm - bnorm ) / expected_bnorm;
   xnorm_err = abs( expected_xnorm - xnorm ) / expected_xnorm;
   print*, "error in 2 norm of x = ",bnorm_err
   print*, "error in 2 norm of b = ",xnorm_err
   if (bnorm_err > err_tol) success = .false.;
   if (xnorm_err > err_tol) success = .false.;
 
   ! Clean up memory (in reverse order)
!  call b%finalize()
!  call x%finalize()
!  call map%finalize()
 
   if (success) then
     print *  
     print *, "End Result: TEST PASSED" 
   else
     print *  
     print *, "End Result: TEST FAILED"
   end if
end program main
