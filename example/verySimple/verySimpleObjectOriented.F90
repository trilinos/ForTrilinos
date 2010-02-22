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

 
  use ,intrinsic :: iso_c_binding ,only : c_int,c_double
  use FEpetra_SerialComm   ,only : epetra_serialcomm
  use FEpetra_Map          ,only : epetra_map
  use FEpetra_Vector       ,only : epetra_vector
  use ForTrilinos_utils    ,only : valid_kind_parameters
  implicit none

  ! Data declarations 
  
  type(epetra_serialcomm) :: communicator
  type(epetra_map)    :: map
  type(epetra_vector) :: x, b
  integer(c_int) :: numGlobalElements_local, numGlobalElements_return
  integer(c_int) :: Index_Base=1
  real(c_double) ,allocatable ,dimension(:) :: bnorm, xnorm
  real(c_double) ,allocatable ,dimension(:) :: err_tol,expected_bnorm,expected_xnorm,bnorm_err,xnorm_err 
  real(c_double) :: two = 2.0, zero = 0.0
  logical        :: success = .true.,zero_initial=.true.
  
  if (.not. valid_kind_parameters()) stop 'C interoperability not supported on this platform.'
  
  ! Executable code
  
! Create a serial comm
  communicator= epetra_serialcomm() 

! Create a map 
  numGlobalElements_local = 4
  map = epetra_map(numGlobalElements_local,Index_Base,communicator)
  numGlobalElements_return = map%NumGlobalElements()
  print *,'NumGlobalElements = ', numGlobalElements_return
  if ( numGlobalElements_local /= numGlobalElements_return ) &
    stop 'In ForTrilinos (verySimpleObjectOriented.F90: return mismatch'
   
  ! Create vectors
  x = epetra_vector(map,zero_initial)
  b = epetra_vector(map,zero_initial)
 
  ! Do some vector operations
  call b%Random()
  call x%Update(two, b, zero) ! /* x = 2*b */
 
  bnorm = b%Norm1()
  xnorm = x%Norm1()
 
   print *, "2 norm of x = ", xnorm 
   print *, "2 norm of b = ", bnorm 

! Test the expected value 
   err_tol = 1.0e-14
   expected_bnorm = [sqrt( 2.0 * 2.0 * numGlobalElements_return )]
   expected_xnorm = [sqrt( 4.0 * 4.0 * numGlobalElements_return )]
   bnorm_err = abs( expected_bnorm - bnorm ) / expected_bnorm
   xnorm_err = abs( expected_xnorm - xnorm ) / expected_xnorm
   print*, "error in 2 norm of x = ",bnorm_err
   print*, "error in 2 norm of b = ",xnorm_err
   if (any(bnorm_err > err_tol)) success = .false.
   if (any(xnorm_err > err_tol)) success = .false.
 
  ! Clean up memory (in reverse order).  This step is not required
  ! with compilers that fupport Fortran 2003 type finalization:
  call b%force_finalization()
  call x%force_finalization()
  call map%force_finalization()
  call communicator%force_finalization()
 
   if (success) then
     print *  
     print *, "End Result: TEST PASSED" 
   else
     print *  
     print *, "End Result: TEST FAILED"
   end if
end program main
