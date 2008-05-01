program main
  use ,intrinsic :: iso_c_binding ,only : c_int,c_double
  use Epetra_Map_Module
  use Epetra_Vector_Module
  !use :: forepetraext
  implicit none

 !
 ! Data declarations 
 !

  integer(c_int)   :: numGlobalElements_loc, numGlobalElements_rtn
  type(Epetra_Map) :: map
  type(Epetra_Vector)  :: x, b
  real(c_double) :: bnorm, xnorm,err_tol,expected_bnorm,expected_xnorm,bnorm_err,xnorm_err 
  real(c_double) :: two = 2.0, zero = 0.0
  logical        :: success = .true.
  
! /*
!  * Executable code
!  */
  
! /* Create a map */
  numGlobalElements_loc = 4;
  map = Create(numGlobalElements_loc);

  numGlobalElements_rtn = NumGlobalElements(map)
  print *,'NumGlobalElements = ', numGlobalElements_rtn
!! assert( numGlobalElements == numGlobalElements_rtn )
   
!! /* Create vectors */
   x = Create(map)
!  b = FEpetra_Vector::Create(map)
!
!! /* Do some vector operations */
!  call PutScalar(b, two)
!  call Update(x, two, b, zero) ! /* x = 2*b */
!
!  bnorm = Norm2(b)
!  xnorm = Norm2(x)
!
!  print *, "2 norm of x = ", xnorm 
!  print *, "2 norm of b = ", bnorm 
!! /* Test the expected value */
!
!  err_tol = 1e-14;
!  expected_bnorm = sqrt( 2.0 * 2.0 * numGlobalElements );
!  expected_xnorm = sqrt( 4.0 * 4.0 * numGlobalElements );
!  bnorm_err = abs( expected_bnorm - bnorm ) / expected_bnorm;
!  xnorm_err = abs( expected_xnorm - xnorm ) / expected_xnorm;
!  print*, "error in 2 norm of x = ",bnorm_err
!  print*, "error in 2 norm of b = ",xnorm_err
!  if (bnorm_err > err_tol) success = .false.;
!  if (xnorm_err > err_tol) success = .false.;
!
!! /* Clean up memory (in reverse order)! */
!  call Destroy(b)
!  call Destroy(x)
!  call Destroy(map)
!
!! /* This should throw an exception and print an error message! */
!! /* NumGlobalElements(map); */
!
!   if (success) then
!    print *  
!    print *, "End Result: TEST PASSED" 
!  else
!    print *  
!    print *, "End Result: TEST FAILED"
!  end if

end program main
