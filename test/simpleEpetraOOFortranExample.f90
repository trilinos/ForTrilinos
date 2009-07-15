program main
  use ,intrinsic :: iso_c_binding ,only : c_int,c_double
  use Epetra_Map_Module           ,only : EPetra_Map
  use Epetra_Vector_Module        ,only : EPetra_Vector
  implicit none

 !
 ! Data declarations 
 !
  type(Epetra_Map)    :: map
  type(Epetra_Vector) :: x, b
  integer(c_int) :: numGlobalElements_local, numGlobalElements_return
  real(c_double) :: bnorm, xnorm,err_tol,expected_bnorm,expected_xnorm,bnorm_err,xnorm_err 
  real(c_double) :: two = 2.0, zero = 0.0
  logical        :: success = .true.
! 
! Executable code
!  
! Create a map 
  numGlobalElements_local = 4;
  call map%Create(numGlobalElements_local);

  numGlobalElements_return = map%NumGlobalElements()
  print *,'NumGlobalElements = ', numGlobalElements_return
!! assert( numGlobalElements == numGlobalElements_return )
   
!! Create vectors
   call x%Create(map)
   call b%Create(map)
!
! Do some vector operations
   call b%PutScalar(two)
   call x%Update(two, b, zero) ! /* x = 2*b */
 
   bnorm = b%Norm2()
   xnorm = x%Norm2()
 
   print *, "2 norm of x = ", xnorm 
   print *, "2 norm of b = ", bnorm 
! Test the expected value 
 
   err_tol = 1e-14;
   expected_bnorm = sqrt( 2.0 * 2.0 * numGlobalElements_return );
   expected_xnorm = sqrt( 4.0 * 4.0 * numGlobalElements_return );
   bnorm_err = abs( expected_bnorm - bnorm ) / expected_bnorm;
   xnorm_err = abs( expected_xnorm - xnorm ) / expected_xnorm;
   print*, "error in 2 norm of x = ",bnorm_err
   print*, "error in 2 norm of b = ",xnorm_err
   if (bnorm_err > err_tol) success = .false.;
   if (xnorm_err > err_tol) success = .false.;
!
!! /* Clean up memory (in reverse order)! */
   call b%Destroy()
   call x%Destroy()
  !call Destroy(map)
!
! This should throw an exception and print an error message!
! call NumGlobalElements(map)

   if (success) then
     print *  
     print *, "End Result: TEST PASSED" 
   else
     print *  
     print *, "End Result: TEST FAILED"
   end if
end program main
