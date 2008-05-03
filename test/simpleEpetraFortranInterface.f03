program main
  use ,intrinsic :: iso_c_binding ,only : c_int,c_double
  use            :: forepetraext
  implicit none

 !
 ! Data declarations 
 !

  integer(c_int) :: numGlobalElements, numGlobalElements_rtn
  integer(c_int) :: mapID    ! TYPE(MapID)
  integer(c_int) :: xID, bID ! TYPE(VectorID)
  real(c_double) :: bnorm, xnorm,err_tol,expected_bnorm,expected_xnorm,bnorm_err,xnorm_err 
  real(c_double) :: two = 2.0, zero = 0.0
  logical        :: success = .true.
  
! /*
!  * Executable code
!  */
  
! /* Create a map */
  numGlobalElements = 4;
  mapID = Epetra_Map_Create(numGlobalElements);

  numGlobalElements_rtn = Epetra_Map_NumGlobalElements(mapID)
  print *,'NumGlobalElements = ', numGlobalElements_rtn
! assert( numGlobalElements == numGlobalElements_rtn )
  
! /* Create vectors */
  xID = Epetra_Vector_Create(mapID)
  bID = Epetra_Vector_Create(mapID)

! /* Do some vector operations */
  call Epetra_Vector_PutScalar(bID, two)
  call Epetra_Vector_Update(xID, two, bID, zero) ! /* x = 2*b */

  bnorm = Epetra_Vector_Norm2(bID)
  xnorm = Epetra_Vector_Norm2(xID)

  print *, "2 norm of x = ", xnorm 
  print *, "2 norm of b = ", bnorm 
! /* Test the expected value */

  err_tol = 1e-14;
  expected_bnorm = sqrt( 2.0 * 2.0 * numGlobalElements );
  expected_xnorm = sqrt( 4.0 * 4.0 * numGlobalElements );
  bnorm_err = abs( expected_bnorm - bnorm ) / expected_bnorm;
  xnorm_err = abs( expected_xnorm - xnorm ) / expected_xnorm;
  print*, "error in 2 norm of x = ",bnorm_err
  print*, "error in 2 norm of b = ",xnorm_err
  if (bnorm_err > err_tol) success = .false.;
  if (xnorm_err > err_tol) success = .false.;

! /* Clean up memory (in reverse order)! */
  call Epetra_Vector_Destroy(bID)
  call Epetra_Vector_Destroy(xID)
  call Epetra_Map_Destroy(mapID)

! /* This should throw an exception and print an error message! */
! /* Epetra_Map_NumGlobalElements(mapID); */
   if (success) then
    print *  
    print *, "End Result: TEST PASSED" 
  else
    print *  
    print *, "End Result: TEST FAILED"
  end if

end program main
