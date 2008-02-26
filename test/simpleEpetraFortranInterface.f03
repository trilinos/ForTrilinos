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
  real(c_double) :: bnorm, xnorm

! /*
!  * Executable code
!  */
  
! /* Create a map */
  numGlobalElements = 4;
  mapID = FEpetra_Map_Create(numGlobalElements);

  numGlobalElements_rtn = FEpetra_Map_NumGlobalElements(mapID)
  print *,'NumGlobalElements = ', numGlobalElements_rtn
! assert( numGlobalElements == numGlobalElements_rtn )
  
! /* Create vectors */
  xID = FEpetra_Vector_Create(mapID)
  bID = FEpetra_Vector_Create(mapID)

! /* Do some vector operations */
  call FEpetra_Vector_Random(bID)
  call FEpetra_Vector_Update(xID,2.0_c_double,bID,0.0_c_double) ! /* x = 2*b */

  bnorm = FEpetra_Vector_Norm2(bID)
  xnorm = FEpetra_Vector_Norm2(xID)

  print *, "2 norm of x = ", xnorm 
  print *, "2 norm of b = ", bnorm 

! /* Clean up memory (in reverse order)! */
  call FEpetra_Vector_Destroy(bID)
  call FEpetra_Vector_Destroy(xID)
  call FEpetra_Map_Destroy(mapID)

! /* This should throw an exception and print an error message! */
! /* FEpetra_Map_NumGlobalElements(mapID); */

end program main
