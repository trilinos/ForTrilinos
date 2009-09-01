program main
  use iso_c_binding ,only : c_int,c_double,c_bool
  use iso_fortran_env ,only : error_unit ,output_unit
  use fortrilinos_utils ,only : valid_kind_parameters
  use forepetra 

  ! This file is the Fortran equivalent of CTrilionos/example/verySimple.c.
  ! As such, it makes direct use of the procedural bindings in forepetra.F90.
  ! This style of accessing Trilinos is recommended for Fortran users whose
  ! compilers do not support the object-oriented features of Fortran 2003.
  
  !/*
  ! * Data declarations 
  ! */

  integer(c_int) numGlobalElements, numGlobalElements_rtn ;
  integer(c_int) junk;

  type(FT_Epetra_SerialComm_ID_t) scommID;
  type(FT_Epetra_Comm_ID_t) commID;

  type(FT_Epetra_Map_ID_t) mapID;
  type(FT_Epetra_BlockMap_ID_t) bmapID;

  type(FT_Epetra_Vector_ID_t) xID, bID;
  type(FT_Epetra_MultiVector_ID_t) mxID, mbID;

  real(c_double) bnorm, xnorm, expected_bnorm, expected_xnorm, bnorm_err, xnorm_err ;
  real(c_double) ,parameter :: err_tol=0.0

  logical(c_bool) :: success = .TRUE.;

  !/*
  ! * Executable code
  ! */

  if (.not. valid_kind_parameters()) &
     stop 'In ForTrilinos (verySimple.F90): C interoperability not supported on this platform.'
  
  !/* Create an Epetra_SerialComm and cast to an Epetra_Comm so that
  ! * it can be passed to functions expecting the latter */

  scommID = Epetra_SerialComm_Create();
  commID = Epetra_Comm_Cast(Epetra_SerialComm_Abstract(scommID));

  !/* Create an Epetra_Map and cast to an Epetra_BlockMap so that
  ! * a) it can be passed to functions expecting the latter and
  ! * b) methods implemented only in BlockMap can be invoked on the Map */

  numGlobalElements = 4;
  !/* use indexBase = 0 unless you know what you're doing! */
  mapID = Epetra_Map_Create(numGlobalElements, indexBase=0, CommID=commID);
  bmapID = Epetra_BlockMap_Cast(Epetra_Map_Abstract(mapID));

  !/* Check the properties of the map */
  numGlobalElements_rtn = Epetra_BlockMap_NumGlobalElements(bmapID);
  print *, "NumGlobalElements = ", numGlobalElements_rtn ;
  if ( numGlobalElements /= numGlobalElements_rtn ) stop 'main: incorrect numGlobalElements_rtn';

  print *,"NumMyElements = ", Epetra_BlockMap_NumMyElements(bmapID);
  
  !/* Create an Epetra_Vector and cast to an Epetra_MultiVector so that
  ! * methods implemented only in MultiVector can be invoked on the Vector */
  xID = Epetra_Vector_Create(bmapID, .TRUE._c_bool); !/* zero this one */
  mxID = Epetra_MultiVector_Cast(Epetra_Vector_Abstract(xID));

  !/* Do the same thing, but do not initialize this one to zero */
  bID = Epetra_Vector_Create(bmapID, .FALSE._c_bool);
  mbID = Epetra_MultiVector_Cast(Epetra_Vector_Abstract(bID));

  !/* Do some vector operations */
  junk = Epetra_MultiVector_PutScalar(mbID, 2.0_c_double);
  junk = Epetra_MultiVector_Update_WithA(mxID, 2.0_c_double, mbID, 0.0_c_double);!/* x = 2*b */

  junk = Epetra_MultiVector_Norm2(mbID, bnorm);
  junk = Epetra_MultiVector_Norm2(mxID, xnorm);

  print *,"2 norm of x = ", xnorm;
  print *,"2 norm of b = ", bnorm;

  !/* Test the expected value */
  expected_bnorm = sqrt( 2.0 * 2.0 * numGlobalElements );
  expected_xnorm = sqrt( 4.0 * 4.0 * numGlobalElements );
  bnorm_err = fabs( expected_bnorm - bnorm ) / expected_bnorm;
  xnorm_err = fabs( expected_xnorm - xnorm ) / expected_xnorm;
  print *,"error in 2 norm of x = ", bnorm_err ;
  print *,"error in 2 norm of b = ", xnorm_err ;

  if (bnorm_err > err_tol) success = .FALSE.;
  if (xnorm_err > err_tol) success = .FALSE.;

  !/* Clean up memory (in reverse order)! */
  call Epetra_MultiVector_Destroy(mxID);
  call Epetra_MultiVector_Destroy(mbID);
  call Epetra_Vector_Destroy(bID);
  call Epetra_Vector_Destroy(xID);

  call Epetra_BlockMap_Destroy(bmapID);
  call Epetra_Map_Destroy(mapID);

  call Epetra_Comm_Destroy(commID);
  call Epetra_SerialComm_Destroy(scommID);

  !/* This should throw an exception and print an error message
  ! * since the object has already been destroyed! */
  !/* Epetra_BlockMap_NumGlobalElements(bmapID); */

  if (success) then
    write(output_unit,*) 
    write(output_unit,fmt='(a)') "End Result: TEST PASSED" ;
  else
    write(error_unit,*) 
    write(error_unit,fmt='(a)') "End Result: TEST FAILED" ;
  end if
  
  ! return ( (success == .TRUE.) ? 0 : 1 );

end program
