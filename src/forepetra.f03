module forepetra
  use ,intrinsic :: iso_c_binding ,only : c_int,c_double
  implicit none

  !
  ! C procedure bindings
  !

  interface 
    integer(c_int) function FEpetra_Map_Create( numGlobalElements ) bind(c,name='Epetra_Map_Create')
      import :: c_int
      integer(c_int) ,value :: numGlobalElements
    end function 
  
    subroutine FEpetra_Map_Destroy( mapID ) bind(c,name='Epetra_Map_Destroy')
      import :: c_int
      integer(c_int) ,value :: mapID
    end subroutine 
  
    integer(c_int) function FEpetra_Map_NumGlobalElements( mapID ) bind(c,name='Epetra_Map_NumGlobalElements')
      import :: c_int
      integer(c_int) ,value :: mapID
    end function 
  end interface
  
  interface 
    integer(c_int) function FEpetra_Vector_Create( mapID ) bind(c,name='Epetra_Vector_Create')
      import :: c_int
      integer(c_int) ,value :: mapID
    end function
  
    subroutine FEpetra_Vector_Destroy( vectorID ) bind(c,name='Epetra_Vector_Destroy')
      import :: c_int
      integer(c_int) ,value :: vectorID
    end subroutine 
  
    subroutine FEpetra_Vector_PutScalar( vectorID, scalarConstant ) bind(c,name='Epetra_Vector_PutScalar')
      import :: c_int,c_double
      integer(c_int) ,value :: vectorID 
      real(c_double) ,value :: scalarConstant
    end subroutine
  
    subroutine FEpetra_Vector_Random( vectorID ) bind(c,name='Epetra_Vector_Random')
      import :: c_int
      integer(c_int) ,value :: vectorID 
    end subroutine 
  
    subroutine FEpetra_Vector_Update(vectorID, alpha, vector2ID, beta ) bind(c,name='Epetra_Vector_Update')
      import :: c_int,c_double
      integer(c_int) ,value :: vectorID ,vector2ID
      real(c_double) ,value :: alpha ,beta
    end subroutine
  
    real(c_double) function FEpetra_Vector_Norm2( vectorID ) bind(c,name='Epetra_Vector_Norm2')
      import :: c_int,c_double
      integer(c_int) ,value :: vectorID 
    end function 
  end interface
end module forepetra
