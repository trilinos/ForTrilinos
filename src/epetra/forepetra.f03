module forepetra  ! Companion to CEpetra_*.h
  use ,intrinsic    :: iso_c_binding ,only : c_int,c_double,c_char,c_bool ! Kind parameters 
  use ,nonintrinsic :: ForTrilinos_enums
  implicit none   ! Prevent implicit typing

  ! This file provides Fortran interface blocks that bind the argument types, return value types, and procedure names to those 
  ! in the C function prototypes in each of the CTrilinos/src/epetra/CEpetra_*.h header files.  The Fortran 2003 standard 
  ! guarantees that the types and names used in these bindings interoperate with a standard-conforming, companion C compiler.

  interface ! Since this file contains only interface bodies, this interface block ends at the bottom of the file.

  ! ____________________ Epetra_Data interface bodies ___________________________________________________________________
  ! C prototype: CT_Epetra_Data_ID_t Epetra_Data_Cast( CTrilinos_Object_ID_t id );
  
  type(FT_Epetra_Data_ID_t) function epetra_data_cast( id ) bind(C,name='Epetra_Data_Cast')
    type(ForTrilinos_Object_ID_t), value :: id 
  end function 

  ! ____________________ Epetra_SerialCommData interface bodies _________________________________________________________
  ! C prototype: CT_Epetra_SerialCommData_ID_t Epetra_SerialCommData_Cast( CTrilinos_Object_ID_t id );
  
  type(FT_Epetra_SerialCommData_ID_t) function epetra_serialCommData_Cast( id ) bind(C,name='Epetra_SerialCommData_Cast')
    type(ForTrilinos_Object_ID_t) ,value :: id 
  end function

  ! ____________________ Epetra_MpiCommData interfaces bodies ___________________________________________________________
  ! C prototype: CT_Epetra_MpiCommData_ID_t Epetra_MpiCommData_Cast( CTrilinos_Object_ID_t id );

  type(FT_Epetra_MpiCommData_ID_t) function epetra_mpiCommData_cast( id ) bind(C,name='Epetra_MpiCommDoata_Cast')
    type(ForTrilinos_Object_ID_t) ,value :: id 
  end function
   
  ! ____________________ Epetra_BlockMapData interface bodies ___________________________________________________________
  ! C prototype: CT_Epetra_BlockMapData_ID_t  Epetra_BlockMapData_Cast( CTrilinos_Object_ID_t id );

  type(FT_Epetra_BlockMapData_ID_t ) function epetra_blockMapData_cast(id) bind(C,name='Epetra_BlockMapData_Cast')
    type(ForTrilinos_Object_ID_t) ,value :: id 
  end function

  ! ____________________ Epetra_CrsGraphData interface bodies ___________________________________________________________
  ! C prototype: CT_Epetra_CrsGraphData_ID_t  Epetra_CrsGraphData_Cast( CTrilinos_Object_ID_t id ); 

  type(FT_Epetra_CrsGraphData_ID_t ) function epetra_crsGraphData_cast( id ) bind(C,name='Epetra_CrsGraphData_Cast')
    type(ForTrilinos_Object_ID_t) ,value :: id 
  end function
   
  ! ____________________ Epetra_Object interface bodies _________________________________________________________________
  ! C prototype: CT_Epetra_Object_ID_t Epetra_Object_Cast( CTrilinos_Object_ID_t id );

  type(FT_Epetra_Object_ID_t) function epetra_object_cast( id ) bind(C,name='Epetra_Object_Cast')
    type(ForTrilinos_Object_ID_t) ,value :: id 
  end function

  ! C prototype: CT_Epetra_Object_ID_t Epetra_Object_Create ( int TracebackModeIn, boolean set_label );

  type(FT_Epetra_Object_ID_t) function epetra_object_create(tracebackModeIn,set_label) bind(C,name='Epetra_Object_Create')
    import :: c_int,c_bool
    integer(c_int)  ,value :: tracebackModeIn
    logical(c_bool) ,value :: set_label 
  end function

  ! C prototype: CT_Epetra_Object_ID_t Epetra_Object_Create_WithLabel ( const char * const Label, int TracebackModeIn );

  type(FT_Epetra_Object_ID_t) function epetra_object_create_withLabel( label, tracebackModeIn ) &
    bind(C,name='Epetra_Object_Create_WithLabel')
    import :: c_int,c_char
    integer(c_char) ,intent(in) :: label
    integer(c_int ) ,value      :: tracebackModeIn
  end function

  ! C prototype: CT_Epetra_Object_ID_t Epetra_Object_Duplicate ( CT_Epetra_Object_ID_t ObjectID );

  type(FT_Epetra_Object_ID_t) function epetra_object_duplicate( objectID ) bind(C,name='Epetra_Object_Duplicate')
    type(FT_Epetra_Object_ID_t) ,value :: objectID 
  end function

  ! C prototype: void Epetra_Object_Destroy ( CT_Epetra_Object_ID_t * selfID );

  subroutine epetra_object_destroy( selfID ) bind(C,name='Epetra_Object_Destroy')
    type(ForTrilinos_Object_ID_t) ,value :: selfID 
  end function

  ! C prototype: void Epetra_Object_SetLabel ( CT_Epetra_Object_ID_t selfID, const char * const Label );

  subroutine epetra_object_setLabel( selfID ,label ) bind(C,name='Epetra_Object_SetLabel')
    import :: c_char
    type(FT_Epetra_Object_ID_t) ,value :: selfID 
    character(c_char) ,intent(in)      :: label 
  end function

  ! C prototype: const char * Epetra_Object_Label ( CT_Epetra_Object_ID_t selfID );

  character(c_char) function epetra_object_label ( selfID ) ,bind(C,name='Epetra_Object_Label')
    import :: c_char
    type(FT_Epetra_Object_ID_t) :: selfID
  end function

  ! ____________________ Epetra_Distributor interface bodies ________________________________________________________

  ! C prototype: CT_Epetra_Distributor_ID_t Epetra_Distributor_Cast( CTrilinos_Object_ID_t id );

  type(FT_Epetra_Distributor_ID_t) function epetra_distributor_cast( id ) ,bind(C,name='Epetra_Distributor_Cast')
    type(ForTrilinos_Object_ID_t) ,value :: id
  end function

  ! C prototype: CT_Epetra_Distributor_ID_t Epetra_Distributor_Clone (  CT_Epetra_Distributor_ID_t selfID );

  type(FT_Epetra_Distributor_ID_t) function epetra_distributor_clone ( selfID ) ,bind(C,name='Epetra_Distributor_Clone ')
    type(FT_Epetra_Distributor_ID_t) ,value :: selfID 
  end function

  ! C prototype: void Epetra_Distributor_Destroy ( CT_Epetra_Distributor_ID_t * selfID );

  subroutine epetra_distributor_destroy ( selfID ) ,bind(C,name='Epetra_Distributor_Destroy ')
    type(FT_Epetra_Distributor_ID_t) :: selfID 
  end subroutine

  ! C prototype: int Epetra_Distributor_CreateFromSends ( 
  !                CT_Epetra_Distributor_ID_t selfID,  const int * NumExportIDs, 
  !                const int * ExportPIDs, boolean Deterministic, 
  !                int * NumRemoteIDs );

  integer(c_int) function epetra_distributor_createFromSends( selfID,numExportIDs,exportPIDs,deterministic,numRemoteIDs ) &
    bind(C,name='Epetra_Distributor_CreateFromSends')
    type(FT_Epetra_Distributor_ID_t) ,value ::  selfID
    integer(c_int) ,intent(in)              :: numExportIDs 
    integer(c_int) ,intent(in)              :: exportPIDs
    logical(c_ bool) ,value                 :: deterministic 
    integer(c_int)                          :: NumRemoteIDs 
   end function

  ! C prototype: int Epetra_Distributor_CreateFromRecvs ( 
  !                CT_Epetra_Distributor_ID_t selfID,  const int * NumRemoteIDs, 
  !                const int * RemoteGIDs, const int * RemotePIDs, 
  !                boolean Deterministic, int * NumExportIDs, int ** ExportGIDs, 
  !                int ** ExportPIDs );

  integer(c_int) function epetra_distributor_createFromRecvs ( ) &
    bind(C,name='Epetra_Distributor_CreateFromRecvs ')
    type(CT_Epetra_Distributor_ID_t) ,value :: selfID
    integer(c_int) ,intent(in)              :: numRemoteIDs
    integer(c_int) ,intent(in)              :: remoteGIDs
    integer(c_int) ,intent(in)              :: remotePIDs, 
    logical(c_bool)                  ,value :: deterministic
    integer(c_int) ,intent(in)              :: numExportIDs
    integer(c_int) ,dimension(*)            :: exportGIDs 
    integer(c_int) ,dimension(*)            :: exportPIDs 
  end function

  ! C prototype: int Epetra_Distributor_Do ( 
  !                CT_Epetra_Distributor_ID_t selfID,  char * export_objs, 
  !                int obj_size, int * len_import_objs, char ** import_objs );

  ! C prototype: int Epetra_Distributor_DoReverse ( 
  !                CT_Epetra_Distributor_ID_t selfID,  char * export_objs, 
  !                int obj_size, int * len_import_objs, char ** import_objs );

  ! C prototype: int Epetra_Distributor_DoPosts ( 
  !                CT_Epetra_Distributor_ID_t selfID,  char * export_objs, 
  !                int obj_size, int * len_import_objs, char ** import_objs );

  ! C prototype: int Epetra_Distributor_DoWaits ( CT_Epetra_Distributor_ID_t selfID );

  ! C prototype: int Epetra_Distributor_DoReversePosts ( 
  !                CT_Epetra_Distributor_ID_t selfID,  char * export_objs, 
  !                int obj_size, int * len_import_objs, char ** import_objs );

  ! C prototype: int Epetra_Distributor_DoReverseWaits ( CT_Epetra_Distributor_ID_t selfID );
  

  ! C prototype: int Epetra_Distributor_Do_VarLen ( 
  !                CT_Epetra_Distributor_ID_t selfID,  char * export_objs, 
  !                int obj_size, int ** sizes, int * len_import_objs, 
  !                char ** import_objs );

  ! C prototype: int Epetra_Distributor_DoReverse_VarLen ( 
  !                CT_Epetra_Distributor_ID_t selfID,  char * export_objs, 
  !                int obj_size, int ** sizes, int * len_import_objs, 
  !                char ** import_objs );

  ! C prototype: int Epetra_Distributor_DoPosts_VarLen ( 
  !                CT_Epetra_Distributor_ID_t selfID,  char * export_objs, 
  !                int obj_size, int ** sizes, int * len_import_objs, 
  !                char ** import_objs );

  ! C prototype: int Epetra_Distributor_DoReversePosts_VarLen ( 
  !                CT_Epetra_Distributor_ID_t selfID,  char * export_objs, 
  !                int obj_size, int ** sizes, int * len_import_objs, 
  !                char ** import_objs );

   
  ! ____________________ Epetra_Map bindings ________________________
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
  
  ! ___________________ Epetra_Vector bindings ____________________

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

  ! ___________________ Epetra_SerialComm bindings ____________________
      integer(c_int) function FEpetra_SerialComm_Create() bind(c,name='Epetra_SerialComm_Create')
        import :: c_int
      end function

      subroutine FEpetra_SerialComm_Destroy(id) bind(c,name='Epetra_SerialComm_Destroy')
        import :: c_int
        integer(c_int) ,value :: id
      end subroutine

  ! ___________________ Epetra_MpiComm bindings ____________________
      integer(c_int) function FEpetra_MpiComm_Create() bind(c,name='Epetra_MpiComm_Create')
        import :: c_int
      end function

      subroutine FEpetra_MpiComm_Destroy(id) bind(c,name='Epetra_MpiComm_Destroy')
        import :: c_int
        integer(c_int) ,value :: id
      end subroutine
  
    end interface

end module forepetra
