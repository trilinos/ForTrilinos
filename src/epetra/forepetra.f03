module forepetra  ! Companion to CEpetra_*.h
  use ,intrinsic     :: iso_c_binding ,only : c_int,c_double,c_char,c_bool ! Kind parameters 
  use ,non_intrinsic :: ForTrilinos_enums
  implicit none   ! Prevent implicit typing

  ! This file provides Fortran interface blocks that bind the argument types, return value types, and procedure names to those 
  ! in the C function prototypes in each of the CTrilinos/src/epetra/CEpetra_*.h header files.  The Fortran 2003 standard 
  ! guarantees that the types and names used in these bindings interoperate with a standard-conforming, companion C compiler.

  interface ! Since this file contains only interface bodies, this interface block ends at the bottom of the file.

  ! ____________________ Epetra_Object interface bodies _________________________________________________________________

  !CT_Epetra_Object_ID_t Epetra_Object_Cast(
  !  CTrilinos_Object_ID_t id );
  
  type(FT_Epetra_Object_ID_t) function Epetra_Object_Cast(id) bind(C,name='Epetra_Object_Cast')
    import :: FT_Epetra_Object_ID_t ,ForTrilinos_Object_ID_t
    type(ForTrilinos_Object_ID_t) ,value :: id 
  end function
  
  !/* Original C++ prototype:
  !   Epetra_Object(int TracebackModeIn = -1, bool set_label = true);
  !*/
  !CT_Epetra_Object_ID_t Epetra_Object_Create ( 
  !  int TracebackModeIn, boolean set_label );
  
  type(FT_Epetra_Object_ID_t) function Epetra_Object_Create(TracebackModeIn,set_label ) bind(C,name='Epetra_Object_Create')
   import :: FT_Epetra_Object_ID_t ,c_int ,c_bool
    integer(c_int)  ,value :: id 
    logical(c_bool) ,value :: set_label
  end function
  
  !CT_Epetra_Object_ID_t Epetra_Object_Create_WithLabel ( 
  !  const char * const Label, int TracebackModeIn );
  !
  !/* Original C++ prototype:
  !   Epetra_Object(const char * const Label, int TracebackModeIn = -1);
  !*/
  
  type(FT_Epetra_Object_ID_t) Epetra_Object_Create_WithLabel( Label, TracebackModeIn )bind(C,name='Epetra_Object_Create_WithLabel')
    import :: FT_Epetra_Object_ID_t ,c_int ,c_char
    character(kind=c_char) ,intent(in) ,dimension(*) :: Label 
    integer(c_int) ,value :: TracebackModeIn
  end function
 
  !/* Original C++ prototype:
  !   Epetra_Object(const Epetra_Object& Object);
  !*/
  !
  !CT_Epetra_Object_ID_t Epetra_Object_Duplicate ( 
  !  CT_Epetra_Object_ID_t ObjectID );
  !
  type(FT_Epetra_Object_ID_t) function Epetra_Object_Duplicate(FT_Epetra_Object_ID_t ObjectID)bind(C,name='Epetra_Object_Duplicate')
    import :: FT_Epetra_Object_ID_t
    type(FT_Epetra_Object_ID_t)  :: ObjectID
  end function 
 
  !/* Original C++ prototype:
  !   virtual ~Epetra_Object();
  !*/
  !void Epetra_Object_Destroy ( CT_Epetra_Object_ID_t * selfID );
  !
  subroutine Epetra_Object_Destroy ( selfID ) bind(C,name='Epetra_Object_Destroy')
    import :: FT_Epetra_Object_ID_t 
    type(FT_Epetra_Object_ID_t) :: selfID
  end subroutine
  
  !/* Original C++ prototype:
  !   virtual void SetLabel(const char * const Label);
  !*/
  !void Epetra_Object_SetLabel ( 
  !  CT_Epetra_Object_ID_t selfID, const char * const Label );
  
  subroutine Epetra_Object_SetLabel( selfID, Label ) bind(C,name='Epetra_Object_SetLabel')
    import :: FT_Epetra_Object_ID_t ,c_char
    type(FT_Epetra_Object_ID_t) ,value :: selfID
    character(kind=c_char) ,intent(in) ,dimension(*) :: Label 
  end subroutine
  
  !  CT_Epetra_Object_ID_t selfID, const char * const Label );
  !
  !/* Original C++ prototype:
  !   virtual const char * Label() const;
  !*/
  !const char * Epetra_Object_Label ( CT_Epetra_Object_ID_t selfID );
 
  function Epetra_Object_Label ( selfID ) bind(C,'Epetra_Object_Label')
    import :: FT_Epetra_Object_ID_t ,c_char
    character(kind=c_char) ,dimension(*) :: Epetra_Object_Label 
    type(FT_Epetra_Object_ID_t) :: selfID 
  end function
 
  ! ____________________ Epetra_Distributor interface bodies ________________________________________________________
  !CT_Epetra_Distributor_ID_t Epetra_Distributor_Cast(
  !  CTrilinos_Object_ID_t id );

   type(FT_Epetra_Distributor_ID_t) function Epetra_Distributor_Cast(id) bind(C,name='Epetra_Distributor_Cast')
    import :: FT_Epetra_Distributor_ID_t
    type(FTrilinos_Object_ID_t) ,value :: id 
   end function
 
  !
  !/* Original C++ prototype:
  !   virtual Epetra_Distributor * Clone() = 0;
  !*/
  !CT_Epetra_Distributor_ID_t Epetra_Distributor_Clone ( 
  !  CT_Epetra_Distributor_ID_t selfID );
 
  type(FT_Epetra_Distributor_ID_t) Epetra_Distributor_Clone( selfID ) bind(C,name='Epetra_Distributor_Clone')
    import :: FT_Epetra_Distributor_ID_t 
    type(CT_Epetra_Distributor_ID_t) ,vaue :: selfID 
   end function

  !/* Original C++ prototype:
  !  virtual ~Epetra_Distributor();
  !*/
  !void Epetra_Distributor_Destroy ( 
  ! CT_Epetra_Distributor_ID_t * selfID );

  subroutine Epetra_Distributor_Destroy ( selfID ) bind(C,name='Epetra_Distributor_Destroy')
    import :: FT_Epetra_Distributor_ID_t
    type(FT_Epetra_Distributor_ID_t)  :: selfID 
  end subroutine
 
  !/* Original C++ prototype:
  !  virtual int CreateFromSends( const int & NumExportIDs, const int * ExportPIDs, bool Deterministic, int & NumRemoteIDs ) = 0;
  !*/
  !int Epetra_Distributor_CreateFromSends ( 
  ! CT_Epetra_Distributor_ID_t selfID, int NumExportIDs, 
  ! const int * ExportPIDs, boolean Deterministic, 
  ! int * NumRemoteIDs );

  ! Assuming that the C prototype will be modified to 'const int * NumExportIDs'
 
  integer(c_int) function Epetra_Distributor_CreateFromSends ( selfID, NumExportIDs, ExportPIDs, Deterministic, NumRemoteIDs ) &
    bind(C,name='Epetra_Distributor_CreateFromSends')
    import :: c_int, FT_Epetra_Distributor_ID_t selfID, c_bool 
    type(FT_Epetra_Distributor_ID_t) ,value :: selfID
    integer(c_int) ,intent(in) :: NumExportIDs
    integer(c_int) ,intent(in) ,dimension(NumExportIDs) :: ExportPIDs
    logical(c_bool) ,value :: Deterministic
    integer(c_int) :: NumRemoteIDs
  end function
 
  !/* Original C++ prototype:
  !   virtual int CreateFromRecvs( const int & NumRemoteIDs, const int * RemoteGIDs, const int * RemotePIDs, 
  !   bool Deterministic, int & NumExportIDs, int *& ExportGIDs, int *& ExportPIDs) = 0;
  !*/
  ! int Epetra_Distributor_CreateFromRecvs ( 
  !   CT_Epetra_Distributor_ID_t selfID, int NumRemoteIDs, 
  !   const int * RemoteGIDs, const int * RemotePIDs, 
  !   boolean Deterministic, int * NumExportIDs, int ** ExportGIDs, 
  !   int ** ExportPIDs );
 
  ! Assuming that the C prototype will be modified to 'const int * NumRemoteIDs'

  integer(c_int) function Epetra_Distributor_CreateFromRecvs ( &
    selfID, NumRemoteIDs, RemoteGIDs, RemotePIDs, Deterministic, NumExportIDs, ExportGIDs, ExportPIDs &
    ) bind(C,name='Epetra_Distributor_CreateFromRecvs') 
    import :: CT_Epetra_Distributor_ID_t ,c_int ,c_bool
    type(CT_Epetra_Distributor_ID_t) ,value :: selfID
    integer(c_int) ,intent(in) :: NumRemoteIDs
    integer(c_int) ,intent(in) ,dimension(NumRemoteIDs) :: RemoteGIDs
    integer(c_int) ,intent(in) ,dimension(NumRemoteIDs) :: RemotePIDs, 
    logical(c_bool),value :: Deterministic
    integer(c_int) :: NumExportIDs
    integer(c_int) ,dimension(NumExportIDs) :: ExportGIDs, 
    integer(c_int) ,dimension(NumExportIDs) :: ExportPIDs 
  end function

  !/* Original C++ prototype:
  !    virtual int Do( char * export_objs, int obj_size, int & len_import_objs, char *& import_objs) = 0;
  !*/
  !int Epetra_Distributor_Do ( 
  !   CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  !   int * len_import_objs, char ** import_objs );

  integer(c_int) function Epetra_Distributor_Do ( selfID, export_objs, obj_size, len_import_objs, import_objs ) &
    bind(C,name='Epetra_Distributor_Do ')
    import :: FT_Epetra_Distributor_ID_t ,c_char ,c_int
    type(FT_Epetra_Distributor_ID_t) :: selfID
    character(kind=c_char) ,dimension(*) :: export_objs
    integer(c_int) ,value :: obj_size, 
    integer(c_int) :: len_import_objs
    charater(kind=c_char) ,dimension(len_import_objs) :: import_objs 
  end function

  !/* Original C++ prototype:
  !   virtual int DoReverse( char * export_objs, int obj_size, int & len_import_objs, char *& import_objs ) = 0;
  !*/
  !int Epetra_Distributor_DoReverse ( 
  !  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  !  int * len_import_objs, char ** import_objs );

  integer(c_int) function Epetra_Distributor_DoReverse ( selfID, export_objs, obj_size, len_import_objs, import_objs ) &
    bind(C,name='Epetra_Distributor_DoReverse') 
    import :: FT_Epetra_Distributor_ID_t ,c_char ,c_int
    type(FT_Epetra_Distributor_ID_t) ,value :: selfID
    character(kind=c_char) ,dimension(*) :: export_objs
    integer(c_int) ,value :: obj_size
    integer(c_int) :: len_import_objs
    character(kind=c_char) ,dimension(len_import_objs) :: import_objs
  end function

  !/* Original C++ prototype:
  !   virtual int DoPosts( char * export_objs, int obj_size, int & len_import_objs, char *& import_objs ) = 0;
  !*/
  !int Epetra_Distributor_DoPosts ( 
  ! CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  ! int * len_import_objs, char ** import_objs );

  integer(c_int) function Epetra_Distributor_DoPosts ( selfID, export_objs, obj_size, len_import_objs, import_objs ) &
    bind(C,name='Epetra_Distributor_DoPosts ')
    import :: FT_Epetra_Distributor_ID_t ,c_char , c_int  
    type(FT_Epetra_Distributor_ID_t) :: selfID
    character(kind=c_char) ,dimension(*) :: export_objs
    integer(c_int) ,value :: obj_size, 
    integer(c_int) :: len_import_objs
    character(c_char) ,dimension(len_import_objs) :: import_objs
  end function
     
  ! /* Original C++ prototype:
  ! virtual int DoWaits() = 0;
  ! */
  ! int Epetra_Distributor_DoWaits ( CT_Epetra_Distributor_ID_t selfID );

  integer(c_int) function Epetra_Distributor_DoWaits ( selfID ) bind(C,name='Epetra_Distributor_DoWaits ')
    import :: FT_Epetra_Distributor_ID_t ,c_int
    type(FT_Epetra_Distributor_ID_t) ,value :: selfID 
  end function
 
  !/* Original C++ prototype:
  !   virtual int DoReversePosts( char * export_objs, int obj_size, int & len_import_objs, char *& import_objs) = 0;
  !*/
  !int Epetra_Distributor_DoReversePosts ( 
  !  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  !  int * len_import_objs, char ** import_objs );

  integer(c_int) function Epetra_Distributor_DoReversePosts(selfID, export_objs, obj_size, len_import_objs, import_objs ) &
    bind(C,name='Epetra_Distributor_DoReversePosts') 
    import :: FT_Epetra_Distributor_ID_t , c_char ,c_int
    FT_Epetra_Distributor_ID_t ,value :: selfID
    character(kind=c_char) ,dimension(*) :: export_objs
    integer(c_int) ,value :: obj_size
    integer(c_int) :: len_import_objs
    character(kind=c_char) ,dimension(len_import_objs) :: import_objs 
  end function
   
  !/* Original C++ prototype:
  !virtual int DoReverseWaits() = 0;
  !*/
  ! int Epetra_Distributor_DoReverseWaits ( 
  !   CT_Epetra_Distributor_ID_t selfID );

  integer(c_int) function Epetra_Distributor_DoReverseWaits ( selfID ) bind(C,name='Epetra_Distributor_DoReverseWaits ')
    import :: FT_Epetra_Distributor_ID_t ,c_int
    type(FT_Epetra_Distributor_ID_t) ,value :: selfID
  end function

  !/* Original C++ prototype:
  !   virtual int Do( char * export_objs, int obj_size, int *& sizes, int & len_import_objs, char *& import_objs) = 0;
  !*/
  !int Epetra_Distributor_Do_VarLen ( 
  !  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  !  int ** sizes, int * len_import_objs, char ** import_objs );

  integer(c_int) function Epetra_Distributor_Do_VarLen ( selfID, export_objs, obj_size, sizes, len_import_objs, import_objs ) &
    bind(C,name='Epetra_Distributor_Do_VarLen') 
    import :: FT_Epetra_Distributor_ID_t ,c_int ,c_char 
    type(FT_Epetra_Distributor_ID_t) ,value :: selfID
    character(c_char) ,dimension(*) ::  export_objs
    integer(c_int) ,value :: obj_size
    integer(c_int) ,dimension(:) :: sizes
    integer(c_int) :: len_import_objs
    character(c_char) ,dimension(len_import_objs) :: import_objs 
  end function
  
  !/* Original C++ prototype:
  !   virtual int DoReverse( char * export_objs, int obj_size, int *& sizes, int & len_import_objs, char *& import_objs) = 0;
  !*/
  !int Epetra_Distributor_DoReverse_VarLen ( 
  !  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  !  int ** sizes, int * len_import_objs, char ** import_objs );

  integer(c_int) function Epetra_Distributor_DoReverse_VarLen(selfID, export_objs, obj_size, sizes, len_import_objs, import_objs)&
    bind(C,name='Epetra_Distributor_DoReverse_VarLen')  
    import ::  FT_Epetra_Distributor_ID_t, c_char , c_int  
    type(FT_Epetra_Distributor_ID_t) ,value :: selfID
    character(kind=c_char) ,dimension(*) ::  export_objs
    integer(c_int) ,value :: obj_size
    integer(c_int) ,dimension(*) :: sizes
    integer(c_int) ::  len_import_objs
    character(c_char) ,dimension(len_import_objs) :: import_objs 
  end function
     
  !/* Original C++ prototype:
  !  virtual int DoPosts( char * export_objs, int obj_size, int *& sizes, int & len_import_objs, char *& import_objs) = 0;
  !*/
  !int Epetra_Distributor_DoPosts_VarLen ( 
  !  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  !  int ** sizes, int * len_import_objs, char ** import_objs );
 
  integer(c_int) function Epetra_Distributor_DoPosts_VarLen( selfID, export_objs, obj_size, sizes, len_import_objs, import_objs ) &
    bind(C,name='Epetra_Distributor_DoPosts_VarLen')
    import :: c_int ,CT_Epetra_Distributor_ID_t ,c_char 
    type(CT_Epetra_Distributor_ID_t) ,value :: selfID
    character(c_char) ,dimension(*) :: export_objs
    integer(c_int) ,value :: obj_size 
    integer(c_int) ,dimension(*) :: sizes
    integer(c_int)  :: len_import_objs
    character(c_char) ,dimension(len_import_objs) :: import_objs
  end function
 
  !/* Original C++ prototype:
  !   virtual int DoReversePosts( char * export_objs, int obj_size, int *& sizes, int & len_import_objs, char *& import_objs) = 0;
  !*/
  !int Epetra_Distributor_DoReversePosts_VarLen ( 
  !  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  !  int ** sizes, int * len_import_objs, char ** import_objs );

  integer(c_int) Epetra_Distributor_DoReversePosts_VarLen ( 
    import :: c_int ,CT_Epetra_Distributor_ID_t ,c_int 
    type(CT_Epetra_Distributor_ID_t) ,value :: selfID
    character(c_char) ,dimension(*) :: export_objs
    integer(c_int) ,value :: obj_size
    integer(c_int) dimension(*) :: sizes
    integer(c_int) :: len_import_objs
    character(c_char) ,dimension(len_import_objs) :: import_objs 
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

