module forepetra  ! Companion to CEpetra_*.h
  use ,intrinsic     :: iso_c_binding ,only : c_int,c_double,c_char,c_bool,c_ptr ,c_long ! Kind parameters 
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
    implicit none
    type(ForTrilinos_Object_ID_t) ,value :: id 
  end function
  
  !/* Original C++ prototype:
  !   Epetra_Object(int TracebackModeIn = -1, bool set_label = true);
  !*/
  !CT_Epetra_Object_ID_t Epetra_Object_Create ( 
  !  int TracebackModeIn, boolean set_label );
  
  type(FT_Epetra_Object_ID_t) function Epetra_Object_Create(TracebackModeIn,set_label ) bind(C,name='Epetra_Object_Create')
    import :: FT_Epetra_Object_ID_t ,c_int ,c_bool
    implicit none
    integer(c_int)  ,value :: TracebackModeIn
    logical(c_bool) ,value :: set_label
  end function
  
  !CT_Epetra_Object_ID_t Epetra_Object_Create_WithLabel ( 
  !  const char * const Label, int TracebackModeIn );
  !
  !/* Original C++ prototype:
  !   Epetra_Object(const char * const Label, int TracebackModeIn = -1);
  !*/
  
  type(FT_Epetra_Object_ID_t) function Epetra_Object_Create_WithLabel(Label,TracebackModeIn) &
    bind(C,name='Epetra_Object_Create_WithLabel')
    import :: FT_Epetra_Object_ID_t ,c_int ,c_char
    implicit none
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
  type(FT_Epetra_Object_ID_t) function Epetra_Object_Duplicate(ObjectID) bind(C,name='Epetra_Object_Duplicate')
    import :: FT_Epetra_Object_ID_t
    implicit none
    type(FT_Epetra_Object_ID_t)  :: ObjectID
  end function 
 
  !/* Original C++ prototype:
  !   virtual ~Epetra_Object();
  !*/
  !void Epetra_Object_Destroy ( CT_Epetra_Object_ID_t * selfID );
  !
  subroutine Epetra_Object_Destroy ( selfID ) bind(C,name='Epetra_Object_Destroy')
    import :: FT_Epetra_Object_ID_t 
    implicit none
    type(FT_Epetra_Object_ID_t) :: selfID
  end subroutine
  
  !/* Original C++ prototype:
  !   virtual void SetLabel(const char * const Label);
  !*/
  !void Epetra_Object_SetLabel ( 
  !  CT_Epetra_Object_ID_t selfID, const char * const Label );
  
  subroutine Epetra_Object_SetLabel( selfID, Label ) bind(C,name='Epetra_Object_SetLabel')
    import :: FT_Epetra_Object_ID_t ,c_char
    implicit none
    type(FT_Epetra_Object_ID_t) ,value :: selfID
    character(kind=c_char) ,intent(in) ,dimension(*) :: Label 
  end subroutine
  
  !  CT_Epetra_Object_ID_t selfID, const char * const Label );
  !
  !/* Original C++ prototype:
  !   virtual const char * Label() const;
  !*/
  !const char * Epetra_Object_Label ( CT_Epetra_Object_ID_t selfID );
 
  type(c_ptr) function Epetra_Object_Label ( selfID ) bind(C,name='Epetra_Object_Label')
    import :: FT_Epetra_Object_ID_t ,c_ptr
    implicit none
    type(FT_Epetra_Object_ID_t) :: selfID 
  end function
 
  ! ____________________ Epetra_Distributor interface bodies ________________________________________________________

  !CT_Epetra_Distributor_ID_t Epetra_Distributor_Cast(
  !  CTrilinos_Object_ID_t id );

   type(FT_Epetra_Distributor_ID_t) function Epetra_Distributor_Cast(id) bind(C,name='Epetra_Distributor_Cast')
    import :: FT_Epetra_Distributor_ID_t ,ForTrilinos_Object_ID_t
    implicit none
    type(ForTrilinos_Object_ID_t) ,value :: id 
   end function
 
  !
  !/* Original C++ prototype:
  !   virtual Epetra_Distributor * Clone() = 0;
  !*/
  !CT_Epetra_Distributor_ID_t Epetra_Distributor_Clone ( 
  !  CT_Epetra_Distributor_ID_t selfID );
 
  type(FT_Epetra_Distributor_ID_t) function Epetra_Distributor_Clone( selfID ) bind(C,name='Epetra_Distributor_Clone')
    import :: FT_Epetra_Distributor_ID_t 
    implicit none
    type(FT_Epetra_Distributor_ID_t) ,value :: selfID 
   end function

  !/* Original C++ prototype:
  !  virtual ~Epetra_Distributor();
  !*/
  !void Epetra_Distributor_Destroy ( 
  ! CT_Epetra_Distributor_ID_t * selfID );

  subroutine Epetra_Distributor_Destroy ( selfID ) bind(C,name='Epetra_Distributor_Destroy')
    import :: FT_Epetra_Distributor_ID_t
    implicit none
    type(FT_Epetra_Distributor_ID_t)  :: selfID 
  end subroutine
 
  !/* Original C++ prototype:
  !  virtual int CreateFromSends( const int & NumExportIDs, const int * ExportPIDs, bool Deterministic, int & NumRemoteIDs ) = 0;
  !*/
  !int Epetra_Distributor_CreateFromSends ( 
  ! CT_Epetra_Distributor_ID_t selfID, int NumExportIDs, 
  ! const int * ExportPIDs, boolean Deterministic, 
  ! int * NumRemoteIDs );

  integer(c_int) function Epetra_Distributor_CreateFromSends ( selfID, NumExportIDs, ExportPIDs, Deterministic, NumRemoteIDs ) &
    bind(C,name='Epetra_Distributor_CreateFromSends')
    import :: c_int ,FT_Epetra_Distributor_ID_t ,c_bool 
    implicit none
    type(FT_Epetra_Distributor_ID_t) ,value :: selfID
    integer(c_int) ,intent(in) ,value :: NumExportIDs
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
 
  integer(c_int) function Epetra_Distributor_CreateFromRecvs ( &
    selfID, NumRemoteIDs, RemoteGIDs, RemotePIDs, Deterministic, NumExportIDs, ExportGIDs, ExportPIDs &
    ) bind(C,name='Epetra_Distributor_CreateFromRecvs') 
    import :: FT_Epetra_Distributor_ID_t ,c_int ,c_bool
    implicit none
    type(FT_Epetra_Distributor_ID_t) ,value :: selfID
    integer(c_int) ,intent(in) ,value :: NumRemoteIDs
    integer(c_int) ,intent(in) ,dimension(NumRemoteIDs) :: RemoteGIDs
    integer(c_int) ,intent(in) ,dimension(NumRemoteIDs) :: RemotePIDs
    logical(c_bool),value :: Deterministic
    integer(c_int) :: NumExportIDs
    integer(c_int) ,dimension(NumExportIDs) :: ExportGIDs
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
    implicit none
    type(FT_Epetra_Distributor_ID_t) :: selfID
    character(kind=c_char) ,dimension(*) :: export_objs
    integer(c_int) ,value :: obj_size
    integer(c_int) :: len_import_objs
    character(kind=c_char) ,dimension(len_import_objs) :: import_objs 
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
    implicit none
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
    implicit none
    type(FT_Epetra_Distributor_ID_t) :: selfID
    character(kind=c_char) ,dimension(*) :: export_objs
    integer(c_int) ,value :: obj_size
    integer(c_int) :: len_import_objs
    character(kind=c_char) ,dimension(len_import_objs) :: import_objs
  end function
     
  ! /* Original C++ prototype:
  ! virtual int DoWaits() = 0;
  ! */
  ! int Epetra_Distributor_DoWaits ( CT_Epetra_Distributor_ID_t selfID );

  integer(c_int) function Epetra_Distributor_DoWaits ( selfID ) bind(C,name='Epetra_Distributor_DoWaits ')
    import :: FT_Epetra_Distributor_ID_t ,c_int
    implicit none
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
    implicit none
    type(FT_Epetra_Distributor_ID_t) ,value :: selfID
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
    implicit none
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
    implicit none
    type(FT_Epetra_Distributor_ID_t) ,value :: selfID
    character(kind=c_char) ,dimension(*) ::  export_objs
    integer(c_int) ,value :: obj_size
    integer(c_int) ,dimension(*) :: sizes
    integer(c_int) :: len_import_objs
    character(kind=c_char) ,dimension(len_import_objs) :: import_objs 
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
    implicit none
    type(FT_Epetra_Distributor_ID_t) ,value :: selfID
    character(kind=c_char) ,dimension(*) ::  export_objs
    integer(c_int) ,value :: obj_size
    integer(c_int) ,dimension(*) :: sizes
    integer(c_int) ::  len_import_objs
    character(kind=c_char) ,dimension(len_import_objs) :: import_objs 
  end function
     
  !/* Original C++ prototype:
  !  virtual int DoPosts( char * export_objs, int obj_size, int *& sizes, int & len_import_objs, char *& import_objs) = 0;
  !*/
  !int Epetra_Distributor_DoPosts_VarLen ( 
  !  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  !  int ** sizes, int * len_import_objs, char ** import_objs );
 
  integer(c_int) function Epetra_Distributor_DoPosts_VarLen( selfID, export_objs, obj_size, sizes, len_import_objs, import_objs ) &
    bind(C,name='Epetra_Distributor_DoPosts_VarLen')
    import :: c_int ,FT_Epetra_Distributor_ID_t ,c_char 
    implicit none
    type(FT_Epetra_Distributor_ID_t) ,value :: selfID
    character(kind=c_char) ,dimension(*) :: export_objs
    integer(c_int) ,value :: obj_size 
    integer(c_int) ,dimension(*) :: sizes
    integer(c_int)  :: len_import_objs
    character(kind=c_char) ,dimension(len_import_objs) :: import_objs
  end function
 
  !/* Original C++ prototype:
  !   virtual int DoReversePosts( char * export_objs, int obj_size, int *& sizes, int & len_import_objs, char *& import_objs) = 0;
  !*/
  !int Epetra_Distributor_DoReversePosts_VarLen ( 
  !  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  !  int ** sizes, int * len_import_objs, char ** import_objs );

  integer(c_int) function Epetra_Distributor_DoReversePosts_VarLen(selfID,export_objs,obj_size,sizes,len_import_objs,import_objs) &
    bind(C,name='Epetra_Distributor_DoReversePosts_VarLen')
    import :: c_int ,FT_Epetra_Distributor_ID_t ,c_char
    implicit none
    type(FT_Epetra_Distributor_ID_t) ,value :: selfID
    character(kind=c_char) ,dimension(*) :: export_objs
    integer(c_int) ,value :: obj_size
    integer(c_int) ,dimension(*) :: sizes
    integer(c_int) :: len_import_objs
    character(kind=c_char) ,dimension(len_import_objs) :: import_objs 
  end function

  ! ____________________ Epetra_Comm interface bodies ________________________________________________________

  !CT_Epetra_Comm_ID_t Epetra_Comm_Cast(
  !  CTrilinos_Object_ID_t id );

  type(FT_Epetra_Comm_ID_t) function Epetra_Comm_Cast( id ) bind(C,name='Epetra_Comm_Cast')
    import :: FT_Epetra_Comm_ID_t, ForTrilinos_Object_ID_t 
    type(ForTrilinos_Object_ID_t) ,value :: id  
  end function

  !/* Original C++ prototype:
  !   virtual Epetra_Comm * Clone() const = 0;
  !*/
  !CT_Epetra_Comm_ID_t Epetra_Comm_Clone ( CT_Epetra_Comm_ID_t selfID );

  type(FT_Epetra_Comm_ID_t) function Epetra_Comm_Clone ( selfID ) bind(C,name='Epetra_Comm_Clone')
    import :: FT_Epetra_Comm_ID_t 
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID
  end function

  !/* Original C++ prototype:
  !  virtual ~Epetra_Comm();
  !*/
  !void Epetra_Comm_Destroy ( CT_Epetra_Comm_ID_t * selfID );

  subroutine Epetra_Comm_Destroy( selfID ) bind(C,name='Epetra_Comm_Destroy') 
    import :: FT_Epetra_Comm_ID_t
    implicit none
    type(FT_Epetra_Comm_ID_t) :: selfID
  end subroutine

  !/* Original C++ prototype:
  !   virtual void Barrier() const = 0;
  !*/
  !void Epetra_Comm_Barrier ( CT_Epetra_Comm_ID_t selfID );

  subroutine Epetra_Comm_Barrier ( selfID ) bind(C,name='Epetra_Comm_Barrier')
    import :: FT_Epetra_Comm_ID_t 
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID 
  end subroutine

  !/* Original C++ prototype:
  !   virtual int Broadcast(double * MyVals, int Count, int Root) const = 0;
  !*/
  !int Epetra_Comm_Broadcast_Double ( 
  !  CT_Epetra_Comm_ID_t selfID, double * MyVals, int Count, int Root );

  integer(c_int) function Epetra_Comm_Broadcast_Double( selfID, MyVals, Count, Root ) bind(C,name='Epetra_Comm_Broadcast_Double')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_double 
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID
    integer(c_int) ,value :: Count
    real(c_double) ,dimension(Count) :: MyVals
    integer(c_int) ,value :: Root 
  end function

  !/* Original C++ prototype:
  !   virtual int Broadcast(int * MyVals, int Count, int Root) const = 0;
  !*/
  !int Epetra_Comm_Broadcast_Int ( 
  !  CT_Epetra_Comm_ID_t selfID, int * MyVals, int Count, int Root );

  integer(c_int) function Epetra_Comm_Broadcast_Int ( selfID, MyVals, Count, Root ) bind(C,name='Epetra_Comm_Broadcast_Int ')
    import :: FT_Epetra_Comm_ID_t , c_int 
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID
    integer(c_int) ,value :: Count
    integer(c_int) ,dimension(Count) :: MyVals
    integer(c_int) ,value ::  Root 
  end function

  !/* Original C++ prototype:
  !   virtual int Broadcast(long * MyVals, int Count, int Root) const = 0;
  !  */
  !int Epetra_Comm_Broadcast_Long ( 
  !  CT_Epetra_Comm_ID_t selfID, long * MyVals, int Count, int Root );

  integer(c_int) function Epetra_Comm_Broadcast_Long( selfID, MyVals, Count, Root ) bind(C,name='Epetra_Comm_Broadcast_Long')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_long
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID
    integer(c_int) ,value :: Count
    integer(c_long) ,dimension(Count) :: MyVals
    integer(c_int) ,value :: Root 
  end function

  !/* Original C++ prototype:
  !   virtual int Broadcast(char * MyVals, int Count, int Root) const = 0;
  !*/
  !int Epetra_Comm_Broadcast_Char ( 
  !  CT_Epetra_Comm_ID_t selfID, char * MyVals, int Count, int Root );

  integer(c_int) function Epetra_Comm_Broadcast_Char ( selfID, MyVals, Count, Root ) bind(C,name='Epetra_Comm_Broadcast_Char')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_char 
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID
    integer(c_int) ,value :: Count
    character(kind=c_char) ,dimension(Count) :: MyVals
    integer(c_int) ,value :: Root 
  end function

  !/* Original C++ prototype:
  !   virtual int GatherAll(double * MyVals, double * AllVals, int Count) const = 0;
  !*/
  !int Epetra_Comm_GatherAll_Double ( 
  !  CT_Epetra_Comm_ID_t selfID, double * MyVals, double * AllVals, 
  !  int Count );

  integer(c_int) function Epetra_Comm_GatherAll_Double( selfID, MyVals, AllVals, Count ) bind(C,name='Epetra_Comm_GatherAll_Double')
    import :: c_int ,FT_Epetra_Comm_ID_t , c_double
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID
    real(c_double) ,dimension(*) :: AllVals 
    integer(c_int) ,value :: Count 
    real(c_double) ,dimension(Count) :: MyVals
  end function

  !/* Original C++ prototype:
  !   virtual int GatherAll(int * MyVals, int * AllVals, int Count) const = 0;
  !*/
  !int Epetra_Comm_GatherAll_Int ( 
  !  CT_Epetra_Comm_ID_t selfID, int * MyVals, int * AllVals, int Count );

  integer(c_int) function Epetra_Comm_GatherAll_Int (selfID, MyVals, AllVals, Count ) bind(C,name='Epetra_Comm_GatherAll_Int')
    import :: c_int ,FT_Epetra_Comm_ID_t 
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID
    integer(c_int) ,dimension(*) ::  AllVals
    integer(c_int) ,value :: Count 
    integer(c_int) ,dimension(Count) :: MyVals
  end function

  !/* Original C++ prototype:
  !   virtual int GatherAll(long * MyVals, long * AllVals, int Count) const = 0;
  !*/
  !int Epetra_Comm_GatherAll_Long ( 
  !  CT_Epetra_Comm_ID_t selfID, long * MyVals, long * AllVals, 
  !  int Count );

  integer(c_int) function Epetra_Comm_GatherAll_Long ( selfID, MyVals, AllVals, Count ) bind(C,name='Epetra_Comm_GatherAll_Long')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_long 
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID
    integer(c_long) ,dimension(*) :: AllVals
    integer(c_int) ,value :: Count
    integer(c_long) ,dimension(Count) :: MyVals
  end function

  !/* Original C++ prototype:
  !   virtual int SumAll(double * PartialSums, double * GlobalSums, int Count) const = 0;
  !*/
  !int Epetra_Comm_SumAll_Double ( 
  !  CT_Epetra_Comm_ID_t selfID, double * PartialSums, 
  !  double * GlobalSums, int Count );

  integer(c_int) function Epetra_Comm_SumAll_Double(selfID,PartialSums,GlobalSums,Count) bind(C,name='Epetra_Comm_SumAll_Double')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_double
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID
    integer(c_int) ,value :: Count
    real(c_double) ,dimension(Count) :: PartialSums 
    real(c_double) ,dimension(Count) :: GlobalSums
  end function

  !/* Original C++ prototype:
  !   virtual int SumAll(int * PartialSums, int * GlobalSums, int Count) const = 0;
  !*/
  !int Epetra_Comm_SumAll_Int ( 
  !  CT_Epetra_Comm_ID_t selfID, int * PartialSums, int * GlobalSums, 
  !  int Count );

  integer(c_int) function Epetra_Comm_SumAll_Int ( selfID, PartialSums, GlobalSums, Count ) bind(C,name='Epetra_Comm_SumAll_Int')
    import :: c_int ,FT_Epetra_Comm_ID_t 
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID
    integer(c_int) ,value :: Count
    integer(c_int) ,dimension(Count) :: PartialSums
    integer(c_int) ,dimension(Count) :: GlobalSums
  end function

  !/* Original C++ prototype:
  !   virtual int SumAll(long * PartialSums, long * GlobalSums, int Count) const = 0;
  !*/
  !int Epetra_Comm_SumAll_Long ( 
  !  CT_Epetra_Comm_ID_t selfID, long * PartialSums, long * GlobalSums, 
  !  int Count );

  integer(c_int) function Epetra_Comm_SumAll_Long ( selfID, PartialSums, GlobalSums, Count ) bind(C,name='Epetra_Comm_SumAll_Long')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_long 
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID
    integer(c_int) ,value :: Count 
    integer(c_long) ,dimension(Count) :: PartialSums
    integer(c_long) ,dimension(Count) :: GlobalSums 
  end function

  !/* Original C++ prototype:
  !  virtual int MaxAll(double * PartialMaxs, double * GlobalMaxs, int Count) const = 0;
  !*/
  !int Epetra_Comm_MaxAll_Double ( 
  !  CT_Epetra_Comm_ID_t selfID, double * PartialMaxs, 
  !  double * GlobalMaxs, int Count );
  
  integer(c_int) function Epetra_Comm_MaxAll_Double ( selfID,PartialMaxs,GlobalMaxs,Count ) bind(C,name='Epetra_Comm_MaxAll_Double')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_double 
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID
    integer(c_int) ,value :: Count
    real(c_double) ,dimension(Count) :: PartialMaxs
    real(c_double) ,dimension(Count) :: GlobalMaxs
  end function

  !/* Original C++ prototype:
  !   virtual int MaxAll(int * PartialMaxs, int * GlobalMaxs, int Count) const = 0;
  !*/
  !int Epetra_Comm_MaxAll_Int ( 
  !  CT_Epetra_Comm_ID_t selfID, int * PartialMaxs, int * GlobalMaxs, 
  !  int Count );

  integer(c_int) function Epetra_Comm_MaxAll_Int ( selfID, PartialMaxs, GlobalMaxs, Count ) bind(C,name='Epetra_Comm_MaxAll_Int')
    import :: c_int ,FT_Epetra_Comm_ID_t
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID
    integer(c_int) ,value :: Count
    integer(c_int) ,dimension(Count) :: PartialMaxs
    integer(c_int) ,dimension(Count) :: GlobalMaxs
  end function
  
  !/* Original C++ prototype:
  !   virtual int MaxAll(long * PartialMaxs, long * GlobalMaxs, int Count) const = 0;
  !*/
  !int Epetra_Comm_MaxAll_Long ( 
  !  CT_Epetra_Comm_ID_t selfID, long * PartialMaxs, long * GlobalMaxs, 
  !  int Count );

  integer(c_int) function Epetra_Comm_MaxAll_Long ( selfID, PartialMaxs, GlobalMaxs, Count ) bind(C,name='Epetra_Comm_MaxAll_Long')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_long
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID
    integer(c_int) ,value :: Count
    integer(c_long) ,dimension(Count) :: PartialMaxs
    integer(c_long) ,dimension(Count) :: GlobalMaxs
  end function
  
  !/* Original C++ prototype:
  !   virtual int MinAll(double * PartialMins, double * GlobalMins, int Count) const = 0;
  !*/
  !int Epetra_Comm_MinAll_Double ( 
  !  CT_Epetra_Comm_ID_t selfID, double * PartialMins, 
  !  double * GlobalMins, int Count );
    
  integer(c_int) function Epetra_Comm_MinAll_Double(selfID, PartialMins, GlobalMins, Count) bind(C,name='Epetra_Comm_MinAll_Double')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_double 
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID
    integer(c_int) ,value :: Count 
    real(c_double) ,dimension(Count) :: PartialMins
    real(c_double) ,dimension(Count) :: GlobalMins
  end function

  !/* Original C++ prototype:
  !   virtual int MinAll(int * PartialMins, int * GlobalMins, int Count) const = 0;
  !*/
  !int Epetra_Comm_MinAll_Int ( 
  !  CT_Epetra_Comm_ID_t selfID, int * PartialMins, int * GlobalMins, 
  !  int Count );
  
  integer(c_int) function Epetra_Comm_MinAll_Int ( selfID, PartialMins, GlobalMins, Count ) bind(C,name='Epetra_Comm_MinAll_Int')
    import :: c_int ,FT_Epetra_Comm_ID_t 
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID
    integer(c_int) ,value :: Count
    integer(c_int) ,dimension(Count) :: PartialMins
    integer(c_int) ,dimension(Count) :: GlobalMins
  end function

  !/* Original C++ prototype:
  !   virtual int MinAll(long * PartialMins, long * GlobalMins, int Count) const = 0;
  !*/
  !int Epetra_Comm_MinAll_Long ( 
  !  CT_Epetra_Comm_ID_t selfID, long * PartialMins, long * GlobalMins, 
  !  int Count );

  integer(c_int) function Epetra_Comm_MinAll_Long ( selfID, PartialMins, GlobalMins, Count ) bind(C,name='Epetra_Comm_MinAll_Long')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_long 
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID
    integer(c_int) ,value :: Count
    integer(c_long) ,dimension(Count) :: PartialMins
    integer(c_long) ,dimension(Count) :: GlobalMins
  end function
  
  !/* Original C++ prototype:
  !   virtual int ScanSum(double * MyVals, double * ScanSums, int Count) const = 0;
  !*/
  !int Epetra_Comm_ScanSum_Double ( 
  !  CT_Epetra_Comm_ID_t selfID, double * MyVals, double * ScanSums, 
  !  int Count );

  integer(c_int) function Epetra_Comm_ScanSum_Double ( selfID, MyVals, ScanSums, Count ) bind(C,name='Epetra_Comm_ScanSum_Double')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_double 
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID
    integer(c_int) ,value :: Count
    real(c_double) ,dimension(Count) :: MyVals
    real(c_double) ,dimension(Count) :: ScanSums
  end function
  
  !/* Original C++ prototype:
  !   virtual int ScanSum(int * MyVals, int * ScanSums, int Count) const = 0;
  !*/
  !int Epetra_Comm_ScanSum_Int ( 
  !CT_Epetra_Comm_ID_t selfID, int * MyVals, int * ScanSums, 
  !  int Count );

  integer(c_int) function Epetra_Comm_ScanSum_Int( selfID, MyVals, ScanSums, Count ) bind(C,name='Epetra_Comm_ScanSum_Int')
    import :: c_int ,FT_Epetra_Comm_ID_t 
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID
    integer(c_int) ,value :: Count
    integer(c_int) ,dimension(Count) :: MyVals
    integer(c_int) ,dimension(Count) :: ScanSums
  end function
  
  !/* Original C++ prototype:
  !   virtual int ScanSum(long * MyVals, long * ScanSums, int Count) const = 0;
  !*/
  !int Epetra_Comm_ScanSum_Long ( 
  !  CT_Epetra_Comm_ID_t selfID, long * MyVals, long * ScanSums, 
  !  int Count );
  
  integer(c_int) function Epetra_Comm_ScanSum_Long( selfID, MyVals, ScanSums, Count ) bind(C,name='Epetra_Comm_ScanSum_Long')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_long
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID
    integer(c_int) ,value :: Count
    integer(c_long) ,dimension(Count) :: MyVals
    integer(c_long) ,dimension(Count) :: ScanSums
  end function
  
  !/* Original C++ prototype:
  !  virtual int MyPID() const = 0;
  !*/
  !int Epetra_Comm_MyPID ( CT_Epetra_Comm_ID_t selfID );

  integer(c_int) function Epetra_Comm_MyPID ( selfID ) bind(C,name='Epetra_Comm_MyPID')
    import :: c_int ,FT_Epetra_Comm_ID_t
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID 
  end function
  
  !/* Original C++ prototype:
  !   virtual int NumProc() const = 0;
  !*/
  !int Epetra_Comm_NumProc ( FT_Epetra_Comm_ID_t selfID );

  integer(c_int) function Epetra_Comm_NumProc ( selfID ) bind(C,name='Epetra_Comm_NumProc')
    import :: c_int ,FT_Epetra_Comm_ID_t 
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID 
  end function
  
  !/* Original C++ prototype:
  !   virtual Epetra_Distributor * CreateDistributor() const = 0;
  !*/
  !CT_Epetra_Distributor_ID_t Epetra_Comm_CreateDistributor ( 
  !  CT_Epetra_Comm_ID_t selfID );
  
  type(FT_Epetra_Distributor_ID_t) function Epetra_Comm_CreateDistributor ( selfID ) bind(C,name='Epetra_Comm_CreateDistributor')
    import :: FT_Epetra_Distributor_ID_t ,FT_Epetra_Comm_ID_t 
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID 
  end function

  !/* Original C++ prototype:
  !   virtual Epetra_Directory * CreateDirectory(const Epetra_BlockMap & Map) const = 0;
  !*/
  !CT_Epetra_Directory_ID_t Epetra_Comm_CreateDirectory ( 
  !  CT_Epetra_Comm_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID );
   
  type(FT_Epetra_Directory_ID_t) function Epetra_Comm_CreateDirectory ( selfID, MapID ) bind(C,name='Epetra_Comm_CreateDirectory') 
    import :: FT_Epetra_Directory_ID_t ,FT_Epetra_Comm_ID_t ,FT_Epetra_BlockMap_ID_t 
    implicit none
    type(FT_Epetra_Comm_ID_t) ,value :: selfID 
    type(FT_Epetra_BlockMap_ID_t) ,value :: MapID 
  end function

  ! ___________________ Epetra_BlockMap interface bodies ____________________________________________


  ! ___________________ Epetra_Map interface bodies ____________________________________________

  !CT_Epetra_Map_ID_t Epetra_Map_Cast(
  !  CTrilinos_Object_ID_t id );

  type(FT_Epetra_Map_ID_t) function Epetra_Map_Cast(id ) bind(C,name='Epetra_Map_Cast')
    import :: FT_Epetra_Map_ID_t ,ForTrilinos_Object_ID_t 
    implicit none
    type(ForTrilinos_Object_ID_t) ,value :: id
  end function 

  !/* Original C++ prototype:
  !   Epetra_Map(int NumGlobalElements, int IndexBase, const Epetra_Comm& Comm);
  !*/
  !CT_Epetra_Map_ID_t Epetra_Map_Create ( 
  !  int NumGlobalElements, int IndexBase, CT_Epetra_Comm_ID_t CommID );

  type(FT_Epetra_Map_ID_t) function Epetra_Map_Create ( NumGlobalElements, IndexBase, CommID ) bind(C,name='Epetra_Map_Create')
    import :: FT_Epetra_Map_ID_t ,c_int ,FT_Epetra_Comm_ID_t 
    implicit none
    integer(c_int) ,value :: NumGlobalElements
    integer(c_int) ,value :: IndexBase
    type(FT_Epetra_Comm_ID_t) ,value :: CommID
  end function 

  !/* Original C++ prototype:
  !   Epetra_Map(int NumGlobalElements, int NumMyElements, int IndexBase, const Epetra_Comm& Comm);
  !*/
  !CT_Epetra_Map_ID_t Epetra_Map_Create_Linear ( 
  !  int NumGlobalElements, int NumMyElements, int IndexBase, 
  !  CT_Epetra_Comm_ID_t CommID );

  type(FT_Epetra_Map_ID_t) function Epetra_Map_Create_Linear ( NumGlobalElements, NumMyElements, IndexBase, CommID ) &
    bind(C,name='Epetra_Map_Create_Linear')
    import :: FT_Epetra_Map_ID_t ,c_int ,FT_Epetra_Comm_ID_t 
    implicit none
    integer(c_int) ,value :: NumGlobalElements
    integer(c_int) ,value :: NumMyElements
    integer(c_int) ,value :: IndexBase
    type(FT_Epetra_Comm_ID_t) ,value :: CommID 
  end function 
  
  !/* Original C++ prototype:
  !   Epetra_Map(int NumGlobalElements, int NumMyElements, const int *MyGlobalElements, int IndexBase, const Epetra_Comm& Comm);
  !*/
  !CT_Epetra_Map_ID_t Epetra_Map_Create_Arbitrary ( 
  !  int NumGlobalElements, int NumMyElements, 
  !  const int * MyGlobalElements, int IndexBase, 
  !  CT_Epetra_Comm_ID_t CommID );
  
  type(FT_Epetra_Map_ID_t) function Epetra_Map_Create_Arbitrary ( &
    NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID ) bind(C,name='Epetra_Map_Create_Arbitrary') 
    import :: FT_Epetra_Map_ID_t ,c_int ,FT_Epetra_Comm_ID_t 
    implicit none
    integer(c_int) ,value :: NumGlobalElements
    integer(c_int) ,value :: NumMyElements
    integer(c_int) ,intent(in) ,dimension(NumMyElements) :: MyGlobalElements
    integer(c_int) ,value :: IndexBase
    type(FT_Epetra_Comm_ID_t) ,value :: CommID 
  end function 

  !/* Original C++ prototype:
  !   Epetra_Map(const Epetra_Map& map);
  !*/
  !CT_Epetra_Map_ID_t Epetra_Map_Duplicate ( CT_Epetra_Map_ID_t mapID );

  type(FT_Epetra_Map_ID_t) function Epetra_Map_Duplicate ( mapID ) bind(C,name='Epetra_Map_Duplicate')
    import :: FT_Epetra_Map_ID_t
    implicit none
    type(FT_Epetra_Map_ID_t) ,intent(in) ,value :: mapID 
  end function 

  !/* Original C++ prototype:
  !   virtual ~Epetra_Map(void);
  !*/
  !void Epetra_Map_Destroy ( CT_Epetra_Map_ID_t * selfID );

  subroutine Epetra_Map_Destroy ( selfID ) bind(C,name='Epetra_Map_Destroy')
    import :: FT_Epetra_Map_ID_t
    implicit none
    type(FT_Epetra_Map_ID_t) :: selfID
  end subroutine
  
  !/* Original C++ prototype:
  !   Epetra_Map & operator=(const Epetra_Map & map);
  !*/
  !void Epetra_Map_Assign ( 
  !  CT_Epetra_Map_ID_t selfID, CT_Epetra_Map_ID_t mapID );

  subroutine Epetra_Map_Assign ( selfID, mapID ) bind(C,name='Epetra_Map_Assign')
    import :: FT_Epetra_Map_ID_t 
    implicit none
    type(FT_Epetra_Map_ID_t) ,value :: selfID
    type(FT_Epetra_Map_ID_t) ,value ,intent(in) :: mapID 
  end subroutine

  end interface
end module forepetra

! !CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Cast(
! !  CTrilinos_Object_ID_t id );

! type(FT_Epetra_BlockMap_ID_t) function Epetra_BlockMap_Cast(id ) bind(C,name='Epetra_BlockMap_Cast')
!   import :: FT_Epetra_BlockMap_ID_t ,ForTrilinos_Object_ID_t 
!   implicit none
!   type(ForTrilinos_Object_ID_t) ,value :: id 
! end function

! !/* Original C++ prototype:
! !   Epetra_BlockMap(int NumGlobalElements, int ElementSize, int IndexBase, const Epetra_Comm& Comm);
! !*/
! !CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create ( 
! !  int NumGlobalElements, int ElementSize, int IndexBase, 
! !  CT_Epetra_Comm_ID_t CommID );

! type(FT_Epetra_BlockMap_ID_t) Epetra_BlockMap_Create ( NumGlobalElements, ElementSize, IndexBase, CommID ) &
!   bind(C,name='Epetra_BlockMap_Create') 
!   import :: FT_Epetra_BlockMap_ID_t ,c_int ,FT_Epetra_Comm_ID_t 
!   implicit none
!   integer(c_int) ,value :: NumGlobalElements
!   integer(c_int) ,value :: ElementSize
!   integer(c_int) ,value :: IndexBase, 
!   type(FT_Epetra_Comm_ID_t) ,value :: CommID 
! end function


! !/* Original C++ prototype:
! !   Epetra_BlockMap(int NumGlobalElements, int NumMyElements, int ElementSize, int IndexBase, const Epetra_Comm& Comm);
! !*/
! !CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Linear ( 
! !  int NumGlobalElements, int NumMyElements, int ElementSize, 
! !  int IndexBase, CT_Epetra_Comm_ID_t CommID );

! type(FT_Epetra_BlockMap_ID_t) function Epetra_BlockMap_Create_Linear( &
!   NumGlobalElements, NumMyElements, ElementSize, IndexBase, CommID) bind(C,name='Epetra_BlockMap_Create_Linear')
!   import :: FT_Epetra_BlockMap_ID_t ,c_int ,FT_Epetra_Comm_ID_t 
!   implicit none
!   integer(c_int) ,value :: NumGlobalElements
!   integer(c_int) ,value :: NumMyElements
!   integer(c_int) ,value :: ElementSize
!   integer(c_int) ,value :: IndexBase
!   type(FT_Epetra_Comm_ID_t) ,value ,intent(in) :: CommID
! end function

! !/* Original C++ prototype:
! !  Epetra_BlockMap(int NumGlobalElements, int NumMyElements, const int *MyGlobalElements, int ElementSize, int IndexBase, const Epetra_Comm& Comm);
! !*/
! !CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Arbitrary ( 
! !  int NumGlobalElements, int NumMyElements, 
! !  const int * MyGlobalElements, int ElementSize, int IndexBase, 
! !  CT_Epetra_Comm_ID_t CommID );

! type(CT_Epetra_BlockMap_ID_t) function Epetra_BlockMap_Create_Arbitrary ( &
!   NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize, IndexBase, CommID ) &
!   bind(C,name='Epetra_BlockMap_Create_Arbitrary') 
!   import :: CT_Epetra_BlockMap_ID_t ,c_int ,CT_Epetra_Comm_ID_t 
!   implicit none
!   integer(c_int) ,value :: NumGlobalElements
!   integer(c_int) ,value :: NumMyElements, 
!   integer(c_int) ,intent(in) ,dimension(NumMyElements) :: MyGlobalElements
!   integer(c_int) ,value :: ElementSize
!   integer(c_int) ,vaue :: IndexBase
!   type(CT_Epetra_Comm_ID_t) ,intent(in) ,value :: CommID
! end function
