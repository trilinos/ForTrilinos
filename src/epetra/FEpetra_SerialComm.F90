module FEpetra_SerialComm
  use ForTrilinos_enums ,only : FT_Epetra_Comm_ID,FT_Epetra_SerialComm_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use FEpetra_Comm      ,only : Epetra_Comm
  use iso_c_binding     ,only : c_int,c_long,c_double,c_char
  use forepetra
#include "ForTrilinos_config.h"
  implicit none
  private                     ! Hide everything by default
  public :: Epetra_SerialComm ! Expose type/constructors/methods

  type ,extends(Epetra_Comm)                 :: Epetra_SerialComm !"shell"
    private
    type(FT_Epetra_SerialComm_ID_t) ,pointer :: SerialComm_id => null()
  contains
     !Constructor
     procedure         :: clone 
     !Developers only
     procedure         :: get_EpetraSerialComm_ID 
     procedure ,nopass :: alias_EpetraSerialComm_ID
     procedure         :: generalize 
     procedure         :: SerialComm_assign => assign_to_Epetra_SerialComm
     procedure         :: Comm_assign => assign_to_Epetra_Comm
#ifdef HAVE_MPI
     procedure         :: MpiComm_assign
#endif
     !Barrier Methods
     procedure         :: barrier
     !Broadcast Methods
     procedure         :: broadcast_double
     procedure         :: broadcast_int
     procedure         :: broadcast_long
     procedure         :: broadcast_char
     !Gather Methods
     procedure         :: gather_double
     !Sum Methods
     !generic :: SumAll=>
     !Max/Min Methods
     !generic :: MaxAll=>
     !generic :: MinAll=>
     !Parallel Prefix Methods
     !Attribute Accessor Methods
     procedure         :: MyPID
     procedure         :: NumProc
     !Gather/catter and Directory Constructors
     !I/O methods
     !Memory Management
     procedure         :: force_finalization 
     final :: finalize
  end type

   interface Epetra_SerialComm ! constructors
     module procedure from_scratch,duplicate,from_struct
   end interface

contains

  type(FT_Epetra_SerialComm_ID_t) function from_struct(id)
     type(FT_Epetra_SerialComm_ID_t) ,intent(in) :: id
     from_struct = id
  end function

  ! Original C++ prototype:
  ! Epetra_SerialComm();
  ! CTrilinos prototype:
  ! CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Create (  );
  
  type(FT_Epetra_SerialComm_ID_t) function from_scratch()
    from_scratch = Epetra_SerialComm_Create()
  end function

  ! Original C++ prototype:
  ! Epetra_SerialComm(const Epetra_SerialComm& Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Duplicate ( CT_Epetra_SerialComm_ID_t CommID );

  type(FT_Epetra_SerialComm_ID_t) function duplicate(this)
    type(Epetra_SerialComm) ,intent(in) :: this 
    duplicate = Epetra_SerialComm_Duplicate(this%SerialComm_id)
  end function

  function clone(this)
    class(Epetra_SerialComm) ,intent(in)  :: this
    class(Epetra_Comm)       ,allocatable :: clone
    class(Epetra_Comm)       ,allocatable :: clone_temp
    type(FT_Epetra_SerialComm_ID_t) :: test
    type(FT_Epetra_Comm_ID_t) :: test1
    allocate(Epetra_SerialComm :: clone) 
    allocate(Epetra_SerialComm :: clone_temp) 
    clone_temp = Epetra_SerialComm_Clone(this%SerialComm_id)
   ! test = clone_temp%SerialComm_id
   ! print *,'clone_temp%serialcom',test%table,test%index
    test1 = clone_temp%get_EpetraComm_ID()
    print *,'clone_temp%comm',test1%table,test1%index
    clone=Epetra_SerialComm(alias_EpetraSerialComm_ID(clone_temp%generalize_EpetraComm()))
    !test = clone%SerialComm_id
    !print *,'clone%serialcomm',test%table,test%index
    test1 = clone%get_EpetraComm_ID()
    print *,'clone%comm',test1%table,test1%index
    call clone_temp%force_finalization_EpetraComm()
  end function

  type(FT_Epetra_SerialComm_ID_t) function get_EpetraSerialComm_ID(this)
   class(Epetra_SerialComm) ,intent(in) :: this 
   if (associated(this%SerialComm_id)) then
    get_EpetraSerialComm_ID=this%SerialComm_id
   else
    stop 'get_EpetraSerialComm_ID: SerialComm_id is unassociated'
   end if
  end function
  
  type(FT_Epetra_SerialComm_ID_t) function alias_EpetraSerialComm_ID(generic_id)
    use iso_c_binding        ,only: c_loc
    use ForTrilinos_enums    ,only: FT_Epetra_SerialComm_ID,ForTrilinos_Universal_ID_t
    use ForTrilinos_table_man,only: CT_Alias
    type(Fortrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(Fortrilinos_Universal_ID_t) ,pointer    :: alias_id
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_SerialComm_ID))
    alias_EpetraSerialComm_ID=degeneralize_EpetraSerialComm(c_loc(alias_id))
    deallocate(alias_id)
  end function


  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding ,only : c_loc
   class(Epetra_SerialComm) ,intent(in) ,target :: this
   generalize = generalize_all( c_loc(this%SerialComm_id) )
   ! ____ Use for ForTrilinos function implementation ______
   
   ! ____ Use for CTrilinos function implementation ______
   ! class(Epetra_SerialComm) ,intent(in) ,target :: this
   ! generalize = Epetra_SerialComm_Generalize ( this%SerialComm_id ) 
   ! ____ Use for CTrilinos function implementation ______
  end function
 
 type(FT_Epetra_SerialComm_ID_t) function degeneralize_EpetraSerialComm(generic_id) bind(C)
   ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_SerialComm_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                     ,value   :: generic_id
    type(FT_Epetra_SerialComm_ID_t) ,pointer :: local_ptr
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraSerialComm = local_ptr
   ! ____ Use for ForTrilinos function implementation ______
   
   ! ____ Use for CTrilinos function implementation ______
   ! type(Fortrilinos_Universal_ID_t) ,intent(in) :: generic_id
   ! degeneralize_EpetraSerialComm = Epetra_SerialComm_Degeneralize( generic_id )
   ! ____ Use for CTrilinos function implementation ______
  end function
 
  subroutine assign_to_Epetra_SerialComm(lhs,rhs)
    class(Epetra_SerialComm)        ,intent(inout) :: lhs
    type(FT_Epetra_SerialComm_ID_t) ,intent(in)    :: rhs
    allocate( lhs%SerialComm_id, source=rhs)
    call lhs%set_EpetraComm_ID(lhs%alias_EpetraComm_ID(lhs%generalize()))
  end subroutine

#ifdef HAVE_MPI 
  subroutine MpiComm_assign(lhs,rhs)
    use ForTrilinos_enums, only : FT_Epetra_MpiComm_ID_t
    class(Epetra_SerialComm),      intent(inout) :: lhs
    type(FT_Epetra_MpiComm_ID_t), intent(in) :: rhs 
    print *, 'MpiComm_assign in FEpetra_SerialComm no-op'
  end subroutine
#endif
  subroutine assign_to_Epetra_Comm(lhs,rhs)
    class(Epetra_SerialComm) ,intent(inout) :: lhs
    class(Epetra_Comm)       ,intent(in)    :: rhs
    select type(rhs)
      class is (Epetra_SerialComm)
        allocate(lhs%SerialComm_id,source=alias_EpetraSerialComm_ID(rhs%generalize()))
        call lhs%set_EpetraComm_ID(lhs%alias_EpetraComm_ID(lhs%generalize()))
        !allocate(Epetra_SerialComm :: lhs)
        !lhs=Epetra_SerialComm(alias_EpetraSerialComm_ID(rhs%generalize()))
     class default
        stop 'assign_to_Epetra_Comm: unsupported class'
     end select
  end subroutine

  subroutine barrier(this)
   class(Epetra_SerialComm) ,intent(in) :: this
   call Epetra_SerialComm_Barrier(this%SerialComm_id)
  end subroutine
 
  subroutine broadcast_double(this,MyVals,count,root)
   class(Epetra_SerialComm)     ,intent(in)    :: this
   real(c_double), dimension(:) ,intent(inout) :: MyVals
   integer(c_int)               ,intent(in)    :: count
   integer(c_int)               ,intent(in)    :: root
   integer(c_int)                              :: error 
   error = Epetra_SerialComm_Broadcast_Double(this%SerialComm_id,MyVals,count,root)
  end subroutine

  subroutine broadcast_int(this,MyVals,count,root)
   class(Epetra_SerialComm)     ,intent(in)    :: this
   integer(c_int), dimension(:) ,intent(inout) :: MyVals
   integer(c_int)               ,intent(in)    :: count
   integer(c_int)               ,intent(in)    :: root
   integer(c_int)                              :: error 
   error = Epetra_SerialComm_Broadcast_Int(this%SerialComm_id,MyVals,count,root)
  end subroutine

  subroutine broadcast_long(this,MyVals,count,root)
   class(Epetra_SerialComm)     ,intent(in)    :: this
   integer(c_long),dimension(:) ,intent(inout) :: MyVals
   integer(c_int)               ,intent(in)    :: count
   integer(c_int)               ,intent(in)    :: root
   integer(c_int)                              :: error 
   error = Epetra_SerialComm_Broadcast_Long(this%SerialComm_id,MyVals,count,root)
  end subroutine
 
  subroutine broadcast_char(this,MyVals,count,root)
   class(Epetra_SerialComm)           ,intent(in)    :: this
   character(kind=c_char),dimension(:),intent(inout) :: MyVals
   integer(c_int)                     ,intent(in)    :: count
   integer(c_int)                     ,intent(in)    :: root
   integer(c_int)                                    :: error 
   error = Epetra_SerialComm_Broadcast_Char(this%SerialComm_id,MyVals,count,root)
  end subroutine

  subroutine gather_double(this,MyVals,AllVals,count)
   class(Epetra_SerialComm)     ,intent(in)    :: this
   !real(c_double), dimension(:) ,intent(inout) :: MyVals
   !real(c_double), dimension(:) ,intent(inout) :: AllVals
   real(c_double), dimension(:)  :: MyVals
   real(c_double), dimension(:)  :: AllVals
   integer(c_int)               ,intent(in)    :: count
   integer(c_int)                              :: error 
   error = Epetra_SerialComm_GatherAll_Double(this%SerialComm_id,MyVals,AllVals,count)
  end subroutine

  integer(c_int) function MyPID(this)
   class(Epetra_SerialComm)     , intent(in) :: this
   MyPID=Epetra_SerialComm_MyPID(this%SerialComm_id)
  end function

  integer(c_int) function NumProc(this)
   class(Epetra_SerialComm)     , intent(in) :: this
   NumProc=Epetra_SerialComm_NumProc(this%SerialComm_id)
  end function

  subroutine finalize(this)
    type(Epetra_SerialComm)     :: this
    call this%force_finalization_EpetraComm()
    call Epetra_SerialComm_Destroy( this%SerialComm_id ) 
    deallocate(this%SerialComm_id)
  end subroutine

  subroutine force_finalization(this)
    class(Epetra_SerialComm) ,intent(inout) :: this
    if (associated(this%SerialComm_id)) then
      call finalize(this) 
    else
      print *,' finalization for Epetra_SerialComm received object with unassociated SerialComm_id'
    end if
  end subroutine
end module 

