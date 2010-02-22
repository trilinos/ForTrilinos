module FEpetra_SerialComm
  use ForTrilinos_enums ,only : FT_Epetra_Comm_ID,FT_Epetra_SerialComm_ID_t,ForTrilinos_Universal_ID_t
  use FEpetra_Comm      ,only : epetra_comm
  use iso_c_binding     ,only : c_int,c_long,c_double,c_char
  use forepetra
  implicit none
  private                     ! Hide everything by default
  public :: epetra_serialcomm ! Expose type/constructors/methods

  type ,extends(epetra_comm)                 :: epetra_serialcomm !"shell"
    private
    type(FT_Epetra_SerialComm_ID_t) ,pointer :: SerialComm_id => null()
  contains
     !Constructor
     procedure         :: clone 
     !Developers only
     procedure         :: get_EpetraSerialComm_ID 
     procedure ,nopass :: alias_EpetraSerialComm_ID
     procedure         :: generalize 
     procedure         :: SerialComm_assign => assign_to_epetra_SerialComm
     procedure         :: Comm_assign => assign_to_epetra_Comm
     !Barrier Methods
     procedure         :: barrier
     !Broadcast Methods
     procedure,private :: broadcast_double
     procedure,private :: broadcast_int
     !procedure,private :: broadcast_long
     procedure,private :: broadcast_char
     !generic :: broadcast=>broadcast_double,broadcast_int,broadcast_long,broadcast_char
     generic :: broadcast=>broadcast_double,broadcast_int,broadcast_char
     !Memory Management
     procedure         :: force_finalization 
     final :: finalize
  end type

   interface epetra_serialcomm ! constructors
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
    type(epetra_serialcomm) ,intent(in) :: this 
    duplicate = Epetra_SerialComm_Duplicate(this%SerialComm_id)
  end function

  function clone(this)
    class(epetra_serialcomm) ,intent(in)  :: this
    class(epetra_comm)       ,allocatable :: clone
    allocate(epetra_serialcomm :: clone) 
    clone = Epetra_SerialComm_Clone(this%SerialComm_id)
    !allocate(clone%SerialComm_id,source=clone%alias_EpetraSerialComm(clone%generalize_Comm()))
  end function

  type(FT_Epetra_SerialComm_ID_t) function get_EpetraSerialComm_ID(this)
   class(epetra_serialcomm) ,intent(in) :: this 
   if (associated(this%SerialComm_id)) then
    get_EpetraSerialComm_ID=this%SerialComm_id
   else
    stop 'get_EpetraSerialComm_ID: SerialComm_id is unassociated'
   end if
  end function
  
  type(FT_Epetra_SerialComm_ID_t) function alias_EpetraSerialComm_ID(generic_id)
    use iso_c_binding     ,only: c_loc
    use ForTrilinos_enums ,only: FT_Epetra_SerialComm_ID,ForTrilinos_Universal_ID_t
    use ForTrilinos_utils ,only: CT_Alias
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
   class(epetra_serialcomm) ,intent(in) ,target :: this
   generalize = generalize_all( c_loc(this%SerialComm_id) )
   ! ____ Use for ForTrilinos function implementation ______
   
   ! ____ Use for CTrilinos function implementation ______
   ! class(epetra_serialcomm) ,intent(in) ,target :: this
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
 
  subroutine assign_to_epetra_SerialComm(lhs,rhs)
    class(epetra_serialcomm)        ,intent(inout) :: lhs
    type(FT_Epetra_SerialComm_ID_t) ,intent(in)    :: rhs
    type(FT_Epetra_Comm_ID_t)                      :: test  ! test line
    allocate( lhs%SerialComm_id, source=rhs)
    PRINT *,'assignment after serial part'
    PRINT *,lhs%SerialComm_id%table,lhs%SerialComm_id%index,lhs%SerialComm_id%is_const
    call lhs%set_EpetraComm_ID(lhs%alias_EpetraComm_ID(lhs%generalize()))
    test=lhs%get_EpetraComm_ID()  ! test line
    PRINT *,test%table,test%index,test%is_const
  end subroutine

  subroutine assign_to_epetra_Comm(lhs,rhs)
    class(epetra_serialcomm) ,intent(inout) :: lhs
    class(epetra_comm)       ,intent(in)    :: rhs
    select type(rhs)
      class is (epetra_serialcomm)
        allocate(lhs%SerialComm_id,source=rhs%SerialComm_id)
        call lhs%set_EpetraComm_ID(lhs%alias_EpetraComm_ID(lhs%generalize()))
     class default
        stop 'assign_to_epetra_Comm: unsupported class'
     end select
  end subroutine

  subroutine barrier(this)
   class(epetra_serialcomm) ,intent(in) :: this
   call Epetra_SerialComm_Barrier(this%SerialComm_id)
  end subroutine
 
  subroutine broadcast_double(this,MyVals,count,root)
   class(epetra_serialcomm)     ,intent(in)    :: this
   real(c_double), dimension(:) ,intent(inout) :: MyVals
   integer(c_int)               ,intent(in)    :: count
   integer(c_int)               ,intent(in)    :: root
   integer(c_int)                              :: error 
   error = Epetra_SerialComm_Broadcast_Double(this%SerialComm_id,MyVals,count,root)
  end subroutine

  subroutine broadcast_int(this,MyVals,count,root)
   class(epetra_serialcomm)     ,intent(in)    :: this
   integer(c_int), dimension(:) ,intent(inout) :: MyVals
   integer(c_int)               ,intent(in)    :: count
   integer(c_int)               ,intent(in)    :: root
   integer(c_int)                              :: error 
   error = Epetra_SerialComm_Broadcast_Int(this%SerialComm_id,MyVals,count,root)
  end subroutine

  !subroutine broadcast_long(this,MyVals,count,root)
  ! class(epetra_serialcomm)     ,intent(in)    :: this
  ! integer(c_long),dimension(:) ,intent(inout) :: MyVals
  ! integer(c_int)               ,intent(in)    :: count
  ! integer(c_int)               ,intent(in)    :: root
  ! integer(c_int)                              :: error 
  ! error = Epetra_SerialComm_Broadcast_Long(this%SerialComm_id,MyVals,count,root)
  !end subroutine
 
  subroutine broadcast_char(this,MyVals,count,root)
   class(epetra_serialcomm)           ,intent(in)    :: this
   character(kind=c_char),dimension(:),intent(inout) :: MyVals
   integer(c_int)                     ,intent(in)    :: count
   integer(c_int)                     ,intent(in)    :: root
   integer(c_int)                                    :: error 
   error = Epetra_SerialComm_Broadcast_Char(this%SerialComm_id,MyVals,count,root)
  end subroutine

  subroutine finalize(this)
    type(epetra_serialcomm) :: this
    call Epetra_SerialComm_Destroy( this%SerialComm_id ) 
    deallocate(this%SerialComm_id)
  end subroutine

  subroutine force_finalization(this)
    class(epetra_serialcomm) ,intent(inout) :: this
    if (associated(this%SerialComm_id)) then
      call finalize(this) 
      deallocate(this%SerialComm_id)
    else
      print *,' finalization for epetra_serialcomm received object with unassociated SerialComm_id'
    end if
  end subroutine
end module 

