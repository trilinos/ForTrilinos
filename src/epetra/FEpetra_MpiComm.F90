module FEpetra_MpiComm
#ifdef HAVE_MPI
  use ForTrilinos_enums ,only: FT_Epetra_MpiComm_ID_t,ForTrilinos_Universal_ID_t
  use FEpetra_Comm      ,only: epetra_comm
  use iso_c_binding     ,only: c_int,c_double,c_long,c_char
  use forepetra
  implicit none
  
  type ,extends(epetra_comm) :: epetra_mpicomm
    private
    type(FT_Epetra_MpiComm_ID_t), pointer :: MpiComm_id  
  contains
    !Constructor
    procedure         :: clone
    !Developers only
    procedure         :: get_EpetraMpiComm_ID
    procedure ,nopass :: alias_EpetraMpiComm_ID
    procedure         :: generalize
    procedure         :: MpiComm_assign => assign_to_epetra_MpiComm 
    procedure         :: Comm_assign => assign_to_epetra_Comm
    !Barrier Method
    procedure         :: barrier
    !Broadcast Method
    procedure,private :: broadcast_double
    procedure,private :: broadcast_int
    procedure,private :: broadcast_long
    procedure,private :: broadcast_char
    generic :: broadcast=>broadcast_double,broadcast_int,broadcast_long,broadcast_char
    !Memory Management 
    procedure :: force_finalization 
    final :: finalize
  end type
  
 interface epetra_mpicomm ! constructors
   procedure from_scratch,duplicate,from_struct
 end interface
  
contains

 type(FT_Epetra_MpiComm_ID_t) function from_struct(id)
     type(FT_Epetra_MpiComm_ID_t) ,intent(in) :: id
     from_struct = id
  end function
 
  ! Original C++ prototype:
  ! Epetra_MpiComm(MPI_Comm comm);
  ! CTrilinos prototype:
  ! CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Create ( MPI_Comm comm );

  type(FT_Epetra_MpiComm_ID_t) function from_scratch(comm)
    include 'mpif.h'
    MPI_Comm ,intent(in) :: comm
    from_scratch = Epetra_MpiComm_Create(comm)
  end function

  ! Original C++ prototype:
  ! Epetra_MpiComm(const Epetra_MpiComm & Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Duplicate ( CT_Epetra_MpiComm_ID_t CommID );

  type(FT_Epetra_MpiComm_ID_t) function duplicate(original)
    type(epetra_mpicomm) ,intent(in) :: original
    duplicate = Epetra_MpiComm_Duplicate(original%MpiComm_id)
  end function

  function clone(this)
    class(epetra_mpicomm)    ,intent(in)  :: this
    class(epetra_comm)       ,allocatable :: clone
    allocate(epetra_mpicomm :: clone)
    clone = Epetra_MpiComm_Clone(this%MpiComm_id)
  end function

 type(FT_Epetra_MpiComm_ID_t) function get_EpetraMpiComm_ID(this)
   class(epetra_mpicomm) ,intent(in) :: this
   if (associated(this%MpiComm_id)) then
    get_EpetraMpiComm_ID=this%MpiComm_id
   else
    stop 'get_EpetraMpiComm_ID: MpiComm_id is unassociated'
   end if
  end function
 
  type(FT_Epetra_MpiComm_ID_t) function alias_EpetraMpiComm_ID(generic_id)
    use iso_c_binding     ,only: c_loc
    use ForTrilinos_enums ,only: FT_Epetra_MpiComm_ID,ForTrilinos_Universal_ID_t
    use ForTrilinos_utils ,only: CT_Alias
    type(Fortrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(Fortrilinos_Universal_ID_t) ,pointer    :: alias_id
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_MpiComm_ID))
    alias_EpetraMpiComm_ID=degeneralize_EpetraMpiComm(c_loc(alias_id))
    deallocate(alias_id)
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding ,only : c_loc
   class(epetra_mpicomm) ,intent(in) ,target :: this
   generalize = generalize_all( c_loc(this%MpiComm_id) )
   ! ____ Use for ForTrilinos function implementation ______
  
   ! ____ Use for CTrilinos function implementation ______
   ! class(epetra_mpicomm) ,intent(in) ,target :: this
   ! generalize = Epetra_MpiComm_Generalize ( this%MpiComm_id )
   ! ____ Use for CTrilinos function implementation ______
  end function

  type(FT_Epetra_MpiComm_ID_t) function degeneralize_EpetraMpiComm(generic_id) bind(C)
   ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_MpiComm_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                     ,value   :: generic_id
    type(FT_Epetra_MpiComm_ID_t) ,pointer :: local_ptr
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraMpiComm = local_ptr
   ! ____ Use for ForTrilinos function implementation ______
  
   ! ____ Use for CTrilinos function implementation ______
   ! type(Fortrilinos_Universal_ID_t) ,intent(in) :: generic_id
   ! degeneralize_EpetraMpiComm = Epetra_MpiComm_Degeneralize( generic_id )
   ! ____ Use for CTrilinos function implementation ______
  end function

  subroutine assign_to_epetra_MpiComm(lhs,rhs)
    class(epetra_mpicomm)     ,intent(inout) :: lhs
    type(FT_Epetra_MpiComm_ID_t) ,intent(in)    :: rhs
    type(FT_Epetra_Comm_ID_t)                   :: test  ! test line
    allocate( lhs%MpiComm_id, source=rhs)
    PRINT *,'assignment after mpi  part'
    PRINT *,lhs%MpiComm_id%table,lhs%MpiComm_id%index,lhs%MpiComm_id%is_const
    call lhs%set_EpetraComm_ID(lhs%alias_EpetraComm_ID(lhs%generalize()))
    test=lhs%get_EpetraComm_ID()  ! test line
    PRINT *,test%table,test%index,test%is_const
  end subroutine

  subroutine assign_to_epetra_Comm(lhs,rhs)
    class(epetra_mpicomm) ,intent(inout) :: lhs
    class(epetra_comm)       ,intent(in)    :: rhs
    select type(rhs)
      class is (epetra_mpicomm)
        allocate(lhs%MpiComm_id,source=rhs%MpiComm_id)
        call lhs%set_EpetraComm_ID(lhs%alias_EpetraComm_ID(lhs%generalize()))
     class default
        stop 'assign_to_epetra_Comm: unsupported class'
     end select
  end subroutine

  subroutine barrier(this)
   class(epetra_mpicomm) ,intent(in) :: this
   call Epetra_MpiComm_Barrier(this%MpiComm_id)
  end subroutine

  subroutine broadcast_double(this,My_Vals,count,root)
    class(epetra_mpicomm)       ,intent(in)    :: this
    real(c_double) ,dimension(:),intent(inout) :: My_Vals
    integer(c_int)              ,intent(in)    :: count,root
    integer(c_int)                             :: error_out
    error_out=Epetra_MpiComm_Broadcast_Double(this%MpiComm_id,My_Vals,count,root)
  end subroutine

  subroutine broadcast_int(this,My_Vals,count,root)
    class(epetra_mpicomm)       ,intent(in)    :: this
    integer(c_int) ,dimension(:),intent(inout) :: My_Vals
    integer(c_int)              ,intent(in)    :: count,root
    integer(c_int)                             :: error_out
    error_out=Epetra_MpiComm_Broadcast_Int(this%MpiComm_id,My_Vals,count,root)
  end subroutine

  subroutine broadcast_long(this,My_Vals,count,root)
    class(epetra_mpicomm)       ,intent(in)    :: this
    real(c_long)   ,dimension(:),intent(inout) :: My_Vals
    integer(c_int)              ,intent(in)    :: count,root
    integer(c_int)                             :: error_out
    error_out=Epetra_MpiComm_Broadcast_Long(this%MpiComm_id,My_Vals,count,root)
  end subroutine

  subroutine broadcast_char(this,My_Vals,count,root)
    class(epetra_mpicomm)              ,intent(in)    :: this
    character(kind=c_char),dimension(:),intent(inout) :: My_Vals
    integer(c_int)                     ,intent(in)    :: count,root
    integer(c_int)                                    :: error_out
    error_out=Epetra_MpiComm_Broadcast_Char(this%MpiComm_id,My_Vals,count,root)
  end subroutine

  subroutine finalize(this)
    type(epetra_mpicomm) :: this
    call Epetra_MpiComm_Destroy ( this%MpiComm_id ) 
    deallocate(this%MpiComm_id)
  end subroutine

  subroutine force_finalize(this)
    class(epetra_mpicomm) ,intent(inout) :: this
    if (associated(this%MpiComm_id)) then
      call finalize(this) 
      deallocate(this%MpiComm_id)
    else
      print *,'finalization for epetra_mpicomm received object with unassociated MpiComm_id'
    end if
  end subroutine
#endif
end module 
