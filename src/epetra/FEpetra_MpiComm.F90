module FEpetra_MpiComm
#include "ForTrilinos_config.h"
#ifdef HAVE_MPI
#include "mpif.h"
  use ForTrilinos_enums ,only: FT_Epetra_MpiComm_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use FEpetra_Comm      ,only: epetra_comm
  use iso_c_binding     ,only: c_int,c_double,c_long,c_char
  use forepetra
  implicit none
  private               ! Hide everything by default
  public :: epetra_mpicomm ! Expose type/methods

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
    procedure         :: SerialComm_assign
    !Barrier Method
    procedure         :: barrier
    !Broadcast Method
    procedure         :: broadcast_double
    procedure         :: broadcast_int
    procedure         :: broadcast_long
    procedure         :: broadcast_char
    !generic :: broadcast=>broadcast_double,broadcast_int,broadcast_char
    !generic :: broadcast=>broadcast_double,broadcast_int,broadcast_long,broadcast_char
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
 !   MPI_Comm ,intent(in) :: comm
    integer ,intent(in) :: comm
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

  !function clone(this)
  !  class(epetra_mpicomm)    ,intent(in)  :: this
  !  class(epetra_comm)       ,allocatable :: clone
  !  allocate(epetra_mpicomm :: clone)
  !  clone = Epetra_MpiComm_Clone(this%MpiComm_id)
  !end function

  function clone(this)
    class(epetra_mpicomm) ,intent(in)  :: this
    class(epetra_comm)       ,allocatable :: clone
    class(epetra_comm)       ,allocatable :: clone_temp
    type(FT_Epetra_MpiComm_ID_t) :: test
    type(FT_Epetra_Comm_ID_t) :: test1
    allocate(epetra_mpicomm :: clone)
    allocate(epetra_mpicomm :: clone_temp)
    clone_temp = Epetra_MpiComm_Clone(this%MpiComm_id)
   ! test = clone_temp%MpiComm_id
   ! print *,'clone_temp%mpicom',test%table,test%index
    test1 = clone_temp%get_EpetraComm_ID()
    print *,'clone_temp%comm',test1%table,test1%index
    clone=epetra_mpicomm(alias_EpetraMpiComm_ID(clone_temp%generalize_EpetraComm()))
    !test = clone%MpiComm_id
   ! test = clone%get_EpetraMpiComm_ID()
   ! print *,'clone%mpicomm',test%table,test%index
    test1 = clone%get_EpetraComm_ID()
    print *,'clone%comm',test1%table,test1%index
    call clone_temp%force_finalization_EpetraComm()
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
    use iso_c_binding        ,only: c_loc
    use ForTrilinos_enums    ,only: FT_Epetra_MpiComm_ID,ForTrilinos_Universal_ID_t
    use ForTrilinos_table_man,only: CT_Alias
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

  subroutine SerialComm_assign(lhs,rhs)
    use ForTrilinos_enums, only : FT_Epetra_SerialComm_ID_t
    class(epetra_mpicomm),      intent(inout) :: lhs
    type(FT_Epetra_SerialComm_ID_t), intent(in) :: rhs 
    print *, 'SerialComm_assign in FEpetra_MpiComm no-op'
  end subroutine


  subroutine assign_to_epetra_Comm(lhs,rhs)
    class(epetra_mpicomm) ,intent(inout) :: lhs
    class(epetra_comm)       ,intent(in)    :: rhs
    select type(rhs)
      class is (epetra_mpicomm)
        allocate(lhs%MpiComm_id,source=alias_EpetraMpiComm_ID(rhs%generalize()))
        call lhs%set_EpetraComm_ID(lhs%alias_EpetraComm_ID(lhs%generalize()))
     class default
        stop 'assign_to_epetra_Comm: unsupported class'
     end select
  end subroutine

  subroutine barrier(this)
   class(epetra_mpicomm) ,intent(in) :: this
   call Epetra_MpiComm_Barrier(this%MpiComm_id)
  end subroutine

  subroutine broadcast_double(this,MyVals,count,root)
    class(epetra_mpicomm)       ,intent(in)    :: this
    real(c_double) ,dimension(:),intent(inout) :: MyVals
    integer(c_int)              ,intent(in)    :: count,root
    integer(c_int)                             :: error_out
    error_out=Epetra_MpiComm_Broadcast_Double(this%MpiComm_id,MyVals,count,root)
  end subroutine

  subroutine broadcast_int(this,MyVals,count,root)
    class(epetra_mpicomm)       ,intent(in)    :: this
    integer(c_int) ,dimension(:),intent(inout) :: MyVals
    integer(c_int)              ,intent(in)    :: count,root
    integer(c_int)                             :: error_out
    error_out=Epetra_MpiComm_Broadcast_Int(this%MpiComm_id,MyVals,count,root)
  end subroutine

  subroutine broadcast_long(this,MyVals,count,root)
    class(epetra_mpicomm)       ,intent(in)    :: this
    integer(c_long),dimension(:),intent(inout) :: MyVals
    integer(c_int)              ,intent(in)    :: count,root
    integer(c_int)                             :: error_out
    error_out=Epetra_MpiComm_Broadcast_Long(this%MpiComm_id,MyVals,count,root)
  end subroutine

  subroutine broadcast_char(this,MyVals,count,root)
    class(epetra_mpicomm)              ,intent(in)    :: this
    character(kind=c_char),dimension(:),intent(inout) :: MyVals
    integer(c_int)                     ,intent(in)    :: count,root
    integer(c_int)                                    :: error_out
    error_out=Epetra_MpiComm_Broadcast_Char(this%MpiComm_id,MyVals,count,root)
  end subroutine

  subroutine finalize(this)
    type(epetra_mpicomm) :: this
    call Epetra_MpiComm_Destroy ( this%MpiComm_id ) 
    deallocate(this%MpiComm_id)
  end subroutine

  subroutine force_finalization(this)
    class(epetra_mpicomm) ,intent(inout) :: this
    if (associated(this%MpiComm_id)) then
      call finalize(this) 
    else
      print *,'finalization for epetra_mpicomm received object with unassociated MpiComm_id'
    end if
  end subroutine
#endif
end module 
