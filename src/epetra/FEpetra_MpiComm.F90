module FEpetra_MpiComm
#include "ForTrilinos_config.h"
#ifdef HAVE_MPI
!#include "mpif.h"
  use ForTrilinos_enums ,only: FT_Epetra_MpiComm_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_external_utils
  use ForTrilinos_hermetic_new, only : hermetic
  use FEpetra_Comm      ,only: Epetra_Comm
  use iso_c_binding     ,only: c_int,c_double,c_long,c_char
  use forepetra
  implicit none
  private               ! Hide everything by default
  public :: Epetra_MpiComm ! Expose type/methods

  type ,extends(Epetra_Comm) :: Epetra_MpiComm
    private
    type(FT_Epetra_MpiComm_ID_t), pointer :: MpiComm_id  
  contains
    !Constructor
    procedure         :: clone
    !Developers only
    procedure, private::remote_dealloc_EpetraMpiComm
    procedure         :: get_EpetraMpiComm_ID
    procedure ,nopass :: alias_EpetraMpiComm_ID
    procedure         :: generalize
    !Barrier Method
    procedure         :: barrier
    !Broadcast Method
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
  end type
  
 interface Epetra_MpiComm ! constructors
   procedure from_scratch,duplicate,from_struct
 end interface
  
contains

 type(Epetra_MpiComm) function from_struct(id)
   type(FT_Epetra_MpiComm_ID_t) ,intent(in) :: id
   from_struct%hermetic = hermetic()
   from_struct%MpiComm_id = id
   call from_struct%set_EpetraComm_ID(from_struct%alias_EpetraComm_ID(from_struct%generalize()))
  end function
 
  ! Original C++ prototype:
  ! Epetra_MpiComm(MPI_Comm comm);
  ! CTrilinos prototype:
  ! CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Create ( MPI_Comm comm );

  type(Epetra_MpiComm) function from_scratch(comm)
    integer(c_int) ,intent(in) :: comm
    from_scratch%hermetic = hermetic()
    from_scratch%MpiComm_id = Epetra_MpiComm_Fortran_Create(comm)
    call from_scratch%set_EpetraComm_ID(from_scratch%alias_EpetraComm_ID(from_scratch%generalize()))
  end function

  ! Original C++ prototype:
  ! Epetra_MpiComm(const Epetra_MpiComm & Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Duplicate ( CT_Epetra_MpiComm_ID_t CommID );

  type(Epetra_MpiComm) function duplicate(original)
   type(Epetra_MpiComm) ,intent(in) :: original
   duplicate%hermetic = hermetic()
   duplicate%MpiComm_id = Epetra_MpiComm_Duplicate(original%MpiComm_id)
   call duplicate%set_EpetraComm_ID(duplicate%alias_EpetraComm_ID(duplicate%generalize()))
  end function

  function clone(this)
    class(Epetra_MpiComm)    ,intent(in)  :: this
    class(Epetra_Comm)       ,allocatable :: clone
    allocate(Epetra_MpiComm :: clone)
    clone=Epetra_MpiComm(Epetra_MpiComm_Clone(this%MpiComm_id))
  end function

  !function clone(this)
   ! class(Epetra_MpiComm) ,intent(in)  :: this
   ! class(Epetra_Comm)       ,allocatable :: clone
    !class(Epetra_Comm)       ,allocatable :: clone_temp
    !type(FT_Epetra_MpiComm_ID_t) :: test
    !type(FT_Epetra_Comm_ID_t) :: test1
    !allocate(Epetra_MpiComm :: clone)
    !allocate(Epetra_MpiComm :: clone_temp)
    !clone_temp = Epetra_MpiComm_Clone(this%MpiComm_id)
   !! test = clone_temp%MpiComm_id
   !! print *,'clone_temp%mpicom',test%table,test%index
    !test1 = clone_temp%get_EpetraComm_ID()
    !print *,'clone_temp%comm',test1%table,test1%index
    !clone=Epetra_MpiComm(alias_EpetraMpiComm_ID(clone_temp%generalize_EpetraComm()))
    !!test = clone%MpiComm_id
   !! test = clone%get_EpetraMpiComm_ID()
   !! print *,'clone%mpicomm',test%table,test%index
    !test1 = clone%get_EpetraComm_ID()
    !print *,'clone%comm',test1%table,test1%index
    !call clone_temp%force_finalization_EpetraComm()
  !end function


 type(FT_Epetra_MpiComm_ID_t) function get_EpetraMpiComm_ID(this)
   class(Epetra_MpiComm) ,intent(in) :: this
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
   class(Epetra_MpiComm) ,intent(in) ,target :: this
   generalize = generalize_all( c_loc(this%MpiComm_id) )
   ! ____ Use for ForTrilinos function implementation ______
  
   ! ____ Use for CTrilinos function implementation ______
   ! class(Epetra_MpiComm) ,intent(in) ,target :: this
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

  subroutine barrier(this)
   class(Epetra_MpiComm) ,intent(in) :: this
   call Epetra_MpiComm_Barrier(this%MpiComm_id)
  end subroutine

  subroutine broadcast_double(this,MyVals,count,root)
    class(Epetra_MpiComm)       ,intent(in)    :: this
    real(c_double) ,dimension(:),intent(inout) :: MyVals
    integer(c_int)              ,intent(in)    :: count,root
    integer(c_int)                             :: error_out
    error_out=Epetra_MpiComm_Broadcast_Double(this%MpiComm_id,MyVals,count,root)
  end subroutine

  subroutine broadcast_int(this,MyVals,count,root)
    class(Epetra_MpiComm)       ,intent(in)    :: this
    integer(c_int) ,dimension(:),intent(inout) :: MyVals
    integer(c_int)              ,intent(in)    :: count,root
    integer(c_int)                             :: error_out
    error_out=Epetra_MpiComm_Broadcast_Int(this%MpiComm_id,MyVals,count,root)
  end subroutine

  subroutine broadcast_long(this,MyVals,count,root)
    class(Epetra_MpiComm)       ,intent(in)    :: this
    integer(c_long),dimension(:),intent(inout) :: MyVals
    integer(c_int)              ,intent(in)    :: count,root
    integer(c_int)                             :: error_out
    error_out=Epetra_MpiComm_Broadcast_Long(this%MpiComm_id,MyVals,count,root)
  end subroutine

  subroutine broadcast_char(this,MyVals,count,root)
    class(Epetra_MpiComm)              ,intent(in)    :: this
    character(kind=c_char),dimension(:),intent(inout) :: MyVals
    integer(c_int)                     ,intent(in)    :: count,root
    integer(c_int)                                    :: error_out
    error_out=Epetra_MpiComm_Broadcast_Char(this%MpiComm_id,MyVals,count,root)
  end subroutine

  subroutine gather_double(this,MyVals,AllVals,count)
   class(Epetra_MpiComm)     ,intent(in)    :: this
   !real(c_double), dimension(:) ,intent(inout) :: MyVals
   !real(c_double), dimension(:) ,intent(inout) :: AllVals
   real(c_double), dimension(:)  :: MyVals
   real(c_double), dimension(:)  :: AllVals
   integer(c_int)               ,intent(in)    :: count
   integer(c_int)                              :: error
   error = Epetra_MpiComm_GatherAll_Double(this%MpiComm_id,MyVals,AllVals,count)
  end subroutine

  integer(c_int) function MyPID(this)
   class(Epetra_MpiComm)     , intent(in) :: this
   MyPID=Epetra_MpiComm_MyPID(this%MpiComm_id)
  end function

  integer(c_int) function NumProc(this)
   class(Epetra_MpiComm)     , intent(in) :: this
   NumProc=Epetra_MpiComm_NumProc(this%MpiComm_id)
  end function

  subroutine remote_dealloc_EpetraMpiComm(this)
    class(Epetra_MpiComm) ,intent(inout) :: this
    call this%Epetra_Comm%remote_dealloc_EpetraComm()
    call Epetra_MpiComm_Destroy(this%MpiComm_id)
  end subroutine
#endif
end module 
