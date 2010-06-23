!*********************************************************************
! ForTrilinos: Object-Oriented Fortran 2003 interface to Trilinos
!                Copyright 2010 Sandia Corporation
!
! Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
! the U.S. Government retains certain rights in this software.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright
!    notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright
!    notice, this list of conditions and the following disclaimer in the
!    documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the Corporation nor the names of the
!    contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! Questions? Contact Karla Morris  (knmorri@sandia.gov) or
!                    Damian Rouson (rouson@sandia.gov)
!*********************************************************************


module FEpetra_BlockMap
  use ForTrilinos_enums ,only: FT_Epetra_Comm_ID_t,FT_Epetra_BlockMap_ID_t,FT_Epetra_Map_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_hermetic,only:hermetic
  use ForTrilinos_universal,only:universal
  use FEpetra_Comm  ,only: Epetra_Comm
  use iso_c_binding ,only: c_int
  use forepetra
  implicit none
  private                   ! Hide everything by default
  public :: Epetra_BlockMap ! Expose type/constructors/methods

  type ,extends(universal)      :: Epetra_BlockMap !"shell"
    private
    type(FT_Epetra_BlockMap_ID_t) :: BlockMap_id 
  contains
     !Developers only
     procedure         :: invalidate_id => invalidate_EpetraBlockMap_ID
     procedure         :: ctrilinos_delete => ctrilinos_delete_EpetraBlockMap
     procedure         :: get_EpetraBlockMap_ID 
     procedure ,nopass :: alias_EpetraBlockMap_ID
     procedure         :: generalize 
     !Local/Global ID accessor methods
     !Size and dimension acccessor functions
     procedure         :: NumGlobalElements
     procedure         :: NumMyElements
     procedure         :: MyGlobalElements
     procedure         :: ElementSize_Const
     procedure         :: ElementSize_LID
     generic :: ElementSize=>ElementSize_Const,ElementSize_LID
     !Miscellaneous boolean tests
     procedure         :: LinearMap
     procedure         :: DistributedGlobal
     !Array accessor functions
     !Miscellaneous
  end type

   interface Epetra_BlockMap ! constructors
     module procedure from_scratch,duplicate,from_struct,from_scratch_linear,from_scratch_arbitrary,from_scratch_variable
   end interface
 
contains
  type(Epetra_BlockMap) function from_struct(id)
     type(FT_Epetra_BlockMap_ID_t) ,intent(in) :: id
     from_struct%BlockMap_id = id
     call from_struct%register_self()
  end function
 
  ! Original C++ prototype:
  ! Epetra_BlockMap(int NumGlobalElements, int ElementSize, int IndexBase, const Epetra_Comm& Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create ( int NumGlobalElements, int ElementSize, int IndexBase, CT_Epetra_Comm_ID_t CommID );

  type(Epetra_BlockMap) function from_scratch(Num_GlobalElements,ElementSize,IndexBase,comm)
   !use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_Map_ID_t
    integer(c_int) ,intent(in) :: Num_GlobalElements
    integer(c_int) ,intent(in) :: ElementSize
    integer(c_int) ,intent(in) :: IndexBase
    class(Epetra_Comm)         :: comm
    type(FT_Epetra_BlockMap_ID_t) :: from_scratch_id
    from_scratch_id = Epetra_BlockMap_Create(Num_GlobalElements,ElementSize,IndexBase,comm%get_EpetraComm_ID())
    from_scratch = from_struct(from_scratch_id)
  end function

! Original C++ prototype:
  ! Epetra_BlockMap(int NumGlobalElements, int NumMyElements, int ElementSize, int IndexBase, const Epetra_Comm& Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Linear ( int NumGlobalElements, int NumMyElements, int
  ! ElementSize, int IndexBase, CT_Epetra_Comm_ID_t CommID );

  type(Epetra_BlockMap) function from_scratch_linear(Num_GlobalElements,Num_MyElements,Element_Size,IndexBase,comm)
   !use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_Map_ID_t
    integer(c_int) ,intent(in) :: Num_GlobalElements
    integer(c_int) ,intent(in) :: Num_MyElements
    integer(c_int) ,intent(in) :: Element_Size
    integer(c_int) ,intent(in) :: IndexBase
    class(Epetra_Comm)         :: comm
    type(FT_Epetra_BlockMap_ID_t) :: from_scratch_linear_id
    from_scratch_linear_id = Epetra_BlockMap_Create_Linear(Num_GlobalElements,Num_MyElements,Element_Size,IndexBase,comm%get_EpetraComm_ID())
    from_scratch_linear = from_struct(from_scratch_linear_id)
  end function

!Original C++ prototype:
  ! Epetra_BlockMap(int NumGlobalElements, int NumMyElements, const int *MyGlobalElements, int ElementSize, int IndexBae, const Epetra_Comm& Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Arbitrary ( int NumGlobalElements, int NumMyElements,
  !const int * MyGlobalElements, int ElementSize, int IndexBase, CT_Epetra_Comm_ID_t CommID );

  type(Epetra_BlockMap) function from_scratch_arbitrary(Num_GlobalElements,Num_MyElements,My_GlobalElements,Element_Size,IndexBase,comm)
   !use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_Map_ID_t
    integer(c_int) ,intent(in) :: Num_GlobalElements
    integer(c_int) ,intent(in) :: Num_MyElements
    integer(c_int) ,intent(in) ,dimension(:) :: My_GlobalElements  
    integer(c_int) ,intent(in) :: Element_Size
    integer(c_int) ,intent(in) :: IndexBase
    class(Epetra_Comm)         :: comm
    type(FT_Epetra_BlockMap_ID_t) :: from_scratch_arbitrary_id
    from_scratch_arbitrary_id = Epetra_BlockMap_Create_Arbitrary(Num_GlobalElements,Num_MyElements,My_GlobalElements,Element_Size,IndexBase,comm%get_EpetraComm_ID())
    from_scratch_arbitrary = from_struct(from_scratch_arbitrary_id)
  end function

! Original C++ prototype:
  ! Epetra_BlockMap(int NumGlobalElements, int NumMyElements, const int *MyGlobalElements, const int *ElementSizeList,
  !int IndexBase, const Epetra_Comm& Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Variable ( int NumGlobalElements, int NumMyElements,
  !const int * MyGlobalElements, const int * ElementSizeList, int IndexBase, CT_Epetra_Comm_ID_t CommID );

  type(Epetra_BlockMap) function from_scratch_variable(Num_GlobalElements,Num_MyElements,My_GlobalElements,Element_SizeList,IndexBase,comm)
   !use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_Map_ID_t
    integer(c_int) ,intent(in) :: Num_GlobalElements
    integer(c_int) ,intent(in) :: Num_MyElements
    integer(c_int) ,intent(in) ,dimension(:) :: My_GlobalElements  
    integer(c_int) ,intent(in) ,dimension(:) :: Element_SizeList    
    integer(c_int) ,intent(in) :: IndexBase
    class(Epetra_Comm)         :: comm
    type(FT_Epetra_BlockMap_ID_t) :: from_scratch_variable_id
    from_scratch_variable_id = Epetra_BlockMap_Create_Variable(Num_GlobalElements,Num_MyElements,My_GlobalElements,Element_SizeList,IndexBase,comm%get_EpetraComm_ID())
    from_scratch_variable = from_struct(from_scratch_variable_id)
  end function

  ! Original C++ prototype:
  ! Epetra_BlockMap(const Epetra_BlockMap& map);
  ! CTrilinos prototype:
  ! CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Duplicate ( CT_Epetra_BlockMap_ID_t mapID );

  type(Epetra_BlockMap) function duplicate(this)
    type(Epetra_BlockMap) ,intent(in) :: this 
    type(FT_Epetra_BlockMap_ID_t) :: duplicate_id
    duplicate_id = Epetra_BlockMap_Duplicate(this%BlockMap_id)
    duplicate = from_struct(duplicate_id)
  end function

  type(FT_Epetra_BlockMap_ID_t) function get_EpetraBlockMap_ID(this)
    class(Epetra_BlockMap) ,intent(in) :: this 
    get_EpetraBlockMap_ID=this%BlockMap_id
  end function
  
  type(FT_Epetra_BlockMap_ID_t) function alias_EpetraBlockMap_ID(generic_id)
    use ForTrilinos_table_man
    use ForTrilinos_enums ,only: ForTrilinos_Universal_ID_t,FT_Epetra_BlockMap_ID
    use iso_c_binding     ,only: c_loc
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,pointer    :: alias_id
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_BlockMap_ID))
    alias_EpetraBlockMap_ID=degeneralize_EpetraBlockMap(c_loc(alias_id))
    deallocate(alias_id)
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding     ,only: c_loc
   class(Epetra_BlockMap) ,intent(in) ,target :: this
   generalize = generalize_all(c_loc(this%BlockMap_ID))
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   !class(Epetra_BlockMap) ,intent(in) ,target :: this
   !generalize = Epetra_BlockMap_Generalize ( this%BlockMap_id)
   ! ____ Use for CTrilinos function implementation ______
  end function

 type(FT_Epetra_BlockMap_ID_t) function degeneralize_EpetraBlockMap(generic_id) bind(C)
   ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_BlockMap_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                   ,value   :: generic_id
    type(FT_Epetra_BlockMap_ID_t) ,pointer :: local_ptr
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraBlockMap = local_ptr
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   !type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
   !degeneralize_EpetraBlockMap = Epetra_BlockMap_Degeneralize(generic_id)
   ! ____ Use for CTrilinos function implementation ______
  end function
 
  integer(c_int) function NumGlobalElements(this)
    class(Epetra_BlockMap) ,intent(in) :: this
    NumGlobalElements=Epetra_BlockMap_NumGlobalElements(this%BlockMap_id)
  end function 

  integer(c_int) function NumMyElements(this)
    class(Epetra_BlockMap) ,intent(in) :: this
    NumMyElements=Epetra_BlockMap_NumMyElements(this%BlockMap_id)
  end function 
 
  function MyGlobalElements(this) result(MyGlobalElementsList)
    class(Epetra_BlockMap)     ,intent(in)    :: this
    integer(c_int),dimension(:),allocatable   :: MyGlobalElementsList
    integer(c_int)                            :: junk
    allocate(MyGlobalElementsList(this%NumMyElements()))
    junk=Epetra_BlockMap_MyGlobalElements_Fill(this%BlockMap_id,MyGlobalElementsList)
  end function 

  integer(c_int) function ElementSize_Const(this)
    class(Epetra_BlockMap) ,intent(in) :: this
    ElementSize_Const=Epetra_BlockMap_ElementSize_Const(this%BlockMap_id)
  end function 

  integer(c_int) function ElementSize_LID(this,L_ID)
    class(Epetra_BlockMap) ,intent(in) :: this
    integer(c_int)         ,intent(in) :: L_ID
    ElementSize_LID=Epetra_BlockMap_ElementSize(this%BlockMap_id,L_ID)
  end function 

  logical function LinearMap(this)
    use ForTrilinos_enums, only:FT_boolean_t,FT_FALSE,FT_TRUE
    class(Epetra_BlockMap),intent(in) :: this
    integer(FT_boolean_t) :: LinearMap_out
    LinearMap_out=Epetra_BlockMap_LinearMap(this%BlockMap_id)
    if (LinearMap_out==FT_FALSE) LinearMap=.false.
    if (LinearMap_out==FT_TRUE) LinearMap=.true.
  end function

  logical function DistributedGlobal(this)
    use ForTrilinos_enums, only:FT_boolean_t,FT_FALSE,FT_TRUE
    class(Epetra_BlockMap),intent(in) :: this
    integer(FT_boolean_t) :: DistributedGlobal_out
    DistributedGlobal_out=Epetra_BlockMap_DistributedGlobal(this%BlockMap_id)
    if (DistributedGlobal_out==FT_FALSE) DistributedGlobal =.false.
    if (DistributedGlobal_out==FT_TRUE) DistributedGlobal=.true.
  end function

  subroutine invalidate_EpetraBlockMap_ID(this)
    class(Epetra_BlockMap),intent(inout) :: this
    this%BlockMap_id%table = FT_Invalid_ID
    this%BlockMap_id%index = FT_Invalid_Index 
    this%BlockMap_id%is_const = FT_FALSE
  end subroutine

  subroutine ctrilinos_delete_EpetraBlockMap(this)
    class(Epetra_BlockMap),intent(inout) :: this
    call Epetra_BlockMap_Destroy( this%BlockMap_id ) 
  end subroutine

end module 

