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


module FEpetra_Map
  use ForTrilinos_enums !,only: FT_Epetra_BlockMap_ID,FT_Epetra_Map_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use FEpetra_Comm       ,only: Epetra_Comm
  use FEpetra_BlockMap   ,only: Epetra_BlockMap
  use iso_c_binding      ,only: c_int
  use forepetra
  private                     ! Hide everything by default
  public :: Epetra_Map!,Epetra_BlockMap ! Expose type/constructors/methods
  implicit none

  type, extends(Epetra_BlockMap)      :: Epetra_Map !"shell"
    private
    type(FT_Epetra_Map_ID_t) ,pointer :: map_id => null()
  contains
     !Developers only
     procedure         :: get_EpetraMap_ID 
     procedure ,nopass :: alias_EpetraMap_ID
     procedure         :: generalize 
     procedure         :: assign_to_Epetra_Map
     procedure         :: EpetraMap_assign_to_EpetraMap  
     generic :: assignment(=) => assign_to_Epetra_Map,EpetraMap_assign_to_EpetraMap
     !Memory Management
     procedure         :: force_finalization 
     final :: finalize
  end type

   interface Epetra_Map ! constructors
     module procedure from_scratch,duplicate,from_struct,from_scratch_linear,from_scratch_arbitrary
   end interface
 
contains

  type(FT_Epetra_Map_ID_t) function from_struct(id)
     type(FT_Epetra_Map_ID_t) ,intent(in) :: id
     from_struct = id
  end function

  ! Original C++ prototype:
  ! Epetra_Map(int NumGlobalElements, int IndexBase, const Epetra_Comm& Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_Map_ID_t Epetra_Map_Create ( int NumGlobalElements, int IndexBase, CT_Epetra_Comm_ID_t CommID );

  type(FT_Epetra_Map_ID_t) function from_scratch(Num_GlobalElements,IndexBase,comm)
    use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_Map_ID_t
    integer(c_int) ,intent(in) :: Num_GlobalElements
    integer(c_int) ,intent(in) :: IndexBase
    class(Epetra_Comm)         :: comm
    from_scratch = Epetra_Map_Create(Num_GlobalElements,IndexBase,comm%get_EpetraComm_ID())
  end function

! Original C++ prototype:
  ! Epetra_Map(int NumGlobalElements, int NumMyElements, int IndexBase, const Epetra_Comm& Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_Map_ID_t Epetra_Map_Create_Linear ( int NumGlobalElements, int NumMyElements, int IndexBase,
  ! CT_Epetra_Comm_ID_t CommID );

  type(FT_Epetra_Map_ID_t) function from_scratch_linear(Num_GlobalElements,Num_MyElements,IndexBase,comm)
    use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_Map_ID_t
    integer(c_int) ,intent(in) :: Num_GlobalElements
    integer(c_int) ,intent(in) :: Num_MyElements
    integer(c_int) ,intent(in) :: IndexBase
    class(Epetra_Comm)         :: comm
    from_scratch_linear = Epetra_Map_Create_Linear(Num_GlobalElements,Num_MyElements,IndexBase,comm%get_EpetraComm_ID())
  end function

! Original C++ prototype:
  ! Epetra_Map(int NumGlobalElements, int NumMyElements, const int *MyGlobalElements, int IndexBase, const
  ! Epetra_Comm& Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_Map_ID_t Epetra_Map_Create_Arbitrary ( int NumGlobalElements, int NumMyElements, const int
  !* MyGlobalElements, int IndexBase, CT_Epetra_Comm_ID_t CommID );

 type(FT_Epetra_Map_ID_t) function from_scratch_arbitrary(Num_GlobalElements,Num_MyElements,My_GlobalElements,IndexBase,comm)
    use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_Map_ID_t
    integer(c_int) ,intent(in)              :: Num_GlobalElements
    integer(c_int) ,intent(in)              :: Num_MyElements
    integer(c_int) ,intent(in) ,dimension(:),allocatable:: My_GlobalElements
    integer(c_int) ,intent(in)              :: IndexBase
    class(Epetra_Comm)                      :: comm
    from_scratch_arbitrary = Epetra_Map_Create_Arbitrary(Num_GlobalElements,Num_MyElements,My_GlobalElements,IndexBase,comm%get_EpetraComm_ID())
  end function

  ! Original C++ prototype:
  ! Epetra_Map(const Epetra_Map& map);
  ! CTrilinos prototype:
  ! CT_Epetra_Map_ID_t Epetra_Map_Duplicate ( CT_Epetra_Map_ID_t mapID );

  type(FT_Epetra_Map_ID_t) function duplicate(original)
    type(Epetra_Map) ,intent(in) :: original
    duplicate = Epetra_Map_Duplicate(original%map_id)
  end function

  type(FT_Epetra_Map_ID_t) function get_EpetraMap_ID(this)
    class(Epetra_Map) ,intent(in) :: this 
    if (associated(this%map_id)) then
     get_EpetraMap_ID=this%map_id
    else
     stop 'get_EpetraMap_ID: map_id is unassociated'
    end if
  end function
 
  type(FT_Epetra_Map_ID_t) function alias_EpetraMap_ID(generic_id)
    use iso_c_binding        ,only: c_loc
    use ForTrilinos_enums    ,only: ForTrilinos_Universal_ID_t, FT_Epetra_Map_ID
    use ForTrilinos_table_man,only: CT_Alias
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,pointer    :: alias_id
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_Map_ID))
    alias_EpetraMap_ID=degeneralize_EpetraMap(c_loc(alias_id))
    deallocate(alias_id)
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all 
   use iso_c_binding     ,only: c_loc
   class(Epetra_Map) ,intent(in) ,target :: this
   generalize =generalize_all(c_loc(this%map_id))
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   !class(Epetra_Map) ,intent(in) ,target :: this
   !generalize = Epetra_Map_Generalize ( this%map_id )
   ! ____ Use for CTrilinos function implementation ______
  
  end function
  
  type(FT_Epetra_Map_ID_t) function degeneralize_EpetraMap(generic_id) bind(C)
   ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_Map_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)              ,value   :: generic_id
    type(FT_Epetra_Map_ID_t) ,pointer :: local_ptr
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraMap = local_ptr
   ! ____ Use for ForTrilinos function implementation ______
   
   ! ____ Use for CTrilinos function implementation ______
   ! type(ForTrilinos_Universal_ID_t) ,intent(in) : generic_id
   ! degeneralize_EpetraMap = Epetra_Map_Degeneralize(generic_id)
   ! ____ Use for CTrilinos function implementation ______
  end function

  subroutine assign_to_Epetra_Map(lhs,rhs)
    class(Epetra_Map)        ,intent(inout):: lhs
    type(FT_Epetra_Map_ID_t) ,intent(in)   :: rhs
    allocate(lhs%map_id,source=rhs)
    lhs%Epetra_BlockMap=Epetra_BlockMap(lhs%alias_EpetraBlockMap_ID(lhs%generalize()))
  end subroutine
  
  subroutine EpetraMap_assign_to_EpetraMap(lhs,rhs)
    class(Epetra_Map) ,intent(inout):: lhs
    class(Epetra_Map) ,intent(in)   :: rhs
    !call Epetra_Map_Assign(lhs%map_id,rhs%map_id)
    !lhs = Epetra_Map(lhs%map_id)
    lhs = Epetra_Map(rhs%alias_EpetraMap_ID(rhs%generalize()))
  end subroutine

  subroutine finalize(this)
    type(Epetra_Map) :: this
    call Epetra_Map_Destroy( this%map_id ) 
    deallocate(this%map_id)
  end subroutine

  subroutine force_finalization(this)
    class(Epetra_Map) ,intent(inout) :: this
    call this%Epetra_BlockMap%force_finalization()
    if (associated(this%map_id)) then
      call finalize(this) 
    else
      print *,' finalization for Epetra_Map received object with unassociated map_id'
    end if
  end subroutine

end module 

