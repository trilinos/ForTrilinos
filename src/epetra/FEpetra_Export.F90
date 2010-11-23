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

module FEpetra_Export
  use ForTrilinos_enums ,only: FT_Epetra_Comm_ID_t,FT_Epetra_Export_ID_t,FT_Epetra_BlockMap_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_universal,only:universal
  use ForTrilinos_error
  use FEpetra_Comm  ,only: Epetra_Comm
  use FEpetra_BlockMap ,only: Epetra_BlockMap
  use iso_c_binding ,only: c_int
  use forepetra
  implicit none
  private                   ! Hide everything by default
  public :: Epetra_Export ! Expose type/constructors/methods

  type ,extends(universal)                 :: Epetra_Export !"shell"
    private
    type(FT_Epetra_Export_ID_t)  :: Export_id 
  contains
     !Developers only
     procedure         :: invalidate_id => invalidate_EpetraExport_ID
     procedure         :: ctrilinos_delete => ctrilinos_delete_EpetraExport
     procedure         :: get_EpetraExport_ID 
     procedure ,nopass :: alias_EpetraExport_ID
     procedure         :: generalize 
     ! Public member functions
     procedure        :: NumSameIDs
     procedure        :: NumPermuteIDs
     procedure        :: NumRemoteIDs
     procedure        :: NumExportIDs
     procedure        :: NumSend
     procedure        :: NumRecv
     procedure        :: SourceMap
     procedure        :: TargetMap
  end type

   interface Epetra_Export ! constructors
     module procedure create,duplicate,from_struct
   end interface
 
contains
  type(Epetra_Export) function from_struct(id)
     type(FT_Epetra_Export_ID_t) ,intent(in) :: id
     from_struct%Export_id = id
     call from_struct%register_self
  end function
 
  ! Original C++ prototype:
  ! Epetra_Export( const Epetra_BlockMap & SourceMap, const Epetra_BlockMap & TargetMap );
  ! CTrilinos prototype:
  ! CT_Epetra_Export_ID_t Epetra_Export_Create ( CT_Epetra_BlockMap_ID_t SourceMapID, CT_Epetra_BlockMap_ID_t TargetMapID );

  type(Epetra_Export) function create(SourceMap,TargetMap)
   !use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_BlockMap_ID_t
    type(Epetra_BlockMap), intent(in) :: TargetMap
    type(Epetra_BlockMap), intent(in) :: SourceMap
    type(FT_Epetra_Export_ID_t) :: create_id
    create_id = Epetra_Export_Create(SourceMap%get_EpetraBlockMap_ID(),TargetMap%get_EpetraBlockMap_ID())
    create = from_struct(create_id)
  end function

  ! Original C++ prototype:
  ! Epetra_Export(const Epetra_Import& Importer);
  ! CTrilinos prototype:
  ! CT_Epetra_Export_ID_t Epetra_Export_Duplicate ( CT_Epetra_Import_ID_t ImporterID );

  type(Epetra_Export) function duplicate(this)
    type(Epetra_Export) ,intent(in) :: this
    type(FT_Epetra_Export_ID_t) :: duplicate_id
    duplicate_id = Epetra_Export_Duplicate(this%Export_id)
    duplicate = from_struct(duplicate_id)
  end function

  type(FT_Epetra_Export_ID_t) function get_EpetraExport_ID(this)
    class(Epetra_Export) ,intent(in) :: this 
    get_EpetraExport_ID=this%Export_id
  end function
  
  type(FT_Epetra_Export_ID_t) function alias_EpetraExport_ID(generic_id)
    use ForTrilinos_table_man
    use ForTrilinos_enums ,only: ForTrilinos_Universal_ID_t,FT_Epetra_Export_ID
    use iso_c_binding     ,only: c_loc,c_int
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,pointer    :: alias_id
    integer(c_int) :: status
    type(error) :: ierr
    if (.not.associated(alias_id)) then
      allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_Export_ID),stat=status)
      ierr=error(status,'FEpetra_Export:alias_EpetraExport_ID')
      call ierr%check_success()
    endif
    alias_EpetraExport_ID=degeneralize_EpetraExport(c_loc(alias_id))
    call deallocate_and_check_error(alias_id,'FEpetra_Export:alias_EpetraExport_ID')
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding     ,only: c_loc
   class(Epetra_Export) ,intent(in) ,target :: this
   generalize = generalize_all(c_loc(this%Export_ID))
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   !class(Epetra_Export) ,intent(in) ,target :: this
   !generalize = Epetra_Export_Generalize ( this%Export_id)
   ! ____ Use for CTrilinos function implementation ______
  end function

 type(FT_Epetra_Export_ID_t) function degeneralize_EpetraExport(generic_id) bind(C)
   ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_Export_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                   ,value   :: generic_id
    type(FT_Epetra_Export_ID_t) ,pointer :: local_ptr
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraExport = local_ptr
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   !type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
   !degeneralize_EpetraExport = Epetra_Export_Degeneralize(generic_id)
   ! ____ Use for CTrilinos function implementation ______
  end function
 
  integer(c_int) function NumSameIDs(this)
    class(Epetra_Export), intent(in) :: this
    NumSameIDs=Epetra_Export_NumSameIDs(this%Export_id)
  end function

  integer(c_int) function NumPermuteIDs(this)
    class(Epetra_Export), intent(in) :: this
    NumPermuteIDs=Epetra_Export_NumPermuteIDs(this%Export_id)
  end function

  integer(c_int) function NumRemoteIDs(this)
    class(Epetra_Export), intent(in) :: this
    NumRemoteIDs=Epetra_Export_NumRemoteIDs(this%Export_id)
  end function

  integer(c_int) function NumExportIDs(this)
    class(Epetra_Export), intent(in) :: this
    NumExportIDs=Epetra_Export_NumExportIDs(this%Export_id)
  end function

  integer(c_int) function NumSend(this)
    class(Epetra_Export), intent(in) :: this
    NumSend=Epetra_Export_NumSend(this%Export_id)
  end function

  integer(c_int) function NumRecv(this)
    class(Epetra_Export), intent(in) :: this
    NumRecv=Epetra_Export_NumRecv(this%Export_id)
  end function

  type(Epetra_BlockMap) function SourceMap(this)
   class(Epetra_Export), intent(in) :: this
   type(FT_Epetra_BlockMap_ID_t) :: SourceMap_id
   SourceMap_id=Epetra_Export_SourceMap(this%Export_id)
   SourceMap=Epetra_BlockMap(SourceMap_id)
  end function

  type(Epetra_BlockMap) function TargetMap(this)
   class(Epetra_Export), intent(in) :: this
   type(FT_Epetra_BlockMap_ID_t) :: TargetMap_id
   TargetMap_id=Epetra_Export_TargetMap(this%Export_id)
   TargetMap=Epetra_BlockMap(TargetMap_id)
  end function

  subroutine invalidate_EpetraExport_ID(this)
    class(Epetra_Export) ,intent(inout) :: this
    this%Export_id%table = FT_Invalid_ID
    this%Export_id%index = FT_Invalid_Index 
    this%Export_id%is_const = FT_FALSE
  end subroutine

  subroutine ctrilinos_delete_EpetraExport(this)
    class(Epetra_Export) ,intent(inout) :: this
    call Epetra_Export_Destroy( this%Export_id ) 
  end subroutine

end module 

