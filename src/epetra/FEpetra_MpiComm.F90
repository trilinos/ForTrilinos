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

module FEpetra_MpiComm
#include "ForTrilinos_config.h"
#ifdef HAVE_MPI
!#include "mpif.h"
  use ForTrilinos_enums ,only: FT_Epetra_MpiComm_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_external_utils
  use ForTrilinos_hermetic, only : hermetic
  use FEpetra_Comm      ,only: Epetra_Comm
  use iso_c_binding     ,only: c_int,c_double,c_long,c_char
  use forepetra
  implicit none
  private               ! Hide everything by default
  public :: Epetra_MpiComm ! Expose type/methods

  type ,extends(Epetra_Comm) :: Epetra_MpiComm
    private
    type(FT_Epetra_MpiComm_ID_t) :: MpiComm_id  
  contains
    !Constructor
    !procedure         :: clone
    !Developers only
    procedure         ::remote_dealloc
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
   from_struct%MpiComm_id = id
   call from_struct%set_EpetraComm_ID(from_struct%alias_EpetraComm_ID(from_struct%generalize()))
   call from_struct%register_self
  end function
 
  ! Original C++ prototype:
  ! Epetra_MpiComm(MPI_Comm comm);
  ! CTrilinos prototype:
  ! CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Create ( MPI_Comm comm );

  type(Epetra_MpiComm) function from_scratch(comm)
   integer(c_int) ,intent(in) :: comm
   type(FT_Epetra_MpiComm_ID_t) :: from_scratch_id
   from_scratch_id = Epetra_MpiComm_Fortran_Create(comm)
   from_scratch=from_struct(from_scratch_id)
  end function

  ! Original C++ prototype:
  ! Epetra_MpiComm(const Epetra_MpiComm & Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Duplicate ( CT_Epetra_MpiComm_ID_t CommID );

  type(Epetra_MpiComm) function duplicate(this)
   type(Epetra_MpiComm) ,intent(in) :: this
    type(FT_Epetra_MpiComm_ID_t) :: duplicate_id
    duplicate_id = Epetra_MpiComm_Duplicate(this%MpiComm_id)
    duplicate = from_struct(duplicate_id)
  end function

  !function clone(this)
  !  class(Epetra_MpiComm)    ,intent(in)  :: this
  !  class(Epetra_Comm)       ,allocatable :: clone
  !  type(Epetra_MpiComm)      :: clone_local
  !  type(FT_Epetra_MpiComm_ID_t) :: clone_mpi_id
  !  type(FT_Epetra_Comm_ID_t) :: clone_comm_id
  !  clone_comm_id=Epetra_MpiComm_Clone(this%MpiComm_id)
  !  call clone_local%set_EpetraComm_ID(clone_comm_id)
  !  allocate(Epetra_MpiComm :: clone)
  !  clone=Epetra_MpiComm(alias_EpetraMpiComm_ID(clone_local%generalize_EpetraComm()))
  !  call clone_local%force_finalize()
  !end function

 type(FT_Epetra_MpiComm_ID_t) function get_EpetraMpiComm_ID(this)
   class(Epetra_MpiComm) ,intent(in) :: this
   get_EpetraMpiComm_ID=this%MpiComm_id
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

  subroutine remote_dealloc(this)
    class(Epetra_MpiComm) ,intent(inout) :: this
    call this%remote_dealloc_EpetraComm()
    call Epetra_MpiComm_Destroy(this%MpiComm_id)
  end subroutine
#endif
end module 
