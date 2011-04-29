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
  use ForTrilinos_enums ,only: FT_Epetra_Comm_ID_t,FT_Epetra_BlockMap_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_universal ,only : universal
  use ForTrilinos_error ,only : error
  use FEpetra_Comm  ,only : Epetra_Comm
  use iso_c_binding ,only : c_int
  use forepetra
  implicit none
  private                   ! Hide everything by default
  public :: Epetra_BlockMap ! Expose type/constructors/methods

  type :: Epetra_BlockMap !,extends(universal)      :: Epetra_BlockMap 
  contains
     !Local/Global ID accessor methods
     !Size and dimension acccessor functions
     procedure         :: NumGlobalElements
     procedure         :: NumMyElements
     procedure         :: IndexBase
     procedure         :: SameAs
     procedure         :: PointSameAs
     procedure         :: MyGlobalElements
     procedure         :: ElementSize_Const
     procedure         :: ElementSize_LID
     generic :: ElementSize=>ElementSize_Const,ElementSize_LID
     !Miscellaneous boolean tests
     procedure         :: LinearMap
     procedure         :: DistributedGlobal
     !Array accessor functions
     !Miscellaneous
     procedure         :: Comm
  end type

 
contains

  !> @name Constructor Functions
  !! @{
 
  !< <BR> Epetra_BlockMap constructor for a Epetra-defined uniform linear distribution of constant size elements.
  type(Epetra_BlockMap) function Epetra_BlockMap(Num_GlobalElements,Element_Size,IndexBase,comm)
    integer(c_int) ,intent(in) :: Num_GlobalElements
    integer(c_int) ,intent(in) :: Element_Size
    integer(c_int) ,intent(in) :: IndexBase
    class(Epetra_Comm)         :: comm
  end function

  !> @name Constructor Functions
  !! @{
  type(Epetra_BlockMap) function Epetra_BlockMap(Num_GlobalElements,Num_MyElements,Element_Size,IndexBase,comm)
    integer(c_int) ,intent(in) :: Num_GlobalElements
    integer(c_int) ,intent(in) :: Num_MyElements
    integer(c_int) ,intent(in) :: Element_Size
    integer(c_int) ,intent(in) :: IndexBase
    class(Epetra_Comm)         :: comm
  end function

  !> @name Constructor Functions
  !! @{
  type(Epetra_BlockMap) function Epetra_BlockMap(Num_GlobalElements,Num_MyElements,&
                                                        My_GlobalElements,Element_Size,IndexBase,comm)
    integer(c_int) ,intent(in) :: Num_GlobalElements
    integer(c_int) ,intent(in) :: Num_MyElements
    integer(c_int) ,intent(in) ,dimension(:) :: My_GlobalElements  
    integer(c_int) ,intent(in) :: Element_Size
    integer(c_int) ,intent(in) :: IndexBase
    class(Epetra_Comm)         :: comm
  end function


  !> @name Constructor Functions
  !! @{
  type(Epetra_BlockMap) function Epetra_BlockMap(Num_GlobalElements,Num_MyElements,&
                                                       My_GlobalElements,Element_SizeList,IndexBase,comm)
    integer(c_int) ,intent(in) :: Num_GlobalElements
    integer(c_int) ,intent(in) :: Num_MyElements
    integer(c_int) ,intent(in) ,dimension(:) :: My_GlobalElements  
    integer(c_int) ,intent(in) ,dimension(:) :: Element_SizeList    
    integer(c_int) ,intent(in) :: IndexBase
    class(Epetra_Comm)         :: comm
  end function

  !> @name Constructor Functions
  !! @{
  type(Epetra_BlockMap) function Epetra_BlockMap(this)
    type(Epetra_BlockMap) ,intent(in) :: this 
  end function

  integer(c_int) function NumGlobalElements(this)
    class(Epetra_BlockMap) ,intent(in) :: this
  end function 

  integer(c_int) function NumMyElements(this)
    class(Epetra_BlockMap) ,intent(in) :: this
  end function 

  integer(c_int) function IndexBase(this)
    class(Epetra_BlockMap) ,intent(in) :: this
  end function 

  logical function  SameAs(lhs,rhs)
    use ForTrilinos_enums, only:FT_boolean_t,FT_FALSE,FT_TRUE
    class(Epetra_BlockMap)        ,intent(in) :: lhs
    class(Epetra_BlockMap)        ,intent(in) :: rhs
  end function SameAs

  logical function  PointSameAs(lhs,rhs)
    use ForTrilinos_enums, only:FT_boolean_t,FT_FALSE,FT_TRUE
    class(Epetra_BlockMap)        ,intent(in) :: lhs
    class(Epetra_BlockMap)        ,intent(in) :: rhs
  end function PointSameAs
 
  function MyGlobalElements(this) result(MyGlobalElementsList)
    class(Epetra_BlockMap)     ,intent(in)    :: this
    integer(c_int),dimension(:),allocatable   :: MyGlobalElementsList
    integer(c_int)                            :: junk
  end function 

  integer(c_int) function ElementSize_Const(this)
    class(Epetra_BlockMap) ,intent(in) :: this
  end function 

  integer(c_int) function ElementSize_LID(this,L_ID)
    class(Epetra_BlockMap) ,intent(in) :: this
    integer(c_int)         ,intent(in) :: L_ID
  end function 

  logical function LinearMap(this)
    use ForTrilinos_enums, only:FT_boolean_t,FT_FALSE,FT_TRUE
    class(Epetra_BlockMap),intent(in) :: this
  end function

  logical function DistributedGlobal(this)
    use ForTrilinos_enums, only:FT_boolean_t,FT_FALSE,FT_TRUE
    class(Epetra_BlockMap),intent(in) :: this
  end function

  subroutine Comm(this,communicator)
    use FEpetra_MpiComm
    use FEpetra_SerialComm
    class(Epetra_BlockMap), intent(in) :: this
    class(Epetra_Comm), allocatable :: communicator
 end subroutine

end module 

