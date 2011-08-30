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


module FEpetra_Import
  use ForTrilinos_enums !,only: FT_Epetra_Comm_ID_t,FT_Epetra_Import_ID_t,FT_Epetra_BlockMap_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_universal ,only:universal
  use ForTrilinos_error
  use FEpetra_Comm  ,only: Epetra_Comm
  use FEpetra_BlockMap ,only: Epetra_BlockMap
  use iso_c_binding ,only: c_int
  use forepetra
  implicit none
  private                   ! Hide everything by default
  public :: Epetra_Import ! Expose type/constructors/methods

  !> <BR> Epetra_Import:This class builds an import object for efficient importing of off-processor elements.

   !> @brief Epetra_Import is used to construct a communication plan that can be called repeatedly by computational classes such the Epetra matrix, vector and multivector classes to efficiently obtain off-processor elements.
   !! This class currently has one constructor, taking two Epetra_Map or Epetra_BlockMap objects. The first map specifies the global IDs of elements that we want to import later. The second map specifies the global IDs that are owned by the calling processor. 

  type Epetra_Import  !,extends(universal)  :: Epetra_Import !"shell"
  contains
     ! Public member functions
     !procedure        :: NumSameIDs
     !procedure        :: NumPermuteIDs
     !procedure        :: PermuteFromLIDs
     !procedure        :: NumRemoteIDs
     !procedure        :: NumExportIDs
     !procedure        :: NumSend
     !procedure        :: NumRecv
     !procedure        :: SourceMap
     !procedure        :: TargetMap
  end type

contains
  ! Original C++ prototype:
  ! Epetra_Import( const Epetra_BlockMap & TargetMap, const Epetra_BlockMap & SourceMap );
  ! CTrilinos prototype:
  ! CT_Epetra_Import_ID_t Epetra_Import_Create ( CT_Epetra_BlockMap_ID_t TargetMapID,
  ! CT_Epetra_BlockMap_ID_t SourceMapID );

  !> @name Constructor Functions
  !! @{

  !> <BR> Epetra_Import Constructor
  !> @brief Constructs a Epetra_Import object from the source and target maps.
  !! This constructor builds an Epetra_Import object by comparing the GID lists of the source and target maps.
  type(Epetra_Import) function Epetra_Import(TargetMap,SourceMap)
   !use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_BlockMap_ID_t
    class(Epetra_BlockMap), intent(in) :: TargetMap &
    !< In Map containing the GIDs from which data should be imported to each processor from the source map whenever an import operation is performed using this importer.
    class(Epetra_BlockMap), intent(in) :: SourceMap &
    !< In Map containing the GIDs that should be used for importing data. 
    !< <ol>
    !> @details Warning: Note that the SourceMap \e must have GIDs uniquely owned, each GID of the source map can occur only once.
    !< <ol>
    !< Builds an import object that will transfer objects built with SourceMap to objects built with TargetMap.
    !< A Epetra_Import object categorizes the elements of the target map into three sets as follows:
    !< <ol>
    !< <li> All elements in the target map that have the same GID as the corresponding element of the source map, starting with the first element in the target map, going up to the first element that is different from the source map. The number of these IDs is returned by NumSameIDs().
    <li> All elements that are local to the processor, but are not part of the first set of elements.  These elements         have GIDs that are owned by the calling processor, but at least the first element of this list is permuted.
         Even if subsequent elements are not permuted, they are included in this list.  The number of permuted elements
         is returned by NumPermutedIDs().  The list of elements (local IDs) in the source map that are permuted can be
found in the list PermuteFromLIDs().  The list of elements (local IDs) in the target map that are the new locations
         of the source elements can be found in the list PermuteToLIDs().
    <li> All remaining elements of the target map correspond to global IDs that are owned by remote processors.  The number 
         of these elements is returned by NumRemoteIDs() and the list of these is returned by RemoteLIDs().
    </ol>

Given the above information, the Epetra_Import constructor builds a list of elements that must be communicated to other
processors as a result of import requests.  The number of exported elements (where multiple sends of the same element
to different processors is counted) is returned by NumExportIDs().  The local IDs to be sent are returned by the list 
ExportLIDs().  The processors to which each of the elements will be sent in returned in a list of the same length by 
ExportPIDs().

The total number of elements that will be sent by the calling processor is returned by NumSend().  The total number of
elements that will be received is returned by NumRecv().


The following example illustrates the basic concepts.

Assume we have 3 processors and 9 global elements with each processor owning 3 elements as follows
\verbatim
 PE 0 Elements |  PE 1 Elements  |  PE 2 Elements
    0  1  2          3  4  5           6  7  8
\endverbatim

The above layout essentially defines the source map argument of the import object.

This could correspond to a 9 by 9 matrix with the first three rows on PE 0, and so on.  Suppose that this matrix
is periodic tridiagonal having the following sparsity pattern:
\verbatim

PE 0 Rows:

  X  X  0  0  0  0  0  0  X
  X  X  X  0  0  0  0  0  0
  0  X  X  X  0  0  0  0  0

PE 1 Rows:

  0  0  X  X  X  0  0  0  0
  0  0  0  X  X  X  0  0  0
  0  0  0  0  X  X  X  0  0

PE 2 Rows:

  0  0  0  0  0  X  X  X  0
  0  0  0  0  0  0  X  X  X
  X  0  0  0  0  0  0  X  X

\endverbatim

To perform a matrix vector multiplication operation y = A*x (assuming that x has the same distribution as the 
rows of the matrix A) each processor will need to import elements of x that
are not local.  To do this, we build a target map on each processor as follows:
\verbatim
    PE 0 Elements    |  PE 1 Elements    |  PE 2 Elements
    0  1  2  3  8       2  3  4  5  6       0  5  6  7  8
\endverbatim
The above list is the elements that will be needed to perform the matrix vector multiplication locally on each processor.
Note that the ordering of the elements on each processor is not unique, but has been chosen for illustration.

With these two maps passed into the Epetra_Import constructor, we get the following attribute definitions:

On PE 0:

\verbatim
NumSameIDs      = 3

NumPermuteIDs   = 0
PermuteToLIDs   = 0
PermuteFromLIDs = 0
NumRemoteIDs    = 2
RemoteLIDs      = [3, 4]

NumExportIDs    = 2
ExportLIDs      = [0, 2]
ExportPIDs      = [1, 2]

NumSend         = 2
NumRecv         = 2

\endverbatim

On PE 1:

\verbatim
NumSameIDs      = 0

NumPermuteIDs   = 3
PermuteFromLIDs = [0, 1, 2]
PermuteToLIDs   = [1, 2, 3]

NumRemoteIDs    = 2
RemoteLIDs      = [0, 4]

NumExportIDs    = 2
ExportLIDs      = [0, 2]
ExportPIDs      = [0, 2]

NumSend         = 2
NumRecv         = 2

\endverbatim

On PE 2:

\verbatim
NumSameIDs      = 0

NumPermuteIDs   = 3
PermuteFromLIDs = [0, 1, 2]
PermuteToLIDs   = [2, 3, 4]

NumRemoteIDs    = 2
RemoteLIDs      = [0, 1]

NumExportIDs    = 2
ExportLIDs      = [0, 2]
ExportPIDs      = [0, 1]

NumSend         = 2
NumRecv         = 2

\endverbatim

<b> Using Epetra_Import Objects </b>
Once a Epetra_Import object has been constructed, it can be used by any of the Epetra classes that support distributed global
objects, namely Epetra_Vector, Epetra_MultiVector, Epetra_CrsGraph, Epetra_CrsMatrix and Epetra_VbrMatrix.  
All of these classes have Import and Export methods that will fill new objects whose distribution is described by
the target map, taking elements from the source object whose distribution is described by the source map.  Details of usage
for each class is given in the appropriate class documentation.

Note that the reverse operation, an export, using this importer is also possible and appropriate in some instances.
For example, if we compute y = A^Tx, the transpose matrix-multiplication operation, then we can use the importer we constructed
in the above example to do an export operation to y, adding the contributions that come from multiple processors.

  end function

  ! Original C++ prototype:
  ! Epetra_Import(const Epetra_Import& Importer);
  ! CTrilinos prototype:
  ! CT_Epetra_Import_ID_t Epetra_Import_Duplicate ( CT_Epetra_Import_ID_t ImporterID );

  !> @name Constructor Functions
  !! @{

  !> <BR> Epetra_Import Constructor
  !> @brief 
  type(Epetra_Import) function duplicate(this)
    type(Epetra_Import) ,intent(in) :: this 
    type(FT_Epetra_Import_ID_t) :: duplicate_id
    duplicate_id = Epetra_Import_Duplicate(this%Import_id)
    duplicate = from_struct(duplicate_id)
  end function

  type(FT_Epetra_Import_ID_t) function get_EpetraImport_ID(this)
    class(Epetra_Import) ,intent(in) :: this 
    get_EpetraImport_ID=this%Import_id
  end function
  
  type(FT_Epetra_Import_ID_t) function alias_EpetraImport_ID(generic_id)
    use ForTrilinos_table_man
    use ForTrilinos_enums ,only: ForTrilinos_Universal_ID_t,FT_Epetra_Import_ID
    use iso_c_binding     ,only: c_loc,c_int
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,allocatable ,target :: alias_id
    integer(c_int) :: status
    type(error) :: ierr
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_Import_ID),stat=status)
    ierr=error(status,'FEpetra_Import:alias_EpetraImport_ID')
    call ierr%check_success()
    alias_EpetraImport_ID=degeneralize_EpetraImport(c_loc(alias_id))
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding     ,only: c_loc
   class(Epetra_Import) ,intent(in) ,target :: this
   generalize = generalize_all(c_loc(this%Import_ID))
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   !class(Epetra_Import) ,intent(in) ,target :: this
   !generalize = Epetra_Import_Generalize ( this%Import_id)
   ! ____ Use for CTrilinos function implementation ______
  end function

 type(FT_Epetra_Import_ID_t) function degeneralize_EpetraImport(generic_id) bind(C)
   ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_Import_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                   ,value   :: generic_id
    type(FT_Epetra_Import_ID_t) ,pointer :: local_ptr=>null()
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraImport = local_ptr
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   !type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
   !degeneralize_EpetraImport = Epetra_Import_Degeneralize(generic_id)
   ! ____ Use for CTrilinos function implementation ______
  end function
 
  integer(c_int) function NumSameIDs(this)
    class(Epetra_Import), intent(in) :: this
    NumSameIDs=Epetra_Import_NumSameIDs(this%Import_id)
  end function

  integer(c_int) function NumPermuteIDs(this)
    class(Epetra_Import), intent(in) :: this
    NumPermuteIDs=Epetra_Import_NumPermuteIDs(this%Import_id)
  end function

 function PermuteFromLIDs(this)
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer,c_int
    class(Epetra_Import), intent(in) :: this
    integer(c_int),dimension(:),allocatable :: PermuteFromLIDs
    type(c_ptr)   :: PermuteFromLIDs_external_ptr 
    integer(c_int),pointer :: PermuteFromLIDs_local_ptr=>null()
    integer(c_int) :: status
    type(error) :: ierr
    allocate(PermuteFromLIDs(this%NumPermuteIDs()),stat=status)
    ierr=error(status,'FEpetra_Import:alias_EpetraImport_ID')
    call ierr%check_success()
    PermuteFromLIDs_external_ptr=Epetra_Import_PermuteFromLIDs(this%Import_id)
    call c_f_pointer (PermuteFromLIDs_external_ptr, PermuteFromLIDs_local_ptr)
    PermuteFromLIDs=PermuteFromLIDs_local_ptr
  end function

  integer(c_int) function NumRemoteIDs(this)
    class(Epetra_Import), intent(in) :: this
    NumRemoteIDs=Epetra_Import_NumRemoteIDs(this%Import_id)
  end function

  integer(c_int) function NumExportIDs(this)
    class(Epetra_Import), intent(in) :: this
    NumExportIDs=Epetra_Import_NumExportIDs(this%Import_id)
  end function

  integer(c_int) function NumSend(this)
    class(Epetra_Import), intent(in) :: this
    NumSend=Epetra_Import_NumSend(this%Import_id)
  end function

  integer(c_int) function NumRecv(this)
    class(Epetra_Import), intent(in) :: this
    NumRecv=Epetra_Import_NumRecv(this%Import_id)
  end function

  type(Epetra_BlockMap) function SourceMap(this)
   class(Epetra_Import), intent(in) :: this
   type(FT_Epetra_BlockMap_ID_t) :: SourceMap_id
   SourceMap_id=Epetra_Import_SourceMap(this%Import_id)
   SourceMap=Epetra_BlockMap(SourceMap_id)
  end function

  type(Epetra_BlockMap) function TargetMap(this)
   class(Epetra_Import), intent(in) :: this
   type(FT_Epetra_BlockMap_ID_t) :: TargetMap_id
   TargetMap_id=Epetra_Import_TargetMap(this%Import_id)
   TargetMap=Epetra_BlockMap(TargetMap_id)
  end function

  subroutine invalidate_EpetraImport_ID(this)
    class(Epetra_Import),intent(inout) :: this
    this%Import_id%table = FT_Invalid_ID
    this%Import_id%index = FT_Invalid_Index 
    this%Import_id%is_const = FT_FALSE
  end subroutine

  subroutine ctrilinos_delete_EpetraImport(this)
    class(Epetra_Import),intent(inout) :: this
    call Epetra_Import_Destroy( this%Import_id ) 
  end subroutine
end module 
