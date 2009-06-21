module ForTrilinos_enums
  use iso_c_binding ,only : c_int ! Kind parameter (precision specifier)
  implicit none                   ! Prevent implicit typing

  ! This file is the Fortran companion to CTrilinos/src/CTrilinos_enums.h.
  ! The definitions here are the nearest Fortran equivalents to the contents of that file.
  ! The Fortran 2003 standard guarantees that the types (e.g., 'integer(c_int)') and 
  ! enumeration values are guaranteed by the Fortran 2003 standard to be interoperable with 
  ! those employed by a companion C compiler.

  enum ,bind(C) 
    enumerator ::
      FT_Invalid_ID, 
      FT_Epetra_Distributor_ID, 
      FT_Epetra_SerialCommData_ID, 
      FT_Epetra_SerialComm_ID, 
      FT_Epetra_BLAS_ID, 
      FT_Epetra_Comm_ID, 
      FT_Epetra_Operator_ID, 
      FT_Epetra_MultiVector_ID, 
      FT_Epetra_OffsetIndex_ID, 
      FT_Epetra_Object_ID, 
      FT_Epetra_Data_ID, 
      FT_Epetra_RowMatrix_ID, 
      FT_Epetra_CompObject_ID, 
      FT_Epetra_Directory_ID, 
      FT_Epetra_BlockMapData_ID, 
      FT_Epetra_Flops_ID, 
      FT_Epetra_SrcDistObject_ID, 
      FT_Epetra_MpiComm_ID, 
      FT_Epetra_MpiCommData_ID, 
      FT_Epetra_CrsMatrix_ID, 
      FT_Epetra_CrsGraph_ID, 
      FT_Epetra_DistObject_ID, 
      FT_Epetra_Export_ID, 
      FT_Epetra_Vector_ID, 
      FT_Epetra_Map_ID, 
      FT_Epetra_CrsGraphData_ID, 
      FT_Epetra_BlockMap_ID, 
      FT_Epetra_Import_ID
  end enum 

  type ,bind(C) :: ForTrilinos_Type_ID_t  
    integer(c_int) :: id_t                 ! Companion to CTrilinos_Type_ID_t  
  end type
  
  type ,bind(C) :: ForTrilinos_Object_ID_t ! interoperable with CTrilinos_Object_ID_t
   !type(ForTrilinos_Type_ID_t)  :: type   ! Could use this form if CTrilinos_Type_ID_t were a struct
    integer(c_int)  :: type                ! Data type of the object (interoperable with CTrilinos_Type_ID_t)
    integer(c_int)  :: index               ! Array index of the object 
  end type 

  ! To keep the function interfaces readable... 
  ! (each type_ID_t component interoperates with a CT_Epetra_*_ID_t type, *=Distributor,SerialCommData,...)
  ! (each 'type FT_*_ID_t type would interoperate directly with CT_Epetra_*_ID_t if the latter were a struct)
  type FT_Epetra_Distributor_ID_t   ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_SerialCommData_ID_t; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_SerialComm_ID_t    ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_BLAS_ID_t          ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_Comm_ID_t          ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_Operator_ID_t      ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_MultiVector_ID_t   ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_OffsetIndex_ID_t   ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_Object_ID_t        ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_Data_ID_t          ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_RowMatrix_ID_t     ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_CompObject_ID_t    ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_Directory_ID_t     ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_BlockMapData_ID_t  ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_Flops_ID_t         ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_SrcDistObject_ID_t ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_MpiComm_ID_t       ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_MpiCommData_ID_t   ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_CrsMatrix_ID_t     ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_CrsGraph_ID_t      ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_DistObject_ID_t    ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_Export_ID_t        ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_Vector_ID_t        ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_Map_ID_t           ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_CrsGraphData_ID_t  ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_BlockMap_ID_t      ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  type FT_Epetra_Import_ID_t        ; type(ForTrilinos_Object_ID_t) :: type_ID_t; end type
  
end module ForTrilinos_enums
