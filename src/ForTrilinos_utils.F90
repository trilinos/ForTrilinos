module fortrilinos_utils
  implicit none
  private
  public :: count
  public :: generalize_all
  public :: valid_kind_parameters
contains
  
  ! Fortran-style strings take two forms: one with no C counterpart and a second that is
  ! an array of character values.  A C string is interoperable with the second type of 
  ! Fortran string with the primary distinction being that Fortran strings carry their
  ! own string-length information, whereas the length of C strings is indicated by a 
  ! terminating null character. The "count" procedure calculates C string length by 
  ! searching for the null terminator.

  integer(c_int) function count(ptr) bind(C,name="forLinkingOnly") 
    use ,intrinsic :: iso_c_binding ,only: c_char ,c_int  ,c_null_char
    implicit none 
    character(kind=c_char) :: ptr(*) 
    count = 0 
    do 
       if(ptr(count+1) == c_null_char) return 
       count = count+1 
    end do 
  end function 

  ! This is a Fortran implementation of the functionality in the Epetra_*_Abstract procedures
  ! in each CTrilinos/src/CEpetra* file.  It effectively casts any given Epetra derived type
  ! to a general type that can represent any of the Epetra derived types in 
  ! ForTrilinos/src/ForTrilinos_enums.F90.

  type(ForTrilinos_Object_ID_t) function generalize_all(object_id) bind(C,name="for_linking_only") 
    use ForTrilinos_enums ,only : ForTrilinos_Object_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr) ,value :: object_id
    type(ForTrilinos_Object_ID_t), pointer :: local_ptr

    call c_f_pointer (object_id, local_ptr)
    generalize_all = local_ptr
  end function

  ! This procedure checks the values of parameters required to interoperate with CTrilinos.  The Fortran 2003 
  ! standard requires that these parameters be defined in the intrinsic module iso_c_binding with values
  ! that communicate meanings specified in the standard.  This procedure's quoted interpretations of these 
  ! values are largely excerpted from the standard.  This procedure returns true if all of the interoperating
  ! Fortran kind  parameters required by ForTrilinos have a corresponding C type defined by the companion 
  ! C processor.  Otherwise, return false.  (For purposes of ForTrilinos, the Fortran standard's use of the
  ! word 'processor' is interpreted as denoting the combination of a compiler, an operating system, and a
  ! hardware architecture.)

  logical function valid_kind_parameters(verbose)
    use ,intrinsic :: iso_c_binding ,only : &
      c_int,c_char,c_double,c_ptr,c_long,c_bool,c_null_char
    use ,intrinsic :: iso_fortran_env ,only : error_unit ,output_unit
    logical ,optional :: verbose
    logical           :: verbose_output


    character(len=*) ,parameter :: no_fortran_kind= &
      'The companion C processor defines the corresponding C type, &
      & but there is no interoperating Fortran processor kind.'
    character(len=*) ,parameter :: no_c_type= &
      'The C processor does not define the corresponding C type.'
    character(len=*) ,parameter :: interoperable = &
      'This interoperating Fortran kind has a corresponding C type.'
    character(len=*) ,parameter :: imprecise= &
      'The C processor type does not have a precision equal to the precision of any of the Fortran processor real kinds.'
    character(len=*) ,parameter :: limited= &
      'The C processor type does not have a range equal to the range of any of the Fortran processor real kinds.'
    character(len=*) ,parameter :: limited_and_imprecise= &
      'The C processor type has neither the precision nor range of any of the Fortran processor real kinds.'
    character(len=*) ,parameter :: not_interoperable_nonspecific = &
      'There is no interoperating Fortran processor kind for unspecified reasons.'

    valid_kind_parameters = .true. ! default return value

    if (present(verbose)) then 
      verbose_output=verbose
    else
      verbose_output=.false.
    end if

    select case(c_long)
      case(-1)
        write(error_unit ,fmt='(2a)') 'c_long error: ',no_fortran_kind
        valid_kind_parameters = .false.
      case(-2)
        write(error_unit ,fmt='(2a)') 'c_long error: ',no_c_type
        valid_kind_parameters = .false.
      case default
        if (verbose_output) write(output_unit,fmt='(2a)') 'c_long: ',interoperable
    end select
  
    select case(c_double)
      case(-1)
        write(error_unit ,fmt='(2a)') 'c_double error: ',imprecise
        valid_kind_parameters = .false.
      case(-2)
        write(error_unit ,fmt='(2a)') 'c_double error: ',limited
        valid_kind_parameters = .false.
      case(-3)
        write(error_unit ,fmt='(2a)') 'c_double error: ',limited_and_imprecise
        valid_kind_parameters = .false.
      case(-4)
        write(error_unit ,fmt='(2a)') 'c_double error: ',not_interoperable_nonspecific
        valid_kind_parameters = .false.
      case default
        if (verbose_output) write(output_unit,fmt='(2a)') 'c_double: ',interoperable
    end select
  
    select case(c_bool)
      case(-1)
        write(error_unit ,fmt='(a)') 'c_bool error: invalid value for a logical kind parameter on the processor.'
        valid_kind_parameters = .false.
      case default
        if (verbose_output) write(output_unit ,fmt='(a)') 'c_bool:  valid value for a logical kind parameter on the processor.'
    end select
  
    select case(c_char)
      case(-1)
        write(error_unit ,fmt='(a)') 'c_char error: invalid value for a character kind parameter on the processor.'
        valid_kind_parameters = .false.
      case default
        if (verbose_output) write(output_unit ,fmt='(a)') 'c_char:  valid value for a character kind type parameter on the processor.'
    end select
  end function

end module fortrilinos_utils
