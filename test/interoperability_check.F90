module interoperability_check

  ! The valid_kind_parameters() procedure in this module checks values of parameters required to interoperate with CTrilinos. 
  ! The Fortran 2003 standard requires that these parameters be defined in the intrinsic module iso_c_binding with values
  ! that communicate meanings specified in the standard.  Much of the quoted text in this module is excerpted from the standard.
 
  use ,intrinsic :: iso_c_binding ,only : &
    c_int,c_char,c_double,c_ptr,c_long,c_bool,c_null_char
  use ,intrinsic :: iso_fortran_env ,only : error_unit ,output_unit

  implicit none
  private
  public :: valid_kind_parameters

  ! For purposes of ForTrilinos, the Fortran standard's use of the word 'processor' is interpreted as denoting the 
  ! combination of a compiler, an operating system, and a hardware architecture.

  character(len=*) ,parameter :: no_fortran_kind= &
    'The companion C processor defines the corresponding C type, &
    & but there is no interoperating Fortran processor kind.'
  character(len=*) ,parameter :: no_c_type= &
    'The C processor does not define the corresponding C type.'
  character(len=*) ,parameter :: interoperable = &
    'An interoperating Fortran kind has a corresponding C type.'
  character(len=*) ,parameter :: imprecise= &
    'The C processor type does not have a precision equal to the precision of any of the Fortran processor real kinds.'
  character(len=*) ,parameter :: limited= &
    'The C processor type does not have a range equal to the range of any of the Fortran processor real kinds.'
  character(len=*) ,parameter :: limited_and_imprecise= &
    'The C processor type has neither the precision nor range of any of the Fortran processor real kinds.'
  character(len=*) ,parameter :: not_interoperable_nonspecific = &
    'There is no interoperating Fortran processor kind for unspecified reasons.'

contains

  ! Return true if all of the interoperating Fortran kind parameters required by ForTrilinos have a corresponding C type
  ! defined by the companion C processor.  Otherwise, return false.

  logical function valid_kind_parameters()

    valid_kind_parameters = .true.

    select case(c_long)
      case(-1)
        write(error_unit ,fmt='(2a)') 'c_long error: ',no_fortran_kind
        valid_kind_parameters = .false.
      case(-2)
        write(error_unit ,fmt='(2a)') 'c_long error: ',no_c_type
        valid_kind_parameters = .false.
      case default
        write(output_unit,fmt='(2a)') 'c_long: ',interoperable
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
        write(output_unit,fmt='(2a)') 'c_double: ',interoperable
    end select
  
    select case(c_bool)
      case(-1)
        write(error_unit ,fmt='(a)') 'c_bool error: invalid value for a logical kind parameter on the processor.'
        valid_kind_parameters = .false.
      case default
        write(output_unit ,fmt='(a)') 'c_bool:  valid value for a logical kind parameter on the processor.'
    end select
  
    select case(c_char)
      case(-1)
        write(error_unit ,fmt='(a)') 'c_char error: invalid value for a character kind parameter on the processor.'
        valid_kind_parameters = .false.
      case default
        write(output_unit ,fmt='(a)') 'c_char:  valid value for a character kind type parameter on the processor.'
    end select
  
  end function

end module interoperability_check
