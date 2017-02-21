!-----------------------------------------------------------------------------!
! \file   parameterlist/test.f90
! \author Seth R Johnson
! \date   Tue Dec 06 18:15:21 2016
! \brief  test module
! \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
!-----------------------------------------------------------------------------!

program main
    use ISO_FORTRAN_ENV

    call test_plist()

contains
subroutine test_plist()
    use ISO_FORTRAN_ENV
    use, intrinsic :: ISO_C_BINDING
    use teuchos
    implicit none
    type(ParameterList) :: plist, sublist
    integer :: val
    integer, dimension(6) :: test_int = (/ -1, 1, 3, 3, 5, 7 /)
    real(C_DOUBLE), dimension(4) :: test_dbl = (/ 0.1d0, 1.9d0, -2.0d0, 4.0d0 /)

    write(0, *) "Constructing..."
    call plist%create("myname")
    ! call load_from_xml(plist, "input_params.xml")

    call plist%set("myint", 4)
    call plist%set("mydouble", 4.0d0)
    call plist%set("intarr", test_int)
    call plist%set("dblarr", test_dbl)
    write(0, *) "Array length: ", plist%get_length("intarr")
    call plist%set("deleteme", 123)
    write(0, *) "Should be a parameter: ", plist%is_parameter("deleteme")
    call plist%remove("deleteme")
    write(0, *) "Should be false: ", plist%is_parameter("deleteme")

    call plist%get("myint", val)
    write(0, *) "Retrieved ", val

    call sublist%create("sublist")
    call sublist%set("anotherval", 4.0d0)
    call sublist%set("stringval", "some string!")

    call plist%set("sublist", sublist)
    call sublist%release()

    write(0, *) "Printing"
    call plist%print()
    write(0, *) "Saving to XML file"
    call save_to_xml(plist, "myparams.xml")

    call plist%release()
end subroutine

end program

!-----------------------------------------------------------------------------!
! end of parameterlist/test.f90
!-----------------------------------------------------------------------------!
