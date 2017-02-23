!-----------------------------------------------------------------------------!
! \file   parameterlist/test.f90
! \author Seth R Johnson
! \date   Tue Dec 06 18:15:21 2016
! \brief  test module
! \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
!-----------------------------------------------------------------------------!

program main

#include "FortranTestMacros.h"

    use ISO_FORTRAN_ENV
    implicit none

    call test_plist()
contains

subroutine test_plist()
    use ISO_FORTRAN_ENV
    use, intrinsic :: ISO_C_BINDING
    use forteuchos
    implicit none

    type(ParameterList) :: plist, sublist
    integer :: val
    logical :: bool_result
    integer, dimension(6) :: test_int = (/ -1, 1, 3, 3, 5, 7 /)
    real(C_DOUBLE), dimension(4) :: test_dbl = (/ 0.1d0, 1.9d0, -2.0d0, 4.0d0 /)

    write(0, *) "Constructing..."
    call plist%create("myname")
    ! call load_from_xml(plist, "input_params.xml")

    ! Get and set a vlaue
    call plist%set("myint", 4)
    call plist%get("myint", val)
    EXPECT_EQ(4, val)

    call plist%set("mydouble", 4.0d0)
    call plist%set("intarr", test_int)
    call plist%set("dblarr", test_dbl)
    
    val = plist%get_length("intarr")
    EXPECT_EQ(6, val)

    call plist%set("deleteme", 123)
    bool_result = plist%is_parameter("deleteme")
    EXPECT_TRUE(bool_result)
    call plist%remove("deleteme")
    bool_result = plist%is_parameter("deleteme")
    EXPECT_FALSE(bool_result)

    call sublist%create("sublist")
    call sublist%set("anotherval", 4.0d0)
    call sublist%set("stringval", "some string!")
    val = sublist%get_length("stringval")
    EXPECT_EQ(12, val)

    call plist%set("sublist", sublist)
    call sublist%release()

    write(0, *) "Printing..."
    call plist%print()
    write(0, *) "Saving to XML file"
    call save_to_xml(plist, "myparams.xml")

    call plist%release()
end subroutine

end program

!-----------------------------------------------------------------------------!
! end of parameterlist/test.f90
!-----------------------------------------------------------------------------!
