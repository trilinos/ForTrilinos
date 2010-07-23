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
! Questions? Contact Karla Morris  (knmorri@sandia.gov)
!                    Damian Rouson (rouson@sandia.gov)
!*********************************************************************



program main

  ! This file is the object-oriented equivalent of verySimple.F90.  In Trilinos 10.4,
  ! this is a snapshot of an unstable (evolving) file expected to become stable in a
  ! subsequent release.  This file exercises the derived types defined in 
  ! ForTrilinos/src/epetra/FEpetra*.F90, which wrap the interface bodies in 
  ! ForTrilinos/src/epetra/forepetra.F90.   
    
  ! This file represents the preferred style for using ForTrilinos and is recommended for 
  ! Fortran users whose compilers support the object-oriented features of Fortran 2003.
  ! As of the Trilinos 10.4 release date, the latest versions of the IBM and Cray compilers 
  ! nominally support the required features.  The Numerical Algorithms Group (NAG) and Intel 
  ! compilers support all features but one: final subroutines.  (In each case, the support is
  ! somewhat immature and buggy.)  

  use iso_c_binding,only:c_int,c_double
  use FEpetra_SerialComm   ,only : Epetra_SerialComm
  use FEpetra_Map          ,only : Epetra_Map
  use FEpetra_Vector       ,only : Epetra_Vector
  use ForTrilinos_utils    ,only : valid_kind_parameters
  use ForTrilinos_error    !,only : error
  implicit none

  ! Data declarations 
  
  type(Epetra_SerialComm) :: communicator
  type(Epetra_Map)    :: map
  type(Epetra_Vector) :: x, b
  type(error)         :: err
  integer(c_int) :: numGlobalElements_local, numGlobalElements_return
  integer(c_int) :: Index_Base=1
  real(c_double) :: bnorm(1), xnorm(1)
  real(c_double) :: err_tol,expected_bnorm,expected_xnorm,bnorm_err,xnorm_err 
  real(c_double) :: two = 2.0, zero = 0.0
  logical        :: success = .true.,zero_initial=.true.
  
  if (.not. valid_kind_parameters()) stop 'C interoperability not supported on this platform.'
  
  ! Executable code
  
! Create a serial comm
  communicator= Epetra_SerialComm() 

! Create a map 
  numGlobalElements_local = 4 
  map = Epetra_Map(numGlobalElements_local,Index_Base,communicator)
  numGlobalElements_return = map%NumGlobalElements()
  print *,'NumGlobalElements = ', numGlobalElements_return
  print *,'NumMyElements=', map%NumMyElements()
  if ( numGlobalElements_local /= numGlobalElements_return ) &
    stop 'In ForTrilinos (verySimpleObjectOriented.F90: return mismatch'
   
  ! Create vectors
  x = Epetra_Vector(map,zero_initial)
  b = Epetra_Vector(map,zero_initial)
 
  ! Do some vector operations
  call b%PutScalar(two)
  call x%Update(two, b, zero) ! /* x = 2*b */
 
  bnorm = b%Norm2(err)
  print *,'Error value from b%Norm2():',err%error_code()
  print *,'Error message from b%Norm2():',err%text()
  xnorm = x%Norm2(err)
  print *,'Error value from x%Norm2():',err%error_code()
  print *,'Error message from x%Norm2():',err%text()
 
  print *, "2 norm of x = ", xnorm(1) 
  print *, "2 norm of b = ", bnorm(1) 

! Test the expected value 
  err_tol = 1.0e-14
  expected_bnorm = sqrt( 2.0 * 2.0 * numGlobalElements_return )
  expected_xnorm = sqrt( 4.0 * 4.0 * numGlobalElements_return )
  bnorm_err = abs( expected_bnorm - bnorm(1) ) / expected_bnorm
  xnorm_err = abs( expected_xnorm - xnorm(1) ) / expected_xnorm
  print *, "error in 2 norm of x = ",bnorm_err
  print *, "error in 2 norm of b = ",xnorm_err
  if (bnorm_err > err_tol) success = .false.
  if (xnorm_err > err_tol) success = .false.
 
  ! Clean up memory (in reverse order).  This step is not required
  ! with compilers that support Fortran 2003 type finalization:
  call b%force_finalize()
  call x%force_finalize()
  call map%force_finalize()
  call communicator%force_finalize()
 
  if (success) then
    print *  
    print *, "End Result: TEST PASSED" 
  else
    print *  
    print *, "End Result: TEST FAILED"
  end if
end program main
