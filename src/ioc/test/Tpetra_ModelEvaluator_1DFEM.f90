module ThyraToTpetraModelEvaluator_class

! Module defines the base model evaluator type.  Derived types must:
!
!  1) call Fortrilinosmodelevaluator%setup(...)
!  2) Provide an implementations of
!     - Modelevaluator_eval_resid
!     - Modelevaluator_eval_jac
!     - Modelevaluator_eval_prec
!
use forteuchos
use fortpetra

#include "ForTrilinos.h"

use ISO_FORTRAN_ENV
use, intrinsic :: ISO_C_BINDING

#ifdef HAVE_MPI
  use mpi
#endif

implicit none

type, abstract :: ThyraToTpetraModelEvaluator
  ! FIXME (TJF May 2018): This class should be replaced by the SWIG generated
  ! parent class
  contains
    procedure, public :: evaluate_model => ThyraToTpetraModelEvaluator_eval_model
    procedure, public :: setup => ThyraToTpetraModelEvaluator_setup

    procedure(eval_resid_interface), deferred :: ModelEvaluator_eval_resid
    generic :: evaluate_residual => Modelevaluator_eval_resid

    procedure(update_x_interface), deferred :: ModelEvaluator_update_x
    generic :: update_x => Modelevaluator_update_x

    procedure, public :: evaluate_jacobian => ModelEvaluator_eval_jac
    procedure, public :: evaluate_preconditioner => ModelEvaluator_eval_prec

end type

interface
  subroutine eval_resid_interface(this, f)
    import
    class(ThyraToTpetraModelEvaluator) :: this
    type(TpetraMultiVector) :: f
  end subroutine eval_resid_interface

  subroutine update_x_interface(this, x)
    import
    class(ThyraToTpetraModelEvaluator) :: this
    type(TpetraMultiVector) :: x
  end subroutine update_x_interface
end interface

contains

! ---------------------------------------------------------------------------- !

subroutine ThyraToTpetraModelEvaluator_setup(this, owned_map, ghosted_map)
  ! -------------------------------------------------------------------------- !
  class(ThyraToTpetraModelEvaluator) :: this
  type(TpetraMap) :: owned_map
  type(TpetraMap), optional :: ghosted_map
  ! -------------------------------------------------------------------------- !
  ! FIXME (TJF May 2018): This procedure will be provided by the actual SWIG
  ! generated parent class.  A dummy implementation is provided so that the
  ! "user" model evaluator can be fully implemented.
end subroutine ThyraToTpetraModelEvaluator_setup

! ---------------------------------------------------------------------------- !

subroutine Modelevaluator_eval_jac(this, op)
  class(ThyraToTpetraModelEvaluator) :: this
  type(TpetraCrsMatrix) :: op
  print*, 'evaluate_jacobian must be implemented by derived class'
  stop 666
end subroutine Modelevaluator_eval_jac

! ---------------------------------------------------------------------------- !

subroutine ModelEvaluator_eval_prec(this, op)
  class(ThyraToTpetraModelEvaluator) :: this
  type(TpetraCrsMatrix) :: op
  print*, 'evaluate_preconditioner must be implemented by derived class'
  stop 666
end subroutine ModelEvaluator_eval_prec

! ---------------------------------------------------------------------------- !

subroutine ThyraToTpetraModelEvaluator_eval_model(this, x, fill_f, fill_J, fill_prec)
  ! -------------------------------------------------------------------------- !
  class(ThyraToTpetraModelEvaluator) :: this
  logical, intent(in) :: fill_f, fill_J, fill_prec
  type(TpetraMultiVector) :: x
  type(TpetraMultiVector) :: f
  type(TpetraCrsMatrix) :: J
  type(TpetraCrsMatrix) :: M_inv
  ! -------------------------------------------------------------------------- !

  call this%update_x(x)

  if (fill_f) then
    call this%evaluate_residual(f)
  end if

  if (fill_J) then
    call this%evaluate_jacobian(J)
  end if

  if (fill_prec) then
    call this%evaluate_preconditioner(M_inv)
  end if

end subroutine ThyraToTpetraModelEvaluator_eval_model

! ---------------------------------------------------------------------------- !

end module ThyraToTpetraModelEvaluator_class


module TpetraModelEvaluator1DFEM_module

use ThyraToTpetraModelEvaluator_class
use forteuchos
use fortpetra

#include "ForTrilinos.h"

use ISO_FORTRAN_ENV
use, intrinsic :: ISO_C_BINDING

#ifdef HAVE_MPI
  use mpi
#endif

implicit none

real(scalar_type), parameter :: zero=0., one=1.

type, extends(ThyraToTpetraModelEvaluator) :: TpetraModelEvaluator1DFEM
  type(TeuchosComm), private :: comm
  type(TpetraMap), private :: x_owned_map, x_ghosted_map, f_owned_map
  type(TpetraImport), private :: importer
  type(TpetraMultiVector), private :: node_coords
  type(TpetraMultiVector), private, allocatable :: u
  type(TpetraMultiVector), private, allocatable :: x
  type(TpetraMultiVector), private, allocatable :: J_diagonal

  contains
    procedure, private :: create_mesh
    procedure, private :: create_graph
    procedure :: ModelEvaluator_update_x => TpetraModelEvaluator1DFEM_update_x
    procedure :: ModelEvaluator_eval_resid => TpetraModelEvaluator1DFEM_eval_resid
    procedure :: ModelEvaluator_eval_jac => TpetraModelEvaluator1DFEM_eval_jac
    procedure :: ModelEvaluator_eval_prec => TpetraModelEvaluator1DFEM_eval_prec

end type

interface TpetraModelEvaluator1DFEM
  procedure new_TpetraModelEvaluator1DFEM
end interface

! Linear FE basis
type :: Linear2NodeFEBasis
  real(scalar_type) :: phi(2), dphide(2)
  real(scalar_type) :: uu, zz, duu, eta, wt
  real(scalar_type) :: dz
  ! These are only needed for transient
  real(scalar_type) :: uuold, duuold;

  contains
    ! Calculates the values of u and x at the specified gauss point
    procedure, public :: compute_basis
end type

interface Linear2NodeFEBasis
  procedure new_Linear2NodeFEBasis
end interface

contains

! ---------------------------------------------------------------------------- !

function new_TpetraModelEvaluator1DFEM(comm, num_global_elems, z_min, z_max) &
    result(this)
  ! -------------------------------------------------------------------------- !
  type(TpetraModelEvaluator1DFEM) :: this
  type(TeuchosComm), intent(in) :: comm
  integer(global_size_type), intent(in) :: num_global_elems
  real(scalar_type), intent(in) :: z_min, z_max
  integer :: i
  integer(global_size_type) :: num_nodes, invalid
  integer(global_ordinal_type) :: min_overlap_GID, gid
  integer(size_type) :: num_overlap_nodes
  integer(global_ordinal_type), allocatable :: node_gids(:)
  type(TpetraCrsGraph) :: graph
  ! -------------------------------------------------------------------------- !

  this%comm = comm
  num_nodes = num_global_elems + 1;

  ! owned space
  this%x_owned_map = TpetraMap(num_nodes, comm)

  ! ghosted space
  if (comm%getSize() == 1) then
    this%x_ghosted_map = this%x_owned_map
  else
    num_overlap_nodes = this%x_owned_map%getNodeNumElements() + 2
    if ((comm%getRank() == 0) .or. (comm%getRank() == (comm%getSize() - 1))) &
      num_overlap_nodes = num_overlap_nodes - 1
    if (comm%getRank() == 0) then
      min_overlap_GID = this%x_owned_map%getMinGlobalIndex()
    else
      min_overlap_GID = this%x_owned_map%getMinGlobalIndex() - 1
    end if

    allocate(node_gids(num_overlap_nodes))
    gid = min_overlap_GID
    do i=1, num_overlap_nodes
      node_gids(i) = gid
      gid = gid + 1
    end do

    invalid = -1
    this%x_ghosted_map = TpetraMap(invalid, node_gids, comm)
  end if

  this%importer = TpetraImport(this%x_owned_map, this%x_ghosted_map)

  ! residual space
  this%f_owned_map = this%x_owned_map

  ! Initialize the graph for W CrsMatrix object
  graph = this%create_graph(this%x_owned_map, this%x_ghosted_map)

  ! Create the nodal coordinates
  this%node_coords = &
    this%create_mesh(this%x_owned_map, z_min, z_max, num_global_elems)

  call this%setup(this%x_owned_map)

end function new_TpetraModelEvaluator1DFEM

! ---------------------------------------------------------------------------- !

function new_Linear2NodeFEBasis() result(this)
  ! -------------------------------------------------------------------------- !
  type(Linear2NodeFEBasis) :: this
  ! -------------------------------------------------------------------------- !
  this%uu = 0
  this%zz = 0
  this%duu = 0
  this%eta = 0
  this%wt = 0
  this%dz = 0
  this%uuold = 0
  this%duuold = 0
end function new_Linear2NodeFEBasis

! ---------------------------------------------------------------------------- !

subroutine TpetraModelEvaluator1DFEM_update_x(this, x)
  ! -------------------------------------------------------------------------- !
  class(TpetraModelEvaluator1DFEM) :: this
  type(TpetraMultiVector) :: x
  integer(size_type) :: num_vecs=1
  ! -------------------------------------------------------------------------- !
  ! Create ghosted objects
  if (.not. allocated(this%u)) then
    allocate(this%u, source=TpetraMultiVector(this%x_ghosted_map, num_vecs))
  end if
  call this%u%doImport(x, this%importer, TpetraREPLACE)

  if (.not. allocated(this%x)) then
    allocate(this%x, source=TpetraMultiVector(this%x_ghosted_map, num_vecs))
    call this%x%doImport(this%node_coords, this%importer, TpetraINSERT)
  end if
end subroutine TpetraModelEvaluator1DFEM_update_x

! ---------------------------------------------------------------------------- !

function create_graph(this, owned_map, ghosted_map) result(graph)
  ! -------------------------------------------------------------------------- !
  class(TpetraModelEvaluator1DFEM) :: this
  type(TpetraCrsGraph) :: graph
  type(TpetraMap) :: owned_map, ghosted_map
  integer(local_ordinal_type) :: ne, i, j
  integer(global_ordinal_type) :: gblrow, gblinds(1)
  integer(size_type) :: num_my_overlap_nodes, num_ent_per_row
  ! -------------------------------------------------------------------------- !

  ! Create the shell for the graph
  num_ent_per_row = 5
  graph = TpetraCrsGraph(owned_map, ghosted_map, num_ent_per_row, TpetraDynamicProfile)

  ! Declare required variables
  num_my_overlap_nodes = ghosted_map%getNodeNumElements()

  ! Loop Over # of Finite Elements on Processor
  do ne=1, num_my_overlap_nodes-1

    ! Loop over nodes in element
    do i=1, 2

      gblrow = ghosted_map%getGlobalElement(ne+i)

      ! Loop over Trial Functions
      do j=1, 2

        ! If this row is owned by current processor, add the index
        if (owned_map%isNodeGlobalElement(gblrow)) then
          gblinds(1) = ghosted_map%getGlobalElement(ne+j)
          call graph%insertGlobalIndices(gblrow, gblinds)
        end if

      end do
    end do
  end do

  call graph%fillComplete()

end function create_graph

! ---------------------------------------------------------------------------- !

function create_mesh(this, owned_map, z_min, z_max, num_elems) result(coords)
  ! -------------------------------------------------------------------------- !
  class(TpetraModelEvaluator1DFEM) :: this
  type(TpetraMultiVector) :: coords
  type(TpetraMap), intent(in) :: owned_map
  real(scalar_type), intent(in) :: z_min, z_max
  integer(global_size_type), intent(in) :: num_elems
  integer(local_ordinal_type) :: lclrow
  integer(size_type) :: num_local_nodes, num_vecs=1
  integer(global_ordinal_type) :: min_GID, col=1
  real(scalar_type) :: dz
  ! -------------------------------------------------------------------------- !
  num_local_nodes = owned_map%getNodeNumElements()
  min_GID = owned_map%getMinGlobalIndex()
  dz = (z_max - z_min) / num_elems
  coords = TpetraMultiVector(owned_map, num_vecs);
  do lclrow=1, num_local_nodes
    call coords%replaceLocalValue(lclrow, col, z_min + dz * (min_GID+lclrow))
  end do
end function create_mesh

! ---------------------------------------------------------------------------- !

subroutine TpetraModelEvaluator1DFEM_eval_resid(this, f)
  ! -------------------------------------------------------------------------- !
  class(TpetraModelEvaluator1DFEM) :: this
  type(TpetraMultiVector) :: f
  type(Linear2NodeFEBasis) :: basis
  integer(local_ordinal_type) :: ne, invalid, num_my_elems, gp, i, lclrow
  integer(size_type) :: col
  integer :: my_rank
  integer(size_type), parameter :: ione=1
  real(scalar_type), dimension(:), pointer :: xdata
  real(scalar_type), dimension(:), pointer :: udata
  real(scalar_type) :: xx(2), uu(2), val
  ! -------------------------------------------------------------------------- !

  call f%putScalar(zero)

  invalid = -1
  my_rank = this%comm%getRank()
  num_my_elems = this%x_ghosted_map%getNodeNumElements()-1

  ! Loop Over # of Finite Elements on Processor
  xdata => this%x%getData(ione)
  udata => this%u%getData(ione)

  do ne = 1, num_my_elems

    ! Get the solution and coordinates at the nodes
    xx(1) = xdata(ne)
    xx(2) = xdata(ne+1)

    uu(1) = udata(ne)
    uu(2) = udata(ne+1)

    basis = Linear2NodeFEBasis()

    ! Loop Over Gauss Points
    do gp=1, 2

      ! Calculate the basis function at the gauss point
      call basis%compute_basis(gp, xx, uu)

      ! Loop over Nodes in Element
      do i = 1, 2
        lclrow = this%x_owned_map%getLocalElement(this%x_ghosted_map%getGlobalElement(ne+i))
        if (lclrow /= invalid) then
          val = basis%wt * basis%dz * &
            (basis%uu * basis%uu * basis%phi(i) + &
            (basis%duu * basis%dphide(i))/(basis%dz * basis%dz))
          col = 1
          call f%sumIntoLocalValue(lclrow, col, val)
        end if
      end do
    end do

    ! Correct for Dirichlet BCs
    if ((my_rank == 0) .and. (ne == 0)) then
      lclrow = 1; col = 1
      val = udata(1) - 1.0
      call f%replaceLocalValue(lclrow, col, val)
    end if

  end do

  nullify(xdata)
  nullify(udata)

end subroutine TpetraModelEvaluator1DFEM_eval_resid

! ---------------------------------------------------------------------------- !

subroutine TpetraModelEvaluator1DFEM_eval_jac(this, Jac)
  ! -------------------------------------------------------------------------- !
  class(TpetraModelEvaluator1DFEM) :: this
  type(TpetraCrsMatrix), intent(inout) :: Jac
  type(Linear2NodeFEBasis) :: basis
  integer(local_ordinal_type) :: ne, invalid, num_my_elems, gp, i, j
  integer(local_ordinal_type) :: lclrow, lclcol, numvalid
  integer(global_ordinal_type) :: gblrow, cols(1)
  real(scalar_type), dimension(:), pointer :: xdata
  real(scalar_type), dimension(:), pointer :: udata
  integer(size_type), parameter :: ione=1
  integer :: my_rank
  real(scalar_type) :: xx(2), uu(2), vals(1)
  ! -------------------------------------------------------------------------- !
  call Jac%setAllToScalar(zero)

  invalid = -1
  my_rank = this%comm%getRank()
  num_my_elems = this%x_ghosted_map%getNodeNumElements()-1

  ! Loop Over # of Finite Elements on Processor
  xdata => this%x%getData(ione)
  udata => this%u%getData(ione)

  do ne=1, num_my_elems

    ! Get the solution and coordinates at the nodes
    xx(1) = xdata(ne)
    xx(2) = xdata(ne+1)

    uu(1) = udata(ne)
    uu(2) = udata(ne+1)

    basis = Linear2NodeFEBasis()

    ! Loop Over Gauss Points
    do gp=1, 2

      ! Calculate the basis function at the gauss point
      call basis%compute_basis(gp, xx, uu)

      ! Loop over Nodes in Element
      do i=1, 2
        lclrow = this%x_owned_map%getLocalElement(this%x_ghosted_map%getGlobalElement(ne+i));
        if (lclrow /= invalid) then
          ! Loop over trial functions
          do j=1, 2
            lclcol = ne + j;
            vals(1) = basis%wt * basis%dz &
              * ((basis%dphide(j)*basis%dphide(i))/(basis%dz*basis%dz) &
              + 2.0*basis%uu*basis%phi(j)*basis%phi(i))
            ! FIXME (TJF: May 2018) Should replace *GlobalValues with ! *LocalValues
            gblrow = this%x_owned_map%getGlobalElement(lclrow)
            cols(1) = this%x_ghosted_map%getGlobalElement(lclcol)
            numvalid = Jac%sumIntoGlobalValues(gblrow, cols, vals)
          end do
        end if
      end do
    end do

    ! Correct for Dirichlet BCs
    if ((my_rank == 0) .and. (ne == 1)) then
      gblrow = 1;
      cols(1) = 1;
      vals(1) = 1.0;
      ! FIXME (TJF: May 2018) Should replace *GlobalValues with ! *LocalValues
      numvalid = Jac%replaceGlobalValues(gblrow, cols, vals)
      cols(1) = 2;
      vals(1) = 0.0;
      ! FIXME (TJF: May 2018) Should replace *GlobalValues with ! *LocalValues
      numvalid = Jac%replaceGlobalValues(gblrow, cols, vals)
    end if
  end do

  nullify(xdata)
  nullify(udata)

  call Jac%fillComplete()

end subroutine TpetraModelEvaluator1DFEM_eval_jac

! ---------------------------------------------------------------------------- !

subroutine TpetraModelEvaluator1DFEM_eval_prec(this, M)
  ! -------------------------------------------------------------------------- !
  class(TpetraModelEvaluator1DFEM) :: this
  type(TpetraCrsMatrix) :: M
  type(TpetraMap) :: row_map, col_map
  type(Linear2NodeFEBasis) :: basis
  integer(local_ordinal_type) :: ne, invalid, num_my_elems, gp, i, j
  integer(local_ordinal_type) :: lclrow, lclcol, numvalid
  integer(global_ordinal_type) :: gblrow, cols(1)
  integer(size_type), parameter :: ione=1
  real(scalar_type), dimension(:), pointer :: xdata
  real(scalar_type), dimension(:), pointer :: udata
  type(TpetraMultiVector) :: diag
  integer(size_type) :: num_vecs=1
  integer :: my_rank
  real(scalar_type) :: xx(2), uu(2), vals(1)
  ! -------------------------------------------------------------------------- !

  if (.not. allocated(this%J_diagonal)) then
    allocate(this%J_diagonal, source=TpetraMultiVector(this%x_owned_map, num_vecs))
  end if
  call M%setAllToScalar(zero)
  call this%J_diagonal%putScalar(zero)

  invalid = -1
  my_rank = this%comm%getRank()
  num_my_elems = this%x_ghosted_map%getNodeNumElements()-1
  row_map = M%getRowMap()
  col_map = M%getColMap()

  ! Loop Over # of Finite Elements on Processor
  xdata => this%x%getData(ione)
  udata => this%u%getData(ione)

  do ne=1, num_my_elems

    ! Get the solution and coordinates at the nodes
    xx(1) = xdata(ne)
    xx(2) = xdata(ne+1)

    uu(1) = udata(ne)
    uu(2) = udata(ne+1)

    basis = Linear2NodeFEBasis()

    ! Loop Over Gauss Points
    do gp=1, 2

      ! Calculate the basis function at the gauss point
      call basis%compute_basis(gp, xx, uu)

      ! Loop over Nodes in Element
      do i=1, 2
        lclrow = this%x_owned_map%getLocalElement(this%x_ghosted_map%getGlobalElement(ne+i))
        if (lclrow /= invalid) then
          ! Loop over trial functions
          do j=1, 2
            lclcol = ne + j
            if (row_map%getGlobalElement(lclrow) == col_map%getGlobalElement(lclcol)) then
              vals(1) = basis%wt * basis%dz &
                * ((basis%dphide(j)*basis%dphide(i))/(basis%dz*basis%dz) &
                + 2.0*basis%uu*basis%phi(j)*basis%phi(i))
              gblrow = this%x_owned_map%getGlobalElement(lclrow)
              cols(1) = this%x_ghosted_map%getGlobalElement(lclcol)
              ! FIXME (TJF: May 2018) Should replace *GlobalValues with ! *LocalValues
              numvalid = M%sumIntoGlobalValues(gblrow, cols, vals)
            end if
          end do
        end if
      end do
    end do

    ! Correct for Dirichlet BCs
    if ((my_rank == 0) .and. (ne == 0)) then
      gblrow = 1;
      cols(1) = 1;
      vals(1) = 1.0;
      ! FIXME (TJF: May 2018) Should replace *GlobalValues with ! *LocalValues
      numvalid = M%replaceGlobalValues(gblrow, cols, vals)
    end if
  end do

  nullify(xdata)
  nullify(udata)

  ! Invert the Jacobian diagonal for the preconditioner
  ! For some reason the matrix must be fill complete before calling rightScale
  call M%fillComplete()
  diag = this%J_diagonal
  ! FIXME (TJF May 2018): getLocalDiagCopy is not implemented!
  ! call M%getLocalDiagCopy(diag);
  call diag%reciprocal(diag)
  ! FIXME (TJF May 2018): rightScale is not implemented!
  ! call M%rightScale(diag)

end subroutine TpetraModelEvaluator1DFEM_eval_prec

! ---------------------------------------------------------------------------- !

subroutine compute_basis(this, gp, z, u, uold)
  ! -------------------------------------------------------------------------- !
  ! Calculates the values of z and u at the specified Gauss point
  ! -------------------------------------------------------------------------- !
  class(Linear2NodeFEBasis) :: this
  integer(local_ordinal_type), intent(in) :: gp
  real(scalar_type), intent(in) :: z(2), u(2)
  real(scalar_type), intent(in), optional :: uold(2)
  real(scalar_type) :: eta, wt
  integer(local_ordinal_type) :: i

  ! -------------------------------------------------------------------------- !
  if (gp==1) then
    eta = -1.0 / sqrt(3.0)
    wt = 1.0
  else if (gp==2) then
    eta = 1.0 / sqrt(3.0)
    wt = 1.0
  end if

  ! Calculate basis function and derivatives at Gauss point
  this%phi(1) = (1.0 - eta) / 2.0
  this%phi(2) = (1.0 + eta) / 2.0
  this%dphide(1) = -0.5;
  this%dphide(2) = -this%dphide(1)

  ! Caculate function and derivative approximations at GP.
  this%dz = 0.5 * (z(2) - z(1))
  this%zz = 0.0
  this%uu = 0.0
  this%duu = 0.0
  this%uuold = 0.0
  this%duuold = 0.0

  do i = 1, 2
    this%zz = this%zz + z(i) * this%phi(i)
    this%uu = this%uu + u(i) * this%phi(i)
    this%duu = this%duu + u(i) * this%dphide(i)
    if (present(uold)) then
      this%uuold = this%uuold + uold(i) * this%phi(i);
      this%duuold = this%duuold + uold(i) * this%dphide(i);
    end if
  end do
end subroutine compute_basis

end module TpetraModelEvaluator1DFEM_module
