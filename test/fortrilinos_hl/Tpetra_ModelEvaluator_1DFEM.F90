! Copyright 2017-2018, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE

! Module defines the base model evaluator type.  Derived types must:
!
!  1) call Fortrilinosmodelevaluator%setup(...)
!  2) Provide an implementations of
!     - Modelevaluator_eval_resid
!     - Modelevaluator_eval_jac
!     - Modelevaluator_eval_prec
!
module TpetraModelEvaluator1DFEM_module
  use forteuchos
  use fortpetra
  use fortrilinos_hl
  implicit none

  real(scalar_type), parameter :: zero=0., one=1.
  type(TpetraMultiVector), private, allocatable :: solnvec

  type, extends(ForModelEvaluator) :: TpetraModelEvaluator1DFEM
    type(TeuchosComm), private :: comm
    type(TpetraMap), private :: x_owned_map, x_ghosted_map, f_owned_map
    type(TpetraCrsGraph), private :: graph
    type(TpetraImport), private :: importer
    type(TpetraMultiVector), private :: node_coords
    type(TpetraMultiVector), private :: x
    type(TpetraMultiVector), private :: J_diagonal

  contains
    procedure, private :: create_mesh
    procedure, private :: create_graph
    procedure :: evaluate_residual => TpetraModelEvaluator1DFEM_eval_resid
    procedure :: evaluate_jacobian => TpetraModelEvaluator1DFEM_eval_jac
    procedure :: evaluate_preconditioner => TpetraModelEvaluator1DFEM_eval_prec
    procedure :: update_solution_vector => TpetraModelEvaluator1DFEM_update_solution_vector
    procedure :: get_x_map => TpetraModelEvaluator1DFEM_get_x_map
    procedure :: get_f_map => TpetraModelEvaluator1DFEM_get_f_map
    procedure :: create_operator => TpetraModelEvaluator1DFEM_create_operator
    procedure :: release => delete_TpetraModelEvaluator1DFEM
  end type

  interface TpetraModelEvaluator1DFEM
    procedure new_TpetraModelEvaluator1DFEM
  end interface

  ! Linear FE basis
  type :: Linear2NodeFEBasis
    real(scalar_type) :: phi(2), dphide(2)
    real(scalar_type) :: uu, zz, duu, wt
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

  ! -------------------------------------------------------------------------- !

  function new_Linear2NodeFEBasis() result(self)
    type(Linear2NodeFEBasis) :: self
    ! ------------------------------------------------------------------------ !
    self%uu = 0
    self%zz = 0
    self%duu = 0
    self%wt = 0
    self%dz = 0
    self%uuold = 0
    self%duuold = 0
  end function new_Linear2NodeFEBasis

  ! -------------------------------------------------------------------------- !

  subroutine compute_basis(self, gp, z, u, uold)
    ! ------------------------------------------------------------------------ !
    ! Calculates the values of z and u at the specified Gauss point
    ! ------------------------------------------------------------------------ !
    class(Linear2NodeFEBasis) :: self
    integer, intent(in) :: gp
    real(scalar_type), intent(in) :: z(2), u(2)
    real(scalar_type), intent(in), optional :: uold(2)
    real(scalar_type) :: eta
    integer :: i

    if (gp==1) then
      eta = -1.0 / sqrt(3.0)
      self%wt = 1.0
    else if (gp==2) then
      eta = 1.0 / sqrt(3.0)
      self%wt = 1.0
    end if

    ! Calculate basis function and derivatives at Gauss point
    self%phi(1) = (1.0 - eta) / 2.0
    self%phi(2) = (1.0 + eta) / 2.0
    self%dphide(1) = -0.5;
    self%dphide(2) = -self%dphide(1)

    ! Caculate function and derivative approximations at GP.
    self%dz = 0.5 * (z(2) - z(1))
    self%zz = 0.0
    self%uu = 0.0
    self%duu = 0.0
    self%uuold = 0.0
    self%duuold = 0.0

    do i = 1, 2
    self%zz = self%zz + z(i) * self%phi(i)
    self%uu = self%uu + u(i) * self%phi(i)
    self%duu = self%duu + u(i) * self%dphide(i)
    if (present(uold)) then
      self%uuold = self%uuold + uold(i) * self%phi(i);
      self%duuold = self%duuold + uold(i) * self%dphide(i);
    end if
    end do
  end subroutine compute_basis

  ! -------------------------------------------------------------------------- !

  function new_TpetraModelEvaluator1DFEM(comm, num_global_elems, z_min, z_max) &
      result(self)
    ! ------------------------------------------------------------------------ !
    type(TpetraModelEvaluator1DFEM) :: self
    type(TeuchosComm), intent(in) :: comm
    integer(global_size_type), intent(in) :: num_global_elems
    real(scalar_type), intent(in) :: z_min, z_max
    integer :: i, num_overlap_nodes
    integer(global_size_type) :: num_nodes
    integer(global_ordinal_type) :: min_overlap_GID, gid
    integer(size_type) :: num_vecs=1
    integer(global_ordinal_type), allocatable :: node_gids(:)
    ! ------------------------------------------------------------------------ !

    self = ForModelEvaluator()

    self%comm = comm
    num_nodes = num_global_elems + 1;

    ! owned space
    self%x_owned_map = TpetraMap(num_nodes, comm)

    ! ghosted space
    if (comm%getSize() == 1) then
      self%x_ghosted_map = self%x_owned_map
    else
      num_overlap_nodes = self%x_owned_map%getNodeNumElements() + 2
      if ((comm%getRank() == 0) .or. (comm%getRank() == (comm%getSize() - 1))) &
        num_overlap_nodes = num_overlap_nodes - 1
      if (comm%getRank() == 0) then
        min_overlap_GID = self%x_owned_map%getMinGlobalIndex()
      else
        min_overlap_GID = self%x_owned_map%getMinGlobalIndex() - 1
      end if

      allocate(node_gids(num_overlap_nodes))
      gid = min_overlap_GID
      do i=1, num_overlap_nodes
        node_gids(i) = gid
        gid = gid + 1
      end do

      self%x_ghosted_map = TpetraMap(TPETRA_GLOBAL_INVALID, node_gids, comm)
      deallocate(node_gids)
    end if

    self%importer = TpetraImport(self%x_owned_map, self%x_ghosted_map)

    ! residual space
    self%f_owned_map = self%x_owned_map

    ! Initialize the graph for W CrsMatrix object
    self%graph = self%create_graph(self%x_owned_map, self%x_ghosted_map)

    ! Create the nodal coordinates
    self%node_coords = &
      self%create_mesh(self%x_owned_map, z_min, z_max, num_global_elems)

    ! Allocate space for domain and solution vectors
    self%x = TpetraMultiVector(self%x_ghosted_map, num_vecs)
    self%J_diagonal = TpetraMultiVector(self%x_owned_map, num_vecs)

    call self%x%doImport(self%node_coords, self%importer, TpetraINSERT)

  end function new_TpetraModelEvaluator1DFEM

  ! -------------------------------------------------------------------------- !

  function create_graph(self, owned_map, ghosted_map) result(graph)
    ! ------------------------------------------------------------------------ !
    class(TpetraModelEvaluator1DFEM) :: self
    type(TpetraCrsGraph) :: graph
    type(TpetraMap) :: owned_map, ghosted_map
    integer :: ne, num_my_overlap_nodes, i, j
    integer(global_ordinal_type) :: gblrow, gblinds(1)
    integer(size_type) :: num_ent_per_row
    ! ------------------------------------------------------------------------ !

    ! Create the shell for the graph
    num_ent_per_row = 5
    graph = TpetraCrsGraph(owned_map, ghosted_map, num_ent_per_row, TpetraStaticProfile)

    ! Declare required variables
    num_my_overlap_nodes = ghosted_map%getNodeNumElements()

    ! Loop Over # of Finite Elements on Processor
    do ne=1, num_my_overlap_nodes-1

      ! Loop over nodes in element
      do i=1, 2

        gblrow = ghosted_map%getGlobalElement(ne+i-1)

        ! Loop over Trial Functions
        do j=1, 2

          ! If self row is owned by current processor, add the index
          if (owned_map%isNodeGlobalElement(gblrow)) then
            gblinds(1) = ghosted_map%getGlobalElement(ne+j-1)
            call graph%insertGlobalIndices(gblrow, gblinds)
          end if

        end do
      end do
    end do

    call graph%fillComplete()

  end function create_graph

  ! -------------------------------------------------------------------------- !

  function create_mesh(self, owned_map, z_min, z_max, num_elems) result(coords)
    ! ------------------------------------------------------------------------ !
    class(TpetraModelEvaluator1DFEM) :: self
    type(TpetraMultiVector) :: coords
    type(TpetraMap), intent(in) :: owned_map
    real(scalar_type), intent(in) :: z_min, z_max
    integer(global_size_type), intent(in) :: num_elems
    integer(size_type) :: num_vecs=1
    integer :: lclrow, num_local_nodes
    integer(global_ordinal_type) :: min_GID, col=1
    real(scalar_type) :: dz
    ! ------------------------------------------------------------------------ !

    num_local_nodes = owned_map%getNodeNumElements()
    min_GID = owned_map%getMinGlobalIndex()
    dz = (z_max - z_min) / num_elems
    coords = TpetraMultiVector(owned_map, num_vecs)
    do lclrow=1, num_local_nodes
      call coords%replaceLocalValue(lclrow, col, z_min + dz * (min_GID+lclrow))
    end do
  end function create_mesh

  ! -------------------------------------------------------------------------- !

  subroutine TpetraModelEvaluator1DFEM_update_solution_vector(self, xp)
    ! ------------------------------------------------------------------------ !
    class(TpetraModelEvaluator1DFEM), intent(in) :: self
    class(TpetraMultiVector), intent(in) :: xp
    integer(size_type) :: num_vecs=1
    ! ------------------------------------------------------------------------ !
    if (.not.allocated(solnvec)) then
      allocate(solnvec, source=TpetraMultiVector(self%x_ghosted_map, num_vecs))
    endif
    call solnvec%doImport(xp, self%importer, TpetraREPLACE)
  end subroutine TpetraModelEvaluator1DFEM_update_solution_vector

  ! -------------------------------------------------------------------------- !

  subroutine TpetraModelEvaluator1DFEM_eval_resid(self, x, f)
    ! ------------------------------------------------------------------------ !
    class(TpetraModelEvaluator1DFEM), intent(in) :: self
    class(TpetraMultiVector), intent(in) :: x
    class(TpetraMultiVector), intent(in) :: f
    type(Linear2NodeFEBasis) :: basis
    integer :: num_my_elems, ne, gp, i, lclrow
    integer(size_type) :: col
    integer :: my_rank
    integer(size_type), parameter :: ione=1
    real(scalar_type), dimension(:), pointer :: xdata
    real(scalar_type), dimension(:), pointer :: udata
    real(scalar_type) :: xx(2), uu(2), val
    ! ------------------------------------------------------------------------ !

    call self%update_solution_vector(x)

    call f%putScalar(zero)

    my_rank = self%comm%getRank()
    num_my_elems = self%x_ghosted_map%getNodeNumElements()-1

    ! Loop Over # of Finite Elements on Processor
    xdata => self%x%getData(ione)
    udata => solnvec%getData(ione)

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
    lclrow = self%x_owned_map%getLocalElement(self%x_ghosted_map%getGlobalElement(ne+i-1))
    if (lclrow /= TPETRA_LOCAL_INVALID) then
      val = basis%wt * basis%dz * &
        (basis%uu * basis%uu * basis%phi(i) + &
        (basis%duu * basis%dphide(i))/(basis%dz * basis%dz))
      col = 1
      call f%sumIntoLocalValue(lclrow, col, val)
    end if
    end do
    end do

    ! Correct for Dirichlet BCs
    if ((my_rank == 0) .and. (ne == 1)) then
      lclrow = 1; col = 1
      val = udata(1) - 1.0
      call f%replaceLocalValue(lclrow, col, val)
    end if

    end do

    nullify(xdata)
    nullify(udata)

  end subroutine TpetraModelEvaluator1DFEM_eval_resid

  ! -------------------------------------------------------------------------- !

  subroutine TpetraModelEvaluator1DFEM_eval_jac(self, x, J)
    ! ------------------------------------------------------------------------ !
    class(TpetraMultiVector), intent(in) :: x
    class(TpetraModelEvaluator1DFEM), intent(in) :: self
    class(TpetraOperator), intent(in) :: J
    type(TpetraCrsMatrix) :: Jmat
    type(Linear2NodeFEBasis) :: basis
    integer :: ne, num_my_elems, gp, i, jj
    integer :: lclrow, lclcol, numvalid
    integer(global_ordinal_type) :: gblrow, cols(1)
    real(scalar_type), dimension(:), pointer :: xdata
    real(scalar_type), dimension(:), pointer :: udata
    integer(size_type), parameter :: ione=1
    integer :: my_rank
    real(scalar_type) :: xx(2), uu(2), vals(1)
    ! ------------------------------------------------------------------------ !

    call self%update_solution_vector(x)

    Jmat = operator_to_matrix(J)
    call Jmat%resumeFill()

    call Jmat%setAllToScalar(zero)

    my_rank = self%comm%getRank()
    num_my_elems = self%x_ghosted_map%getNodeNumElements()-1

    ! Loop Over # of Finite Elements on Processor
    xdata => self%x%getData(ione)
    udata => solnvec%getData(ione)

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
    lclrow = self%x_owned_map%getLocalElement(self%x_ghosted_map%getGlobalElement(ne+i-1));
    if (lclrow /= TPETRA_LOCAL_INVALID) then
      ! Loop over trial functions
      do jj=1, 2
      lclcol = ne + jj - 1;
      vals(1) = basis%wt * basis%dz &
        * ((basis%dphide(jj)*basis%dphide(i))/(basis%dz*basis%dz) &
        + 2.0*basis%uu*basis%phi(jj)*basis%phi(i))
      ! FIXME (TJF: May 2018) Should replace *GlobalValues with ! *LocalValues
      gblrow = self%x_owned_map%getGlobalElement(lclrow)
      cols(1) = self%x_ghosted_map%getGlobalElement(lclcol)
      numvalid = Jmat%sumIntoGlobalValues(gblrow, cols, vals)
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
      numvalid = Jmat%replaceGlobalValues(gblrow, cols, vals)
      cols(1) = 2;
      vals(1) = 0.0;
      ! FIXME (TJF: May 2018) Should replace *GlobalValues with ! *LocalValues
      numvalid = Jmat%replaceGlobalValues(gblrow, cols, vals)
    end if
    end do

    nullify(xdata)
    nullify(udata)

    call Jmat%fillComplete()

  end subroutine TpetraModelEvaluator1DFEM_eval_jac

  ! -------------------------------------------------------------------------- !

  subroutine TpetraModelEvaluator1DFEM_eval_prec(self, x, M)
    ! ------------------------------------------------------------------------ !
    class(TpetraModelEvaluator1DFEM), intent(in) :: self
    class(TpetraMultiVector), intent(in) :: x
    class(TpetraOperator), intent(in) :: M
    type(TpetraCrsMatrix) :: Mmat
    type(TpetraMap) :: row_map, col_map
    type(Linear2NodeFEBasis) :: basis
    integer :: ne, num_my_elems, gp, i, j
    integer :: lclrow, lclcol, numvalid
    integer(global_ordinal_type) :: gblrow, cols(1)
    integer(size_type), parameter :: ione=1
    real(scalar_type), dimension(:), pointer :: xdata
    real(scalar_type), dimension(:), pointer :: udata
    integer :: my_rank
    real(scalar_type) :: xx(2), uu(2), vals(1)
    ! ------------------------------------------------------------------------ !

    call self%update_solution_vector(x)

    Mmat = operator_to_matrix(M)
    call Mmat%resumeFill()

    call Mmat%setAllToScalar(zero)
    ! FIXME
    ! call self%J_diagonal%putScalar(zero)

    my_rank = self%comm%getRank()
    num_my_elems = self%x_ghosted_map%getNodeNumElements()-1
    row_map = Mmat%getRowMap()
    col_map = Mmat%getColMap()

    ! Loop Over # of Finite Elements on Processor
    xdata => self%x%getData(ione)
    udata => solnvec%getData(ione)

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
    lclrow = self%x_owned_map%getLocalElement(self%x_ghosted_map%getGlobalElement(ne+i-1))
    if (lclrow /= TPETRA_LOCAL_INVALID) then
      ! Loop over trial functions
      do j=1, 2
      lclcol = ne + j
      if (row_map%getGlobalElement(lclrow) == col_map%getGlobalElement(lclcol)) then
        vals(1) = basis%wt * basis%dz &
          * ((basis%dphide(j)*basis%dphide(i))/(basis%dz*basis%dz) &
          + 2.0*basis%uu*basis%phi(j)*basis%phi(i))
        gblrow = self%x_owned_map%getGlobalElement(lclrow)
        cols(1) = self%x_ghosted_map%getGlobalElement(lclcol)
        ! FIXME (TJF: May 2018) Should replace *GlobalValues with ! *LocalValues
        numvalid = Mmat%sumIntoGlobalValues(gblrow, cols, vals)
      end if
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
      numvalid = Mmat%replaceGlobalValues(gblrow, cols, vals)
    end if
    end do

    nullify(xdata)
    nullify(udata)

    ! Invert the Jacobian diagonal for the preconditioner
    ! For some reason the matrix must be fill complete before calling rightScale
    call Mmat%fillComplete()
    ! FIXME
    ! diag = self%J_diagonal
    ! FIXME (TJF May 2018): getLocalDiagCopy is not implemented!
    ! call Mmat%getLocalDiagCopy(diag);
    ! call diag%reciprocal(diag)
    ! FIXME (TJF May 2018): rightScale is not implemented!
    ! call Mmat%rightScale(diag)
    ! call Mmat%rightScale(diag)

  end subroutine TpetraModelEvaluator1DFEM_eval_prec

  ! -------------------------------------------------------------------------- !

  function TpetraModelEvaluator1DFEM_get_x_map(self) result(map)
    class(TpetraModelEvaluator1DFEM), intent(in) :: self
    type(TpetraMap) :: map

    map = self%x_owned_map
  end function TpetraModelEvaluator1DFEM_get_x_map

  ! -------------------------------------------------------------------------- !

  function TpetraModelEvaluator1DFEM_get_f_map(self) result(map)
    class(TpetraModelEvaluator1DFEM), intent(in) :: self
    type(TpetraMap) :: map

    map = self%x_owned_map
  end function TpetraModelEvaluator1DFEM_get_f_map

  ! -------------------------------------------------------------------------- !

  function TpetraModelEvaluator1DFEM_create_operator(self) result(op)
    class(TpetraModelEvaluator1DFEM), intent(in) :: self
    type(TpetraCrsMatrix) :: matrix
    type(TpetraOperator) :: op

    matrix = TpetraCrsMatrix(self%graph)

    op = matrix_to_operator(matrix)

  end function TpetraModelEvaluator1DFEM_create_operator

  ! -------------------------------------------------------------------------- !

  subroutine delete_TpetraModelEvaluator1DFEM(self)
    class(TpetraModelEvaluator1DFEM), intent(inout) :: self

    call self%comm%release()
    call self%x_owned_map%release()
    call self%x_ghosted_map%release()
    call self%f_owned_map%release()
    call self%graph%release()
    call self%importer%release()
    call self%node_coords%release()
    call self%x%release()
    call self%J_diagonal%release()

#ifdef __GNUC__
    ! FIXME This segfaults with Flang
    ! Call base class release()
    call self%ForModelEvaluator%release()
#endif
  end subroutine

end module TpetraModelEvaluator1DFEM_module
