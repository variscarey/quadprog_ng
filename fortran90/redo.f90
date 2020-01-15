module quadprog_ng_redo
  implicit none
contains
  subroutine do_cholesky_and_inverse(rank_A, in_mat_A, out_mat_L, out_mat_Inv)
    implicit none
    integer, intent(in) :: rank_A
    real(8), allocatable, intent(in) :: in_mat_A(:,:)
    real(8), allocatable, intent(out) :: out_mat_L(:,:), out_mat_Inv(:,:)

    integer :: irow, icol = 0
    real(8), allocatable :: mat_U(:,:)

    integer :: ierr

    allocate(mat_U(rank_A, rank_A))

    if (.not. allocated(out_mat_L)) then
      allocate(out_mat_L(rank_A, rank_A))
      out_mat_L = 0
    endif

    if (.not. allocated(out_mat_Inv)) then
      allocate(out_mat_Inv(rank_A, rank_A))
      out_mat_Inv = 0
    endif

    out_mat_L = in_mat_A
    mat_U = in_mat_A

    call dpotrf('L', rank_A, out_mat_L, rank_A, ierr)
    call dpotrf('U', rank_A, mat_U, rank_A, ierr)

    out_mat_Inv = out_mat_L

    call dpotri('L', rank_A, out_mat_Inv, rank_A, ierr)
    call dpotri('U', rank_A, mat_U, rank_A, ierr)

    !! Zero out bad entries in upper and lower
    do icol=1,rank_A
        do irow=1,rank_A
            if (irow .lt. icol) then
                out_mat_Inv(irow, icol) = 0
                out_mat_L(irow, icol) = 0
            endif
        enddo
    enddo

    do icol=1,rank_A
        do irow=1,rank_A
            if (irow .ge. icol) then
                mat_U(irow, icol) = 0
            endif
        enddo
    enddo

    out_mat_Inv = out_mat_Inv + mat_U

    return
  end subroutine

  subroutine get_inverse(rank_A, in_mat_A, out_mat_A_Inv)
    implicit none
    integer, intent(in) :: rank_A
    real(8), allocatable, intent(in) :: in_mat_A(:,:)
    real(8), allocatable, intent(out) :: out_mat_A_Inv(:,:)

    integer, allocatable :: ipiv(:)
    real(8), allocatable :: work(:)

    integer :: ierr
    integer :: lwork

    allocate(ipiv(rank_A))

    if (.not. allocated(out_mat_A_Inv)) then
      allocate(out_mat_A_Inv(rank_A, rank_A))
    endif

    out_mat_A_Inv = in_mat_A

    call dgetrf(rank_A, rank_A, out_mat_A_Inv, rank_A, ipiv, ierr)

    lwork = ceiling(1.5 * rank_A)
    allocate(work(lwork))

    call dgetri(rank_A, out_mat_A_Inv, rank_A, ipiv, work, lwork, ierr)

    deallocate(ipiv)
    deallocate(work)
  end subroutine

  subroutine get_qr(nrow, ncol, in_mat_A, mat_Q, mat_R)
    implicit none
    integer, intent(in) :: nrow, ncol
    real(8), allocatable, intent(in) :: in_mat_A(:,:)
    real(8), allocatable, intent(out) :: mat_Q(:,:), mat_R(:,:)

    real(8), allocatable :: work(:), tau(:), temp(:), temp_R(:,:)
    integer :: lwork, ierr, irow, icol

    allocate(tau(min(nrow, ncol)))
    tau = 0
    allocate(temp(1))
    temp = 0

    print *, "doing allocations"

    if (.not. allocated(mat_Q)) then
      allocate(mat_Q(nrow, nrow))
      mat_Q = 0
    endif

    print *, "allocated q"

    if (.not. allocated(mat_R)) then
      allocate(mat_R(nrow, ncol))
      mat_R = 0
    endif

    print *, "allocated r"

    mat_R = 0
    mat_R(1:nrow, 1:ncol) = in_mat_A(1:nrow, 1:ncol)

    !! Do a dummy call to find optimal lwork value
    call dgeqrf(nrow, ncol, mat_R, nrow, tau, temp, -1, ierr)

    lwork = int(temp(1))
    allocate(work(lwork))

    print *, "starting R"

    !! Form R
    call dgeqrf(nrow, ncol, mat_R, nrow, tau, work, lwork, ierr)

    print *, "doing temp move"

    allocate(temp_R(nrow, max(nrow, ncol)))
    temp_R = 0
    temp_R(1:nrow,1:ncol) = mat_R(1:nrow, 1:ncol)

    print *, "start Q"

    !! Get Q back from it
    call dorgqr(nrow, nrow, nrow, temp_R, nrow, tau, work, lwork, ierr)

    mat_Q = 0
    mat_Q(1:nrow, 1:ncol) = temp_R(1:nrow, 1:nrow)

    print *, "done Q"

    !! zero out bad entries to make R upper triangular
    do icol=1,ncol
        do irow=1,nrow
            if (irow .gt. icol) then
                mat_R(irow, icol) = 0
            endif
        enddo
    enddo

    deallocate(temp_R)
    deallocate(work)
    deallocate(tau)
    deallocate(temp)
  end subroutine


  subroutine solve_qp(quadr_coeff_G, linear_coeff_a, &
                      n_ineq, ineq_coef_C, ineq_vec_d, &
                      m_eq, eq_coef_A, eq_vec_b, &
                      nvars, sol, ierr)
    !
    implicit none

    !!
    !! Externals
    !! 

    ! quadratic coefficient matrix G in 
    real(8), allocatable, intent(in) :: quadr_coeff_G(:,:)
    real(8), allocatable, intent(in) :: linear_coeff_a(:)
    
    integer, intent(in) :: n_ineq
    real(8), allocatable, intent(in) :: ineq_coef_C(:,:)
    real(8), allocatable, intent(in) :: ineq_vec_d(:)

    integer, intent(in) :: m_eq
    real(8), allocatable, intent(in) :: eq_coef_A(:,:)
    real(8), allocatable, intent(in) :: eq_vec_b(:)

    integer, intent(in) :: nvars
    integer, intent(inout) :: ierr

    ! the solution iterate
    real(8), allocatable, intent(inout) :: sol(:)

    !!~~~~~~~~~~ Internals ~~~~~~~~~~!!
    logical :: DONE, FULL_STEP = .false.
    logical :: FIRST_PASS = .true. 

    integer :: irow, icol

    real(8), allocatable :: G_inv(:,:)
    real(8), allocatable :: L_chol(:,:), L_inv(:,:)

    real(8) :: t, t1, t2
    real(8) :: MAX_DOUBLE = huge(t)

    !!~~~~~~~~ Allocations ~~~~~~~~!!
    allocate(L_chol(nvars,nvars))
    allocate(L_inv(nvars,nvars))
    allocate(G_inv(nvars,nvars))
    
    call do_cholesky_and_inverse(nvars, quadr_coeff_G, L_chol, G_inv)
    call get_inverse(nvars, L_chol, L_inv)

    sol = (-1) * matmul(G_inv, linear_coeff_a)

    deallocate(L_chol)
    deallocate(G_inv)

  end subroutine solve_qp
end module


program test
  use quadprog_ng_redo
  implicit none
  real(8), allocatable :: G(:,:), lin_vec(:), C(:,:), d(:), A(:,:), b(:), sol(:)
  integer :: irow, icol
  integer :: m_eq, n_ineq, ierr

  integer, dimension(2) :: mat_dim

  G = transpose(reshape((/4, -2, -2, 4/),(/2,2/)))

  lin_vec = (/6, 0/)

  C = transpose(reshape((/1, 0, 1, 0, 1, 1/),(/3,2/)))

  d = (/0, 0, 2/)

  allocate(A(1,1))
  A = 0

  allocate(b(1))
  b = 1

  m_eq = 0
  n_ineq = 3

  allocate(sol(2))

  call solve_qp(G, lin_vec, &
                      n_ineq, C, d, &
                      m_eq, A, b, &
                      2, sol, ierr)

  print *, sol

end program test