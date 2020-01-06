module quadprog_ng
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

    lwork = floor(1.5 * rank_A)
    allocate(work(lwork))

    call dgetri(rank_A, out_mat_A_Inv, rank_A, ipiv, work, lwork, ierr)

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
    ! 
    real(8), allocatable, intent(in) :: quadr_coeff_G(:,:)
    real(8), allocatable, intent(in) :: linear_coeff_a(:)
    
    integer, intent(in) :: n_ineq
    real(8), allocatable, intent(in) :: ineq_coef_C(:,:)
    real(8), allocatable, intent(in) :: ineq_vec_d(:)

    integer, intent(in) :: m_eq
    real(8), allocatable, intent(in) :: eq_coef_A(:,:)
    real(8), allocatable, intent(in) :: eq_vec_b(:)

    integer, intent(in) :: nvars

    ! the solution iterate
    real(8), allocatable, intent(inout) :: sol(:)

    ! If ierr is set to anything except for zero, a problem happened
    integer, intent(inout) :: ierr 
    integer :: info

    !!
    !! Internals
    !!
    logical :: DONE = .false.
    logical :: FULL_STEP = .false. 
    logical :: ADDING_EQ_CONSTRAINTS = .false.

    logical :: first_pass = .true. 

    integer :: status = 0
    integer :: irow, icol, iactive_set, icopy_idx = 1

    real(8), allocatable :: G_inv(:,:)

    real(8), allocatable :: ineq_prb(:)
    real(8), allocatable :: eq_prb(:)

    !! Cholesky decomp of quadr_coeff_G
    real(8), allocatable :: chol_L(:,:)
    real(8), allocatable :: inv_chol_L(:,:)

    !! QR factorization of B = L^{-1} N
    real(8), allocatable :: Q_mat(:,:)
    real(8), allocatable :: R(:,:)
    real(8), allocatable :: R_inv(:,:)

    real(8), allocatable :: J(:,:)
    real(8), allocatable :: B(:,:)

    integer, allocatable :: active_set(:)
    integer, allocatable :: n_p(:)
    integer, allocatable :: copy_integer(:)
    integer :: p, q = 0 

    !! lagrangian, each step in dual space
    real(8), allocatable :: u(:)
    !! lagrangian for each constraint in the active set
    real(8), allocatable :: lagr(:)

    integer :: k_dropped, j_dropped, k = 0

    real(8), allocatable :: z_step(:), r_step(:)

    real(8) :: t1, t2, t
    real(8) :: MAX_DOUBLE = huge(t1)

    !! intermediate matrices for doing the inversions
    real(8), allocatable :: work(:), tau(:)
    integer, allocatable :: ipiv(:)

    !! temps for holding shape information
    integer, dimension(2) :: mat_dim

    !!~~~ Allocations & Initializations ~~~!!
    if (m_eq .eq. 0) then
      ADDING_EQ_CONSTRAINTS = .false.
    else
      ADDING_EQ_CONSTRAINTS = .true. 
    endif

    if (.not. (allocated(sol))) then
      allocate(sol(nvars))
    endif

    !! Begin chol factorization and use the factors to get
    !! the inverse of the matrix G

    allocate(chol_L(nvars, nvars))
    allocate(inv_chol_L(nvars, nvars))


    !! HAD OLD CHOLESKY STUFF HERE
    print *, "WARNING 1:> ADD CHOLESKY STUFF BACK HERE"


    !! Allocate the other internals
    allocate(n_p(nvars))
    allocate(active_set(m_eq + n_ineq))  
    allocate(lagr(m_eq + n_ineq))
    allocate(u(m_eq + n_ineq))

    allocate(z_step(nvars))
    allocate(r_step(m_eq + n_ineq))

    allocate(copy_integer(m_eq + n_ineq))
    copy_integer = 0

    allocate(Q_mat(nvars,nvars))
    allocate(R(m_eq + n_ineq, m_eq + n_ineq))
    allocate(R_inv(m_eq + n_ineq, m_eq + n_ineq))

    allocate(J(nvars,nvars))
    J = 0

    allocate(B(nvars, m_eq + n_ineq))
    B = 0

    !!~~~ Begin Processing ~~~!!
    !! Solution iterate
    !! Start soln iterate at unconstrained minimum
    sol = (-1) * matmul(G_inv, linear_coeff_a)

    !! Main loop
    !! ###~~~~~~~~ Step 1 ~~~~~~~~###
    do while (.not. DONE)
      ineq_prb = matmul(transpose(ineq_coef_C),sol) - ineq_vec_d

      if (ADDING_EQ_CONSTRAINTS) then
        eq_prb = matmul(transpose(eq_coef_A),sol) - eq_vec_b
      endif

      if (q .ge. m_eq) then
        ADDING_EQ_CONSTRAINTS = .false. 
      endif

      if (any(ineq_prb < 0) .or. ADDING_EQ_CONSTRAINTS) then
        if (ADDING_EQ_CONSTRAINTS) then
          p = q+1
          n_p = eq_coef_A(:,p)
        else
          do icol=1,n_ineq
            if (ineq_prb(icol) .lt. 0) then
              p = icol
              exit
            endif
          enddo
          n_p = ineq_coef_C(:,p)
        endif

        if (q .eq. 0) then
          u = 0
        endif

        lagr = 0
        lagr(1:q) = u(1:q)

        FULL_STEP = .false. 

        !!###~~~~~~~~ Step 2 ~~~~~~~~###      
        do while (.not. FULL_STEP)
          ineq_prb = matmul(transpose(ineq_coef_C),sol) - ineq_vec_d

          if (ADDING_EQ_CONSTRAINTS) then
            ineq_prb = matmul(transpose(ineq_coef_C),sol) - ineq_vec_d
          endif

          !!###~~~~~~~~ Step 2(a) ~~~~~~~~###
          !!## Calculate step directions
          if (first_pass) then
            z_step = matmul(G_inv, n_p)

            first_pass = .false.
          else
            z_step = matmul(matmul(J(:,q+1:nvars), transpose(J(:,q+1:nvars))), n_p)

            if (q .gt. 0) then
              allocate(work(nvars))
              allocate(ipiv(nvars))

              info = 0
              R_inv(1:q,1:q) = R(1:q,1:q)
              
              call dgetrf(q,q,R_inv,q,ipiv,info)
              call dgetri(q,R_inv,q,ipiv,work,nvars,info)

              deallocate(ipiv)
              deallocate(work)

              r_step = matmul(matmul(R_inv(1:q,1:q), transpose(J(:,1:q))), n_p)
            endif
          endif

          !!###~~~~~~~~ Step 2(b) ~~~~~~~~###
          !!# partial step length t1 - max step in dual space
          if ( ((q .eq. 0) .or. all(r_step .le. 0)) .or. ADDING_EQ_CONSTRAINTS) then
            t1 = MAX_DOUBLE 
          else
            t1 = MAX_DOUBLE
            k_dropped = 0

            do iactive_set=m_eq, q
              k = active_set(iactive_set)
              if ((r_step(iactive_set) .gt. 0) .and. ((lagr(iactive_set) / r_step(iactive_set)) .lt. t1)) then
                t1 = lagr(iactive_set) / r_step(iactive_set)
                k_dropped = k
                j_dropped = iactive_set + m_eq
              endif
            enddo
          endif

          !!# full step length t2 - min step in primal space
          if (all(z_step .eq. 0)) then
            t2 = MAX_DOUBLE
          else
            if (ADDING_EQ_CONSTRAINTS) then
              t2 = (-1) * eq_prb(p) / (dot_product(z_step, n_p))
            else
              t2 = (-1) * ineq_prb(p) / ((dot_product(z_step, n_p)))
            endif
          endif

          !!# current step length
          t = min(t1, t2)

          !!###~~~~~~~~ Step 2(c) ~~~~~~~~###
          if(t .eq. MAX_DOUBLE) then
            print *, "infeasible! Stop here!"
            FULL_STEP = .true.
            DONE = .true.
            ierr = 420
            return
          endif

          !!# If t2 is infinite, then we took a partial step in the dual space.
          if (t2 .eq. MAX_DOUBLE) then
            !update lagrangian
            r_step(q+1) = -1
            lagr = lagr + (-1 * t * r_step)

            !! remove dropped constraint
            active_set(j_dropped) = 0
            icopy_idx = 1
            copy_integer = 0
            do iactive_set=1,q
              if (active_set(iactive_set) .gt. 0) then
                copy_integer(icopy_idx) = active_set(iactive_set)
                icopy_idx = icopy_idx + 1
              endif
            enddo
            active_set = copy_integer

            !TODO:> ADD QR UPDATE
            B(:,j_dropped) = 0
            B(:,j_dropped:q) = B(:,j_dropped+1:q+1)

            q = q - 1

            R = 0
            Q_mat = 0
            Q_mat(1:nvars, 1:q) = B(1:nvars,1:q)

            allocate(tau(nvars))
            allocate(work(nvars))

            call dgeqrf(nvars, q, Q_mat(1:nvars,1:q), nvars, tau, work, nvars, info) !!TODO:> figure out what should happen when q = 0

            R(1:q,1:q) = Q_mat(1:q, 1:q)

            call dorgqr(nvars, nvars, nvars, Q_mat(1:q,1:q), nvars, tau, work, nvars, info)

            deallocate(tau)
            deallocate(work)

            !! zero out lower entries of R
            do icol=1,q
              do irow=1,q
                if (irow .gt. icol) then
                  R(irow,icol) = 0
                endif
              enddo
            enddo

            J = 0
            J = matmul(transpose(inv_chol_L), Q_mat)

            !# go back to step 2(a)
            cycle
          endif

          sol = sol + t * z_step
          r_step(q+1) = -1
          lagr = lagr + (-1 * t * r_step)

          !!# if we took a full step
          if (t .eq. t2) then
            q = q+1
            active_set(q) = p

            u = 0
            u(1:q) = lagr(1:q)

            !TODO:> ADD QR_UPDATE
            B(:,q) = matmul(inv_chol_L, n_p)

            R = 0
            Q_mat = 0
            Q_mat(1:nvars, 1:q) = B(1:nvars,1:q)

            allocate(tau(nvars))
            allocate(work(nvars))

            call dgeqrf(nvars, q, Q_mat(1:nvars,1:q), nvars, tau, work, nvars, info)

            R(1:q,1:q) = Q_mat(1:q, 1:q)

            call dorgqr(q, q, q, Q_mat, nvars, tau, work, nvars, info)

            deallocate(tau)
            deallocate(work)

            !! zero out lower entries of R
            do icol=1,q
              do irow=1,q
                if (irow .gt. icol) then
                  R(irow,icol) = 0
                endif
              enddo
            enddo

            J = 0
            J = matmul(transpose(inv_chol_L), Q_mat)

            FULL_STEP = .true.
            exit
          endif

          if (t .eq. t1) then
            !! remove dropped constraint
            active_set(j_dropped) = 0
            icopy_idx = 1
            copy_integer = 0
            do iactive_set=1,q+1
              if (active_set(iactive_set) .gt. 0) then
                copy_integer(icopy_idx) = active_set(iactive_set)
                icopy_idx = icopy_idx + 1
              endif
            enddo
            active_set = 0
            active_set = copy_integer

            !TODO:> ADD QR UPDATE
            B(:,j_dropped) = 0
            B(:,j_dropped:q) = B(:,j_dropped+1:q+1)

            q = q - 1

            R = 0
            Q_mat = 0
            Q_mat(1:nvars, 1:q) = B(1:nvars,1:q)

            allocate(tau(nvars))
            allocate(work(nvars))

            call dgeqrf(nvars, q, Q_mat(1:nvars,1:q), nvars, tau, work, nvars, info)

            R(1:q,1:q) = Q_mat(1:q, 1:q)

            call dorgqr(nvars, nvars, nvars, Q_mat(1:q,1:q), nvars, tau, work, nvars, info)

            deallocate(tau)
            deallocate(work)

            !! zero out lower entries of R
            do icol=1,q
              do irow=1,q
                if (irow .gt. icol) then
                  R(irow,icol) = 0
                endif
              enddo
            enddo

            J = 0
            J = matmul(transpose(inv_chol_L), Q_mat)

            cycle
          endif

        enddo
      else
        DONE = .true.
      endif
    enddo 

    deallocate(n_p)
    deallocate(active_set)  
    deallocate(lagr)
    deallocate(u)
    deallocate(z_step)
    deallocate(r_step)
    deallocate(copy_integer)
    deallocate(Q_mat)
    deallocate(R)
    deallocate(R_inv)
    deallocate(J)
    deallocate(B)

    return

  end subroutine solve_qp
end module


program test
  use quadprog_ng
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