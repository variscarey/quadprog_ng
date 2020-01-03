module quadprog_ng
  implicit none
contains
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

    logical :: first_pass = .false. 

    integer :: status = 0
    integer :: irow, icol, iactive_set, icopy_idx = 1

    real(8), allocatable :: G_inv(:,:)
    real(8), allocatable :: U_work(:,:)

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
      ADDING_EQ_CONSTRAINTS = .true.
    else
      ADDING_EQ_CONSTRAINTS = .false. 
    endif

    if (.not. (allocated(sol))) then
      allocate(sol(nvars))
    endif

    !! Begin chol factorization and use the factors to get
    !! the inverse of the matrix G

    allocate(chol_L(nvars, nvars))
    allocate(inv_chol_L(nvars, nvars))

    chol_L = quadr_coeff_G

    ! Lower triangular cholesky 
    call dpotrf('L', nvars, chol_L, nvars, status)

    ! Set all non-lower-triangular entries to 0
    do icol=1,nvars
        do irow=1,nvars
            if (irow .lt. icol) then
                chol_L(irow, icol) = 0
            endif
        enddo
    enddo

    allocate(G_inv(nvars,nvars))

    G_inv = chol_L

    ! calc G^{-1}
    call dpotri('L', nvars, G_inv, nvars, status)

    allocate(U_work(nvars,nvars))

    U_work = quadr_coeff_G

    call dpotrf('U', 4, U_work, 4, info)
    call dpotri('U', 4, U_work, 4, info)

    do icol=1,nvars
        do irow=1,nvars
            if (irow .ge. icol) then
                U_work(irow, icol) = 0
            endif
        enddo
    enddo

    !! Now hold onto L inverse, use LU factorization on it
    allocate(tau(nvars))
    allocate(ipiv(nvars))

    inv_chol_L = chol_L
    call dgetrf(nvars,nvars,inv_chol_L,nvars,ipiv,info)
    call dgetri(nvars,inv_chol_L,nvars,ipiv,work,nvars,info)

    deallocate(chol_L)  ! Bounce the lower triangular, don't need it
    deallocate(ipiv)
    deallocate(work)

    !! Begin adding the lower and upper terms
    do icol=1,nvars
      do irow=1,nvars
        G_inv(irow, icol) = G_inv(irow,icol) + U_work(irow,icol)
      enddo
    enddo

    !! Bounce the U_work matrix so we can re-use it later
    deallocate(U_work)

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
          p = q
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
              allocate(work(q))
              allocate(ipiv(q))

              R_inv(1:q,1:q) = R(1:q,1:q)
              call dgetrf(q,q,R_inv(1:q,1:q),q,ipiv,info)
              call dgetri(q,R_inv(1:q,1:q),q,ipiv,work,q,info)

              deallocate(ipiv)
              deallocate(work)

              r_step = matmul(matmul(R_inv(1:q,1:q), transpose(J(:,1:q))), n_p)

              deallocate(R_inv)
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

            allocate(tau(q))
            allocate(work(q))

            call dgeqrf(nvars, q, Q_mat(1:nvars,1:q), nvars, tau, work, q, info)

            R(1:q,1:q) = Q_mat(1:q, 1:q)

            call dorgqr(q, q, q, Q_mat(1:q,1:q), q, tau, work, q, info)

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

            allocate(tau(q))
            allocate(work(q))

            call dgeqrf(nvars, q, Q_mat(1:nvars,1:q), nvars, tau, work, q, info)

            R(1:q,1:q) = Q_mat(1:q, 1:q)

            call dorgqr(q, q, q, Q_mat(1:q,1:q), q, tau, work, q, info)

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

            allocate(tau(q))
            allocate(work(q))

            call dgeqrf(nvars, q, Q_mat(1:nvars,1:q), nvars, tau, work, q, info)

            R(1:q,1:q) = Q_mat(1:q, 1:q)

            call dorgqr(q, q, q, Q_mat(1:q,1:q), q, tau, work, q, info)

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