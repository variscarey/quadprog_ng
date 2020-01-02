module quadprog_ng
  implicit none
contains
  subroutine solve_qp(quadr_coeff_G, linear_coeff_a,
                      n_ineq, ineq_coef_C, ineq_vec_d,
                      m_eq, eq_coef_A, eq_vec_b,
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
    real(8), allocatable, intent(in) :: eq_coef_A
    real(8), allocatable, intent(in) :: eq_vec_b

    integer, intent(in) :: nvars

    ! the solution iterate
    real(8), allocatable, intent(inout) :: sol(:)

    ! If ierr is set to anything except for zero, a problem happened
    integer, intent(inout) :: ierr = 0    

    !!
    !! Internals
    !!
    logical DONE = .FALSE.
    logical FULL_STEP = .FALSE. 
    logical ADDING_EQ_CONSTRAINTS = .TRUE.

    logical first_pass = .TRUE. 

    integer status = 0
    integer irow, icol, iactive_set, icopy_idx = 1

    real(8), allocatable :: G_inv(:,:)
    real(8), allocatable :: U_work(:,:)

    real(8), allocatable :: ineq_prb(:)
    real(8), allocatable :: eq_prb(:)

    !! Cholesky decomp of quadr_coeff_G
    real(8), allocatable :: chol_L(:,:)
    real(8), allocatable :: inv_chol_L(:,:)

    !! QR factorization of B = L^{-1} N
    real(8), allocatable :: Q(:,:)
    real(8), allocatable :: R(:,:)
    real(8), allocatable :: R_inv(:,:)

    real(8), allocatable :: J(:,:)

    integer, allocatable :: active_set(:)
    integer, allocatable :: n_p(:)
    integer, allocatable :: copy_integer(:)
    integer :: p = 0, &
               q = 0 

    !! lagrangian, each step in dual space
    real(8), allocatable :: u(:)
    !! lagrangian for each constraint in the active set
    real(8), allocatable :: lagr(:)
    real(8), allocatable :: copy_real(:)

    integer :: k_dropped = 0, &
               j_dropped = 0, &
               k = 0

    real(8), allocatable :: z(:), &  ! Step direction in primal space
                            r(:)     ! Step dir, dual space

    real(8) :: t1, t2, t
    real(8) :: MAX_DOUBLE = huge(t1)

    !! intermediate matrices for doing the inversions
    real(8), allocatable :: work(:)
    integer, allocatable :: ipiv(:)

    !! temps for holding 
    integer, dimension(2) :: R_dim

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

    ! Lower triangular cholesky 
    call dpotrf('L', nvars, L, nvars, status)

    ! Set all non-lower-triangular entries to 0
    do icol=1,nvars
        do irow=1,nvars
            if (irow .lt. icol) then
                L(irow, icol) = 0
            endif
        enddo
    enddo

    allocate(G_inv(nvars,nvars))

    G_inv = L

    ! calc G^{-1}
    call dpotri('L', nvars, G_inv, nvars, status)

    allocate(U_work(nvars,nvars))

    U_work = G

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

    allocate(z(nvars))
    allocate(r(m_eq + n_ineq))

    allocate(copy_integer(m_eq + n_ineq))
    copy_integer = 0
    allocate(copy_real(m_eq + n_ineq))
    copy_real = 0

    allocate(J(m_eq + n_ineq, m_eq + n_ineq))

    !!~~~ Begin Processing ~~~!!
    !! Solution iterate
    !! Start soln iterate at unconstrained minimum
    sol = (-1) * matmul(G_inv, linear_coeff_a)

    !! Main loop
    !! ###~~~~~~~~ Step 1 ~~~~~~~~###
    do while (.not. DONE)
      ineq_prb = matmul(transpose(ineq_coef_C),sol) - ineq_vec_d

      if (ADDING_EQ_CONSTRAINTS) then
        ineq_prb = matmul(transpose(ineq_coef_C),sol) - ineq_vec_d
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
        do irow=1,q
          lagr(irow) = u(irow)
        enddo

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
            z = matmul(G_inv, n_p)

            first_pass = .false.
          else
            z = matmul(matmul(J2, transpose(J2)), n_p)

            if (q .gt. 0) then
              R_dim = shape(R)
              allocate(R_inv(R_dim(1), R_dim(2)))
              allocate(tau(q))
              allocate(ipiv(q))

              R_inv = R
              call dgetrf(q,q,inv_chol_L,q,ipiv,info)
              call dgetri(q,inv_chol_L,q,ipiv,work,q,info)

              deallocate(ipiv)
              deallocate(work)

              r = matmul(matmul(R_inv(1:q,1:q), transpose(J1)), n_p)

              deallocate(R_inv)
            endif
          endif

          !!###~~~~~~~~ Step 2(b) ~~~~~~~~###
          !!# partial step length t1 - max step in dual space
          if ((q .eq. 0) .or. (r .le. 0) .or. ADDING_EQ_CONSTRAINTS) then
            t1 = MAX_DOUBLE 
          else
            t1 = MAX_DOUBLE
            k_dropped = 0

            do iactive_set=meq+1, q
              k = active_set(iactive_set)
              if ((r(iactive_set) .gt. 0) .and. ((lagr(iactive_set) / r(iactive_set)) .lt. t1)) then
                t1 = lagr(iactive_set) / r(iactive_set)
                k_dropped = k
                j_dropped = iactive_set
              endif
            enddo
          endif

          !!# full step length t2 - min step in primal space
          if (all(z .eq. 0)) then
            t2 = MAX_DOUBLE
          else
            if (ADDING_EQ_CONSTRAINTS) then
              t2 = (-1) * eq_prb(p) / (matmul(transpose(z), n_p))
            else
              t2 = (-1) * ineq_prb / (matmul(transpose(z), n_p))
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
            r(q+1) = -1
            lagr = lagr + (-1 * t * r)

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
            active_set = copy_integer

            q = q - 1

            !TODO:> ADD QR UPDATE


            !# go back to step 2(a)
            cycle
          endif

          sol = sol + t * z
          r(q+1) = -1
          lagr = lagr + (-1 * t * r)

          !!# if we took a full step
          if (t .eq. t2) then
            q = q+1
            active_set(q) = p

            u = 0
            icopy_idx = 1
            do iactive_set=1,q
              if(lagr(iactive_set) .gt. 0) then
                u(icopy_idx) = lagr(iactive_set)
                icopy_idx = icopy_idx + 1
              endif
            enddo

            !TODO:> ADD QR_UPDATE
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
            active_set = copy_integer

            q = q - 1

            !TODO:> ADD QR UPDATE
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

end program test