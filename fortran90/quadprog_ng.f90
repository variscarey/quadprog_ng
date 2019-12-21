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

    integer :: intent(in) :: nvars

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

    !! Cholesky decomp of quadr_coeff_G
    real(8), allocatable :: chol_L(:,:)
    real(8), allocatable :: inv_chol_L(:,:)

    !! QR factorization of B = L^{-1} N
    real(8), allocatable :: Q(:,:)
    real(8), allocatable :: R(:,:)

    !! J = L^{-T} Q, inverse transpose of L by columns of Q = [Q1 | Q2]
    !! Q1 has the columns corresponding to active constraints 
    real(8), allocatable :: J1(:,:)
    real(8), allocatable :: J2(:,:)

    integer, allocatable :: active_set(:)
    integer, allocatable :: n_p(:)
    integer :: p, &
               q 

    real(8), allocatable :: u(:)
    real(8), allocatable :: lagr(:)

    integer :: k_dropped, &
               j_dropped

    real(8), allocatable :: z(:), & 
                            r(:)

    !!~~~ Allocations ~~~!!

    !!~~~ Begin Processing ~~~!!
    !! Solution iterate


  end subroutine solve_qp
end module


program test

end program test