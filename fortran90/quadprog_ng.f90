module quadprog_ng
  implicit none
contains
  subroutine solve_qp(quadr_coeff_G, linear_coeff_a,
                      n_ineq, ineq_coef_C, ineq_vec_d,
                      m_eq, eq_coef_A, eq_vec_b,
                      sol, ierr)
    !
    implicit none
    !
    ! Externals
    ! 
    real(8), allocatable, intent(in) :: quadr_coeff_G(:,:)
    real(8), allocatable, intent(in) :: linear_coeff_a(:)
    
    integer, intent(in) :: n_ineq
    real(8), allocatable, intent(in) :: ineq_coef_C(:,:)
    real(8), allocatable, intent(in) :: ineq_vec_d(:)

    integer, intent(in) :: m_eq
    real(8), allocatable, intent(in) :: eq_coef_A
    real(8), allocatable, intent(in) :: eq_vec_b

    real(8), allocatable, intent(out) :: sol(:)
    integer, intent(out) :: ierr = 0    
    ! 
    ! Internals
    !
    logical DONE = .FALSE.
    logical FULL_STEP = .FALSE. 
    logical ADDING_EQ_CONSTRAINTS = .TRUE.

    logical first_pass = .TRUE. 

    !!~~~ Begin Processing ~~~!!
    !! Solution iterate

  end subroutine solve_qp
end module


program test

end program test