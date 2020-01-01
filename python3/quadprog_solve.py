#!/usr/bin/python3
import numpy as np
from scipy.linalg import qr_insert, qr_delete

###-------------------------------------###
def quadprog_solve(quadr_coeff_G, linear_coeff_a, 
                   n_ineq, ineq_coef_C, ineq_vec_d,
                   m_eq, eq_coef_A, eq_vec_b):
    DONE = False
    FULL_STEP = False

    if eq_coef_A is not None:
        ADDING_EQ_CONSTRAINTS = True
    else:
        ADDING_EQ_CONSTRAINTS = False

    first_pass = True

    ###-------------------------------------###
    ## Solution iterate
    sol = (-1) * np.linalg.inv(quadr_coeff_G) * linear_coeff_a   

    ## We only need to keep ahold of L^{-1} for the implementation. 
    ## L is the triangular result of 
    L = np.linalg.cholesky(quadr_coeff_G)
    Linv = np.linalg.inv(L)

    ##   Need to adopt the conventions used for directly computing
    ## the values of ``z`` and ``r``, we need to hold onto 
    ## the values of the factorization QR = B = L^{-1} N. 
    Q = None
    R = None
    J1 = None
    J2 = None

    ## indices of constraints being considered 
    active_set = np.array([], dtype=(np.dtype(int))) 

    ## normal vector of a given constraint
    n_p = None
    ## index of n_p in the choice of equality or inequality constraint
    p = None

    ## number of considered constraints in active set
    q = 0

    u = None
    lagr = None

    # index out of all constraints which was dropped
    k_dropped = None
    # index of active constraint which was dropped
    j_dropped = None

    z = None ## Step direction in `primal' space
    r = 0    ## Step direction in `dual' space

    ###-------------------------------------###

    ###~~~~~~~~ Step 1 ~~~~~~~~###
    while not DONE:
        ineq = np.ravel((ineq_coef_C.T * sol) - ineq_vec_d)

        if ADDING_EQ_CONSTRAINTS:
            eq_prb = np.ravel((eq_coef_A.T * sol) - eq_vec_b)

        if (len(active_set) == m_eq):
                    ADDING_EQ_CONSTRAINTS = False

        if (np.any(ineq < 0) or ADDING_EQ_CONSTRAINTS):
            if (ADDING_EQ_CONSTRAINTS):
                p = len(active_set)
                n_p = eq_coef_A[:,p]
            else:
                # Choose a violated constraint not in active set.
                #  This is the most naive way, can be improved. 
                violated_constraints = np.ravel(np.where(ineq < 0))
                v = [x for x in violated_constraints if x not in active_set]

                # Pick the first violated constraint.
                p = v[0]

                # normal vector for each constraint, vector normal to the plane.
                n_p = ineq_coef_C[:,p]

            if q == 0:
                u = 0
            
            lagr = np.hstack((u, 0))

            ###~~~~~~~~ Step 2 ~~~~~~~~###
            FULL_STEP = False

            while not FULL_STEP:
                # algo as writ will cycle back here after taking a step
                # in dual space, update inequality portion
                ineq = np.ravel((ineq_coef_C.T * sol) - ineq_vec_d) 
                if ADDING_EQ_CONSTRAINTS:
                    eq_prb = np.ravel((eq_coef_A.T * sol) - eq_vec_b)

                ###~~~~~~~~ Step 2(a) ~~~~~~~~###
                ## Calculate step directions
                if first_pass:
                    z = np.linalg.inv(quadr_coeff_G) * n_p
                    first_pass = False
                else:
                    # step direction in the primal space
                    z = J2 * J2.T * n_p

                    if (q > 0):
                        # negative of step direction in the dual space
                        # r will have num_rows = len(active_set)
                        r = np.linalg.inv(R[0:q, 0:q]) * J1.T * n_p

                ###~~~~~~~~ Step 2(b) ~~~~~~~~###
                # partial step length t1 - max step in dual space
                if ((q == 0) or (r <= 0) or ADDING_EQ_CONSTRAINTS):
                    t1 = np.inf
                else:
                    t1 = np.inf
                    k_dropped = None

                    for j in range(m_eq+1, len(active_set)):
                        k = active_set[j]
                        if (r[j] > 0) and (lagr[j] / r[j]) < t1:
                            t1 = lagr[j]/r[j]
                            k_dropped = k
                            j_dropped = j

                    t1 = np.ravel(t1)[0]


                # full step length t2 - min step in primal space
                if (np.all(z == 0)):
                    # If no step in primal space
                    t2 = np.inf
                else:
                    if ADDING_EQ_CONSTRAINTS:
                        t2 = (-1) * eq_prb[p] / (z.T * n_p)
                    else:
                        t2 = (-1) * ineq[p] / (z.T * n_p)
                    
                    t2 = np.ravel(t2)[0]

                # current step length
                t = np.min([t1, t2])


                ###~~~~~~~~ Step 2(c) ~~~~~~~~###
                if(t == np.inf):
                    print("infeasible! Stop here!")
                    FULL_STEP = True
                    DONE = True
                    return None
                    #break


                # If t2 is infinite, then we took a partial step in the dual space.
                if(t2 == np.inf):
                    #print("t2 infinite")
                    # Update lagrangian
                    lagr = lagr + t * np.hstack((np.ravel(-1 * r), 1))

                    # Drop the constraint which minimized the step we took at that
                    # point.
                    active_set = np.delete(active_set, j_dropped)

                    q = q - 1

                    Q,R = qr_delete(Q, R, j_dropped, 1, 'col')
                    J1 = Linv.T * Q[:,[x for x in range(0, q)]]
                    J2 = Linv.T * Q[:,[x for x in range(q, Q.shape[1])]]

                    # go back to step 2(a)
                    continue


                # Update iterate for x, and the lagrangian
                sol = sol + t * z
                lagr = lagr + t * np.hstack((np.ravel(-1 * r), 1))

                # if we took a full step
                if (t == t2):
                    #print("full step")
                    ## Add the constraint to the active set
                    active_set = np.hstack((active_set, p))
                    q = q + 1
                    u = lagr[-q:]

                    if Q is None:
                        Q,R = np.linalg.qr(Linv * n_p, mode="complete")
                        J1 = Linv.T * Q[:,[x for x in range(0, q)]]
                        J2 = Linv.T * Q[:,[x for x in range(q, Q.shape[1])]]
                    else:
                        Q,R = qr_insert(Q, R, (Linv * n_p) ,len(active_set)-1, 'col')
                        J1 = Linv.T * Q[:,[x for x in range(0, q)]]
                        J2 = Linv.T * Q[:,[x for x in range(q, Q.shape[1])]]

                    # Exit current loop for Step 2, go back to Step 1
                    FULL_STEP = True
                    break

                # if we took a partial step
                if (t == t1):
                    #print("partial step")
                    # Drop constraint k
                    active_set = np.delete(active_set, j_dropped)
                    q = q - 1

                    Q,R = qr_delete(Q, R, j_dropped, 1, 'col')
                    J1 = Linv.T * Q[:,[x for x in range(0, q)]]
                    J2 = Linv.T * Q[:,[x for x in range(q, Q.shape[1])]]

                    # Go back to step 2(a)
                    continue

        else:
            DONE = True

    #print(ineq)
    #print(sol)
    return sol


if __name__ == "__main__":
    ## Test 00
    #### https://tinyurl.com/ux8x4dv
    print("Running Test 00 :: https://tinyurl.com/ux8x4dv")
    nconstraint = 3

    G = np.matrix([[4, -2],
                   [-2, 4]])

    a = np.matrix([[6], [0]])

    C = np.matrix([[1, 0, 1], 
                   [0, 1, 1]])

    d = np.matrix([[0],[0],[2]])

    m_eq = 0

    truth = np.matrix([[0.5],
                       [1.5]])

    est = quadprog_solve(G, a, 
                         nconstraint, C, d,
                         m_eq, None, None)

    if (np.allclose(truth, est)):
        print("Test 00 successful!")
    else:
        print("test failed... whoops.")


    ## Test 01
    #### https://tinyurl.com/soytugm
    print("Running Test 01 :: https://tinyurl.com/soytugm")
    nconstraint = 2
    C = np.matrix([[1, 0], 
                   [0, 1]])
    d = np.matrix([[1],[0]])

    A = np.matrix([[1],
                   [1]])
    b = np.matrix([[6]])

    m_eq = 1

    truth = np.matrix([[5/2],
                       [7/2]])

    est = quadprog_solve(G, a,
                         nconstraint, C, d,
                         m_eq, A, b)

    if (np.allclose(truth, est)):
        print("Test 01 successful!")
    else:
        print("test failed... whoops.")
