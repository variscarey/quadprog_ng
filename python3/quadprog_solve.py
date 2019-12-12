#!/usr/bin/python3
import numpy as np

###-------------------------------------###
def quadprog_solve(G, a, nconstraint, C, b):
    DONE = False
    FULL_STEP = False

    ###-------------------------------------###
    sol = (-1) * np.linalg.inv(G) * a  # Solution iterate 

    ## We only need to keep ahold of L^{-1} for the implementation. 
    ## L is the triangular result of 
    L = np.linalg.cholesky(G)
    Linv = np.linalg.inv(L)

    ##   Need to adopt the conventions used for directly computing
    ## the values of ``z`` and ``r``, we need to hold onto 
    ## the values of the factorization QR = B = L^{-1} N. 
    Q = None
    R = None
    J1 = None
    J2 = None

    first_pass = True

    ## indices of constraints being considered 
    active_set = np.array([], dtype=(np.dtype(int))) 

    ## number of considered constraints in active set
    q = 0

    u = None
    lagr = None

    k_dropped = None

    z = None ## Step direction in `primal' space
    r = 0    ## Step direction in `dual' space

    ###-------------------------------------###

    ###~~~~~~~~ Step 1 ~~~~~~~~###
    while not DONE:
        ineq = np.ravel((C.T * sol) - b)

        if (np.any(ineq < 0)):
            # Choose a violated constraint not in active set.
            #  This is the most naive way, can be improved. 
            violated_constraints = np.ravel(np.where(ineq < 0))
            v = [x for x in violated_constraints if x not in active_set]

            # Pick the first violated constraint.
            p = v[0]

            # normal vector for each constraint, vector normal to the plane.
            n_p = C[:,p]

            if q == 0:
                u = 0
            
            lagr = np.hstack((u, 0))

            ###~~~~~~~~ Step 2 ~~~~~~~~###
            FULL_STEP = False

            while not FULL_STEP:
                # algo as writ will cycle back here after taking a step
                # in dual space, update inequality portion
                ineq = np.ravel((C.T * sol) - b) 

                ###~~~~~~~~ Step 2(a) ~~~~~~~~###
                ## Calculate step directions
                if first_pass:
                    z = np.linalg.inv(G) * n_p
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
                if ((q == 0) or (r <= 0)):
                    t1 = np.inf
                else:
                    t1 = np.inf
                    k_dropped = None

                    for j in range(0, len(active_set)):
                        k = active_set[j]
                        if (r[j] > 0) and (lagr[j] / r[j]) < t1:
                            t1 = lagr[j]/r[j]
                            k_dropped = k

                    t1 = np.ravel(t1)[0]


                # full step length t2 - min step in primal space
                if (np.all(z == 0)):
                    # If no step in primal space
                    t2 = np.inf
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
                    active_set = np.delete(active_set,
                                           np.where(active_set == k_dropped))

                    q = q - 1

                    # Update H and N*
                    N = np.matrix(C[:,active_set]) 
                    B = Linv * N

                    Q,R = np.linalg.qr(B, mode='complete')
                    Q1 = Q[:,[x for x in range(0, q)]]
                    Q2 = Q[:,[x for x in range(q, Q.shape[1])]]
                    J1 = Linv.T * Q1
                    J2 = Linv.T * Q2

                    Rsquare = R[0:q, 0:q]

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

                    # Update Q,R,J's,etc. 
                    N = np.matrix(C[:,active_set]) 
                    B = Linv * N

                    Q,R = np.linalg.qr(B, mode='complete')
                    Q1 = Q[:,[x for x in range(0, q)]]
                    Q2 = Q[:,[x for x in range(q, Q.shape[1])]]
                    J1 = Linv.T * Q1
                    J2 = Linv.T * Q2

                    Rsquare = R[0:q, 0:q]

                    # Exit current loop for Step 2, go back to Step 1
                    FULL_STEP = True
                    break

                # if we took a partial step
                if (t == t1):
                    #print("partial step")
                    # Drop constraint k
                    active_set = np.delete(active_set,
                                           np.where(active_set == k_dropped))
                    q = q - 1

                    # Update H and N*
                    N = np.matrix(C[:,active_set]) 
                    B = Linv * N

                    Q,R = np.linalg.qr(B, mode='complete')
                    Q1 = Q[:,[x for x in range(0, q)]]
                    Q2 = Q[:,[x for x in range(q, Q.shape[1])]]
                    J1 = Linv.T * Q1
                    J2 = Linv.T * Q2

                    # Go back to step 2(a)
                    continue

        else:
            DONE = True

    #print(ineq)
    #print(sol)
    return sol


if __name__ == "__main__":
    nconstraint = 3

    G = np.matrix([[4, -2],
                   [-2, 4]])

    C = np.matrix([[1, 0, 1], 
                   [0, 1, 1]])

    b = np.matrix([[0],[0],[2]])

    a = np.matrix([[6], [0]])

    truth = np.matrix([[0.5],
                       [1.5]])

    est = quadprog_solve(G, a, nconstraint, C, b)

    if (np.allclose(truth, est)):
        print("Test successful!")