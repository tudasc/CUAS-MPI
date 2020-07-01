import numpy as np
from scipy.sparse import coo_matrix
cimport numpy as np
cimport cython

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t


cdef DTYPE_t TINY = np.finfo(np.float32).eps


@cython.cdivision # dont check for division through 0
cdef inline DTYPE_t hmean(DTYPE_t x1, DTYPE_t x2):
    return 2.0 * x1 * x2 / (x1 + x2 + TINY)

cdef inline int m(int i, int j, int Nx):
    return j*Nx + i

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.cdivision # dont check for division through 0
def fill_matrix_coo(int Nx, int Ny,
        np.ndarray[DTYPE_t, ndim=2] S,
        np.ndarray[DTYPE_t, ndim=2] T,
        DTYPE_t dx,
        DTYPE_t dt,
        DTYPE_t theta,
        np.ndarray[DTYPE_t, ndim=2] u_n,
        np.ndarray[DTYPE_t, ndim=2] Q,
        np.ndarray[DTYPE_t, ndim=2] dirichlet_values,
        np.ndarray[np.uint8_t, ndim=2] dirichlet_mask,
        ):
    cdef int N = Nx*Ny

    cdef Py_ssize_t i, j, k, p
    cdef DTYPE_t S_P, d_N, d_S, d_W, d_E, d_P
    cdef DTYPE_t A_N, A_S, A_W, A_E, A_P

    # cdef np.ndarray[int, ndim=1] rows = np.zeros(N*5, dtype=np.int32)
    # cdef np.ndarray[int, ndim=1] cols = np.zeros(N*5, dtype=np.int32)
    # cdef np.ndarray[DTYPE_t, ndim=1] data = np.zeros(N*5, dtype=DTYPE)
    Ndirichlet = np.count_nonzero(dirichlet_mask)
    Ndirichlet = Ndirichlet + 2*Nx + 2*Ny ## this is only necessary, because we dont do the boundary points at the edge of the grid double...
    cdef np.ndarray[int, ndim=1] rows = np.zeros(N*5+Ndirichlet, dtype=np.int32)
    cdef np.ndarray[int, ndim=1] cols = np.zeros(N*5+Ndirichlet, dtype=np.int32)
    cdef np.ndarray[DTYPE_t, ndim=1] data = np.zeros(N*5+Ndirichlet, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] b = np.zeros(N, dtype=DTYPE)

    #############################


    k = 0

    for j in range(Ny):
        for i in range(Nx):
            if dirichlet_mask[i,j]:
                # check if it is a dirichlet point and if yes, then set the
                # diagonal entry to 1 and the RHS to the corresponding dirichlet
                # value.
                p = m(i,j, Nx)
                rows[k] = p
                cols[k] = p
                data[k] = 1
                b[p] = dirichlet_values[i,j]
                k = k+1

            else:
                S_P = S[i,j]

                d_N = hmean(T[i,j], T[i,j+1]) / dx**2
                d_S = hmean(T[i,j], T[i,j-1]) / dx**2
                d_W = hmean(T[i,j], T[i-1,j]) / dx**2
                d_E = hmean(T[i,j], T[i+1,j]) / dx**2
                d_P = -(d_N + d_S + d_W + d_E)

                A_N = - theta * dt/S_P * d_N
                A_S = - theta * dt/S_P * d_S
                A_W = - theta * dt/S_P * d_W
                A_E = - theta * dt/S_P * d_E
                A_P = 1 - theta * dt/S_P * d_P

                p = m(i,j, Nx)
                ## fill A matrix
                # A[p, m(i,j-1)] = A_S
                rows[k] = p
                cols[k] = m(i,j-1, Nx)
                data[k] = A_S
                k = k+1

                # A[p, m(i-1,j)] = A_W
                rows[k] = p
                cols[k] = m(i-1,j, Nx)
                data[k] = A_W
                k = k+1
                # A[p, p]        = A_P
                rows[k] = p
                cols[k] = p
                data[k] = A_P
                k = k+1
                # A[p, m(i+1,j)] = A_E
                rows[k] = p
                cols[k] = m(i+1,j, Nx)
                data[k] = A_E
                k = k+1
                # A[p, m(i,j+1)] = A_N
                rows[k] = p
                cols[k] = m(i,j+1, Nx)
                data[k] = A_N
                k = k+1

                ## fill b
                b[p] = u_n[i,j] + \
                  (1-theta)* \
                  dt/S_P * (d_N*u_n[i,j+1] + d_S*u_n[i,j-1] +\
                  d_P*u_n[i,j] + d_W*u_n[i-1,j] +\
                  d_E*u_n[i+1,j])+\
                  dt/S_P * Q[i,j]
    return coo_matrix((data, (rows, cols)), shape=(Nx*Ny, Nx*Ny)), b
