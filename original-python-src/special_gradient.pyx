import numpy as np
from scipy.sparse import coo_matrix
cimport numpy as np
cimport cython

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t


cdef DTYPE_t TINY = np.finfo(np.float32).eps



@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.cdivision # dont check for division through 0
def gradient2(int Nx, int Ny,
        np.ndarray[DTYPE_t, ndim=2] h,
        DTYPE_t dx,
        ):

    cdef Py_ssize_t i, j
    cdef DTYPE_t hgrad_x, hgrad_y
    cdef np.ndarray[DTYPE_t, ndim=2] hgrad = np.zeros((Nx, Ny), dtype=DTYPE)

    #############################

    for j in range(1, Ny-1):
        for i in range(1, Nx-1):
                        # gradient_x(i,j) = (0.5 * (((array(i,j) - array(i-1,j)) / dx)**2+((array(i+1,j) - array(i,j)) / dx)**2))
            hgrad_x = (0.5 * (((h[i,j] - h[i-1,j]) / dx)**2+((h[i+1,j] - h[i,j]) / dx)**2))
            hgrad_y = (0.5 * (((h[i,j] - h[i,j-1]) / dx)**2+((h[i,j+1] - h[i,j]) / dx)**2))
            hgrad[i,j] = hgrad_x + hgrad_y

    return hgrad

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.cdivision # dont check for division through 0
def gradient2_central(int Nx, int Ny,
        np.ndarray[DTYPE_t, ndim=2] h,
        DTYPE_t dx,
        ):

    cdef Py_ssize_t i, j
    cdef DTYPE_t hgrad_x, hgrad_y
    cdef np.ndarray[DTYPE_t, ndim=2] hgrad = np.zeros((Nx, Ny), dtype=DTYPE)

    #############################

    for j in range(1, Ny-1):
        for i in range(1, Nx-1):
            hgrad_x = (h[i+1,j] - h[i-1,j]) / (2*dx)
            hgrad_y = (h[i,j+1] - h[i,j-1]) / (2*dx)
            hgrad[i,j] = hgrad_x**2 + hgrad_y**2

    return hgrad

