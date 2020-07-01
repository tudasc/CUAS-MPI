#!/usr/bin/env python3
"""
Relative sea level forcing format is e.g.:
    netcdf GPS13Vertical_Displacement {
      dimensions:
        time = 114512 ;
      variables:
    	 double time(time) ;
    	   time:units = "seconds since 2017-07-05 17:07:57" ;
    	 float height(time) ;
    		height:units = "m" ;
    }
"""
import numpy as np
import scipy.sparse.linalg as spla
from timeit import default_timer as timer
import datetime
import argparse
from netCDF4 import Dataset
from tqdm import tqdm
import scipy.ndimage.morphology as morph
from scipy import ndimage
from scipy.interpolate import interp1d
# CUAS stuff:
from CUAS import __version__
from CUAS import helpers
from CUAS import timeparse
from CUAS.fill_matrix_fast import fill_matrix_coo
from CUAS.special_gradient import gradient2
import CUAS.physical_constants as pc

from CUAS import model as md

stencil = np.array([[0, 1, 0],
                    [1, 4, 1],
                    [0, 1, 0]]) / 8



parser = argparse.ArgumentParser(description='Run CUAS on an input geometry.')

parser.add_argument('netcdf', metavar='Netcdf_input_file', type=str)
parser.add_argument('output', metavar='Netcdf_output_file', type=str)

parser.add_argument('--totaltime', type=str, default='10 years',
                    help="""Total time to run model. 
                    Example: --totaltime '4 weeks', --totaltime '3 years 6 months' or --totaltime 50year.""")
parser.add_argument('--dt', type=str, default='12 hours',
                    help="""Time step length. 
                    Example: --dt '12 hours', --dt 12hours or --dt 1day.""")
parser.add_argument('--saveEvery', type=int, default=100, help="Save every nth timestep to netcdf.")
parser.add_argument('--conductivity', metavar='K', type=float, default=10, help="Conductivity of layer.")
parser.add_argument('--dochannels', action='store_true', help="Evolve channels?")
parser.add_argument('--disableUnconfined', action='store_true', help="Disable onconfined aquifer case.")
parser.add_argument('--Tmax', type=float, default=20.0, help="Maximum T to be allowed in the evolution")
parser.add_argument('--Tmin', type=float, default=0.0000001, help="Minimum T to be allowed in the evolution")
parser.add_argument('--flowConstant', type=float, default=5e-25, help="Ice Flow Constant A")
parser.add_argument('--roughnessFactor', type=float, default=1.0, help="Roughness factor for opening term")
parser.add_argument('--supplyMultiplier', type=float, default=1, help="Multiplier for supply")
parser.add_argument('--layerThickness', type=float, default=0.1, help="Water layer thickness (m)")
parser.add_argument('--unconfSmooth', type=float, default=0.0, help="Unconfined confined transition (m)")
parser.add_argument('--restart', type=str, default=None, help="Restart from this file.")
parser.add_argument('--Ssmulti', type=float, default=1.0, help="Multiplier for specific storage Ss.")
parser.add_argument('--Sy', type=float, default=0.4, help="Specific yield for unconfined.")
parser.add_argument('--Texp', type=float, default=1, help="Exponent of T.")
parser.add_argument('--noSmoothMelt', action='store_true',
                    help="Smooth melt term before computing change in T?")
parser.add_argument('--loopForcing', action='store_true',
                    help="""Loop the forcing when total time is longer than forcing. 
                    Otherwise the last step of the forcing is used""")
parser.add_argument('--basalVelocityIce', type=float, default=1e-6,
                    help="Basal velocity of the ice (m/s)")
parser.add_argument('--cavityBeta', type=float, default=5e-4,
                    help="cavity opening parameter")
parser.add_argument('--initialHead', default='Nzero', choices=['zero', 'Nzero', 'topg'],
                    help="""Initial value for head. 
                    Nzero means, that head is set such, that the effective pressure N is zero.""" )
parser.add_argument('--tempResults', type=str, default=None,
                    help="Save temporary results to this file(s) to later restart from them.")
parser.add_argument('--version', action='version',
                    version='%(prog)s {version}'.format(version=__version__))
parser.add_argument('--seaLevelForcing', type=str, default=None,
                    help="Apply sea level forcing from netcdf scalar time series file.")
parser.add_argument('-v', '--verbose', action='store_true',
                    help="Verbose output. Disables Progressbar")

def main():
    ##############################################
    # config stuff
    ##############################################
    args = parser.parse_args()
    verbose = args.verbose
    Tmax = args.Tmax
    Tmin = args.Tmin
    outfile = args.output
    save_every = args.saveEvery
    if args.disableUnconfined:
        enable_unconfined = False
    else:
        enable_unconfined = True
    ########################################
    totaltime_secs = timeparse.timeparse(args.totaltime)
    time_scaling, time_units = helpers.time_display_units(totaltime_secs)
    dt_secs = timeparse.timeparse(args.dt)
    Nt = int(np.rint(totaltime_secs / dt_secs))
    ########################################
    if args.netcdf == "NoData":
        geom = md.geometry.no_data_opt()
    else:
        geom = md.geometry.from_netcdf(args.netcdf)
    ########################################

    bt = args.layerThickness
    Sp = np.ones((geom.Ny, geom.Nx)) * 0
    K = np.ones((geom.Ny, geom.Nx)) * args.conductivity
    T = K*bt
    T[:,:] = args.Tmin  # initial T
    T[:,:] = 0.2
    # T = np.random.rand(Nx, Ny, dtype=np.float32) * 1
    S = np.ones((geom.Ny, geom.Nx)) * pc.Ss * bt * args.Ssmulti

    # neumann bnd conditions as very low T/K
    noflow_mask = geom.bnd_mask == pc.NOFLOW_FLAG
    K[noflow_mask] = pc.NOFLOW_VALUE
    T[noflow_mask] = pc.NOFLOW_VALUE
    T_n = T.copy()

    Q = geom.Q / pc.SPY * args.supplyMultiplier  # in seconds again
    Q = Q.astype(np.float64)

    ##################################
    # boundary points
    ##################################
    dirichlet_mask = geom.bnd_mask == pc.DIRICHLET_FLAG
    dirichlet_mask = np.logical_or(dirichlet_mask, geom.bnd_mask == pc.DIRICHLET_LAKE_FLAG)
    dirichlet_mask = np.logical_or(dirichlet_mask, noflow_mask)
    dirichlet_values = helpers.pressure2head(geom.p_ice, geom.topg)
    dirichlet_values = dirichlet_values.astype(np.float64)

    sea_level_forcing_mask = geom.bnd_mask == pc.DIRICHLET_FLAG # connected to ocean
    # np.savetxt('data.csv', sea_level_forcing_mask, delimiter=',')

    ###################################
    # mask for gradient
    grad_mask = morph.binary_dilation(noflow_mask)

    # time dependent forcing?
    if len(Q.shape) == 3:
        interpolated_forcing = interp1d(geom.time_forcing, Q, axis=0)

        if args.loopForcing:
            def get_current_Q(time):
                time_looped = time % geom.time_forcing[-1]
                return interpolated_forcing(time_looped)
        else:
            def get_current_Q(time):
                if time > geom.time_forcing[-1]:
                    return interpolated_forcing(geom.time_forcing[-1])
                else:
                    return interpolated_forcing(time)
    else:
        def get_current_Q(time):
            return Q

    # time dependent tidal (scalar) forcing
    if args.seaLevelForcing:
        ## read data file
        if verbose:
            print("Read sea level forcing from file <%s> ..."%args.seaLevelForcing)
        ncin = Dataset(args.seaLevelForcing, 'r', format='NETCDF4')
        time_rsl_var = ncin.variables["time"]
        height_rsl_var = ncin.variables["height"]
        time_rsl = time_rsl_var[:]
        height_rsl = height_rsl_var[:]
        ncin.close()
        if verbose:
            print("... Done (found nTimes = %d within [%f, %f] %s)." %
                  (len(time_rsl),time_rsl[0]/time_scaling, time_rsl[-1]/time_scaling, time_units))

        # This is the 1d (time) interpolation context.
        # No errors if time is out of bounds (return zero as fill value in those cases)
        rsl_forcing_func = interp1d(time_rsl, height_rsl, bounds_error=False, fill_value={0.0, 0.0})

        def get_current_rsl(time):
            return rsl_forcing_func(time)

    else:
        def get_current_rsl(time):
            return 0.0




    ##################################
    # SOLVER PREPARATION
    ##################################
    It = np.arange(Nt+1, dtype=np.int)
    u = np.zeros((geom.Ny, geom.Nx))   # unknown u at new time level
    u_n = np.zeros((geom.Ny, geom.Nx))   # u at the previous time level
    theta = 1  # 1 means fully implicit, 0 means fully explicit, 0.5 is Crank-Nicholson
    b = np.zeros(geom.Ny*geom.Nx)

    # initial head
    def select_initial_head(choice, p_ice, topg):
        return {
            'zero': 0,
            'Nzero': helpers.pressure2head(geom.p_ice, topg),
            'topg': topg,
        }[choice]

    u_n[:,:] = select_initial_head(args.initialHead, geom.p_ice, geom.topg)

    # restart from previous file
    if args.restart:
        if verbose:
            print("Read restart from file <%s> ..."%args.restart)
        data = Dataset(args.restart, mode='r')
        u_n[:,:] = data.variables['head'][:]
        T = data.variables['transmissivity'][:]
        T_n = T.copy()
        u_n = u_n.astype(np.float64)
        T_n = T_n.astype(np.float64)
        T = T.astype(np.float64)

        # check params in restart
        print("CUAS parameters used in restart file:")
        for key, value in vars(args).items():
            if key in ['conductivity', 'Tmax', 'Tmin', 'flowConstant', 'roughnessFactor'
                       'supplyMultiplier', 'layerThickness', 'unconfSmooth', 'Ssmulti',
                       'Sy', 'Texp', 'basalVelocityIce', 'cavityBeta']:
                print("\t", key, value)

        data.close()

    ##################################
    # CAVITY OPENING
    ##################################

    def cavity_open_b(beta, v_b, K):
        """From Werder2013 / summers2018:
        beta = (b r − b)/l r for b < b r , beta = 0 for b ≥ b r
        """
        # if limit:
        #     beta = (b_r*K - b*K) / (l_r*K)
        #     beta[b >= b_r*K] = 0
        #     return beta * u_b
        # elif multiK:
        #     beta = b_r/l_r
        #     return beta * u_b * K
        # else:
        return beta * v_b * K

    ################################################
    # like deFleurian2016
    def compute_melt(r, g, rho_w, T, K, gradh2, rho_i, L, bt):
        """avoid gradh**2 and use gradh2 = gradh**2 directly from gradient2()"""
        return r * g * rho_w * T*K * gradh2 / (rho_i * L)




    ##################################
    # save timedependent values
    ##################################
    Ntsaved, r = np.divmod(It.size, save_every)
    if r > 0: Ntsaved +=1

    if verbose:
        print("runtime = %f, time step = %f, Ntsaved = %d for save_every = %d"
              %(totaltime_secs, dt_secs, Ntsaved, save_every))

    solution = md.solution(geom.Nx, geom.Ny, Ntsaved)
    # set initial conditions as first timestep
    solution.setInitialConditions(head=u_n, T=T_n, melt=0, creep=0, Q=0)
    ##################################
    time_current = 0
    start = timer()
    for n in tqdm(It[1:], disable=verbose):
        time_current = time_current + dt_secs # todo tkleiner: avoid accumulation of floating point errors

        current_Q = get_current_Q(time_current)

        if args.seaLevelForcing:
            rsl = get_current_rsl(time_current)
            head_with_SL = helpers.pressure2head(geom.p_ice, geom.topg, sea_level=rsl).astype(np.float64)
            dirichlet_values = np.where(sea_level_forcing_mask, head_with_SL, dirichlet_values)

        Teff = T.copy()
        if enable_unconfined:
            Psi = u_n - geom.topg
            Psi[Psi < 0] = 0.01
            Teff = np.where(Psi < bt, K*Psi, T_n)

            Sp[:,:] = 0
            # TODO use specific YIELD!!
            Sp = np.where(Psi < (bt-args.unconfSmooth), 0.4, 0)
            # Sp = np.where(Psi>=(bt-args.unconfSmooth and Psi<bt), 0.4/args.unconfSmooth*bt-Psi, 0)
        else:
            Sp[:,:] = 0

        Se = S + Sp # effective storage

        # cython can not directly handle bools, so give the mask as uint8:
        A, b = fill_matrix_coo(geom.Ny,
                               geom.Nx,
                               Se,
                               Teff**args.Texp,
                               geom.dx,
                               dt_secs,
                               theta,
                               u_n,
                               current_Q,
                               dirichlet_values,
                               dirichlet_mask.astype(np.uint8))
        A = A.tocsc()
        # direct solver:
        c = spla.spsolve(A, b)
        u = c.reshape(geom.Ny, geom.Nx, order='F')

        # ensure h >= topg:
        # u = np.maximum(u, topg)
        # u = np.where(u_n-topg<0, topg, u)

        p_w = (u_n-geom.topg - bt) * pc.RHO_WATER * pc.GRAVITY #TODO: is the use of u_n (last solution) intendet here?
        N = geom.p_ice - p_w
        if args.dochannels:
            # grady, gradx = np.gradient(u_n)
            # maggrad = np.sqrt(gradx**2 + grady**2) / dx
            # maggrad2 = maggrad**2
            maggrad2 = gradient2(geom.Ny, geom.Nx, u_n, geom.dx)
            maggrad2[grad_mask] = 0

            creep = 2 * args.flowConstant * (N/3)**3 * T

            melt = compute_melt(args.roughnessFactor,
                                pc.GRAVITY,
                                pc.RHO_WATER,
                                T_n**args.Texp,
                                K, maggrad2,
                                pc.RHO_ICE,
                                pc.LATENT_HEAT,
                                bt)
            if not args.noSmoothMelt:
                melt = ndimage.convolve(melt, stencil, mode='nearest')
            cavity_opening = cavity_open_b(args.cavityBeta,
                                           args.basalVelocityIce,
                                           K )
            T = T_n + (melt + cavity_opening - creep) * dt_secs
            np.clip(T, Tmin, Tmax, out=T)
            T[noflow_mask] = pc.NOFLOW_VALUE
        else:
            melt = 0
            cavity_opening = 0
            creep = 0

        # update u_n before next step
        u_n, u = u, u_n
        T_n, T = T, T_n
        # print(n,'/',It[-1]) ## now handled by tqdm
        #
        if n % save_every == 0:
            if verbose: print("time(%06d/%06d) = %f (%s)" %(n, Nt, time_current/time_scaling, time_units))
            eps_head = np.max(np.abs(u - u_n))
            eps_T = np.max(np.abs(T - T_n))
            solution.saveTimestep(
                    time_current,
                    head=u,
                    T=T,
                    melt=melt+cavity_opening,
                    creep=creep,
                    Q=current_Q,
                    eps_head=eps_head,
                    eps_T=eps_T)

    end = timer()
    print('computation took: '+str(datetime.timedelta(seconds=end-start)))
    print('Tmax:', Tmax)
    print('Tmin:', Tmin)
    print('b:', bt)

    ###############################################
    ## save to netcdf
    ###############################################
    config = vars(args)
    config['cputime'] = end-start
    config['cputime_human'] = str(datetime.timedelta(seconds=end-start))
    config['dx'] = geom.dx

    solution.saveAsNetcdf(outfile, geom, geom.bnd_mask, config)

