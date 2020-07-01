import numpy as np
from netCDF4 import Dataset
import CUAS.physical_constants as pc
import CUAS.IO as IO
from timeit import default_timer as timer
import datetime
import scipy.ndimage.morphology as morph
from scipy.interpolate import interp1d


class geometry:
    """Model geometry
    """

    def __init__(self, x, y, usurf, topg, thk, bnd_mask,
                 Q, time_forcing):
        self.x = x
        self.y = y
        self.usurf = usurf
        self.topg = topg
        self.thk = thk
        self.Q = Q
        self.bnd_mask = bnd_mask
        self.time_forcing = time_forcing

        self.dx = x[1] - x[0]
        self.Ny, self.Nx = usurf.shape

        # not sure if ice pressure should really be in geometry
        self.p_ice = self.thk * pc.RHO_ICE * pc.GRAVITY

    @classmethod
    def from_netcdf(cls, ncName):
        print("Loading netCDF...")
        start = timer()
        data = Dataset(ncName, mode='r')
        x = data.variables['x'][:]
        y = data.variables['y'][:]
        usurf = data.variables['usurf'][:]
        topg = data.variables['topg'][:]
        thk = data.variables['thk'][:]
        Q = data.variables['bmelt'][:]
        bnd_mask = data.variables['bnd_mask'][:]

        # hack: alter bnd_mask if topg is not given (e.g. old BedMachine data)
        # todo: check for fill value or missing value
        bnd_mask[((topg < -1.0e4) | (topg > 1.0e4))] = pc.NOFLOW_FLAG

        # todo: Ensure we do not have active grid points at the outer margin
        #       in better python. This is more a hack.
        edge = bnd_mask[:, 0]
        edge[edge == pc.COMPUTE_FLAG] = pc.NOFLOW_FLAG
        bnd_mask[:, 0] = edge

        edge = bnd_mask[:, -1]
        edge[edge == pc.COMPUTE_FLAG] = pc.NOFLOW_FLAG
        bnd_mask[:, -1] = edge

        edge = bnd_mask[0, :]
        edge[edge == pc.COMPUTE_FLAG] = pc.NOFLOW_FLAG
        bnd_mask[0, :] = edge

        edge = bnd_mask[-1, :]
        edge[edge == pc.COMPUTE_FLAG] = pc.NOFLOW_FLAG
        bnd_mask[-1, :] = edge

        try:
            time_forcing = data.variables['time'][:]
        except KeyError:
            print('Time variable not found in netcdf.')
            print('Assuming constant forcing.')
            time_forcing = None
        data.close()
        end = timer()
        print('NetCDF loading took: ' + str(datetime.timedelta(seconds=end - start)))
        return cls(x, y, usurf, topg, thk, bnd_mask, Q, time_forcing)

    @classmethod
    def no_data_opt(cls):
        """Simple geometry without loading a file"""
        start = timer()
        x = np.arange(20) * 1000.0
        y = np.arange(10) * 1000.0
        X, Y = np.meshgrid(x, y)
        usurf = 2000.0 - 0.1 * X
        topg = np.zeros_like(usurf)
        thk = usurf - topg
        bnd_mask = np.zeros_like(usurf)
        bnd_mask[:, 0] = pc.NOFLOW_FLAG
        bnd_mask[:, -1] = pc.DIRICHLET_FLAG
        bnd_mask[0, :] = pc.NOFLOW_FLAG
        bnd_mask[-1, :] = pc.NOFLOW_FLAG
        bmelt = np.ones_like(usurf)
        time_forcing = None

        return cls(x, y, usurf, topg, thk, bnd_mask, bmelt, time_forcing)

    def saveNc(self, ncName):
        root_grp = Dataset(ncName, 'w', format='NETCDF4')
        root_grp.description = 'CUAS geometry'
        root_grp.createDimension('x', self.Nx)
        root_grp.createDimension('y', self.Ny)
        # variables
        nc_x = root_grp.createVariable('x', 'd', ('x',))
        nc_x.units = 'meters'
        nc_x.axis = 'X'
        nc_x.long_name = 'X-coordinate in Cartesian system'
        nc_x.standard_name = 'projection_x_coordinate'
        nc_y = root_grp.createVariable('y', 'd', ('y',))
        nc_y.units = 'meters'
        nc_y.axis = 'Y'
        nc_y.long_name = 'Y-coordinate in Cartesian system'
        nc_y.standard_name = 'projection_y_coordinate'

        nc_x[:] = self.x
        nc_y[:] = self.y

        nc_usurf = root_grp.createVariable('usurf', 'd', ('y', 'x'))
        nc_usurf.units = 'meters'
        nc_usurf.standard_name = 'usurf'
        nc_usurf.long_name = 'surface elevation'
        nc_usurf[:] = self.usurf
        nc_topg = root_grp.createVariable('topg', 'd', ('y', 'x'))
        nc_topg.units = 'meters'
        nc_topg.standard_name = 'topg'
        nc_topg.long_name = 'base elevation'
        nc_topg[:] = self.topg
        nc_thk = root_grp.createVariable('thk', 'd', ('y', 'x'))
        nc_thk.units = 'meters'
        nc_thk.standard_name = 'ice thickness'
        nc_thk.long_name = 'ice thickness'
        nc_thk[:] = self.thk
        nc_bmelt = root_grp.createVariable('bmelt', 'd', ('y', 'x'))
        nc_bmelt.units = 'm/yr'
        nc_bmelt.standard_name = 'basal melt'
        nc_bmelt.long_name = 'basal melt'
        nc_bmelt[:] = self.Q
        nc_bnd_mask = root_grp.createVariable('bnd_mask', 'i', ('y', 'x'))
        nc_bnd_mask.units = ''
        nc_bnd_mask.standard_name = 'boundary condition mask'
        nc_bnd_mask.long_name = 'boundary condition mask'
        nc_bnd_mask[:] = self.bnd_mask
        # nc_time_forcing = root_grp.createVariable('time', 'i', ('y', 'x'))
        # nc_time_forcing.units = ''
        # nc_time_forcing.standard_name = 'boundary condition mask'
        # nc_time_forcing.long_name = 'boundary condition mask'
        # nc_time_forcing[:] = self.time_forcing
        root_grp.close()


class boundaryConditions:
    """Boundary conditions"""

    def __init__(self, bnd_mask):
        self.noflow_mask = bnd_mask == pc.NOFLOW_FLAG
        self.dirichlet_mask = bnd_mask == pc.DIRICHLET_FLAG
        self.dirichlet_mask = np.logical_or(self.dirichlet_mask, bnd_mask == pc.DIRICHLET_LAKE_FLAG)

        # dirichlet values are set as ?
        self.dirichlet_values = 0
        self.gradient_mask = morph.binary_dilation(self.noflow_mask)


class initialConditions:
    """Initial conditions"""

    def __init__(self, head, T):
        pass

    def select_initial_head(choice, p_ice, topg):
        return {
            'zero': 0,
            # 'Nzero': helpers.pressure2head(geom.p_ice, topg),
            'topg': topg,
        }[choice]


class forcing:
    """forcing"""

    def __init__(self, geom):
        # time dependent forcing?
        self.Q = geom.Q
        self.time_forcing = geom.time_forcing
        if len(self.Q.shape) == 3:
            self.interpolated_forcing = interp1d(self.time_forcing, self.Q, axis=0)

    def get_current_Q(self, time):
        if len(self.Q.shape) == 3:
            time_looped = time % self.time_forcing[-1]
            return self.interpolated_forcing(time_looped)
        else:
            return self.Q


class config:
    """config"""

    def __init__(self, args):
        self.Tmax = args.Tmax
        self.Tmin = args.Tmin

        self.save_every = args.saveEvery
        self.enable_unconfined = ~args.disableUnconfined

        self.layerThickness = args.layerThickness
        self.conductivity = args.conductivity
        self.specific_yield = args.Sy

        # self.cavity_opening =


class solution:
    """Model solution
    """

    def __init__(self, Nx, Ny, Ntsaved):
        self.head = np.empty((Ntsaved, Ny, Nx))
        self.T = np.empty((Ntsaved, Ny, Nx))
        self.melt = np.empty((Ntsaved, Ny, Nx))
        self.creep = np.empty((Ntsaved, Ny, Nx))
        self.Q = np.empty((Ntsaved, Ny, Nx))
        self.time = np.zeros(Ntsaved)
        self.eps_head = np.empty(Ntsaved)
        self.eps_T = np.empty(Ntsaved)
        self.nsaved = 0

    def setInitialConditions(self, head, T, melt, creep, Q):
        """set initial conditions for first timestep"""
        self.head[0, :, :] = head
        self.T[0, :, :] = T
        self.melt[0, :, :] = melt
        self.creep[0, :, :] = creep
        self.Q[0, :, :] = Q
        self.eps_head[0] = -9999.0  # fill_value in IO.py
        self.eps_T[0] = -9999.0  # fill_value in IO.py
        self.nsaved = 1

    def saveTimestep(self, time_current, head, T, melt, creep, Q, eps_head, eps_T):
        """Save data for timestep"""
        i = self.nsaved
        self.head[i, :, :] = head
        self.T[i, :, :] = T
        self.melt[i, :, :] = melt
        self.creep[i, :, :] = creep
        self.Q[i, :, :] = Q
        self.time[i] = time_current
        self.eps_head[i] = eps_head
        self.eps_T[i] = eps_T
        self.nsaved = self.nsaved + 1

    def saveAsNetcdf(self, filename, geom, mask, config):
        """Save the complete solution to netcdf output"""
        pressure_water = (self.head - geom.topg) * pc.RHO_WATER * pc.GRAVITY
        pressure_effective = geom.p_ice - pressure_water

        nc_head = IO.ncvar(self.head[-1, :, :], 'head')
        nc_htime = IO.ncvar(self.head, 'head', timedep=True)
        nc_N = IO.ncvar(pressure_effective[-1, :, :], 'peffective')
        nc_Ntime = IO.ncvar(pressure_effective, 'peffective', timedep=True)
        nc_T = IO.ncvar(self.T[-1, :, :], 'transmissivity')
        nc_Ttime = IO.ncvar(self.T, 'transmissivity', timedep=True)
        nc_melttime = IO.ncvar(self.melt, 'opening', timedep=True)
        nc_creeptime = IO.ncvar(self.creep, 'closure', timedep=True)
        nc_Qtime = IO.ncvar(self.Q, 'watersource', timedep=True)

        nc_noflow_mask = IO.ncvar(mask, 'noflow_mask')
        nc_topg = IO.ncvar(geom.topg, 'topg')
        nc_thk = IO.ncvar(geom.thk, 'thk')

        # flux
        q_solution = np.empty((self.nsaved, geom.Ny, geom.Nx))
        q_solution[:] = np.nan
        for i in range(self.nsaved):
            h_area = self.head[i, 1:-1, 1:-1]  # ignore border
            T_area = self.T[i, 1:-1, 1:-1]  # ignore border
            grad_h_x = np.gradient(h_area)[1] / geom.dx
            grad_h_y = np.gradient(h_area)[0] / geom.dx
            qx = T_area * grad_h_x
            qy = T_area * grad_h_y
            q_mag = np.sqrt(qx ** 2 + qy ** 2)
            q_solution[i, 1:-1, 1:-1] = q_mag
        nc_qtime = IO.ncvar(q_solution, 'flux', timedep=True)

        IO.write_netcdf(filename, config, geom.x, geom.y, self.time, self.eps_head, self.eps_T,
                        nc_head,
                        nc_htime,
                        nc_N,
                        nc_Ntime,
                        nc_noflow_mask,
                        nc_topg,
                        nc_thk,
                        nc_T,
                        nc_Ttime,
                        nc_melttime,
                        nc_creeptime,
                        nc_Qtime,
                        nc_qtime,
                        )

