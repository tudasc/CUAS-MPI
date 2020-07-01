from netCDF4 import Dataset

vardict = {
    'head': {
        'ncname': 'head',
        'units': 'm',
        'stdname': 'hydraulic_head',
        'lname': 'Hydraulic head'
    },
    'topg': {
        'ncname': 'topg',
        'units': 'm',
        'stdname': 'bedrock_elevation',
        'lname': 'Bedrock elevation'
    },
    'thk': {
        'ncname': 'thk',
        'units': 'm',
        'stdname': 'ice_thickness',
        'lname': 'Ice thickness'
    },
    'noflow_mask': {
        'ncname': 'noflow_mask',
        'units': '',
        'stdname': 'noflow_mask',
        'lname': 'Noflow mask'
    },
    'pwater': {
        'ncname': 'pwater',
        'units': 'Pa',
        'stdname': 'water_pressure',
        'lname': 'Water pressure'
    },
    'pice': {
        'ncname': 'pice',
        'units': 'Pa',
        'stdname': 'ice_pressure',
        'lname': 'Ice pressure'
    },
    'peffective': {
        'ncname': 'peffective',
        'units': 'Pa',
        'stdname': 'effective_pressure',
        'lname': 'Effective pressure'
    },
    'transmissivity': {
        'ncname': 'transmissivity',
        'units': 'm^2 s^-1',
        'stdname': 'transmissivity',
        'lname': 'Layer transmissivity'
    },
    'opening': {
        'ncname': 'opening',
        'units': ' ',
        'stdname': 'opening',
        'lname': 'opening rate'
    },
    'closure': {
        'ncname': 'closure',
        'units': ' ',
        'stdname': 'closure',
        'lname': 'closure rate'
    },
    'flux': {
        'ncname': 'flux',
        'units': 'm^2/s m^3/s',
        'stdname': 'water_flux',
        'lname': 'Water flux'
    },
    'watersource': {
        'ncname': 'watersource',
        'units': 'm/s',
        'stdname': 'water_source',
        'lname': 'Water input into model'
    }
}


class ncvar:
    def __init__(self, data, name, timedep=False):
        self.data = data
        self.nc_names = vardict[name]
        self.timedep = timedep

    def __repr__(self):
        return ("data: " + str(self.data) + "\n" +
                "nc_names: " + str(self.nc_names) + "\n" +
                "time_dep: " + str(self.timedep))


def write_netcdf(filename, config, x, y, t, eps_inf, Teps_inf, *args):
    '''Write the data to netcdf file'''
    nx = x.size
    ny = y.size
    root_grp = Dataset(filename, 'w', format='NETCDF4')
    root_grp.description = 'CUAS result for basal hydrology.'
    root_grp.createDimension('x', nx)
    root_grp.createDimension('y', ny)
    root_grp.createDimension('time', None)

    # variables
    nc_x = root_grp.createVariable('x', 'f4', ('x',))
    nc_x.units = 'm'
    nc_x.axis = 'X'
    nc_x.long_name = 'X-coordinate in Cartesian system'
    nc_x.standard_name = 'projection_x_coordinate'

    nc_y = root_grp.createVariable('y', 'f4', ('y',))
    nc_y.units = 'm'
    nc_y.axis = 'Y'
    nc_y.long_name = 'Y-coordinate in Cartesian system'
    nc_y.standard_name = 'projection_y_coordinate'

    nc_time = root_grp.createVariable('time', 'f4', ('time',))
    nc_time.units = 'seconds since 01-01-01 00:00:00'  # arbitrary reference date and time. Could be changed later.
    nc_time.calendar = '365_day'  # no other option as we use timeparse.py to convert years to seconds
    nc_time.standard_name = 'time'
    nc_time.axis = "T"

    nc_eps_inf = root_grp.createVariable('eps_inf', 'f4', ('time',), fill_value=-9999.0)
    nc_eps_inf.long_name = 'rate of change in head inf-norm: eps_inf = max(|h^n-h^(n-1)|)/dt'
    nc_eps_inf.units = 'm/s'

    nc_Teps_inf = root_grp.createVariable('Teps_inf', 'f4', ('time',), fill_value=-9999.0)
    nc_Teps_inf.long_name = 'rate of change in transmissivity inf-norm'
    nc_Teps_inf.units = 'm^2 s^-2'

    nc_x[:] = x
    nc_y[:] = y
    nc_time[:] = t
    nc_eps_inf[:] = eps_inf
    nc_Teps_inf[:] = Teps_inf

    for arg in args:
        if isinstance(arg, ncvar):
            if arg.timedep:
                nc_variable = root_grp.createVariable(arg.nc_names['ncname'] + '_t', 'f4', ('time', 'y', 'x'),
                                                      zlib=True)
            else:
                nc_variable = root_grp.createVariable(arg.nc_names['ncname'], 'f4', ('y', 'x'), zlib=True)
            nc_variable.units = arg.nc_names['units']
            nc_variable.standard_name = arg.nc_names['stdname']
            nc_variable.long_name = arg.nc_names['lname']
            nc_variable[:] = arg.data
        elif isinstance(arg, tuple) and len(arg) == 2:
            # tuple of (data, name)
            nc_variable = root_grp.createVariable(arg[1], 'f4', ('y', 'x'), zlib=True)
            nc_variable[:] = arg[0]

    # write config
    for key, value in config.items():
        if value is None:
            value = 'None'
        if value is True:
            value = 'True'
        if value is False:
            value = 'False'
        root_grp.setncattr(key, value)
    # root_grp.setncatts(config)

    # write version
    from . import __version__
    root_grp.setncattr('version', __version__)

    root_grp.close()
