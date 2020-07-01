"""Helper scripts for CUAS
"""
import CUAS.physical_constants as pc
from CUAS import timeparse


def time_display_units(totaltime_secs):
    """Convert time units from seconds to hours, days or years
    depending on the total run length"""
    scaling = 1
    units = 'seconds'
    if totaltime_secs > timeparse.timeparse('10 hours'):
        scaling = timeparse.timeparse('1 hour')
        units = 'hours'
    if totaltime_secs > timeparse.timeparse('10 days'):
        scaling = timeparse.timeparse('1 day')
        units = 'days'
    if totaltime_secs > timeparse.timeparse('10 years'):
        scaling = timeparse.timeparse('1 year')
        units = 'years'

    return scaling, units


def head2pressure(head, bed_elevation, sea_level=0.0):
    """Converts hydraulic head to water pressure
    bed_elevation and sea_level are relative to the
    mean sea level = 0.0
    """
    effective_bed_elevation = bed_elevation - sea_level
    return pc.RHO_WATER * pc.GRAVITY * (head - effective_bed_elevation)


def overburdenPressure(thk):
    """Compute ice overburden pressure
    """
    return thk * pc.RHO_ICE * pc.GRAVITY


def pressure2head(pressure, bed_elevation, sea_level=0.0):
    """Convert water pressure to hydraulic head
    bed_elevation and sea_level are relative to the
    mean sea level = 0.0
    pressure is ice pressure
    """
    effective_bed_elevation = bed_elevation - sea_level
    return pressure / (pc.RHO_WATER * pc.GRAVITY) + effective_bed_elevation
