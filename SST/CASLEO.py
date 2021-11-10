from astropy.coordinates import solar_system_ephemeris, EarthLocation
from astropy import units as u
from astropy.utils import iers

def Observatory_Coordinates():
    #
    # This is needed because IERS data repository is not
    # accepting downloads. The precision is lower, but
    # it' enough for our computations
    #
    iers.conf.auto_download = False
    # CASLEO coordinates
    lon = -69.29669444 * u.degree
    lat = -31.79897222 * u.degree
    height = 2.491e+03 * u.m
    return EarthLocation.from_geodetic(lon,lat,height)

