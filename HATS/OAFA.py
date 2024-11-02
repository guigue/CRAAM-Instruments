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
    # OAFA coordinates
    # From https://en.wikipedia.org/wiki/F%C3%A9lix_Aguilar_Observatory 2023-04-18T1501BRT
    lon = -69.3265 * u.degree
    lat = -31.8023 * u.degree
    height = 2.420e+03 * u.m
    return EarthLocation.from_geodetic(lon,lat,height)

