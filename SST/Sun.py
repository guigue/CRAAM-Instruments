class Sun(object):

    def __init__(self):
        self.quiet_sun_temp = pd.dataFrame({'212': 5.9E+03, '405':5.1E+03})

    def Get_Sun_Coordinates(self,ctime):
        #
        #  Get_Sun_Coordinates
        #
        #     This method obtains the Sun Coordinates both in (Ra,Dec) and (Alt, Az)
        #     While astropy may use very accurate calculations and corrections from
        #     IERS, for some reason it cannot download the data. So, I decided to
        #     work off-line with a lesser precision, although enoug for us.
        #
        #  @guiguesp - 2020-04-05
        #
        #_____________________________________________________________________________


        # CASLEO coordinates
            lon = -69.29669444 * u.degree
            lat = -31.79897222 * u.degree
            height = 2.550e+03 * u.m
            cg=EarthLocation.from_geodetic(lon,lat,height)

            t = Time(self.ctime)
            #
            # This is needed because IERS data repository is not
            # accepting downloads. The precision is lower, but
            # it' enough for our computations
            #
            iers.conf.auto_download = False

            with solar_system_ephemeris.set('builtin'):
                sun_radec=get_body('sun',t,cg)

            aa = AltAz(location=cg,obstime=t)
            sun_altaz = sun_radec.transform_to(aa)

            return sun_radec,sun_altaz,cg
