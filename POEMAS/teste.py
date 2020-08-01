import timeit
import poemas
import time

start_time = time.time()

trkf = poemas.openP("TestData/SunTrack_120127_120001.TRK", ms=1)
trkf2 = poemas.openP("TestData/SunTrack_120127_120001.TRK")


trkf.to_fits(name="SunTrack_120127_120001_Integration.fits")
trkf2.to_fits(name="SunTrack_120127_120001_Subintegration.fits")

end_time = time.time()
total = end_time-start_time

print(total)