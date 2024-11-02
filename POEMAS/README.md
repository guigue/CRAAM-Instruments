# POEMAS
  
Python modules for working with POEMAS data.  
* poemas - Convert POEMAS Raw Binary Data to FITS.
  
# How to use

```python
import poemas

trkf = poemas.open("TestData/SunTrack_120126_105045.TRK")

#Get info
print(trkf.type) #TRK
print(trkf.date) #2012-01-26
print(trkf.time) #10:50:45
print(trkf.filename) #SunTrack_120126_105045.TRK
print(trkf.headercolumns) # ('Code', 'NRS', 'FreqNo', 'Freq1', 'Freq2', 'BRTMin', 'BRTMAX')
print(trkf.columns) # ('sec', 'ele_ang', 'zi_ang', 'TBL_45', 'TBR_45', 'TBL_90', 'TBR_90')
print(trkf.headerdata) #[(5566441, 552, 2, 45., 90., 488.25, 620.55)]
print(trkf.data[0]) #(349267845, 10.06, 106.020004, [559.1 , 579.75, 489.75, 500.35, 559., ....])

#Write data to a FITS file
trkf.to_fits()
```

# Compatibility
Python >= 3.7
Numpy = 1.18.1
