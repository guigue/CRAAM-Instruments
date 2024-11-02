# routines_IDL_30THz

There are three routines (and must be used in the following order): 
1 - transform_fpf_to_fits_30THz.pro: required to transform the files from the format .fpf to .fits 
                                     Requires: "read_fpf.pro"
2 - obtain_flat_30THz.pro:           makes the flat necessary for the processing of the data
3 - calibration_30THz.pro:           performs the calibration 
