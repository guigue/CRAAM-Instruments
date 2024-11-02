#####################################################################
#
# testFITS.py
#
# This is the first phase of the conversion tests for the SST RBD
# files. It creates fits files for every different RBD format. It also
# checks come functionalities, like writing in a different directory,
# or giving a non-standard name to the fits file.
#
# Remember to set the environment variable $RBDXMLPATH where the
# xml files are stored in the system.
#
# Change the path where the actual data are, also.
#
# After running this script, you can run the IDL pro testFITS.pro
#
# Use:
#  python testFITS.py
#
# Written
#  Guigue@gcastro.net - 2017-08-31
#
##################################################################

import oRBD

d=oRBD.RBD(InputPath='~/Programming/CRAAM-Instruments/SST/TestData',
           OutputPath='~/Programming/CRAAM-Instruments/SST/TestData')
d.readRBDinDictionary('bi1010822')
d.writeFITS()

d=oRBD.RBD(InputPath='~/Programming/CRAAM-Instruments/SST/TestData',
           OutputPath='~/Programming/CRAAM-Instruments/SST/TestData')

d.readRBDinDictionary('bi1021019')
d.writeFITS()

d=oRBD.RBD(InputPath='~/Programming/CRAAM-Instruments/SST/TestData',
           OutputPath='~/Programming/CRAAM-Instruments/SST/TestData')

d.readRBDinDictionary('bi1021202')
d.writeFITS()

d=oRBD.RBD(InputPath='~/Programming/CRAAM-Instruments/SST/TestData',
           OutputPath='~/Programming/CRAAM-Instruments/SST/TestData')

d.readRBDinDictionary('bi1021221')
d.writeFITS()

d=oRBD.RBD(InputPath='~/Programming/CRAAM-Instruments/SST/TestData',
           OutputPath='~/Programming/CRAAM-Instruments/SST/TestData')

d.readRBDinDictionary('rs1020715.1300')
d.writeFITS()

d=oRBD.RBD(InputPath='~/Programming/CRAAM-Instruments/SST/TestData',
           OutputPath='~/Programming/CRAAM-Instruments/SST/TestData')

d.readRBDinDictionary('rs1021205.2200')
d.writeFITS()

d=oRBD.RBD(InputPath='~/Programming/CRAAM-Instruments/SST/TestData',
           OutputPath='~/Programming/CRAAM-Instruments/SST/TestData')

d.readRBDinDictionary('rs1061206.2100')
d.writeFITS()

d=oRBD.RBD(InputPath='~/Programming/CRAAM-Instruments/SST/TestData',
           OutputPath='~/Programming/CRAAM-Instruments/SST/TestData')

d.readRBDinDictionary('rs990909.1700')
d.writeFITS()
d.reduced()
d.writeFITSwithName('rs990909.1700.red.fits')

d1=oRBD.RBD(InputPath='~/Programming/CRAAM-Instruments/SST/TestData',
           OutputPath='~/Programming/CRAAM-Instruments/SST/TestData')
d2=oRBD.RBD(InputPath='~/Programming/CRAAM-Instruments/SST/TestData',
           OutputPath='~/Programming/CRAAM-Instruments/SST/TestData')

d1=oRBD.RBD(InputPath='~/Programming/CRAAM-Instruments/SST/TestData')
d2=oRBD.RBD(InputPath='~/Programming/CRAAM-Instruments/SST/TestData')

d1.readRBDinDictionary('rs1150621.1700')
d2.readRBDinDictionary('rs1150621.1800')

d=oRBD.RBD(InputPath='~/Programming/CRAAM-Instruments/SST/TestData',
           OutputPath='~/Programming/CRAAM-Instruments/SST/TestData')

d = oRBD.RBD(OutputPath='~/Programming/CRAAM-Instruments/SST/TestData')
d.concat([d1,d2])
d.writeFITSwithName('rs1150621.17-18-concat.fits')

d1.reduced()
d2.reduced()
d = oRBD.RBD(OutputPath='~/Programming/CRAAM-Instruments/SST/TestData')
d.concat([d1,d2])
d.writeFITSwithName('rs1150621.17-18-red-concat.fits')










