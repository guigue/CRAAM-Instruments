# Makefile for the FLIR Programs
#
# Guigue - 2017-10-24 @Sampa
#
#
########################################################################

#CFLAGS = -g -L/usr/local/lib
CFLAGS =  -g -L/usr/local/lib -I/usr/local/include
CLIBS = -lCCfits -lcfitsio
CC = g++

MAIN_PROGS = read_fpf fpf2fits
READ_FPF = read_fpf.o oFPF.o
FPF2FITS = fpf2fits.o oFPF.o

#########################################################################
all: compile

oFPF.o: oFPF.h oFPF.cpp
	$(CC) $(CFLAGS) -c -o oFPF.o oFPF.cpp

read_fpfp.o: oFPF.h read_fpf.cpp
	$(CC) $(CFLAGS) -c -o read_fpf.o read_fpf.cpp

fpf2fits.o: oFPF.h fpf2fits.cpp oFPF.cpp
	$(CC) $(CFLAGS) -c -o fpf2fits.o fpf2fits.cpp

read_fpf: $(READ_FPF) 
	$(CC) $(CFLAGS) $(READ_FPF) -o read_fpf $(CLIBS)

fpf2fits: $(FPF2FITS)
	$(CC) $(CFLAGS) $(FPF2FITS) -o fpf2fits $(CLIBS)

compile: $(MAIN_PROGS)

clean:
	rm -f *.o $(MAIN_PROGS)




