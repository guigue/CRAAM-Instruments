/*****************************************************************

  HATS_fft

  @guiguesp - Sao Paulo - 2021-12-17 - 

  Compile with:
   gcc -Wall -g windowed_dft.c HATS_fft.c -lm -lfftw3 -o HATS_fft

  usage:
    HATS_fft [window_size steps]
             where window_size = number of records taken to get the FFT
                   steps = shift between windows
                   steps < window_size

******************************************************************/

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <semaphore.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <complex.h>

#include "HATS_fft.h"
#include "windowed_dft.h"

unsigned long bin_size(FILE * fp)
{
  fseek(fp, 0L, SEEK_END);
  return ftell(fp);
}

void fft(long int window_size,
	 long int steps,
	 long int num_windows,
	 double * signal,
	 unsigned long * husec,
	 unsigned long * m_husec,
	 double * fftCoeff)
{
  int i;
  for(i = 0; i < num_windows; i++)
    {
      fftCoeff[i] = goertzel_amplitude(&signal[i*steps]     ,
				       TARGET_FREQUENCY     ,
				       window_size          ,
				       SAMPLING_FREQUENCY   ,
				       FLATTOP) / (window_size)      ;
      m_husec[i] = husec[i*steps+window_size/2];    
    }
}

int main(int argc, char ** argv){
  FILE *fpi, *fpo;
  long int num_windows, nrec, window_size, steps;

  if (argc != 3)
    {
      window_size = WINDOW_SIZE;
      steps = STEPS;
    }
  else
    {
      window_size = atoi(argv[1]);
      steps = atoi(argv[2]);
    }

  if (steps > window_size) steps = window_size;
  
  if  (NULL == (fpi = fopen(DATA_FILE,"rb"))) exit(1);
  nrec = bin_size(fpi) / sizeof(double) ;
  fseek(fpi,0L,SEEK_SET);  
  //int * adc = malloc(sizeof(int)*nrec);
  unsigned long int * husec = malloc(sizeof(unsigned long int) * nrec);
  double * signal = malloc(sizeof(double) * nrec);
  
  num_windows = floor((nrec - window_size + 1) / steps);
  unsigned long int * m_husec = malloc(sizeof(unsigned long int) * num_windows);
  double * fftAmp = malloc(sizeof(double) * num_windows) ;
  
  fread(signal, sizeof(double), nrec, fpi);
  fclose(fpi);
  //  for (i=0;i<nrec;i++) signal[i] = (double) adc[i];
  
  if  (NULL == (fpi = fopen(HUSEC_FILE,"rb"))) exit(1);
  fread(husec, sizeof(unsigned long int), nrec, fpi);
  fclose(fpi);
      
  fft(window_size,steps,num_windows, signal, husec, m_husec, fftAmp) ;

  if  (NULL == (fpo = fopen(HUSEC_FILE,"wb"))) exit(1);
  fwrite(m_husec, sizeof(long unsigned int), num_windows, fpo);
  fclose(fpo);

  if  (NULL == (fpo = fopen(DATA_FILE,"wb"))) exit(1);
  fwrite(fftAmp, sizeof(double), num_windows, fpo);
  fclose(fpo);

}
