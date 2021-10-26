/*****************************************************************

  HATS_fft

  @guiguesp - Sao Paulo - 2021-10-12 - 

  Compile with:
   gcc -Wall -g windowed_dft.c HATS_fft.c -lm -lfftw3 -o HATS_fft

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

void fft(int num_windows, double * signal, unsigned long * husec, unsigned long * m_husec, double * fftCoeff)
{
  int i;
  for(i = 0; i < num_windows; i++)
    {
      fftCoeff[i] = goertzel_amplitude(&signal[i*STEPS]     ,
				       TARGET_FREQUENCY     ,
				       WINDOW_SIZE          ,
				       SAMPLING_FREQUENCY   ,
				       FLATTOP) / (WINDOW_SIZE)      ;
      m_husec[i] = husec[i*STEPS+WINDOW_SIZE/2];    
    }
}

int main(){
  FILE *fpi, *fpo;
  long num_windows,nrec,i;

  if  (NULL == (fpi = fopen(DATA_FILE,"rb"))) exit(1);
  nrec = bin_size(fpi) / sizeof(int) ;
  fseek(fpi,0L,SEEK_SET);  
  int * adc = malloc(sizeof(int)*nrec);
  unsigned long int * husec = malloc(sizeof(unsigned long int) * nrec);
  double * signal = malloc(sizeof(double) * nrec);
  
  num_windows = floor((nrec - WINDOW_SIZE + 1) / STEPS);
  unsigned long int * m_husec = malloc(sizeof(unsigned long int) * num_windows);
  double * fftAmp = malloc(sizeof(double) * num_windows) ;
  
  fread(adc, sizeof(int), nrec, fpi);
  fclose(fpi);
  for (i=0;i<nrec;i++) signal[i] = (double) adc[i];
  
  if  (NULL == (fpi = fopen(HUSEC_FILE,"rb"))) exit(1);
  fread(husec, sizeof(unsigned long int), nrec, fpi);
  fclose(fpi);
      
  fft(num_windows, signal, husec, m_husec, fftAmp) ;

  if  (NULL == (fpo = fopen(HUSEC_FILE,"wb"))) exit(1);
  fwrite(m_husec, sizeof(long unsigned int), num_windows, fpo);
  fclose(fpo);

  if  (NULL == (fpo = fopen(DATA_FILE,"wb"))) exit(1);
  fwrite(fftAmp, sizeof(double), num_windows, fpo);
  fclose(fpo);

}
