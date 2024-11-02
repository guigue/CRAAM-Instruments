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

#include "getPos.h"

#define DATA_FILE "/data/HATS/hats-2021-10-27T0200.aux"
  
int main(){
	   
  FILE *fpi;
  pos_data_type data;
  char * buffer;

  buffer = (char *) &data;

  printf("\n\n size of data = %d\n\n",sizeof(data));  
  fpi = fopen(DATA_FILE,"rb");
  fread(buffer, sizeof(pos_data_type), sizeof(pos_data_type), fpi);

  fclose(fpi);
}



