#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdlib.h>

#include "oFPF.h"

int main( int argc, char* argv[]) 
{
  std::string fname, fitsname;
  bool noshow=false;
  int f=0,c=0;
  ifstream file_p;
  oFPF fpf_obj;
  
  if (argc == 1) {
    std::cout << "Exit..." << endl;
    exit(0);
  }
  
  if (argc == 2) {
    fname=(string) argv[1]+".fpf";
    fitsname=(string) argv[1]+"fits";
  }
  
  if (! fpf_obj.OpenImage(fname, file_p)) exit(1);
  if (! fpf_obj.ReadImageSize(file_p)) exit(1) ;
  if (! fpf_obj.DimensionDataArray(noshow)) exit(1);
  if (! fpf_obj.ReadImage(file_p)) exit(1) ;
  fpf_obj.CloseImage(file_p);

  oFPF fpf_obj_2 = fpf_obj;

  //fpf_obj.ShowHeader();

  if (argc == 2) {
    string altFITSname="test.fits" ;
    if (! fpf_obj.writeFITS(fitsname)) exit(1);
    if (! fpf_obj_2.writeFITS(altFITSname)) exit(1);
  }

  exit(0);
    
}
