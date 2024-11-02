#include <string>
#include <iostream>

#include "oFPF.h"

int main( int argc, char* argv[]) 
{
  std::string fpfname,fitsname;
  bool noshow=false;
  int f=0,c=0;
  ifstream file_p;
  oFPF fpf_obj;

  if (argc < 2) {
    std::cout << endl << endl << "Usage:" << endl ;
    std::cout << "    " << argv[0] << " FPFfilename (without extension)" << endl;
    std::cout << "    " << "The FITS file will have the same name." << endl << endl;
    exit(0);
      }
  
  fpfname  = (string) argv[1]+".fpf";
  fitsname = (string) argv[1]+".fits";
  
  fpf_obj.OpenImage(fpfname, file_p);
  fpf_obj.ReadImageSize(file_p);
  fpf_obj.DimensionDataArray(noshow);
  fpf_obj.ReadImage(file_p);
  fpf_obj.CloseImage(file_p);
  fpf_obj.writeFITS(fitsname);

}
