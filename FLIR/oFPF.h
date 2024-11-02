/*--------------------------------------------------------------------------

  oFPF :  Flir Public Image Format Object Manipulation
          A series of routines to read and convert images recorded with the
          FLIR camera and saved in FPF format. 

          FPF format is the most extensive format delivered by the FLIR
          camera, besides SEQ(FLIR). The main difference: FPF is for 
          images (frames) while SEQ(FLIR) is a format for movies. 

          SEQ(FLIR) is a reformulation of the original open format SEQ, 
          which is proprietary and can be read, and manipulated using the 
          FLIR software only.

          FPF format brings all the header information (meta-data, auxiliary data)
          which is found in the SEQ(FLIR) files. None of the other formats
          give as much information as FPF.

          The class oFPF can read flir FPF files, including the header 
          and convert to FITS, creating a simple image fits format. 

          It needs the CCfits and cfitsio libraries. They are normally included
          with Ubuntu distros. 

          Guigue - 2017-10-24 @Sampa

  ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <valarray>
#include <cmath>
#include <CCfits/CCfits>

#ifndef oFPF_H_
#define oFPF_H_

using namespace std ;
typedef int FPFlong              ;   // signed 4 bytes integer
typedef unsigned int FPFulong    ;   // unsigned 4 bytes integer
typedef unsigned short FPFuint   ;   // unsigned 2 bytes integer
typedef short FPFint             ;   // signed 2 bytes integer
typedef float FPFfloat           ;   // 4 bytes float
typedef char FPFchar             ;   // 1 byte

const short FPF_DESCR_LEN = 32 ;
const short FPF_SPARE = 16     ;
const long MaxfpfDataSize = 1000000 ;

// Start description of the different header sections
// We define every header section as an structure
struct  FPF_IMAGE_DATA_T
{
  FPFchar fpfID[FPF_DESCR_LEN] ; // "FLIR Public Image Format" 
  FPFulong version             ;             
  FPFulong pixelOffset         ; // Offset to pixel values from start of fpfID.   
  FPFuint ImageType            ; // Temperature = 0, 
                                 // Diff Temp = 2, 
                                 // Object Signal = 4,
                                 // Diff Object Signal = 5, etc 
  FPFuint pixelFormat          ; // 0 = short integer = 2 bytes
                                 // 1 = long integer = 4 bytes 
                                 // 2 = float (single precision)  = 4 bytes
                                 // 3 = double (double precision) = 8 bytes 
  FPFuint xSize                ;
  FPFuint ySize                ;
  FPFulong trig_count          ; // external trig counter 
  FPFulong frame_count         ; // frame number in sequence 
} ;

struct FPF_CAMDATA_T {
  FPFchar  camera_name[FPF_DESCR_LEN]  ;
  FPFchar  camera_partn[FPF_DESCR_LEN] ;
  FPFchar  camera_sn[FPF_DESCR_LEN]    ;
  FPFfloat camera_range_tmin           ;
  FPFfloat camera_range_tmax           ;
  FPFchar  lens_name[FPF_DESCR_LEN]    ;
  FPFchar  lens_partn[FPF_DESCR_LEN]   ;
  FPFchar  lens_sn[FPF_DESCR_LEN]      ;
  FPFchar  filter_name[FPF_DESCR_LEN]  ;
  FPFchar  filter_partn[FPF_DESCR_LEN] ;
  FPFchar  filter_sn[FPF_DESCR_LEN]    ;
};

struct FPF_OBJECT_PAR_T
{
  FPFfloat emissivity     ; // 0 - 1 
  FPFfloat objectDistance ; // Meters 
  FPFfloat ambTemp        ; // Reflected Ambient temperature in Kelvin 
  FPFfloat atmTemp        ; // Atmospheric temperature in Kelvin 
  FPFfloat relHum         ; // 0 - 1 
  FPFfloat compuTao       ; // Computed atmospheric transmission 0 - 1
  FPFfloat estimTao       ; // Estimated atmospheric transmission 0 - 1
  FPFfloat refTemp        ; // Reference temperature in Kelvin 
  FPFfloat extOptTemp     ; // Kelvin 
  FPFfloat extOptTrans    ; // 0 - 1 
} ;

struct FPF_DATETIME_T {
   FPFlong Year          ;
   FPFlong Month         ;
   FPFlong Day           ;
   FPFlong Hour          ;
   FPFlong Minute        ;
   FPFlong Second        ;
   FPFlong MilliSecond   ;
} ;

struct FPF_SCALING_T {
  FPFfloat tMinCam       ; //  Camera scale min, in current output 
  FPFfloat tMaxCam       ; //  Camera scale max 
  FPFfloat tMinCalc      ; //  Calculated min (almost true min) 
  FPFfloat tMaxCalc      ; //  Calculated max (almost true max) 
  FPFfloat tMinScale     ; //  Scale min 
  FPFfloat tMaxScale     ; //  Scale max 
} ;

// We assemble all header sections in one unique structure
struct FPFHEADER_T {
  FPF_IMAGE_DATA_T image_data    ;
  FPF_CAMDATA_T    camdata       ;
  FPF_OBJECT_PAR_T object_par    ;
  FPF_DATETIME_T   datetime      ;
  FPF_SCALING_T    scaling       ;
  string           fpf_fname     ;
} ;

class oFPF {
 private:
  struct image_t {
    FPFHEADER_T header  ;
    vector <float> data ; // Data is a 1D vector of type float.
                          // No dimension is given when created.
  } image ;
 public:
  oFPF()                                         ; // Default
  oFPF(const std::string)                        ; // Set the FPF file name
  oFPF(const oFPF &)                             ; // Copy
  oFPF & operator=(const oFPF &)                 ; // Assignment
 ~oFPF()                                         ;
  bool OpenImage (std::string , std::ifstream & );
  bool OpenImage (ifstream &)                    ;
  void CloseImage (std::ifstream & )             ;
  bool ReadImageSize(std::ifstream &)            ;
  bool ReadImage(std::ifstream & )               ;
  bool DimensionDataArray(bool)                  ;
  void ShowHeader() const                        ;
  float getElement(short , short ) const         ;
  bool  writeFITS(std::string & ) const           ;

} ;

#endif
  
