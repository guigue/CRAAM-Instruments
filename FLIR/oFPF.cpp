/*----------------------------------------------------------------------
 *
 * oFPF.cpp : The Class Interface
 *            More comments on the header file and in the lines below
 *
 * Guigue - 2017-10-24 @Sampa
 *
 *---------------------------------------------------------------------- */

#include "oFPF.h"

// Default creator. It creates an empty object.
// The size of data is set after reading the
// FPF header.
oFPF::oFPF() {};

// A constructor to set the name of the FPF file
oFPF::oFPF(const std::string fname) {
  image.header.fpf_fname=fname;
}

// The copy constructor
oFPF::oFPF(const oFPF & fpf){
  image.header = fpf.image.header;
  image.data = fpf.image.data;
}

// The assignment constructor
oFPF & oFPF::operator=(const oFPF & fpf) {
  if (this == &fpf) return *this;
  image.header = fpf.image.header;
  image.data = fpf.image.data;
  return *this;
}
  


// Default Destructor
oFPF::~oFPF() {} ;

// Open the FPP file and set the name in the header structure
bool oFPF::OpenImage (const std::string fname, ifstream &fin ){
  
  fin.exceptions(std::ifstream::failbit);
  try
    {
      fin.open(fname.c_str(), std::ifstream::in | std::ifstream::binary);
    }
  catch (std::ifstream::failure e) 
    {
      cout << endl << endl << "Cannot open file " << fname << endl << endl << endl ;
      return false;
    }      
  image.header.fpf_fname=fname;
  return true ;
}

// The name is already defined in the header. Open the FPP file.
bool oFPF::OpenImage (ifstream &fin ){  
  fin.exceptions(std::ifstream::failbit);
  try
    {
      fin.open(image.header.fpf_fname.c_str(), std::ifstream::in | std::ifstream::binary) ;
    }
  catch (std::ifstream::failure e) 
    {
      cout << endl << endl << "Cannot open file " << image.header.fpf_fname << endl << endl << endl ;
      return false;
    }      
  return true;
}


// Close the FPF file
void oFPF::CloseImage (ifstream & fin){
  fin.close();
  return;
}

/*
 The Image Size is included in the FPF header.
 So we first read the first header section, where
 X and Y sizes are written.  

 The way to read from a file is using a buffer of adequate
 size, and then cast its content to the right format.

*/
bool oFPF::ReadImageSize(ifstream & fin) {
  FPFulong ulbuffer;
  FPFuint  uibuffer;
  FPFchar  buffer[FPF_SPARE*sizeof(FPFlong)];

  fin.exceptions (std::ifstream::failbit | std::ifstream::badbit);

  try
    {
      fin.read(image.header.image_data.fpfID, FPF_DESCR_LEN);

      fin.read((char *) &ulbuffer,sizeof(FPFulong)); 
      image.header.image_data.version = (FPFulong) (ulbuffer) ;

      fin.read((char *) &ulbuffer,sizeof(FPFulong));
      image.header.image_data.pixelOffset = (FPFulong) (ulbuffer);

      fin.read((char *) &uibuffer,sizeof(FPFuint));
      image.header.image_data.ImageType = (FPFuint) (uibuffer);

      fin.read((char *) &uibuffer,sizeof(FPFuint));
      image.header.image_data.pixelFormat = (FPFuint) (uibuffer);

      fin.read((char *) &uibuffer,sizeof(FPFuint));
      image.header.image_data.xSize = (FPFuint) (uibuffer);

      fin.read( (char *) &uibuffer,sizeof(FPFuint));
      image.header.image_data.ySize = (FPFuint) (uibuffer);

      fin.read((char *) &ulbuffer,sizeof(FPFulong));
      image.header.image_data.trig_count = (FPFulong) (ulbuffer);

      fin.read((char *) &ulbuffer,sizeof(FPFulong));
      image.header.image_data.frame_count = (FPFulong) (ulbuffer);

      fin.read(buffer,FPF_SPARE*sizeof(FPFlong));
    }
  catch (std::ifstream::failure e) 
    {
      cout << endl << endl << "Cannot read file " << image.header.fpf_fname << endl << endl << endl ;
      return false;
    }   
  return true;
}

// Resize the 1D vector to host the data

bool oFPF::DimensionDataArray(bool noShow){
  long fpfDataSize ; 
  if  (image.header.image_data.xSize > 0  || image.header.image_data.ySize > 0) {
      fpfDataSize = image.header.image_data.xSize * image.header.image_data.ySize;
      if (fpfDataSize < MaxfpfDataSize)
	{
	image.data.resize(fpfDataSize);
	if (!noShow) cout << endl << "Size of Data Image = " << image.data.size() << endl;
	}
      else
	{
	  cout << endl<< endl<< "Image Sizee bigger than Maximum allowed size = " << MaxfpfDataSize << endl << endl;
	  return false;
	}
    }
  else
    {
      cout << endl << endl << "Problems reading FPF Data Size" << endl << endl ;
      return false;
    }
  return true;
}

/*
  ReadImage: It actually reads the rest of the header, and, finnally
             the image data.

*/
bool oFPF::ReadImage(ifstream & fin) {

  FPFchar null_buffer[FPF_SPARE*sizeof(FPFlong)];

  fin.exceptions (std::ifstream::failbit | std::ifstream::badbit);

  try
    {
  /*******************************************************************
   The Camera Data Header 
  ******************************************************************/
      {
	FPFchar  buffer[FPF_DESCR_LEN];
	FPFfloat fbuffer[2];
    
	fin.read(image.header.camdata.camera_name, FPF_DESCR_LEN)  ;
	fin.read(image.header.camdata.camera_partn,FPF_DESCR_LEN)  ;
	fin.read(image.header.camdata.camera_sn, FPF_DESCR_LEN)    ;
  
	fin.read((char *) &fbuffer,sizeof(FPFfloat)*2)                ;
	image.header.camdata.camera_range_tmin = (FPFfloat) fbuffer[0];
	image.header.camdata.camera_range_tmax= (FPFfloat) fbuffer[1] ;
  
	fin.read(image.header.camdata.lens_name,FPF_DESCR_LEN)     ;
	fin.read(image.header.camdata.lens_partn,FPF_DESCR_LEN)    ;
	fin.read(image.header.camdata.lens_sn,FPF_DESCR_LEN)       ;
	fin.read(image.header.camdata.filter_name,FPF_DESCR_LEN)   ;
	fin.read(image.header.camdata.filter_partn,FPF_DESCR_LEN)  ;
	fin.read(image.header.camdata.filter_sn,FPF_DESCR_LEN)     ;
	fin.read(null_buffer,FPF_SPARE*sizeof(FPFlong))            ; //blank data
      }
  
  /*******************************************************************
   The Object Parameters Header 
  ******************************************************************/
      {
	FPFfloat fbuffer[10];
    
	fin.read((char *) &fbuffer,sizeof(FPFfloat)*10)                ;
	image.header.object_par.emissivity     = (FPFfloat) fbuffer[0] ;
	image.header.object_par.objectDistance = (FPFfloat) fbuffer[1] ;
	image.header.object_par.ambTemp        = (FPFfloat) fbuffer[2] ;
	image.header.object_par.atmTemp        = (FPFfloat) fbuffer[3] ;
	image.header.object_par.relHum         = (FPFfloat) fbuffer[4] ;
	image.header.object_par.compuTao       = (FPFfloat) fbuffer[5] ;
	image.header.object_par.estimTao       = (FPFfloat) fbuffer[6] ;
	image.header.object_par.refTemp        = (FPFfloat) fbuffer[7] ;
	image.header.object_par.extOptTemp     = (FPFfloat) fbuffer[8] ;
	image.header.object_par.extOptTrans    = (FPFfloat) fbuffer[9] ;
	
	fin.read(null_buffer,FPF_SPARE*sizeof(FPFlong))                ; //blank data
      }
  
  /*****************************************************************
   The Date Time Header 
  ******************************************************************/
      {
	FPFlong ibuffer[7] ;
    
	fin.read((char *) &ibuffer,sizeof(FPFlong) * 7)          ;
	image.header.datetime.Year        = (FPFlong) ibuffer[0] ;
	image.header.datetime.Month       = (FPFlong) ibuffer[1] ;
	image.header.datetime.Day         = (FPFlong) ibuffer[2] ;
	image.header.datetime.Hour        = (FPFlong) ibuffer[3] ;
	image.header.datetime.Minute      = (FPFlong) ibuffer[4] ;
	image.header.datetime.Second      = (FPFlong) ibuffer[5] ;
	image.header.datetime.MilliSecond = (FPFlong) ibuffer[6] ;
	fin.read(null_buffer,FPF_SPARE*sizeof(FPFlong))          ; //blank data

      }

  /*******************************************************************
   The Scaling Header 
  ******************************************************************/
      {
	FPFfloat fbuffer[6];
	fin.read((char *) &fbuffer,sizeof(FPFfloat) * 6)    ;
	image.header.scaling.tMinCam   = (FPFfloat) fbuffer[0] ;
	image.header.scaling.tMaxCam   = (FPFfloat) fbuffer[1] ;
	image.header.scaling.tMinCalc  = (FPFfloat) fbuffer[2] ;
	image.header.scaling.tMaxCam   = (FPFfloat) fbuffer[3] ;
	image.header.scaling.tMinScale = (FPFfloat) fbuffer[4] ;
	image.header.scaling.tMaxScale = (FPFfloat) fbuffer[5] ;
	fin.read(null_buffer,FPF_SPARE*sizeof(FPFlong))        ;  // blank data
      }

  /***************************************************************
   The actual Data
  ****************************************************************/
      {
	fin.read(null_buffer,FPF_SPARE * sizeof(FPFlong)) ; // Blank Data
	fin.read(null_buffer,FPF_SPARE * sizeof(FPFlong)) ; // Blank Data
	fin.read(reinterpret_cast<char*>(image.data.data()), image.data.size()*sizeof(FPFfloat));
      }
    }
  catch (std::ifstream::failure e)
    {
      cout << endl << endl << "Error while reading file " << image.header.fpf_fname << endl << endl ;
      return false;
    }
  return true;
}

// Simple look at the Header Sections and its contents  
void oFPF::ShowHeader() const{

  cout << endl << endl << "   Image Header  " << endl ;

  cout << endl << "Image Data" << endl ;
  cout << "fpfID = " << image.header.image_data.fpfID << endl ;
  cout << "version = " << image.header.image_data.version << endl ;
  cout << "pixelOffset = " << image.header.image_data.pixelOffset << endl;
  cout << "ImageType = " << image.header.image_data.ImageType << endl;
  cout << "pixelFormat = " << image.header.image_data.pixelFormat << endl;
  cout << "xSize = " << image.header.image_data.xSize << endl;
  cout << "ySize = " << image.header.image_data.ySize << endl;
  cout << "trig_count = " << image.header.image_data.trig_count << endl;
  cout << "frame_count = " << image.header.image_data.frame_count << endl;
  
  cout << endl << "Camera Data" << endl;
  cout << "Camera Name = " << image.header.camdata.camera_name << endl;
  cout << "Camera Part N = " << image.header.camdata.camera_partn << endl;
  cout << "Camera SN = " << image.header.camdata.camera_sn << endl;
  cout << "Camera Tmin = " << image.header.camdata.camera_range_tmin << endl;
  cout << "Camera Tmax = " << image.header.camdata.camera_range_tmax << endl;
  cout << "Camera Lens name = " << image.header.camdata.lens_name << endl;
  cout << "Camera Lens PartN = " << image.header.camdata.lens_partn << endl;
  cout << "Camera Lens SN = " << image.header.camdata.lens_sn << endl;
  cout << "Camera Filter Name = " << image.header.camdata.filter_name << endl;
  cout << "Camera Filter PartN = " << image.header.camdata.filter_partn << endl;
  cout << "Camera Filter SN = " << image.header.camdata.filter_sn << endl;

  cout << endl << "Object Parameters" << endl;
  cout << "Emisivity = " << image.header.object_par.emissivity << endl ;
  cout << "Object Distance = " << image.header.object_par.objectDistance << endl ;
  cout << "ambTemp = " << image.header.object_par.ambTemp << endl ;
  cout << "atmTemp = " << image.header.object_par.atmTemp << endl ;
  cout << "relHum = " << image.header.object_par.relHum << endl ;
  cout << "compuTao = " << image.header.object_par.compuTao << endl ;
  cout << "estimTao = " << image.header.object_par.estimTao << endl ;
  cout << "refTemp = " << image.header.object_par.refTemp << endl ;
  cout << "extOptTemp = " << image.header.object_par.extOptTemp << endl ;
  cout << "extOptTrans = " << image.header.object_par.extOptTrans << endl ;

  cout << endl << "Date Time Data" << endl ;
  cout << "Year = " << image.header.datetime.Year << " Month = " <<
    image.header.datetime.Month << " Day = " << image.header.datetime.Day << endl;
  cout << "Hour = " << image.header.datetime.Hour << " Minutes = " << image.header.datetime.Minute
       << " Seconds = " << image.header.datetime.Second << " Milliseconds = "
       << image.header.datetime.MilliSecond << endl;

  cout << endl << "Scaling Header " << endl ;
  cout << "tMinCam = " << image.header.scaling.tMinCam << endl ;
  cout << "tMaxCam = " << image.header.scaling.tMaxCam << endl ;
  cout << "tMinCalc = " << image.header.scaling.tMinCalc << endl;
  cout << "tMaxCam = " << image.header.scaling.tMaxCam << endl ;
  cout << "tMinScale = " << image.header.scaling.tMinScale << endl ;
  cout << "tMaxScale = " << image.header.scaling.tMaxScale << endl;

}

// Since we are using a 1D array to represent a 2D Matrix,
// this routine returns the Matrix(f,c) element of the Data Image 
float oFPF::getElement(short f, short c) const
{
  unsigned long image_index;
  
  if (f < image.header.image_data.ySize && f >= 0 &&
      c < image.header.image_data.xSize && c >= 0 )
    {
      image_index = c + f * image.header.image_data.xSize ;
      return image.data[image_index] ;
    }
   else
     {
       return -999;
     } 

}

/*
  writeFITS: Creates a FITS file from the image read. The FITS file is a simple
             image file containg just a Primary Data Header Unit (PHU) without 
             extensions.
             
             The PHU is created when CCfits::FITS is invoked and includes some
             fixed tags and comments. 

             Afterwards we add some optional tags like ORIGIN and T_START, etc.
             Finnally the rest of the image header is included, in separate sections
             identified by comments.

 */
bool oFPF::writeFITS(string & fitsname ) const
{
  char buff1[10], buff2[25]       ;
  string fitsAname="!"+fitsname   ;         // fits Alternate filename, to overwrite in case it exist
  string Date_Obs, T_Start        ;
  long naxis = 2                  ;
  long nelements=0, first_pixel=1 ;
  long naxes[2] = {image.header.image_data.xSize , image.header.image_data.ySize};

  // CCfits uses std::valarray as data, and I prefer to use std::vector
  // The temporary valarray is used to transfer the data to the fits file.
  // There is no way to cast a <vector> to a <valarray> , sniff...
  valarray<float> temp_data(image.data.data(),image.data.size());
    
  /*
     From the Manual of CCfits-2.5

     Declare auto_pointer to FITS at function scope. 
     Ensures no resources leaked if something fails 
     in dynamic allocation.
  */  
  std::auto_ptr<CCfits::FITS> pFits(0);

  nelements = naxes[0]*naxes[1] ;

  // The DATE-OBS fomatted keyword
  sprintf(buff1,"%4d-%2.2d-%2.2d",
	  image.header.datetime.Year,
	  image.header.datetime.Month,
	  image.header.datetime.Day);
  Date_Obs=buff1;

  // The T_START formated keyword
  sprintf(buff2,"%4d-%2.2d-%2.2dT%2.2d:%2.2d:%2.2d.%3.3dZ",
	  image.header.datetime.Year,
	  image.header.datetime.Month,
	  image.header.datetime.Day,
	  image.header.datetime.Hour,
	  image.header.datetime.Minute,
	  image.header.datetime.Second,
	  image.header.datetime.MilliSecond);
  T_Start=buff2;

  try
    {
  /* 
    From the Manual of CCfits-2.5
     Create a new FITS object, specifying the data 
     type and axes for the primary image. Simultaneously 
     create the corresponding file. This image is float
     data, demonstrating the cfitsio extension to the 
     FITS standard.

     The fitsAname (starting with a !) ensures the creation 
     of the file even if it already exists.

  */
      pFits.reset(new CCfits::FITS(fitsAname, FLOAT_IMG, naxis, naxes)) ;
    }
  catch (CCfits::FITS::CantCreate)
    {
      cout << endl << endl << "Cannot create fits file." << endl << endl ;
      return false;
    }

  /* 
       Kind of Subroutine. It writes the header to the fits file.
       Boring... 
  */

  // These are keywords considered mandatory.
  pFits->pHDU().addKey("ORIGIN","CRAAM/UPM"," ");
  pFits->pHDU().addKey("TELESCOP","Mid-IR","  ");
  pFits->pHDU().addKey("FREQ","30","THz");
  pFits->pHDU().addKey("ORIGIFIL",image.header.fpf_fname,"Original FPF File Name");
  pFits->pHDU().addKey("DATE-OBS", Date_Obs," ");
  pFits->pHDU().addKey("T_START", T_Start," ");
  pFits->pHDU().addKey("LVL_NUM","0.0","Raw Data");

  // These are the values obtained from the camera.
  pFits->pHDU().writeComment(" ");
  pFits->pHDU().writeComment(" ");
  pFits->pHDU().writeComment("--------- Original Camera Header Cards ---------");
  pFits->pHDU().writeComment(" ");
  pFits->pHDU().writeComment("Image Header");
  pFits->pHDU().addKey("fpfID", image.header.image_data.fpfID," ");
  pFits->pHDU().addKey("version", image.header.image_data.version," ");
  pFits->pHDU().addKey("pxOffset",image.header.image_data.pixelOffset," ");
  pFits->pHDU().addKey("ImgType",image.header.image_data.ImageType," ");
  pFits->pHDU().addKey("pxFormat",image.header.image_data.pixelFormat," ");
  pFits->pHDU().addKey("xSize",image.header.image_data.xSize, " ");
  pFits->pHDU().addKey("ySize",image.header.image_data.ySize," ");
  pFits->pHDU().addKey("trgCount",image.header.image_data.trig_count," ");
  pFits->pHDU().addKey("frmCount",image.header.image_data.frame_count," ");
  
  pFits->pHDU().writeComment(" ");
  pFits->pHDU().writeComment("Camera Data");
  pFits->pHDU().addKey("Name",image.header.camdata.camera_name,"  ");
  pFits->pHDU().addKey("Part_N",image.header.camdata.camera_partn," ");
  pFits->pHDU().addKey("SN",image.header.camdata.camera_sn,"  ");
  pFits->pHDU().addKey("Tmin",image.header.camdata.camera_range_tmin,"  ");
  pFits->pHDU().addKey("Tmax",image.header.camdata.camera_range_tmax,"  ");
  pFits->pHDU().addKey("L_Name",image.header.camdata.lens_name,"Lens Name");
  pFits->pHDU().addKey("L_Part",image.header.camdata.lens_partn,"Lens Part");
  pFits->pHDU().addKey("L_SN",image.header.camdata.lens_sn,"Lens SN");
  pFits->pHDU().addKey("F_Name",image.header.camdata.filter_name,"Filter Name");
  pFits->pHDU().addKey("F_PartN",image.header.camdata.filter_partn,"Filter Part N");
  pFits->pHDU().addKey("F_SN",image.header.camdata.filter_sn,"Filter SN");

  pFits->pHDU().writeComment(" ");
  pFits->pHDU().writeComment("Object Parameters");
  pFits->pHDU().addKey("Emisivit",image.header.object_par.emissivity,"  ");
  pFits->pHDU().addKey("ObjDist",image.header.object_par.objectDistance,"Object Distance");
  pFits->pHDU().addKey("ambTemp",image.header.object_par.ambTemp,"  ");
  pFits->pHDU().addKey("atmTemp",image.header.object_par.atmTemp,"  ");
  pFits->pHDU().addKey("relHum",image.header.object_par.relHum,"  ");
  pFits->pHDU().addKey("compuTao",image.header.object_par.compuTao,"  ");
  pFits->pHDU().addKey("estimTao",image.header.object_par.estimTao,"  ");
  pFits->pHDU().addKey("refTemp",image.header.object_par.refTemp,"  ");
  pFits->pHDU().addKey("eOptTemp",image.header.object_par.extOptTemp,"Originally: extOptTemp");
  pFits->pHDU().addKey("eOptTrns",image.header.object_par.extOptTrans,"Originally: extOptTrans  ");

  pFits->pHDU().writeComment(" ");
  pFits->pHDU().writeComment("Scaling Header");
  pFits->pHDU().addKey("tMinCam",image.header.scaling.tMinCam," ");
  pFits->pHDU().addKey("tMaxCam",image.header.scaling.tMaxCam," ");
  pFits->pHDU().addKey("tMinCalc",image.header.scaling.tMinCalc," ");
  pFits->pHDU().addKey("tMaxCam",image.header.scaling.tMaxCam," ");
  pFits->pHDU().addKey("tMinScal",image.header.scaling.tMinScale," ");
  pFits->pHDU().addKey("tMaxScal",image.header.scaling.tMaxScale," ");

  // NOW... the actual Data
  pFits->pHDU().write(first_pixel,nelements,temp_data);

  return true;
  
}
