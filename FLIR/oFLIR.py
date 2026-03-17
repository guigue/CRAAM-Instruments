# Requirements:
# pip install astropy
# pip install install/FileSDK-2024.7.1-cp312-cp312-linux_x86_64.whl 
# (corresponding python version)

"""
Running the script: 

 python seq2fits_30THz path_to_seq_files telescope outdir
 
 path_to_seq_files = the path where the .seq are located (e.g. /home/user/data/seq/)
 telescope = AR30T or SP30T according to the corresponding .seq (do not mix SP30T and AR30T .seq files)
 outdir = path to save the .fits files
 
 Return: 
 - The code will extract the .fits files and create header for every frame (different frames will be created for AR30T and SP30T)
 - Using outdir as base the script will create (if not existing) the following directories: outdir + /telescope/yyyymmdd/level00/
 - filename examples: AR30T_20250612T123030.870_level00.fits ; SP30T_20250612T123030_level00.fits

Created by Ian Grinkraud (CRAAM-Universidade Presbiteriana Mackenzie)
Modified by Fernando Lopez (Universidad de Mendoza-CONICET & Observatorio Astronómico Félix Aguilar) 
Last Modification: 2025-12-02

"""
from datetime import datetime, timezone
from astropy.coordinates import get_sun, EarthLocation
from astropy.time import Time
import astropy.units as u
from astropy.constants import R_sun
import os
import glob


from astropy.io import fits
import numpy as np
import datetime
import fnv
import fnv.reduce
import fnv.file
import sys
import os
import shlex

class Telescope:
    
    def __init__(self, name):
        self.name = name
        if name == 'AR30T':
            self.instrument = 'FLIR A645sc'
            self.origin = 'OAFA/UNSJ - CRAAM/UPM'
            self.observatory = 'Obs. Astronomico Felix Aguilar'
            self.place = 'El Leoncito - San Juan, Argentina'
            self.longitude = -69.318889
            self.latitude = -31.802222
            self.elevation = 2430
            self.comments = ['Level 0.0: uncalibrated image','COPYRIGHT. Grant of use.',
                             'These data are property of Universidad Presbiteriana Mackenzie and Observatorio Astronomico Felix Aguilar.',
                             'The Centro de Radio Astronomia e Astrofisica Mackenzie',
                             'and Observatorio Astronomico Felix Aguilar are responsable for their distribution.',
                             'Grant of use permission is given for Academic purposes only.',
                             'Contact:guigue@craam.mackenzie.br; fernando.lopez@um.edu.ar' ]
        elif name == 'SP30T':
            self.instrument = 'FLIR A20'
            self.origin = 'CRAAM/UPM'
            self.observatory = 'Centro de Radio Astronomia e Astrofisica Mackenzie'
            self.place = 'Sao Paulo, Brazil'
            self.longitude = -46.633308
            self.latitude = -23.550520
            self.elevation = 760
            self.comments = ['Level 0.0: uncalibrated image','COPYRIGHT. Grant of use.',
                             'These data are property of Universidad Presbiteriana Mackenzie',
                             'The Centro de Radio Astronomia e Astrofisica Mackenzie', 
                             'is responsable for their distribution.',
                             'Grant of use permission is given for Academic purposes only.',
                             'Contact:guigue@craam.mackenzie.br' ]
        else:
            raise ValueError("Unsupported telescope. Choose 'AR30T' or 'SP30T'.")

        return

    def create_fits_header(self,frameInfo, imageRoi, objPar, width, height, in_file, telescope, fits_filename=None):
   
        self.hdr = fits.Header()

        # General information
        hdr['SIMPLE']   = True
        hdr['BITPIX']   = -32
        hdr['NAXIS']    = 2
        hdr['NAXIS1']   = (width, ' ')
        hdr['NAXIS2']   = (height, ' ')
        hdr['EXTEND']   = (True, ' ')
        hdr['ORIGIFIL'] = (os.path.basename(in_file), ' ')
        hdr['TELESCOP'] = (telescope.name, ' ')
        hdr['INSTRUME'] = (telescope.instrument, ' ')
        hdr['ORIGIN']   = (telescope.origin, ' ')
        hdr['DATE']     = (datetime.now(timezone.utc).strftime('%Y-%m-%d'), ' ')
        hdr['DATE-OBS'] = (frameInfo.time.isoformat(), ' ')
        hdr['LVL_NUM']  = ('0.0', ' ')
        hdr['EXPTYPE']  = ('EXPOSURE', ' ')
        hdr['WAVELNTH'] = (10, 'center wavelength of bandpass filter')
        hdr['WAVEUNIT'] = ('micrometers', ' ')
        hdr['OBSERVAT'] = (telescope.observatory, ' ')
        hdr['PLACE']    = (telescope.place, ' ')
        hdr['LONGITUD'] = (telescope.longitude, ' ')
        hdr['LATITUDE'] = (telescope.latitude, ' ')
        hdr['ELEVATIO'] = (telescope.elevation, 'Altitude in meters')
        hdr['CUNIT1']   = ('arcsec', ' ')
        hdr['CRVAL1']   = (0.0, ' ')
        hdr['CRPIX1']   = (width / 2, ' ')
        hdr['CUNIT2']   = ('arcsec', ' ')
        hdr['CRVAL2']   = (0.0, ' ')
        hdr['CRPIX2']   = (height / 2, ' ')
        
        # Original header from .seq file
        hdr['PXOFFSET'] = (892, 'image_pixel_offset ')
        hdr['IMGTYPE']  = (frameInfo.preset, 'image_type ')
        hdr['PXFORMAT'] = (2, 'image_pixel_format ')
        hdr['TRGCOUNT'] = (0, 'image_trig_count')
        hdr['FRMCOUNT'] = (frameInfo.frame, 'image_frame_count')
        hdr['FRM_TIME'] = (frameInfo.time.isoformat() + 'Z', 'frame_time')
        
        # Object parameters
        hdr['EMISIVIT'] = (round(objPar.emissivity, 2), 'object emissivity')
        hdr['OBJDIST']  = (round(objPar.distance, 2), 'object distance ')
        hdr['AMBTEMP']  = (round(objPar.atmosphere_temp, 2), 'amb_temperature ')
        hdr['ATMTEMP']  = (round(objPar.atmosphere_temp, 2), 'atm_temperature ')
        hdr['RELHUM']   = (round(objPar.relative_humidity, 2), 'rel_humidity ')
        hdr['COMPUTAO'] = (round(objPar.atmospheric_transmission, 6), 'compu_tao ')
        hdr['ESTIMTAO'] = (round(objPar.est_atmospheric_transmission, 2), 'estim_tao ')
        hdr['REFTEMP']  = (round(objPar.reflected_temp, 2), 'ref_temp ')
        hdr['EOPTTEMP'] = (round(objPar.ext_optics_temp, 2), 'ext_opt_temp ')
        hdr['EOPTTRNS'] = (round(objPar.ext_optics_transmission, 2), 'ext_opt_trans ')
        
        # Image scales
        hdr['TMINCAM']  = (round(imageRoi.min_value, 4), 't_min_cam')
        hdr['TMAXCAM']  = (round(imageRoi.max_value, 4), 't_max_cam')
        hdr['TMINCALC'] = (round(imageRoi.min_value, 4), 't_min_calc')
        hdr['TMAXSCAL'] = (round(imageRoi.max_value, 4), 't_max_scale')
        hdr['TMINSCAL'] = (round(imageRoi.min_value, 4), 't_min_scale')
        
        # Add comments
        hdr.add_comment(telescope.comments[0])
        hdr.add_comment(telescope.comments[1])
        hdr.add_comment(telescope.comments[2])
        hdr.add_comment(telescope.comments[3])
        hdr.add_comment(telescope.comments[4])
        hdr.add_comment(telescope.comments[5])
    
        #hdr.add_comment('Level 0.0: uncalibrated image')
        #hdr.add_comment('COPYRIGHT. Grant of use.')
        #hdr.add_comment('These data are property of Universidade Presbiteriana Mackenzie')
        #hdr.add_comment('and Observatorio Astronomico Felix Aguilar.')
        #hdr.add_comment('The Centro de Radio Astronomia e Astrofisica Mackenzie and Observatorio Astronomico Felix Aguilar are responsable for their distribution.')
        #hdr.add_comment('Grant of use permission is given for Academic purposes only.')
        #hdr.add_comment('Contact: guigue@craam.mackenzie.br; fernando.lopez@um.edu.ar')

        return

    def define_fits_filename(self,base_dir):
   
        # Extract necessary information from the header
        telescope = header['TELESCOP']
        date_obs = header['DATE-OBS']  # Format: YYYY-MM-DDTHH:MM:SS
        level = header['LVL_NUM'].replace('.', '')
        
        # Parse the observation date and time
        date_part, time_part = date_obs.split('T')
        yyyymmdd = date_part.replace('-', '')
        hhmmss_sss = time_part.replace(':', '')
        hhmmss_sss = f"{float(hhmmss_sss):.3f}"
        
        # Define the directory structure
        output_dir = os.path.join(base_dir, telescope, yyyymmdd, f"level{level}")
        os.makedirs(output_dir, exist_ok=True)

        # Construct the filename
        fits_filename = os.path.join(output_dir, f"{telescope}_{yyyymmdd}T{hhmmss_sss}_level{level}.fits")

        return fits_filename


    def seq2fits(self,seq_directory, outdir):
        
    """
    Process all .seq files in the given directory and extract .fits files.

    Parameters:
        seq_directory (str): Path to the directory containing .seq files.
        telescope (Telescope): Telescope object containing metadata.
    """

        self.input_directory = seq_directory
        self.outdir = outdir
        
        # Find all .seq files in the directory
        seq_files = glob.glob(f"{seq_directory}/*.seq")
        seq_files.sort()

        for seq_file in seq_files:
            print(f"Processing file: {seq_file}")

            # Open the .seq file
            im = fnv.file.ImagerFile(seq_file)
            frame_number = im.first_frame_number(fnv.Preset.ANY)[0]

            while frame_number is not None:
                if im.get_frame(frame_number):
                    frameInfo = im.frame_info
                    imageRoi = im.rois[0]

                    width = im.width
                    height = im.height
                    final = im.final
                    objPar = im.object_parameters

                    # Get data for .fits
                    data = np.array(final, dtype=np.float32).reshape((height, width))

                    # Create header for each .fits image
                    #fits_filename = define_fits_filename(hdr)
                    hdr = create_fits_header(self,frameInfo, imageRoi, objPar, width, height, seq_file)

                    hdu = fits.PrimaryHDU(data=data, header=hdr)
                    hdul = fits.HDUList([hdu])

                    # Define filename and save .fits files
                    fits_filename = define_fits_filename(hdr, base_dir)
                    hdul.writeto(fits_filename, overwrite=True)

                    print(f"fits file {frame_number} saved in: {fits_filename}")
                    print(f"Full path of saved .fits file: {os.path.abspath(fits_filename)}")

                    # Go to the next frame
                    frame_number = im.next_frame_number(frame_number, fnv.Preset.ANY)[0]

                    # Close the .seq file
                    im.close()

                    print("Processing of all .seq files completed.")

            return
        

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python seq2fits_AR30T.py <input_directory> <telescope> <base_dir>")
        sys.exit(1)

    seq_directory = sys.argv[1]  # Get the input directory from the command-line argument
    telescope_name = sys.argv[2]  # Get the telescope name from the command-line argument
    base_dir = sys.argv[3]  # Get the base directory from the command-line argument
    if not os.path.exists(base_dir):
        print(f"Error: The base directory '{base_dir}' does not exist.")
        sys.exit(1)

    try:
        telescope = Telescope(telescope_name)
    except ValueError as e:
        print(e)
        sys.exit(1)

    # Process the .seq files in the provided directory
    process_seq_files(seq_directory, telescope)
