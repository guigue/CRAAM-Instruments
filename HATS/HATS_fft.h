#ifndef HATS_FFT
#define HATS_FFT

#define DATA_FILE "hats_data_rbd.bin"
#define HUSEC_FILE "hats_husec.bin"
#define OUTPUT_FILE "hats_fft_amplitudes.bin"

#define TARGET_FREQUENCY 20                     // frequency that we are interested 
#define WINDOW_SIZE 128                         // number of sample point that compose each DFT analysis
#define STEPS 32                                // number of "walked" sample points for each window
#define SAMPLING_INTERVAL 0.001
#define	SAMPLING_FREQUENCY 1/SAMPLING_INTERVAL
#define SIGNAL_LENGTH 8192
#define NREAD 64 

#endif
