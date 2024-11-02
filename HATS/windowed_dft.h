#ifndef WINDOWED_DFT
#define WINDOWED_DFT

#include <complex.h> //to handle the complex arrays
#include <math.h> //sine, cosine, and fabs functions in the windows filters
#include <fftw3.h> //for fft functions

#define RECTANGULAR 0
#define FLATTOP 1
#define CHEBYSHEV 2
#define HANN 3
#define TRIANGULAR 4

/* Functions here! */

void windowed_dft(double *signal, complex *fft_out, double *windowed_signal, double *mod, double *arg, int window_type, int saved_files, int signal_lenght);

double goertzel_amplitude(double *signal, double target_frequency, int signal_lenght, int sampling_frequency, int window_type);

void mod_array(complex *in, double *out, int signal_lenght);

void arg_array(complex *in, double *out, int signal_lenght);

void rectangularwin(double *signal, double *windowed_signal, int signal_lenght);

void flattopwin(double *in, double *out, int signal_lenght);

//void dolphwin(double *in, double *out, int signal_lenght, int attenuation); //Not properly implemented

void hannwin(double *in, double *out, int signal_lenght);

void triangwin(double *in, double *out, int signal_lenght);

#endif
