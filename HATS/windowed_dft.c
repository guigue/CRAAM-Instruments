#include "windowed_dft.h"

/* Windowed DFT for HATS using FFTW. 

   Based on Tiago Giorgetti's implementation in Matlab, 
   "Influência das janelas de ponderação na  obtenção da amplitude de sinais periódicos no domínio da frequência".
   Goertzel Algorithm taken from references/Special_fft_cases.pdf

Description: 
	Calculates the DFT using FFTW algorithm after applying a weight on the whole signal to attenuate spectral leakage. 
	Saves on <fft_out>, <windowed_signal>, <mod>, and <arg> arrays 
Obs:
	Please note that this algorithm assumes that the window size IS the whole signal passed to the function.
Recieves:
	<signal> pointer to double: is an array with lenght <signal_lenght>
	<fft_out> pointer to complex: is an array with lenght floor(<signal_lenght>/2)+1. (This is due to periodic properties of the real valued dft)
	<windowed_signal> pointer to double: is an array with size <signal_lenght> and will store the filtered, or windowed, signal
	<mod> pointer to double: array that will be stored the modulus of the DFT <signal_lenght> (This is easily optimized, it is odd symmetric) 
	<arg> pointer to double: array that will be store the principal argument of the DFT sized <signal_lenght>
	<window_type> int: defines the type of filter to be applyied to cut the signal
		-> 0, Rectangular
		-> 1, Flat top
		-> 2, Dolph-Chebyshev (Not implemented)
		-> 3, Hanning
		-> 4, Triangular
	<saved_files> int: defines which data is calculated 
		-> 0, dft complex array + the amplitude or modulus
		-> 1, dft complex array + amplitude + principal argument or phase
	<signal_lenght> int: number of data points in the window
Returns:
	None
Programmed by Manuel Giménez de Castro in 2021.*/

void windowed_dft(double *signal,
		  complex *fft_out,
		  double *windowed_signal,
		  double *mod,
		  double *arg,
		  int window_type,
		  int saved_files,
		  int signal_lenght){
  int i;
  fftw_plan p;
  switch(window_type){
  case 0:
    rectangularwin(signal, windowed_signal, signal_lenght);
    break;
  case 1:
    flattopwin(signal, windowed_signal, signal_lenght);
    break;
    /*case 2:
      Is way different. This is defined in terms of the frequency response. 
         We need to calculate the inverse dft and of an constant signal(?) 
         then normalize it by its max norm and then apply it to the signal. Weird......

      dolphwin(signal, windowed_signal, signal_lenght, 100); 
      break;
    */
  case 3:
    hannwin(signal, windowed_signal, signal_lenght);
    break;
  case 4:
    triangwin(signal, windowed_signal, signal_lenght);
    break;
  default:
    //just return if a non valid type of window is selected
    return;
    break;
  }
  p = fftw_plan_dft_r2c_1d(signal_lenght, windowed_signal, fft_out, FFTW_ESTIMATE); //apply the DFT to the windowed signal
  fftw_execute(p); /* repeat as needed */
  fftw_destroy_plan(p);
  switch (saved_files){
  case 0:
    mod_array(fft_out, mod, signal_lenght); //calculates the amplitude or modulus data array
    break;
  case 1:	
    mod_array(fft_out, mod, signal_lenght); 
    arg_array(fft_out, arg, signal_lenght); //calculates the argument or phase
    break;
  default:
    return;
    break;		
  }  
  return;
}

/* goertzel_amplitude

   Description: Returns only the amplitude of <target_frequency> of the <signal> array using the goertzel algorithm.
   Recieves:
	<signal>              pointer to complex: is an array with the dft of the signal with lenght <signal_lenght>
	<target_frequency>    pointer to double: is an array with lenght floor(<signal_lenght>/2)+1. 
                              (This is due to periodic properties of the real valued dft)
	<signal_lenght>       int number of data points in the signal
	<sampling_frequency>  int frequency of the samples
	<window_type>         int defines the type of filter to be applyied to cut the signal
		-> 0, Rectangular
		-> 1, Flat top
		-> 2, Dolph-Chebyshev (Not implemented)
		-> 3, Hanning
		-> 4, Triangular
   Returns: Double with the amplitude of the <target_frequency_k> frequency*/

double goertzel_amplitude(double *signal,
			  double target_frequency,
			  int signal_lenght,
			  int sampling_frequency,
			  int window_type){
  int i      ; 
  double Q[3];                                                 //coeficients of our expansion
  double windowed_signal[signal_lenght];                       //filtered signal
                                                               
  const int target_frequency_k = (int) floor(target_frequency*signal_lenght/sampling_frequency); // Frequency to be calculated
  const double A = 2 * cos(2 * M_PI * target_frequency_k / signal_lenght);
  
  //apply filter
  switch (window_type){
  case 0:
    rectangularwin(signal, windowed_signal, signal_lenght);
    break;
  case 1:
    flattopwin(signal, windowed_signal, signal_lenght);
    break;
    case 2:
      return 0;
      /* Is way different. This is defined in terms of the frequency response. 
	 We need to calculate the inverse dft and of an constant signal(?) 
	 then normalize it by its max norm and then apply it to the signal. Weird......

	 dolphwin(signal, windowed_signal, signal_lenght, 100); //With attenuation parameter of -100 dB.
      */
      break;
  case 3:
    hannwin(signal, windowed_signal, signal_lenght);
    break;
  case 4:
    triangwin(signal, windowed_signal, signal_lenght);
    break;
  default:
    //just return 0 if a non valid type of window is selected
    return 0;
    break;
  }
  //Initial values of the coeficients
  Q[0] = 0;
  Q[1] = windowed_signal[0];
  for (i = 1; i < signal_lenght; i++){
    Q[(i+1)%3] = windowed_signal[i] + A*Q[i%3] - Q[(i-1)%3];
  }
  return sqrt(pow(Q[(i-1)%3],2) + pow(Q[(i-2)%3],2) - A*Q[(i-1)%3]*Q[(i-2)%3]);
}

/* mod_array

Description: Calculates the modulus or amplitude of an array of complex numbers.
Recieves:
	<in> pointer to complex: is an array with the dft of the signal with lenght <signal_lenght>
	<out> pointer to double: is an array with lenght floor(<signal_lenght>/2)+1. (This is due to periodic properties of the real valued dft)
	<signal_lenght> int: number of data points in the signal
Returns:
	None

*/

void mod_array(complex *in, double *out, int signal_lenght){
	int i, half;
	//rounded integral of the halved signal lenght. This is because we don't have enough sampling to work with those higher frequencies 
	half = (int) floor(signal_lenght/2)+1;	
	for(i = 0; i < half; i++){
		out[i] = cabs(in[i]);
	}
}

/* arg_array

Description: Calculates the phase or principal argument of an array of complex numbers. Returns an array of real numbers.
Recieves:
	<in> pointer to complex: is an array with the dft of the signal with lenght <signal_lenght>
	<out> pointer to double: is an array with lenght floor(<signal_lenght>/2)+1. (This is due to periodic properties of the real valued dft)
	<signal_lenght> int: number of data points in the signal
Returns:
	None

*/

void arg_array(complex *in, double *out, int signal_lenght){
	int i, half;
	//rounded integral of the halved signal lenght. This is because we don't have enough sampling to work with those higher frequencies 
	half = (int) floor(signal_lenght/2)+1;
	
	for(i = 0; i < half; i++){
		out[i] = cproj(in[i]);
	}
	return;
}

/* rectangularwin

Description: Applyies the rectangular filter to the signal.
Recieves:
	<in> pointer to double: is the array with signal of lenght <signal_lenght>
	<out> pointer to double: is the array after applying the filter with lenght <signal_lenght>.
	<signal_lenght> int: number of data points in the signal
Returns:
	None

*/

void rectangularwin(double *in, double *out, int signal_lenght){
	int i;
	const double correction_param = 2;   //multiplying by the amplitude calibration coeficient.
	                                     // This makes sure that the sinusoidal DFT windowed by this will be accurate in its amplitude 
	for(i=0; i<signal_lenght; i++)  out[i] = correction_param*in[i]; 
}

/* flattopwin

Description: Applyies the Symmetric Flat top window to the signal. Taken from "type 0", ISO 18431-1 (Chique, né?)
	     https://www.dadisp.com/webhelp/mergedProjects/refman2/FncrefFK/FLATTOP.htm
Recieves:
	<in> pointer to double: is the array with signal of lenght <signal_lenght>
	<out> pointer to double: is the array after applying the filter with lenght <signal_lenght>.
	<signal_lenght> int: number of data points in the signal
Returns:
	None
*/

void flattopwin(double *in, double *out, int signal_lenght){
	int i;
	double weight;
	const double correction_param = 2.000122419602196; //with the corrected coeficient by amplitude.
	                                                   // Differs from Tiago Giorgetti and I don't know why.
	for (i=0; i<signal_lenght;i++){
		weight = 1.0 - 1.9330 * cos(2*M_PI*i/(signal_lenght-1)) +
		  1.2860*cos(4*M_PI*i/(signal_lenght-1)) -
		  0.388*cos(6*M_PI*i/(signal_lenght-1))+0.0322*cos(8*M_PI*i/(signal_lenght-1));
		out[i] = correction_param * weight*in[i]; 	
	}		
}


/* hannwin

Description: Applyies the Hanning window to the signal. Taken from Lizhe Tan and Jean Jiang, 
             "Chapter 4 - Discrete Fourier Transform and Signal Spectrum" 
             and from https://www.dadisp.com/webhelp/mergedProjects/refman2/FncrefFK/HANNING.htm
Recieves:
	<in> pointer to double: is the array with signal of lenght <signal_lenght>
	<out> pointer to double: is the array after applying the filter with lenght <signal_lenght>.
	<signal_lenght> int: number of data points in the signal
Returns:
	None
*/

void hannwin(double *in, double *out, int signal_lenght){
	int i;
	double weight;
	const double correction_param = 4.0002441565037525; //This coeficient makes sure that a "pure" sinusoidal will have its amplitude preserved
	for(i = 0; i<signal_lenght;i++){
		weight = 0.5-0.5*cos(2*M_PI*i/(signal_lenght-1));
		out[i] = correction_param*weight*in[i]; 	
	}
}

/* triangwin 

Description: Applyies the Triangular window to the signal. Taken from Lizhe Tan and Jean Jiang, 
             "Chapter 4 - Discrete Fourier Transform and Signal Spectrum"
Recieves:
	<in> pointer to double: is the array with signal of lenght <signal_lenght>
	<out> pointer to double: is the array after applying the filter with lenght <signal_lenght>.
	<signal_lenght> int: number of data points in the signal
Returns:
	None

*/

void triangwin(double *in, double *out, int signal_lenght){
	int i;
	double weight;
	const double correction_param = 4.000244170177297; // This coeficient makes sure that a
	                                                   // "pure" sinusoidal will have its amplitude preserved
	for(i = 0; i<signal_lenght; i++){
		weight = 1-fabs(2*i-signal_lenght+1)/(signal_lenght-1);
		out[i] = correction_param*weight*in[i];    //coeficient to correct the amplitude of a sinusoidal
	}
}


/*-------------------------------------------------------------------------------------------------- 

dolphwin

Description: 
	Applyies the Dolph-Chebyshev window to the signal. 
Obs:
	Completely different filter type. It is defined by its frequency response. We need to calculate the inverse dft, normalize it, 
	then apply it to the signal.

	Completely broken.
Recieves:
	<in> pointer to double: is the array with signal of lenght <signal_lenght>
	<out> pointer to double: is the array after applying the filter with lenght <signal_lenght>.
	<signal_lenght> int: number of data points in the signal
	<attenuation> double: attenuation parameter in dB 
Returns:
	None

void dolphwin(double *in, double *out, int signal_lenght, int attenuation){
	int i, j; //iterators
	const int order = signal_lenght - 1;
	double max_weight = 0; //max norm of the weight array. This is used to force a unitary max.
	const double beta = cosh(acosh(pow(10,attenuation/20))/order);
	const double r = 1/pow(10,attenuation/20);
	const double correction_param = 1; //correction value so that sine signal has the correct amplitude 
	const double frequency_den = cosh(order*acosh(beta)); //same denominator in every frequency response point
	const int half = floor(signal_lenght)/2+1;
	double weight[signal_lenght]; //array with weight function
	double frequency_response[signal_lenght]; //array with the frequency response. The Dolph-Chebyshev window is defined in frequency
	double sum, x;
	//calculates the frequency response using chebyshev polynomials   
	for(i = 0; i < signal_lenght; i++){
		x = beta*cos(M_PI*((double)i/order));
		if (fabs(x) <= 1){
			//frequency_response[i] = cos(order*acos(x))/frequency_den;
			frequency_response[i] = cos(order*acos(x));
		} else {
			//frequency_response[i] = cosh(order*acosh(x))/frequency_den;
			frequency_response[i] = cosh(order*acosh(x));
		}
	}
	//inverse fft via direct sum
	for(i = 0; i < signal_lenght; i++){
		sum = 0;
		if (i <= half) {	
			for (j = 0; j < half; j++){
				x = beta*cos(M_PI*((double)j/order));
				if (fabs(x) <= 1) {
					//sum += cos(2*j*acos(x))*cos(2*M_PI*((double)j*i/order));
					sum += cos(j*acos(x))*cos(2*M_PI*((double)j*i/order));
				} else {
					//sum += cosh(2*j*acosh(x))*cos(2*M_PI*((double)j*i/order));
					sum += cosh(j*acosh(x))*cos(2*M_PI*((double)j*i/order));
				}
			}
		} else {
			weight[i] = weight[half-i];
		}
		weight[i] = ((double)1/order)*(1+2*r*sum);
		if (max_weight < weight[i]){
			max_weight = weight[i];
		}
	}
	
	//applyies the filter to the signal
	for (i = 0; i < signal_lenght; i++){
		out[i] = correction_param*weight[i]*in[i]/max_weight;
	}
	return;
}

*/
