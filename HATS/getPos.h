#ifndef GETPOS
#define GETPOS

//Parametros de conexao
//#define IP_SERVER "10.0.92.12"
//#define TCP_PORT 3040
//#define RCV_BUFFER_SIZE 2000
//#define TX_DELAY 0		//Microseconds


//Parametros para Husec_time()
//#define SEC2HUSEC     10000L
//#define HUSEC2NSEC   100000L
//#define ONEDAY        86400L
//#define HUSECS2HRS 36000000L
//#define MIN2HUSEC    600000L
//#define MIN              60L


//Parametros para Shared Memory
//#define BackingFile "HatsRingBuffer"
//#define SemaphoreName "HatsSemaphore"
//#define AccessPerms 0664


//Parametros do Ring Buffer
//#define RINGSIZE 10

//#define DATAFILENAME  "getposition_file.bin"

typedef struct
{
  unsigned long long time_Husec 	;	//Time: Husec
  double time_JD			;	//Time: Julian Date
  double time_Sid			;	//Time: Sideral
  //--------------------------------------------------------
  double pos_tele_alt		;	//Telescope Position: Altitude
  double pos_tele_az		;	//Telescope Position: Azimute
  double pos_tele_ra		;	//Telescope Position: Right Ascension
  double pos_tele_dec		;	//Telescope Position: Declination
  //---------------------------------------------------------
  double rate_ObjId_ra		;	//Tracking rate: Right Ascension diff from sidereal rate
  double rate_ObjId_dec		;	//Tracking rate: Declination diff from sidereal rate
  //---------------------------------------------------------
  int object			;	//Operation mode
  int opmode			;  	//Operation mode

} pos_data_type ;


//----------Begin--Of--Javascripts--------------------------------------------

int print_usage()
{
	printf("\n\n  getPos - Get Position program.\n\n\n");
	
	
	printf("SYNOPSIS\n\n");
	
	printf(" getPos [OPTIONS]\n\n");
	
	printf("DESCRIPTION\n\n");
	
	printf(" This program performs, in a infinite loop, data aquisition from TheSkyX server collecting data     \n");
	printf(" about current position of the telescope in equatorial and horizontal coordinates as well as their  \n");
	printf(" respective tracking rates. Collects also time information for use in timestamps in julian date,    \n");
	printf(" hundred of micro seconds (Husec) and sideral time. And collects the operation mode of telescope,   \n");
	printf(" by accessing its shared memory. The data structure is ordenated as described below:		    \n");
	printf(" ________________________________________________________________________________________________________ \n");
	printf(" Husec| Julian Date | Sideral Time | Altitude | Azimute | RA | DEC | RA rate | DEC rate | Object | OpMode \n");
	printf("   |	    |		   |		|	   |	  |     | 	|	  |	    |	     |	  \n");
	printf("   |	  double         double	      double	   |	  |	|     double	double 	   int      int	  \n");				   
	printf("   |		   	  			double	  |   double  					  \n");
	printf("  unsigned long long					double					          \n");
	printf(" ________________________________________________________________________________________________________ \n");
	printf(" Was built to working in background, but for tests, is possible some options as described below.  \n\n");
	
	printf("OPTIONS\n\n");
	
	printf(" -h, --help             Print this help and exit.\n");
	printf(" -v, --verbose          Verbose option. Basic output informations.\n");
	printf(" -c, --clock		Time measure of loops to estimate ring buffer size.\n");
	printf(" -d, --daemon		Run getPos as a daemon executing in background.\n");
	printf("			This option have safe instruction using SIGTERM, by\n");
	printf("			closing all descriptors and memories before be killed.\n");
	printf("			This option must be used without another option.\n");
	printf(" --debug                Full output information from getPos and log.\n");
	printf(" --version              Show getPos's version and exit.\n\n");
	
	printf("AUTHOR\n\n");
	
	printf(" Written by Tiago Giorgetti and Guillermo de Castro\n\n");
	
	printf("REPORTING BUGS\n\n");
	
	printf(" Please send an email to tgiorgetti@gmail.com or guigue@craam.mackenzie.br\n\n");
	
	printf("COPYRIGHT\n\n");
	
	printf(" Copyright © 2021 Free Software Foundation, Inc. License GPLv3+: GNU GPL version 3  or later\n");
	printf(" <https://gnu.org/licenses/gpl.html>.\n");
	printf(" This  is free software: you are free to change and redistribute it.\n");
	printf(" There is NO WARRANTY, to the extent permitted by law.\n\n");
        
	printf("SEE ALSO\n\n");
	
	printf("Full documentation at: <http://github.com>\n\n");
	
	exit(0);
}




















#endif
