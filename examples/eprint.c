/*===========================eprint.c=================================== 
   
      PROGRAMMER:  Anthony Villano 09/12/14

      UPDATES:

      PURPOSE:  A first program to utilize the partgen library insipred by
		Allison Kennedy.  This program reads in a root file on
		input, which contains a single histogram named
		"returnedhist." The program then prints the generated
		energies in the distribution specified by the histogram.
              
      INPUT:      

      OUTPUT:   
              
======================================================================*/
#include <iostream> 
#include <iomanip> 
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <string> 
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

//root libraries
#include "TRandom.h"

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

//partgen library
#include "partgen.h"

//The name of this program. 
const char* program_name;

//Prints usage information for this program to STREAM (typically
//stdout or stderr), and exit the program with EXIT_CODE.  Does not
//return. 

void print_usage (FILE* stream, int exit_code)
{
  fprintf (stream, "Usage:  %s options [ inputfile(s) ]\n", program_name);
  fprintf (stream,
	   //"\n"
           "  -d, --histfile      <filename>    name the input hist file file \n"
           "  -h, --help                        display this usage information\n"
           "  -n, --ngen                        number of generated events\n"
           "  -v, --verbose       <level>       Print verbose messages at level <level>\n"
           "  -V, --version                     print version and exit\n");

  exit (exit_code);
}

#define no_argument 0
#define required_argument 1 
#define optional_argument 2

using namespace std;

int main(int argc, char *argv[])
{
   program_name = argv[0];

   string rep; 
   string histfile="muenergies_1.root";
   uint verbosity=0;
   uint quietness=0;
   int  ntrig=10;
   bool quiet=false;

   //***********Begin get input*********************//
 
   
   const struct option longopts[] =
   {
     {"histfile",  required_argument,  0, 'd'},
     {"help",      no_argument,        0, 'h'},
     {"ngen",     required_argument,  0, 'n'},
     {"verbose",   optional_argument,  0, 'v'},
     {"version",   no_argument,        0, 'V'},
     {0,0,0,0},
   };

   int index;
   int iarg=0;

   //turn off getopt error message
   opterr=1; 

   while(iarg != -1)
   {
     iarg = getopt_long(argc, argv, "+d:hn:v::V", longopts, &index);

     switch (iarg)
     {


       case 'd':
         histfile = optarg;
         break;

       case 'h':
         print_usage(stdout,0);
         break;

       case 'n':
         ntrig = atoi(optarg);
	 break;

       case 'v':
         if(optarg)
           verbosity = atoi(optarg);
	 else
	   verbosity = 1;
         break;

       case 'V':
         printf("Version: %s\n", __GIT_VERSION);
         break;

       //case for long only options
       case 0:
         break;

       case '?':
         print_usage(stderr,1);
         break;
     }
   } 

   //***********End get input*********************//

   //encode the histogram
   std::pair<double *,double *> dist = EGen_dist(histfile.c_str(),"returnedhist");

   //seed a random number
   TRandom frand(time(0));

   for(int i=0;i<ntrig;i++){

     std::cout << EGen_rand(dist.first,dist.second,&frand) << std::endl;

   }

   return 0;
}
