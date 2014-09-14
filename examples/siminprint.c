/*===========================siminprint.c=============================== 
   
      PROGRAMMER:  Anthony Villano 09/12/14

      UPDATES:

      PURPOSE:  A program to utilize the partgen library insipred by
		Allison Kennedy.  This program reads in a root file on
		input, which contains a single histogram named
		"returnedhist." The program then prints the generated
		lines which specify thrown particles with the energy
		spectrum that was read in, and a uniform angular
		distribution generated on a hemisphere.
              
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
           "  -g, --geomfile      <filename>    name plane geometry file\n"
           "  -h, --help                        display this usage information\n"
           "  -n, --ngen                        number of generated events\n"
           "  -u, --normal        <direction>   direction of normal to rad emission\n"
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
   string geomfile="null";
   string normal="-z";
   uint verbosity=0;
   uint quietness=0;
   int  ntrig=10;
   bool quiet=false;

   //***********Begin get input*********************//
 
   
   const struct option longopts[] =
   {
     {"histfile",  required_argument,  0, 'd'},
     {"geomfile",  required_argument,  0, 'g'},
     {"help",      no_argument,        0, 'h'},
     {"ngen",     required_argument,  0, 'n'},
     {"normal",   required_argument,  0, 'u'},
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
     iarg = getopt_long(argc, argv, "+d:g:hn:u:v::V", longopts, &index);

     switch (iarg)
     {


       case 'd':
         histfile = optarg;
         break;

       case 'g':
         geomfile = optarg;
         break;

       case 'h':
         print_usage(stdout,0);
         break;

       case 'n':
         ntrig = atoi(optarg);
	 break;

       case 'u':
         normal = optarg;
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
   
   //read the number of geometrical planes for generation
   double zmid,xmin,xmax,ymin,ymax,zthk;
    std::ifstream fileStream(geomfile.c_str(),std::ios::in);
    if(!fileStream){
      std::cerr << "ERROR! (simprint): geometry file not found." << std::endl;
      exit(1);
    }

    vector<vector<double> > planes;
    fileStream >> zmid >> xmin >> xmax >> ymin >> ymax >> zthk;
    while(!fileStream.eof()){

      vector<double> oneplane;
      oneplane.push_back(zmid);
      oneplane.push_back(xmin);
      oneplane.push_back(xmax);
      oneplane.push_back(ymin);
      oneplane.push_back(ymax);
      oneplane.push_back(zthk);

      planes.push_back(oneplane);

      fileStream >> zmid >> xmin >> xmax >> ymin >> ymax >> zthk;
    }

   //encode the energy histogram
   std::pair<double *,double *> dist = EGen_dist(histfile.c_str(),"returnedhist");

   //seed a random number
   TRandom frand(time(0));

   //set up necessary vars
   double x,y,z,px,py,pz,E;
   uint evnum=0;
   int partnum = 22; //gammas for now

   string space="\t";
   std::cout.precision(5);
   //std::cout << std::dec << std::setw(20) << std::setfill('0');
   std::cout << "\t\t" << 0 << "\t" << 0 << "\t" << 0 << endl;
   for(int i=0;i<planes.size();i++){
     for(int j=0;j<ntrig;j++){

       vector<double> linespec;
       E = EGen_rand(dist.first,dist.second,&frand);
       //FIXME
       E = sqrt(E*E); //make positive, sometimes is negative not sure why
       genFlatDist(linespec,&frand,planes[i]);
       x = linespec[2];
       y = linespec[3];
       z = linespec[4];

       //use the normal to decide the px,py,pz normal is always along thet=0
       if(normal=="-z"){
         pz = -cos(linespec[0]);
	 px = sin(linespec[0])*cos(linespec[1]);
	 py = sin(linespec[0])*sin(linespec[1]);
       }
       else if(normal=="+z"){
         pz = cos(linespec[0]);
	 px = sin(linespec[0])*cos(linespec[1]);
	 py = sin(linespec[0])*sin(linespec[1]);
       }
       else if(normal=="-x"){
         px = -cos(linespec[0]);
	 py = sin(linespec[0])*cos(linespec[1]);
	 pz = sin(linespec[0])*sin(linespec[1]);
       }
       else if(normal=="+x"){
         px = cos(linespec[0]);
	 py = sin(linespec[0])*cos(linespec[1]);
	 pz = sin(linespec[0])*sin(linespec[1]);
       }
       else if(normal=="-y"){
         py = -cos(linespec[0]);
	 pz = sin(linespec[0])*cos(linespec[1]);
	 px = sin(linespec[0])*sin(linespec[1]);
       }
       else if(normal=="+y"){
         py = cos(linespec[0]);
	 pz = sin(linespec[0])*cos(linespec[1]);
	 px = sin(linespec[0])*sin(linespec[1]);
       }
       else{
         pz = -cos(linespec[0]);
	 px = sin(linespec[0])*cos(linespec[1]);
	 py = sin(linespec[0])*sin(linespec[1]);
       }


       std::cout << space 
          << std::setw(10) << (evnum+1) << space 
	  << std::setw(10) << partnum << space 
          << std::setw(10) << E << space
	  << std::setw(10) << x << space
          << std::setw(10) << y << space 
	  << std::setw(10) << z << space 
	  << std::setw(10) << px << space 
	  << std::setw(10) << py << space 
	  << std::setw(10) << pz << std::endl;

       evnum++;

     }
   }

   return 0;
}
