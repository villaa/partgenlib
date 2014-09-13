/*=========================partgen.h==================================== 
   
      PROGRAMMER:  Anthony Villano/Allison Kennedy 09/12/14

      UPDATES:      
       

      PURPOSE: a library to store some functions for generating particle
               energy/angle distributions.
              
======================================================================*/
#include "Math/Minimizer.h"
#include "TRandom.h"
#include <TFile.h>
#include <TTree.h>
#include <Riostream.h> 
#include <iostream>
#include <fstream>
#include "TH1.h"
#include <algorithm>
#include <map>
#include <vector>
#include <string>
#include <stdio.h>



#ifndef partgen_h
#define partgen_h 1

std::pair<double *, double *> EGen_dist(const char* filename, const char* histname);
double EGen_rand(double* bins, double* intbinamp, TRandom *frand);
bool genFlatDist(vector<double> &linespec, TRandom *rand, vector<double> refplane);
bool genCosDist(vector<double> &linespec, TRandom *rand, vector<double> refplane);

#endif
