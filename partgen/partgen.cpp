/*=========================partgen.h==================================== 
   
      PROGRAMMER:  Anthony Villano/Allison Kennedy 09/12/14

      UPDATES:      
       

      PURPOSE: a library to store some functions for generating particle
               energy/angle distributions.
              
              
======================================================================*/

#define PI 3.141592653589793
#include "partgen.h"

//Functions to read histogram energy distribution [ABK]
std::pair<double *, double *> EGen_dist(const char* filename, const char* histname)
	{

	Double_t Energy=0.0;

	const char * suffix;
	string rootfiletest=".root";
	string txtfiletest=".txt";
  	suffix = strstr (filename,".");

///////////////////////////////// For text files /////////////////////////////////

	if (string(suffix)==txtfiletest) 
		{
		ifstream in_c1;
		in_c1.open(filename);
		Double_t c1_temp;
		vector<double> c1;
		int i=0;

		while (!in_c1.eof())
			{
			in_c1 >> c1_temp;
			c1.push_back(c1_temp);
			i++;
			//cout << "i: " << i << endl;
			}
		Int_t nevents=i;
		//Double_t *c1_arr=&c1[0];
		Double_t maxE=*std::max_element(c1.begin(),c1.end());
		Double_t minE=*std::min_element(c1.begin(),c1.end());
		Int_t nbins=2000;
		vector <double> bins;
		vector <double> binamp;
		vector <double> intbinamp;

		Double_t binwidth=(maxE-minE)/(double)nbins;

		//cout << "minE: " << minE << " maxE: " << maxE << " binwidth: " << binwidth << endl;
		for (int i=0;i<nbins;i++)
			{
			bins.push_back(minE+i*binwidth); //define minimum value for each bin
			//cout << "bins: " << bins[i] << endl;
			binamp.push_back(0.0);	// define the amplitude for each bin value
			intbinamp.push_back(0.0);	// define the integrate amplitude for each bin value
			}

		//cout << "GOT TO HERE" << endl;

		for (int i=0;i<nevents;i++)
			{
			int j=(int) floor((c1[i]-minE)/(double) binwidth); // for each event in a bind width, increase binamp by 1
			//cout << "j: " << (c1[i]-minE)/(double) binwidth << endl;
			binamp[j]=binamp[j]+1; //((double) nevents);
			//cout << "j: " << j << endl;
			}
	
		for (int i=0; i<nbins; i++)
			{
			//cout << "binamp: " << binamp[i] << endl;
			if (i==0) {intbinamp[i]=binamp[i];}
			else
				{
				intbinamp[i]=intbinamp[i-1]+binamp[i]; //integrated bin amp is sum of all previous amplitudes
				}
			//cout << "intbinamp: " << intbinamp[i] << endl;
			}
		}


///////////////////////////////// For Root files /////////////////////////////////

		vector<double> bins;
		vector<double> binamp;
		vector<double> intbinamp;

	if (string(suffix)==rootfiletest) 
		{

		TFile *f = new TFile(filename);  
        	TH1F *h1 = (TH1F*)f->FindObjectAny(histname); // This will use the histogram you enetered at the top.


		Int_t nbins=h1->GetNbinsX();
		Double_t maxE=h1->GetBinLowEdge(nbins);
		Double_t minE=h1->GetBinLowEdge(1);

		bins.clear();
		binamp.clear();
		intbinamp.clear();
		
		bins.push_back(nbins);   	//Store bin number as first part of arrays
		intbinamp.push_back(nbins);


		//cout << "minE: " << minE << " maxE: " << maxE << " nbins: " << nbins << endl;
		for (int i=0;i<nbins;i++)
			{
			bins.push_back(h1->GetBinLowEdge(i)); //define minimum value for each bin
			binamp.push_back((h1->GetBinContent(i)));	// define the amplitude for each bin value
			intbinamp.push_back(0.0);	// define the integrate amplitude for each bin value
			//cout << "bins: " << bins[i] << " " << binamp[i] << endl;
			}

	
		for (int i=1; i<=nbins; i++)
			{
			//cout << "binamp: " << binamp[i] << endl;
			if (i==1) {intbinamp[i]=binamp[i];}
			else
				{
				intbinamp[i]=intbinamp[i-1]+binamp[i]; //integrated bin amp is sum of all previous amplitudes
				}
			//cout << "intbinamp: " << intbinamp[i] << endl;
			}
		//cout << ""; //Not sure why this HAS to be here, but there ya go...
		f->Close();
		}

///////////////////////////////// Output Bin Values  and Cumulative Bin Amplitude /////////////////////////////////

    
	double * intbinamp_Arr = (double*) malloc(intbinamp.size()*sizeof(double));
	double * bins_Arr = (double*) malloc(bins.size()*sizeof(double));

	for(int i=0;i<(int)bins.size();i++){
          *(intbinamp_Arr + i) = intbinamp[i];
	  *(bins_Arr + i) = bins[i];
	}

	return std::make_pair(bins_Arr,intbinamp_Arr);
	}
	



double EGen_rand(double* bins, double* intbinamp, TRandom *frand)
	{

	Double_t Energy=0.0;
	        //Use the root random generator, note it's not a great one like the M. Twister, 
		//but I think it's OK for this application. 
		//srand (time(0);
		//srand (time(0)*randseed);
		//Double_t rand1=((double) rand() / (RAND_MAX));
		Double_t rand1 = frand->Uniform();
		int nbins=bins[0];  			//Number of Bins in the original histogram.
		//cout << "Random: " << rand1 <<   " " << nbins << endl;
		int eventnum=intbinamp[nbins-2];
		//cout << eventnum << endl;


                Double_t newEmax,newEmin,deltaE;
		for( int i=1; i<=nbins; i++)
			{
			//cout << "intbinamp: " << intbinamp[i] << endl;
			if (rand1> intbinamp[i]/((double) eventnum) && rand1<intbinamp[i+1]/((double) eventnum)) //Find new random within selcted bin
				{
				newEmax=bins[i+1]; //Energy will rest between these two bins
				newEmin=bins[i];
				Double_t binwidth=newEmax-newEmin; // Find the width between the two bins
				deltaE=(rand1*((double) eventnum)-intbinamp[i])*binwidth/(intbinamp[i+1]-intbinamp[i]);
			//cout << "Emax " << newEmax << " Emin " << newEmin << " Delta " << deltaE << endl;
		//cout << "rand1 " << rand1 << " i " << intbinamp[i] << " i+1 " << intbinamp[i+1] << endl;
				}
			}

		Energy=newEmin+deltaE;

	return Energy;
	}

bool genFlatDist(vector<double> &linespec, TRandom *rand, vector<double> refplane)
{

  //random number should already be seeded. 
  linespec.clear();
  if(refplane.size() != 6)
  {
     cerr << "No location specification for random generation" << endl;
     return 0;
  }
  double costhet;
  costhet=rand->Uniform();  //only generate over the upper hemisphere
  linespec.push_back(acos(costhet));
  linespec.push_back(2.0*PI*rand->Uniform());
  linespec.push_back((refplane[2]-refplane[1])*rand->Uniform()+refplane[1]);
  linespec.push_back((refplane[4]-refplane[3])*rand->Uniform()+refplane[3]);
  linespec.push_back((refplane[5])*rand->Uniform()+(refplane[0]-(refplane[5]/2.0)));
  
  return 1;

}
bool genCosDist(vector<double> &linespec, TRandom *rand, vector<double> refplane)
{

  //random number should already be seeded. 
  linespec.clear();
  if(refplane.size() != 6)
  {
     cerr << "No location specification for random generation" << endl;
     return 0;
  }
  double costhet,prob;
  do
  {
    prob = rand->Uniform();
    costhet=rand->Uniform();  //only generate over the upper hemisphere
  }
  while(prob >= pow(costhet,1.3)/0.9186);
    
  
  linespec.push_back(acos(costhet));
  linespec.push_back(2.0*PI*rand->Uniform());
  linespec.push_back((refplane[2]-refplane[1])*rand->Uniform()+refplane[1]);
  linespec.push_back((refplane[4]-refplane[3])*rand->Uniform()+refplane[3]);
  linespec.push_back((refplane[5])*rand->Uniform()+(refplane[0]-(refplane[5]/2.0)));
  
  return 1;

}
