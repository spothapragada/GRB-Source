#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define _USESTDVECTOR_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <NR3/nr3.h>
#include <NR3/gaussj.h>
#include <NR3/fitmrq.h>


using namespace std;

double Epiv = 1000.0;
double Gamma = 100.0;

void BandFn(const Doub En, VecDoub_I &a, Doub &y, VecDoub_O &dyda) ;

int main(int argc, char *argv[])
{
	ifstream readdata,readdata2;
    
 	readdata.open(argv[1]);
	readdata2.open(argv[2]);

	fstream tabledata;
 	tabledata.open(argv[3],ios::out);

	double time_left, time_right,Ener, flux;

	double win_left    = atof(argv[4]);
	double win_right   = atof(argv[5]);
	double err_frac    = atof(argv[6]);
	
	int spectra_ctr = 0;
    double fitparam_init[4] = {1.0E28, 2.0, -0.5, 5000.0};
	
	VecDoub guess(fitparam_init,fitparam_init + 4);
	
	//int n = 240;
	
	//VecDoub E (n, 0.0);
	//VecDoub fE(n, 0.0);
	//VecDoub error(n, 0.0);
	//VecDoub zeros(n, 0.0);
	
	
	vector<double> E;
	vector<double> fE;
	vector<double> error;
	
	cout << "I am here" << endl;
	
	while(readdata)
	{
		//double fitparam_init[4] = {1E30, 1.0/3.0, -2.5, 5000.0};
		
		for(int ctr = 0; ctr < 240; ctr ++)
		{
			
				
			//readdata >> time_left >> time_right >> E[ctr] >> fE[ctr] ;
			readdata >> time_left >> time_right >> Ener >> flux ;
			
			//E[ctr] = Ener;
			//fE[ctr] = Ener*flux;
						
			//error[ctr] = err_frac*fE[ctr];
			
			if (Ener >= win_left && Ener <= win_right) 
			{
				E.push_back(Ener);
				fE.push_back(flux);
				error.push_back(err_frac*flux);
			}
		}
			
		
		int n = E.size();
		
		Fitmrq GRBmrq(E,fE,error,guess,BandFn,1e-7);
		
		GRBmrq.fit();
		
		
		Doub funcvalue = GRBmrq.a[0]* pow(GRBmrq.a[3]/Epiv, GRBmrq.a[1])*exp((GRBmrq.a[3]*(GRBmrq.a[2]-GRBmrq.a[1]))/GRBmrq.a[3]);
		
		
		double F_bol;
		
		readdata2 >> F_bol >> F_bol >> F_bol >> F_bol;
	
		{
		
			
			tabledata << setprecision(8) << setw(15) << time_left << setw(15) << time_right << setw(15) << GRBmrq.a[0] << setw(15) << GRBmrq.a[1] << setw(15) << GRBmrq.a[2] << setw(15) << GRBmrq.a[3] << setw(15) << funcvalue << setw(15) << GRBmrq.chisq << setw(15) << F_bol << endl;

			cout << "Fitting Spectra # = " << setw(5) << spectra_ctr << "\n";
			spectra_ctr++;
		}
		
	
		
		guess = GRBmrq.a; //feedback the guesses
		
		//E     = zeros;
		//fE    = zeros;
		//error = zeros;
		
		E.clear();
		fE.clear();
		error.clear();
		

	}
	cout << endl; 
	readdata.close();
	tabledata.close();
	
	return 0;
		
		

}

void BandFn(const Doub En, VecDoub_I &a, Doub &y, VecDoub_O &dyda) {
	
	Doub a_sig = 3.0;
	Doub A = a[0];	//A
	Doub B = a[1];	//Alpha
	Doub C = a[2];	//Beta
	Doub D = a[3];	//Ebreak
	
	
	Doub f1E  = pow(En/Epiv, B)*exp((En*(C-B))/D);
	Doub f2E  = pow(D/Epiv, B-C)* pow(En/Epiv, C) * exp(C-B);
	Doub sig1 = 0.5*(1.0 - tanh(a_sig*(En - D)));
	Doub sig2 = 0.5*(1.0 + tanh(a_sig*(En - D)));
	
	y = (A*f1E*sig1) + (A*f2E*sig2);
	
	
	dyda[0] = (f1E*sig1) + (f2E*sig2);
	dyda[1] = (A*f1E*sig1*(log(En/Epiv)-(En/D))) + (A*f2E*sig2*(log(En/Epiv)-1.0));      
	dyda[2] = (A*f1E*sig1*(En/D)) + (A*f2E*sig2*(log(En/D)+1.0));
	dyda[3] = (A*f1E*((sig1*En*(B-C)/(D*D))+0.5*a_sig*pow(cosh(En-D),-2.0)) + (A*f2E*(((B-C)*sig2/D)-(0.5*a_sig*pow(cosh(En-D),-2.0)))));
	
}


