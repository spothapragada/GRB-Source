/*
 *  BinData.cpp
 *  
 *
 *  Created by Sriharsha Pothapragada on 10/13/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include<iomanip>
#include<iostream>
#include<math.h>
#include<cmath>
#include<fstream>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<stdlib.h>

using namespace std;

int main(int argc, char *argv[])
{
	
	fstream readdata,readtimes,writedata;

	double t, F, E;
	double dE = 10.0;

	double F_bol_t = 0.0;
	double bol_norm = 0.0;
	double win_left = atof(argv[2]);
	double win_right = atof(argv[3]);
	
	readdata.open(argv[1]);						//Read the synthetic GRB file FakeGRB_sharp.txt
  	writedata.open(argv[4],ios::out);			//Open the Bolometric Flux File t vs F_bol(t)
	
	while(readdata)
	{
		double F_bol_t = 0.0;
		for(int i = 0; i < 240; i++)
		{
			readdata >> t >> E >> F;
			if(E >= win_left && E <= win_right)
			{
				//F_bol_t = F_bol_t + F*pow(dE,0.025*i); 
				F_bol_t = F_bol_t + (F*(E*pow(10.0,0.025) - E));
			}
			//writedata  << setw(10) << t << "\t" << F_bol_t <<endl;
			
			 
		}
		writedata  << setw(10) << t << "\t" << F_bol_t << endl; 
// 		writedata << endl ;
		if(F_bol_t > bol_norm) {bol_norm = F_bol_t;}
		cout << "Computing Bolometric Flux for t = " << t << "\r";
	}

	cout << "Normalization for Bolometric Flux is = " << bol_norm << endl;
	cout << endl << "Finished Computing Bolometric Flux!!" << endl;
	//cout << endl;

	readdata.close();
	writedata.close();

 	readdata.open(argv[4]);					// Open the Bolometric Flux file
	
	double F_bol_read, t_bol_read;
	double F_bol_cum = 0.0;

 	//double threshold = atof(argv[5]) ; //  THRESHOLD   // Obtain the threshold for binning

	double t_left_bin = 0.0, t_right_bin;
	double t_step = 1e-4; //Same as the grid resolution for the FakeGRB  
	
	int num_of_times  = 0;

	while(readdata)
	{
		readdata >> t_bol_read >> F_bol_read;			// Get the number of time ticks in the bolometric flux

		num_of_times++;
	}
	
	readdata.close();

	cout << endl << "Number of times are = " << num_of_times << endl;
	double t_bol[num_of_times], F_bol[num_of_times];

	for(int i = 0 ;i < num_of_times; i++)

	{
		t_bol[i] = 0.0;
		F_bol[i] = 0.0;
	}

	
	readdata.open(argv[4]);

	double temp1, fluenceTotal = 0;
	for(int i = 0; i < num_of_times ; i++)
	{

		//readdata >> t_bol[i] >> temp1;				// Read in the times and the Normalised bolometric fluxes into thir arrays
 		//F_bol[i] = temp1/bol_norm;
		readdata >> t_bol[i] >> F_bol[i];
		fluenceTotal = fluenceTotal + F_bol[i] * t_step;
	}

	//double avgFBol = fluenceTotal / (num_of_times * t_step);
	
	readdata.close();

	double threshold = atof(argv[5])*fluenceTotal ; //  THRESHOLD   // Obtain the threshold for binning fraction of bolometric fluence
	
	//double threshold = atof(argv[5]) * avgFBol;
	
	readdata.open(argv[4]);							// Reopen the bolometric flux file
  	writedata.open(argv[6],ios::out);				// Binned Bolometric Flux

	int bin_ctr = 0;
	for(int i = 0;i < num_of_times; i++)
	{
		
		cout << "Binning Normalized Bolometric Flux.Please be Patient" << "\r";

		//F_bol_cum = F_bol_cum + F_bol[i];
		F_bol_cum = F_bol_cum + F_bol[i] * t_step;
		//F_bol_cum = F_bol_cum + F_bol[i]/(num_of_times);
		
		if(F_bol_cum >= threshold)
		{
			t_right_bin = t_bol[i+1];
			writedata << t_left_bin << "\t" << t_right_bin << "\t" << 0.5*(t_left_bin + t_right_bin) <<  "\t" << (F_bol_cum*t_step)/((t_right_bin - t_left_bin)*bol_norm) << endl;
			//writedata << t_left_bin << "\t" << t_right_bin << "\t" << 0.5*(t_left_bin + t_right_bin) <<  "\t" << (F_bol_cum) << endl;
			//bin_ctr++;
			
			F_bol_cum = 0.0;
			t_left_bin = t_right_bin;
		}

	}

	cout << endl << "Done Binning!" << endl;

	readdata.close();
	writedata.close();

 	readdata.open(argv[1]);

	double t_read, E_read, F_read;
	//double Ener[240];
	//double F_spec[9000][240];
	//double time[9000];
	
	int n_time = num_of_times - 2;
	
	gsl_vector *Ener   = gsl_vector_alloc(240        );
	gsl_vector *time   = gsl_vector_alloc(n_time     );
	gsl_matrix *F_spec = gsl_matrix_alloc(n_time, 240);
	
	gsl_matrix_set_zero(F_spec);
	gsl_vector_set_zero(time);
	gsl_vector_set_zero(Ener);
	
	
	
	for(int i = 0; i < n_time; i++)
	{
		for(int j = 0; j < 240 ; j++)
		{
			//readdata >> time[i] >> Ener[j] >> F_spec[i][j];	
			readdata >> t >> E >> F;
			gsl_matrix_set(F_spec, i, j, F);
			gsl_vector_set(time, i, t);
			gsl_vector_set(Ener, j, E);
			
		}
	}
	readdata.close();

	readdata.open(argv[6]);

	double t_l_ctr,t_r_ctr,t_mid_ctr, F_ctr; 
	while(readdata)
	{
		readdata >> t_l_ctr >> t_r_ctr >> t_mid_ctr >>  F_ctr;
		bin_ctr++;
	}

	cout << "Number of Bins is = " << bin_ctr << endl;
	readdata.close();

	
	gsl_matrix *F_avg = gsl_matrix_alloc(bin_ctr, 240);
	gsl_matrix_set_zero(F_avg);
	
	
	//double F_avg[bin_ctr][240];
	/*for(int l = 0; l < bin_ctr; l++)
		{
			for(int j = 0; j < 240; j++)
			{
				F_avg[l][j] = 0.0;
			}
		}
	 */
  	readtimes.open(argv[6]);

 	writedata.open(argv[7], ios::out);

	
	int bin_mark = 0;
	double t_mid, F_bin;
 	while(readtimes)
 	{
		readtimes >> t_left_bin >> t_right_bin >> t_mid >> F_bin;

		for(int i = 0; i < n_time; i++)
		{
			if(gsl_vector_get(time,i) >= t_left_bin && gsl_vector_get(time,i) <= t_right_bin)
			{	
				for(int j = 0; j < 240; j++)
				{
					//F_avg[bin_mark][j] = F_avg[bin_mark][j] + F_spec[i][j]*t_step/(t_right_bin - t_left_bin);
					double temp = gsl_matrix_get(F_avg, bin_mark,j) 
									+ gsl_matrix_get(F_spec, i,j)*t_step/(t_right_bin - t_left_bin);
					
					gsl_matrix_set(F_avg,bin_mark,j,temp);
					
				}
			}
		}
		
		for(int k = 0; k < 240 ;k++)
		{
			writedata << t_left_bin << "\t" << t_right_bin << "\t" << gsl_vector_get(Ener,k) << "\t" << gsl_matrix_get(F_avg,bin_mark,k) << endl; 
		}
		writedata << endl;
		cout << "Wrote Spectra for Bin # = " << bin_mark+1 << "\r";		
		bin_mark++;		
		
		 		
 	}
	cout << endl;
	readtimes.close();
	writedata.close();
	
	gsl_matrix_free(F_spec);
	gsl_matrix_free(F_avg);
	gsl_vector_free(Ener);
	gsl_vector_free(time);
	return 0;
}


