/*
 *  FakeGRBredux.cpp
 *  
 *
 *  Created by Sriharsha Pothapragada on 2/4/11.
 *  Copyright 2011 University of Kansas. All rights reserved.
 *
 */

#include<iomanip>
#include<iostream>
#include<math.h>
#include<cmath>
#include<fstream>
#include<gsl/gsl_rng.h>
#include</usr/include/time.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>

using namespace std;

double LinearInterpolate(int x, int x1, int x2, double y1, double y2);

int main(int argc, char *argv[])
{
	fstream readgrbpiece, writegrbfull, readconstants;
	readconstants.open(argv[1]);
	
	int N_spikes = atoi(argv[2]) + 1;
	cout << "The Number of GRB spikes to be incorporated is = " << N_spikes << endl;
	
	double t_o, F_o;
	double tTmp, vTmp, val, dummyAngle;
	int tCtr = 0;
	double factorDilation = atof(argv[3]);
	
		
	readgrbpiece.open(argv[4]);
	while(readgrbpiece)
	{
		for(int i = 0; i < 240; i++)
		{
			readgrbpiece >> tTmp >> vTmp >> val >> dummyAngle;
		}
		tCtr++;
	}
	readgrbpiece.close();
	
	double spike_end = tTmp;
	cout << "Each spike ends at a time of = " << spike_end << endl;
	
	gsl_matrix *F_piece = gsl_matrix_alloc(tCtr, 240);
	gsl_vector *tList = gsl_vector_alloc(tCtr);
	
	gsl_matrix_set_zero(F_piece);
	gsl_vector_set_zero(tList);
	
	readgrbpiece.open(argv[4]);
	int i = 0;
	int j = 0;
	while(readgrbpiece)
	{
		if(j == 240) 
		{
			j = 0;
			i++;
		}
		
		readgrbpiece >> tTmp >> vTmp >> val >> dummyAngle;
		cout << "Inside while Loop doing the " << i+1 << "th iteration"<< "\r";
		gsl_matrix_set(F_piece, i, j, val);
		gsl_vector_set(tList, i, tTmp);
		
		j++;
		
	}
	readgrbpiece.close();
	cout << endl;
	int n_time;
	double dt_grid, tFull, valSpike, F_tmp, F_full;
	
	dt_grid = 1e-4;												//Resolution of time grid
	
	n_time = int((spike_end*factorDilation + 2.0)/dt_grid)  ;  // To calculate the number of grid points to accomodate the full GRB. 
	//n_time = 4100;											// 2.0 seconds is the durations of the longest test LC. 
	//n_time = 150000;
	//n_time = 2000;
	//n_time = 18000;
	//n_time = 140000;
	
	
	int nOld, nNew;
	int tOld, tNew, tGrd;

	/*for (int m = 0; m < 1; m++) 
	{
		for (int n = 0; n < 240; n++)
		{
			cout << setw(10) << gsl_vector_get(tList, m) 
				 << setw(10) << n	   
				 << setw(20) << gsl_matrix_get(F_piece, m, n)
				 << endl;
		}
		cout << endl;
	}*/
	
	
	gsl_matrix *F_sum  = gsl_matrix_alloc(n_time, 240);
	gsl_matrix *F_grid = gsl_matrix_alloc(n_time, 240);	
	
	gsl_matrix_set_zero(F_sum );
	gsl_matrix_set_zero(F_grid);
	
	cout << endl << "Initializations OK " << endl;
	cout << "I have started working... Please be Patient...." << endl; 
	
	for(int ctrSpike = 0; ctrSpike < N_spikes; ctrSpike++)
	{
		cout << "N_spike started = " << ctrSpike+1 << endl;
		
		gsl_matrix_set_zero(F_grid);
		readconstants >> t_o >> F_o;
		
		//cout << "nOld is = " << nOld << endl;
		//cout << "Multiplication factor is = " << F_o << endl;
		
		
		for(int v_step = 0 ; v_step < 240; v_step++)
		{
			nOld = int(t_o/dt_grid);
			for (int i = 0; i < tCtr-1; i++) 
			{
				tFull    = (factorDilation * gsl_vector_get(tList,i)) + t_o;
				nNew     = int(tFull/dt_grid);
				valSpike = F_o * gsl_matrix_get(F_piece, i, v_step);
				//valSpike = F_o;
				
				for (int nGrd = nOld+1; nGrd <= nNew; nGrd++) 
				{
					tGrd     = dt_grid * double(nGrd);
					F_tmp = LinearInterpolate(nGrd, nOld, nNew, 
											  gsl_matrix_get(F_grid, nOld, v_step),
											  valSpike); 
					gsl_matrix_set(F_grid, nGrd, v_step, F_tmp);
					//cout << " nGrd = " << nGrd   << 
					//		" nOld = " << nOld   << 
					//		" nNew = " << nNew   <<
					//		" Ftmp = " << F_tmp  << 
					//		endl; 
					
				}
				
				gsl_matrix_set(F_grid, nNew, v_step, F_tmp);
				nOld = nNew;
			}
		}
		
		for (int m = 0; m < n_time; m++) 
		{
			for (int n = 0; n < 240; n++)
			{
				F_full = gsl_matrix_get(F_sum, m, n) + (gsl_matrix_get(F_grid, m, n));
				gsl_matrix_set(F_sum, m, n, F_full);
			}
		}
	}
	
	writegrbfull.open(argv[5],ios::out);
	if(!writegrbfull) {cout << "Cant Open File to write Full GRB" << endl;}
	
	double v_fake = 1.0;
	for(int i = 0; i < n_time; i++)
	{
		for (int j = 0; j < 240; j++)
		{
			writegrbfull << setw(20) << dt_grid*i << setw(20) << v_fake << setw(20) << gsl_matrix_get(F_sum,i, j) << endl;
			v_fake = v_fake*pow(10.0, 0.025);
		}
		cout << "Writing for time step = " << i << " out of " << n_time << "\r";
		v_fake = 1.0;
		//writegrbfull << endl << endl;
		writegrbfull << endl;
	}
	
	cout << endl;
	
	writegrbfull.close();
	readconstants.close();
	
	gsl_matrix_free(F_sum  );
	gsl_matrix_free(F_piece);
	gsl_matrix_free(F_grid );
	gsl_vector_free(tList  );
	

	return 0;

}

double LinearInterpolate(int x, int x1, int x2, double y1, double y2)
{
	double yIntpl;
	double delxx1  = 1.0*(x  - x1);
	double delx2x1 = 1.0*(x2 - x1);
	
	yIntpl = y1 + (delxx1*(y2 - y1)/delx2x1);
	return yIntpl;
}