#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <fstream>
#include <iostream>
#include <iostream>
#include <iomanip>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <malloc.h>
#include "nrutilc.h"
#include <time.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
//include's de c++
#include <string>
#include <sstream>

#define PI 3.14159

using namespace std;

typedef NRVec<int> VecInt, VecInt_O, VecInt_IO;
typedef const NRVec<int>VecInt_I;

typedef NRVec<double> VecDoub, VecDoub_O, VecDoub_IO;
typedef const NRVec<double> VecDoub_I;

typedef NRMat<int> MatInt, MatInt_O, MatInt_IO;
typedef const NRMat<int> MatInt_I;

typedef NRMat<double> MatDoub, MatDoub_O, MatDoub_IO;
typedef const NRMat<double> MatDoub_I;

void GeneraPatron(VecDoub &z);
double Prod(VecDoub &z,MatDoub &sigma);

int main(){
	
	int N_IN=4;
	/*printf("Escriba el numero de neuronas de entrada.\n");
	scanf("%d",&N_IN);*/
	double eta=0.01;

	VecDoub J(N_IN);
	VecDoub z(N_IN);
	MatDoub sigma(N_IN,N_IN);
	double v;

	srand(time(0));
	for(int i=0;i<J.size();i++)
		J[i]=rand()*1.0/RAND_MAX;
	GeneraPatron(z);
	
	double t=0;
	while(t<10000){
		v=0;
		for(int i=0;i<N_IN;i++)
			v+=J[i]*z[i];
		for(int i=0;i<N_IN;i++){
			J[i]+=eta*v*(z[i]-v*J[i]);
			}
		t++;
	}
	for(int i=0;i<N_IN;i++)
		cout << J[i] << endl;		
	return 0;
	}
	
void GeneraPatron(VecDoub &z){
    
    MatDoub inv_sigma(4,4);
	for(int i=0;i<inv_sigma.nrows();i++){
		for(int j=0;j<inv_sigma.ncols();j++){
			inv_sigma[i][j]=-0.2;	
			}
		inv_sigma[i][i]+=1;	
		}
				
    double det=5;
    int patron = 0;
    while(patron == 0){
		
		for(int i=0;i<z.size();i++){
			z[i]=rand()*1.0/RAND_MAX;
			cout << "z[i]=" << z[i] << endl;
			}
		double p;
		p = exp(-0.5*Prod(z,inv_sigma))/(4*PI*PI*sqrt(det));
		cout << "p=" << p << endl;
		double aleatorio = rand()*1.0/RAND_MAX;
		if(aleatorio < p)
			patron=1;
		}
	}
	
double Prod(VecDoub &z,MatDoub &sigma){
	
	double prod=0;
	VecDoub aux(z.size());
	
	for(int k=0;k<z.size();k++){
		for(int i=0;i<sigma.nrows();i++)
			for(int j=0;j<sigma.ncols();j++)
				aux[i]+=sigma[i][j]*z[i];
		prod+=aux[k]*z[k];
		}

	return prod;
	}
