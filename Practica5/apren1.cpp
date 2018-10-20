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
//include's de c++
#include <string>
#include <sstream>
using namespace std;

typedef NRVec<int> VecInt, VecInt_O, VecInt_IO;
typedef const NRVec<int>VecInt_I;

typedef NRVec<double> VecDoub, VecDoub_O, VecDoub_IO;
typedef const NRVec<double> VecDoub_I;

typedef NRMat<int> MatInt, MatInt_O, MatInt_IO;
typedef const NRMat<int> MatInt_I;

typedef NRMat<double> MatDoub, MatDoub_O, MatDoub_IO;
typedef const NRMat<double> MatDoub_I;

#define PI 3.14159

void GeneraPatron(VecDoub &z);
void ProdMatVec(MatDoub &sigma,VecDoub &z);

int main(){
	
	int N_IN=4;
	double eta=0.01;

	VecDoub J(N_IN);
	VecDoub z(N_IN);
	MatDoub sigma(N_IN,N_IN);
	double v;

	for(int i=0;i<sigma.nrows();i++){
		for(int j=0;j<sigma.ncols();j++)
			sigma[i][j]=0.309;
		sigma[i][i]+=1;
		}

	srand(time(0));
	for(int i=0;i<J.size();i++)
		J[i]=rand()*1.0/RAND_MAX;

	double t=0;
	while(t<10000){
		v=0;
		GeneraPatron(z);
		ProdMatVec(sigma,z);
		for(int i=0;i<N_IN;i++)
			v+=J[i]*z[i];
		for(int i=0;i<N_IN;i++){
			J[i]+=eta*v*(z[i]-v*J[i]);
			}
		cout << J[0] << " " << J[1] << " " << J[2] << " " << J[3] << endl;
		t++;
		}

	return 0;
	}
	
void GeneraPatron(VecDoub &z){
	
	double p;
	double sigma=1.0;
	double aleat;
	double r;
	int prueba;
	int n=z.size();
	
	for(int i=0;i<n;i++){
		prueba=0;
		while (prueba==0){
			aleat = rand()*1.0/RAND_MAX;
			p=exp(-aleat*aleat/2*sigma*sigma)/(sigma*2*PI);
			r=rand()*1.0/RAND_MAX;
			if(r<p)
				prueba=1;
			}
		z[i]=aleat;
		}
	}

void ProdMatVec(MatDoub &sigma,VecDoub &z){
	
	int n=z.size();
	VecDoub aux(n);
	for(int i=0;i<n;i++){
		aux[i]=0;
		for(int j=0;j<n;j++)
			aux[i]+=sigma[i][j]*z[j];
			}
	
	for(int i=0;i<n;i++)
		z[i]=aux[i];
	}
