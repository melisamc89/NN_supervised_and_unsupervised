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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
//include's de c++
#include <string>
#include <sstream>
using namespace std;

template<class T>
inline string to_string(const T& t){
    stringstream ss;
    ss << t;
    return ss.str();
}

typedef NRVec<int> VecInt, VecInt_O, VecInt_IO;
typedef const NRVec<int>VecInt_I;

typedef NRVec<double> VecDoub, VecDoub_O, VecDoub_IO;
typedef const NRVec<double> VecDoub_I;

typedef NRMat<int> MatInt, MatInt_O, MatInt_IO;
typedef const NRMat<int> MatInt_I;

typedef NRMat<double> MatDoub, MatDoub_O, MatDoub_IO;
typedef const NRMat<double> MatDoub_I;

void GeneraPatronUniforme(VecDoub &z);
double ProbPolar(VecDoub &z);
void InitWeigh(MatDoub &matrix);
void CalculoH(MatDoub &J,VecDoub &z,VecDoub &h);
int MaxVectorPosition(VecDoub &vec);
double Lambda(int i1,int i2,double sigma);
void CalculoDeltaJ(double eta,MatDoub &J,VecDoub &z,MatDoub &AJ,int winner,double sigma);
void ModificoPesos(MatDoub &pesos,MatDoub &delta);
	
int main(){

	int N_IN=2,N_OUT=10;
	double eta=0.01;
	double sigma=1;
	
	MatDoub J(N_OUT,N_IN);
	VecDoub z(N_IN);
	MatDoub AJ(N_OUT,N_IN);
	VecDoub h(N_OUT);
	double prob=0;
	
	srand(time(0));
	InitWeigh(J);
	
	while(prob==0){
		GeneraPatronUniforme(z);
		prob=ProbPolar(z);
		//cout << prob << endl;
	}
	FILE *datos;
	int winner;
	int tiempo=0;
	int tmax=10000;
	while(tiempo<tmax){
		CalculoH(J,z,h);
		if(tiempo%100==0){
		string output=to_string(tiempo)+to_string("_")+to_string(sigma)+to_string(".dat");
        if((datos = fopen(output.c_str(), "w")) == NULL){
            printf("No puedo abrir el archivo %s.\n", output.c_str());
            exit(1);
        }
		for(int i=0;i<J.nrows();i++){
			for(int j=0;j<J.ncols();j++)
				fprintf(datos,"%lf\t", J[i][j]);
				fprintf(datos,"\n");
				}
			}
		winner = MaxVectorPosition(h);
		CalculoDeltaJ(eta,J,z,AJ,winner,sigma);
		ModificoPesos(J,AJ);
	tiempo++;
	}
	
	return 0;

}

void GeneraPatronUniforme(VecDoub &z){
	int n=z.size();
    for(int i=0;i<n;i++)
		z[i]=rand()*1.0/RAND_MAX;
	}
	
double ProbPolar(VecDoub &z){
	
	int n=z.size();
	double r=0,theta;
	double p=0;
	
	for(int i=0;i<n;i++)
		r+=z[i]*z[i];
	r=sqrt(r);
	
	theta=atan(z[1]/z[0]);
	
	if(r>0.9 && r<1.1 && theta>0 && theta<atan(1))
		p=1;
	
	}
	
void InitWeigh(MatDoub &matrix){
	
	int n = matrix.nrows();
	int m = matrix.ncols();
	
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			matrix[i][j] = rand()*1.0/RAND_MAX;
			}
		}
	}
	
void CalculoH(MatDoub &J,VecDoub &z,VecDoub &h){
	
	int n=J.nrows();
	//cout << "n=" << n << endl;
	int m=J.ncols();
	//cout << "m=" << m << endl;
	
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			h[i]+=J[i][j]*z[j];
			}
		}
	}
	
int MaxVectorPosition(VecDoub &vec){
	
	int pos;
	int n=vec.size();
	double max;
	
	pos=0;
	max=vec[0];
	
	for(int i=0;i<n;i++){
		if(vec[i]>max){
			pos=i;
			max=vec[i];
			}
		}
	return pos;
	}

double Lambda(int i1,int i2,double sigma){
	
	double resultado;
	double exponente;
	
	exponente=(i1-i2)*(i1-i2);
	exponente=sqrt(exponente);
	
	resultado = exp(-exponente/(2*sigma*sigma));
	
	return resultado;
	}
	
void CalculoDeltaJ(double eta,MatDoub &J,VecDoub &z,MatDoub &AJ,int winner,double sigma){
	
	int n=J.nrows();
	int m=J.ncols();
	
	for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
			AJ[i][j]=eta*Lambda(i,winner,sigma)*(z[j]-J[i][j]);
	}
	
void ModificoPesos(MatDoub &pesos,MatDoub &delta){
	
	int n=pesos.nrows();
	int m=pesos.ncols();
	
	for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
			pesos[i][j]+=delta[i][j];

	}
