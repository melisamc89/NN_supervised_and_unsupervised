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

void GeneraCoordenadas(MatDoub &mat);
void InitWeigh(MatDoub &matrix);
int MaxVectorPosition(VecDoub &vec);
void OrdenaVector(VecDoub &vector,MatDoub &mat);
double CalculoDistancia(const MatDoub &mat);
	
int main(){
	
	int N_IN=2,N_OUT;
	double alpha=0.8;
	
	for(N_OUT=50;N_OUT<=100;N_OUT++){
	
		double d = N_OUT;
		MatDoub J(N_OUT,N_IN);
		VecDoub z(N_IN);
		VecDoub angulo(N_OUT);
		VecDoub ord_ciudades(N_OUT);
		
		MatDoub ciudades(N_OUT,N_IN);
		GeneraCoordenadas(ciudades);
		
		for(int i=0;i<N_OUT;i++)
			angulo[i]=i*4*1.0*asin(1)/N_OUT;
		
		/*FILE *coordenadas;
		coordenadas=fopen("coordenadas.dat","w");
		for(int i=0;i<N_OUT;i++)
			fprintf(coordenadas,"%i\t%lf\t%lf\n",i,ciudades[i][0],ciudades[i][1]);*/
		
		srand(time(0));
		InitWeigh(J);
		
		int winner;
		int tiempo=0;	
		int tmax=1;	
		double d0=d;
		while(tiempo<=tmax){
			for(int i=0;i<ciudades.nrows();i++){
				VecDoub h(N_OUT);
				for(int k=0;k<N_OUT;k++){
					for(int j=0;j<N_IN;j++)
						h[k]+=J[k][j]*ciudades[i][j];
						/*if(tiempo==tmax && i==2)
							cout << h[k] << endl;*/
						}			
				winner = MaxVectorPosition(h);
				for(int k=0;k<N_OUT;k++){
					double dist = ((k-winner)%N_OUT);
					if(dist<=d){
						for(int j=0;j<N_IN;j++)
							J[k][j]= J[k][j] + alpha*(ciudades[i][j]-J[k][j]);
						}
					}
			if(tiempo==tmax){
				double prom_cos=0,prom_sin=0,prom=0;
				for(int j=0;j<N_OUT;j++){
					prom_cos+=h[j]*cos(angulo[j]);
					prom_sin+=h[j]*sin(angulo[j]);
					}	
				prom=atan2(prom_sin,prom_cos);
				if(prom<0)
					prom+=4*asin(1);
				ord_ciudades[i]=prom;
				}	
					
			}
			alpha-=0.01;
			d=d-(d0-1)/tmax;
			tiempo++;
		}
		
		OrdenaVector(ord_ciudades,ciudades);
		double distancia = CalculoDistancia(ciudades);
		cout << N_OUT << " " << distancia << endl;
	}

	return 0;
}

void GeneraCoordenadas(MatDoub &mat){
	
	int n=mat.nrows();
	int m=mat.ncols();
	
    for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
			mat[i][j]=rand()*1.0/RAND_MAX;
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
	
void OrdenaVector(VecDoub &vector,MatDoub &mat){
	
	double a,x,y;
	int n=vector.size();
	
	for(int j=1;j<n;j++)
		for(int i=0;i<n-j;i++)
			if(vector[i]>vector[i+1]){
				a=vector[i+1];
				vector[i+1]=vector[i];
				vector[i]=a;
				x=mat[i+1][0];
				mat[i+1][0]=mat[i][0];
				mat[i][0]=x;
				y=mat[i+1][1];
				mat[i+1][1]=mat[i][1];
				mat[i][1]=y;
			}
	}
	
double CalculoDistancia(const MatDoub &mat){
	
	int n=mat.nrows();
	double dist=0;
	for(int i=0;i<n-1;i++)
		dist+=sqrt((mat[i][0]-mat[i+1][0])*(mat[i][0]-mat[i+1][0])+(mat[i][1]-mat[i+1][1])*(mat[i][1]-mat[i+1][1]));

	return dist;
	}
