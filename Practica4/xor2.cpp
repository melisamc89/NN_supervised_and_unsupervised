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

void InitWeigh(MatDoub &matrix);
void XOR_patrones(MatDoub &matrix);
	
#define N_IN 2
#define N_HIDE1 1
#define N_OUT 1
#define P 4
#define eta 0.01

int main(){
	
	MatDoub patrones(P,N_IN);
	MatDoub J(N_HIDE1,N_IN),W(N_OUT,N_HIDE1+N_IN);
	VecDoub umbral1(N_HIDE1);
	double umbral2;
	int estado_umbral=-1;
	double salida_deseada,o;
	VecDoub v(N_HIDE1);	
	double error,delta;
	VecDoub delta2(N_HIDE1+N_IN);

	srand(time(0));
	
	InitWeigh(J);
	InitWeigh(W);
	for(int i=0;i<N_HIDE1;i++)
		umbral1=rand()*1.0/RAND_MAX;
	umbral2 = rand()*1.0/RAND_MAX;
	XOR_patrones(patrones);

	double error_global=1;
	int it=0;
	while(error_global>0.001){
		error_global=0;
		for(int i=0;i<P;i++){
			
			salida_deseada=patrones[i][0]*patrones[i][1];
			
			for(int k=0;k<N_HIDE1;k++){
				v[k]=0;
				for(int j=0;j<N_IN;j++)
					v[k]+=J[k][j]*patrones[i][j];
				v[k]+=umbral1[k]*estado_umbral;			
				v[k]=tanh(v[k]);
				}
			o=0;	
			for(int k=0;k<N_HIDE1;k++)
				o+=W[0][k]*v[k];
			for(int j=0;j<N_IN;j++)
				o+=W[0][N_HIDE1+j]*patrones[i][j];
			o+=umbral2*estado_umbral;
			o=tanh(o);
			//cout << salida_deseada << " " << o << endl;
			error = 0.5*(salida_deseada-o)*(salida_deseada-o);	
			error_global+=error;
			delta = (salida_deseada-o)*(1-o*o);

			for(int k=0;k<N_HIDE1;k++)
				delta2[k]=delta*W[0][k]*(1-v[k]*v[k]);
				
			for(int k=0;k<N_HIDE1;k++)
				W[0][k]+=eta*delta*v[k];
			for(int j=0;j<N_IN;j++)
				W[0][N_HIDE1+j]+=eta*delta*patrones[i][j];
			umbral2+=eta*delta*estado_umbral;	

			for(int k=0;k<N_HIDE1;k++){
				for(int j=0;j<N_IN;j++)
					J[k][j]+=eta*delta2[k]*patrones[i][j];
				umbral1[k]+=eta*delta2[k]*estado_umbral;
				}
			}
		//cout << it << " " << error_global << endl;
		it++;
		}	
	
	for(int i=0;i<N_HIDE1;i++){
		for(int j=0;j<N_IN;j++)
			cout << J[i][j] << endl;
		cout << umbral1[i] << endl;
		cout << endl;
		}
	
	for(int i=0;i<N_OUT;i++){
		for(int j=0;j<N_HIDE1;j++)
			cout << W[i][j] << endl;
		cout << umbral2 << endl;
		cout << endl;
		}
	
		
	return 0;
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

void XOR_patrones(MatDoub &matrix){
	
	int n = matrix.nrows();
	int m = matrix.ncols();
	
	matrix[0][0]=1;
	matrix[0][1]=1;
	matrix[1][0]=-1;
	matrix[1][1]=-1;
	matrix[2][0]=1;
	matrix[2][1]=-1;
	matrix[3][0]=-1;
	matrix[3][1]=1;
	
	}
