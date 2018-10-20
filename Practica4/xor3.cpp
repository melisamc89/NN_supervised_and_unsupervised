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
void GeneraPatrones(MatDoub &matrix);
	
int main(){

	int N_IN,N_HIDE1,N_OUT,P;
	printf("Escriba el numero de neuronas de entrada.\n");
	scanf("%d",&N_IN);
	printf("Escriba el numero de neuronas de la capa oculta.\n");
	scanf("%d",&N_HIDE1);
	printf("Escriba el numero neuronas de salida\n");
	scanf("%d",&N_OUT);
	printf("Escriba el numero de patrones\n");
	scanf("%d",&P);
	
	FILE *archivo1,*archivo2;
	string output1=to_string("error_xor")+to_string("_")+to_string(N_IN)+to_string("_")+to_string(N_HIDE1)+to_string(".dat");
	if((archivo1 = fopen(output1.c_str(), "w")) == NULL){
            printf("No puedo abrir el archivo %s.\n", output1.c_str());
            exit(1);
        }
        
	string output2=to_string("salidas_xor")+to_string("_")+to_string(N_IN)+to_string("_")+to_string(N_HIDE1)+to_string(".dat");
	if((archivo2 = fopen(output2.c_str(), "w")) == NULL){
            printf("No puedo abrir el archivo %s.\n", output1.c_str());
            exit(1);
        }
     
    double eta = 0.01;   
	MatDoub patrones(P,N_IN);
	MatDoub J(N_HIDE1,N_IN),W(N_OUT,N_HIDE1);
	VecDoub umbral1(N_HIDE1);
	double umbral2;
	int estado_umbral=-1;
	double o;
	VecDoub v(N_HIDE1);	
	double error,delta;
	VecDoub delta2(N_HIDE1);

	srand(time(0));
	
	InitWeigh(J);
	InitWeigh(W);
	for(int i=0;i<N_HIDE1;i++)
		umbral1=rand()*1.0/RAND_MAX;
	umbral2 = rand()*1.0/RAND_MAX;
	GeneraPatrones(patrones);
	VecDoub salida_deseada(P);
	for(int i=0;i<P;i++){
		double prod=1;
		for(int j=0;j<N_IN;j++){
			prod=prod*patrones[i][j];
			}
		salida_deseada[i]=prod;
		}

	double error_global=1;
	int it=0;
	while(error_global>0.001){
		error_global=0;
		for(int i=0;i<P;i++){
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
			o+=umbral2*estado_umbral;
			o=tanh(o);
			cout << salida_deseada[i] << " " << o << endl;
			error = 0.5*(salida_deseada[i]-o)*(salida_deseada[i]-o);	
			error_global+=error;
			delta = (salida_deseada[i]-o)*(1-o*o);

			for(int k=0;k<N_HIDE1;k++)
				delta2[k]=delta*W[0][k]*(1-v[k]*v[k]);

			for(int k=0;k<N_HIDE1;k++)
				W[0][k]+=eta*delta*v[k];	
			umbral2+=eta*delta*estado_umbral;	

			for(int k=0;k<N_HIDE1;k++){
				for(int j=0;j<N_IN;j++)
					J[k][j]+=eta*delta2[k]*patrones[i][j];
				umbral1[k]+=eta*delta2[k]*estado_umbral;
				}
			}
		fprintf(archivo1,"%i\t%lf\n",it,error_global);
		it++;
		}
		
	for(int i=0;i<N_HIDE1;i++){
		for(int j=0;j<N_IN;j++)
			fprintf(archivo2,"%lf\n",J[i][j]);
		fprintf(archivo2,"%lf\n",umbral1[i]);
		}
	
	for(int i=0;i<N_OUT;i++){
		for(int j=0;j<N_HIDE1;j++)
			fprintf(archivo2,"%lf\n",W[i][j]);
		fprintf(archivo2,"%lf\n",umbral2);
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

		
void GeneraPatrones(MatDoub &matrix){
	
	int n = matrix.nrows();
	int m = matrix.ncols();
	
	double aleatorio;
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			matrix[i][j]=1;
			aleatorio = rand()*1./RAND_MAX;
			if(aleatorio < 0.5)
				matrix[i][j] = -1;
			}
		}
	}
