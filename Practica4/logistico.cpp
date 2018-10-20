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

#include "funciones.h"

template<class T>
inline string to_string(const T& t){
    stringstream ss;
    ss << t;
    return ss.str();
}

int main(){
	
	int N_IN=1,N_HIDE1=5,N_OUT=1,P=100;
	/*printf("Escriba el numero de neuronas de entrada.\n");
	scanf("%d",&N_IN);
	printf("Escriba el numero de neuronas de la capa oculta.\n");
	scanf("%d",&N_HIDE1);
	printf("Escriba el numero neuronas de salida\n");
	scanf("%d",&N_OUT);
	printf("Escriba el numero de patrones\n");
	scanf("%d",&P);*/
	
	double eta = 0.001;
	
	VecDoub J(N_HIDE1),W(N_HIDE1);
	double umbral,deltaumbral;
	VecDoub salida_deseada(P);
	VecDoub z(P),o(P);
	VecDoub v(N_HIDE1),delta2(N_HIDE1);
	double delta1;
	
	srand(time(0));
	
	for(int i=0;i<N_HIDE1;i++){
		J[i]=rand()*1.0/RAND_MAX;
		W[i]=rand()*1.0/RAND_MAX;		
		}
	umbral = rand()*1.0/RAND_MAX;
	for(int k=0;k<P;k++)
		z[k]=k*1.0/P;

	double error_patrones=1;
	int it=0;
	while(error_patrones > 0.05){
		error_patrones=0;
		for(int k=0;k<P;k++){	
			salida_deseada[k] = 4*z[k]*(1-z[k]);
			for(int i=0;i<N_HIDE1;i++)
				v[i]=tanh(J[i]*z[k]);
		
			o[k]=umbral*z[k];	
			for(int i=0;i<N_HIDE1;i++)
				o[k]+=W[i]*v[i];			
			
			delta1 = (salida_deseada[k]-o[k]);
			for(int i=0;i<N_HIDE1;i++)
				delta2[i]=W[i]*delta1*(1-v[i]*v[i]);
				
			for(int i=0;i<N_HIDE1;i++){
				W[i]+=eta*delta1*v[i];
				J[i]+=eta*delta2[i]*z[k];
				}
			umbral+=eta*delta1*z[k];
			error_patrones+=delta1*delta1*0.5;
		}
		//cout << it << " " << error_patrones << endl;
		it++;
		}
	
	double error=0;		
	for(int i=0;i<P;i++)
		error+=0.5*(salida_deseada[i]-o[i])*(salida_deseada[i]-o[i]);
	error/=P;
		cout << "Error de aprendizaje = " << error << endl;
		//cout << z[i] << " " << salida_deseada[i] << " " << o[i] << endl;
		
	double s,salida;
	it=0;
	error=0;
	double out;
	while(it<1000){	
		s=rand()*1.0/RAND_MAX;
		out=4*s*(1-s);
		VecDoub v(N_HIDE1);
		for(int i=0;i<N_HIDE1;i++)
			v[i]=tanh(s*J[i]);	
		salida=umbral*s;
		for(int i=0;i<N_HIDE1;i++)	
			salida+=v[i]*W[i];
		//cout << s << " " << salida << endl;;
		error+=0.5*(out-salida)*(out-salida);
		it++;	
		}
	error/=1000;
	cout << "Error de generalizacion = " << error << endl;
		
	return 0;
	}
