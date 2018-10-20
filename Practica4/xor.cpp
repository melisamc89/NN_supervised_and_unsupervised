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
	
	int N_IN,N_HIDE1,N_OUT,P;
	printf("Escriba el numero de neuronas de entrada.\n");
	scanf("%d",&N_IN);
	printf("Escriba el numero de neuronas de la capa oculta.\n");
	scanf("%d",&N_HIDE1);
	printf("Escriba el numero neuronas de salida\n");
	scanf("%d",&N_OUT);
	printf("Escriba el numero neuronas patrones\n");
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
        
	double eta = 0.05;
	
	MatDoub J(N_HIDE1,N_IN),W(N_OUT,N_HIDE1);
	int estado_umbral=-1;
	VecDoub umbral1(N_HIDE1),umbral2(N_OUT);
	MatDoub z(P,N_IN),v(P,N_HIDE1),h1(P,N_HIDE1),salida_deseada(P,N_OUT),o(P,N_OUT),error(P,N_OUT),h2(P,N_OUT);

	srand(time(0));
	
	InitWeigh(J);
	InitWeigh(W);
	InitWeighUmbral(umbral1);
	InitWeighUmbral(umbral2);
	//XOR_patrones(z);
	GeneraPatrones(z);
	ExpectedAnswer_XOR(z,salida_deseada);

	double er=1;
	int it=0;
	while(er > 0.0001){
					
		Calculo(z,J,h1,estado_umbral,umbral1);
		FuntoMat(h1,v,funtanh);
		Calculo(v,W,h2,estado_umbral,umbral2);
		FuntoMat(h2,o,funtanh);

		CalculoError(salida_deseada,o,error);
		er=SumaMat2(error);
		cout << it << " " << er << endl;
		fprintf(archivo1,"%i\t%lf\n",it,er);
		//er=0;
		
		MatDoub delta1(P,N_OUT);
		CalculoDelta(delta1,error,h2,inv_cosh2);

		MatDoub AW(N_OUT,N_HIDE1);
		VecDoub AUmbral2(N_OUT);
		CalculoAWeigh(eta,delta1,v,AW,estado_umbral,AUmbral2);
	
		ModificoPesos(W,AW,umbral2,AUmbral2);
		
		MatDoub delta2(P,N_HIDE1);
		CalculoDelta2(delta2,delta1,W,h1,inv_cosh2);

		MatDoub AJ(N_HIDE1,N_IN);
		VecDoub AUmbral1(N_HIDE1);
		CalculoAWeigh(eta,delta2,z,AJ,estado_umbral,AUmbral1);

		ModificoPesos(J,AJ,umbral1,AUmbral1);
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
		fprintf(archivo2,"%lf\n",umbral2[i]);
		}
		
	return 0;
	}
