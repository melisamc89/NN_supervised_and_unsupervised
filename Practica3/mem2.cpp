#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <fstream>
#include <iostream>
#include <iostream>
#include <iomanip>
#include "nrutilc.h"
#include "funciones.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <malloc.h>

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

#define ITMAX 10

int main(){

	int N,P;
	printf("Escriba el numero de neuronas de la red.\n");
	scanf("%d",&N);
	printf("Escriba el numero de patrones que desea almacenar\n");
	scanf("%d",&P);
		
	FILE *over;
	string red_overlap = to_string("Hopfield2_overlapmedio_") + to_string(N)+ to_string("_")+to_string(P)+to_string(".dat");
	if((over = fopen(red_overlap.c_str(), "w"))==NULL){
			printf("Problemas para abrir archivo overlap.dat");
			}
	
	MatInt patrones(N,P);
	srand(time(0));
	GeneraPatrones(patrones);
	MatDoub conexiones(N,N);
	MatConexiones(patrones,conexiones);

	VecDoub estados(N),estado_inicial(N);
	VecDoub overlap(P);
	VecDoub aux(N);
	VecDoub h(N);
	MatDoub prob(N,2);
	double beta;
	double TMAX=2;
	double temperature=0.1;
	double overlap_medio;
	
	while(temperature<=TMAX){
		for(int p=0;p<=P;p++){
			for(int i=0;i<N;i++){
				estados[i]=patrones[i][p];
				estado_inicial[i]=patrones[i][p];
				}
			int iteraciones=0;
			while(iteraciones < ITMAX){
				CalculaH(conexiones,estados,h);
				beta=1./temperature;
				CalculaProb(h,prob,beta);
				DinamicaProb(prob,estados);
				iteraciones++;
			}
			overlap[p] = CalculaOverlap(estados,estado_inicial);
			}
		overlap_medio=PromVector(overlap);	
		fprintf(over,"%lf\t%lf\n", temperature ,overlap_medio);
		cout << temperature << " " << overlap_medio << endl;
		temperature+=0.1;
		}
	
	fclose(over);
	
	return 0;
	}
