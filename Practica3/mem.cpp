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

int main(){

	int N,P;
	printf("Escriba el numero de neuronas de la red.\n");
	scanf("%d",&N);
	printf("Escriba el numero de patrones que desea almacenar\n");
	scanf("%d",&P);
		
	FILE *over;
	string red_overlap = to_string("Hopfield_overlap_") + to_string(N)+ to_string("_")+to_string(P)+to_string(".dat");
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
	
	for(int p=0;p<=P;p++){
		for(int i=0;i<N;i++){
			estados[i]=patrones[i][p];
			estado_inicial[i]=patrones[i][p];
			//cout << estados[i] << "\t" << estado_inicial[i] << endl;
			}
		/*FILE *archivo;
		if((archivo=fopen("Evolucion.dat","w"))==NULL){
			printf("Problemas para abrir archivo evolucion.dat");
			}*/
		int conv=1;
		int convtime=0;
		while(conv != 0 && convtime<500){
			convtime++;
			conv=0;
			for(int i=0;i<N;i++)
				aux[i]=estados[i];
			Dinamica(estados,conexiones);
			for(int i=0;i<N;i++){
				if(estados[i]*aux[i]<0)
					conv=1;
				}
			/*for(int i=0;i<N;i++)
				fprintf(archivo,"%lf\t",estados[i]);
			fprintf(archivo,"\n");*/
		}
		
		overlap[p] = CalculaOverlap(estados,estado_inicial);
		cout << p << " " << overlap[p] << " " << convtime << endl;
		fprintf(over,"%lf\n", overlap[p]);
		
	}
	
	fclose(over);
	
	return 0;
	}
