#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <string>
#include <fstream>
#include <iostream>
#include "nrutilc.h"
#include "LUdcmp.h"
#include "jacobi.h"

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
	
		int n,elemento;
		
		printf("Escriba el numero de filas.\n");
		scanf("%d",&n);
		MatDoub matrix(n,n);
		
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				printf("%d %d ", i,j);
				scanf("%d",&elemento);
				matrix[i][j]=elemento;
			}
		}

		double det;
		//------------------------------------------CALCULO DETERMINANTE
		
		//cout << "Descomposicion LU" << endl;
		LUdcmp auxiliar=LUdcmp(matrix);
		det = auxiliar.det();
		cout << "determinante  = " << det << endl;
		//-------------------------------------------CALCULO AUTOVALORES
		
		//cout << "Jacobi" << endl;
		//Jacobi jac(diagonalizar);
		Jacobi jac(matrix);
		eigsrt(jac.d,&jac.v);
		
		//--------------------------ESCRITURA AUTOVALORES Y AUTOVECTORES
		
		cout << "Escritura de autovalores" << endl;
		for(int i=0;i<n;i++)
			cout << jac.d[i] << endl;
			
		cout << "Escritura de autovectores" << endl;	
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				cout << jac.v[i][j] << " ";
			}
			cout << endl;
				}
	
		
return 0;
}
