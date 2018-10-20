#ifndef FUNCIONESINT_H
#define FUNCIONESINT_H

typedef NRVec<int> VecInt, VecInt_O, VecInt_IO;
typedef const NRVec<int>VecInt_I;

typedef NRVec<double> VecDoub, VecDoub_O, VecDoub_IO;
typedef const NRVec<double> VecDoub_I;

typedef NRMat<int> MatInt, MatInt_O, MatInt_IO;
typedef const NRMat<int> MatInt_I;

typedef NRMat<double> MatDoub, MatDoub_O, MatDoub_IO;
typedef const NRMat<double> MatDoub_I;

//------------------DECLARACION DE FUNCIONES----------------------------
void GeneraPatrones(MatDoub &matrix);
void InitWeigh(MatDoub &matrix);
void InitWeighUmbral(VecDoub &vec);
void XOR_patrones(MatDoub &matrix);
int ExpectedAnswer_XOR(MatDoub &entrada,MatDoub &salida);
void Calculo(MatDoub &inicial,MatDoub &pesos,MatDoub &final,int estado_umbral,VecDoub &umbral);
void Calculob(MatDoub &patrones,MatDoub &inicial,MatDoub &pesos,MatDoub &final,int estado_umbral,VecDoub &umbral);
void FuntoMat(const MatDoub &entrada,MatDoub &salida,double (*fun)(double));
int Signo(double x);
void CalculoError(const MatDoub &esperado,const MatDoub &salida,MatDoub &error);
double inv_cosh2(double x);
void CalculoDelta(MatDoub &delta,const MatDoub &error,const MatDoub &h,double (*fun)(double));
void CalculoDelta2(MatDoub &delta2,MatDoub &delta1,MatDoub &pesos,const MatDoub &h,double (*fun)(double));
void CalculoAWeigh(double eta,MatDoub &delta,MatDoub &v,MatDoub &A,int estado_umbral,VecDoub &AUmbral);
void CalculoAWeighb(double eta,MatDoub &delta,MatDoub &z,MatDoub &v,MatDoub &A,int estado_umbral,VecDoub &AUmbral);
void ModificoPesos(MatDoub &pesos,MatDoub &delta,VecDoub &umbral,VecDoub &deltaumbral);
int Comparacion(const MatDoub &salida,const MatDoub salida_deseada);
double SumaMat2(const MatDoub &mat);
	
//------------------------FUNCIONES-----------------------------	
		
void GeneraPatrones(MatDoub &matrix){
	
	int p = matrix.ncols();
	int n = matrix.nrows();
	
	double aleatorio;
	for(int i=0;i<n;i++){
		for(int j=0;j<p;j++){
			aleatorio = rand()*1./RAND_MAX;
			if(aleatorio > 0.5)
				matrix[i][j] = 1;
			else
				matrix[i][j] = -1;
			}
		}
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
	
void InitWeighUmbral(VecDoub &vec){
	
	int n = vec.size();
	
	for(int i=0;i<n;i++){
		vec[i] = rand()*1.0/RAND_MAX;
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

void LOGISTICO_patrones(MatDoub &matrix){
	
	int n = matrix.nrows();
	int m = matrix.ncols();
	
	for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
			matrix[i][j]= rand()*1.0/RAND_MAX;
	
	}

int ExpectedAnswer_XOR(MatDoub &entrada,MatDoub &salida){
	
	int n = entrada.nrows();
	int m = entrada.ncols();
	
	int l=salida.nrows();

	if(n!=l)
		cout << "problemas de dimension" << endl;
	
	VecInt prod(n);
	for(int i=0;i<n;i++)
		prod[i]=1;
	
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			prod[i]=prod[i]*entrada[i][j];
			}
		}
		
	for(int i=0;i<n;i++)
		salida[i][0]=prod[i];
	}
	
int ExpectedAnswer_LOGISTICO(MatDoub &entrada,MatDoub &salida){
	
	int n = entrada.nrows();
	int m = entrada.ncols();
	
	int l=salida.nrows();

	if(n!=l)
		cout << "problemas de dimension" << endl;
	
	VecInt out(n);
	
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			out[i]=4*entrada[i][j]*(1-entrada[i][j]);
			}
		}
		
	for(int i=0;i<n;i++)
		salida[i][0]=out[i];
	}
	
int Signo(double x){
	if(x>0)
		return 1;
	else
		return -1;
	}
	
void Calculo(MatDoub &inicial,MatDoub &pesos,MatDoub &final,int estado_umbral,VecDoub &umbral){
	
	int n=inicial.nrows();
	int m=inicial.ncols();
	
	int l=final.nrows();
	int s=final.ncols();
	
	double suma;
	
	for(int i=0;i<l;i++){
		for(int j=0;j<s;j++){
			suma=0;
			for(int k=0;k<m;k++)
				suma+=pesos[j][k]*inicial[i][k];	
			final[i][j]=suma+estado_umbral*umbral[j];
			}
		}	
	}
	
void Calculob(MatDoub &patrones,MatDoub &inicial,MatDoub &pesos,MatDoub &final,int estado_umbral,VecDoub &umbral){
	
	int n=inicial.nrows();
	int m=inicial.ncols();
	
	int l=final.nrows();
	int s=final.ncols();
	double suma;
	
	for(int i=0;i<l;i++){
		for(int j=0;j<s;j++){
			suma=0;	
			for(int k=0;k<patrones.ncols();k++)
				suma+=pesos[j][k]*patrones[i][k];
			for(int k=0;k<m;k++)
				suma+=inicial[i][k]*pesos[j][k];
			final[i][j]+=suma+estado_umbral*umbral[j];
			}
		}
		
	}
	
void FuntoMat(const MatDoub &entrada,MatDoub &salida,double (*fun)(double)){
	
	int n=entrada.nrows();
	int m=entrada.ncols();
	
	for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
			salida[i][j]=fun(entrada[i][j]);
	
	}	
	
void CalculoError(const MatDoub &esperado,const MatDoub &salida,MatDoub &error){
	
	int m=salida.ncols();
	int n=salida.nrows();

	for(int i=0;i<n;i++)
		for (int j=0;j<m;j++)
			error[i][j]=(esperado[i][j]-salida[i][j]);
		
	}
	
double SumaMat2(const MatDoub &mat){
	
	int n=mat.nrows();
	int m=mat.ncols();
	double suma=0;
	
	for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
			suma+=mat[i][j]*mat[i][j];
			
	return 0.5*suma;
	}
	
double inv_cosh2(double x){
	double s;
	s = 1.1/(cosh(x)*cosh(x));
	return s;
	}
	
double funtanh(double x){
	return 1.1*tanh(x);
	}
	
double fun_lineal(double x){
	return x;
	}

double fun_uno(double x){
	return 1;
	}

void CalculoDelta(MatDoub &delta,const MatDoub &error,const MatDoub &h,double (*fun)(double)){
	
	int n= delta.nrows();
	int m= delta.ncols();
	
	MatDoub mataux(h.nrows(),h.ncols());
	FuntoMat(h,mataux,fun);
	
	for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
			delta[i][j]=mataux[i][j]*error[i][j];
			
	}
	
void CalculoAWeigh(double eta,MatDoub &delta,MatDoub &v,MatDoub &A,int estado_umbral,VecDoub &AUmbral){
	
	int n=A.nrows();
	int m=A.ncols();
	int l=delta.nrows();
	int s=AUmbral.size();
	double suma;
	
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			suma=0;
			for(int k=0;k<l;k++)
				suma+=delta[k][i]*v[k][j];
			A[i][j]=suma*eta;
			}
		}
		
	for(int i=0;i<s;i++){
		suma=0;
		for(int k=0;k<l;k++)
			suma+=delta[k][i]*estado_umbral;
		AUmbral[i]=suma*eta;
		}
	}
	
void CalculoAWeighb(double eta,MatDoub &delta,MatDoub &z,MatDoub &v,MatDoub &A,int estado_umbral,VecDoub &AUmbral){
	
	int n=A.nrows();
	int m=A.ncols();
	int l=delta.nrows();
	int s=AUmbral.size();
	double suma;
	
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			suma=0;
			for(int k=0;k<l;k++){
				suma+=delta[k][i]*(v[k][j]+z[k][i]);
				}
			A[i][j]=suma*eta;
			}
		}
		
	for(int i=0;i<s;i++){
		suma=0;
		for(int k=0;k<l;k++)
			suma+=delta[k][i]*estado_umbral;
		AUmbral[i]=suma*eta;
		}
	}

void CalculoDelta2(MatDoub &delta2,MatDoub &delta1,MatDoub &pesos,const MatDoub &h,double (*fun)(double)){
	
	int l = delta1.ncols();
	int n=delta2.nrows();
	int m=delta2.ncols();
	double suma;
	
	MatDoub mataux(h.nrows(),h.ncols());
	FuntoMat(h,mataux,fun);
	
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			suma=0;
			for(int k=0;k<l;k++)
				suma+=delta1[i][k]*pesos[k][j];
			delta2[i][j]=mataux[i][j]*suma;
			}
		}
	}
	
void ModificoPesos(MatDoub &pesos,MatDoub &delta,VecDoub &umbral,VecDoub &deltaumbral){
	
	int n=pesos.nrows();
	int m=pesos.ncols();
	int l=umbral.size();
	
	for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
			pesos[i][j]+=delta[i][j];
	
	for(int i=0;i<l;i++)
		umbral[i]+=deltaumbral[i];
	
	}
	
int Comparacion(const MatDoub &salida,const MatDoub salida_deseada){
	
	int n=salida.nrows();
	int m=salida.ncols();
	int comp=0;
	
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			if(salida[i][j]!=salida_deseada[i][j]);
				comp++;
		}
	}
			
	return comp;
	}
	
#endif
