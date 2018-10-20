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

void GeneraPatrones(MatInt &matrix);
void MatConexiones(MatInt &matrix,MatDoub &Conexiones);
int Signo(double x);
void Dinamica(VecDoub &estados,MatDoub &conexiones);
double CalculaOverlap(VecDoub &efinal,VecDoub &einicial);

void GeneraPatrones(MatInt &matrix){
	
	int p=matrix.ncols();
	int n=matrix.nrows();
	
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

void MatConexiones(MatInt &matrix,MatDoub &Conexiones){
	
	int n = matrix.nrows();
	int m = matrix.ncols();
	
	for(int i=0;i<n;i++){
		for(int j=0;j<i;j++){
			for(int k=0;k<m;k++)
				Conexiones[i][j]+=matrix[i][k]*matrix[j][k];
			Conexiones[i][j]/=n;
			}
		}
		
	for(int i=0;i<n;i++){
		for(int j=i;j<n;j++){
			Conexiones[i][j]=Conexiones[j][i];
			}
		}
	
	for(int i=0;i<n;i++){
		Conexiones[i][i]=0;
		}
	}
	
int Signo(double x){
	if(x>0)
		return 1;
	else
		return -1;
	}
	
void Dinamica(VecDoub &estados,MatDoub &conexiones){
	
	int n = estados.size();
	
	VecDoub auxiliar(n);
	
	for(int i=0;i<n;i++){
		auxiliar[i]=estados[i];
		estados[i]=0;	
		}
	
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++)
			estados[i]+=conexiones[i][j]*auxiliar[j];
		estados[i]=Signo(estados[i]);
		}	
	}
	
double CalculaOverlap(VecDoub &efinal,VecDoub &einicial){
	
	double mu,suma=0;
	int n=efinal.size();
	
	for(int i=0;i<n;i++)
		suma+=efinal[i]*einicial[i];
	mu = suma*1.0/n;
	return mu;
	}
	
void CalculaH(const MatDoub &conexiones,const VecDoub &estados,VecDoub &h){
	
	int n=h.size();
	
	for(int i=0;i<n;i++)
		h[i]=0;
	
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			h[i]+=conexiones[i][j]*estados[j];
		}
		//cout << h[i] << endl;
	}

	}

void CalculaProb(const VecDoub &h,MatDoub &prob,double beta){
	
	int n=h.size();
	
	for(int i=0;i<n;i++){
		prob[i][0]=exp(beta*h[i])/(exp(beta*h[i])+exp(-beta*h[i]));
		prob[i][1]=exp(-beta*h[i])/(exp(beta*h[i])+exp(-beta*h[i]));		
		}

	}
	
void DinamicaProb(const MatDoub &prob,VecDoub &estados){
	
	int n=estados.size();
	
	for(int i=0;i<n;i++){
		estados[i]=1*prob[i][0]+(-1)*prob[i][1];
		}
	}

double PromVector(const VecDoub &vector){
	
	int n=vector.size();
	double prom=0;
	
	for(int i=0;i<n;i++)
		prom+=vector[i];
	prom/=n;
	return prom;
	}


#endif
