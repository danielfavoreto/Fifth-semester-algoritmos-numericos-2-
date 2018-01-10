#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "GMRES.h"

void subtract(double* uc,double beta,double* u,int tamU) 
{

	int i;
	
	for (i = 0; i < tamU; i++) 
	{
	
		uc[i] -= beta*u[i];
		
	}
	
}

void residuo(MAT* A,double* x,double* Fp,double* u,int tamX) 
{

	int i;
	
	double* r = CSR_MATRIX_ARRAY_Product(A,x,tamX);
	
	for (i = 0; i < tamX; i++) 
	{

		u[i] = Fp[i] - r[i];
		
	}
	
	free(r);
	
}

double norm(double* v,int tamV) 
{

	double sum = 0;
	
	int i;
	
	for (i = 0; i < tamV; i++) 
	{
	
		sum += pow(v[i],2);
		
	}
	
	return sqrt(sum);
	
}

double prod_interno(double* a,double* b,int k) 
{

	int i;
	
	double r = 0;
	
	for (i = 0; i < k; i++) 
	{
	
		r += a[i]*b[i];
		
	}
	
	return r;
	
}

double* GMRES(MAT* Ap,MAT* L,MAT* U,double* x,double* Fp,double tol,int nVetBase,int maxIter,int tamSol) 
{

	double epsilon = tol*norm(Fp,tamSol);
	
	double** u,*e,rol,**H,*c,*s,r,*y;
	
	double *p,*q,*w;
	
	int iter = 1,I,j;

	e = (double*) calloc((nVetBase+1),sizeof(double));
	u = (double**) calloc((nVetBase+1),sizeof(double*));
	H = (double**) calloc((nVetBase+1),sizeof(double*));
	c = (double*) calloc(nVetBase,sizeof(double));
	s = (double*) calloc(nVetBase,sizeof(double));
	y = (double*) calloc(nVetBase,sizeof(double));

	for (I = 0; I < nVetBase + 1; I++) {
	
		u[I] = (double*) calloc(tamSol,sizeof(double));
		
		H[I] = (double*) calloc(nVetBase,sizeof(double));
		
	}

	do {
	
		int i = 0;

		//Ação do precondicionador
		p = (double*)calloc(tamSol,sizeof(double));
		
		residuo(Ap,x,Fp,p,tamSol);
		
		q = solve_system_L(L,p,tamSol);
		
		w = solve_system_U(U,q,tamSol);
		
		for (I = 0; I < tamSol; I++) 
		{

			u[i][I] = w[I];
			
		}
		
		free(p); free(q); free(w);

		e[i] = norm(u[i],tamSol);

		for (I = 0; I < tamSol; I++) 
		{
			
			u[i][I] = u[i][I]/e[i];
			
		}
		
		rol = e[i];

		while (rol > epsilon && i < nVetBase) 
		{
		
			//Ações do precondicionador
			
			p = CSR_MATRIX_ARRAY_Product(Ap,u[i],tamSol);
			q = solve_system_L(L,p,tamSol);
			w = solve_system_U(U,q,tamSol);
			
			for (I = 0; I < tamSol; I++)
			{
			
				u[i+1][I] = w[I];
				
			}
			
			free(p); free(q); free(w);

			//Gram-Schmidt
			for (j = 0; j <= i; j++)
			{
			
				H[j][i] = prod_interno(u[i+1],u[j],tamSol);
				
				subtract(u[i+1],H[j][i],u[j],tamSol);
				
			}

			H[i+1][i] = norm(u[i+1],tamSol);

			for (I = 0; I < tamSol; I++) 
			{
			
				u[i+1][I] = u[i+1][I]/H[i+1][i];
				
			}

			//Algoritmo QR
			for (j = 0; j <= i-1; j++)
			{
			
				double temp = H[j][i];
				
				H[j][i] = c[j]*H[j][i] + s[j]*H[j+1][i];
				
				H[j+1][i] = -s[j]*temp + c[j]*H[j+1][i];
				
			}

			r = sqrt(pow(H[i][i],2) + pow(H[i+1][i],2));
			
			c[i] = H[i][i]/r;
			
			s[i] = H[i+1][i]/r;
			
			H[i][i] = r;
			
			H[i+1][i] = 0;
			
			e[i+1] = -s[i]*e[i];
			
			e[i] = c[i]*e[i];
			
			rol = fabs(e[i+1]);
			
			i++;
		}
		
		i--;

		for (j = i; j >= 0; j--) 
		{
		
			int l;
			
			double sum = 0;
			
			for (l = j+1; l <= i; l++) 
			{
			
				sum += H[j][l]*y[l];
				
			}
			
			y[j] = (e[j]-sum)/H[j][j];
			
		}

		for (j = 0; j <= i; j++) 
		{

			for (I = 0; I < tamSol; I++) 
			{
				
				x[I] += y[j]*u[j][I];
				
			}
		}

		iter++;
		
	} while ((rol >= epsilon) && (iter < maxIter));

	free(e); 
	free(c);
	free(s); 
	free(y); 

	for (I = 0; I < nVetBase + 1; I++) {
	
		free(u[I]);
		
		free(H[I]);
		
	}
	
	free(u); free(H);

	printf("Número de iterações: %d\n",iter);

	return x;
}

double* CSR_MATRIX_ARRAY_Product(MAT* mat,double* arr,int sizeArray){

	double* result = (double*)calloc(sizeArray,sizeof(double));
	
	int i,j;
	
	for (i = 0; i < sizeArray; i++) {

		for (j = mat->IA[i]; j < mat->IA[i+1]; j++) {

			result[i] += mat->AA[j]*arr[mat->JA[j]];
			
		}
	}
	
	return result;
}

double* solve_system_L(MAT* L,double* arr,int size){

	int i,j;
	
	double* r = (double*)calloc(size,sizeof(double));
	
	double s = 0;
	
	r[0] = arr[0]/L->D[0];
	
	for (i = 1; i < size; i++) {

		s = 0;
		
		for (j = L->IA[i]; j < L->IA[i+1];j++) {
		
			s += r[L->JA[j]]*L->AA[j];
			
		}
		
		r[i] = (arr[i] - s)/L->D[i];
		
	}
	return r;
}

double* solve_system_U(MAT* U,double* arr,int size)
{

	int i,j;
	
	double* r = (double*)calloc(size,sizeof(double));
	
	double s = 0;
	
	r[size -1] = arr[size -1]/U->D[size -1];
	
	for (i = size - 2; i >=0; i--)
	{
	
		s = 0;
		
		for (j = U->IA[i]; j < U->IA[i+1];j++)
		{
		
			s += r[U->JA[j]]*U->AA[j];
			
		}
		
		r[i] = (arr[i] - s)/U->D[i];
		
	}
	
	return r;
}


double* GMRES_No_Preconditioner(MAT* Ap,double* x,double* Fp,double tol,int nVetBase,int maxIter)
{

	int tamSol = Ap->n;
	
	double epsilon = tol*norm(Fp,tamSol);
	
	double** u,*e,rol,**H,*c,*s,r,*y;
	
	double *p;
	
	int iter = 1,I,j;
	
	e = (double*) calloc((nVetBase+1),sizeof(double));
	u = (double**) calloc((nVetBase+1),sizeof(double*));
	H = (double**) calloc((nVetBase+1),sizeof(double*));
	c = (double*) calloc(nVetBase,sizeof(double));
	s = (double*) calloc(nVetBase,sizeof(double));
	y = (double*) calloc(nVetBase,sizeof(double));

	for (I = 0; I < nVetBase + 1; I++)
	{
	
		u[I] = (double*)calloc(tamSol,sizeof(double));
		
		H[I] = (double*)calloc(nVetBase,sizeof(double));
		
	}

	do {
	
		int i = 0;

		/*Ação Precondicionador:
		*  u[i] = b - A*x
		*  u[i] = L*U*b - L*U*A*x
		*  u[i] = L*b1 - L*U*p
		*  u[i] = b2 - L*q
		*/
		
		p = CSR_MATRIX_ARRAY_Product(Ap,x,tamSol);

		for (I = 0; I < tamSol; I++) {
		
			u[i][I] = Fp[I] - p[I];
			
		}
		
		free(p);
		
		e[i] = norm(u[i],tamSol);

		for (I = 0; I < tamSol; I++)
		{
		
			u[i][I] = u[i][I]/e[i];
			
		}
		
		rol = e[i];

		while (rol > epsilon && i < nVetBase)
		{
			/*Ação do precondicionador:
			*  u[i+1] = A*u[i]
			*  u[i+1] = L*U*A*u[i]
			*  u[i+1] = L*U*p
			*  u[i+1] = L*q
			*/
			p = CSR_MATRIX_ARRAY_Product(Ap,u[i],tamSol);
			
			for (I = 0; I < tamSol; I++)
			{
			
				u[i+1][I] = p[I];
				
			}
			
			free(p);

			//Gram-Schmidt
			for (j = 0; j <= i; j++) {
			
				H[j][i] = prod_interno(u[i+1],u[j],tamSol);
				
				subtract(u[i+1],H[j][i],u[j],tamSol);
				
			}

			H[i+1][i] = norm(u[i+1],tamSol);

			for (I = 0; I < tamSol; I++)
			{
			
				u[i+1][I] = u[i+1][I]/H[i+1][i];
				
			}

			//Algoritmo QR
			for (j = 0; j <= i-1; j++)
			{
			
				double temp = H[j][i];
				
				H[j][i] = c[j]*H[j][i] + s[j]*H[j+1][i];
				
				H[j+1][i] = -s[j]*temp + c[j]*H[j+1][i];
				
			}

			r = sqrt(pow(H[i][i],2) + pow(H[i+1][i],2));
			
			c[i] = H[i][i]/r;
			
			s[i] = H[i+1][i]/r;
			
			H[i][i] = r;
			
			H[i+1][i] = 0;
			
			e[i+1] = -s[i]*e[i];
			
			e[i] = c[i]*e[i];
			
			rol = fabs(e[i+1]);
			
			i++;
			
		}
		
		i--;

		for (j= i; j >= 0; j--) {
		
			int l;
			
			double sum = 0;
			
			for (l = j+1; l <= i; l++) {

				sum += H[j][l]*y[l];
				
			}
			
			y[j] = (e[j]-sum)/H[j][j];
			
		}

		for (j = 0; j <= i; j++){
		
			for (I = 0; I < tamSol; I++) {
			
				x[I] += y[j]*u[j][I];
				
			}
		}

		iter++;
		
	} while ((rol >= epsilon) && (iter < maxIter));

	free(e); 
	free(c);
	free(s); 
	free(y); 

	for (I = 0; I < nVetBase + 1; I++){
	
		free(u[I]);
		
		free(H[I]);
		
	}
	
	free(u); free(H);

	printf("Número de iterações: %d\n",iter);

	return x;
}
