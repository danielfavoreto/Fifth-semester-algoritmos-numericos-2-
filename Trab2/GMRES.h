#include "CommonFiles/heads.h"

void subtract(double* uc,double beta,double* u,int tamU);
void residuo(MAT* A,double* x,double* Fp,double* u,int tamX);
double norm(double* v,int tamV);
double prod_interno(double* a,double* b,int k);
double* GMRES(MAT* Ap,MAT* L,MAT* U,double* x,double* Fp,double tol,int nVetBase,int maxIter,int tamSol);
double* CSR_MATRIX_ARRAY_Product(MAT* mat,double* arr,int sizeArray);
double* solve_system_L(MAT* L,double* arr,int size);
double* solve_system_U(MAT* U,double* arr,int size);
double* GMRES_No_Preconditioner(MAT* Ap,double* x0,double* b,double tol,int kmax,int lmax);
