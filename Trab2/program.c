#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "./CommonFiles/protos.h"
#include "GMRES.h"

#define FREE(p){if(p != NULL) free(p);}

double get_time ()
{

	struct timeval tv; gettimeofday(&tv, 0);
	
	return (double)(tv.tv_sec * 100.0 + tv.tv_usec / 10000.0);
	
}


int main (int argc, char* argv[])
{
	double time;
  
	if (argc != 2)
	{
	
		printf("\n Erro! Sem arquivo da matriz (.mtx)"); 
		
		printf("\n Modo de usar: ./program <nome_da_matriz> Saindo... [main]\n\n");
		
		return 0;
		
	}
	
	MAT *A = (MAT*) malloc (sizeof(MAT));						
	
	MATRIX_readCSR (A,argv[1]);

	double *x = (double*)calloc(A->n,sizeof(double));
	
	int i,j;
	
	double *b = (double*)calloc(A->n,sizeof(double));
	
	for (i = 0; i < A->n; i++)
	{
	
		for (j = A->IA[i]; j < A->IA[i+1]; j++)
		{
		
			b[i] += A->AA[j];
			
		}
	}

	double* b_permutado = (double*)calloc(A->n,sizeof(double));		

	int choiceRCM,choiceILU;
	
	printf("Deseja usar o reordenamento RCM? Digite 1 para SIM e 0 para NAO\n");
	
	scanf("%d",&choiceRCM);
	
	/*---------------------------------------------*/
	/*---COMO USAR O REORDENAMENTO RCM-------------*/
	/*---------------------------------------------*/
	
	int *p;									// Vetor de permutação
	int  bandwidth;

	if(choiceRCM)
	{
	
		bandwidth = (int) MATRIX_bandwidth(A);					// Calcula Largura de Banda da matriz original
		
		printf("\n  [ REORDENANDO com RCM ]\n");
		
		printf("  - Largura de Banda inicial : %d\n", bandwidth);
		
		/*---START TIME---------------> */ 
		time = get_time(); 
		
		REORDERING_RCM_opt(A,&p);						// Aplica o reordenamento RCM na matriz A
		
		MATRIX_permutation(A,p); 						// Aplica a permutação em A para trocar linhas e colunas

		for (i = 0; i < A->n; ++i) 
		{
		
			b_permutado[i] = b[p[i]];
			
		}
		/*---FINAL TIME---------------> */ 
		time = (get_time() - time)/100.0;
		
		bandwidth = (int) MATRIX_bandwidth(A);					// Calcula Largura de Banda da matriz reordenada
		
		printf("  - Largura de Banda final   : %d\n", bandwidth);
		
		printf("  - Tempo total              : %.6f sec\n\n", time);
		
	}

	printf("Deseja usar precondicionamento? Digite 1 para SIM e 0 para NAO\n");
	
	scanf("%d",&choiceILU);
	
	/*---------------------------------------------*/
	/*---COMO USAR O ALGORITMO ILUP----------------*/
	/*---------------------------------------------*/

	MAT *L = (MAT*) malloc(sizeof(MAT));						// Alocando matriz L
	
	MAT *U = (MAT*) malloc(sizeof(MAT));						// Alocando matriz U
	
	if (choiceILU) {
	
		SparMAT* mat = (SparMAT*) malloc(sizeof(SparMAT));				// Alocando estruturas para o ILU(p)
		
		SparILU* lu  = (SparILU*) malloc(sizeof(SparILU));
		
		int p;
		
		printf("Digite a quantidade de fill-in desejada, de 0 para cima.\n");
		
		scanf("%d",&p);
		
		printf("\n  [ CALCULANDO PRECONDICIONADOR ILU ]\n");
		
		/*---START TIME---------------> */ 
		
		time = get_time(); 
		
		CSRto_SPARMAT (A,mat);								// Convertendo CSR para estrutura especial
		
		ILUP          (mat,lu,p);							// Algoritmo ILU(p)
		
		SPARILU_toCSR (lu,L,U);								// Convertendo estrutura especial para CSR
		
		/*---FINAL TIME---------------> */ 
		
		time = (get_time() - time)/100.0;
		
		printf("  - Tempo total              : %.6f sec\n", time);
		
		SPARILU_clean (lu);								// Liberando memória da estrutura lu
		
		SPARMAT_clean (mat);								// Liberando memória da estrutura mat
	
		/* L contém a parte estritamente inferior de M / L->D contém a diagonal = 1.0 */
		/* U contém a parte estritamente superior de M / U->D contém a diagonal       */
		//MATRIX_printLU (A,L,U);	
	}
	
	double tol;
	
	printf("Digite a tolerância desejada\n");
	
	scanf("%lf",&tol);
	
	int maxIter;
	
	printf("Digite o número máximo de iterações\n");
	
	scanf("%d",&maxIter);
	
	int nVetBase;
	
	printf("Digite o número de iterações para reinício\n");
	
	scanf("%d",&nVetBase);
	
	time = get_time(); 

	if (choiceRCM)
	{
		if (choiceILU)
		{
		
			x = GMRES(A,L,U,x,b_permutado,tol,nVetBase,maxIter,A->n);
			
		}
		else
		{
		
			x = GMRES_No_Preconditioner(A,x,b_permutado,tol,nVetBase,maxIter);
			
		}
	}
	else {
	
		if (choiceILU)
		{
		
			x = GMRES(A,L,U,x,b,tol,nVetBase,maxIter,A->n);
			
		}
		else
		{
		
			x = GMRES_No_Preconditioner(A,x,b,tol,nVetBase,maxIter);
			
		}
	}
	
	time = (get_time() - time)/100.0;
	
	printf("  - Tempo total              : %.6f sec\n\n", time);

	double* x_permutado = (double*)calloc(A->n,sizeof(double));	
	
	double* r = (double*)calloc(A->n,sizeof(double));

	if (choiceRCM)
	{	
		for (i = 0; i < A->n; ++i)
		{
		
			x_permutado[p[i]] = x[i];
			
		}
		
		residuo(A,x,b_permutado,r,A->n);
		
	}
	else
	{
	
		residuo(A,x,b,r,A->n);
		
	}
	
	FILE* fout = fopen("file.txt","w");
	
	for (j = 0; j < A->n; j++)
	{
	
		fprintf(fout,"%.20lf\n",x[j]);
		
	}
	
	fclose(fout);

	double r_r = norm(r,A->n);
	
	free(r);

	printf("norma residuo: %.20lf\n",r_r);

	if (choiceRCM)
	{
	
		free(p);
		
	}
	
	MATRIX_clean(A);

	if (choiceILU)
	{
	
		MATRIX_clean(L);
		
		MATRIX_clean(U);
		
	}
	else
	{
	
		free(L);
		
		free(U);
		
	}	

	free(x);
	free(x_permutado);
	free(b);
	free(b_permutado);
	
	return 0;
}
