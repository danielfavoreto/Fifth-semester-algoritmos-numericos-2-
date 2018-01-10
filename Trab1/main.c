#include <stdio.h>
#include <time.h>
#include "Vector.h"

#define M_E  2.71828182845904523536
#include <math.h>
#include <float.h>
Vector resolverPVCValidacao1(Vector b, double alpha, double omega, int m, int n,double k, int iterationMax, int valorPrescritoAbaixo, int valorPrescritoAcima, int valorPrescritoDireita, 
                             int valorPrescritoEsquerda);
Vector resolverPVCValidacao2(Vector b, double alpha, double omega, int m, int n,double k, int iterationMax, int valorPrescritoAbaixo, int valorPrescritoAcima, int valorPrescritoDireita, 
                             int valorPrescritoEsquerda);
Vector resolverPVCAplicacao1(Vector b, double alpha, double omega, int m, int n,double k, int iterationMax, int valorPrescritoAbaixo, int valorPrescritoAcima, int valorPrescritoDireita, 
                             int valorPrescritoEsquerda);

int main(void) {

    // validacao 1
    
    int n = 40,m = 250 ;

    clock_t begin, end;
    double time_spent;

    begin = clock();

    /* descomente a linha para executar a validacao 1
    Vector bValidacao1 = buildVectorWithValue(m*n,0);
    
    Vector vetorResultado1 = resolverPVCValidacao1(bValidacao1,0.00001,1.6,m,n,1.0,10000,2,2,2,2);
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
     printf ("Tempo execucao: %f\n",time_spent);
    
    FILE * temp = fopen("RESULTADO.m", "w");

    int i;

    fprintf(temp,"vetor = [");
    for (i=0; i < m*n; i++)
    {
        fprintf(temp, "%lf ", vetorResultado1.v[i]); 
    }
    fprintf(temp,"]\n");
    fprintf(temp, "matriz = reshape(vetor, %d, %d)\n", m, n);
    fprintf(temp, "surfc(matriz)\n");
    fprintf(temp, "pause\n");
    fclose(temp);

    deleteVector(bValidacao1);
    deleteVector(vetorResultado1);
    */


    /* validacao 2

    Vector xSolucaoExata = buildVectorWithValue(m*n,0);
    Vector bValidacao2 = buildVectorWithValue(m*n,0);

    double hx = (1.0 - 0.0)/(n - 1);
    double hy = (1.0 - 0.0)/(m - 1);
    double gradiente = 0.0;
    double derivadaParcialUx = 0.0;
    double derivadaParcialUy = 0.0;
    double maxDiferenca = 0.0;
    double diferencaAtual = 0.0;

    for (int i = 0; i < m*n; i++){

        double x1 = (i%n) * hx;            
        double y1 = (i/n) * hy;
        double u = 10*x1*y1*(1-x1)*(1-y1)*pow(M_E,pow(x1,4.5));

        xSolucaoExata.v[i] = u;

        gradiente = -(pow(M_E,pow(x1,4.5)) * (90*pow(x1,4.5)*(y1-1)*y1 + 247.5*(x1-1)*pow(x1,3.5)*(y1-1)*y1 + 202.5*(x1-1)*pow(x1,8)*(y1-1)*y1 + 20*(x1-1)*x1 + 20*(y1-1)*y1));

        derivadaParcialUx = -45*pow(M_E,pow(x1,4.5)) * ((-pow(x1,5.5)) + pow(x1,4.5) - 0.444444*x1 + 0.222222)*(y1-1)*y1;

        derivadaParcialUy = 10*pow(M_E,pow(x1,4.5)) * (x1-1)*x1*((2*y1) - 1);

        bValidacao2.v[i] = gradiente + derivadaParcialUx + 20*y1*derivadaParcialUy + u;

    }
    
    Vector vetorResultado2 = resolverPVCValidacao2(bValidacao2,0.00001,1.6,m,n,1.0,10000,0,0,0,0);

    for (int i = 0; i < m*n; i++){
        
        diferencaAtual = fabs(vetorResultado2.v[i] - xSolucaoExata.v[i]);
        
        if (diferencaAtual > maxDiferenca){

            maxDiferenca = diferencaAtual;
            
        }

    }

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    printf ("Erro Validacao 2: %f\n Tempo execucao: %f\n", maxDiferenca,time_spent);

    FILE * temp = fopen("RESULTADO.m", "w");

    int i;

    fprintf(temp,"vetor = [");
    for (i=0; i < m*n; i++)
    {
        fprintf(temp, "%lf ", vetorResultado2.v[i]); 
    }
    fprintf(temp,"]\n");
    fprintf(temp, "matriz = reshape(vetor, %d, %d)\n", m, n);
    fprintf(temp, "surfc(matriz)\n");
    fprintf(temp, "pause\n");
    fclose(temp);
    deleteVector(vetorResultado2);
*/ 

    Vector bAplicacao1 = buildVectorWithValue(m*n,70);
    Vector vetorResultado3 = resolverPVCAplicacao1(bAplicacao1,0.00001,1.6,m,n,1.0,1000,70,70,0,200);

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
     printf ("Tempo execucao: %f\n",time_spent);
    
    FILE * temp = fopen("RESULTADO.m", "w");

    int i;

    fprintf(temp,"vetor = [");
    for (i=0; i < m*n; i++)
    {
        fprintf(temp, "%lf ", vetorResultado3.v[i]); 
    }
    fprintf(temp,"]\n");
    fprintf(temp, "matriz = reshape(vetor, %d, %d)\n", m, n);
    fprintf(temp, "surfc(matriz)\n");
    fprintf(temp, "pause\n");
    fclose(temp);

    deleteVector(bAplicacao1);
    deleteVector(vetorResultado3);
    return(0);

}

// Aplicacao 1
Vector resolverPVCAplicacao1(Vector b, double alpha, double omega, int m, int n,double k, int iterationMax, int valorPrescritoAbaixo, int valorPrescritoAcima, int valorPrescritoDireita, 
                             int valorPrescritoEsquerda) {

    double ax = 0.0,bx = 1.0,cx = 0.0,dx = 1.0;

    int i,j,index;
    int iterationNumber = 0;
    int valor_prescrito_na_linha_de_cima = 1;
    int valor_prescrito_na_linha_de_baixo = 1;
    int valor_prescrito_na_linha_da_esquerda = 1;
    int valor_prescrito_na_linha_da_direita = 0;
    
    double maxDiferenca = -1;
    double maxValor = -1;
    double erro = 0.0;
    double valorAtualizado = 0.0;
    double diferencaAtual = 0.0;
    int nao_convergir = 0;
    double derivadaNormalY = 0.0;
    double derivadaNormalX = 0.0;
    double hx = (bx - ax)/(n - 1);
    double hy = (dx - cx)/(m - 1);
    double novoA1 = 0.0;
    double kSobreHyQuadrado = (-k/(hy*hy));
    
    Vector X = buildVectorWithValue(m*n,0);
    Vector xAntigo = buildVectorWithValue(m*n,FLT_MAX);
        
    // linha baixo
    if (valor_prescrito_na_linha_de_baixo) {
        
        for (i = 0; i < n; i++) {
            
            X.v[i] = valorPrescritoAbaixo;
            b.v[i] = valorPrescritoAbaixo;
            
        }
    }

    // linha cima
    if (valor_prescrito_na_linha_de_cima) {
        
        for (i = 0; i < n; i++) {
            
            // ????
            index = m*n + (i - n);
            
            X.v[index] = valorPrescritoAcima;
            b.v[index] = valorPrescritoAcima;
            
        }
    }
    // linha esquerda
       if (valor_prescrito_na_linha_da_esquerda) {
        
        for (i = 0; i < m; i++) {

			index = i*n;

            X.v[index] = valorPrescritoEsquerda;
            b.v[index] = valorPrescritoEsquerda;
        }
    }
    // linha direita    
       if (valor_prescrito_na_linha_da_direita) {
        
        for (i = 1; i < m - 1; i++) {

            index = i*n + (n-1);

            X.v[index] = valorPrescritoDireita;
            b.v[index] = valorPrescritoDireita;
            
        }
    }    
//    showVector(b);
    
    double a1, b1, c1, d1, e1;

    printf ("Aplicação 1:\n");

    // validacao 2
    a1 = 1.0 + (2*k * (1.0/(hx*hx) + (1.0/(hy*hy))));
    // e1 sendo calculado dentro do for
    c1 = (-k/(hx*hx)) - (1.0/(2*hx));
    b1 = (-k/(hx*hx)) + (1.0/(2*hx));
    d1 = 0.0;
    e1 = 0.0;
    // d1 sendo calculado dentro do for

             
    // fim validacao 2

    do {

        iterationNumber++;
        
        // canto inferior esquerdo 
        if (!valor_prescrito_na_linha_da_esquerda && !valor_prescrito_na_linha_de_baixo) {

            e1 = kSobreHyQuadrado;
            d1 = kSobreHyQuadrado;
            // validacao 1 e validacao 2
            novoA1 = a1 + c1 + e1;
             
            valorAtualizado = (e1*(-derivadaNormalY)*hy) + c1*(-derivadaNormalX)*hx;
             
            X.v[0] = X.v[0] + omega * (((b.v[0] - valorAtualizado) - (b1*X.v[1] + d1*X.v[n]))/novoA1 - X.v[0]); 
            
        }
        
        // linha de baixo
        if (!valor_prescrito_na_linha_de_baixo) {
            
             e1 = kSobreHyQuadrado;
             d1 = kSobreHyQuadrado;

             novoA1 = a1 + e1;
             valorAtualizado = e1*(-derivadaNormalY)*hy;
            
             for (j = 1; j < n-1; j++) {

                //index = j;
             
                X.v[j] = X.v[j] + omega * (((b.v[j] - valorAtualizado) - (c1*X.v[j - 1] + b1*X.v[j + 1] + d1*X.v[j + n]))/novoA1 - X.v[j]);
           
             }
        }
        
        // canto inferior direito 
        if (!valor_prescrito_na_linha_da_direita && !valor_prescrito_na_linha_de_baixo) {
            
            index = n-1;

            e1 = kSobreHyQuadrado;
            d1 = kSobreHyQuadrado;

            // validacao 1 e validacao 2
            novoA1 = a1 + b1 + e1;
             
            valorAtualizado = (e1*(-derivadaNormalY)*hy) + b1*(derivadaNormalX)*hx;
             
            X.v[index] = X.v[index] + omega * ( ((b.v[index] - valorAtualizado) - (c1*X.v[index - 1] + d1*X.v[index + n]))/novoA1 - X.v[index]);
            
        }
        
        for (i = 1; i < m - 1; i++){
                        
               double parcela = (20*(i*hy)/(2*hy));

               d1 = kSobreHyQuadrado + parcela;
               e1 = kSobreHyQuadrado - parcela;

            // linha esquerda 
            if (!valor_prescrito_na_linha_da_esquerda){
                
                index = i*n;
               

                novoA1 = a1 + c1;
                valorAtualizado = c1*(-derivadaNormalX)*hx;
                
                X.v[index] = X.v[index] + omega * (((b.v[index] - valorAtualizado) - (e1*X.v[index - n] + b1*X.v[index + 1] + d1*X.v[index + n]))/novoA1 - X.v[index]);

            }
 

           
            for (j = 1; j < n-1; j++) {
                
                index = i*n + j;

               X.v[index] = X.v[index] + omega * ((b.v[index] - (e1*X.v[index - n] + c1*X.v[index - 1] + b1*X.v[index + 1] + d1*X.v[index + n]))/a1 - X.v[index]);


            }
            
            // linha direita
            if (!valor_prescrito_na_linha_da_direita){
                
                novoA1 = a1 + b1-(b1*hx);
//                novoA1 = a1 + b1;

                valorAtualizado = b1*70*hx; 
//                valorAtualizado = b1*derivadaNormalX*hx;
								
                index = i*n + (n-1);

                X.v[index] = X.v[index] + omega * (((b.v[index] - valorAtualizado) - (e1*X.v[index - n] + c1*X.v[index - 1] + d1*X.v[index + n]))/novoA1 - X.v[index]);

            }
			
        }
        
            d1 = kSobreHyQuadrado + (20*((m-1)*hy)/(2*hy));
            e1 = kSobreHyQuadrado - (20*((m-1)*hy)/(2*hy));                


        // canto superior esquerdo 
        if (!valor_prescrito_na_linha_da_esquerda && !valor_prescrito_na_linha_de_cima) {
            
            index = m*n - n;

            // validacao 1 e validacao 2
            novoA1 = a1 + c1 + d1;  
             
            valorAtualizado = (d1*derivadaNormalY*hy) + c1*(-derivadaNormalX)*hx;
             
            X.v[index] = X.v[index] + omega * ( ((b.v[index] - valorAtualizado) - (b1*X.v[index + 1] + e1*X.v[index - n]))/novoA1 - X.v[index]); 
            
        }
        
        
            
        // linha cima
        if (!valor_prescrito_na_linha_de_cima) {

            for (j = m*n - (n-1); j < m*n - 1; j++){        
    
                novoA1 = a1 + d1;
                valorAtualizado = d1*derivadaNormalY*hy;
                  
                X.v[j] = X.v[j] + omega * (((b.v[j] - valorAtualizado) - (e1*X.v[j - n] + b1*X.v[j + 1] + c1*X.v[j -1]))/novoA1 - X.v[j]);
            }

        }
        
        // canto superior direito 
        if (!valor_prescrito_na_linha_da_direita && !valor_prescrito_na_linha_de_cima) {
            
            index = m*n - 1;
            
            // validacao 1 e validacao 2
            novoA1 = a1 + b1 + d1;
             
            valorAtualizado = (d1*derivadaNormalY*hy) + b1*derivadaNormalX*hx;
             
            X.v[index] = X.v[index] + omega * ( ((b.v[index] - valorAtualizado) - (c1*X.v[index - 1] + e1*X.v[index - n]))/novoA1 - X.v[index]); 
            
        }
 
        for (i = 0; i < m*n; i++){
            
            diferencaAtual = fabs(X.v[i] - xAntigo.v[i]);
            
            if (diferencaAtual > maxDiferenca){

                maxDiferenca = diferencaAtual;
                
            }
            
            if (fabs(X.v[i]) > maxValor){
                
                maxValor = fabs(X.v[i]);
                
            }

          xAntigo.v[i] = X.v[i];   

        }
        
        erro = maxDiferenca/maxValor;
    
        if (erro > alpha)
            nao_convergir = 1;
        else 
            nao_convergir = 0;

        maxDiferenca = -1;
        maxValor = -1;
    
    } while (nao_convergir && iterationNumber < iterationMax); 

    deleteVector(xAntigo);

    printf ("Numero iteracoes: %d\n ", iterationNumber);

    return X;
}


// validacao 2
Vector resolverPVCValidacao2(Vector b, double alpha, double omega, int m, int n,double k, int iterationMax, int valorPrescritoAbaixo, int valorPrescritoAcima, int valorPrescritoDireita, 
                             int valorPrescritoEsquerda) {

    double ax = 0.0,bx = 1.0,cx = 0.0,dx = 1.0;

    int i,j,index;
    int iterationNumber = 0;
    int valor_prescrito_na_linha_de_cima = 1;
    int valor_prescrito_na_linha_de_baixo = 1;
    int valor_prescrito_na_linha_da_esquerda = 1;
    int valor_prescrito_na_linha_da_direita = 1;
    
    double maxDiferenca = -1;
    double maxValor = -1;
    double erro = 0.0;
    double valorAtualizado = 0.0;
    double diferencaAtual = 0.0;
    int nao_convergir = 0;
    double derivadaNormalY = 0.0;
    double derivadaNormalX = 0.0;
    double hx = (bx - ax)/(n - 1);
    double hy = (dx - cx)/(m - 1);
    double novoA1 = 0.0;
    double kSobreHyQuadrado = (-k/(hy*hy));
    
    Vector X = buildVectorWithValue(m*n,0);
    Vector xAntigo = buildVectorWithValue(m*n,FLT_MAX);
        
    // linha baixo
    if (valor_prescrito_na_linha_de_baixo) {
        
        for (i = 0; i < n; i++) {
            
            X.v[i] = valorPrescritoAbaixo;
            b.v[i] = valorPrescritoAbaixo;
            
        }
    }

    // linha cima
    if (valor_prescrito_na_linha_de_cima) {
        
        for (i = 0; i < n; i++) {
            
            index = m*n + (i - n);
            
            X.v[index] = valorPrescritoAcima;
            b.v[index] = valorPrescritoAcima;
            
        }
    }
    // linha esquerda
       if (valor_prescrito_na_linha_da_esquerda) {
        
        for (i = 0; i < m; i++) {

			index = i*n;

            X.v[index] = valorPrescritoEsquerda;
            b.v[index] = valorPrescritoEsquerda;
        }
    }
    // linha direita    
       if (valor_prescrito_na_linha_da_direita) {
        
        for (i = 1; i < m - 1; i++) {

            index = i*n + (n-1);

            X.v[index] = valorPrescritoDireita;
            b.v[index] = valorPrescritoDireita;
            
        }
    }    
    printf ("Validacao 2:\n");
//    showVector(b);
    
    double a1, b1, c1, d1, e1;

    // validacao 2
    a1 = 1.0 + (2*k * (1.0/(hx*hx) + (1.0/(hy*hy))));
    // e1 sendo calculado dentro do for
    c1 = (-k/(hx*hx)) - (1.0/(2*hx));
    b1 = (-k/(hx*hx)) + (1.0/(2*hx));
    d1 = 0.0;
    e1 = 0.0;
    // d1 sendo calculado dentro do for

             
    // fim validacao 2

    do {

        iterationNumber++;
        
        // canto inferior esquerdo 
        if (!valor_prescrito_na_linha_da_esquerda && !valor_prescrito_na_linha_de_baixo) {

            e1 = kSobreHyQuadrado;
            d1 = kSobreHyQuadrado;
            // validacao 1 e validacao 2
            novoA1 = a1 + c1 + e1;
             
            valorAtualizado = (e1*(-derivadaNormalY)*hy) + c1*(-derivadaNormalX)*hx;
             
            X.v[0] = X.v[0] + omega * (((b.v[0] - valorAtualizado) - (b1*X.v[1] + d1*X.v[n]))/novoA1 - X.v[0]); 
            
        }
        
        // linha de baixo
        if (!valor_prescrito_na_linha_de_baixo) {
            
             e1 = kSobreHyQuadrado;
             d1 = kSobreHyQuadrado;

             novoA1 = a1 + e1;
             valorAtualizado = e1*(-derivadaNormalY)*hy;
            
             for (j = 1; j < n-1; j++) {

                //index = j;
             
                X.v[j] = X.v[j] + omega * (((b.v[j] - valorAtualizado) - (c1*X.v[j - 1] + b1*X.v[j + 1] + d1*X.v[j + n]))/novoA1 - X.v[j]);
           
             }
        }
        
        // canto inferior direito 
        if (!valor_prescrito_na_linha_da_direita && !valor_prescrito_na_linha_de_baixo) {
            
            index = n-1;

            e1 = kSobreHyQuadrado;
            d1 = kSobreHyQuadrado;

            // validacao 1 e validacao 2
            novoA1 = a1 + b1 + e1;
             
            valorAtualizado = (e1*(-derivadaNormalY)*hy) + b1*(derivadaNormalX)*hx;
             
            X.v[index] = X.v[index] + omega * ( ((b.v[index] - valorAtualizado) - (c1*X.v[index - 1] + d1*X.v[index + n]))/novoA1 - X.v[index]);
            
        }
        
        for (i = 1; i < m - 1; i++){
                        
               double parcela = (20*(i*hy)/(2*hy));

               d1 = kSobreHyQuadrado + parcela;
               e1 = kSobreHyQuadrado - parcela;

            // linha esquerda 
            if (!valor_prescrito_na_linha_da_esquerda){
                
                index = i*n;
               

                novoA1 = a1 + c1;
                valorAtualizado = c1*(-derivadaNormalX)*hx;
                
                X.v[index] = X.v[index] + omega * (((b.v[index] - valorAtualizado) - (e1*X.v[index - n] + b1*X.v[index + 1] + d1*X.v[index + n]))/novoA1 - X.v[index]);

            }
 
            for (j = 1; j < n-1; j++) {
                
                index = i*n + j;

               X.v[index] = X.v[index] + omega * ((b.v[index] - (e1*X.v[index - n] + c1*X.v[index - 1] + b1*X.v[index + 1] + d1*X.v[index + n]))/a1 - X.v[index]);


            }
            
            // linha direita
            if (!valor_prescrito_na_linha_da_direita){
                
                novoA1 = a1 + b1;

                valorAtualizado = b1*derivadaNormalX*hx;
								
                index = i*n + (n-1);

                X.v[index] = X.v[index] + omega * (((b.v[index] - valorAtualizado) - (e1*X.v[index - n] + c1*X.v[index - 1] + d1*X.v[index + n]))/novoA1 - X.v[index]);

            }
			
        }
        
            d1 = kSobreHyQuadrado + (20*((m-1)*hy)/(2*hy));
            e1 = kSobreHyQuadrado - (20*((m-1)*hy)/(2*hy));                

        // canto superior esquerdo 
        if (!valor_prescrito_na_linha_da_esquerda && !valor_prescrito_na_linha_de_cima) {
            
            index = m*n - n;

            // validacao 1 e validacao 2
            novoA1 = a1 + c1 + d1;  
             
            valorAtualizado = (d1*derivadaNormalY*hy) + c1*(-derivadaNormalX)*hx;
             
            X.v[index] = X.v[index] + omega * ( ((b.v[index] - valorAtualizado) - (b1*X.v[index + 1] + e1*X.v[index - n]))/novoA1 - X.v[index]); 
            
        }
                
        // linha cima
        if (!valor_prescrito_na_linha_de_cima) {

            for (j = m*n - (n-1); j < m*n - 1; j++){        
    
                novoA1 = a1 + d1;
                valorAtualizado = d1*derivadaNormalY*hy;
                  
                X.v[j] = X.v[j] + omega * (((b.v[j] - valorAtualizado) - (e1*X.v[j - n] + b1*X.v[j + 1] + c1*X.v[j -1]))/novoA1 - X.v[j]);
            }

        }
        
        // canto superior direito 
        if (!valor_prescrito_na_linha_da_direita && !valor_prescrito_na_linha_de_cima) {
            
            index = m*n - 1;
            
            // validacao 1 e validacao 2
            novoA1 = a1 + b1 + d1;
             
            valorAtualizado = (d1*derivadaNormalY*hy) + b1*derivadaNormalX*hx;
             
            X.v[index] = X.v[index] + omega * ( ((b.v[index] - valorAtualizado) - (c1*X.v[index - 1] + e1*X.v[index - n]))/novoA1 - X.v[index]); 
            
        }
 
        for (i = 0; i < m*n; i++){
            
            diferencaAtual = fabs(X.v[i] - xAntigo.v[i]);
            
            if (diferencaAtual > maxDiferenca){

                maxDiferenca = diferencaAtual;
                
            }
            
            if (fabs(X.v[i]) > maxValor){
                
                maxValor = fabs(X.v[i]);
                
            }

          xAntigo.v[i] = X.v[i];   

        }
        
        erro = maxDiferenca/maxValor;
    
        if (erro > alpha)
            nao_convergir = 1;
        else 
            nao_convergir = 0;

        maxDiferenca = -1;
        maxValor = -1;
    
    } while (nao_convergir && iterationNumber < iterationMax); 

    deleteVector(xAntigo);

    printf ("Numero iteracoes: %d\n ", iterationNumber);

    return X;
}




// validacao 1
Vector resolverPVCValidacao1(Vector b, double alpha, double omega, int m, int n,double k, int iterationMax, int valorPrescritoAbaixo, int valorPrescritoAcima, int valorPrescritoDireita, 
                             int valorPrescritoEsquerda) {

    double ax = 0.0,bx = 1.0,cx = 0.0,dx = 1.0;

    int i,j,index;
    int iterationNumber = 0;
    int valor_prescrito_na_linha_de_cima = 1;
    int valor_prescrito_na_linha_de_baixo = 1;
    int valor_prescrito_na_linha_da_esquerda = 1;
    int valor_prescrito_na_linha_da_direita = 1;
    
    double maxDiferenca = -1;
    double maxValor = -1;
    double erro = 0.0;
    double valorAtualizado = 0.0;
    double diferencaAtual = 0.0;
    int nao_convergir = 0;
    double derivadaNormalY = 0.0;
    double derivadaNormalX = 0.0;
    double hx = (bx - ax)/(n - 1);
    double hy = (dx - cx)/(m - 1);
    double novoA1 = 0.0;
    
    
    Vector X = buildVectorWithValue(m*n,0);
    Vector xAntigo = buildVectorWithValue(m*n,FLT_MAX);
        
    // linha baixo
    if (valor_prescrito_na_linha_de_baixo) {
        
        for (i = 0; i < n; i++) {
            
            X.v[i] = valorPrescritoAbaixo;
            b.v[i] = valorPrescritoAbaixo;
            
        }
    }

    // linha cima
    if (valor_prescrito_na_linha_de_cima) {
        
        for (i = 0; i < n; i++) {
            
            // ????
            index = m*n + (i - n);
            
            X.v[index] = valorPrescritoAcima;
            b.v[index] = valorPrescritoAcima;
            
        }
    }
    // linha esquerda
       if (valor_prescrito_na_linha_da_esquerda) {
        
        for (i = 0; i < m; i++) {

			index = i*n;
            X.v[index] = valorPrescritoEsquerda;
            b.v[index] = valorPrescritoEsquerda;
        }
    }
    // linha direita    
       if (valor_prescrito_na_linha_da_direita) {
        
        for (i = 1; i < m - 1; i++) {

            index = i*n + (n-1);
            X.v[index] = valorPrescritoDireita;
            b.v[index] = valorPrescritoDireita;
            
        }
    }    
    printf ("Validacao 1\n");
//    showVector(b);
    
    double a1, b1, c1, d1, e1;

    // validacao 1
    a1 = (2*k * (1.0/(hx*hx) + (1.0/(hy*hy))));
    e1 = (-k/(hy*hy));
    c1 = (-k/(hx*hx));
    b1 = (-k/(hx*hx));
    d1 = (-k/(hy*hy));
    printf ("------------------VALORES-\n");
    printf ("%f,%f,%f,%f,%f \n",a1,b1,c1,d1,e1);
             
    // fim validacao 1

    
    
    do {

		
        iterationNumber++;
        
        // canto inferior esquerdo 
        if (!valor_prescrito_na_linha_da_esquerda && !valor_prescrito_na_linha_de_baixo) {
            
            
            // validacao 1 e validacao 2
            novoA1 = a1 + c1 + e1;
             
			
            valorAtualizado = (e1*(-derivadaNormalY)*hy) + c1*(-derivadaNormalX)*hx;
             
            X.v[0] = X.v[0] + omega * (((b.v[0] - valorAtualizado) - (b1*X.v[1] + d1*X.v[n]))/novoA1 - X.v[0]); 
            
        }
        
        // linha de baixo
        if (!valor_prescrito_na_linha_de_baixo) {
            
             novoA1 = a1 + e1;
             valorAtualizado = e1*(-derivadaNormalY)*hy;
            
             for (j = 1; j < n-1; j++) {

                //index = j;
             
                X.v[j] = X.v[j] + omega * (((b.v[j] - valorAtualizado) - (c1*X.v[j - 1] + b1*X.v[j + 1] + d1*X.v[j + n]))/novoA1 - X.v[j]);
           
             }
        }
        
        // canto inferior direito 
        if (!valor_prescrito_na_linha_da_direita && !valor_prescrito_na_linha_de_baixo) {
            
            index = n-1;
            
            // validacao 1 e validacao 2
            novoA1 = a1 + b1 + e1;
             
            valorAtualizado = (e1*(-derivadaNormalY)*hy) + b1*(derivadaNormalX)*hx;
             
            X.v[index] = X.v[index] + omega * ( ((b.v[index] - valorAtualizado) - (c1*X.v[index - 1] + d1*X.v[index + n]))/novoA1 - X.v[index]);
            
        }
        
        for (i = 1; i < m - 1; i++){
                        
            // linha esquerda 
            if (!valor_prescrito_na_linha_da_esquerda){
                
                index = i*n;
                
                novoA1 = a1 + c1;
                valorAtualizado = c1*(-derivadaNormalX)*hx;
                
                X.v[index] = X.v[index] + omega * (((b.v[index] - valorAtualizado) - (e1*X.v[index - n] + b1*X.v[index + 1] + d1*X.v[index + n]))/novoA1 - X.v[index]);

            }
            
            for (j = 1; j < n-1; j++) {
                
               index = i*n + j;

               X.v[index] = X.v[index] + omega * ((b.v[index] - (e1*X.v[index - n] + c1*X.v[index - 1] + b1*X.v[index + 1] + d1*X.v[index + n]))/a1 - X.v[index]);


            }
            
            // linha direita
            if (!valor_prescrito_na_linha_da_direita){
                
                novoA1 = a1 + b1;
                valorAtualizado = b1*derivadaNormalX*hx;
								
                index = i*n + (n-1);

                X.v[index] = X.v[index] + omega * (((b.v[index] - valorAtualizado) - (e1*X.v[index - n] + c1*X.v[index - 1] + d1*X.v[index + n]))/novoA1 - X.v[index]);

            }
			
        }
        
                
        // canto superior esquerdo 
        if (!valor_prescrito_na_linha_da_esquerda && !valor_prescrito_na_linha_de_cima) {
            
            index = m*n - n; 
            
            // validacao 1 e validacao 2
            novoA1 = a1 + c1 + d1;  
             
            valorAtualizado = (d1*derivadaNormalY*hy) + c1*(-derivadaNormalX)*hx;
             
            X.v[index] = X.v[index] + omega * ( ((b.v[index] - valorAtualizado) - (b1*X.v[index + 1] + e1*X.v[index - n]))/novoA1 - X.v[index]); 
            
        }
        
        
            
        // linha cima
        if (!valor_prescrito_na_linha_de_cima) {
            
            for (j = m*n - (n-1); j < m*n - 1; j++){        
    
                novoA1 = a1 + d1;
                valorAtualizado = d1*derivadaNormalY*hy;
  
                
                X.v[j] = X.v[j] + omega * (((b.v[j] - valorAtualizado) - (e1*X.v[j - n] + b1*X.v[j + 1] + c1*X.v[j -1]))/novoA1 - X.v[j]);
            }
        }
        
        // canto superior direito 
        if (!valor_prescrito_na_linha_da_direita && !valor_prescrito_na_linha_de_cima) {
            
            index = m*n - 1;
            
            // validacao 1 e validacao 2
            novoA1 = a1 + b1 + d1;
             
            valorAtualizado = (d1*derivadaNormalY*hy) + b1*derivadaNormalX*hx;
             
            X.v[index] = X.v[index] + omega * ( ((b.v[index] - valorAtualizado) - (c1*X.v[index - 1] + e1*X.v[index - n]))/novoA1 - X.v[index]); 
            
        }
 
        for (i = 0; i < m*n; i++){
            
            diferencaAtual = fabs(X.v[i] - xAntigo.v[i]);
            
            if (diferencaAtual > maxDiferenca){

                maxDiferenca = diferencaAtual;
                
            }
            
            if (fabs(X.v[i]) > maxValor){
                
                maxValor = fabs(X.v[i]);
                
            }

          xAntigo.v[i] = X.v[i];   

        }
        
        erro = maxDiferenca/maxValor;
         printf ("erro validacao 1: %f\n",erro);
    
        if (erro > alpha)
            nao_convergir = 1;
        else 
            nao_convergir = 0;
    
        maxDiferenca = -1;
        maxValor = -1;
    
    } while (nao_convergir && iterationNumber < iterationMax); 

    deleteVector(xAntigo);

    printf ("Numero iteracoes: %d\n ", iterationNumber);

    return X;
}
