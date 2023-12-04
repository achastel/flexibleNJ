#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
//#include <iostream.h>

struct p_node {
	long int d;
	int i,j;
};

struct otuName {
   char *name;
   int ordem;
};

typedef struct p_node p_node_type;


int father(int x){
    return x/2; }

void heap_fica(p_node_type vet[], int i, int qtde){
   int g,j=i;
   p_node_type aux;
   while( (2*j)<qtde ){
      g=2*j;
      if ((g<qtde) && (vet[g].d > vet[g+1].d))
         g=g+1;
      if (vet[j].d <= vet[g].d)
         j=qtde;
      else {
         aux = vet[j];
         vet[j] = vet[g];
         vet[g] = aux;
         j=g;
      }
   }
   return;
}
       
        
void matrizQ(long int **nm, int sz, long int **d, long int D[]){
   int i,j;
   for(i=0;i<sz;i++){
      #pragma omp parallel for
      for(j=0;j<i;j++){
	        nm[i][j]=((long int)(sz - 2) * d[i][j]) - (D[i] + D[j]);
	 }
   }
   return;
}//matrizQ()

void minimoModificado(p_node_type R[], int size, int nelem, long int **mq, double p){
        // esse metodo obtem uma lista R de n/p  pares de OTUs que possuem menor valor Q
        // carregar a matriz Q em uma lista L = (i, j, Q[i][j]) onde i, j sao indices de Q
        //printf("minimoModificado size=%d n.elem=%d\n",size,nelem);
        p_node_type aux;
        p_node_type *L;
        L=(p_node_type *) malloc((nelem+1)*sizeof(p_node_type));
        //printf("1 ");
        int indices[size];
        int k, lim, j=0, ind=0,x;
        lim = size*p;
//printf("2 ");

        #pragma omp parallel for
        for (int i=0; i<size; i++) indices[i]=0;
        //printf("3 ");
	
        #pragma omp parallel for
	for (int i=0; i<=nelem; i++){
           L[i].i=-1;
           L[i].j=-1;
           L[i].d=0;
	}
        //printf ("\nModificado size = %d p=%f nelem=%d lim=%d size*p=%f\n",size,p,nelem,lim,size*p);
        k=1;
        for (int i=0; i<size; i++)		     //monta uma lista L com as dist. em m
            for (int j=0; j<i; j++){ 
                L[k].i = i; L[k].j=j; L[k++].d=mq[i][j];
            }
        k--;
        // transformar L em um heap
        x=k/2;//floor((k-1)/2);
        //printf("\nX=%d   K=%d L[x].i=%d L[k].i=%d \n",x,k,L[x].i,L[k].i);
                
        for (int h=x; h>0; h--){
           heap_fica(L,h,k);
        }
        //printf("heap montado \n");
        //for (int h=1; h<k; h++) printf("L[%d]=(%d,%d)=%ld \n",h,L[h].i,L[h].j,L[h].d);
                
        // obter a lista de nos l[u] = True se o no u foi escolhido; False caso contrario.
        if (lim < 1) lim=1 ;
        //printf("\np=%f lim=%d   size=%d  nelem=%d  k=%d\n",p,lim,size,nelem,k);//<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        int i=0;	
        if(lim==1){
           R[0]=L[1];}
           //printf("\nR[0]=(%d,%d)=%ld\n",R[0].i,R[0].j,R[0].d);}
        else{
           i=0;
           while((i<lim)&&(nelem>0)){					//selecioa os menores o topo do heap 
              //printf("i=%d nelem=%d L[1].i=%d L[1].j=%d \n",i,nelem,L[1].i,L[1].j);
              if(indices[L[1].i]==0 && indices[L[1].j]==0){
                 R[i++]=L[1];
                 indices[L[1].i]=1; 
                 indices[L[1].j]=1; }
              /*remover (i+1) do heap*/
              //printf("0,");
              L[1]=L[nelem];
              nelem--;
              //printf("1,");
              heap_fica(L, 1, nelem);
              //printf("2\n");
           }//while()
        }
        //printf("\nMinimos selecionados= ");
        //for (int y=0;y<lim;y++) printf("R[%d]=(%d,%d)=%ld , ",y,R[y].i, R[y].j,R[y].d);
        //printf("\n");
        if(L) free(L);
        return;
}//minimo modificado


// Driver Code
int main(int argc, char *argv[])
{
    // *************************************ABRINDO ARQUIVO COM MATRIZ DE DISTÂNCIA
	FILE* ptr = fopen(argv[1], "r");
	double f = strtod(argv[2],NULL);

	if (ptr == NULL) {
		printf("no such file.");
		return 0;
	}
	
	int  size, i=0,j=0,k=1, ind=0;
	fscanf(ptr,"%d", &size);									// ****LENDO No. OTUs
	printf("N. OTUs=%d f=%f\n",size,f);
	float s = 0.0;
	s = (size/2);
	s = (size+1)*s;										//******DETERMINANDO TAMANHO MATRIZ TRIANGULAR
	//printf("%.2f\n",s);
	int t=s-size;
     
     struct otuName otus[size];
	long int number, **m;										// MATRIZ DE DISTANCIAS
	m=(long int**)malloc(size*sizeof(long int));
	long int D[size];											// SOMA DAS DISTANCIAS DE UMA OTU(linha)
	int z=0;
	
	long int *Q[size];
	for(int i=1;i<size;i++)
	   Q[i]= (long int *) malloc(i*sizeof(long int));
	
	float x= size*f;
	int limite = size*f;
	if(limite==0 && x>0.0) limite=1;
	p_node_type h[limite];
	int menor,maior;
	long int mij, Di, Dj;
    
     printf("Tam. matriz triangular t=%d\n",t);     		
	for(i=0;i<size;i++){										// LAÇO ALOCAÇÃO LINHAS DA MATRIZ E
	       													// INICIALIZAÇÃO DE VALORES
	   otus[i].name = (char *) malloc(8*sizeof(char));
	   D[i]=0;
	   m[i]=(long int*)malloc((i+1)*sizeof(long int));
	   for(j=0;j<i;j++) {m[i][j]=-1; z++;}
	   sprintf(otus[i].name,"%d",i);
	   otus[i].ordem=i;
	}
	printf("Matriz Triangular:\n");
	i=0;k=1;
	while (i<size){											// LENDO DISTANCIAS DO ARQUIVO...
		j=0;
		while(j<=i){
	   		fscanf(ptr, "%ld ", &m[i][j]);		   //..inserindo na matriz
	   		j++;
	   	}
	   	i++;
	}
	//printf("\n Nomes- INICIAL z=%d k=%d\n",z,k);
	//for(int k=0;k<size;k++) printf("%d-%s\n",otus[k].ordem,otus[k].name);
	//    for(i=0;i<size;i++){
     //       for(j=0;  j<=i;   j++)
     //          printf("%ld ",m[i][j]);
     //       printf("\n");
     //    }
     //    printf("\n");	   	
	//printf("\n Lista Prioridade\n");
     //for(k=1;k<t-size;k++) printf("k=%d -> %ld ->[%d,%d]\n",k,prior[k].d, prior[k].i, prior[k].j);

     for(i=0;i<size;i++){										// calculando soma das distancias das OTUS - D[i]
        for(j=0;  j<i;   j++) D[i]=D[i]+m[i][j];
        for(j=i+1;j<size;j++) D[i]=D[i]+m[j][i];
     }
     //printf("Soma das Distacias na liha i:\n");
     //for(i=0;i<size;i++) printf("D[%d] = %ld \n",i,D[i]);

	//**************************************** LAÇO PRINCIPAL
	clock_t tempo=clock();
	
	while(size>=2){
        if (size==2){
	  //char aux[(strlen(otus[0].name)+strlen(otus[1].name)+3)];
          // 
          // if(otus[0].ordem<otus[1].ordem)
	  //       sprintf(aux,"(%s,%s)",otus[0].name,otus[1].name);
	  //    else
	  //       sprintf(aux,"(%s,%s)",otus[1].name,otus[0].name);
	  //    //sprintf(otus[0].name,"%s",aux);
	  //    printf("ÁRVORE = %s\n",aux); //otus[0].name);
	      size--;//exit(1);
	  //    //free(aux);
        }
        else
        {
           s = 0.0;
           s = size/2.0;
   	      s = s*(size+1);					//******DETERMINANDO TAMANHO MATRIZ TRIANGULAR
   	      limite = size*f;					// número de pares a serem selecionados
  	      if (limite==0) limite=1;
   	      t = s-size;
           //printf("Calculando Matriz Q %d\n",size);
	      matrizQ(Q, size, m, D);								// CALCULANDO MATRIZ Q
	      //printf("Selecionando %d menores pares\n",limite);
  	      minimoModificado(h,size, t, Q, f);					// Seleciona Menores

  	      //for(int z=1;z<size;z++) if(Q[z]) free(Q[z]);
  	      //printf("Juntando Pares de OTU's\n");
	      for(int c=0; c<limite; c++) {
              if (size > 2) {
                 double d_i_novo,d_j_novo=0.0;
                 double t1, t2;
              
                 if (h[c].i > h[c].j) {
	               mij=m[h[c].i][h[c].j]; Di=D[h[c].j]; Dj=D[h[c].i];
                    menor=h[c].j;
                    maior=h[c].i; } 
                 else {
                    mij = m[h[c].j][h[c].i]; Di = D[h[c].i]; Dj = D[h[c].j];
                    menor=h[c].i;
                    maior=h[c].j; }

                 //char *aux=(char *) malloc((strlen(otus[menor].name)+strlen(otus[maior].name))*sizeof(char));
                 int x = (strlen(otus[menor].name)+ strlen(otus[maior].name)+3)%8;
                 if (x>0) x = 8-x; 
                 char aux[(strlen(otus[menor].name)+strlen(otus[maior].name)+3)];
                 
                 if(otus[menor].ordem  < otus[maior].ordem)
                    sprintf(aux,"(%s,%s)",otus[menor].name,otus[maior].name);
                 else{
                    sprintf(aux,"(%s,%s)",otus[maior].name,otus[menor].name);
                    otus[menor].ordem=otus[maior].ordem;}
                    
                 //printf("\n%ld ",(strlen(otus[menor].name)+strlen(otus[maior].name)+3));
                 //printf("name = %s tamenho=%ld menor=%d\n",otus[menor].name,strlen(otus[menor].name),menor);
                 //printf("name = %s tamenho=%ld\n",otus[maior].name,strlen(otus[maior].name));
                 //printf("aux = %s tamanho=%ld\n",aux,strlen(aux));
                 
                 if (otus[menor].name != NULL) {
                    char *temp;
                    temp = otus[menor].name;
                    otus[menor].name=NULL;
                    //printf("LIBEROU*****\n");
                 }
                 //else printf("Não Liberou*******\n\n");
                 //if(strlen(aux)>strlen(otus[menor].name)){  printf("R ");
                    int y;
                    y = 8 - (strlen(aux)%8) + 1;
                    //printf("%d ",y);
                   // otus[menor].name=(char*)realloc(otus[menor].name,(strlen(aux)+y));
                   if(!(otus[menor].name=(char*) malloc((strlen(aux)+y)*sizeof(char))))
                      printf("Não alocou*************\n");
                   //else
                   //   printf("ALOCOU************\n");

	            //if (strlen(aux)%8>0) aux[strlen(aux)+1]='\x0';
	            sprintf(otus[menor].name,"%s",aux);
	            //otus[menor].name[strlen(otus[menor].name)+1]='\x0';
	            
	            //free(aux);
	            //otus[menor].name=aux;
                 //printf("name2= %s tamanho=%ld \n",otus[menor].name,strlen(otus[menor].name));
                 
                 //printf("menor=%d maior=%d size=%d\n",menor,maior,size);
                 //printf("\n Nomes\n");
	            //for(int k=0;k<size;k++) printf("%d- %s\n",otus[k].ordem,otus[k].name);
                 //printf("Calculando distancias da nova OTU para as demais.\n"); 
                 t1= (size - 2) * 2;
                 t2= (Di - Dj);
                 t2 = t2/t1;
                 //printf("t1 = %f  t2=%f \n",t1,t2);
                 d_i_novo = (0.5*mij);
                 d_i_novo = d_i_novo + t2 ;
                 d_j_novo = mij - d_i_novo;
		 
                 //calculando distancia do novo aos demais ---- calculando a linha-coluna
                 for(int k=0;k<size;k++){
                    if ((k!=menor)&&(k!=maior))
                       if (k<menor){ 
                          m[menor][k] = ((m[menor][k]+m[maior][k])-mij)*0.5;
                          //printf("1-(%d,%d)=((%ld+%ld)-%ld)*0,5=%ld\n",menor,k,m[menor][k],m[maior][k],mij,m[menor][k]);
                       }
                       else
                          if (k<maior){
                             m[k][menor]=((m[k][menor]+m[maior][k])-mij)*0.5;
                           //printf("2-(%d,%d)=((%ld+%ld)-%ld)*0,5=%ld\n",k,menor,m[k][menor],m[maior][k],mij,m[k][menor]);
                          }
                          else{
                             t1=(((m[k][menor]+m[k][maior])-mij));
                             t1=t1*.5;
                             //printf("3-(%d,%d)=((%ld+%ld)-%ld)*0,5=%f\n",k,menor,m[k][menor],m[k][maior],mij,t1);
                             m[k][menor]=t1;
                             t1=0.0;
                          }
                 }
                 
                 //printf("Eliminando linha-coluna com cópia.\n");
                 for(int k=0;k<size;k++) 
                    if (k>maior) m[k][maior]=m[size-1][k];
                    else m[maior][k]=m[size-1][k];
                 char *temp;
                 temp=otus[maior].name;
                 otus[maior].name = (char *) malloc(strlen(otus[size-1].name)*sizeof(char));
                 strcpy(otus[maior].name,otus[size-1].name);
                 otus[maior].ordem=otus[size-1].ordem;
                 
                 //printf("Corrigindo coordenadas.\n");
                 for(int e=c+1; e<limite; e++){
                    if(h[e].i==size-1)
                       if (h[e].j > maior) {
                          h[e].i=h[e].j;
                          h[e].j=maior;}
                      else h[e].i=maior;
                 
                    if (h[e].j==size-1)
                       h[e].j =maior;
                 }   
              
                 //printf("Recalculando somas das linhas.\n\n");
		 #pragma omp parallel for
                 for(int k=0; k<size-1; k++) D[k]=0;
		 #pragma omp parallel for
                 for(int k=0; k<size-1; k++) {
                    for (int l = 0; l<(size-(size-k)); l++) {
                       D[k]=D[k]+m[k][l];
                       D[l]=D[l]+m[k][l];}
                 }
                 
                //for(int k=0; k<size-1; k++) printf("D[%d]=%ld\n",k,D[k]);
              }//if (size > 2)	
	         
	         size--;
	         //printf("\n++++++>>>>>>>>>>>SIZE=%d\n",size);
	      }//for
	      //free(h);
	   }
	   
	}//while(size>1)
     float tempo2 = (clock() - tempo) / (double)CLOCKS_PER_SEC;
     printf("Tempo(s) = %f \n",tempo2);
    
	return 0;
}


