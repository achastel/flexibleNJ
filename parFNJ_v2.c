#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

struct p_node {
	long int d;
	int i,j;
};

struct otuName {
   char *name;
   int ordem;
};

typedef struct p_node p_node_type;

// minHeapify() definition
// s: start of array
// t_last: last position in the thread 
void minHeapify(p_node_type h[], int index, int s, int t_last)
{
   int left = 2*index - s + 1; // relative positions in h
   int right = 2*index - s  + 2;
   int min = index;
 
   // Checking whether our left or child element
   // is at right index or not to avoid index error
   if (left > t_last || left < s)
      left = -1;
   if (right > t_last || right < s)
      right = -1; 
   // store left or right element in min if
   // any of these is smaller that its parent
   if (left != -1 && h[left].d < h[index].d)
       min = left;
   if (right != -1 && h[right].d < h[min].d)
       min = right;
 
   // Swapping the nodes
   if (min != index) {
      p_node_type temp = h[min];
      h[min] = h[index];
      h[index] = temp;
      // recursively calling for their child elements
      // to maintain min heap
      minHeapify(h, min, s, t_last);
   }
}

void create_heap(p_node_type h[], int s, int t_last){
  
   for(int i = (s+t_last-1)/2; i >= s; i--){
      minHeapify(h, i, s, t_last);
  }
}

void heap_fica(p_node_type vet[], int i, int qtde){
   int g, min=i;
   p_node_type aux;
   while( (2*min) <= qtde ){
      g = 2*min;
      if ((g<qtde) && (vet[g].d > vet[g+1].d))
         g = g+1;
      if (vet[min].d <= vet[g].d)
         min = qtde;
      else {
         aux = vet[min];
         vet[min] = vet[g];
         vet[g] = aux;
         min = g;
      }
   }
   return;
}
       
//Calculate Q Matrix        
void matrizQ(long int **nm, int sz, long int **d, long int D[]){
   int i,j;
   for(i=0;i<sz;i++){
      #pragma omp  for nowait
      for(j=0;j<i;j++){
	        nm[i][j]=((long int)(sz - 2) * d[i][j]) - (D[i] + D[j]);
	 }
   }
   return;
}//matrizQ()

void minimoModificado(p_node_type R[], int size, int nelem, long int **mq, double p){
   // This method obtains a list R with n/p pairs of OTUs that have the lowest Q values
   // loads the matrix Q into a list L = (i, j, Q[i][j]) where i,j are indices of Q
   p_node_type aux;
   p_node_type *L;
   L=(p_node_type *) malloc((nelem+1)*sizeof(p_node_type));
  
   int *indices;		//used to mark already selected OTUs in pairs
   indices = (int *) malloc(size*sizeof(int));
   int k, lim;
   lim = size*p;

   for (int i=1; i<size; i++)  
     for(int j=0; j<i; j++){
       k = i*(i-1)/2 + j +1;
       L[k].i = i; L[k].j=j; L[k].d=mq[i][j];
     }

   // #pragma omp parallel
   //{
   //#pragma omp  for 
   for (int i=0; i<size; i++) indices[i]=0;
    
   // turn L into a heap L
   // #pragma omp single
   //{ 
   for (int h=nelem/2; h>0; h--){
    heap_fica(L,h,nelem);
   }
                  
   if(lim > 1){
    int i=0;
    //select the smallest ones - the top of the heap
    while((i<lim)&&(nelem>0)){
       if(indices[L[1].i]==0 && indices[L[1].j]==0){
	     R[i++]=L[1];
	     indices[L[1].i]=1; 
	     indices[L[1].j]=1; 
	  }
       //remove (i+1) from heap */
       L[1]=L[nelem];
       nelem--;
     
       heap_fica(L, 1, nelem);
    }//while()
  }
  else
    R[0]=L[1];
  //}
  //}
  if(L) free(L);
  free(indices);			   
  return;
}//minimo modificado


// Driver Code
int main(int argc, char *argv[]){
  
   // OPENING FILE WITH DISTANCE MATRIX
   FILE* ptr = fopen(argv[1], "r");
   double f = strtod(argv[2],NULL);
  
   if (ptr == NULL) {
     printf("no such file.");
     return 0;
   }
	
   int  size, i=0,j=0,k=1, ind=0;
   fscanf(ptr,"%d", &size);	        	//Reading No. OTUs
   printf("N. OTUs=%d f=%f\n",size,f);
   int nelem = size*(size-1)/2;	    //Triangular matrix size
	 
   struct otuName otus[size];
   long int number, **m;				// DISTANCE MATRIX
   m=(long int**)malloc(size*sizeof(long int));
   long int D[size];				// sum of the distances of an OTU(line)
  
   long int *Q[size];				// Matrix Q
   for(int i=1;i<size;i++)
      Q[i]= (long int *) malloc(i*sizeof(long int));

   int th;	//threads	
   float x=size*f;
   int limite = size*f;	//n. of pairs to be selected
   if(limite==0 && x>0.0) limite=1;

   p_node_type h[limite];		//HEAP
   p_node_type L[nelem+1];		// auxiliar List
   int indices[size];			
   
   int menor,maior, first[size],last[size];
   long int mij, Di, Dj;
  
   //printf("Tam. matriz triangular t=%d\n",nelem);     		
   for(i=0;i<size;i++){	// LAÇO ALOCAÇÃO LINHAS DA MATRIZ E INICIALIZAÇÃO DE VALORES
      otus[i].name = (char *) malloc(8*sizeof(char));
      D[i]=0;
      m[i]=(long int*)malloc((i+1)*sizeof(long int));
      for(j=0;j<i;j++) {
        m[i][j]=-1;
      }
      sprintf(otus[i].name,"%d",i);
      otus[i].ordem=i;
   }

   //printf("Matriz Triangular:\n");
   i=0;k=1;
   while (i<size){						// READING FILE DISTANCES...
      j=0;
      while(j<=i){
         fscanf(ptr, "%ld ", &m[i][j]); 	// inserting into matrix
         j++;
      }
      i++;
   }
    
   double tempo  = omp_get_wtime();
   long int soma=0;
   //#pragma omp parallel
   //{
   //int tid = omp_get_thread_num(); 
   //int nth = omp_get_num_threads();
   //#pragma omp parallel for
   for(i=0;i<size;i++){ 	       // calculating the sum of the OTUS distances - D[i]
      soma=0;
      #pragma omp parallel
      {
          #pragma omp  for reduction(+ : D[i])
          for(j=0;  j<i;   j++) D[i]=D[i]+m[i][j];
          #pragma omp  for reduction(+ : soma)
          //#pragma omp for nowait
          for(j=i+1;j<size;j++) soma=soma+m[j][i]; //D[i]=D[i]+m[j][i];
      }
      D[i]=D[i]+soma;
   }
     
   //********************************** MAIN LOOP
   while(size>2){
      limite = size*f;		 // number of pairs to be selected
      if (limite==0) limite=1;
      nelem = size*(size-1)/2;
      //printf("Calculando Matriz Q %d\n",size);
      #pragma omp parallel 
      {
         th = omp_get_num_threads();		// N. threads
         int tam = nelem/th;				// numero de elementos por thread
         int tid = omp_get_thread_num();	// identificação thread
   
         matrizQ(Q, size, m, D);	// CALCULATING Q MATRIX

         #pragma omp for nowait
         for (int i=0; i<size; i++) indices[i]=0;

	    for (int i=1; i<size; i++)  //monta uma lista L com as dist. em m
	       #pragma omp for 
	       for(int j=0; j<i; j++){
	          int k = i*(i-1)/2 + j; 
	          L[k].i = i; L[k].j=j; L[k].d=Q[i][j];
	       }

         //creates local heaps
         #pragma omp for 
         for (int i=0; i < th; i++){
	       first[tid] = tid*tam;
	       last[tid] = first[tid]+tam-1;
	       create_heap(L, first[tid], last[tid]);
         }
         
         //selects smaller among local heaps
         #pragma omp master
	    {	 
	       i = 0;
	       while(i < limite && limite > 1 && nelem > 0){
	          int min = 0;
	          for(int j=1; j<th && j < size ; j++)
	             if (L[first[j]].d < L[first[min]].d) 
	                min = j;
	     
	             if(indices[L[first[min]].i]==0 && indices[L[first[min]].j]==0){
	                h[i]=L[first[min]];
	                indices[L[first[min]].i]=1;
	                indices[L[first[min]].j]=1;
	                i++;
	             }
	             L[first[min]]=L[last[min]];
	             nelem--;
	             last[min]--;
	             minHeapify(L, first[min], first[min], last[min]);
	       }
	       if (limite == 1) h[0]=L[0];
	    }
      }//#pragma

      //construindo arestas da árvore
      for(int c=0; c<limite; c++) {
         if (size > 2) {
           double d_i_novo=0.0, d_j_novo=0.0;
           double t1, t2;
           #pragma omp parallel 
           {
	         #pragma omp single
	         {
	             if (h[c].i > h[c].j) {
	                mij=m[h[c].i][h[c].j]; Di=D[h[c].j]; Dj=D[h[c].i];
	                menor=h[c].j;
	                maior=h[c].i;
	             } 
	             else {
	                mij = m[h[c].j][h[c].i]; Di = D[h[c].i]; Dj = D[h[c].j];
	                menor=h[c].i;
	                maior=h[c].j;
	             }
	             char aux[(strlen(otus[menor].name)+strlen(otus[maior].name)+4)];
	             if(otus[menor].ordem  < otus[maior].ordem)
	                sprintf(aux,"(%s,%s)",otus[menor].name,otus[maior].name);
	             else{
	                sprintf(aux,"(%s,%s)",otus[maior].name,otus[menor].name);
	                otus[menor].ordem=otus[maior].ordem;
	             }
	             free(otus[menor].name);

	             int y;
	             y = 8 - (strlen(aux)%8) + 1;
	             if(!(otus[menor].name=(char*) malloc((strlen(aux)+y)*sizeof(char))))
	                printf("Não alocou*************\n");
	 
	             //sprintf(otus[menor].name,"%s",aux);
	             strcpy(otus[menor].name,aux);   
	             //printf("Calculando distancias da nova OTU para as demais.\n"); 
	             t1= (size - 2) * 2;
	             t2= (Di - Dj);
	             t2= t2/t1;
	             d_i_novo = (0.5*mij);
	             d_i_novo = d_i_novo + t2; //d_i_novo =    DISTANCIAS ARESTAS
	             d_j_novo = mij - d_i_novo; //d_j_novo =
	         }
	         //calculating distance from the new OTU to the others - calculating the row-column
              #pragma omp for 
	            for(int k=0;k<menor;k++){
  		          D[menor] = D[menor] - (m[menor][k]);
		          D[k] = D[k] - (m[menor][k]+m[maior][k]);
		          m[menor][k] = (int)((m[menor][k]+m[maior][k])-mij)*0.5;
		          D[k] = D[k] + m[menor][k];
		          D[menor] = D[menor] + m[menor][k];
	            }
	            #pragma omp for 
	             for(int k=menor+1;k<maior;k++){
                     D[menor] = D[menor]-m[k][menor]; //-----
		           D[k] = D[k] - (m[k][menor]+m[maior][k]); //+m[maior][k]);
		           m[k][menor]= (int)((m[k][menor]+m[maior][k])-mij)*0.5;
		           D[k] = D[k] + m[k][menor];
		           D[menor]=D[menor]+m[k][menor];//---
	             }
	  	        #pragma omp for 
	             for(int k=maior;k<size;k++){
		           D[menor] = D[menor]-m[k][menor]; //-----
		           D[k] = D[k] - (m[k][menor]+m[k][maior]);
		           m[k][menor] = (int)((m[k][menor]+m[k][maior])-mij)*0.5;
		           D[k] = D[k] + m[k][menor];
		           D[menor]=D[menor]+m[k][menor]; //---
	             }

                  #pragma omp single
	             {
	                 D[maior]=D[size-1];
	                 free(otus[maior].name);
	                 otus[maior].name = otus[size-1].name;  
	                 otus[maior].ordem=otus[size-1].ordem;
	    
	                 free(m[maior]);
	                 m[maior]=m[size-1];
	            }
	            #pragma omp for 
                 for(int k=maior+1;k<size-1;k++) 
	              m[k][maior]=m[maior][k];
	  
	            //printf("Corrigindo coordenadas.\n");
	            #pragma omp for   
	            for(int e=c+1; e<limite; e++){
	               if(h[e].i==size-1)
	                  if (h[e].j > maior) {
		                h[e].i=h[e].j;
		                h[e].j=maior;
	                  }
	                  else h[e].i=maior;
	    
	               if (h[e].j==size-1)
	                   h[e].j =maior;
	            }
	      }// omp parallel
        }//if (size > 2)	
        size--;
     }//for(c=0 ...
   }//while(size>1)

   // impress
   char aux[(strlen(otus[0].name)+strlen(otus[1].name)+4)];
   if(otus[0].ordem  < otus[1].ordem)
     sprintf(aux,"(%s,%s)",otus[0].name,otus[1].name);
   else
     sprintf(aux,"(%s,%s)",otus[1].name,otus[0].name);
  
   printf("%s\n",aux);
   printf("Time(seconds) = %f \n",omp_get_wtime()-tempo);
  
  return 0;
}


