#include<stdio.h>
#include<math.h>
#include<gsl/gsl_blas.h>
void print_matrix(const gsl_matrix *m) {
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			printf("%.20f\t", gsl_matrix_get(m, i, j));


		}
		printf("\n");
	}
}


void range(gsl_matrix *A){
	printf("ENTROOOO\n");
	double e = 1E-10;
	int range = A->size1;
	int cero = 0;
	for(int i = 0;i<A->size1 ;i++){
		cero = 0;
		for(int j = 0;j<A->size2;j ++){
			
			if(gsl_matrix_get(A,i,j) < e){
				cero ++;			
			}
			
		}
		if(cero ==A->size2){
			range--;		
			
		}

	}
	
	printf("RANGO de la matriz es: %d \n",range);


}

void calculate_range( gsl_matrix *A,gsl_matrix *B){
	printf("Calculate Range");

	int i = 0;
	int j = 0;
	double c,m;
	int cont=0;
	int counter = 0;
	int Range  = A->size1 ;
	gsl_vector* v = gsl_vector_alloc(A->size1);
	gsl_vector* v1 = gsl_vector_alloc(A->size1);
	gsl_vector* column =gsl_vector_alloc(A->size2);
	gsl_vector* column1 =gsl_vector_alloc(A->size2);
		
	for (int row = 0; row < Range; row++) {
		printf("Interation i = %d\n", i);
		for (int col = 0; col < A->size2; col++) {
			printf("interation j %d\n ", col);
			//printf("%f\t", gsl_matrix_get(m, i, j));
			if(col<row){	
				if (gsl_matrix_get(A,col,col) == 0){
					
					//caso1: evaluar elemento no nulo en la columna
					int caso = 0;
					//recorrer columna 
					for(int i = row ; i < A->size1; i++ ){		
						printf("Interacion (%d,%d)\n", i, j);	
						if(gsl_matrix_get(A,i,col)){		
							//swap rows
							printf("------Caso1------\n");
							//extract zero row

							gsl_matrix_swap_rows(A,i,col);
							caso = 1;
							printf("row swap (%d,%d)\n", i, col);	
							print_matrix(A);
							col = -1;
							row =-1;
							break;
						}
					}
					if(caso == 0){
						printf("------Caso2------\n");
						//printf("pivot zero\n");
						//swap columns
						

						if(col= A->size2 -1){
							break;
							} 
						for(int j = col+1; j < Range; j++){
							if(gsl_matrix_get(A,row-1,j)){
								gsl_matrix_swap_columns(A, j, col);
							}
							
							Range --;
							printf("column swap (%d,%d)\n", row, col);
							print_matrix(A);


						}
						
						

					}
					
						
					
				}else{	
					if(row== A->size1 ){
						break;
					} 
					if(col== A->size2 ){
						break;
					} 
					
					printf("caso normal (%d,%d)\n", row, col);
					c = gsl_matrix_get(A,row,col)/gsl_matrix_get(A,col,col);
					//printf("constatn %f/n", c);
					for(int k =0;k< A->size2;k++){
						m = gsl_matrix_get(A,row,k) - c*gsl_matrix_get(A,col,k);
						gsl_matrix_set(A,row,k,m);				
					}
					print_matrix(A);
				}			
			}
		printf("###########################################\n");
		//print_matrix(A);
		//printf("interation j %d\n ", j);
		printf("###########################################\n");

		}
		printf("\n");
	}
	
	//revisar filas nulas

	printf("Calculando el rango......\n");
	range(A);	
	
}






int main(){
	
	/*
	const unsigned int M=5,N=5;
	
	gsl_matrix *A = gsl_matrix_alloc(M,N);
	
	for(int i=0; i < A->size1; i++){
	
		gsl_matrix_set(A,i,i,1);	
	}		

	for(int i=0; i<A->size1; i++){
		
		for(int j=0; j<A->size2; j++){
			printf("%f\t", gsl_matrix_get(A, i,j));
		}
		printf("\n");				
	}

	
	*/
	

	int m,n;            // Número de filas y columnas de la matriz de usuario
  	printf("Introduzca el número de filas de la matriz: ");
  	scanf("%d",&m);
	const unsigned int M = m;

  	printf("Introduzca el número de columnas de la matriz: ");
  	scanf("%d",&n);
	const unsigned int N = n;

	//const double rcond = 1E-15;
	gsl_matrix *A = gsl_matrix_alloc(M,N);  
	gsl_matrix *A_new = gsl_matrix_alloc(M,N);
	printf("Introduzca los elementos de la matriz A(%dx%d)\n",m,n);
	double elem;
	for(int i=0; i<M; i++){
    		for(int j=0; j<N; j++){
      			printf("a%d,%d: ", i+1,j+1);
      			scanf("%lf", &elem);
			gsl_matrix_set(A, i, j, elem);
    		}
  	}

	print_matrix(A);

	calculate_range(A,A_new);
	
	print_matrix(A);
	return 0;





}
