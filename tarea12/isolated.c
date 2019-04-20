#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//structure definition to handle Matrix
typedef struct Matrix{
	double *mat;
	int m;
	int n;
}Matrix,*Matrix_ptr;


typedef struct LDR_Matrix{
	Matrix_ptr L;
	Matrix_ptr R;
	Matrix_ptr D;
}LDR_Matrix, *LDR_Matrix_ptr;



Matrix_ptr matrix_alloc(int m, int n){
	//allocate memory for the matrix structure
	Matrix_ptr M = (Matrix_ptr)malloc(sizeof(Matrix));
	//allocate memory for the matrix
	double *matriz = (double *)calloc(m*n,sizeof(double));
	//assign pointers to the matrix pointer
	M->mat = matriz;
	M->m = m;
	M->n = n;
	return M;
}

//set value val to the matrix A in the position i,j
void set_matrix_value(Matrix_ptr mat, int i, int j, double val){

	*(mat->mat + (j + (i * mat->n))) = val;
}

//get value i,j from the matrix A
double get_matrix_value(Matrix_ptr mat, int i, int j){

	return *(mat->mat + (j + (i * mat->n)));
}

//matrix multiplication A*B
Matrix_ptr mult_matrix(Matrix_ptr A, Matrix_ptr B){

	double val,a,b;
	int m = A->m;
	int n = B->n;

	if(A->n != B->m){
		return -1;
	}
	//allocate memory for resulting matrix
	Matrix_ptr R =  matrix_alloc(m,n);
	//perform dot product
	for(int i = 0; i< m; i++){
		for(int j = 0; j< n; j++){
			val = 0;
			for(int k = 0; k < m;k++){
					a = get_matrix_value(A,i,k);
					b = get_matrix_value(B,k,j);
					val = val + a*b;
					set_matrix_value(R, i,j,val);
			}
		}
	}
	return R;
}


void print_matrix(Matrix_ptr M){

	for (int i = 0; i<M->m; i++){
		for (int j = 0; j<M->n; j++){
			printf("%f\t",get_matrix_value(M, i,j));
		}
		printf("\n");
	}
}
//add matrix

//transpose a matrix
Matrix_ptr matrix_transpose(Matrix_ptr M){
		int n = M->n;
		int m = M->m;
		double val;

		Matrix_ptr Mt = matrix_alloc(n,m);

		for(int i = 0; i < Mt->m; i++){
			for(int j = 0; j < Mt->n; j++){
				val = get_matrix_value(M, j,i);
				set_matrix_value(Mt, i,j, val);
			}
		}

		return Mt;
}



Matrix_ptr user_request_matrix(){
	int m,n;            // Número de filas y columnas de la matriz de usuario
	  	printf("Introduzca el número de filas de la matriz: ");
	  	scanf("%d",&m);
		const unsigned int M = m;

	  	printf("Introduzca el número de columnas de la matriz: ");
	  	scanf("%d",&n);
		const unsigned int N = n;
	Matrix_ptr A = matrix_alloc(M,N);
	printf("Introduzca los elementos de la matriz A(%dx%d)\n",m,n);
	double elem;
	for(int i=0; i<M; i++){
    		for(int j=0; j<N; j++){
      			printf("a%d,%d: ", i+1,j+1);
      			scanf("%lf", &elem);
						set_matrix_value(A, i, j, elem);
    		}
  	}

	printf("%s\n","Matriz ingresada:" );
	print_matrix(A);
	return A;

}

int check_error(Matrix_ptr E, double err){
	int n = E->n;
	int m = E->m;
	double val;
	for(int i = 0; i<m ; i++){
		for(int j = 0; j<n ; j++){

			val = get_matrix_value(E,i,j);
			if(val > err){
				return 0;
			}

		}
	}

	return 1;

}

void liebman_iteration(Matrix_ptr T){
	int n = T->n;
	int m = T->m;
	double t,a,b,c,d,t_ij;
	double e;
	Matrix_ptr E = matrix_alloc(m-1,n-2);
	int cont = 0;
	int error = 0;
	double lambda = 1.5;


	do{
		for(int i = 1; i< m; i++){
			for (int j =1; j < n-1; j++)
			{
				/* iteration method */
				
				//liebman method
				//if is on the isolated edge
				if(i == m-1){
					a = get_matrix_value(T, i, j+1);
					b = get_matrix_value(T, i, j-1);
					c = get_matrix_value(T, i-1, j);
					t = get_matrix_value(T, i,j);
					t_ij = (a+b+(2*c))/4;
				}
				else{
					a = get_matrix_value(T, i-1, j);
					b = get_matrix_value(T, i+1, j);
					c = get_matrix_value(T, i, j+1);
					d = get_matrix_value(T, i, j-1);
					t = get_matrix_value(T, i,j);
					t_ij = (a+b+c+d)/4; 
				}

				//acelleration method
				t_ij = (lambda *t_ij) + (1- lambda)*t;

				e = fabs((t_ij-t)/t_ij);
				set_matrix_value(E,i-1,j-1,e);
				//set value in the matrix
				set_matrix_value(T,i,j,t_ij);
				
			}
		}
		//set isoleated edge temperature



		error  = check_error(E, 0.0001);


		printf("%s\n","############ERROR#############" );
		print_matrix(E);
		printf("%s\n","-------------------------" );
		printf("%s\n","New Matrix " );
		print_matrix(T);
		printf("%s\n","#########################" );
		cont ++;
	}
	while(error == 0);
		
	printf("Numero de Iteraciones: %d\n", cont);



}


int main(){

	//create matrix to represent the plate
	Matrix_ptr T = matrix_alloc(5,5);
	//set initial condition matrix
	set_matrix_value(T,0,0,0);
	set_matrix_value(T,0,1,100);
	set_matrix_value(T,0,2,100);
	set_matrix_value(T,0,3,100);
	set_matrix_value(T,0,4,0);

	set_matrix_value(T,1,0,75);
	set_matrix_value(T,1,1,0);
	set_matrix_value(T,1,2,0);
	set_matrix_value(T,1,3,0);
	set_matrix_value(T,1,4,50);
	
	set_matrix_value(T,2,0,75);
	set_matrix_value(T,2,1,0);
	set_matrix_value(T,2,2,0);
	set_matrix_value(T,2,3,0);
	set_matrix_value(T,2,4,50);
	
	set_matrix_value(T,3,0,75);
	set_matrix_value(T,3,1,0);
	set_matrix_value(T,3,2,0);
	set_matrix_value(T,3,3,0);
	set_matrix_value(T,3,4,50);
	
	set_matrix_value(T,4,0,75);
	set_matrix_value(T,4,1,0);
	set_matrix_value(T,4,2,0);
	set_matrix_value(T,4,3,0);
	set_matrix_value(T,4,4,50);


	printf("Matriz inicial\n");
	print_matrix(T);

	liebman_iteration(T);

}