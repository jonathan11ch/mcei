#include <stdio.h>
#include <stdlib.h>


typedef struct Matrix{
	double *mat;
	int m;
	int n;
}Matrix,*Matrix_ptr;

Matrix_ptr matrix_alloc(int m, int n){
	Matrix_ptr M = (Matrix_ptr)malloc(sizeof(Matrix));

	//develop allocation
	double *matriz = (double *)calloc(m*n,sizeof(double));
	M->mat = matriz;
	M->m = m;
	M->n = n;
	return M;
}

void set_matrix_value(Matrix_ptr mat, int i, int j, double val){

	*(mat->mat + (j + (i * mat->n))) = val;
}


double get_matrix_value(Matrix_ptr mat, int i, int j){

	return *(mat->mat + (j + (i * mat->n)));
}


Matrix_ptr mult_matrix(Matrix_ptr A, Matrix_ptr B){

	double val,a,b;
	int m = A->m;
	int n = B->n;
	Matrix_ptr R =  matrix_alloc(m,n);

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


int main()
{
	/* code */
	int n = 3;

	Matrix_ptr A,B,C, At;
	printf("%s\n","#############1###########" );
	A = matrix_alloc(3,2);
	printf("%s\n","##############3##########" );
	B = matrix_alloc(2,4);
	printf("%s\n","#############4###########" );
	//C = matrix_alloc(3,4);
	printf("%s\n","########################" );




	for (int i = 0; i<A->m; i++){
		for (int j = 0;j<A->n; j++){
			set_matrix_value(A, i,j, i*4);
			//set_matrix_value(B, i,j, j);
		}

	}
	for (int i = 0; i<B->m; i++){
		for (int j = 0;j<B->n; j++){
			//set_matrix_value(A, i,j, i*4);
			set_matrix_value(B, i,j, j+1);
		}

	}

	//set_matrix_value(A,1,2,1);
	print_matrix(A);
	printf("%s\n","########################" );
	print_matrix(B);

	//set_matrix_value(A, 2,2, 69);
	C = mult_matrix(A,B);
	printf("%s\n","########################" );
	print_matrix(C);
	printf("%s\n","###########transpose#############" );
	At = matrix_transpose(A);
	print_matrix(At);

	Matrix_ptr Z = user_request_matrix();

	return 0;
}
