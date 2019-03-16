#include <stdio.h>
#include <stdlib.h>

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



LDR_Matrix_ptr LDR_decomposition(Matrix_ptr M){

	int m = M->m;
	int n = M->n;
	double val;
	//allocate memory for the ldr structure
	LDR_Matrix_ptr ldr = (LDR_Matrix_ptr)malloc(sizeof(LDR_Matrix));
	//allocate memory for the three matrices L,R and def
	ldr->L = matrix_alloc(m,n);
	ldr->D = matrix_alloc(m,n);
	ldr->R = matrix_alloc(m,n);

	for(int i = 0; i<m; i++){
		for(int j = 0; j< n; j++){
			val = get_matrix_value(M,i,j);
			if(i == j){
				set_matrix_value(ldr->D,i,j,val);
			}
			else if(i > j){
				set_matrix_value(ldr->L, i,j,val);
			}
			else{
				set_matrix_value(ldr->R,i,j,val);
			}
		}

	}
	return ldr;
}


double det3x3(Matrix_ptr M){
	double a11,a22,a33,a12,a13,a21,a23,a31,a32;
	double det;
	a11 = get_matrix_value(M,1,1);
	a12 = get_matrix_value(M,1,2);
	a13 = get_matrix_value(M,1,3);
	a21 = get_matrix_value(M,2,1);
	a22 = get_matrix_value(M,2,2);
	a23 = get_matrix_value(M,2,3);
	a31 = get_matrix_value(M,3,1);
	a32 = get_matrix_value(M,3,2);
	a33 = get_matrix_value(M,3,3);
	det = (a11*a22*a33) + (a12*a23*a31) + (a13*a21*a32) - ((a31*a22*a13)+(a32*a23*a11)+(a33*a21*a12));
	return det;
}

Matrix_ptr adjoint_3x3(Matrix_ptr M){
	double a11,a22,a33,a12,a13,a21,a23,a31,a32;
	int m = M->m;
	int n = M->n;
	Matrix_ptr Adj = matrix_alloc(m,n);
	a11 = get_matrix_value(M,1,1);
	a12 = get_matrix_value(M,1,2);
	a13 = get_matrix_value(M,1,3);
	a21 = get_matrix_value(M,2,1);
	a22 = get_matrix_value(M,2,2);
	a23 = get_matrix_value(M,2,3);
	a31 = get_matrix_value(M,3,1);
	a32 = get_matrix_value(M,3,2);
	a33 = get_matrix_value(M,3,3);
	set_matrix_value(Adj, 1,1, a22*a33-a32*a23);
	set_matrix_value(Adj, 1,2, a21*a33-a31*a23);
	set_matrix_value(Adj, 1,3, a21*a32-a31*a22);
	set_matrix_value(Adj, 2,1, a12*a33-a32*a13);
	set_matrix_value(Adj, 2,2, a11*a33-a31*a13);
	set_matrix_value(Adj, 2,3, a11*a32-a12*a31);
	set_matrix_value(Adj, 3,1, a12*a23-a22*a13);
	set_matrix_value(Adj, 3,2, a11*a23-a21*a13);
	set_matrix_value(Adj, 3,3, a11*a22-a21*a12);

	return Adj;

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
set_matrix_value(Adj, 1,1,);



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


	LDR_Matrix_ptr ldr =  LDR_decomposition(Z);
	printf("%s\n","########################" );
	print_matrix(ldr->L);
	printf("%s\n","########################" );
	print_matrix(ldr->D);
	printf("%s\n","########################" );
	print_matrix(ldr->R);


	
	return 0;
}
