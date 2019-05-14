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




/*
Matrix_ptr get_nxn_matrix(int n, Matrix_ptr A){
  Matrix_ptr n_matrix = matrix_alloc(n,n);
  double val;
  for(int i = 0; i<n; i++){
    for(int j = 0; j<n; j++){
      val = get_matrix_value(A,i,j);
      set_matrix_value(n_matrix, i,j,val);
    }
  }

  return n_matrix;
}



*/


double u_function(double x, double t, double a){

	if(x <= t){
		return 1.0;
	}
	else if((x > t) &(x<=a)){
		return (a-x)/(a-t);
	}
	else{
		return 0.0;
	}

}

void pde_burgers(int t_run){
	double lambda = 0;
	double t = 0;
	double x_min = -0.5;
	double x_max  = 2.0;
	double dx = 0.1;
	double dt = 0.25;
	double a = 0.8;

	double x_run = 2 + (x_max - x_min)/dx;
	double Ui, Ui_1, U_k;
	Matrix_ptr U =  matrix_alloc(x_run, t_run);
	Matrix_ptr X =  matrix_alloc(x_run,1);
	Matrix_ptr T =  matrix_alloc(t_run,1);


	printf("%d\n", U->m);


	for (int i = 0; i < X->m; i++) {
		/* code */
		set_matrix_value(X, i+1	,0, x_min + (dx*i));
	}
	//print_matrix(X);

	for (int i = 0; i < T->m; i++) {
		/* code */
		set_matrix_value(T, i	,0, (dt*i));
	}

	//set initial condition
	for (int i = 0; i < X->m; i++) {
		/* code */
		set_matrix_value(U,i,0,u_function(get_matrix_value(X,i,0), 0, a));
	}

	//print_matrix(U);
	double U_new_i;
	for (int k = 0; k < T->m -1; k++) {
	  /* code */
	  //space
	  for(int i = 1; i <  X->m; i++) {
	    /* code */
			//double Ui  = get_matrix_value(U,i,k);
			double Ui = u_function(get_matrix_value(X,i,0), get_matrix_value(T,k,0), a);
			//double U_i = get_matrix_value(U,i-1,k);
			double U_i = u_function(get_matrix_value(X,i-1,0),get_matrix_value(T,k,0),a);
			double u_x = u_function(get_matrix_value(X,i,0),get_matrix_value(T,k,0),a);

	    //double U_new_i = Ui - (dt/dx) *u_x*(Ui-U_i);
			double U_new_i = Ui - (dt/dx)* (1/2) * ((Ui*Ui)-(U_i*U_i));
			set_matrix_value(U, i,k+1, U_new_i);
			//print_matrix(U);
	  }
	}
	print_matrix(U);






}




int main(){

	//create matrix to represent the plate
	pde_burgers(5);

}
