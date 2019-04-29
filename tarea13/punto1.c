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






Matrix_ptr calculate_cholesky(Matrix_ptr A){
  int n = A->m;
  Matrix_ptr L = matrix_alloc(n,n);
  double val, sum,a;
  double l_jk,l_ik, l_jj;
  for(int i = 0; i<n; i++){
    for(int j = 0; j<n; j++){

      //diagonal
      //printf("Coordenada (%d,%d)\n",i, j );
      if(i==j){
        l_ik = 0;
        l_jk = 0;
        a = 0;
        l_jj = 1;

        //get a_ii
        a = get_matrix_value(A,j,j);
        //summatory of l_jk^2
        sum = 0;
        //printf("%s\n","************************" );
        //print_matrix(L);
        for(int k = 0;k < j;k++){
          l_jk = get_matrix_value(L,j,k);
          sum = sum + pow(l_jk,2);
        }
        //l_jj = sqrt(a_jj - sum(l_jk^2))
        val = sqrt(a-sum);
        //printf("%f,%f,%f,%f\n", l_jj,a,l_ik,l_jk);
        //printf("val : %f\n", val);
        //printf("%s\n","************************" );
        set_matrix_value(L, i,j,val);

      }
      //lower triangular matrix
      else if(i > j){
        l_ik = 0;
        l_jk = 0;
        a = 0;
        l_jj = 1;

        a = get_matrix_value(A,i,j);
        l_jj = get_matrix_value(L,j,j);
        sum = 0;
        //printf("%s\n","************************" );
        //print_matrix(L);
        for(int k =0;k <j;k++){

          l_ik = get_matrix_value(L,i,k);
          l_jk = get_matrix_value(L,j,k);

          sum = sum + (l_ik*l_jk);
        }


        val = (1/l_jj)*(a-sum);
        //printf("%f,%f,%f,%f\n", l_jj,a,l_ik,l_jk);
        //printf("val : %f\n", val);
        set_matrix_value(L, i,j,val);
        //printf("%s\n","************************" );

      }

    }
  }
  return L;
}



Matrix_ptr solve_linear(Matrix_ptr L, Matrix_ptr b){
  Matrix_ptr x = matrix_alloc(L->m, 1);
  Matrix_ptr y = matrix_alloc(L->m, 1);
  Matrix_ptr Lt = matrix_alloc(L->m, L->n);
  int i = 0;
  int n = L->m;
  int j = 0;
  double c = 0;
  //gsl_matrix_set_zero(x);
  //gsl_matrix_set_zero(y);

  //Sustitucion progresiva
  for (i = 0; i < n; i++ ){
      c = get_matrix_value(b, i, 0);
      for (j = 0; j < i; j ++ ){
        c = c - get_matrix_value(L, i, j)*get_matrix_value(y, j, 0);
      }
      c = c/get_matrix_value(L, i, i);
      set_matrix_value(y, i, 0, c);
  }
  //gsl_matrix_transpose_memcpy (Lt, L);
  Lt = matrix_transpose(L);

  //Sustitucion regresiva
  for (i = n-1; i >= 0; i-- ){
    c = get_matrix_value(y, i, 0);
    for (j = n-1; j > i; j--){
      c =  c - get_matrix_value(Lt, i, j)*get_matrix_value(x, j, 0);
    }
    c = c/get_matrix_value(Lt, i, i);
    set_matrix_value(x, i, 0, c);
  }
  return x;
}


Matrix_ptr tridiagonal_solution(Matrix_ptr A, Matrix_ptr b){
	int m = A->m;
	int n =A->n;

	Matrix_ptr L = calculate_cholesky(A);

	Matrix_ptr R = solve_linear(L, b);

	return R;
}


void temperature(){
	//initial parameters
	double lambda;
	double dx = 2; //cm
	double dt = 0.1; //s
	double alpha = 0.835; //cm2/s

	double f_0 = 100;
	double f_m = 50;
	int m = 4;
	int n = 4;

	/////////////////////////////////////////
	lambda = alpha * (dt /(dx*dx));

	Matrix_ptr T_l = matrix_alloc(4,1);
	Matrix_ptr T_new = matrix_alloc(4,1);


	printf("Lambda parameter: %f\n", lambda);
	/////////////////////////////////////////

	Matrix_ptr T = matrix_alloc(4,4);
	printf("%s\n", "initial matrix" );
	for(int i = 0;i< m;i++){
		for(int j = 0;j< n;j++){
			if(i == j){
				set_matrix_value(T,i,j,2*(1+lambda));

			}
			else if((i+1 == j) || (i-1 == j)){
				set_matrix_value(T,i,j,-lambda);
			}

		}
	}
	/////////////


	int simulation = 100;
	for(int time = 0;time < simulation ;time ++){
		//set vector
		//initial border
		printf("%s\n", "");
		printf("%s\n", "");
		printf("%s\n", "");
		double t1 = get_matrix_value(T_l,0,0);
		double t2 = get_matrix_value(T_l,1,0);
		double val = 2*lambda*f_0 + 2*(1-lambda)*t1 + lambda *t2;
		set_matrix_value(T_new, 0,0, val);

		double t3;
		//inner points
		for(int i = 1;i < n-1;i++){
			t1 = get_matrix_value(T_l,i-1,0);
			t2 = get_matrix_value(T_l,i ,0);
			t3 = get_matrix_value(T_l,i+1,0);
			val = lambda*t1 + 2*(1-lambda)*t2 + lambda*t3;
			set_matrix_value(T_new, i,0, val);
		}
		//final border
		double tm = get_matrix_value(T_l,3,0);
		double tm_1 = get_matrix_value(T_l,2,0);
		val = 2*lambda*f_m + 2*(1-lambda)*tm + lambda*tm_1;
		set_matrix_value(T_new, 3,0, val);


		printf("temperature at each point at time:   %f s\n ", dt*time );
		//print_matrix(T_new);

		
		T_l = T_new;
		

		T_l = tridiagonal_solution(T,T_l);
		print_matrix(T_l);


	}
		


	
	/////////////////////////////////////////

	/////////////////////////////////////////

	/////////////////////////////////////////

}


int main(){

	//create matrix to represent the plate
	temperature();

}

