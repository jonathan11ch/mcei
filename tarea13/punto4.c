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


Matrix_ptr set_initial_condition(Matrix_ptr T){
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
	
	set_matrix_value(T,4,0,0);
	set_matrix_value(T,4,1,0);
	set_matrix_value(T,4,2,0);
	set_matrix_value(T,4,3,0);
	set_matrix_value(T,4,4,0);

	return T;

}

void temperature_surface(int count){
	int n = 5;
	double dt = 10;
	double dx = 10;
	double k = 0.835;
	double lambda = k*(dt/(dx*dx));
	printf("Lambda constant: %f\n", lambda);
	/////////////////////////////////
	Matrix_ptr S = matrix_alloc(n,n);
	S = set_initial_condition(S);
	printf("%s\n", "intial surface state" );
	print_matrix(S);
	printf("\n");

	///////////create matrix

	Matrix_ptr T = matrix_alloc(3,3);
	printf("%s\n", "initial matrix" );
	for(int i = 0;i< T->m;i++){
		for(int j = 0;j< T->n;j++){
			if(i == j){
				set_matrix_value(T,i,j,2*(1+lambda));

			}
			else if((i+1 == j) || (i-1 == j)){
				set_matrix_value(T,i,j,-lambda);
			}

		}
	}

	print_matrix(T);

	Matrix_ptr T_l = matrix_alloc(3,1);
	double val, t1,t2,t3;
	/////////////////////////////////
	for (int con = 0; con < count; ++con)
	{
		/* code */
		//first direction
		for (int i = 1; i < (S->m)-1; i++)
		{
			for (int j = 1; j < (S->n)-1; j++)
			{
				/* code */
				//printf("%d\n", j);


				t1 = get_matrix_value(S,i-1,j);
				t2 = get_matrix_value(S,i,j);
				t3 = get_matrix_value(S,i+1,j);

				if(j == 1){
					val = lambda*t1 + 2*(1-lambda)*t2 + lambda*t3 + lambda*get_matrix_value(S,i,j-1);
				}else if(j == 3){
					val = lambda*t1 + 2*(1-lambda)*t2 + lambda*t3 + lambda*get_matrix_value(S,i,j+1);

				}else{
					val = lambda*t1 + 2*(1-lambda)*t2 + lambda*t3;
				}
				set_matrix_value(T_l,j-1,0,val);

			}
			printf("%s\n", "fila");
			print_matrix(T_l);
			printf("%s\n", "solucion");
			T_l = tridiagonal_solution(T,T_l);
			print_matrix(T_l);



			//UPDATE S SURFACE
			for (int k = 1;k < 4 ;k++){
				set_matrix_value(S,i,k,get_matrix_value(T_l,k-1,0));
			}

			


		}

		
		for (int j = 1; j< (S->n)-1; j++)
		{
			for (int i = 1; i < (S->m)-1; i++)
			{
			
				//printf("%d\n", j);
				t1 = get_matrix_value(S,i,j-1);
				t2 = get_matrix_value(S,i,j);
				t3 = get_matrix_value(S,i,j+1);

				if(i == 1){
					val = lambda*t1 + 2*(1-lambda)*t2 + lambda*t3 + lambda*get_matrix_value(S,i-1,j);
				}else if(i == 3){
					val = lambda*t1 + 2*(1-lambda)*t2 + lambda*t3 + lambda*get_matrix_value(S,i+1,j);

				}else{
					val = lambda*t1 + 2*(1-lambda)*t2 + lambda*t3;
				}

				//val = lambda*t1 + 2*(1-lambda)*t2 + lambda*t3;
				set_matrix_value(T_l,i-1,0,val);

			}
			printf("%s\n", "fila");
			print_matrix(T_l);
			printf("%s\n", "solucion");
			T_l = tridiagonal_solution(T,T_l);
			print_matrix(T_l);



			//UPDATE S SURFACE
			for (int k = 1;k < 4 ;k++){
				set_matrix_value(S,j,k,get_matrix_value(T_l,k-1,0));
			}

			

		}
		printf("%s\n", "update");
		print_matrix(S);
			//
		printf("%s\n","" );
		
		

	}
	

	//update matriz s





}

int main(){


	temperature_surface(100);

}
