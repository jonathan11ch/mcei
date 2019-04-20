//eigen values and eigen vectors
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


typedef struct EigenStruct{
  Matrix_ptr v;
  double lambda;
}EigenStruct, *EigenStruct_ptr;


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






Matrix_ptr scalar_mult(Matrix_ptr M, double s){

	int m = M->m;
	int n = M->n;
	double val;
	Matrix_ptr Mult = matrix_alloc(m,n);

	for(int i =0; i < m ;i++){
		for(int j= 0; j < n ;j++){
			val = get_matrix_value(M,i,j) * s;
			set_matrix_value(Mult,i,j,val);
		}
	}

	return Mult;
}



//calculate F norm
double matrix_norm(Matrix_ptr M){
	double val = 0;
	double a;
	int n = M->n;
	int m = M->m;

	for(int i= 0; i<m;i++){
		for(int j = 0 ;j<n; j++){
			a = get_matrix_value(M,i,j);
			val = val + fabs(pow(a,2));
		}

	}

	val = sqrt(val);
	return val;
}

Matrix_ptr get_identity_matrix(int m, int n){

	Matrix_ptr I = matrix_alloc(m,n);

	for(int i= 0; i<m;i++){
		for(int j = 0 ;j<n; j++){
			if (i ==j){
				set_matrix_value(I,i,j,1);
			}
		}

	}

	return I;

}


Matrix_ptr sum_matrix(Matrix_ptr A, Matrix_ptr B){

	int n = A->n;
	int m = A->m;
	double a,b,c;
	Matrix_ptr C = matrix_alloc(m,n);

	for(int i= 0; i<m;i++){
		for(int j = 0 ;j<n; j++){
			a = get_matrix_value(A,i,j);
			b = get_matrix_value(B,i,j);
			c = a+b;
			set_matrix_value(C,i,j,c);
		}
	}

	return C;


}

double get_out_condition(Matrix_ptr A){
	double err;
	int m =A->m;
	int n =A->n;
	double val;

	for(int i = 0; i < m; i++){
		for (int j = 0; j < n; j++)
		{
			/* code */
			if(i!=j){
				val = get_matrix_value(A, i,j);
				val = pow(val,2.0);
				err += val;	
			}
			
		}
	}
	return err;
}

int jacobi_clasic_method(Matrix_ptr A, double e){
	int m =A->m;
	int n =A->n;
	
	int i_n,j_n;
	double N = 0;
	double err;
	/////////
	double Tan2, Cos2, Cos, Sen;
	double a_pq, a_pp,a_qq;
	/////////
	int cont = 0;
	
	//
	Matrix_ptr A_new;
	Matrix_ptr V = get_identity_matrix(m,n);

	do{
		printf("%s\n","-------------------------------------------------" );
		//get elimination element
		for(int i = 0; i < m; i++){
			for (int j = 0; j < n; j++)
			{
			/* code */
				if(j>i){
					double val = get_matrix_value(A,i,j);
					if(fabs(val)> N){
						N = fabs(val);
						i_n = i;
						j_n = j;
					}
				}
			}
			
		}
		N = 0;
		printf("Indexxx %d,%d\n", i_n, j_n);
		a_pq = get_matrix_value(A,i_n, j_n);
		a_pp = get_matrix_value(A,i_n, i_n);
		a_qq = get_matrix_value(A,j_n, j_n);

		printf("valores: %f,%f,%f\n",a_pq, a_pp, a_qq);
		//get rotation sin and cos
		if(a_pp != a_qq){
			printf("%s\n","caso 1" );
			Tan2 = 2*a_pq/(a_pp - a_qq);
			Cos2 = 1.0/(sqrt(pow(Tan2,2.0) + 1 ));
			Cos = sqrt(0.5*(1+Cos2));
			Sen	= sqrt(0.5*(1-Cos2));		
		}
		else{
			printf("%s\n","caso 2" );
			Cos = cos(3.1416/4);
			Sen = sin(3.1416/4);
		}
		
		//get J and Jt template
		Matrix_ptr J = get_identity_matrix(m,n);
		Matrix_ptr J_t = get_identity_matrix(m,n);
		//build rotation matrix
		set_matrix_value(J, i_n,i_n,  Cos);
		set_matrix_value(J, j_n,j_n,  Cos);
		set_matrix_value(J, i_n,j_n, -Sen);
		set_matrix_value(J, j_n,i_n,  Sen);
		//J_transpose
		set_matrix_value(J_t, i_n,i_n,  Cos);
		set_matrix_value(J_t, j_n,j_n,  Cos);
		set_matrix_value(J_t, i_n,j_n,  Sen);
		set_matrix_value(J_t, j_n,i_n, -Sen);

		//printf("%s\n","J matrix");
		//print_matrix(J);
		//printf("%s\n","J transpose matrix");
		//print_matrix(J_t);
		V = mult_matrix(V,J);

		//operate iteration
		Matrix_ptr A_new = mult_matrix(J_t, A);
		A_new = mult_matrix(A_new,J);
		printf("%s\n","····A_new····" );
		print_matrix(A_new);
		printf("%s\n","····A_new····" );
		//calculate out condition
		err = get_out_condition(A_new);
		printf("Error obtained %f\n", err );
		A = A_new;
		cont ++;
		printf("%s\n","-------------------------------------------------" );
	}	
	while(err > e);
	printf("%s\n", "######RESULTS#######" );
	printf("%s\n", "######A#######" );
	print_matrix(A);
	printf("%s\n", "######V#######" );
	print_matrix(V);
	printf("%s\n", "##############" );

	printf("---------------Valores Propios -------------\n");
	for(int i= 0;i<m ;i++){
		double val = get_matrix_value(A, i,i);
		printf("Lamba(%d): %f\n",i, val);
	}
	printf("---------------Vectores Propios -------------\n");
		
	for(int i= 0;i<m ;i++){
		printf("V(%d):\n",i );
		for (int j = 0; j < n; j++)
		{
			/* code */
			double val = get_matrix_value(V, i,j);
		    printf("%f\t", val);
		}
		printf("\n");
	}
	printf("Cantidad de Iteraciones: %d\n", cont);
	printf("Cantidad de Iteraciones: %d\n", cont);
	printf("--------------------------------------------\n");

	return cont;
	
}


int jacobi_cyclic_method(Matrix_ptr A, double e){
	int m =A->m;
	int n =A->n;
	
	int i_n,j_n;
	double N = 0;
	double err;
	/////////
	double Tan2, Cos2, Cos, Sen;
	double a_pq, a_pp,a_qq;
	/////////
	int cont = 0;
	
	//
	Matrix_ptr A_new;
	Matrix_ptr V = get_identity_matrix(m,n);
	i_n = 0;
	j_n = 0;

	do{
		printf("%s\n","-------------------------------------------------" );
		//get elimination element
		if(j_n < n -1 ){
			j_n++;
		}
		else if((i_n < m -2)){
			i_n ++;
			j_n = i_n + 1;
		}else{
			i_n = 0;
			j_n = 1;
		}

		printf("Indexxx %d,%d\n", i_n, j_n);
		a_pq = get_matrix_value(A,i_n, j_n);
		a_pp = get_matrix_value(A,i_n, i_n);
		a_qq = get_matrix_value(A,j_n, j_n);

		printf("valores: %f,%f,%f\n",a_pq, a_pp, a_qq);
		//get rotation sin and cos
		if(a_pp != a_qq){
			printf("%s\n","caso 1" );
			Tan2 = 2*a_pq/(a_pp - a_qq);
			Cos2 = 1.0/(sqrt(pow(Tan2,2.0) + 1 ));
			Cos = sqrt(0.5*(1+Cos2));
			Sen	= sqrt(0.5*(1-Cos2));		
		}
		else{
			printf("%s\n","caso 2" );
			Cos = cos(3.1416/4);
			Sen = sin(3.1416/4);
		}
		
		//get J and Jt template
		Matrix_ptr J = get_identity_matrix(m,n);
		Matrix_ptr J_t = get_identity_matrix(m,n);
		//build rotation matrix
		set_matrix_value(J, i_n,i_n,  Cos);
		set_matrix_value(J, j_n,j_n,  Cos);
		set_matrix_value(J, i_n,j_n, -Sen);
		set_matrix_value(J, j_n,i_n,  Sen);
		//J_transpose
		set_matrix_value(J_t, i_n,i_n,  Cos);
		set_matrix_value(J_t, j_n,j_n,  Cos);
		set_matrix_value(J_t, i_n,j_n,  Sen);
		set_matrix_value(J_t, j_n,i_n, -Sen);

		//printf("%s\n","J matrix");
		//print_matrix(J);
		//printf("%s\n","J transpose matrix");
		//print_matrix(J_t);
		V = mult_matrix(V,J);

		//operate iteration
		Matrix_ptr A_new = mult_matrix(J_t, A);
		A_new = mult_matrix(A_new,J);
		printf("%s\n","····A_new····" );
		print_matrix(A_new);
		printf("%s\n","····A_new····" );
		//calculate out condition
		err = get_out_condition(A_new);
		printf("Error obtained %f\n", err );
		A = A_new;
		cont ++;
		printf("%s\n","-------------------------------------------------" );
	}	
	while(err > e);
	printf("%s\n", "######RESULTS#######" );
	printf("%s\n", "######A#######" );
	print_matrix(A);
	printf("%s\n", "######V#######" );
	print_matrix(V);
	printf("%s\n", "##############" );

	printf("---------------Valores Propios -------------\n");
	for(int i= 0;i<m ;i++){
		double val = get_matrix_value(A, i,i);
		printf("Lamba(%d): %f\n",i, val);
	}
	printf("---------------Vectores Propios ---el ----------\n");

	for(int i= 0;i<m ;i++){
		printf("V(%d):\n", i);
		for (int j = 0; j < n; j++)
		{
			/* code */
			double val = get_matrix_value(V, i,j);
		    printf(" %f\t", val);
		}
		printf("\n");
	}
	printf("Cantidad de Iteraciones: %d\n", cont);

	printf("--------------------------------------------\n");

	return cont;		
	
}

int main(){
  printf("%s\n", "INGRESE LA MATRIZ");
  Matrix_ptr A = user_request_matrix();
  int a,b;
  double err = 0.000001;
  printf("%s\n", "INGRESE LA PRECISION DESEADA");
  scanf("%lf", &err);

  
  printf("%s\n","------------------ Metodo Ciclico -----------------" );
  a = jacobi_cyclic_method(A, err);
  printf("%s\n","------------------ Metodo Clasico -----------------" );
  b = jacobi_clasic_method(A, err);

  if(a>b){
  	printf("El metodo que menos iteraciones requiere es el método CLASICO (%d).\n", b);
  }else if(a < b){
  	printf("El metodo que menos iteraciones requiere es el método CICLICO (%d).\n",a);
  }

}
