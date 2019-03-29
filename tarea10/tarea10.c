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

typedef struct SouriauMatrix{
	Matrix_ptr Bn;
	Matrix_ptr qns;
}SouriauMatrix, *SouriauMatrix_ptr;

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
		printf("      ");
		for (int j = 0; j<M->n; j++){
			printf("%f\t",get_matrix_value(M, i,j));
		}
		printf("\n");
	}
}

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
	  	printf("      Introduzca el número de filas de la matriz: ");
	  	scanf("%d",&m);
		const unsigned int M = m;

	  	printf("      Introduzca el número de columnas de la matriz: ");
	  	scanf("%d",&n);
		const unsigned int N = n;
	Matrix_ptr A = matrix_alloc(M,N);
	printf("      Introduzca los elementos de la matriz (%dx%d)\n",m,n);
	double elem;
	for(int i=0; i<M; i++){
    		for(int j=0; j<N; j++){
      			printf("      %d,%d: ", i+1,j+1);
      			scanf("%lf", &elem);
						set_matrix_value(A, i, j, elem);
    		}
  	}

	printf("      %s\n","Matriz ingresada:" );
	print_matrix(A);
	printf("\n" );
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
	a11 = get_matrix_value(M,0,0);
	a12 = get_matrix_value(M,0,1);
	a13 = get_matrix_value(M,0,2);
	a21 = get_matrix_value(M,1,0);
	a22 = get_matrix_value(M,1,1);
	a23 = get_matrix_value(M,1,2);
	a31 = get_matrix_value(M,2,0);
	a32 = get_matrix_value(M,2,1);
	a33 = get_matrix_value(M,2,2);
	det = (a11*a22*a33) + (a12*a23*a31) + (a13*a21*a32) - ((a31*a22*a13)+(a32*a23*a11)+(a33*a21*a12));
	return det;
}

Matrix_ptr adjoint_3x3(Matrix_ptr M){
	double a11,a22,a33,a12,a13,a21,a23,a31,a32;

	int m = M->m;
	int n = M->n;
	double val;
	Matrix_ptr Ad = matrix_alloc(m,n);
	a11 = get_matrix_value(M,0,0);
	a12 = get_matrix_value(M,0,1);
	a13 = get_matrix_value(M,0,2);
	a21 = get_matrix_value(M,1,0);
	a22 = get_matrix_value(M,1,1);
	a23 = get_matrix_value(M,1,2);
	a31 = get_matrix_value(M,2,0);
	a32 = get_matrix_value(M,2,1);
	a33 = get_matrix_value(M,2,2);


	set_matrix_value(Ad, 0,0, (a22*a33)-(a32*a23));
	set_matrix_value(Ad, 0,1, -(a21*a33)+(a31*a23));
	set_matrix_value(Ad, 0,2, (a21*a32)-(a31*a22));
	set_matrix_value(Ad, 1,0, -(a12*a33)+(a32*a13));
	set_matrix_value(Ad, 1,1, (a11*a33)-(a31*a13));
	set_matrix_value(Ad, 1,2, -(a11*a32)+(a12*a31));
	set_matrix_value(Ad, 2,0, (a12*a23)-(a22*a13));
	set_matrix_value(Ad, 2,1, -(a11*a23)+(a21*a13));
	set_matrix_value(Ad, 2,2, (a11*a22)-(a21*a12));
	Ad = matrix_transpose(Ad);
	return Ad;

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

Matrix_ptr matrix_inverse_3x3(Matrix_ptr M ){

	int m = M->m;
	int n = M->n;

	Matrix_ptr Inv = matrix_alloc(m,n);

	Matrix_ptr ad = adjoint_3x3(M);
	double det = det3x3(M);

	if(det == 0 ){
		return -1;
	}

	Inv = scalar_mult(ad, (1/det));

	free(ad);

	return Inv;
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

//Ax = b
Matrix_ptr general_interation_method(Matrix_ptr A, Matrix_ptr B, Matrix_ptr b, Matrix_ptr x_0){

	int n = A->n;
	int m = A->m;
	double error = 1e-5;
	double norm_err;
	int iteration = 0;
	//allox x as initial condition (0 0 0)
	if(x_0 == NULL){
		Matrix_ptr x_0 = matrix_alloc(m, 1);
	}

	Matrix_ptr x;
	//get Identity matrix
	Matrix_ptr B_inv =  matrix_inverse_3x3(B);
	Matrix_ptr I = get_identity_matrix(m,n);
	//Matrix_ptr G = matrix_alloc(m,n);

	//generate iteration matrix G
	Matrix_ptr G = mult_matrix(B_inv, A);
	G = scalar_mult(G, -1);
	G = sum_matrix(I, G);


	do{
		//iteration x =Gx_0 + B_inv*b
		x = mult_matrix(G,x_0);
		Matrix_ptr b_product = mult_matrix(B_inv, b);
		x = sum_matrix(x, b_product);
		//substract x-x_0
		x_0 = scalar_mult(x_0, -1);
		Matrix_ptr err = sum_matrix(x,x_0);

		norm_err = matrix_norm(err);
		x_0 = x;
		//printf("Iteration error: %f, vs expected error: %f\n", norm_err, error );
		free(err);
		free(b_product);
		iteration ++;
	}
	while(norm_err > error);

	printf("Numero total de iteraciones: %d\n",iteration );

	return x;
}

//Ax = b
int general_interation_method_2(Matrix_ptr A, Matrix_ptr B, Matrix_ptr b, Matrix_ptr x_0, double error, double w){

	int n = A->n;
	int m = A->m;
	//double error = 1e-5;
	double norm_err;
	int iteration = 0;
	//allox x as initial condition (0 0 0)
	if(x_0 == NULL){
		Matrix_ptr x_0 = matrix_alloc(m, 1);
	}
	Matrix_ptr x;
	//get Identity matrix
	Matrix_ptr B_inv =  matrix_inverse_3x3(B);
	B_inv = scalar_mult(B_inv, w);
	Matrix_ptr I = get_identity_matrix(m,n);
	//Matrix_ptr G = matrix_alloc(m,n);
	//generate iteration matrix G
	Matrix_ptr G = mult_matrix(B_inv, A);
	G = scalar_mult(G, -1);
	G = sum_matrix(I, G);
	do{
		//iteration x =Gx_0 + B_inv*b
		x = mult_matrix(G,x_0);
		Matrix_ptr b_product = mult_matrix(B_inv, b);
		x = sum_matrix(x, b_product);
		//substract x-x_0
		x_0 = scalar_mult(x_0, -1);
		Matrix_ptr err = sum_matrix(x,x_0);
		norm_err = matrix_norm(err);
		x_0 = x;
		//printf("Iteration error: %f, vs expected error: %f\n", norm_err, error );
		free(err);
		free(b_product);
		iteration ++;
	}
	while(norm_err > error);
	//printf("Numero total de iteraciones: %d\n",iteration );
	return iteration;
}

Matrix_ptr jacobi_iteration_method(Matrix_ptr A, Matrix_ptr b, Matrix_ptr x_0){
	//this method sets the B matrix as Ad
	//1. get ldr decomposition
	LDR_Matrix_ptr ldr=  LDR_decomposition(A);
	//set B as the D matrix;
	Matrix_ptr B = ldr->D;

	Matrix_ptr x= general_interation_method(A,B,b, x_0);
	printf("Solucion de x por Jacobi\n");
	print_matrix(x);
	return x;

}

Matrix_ptr gauss_seidel_iteration_method(Matrix_ptr A, Matrix_ptr b, Matrix_ptr x_0){
	//this method set the B matrix as D+L
	//1. get ldr decomposition
	LDR_Matrix_ptr ldr=  LDR_decomposition(A);
	//add R and L matrices
	Matrix_ptr B = sum_matrix(ldr->D, ldr->L);

	Matrix_ptr x= general_interation_method(A,B,b, x_0);
	printf("Solucion de x por Gauss-seidel\n");
	print_matrix(x);
	return x;

}

Matrix_ptr generate_hilbert_matrix(int m){

	Matrix_ptr H = matrix_alloc(m,m);
	double val;

	for(int i =  0; i < m ;i++){
		for(int j = 0; j < m; j++){
			val = 1.0/(i+1+j+1-1);
			//printf("%f\n", val);
			set_matrix_value(H, i,j, val);
		}
	}

	printf("HILBERT Matrix\n");
	print_matrix(H);
	return H;

}

double get_max(Matrix_ptr a){
  int m = a->m;
  int n = a->n;
  //printf("%s\n","1" );
  double max = fabs(get_matrix_value(a,0,0));
  //printf("%s\n","2" );
  for(int i = 0;i<m;i++){
    for(int j = 0;j<n; j++){
      //printf("%d\n",i );
      double item = fabs(get_matrix_value(a,i,j));
      if (max < item){
        max = item;
      }
    }

  }


  return max;
}

double get_min(Matrix_ptr a){
  int m = a->m;
  int n = a->n;
  double min = fabs(get_matrix_value(a,0,0));

  for(int i = 0;i<m;i++){
    for(int j = 0;j<n; j++){

      double item = fabs(get_matrix_value(a,i,j));
      if (min > item){
        min = item;
      }

    }

  }
  return min;
}

EigenStruct_ptr get_eigen_values(Matrix_ptr A, Matrix_ptr x, double prec ){
  //get matrix size
  int m = A->m;
  int n = A->n;
  //alloc memory for eigen struct
  EigenStruct_ptr eigen = (EigenStruct_ptr)malloc(sizeof(EigenStruct));
  Matrix_ptr y = matrix_alloc(m,1);
  Matrix_ptr x_new = matrix_alloc(m,1);
  double C;
  double err;
  double e = 10e-5;

  int cont = 0;
  do{
    /* code */
    y = mult_matrix(A,x);
    //get max of y
    C = get_max(y);
    x = scalar_mult(x,-1);
    x_new = scalar_mult(y, 1/C);
    Matrix_ptr diff = sum_matrix(x_new, x);
    err = matrix_norm(diff);
    x = x_new;

    printf("Error obtained %f, Error desired %f\n",err,e );
    printf(" C variable %f\n",C );
    cont ++;

  }while(err > e);

  eigen->v = x_new;
  eigen->lambda = C;

  printf("      %s\n","################################" );
  printf("      Max eigen value: %f\n", eigen->lambda);
  printf("      %s\n", "Eigen vector: ");
  print_matrix(eigen->v);
  printf("      Numero de iteraciones: %d\n", cont);
  printf("      %s\n\n","################################" );
  return eigen;
}

double relaxation_method(Matrix_ptr A, Matrix_ptr b, Matrix_ptr x_0, double error){
	double limit_max = 2.0;
	double limit_min = 0.1;
	double w = 0.1;
	double w_min = 0.1;
	int ext = 0;
	int iter_min = 0;
	int iter = 0;
	LDR_Matrix_ptr ldr = LDR_decomposition(A);
	Matrix_ptr B = ldr->D;
	do {
		while (w <= limit_max) {
				iter = general_interation_method_2(A, B, b, x_0, error, w);
				printf("\n Iter %d ----- w %lf ", iter, w);
				if (iter_min == 0 || iter < iter_min) {
					iter_min = iter;
					w_min = w;
				}
				w = w + 0.1;
		}
		if ((w_min != limit_min && w_min != limit_max) || w_min == limit_min){
			ext = 1;
		} else {
			limit_min = limit_max;
			limit_max = limit_min + 2.0;
			w = limit_min;
			w_min = w;
			iter_min = 0;
			iter = 0;
		}
	} while (ext != 1);
	printf("\n%s\n", "      Utilizando metodo de relajacion" );
	printf("      w optimo: %lf\n", w_min );
	printf("      iteraciones: %d\n", iter_min );
	//printf("\n\n   W minimo: %f", w_min);
	//printf("\n   Iteraciones minimas: %d \n\n", iter_min);
	return w_min;
}

double over_relaxation_method(Matrix_ptr A, Matrix_ptr b, Matrix_ptr x_0, double error){
	double limit_max = 2.0;
	double limit_min = 0.1;
	double w = 0.1;
	double w_min = 0.1;
	int ext = 0;
	int iter_min = 0;
	int iter = 0;
	LDR_Matrix_ptr ldr = LDR_decomposition(A);
	Matrix_ptr B = sum_matrix(ldr->D, ldr->L);
	do {
		while (w <= limit_max) {
				iter = general_interation_method_2(A, B, b, x_0, error, w);
				printf("\n\n Iter %d ----- w %lf \n", iter, w);
				if (iter_min == 0 || iter < iter_min) {
					iter_min = iter;
					w_min = w;
				}
				w = w + 0.1;
		}
		if ((w_min != limit_min && w_min != limit_max) || w_min == limit_min){
			ext = 1;
		} else {
			limit_min = limit_max;
			limit_max = limit_min + 2.0;
			w = limit_min;
			w_min = w;
			iter_min = 0;
			iter = 0;
		}
	} while (ext != 1);
	printf("\n%s\n", "      Utilizando metodo de sobre-relajacion" );
	printf("      w optimo: %lf\n", w_min );
	printf("      iteraciones: %d\n", iter_min );
	//printf("\n\n   W minimo: %f", w_min);
	//printf("\n   Iteraciones minimas: %d", iter_min);
	return w_min;
}

double get_trace(Matrix_ptr A){
	int m = A->m;
	int n = A->n;
	int i = 0;
	double accum = 0;
	for(int i= 0; i<m;i++){
		for(int j = 0 ;j<n; j++){/* code */
			if (i ==j){
				accum = accum + get_matrix_value(A,i,j);
			}
		}
	}
	return accum;
}

Matrix_ptr pow_matrix(Matrix_ptr A, int power){
	int m = A->m;
	int n = A->n;
	int i = 0;
	Matrix_ptr R = matrix_alloc(m,n);
	R = get_identity_matrix(m,n);
	if (power != 0){
		for (i=1; i<=power; i++ ){
			R = mult_matrix(R, A);
		}
	}
	return R;
}

Matrix_ptr matrix_zeros(int m, int n){
	Matrix_ptr A = matrix_alloc(m,n);
	for(int i= 0; i<m;i++){
		for(int j = 0 ;j<n; j++){
			set_matrix_value(A, i, j, 0);
		}
	}
	return A;
}


double leverrier_aux(int index_p, int index_s, int k, Matrix_ptr pol, Matrix_ptr A){
	double p_value, s_value;
	if (index_p==k || index_s==0){
		return 0;
	} else {
		p_value = get_matrix_value(pol, 0, index_p);
		s_value = get_trace(pow_matrix(A, index_s));
		return p_value*s_value + leverrier_aux(index_p+1, index_s-1, k, pol, A);
	}
}

Matrix_ptr leverrier_algorithm(Matrix_ptr A){
	int m = A->m;
	int n = A->n;
	int k = 1;
	double pk = 0;
	//allocate memory for the ldr structure
	Matrix_ptr pol = matrix_alloc(1,n+1);
	set_matrix_value(pol, 0, 0, 1);
	//allocate memory for the B and I
	Matrix_ptr Ak = matrix_alloc(m,n);
	for (k=1; k<=n; k++){
		Ak = pow_matrix(A, k);
		pk = (-1.0/k)*(get_trace(Ak) + leverrier_aux(1, k-1, k, pol, A));
		set_matrix_value(pol, 0, k, pk);
	}
	return pol;
}

SouriauMatrix_ptr souriau_method(Matrix_ptr A){
	int m = A->m;
	int n = A->n;
	int k = 1;
	double qn = 0;
	//allocate memory for the ldr structure
	Matrix_ptr pol = matrix_alloc(1,n+1);
	set_matrix_value(pol, 0, 0, 1);
	//allocate structure
	SouriauMatrix_ptr souriau = (SouriauMatrix_ptr)malloc(sizeof(SouriauMatrix));
	souriau->Bn = matrix_alloc(m,n);
	souriau->qns = matrix_alloc(1,n+1);
	//allocate memory for the B and I
	Matrix_ptr An = matrix_alloc(m,n);
	Matrix_ptr Bn = matrix_alloc(m,n);
	Matrix_ptr Bn_1 = matrix_alloc(m,n);
	Matrix_ptr I = get_identity_matrix(m,n);

	Bn = I;
	for (k=1; k<=n; k++){
		An = mult_matrix(A, Bn);
		qn = (-1.0/k)*get_trace(An);
		Bn_1 = Bn;
		Bn = sum_matrix(An, scalar_mult(I, qn));
		set_matrix_value(pol, 0, k, qn);
	}
	souriau->Bn = Bn_1;
	souriau->qns = pol;
	return souriau;
}

Matrix_ptr matrix_inverse_leverrier(Matrix_ptr pol, Matrix_ptr A){
	int m = A->m;
	int n = A->n;
	Matrix_ptr A_inv = matrix_alloc(m, n);
	A_inv = matrix_zeros(m, n);
	for (int i=1; i<=n; i++){
		A_inv = sum_matrix(A_inv, scalar_mult(pow_matrix(A, n-i), get_matrix_value(pol, 0, i-1)));
	}
	A_inv = scalar_mult(A_inv, (-1.0/get_matrix_value(pol, 0, n)));
	return A_inv;
}

Matrix_ptr matrix_inverse_souriau(Matrix_ptr Bn_1, double qn){
	Matrix_ptr A_inv = scalar_mult(Bn_1, (-1.0/qn));
	return A_inv;
}

void print_polynomial(Matrix_ptr pol, char* var) {
	for (int i = 0; i<pol->m; i++){
		printf("      ");
		for (int j = 0; j<pol->n; j++){
			if ( pol->n - j -1 > 0){
				printf("%.2f%s^%d \t",get_matrix_value(pol, i,j), var, pol->n -j-1 );
			} else {
				printf("%.2f\t",get_matrix_value(pol, i,j));
			}

		}
		printf("\n");
	}
}

int main(){
	int opcion;
	double error, w1, w2;
	printf("\n    ---------------------------------------------" );
	printf("\n   |        Metodos Computaionales (MCEI)        |" );
	printf("\n   |                  Tarea 10                   |" );
	printf("\n    ---------------------------------------------" );
	printf("\n   |                                             |" );
	printf("\n   |  1. Polinomio Caracteristico e Inversa      |" );
	printf("\n   |  2. Metodo QR                               |" );
	printf("\n   |  3. Salir                                   |" );
	printf("\n   |                                             |" );
	printf("\n    ---------------------------------------------" );
	printf("\n\n      Seleccione una opcion (1 - 3): " );
	scanf("%d", &opcion );

	switch (opcion) {
		case 1:
			printf("\n\n      Se ha seleccionado la opcion (1)\n" );
			printf("%s\n", "      Ingrese la matriz A\n");
			Matrix_ptr A = user_request_matrix();
			printf("\n%s", "      ---------------------------------------" );
			printf("\n%s\n", "      Utilizando Algoritmo de Leverrier:" );
			printf("%s\n", "      Polinomio Caracteristico:" );
			Matrix_ptr pol1 = leverrier_algorithm(A);
			print_polynomial(pol1, "x");
			Matrix_ptr A_inv1 = matrix_inverse_leverrier(pol1, A);
			printf("%s\n", "      Inversa:" );
			print_matrix(A_inv1);
			printf("\n\n%s", "      ---------------------------------------" );
			printf("\n%s\n", "      Utilizando Metodo de Souriau" );
			printf("%s\n", "      Polinomio Caracteristico:" );
			SouriauMatrix_ptr sou = souriau_method(A);
			print_polynomial(sou->qns, "x");
			Matrix_ptr A_inv2 = matrix_inverse_souriau(sou->Bn, get_matrix_value(sou->qns, 0, sou->qns->n-1));
			printf("%s\n", "      Inversa:" );
			print_matrix(A_inv2);
			break;
		case 2:
			printf("\n\n Se ha seleccionado la opcion (2)" );
			/*
			TODO:
			QR Method*/
			break;
		case 3:
			printf("\n      Finalizando programa\n\n" );
			break;
		default:
			printf("\n      Por favor ingrese una opcion valida \n" );
			main();
	}
	return 0;


}
