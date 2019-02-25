#include <stdio.h>
#include <gsl/gsl_blas.h>


typedef unsigned int bool;
#define false		0
#define true 		1

//gcc -Wall -02 mpenrose.c -o penrose -lgsl -lgslcblas -lm

// Función que muestra en consola una matriz (estructura propia de libreria gsl)
void mostrar_matriz(const gsl_matrix *m) {
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			printf("%f\t", gsl_matrix_get(m, i, j));
		}
		printf("\n");
	}
}

void mostrar_vector(const gsl_vector *v){
	printf("vector:\n");
	for(int i = 0 ;i<v->size ; i++){
			printf("%f\t", gsl_vector_get(v, i));
	}
	printf("\n");


}

gsl_matrix* penrose(gsl_matrix* A){
	unsigned int n = A->size2;
	unsigned int m = A->size1;


	gsl_matrix *V, *Sigma_inv, *U, *A_pinv;
	gsl_vector *_tmp_vec, *s;
	gsl_matrix *_tmp_mat;

	bool transpose = false;
	


	//validate size of matrix
	if(n > m){
		transpose = true;
		_tmp_mat = gsl_matrix_alloc(n, m);
		gsl_matrix_transpose_memcpy(_tmp_mat, A);
		
		n = A->size1;
		m = A->size2;
		printf("indices m: %d, n: %d \n",m,n );
		A = _tmp_mat;
		mostrar_matriz(A);
		

	}


	V = gsl_matrix_alloc(n, n);
	U = gsl_matrix_alloc(m, m);
	s = gsl_vector_alloc(n);

	gsl_matrix_set_zero(U);
	_tmp_vec = gsl_vector_alloc(n);

	gsl_linalg_SV_decomp(A, V, s, _tmp_vec);

	printf("\n Matriz V\n");
	mostrar_matriz(V);
	printf("\n Matriz S\n");
	mostrar_vector(s);
	printf("\n Matriz U\n");
	mostrar_matriz(A);

	//llena la matriz U 
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			gsl_matrix_set(U, i, j, gsl_matrix_get(A, i, j));
		}
	}
 

	//create sigma inverse function
	Sigma_inv = gsl_matrix_alloc(n,m);
	gsl_matrix_set_zero(Sigma_inv);

	for(int i = 0;i<s->size;i++){
		if(gsl_vector_get(s,i)>0){
			gsl_matrix_set(Sigma_inv,i,i,1/ gsl_vector_get(s,i));
		}
		else{
			gsl_matrix_set(Sigma_inv,i,i,0.0);
		}

	}
	printf("Matriz Sigma\n");
	mostrar_matriz(Sigma_inv);
	printf("\n Matriz U\n");
	mostrar_matriz(U);

	//A_pinv= V * Sigma_inv * U_transpose
 
	return	A;

}





int main()
{
	int m,n;            // Número de filas y columnas de la matriz de usuario
	printf("Introduzca el número de filas de la matriz: ");
	scanf("%d",&m);
	const unsigned int M = m;

  	printf("Introduzca el número de columnas de la matriz: ");
  	scanf("%d",&n);
	const unsigned int N = n;

	const double rcond = 1E-15;
	gsl_matrix *A = gsl_matrix_alloc(M,N);    // Matriz de usuario
	gsl_matrix *A_pinv;											  // Matriz de pseudo-inversa

	printf("Introduzca los elementos de la matriz A(%dx%d)\n",m,n);
	double elem;
	
	for(int i=0; i<M; i++){
    	
    	for(int j=0; j<N; j++){
	      	
	      	printf("a%d,%d: ", i+1,j+1);
	      	scanf("%lf", &elem);
			gsl_matrix_set(A, i, j, elem);
    	
    	}
  	
  	}

	printf("\n     Matriz de usuario\n");
	mostrar_matriz(A);
	printf("\nCalculando Inversa de Moore-Penrose\n");
	A_pinv = penrose(A);
	//mostrar_matriz(A_pinv);

	printf("Fin.\n");


	return 0;
}