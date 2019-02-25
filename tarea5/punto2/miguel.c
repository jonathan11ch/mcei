/**
 * Programa que determina la pseudo-inversa de Moore-Penrose de una matriz rectangular
 *
 * Comando de compilación:
 *
 *     gcc -Wall -02 moore_penrose_inv.c -o moore_penrose -lgsl -lgslcblas -lm
 *
 * Dependencias que deben instalarse:
 * - libgsl (GNU Scientific Library)
 * - libblas (Basic Linear Algebra Subprograms)
**/

#include <stdio.h>
#include <gsl/gsl_blas.h>

#define max(a,b)		((a) > (b) ? (a) : (b))
#define min(a,b)		((a) < (b) ? (a) : (b))

typedef unsigned int bool;
#define false		0
#define true 		1

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
/**
 * Cómputo de la pseudo-inversa de Moore-Penrose.
 *
 * Si se tiene la descomposición de valores singulares (SVD) de A = USVᵀ, entonces
 * la pseudo-inversa se define como VS⁺Uᵀ. Elementos más pequeños que "rcond"
 * veces el valor singular más grande se consideran cero (0).
 *
 * Parámetros de entrada:
 * A:	Matriz de entrada.
 * rcond:	Un número real que especifica el umbral de valor singular para inclusión.
 * Valor por defecto: 1E-15.
 *
 * Variables de salida:
 * A_pinv:	Matriz que contiene la pseudo-inversa resultante.
**/
gsl_matrix* moore_penrose_pinv(gsl_matrix *A, const double rcond) {
	gsl_matrix *V, *Sigma_mas, *U, *A_pinv;
	gsl_matrix *_tmp_mat = NULL;
	gsl_vector *_tmp_vec, *u;
	double x, umbral;
	int i, j;
	unsigned int m = A->size1;
	unsigned int n = A->size2;
	bool transpuesta = false;

	if (n > m) {
		/* La librería libgsl SVD sólo puede manejar el caso en que n <= m, por lo que
		la matriz A se transpone y almacena en una matriz temporal en caso contrario*/
		transpuesta = true;
		_tmp_mat = gsl_matrix_alloc(n, m);
		/** Función: int gsl_matrix_transpose_memcpy (gsl_matrix * dest, const gsl_matrix * src)
		/* Esta función hace de la matriz "dest" la transpuesta de la matriz "src"
		/* al copiar sus elementos a "dest".
		*/
		gsl_matrix_transpose_memcpy(_tmp_mat, A);
		A = _tmp_mat;
		i = n;
		n = m;
		m = i;
	}

	// Realizando la descomposición por valores singulares (SVD)
	V = gsl_matrix_alloc(n, n);
	u = gsl_vector_alloc(n);
	_tmp_vec = gsl_vector_alloc(n);
	/* Función: int gsl_linalg_SV_decomp (gsl_matrix * A, gsl_matrix * V, gsl_vector * S, gsl_vector * work)
	/* Esta función factoriza la matriz A de dimensiones MxN en una descomposición
	/* de valores singulares de la forma A = U S Vᵀ para M >= N. La matriz A es reemplazada
	/* a la salida por U. Los elementos de la diagonal de la matriz de valores singulares S
	/* son almacenados en el vector S. Los valores singulares son positivos y forman una
	/* secuencia decreciente de s_1 a s_n. La matriz V debe contener los elementos de V
	/* en forma no transpuesta. Esta función usa el algoritmo Golub-Reinsch.
	*/
	gsl_linalg_SV_decomp(A, V, u, _tmp_vec);
	printf("\n Matriz V\n");
	mostrar_matriz(V);
	printf("\n Matriz S\n");
	mostrar_vector(u);
	printf("\n Matriz U\n");
	mostrar_matriz(A);
	gsl_vector_free(_tmp_vec);

	// Cálculo de S⁺
	Sigma_mas = gsl_matrix_alloc(n, m);
	gsl_matrix_set_zero(Sigma_mas);
	/* Función: double gsl_vector_max (const gsl_vector * v)
	/* Esta función devuelve el valor máximo del vector v.
	*/
	umbral = rcond * gsl_vector_max(u);

	for (i = 0; i < n; ++i) {
		if (gsl_vector_get(u, i) > umbral) {
			x = 1. / gsl_vector_get(u, i);
		}
		else {
			x = 0.;
		}
		gsl_matrix_set(Sigma_mas, i, i, x);
	}
	printf("············sigma··············\n");
	mostrar_matriz(Sigma_mas);
	/* La librería libgsl SVD provee la matriz U sin las filas o columnas de ceros,
	/* por lo que hay que añadírselas
	*/
	U = gsl_matrix_alloc(m, m);
	gsl_matrix_set_zero(U);

	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			gsl_matrix_set(U, i, j, gsl_matrix_get(A, i, j));
		}
	}

	if (_tmp_mat != NULL) {
		gsl_matrix_free(_tmp_mat);
	}

	_tmp_mat = gsl_matrix_alloc(n, m);
	// Se realizan dos productos punto para obtener la matriz pseudo-inversa
	/* Función: int gsl_blas_dgemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, double alpha, const gsl_matrix * A, const gsl_matrix * B, double beta, gsl_matrix * C)
	/* Esta función calcula el producto matricial y la suma C = \alpha op(A) op(B) + \beta C
	/* donde op(A) = A, A^T, A^H para TransA = CblasNoTrans, CblasTrans, CblasConjTrans
	/* y de forma similar para el parámetro TransB.
	*/
	// Cálculo de V*S⁺ y almacenamiento en matriz temporal (_tmp_mat)
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., V, Sigma_mas, 0., _tmp_mat);

	// La pseudo-inversa se calcula dependiendo de si la matriz fue transpuesta o no
	if (transpuesta) {
		A_pinv = gsl_matrix_alloc(m, n);
		// Cálculo de U * (V * S⁺)ᵀ para matriz A transpuesta
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., U, _tmp_mat, 0., A_pinv);
	}
	else {
		A_pinv = gsl_matrix_alloc(n, m);
		// Cálculo de V * S⁺ * Uᵀ para matriz A no transpuesta
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., _tmp_mat, U, 0., A_pinv);
	}
	printf("sigma\n");
	mostrar_matriz(Sigma_mas);

	printf("U MATRIX\n");
	mostrar_matriz(U);
	printf("V MATRIX\n");
	mostrar_matriz(V);

	gsl_matrix_free(_tmp_mat);
	gsl_matrix_free(U);
	gsl_matrix_free(Sigma_mas);
	gsl_vector_free(u);
	gsl_matrix_free(V);

	return A_pinv;
}

int main() {
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
	A_pinv = moore_penrose_pinv(A, rcond);
	printf("\n     Pseudo-inversa de A:\n");
	mostrar_matriz(A_pinv);

	// Liberando las matrices creadas
	gsl_matrix_free(A);
	gsl_matrix_free(A_pinv);

	return 1;
}

