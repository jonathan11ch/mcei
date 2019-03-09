#include<stdio.h>
#include<math.h>
#include<gsl/gsl_blas.h>


void print_matrix(const gsl_matrix *m) {
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			printf("%.20f\t", gsl_matrix_get(m, i, j));
		}
		printf("\n");
	}
}

void print_vector(const gsl_vector *v){
  for (int j = 0; j < v->size; j++) {
    printf("%.20f\t", gsl_vector_get(v, j));
  }

}

double det_matrix_R(gsl_matrix *r){
	double det = 1;
	double val;
	int M =  r->size1;
	int N = r->size2;

	//calculate multiplication of main diagonal
	for(int i = 0; i < M ;i++){
			for(int j = 0; j < N ;j++){
					if(i == j){
						val = gsl_matrix_get(r,i,j);
						det = det * val;
					}
			}
	}

	return det;

}

void calculate_gauss(gsl_matrix *A){
    printf("%s\n","Calculate gauss of Matrix" );
    int M = A->size1;
    int N = A->size2;
    int index;
    double det, val;
    double pivot,c,max;
    gsl_vector *col =gsl_vector_alloc(M);
    gsl_vector *row =gsl_vector_alloc(N);
    gsl_vector *row_1 =gsl_vector_alloc(N);
    gsl_matrix *L = gsl_matrix_alloc(M,N);
    gsl_matrix *R = gsl_matrix_alloc(M,N);
    gsl_matrix *P = gsl_matrix_alloc(M,N);
    gsl_matrix *P_NEW = gsl_matrix_alloc(M,N);
    gsl_matrix *PP = gsl_matrix_alloc(M,N);

    gsl_matrix_set_zero(L);
    gsl_matrix_set_zero(R);
    //copy A to R
    gsl_matrix_memcpy(R,A);

    printf("%s\n","set diagonals" );
    //set L main diagonal with ones
    for(int i = 0; i < M ;i++){
        for(int j = 0; j < N ;j++){
            if(i == j){
              //gsl_matrix_set(L,i,j,1);
              gsl_matrix_set(P,i,j,1);
              //gsl_matrix_set(P_NEW,i,j,1);
            }
        }
    }

    //start gauss
    printf("%s\n","Start Gauss calculation");
    print_matrix(P);
    //go column by column
    for(int j = 0; j < N ;j++){
        printf("%s\n","get col" );
        //check maximun value of the row
        gsl_matrix_get_col(col,A,j);
        max = gsl_vector_get(col,j);
        index = j;
        printf("%s\n","get max index");
        for (int x = j; x < M; x++){
          if(gsl_vector_get(col,x) > max){
            index = x;
            max = gsl_vector_get(col,x);
          }
        }
        //get index
        //index = gsl_vector_max_index(col);
        printf("Max value located in row: %d\n", index );
        if (index != j){
          //must perform permutation
          gsl_matrix_swap_rows(A,j,index);
          gsl_matrix_swap_rows(P,j,index);
          //gsl_matrix_swap_rows(L,j,index);

          //gsl_blas_dgemm(CblasNoTranblasNoTrans,1., P_NEW,P,0.,PP);
          //gsl_matrix_memcpy(P, PP);
          printf("%s\n","SWAPPPPPPPPP[PPP]" );
          print_matrix(A);
        }

        pivot = gsl_matrix_get(A,j,j);
        printf("Pivot: %f\n",pivot );
        if(pivot == 0){
          //null row, must change
          printf("%s\n", "null pivot, move to the next column");
        }
        else{
          //apply gauss
            for(int i = j+1; i < M ;i++){

                c = gsl_matrix_get(A,i,j);
                printf("Element i, j %d,%d: %f\n",i,j,c );
                c = c/pivot;
                c = c*-1;
                //set value in L
                printf("indices para settear L %d,%d\n",i,j );
                gsl_matrix_set(L,i,j,-c);
                printf("Constant C %f\n",c );
                print_matrix(L);
                printf("%s\n","--------------------------" );
                gsl_matrix_get_row(row, A, j);
                gsl_vector_scale(row,c);
                gsl_matrix_get_row(row_1,A,i);
                gsl_vector_add(row_1,row);
                //printf("%s\n","Vector" );
                //print_vector(row);
                //printf("%s\n","Vector1" );
                //print_vector(row_1);
                gsl_matrix_set_row(A,i,row_1);
                printf("%s\n", "#####################################");
                print_matrix(A);
                printf("%s\n", "#####################################");
            }

        }


    }
    for(int i = 0; i < M ;i++){
        for(int j = 0; j < N ;j++){
            if(i == j){
              gsl_matrix_set(L,i,j,1);
              //gsl_matrix_set(P,i,j,1);
              //gsl_matrix_set(P_NEW,i,j,1);
            }
        }
    }


    printf("-----------------Matrix L ------------------\n");
    print_matrix(L);
    printf("-----------------Matrix R ------------------\n");
    print_matrix(A);
    printf("-----------------Matrix P ------------------\n");
    print_matrix(P);
    printf("-----------------Matrix A ------------------\n");
    print_matrix(R);

    //calculate det of A
		det = det_matrix_R(A);
		printf("Calculated determinant %f\n", det);


    det = 1;
    for(int i = 0; i < M ;i++){
        for(int j = 0; j < N ;j++){
            if(i == j){
              val = gsl_matrix_get(A,i,j);
              det = det * val;
            }

        }
    }
    printf("Determinant of Matrix A: %f\n", det);


}





int main(){

  gsl_vector *test = gsl_vector_alloc(3);
  int PIVOT;
	int m,n;            // Número de filas y columnas de la matriz de usuario
  	printf("Introduzca el número de filas de la matriz: ");
  	scanf("%d",&m);
	const unsigned int M = m;

  	printf("Introduzca el número de columnas de la matriz: ");
  	scanf("%d",&n);
	const unsigned int N = n;

	//const double rcond = 1E-15;
	gsl_matrix *A = gsl_matrix_alloc(M,N);
	gsl_matrix *A_new = gsl_matrix_alloc(M,N);
	printf("Introduzca los elementos de la matriz A(%dx%d)\n",m,n);
	double elem;
	for(int i=0; i<M; i++){
    		for(int j=0; j<N; j++){
      			printf("a%d,%d: ", i+1,j+1);
      			scanf("%lf", &elem);
			gsl_matrix_set(A, i, j, elem);
    		}
  	}

  gsl_vector_set(test,0,1);
  gsl_vector_set(test,1,1);
  gsl_vector_set(test,2,1);
  PIVOT = gsl_vector_max_index(test);
  printf("Pivot test .. %d\n", PIVOT);
  calculate_gauss(A);

  printf("--------------------- Result -------------------\n");
	print_matrix(A);

	return 0;





}
