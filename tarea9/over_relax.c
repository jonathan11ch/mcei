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
	return w_min;
}
