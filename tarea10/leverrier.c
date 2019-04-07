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
