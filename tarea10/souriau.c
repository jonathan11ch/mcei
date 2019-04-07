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

Matrix_ptr matrix_inverse_souriau(Matrix_ptr Bn_1, double qn){
	Matrix_ptr A_inv = scalar_mult(Bn_1, (-1.0/qn));
	return A_inv;
}
