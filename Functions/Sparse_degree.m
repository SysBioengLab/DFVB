function degree = Sparse_degree(A)

num_zeros = nnz(A(:) == 0);
num_total = numel(A);
degree    = num_zeros / num_total;