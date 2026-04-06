function degree = Sparse_degree(A)
% --------------------------------------------------------------------------
% Computes the sparsity degree of a matrix A.
% Inputs:
%   A - Input matrix.
% Outputs:
%   degree - Sparsity degree of the matrix, defined as the ratio of zero elements to total elements.
% --------------------------------------------------------------------------

num_zeros = nnz(A(:) == 0);
num_total = numel(A);
degree    = num_zeros / num_total;