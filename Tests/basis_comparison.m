load yeast9-GEM.mat
tol = 1e-9;
A = model.S;
% Dense
N_dense = null(full(A));
N_dense(abs(N_dense) < tol) = 0; 
nnz_dense = nnz(N_dense);
fill_dense = (nnz_dense / numel(N_dense)) * 100;

any
% Sparse
[N_sparse, rankS] = sparseNull(A, tol);
nnz_sparse = nnz(N_sparse);
fill_sparse = (nnz_sparse / numel(N_sparse)) * 100;

%% Plot
figure('Name', 'Sparsity Pattern: A, Nativo vs LUQ', 'Position', [100, 100, 1200, 450]);

subplot(1, 3, 1);
spy(A, 'k');
title('Stoichiometric Matrix')
axis square;
set(gca, 'FontSize', 11);

subplot(1, 3, 2);
spy(N_dense, 'r');
title('SVD basis')
axis square; 
set(gca, 'FontSize', 11);

subplot(1, 3, 3);
spy(N_sparse, 'b');
title('RRLU basis')
axis square;
set(gca, 'FontSize', 11);