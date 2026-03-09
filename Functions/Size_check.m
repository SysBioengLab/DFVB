function sizes = Size_check(model,iterations,correct_size)

if nargin < 3
    correct_size = size(NullS(model,'ep_vTol',1e-4,'vTol',1e-8));
end

wrong = 0;
right = 0;
sizes = zeros(iterations,2);

for i = 1:iterations
    base   = NullS(model,'ep_vTol',1e-4,'vTol',1e-8);

    sizes(i,:) = size(base);

    if isequal(size(base),correct_size)
        right = right + 1;
    else
        wrong = wrong + 1;
    end
end

fprintf('\nResults of %d iterations \n',iterations)
disp('-------------------------')
fprintf('Same size base: %d \n',right)
fprintf('Different size base: %d \n',wrong)

if wrong == 0
    fprintf('\nCheck Rank of basis\n')
    disp('-------------------------')
    fprintf('Rank of basis: %d \n',rank(base))
    fprintf('Col number of basis: %d \n',size(base,2))
end