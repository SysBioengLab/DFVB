function [data_b, idx] = resampleAlignedWithReplacement(data_in, includeFirst)
%---------------------------------------------------------------------------------
% Performs bootstrap resampling with replacement of the same rows in all fields:
%   - data_in.conc   (all columns)
%   - data_in.t_exp
%
% Inputs:
%   data_in      : struct with fields .conc [n x m] and .t_exp [n x 1]
%   includeFirst : if true, guarantees that index 1 appears in the sample
% Outputs:
%   data_b : struct with the same fields as data_in but with resampled rows
%   idx    : indices of the selected rows (with replacement)
%---------------------------------------------------------------------------------
    conc = data_in.conc;
    t    = data_in.t_exp(:);
    n    = size(conc,1);
    % Muestreo con reemplazo de n índices
    idx = randi(n, n, 1);
    if includeFirst
        % Forzar la presencia del primer tiempo
        idx(1) = 1;
        % Mezclar para no fijar la posición del tiempo 1
        idx = idx(randperm(n));
    end
    % Aplicar índices a todas las variables (alineado)
    data_b.conc  = conc(idx, :);
    data_b.t_exp = t(idx);
end