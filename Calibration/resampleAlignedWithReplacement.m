function [data_b, idx] = resampleAlignedWithReplacement(data_in, includeFirst)
% Muestrea con reemplazo EXACTAMENTE las mismas filas para:
%   - data_in.conc (todas las columnas)
%   - data_in.t_exp
% Parámetros:
%   data_in      : struct con campos .conc [n x m], .t_exp [n x 1]
%   includeFirst : si true, asegura que el índice 1 esté incluido.
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