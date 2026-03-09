function [intN, lb_int, ub_int, nullS, lb_ex, ub_ex, c, zero_rows_int, zero_rows_ex, rxn_names]  = forced_idx(model)

l_int = size(model.intN,1);
l_ex  = size(model.nullS,1);

% Extract boundaries of internal and external fluxes
lb_int  = model.rhs(l_ex+l_int+1:end);
ub_int  = model.rhs(l_ex+1:l_int+l_ex);

ub_ex   = model.ub(1:l_ex);
lb_ex   = model.lb(1:l_ex);

c = model.obj(1:l_ex);

% Clean Internal matrix and boundaries
zero_rows_int  = unique(find(all(model.intN == 0, 2))');
non_zero_rows  = setdiff(1:size(model.intN,1), zero_rows_int);

intN     = model.intN(non_zero_rows,:);
lb_int   = lb_int(non_zero_rows);
ub_int   = ub_int(non_zero_rows);


% Clean external matrix and boundaries
zero_rows_ex      = find(all(model.nullS == 0, 2))';
non_zero_rows_ex  = setdiff(1:size(model.nullS,1), zero_rows_ex);

nullS = model.nullS(non_zero_rows_ex,:);
lb_ex = lb_ex(non_zero_rows_ex);
ub_ex = ub_ex(non_zero_rows_ex);
c     = c(non_zero_rows_ex);

rxn_names = model.rxns(non_zero_rows_ex);
end