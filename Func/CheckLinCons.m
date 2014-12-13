function bool = CheckLinCons(x,LHS,RHS)
tol = 1e-9;

bool = isempty(find(LHS*x-RHS>tol,1));