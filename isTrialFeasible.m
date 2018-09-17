function feasible = isTrialFeasible(X,Aineq,bineq,Aeq,beq,lb,ub,tol)
%isTrialFeasible Checks if X is feasible w.r.t. linear constraints. 
% 	

%   Copyright 2003-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2012/08/21 00:24:22 $

feasibleIneq = true;
feasibleEq = true;

% Inequality constraints
if ~isempty(Aineq)
    feasibleIneq = max(Aineq*X-bineq) <= tol;
end
% Upper bounds
argub = ~eq(ub,inf);
if feasibleIneq && any(argub) && ~isempty(argub)
    feasibleIneq = all(X(argub) <= ub(argub));
end
% Lower bounds
arglb = ~eq(lb,-inf);
if feasibleIneq && any(arglb) && ~isempty(arglb)
    feasibleIneq = all(X(arglb) >= lb(arglb));
end
% Equality constraints
if feasibleIneq && ~isempty(Aeq) 
    constrViolation = (Aeq*X-beq);
    feasibleEq   = all(abs(constrViolation(~isinf(constrViolation))) <= tol);
end
feasible = feasibleIneq && feasibleEq;
