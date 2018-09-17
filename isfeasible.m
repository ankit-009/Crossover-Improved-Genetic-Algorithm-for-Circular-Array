function feasible = isfeasible(X,A,L,U,tol,IndEqcstr)
%isfeasible Checks to see if current iterate X is feasible or not. 
% 	
% 	X: is the current point which for which we check feasibility.
% 	A,L,U: Defines the feasible region in case of linear/bound constraints.
% 	L<=X<=U.
% 	
% 	TOL: Tolerance used to determine feasibility.
% 	
% 	IndEqcstr: Logical indices of equality constraints. A(IndEqcstr), LB(IndEqcstr)
% 	UB(IndEqcstr) represents equality constraints.

%   Copyright 2003-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2012/08/21 00:24:25 $


%If there are no equality constraints
if ~any(IndEqcstr)
    feasible = all(A*X-L>=0 & A*X-U<=0);
    return;
end

if ~isempty(IndEqcstr)
    feasIneq = all(A(~IndEqcstr,:)*X-L(~IndEqcstr,:)>=0) && all(A(~IndEqcstr,:)*X-U(~IndEqcstr,:)<=0);
    lower    = (A(IndEqcstr,:)*X-L(IndEqcstr,:));
    upper    = (A(IndEqcstr,:)*X-U(IndEqcstr,:));
    feasEq   = all(abs(lower(~isinf(lower)))<=tol ) && all(abs(upper(~isinf(upper)))<=tol);
    feasible = feasIneq && feasEq;
else
    feasible = all(A*X-L>=0 & A*X-U<=0);
end