function [feasible,XOUT] = getFeasiblePoints(X,Aineq,bineq,Aeq,beq,lb,ub,tol)
%GETFEASIBLEPOINTS returns all feasible points in point matrix 'X' w.r.t. 
%   linear/bound constraints within the tolerance 'tol'.
%
%   Output arguments 'XOUT' is a matrix of all the feasible points and 'feasible'
%   is a logical array of length size(X,2) indicating the point is feasible (TRUE),
%   or infeasible (FALSE).
%
%   Example:
%     If there are 4 points in 2 dimension space, [2;-4],[1;5],[9;0]
%     and [-2;1] then
%
%     X  =   [2  1 9 -2
%            -4  5 0  1 ]
%
%     and if Aineq = diag([-2,2]), bineq = zeros(2,1);
%     tol = 1e-6; we obtain [Xout,feasible] = getFeasiblePoints(X,Aineq,bineq,[],[],[],[],tol)
%
%     Xout = [-2;1]
%

%   Copyright 2003-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2012/08/21 00:24:13 $


feasible = true(size(X,2),1);
XOUT=[];
for i = 1:size(X,2)
    feasible(i) = isTrialFeasible(X(:,i),Aineq,bineq,Aeq,beq,lb,ub,tol);
    if feasible(i)
        XOUT(:,end+1) = X(:,i);
    end
end
