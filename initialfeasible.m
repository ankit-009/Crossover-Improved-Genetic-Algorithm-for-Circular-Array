function [XIN,success] = initialfeasible(XIN,n,Aineq,bineq,Aeq,beq,lb,ub)
%INITIALFEASIBLE Finds a feasible point subject to linear constraints.
% 	This function sets up an LP in order to find initial point feasible w.r.t.
% 	the bounds and linear constraints.
%
% 	XIN: Starting guess (may be modified to satisfy the box constraints)
%
% 	SUCCESS: Indicates success or failure in finding a feasible point.

%   Copyright 2003-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2012/08/21 00:24:18 $

% Use linear programming to get initial feasible point
% that is closest to XIN

ineq_ind  = ~isinf(bineq);
nineqcstr = nnz(ineq_ind); % Remove infinite rows.
eq_ind    = ~isinf(beq);
neqcstr   = nnz(eq_ind);

f = [zeros(n,1); 1];
% Combine Aineq and bounds
Aineq2 = [ eye(n)        , -1*ones(n,1);
          -eye(n)        , -1*ones(n,1);
      Aineq(ineq_ind,:)  ,  zeros(nineqcstr,1)];
bineq2 = [XIN; -XIN; bineq(ineq_ind)];
% We want extra variable to be >= 0
lb2 = [lb(:); 0.0];
ub2 = [ub(:); Inf];

Aeq2 = [Aeq(eq_ind,:) zeros(neqcstr,1)];
beq2 = beq(eq_ind);


% LINPROG
[XOUT,~,success] = linprog(f,Aineq2,bineq2,Aeq2,beq2,lb2,ub2,[],optimset('Display','off'));
XIN = XOUT(1:end-1); % Last element is a slack variable

