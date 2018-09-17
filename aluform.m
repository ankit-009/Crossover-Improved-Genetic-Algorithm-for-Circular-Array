function [Iterate,A,L,U,nineqcstr,neqcstr,ncstr,IndIneqcstr,IndEqcstr,msg,exitflag] = ...
    aluform(XOUT,Aineq,Bineq,Aeq,Beq,LB,UB,numberOfVariables,type,verbosity,nonlcon,X,conFcnArg)
%ALUFORM  private to PATTERNSEARCH and GA.

%   Copyright 2003-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2012/10/24 04:17:59 $

% Initialize
Xin = XOUT;
Iterate.x = XOUT(:);
msg = '';
exitflag = 1;
tol = sqrt(eps);
% if subtype is unconstrained
if strcmpi(type,'unconstrained')
    A = [];L = [];U = []; nineqcstr = 0;neqcstr = 0; ncstr = 0;
    IndIneqcstr = [];IndEqcstr = [];
    % Evaluate nonlinear constraints for the first time
    if nargin >= 11 && ~isempty(nonlcon)
        try
            [cineq,ceq] = feval(nonlcon,reshapeinput(X,Iterate.x),conFcnArg{:});
             Iterate.cineq = zeros(numel(cineq),1);
             Iterate.ceq = zeros(numel(ceq),1);
             Iterate.cineq(:) = cineq;
             Iterate.ceq(:) = ceq;
        catch userFcn_ME
            gads_ME = MException('globaloptim:aluform:confunCheck', ...
                'Failure in initial user-supplied nonlinear constraint function evaluation.');
            userFcn_ME = addCause(userFcn_ME,gads_ME);
            rethrow(userFcn_ME)
        end
        c = [Iterate.cineq;Iterate.ceq];
        if ~all(isreal(c) & isfinite(c))
            error(message('globaloptim:aluform:confunNotReal'));
        end
    end
    return;
end

% If called from PFMINBND, do this bound check only
if strcmpi(type,'boundconstraints')
    A = eye(numberOfVariables);
    L = LB; U = UB;
    % Check the box constraints (bounds) first.
    XOUT = Iterate.x;
    lbound   = A*XOUT-LB >= 0;
    ubound   = A*XOUT-UB <= 0;
    feasible = all(lbound) && all(ubound);
    if ~feasible
        Iterate.x(~lbound) = LB(~lbound);
        Iterate.x(~ubound) = UB(~ubound);
        if verbosity > 2
            warning(message('globaloptim:aluform:infeasibleInitialPointBoundsOnly'));
        end
    end
    % Initialize these output arguments because this function might get
    % called with all output arguments (pfmincon)
    nineqcstr = 0;neqcstr = 0; ncstr = 0;
    IndIneqcstr = [];IndEqcstr = [];
    % Evaluate nonlinear constraints for the first time
    if nargin >= 11 && ~isempty(nonlcon)
        try
            [cineq,ceq] = feval(nonlcon,reshapeinput(X,Iterate.x),conFcnArg{:});
             Iterate.cineq = zeros(numel(cineq),1);
             Iterate.ceq = zeros(numel(ceq),1);
             Iterate.cineq(:) = cineq;
             Iterate.ceq(:) = ceq;
        catch userFcn_ME
            gads_ME = MException('globaloptim:aluform:confunCheck', ...
                'Failure in initial user-supplied nonlinear constraint function evaluation.');
            userFcn_ME = addCause(userFcn_ME,gads_ME);
            rethrow(userFcn_ME)
        end
        c = [Iterate.cineq;Iterate.ceq];
        if ~all(isreal(c) & isfinite(c))
            error(message('globaloptim:aluform:confunNotReal'));
        end
    end
    return;
end

% We allow row or column vectors for Beq and Bineq
Bineq = Bineq(:);
Beq   = Beq(:);
% Set the constraints up: defaults and check size
[nineqcstr,n] = size(Aineq);
[neqcstr,m]   = size(Aeq);
if ~isempty(Aineq)
    if ~isequal(length(Bineq),nineqcstr)
        error(message('globaloptim:aluform:inconsistentAineqAndBineq'))
    elseif ~isequal(numberOfVariables,n)
        error(message('globaloptim:aluform:inconsistentAineqAndX0'))
    end
elseif ~isempty(Bineq)
    error(message('globaloptim:aluform:emptyAineqNotBineq'))
end
if ~isempty(Aeq)
    if ~isequal(length(Beq),neqcstr)
        error(message('globaloptim:aluform:inconsistentAeqAndBeq'))
    elseif ~isequal(numberOfVariables,m)
        error(message('globaloptim:aluform:inconsistentAeqAndX0'))
    end
elseif ~isempty(Beq)
    error(message('globaloptim:aluform:emptyAeqNotBeq'))
end
% Remove dependent constraint, if any and find a basic solution w.r.t. A,L,U
[Iterate.x,Aineq,Bineq,Aeq,Beq,msg,how,exitflag]= eqnsolv(Iterate.x,Aineq,Bineq,Aeq,Beq,LB,UB,verbosity);

% Reinitialize these
nineqcstr = size(Aineq,1);
neqcstr   = size(Aeq,1);
ncstr     = nineqcstr + neqcstr;

% Create A, as described in L & T  L <= AX <=U.
Abox = eye(numberOfVariables);
A    = full([Aeq;Aineq;Abox]);
U    = full([Beq;Bineq;UB]);
L    = full([Beq;repmat(-Inf,nineqcstr,1);LB]);

% Logical indices of constraints.
IndIneqcstr = false(size(A,1),1);
IndEqcstr   = false(size(A,1),1);
IndEqcstr(1:neqcstr) = 1;
IndIneqcstr(neqcstr+1:ncstr) = 1;

% Is initial point feasible?
if strcmp(how,'infeasible')
    % Equalities are inconsistent, so return original X
    Iterate.x = Xin;
    return
end

% Check the box constraints (bounds).
XOUT = Iterate.x;
lbound   = Abox*XOUT-LB>=0;
ubound   = Abox*XOUT-UB<=0;
feasible = all(lbound) && all(ubound);
if ~feasible
    Iterate.x(~lbound) = LB(~lbound);
    Iterate.x(~ubound) = UB(~ubound);
    XOUT = Iterate.x;
    if verbosity > 2
        warning(message('globaloptim:aluform:infeasibleInitialPointBounds'));
    end
end

% Now add the linear constraints too and check it.
 feasible = isfeasible(XOUT,A,L,U,tol,IndEqcstr);

% Not a feasible point? find an initial feasible point using LP
if ~feasible
    if verbosity > 2
        warning(message('globaloptim:aluform:infeasibleInitialPointLinearConstr'));
    end
    % Find a feasible initial point using linprog-active-set algorithm
    [Iterate.x,~,success] = linprog([],Aineq,Bineq,Aeq,Beq,LB,UB,XOUT, ...
        optimset('Algorithm','active-set','Display','off'));
   feasible = isfeasible(Iterate.x,A,L,U,tol,IndEqcstr);
    if success <=0 || ~feasible
        % Add a slack variable and find a feasible initial point
        [Iterate.x,success] = initialfeasible(Iterate.x,A,L,U);
    end
    if success <= 0
        Iterate.x = Xin;
        msg = sprintf('Could not find a feasible initial point.');
        exitflag = -2;
        return;
    end
end
% Evaluate nonlinear constraints for the first time
if nargin >= 11 && ~isempty(nonlcon)
        try
            [cineq,ceq] = feval(nonlcon,reshapeinput(X,Iterate.x),conFcnArg{:});
             Iterate.cineq = zeros(numel(cineq),1);
             Iterate.ceq = zeros(numel(ceq),1);
             Iterate.cineq(:) = cineq;
             Iterate.ceq(:) = ceq;
        catch userFcn_ME
            gads_ME = MException('globaloptim:aluform:confunCheck', ...
                'Failure in initial user-supplied nonlinear constraint function evaluation.');
            userFcn_ME = addCause(userFcn_ME,gads_ME);
            rethrow(userFcn_ME)
        end
        c = [Iterate.cineq;Iterate.ceq];
        if ~all(isreal(c) & isfinite(c))
            error(message('globaloptim:aluform:confunNotReal'));
        end
end
