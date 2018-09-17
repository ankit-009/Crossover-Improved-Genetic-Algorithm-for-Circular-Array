function [XOUT,Aineq,Bineq,Aeq,Beq,LB,UB,msg,exitflag] = ...
    preProcessLinearConstr(XOUT,Aineq,Bineq,Aeq,Beq,LB,UB,numberOfVariables,type,verbosity)
%PREPROCESSLINEARCONSTR validates dimension of constraint matrices, removes
%   redundancy in them, and finds initial feasible point.
%
%   private to PATTERNSEARCH and GA.

%   Copyright 2003-2012 The MathWorks, Inc.
%   $Revision: 1.1.6.6 $  $Date: 2012/10/24 04:18:02 $

% Initialize
msg = '';
exitflag = 1;

% If problem type is unconstrained return here
if strcmpi(type,'unconstrained')
    return;
end

% If called from PFMINBND, do this bound check only
if strcmpi(type,'boundconstraints')
    % Check the box constraints (bounds) first.
    lbound   = XOUT-LB >= 0;
    ubound   = XOUT-UB <= 0;
    feasible = all(lbound) && all(ubound);
    if ~feasible
        XOUT(~lbound) = LB(~lbound);
        XOUT(~ubound) = UB(~ubound);
    end
    return;
end

% ------------Only linear constrained problems reach here--------------
% Initialize
tol = sqrt(eps);
Xin = XOUT;

% We allow row or column vectors for Beq and Bineq
Bineq = Bineq(:);
Beq   = Beq(:);
% Set the constraints up: defaults and check size
[nineqcstr,n] = size(Aineq);
[neqcstr,m]   = size(Aeq);
if ~isempty(Aineq)
    if ~isequal(length(Bineq),nineqcstr)
        error(message('globaloptim:preprocesslinearconstr:inconsistentAineqAndBineq'))
    elseif ~isequal(numberOfVariables,n)
        error(message('globaloptim:preprocesslinearconstr:inconsistentAineqAndX0'))
    end
    % Remove empty constraints.
    zeroConstrID = all(Aineq == 0,2); % get rows with all zeros
    zeroConstrID = zeroConstrID & (Bineq(:) == 0);
    if ~isempty(zeroConstrID)
        Aineq(zeroConstrID,:) = [];
        Bineq(zeroConstrID) = [];
    end
elseif ~isempty(Bineq)
    error(message('globaloptim:preprocesslinearconstr:emptyAineqNotBineq'))
end
if ~isempty(Aeq)
    if ~isequal(length(Beq),neqcstr)
        error(message('globaloptim:preprocesslinearconstr:inconsistentAeqAndBeq'))
    elseif ~isequal(numberOfVariables,m)
        error(message('globaloptim:preprocesslinearconstr:inconsistentAeqAndX0'))
    end
    % Remove empty constraints.
    zeroConstrID = all(Aeq == 0,2); % get rows with all zeros
    zeroConstrID = zeroConstrID & (Beq(:) == 0);
    if ~isempty(zeroConstrID)
        Aeq(zeroConstrID,:) = [];
        Beq(zeroConstrID) = [];
    end
elseif ~isempty(Beq)
    error(message('globaloptim:preprocesslinearconstr:emptyAeqNotBeq'))
end

% Remove dependent constraint, if any and find a basic solution
[XOUT,Aineq,Bineq,Aeq,Beq,msg,how,exitflag]= eqnsolv(XOUT,Aineq,Bineq,Aeq,Beq,LB,UB,verbosity);

% Turn constraints into right size if they are empty after pre-processing.
if isempty(Aineq)
    Aineq = zeros(0,numberOfVariables);
end
if isempty(Bineq)
    Bineq = zeros(0,1);
end
if isempty(Aeq)
    Aeq = zeros(0,numberOfVariables); 
end
if isempty(Beq)
    Beq = zeros(0,1);
end


% Is initial point feasible?
if strcmp(how,'infeasible')
    % Equalities are inconsistent, so return original X
    XOUT = Xin;
    return
end

%------------- Find a feasible point-------------
% Check the bound constraints
lbound   = XOUT - LB >= 0;
ubound   = XOUT - UB <= 0;
feasible = all(lbound) && all(ubound);
if ~feasible
    XOUT(~lbound) = LB(~lbound);
    XOUT(~ubound) = UB(~ubound);
end

% Now add the linear constraints too and check it
feasible = isTrialFeasible(XOUT,Aineq,Bineq,Aeq,Beq,[],[],tol);

% Not a feasible point? find an initial feasible point using LP (Two
% approaches)
if ~feasible
    % Find a feasible initial point using linprog active-set algorithm
    [XOUT,~,success] = linprog([],Aineq,Bineq,Aeq,Beq,LB,UB,XOUT, ...
        optimset('Algorithm','active-set','Display','off'));
    feasible = isTrialFeasible(XOUT,Aineq,Bineq,Aeq,Beq,LB,UB,tol);
    if success <=0 || ~feasible
        % Add a slack variable and find a feasible initial point
        [XOUT,success] = initialfeasible(XOUT,numberOfVariables,Aineq,Bineq,Aeq,Beq,LB,UB);
        if success > 0
            feasible = isTrialFeasible(XOUT,Aineq,Bineq,Aeq,Beq,LB,UB,tol);
        end
    end
    
    % Add additional linprog algorithm
    if success <= 0 || ~feasible
        % Find a feasible initial point using lipsol algorithm
        [XOUT,~,success] = linprog([],Aineq,Bineq,Aeq,Beq,LB,UB,XOUT, ...
            optimset('Display','off'));
        if success > 0
            feasible = isTrialFeasible(XOUT,Aineq,Bineq,Aeq,Beq,LB,UB,tol);
        end
    end
    if success <= 0 || ~feasible
        % Find a feasible initial point using Simplex algorithm
        [XOUT,~,success] = linprog([],Aineq,Bineq,Aeq,Beq,LB,UB,XOUT, ...
            optimset('Algorithm','simplex','Display','off'));
        if success > 0        
            feasible = isTrialFeasible(XOUT,Aineq,Bineq,Aeq,Beq,LB,UB,tol);
        end
    end
    
    % Quit now
    if success <= 0 || ~feasible
        XOUT = Xin;
        msg = sprintf('Could not find a feasible initial point.');
        exitflag = -2;
    end
end

