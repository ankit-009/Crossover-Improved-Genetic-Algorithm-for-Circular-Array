function [startPointsToRun, numFunEvals] = ...
    filterStartPoints(RunStartPoints,startPointsToRun,problem,numNlinIneqCon)
%filterStartPoints Filter the supplied start points.
%
%   [FILTPTS, NUMFUNEVALS] = filterStartPoints(RUNSTARTPTS, PTS, PROBLEM,
%   NUMNLININEQCON) filters the supplied start points, PTS, according to
%   the specified filter, RUNSTARTPTS. Possible values of RUNSTARTPTS are:-
%  
%   RUNSTARTPTS    |  FILTER ACTION
%   ---------------+---------------
%   'all'          |  All start points are allowed by the filter
%   'bounds'       |  Start points which satisfy the bounds are allowed by 
%                  |  the filter
%   'bounds-ineqs' |  Start points which satisfy the bounds, linear
%                  |  inequality and nonlinear inequality constraints are 
%                  |  allowed by the filter
%
%   See also FMULTISTART, GLOBALSEARCHNLP

%   Copyright 2009 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2012/08/21 00:23:36 $

% Initialize numFunEvals
numFunEvals = 0;

% No need to filter the start points when 'all' is specified. Just return.
if strcmp(RunStartPoints,'all')
    return
end

% If a filter has been specified, filter the supplied points
if any(strcmp(RunStartPoints,{'bounds','bounds-ineqs'}))

    numPoints = size(startPointsToRun,2);
    numVars = numel(problem.x0);

    % Evaluate bound constraints
    if isfield(problem.options,'TolCon')
        TOLCON = optimget(problem.options,'TolCon',1e-6);
    else
        TOLCON = 1e-6;
    end
    if ~isfield(problem,'lb') || isempty(problem.lb)
        lBnd = zeros(0, numPoints);
    else
        LB = problem.lb(:);
        lenLB = length(LB);
        if lenLB > numVars
            LB = LB(1:numVars);
        elseif lenLB < numVars
            LB = [LB; -inf*ones(numVars-lenLB,1)];
        end
        rLB = LB(:, ones(1, numPoints));
        lBnd = rLB - startPointsToRun;
    end
    if ~isfield(problem,'ub') || isempty(problem.ub)
        uBnd = zeros(0, numPoints);
    else
        UB = problem.ub(:);        
        lenUB = length(UB);        
        if lenUB > numVars
            UB = UB(1:numVars);
        elseif lenUB < numVars
            UB = [UB; inf*ones(numVars-lenUB,1)];
        end        
        rUB = UB(:, ones(1, numPoints));
        uBnd = startPointsToRun - rUB;
    end
    allCon = [lBnd;uBnd];
    feasIdx = allCon < TOLCON;
    
    % Evaluate inequality constraints
    if strcmp(RunStartPoints,'bounds-ineqs')
        
        % Evaluate the linear inequality constraints
        if ~isfield(problem,'Aineq') || isempty(problem.Aineq)
            linIneq = [];
        else
            rB = problem.bineq(:, ones(1, numPoints));
            linIneq = problem.Aineq*startPointsToRun - rB;
        end
        
        % Evaluate the nonlinear inequality constraints
        if ~isfield(problem,'nonlcon') || isempty(problem.nonlcon)
            nlinIneq = [];
        else
            nlinIneq = zeros(numNlinIneqCon, numPoints);
            for i = 1:numPoints
                thisNlinIneq = i_evaluateNonlcon(problem.nonlcon, ...
                    startPointsToRun(:, i), numNlinIneqCon);
                nlinIneq(:, i) = thisNlinIneq(:);
            end
            numFunEvals = numPoints;
        end
        allIneqCon = [linIneq;nlinIneq];
        % Concatenate the feasibility status of the inequality constraints
        % for each point with the feasibility status of the bounds.
        feasIdx = [feasIdx; allIneqCon < TOLCON];
    end
    
    % Filter the points
    keepIdx = all(feasIdx,1);
    startPointsToRun = startPointsToRun(:,keepIdx);
    
end

function c = i_evaluateNonlcon(nonlcon, x, numCon)

try
    % Try to evaluate the nonlinear inequality constraints. If this
    % evaluation encounters an error, or they evaluate to a complex number,
    % we consider the constraints to be undefined at this point.
    c = nonlcon(x);
    if ~isreal(c)
        error(message('globaloptim:filterStartPoints:UsrNonlConstrUndefAtX0'));
    end
catch ME
    % Set constraint evaluation to be Inf for this point
    c = Inf(numCon, 1);
end