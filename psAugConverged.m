function [X,FVAL,EXITFLAG,maxConstr,msg,run] = psAugConverged(Iterate,X,optimState, ...
    psAugParam,options,numNonlinCstr,numNonlinIneqcstr,EXITFLAG,run)
%PSAUGCONVERGED Augmented Lagrangian barrier convergence test.
% Private to PATTERNSEARCH

%   Copyright 2005-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2012/08/21 00:24:53 $

verbosity = options.Verbosity;
Iter = optimState.Iter;
FunEval = optimState.FunEval;
how = optimState.how;
deltaX = optimState.deltaX;
TolX = options.TolX;
currentTolMesh = psAugParam.currentTolMesh;
minMesh = options.TolMesh;
infMessage = optimState.infMessage;
TolCon = options.TolCon;

X(:) = Iterate.x;
FVAL = Iterate.f;
comp_slack = 0;
maxConstr = 0;
msg = '';
if numNonlinIneqcstr
    comp_slack = norm(Iterate.cineq.*psAugParam.lambdabar(1:numNonlinIneqcstr));
    maxConstr = max([maxConstr;Iterate.cineq(:)]);
end
if numNonlinCstr > numNonlinIneqcstr
    maxConstr = max([maxConstr;abs(Iterate.ceq(:))]);
end
% Print some iterative information
if verbosity > 1
    fprintf('%5.0f %5.0f %12.6g %12.4g %12.4g   %s', ...
        Iter,FunEval,FVAL,maxConstr,currentTolMesh,how);
    fprintf('\n');
end
stallTol = min(minMesh,eps);

% Check mesh size tolerance and complementary slackness
if currentTolMesh <= minMesh && comp_slack <= sqrt(TolCon) && maxConstr <= TolCon
    EXITFLAG = 1;
    run  = false;
    msg = sprintf('%s','Optimization terminated: ');
    msg = [msg,sprintf('%s', 'mesh size less than options.TolMesh')];
    msg = [msg,sprintf('\n%s', ' and constraint violation is less than options.TolCon.')];
    if verbosity > 0
        fprintf('%s\n',msg);
    end
    return;
end

% Constraints are satisfied, deltaX is small, and mesh size is small enough
if deltaX <= TolX && currentTolMesh <= TolX && maxConstr <= TolCon
    EXITFLAG = 2;
    run  = false;
    msg = sprintf('%s','Optimization terminated: Change in X less than options.TolX');
    msg = [msg, sprintf('\n%s', ' and constraints violation is less that options.TolCon.')];
    if verbosity > 0
        fprintf('%s\n',msg);
    end
    return;
end

% Constraints are satisfied but step is too small
if currentTolMesh <= stallTol && maxConstr <= TolCon
    EXITFLAG = 4;
    run  = false;
    msg = sprintf('%s %g','Optimization terminated: norm of the step is less than ',stallTol);
    msg = [msg, sprintf('\n%s', ' and constraints violation is less that options.TolCon.')];
    if verbosity > 0
        fprintf('%s\n',msg);
    end
    return;
end

% Mesh size below tolerance and no improvement in X
if (deltaX <= eps && currentTolMesh <= minMesh)
    run  = false;
    msg = sprintf('%s %g','Optimization terminated: norm of the step is less than ',stallTol);
    % Check feasibility
    if  maxConstr <= TolCon
        EXITFLAG = 4;
        msg = [msg, sprintf('\n%s', ' and constraints violation is less that options.TolCon.')];
    else % Stall in infeasible region
        EXITFLAG = -2;
        msg = sprintf('%s\n','Optimization terminated: no feasible point found.');
    end
    if verbosity > 0
        fprintf('%s\n',msg);
    end
    return;
end

% fmincon enocntered NaN or Inf and could not continue; error here
if ~isempty(infMessage) & strmatch('optimlib:optimfcnchk',infMessage)
    error(message('globaloptim:psAugConverged:NaNFval'));
end
% No feasible solution or stall
if strcmpi(psAugParam.step,'Infeasible')
    % Check feasibility
    if  (maxConstr <= TolCon)
         return; % This will stop in later iteration due to TolX
    else % Stall in infeasible region
        run = false;
        EXITFLAG = -2;
        msg = sprintf('%s\n','Optimization terminated: no feasible point found.');
    end
    if verbosity > 0
        fprintf('%s\n',msg);
    end
    return;
end

% User interruption
if optimState.stopOutput || optimState.stopPlot
    EXITFLAG = -1;
    run = false;
    msg = sprintf('%s','Stop requested.');
    if verbosity > 0
        fprintf('%s\n',msg);
    end
    return;
end

% Maxiter test
if Iter > options.MaxIter
    EXITFLAG = 0;
    run  = false;
    msg = sprintf('%s', 'Maximum number of iterations exceeded: ');
    msg = [msg,sprintf('%s', 'increase options.MaxIter.')];
    if verbosity > 0
        fprintf('%s\n',msg);
    end
    return;
end
% max Fun Evaluation test
if FunEval >=  options.MaxFunEvals
    EXITFLAG = 0;
    run  = false;
    msg = sprintf('%s', 'Maximum number of function evaluations exceeded: ');
    msg = [msg, sprintf('%s', 'increase options.MaxFunEvals.')];
    if verbosity > 0
        fprintf('%s\n',msg);
    end
    return;
end

% Max time limit test
if(cputime- optimState.StartTime) > options.TimeLimit
    EXITFLAG = 0;
    run  = false;
    msg = sprintf('%s', 'Time limit reached: ');
    msg = [msg, sprintf('%s', 'increase options.TimeLimit.')];
    if verbosity > 0
        fprintf('%s\n',msg);
    end
    return;
end
% Setup display header every twenty iterations
if verbosity > 1 && rem(Iter,20)== 0 && Iter > 0
    fprintf('\n                                  max\n');
    fprintf('  Iter   f-count      f(x)      constraint   MeshSize      Method\n');
end
