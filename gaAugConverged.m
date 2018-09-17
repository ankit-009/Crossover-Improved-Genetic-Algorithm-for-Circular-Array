function [X,Fval,exitFlag,maxConstr,reasonToStop,run] = gaAugConverged(options,state,gaAugParam,Iterate, ...
    numNonlinCstr,numNonlinIneqcstr,X,run,infMessage)
%GAAUGCONVERGED Augmented Lagrangian convergence test.
%   Private to GA

%   Copyright 2005-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2012/08/21 00:23:41 $


X(:) = Iterate.x;
Fval = Iterate.f;
exitFlag = [];
reasonToStop = '';
comp_slack = 0;
maxConstr = 0;
Gen = state.Generation;
stallTol = min(options.TolFun,eps);
% Calculate complementary condition and constraint violation
if numNonlinIneqcstr
    comp_slack = norm(Iterate.cineq.*gaAugParam.lambdabar(1:numNonlinIneqcstr));
    maxConstr = max([maxConstr;Iterate.cineq(:)]);
end
if numNonlinCstr > numNonlinIneqcstr
    maxConstr = max([maxConstr;abs(Iterate.ceq(:))]);
end
% Print iterative information
if  options.Verbosity > 1 && Gen > 0
    FunEval  = state.FunEval;
    StallGen = Gen - state.LastImprovement;
    fprintf('%5.0f       %5.0f  %12.6g %12.4g    %3.0f\n', ...
        Gen, FunEval, Fval, maxConstr, StallGen);
end
% Converged at an optimum
if gaAugParam.currentTolFun < options.TolFun && maxConstr <= options.TolCon && ...
        comp_slack <= sqrt(options.TolCon)
    run = false;
    reasonToStop = sprintf('%s','Optimization terminated: ');
    reasonToStop = [reasonToStop,sprintf('average change in the fitness value less than options.TolFun')];
    % Check if linear constraints are satisfied
    linCon = options.LinearConstr;
    if ~isempty(linCon) && linCon.linconCheck && ~feasibleLinearConstraints
        reasonToStop = [reasonToStop,sprintf('\n%s',' but linear constraints are not satisfied.')];
        exitFlag = -2;
    else
        reasonToStop = [reasonToStop,sprintf('\n%s', ' and constraint violation is less than options.TolCon.')];
        exitFlag = 1;
    end
    if options.Verbosity > 0
        fprintf('%s\n',reasonToStop);
    end
    return;
end
% Stall but constraints are satisfied
if gaAugParam.currentTolFun < stallTol && maxConstr <= options.TolCon
    run = false;
    reasonToStop = sprintf('%s %6.5g','Optimization terminated: norm of the step is less than ',eps);
    % Check if linear constraints are satisfied
    linCon = options.LinearConstr;
    if ~isempty(linCon) && linCon.linconCheck && ~feasibleLinearConstraints
        reasonToStop = [reasonToStop,sprintf('\n%s',' but linear constraints are not satisfied.')];
        exitFlag = -2;
    else
        reasonToStop = [reasonToStop, sprintf('\n%s', ' and constraint violation is less than options.TolCon.')];
        exitFlag = 4;
    end
    if options.Verbosity > 0
        fprintf('%s\n',reasonToStop);
    end
    return;
end
% fmincon enocntered NaN or Inf and could not continue; error here
if ~isempty(infMessage) && ~isempty(strmatch('optimlib:optimfcnchk',infMessage))
    error(message('globaloptim:GAAUGCONVERGED:NaNFval', 'Constraint function returned non-real value;can not continue.'));
end
% Stalls and constraints are not satisfied
if strcmpi(gaAugParam.step,'Infeasible')
    run  = false;
    reasonToStop = sprintf('Optimization terminated: no feasible point found.');
    exitFlag = -2;
    if options.Verbosity > 0
        fprintf('%s\n',reasonToStop);
    end
    return;
end
% Maximum generations
if(Gen >= options.Generations)
    reasonToStop = sprintf('Optimization terminated: maximum number of generations exceeded.');
    exitFlag = 0;
    run  = false;
    if options.Verbosity > 0
        fprintf('%s\n',reasonToStop);
    end
    return;
end
% Maximum time limit
if((cputime-state.StartTime) > options.TimeLimit)
    reasonToStop = sprintf('Optimization terminated: time limit exceeded.');
    exitFlag = -5;
    run  = false;
    if options.Verbosity > 0
        fprintf('%s\n',reasonToStop);
    end
    return;
end
% Stall time limit
if((cputime-state.LastImprovementTime) > options.StallTimeLimit)
    reasonToStop = sprintf('Optimization terminated: stall time limit exceeded.');
    exitFlag = -4;
    run  = false;
    if options.Verbosity > 0
        fprintf('%s\n',reasonToStop);
    end
    return;
end
% Stall generation limit
if((state.Generation  - state.LastImprovement) > options.StallGenLimit)
    reasonToStop = sprintf('Optimization terminated: stall generations limit exceeded');
    if maxConstr <= options.TolCon
        % Next, check if linear constraints are satisfied
        linCon = options.LinearConstr;
        if ~isempty(linCon) && linCon.linconCheck && ~feasibleLinearConstraints
            reasonToStop = [reasonToStop,sprintf('\n%s',' but linear constraints are not satisfied.')];
            exitFlag = -2;
        else
            reasonToStop = [reasonToStop, sprintf('\n%s', ' and constraint violation is less than options.TolCon.')];
            exitFlag = 3;
        end
    else
        reasonToStop = [reasonToStop,sprintf('\n%s',' but constraints are not satisfied.')];
        exitFlag = -2;
    end
    run  = false;
    if options.Verbosity > 0
        fprintf('%s\n',reasonToStop);
    end
    return;
end
% Minimum fitness limit
if (Fval <= options.FitnessLimit )
    reasonToStop = sprintf('Optimization terminated: minimum fitness limit reached');
    if maxConstr <= options.TolCon
        % Next, check if linear constraints are satisfied
        linCon = options.LinearConstr;
        if ~isempty(linCon) && linCon.linconCheck && ~feasibleLinearConstraints
            reasonToStop = [reasonToStop,sprintf('\n%s',' but linear constraints are not satisfied.')];
            exitFlag = -2;
        else
            reasonToStop = [reasonToStop, sprintf('\n%s', ' and constraint violation is less than options.TolCon.')];
            exitFlag = 5;
        end
    else
        reasonToStop = [reasonToStop,sprintf('\n%s',' but constraints are not satisfied.')];
        exitFlag = -2;
    end
    run  = false;
    if options.Verbosity > 0
        fprintf('%s\n',reasonToStop);
    end
    return;
end
% Stop requested from user
if(~isempty(state.StopFlag))
    reasonToStop = sprintf('Optimization terminated: %s',state.StopFlag);
    exitFlag = -1;
    run  = false;
    if options.Verbosity > 0
        fprintf('%s\n',reasonToStop);
    end
    return;
end

% Print header again
if options.Verbosity > 1 && rem(Gen,20)==0 && Gen > 0
    fprintf('\n                           Best       max        Stall\n');
    fprintf('Generation  f-count        f(x)     constraint  Generations\n');
end
%---------------------------------------------------------
    function feasible  = feasibleLinearConstraints
        % Function to check if linear constraints are satisfied at final point
        % If it is a constrained problem and we want to check linear constraints
        tol = max(sqrt(eps),sqrt(options.TolCon));
        feasible = isTrialFeasible(Iterate.x,linCon.Aineq,linCon.bineq,linCon.Aeq, ...
        linCon.beq,linCon.lb,linCon.ub,tol);
    end % End of feasibleLinearConstraints
end
