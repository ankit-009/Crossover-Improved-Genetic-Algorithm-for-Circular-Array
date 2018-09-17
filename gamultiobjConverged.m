function [state,exitFlag,reasonToStop] = gamultiobjConverged(options,state)
%gamultiobjConverged Check to see if any of the stopping criteria have been met.


%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2012/08/21 00:24:06 $


Gen = state.Generation;

if options.Verbosity > 1
    fprintf('%5.0f      %5.0f      %12.6g      %12.6g\n',Gen,state.FunEval,mean(state.AverageDistance),mean(state.Spread(end,:)));
end

% Window used to get best fval
 Window = options.StallGenLimit + 1;
 spreadChange = Inf;
 % Compute change in fval and individuals in last 'Window' generations
 if Gen > Window
     if size(state.Spread,2) > 1
        Spread =  mean(state.Spread((Gen - Window):end,:),2);
     else
         Spread = state.Spread((Gen - Window):end);
     end
     meanSpread = mean(Spread);
     spreadChange = 0;
     Weight = 0.5;
     for i = 1:Window-1
        spreadChange = spreadChange + (Weight)^(Window- i)*abs(Spread(i+1) - Spread(i))/(1+Spread(i));
     end
     % Take an average of function value change
     spreadChange = spreadChange/Window;
 end

reasonToStop = '';
exitFlag = [];
if(state.Generation >= options.Generations)
    reasonToStop = sprintf('Optimization terminated: maximum number of generations exceeded.');
    exitFlag = 0;
elseif((cputime-state.StartTime) > options.TimeLimit)
    reasonToStop = sprintf('Optimization terminated: time limit exceeded.');
    exitFlag = -5;
elseif(~isempty(state.StopFlag))
    reasonToStop = sprintf('Optimization terminated: %s',state.StopFlag);
    exitFlag = -1;
 elseif spreadChange <= options.TolFun && meanSpread >= Spread(end)
     reasonToStop = sprintf('Optimization terminated: average change in the spread of Pareto solutions less than options.TolFun.');
     exitFlag = 1;
end

% If it is a constrained problem and we want to check linear
% constraints (when GAMULIOBJ is terminating)
if ~isempty(reasonToStop) && isfield(options,'LinearConstr')
    linConstr = options.LinearConstr;
    if linConstr.linconCheck
        feasible = feasibleLinearConstraints(state,linConstr,sqrt(options.TolCon));
        if  ~all(feasible)
            if options.Verbosity > 0
                fprintf('%s\n','Infeasible individuals are present in the final population.');
            end
             state.Rank(~feasible) = -inf;
        end
        if ~all(isfinite(state.Rank))
            reasonToStop = [reasonToStop,...
                sprintf('\n%s','Linear constraints are not satisfied within constraint tolerance.')];
            exitFlag = -2;
        end
    end
end

if ~isempty(reasonToStop) && options.Verbosity > 0
    fprintf('%s\n',reasonToStop);
    return;
end
% Print header again
if options.Verbosity > 1 && rem(Gen,30)==0 && Gen > 0
    fprintf('\n                           Average            Average\n');
    fprintf('Generation   f-count    Pareto distance    Pareto spread\n');
end

%------------------------------------------------
function feasible  = feasibleLinearConstraints(state,linCon,tol)
% Function to check if linear constraints are satisfied at final
% population
popSize = size(state.Population,1);
feasible = false(popSize,1);
for i = 1:popSize
    feasible(i)  = isTrialFeasible(state.Population(i,:)',linCon.Aineq,linCon.bineq,linCon.Aeq, ...
        linCon.beq,linCon.lb,linCon.ub,tol);
end


