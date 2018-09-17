function [X,FVAL,EXITFLAG,output,POPULATION,SCORES] = gacon(FitnessFcn,GenomeLength, ...
    Aineq,bineq,Aeq,beq,lb,ub,NonconFcn,options,output,Iterate,subtype)
%GACON Genetic algorithm generalized constrained solver.
%   Private function to GA
%
%  GACON solves problems of the form:
%           min F(X)    subject to:    A*X <= b
%            X                         Aeq*X = beq
%                                      LB <= X <= UB
%                                      C(X) <= 0
%                                      Ceq(X) = 0
%

%   Copyright 2005-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.8 $  $Date: 2012/08/21 00:23:46 $


LinearConstr = options.LinearConstr;
% Number of linear constraints
nineqcstr = size(Aineq,1);
neqcstr   = size(Aeq,1);
% Initialize augmented Lagrangian related parameters
[gaAugParam,numNonlinIneqcstr,numNonlinEqcstr,numNonlinCstr] = gaAugInit(Iterate,options);
Iterate.f = FitnessFcn(Iterate.x');

% Parameter to track if nonlcon returns non-real value to fmincon
infMessage = [];
% Output structure initialized for sub-problem solvers
innerOutput.problemtype = subtype;
run = true;

% Fitness function formulation as augmented Lagrangian
SubFitness = createAnonymousFcn(@augLagFun,{FitnessFcn,NonconFcn,1,GenomeLength, ...
    gaAugParam,numNonlinIneqcstr,numNonlinEqcstr,numNonlinCstr, ...
    options.TolCon,options.FitnessFcnArgs,options.NonconFcnArgs});

% Create initial state: population, scores, status data. State is merely
% used to pass some values to output/plot functions. It is not used in
% ALGA.
state = makeState(GenomeLength,SubFitness,Iterate,subtype,options);
innerPopulation = state.Population;
state.Best(1) = state.Score(1);
X = Iterate.x';
% Constraint information in state
state.NonlinIneq = Iterate.cineq;
state.NonlinEq =   Iterate.ceq;

% Give the plot/output Fcns a chance to do any initialization they need.
state = gadsplot(options,state,'init','Genetic Algorithm');
[state,options] = gaoutput(FitnessFcn,options,state,'init');
% Setup display header
if  options.Verbosity > 1
    fprintf('\n                           Best       max        Stall\n');
    fprintf('Generation  f-count        f(x)     constraint  Generations\n');
end

% Outer loop setup the augmented Lagrangian formulation and inner loop
% solves the sub-problem
while run
    % Check for augmented Lagrangian convergence
    [X,FVAL,EXITFLAG,maxConstr,output.message,run] = gaAugConverged(options,state,gaAugParam,Iterate, ...
        numNonlinCstr,numNonlinIneqcstr,X,run,infMessage);
    if ~run, break; end

    % Fitness function formulation as augmented Lagrangian
    SubFitness = createAnonymousFcn(@augLagFun,{FitnessFcn,NonconFcn,1,GenomeLength, ...
        gaAugParam,numNonlinIneqcstr,numNonlinEqcstr,numNonlinCstr, ...
        options.TolCon,options.FitnessFcnArgs,options.NonconFcnArgs});    
    
    % Update Iterate with value from sub-problem
    Iterate.f = SubFitness(Iterate.x');

    % Reset some parameters before solving sub-problem
    innerOptions = gaAugReset(innerPopulation,options,state,gaAugParam);

    % Solve the sub-problem
    if ~strcmp(gaAugParam.step,'Infeasible')
        if any(strcmp(subtype,{'boundconstraints','linearconstraints'}))
            [~,Iterate.f,innerExitFlag,innerOutput,innerPopulation,innerScore] = ...
                galincon(SubFitness,GenomeLength,Aineq,bineq,Aeq,beq,lb,ub,innerOptions,innerOutput,Iterate);
        elseif strcmp(subtype,'unconstrained')
            [~,Iterate.f,innerExitFlag,innerOutput,innerPopulation,innerScore] = ...
                gaunc(SubFitness,GenomeLength,innerOptions,innerOutput,Iterate);
        end
    end
    % Which individual is the best (same as innerX)
    [Iterate.f,bestIndividual] = min(innerScore);
    % Find the best solution found for sub-problem
    Iterate.x = innerPopulation(bestIndividual,:)';

    % Evaluate constraints at the intermediate solution
    [Iterate.cineq(:),Iterate.ceq(:)] = NonconFcn(Iterate.x');
    % Update Lagrange multipliers and penalty parameters (based on
    % constraint values only)
    [gaAugParam,state] = gaAugUpdate(Iterate,gaAugParam,state,options,numNonlinIneqcstr,numNonlinCstr);
    % Make sure that intermediate point is feasible for next sub-problem
    if numNonlinIneqcstr
        shiftedConstr = gaAugParam.shift(1:numNonlinIneqcstr) - Iterate.cineq + options.TolCon;
        % Infeasible next point?
        if any(shiftedConstr <= 0)
            e = -2; foundFeasible = false; infMessage = [];
            infeasConfunSetup = @(x) infeasConstr(x,gaAugParam.shift);
            % Solve a problem to find a feasible solution
            warnstate = warning;
            warning off;
            try
                [x,f,e,o] = fmincon(@(x) x(end),[Iterate.x;options.PenaltyFactor], ...
                    [Aineq zeros(nineqcstr,1)],bineq,[Aeq zeros(neqcstr,1)],beq,[lb;0],[ub;Inf], ...
                    infeasConfunSetup,optimset('OutputFcn',@infeasOutputFcn,'Algorithm','active-set', ...
                    'FunValCheck','on','Display','off'));
                state.FunEval = state.FunEval + o.funcCount;
            catch fminconME
                infMessage = fminconME.message;
            end
            warning(warnstate);
            % Accept this condition as feasible iterate for next
            % iteration if f is close to 0
            if e >= -1 && f < 1
                foundFeasible = true;
            end

            % If exitflag is negative or 'f' is not close to zero,
            % then the problem does not have a feasible solution
            if e < -1 || ~foundFeasible
                gaAugParam.step = 'Infeasible';
                Iterate.f = FitnessFcn(Iterate.x');
            else
                Iterate.x(:) = x(1:end-1);
            end
        end
    end
    if isempty(gaAugParam.step) % The step is either empty or 'Infeasible'
        % Compute the constraints at the intermediate point and function value
        [Iterate.cineq(:),Iterate.ceq(:)] = NonconFcn(Iterate.x');
        % This should not happen
        if numNonlinIneqcstr
            shiftedConstr = gaAugParam.shift(1:numNonlinIneqcstr) - Iterate.cineq + options.TolCon;
            if any (shiftedConstr <= 0)
                % Adjust shift
                gaAugParam.shift(1:numNonlinIneqcstr) = max(0,Iterate.cineq) + 1e-4;
                state.how = 'Infeasible point';
            end
        end
        Iterate.f = FitnessFcn(Iterate.x');
    else % Restore the old Iterate
        Iterate.f = FVAL;
        Iterate.x = X(:);
        Iterate.cineq = state.NonlinIneq;
        Iterate.ceq = state.NonlinEq;
    end


    % Update the state
    innerPopulation(bestIndividual,:) = Iterate.x';
    state.Population = innerPopulation;
    state.Score = getFitnessValue(innerPopulation);
    best = Iterate.f;
    generation = state.Generation;
    state.Best(generation) = best;
    state.NonlinIneq = Iterate.cineq;
    state.NonlinEq =   Iterate.ceq;
    state.FunEval = state.FunEval + innerOutput.funccount;
    % Keep track of improvement in the best
    if((generation > 1) && isfinite(best))
        if(state.Best(generation-1) ~= best)
            state.LastImprovement = generation;
            state.LastImprovementTime = cputime;
        end
    end
    % Special case: If solver is stopped by plot/output function in the
    % sub-problem we terminate
    if innerExitFlag == -1
        state.StopFlag = 'Stopped in inner iterations.';
        continue;
    end
    % Update plot and output functions
    state = gadsplot(options,state,'iter','Genetic Algorithm');
    [state,options,optchanged] = gaoutput(FitnessFcn,options,state,'iter');
    if optchanged
        options.LinearConstr = LinearConstr;
    end

end % End outer while loop

% Update output structure
output.generations = state.Generation;
output.funccount   = state.FunEval;
output.maxconstraint = maxConstr;
POPULATION = state.Population;
SCORES = state.Score;


% Call hybrid function
if ~isempty(options.HybridFcn)
    if strcmpi(options.PopulationType,'doubleVector')
        [X,FVAL] = callHybridFunction;
    else
        warning(message('globaloptim:gacon:notValidHybrid'));
    end
end
% Give the Output functions a chance to finish up
gadsplot(options,state,'done','Genetic Algorithm');
gaoutput(FitnessFcn,options,state,'done');

%---------------------------------------------------------------
% Constraint formulation to handle infeasiblity problem
    function [cin,ceq] = infeasConstr(input,oldShift)
        x = input(1:end-1);
        cin = zeros(numNonlinIneqcstr,1);
        [cin(:),ceq] = NonconFcn(x');
        cin = cin - input(end)*oldShift;
    end % End of infeasConstr
%---------------------------------------------------------------
% Output function to stop the infeasibility problem when fitness is less
% than one.
    function stop = infeasOutputFcn(~,optimValues,flag)
        stop = false;
        switch flag
            case {'init','iter'}
                if optimValues.fval < 1 && ...
                        optimValues.constrviolation <= options.TolCon
                    stop = true;
                end
            otherwise
        end
    end
%-----------------------------------------------------------------
    function obj = getFitnessValue(pop)
        % Calculate fitness of population
        if strcmpi(options.Vectorized, 'off')
            obj = fcnvectorizer(pop,FitnessFcn,1,options.SerialUserFcn);
        else
            obj = FitnessFcn(pop);
        end
    end
%-----------------------------------------------------------------
% Hybrid function
    function [xhybrid,fhybrid] = callHybridFunction
        xhybrid = X;
        fhybrid = FVAL;
        % Who is the hybrid function
        if isa(options.HybridFcn,'function_handle')
            hfunc = func2str(options.HybridFcn);
        else
            hfunc = options.HybridFcn;
        end
        % Inform about hybrid scheme
        if  options.Verbosity > 1
            fprintf('%s%s%s\n','Switching to the hybrid optimization algorithm (',upper(hfunc),').');
        end
        % Create functions handle to be passed to hybrid function
        FitnessHybridFcn = @(x) FitnessFcn(x,options.FitnessFcnArgs{:});
        if ~isempty(NonconFcn)
            ConstrHybridFcn  = @(x) NonconFcn(x,options.NonconFcnArgs{:});
        else
            ConstrHybridFcn = [];
        end
        if ~any(strcmpi(hfunc,{'fmincon', 'patternsearch'}))
            warning(message('globaloptim:gacon:unconstrainedHybridFcn', ...
                upper(hfunc),'FMINCON'));
            hfunc = 'fmincon';
        end
        [x_temp,f_temp,funccount,theMessage,conviol_temp] = ...
            callHybrid(hfunc,FitnessHybridFcn,X,options.HybridFcnArgs,Aineq,bineq,Aeq,beq,lb,ub,ConstrHybridFcn);
        output.funccount = output.funccount + funccount;
        output.message   = sprintf([output.message '\n', theMessage '\n']);
        hybridPointFeasible = isHybridPointFeasible(conviol_temp, ...
            hfunc, options.HybridFcnArgs{:}); 
        gaPointFeasible = (output.maxconstraint < options.TolCon);
        if hybridPointFeasible && (~gaPointFeasible || f_temp < fhybrid)
            fhybrid = f_temp;
            xhybrid = x_temp;
        end
        % Inform about hybrid scheme termination
        if  options.Verbosity > 1
            fprintf('%s%s\n',upper(hfunc), ' terminated.');
        end
    end % End of callHybridFunction
end  % End of gaconstr.m

