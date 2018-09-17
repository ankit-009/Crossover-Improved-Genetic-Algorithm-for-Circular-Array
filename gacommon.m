function [x,fval,exitFlag,output,population,scores,FitnessFcn,nvars,Aineq,bineq,Aeq,beq,lb,ub, ...
    NonconFcn,options,Iterate,type] = gacommon(nvars,fun,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,intcon,options,user_options,output)
%gacommon Common validation tasks for GA and GAMULTIOBJ
%
%   This function is private to GA and GAMULTIOBJ

%   Copyright 2007-2012 The MathWorks, Inc.
%   $Revision: 1.1.6.6 $  $Date: 2012/08/21 00:23:45 $

% Determine the population type
stringSet('PopulationType',options.PopulationType,{'doubleVector','custom','bitString'});

% Determine the 'type' of the problem
if ~isempty(nonlcon)
    type = 'nonlinearconstr';
    % Determine the sub-problem type for the constrained problem (used in ALPS)
    if (~isempty(Aeq)   && ~all(Aeq(:) == 0))   || ~isempty(beq) || ... 
       (~isempty(Aineq) && ~all(Aineq(:) == 0)) || ~isempty(bineq)
        subtype = 'linearconstraints';
    elseif (~isempty(lb) && any(isfinite(lb))) || ...
           (~isempty(ub) && any(isfinite(ub)))
        subtype = 'boundconstraints';
    else
        subtype = 'unconstrained';
    end
    % If Aeq or Aineq is not empty, then problem has linear constraints.
elseif (~isempty(Aeq)   && ~all(Aeq(:) == 0))   || ~isempty(beq) || ... 
       (~isempty(Aineq) && ~all(Aineq(:) == 0)) || ~isempty(bineq)
       type = 'linearconstraints';
    % This condition satisfies bound constraints
elseif (~isempty(lb) && any(isfinite(lb))) || ...
       (~isempty(ub) && any(isfinite(ub)))
    type = 'boundconstraints';
    % If all constraints fields are empty then it is unconstrained
else
    type = 'unconstrained';
end
% Make sure 'custom' or 'bitString' don't bring in unsupported constraints
% with them.
if strcmpi(options.PopulationType,'bitString') && ...
        ~strcmpi(type,'unconstrained')
    warning(message('globaloptim:gacommon:bitStringConstraintsSupport'));
    type = 'unconstrained';
    Aineq = []; bineq = []; Aeq = []; beq = []; lb = []; ub = [];
elseif strcmpi(options.PopulationType,'custom') && ...
        ~any(strcmpi(type,{'unconstrained','boundconstraints'}))
    warning(message('globaloptim:gacommon:customConstraintsSupport'));
    if (~isempty(lb) || ~isempty(ub))
        type = 'boundconstraints';
        Aineq = []; bineq = []; Aeq = []; beq = [];
    else
        type = 'unconstrained';
        Aineq = []; bineq = []; Aeq = []; beq = []; lb = []; ub = [];
    end
end
output.problemtype = type;
% If nonlinear constraints, then subtype is needed to process linear
% constraints (see function preProcessLinearConstr)
if strcmp(type,'nonlinearconstr')
    type = subtype;
end

% Remember the random number state used
dflt = RandStream.getGlobalStream;
output.rngstate = struct('state',{dflt.State}, 'type',{dflt.Type});
% Initialize other fields of output
output.generations = 0;
output.funccount   = 0;
output.message   = '';

% Validate options and fitness function
[options,nvars,FitnessFcn,NonconFcn] = validate(options,type,nvars,fun,nonlcon,user_options);

% Perform check on initial population, score, and range
options = initialPopulationCheck(options);

if ~strcmp(output.problemtype,'unconstrained') 
    % Determine a start point
    if ~isempty(options.InitialPopulation)
        x = options.InitialPopulation(1,:);
    else
        x = randn(1,nvars);
    end
    Iterate.x = x(:);
else
    Iterate.x = [];
end
% Initialize output
fval = [];
x = [];
population = [];
scores = [];

% Bound correction
[lb,ub,msg,exitFlag] = checkbound(lb,ub,nvars);
if exitFlag < 0
    output.message = msg;
    if options.Verbosity > 0
        fprintf('%s\n',msg)
    end
    return;
end
% Linear constraints correction
[Iterate.x,Aineq,bineq,Aeq,beq,lb,ub,msg,exitFlag] = ...
    preProcessLinearConstr(Iterate.x,Aineq,bineq,Aeq,beq,lb,ub,nvars,type,options.Verbosity);
if exitFlag < 0
    output.message = msg;
    if options.Verbosity > 0
        fprintf('%s\n',msg)
    end
    return;
end

% Verify that individuals in InitialPopulation are feasible
if ~isempty(options.InitialPopulation) && ~strcmp(type,'unconstrained')
    pop = size(options.InitialPopulation,1);
    feasible = true(pop,1);
    for i = 1:pop
        feasible(i) = isTrialFeasible(options.InitialPopulation(i,:)',Aineq,bineq,Aeq,beq,lb,ub,options.TolCon);
    end
    options.InitialPopulation(~feasible,:) = [];
    try % InitialScores may not be present
        options.InitialScores(~feasible,:) = [];
    catch
    end
end
% If initial population is empty at this point then we stuff the feasible point
% found before. Don't do this if integer constraints are present - we have
% found that this degrades the performance of some of the integer
% constrained problems.
if isempty(options.InitialPopulation) && ~isempty(Iterate.x) && isempty(intcon)
    options.InitialPopulation(1,:) = Iterate.x';
end

% Validate nonlinear constraints
[LinearConstr, Iterate,nineqcstr,neqcstr,ncstr] = constrValidate(NonconFcn, ...
    Iterate,Aineq,bineq,Aeq,beq,lb,ub,intcon,type,options);
options.LinearConstr = LinearConstr;

% Add integer constraint field to the options structure.
options.IntegerVars = [];

% Make sure that bounds and PopInitRange are consistent 
options.PopInitRange = checkPopulationInitRange(lb,ub,options.PopInitRange);

% Print some diagnostic information if asked for
if options.Verbosity > 2
    % Find numObj (for gamultiobj) by making one call to fitness function
    % (if Initial Score is not present)
    % Determine who is the caller
    callStack = dbstack;
    caller = callStack(2).file(1:end-2);
    if strcmp(caller,'gamultiobj')
        if isempty(options.InitialScores)
            if isempty(options.InitialPopulation)
                % Do we have a feasible point
                options.InitialPopulation(1,:) =  randn(1,nvars);
            end
            % Evaluate InitialPopulation(1,:)
             score = FitnessFcn(options.InitialPopulation(1,:));
             options.InitialScores = score(:)';
        end
        numObj = size(options.InitialScores,2);
    else
        numObj = 1;
    end
    gadiagnose(fun,nonlcon,nvars,nineqcstr,neqcstr,ncstr,numObj,user_options,caller);
end

