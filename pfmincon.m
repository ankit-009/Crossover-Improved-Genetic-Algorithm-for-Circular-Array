function [X,FVAL,EXITFLAG,OUTPUT] = pfmincon(FUN,objFcnArg,initialX,numberOfVariables,Iterate, ...
    Aineq,bineq,Aeq,beq,NullBasisAeq,lb,ub,nonlcon,conFcnArg,options,defaultopt,OUTPUT,subproblemtype)
%PFMINCON Finds a constrained minimum of a function.
%   PFMINCON solves problems of the form:
%           min F(X)    subject to:    A*X <= b
%            X                         Aeq*X = beq
%                                      LB <= X <= UB
%                                      C(X) <= 0
%                                      Ceq(X) = 0
%

%   Copyright 2005-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2012/08/21 00:24:40 $



% Number of linear constraints
nineqcstr = size(Aineq,1);
neqcstr   = size(Aeq,1);
[row,col] = size(initialX);

% Get some initial values
[optimState,nextIterate,unused,EXITFLAG,run] = ...
    getinitial(Iterate,numberOfVariables,neqcstr,lb,ub,options);

X = initialX;
X(:) = Iterate.x;
FVAL = Iterate.f;

% Initialize augmented Lagrangian related parameters
[psAugParam,numNonlinIneqcstr,numNonlinEqcstr,numNonlinCstr] = psAugInit(Iterate,options);

% Output structure initialized for sub-problem solvers
innerOutput = struct('function',FUN,'problemtype',subproblemtype);

% Set up output function plot
if options.OutputTrue || options.PlotTrue
    callOutputPlotFunctions('init');
end

% Setup display header
if  options.Verbosity>1
    fprintf('\n                                  max\n');
    fprintf('  Iter   f-count      f(x)      constraint   MeshSize      Method\n');
end
% Setup and update sub-problem; call private solvers to solve sub-problem
while run
    % Check for augmented Lagrangian convergence
    [X,FVAL,EXITFLAG,maxConstr,msg,run] = psAugConverged(Iterate,X,optimState,psAugParam, ...
        options,numNonlinCstr,numNonlinIneqcstr,EXITFLAG,run);
    if ~run, break; end

    % Objective function formulation as augmented Lagrangian
     SubObjective = createAnonymousFcn(@augLagFun, {FUN,nonlcon,row,col,psAugParam,numNonlinIneqcstr,numNonlinEqcstr, ...
        numNonlinCstr,options.TolCon,objFcnArg,conFcnArg});
    
    % Reset some parameters before solving sub-problem
    innerOptions = psAugReset(options,optimState,psAugParam);
    % Function value is calculated as sub-problem function value
    Iterate.f = SubObjective(reshapeinput(X,Iterate.x));
    %------Solve the sub-problem---------%
    if ~strcmp(psAugParam.step,'Infeasible')
        if strcmp(subproblemtype,'linearconstraints')
            [nextIterate.x(:),nextIterate.f,innerExitFlag,innerOutput] = pfminlcon(SubObjective,{},X,numberOfVariables,Iterate, ...
                Aineq,bineq,Aeq,beq,NullBasisAeq,lb,ub,innerOptions,defaultopt,innerOutput);
        elseif strcmp(subproblemtype,'boundconstraints')
            [nextIterate.x(:),nextIterate.f,innerExitFlag,innerOutput] = pfminbnd(SubObjective,{},X,numberOfVariables,Iterate, ...
                lb,ub,innerOptions,defaultopt,innerOutput);
        elseif strcmp(subproblemtype,'unconstrained')
            [nextIterate.x(:),nextIterate.f,innerExitFlag,innerOutput] = pfminunc(SubObjective,{},X,numberOfVariables,Iterate, ...
                innerOptions,defaultopt,innerOutput);
        end
    end
    % Update function evaluation counter
    optimState.FunEval = optimState.FunEval + innerOutput.funccount;
    % Compute the constraints at the intermediate point but no need to
    % calculate function value
    [nextIterate.cineq(:),nextIterate.ceq(:)] = feval(nonlcon,reshapeinput(X,nextIterate.x),conFcnArg{:});
    % Update Lagrange multipliers and penalty parameters (based on
    % constraint values only)
    [psAugParam,optimState] = psAugUpdate(nextIterate,psAugParam,optimState,options, ...
        numNonlinIneqcstr,numNonlinCstr);

    %------Make sure the intermediate solution is feasible------
    if numNonlinIneqcstr
        shiftedConstr = psAugParam.shift(1:numNonlinIneqcstr) - nextIterate.cineq + options.TolCon;
        % Infeasible point for next sub-problem?
        if any(shiftedConstr <= 0)
            e = -2; foundFeasible = false; optimState.infMessage = [];
           infeasConfunSetup = @(x) infeasConstr(x,psAugParam.shift);
            % Find a feasible solution
            warnstate = warning;
            warning off;
            try
                [x,f,e,o] = fmincon(@(x) x(end),[nextIterate.x;options.PenaltyFactor], ...
                    [Aineq zeros(nineqcstr,1)],bineq,[Aeq zeros(neqcstr,1)],beq,[lb;0],[ub;Inf], ...
                    infeasConfunSetup,optimset('OutputFcn',@infeasOutputFcn,'Algorithm','active-set', ...
                    'FunValCheck','on','Display','off','MaxFunEvals', ...
                            options.MaxFunEvals - optimState.FunEval,'TolCon', options.TolCon));
                optimState.FunEval = optimState.FunEval + o.funcCount;
            catch optim_ME
                optimState.infMessage = optim_ME.message;
                %[unused,optimState.infMessage] = lasterr;
            end
            warning(warnstate);
            % Even if f is significantly greater than zero, accept
            % this iterate temporarily. Adjust shift s so that (s - c)
            % is greater than zero
            if e >= -1 && f < 1
                foundFeasible = true;
            end
            % If exitflag is less than -1 or 'f' is greater than 1,
            % then the next sub-problem is feasible
            if e < -1 || ~foundFeasible
                psAugParam.step = 'Infeasible';
                nextIterate.f = feval(FUN,reshapeinput(X,nextIterate.x),objFcnArg{:});
            else % Update the current Iterate
                nextIterate.x = x(1:end-1);
            end
        end
    end
    
    %-------- Update the iterate only if the step is feasible-----------
    if isempty(psAugParam.step) % The step is either empty or 'infeasible'
        % Compute the constraints at the intermediate point (which should be feasible)
        [nextIterate.cineq(:),nextIterate.ceq(:)] = feval(nonlcon,reshapeinput(X,nextIterate.x),conFcnArg{:});
        % This should not happen (sanity check)
        if numNonlinIneqcstr
            shiftedConstr  =  psAugParam.shift(1:numNonlinIneqcstr) - nextIterate.cineq + options.TolCon;
            if any (shiftedConstr <= 0)
                % Adjust shift
                psAugParam.shift(1:numNonlinIneqcstr) = max(0,nextIterate.cineq) + 1e-4;
                optimState.how = [optimState.how '/Infeasible'];
            end
        end
        % Compute the objective function value
        nextIterate.f = feval(FUN,reshapeinput(X,nextIterate.x),objFcnArg{:});
        % Compute deltaX
        optimState.deltaX = norm(nextIterate.x - Iterate.x);
        optimState.deltaF = abs(nextIterate.f - Iterate.f);
        % Prepare Iterate for next iteration
        Iterate = nextIterate;
    else
        Iterate.f = FVAL;
    end
    % Special case: If solver is stopped by plot/output function in the
    % sub-problem we terminate 
    if innerExitFlag == -1
        optimState.stopOutput = true;
        continue;
    end
    % Call output and plot functions
    if options.OutputTrue || options.PlotTrue
        callOutputPlotFunctions('iter')
    end
end % End of while loop

% Call output/plot functions
if options.OutputTrue || options.PlotTrue
    callOutputPlotFunctions('done')
end

OUTPUT.pollmethod = options.PollMethod;
OUTPUT.searchmethod = options.SearchMethod;
OUTPUT.iterations = optimState.Iter;
OUTPUT.funccount = optimState.FunEval;
OUTPUT.meshsize = psAugParam.currentTolMesh;
OUTPUT.maxconstraint = maxConstr;
OUTPUT.message = msg;
%---------------------------------------------------------------
% Constraint formulation to handle infeasiblity problem
    function [cin,ceq] = infeasConstr(input,oldShift)
        x = input(1:end-1);
        cin = zeros(numNonlinIneqcstr,1);
        [cin(:),ceq] = feval(nonlcon,reshapeinput(X,x),conFcnArg{:});
        cin = cin - input(end)*oldShift;
    end % End of infeasConstr
%---------------------------------------------------------------
% Output function to stop fmincon when objective function is less than one 
    function stop = infeasOutputFcn(unused1,optimValues,flag)
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
%---------------------------------------------------------------
% Nested function to call output/plot functions
    function callOutputPlotFunctions(state)
        optimvalues.x = X; 
        optimvalues.iteration = optimState.Iter; 
        optimvalues.fval = Iterate.f;
        optimvalues.meshsize = psAugParam.currentTolMesh;
        optimvalues.funccount = optimState.FunEval;
        optimvalues.method = optimState.how; 
        optimvalues.TolFun = optimState.deltaF; 
        optimvalues.TolX = optimState.deltaX;
        optimvalues.problemtype = OUTPUT.problemtype;
        optimvalues.nonlinineq = Iterate.cineq; 
        optimvalues.nonlineq = Iterate.ceq;

        solverName = 'Pattern Search';
        switch state
            case {'init','iter'}
                if(options.PlotTrue)
                    optimState.stopPlot = gadsplot(options,optimvalues,state,solverName);
                end                
                if (options.OutputTrue)
                    [optimState.stopOutput,options,optchanged] = psoutput(options.OutputFcns,options.OutputFcnsArg, ...
                        optimvalues,options,state);
                    if optchanged % Check options
                        options = checkoptions(options,defaultopt,numberOfVariables);
                    end
                end
            case 'done'
                if(options.PlotTrue)
                    gadsplot(options,optimvalues,state,solverName);
                end                
                if (options.OutputTrue)
                    psoutput(options.OutputFcns,options.OutputFcnsArg,optimvalues,options,state);
                end
        end
    end % End of callOutputPlotFunctions
%------------------------------------------------------------------
end  % End of pfmincon


