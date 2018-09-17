function [X,FVAL,EXITFLAG,OUTPUT] = pfminlcon(FUN,objFcnArg,initialX,numberOfVariables,Iterate, ...
            Aineq,bineq,Aeq,beq,NullBasisAeq,lb,ub,options,defaultopt,OUTPUT)
%PFMINLCON Finds a linearly constrained minimum of a function.
%   PFMINLCON solves problems of the form:
%           min F(X)    subject to:      A*x <= b
%            X                          Aeq*x = beq
%                                      LB <= X <= UB
%
%   Private to PATTERNSEARCH

%   Copyright 2003-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2012/08/21 00:24:41 $

neqcstr   = size(Aeq,1);

% Get some initial values
[optimState,nextIterate,MeshSize,EXITFLAG,run] = ...
    getinitial(Iterate,numberOfVariables,neqcstr,lb,ub,options);

X = initialX;
X(:) = Iterate.x;
FVAL = Iterate.f;
% Determine who is the caller
callStack = dbstack;
caller = callStack(2).file(1:end-2);

% Call output and plot functions
if options.OutputTrue || options.PlotTrue
    % Set state for plot and output functions (only pfmincon will have
    % 'interrupt' state)
    if ~strcmp(caller,'pfmincon')
        currentState = 'init';
    else
        currentState = 'interrupt';
    end
    callOutputPlotFunctions(currentState)
end

% Setup display header
if  options.Verbosity>1
    fprintf('\n\nIter     f-count          f(x)      MeshSize     Method\n');
end
% Set state for plot and output functions (only pfmincon will have
% 'interrupt' state)
if ~strcmp(caller,'pfmincon')
    currentState = 'iter';
else
    currentState = 'interrupt';
end
while run
    % Check for convergence
    [X,EXITFLAG,FVAL,msg,run] = isconverged(optimState,options,MeshSize, ...
        nextIterate,X,EXITFLAG,run);

    if ~run
        continue;
    end
    % SEARCH.
    [successSearch,nextIterate,optimState] = search(FUN,X,Iterate,MeshSize,Aineq,bineq, ...
            Aeq,beq,NullBasisAeq,lb,ub,OUTPUT.problemtype,objFcnArg,optimState,options);
    % POLL
    if ~successSearch  % Unsuccessful search
        [successPoll,nextIterate,optimState] = poll(FUN,X,Iterate,MeshSize,Aineq,bineq, ...
            Aeq,beq,NullBasisAeq,lb,ub,OUTPUT.problemtype,objFcnArg,optimState,options);
    else
        successPoll =0; % Reset this parameter because this is used to update meshsize
    end

    % Scale the variables in every iterations
    if any(strcmpi(options.ScaleMesh,{'dynamic','on'})) && ~neqcstr
        optimState.scale = logscale(lb,ub,mean(Iterate.x,2));
    end

    % Update
    [MeshSize,Iterate,X,optimState] = updateparam(successPoll,successSearch, ...
        MeshSize,nextIterate,Iterate,X,optimState,options);
    
    % Call output and plot functions
    if options.OutputTrue || options.PlotTrue
        callOutputPlotFunctions(currentState)
    end
end % End of while

% Call output and plot functions
if options.OutputTrue || options.PlotTrue
    % Set state for plot and output functions (only pfmincon will have
    % 'interrupt' state)
    if ~strcmp(caller,'pfmincon')
        currentState = 'done';
    else
        currentState = 'interrupt';
    end
    callOutputPlotFunctions(currentState)
end

% Update values of OUTPUT structure
OUTPUT.pollmethod = options.PollMethod; % This might change via output function
OUTPUT.searchmethod = options.SearchMethod; % This might change via output function
OUTPUT.iterations = optimState.Iter;
OUTPUT.funccount = optimState.FunEval;
OUTPUT.meshsize = MeshSize;
OUTPUT.maxconstraint = norm([Aeq*X(:)-beq; max([Aineq*X(:) - bineq;X(:) - ub(:);lb(:) - X(:)],0)],Inf);
OUTPUT.message = msg;

%-----------------------------------------------------------------
% Nested function to call output/plot functions
    function callOutputPlotFunctions(state)
        optimvalues.x = X; 
        optimvalues.iteration = optimState.Iter;
        optimvalues.fval = Iterate.f;
        optimvalues.problemtype = OUTPUT.problemtype;
        optimvalues.meshsize = MeshSize;
        optimvalues.funccount = optimState.FunEval;
        optimvalues.method = optimState.how;
        optimvalues.TolFun = optimState.deltaF;
        optimvalues.TolX = optimState.deltaX;
        solverName = 'Pattern Search';
        switch state
            case {'init', 'iter'}
                if options.PlotTrue
                    optimState.stopPlot = gadsplot(options,optimvalues,state,solverName);
                end                
                if options.OutputTrue
                    [optimState.stopOutput,options,optchanged] = psoutput(options.OutputFcns,options.OutputFcnsArg, ...
                        optimvalues,options,state);
                    if optchanged % Check options
                        options = checkoptions(options,defaultopt,numberOfVariables);
                    end
                end
            case 'interrupt'
                if options.PlotTrue
                    optimState.stopPlot = gadsplot(options,optimvalues,state,solverName);
                end
                if options.OutputTrue
                    optimState.stopOutput = psoutput(options.OutputFcns,options.OutputFcnsArg, ...
                        optimvalues,options,state);
                end
           case 'done'
                optimvalues.x = X; optimvalues.iteration = optimState.Iter; optimvalues.fval = Iterate.f;
                optimvalues.problemtype = OUTPUT.problemtype; optimvalues.meshsize = MeshSize; optimvalues.funccount = optimState.FunEval;
                optimvalues.method = optimState.how; optimvalues.TolFun = optimState.deltaF; optimvalues.TolX = optimState.deltaX;
                if options.PlotTrue
                    gadsplot(options,optimvalues,state,solverName);
                end
                if options.OutputTrue
                    psoutput(options.OutputFcns,options.OutputFcnsArg,optimvalues,options,state);
                end
        end
    end % End of callOutputPlotFunctions
%------------------------------------------------------------------
end  % End of pfminlcon

