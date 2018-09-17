function state  = gamultiobjHybrid(FitnessFcn,x,fval,state,Aineq,bineq,Aeq,beq,lb,ub,options)
%gamultiobjHybrid setup arguments for the hybrid function and run it
%
%   This function is private to GAMULTIOBJ

%   Copyright 2007-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.6 $  $Date: 2012/08/21 00:24:07 $

% Who is the hybrid function
if isa(options.HybridFcn,'function_handle')
    hfunc = func2str(options.HybridFcn);
else
    hfunc = options.HybridFcn;
end

% Create anonymous function to be passed to the hybrid function
FitnessHybridFcn = @(x) FitnessFcn(x,options.FitnessFcnArgs{:});
% Inform about hybrid scheme
if  options.Verbosity > 1
    fprintf('%s%s%s\n','Switching to the hybrid optimization solver (',upper(hfunc),').');
end
if ~strcmpi(hfunc,'fgoalattain')
    warning(message('globaloptim:gamultiobjHybrid:notMultiObjHybridFcn', ... 
        upper(hfunc),'FGOALATTAIN'));
    hfunc = 'fgoalattain';
end
% Create room for new population from the hybrid solver
xHybrid = zeros(size(x));
fHybrid = zeros(size(fval));
funccount = zeros(1,size(x,1));
% Local variable to use inside parfor loop
args = options.HybridFcnArgs;
% Use for if SerialUserFcn is 'true'; the if-and-else part should
% be exactly same other than parfor-for syntax difference.
if options.SerialUserFcn
    for i = 1:size(x,1)
        % Do not call hybrid function if fval(i,:) have Inf or NaN
        if any(~isfinite(fval(i,:)))
            acceptHybridSolution = false;
        else
            [pweight,pgoal] = pseudoWeightAndGoal(fval(i,:),fval);
            try
                [xHybrid(i,:),fHybrid(i,:),funccount(i)] = callHybrid(hfunc,FitnessHybridFcn,x(i,:),args,Aineq,bineq,Aeq,beq,lb,ub,[],pgoal,pweight);
                % For fgoalattain, we need to pass the final point and the
                % problem constraints to isHybridPointFeasible.
                hybridConInfo = {xHybrid(i,:)',Aineq,bineq,Aeq,beq,lb,ub};
                hybridPointFeasible = isHybridPointFeasible(hybridConInfo, 'fgoalattain', args{:}); 
                if hybridPointFeasible
                    acceptHybridSolution = true;
                else
                    acceptHybridSolution = false;
                end
            catch optim_ME
                acceptHybridSolution = false;
                if ~isempty(findstr(optim_ME.identifier,'optimlib:optimfcnchk:checkfun'))
                    % do nothing
                else
                    rethrow(optim_ME);
                end
            end
        end
        if ~acceptHybridSolution
            xHybrid(i,:) = x(i,:); % Do not change x
            fHybrid(i,:) = fval(i,:);
        end
    end
else % Run fgoalattain in parallel using parfor
    parfor (i = 1:size(x,1))
        % Do not call hybrid function if fval(i,:) have Inf or NaN
        if any(~isfinite(fval(i,:)))
            acceptHybridSolution = false;
        else
            [pweight,pgoal] = pseudoWeightAndGoal(fval(i,:),fval);
            acceptHybridSolution = false;
            try
                [xHybrid(i,:),fHybrid(i,:),funccount(i)] = callHybrid(hfunc,FitnessHybridFcn,x(i,:),args,Aineq,bineq,Aeq,beq,lb,ub,[],pgoal,pweight);
                % For fgoalattain, we need to pass the final point and the
                % problem constraints to isHybridPointFeasible.
                hybridConInfo = {xHybrid(i,:)',Aineq,bineq,Aeq,beq,lb,ub};
                hybridPointFeasible = isHybridPointFeasible(hybridConInfo, 'fgoalattain', args{:}); 
                if hybridPointFeasible
                    acceptHybridSolution = true;
                end
            catch optim_ME
                if ~isempty(findstr(optim_ME.identifier,'optimlib:optimfcnchk:checkfun'))
                    % do nothing
                else
                    rethrow(optim_ME);
                end
            end
        end
        if ~acceptHybridSolution
            xHybrid(i,:) = x(i,:); % Do not change x
            fHybrid(i,:) = fval(i,:);
        end
    end
end
state.FunEval = state.FunEval + sum(funccount);

% Merge xHybrid with current population and find a new non-dominated front
% of size(x,1)
[state.Population,state.Score,state.Rank,state.Distance]  = ...
    rankAndDistance([state.Population; xHybrid],[state.Score; fHybrid],options,size(state.Score,1));


% Inform about hybrid scheme termination
if  options.Verbosity > 1
    fprintf('%s %s\n',upper(hfunc), 'terminated.');
end
