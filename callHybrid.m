function [x,fval,funccount,theMessage,conviol] = callHybrid(hfunc,FitnessHybridFcn,x0,HybridFcnArgs,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,goal,weight)
%callHybrid Calls hybrid function used in GA and GAMULTIOBJ
%
%   This function is private to GA and GAMULTIOBJ

%   Copyright 2007-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2012/08/21 00:23:21 $

if nargin < 13, weight = [];
    if nargin < 12, goal = [];
        if nargin < 11,  nonlcon = [];
            if nargin < 10, ub = [];
                if nargin < 9, lb = [];
                    if nargin < 8, beq = [];
                        if nargin < 7, Aeq = [];
                            if nargin < 6, bineq = [];
                                if nargin < 5, Aineq = [];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% We want to set FunValCheck to 'on'  and Display to 'off' by default
if isempty(HybridFcnArgs)
    HybridFcnArgs = {struct('Display','off','FunValCheck','on')};
end

[lastmsg, lastid] = lastwarn;
warnstate = warning;
warning off;

try
    % Determine which solver to call
    switch hfunc
        case 'fminsearch'
            [x,fval,~,output] = fminsearch(FitnessHybridFcn,x0,HybridFcnArgs{:});
            funccount = output.funcCount;
            theMessage   = sprintf('FMINSEARCH: %s', output.message);
            % fminsearch only applies to unconstrained problems
            conviol = [];
        case 'patternsearch'
            [x,fval,~,output] = patternsearch(FitnessHybridFcn,x0,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,HybridFcnArgs{:});
            funccount =  output.funccount;
            theMessage   = sprintf('PATTERNSEARCH: %s', output.message);
            if isfield(output, 'maxconstraint')
                conviol = output.maxconstraint;
            else
                conviol = [];
            end
        case 'fminunc'
            [x,fval,~,output] = fminunc(FitnessHybridFcn,x0,HybridFcnArgs{:});
            funccount = output.funcCount;
            theMessage   = sprintf('FMINUNC: %s', output.message);
            % fminunc only applies to unconstrained problems
            conviol = [];            
        case 'fmincon'
            [x,fval,~,output] = fmincon(FitnessHybridFcn,x0,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,HybridFcnArgs{:});
            funccount =  output.funcCount;
            theMessage   = sprintf('FMINCON: %s', output.message);
            conviol = output.constrviolation;
        case 'fgoalattain'
            [x,fval,~,~,output] = fgoalattain(FitnessHybridFcn,x0,goal,weight,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,HybridFcnArgs{:});
            funccount =  output.funcCount;
            theMessage   = sprintf('FGOALATTAIN: %s', output.message);
            conviol = output.constrviolation;
    end   
    warning(warnstate);
    lastwarn(lastmsg,lastid);
catch optimFcn_ME % Restore warnings in case optim functions error out.
    warning(warnstate);
    lastwarn(lastmsg,lastid);
    rethrow(optimFcn_ME);
end
