function [LinearConstr,Iterate,nineqcstr,neqcstr,ncstr] = constrValidate(nonlcon,Iterate,Aineq,bineq,Aeq,beq,lb,ub,intcon,type,options)
%constrValidate validate nonlinear constraint function and create LinearConstr
%   structure for linear constraints
%
%   Private to GA

%   Copyright 2005-2012 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2012/08/21 00:23:28 $

% Number of constraints
nineqcstr = size(Aineq,1);
neqcstr   = size(Aeq,1);
ncstr     = nineqcstr + neqcstr;

% Validate nonlinear constraint function. When integer constraints are
% present, nonlinear constraints are validated in gaminlpintconstrvalidate.
if ~isempty(nonlcon) && isempty(intcon)
    % Evaluate nonlinear constraints for the first time
    try
        if strcmpi(options.Vectorized,'on')
            % Evaluate the nonlinear constraints at two points, [Iterate.x Iterate.x]'.
            [tmpCineq,tmpCeq] = nonlcon([Iterate.x Iterate.x]');
        else
            [cineq,ceq] = nonlcon(Iterate.x');
        end
    catch userFcn_ME
        gads_ME = MException('globaloptim:constrvalidate:confunCheck', ...
            'Failure in initial user-supplied nonlinear constraint function evaluation.');
            userFcn_ME = addCause(userFcn_ME,gads_ME);
            rethrow(userFcn_ME)
    end
    if strcmpi(options.Vectorized,'on')
        % Check whether row dimension of the constraint matrix is
        % two, corresponding to two points.
        % We don't perform any additional checks for the two-constraint case.
        if ( ~isempty(tmpCineq) && ( size(tmpCineq,1) ~= 2 ) ) || ...
                ( ~isempty(tmpCeq) && ( size(tmpCeq,1) ~= 2 ) )
            error(message('globaloptim:constrvalidate:confunCheckVectorized', ...
                'Vectorized','on'));
            
        end
        
        cineq = [];
        ceq = [];        
        if ~isempty(tmpCineq)
            cineq = tmpCineq(1,:);
        end
        if ~isempty(tmpCeq)
            ceq = tmpCeq(1,:);
        end
        
    end
    
    Iterate.cineq = zeros(numel(cineq),1);
    Iterate.ceq = zeros(numel(ceq),1);
    Iterate.cineq(:) = cineq;
    Iterate.ceq(:) = ceq;
    c = [Iterate.cineq;Iterate.ceq];
    if ~all(isreal(c) & isfinite(c))
        error(message('globaloptim:constrvalidate:confunNotReal'));
    end
end

if ~strcmpi(type,'unconstrained')
    % Check linear constraint satisfaction at stopping?
    linconCheck = false;
    LinearConstr = struct('Aineq',Aineq,'bineq',bineq,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub, ...
    'feasibleX',Iterate.x,'type',type);
    % Mutation function for constrained GA
    mutationFcn = func2str(options.MutationFcn);
    if ~strcmpi(mutationFcn,'mutationadaptfeasible')
        linconCheck = true;
        % Warn for using 'mutationuniform' or 'mutationgaussian' mutation
        % functions. Do not warn if we're solving a mixed integer problem,
        % as we will set the mutation operator.
        if any(strcmpi(mutationFcn,{'mutationuniform','mutationgaussian'})) && isempty(intcon)
           warning(message('globaloptim:constrvalidate:unconstrainedMutationFcn', mutationFcn))
        end
    end
    LinearConstr.linconCheck = linconCheck;
else
    LinearConstr = struct('Aineq',Aineq,'bineq',bineq,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub, ...
    'feasibleX',Iterate.x,'type',type,'linconCheck',false);
end
