function [gaAugParam,state] = gaAugUpdate(Iterate,gaAugParam,state,options,numNonlinIneqcstr,numNonlinCstr)
%GAAUGUPDATE Updates values of parameters for augmented Lagrangian GA 
%   Private to GA

%   Copyright 2005-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2012/08/21 00:23:44 $


state.how = ' ';
 state.Generation = state.Generation + 1;
stallTol = min(options.TolFun,eps);

lambda = gaAugParam.lambda;
penalty = gaAugParam.penalty;
shift = gaAugParam.shift;

% Initialize new values of lambda
lambdabar = zeros(numNonlinCstr,1);
C1 = 0;

% Calculate switching condition which will help us decide whether to
% update multipliers or increase penalty. We also update lambda.

if numNonlinIneqcstr % Nonlinear inequality constraints
    shiftedConstr = -Iterate.cineq + shift(1:numNonlinIneqcstr);
    lambdabar(1:numNonlinIneqcstr) = (lambda(1:numNonlinIneqcstr).*shift(1:numNonlinIneqcstr))./shiftedConstr;
    if numNonlinCstr > numNonlinIneqcstr % Equalities too
    lambdabar(numNonlinIneqcstr+1:numNonlinCstr) = ...
        lambda(numNonlinIneqcstr+1:numNonlinCstr) + Iterate.ceq.*penalty;
    end
    for i = 1:numNonlinIneqcstr
        if lambda(i) > eps
          C1 =  C1 + (Iterate.cineq.*lambdabar(i))./(lambda(i).^gaAugParam.alphaL);
        end
    end
    C1 = norm(C1);
    C2 = norm(Iterate.ceq);
else % Only nonlinear equality constraints
    lambdabar(numNonlinIneqcstr+1:numNonlinCstr) = ...
        lambda(numNonlinIneqcstr+1:numNonlinCstr) + Iterate.ceq.*penalty;
C2 = norm(Iterate.ceq);
end
% Update multipliers?
updateLang = max(C1 , C2) <= gaAugParam.currentTolCon;
if updateLang % Update Lagrange multipliers estimate
    state.how = 'Update multipliers';
    invPenalty = max(stallTol,1/penalty);
    % Update these three quantities.
    gaAugParam.lambda = lambdabar;
    gaAugParam.currentOmega = min(1e-1,gaAugParam.currentOmega*(invPenalty^gaAugParam.betafun));
    gaAugParam.currentTolFun = gaAugParam.currentOmega;
    gaAugParam.currentTolCon  = min(1e-1,gaAugParam.currentTolCon*(invPenalty^gaAugParam.betaconstr));
else % Increase penalty
    state.how = 'Increase penalty';
    gaAugParam.penalty = options.PenaltyFactor*penalty;
    invPenalty = max(stallTol,1/gaAugParam.penalty);
    gaAugParam.currentOmega = min(1e-1,gaAugParam.startOmega*(invPenalty^gaAugParam.alphafun));
    gaAugParam.currentTolFun = gaAugParam.currentOmega;
    gaAugParam.currentTolCon = min(1e-1,gaAugParam.startTolCon*(invPenalty^gaAugParam.alphaconstr));
end
% If cuurentTolFun is too low then use penalty values to limit accuracy
% of inner sub-problem's solution and constraint tolerance
if gaAugParam.currentTolFun <= options.TolFun
    gaAugParam.currentTolFun = max(gaAugParam.currentOmega,gaAugParam.currentTolFun);
end
if gaAugParam.currentTolCon <= options.TolCon
    gaAugParam.currentTolCon = max(options.TolCon,gaAugParam.currentTolCon);
end
 
% Compute shift
gaAugParam.shift = invPenalty*(gaAugParam.lambda(1:numNonlinIneqcstr).^gaAugParam.alphaL);
gaAugParam.lambdabar = lambdabar;
