function [psAugParam,optimState] = psAugUpdate(Iterate,psAugParam,optimState,options, ...
        numNonlinIneqcstr,numNonlinCstr)
%PSAUGUPDATE Updates values of parameters for augmented Lagrangian PS 
%   Private to PATTERNSEARCH

%   Copyright 2005-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2012/08/21 00:24:56 $

optimState.how = ' ';
optimState.Iter = optimState.Iter + 1;
stallTol = min(options.TolMesh,eps);

lambda = psAugParam.lambda;
penalty = psAugParam.penalty;
shift = psAugParam.shift;

% Initialize new values of lambda
lambdabar = zeros(numNonlinCstr,1);
C1 = 0;
% Calculate switching condition which will decide whether to
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
          C1 =  C1 + (Iterate.cineq.*lambdabar(i))./(lambda(i).^psAugParam.alphaL);
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
updateLang = max(C1 , C2) <= psAugParam.currentTolCon;
if updateLang % Update Lagrange multipliers estimate
    optimState.how = 'Update multipliers';
    invPenalty = max(stallTol,1/penalty);
    % Update these three quantities.
    psAugParam.lambda = lambdabar;
    psAugParam.currentOmega = min(1e-1,psAugParam.currentOmega*(invPenalty^psAugParam.betamesh));
    psAugParam.currentTolMesh = psAugParam.currentOmega;
    psAugParam.currentTolCon  = min(1e-1,psAugParam.currentTolCon*(invPenalty^psAugParam.betaconstr));
else % Increase penalty
    optimState.how = 'Increase penalty';
    psAugParam.penalty = options.PenaltyFactor*penalty;
    invPenalty = max(stallTol,1/psAugParam.penalty);
    psAugParam.currentOmega = min(1e-1,psAugParam.startOmega*(invPenalty^psAugParam.alphamesh));
    psAugParam.currentTolMesh = psAugParam.currentOmega; % Lewis & Torczon
    psAugParam.currentTolCon = min(1e-1,psAugParam.startTolCon*(invPenalty^psAugParam.alphaconstr));
end
% If currentTolCon is too low then use penalty values to
% limit accuracy of inner sub-problem's solution and constraint tolerance
if psAugParam.currentTolCon <= options.TolCon
    psAugParam.currentTolCon = max(options.TolCon,psAugParam.currentTolCon);
end

% Compute shift
psAugParam.shift = invPenalty*(psAugParam.lambda(1:numNonlinIneqcstr).^psAugParam.alphaL);
psAugParam.lambdabar = lambdabar;
