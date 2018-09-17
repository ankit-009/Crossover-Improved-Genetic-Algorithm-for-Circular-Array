function [psAugParam,numNonlinIneqcstr,numNonlinEqcstr,numNonlinCstr] = psAugInit(Iterate,options)
%PSAUGINIT Initial values of parameters for augmented Lagrangian PS 
%   Private to PATTERNSEARCH

%   Copyright 2005-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2012/08/21 00:24:54 $


% How many nonlinear constraints?
numNonlinIneqcstr = length(Iterate.cineq);
numNonlinEqcstr   = length(Iterate.ceq);
numNonlinCstr = numNonlinIneqcstr + numNonlinEqcstr;
% Initialize shift
psAugParam.shift = zeros(numNonlinIneqcstr,1);
% Initialize Lagrange estimates
psAugParam.lambda = zeros(numNonlinCstr,1);
% Initialize Lagrange estimate updates (to be used in AugUpdate)
psAugParam.lambdabar = zeros(numNonlinCstr,1);
% Parameter used to update currentTolMesh (successful subproblem)
psAugParam.betamesh  = 1.01; % 1.0 0.9 0.25
% Parameter used to update currentTolCon (successful subproblem)
psAugParam.betaconstr  = 0.25; %0.9, 0.25;
% Parameter used to update currentTolMesh (unsuccessful subproblem)
psAugParam.alphamesh = 1.0; % 1.0, 0.75, 0.1
% Parameter used to update currentTolCon (unsuccessful subproblem)
psAugParam.alphaconstr = 0.75; %1.0, 0.75, 0.1;
% Parameter used for constraint satisfaction after a sub-problem is solved
psAugParam.alphaL = min(1,psAugParam.alphaconstr/(1-psAugParam.alphaconstr)); % <= 1
% Constants used for tolerances
psAugParam.startTolCon  = 1.0;  % etas
psAugParam.startOmega   = 1.0;  % omega
% Initial mesh size should not be too small (multiply by this factor)
meshConstant = 1e3;
% Positive constant added to shift
shiftConstant = 1e-4;

% Initial value for penalty parameter
penalty = psoptimget(options,'InitialPenalty',psoptimset,'fast');

% Inverse of penalty (to be used in equations)
invPenalty = 1/penalty;
% Shifts s.t.  s_i - c_i(x) > 0; s_i > 0
psAugParam.shift = max(0,Iterate.cineq) + shiftConstant;
%psAugParam.shift(numNonlinIneqcstr+1:numNonlinCstr) = max(shiftConstant,Iterate.ceq);
%psAugParam.lambda = (psAugParam.shift./invPenalty).^(1/psAugParam.alphaL);

psAugParam.lambda(1:numNonlinIneqcstr) = (psAugParam.shift./invPenalty).^(1/psAugParam.alphaL);
% Initial value of lambda for equality constraints (for backward compatibility reasons)
temp = max(shiftConstant,Iterate.ceq); 
psAugParam.lambda(numNonlinIneqcstr+1:numNonlinCstr) = (temp./invPenalty).^(1/psAugParam.alphaL);


% Parameters to be changed in outer loop (main problem)
psAugParam.currentTolCon  = min(1e-1,psAugParam.startTolCon* ...
    (invPenalty^psAugParam.alphaconstr));
psAugParam.currentOmega   = min(1e-1,psAugParam.startOmega* ...
    (invPenalty^psAugParam.alphamesh));
psAugParam.currentTolMesh = max(options.TolMesh,meshConstant* ...
    psAugParam.currentOmega/(1 + norm(psAugParam.lambda) + penalty));
psAugParam.currentTolMesh = min(psAugParam.currentTolMesh,options.InitialMeshSize);
psAugParam.startTolMesh = psAugParam.currentTolMesh;
psAugParam.step = '';
psAugParam.penalty = penalty;
