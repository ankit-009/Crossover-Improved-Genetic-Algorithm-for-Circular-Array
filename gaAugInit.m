function [gaAugParam,numNonlinIneqcstr,numNonlinEqcstr,numNonlinCstr] = gaAugInit(Iterate,options)
%GAAUGINIT Initial values of parameters for augmented Lagrangian GA 
%   Private to GA

%   Copyright 2005-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2012/08/21 00:23:42 $


% How many nonlinear constraints?
numNonlinIneqcstr = length(Iterate.cineq);
numNonlinEqcstr   = length(Iterate.ceq);
numNonlinCstr = numNonlinIneqcstr + numNonlinEqcstr;
% Initialize shift 
gaAugParam.shift = zeros(numNonlinIneqcstr,1);
% Initialize Lagrange estimates
gaAugParam.lambda = zeros(numNonlinCstr,1);
% Initialize Lagrange estimate updates (to be used in AugUpdate)
gaAugParam.lambdabar = zeros(numNonlinCstr,1);
% Parameter used to update currentTolFun (successful subproblem)
gaAugParam.betafun  = 1.01; % 1.0 0.9 0.25
% Parameter used to update currentTolCon (successful subproblem)
gaAugParam.betaconstr  = 0.25; %0.9, 0.25;
% Parameter used to update currentTolFun (unsuccessful subproblem)
gaAugParam.alphafun = 1.0; % 1.0 0.75 0.1
% Parameter used to update currentTolCon (unsuccessful subproblem)
gaAugParam.alphaconstr = 0.75; %0.1, 0.75;
% Used to check constraint satisfaction after a sub-problem is solved
gaAugParam.alphaL = min(1,gaAugParam.alphaconstr/(1-gaAugParam.alphaconstr)); % <= 1
% Initial value of tolerances
gaAugParam.startTolCon  = 1.0;  % etas
gaAugParam.startOmega   = 1.0;  % omega
% Initial tolerance on function value should not be too small (multiply by this factor)
funtolConstant = 1e3;
% Positive constant added to shift
shiftConstant = 1e-4;

% Initial values for changing parameters
% Recommended values are 1, 5, and 10  < 1 (rho)
penalty = gaoptimget(options,'InitialPenalty',gaoptimset,'fast');

% Inverse of penalty (to be used in equations)
invPenalty = 1/penalty;
% Shifts s.t. c_i(x) - s_i <0; s_i > 0
gaAugParam.shift = max(0,Iterate.cineq) + shiftConstant;
%gaAugParam.shift(numNonlinIneqcstr+1:numNonlinCstr) = max(shiftConstant,Iterate.ceq); 
gaAugParam.lambda(1:numNonlinIneqcstr) = (gaAugParam.shift./invPenalty).^(1/gaAugParam.alphaL);
% Initial value of lambda for equality constraints (for backward compatibility reasons)
temp = max(shiftConstant,Iterate.ceq); 
gaAugParam.lambda(numNonlinIneqcstr+1:numNonlinCstr) = (temp./invPenalty).^(1/gaAugParam.alphaL);
% Parameters to be changed in outer loop (main problem)
gaAugParam.currentTolCon  = gaAugParam.startTolCon*(invPenalty^gaAugParam.alphaconstr);
gaAugParam.currentOmega   = gaAugParam.startOmega*(invPenalty^gaAugParam.alphafun);
% Lower bound on tolerance for first generation
gaAugParam.currentTolFun =  max(options.TolFun,gaAugParam.currentOmega*funtolConstant);
% Upper bound on tolerance for first generation
gaAugParam.currentTolFun = min(1,gaAugParam.currentTolFun); 
gaAugParam.startTolFun = gaAugParam.currentTolFun;

gaAugParam.step = '';
gaAugParam.penalty = penalty;
