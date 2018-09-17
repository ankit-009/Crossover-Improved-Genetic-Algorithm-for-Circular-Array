function [Basis,TangentCone] = lcondirections(AdaptiveMesh,MeshSize,x,Aineq,bineq,Aeq,lb,ub,tol)
%LCONDIRECTIONS finds search vectors when linear constraints and bounds are
%   present. AdaptiveMesh is a boolean value which is true if MADS is used to poll
%   and X is the current point. Input arguments 'Aineq','bineq','Aeq','lb', and
%   'ub' are linear and bound constraint information and 'tol' is the binding tolerance
%   used for determining whether constraints are active or not. Note that 'beq' is not 
%   used to calculate tangent cone.

%   Copyright 2003-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2012/08/21 00:24:27 $



% Initialization
TangentCone = [];
vars = length(x);
% N linearly independent vectors
if ~AdaptiveMesh
    Basis = eye(vars);
else
    % Create random generating directions
    pollParam = 1/sqrt(MeshSize);
    lowerT = tril((round((pollParam+1)*rand(vars)-0.5)),-1);
    diagtemp = pollParam*sign(rand(vars,1)-0.5);
    diagtemp(~diagtemp) = pollParam*sign(0.5-rand);
    diagT  = diag(diagtemp);
    Basis = lowerT + diagT;
    order = randperm(vars);
    Basis = Basis(order,order);
end

Normals = zeros(1,vars);
tolDep = 100*vars*eps;
% The normal cone generators for minimum epsilon is in active set
% (Lewis & Torczon section 8.2)
if ~isempty(Aineq)
    while rank(Normals) ~= min(size(Normals))
        if tol < tolDep
            error(message('globaloptim:lcondirections:degenconstr'));
        end
        active = abs(Aineq*x - bineq) <= tol;
        Normals = Aineq(active,:);
        tol = tol/2;
    end
end
% Add equality constraints to the normal cone
Normals = [Aeq; Normals]';
% Lewis & Torczon section 8.2. T = V*inv(V'V), which is computed using QR
% decomposition
if ~isempty(Normals)
     Normals(:,(~any(Normals))) = [];
    [Q,R] = qr(Normals,0);
    TangentCone = Q/R';
    Basis = Basis - TangentCone*Normals';
end

% Add active bounds in tangent cone as well
if any(isfinite(ub) | isfinite(lb))
    I = eye(vars);
    % Check which bounds are active for lb <= x <= ub at 'x'
    active = abs(x - lb) <=tol | abs(x - ub) <=tol;
    % Include all directions parallel to active bounds
    TangentCone = [TangentCone, I(:,active)];
end
