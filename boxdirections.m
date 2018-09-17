function [Basis,TangentCone] = boxdirections(AdaptiveMesh,MeshSize,x,lb,ub,tol)
%BOXDIRECTIONS finds search vectors when bound constraints are present.
%   AdaptiveMesh is a boolean value which is true if MADS is used to poll
%   and X is the current point. Input arguments 'lb' and 'ub' are lower 
%   and upper bounds respectively and 'tol' is the binding tolerance used
%   to determine whether bounds are active or not.

%   Copyright 2003-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2012/08/21 00:23:20 $


vars = length(x);
I = eye(vars);
% Check which bounds are active for lb <= x <= ub at 'x'
active = abs(x - lb) <=tol | abs(x - ub) <=tol;
% Include all directions parallel to active bounds
TangentCone = I(:,active);

% N linearly independent vectors 
if ~AdaptiveMesh
   Basis = I(:,~active);
else
    pollParam = 1/sqrt(MeshSize);
    lowerT = tril((round((pollParam+1)*rand(vars)-0.5)),-1);
    diagtemp = pollParam*sign(rand(vars,1)-0.5);
    diagtemp(~diagtemp) = pollParam*sign(0.5-rand);
    diagT  = diag(diagtemp);
    Basis = lowerT + diagT;
    order = randperm(vars);
    Basis = Basis(order,order);
end


