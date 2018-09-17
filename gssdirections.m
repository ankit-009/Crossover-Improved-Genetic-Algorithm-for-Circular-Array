function [Basis,TangentCone,AugmentedDirs] = gssdirections(MeshSize,XCurrent,Aineq,bineq,lb,ub,NullBasisAeq,epsBall)
%GSSDIRECTIONS finds search vectors when linear constraints and bounds are
% present by explicitly calculating the Tangent Cone (of inequality 
% constraints) in the Null Space of the linear equality constraints. 
% 'MeshSize' is the current step size and 'XCurrent' is the current point. 
% 'Aineq','bineq','Aeq','lb', and 'ub' are linear and bound constraint 
%  information. 
% 'NullBasisAeq' is a set of basis vectors that spans the Null Space of the
%  linear equality constraints in 'Aeq'.
% 'epsBall' is the binding tolerance used for determining whether constraints
%  are active or not.  
% 'Basis' is the set of core search directions spanning the unconstrained space
% 'TangentCone' is the set of core search directions spanning the 
%  space of linear inequality constraints
% 'AugmentedDirs' contains additional search directions besides the ones
% contained in 'Basis' and 'TangentCone'. These are not required to ensure 
% convergence but might speed up the search process.
%
% Private to POLL
%
% See also LCONDIRECTIONS
%   
% Copyright 2009 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2012/08/21 00:24:16 $

% Initialization
nVars = length(XCurrent);
Normals = zeros(1);
tolDep = 100*nVars*eps;

% Find the active inequality constraints at the given iterate 
% in an epsilon ball of starting size = tolBind.
% If resulting set is dependent, reduce epsilon ball by a factor of 
% 0.5 each time till either an independ set is found. If independent set
%  is not found, error out.
if ~isempty(Aineq)
  while rank(Normals) ~= min(size(Normals))
    if epsBall < tolDep
      % This direction generator will not work for dependent active constraints 
      error(message('globaloptim:gssdirections:degenconstr'));
    end
    active = abs(Aineq*XCurrent - bineq) <= epsBall;
    Normals = Aineq(active,:);
    epsBall = epsBall/2;
  end
end

if (isempty(Normals) || isscalar(Normals))
  % If there are inequalities but no active ones, Normals is empty.
  % This is also the code path taken when only equality constraints are present
  % In the absence of active inequality constraints  no Augmented Directions 
  % and TangentCone is empty too.
  Basis = NullBasisAeq;
  TangentCone = zeros(nVars,0);
  AugmentedDirs = zeros(nVars,0); 
else
  ActiveIneq = Normals;
  % Project active inequality constraints to Null space of  equality constraints
  ProjIneq = (ActiveIneq*NullBasisAeq)';
  
  % Error if working set of constraints is degenerate
  if rank(ProjIneq) ~= min(size(ProjIneq))
    % This direction generator will not work for dependent active constraints 
    error(message('globaloptim:gssdirections:degenerateWSet'));
  end
  
  % First set of generators : Null space of ProjIneq'
  % space perpendicular to Col space of ProjIneq
  [numRows numCols] = size(ProjIneq);
  [QIneq,RIneq] = qr(ProjIneq);
  
  if (numRows > numCols)
    Basis = NullBasisAeq*QIneq(:,numCols+1:end);
  else
    Basis = zeros(nVars,0);
  end
  
  % Second set of generators, Tangent Cone : Right inverse of ProjIneq'
  % see  Lewis et al, Implementing Generating Set Search Methods 
  % For Linearly Constrained Minimization, SIAM Journal on Scientific 
  % Computing, 29(2007).
  rangeDim = min(numCols,size(QIneq,2));
  QIneq = QIneq(:,1:rangeDim);
  RIneq = RIneq(1:rangeDim,1:rangeDim);
  RightInv = QIneq/RIneq';
  
  % Project back to original space, positive spanning generators
  TangentCone = -NullBasisAeq*RightInv;
  
  % NormalCone of linear constraints, Constraint vectors are positive 
  % spanning set of Normal Cone. 
  % Effectively, AugmentedDirs = NullBasisAeq* NullBasisAeq'*ActiveIneq' 
  % These are not required for convergence but only here to supplement core directions
  % in 'Basis' and 'TangentCone'
  AugmentedDirs = NullBasisAeq*ProjIneq;
end

% Note  bounds are handled separately from inequality constraints because its 
% trivial to find the tangent cone of bounds.
if any(isfinite(lb) | isfinite(ub))
  I = eye(nVars);
  % Check which bounds are active for lb <= XCurrent <= ub at 'XCurrent'
  activeLB = abs(XCurrent - lb) <= epsBall;
  activeUB = abs(XCurrent - ub) <= epsBall;
  % Include all directions parallel to active bounds
  TangentCone = [TangentCone I(:,activeLB) -I(:,activeUB)];
end
