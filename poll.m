function  [msg,nextIterate,optimState] = poll(ObjFunc,Xin,Iterate, ...
    MeshSize,Aineq,bineq,Aeq,beq,NullBasisAeq,lb,ub,problemtype,objFcnArg,optimState,options)
%POLL Performs the poll step in GPS
%
%  OBJFUNC: The objective function on which POLL step is implemented.
%
%  POLLMETHOD: Search directions in POLL step is obtained according to
%  different pollmethod.
%
%  COMPLETEPOLL: If 'off' indicates that POLL can be called off as soon
%  as a better point is found i.e. no sufficient decrease condition is imposed;
%  If 'on' then ALL the points are evaluated and point with least function value
%  is returned. Default is 'off'. If function is expensive, make this 'off'.
%
%  POLLORDER: Ordering of poll directions.
%
%  ITERATE: Incumbent point around which polling will be done. Iterate Stores
%  the current point 'x' and function value 'f' at this point.
%
%  SUCCESSDIR: Last successful POLL/SEARCH direction. This information can be used
%  by the POLL step in ordering the search direction (last successful
%  direction is polled first).
%
%  MESHSIZE: Current mesh size used in POLL step.
%
%  SCALE: Scale factor used to scale the design points.
%
%  TOL: Tolerance used for determining whether constraints are active or not.
%
%  PROBLEMTYPE: This flag is passed to the SEARCH routines, indicating that the
%  problem is 'unconstrained', 'boundconstraints', 'linearconstraints'.
%
%  NEXTITERATE: Successful iterate after polling is done. If POLL is NOT
%  successful, NEXTITERATE is same as ITERATE.
%
%  MSG: A binary flag indicating, POLL is successful or not.

%   Copyright 2003-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2012/08/21 00:24:43 $


constr = true;
tol = options.TolBind;
Successdir = optimState.Successdir;
scale = optimState.scale;
TangentCone = [];
AugmentedDirs = [];
AdaptiveMesh = any(strcmpi(options.PollMethod,{'madspositivebasisnp1','madspositivebasis2n'}));

% Get the directons which forms the positive basis(minimal or maximal basis)
if  strcmpi(problemtype, 'unconstrained')
  DirVector = uncondirections(AdaptiveMesh,MeshSize,Iterate.x);
  constr = false;
elseif strcmpi(problemtype,'boundconstraints')    % Only Box constraints
  [DirVector,TangentCone] = boxdirections(AdaptiveMesh,MeshSize,Iterate.x,lb,ub,tol);
elseif any(strcmpi(problemtype, {'linearconstraints','nonlinearconstr'}))  % Constraints
  if any(strcmpi(options.PollMethod,{'gsspositivebasisnp1','gsspositivebasis2n'}))
    [DirVector,TangentCone,AugmentedDirs] = gssdirections(MeshSize,Iterate.x, ...
      Aineq,bineq,lb,ub,NullBasisAeq,tol);
  else
    [DirVector,TangentCone] = lcondirections(AdaptiveMesh,MeshSize,Iterate.x, ...
      Aineq,bineq,Aeq,lb,ub,tol);
  end
end

% If the point is on the constraint boundary (nonempty TangentCone)
% we use scale = 1
if ~isempty(TangentCone)
  scale = 1;
  TangentCone(:,(~any(TangentCone))) = [];
end

nDirTan = size(TangentCone,2);
% Form directions that positively span the tangent cone at x
switch lower(options.PollMethod)
  case {'positivebasisnp1','gpspositivebasisnp1','madspositivebasisnp1'}
    % Add n+1st vector to the Basis
    DirVector = [-sum(DirVector,2) DirVector];
    nDirVector = size(DirVector,2);
    % Add tangent cone to the direction vectors
    DirVector = [DirVector TangentCone];
    % Total number of search directions
    nDirTotal = nDirTan + nDirVector;
    % Make the index vector to be used to access directions
    indexVec = [1:nDirVector (nDirVector+1):nDirTotal (nDirVector+1):nDirTotal];
    % Vector to take care of sign of directions
    dirSign = [ones(1,nDirVector) ones(1,nDirTan) -ones(1,nDirTan)];
  case {'positivebasis2n','gpspositivebasis2n','madspositivebasis2n'}
    nDirVector = size(DirVector,2);
    % Add tangent cone to the direction vectors
    DirVector = [DirVector TangentCone];
    % Total number of search directions
    nDirTotal = nDirTan + nDirVector;
    % Make the index vector to be used to access directions
    indexVec = [1:nDirVector 1:nDirVector (nDirVector+1):nDirTotal (nDirVector+1):nDirTotal];
    % Vector to take care of sign of directions
    dirSign = [ones(1,nDirVector) -ones(1,nDirVector) ones(1,nDirTan) -ones(1,nDirTan)];
  case {'gsspositivebasisnp1'}
    % Handle the empty DirVector case below in a special if loop
    if isempty(DirVector)
      nDirVector = 0;
    else
      % Add n+1st vector to the Basis
      DirVector = [-sum(DirVector,2) DirVector];
      nDirVector = size(DirVector,2);
    end
    % Add tangent cone to the direction vectors
    DirVector = [DirVector TangentCone];
    % Total number of search directions
    nDirTotal = nDirTan + nDirVector;
    % Make the index vector to be used to access directions
    indexVec = [1:nDirVector (nDirVector+1):nDirTotal];
    % Vector to take care of sign of directions
    dirSign = [ones(1,nDirVector) ones(1,nDirTan)];
  case {'gsspositivebasis2n'}
    nDirVector = size(DirVector,2);
    % Add tangent cone to the direction vectors
    DirVector = [DirVector TangentCone];
    % Total number of search directions
    nDirTotal = nDirTan + nDirVector;
    % Make the index vector to be used to access directions
    indexVec = [1:nDirVector 1:nDirVector (nDirVector+1):nDirTotal];
    % Vector to take care of sign of directions
    dirSign = [ones(1,nDirVector) -ones(1,nDirVector) ones(1,nDirTan)];
  otherwise
    error(message('globaloptim:poll:pollmethod'));
end

% Get order of direction vector based on polling order choice. (Mark Abramson)
pollorder = options.PollingOrder;
OrderVec = 1:length(indexVec);

% MADS algorithms are not affected by poll order
if AdaptiveMesh
  % Do nothing
elseif strcmpi(pollorder,'consecutive')
  % Do nothing
elseif strcmpi(pollorder, 'success')
  % If using 'success' then it is assumed that CompletePoll is 'off' so
  % duplicating search direction (if any) is okay.
  if ~isempty(Successdir)
    DirVector = [DirVector Successdir];
    indexVec  = [nDirTotal+1 indexVec];
    dirSign  =  [1 dirSign];
    OrderVec = 1:length(indexVec);
  end
elseif strcmpi(pollorder, 'random')
  OrderVec = randperm(length(indexVec));
else
  warning(message('globaloptim:poll:invalidPollOrder'));
end
% Total number of trial points
numberOfXtrials = length(OrderVec);

%Get the trial points along the direction vectors;
sites = struct('x',cell(numberOfXtrials,1,1),'f',cell(numberOfXtrials,1,1));
for k = 1:numberOfXtrials
  direction = dirSign(k).*DirVector(:,indexVec(OrderVec(k)));
  sites(k).x = Iterate.x + MeshSize*scale.*direction;
end

%Find an iterate with lower objective
[msg,nextIterate,dirNum,optimState.FunEval] = ...
  psnextfeasible(ObjFunc,Xin,sites,Iterate,Aineq,bineq,Aeq,beq,lb,ub,tol,constr,objFcnArg,optimState.FunEval,options);
if (msg)
  Successdir = dirSign(indexVec(OrderVec(dirNum))).*DirVector(:,indexVec(OrderVec(dirNum)));
  if  AdaptiveMesh
    Iterate = nextIterate;
    % Search again along the successful direction
    site.x = Iterate.x + min(1,4*MeshSize)*scale.*Successdir;
    % Find an iterate with lower objective
    [unused1,nextIterate,unused2,optimState.FunEval] =  ...
      psnextfeasible(ObjFunc,Xin,sites,Iterate,Aineq,bineq,Aeq,beq,lb,ub,tol,constr,objFcnArg,optimState.FunEval,options);
  end
elseif ~isempty(AugmentedDirs)
  % Poll Augmented directions
  numberOfXtrials = size(AugmentedDirs,2);
  sites = struct('x',cell(numberOfXtrials,1,1),'f',cell(numberOfXtrials,1,1));
  for k = 1:numberOfXtrials
    direction = AugmentedDirs(:,k); % Always positively generated directions
    sites(k).x = Iterate.x + MeshSize*scale.*direction;
  end
  [msg,nextIterate,dirNum,optimState.FunEval] = ...
    psnextfeasible(ObjFunc,Xin,sites,Iterate,Aineq,bineq,Aeq,beq,lb,ub,tol,constr,objFcnArg,optimState.FunEval,options);
  Successdir = AugmentedDirs(:,dirNum); % dirNum is [] if msg = false
else
  Successdir = [];
end
% Update successful direction in optimState
optimState.Successdir = Successdir;

