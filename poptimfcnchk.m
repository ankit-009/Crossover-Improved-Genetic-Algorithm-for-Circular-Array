function [Iterate,count] = poptimfcnchk(FUN,nonlcon,Xin,Iterate,Vectorized,objFcnArg,conFcnArg)
%POPTIMFCNCHK Calls objective function 'FUN' and nonlinear constraint.
%   function 'nonlcon' for the first time at the start point 'Iterate.x'. If
%   'Vectorized' option in 'on' the objective function is called with two
%   points.

%   Private to PATTERNSEARCH.


%   Copyright 2003-2012 The MathWorks, Inc.
%   $Revision: 1.1.6.6 $  $Date: 2012/08/21 00:24:44 $

y = NaN;
count = 0;
X = Iterate.x;
% Check the objective function. 
FUN = fcnchk(FUN);
if strcmpi(Vectorized,'off')  % The function is not vectorized. 
    try
        [y,count] = funevaluate(FUN,Xin,X,'init',[],[],objFcnArg{:});
    catch userFcn_ME
        gads_ME = MException('globaloptim:poptimfcnchk:objfunCheck', ...
            'Failure in initial user-supplied objective function evaluation. PATTERNSEARCH cannot continue.');
        userFcn_ME = addCause(userFcn_ME,gads_ME);
        rethrow(userFcn_ME)
   end
    if numel(y) ~=1
        error(message('globaloptim:poptimfcnchk:objfunCheck'));
    end
elseif strcmpi(Vectorized,'on') % If vectorized is 'on', Completepoll MUST be 'on' too
    X2 = [X, X];
    try
        [f,count] = funevaluate(FUN,Xin,X2,'init',[],[],objFcnArg{:});
    catch userFcn_ME
        gads_ME = MException('globaloptim:poptimfcnchk:objfunCheck', ...
            'Failure in initial user-supplied objective function evaluation (hint: ''vectorized'' option is ''on''). PATTERNSEARCH cannot continue.');
        userFcn_ME = addCause(userFcn_ME,gads_ME);
        rethrow(userFcn_ME)
    end
    if 2 ~=numel(f)
        error(message('globaloptim:poptimfcnchk:objfunCheckVectorized', ... 
            'Vectorized','on'));
    end
    y = f(1);
end
% We want the function value to be real
if isnan(y)
    error(message('globaloptim:poptimfcnchk:objfunNaN'));
end
Iterate.f = y;

% Evaluate nonlinear constraints for the first time
if  ~isempty(nonlcon)
    try
        if strcmpi(Vectorized,'on')
            % Evaluate the nonlinear constraints at two points; X2 = [X,X].
            [tmpCineq,tmpCeq] = feval(nonlcon,reshapeinput(Xin,X2),conFcnArg{:});            
        else
            [cineq,ceq] = feval(nonlcon,reshapeinput(Xin,X),conFcnArg{:});
        end
    catch userFcn_ME
        gads_ME = MException('globaloptim:poptimfcnchk:confunCheck', ...
            'Failure in initial user-supplied nonlinear constraint function evaluation. PATTERNSEARCH cannot continue.');
        userFcn_ME = addCause(userFcn_ME,gads_ME);
        rethrow(userFcn_ME)
    end
    if strcmpi(Vectorized,'on')
        % Is the initial point a column or row vector? A scalar is defined
        % as row vector. Note that we want Xin to be a vector.
        isXinColumnVector = size(Xin,1) > size(Xin,2);
        % Check whether row/column dimension of the constraint matrix is
        % two, corresponding to two points.         
        % We don't perform any additional checks for the two-constraint
        % case.
        % If Xin is a row vector, 2^isXinColumnVector) = 1, therefore the
        % number of rows in the constraint matrices should be 2.
        % If Xin is a column vector, 2^isXinColumnVector) = 2, therefore the
        % number of columns in the constraint matrices should be 2.                
        if ( ~isempty(tmpCineq) && ( size(tmpCineq,2^isXinColumnVector) ~= 2 ) ) || ...
                ( ~isempty(tmpCeq) && ( size(tmpCeq,2^isXinColumnVector) ~= 2 ) )            
            if isXinColumnVector
                error(message('globaloptim:poptimfcnchk:x0ColConfunCheckVectorized', ...
                    'Vectorized','on'));
            else
                error(message('globaloptim:poptimfcnchk:x0RowConfunCheckVectorized', ...
                    'Vectorized','on'));
            end
        end
        cineq = [];
        ceq = [];
        if isXinColumnVector
            if ~isempty(tmpCineq)
                cineq = tmpCineq(:,1);
            end
            if ~isempty(tmpCeq)
                ceq = tmpCeq(:,1);
            end                        
        else
            if ~isempty(tmpCineq)
                cineq = tmpCineq(1,:);
            end
            if ~isempty(tmpCeq)             
                ceq = tmpCeq(1,:);
            end                                    
        end
    end
    Iterate.cineq = zeros(numel(cineq),1);
    Iterate.ceq = zeros(numel(ceq),1);
    Iterate.cineq(:) = cineq;
    Iterate.ceq(:) = ceq;
    c = [Iterate.cineq;Iterate.ceq];
    if ~all(isreal(c) & isfinite(c))
        error(message('globaloptim:poptimfcnchk:confunNotReal'));
    end
end