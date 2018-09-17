function [X,Aineq,Bineq,Aeq,Beq,msg,how,exitflag]= eqnsolv(X,Aineq,Bineq,Aeq,Beq,LB,UB,verbosity)
%EQNSOLV Helper function for PFMINLCON. 
%   Removes dependent equality constraints.

%   Copyright 2003-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2012/08/21 00:23:34 $

how = '';
msg = '';
exitflag = 1;
numberOfVariables = length(X);
tolDep = 100*numberOfVariables*eps;
% Make sure that the constraints are not sparse; 
if issparse(Aeq)
    warning(message('globaloptim:eqnsolv:eqSparseToFull'));
    Aeq = full(Aeq);
    Beq = full(Beq);
end
if issparse(Aineq)
    warning(message('globaloptim:eqnsolv:ineqSparseToFull'));
    Aineq = full(Aineq);
    Bineq  = full(Bineq);
end
A = [Aeq;Aineq];
B = [Beq;Bineq];
% Total number of constraints
ncstr   = size(A,1);
% Normalize constraints
for i=1:ncstr
    normRow = norm(A(i,:));
    if (normRow >= tolDep)
        A(i,:) = A(i,:)/normRow;
        B(i) = B(i)/normRow;
    end
end

% Handle bounds as linear constraints
arglb = ~eq(LB,-inf);
lenlb=length(LB); 
if nnz(arglb) > 0     
    lbmatrix = -eye(lenlb,numberOfVariables);
    A=[A; lbmatrix(arglb,1:numberOfVariables)]; % select finite bounds
    B=[B;-LB(arglb)];
end
argub = ~eq(UB,inf);
lenub=length(UB);
if nnz(argub) > 0
    ubmatrix = eye(lenub,numberOfVariables);
    A=[A; ubmatrix(argub,1:numberOfVariables)];
    B=[B; UB(argub)];
end

eqix = 1:size(Aeq,1);
% Form the active set index (active region is equality constraints
ACTIND = eqix;
% Number of equality constraints
neqcstr = length(eqix);
% Total number of constraints will also include the bounds
ncstr   = size(A,1);

% Index of equality constraints is 'eqix'
% A logical variable 'aix' to denote equality constraints
aix = zeros(ncstr,1);
% We set indices of equality constraints to true (ONE).
aix(ACTIND) = 1;
ACTCNT = length(ACTIND);
ACTSET = A(ACTIND,:);

remove =[];

% First see if the equality constraints form a consistent system.
[Qa,Ra,Ea]=qr(A(eqix,:)); % Using QR for full matrix A so we need the 3rd output argument

% Form vector of dependent indices.
if min(size(Ra))==1 % Make sure Ra isn't a vector
    depInd = find( abs(Ra(1,1)) < tolDep);
else
    depInd = find( abs(diag(Ra)) < tolDep );
end
if neqcstr > numberOfVariables
    depInd = [depInd; ((numberOfVariables+1):neqcstr)'];
end      

if ~isempty(depInd)    % equality constraints are dependent
    msg = sprintf('The equality constraints are dependent.\n');
    if verbosity > 2
        fprintf(msg);
    end
    how='dependent';
    exitflag = 1;
    bdepInd =  abs(Qa(:,depInd)'*B(eqix)) >= tolDep ;
    
    if any( bdepInd ) % Not consistent
        how='infeasible';   
        exitflag = -2;
        msg = [msg,sprintf('The system of equality constraints is not consistent.\n')];
        if ncstr > neqcstr
            msg = [msg sprintf('The inequality constraints may or may not be satisfied.\n')];
        end
        msg = [msg, sprintf('There is no feasible solution.\n')];
   else % the equality constraints are consistent
        % Delete the redundant constraints
        % By QR factoring the transpose, we see which columns of A'
        %   (rows of A) move to the end
        [Qat,Rat,Eat]=qr(A(eqix,:)');        
        [i,j] = find(Eat); % Eat permutes the columns of A' (rows of A)
        remove = i(depInd);
        numDepend = nnz(remove);
        if verbosity > 2
            fprintf('The system of equality constraints is consistent.\n');
            fprintf('Removing the following dependent constraints before continuing:\n');
            disp(remove)
        end
        A(eqix(remove),:)=[];
        B(eqix(remove))=[];
        %Taking care of equality
        Aeq(eqix(remove),:)=[];
        Beq(eqix(remove))=[];
        neqcstr = neqcstr - numDepend;
        ncstr = ncstr - numDepend;
        eqix = 1:neqcstr;
        aix(remove) = [];
        ACTIND(1:numDepend) = [];
        ACTIND = ACTIND - numDepend;      
        ACTSET = A(ACTIND,:);
        ACTCNT = ACTCNT - numDepend;
    end % consistency check
end % dependency check

% Now that we have done all we can to make the equality constraints
% consistent and independent we will check the inequality constraints
% in the working set.  First we want to make sure that the number of 
% constraints in the working set is only greater than or equal to the
% number of variables if the number of (non-redundant) equality 
% constraints is greater than or equal to the number of variables.
if ACTCNT >= numberOfVariables
    ACTCNT = max(neqcstr, numberOfVariables-1);
    ACTIND = ACTIND(1:ACTCNT);
    ACTSET = A(ACTIND,:);
    aix = zeros(ncstr,1);
    aix(ACTIND) = 1;
end

% Now check to see that all the constraints are linearly independent.
if ACTCNT > neqcstr
    [Qat,Rat,Eat]=qr(ACTSET');
    
    % Form vector of dependent indices.
    if min(size(Rat))==1 % Make sure Rat isn't a vector
        depInd = find( abs(Rat(1,1)) < tolDep);
    else
        depInd = find( abs(diag(Rat)) < tolDep );
    end
    
    if ~isempty(depInd)
        [i,j] = find(Eat); % Eat permutes the columns of A' (rows of A)
        remove2 = i(depInd);
        removeEq   = remove2(find(remove2 <= neqcstr));
        removeIneq = remove2(find(remove2 > neqcstr));
        
        if ~isempty(removeEq)
            % Just take equalities as initial working set.
            ACTIND = 1:neqcstr; 
        else
            % Remove dependent inequality constraints.
            ACTIND(removeIneq) = [];
        end
        aix = zeros(ncstr,1);
        aix(ACTIND) = 1;
        ACTSET = A(ACTIND,:);
        ACTCNT = length(ACTIND);
    end  
end

[Q,R,E]=qr(ACTSET');

if ~strcmp(how,'infeasible') && ACTCNT > 0
    % Find point closest to the given initial X which satisfies
    % working set constraints.
    minnormstep = Q(:,1:ACTCNT) * ...
        ((R(1:ACTCNT,1:ACTCNT)') \ (B(ACTIND) - ACTSET*X));
    X = X + minnormstep; 
    % Sometimes the "basic" solution satisfies Aeq*x= Beq 
    % and A*X < B better than the minnorm solution. Choose the one
    % that the minimizes the max constraint violation.
    err = A*X - B;
    err(eqix) = abs(err(eqix));
    if any(err > eps)
        Xbasic = ACTSET\B(ACTIND);
        errbasic = A*Xbasic - B;
        errbasic(eqix) = abs(errbasic(eqix));
        if max(errbasic) < max(err) 
            X = Xbasic;
        end
    end
end

% End of eqnsolv.m


