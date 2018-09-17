function [X, A, B, Aeq, Beq, LB, UB] = checkglobalsearchnlpinputs(X, A, B, Aeq, Beq, LB, UB)
%CHECKGLOBALSEARCHNLPINPUTS Check inputs to globalsearchnlp.
%
%   [X, A, B, AEQ, BEQ, LB, UB] = CHECKGLOBALSEARCHNLPINPUTS(X, A, B, AEQ,
%   BEQ, LB, UB) performs checks on the inputs to GLOBALSEARCHNLP. This
%   function provides the following input checks:
%
%   1. All variables are real.
%   2. The linear constraints are specified as 2-d matrices.
%   3. Bounds are sensible.
% 
%   This function is private to the Global Optimization Toolbox.

%   Copyright 2009-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2012/08/21 00:23:26 $

try

    % Ensure variables are real
    i_checkReal(X, 'X');
    i_checkReal(A, 'A');
    i_checkReal(B, 'B');
    i_checkReal(Aeq, 'Aeq');
    i_checkReal(Beq, 'Beq');
    i_checkReal(LB, 'LB');
    i_checkReal(UB, 'UB');
    
    % Ensure linear constraint matrices are 2-d only. This check should be in
    % fmincon.
    i_check2dMatrix(A, 'A');
    i_check2dMatrix(B, 'B');
    i_check2dMatrix(Aeq, 'Aeq');
    i_check2dMatrix(Beq, 'Beq');
       
    % Ensure bounds are sensible.
    try
        [X,LB,UB,msg] = checkbounds(X,LB,UB,numel(X));
    catch ME
       msg = ME.message;
    end
    if ~isempty(msg)
        error(message('globaloptim:globalsearchnlp:IncorrectBounds'));
    end
    
catch ME
    throwAsCaller(ME);
end

function i_check2dMatrix(X, name)

if ndims(X) ~= 2
    error(message('globaloptim:globalsearchnlp:InvalidArgumentNot2D', name));
end

function i_checkReal(X, name)

if ~isreal(X)
   error(message('globaloptim:globalsearchnlp:InvalidArgumentNotReal', name));
end 

