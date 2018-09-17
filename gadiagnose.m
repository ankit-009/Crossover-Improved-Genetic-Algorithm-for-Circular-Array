function gadiagnose(FUN,nonlcon,GenomeLength,nineqcstr,neqcstr,ncstr,nObj,options,caller)
%GADIAGNOSE prints some diagnostic information about the problem
%
%   This function is private to GA

%   Copyright 2004-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2012/08/21 00:23:47 $

properties =  fieldnames(gaoptimset);
defaultOpt = gaoptimset(caller);
Output_String = sprintf('\nDiagnostic information.');

Output_String = [Output_String sprintf('\n\tFitness function = %s',value2RHS(FUN))];
if ~isempty(GenomeLength)
    Output_String = [Output_String sprintf('\n\tNumber of variables = %d',GenomeLength)];
end
if nObj > 1
    Output_String = [Output_String sprintf('\n\tNumber of objectives = %d',nObj)];
end

%print some information about constraints
if ~isempty(nonlcon)
    Output_String = [Output_String sprintf('\n\tnonlinear constraint function = %s',value2RHS(nonlcon))];
end
if ~isempty(nineqcstr)
    Output_String = [Output_String sprintf('\n\t%d Inequality constraints',nineqcstr)];
end
if ~isempty(neqcstr)
    Output_String = [Output_String sprintf('\n\t%d Equality constraints',neqcstr)];
end
if ~isempty(ncstr)
    Output_String = [Output_String sprintf('\n\t%d Total number of linear constraints\n',ncstr)];
end

Output_String = [Output_String sprintf('\n%s','Modified options:')];
for i = 1:length(properties)
    prop = properties{i};
    if (~isempty(prop)) && isfield(options,prop) % the property list has blank lines, ignore them
        value = options.(prop);
        if ~(isequal(value,defaultOpt.(prop)) || isempty(value)) 
            Output_String = [Output_String sprintf('\n\toptions.%s = %s',prop,value2RHS(value))];
        end
    end
end
Output_String = [Output_String sprintf('\nEnd of diagnostic information.')];
fprintf('%s\n',Output_String)
