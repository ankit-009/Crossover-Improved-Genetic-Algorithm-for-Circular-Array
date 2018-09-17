function [x,fval,exitFlag,output,population,scores] = gaunc(FitnessFcn,GenomeLength,options,output,Iterate)
%GAUNC Genetic algorithm unconstrained solver.
%   X = GAUNC(FITNESSFCN,NVARS) finds the minimum of FITNESSFCN using
%   GAUNC. NVARS is the dimension (number of design variables) of the
%   FITNESSFCN. FITNESSFCN accepts a vector X of size 1-by-NVARS,
%   and returns a scalar evaluated at X. 
% 		
%   X = GAUNC(FITNESSFCN,NVARS,OPTIONS) finds the minimum for
%   FITNESSFCN with the default optimization parameters replaced by values
%   in the structure OPTIONS. OPTIONS can be created with the GAOPTIMSET
%   function.
% 		
%   [X, FVAL] = GAUNC(FITNESSFCN, ...) returns FVAL, the value of the fitness
%   function FITNESSFCN at the solution X.
% 		
%   [X,FVAL,REASON] = GAUNC(FITNESSFCN, ...) returns the REASON for stopping.
%
%   [X,FVAL,REASON,OUTPUT] = GAUNC(FITNESSFCN, ...) returns a
%   structure OUTPUT with the following information: 
%             rngstate: <State of the random number generator before GAUNC started>
%          generations: <Total generations, excluding HybridFcn iterations>
%            funccount: <Total function evaluations>
%              message: <GAUNC termination message>
%
%   [X,FVAL,REASON,OUTPUT,POPULATION] = GAUNC(FITNESSFCN, ...) returns the final
%   POPULATION at termination.
% 		
%   [X,FVAL,REASON,OUTPUT,POPULATION,SCORES] = GAUNC(FITNESSFCN, ...) returns the
%   SCORES of the final POPULATION.
% 		
%   See also GAOPTIMSET, FITNESSFUNCTION, PATTERNSEARCH, @.

%   Copyright 2005-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.7 $  $Date: 2012/08/21 00:24:11 $



exitFlag=[];

% Create initial state: population, scores, status data
state = makeState(GenomeLength,FitnessFcn,Iterate,output.problemtype,options);
% Determine who is the caller
callStack = dbstack;
caller = callStack(2).file(1:end-2);

% Set state for plot and output functions (only gacon will have
% 'interrupt' state)
if ~strcmp(caller,'gacon')
    currentState = 'init';
else
    currentState = 'interrupt';
end
% Give the plot/output Fcns a chance to do any initialization they need.
state = gadsplot(options,state,currentState,'Genetic Algorithm');
[state,options] = gaoutput(FitnessFcn,options,state,currentState);

% Setup display header 
if  options.Verbosity > 1
    fprintf('\n                               Best           Mean      Stall\n');
    fprintf('Generation      f-count        f(x)           f(x)    Generations\n');
end
% Set state for plot and output functions (only gacon will have
% 'interrupt' state)
if ~strcmp(caller,'gacon')
    currentState = 'iter';
else
    currentState = 'interrupt';
end
% Run the main loop until some termination condition becomes true
while isempty(exitFlag)
        state.Generation = state.Generation + 1;
        % Repeat for each subpopulation (element of the populationSize vector)
        offset = 0;
        totalPop = options.PopulationSize;
        % Each sub-population loop
        for pop = 1:length(totalPop)
            populationSize =  totalPop(pop);
            thisPopulation = 1 + (offset:(offset + populationSize - 1));
            population = state.Population(thisPopulation,:);
            score = state.Score( thisPopulation );
            % Empty population is also possible
            if isempty(thisPopulation)
                continue;
            end
            [score,population,state] = stepGA(score,population,options,state,GenomeLength,FitnessFcn);
            
            % Store the results for this sub-population
            state.Population(thisPopulation,:) = population;
            state.Score(thisPopulation) = score;
            offset = offset + populationSize;
        end 
        
        % Remember the best score
        best = min(state.Score);
        generation = state.Generation;
        state.Best(generation) = best;
        
        % Keep track of improvement in the best
        if (generation > 1) && isfinite(best)
            if state.Best(generation-1) > best
                state.LastImprovement = generation;
                state.LastImprovementTime = cputime;
            end
        end
        
        % Do any migration
        state = migrate(FitnessFcn,GenomeLength,options,state);
        % Update the Output
        state = gadsplot(options,state,currentState,'Genetic Algorithm');
        [state,options] = gaoutput(FitnessFcn,options,state,currentState);

        % check to see if any stopping criteria have been met
        [exitFlag, output.message] = isItTimeToStop(options,state);
end % End while loop

% Find and return the best solution
[fval,best] = min(state.Score);
x = state.Population(best,:);

% Update output structure
output.generations = state.Generation;
output.funccount   = state.FunEval;

population = state.Population;
scores = state.Score;

% A hybrid scheme. Try another minimization method if there is one.
if ~isempty(options.HybridFcn)
    if strcmpi(options.PopulationType,'doubleVector')
        [x,fval] = callHybridFunction;
    else
        warning(message('globaloptim:gaunc:notValidHybrid'));
    end
end
% Set state for plot and output functions (only gacon will have
% 'interrupt' state)
if ~strcmp(caller,'gacon')
    currentState = 'done';
else
    currentState = 'interrupt';
end
% Give the Output functions a chance to finish up
gadsplot(options,state,currentState,'Genetic Algorithm');
gaoutput(FitnessFcn,options,state,currentState);

    function [xhybrid,fhybrid] = callHybridFunction
        xhybrid = x;
        fhybrid = fval;
        % Who is the hybrid function
        if isa(options.HybridFcn,'function_handle')
            hfunc = func2str(options.HybridFcn);
        else
            hfunc = options.HybridFcn;
        end
        % Create functions handle to be passed to hybrid function
        FitnessHybridFcn = @(x) FitnessFcn(x,options.FitnessFcnArgs{:});
        % Inform about hybrid scheme
        if  options.Verbosity > 1
            fprintf('%s%s%s\n','Switching to the hybrid optimization algorithm (',upper(hfunc),').');
        end
        if ~any(strcmpi(hfunc,{'fminsearch','fminunc','patternsearch'}))
             warning(message('globaloptim:gaunc:constrainedHybridFcn',upper(hfunc),'FMINUNC'));
            hfunc = 'fminunc';
        end
        [x_temp,f_temp,funccount,theMessage] = callHybrid(hfunc,FitnessHybridFcn,x,options.HybridFcnArgs);
        output.funccount = output.funccount + funccount;
        output.message   = sprintf([output.message '\n', theMessage '\n']);

        % Inform about hybrid scheme termination
        if  options.Verbosity > 1
            fprintf('%s%s\n',upper(hfunc), ' terminated.');
        end
        % The solution returned from the hybrid function will always be
        % feasible, so we can accept it if it has a lower function value.
        if f_temp < fval
            fhybrid = f_temp;
            xhybrid = x_temp;
        end
    end
end
