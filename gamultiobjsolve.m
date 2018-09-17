function [x,fval,exitFlag,output,population,scores] = gamultiobjsolve(FitnessFcn,GenomeLength, ...
     Aineq,bineq,Aeq,beq,lb,ub,options,output)
%GAMULTIOBJSOLVE Genetic algorithm multi-objective solver.

%   Copyright 2007-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2012/08/21 00:24:09 $



% Create initial state: population, scores, status data
state = gamultiobjMakeState(GenomeLength,FitnessFcn,output.problemtype,options);

currentState = 'init';
% Give the plot/output Fcns a chance to do any initialization they need.
state = gadsplot(options,state,currentState,'Genetic Algorithm');
[state,options] = gaoutput(FitnessFcn,options,state,currentState);

% Setup display header 
if  options.Verbosity > 1
    fprintf('\n                           Average            Average\n');
    fprintf('Generation   f-count    Pareto distance    Pareto spread\n');
end

currentState = 'iter';
% Run the main loop until some termination condition becomes true
exitFlag = [];
while true
       state.Generation = state.Generation + 1;
        % check to see if any stopping criteria have been met
       [state,exitFlag,reasonToStop] = gamultiobjConverged(options,state);
       if ~isempty(exitFlag)
           break;
       end
       
        % Repeat for each sub-population (element of the PopulationSize vector)
        offset = 0;
        totalPop = options.PopulationSize;
        % Each sub-population loop
        for pop = 1:length(totalPop)
            populationSize =  totalPop(pop);
            thisPopulation = 1 + (offset:(offset + populationSize - 1));
            % Empty population is also possible
            if isempty(thisPopulation)
                continue; 
            end
            state = stepgamultiobj(pop,thisPopulation,options,state,GenomeLength,FitnessFcn);
            offset = offset + populationSize;
        end 
        
        % Migration
        state = migrate(FitnessFcn,GenomeLength,options,state);

        % Output and plot functions
        state = gadsplot(options,state,currentState,'Genetic Algorithm');
        [state,options] = gaoutput(FitnessFcn,options,state,currentState);
end % End while loop
% Update output structure
output.generations = state.Generation;
output.message = reasonToStop;

% If sub-population model is used, merge all sub-population and perform
% another non-dominated sorting
if length(options.PopulationSize) > 1
    [state.Population,state.Score,state.Rank,state.Distance]  = ...
        rankAndDistance(state.Population,state.Score,options);
    % Calculate average distance and spread
    [output.averagedistance,output.spread] = distanceAndSpread(state.Distance, ...
        state.Rank,state.Score,state.Score);
else
    % Calculate front statistics for output structure
    output.averagedistance = state.AverageDistance;
    output.spread = state.Spread(end);
end

% Find and return the solutions on Pareto front
fval = state.Score(state.Rank == 1,:);
x = state.Population((state.Rank == 1),:);


% A hybrid scheme; try another minimization method
if ~isempty(options.HybridFcn)
    if strcmpi(options.PopulationType,'doubleVector')
        state  = gamultiobjHybrid(FitnessFcn,x,fval,state,Aineq,bineq,Aeq,beq,lb,ub,options);
        % Calculate front statistics for output structure
        [output.averagedistance,output.spread] = distanceAndSpread(state.Distance,state.Rank,state.Score,state.Score);
        % Find and return the solutions on Pareto front
        fval = state.Score(state.Rank == 1,:);
        x = state.Population((state.Rank == 1),:);
    else
        warning(message('globaloptim:gamultiobjsolve:notValidHybrid'));
    end
end

output.funccount   = state.FunEval;
population = state.Population;
scores = state.Score;

currentState = 'done';
% Give the Output functions a chance to finish up
gadsplot(options,state,currentState,'Genetic Algorithm');
gaoutput(FitnessFcn,options,state,currentState);

