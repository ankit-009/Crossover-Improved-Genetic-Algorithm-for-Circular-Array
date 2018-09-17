function [state,options,optchanged] = gaoutput(~,options,state,flag)
%GAOUTPUT Helper function that manages the output functions for GA.
%
%   [STATE, OPTIONS, OPTCHANGED] = GAOUTPUT(~, OPTIONS, STATE, FLAG)
%   runs each of the display functions in the options.OutputFcns cell
%   array.
%
%   This is a helper function called by ga between each generation, and is
%   not typically called directly.

%   Copyright 2003-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2012/10/24 04:18:01 $


% get the functions and return if there are none
optchanged = false;
functions = options.OutputFcns;
if(isempty(functions))
    return
end

% call each output function
stopFlag = [];
args = options.OutputFcnsArgs;
for i = 1:length(functions)
    % Always clear state.StopFlag before calling OutputFcns so that each
    % function can independently set a termination message.
    state.StopFlag = '';
    [state,optnew,changed] = feval(functions{i},options,state,flag,args{i}{:});
    if ~isempty(state.StopFlag)
        stopFlag = [stopFlag state.StopFlag ';'];
    end
    if changed %If changes are not duplicates, we will get all the changes
       % Keep LinearConstr out of options and accept new options
       LinearConstr = options.LinearConstr;
       type = LinearConstr.type;
       gLength = size(state.Population,2);
       % Store integer constraints information
       [intcon, UserSpecPopInitRange, UserVectorized] = ...
           i_storeIntegerInfo(options);
       options = gaoptimset(options,optnew);
       options = validate(options,type,gLength,[],[]);
       options.LinearConstr = LinearConstr;
       % Perform integer constraint specific tasks
       options = i_resetIntegerInfo(options, intcon, ...
           UserSpecPopInitRange, UserVectorized, gLength, LinearConstr);
       optchanged = true;
    end
end

state.StopFlag = stopFlag;

function [intcon, UserSpecPopInitRange, UserVectorized] = ...
    i_storeIntegerInfo(options)

intcon = options.IntegerVars;
if isempty(intcon)
    UserSpecPopInitRange = [];
    UserVectorized = [];
else
    % Keep UserSpecPopInitRange and UserVectorized out of options
    UserSpecPopInitRange = options.UserSpecPopInitRange;
    UserVectorized = options.UserVectorized;
end

function options = i_resetIntegerInfo(options, intcon, ...
    UserSpecPopInitRange, UserVectorized, gLength, LinearConstr)

% Add IntegerVars field back to options structure
options.IntegerVars = intcon;

% If there are any integer constraints, need to perform extra
% validation checks
if ~isempty(intcon)
    % Add UserSpecPopInitRange back to options
    options.UserSpecPopInitRange = UserSpecPopInitRange;
    % Validate mixed integer GA options
    gaminlpvalidateoptions(options);
    % It is possible for users to alter the Vectorized option in an output
    % function. Say the user originally set the Vectorized option to 'off'
    % and they set it to 'on' in their output function. Because we set
    % options.Vectorized = 'on' for mixed integer GA, it will look like
    % this option has not changed. As such, we cannot tell whether or not
    % the user changed options.Vectorized to 'on' in an output function.
    % Note that we can determine if a user sets options.Vectorized to
    % 'off'.
    %
    % So, for mixed integer GA, we ignore any change that a user makes to
    % the Vectorized option in an output function. We will use
    % gaminlpoverwriteoptions to overwrite the Vectorized option (to 'on')
    % and set the internal UserVectorized option to reflect the Vectorized
    % status when the user called GA.
    % 
    % gaminlpoverwriteoptions sets up UserVectorized from the value in the
    % Vectorized option, so we set options.Vectorized to reflect the value
    % in UserVectorized here.
    if UserVectorized
        options.Vectorized = 'on';
    else
        options.Vectorized = 'off';
    end
    % The call to gaoptimset removes the mixed integer GA specific
    % options. We need to recreate them here.
    options = gaminlpoverwriteoptions(options, gLength, ...
        LinearConstr.lb, LinearConstr.ub, intcon);
end
