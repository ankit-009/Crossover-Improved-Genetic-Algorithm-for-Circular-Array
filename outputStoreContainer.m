function fHdls = outputStoreContainer
%OUTPUTSTORECONTAINER Container to store output information from a local
%solver.
%
%   FHDLS = OUTPUTSTORECONTAINER creates a structure of nested function
%   handles. The nested function workspace of this function holds a partial
%   optimization toolbox output structure. Specifically, this workspace
%   holds a structure, S, containing two fields, 'iterations' and
%   'funcCount'.
%
%   Once the container is created, the nested function handle structure can
%   be used to act on S. In particular, S can be updated during the
%   progress of a local solver call via an output function. The function
%   handle structure, FHDLS, contains three functions that can be used to
%   update S:
%
%   CURR_S = FHDLS.GETOUTPUTSTORE returns the current S in the nested
%   function workspace.
%
%   FHDLS.SETOUTPUTSTORE(NEW_S) replaces the current S in the nested
%   function workspace with NEW_S.
%
%   FHDLS.RESETOUTPUTSTORE resets S to be an empty structure with fields
%   'iterations' and 'funcCount'

%   Copyright 2009 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2012/08/21 00:24:36 $

OutputStore = i_createOutputStore;

fHdls.getOutputStore = @getOutputStore;
fHdls.setOutputStore = @setOutputStore;
fHdls.resetOutputStore = @resetOutputStore;

    function setOutputStore(NewOutputStore)
        OutputStore = NewOutputStore;
    end     

    function CurrentOutputStore = getOutputStore
        CurrentOutputStore = OutputStore;
    end

    function resetOutputStore
        OutputStore = i_createOutputStore;
    end
end


function os = i_createOutputStore

os = struct('iterations', 0, ...
    'funcCount', 0);

end