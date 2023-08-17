% CPM Class CPMAAPOSimulator < CPMPCPSimulator.
% Written by Satoshi Yamashita.
% Controller class managing a simulation with adjusted Ao and Po values.

classdef CPMAAPOSimulator < CPMPCPSimulator
properties
    PoProjectorArray = nan; % SYArray.
    AoProjectorArray = nan; % SYArray.
end

methods
function obj = CPMAAPOSimulator
% Controller class managing a simulation.
% obj = CPMAAPOSimulator
    obj.mapClass = 'CPMPCPMap';
    obj.cellClass = 'CPMAAPOCell';
end

function dest = copy(obj,dest)
% Method to make a shallow copy.
% dest = copy(obj,dest)
% Return value is an instance of the class.
    if nargin < 2
        dest = CPMAAPOSimulator;
    end
    copy@CPMPCPSimulator(obj,dest);

    dest.PoProjectorArray = obj.PoProjectorArray;
    dest.AoProjectorArray = obj.AoProjectorArray;
end

function readHint(obj)
% Method to import parameters from the hint.
% readHint(obj)
% The hint includes:
%   CPMPoProjectorArrayKey: (SYArray) [SYData] matrix whose row corresponds
%       to a perimeter segment and i-th column represents a coefficient for
%       a term of i-1 degree.
%   CPMAoProjectorArrayKey: (SYArray) [SYData] matrix whose row corresponds
%       to a cell total area or subcellular component and i-th column
%       represents a coefficient for a term of i-1 degree.
    readHint@CPMPCPSimulator(obj);
    
    obj.PoProjectorArray = obj.hint.objectForKey("CPMPoProjectorArrayKey");
    obj.AoProjectorArray = obj.hint.objectForKey("CPMAoProjectorArrayKey");
end

end
end
