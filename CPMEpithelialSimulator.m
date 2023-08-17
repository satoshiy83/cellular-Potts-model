% CPM Class CPMEpithelialSimulator < CPMPolarizedSimulator.
% Written by Satoshi Yamashita.
% Controller class managing a simulation of epithelial tissue.

classdef CPMEpithelialSimulator < CPMPolarizedSimulator
properties
    epitheliaSwitch = nan; % SYData.
end

methods
function obj = CPMEpithelialSimulator
% Controller class managing a simulation.
% obj = CPMEpithelialSimulator
    obj.mapClass = 'CPMMap';
    obj.cellClass = 'CPMEpithelialCell';
end

function dest = copy(obj,dest)
% Method to make a shallow copy.
% dest = copy(obj,dest)
% Return value is an instance of the class.
    if nargin < 2
        dest = CPMEpithelialSimulator;
    end
    copy@CPMPolarizedSimulator(obj,dest);
end

function readHint(obj)
% Method to import parameters from the hint.
% readHint(obj)
% The hint includes:
%   CPMEpitheliaSwitchKey: (SYData) {n,2} matrix specifying cell types
%       either epithelial or not.
    readHint@CPMPolarizedSimulator(obj);
    
    obj.epitheliaSwitch = obj.hint.objectForKey("CPMEpitheliaSwitchKey");
end

function markSubcellularLocations(obj)
% Method to make cells mark subcellular locations.
% markSubcellularLocations(obj)
    markSubcellularLocations@CPMPolarizedSimulator(obj);
    
    obj.cellArray.makeObjectsPerformSelector(@markAdherensJunction);
end

function updateCellsInList(obj)
% Method to make cells around the changed point receive a change.
% updateCellsInList(obj)
    updateCellsInList@CPMPolarizedSimulator(obj);
    
    for i = obj.updateList
        cell = obj.cellWithLabel(i);
        cell.markAdherensJunction;
    end
end

end
end
