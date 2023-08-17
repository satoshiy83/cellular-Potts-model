% CPM Class CPMPCPSimulator < CPMEpithelialSimulator.
% Written by Satoshi Yamashita.
% Controller class managing a simulation of planarly polarized tissue.

classdef CPMPCPSimulator < CPMEpithelialSimulator
properties
    PCPSwitch = nan; % SYData.
end

methods
function obj = CPMPCPSimulator
% Controller class managing a simulation.
% obj = CPMPCPSimulator
    obj.mapClass = 'CPMPCPMap';
    obj.cellClass = 'CPMPCPCell';
end

function dest = copy(obj,dest)
% Method to make a shallow copy.
% dest = copy(obj,dest)
% Return value is an instance of the class.
    if nargin < 2
        dest = CPMPCPSimulator;
    end
    copy@CPMEpithelialSimulator(obj,dest);

    dest.PCPSwitch = obj.PCPSwitch;
end

function readHint(obj)
% Method to import parameters from the hint.
% readHint(obj)
% The hint includes:
%   CPMPCPSwitchKey: (SYData) {n,2} matrix specifying cell types either
%       planarly polarized or not.
    readHint@CPMEpithelialSimulator(obj);
    
    obj.PCPSwitch = obj.hint.objectForKey("CPMPCPSwitchKey");
end

function markSubcellularLocations(obj)
% Method to make cells mark subcellular locations.
% markSubcellularLocations(obj)
    markSubcellularLocations@CPMEpithelialSimulator(obj);
    
    obj.map.setUpstreamBoundary;
    obj.propagatePCP;
end
function propagatePCP(obj)
% Method to set PCP.
% propagatePCP(obj)
    flag = true;
    while flag
        flag = false;
        for i = 1:obj.cellArray.count
            cell = obj.cellArray.objectAtIndex(i);
            if cell.polarizeCell
                flag = true;
            end
        end
    end
end

function updateCellsInList(obj)
% Method to make cells around the changed point receive a change.
% updateCellsInList(obj)
    updateCellsInList@CPMEpithelialSimulator(obj);
    
    for i = obj.updateList
        cell = obj.cellWithLabel(i);
        cell.clearPCP;
        cell.clearPCPRim;
    end
    
    flag = true;
    while flag
        flag = false;
        for i = obj.updateList
            cell = obj.cellWithLabel(i);
            if cell.polarizeCell
                flag = true;
            end
        end
    end
end

end
end
