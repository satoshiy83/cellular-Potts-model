
classdef CPMPolarizedSimulator < CPMSimulator
properties
    polaritySwitch = nan; % SYData.
end

methods
function obj = CPMPolarizedSimulator
% Controller class managing a simulation.
% obj = CPMPolarizedSimulator
    obj.mapClass = 'CPMMap';
    obj.cellClass = 'CPMPolarizedCell';
end

function dest = copy(obj,dest)
% Method to make a shallow copy.
% dest = copy(obj,dest)
% Return value is an instance of the class.
    if nargin < 2
        dest = CPMPolarizedSimulator;
    end
    copy@CPMSimulator(obj,dest);

    dest.polaritySwitch = obj.polaritySwitch;
end

function readHint(obj)
% Method to import parameters from the hint.
% readHint(obj)
% The hint includes:
%   CPMPolaritySwitchKey: (SYData) {n,2} matrix specifying cell types
%       either apico-basally polarized or not.
    readHint@CPMSimulator(obj);
    
    obj.polaritySwitch = obj.hint.objectForKey("CPMPolaritySwitchKey");
end

function markSubcellularLocations(obj)
% Method to make cells mark subcellular locations.
% markSubcellularLocations(obj)
    markSubcellularLocations@CPMSimulator(obj);
    
    obj.markApicalSurface;
end
function markApicalSurface(obj)
% Method to mark tissue apical surface.
% markApicalSurface(obj)
    bitmap = zeros(obj.map.frameSize);
    citmap = zeros(obj.map.frameSize);
    
    for i = 1:obj.cellArray.count
        pc = obj.cellArray.objectAtIndex(i);
        
        n = find(any(obj.polaritySwitch.var == pc.cellType,1),1);
        if n == 1
            continue;
        end
        
        pc.markLateral;
        pc.markApical(obj.map.labelDict.objectForKey("EC-apical"));
        
        bitmap(pc.apicoBasalRim) = 1;
        citmap(pc.apicalRim) = 1;
    end
    
    bitmap = IPConnectedComponents.connectedBinaryComponents(bitmap,8);
    for i = 1:max(bitmap(:))
        ditmap = (bitmap == i) & citmap;
        if ~any(ditmap(:))
            bitmap(bitmap == i) = 0;
        end
    end
    
    maskData = SYData(bitmap > 0);
    for i = 1:obj.cellArray.count
        pc = obj.cellArray.objectAtIndex(i);
        
        n = find(any(obj.polaritySwitch.var == pc.cellType,1),1);
        if n == 2 % polarized cell.
            pc.readApicalRim(maskData);
        end
    end
end

function updateCellsInList(obj)
% Method to make cells around the changed point receive a change.
% updateCellsInList(obj)
    updateCellsInList@CPMSimulator(obj);
    
    for i = obj.updateList
        cell = obj.cellWithLabel(i);
        cell.updateRim;
    end
end

end
end
