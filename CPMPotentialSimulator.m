% CPM Class CPMPotentialSimulator < CPMAAPOSimulator.
% Written by Satoshi Yamashita.
% Controller class managing a simulation with potential energy.

classdef CPMPotentialSimulator < CPMAAPOSimulator
% properties (Access = private)
%     cellClass = 'CPMPotentialCell';
% end
properties
    potentialComponentArray = nan; % SYArray.
end

methods
function obj = CPMPotentialSimulator
% Controller class managing a simulation.
% obj = CPMPotentialSimulator
    obj.mapClass = 'CPMPotentialMap';
    obj.cellClass = 'CPMPotentialCell';
end
function obj = initWithMap(obj,map,hint)
% Initialization method with CPMPotentialMap instance and hint.
% obj = initWithMap(obj,map,hint)
% Argument hint includes:
%   CPMPotentialComponentArrayKey: (SYArray) [SYArray [SYData]] matrices
%       specifying subcellular locations and components and slice of the
%       potential energy field stack.
    if isnumeric(map)
        map = CPMPotentialMap(map(:,:,1),map(:,:,2), ...
            map(:,:,3),map(:,:,4:end),hint);
    end
    initWithMap@CPMAAPOSimulator(obj,map,hint);
end

function dest = copy(obj,dest)
% Method to make a shallow copy.
% dest = copy(obj,dest)
% Return value is an instance of the class.
    if nargin < 2
        dest = CPMPotentialSimulator;
    end
    copy@CPMAAPOSimulator(obj,dest);
end

function readHint(obj)
% Method to import parameters from the hint.
% readHint(obj)
% Argument hint includes:
%   CPMPotentialComponentArrayKey: (SYArray) [SYArray [SYData]] matrices
%       specifying subcellular locations and components and slice of the
%       potential energy field stack.
    readHint@CPMAAPOSimulator(obj);
    
    obj.potentialComponentArray = ...
        obj.hint.objectForKey("CPMPotentialComponentArrayKey");
end
% function convertNameToIndexInPotentialComponentsArray(obj)
%     for i = 1:obj.potentialComponentArray.count
%         pMatrixArray = obj.potentialComponentArray.objectAtIndex(i);
%         pMatrix = pMatrixArray.objectAtIndex(1);
%         if isa(pMatrix,'SYDictionary')
%             data = obj.convertNameToPMatrixSubcellularComponents(pMatrix);
%             pMatrixArray.replaceObjectAtIndex(data,1);
%         end
%         pMatrix = pMatrixArray.objectAtIndex(2);
%         if isa(pMatrix,'SYDictionary')
%             data = obj.convertNameToPMatrixSubcellularLocations(pMatrix);
%             pMatrixArray.replaceObjectAtIndex(data,2);
%         end
%     end
% end
% function result = convertNameToPMatrixSubcellularComponents(obj,dict)
%     if dict.count < 1
%         result = SYData;
%         return;
%     end
%     keys = dict.allKeys;
%     for i = 1:keys.count
%         key = keys.objectAtIndex(i);
%         
%     end
% end
% function result = convertNameToPMatrixSubcellularLocations(obj,dict)
%     
% end

end
end
