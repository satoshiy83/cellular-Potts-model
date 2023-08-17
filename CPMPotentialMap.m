% CPM Class CPMPotentialMap < CPMPCPMap.
% Written by Satoshi Yamashita.
% Data class representing a lattice with potential energy.

classdef CPMPotentialMap < CPMPCPMap
properties
    potentialMap = nan; % double[h,w,n].
    potentialMapsCount = nan; % double.
end

methods
function obj = ...
        CPMPotentialMap(cellMap,subcellMap,cellTypeMap,potentialMap,hint)
% Data class representing a lattice with potential energy.
% obj = CPMPotentialMap(cellMap,subcellMap,cellTypeMap,potentialMap,hint)
% Argument cellMap is a bitmap of cell labels.
% Argument subcellMap is a bitmap of subcellular components labels.
% Argument cellTypeMap is a bitmap of cell type labels.
% Argument potentialMap is a stack of potential energy fields.
% Argument hint is an SYDictionary instance holding parameters.
% Parameters:
%   CPMLabelDictKey: (SYDictionary) dictinary of labels refered in
%   subclass.
    if nargin < 6
        return;
    end
    obj.initWithMaps(cellMap,subcellMap,cellTypeMap,potentialMap,hint);
end
function obj = ...
        initWithMaps(obj,cellMap,subcellMap,cellTypeMap,potentialMap,hint)
% Initialization method with maps.
% obj = initWithMaps(obj,cellMap,subcellMap,cellTypeMap,potentialMap,hint)
% Argument cellMap is a bitmap of cell labels.
% Argument subcellMap is a bitmap of subcellular components labels.
% Argument cellTypeMap is a bitmap of cell type labels.
% Argument potentialMap is a stack of potential energy fields.
% Argument hint is an SYDictionary instance holding parameters.
% Parameters:
%   CPMLabelDictKey: (SYDictionary) dictinary of labels refered in
%   subclass.
    if nargin == 5
        hint = potentialMap;
    end
    initWithMaps@CPMPCPMap(obj,cellMap,subcellMap,cellTypeMap,hint);

    if nargin < 6
        return;
    end
    
    obj.potentialMap = potentialMap;
    obj.potentialMapsCount = size(potentialMap,3);
end
function obj = initWithData(obj,data)
% Initialization method with SYData instance.
% obj = initWithData(obj,data)
    obj = initWithData@CPMPCPMap(obj,data);
    
    s = data.var;
    obj.potentialMap = s.potentialMap_;
    obj.potentialMapsCount = size(obj.potentialMap,3);
end

function dest = copy(obj,dest)
% Method to make a shallow copy.
% dest = copy(obj,dest)
% Return value is an instance of the class.
    if nargin < 2
        dest = CPMPotentialMap;
    end
    copy@CPMPCPMap(obj,dest);
    
    dest.potentialMap = obj.potentialMap;
    dest.potentialMapsCount = obj.potentialMapsCount;
end

function result = data(obj)
% Method to convert the instance to an SYData instance.
% result = data(obj)
% Return value is an SYData instance.
    data_ = data@CPMPCPMap(obj);
    s = data_.var;
    s.potentialMap_ = obj.potentialMap;
    
    data_ = SYData(s);
    result = data_;
end

end
end
