% CPM Class CPMPCPMap < CPMMap.
% Written by Satoshi Yamashita.
% Data class representing a lattice with PCP map.

classdef CPMPCPMap < CPMMap
properties
    PCPMap = nan; % double[h,w];
end

methods
function obj = initWithMaps(obj,cellMap,subcellMap,cellTypeMap,hint)
% Data class representing a lattice.
% obj = initWithMaps(obj,cellMap,subcellMap,cellTypeMap,hint)
%     initWithMaps@CPMMap(obj,cellMap,subcellMap,cellTypeMap,hint);
% Argument cellMap is a bitmap of cell labels.
% Argument subcellMap is a bitmap of subcellular components labels.
% Argument cellTypeMap is a bitmap of cell type labels.
% Argument hint is an SYDictionary instance holding parameters.
% Parameters:
%   CPMLabelDictKey: (SYDictionary) dictinary of labels refered in
%   subclass.
    initWithMaps@CPMMap(obj,cellMap,subcellMap,cellTypeMap,hint);
    
    obj.PCPMap = zeros(obj.frameSize);
    obj.setUpstreamBoundary;
end
function obj = initWithData(obj,data)
% Initialization method with SYData instance.
% obj = initWithData(obj,data)
    initWithData@CPMMap(obj,data);
    
    obj.PCPMap = zeros(obj.frameSize);
    obj.setUpstreamBoundary;
    s = data.var;
    obj.PCPMap = s.PCPMap;
end

function dest = copy(obj,dest)
% Method to make a shallow copy.
% dest = copy(obj,dest)
% Return value is an instance of the class.
    if nargin < 2
        dest = CPMPCPMap;
    end
    copy@CPMMap(obj,dest);
    
    dest.PCPMap = obj.PCPMap;
end

function result = data(obj)
% Method to convert the instance to an SYData instance.
% result = data(obj)
% Return value is an SYData instance.
    data_ = data@CPMMap(obj);
    s = data_.var;
    
    s.PCPMap = obj.PCPMap;
    
    data_ = SYData(s);
    result = data_;
end

function setUpstreamBoundary(obj)
% Method to mark PCP upstream boundary.
% setUpstreamBoundary(obj)
    label = obj.labelDict.objectForKey("PCP-upstream_boundary");
    mask = any(obj.subcellMap == label,3);
    obj.PCPMap(mask) = 1;
end

function obj = revert(obj)
% Method to reject an update.
% obj = revert(obj)
    revert@CPMMap(obj);
    
    pixel = obj.pointUpdated;
    obj.PCPMap(pixel) = 0;
end

end
end
