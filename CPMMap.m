% CPM Class CPMMap < SYObject.
% Written by Satoshi Yamashita.
% Data class representing a lattice.

classdef CPMMap < SYObject
properties
    cellMap = nan; % uint16[h,w].
    subcellMap = nan; % uint16[h,w].
    cellTypeMap = nan; % uint16[h,w].
    
    hint = nan; % SYDictionary.
    labelDict = nan; % SYDictionary.
    
    frameSize = nan; % double[2].
    distanceMap = nan; % double[h,w].
    dmCenter = nan; % double[2].
    
    pointUpdated = nan; % double[2].
    labelUpdated = nan; % double[3].
end

methods
function obj = CPMMap(cellMap,subcellMap,cellTypeMap,hint)
% Data class representing a lattice.
% obj = CPMMap(cellMap,subcellMap,cellTypeMap,hint)
% Argument cellMap is a bitmap of cell labels.
% Argument subcellMap is a bitmap of subcellular components labels.
% Argument cellTypeMap is a bitmap of cell type labels.
% Argument hint is an SYDictionary instance holding parameters.
% Parameters:
%   CPMLabelDictKey: (SYDictionary) dictinary of labels refered in
%   subclass.
    if nargin < 4
        return;
    end
    obj.initWithMaps(cellMap,subcellMap,cellTypeMap,hint);
end
function obj = initWithMaps(obj,cellMap,subcellMap,cellTypeMap,hint)
% Initialization method with maps.
% obj = initWithMaps(obj,cellMap,subcellMap,cellTypeMap,hint)
% Argument cellMap is a bitmap of cell labels.
% Argument subcellMap is a bitmap of subcellular components labels.
% Argument cellTypeMap is a bitmap of cell type labels.
% Argument hint is an SYDictionary instance holding parameters.
% Parameters:
%   CPMLabelDictKey: (SYDictionary) dictinary of labels refered in
%   subclass.
    if nargin < 5
        return;
    end
    
    obj.cellMap = uint16(cellMap);
    obj.subcellMap = uint16(subcellMap);
    obj.cellTypeMap = uint16(cellTypeMap);
    
    obj.hint = hint.copy;
    obj.labelDict = hint.objectForKey("CPMLabelDictKey");
    
    obj.frameSize = [size(cellMap,1),size(cellMap,2)];
    
    rc = round(obj.frameSize(1) / 2);
    cc = round(obj.frameSize(2) / 2);
    obj.dmCenter = [rc,cc];
    bitmap = zeros(obj.frameSize(1),obj.frameSize(2),2);
    bitmap(1:rc,:,1) = repmat((rc - 1:-1:0)',1,obj.frameSize(2));
    bitmap(rc + 1:end,:,1) = ...
        repmat((1:obj.frameSize(1) - rc)',1,obj.frameSize(2));
    bitmap(:,1:cc,2) = repmat(cc - 1:-1:0,obj.frameSize(1),1);
    bitmap(:,cc + 1:end,2) = ...
        repmat(1:obj.frameSize(2) - cc,obj.frameSize(1),1);
    bitmap = sqrt(sum(bitmap .^ 2,3));
    obj.distanceMap = bitmap;
end
function obj = initWithData(obj,data)
% Initialization method with SYData instance.
% obj = initWithData(obj,data)
    s = data.var;
    cellMap_ = s.cellMap_;
    subcellMap_ = s.subcellMap_;
    cellTypeMap_ = s.cellTypeMap_;
    data_ = SYData(s.hint_var);
    hint_ = SYDictionary;
    hint_.initWithData(data_);
    
    obj = obj.initWithMaps(cellMap_,subcellMap_,cellTypeMap_,hint_);
end
function dest = copy(obj,dest)
% Method to make a shallow copy.
% dest = copy(obj,dest)
% Return value is an instance of the class.
    if nargin < 2
        dest = CPMMap;
    end
    copy@SYObject(obj,dest);
    
    dest.cellMap = obj.cellMap;
    dest.subcellMap = obj.subcellMap;
    dest.cellTypeMap = obj.cellTypeMap;
    dest.hint = obj.hint;
    dest.labelDict = obj.labelDict;
    dest.frameSize = obj.frameSize;
    dest.distanceMap = obj.distanceMap;
    dest.dmCenter = obj.dmCenter;
    dest.pointUpdated = obj.pointUpdated;
    dest.labelUpdated = obj.labelUpdated;
end

function delete(obj)
% Method to clean the instance.
% delete(obj)
    delete@SYObject(obj);
end

function result = data(obj)
% Method to convert the instance to an SYData instance.
% result = data(obj)
% Return value is an SYData instance.
    s.cellMap_ = obj.cellMap;
    s.subcellMap_ = obj.subcellMap;
    s.cellTypeMap_ = obj.cellTypeMap;
    data = obj.hint.data;
    s.hint_var = data.var;
    
    data = SYData(s);
    result = data;
end

function drawCellAtPoint(obj,bitmap,point)
% Method to draw a cell into the map.
% drawCellAtPoint(obj,bitmap,point)
% Argument bitmap is a stack of images for cell label, subcellular
%   components label, and cell type label.
% Argument point specifies an origin in the map where the cell is drawn.
    mask = bitmap(:,:,1) > 0;
    array = unique(obj.cellMap(:));
    array(array < 1) = [];
    array = setdiff(1:length(array) + 1,array);
    
    indices = (1:size(bitmap,1)) + point(1) - 1;
    jndices = (1:size(bitmap,2)) + point(2) - 1;
    kndices = obj.indices(indices,jndices);
    
    % draw cell.
    citmap = obj.cellMap(kndices);
    citmap = citmap .* uint16(~mask) + uint16(array(1) .* mask);
    obj.cellMap(kndices) = citmap;
    
    % draw subcellular components.
    citmap = obj.subcellMap(kndices);
    citmap = citmap .* uint16(~mask) + bitmap(:,:,2) .* uint16(mask);
    obj.subcellMap(kndices) = citmap;
    
    % draw cell type.
    citmap = obj.cellTypeMap(kndices);
    citmap = citmap .* uint16(~mask) + bitmap(:,:,3) .* uint16(mask);
    obj.cellTypeMap(kndices) = citmap;
end

function obj = updatePixel(obj,pixel,label)
% Method to update a pixel.
% obj = updatePixel(obj,pixel,label)
% Argument pixel is a point to be updated.
% Argument label is an array of cell label, subcellular component label,
%    and cell type label.
    pixel = obj.indices(pixel(1),pixel(2));
    obj.pointUpdated = pixel;
    obj.labelUpdated = [obj.cellMap(pixel), ...
        obj.subcellMap(pixel), ...
        obj.cellTypeMap(pixel)];
    
    obj.cellMap(pixel) = label(1);
    obj.subcellMap(pixel) = label(2);
    obj.cellTypeMap(pixel) = label(3);
end
function obj = revert(obj)
% Method to reject an update.
% obj = revert(obj)
    pixel = obj.pointUpdated;
    label = obj.labelUpdated;
    obj.cellMap(pixel) = label(1);
    obj.subcellMap(pixel) = label(2);
    obj.cellTypeMap(pixel) = label(3);
end

function result = indices(obj,rowIndices,columnIndices)
% Method to convert  row and column indices to a matrix of linear index.
% result = indices(obj,rowIndices,columnIndices)
% Argument rowIndices is an array of index for the row.
% Argument columnIndices is an array of index for the column.
% Return value is a matrix whose elements are index to the point specified
%   by the rowIndices and columnIndices.
    result = obj.indicesInTorus(rowIndices,columnIndices);
end
function result = indicesInTorus(obj,rowIndices,columnIndices)
% Method to convert row and column indices to a matrix of linear index.
% result = indicesInTorus(obj,rowIndices,columnIndices)
% Argument rowIndices is an array of index for the row.
% Argument columnIndices is an array of index for the column.
% Return value is a matrix whose elements are index to the point specified
%   by the rowIndices and columnIndices, assuming the indices a modulus and
%   a space of the map is torus.
    m = zeros(obj.frameSize(1),obj.frameSize(2));
    m(:) = 1:length(m(:));
    
    rowIndices = mod(rowIndices - 1,obj.frameSize(1)) + 1;
    columnIndices = mod(columnIndices - 1,obj.frameSize(2)) + 1;
    n = m(rowIndices,columnIndices);
    
    result = n;
end
function result = indicesFromSub(obj,array)
% Method to convert subscripts to linear indices.
% result = indicesFromSub(obj,array)
% Argument array is a {n,2} array of points, where the first column
%   specifies the row and the second column specifies the column.
% Return value is an array of linear index.
    result = obj.indicesFromSubInTorus(array);
end
function result = indicesFromSubInTorus(obj,array)
% Method to convert subscripts to linear indices.
% result = indicesFromSubInTorus(obj,array)
% Argument array is a {n,2} array of points, where the first column
%   specifies the row and the second column specifies the column.
% Return value is an array of linear index, assuming the indices a modulus 
%   and a space of the map is torus.
    array(:,1) = mod(array(:,1) - 1,obj.frameSize(1)) + 1;
    array(:,2) = mod(array(:,2) - 1,obj.frameSize(2)) + 1;
    
    result = sub2ind(obj.frameSize,array(:,1),array(:,2));
end
function result = subscrptsFromInd(obj,array)
% Method to convert index to subscripts.
% result = subscrptsFromInd(obj,array)
% Argument array is an array of index.
% Return value is a {n,2} array of points, where the first column
%   specifies the row and the second column specifies the column.
    result = obj.subscriptsFromIndInTorus(array);
end
function result = subscriptsFromIndInTorus(obj,array)
% Method to convert index to subscripts.
% result = subscriptsFromIndInTorus(obj,array)
% Argument array is an array of index.
% Return value is a {n,2} array of points, where the first column
%   specifies the row and the second column specifies the column, assuming
%   the indices a modulus and a space of the map is torus.
    [r,c] = ind2sub(obj.frameSize,array(:));
    result = [r,c];
end

function result = labelOf(obj,str)
% Method to get labels for a category.
% result = labelOf(obj,str)
% Argument str is a key string specifying the category.
% Return value is an array of labels in the category.
    result = obj.labelDict.objectForKey(str);
end
function result = mapOfLabel(obj,label)
% Method to get a mask for labels.
% result = mapOfLabel(obj,label)
% Argument label is an array of labels.
% Return value is a mask covering the labels.
    if ischar(label)
        label = obj.labelDict.objectForKey(label);
    end
    a = zeros(1,1,length(label(:)));
    a(:) = label(:);
    bitmap = any(obj.subcellMap == a,3);
    
    result = bitmap;
end

function result = distanceMapAround(obj,point)
% Method to get a map showing distance from a point.
% result = distanceMapAround(obj,point)
% Argument point specifies an origin.
% Return value is a map of distance from the origin.
    result = obj.distanceMapOnTorusAround(point);
end
function result = distanceMapOnClosedPlaneAround(obj,point)
    d_r = (1:obj.frameSize(1))' - point(1);
    d_c = (1:obj.frameSize(2)) - point(2);
    result = sqrt(d_r .^ 2 + d_c .^ 2);
end
function result = distanceMapOnTorusAround(obj,point)
% Method to get a map showing distance from a point.
% result = distanceMapOnTorusAround(obj,point)
% Argument point specifies an origin.
% Return value is a map of distance from the origin, assuming the indices
%   a modulus and a space of the map is torus.
    u = obj.dmCenter(1) - round(point(1));
    l = obj.dmCenter(2) - round(point(2));
    indices = mod(u:u + obj.frameSize(1) - 1,obj.frameSize(1)) + 1;
    jndices = mod(l:l + obj.frameSize(2) - 1,obj.frameSize(2)) + 1;
    
    result = obj.distanceMap(indices,jndices);
end
function result = distanceMapFromComponent(obj,component,isPlaneClosed)
    if isa(component,"string")
        labelArray = obj.hint.objectForKey("CPMLabelArrayKey");
        subcellIndices = obj.hint.objectForKey("CPMSubcellIndicesKey");
        component = ...
            subcellIndices.var(labelArray.indexOfObject(component));
    end

    ditmap = [];
    
    citmap = obj.subcellMap == component;
    indices = find(citmap);
    for i = indices'
        [r,c] = ind2sub(obj.frameSize,i);
        point = [r,c];
        if isPlaneClosed
            bitmap = obj.distanceMapOnClosedPlaneAround(point);
        else
            bitmap = obj.distanceMapAround(point);
        end
        ditmap = min(cat(3,ditmap,bitmap),[],3);
    end

    result = ditmap;
end

end
end
