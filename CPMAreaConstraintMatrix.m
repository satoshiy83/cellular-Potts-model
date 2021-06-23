% CPM Class CPMAreaConstraintMatrix < SYObject.
% Written by Satoshi Yamashita.
% Data class for area vulk modulus.

classdef CPMAreaConstraintMatrix < SYObject
properties
    labelArray = nan; % SYArray.
    
    celltypeArray = nan; % SYArray.
    lambdaArray = nan; % SYArray.
end

methods
function obj = initWithData(obj,data)
% Initialization method with an SYData instance.
% obj = initWithData(obj,data)
    s = data.var;
    array = SYArray;
    array.initWithData(SYData(s.labelArray_));
    obj.labelArray = array;
    array = SYArray;
    array.initWithData(SYData(s.celltypeArray_));
    obj.celltypeArray = array;
    array = SYArray;
    array.initWithData(SYData(s.lambdaArray_));
    obj.lambdaArray = array;
end
function result = data(obj)
% Method to convert the instance to an SYData instance.
% result = data(obj)
% Return value is an SYData instance.
    s.labelArray_ = obj.labelArray.data.var;
    s.celltypeArray_ = obj.celltypeArray.data.var;
    s.lambdaArray_ = obj.lambdaArray.data.var;
    data = SYData(s);
    result = data;
end

function result = ...
        areaVModulusArrayWithCellTypesAndLabels(obj,cellTypes,labels)
% Method to get an array of area vulk modulus.
% result = areaVModulusArrayWithCellTypesAndLabels(obj,cellTypes,labels)
% Argument cellTypes is an array of cell types.
% Argument labels is an array of subcellular components.
% Return value is an array of the vulk modulus of cell and subcellular
%   components for the cell types.
    indices = zeros(1,labels.count);
    for i = 1:labels.count
        index = obj.labelArray.indexOfObject(labels.objectAtIndex(i));
        indices(i) = index + 1;
    end
    if length(indices) ~= labels.count
        disp('Labels mismatch.');
        result = nan;
        return;
    end
    
    array = SYArray;
    for i = 1:cellTypes.count
        cellType = cellTypes.objectAtIndex(i);
        index = obj.celltypeArray.indexOfObject(cellType);
        if ~isempty(index)
            lambda = obj.lambdaArray.objectAtIndex(index);
            array.addObject(SYData(lambda([1,indices])));
        else
            array.addObject(SYData(zeros(labels.count + 1,1)));
        end
    end
    result = array;
end

function result = aLambdaArrayWithCellTypesAndLabels(obj,cellTypes,labels)
% depreted method.
    indices = zeros(1,labels.count);
    for i = 1:labels.count
        index = obj.labelArray.indexOfObject(labels.objectAtIndex(i));
        indices(i) = index + 1;
    end
    if length(indices) ~= labels.count
        disp('Labels mismatch.');
        result = nan;
        return;
    end
    
    array = SYArray;
    for i = 1:cellTypes.count
        cellType = cellTypes.objectAtIndex(i);
        index = obj.celltypeArray.indexOfObject(cellType);
        if ~isempty(index)
            lambda = obj.lambdaArray.objectAtIndex(index);
            array.addObject(lambda([1,indices]));
        else
            array.addObject(zeros(labels.count + 1,1));
        end
    end
    result = array;
end

end
end
