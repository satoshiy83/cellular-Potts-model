% CPM Class CPMConverterMatrix < SYObject.
% Written by Satoshi Yamashita.
% Data class for label conversion.

classdef CPMConverterMatrix < SYObject
properties
    labelArray = nan; % SYArray.
    converterMatrix = nan; % double[n,n,2];
end

methods
function obj = initWithData(obj,data)
% Initialization method with an SYData instance.
% obj = initWithData(obj,data)
    s = data.var;
    array = SYArray;
    array.initWithData(SYData(s.labelArray_));
    obj.labelArray = array;
    obj.converterMatrix = s.converterMatrix_;
end
function result = data(obj)
% Method to convert the instance to an SYData instance.
% result = data(obj)
% Return value is an SYData instance.
    s.labelArray_ = obj.labelArray.data.var;
    s.converterMatrix_ = obj.converterMatrix;
    data = SYData(s);
    result = data;
end

function result = converterWithLabels(obj,labels)
% Method to get an label converter.
% result = converterWithLabels(obj,labels)
% Argument labels is an array of subcellular components.
% Return value is a matrix defining conversion of labels.
    indices = zeros(labels.count,1);
    for i = 1:labels.count
        indices(i) = obj.labelArray.indexOfObject(labels.objectAtIndex(i));
    end
    result = SYData(obj.converterMatrix(indices,indices,:));
end

function result = matrixWithLabels(obj,labels)
% depreted method.
    indices = zeros(labels.count,1);
    for i = 1:labels.count
        indices(i) = obj.labelArray.indexOfObject(labels.objectAtIndex(i));
    end
    result = obj.converterMatrix(indices,indices,:);
end

function result = indexOfLabel(obj,label)
% depreted method.
    result = obj.labelArray.indexOfObject(label);
end

end
end
