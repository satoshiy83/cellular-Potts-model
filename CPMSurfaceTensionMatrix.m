% CPM Class CPMSurfaceTensionMatrix < SYObject.
% Written by Satoshi Yamashita.
% Data class for surface elasticity and contractility.

classdef CPMSurfaceTensionMatrix < SYObject
properties
    labelArray = nan; % SYArray.
    celltypeArray = nan; % SYArray.
    
    contactEnergyMatrix = nan; % double[n,n,2].
    baseTension = nan; % double[2].
    
    surfaceYModulusArray = nan; % SYArray.
    
%     pContrArray = nan; % SYArray.
%     pLambdaArray = nan; % SYArray.
    
    contactEnergyGradientMatrix = nan; % double[n,m].
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
    obj.contactEnergyMatrix = s.contactEnergyMatrix_;
    obj.baseTension = s.baseTension_;
    array = SYArray;
    obj.surfaceYModulusArray = ...
        array.initWithData(SYData(s.surfaceYModulusArray_));
%     array = SYArray;
%     obj.pContrArray = array.initWithData(SYData(s.pContrArray_));
%     array = SYArray;
%     obj.pLambdaArray = array.initWithData(SYData(s.pLambdaArray_));
    obj.contactEnergyGradientMatrix = s.contactEnergyGradientMatrix_;
end
function result = data(obj)
% Method to convert the instance to an SYData instance.
% result = data(obj)
% Return value is an SYData instance.
    s.labelArray_ = obj.labelArray.data.var;
    s.celltypeArray_ = obj.celltypeArray.data.var;
    s.contactEnergyMatrix_ = obj.contactEnergyMatrix;
    s.baseTension_ = obj.baseTension;
    s.surfaceYModulusArray_ = obj.surfaceYModulusArray.data.var;
%     s.pContrArray_ = obj.pContrArray.data.var;
%     s.pLambdaArray_ = obj.pLambdaArray.data.var;
    s.contactEnergyGradientMatrix_ = obj.contactEnergyGradientMatrix;
    
    data = SYData(s);
    result = data;
end

function result = contactEnergyWithLabels(obj,labels)
% Method to get a matrix of contact energy.
% result = contactEnergyWithLabels(obj,labels)
% Argument labels is an array of subcellular components.
% Return value is a matrix of contact energy between the subcellular
%   components.
    m = zeros(labels.count,labels.count,2);
    m(:,:,1) = obj.baseTension(1);
    m(:,:,2) = obj.baseTension(2);
    for i = 1:obj.labelArray.count
        index = labels.indexOfObject(obj.labelArray.objectAtIndex(i));
        if isempty(index)
            continue;
        end
        for j = i:obj.labelArray.count
            jndex = labels.indexOfObject(obj.labelArray.objectAtIndex(j));
            if isempty(jndex)
                continue;
            end
            m(index,jndex,:) = obj.contactEnergyMatrix(i,j,:);
            m(jndex,index,:) = obj.contactEnergyMatrix(i,j,:);
        end
    end
    result = SYData(m);
end
function result = surfaceYModulusArrayWithCellTypes(obj,cellTypes)
% Method to get an array of surface Young's modulus.
% result = surfaceYModulusArrayWithCellTypes(obj,cellTypes)
% Argument cellTypes is an array of cell types.
% Return value is an array of the surface Young's modulus for the cell
%   types.
    array = SYArray;
    for i = 1:cellTypes.count
        index = obj.celltypeArray.indexOfObject(cellTypes.objectAtIndex(i));
        if ~isempty(index)
            array.addObject(obj.surfaceYModulusArray.objectAtIndex(index));
        else
            array.addObject(SYData([]));
        end
    end
    result = array;
end
function result = CEGradArrayWithLabels(obj,labels)
% Method to get a matrix of contact energy for gradient components.
% result = CEGradArrayWithLabels(obj,labels)
% Argument labels is an array of subcellular components.
% Return value is a contact energy between the subcellular components and
%   gradient components.
    m = zeros(labels.count,1,size(obj.contactEnergyGradientMatrix,2));
    for i = 1:labels.count
        label = labels.objectAtIndex(i);
        index = obj.indexOfLabel(label);
        m(i,1,:) = obj.contactEnergyGradientMatrix(index,:);
    end
    result = SYData(m);
end

% function result = matrixWithLabels(obj,labels)
%     m = zeros(labels.count,labels.count,2);
%     m(:,:,1) = obj.baseTension(1);
%     m(:,:,2) = obj.baseTension(2);
%     for i = 1:obj.labelArray.count
%         index = labels.indexOfObject(obj.labelArray.objectAtIndex(i));
%         if isempty(index)
%             continue;
%         end
%         for j = i:obj.labelArray.count
%             jndex = labels.indexOfObject(obj.labelArray.objectAtIndex(j));
%             if isempty(jndex)
%                 continue;
%             end
%             m(index,jndex,:) = obj.contactEnergyMatrix(i,j,:);
%             m(jndex,index,:) = obj.contactEnergyMatrix(i,j,:);
%         end
%     end
%     result = m;
% end
% function result = pContrArrayWithCellTypes(obj,cellTypes)
%     array = SYArray;
%     for i = 1:cellTypes.count
%         index = obj.celltypeArray.indexOfObject(cellTypes.objectAtIndex(i));
%         if ~isempty(index)
%             array.addObject(obj.pContrArray.objectAtIndex(index));
%         else
%             array.addObject([]);
%         end
%     end
%     result = array;
% end
% function result = pLambdaArrayWithCellTypes(obj,cellTypes)
%     array = SYArray;
%     for i = 1:cellTypes.count
%         index = obj.celltypeArray.indexOfObject(cellTypes.objectAtIndex(i));
%         if ~isempty(index)
%             array.addObject(obj.pLambdaArray.objectAtIndex(index));
%         else
%             array.addObject([]);
%         end
%     end
%     result = array;
% end
% function result = coefGradArrayWithLabels(obj,labels)
%     m = zeros(labels.count,1,size(obj.contactEnergyGradientMatrix,2));
%     for i = 1:labels.count
%         label = labels.objectAtIndex(i);
%         index = obj.indexOfLabel(label);
%         m(i,1,:) = obj.contactEnergyGradientMatrix(index,:);
%     end
%     result = m;
% end

function result = indexOfLabel(obj,label)
    result = obj.labelArray.indexOfObject(label);
end

end 
end