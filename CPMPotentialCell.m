% CPM Class CPMPotentialCell < CPMAAPOCell.
% Written by Satoshi Yamashita.
% Data class representing a cell with potential energy.

classdef CPMPotentialCell < CPMAAPOCell
properties
    % Shared property.
    potentialComponentArray = nan; % SYArray.
    pMatrixArray = nan; % SYArray.
end

methods

function obj = initWithData(obj,data)
% Initialization method with an SYData instance.
% obj = initWithData(obj,data)
    initWithData@CPMAAPOCell(obj,data);
    
    s = data.var;
    array = SYArray;
    array.initWithData(SYData(s.pMatrixArray_));
    obj.pMatrixArray = array;
end

function dest = copy(obj,dest)
% Method to make a shallow copy.
% dest = copy(obj,dest)
% Return value is an instance of the class.
    if nargin < 2
        dest = CPMPotentialCell;
    end
    copy@CPMAAPOCell(obj,dest);
    
    dest.potentialComponentArray = obj.potentialComponentArray;
    dest.pMatrixArray = obj.pMatrixArray;
end

function importFromOwner(obj)
% Method method to import shared variables from the owner.
% importFromOwner(obj)
    importFromOwner@CPMAAPOCell(obj);
    
    obj.potentialComponentArray = obj.owner.potentialComponentArray;
end

function result = data(obj)
% Method to convert the instance to an SYData instance.
% result = data(obj)
% Return value is an SYData instance.
    data_ = data@CPMAAPOCell(obj);
    s = data_.var;
    
    eata_ = obj.pMatrixArray.data;
    s.pMatrixArray_ = eata_.var;
    
    data_ = SYData(s);
    result = data_;
end

function readCellType(obj)
% Method to import cell type from map.
% readCellType(obj)
    readCellType@CPMAAPOCell(obj);
    
    obj.pMatrixArray = ...
        obj.potentialComponentArray.objectAtIndex(obj.cellType);
end

function result = getHamiltonian(obj)
% Method to get energy of the cell.
% result = getHamiltonian(obj)
% Return value is the energy of the cell.
    n = find(any(obj.HamiltonianSwitch.var == obj.cellType,1),1);
    if n == 1 % media.
        result = 0;
    else % cell.
        H = getHamiltonian@CPMAAPOCell(obj);
        
        H = H + obj.potentialEnergy;
        
        result = H;
    end
end
function result = potentialEnergy(obj)
% Method to get potential energy of the cell.
% result = potentialEnergy(obj)
% Return value is the potential energy, where the potential energy was
% averaged among sites labeled with a subcellular component or location
% assinged the potential energy.
    if obj.pMatrixArray.count < 2
        result = 0;
        return;
    end
    
    H = 0;
    
    indices = obj.frame(1) + 1:obj.frame(2) - 1;
    jndices = obj.frame(3) + 1:obj.frame(4) - 1;
    kndices = obj.map.indices(indices,jndices);
    bitmap = obj.map.cellMap(kndices);
    mask = bitmap == obj.label;
    bitmap = double(obj.map.subcellMap(kndices));
    citmap = bitmap .* mask;
    ditmap = obj.map.potentialMap(kndices);
    if obj.map.potentialMapsCount > 1
        l = obj.map.frameSize(1) * obj.map.frameSize(2);
        kndices = kndices + l;
        for i = 2:obj.map.potentialMapsCount
            ditmap = cat(3,ditmap,obj.map.potentialMap(kndices));
        end
    end
    
    % Potential energy of subcellular components.
    pMatrix = obj.pMatrixArray.objectAtIndex(1);
    if ~isempty(pMatrix.var)
        for i = 1:size(pMatrix.var,1)
            bitmap(:) = 0;
            label = pMatrix.var(i,1);
            index = pMatrix.var(i,2);
            nask = citmap == label;
            bitmap(nask) = 1;
            bitmap = bitmap .* ditmap(:,:,index);
            
            n = sum(nask(:));
            if n > 0
                H = H + sum(bitmap(:)) / n;
            end
        end
    end
    
    % Potential energy of subcellular locations.
    pMatrix = obj.pMatrixArray.objectAtIndex(2);
    if ~isempty(pMatrix.var)
        for i = 1:size(pMatrix.var,1)
            switch pMatrix.var(i,1)
                case 1 % cell body.
                    indices = find(mask);
                case 2 % cell perimeter.
                    indices = obj.cellRim;
                case 3 % apical perimeter.
                    indices = obj.ind2indInFrame(obj.apicalRim);
                case 4 % lateral perimeter.
                    indices = obj.ind2indInFrame(obj.lateralRim);
                case 5 % basal perimeter.
                    indices = obj.ind2indInFrame(obj.basalRim);
                case 6 % adherens junction.
                    indices = obj.ind2indInFrame(obj.AJPixels);
                case 7 % upstream perimeter.
                    indices = obj.ind2indInFrame(obj.upstreamRim);
                case 8 % downstream perimeter.
                    indices = obj.ind2indInFrame(obj.downstreamRim);
                case 9 % upstream adherens junction.
                    indices = obj.ind2indInFrame(obj.upstreamAJPixels);
                case 10 % downstream adherens junction.
                    indices = obj.ind2indInFrame(obj.downstreamAJPixels);
                otherwise
                    continue;
            end
            
            n = length(indices);
            if n > 0
                bitmap(:) = 0;
                index = pMatrix.var(i,2);
                bitmap(indices) = 1;
                bitmap = bitmap .* ditmap(:,:,index);
            
                H = H + sum(bitmap(:)) / n;
            end
        end
    end
    
    result = H;
end

end
end
