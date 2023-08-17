% CPM Class CPMepithelialCell < CPMPolarizedCell.
% Written by Satoshi Yamashita.
% Data class representing an epithelial cell with apical-lateral adherens
% junction.


classdef CPMEpithelialCell < CPMPolarizedCell
properties
    % Shared property.
    epitheliaSwitch = nan; % SYData.
    
    % Instance variables.
    AJPixels = nan; % double[n]. Indices in ex-frame.
end

methods

function obj = initWithData(obj,data)
% Initialization method with an SYData instance.
% obj = initWithData(obj,data)
    initWithData@CPMPolarizedCell(obj,data);
    
    s = data.var;
    obj.AJPixels = s.AJPixels;
end
function dest = copy(obj,dest)
% Method to make a shallow copy.
% dest = copy(obj,dest)
% Return value is an instance of the class.
    if nargin < 2
        dest = CPMEpithelialCell;
    end
    copy@CPMPolarizedCell(obj,dest);

    dest.epitheliaSwitch = obj.epitheliaSwitch;
    dest.AJPixels = obj.AJPixels;
end

function importFromOwner(obj)
% Method method to import shared variables from the owner.
% importFromOwner(obj)
    importFromOwner@CPMPolarizedCell(obj);
    
    obj.epitheliaSwitch = obj.owner.epitheliaSwitch;
end

function result = data(obj)
% Method to convert the instance to an SYData instance.
% result = data(obj)
% Return value is an SYData instance.
    data_ = data@CPMPolarizedCell(obj);
    s = data_.var;
    
    s.AJPixels = obj.AJPixels;
    
    data_ = SYData(s);
    result = data_;
end

function markAdherensJunction(obj)
% Method to mark adherens junction sites.
% markAdherensJunction(obj)
% Lateral sites are marked as adherens junction sites when it is adjacent
% to an apical site of any cell and adjacent to a lateral site of different
% cell.
    n = find(any(obj.epitheliaSwitch.var == obj.cellType,1),1);
    if n == 1 % not epithelial cell.
        return;
    end
    
    array = zeros(length(obj.lateralRim),1);
    c = 0;
    
    point = [0,0];
    m = [-1,-1; -1,0; -1,1; 0,-1; 0,1; 1,-1; 1,0; 1,1;];
    for i = 1:length(obj.lateralRim)
        [point(1),point(2)] = ind2sub(obj.map.frameSize,obj.lateralRim(i));
        flags = [false,false];
        for j = 1:size(m,1)
            qoint = point + m(j,:);
            qoint = obj.map.indicesFromSub(qoint);
            
            label = obj.map.cellMap(qoint);
            pc = obj.owner.cellWithLabel(label);
            % Check if adjacent to apical pixel.
            if ~flags(1) && pc.isPixelApical(qoint)
                flags(1) = true;
            end
            
            % Check if adjacent to lateral pixel of different cell.
            if ~flags(2) && label ~= obj.label && pc.isPixelLateral(qoint)
                flags(2) = true;
            end
            
            if all(flags)
                c = c + 1;
                array(c) = obj.lateralRim(i);
                break;
            end
        end
    end
    
    if c > 0
        obj.AJPixels = array(1:c);
    else
        obj.AJPixels = [];
    end
end

function result = imageOfSubcellularLocations(obj)
% Method to get an image of subcellular locations.
% result = imageOfSubcellularLocations(obj)
% Return value is an image of subcellular locations of the cell.
    siz = obj.frameSize - 2;
    bitmap = imageOfSubcellularLocations@CPMPolarizedCell(obj);
    
    n = find(any(obj.drawingSwitch.var == obj.cellType,1),1);
    if n == 1 % media
        result = bitmap;
        return;
    end
    
    n = find(any(obj.epitheliaSwitch.var == obj.cellType,1),1);
    if n == 2 % epithelial cell.
        % adherens junction.
        color = obj.lut.var(end - 4,:);
        mask = zeros(siz(1),siz(2),'logical');
        indices = obj.ind2indInFrame(obj.AJPixels);
        mask(indices) = true;
        bitmap(mask(:,:,[1,1,1])) = repmat(color,sum(mask(:)),1);
    end
    
    result = bitmap;
end

end
end
