% CPM Class CPMAAPOCell < CPMPCPCell.
% Written by Satoshi Yamashita.
% Data class representing a cell with adjustable Ao and Po values.

classdef CPMAAPOCell < CPMPCPCell
properties
    
end

methods
function obj = CPMAAPOCell

end

function dest = copy(obj,dest)
% Method to make a shallow copy.
% dest = copy(obj,dest)
% Return value is an instance of the class.
    if nargin < 2
        dest = CPMAAPOCell;
    end
    copy@CPMPCPCell(obj,dest);
end

function makePoCurrentLength(obj)
% Method to initialize reference lengths of cell surface.
% makePoCurrentLength(obj)
% Po are calculated with a polynomial with respect to current lengths.
% Constants and coefficients are stored in PoProjector/PoProjectorArray.
    makePoCurrentLength@CPMPCPCell(obj);
    
    data = obj.owner.PoProjectorArray.objectAtIndex(obj.cellType);
    for i = 1:length(obj.Po)
        p = 0;
        for j = 1:size(data.var,2)
            p = p + data.var(i,j) * obj.Po(i) ^ (j - 1);
        end
        obj.Po(i) = p;
    end
end
function makeAoCurrentVolume(obj)
% Method to initialize reference volume of the cell and subcellular
%   components in it.
% makeAoCurrentVolume(obj)
% Ao are calculated with a polynomial with respect to current lengths.
% Constants and coefficients are stored in AoProjector/AoProjectorArray.
    makeAoCurrentVolume@CPMPCPCell(obj);
    
    data = obj.owner.AoProjectorArray.objectAtIndex(obj.cellType);
    for i = 1:length(obj.Ao)
        a = 0;
        for j = 1:size(data.var,2)
            a = a + data.var(i,j) * obj.Ao(i) ^ (j - 1);
        end
        obj.Ao(i) = a;
    end
end

end
end
