% CPM Class CPMCACell < CPMCell.
% Written by Satoshi Yamashita.
% Subclass of CPMCell enable to import reference volume and surface area
% from an array prepared by ss_export_AoPo().

classdef CPMCACell < CPMCell
properties
end

methods

function importAoPo(obj,labels,array)
% Method to import reference volume and surface area from an array.
% importAoPo(obj,labels,array)
% Argument labels is an array of cell labels.
% Argument array is an SYArray instance carrying the reference values
%   prepared by ss_export_AoPo().
    index = find(labels == obj.label);
    if isempty(index)
        disp('Unable to find AoPo.');
        return;
    end
    AoPo = array.objectAtIndex(index);
    Ao = AoPo.objectAtIndex(1);
    obj.Ao = Ao.var;
    Po = AoPo.objectAtIndex(2);
    obj.Po = Po.var;
    
    obj.registerEnergy;
end

end
end
