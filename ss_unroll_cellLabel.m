% Short script to draw cell labels in a frame.
% Argument map is a CPMMap instance.
% Argument cells is an SYArray instance containing cells.
% Argument labels is an array of cell labels.
% Argument upperCells is a boolean array marking cells in upper half.
% Argument leftCells is a boolean array marking cells in left half.
% Argument frame specifies canvas size.
% Return value is a bitmap of cell labels.

function result = ss_unroll_cellLabel(map,cells,labels,upperCells, ...
    leftCells,frame)
    citmap = zeros(frame,'uint16');
    
    u = round((frame(1) - map.frameSize(1)) / 2);
    l = round((frame(2) - map.frameSize(2)) / 2);
    for i = 1:cells.count
        cel = cells.objectAtIndex(i);
        index = labels == cel.label;
        
        indices = cel.frame(1) + 1:cel.frame(2) - 1;
        jndices = cel.frame(3) + 1:cel.frame(4) - 1;
        kndices = map.indices(indices,jndices);
        mask = map.cellMap(kndices) == cel.label;
        
        if indices(end) > map.frameSize(1) && upperCells(index)
            indices = indices - map.frameSize(1);
        end
        indices = indices + u;
        if jndices(end) > map.frameSize(2) && leftCells(index)
            jndices = jndices - map.frameSize(2);
        end
        jndices = jndices + l;
        
        kndices = indices < frame(1) & indices > 1;
        lndices = jndices < frame(2) & jndices > 1;
        indices = indices(kndices);
        jndices = jndices(lndices);
        mask = mask(kndices,lndices);
        
        scope = citmap(indices,jndices);
        scope(mask) = cel.label;
        citmap(indices,jndices) = scope;
    end
    
    result = citmap;
end
