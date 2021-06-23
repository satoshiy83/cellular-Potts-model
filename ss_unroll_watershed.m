% Short script to draw segmented cells in a frame.
% Argument map is a CPMMap instance.
% Argument cells is an SYArray instance containing cells.
% Argument labels is an array of cell labels.
% Argument upperCells is a boolean array marking cells in upper half.
% Argument leftCells is a boolean array marking cells in left half.
% Argument watershed is a bitmap of watershed.
% Argument frame specifies canvas size.
% Return value is a bitmap of watershed.

function result = ss_unroll_watershed(map,cells,labels,upperCells, ...
    leftCells, watershed,frame)
    witmap = ones(frame);
    citmap = zeros(frame,'uint16');
    
    u = round((frame(1) - map.frameSize(1)) / 2);
    l = round((frame(2) - map.frameSize(2)) / 2);
    for i = 1:cells.count
        cel = cells.objectAtIndex(i);
        index = labels == cel.label;
        
        indices = cel.frame(1) + 1:cel.frame(2) - 1;
        jndices = cel.frame(3) + 1:cel.frame(4) - 1;
        kndices = map.indices(indices,jndices);
        bitmap = watershed(kndices);
        mask = map.cellMap(kndices) == cel.label;
        bitmap = bitmap .* double(mask);
        
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
        bitmap = bitmap(kndices,lndices);
        mask = mask(kndices,lndices);
        
        scope = witmap(indices,jndices);
        scope(mask) = 0;
        scope = scope + bitmap;
        witmap(indices,jndices) = scope;
        
        scope = citmap(indices,jndices);
        scope(mask) = cel.label;
        citmap(indices,jndices) = scope;
    end
    
    mask = citmap > 0;
    mask_ = mask;
    mask_(1:end - 1,:) = mask_(1:end - 1,:) | mask(2:end,:);
    mask_(:,1:end - 1) = mask_(:,1:end - 1) | mask(:,2:end);
    mask_(2:end,:) = mask_(2:end,:) | mask(1:end - 1,:);
    mask_(:,2:end) = mask_(:,2:end) | mask(:,1:end - 1);
    rim = mask_ & ~mask;
    
    indices = find(rim(:));
    jndices = [2,4,6,8]';
    for i = indices'
        [r,c] = ind2sub(frame,i);
        kndices = r - 1:r + 1;
        kndices(kndices < 1) = 1; kndices(kndices > frame(1)) = frame(1);
        lndices = c - 1:c + 1;
        lndices(lndices < 1) = 1; lndices(lndices > frame(2)) = frame(2);
        scope = mask(kndices,lndices) & logical(witmap(kndices,lndices));
        if any(scope(jndices))
            witmap(i) = 0;
        end
    end
    
    result = witmap;
end
