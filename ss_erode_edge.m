% Short script to make 2 pixel width edge of Potts cells to 1 pixel width.
% Argument map is a CPMMap instance.
% Argument cells is an SYArray instance containing cells.
% Return value is a 1/0 float bitamp data.

function result = ss_erode_edge(map,cells)
    % Weight edge pixels.
    witmap = lf_weight_edge(map,cells);
    witmap = lf_dig_edge(map,cells,witmap);
    
    % Check connectedness.
    citmap = lf_channel_basin(map,cells);
    witmap(citmap) = 0;
    
    % Erode edge.
    witmap = lf_erode_edge(map,witmap);
    
    result = double(~witmap);
end

function result = lf_dilate_binary(bitmap)
    if isnumeric(bitmap)
        bitmap = bitmap > 0;
    end
    
    bitmap_ = bitmap;
    bitmap_(1:end - 1,:) = bitmap_(1:end - 1,:) | bitmap(2:end,:);
    bitmap_(:,1:end - 1) = bitmap_(:,1:end - 1) | bitmap(:,2:end);
    bitmap_(2:end,:) = bitmap_(2:end,:) | bitmap(1:end - 1,:);
    bitmap_(:,2:end) = bitmap_(:,2:end) | bitmap(:,1:end - 1);
    
    result = bitmap_;
end

function result = lf_weight_edge(map,cells)
    % Measure distance from basin to edge.
    ditmap = zeros(map.frameSize);
    for i = 1:cells.count
        cel = cells.objectAtIndex(i);
        siz = cel.frameSize - 2;
        rim = zeros(siz,'logical');
        rim(cel.cellRim) = true;
        indices = cel.frame(1) + 1:cel.frame(2) - 1;
        jndices = cel.frame(3) + 1:cel.frame(4) - 1;
        kndices = map.indices(indices,jndices);
        citmap = map.cellMap(kndices) == cel.label;
        inner = citmap;
        inner(rim) = false;
        scope = ditmap(kndices);
        d = 1;
        while true
            inner = lf_dilate_binary(inner);
            inner = inner & citmap;
            mask = rim & inner;
            if ~any(mask(:))
                scope(rim) = d;
                break;
            end
            scope(mask) = d;
            d = d + 1;
            rim(mask) = false;
        end
        ditmap(kndices) = scope;
    end
    
    result = ditmap;
end
function result = lf_dig_edge(map,cells,witmap)
    % Set pixel value to maximum in its neighbor.
    xitmap = zeros(map.frameSize);
    nask = [2,4,5,6,8]';
    for i = 1:cells.count
        cel = cells.objectAtIndex(i);
        indices = cel.frame(1):cel.frame(2);
        jndices = cel.frame(3):cel.frame(4);
        kndices = map.indices(indices,jndices);
        mask = double(map.cellMap(kndices) == cel.label);
        
        scope = witmap(kndices) .* mask;
        tcope = xitmap(kndices);
        indices = find(scope(:));
        for j = indices'
            [r,c] = ind2sub(size(scope),j);
            neig = scope(r - 1:r + 1,c - 1:c + 1);
            neig = neig(nask);
            tcope(r,c) = max(neig);
        end
        xitmap(kndices) = tcope;
    end
    
    result = xitmap;
end

function result = lf_channel_basin(map,cells)
    % Mark path between disconnected components for each basin if possible.
    chtmap = zeros(map.frameSize,'logical');
    for i = 1:cells.count
        cel = cells.objectAtIndex(i);
        indices = cel.frame(1):cel.frame(2);
        jndices = cel.frame(3):cel.frame(4);
        kndices = map.indices(indices,jndices);
        citmap = map.cellMap(kndices) == cel.label;
        bitmap_ = citmap(2:end - 1,2:end - 1);
        bitmap_(cel.cellRim) = false;
        bitmap = zeros(cel.frameSize,'logical');
        bitmap(2:end - 1,2:end - 1) = bitmap_;
        nitmap = chtmap(kndices);
        
        nitmap = lf_dilate_binary(nitmap);
        
        flag = false;
        while true
            sitmap = ...
                IPConnectedComponents.connectedBinaryComponents(...
                double(bitmap),4);
            if max(sitmap(:)) == 1
                break;
            end
            
            b1 = sitmap == 1;
            bo = sitmap > 1;
            ditmap = double(b1);
            d = 2;
            while true
                mask = ditmap > 0;
                rim = lf_dilate_binary(mask);
                rim = rim & citmap;
                rim(mask) = false;
                rim(nitmap) = false;
                if ~any(rim(:))
                    flag = true;
                    break;
                end
                
                ditmap(rim) = d;
                d = d + 1;
                
                bo_ = bo & rim;
                if any(bo_(:))
                    break;
                end
            end
            if flag
                break;
            end
            
            rp = bo_;
            d = d - 2;
            while true
                rp_ = lf_dilate_binary(rp);
                rp = rp | (rp_ & (ditmap == d));
                d = d - 1;
                if d < 1
                    break;
                end
            end
            bitmap(rp) = true;
        end
        
        scope = chtmap(kndices);
        scope(bitmap) = true;
        chtmap(kndices) = scope;
    end
    
    result = chtmap;
end

function result = lf_erode_edge(map,witmap)
    % Mark edge off in an order of weight as long as it won't make
    % basins adjacent.
    mask = [1,3,7,9];
    while true
        bitmap = witmap == 0;
        bitmap = double(lf_dilate_full(bitmap));
        ritmap = witmap .* bitmap;
        d = max(ritmap(:));
        indices = find(ritmap(:) == d);
        indices = indices(randperm(length(indices)));
        flag = true;
        for i = indices'
            [r,c] = ind2sub(map.frameSize,i);
            indices = r - 1:r + 1;
            jndices = c - 1:c + 1;
            kndices = map.indices(indices,jndices);
            citmap = map.cellMap(kndices);
            citmap = citmap ~= citmap(5);
            citmap(witmap(kndices) > 0) = false;
            citmap(mask) = false;
            if any(citmap(:))
                continue;
            end
            
            witmap(i) = 0;
            flag = false;
        end
        if flag
            if d < 2
                break;
            end
            
            ritmap = witmap .* bitmap;
            indices = ritmap(:) == d;
            witmap(indices) = d - 1;
        end
    end
    
    result = witmap > 0;
end
function result = lf_dilate_full(bitmap)
    bitmap_ = bitmap;
    indices = [size(bitmap,1),1:size(bitmap,1) - 1];
    bitmap_ = bitmap_ | bitmap(indices,:);
    indices = [2:size(bitmap,1),1];
    bitmap_ = bitmap_ | bitmap(indices,:);
    indices = [size(bitmap,2),1:size(bitmap,2) - 1];
    bitmap_ = bitmap_ | bitmap(:,indices);
    indices = [2:size(bitmap,2),1];
    bitmap_ = bitmap_ | bitmap(:,indices);
    
    result = bitmap_;
end
