% CPM Class CPMPCPCell < CPMEpithelialCell.
% Written by Satoshi Yamashita.
% Data class representing a cell with PCP.

classdef CPMPCPCell < CPMEpithelialCell
properties
    % Shared property.
    PCPSwitch = nan; % SYData.
    
    % Instances variables.
    upstreamRim = nan; % double[n]: subset of lateralRim.
    downstreamRim = nan; % double[n]: subset of lateralRim.
    upstreamAJPixels = nan; % doble[n]: subset of AJPixels.
    downstreamAJPixels = nan; % double[n]: subset of AJPixels.
    ppolarized = nan; % bool.
end

methods
function obj = initWithData(obj,data)
    initWithData@CPMEpithelialCell(obj,data);
% Initialization method with an SYData instance.
% obj = initWithData(obj,data)
    
    s = data.var;
    obj.upstreamRim = s.upstreamRim;
    obj.downstreamRim = s.downstreamRim;
    obj.upstreamAJPixels = s.upstreamAJPixels;
    obj.downstreamAJPixels = s.downstreamAJPixels;
    obj.ppolarized = s.ppolarized;
end

function dest = copy(obj,dest)
% Method to make a shallow copy.
% dest = copy(obj,dest)
% Return value is an instance of the class.
    if nargin < 2
        dest = CPMPCPCell;
    end
    copy@CPMEpithelialCell(obj,dest);

    dest.upstreamRim = obj.upstreamRim;
    dest.downstreamRim = obj.downstreamRim;
    dest.upstreamAJPixels = obj.upstreamAJPixels;
    dest.downstreamAJPixels = obj.downstreamAJPixels;
    dest.ppolarized = obj.ppolarized;
end

function importFromOwner(obj)
% Method method to import shared variables from the owner.
% importFromOwner(obj)
    importFromOwner@CPMEpithelialCell(obj);
    
    obj.PCPSwitch = obj.owner.PCPSwitch;
end

function result = data(obj)
% Method to convert the instance to an SYData instance.
% result = data(obj)
% Return value is an SYData instance.
    data_ = data@CPMEpithelialCell(obj);
    s = data_.var;
    
    s.upstreamRim = obj.upstreamRim;
    s.downstreamRim = obj.downstreamRim;
    s.upstreamAJPixels = obj.upstreamAJPixels;
    s.downstreamAJPixels = obj.downstreamAJPixels;
    s.ppolarized = obj.ppolarized;
    
    data_ = SYData(s);
    result = data_;
end

function makePoCurrentLength(obj)
% Method to set apical, lateral, and basal Po current length.
% makePoCurrentLength(obj)
    makePoCurrentLength@CPMEpithelialCell(obj);
    
    n = find(any(obj.PCPSwitch.var == obj.cellType,1),1);
    if n == 1 || ~obj.ppolarized
        obj.Po = cat(1,obj.Po,0,0);
        return;
    end
    
    indices = obj.frame(1):obj.frame(2);
    jndices = obj.frame(3):obj.frame(4);
    kndices = obj.map.indices(indices,jndices);
    bitmap = obj.map.cellMap(kndices);
    bitmap = bitmap == obj.label;
    
    siz = obj.frameSize - 2;
    citmap = zeros(siz,'uint8');
    indices = obj.ind2indInFrame(obj.upstreamRim);
    citmap(indices) = 1;
    
    m = [-1,-1; -1,0; -1,1; 0,-1; 0,1; 1,-1; 1,0; 1,1;];
    indices = 2:size(bitmap,1) - 1;
    jndices = 2:size(bitmap,2) - 1;
    ditmap = zeros(obj.frameSize,'uint8');
    for i = 1:size(m,1)
        kndices = indices + m(i,1);
        lndices = jndices + m(i,2);
        ditmap(kndices,lndices) = ditmap(kndices,lndices) + citmap;
    end
    ditmap(bitmap) = 0;
    u = sum(ditmap(:));
    
    citmap(:) = 0;
    kndices = obj.ind2indInFrame(obj.downstreamRim);
    citmap(kndices) = 1;
    ditmap(:) = 0;
    for i = 1:size(m,1)
        kndices = indices + m(i,1);
        lndices = jndices + m(i,2);
        ditmap(kndices,lndices) = ditmap(kndices,lndices) + citmap;
    end
    ditmap(bitmap) = 0;
    d = sum(ditmap(:));
    
    obj.Po = cat(1,obj.Po,u,d);
end

function result = polarizeCell(obj)
% Method to polarize cell planarly.
% result = polarizeCell(obj)
    n = find(any(obj.PCPSwitch.var == obj.cellType,1),1);
    if n ~= 2
        obj.ppolarized = false;
        result = false;
        return;
    end
    
    indices = obj.frame(1):obj.frame(2);
    jndices = obj.frame(3):obj.frame(4);
    kndices = obj.map.indices(indices,jndices);
    bitmap = obj.map.cellMap(kndices);
    citmap = obj.map.PCPMap(kndices) == 1;
    
    siz = obj.frameSize - 2;
    ditmap = zeros(siz);
    indices = obj.ind2indInFrame(obj.lateralRim);
    ditmap(indices) = 1;
%     connectedComponents = IPConnectedComponents;
    eitmap = IPConnectedComponents.connectedBinaryComponents(ditmap,8);
    if max(eitmap(:)) == 1
        ditmap = obj.cutPolarizedRim;
        eitmap = IPConnectedComponents.connectedBinaryComponents(ditmap,8);
    end
    
    mask = bitmap == obj.label;
    citmap(mask) = 0;
    fitmap = zeros(siz,'logical');
    m = [-1,-1; -1,0; -1,1; 0,-1; 0,1; 1,-1; 1,0; 1,1];
    indices = 2:obj.frameSize(1) - 1;
    jndices = 2:obj.frameSize(2) - 1;
    for i = 1:size(m,1)
        fitmap = fitmap | citmap(indices + m(i,1),jndices + m(i,2));
    end
    
    for i = 1:max(eitmap(:))
        gitmap = eitmap == i;
        if any(fitmap(:) & gitmap(:))
            ditmap(gitmap) = 2;
        end
    end
    indices = find(ditmap == 1);
    jndices = find(ditmap == 2);
    if length(jndices) < 1
        obj.upstreamRim = nan;
        obj.downstreamRim = nan;
        obj.ppolarized = false;
        result = false;
        return;
    end
    indices = obj.ind2indExFrame(indices);
    jndices = obj.ind2indExFrame(jndices);
    
    result = false;
    if length(indices) ~= length(obj.upstreamRim) || ...
            ~isempty(setdiff(indices,obj.upstreamRim))
        result = true;
    elseif length(jndices) ~= length(obj.downstreamRim) || ...
            ~isempty(setdiff(jndices,obj.downstreamRim))
        result = true;
    end
    
    obj.upstreamRim = indices;
    obj.downstreamRim = jndices;
    obj.ppolarized = true;
    
    obj.upstreamAJPixels = intersect(obj.AJPixels,obj.upstreamRim);
    obj.downstreamAJPixels = intersect(obj.AJPixels,obj.downstreamRim);
    
    obj.map.PCPMap(obj.upstreamRim) = 1;
    obj.map.PCPMap(obj.downstreamRim) = 2;
end
function clearPCP(obj)
% Method to clear PCP labels in the map.
% clearPCP(obj)
    n = find(any(obj.PCPSwitch.var == obj.cellType,1),1);
    if n ~= 2
        return;
    end
    
    indices = obj.frame(1) + 1:obj.frame(2) - 1;
    jndices = obj.frame(3) + 1:obj.frame(4) - 1;
    kndices = obj.map.indices(indices,jndices);
    bitmap = obj.map.cellMap(kndices);
    citmap = obj.map.PCPMap(kndices);
    
    mask = bitmap == obj.label;
    citmap(mask) = 0;
    obj.map.PCPMap(kndices) = citmap;
end
function clearPCPRim(obj)
% Method to clear PCP record in the cell.
% clearPCPRim(obj)
    obj.upstreamRim = nan;
    obj.downstreamRim = nan;
    obj.upstreamAJPixels = nan;
    obj.downstreamAJPixels = nan;
end
function result = checkPCP(obj)
% Method to check PCP marking validity.
% result = checkPCP(obj)
% Return value is boolean.
    n = find(any(obj.PCPSwitch.var == obj.cellType,1),1);
    if n ~= 2
        result = true;
        return;
    end
    
    if ~obj.checkFrame
        disp('break!');
    end
    
    indices = obj.frame(1) + 1:obj.frame(2) - 1;
    jndices = obj.frame(3) + 1:obj.frame(4) - 1;
    kndices = obj.map.indices(indices,jndices);
    mask = obj.map.cellMap(kndices) == obj.label;
    
    bitmap = obj.map.PCPMap(kndices);
    bitmap(~mask) = 0;
    bitmap = bitmap > 0;
    
    siz = obj.frameSize - 2;
    citmap = zeros(siz,'logical');
    indices = obj.ind2indInFrame(obj.lateralRim);
    citmap(indices) = true;
    
    ditmap = zeros(siz,'logical');
    indices = obj.ind2indInFrame(obj.upstreamRim);
    ditmap(indices) = true;
    indices = obj.ind2indInFrame(obj.downstreamRim);
    ditmap(indices) = true;
    
    if ~obj.checkFrame
        disp('break!');
    end
    if any(size(bitmap) ~= size(citmap)) || any(size(bitmap) ~= size(ditmap))
        disp('break!');
    end
    if any(bitmap(:) ~= ditmap(:))
        result = false;
    else
        result = true;
    end
end
function result = cutPolarizedRim(obj)
% Method to divide lateral rim to segments adjacent to different cells.
% result = cutPolarizedRim(obj)
% Return value is an in-frame bitmap in which the segments were drawn by 1.
    indices = obj.frame(1):obj.frame(2);
    jndices = obj.frame(3):obj.frame(4);
    kndices = obj.map.indices(indices,jndices);
    bitmap = obj.map.cellMap(kndices);
    citmap = obj.map.cellTypeMap(kndices);
    
    mask = bitmap == obj.label;
    bitmap(mask) = 0;
    
    labels = permute(obj.polaritySwitch.var(:,2),[2,3,1]);
    nask = any(citmap == labels,3);
    bitmap(~nask) = 0;
    
    array = unique(bitmap(:))';
    array(array == 0) = [];
    siz = obj.frameSize - 2;
    ditmap = zeros(siz,'logical');
    eitmap = zeros(siz);
    indices = 2:obj.frameSize(1) - 1;
    jndices = 2:obj.frameSize(2) - 1;
    m = [-1,-1; -1,0; -1,1; 0,-1; 0,1; 1,-1; 1,0; 1,1];
    for i = array
        ditmap(:) = false;
        fitmap = bitmap == i;
        for j = 1:size(m,1)
            ditmap = ditmap | fitmap(indices + m(j,1),jndices + m(j,2));
        end
        eitmap = eitmap + double(ditmap);
    end
    ditmap = zeros(siz);
    indices = obj.ind2indInFrame(obj.lateralRim);
    ditmap(indices) = 1;
    
    mask = eitmap > 1;
    ditmap(mask) = 0;
    
    result = ditmap;
end

function result = surfaceElasticity(obj)
% Method to get an energy by the surface elasticity.
% result = surfaceElasticity(obj)
% Return value is the elastic energy of the cell surface.
    n = find(any(obj.PCPSwitch.var == obj.cellType,1),1);
    if n ~= 2 || ~obj.ppolarized
        result = surfaceElasticity@CPMEpithelialCell(obj);
        return;
    end
    
    talb = surfaceElasticity@CPMEpithelialCell(obj);
    
    indices = obj.frame(1):obj.frame(2);
    jndices = obj.frame(3):obj.frame(4);
    kndices = obj.map.indices(indices,jndices);
    bitmap = obj.map.cellMap(kndices);
    bitmap = bitmap == obj.label;
    
    siz = obj.frameSize - 2;
    citmap = zeros(siz,'uint8');
    indices = obj.ind2indInFrame(obj.upstreamRim);
    citmap(indices) = 1;
    
    m = [-1,-1; -1,0; -1,1; 0,-1; 0,1; 1,-1; 1,0; 1,1;];
    indices = 2:size(bitmap,1) - 1;
    jndices = 2:size(bitmap,2) - 1;
    ditmap = zeros(obj.frameSize,'uint8');
    for i = 1:size(m,1)
        kndices = indices + m(i,1);
        lndices = jndices + m(i,2);
        ditmap(kndices,lndices) = ditmap(kndices,lndices) + citmap;
    end
    ditmap(bitmap) = 0;
    u = sum(ditmap(:));
    
    citmap(:) = 0;
    kndices = obj.ind2indInFrame(obj.downstreamRim);
    citmap(kndices) = 1;
    ditmap(:) = 0;
    for i = 1:size(m,1)
        kndices = indices + m(i,1);
        lndices = jndices + m(i,2);
        ditmap(kndices,lndices) = ditmap(kndices,lndices) + citmap;
    end
    ditmap(bitmap) = 0;
    d = sum(ditmap(:));
    
    P = [u;d];
    result = talb + ...
        sum(obj.surfaceYModulus.var(5:6) .* (P - obj.Po(5:6)) .^ 2);
end

function revert(obj)
% Method to reject the update.
% revert(obj)
    revert@CPMEpithelialCell(obj);
    
    n = find(any(obj.PCPSwitch.var == obj.cellType,1),1);
    if n ~= 2 || ~obj.ppolarized
        return;
    end
    
    obj.clearPCP;
    obj.map.PCPMap(obj.upstreamRim) = 1;
    obj.map.PCPMap(obj.downstreamRim) = 2;
end

function result = imageOfSubcellularLocations(obj)
% Method to get an image of subcellular locations.
% result = imageOfSubcellularLocations(obj)
% Return value is an image of subcellular locations of the cell.
    siz = obj.frameSize - 2;
    bitmap = imageOfSubcellularLocations@CPMEpithelialCell(obj);
    
    n = find(any(obj.drawingSwitch.var == obj.cellType,1),1);
    if n == 1 % media
        result = bitmap;
        return;
    end
    
    n = find(any(obj.PCPSwitch.var == obj.cellType,1),1);
    if n == 2 && obj.ppolarized % planar polarized cell.
        color = obj.lut.var(end - 5,:);
        mask = zeros(siz(1),siz(2),'logical');
        indices = obj.ind2indInFrame(obj.upstreamRim);
        mask(indices) = true;
        bitmap(mask(:,:,[1,1,1])) = repmat(color,sum(mask(:)),1);
        color = obj.lut.var(end - 6,:);
        mask(:) = false;
        indices = obj.ind2indInFrame(obj.downstreamRim);
        mask(indices) = true;
        bitmap(mask(:,:,[1,1,1])) = repmat(color,sum(mask(:)),1);
        color = obj.lut.var(end - 7,:);
        mask(:) = false;
        indices = obj.ind2indInFrame(obj.upstreamAJPixels);
        mask(indices) = true;
        bitmap(mask(:,:,[1,1,1])) = repmat(color,sum(mask(:)),1);
        color = obj.lut.var(end - 8,:);
        mask(:) = false;
        indices = obj.ind2indInFrame(obj.downstreamAJPixels);
        mask(indices) = true;
        bitmap(mask(:,:,[1,1,1])) = repmat(color,sum(mask(:)),1);
    end
    
    result = bitmap;
end

end
end
