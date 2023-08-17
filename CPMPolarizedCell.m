
classdef CPMPolarizedCell < CPMCell
properties
    % Shared property.
    polaritySwitch = nan; % SYData.
    
    % Instance variables.
    apicalRim = nan; % double[n]. Indices in ex-frame.
    basalRim = nan; % double[n]. Indices in ex-frame.
    lateralRim = nan; % double[n]. Indices in ex-frame.
    neighborCells = nan; % double[n].
end

methods

function obj = initWithData(obj,data)
% Initialization method with an SYData instance.
% obj = initWithData(obj,data)
    initWithData@CPMCell(obj,data);
    
    s = data.var;
    obj.apicalRim = s.apicalRim;
    obj.basalRim = s.basalRim;
    obj.lateralRim = s.lateralRim;
end
function dest = copy(obj,dest)
% Method to make a shallow copy.
% dest = copy(obj,dest)
% Return value is an instance of the class.
    if nargin < 2
        dest = CPMPolarizedCell;
    end
    copy@CPMCell(obj,dest);

    dest.apicalRim = obj.apicalRim;
    dest.basalRim = obj.basalRim;
    dest.lateralRim = obj.lateralRim;
end

function importFromOwner(obj)
% Method method to import shared variables from the owner.
% importFromOwner(obj)
    importFromOwner@CPMCell(obj);
    
    obj.polaritySwitch = obj.owner.polaritySwitch;
end

function result = data(obj)
% Method to convert the instance to an SYData instance.
% result = data(obj)
% Return value is an SYData instance.
    data_ = data@CPMCell(obj);
    s = data_.var;
    
    s.apicalRim = obj.apicalRim;
    s.basalRim = obj.basalRim;
    s.lateralRim = obj.lateralRim;
    
    data_ = SYData(s);
    result = data_;
end

% function obj = update(obj)
%     update@CPMCell(obj);
%     
%     obj.updateRim;
% end

function markLateral(obj)
% Method to mark lateral sites.
% markLateral(obj)
% Sites adjacet to other polarized cells are marked as lateral.
    indices = obj.frame(1):obj.frame(2);
    jndices = obj.frame(3):obj.frame(4);
    kndices = obj.map.indices(indices,jndices);
    bitmap = obj.map.cellMap(kndices);
    mask = obj.map.cellTypeMap(kndices);
%     labels = obj.map.labelOf('celltype_media');
    labels = permute(obj.polaritySwitch.var(:,2),[2,3,1]);
%     mask = any(mask == labels,3);
    mask = ~any(mask == labels,3);
    bitmap(mask) = 0;
    
    indices = 2:obj.frameSize(1) - 1;
    jndices = 2:obj.frameSize(2) - 1;
    citmap = bitmap(indices,jndices);
    m = [-1,-1; -1,0; -1,1; 0,-1; 0,1; 1,-1; 1,0; 1,1;];
    for i = 1:size(m,1)
        citmap = cat(3,citmap,bitmap(indices + m(i,1),jndices + m(i,2)));
    end
    ditmap = any(citmap ~= obj.label,3) & (citmap(:,:,1) == obj.label);
    eitmap = ditmap & any(citmap == 0,3);
    ditmap = xor(ditmap,eitmap);
    
    indices = find(ditmap);
    obj.lateralRim = obj.ind2indExFrame(indices);
end
function result = apicoBasalRim(obj)
% Method to return apical and basal rim.
% result = apicoBasalRim(obj)
    if isempty(obj.lateralRim)
        result = obj.ind2indExFrame(obj.cellRim);
    else
        indices = obj.ind2indExFrame(obj.cellRim);
        result = setdiff(indices,obj.lateralRim);
    end
end
function markApical(obj,cuticleLabel)
% Method to mark apical sites.
% markApical(obj,cuticleLabel)
% Argument cuticleLabel is an array of aubcellular component labels.
% Sites adjacent to a site with the cuticleLabel are marked as apical.
    indices = obj.frame(1):obj.frame(2);
    jndices = obj.frame(3):obj.frame(4);
    kndices = obj.map.indices(indices,jndices);
    bitmap = obj.map.subcellMap(kndices);
    
    array = obj.ind2indInFrame(obj.apicoBasalRim);
    brray = array;
    siz = obj.frameSize - 2;
    m = [-1,-1; -1,0; -1,1; 0,-1; 0,1; 1,-1; 1,0; 1,1;];
    for i = 1:length(array)
        [r,c] = ind2sub(siz,array(i));
        r = r + 1; c = c + 1;
        for j = 1:size(m,1)
            p = [r,c] + m(j,:);
            if any(bitmap(p(1),p(2)) == cuticleLabel(:))
                array(i) = 0;
                break;
            end
        end
    end
    brray = setdiff(brray,array);
    
    obj.apicalRim = obj.ind2indExFrame(brray);
end
function readApicalRim(obj,maskData)
% Method to mark apical rim.
% readApicalRim(obj,maskData)
% Argument maskData is an SYData instance of bitmap marked apical.
    indices = obj.frame(1) + 1:obj.frame(2) - 1;
    jndices = obj.frame(3) + 1:obj.frame(4) - 1;
    kndices = obj.map.indices(indices,jndices);
    
    bitmap = obj.map.cellMap(kndices);
    citmap = maskData.var(kndices);
    bitmap = citmap & (bitmap == obj.label);
    
    lndices = find(bitmap);
    obj.apicalRim = obj.ind2indExFrame(lndices);
    obj.basalRim = setdiff(obj.apicoBasalRim,obj.apicalRim);
end
function updateRim(obj)
% Method to update rim marking.
% updateRim(obj)
    n = find(any(obj.polaritySwitch.var == obj.cellType,1),1);
    if n == 2 % polarized cell.
        obj.markLateral;

        siz = obj.frameSize - 2;
        bitmap = zeros(siz);
        indices = obj.ind2indInFrame(obj.apicoBasalRim);
        bitmap(indices) = 1;
        bitmap = IPConnectedComponents.connectedBinaryComponents(bitmap,8);
        indices = find(bitmap);
        indices = obj.ind2indExFrame(indices);

        for i = 1:max(bitmap(:))
            jndices = find(bitmap == i);
            jndices = obj.ind2indExFrame(jndices);
            if isempty(intersect(obj.apicalRim,jndices))
                indices = setdiff(indices,jndices);
            end
        end

        obj.apicalRim = indices;
        obj.basalRim = setdiff(obj.apicoBasalRim,obj.apicalRim);
    end
end
function result = isPixelApical(obj,pixel)
% Method to check a site if apical.
% result = isPixelApical(obj,pixel)
% Argument pixel is a position of the site: int[2]
% Return value is a boolean.
    if length(pixel) > 1
        pixel = sub2ind(obj.map.frameSize,pixel(1),pixel(2));
    end
    
    if any(obj.apicalRim == pixel)
        result = true;
    else
        result = false;
    end
end
function result = isPixelLateral(obj,pixel)
% Method to check a site if lateral.
% result = isPixelLateral(obj,pixel)
% Argument pixel is a position of the site: int[2]
% Return value is a boolean.
    if length(pixel) > 1
        pixel = sub2ind(obj.map.frameSize,pixel(1),pixel(2));
    end
    
    if any(obj.lateralRim == pixel)
        result = true;
    else
        result = false;
    end
end
function result = isPixelBasal(obj,pixel)
% Method to check a site if basal.
% result = isPixelBasal(obj,pixel)
% Argument pixel is a position of the site: int[2]
% Return value is a boolean.
    if length(pixel) > 1
        pixel = sub2ind(obj.map.frameSize,pixel(1),pixel(2));
    end
    
    if any(obj.basalRim == pixel)
        result = true;
    else
        result = false;
    end
end

function makePoCurrentLength(obj)
% Method to set apical, lateral, and basal Po current length.
% makePoCurrentLength(obj)
    makePoCurrentLength@CPMCell(obj);
    
    n = find(any(obj.polaritySwitch.var == obj.cellType,1),1);
    if n == 1
        obj.Po = cat(1,obj.Po,0,0,0);
        return;
    end
    
    indices = obj.frame(1):obj.frame(2);
    jndices = obj.frame(3):obj.frame(4);
    kndices = obj.map.indices(indices,jndices);
    bitmap = obj.map.cellMap(kndices);
    bitmap = bitmap == obj.label;
    
    siz = obj.frameSize - 2;
    citmap = zeros(siz,'uint8');
    indices = obj.ind2indInFrame(obj.apicalRim);
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
    a = sum(ditmap(:));
    
    citmap(:) = 0;
    kndices = obj.ind2indInFrame(obj.lateralRim);
    citmap(kndices) = 1;
    ditmap(:) = 0;
    for i = 1:size(m,1)
        kndices = indices + m(i,1);
        lndices = jndices + m(i,2);
        ditmap(kndices,lndices) = ditmap(kndices,lndices) + citmap;
    end
    ditmap(bitmap) = 0;
    l = sum(ditmap(:));
    
    citmap(:) = 0;
    kndices = obj.ind2indInFrame(obj.basalRim);
    citmap(kndices) = 1;
    ditmap(:) = 0;
    for i = 1:size(m,1)
        kndices = indices + m(i,1);
        lndices = jndices + m(i,2);
        ditmap(kndices,lndices) = ditmap(kndices,lndices) + citmap;
    end
    ditmap(bitmap) = 0;
    b = sum(ditmap(:));
    
    obj.Po = cat(1,obj.Po,a,l,b);
end

function result = surfaceElasticity(obj)
% Method to get an energy by the surface elasticity.
% result = surfaceElasticity(obj)
% Return value is the elastic energy of the cell surface.
    n = find(any(obj.polaritySwitch.var == obj.cellType,1),1);
    if n ~= 2
        result = surfaceElasticity@CPMCell(obj);
        return;
    end
    
    t = surfaceElasticity@CPMCell(obj);
    
    indices = obj.frame(1):obj.frame(2);
    jndices = obj.frame(3):obj.frame(4);
    kndices = obj.map.indices(indices,jndices);
    bitmap = obj.map.cellMap(kndices);
    bitmap = bitmap == obj.label;
    
    siz = obj.frameSize - 2;
    citmap = zeros(siz,'uint8');
    indices = obj.ind2indInFrame(obj.apicalRim);
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
    a = sum(ditmap(:));
    
    citmap(:) = 0;
    kndices = obj.ind2indInFrame(obj.lateralRim);
    citmap(kndices) = 1;
    ditmap(:) = 0;
    for i = 1:size(m,1)
        kndices = indices + m(i,1);
        lndices = jndices + m(i,2);
        ditmap(kndices,lndices) = ditmap(kndices,lndices) + citmap;
    end
    ditmap(bitmap) = 0;
    l = sum(ditmap(:));
    
    citmap(:) = 0;
    kndices = obj.ind2indInFrame(obj.basalRim);
    citmap(kndices) = 1;
    ditmap(:) = 0;
    for i = 1:size(m,1)
        kndices = indices + m(i,1);
        lndices = jndices + m(i,2);
        ditmap(kndices,lndices) = ditmap(kndices,lndices) + citmap;
    end
    ditmap(bitmap) = 0;
    b = sum(ditmap(:));
    
    P = [a; l; b];
    result = t + sum(obj.surfaceYModulus.var(2:4) .* (P - obj.Po(2:4)) .^ 2);
end

function result = imageOfSubcellularLocations(obj)
% Method to get an image of subcellular locations.
% result = imageOfSubcellularLocations(obj)
% Return value is an image of subcellular locations of the cell.
    siz = obj.frameSize - 2;
    bitmap = imageOfSubcellularLocations@CPMCell(obj);
    
    n = find(any(obj.drawingSwitch.var == obj.cellType,1),1);
    if n == 1 % media
        result = bitmap;
        return;
    end
    
    n = find(any(obj.polaritySwitch.var == obj.cellType,1),1);
    if n == 2 % polarized cell.
        %  lateral rim.
        color = obj.lut.var(end - 1,:);
        mask = zeros(siz(1),siz(2),'logical');
        indices = obj.ind2indInFrame(obj.lateralRim);
        mask(indices) = true;
        bitmap(mask(:,:,[1,1,1])) = repmat(color,sum(mask(:)),1);
        %  basal rim.
        color = obj.lut.var(end - 2,:);
        mask(:) = false;
        indices = obj.ind2indInFrame(obj.basalRim);
        mask(indices) = true;
        bitmap(mask(:,:,[1,1,1])) = repmat(color,sum(mask(:)),1);
        %  apical rim.
        color = obj.lut.var(end - 3,:);
        mask(:) = false;
        indices = obj.ind2indInFrame(obj.apicalRim);
        mask(indices) = true;
        bitmap(mask(:,:,[1,1,1])) = repmat(color,sum(mask(:)),1);
    end
    
    result = bitmap;
end

end
end
