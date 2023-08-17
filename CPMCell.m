% CPM Class CPMCell < SYObject.
% Written by Satoshi Yamashita.
% Data class representing a cell.

classdef CPMCell < SYObject
properties
    % Shared property.
    owner = nan; % CPMSimulator.
    map = nan; % CPMMap.
    subcellIndices = nan; % SYData.
    subcellConverter = nan; % SYData.
    rimCounts = nan; % SYData.
    cellEnergy = nan; % SYData.
    imageData = nan; % SYData.
    drawingSwitch = nan; % SYData.
    lut = nan; % SYData.
    % % Parameters.
    contactEnergy = nan; % SYData.
    surfaceYModulus = nan; % SYData.
    areaBModulus = nan; % SYData.
    HamiltonianSwitch = nan; % SYData.
    
    % Instance variables.
    cellIndex = nan; % double: index of cell in simulator.cellArray.
    label = nan; % double : label of cell on map.cellMap.
    cellType = nan; % double.
    frame = nan; % double(4).
    frameSize = nan; % double(2);
    cellRim = nan; % double[n]. Indices in in-frame.
    subcellRim = nan; % double[n]. subcellRim \supset cellRim
    Ao = nan; % double.
    Po = nan; % double[n].
    
    oldData = nan; % SYData.
    newEnergy = nan; % double.
end

methods
function obj = initWithData(obj,data)
% Initialization method with an SYData instance.
% obj = initWithData(obj,data)
    if isnan(obj.owner)
        disp('CPMCell instance needs owner before initialization.');
        return;
    end
    
    s = data.var;
    obj.cellIndex = s.cellIndex;
    obj.label = s.label;
    obj.cellType = s.cellType;
    obj.frame = s.frame;
    obj.frameSize = s.frameSize;
    obj.cellRim = s.cellRim;
    obj.subcellRim = s.subcellRim;
    obj.Ao = s.Ao;
    obj.Po = s.Po;
    obj.surfaceYModulus = ...
        obj.owner.surfaceYModulusArray.objectAtIndex(obj.cellType);
    obj.areaBModulus = ...
        obj.owner.areaBModulusArray.objectAtIndex(obj.cellType);
end

function dest = copy(obj,dest)
% Method to make a shallow copy.
% dest = copy(obj,dest)
% Return value is an instance of the class.
    if nargin < 2
        dest = CPMCell;
    end
    copy@SYObject(obj,dest);
    
    dest.owner = obj.owner;
    dest.initWithData(obj.data);
end

% function delete(obj)
%     delete@SYObject(obj);
% end

function set.owner(obj,newOwner)
    obj.owner = newOwner;
    obj.setOwner(newOwner);
end
function setOwner(obj,~)
% Bypass method from set.owner.
% Do not call directly.
    obj.importFromOwner;
end
function importFromOwner(obj)
% Method method to import shared variables from the owner.
% importFromOwner(obj)
    obj.map = obj.owner.map;
    obj.subcellIndices = obj.owner.subcellIndices;
    obj.subcellConverter = obj.owner.subcellConverter;
    obj.rimCounts = obj.owner.rimCounts;
    obj.cellEnergy = obj.owner.cellEnergy;
    obj.imageData = obj.owner.imageData;
    obj.drawingSwitch = obj.owner.drawingSwitch;
    obj.lut = obj.owner.lut;
    obj.contactEnergy = obj.owner.contactEnergy;
    obj.HamiltonianSwitch = obj.owner.HamiltonianSwitch;
end

function result = data(obj)
% Method to convert the instance to an SYData instance.
% result = data(obj)
% Return value is an SYData instance.
    s.cellIndex = obj.cellIndex;
    s.label = obj.label;
    s.cellType = obj.cellType;
    s.frame = obj.frame;
    s.frameSize = obj.frameSize;
    s.cellRim = obj.cellRim;
    s.subcellRim = obj.subcellRim;
    s.Ao = obj.Ao;
    s.Po = obj.Po;
    
    data = SYData(s);
    result = data;
end

% Depleted because of frequent call from initWithData().
% function set.cellType(obj,newCellType)
%     obj.cellType = newCellType;
%     obj.setCellType(newCellType);
% end
% function setCellType(~,~)
%     % (^_^).
% end

function set.subcellRim(obj,newRim)
    obj.subcellRim = newRim;
    obj.setSubcellRim(newRim);
end
function setSubcellRim(obj,newRim)
% Bypass method from set.subcellRim.
% Do not call directly.
    obj.rimCounts.var(obj.cellIndex) = length(newRim);
end

function enframe(obj)
% Method to find an enclosing frame.
% enframe(obj)
    y = any(obj.map.cellMap == obj.label,2);
    x = any(obj.map.cellMap == obj.label,1);
    f = [find(y,1),find(y,1,'last'),find(x,1),find(x,1,'last')];
    if f(1) == 1 && f(2) == obj.map.frameSize(1) && any(~y)
        f(1) = find(~y,1,'last') + 1;
        f(2) = find(~y,1) - 1 + obj.map.frameSize(1);
    end
    if f(3) == 1 && f(4) == obj.map.frameSize(2) && any(~x)
        f(3) = find(~x,1,'last') + 1;
        f(4) = find(~x,1) - 1 + obj.map.frameSize(2);
    end
    obj.frame = f + [-1,1,-1,1];
    obj.frameSize = [obj.frame(2) - obj.frame(1) + 1, ...
        obj.frame(4) - obj.frame(3) + 1];
    
%     if ~obj.checkFrame
%         disp('break!');
%     end
end
function result = checkFrame(obj)
% Method to check frame validity.
% result = checkFrame(obj)
% Return value is a boolean.
    indices = obj.frame(1) + 1:obj.frame(2) - 1;
    jndices = obj.frame(3) + 1:obj.frame(4) - 1;
    kndices = obj.map.indices(indices,jndices);
    siz = obj.frameSize - 2;
    if any(size(kndices) ~= siz)
        result = false;
    else
        result = true;
    end
end
function readCellType(obj)
% Method to import cell type from map.
% readCellType(obj)
    indices = obj.frame(1):obj.frame(2);
    jndices = obj.frame(3):obj.frame(4);
    kndices = obj.map.indices(indices,jndices);
    bitmap = obj.map.cellMap(kndices);
    mask = bitmap == obj.label;
    bitmap = obj.map.cellTypeMap(kndices);
    array = unique(bitmap(mask));
    obj.cellType = array(1);
end
function markSurface(obj)
% Method to mark cell surface to be updated.
% markSurface(obj)
    indices = obj.frame(1):obj.frame(2);
    jndices = obj.frame(3):obj.frame(4);
    kndices = obj.map.indices(indices,jndices);
    bitmap = obj.map.cellMap(kndices);
    mask = bitmap == obj.label;
    
    citmap = zeros(length(indices) - 2,length(jndices) - 2,'logical');
    ditmap = mask(2:end - 1,2:end - 1);
    m = [-1,-1; -1,0; -1,1; 0,-1; 0,1; 1,-1; 1,0; 1,1;];
    for i = 1:size(m,1)
        eitmap = ditmap ~= mask((2:end - 1) + m(i,1), ...
            (2:end - 1) + m(i,2));
        citmap = citmap | eitmap;
    end
    eitmap = citmap & mask(2:end - 1,2:end - 1);
    obj.cellRim = find(eitmap);
end
function makePoCurrentLength(obj)
% Method to initialize reference lengths of cell surface.
% makePoCurrentLength(obj)
    indices = obj.frame(1):obj.frame(2);
    jndices = obj.frame(3):obj.frame(4);
    kndices = obj.map.indices(indices,jndices);
    bitmap = obj.map.cellMap(kndices);
    bitmap = bitmap == obj.label;
    
    m = [-1,-1; -1,0; -1,1; 0,-1; 0,1; 1,-1; 1,0; 1,1;];
    indices = 2:size(bitmap,1) - 1;
    jndices = 2:size(bitmap,2) - 1;
    citmap = zeros(size(bitmap,1),size(bitmap,2),'uint8');
    for j = 1:size(m,1)
        kndices = indices + m(j,1);
        lndices = jndices + m(j,2);
        citmap(kndices,lndices) = ...
            citmap(kndices,lndices) + uint8(bitmap(indices,jndices));
    end
    citmap(bitmap) = 0;
    obj.Po = sum(citmap(:));
    
    obj.surfaceYModulus = ...
        obj.owner.surfaceYModulusArray.objectAtIndex(obj.cellType);
end
function collectRim(obj)
% Method to mark subcellular components surface to be updated.
% collectRim(obj)
    indices = obj.frame(1):obj.frame(2);
    jndices = obj.frame(3):obj.frame(4);
    kndices = obj.map.indices(indices,jndices);
    bitmap = obj.map.cellMap(kndices);
    mask = bitmap == obj.label;
    bitmap = obj.map.subcellMap(kndices);
    bitmap(~mask) = 0;
    citmap = bitmap(2:end - 1,2:end - 1);
    ditmap = zeros(size(citmap,1),size(citmap,2),'logical');
    m = [-1,0; 0,-1; 0,1; 1,0];
    for i = 1:size(m,1)
        eitmap = citmap ~= bitmap((2:end - 1) + m(i,1), ...
            (2:end - 1) + m(i,2));
        ditmap = ditmap | eitmap;
    end
    ditmap = ditmap & mask(2:end - 1,2:end - 1);
    obj.subcellRim = find(ditmap);
end
function makeAoCurrentVolume(obj)
% Method to initialize reference volume of the cell and subcellular
%   components in it.
% makeAoCurrentVolume(obj)
    obj.Ao = zeros(length(obj.subcellIndices.var) + 1,1);
    
    indices = obj.frame(1) + 1:obj.frame(2) - 1;
    jndices = obj.frame(3) + 1:obj.frame(4) - 1;
    kndices = obj.map.indices(indices,jndices);
    bitmap = obj.map.cellMap(kndices);
    mask = bitmap == obj.label;
    
    obj.Ao(1) = sum(mask(:));
    
    bitmap = obj.map.subcellMap(kndices);
    bitmap(~mask) = 0;
    array = unique(bitmap(:))';
    array(array < 1) = [];
    A = zeros(length(obj.subcellIndices.var),1);
    for i = array
        a = obj.subcellIndices.var == i;
        b = bitmap == i;
        b = sum(b(:));
        A(a) = b;
    end
    obj.Ao(2:end) = A;
    
    obj.areaBModulus = ...
        obj.owner.areaBModulusArray.objectAtIndex(obj.cellType);
end
function registerEnergy(obj)
% Method to register the energy of the cell.
% registerEnergy(obj)
    H = obj.getHamiltonian;
    obj.cellEnergy.var(obj.cellIndex) = H;
end

function result = ind2indInFrame(obj,indices)
% Method to convert indices in map to indices in cell frame.
% result = ind2indInFrame(obj,indices)
% Argument indices is an array of indices.
% Return value is an array of indices.
% Note that the cell frame is frameSize - 2.
    if isempty(indices)
        result = [];
        return;
    end
    
    siz = obj.frameSize - 2;
    sub = obj.ind2subInFrame(indices);
    result = sub2ind(siz,sub(:,1),sub(:,2));
end
function result = ind2subInFrame(obj,indices)
% Method to convert indices in map to subscripts in cell frame.
% result = ind2subInFrame(obj,indices)
% Argument indices is an array of indices.
% Return value is a {n,2} array of subscripts.
% Note that the cell frame is frameSize - 2.
    if isempty(indices)
        result = [];
        return;
    end
    
    [r,c] = ind2sub(obj.map.frameSize,indices);
    r = mod(r - obj.frame(1) - 1,obj.map.frameSize(1)) + 1;
    c = mod(c - obj.frame(3) - 1,obj.map.frameSize(2)) + 1;
    result = [r,c];
end
function result = ind2indExFrame(obj,indices)
% Method to convert indices in cell frame to indices in map.
% result = ind2indExFrame(obj,indices)
% Argument indicies is an array of indices.
% Note that the cell frame is frameSize - 2.
    if isempty(indices)
        result = [];
        return;
    end
    
    sub = obj.ind2subExFrame(indices);
    
    result = obj.map.indicesFromSub(sub);
end
function result = ind2subExFrame(obj,indices)
% Method to convert indices in cell frame to subscripts in map.
% result = ind2subExFrame(obj,indices)
% Argument indices is an array of indices.
% Note that the cell frame is frameSize - 2.
    if isempty(indices)
        result = [];
        return;
    end
    
    siz = obj.frameSize - 2;
    [r,c] = ind2sub(siz,indices);
    r = r + obj.frame(1);
    c = c + obj.frame(3);
    result = [r,c];
end

function result = getHamiltonian(obj)
% Method to get energy of the cell.
% result = getHamiltonian(obj)
% Return value is the energy of the cell.
    n = find(any(obj.HamiltonianSwitch.var == obj.cellType,1),1);
    if n == 1 % media.
        result = 0;
    else % cell.
        H = 0;
        % Surface tension.
        H = H + obj.surfaceContractility;
        H = H + obj.surfaceElasticity;
        % Area constraint.
        H = H + obj.areaElasticity;

        result = H;
    end
end
function result = surfaceContractility(obj)
% Method to get a surface contractilty of subcellular components in the
%   cell.
% result = surfaceContractility(obj)
% Return value is a total contact energy of the all components in the cell.
    indices = obj.frame(1):obj.frame(2);
    jndices = obj.frame(3):obj.frame(4);
    kndices = obj.map.indices(indices,jndices);
    bitmap = obj.map.cellMap(kndices);
    mask = bitmap == obj.label;
    bitmap = double(obj.map.subcellMap(kndices));
    citmap = bitmap .* mask;
    
    tau = 0;
    
    % Subcellular components surface contractility.
    m = [-1,-1; -1,0; -1,1; 0,-1; 0,1; 1,-1; 1,0; 1,1;];
    indices = 2:size(bitmap,1) - 1;
    jndices = 2:size(bitmap,2) - 1;
    array = unique(citmap(:))';
    array(array < 1) = [];
    stack = zeros(size(citmap,1),size(citmap,2),length(array));
    for i = 1:length(array)
        ditmap = citmap(indices,jndices) == array(i);
        for j = 1:size(m,1)
            kndices = indices + m(j,1);
            lndices = jndices + m(j,2);
            eitmap = citmap(kndices,lndices) ~= array(i);
            stack(kndices,lndices,i) = ...
                stack(kndices,lndices,i) + (ditmap & eitmap);
        end
    end
    brray = unique(bitmap(:))';
    brray(brray < 1) = [];
    for i = 1:length(array)
        r = obj.subcellIndices.var == array(i);
        for j = brray
            c = obj.subcellIndices.var == j;
            ditmap = bitmap == j;
            eitmap = stack(:,:,i) .* ditmap .* mask;
            tau = tau + sum(eitmap(:)) * obj.contactEnergy.var(r,c,1);
            eitmap = stack(:,:,i) .* ditmap .* ~mask;
            tau = tau + sum(eitmap(:)) * obj.contactEnergy.var(r,c,2);
        end
    end
    
    result = tau;
end
function result = surfaceElasticity(obj)
% Method to get an energy by the surface elasticity.
% result = surfaceElasticity(obj)
% Return value is the elastic energy of the cell surface.
    indices = obj.frame(1):obj.frame(2);
    jndices = obj.frame(3):obj.frame(4);
    kndices = obj.map.indices(indices,jndices);
    bitmap = obj.map.cellMap(kndices);
    bitmap = bitmap == obj.label;
    
    m = [-1,-1; -1,0; -1,1; 0,-1; 0,1; 1,-1; 1,0; 1,1;];
    indices = 2:size(bitmap,1) - 1;
    jndices = 2:size(bitmap,2) - 1;
    citmap = zeros(size(bitmap,1),size(bitmap,2),'uint8');
    for j = 1:size(m,1)
        kndices = indices + m(j,1);
        lndices = jndices + m(j,2);
        citmap(kndices,lndices) = ...
            citmap(kndices,lndices) + uint8(bitmap(indices,jndices));
    end
    citmap(bitmap) = 0;
    P = sum(citmap(:));
    
    result = obj.surfaceYModulus.var(1) * (P - obj.Po(1)) ^ 2;
end
function result = areaElasticity(obj)
% Method to get an energy by the volume elasticity.
% result = areaElasticity(obj)
% Return value is the elastic energy of the cell and subcellular
%   components.
    H = 0;
    indices = obj.frame(1) + 1:obj.frame(2) - 1;
    jndices = obj.frame(3) + 1:obj.frame(4) - 1;
    kndices = obj.map.indices(indices,jndices);
    bitmap = obj.map.cellMap(kndices);
    mask = bitmap == obj.label;
    
    A = sum(mask(:));
    H = H + obj.areaBModulus.var(1) * (A - obj.Ao(1)) ^ 2;
    
    bitmap = obj.map.subcellMap(kndices);
    bitmap(~mask) = 0;
    array = unique(bitmap(:))';
    array(array < 1) = [];
    for i = array
        j = find(obj.subcellIndices.var == i,1) + 1;
        citmap = bitmap == i;
        A = sum(citmap(:));
        H = H + obj.areaBModulus.var(j) * (A - obj.Ao(j)) ^ 2;
    end
    
    result = H;
end

function result = tryUpdate(obj)
% Method to raies a candidate of change.
% result = tryUpdate(obj)
% Return value is a boolean indicating if the candidate was raised.
    
    % select a pixel.
    i = ceil(rand() * obj.rimCounts.var(obj.cellIndex));
    if i == 0
        i = 1;
    end
    p = obj.ind2subExFrame(obj.subcellRim(i));
    
    indices = obj.map.indices(p(1) - 1:p(1) + 1,p(2) - 1:p(2) + 1);
    bitmap = obj.map.cellMap(indices);
    citmap = bitmap ~= bitmap(2,2);
    ditmap = obj.map.subcellMap(indices);
    eitmap = obj.map.cellTypeMap(indices);
    
    % check connectedness.
    array = citmap([1,2,3,6,9,8,7,4]);
    brray = citmap([2,3,6,9,8,7,4,1]);
    if sum(array ~= brray) > 2
        connectedFlag = false;
    else
        connectedFlag = true;
    end
    
    % select a neighbor.
    index = randi(4) * 2;
    r = obj.subcellIndices.var == ditmap(2,2);
    c = obj.subcellIndices.var == ditmap(index);
    if ~citmap(index) % the neighbor is in the same cell.
        l = obj.subcellConverter.var(r,c,1);
    else % the neighbor is in a different cell.
        if ~connectedFlag
            result = false;
            return;
        end
        
        l = obj.subcellConverter.var(r,c,2);
    end
    if isnan(l)
        result = false;
        return;
    end
    
    l = [bitmap(index),l,eitmap(index)];
    neighbors = unique(bitmap(:));
    neighbors(neighbors < 1) = [];
    obj.owner.updatePixel(p,l,neighbors);
    
    result = true;
end
function update(obj)
% Method to update the cell tentatively.
% update(obj)
    obj.oldData = obj.data;
    
    obj.enframe;
    obj.markSurface;
    obj.collectRim;
end
function result = dH(obj)
% Method to get a change in the energy.
% result = dH(obj)
% Return value is a change in the cell energy before and after the update.
    obj.newEnergy = obj.getHamiltonian;
    result = obj.newEnergy - obj.cellEnergy.var(obj.cellIndex);
end
function fixUpdate(obj)
% Method to accept the update.
% fixUpdate(obj)
    obj.cellEnergy.var(obj.cellIndex) = obj.newEnergy;
end
function revert(obj)
% Method to reject the update.
% revert(obj)
    obj.initWithData(obj.oldData);
end

function drawCell(obj)
% Method to draw an image of the cell on the simulator image viewer.
% drawCell(obj)
    bitmap = uint8(obj.imageOfCell);
    
    indices = obj.frame(1) + 1:obj.frame(2) - 1;
    jndices = obj.frame(3) + 1:obj.frame(4) - 1;
    kndices = obj.map.indices(indices,jndices);
    mask = obj.map.cellMap(kndices) == obj.label;
    bitmap(~mask(:,:,[1,1,1])) = 0;
    
    l = obj.map.frameSize(1) * obj.map.frameSize(2);
    lndices = cat(3,kndices,kndices + l,kndices + l * 2);
    citmap = obj.imageData.var(lndices);
    citmap(mask(:,:,[1,1,1])) = 0;
    
    citmap = citmap + bitmap;
    obj.imageData.var(lndices) = citmap;
end
function result = imageOfCell(obj)
% Method to get an image of the cell.
% result = imageOfCell(obj)
% Return value is an image of the cell in the cell frame. Note that the
% frame is frameSize - 2.
    bitmap = obj.imageOfSubcellularLocations;
    citmap = obj.imageOfSubcellularComponents;
    bitmap = cat(4,bitmap,citmap);
    
    bitmap = min(bitmap,[],4);
    result = bitmap;
end
function result = imageOfSubcellularLocations(obj)
% Method to get an image of subcellular locations.
% result = imageOfSubcellularLocations(obj)
% Return value is an image of subcellular locations of the cell.
    siz = obj.frameSize - 2;
    bitmap = zeros(siz(1),siz(2),3);
    bitmap(:) = 255;
    
    mask = zeros(siz(1),siz(2),'logical');
    mask(obj.cellRim) = true;
    n = find(any(obj.drawingSwitch.var == obj.cellType,1),1);
    if n == 1 % media
        color = obj.lut.var(end,:);
    else % cell.
        color = obj.lut.var(end - 1,:);
    end
    bitmap(mask(:,:,[1,1,1])) = repmat(color,sum(mask(:)),1);
    
    result = bitmap;
end
function result = imageOfSubcellularComponents(obj)
% Method to get an image of subcellular components.
% result = imageOfSubcellularComponents(obj)
% Return value is an image of subcellular components in the cell.
    siz = obj.frameSize - 2;
    bitmap = zeros(siz(1),siz(2),3);
    bitmap(:) = 255;
    
    indices = obj.frame(1) + 1:obj.frame(2) - 1;
    jndices = obj.frame(3) + 1:obj.frame(4) - 1;
    kndices = obj.map.indices(indices,jndices);
    mask = obj.map.cellMap(kndices) == obj.label;
    citmap = double(obj.map.subcellMap(kndices)) .* double(mask);
    
    array = unique(citmap(:))';
    array(array < 1) = [];
    for i = array
        color = obj.lut.var(i,:);
        mask = citmap == i;
        bitmap(mask(:,:,[1,1,1])) = repmat(color,sum(mask(:)),1);
    end
    
    result = bitmap;
end

end
end
