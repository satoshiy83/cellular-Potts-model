% CPM Class CPMSimulator < SYObject.
% Written by Satoshi Yamashita.
% Controller class managing a simulation.

classdef CPMSimulator < SYObject
properties
    % Map and cell classes specifier.
    mapClass = nan; % string.
    cellClass = nan; % string.
    
    % Data representations.
    map = nan; % CPMMap instance.
    hint = nan; % SYDictionary.
    
    % Simulation data managers.
    cellLabels = nan; % double[n,1].
    cellArray = nan; % SYArray.
    rimCounts = nan; % SYData.
    updateList = nan; % double[1,n].
    cellEnergy = nan; % SYData.
    subcellIndices = nan; % SYData.
    subcellConverter = nan; % SYData.
    intervalCounter = nan; % double.
    frameCounter = nan; % double.
    stopFlag = nan; % bool.
    endFlag = nan; % bool.
    dH = nan; % double.
    p = nan; % double.
    
    % Parameters.
    contactEnergy = nan; % SYData.
    surfaceYModulusArray = nan; % SYArray.
    areaBModulusArray = nan; % SYArray.
    HamiltonianSwitch = nan; % SYData.
    T = nan; % double.
    
    % Display.
    imageData = nan; % SYData.
    drawingSwitch = nan; % SYData.
    lut = nan; % SYData.
    % Interface.
    window = nan; % figure panel handle.
    textField_progress = nan; % text field on the window.
    imageView = nan; % image handle on the window.
    
    % Output.
    frameInterval = nan; % double.
    stackLimit = nan; % double.
    outputName = nan; % string.
end

methods
function obj = CPMSimulator
% Controller class managing a simulation.
% obj = CPMSimulator
% Subclass of CPMSimulator may define map class and cell class here.
    obj.mapClass = 'CPMMap';
    obj.cellClass = 'CPMCell';
end

function obj = initWithMap(obj,map,hint)
% Initialization method with CPMMap instance and SYDictionary hint.
% obj = initWithMap(obj,map,hint)
% Argument map is a CPMMap family instance.
% Argument hint is an SYDictionary instance holding parameters.
% Parameters:
%   CPMMapClassKey: name of map class. Default value is set at allocation.
%   CPMCellClassKey: name of cell class. Default value is set at allcation.
%   CPMSubcellIndicesKey: (SYData) array of subcellular components labels.
%   CPMSubcellConverterKey: (SYData) {n,n,2} matrix defining conversion of
%       subcellular components label.
%   CPMContactEnergyKey: (SYData) {n,n,2} matrix defining contact  energy 
%       between subcellular components.
%   CPMSurfaceYModulusArrayKey: (SYArray) array containing surface 
%       elasticity modulus of subcellular locations for cell types.
%   CPMAreaBModulusArrayKey: (SYArray) array containing area  elasticity
%       modulus of cell and subcellular components for cell types.
%   CPMHamiltonianSwitchKey: (YSData) {n,2} matrix  specifying cell types
%       either having energy or not.
%   CPMTKey: (double) temperature.
%   CPMDrawingSwitchKey:  (SYData) {n,2}  matrix specifying cell types
%       either  media or cell.
%   CPMLutKey: (SYData) {n,3} look up table.
%   CPMOutputNameKey: (string) path to save results.
%   CPMFrameIntervalKey: (double) number of updates between frames.
%   CPMStackLimitKey: (double) number of slices in a stack to be saved.

    obj.init;
    
    if isnumeric(map)
        map = CPMMap(map(:,:,1),map(:,:,2),map(:,:,3),hint);
    end
    obj.map = map;
    obj.hint = hint;
    
    obj.readHint;
    
    indices = unique(map.cellMap(:));
    indices(indices < 1) = [];
    obj.cellLabels = indices;
    c = length(indices);
    
    obj.cellArray = SYArray;
    obj.rimCounts = SYData(zeros(c,1));
    obj.updateList = [];
    obj.cellEnergy = SYData(zeros(c,1));
    range = [1,size(obj.lut.var,1)];
    bitmap = map.subcellMap;
    context = SYGraphicsContext(SYData(bitmap), ...
        size(bitmap,2),size(bitmap,1),8, ...
        SYGraphicsContext.CompositeModeOver, ...
        SYGraphicsContext.ColorSpaceIndexed,obj.lut.var,range);
    image = SYImage(context);
    bitmap = image.drawBitmapRep(nan);
    
    obj.imageData = SYData(bitmap);
    
    obj.prepareCells(indices);
end
function obj = initWithData(obj,data)
% Initialization method with SYData instance.
% obj = initWithData(obj,data)
    hintData = SYData(data.var.hintData);
    hint_ = SYDictionary;
    hint_.initWithData(hintData);
    obj.hint = hint_;
    
    mapData = SYData(data.var.mapData);
    fh = str2func(obj.mapClass);
    map_ = fh();
    map_.initWithData(mapData);
    obj.map = map_;
    
    obj.readHint;
    
    indices = unique(map_.cellMap(:));
    indices(indices < 1) = [];
    obj.cellLabels = indices;
    c = length(indices);
    
    obj.cellArray = SYArray;
    obj.rimCounts = SYData(zeros(c,1));
    obj.updateList = [];
    obj.cellEnergy = SYData(zeros(c,1));
    range = [1,size(obj.lut.var,1)];
    bitmap = map_.subcellMap;
    context = SYGraphicsContext(SYData(bitmap), ...
        size(bitmap,2),size(bitmap,1),8, ...
        SYGraphicsContext.CompositeModeOver, ...
        SYGraphicsContext.ColorSpaceIndexed,obj.lut.var,range);
    image = SYImage(context);
    bitmap = image.drawBitmapRep(nan);
    
    obj.imageData = SYData(bitmap);
    
    cellsData = SYData(data.var.cellsData);
    array = SYArray;
    array.initWithData(cellsData);
    fh = str2func(obj.cellClass);
    for i = 1:array.count
        c = fh();
        c.owner = obj;
        c.initWithData(array.objectAtIndex(i));
        obj.cellArray.addObject(c);
    end
    obj.makeCellsRegisterEnergy;
    obj.cellArray.makeObjectsPerformSelector(@drawCell);
end

function dest = copy(obj,dest)
% Method to make a shallow copy.
% dest = copy(obj,dest)
% Return value is an instance of the class.
    if nargin < 2
        dest = CPMSimulator;
    end
    copy@SYObject(obj,dest);
    
    dest.initWithMap(obj.map.copy,obj.hint.copy);
    
    dest.cellLabels = obj.cellLabels;
    dest.cellArray = obj.cellArray.copy;
    dest.rimCounts = obj.rimCounts.copy;
    dest.updateList = obj.updateList;
    dest.cellEnergy = obj.cellEnergy.copy;
    
    dest.imageData = obj.imageData.copy;
    
    if ishandle(obj.window)
        dest.window = figure;
    end
end

function delete(obj)
% Method to clean the instance.
% delete(obj)
    if ishandle(obj.window)
        delete(obj.window);
    end
    
    delete@SYObject(obj);
end

function result = data(obj)
% Method to convert the instance to an SYData instance.
% result = data(obj)
% Return value is an SYData instance.
    s.hintData = obj.hint.data.var;
    s.mapData = obj.map.data.var;
    
    array = SYArray;
    for i = 1:obj.cellArray.count
        c = obj.cellArray.objectAtIndex(i);
        array.addObject(c.data);
    end
    s.cellsData = array.data.var;
    
    result = SYData(s);
end

function readHint(obj)
% Method to import parameters from the hint.
% readHint(obj)
    if ~obj.hint.isNanForKey("CPMMapClassKey")
        obj.mapClass = obj.hint.objectForKey("CPMMapClassKey");
    end
    if ~obj.hint.isNanForKey("CPMCellClassKey")
        obj.cellClass = obj.hint.objectForKey("CPMCellClassKey");
    end
    
    obj.subcellIndices = obj.hint.objectForKey("CPMSubcellIndicesKey");
    obj.subcellConverter = obj.hint.objectForKey("CPMSubcellConverterKey");
    
    obj.contactEnergy = obj.hint.objectForKey("CPMContactEnergyKey");
    obj.surfaceYModulusArray = ...
        obj.hint.objectForKey("CPMSurfaceYModulusArrayKey");
    obj.areaBModulusArray = ...
        obj.hint.objectForKey("CPMAreaBModulusArrayKey");
    obj.HamiltonianSwitch = ...
        obj.hint.objectForKey("CPMHamiltonianSwitchKey");
    obj.T = obj.hint.objectForKey("CPMTKey");
    
    obj.drawingSwitch = obj.hint.objectForKey("CPMDrawingSwitchKey");
    obj.lut = obj.hint.objectForKey("CPMLutKey");
    
    obj.outputName = obj.hint.objectForKey("CPMOutputNameKey");
    if isstring(obj.outputName)
        obj.outputName = convertStringsToChars(obj.outputName);
    end
    if ~isnan(obj.outputName)
        obj.frameInterval = obj.hint.objectForKey("CPMFrameIntervalKey");
        obj.stackLimit = obj.hint.objectForKey("CPMStackLimitKey");
    else
        obj.frameInterval = -1;
    end
end
function prepareCells(obj,labels)
% Method to initialize cells.
% prepareCells(obj,labels)
% Argument labels is an array of cell labels corresponding to labels in
%   cell map.
    % Allocate cells.
    obj.makeCells(labels);
    
    % Fit a frame to a cell.
    obj.enframeCells;
    
    % Set cell type and variables dependent on the cell types.
    obj.makeCellsReadCelltype;
    
    % Mark surface.
    obj.markSubcellularLocations;
    obj.prepareSurfaceElasticity;
    
    % Mark area.
    obj.markSubcellularComponents;
    obj.prepareAreaConstraint;
    
    % Get energy of the cells.
    obj.makeCellsRegisterEnergy;
    
    % Draw cells.
    obj.cellArray.makeObjectsPerformSelector(@drawCell);
end
function makeCells(obj,labels)
% Method to allocate cells.
% makeCells(obj,labels)
% Argument labels is an array of cell labels corresponding to labels in
%   the cellMap.
    fh = str2func(obj.cellClass);
    for i = 1:length(labels)
        cell = fh();
        cell.owner = obj;
        cell.cellIndex = i;
        cell.label = labels(i);
        
        obj.cellArray.addObject(cell);
    end
end
function enframeCells(obj)
% Method to make cells find enclosing frame.
% enframeCells(obj)
    obj.cellArray.makeObjectsPerformSelector(@enframe);
end
function makeCellsReadCelltype(obj)
% Method to make cells set cell type.
% makeCellsReadCelltype(obj)
    obj.cellArray.makeObjectsPerformSelector(@readCellType);
end
function markSubcellularLocations(obj)
% Method to make cells mark subcellular locations.
% markSubcellularLocations(obj)
    obj.cellArray.makeObjectsPerformSelector(@markSurface);
end
function prepareSurfaceElasticity(obj)
% Method to make cells initialize surface length reference values.
% prepareSurfaceElasticity(obj)
    obj.cellArray.makeObjectsPerformSelector(@makePoCurrentLength);
end
function markSubcellularComponents(obj)
% Method to make cells mark pixels to be updated.
% markSubcellularComponents(obj)
    obj.cellArray.makeObjectsPerformSelector(@collectRim);
end
function prepareAreaConstraint(obj)
% Method to make cells initialize area reference values.
% prepareAreaConstraint(obj)
    obj.cellArray.makeObjectsPerformSelector(@makeAoCurrentVolume);
end
function makeCellsRegisterEnergy(obj)
% Method to make cells calculate and record energy.
% makeCellsRegisterEnergy(obj)
    obj.cellArray.makeObjectsPerformSelector(@registerEnergy);
end

function result = cellWithLabel(obj,label)
% Method to get a cell by a label.
% result = cellWithLabel(obj,label)
% Argument label is a value representing the cell in the cellMap.
% Return value is an instance of CPMCell family.
    index = find(obj.cellLabels == label,1);
    if isempty(index)
        result = nan;
    else
        result = obj.cellArray.objectAtIndex(index);
    end
end

function runSimulation(obj)
% Method to run a simulation.
% runSimulation(obj)
    obj.prepareSimulation;
    obj.endFlag = false;
    
    obj.startSimulation;
end
function prepareSimulation(obj)
% Method to prepare a simulation.
% prepareSimulation(obj)
% This method is called by runSimulation() before starting the simulation.
    if ~isnan(obj.outputName)
        data = obj.hint.data;
        data.writeToFile([obj.outputName,'.hint.mat']);
        text_ = obj.hint.description;
        fid = fopen([obj.outputName,'.hint.txt'],'w');
        fprintf(fid,text_);
        fclose(fid);
        image = SYImage(obj.imageData.var);
        image.writeToFile([obj.outputName,'.tif'],false);
    end
    
    obj.intervalCounter = obj.frameInterval;
    obj.frameCounter = 0;
    
    obj.prepareWindow;
end
function startSimulation(obj)
% Method to start simulating.
% startSimulation(obj)
% This method is called by runSimulation() and implement a loop of updates
%   until end.
    obj.stopFlag = false;
    
    while true
        if obj.tryUpdate
            obj.advanceCounters;
        end
        
        if obj.stopFlag
            break;
        end
    end
    
    if obj.endFlag
        obj.endSimulation;
    end
end
function pauseSimulation(obj)
% Method to pause a running simulation.
% pauseSimulation(obj)
    obj.stopFlag = true;
end
function endSimulation(obj)
% Method to terminate a simulation and write its result.
% endSimulation(obj)
    if ~isnan(obj.outputName)
        str = [obj.outputName,'.end.mat'];
        data = obj.data;
        data.writeToFile(str);
    end
    
    if ishandle(obj.window)
        delete(obj.window);
    end
end
function rerunSimulation(obj)
% Method to rerun a simulation.
% rerunSimulation(obj)
    str = [obj.outputName,'.end.mat'];
    if exist(str,'file') == 2
        ttr = [obj.outputName,'midway.mat'];
        if exist(ttr,'file') == 2
            delete(ttr)
        end
        movefile(str,ttr);
    end

    obj.intervalCounter = obj.frameInterval;
    obj.frameCounter = 0;
    obj.endFlag = false;

    obj.startSimulation;
end

function prepareWindow(obj)
% Method to initialize window.
% prepareWindow(obj)
    obj.window = figure;
    obj.textField_progress = uicontrol(obj.window,'Style','text', ...
        'Position',[5 35 480 20],'FontSize',14);
    obj.imageView = imshow(zeros(obj.map.frameSize));
    obj.window.Pointer = 'crosshair';
end

function result = tryUpdate(obj)
% Method to implement Metropolis-step.
% result = tryUpdate(obj)
% Return value is a bool indicating if the update was accepted.
    obj.selectPoint;
    
    obj.updateCellsInList;
    
    obj.dH = obj.dHofCells;
    obj.p = exp(-obj.dH / obj.T);
    if obj.p > rand()
        obj.fixUpdate;
        obj.updateView;
        
        result = true;
    else
        obj.revert;
        
        result = false;
    end
end
function selectPoint(obj)
% Method to make a cell raise a candidate of change.
% selectPoint(obj)
    while true
        % Select a boundary point randomly.
        array = cumsum(obj.rimCounts.var);
        i = ceil(rand() * array(end));
        j = find(array >= i,1);
        
        % Select a label to update.
        cell = obj.cellArray.objectAtIndex(j);
        if cell.tryUpdate
            break;
        end
    end
end
function fixUpdate(obj)
% Method to make cells accept the update.
% fixUpdate(obj)
    for i = obj.updateList
        cell = obj.cellWithLabel(i);
        cell.fixUpdate;
        cell.drawCell;
    end
end
function updateView(obj)
% Method to update window for an update.
% updateView(obj)
    % Draw map.
    if ishandle(obj.window)
        figure(obj.window);
        obj.imageView.CData = obj.imageData.var;
        H = sum(obj.cellEnergy.var);
        obj.textField_progress.String = ['H:',num2str(H,'%.4g'),...
            ' dH:',num2str(obj.dH,'%.4g'),...
            ' P:',num2str(obj.p,'%.4g'),...
            ' fi:',num2str(obj.intervalCounter)];
        drawnow;
    end
end
function advanceCounters(obj)
% Method to count updates.
% advanceCounters(obj)
    if obj.intervalCounter > 0
        obj.intervalCounter = obj.intervalCounter - 1;
    elseif obj.intervalCounter == 0
        str = [obj.outputName,'.tif'];
        image = SYImage(str);
        image.addRepresentation(obj.imageData.var);
        image.writeToFile(str,false);
        
        obj.intervalCounter = obj.frameInterval;
        
        obj.frameCounter = obj.frameCounter + 1;
        if obj.frameCounter >= obj.stackLimit
            obj.stopFlag = true;
            obj.endFlag = true;
        end
    end
end
function revert(obj)
% Method to reject an update.
% revert(obj)
    for i = obj.updateList
        cell = obj.cellWithLabel(i);
        cell.revert;
    end
    obj.map.revert;
end
function result = dHofCells(obj)
% Method to get change of energy.
% result = dHofCells(obj)
% Return value is the change of energy.
    dH_ = 0;
    for i = obj.updateList
        cell = obj.cellWithLabel(i);
        dH_ = dH_ + cell.dH;
    end
    result = dH_;
end
function updatePixel(obj,pixel,label,neighbors)
% Method to receive a candidate of change.
% updatePixel(obj,pixel,label,neighbors)
% Argument pixel is a point to be updated.
% Argument label is an arry of cell label, subcellular component label,
%   and cell type label to be pasted.
% Argument neighbors is an array of cell labels around the pixel.
    obj.map.updatePixel(pixel,label);
    
    array = zeros(1,length(neighbors(:)));
    array(:) = neighbors(:);
    
    obj.updateList = array;
end
function updateCellsInList(obj)
% Method to make cells around the changed point receive a change.
% updateCellsInList(obj)
    for i = obj.updateList
        cell = obj.cellWithLabel(i);
        cell.update;
    end
end

end
end
