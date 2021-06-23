# Cellular Potts model simulator.
This project provides a framework for cellular Potts model simulation.

## Installation
Download files and put them in a folder with a suitable name. Go to Matlab command line, enter "addpath" + a full path of the folder, and enter "savepath".

## Requirement
This project requires no Matlab toolbox, but a custom framework of objective classes SYObject family which is available at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3625631.svg)](https://doi.org/10.5281/zenodo.3625631).

## Example
Below is an example of simulation in which a tissue with two types of cells is compressed. The cells have different surface tensions. It first sets parameters, prepares an initial configuration of cells, and then simulates motion of the cells.

### Parameters
The parameters are set in a dictionary ‘hint’ and saved in a folder ‘properties’.
To keep reference surface area before and after the compression, CPMCACell its used.
It defines three cell types and three subcellular components.
The variable subcellIndices is an numerical representation of the subcellular components.
```
mkdir('properties');
 
hint = SYDictionary;
 
hint.setObjectForKey('CPMCellClassKey','CPMCACell');
 
%% Prepare labels.
celltypeArray = SYArray('media','softCell','stiffCell');
labelArray = SYArray('media','cytosol_soft','cytosol_stiff');
subcellIndices = SYData((1:labelArray.count)');
 
hint.setObjectForKey('CPMCelltypeArrayKey',celltypeArray);
hint.setObjectForKey('CPMLabelArrayKey',labelArray);
hint.setObjectForKey('CPMSubcellIndicesKey',subcellIndices);
```

For a Metropolis’s step, a label is not copied directly, but copied from a matrix by looking up labels on a pixel to be updated and a neighbor pixel. It also depends on whether the neighbor pixel is in the same cell or not.
```
%% Prepare label converter.
c = labelArray.count;
cc = nan(c,c,2);
% intracellular.
z = 1;
% media : ~.
label = 'media'; % Label on a pixel to be updated.
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject('media'),z) = ...
%     labelArray.indexOfObject('media');
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('cytosol_soft'),z) = ... Label on a pixel to be referred.
    labelArray.indexOfObject('cytosol_soft'); % Label to be pasted.
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('cytosol_stiff'),z) = ...
    labelArray.indexOfObject('cytosol_stiff');
% cytosol_soft : ~.
label = 'cytosol_soft';
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('media'),z) = ...
    labelArray.indexOfObject('media');
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject('cytosol_soft'),z) = ...
%     labelArray.indexOfObject('cytosol_soft');
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('cytosol_stiff'),z) = ...
    labelArray.indexOfObject('cytosol_stiff');
% cytosol_stiff : ~.
label = 'cytosol_stiff';
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('media'),z) = ...
    labelArray.indexOfObject('media');
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('cytosol_soft'),z) = ...
    labelArray.indexOfObject('cytosol_soft');
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject('cytosol_stiff'),z) = ...
%     labelArray.indexOfObject('cytosol_stiff');
 
% intercellular.
z = 2;
% media : ~.
label = 'media';
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('media'),z) = ...
    labelArray.indexOfObject('media');
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('cytosol_soft'),z) = ...
    labelArray.indexOfObject('cytosol_soft');
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('cytosol_stiff'),z) = ...
    labelArray.indexOfObject('cytosol_stiff');
% cytosol_soft : ~.
label = 'cytosol_soft';
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('media'),z) = ...
    labelArray.indexOfObject('media');
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('cytosol_soft'),z) = ...
    labelArray.indexOfObject('cytosol_soft');
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('cytosol_stiff'),z) = ...
    labelArray.indexOfObject('cytosol_stiff');
% cytosol_stiff : ~.
label = 'cytosol_stiff';
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('media'),z) = ...
    labelArray.indexOfObject('media');
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('cytosol_soft'),z) = ...
    labelArray.indexOfObject('cytosol_soft');
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('cytosol_stiff'),z) = ...
    labelArray.indexOfObject('cytosol_stiff');
 
%
cMatrix = CPMConverterMatrix;
cMatrix.labelArray = labelArray;
cMatrix.converterMatrix = cc;
data = cMatrix.data;
data.writeToFile('properties/cMatrix.mat');
 
cc = cMatrix.converterWithLabels(labelArray);
hint.setObjectForKey('CPMSubcellConverterKey',cc);
```

The area constraint is counted for a total area of the cell and each subcellular component in it.
```
%% Prepare area constraint.
array = SYArray;
% area constraint arrays.
% Vulk modulus for [total area; media; cytosol_soft; cytosol_stiff].
% media.
lambda = [1; 0; 0; 0];
array.addObject(lambda);
% soft cell.
lambda = [1; 0; 0; 0];
array.addObject(lambda);
% stiff cell.
lambda = [1; 0; 0; 0];
array.addObject(lambda);
 
acMatrix = CPMAreaConstraintMatrix;
acMatrix.labelArray = labelArray;
acMatrix.celltypeArray = celltypeArray;
acMatrix.lambdaArray = array;
data = acMatrix.data;
data.writeToFile('properties/acMatrix.mat');
 
av = acMatrix.areaVModulusArrayWithCellTypesAndLabels(celltypeArray, ...
    labelArray);
hint.setObjectForKey('CPMAreaBModulusArrayKey',av);
```

The surface tension includes a surface contractility and surface elasticity.
The surface contractility is defined by contact energy between subcellular components either inside or between cells.
```
%% Prepare surface tension.
% Perimeter contractility.
baseTension = [0.0,1.0];
 
c = labelArray.count;
cm = zeros(c,c,2);
% intracellular.
z = 1;
% media : ~.
label = 'media'; % Label on a pixel of interest.
% cm(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject('media'),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('cytosol_soft'),z) = 1.0; % Label on an adjacent pixel.
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('cytosol_stiff'),z) = 1.0;
% cytosol_soft : ~.
label = 'cytosol_soft';
% cm(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject('cytosol_soft'),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('cytosol_stiff'),z) = 4.0;
% cytosol_stiff : ~.
% label = 'cytosol_stiff';
% cm(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject('cytosol_stiff'),z) = 1.0;
 
% intercellular.
z = 2;
% media : ~.
label = 'media';
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('media'),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('cytosol_soft'),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('cytosol_stiff'),z) = 1.0;
% cytosol_soft : ~.
label = 'cytosol_soft';
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('cytosol_soft'),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('cytosol_stiff'),z) = 4.0;
% cytosol_stiff : ~.
label = 'cytosol_stiff';
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject('cytosol_stiff'),z) = 4.0;
 
cm = ss_symmetrize_upper_triangle(cm);
```
 
The surface elasticity is defined for cell types.
```
% perimeter elasticity.
array = SYArray;
% [total].
% media.
lambda = [0];
array.addObject(SYData(lambda));
% softCell.
lambda = [0];
array.addObject(SYData(lambda));
% stiffCell.
lambda = [0];
array.addObject(SYData(lambda));
%
surfaceYModulusArray = array;
 
stMatrix = CPMSurfaceTensionMatrix;
stMatrix.labelArray = labelArray;
stMatrix.celltypeArray = celltypeArray;
stMatrix.baseTension = baseTension;
stMatrix.contactEnergyMatrix = cm;
stMatrix.surfaceYModulusArray = surfaceYModulusArray;
data = stMatrix.data;
data.writeToFile('properties/stMatrix.polar.mat');
 
ce = stMatrix.contactEnergyWithLabels(labelArray);
hint.setObjectForKey('CPMContactEnergyKey',ce);
sm = stMatrix.surfaceYModulusArrayWithCellTypes(celltypeArray);
hint.setObjectForKey('CPMSurfaceYModulusArrayKey',sm);
```

Switches categorize cell types for various process such as drawing.
```
%% Prepare swiths.
% Hamiltonian switch.
m = zeros(2,2);
%  No energy.
m(1,1) = celltypeArray.indexOfObject('media');
%  With energy.
m(1,2) = celltypeArray.indexOfObject('softCell');
m(2,2) = celltypeArray.indexOfObject('stiffCell');
HamiltonianSwitch = SYData(m);
 
% drawing switch.
m = zeros(2,2);
%  Gray rim.
m(1,1) = celltypeArray.indexOfObject('media');
%  Colored rim.
m(1,2) = celltypeArray.indexOfObject('softCell');
m(2,2) = celltypeArray.indexOfObject('stiffCell');
drawingSwitch = SYData(m);
 
 
array = SYArray(HamiltonianSwitch,drawingSwitch);
data = array.data;
data.writeToFile('properties/switches.mat');
 
hint.setObjectForKey('CPMHamiltonianSwitchKey',array.objectAtIndex(1));
hint.setObjectForKey('CPMDrawingSwitchKey',array.objectAtIndex(2));
```

It needs also several parameters such as temperature and duration.
```
% tempareture.
hint.setObjectForKey('CPMTKey',2);
 
% color space.
lut = [255,255,255; ... media
    220,220,220; ... cytosol_soft
    220,220,0; ... cytosol_stiff
    
    255,255,255; ... clear-color.
    
    126,245,126; ... lateral rim. Second from bottom.
    126,126,126 ... media rim. Fist from bottom.
    ];
hint.setObjectForKey('CPMLutKey',SYData(lut));
 
% frame interval.
hint.setObjectForKey('CPMFrameIntervalKey',5000);
 
% stack limit.
hint.setObjectForKey('CPMStackLimitKey',12);
 
% export.
data = hint.data;
data.writeToFile('properties/hint.mat');
txt = hint.description;
fid = fopen('properties/hint.txt','w');
fprintf(fid,txt);
fclose(fid);
```

### Initial configuration
The initial configuration is prepared by Voronoi tessellation from random sparse points.
```
%% Prepare a map of compressed tissue.
frame = [270,480];
c_max = 600;
r_min = 10;
comp_rate = 1.8;
margin = 40;

Frame = [frame(1) / comp_rate,frame(2) * comp_rate];
seeds = ss_random_sparse(Frame,r_min,c_max);
data = SYData(seeds);
data.writeToFile('seeds.mat');

upperCells = seeds(:,2) < Frame(1) / 2;
leftCells = seeds(:,1) < Frame(2) / 2;

bitmap = uint16(ss_voronoi_torus(seeds,Frame));
citmap = zeros(Frame,'uint16');
ditmap = zeros(Frame,'uint16');

label = subcellIndices.var(labelArray.indexOfObject('cytosol_soft'));
citmap(:) = label;
celltype = celltypeArray.indexOfObject('softCell');
ditmap(:) = celltype;

label = subcellIndices.var(labelArray.indexOfObject('cytosol_stiff'));
celltype = celltypeArray.indexOfObject('stiffCell');
indices = find(seeds(:,1) > Frame(2) / 2);
for i = indices(:)'
    mask = bitmap == i;
    citmap(mask) = label;
    ditmap(mask) = celltype;
end

map = CPMMap;
map.initWithMaps(bitmap,citmap,ditmap,hint);
data = map.data;
data.writeToFile('map0.mat');
```

Then the tissue is compressed.
``` 
% Reference surface area and volume are saved.
simulator = CPMSimulator;
simulator.initWithMap(map,hint);
cells = simulator.cellArray;
[labels,AoPo] = ss_export_AoPo(cells);

% Images are resized.
image = SYImage(SYData(map.cellMap));
image.frameSize = frame;
cellMap = image.drawBitmapRep;

image = SYImage(SYData(map.subcellMap));
image.frameSize = frame;
subcellMap = image.drawBitmapRep;

image = SYImage(SYData(map.cellTypeMap));
image.frameSize = frame;
celltypeMap = image.drawBitmapRep;

map = CPMMap(cellMap,subcellMap,celltypeMap,hint);

% Saved reference value are assigned to cells.
simulator = CPMSimulator;
simulator.initWithMap(map,hint);
cells = simulator.cellArray;
for j = 1:cells.count
    cel = cells.objectAtIndex(j);
    cel.importAoPo(labels,AoPo);
end
```

### Running simulation
Finally, it simulates relaxation of the cells in the compressed tissue. Images of cells and cell segmentation are exported when runSimulation() terminates. Therefore it iterates for several times.
```
mkdir('segmentations');
mkdir('cells');

iterN = 15;
for i = 1:iterN
    str = ['relaxation.t',num2str(i)];
    hint.setObjectForKey('CPMOutputNameKey',str);
    simulator.readHint;
    
    simulator.runSimulation;

    citmap = ss_unroll_cellLabel(map,cells,labels,upperCells, ...
        leftCells,frame + margin);
    citmap = citmap + 1;
    witmap = ss_erode_edge(map,cells);
    witmap = ss_unroll_watershed(map,cells,labels,upperCells, ...
        leftCells,witmap,frame + margin);
    citmap(witmap == 0) = 0;
    image = SYImage(SYData(uint16(citmap)));
    image.writeToFile(['cells/',str,'.tif'],false);
    image = SYImage(SYData(logical(witmap)));
    image.writeToFile(['segmentations/',str,'.png'],false);
end
```
