# Cellular Potts model simulator for an epithelium
This project provides a framework for cellular Potts model simulation of an epithelial tissue.

## Installation
Download files and put them in a folder with a suitable name. Go to Matlab command line, enter "addpath" + a full path of the folder, and enter "savepath".

## Requirement
This project requires no Matlab toolbox, but a custom framework of objective classes SYObject family which is available at [![DOI](https://zenodo.org/badge/235579182.svg)](https://zenodo.org/badge/latestdoi/235579182).

## Example
Below is an example of simulation in which the epithelium was invaginated by a modulated surface elasticity and a supracellular myosin cable.
In this simulation framework, a cell is assigned a cell type and composed of materials. The cell occupies a connected area, while the materials inside the cell can be dispersed disconnected.

### Setting parameters
The parameters are stored in a dictionary ‘hint’.
First, labels representing cell types and materials are prepared.
The variable subcellIndices is a numerical representation of the materials.
```
hint = SYDictionary;

%% prepare labels
celltypeArray = SYArray("media","constrictingCell","softCell", ...
    "stiffCell","cableCell","staticCell");
labelArray = SYArray("bodyFluid","extraBodyFluid","apical_ECM", ...
    "cytosol_const","cytosol_soft","cytosol_stiff","cytosol_static");
subcellIndices = SYData((1:labelArray.count)');

hint.setObjectForKey("CPMCelltypeArrayKey",celltypeArray);
hint.setObjectForKey("CPMLabelArrayKey",labelArray);
hint.setObjectForKey("CPMSubcellIndicesKey",subcellIndices);

%% hints for map.
% label dict.
labelDict = SYDictionary;
indices = [labelArray.indexOfObject("extraBodyFluid"), ...
           labelArray.indexOfObject("apical_ECM")];
labelDict.setObjectForKey("EC-apical",subcellIndices.var(indices));
labelDict.setObjectForKey("celltype_media", ...
    celltypeArray.indexOfObject("media")); % dimension in [1,1,n].
hint.setObjectForKey("CPMLabelDictKey",labelDict);
```

For an update of labeling in the Metropolis’s step, a label is not copied directly from a neighbor, but copied from a matrix by looking up labels on a site to be updated and the neighbor site. It also depends on whether the neighbor site is in the same cell or not.
```
%% prepare converter
c = labelArray.count;
cc = nan(c,c,2);
%% intracellular.
z = 1;
% body fluid : ~.
label = "bodyFluid";
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("bodyFluid"),z) = ...
%     labelArray.indexOfObject("bodyFluid");
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("extraBodyFluid"),z) = ...
%     labelArray.indexOfObject("extraBodyFluid");
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("apical_ECM"),z) = ...
%     labelArray.indexOfObject("apical_ECM");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_const"),z) = ...
    labelArray.indexOfObject("cytosol_const");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_soft"),z) = ...
    labelArray.indexOfObject("cytosol_soft");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_stiff"),z) = ...
    labelArray.indexOfObject("cytosol_stiff");
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("cytosol_static"),z) = ...
%     labelArray.indexOfObject("cytosol_static");
% external body fluid : ~.
label = "extraBodyFluid";
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("bodyFluid"),z) = ...
%     labelArray.indexOfObject("bodyFluid");
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("extraBodyFluid"),z) = ...
%     labelArray.indexOfObject("extraBodyFluid");
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("apical_ECM"),z) = ...
%     labelArray.indexOfObject("apical_ECM");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_const"),z) = ...
    labelArray.indexOfObject("cytosol_const");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_soft"),z) = ...
    labelArray.indexOfObject("cytosol_soft");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_stiff"),z) = ...
    labelArray.indexOfObject("cytosol_stiff");
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("cytosol_static"),z) = ...
%     labelArray.indexOfObject("cytosol_static");
% apical ECM : ~.
% % static.
% cytosol constricting : ~.
label = "cytosol_const";
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("bodyFluid"),z) = ...
    labelArray.indexOfObject("bodyFluid");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("extraBodyFluid"),z) = ...
    labelArray.indexOfObject("extraBodyFluid");
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("apical_ECM"),z) = ...
%     labelArray.indexOfObject("apical_ECM");
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("cytosol_const"),z) = ...
%     labelArray.indexOfObject("cytosol_const");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_soft"),z) = ...
    labelArray.indexOfObject("cytosol_soft");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_stiff"),z) = ...
    labelArray.indexOfObject("cytosol_stiff");
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("cytosol_static"),z) = ...
%     labelArray.indexOfObject("cytosol_static");
% cytosol soft : ~.
label = "cytosol_soft";
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("bodyFluid"),z) = ...
    labelArray.indexOfObject("bodyFluid");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("extraBodyFluid"),z) = ...
    labelArray.indexOfObject("extraBodyFluid");
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("apical_ECM"),z) = ...
%     labelArray.indexOfObject("apical_ECM");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_const"),z) = ...
    labelArray.indexOfObject("cytosol_const");
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("cytosol_soft"),z) = ...
%     labelArray.indexOfObject("cytosol_soft");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_stiff"),z) = ...
    labelArray.indexOfObject("cytosol_stiff");
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("cytosol_static"),z) = ...
%     labelArray.indexOfObject("cytosol_static");
% cytosol stiff : ~.
label = "cytosol_stiff";
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("bodyFluid"),z) = ...
    labelArray.indexOfObject("bodyFluid");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("extraBodyFluid"),z) = ...
    labelArray.indexOfObject("extraBodyFluid");
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("apical_ECM"),z) = ...
%     labelArray.indexOfObject("apical_ECM");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_const"),z) = ...
    labelArray.indexOfObject("cytosol_const");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_soft"),z) = ...
    labelArray.indexOfObject("cytosol_soft");
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("cytosol_stiff"),z) = ...
%     labelArray.indexOfObject("cytosol_stiff");
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("cytosol_static"),z) = ...
%     labelArray.indexOfObject("cytosol_static");
% cytosol static : ~.
% % static.
% PCP upstream boundary : ~.
% % static.

%% intercellular.
z = 2;
% body fluid : ~.
label = "bodyFluid";
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("bodyFluid"),z) = ...
    labelArray.indexOfObject("bodyFluid");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("extraBodyFluid"),z) = ...
    labelArray.indexOfObject("extraBodyFluid");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("apical_ECM"),z) = ...
    labelArray.indexOfObject("extraBodyFluid");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_const"),z) = ...
    labelArray.indexOfObject("cytosol_const");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_soft"),z) = ...
    labelArray.indexOfObject("cytosol_soft");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_stiff"),z) = ...
    labelArray.indexOfObject("cytosol_stiff");
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("cytosol_static"),z) = ...
%     labelArray.indexOfObject("cytosol_static");
% external body fluid : ~.
label = "extraBodyFluid";
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("bodyFluid"),z) = ...
    labelArray.indexOfObject("bodyFluid");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("extraBodyFluid"),z) = ...
    labelArray.indexOfObject("extraBodyFluid");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("apical_ECM"),z) = ...
    labelArray.indexOfObject("extraBodyFluid");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_const"),z) = ...
    labelArray.indexOfObject("cytosol_const");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_soft"),z) = ...
    labelArray.indexOfObject("cytosol_soft");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_stiff"),z) = ...
    labelArray.indexOfObject("cytosol_stiff");
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("cytosol_static"),z) = ...
%     labelArray.indexOfObject("cytosol_static");
% apical ECM : ~.
% % static.
% cytosol constricting : ~.
label = "cytosol_const";
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("bodyFluid"),z) = ...
    labelArray.indexOfObject("bodyFluid");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("extraBodyFluid"),z) = ...
    labelArray.indexOfObject("extraBodyFluid");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("apical_ECM"),z) = ...
    labelArray.indexOfObject("extraBodyFluid");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_const"),z) = ...
    labelArray.indexOfObject("cytosol_const");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_soft"),z) = ...
    labelArray.indexOfObject("cytosol_soft");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_stiff"),z) = ...
    labelArray.indexOfObject("cytosol_stiff");
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("cytosol_static"),z) = ...
%     labelArray.indexOfObject("cytosol_static");
% cytosol soft : ~.
label = "cytosol_soft";
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("bodyFluid"),z) = ...
    labelArray.indexOfObject("bodyFluid");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("extraBodyFluid"),z) = ...
    labelArray.indexOfObject("extraBodyFluid");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("apical_ECM"),z) = ...
    labelArray.indexOfObject("extraBodyFluid");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_const"),z) = ...
    labelArray.indexOfObject("cytosol_const");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_soft"),z) = ...
    labelArray.indexOfObject("cytosol_soft");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_stiff"),z) = ...
    labelArray.indexOfObject("cytosol_stiff");
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("cytosol_static"),z) = ...
%     labelArray.indexOfObject("cytosol_static");
% cytosol stiff : ~.
label = "cytosol_stiff";
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("bodyFluid"),z) = ...
    labelArray.indexOfObject("bodyFluid");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("extraBodyFluid"),z) = ...
    labelArray.indexOfObject("extraBodyFluid");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("apical_ECM"),z) = ...
    labelArray.indexOfObject("extraBodyFluid");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_const"),z) = ...
    labelArray.indexOfObject("cytosol_const");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_soft"),z) = ...
    labelArray.indexOfObject("cytosol_soft");
cc(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_stiff"),z) = ...
    labelArray.indexOfObject("cytosol_stiff");
% cc(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("cytosol_static"),z) = ...
%     labelArray.indexOfObject("cytosol_static");
% cytosol static : ~.
% % static.

hint.setObjectForKey("CPMSubcellConverterKey", SYData(cc));
```

The surface tension includes a surface contractility and a surface elasticity.
The surface contractility is defined by a contact energy between materials, and it also depends on either the materials are in a cell or different cells.
The surface elasticity is defined by an array of elastic moduli for each cell type, and the moduli are for total, apical, lateral, basal, upstream, and downstream perimeters.
The upstream and downstream perimeter are defined for a cell with planar cell polarity, but are not used in this example.
```
%% prepare contact energy
c = labelArray.count;
cm = zeros(c,c,2);
%% intracellular.
z = 1;
% body fluid : ~.
label = "bodyFluid";
% cm(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("bodyFluid"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("extraBodyFluid"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("apical_ECM"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_const"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_soft"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_stiff"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_static"),z) = 1.0;
% external body fluid : ~.
label = "extraBodyFluid";
% cm(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("extraBodyFluid"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("apical_ECM"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_const"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_soft"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_stiff"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_static"),z) = 1.0;
% apical ECM : ~.
label = "apical_ECM";
% cm(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("apical_ECM"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_const"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_soft"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_stiff"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_static"),z) = 1.0;
% cytosol constricting : ~.
label = "cytosol_const";
% cm(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("cytosol_const"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_soft"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_stiff"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_static"),z) = 1.0;
% cytosol soft : ~.
label = "cytosol_soft";
% cm(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("cytosol_soft"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_stiff"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_static"),z) = 1.0;
% cytosol stiff : ~.
label = "cytosol_stiff";
% cm(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("cytosol_stiff"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_static"),z) = 1.0;
% cytosol static : ~.
label = "cytosol_static";
% cm(labelArray.indexOfObject(label), ...
%     labelArray.indexOfObject("cytosol_static"),z) = 1.0;


%% intercellular.
z = 2;
% body fluid : ~.
label = "bodyFluid";
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("bodyFluid"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("extraBodyFluid"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("apical_ECM"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_const"),z) = 4.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_soft"),z) = 4.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_stiff"),z) =4.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_static"),z) = 1.0;
% external body fluid : ~.
label = "extraBodyFluid";
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("extraBodyFluid"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("apical_ECM"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_const"),z) = 4.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_soft"),z) = 4.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_stiff"),z) = 4.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_static"),z) = 1.0;
% apical ECM : ~.
label = "apical_ECM";
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("apical_ECM"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_const"),z) = 4.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_soft"),z) = 4.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_stiff"),z) = 4.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_static"),z) = 1.0;
% cytosol constricting : ~.
label = "cytosol_const";
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_const"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_soft"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_stiff"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_static"),z) = 1.0;
% cytosol soft : ~.
label = "cytosol_soft";
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_soft"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_stiff"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_static"),z) = 1.0;
% cytosol stiff : ~.
label = "cytosol_stiff";
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_stiff"),z) = 1.0;
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_static"),z) = 1.0;
% cytosol static : ~.
label = "cytosol_static";
cm(labelArray.indexOfObject(label), ...
    labelArray.indexOfObject("cytosol_static"),z) = 1.0;

cm = ss_symmetrize_upper_triangle(cm);
hint.setObjectForKey("CPMContactEnergyKey", SYData(cm));


%% perimeter elascticity
array = SYArray;
% [total, apical, lateral, basal, upsteram, downstream].
% media.
lambda = [0; 0; 0; 0; 0; 0];
array.addObject(SYData(lambda));
% constrictingCell.
lambda = [0; 5; 5; 5; 0; 0];
array.addObject(SYData(lambda));
% softCell.
lambda = [0; 0.5; 0.5; 0.5; 0; 0];
array.addObject(SYData(lambda));
% stiffCell.
lambda = [0; 5; 5; 5; 0; 0];
array.addObject(SYData(lambda));
% cableCell.
lambda = [0; 5; 5; 5; 0; 0];
array.addObject(SYData(lambda));
% staticCell.
lambda = [0; 0; 0; 0; 0; 0];
array.addObject(SYData(lambda));

hint.setObjectForKey("CPMSurfaceYModulusArrayKey", array);
```

The area constraint is defined for each cell type, and is counted for a total area of a cell and each material in the cell.
```
%% area constraint
array = SYArray;
% area constraint arrays.
% [total area, bodyFluid, extraBodyFluid, apical_ECM, cytosol_const, 
%   cytosol_soft, cytosol_stiff, cytosol_static].
% media.
lambda = [0; 0; 0; 0; 0; 0; 0; 0];
array.addObject(SYData(lambda));
% constrictingCell.
lambda = [0; 0; 0; 0; 1; 1; 1; 0];
array.addObject(SYData(lambda));
% softCell.
lambda = [0; 0; 0; 0; 1; 1; 1; 0];
array.addObject(SYData(lambda));
% stiffCell.
lambda = [0; 0; 0; 0; 1; 1; 1; 0];
array.addObject(SYData(lambda));
% cableCell.
lambda = [0; 0; 0; 0; 1; 1; 1; 0];
array.addObject(SYData(lambda));
% staticCell.
lambda = [0; 0; 0; 0; 0; 0; 0; 0];
array.addObject(SYData(lambda));

hint.setObjectForKey("CPMAreaBModulusArrayKey", array);
```

Switches categorize cell types for various processes.
Hamiltonian switch indicates if a cell type is subject to a calculation of energy.
Drawing switch specifies how a cell is drawn.
Polarity switch indicates if a cell is apico-basally polarized or not.
Epithelial switch indicates if a cell type is epithelial or not.
PCP switch indicates if a cell type is plenarily polarized or not.
```
%% cell types categorization
% Hamiltonian switch.
m = zeros(4,2);
%  No energy.
m(1,1) = celltypeArray.indexOfObject("media");
m(2,1) = celltypeArray.indexOfObject("staticCell");
%  With energy.
m(1,2) = celltypeArray.indexOfObject("constrictingCell");
m(2,2) = celltypeArray.indexOfObject("softCell");
m(3,2) = celltypeArray.indexOfObject("stiffCell");
m(4,2) = celltypeArray.indexOfObject("cableCell");
hint.setObjectForKey("CPMHamiltonianSwitchKey", SYData(m));

% drawing switch.
m = zeros(4,2);
%  Gray rim.
m(1,1) = celltypeArray.indexOfObject("media");
m(2,1) = celltypeArray.indexOfObject("staticCell");
%  Colored rim.
m(1,2) = celltypeArray.indexOfObject("constrictingCell");
m(2,2) = celltypeArray.indexOfObject("softCell");
m(3,2) = celltypeArray.indexOfObject("stiffCell");
m(4,2) = celltypeArray.indexOfObject("cableCell");
hint.setObjectForKey("CPMDrawingSwitchKey", SYData(m));

% polarity switch.
m = zeros(5,2);
%  non-polarized cell.
m(1,1) = celltypeArray.indexOfObject("media");
%  polarized cell.
m(1,2) = celltypeArray.indexOfObject("constrictingCell");
m(2,2) = celltypeArray.indexOfObject("softCell");
m(3,2) = celltypeArray.indexOfObject("stiffCell");
m(4,2) = celltypeArray.indexOfObject("cableCell");
m(5,2) = celltypeArray.indexOfObject("staticCell");
hint.setObjectForKey("CPMPolaritySwitchKey", SYData(m));

% epithelial switch.
m = zeros(4,2);
%  non-epithelial cell.
m(1,1) = celltypeArray.indexOfObject("media");
m(2,1) = celltypeArray.indexOfObject("staticCell");
%  epithelial cell.
m(1,2) = celltypeArray.indexOfObject("constrictingCell");
m(2,2) = celltypeArray.indexOfObject("softCell");
m(3,2) = celltypeArray.indexOfObject("stiffCell");
m(4,2) = celltypeArray.indexOfObject("cableCell");
hint.setObjectForKey("CPMEpitheliaSwitchKey", SYData(m));

% PCP switch.
m = zeros(6,3);
%  non-PCP cell.
m(1,1) = celltypeArray.indexOfObject("media");
m(2,1) = celltypeArray.indexOfObject("constrictingCell");
m(3,1) = celltypeArray.indexOfObject("softCell");
m(4,1) = celltypeArray.indexOfObject("stiffCell");
m(5,1) = celltypeArray.indexOfObject("cableCell");
m(6,1) = celltypeArray.indexOfObject("staticCell");
PCPSwitch = SYData(m);
hint.setObjectForKey("CPMPCPSwitchKey", SYData(m));
```

Po and Ao projectors are set for each cell type and define how to calculate reference values for the surface elasticity and the area constraint.
It calculates the reference value by a polynomial with respect to an initial value of the perimeter or the area. In a row vector, i-th element represents a coefficient for a term of i - 1 degree. To set the reference value equal to the initial value, the vector shall be (0, 1).
Here the apical surface reference value of the constricting cell is set 0.
```
%% Po and Ao projectors
array = SYArray;
% [total, apical, lateral, basal, upstream, downstream].
% media
f = [0; 0; 0; 0; 0; 0];
data = SYData(f);
array.addObject(data);
% constrictingCell
f = [-43,1; 0,0; 0,1; 0,1; 0,1; 0,1];
data = SYData(f);
array.addObject(data);
% softCell
f = [0,1; 0,1; 0,1; 0,1; 0,1; 0,1];
data = SYData(f);
array.addObject(data);
% stiffCell
f = [0,1; 0,1; 0,1; 0,1; 0,1; 0,1];
data = SYData(f);
array.addObject(data);
% cableCell
f = [0,1; 0,1; 0,1; 0,1; 0,1; 0,1];
data = SYData(f);
array.addObject(data);
% staticCell
f = [0,1; 0,1; 0,1; 0,1; 0,1; 0,1];
data = SYData(f);
array.addObject(data);
hint.setObjectForKey("CPMPoProjectorArrayKey", array);

array = SYArray;
% [total area, bodyFluid, extraBodyFluid, apical_ECM, cytosol_const, 
%   cytosol_soft, cytosol_stiff, cytosol_static].
% media
f = [0; 0; 0; 0; 0; 0; 0; 0];
data = SYData(f);
array.addObject(data);
% constrictingCell
f = [0,1; 0,1; 0,1; 0,1; 0,1; 0,1; 0,1; 0,1];
data = SYData(f);
array.addObject(data);
% softCell
f = [0,1; 0,1; 0,1; 0,1; 0,1; 0,1; 0,1; 0,1];
data = SYData(f);
array.addObject(data);
% stiffCell
f = [0,1; 0,1; 0,1; 0,1; 0,1; 0,1; 0,1; 0,1];
data = SYData(f);
array.addObject(data);
% cableCell
f = [0,1; 0,1; 0,1; 0,1; 0,1; 0,1; 0,1; 0,1];
data = SYData(f);
array.addObject(data);
% staticCell
f = [0; 0; 0; 0; 0; 0; 0; 0];
data = SYData(f);
array.addObject(data);
hint.setObjectForKey("CPMAoProjectorArrayKey", array);
```

Potential components defines which parts to be assigned with a potential energy for each cell type.
The potential energy is represented by a scalar field, and a simulation can be set with multiple fields.
The potential energy can be assigned to both the materials and the marked surface parts. A vector v = (m, i) indicates that a material m (value in subcellIndices) is assigned a potential energy represented by i-th field. A vector w = (n, j) indicates that n-th marked part is assigned j-th field.
For parts labeling, 1: cell total body, 2: cell total perimeter, 3: apical surface, 4: lateral cell-cell junction, 5: basal surface, 6: adherens junction, 7: PCP upstream cell-cell junction, 8: PCP downstream cell-cell junction, 9: PCP upstream adherens junction, 10: PCP downstream adherens junction.
For a cell type, the potential energies are set by an array containing a first matrix whose row vectors are v and a second matrix whose row vectors are w.
```
%% potential components
array = SYArray;
% media.
array.addObject(SYArray);
% constrictingCell.
array.addObject(SYArray);
% softCell.
array.addObject(SYArray);
% stiffCell.
array.addObject(SYArray);
% cableCell.
sc = SYData([]);
sl = SYData([6,1]);
array.addObject(SYArray(sc,sl));
% staticCell.
array.addObject(SYArray);
hint.setObjectForKey("CPMPotentialComponentArrayKey", array);
```

It also needs several parameters including fluctuation allowance, color table, and duration.
```
%% fluctuation allowance.
hint.setObjectForKey("CPMTKey", 240);


%% drawing color
lut = [255,255,255; ... body fluid
    220,220,220; ... extra-body fluid
    220,220,0; ... apical_ECM.
    245,245,255; ... cytosol constricting
    255,245,245; ... cytosol soft
    245,255,245; ... cytosol stiff
    200,200,200; ... cytosol static
    
    255,255,255; ... clear-color.
    
    30,235,126; ... PCP downstream adherens junction.
    235,126,30; ... PCP upstream adherens junction.
    30,200,126; ... PCP downstream rim.
    200,126,30; ... PCP upstream rim.
    222,222,12; ... adherens junction.
    126,126,245; ... apical rim.
    245,126,126; ... basal rim.
    126,245,126; ... lateral rim.
    126,126,126 ... media rim.
    ];
hint.setObjectForKey("CPMLutKey",SYData(lut));


%% stack setting
hint.setObjectForKey("CPMFrameIntervalKey", 10000);
hint.setObjectForKey("CPMStackLimitKey", 200);

%% exporting hint
data = hint.data;
data.writeToFile("hint.mat");
txt = hint.description;
fid = fopen("hint.txt", "w");
fprintf(fid, txt);
fclose(fid);
```

### Adjusting parameters
The parameters can be easily edited.
In below the contact energy was edited according to an initial cell size and compression ratio.
```
%% map parameters
w = 13;
h = 18;
compression = 1.01;
lambda = 1;
bl = 2;
g = 1000;

%% adjusting Ao projector and contact energy
array = hint.objectForKey("CPMAoProjectorArrayKey");
cell_types = ["constrictingCell","softCell","stiffCell","cableCell"];
cell_soft = "softCell";
for cell_type = cell_types
    index = celltypeArray.indexOfObject(cell_type);
    proj = array.objectAtIndex(index);
    proj.var(:,2) = compression;
end

J_l = lambda * w * h * (compression - 1) * w;
data = hint.objectForKey("CPMContactEnergyKey");
cytosol_array = ["cytosol_const","cytosol_soft","cytosol_stiff", ...
    "cytosol_static", "PCP_upstream_boundary"];
for cytosol = cytosol_array
    index = labelArray.indexOfObject(cytosol);
    for dytosol = cytosol_array
        jndex = labelArray.indexOfObject(dytosol);
        data.var(index,jndex,2) = J_l;
        data.var(jndex,index,2) = J_l;
    end
end

J_a = J_l * bl;
J_b = J_l * bl;
apical_EC_array = ["extraBodyFluid","apical_ECM"];
basal_EC_array = ["bodyFluid"];
cytosol_const = "cytosol_const";
for cytosol = cytosol_array
    index = labelArray.indexOfObject(cytosol);
    for apical_EC = apical_EC_array
        jndex = labelArray.indexOfObject(apical_EC);
        data.var(index,jndex,2) = J_a;
        data.var(jndex,index,2) = J_a;
    end
    for basal_EC = basal_EC_array
        jndex = labelArray.indexOfObject(basal_EC);
        data.var(index,jndex,2) = J_b;
        data.var(jndex,index,2) = J_b;
    end
end
```

### Initial configuration
The initial configuration is prepared by drawing cells one by one into a map.
```
%% drawing map
frame = [50,w * 16];
map = CPMPotentialMap;
bitmap = zeros(frame,'uint16'); % cell map.
citmap = zeros(frame,'uint16'); % subcellular map.
ditmap = zeros(frame,'uint16'); % cell type map.
eitmap = zeros(frame,'double'); % potential map.
% clear map with media.
bitmap(:) = 1;
citmap(:) = subcellIndices.var(labelArray.indexOfObject("bodyFluid"));
ditmap(:) = celltypeArray.indexOfObject("media");
% draw cuticle.
ecm = zeros(frame,'logical');
ecm(1:4,:) = true;
citmap(ecm) = subcellIndices.var(labelArray.indexOfObject("apical_ECM"));
% draw potential.
half_w = w * 8;
eitmap(:,1:half_w) = repmat(half_w:-1:1,frame(1),1);
eitmap(:,half_w + 1:end) = repmat(1:half_w,frame(1),1);
eitmap = eitmap .* g;
%
map.initWithMaps(bitmap,citmap,ditmap,eitmap,hint);

% draw cells.
% apical-constricting cells.
siz = [h,w];
bitmap = zeros(siz,'uint16');
citmap = zeros(siz,'uint16');
ditmap = zeros(siz,'uint16');
% subcellular components.
% % cytosol.
mask = zeros(siz,'logical');
mask(:) = true;
citmap(mask) = ...
    subcellIndices.var(labelArray.indexOfObject("cytosol_const"));
% cell mask.
mask = citmap > 0;
bitmap(mask) = 1;
% cell type.
ditmap(mask) = celltypeArray.indexOfObject("constrictingCell");
% stack.
stack = cat(3,bitmap,citmap,ditmap);
% 
for i = 5:9
    map.drawCellAtPoint(stack,[5,w * i + 7]);
end

% cable cells.
siz = [h,w];
bitmap = zeros(siz,'uint16');
citmap = zeros(siz,'uint16');
ditmap = zeros(siz,'uint16');
% subcellular components.
% % cytosol.
mask = zeros(siz,'logical');
mask(:) = true;
citmap(mask) = ...
    subcellIndices.var(labelArray.indexOfObject("cytosol_stiff"));
% cell mask.
mask = citmap > 0;
bitmap(mask) = 1;
% cell type.
ditmap(mask) = celltypeArray.indexOfObject("cableCell");
% stack.
stack = cat(3,bitmap,citmap,ditmap);
% 
for i = [4,10]
    map.drawCellAtPoint(stack,[5,w * i + 7]);
end

% stiff cells.
siz = [h,w];
bitmap = zeros(siz,'uint16');
citmap = zeros(siz,'uint16');
ditmap = zeros(siz,'uint16');
% subcellular components.
% % cytosol.
mask = zeros(siz,'logical');
mask(:) = true;
citmap(mask) = ...
    subcellIndices.var(labelArray.indexOfObject("cytosol_stiff"));
% cell mask.
mask = citmap > 0;
bitmap(mask) = 1;
% cell type.
ditmap(mask) = celltypeArray.indexOfObject("stiffCell");
% stack.
stack = cat(3,bitmap,citmap,ditmap);
% 
for i = [1:3,11:13]
    map.drawCellAtPoint(stack,[5,w * i + 7]);
end

% soft cells.
siz = [h,w];
bitmap = zeros(siz,'uint16');
citmap = zeros(siz,'uint16');
ditmap = zeros(siz,'uint16');
% subcellular components.
% % cytosol.
mask = zeros(siz,'logical');
mask(:) = true;
citmap(mask) = ...
    subcellIndices.var(labelArray.indexOfObject("cytosol_soft"));
% cell mask.
mask = citmap > 0;
bitmap(mask) = 1;
% cell type.
ditmap(mask) = celltypeArray.indexOfObject("softCell");
% stack.
stack = cat(3,bitmap,citmap,ditmap);
% 
for i = [0,14]
    map.drawCellAtPoint(stack,[5,w * i + 7]);
end

% boundary cell.
siz = [h,w];
bitmap = zeros(siz,'uint16');
citmap = zeros(siz,'uint16');
ditmap = zeros(siz,'uint16');
% subcellular components.
% % cytosol.
mask = zeros(siz,'logical');
mask(:) = true;
citmap(mask) = ...
    subcellIndices.var(labelArray.indexOfObject("cytosol_static"));
% cell mask.
mask = citmap > 0;
bitmap(mask) = 1;
% cell type.
ditmap(mask) = celltypeArray.indexOfObject("staticCell");
% stack.
stack = cat(3,bitmap,citmap,ditmap);
% 
for i = 15
    map.drawCellAtPoint(stack,[5,w * i + 7]);
end

% data = map.data;
% data.writeToFile("map.mat");
```

### Running simulation
The cells and tissue deformation is simulated with the parameters and initial configuration.
To run multiple trials, the hint and map are first converted to data, and re-instantiated for each trial.
```
%% running simulation
hintData = hint.data;
mapData = map.data;
N = 3;
for n = 1:N
    str = "c01_bl2_g1000_n" + string(n);

    hint = SYDictionary;
    hint.initWithData(hintData);
    map = CPMPotentialMap;
    map.initWithData(mapData);

    hint.setObjectForKey("CPMOutputNameKey", str);

    simulator = CPMPotentialSimulator;
    simulator.initWithMap(map,hint);
    simulator.runSimulation;
end
```

### Invagination evaluation
Distance between the cells and the apical ECM is measured and exported.
```
%% measuring depth
ECM_label = "apical_ECM";
celltype_cc = "constrictingCell";
celltype_sc = "stiffCell";

fileName = "depth_cell.2.txt";

txt = "";
for n = 1:N
    str = "c01_bl2_g1000_n" + string(n);
    txt = txt + str;
    
    data = SYData;
    data.initWithContentsOfFile(str + ".end.mat");
    simulator = CPMPotentialSimulator;
    simulator.initWithData(data);

    txt = txt + "\tconstricting cells";
    d = ss_measure_distance_to_cell(simulator,ECM_label,celltype_cc,false);
    for i = 1:length(d)
        txt = txt + "\t" + string(d(i));
    end
    
    txt = txt + "\tsurrounding cells";
    d = ss_measure_distance_to_cell(simulator,ECM_label,celltype_sc,false);
    for i = 1:length(d)
        txt = txt + "\t" + string(d(i));
    end

    txt = txt + "\n";
end

fid = fopen(fileName, "w");
fprintf(fid, txt);
fclose(fid);
```
