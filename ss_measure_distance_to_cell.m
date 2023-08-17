

function result = ss_measure_distance_to_cell(simulator,ECM_label, ...
    celltype,maskDelaminated)

if isa(celltype,"string")
    celltypeArray = simulator.hint.objectForKey("CPMCelltypeArrayKey");
    celltype = celltypeArray.indexOfObject(celltype);
end

map = simulator.map;
ditmap = map.distanceMapFromComponent(ECM_label,true);

cellArray = simulator.cellArray;
d = nan([1,cellArray.count]);
m = zeros([1,cellArray.count],'logical');

for i = 1:cellArray.count
    cel = cellArray.objectAtIndex(i);
    if cel.cellType ~= celltype
        continue
    end

    m(i) = true;

    if maskDelaminated && (isempty(cel.apicalRim) || ...
            (length(cel.apicalRim) == 1 && isnan(cel.apicalRim)))
        continue
    end

    indices = cel.frame(1) + 1:cel.frame(2) - 1;
    jndices = cel.frame(3) + 1:cel.frame(4) - 1;
    kndices = map.indices(indices,jndices);
    scope = ditmap(kndices);
    mask = map.cellMap(kndices) == cel.label;
    d(i) = min(scope(mask));
end

result = d(m);
end
