% Short script to export Ao and Po of cells.
% Written by Satoshi Yamashita.

function [result1,result2] = ss_export_AoPo(cells,name)
% Function to save reference volume and surface are of cells.
% result = ss_export_AoPo(cells,name)
% Argument cells is an SYArray instance carrying cells.
% Optional argument name is a path to which a result is saved.
% Return value is an SYData instance of labels and an array of the
%   reference values if the number of output is 1.
% Return values are labels and an array of the reference values if the
%   number of outputs is 2.

    labels = zeros(cells.count,1);
    array = SYArray;
    for i = 1:cells.count
        c = cells.objectAtIndex(i);
        labels(i) = c.label;
        array.addObject(SYArray(SYData(c.Ao),SYData(c.Po)));
    end
    s.labels = labels;
    s.AoPo = array.data.var;
    
    data = SYData(s);
%     result = data;
    
    if nargin > 1 && ~isnan(name)
        data.writeToFile(name);
    end
    
    if nargout > 1
        result1 = labels;
        result2 = array;
    else
        result1 = data;
    end
end
