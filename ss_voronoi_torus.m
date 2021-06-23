% Short script to tessellate torus plane with given seeds.
% Written by Satoshi Yamashita.

function result = ss_voronoi_torus(seeds,frame)
% Function to get an image of Voronoi tessellation.
% result = ss_voronoi_torus(seeds,frame)
% Argument seeds is a matrix where its row represents a seed, its first and
%   second column represent x and y coordinates of the seed.
% Argument frame specify a size of the plane (rows, columns).
% Return value is a bitmap tesselated by cells with their label.

% parameters.
iterM = 50;

% initialize.
bitmap = zeros(frame(1), frame(2));
hint = SYDictionary;
map = CPMMap(bitmap,bitmap,bitmap,hint);

% mark seeds.
for i = 1:size(seeds,1)
    p = seeds(i,:);
    bitmap(p(2),p(1)) = i;
end

% expand cells from seeds.
c = size(seeds,1);
c_m = frame(1) * frame(2);
for r = 1:iterM
    % expand cells.
    mask = lf_circle(r);
    for i = 1:size(seeds,1)
        p = seeds(i,:);
        indices = p(2) - r:p(2) + r;
        jndices = p(1) - r:p(1) + r;
        kndices = map.indices(indices,jndices);
        scope = bitmap(kndices);
        nask = mask & (scope == 0);
        scope(nask) = i;
        bitmap(kndices) = scope;
        
        c = c + sum(nask(:));
        if c >= c_m
            break;
        end
    end
end

result = bitmap;
end

function result = lf_circle(r)
a = (-r:r).^2;
m = sqrt(a + a');
result = round(m) <= r;
end
