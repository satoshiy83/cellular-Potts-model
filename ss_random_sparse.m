% Short script to spread points randomly and sparsely.
% Written by Satoshi Yamashita.

function result = ss_random_sparse(frame,r_min,c_max)
% Function to get an array of points randmly and sparsely distributed.
% result = ss_random_sparse(frame,r_min,c_max)
% Argument frame specifies a space in which points will be spread.
% Argument r_min specifies a minimum distance between points.
% Argument c_max specifies a maximum number of points to be spread.
% Return value is a {n,2} array of points, where the first column
%   represents x and the second column represents y.

% parameters.
n_add = 20;
iterM = 20000;

% initialize.
bitmap = zeros(frame(1),frame(2),'uint16');
hint = SYDictionary;
map = CPMMap(bitmap,bitmap,bitmap,hint);

mask_m = lf_circle(r_min);
bitmap = ones(frame(1),frame(2),'logical');

seeds = [];

% spread points.
count = 0;
for i = 1:iterM
    % choose candidates.
    array = find(bitmap);
    if length(array) < 1
        break;
    elseif length(array) < n_add
        n = length(array);
        n_add = round(n_add / 2);
    else
        n = n_add;
    end
    [r,c] = ind2sub(frame,array(randperm(length(array),n)));
    
%     for j = 1:length(r)
%         if ~bitmap(r,c)
%             disp('wrong candidate!');
%         end
%     end
    
    dx = c - c';
    dy = r - r';
    d = sqrt(dx .^ 2 + dy .^ 2);
    array = find(~any(triu(d < r_min,1),2));
    if isempty(array)
        continue;
    end
    
    if count + length(array) > c_max
        array = array(1:c_max - count);
    end
    
    % mask around.
    seeds = cat(1,seeds,[c(array),r(array)]);
    
%     dx = seeds(:,1) - seeds(:,1)';
%     dy = seeds(:,2) - seeds(:,2)';
%     d = sqrt(dx .^ 2 + dy .^ 2);
%     if any(any(triu(d < r_min,1),2))
%         disp('close points!');
%     end
    
    for j = array(:)'
        indices = r(j) - r_min:r(j) + r_min;
        jndices = c(j) - r_min:c(j) + r_min;
        kndices = map.indices(indices,jndices);
        scope = bitmap(kndices);
        scope(mask_m) = false;
        bitmap(kndices) = scope;
    end
    
    count = count + length(array);
    if count >= c_max
        break;
    end
end

result = seeds;
end

function result = lf_circle(r)
a = (-r:r).^2;
m = sqrt(a + a');
result = round(m) <= r;
end
