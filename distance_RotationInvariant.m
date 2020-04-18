function [dist,shift] = distance_RotationInvariant(WECT1,WECT2)

% Input: two WECTs or SWECTs of the same size. This assumes that the
% directions used to compute each are the same.

% Output: check every cyclic permutation of the second WECT. dist is the
% L2 distance between WECT1 and the permuted WECT2. shift is the best
% shift.

num_directions = size(WECT1,2);

distances = zeros(1,num_directions);

for d = 1:num_directions
    WECT2_shifted = circshift(WECT2,-d,2);
    distances(d) = norm(WECT1 - WECT2_shifted,'fro');
end

[dist,shift] = min(distances);