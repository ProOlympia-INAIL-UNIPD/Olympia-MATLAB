function [PointsLocal, PointsLocalMean] = points2local(points, gndT_loc)
% POINTS2LOCAL is a function that given N points stacked in the matrix
% POINT (first input) and the pose of the local coordinate system with
% respect to a ground, gives back the points in the local frame averaged
% over time.
% 
% INPUTS:  - POINTS [nF, 3, nP]
%          - gndT_loc [4, 4, nF]
%
% OUTPUT:  - POINTSLOCAL [3 x nP]

points = permute(points, [2 3 1]);
points(4,:,:) = 1;

PointsLocal = nan(size(points));

locT_gnd = Tinv(gndT_loc);
for iF = 1:size(points,3)
    PointsLocal(:,:,iF) = locT_gnd(:,:,iF)*points(:,:,iF);
end

PointsLocal(4,:, :) = [];
PointsLocalMean = nanmean(PointsLocal, 3);
PointsLocal = permute(PointsLocal, [3,1,2]);


end