function [cycle,t_new]=time2cycle(t,y,npoints,dim,NameValue)
% [cycle,t_new]=time2cycle(t,y,npoints,dim,NameValue)
% Returns the interpolated cycle starting from the data in y. if y has more than 1-D, 
% the function interpolates along the longest dimension of y
% INPUTS:
% -t: time base of original data 
% -y: original data
% -npoints: number of target points
% -dim: dimension over which interpolate data
% -NameValue controls missing value
% OUTPUTS:
% -cycle: interpolated data
% -t_new: new time vector associated with cycle
% Author: Giuseppe Zullo
% Date: 26-06-2024
% Version: 2.0
% Changelog: added nan behavior for tail and long gaps
arguments
    t (:,1) double
    y (:,:,:) double
    npoints (1,1) double
    dim (1,1) double =0
    NameValue.MaxGap=20/100; % 20% max gap is allowed
    NameValue.AllowTailPrediction=false; %tails are removed
end

if dim==0
    [~,dim]=max(size(y));
end


if isempty(t)
    t=1:size(y,dim);
end

t_new=linspace(t(1),t(end),npoints);
nv=namedargs2cell(NameValue);
if dim==1 %data is column oriented
    cycle=nanspline(t,y',t_new,nv{:})';
    
else      %data is row oriented or N-D
    cycle=nanspline(t,y,t_new,nv{:});
end
end
