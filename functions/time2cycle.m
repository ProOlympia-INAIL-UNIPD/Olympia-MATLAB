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
% Version: 2.0.1
% Changelog: 
%            - 2.0 added nan behavior for tail and long gaps
% 12/05/2025 - 2.0.1 buxfixing in splinefix
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

if dim==1 %data is column oriented
    cycle=splinefix(t,y',t_new,NameValue)';
    
else      %data is row oriented or N-D
    cycle=splinefix(t,y,t_new,NameValue);
end
end

function cycle=splinefix(t,y,t_new,NameValue)
nans=any(isnan(y),1); %search nans
discard=false(1,length(t_new)); %create array to mask unwanted values
i_max=floor(length(t)*NameValue.MaxGap); %maximum number of points of gap
if all(nans)||sum(nans)==0
else
nans=[false nans false]; %to consider also tails in searching for gaps
st_gap=find(diff(nans)==1); %begin of a gap
end_gap=find(diff(nans)==-1)-1; %end of a gap

   if NameValue.AllowTailPrediction==false %if no tailprediction
      if nans(2) %if tip is missing
         discard(t_new<t(end_gap(1)))=true; %mask the output
         st_gap(1)=[];                      %remove this gap
         end_gap(1)=[];
      end
      if nans(end-1)
         discard(t_new>t(st_gap(end)))=true;
         st_gap(end)=[];
         end_gap(end)=[];
      end
   end

   for i=1:length(st_gap)
        if end_gap(i)-st_gap(i)>i_max %check that the i-th gap is shorter than allowed or mask output
            discard((t(st_gap(i))<t_new)&(t_new<t(end_gap(i))))=true;
        end
   end
end
warning off
try
   cycle=spline(t,y,t_new);
catch
   cycle=y(:,1:length(t_new));
   warning on
   return
end
warning on
cycle(:,discard)=nan;
end
