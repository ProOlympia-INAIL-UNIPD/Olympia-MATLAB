function cycle=nanspline(t,y,t_new,NameValue)
arguments
    t (1,:) double
    y (:,:) double
    t_new (1,:) double
    NameValue.MaxGap=20/100; % 20% max gap is allowed
    NameValue.AllowTailPrediction=false; %tails are not interpolated
end

nans=any(isnan(y),1); %search nans
discard=false(1,length(t_new)); %create array to mask unwanted values
i_max=floor(length(t)*NameValue.MaxGap); %maximum number of points of gap
if all(nans)
   cycle=nan(size(y,1),length(t_new));
   return
elseif sum(nans)==0
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