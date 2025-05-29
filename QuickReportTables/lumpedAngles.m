function angles=lumpedAngles(markers,markerlist,DP,FS,FO)
%function angles=lumpedAngles(markers,markerlist,DP,FS,FO)
%
% computes 1D angles with lumped approach between points in markers defined by markerlist
% and distal point DP given as input for each step between FootStrike (FS)
% and FootOff(FO) events
% output is a structure containing the segmented track, the normalized
% cycles, the average cycle and its standard deviation
% 
% all data must be given in X (medio-lateral), Y (antero-posterior), Z
% (vertical) coordinates.
% example of proximal points (markers)
% COM: Center of Mass
% GT: Great Trocanter
% HJC: Hip Joint Center
% KJC: Knee Joint Center/Lateral Femoural Epicondyle
% AJC: Ankle Joint Center/Lateral Malleoulus
% 
% example of distal points
% FH:  Foot Head (e.g., V Metatarsal Head) for lambda
% COP: Center of Pressure for theta
% 
% author: Giuseppe Zullo
% date: 24/01/2024
% version: 1.1
% supports single point (Nx3) as input, markerlist assign name to that
% point
markerlist={markerlist};
if isnumeric(markers)
markers=struct(markerlist{1},markers);
end

%% angle computation
newlist=markerlist; %markerlist contains side (L/R) which is unnecessary
for j=1:length(markerlist)
        if ismember(markerlist{j}(1),{'L','R'}); newlist{j}(1)=[];end %remove side from markername
end
for i=1:length(FS) %repeat for each step
    for j=1:length(newlist)
v=markers.(markerlist{j})(FS(i):FO(i),:)-DP(FS(i):FO(i),:); %create vector
a.(newlist{j}){i}=atan2d(v(:,2),v(:,3)); %calculate 1D angle
    end
end

%% rearrange datastruct
fns=fieldnames(a);
for i=1:length(fns)
    if isnumeric(a.(fns{i})) %method for scalar values
    angles.(fns{i}).all=a.(fns{i});
    angles.(fns{i}).mean=mean(a.(fns{i}));
    angles.(fns{i}).std=std(a.(fns{i}));
    elseif iscell(a.(fns{i})) %method for Nxm  arrays (N: number of samples)
    angles.(fns{i}).all=a.(fns{i});
    for j=length(a.(fns{i})):-1:1 % normalize j-th step
    tmp(:,:,j)=time2cycle([],a.(fns{i}){j},101);
    end
    angles.(fns{i}).cycle=squeeze(tmp); % removes third dimension if m=1
    angles.(fns{i}).mean=mean(tmp,3);   % average steps
    angles.(fns{i}).std=std(tmp,1,3);   % calculate std
    clear tmp
    end
end


end
   