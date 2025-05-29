function stiffness=lumpedStiffness(PP,COP,GRF,FS,FO,mass)
%function stiffness=lumpedStiffness(PP,COP,GRF,FS,FO,mass)
%
% computes leg stiffness and deltaL with lumped approach using a selected
% Proximal Point (PP) and the Center of Pressure (COP) of the Ground
% Reaction Force (GRF)
% for each step between FootStrike (FS) and FootOff(FO) events
% output is a structure containing the segmented track, the normalized
% cycles, the average cycle and its standard deviation
% 
% all data must be given in X (medio-lateral), Y (antero-posterior), Z
% (vertical) coordinates.
%
% author: Giuseppe Zullo
% date: 24/01/2024
% version: 1.0

for i=length(FS):-1:1
    while isnan(COP(FS(i),1)) %solve potential NaN issue with COP and GRF
        FS(i)=FS(i)+1;
    end
    while isnan(COP(FO(i),1))
        FO(i)=FO(i)-1;
    end

L{i}=PP(FS(i):FO(i),:)- COP(FS(i):FO(i),:); %vector PP-COP
Lnorm{i}=vecnorm(L{i},2,2);                 % ||PP-COP||
Lvers{i}=L{i}./Lnorm{i};                    %versor PP-COP
GRF_stance{i}=GRF(FS(i):FO(i),:);           %GRF track during stance
end
%% Parameters Calculation
for i=1:length(FS) %repeat for each step
    s.L{i}=Lnorm{i};                            % L during stance
    s.L0(i)=Lnorm{i}(1);                        % L0 at FS
    s.Lfin(i)=Lnorm{i}(end);                    % L at FO
    s.Lmin(i)=min(Lnorm{i});                    % minimum L during step
    s.deltaAPprox(i)=PP(FO(i),2)-PP(FS(i),2);   % forward displacement of proximal point relative to initial position
    s.deltaL{i}=Lnorm{i}(1)-Lnorm{i};           % dL=L0-L
    s.V0(i)=PP(FS(i),3);                        % vertical position of PP at FS
    s.Vfin(i)=PP(FO(i),3);                      % vertical position of PP at FS
    s.Vmin(i)=min(PP(FS(i):FO(i),3));           % minimum vertical position of PP during stance 
    s.deltaVfin(i)=PP(FS(i),3)-PP(FO(i),3);     % vertical trajectory of PP during stance
    s.deltaV{i}=PP(FS(i):FO(i),3)-PP(FS(i),3);  % vertical trajectory of PP relative to initial position
    s.modGRF{i}=dot(Lvers{i},GRF_stance{i},2);  % magnitude of GRF projection on leg versor
    s.GRFcomp{i}=s.modGRF{i}.*Lvers{i};         % components GRF projection on leg versor
    [GRFmax, Imax]=max(s.modGRF{i});
    s.GRFmax(i)=GRFmax;                         % maximum GRF during stance
    s.K(i)=GRFmax/s.deltaL{i}(Imax);            % leg stiffness at GRFmax [N/mm]
    s.Knorm(i)=s.K(i)/(9.81*mass);              % leg stiffness at GRFmax [(N/BW)/mm]
    s.Kist{i}=(s.modGRF{i})./s.deltaL{i};       % istantaneous leg stiffness [N/mm]
    s.Kist{1,i}(s.Kist{1,i}==Inf)=nan;

    s.PP{i}=PP(FS(i):FO(i),:);
end
%% arrange data for DataBase
fns=fieldnames(s);
for i=1:length(fns)
    if isnumeric(s.(fns{i})) %method for scalar values
    stiffness.(fns{i}).all=s.(fns{i});
    stiffness.(fns{i}).mean=mean(s.(fns{i}));
    stiffness.(fns{i}).std=std(s.(fns{i}));
    elseif iscell(s.(fns{i})) %method for Nxm  arrays (N: number of samples)
    stiffness.(fns{i}).all=s.(fns{i});
    for j=length(s.(fns{i})):-1:1 % normalize j-th step
    tmp(:,:,j)=time2cycle([],s.(fns{i}){j},101);
    end
    stiffness.(fns{i}).cycle=squeeze(tmp); % removes third dimension if m=1
    stiffness.(fns{i}).mean=mean(tmp,3);   % average steps
    stiffness.(fns{i}).std=std(tmp,1,3);   % calculate std
    clear tmp
    end
end