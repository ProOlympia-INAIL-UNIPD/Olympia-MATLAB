function [out]=getParamGRFv2(forceplat,events,mass)
%getParamGRF has in input GRF in the following order X (medial lateral
%axis), Y (anterior posterior axis) and Z (vertical axis). input GRF are NOT
%normalized to body weight. FS-foot strike, FO-foot off, strike and frequency (force plate acquisition rate).
%NO CHANGE in the acquisizion rate of kinematics is either performed.
%Output is a struct containing impulses normalized along vertical and
%anterior-posterior axes and GRF_param struct with GRFs parameters (min,
%max, mean of all steps of each limb) and curves of GRFs normalized on stance time.
%AUTHOR Breban Samira Giuliana
%25-05-2023
%edit on 22/01/2024 by Giuseppe Zullo
%added additional (optional, required for kinside) input GRM, COP. Formatted accordingly to GRF.
%added additional output kinside which returns a structure to be appended
%in the DB in the form ....kinetics.filtered.(side)=kinside 
%edit on 17/04/2024 by Giuseppe Zullo
%replaced GeneralInfo with mass
%removed strike
arguments
    forceplat forcePlatformType2
    events Event
    mass double {mustBePositive(mass)}
end
% components
V=3;
AP=2;
ML=1;

GRF=forceplat.GRF/(9.81*mass);
GRM=forceplat.GRM/(9.81*mass);
COP=forceplat.COP/(9.81*mass);
f_force=forceplat.SampleRate;
time=(0:length(GRF(:,3))-1)/f_force;
footContact=0*time;
events=events.exportEvents('analog',false);
for s=["Left","Right"]
    try
    FS=events.(s).Foot_Strike;
    FO=events.(s).Foot_Off;
for j=length(FS):-1:1
    fc=FS(j);
    fo=FO(j);
    footContact(fc:fo)=1;
    contactTime(j)=(fo-fc)/f_force;
    %% GRFs normalized along stance time


    %needed for DB
    GRF_stance(:,:,j)=time2cycle([],GRF(fc:fo,:),101);
    GRM_stance(:,:,j)=time2cycle([],GRM(fc:fo,:),101);
    COP_stance(:,:,j)=time2cycle([],COP(fc:fo,:),101);

    %needed for impulse and GRF_param
    F_AP_stance(:,j)=GRF_stance(:,AP,j);
    F_V_stance(:,j)=GRF_stance(:,V,j);
    F_ML_stance(:,j)=GRF_stance(:,ML,j);
    %% GRFs mean values at each step of the same limb
    GRF_AP_mean(j)=mean(GRF(fc:fo,AP));
    GRF_V_mean(j)=mean(GRF(fc:fo,3));

    %% GRFs max values at each step of the same limb
    GRF_AP_max(j)=max(GRF(fc:fo,AP));
    GRF_V_max(j)=max(GRF(fc:fo,V));

    %% GRFs min values at each step of the same limb
    GRF_AP_min(j)=min(GRF(fc:fo,AP));
    GRF_V_min(j)=min(GRF(fc:fo,V));

    %% Impulses Y at each step of the same limb
    I_V(j)=trapz(time(fc:fo),GRF(fc:fo,V)-1); %(Fy-BW)/BW= Fy/BW-1
    %% Braking/Propulsive impulses X for each step of the same limb
    GRF_AP_e=GRF(fc:fo,AP);
    zerocross=[ne(diff(sign(GRF_AP_e)),0);false]; %find zero crossing
    Ic=cumtrapz(GRF_AP_e)/f_force;
    Ic=[Ic(zerocross)', Ic(end)];
    In=Ic-[0,Ic(1:end-1)];
    I_APbrak(j)=sum(In(In<0));
    I_APprop(j)=max(In);
    I_AP(j)=sum(In);
    I_APsp(j)=I_AP(j)-(I_APbrak(j)+I_APprop(j)); %initial propulsion
end
%% MEAN IMPULSES (net, braking, propulsion) along Y and X direction for all steps of the same limb
impulses.AP_mean=mean(I_AP);
impulses.AP_step=I_AP;
impulses.AP_std=std(I_AP);
impulses.V_mean=mean(I_V);
impulses.V_step=I_V;
impulses.V_std=std(I_V);

impulses.APprop_mean=mean(I_APprop);
impulses.APprop_step=I_APprop;
impulses.APprop_std=std(I_APprop);
impulses.APbrak_mean=mean(I_APbrak);
impulses.APbrak_step=I_APbrak;
impulses.APbrak_std=std(I_APbrak);
impulses.APsp_mean=mean(I_APsp);
impulses.APsp_step=I_APsp;
impulses.APsp_std=std(I_APsp);
%% IMPULSES along the last 20% of the stance time
%{
FC_80=round(FS+(FO-FS)*0.80);
for j=1:length(FS)
    I_APprop80(j)=trapz(time(FC_80(j):FO(j)),GRF_comb(FC_80(j):FO(j),2));
    I_V80(j)=trapz(time(FC_80(j):FO(j)),GRF_comb(FC_80(j):FO(j),3));
end
impulses.APprop80_mean=mean(I_APprop80);
impulses.APprop80_step=I_APprop80;
impulses.APprop80_std=std(I_APprop80);
impulses.V80_mean=mean(I_V80);
impulses.V80_step=I_V80;
impulses.V80_std=std(I_V80);
%}
%% MEAN, MAX, MIN GRFs FOR ALL STEPS of the same limb
GRF_param.AP_mean=mean(GRF_AP_mean);
GRF_param.AP_std=std(GRF_AP_mean);
GRF_param.V_mean=mean(GRF_V_mean);
GRF_param.V_std=std(GRF_V_mean);

GRF_param.AP_max_mean=mean(GRF_AP_max);
GRF_param.AP_max_step=GRF_AP_max;
GRF_param.AP_max_std=std(GRF_AP_max);
GRF_param.V_max_mean=mean(GRF_V_max);
GRF_param.V_max_step=GRF_V_max;
GRF_param.V_max_std=std(GRF_V_max);

GRF_param.AP_min_mean=mean(GRF_AP_min);
GRF_param.AP_min_step=GRF_AP_min;
GRF_param.AP_min_std=std(GRF_AP_min);
GRF_param.V_min_mean=mean(GRF_V_min);
GRF_param.V_min_step=GRF_V_min;
GRF_param.V_min_std=std(GRF_V_min);

F_AP_MEAN=mean(F_AP_stance,2);
F_V_MEAN=mean(F_V_stance,2);
F_ML_MEAN=mean(F_ML_stance,2);
%% Module of GRF where Fx is maximum (propulsion side) of all steps of the same limb
GRF_abs=sqrt(F_AP_MEAN.^2+F_V_MEAN.^2+F_ML_MEAN.^2);
[~, pos]=(max(F_AP_MEAN));
GRF_param.GRF_abs_max=GRF_abs(pos);

%% GRF_stance
GRF_param.AP_stance=F_AP_stance;
GRF_param.V_stance=F_V_stance;
GRF_param.ML_stance=F_ML_stance;

%% make DB entry
GRFdb=GRF;
GRMdb=GRM;
COPdb=COP;
GRFdb(~footContact,:)=nan;
GRMdb(~footContact,:)=nan;
COPdb(~footContact,:)=nan;
ax={'x','y','z'};
out.(s).GRF.track=GRFdb;
out.(s).GRF.cycle=GRF_stance;
out.(s).GRF.mean=mean(GRF_stance,3);
out.(s).GRF.std=std(GRF_stance,1,3);

out.(s).GRM.track=GRMdb;
out.(s).GRM.cycle=GRM_stance;
out.(s).GRM.mean=mean(GRM_stance,3);
out.(s).GRM.std=std(GRM_stance,1,3);

out.(s).COP.track=COPdb;
out.(s).COP.cycle=COP_stance;
out.(s).COP.mean=mean(COP_stance,3);
out.(s).COP.std=std(COP_stance,1,3);


impdirDBID={'vertical','horizontal','braking','propulsive','startprop'};
impdirID=  {'V','AP','APbrak','APprop','APsp'};

for j=1:length(impdirID)
out.(s).impulse.(impdirDBID{j}).cycle=impulses.([impdirID{j} '_step']);
out.(s).impulse.(impdirDBID{j}).mean=impulses.([impdirID{j} '_mean']);
out.(s).impulse.(impdirDBID{j}).std=impulses.([impdirID{j} '_std']);
end

paramdirDBID={'vertical','horizontal','mediolateral'};
idx=[V, AP, ML];
for j=1:length(paramdirDBID)
out.(s).parameters.(paramdirDBID{j}).max.cycle=squeeze(max(GRF_stance(:,idx(j),:)))';
out.(s).parameters.(paramdirDBID{j}).min.cycle=squeeze(min(GRF_stance(:,idx(j),:)))';
out.(s).parameters.(paramdirDBID{j}).mean.cycle=squeeze(mean(GRF_stance(:,idx(j),:)))';
for k={'max','min','mean'}
out.(s).parameters.(paramdirDBID{j}).(k{:}).mean=mean(out.(s).parameters.(paramdirDBID{j}).(k{:}).cycle);
out.(s).parameters.(paramdirDBID{j}).(k{:}).std=std(out.(s).parameters.(paramdirDBID{j}).(k{:}).cycle);
end
end
out.(s).ContactTime=contactTime;
clear contactTime GRF_AP_mean GRF_AP_min GRF_AP_max GRF_V_max GRF_V_min GRF_V_mean
clear I_V I_AP I_APbrak I_APprop I_APsp 
clear F_ML_stance F_V_stance F_AP_stance GRF_stance GRM_stance COP_stance
    catch ME
        warning(ME.message)
        out.(s)=[];
    end
end