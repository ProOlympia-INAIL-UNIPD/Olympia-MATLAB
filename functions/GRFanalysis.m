function [out]=GRFanalysis(trial,normalize)
% this function inherits the structure from the getGRFparam and
% getGRFparamv2, adapting inputs to the class Trial and the outputs to the
% desired functionality.
% The function embeds the PCA algorithm to automatically detect the
% AnteroPosterior, MedioLateral
% and Vertical directions, aligning the GRF to the actual running
% directions.
arguments
    trial (1,1) Trial
    normalize string {mustBeMember(normalize,["no","bodymass","bodyweight"])}="bodyweight";
end


forceplat=trial.ForcePlatform;
if not(isscalar(forceplat))
    warning('Number of active forceplates in the trial (%i) is not compatible with GRF analysis',numel(forceplat));
end

events=trial.Events;

try %import mass from trial and adjust units
mass=trial.Metadata.PROCESSING.mass;
mustBePositive(mass);
switch normalize
    case "no"
        mass=1;
        unit="N%s";
    case "bodymass"
        mass=mass;
        unit="N%s/kg";
    case "bodyweight"
        mass=mass*9.81;
        unit="N%s/N";
end
catch
    mass=1;
    unit="N%s";
end

   
GRF=sum(cat(3,forceplat.GRF),3)/mass;
if isscalar(forceplat)
COP=forceplat.COP;
% align the COP to the principal directions
R=pca(COP);
else
    notok=true;
    while notok
    xyz=inputdlg(["AnteroPosterior:","MedioLateral:","Vertical:"],"Enter trial direction",[1 30],["+Y","+X","+Z"]);
    xyz=char(xyz);
    s=xyz(:,1);
    s=s=='-';
    xyz=double(xyz(:,end)-'X')+1;
    R=eye(3);
    R=R(:,xyz);
    R(:,s)=-R(:,s);
    if det(R)~=1
       ans=questdlg("Invalid set of axes, Try again?");
    switch ans
    case "Yes"
        continue
    case "No"   
        return
    otherwise
        return
    end
    end
    notok=false;
    end
end

    GRF=GRF*R';

AP=1; % now the first direction has the most variability (i.e., is the running direction)
ML=2;
V=3;  % in a typical flat laboratory, the vertical should change least (i.e., is vertical)


f_force=forceplat.SampleRate;

events=exportEvents(events,'analog',false);
events=rmfield(events,'units');

for ctx=string(fieldnames(events))' %run the analysis for each context of events
    try
    FS=events.(ctx).Foot_Strike;
    FO=events.(ctx).Foot_Off;
for j=length(FS):-1:1
    fc=FS(j);
    fo=FO(j);
    
    %% MAX, MIN, MEAN
    % anteroposterior
    out.GRF_PARAM.Units=sprintf(unit,"");
    out.GRF_PARAM.(ctx+"ContactTime")(j)=(fo-fc)/f_force;
    out.GRF_PARAM.(ctx+"AnteroPosteriorMax")(j)=mean(GRF(fc:fo,AP));
    out.GRF_PARAM.(ctx+"AnteroPosteriorMin")(j)=min(GRF(fc:fo,AP));
    out.GRF_PARAM.(ctx+"AnteroPosteriorMaxMean")(j)=max(GRF(fc:fo,AP));
    % vertical
    out.GRF_PARAM.(ctx+"VerticalMax")(j)=max(GRF(fc:fo,V));
    out.GRF_PARAM.(ctx+"VerticalMean")(j)=mean(GRF(fc:fo,V));
    %% IMPULSES
    out.GRF_IMPULSES.Units=sprintf(unit,"s");
    out.GRF_IMPULSES.(ctx+"Vertical")(j)=trapz(GRF(fc:fo,V)-1)/f_force; %(Fy-BW)/BW= Fy/BW-1
    %% Braking/Propulsive impulses X for each step of the same limb
    GRF_AP_e=GRF(fc:fo,AP);
    zerocross=[ne(diff(sign(GRF_AP_e)),0);false]; %find zero crossing
    Ic=cumtrapz(GRF_AP_e)/f_force;
    Ic=[Ic(zerocross)', Ic(end)];
    In=Ic-[0,Ic(1:end-1)];
    
    out.GRF_IMPULSES.(ctx+"AnteroPosteriorStartPropulsion")(j)=sum(In)-(sum(In(In<0))+max(In)); %initial propulsion
    out.GRF_IMPULSES.(ctx+"AnteroPosteriorBraking")(j)=sum(In(In<0));
    out.GRF_IMPULSES.(ctx+"AnteroPosteriorPropulsive")(j)=max(In);
    out.GRF_IMPULSES.(ctx+"AnteroPosteriorNet")(j)=sum(In);
    
end
    catch ME
        warning(ME.message)
    end
end