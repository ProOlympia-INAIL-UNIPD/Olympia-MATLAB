function [Tfile,Tab]=tableCoreProcessor(c3dfile,Info,tablepath,exmode)
arguments
    c3dfile char
    Info struct
    tablepath char
    exmode char='';
end
    mass=Info.Athlete.Mass;

    % read c3d and create platforms and events
    trial=Trial(c3dfile);
    fp=trial.ForcePlatform;
    markers=trial.Points.PointStruct;
    fs_mkr=trial.Metadata.POINT.RATE;
    ev=trial.Events;
    
    if ev.getEventCount==0
        error("No events in the C3D!");
    end
    ev=ev.selectFootContacts;
    
    %fix GRF
    fp=fp.cleanSignal;
    fp=fp.combineFP();
    fp=fp.getGRF;

    % export events
    events=ev.exportEvents('analog');
    
%% ORGANIZE STEP ORDER
stepside={};
time=[];
try
    fs=events.Left.Foot_Strike;
    fo=events.Left.Foot_Off;
    nsteps=min(numel(fs),numel(fo));
    stepside(1:nsteps)={'Left'};
    time=fs(1:nsteps);
catch
stepside={};
end
try
    fs=events.Right.Foot_Strike;
    fo=events.Right.Foot_Off;
    nsteps=min(numel(fs),numel(fo));
    stepside(end+1:end+nsteps)={'Right'};
    time=[time fs(1:nsteps)];
catch
end
[~,o]=sort(time);
stepside=stepside(o);

%% FILL THE TABLE
% GRFs
try
[kinDB]=getParamGRFv2(fp,ev,mass);
catch ME
    kinDB=[];
    [~,f]=fileparts(c3dfile);
    warning("%s: %s",f ,ME.message);
end
% Tab=makeQuickReportTable(Info,stepside,kinDB,[],[]);
% Tfile=writeSessionTableXLSX(Tab,fullfile(tablepath),exmode);

lumpOK=false;

% Lumped kinematics
fp=fp.resample(fs_mkr);
events=ev.exportEvents('point');
if isfield(markers,'LGT') && isfield(events,'Left') %right side
    lumped.Stiffness.Left=lumpedStiffness(markers.LGT,fp.COP,fp.GRF,events.Left.Foot_Strike,events.Left.Foot_Off,mass);
    lumped.Theta.Left=lumpedAngles(markers,'LGT',fp.COP,events.Left.Foot_Strike,events.Left.Foot_Off);
    lumpOK=true;
    sm=markers.LGT;
    sm(sm==0)=nan;
    speed(1)=mean(vecnorm(diff(sm/1000)*fs_mkr,2,2),'omitnan');
else
    speed(1)=nan;
    warning('LGT trajectory is missing')
end

if isfield(markers,'RGT') && isfield(events,'Right') %left side
    lumped.Stiffness.Right=lumpedStiffness(markers.RGT,fp.COP,fp.GRF,events.Right.Foot_Strike,events.Right.Foot_Off,mass);
    lumped.Theta.Right=lumpedAngles(markers,'RGT',fp.COP,events.Right.Foot_Strike,events.Right.Foot_Off);
    lumpOK=true;
    sm=markers.RGT;
    sm(sm==0)=nan;
    speed(2)=mean(vecnorm(diff(sm/1000)*fs_mkr,2,2),'omitnan');
else
    speed(2)=nan;
    warning('RGT trajectory is missing')
end
speed=mean(speed,'omitnan');
if lumpOK %fill kinematics table
    Tab=makeQuickReportTable(Info,stepside,kinDB,lumped.Theta,lumped.Stiffness,speed);
    Tfile=writeSessionTableXLSX(Tab,fullfile(tablepath),exmode);
end
end