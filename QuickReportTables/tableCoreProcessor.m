function [Tfile,Tab]=tableCoreProcessor(c3dfile,Info,tablepath,exmode)
arguments
    c3dfile char
    Info struct
    tablepath char
    exmode char='';
end
    mass=Info.Athlete.Mass;

    % read c3d and create platforms and events
    H=btkReadAcquisition(c3dfile);
    fp=forcePlatformType2(c3dfile);
    markers=btkGetMarkers(H);
    fs_mkr=btkGetPointFrequency(H);
    ev=Event(H);
    if ev.getEventCount==0
        error("No events in the C3D!");
    end
    btkCloseAcquisition(H);
    
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
stepside(1:length(events.Left.Foot_Strike))={'Left'};
time=events.Left.Foot_Strike;
catch
stepside={};
end
try
stepside(end+1:end+length(events.Right.Foot_Strike))={'Right'};
time=[time events.Right.Foot_Strike];
catch
end
[~,o]=sort(time);
stepside=stepside(o);

%% FILL THE TABLE
% GRFs
try
[kinDB]=getParamGRFv2(fp,ev,mass);
catch ME
    [~,f]=fileparts(c3dfile);
    warning("%s: %s",f ,ME.message);
end
Tab=makeQuickReportTable(Info,stepside,kinDB,[],[]);
Tfile=writeSessionTableXLSX(Tab,fullfile(tablepath),exmode);

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
    speed(1)=mean(diff(sm/1000)*fs_mkr,'omitnan');
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
    speed(2)=mean(diff(sm/1000)*fs_mkr,'omitnan');
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