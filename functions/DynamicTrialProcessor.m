function [trial,stat]=DynamicTrialProcessor(file,NameValue)
%This function executes the analysis on a dynamic file using a predefined
%set of instructions
arguments
    file char=[];
    %general
    NameValue.Configuration char=[];
    NameValue.StaticTrial char=[];
    NameValue.TrialUnits string {mustBeMember(NameValue.TrialUnits,["m","mm"])} ="mm";
    NameValue.RotateGround=eye(4);
    NameValue.GravityVector=[0 0 -9.806];
    %forceplat
    NameValue.MaxRadius=200;      %expected half-duration of the contact (in samples) 
    NameValue.FilterCutOff=100;   %cut-off frequency for low pass filter
    NameValue.FilterOrder=2;      %low pass filter order (will be doubled in filtfilt)
    NameValue.ActiveThreshold=100;%minimum threshold to accept FP as used during the acquisition
    NameValue.MaxNumContacts=100;
    NameValue.CombineFP='all';
    NameValue.Reflect=true;
    %kinematic
    NameValue.ClearUnlabeled='C_';
    NameValue.KinematicFilterOrder double=2;
    NameValue.KinematicFilterCutOff double=20;
    NameValue.SVDmode {mustBeMember(NameValue.SVDmode,["replace","ignore"])}="replace";
    %other
    NameValue.EditEvents=false;
    NameValue.EventTreshold=20;
    % NameValue.
end
if isempty(file)
   [f,p]=uigetfile('.c3d');
   file=fullfile(p,f);
end
if isempty(NameValue.Configuration)
   trial=Trial(file);
else
   trial=Trial(file,NameValue.Configuration);
end
if NameValue.EditEvents
   trial.Events=trial.Events.detectfromForcePlates("overwrite",NameValue.EventTreshold);
end
trial=trial.setUnits(NameValue.TrialUnits);

trial=trial.changeCoordinates(NameValue.RotateGround);

trial.ForcePlatform=trial.ForcePlatform.cleanSignal(...
    "FilterCutOff",NameValue.FilterCutOff,"FilterOrder",NameValue.FilterOrder,...
    "ActiveThreshold",NameValue.ActiveThreshold,"MaxNumContacts",...
    NameValue.MaxNumContacts,"MaxRadius",NameValue.MaxRadius,"Reflect",NameValue.Reflect);

if isempty(NameValue.CombineFP)
elseif isequal(NameValue.CombineFP,'all')
   trial.ForcePlatform=trial.ForcePlatform.combineFP;
else
   trial.ForcePlatform=trial.ForcePlatform.combineFP(NameValue.CombineFP);
end

if isempty(NameValue.ClearUnlabeled)
else
    trial.Points=trial.Points.clearUnlabeled(NameValue.ClearUnlabeled);
end

if isempty(NameValue.KinematicFilterCutOff)||NameValue.KinematicFilterOrder==0
else
   [b,a]=butter(NameValue.KinematicFilterOrder,NameValue.KinematicFilterCutOff/(trial.Metadata.POINT.RATE/2));
   trial.Points=trial.Points.filtfilt(b,a);
end

if isempty(NameValue.StaticTrial)
else
   stat=Trial(NameValue.StaticTrial,NameValue.Configuration);
   stat=stat.setUnits(NameValue.TrialUnits);
   stat=stat.buildSkeleton;
   %stat=stat.scaleInertialProp;
   stat=readC3DInertialProperties(stat);
   stat=stat.setUnits(NameValue.TrialUnits);
   trial=trial.segmentSolidification(stat,NameValue.SVDmode);
end


if not(isempty(trial.ConfigFile))
trial=trial.buildSkeleton;
trial=trial.moveInertialProperties(stat);
trial=trial.inverseDynamics(1,NameValue.GravityVector);
trial=lumpedAnalysis(trial,"LGT","RGT","LHJC","RHJC");
end
trial=trial.GRFAnalysis;
