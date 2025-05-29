function [stat]=StaticTrialProcessor(file,config,NameValue)
%This function executes the analysis on a static file using a predefined
%set of instructions
arguments
    file char=[];
    config =[];
    %general
    NameValue.TrialUnits="mm";
    NameValue.ChangeCoordinates=eye(4);
    NameValue.ExportSuffix='_PROCESSED';
    NameValue.AverageTrial=true;
    NameValue.XMLatt="_att";
    %forceplat
    NameValue.CombineFP=nan;
    %kinematic
    NameValue.ClearUnlabeled='C_';
end

if isempty(config)
stat=Trial(file);
elseif isstruct(config)
stat=Trial(file);
stat.XMLatt=NameValue.XMLatt;
stat.applyConfiguration(config);


elseif isfile(config)
stat=Trial(file,config);
stat.XMLatt=NameValue.XMLatt;
end
if NameValue.AverageTrial
   stat=mean(stat);
end

stat=stat.changeCoordinates(NameValue.ChangeCoordinates);
if isequal(NameValue.CombineFP,nan)
else
stat.ForcePlatform=stat.ForcePlatform.combineFP([],[],NameValue.CombineFP);
end
stat.Points=stat.Points.clearUnlabeled(NameValue.ClearUnlabeled);

stat=staticMarkerReconstruction(stat);
stat=stat.setUnits(NameValue.TrialUnits);
stat=stat.buildSkeleton;
stat=stat.setSubjectAntropometry();
stat=stat.scaleInertialProp;
end