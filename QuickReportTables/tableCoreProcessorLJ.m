function [Tfile]=tableCoreProcessorLJ(c3dfile,Info,tablepath,exmode)
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
    btkCloseAcquisition(H);
    
    %fix GRF
    fp=fp.cleanSignal;
    fp=fp.combineFP();
    fp=fp.getGRF;

try
    [kinDB]=getParamGRFv2(fp,ev,mass);
catch ME
    [~,f]=fileparts(c3dfile);
    warning("%s: %s",f ,ME.message);
end

events=ev.exportEvents('point',false);
try
JT.deltaD=deltaDistance_computation(markers, Info.Athlete.AmputationSide, takeoff_leg, events.(takeoff_leg).Foot_Off(end));
%angles computation
[JT.trunk_resampled,JT.trunk_angle_TD,JT.trunk_angle_TO,JT.thigh_angle_TD, JT.shank_angle_TD, JT.GT_tip_angle_TD, JT.GT_tip_angle_TO, JT.angleGTsec] = angles_computation(markers, Info.Athlete.AmputationSide, Info.Athlete.AmpLevel, takeoff_leg, events.(takeoff_leg).Foot_Strike(end), events.(takeoff_leg).Foot_Off(end));
%heights computation
[JT.H_GT_TD, JT.H_GT_TO, JT.H_RGT_max, JT.H_LGT_max, JT.H_GT_mean_2ndLast, JT.GTy, JT.GTz, JT.GTy_resampled, JT.GTz_resampled] = computation_height(markers, takeoff_leg, events);
%speed
[JT.Vh_TD, JT.Vh_TO, JT.Vv_TD, JT.Vv_TO, JT.Vh_mean_2ndLast, JT.Vv_mean_2ndLast] = computation_speed(markers, takeoff_leg, events, fs_mkr);
%Lh GT
JT.Lh_GT = computation_Lh(markers, takeoff_leg, events);
catch
end
Tab=makeQuickReportTableLJ(Info,kinDB,[],[],JT);
Tfile=writeSessionTableXLSX(Tab,fullfile(tablepath),exmode,'JUMP');

fp=fp.resample(fs_mkr);
events=ev.exportEvents('point');
%% GT-COP filling
if isfield(markers,[takeoff_leg(1) 'GT'])
lumped.Stiffness.(takeoff_leg)=lumpedStiffness(markers.([takeoff_leg(1) 'GT']),fp.COP,fp.GRF,events.(takeoff_leg).Foot_Strike,events.(takeoff_leg).Foot_Off,mass);
lumped.Theta.(takeoff_leg)=lumpedAngles(markers,[takeoff_leg(1) 'GT'],fp.COP,events.(takeoff_leg).Foot_Strike,events.(takeoff_leg).Foot_Off);

Tab=makeQuickReportTableLJ(Info,[],lumped.Theta,lumped.Stiffness);
Tfile=writeSessionTableXLSX(Tab,fullfile(tablepath),exmode,'JUMP');
else
warning('%s trajectory is missing',[takeoff_leg(1) 'GT'])
end
end



%% subfuns
