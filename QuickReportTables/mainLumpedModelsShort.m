[fileName, pathName ] = uigetfile( [TrialPath sepWinMac 'Task' sepWinMac 'Solidified' sepWinMac '*.c3d'], 'Choose c3d file for lumped parameters model');
h=btkReadAcquisition([pathName fileName]);
[markersData, ~]=btkGetMarkers(h);
btkCloseAcquisition(h);
trialType=split(GeneralInfo.trial_info,';');
switch trialType{2}
    case 'SSR'
        force=F_gnd(:,:,10);
    case'longJumping'
        force=F_gnd;
end
%load required markers
MIlump = fileReader([TrialPath sepWinMac 'ConfigFiles' sepWinMac 'MarkersInfo.txt']);
lumpK_L = MIlump.markers_lumpk_left;
lumpK_R = MIlump.markers_lumpk_right;
lumpAng_L = MIlump.markers_lumpang_left; %last is distal for theta
lumpAng_R = MIlump.markers_lumpang_right;%last is distal for theta
mass=str2double(GeneralInfo.session_info{2});%mass for normalized stiffness

%events from kinetic to kinematics and GRF downsampled to match kinematics
eventsKin=resampleEvents(eventsForce,Fs_analog,Fs_mkr);
GRFkin.F=force(1:Fs_ratio:end,:);
GRFkin.COP=COP_gnd_filt(1:Fs_ratio:end,:);
switch trialType{2}
    case 'SSR'
        try
            for i=1:length(lumpK_R)
                lumped(1).Stiffness.R.(lumpK_R{i}(2:end))=lumpedStiffness(markersData.(lumpK_R{i}),GRFkin.COP,GRFkin.F,eventsKin.Right_Foot_Strike,eventsKin.Right_Foot_Off,mass);
            end
        catch
        end
        try
            lumped(1).Stiffness.L.(lumpK_L{1}(2:end))=lumpedStiffness(markersData.(lumpK_L{1}),GRFkin.COP,GRFkin.F,eventsKin.Left_Foot_Strike,eventsKin.Left_Foot_Off,mass);
            lumped(1).Theta.R=lumpedAngles(markersData,lumpAng_R(1),GRFkin.COP,eventsKin.Right_Foot_Strike,eventsKin.Right_Foot_Off);
            lumped(1).Theta.L=lumpedAngles(markersData,lumpAng_L(1),GRFkin.COP,eventsKin.Left_Foot_Strike,eventsKin.Left_Foot_Off);
        catch
        end
        try
            writeSessionTable(GeneralInfo,eventsForce,[],lumped(1),athletePath)
        catch
        end
    case 'longJumping'
        if(strcmpi(leg_takeoff, 'R'))
            lumped(1).Theta.R=lumpedAngles(markersData,lumpAng_R(1),GRFkin.COP,eventsKin.Right_Foot_Strike(end),eventsKin.Right_Foot_Off(end));
            for i=1:length(lumpK_R)
                lumped(1).Stiffness.R.(lumpK_R{i}(2:end))=lumpedStiffness(markersData.(lumpK_R{i}),GRFkin.COP,GRFkin.F,eventsKin.Right_Foot_Strike(end),eventsKin.Right_Foot_Off(end),mass);
            end
        else
            for i=1:length(lumpK_L)
                lumped(1).Stiffness.L.(lumpK_L{i}(2:end))=lumpedStiffness(markersData.(lumpK_L{i}),GRFkin.COP,GRFkin.F,eventsKin.Left_Foot_Strike(end),eventsKin.Left_Foot_Off(end),mass);
            end
            lumped(1).Theta.L=lumpedAngles(markersData,lumpAng_L(1),GRFkin.COP,eventsKin.Left_Foot_Strike(end),eventsKin.Left_Foot_Off(end));
        end
        writeSessionTableLJ(GeneralInfo,[],lumped(1),[],'',athletePath)
end