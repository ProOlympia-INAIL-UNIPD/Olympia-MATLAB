%% Height GT

function [H_GT_FS, H_GT_FO, H_RGT_max, H_LGT_max, H_GT_mean_beforeLast, GTy, GTz, GTy_cycle, GTz_cycle]= computation_height(markersData, takeoff_leg, eventsKin)

debugplot=false;
% OUTPUTS:
% H_GT_FS --> GT height at touchdown
% H_GT_FO --> GT height at take-off
% H_RGT_max --> RGT max height
% H_LGT_max --> LGT max height
%GTy_resampled --> GTy values during the stance phase
%GTy_resampled --> GTy values during the stance phase
%GTy --> GTy values for all the frames
%GTz --> GTz values for all the frames
%Outputs expressed in metres

LGT_heights=markersData.LGT(:,3) * 1e-3;  %[m]

RGT_heights=markersData.RGT(:,3) * 1e-3;  %[m]

H_RGT_max=max(RGT_heights);  
H_LGT_max=max(LGT_heights);  

%Se ho più eventi di foot strike e foot off scelgo l'ultimo che
%corrisponmde allo stacco
FS=max(eventsKin.(takeoff_leg).Foot_Strike);
FO=max(eventsKin.(takeoff_leg).Foot_Off);

leg=["Left","Right"];
secondLast=leg(not(leg==takeoff_leg));

FS_beforeLast =max(eventsKin.(secondLast).Foot_Strike);
FO_beforeLast =max(eventsKin.(secondLast).Foot_Off);

if(strcmpi(takeoff_leg, 'Right'))
    H_GT_FO = RGT_heights(FO);  
    H_GT_FS = RGT_heights(FS);
    H_GT_mean_beforeLast=(LGT_heights(FO_beforeLast)+LGT_heights(FS_beforeLast))/2;
    

    GTy=markersData.RGT(:,2) * 1e-3; %[m]
    GTz=RGT_heights;
else
    H_GT_FO = LGT_heights(FO);  
    H_GT_FS = LGT_heights(FS); 
    H_GT_mean_beforeLast=(RGT_heights(FO_beforeLast)+RGT_heights(FS_beforeLast))/2;

    GTy=markersData.LGT(:,2) * 1e-3; %[m]
    GTz=LGT_heights;
end


GTy_FSFO=GTy(FS:FO);
GTy_cycle=time2cycle([],GTy_FSFO,101);

GTz_FSFO=GTz(FS:FO);
GTz_cycle=time2cycle([],GTz_FSFO,101);

if debugplot
figure(), 
subplot(221), plot(GTy, 'LineWidth', 2)
hold on
plot(FS, GTy(FS), 'r*', 'MarkerSize', 10, 'LineWidth', 2)
plot(FO, GTy(FO), 'g*', 'MarkerSize', 10, 'LineWidth', 2)
xlabel('Number of frames')
ylabel('Distance [m]')
legend('GTy trend', 'Foot Strike', 'Foot Off', 'Location', 'southeast')
title('Anteroposterior GT trend of the takeoff leg in long jumping')
xlim([1 length(GTy)])

subplot(223), plot(GTy_100, 'LineWidth', 2)
hold on
plot(1, GTy_100(1), 'r*', 'MarkerSize', 10, 'LineWidth', 2);
plot(101, GTy_100(end), 'g*', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('% stance');
ylabel('Distance [m]');
xlim([1 101])
legend('GTy trend', 'Foot Strike', 'Foot Off', 'Location', 'southeast');
title('AP GT trend of the takeoff leg during the stance phase');

subplot(222), plot(GTz, 'LineWidth', 2)
hold on
plot(FS, H_GT_FS, 'r*', 'MarkerSize', 10, 'LineWidth', 2)
plot(FO, H_GT_FO, 'g*', 'MarkerSize', 10, 'LineWidth', 2)
xlabel('Number of frames')
ylabel('Height [m]')
legend('GTz trend', 'Foot Strike', 'Foot Off', 'Location', 'southeast')
title('Vertical GT trend of the takeoff leg in long jumping')
xlim([1 length(GTz)])

subplot(224), plot(GTz_cycle, 'LineWidth', 2)
hold on
plot(1, GTz_cycle(1), 'r*', 'MarkerSize', 10, 'LineWidth', 2)
plot(101, GTz_cycle(end), 'g*', 'MarkerSize', 10, 'LineWidth', 2)
xlabel('% stance')
ylabel('Height [m]')
xlim([1 101])
legend('GTz trend', 'Foot Strike', 'Foot Off', 'Location', 'southeast')
title('Vertical GT trend of the takeoff leg during the stance phase')
end


