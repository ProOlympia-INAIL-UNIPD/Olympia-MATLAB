%% Speed GT

function [Vh_FS, Vh_FO, Vv_FS, Vv_FO, Vh_mean_beforeLast, Vv_mean_beforeLast]= computation_speed(markersData, takeoff_leg, eventsKin, Fs_mkr)

debugplot=false
% OUTPUTS:
% Vh_FS -> horizontal velocity at touchdown
% Vh_FO -> horizontal velocity at takeoff
% Vv_FS -> vertical velocity at touchdown
% Vv_FO -> vertical velocity at takeoff
% Vh_mean_beforeLast -> media della velocità orizzontale tra FS e FO del
% penultimo passo
% Vv_mean_beforeLast -> media della velocità verticale tra FS e FO del
% penultimo passo
%Outputs expressed in [m/s]


RGT_pos_h=markersData.RGT(:,2) * 1e-3;  %[m], y is the anteroposterior direction
LGT_pos_h=markersData.LGT(:,2) * 1e-3;  %[m],

RGT_pos_v=markersData.RGT(:,3) * 1e-3;  %[m], z is the vertical direction
LGT_pos_v=markersData.LGT(:,3) * 1e-3;  %[m]

RGT_Vh=diff([0; RGT_pos_h])*Fs_mkr; %[m/s]
RGT_Vv=diff([0; RGT_pos_v])*Fs_mkr; 

LGT_Vh=diff([0; LGT_pos_h])*Fs_mkr;
LGT_Vv=diff([0; LGT_pos_v])*Fs_mkr;

%Per La Barbera le velocità del penultimo passo sono state calcolate 
%senza lo 0 prima di inserire il vettore posizione

%Se ho più eventi di foot strike e foot off scelgo l'ultimo che
%corrisponde allo stacco
FS=max(eventsKin.(takeoff_leg).Foot_Strike);
FO=max(eventsKin.(takeoff_leg).Foot_Off);

leg=["Left","Right"];
secondLast=leg(not(leg==takeoff_leg));

try
FS_beforeLast =max(eventsKin.(secondLast).Foot_Strike);
FO_beforeLast =max(eventsKin.(secondLast).Foot_Off);
sl=true;
catch
sl=false;
end

if(strcmpi(takeoff_leg, 'R'))
    
    Vh_FS = RGT_Vh(FS);
    Vh_FO = RGT_Vh(FO);
    Vv_FS = RGT_Vv(FS);
    Vv_FO = RGT_Vv(FO);
    if sl
    Vh_mean_beforeLast=(LGT_Vh(FO_beforeLast)+LGT_Vh(FS_beforeLast))/2;
    Vv_mean_beforeLast=(LGT_Vv(FO_beforeLast)+LGT_Vv(FS_beforeLast))/2;
    else
    Vh_mean_beforeLast=nan;
    Vv_mean_beforeLast=nan;
    end

    Vh=RGT_Vh;
    Vv=RGT_Vv;

else
    Vh_FS = LGT_Vh(FS);
    Vh_FO = LGT_Vh(FO);
    Vv_FS = LGT_Vv(FS);
    Vv_FO = LGT_Vv(FO);
    if sl
    Vh_mean_beforeLast=(RGT_Vh(FO_beforeLast)+RGT_Vh(FS_beforeLast))/2;
    Vv_mean_beforeLast=(RGT_Vv(FO_beforeLast)+RGT_Vv(FS_beforeLast))/2;
    else
    Vh_mean_beforeLast=nan;
    Vv_mean_beforeLast=nan;
    end
    Vh=LGT_Vh;
    Vv=LGT_Vv;
end

Vh_FSFO=Vh(FS:FO);
Vh_cycle=time2cycle([],Vh_FSFO,101);

Vv_FSFO=Vv(FS:FO);
Vv_cycle=time2cycle([],Vv_FSFO,101);

if debugplot
figure(), subplot(221), plot(Vh(2:end), 'LineWidth', 2)
hold on
plot(FS, Vh_FS, 'r*', 'MarkerSize', 10, 'LineWidth', 2);
plot(FO, Vh_FO, 'g*', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Number of frames');
ylabel('Speed [m/s]');
legend('GT_{h} speed', 'Foot Strike', 'Foot Off', 'Location', 'southeast');
title('GT horizontal speed of the takeoff leg in long jumping');
xlim([1 length(Vh)])

subplot(223), plot(Vh_cycle, 'LineWidth', 2)
hold on
plot(1, Vh_cycle(1), 'r*', 'MarkerSize', 10, 'LineWidth', 2);
plot(101, Vh_cycle(end), 'g*', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('% stance');
ylabel('Speed [m/s]');
xlim([1 101])
legend('GT_{h} speed', 'Foot Strike', 'Foot Off', 'Location', 'southeast');
title('GT horizontal speed of the takeoff leg during the stance phase');

subplot(222), plot(Vv(2:end), 'LineWidth', 2)
hold on
plot(FS, Vv_FS, 'r*', 'MarkerSize', 10, 'LineWidth', 2);
plot(FO, Vv_FO, 'g*', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Number of frames');
ylabel('Speed [m/s]');
legend('GT_{v} speed', 'Foot Strike', 'Foot Off', 'Location', 'southeast');
title('GT vertical speed of the takeoff leg in long jumping');
xlim([1 length(Vv)])

subplot(224), plot(Vv_cycle, 'LineWidth', 2)
hold on
plot(1, Vv_cycle(1), 'r*', 'MarkerSize', 10, 'LineWidth', 2);
plot(101, Vv_cycle(end), 'g*', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('% stance');
ylabel('Speed [m/s]');
xlim([1 101])
legend('GT_{v} speed', 'Foot Strike', 'Foot Off', 'Location', 'southeast');
title('GT vertical speed of the takeoff leg during the stance phase');

end

