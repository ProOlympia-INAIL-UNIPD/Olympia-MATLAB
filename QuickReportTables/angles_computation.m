%% Angles computation
function [trunk_resampled,trunk_angle_FS,trunk_angle_FO,thigh_angle_FS,shank_angle_FS,GT_tip_angle_FS,GT_tip_angle_FO,sec_gitt_angle] = angles_computation(markersData, amputationSide, amputationLevel, takeoffLeg, FS_kin, FO_kin)
arguments
    markersData struct
    amputationSide char {mustBeMember(amputationSide,{'Left','Right','none'})}
    amputationLevel char {mustBeMember(amputationLevel,{'TT','TF','none'})}
    takeoffLeg char {mustBeMember(takeoffLeg,{'Left','Right'})}
    FS_kin double
    FO_kin double
end
debugplot=false;
% markersData -> una struct contenente i dati di tutti i marker
% amputationSide -> stringa indicante il lato amputato, può essere 'none', 'Right', 'Left'
% amputationLevel -> stringa indicante il lato amputato, può essere 'none', 'TT', 'TF'
% takeoffLeg -> gamba di stacco, può essere 'Right' o 'Left'
% FO_kin -> foot off, corrisponde al frame in cui avviene il take off
% FS_kin -> foot strike, corrisponde al frame in cui avviene il touchdown

%se ho gli eventi di più passi selezioni il frame più grande che
%corrisponde allo stacco
FO_kin=max(FO_kin);
FS_kin=max(FS_kin);

%% TRUNK
%segmento che unisce il punto medio di C7-IJ con il punto medio di T8-PX
try
    %coordinate traiettoria punto medio T8-PX
    y_T8_PX = (markersData.T8(:,2) + markersData.PX(:,2)) / 2;
    z_T8_PX = (markersData.T8(:,3) + markersData.PX(:,3)) / 2;

    %coordinate traiettoria punto medio C7-IJ
    y_C7_IJ = (markersData.C7(:,2) + markersData.IJ(:,2)) / 2;
    z_C7_IJ = (markersData.C7(:,3) + markersData.IJ(:,3)) / 2;

    %proiezioni verticale e orizzontale per ogni frame rispetto cui devo calcolare
    %l'inclinazione
    vertical_tr_length = z_C7_IJ-z_T8_PX;
    horizontal_tr_length= y_C7_IJ-y_T8_PX;

    %TRUNK INCLINATION at touchdown (FS) and takeoff (FO)
    trunk_inclination=atan2d(horizontal_tr_length,vertical_tr_length);
    trunk_angle_FS = trunk_inclination(FS_kin); %touchdown
    trunk_angle_FO = trunk_inclination(FO_kin); %take-off

    trunk_angles_FSFO=trunk_inclination(FS_kin:FO_kin);
    trunk_resampled=time2cycle([],trunk_angles_FSFO,101);
    if debugplot
    figure(), subplot(211), plot(trunk_inclination, 'LineWidth', 2)
    hold on
    plot(FS_kin, trunk_angle_FS, 'r*', 'MarkerSize', 10, 'LineWidth', 2);
    plot(FO_kin, trunk_angle_FO, 'g*', 'MarkerSize', 10, 'LineWidth', 2);
    xlabel('Number of frames');
    ylabel('Angles [°]');
    legend('Trunk inclination', 'Foot Strike', 'Foot Off');
    title('Trunk inclination in long jumping');
    xlim([1 length(trunk_inclination)])

    subplot(212), plot(trunk_resampled, 'LineWidth', 2)
    hold on
    plot(1, trunk_resampled(1), 'r*', 'MarkerSize', 10, 'LineWidth', 2);
    plot(101, trunk_resampled(end), 'g*', 'MarkerSize', 10, 'LineWidth', 2);
    xlabel('% stance');
    ylabel('Angles [°]');
    xlim([1 101])
    legend('Trunk inclination', 'Foot Strike', 'Foot Off');
    title('Trunk inclination during the stance phase in long jumping');
    end
catch
    trunk_resampled=NaN;
    trunk_angle_FO=NaN;
    trunk_angle_FS=NaN;
end
clearvars 'marker'
%% THIGH
%none and TT -> calcolo dell'angolo del segmento GT-FLE rispetto alla
%verticale
%TF -> angolo del segmento che unisce GT a LAD rispetto alla
%verticale
try 
    if(strcmpi(amputationLevel, 'none') || strcmpi(amputationLevel, 'TT')) %in questi due casi non ce ne importa niente del lato amputato perché i marker fino alla coscia hanno lo stesso nome
        if(strcmpi(takeoffLeg, 'Right'))
            vertical_th_length = markersData.RGT(:,3)-markersData.RFLE(:,3);
            horizontal_th_length = markersData.RGT(:,2)-markersData.RFLE(:,2);
        else
            vertical_th_length = markersData.LGT(:,3)-markersData.LFLE(:,3);
            horizontal_th_length = markersData.LGT(:,2)-markersData.LFLE(:,2);
        end

    else %TF
        if(strcmpi(takeoffLeg, 'Right') && strcmpi(amputationSide, 'Right'))

            vertical_th_length = markersData.RGT(:,3)-markersData.RLAD(:,3);
            horizontal_th_length = markersData.RGT(:,2)-markersData.RLAD(:,2);

        elseif(strcmpi(takeoffLeg, 'Left') && strcmpi(amputationSide, 'Left'))

            vertical_th_length = markersData.LGT(:,3)-markersData.LLAD(:,3);
            horizontal_th_length = markersData.LGT(:,2)-markersData.LLAD(:,2);

        elseif(strcmpi(takeoffLeg, 'Left') && strcmpi(amputationSide, 'Right'))

            vertical_th_length = markersData.LGT(:,3)-markersData.LFLE(:,3);
            horizontal_th_length = markersData.LGT(:,2)-markersData.LFLE(:,2);

        elseif(strcmpi(takeoffLeg, 'Right') && strcmpi(amputationSide, 'Left'))

            vertical_th_length = markersData.RGT(:,3)-markersData.RFLE(:,3);
            horizontal_th_length = markersData.RGT(:,2)-markersData.RFLE(:,2);
        end

    end

    %Thigh inclination at touchdown
    thigh_angle=atan2d(horizontal_th_length,vertical_th_length);
    thigh_angle_FS=thigh_angle(FS_kin);
    %thigh_angles_FSFO=thigh_angle(FS_kin:FO_kin);
    %thigh_resampled=time2cycle([],thigh_angles_FSFO,1001);
    if debugplot
    figure(), plot(thigh_angle, 'LineWidth', 2)
    hold on
    plot(FS_kin, thigh_angle_FS, 'r*', 'MarkerSize', 10, 'LineWidth', 2);
    plot(FO_kin, thigh_angle(FO_kin), 'g*', 'MarkerSize', 10, 'LineWidth', 2);
    xlabel('Number of frames');
    ylabel('Angles [°]');
    legend('Thigh inclination', 'Foot Strike', 'Foot Off');
    title('Thigh inclination of the takeoff leg in long jumping');
    xlim([1 length(thigh_angle)])

    % subplot(212), plot(thigh_resampled, 'LineWidth', 2)
    % hold on
    % plot(0, thigh_resampled(1), 'r*', 'MarkerSize', 10, 'LineWidth', 2);
    % plot(1001, thigh_resampled(end), 'g*', 'MarkerSize', 10, 'LineWidth', 2);
    % xlabel('% stance');
    % ylabel('Angles [°]');
    % xlim([0 1001])
    % legend('Thigh inclination', 'Foot Strike', 'Foot Off');
    % title('Thigh inclination during the stance phase in long jumping');
    end
catch
    thigh_angle_FS=NaN;
end
%% SHANK

%none -> calcolo shank come FLE-LM
%TF -> calcolo angolo del clamp rispetto alla verticale
%TT -> uso marker del socket rispetto alla verticale
try
    if(strcmpi(amputationLevel, 'none'))
        if(strcmpi(takeoffLeg, 'Right'))
            vertical_sh_length = markersData.RFLE(:,3)-markersData.RLM(:,3);
            horizontal_sh_length =markersData.RFLE(:,2)-markersData.RLM(:,2);
        else
            vertical_sh_length = markersData.LFLE(:,3)-markersData.LLM(:,3);
            horizontal_sh_length= markersData.LFLE(:,2)-markersData.LLM(:,2);
        end
        shank_angle=atan2d(horizontal_sh_length, vertical_sh_length);


    elseif(strcmpi(amputationLevel, 'TT'))

        if(strcmpi(takeoffLeg, 'Right') && strcmpi(amputationSide, 'Right'))

            vertical_sh_length = markersData.RLAP(:,3)-markersData.RLAD(:,3);
            horizontal_sh_length = markersData.RLAP(:,2)-markersData.RLAD(:,2);

        elseif(strcmpi(takeoffLeg, 'Left') && strcmpi(amputationSide, 'Left'))

            vertical_sh_length = markersData.LLAP(:,3)-markersData.LLAD(:,3);
            horizontal_sh_length = markersData.LLAP(:,2)-markersData.LLAD(:,2);

        elseif(strcmpi(takeoffLeg, 'Left') && strcmpi(amputationSide, 'Right'))

            vertical_sh_length = markersData.LFLE(:,3)-markersData.LLM(:,3);
            horizontal_sh_length = markersData.LFLE(:,2)-markersData.LLM(:,2);

        elseif(strcmpi(takeoffLeg, 'Right') && strcmpi(amputationSide, 'Left'))

            vertical_sh_length = markersData.RFLE(:,3)-markersData.RLM(:,3);
            horizontal_sh_length=markersData.RFLE(:,2)-markersData.RLM(:,2);
        end

        shank_angle=atan2d(horizontal_sh_length, vertical_sh_length);

    else %TF -> calcolo clamp rispetto a ORIZZONTALE
        if(strcmpi(takeoffLeg, 'Right') && strcmpi(amputationSide, 'Right'))

            y_cl_prox_medio = (markersData.RCLAM(:,2) + markersData.RCLAL(:,2)) / 2;
            z_cl_prox_medio = (markersData.RCLAM(:,3) + markersData.RCLAL(:,3)) / 2;

            y_cl_dist_medio = (markersData.RCLPM(:,2) + markersData.RCLPL(:,2)) / 2;
            z_cl_dist_medio = (markersData.RCLPM(:,3) + markersData.RCLPL(:,3)) / 2;

            %shank_length -> vettore che collega il punto medio tra i marker prossimali del clamp e il punto medio dei marker distali del clamp
            vertical_sh_length = z_cl_prox_medio-z_cl_dist_medio;
            horizontal_sh_length = y_cl_prox_medio-y_cl_dist_medio;

            shank_angle=atan2d(vertical_sh_length, horizontal_sh_length); %rispetto a orizzontale

        elseif(strcmpi(takeoffLeg, 'Left') && strcmpi(amputationSide, 'Left'))

            y_cl_prox_medio = (markersData.LCLAM(:,2) + markersData.LCLAL(:,2)) / 2;
            z_cl_prox_medio = (markersData.LCLAM(:,3) + markersData.LCLAL(:,3)) / 2;

            y_cl_dist_medio = (markersData.LCLPM(:,2) + markersData.LCLPL(:,2)) / 2;
            z_cl_dist_medio = (markersData.LCLPM(:,3) + markersData.LCLPL(:,3)) / 2;

            %shank_length -> vettore che collega il punto medio tra i marker prossimali del clamp al punto medio dei marker distali del clamp
            vertical_sh_length = z_cl_prox_medio-z_cl_dist_medio;
            horizontal_sh_length = y_cl_prox_medio-y_cl_dist_medio;

            shank_angle=atan2d(vertical_sh_length, horizontal_sh_length); %rispetto a orizzontale

        elseif(strcmpi(takeoffLeg, 'Left') && strcmpi(amputationSide, 'Right'))

            vertical_sh_length = markersData.LFLE(:,3)-markersData.LLM(:,3);
            horizontal_sh_length = markersData.LFLE(:,2)-markersData.LLM(:,2);

            shank_angle=atan2d(horizontal_sh_length, vertical_sh_length);

        elseif(strcmpi(takeoffLeg, 'Right') && strcmpi(amputationSide, 'Left'))

            vertical_sh_length = markersData.RFLE(:,3)-markersData.RLM(:,3);
            horizontal_sh_length = markersData.RFLE(:,2)-markersData.RLM(:,2);

            shank_angle=atan2d(horizontal_sh_length, vertical_sh_length);
        end

    end

    shank_angle_FS=shank_angle(FS_kin);
    %shank_angles_FSFO=shank_angle(FS_kin:FO_kin);
    %shank_resampled=time2cycle([],shank_angles_FSFO,1001);
    if debugplot
    figure(), plot(shank_angle, 'LineWidth', 2)
    hold on
    plot(FS_kin, shank_angle_FS, 'r*', 'MarkerSize', 10, 'LineWidth', 2);
    plot(FO_kin, shank_angle(FO_kin), 'g*', 'MarkerSize', 10, 'LineWidth', 2);
    xlabel('Number of frames');
    ylabel('Angles [°]');
    legend('Shank inclination', 'Foot Strike', 'Foot Off');
    title('Shank inclination of the takeoff leg in long jumping');
    xlim([1 length(shank_angle)])
    end
catch
    shank_angle_FS=NaN;
end
%% GT-Tip protesi
try
    if(strcmpi(takeoffLeg, 'Left') && strcmpi(amputationSide, 'Left'))
        vertical_length = markersData.LGT(:,3)-markersData.LFD1(:,3);
        horizontal_length = markersData.LGT(:,2)-markersData.LFD1(:,2);

    elseif(strcmpi(takeoffLeg, 'Right') && strcmpi(amputationSide, 'Right'))

        vertical_length = markersData.RGT(:,3)-markersData.RFD1(:,3);
        horizontal_length=markersData.RGT(:,2)-markersData.RFD1(:,2);
    end
    GT_tip_angle_FS=atan2d(horizontal_length(FS_kin), vertical_length(FS_kin));
    GT_tip_angle_FO=atan2d(horizontal_length(FO_kin), vertical_length(FO_kin));
catch
    GT_tip_angle_FS=NaN;
    GT_tip_angle_FO=NaN;
end
%% angolo GT allo stacco-punto max gittata GT
try
    if strcmpi(takeoffLeg, 'Left')
        [~,pos]=max(markersData.LGT(:,3));
        vec=markersData.LGT(pos,:)-markersData.LGT(FO_kin,:);
        sec_gitt_angle=atan2d(vec(3),vec(2));
    elseif strcmpi(takeoffLeg, 'Right')
        [~,pos]=max(markersData.RGT);
        vec=markersData.RGT(pos,:)-markersData.RGT(FO_kin,:);
        sec_gitt_angle=atan2d(vec(3),vec(2));
    end
catch
    sec_gitt_angle=NaN;
end



