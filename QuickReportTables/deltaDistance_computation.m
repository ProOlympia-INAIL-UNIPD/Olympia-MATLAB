%% IDEA: calcolare automaticamente la distanza dal foot off a dove lascia il piede
function distance = deltaDistance_computation(markersData, amputationSide, leg_takeoff, FO_kin)

% markersData -> una struct contenente i dati di tutti i marker
% amputationSide -> stringa indicante il lato amputato, può essere 'none', 'Right' 'Left''
% leg_takeoff -> tipo char, 'Right' o 'Left''
% FO_kin -> foot off, corrisponde al frame in cui avviene il foot off
% ouput -> distance between the foot off and the end of the force platform,
% calculated in [m]

end_point=900; %[mm], indica la fine della pedana

%considero diversi casi in base al lato amputato e al piede di stacco
%considero la punta della protesi per il lato protesico e la seconda testa metatarsale per il lato sano

FO_kin=max(FO_kin); %considero solo il passo di stacco

if strcmpi(amputationSide, 'Right') && strcmpi(leg_takeoff, 'Right')
    try
        start_point = markersData.RFD1(FO_kin,2);
    catch
        start_point = (markersData.RFD2(FO_kin,2)+markersData.RFD3(FO_kin,2))/2 + 50; %se non c'è FD1 nella prova, ne calcolo la posizione come punto medio tra FD2 e FD3 +5cm
    end

elseif strcmpi(amputationSide, 'Left') && strcmpi(leg_takeoff, 'Left')
    try
        start_point = markersData.LFD1(FO_kin,2);
    catch
        start_point = (markersData.LFD2(FO_kin,2)+markersData.LFD3(FO_kin,2))/2 + 50;
    end

elseif (strcmpi(amputationSide, 'Left') && strcmpi(leg_takeoff, 'Right')) || (strcmpi(amputationSide, 'none') && strcmpi(leg_takeoff, 'Right'))
    try
        start_point = markersData.RIIMH(FO_kin,2);
    catch
        start_point = markersData.RIMH(FO_kin,2);
    end

elseif (strcmpi(amputationSide, 'Right') && strcmpi(leg_takeoff, 'Left')) || (strcmpi(amputationSide, 'none') && strcmpi(leg_takeoff, 'Left'))

    try
        start_point = markersData.LIIMH(FO_kin,2);
    catch
        start_point = markersData.LIMH(FO_kin,2);
    end
end



distance=(end_point-start_point) * 1e-3;


