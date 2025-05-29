function Lh_GT= computation_Lh(markersData, takeoff_leg, eventsKin)

% OUTPUTS:
%Horizontal distance covered by GT during the stance phase of the last step
%Output expressed in metres

%Estraggo i dati in direzione anteroposteriore dei due GT
LGT_ap=markersData.LGT(:,2) * 1e-3;  %[m]
RGT_ap=markersData.RGT(:,2) * 1e-3;  %[m]

%Se ho più eventi di foot strike e foot off scelgo l'ultimo che
%corrisponmde allo stacco
FS=max(eventsKin.(takeoff_leg).Foot_Strike);
FO=max(eventsKin.(takeoff_leg).Foot_Off);

if(strcmpi(takeoff_leg, 'Right'))
    GT_FO = RGT_ap(FO);  
    GT_FS = RGT_ap(FS);
    
else
    GT_FO = LGT_ap(FO);  
    GT_FS = LGT_ap(FS);
    
end

Lh_GT=GT_FO-GT_FS;
