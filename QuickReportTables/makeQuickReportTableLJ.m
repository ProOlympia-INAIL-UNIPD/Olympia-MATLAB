function out=makeQuickReportTableLJ(Info,kinDB,theta,stiffness,jumpdata,speed)
%Tfile=writeSessionTableV2(Info,kinDB,lumped,stepside,Tfile,existingdatamode,speed)
%creates/updates a Session Excel Spreadsheet containing the info, kinetics and/or lumped results
%for each step existing in the input trial
% the function splits the existing template table into three blocks: INFO, AMTI,
% LUMPED which can be filled/overwrited depending on user inputs
% Inputs (Required):
% Info: structure with fields: Athlete, Session, Trial
%       Athlete contains ID, mass, Amputation Side
%       Session contains ID, Date
%       Trial contains ID of the run and ID of the trial
% kinDB: struct containing forces, impulses and other kinetic measurements
% lumped: struct containing 1-D leg lumped model results
% stepside: cell array containing the order of step laterality
% Inputs (Optional):
% Tfile: directory in to which save the .xlsx file (defaults to
% current directory), or specific destination file
% existingdatamode: defines how the function behaves when user tries to
%                   write data of a step which was already in the table 
%                   -Ignore: skip existing data entries, write missing
%                            entries (if any new)
%                   -Overwrite: overwrite the FULL row in which data exist
%                   -Append: move to the first empty row and write again
%                            the data
% speed: the average speed during the trial
arguments
Info struct   =[];
kinDB struct  =[];
theta struct  =[];
stiffness struct  =[];
jumpdata struct  =[];
speed double =nan;
end

%% initialize missing input to empty
if isempty(kinDB)
   kinDB=struct('Left',[],'Right',[]);
end

if isempty(theta)
   theta=struct('Left',[],'Right',[]);
end

if isempty(stiffness)
   stiffness=struct('Left',[],'Right',[]);
end
cycles=nan(101,7);
%% Initialize Table 
[Tathl, Tprost, Ttrial]=infoTable(Info);

out.Parameters=table();    
side=Info.Trial.TakeOffLeg;
Ttrial.StepSide=string(side);
%% video+AMTI
TvideoAMTI=table(nan,nan,nan,'VariableNames',{'FP2land','TakeOff2FP','JumpLength'});
%FP2land: distanza (misurata in campo) tra bordo anteriore pedana e atterraggio
%Toff2FP: distanza tra COP al takeoff e bordo anteriore pedana
%JumpLength: lunghezza effettiva salto
if not(isempty(jumpdata))
    TvideoAMTI.FP2land=nan;
    TvideoAMTI.TakeOff2FP=jumpdata.deltaD;
    TvideoAMTI.JumpLength=nan;
end
%% AMTI
Tamti=table(nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,'VariableNames',{'meanGRFbrak','meanGRFprop',...
    'meanGRFvert','maxGRFbrak','maxGRFprop','maxGRFvert','Istprop','Ibrak','Iprop','Ivert'});

if not(isempty(kinDB.(side)))
    Tamti.meanGRFbrak=nan;
    Tamti.meanGRFprop=nan;
    Tamti.meanGRFvert=kinDB.(side).parameters.vertical.mean.cycle(end);

    Tamti.maxGRFbrak=kinDB.(side).parameters.horizontal.min.cycle(end);
    Tamti.maxGRFprop=kinDB.(side).parameters.horizontal.max.cycle(end);
    Tamti.maxGRFvert=kinDB.(side).parameters.vertical.max.cycle(end);

    Tamti.Istprop=kinDB.(side).impulse.braking.cycle(end);
    Tamti.Ibrak=kinDB.(side).impulse.braking.cycle(end);
    Tamti.Iprop=kinDB.(side).impulse.propulsive.cycle(end);
    Tamti.Ivert=kinDB.(side).impulse.vertical.cycle(end);
    cycles(:,1)=kinDB.(side).GRF.cycle(:,2,ii);
    cycles(:,2)=kinDB.(side).GRF.cycle(:,3,ii);
    cycles(:,3)=vecnorm(cycles(:,1:2),2,2);
    cycles(:,4:5)=kinDB.(side).COP.cycle(:,2:-1:1,ii);
end
%% VICON
Tvicon=table(nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,'VariableNames',{'GTh2ndlast','GThTD','GThTO','GTv2ndlast','GTvTD','GTvTO','thTrTD','thTrTO','thThTD','thShTD','hGT2ndlast','hGTTD','hGTTO','hLGTmax','hRGTmax','secGT'});
%'GTh2ndlast: media velocità del GT orizzontale nel penultimo passo
% GThTD: velocità orizzontale GT al touchdown
% GThTO: velocità orizzontale GT al TakeOff
% GTv2ndlast: come sopra ma verticale
% GTvTD','
% GTvTO','
% thTrTD: angolo Tronco al Touchdown
% thTrTO: angolo Tronco al TakeOff
% thThTD: angolo Tronco al Touchdown
% thShTD: angolo Tronco al Touchdown
% hGT2ndlast: altezza GT media nel penultimo passo
% hGTTD: altezza GT al Touchdown
% hGTTO
% hLGTmax
% hRGTmax
if not(isempty(jumpdata))
Tvicon.GTh2ndlast=jumpdata.Vh_mean_2ndLast;% media velocità del GT orizzontale nel penultimo passo
Tvicon.GThTD=jumpdata.Vh_TD;% velocità orizzontale GT al touchdown
Tvicon.GThTO=jumpdata.Vh_TO;% velocità orizzontale GT al TakeOff
Tvicon.GTv2ndlast=jumpdata.Vv_mean_2ndLast;% come sopra ma verticale
Tvicon.GTvTD=jumpdata.Vv_TD;%
Tvicon.GTvTO=jumpdata.Vv_TO;%
Tvicon.thTrTD=jumpdata.trunk_angle_TD;% angolo Tronco al Touchdown
Tvicon.thTrTO=jumpdata.trunk_angle_TO;% angolo Tronco al TakeOff
Tvicon.thThTD=jumpdata.thigh_angle_TD;% angolo Coscia al Touchdown
Tvicon.thShTD=jumpdata.shank_angle_TD;% angolo Gamba al Touchdown
Tvicon.hGT2ndlast=jumpdata.H_GT_mean_2ndLast;% altezza GT media nel penultimo passo
Tvicon.hGTTD=jumpdata.H_GT_TD;% altezza GT al Touchdown
Tvicon.hGTTO=jumpdata.H_GT_TO;%
Tvicon.hLGTmax=jumpdata.H_LGT_max;%
Tvicon.hRGTmax=jumpdata.H_RGT_max;%
Tvicon.secGT=jumpdata.angleGTsec;
cycles(:,9)=jumpdata.trunk_resampled;
end
%% VICON+AMTI
Tlumped=table(nan,nan,nan,nan,nan,nan,'VariableNames',{'lGT_TD','lGT_min','lGT_TO','theta_TD','theta_TO','LhGT'});

if not(isempty(theta.(side)) && isempty(stiffness.(side)))
    Tlumped.lGT_TD=stiffness.(side).L0.all(end);
    Tlumped.lGT_min=stiffness.(side).Lmin.all(end);
    Tlumped.lGT_TO=stiffness.(side).Lfin.all(end);
    Tlumped.theta_TD=theta.(side).GT.cycle(1,end);
    Tlumped.theta_TO=theta.(side).GT.cycle(end,end);
    Tlumped.LhGT=stiffness.(side).deltaAPprox.all(end);
    cycles(:,6:7)=stiffness.(side).PP.cycle(:,[2 3],ii);
end

    b=cycles(:,4:5)-cycles(:,6:7);
    b(:,3)=0;
    F=cycles(:,1:2);
    F(:,3)=0;
    Mhip=cross(b,F);
    cycles(:,8)=Mhip(:,3);

%% creating output struct
cyclesvarnames={'GRFh','GRFv','GRFmod','COPh','COPl','GTh','GTv','Mhip','TH_trunk'};
Tcycles=array2table(cycles,'VariableNames',cyclesvarnames);
sname=[Ttrial.RunID{:} '_Cycles'];
out.Parameters=[cell2table({Info.Trial.ID},'VariableNames',"TrialID") Tathl Tprost Ttrial TvideoAMTI Tvicon Tamti Tlumped];
out.(sname)=[array2table((1:101)','VariableNames',{'% Stance'}) Tcycles];

end 

%% SUBFUNS
function [Tathl, Tprost, Ttrial]=infoTable(Info)
%initialize tables
Tathl=table({'n/a'},nan,{'n/a'},nan,'VariableNames',{'Name','mass','AffSide','UnaffGTheight'});
Tprost=table({'n/a'},nan,nan,{'n/a'},nan,nan,'VariableNames',{'Model','Cat','AffGTheight','SocketConfigID','TAP','SAP'});
Ttrial=table({'n/a'},{'n/a'},{'n/a'},{'n/a'},nan,nan,'VariableNames',{'Date','RunID','StepID','StepSide','NSteps','RunLength'});


%% fill table: common
% athlete
try
Tathl.Name={Info.Athlete.ID};
Tathl.mass=double(Info.Athlete.Mass);
Tathl.AffSide={Info.Athlete.AmputationSide};
catch
    warning('No Athlete Info Available')
end

%med device, if any
if not(isempty(Info.Athlete.AmputationSide))
try
Tprost.Model=Info.MedDev.Foot.Model;
Tprost.Cat=Info.MedDev.Foot.CAT;
Tprost.SocketConfigID=Info.MedDev.Socket.ID;
Tprost.TAP=Info.MedDev.Foot.TAP;
Tprost.SAP=Info.MedDev.Socket.SAP;
catch
    warning('No prostesis data detected!')
end
end

%trial
Ttrial.Date={Info.Session.Date};
Ttrial.RunID={Info.Trial.RunID};
Ttrial.StepID="Takeoff";

end